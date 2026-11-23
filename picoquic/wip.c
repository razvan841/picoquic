/*
* Author: Christian Huitema
* Copyright (c) 2019, Private Octopus, Inc.
* All rights reserved.
*
* Permission to use, copy, modify, and distribute this software for any
* purpose with or without fee is hereby granted, provided that the above
* copyright notice and this permission notice appear in all copies.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Private Octopus, Inc. BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
* "wip": work-in-progress congestion control for the thesis.
*
* For now this is a copy of picoquic's Cubic (cubic.c, without the dcubic
* variant), with a single change: in slow start the congestion window grows
* by four times per RTT instead of doubling. It exists to check that a new
* algorithm is registered and actually used ("-G wip" in picoquicdemo,
* "main_cc_algo: wip" in pico_sim).
*/

#include "picoquic_internal.h"
#include <stdlib.h>
#include <string.h>
#include "cc_common.h"
#include "picoquic_wip.h"

#define PICOQUIC_WIP_C 0.4
#define PICOQUIC_WIP_BETA_ECN (7.0 / 8.0)
#define PICOQUIC_WIP_BETA (3.0 / 4.0)

/* Growth of the window per RTT in slow start. Standard slow start adds one
 * byte to the window for each byte acknowledged, which doubles the window
 * every RTT. Adding (FACTOR - 1) bytes per acknowledged byte multiplies the
 * window by FACTOR every RTT. */
#define PICOQUIC_WIP_SLOW_START_FACTOR 4

typedef enum {
    picoquic_wip_alg_slow_start = 0,
    picoquic_wip_alg_recovery,
    picoquic_wip_alg_congestion_avoidance
} picoquic_wip_alg_state_t;

typedef struct st_picoquic_wip_state_t {
    picoquic_wip_alg_state_t alg_state;
    uint64_t recovery_sequence;
    uint64_t start_of_epoch;
    uint64_t previous_start_of_epoch;
    uint64_t previous_alg_state;
    uint64_t previous_ssthresh;
    uint64_t previous_cwin;
    double K;
    double W_max;
    double W_last_max;
    double W_reno;
    uint64_t ssthresh;
    picoquic_min_max_rtt_t rtt_filter;
} picoquic_wip_state_t;

static void wip_reset(picoquic_wip_state_t* wip_state, picoquic_path_t* path_x, uint64_t current_time) {
    memset(&wip_state->rtt_filter, 0, sizeof(picoquic_min_max_rtt_t));
    memset(wip_state, 0, sizeof(picoquic_wip_state_t));
    path_x->cwin = PICOQUIC_CWIN_INITIAL;
    wip_state->alg_state = picoquic_wip_alg_slow_start;
    wip_state->ssthresh = UINT64_MAX;
    wip_state->W_last_max = (double)wip_state->ssthresh / (double)path_x->send_mtu;
    wip_state->W_max = wip_state->W_last_max;
    wip_state->start_of_epoch = current_time;
    wip_state->previous_start_of_epoch = current_time;
    wip_state->previous_alg_state = wip_state->alg_state;
    wip_state->previous_cwin = path_x->cwin;
    wip_state->previous_ssthresh = UINT64_MAX;
    wip_state->W_reno = PICOQUIC_CWIN_INITIAL;
    wip_state->recovery_sequence = 0;
}

static void wip_init(picoquic_path_t* path_x, char const* UNUSED(option_string), uint64_t current_time)
{
    /* Initialize the state of the congestion control algorithm */
    picoquic_wip_state_t* wip_state = (picoquic_wip_state_t*)malloc(sizeof(picoquic_wip_state_t));
#ifdef _WINDOWS
    UNREFERENCED_PARAMETER(option_string);
#endif
    path_x->congestion_alg_state = (void*)wip_state;
    if (wip_state != NULL) {
        wip_reset(wip_state, path_x, current_time);
    }
}

static double wip_cubic_root(double x)
{
    /* First find an approximation */
    double v = 1;
    double y = 1.0;
    double y2;
    double y3;

    while (v > x * 8) {
        v /= 8;
        y /= 2;
    }

    while (v < x) {
        v *= 8;
        y *= 2;
    }

    for (int i = 0; i < 3; i++) {
        y2 = y * y;
        y3 = y2 * y;
        y += (x - y3) / (3.0*y2);
    }

    return y;
}

/* Compute W_cubic(t) = C * (t - K) ^ 3 + W_max */
static double wip_W_cubic(
    picoquic_wip_state_t* wip_state,
    uint64_t current_time)
{
    double delta_t_sec = ((double)(current_time - wip_state->start_of_epoch) / 1000000.0) - wip_state->K;
    double W_cubic = (PICOQUIC_WIP_C * (delta_t_sec * delta_t_sec * delta_t_sec)) + wip_state->W_max;

    return W_cubic;
}

/* On entering congestion avoidance, need to compute the new coefficients
 * of the cubic curve */
static void wip_enter_avoidance(
    picoquic_wip_state_t* wip_state,
    uint64_t current_time)
{
    wip_state->K = wip_cubic_root(wip_state->W_max*(1.0 - PICOQUIC_WIP_BETA_ECN) / PICOQUIC_WIP_C);
    wip_state->alg_state = picoquic_wip_alg_congestion_avoidance;
    wip_state->start_of_epoch = current_time;
    wip_state->previous_start_of_epoch = wip_state->start_of_epoch;
}

/* RFC 9002 7.6.2: collapse to the minimum window, and slow start up to ssthresh. */
static void wip_enter_persistent_congestion(picoquic_path_t* path_x, picoquic_wip_state_t* wip_state, uint64_t current_time)
{
    path_x->cwin = PICOQUIC_CWIN_MINIMUM;
    wip_state->previous_start_of_epoch = wip_state->start_of_epoch;
    wip_state->start_of_epoch = current_time;
    wip_state->alg_state = picoquic_wip_alg_slow_start;
}

/* The recovery state last 1 RTT, during which parameters will be frozen
 */
static void wip_enter_recovery(picoquic_cnx_t * cnx,
    picoquic_path_t* path_x,
    picoquic_congestion_notification_t notification,
    picoquic_wip_state_t* wip_state,
    uint64_t current_time)
{
    wip_state->recovery_sequence = picoquic_cc_get_sequence_number(cnx, path_x);
    wip_state->previous_start_of_epoch = wip_state->start_of_epoch;
    wip_state->previous_alg_state = wip_state->alg_state;
    wip_state->previous_ssthresh = wip_state->ssthresh;
    wip_state->previous_cwin = path_x->cwin;
    /* Update similar to new reno, but different beta */
    wip_state->W_max = (double)path_x->cwin / (double)path_x->send_mtu;
    /* Apply fast convergence */
    if (wip_state->W_max < wip_state->W_last_max) {
        wip_state->W_last_max = wip_state->W_max;
        wip_state->W_max = wip_state->W_max * PICOQUIC_WIP_BETA_ECN;
    }
    else {
        wip_state->W_last_max = wip_state->W_max;
    }
    /* Compute the new ssthresh */
    wip_state->ssthresh = (uint64_t)(wip_state->W_max * PICOQUIC_WIP_BETA_ECN * (double)path_x->send_mtu);
    if (wip_state->ssthresh < PICOQUIC_CWIN_MINIMUM) {
        /* If things are that bad, fall back to slow start */

        wip_state->alg_state = picoquic_wip_alg_slow_start;
        if (notification == picoquic_congestion_notification_persistent_congestion) {
            /* A finite ssthresh keeps the bandwidth estimate from undoing the collapse. */
            wip_state->ssthresh = PICOQUIC_CWIN_MINIMUM;
        }
        else {
            wip_state->ssthresh = UINT64_MAX;
            path_x->is_ssthresh_initialized = 0;
        }
        wip_state->previous_start_of_epoch = wip_state->start_of_epoch;
        wip_state->start_of_epoch = current_time;
        wip_state->W_reno = PICOQUIC_CWIN_MINIMUM;
        path_x->cwin = PICOQUIC_CWIN_MINIMUM;
    }
    else {
        if (notification == picoquic_congestion_notification_persistent_congestion) {
            wip_enter_persistent_congestion(path_x, wip_state, current_time);
        }
        else {
            /* Enter congestion avoidance immediately */
            wip_enter_avoidance(wip_state, current_time);
            /* Compute the initial window for both Reno and Cubic */
            double W_cubic = wip_W_cubic(wip_state, current_time);
            uint64_t win_cubic = (uint64_t)(W_cubic * (double)path_x->send_mtu);
            wip_state->W_reno = ((double)path_x->cwin) / 2.0;

            /* The formulas that compute "W_cubic" at the beginning of congestion avoidance
            * guarantee that "w_cubic" is larger than "w_reno" even if "fast convergence"
            * is applied as long as "beta_cubic" is greater than
            * (-1 + sqrt(1+4))/2, about 0.618033988749895.
            * Since beta_cubic is set to 3/4, we do not need to compare "w_cubic" and
            * "w_reno" to pick the largest. */
            path_x->cwin = win_cubic;
        }
    }
}

/* On spurious repeat notification, restore the previous congestion control.
 * If the previous state was slow start, we get back to slow start with
 * exactly the initial parameters. Otherwise,
 * we assume that K is still valid -- we only update it after exiting recovery.
 * Set cwin to the value of W_max before the recovery event
 * Set W_max to W_max_last, i.e. the value before the recovery event
 * Set the epoch back to where it was, by computing the inverse of the
 * W_cubic formula */
static void wip_correct_spurious(picoquic_path_t* path_x,
    picoquic_wip_state_t* wip_state,
    uint64_t current_time)
{
    if (wip_state->ssthresh != UINT64_MAX) {
        wip_state->W_max = wip_state->W_last_max;
        wip_state->start_of_epoch = wip_state->previous_start_of_epoch;
        wip_state->alg_state = wip_state->previous_alg_state;
        if (wip_state->alg_state != picoquic_wip_alg_slow_start) {
            wip_enter_avoidance(wip_state, wip_state->previous_start_of_epoch);
            double W_cubic = wip_W_cubic(wip_state, current_time);
            wip_state->W_reno = W_cubic * (double)path_x->send_mtu;
            wip_state->ssthresh = (uint64_t)(wip_state->W_max * PICOQUIC_WIP_BETA * (double)path_x->send_mtu);
            path_x->cwin = (uint64_t)wip_state->W_reno;
        }
        else {
            wip_state->ssthresh = wip_state->previous_ssthresh;
            path_x->cwin = wip_state->previous_cwin;
            wip_state->W_reno = (double)wip_state->previous_cwin;
        }
    }
}

static void wip_notify(
    picoquic_cnx_t* cnx, picoquic_path_t* path_x,
    picoquic_congestion_notification_t notification,
    picoquic_per_ack_state_t * ack_state,
    uint64_t current_time)
{
    picoquic_wip_state_t* wip_state = (picoquic_wip_state_t*)path_x->congestion_alg_state;
    path_x->is_cc_data_updated = 1;

    if (wip_state != NULL) {
        switch (notification) {
            /* RTT measurements will happen before acknowledgement is signalled */
            case picoquic_congestion_notification_acknowledgement:
                switch (wip_state->alg_state) {
                    case picoquic_wip_alg_slow_start:
                        if (wip_state->ssthresh == UINT64_MAX) {
                            /* Increase cwin based on bandwidth estimation. */
                            path_x->cwin = picoquic_cc_update_target_cwin_estimation(path_x);
                        }

                        if (path_x->last_time_acked_data_frame_sent > path_x->last_sender_limited_time) {
                            /* The only change from Cubic: the increase is multiplied by
                             * (FACTOR - 1), so the window quadruples every RTT instead of
                             * doubling. The increase is still zero when the path was not
                             * cwin limited (see picoquic_cc_slow_start_increase). */
                            path_x->cwin += (PICOQUIC_WIP_SLOW_START_FACTOR - 1) *
                                picoquic_cc_slow_start_increase_ex(path_x, ack_state->nb_bytes_acknowledged, 0);

                            /* if cnx->cwin exceeds SSTHRESH, exit and go to CA */
                            if (path_x->cwin >= wip_state->ssthresh) {
                                wip_state->W_reno = ((double)path_x->cwin) / 2.0;
                                path_x->is_ssthresh_initialized = 1;
                                wip_enter_avoidance(wip_state, current_time);
                            }
                        }
                        break;
                    case picoquic_wip_alg_recovery:
                        /* exit recovery, move to CA or SS, depending on CWIN */
                        wip_state->alg_state = picoquic_wip_alg_slow_start;
                        path_x->cwin += ack_state->nb_bytes_acknowledged;
                        /* if cnx->cwin exceeds SSTHRESH, exit and go to CA */
                        if (path_x->cwin >= wip_state->ssthresh) {
                            wip_state->alg_state = picoquic_wip_alg_congestion_avoidance;
                        }
                        break;
                    case picoquic_wip_alg_congestion_avoidance:
                        if (path_x->last_time_acked_data_frame_sent > path_x->last_sender_limited_time) {
                            double W_cubic;
                            uint64_t win_cubic;
                            /* Protection against limited senders. */
                            if (wip_state->start_of_epoch < path_x->last_sender_limited_time) {
                                wip_state->start_of_epoch = path_x->last_sender_limited_time;
                            }
                            /* Compute the cubic formula */
                            W_cubic = wip_W_cubic(wip_state, current_time);
                            win_cubic = (uint64_t)(W_cubic * (double)path_x->send_mtu);
                            /* Also compute the Reno formula */
                            wip_state->W_reno += ((double)ack_state->nb_bytes_acknowledged) * ((double)path_x->send_mtu) / wip_state->W_reno;

                            /* Pick the largest */
                            if ((double)win_cubic > wip_state->W_reno) {
                                /* if cubic is larger than threshold, switch to cubic mode */
                                path_x->cwin = win_cubic;
                            }
                            else {
                                path_x->cwin = (uint64_t)wip_state->W_reno;
                            }
                        }
                        break;
                }
                break;
            case picoquic_congestion_notification_repeat:
            case picoquic_congestion_notification_timeout:
            case picoquic_congestion_notification_ecn_ec:
                switch (wip_state->alg_state) {
                    case picoquic_wip_alg_slow_start:
                        /* For compatibility with Linux-TCP deployments, we implement a filter so
                         * Cubic will only back off after repeated losses, not just after a single loss.
                         */
                        if ((notification == picoquic_congestion_notification_ecn_ec ||
                            picoquic_cc_hystart_loss_test(&wip_state->rtt_filter, notification, ack_state->lost_packet_number, PICOQUIC_SMOOTHED_LOSS_THRESHOLD)) &&
                            (current_time - wip_state->start_of_epoch > path_x->smoothed_rtt ||
                                wip_state->recovery_sequence <= picoquic_cc_get_ack_number(cnx, path_x))) {
                            path_x->is_ssthresh_initialized = 1;
                            wip_enter_recovery(cnx, path_x, notification, wip_state, current_time);
                        }
                        break;
                    case picoquic_wip_alg_recovery:
                    case picoquic_wip_alg_congestion_avoidance:
                        /* For compatibility with Linux-TCP deployments, we implement a filter so
                         * Cubic will only back off after repeated losses, not just after a single loss.
                         */
                        if (ack_state->lost_packet_number >= wip_state->recovery_sequence &&
                            (notification == picoquic_congestion_notification_ecn_ec ||
                                picoquic_cc_hystart_loss_test(&wip_state->rtt_filter, notification, ack_state->lost_packet_number, PICOQUIC_SMOOTHED_LOSS_THRESHOLD))) {
                            /* Re-enter recovery */
                            wip_enter_recovery(cnx, path_x, notification, wip_state, current_time);
                        }
                        break;
                }
                break;
            case picoquic_congestion_notification_persistent_congestion:
                /* Not filtered: if the losses of the same ACK entered recovery, keep their ssthresh and only collapse the window. */
                if (ack_state->lost_packet_number >= wip_state->recovery_sequence) {
                    path_x->is_ssthresh_initialized = 1;
                    wip_enter_recovery(cnx, path_x, notification, wip_state, current_time);
                }
                else {
                    wip_enter_persistent_congestion(path_x, wip_state, current_time);
                }
                break;
            case picoquic_congestion_notification_spurious_repeat:
                /* Reset CWIN based on ssthresh, not based on current value. */
                wip_correct_spurious(path_x, wip_state, current_time);
                break;
            case picoquic_congestion_notification_rtt_measurement:
                if (wip_state->alg_state == picoquic_wip_alg_slow_start &&
                    wip_state->ssthresh == UINT64_MAX) {

                    /* HyStart. */
                    /* Using RTT increases as signal to get out of initial slow start */
                    if (picoquic_cc_hystart_test(&wip_state->rtt_filter, (cnx->is_time_stamp_enabled) ? ack_state->one_way_delay : ack_state->rtt_measurement,
                            cnx->path[0]->pacing.packet_time_microsec, current_time, cnx->is_time_stamp_enabled)) {
                        /* RTT increased too much, get out of slow start! */

                        if (wip_state->rtt_filter.rtt_filtered_min > PICOQUIC_TARGET_RENO_RTT){
                            double correction;
                            if (wip_state->rtt_filter.rtt_filtered_min > PICOQUIC_TARGET_SATELLITE_RTT) {
                                correction = (double)PICOQUIC_TARGET_SATELLITE_RTT / (double)wip_state->rtt_filter.rtt_filtered_min;
                            }
                            else {
                                correction = (double)PICOQUIC_TARGET_RENO_RTT / (double)wip_state->rtt_filter.rtt_filtered_min;
                            }
                            uint64_t base_window = (uint64_t)(correction * (double)path_x->cwin);
                            uint64_t delta_window = path_x->cwin - base_window;
                            path_x->cwin -= (delta_window / 2);
                        }
                        else {
                            /* In the general case, compensate for the growth of the window after the acknowledged packet was sent. */
                            path_x->cwin /= 2;
                        }

                        wip_state->ssthresh = path_x->cwin;
                        wip_state->W_max = (double)path_x->cwin / (double)path_x->send_mtu;
                        wip_state->W_last_max = wip_state->W_max;
                        wip_state->W_reno = ((double)path_x->cwin);
                        path_x->is_ssthresh_initialized = 1;
                        /* enter recovery to ignore the losses expected if the window grew
                        * too large after the acknowleded packet was sent. */
                        wip_enter_recovery(cnx, path_x, notification, wip_state, current_time);
                        /* apply a correction to enter the test phase immediately */
                        uint64_t K_micro = (uint64_t)(wip_state->K * 1000000.0);
                        if (K_micro > current_time) {
                            wip_state->K = ((double)current_time) / 1000000.0;
                            wip_state->start_of_epoch = 0;
                        }
                        else {
                            wip_state->start_of_epoch = current_time - K_micro;
                        }
                    }
                }
                break;
            case picoquic_congestion_notification_seed_cwin:
                if (wip_state->alg_state == picoquic_wip_alg_slow_start) {
                    if (wip_state->ssthresh == UINT64_MAX) {
                        if (path_x->cwin < ack_state->nb_bytes_acknowledged) {
                            path_x->cwin = ack_state->nb_bytes_acknowledged;
                        }
                        wip_state->ssthresh = ack_state->nb_bytes_acknowledged;
                        wip_state->W_max = (double)path_x->cwin / (double)path_x->send_mtu;
                        wip_state->W_last_max = wip_state->W_max;
                        wip_state->W_reno = ((double)path_x->cwin);
                        path_x->is_ssthresh_initialized = 1;
                        wip_enter_avoidance(wip_state, current_time);
                    }
                }
                break;
            case picoquic_congestion_notification_reset:
                wip_reset(wip_state, path_x, current_time);
                break;
            default:
                break;

        }

        /* Compute pacing data */
        picoquic_update_pacing_data(path_x, wip_state->alg_state == picoquic_wip_alg_slow_start &&
            wip_state->ssthresh == UINT64_MAX);
    }
}

/* Release the state of the congestion control algorithm */
static void wip_delete(picoquic_path_t* path_x)
{
    if (path_x->congestion_alg_state != NULL) {
        free(path_x->congestion_alg_state);
        path_x->congestion_alg_state = NULL;
    }
}

/* Observe the state of congestion control */
static void wip_observe(picoquic_path_t* path_x, uint64_t* cc_state, uint64_t* cc_param)
{
    picoquic_wip_state_t* wip_state = (picoquic_wip_state_t*)path_x->congestion_alg_state;
    *cc_state = (uint64_t)wip_state->alg_state;
    *cc_param = (uint64_t)wip_state->W_max;
}

/* Definition record for the wip algorithm */

#define picoquic_wip_ID "wip"

picoquic_congestion_algorithm_t picoquic_wip_algorithm_struct = {
    picoquic_wip_ID, PICOQUIC_CC_ALGO_NUMBER_WIP, PICOQUIC_ECN_ECT_0,
    wip_init,
    wip_notify,
    wip_delete,
    wip_observe
};

picoquic_congestion_algorithm_t* picoquic_wip_algorithm = &picoquic_wip_algorithm_struct;
