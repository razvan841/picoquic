/*
* Work-in-progress congestion control for the thesis (see wip.c).
* Same licence terms as the rest of picoquic.
*/

#ifndef PICOQUIC_WIP_H
#define PICOQUIC_WIP_H

#include "picoquic.h"

#ifdef __cplusplus
extern "C" {
#endif

extern picoquic_congestion_algorithm_t* picoquic_wip_algorithm;

#ifdef __cplusplus
}
#endif
#endif
