/*
 * Implementation of FMEDS structure and its related routines
 */
#ifndef FMEDS_H
#define FMEDS_H

#include "FME_Support_linAlg.h"
#include "FME_Support_inequality.h"
#include "FME_Support_unrolledll.h"
#include "FME_Support_inequalityOperations.h"

/*
 * This data structure contains decomposition of a linear inequality
 * system to facilitate Fourier-Motzkin elimination.
 */

typedef struct
{
	int dimension;
    unrolledLl_t * posSubSys;
    unrolledLl_t * negSubSys;
    unrolledLl_t * zerSubSys;
    inequalityNode_t ** currentPos;
    inequalityNode_t ** currentNeg;
    inequalityNode_t ** currentZer;
    int * posIndex;
	int * negIndex;
	int * zerIndex;
} FMEDS_t;

void allocFMEDS(int varNum , FMEDS_t * newds);

void FMEDSRewind(FMEDS_t * a);

void addToFMEDS(FMEDS_t * a , inequality_t * newIneq , int varIndex);

int getFMEDSSize(FMEDS_t * a);

void FMEDSListRewind(FMEDS_t * a , char n);

void freeFMEDS(FMEDS_t * inputFMEDS);

#endif

