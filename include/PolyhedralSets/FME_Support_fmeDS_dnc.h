/*
 * Implementation of FMEDS structure and its related routines
 */
#ifndef FMEDS_DNC_H
#define FMEDS_DNC_H

#include "FME_Support_linAlg.h"
#include "FME_Support_inequality.h"
#include "FME_Support_unrolledll_dnc.h"
#include "FME_Support_redundancyTest.h"
#include "FME_Support_inequalityOperations.h"

/*
 * This data structure contains decomposition of a linear inequality
 * system to facilitate Fourier-Motzkin elimination.
 */

typedef struct
{
    int dimension;
    unrolledLl_dnc_t * posSubSys;
    unrolledLl_dnc_t * negSubSys;
    unrolledLl_dnc_t * zerSubSys;
    inequalityNode_dnc_t ** currentPos;
    inequalityNode_dnc_t ** currentNeg;
    inequalityNode_dnc_t ** currentZer;
    int * posIndex;
    int * negIndex;
    int * zerIndex;
} FMEDS_dnc_t;

void allocFMEDS_dnc(int varNum , FMEDS_dnc_t * newds);

void FMEDSRewind_dnc(FMEDS_dnc_t * a);

void addToFMEDS_dnc(FMEDS_dnc_t * a , inequality_t * newIneq , int varIndex);

int getFMEDSSize_dnc(FMEDS_dnc_t * a);

void FMEDSListRewind_dnc(FMEDS_dnc_t * a , char n);

void freeFMEDS_dnc(FMEDS_dnc_t * inputFMEDS);

void addZeros_dnc(FMEDS_dnc_t * input, int varIndex, FMEDS_dnc_t * output, 
		matrix_t * testCone);
		
void divideFMEDS_dnc(FMEDS_dnc_t * input, FMEDS_dnc_t * out1, FMEDS_dnc_t * out2);

void mergeFMEDS_dnc(FMEDS_dnc_t * out1, FMEDS_dnc_t * out2, FMEDS_dnc_t * out3, FMEDS_dnc_t * out4,
		FMEDS_dnc_t * output);

#endif

