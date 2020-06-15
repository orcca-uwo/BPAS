#ifndef FME_DNC_H
#define FME_DNC_H

#include "FME_Support_linAlg.h"
#include "FME_Support_inequality.h"
#include "FME_Support_unrolledll_dnc.h"
#include "FME_Support_inequalityOperations.h"
#include "FME_Support_fmeDS_dnc.h"
#include "FME_Support_redundancyTest.h"


void oneStepFME__dncSupp(FMEDS_dnc_t * input, int varIndex, FMEDS_dnc_t * output,
		matrix_t * testCone);
		
void oneStepFME_dnc(FMEDS_dnc_t * input, int varIndex, FMEDS_dnc_t * output,
		matrix_t * testCone, int thr);
#endif
