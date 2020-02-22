
/*
 * Implementation of one Fourier-step.
 */

#ifndef FME_H
#define FME_H

#include "FME_Support_linAlg.h"
#include "FME_Support_inequality.h"
#include "FME_Support_unrolledll.h"
#include "FME_Support_inequalityOperations.h"
#include "FME_Support_fmeDS.h"
#include "FME_Support_redundancyTest.h"

/*
 * Eliminate the variable in 'varIndex' position from the 'input' structure,
 * place the result in the 'output' structure, by running FME and eliminating
 * ALL redundant inequalities with Balas' method, using extreme rays of the cone
 * defiend by 'testCone'.
 */

void oneStepFME(FMEDS_t * input, int varIndex, FMEDS_t * output, matrix_t * testCone);

#endif

