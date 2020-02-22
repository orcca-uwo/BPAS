/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __MultiDFFT_h
#define __MultiDFFT_h 

#include "Types.h"
#include "FMUL.h"
#include "matrix.h"
#include "generalFuncs.h"
#include <stdlib.h>            /* rand() and srand() functions               */
#include <assert.h>            /* assert()                                   */



// dgs1[0], dgs2[0] are useless.
// ccumx[i] keeps the base before ith-dimension.
void plainMultiDMul(sfixn N, sfixn * ccum, sfixn * res, sfixn * ccum1, sfixn * dgs1, sfixn * ccum2, sfixn * dgs2, sfixn * coeffs1, sfixn * coeffs2, MONTP_OPT2_AS_GENE * pPtr);


void fftMultiD_square_test(sfixn * coeffs1, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);

void fftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);



void fromtofftRepMultiD(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs);

void fftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);

void tftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr);




#endif
