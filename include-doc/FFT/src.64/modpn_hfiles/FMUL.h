/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __FMUL_h
#define __FMUL_h 

#include "Types.h"
#include "generalFuncs.h"
#include "AS.h"
#include <stdlib.h>            /* rand() and srand() functions               */
#include <assert.h>            /* assert()                                   */

#define CutOffPlainFFTMul 138



extern sfixn BASE;




//==========================================================
// Exported functions:



void EX_Mont_TFTMul_OPT2_AS_GENE(sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);


void EX_Mont_PlainMul_OPT2_AS_GENE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void EX_MontP_Init_OPT2_AS_GENE(MONTP_OPT2_AS_GENE * pPtr, sfixn p);

void EX_MontP_Free_OPT2_AS_GENE(MONTP_OPT2_AS_GENE * pPtr);

void EX_MontP_Print_OPT2_AS_GENE(MONTP_OPT2_AS_GENE * pPtr);


// This get roots can't used by other FFTs (ONLY for OPT2_AS_GENE), since the roots consists of a factor of R and shifted left BASE_Rpow bits.
void EX_Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr);

void EX_Mont_GetNthRoots_OPT2_AS_GENE_RAND(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr);


// all function with suffix _R suppose the input or ouput contains a R^-1 besides the correct result.
void 
EX_Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

void 
EX_Mont_PairwiseMul_OPT2_AS(sfixn n, sfixn * APtr, sfixn * BPtr, sfixn p);



void * BlockPairwiseMul(void * PTR);

void 
EX_Mont_DFT_OPT2_AS_GENE ( sfixn n, sfixn power, sfixn * rootsPtr,  sfixn * tmpVecPtr, sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr);
void 
EX_Mont_INVDFT_OPT2_AS_GENE_R (sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, sfixn degRes, sfixn * ResPtr, MONTP_OPT2_AS_GENE * pPtr );
void 
EX_Mont_INVDFT_OPT2_AS_GENE (sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, sfixn degRes, sfixn * ResPtr, MONTP_OPT2_AS_GENE * pPtr );

void
EX_Mont_DFT_OPT2_AS_GENE_1 ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, MONTP_OPT2_AS_GENE * pPtr );

void 
EX_Mont_INVDFT_OPT2_AS_GENE_R_1 ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, MONTP_OPT2_AS_GENE * pPtr);

void
EX_Mont_INVDFT_OPT2_AS_GENE_1 ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, MONTP_OPT2_AS_GENE * pPtr);



void EX_Mont_FFTMul_OPT2_AS_GENE(sfixn n, sfixn e, sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

void EX_Mont_FFTSQUARE_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr);

void EX_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

void EX_Mont_FFTMul(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);


void EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);


void EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr);


void EX_Mont_FFTMul_OPT2_AS_GENE_1_2(sfixn n, sfixn e, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

void Mont_TFTMul_OPT2_AS_GENE(sfixn * resPtr, sfixn d1, sfixn * v1Ptr, sfixn d2, sfixn * v2Ptr, MONTP_OPT2_AS_GENE * pPtr);

void Mont_TFTMul_OPT2_AS_GENE_SPE(sfixn * resPtr, sfixn d1, sfixn * v1Ptr, sfixn d2, sfixn * v2Ptr, MONTP_OPT2_AS_GENE * pPtr);


void EX_Mont_TFTMul_OPT2_AS_GENE_WHOLE(sfixn * resPtr, sfixn d1, sfixn * v1Ptr, sfixn d2, sfixn * v2Ptr, MONTP_OPT2_AS_GENE * pPtr);


void 
EX_Mont_TDFT_OPT2_AS_GENE_1 ( sfixn l, sfixn * rootsPtr, sfixn * tmpVecPtr, MONTP_OPT2_AS_GENE * pPtr);


void 
EX_Mont_INVTDFT_OPT2_AS_GENE_1 ( sfixn l, sfixn * rootsPtr, sfixn * tmpVecPtr, MONTP_OPT2_AS_GENE * pPtr);


void 
EX_Mont_INVTDFT_OPT2_AS_GENE_R_1 ( sfixn l, sfixn * rootsPtr, sfixn * tmpVecPtr, MONTP_OPT2_AS_GENE * pPtr);

sfixn *EX_Mont_Mul(sfixn *degResAddr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);


#endif
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


