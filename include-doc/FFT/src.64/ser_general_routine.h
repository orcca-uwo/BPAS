// basic-routine.h

#ifndef __ser_general_routine_h
#define __ser_general_routine_h

#include "modpn.h"

namespace SERBPAS {

  void checkFrourierPrime(sfixn e, MONTP_OPT2_AS_GENE * pPtr);

  //compute primitive roots
  void EX_Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr);
  
  void Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr);
  
  void Mont_GetNthRoots_OPT2_AS_GENE_SPE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr);


  //pairwise multiplication modulo a prime
  void EX_Mont_PairwiseMul_OPT2_AS(sfixn n, sfixn * APtr, sfixn * BPtr, sfixn p);

  void EX_Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void Mont_PairwiseMul_OPT2_AS_SPE_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  int forcutoff_Multi_TFT_FFT(sfixn N, sfixn *dgs1, sfixn *dgs2, int cutoff);

  void transpose(sfixn *A, sfixn lda, sfixn *B, sfixn ldb,
				 sfixn i0, sfixn i1, sfixn j0, sfixn j1);

  void sqtranspose(sfixn *A, sfixn lda,
				   sfixn i0, sfixn i1,
				   sfixn j0, sfixn dj0, sfixn j1 /*, int dj1 = 0 */);

  //for data transposition form dimension dm to dimension 1 
  void multi_mat_transpose(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data);
  void decompose(sfixn N, sfixn dm, sfixn * accum, sfixn * dims, sfixn * tmpPtr, sfixn * data);

  //MPMMTS.c
  sfixn getDenseSiz(sfixn N, preFFTRep *poly, sfixn buszsN, sfixn *data, sfixn cumN);

  void InitKroFFTRep(KroFFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr);

  void freeKroFFTRep(KroFFTRep * x);

  void InitKroTFTRep(KroTFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr);

  void freeKroTFTRep(KroTFTRep * x);

  void InitResPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs, sfixn * p2dgs);

  void copyVec_0_to_d(int d, sfixn * desV, sfixn * srcV);

  void copyVec_1_to_n(int n, sfixn * desV, sfixn * srcV);

  int compareVec(sfixn deg, sfixn * vec1, sfixn * vec2);

  int EX_IsEqualPoly(preFFTRep * Ptr1, preFFTRep * Ptr2);

  void fromtofftRepMultiD(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs);

  sfixn * reverseUni(sfixn deg, sfixn * vec1, sfixn * vec2);

  void cleanVec(sfixn deg, sfixn * cof);

  void InitOneRevInv(sfixn N, preFFTRep * tRIPtr, sfixn * bounds, sfixn di);

  sfixn *reverseMulti(sfixn deg, sfixn sizOfCoef, sfixn * outVec, sfixn * inVec);

  void PolyCleanData(preFFTRep * prePtr);

  void decreaseKroFFTRep(KroFFTRep * kPtr, sfixn * resDgs);

  void KroneckerCleanData(KroFFTRep * kPtr);

  void fromtofftRep(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs);

  void InitOneReducedPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs);

  void initCopyOneRevInv(preFFTRep * outPtr, preFFTRep * inPtr);

  sfixn *reverseMulti_1(sfixn deg, sfixn sizOfCoef, sfixn * vec);

  void subPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);

  void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p);

  void subPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);

  void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p);

  void subEqDgPoly_inner_1(sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int selector);

  void subEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p, int selector);
  
  void addEqDgPoly_inner(sfixn N, sfixn *dgs, sfixn *accum, sfixn *data1, sfixn *data2, sfixn p);

  void addEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);

  void copyPolyPointers(preFFTRep *D, preFFTRep *S);

  preFFTRep *
	getCoefMulti(sfixn N, preFFTRep* co, preFFTRep* f, sfixn i);

  void initTriRevInvSet(sfixn *dgs, sfixn N, TriRevInvSet *tRevInvPtr, TriSet * tPtr);

  TriSet * Ex_initRandomMonicDenseTriSet(sfixn N, sfixn * dgs, MONTP_OPT2_AS_GENE * pPtr);

  void initRandomDenseTriSet( sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);

  TriSet *
	EX_initRandomTriSet( sfixn N, sfixn dgbound, MONTP_OPT2_AS_GENE * pPtr);
  
  void initRandomTriSet( sfixn N, sfixn dgbound, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);

  void
	InitOneMonicRandomPolys(sfixn * bounds, sfixn N, preFFTRep * p1Ptr,  MONTP_OPT2_AS_GENE * pPtr, sfixn seed);

  void randomMonicPoly(preFFTRep * Ptr, sfixn p);

  void randomMonicPoly_inner(sfixn N, sfixn N1, sfixn * dgs, sfixn * accum, sfixn * data, sfixn p);

  preFFTRep *
	EX_randomPoly(sfixn N, sfixn * dgs, sfixn p);
  
  preFFTRep *
	EX_InitOnePoly(sfixn N, sfixn * dgs);

  void InitOnePoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs);

  void randomPoly(preFFTRep * Ptr, sfixn p);
  
  void randomPoly_inner(sfixn N, sfixn * dgs, sfixn * accum, sfixn * data, sfixn p);

  void CopyOnePoly(preFFTRep * rPtr, preFFTRep * fPtr );

  preFFTRep * CopyOnePoly_sameData(preFFTRep * fPtr );
  
  void freePoly_not_data(preFFTRep * x);
}

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


