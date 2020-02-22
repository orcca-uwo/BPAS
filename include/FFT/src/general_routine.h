// basic-routine.h

#ifndef __general_routine_h
#define __general_routine_h

#include "modpn.h"
#if FURER
#include "fft_furer1.h"
#include "fft_furer2.h"
#else
#include "fft_iter1.h"
#include "fft_iter2.h"
#endif


namespace PBPAS {

  void checkFrourierPrime(sfixn e, MONTP_OPT2_AS_GENE * pPtr);

  int forcutoff_Multi_TFT_FFT(sfixn N, sfixn *dgs1, sfixn *dgs2, int cutoff);


  //for data transposition form dimension dm to dimension 1 
  void multi_mat_transpose(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data);
  void decompose(sfixn N, sfixn dm, sfixn * accum, sfixn * dims, sfixn * tmpPtr, sfixn * data);

  //MPMMTS.c
  sfixn getDenseSiz(sfixn N, preFFTRep *poly, sfixn buszsN, sfixn *data, sfixn cumN);

  void InitKroFFTRep(KroFFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr);

  void freeKroFFTRep(KroFFTRep * x);

  void InitKroTFTRep(KroTFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr);

  void InitKroTFTRep_multiplier_1(KroTFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M, MONTP_OPT2_AS_GENE * pPtr);

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

  void subEqDgPoly_inner_1(sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int selector);

  void subEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p, int selector);

  void addEqDgPoly_inner(sfixn N, sfixn *dgs, sfixn *accum, sfixn *data1, sfixn *data2, sfixn p);

  void addEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);

  void copyPolyPointers(preFFTRep *D, preFFTRep *S);

  preFFTRep *
	getCoefMulti(sfixn N, preFFTRep* co, preFFTRep* f, sfixn i);

  void initTriRevInvSet(sfixn *dgs, sfixn N, TriRevInvSet *tRevInvPtr, TriSet * tPtr);

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

  sfixn noExtensionBalancedBivar(sfixn N, sfixn* ls);

  void ReductionToBalancedBivar(sfixn N, sfixn* ls, sfixn *factm);

  void prime_factorization(sfixn x, sfixn fact[]);
  
  sfixn V1to2V_pivot(sfixn d1, sfixn d2);
  sfixn multi_V1to2V_pivot(sfixn N, sfixn * degs1, sfixn *degs2);
  
  sfixn multi_V1to2V_pivot_opt(sfixn sigma, sfixn q1f, sfixn q1g);

  void InitKroTFTRep_1V2V(KroTFTRep * kPtr, sfixn b, sfixn d1, sfixn d2);

  void InitKroTFTRep_1V2V_multi(KroTFTRep *kPtr, sfixn b, sfixn sigma, sfixn q1f, sfixn q1g );

  void from1Vto2VTFTRep(sfixn b, sfixn d, sfixn * coeffs, sfixn rccum2, sfixn * data);
  void fromMulti_1Vto2VTFTRep_1(sfixn N, sfixn b, sfixn q1, sfixn *pdegs, sfixn *M1, sfixn num_boxes_B, sfixn num_blocks_B, sfixn *B);

  void fromMulti_1Vto2VTFTRep(sfixn n, sfixn b, sfixn q1f, sfixn q1g, sfixn *degs_M1, sfixn *M1, sfixn *degs_M2, sfixn *B);

  void from2VTFTRepto1V(sfixn b, sfixn s1, sfixn s2, sfixn * data, sfixn d, sfixn * coeffs, MONTP_OPT2_AS_GENE *pPtr);

  void from2VTFTRepto1V_2Traversal(sfixn b, sfixn s1, sfixn s2, sfixn * data, sfixn siz, sfixn * coeffs, MONTP_OPT2_AS_GENE *pPtr);
  
  void from2VTFTReptoMultiV(sfixn b, sfixn sigma, sfixn num_block_B, sfixn *B, sfixn box_size_M, sfixn *M, MONTP_OPT2_AS_GENE *pPtr);

  void nextPosition(sfixn * exp, sfixn *deg, sfixn n, sfixn * size);

  void multi_mat_transpose_2DTran_interp(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data);

  void multi_mat_transpose_2DTran_eval(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data);


  void EX_Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr,
                                        MONTP_OPT2_AS_GENE * pPtr);

  void Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr,
                                     MONTP_OPT2_AS_GENE * pPtr);

  void Mont_GetNthRoots_OPT2_AS_GENE_SPE(sfixn e, sfixn n, sfixn * rootsPtr,
                                         MONTP_OPT2_AS_GENE * pPtr);

  void EX_Mont_PairwiseMul_OPT2_AS(sfixn n, sfixn * APtr, sfixn * BPtr, sfixn p);

  void EX_Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr,
                                     MONTP_OPT2_AS_GENE * pPtr);

  void Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr,
                                  MONTP_OPT2_AS_GENE * pPtr);

  void Mont_PairwiseMul_OPT2_AS_SPE_R(sfixn n, sfixn * APtr, sfixn * BPtr,
                                      MONTP_OPT2_AS_GENE * pPtr);

  void transpose(sfixn *A, sfixn lda, sfixn *B, sfixn ldb,
                 sfixn i0, sfixn i1, sfixn j0, sfixn j1);

  void sqtranspose(sfixn *A, sfixn lda,
                   sfixn i0, sfixn i1,
                   sfixn j0, sfixn dj0, sfixn j1 /*, int dj1 = 0 */);

  int ceil_log2_long(unsigned long long x);

  //new DFT
  void DFT_eff(int n, int r,
               sfixn *A,
               sfixn *W,
               MONTP_OPT2_AS_GENE *pPtr,
               int H, int *RevBidMap,
               sfixn *B,sfixn num);

  void DFT_rec(int n, int r,
               sfixn *A,
               sfixn *W,
               MONTP_OPT2_AS_GENE *pPtr,
               int H, int *RevBidMap,
               sfixn *B);

  void DFT_iter(int n, int r,
                sfixn *A,
                sfixn *W,
                MONTP_OPT2_AS_GENE *pPtr,
                sfixn *B);

  void Shuffle(int n, sfixn *A, sfixn *B);

  void ArrayBitReversal(int n, sfixn *A,
                        int *RevBidMap);

  void RevBitIncr( int *n, int bit );

  void RevBitInd(int max, int *A);

  void RootsTable(int n, int r,
                  sfixn *T,
                  MONTP_OPT2_AS_GENE *pPtr);

  void InverseRootsTable(int n, sfixn *T, sfixn *Ti);
  void RootsTableSpe(int n, int r,
                  sfixn *T,MONTP_OPT2_AS_GENE *p,
                  int m);

  void RootsTable2(int n, int r,
                   sfixn *T, sfixn *Ti,
                   MONTP_OPT2_AS_GENE *pPtr);

  void Weight_RootsTable2(int n, int r, int K,
                          sfixn *Th, sfixn *T, sfixn *Ti,
                          MONTP_OPT2_AS_GENE *pPtr);

  void InvDFT_eff(int n, int r,
                  sfixn *A,
                  sfixn *W,
                  MONTP_OPT2_AS_GENE *pPtr,
                  int H, int *RevBidMap,
                  sfixn *B,sfixn invn,sfixn num);
  void InvDFT_eff_keepMontgomery(int n, int r,
                  sfixn *A,
                  sfixn *W,
                  MONTP_OPT2_AS_GENE *pPtr,
                  int H, int *RevBidMap,
                  sfixn *B,sfixn invn,sfixn num);

}

#endif
