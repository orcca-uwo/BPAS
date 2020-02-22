// basic-routine.h

#ifndef __basic_routine_h
#define __basic_routine_h

#include "modpn.h"
#include "modpn_export.h"
#include "general_routine.h"
#if FURER
#include "fft_furer1.h"
#else
#include "fft_iter1.h"
#endif


namespace PBPAS {

  /**MultiDFFT.c
   *multiply two polys (defined by their coeffs arrays) by multi-FFT
   *in place
   * seperated evaluation by FFT
   *result is in coeffs1
   **/
  void fftMultiD_test_1(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr);
  
  sfixn * MultiEvalByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr);
  
  void InterpolByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr1);
  
  /**MultiDFFT.c
   *multiply two polys (defined by their coeffs arrays) by multi-TFT
   *use a temp vector for TFT
   *result in coeffs1
   *seperated evaluation by TFT
   **/
  void tftMultiD_test_1(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr);

  void EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(sfixn * coeffs1j, sfixn maxdim, sfixn lsi, sfixn * tmprootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  sfixn * MultiEvalByTFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr);

  void EX_Mont_TDFT_OPT2_AS_GENE_1_par(sfixn lsi, sfixn * tmprootsPtr, sfixn maxdim, sfixn * coeffs1j, MONTP_OPT2_AS_GENE *pPtr);

  void InterpolByTFT(sfixn *coeffs1, sfixn N, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr);


  void MultiplyByFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_RBB(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_3D2(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_2D2(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);


  void MultiplyByKroneckerFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  //----plain multiply
  void plainMultiDMul(sfixn N, sfixn * ccum, sfixn * res, sfixn * ccum1, sfixn * dgs1, sfixn * ccum2, sfixn * dgs2, sfixn * coeffs1, sfixn * coeffs2, MONTP_OPT2_AS_GENE * pPtr);

  void decomposePoly(sfixn N, sfixn * ccum, sfixn * res, sfixn N1, sfixn * dgs1, sfixn * coeffs1, sfixn * ccum1, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr, sfixn R, sfixn SFT);

  void decomposePoly2(sfixn N, sfixn * ccum, sfixn * res, sfixn num, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr);

  void UniFastMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, sfixn Lbound, sfixn * BRevInvPtr,  MONTP_OPT2_AS_GENE * pPtr);

  void UniPlainMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

  void EX_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void EX_Mont_FFTMul_OPT2_AS_GENE(sfixn n, sfixn e, sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void getRevInvTiSet(sfixn *dgs, sfixn N, TriRevInvSet * tris, TriSet * ts,  MONTP_OPT2_AS_GENE * pPtr);

  void NewtonRevInverse(sfixn N, preFFTRep * tRIPtrtmp, preFFTRep * tRIPtr, preFFTRep * tPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);


  sfixn * modularInvPM(sfixn degG, // degG=n-1;
				 sfixn * GPtr, sfixn degF, sfixn * FPtr, 
				 sfixn r, sfixn n,
				 MONTP_OPT2_AS_GENE * pPtr);
 
  void squarePoly_FFT(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep *p1Ptr,  MONTP_OPT2_AS_GENE * pPtr);
 
  void fftMultiD_square_test_1(sfixn * coeffs1, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);

  void EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr);

  void
	MultiMod(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void
	MultiMod_DF(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiMod_1_par(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiUniFastMod_1(sfixn N, preFFTRep * tmpIn,  sfixn n, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void EX_mulPoly(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);

  void EX_mulPoly_TFTFFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr);

  void EX_mulPoly_FFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);
  
  void mulPoly_FFT(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr);

  void polyMul_TFT(sfixn N, KroTFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr);

  void
	MultiMod_BULL(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiMod_1_BULL_par(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiUniFastMod_1_TFT(sfixn N, preFFTRep * tmpIn,  sfixn n, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  preFFTRep *
	MultiUniPlainMod_1(sfixn N, sfixn d1, preFFTRep * FPtr, sfixn d2, preFFTRep * GPtr, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);
  
  sfixn
	fmecgEqDg_1(sfixn N, preFFTRep * f1, sfixn e, preFFTRep * co, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

  preFFTRep * 
	mul_Reduced(preFFTRep * resPtr, sfixn N, preFFTRep * out, preFFTRep * f1, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);
  
  void EX_Mont_PlainMul_OPT2_AS_GENE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

  void
	Mont_PlainMul_OPT2_AS_GENE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

  void
	Mont_PlainMul_OPT2_AS_GENE_SPE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

  void reduceCoeffs(sfixn N, preFFTRep * toPtr, preFFTRep * fromPtr, TriSet * ts, TriRevInvSet * tris,   MONTP_OPT2_AS_GENE * pPtr);

  TriRevInvSet *
	EX_getRevInvTriSet(sfixn N,  TriSet * ts,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiplyByTFT_1V2V(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_RBBnoE(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByKronecker1DTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_RBB_KN1D_to_2D(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_1Vto2V_multiV_2V(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void tft2D_no_transp(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE * pPtr);

  void InterpolBy2DTFT_no_transp(sfixn *coeffs1, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr);

  void inv_TDFT_dim2_no_transp(sfixn d, sfixn ls2, sfixn * tmprootsPtr, sfixn maxdim, sfixn * coeffs1, MONTP_OPT2_AS_GENE *pPtr);

  void MultiEvalBy2DTFT_no_transp(sfixn * coeffs1, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr);

  void TDFT_dim2_no_transp(sfixn d, sfixn ls2, sfixn * tmprootsPtr, sfixn maxdim, sfixn * coeffs1, MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyBy2DTFT_no_transp(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void bivarMultiplyBy2DTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void  bivarEvalBy2DTFT(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn pd1, sfixn pd2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void bivarInterpolBy2DTFT(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void V1to2VEvalBy2DTFT(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn d, sfixn b, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void bivarMultiplyBy2DTFT_spawn(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void bivarMultiplyBy2DTFT_for_spawn(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void bivarInterpolBy2DTFT_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void bivarInterpolBy2DTFT_for_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void copyDataOutInvTFT(sfixn * data, sfixn dims, size_t size, sfixn ls, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void bivarEvalBy2DTFT_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn * pdegs, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void bivarEvalBy2DTFT_for_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn * pdegs, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void copyDataOutTFT(sfixn * indata, sfixn * outdata, sfixn dims, size_t insize, size_t outsize, sfixn ls, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void InterpolByTFT_2DTran(sfixn *coeffs1, sfixn N, sfixn n, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr);
  
  sfixn * MultiEvalByTFT_2DTran(sfixn * coeffs1, sfixn N, sfixn n, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr);

  void tftMultiD_test_1_2DTran(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiplyByTFT_2DTran(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_2DTran_noKroRep(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  sfixn * MultiEvalByTFT_2DTran_saveOn1stVar(sfixn * resCoeffs, sfixn N, sfixn n, sfixn *dims, sfixn * ls, sfixn * rccum, sfixn * coeffs, sfixn n1, sfixn * pdegs, sfixn * ccum, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void eval1stVar(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum, sfixn * dgs, sfixn * coeffs, sfixn dim1, sfixn ls1, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void MultiplyByTFT_RBBnoE_saveOnDim1(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void MultiEvalByTFT_RBBnoE_saveOnDim1(sfixn N, sfixn t, sfixn *resCoeffs, sfixn *rccum, sfixn ls1, sfixn ls2, sfixn dim1, sfixn dim2, sfixn *coeffs, sfixn *dgs, sfixn *ccum, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void evalDim1(sfixn N, sfixn t, sfixn *resCoeffs, sfixn *rccum, sfixn ls1, sfixn dim1, sfixn *coeffs, sfixn *dgs, sfixn *ccum, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  sfixn * IntPolyMulMod_2DTFT(sfixn * f1, sfixn * f2, sfixn d1, sfixn d2, sfixn k, sfixn p);

  void WeightVectorEval(sfixn *thetaPtr, sfixn *A, sfixn *AP, sfixn *res, 		
			sfixn K, sfixn d, sfixn dims2, sfixn ls2, 
			sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);


// In the above function, the last argument "num"
// can take values 1, 2 or any another sfixn value.
// The values 1 and 2 correspond to the first and second
// prime numbers defined in the configuration file
// located in the directory where the 1D FFT code generator resides.
// When these primes are selected, the 1D FFT code which is run 
// is that generated when the BPAS is compiled. In this
// automatically generated code, many code optimization techniques
// are used and, in particular, the prime number is hard cooded.
// When "num" has value different from 1 and 2, the prime
// number is defined by the data structure pPtr; in this case
// the 1D FFT code is not automatically generated and thus it
// is not specialized at a particular prime.
    sfixn *TwoConvolutionModNew(sfixn *Ap, sfixn *Bp,
                           sfixn d1, sfixn d2, sfixn K,
                           MONTP_OPT2_AS_GENE *pPtr,
                           int H, int *RevBidMap,sfixn num);

  void CyclicConvolution(sfixn *res, sfixn s,
                         sfixn es1, sfixn es2,
                         sfixn K, sfixn dims2,
                         sfixn *A, sfixn *B, sfixn dA, sfixn dB,
                         sfixn *KRT, sfixn *dRT,
                         sfixn *invKRT, sfixn *invdRT,
                         MONTP_OPT2_AS_GENE *pPtr,
                         int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num);

  void AdaptiveEvaluation(sfixn *res, sfixn es1, sfixn es2,
                          sfixn K, sfixn dims2,
                          sfixn *A, sfixn dA,
                          sfixn *KRT, sfixn *dRT,
                          MONTP_OPT2_AS_GENE *pPtr,
                          int H, int *RevBidMap,sfixn num);

  void AdaptiveInterpolation(sfixn *res, sfixn K, sfixn es1,
                             sfixn es2, sfixn dims2,
                             sfixn *invKRT, sfixn *invdRT,
                             MONTP_OPT2_AS_GENE *pPtr,
                             int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num);

  void NegacyclicConvolution(sfixn *res, sfixn s,
                             sfixn es1, sfixn es2,
                             sfixn K, sfixn dims2,
                             sfixn *A, sfixn *B,
                             sfixn dA, sfixn dB, sfixn K2,
                             sfixn *thetaPtr,
                             sfixn *KRT, sfixn *dRT,
                             sfixn *invKRT, sfixn *invdRT,
                             MONTP_OPT2_AS_GENE *pPtr,
                             int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num);

  void WeightVectorAdaptiveEvaluation(sfixn *res, sfixn es1,
                                      sfixn es2, sfixn K,
                                      sfixn dims2,
                                      sfixn *A, sfixn dA,
                                      sfixn *thetaPtr,
                                      sfixn *KRT, sfixn *dRT,
                                      MONTP_OPT2_AS_GENE *pPtr,
                                      int H, int *RevBidMap,sfixn num);

  void AdaptiveInterpolationWeightVector(sfixn *res, sfixn K,
                                         sfixn es1, sfixn es2,
                                         sfixn dims2,
                                         sfixn *thetaPtr,
                                         sfixn *invKRT,
                                         sfixn *invdRT,
                                         MONTP_OPT2_AS_GENE *pPtr,
                                         int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num);


  //old-------------------------------------------------
  sfixn * TwoConvolutionMod(sfixn *Ap, sfixn *Bp,
                            sfixn d1, sfixn d2,
                            sfixn p, sfixn K);

  void AdaptiveEvaluation1(sfixn *res, sfixn es1, sfixn es2,
                          sfixn K, sfixn dims2, sfixn ls2,
                          sfixn *A, sfixn dA,
                          sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void AdaptiveInterpolation1(sfixn *res, sfixn K, sfixn es1,
                             sfixn es2, sfixn dims2, sfixn ls2,
                             sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void CyclicConvolution1(sfixn *res, sfixn s,
                         sfixn es1, sfixn es2, sfixn K, sfixn dims2, sfixn ls2,
                         sfixn *A, sfixn *B, sfixn dA, sfixn dB,
                         sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void NegacyclicConvolution1(sfixn *res, sfixn s,
                             sfixn es1, sfixn es2,
                             sfixn K, sfixn dims2, sfixn ls2,
                             sfixn *A, sfixn *B, sfixn dA, sfixn dB, sfixn K2,
                             sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void WeightVectorAdaptiveEvaluation1(sfixn *res, sfixn es1, sfixn es2,
                                      sfixn K, sfixn K2,
                                      sfixn dims2, sfixn ls2,
                                      sfixn *A, sfixn dA,
                                      sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

  void AdaptiveInterpolationWeightVector1(sfixn *res, sfixn K, sfixn es1,
                                         sfixn es2, sfixn dims2, sfixn ls2,
                                         sfixn *rootsPtr,
                                         MONTP_OPT2_AS_GENE *pPtr);


  /*
      C = A*B mod p
      p: prime number
      N: number of variables
      C, A, B: dense recursive representation of C, A and B
      dc, da, db: partial degree vectors for C, A and B
      C should be initialized to zero
      default wp (whichprime) =0 to call the fft for a general Frouier prime
     */
    bool DMPMul(sfixn p,  int N,
                sfixn *C, int *dc,
                sfixn *A, int *da,
                sfixn *B, int *db,
                int wp);

/*
      C = A*B mod pPtr->P
      two variables
      C, A, B: dense recursive representation of C, A and B
      dc, da, db: partial degree vectors for C, A and B
      dc[0] = partial degree of variable 1
      pPtr: struct for Montgmery mod arithmetic
      C should be initialized to zero
     */
    void BivarMul(sfixn *C, int *dc,
                  sfixn *A, int *da,
                  sfixn *B, int *db,
                  MONTP_OPT2_AS_GENE *pPtr,
                  int wp);

/* 
       res should be initialized to zero
       s1 = partial degree of v1 of A + 1
       s2 = partial degree of v2 of A + 1
    */
    void Evaluation2D(sfixn *res, int es1, int es2,
                      int dims1, int dims2,
                      int ls1, int ls2,
                      sfixn *A, int s1, int s2,
                      sfixn *RT1, sfixn *RT2,
                      MONTP_OPT2_AS_GENE *pPtr,
                      int H, int *RevBidMap, int wp);

    void Interpolation2D(sfixn *res, int es1, int es2,
                         int dims1, int dims2,
                         int ls1, int ls2,
                         sfixn *invRT1, sfixn *invRT2,
                         MONTP_OPT2_AS_GENE *pPtr,
                         int H, int *RevBidMap, int wp);

	/**
	 * FFT-based multiplication
	 * For the particular prime 4179340454199820289
	 *
	 * @param l: log_2(s)
	 * @param s: size of the f and g (power of 2)
	 * @param f: an input vector
	 * @param g: another input vector
	 * @param p: the prime
	 * @param basezie: the base size (=1024)
	 **/
	bool ks_mul (int l, int s, sfixn* f, sfixn* g, sfixn p, int basesize);
} //PBPAS

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


