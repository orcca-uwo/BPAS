
// ser_basic_routine.h

#ifndef __ser_basic_routine_h
#define __ser_basic_routine_h


#include "modpn.h"
#include "modpn_export.h"
#include "ser_general_routine.h"

namespace SERBPAS {

  void MultiplyByTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void tftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiplyByFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);

  void fftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr);

  sfixn * MultiEvalByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr);

  void InterpolByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr);

  void MultiplyByKroneckerFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr);
  
  void EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void plainMultiDMul(sfixn N, sfixn * ccum, sfixn * res, sfixn * ccum1, sfixn * dgs1, sfixn * ccum2, sfixn * dgs2, sfixn * coeffs1, sfixn * coeffs2, MONTP_OPT2_AS_GENE * pPtr);

  void decomposePoly(sfixn N, sfixn * ccum, sfixn * res, sfixn N1, sfixn * dgs1, sfixn * coeffs1, sfixn * ccum1, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr, sfixn R, sfixn SFT);

  void decomposePoly2(sfixn N, sfixn * ccum, sfixn * res, sfixn num, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr);
  
  void UniFastMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, sfixn Lbound, sfixn * BRevInvPtr,  MONTP_OPT2_AS_GENE * pPtr);

  void 
	UniPlainMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );


  void EX_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void EX_Mont_FFTMul_OPT2_AS_GENE(sfixn n, sfixn e, sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

  void getRevInvTiSet(sfixn *dgs, sfixn N, TriRevInvSet * tris, TriSet * ts,  MONTP_OPT2_AS_GENE * pPtr);

  void NewtonRevInverse(sfixn N, preFFTRep * tRIPtrtmp, preFFTRep * tRIPtr, preFFTRep * tPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);


  sfixn *
	modularInvPM(sfixn degG, // degG=n-1;
				 sfixn * GPtr, sfixn degF, sfixn * FPtr, 
				 sfixn r, sfixn n,
				 MONTP_OPT2_AS_GENE * pPtr);
 
  void squarePoly_FFT(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep *p1Ptr,  MONTP_OPT2_AS_GENE * pPtr);
 
  void fftMultiD_square_test(sfixn * coeffs1, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);

  void EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr);

  void
	MultiMod(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void
	MultiMod_DF(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiMod_1(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiUniFastMod_1(sfixn N, preFFTRep * tmpIn,  sfixn n, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void EX_mulPoly(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);

  void EX_mulPoly_TFTFFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr);

  void EX_mulPoly_FFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);
  
  void mulPoly_FFT(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr);

  void polyMul_TFT(sfixn N, KroTFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr);

  void
	MultiMod_BULL(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  void MultiMod_1_BULL(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);


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
  

  
} //end of SERBPAS
#endif
