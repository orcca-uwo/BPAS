/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __MPMMTS_h
#define __MPMMTS_h 



#include "Types.h"
#include "MultiDFFT.h"
#include "generalFuncs.h"
#include "GCD.h"
#include "HGCD.h"
#include "FDIV.h"
#include "MapleCConverter.h"
#include <stdlib.h>            /* rand() and srand() functions               */
#include <assert.h>            /* assert()                                   */




extern char letters[];


preFFTRep * GcdAsUni(preFFTRep *f1, preFFTRep *f2, MONTP_OPT2_AS_GENE * pPtr);

int IsAllNumberCoeffs(preFFTRep *poly);

preFFTRep * QuoAsUni(preFFTRep *f1, preFFTRep *f2, MONTP_OPT2_AS_GENE * pPtr);

void getRevInvTiSet(sfixn *dgs, sfixn N, TriRevInvSet * tris, TriSet * ts,  MONTP_OPT2_AS_GENE * pPtr);

preFFTRep * 
EX_mul_Coef_Reduced(sfixn N, preFFTRep * out, preFFTRep * f1, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

sfixn 
MultiRecip(sfixn N, preFFTRep * invPtr, preFFTRep * cPtr, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

void 
copyPolyPointers(preFFTRep *D, preFFTRep *S);


void printPolyStrut(preFFTRep * Ptr);


void initRandomDenseTriSet( sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);


void initRandomTriSet( sfixn N, sfixn dgbound, TriSet * tPtr, MONTP_OPT2_AS_GENE * pPtr);

void InitOnePoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs);


void CopyOnePoly(preFFTRep * rPtr, preFFTRep * fPtr );

void freeTriSet(TriSet * tPtr);

void InitOneRevInv(sfixn N, preFFTRep * tRIPtr, sfixn * bounds, sfixn di);

void initTriRevInvSet(sfixn *dgs, sfixn N,  TriRevInvSet * tRevInvPtr, TriSet * tPtr);

void initCopyOneRevInv(preFFTRep * outPtr, preFFTRep * inPtr);

void freeTriRevInvSet(TriRevInvSet * tRIPtr);


void printPoly(preFFTRep * Ptr);

void printTriSet(TriSet * tPtr);

void printTriRevInvSet(TriRevInvSet * triPtr);

void freePoly(preFFTRep * x);

void freeKroFFTRep(KroFFTRep * x);

void InitOneReducedPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs);

void
InitOneMonicRandomPolys(sfixn * , sfixn, preFFTRep *, MONTP_OPT2_AS_GENE * , sfixn);

void initRandomTriSet_deg( sfixn N, sfixn dgbound, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);
void randomPoly(preFFTRep * Ptr, sfixn p);

void randomSparsePoly(preFFTRep * Ptr, sfixn p);


void PolyCleanData(preFFTRep * prePtr);


void randomPoly_allOne(preFFTRep * Ptr, sfixn p);

void randomMonicPoly(preFFTRep * Ptr, sfixn p);

void randomMonicPoly_allOne(preFFTRep * Ptr, sfixn p);


void InitKroFFTRep(KroFFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,  MONTP_OPT2_AS_GENE * pPtr);

void InitResPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs, sfixn * p2dgs);

void
InitOneRandomReducedInputPoly(sfixn * bounds, sfixn N, preFFTRep * p1Ptr, MONTP_OPT2_AS_GENE * pPtr);


void initRandomTriSet_dgs_allOne( sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);


void InitTwoRandomReducedInputPolys(sfixn * dgs, sfixn N, preFFTRep * p1Ptr, preFFTRep * p2Ptr, MONTP_OPT2_AS_GENE * pPtr);

void MultiD_KroneckerMul(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr, MONTP_OPT2_AS_GENE * pPtr);

void MultiD_KroneckerSquare(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep * 
EX_mul_Reduced(sfixn N, preFFTRep * out, preFFTRep * f1, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

sfixn
onePolyp(preFFTRep * polyPtr);

void fprintPoly(FILE *file, preFFTRep * Ptr);

void MultiMod(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

sfixn
MultiUniEuclidean(sfixn N, preFFTRep * gcd, preFFTRep * FPtr, preFFTRep * GPtr, 
		  TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

void reduceCoeffs(sfixn N, preFFTRep * toPtr, preFFTRep * fromPtr, TriSet * ts, TriRevInvSet * tris,   MONTP_OPT2_AS_GENE * pPtr);

int
comparePolyData(preFFTRep * p1, preFFTRep * p2);


void addPoly(sfixn N, preFFTRep * Ptr,  preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);

void subPoly(sfixn N, preFFTRep * Ptr, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);


preFFTRep *
addMulPoly_1(preFFTRep *f, preFFTRep *f1, preFFTRep *f2, sfixn N, 
             TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep *
subMulPoly_1(preFFTRep *f, preFFTRep *f1, preFFTRep *f2, sfixn N, 
             TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep *
powPoly_1(preFFTRep *f, sfixn e, sfixn N, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

void initRandomDenseTriSetForLifting(sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);

void LiftTinTriSet( sfixn N, TriSet * tPtr, int i, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
scalarMulJMatrix_1(sfixn r, POLYMATRIX * mat, MONTP_OPT2_AS_GENE * pPtr);

void
MultiMod_BULL(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

void * BlockBottom_UniDiv_1(void * Ptr);


void * Middle_UniDiv_BULL_PARA1 (void * pptr);

void
MultiMod_BULL_PARA1(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

  
void * reduceCoeffs_Inner_BULL_Para1(void * Ptr);

void init_example_1_DenseTriSetForLifting_y0_10( sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);

preFFTRep * 
EX_mul_Reduced_1(sfixn N, preFFTRep * f1, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

void EX_mulPoly_FFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr);


sfixn *EX_getPartialDegVec(preFFTRep *Ptr);


void randomSparsePoly(preFFTRep * Ptr, sfixn p);

int  EX_getPartialDegVec_inner(sfixn *coefVecSizAddr, sfixn N, sfixn *dgs, sfixn *accum, sfixn *data, sfixn *locAddr, sfixn *pdegVec);


sfixn *RemoveMinusTwo(sfixn *vec);

void RemoveLCMinusOneAndFillCoefs(sfixn base, sfixn *coefSizeFinalAddr, preFFTRep *Ptr, sfixn *coefVec, sfixn *locAddr, sfixn level,  sfixn *vec);

void EX_mulPoly_TFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr);


void EX_mulPoly_FFT_select(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr, int switcher);

void EX_mulPoly_TFTFFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr);

void EX_mulPoly_TFTFFT_Bench(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr, int fftNOTtft);

void initDenseTriSet( sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr);


void
InitOneMonicPolys(sfixn * bounds, sfixn N, preFFTRep * p1Ptr,  MONTP_OPT2_AS_GENE * pPtr, sfixn seed);

preFFTRep *EX_CopyOnePoly(preFFTRep * fPtr );

void EX_FreeOnePoly(preFFTRep * rPtr);

TriSet *EX_CopyOneTriSet(TriSet * srcts);

void EX_freeTriSet(TriSet * tPtr);

TriSet *
Ex_InitDenseTriSet( sfixn N, sfixn * dgs, MONTP_OPT2_AS_GENE * pPtr);


sfixn *bounds2dgs(TriSet *ts);

TriRevInvSet *EX_initTriRevInvSet(sfixn *dgs, sfixn N, TriSet * tPtr);

void EX_freeTriRevInvSet(TriRevInvSet * tRIPtr);

preFFTRep *EX_InitOnePoly(sfixn N, sfixn * dgs);

preFFTRep *EX_randomPoly(sfixn N, sfixn * dgs, sfixn p);

int
MonicizePoly_1(sfixn N, preFFTRep *inPoly, TriSet *ts, 
             TriRevInvSet *tris, MONTP_OPT2_AS_GENE *pPtr);

TriSet *
EX_initRandomNonMonicTriSet(sfixn N, sfixn dgbound, MONTP_OPT2_AS_GENE * pPtr);

int
MonicizeTriSet_1(sfixn N, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);


sfixn 
MultiRecip_Lift(sfixn N, preFFTRep * invPtr, preFFTRep * cPtr, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

void EX_mulPoly(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);

void
MultiMod_ForLifting(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr, sfixn size);

preFFTRep * 
EX_mul_Reduced_ForLifting(sfixn N, preFFTRep * out, preFFTRep * f1, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep * 
EX_ShiftPoly(preFFTRep * fPtr, int m);

preFFTRep *
direvativeMulti(preFFTRep * poly,  MONTP_OPT2_AS_GENE * pPtr);

void NormalizeTriSetBDS(sfixn N, TriSet * ts);


preFFTRep * 
MultiCoefPolyMul_1(sfixn N, preFFTRep * co, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr);

void
MultiMod_OPT(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr, int opt);


void
MultiMod_DF(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

TriSet *
EX_initRandomTriSet(sfixn N, sfixn dgbound, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep *
EX_MonicMultiPlainDivide(sfixn N, preFFTRep *FPtr, preFFTRep *GPtr, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);


preFFTRep *
EX_randomMonicPoly(sfixn N, sfixn * dgs, sfixn p);

preFFTRep * 
Ex_ReduceCoeffs(sfixn N, preFFTRep *fromPtr, TriSet *ts,  TriRevInvSet *tris, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep * 
EX_ReduceCoeffs(sfixn N, preFFTRep *fromPtr, TriSet *ts,  TriRevInvSet *tris, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep **EX_CopyOnePolyList(sfixn no, preFFTRep **PolyList );

preFFTRep *EX_getInitial(preFFTRep *poly);

preFFTRep *shrinkOneDim(preFFTRep *inPoly);

preFFTRep *EX_NormalizePoly(preFFTRep *inPoly);

preFFTRep *EX_NormalForm(sfixn N, preFFTRep *poly, TriSet *ts, TriRevInvSet *tris, MONTP_OPT2_AS_GENE *pPtr);

sfixn *EX_getDgsForNormalForm(preFFTRep *poly, TriSet *ts, sfixn e);

preFFTRep *EX_EY_ForNormalForm(preFFTRep *poly, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);




int
EX_IsEqualPoly(preFFTRep * Ptr1, preFFTRep * Ptr2);

preFFTRep *
EX_EY_Normalize(preFFTRep *poly, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);

sfixn realDegPoly(preFFTRep *poly, sfixn M);

preFFTRep *CreateZeroPoly();

preFFTRep *CreateConsPoly(sfixn cons);

preFFTRep *CreateUniPoly(sfixn dg, sfixn *vec);

TriSet * EX_ExchangeOnePoly(preFFTRep *poly, TriSet *ints, MONTP_OPT2_AS_GENE *pPtr);

TriSet * EX_MergeTriSet(sfixn start, sfixn index, preFFTRep *poly, TriSet *ts_top, TriSet *ts_under, MONTP_OPT2_AS_GENE *pPtr);

TriSet *EX_getLowerTriSet(sfixn M, TriSet * srcts);

TriRevInvSet *
EX_getRevInvTriSet(sfixn N,  TriSet *ts,  MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *
EX_InitResPoly(sfixn N, sfixn *p1dgs, sfixn *p2dgs);

preFFTRep *
EX_EX_mulPoly(sfixn N,  preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *
EX_EX_PlainMulPoly(sfixn N,  preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);


preFFTRep *
EX_EX_TFTFFTMulPoly(sfixn N,  preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *EX_GetPolyTail(preFFTRep * inPoly);


#endif
