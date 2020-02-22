/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __generalFuncs_h
#define __generalFuncs_h 


#include "Types.h"
#include "AS.h"
#include "inlineFuncs.h"

extern sfixn BASE;


void *my_calloc(size_t size, size_t thesizeof);
void *my_malloc(size_t totsize);
void my_free(void *ptr);


void CDECL catch_intr( int32 signo );

// Create random vector with random coefficients.
sfixn * randomVec(sfixn *, sfixn, sfixn);

void
printPolyUni(int32 deg, sfixn * coeffs);

sfixn * randomMonicVec(sfixn * vec, sfixn s, sfixn p);

void qr_seprator(int32 *nloop, int32 *q, int32 *r, int32 size, int32 n);

void qr_seprator2(int32 *nloop, int32 *q, int32 *r, int32 size, int32 n, int32 co);

void printVec_1_to_n_double(sfixn n, double * data);

void copyVec_1_to_n(int32 n, sfixn * desV, sfixn * srcV);

void copyVec_0_to_d(int32 d, sfixn * desV, sfixn * srcV);

sfixn * allOneVec(sfixn *, sfixn);
// Create random vector with random coefficients.
sfixn * randomVecSeed(sfixn, sfixn *, sfixn, sfixn, unsigned);

// Print32 the polynomial.
//void printPoly(int32, sfixn *);

double gettime();

unsigned getSeed();

sfixn * reverse1(sfixn, sfixn *);

sfixn * reverseUni(sfixn, sfixn *, sfixn *);

sfixn logceiling(sfixn );

int32 compareVec(sfixn, sfixn *, sfixn *);

sfixn * reverseMulti_1(sfixn deg, sfixn sizOfCoef, sfixn * vec);

sfixn * reverseMulti(sfixn deg, sfixn sizOfCoef, sfixn * vec1, sfixn * vec2);

void printVec(sfixn deg, sfixn * data);

void fprintVec(FILE * F, sfixn deg, sfixn * data);


void printVecST(sfixn S, sfixn T, sfixn * data);

void printException(int32 menu);

void printVecLong(sfixn deg, longfixnum * data);

void printVecFrom1(sfixn n, sfixn * data);


void fprintVecDoubleDoubleBlock(FILE * F, sfixn Xst, sfixn Xed, sfixn Xinterval, sfixn Yst, sfixn Yed,sfixn Yinterval,  double ** data, sfixn deg, sfixn * vec);


void printVecDoubleDoubleBlock(sfixn Xst, sfixn Xed, sfixn Xinterval, sfixn Yst, sfixn Yed, sfixn Yinterval, double ** data, sfixn deg, sfixn * vec);


void fprintVecDoubleDoubleBlockGnuplot(FILE * F, sfixn Xst, sfixn Xed, sfixn Xinterval, sfixn Yst, sfixn Yed,sfixn Yinterval,  double ** data, sfixn deg, sfixn * vec);

void cleanVec(sfixn deg, sfixn * cof);

void cleanVecft(sfixn from, sfixn to, sfixn * cof);

void cleanVecINT(sfixn deg, int32 * cof);

void cleanVecDOUBLE(int32 deg, double * cof);

//sfixn inverseMod(sfixn n, sfixn PPP);

//void egcd (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo);

void aborting(const char * );


// 0 means NO, it is NOT a constatn vector.
int32 isConstVec(sfixn n, sfixn * data);

int32 isZeroVec(int32 siz, sfixn *vec);

void freeVecVec(sfixn m, sfixn **ptsPtr);

sfixn*  EX_copyVec_0_to_d(int32 d, sfixn * srcV);

void
EX_printPolyUni(int32 deg, sfixn * coeffs, char var);


void printVecAndIndex(sfixn deg, sfixn * data);

int32 isVecContainsZero(int32 siz, sfixn *vec);

sfixn
shiftBigger_1(int32 d, sfixn * vec, int32 m);

int32 isZeroPolyUni(sfixn d, sfixn *vec);

void
fprintPolyUni(FILE *file, int32 deg, sfixn * coeffs);

sfixn *EX_RandomUniPolyCoeffsVec(sfixn d, sfixn p);
sfixn * EX_RandomUniPolyCoeffsVecMonic(sfixn d, sfixn p);

#endif
