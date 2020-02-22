#ifndef __generalFuncs_h
#define __generalFuncs_h 


#include "Types.h"
#include "AS.h"
#include "inlineFuncs.h"

extern sfixn BASE;


void *my_calloc(size_t size, size_t thesizeof);
void *my_malloc(size_t totsize);
void my_free(void *ptr);


void CDECL catch_intr( int signo );

// Create random vector with random coefficients.
sfixn * randomVec(sfixn *, sfixn, sfixn);

void
printPolyUni(int deg, sfixn * coeffs);

sfixn * randomMonicVec(sfixn * vec, sfixn s, sfixn p);

void qr_seprator(int *nloop, int *q, int *r, int size, int n);

void qr_seprator2(int *nloop, int *q, int *r, int size, int n, int co);

void printVec_1_to_n_double(sfixn n, double * data);

void copyVec_1_to_n(int n, sfixn * desV, sfixn * srcV);

void copyVec_0_to_d(int d, sfixn * desV, sfixn * srcV);

sfixn * allOneVec(sfixn *, sfixn);
// Create random vector with random coefficients.
sfixn * randomVecSeed(sfixn, sfixn *, sfixn, sfixn, unsigned);

// Print the polynomial.
//void printPoly(int, sfixn *);

double gettime();

unsigned getSeed();

sfixn * reverse1(sfixn, sfixn *);

sfixn * reverseUni(sfixn, sfixn *, sfixn *);

sfixn logceiling(sfixn );

int compareVec(sfixn, sfixn *, sfixn *);

sfixn * reverseMulti_1(sfixn deg, sfixn sizOfCoef, sfixn * vec);

sfixn * reverseMulti(sfixn deg, sfixn sizOfCoef, sfixn * vec1, sfixn * vec2);

void printVec(sfixn deg, sfixn * data);

void fprintVec(FILE * F, sfixn deg, sfixn * data);


void printVecST(sfixn S, sfixn T, sfixn * data);

void printException(int menu);

void printVecLong(sfixn deg, longfixnum * data);

void printVecFrom1(sfixn n, sfixn * data);


void fprintVecDoubleDoubleBlock(FILE * F, sfixn Xst, sfixn Xed, sfixn Xinterval, sfixn Yst, sfixn Yed,sfixn Yinterval,  double ** data, sfixn deg, sfixn * vec);


void printVecDoubleDoubleBlock(sfixn Xst, sfixn Xed, sfixn Xinterval, sfixn Yst, sfixn Yed, sfixn Yinterval, double ** data, sfixn deg, sfixn * vec);


void fprintVecDoubleDoubleBlockGnuplot(FILE * F, sfixn Xst, sfixn Xed, sfixn Xinterval, sfixn Yst, sfixn Yed,sfixn Yinterval,  double ** data, sfixn deg, sfixn * vec);

void cleanVec(sfixn deg, sfixn * cof);

void cleanVecft(sfixn from, sfixn to, sfixn * cof);

void cleanVecINT(sfixn deg, int * cof);

void cleanVecDOUBLE(int deg, double * cof);

//sfixn inverseMod(sfixn n, sfixn PPP);

//void egcd (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo);

void aborting(const char * );


// 0 means NO, it is NOT a constatn vector.
int isConstVec(sfixn n, sfixn * data);

int isZeroVec(int siz, sfixn *vec);

void freeVecVec(sfixn m, sfixn **ptsPtr);

sfixn*  EX_copyVec_0_to_d(int d, sfixn * srcV);

void
EX_printPolyUni(int deg, sfixn * coeffs, char var);


void printVecAndIndex(sfixn deg, sfixn * data);

int isVecContainsZero(int siz, sfixn *vec);

sfixn
shiftBigger_1(int d, sfixn * vec, int m);

int isZeroPolyUni(sfixn d, sfixn *vec);

void
fprintPolyUni(FILE *file, int deg, sfixn * coeffs);

sfixn *EX_RandomUniPolyCoeffsVec(sfixn d, sfixn p);
sfixn * EX_RandomUniPolyCoeffsVecMonic(sfixn d, sfixn p);

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


