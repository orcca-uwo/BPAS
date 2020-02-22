
#ifndef _LIN_ALG_IML_H_
#define _LIN_ALG_IML_H_

#include <stdlib.h>
#include <gmp.h>
#include "../iml/include/iml.h"
#include "../iml/include/basisop.h"
#include "../iml/include//RNSop.h"
#include "../iml/include/certsolve.h"
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Attempt to solve a system of rational number linear equations.
 * A is an n x n matrix stored in row-major format.
 * b is an n x 1 matrix stored in row-major format.
 * n is the dimension of the system. 
 *
 * returns r, the rank of matrix A.
 * 
 * The solution vector is returned in X, if a solution can be found. 
 * X should be pre-allocated and initialized.
 *
 * rowList returns a list of indices whereby the first r indices represent good, 
 * linearly independent rows of A. The remaining n-r indices indicate rows that
 * are dependent and defective. 
 */
long solveMPQSystem(const mpq_t* A, const mpq_t* b, long n, mpq_t* X, long** rowList);

/**
 * Attempts to solve a system of integer linear equations with a rational number 
 * result.
 * A is an n x n matrix stored in row-major format.
 * b is an n x 1 matrix stored in row-major format.
 * n is the dimension of the system. 
 *
 * returns r, the rank of matrix A.
 * 
 * The solution vector is returned in mp_N, if a solution can be found.
 * A common denominator of the solution vector is returned in mp_D. 
 * mp_N should be pre-allocated and initialized. mp_D should be initialized prior.
 * 
 * rowList returns a list of indices whereby the first r indices represent good, 
 * linearly independent rows of A. The remaining n-r indices indicate rows that
 * are dependent and defective. 
 */
long solveMPZSystem(const mpz_t* mp_A, const mpz_t* mp_b, long n, mpz_t* mp_N, mpz_t mp_D, long** rowList);

/* Generate a n x m random dense matrix of type mpz_t with entries lying in
 * (-2^bd, 2^bd).
 * Returns the newly created matrix, stored in row-major order.
 */
mpz_t * randomMPZMat (const long n, const long m, const long bits, unsigned long seed);

/* Generate a n x m random dense matrix of type mpq_t with entries lying in
 * (-2^bd, 2^bd) 
 * Returns the newly created matrix, stored in row-major order.
 */
mpq_t * randomMPQMat (const long n, const long m, const long bits, unsigned long seed);

/** 
 * Generate an n x n ration number identity matrix.
 * Returns the newly created matrix, stored in row-major order. 
 */
mpq_t* getMPQIdentity(int n);

/**
 * Print an Integer matrix of dimensions n x m.
 */
void printMPZMat(const mpz_t* A_mp, int n, int m);

/**
 * Print a Rational number matrix of dimensions n x m.
 */
void printMPQMat(const mpq_t* A_mp, int n, int m);

/**
 * Free an n x m rational number matrix.
 */
void freeMPQMat(mpq_t* A, int n, int m);

/**
 * Free an n x m integer matrix.
 */
void freeMPZMat(mpz_t* A, int n, int m);




#ifdef __cplusplus
}
#endif

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


