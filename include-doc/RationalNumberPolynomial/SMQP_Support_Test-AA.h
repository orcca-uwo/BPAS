
#ifndef _SMQP_SUPPORT_TEST_AA_H_
#define _SMQP_SUPPORT_TEST_AA_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "SMQP_Support-AA.h"
#include <math.h>
#include <time.h>

/** 
 * Get the next degree_t for a polynomial given the previous term's degree_t
 * and a "step" to take. In the univariate case this is prev+step. 
 * In the multivatirate case we consider an integer of radix maxUniDeg with 
 * coefficients described our degrees_t. We step such that the returned value
 * is equal to prev + step in this radix maxUniDeg representation.
 * e.g: prev = [1,2,7], step = 5, maxUniDeg = 10. Then next is [1,3,2];
 */
degrees_t getNextDegrees (degrees_t prev, degree_t step, degree_t maxUniDeg, int nvar, int* sizes, unsigned long long int* masks);

/** 
 * Build a random polynomial given the the number of variables, nvar,
 * the number of terms, nterms, an (exclusive) upper bound on the absolute value
 * of the cofficients and a sparsity factor. 
 *
 * The sparsity factor is such that the difference in degree_t between sucsessive terms 
 * in the generated polynomial is 1 <= diff < sparsity;
 *
 */
Node* buildRandomPoly(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg);

/**
 * Build a random polynomial given the number of variables, nvar, 
 * and the maximum degree of each variable, as the maxDegs arrays. 
 *
 * coefBound is the maximum number of bits in the coefficients.
 *
 * if includeNeg == 0 then all cofficients will be positive, otherwise randomly
 * negative.
 *
 * sparsity is a percentage of zero terms between the term with monomial of maxDegs
 * and the constant term. A sparsity of 0 produces a dense polynomial, a sparisty of 1
 * produces a polynomial of only one term, the one whose monomial is maxDegs.
 *
 * returns the randomly generated polynomial.
 */
AltArr_t* buildRandomPolyFromMax(int nvar, const int* maxDegs, unsigned long int coefBound, float sparsity, int includeNeg);
AltArr_t* buildRandomPolyFromMax_seeded(time_t seed, int nvar, const int* maxDegs, unsigned long int coefBound, float sparsity, int includeNeg);


/****************
* Random Polynomial Algorithm with Maximum Degree
*****************/

/**
 * @param[in] lastDegs Previous degree
 * @param[in] step The distance between the previous and returned degree
 * @param[in] nvar The number of variables
 * \example if the lastDegs = (2 2 4) and step = 3 then output is (1 2 4), and
 *  the sequence of in/out-puts is -> (0 2 4) -> (0 0 4) -> (0 0 3) -> (0 0 2).
 */
degrees_t getNextMaxDegs(degrees_t lastDegs, degrees_t step, int nvar, int* sizes, unsigned long long int* masks);
  
/**
 * Random Polynomial Algorithm with MaxDegs
 * @param[in] nvar The number of variables
 * @param[in] maxDegs The Maximum degrees of the output polynomial
 * @param[in] sparsity
 * @param[in] coefBound
 * @param[in] includeNeg
 * \brief This algorithm generates a random polynomial where the maximum degree is maxDegs.
 */
Node* randomMaxDegsPoly(int nvar, degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg);

/***************
* Random Triangular Set
***************/

/**
 * Convert the exponents of a polynomial
 * @param[in] inPoly The polynomial
 * @param[in] inNvar The number of variables of inPoly
 * @param[in] outNvar The number of variables of output polynomial.
 */
//void convertExponentPoly(Node* inPoly, int inNvar, int outNvar);

/** 
 * Build a Random Triangular Set based on maxDegs
 * @param[in] T A pointer to the triangular set
 * @param[in] outNvar The number of variables of output polynomial
 * @param[in] maxDegs The Maximum degree of polynomials
 * @param[in] sparsity
 * @param[in] coefBound
 * @param[in] includeNeg
 * \brief Output is a Lazard or general triangular set (w.r.t the lazard_flag)
 * such that if lazard_flag != 0, then the output set is Lazard triangular set.
 */
void randomTriangularSet (Node** T, int outNvar,degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg, int lazard_flag);
void randomTriangularSet_AA (AltArr_t** T, int outNvar,degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg, int lazard_flag);



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


