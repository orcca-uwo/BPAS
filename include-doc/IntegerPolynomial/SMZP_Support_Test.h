
#ifndef _SMZP_SUPPORT_TEST_AA_H_
#define _SMZP_SUPPORT_TEST_AA_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "SMZP_Support.h"
#include <math.h>
#include <time.h>

/** 
 * Build a random polynomial given the the number of variables, nvar,
 * the number of terms, nterms, an (exclusive) upper bound on the absolute value
 * of the cofficients and a sparsity factor. 
 *
 * The sparsity fa-ctor is such that the difference in degree_t between sucsessive terms 
 * in the generated polynomial is 1 <= diff < sparsity;
 *
 */
Node* buildRandomZPoly(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg);

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
AltArrZ_t* buildRandomZPolyFromMax(int nvar, const int* maxDegs, unsigned long int coefBound, float sparsity, int includeNeg);


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


