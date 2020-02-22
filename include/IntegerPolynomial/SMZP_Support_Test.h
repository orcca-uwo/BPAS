
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
AltArrZ_t* buildRandomPoly_AAZ_unpk(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg);
AltArrZ_t* buildRandomSeededPoly_AAZ_unpk(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg, time_t seed);

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
