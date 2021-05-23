


#ifndef _SMZP_FACTORING_SUPPORT_H_
#define _SMZP_FACTORING_SUPPORT_H_

#include "SMZP_Support.h"

#ifdef __cplusplus
extern "C" {
#endif

// /**
//  * Choose a "random" evaluation point
//  */
// void chooseRandomEvaluationPoint_factAAZ(mpz_t* values, AltArrZ_t** U0, mpz_t delta, const AltArrZ_t* nonzero, const AltArrZ_t* f_in, mpz_t bound);

/**
 * Choose a "random" evaluation point for Hensel.
 * This chooses nvar-1 values, each within bound, such that
 * the polynomial nonzero does not evaluate to 0 at that point.
 * The point is chosen such that the chosen values are pair-wise co-prime.
 *
 * @param nonzero a polynomial which should not evaluate to 0 under the chosen point
 * @param bound the bound for each component of the point
 * @param[out] values a vector of size nvar which stores to the point at indices 1..nvar-1
 * @param nvar the number of variables in nonzero
 */
void chooseRandomEvaluationPointHensel_AAZ(AltArrZ_t const*const* nonzero, int nnzero, mpz_t bound, mpz_t* values, int nvar);


/**
 * Algorithm N "Nondivisors" from Wang, 1978.
 * Used as part of leading coefficient reconsutrction.
 * Returns 1 iff the evaluation point is sufficient to reconstruct the leading coefficients.
 *
 * @param[out] d an array of size k to store "special integers" and as working space for the algorithm.
 * @param Omega the content of the leading coefficient of the polynomial to be factored
 * @param F_tilde an array of the evaluations of the factors of the l.c. at the evaluation point
 * @param k the number of l.c. factors
 * @param delta the content of the polynomial to be factored when evaluated at the evaluation point
 * @return 1 if the evaluation point should be successful for leading coefficient reconstruction
 */
int checkEvaluationPoint_factAAZ(mpz_t* d, const mpz_t Omega, const mpz_t* F_tilde, long k, const mpz_t delta);


/**
 * "Reconstruct"--obtain--the correct multivariate leading coefficients
 * of the factors (univariate, to be lifted) of f. Given a factorization
 * of the leading coefficient of f as lc_factors, distribute them and the integer content
 * of f across the univariate images of the factors of f in preparation for hensel lifitng.
 *
 * @param[out] lcs an array of size nufs to the true leading coefficients of the factors of f
 * @param f the polynomial being factored; must be non-const since it is modified as a last resort
 * @param lc_factors the factors of the leadng coefficient of f
 * @param lc_exps the exponents on the factors of the leading coefficient of f
 * @param nlcfs the number of leading coefficient factors
 * @param F_tilde an array of the evaluations of the factors of the l.c. at the evaluation point
 * @param U0_factors the univariate images of the factors of f
 * @param delta_in the content of the univariate image of f
 * @param values the evaluation point used to obtain the univariate images
 * @return 1 iff the true leading coefficients of the factors were determined.
 */
int reconstructLeadingCoefficients_factAAZ(AltArrZ_t** lcs, AltArrZ_t* f, AltArrZ_t const*const* lc_factors, const int* lc_exps, int nlcfs,
	const mpz_t* F_tilde, DUZP_t** U0_factors, int nufs, const mpz_t delta_in, const mpz_t* values);



#ifdef __cplusplus
}
#endif

#endif