

/******************************
 * A collection of routines for hensel lifting of integer polynomials. 
 * NOTE: functions herein are still experimental, may not be bug-free, 
 *       and may change in future releases.
 ******************************/

#ifndef _SMZP_HENSEL_
#define _SMZP_HENSEL_

#include "SMZP_Support.h"
#include "DUZP_Support.h"
#include "DUZP_Hensel.h"


#ifdef __cplusplus
extern "C" {
#endif

/**
 * Convert the AltArrZ_t polynomial a, which should be univariate, 
 * to a duspoly_t over the prime field described by Pptr.
 * That is, returns a mod Pptr->prime.
 *
 * @param a the polynomial to convert
 * @param Pptr the Prime_ptr describing the finite field
 * 
 * @return the resulting modular polynomial.
 *
 */
duspoly_t* univarToPrimeField_AAZ(const AltArrZ_t* a, const Prime_ptr* Pptr);


/**
 * Efficiently expand the simple binomial (x - a) raised to the power of n.
 * For nvar > 1, x is assumed to be the least significant variable.
 *
 * @param a the constant term of the binomial.
 * @param n the power with which to raise the binomial.
 * @param nvar the number of variables in the final polynomial 
 *
 * @return the binomial power as an expanded polynomial.
 *
 */
AltArrZ_t* binomialExpansion(const mpz_t a, unsigned int n, int nvar);

/**
 * Lift the r univariate polynomials pointed to by fs to be bivariate polynomials
 * satisfying a = prod(liftedF) where a(x_2 = x) = prod(fs).
 * In this restricted format, a should be monic
 * and the lifting occurs w.r.t the second variable of a.
 * The lifted polynomials are returned in liftedF.
 *
 * @param a a bivariate polynomial.
 * @param fs the list of univariate image polynomials to lift.
 * @param r the number of polynomials in the list fs.
 * @param x the evaluation point for the second variable in a.
 * @param Pptr the prime modulo.
 * @param[out] liftedF the pre-allocated list where lifted images are stored.
 *
 * @return 1 iff the lift was successful
 */ 
int monicBivarHenselLift_AAZ(const AltArrZ_t* a, AltArrZ_t const*const* fs, unsigned int r, const mpz_t x, const Prime_ptr* Pptr, AltArrZ_t** liftedF);

/**
 * Lift the r univariate polynomials pointed to by fs to be bivariate polynomials
 * satisfying a = prod(liftedF) where a(x_2 = x) = prod(fs).
 * The lifted polynomials are returned in liftedF.
 * Both a and the polynomials in lcF should be in the same number of variables, 
 * but the partial degree of lcF in x_1 should be 0.
 *
 * @param a a bivariate polynomial.
 * @param lcF the true leading coefficients of factors to lift.
 * @param fs the list of univariate image polynomials to lift.
 * @param r the number of polynomials in the list fs.
 * @param x the evaluation point for the second variable in a.
 * @param Pptr the prime modulo.
 * @param[out] liftedF the pre-allocated list where lifted images are stored.
 * 
 * @return 1 iff the lift was successful
 */ 
int bivarHenselLift_AAZ(const AltArrZ_t* a, DUZP_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t x, const Prime_ptr* Pptr, AltArrZ_t** liftedF);


/**
 * Lift the r univariate polynomials pointed to by fs to be trivariate polynomials
 * satisfying a = prod(liftedF) where a(x_2 = x) = prod(fs).
 * The lifted polynomials are returned in liftedF.
 * Both a and the polynomials in lcF should be in the same number of variables, 
 * but the partial degree of lcF in x_1 should be 0.
 *
 * @param a A trivariate polynomial.
 * @param lcF the true leading coefficients of factors to lift.
 * @param fs the list of univariate image polynomials to lift.
 * @param r the number of polynomials in the list fs.
 * @param x the evaluation point for the second variable in a.
 * @param Pptr the prime modulo.
 * @param[out] liftedF the pre-allocated list where lifted images are stored.
 * 
 * @return 1 iff the lift was successful
 */ 
int trivarHenselLift_AAZ(const AltArrZ_t* a, AltArrZ_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t* x, const Prime_ptr* Pptr, AltArrZ_t** liftedF);


/**
 * Lift the r univariate polynomials pointed to by fs to be n-variable polynomials
 * satisfying a = prod(liftedF) where a(x_2=x[1],...x_n=x[n-1]) = prod(fs).
 * The lifted polynomials are returned in liftedF.
 * Both a and the polynomials in lcF should be in the same number of variables, 
 * but the partial degree of lcF in x_1 should be 0.
 *
 * @param a a multivariate polynomial with positive degree in at least two variables.
 * @param lcF the true leading coefficients of factors to lift.
 * @param fs the list of univariate image polynomials to lift.
 * @param r the number of polynomials in the list fs.
 * @param x the evaluation point for the other variables in a used to determine fs.
 * @param Pptr the prime modulo.
 * @param[out] liftedF the pre-allocated list where lifted images are stored.
 * 
 * @return 1 iff the lift was successful
 */ 
int henselLift_AAZ(const AltArrZ_t* a, AltArrZ_t** lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t* x, const Prime_ptr* Pptr, AltArrZ_t** liftedF);

/**
 * Multi-term bivariate diophantine solver. 
 * Solve the multi-term bivariate diophantine problem for sigma_i in 
 * b[0]*sigma[0] + ... + b[r-1]*sigma[r-1] = c mod p 
 * where b[j] = prod us[i] for i != j.
 * Each sigma_i returned in the pre-allocated array sigmas.
 *
 *
 */
int multiBDP_AAZ(const AltArrZ_t* c, AltArrZ_t const* const* us, AltArrZ_t** sigmas, unsigned int r, const Prime_ptr* Pptr, const mpz_t mpPrime);

/**
 * Multi-term multivariate diophantine solver. 
 * Solve the multi-term multivariate diophantine problem for sigma_i in 
 * b[0]*sigma[0] + ... + b[r-1]*sigma[r-1] = c mod p 
 * where b[j] = prod us[i] for i != j.
 * Each sigma_i returned in the pre-allocated array sigmas.
 * Number of variables of each us[i] and c must be the same and >= 3.
 * Problem is solved over the finite field described by Pptr.
 *
 */
int multiMDP_AAZ(const AltArrZ_t* c, AltArrZ_t const* const* us, AltArrZ_t** sigmas, unsigned int r, const Prime_ptr* Pptr, const mpz_t mpPrime);





//TODO move to right header file


/**
 * Incrementally interpolate the jth point (x[j],c) using the polynomial pointed to
 * by a_ptr as the existing interpolant. Returns a 1 on successful interpolation.
 * Failure may occur if the points being interpolated cannot be fitted by an integer
 * polynomial. 
 *
 *
 * @param a_ptr a pointer to the existing interpolating polynomial to be updated
 * @param x the list of interpolation nodes
 * @param j the index of the new point to interpolate within x. 
 * @param bj the value corresponding with x[j] forming the point to interpolate.  
 *
 * @return 1 if the interpolation was successful, 
 *         0 if the interpolation did not add degree to the interpolant,
 *        -1 if the interpolation was not successful.
 */
int univarIncrNewtonInterp_AAZ(AltArrZ_t** a_ptr, const mpz_t* x, int j, const mpz_t bj);

/**
 * Incrementally interpolate the jth point (x[j],c) using the polynomial pointed to
 * by a_ptr as the existing interpolant. Returns 1 on successful interpolation.
 * Failure may occur if the interpolation nodes are not unique.
 * The interpolation occurs modulo Pptr->prime.
 *
 *
 * @param a_ptr a pointer to the existing interpolating polynomial to be updated
 * @param x the list of interpolation nodes
 * @param j the index of the new point to interpolate within x. 
 * @param bj the value corresponding with x[j] forming the point to interpolate.  
 * @param Pptr the modulus 
 *
 * @return 1 if the interpolation was successful, 
 *         0 if the interpolation did not add degree to the interpolant,
 *        -1 if the interpolation was not successful.
 */
int univarIncrNewtonInterpModP_AAZ(AltArrZ_t** a_ptr, const mpz_t* x, int j, const mpz_t bj, const Prime_ptr* Pptr);

/** 
 * Evaluate all variables except the main variable of the input polynomial,
 * returning the result as a univariate DUZP polynomial.
 * 
 * @param aa the polynomial to evaluate
 * @param vals a list of aa->nvar-1 values at which to evaluate aa.
 *
 * @return the resulting evaluated polynomial in a dense representation DUZP.  
 */
DUZP_t* evaluatePolyToDUZP_AAZ(const AltArrZ_t* aa, const mpz_t* vals);

/**
 * Evaluate all variables except the main variable of the input polynomial,
 * returning the result as a univariate DUZP polynomial.
 * 
 * @param aa the polynomial to evaluate
 * @param vals a list of aa->nvar-1 values at which to evaluate aa.
 * @param res_p a pointer to the (possibly pre-allocated) duspoly_t to be returned after evaluation.
 * @param Pptr the Prime_ptr describing the finite field.
 *
 * @return the resulting evaluated polynomial in a dense representation DUZP.  
 */
void evaluatePolyToDUSP_AAZ(const AltArrZ_t* aa, const mpz_t* vals, duspoly_t** res_p, const Prime_ptr* Pptr);

/**
 * Combine many polynomials which represent the dense representation of
 * a polynomial with polynomial coefficients in a sparse polynomial in one more variable.
 * That is the index of the polynomial in coefs represents the degree of
 * a new variable corresponding to that coefficient. 
 *
 * The polynomials in coefs must all have the same number of variables.
 *
 * @param coefs the dense list of polynomial coefficients
 * @param n the number of coefficients in coefs.
 * 
 * @return a sparse polynomial in coefs[0]->nvar + 1 variables
 *         whose terms are built from coefs. 
 */
AltArrZ_t* combinePolyCoefs(AltArrZ_t const*const* coefs, int n);

/**
 * Combine many univariate polynomials which represent the dense representation of
 * a bivariate polynomial with.
 * That is, the index of the polynomial in coefs represents the degree of
 * a new variable corresponding to that coefficient. 
 *
 * The polynomials in coefs must all have the same number of variables.
 *
 * @param coefs the dense list of polynomial coefficients
 * @param n the number of coefficients in coefs.
 * 
 * @return a sparse polynomial in coefs[0]->nvar + 1 variables
 *         whose terms are built from coefs. 
 */
AltArrZ_t* combineUnivarPolyCoefs(AltArrZ_t const*const* coefs, int n);


/**
 * Replace the leading coefficient of p with lc, where p is viewed 
 * recursively as a univariate polynomial in its first variable.
 *
 * Here, p and lc should have the same number of variables.
 * Naturally, though, lc should have degree 0 in the first variable.
 *
 * @param p A pointer to the polynomial to modify.
 * @param lc the polynomial to set as p's leading coefficient.
 *
 */ 
void replaceLC_AAZ_inp(AltArrZ_t** p, const AltArrZ_t* lc);

#ifdef __cplusplus
}
#endif


#endif
