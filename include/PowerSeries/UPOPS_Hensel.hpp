

#ifndef _UPOPS_HENSEL_H_
#define _UPOPS_HENSEL_H_

#include "UnivariatePolynomialOverPowerSeries.h"

/**
 * Given a UPOPS, evaluate power series variables to 0 and
 * factor the resulting univariate polynomial into monic linear factors.
 * Those factors are returned as an array of rational numbers c
 * and an array of exponents k, such that f = \prod_i (x-c_i)^{k_i}.
 *
 * @note It is assumed that f is monic and factors into linear factors over Q.
 *
 * @param f: the upops to evaluate and factor
 * @param c[out]: a pointer to an array of mpq_t.
 *                If the pointer points to NULL, an array is allocated, otherwise
 *                it is assumed that the array is of size at least f->deg + 1,
 *                and that each mpq_t is initialized.
 * @param k[out]: a pointer to an array of ints.
 *                If the pointer points to NULL, an array is allocated, otherwise
 *                it is assumed that the array is of size at least f->deg + 1.
 * @param n[out]: the number of factors is returned within.
 */
void evalToZeroAndFactor_UPOPS(Upops_t* f, mpq_t** c, int** k, int* n);



/**
 * Perform a Taylor shift (a translation) of f
 * and return the result constructed lazily.
 *
 * @param f : the input upops
 * @param c_r : the rational number for the translation
 * @return the translated upops
 */
Upops_t* taylorShift_UPOPS(Upops_t* f, const mpq_t c_r);


#if defined(WITH_NTL) && WITH_NTL
/**
 * Factorize upops f using Hensel's lemma.
 * Returns an array of size n of the upops factors of f
 * found using Hensel's lemma and Weierstrass preparation.
 * 
 * @Note: f is assumed to be monic in its polynomial variable.
 *
 * @Note: it is assumed that \bar{f} factorizes into linear
 * factors over Q, with \bar{f} being the polynomial resulting
 * from evaluating all power series coefficeints of f at the origin.
 *
 * @Note requires to built with NTL
 *
 * If either factrs or n is NULL, this function does nothing.
 *
 * @param f, the input upops
 * @param facts[out] a pointer to an array of size n of the upops factors of f,
                     if facts does not point to NULL, it is assumed to be
                     pre-allocated with a large enough allocation.
 * @param n[out], a pointer to the number of roots.
 */
void HenselFactorization_UPOPS(Upops_t* f, Upops_t*** facts, int* n);
#endif


/**
 * Factorize upops f using Hensel's lemma.
 * Returns an array of size n of the upops factors of f
 * found using Hensel's lemma and Weierstrass preparation.
 *
 * C is a list of roots of size n of \bar{f}, with
 * \bar{f} being the polynomial resulting from evaluating
 * all power series coefficeints of f at the origin.
 *
 * @Note: f is assumed ot be monic in its polynomial variable. 
 *
 * @Note: it is assumed that \bar{f} factorizes into linear
 * factors over Q whose roots are given by C.
 *
 * @param f, the input upops
 * @param C, the list of roots of size of of \bar{f}.
 * @param n, the number of roots.
 * @param facts[out] a pointer to an array of size n of the upops factors of f,
                     if facts does not point to NULL, it is assumed to be
                     pre-allocated with a large enough allocation.
 */
void HenselFactorization_UPOPS(Upops_t* f, mpq_t* C, int n, Upops_t*** facts);





#endif
