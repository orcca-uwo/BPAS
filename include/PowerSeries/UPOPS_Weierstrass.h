#ifndef _WEIERSTRASS_H_
#define _WEIERSTRASS_H_


#include "UnivariatePolynomialOverPowerSeries.h"
#include "RationalNumberPolynomial/SMQP_Support-AA.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Given power series F, G, H, with F = GH, compute the homogeneous part of G of degree r.
 * This assumes that F is known to at least degree r,
 * G has 0 as a constant term and is known to degree r-1,
 * and H is known to at least degree r-1.
 *
 * @param F, G, H: power series satisfying F = GH
 * @param r: the degree of the homogeneous part of G to compute
 *
 */
Poly_ptr lemmaForWeiestrass_UPOPS(PowerSeries_t* F, PowerSeries_t* G, PowerSeries_t* H, int r);

/**
 * Update all power series coefficients of p and alpha by one additional degree of precision,
 * given upops = p*alpha
 * @param upops : the input upops
 * @param p : the upops p
 * @param alpha : the upops alpha
 */
void weierstrassUpdate_UPOPS(Upops_t* p, Upops_t* alpha);


/**
 * Given an input upops, f, apply Weierstrass preparation to obtain f = p*alpha.
 * F must have at least one power series coefficient being a unit.
 * p and alpha are computed lazily.
 * p is of degree d, where d is the smallest integer such that the d'th coefficient
 * of f is a unit.
 *
 * @note It is invalid to destory p if there are any references to its coefficients
 * by other power series. E.g. when a product is computed for p*alpha.
 * If the returned p happens to be destroy destroyed, then
 * its power series coefficients will be truncated,
 * despite any reference to them from any other ancestry. Below explains.
 *
 * The returned p and alpha, and the power series coefficients contained therein,
 * are tightly coupled. Both are needed to updated a single one of them.
 * Enough information is captured internally to destroy the resulting alpha,
 * but the same is not true for p. Destroying the returned p
 * (without first reserving it or copying it) will result in its
 * power series coefficients being truncated (if referenced by another
 * power series ancestry) or free'd completely.
 * Indeed, since there is a circular reference between p and it's power series coefficeints
 * (the upops p is a generator parameter of the power series coefficients),
 * then the reference from the coefficients to the upops p is implemented as a
 * *weak reference*. Thus, p will be fully destroyed if no strong reference remains to it.
 * At that time, p's power series coefficients will be truncated or free'd.
 *
 * @param f : the input upops to factorize
 * @p_out : a pointer to p
 * @alpha_out : a pointer to alpha
 */
void weierstrassPreparation_UPOPS(Upops_t* f, Upops_t** p_out, Upops_t** alpha_out);


#ifdef __cplusplus
}
#endif

#endif
