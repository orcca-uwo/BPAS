
#ifndef _DUZP_HENSEL_
#define _DUZP_HENSEL_

#ifdef __cplusplus
extern "C" {
#endif


#include "IntegerPolynomial/SMZP_Support.h"
#include "IntegerPolynomial/DUZP_Support.h"
#include "ModularPolynomial/DUSP_Support.h"
#include "ModularPolynomial/DUSP_Support_Test.h"
#include "Utils/Unix_Timer.h"




////////////////////////
// Diophantine Solvers
////////////////////////

/**
 * Solve the two-term diophantine equation u*sigma + w*tau = c
 * for sigma nad tau. 
 *
 * sigma and tau are returned by pointer.
 *
 * @param u a polynomial over a finite field
 * @param w a polynomial over a finite field
 * @param c a polynomial over a finite field, the rhs of the diopahntine equation
 * @param[out] sigma a polynomial over the finite field satisfying the diophantine equation
 * @param[out] tau a polynomial over the finite field satisfying the diophantine equation
 * @param Pptr the prime of the finite field
 *
 *
 * @return 1 iff the equation can be solved (i.e. gcd(u,w) does not divide c )
 */
int UDP_spX(const duspoly_t* u, const duspoly_t* w, const duspoly_t* c, duspoly_t** sigma, duspoly_t** tau, const Prime_ptr* Pptr);


/**
 * Solve the multi-term univariate diophantine problem for sigma_i in 
 * b[0]*sigma[0] + ... + b[r-1]*sigma[r-1] = c mod p 
 * where b[j] = prod f[i] for i != j.
 * Each sigma_i returned in the pre-allocated array sigmas.
 *
 * @param c the rhs of the diophantine equation
 * @param fs an array of the factors of the b polynomials in the diophantine equation
 * @param r the size of the fs array
 * @param[out] sigmas a pre-allocated of size r to store the sigmas satisfying the equation
 * @param Pptr the prime of the finite field
 * 
 * @return 1 iff the problem could be solved (i.e. if gcd(f[i],f[j]) == 1 for all i,j
 * and p does not divide l.c. of f[i] for all i).
 */
int multiUDP_spX(const duspoly_t* c, duspoly_t const*const* fs, int r, duspoly_t** sigmas, const Prime_ptr* Pptr);





////////////////////////
// Two Term Lifts
////////////////////////

/**
 * Lift the factorization congruence a = uw mod p to obtain u and w over the integers.
 * This function is specialized to the case where a is monic and all coefficients of
 * it and its possible factors are positive.
 *
 * @param a a polynomial over the integers whose factors are to be lifted.
 * @param u a polynomial over the finite field which is a factor of a modulo the prime.
 * @param w a polynomial over the finite field which is a factor of a modulo the prime.
 * @param[out] liftedU the lifted factor corresponding to u.
 * @param[out] liftedW the lifted factor corresponding to w.
 * @param Pptr the prime of the finite field
 * 
 */
int monicPositiveTwoTermPadicLift(const DUZP_t* a, const duspoly_t* u, const duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const Prime_ptr* Pptr);


/**
 * Lift the factorization congruence a = uw mod p to obtain u and w over the integers.
 * This function is specialized to the case where a is monic.
 *
 * @param a a polynomial over the integers whose factors are to be lifted.
 * @param u a polynomial over the finite field which is a factor of a modulo the prime.
 * @param w a polynomial over the finite field which is a factor of a modulo the prime.
 * @param[out] liftedU the lifted factor corresponding to u.
 * @param[out] liftedW the lifted factor corresponding to w.
 * @param Pptr the prime of the finite field
 * 
 */
int monicTwoTermPadicLift(const DUZP_t* a, const duspoly_t* u, const duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const Prime_ptr* Pptr);


/**
 * Lift the factorization congruence a = uw mod p to obtain u and w over the integers.
 * This function is generic to handle when a is not monic and with possibly negative 
 * coefficients.
 *
 * @note the polynomials u and w are modified for efficiency.
 * 
 * @param a a primtive polynomial over the integers whose factors are to be lifted.
 * @param u a polynomial over the finite field which is a factor of a modulo the prime.
 * @param w a polynomial over the finite field which is a factor of a modulo the prime.
 * @param[out] liftedU the lifted factor corresponding to u.
 * @param[out] liftedW the lifted factor corresponding to w.
 * @param gamma a number is which known to be a multiple of the leading coefficient of u after lifting.
 * @param Pptr the prime of the finite field
 * 
 */
int twoTermPadicLift(const DUZP_t* aa, duspoly_t* u, duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const mpz_t gamma, const Prime_ptr* Pptr);


/**
 * Lift the factorization congruence a = uw mod p to obtain u and w over the integers.
 * This function is specialized to handle when a is not monic but with only positive coefficeints
 * in it and any of its possible factors.
 *
 * @note the polynomials u and w are modified for efficiency.
 *
 * @param a a primtive polynomial over the integers whose factors are to be lifted.
 * @param u a polynomial over the finite field which is a factor of a modulo the prime.
 * @param w a polynomial over the finite field which is a factor of a modulo the prime.
 * @param[out] liftedU the lifted factor corresponding to u.
 * @param[out] liftedW the lifted factor corresponding to w.
 * @param gamma a number is which known to be a multiple of the leading coefficient of u after lifting.
 * @param Pptr the prime of the finite field
 * 
 */
int positiveTwoTermPadicLift(const DUZP_t* aa, duspoly_t* u, duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const mpz_t gamma, const Prime_ptr* Pptr);





////////////////////////
// Multi Term Lifts
////////////////////////

/**
 * Given a monic polynomial, a, a = f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)liftedF[0]*...liftedF[nf-1] over Z.
 * This process is iterative in that it solves the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 *
 * @param a a monic polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials ove ra finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param nf the number of factors in f to be lifted.
 * @param Pptr the prime of the finite field.
 * 
 * @return the exponent of the modulus where lifting terminated. That is, returns e where a = lc(a)*liftedF[0]*...liftedF[nf-1] mod Pptr->prime^e.
 */ 
int monicMultiTermPadicLift_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr);


/**
 * Given a monic polynomial, a, a = f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)liftedF[0]*...liftedF[nf-1] over Z.
 * This process is optimized and iterative in that it solves the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 *
 * @param a a monic polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials ove ra finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param nf the number of factors in f to be lifted.
 * @param Pptr the prime of the finite field.
 * 
 */ 
int monicMultiTermPadicLiftOpt_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr);


/**
 * Given a primtive polynomial, a, a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)*liftedF[0]*...liftedF[nf-1] up to some much higher modulus.
 * That modulus is a power of Pptr->prime, where that power is the return of this function.
 * This process is optimized and iterative in that it solves the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 *
 * @param a a polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials over a finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param[out] sigmas the nf coefficients which solve the diophantine equation involved in the last hensel step.
 *             These coefficients are needed to pass to multiTermQuadraticPadicLiftResume in order to resume lifting.
 * @param nf the number of factors in f to be lifted; also the size of the arrays liftedF and sigmas.
 * @param Pptr the prime of the finite field.
 * 
 * @note The array sigmas is only filled if the returned exponent is greater than 1.
 * 
 * @return the exponent of the modulus where lifting terminated. That is, returns e where a = lc(a)*liftedF[0]*...liftedF[nf-1] mod Pptr->prime^e.
 * 
 */ 
int multiTermPadicLiftStart(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, duspoly_t*** sigmas, short nf, const mpz_t bound, const Prime_ptr* Pptr);


/**
 * Given a primtive polynomial, a, a = lc(a)*liftedF[0]*liftedF[1]*...*liftedF[nf-1] mod (Pptr->prime)^e, 
 * continue the lifting process so that a = lc(a)*liftedF[0]*liftedF[1]*...*liftedF[nf-1] mod (Pptr->prime)^m
 * where m is the smallest exponent such that Pptr->prime^m is greater than targetBound.
 * If the factors given in liftedF are not already lifted to Ppre->prime^e then no work is 
 * performed and -1 is returned to indicate an error. Otherwise, m is returned.
 *
 *
 * @param a a polynomial over the integers whose factors are to be lifted.
 * @param[in,out] liftedF an array of nf co-prime integer polynomials which are the factors of a modulo Pptr->prime^e. 
 *                These polynomials are modified in place so that on return they are factors of a 
 *                modulo Pptr->prime^m.
 * @param nf the number (> 1) of polynomials in the array liftedF.
 * @param[in,out] sigmas the coefficients solving the diophantine equation arising in the previous hensel lift.
 *                These polynomials are modified in place so that they may be used to resume lifting again.
 * @param Pptr the prime used in the modular equivalence.
 * @param e the exponent of the prime used in the modular equivalence of factors on input.
 * @param targetBound the minimum modulus to which factors are lifted. 
 * 
 * @return m, the smalled exponent such that Pptr->prime^m is greater than targetBound,
 *         or -1 if input polynomials do not satisfy input conditions. 
 *                
 */
int multiTermPadicLiftResume(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, duspoly_t const*const* sigmas, short nf, int e, const mpz_t bound, const Prime_ptr* Pptr);

/**
 * Given a primtive polynomial, a, a = lc(a)f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)*liftedF[0]*...liftedF[nf-1] over Z.
 * This process is optimized and iterative in that is solved the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 *
 * @note the polynomials in f are temporarily modified during this function, but returned to their original state before returning.
 * 
 * @param a a polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials ove ra finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param nf the number of factors in f to be lifted.
 * @param Pptr the prime of the finite field.
 * 
 * @return the exponent of the modulus where lifting terminated. That is, returns e where a = lc(a)*liftedF[0]*...liftedF[nf-1] mod Pptr->prime^e.
 */ 
static inline int multiTermPadicLiftOpt_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);

	int ret = multiTermPadicLiftStart(a, f, liftedF, NULL, nf, bound, Pptr);

	mpz_clear(bound);
	return ret;
}


/**
 * Given a monic polynomial, a, a = f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = liftedF[0]*...liftedF[nf-1] over Z.
 * This process is optimized and iterative in that it solves the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 * The lifting here is done quadratically, that is, each step lifts the factors to the modulus squared.
 *
 * @param a a monic polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials ove ra finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param nf the number of factors in f to be lifted.
 * @param Pptr the prime of the finite field.
 * 
 * @return the exponent of the modulus where lifting terminated. That is, returns e where a = lc(a)*liftedF[0]*...liftedF[nf-1] mod Pptr->prime^e.
 */ 
int monicMultiTermQuadraticPadicLiftOpt_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr);


/**
 * Given a primtive polynomial, a, a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)*liftedF[0]*...liftedF[nf-1] up to some much higher modulus.
 * That modulus is a power of Pptr->prime, where that power is the return of this function.
 * This process is optimized and iterative in that it solves the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 * The lifting here is done quadratically, that is, each step lifts the factors to the modulus squared.
 *
 * @param a a polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials over a finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param[out] sigmas the nf coefficients which solve the diophantine equation involved in the last hensel step.
 *             These coefficients are needed to pass to multiTermQuadraticPadicLiftResume in order to resume lifting.
 * @param nf the number of factors in f to be lifted; also the size of the arrays liftedF and sigmas.
 * @param Pptr the prime of the finite field.
 * 
 * @return the exponent of the modulus where lifting terminated. That is, returns e where a = lc(a)*liftedF[0]*...liftedF[nf-1] mod Pptr->prime^e.
 * 
 */ 
int multiTermQuadraticPadicLiftStart(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, DUZP_t*** sigmas, short nf, const mpz_t bound, const Prime_ptr* Pptr);


/**
 * Given a primtive polynomial, a, a = lc(a)*liftedF[0]*liftedF[1]*...*liftedF[nf-1] mod (Pptr->prime)^e, 
 * continue the lifting process so that a = lc(a)*liftedF[0]*liftedF[1]*...*liftedF[nf-1] mod (Pptr->prime)^m
 * where m is the smallest exponent such that Pptr->prime^m is greater than targetBound.
 * If the factors given in liftedF are not already lifted to Ppre->prime^e then no work is 
 * performed and -1 is returned to indicate an error. Otherwise, m is returned.
 *
 *
 * @param a a polynomial over the integers whose factors are to be lifted.
 * @param[in,out] liftedF an array of nf co-prime integer polynomials which are the factors of a modulo Pptr->prime^e. 
 *                These polynomials are modified in place so that on return they are factors of a 
 *                modulo Pptr->prime^m.
 * @param nf the number (> 1) of polynomials in the array liftedF.
 * @param[in,out] sigmas the coefficients solving the diophantine equation arising in the previous hensel lift.
 *                These polynomials are modified in place so that they may be used to resume lifting again.
 * @param Pptr the prime used in the modular equivalence.
 * @param e the exponent of the prime used in the modular equivalence of factors on input.
 * @param targetBound the minimum modulus to which factors are lifted. 
 * 
 * @return m, the smalled exponent such that Pptr->prime^m is greater than targetBound,
 *         or -1 if input polynomials do not satisfy input conditions. 
 *                
 */
int multiTermQuadraticPadicLiftResume(const DUZP_t* a, DUZP_t** liftedF, DUZP_t** sigmas, short nf, int e, const mpz_t targetBound, const Prime_ptr* Pptr);


/**
 * Given a primtive polynomial, a, a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)*liftedF[0]*...liftedF[nf-1] (hopefully) over Z.
 * This process is optimized and iterative in that it solves the multi-term diophantine equation directly
 * rather than a recursive scheme which uses only two-term diophantine equations.
 * The lifting here is done quadratically, that is, each step lifts the factors to the modulus squared.
 *
 * @param a a polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials over a finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param nf the number of factors in f to be lifted.
 * @param Pptr the prime of the finite field.
 * 
 * @return the exponent of the modulus where lifting terminated. That is, returns e where a = lc(a)*liftedF[0]*...liftedF[nf-1] mod Pptr->prime^e.
 */ 
static inline int multiTermQuadraticPadicLiftOpt_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {
	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);

	int ret = multiTermQuadraticPadicLiftStart(a, f, liftedF, NULL, nf, bound, Pptr);
	
	mpz_clear(bound);
	return ret;
}


/**
 * Given a polynomials a, u, and w, such that a = uw mod m and su + tw = 1 mod m, 
 * lift both congruences to mod m^2. That is, we compute 
 * u', w', s', t', such that a = u'w' mod m^2 and u's' + w't' = 1 mod m^2.
 * 
 * @note w must be monic.
 *
 * The lifted polynomials u', w', s', t', are returned in uu, ww, ss, tt, respectively.
 *
 * @param a a polynomial over the integers whose factors are to be lifted.
 * @param u a polynomial which is a factor of a modulo m.
 * @param w a monic polynomial which is a factor of a modulo m.
 * @param[in,out] m the modulo; on function exit it is set to m^2. 
 * @param[out] uu the lifted factor u' is returned in this pointer.
 * ...
 *
 * @return 1 iff the hensel step was effective i.e. the error was non-zero and u, w, s, and t, were updated.
 *
 */
// int quadraticPadicHenselStep_DUZP(const DUZP_t* a, const DUZP_t* u, const DUZP_t* w, const DUZP_t* s, const DUZP_t* t, const mpz_t m, DUZP_t** uu, DUZP_t** ww, DUZP_t** ss, DUZP_t** tt);





////////////////////////
// Lfits for GCDs
////////////////////////

/**
 * Compute the GCD of a and b using a modular method and reconstruct using Hensel lifting.
 * Input polynomials should both be monic.
 *
 * @param a a monic polynomial over the integers.
 * @param b a monic polynomial over the integers.
 * @param Pptr the prime to use for the modular method and subsequent lifting.
 *
 * @return the gcd of a and b over the integers.
 */
DUZP_t* monicUnivarHenselGCD_DUZP(const DUZP_t* a, const DUZP_t* b, const Prime_ptr* Pptr);


/**
 * Compute the GCD of a and b using a modular method and reconstruct using Hensel lifting.
 *
 * @param a a polynomial over the integers.
 * @param b a polynomial over the integers.
 * @param Pptr the prime to use for the modular method and subsequent lifting.
 *
 * @return the gcd of a and b over the integers.
 */
DUZP_t* univarHenselGCD_DUZP(const DUZP_t* a, const DUZP_t* b, const Prime_ptr* Pptr);


/**
 * Compute the GCD of a and b using a modular method and reconstruct using Hensel lifting.
 *
 * @param a a polynomial over the integers.
 * @param b a polynomial over the integers.
 * @param Pptr the prime to use for the modular method and subsequent lifting.
 *
 * @return the gcd of a and b over the integers.
 */
DUZP_t* univarQuadraticHenselGCD_DUZP(const DUZP_t* a, const DUZP_t* b, const Prime_ptr* Pptr);


#ifdef __cplusplus
}
#endif


#endif
