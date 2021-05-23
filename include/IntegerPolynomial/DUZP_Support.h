
#ifndef _DUZP_SUPPORT_
#define _DUZP_SUPPORT_

#ifdef __cplusplus
    extern "C" {
#endif

#include "../ModularPolynomial/DUSP_Support.h"
#include "../ModularPolynomial/DBSP_Support.h"
#include "../ModularPolynomial/DUSP_Support_Test.h"
#include "../FiniteFields/prime64_constants.h"
#include "../RationalNumberPolynomial/SMQP_Support-AA.h"
#include "SMZP_Support.h"
#include <gmp.h>

#define MAX(a, b) (a > b) ? a : b
#define MIN(a, b) (a < b) ? a : b


/**
 * Dense Univariate Integer Polynomial C struct.
 * Dense array of a coefficients, coefs.
 * coefs allocation size, alloc
 * last non-zero term, lt.
 */
typedef struct duzp {
	mpz_t* coefs;
	polysize_t alloc;
	polysize_t lt;
} DUZP_t;

/**
 * Allocate a DUZP with alloc number of coefficients.
 * i.e., alloc = maximum degree - 1, and the polynomial set to 0.
 * returns the newly allocated DUZP.
 */
static inline DUZP_t* makePolynomial_DUZP(int alloc) {
	if (alloc <= 0) {
		alloc = 1;
	}
	DUZP_t* p = (DUZP_t*) malloc(sizeof(DUZP_t));
	mpz_t* coef = (mpz_t*) malloc(sizeof(mpz_t)*alloc);
	mpz_init(coef[0]);
	p->coefs = coef;
	p->alloc = alloc;
	p->lt = 0;
	return p;
}

/**
 * Allocate a DUZP with alloc number of coefficients.
 * This polynomial has no initialized coeffiicents and it's degree is -1.
 * @return the newly allocated DUZP
 */
static inline DUZP_t* allocPolynomial_DUZP(int alloc) {
	if (alloc <= 0) {
		alloc = 1;
	}
	DUZP_t* p = (DUZP_t*) malloc(sizeof(DUZP_t));
	mpz_t* coef = (mpz_t*) malloc(sizeof(mpz_t)*alloc);
	p->coefs = coef;
	p->alloc = alloc;
	p->lt = -1;
	return p;
}

static inline DUZP_t* makeConstPolynomial_DUZP(int alloc, signed long z) {
	DUZP_t* ret = makePolynomial_DUZP(alloc);
	mpz_set_si(ret->coefs[0], z);
	return ret;
}

static inline DUZP_t* makeBigConstPolynomial_DUZP(int alloc, const mpz_t z) {
	DUZP_t* ret = makePolynomial_DUZP(alloc);
	mpz_set(ret->coefs[0], z);
	return ret;
}

/**
 * Free a DUZP.
 */
static inline void freePolynomial_DUZP(DUZP_t* p) {
	if (p != NULL) {
		polysize_t lt = p->lt;
		mpz_t* coefs = p->coefs;
		for (int i = 0; i <= lt; ++i) {
			mpz_clear(coefs[i]);
		}
		free(p->coefs);
		free(p);
	}
}

static inline void resizePolynomial_DUZP(DUZP_t* p, polysize_t allocSize) {
	if (p == NULL) {
		return;
	}

	if (allocSize == p->alloc) {
		return;
	}

	if (allocSize <= p->lt) {
		for (polysize_t i = p->lt; i >= allocSize; --i) {
			mpz_clear(p->coefs[i]);
		}
		p->lt = allocSize-1;
	}
	p->coefs = (mpz_t*) realloc(p->coefs, sizeof(mpz_t)*allocSize);
	p->alloc = allocSize;
}

/**
 * Evaluate a DUZP at a constant c
 * returns the resulting constant as val
 */
void evaluate_DUZP(const DUZP_t* p, const mpz_t c, mpz_t val);

/**
 * Deep copy a DUZP.
 * returns the new DUZP.
 */
DUZP_t* deepCopyPolynomial_DUZP(const DUZP_t* p);

/**
 * Create a copy of the polynomial a and store it in the (possibly)
 * pre-allocated DUZP_t* pointed to by bb.
 *
 * @param a the polynomial to be copied
 * @param bb a pointer to a (possibly) pre-allocated DUZP_t to store the copy.
 */
void deepCopyPolynomial_DUZP_inp(const DUZP_t* a, DUZP_t** bb);

/**
 * Convert a DUZP to a Dense Univariate Small Prime polynomial
 * using the prime contained in pptr.
 * returns the newly created DUSP.
 */
duspoly_t* convertToPrimeField_DUZP(const DUZP_t* p, const Prime_ptr* pptr);

/**
 * Convert a DUZP to a univariate polynomial over a finite field.
 * That is, obtain p mod pptr->prime.
 * The result is stored in the pointer pp, which may point to a pre-allocated
 * duspoly_t to re-use allocated space.
 *
 * @param p the DUZP to mod
 * @param pptr the prime pointer describing the prime field
 * @param[out] pp a pointer to the (possibly) pre-allocated duspoly_t to store the result.
 */
void convertToPrimeField_DUZP_inp(const DUZP_t* p, const Prime_ptr* pptr, duspoly_t** pp);

/**
 * Convert a DUSP directly to DUZP. Coffecicients are not
 * properly lifted but rather copied explicitly.
 * returns the newly created DUZP.
 */
DUZP_t* deepCopyPolynomial_DUZPFromspX(const duspoly_t* p, const Prime_ptr* pptr);

/**
 * Convert a DUSP to DUZP using symmetric range of finite field.
 * That is, if p->coefs[i] > prime / 2 then the returned coef
 * is p->coefs[i] - prime.
 * Coefficients are not lifted.
 *
 * returns the newly created DUZP.
 */
DUZP_t* deepCopyPolynomial_DUZPFromspXSymmetric(const duspoly_t* p, const Prime_ptr* pptr);

/**
 * Reduce the coefficients of the input polynomial by applying
 * the modulo mod. The polynomial is reduced in place.
 *
 * @param p the polynomial to be reduced.
 * @param mod the modulo to apply.
 */
void applyModulo_DUZP_inp(DUZP_t* p, const mpz_t mod);


/**
 * Reduce the coefficients of the input polynomial by applying
 * the modulo mod. The polynomial is reduced in place.
 * The mod is applied in a symmetric range so that coefficients
 * after being reduced are between -mod/2 and mod/2.
 *
 * @param p the polynomial to be reduced.
 * @param mod the modulo to apply.
 */
void applyModuloSymmetric_DUZP_inp(DUZP_t* p, const mpz_t mod);




/**
 * Determine if the input polynomial is 0.
 *
 * returns 1 iff a is 0.
 */
static inline int isZero_DUZP(const DUZP_t* a) {
	return (a == NULL) || (a->alloc == 0) || (a->lt == 0 && mpz_sgn(a->coefs[0]) == 0);
}

/**
 * Determine if two DUZP polynomials are equal.
 * returns non-zero if equal.
 */
int isEqual_DUZP(const DUZP_t* a, const DUZP_t* b);

/**
 * Print a DUZP using sym as the variable symbol.
 */
void printPoly_DUZP(const DUZP_t* poly, char* sym);


/**
 * Print a DUZP using sym as the variable symbol in maple format
 * where maplevar is the variable in maple to be assigned the value of poly.
 * i.e.: "maplevar:= poly;"
 */
void printPoly_DUZP_maple(DUZP_t* poly, char* sym, char* maplevar);

/**
 * Print a polynomial to the file pointer fp.
 *
 * @param fp the file pointer
 * @param poly the polynomial to print
 * @param sym the symbol for the polynomial variable.
 *
 */
void printPolyToFile_DUZP(FILE* fp, DUZP_t* poly, char* sym);

/**
 * Build a random DUZP given the maxDeg, bound of coefficent size in bits,
 * sparsity percentage, and boolean for in coefficients should be all positive
 * or randomly positive or negative.
 * returns the newly created DUZP.
 */
DUZP_t* buildRandomPoly_DUZP(polysize_t maxDeg, int coefBoundBits, float sparsity, int includeNeg);

/**
 * Get the infinity norm of the polynomial a. That is,
 * absolute value of the coefficient of a with maximum magnitude.
 */
static inline void infinityNorm_DUZP (const DUZP_t* a, mpz_t c) {
	if (a == NULL) {
		mpz_set_ui(c, 0ul);
		return;
	}

	mpz_abs(c, a->coefs[0]);
	for (int i = 1; i <= a->lt; ++i) {
		if (mpz_cmpabs(c, a->coefs[i]) < 0) {
			mpz_abs(c, a->coefs[i]);
		}
	}
}


/**
 * Get the content of a DUZP and return it in c.
 */
void content_DUZP(const DUZP_t* p, mpz_t c);

/**
 * Get the primitive part of a DUZP, p.
 * returns the primitive part of p.
 */
DUZP_t* primitivePart_DUZP(const DUZP_t* p);

/**
 * Convert the DUZP p to its primitive part in place.
 */
void primitivePart_DUZP_inp(DUZP_t* p);

/**
 * Convert the DUZP p to its primitive part in place.
 * Returns the content of p in cont.
 */
void primitivePartAndContent_DUZP_inp(DUZP_t* p, mpz_t cont);

/**
 * Get the content and primitive part of a DUZP p.
 * Content is returned in cont while function returns the primitive part.
 * returns the primitive part of p as a DUZP.
 */
DUZP_t* primitivePartAndContent_DUZP(const DUZP_t* p, mpz_t cont);

/**
 * Divive a DUZP, in place, by the integer z where division must be exact.
 */
static inline void divideByIntegerExact_DUZP_inp(DUZP_t* pp, const mpz_t z) {
	mpz_t* coefs = pp->coefs;
	polysize_t lt = pp->lt;
	for (polysize_t i = 0; i <= lt; ++i) {
		if (!mpz_divisible_p(coefs[i], z)) {
			fprintf(stderr, "not divisible yo\n");
			free((int*) -1);
			exit(1);
		}
		mpz_divexact(coefs[i], coefs[i], z);
	}
}

static inline void negatePolynomial_DUZP_inp(const DUZP_t* p) {
	mpz_t* coefs = p->coefs;
	polysize_t lt = p->lt;
	for (polysize_t i = 0; i <= lt; ++i) {
		mpz_neg(coefs[i], coefs[i]);
	}
}

/**
 * Multiply a DUZP, in place, by the gmp integer z.
 */
static inline void multiplyByInteger_DUZP_inp(DUZP_t* pp, const mpz_t z) {
	mpz_t* coefs = pp->coefs;
	polysize_t lt = pp->lt;
	for (polysize_t i = 0; i <= lt; ++i) {
		mpz_mul(coefs[i], coefs[i], z);
	}
}

/**
 * Multiply a DUZP, in place, by the unsigned long integer z.
 */
static inline void multiplyByIntegerUI_DUZP_inp(DUZP_t* pp, unsigned long z) {
	mpz_t* coefs = pp->coefs;
	polysize_t lt = pp->lt;
	for (polysize_t i = 0; i <= lt; ++i) {
		mpz_mul_ui(coefs[i], coefs[i], z);
	}
}

/**
 * Multiply a DUZP by the integer z.
 * returns the product of pp and z.
 */
static inline DUZP_t* multiplyByInteger_DUZP(const DUZP_t* pp, const mpz_t z) {
	polysize_t lt = pp->lt;
	DUZP_t* ret = makePolynomial_DUZP(lt+1);
	ret->lt = lt;
	mpz_t* ret_coefs = ret->coefs;
	mpz_t* coefs = pp->coefs;
	mpz_mul(ret_coefs[0], coefs[0], z);
	for (polysize_t i = 1; i <= lt; ++i) {
		mpz_init(ret_coefs[i]);
		mpz_mul(ret_coefs[i], coefs[i], z);
	}
	return ret;
}



/************************
 * Arithmetic Functions
 ************************/

/**
 * Subtract a duspoly_t from a DUZP, viewing the finite field poly as being
 * over the integers.
 *
 * returns the difference between a and b.
 */
DUZP_t* subtractDUSP_DUZP(DUZP_t* a, duspoly_t* b, const Prime_ptr* Pptr);

/**
 * Add two DUZP polynomials.
 *
 * returns the sum of a and b.
 */
DUZP_t* addPolynomials_DUZP(DUZP_t* a, DUZP_t* b);

/**
 * Add two DUZP polynomials, resuing the allocation
 * of a to store the sum in place of a.
 */
void addPolynomials_DUZP_inp(DUZP_t** a, DUZP_t* b);

/**
 * Add two DUZP polynomials together, storing the sum
 * in the pre-allocated polynomial pointed to by sump.
 * The polynomial pointed to by sump does not necessarily need to be
 * pre-allocated but should be.
 *
 * @param a one of the polynomials to add together.
 * @param b one of the polynomials to add together.
 * @param[out] sump a pointer to a pre-allocated polynomial to store the sum.
 */
void addPolynomialsPreAlloc_DUZP(DUZP_t* a, DUZP_t* b, DUZP_t** sump);

/**
 * Subtract two DUZP polynomials.
 *
 * returns the difference between a and b.
 */
DUZP_t* subtractPolynomials_DUZP(DUZP_t* a, DUZP_t* b);

/**
 * Subtract two DUZP polynomials, resuing the allocation
 * of a to store the difference in place of a.
 */
void subtractPolynomials_DUZP_inp(DUZP_t** a, DUZP_t* b);

/**
 * Subtract two polynomial, b from a, putting the difference back in b.
 */
void subtractPolynomials_DUZP_inpRHS(const DUZP_t* a, DUZP_t** b);

/**
 * Multiply a DUZP and a duspoly then apply mod p to return the
 * product as a duspoly_t.
 *
 * returns the product of a and b after reducing mod p.
 */
duspoly_t* multiplyAndModDUSP_DUZP(DUZP_t* a, duspoly_t* b, const Prime_ptr* Pptr);

/**
 * Multiply a DUZP and a duspoly then apply mod p to return the
 * product as a duspoly_t.
 * The product is returned in b.
 *
 */
void multiplyAndModDUSP_DUZP_inp(DUZP_t* a, duspoly_t** b, const Prime_ptr* Pptr);


/**
 * Multiply two DUZP polynomials together.
 * returns the product of a and b.
 */
DUZP_t* multiplyPolynomials_DUZP(const DUZP_t* a, const DUZP_t* b);


/**
 * Multiply two DUZP polynomials together, reusing
 * the allocation of bb to store the product.
 */
void multiplyPolynomials_DUZP_inp(const DUZP_t* a, DUZP_t** bb);


/**
 * Multiply all the polynomials in the list polys, of size n,
 * into a single product.
 *
 * returns the product of polynomials in the list polys.
 *
 */
DUZP_t* multiplyManyPolynomials_DUZP(DUZP_t const* const* polys, int n);


/**
 * Multiply all the polynomials in the list polys, of size n,
 * into a single product. The DUZP_t prod should point
 * to a pre-allocated DUZP_t to store the product.
 * The product is returned in prod.
 *
 */
void multiplyManyPolynomialsPreAlloc_DUZP(DUZP_t const* const* polys, int n, DUZP_t** prod);


/**
 * Multiply polynomials a and b together, using (possibly pre-allocated) polynomial c
 * to store the result.
 */
void multiplyPolynomialsPreAlloc_DUZP(const DUZP_t* a, const DUZP_t* b, DUZP_t** c);

/**
 * Multiply together all the polynomials in the list polys, of size n,
 * exlucding the polynomial at index idx.
 * The single DUZP_t product is returned.
 *
 * @param polys an array of DUZP_t polys of size n
 * @param n the size of the array polys
 * @param idx the index of the polynomial to exclude.
 */
DUZP_t* multiplyAllButOnePolynomials_DUZP(DUZP_t const*const* polys, int n, int idx);

/**
 * Multiply together all the polynomials in the list polys, of size n,
 * exlucding the polynomial at index idx. The product is returned
 * in prod, which should be, but not neccesarily be, pre-alloccated.
 *
 * @param polys an array of DUZP_t polys of size n.
 * @param n the size of the array polys.
 * @param idx the index of the polynomial to exclude.
 * @param[out] prod the product of the polys excluding index idx is returned in this pointer.
 */
void multiplyAllButOnePolynomialsPreAlloc_DUZP(DUZP_t const*const* polys, int n, int idx, DUZP_t** prod);

/**
 * Multiply the polynomial pointed to by a_ptr by the monic binomial (x - b).
 * This operation is done in-place, eseentially by in-place arithmetic and a shift.
 *
 * @param a_ptr a pointer to the polynomial to be multiplied.
 * @param b the constant term of the monic binomial multiplier.
 *
 */
void multiplyByBinomial_DUZP_inp(DUZP_t** a_ptr, const mpz_t b);

/**
 * Divide the polynomial (a) by a monic linear polynomial (x - b).
 *
 * @param a the dividend
 * @param b the constant term of the divisor, (x-b);
 * @param[out] q a pointer to the resulting quotient; may be NULL;
 * @param[out] rem a pointer to the resulting remainder; may be NULL.
 *
 * @return 1 if the division was exact,
 *         0 if the division was not exact,
 *        -1 if an error occurred (such as invalid parameters;
 */
int divideByMonicLinear_DUZP(const DUZP_t* a, const mpz_t b, DUZP_t** q, mpz_t* rem);

/**
 * Divide polynomial a by polynomial b. Quotient is returned in q
 * and remainder returned in r. This functions returns 1 or 0
 * to indicate if the division works over the integers.
 * Specifically, if the leading coefficient of b is a divisor
 * of contents of both a and b. If this function returns 0
 * then q and r both point to NULL.
 *
 * returns 1 iff the division is possible over the integers.
 */
int dividePolynomials_DUZP(const DUZP_t* a, const DUZP_t* b, DUZP_t** q, DUZP_t** r);

/**
 * Determine if DUZP a is divisible by DUZP b.
 * If the division is exact, the quotient is returned in q if non-NULL.
 * returns non-zero if a is divisible by b.
 */
int divideTest_DUZP(const DUZP_t* a, const DUZP_t* b, DUZP_t** q);

/**
 * Divide the polynomial pointed to by aa by b.
 * The remainder is returned in aa in-place.
 * If the division could not be carried out over the integers
 * then 0 is returned and the polynomial pointed to by aa is
 * left in a valid yet undetermined state.
 *
 * @param[in,out] aa a pointer to the dividend and location of remainder
 * @param b the divisor polynomial.
 *
 * @return 1 iff the division could be carried out over the integers
 */
int remainder_DUZP_inp(DUZP_t** aa, DUZP_t* b);

int remainderMod_DUZP_inp(DUZP_t** aa, const DUZP_t* b, const mpz_t mod);

/**
 * Get the GCD of DUZP polynomials a and b.
 * Performs the GCD by many small primes and CRT.
 * returns the gcd of a and b as a DUZP.
 */
DUZP_t* GCD_DUZP(DUZP_t* a, DUZP_t* b);

/**
 * Get the GCD of two DUZP primitive polynomials.
 * returns the gcd of a and b as a DUZP.
 */
DUZP_t* primitiveGCD_DUZP(const DUZP_t* a, const DUZP_t* b);

/**
 * Get the derivative of the polynomial pointed to by a.
 * The derivative is returned in a in place.
 */
void differentiate_DUZP_inp(DUZP_t** a);

/**
 * Get the derivative of the polynomial pointed to by a.
 * returns the derivative of a as a DUZP.
 */
inline DUZP_t* differentiate_DUZP(const DUZP_t* a) {
	DUZP_t* b = deepCopyPolynomial_DUZP(a);
	differentiate_DUZP_inp(&b);
	return b;
}


/**
 * Modular Resultant over Z
 * @note non-determinisitic algorithm (determinisitic=0) works faster in practice.
 *
 * @param a polynomial
 * @param b polynomial
 * @param uBound is returned upper-bound used to choose 64-bits primes
 * @param deterministic has two states: 0, 1
 * @param hgcd if 1, it uses hgcd-based subresultant chain algoirthm [see MCA, chapter 13]
 * @return resultant (a, b)
 *
 */
DUZP_t* modularResultant_DUZP (const DUZP_t* a, const DUZP_t* b, mpz_t uBound, int deterministic, int hgcd);

/**
 * Modular Resultant over Z
 * @note this is a deterministic algorithm.
 *
 * @param a polynomial
 * @param b polynomial
 * @param chain_size is returned the size of subresultantChain
 * @param uBound is returned upper-bound used to choose 64-bits primes
 * @param k returns k-th subresultant chain
 * @param hgcd if 1, it uses hgcd-based subresultant chain algoirthm [see MCA, chapter 13]
 * @return subresultantChain (a, b)
 */
DUZP_t** modularSubresultantChain_DUZP (const DUZP_t* a, const DUZP_t* b, int* chain_size, mpz_t uBound, int deterministic);
DUZP_t** modularSubresultantAtDeg_DUZP (const DUZP_t* a, const DUZP_t* b, int k, int* chain_size, mpz_t uBound, int deterministic, int hgcd);


/**
 * Speculative Subresultant algorithm over Z[x]
 * 
 * @param a polynomial
 * @param b polynomial
 * @param k in range(0, deg(b))
 * @param chain_size is returned the size of the chain 
 * @param uBound is returned upper-bound used to choose 64-bits primes 
 * @param uspecInfo a list of specAQR_spX_t* containging speculative information [A, Q, R]
 * @param inf_size is returned the size of uspecInfo 
 * @param results_mode if it's 0 : computes the entire subresultant, otherwise : computes only its principle coeff.
 */
DUZP_t** regularGCDUnivariateSpeculativeSRC_DUZP (const DUZP_t* a, const DUZP_t* b, int k, int* chain_size, 
											polysize_t *degs, mpz_t uBound, specAQRArray_spX_t **uspecInfoArray, 
											int *info_size, int results_mode);

/************************
 * Conversion Functions
 ************************/

/**
 * Convert an AltArr_t to a DUZP, failing
 * if not univariate or not all integer coefficeints.
 * returns the newly created DUZP.
 */
DUZP_t* convertFromAltArr_DUZP(const AltArr_t* aa);

/**
 * Convert an AltArrZ_t to a DUZP, failing
 * if not univariate.
 * returns the newly created DUZP.
 */
DUZP_t* convertFromAltArrZ_DUZP(const AltArrZ_t* aa);

/**
 * Convert a duspoly_t to a DUZP.
 *
 * returns the newly created DUZP.
 */
DUZP_t* convertFromDUSP_DUZP(const duspoly_t* a, const Prime_ptr* p);

DUZP_t* convertFromDUSP_DUZP_inRange (const duspoly_t* a, const Prime_ptr* Pptr);


/**
 * Convert a duspoly_t to a DUZP using the symmetric range of the
 * finite field.
 *
 * returns the newly created DUZP.
 */
static inline DUZP_t* convertFromDUSPSymmetric_DUZP(const duspoly_t* a, const Prime_ptr* p) {
	return deepCopyPolynomial_DUZPFromspXSymmetric(a, p);
}

/**
 * Apply prime and convert to DUSP.
 *
 * returns the newly created duspoly_t.
 */
duspoly_t* convertToDUSP_DUZP(const DUZP_t* a, const Prime_ptr* p);

/**
 * Convert a DUZP to an AltArr_t.
 * returns the newly created AltArr_t.
 */
AltArr_t* convertToAltArr_DUZP(const DUZP_t* poly);

/**
 * Convert a DUZP to an AltArrZ_t.
 * returns the newly created AltArrZ_t.
 */
AltArrZ_t* convertToAltArrZ_DUZP(const DUZP_t* poly);

#ifdef __cplusplus
}
#endif


#endif
