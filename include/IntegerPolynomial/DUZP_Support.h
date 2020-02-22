
#ifndef _DUZP_SUPPORT_
#define _DUZP_SUPPORT_

#ifdef __cplusplus
    extern "C" {
#endif

#include "../ModularPolynomial/DUSP_Support.h"
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
typedef struct {
	mpz_t* coefs;
	polysize_t alloc;
	polysize_t lt;
} DUZP_t;


/**
 * Allocate a DUZP with alloc number of coefficients. 
 * i.e., alloc = maximum degree - 1.
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
 * Convert a DUZP to a Dense Univariate Small Prime polynomial
 * using the prime contained in pptr.
 * returns the newly created DUSP.
 */
duspoly_t* convertToPrimeField_DUZP(const DUZP_t* p, const Prime_ptr* pptr);

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
void content_DUZP(DUZP_t* p, mpz_t c);

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



/***********************
 * DBZP Polynomials 
 **********************/

/** 
 * Dense Bivariate Integer Polynomial C struct Z[x > y].
 * Dense array of DUZP_t polynomials over Z[y], coefsY.
 * partial leading term based on X, ltX and coefsY allocation size, allocX. 
 */
typedef struct {
	polysize_t ltX;
	polysize_t allocX;
	DUZP_t** coefsY;
} DBZP_t;


/**
 * Allocate a DBZP with alloc number of coefficients. 
 * i.e., alloc = partial maximum degree - 1.
 * returns the newly allocated DBZP.
 */
static inline DBZP_t* makePolynomial_DBZP(int allocX) {
	if (allocX <= 0) {
		allocX = 1; 
	}
	DBZP_t* p = (DBZP_t*) malloc(sizeof(DBZP_t));
	p->coefsY = (DUZP_t**) malloc (allocX*sizeof(DUZP_t*));
	for (polysize_t i = 0; i < allocX; i++) {
		p->coefsY[i] = NULL;
	}
	p->allocX = allocX;
	p->ltX = 0;
	return p;
}

/**
 * Make a constant DBZP_t polynomial
 * 
 * @param allocX allocation size
 * @param z coef (signed long)
 * @return 
 */
static inline DBZP_t* makeConstPolynomial_DBZP(int allocX, signed long z) {
	DBZP_t* ret = makePolynomial_DBZP(allocX);
	ret->coefsY[0] = makeConstPolynomial_DUZP (1, z);
	return ret;
}

/**
 * Make a constant DBZP_t polynomial
 * 
 * @param allocX allocation size
 * @param z coef (mpz_t)
 * @return 
 */
static inline DBZP_t* makeBigConstPolynomial_DBZP(int allocX, const mpz_t z) {
	DBZP_t* ret = makePolynomial_DBZP(allocX);
	ret->coefsY[0] = makeBigConstPolynomial_DUZP (1, z);
	return ret;
}

/**
 * Free a DBZP.
 * @note
 */
static inline void freePolynomial_DBZP(DBZP_t* p) {
	if (p != NULL) {
		polysize_t lt = p->ltX;
		for (polysize_t i = 0; i <= lt; i++) {
			freePolynomial_DUZP (p->coefsY[i]);
		}
		free(p->coefsY);
		free(p);
	}
}

/**
 * Resize a DBZP polynomial 
 * @note this works in-place
 * 
 * @param p DBZP_t polynomial 
 * @param allocSize new allocation size
 */
static inline void resizePolynomial_DBZP(DBZP_t* p, polysize_t allocSize) {
	if (p == NULL) {
		return;
	}

	if (allocSize == p->allocX) {
		return;
	}

	if (allocSize <= p->ltX) {
		for (polysize_t i = p->ltX; i >= allocSize; --i) {
			freePolynomial_DUZP (p->coefsY[i]);
		}
		p->ltX = allocSize-1;
	}
	p->coefsY = (DUZP_t**) realloc (p->coefsY, sizeof(DUZP_t*)*allocSize);
	for (polysize_t i = p->ltX+1; i < allocSize; i++) {
		p->coefsY[i] = NULL;
	}
	p->allocX = allocSize;
}

/**
 * Evaluate coefsY of a DBZP at a constant c
 * 
 * @param p DBZP_t polynomial 
 * @param c a evaluating point
 * @return a bivariate polynomial (DBZP_t) over Z[x].
 */
static inline DBZP_t* partialEvaluate_DBZP(const DBZP_t* p, const mpz_t c) {
	if (p == NULL) {
		return NULL;
	}

	polysize_t ltX = p->ltX;
	DBZP_t* px = makePolynomial_DBZP (ltX+1);
	px->ltX = ltX;
	mpz_t vals;
	for (polysize_t i = 0; i <= ltX; i++) {
		mpz_init (vals);
		evaluate_DUZP (p->coefsY[i], c, vals);
		px->coefsY[i] = makeBigConstPolynomial_DUZP (1, vals);
		mpz_clear (vals);
	}
	
	polysize_t newlt = ltX;
	while (newlt >= 0 && px->coefsY[newlt] == NULL) {
		--newlt;
	}
	if (newlt < 0) {
		freePolynomial_DBZP (px);
		return NULL;
	}
	resizePolynomial_DBZP (px, newlt+1);
	return px;
}

/**
 * Evaluate a DBZP \in Z[x,y] at (x,y) = (c,d)
 * 
 * @param p DBZP_t polynomial 
 * @param c a X evaluating point
 * @param d a Y evaluating point
 * @param var result of evaluation
 */
static inline void evaluate_DBZP(const DBZP_t* p, const mpz_t c, const mpz_t d, mpz_t var) {
	if (p == NULL) {
		return;
	}

	polysize_t ltX = p->ltX;
	DBZP_t* px = makePolynomial_DBZP (ltX+1);
	px->ltX = ltX;
	mpz_t vals;
	// evaluate at y=d
	for (polysize_t i = 0; i <= ltX; i++) {
		mpz_init (vals);
		evaluate_DUZP (p->coefsY[i], d, vals);
		px->coefsY[i] = makeBigConstPolynomial_DUZP (1, vals);
		mpz_clear (vals);
	}
	
	polysize_t newlt = ltX;
	while (newlt >= 0 && px->coefsY[newlt] == NULL) {
		--newlt;
	}
	if (newlt < 0) {
		freePolynomial_DBZP (px);
		return;
	}

	// evaluate at x=c
	mpz_t result;
	mpz_init_set (result, px->coefsY[newlt]->coefs[0]);
	for (polysize_t i = newlt-1; i >= 0; i++) {
		mpz_mul (result, result, c); 
		if (px->coefsY[i] != NULL) {
			mpz_add (result, result, px->coefsY[i]->coefs[0]);
		}
	}

	mpz_init_set (var, result);
}

/**
 * Deep copy a DBZP.
 * @return the new DBZP.
 */
static inline DBZP_t* deepCopyPolynomial_DBZP(const DBZP_t* p) {
	if (p == NULL) {
		return NULL;
	}

	polysize_t ltX = p->ltX;
	DBZP_t* dp = makePolynomial_DBZP (ltX+1);
	for (polysize_t i = 0; i < ltX; i++) {
		dp->coefsY[i] = deepCopyPolynomial_DUZP (p->coefsY[i]);
	}
	dp->ltX = ltX;
	return dp;
}


/**
 * (distributed) Bivariate Subresultant Chain Data-Type over Z 
 * @note each subresultant-chain is stored in one C-array
 * @note using distributed-data-type led to a better performance for modularSubres over Z[x>y].
 */
typedef struct {
    polysize_t n; // Number of polynomials in a subresultant-chain, at least 2
    polysize_t** deg; // Partial degree of each polynomial (bi-var: deg[2])
    polysize_t* size; // size of each polynomial array
    mpz_t** coefs; // coefficients of polynomials 
} biSubresZ_t;

static inline biSubresZ_t* makeBiSubresultant_DBZP (polysize_t n) {
    if (n < 0) {
        return NULL;
    }
    biSubresZ_t* res = (biSubresZ_t*) malloc (sizeof(biSubresZ_t));
    res->n = n;
    if (n > 0) {
        if (n > 1) {
            res->deg = (polysize_t**) malloc (n*sizeof(polysize_t*));
            res->coefs = (mpz_t**) malloc (n*sizeof(mpz_t*)); 
        }
        res->size = (polysize_t*) calloc (n, sizeof (polysize_t));
        for (polysize_t i = 0; i < n; i++) {
            res->deg[i] = (polysize_t*) calloc (2, sizeof(polysize_t));  
        }
    } else {
        res->n = 0;
        res->deg = NULL;
        res->size = NULL;
        res->coefs = NULL;
    }
	return res;
}

static inline void freeBiSubresultant_DBZP (biSubresZ_t* s) {
    if (s == NULL) {
        return;
    }
    if (s->n) {
        for (polysize_t i = 0; i < s->n; i++) {
            free (s->deg[i]);
            if (s->size[i]) {
				for (polysize_t j = 0; j < s->size[i]; j++) {
					mpz_clear(s->coefs[i][j]);
				}
				free (s->coefs[i]);
			}	
        }
        free (s->deg);
        free (s->size);
		free (s->coefs);
    }
    free (s);
}

static inline biSubresZ_t* deepCopyBiSubresultant_DBZP (biSubresZ_t* a) {
	if (a == NULL) {
		return NULL;
	}

	biSubresZ_t* ac = makeBiSubresultant_DBZP (a->n);
	if (!a->n) {
		return ac;
	}
	for (polysize_t i = 0; i < a->n; i++) {
		ac->deg[i][0] = a->deg[i][0];
		ac->deg[i][1] = a->deg[i][1];
		ac->size[i]   = a->size[i];
		if (a->size[i]) {
			for (polysize_t j = 0; j < a->size[i]; j++) {
				mpz_set (ac->coefs[i][j], a->coefs[i][j]);
			}
		}
	}
	return ac;
}

static inline int isBiSubresultantEquals_DBZP (biSubresZ_t* a, biSubresZ_t* b) {
	if (a == NULL && b == NULL) {
		return 1;
	} else if (a == NULL) {
		return 0;
	} else if (b == NULL) {
		return 0;
	}
	if (a->n != b->n) {
		// fprintf (stderr, "[DBZP] a->n != b->n\n");
		return 0;
	}
	for (polysize_t i = 0; i < a->n; i++) {
		if (a->deg[i][0] != b->deg[i][0] || a->deg[i][1] != b->deg[i][1]) {
			// fprintf (stderr, "[DBZP] a->deg[i][0] != b->deg[i][0] || a->deg[i][1] != b->deg[i][1]\n");
			return 0;
		} else if (a->size[i] != b->size[i]) {
			// fprintf (stderr, "[DBZP] a->size[i] != b->size[i]\n");
			return 0;
		} 
		for (polysize_t j = 0; j < a->size[i]; j++) {
			if (mpz_cmp (a->coefs[i][j], b->coefs[i][j]) ) {
				// fprintf (stderr, "[DBZP] a->coefs[i][j] != b->coefs[i][j]\n");
				return 0;
			}
		}
	}
	return 1;
}

static inline AltArrZ_t* convertToAltArrZ_DBZP (const mpz_t* a, polysize_t* apd) 
{
	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	AltArrZ_t* aZZ = NULL;
	
	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t));
	int n = apd[0] + 1;
	for (polysize_t i = 0; i <= apd[1]; i++) {
		for (polysize_t j = 0; j <= apd[0]; j++) {
			if (mpz_sgn(a[n*i+j])) {
				// fprintf (stderr, "deg[0] = %ld , deg[1] = %ld\n", j, i);
				// gmp_fprintf (stderr, " a[n*i+j] = %Zd\n", a[n*i+j]);
				// a true monomial with [j][i]
				degs[0] = i; degs[1] = j;
				if (aZZ == NULL) {
					aZZ = makeConstPolynomial_AAZ(1, 2, a[n*i+j]);
					setDegrees_AAZ_inp(aZZ, 0, degs, 2);
				} else {
					setCoefficient_AAZ (aZZ, degs, 2, a[n*i+j]);
				}
			}
		}
	}
	free (degs);
	if (aZZ != NULL) {
		mergeSortPolynomial_AAZ (aZZ);
	}
	return aZZ;
}


static inline AltArrZ_t* convertToAltArrZ_DBSP (const elem_t* a, polysize_t* apd, Prime_ptr* Pptr) 
{
	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	AltArrZ_t* aZZ = NULL;
	
	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t));
	mpz_t tmp;
	mpz_init (tmp);
	int n = apd[0] + 1;
	for (polysize_t i = 0; i <= apd[1]; i++) {
		for (polysize_t j = 0; j <= apd[0]; j++) {
			if (a[n*i+j]) {
				degs[0] = i; degs[1] = j;
				// gmp_fprintf (stderr, " a[n*i+j] = %Zd\n", a[n*i+j]);
				mpz_set_si (tmp, smallprimefield_convert_out (a[n*i+j], Pptr));
				if (aZZ == NULL) {
					aZZ = makeConstPolynomial_AAZ(1, 2, tmp);
					setDegrees_AAZ_inp(aZZ, 0, degs, 2);
				} else {
					setCoefficient_AAZ (aZZ, degs, 2, tmp);
				}
			}
		}
	}
	mpz_clear (tmp);
	free (degs);
	if (aZZ != NULL) {
		mergeSortPolynomial_AAZ (aZZ);
	}	
	return aZZ;
}

// don't need to init mpz_t* a 
static inline mpz_t* convertFromAltArrZ_DBZP (AltArrZ_t* aZ, polysize_t* apd) 
{
	if (isZero_AAZ (aZ) || aZ->nvar != 2) {
		return NULL;
	}

	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	if (apd == NULL) {
		fprintf (stderr, "DUZP Error: apd must be initialized\n");
		exit (1);
	}
	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t)); 
	partialDegrees_AAZ (aZ, degs);
	apd[0] = degs[1]; 
	apd[1] = degs[0]; 

	int allocSz = (degs[1]+1) * (degs[0]+1);
	// fprintf(stderr, "allocSz = %d\n", allocSz);
	mpz_t* a = (mpz_t*) malloc (allocSz *sizeof(mpz_t));
	for (int k = 0; k < allocSz; k++) {
		mpz_init_set_d (a[k], 0);
	}
	int n = apd[1]+1;
	for (int k = 0; k < aZ->size; k++) {
		partialDegreesTerm_AAZ (aZ, k, degs);
		// fprintf (stderr, "a[%d] set!\n", n*degs[1]+degs[0]);
		mpz_set (a[n*degs[1]+degs[0]], aZ->elems[k].coef);
		degs[0] = degs[1] = 0;
	}
	// for (int i = 0; i < allocSz; i++) {
	// 	gmp_fprintf (stderr, "a[%d] = %Zd\n", i, a[i]);
	// }
	free (degs);
	return a;
}

static inline elem_t* convertFromAltArrZ_DBSP (AltArrZ_t* aZ, polysize_t* apd, Prime_ptr* Pptr) 
{
	if (isZero_AAZ (aZ) || aZ->nvar != 2) {
		return NULL;
	}

	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	if (apd == NULL) {
		fprintf (stderr, "DUZP Error: apd must be initialized\n");
		exit (1);
	}
	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t)); 
	partialDegrees_AAZ (aZ, degs);
	apd[0] = degs[1]; 
	apd[1] = degs[0]; 

	int allocSz = (degs[1]+1) * (degs[0]+1);
	// fprintf(stderr, "allocSz = %d\n", allocSz);
	elem_t* a = (elem_t*) calloc (allocSz, sizeof(elem_t));
	mpz_t tmp;
	mpz_init (tmp);
	int n = apd[1]+1;
	for (int k = 0; k < aZ->size; k++) {
		partialDegreesTerm_AAZ (aZ, k, degs);
		// fprintf (stderr, "a[%d] set!\n", n*degs[1]+degs[0]);
		mpz_mod_ui (tmp, aZ->elems[k].coef, (unsigned long) Pptr->prime);
		a[n*degs[1]+degs[0]] = mpz_get_si (aZ->elems[k].coef);
		a[n*degs[1]+degs[0]] = smallprimefield_convert_in (a[n*degs[1]+degs[0]], Pptr);
		degs[0] = degs[1] = 0;
	}
	// for (int i = 0; i < allocSz; i++) {
	// 	fprintf (stderr, "a[%d] = %Zd\n", i, a[i]);
	// }
	free (degs);
	return a;
}


static inline void printDBZP_AAZ (mpz_t* a, polysize_t* apd)
{
  if (apd[0]==0 || apd[1]==0) {
    fprintf (stderr, "0;\n");
  }

  polysize_t ad = apd[1] + 1;
  for (int i = 0; i <= apd[0]; i++) {
    for (int j = 0; j <= apd[1]; j++) {
      gmp_fprintf (stderr, "%Zd\t",a[ad*i+j]);
    }
    fprintf (stderr, "\n");
  }
}

static inline void printDBSP_AAZ (elem_t* a, polysize_t* apd, Prime_ptr* Pptr)
{
  if (apd[0]==0 || apd[1]==0) {
    fprintf (stderr, "0;\n");
  }

  polysize_t ad = apd[1] + 1;
  for (int i = 0; i <= apd[0]; i++) {
    for (int j = 0; j <= apd[1]; j++) {
      fprintf (stderr, "%lld\t", smallprimefield_convert_out (a[ad*i+j], Pptr));
    }
    fprintf (stderr, "\n");
  }
}

/**
 * Modular Bivariate Subresultant Chain working with a prime-set 
 * @note see fourier_primes_u64_table
 * 
 * @param g a distributed 2D polynomial over Z
 * @param gpd partial degrees gpd[0], gpd[1]
 * @param f a distributed 2D polynomial over Z s.t. deg (g, x) > deg (f, x)
 * @param fpd partial degrees fpd[0], fpd[1]
 * @param m returned-middle-computation in CRT
 * @param startIdx starting from prime64_ptr[startIdx] in CRT  
 * @param primeSz size of 'primes' list: {startIdx+0, startIdx+1, ..., startIdx+primeSz-1}
 * @return biSubresZ_t mod \prod_{i} primes[i] with CRT properties
 */
biSubresZ_t* modularSetBiSubresultant_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd, mpz_t m, int startIdx, int primeSz);

/**
 * Modular Bivariate Subresultant Chain over Z
 * @note see fourier_primes_u64_table
 * 
 * @param g a distributed 2D polynomial over Z
 * @param gpd partial degrees gpd[0], gpd[1]
 * @param f a distributed 2D polynomial over Z s.t. deg (g, x) > deg (f, x)
 * @param fpd partial degrees fpd[0], fpd[1]
 * @return biSubresZ_t mod \prod_{i} primes[i] with CRT properties
 */
biSubresZ_t* modularBiSubresultant_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd);


#ifdef __cplusplus
}
#endif


#endif
