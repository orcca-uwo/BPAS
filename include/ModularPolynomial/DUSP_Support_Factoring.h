
#ifndef DUSP_SUPPORT_FACTORING_H
#define DUSP_SUPPORT_FACTORING_H


#ifdef __cplusplus
extern "C" {
#endif

#include "./DUSP_Support.h"
#include "./DUSP_Support_Test.h"
#include <time.h>
#include <math.h>

#if TIMER_SUBPROGRAMS
#include "Unix_Timer.h"
#endif 

/**
 * Initial factors_t data-type
 * Given the alloc size and return an array of zero polys. 
 */
static inline factors_t* initFactors_spX (polysize_t alloc) 
{
	if (alloc == 0) {
		return NULL;
	}

	factors_t* f = (factors_t*) malloc (sizeof(factors_t));
	f->polys = (duspoly_t**) calloc (alloc, sizeof(duspoly_t*));
	f->exps  = (polysize_t*) calloc (alloc, sizeof(polysize_t));
	f->alloc = alloc;

	return f;
}

/**
 * Free factors_t data-type 
 */
static inline void freeFactors_spX (factors_t** f)
{
	if (*f == NULL){
		return;
	}

	for (polysize_t i = 0; i < (*f)->alloc; i++){
		freePolynomial_spX (&((*f)->polys[i]));
	}
	free ((*f)->polys);
	free ((*f)->exps);
	(*f)->alloc = 0;

	free (*f);
	*f = NULL;
}


/**
 * Initial factsll_t data-type
 * Given the exponent value, return a node with zero poly. 
 */
static inline factsll_t* initFactsll_spX (polysize_t exp)
{
	factsll_t* node = (factsll_t*) malloc (sizeof (factsll_t));
	node->exp = exp;
	node->poly = NULL;
	node->next = NULL;

	return node;
}


/**
 * Free factsll_t data-type
 */
static inline void freeFactsll_spX (factsll_t** f)
{
	if (*f == NULL){
		return;
	}

	factsll_t* tmp = NULL;

	while (*f != NULL) {
		tmp = *f;
		*f = (*f)->next;

		freePolynomial_spX (&(tmp->poly));
		free(tmp);
	}

	free (*f);
	*f = NULL;
}

/**
 * Given a factors_t of factors, 
 * return 1, if it is a empty array of factors and 0, otherwise
 */
static inline int isZeroFactors_spX (factors_t* f) 
{
	if (f == NULL){
		return 1;
	}

	if (f->alloc == 0) {
		return 1;
	}

	return 0;
}

/**
 * Given a factsll_t of factors, 
 * return 1, if it is a empty linked-list of factors and 0, otherwise
 */
static inline int isZeroFactsll_spX (factsll_t* f) 
{
	if (f == NULL){
		return 1;
	}

	return 0;
}

/** 
 * Deep copy of a factors_t data-type of factors 
 */
static inline factors_t* deepCopyFactors_spX (factors_t* f)
{
	if (isZeroFactors_spX(f)) {
		return NULL;
	}

	factors_t* f2 = (factors_t*) malloc (sizeof(factors_t));
	f2->polys = (duspoly_t**) malloc (sizeof(duspoly_t*)*(f->alloc));
	f2->exps = (polysize_t*) malloc (sizeof(polysize_t)*(f->alloc));
	f2->alloc = f->alloc;

	for (polysize_t i = 0; i < f->alloc; i++) {
		f2->polys[i] = deepCopyPolynomial_spX (f->polys[i]);
		f2->exps[i] = f->exps[i];
	}

	return f2;
}

/** 
 * Deep copy of a factsll_t data-type of factors 
 */
static inline factsll_t* deepCopyFactsll_spX (factsll_t* f)
{
	if (isZeroFactsll_spX(f)) {
		return NULL;
	}

	factsll_t* node = (factsll_t*) malloc (sizeof (factsll_t));

	factsll_t* head = NULL;
	factsll_t* tail = NULL; 


	factsll_t* cur = f;
	while (cur != NULL) {
		
		node->exp  = cur->exp;
		node->poly = deepCopyPolynomial_spX (cur->poly);
		node->next = NULL;

		if (head == NULL) {
			head = node;
			tail = node;
		} else {
			tail->next = node;
			tail = tail->next;
		}

		cur = cur->next;

		if (cur != NULL) {
			node = (factsll_t*) malloc (sizeof (factsll_t));
		}
	}

	return head;
}

/** 
 * Determine if two factors_t of factors are equal or not.
 * Note this function checks not only polynoimals but also exponents in order.
 * Note it is a lazy implementation of isEqualFacotrs_spX !
 */
static inline int isEqualFactors_spX (factors_t* f1, factors_t* f2)
{
	if (isZeroFactors_spX (f1)) {
		if (isZeroFactors_spX (f2)) {
			return 1;
		} else {
			return 0;
		}
	} else if (isZeroFactors_spX (f2)) {
		return 0;
	}

	if (f1->alloc != f2->alloc) {
		return 0;
	}

	polysize_t i;
	for (i = 0; i < f1->alloc; i++) {
		if (!isEqual_spX (f1->polys[i], f2->polys[i]) || (f1->exps[i] != f2->exps[i]) ) {
			return 0;
		}
	}

	return 1;
}


/** 
 * Determine if two factsll_t of factors are equal or not.
 * Note this function checks not only polynoimals but also exponents in order.
 * Note it is a lazy implementation of isEqualFactsll_spX !
 */
static inline int isEqualFactsll_spX (factsll_t* f1, factsll_t* f2)
{
	if (isZeroFactsll_spX (f1)) {
		if (isZeroFactsll_spX (f2)) {
			return 1;
		} else {
			return 0;
		}
	} else if (isZeroFactsll_spX (f2)) {
		return 0;
	}

	factsll_t* cur1 = f1;
	factsll_t* cur2 = f2;

	while (cur1 != NULL && cur2 != NULL) {
		if (!isEqual_spX (cur1->poly, cur2->poly) || cur1->exp != cur2->exp) {
			return 0;
		}
		
		cur1 = cur1->next;
		cur2 = cur2->next;
	}

	if (cur1 != NULL || cur2 != NULL) {
		return 0;
	}

	return 1;
}

/**
 * Return size of a factors_t of factors 
 * Note this function doesn't consider the allocation size, and
 * return the size of non-zero factors.
 */
static inline polysize_t lenFactors_spX (factors_t* f) 
{
	if (isZeroFactors_spX (f)) {
		return 0;
	}

	polysize_t len = f->alloc;
	while (f->polys[len-1] == NULL) {
		len--;
	}

	return len;
}

/**
 * Return size of a factsll_t of factors 
 */
static inline polysize_t lenFactsll_spX (factsll_t* f) 
{
	if (isZeroFactsll_spX (f)) {
		return 0;
	}

	factsll_t* cur = f;
	polysize_t len = 0;
	while (cur != NULL) {
		len++;
		cur = cur->next;
	}

	return len;
}

/** 
 * Concatenate two Factsll_t Blindly! 
 */
static inline factsll_t* concatFactsll_spX (factsll_t* f1, factsll_t* f2)
{
	if (isZeroFactsll_spX (f1)) {
		return deepCopyFactsll_spX (f2);
	} else if (isZeroFactsll_spX (f2)) {
		return deepCopyFactsll_spX (f1);
	}

	factsll_t* cur = f1;

	factsll_t* node = NULL;
	factsll_t* f1f2 = NULL;
	factsll_t* tail = NULL;

	while (cur != NULL) {

		node = initFactsll_spX (cur->exp);
		node->poly = deepCopyPolynomial_spX (cur->poly);

		if (f1f2 == NULL) {
			f1f2 = node;
			tail = node;
		} else {
			tail->next = node;
			tail = node;
		}

		cur = cur->next;
	}

	cur = NULL; // TEST

	cur = f2;
	while (cur != NULL) {

		node = initFactsll_spX (cur->exp);
		node->poly = deepCopyPolynomial_spX (cur->poly);

		if (f1f2 == NULL) {
			f1f2 = node;
			tail = node;
		} else {
			tail->next = node;
			tail = node;
		}

		cur = cur->next;
	}

	return f1f2;
}

/**
 * Getting a factor with giving exp value 
 */
static inline duspoly_t* getFactorsAtExp_spX (factors_t* f, polysize_t exp)
{
	fprintf(stderr, "DUSP: sorry! getFactorsAtExp_spX is not implemented.\n");
	exit(1);
}


/**
 * Getting a factor with giving exp value 
 * Note this implementation is lazy! 
 */
static inline duspoly_t* getFactorsAtExpll_spX (factsll_t* f, polysize_t exp)
{
	if (isZeroFactsll_spX (f)) {
		return NULL;
	}

	factsll_t* cur = f;
	while (cur != NULL) {
		if (cur->exp == exp) {
			return deepCopyPolynomial_spX (cur->poly);
		}
		cur = cur->next;
	}

	return NULL;
}

/**
 * Compute a^{1/p} where a is the input polynomial and
 * p is Pptr->prime. 
*/
void computePolyToPowerPInv_spX (const duspoly_t* a, duspoly_t** ap, const Prime_ptr* Pptr);
void computePolyToPowerPInv_spX_inp (duspoly_t** a, const Prime_ptr* Pptr);

/**
 * Square Free Factorization 
 * Given a monic polynomial, a, and return its square free decomposition.
 */
void squareFreeFactorizationInForm_spX (const duspoly_t* a, factors_t** f, const Prime_ptr* Pptr);
void squareFreeFactorizationInFormll_spX (const duspoly_t* a, factsll_t** F, const Prime_ptr* Pptr);

/** 
 * Compute h^{p} and h^{p} - x mod f  
 * If h is NULL polynomial, this function computes h = x^{p} mod f
 * otherwise, h = h^{p} mod f and p is Pptr->prime.
 * It returns h^{p} - x mod f.
 */
duspoly_t* expXModSubInForm_spX (polysize_t q, const duspoly_t* f, duspoly_t** h, const Prime_ptr* Pptr);

/** 
 * Distinct Degree Factoriztion 
 * Given a monic square-free polynomial, f, and return distinct-degree factors 
 */
void distinctDegFactorizationInFormll_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr);

/**
* Mod Powering specialized for EqualDegreeFactorization 
*
* A faster way to compute a**((p**d - 1)/2) mod f where p is the field characteristic. 
* we use this fact that (p**d - 1)/2 = (p**(d-1) + ... + p + 1)*(q-1)/2
*
* @param a is a polynomial
* @param d is a bound for computation 
* @param f is modulus polynomial 
* @return mul_pows is the result which can be used directly in equalDegreeFactorization. 
*
*/
void modPoweringForEDF_spX (const duspoly_t* a, polysize_t d, const duspoly_t* f, duspoly_t** mul_pows, const Prime_ptr* Pptr);

/** 
 * Equal Degree Splitting 
 * Given monic polynomial, f of degree n where d|n and there is 
 * a set of irreducable polynomials g_1, ..., g_s so that g = g_1 *...* g_s | f ,
 * this function computes g with following a probabilistic approach.
 * If this algorithm couldn't find a proper g, it returns 0, and you need to call it again.
 * Note the results might be different at each call w.r.t the randomPolynomialInForm_spX 
 */
int equalDegSplittingInForm_spX (const duspoly_t* f, polysize_t d, duspoly_t** g, const Prime_ptr* Pptr);

/** 
 * Equal Degree Factorization 
 * Given a monic square-free polynomial, f, of degree n where d|n and it has 
 * irreducable factors of degree d. Return all of equal-degree factors of f 
 * by calling the equalDegSplittingInForm_spX algorithm and in a divide and conquer manner.
 */
void equalDegFactorizationInFormll_spX (const duspoly_t* f, polysize_t d, factsll_t** G, const Prime_ptr* Pptr);

/**
 * Equal Degree Factorization for a list of monic square-free polynomials, D.
 * And G is a list of all factors which is concatinated.
 * Note the degree of irreducable factors for each polynomial in D, has to be stored 
 * in 'exp' member of factsll_t data-type.
 */
void equalDegFactorsInFormll_spX (factsll_t* D, factsll_t** G, const Prime_ptr* Pptr);

/**
 * Modular Factorization in Montogmery form
 * Given a polynomial f, Return all of factors
 * Note this algorithm is using a better approach for bonding the square-free and
 * distinct-degree factorization with equal-degree factorization algorithm. 
 * See M. Asadi slides in "Make BPAS Great" presentations
 */
void modFactorizationInFormll_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr);

/**
 * Verify the correctness of modFactorizationInFormll_spX
 */
int modFactorVerificationInFormll_spX (const duspoly_t* f, factsll_t* G, const Prime_ptr* Pptr);

/**
 * Convert the output of the modFactorizationInFormll_spX (factsll_t) to the list of polynomials.
 * see factors_t in DUSP_Support.h 
 */ 
factors_t* convertToFactorsList_spX (factsll_t* F, const Prime_ptr* Pptr);


#ifdef __cplusplus
}
#endif

#endif // DUSP_SUPPORT_FACTORING_H




// // Temporary Factoring Interface:
// void squareFreeFactorizationInFormll_spX (const duspoly_t* a, factsll_t** F, const Prime_ptr* Pptr) {
// 	squareFreeFactorizationInFormll1_spX (a, F, Pptr); // TODO: after fixing DUSP SSMul and SSDiv
// }

// void exponentiatePolynomialInForm_spX (const duspoly_t* a, polysize_t n, duspoly_t** an, const Prime_ptr* Pptr) {
// 	exponentiatePolynomialInForm1_spX (a,n,an,Pptr); // TODO: after fixing DUSP SSMul
// }

// void distinctDegFactorizationInFormll_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr) {
// 	distinctDegFactorizationInFormll1_spX (f,G,Pptr);
// }

// void equalDegFactorizationInFormll_spX (const duspoly_t* f, polysize_t d, factsll_t** G, const Prime_ptr* Pptr) {
// 	equalDegFactorizationInFormll1_spX (f,d,G,Pptr);
// }

// void equalDegFactorsInFormll_spX (factsll_t* D, factsll_t** G, const Prime_ptr* Pptr) {
// 	equalDegFactorsInFormll1_spX (D,G,Pptr);
// }

// void modFactorizationInFormll_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr) {
// 	modFactorizationInFormll1_spX (f,G,Pptr);
// }
