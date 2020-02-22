
/*****
 * Supporting methods for sparse multivariate ratioanl polynomials written in 
 * pure C. 
 *
 * Throughout this it is assumed that polynomials are compatiable. They must have
 * the same number of variables and the same variable ordering. 
 *
 * Lexicographical ordering is used for term ordering throughout. 
 *
 * Functions assume that the number of variables is at least 1. But, in the most
 * basic case, a NULL exponent vector encodes nvar = 0. Similar to how Node* = NULL
 * encodes the zero polynomial.
 *****/

#ifndef _SMZP_SUPPORT_H_
#define _SMZP_SUPPORT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>

#include "../RationalNumberPolynomial/ExponentPackingDefines.h"
#include "../RationalNumberPolynomial/SMQP_Support-AA.h"



/*****************
 * Alternating Array definition and helpers
 *****************/

typedef struct AAZElem {
	mpz_t coef;
	degrees_t degs;
} AAZElem_t;

typedef struct AltArrZ {
	int size;
	int alloc;
	int nvar;
	AAZElem_t* elems;
} AltArrZ_t;

#define AA_SIZE(A) ((A)->size)
#define AA_ALLOC(A) ((A)->alloc)

static inline void freePolynomial_AAZ(AltArrZ_t* aa) {
	if (aa != NULL) {
		int size = aa->size;
		AAZElem_t* elems = aa->elems;
		for (int i = 0; i < size; ++i) {
			mpz_clear(elems[i].coef);
		}
		free(aa->elems);
		free(aa);
	}
}

typedef struct AAZElem_DegList {
	mpz_t coef;
	degree_t* degs;
} AAZElem_DegList_t;

typedef struct AltArrZDegList {
	int size;
	int alloc;
	int nvar;
	AAZElem_DegList_t* elems;
} AltArrZDegList_t;

static inline void freePolynomial_AAZDL(AltArrZDegList_t* aa) {
	if (aa != NULL) {
		int size = aa->size;
		AAZElem_DegList_t* elems = aa->elems;
		for (int i = 0; i < size; ++i) {
			mpz_clear(elems[i].coef);
			free(elems[i].degs);
		}
		free(aa->elems);
		free(aa);
	}
}


/**
 * Create a new polynomial alternating array with a specified allocation size.
 * The array is not initialized.
 */
static inline AltArrZ_t* makePolynomial_AAZ(int allocSize, int nvar) {
	if (allocSize < 1) {
		return NULL;
	}
	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	newAA->size = 0; 
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*allocSize);
	return newAA;
}

static inline AltArrZDegList_t* makePolynomial_AAZDL(int allocSize, int nvar) {
	if (allocSize < 1) {
		return NULL;
	}
	AltArrZDegList_t* newAA = (AltArrZDegList_t*) malloc(sizeof(AltArrZDegList_t));
	newAA->size = 0; 
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->elems = (AAZElem_DegList_t*) malloc(sizeof(AAZElem_DegList_t)*allocSize);
	return newAA;	
}

static inline AltArrZ_t* makeConstPolynomial_AAZ(int allocSize, int nvar, const mpz_t coef) {
	if (allocSize < 1) {
		return NULL;
	}

	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	newAA->size = 1;
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*allocSize);
	mpz_init(newAA->elems->coef);
	mpz_set(newAA->elems->coef, coef);
	newAA->elems->degs = 0;
	return newAA;
}

/**
 * Re-pack the exponent vectors of aa so that they have newNvar number of variables.
 * Expansion is done to the right. That is, exponents are packed into indices 
 * in the new exponent vector starting from 0.
 */
void expandNumVars_AAZ(AltArrZ_t* aa, int newNvar);

/**
 * Determine if two polynomials are exactly equal. 
 * This function does NOT check whether two polynomials are mathematically equal
 * for instance, if they have a different number of variables. The polynomials
 * must match exactly.
 *
 * returns 1 iff they match exactly.
 */
int isExactlyEqual_AAZ(AltArrZ_t* a, AltArrZ_t* b);
/**
 * Re-pack the exponent vectors of aa so that they have newNvar number of variables.
 * Expansion is done to the left. That is, exponents are packed so that the the 
 * leading exponents in the new exponent vector are 0.
 */
void expandNumVarsLeft_AAZ(AltArrZ_t* aa, int newNvar);

/**
 * Re-pack the exponent vectors of aa so they have one less variable. 
 * This re-packing is done such that the exponent at index idx is discarded completely. 
 * NOTE sorting may be needed after a call to this function, depending on the circumstances.
 * No combination of like-terms or re-ordering is done here.
 */
void shrinkNumVarsAtIdx_AAZ(AltArrZ_t* aa, int idx);

/**
 * Like reorderVars except that the varMap can map to indices less than 0. That is, 
 * if varMap[i] < 0 then variable at index i is removed from the exponent vector.
 *
 * It is assumed that the polynomial aa has zero exponents for all such i. 
 */
void shrinkAndReorderVars_AAZ(AltArrZ_t* aa, int* varMap, int varmapSize);

/**
 * Given a mapping between current exponent indices and new ones (such that
 * exponent at index i is moved to index varMap[i]), reorder the exponent
 * vectors of aa.
 *
 * For variables of index >= varMapSize, it is assumed they are all zero and so
 * their mapped position is not important.
 */
void reorderVars_AAZ(AltArrZ_t* aa, int* varMap, int varMapSize);

degrees_t calculateMaxDegs_AAZ(AltArrZ_t* aa);

static inline void resizePolynomial_AAZ(AltArrZ_t* aa, int allocSize) {
	if (allocSize < aa->size) {
		for (int i = aa->size; i >= allocSize; --i) {
			mpz_clear(aa->elems[i].coef);
		}
		aa->size = allocSize;
	}
	aa->elems = (AAZElem_t*) realloc(aa->elems, sizeof(AAZElem_t)*allocSize);
	aa->alloc = allocSize;
}

/**
 * Construct an alternating array polynomial representation from a Node rep.
 *
 * returns the new alternating array pointer. 
 */
AltArrZ_t* deepCopyPolynomial_AAZFromNode(Node* a, int nvar);

/**
 * Construct an linked list polynomial representation from an alternating array rep.
 *
 * returns the head node pointer of the linked list. 
 */
Node* deepCopyPolynomial_NodeFromAAZ(AltArrZ_t* aa);

/**
 * Construct an alternating array with un-packed exponents from an alternating
 * array with packed exponents. 
 */
AltArrZDegList_t* deepCopyPolynomial_AAZDegListFromAA(AltArrZ_t* aa);

/** 
 * Create a deep copy of the polynomial represented by the alternating array aa.
 * 
 * returns a pointer to the new alternating array.
 */
AltArrZ_t* deepCopyPolynomial_AAZ(AltArrZ_t* aa);

/**
 * Create a deep copy of the rational number polynomial aa.
 *
 * returns a pointer to the new integer alternating array.
 */
AltArrZ_t* deepCopyPolynomial_AAZFromAA(AltArr_t* aa);

/**
 * Create a deep copy of this integer polynomial and return it as a rational
 * number polynomial.
 *
 * returns a pointer to the new rational number alternating array. 
 */
AltArr_t* deepCopyPolynomial_AAFromAAZ(AltArrZ_t* aa);

/**
 * Sort a polynomial in place represented by an alternating array.
 *
 * returns the same pointer given.
 */
AltArrZ_t* sortPolynomial_AAZ(AltArrZ_t* aa);

/**
 * Sort a polynomial in place using merge sort.
 */
void mergeSortPolynomial_AAZ(AltArrZ_t* aa);

/**
 * Given a polynomial in sorted order but with possible monomial duplicates, 
 * combine like terms and condense the polynomial.
 */
void condensePolyomial_AAZ(AltArrZ_t* aa);

/**
 * Negate a polynomial in place.
 */
void negatePolynomial_AAZ(AltArrZ_t* a);

/**
 * Evaluate a polynomial represented as an alternating array.
 * This method returns another polynomial as not all variables need
 * to have values supplied. But, if they do, a constant term will be returned. 
 * both active and vals are arrays of size nvar. active[i] determines if 
 * the variable at degs[i] is to be evaluated using vals[i] as value.
 */
AltArrZ_t* evaluatePoly_AAZ(AltArrZ_t* aa, int* active, mpz_t* vals, int nvar);


AltArrZ_t* convertFromAAZElemToAAZ (AAZElem_t* coef, int coefSize, int nvar);


AltArrZ_t* swappingExponents_AAZ (AltArrZ_t* aa, int idx1, int idx2);

/**
 * (right) shift main variable of polynomial aa of size n.
 */ 
AltArrZ_t* mainLShiftPolynomial_AAZ (AltArrZ_t* aa, int n);
AltArrZ_t* mainLShiftPolynomial_AAZ_inp (AltArrZ_t* aa, int n);   
/**
 * Given polynomial a, return the leading term
 */
AltArrZ_t* leadingTerm_AAZ (AltArrZ_t* aa, int nvar);

/**
 * Get the leading variable, that is, the highest-order variable with positive degree
 * of this polynomial.
 * returns the postion of leading variable, -1 (for constant polynomials) and -2 (for zero polynomials).
 */
int leadingVariable_AAZ (AltArrZ_t* aa);

/**
 * Main Leading Degree
 * given a sorted multivariable polynomial aa, this function returns 
 * the maximum value of the main variable of aa.
 */ 
int mainLeadingDegree_AAZ (AltArrZ_t* aa);
    
/**
 * Main Leading Coefficient
 * given a multivariable polynomial aa, this function returns 
 * the leading coefficient of aa w.r.t the main variable.
 * NOTE: This function is NOT inplace! 
 */ 
AltArrZ_t* mainLeadingCoefficient_AAZ (AltArrZ_t* aa);

/**
 * Main Leading Coefficient w.r.t the e-th variable
 * given a multivariable polynomial aa and index of a variable,
 *  this function returns the leading coefficient of aa w.r.t 
 * the special variable.
 * NOTE: This function is NOT inplace! 
 */ 
AltArrZ_t* mainCoefficientAtIdx_AAZ (AltArrZ_t* aa, int e);

/**
 * given two polynomials, return polynomial with bigger leading term
 * Note: return polynomial is a deep copy of the polynomial
 */
AltArrZ_t* maxPolynomials_AAZ (AltArrZ_t* a, AltArrZ_t* b);

/**
 * given two polynomials, return polynomial with bigger leading term
 * Note: Is is done INPLACE w.r.t. the first input, a.
 */
AltArrZ_t* maxPolynomials_AAZ_inp (AltArrZ_t* a, AltArrZ_t* b);

    
/*****************
 * SMQP Addition
 *****************/

/**
 * Add two polynomials given their alternating arrays, a and b.
 * nvar: number of variables in the polynomials.
 * returns a pointer to the sum alternating array.
 */
AltArrZ_t* addPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar);

/**
 * Add two polynomials given their alternating arrays, a and b.
 * This is done in-place wrt a. That is, using as much memory of a
 * in the process of building the sum. The pointer a is invalidated
 * in this process and the return of this function should be reassigned to
 * the original a.
 */
AltArrZ_t* addPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar);

/**
 * Subtract two polynomials given their alternating arrays, a and b.
 * nvar: number of variables in the polynomials.
 * returns a pointer to the difference alternating array.
 */
AltArrZ_t* subPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar);

/**
 * Subtract two polynomials given their alternating arrays, a and b.
 * This is done in-place wrt a. That is, using as much memory of a
 * in the process of building the difference. The pointer a is invalidated
 * in this process and the return of this function should be reassigned to
 * the original a.
 */
AltArrZ_t* subPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar);

/**
 * the following subPolynomial is used in DucosSubresultantChainZ 
 * implemented in SMZP_Support_Recursive.c
 */
AltArrZ_t* CFDucosOptZ_subPolynomials_AAZ_inp (AltArrZ_t* a, AltArrZ_t* b, int nvar,
					       AltArrZ_t** Pe, int e);


/*****************
 * SMQP Multiplication & Helpers
 *****************/

typedef struct ProductHeapChain_AAZ {
#if SMQP_INT_PRODHEAP
	int a_i;
	int b;
#else
	AAZElem_t* a_i;
	AAZElem_t* b;
#endif
	struct ProductHeapChain_AAZ* next;
} ProductHeapChain_AAZ;

typedef struct ProductHeapElem_AAZ {
	degrees_t degs;
	ProductHeapChain_AAZ* chain;
} ProductHeapElem_AAZ;

typedef struct {
	ProductHeapElem_AAZ* elements;
	polysize_t heapSize;
	polysize_t maxHeapSize;
	int nvar;
#if SMQP_INT_PRODHEAP
	int lastB;
#else
	AAZElem_t* lastB;
#endif
} ProductHeap_AAZ;


/**
 * Make an element for the product heap, combining a and b as the 
 * element's product.
 */
#if SMQP_INT_PRODHEAP
static inline ProductHeapChain_AAZ* prodheapMakeChain_AAZ(int a, int b, ProductHeapChain_AAZ* next) {
#else
static inline ProductHeapChain_AAZ* prodheapMakeChain_AAZ(AAZElem_t* a, AAZElem_t* b, ProductHeapChain_AAZ* next) {
#endif
	ProductHeapChain_AAZ* chain = (ProductHeapChain_AAZ*) malloc(sizeof(ProductHeapChain_AAZ));
	chain->a_i = a;
	chain->b = b;
	chain->next = next;
	return chain;
}

/**
 * Free an entire product heap chain.
 */
static inline void prodheapFreeChain_AAZ(ProductHeapChain_AAZ* chain) {
	ProductHeapChain_AAZ* next;
	while(chain != NULL) {
		next = chain->next;
		free(chain);
		chain = next;
	}
}

/**
 * Cleanup memory initialized for the heap
 */
static inline void prodheapFree_AAZ(ProductHeap_AAZ* h) {
	ProductHeapElem_AAZ* elems = h->elements;
	polysize_t s = h->heapSize;
	for (polysize_t i = 0; i < s; ++i) {
		prodheapFreeChain_AAZ(elems[i].chain);
	}
	free(h->elements);
	free(h); 
}

/**
 * Create an empty product heap. 
 */
static inline ProductHeap_AAZ* prodheapCreate_AAZ(int nvar) {
	ProductHeap_AAZ* h = (ProductHeap_AAZ*) malloc(sizeof(ProductHeap_AAZ));
	h->elements = NULL;
	h->heapSize = 0;
	h->maxHeapSize = 0;
	h->nvar = nvar;
	return h;
}

/**
 * Initialize the product heap with the two polynomials to multiply, a and b.
 *
 * We know the maximum heap size is numTerms(a) as at most one of a_i*b
 * is in the heap at once. 
 */
ProductHeap_AAZ* prodheapInit_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar);

/**
 * Increase the capacity of the heap to newAllocSize. 
 */
static inline void prodheapResize_AAZ(ProductHeap_AAZ* h, int newAllocSize) {
	h->elements = (ProductHeapElem_AAZ*) realloc(h->elements, sizeof(ProductHeapElem_AAZ)*newAllocSize);
	h->maxHeapSize  = newAllocSize;
}

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
#if SMQP_INT_PRODHEAP
void prodheapInsert_AAZ(ProductHeap_AAZ* h, ProductHeapChain_AAZ* chain, register degrees_t degs);
#else
void prodheapInsert_AAZ(ProductHeap_AAZ* h, ProductHeapChain_AAZ* elem);
#endif

/**
 * Peak into the heap and get the exponent vector of the product of the max elem.
 * NOTE that a non-NULL pointer is invalidated after calling prodheapRemoveMax on h.
 */
static inline degrees_t* prodheapPeek_AAZ(ProductHeap_AAZ* h) {
	if (h->heapSize > 0) {
		return &(h->elements->degs);
	}
	return NULL;
}

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * NOTE it is only valid ot call RemoveMax if prodheapPeek_AAZ returns non-NULL.
 */
ProductHeapChain_AAZ* prodheapRemoveMax_AAZ(ProductHeap_AAZ* h);

/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure. 
 * It is assumed that both a and b have compatible exponent vectors.
 * 
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
AltArrZ_t* multiplyPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar);

/**
 * Multiply two polynomials in-place wrt the first polynomial.
 * That is, the multiplication occurs, reusing elements of a as much as possible
 * and the previous content of a is destroyed in the process.
 * It is assumed that a and b have the same sized exponent vectors. 
 *
 *
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the new head node of a.
 *
 */
AltArrZ_t* multiplyPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar);



/*****************
 * Polynomial exponentiation 
 *****************/


/**
 * Given a polynomial, a, compute a^n.
 *
 * n: a positive integer
 * nvar: the number of variables of a.
 *
 */
AltArrZ_t* exponentiatePoly_AAZ(AltArrZ_t* a, unsigned int n, int nvar);



/*****************
 * Polynomial division
 *****************/


/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree_t
 */
#if SMQP_INT_PRODHEAP 
void divisionGetNextTerm_AAZ(ProductHeap_AAZ* h, const AAZElem_t* __restrict__ aElems, const AAZElem_t* __restrict__ bElems, mpz_t* retCoef);
#else 
void divisionGetNextTerm_AAZ(ProductHeap_AAZ* h, mpz_t* retCoef);
#endif

/** 
 * Given a polynomial, c, and a term, b, determine polynomials a and r
 * such that c = b*a + r.
 * a and r are returned in res_a and res_r, respectively. 
 */
void divideBySingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar);

/** 
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r 
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar);

 void exactDividePolynomials_AAZ (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, register int nvar);
/**
 * Pseudo division for univariate polynomials.
 */
void univariatePseudoDividePolynomials_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy);

/**
 * Determine if a polynomial, c, is exactly divisble by another, b.
 * If the division is exact, returns the quoteint in res_a. 
 * Otherwise, the value of res_a is undefined.
 * returns 1 iff division is exact. 0 otherwise.
 */
int divideTest_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar);



/*****************
 * Derivative / Integral
 *****************/

/**
 * Get the k-th partial derivative of aa for the variable at index idx.
 *
 * returns the partial derivative
 */
AltArrZ_t* derivative_AAZ(AltArrZ_t* aa, int idx, int k);


/**
 * Get the k-th partial integral of aa for the variable at index idx;
 *
 * Note: constants of integration are not included.
 *
 * returns the partial integral.
 */
AltArr_t* integral_AAZ(AltArrZ_t* aa, int idx, int k);



/*****************
 * Content, PrimitivePart, etc.
 *****************/

/**
 * Get the "integral" content of a polynomial.
 * That is, the rational number content whose corresponding
 * primitive part is a primitive integer polynomial.
 *
 * returns the content through ret parameter.
 */
void integralContent_AAZ(AltArrZ_t* aa, mpz_t ret);

/**
 * Get the primitive part, a primitive integer polynomial, of aa.
 *
 * returns the primitive part.
 */
AltArrZ_t* primitivePart_AAZ(AltArrZ_t* aa);

/**
 * Get the primitive part, a primitive integer polynomial, of aa.
 * Returns the integral content in cont.
 *
 * returns the primitive part.
 */
AltArrZ_t* primitivePartAndContent_AAZ(AltArrZ_t* aa, mpz_t cont);

/**
 * Convert the input aa to its primitive part in place.
 */
void primitivePart_AAZ_inp(AltArrZ_t* aa);

/**
 * Given to polynomials, both with 1 variable, get their GCD as a primitive polynomial.
 *
 * returns the GCD.
 */
AltArrZ_t* univariateGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b);

/**
 * Test for if the polynomial aa is actually an integer polynomial.
 * If so, returns the integral content in mpzG. Otherwise, mpzG is set to 0.
 */
void integerPolynomialTest_AAZ(AltArrZ_t* aa, mpz_t mpzG);

/**
 * Given a polynomial, find a monic factor common to all terms of the polynomial.
 * If factored is not NULL then the input polynomial with common factor removed
 * is returned in factored.
 *
 * returns the common factor. 
 */
AltArrZ_t* commonFactor_AAZ(AltArrZ_t* a, AltArrZ_t** factored);



/*****************
 * Interpolation / Evaluation
 *****************/

/**
 * Evaluate a univariate polynomial at the point, point. The result is 
 * returned as a multi-precision rational number, res.
 */
void univarEvaluate_AAZ(AltArrZ_t* aa, const mpz_t point, mpz_t res);



/****************
* Multi-Divisor Division
*****************/

/** Multi-Divisor Division (MDD) using Heap (Johnson's) Division
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void heapMDD_AAZ(AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);

/** Multi-Divisor Division (MDD) for Triangular Sets
 * @param[in] f The dividend
 * @param[in] G The triangular-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void triangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);

/** Multi-Divisor Division (MDD) for Triangular Sets using Primitive Factorization
 * @param[in] f The dividend
 * @param[in] G The triangular-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void primitiveFactorTriangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);

/** Multi-Divisor Division (MDD)
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[in] type The type of MDD (type = 0 ? HeapMDD : (type = 1 ? triangularSetMDD : primitiveFactorTriangularSetMDD))
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void multiDivisorDivision_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar, int type);

int multiDivisorDivisionVerification_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t* r, AltArrZ_t* hPow, int nSet, int tnvar);



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


