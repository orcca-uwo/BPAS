
/*****
 * Supporting methods for sparse univariate rational polynomials written in 
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

#ifndef _SUQP_SUPPORT_H_
#define _SUQP_SUPPORT_H_

#ifdef __cplusplus
extern "C" { 
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>

#include "SMQP_Support-AA.h"
#include <math.h>
#include "../../include/IntegerPolynomial/primes.h"

#include "FiniteFields/SmallPrimeField_Support.h"


/*****************
 * Data types used throughout 
 *****************/

#define compareExponentVectors(A, B) ((A) < (B) ? -1 : ((A) > (B) ? 1 : 0) )

/*****************
 * Alternating Array definition and helpers
 *****************/
// the first Array is the term element and the second is the polymonial:  
 typedef struct AAElemU {
	ratNum_t coef;
	degree_t deg;
} AAElemU_t;

// the second one is the polynomial as we mentioned
typedef struct AltArrU {
	                          // size is number of terms and alloc is space that we totally use
	int size;
	int alloc;
	AAElemU_t* elems;
} AltArrU_t;

#define AA_SIZE(A) ((A)->size)
#define AA_ALLOC(A) ((A)->alloc)

static inline void freePolynomial_AAU(AltArrU_t* aa) {
	if (aa != NULL) {
		int size = aa->size;
		AAElemU_t* elems = aa->elems;
		for (int i = 0; i < size; ++i) {
			mpq_clear(elems[i].coef);
		}
		free(aa->elems);
		free(aa);
	}
}

/**
 * Create a new polynomial alternating array with a specified allocation size.
 * The array is not initialized.
 */


static inline AltArrU_t* makePolynomial_AAU(int allocSize) {
	if (allocSize < 1) {
		return NULL;
	}
	AltArrU_t* newAA = (AltArrU_t*) malloc(sizeof(AltArrU_t));
	newAA->size = 0; 
	newAA->alloc = allocSize;
	newAA->elems = (AAElemU_t*) malloc(sizeof(AAElemU_t)*allocSize);
	return newAA;
}

static inline AltArrU_t* makeConstPolynomial_AAU(int allocSize,  const mpq_t coef) {
	if (allocSize < 1) {
		return NULL;
	}

	AltArrU_t* newAA = (AltArrU_t*) malloc(sizeof(AltArrU_t));
	newAA->size = 1;
	newAA->alloc = allocSize;
	newAA->elems = (AAElemU_t*) malloc(sizeof(AAElemU_t)*allocSize);
	mpq_init(newAA->elems->coef);
	mpq_set(newAA->elems->coef, coef);
	newAA->elems->deg = 0;
	return newAA;
}


static inline void resizePolynomial_AAU(AltArrU_t* aa, int allocSize) {
	if (allocSize < aa->size) {
		for (int i = aa->size; i >= allocSize; --i) {
			mpq_clear(aa->elems[i].coef);
		}
		aa->size = allocSize;
	}
	aa->elems = (AAElemU_t*) realloc(aa->elems, sizeof(AAElemU_t)*allocSize);
	aa->alloc = allocSize;
}

/** 
 * Create a deep copy of the polynomial represented by the alternating array aa.
 * 
 * returns a pointer to the new alternating array.
 */
AltArrU_t* deepCopyPolynomial_AAU(AltArrU_t* aa);

// print polynomial

void printAAU(AltArrU_t* aa);

/**
 * Sort a polynomial in place represented by an alternating array.
 *
 * returns the same pointer given.
 */
AltArrU_t* sortPolynomial_AAU(AltArrU_t* aa);




void condensePolyomial_AAU(AltArrU_t* aa);
/**
 * Add a term to the end of the alternating array with degrees d and coefficient coef.
 * This method does not reallocate the array if the term to add exceeds allocated space.
 * In this case, behaviour is undefined and erroneous. 
 */
static inline void addTerm_AAU(AltArrU_t* aa, degree_t d, const ratNum_t coef) {
	mpq_init(aa->elems[aa->size].coef);
	mpq_set(aa->elems[aa->size].coef, coef);
	aa->elems[aa->size].deg = d;
	++(aa->size);
}


/**
 * Add a term to the end of this polynomial, swapping in the supplied coefficient 
 * for effective memory usage. The passed-in coef becomes 0/1.
 * This method is unsafe as addTerm_AA.
 */
static inline void addTermSwap_AAU(AltArrU_t* aa, degree_t d, ratNum_t coef) {
	mpq_init(coef);
	mpq_swap(aa->elems[aa->size].coef, coef);
	aa->elems[aa->size].deg = d;
	++(aa->size);
}

/**
 * Add a term to the end of this polynomial, reallocating more space as needed.
 */
void addTermSafe_AAU(AltArrU_t* aa, degree_t d, const ratNum_t coef);

/**
 * Multiply two polynomial terms together to create another polynomial term.
 */
AAElemU_t* multiplyTerms_AAU(AAElemU_t* a, AAElemU_t* b);
/**
 * Negate a polynomial in place.
 */

void negatePolynomial_AAU(AltArrU_t* aa);
/*
 * Add two polynomials given their alternating arrays, a and b.
 * nvar: number of variables in the polynomials.
 * returns a pointer to the sum alternating array.
 */
 
 AltArrU_t* addPolynomials_AAU( AltArrU_t* a,  AltArrU_t* b);


 AltArrU_t* subPolynomials_AAU( AltArrU_t* a,  AltArrU_t* b);



/**
 * Add polynomial a and b, putting the sum back into a.=inplace=inp
 * 
 * nvar: size of exponent vectors in a and b.
 *
 * returns a pointer to the new head of a.
 */
AltArrU_t* addPolynomials_AAU_inp(AltArrU_t* a, AltArrU_t* b);
AltArrU_t* subPolynomials_AAU_inp(AltArrU_t* a, AltArrU_t* b);




/*****************
 * SMQP Multiplication & Helpers
 *****************/

/**
 * Data element for the ProductHeap data structure
 */

typedef struct ProductHeapChain_AAU {

	int a_i;
	int b;

	struct ProductHeapChain_AAU* next;
} ProductHeapChain_AAU;
//////////////////////////////////////elements heap

typedef struct ProductHeapElem_AAU {
	degree_t deg;
	ProductHeapChain_AAU* chain;
} ProductHeapElem_AAU;

////////////   answer=C

typedef struct {
	ProductHeapElem_AAU* elements;
	polysize_t heapSize;
	polysize_t maxHeapSize;
	
// size(multiplication a*b^-=a/b)
	int lastB;

} ProductHeap_AAU;
////////////////////////////////////
void prodheapPrint_AAU(ProductHeap_AAU* h);

 //////////////////////////////////////////////////////////
/**
 * Initialize the product heap with the two polynomials to multiply, a and b.
 *
 * We know the maximum heap size is numTerms(a) as at most one of a_i*b
 * is in the heap at once. 
 */
/**
 * Cleanup memory initialized for the heap
 */


static inline void prodheapFreeChain_AAU(ProductHeapChain_AAU* chain) {
	ProductHeapChain_AAU* next;
	while(chain != NULL) {
		next = chain->next;
		free(chain);
		chain = next;
	}
}

/**
 * Cleanup memory initialized for the heap
 *//**
 * Make an element for the product heap, combining nodes a and b as the 
 * element's product.
 */
static inline void prodheapFree_AAU(ProductHeap_AAU* h) {
	ProductHeapElem_AAU* elems = h->elements;
	polysize_t s = h->heapSize;
	for (polysize_t i = 0; i < s; ++i) {
		prodheapFreeChain_AAU(elems[i].chain);
	}
	free(h->elements);
	free(h); 
}

/**
 * Create an empty product heap. 
 */
static inline ProductHeap_AAU* prodheapCreate_AAU() {
	ProductHeap_AAU* h = (ProductHeap_AAU*) malloc(sizeof(ProductHeap_AAU));
	h->elements = NULL;
	h->heapSize = 0;
	h->maxHeapSize = 0;
	return h;
}

/**
 * Make an chain for the product heap.
 * 
 */
static inline ProductHeapChain_AAU* prodheapMakeChain_AAU(int a, int b, ProductHeapChain_AAU* next) {

	ProductHeapChain_AAU* chain = (ProductHeapChain_AAU*) malloc(sizeof(ProductHeapChain_AAU));
	chain->a_i = a;
	chain->b = b;
	chain->next = next;
	return chain;
}
 ProductHeap_AAU* prodheapInit_AAU(AltArrU_t* a, AltArrU_t* b);

/**
 * Increase the capacity of the heap to newAllocSize. 
 */
static inline void prodheapResize_AAU(ProductHeap_AAU* h, int newAllocSize) {
	h->elements = (ProductHeapElem_AAU*) realloc(h->elements, sizeof(ProductHeapElem_AAU)*newAllocSize);
	h->maxHeapSize  = newAllocSize;
}

void prodheapInsert_AAU(ProductHeap_AAU* h, ProductHeapChain_AAU* chain, register degree_t deg) ;

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
 /**
 * Peak into the heap and get the exponent vector of the product of the max elem.
 * NOTE that a non-NULL pointer is invalidated after calling prodheapRemoveMax on h.
 */
static inline degree_t* prodheapPeek_AAU(ProductHeap_AAU* h) {
	if (h->heapSize > 0) {
		return &(h->elements->deg);
	}
	return NULL;
}




ProductHeapChain_AAU* prodheapRemoveMax_AAU(ProductHeap_AAU* h) ;

AltArrU_t* multiplyPolynomials_AAU(AltArrU_t* __restrict__ a, AltArrU_t* __restrict__ b);
/**
 * Multiply two polynomials in-place wrt the first polynomial.
 * That is, the multiplication occurs, reusing elements of a as much as possible
 * and the previous content of a is destroyed in the process.
 * It is assumed that a and b have the same sized exponent vectors. 
 *
 *
 *  
 * returns a pointer to the new head node of a.
 *
 */
AltArrU_t* multiplyPolynomials_AAU_inp(AltArrU_t* a, AltArrU_t* b) ;

/*****************
 * Polynomial division
 *****************/
/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients 
 * 
 */
static inline int monomialDivideTestU(degree_t adegs, degree_t bdegs) {
	return (adegs>=bdegs);
}
/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients 
 * 
 */

void divisionGetNextTerm_AAU(ProductHeap_AAU* h, const AAElemU_t* __restrict__ aElems, const AAElemU_t* __restrict__ bElems, mpq_t* retCoef) ;
/** 
 * Given a polynomial, c, and a term, b, determine polynomials a and r
 * such that c = b*a + r.
 * a and r are returned in res_a and res_r, respectively. 
 */
void divideBySingleTerm_AAU(AltArrU_t* c, AltArrU_t* b, AltArrU_t** res_a, AltArrU_t** res_r);

/** 
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r 
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials_AAU(AltArrU_t* c, AltArrU_t* b, AltArrU_t** res_a, AltArrU_t** res_r);
/**
 * Determine if a polynomial, c, is exactly divisble by another, b.
 * If the division is exact, returns the quoteint in res_a. 
 * Otherwise, the value of res_a is undefined.
 * returns 1 iff division is exact. 0 otherwise.
 */
int divideTest_AAU(AltArrU_t* c, AltArrU_t* b, AltArrU_t** res_a);



// exponentiation
AltArrU_t* exponentiatePoly_AAU(AltArrU_t* a, unsigned int n) ;

/*****************
 * Content, PrimitivePart, etc.
 *****************/

void integralContent_AAU(AltArrU_t* aa, mpq_t ret) ;

AltArrU_t* primitivePart_AAU(AltArrU_t* aa);               

/*****************
 * Derivative / Integral
 *****************/

AltArrU_t* derivative_AAU(AltArrU_t* aa, int k);
AltArrU_t* integral_AAU(AltArrU_t* aa, int k);
/*****************
 * evaluation using horner method
 *****************/
void evaluatePoly_AAU(AltArrU_t* aa, const ratNum_t r, ratNum_t mpqRet);
///////////////////////////////////////////// pomopo /////////////////////////////////
AltArrU_t* pomopoU (AltArrU_t* aa, AltArrU_t* c, AltArrU_t* t, AltArrU_t* b);
///////////////////////////////////////////////////////// lazy and pseudo division
void multiplyByRational_AAU_inp(AltArrU_t* aa, const mpq_t z);
void univariatePseudoDividePolynomials_AAU(AltArrU_t* c, AltArrU_t* b, AltArrU_t** res_a, AltArrU_t** res_r, int* e, int lazy) ;
////////////////////////////////////////////////////////////////////////////////////// Iszero
int  isZero(AltArrU_t* a);
///////////////////////////////////////////////////////////////////////////////////// GCD
AltArrU_t* univariateGCD_AAU(AltArrU_t* a, AltArrU_t* b) ;
///////////////////////////////////////////////////////////////////////////////////// Squar free
AltArrU_t** squareFree_AAU(AltArrU_t * A,int* factsize);
////////////////////////////////////////////////////////////////////////////////////  Random poly
AltArrU_t* buildRandomPolyU(time_t seed, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) ;
/////////////////////////////////////////////////////////////////////////////// ConvertfromAlttoAr
mpz_t * ConvertfromAlttoAr(AltArrU_t* F, int * deg);
AltArrU_t * ConvertfromArtoAlt(mpz_t *Ar,int deg);
/////////////////////////////////////////////////////////////////////////// modularGCD

 AltArrU_t*  modularGCD (AltArrU_t *F,AltArrU_t *G);
//////////////////////////////////////////////////////////  arrayDegree

int arrayDegree(long long int *gp, int curdg);
///////////////////////////////////////////////////////////////////isDivide

int  isDivide(AltArrU_t* F,AltArrU_t*G,AltArrU_t* ret);
//////////////////////////////////////////////////////////////////////content
void content_U(mpz_t * h,mpz_t *ret, int deg) ;

long long int nextprime(int *s, mpz_t m);

void monicGCD_U(long long int *gp, int *gd, long long int *fp, int fd, long long int p, long long int lc);
static  sfixn MulMod_U(sfixn a, sfixn b, sfixn n);
static  sfixn MulMod_U(sfixn a, sfixn b, sfixn n);
static sfixn SubMod_U(sfixn a, sfixn b, sfixn p);
 static sfixn inverseMod_U(sfixn n, sfixn p);
////////////////////////////////////////////////   


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


