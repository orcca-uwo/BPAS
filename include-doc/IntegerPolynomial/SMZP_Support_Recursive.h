
#ifndef _SMZP_SUPPORT_RECURSIVE_H_
#define _SMZP_SUPPORT_RECURSIVE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "SMZP_Support.h"
#include "../RationalNumberPolynomial/SMQP_Support_Recursive-AA.h"

typedef struct {
	int coefSize;
	degree_t exp;
	AAZElem_t* coef;
} RecArrElemZ_t;

typedef struct {
	int alloc;
	int size;
	AltArrZ_t* origAA;
	RecArrElemZ_t* elems;
} RecArrZ_t;

/**
 * Free a RecArrZ_t and it's RecArrElemZ_t array. but not the underlying coef data.
 */
static inline void freeRecArrayZ(RecArrZ_t* rArr) {
	if (rArr->alloc > 0) {
		free(rArr->elems);
	}
	free(rArr);
}

/** 
 * Build a random recursive polynomial.
 * Nterms is the number of terms wrt the main variable. 
 * coefNTerms is the number of terms in the coefficients (when viewed recursively, and nvar > 1);
 * coefBound is the absolute upper limit on the rational numbers in the polynomial.
 * sparsity is the upper limit on the difference between degrees of successive recursive terms. 
 * sparsity is also passed when constructing coefficient polynomials.
 */
RecArrZ_t* buildRandomRecArrZPoly(int nvar, int nterms, int coefNterms, unsigned long int coefBound, degree_t sparsity, int includeNeg);


/***********
 * Conversion to and from recursive view.
 ***********/

/**
 * Convert a polynomial to a recursively viewed polynomial using the variable
 * of highest precedence (index 0) as the main variable.
 *
 * Note this conversion is done INPLACE. And therefore the input Node* is 
 * invalidated by this operation.
 *
 * returns teh recursively viewed polynomial.
 */
RecArrZ_t* convertToRecursiveArrayZ(AltArrZ_t* aa);

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done INPLACE and therefore the original recursively
 * viewed polynomial is then invalidated.
 * 
 * returns the distributed polynomial.
 */
AltArrZ_t* convertFromRecursiveArrayZ(RecArrZ_t* poly, int nvar);

/**
 * Given a recursive polynomial whose head is poly, make a deep copy of it.
 *
 * returns the head to the new polynomial.
 */
RecArrZ_t* deepCopyRecArrayZPolynomial(RecArrZ_t* poly, int nvar);


/**
 * Convert a polynomial to a recursively viewed polynomial given the index
 * within the exponent vector of the variable to become the "main variable"
 * of the recursive view.
 *
 * Note this conversion is not done INPLACE. And therefore the output polynomial is 
 * independent of the input one.
 *
 * returns the recursively viewed polynomial. 
 */
RecArrZ_t* convertToRecursiveArrayZAtIdx (AltArrZ_t* aa, int idx);

/**
 * Convert a recursive viewed polynomial to a AltArr_t* polynomial given the index
 * within the exponent vector of the variable to become the "main variable"
 * of the recursive view.
 *
 * Note this conversion is not done INPLACE. And therefore the output polynomial is 
 * independent of the input one.
 *
 * returns the recursively viewed polynomial. 
 */
AltArrZ_t* convertFromRecursiveArrayZAtIdx (RecArrZ_t* recPoly, int idx, int nvar);


    

/***********
 * Pesudo Division
 ***********/

// Since we are using integer indices, we can use SMQP RecProdHeap definition. 
// The below is just kept for reference. We will use the SMQP methodf.
//

// typedef struct RecProdHeapChain_AA {
// 	int a_i;
// 	int b;
// 	struct RecProdHeapChain_AA* next;
// } RecProdHeapChain_AA;

// typedef struct RecProdHeapElem_AA {
// 	degree_t exp;
// 	RecProdHeapChain_AA* chain;
// } RecProdHeapElem_AA;

// typedef struct {
// 	RecProdHeapElem_AA* elements;
// 	int heapSize;
// 	int maxHeapSize;
// 	int nvar;
// 	int lastB;
// } RecProdHeap_AA;

// /**
//  * Create a new element to be inserted into the heap. 
//  * Takes one term of the multiplied and the entire multiplicand. 
//  *
//  * Returns an element ready to be inserted in the RecProdHeap. 
//  */
// static inline RecProdHeapChain_AA* recProdHeapMakeChain_AA(int a , int b, RecProdHeapChain_AA* next) {
// 	RecProdHeapChain_AA* chain = (RecProdHeapChain_AA*) malloc(sizeof(RecProdHeapChain_AA));
// 	chain->a_i = a;
// 	chain->b = b;
// 	chain->next = next;
// 	return chain;
// }

// static inline void recProdHeapFreeChain_AA(RecProdHeapChain_AA* chain) {
// 	RecProdHeapChain_AA* next;
// 	while(chain != NULL) {
// 		next = chain->next;
// 		free(chain);
// 		chain = next;
// 	}
// }

// /**
//  * Create an empty product heap. 
//  */
// static inline RecProdHeap_AA* recProdHeapCreate_AA(int nvar) {
// 	RecProdHeap_AA* h = (RecProdHeap_AA*) malloc(sizeof(RecProdHeap_AA));
// 	h->elements = NULL;
// 	h->heapSize = 0;
// 	h->maxHeapSize = 0;
// 	h->nvar = nvar;
// 	return h;
// }

// *
//  * Allocate space for a maximum of newAllocSize elements in the heap
  
// static inline void recProdheapResize_AA(RecProdHeap_AA* h, int newAllocSize) {
// 	h->elements = (RecProdHeapElem_AA*) realloc(h->elements, sizeof(RecProdHeapElem_AA)*newAllocSize);
// 	h->maxHeapSize  = newAllocSize;
// }

// static inline void recProdHeapFree_AA(RecProdHeap_AA* h) {
// 	RecProdHeapElem_AA* elems = h->elements;
// 	polysize_t s = h->heapSize;
// 	for (polysize_t i = 0; i < s; ++i) {
// 		recProdHeapFreeChain_AA(elems[i].chain);
// 	}
// 	free(h->elements);
// 	free(h); 
// }

// /**
//  * Insert an RecProdHeap element into the heap h.  
//  */
// void recProdHeapInsert_AA(RecProdHeap_AA* h, RecProdHeapChain_AA* chain, register degree_t insertExp);

// /**
//  * Peak into the heap and get the exponent vector of the product of the max elem.
//  */
// static inline degree_t recProdHeapPeek_AA(RecProdHeap_AA* h) {
// 	if (h->heapSize > 0) {
// 		return h->elements->exp;
// 	}
// 	return -1;
// }

// /**
//  * Remove the maximum chain from the product heap.
//  */
// RecProdHeapChain_AA* recProdHeapRemoveMax_AA(RecProdHeap_AA* h);

/**
 * Get the coefficient of the next term to be produced in the multiplication.
 * Terms are produced in decreasing order by degree.
 * Note that this may return NULL to indicate that term whose degree was 
 * identified as maximum by recProdHeapPeek, ended up being cancelled to 0
 * when combining like terms.
 *
 * returns the coefficient of the next product term, or NULL, as above.
 */
AltArrZ_t* recProdHeapGetNextCoef_AAZ(RecProdHeap_AA* h, const RecArrElemZ_t* __restrict__ aElems, const RecArrElemZ_t* __restrict__ bElems);

/** 
 * Do the pesudo division of c by b. Quotient is returned in res_a, and remainder
 * in res_r. 
 *
 * If e is not NULL then it returns the exact number of division steps which occurred.
 *
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 */
void pesudoDivide_RecArrayZ(RecArrZ_t* c, RecArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, AltArrZ_t** hPow, int nvar, int lazy);


//////////////////////////////////
// Multi-divisor Pseudo-division
//////////////////////////////////

/** 
 * Do the pesudo division of c by b when the idx-th member
 * of exponent vectors is the main variable. 
 * Quotient is returned in res_a, and remainder in res_r. 
 * If e is not NULL then it returns the exact number of division steps which occurred.
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 */
void pesudoDivideAtIdx_AAZ (int idx, AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, AltArrZ_t** hPow, int nvar, int lazy);

void naiveMultiDivisorPseudoDivision_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet);

void multiDivisorPseudoDivision_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet);

void recursiveMDD_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, int nvar, int nSet);

//////////////////////////
// Subresultant Chains
//////////////////////////
/*
 *a linked list of AltArrZ_t 
 */
typedef struct AltArrsZ {
    AltArrZ_t* poly;
    struct AltArrsZ* next;
} AltArrsZ_t;

inline int Log2n (int n)
{
    return (n > 1) ? 1 + Log2n (n/2) : 0;
}

/** 
 * free an AltArrsZ_t 
 */
void freeAltArrsZ (AltArrsZ_t* AAs);
    
/**
 * Lazard Optimization (L. Docus, Optimizations of the Subresultant Algorithm, p4)  
 * Input Sd, Sdm
 * Output Se
 */ 
    AltArrZ_t* LazardOptZ (AltArrZ_t* Sd, AltArrZ_t* Sdm, AltArrZ_t* s);

/**
 * Ducos Optimization (L. Docus, Optimizations of the Subresultant Algorithm, p7)  
 * Input A, Sdm, Se, sd
 * Output Sem
*/
AltArrZ_t* DucosOptZ (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd);

/**
 * Cache friendly Ducos Optimization
 * Input A, Sdm, Se, sd
 * Output Sem
*/
AltArrZ_t* CFDucosOptZ (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd);

/**
 * Ducos SubresultantChain (L. Docus, Optimizations of the Subresultant Algorithm, p8)  
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * and store them on linked-list, AltArrsZ_t.
 * SC is the resultantchain and len is the size.
*/
void DucosSubresultantChainZ_rev (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len, int type);

void DucosSubresultantChainZ (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len);

    
/**
 * Resultant of two polynomials by DucosSubresultantChain
 * given two polynomials P and Q, this algorithm computes the resultant
 * of them by calling DucosSubresultantChain algorithm.
*/
AltArrZ_t* DucosResultantZ (AltArrZ_t* P, AltArrZ_t* Q);

/**
 * GCD of two polynomials by DucosSubresultantChain
 * given two polynomials P and Q, this algorithm computes the gcd
 * of them by calling DucosSubresultantChain algorithm.
*/
AltArrZ_t* DucosGCDZ (AltArrZ_t* P, AltArrZ_t* Q);

/**
 * GCD of two constant polynomials
 * given two constant polynomials P and Q, this algorithm computes the gcd
 * of them and return a constant polynomial.
*/
AltArrZ_t* integralGCD_AAZ_polyOut (AltArrZ_t* P, AltArrZ_t* Q);

/**
 * GCD of a polynomial and an Integer
 * given  polynomial P and integer c, this algorithm computes the gcd
 * of them and return a constant polynomial.
*/
AltArrZ_t* gcd_AAZ_Z (AltArrZ_t* P, mpz_t c);

/**
 * GCD of two polynomials
 * given two polynomials P and Q with the same main variable, 
 * this algorithm computes the gcd of them w.r.t the main 
 * variable and return a polynomial.
*/
AltArrZ_t* gcd_AAZ (AltArrZ_t* P, AltArrZ_t* Q);

/**
 * Primitive Factorization w.r.t the main variable
 * given a polynomial, P, this algorithm computes the primitive part
 * and the content of P w.r.t the main variable.
*/
AltArrZ_t* mainPrimitiveFactorization_AAZ (AltArrZ_t* P, AltArrZ_t** cont);

/**
 * Primitive Part w.r.t the main variable
 * given a polynomial, P, this algorithm computes the primitive part
 * of P w.r.t the main variable.
*/
AltArrZ_t* mainPrimitivePart_AAZ (AltArrZ_t* P, int mvar);


/**
 * Content w.r.t the main variable
 * given a polynomial, P, this algorithm computes the content
 * of P w.r.t the main variable.
*/
AltArrZ_t* mainContent_AAZ (AltArrZ_t* P);
    
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


