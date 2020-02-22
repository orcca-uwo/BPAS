
#ifndef _SMQP_SUPPORT_RECURSIVE_AA_H_
#define _SMQP_SUPPORT_RECURSIVE_AA_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "SMQP_Support-AA.h"

#define PDIVIDE_DIVISBLE_CHECK 1

/**
 * A recursively viewed term of a polynomial. 
 * In this representation, the coefficient is itself a polynomial and a single 
 * exponent is encoded. 
 * Nodes can be strung together in a linked list to produce a polynomial.
 * A NULL RecNode is considered to be 0.
 */
typedef struct RecNode {
	int coefSize;
	degree_t exp;
	AAElem_t* coef;
	// int powersOfH;
	struct RecNode* next;
} RecNode_t;

/**
 * Free a recursively viewed polynomial node 
 * However, it does not free the "next" RecNode_t in its chain nor the coef data. 
 */
static inline void freeRecNode(RecNode_t* rNode) {
	// if (rNode->coef != NULL) {
	// 	freeNode(rNode->coef);
	// }
	free(rNode);
}

/**
 * Creates a new RecNode_t and appends to the linked list whose tail is given 
 * as tail. 
 * 
 * returns the new tail of the list of RecNode_t nodes.
 */
static inline RecNode_t* addRecTerm(RecNode_t* tail, AAElem_t* coef, int coefSize, int exp) {
	RecNode_t* node = (RecNode_t*) malloc(sizeof(RecNode_t));
	node->coef = coef;
	node->coefSize = coefSize;
	node->exp = exp;
	// node->powersOfH = 0;
	node->next = NULL;
	if (tail != NULL) {
		tail->next = node;
	}
	return node;
}

typedef struct {
	int coefSize;
	degree_t exp;
	AAElem_t* coef;
} RecArrElem_t;

typedef struct {
	int alloc;
	int size;
	AltArr_t* origAA;
	RecArrElem_t* elems;
} RecArr_t;

/**
 * Free a RecArr_t and it's RecArrElem_t array. but not the underlying coef data.
 */
static inline void freeRecArray(RecArr_t* rArr) {
	if (rArr != NULL) {
		if (rArr->alloc > 0) {
			free(rArr->elems);
		}
		free(rArr);
	}
}

/** 
 * Build a random recursive polynomial.
 * Nterms is the number of terms wrt the main variable. 
 * coefNTerms is the number of terms in the coefficients (when viewed recursively, and nvar > 1);
 * coefBound is the absolute upper limit on the rational numbers in the polynomial.
 * sparsity is the upper limit on the difference between degrees of successive recursive terms. 
 * sparsity is also passed when constructing coefficient polynomials.
 */
RecNode_t* buildRandomRecPoly(int nvar, int nterms, int coefNterms, unsigned long int coefBound, degree_t sparsity, int includeNeg);
RecArr_t* buildRandomRecArrPoly(int nvar, int nterms, int coefNterms, unsigned long int coefBound, degree_t sparsity, int includeNeg);


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
// RecNode_t* convertToRecursiveNode(Node* poly);
RecNode_t* convertToRecursiveNode(AltArr_t* aa);
RecArr_t* convertToRecursiveArray(AltArr_t* aa);

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
// RecNode_t* convertToRecursiveNodeAtIdx(Node* poly, int idx);
RecNode_t* convertToRecursiveNodeAtIdx(AltArr_t* aa, int idx);
RecArr_t* convertToRecursiveArrayAtIdx (AltArr_t* aa, int idx);

/**
 * Convert a polynomial to a recurisvely viewed polynomial using the varibale
 * of highest precdence (index 0) as the main variable.
 *
 * Note that this conversion is out-of-place. And therefore input Node* is 
 * unaffected. 
 *
 * returns the recursively viewed polynomial.
 */
RecNode_t* convertToRecursiveNodeCopy(Node* poly, int nvar);

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done INPLACE and therefore the original recursively
 * viewed polynomial is then invalidated.
 * 
 * nvar is the number of variables the returning AltArr_t should have. That is, 
 * the number of variables in the RecArr's coefficients + 1.
 *
 * returns the distributed polynomial.
 */
// Node* convertFromRecursiveNode(RecNode_t* poly);
// AltArr_t* convertFromRecursiveNode(RecNode_t* poly);
AltArr_t* convertFromRecursiveArray(RecArr_t* poly, int nvar);

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done out of place and therefore the original recursively
 * viewed polynomial is not affected.
 * 
 * returns the distributed polynomial.
 */
// Node* convertFromRecursiveNodeCopy(RecNode_t* poly, int nvar);

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
AltArr_t* convertFromRecursiveArrayAtIdx (RecArr_t* recPoly, int idx, int nvar);


/**
 * Given a recursive polynomial whose head is poly, make a deep copy of it.
 *
 * returns the head to the new polynomial.
 */
// RecNode_t* deepCopyRecPolynomial(RecNode_t* poly, int nvar);
RecArr_t* deepCopyRecArrayPolynomial(RecArr_t* poly, int nvar);



/***********
 * Pesudo Division
 ***********/

typedef struct RecProdHeapChain_AA {
	int a_i;
	int b;
	struct RecProdHeapChain_AA* next;
} RecProdHeapChain_AA;

typedef struct RecProdHeapElem_AA {
	degree_t exp;
	RecProdHeapChain_AA* chain;
} RecProdHeapElem_AA;

typedef struct {
	RecProdHeapElem_AA* elements;
	int heapSize;
	int maxHeapSize;
	int nvar;
	int lastB;
} RecProdHeap_AA;

/**
 * Create a new chain to be inserted into the heap. 
 * Takes integers representing the index of polynomial terms in 
 * the multiplier and the multiplicand. 
 *
 * returns a chain ready to be inserted in the heap. 
 */
static inline RecProdHeapChain_AA* recProdHeapMakeChain_AA(int a , int b, RecProdHeapChain_AA* next) {
	RecProdHeapChain_AA* chain = (RecProdHeapChain_AA*) malloc(sizeof(RecProdHeapChain_AA));
	chain->a_i = a;
	chain->b = b;
	chain->next = next;
	return chain;
}

/**
 * Free a chain, and anything chained to it. 
 */
static inline void recProdHeapFreeChain_AA(RecProdHeapChain_AA* chain) {
	RecProdHeapChain_AA* next;
	while(chain != NULL) {
		next = chain->next;
		free(chain);
		chain = next;
	}
}

/**
 * Create an empty product heap. 
 */
static inline RecProdHeap_AA* recProdHeapCreate_AA(int nvar) {
	RecProdHeap_AA* h = (RecProdHeap_AA*) malloc(sizeof(RecProdHeap_AA));
	h->elements = NULL;
	h->heapSize = 0;
	h->maxHeapSize = 0;
	h->nvar = nvar;
	return h;
}

/**
 * Allocate space for a maximum of newAllocSize elements in the heap
 */ 
static inline void recProdheapResize_AA(RecProdHeap_AA* h, int newAllocSize) {
	h->elements = (RecProdHeapElem_AA*) realloc(h->elements, sizeof(RecProdHeapElem_AA)*newAllocSize);
	h->maxHeapSize  = newAllocSize;
}

static inline void recProdHeapFree_AA(RecProdHeap_AA* h) {
	RecProdHeapElem_AA* elems = h->elements;
	polysize_t s = h->heapSize;
	for (polysize_t i = 0; i < s; ++i) {
		recProdHeapFreeChain_AA(elems[i].chain);
	}
	free(h->elements);
	free(h); 
}

/**
 * Insert aa RecProdHeapChain into the heap h, given that the chain represents 
 * a degree equal to insertExp.  
 */
void recProdHeapInsert_AA(RecProdHeap_AA* h, RecProdHeapChain_AA* chain, register degree_t insertExp);

/**
 * Peak into the heap and get the exponent vector of the product of the max elem.
 */
static inline long long int recProdHeapPeek_AA(RecProdHeap_AA* h) {
	if (h->heapSize > 0) {
		return h->elements->exp;
	}
	return -1;
}

/**
 * Remove the maximum chain from the product heap.
 */
RecProdHeapChain_AA* recProdHeapRemoveMax_AA(RecProdHeap_AA* h);

/**
 * Get the coefficient of the next term to be produced in the multiplication.
 * Terms are produced in decreasing order by degree.
 * Note that this may return NULL to indicate that term whose degree was 
 * identified as maximum by recProdHeapPeek, ended up being cancelled to 0
 * when combining like terms.
 *
 * Since the chains only store integers, we must pass in the array of polynomial
 * terms for each of the multiplier and multiplicand as aElems and bElems.
 *
 * returns the coefficient of the next product term, or NULL, as above.
 */
AltArr_t* recProdHeapGetNextCoef_AA(RecProdHeap_AA* h, const RecArrElem_t* __restrict__ aElems, const RecArrElem_t* __restrict__ bElems);

/** 
 * Do the pesudo division of c by b. Quotient is returned in res_a, and remainder
 * in res_r. 
 *
 * If e is not NULL then it returns the exact number of division steps which occurred.
 *
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 */
void pesudoDivide_RecArray(RecArr_t* c, RecArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy);


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
void pesudoDivideAtIdx_AA(int idx, AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy);

void naiveMultiDivisorPseudoDivision_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet);

void multiDivisorPseudoDivision_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet);

void recursiveMDD_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, int nvar, int nSet);

//////////////////////////
// Subresultant Chains
//////////////////////////
/*
 *a linked list of AltArr_t 
 */
typedef struct AltArrs {
    AltArr_t* poly;
    struct AltArrs* next;
} AltArrs_t;

/** 
 * free a AltArrs_t 
*/
void freeAltArrs (AltArrs_t* AAs);
    
/**
 * Lazard Optimization (L. Docus, Optimizations of the Subresultant Algorithm, p4)  
 * Input Sd, Sdm
 * Output Se
*/ 
AltArr_t* LazardOpt(AltArr_t* Sd, AltArr_t* Sdm);

/**
 * Ducos Optimization (L. Docus, Optimizations of the Subresultant Algorithm, p7)  
 * Input A, Sdm, Se, sd
 * Output Sem
*/
AltArr_t* DucosOpt(AltArr_t* A, AltArr_t* Sdm, AltArr_t* Se, AltArr_t* sd);

/**
 * Ducos SubresultantChain (L. Docus, Optimizations of the Subresultant Algorithm, p8)  
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * and store them on linked-list, AltArrs_t.
 * SC is the resultantchain and len is the size.
*/
void DucosSubresultantChain_inv (AltArr_t* P, AltArr_t* Q, AltArrs_t** SC, int* len);
void DucosSubresultantChain (AltArr_t* P, AltArr_t* Q, AltArrs_t** SC, int* len);

    
/**
 * Resultant of two polynomials by DucosSubresultantChain
 * given two polynomials P and Q, this algorithm computes the resultant
 * of them by calling DucosSubresultantChain algorithm.
*/
AltArr_t* DucosResultant (AltArr_t* P, AltArr_t* Q);

/**
 * GCD of two polynomials by DucosSubresultantChain
 * given two polynomials P and Q, this algorithm computes the gcd
 * of them by calling DucosSubresultantChain algorithm.
*/
AltArr_t* DucosGCD (AltArr_t* P, AltArr_t* Q);
    
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


