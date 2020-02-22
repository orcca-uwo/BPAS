#ifndef _SMQP_SUPPORT_RECURSIVE_H_
#define _SMQP_SUPPORT_RECURSIVE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "SMQP_Support.h"

/**
 * A recursively viewed term of a polynomial. 
 * In this representation, the coefficient is itself a polynomial and a single
 * exponent is encoded. 
 * Nodes can be strung together in a linked list to produce a polynomial.
 * A NULL RecNode is considered to be 0.
 */
typedef struct RecNode {
	Node* coef;
	degree_t exp;
	struct RecNode* next;
}RecNode_t;

/**
 * Free a recursively viewed polynomial node, and its constituents. 
 * However, it does not free the "next" RecNode_t in its chain.
*/
static inline void freeRecNode(RecNode_t* rNode) {
	if (rNode->coef != NULL) {
		freeNode(rNode->coef);
	}
	free(rNode);
}

/**
 * Creates a new RecNode_t and appends to the linked list whose tail is given 
 * as tail. 
 * 
 * returns the new tail of the list of RecNode_t nodes.
 */
static inline RecNode_t* addRecTerm(RecNode_t* tail, Node* coef, int exp) {
	RecNode_t* node = (RecNode_t*) malloc(sizeof(RecNode_t));
	node->coef = coef;
	node->exp = exp;
	node->next = NULL;
	if (tail != NULL) {
		tail->next = node;
	}
	return node;
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

/***********
 * Conversion to and from recursive view.
 ***********/

/**
 * Convert a polynomial to a recursively viewed polynomial given the index
 * within the exponent vector of the variable to become the "main variable"
 * of the recursive view.
 *
 * Note this conversion is done INPLACE. And therefore the input Node* is 
 * invalidated by this operation.
 *
 * Note: considered that the previous values before idx in exponent vectors
 * are zero.
 *
 * returns the recursively viewed polynomial. 
 */
RecNode_t* convertToRecursiveNodeAtIdx(Node* poly, int idx);

/**
 * Convert a polynomial to a recursively viewed polynomial using the variable
 * of highest precedence (index 0) as the main variable.
 *
 * Note this conversion is done INPLACE. And therefore the input Node* is 
 * invalidated by this operation.
 *
 * returns teh recursively viewed polynomial.
 */
RecNode_t* convertToRecursiveNode(Node* poly);

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
 * Convert representation of polynomial poly from recursive to sparse,
 * such that the main degree of input polynomial is idx.
 *
 * If idx = 0, then this algorithm is equal to convertFromRecursiveNode.
 *
 * Note: considered that the previous values before idx in exponent
 * vectors are zero.
 */
Node* convertFromRecursiveNodeAtIdx (RecNode_t* poly, int idx);

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done INPLACE and therefore the original recursively
 * viewed polynomial is then invalidated.
 * 
 * returns the distributed polynomial.
 */
Node* convertFromRecursiveNode(RecNode_t* poly);

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done out of place and therefore the original recursively
 * viewed polynomial is not affected.
 * 
 * returns the distributed polynomial.
 */
Node* convertFromRecursiveNodeCopy(RecNode_t* poly, int nvar);

/**
 * Given a recursive polynomial whose head is poly, make a deep copy of it.
 *
 * returns the head to the new polynomial.
 */
RecNode_t* deepCopyRecPolynomial(RecNode_t* poly, int nvar);

/***********
 * Pesudo Division
 ***********/

/**
 * An element of the product heap for recursively viewed polynomials.
 * Holds pointers to a_i, the current element of a being distributed over b, 
 * and to b, the multiplicand polynomial. 
 */
typedef struct RecProdHeapElem {
	RecNode_t* a_i;
	RecNode_t* b;
	degree_t prodExp;
	struct RecProdHeapElem* next; //chained heaps;
}RecProdHeapElem_t;

/** 
 * The heap used for multiplication between recursively viewed polynomials.
 * This is a chained heap, each element is a linked list of terms with same degree.
 */
typedef struct {
	RecProdHeapElem_t** elements;
	polysize_t heapSize;
	polysize_t maxHeapSize;
	int nvar;
}RecProdHeap_t;

/**
 * Create a new element to be inserted into the heap. 
 * Takes one term of the multiplied and the entire multiplicand. 
 *
 * Returns an element ready to be inserted in the RecProdHeap. 
 */
RecProdHeapElem_t* recProdHeapMakeElement(RecNode_t* a_i, RecNode_t* b, int nvar);

/** 
 * Create a RecProdHeap.
 */
RecProdHeap_t* recProdHeapCreate(int nvar);

/**
 * Insert an RecProdHeap element into the heap h.  
 */
void recProdHeapInsert(RecProdHeap_t* h, RecProdHeapElem_t* elem);

/**
 * Peek at the maximum degree product term in the heap.
 *
 * returns the degree of the product whose degree is maximum in the heap.
 */
degree_t recProdHeapPeek(RecProdHeap_t* h);

/**
 * Get the coefficient of the next term to be produced in the multiplication.
 * Terms are produced in decreasing order by degree.
 * Note that this may return NULL to indicate that term whose degree was 
 * identified as maximum by recProdHeapPeek, ended up being cancelled to 0
 * when combining like terms.
 *
 * returns the coefficient of the next product term, or NULL, as above.
 */
Node* recProdHeapGetNextCoef(RecProdHeap_t* h);

/** 
 * Do the pesudo division of c by b. Quotient is returned in res_a, and remainder
 * in res_r. 
 *
 * If e is not NULL then it returns the exact number of division steps which occurred.
 *
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 */
void pesudoDivide(RecNode_t* c, RecNode_t* b, RecNode_t** res_a, RecNode_t** res_r, int* e, Node** hPow, int nvar, int lazy);

/*
 * Addition
 */

/**
 * Add two recursive polynomials given their head nodes, a and b.
 * nvar: number of variables in the polynomial.
 * returns a pointer to the head Node of the sum.
 */
RecNode_t* addRecPolynomials(RecNode_t *a, RecNode_t *b, int nvar);

/**
 * Compute the multi divisor pseudo-division polynomial f by triangular set T
 */
void multiDivisorPseudoDivision (RecNode_t* f, RecNode_t** T, RecNode_t** quoSet, RecNode_t** rem, Node** H, int nvar, int nSet, int lazy);

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


