
/*****
 * Supporting methods for sparse multivariate rational polynomials written in 
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

#ifndef _SMQP_SUPPORT_H_
#define _SMQP_SUPPORT_H_

#ifdef __cplusplus
extern "C" { 
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>

/**
 * Defines whether to use chained heaps over basic heaps.
 */
#define SMQP_SUPPORT_CHAINED_HEAP 1

/**
 * Compile methods for debugging SMQP Support.
 */
#define SMQP_SUPPORT_DEBUG 0

/** 
 * Define if SMQP should count number of comparisons in a global variable
 */
#define SMQP_COUNT_COMPARISONS 0

/**
 * Define if SMQP should count heap chains and inserts.
 */
#define SMQP_COUNT_CHAINS 0



/*****************
 * Data types used throughout 
 *****************/

typedef unsigned long int polysize_t;

typedef polysize_t degree_t;

typedef degree_t* degrees_t; //degrees_t is an exponent vector

typedef mpq_t ratNum_t;

/*****************
 * Exponent Vector functions.
 *****************/

#if SMQP_COUNT_COMPARISONS
	static unsigned long long smqpCompareCount = 0;
#endif

/**
 * Compare two exponent vectors for lexicographical term order.
 * a: int array representing exponent vector
 * b: int array representing exponent vector
 * v: size of the arrays
 * returns -1 if a < b, 0 if equal, 1 if a > b
 */
int compareExponentVectors(degrees_t a, degrees_t b, int v);

/**
 * Check if an exponent vector has all 0 exponents.
 * a: the exponent vector to check
 * v: size of a
 */ 
int isZeroExponentVector(degrees_t a, int v);

/**
 * Add two exponent vectors, a and b, and store the result in c.
 * a,b,c: degree_tarrays 
 * v: the number of elements in the arrays
 */
static inline void addExponentVectors(degrees_t a, degrees_t b, degrees_t c, int v){
	for (int i = 0; i < v; ++i) {
		c[i] = a[i] + b[i];				
	}
}

/**
 * Subtract two exponent vectors, b from a, and store the result in c.
 * It is assumed that this operation is valid.
 * see also, monomialDivideTest.
 *
 * a,b,c: degree_tarrays
 * v: the number of elements in the arrays
 */
static inline void subtractExponentVectors(degrees_t a, degrees_t b, degrees_t c, int v) {
	for (int i = 0; i < v; ++i) {
		c[i] = a[i] - b[i];
	}

}



/*****************
 * Node definition and helper functions.
 *****************/

/**
 * A Node struct containing a single term of a polynomial, that is, 
 * the coefficient and the monomial's exponent vector. It also contains
 * a pointer to the next Node in the polynomial. In this way, a Node represents
 * both a single term and a polynomial. 
 *
 * Note that a NULL Node is considered to be 0. 
 */

typedef struct Node {
	ratNum_t coef;
	degrees_t degs;
	// a pointer connecting to the next term
	struct Node* next;
} Node;

/**
 * Free a node and its exponent vector
 */
static inline void freeNode(Node* node) {
	if (node != NULL) {
		mpq_clear(node->coef);
		free(node->degs);
		free(node);
	}
}

/**
 * Create a copy of supplied Node, except for the next pointer.
 * returns this new copy.
 */
static inline Node* deepCopyNode(Node* node, int nvar) {
	Node* newNode = (Node*) calloc(1, sizeof(Node));
	newNode->degs = (degrees_t) malloc(sizeof(degree_t) * nvar);
	memcpy(newNode->degs, node->degs, sizeof(degree_t) * nvar);
	mpq_init(newNode->coef);
	mpq_set(newNode->coef, node->coef);
	newNode->next = NULL;

	return newNode;
}

/**
 * Free a node, and all it's successors. 
 */
void freePolynomial(Node* node);

/**
 * Deep copy the polynomial whose head in the input Node*.
 */
Node* deepCopyPolynomial(Node* node, int nvar);

/**
 * Determine the number of terms in a polynomial given the head node 
 * of the linked-list representation of the polynomial.
 * returns the number of terms including and following the supplied node 
 */
polysize_t numberOfTermsNode(Node* node);

/**
 * Given a linked-list of nodes, whose head is poly, 
 * sort the list based on the ordering of compareExponentVectors.
 * Updates the polynomial in place by modifying *poly to be the new head node.
 * Also returns a pointer which is the new head node. 
 */
Node* sortPolynomial(Node** poly, int nvar);

/**
 * Given a polynomial in sorted order but with possible duplicates, 
 * condense the list such that like-terms are combined into a single node.
 */
void condensePolyomial(Node* poly, int nvar);

/**
 * Add a term to the end of the polynomial linked list given the tail node,
 * trailingTerm, the exponent vector for the new term, d, and the coefficient, 
 * ratNum_t. A copy of the input coef is made and stored in the node. But the
 * degrees_t is stored directly. 
 *
 * returns a pointer to the new tail of the list (the node created)
 */
Node* addTerm(Node* trailingTerm, degrees_t d, const ratNum_t coef);

/**
 * Add a term to the end of linked list with tail trailingTerm. The new term
 * will have degs d and a coef of 0.
 * 
 * returns a pointer to the new tail of the list (the node created)
 */
Node* addZeroTerm(Node* trailingTerm, degrees_t d);

/** 
 * Given two terms, as nodes a and b, multiply these two and return the
 * single result as a node.
 */
Node* multiplyTerms(Node* a, Node* b, int nvar);

/**
 * Negate a polynomial. 
 * Given the head of a polynomial, negate all coefficients in place.
 */
void negatePolynomial(Node* a);

/**
 * Evaluate a polynomial whose head is given by a. 
 * This method returns another polynomial as not all variables need
 * to have values supplied. But, if they do, a constant term will be returned. 
 * both active and vals are arrays of size nvar. active[i] determines if 
 * the variable at degs[i] is to be evaluated using vals[i] as value.
 */
Node* evaluatePoly(Node* a, int* active, mpq_t* vals, int nvar);


/**
 * Unchecked!
 *
 * Given polynomial a, return the leading term
 */ 
Node* leadingTerm (Node* a);



/*****************
 * SMQP Addition and Subtraction
 *****************/

/**
 * Add two polynomials given their head nodes, a and b.
 * nvar: number of variables in the polynomials.
 * returns a pointer to the head Node of the sum.
 */
Node* addPolynomials(Node *a, Node *b, int nvar);


/**
 * Add polynomial a and b, putting the sum back into a.
 * 
 * nvar: size of exponent vectors in a and b.
 *
 * returns a pointer to the new head of a.
 */
Node* addPolynomials_inp(Node* a, Node* b, int nvar);


/**
 * Subtract b from a, putting the difference back into a.
 * 
 * nvar: size of exponent vectors in a and b.
 *
 * returns a pointer to the new head of a.
 */
Node* subPolynomials_inp(Node* a, Node* b, int nvar);




/*****************
 * SMQP Multiplication & Helpers
 *****************/

/**
 * Data element for the ProductHeap data structure
 */
typedef struct productHeapElem {
	Node* a_i;
	Node* b;
	Node* product;
#if SMQP_SUPPORT_CHAINED_HEAP
	struct productHeapElem* next;
#endif
} productHeapElem;

/**
 * Data structure 'heap' holding the intermediate terms produced during multiplication.
 */
typedef struct {
	productHeapElem** elements;
	polysize_t heapSize;
	polysize_t maxHeapSize;
	int nvar;
} ProductHeap;

#if SMQP_SUPPORT_DEBUG
/**
 * Print the product degrees_t currently in the heap.
 */
void prodheapPrint(ProductHeap* h, int nvar);
#endif

/**
 * Make an element for the product heap, combining nodes a and b as the 
 * element's product.
 */
productHeapElem* prodheapMakeElement(ProductHeap* h, Node* a, Node* b);

/**
 * Free a product heap element and its product node. 
 */
static inline void prodheapFreeElement(productHeapElem* elem) {
	freeNode(elem->product);
	free(elem);
}

/**
 * Cleanup memory initialized for the heap
 */
void prodheapFree(ProductHeap* h);

/**
 * Create an empty product heap. 
 */
ProductHeap* prodheapCreate(int nvar);

/**
 * Initialize the product heap with the two polynomials to multiply, a and b.
 *
 * We know the maximum heap size is numTerms(a) as at most one of a_i*b
 * is in the heap at once. 
 */
ProductHeap* prodheapInit(Node* a, Node* b, int nvar);

/**
 * Given an index, i, into the heap, swim that node up so that the heap 
 * is properly ordered. 
 */
void prodheapSwim(ProductHeap* h, int index);

/**
 * Given an index, i, into the heap, sink that node down so that the heap 
 * is properly ordered. 
 * Note: This is not used in the chained heap case
 */
void prodheapSink(ProductHeap* h, int i);

#if SMQP_COUNT_CHAINS
	static unsigned long long int SMQP_CHAINS = 0;
	static unsigned long long int SMQP_INSERTS = 0;
#endif 

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
void prodheapInsert(ProductHeap* h, productHeapElem* elem);

/**
 * Peek into the heap, h, to get the exponent vector of the product
 * of the max element.
 */
degrees_t prodheapPeek(ProductHeap* h);

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
productHeapElem* prodheapRemoveMax(ProductHeap* h);

/**
 * Extract the maximum product from the heap, returned as the Node*. 
 * This automatically adds the next term in the stream of the chosen maximum
 * if such a term exists.
 */
Node* prodheapExtract(ProductHeap* h);

/**
 * Find the next product term of highest degree_t
 * This searches through each i in I..n streams representing the expanded multiplcation
 * of node a_i with b. Picking the maximal a_i*b_fi.
 * 
 * returns the terms to multiply as a_s and b_fs and the value s.
 *
 * This search is done linearly, for use in the multiplication algorithm which does
 * not use heaps.
 */
void findNextProductTerm(Node* a_I, Node** b_fi, polysize_t numTermsA, polysize_t I, int* s, Node** a_s, Node** b_fs, int nvar);

/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure. 
 * It is assumed that both a and b have compatible exponent vectors.
 * 
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
Node* multiplyPolynomials(Node* a, Node* b, int nvar);


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
Node* multiplyPolynomials_inp(Node* a, Node* b, int nvar);


/**
 * Struct to store the Multiplication Producer's state
 */
typedef struct {
	Node* a_I;
	Node** f;
	polysize_t nTerms; 
	int nvar;
	polysize_t I;
	ratNum_t lastCoef;
	degrees_t lastDeg;
} MultState;

/** 
 * Struct to return from producer. The updated state and the result
 */
typedef struct {
	MultState* state;
	Node* nextRes;
} MultResult; 

/** 
 * Initialize the multiplication producer's state given the 
 * two polynomials to multiply and the number of variables of
 * the polynomials.
 */
MultState* initProducerState(Node* a, Node* b, int nvar);

/** 
 * Cleanup the producer state memory allocation.
 */
void releaseProducerState(MultState* state);

/**
 * A 'producer' for the multiplication of two polynomials. 
 * This producer, using the input state, will produce a single
 * term of the product and then return. 
 *
 * Drivers must repeatedly call this to obtain terms, one after the other.
 * Terms are produced in descending order with like-terms already combined.
 * returns a NULL result if no more terms can be produced.
 *
 * Note the underlying multiplication algorithm is not efficient but this can 
 * be handy when only a few leading terms of the product are needed.
 */
MultResult multPolyProducer(MultState* state);



/*****************
 * Polynomial exponentiation 
 *****************/

/**
 * Given a polynomial, a, compute a^n.
 * 
 * n: a positive integer
 * nvar: the number of variables of a
 * 
 */
Node* exponentiatePoly(Node* a, unsigned int n, int nvar);




/*****************
 * Polynomial division
 *****************/

/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients 
 * nvar: number of variables of monomials a and b
 */
int monomialDivideTest(degrees_t adegs, degrees_t bdegs, int nvar);

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree_t
 */
Node* divisionGetNextTerm(ProductHeap* h);

/** 
 * Given a polynomial, c, and a term, b, determine polynomials a and r
 * such that c = b*a + r.
 * a and r are returned in res_a and res_r, respectively. 
 */
void divideBySingleTerm(Node* c, Node* b, Node** res_a, Node** res_r, int nvar);

/** 
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r 
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials(Node* c, Node* b, Node** res_a, Node** res_r, int nvar);


/****************
* Primitive Factorization
*****************/
typedef struct mpzCoefs_t 
{
  mpz_t mpzCoef; 
  struct mpzCoefs_t* next;
}mpzCoefs_t;

typedef struct content_t 
{
  mpz_t cont; 
}content_t;

/**
 * @param[in] nvar The number of variable
 * @param[out] seq The recursive loop
 * \brief Return content of a where a = content(a) * primitivePart(a)
 */
int* recursiveLoop (int nvar);

/** 
 * @param[in] a The polynomial
 * @param[out] cont The content of a
 * \brief Return content of a where a = content(a) * primitivePart(a)
 */
content_t* content (Node* a);

/**
 * @param[in] a The polynomial
 * @param[in] nvarThe number of variables
 * @param[out] returnContent The content of a
 * @param[out] head The primitive part of a
 * \brief Compute the primitive part of a and relevant content in returnContent.
 * For example, let  a = 4*x^2y^3 + 8*x*y^2 + 16 the primitivePart(a, returnContent, 2) 
 * returns x^2y^3 + 2*x*y^2 + 4 and the returnContent becomes 4.  
 */
Node* primitivePart (Node* a, content_t** returnContent,  int nvar);


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
void heapMDD (Node* f, Node** G, Node** Q, Node** r, int s, int nvar);

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
void triangularSetMDD (Node* f, Node** G, Node** Q, Node** r, int s, int nvar);

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
void primitiveFactorTriangularSetMDD (Node* f, Node** G, Node** Q, Node** r, int s, int nvar);

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
void multiDivisorDivision (Node* f, Node** G, Node** Q, Node** r, int s, int nvar, int type);

/***************
* MDD Test and Verification
***************/

/**
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[out] 0-1 Return 1, unless there is a NULL polynomial
 */
int multiDivisorDivisionIsEmpty (Node* f, Node** G, int s);


/** Multi-Divisor Division Verification of Definition:
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[in] rr The remainder
 * @param[in] Q The quotient-set
 * @param[out] 0-1 Return 1, unless Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r != f
 * \brief This algorithm tests the Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r -f = 0.
 */
int multiDivisorDivisionVerif (Node* ff, Node** G, Node** Q, Node* rr, int s, int nvar);

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


