
#ifndef _SMQP_SUPPORT_RECURSIVE_AA_H_
#define _SMQP_SUPPORT_RECURSIVE_AA_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "SMQP_Support-AA.h"
#include "SMQP_Support_Recursive_Unpacked.h"

//NOTE: SMZP currently broken if this is off.
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

typedef struct RecArrElem {
	int coefSize;
	degree_t exp;
	AAElem_t* coef;
} RecArrElem_t;

typedef struct RecArr {
	int alloc;
	int size;
	int unpacked;
	AltArr_t* origAA;
	RecArrElem_t* elems;
} RecArr_t;

/**
 * Free a RecArr_t and it's RecArrElem_t array. but NOT the underlying coef data.
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
 * Free a RecArr_t and it's RecArrElem_t array. as well as the underlying coef data.
 */
static inline void freeRecArrayAndCoef(RecArr_t* rArr) {
	if (rArr != NULL) {
		if (rArr->alloc > 0) {
			if (rArr->unpacked) {
				free( (degree_t*) rArr->elems->coef->degs);
			}
			free(rArr->elems->coef);
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

typedef struct RecProdHeap_AA {
	RecProdHeapElem_AA* elements;
	int heapSize;
	int maxHeapSize;
	int nvar;
	int lastB;
	int unpacked;
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
	h->unpacked = 0;
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
// AltArr_t* recProdHeapGetNextCoef_AA(RecProdHeap_AA* h, const RecArrElem_t* __restrict__ aElems, const RecArrElem_t* __restrict__ bElems);
AltArr_t* recProdHeapGetNextCoef_AA(RecProdHeap_AA* h, const RecArrElem_t* __restrict__ aElems, const RecArrElem_t* __restrict__ bElems, const int* aUnpacked, int bUnpacked);

void pesudoDivideOneTerm_RecArray(RecArr_t* c, RecArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy);

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
// Multi-divisor (Pseudo-)division
//////////////////////////////////

/** 
 * Do the pseudo division of c by b when the idx-th member
 * of exponent vectors is the main variable. 
 * Quotient is returned in res_a, and remainder in res_r. 
 * If e is not NULL then it returns the exact number of division steps which occurred.
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 */
void pesudoDivideAtIdx_AA(int idx, AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy);

/** 
 * Do the pseudo division of c by the triangular-set of B
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 *
 * Note this algorithm is based on both naive and recursive principles.
 */
void multiDivisorPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet);

/** 
 * Do the pseudo division of c by the triangular-set of B in the naive principle
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 */
void naiveMultiDivisorPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet);


/** 
 * Do the pseudo division of c by the normalized triangular-set of B in the recursive principle 
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 */
void normalizedTriangularSetPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t*** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet);

/** 
 * Multi-Divisor Division (MDD) where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
void recursiveTriangularSetMDD_AA (AltArr_t* c, AltArr_t** B, AltArr_t*** quoSet, AltArr_t** rem, int nvar, int nSet);

	
//////////////////////////
// Subresultant Chains
//////////////////////////

/**
 * a linked list of AltArr_t 
 */
typedef struct AltArrs {
    AltArr_t* poly;
    struct AltArrs* next;
} AltArrs_t;


/** 
 * free a AltArrs_t 
*/
static inline void freeAltArrs (AltArrs_t* AAs)
{
	while (AAs == NULL){
		return;
	}
	AltArrs_t* cur = AAs;
	while (cur != NULL){
		if (cur->poly != NULL){
			freePolynomial_AA(cur->poly);
		}
		cur = cur->next;
	}
}

static inline int Log2_AA (int n){
    return (n > 1) ? 1 + Log2_AA (n/2) : 0;
}



AltArr_t* LazardOpt (AltArr_t* Sd, AltArr_t* Sdm, AltArr_t* s);

/**
 * Ducos Optimization (L. Docus, Optimizations of the Subresultant Algorithm, p7)  
 * Input A, Sdm, Se, sd
 * Output Sem
*/
AltArr_t* DucosOpt (AltArr_t* A, AltArr_t* Sdm, AltArr_t* Se, AltArr_t* sd);

/**
 * Cache friendly Ducos Optimization
 * Input A, Sdm, Se, sd
 * Output Sem
*/
AltArr_t* CFDucosOpt (AltArr_t* A, AltArr_t* Sdm, AltArr_t* Se, AltArr_t* sd);

/**
 * Ducos SubresultantChain (L. Docus, Optimizations of the Subresultant Algorithm, p8)  
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * and store them on linked-list, AltArrs_t.
 * SC is the resultantchain and len is the size.
*/
void DucosSubresultantChain_rev (AltArr_t* P, AltArr_t* Q, AltArrs_t** SC, int* len, int type);
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

/** 
 * Return the last non-zero polynomial in subresultant chain of two polynomials P and Q.
*/	
AltArr_t* lastNonZeroChain_AA (AltArr_t* P, AltArr_t* Q);


///////////////////////////////
// Extended Subresultant Chain 
///////////////////////////////

/**
 * To store the Bezout Coefficients and partial gcds
 */
typedef struct exgcds {
	AltArr_t* r;
	AltArr_t* a;
	AltArr_t* b;
	struct exgcds* next;
} exgcds_t;

/**
 * Free exgcds_t
 */
static inline void freeExgcds_AA (exgcds_t* lg) {
	if (lg != NULL) {
		exgcds_t* cur = lg;
		while (cur != NULL) {
	 		if (cur->r != NULL) {
 				freePolynomial_AA (cur->r);
 			}
	 		if (cur->a != NULL) {
 				freePolynomial_AA (cur->a);
 			}
		 	if (cur->b != NULL) {
	 			freePolynomial_AA (cur->b);
	 		}
	 		cur = cur->next;
	 	}
	}
}

/**
 * The semi-lazard optimization based on "L. Docus, Optimizations of the Subresultant Algorithm".
 * used in exLazardOpt_AA algorithm
 */
AltArr_t* semiLazardOpt_AA  (AltArr_t* Sd, AltArr_t* Sdm, AltArr_t* s);
AltArr_t* exLazardOpt_AA (AltArr_t* Sd, exgcds_t* VSdm, AltArr_t* s, AltArr_t** hc, AltArr_t** qc);

/**
 * Extended Ducos Optimization  
 * Input A, VSdm(= [Sdm, VSdm->a, VSdm->b]), Se, sd
 * Output Sem, and return by reference hc and qc. 
*/
AltArr_t* exDucosOpt_AA (exgcds_t* VSd, exgcds_t* VSdm, AltArr_t* Se, AltArr_t* sd, AltArr_t** hb, AltArr_t** qb);

/**
 * Extended Ducos SubresultantChain  
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * with all partial Bezout Coefficients and store them on linked-list, exgcds_t.
 * SC is the resultantchain and len is the size.
*/
void exDucosSubresultantChain_rev_AA (AltArr_t* P, AltArr_t* Q, exgcds_t** SC, int* len, int type);
void exDucosSubresultantChain_AA (AltArr_t* P, AltArr_t* Q, exgcds_t** SC, int* len);

/**
 * Extended Resultant of two polynomials by DucosSubresultantChain
 * given two polynomials P and Q, this algorithm computes the resultant
 * with its Bezout Coefficients by calling exDucosSubresultantChain_AA algorithm.
*/
AltArr_t* exDucosResultant_AA (AltArr_t* P, AltArr_t* Q, AltArr_t** a, AltArr_t** b);


////////////////////////////
/// Mod Subresultant Chain 
////////////////////////////
/**
 * Ducos SubresultantChain Mod alpha
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * and store them on linked-list, AltArrs_t.
 * SC is the resultantchain and len is the size.
 * Note the main variable of alpha should be less than the main variable of p and q.
 * Ex: For bivariate p(x,y) and q(x,y), alpha must satisfy deg(alpha, x) = 0.  
*/
void modDucosSubresultantChain_rev_AA (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, AltArrs_t** SC, int* len, int type, int isDeepMod);
void modDucosSubresultantChain_AA (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, AltArrs_t** SC, int* len, int isDeepMod);


/**
 * Resultant of two polynomials Mod alpha by modDucosSubresultantChain 
 * given two polynomials P and Q, this algorithm computes the resultant
 * of them by calling DucosSubresultantChain algorithm.
*/
AltArr_t* modDucosResultant_AA  (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, int isDeepMod);

/** 
 * Return the last non-zero polynomial in subresultant chain of two polynomials P and Q Mod alphas.
*/	
AltArr_t* modLastNonZeroChain_AA (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, int isDeepMod);

// Note: k == -1 returns the last non-zero subresutlant.
AltArr_t* modKthSubresultantChain_AA (degree_t k, AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, int isDeepMod);


///////////////////////////////////
///// Triangular Set Normalization
///////////////////////////////////

/**
 * Assume that the iterated initials of P are regular w.r.t T,
 * T is a zero-dimensional triangular set.
 * This algorithm returns (r, q) s.t. q is invertible module <T>,
 * qp = r mod <T>, and r is normalized w.r.t T.
 * r is normalized w.r.t T if and only if r \in Q or 
 * mvar(r) not \in mvar(T) and init(r) is normalized w.r.t <T>
 */
AltArr_t* normalizePolynomial_AA (AltArr_t* P, AltArr_t** T, AltArr_t** A, int s, int nvar);


//////////////
// GCD 
//////////////

/**
 * GCD of two constant polynomials
 * given two constant polynomials P and Q, this algorithm computes the gcd
 * of them and return a constant polynomial.
*/
AltArr_t* integralGCD_AA_polyOut (AltArr_t* P, AltArr_t* Q);

/**
 * GCD of a polynomial and an Integer
 * given  polynomial P and integer c, this algorithm computes the gcd
 * of them and return a constant polynomial.
*/
AltArr_t* gcd_AA_Q (AltArr_t* P, mpq_t c);

/**
 * GCD of two polynomials
 * given two polynomials P and Q with the same main variable, 
 * this algorithm computes the gcd of them w.r.t the main 
 * variable and return a polynomial.
*/
AltArr_t* gcd_AA (AltArr_t* P, AltArr_t* Q);

/**
 * Primitive Factorization w.r.t the main variable
 * given a polynomial, P, this algorithm computes the primitive part
 * and the content of P w.r.t the main variable.
*/
AltArr_t* mainPrimitiveFactorization_AA (AltArr_t* P, AltArr_t** cont);

/**
 * Primitive Part w.r.t the main variable
 * given a polynomial, P, this algorithm computes the primitive part
 * of P w.r.t the main variable.
*/
AltArr_t* mainPrimitivePart_AA (AltArr_t* P, int mvar);

/**
 * Content w.r.t the main variable
 * given a polynomial, P, this algorithm computes the content
 * of P w.r.t the main variable.
*/
AltArr_t* mainContent_AA (AltArr_t* P);

	
/******************
 * Square Free 
 *****************/

/**
 * Computes the square free part of a polynomial.
 */
AltArr_t* squareFreePart_AA (AltArr_t* aa, int nvar);

	
#ifdef __cplusplus
}
#endif


#endif

