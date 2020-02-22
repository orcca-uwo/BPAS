
#ifndef _SMZP_SUPPORT_RECURSIVE_H_
#define _SMZP_SUPPORT_RECURSIVE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "SMZP_Support.h"
#include "SMZP_Support_Recursive_Unpacked.h"
#include "../RationalNumberPolynomial/SMQP_Support_Recursive-AA.h"
#include "DUZP_Support.h"
#include "../Utils/Unix_Timer.h"

typedef struct RecArrElemZ {
	int coefSize;
	degree_t exp;
	AAZElem_t* coef;
} RecArrElemZ_t;

typedef struct RecArrZ {
	int alloc;
	int size;
	int unpacked;
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
 * Free a RecArr_t and it's RecArrElem_t array. as well as the underlying coef data.
 */
static inline void freeRecArrayZAndCoef(RecArrZ_t* rArr) {
	if (rArr != NULL) {
		if (rArr->alloc > 0) {
			if (rArr->unpacked) {
				free( (degree_t*) rArr->elems->coef->degs);
			}
			free(rArr->elems->coef);
			free(rArr->elems);
		}
		free(rArr);
		if (rArr->origAA != NULL) {
			free(rArr->origAA);
		}
	}
}


AltArrZ_t* buildRandomAltArrZPoly_MvarSparse(int univarSparsity, int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg, time_t seed);

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
AltArrZ_t* recProdHeapGetNextCoef_AAZ(RecProdHeap_AA* h, const RecArrElemZ_t* __restrict__ aElems, const RecArrElemZ_t* __restrict__ bElems, const int* aUnpacked, int bUnpacked);

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
// Multi-divisor (Pseudo-)division
//////////////////////////////////

/** 
 * Do the pseudo division of c by b when the idx-th member
 * of exponent vectors is the main variable. 
 * Quotient is returned in res_a, and remainder in res_r. 
 * If e is not NULL then it returns the exact number of division steps which occurred.
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 */
void pesudoDivideAtIdx_AAZ (int idx, AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, AltArrZ_t** hPow, int nvar, int lazy);

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
void multiDivisorPseudoDivide_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet);

/** 
 * Do the pseudo division of c by the triangular-set of B in the naive principle
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 */
void naiveMultiDivisorPseudoDivide_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet);


/** 
 * Do the pseudo division of c by the normalized triangular-set of B in the recursive principle 
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 */
void normalizedTriangularSetPseudoDivide_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t*** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet);

/** 
 * Multi-Divisor Division (MDD) where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
void recursiveTriangularSetMDD_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t*** quoSet, AltArrZ_t** rem, int nvar, int nSet);


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

static inline int Log2_AAZ (int n)
{
    return (n > 1) ? 1 + Log2_AAZ (n/2) : 0;
}

/** 
 * free an AltArrsZ_t 
 */
static inline void freeAltArrsZ (AltArrsZ_t* AAs)
{
	while (AAs == NULL){
		return;
	}
	AltArrsZ_t* cur = AAs;
	while (cur != NULL){
		if (cur->poly != NULL){
			freePolynomial_AAZ (cur->poly);
		}
		cur = cur->next;
	}
	free (AAs);
}
    
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
AltArrZ_t* CFDucosOptZ (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd, int flag);
AltArrZ_t* CFDucosOptZ_new (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd, int flag);

AltArrZ_t* mainTailPolynomial_AAZ (AltArrZ_t* a, AltArrZ_t** lc);

/**
 * Ducos SubresultantChain (L. Docus, Optimizations of the Subresultant Algorithm, p8)  
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * and store them on linked-list, AltArrsZ_t.
 * SC is the resultantchain and len is the size.
*/
void DucosSubresultantChainZ_rev (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len, int type);
void DucosSubresultantChainZ (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len);

/**
 * Ducos SubresultantChain At Idx
 * given two polynomials P and Q, and idx, this algorithm computes the idx-th and 
 * idx+1-th resultantchain
 * Ex: if idx=0: SC_idx = resultant(P,Q) and SC_idx1 = upper polynomial in the chain.
 * @note Our motivation is making use of speculative subresultant chain algorithms to get better performance. 
 */
void DucosSubresultantChainAtIdxZ (AltArrZ_t* P, AltArrZ_t* Q, int idx, AltArrZ_t** SC_idx, AltArrZ_t** SC_idx1);

/**
 * Ducos SubresultantChain At Idx
 * given two polynomials P and Q, and idx, this algorithm computes the idx-th and 
 * idx+1-th resultantchain
 * Ex: if idx=0: SC_idx = resultant(P,Q) and SC_idx1 = upper polynomial in the chain.
 * @note Our motivation is making use of speculative subresultant chain algorithms to get better performance. 
 * @note this function returns a list of mainLeadningCoefs of subresultants before idx, stored in principle_coefs of size pcSize.
 */
void DucosSubresultantChainAtIdxZ_withPrincipleCoefs (AltArrZ_t* P, AltArrZ_t* Q, int idx, AltArrZ_t** SC_idx, AltArrZ_t** SC_idx1, AltArrsZ_t** principle_coefs, int* pcSize);


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
 * Return the last non-zero polynomial in subresultant chain of two polynomials P and Q.
*/	
AltArrZ_t* lastNonZeroChain_AAZ (AltArrZ_t* P, AltArrZ_t* Q);


///////////////////////////////
// Extended Subresultant Chain 
///////////////////////////////

/**
 * To store the Bezout Coefficients and partial gcds
 */
typedef struct exgcdsZ {
	AltArrZ_t* r;
	AltArrZ_t* a;
	AltArrZ_t* b;
	struct exgcdsZ* next;
} exgcdsZ_t;

/**
 * Free exgcdsZ_t
 */
static inline void freeExgcds_AAZ (exgcdsZ_t* lg) {
	if (lg != NULL) {
		exgcdsZ_t* cur = lg;
		while (cur != NULL) {
	 		if (cur->r != NULL) {
 				freePolynomial_AAZ (cur->r);
 			}
	 		if (cur->a != NULL) {
 				freePolynomial_AAZ (cur->a);
 			}
		 	if (cur->b != NULL) {
	 			freePolynomial_AAZ (cur->b);
	 		}
	 		cur = cur->next;
	 	}
	}
}

/**
 * The semi-lazard optimization based on "L. Docus, Optimizations of the Subresultant Algorithm".
 * used in exLazardOpt_AA algorithm
 */
AltArrZ_t* semiLazardOpt_AAZ  (AltArrZ_t* Sd, AltArrZ_t* Sdm, AltArrZ_t* s);
AltArrZ_t* exLazardOpt_AAZ (AltArrZ_t* Sd, exgcdsZ_t* VSdm, AltArrZ_t* s, AltArrZ_t** hc, AltArrZ_t** qc);

/**
 * Extended Ducos Optimization  
 * Input A, VSdm(= [Sdm, VSdm->a, VSdm->b]), Se, sd
 * Output Sem, and return by reference hc and qc. 
*/
AltArrZ_t* exDucosOpt_AAZ (exgcdsZ_t* VSd, exgcdsZ_t* VSdm, AltArrZ_t* Se, AltArrZ_t* sd, AltArrZ_t** hb, AltArrZ_t** qb);

/**
 * Extended Ducos SubresultantChain  
 * given two polynomials P and Q, this algorithm computes the resultantchain
 * with all partial Bezout Coefficients and store them on linked-list, exgcdsZ_t.
 * SC is the resultantchain and len is the size.
*/
void exDucosSubresultantChain_rev_AAZ (AltArrZ_t* P, AltArrZ_t* Q, exgcdsZ_t** SC, int* len, int type);
void exDucosSubresultantChain_AAZ (AltArrZ_t* P, AltArrZ_t* Q, exgcdsZ_t** SC, int* len);

/**
 * Extended Resultant of two polynomials by DucosSubresultantChain
 * given two polynomials P and Q, this algorithm computes the resultant
 * with its Bezout Coefficients by calling exDucosSubresultantChain_AA algorithm.
*/
AltArrZ_t* exDucosResultant_AAZ (AltArrZ_t* P, AltArrZ_t* Q, AltArrZ_t** a, AltArrZ_t** b);


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
AltArrZ_t* normalizePolynomial_AAZ (AltArrZ_t* P, AltArrZ_t** T, AltArrZ_t** A, int s, int nvar);

////////
// GCD
////////

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

    
/******************
 * Square Free 
 *****************/

/**
 * Computes the square free part of a polynomial.
 */
AltArrZ_t* squareFreePart_AAZ (AltArrZ_t* aa, int nvar);


/** 
 * Compute the square free factorization of aa.
 * The factorization is returned in two parts, u is the content of aa
 * and facts is a list of the factors such that facts[i] is a primitive square free factor
 * of aa with exponent exps[i].
 * If facts_p points to NULL then an array of factors is allocated, otherwise, it 
 * is assumed to be a pre-allocated array with size *nfacts; the same is true for exps_p.
 *
 * @note the exponents returned may not be unique since partial and trivial factorizations are
 *       performed here as well.
 *
 * @param aa the polynomial to factorization
 * @param[out] u the "content" of aa
 * @param[out] facts_p a pointer in which the array of factors is returned. 
 * @param[out] exps_p a pointer in which the array of exponents is returned.
 * @param[in,out] nfacts a pointer in which the number of factors/exponents is returned, 
 *                on entry if *facts_p or *exps_p is not NULL, *nfacts is the size of the pre-allocated array.
 *
 */
void squareFree_AAZ(const AltArrZ_t* aa, mpz_t u, AltArrZ_t*** facts_p, degree_t** exps_p, int* nfacts);


void biModularSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size);

#ifdef __cplusplus
}
#endif


#endif

