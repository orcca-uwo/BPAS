

#ifndef _SMZP_SUPPORT_UNPK_H_
#define _SMZP_SUPPORT_UNPK_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "../RationalNumberPolynomial/ExponentPackingDefines.h"
#include "../RationalNumberPolynomial/SMQP_Support-AA.h"
#include <gmp.h>

//forward declares
typedef struct AAZElem AAZElem_t;
typedef struct AltArrZ AltArrZ_t;
typedef struct AltArrZDegList AltArrZDegList_t;
typedef struct ProductHeap_AAZ ProductHeap_AAZ;
typedef struct ProductHeapChain_AAZ ProductHeapChain_AAZ;

/********************************
 * Declarations of all unpacked methods for AltArrZ_t.
 * One should avoid inline definitions here due to
 * the intermingling of code here and SMQP_Support.h
 * and SMQP_Support_Packed.h
 ********************************/

/********************************
 * This stuff is very tricky...
 * Please note the following assumptions:
 *  - the degs array is always allocated to match aa->nvar*aa->alloc
 *  - aa->elems[i].degs actually points to the degs array only if i < aa->size;
 *  - if aa->nvar = 0 then aa->elems->degs CAN be NULL. For positive nvar it should NEVER be NULL
 *    the single exception is when aa->alloc = 0 (but if this happens we have other problems).
 ********************************/

/**
 * Make an unpacked AltArrZ_t with allocSize elements.
 * The array of unpacked exponent vectors is returned in elems[0].degs. 
 * NOTE: the following elems do no have their degree pointers set.
 */
AltArrZ_t* makePolynomial_AAZ_unpk(int allocSize, int nvar);

void freePolynomial_AAZ_unpk(AltArrZ_t* aa);

void unpackExponentVectors_AAZ_inp(AltArrZ_t* aa);

void tryPackExponentVectors_AAZ_inp(AltArrZ_t* aa);

void packExponentVectors_AAZ_inp(AltArrZ_t* aa);

int monomialDivideTest_AAZ_unpk(AltArrZ_t* a, int idxa, AltArrZ_t* b, int idxb);

void printDegs_AAZ_unpk(FILE* fp, degrees_t degs, const char** syms, int nvar);

void printPoly_AAZ_unpk(FILE* fp, const AltArrZ_t* aa, const char** syms, int nvar);

int isExactlyEqual_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b);

void nonZeroVariables_AAZ_unpk(AltArrZ_t* a, int* vars);

degree_t totalDegree_AAZ_unpk(const AltArrZ_t* aa);

degree_t partialDegree_AAZ_unpk(const AltArrZ_t* aa, int k);

void partialDegrees_AAZ_unpk(const AltArrZ_t* aa, degree_t* degsList);

degree_t mainDegree_AAZ_unpk(AltArrZ_t* aa);

int mainVariable_AAZ_unpk(AltArrZ_t* aa);

void coefficient_AAZ_unpk(AltArrZ_t* aa, const degree_t* degs, int nvar, mpz_t retCoef);

void setCoefficient_AAZ_unpk(AltArrZ_t* aa, const degree_t* degs, int nvar, const mpz_t coef);

int isEqualWithVariableOrdering_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b, const int* xs, int xsSize);

void expandNumVars_AAZ_unpk(AltArrZ_t* aa, int newNvar);

void expandNumVarsLeft_AAZ_unpk(AltArrZ_t* aa, int newNvar);

/**
 * Re-pack the exponent vectors of aa so they have one less variable. 
 * This re-packing is done such that the exponent at index idx is discarded completely. 
 * NOTE sorting may be needed after a call to this function, depending on the circumstances.
 * No combination of like-terms or re-ordering is done here.
 */
void shrinkNumVarsAtIdx_AAZ_unpk(AltArrZ_t* aa, int idx);

void shrinkAndReorderVars_AAZ_unpk(AltArrZ_t* aa, int* varMap, int varmapSize);

void reorderVars_AAZ_unpk(AltArrZ_t* aa, int* varMap, int varMapSize);

void setDegrees_AAZ_inp_unpk(AltArrZ_t* aa, int idx, const degree_t* degsList, int nvar);

void setExponentTerm_AAZ_inp_unpk(AltArrZ_t* aa, int idx, degree_t deg, int k);


/**
 * Resizes the polynomial and the unpacked exponent vector array.
 */
void resizePolynomial_AAZ_unpk(AltArrZ_t* aa, int allocSize);

/**
 * Deep copy an AltArrZ_t and unpack if not already.
 */
AltArrZ_t* deepCopyPolynomial_AAZ_unpk(const AltArrZ_t* aa);

void deepCopyPolynomial_AAZ_inp_unpk(const AltArrZ_t* aa, AltArrZ_t** bb);

AltArrZDegList_t* deepCopyPolynomial_AAZDegListFromAA_unpk(AltArrZ_t* aa);

AltArr_t* deepCopyPolynomial_AAFromAAZ_unpk(AltArrZ_t* aa);

AltArrZ_t* deepCopyPolynomial_AAZFromAA_unpk(AltArr_t* aa);

void mergeSortPolynomial_AAZ_unpk(AltArrZ_t* aa);

AltArrZ_t* sortPolynomial_AAZ_unpk(AltArrZ_t* aa);

void condensePolynomial_AAZ_unpk(AltArrZ_t* aa);

void removeZeroTerms_AAZ_unpk(AltArrZ_t* aa);

/**
 * Given some polynomial which is a collection of coefficient-monomial pairs, 
 * transform it into canonical form. This includes removing pairs with
 * a zero coefficeint, sorting the terms by monomial, and ensuring
 * each monomial is unique.
 */
static inline void canonicalizePolynomial_AAZ_unpk(AltArrZ_t* aa) {
	removeZeroTerms_AAZ_unpk(aa);
	mergeSortPolynomial_AAZ_unpk(aa);
	condensePolynomial_AAZ_unpk(aa);
}

extern void negatePolynomial_AAZ(AltArrZ_t* aa);
static inline void negatePolynomial_AAZ_unpk(AltArrZ_t* aa) {
	negatePolynomial_AAZ(aa); //this doesn't touch exponents.
}

extern void multiplyByRational_AAZ_inp(AltArrZ_t* aa, const mpz_t z);
static inline void multiplyByRational_AAZ_inp_unpk(AltArrZ_t* aa, const mpz_t z) {
	multiplyByRational_AAZ_inp(aa, z);
}

void evalPolyToVal_AAZ_unpk(const AltArrZ_t* aa, const mpz_t* vals, short nvar, mpz_t res);

AltArrZ_t* evaluatePoly_AAZ_unpk(const AltArrZ_t* aa, const int* active, const mpz_t* vals, short nvar);




/**
 * (right) shift main variable of polynomial aa of size n.
 */ 
AltArrZ_t* mainLShiftPolynomial_AAZ_unpk (AltArrZ_t* aa, int n);
AltArrZ_t* mainLShiftPolynomial_AAZ_inp_unpk (AltArrZ_t* aa, int n);   

/**
 * Given polynomial a, return the leading term
 */
AltArrZ_t* leadingTerm_AAZ_unpk (AltArrZ_t* aa, int nvar);

/**
 * Get the leading variable, that is, the highest-order variable with positive degree
 * of this polynomial.
 * returns the postion of leading variable, -1 (for constant polynomials) and -2 (for zero polynomials).
 */
int leadingVariable_AAZ_unpk (AltArrZ_t* aa);

/**
 * Main Leading Degree
 * given a sorted multivariable polynomial aa, this function returns 
 * the maximum value of the main variable of aa.
 */ 
int mainLeadingDegree_AAZ_unpk (AltArrZ_t* aa);

/**
 * Main Leading Coefficient
 * given a multivariable polynomial aa, this function returns 
 * the leading coefficient of aa w.r.t the main variable.
 * NOTE: This function is NOT inplace! 
 */ 
AltArrZ_t* mainLeadingCoefficient_AAZ_unpk (AltArrZ_t* aa);

/**
 * Main Leading Coefficient w.r.t the e-th variable
 * given a multivariable polynomial aa and index of a variable,
 * this function returns the leading coefficient of aa w.r.t 
 * the special variable.
 * NOTE: This function is NOT inplace! 
 */ 
AltArrZ_t* mainCoefficientAtIdx_AAZ_unpk (AltArrZ_t* aa, int e);

/**
 * given two polynomials, return polynomial with bigger leading term
 * Note: return polynomial is a deep copy of the polynomial
 */
AltArrZ_t* maxPolynomials_AAZ_unpk (AltArrZ_t* a, AltArrZ_t* b);

/**
 * given two polynomials, return polynomial with bigger leading term
 * Note: Is is done INPLACE w.r.t. the first input, a.
 */
AltArrZ_t* maxPolynomials_AAZ_inp_unpk (AltArrZ_t* a, AltArrZ_t* b);





///     Arithmetic                    /////////////////////

void addInteger_AAZ_inp_unpk(AltArrZ_t* aa, const mpz_t coef);

AltArrZ_t* addPolynomials_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar);

AltArrZ_t* subPolynomials_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar);
    
AltArrZ_t* addPolynomials_AAZ_inp_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar);

AltArrZ_t* subPolynomials_AAZ_inp_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar);

void subPolynomials_AAZ_inpRHS_unpk(const AltArrZ_t* a, AltArrZ_t** bb);  

AltArrZ_t* CFDucosOptZ_subPolynomials_AAZ_inp_unpk (AltArrZ_t* a, AltArrZ_t* b, int nvar, AltArrZ_t** Pe, int e);


///     Multiplication                /////////////////////

void prodheapInsert_AAZ_unpk(ProductHeap_AAZ* h, ProductHeapChain_AAZ* chain, register degrees_t degs);

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * NOTE it is only valid ot call RemoveMax if prodheapPeek_AAZ returns non-NULL.
 */
ProductHeapChain_AAZ* prodheapRemoveMax_AAZ_unpk(ProductHeap_AAZ* h);

AltArrZ_t* multiplyPolynomials_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b, int nvar);

void multiplyPolynomialsPreAlloc_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b, AltArrZ_t** cc);


///     Division                     /////////////////////

void divideBySingleTerm_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar);

void dividePolynomials_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, register int nvar);

void exactDividePolynomials_AAZ_unpk (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, register int nvar);

void univariatePseudoDivideBySingleTerm_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy);

void univariatePseudoDividePolynomials_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy);

int divideTestSingleTerm_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar);

/**
 * Determine if a polynomial, c, is exactly divisble by another, b.
 * If the division is exact, returns the quoteint in res_a. 
 * Otherwise, the value of res_a is undefined.
 * returns 1 iff division is exact. 0 otherwise.
 */
int divideTest_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar);

/**
 * Given two polynomials, c and b, to compute the remainder and quotient
 * of leading terms of c and b such that lt(c) = lt(b)*res_a + res_r.   
 */
void divideByLeadingTerms_AAZ_unpk (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar);


/*****************
 * Derivative / Integral
 *****************/

/**
 * Get the k-th partial derivative of aa for the variable at index idx.
 *
 * returns the partial derivative
 */
AltArrZ_t* derivative_AAZ_unpk(const AltArrZ_t* aa, int idx, int k);

/**
 * Turn the polynomial pointed to by aa into its k'th partial derivative 
 * with respect to the varibale at index idx.
 *
 * @param[in,out] aa a pointer to the polynomial whose partial derivative is to be computed, also the place where the resulting derivative is stored.
 * @param idx the index of the variable with respect to which the partial derviative is taken.
 * @param k the order of the derivative. 
 *
 */
void derivative_AAZ_inp_unpk(AltArrZ_t** aa, int idx, int k);



/*****************
 * Content, PrimitivePart, etc.
 *****************/

AltArrZ_t* primitivePartAndContent_AAZFromAA_unpk(AltArr_t* aa, mpq_t cont);

AltArrZ_t* univariateGCD_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b);


/**
 * Given a polynomial, find a monic factor common to all terms of the polynomial.
 * If factored is not NULL then the input polynomial with common factor removed
 * is returned in factored.
 *
 * returns the common factor. 
 */
AltArrZ_t* commonFactor_AAZ_unpk(const AltArrZ_t* a, AltArrZ_t** factored);

void univarEvaluate_AAZ_unpk(AltArrZ_t* aa, const mpz_t point, mpz_t res);


#ifdef __cplusplus
}
#endif

#endif

