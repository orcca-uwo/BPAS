

#ifndef _SMQP_SUPPORT_UNPK_H_
#define _SMQP_SUPPORT_UNPK_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ExponentPackingDefines.h"
#include <gmp.h>

//forward declares
typedef struct AAElem AAElem_t;
typedef struct AltArr AltArr_t;
typedef struct AltArrDegList AltArrDegList_t;
typedef mpq_t ratNum_t;
typedef struct ProductHeap_AA ProductHeap_AA;
typedef struct ProductHeapChain_AA ProductHeapChain_AA;

/********************************
 * Declarations of all unpacked methods for AltArr_t.
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
 * Make an unpacked AltArr_t with allocSize elements.
 * The array of unpacked exponent vectors is returned in elems[0].degs. 
 * NOTE: the following elems do no have their degree pointers set.
 */
AltArr_t* makePolynomial_AA_unpk(int allocSize, int nvar);

void freePolynomial_AA_unpk(AltArr_t* aa);

void unpackExponentVectors_AA_inp(AltArr_t* aa);

void unpackExponentVector(degrees_t packedExp, degree_t* unpackedExp, int nvar);

void tryPackExponentVectors_AA_inp(AltArr_t* aa);

void packExponentVectors_AA_inp(AltArr_t* aa);

int monomialDivideTest_AA_unpk(AltArr_t* a, int idxa, AltArr_t* b, int idxb);

void printDegs_AA_unpk(FILE* fp, degrees_t degs, const char** syms, int nvar);

void printPoly_AA_unpk(FILE* fp, const AltArr_t* aa, const char** syms, int nvar);

int isExactlyEqual_AA_unpk(AltArr_t* a, AltArr_t* b);

void nonZeroVariables_AA_unpk(AltArr_t* a, int* vars);

degree_t totalDegree_AA_unpk(AltArr_t* aa);

degree_t partialDegree_AA_unpk(const AltArr_t* aa, int k);

void partialDegrees_AA_unpk(const AltArr_t* aa, degree_t* degsList);

degree_t mainDegree_AA_unpk(AltArr_t* aa);

int mainVariable_AA_unpk(AltArr_t* aa);

void coefficient_AA_unpk(AltArr_t* aa, const degree_t* degs, int nvar, mpq_t retCoef);

void setCoefficient_AA_unpk(AltArr_t* aa, const degree_t* degs, int nvar, const mpq_t coef);

/**
 * Append a new term ot the end of the AltArr_t aa. 
 * This operations is potentially unsafe; no re-ordering is performed, 
 * it is only a simple append.
 */
void addTerm_AA_unpk(AltArr_t* aa, const degree_t* degs, int nvar, const mpq_t coef);

int isEqualWithVariableOrdering_AA_unpk(AltArr_t* a, AltArr_t* b, const int* xs, int xsSize);

void expandNumVars_AA_unpk(AltArr_t* aa, int newNvar);

void expandNumVarsLeft_AA_unpk(AltArr_t* aa, int newNvar);

/**
 * Re-pack the exponent vectors of aa so they have one less variable. 
 * This re-packing is done such that the exponent at index idx is discarded completely. 
 * NOTE sorting may be needed after a call to this function, depending on the circumstances.
 * No combination of like-terms or re-ordering is done here.
 */
void shrinkNumVarsAtIdx_AA_unpk(AltArr_t* aa, int idx);

void shrinkAndReorderVars_AA_unpk(AltArr_t* aa, int* varMap, int varmapSize);

void reorderVars_AA_unpk(AltArr_t* aa, int* varMap, int varMapSize);

void setDegrees_AA_inp_unpk(AltArr_t* aa, int idx, const degree_t* degsList, int nvar);

/**
 * Resizes the polynomial and the unpacked exponent vector array.
 */
void resizePolynomial_AA_unpk(AltArr_t* aa, int allocSize);

/**
 * Deep copy an AltArr_t and unpack if not already.
 */
AltArr_t* deepCopyPolynomial_AA_unpk(AltArr_t* aa);

AltArrDegList_t* deepCopyPolynomial_AADegListFromAA_unpk(AltArr_t* aa);

void mergeSortPolynomial_AA_unpk(AltArr_t* aa);

AltArr_t* sortPolynomial_AA_unpk(AltArr_t* aa);

void condensePolyomial_AA_unpk(AltArr_t* aa);

extern void negatePolynomial_AA(AltArr_t* aa);
static inline void negatePolynomial_AA_unpk(AltArr_t* aa) {
	negatePolynomial_AA(aa); //this doesn't touch exponents.
}

extern void multiplyByRational_AA_inp(AltArr_t* aa, const mpq_t z);
static inline void multiplyByRational_AA_inp_unpk(AltArr_t* aa, const mpq_t z) {
	multiplyByRational_AA_inp(aa, z);
}

void evalPolyToVal_AA_unpk(const AltArr_t* aa, ratNum_t* vals, short nvar, ratNum_t res);

AltArr_t* evaluatePoly_AA_unpk(AltArr_t* aa, int* active, ratNum_t* vals, short nvar);




/**
 * (right) shift main variable of polynomial aa of size n.
 */ 
AltArr_t* mainLShiftPolynomial_AA_unpk (AltArr_t* aa, int n);
AltArr_t* mainLShiftPolynomial_AA_inp_unpk (AltArr_t* aa, int n);   

/**
 * Given polynomial a, return the leading term
 */
AltArr_t* leadingTerm_AA_unpk (AltArr_t* aa, int nvar);

/**
 * Get the leading variable, that is, the highest-order variable with positive degree
 * of this polynomial.
 * returns the postion of leading variable, -1 (for constant polynomials) and -2 (for zero polynomials).
 */
int leadingVariable_AA_unpk (AltArr_t* aa);

/**
 * Main Leading Degree
 * given a sorted multivariable polynomial aa, this function returns 
 * the maximum value of the main variable of aa.
 */ 
int mainLeadingDegree_AA_unpk (AltArr_t* aa);

/**
 * Main Leading Coefficient
 * given a multivariable polynomial aa, this function returns 
 * the leading coefficient of aa w.r.t the main variable.
 * NOTE: This function is NOT inplace! 
 */ 
AltArr_t* mainLeadingCoefficient_AA_unpk (AltArr_t* aa);

/**
 * Main Leading Coefficient w.r.t the e-th variable
 * given a multivariable polynomial aa and index of a variable,
 * this function returns the leading coefficient of aa w.r.t 
 * the special variable.
 * NOTE: This function is NOT inplace! 
 */ 
AltArr_t* mainCoefficientAtIdx_AA_unpk (AltArr_t* aa, int e);

/**
 * given two polynomials, return polynomial with bigger leading term
 * Note: return polynomial is a deep copy of the polynomial
 */
AltArr_t* maxPolynomials_AA_unpk (AltArr_t* a, AltArr_t* b);

/**
 * given two polynomials, return polynomial with bigger leading term
 * Note: Is is done INPLACE w.r.t. the first input, a.
 */
AltArr_t* maxPolynomials_AA_inp_unpk (AltArr_t* a, AltArr_t* b);





///     Arithmetic                    /////////////////////

void addRationalNumber_AA_inp_unpk(AltArr_t* aa, const mpq_t coef);

AltArr_t* addPolynomials_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar);

AltArr_t* subPolynomials_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar);
    
AltArr_t* addPolynomials_AA_inp_unpk(AltArr_t* a, AltArr_t* b, int nvar);

AltArr_t* subPolynomials_AA_inp_unpk(AltArr_t* a, AltArr_t* b, int nvar);  


///     Multiplication                /////////////////////

void prodheapInsert_AA_unpk(ProductHeap_AA* h, ProductHeapChain_AA* chain, register degrees_t degs);

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * NOTE it is only valid ot call RemoveMax if prodheapPeek_AA returns non-NULL.
 */
ProductHeapChain_AA* prodheapRemoveMax_AA_unpk(ProductHeap_AA* h);

AltArr_t* multiplyPolynomials_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar);



///     Division                     /////////////////////

void divideBySingleTerm_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int nvar);

void dividePolynomials_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, register int nvar);

void exactDividePolynomials_AA_unpk (AltArr_t* c, AltArr_t* b, AltArr_t** res_a, register int nvar);

int divideTestSingleTerm_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, int nvar);

/**
 * Determine if a polynomial, c, is exactly divisble by another, b.
 * If the division is exact, returns the quoteint in res_a. 
 * Otherwise, the value of res_a is undefined.
 * returns 1 iff division is exact. 0 otherwise.
 */
int divideTest_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, int nvar);

/**
 * Given two polynomials, c and b, to compute the remainder and quotient
 * of leading terms of c and b such that lt(c) = lt(b)*res_a + res_r.   
 */
void divideByLeadingTerms_AA_unpk (AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int nvar);


/*****************
 * Derivative / Integral
 *****************/

/**
 * Get the k-th partial derivative of aa for the variable at index idx.
 *
 * returns the partial derivative
 */
AltArr_t* derivative_AA_unpk(AltArr_t* aa, int idx, int k);


/**
 * Get the k-th partial integral of aa for the variable at index idx;
 *
 * Note: constants of integration are not included.
 *
 * returns the partial integral.
 */
AltArr_t* integral_AA_unpk(AltArr_t* aa, int idx, int k);

AltArr_t* integrateExpand_AA_unpk(AltArr_t* aa, int k);

/**
 * Integrate a polynomial in place assuming it's exponent
 * vector has already been expanded to accomodate the new variable.
 */
void integrateExpand_AA_inp_unpk(AltArr_t* aa, int k);



/*****************
 * Content, PrimitivePart, etc.
 *****************/

AltArr_t* univariateGCD_AA_unpk(AltArr_t* a, AltArr_t* b);

/**
 * Given a polynomial, find a monic factor common to all terms of the polynomial.
 * If factored is not NULL then the input polynomial with common factor removed
 * is returned in factored.
 *
 * returns the common factor. 
 */
AltArr_t* commonFactor_AA_unpk(AltArr_t* a, AltArr_t** factored);

void univarEvaluate_AA_unpk(AltArr_t* aa, const mpq_t point, mpq_t res);

#ifdef __cplusplus
}
#endif

#endif

