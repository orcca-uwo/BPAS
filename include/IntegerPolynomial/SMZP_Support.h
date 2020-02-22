
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
#include "SMZP_Support_Unpacked.h"

#define DUZP_INFNORM_SUBRES_THRESHOLD 256
#define DUZP_MINDEG_SUBRES_THRESHOLD 50
#define DBZP_INFNORM_SUBRES_THRESHOLD 5
#define DBZP_MINDEGY_SUBRES_THRESHOLD 5
// #define PRINT_WHICH_SUBRES 1
#define SUBRES_TIME_DEBUG 0

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
	int unpacked;
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

		if (aa->unpacked) {
			degree_t* degs_unpk = (degree_t*) aa->elems[0].degs;
			free(degs_unpk);
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
	newAA->unpacked = 0;
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
	newAA->unpacked = 0;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*allocSize);
	mpz_init(newAA->elems->coef);
	mpz_set(newAA->elems->coef, coef);
	newAA->elems->degs = 0;
	return newAA;
}

static inline AltArrZ_t* makeConstIntPolynomial_AAZ(int allocSize, int nvar, long int coef) {
	if (allocSize < 1) {
		return NULL;
	}

	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	newAA->size = 1;
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->unpacked = 0;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*allocSize);
	mpz_init_set_si(newAA->elems->coef, coef);
	newAA->elems->degs = 0ll;
	return newAA;
}

static inline int isZero_AAZ(const AltArrZ_t* aa) {
    if (aa == NULL || aa->size == 0) {
    	return 1;
    }

    int ret = mpz_sgn(aa->elems[0].coef) == 0;
    if (aa->unpacked) {
    	ret = ret && isZeroExponentVector_unpk(aa->elems[0].degs, aa->nvar);
    } else {
    	ret = ret && isZeroExponentVector(aa->elems[0].degs);
	}

	return ret;
}

static inline int isOne_AAZ(AltArrZ_t* aa) {
	if (aa != NULL) {
		int ret = aa->unpacked ? isZeroExponentVector_unpk(aa->elems[0].degs, aa->nvar) : isZeroExponentVector(aa->elems[0].degs);
		ret = ret && mpz_cmp_ui(aa->elems[0].coef, 1ul) == 0;
		return ret;
    }
    return 0;
}

static inline int isNegativeOne_AAZ(AltArrZ_t* aa) {
	if (aa != NULL) {
		int ret = aa->unpacked ? isZeroExponentVector_unpk(aa->elems[0].degs, aa->nvar) : isZeroExponentVector(aa->elems[0].degs);
		ret = ret && mpz_cmp_si(aa->elems[0].coef, -1l) == 0;
		return ret;
    }
    return 0;
}

static inline int isConstant_AAZ(const AltArrZ_t* aa) {
    if (aa == NULL || aa->size == 0) {
        return 1;
    }
    int isZeroDegs = aa->unpacked ? isZeroExponentVector_unpk(aa->elems->degs, aa->nvar) : isZeroExponentVector(aa->elems->degs);
    if (isZeroDegs) {
        if (mpz_sgn(aa->elems->coef) >= 0) {
            return 1;
        } else {
            return -1;
        }
    }
    return 0;
}

/**
 * Get the partial degree of the k-th variable of the term at index idx of aa.
 * returns that partial degree or -1 if k >= nvar or idx > aa->size
 */
static inline degree_t partialDegreeTerm_AAZ(const AltArrZ_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return -1;
	}
	if (idx >= aa->size || k >= aa->nvar) {
		return -1;
	}
	if (aa->unpacked) {
		return ((degree_t*) aa->elems[idx].degs)[k];
	} else {
		const degrees_t* masks = getExpMaskArray(aa->nvar);
		const int* sizes = getExpOffsetArray(aa->nvar);
		degree_t ret = GET_NTH_EXP(aa->elems[idx].degs, masks[k], sizes[k]);
		return ret;
	}
}

static inline void partialDegreesTerm_AAZ(const AltArrZ_t* aa, int idx, degree_t* degs) {
	if (aa == NULL || aa->size == 0 || idx >= aa->size) {
		return;
	}
	if (aa->unpacked) {
		for (int k = 0; k < aa->nvar; ++k ) {
			degs[k] = ((degree_t*) aa->elems[idx].degs)[k];
		}
	} else {
		const degrees_t* masks = getExpMaskArray(aa->nvar);
		const int* sizes = getExpOffsetArray(aa->nvar);
		for (int k = 0; k < aa->nvar; ++k ) {
			degs[k] = GET_NTH_EXP(aa->elems[idx].degs, masks[k], sizes[k]);
		}
	}
}

void printPoly_AAZ(FILE* fp, const AltArrZ_t* aa, const char** syms, int nvar);

/**
 * Determine if the AAZ aa is in canonical form.
 * This is helpful for debugging.
 * returns non-zero iff in canonical form.
 */
int isInOrder_AAZ(const AltArrZ_t* aa);

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
int isExactlyEqual_AAZ(const AltArrZ_t* a, const AltArrZ_t* b);

/**
 * Get the number of non-zero variabes in aa.
 */
void nonZeroVariables_AAZ(AltArrZ_t* aa, int* vars);

/**
 * Get the total degree of aa.
 */
degree_t totalDegree_AAZ(const AltArrZ_t* aa);

/** 
 * Get the partial degree of the kth-variable.
 * returns 0 if aa is zero or k >= nvar.
 */
degree_t partialDegree_AAZ(const AltArrZ_t* aa, int k);

/**
 * Get the partial degrees of each variable in aa. 
 * The partial degree of variable i is returned in degs[i].
 */
void partialDegrees_AAZ(const AltArrZ_t* aa, degree_t* degs);

/**
 * Get the main degree of aa. In particular, the 
 * partial degree of the first variable whose partial
 * degree is positive.
 */
degree_t mainDegree_AAZ(AltArrZ_t* aa);

/**
 * Get the index of the first non-zero variable in aa, the main variable.
 */
int mainVariable_AAZ(AltArrZ_t* aa);

/**
 * Get the coefficient of aa corresponding to the monomial with partial degrees degs.
 * The coefficient is returned in retCoef.
 */
void coefficient_AAZ(AltArrZ_t* aa, const degree_t* degs, int nvar, mpz_t retCoef);

/**
 * Get the coefficient of aa corresponding to the monomial with partial degrees degs.
 * The coefficient is returned in retCoef.
 */
void setCoefficient_AAZ(AltArrZ_t* aa, const degree_t* degs, int nvar, const mpz_t coef);

/**
 * Determine if two polynomials are equal given a mapping between thier variables in xs.
 * variable xs[2*j]-1 in a matches variable xs[2*j+1]-1 in b. If xs[2*j] is 0 then 
 * the variable in b of xs[2*j+1] does not appear in a, and vice versa.
 */
int isEqualWithVariableOrdering_AAZ(AltArrZ_t* a, AltArrZ_t* b, const int* xs, int xsSize);

/**
 * Get the term of a at index idx. Returns the term as a single term polynomial.
 */
AltArrZ_t* termAtIdx_AAZ(AltArrZ_t* a, int idx);

/**
 * Determine if the constant term of a is zero.
 * returns 1 iff the constant term is 0.
 */
int isConstantTermZero_AAZ(AltArrZ_t* a);

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

/**
 * Given a polynomial, aa, compute a variable map which removes all variables
 * with 0 partial degree.
 *
 * @param aa, the polynomial whose partial degrees are checked.
 * @param varMap, the variable map, pre-allocated, of size at least aa->nvar.
 * @return, the number of variables with positive partial degree.
 */
static inline int computeShrinkMap_AAZ(const AltArrZ_t* aa, int* varMap) {
	if (varMap == NULL) {
		return 0;
	}

	int trueNvar = 0;
	degree_t degs[aa->nvar];
	partialDegrees_AAZ(aa, degs);
	for (int i = 0; i < aa->nvar; ++i) {
		if (degs[i] == 0) {
			varMap[i] = -1;
		} else {
			varMap[i] = trueNvar;
			++trueNvar;
		}
	}
	return trueNvar;
}

/**
 * Given two polynomials, aa and bb, compute a variable map which removes all variables
 * with 0 partial degree in both polynomials. Is is assumed aa and bb have the same number
 * of variables.
 *
 * @param aa, one polynomial whose partial degrees are checked.
 * @param bb, one polynomial whose partial degrees are checked.
 * @param varMap, the variable map, pre-allocated, of size at least aa->nvar.
 * @return, the number of variables with positive partial degree.
 */
static inline int computeShrinkMapPolys_AAZ(const AltArrZ_t* aa, const AltArrZ_t* bb, int* varMap) {
	if (varMap == NULL) {
		return 0;
	}
	if (aa == NULL || bb == NULL || aa->nvar != bb->nvar) {
		return 0;
	}

	int trueNvar = 0;
	degree_t degsA[aa->nvar];
	degree_t degsB[aa->nvar];
	partialDegrees_AAZ(aa, degsA);
	partialDegrees_AAZ(bb, degsB);
	for (int i = 0; i < aa->nvar; ++i) {
		if (degsA[i] == 0 && degsB[i] == 0) {
			varMap[i] = -1;
		} else {
			varMap[i] = trueNvar;
			++trueNvar;
		}
	}
	return trueNvar;

}

/**
 * Given a variable mapping, varMap, which was used exclusively to remove variables, 
 * not re-order variables, compute its reverse map.
 *
 * @param originalSize, the number of variables before shrinking
 * @param shrunkSize, the number of variables after shrinking
 * @param varMap, the variable mapping which performed the shrinking, preallocated.
 * @param revMap, the reverse variable map to return, preallocated.
 *
 */
static inline void reverseShrinkMap_AAZ(int originalSize, int shrunkSize, int* varMap, int* revMap) {
	reverseShrinkMap_AA(originalSize, shrunkSize, varMap, revMap);
}

/**
 * Given a polynomial, aa, remove all variables which have 0 partial degree.
 * A variable map is returned describing how variables were moved.
 * varMap[i] = -1 means variable at index i is removed.
 *
 * @param aa, the polynomial from which to remove variables.
 * @param[out] varMap, the variable map causing the shrink, preallocated.
 * 
 * @return the new number of variables.
 *
 */
int tryShrinkVariables_AAZ_inp(AltArrZ_t* aa, int* varMap);

/**
 * Given a polynomial, aa, and a variable map which was used to shrink (but not re-order)
 * its variables, reverse the mapping to expand the number of variables in aa.
 *
 * @param aa, the polynomial whose variables should be expanded.
 * @param originalSize, the number of variables aa should have after expanding.
 * @param varMap, the original variable map used to shrink the number of variables.
 * 
 * @return the new number of variables.
 *
 */
int reverseShrinkVariables_AAZ_inp(AltArrZ_t* aa, int originalSize, int* varMap);

/**
 * Set the exponent vector of aa's term at index idx to the degrees in 
 * degsList.
 * Note that no sorting or canonicalization is performed.
 *
 * @param aa the polynomial so modify
 * @param idx the index of the term to modify
 * @param degsList the list of new exponents.
 * @param nvar the size of degsList which should be the same as aa->nvar.
 */
//TODO rename to setExponentsTerm
void setDegrees_AAZ_inp(AltArrZ_t* aa, int idx, const degree_t* degsList, int nvar);

/**
 * Set the exponent at index k of aa's term at index idx to the value deg.
 * Note that no sorting or canonicalization is performed.
 *
 * @param aa the polynomial so modify.
 * @param idx the index of the term to modify.
 * @param deg the new exponent value.
 * @param k the index of the exponent to modify.
 */
void setExponentTerm_AAZ_inp(AltArrZ_t* aa, int idx, degree_t deg, int k);

static inline void resizePolynomial_AAZ(AltArrZ_t* aa, int allocSize) {
	if (aa != NULL && aa->unpacked) {
		return resizePolynomial_AAZ_unpk(aa, allocSize);
	}
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
AltArrZ_t* deepCopyPolynomial_AAZ(const AltArrZ_t* aa);

/**
 * Create a copy of the polynomial aa and store it in the (possibly)
 * pre-allocated AltArrZ_t pointed to by bb.
 *
 * @param aa the polynomial to be copied
 * @param bb a pointer to a (possibly) pre-allocated AltArrZ_t to store the copy.
 */
void deepCopyPolynomial_AAZ_inp(const AltArrZ_t* aa, AltArrZ_t** bb);

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
void condensePolynomial_AAZ(AltArrZ_t* aa);

/**
 * Remove all terms with a zero coefficient from aa, shrinking its
 * size as needed. If the polynomial is not null but identically 0, 
 * one zero term will remain.
 */
void removeZeroTerms_AAZ(AltArrZ_t* aa);

/**
 * Given some polynomial which is a collection of coefficient-monomial pairs, 
 * transform it into canonical form. This includes removing pairs with
 * a zero coefficeint, sorting the terms by monomial, and ensuring
 * each monomial is unique.
 */
static inline void canonicalizePolynomial_AAZ(AltArrZ_t* aa) {
	removeZeroTerms_AAZ(aa);
	mergeSortPolynomial_AAZ(aa);
	condensePolynomial_AAZ(aa);
}

/**
 * Negate a polynomial in place.
 */
void negatePolynomial_AAZ(AltArrZ_t* a);

static inline void multiplyByInteger_AAZ_inp(AltArrZ_t* aa, const mpz_t z) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	for (int i = 0; i < aa->size; ++i) {
		mpz_mul(aa->elems[i].coef, aa->elems[i].coef, z);
	}
}

static inline void divideByIntegerExact_AAZ_inp(AltArrZ_t* aa, const mpz_t z) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	for (int i = 0; i < aa->size; ++i) {
		mpz_divexact(aa->elems[i].coef, aa->elems[i].coef, z);
	}
}


AltArrZ_t* multiplyByInteger_AAZ(AltArrZ_t* a, mpz_t div);

/**
 * Divide a through by integer divisor d. It is assumed that
 * div divides content(a). 
 */
AltArrZ_t* divideByInteger_AAZ(AltArrZ_t* a, mpz_t div);

/**
 * Apply a modulo to all of the integer coefficients of the polynomial p,
 * modifying the polynomial in place. 
 * 
 * @param p the polynomial which to apply the modulo.
 * @param mod the modulo.
 *
 */
void applyModuloSymmetric_AAZ_inp(AltArrZ_t* p, const mpz_t mod);
 
/**
 * Evaluate a polynomial represented as an alternating array.
 * This method returns the rational number vlaue in res.
 * vals[i] corresponds to the value for variable i in aa.
 */
void evalPolyToVal_AAZ(const AltArrZ_t* aa, const mpz_t* vals, short nvar, mpz_t res);

/**
 * Evaluate a polynomial represented as an alternating array.
 * This method returns another polynomial as not all variables need
 * to have values supplied. But, if they do, a constant term will be returned. 
 * both active and vals are arrays of size nvar. active[i] determines if 
 * the variable at degs[i] is to be evaluated using vals[i] as value.
 */
AltArrZ_t* evaluatePoly_AAZ(const AltArrZ_t* aa, const int* active, const mpz_t* vals, short nvar);


AltArrZ_t* convertFromAAZElemToAAZ (AAZElem_t* coef, int coefSize, int nvar, int unpacked);


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
 * Make a List of Leading Coefficients w.r.t the idx-th variable
 * given a multivariable polynomial aa and index of a variable,
 * this function returns the list of leading coefficients of aa w.r.t 
 * the special variable.
 * NOTE: This function is NOT inplace! 
 */ 
void mainCoefficientListAtIdx_AAZ (AltArrZ_t* aa, int idx, AltArrZ_t*** cList, int *sz);

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

/**
 * Get maximum absolute coefficient in a. Inifinity norm. 
 */
static inline void infinityNorm_AAZ (const AltArrZ_t* a, mpz_t c) {
	if (a == NULL || a->size == 0) {
		mpz_set_ui(c, 0ul);
		return;
	}
	int size = a->size;
	AAZElem_t* elems = a->elems;
	mpz_abs(c, elems[0].coef);
	for (int i = 1; i < size; ++i) {
		if (mpz_cmpabs(c, elems[i].coef) < 0) {
			mpz_abs(c, elems[i].coef);
		}
	}
}

    
/*****************
 * SMQP Addition
 *****************/

void addInteger_AAZ_inp(AltArrZ_t* aa, const mpz_t coef);

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
 * Subtract the polynomial contained in bb from the polynomial a, 
 * returning the difference back in the polynomial contained in bb.
 * This operation re-uses the allocation of the polynomial
 * pointed to by bb to store the difference. 
 *
 * @param a the polynomial to subtract from
 * @param[in,out] a pointer to the polynomial which is subtracted from a, as well as where the difference is stored.
 * @param nvar the number of variables in both a and b.
 * 
 */
void subPolynomials_AAZ_inpRHS(const AltArrZ_t* a, AltArrZ_t** bb);

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

typedef struct ProductHeap_AAZ {
	ProductHeapElem_AAZ* elements;
	polysize_t heapSize;
	polysize_t maxHeapSize;
	int nvar;
#if SMQP_INT_PRODHEAP
	int lastB;
#else
	AAZElem_t* lastB;
#endif
	degree_t* unpackedDegs;
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
	h->unpackedDegs = NULL;
	return h;
}

/**
 * Initialize the product heap with the two polynomials to multiply, a and b.
 *
 * We know the maximum heap size is numTerms(a) as at most one of a_i*b
 * is in the heap at once. 
 */
ProductHeap_AAZ* prodheapInit_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, int nvar);

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
 * Multiply two polynomials a and b, returning their product.
 * This algorithm makes use of heaps as an efficient search data structure. 
 * It is assumed that both a and b have compatible exponent vectors (i.e. the same number of variables).
 * 
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
AltArrZ_t* multiplyPolynomials_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, int nvar);


/**
 * Multiply two polynomials a and b, returning their product in the pointer cc.
 * The pointer cc should (but not necessarily) point to a pre-allocated polynomial
 * whose space is used to store the product.
 * Polynomials a and b should have compatible exponent vectors (i.e. same number of variables).
 * 
 * @note cc should not point to one of a or b.
 *
 * @param a the multiplicand polynomial.
 * @param b the multiplier polynomial.
 * @param[out] cc a pointer to (a pre-allocated polynomial to store) the product polynomial.
 */
void multiplyPolynomialsPreAlloc_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, AltArrZ_t** cc);

/**
 * Multiply two polynomials in-place wrt the first polynomial.
 * That is, the multiplication occurs, reusing elements of a as much as possible
 * and the previous content of a is destroyed in the process.
 * It is assumed that a and b have the same sized exponent vectors. 
 *
 *
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the new location of parameter a.
 *
 */
AltArrZ_t* multiplyPolynomials_AAZ_inp(AltArrZ_t* a, const AltArrZ_t* b, int nvar);

/**
 * Multiply all the polynomials in the list polys, of size n,
 * into a single product. 
 * 
 * @param polys the list of polynomials to multiply.
 * @param n the size of the list polys.
 *
 * @return the product of polynomials in the list polys.
 *
 */
AltArrZ_t* multiplyManyPolynomials_AAZ(AltArrZ_t const * const * polys, int n);

/**
 * Multiply all the polynomials in the list polys, of size n,
 * into a single product. Reuse the space pointed to by prod
 * to store the resulting product.
 * 
 * @param polys the list of polynomials to multiply.
 * @param n the size of the list polys.
 * @param[out] prod a pointer to a potentially allocated polynomial to store the product. 
 *
 */
void multiplyManyPolynomialsPreAlloc_AAZ(AltArrZ_t const* const* polys, int n, AltArrZ_t** prod);

/**
 * Multiply together all the polynomials in the list polys, of size n,
 * exlucding the polynomial at index idx. 
 * The single AltArrZ_t product polynomial is returned.
 *
 * @param polys an array of DUZP_t polys of size n.
 * @param n the size of the array polys.
 * @param idx the index of the polynomial to exclude.
 *
 * @return the product polynomial.
 */
AltArrZ_t* multiplyAllButOnePolynomials_AAZ(AltArrZ_t const*const* polys, int n, int idx);

/** 
 * Multiply polynomial, aa, by X_{idx}^n where 0 <= idx < nvar. 
 * Note this works inplace.
 */ 
void multiplyPolynomialAtIdxByXn_AAZ_inp (AltArrZ_t* aa, int idx, int n, int nvar);


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

void univariatePseudoDivideBySingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy);

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
int divideTestSingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar);
 
/**
 * Given two polynomials, c and b, to compute the remainder and quotient
 * of leading terms of c and b such that lt(c) = lt(b)*res_a + res_r.   
 */
void divideByLeadingTerms_AAZ (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar);


/*****************
 * Derivative / Integral
 *****************/

/**
 * Get the k-th partial derivative of aa for the variable at index idx.
 *
 * returns the partial derivative
 */
AltArrZ_t* derivative_AAZ(const AltArrZ_t* aa, int idx, int k);

/**
 * Turn the polynomial pointed to by aa into its k'th partial derivative 
 * with respect to the varibale at index idx.
 *
 * @param[in,out] aa a pointer to the polynomial whose partial derivative is to be computed, also the place where the resulting derivative is stored.
 * @param idx the index of the variable with respect to which the partial derviative is taken.
 * @param k the order of the derivative. 
 *
 */
void derivative_AAZ_inp(AltArrZ_t** aa, int idx, int k);

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
void integralContent_AAZ(const AltArrZ_t* aa, mpz_t ret);

/**
 * Get the primitive part, a primitive integer polynomial, of aa.
 *
 * returns the primitive part.
 */
AltArrZ_t* primitivePart_AAZ(const AltArrZ_t* aa);

/**
 * Get the primitive part, a primitive integer polynomial, of aa.
 * Returns the integral content in cont.
 *
 * returns the primitive part.
 */
AltArrZ_t* primitivePartAndContent_AAZ(const AltArrZ_t* aa, mpz_t cont);

/**
 * Get the primitive part, a primitive integer polynomial, 
 * from a rational number polynomial aa.
 * Returns the content in cont.
 *
 * returns the primitive part.
 *
 */
AltArrZ_t* primitivePartAndContent_AAZFromAA(AltArr_t* aa, mpq_t cont);

/**
 * Convert the input aa to its primitive part in place.
 */
void primitivePart_AAZ_inp(AltArrZ_t* aa);

/**
 * Convert the input aa to its primitive part in place and
 * return its content in the mpz_t content.
 */
void primitivePartContent_AAZ_inp(AltArrZ_t* aa, mpz_t content);

/**
 * Given two polynomials, both with 1 variable, get their GCD using a heuristic method.
 * Returns NULL if the heuristic fails or if both inputs encode 0.
 *
 * returns the GCD or NULL
 */
AltArrZ_t* univariateHeuristicGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b);

/**
 * Given two polynomials, both with 1 variable and both primitive, 
 * get their GCD using a heuristic method.
 * Returns NULL if the heuristic fails or if both inputs encode 0.
 *
 * returns the GCD or NULL
 */
AltArrZ_t* univariateHeuristicPrimitiveGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b);

/**
 * Given two polynomials, both with 1 variable, get their GCD.
 *
 * returns the GCD.
 */
AltArrZ_t* univariateGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b);

/** 
 * Given two polynomials, both with 1 variable, get their GCD as a primitive polynomial.
 * Is s or t is non-null then returns the corresponding bezout coefficient as well such that
 * s*a + t*b = g.
 *
 * returns the gcd, g.
 */
// AltArrZ_t* extendedEuclidean(AltArrZ_t* a, AltArrZ_t* b, AltArrZ_t** s, AltArrZ_t** t);

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
AltArrZ_t* commonFactor_AAZ(const AltArrZ_t* a, AltArrZ_t** factored);



/*****************
 * Interpolation / Evaluation
 *****************/

/**
 * Evaluate a univariate polynomial at the point, point. The result is 
 * returned as a multi-precision rational number, res.
 */
void univarEvaluate_AAZ(AltArrZ_t* aa, const mpz_t point, mpz_t res);

 


/****************
* Normal Form (Multi-Divisor Division (MDD))
*****************/

/** 
 * Normal Form (or Multi-Divisor Division (MDD)) 
 * Given the dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}},
 * such that f = q_0*g_0 + ... + q_{s-1}*g_{s-1} + r.
 * Note the variable type is determining the type of normal form algorithm such that
 *  if type = 0 then using the naive normal form algorithm which works with general divisor-sets
 *  if type > 0 then using the normal form algorithm which is specialized for triangular-sets works recursively
 * Note the triangular-set normal form algorithm over recursive representation of polynoimials 
 *  is implemented in SMQP_Support_Recurisve-AA.h
 */
void multiDivisorDivision_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar, int type);

/**
 * The specific Normal Form algorithm to compute only the normal form of f w.r.t G
 */
AltArrZ_t* onlyNormalForm_AAZ (AltArrZ_t* f, AltArrZ_t** G, int s, int nvar); 

/** 
 * Multi-Divisor Division (MDD) using Heap (Johnson's) Division
 * Given the dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
void normalForm_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);
void heapMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);

/** 
 * Multi-Divisor Division (MDD) where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
void triangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);
void recTriangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);

/** 
 * Multi-Divisor Division (MDD) where the divisor-set is a Triangular Set
 * with using primitive factorization techniques 
 * Given the dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
void primitiveFactorTriangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar);

/**
 * given the number of variables, nvar, and return a recursive loop of the routine doing 
 * in triangular-set normal form algorithm  
 */
static inline int* recursiveLoopZ (int nvar)
{
  int size = (1 << nvar);
  int* seq = (int*) calloc( size-1, sizeof (int));

  if (nvar == 1)
      return seq;
  else {
	  int* seq_part = recursiveLoopZ (nvar-1);
	  for (int i = 0; i < (size/2)-1; i++) {
    	  seq[i] = seq_part[i];
    	  seq[(size/2) +i] = seq_part[i];
      }

      free(seq_part);
      seq[(size/2)-1] = nvar-1;
      return seq;
    }
}

/** 
 * Check the input triangular set is normalized or not, 
 * s is the size of the triangular set and nvar is the number of variables 
 */
int isNormalizedTriangularSet_AAZ (AltArrZ_t** G, int s, int nvar);

/** 
 * The verification algorithm to test the correctness of normal form algorithms.
 */
int multiDivisorDivisionVerification_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t* r, AltArrZ_t* hPow, int nSet, int tnvar);

#ifdef __cplusplus
}
#endif

#endif
