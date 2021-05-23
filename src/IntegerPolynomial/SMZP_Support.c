
#include "IntegerPolynomial/SMZP_Support.h"
#include "IntegerPolynomial/SMZP_Hensel.h"
#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "IntegerPolynomial/SMZP_Factoring_Support.h"
#include "IntegerPolynomial/DUZP_Support.h"
#include "Utils/MacroHelpers.h"
#include "Utils/C_To_Cpp.h"
#include "Parser/bpas_parser.h"

#define HEURISTIC_GCD_CUTOFF_INFN    60
#define HEURISTIC_GCD_CUTOFF_SIZE 25000
#define HEURISTIC_GCD_ATTEMPTS        5

/*****************
 * Alternating Array helpers
 *****************/
static void printAAZ(AltArrZ_t* aa) {
	if (aa == NULL) {
		fprintf(stderr, "0\n");
		return;
	}
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		gmp_fprintf(stderr, "%Zd*%llx + ", aa->elems[i].coef, aa->elems[i].degs);
	}
	fprintf(stderr, "\n");
}

void printPolyOptions_AAZ(FILE* fp, const AltArrZ_t* aa, const char** syms, int nvar, int positiveDegsOnly) {
	if (isZero_AAZ(aa)) {
		fprintf(fp, "0");
		return;
	}
	if (isConstant_AAZ(aa)) {
		gmp_fprintf(fp, "%Zd", aa->elems[0].coef);
		return;
	}

	if (aa->unpacked) {
		printPoly_AAZ_unpk(fp, aa, syms, nvar);
		return;
	}

	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	gmp_fprintf(fp, "%Zd", aa->elems[0].coef);
	if (!isZeroExponentVector(aa->elems[0].degs)) {
		fprintf(fp, "*");
		if (positiveDegsOnly) {
			printPositiveDegs_AA(fp, aa->elems[0].degs, syms, nvar, masks, sizes);
		} else {
			printDegs_AA(fp, aa->elems[0].degs, syms, nvar, masks, sizes);
		}
	}
	for (int i = 1; i < AA_SIZE(aa)-1; ++i) {
		if (mpz_sgn(aa->elems[i].coef) > 0) {
			gmp_fprintf(fp, " + %Zd", aa->elems[i].coef);
		} else {
			mpz_neg(aa->elems[i].coef, aa->elems[i].coef);
			gmp_fprintf(fp, " - %Zd", aa->elems[i].coef);
			mpz_neg(aa->elems[i].coef, aa->elems[i].coef);
		}
		if (!isZeroExponentVector(aa->elems[i].degs)) {
			fprintf(fp, "*");
			if (positiveDegsOnly) {
				printPositiveDegs_AA(fp, aa->elems[i].degs, syms, nvar, masks, sizes);
			} else {
				printDegs_AA(fp, aa->elems[i].degs, syms, nvar, masks, sizes);
			}
		}
	}
	if (aa->size > 1) {
		if (mpz_sgn(aa->elems[aa->size-1].coef) > 0) {
			gmp_fprintf(fp, " + %Zd", aa->elems[aa->size-1].coef);
		} else {
			mpz_neg(aa->elems[aa->size-1].coef, aa->elems[aa->size-1].coef);
			gmp_fprintf(fp, " - %Zd", aa->elems[aa->size-1].coef);
			mpz_neg(aa->elems[aa->size-1].coef, aa->elems[aa->size-1].coef);
		}
		if (!isZeroExponentVector(aa->elems[aa->size-1].degs)) {
			fprintf(fp, "*");
			if (positiveDegsOnly) {
				printPositiveDegs_AA(fp, aa->elems[aa->size-1].degs, syms, nvar, masks, sizes);
			} else {
				printDegs_AA(fp, aa->elems[aa->size-1].degs, syms, nvar, masks, sizes);
			}
		}
	}

	//fprintf(fp, "\n");
}


void calculateMaxDegsList_AAZ(const AltArrZ_t* aa, degree_t* maxList) {
	if (maxList == NULL) {
		return;
	}

	if (isZero_AAZ(aa)) {
		if (aa != NULL) {
			int nvar = aa->nvar;
			for (int i = 0; i < nvar; ++i) {
				maxList[i] = 0;
			}
		}
		return;
	}


	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	for (int i = 0; i < nvar; ++i) {
		maxList[i] = 0;
	}
	degree_t deg;
	for (register int i = 0; i < size; ++i) {
		for (register int j = 0; j < nvar; ++j) {
			deg = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
			maxList[j] =  deg > (maxList[j]) ? deg : maxList[j];
		}
	}
}

degrees_t calculateMaxDegs_AAZ(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0 || aa->nvar == 0) {
		return 0;
	}

	register unsigned int size = aa->size;
	register unsigned short nvar = aa->nvar;
	degrees_t maxList[nvar];
	AAZElem_t* elems = aa->elems;
	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	register unsigned int i;
	register unsigned int j;
	for (j = 0; j < nvar; ++j) {
		maxList[j] = 0;
	}
	for (i = 0; i < size; ++i) {
		for (j = 0; j < nvar; ++j) {
			maxList[j] = (elems[i].degs & masks[j]) > (maxList[j]) ? (elems[i].degs & masks[j]) : maxList[j];
		}
	}

	degrees_t max = 0ll;
	for (register int j = 0; j < nvar; ++j) {
		max |= maxList[j];
	}

	return max;
}

int isInOrder_AAZ(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size < 2) {
		return 1;
	}
	for (int i = 1; i < aa->size; ++i) {
		if (!isGreaterExponentVectors(aa->elems[i-1].degs, aa->elems[i].degs)) {
			fprintf(stderr, "\n\nIndex: %d is not in order!\n", i);
			return 0;
		}
	}
	return 1;
}

int isExactlyEqual_AAZ(const AltArrZ_t* a, const AltArrZ_t* b) {
	if (a == NULL) {
		return (b == NULL);
	}
	if (b == NULL) {
		return 0;
	}

	if (a->size != b->size) {
		return 0;
	}

	if (a->unpacked || b->unpacked) {
		return isExactlyEqual_AAZ_unpk(a, b);
	}

	for (int i = 0; i < a->size; ++i) {
		if (a->elems[i].degs != b->elems[i].degs) {
			return 0;
		}
		if (mpz_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
			return 0;
		}
	}

	return 1;
}

void nonZeroVariables_AAZ(const AltArrZ_t* aa, int* foundVar) {
	if (aa == NULL || aa->nvar == 0 || isConstant_AAZ(aa)) {
		return;
	}

	if (aa->unpacked) {
		nonZeroVariables_AAZ_unpk(aa, foundVar);
		return;
	}

	int nvar = aa->nvar;
	degrees_t degs = 0;
	const degrees_t* masks = getExpMaskArray(nvar);

	for (int i = 0; i < nvar; ++i) {
		foundVar[i] = 0;
	}
	foundVar[0] = ( (aa->elems[0].degs & masks[0]) > 0);
	int searchIdx = 1;
	for (int i = 0; i < aa->size && searchIdx < nvar; ++i) {
		degs = aa->elems[i].degs;
		if ( (degs & masks[searchIdx]) > 0) {
			foundVar[searchIdx] = 1;
			while (searchIdx < nvar && foundVar[searchIdx] == 1) {
				++searchIdx;
			}
		}

		for (int j = searchIdx; j < nvar; ++j) {
			if ((degs & masks[j]) > 0) {
				foundVar[j] = 1;
			}
		}
	}
}

degree_t totalDegree_AAZ(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0 || isConstant_AAZ(aa)) {
		return 0;
	}

	if (aa->unpacked) {
		return totalDegree_AAZ_unpk(aa);
	}

	int nvar = aa->nvar;
	int size = aa->size;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	degree_t total = 0, totalMax = 0;
	degrees_t degs;
	for (int i = 0; i < size; ++i) {
		total = 0;
		degs = aa->elems[i].degs;
		for (int j = 0; j < nvar; ++j) {
			total += GET_NTH_EXP(degs, masks[j], sizes[j]);
		}
		if (total > totalMax) {
			totalMax = total;
		}
	}

	return totalMax;
}

degree_t partialDegree_AAZ(const AltArrZ_t* aa, int k) {
	if (aa == NULL || aa->size == 0 || k >= aa->nvar) {
		return -1;
	}

	if (aa->unpacked) {
		return partialDegree_AAZ_unpk(aa, k);
	}

	int nvar = aa->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	degree_t dMax = 0;
	if (k == 0) {
		dMax = GET_NTH_EXP(aa->elems->degs, masks[k], sizes[k]);
	} else {
		degree_t d = 0;
		for (int i = 0; i < aa->size; ++i) {
			d = GET_NTH_EXP(aa->elems[i].degs, masks[k], sizes[k]);
			if (d > dMax) {
				dMax = d;
			}
		}
	}

	return dMax;
}

void partialDegrees_AAZ(const AltArrZ_t* aa, degree_t* degsList) {
	if (aa == NULL || aa->size == 0 || aa->nvar == 0 || degsList == NULL) {
		return;
	}

	if (aa->unpacked) {
		partialDegrees_AAZ_unpk(aa, degsList);
		return;
	}

	calculateMaxDegsList_AAZ(aa, degsList);
}

void setPartialDegreeTerm_AAZ(AltArrZ_t* aa, int idx, int k, degree_t deg) {
	if (isZero_AAZ(aa)) {
		return;
	}
	if (idx >= aa->size || k >= aa->nvar) {
		return;
	}
	if (aa->unpacked) {
		((degree_t*) aa->elems[idx].degs)[k] = deg;
	} else {
		const int* sizes = getExpOffsetArray(aa->nvar);
		const degrees_t* maxDegs = getMaxExpArray(aa->nvar);
		const degrees_t* masks = getExpMaskArray(aa->nvar);
		if (deg > maxDegs[k]) {
			unpackExponentVectors_AAZ_inp(aa);
			setPartialDegreeTerm_AAZ(aa, idx, k, deg);
			return;
		}
		aa->elems[idx].degs &= ~(masks[k]);
		aa->elems[idx].degs |= (((degrees_t) deg) << sizes[k]);
	}
	// mergeSortPolynomial_AA(aa);
}

void lowDegrees_AAZ(const AltArrZ_t* aa, degree_t* degs) {
	if (degs == NULL) {
		return;
	}

	if (isZero_AAZ(aa)) {
		if (aa != NULL) {
			int nvar = aa->nvar;
			for (int i = 0; i < nvar; ++i) {
				degs[i] = 0;
			}
		}
		return;
	}

	if (aa->unpacked) {
		lowDegrees_AAZ_unpk(aa, degs);
		return;
	}


	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	degree_t allZero = 0;
	for (int i = 0; i < nvar; ++i) {
		degs[i] = GET_NTH_EXP(elems[size-1].degs, masks[i], sizes[i]);
		allZero |= degs[i];
	}

	degree_t deg;
	//note the && allZero > 0 for fast bail-out;
	for (register int i = size-2; i >= 0 && allZero > 0; --i) {
		allZero = 0;
		for (register int j = 0; j < nvar; ++j) {
			deg = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
			degs[j] =  deg < (degs[j]) ? deg : degs[j];
			allZero |= degs[j];
		}
	}
}

void removeLowDegrees_AAZ_inp(AltArrZ_t* aa, degree_t* degs) {
	if (isZero_AAZ(aa) || aa->nvar == 0) {
		if (aa != NULL && degs != NULL) {
			int nvar = aa->nvar;
			for (int i = 0; i < nvar; ++i) {
				degs[i] = 0;
			}
		}
		return;
	}

	if (aa->unpacked) {
		removeLowDegrees_AAZ_inp_unpk(aa, degs);
		return;
	}

	degree_t stackDegs[aa->nvar];
	if (degs == NULL) {
		//Note that doing this will not impact calling scope.
		degs = stackDegs;
	}

	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	degrees_t packedLDegs = elems[size-1].degs;
	for (int i = 0; i < nvar; ++i) {
		degs[i] = GET_NTH_EXP(elems[size-1].degs, masks[i], sizes[i]);
	}

	degree_t deg;
	//note the && packedLDegs > 0 for fast bail-out;
	for (register int i = size-2; i >= 0 && packedLDegs > 0; --i) {
		packedLDegs = 0ll;
		for (register int j = 0; j < nvar; ++j) {
			deg = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
			degs[j] =  deg < (degs[j]) ? deg : degs[j];
			packedLDegs |= (((degrees_t) degs[j]) << sizes[j]);
		}
	}

	if (packedLDegs > 0) {
		for (register int i = 0; i < size; ++i) {
			subtractExponentVectors(elems[i].degs, packedLDegs, elems[i].degs);
		}
	}
}



degree_t mainDegree_AAZ(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return -1;
	}

	if (isConstant_AAZ(aa)) {
		return 0;
	}

	if (aa->unpacked) {
		return mainDegree_AAZ_unpk(aa);
	}

	int nvar = aa->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	int ret = 0;
	for (int j = 0; j < nvar; ++j) {
		if ((aa->elems->degs & masks[j]) != 0) {
			ret = GET_NTH_EXP(aa->elems->degs, masks[j], sizes[j]);
			break;
		}
	}
	return ret;
}

int mainVariable_AAZ(AltArrZ_t* aa) {
	if (aa == NULL || aa->nvar == 0) {
		return -1;
	}

	if (aa->unpacked) {
		return mainVariable_AAZ_unpk(aa);
	}

	int lvIdx = -1;
	int nvar = aa->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
	for (int i = 0; i < nvar; ++i) {
		if ((aa->elems->degs & masks[i]) != 0) {
			lvIdx = i;
			break;
		}
	}

	return lvIdx;
}

void coefficient_AAZ(AltArrZ_t* aa, const degree_t* degsList, int nvar, mpz_t retCoef) {
	if (isZero_AAZ(aa)) {
		mpz_set_si(retCoef, 0l);
		return;
	}
	if (nvar == 0 && aa->nvar == 0) {
		mpz_set(retCoef, aa->elems->coef);
		return;
	}

	if (aa->unpacked) {
		return coefficient_AAZ_unpk(aa, degsList, nvar, retCoef);
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	degrees_t degs = 0;
	for (int j = 0; j < nvar; ++j) {
		if (degsList[j] > maxExps[j]) {
			//if we are not unpacked and we are requesting a coefficent for a monomial
			//which cannot be packed then obviously this coefficient does not exist;
			mpz_set_ui(retCoef, 0l);
			return;
		}
		degs |= ((degrees_t)degsList[j] << sizes[j]);
	}

	mpz_set_ui(retCoef, 0l);
	for (int i = 0; i < aa->size; ++i) {
		if (isEqualExponentVectors(degs, aa->elems[i].degs)) {
			mpz_set(retCoef, aa->elems[i].coef);
			break;
		}
		if (isLessExponentVectors(aa->elems[i].degs, degs)) {
			//due to ordering of terms, none will have degs as monomial.
			break;
		}
	}
}

void setCoefficient_AAZ(AltArrZ_t* aa, const degree_t* degsList, int nvar, const mpz_t coef) {
	if (aa == NULL || nvar != aa->nvar) {
		return;
	}

	if (aa->size == 0 || nvar == 0) {
		if (aa->alloc < 1) {
			resizePolynomial_AAZ(aa, 1);
		}
		mpz_set(aa->elems->coef, coef);
		if (nvar != 0) {
			setDegrees_AAZ_inp (aa, 0, degsList, nvar);
		}
		return;
	}

	if (aa->unpacked) {
		setCoefficient_AAZ_unpk(aa, degsList, nvar, coef);
		return;
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	degrees_t degs = 0;
	for (int j = 0; j < nvar; ++j) {
		if (degsList[j] > maxExps[j]) {
			unpackExponentVectors_AAZ_inp(aa);
			setCoefficient_AAZ_unpk(aa, degsList, nvar, coef);
			return;
		}
		degs |= ((degrees_t) degsList[j] << sizes[j]);
	}

	for (int i = 0; i < aa->size; ++i) {
		if (isEqualExponentVectors(aa->elems[i].degs, degs)) {
			mpz_set(aa->elems[i].coef, coef);
			return;
		}
		if (isLessExponentVectors(aa->elems[i].degs, degs)) {
			//shift everything to the right 1 and insert at i.
			if (aa->size >= aa->alloc) {
				resizePolynomial_AAZ(aa, aa->size + 1);
			}
			mpz_init(aa->elems[aa->size].coef);
			for (int j = aa->size; j > i; --j) {
				mpz_swap(aa->elems[j].coef, aa->elems[j-1].coef);
				aa->elems[j].degs = aa->elems[j-1].degs;
			}
			aa->elems[i].degs = degs;
			mpz_set(aa->elems[i].coef, coef);
			++(aa->size);
			return;
		}
	}

	//if we get here then we need to insert at the end of the array;

	if (aa->size >= aa->alloc) {
		resizePolynomial_AAZ(aa, aa->size + 1);
	}
	mpz_init(aa->elems[aa->size].coef);
	mpz_set(aa->elems[aa->size].coef, coef);
	aa->elems[aa->size].degs = degs;
	++(aa->size);
}

int isEqualWithVariableOrdering_AAZ(AltArrZ_t* a, AltArrZ_t*b, const int* xs, int xsSize) {
	if (isZero_AAZ(a)) {
		if (isZero_AAZ(b)) {
			return 1;
		} else {
			return 0;
		}
	}
	if (isZero_AAZ(b)) {
		return 0;
	}

	if (a->size != b->size) {
		return 0;
	}

	if (a->unpacked || b->unpacked) {
		return isEqualWithVariableOrdering_AAZ_unpk(a, b, xs, xsSize);
	}

	int nvar = a->nvar;
	int bnvar = b->nvar;
	int v = xsSize / 2;

	const degrees_t* aMasks = getExpMaskArray(nvar);
	const degrees_t* bMasks = getExpMaskArray(bnvar);

	const int* aSizes = getExpOffsetArray(nvar);
	const int* bSizes = getExpOffsetArray(bnvar);

	int ret = 1;
	degree_t adeg, bdeg;
	for (int i = 0; i < a->size; ++ i){
		if (mpz_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
			ret = 0;
			break;
		}
		for (int j = 0; j < v; ++j) {
			adeg = GET_NTH_EXP(a->elems[i].degs, aMasks[xs[2*j]-1], aSizes[xs[2*j]-1]);
			bdeg = GET_NTH_EXP(b->elems[i].degs, bMasks[xs[2*j+1]-1], bSizes[xs[2*j+1]-1]);
			if (xs[2*j] && xs[2*j+1] && (adeg != bdeg)) {
				ret = 0;
				break;
			} else if (!xs[2*j] && xs[2*j+1] && (bdeg != 0) ) {
				ret = 0;
				break;
			} else if (xs[2*j] && !xs[2*j+1] && (adeg != 0) ) {
				ret = 0;
				break;
			}
		}

	}

	return ret;
}

AltArrZ_t* termAtIdx_AAZ(AltArrZ_t* a, int idx) {
	if (a == NULL || idx >= a->size) {
		return NULL;
	}

	AltArrZ_t localA = *a;
	localA.elems = localA.elems + idx;
	localA.size = 1;
	return deepCopyPolynomial_AAZ(&localA);
}

AltArrZ_t* homogeneousPart_AAZ(AltArrZ_t* a, int tdeg) {
	if (isZero_AAZ(a)) {
		return NULL;
	}

	int curAlloc = 10;
	int curIdx = 0;
	int nvar = a->nvar;
	AltArrZ_t* ret = makePolynomial_AAZ(curAlloc, nvar);
	int size = a->size;

	//try to do this in a way which is agnostic to packed or unpacked.
	degree_t degs[a->nvar];
	int k;
	degree_t termtdeg;
	for (int i = 0; i < size; ++i) {
		partialDegreesTerm_AAZ(a, i, degs);
		termtdeg = 0;
		for (k = 0; k < nvar; ++k) {
			termtdeg += degs[k];
		}
		if (termtdeg == tdeg) {
			if (curIdx >= curAlloc) {
				curAlloc *= 2;
				resizePolynomial_AAZ(ret, curAlloc);
			}
			setDegrees_AAZ_inp(ret, curIdx, degs, nvar);
			mpz_init(ret->elems[curIdx].coef);
			mpz_set(ret->elems[curIdx].coef, a->elems[i].coef);
			++curIdx;
		}
	}


	if (curIdx == 0) {
		freePolynomial_AAZ(ret);
		return NULL;
	}

	ret->size = curIdx;
	mergeSortPolynomial_AAZ(ret);
	return ret;
}

int isHomogeneous_AAZ(const AltArrZ_t* a) {
	if (isZero_AAZ(a) || isConstant_AAZ(a)) {
		return 1;
	}

	int size = a->size;
	degree_t degs[a->nvar];
	int k;
	int nvar = a->nvar;
	degree_t tdeg = 0;
	partialDegreesTerm_AAZ(a, 0, degs);
	for (k = 0; k < nvar; ++k) {
		tdeg += degs[k];
	}

	degree_t termtdeg;
	for (int i = 1; i < size; ++i) {
		partialDegreesTerm_AAZ(a, i, degs);
		termtdeg = 0;
		for (k = 0; k < nvar; ++k) {
			termtdeg += degs[k];
		}
		if (termtdeg != tdeg) {
			return 0;
		}
	}
	return 1;

}

int isConstantTermZero_AAZ(AltArrZ_t* a) {
	if (a == NULL || a->size == 0) {
		return 1;
	}

	if (a->unpacked) {
		return !isZeroExponentVector_unpk(a->elems[a->size-1].degs, a->nvar);
	} else {
		return !isZeroExponentVector(a->elems[a->size-1].degs);
	}
}

void expandNumVars_AAZ(AltArrZ_t* aa_in, int newNvar) {
	if (aa_in == NULL) {
		return;
	}
	if (aa_in->size == 0 || aa_in->alloc == 0) {
		aa_in->nvar = newNvar;
		return;
	}

	if (newNvar <= aa_in->nvar) {
		return;
	}

	if (aa_in->unpacked) {
		expandNumVars_AAZ_unpk(aa_in, newNvar);
		return;
	}

	if (aa_in->nvar == 0) {
		aa_in->nvar = newNvar;
		return;
	}

	AltArrZ_t* aa = deepCopyPolynomial_AAZ(aa_in);

	const degrees_t* __restrict__ oldMasks = getExpMaskArray(aa->nvar);

	const int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	const int* __restrict__ newSizes = getExpOffsetArray(newNvar);

	const degrees_t* maxExps = getMaxExpArray(newNvar);

	degrees_t degs;
	degrees_t curDeg;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < aa->nvar; ++j) {
			curDeg = GET_NTH_EXP(degs, oldMasks[j], oldSizes[j]);
			if (curDeg > maxExps[j]) {
				freePolynomial_AAZ(aa);
				unpackExponentVectors_AAZ_inp(aa_in);
				expandNumVars_AAZ_unpk(aa_in, newNvar);
				return;

				// int diff = newNvar - aa->nvar;
				// fprintf(stderr, "SMZP ERROR: Overflow in exponent packing for expand at index %d; %lld > %lld.\n", j+diff, curDeg, maxExps[j+diff]);
				// exit(1);
			}
			aa->elems[i].degs |= (curDeg << newSizes[j]);
		}
	}

	aa_in->nvar = newNvar;
	AAZElem_t* tempElems = aa_in->elems;
	aa_in->elems = aa->elems;
	aa->elems = tempElems;

	freePolynomial_AAZ(aa);
}

void expandNumVarsLeft_AAZ(AltArrZ_t* aa_in, int newNvar) {

	if (aa_in == NULL) {
		return;
	}
	if (newNvar <= aa_in->nvar) {
		return;
	}
	if (aa_in->size == 0 || aa_in->alloc == 0) {
		aa_in->nvar = newNvar;
		return;
	}

	if (aa_in->unpacked) {
		expandNumVarsLeft_AAZ_unpk(aa_in, newNvar);
		return;
	}

	if (aa_in->nvar == 0) {
		aa_in->nvar = newNvar;
		return;
	}

	AltArrZ_t* aa = deepCopyPolynomial_AAZ(aa_in);

	const degrees_t* __restrict__ oldMasks = getExpMaskArray(aa->nvar);
	const int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	const int* __restrict__ newSizes = getExpOffsetArray(newNvar);
	const degrees_t* maxExps = getMaxExpArray(newNvar);

	degrees_t degs;
	degrees_t curDeg;
	int diff = newNvar - aa->nvar;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < aa->nvar; ++j) {
			curDeg = GET_NTH_EXP(degs, oldMasks[j], oldSizes[j]);
			if (curDeg > maxExps[j+diff]) {
				freePolynomial_AAZ(aa);
				unpackExponentVectors_AAZ_inp(aa_in);
				expandNumVarsLeft_AAZ_unpk(aa_in, newNvar);
				return;
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for expand at index %d; %lld > %lld.\n", j+diff, curDeg, maxExps[j+diff]);
				// exit(1);
			}
			aa->elems[i].degs |= (curDeg << newSizes[j+diff]);
		}
	}

	aa_in->nvar = newNvar;
	AAZElem_t* tempElems = aa_in->elems;
	aa_in->elems = aa->elems;
	aa->elems = tempElems;

	freePolynomial_AAZ(aa);
}

void shrinkNumVarsAtIdx_AAZ(AltArrZ_t* aa, int idx) {
	if (aa == NULL || aa->nvar < 1) {
		return;
	}
	if (aa->size < 1) {
		--(aa->nvar);
		return;
	}
	if (aa->nvar == 1) {
		for (int i = 0; i < aa->size; ++i) {
			aa->elems[i].degs = 0;
		}
		aa->nvar = 0;
		return;
	}

	if (aa->unpacked) {
		shrinkNumVarsAtIdx_AAZ_unpk(aa, idx);
		return;
	}

	int nvar = aa->nvar;
	int newNvar = aa->nvar - 1;

	AAZElem_t* elems = aa->elems;

	const degrees_t* __restrict__ oldMasks = getExpMaskArray(aa->nvar);
	const int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	const int* __restrict__ newSizes = getExpOffsetArray(newNvar);

	degrees_t newDegs;
	degrees_t curDegs;
	degrees_t deg;
	int j;
	for (int i = 0; i < aa->size; ++i ) {
		curDegs = elems[i].degs;
		newDegs = 0;
		//iterate this current nvar, skipping j = idx, repacking exponents.
		//when j > idx, exponents are shifted to be at an index one less than
		//originally.
		for (j = 0; j < idx; ++j) {
			deg = GET_NTH_EXP(curDegs, oldMasks[j], oldSizes[j]);
			newDegs |= (deg << newSizes[j]);
		}
		for (j = idx+1; j < nvar; ++j) {
			deg = GET_NTH_EXP(curDegs, oldMasks[j], oldSizes[j]);
			newDegs |= (deg << newSizes[j-1]);
		}

		elems[i].degs = newDegs;
	}

	aa->nvar = newNvar;
}

void shrinkAndReorderVars_AAZ(AltArrZ_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->nvar < 1) {
		return;
	}
	if (varmapSize > aa->nvar) {
		return;
	}

	if (aa->unpacked) {
		shrinkAndReorderVars_AAZ_unpk(aa, varMap, varmapSize);
		return;
	}

	int newNvar = 0;
	int needSort = 0;
	int maxSoFar = -1;
	for (int i = 0; i < varmapSize; ++i) {
		if (varMap[i] >= 0) {
			++newNvar;

			//the below checks for reordering, setting needSort if reorder occurs.
			if (varMap[i] < maxSoFar) {
				needSort = 1;
			} else {
				maxSoFar = varMap[i];
			}
		}
	}

	if (aa->size < 1) {
		aa->nvar = newNvar;
		return;
	}

	if (newNvar == 0) {
		for (int i = 1; i < aa->size; ++i) {
			mpz_clear(aa->elems[i].coef);
		}
		aa->elems->degs = 0ll;
		aa->size = 1;
		aa->nvar = 0;
		return;
	}

	const degrees_t* oldMasks = getExpMaskArray(aa->nvar);
	const int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	const int* __restrict__ newSizes = getExpOffsetArray(newNvar);
	const degrees_t* __restrict__ maxExps = getMaxExpArray(newNvar);

	AltArrZ_t* tmpAA = deepCopyPolynomial_AAZ(aa);
	AAZElem_t* elems = tmpAA->elems;
	degrees_t newDegs, oldDegs, curDeg;
	int j;
	for (int i = 0; i < aa->size; ++i) {
		oldDegs = elems[i].degs;
		newDegs = 0;
		for (j = 0; j < varmapSize; ++j) {
			if (varMap[j] < 0) {
				continue;
			}
			curDeg = GET_NTH_EXP(oldDegs, oldMasks[j], oldSizes[j]);
			if (curDeg > maxExps[varMap[j]]) {
				freePolynomial_AAZ(tmpAA);
				unpackExponentVectors_AAZ_inp(aa);
				shrinkAndReorderVars_AAZ_unpk(aa, varMap, varmapSize);
				return;
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for shrink and reorder vars. At new index %d.", varMap[j]);
				// exit(1);
			}
			newDegs |= (curDeg << newSizes[varMap[j]]);
		}
		elems[i].degs = newDegs;
	}

	AAZElem_t* tmpElems = aa->elems;
	aa->elems = elems;
	tmpAA->elems = tmpElems;
	freePolynomial_AAZ(tmpAA);

	aa->nvar = newNvar;

	if (needSort) {
		mergeSortPolynomial_AAZ(aa);
	}
}

//     degs[varmap[i]] = curDegs[i];
void reorderVars_AAZ(AltArrZ_t* aa_in, int* varMap, int varmapSize) {
	if (aa_in == NULL || aa_in->size == 0) {
		return;
	}

	if (aa_in->unpacked) {
		reorderVars_AAZ_unpk(aa_in, varMap, varmapSize);
		return;
	}

	int nvar = aa_in->nvar;

	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	const int* __restrict__ sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	AltArrZ_t* aa = deepCopyPolynomial_AAZ(aa_in);

	degrees_t degs;
	degrees_t curDeg;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < varmapSize; ++j) {
			curDeg = GET_NTH_EXP(degs, masks[j], sizes[j]);
			if (curDeg > maxExps[varMap[j]]) {
				freePolynomial_AAZ(aa);
				unpackExponentVectors_AAZ_inp(aa_in);
				reorderVars_AAZ_unpk(aa_in, varMap, varmapSize);
				return;
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for reorder at index %d; %lld > %lld.\n", varMap[j], curDeg, maxExps[varMap[j]]);
				// exit(1);
			}
			aa->elems[i].degs |= (curDeg << sizes[varMap[j]]);
		}
	}

	AAZElem_t* tmpElems = aa_in->elems;
	aa_in->elems = aa->elems;
	aa->elems = tmpElems;
	freePolynomial_AAZ(aa);

	mergeSortPolynomial_AAZ(aa_in);
}

int tryShrinkVariables_AAZ_inp(AltArrZ_t* aa, int* varMap) {
	int trueNvar = computeShrinkMap_AAZ(aa, varMap);
	if (trueNvar != aa->nvar) {
		shrinkAndReorderVars_AAZ(aa, varMap, aa->nvar);
	}
	return trueNvar;
}

int reverseShrinkVariables_AAZ_inp(AltArrZ_t* aa, int varmapSize, int* varMap) {
	if (aa == NULL || varMap == NULL || aa->nvar == varmapSize) {
		return aa == NULL ? 0 : aa->nvar;
	}

	int revMap[varmapSize];
	reverseShrinkMap_AAZ(varmapSize, aa->nvar, varMap, revMap);

	expandNumVars_AAZ(aa, varmapSize);
	reorderVars_AAZ(aa, revMap, varmapSize);
	return varmapSize;
}

void setDegrees_AAZ_inp(AltArrZ_t* aa, int idx, const degree_t* degsList, int nvar) {
	if (aa == NULL || idx >= aa->alloc) {
		return;
	}

	if (aa->unpacked) {
		setDegrees_AAZ_inp_unpk(aa, idx, degsList, nvar);
		return;
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxDegs = getMaxExpArray(nvar);

	aa->elems[idx].degs = 0;

	for (int k = 0; k < nvar; ++k) {
		if (degsList[k] > maxDegs[k]) {
			unpackExponentVectors_AAZ_inp(aa);
			setDegrees_AAZ_inp_unpk(aa, idx, degsList, nvar);
			break;

			// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for setDegrees_AA at index %d; %d > %lld.\n", k, degsList[k], maxDegs[k]);
			// exit(1);
		}

		aa->elems[idx].degs |= ((degrees_t) degsList[k] << sizes[k]);
	}
}

void setExponentTerm_AAZ_inp(AltArrZ_t* aa, int idx, degree_t deg, int k) {
	if (aa == NULL || idx >= aa->alloc || k >= aa->nvar) {
		return;
	}

	if (aa->unpacked) {
		setExponentTerm_AAZ_inp_unpk(aa, idx, deg, k);
		return;
	}

	int nvar = aa->nvar;
	const degrees_t* maxDegs = getMaxExpArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	if (deg > maxDegs[k]) {
		unpackExponentVectors_AAZ_inp(aa);
		setExponentTerm_AAZ_inp_unpk(aa, idx, deg, nvar);
	} else {
		aa->elems[idx].degs = aa->elems[idx].degs & ~(masks[k]);
		aa->elems[idx].degs |= ((degrees_t) deg << sizes[k]);
	}
}


AltArrZ_t* deepCopyPolynomial_AAZFromNode(Node* a, int nvar) {
	if (a == NULL) {
		return NULL;
	}
	polysize_t asize = numberOfTermsNode(a);
	AltArrZ_t* newAA = makePolynomial_AAZ(asize, nvar);
	AA_SIZE(newAA) = asize;

	int size = AA_SIZE(newAA);
	AAZElem_t* elems = newAA->elems;
	for (int i = 0; i < size && a != NULL; ++i) {
		mpz_init(elems[i].coef);
		mpz_set(elems[i].coef, mpq_numref(a->coef));
		elems[i].degs = a->degs;
		a = a->next;
	}

	return newAA;
}

Node* deepCopyPolynomial_NodeFromAAZ(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		fprintf(stderr, "Cannot make Node from unpacked AltArr_t\n");
		exit(1);
	}

	int asize = AA_SIZE(aa);

	AAZElem_t* elems = aa->elems;
	mpq_t qCoef;
	mpq_init(qCoef);
	mpz_set(mpq_numref(qCoef), elems[0].coef);
	Node* head = addTerm(NULL, elems[0].degs, qCoef);
	Node* tail = head;
	for (int i = 1; i < asize; ++i) {
		mpz_set(mpq_numref(qCoef), elems[i].coef);
		tail = addTerm(tail, elems[i].degs, qCoef);
	}
	mpq_clear(qCoef);
	return head;
}

AltArr_t* deepCopyPolynomial_AAFromAAZ(AltArrZ_t* aa) {
	if (isZero_AAZ(aa)) {
		return NULL;
	}

	if (aa->unpacked) {
		return deepCopyPolynomial_AAFromAAZ_unpk(aa);
	}

	AltArr_t* newAA = makePolynomial_AA(aa->size, aa->nvar);
	AAElem_t* elems = newAA->elems;
	AAZElem_t* oldelems = aa->elems;
	int size = AA_SIZE(newAA) = AA_SIZE(aa);
	for (int i = 0; i < size; ++i) {
		mpq_init(elems[i].coef);
		mpz_set(mpq_numref(elems[i].coef), oldelems[i].coef);
		elems[i].degs = oldelems[i].degs;
	}


	//TODO temporaily for testing unpacked SMQP;
	// unpackExponentVectors_AA_inp(newAA);
	// if (newAA->nvar > 0 && newAA->elems->degs == 0) {
	// 	fprintf(stderr, "from AAZ degs is NULL!\n" );
	// 	free((int*)-1);
	// }

	return newAA;
}

AltArrZ_t* deepCopyPolynomial_AAZFromAA(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}


	if (aa->unpacked) {
		return deepCopyPolynomial_AAZFromAA_unpk(aa);
	}

	int size = aa->size;
	int nvar = aa->nvar;
	AltArrZ_t* newAA = makePolynomial_AAZ(aa->alloc, nvar);
	AA_SIZE(newAA) = AA_SIZE(aa);
	AAZElem_t* elems = newAA->elems;
	AAElem_t* oldelems = aa->elems;
	for (int i = 0; i < size; ++i) {
		if (mpz_cmp_si(mpq_denref(oldelems[i].coef), 1l) != 0) {
			fprintf(stderr, "SMZP ERROR: Failed to convert a rational number polynomial to an integer one.\n");
			exit(1);
		}
		mpz_init(elems[i].coef);
		mpz_set(elems[i].coef, mpq_numref(oldelems[i].coef));
		elems[i].degs = oldelems[i].degs;
	}

	return newAA;

}

AltArrZDegList_t* deepCopyPolynomial_AAZDegListFromAA(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return deepCopyPolynomial_AAZDegListFromAA_unpk(aa);
	}

	int nvar = aa->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	AltArrZDegList_t* ret = makePolynomial_AAZDL(aa->alloc, nvar);
	AAZElem_DegList_t* rElems = ret->elems;
	AAZElem_t* elems = aa->elems;

	for (int i = 0; i < aa->size; ++i) {
		mpz_init(rElems[i].coef);
		mpz_set(rElems[i].coef, elems[i].coef);
		rElems[i].degs = malloc(sizeof(degree_t)*nvar);
		for (int j = 0; j < nvar; ++j) {
			rElems[i].degs[j] = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
		}
	}

	ret->size = aa->size;

	return ret;
}

AltArrZ_t* deepCopyPolynomial_AAZ(const AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return deepCopyPolynomial_AAZ_unpk(aa);
	}

	AltArrZ_t* newAA = makePolynomial_AAZ(aa->alloc, aa->nvar);
	int size = AA_SIZE(newAA) = AA_SIZE(aa);
	AAZElem_t* elems = newAA->elems;
	AAZElem_t* oldelems = aa->elems;
	for (int i = 0; i < size; ++i) {
		mpz_init(elems[i].coef);
		mpz_set(elems[i].coef, oldelems[i].coef);
		elems[i].degs = oldelems[i].degs;
	}

	return newAA;
}


void deepCopyPolynomial_AAZ_inp(const AltArrZ_t* aa, AltArrZ_t** bb){
	if (bb == NULL) {
		return;
	}

	if (isZero_AAZ(aa)) {
		if (!isZero_AAZ(*bb)) {
			freePolynomial_AAZ(*bb);
		}
		return;
	}

	if (aa->unpacked || (*bb != NULL && (*bb)->unpacked)) {
		deepCopyPolynomial_AAZ_inp_unpk(aa, bb);
		return;
	}

	AltArrZ_t* b = *bb;
	if (b == NULL) {
		b = deepCopyPolynomial_AAZ(aa);
		*bb = b;
		return;
	} else if (b->alloc < aa->size) {
		resizePolynomial_AAZ(b, aa->size);
	}

	polysize_t i;
	if (aa->size < b->size) {
		for (i = 0; i < aa->size; ++i) {
			mpz_set(b->elems[i].coef, aa->elems[i].coef);
			b->elems[i].degs = aa->elems[i].degs;
		}
		for ( ; i < b->size; ++i) {
			mpz_clear(b->elems[i].coef);
		}
	} else {
		for (i = 0; i < b->size; ++i) {
			mpz_set(b->elems[i].coef, aa->elems[i].coef);
			b->elems[i].degs = aa->elems[i].degs;
		}
		for ( ; i < aa->size; ++i ) {
			mpz_init_set(b->elems[i].coef, aa->elems[i].coef);
			b->elems[i].degs = aa->elems[i].degs;
		}
	}

	b->size = aa->size;
	b->nvar = aa->nvar; //since not packed this is fine.
}

void convertFromDBZP_AAZ_inp(const DBZP_t* p, AltArrZ_t** aa) {
	if (aa == NULL) {
		return;
	}

	if (p == NULL) {
		if (*aa != NULL) {
			freePolynomial_AAZ(*aa);
			*aa = NULL;
		}
		return;
	}

	int requiredSize = 0;
	for (int i = 0; i <= p->ltX; ++i) {
		if (isZero_DUZP(p->coefsY[i])) {
			continue;
		}
		requiredSize += p->coefsY[i]->lt + 1;
	}

	AltArrZ_t* a;
	int nvar = 2;
	if (*aa == NULL) {
		a = makePolynomial_AAZ(requiredSize, nvar);
		*aa = a;
	} else {
		a = *aa;
		if (a->size < requiredSize) {
			resizePolynomial_AAZ(a, requiredSize);
		}
	}

	int i, j, curSize = 0;
	for (i = a->size; i < requiredSize; ++i) {
		mpz_init(a->elems[i].coef);
	}
	for (i = requiredSize; i < a->size; ++i) {
		mpz_clear(a->elems[i].coef);
	}
	for (i = p->ltX; i >= 0; --i) {
		if (isZero_DUZP(p->coefsY[i])) {
			continue;
		}
		for (j = p->coefsY[i]->lt; j >= 0; --j) {
			if (mpz_sgn(p->coefsY[i]->coefs[j]) == 0) {
				continue;
			}
			mpz_set(a->elems[curSize].coef, p->coefsY[i]->coefs[j]);
			a->elems[curSize].degs = (((degrees_t) i) << EXP_OFFSET_1_V2) | ((degrees_t) j);
			++curSize;
		}
	}
	for (i = curSize; i < requiredSize; ++i) {
		mpz_clear(a->elems[i].coef);
	}
	a->size = curSize;
}

// Return a positive value if a > b, zero if a == b, or a negative value if a < b.
// Ordering is based on lexicographical ordering of leading monomials, then leading coefficient,
// then recursively comparing the ordering of the reductum of each of a and b.
int comparePolynomials_AAZ(const AltArrZ_t* a, const AltArrZ_t* b) {
	if (isZero_AAZ(a)) {
		if (isZero_AAZ(b)) {
			return 0;
		}
		return -1;
	}
	if (isZero_AAZ(b)) {
		return 1;
	}

	if (a->nvar != b->nvar) {
		fprintf(stderr, "BPAS ERROR in comparePolynomials_AAZ: polynomials do not belong to the same polynomial ring.");
		exit(1);
	}

	if (a->unpacked || b->unpacked) {
		const AltArrZ_t* aa;
		if (!a->unpacked) {
			aa = deepCopyPolynomial_AAZ(a);
			unpackExponentVectors_AAZ_inp((AltArrZ_t*)aa);
		} else {
			aa = a;
		}
		const AltArrZ_t* bb;
		if (!b->unpacked) {
			bb = deepCopyPolynomial_AAZ(b);
			unpackExponentVectors_AAZ_inp((AltArrZ_t*)bb);
		} else {
			bb = b;
		}

		int ret = 0;
		int size = a->size < b->size ? a->size : b->size;
		int nvar = a->nvar;
		for (int i = 0; i < size && ret == 0; ++i) {
			ret = compareExponentVectors_unpk(aa->elems[i].degs, bb->elems[i].degs, nvar);
			ret = ret == 0 ? mpz_cmp(aa->elems[i].coef, bb->elems[i].coef) : ret;
		}

		if (ret == 0 && a->size < b->size) {
			//a < b
			ret = -1;
		}
		if (ret == 0 && a->size > b->size) {
			//a > b
			ret = 1;
		}

		if (!a->unpacked) {
			freePolynomial_AAZ((AltArrZ_t*)aa);
		}
		if (!b->unpacked) {
			freePolynomial_AAZ((AltArrZ_t*)bb);
		}

		return ret;
	}

	int ret = 0;
	int size = a->size < b->size ? a->size : b->size;
	for (int i = 0; i < size && ret == 0; ++i) {
		ret = compareExponentVectors(a->elems[i].degs, b->elems[i].degs);
		ret = ret == 0 ? mpz_cmp(a->elems[i].coef, b->elems[i].coef) : ret;
	}

	if (ret == 0 && a->size < b->size) {
		//a < b
		ret = -1;
	}
	if (ret == 0 && a->size > b->size) {
		//a > b
		ret = 1;
	}

	return ret;
}


void sortPolynomialList_AAZ(AltArrZ_t** polys, int n) {
	if (polys == NULL || n == 0) {
		return;
	}

	AltArrZ_t* swapPoly;
	for (int i = 1; i < n; ++i) {
		for (int j = i; j > 0 && comparePolynomials_AAZ(polys[j-1], polys[j]) > 0; --j) {
			swapPoly = polys[j-1];
			polys[j-1] = polys[j];
			polys[j] = swapPoly;
		}
	}
}

void sortPolynomialAndExponentList_AAZ(AltArrZ_t** polys, int* exps, int n) {
	if (polys == NULL || n == 0) {
		return;
	}

	AltArrZ_t* swapPoly;
	int swapExp;
	for (int i = 1; i < n; ++i) {
		for (int j = i; j > 0 && comparePolynomials_AAZ(polys[j-1], polys[j]) > 0; --j) {
			swapPoly = polys[j-1];
			polys[j-1] = polys[j];
			polys[j] = swapPoly;

			swapExp = exps[j-1];
			exps[j-1] = exps[j];
			exps[j] = swapExp;
		}
	}
}

AltArrZ_t* sortPolynomial_AAZ(AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return sortPolynomial_AAZ_unpk(aa);
	}

	AAZElem_t* elems = aa->elems;
	int size = AA_SIZE(aa);

	degrees_t swapDegs;
	for (int i = 1; i < size; ++i) {
		for (int j = i; j > 0 && compareExponentVectors(elems[j-1].degs, elems[j].degs) < 0; --j) {
			mpz_swap(elems[j-1].coef, elems[j].coef);
			swapDegs = elems[j-1].degs;
			elems[j-1].degs = elems[j].degs;
			elems[j].degs = swapDegs;
		}
	}

	condensePolynomial_AAZ(aa);

	return aa;
}

static void mergeAAZElems(AAZElem_t* __restrict__ a, AAZElem_t* __restrict__ endA, AAZElem_t* __restrict__ b , AAZElem_t* __restrict__ endB, AAZElem_t* __restrict__ sorted) {
	int i = 0;
	while (a < endA && b < endB) {
		if (isGreaterExponentVectors(a->degs, b->degs)) {
			sorted[i] = *a;
			++a;
		} else {
			sorted[i] = *b;
			++b;
		}
		++i;
	}

	while (a < endA) {
		sorted[i] = *a;
		++a;
		++i;
	}

	while (b < endB) {
		sorted[i] = *b;
		++b;
		++i;
	}
}

void mergeSortPolynomial_AAZ(AltArrZ_t* aa) {
	if (aa->size < 9) {
		sortPolynomial_AAZ(aa);
	}

	if (aa->unpacked) {
		mergeSortPolynomial_AAZ_unpk(aa);
		return;
	}

	if (aa->size < 8) {
		sortPolynomial_AAZ(aa);
		return;
	}

	int size = aa->size;
	AAZElem_t* tempElems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*aa->alloc);
	AAZElem_t* elems = aa->elems;
	int end1 = 0;
	int end2 = 0;
	for (int window = 1; window < size; window <<= 1) {
		//merge elems[i]:elems[i+window] with elems[i+window]:elems[i+2*window]
		for (int i = 0; i < size; i += 2*window) {
			end1 = i + window < size ? i + window : size;
			end2 = i + 2*window < size ? i + 2*window : size;
			mergeAAZElems(elems+i, elems+end1, elems+end1, elems+end2, tempElems+i);
		}

		AAZElem_t* temp = tempElems;
		tempElems = elems;
		elems = temp;
	}

	aa->elems = elems;

	free(tempElems);
}

void condensePolynomial_AAZ(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) <= 1) {
		return;
	}

	if (aa->unpacked) {
		condensePolynomial_AAZ_unpk(aa);
		return;
	}

	int size = AA_SIZE(aa);
	AAZElem_t* elems = aa->elems;
	int insertIdx = 0;
	int compareIdx = 1;
	while (compareIdx < size) {
		if (compareExponentVectors(elems[insertIdx].degs, elems[compareIdx].degs) == 0) {
			mpz_add(elems[insertIdx].coef, elems[insertIdx].coef, elems[compareIdx].coef);
		} else if(compareIdx - insertIdx > 1) {
			if (mpz_sgn(elems[insertIdx].coef) != 0) {
				++insertIdx;
			}
			elems[insertIdx].degs = elems[compareIdx].degs;
			mpz_swap(elems[insertIdx].coef, elems[compareIdx].coef);
		} else if (mpz_sgn(elems[insertIdx].coef) != 0) {
			++insertIdx;
		}
		++compareIdx;
	}

	if (mpz_sgn(elems[insertIdx].coef) != 0) {
		++insertIdx;
	}

	for (int i = insertIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
	}
	AA_SIZE(aa) = insertIdx;
}

void removeZeroTerms_AAZ(AltArrZ_t* aa) {
	if (aa == NULL || aa->size <= 1) {
		return;
	}

	if (aa->unpacked) {
		removeZeroTerms_AAZ_unpk(aa);
		return;
	}

	int size = AA_SIZE(aa);
	AAZElem_t* elems = aa->elems;
	int curIdx = 0;
	for (int i = 0; i < size; ++i) {
		if (mpz_cmp_ui(elems[i].coef, 0ul) != 0) {
			if (i != curIdx) {
				mpz_swap(elems[curIdx].coef, elems[i].coef);
				elems[curIdx].degs = elems[i].degs;
			}
			++curIdx;
		}
	}

	//keep one zero around if aa is not null and indentically 0.
	if (curIdx == 0) {
		elems[curIdx].degs = 0ll;
		++curIdx;
	}

	for (int i = curIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
	}

	aa->size = curIdx;
}

void negatePolynomial_AAZ(AltArrZ_t* aa) {
	if (isZero_AAZ (aa)) {
		return;
	}
	int size = AA_SIZE(aa);
	AAZElem_t* elems = aa->elems;
	for (int i = 0; i < size; ++i) {
		mpz_neg(elems[i].coef, elems[i].coef);
	}
}

AltArrZ_t* mainLShiftPolynomial_AAZ (AltArrZ_t* aa, int n)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}
	if (n < 1){
		return deepCopyPolynomial_AAZ (aa);
	}

	if (aa->unpacked) {
		return mainLShiftPolynomial_AAZ_unpk(aa, n);
	}

	int nvar = aa->nvar;
	int mvarDegOffset = getMVarExpOffset (nvar);

	degrees_t nDegs = (degrees_t) n << mvarDegOffset;
	int fitsPacked = monomialMultFitsPacked(aa->elems->degs, nDegs, nvar);
	if (!fitsPacked) {
		unpackExponentVectors_AAZ_inp(aa);
		AltArrZ_t* ret = mainLShiftPolynomial_AAZ_unpk(aa, n);
		packExponentVectors_AAZ_inp(aa);
		return ret;
	}
	AltArrZ_t* shifted = (AltArrZ_t*) malloc (sizeof (AltArrZ_t));
	shifted->nvar = aa->nvar;
	shifted->size = aa->size;
	shifted->alloc = aa->alloc;
	shifted->unpacked = 0;
	shifted->elems = (AAZElem_t*) malloc (sizeof (AAZElem_t)*(aa->size));
	AAZElem_t* selems = shifted->elems;
	AAZElem_t* aelems = aa->elems;

	for (int i = 0; i < shifted->size; ++i){
		mpz_init (selems[i].coef);
		mpz_set (selems[i].coef, aelems[i].coef);
		selems[i].degs = aelems[i].degs;
		selems[i].degs += ((degrees_t) n << mvarDegOffset);
	}
	return shifted;
}

AltArrZ_t* mainLShiftPolynomial_AAZ_inp (AltArrZ_t* aa, int n)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}
	if (n < 1){
		return aa;
	}

	if (aa->unpacked) {
		return mainLShiftPolynomial_AAZ_inp_unpk(aa, n);
	}

	int nvar = aa->nvar;
	int mvarDegOffset = getMVarExpOffset (nvar);

	degrees_t nDegs = (degrees_t) n << mvarDegOffset;
	int fitsPacked = monomialMultFitsPacked(aa->elems->degs, nDegs, nvar);
	if (!fitsPacked) {
		unpackExponentVectors_AAZ_inp(aa);
		return mainLShiftPolynomial_AAZ_inp_unpk(aa, n);
	}

	AAZElem_t* aelems = aa->elems;

	for (int i = 0; i < aa->size; ++i){
		aelems[i].degs += ((degrees_t) n << mvarDegOffset);
	}
	return aa;
}

AltArrZ_t* multiplyByInteger_AAZ(AltArrZ_t* a, const mpz_t mult) {
	if (isZero_AAZ(a)) {
		return NULL;
	}

	AltArrZ_t* ret = deepCopyPolynomial_AAZ(a);
	AAZElem_t* elems = ret->elems;
	for (int i = 0; i < ret->size; ++i) {
		mpz_mul(elems[i].coef, elems[i].coef, mult);
	}
	return ret;
}

AltArrZ_t* divideByInteger_AAZ(AltArrZ_t* a, mpz_t div) {
	if (isZero_AAZ(a)) {
		return NULL;
	}

	AltArrZ_t* ret = deepCopyPolynomial_AAZ(a);
	divideByInteger_AAZ_inp(ret, div);
	return ret;
}


void divideByInteger_AAZ_inp(AltArrZ_t* ret, mpz_t div) {
	if (isZero_AAZ(ret)) {
		return;
	}

	AAZElem_t* retElems = ret->elems;
	for (int i = 0; i < ret->size; ++i) {
		if (!mpz_divisible_p(retElems[i].coef, div)) {
			fprintf(stderr, "SMZP divideByInteger ERROR!\n");
			exit(1);
		}
		mpz_div(retElems[i].coef, retElems[i].coef, div);
	}
}

void multiplyByMonomial_AAZ_inp(AltArrZ_t* aa, degree_t* __restrict__ degs) {
	if (isZero_AAZ(aa) || degs == NULL) {
		return;
	}

	register int size = aa->size;
	register int nvar = aa->nvar;
	register int i;
	register int j;

	int done = 0;
	if (!aa->unpacked) {
		if (isPackableExponentVector(degs, nvar)) {
			const int* sizes = getExpOffsetArray(nvar);
			degrees_t packedDegs = 0ll;
			for (j = 0; j < nvar; ++j) {
				packedDegs |= ((degrees_t) degs[j]) << sizes[j];
			}
			degrees_t maxDegs = calculateMaxDegs_AAZ(aa);
			if (monomialMultFitsPacked(maxDegs, packedDegs, nvar)) {
				for (i = 0; i < size; ++i) {
					addExponentVectors(aa->elems[i].degs, packedDegs, aa->elems[i].degs);
				}
				done = 1;
			} else {
				unpackExponentVectors_AAZ_inp(aa);
			}
		} else {
			unpackExponentVectors_AAZ_inp(aa);
		}
	}

	if(!done) {
		degree_t* __restrict__ aDegs = (degree_t*) aa->elems->degs;
		for (i = 0; i < size; ++i) {
			for (j = 0; j < nvar; ++j) {
				aDegs[i*nvar + j] += degs[j];
			}
		}
	}
}

void applyModuloSymmetric_AAZ_inp(AltArrZ_t* p, const mpz_t mod) {
	if (isZero_AAZ(p)) {
		return;
	}

	polysize_t newSize = 0;
	mpz_t halfMod;
	mpz_init_set(halfMod, mod);
	mpz_div_ui(halfMod, halfMod, 2ul);

	int unpacked = p->unpacked;
	int nvar = p->nvar;
	polysize_t i;
	polysize_t insertIdx = 0;
	polysize_t size = p->size;
	AAZElem_t* elems = p->elems;
	for (i = 0; i < size; ++i) {
		// fprintf(stderr, "moding poly %p at %ld\n",p, i);
		mpz_mod(elems[insertIdx].coef, elems[i].coef, mod);
		if(mpz_cmp(elems[insertIdx].coef, halfMod) > 0) {
			mpz_sub(elems[insertIdx].coef, elems[insertIdx].coef, mod);
		}

		if (unpacked) {
			 setExponentVector_unpk(elems[insertIdx].degs, elems[i].degs, nvar);
		} else {
			elems[insertIdx].degs = elems[i].degs;
		}

		if (mpz_sgn(elems[insertIdx].coef) != 0) {
			++insertIdx;
		}
	}
	for (i = insertIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
	}

	p->size = insertIdx;
	mpz_clear(halfMod);
}

void evalPolyToVal_AAZ(const AltArrZ_t* aa, const mpz_t* vals, short nvar, mpz_t res) {
	if (aa == NULL || nvar != aa->nvar) {
		mpz_set_ui(res, 0ul);
		return;
	}

	if (nvar == 0 || aa->nvar == 0) {
		mpz_set(res, aa->elems->coef);
		return;
	}

	if (aa->unpacked) {
		evalPolyToVal_AAZ_unpk(aa, vals, nvar, res);
		return;
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);

	mpz_t* valList[nvar];
	int valListSize[nvar];

	degrees_t maxDegs = calculateMaxDegs_AAZ(aa);
	for (int j = 0; j < nvar; ++j) {
		degree_t deg = GET_NTH_EXP(maxDegs, masks[j], sizes[j]);
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(deg+1));
		valListSize[j] = deg+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	mpz_set_ui(res, 0ul);
	mpz_t acc;
	mpz_init(acc);

	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpz_set(acc, aa->elems[i].coef);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = GET_NTH_EXP(aa->elems[i].degs, masks[j], sizes[j]);
			mpz_mul(acc, acc, valList[j][deg]);
		}
		mpz_add(res, res, acc);
	}

	mpz_clear(acc);

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}
}

AltArrZ_t* evaluatePoly_AAZ(const AltArrZ_t* aa, const int* active, const mpz_t* vals, short nvar) {
	if (aa == NULL) {
		return NULL;
	}

	if (nvar == 0 || aa->nvar == 0) {
		return deepCopyPolynomial_AAZ(aa);
	}

	if (aa->unpacked) {
		return evaluatePoly_AAZ_unpk(aa, active, vals, nvar);
	}

	int i;
	int newNvar = nvar;
	for (i = 0; i < nvar; ++i) {
		if (active[i]) {
			--newNvar;
		}
	}

	if (newNvar == 0) {
		mpz_t res;
		mpz_init(res);
		evalPolyToVal_AAZ(aa, vals, nvar, res);
		AltArrZ_t* ret = makeConstPolynomial_AAZ(1, 0, res);
		mpz_clear(res);
		return ret;
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);

	mpz_t* valList[nvar];
	int valListSize[nvar];

	degrees_t maxDegs = calculateMaxDegs_AAZ(aa);
	for (int j = 0; j < nvar; ++j) {
		degree_t deg = GET_NTH_EXP(maxDegs, masks[j], sizes[j]);
		if (!active[j]) {
			valList[j] = NULL;
			valListSize[j] = 0;
			continue;
		}
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(deg+1));
		valListSize[j] = deg+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	AltArrZ_t* res = makePolynomial_AAZ(aa->size, newNvar);
	AAZElem_t* resElems = res->elems;
	res->size = aa->size;

	const int* newSizes = getExpOffsetArray(newNvar);
	int size = aa->size;
	int k = 0;
	degrees_t newDegs = 0;
	for (int i = 0; i < size; ++i) {
		mpz_init(resElems[i].coef);
		mpz_set(resElems[i].coef, aa->elems[i].coef);
		for (int j = 0; j < nvar; ++j) {
			degrees_t deg = GET_NTH_EXP(aa->elems[i].degs, masks[j], sizes[j]);
			if (valList[j] == NULL) {
				newDegs |= (deg << newSizes[k]);
				++k;
			} else {
				mpz_mul(resElems[i].coef, resElems[i].coef, valList[j][deg]);
			}
		}
		resElems[i].degs = newDegs;
		k = 0;
		newDegs = 0;
	}

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}

	canonicalizePolynomial_AAZ(res);

	return res;
}


AltArrZ_t* convertFromAAZElemToAAZ (AAZElem_t* coef, int coefSize, int nvar, int unpacked)
{
	if (coef == NULL || coefSize == 0) {
	  return NULL;
	}

	AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	poly->size = coefSize;
	poly->alloc = coefSize;
	poly->elems = coef;
	poly->nvar = nvar;
	poly->unpacked = unpacked;
	return poly;
}

AltArrZ_t* swappingExponents_AAZ (AltArrZ_t* aa, int idx1, int idx2)
{
  if (aa == NULL)
	  return NULL;


  AltArrZ_t* cPoly = deepCopyPolynomial_AAZ (aa);
  if (idx1 == idx2){
	  return cPoly;
  }

  int varMap[cPoly->nvar];
  for (int i = 0; i < cPoly->nvar; ++i){
	varMap[i] = i;
  }
  varMap[idx1] = idx2;
  varMap[idx2] = idx1;

  reorderVars_AAZ (cPoly, varMap, cPoly->nvar);
  return cPoly;
}


AltArrZ_t* leadingTerm_AAZ (AltArrZ_t* aa, int nvar)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	if (aa->unpacked) {
		return leadingTerm_AAZ_unpk(aa, nvar);
	}

	AltArrZ_t* lt = makeConstPolynomial_AAZ (1, aa->nvar, aa->elems->coef);
	lt->elems->degs = aa->elems->degs;

	return lt;
}

int leadingVariable_AAZ (const AltArrZ_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return -2;
	}

	if (aa->unpacked) {
		return leadingVariable_AAZ_unpk(aa);
	}

	if (aa->nvar == 0) {
		return -1;
	}

	const degrees_t* masks = getExpMaskArray(aa->nvar);
	for (int i = 0; i < aa->nvar; ++i){
		if ((aa->elems[0].degs & masks[i]) != 0) {
			return i;
		}
	}

	return -1;
}


int mainLeadingDegree_AAZ (AltArrZ_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return 0;
	}

	if (aa->unpacked) {
		return mainLeadingDegree_AAZ_unpk(aa);
	}

	int mvarDegOffset = getMVarExpOffset (aa->nvar);
	degrees_t mvarMask = getMVarExpMask (aa->nvar);

	return GET_NTH_EXP(aa->elems->degs, mvarMask, mvarDegOffset);
}


AltArrZ_t* mainLeadingCoefficient_AAZ (const AltArrZ_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	if (aa->unpacked) {
		return mainLeadingCoefficient_AAZ_unpk(aa);
	}

	int mvarDegOffset = getMVarExpOffset (aa->nvar);
	degrees_t mvarMask = getMVarExpMask(aa->nvar);

	AAZElem_t* elems = aa->elems;
	degree_t mvarDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
	degree_t curDeg = mvarDeg;

	int polyElemsAlloc = 10;
	int polyElemsSize = 0;
	AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);

	for (int i = 0; (mvarDeg == curDeg) && i < AA_SIZE (aa); ++i){

		if (polyElemsSize + 1 > polyElemsAlloc){
			polyElemsAlloc += 10;
			polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t) * polyElemsAlloc);
		}

		mpz_init (polyElems[i].coef);
		mpz_set (polyElems[i].coef, elems[i].coef);
		polyElems[i].degs = (elems[i].degs & (~mvarMask));
		++polyElemsSize;

		if (i+1 < AA_SIZE(aa)){
			curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
		} else {
			curDeg = -1;
		}
	}

	AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	poly->alloc = polyElemsAlloc;
	poly->size = polyElemsSize;
	poly->elems = polyElems;
	poly->unpacked = 0;
	poly->nvar = aa->nvar;

	return poly;
}

AltArrZ_t* mainCoefficientAtIdx_AAZ (AltArrZ_t* aa, int e)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}
	if (e < 0){
		return NULL;
	}

	if (aa->unpacked) {
		return mainCoefficientAtIdx_AAZ_unpk(aa, e);
	}

	int mvarDegOffset = getMVarExpOffset (aa->nvar);
	degrees_t mvarMask = getMVarExpMask(aa->nvar);

	AAZElem_t* elems = aa->elems;
	degree_t mvarDeg = e;
	degree_t curDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
	int begin = -1;

	if (mvarDeg > curDeg){
		return NULL;
	}

	int polyElemsAlloc = 10;
	int polyElemsSize = 0;
	AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);

	for (int i = 0; i < AA_SIZE (aa); ++i){
		if (mvarDeg == curDeg){
			if (polyElemsSize + 1 > polyElemsAlloc) {
				polyElemsAlloc += 10;
				polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t)*polyElemsAlloc);
			}

			mpz_init (polyElems[polyElemsSize].coef);
			mpz_set (polyElems[polyElemsSize].coef, elems[i].coef);
			polyElems[polyElemsSize].degs = elems[i].degs & (~mvarMask);
			++polyElemsSize;

			if (begin == -1) {
				begin = i;
			}
		}

		if (begin != -1 && mvarDeg != curDeg){
			break;
		}

		if (i+1 < AA_SIZE(aa)){
			curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
		} else {
			curDeg = -1;
		}
	}

	AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	poly->alloc = polyElemsAlloc;
	poly->size = polyElemsSize;
	poly->elems = polyElems;
	poly->unpacked = 0;
	poly->nvar = aa->nvar;

	return poly;
}

void mainCoefficientListAtIdx_AAZ (AltArrZ_t* aa, int idx, AltArrZ_t*** cList, int *sz)
{
	if (aa == NULL || aa->size == 0) {
		*sz = 0;
		*cList = NULL;
		return;
	}

	if (idx < 0) {
		*sz = 0;
		*cList = NULL;
		return;
	}

	if (aa->unpacked) {
		// TODO
		fprintf (stderr, "mainCoefficientListAtIdx_AAZ_unpk is not implemented yet!\n");
		exit(1);
	}


	AltArrZ_t* AA = aa;
	int isExpand = 0;
	int nnvar;

	if (idx > 0) {
		AA = deepCopyPolynomial_AAZ (aa);

		nnvar = AA->nvar + 1;
		expandNumVarsLeft_AAZ (AA, nnvar);

		int varmap [nnvar];
		for (int i = 0; i < nnvar; i++) {
			varmap[i] = i;
		}
		varmap[0] = idx+1;
		varmap[idx+1] = 0;
		reorderVars_AAZ (AA, varmap, nnvar);

		isExpand = 1;
	} else {
		nnvar = AA->nvar;
	}

	int mvarDegOffset = getMVarExpOffset (nnvar);
	degrees_t mvarMask = getMVarExpMask(nnvar);

	AAZElem_t* elems = AA->elems;
	degree_t mvarDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
	degree_t curDeg =  mvarDeg;
	degree_t expDeg;

	AltArrZ_t** mCoefs = (AltArrZ_t**) calloc (curDeg+1, sizeof (AltArrZ_t*));
	int len = 0;

	AltArrZ_t* poly;
	int polyElemsAlloc;
	int polyElemsSize;
	AAZElem_t* polyElems;
	int start = 0;

	for (expDeg = mvarDeg; expDeg >= 0; expDeg--) {

		if (expDeg != curDeg) {
			mCoefs[expDeg] = NULL;
			len++;

			continue;
		}

		/* fprintf (stderr, "[SMQP] expDeg = %d\n", expDeg); // TEST */

		polyElemsAlloc = 10;
		polyElemsSize = 0;
		polyElems = (AAZElem_t*) malloc (polyElemsAlloc * sizeof(AAZElem_t)); // TODO
		for (int i = start; (curDeg == expDeg) && i < AA_SIZE(AA); i++) {

			if (polyElemsSize + 1 > polyElemsAlloc){
				polyElemsAlloc += 10;
				polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t) * polyElemsAlloc);
			}

			/* fprintf (stderr, "[SMQP] i = %d\n", i); // TEST */
			mpz_init (polyElems[polyElemsSize].coef);
			mpz_set (polyElems[polyElemsSize].coef, elems[i].coef);
			polyElems[polyElemsSize].degs = (elems[i].degs & (~mvarMask));
			++polyElemsSize;

			if (i+1 < AA_SIZE(AA)){
				curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
			} else {
				curDeg = -1;
			}

			start++;
		}

		poly = (AltArrZ_t*) malloc (sizeof (AltArrZ_t));
		poly->alloc = polyElemsAlloc;
		poly->size = polyElemsSize;
		poly->elems = polyElems;
		poly->unpacked = 0;
		poly->nvar = nnvar;

		if (isExpand) {
			shrinkNumVarsAtIdx_AAZ (poly, 0);
		}

		/* printAAZ (poly); // TEST		 */
		mCoefs[expDeg] = poly;
		len++;
	}

	*cList = mCoefs;
	*sz = len;

}

// TODO: used in recursive MDPD
AltArrZ_t* maxPolynomials_AAZ (AltArrZ_t* a, AltArrZ_t* b)
{
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AAZ (b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ (a);
	}

	if (a->unpacked || b->unpacked) {
		return maxPolynomials_AAZ_unpk(a,b);
	}

	int cmp = compareExponentVectors(a->elems[0].degs, b->elems[0].degs);

	if (cmp > 0){
		return deepCopyPolynomial_AAZ (a);
	} else if (cmp < 0){
		return deepCopyPolynomial_AAZ (b);
	} else {
		if (mpz_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return deepCopyPolynomial_AAZ (a);
		} else {
			return deepCopyPolynomial_AAZ (b);
		}
	}
	return NULL;
}

AltArrZ_t* maxPolynomials_AAZ_inp (AltArrZ_t* a, AltArrZ_t* b)
{
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AAZ (b);
	}
	if (b == NULL || b->size == 0){
		return a;
	}

	if (a->unpacked || b->unpacked) {
		return maxPolynomials_AAZ_inp_unpk(a,b);
	}

	int cmp = compareExponentVectors(a->elems[0].degs, b->elems[0].degs);

	if (cmp > 0){
		return a;
	} else if (cmp < 0){
		return deepCopyPolynomial_AAZ (b);
	} else {
		if (mpz_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return a;
		} else {
			return deepCopyPolynomial_AAZ (b);
		}
	}
	return NULL;
}



/*****************
 * SMZP Addition
 *****************/

void addInteger_AAZ_inp(AltArrZ_t* aa, const mpz_t coef) {
	if (aa == NULL) {
		return;
	}

	if (aa->size == 0) {
		if (aa->alloc < 1) {
			resizePolynomial_AAZ(aa, 1);
		}
		mpz_set(aa->elems->coef, coef);
		return;
	}

	if (aa->unpacked) {
		addInteger_AAZ_inp_unpk(aa, coef);
		return;
	}

	if (isZeroExponentVector(aa->elems[aa->size-1].degs)) {
		mpz_add(aa->elems[aa->size-1].coef, aa->elems[aa->size-1].coef, coef);
		return;
	}

	if (aa->size >= aa->alloc) {
		resizePolynomial_AAZ(aa, aa->alloc+10);
	}

	mpz_init(aa->elems[aa->size].coef);
	mpz_set(aa->elems[aa->size].coef, coef);
	aa->elems[aa->size].degs = 0;
	++(aa->size);
}

AltArrZ_t* addPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AAZ(a);
	}

	if (a->unpacked || b->unpacked) {
		return addPolynomials_AAZ_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			mpz_init(cElems[k].coef);
			mpz_add(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, aElems[i].coef);
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);
	if (isZero_AAZ(c)) {
		freePolynomial_AAZ(c);
		c = NULL;
	}

	return c;
}


AltArrZ_t* addPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		if (a != NULL && a->size == 0) {
			freePolynomial_AAZ(a);
		}
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0) {
		return a;
	}

	if (a->unpacked || b->unpacked) {
		return addPolynomials_AAZ_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;
	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = aElems[i];
			mpz_add(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			cElems[k] = aElems[i];
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		k += (asize - i);
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;

	resizePolynomial_AAZ(c, k);
	if (isZero_AAZ(c)) {
		freePolynomial_AAZ(c);
		c = NULL;
	}


	return c;
}

AltArrZ_t* subPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (isZero_AAZ (a) && isZero_AAZ (b)) {
		return NULL;
	}
	if (isZero_AAZ (a)) {
		AltArrZ_t* ret = deepCopyPolynomial_AAZ(b);
		negatePolynomial_AAZ(ret);
		return ret;
	}
	if (isZero_AAZ (b)) {
		return deepCopyPolynomial_AAZ(a);
	}

	if (a->unpacked || b->unpacked) {
		return subPolynomials_AAZ_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			mpz_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, aElems[i].coef);
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);
	if (isZero_AAZ(c)) {
		freePolynomial_AAZ(c);
		c = NULL;
	}

	return c;
}

AltArrZ_t* subPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (isZero_AAZ (a) && isZero_AAZ (b)) {
		return NULL;
	}
	if (isZero_AAZ (a)) {
		if (a != NULL && a->size == 0) {
			freePolynomial_AAZ(a);
		}
		AltArrZ_t* ret = deepCopyPolynomial_AAZ(b);
		negatePolynomial_AAZ(ret);
		return ret;
	}
	if (isZero_AAZ (b)) {
		return a;
	}

	if (a->unpacked || b->unpacked) {
		return subPolynomials_AAZ_inp_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;

	resizePolynomial_AAZ(c, k);
	if (isZero_AAZ(c)) {
		freePolynomial_AAZ(c);
		c = NULL;
	}

	return c;
}


void subPolynomials_AAZ_inpRHS(const AltArrZ_t* a, AltArrZ_t** bb) {
	if (bb == NULL) {
		return;
	}

	if (isZero_AAZ(a)) {
		negatePolynomial_AAZ(*bb);
		return;
	}

	if (isZero_AAZ(*bb)) {
		deepCopyPolynomial_AAZ_inp(a, bb);
	}

	if (a->nvar != (*bb)->nvar) {
		fprintf(stderr, "Polynomials are not compatible in subPolynomials_AAZ_inpRHS\n" );
		exit(1);
	}

	AltArrZ_t* b = *bb;
	if (a->unpacked || b->unpacked) {
		subPolynomials_AAZ_inpRHS_unpk(a, bb);
		return;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AAZElem_t* cElems = (AAZElem_t*) malloc(sizeof(AAZElem_t) * (asize + bsize));

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;

	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			cElems[k] = bElems[j];
			mpz_neg(cElems[k].coef, cElems[k].coef);
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = bElems[j];
			// mpq_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, aElems[i].coef, cElems[k].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is b[j], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, aElems[i].coef);
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if(j < bsize) {
		memcpy(cElems + k, bElems + j, sizeof(AAZElem_t)*(bsize - j));
		while ( j < bsize ) {
			mpz_neg(cElems[k].coef, cElems[k].coef);
			++j;
			++k;
		}
	}
	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(bElems);
	b->elems = cElems;
	b->size = k;
	b->alloc = asize + bsize;

	resizePolynomial_AAZ(b, k);

	//TODO make sure users of this method are okay with the result somtimes being NULL.
	// if (isZero_AAZ(b)) {
	// 	freePolynomial_AAZ(b);
	// 	b = NULL;
	// }
	// *bb = NULL;

}



AltArrZ_t* CFDucosOptZ_subPolynomials_AAZ_inp (AltArrZ_t* a, AltArrZ_t* b, int nvar, AltArrZ_t** Pe, int e)
{
	if (a == NULL && b == NULL) {
		*Pe = NULL;
		return NULL;
	}

	if (a->unpacked || b->unpacked) {
		return CFDucosOptZ_subPolynomials_AAZ_inp_unpk(a, b, nvar, Pe, e);
	}

	register int asize = a == NULL ? 0 : AA_SIZE(a);
	register int bsize = b == NULL ? 0 : AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a == NULL ? NULL : a->elems;
	AAZElem_t* bElems = b == NULL ? NULL : b->elems;
	AAZElem_t* cElems = c->elems;

	int mvarDegOffset = getMVarExpOffset (a->nvar);
	degrees_t mvarMask = getMVarExpMask(a->nvar);

	int polyElemsAlloc = 10;
	int polyElemsSize = 0;
	AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
		//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
		  // a==b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;

	resizePolynomial_AAZ(c, k);

	if (c == NULL || c->size == 0 || e < 0){
		*Pe = NULL;
	}

	AAZElem_t* elems = c->elems;
	degree_t mvarDeg = e;
	degree_t curDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
	int begin = -1;

	for (int i = 0; i < AA_SIZE (c); ++i){
		if (mvarDeg == curDeg){
			if (polyElemsSize + 1 > polyElemsAlloc) {
				polyElemsAlloc += 10;
				polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t)*polyElemsAlloc);
			}

			mpz_init (polyElems[polyElemsSize].coef);
			mpz_set (polyElems[polyElemsSize].coef, elems[i].coef);
			polyElems[polyElemsSize].degs = elems[i].degs & (~mvarMask);
			++polyElemsSize;

			if (begin == -1) {
				begin = i;
			}
		}

		if (begin != -1 && mvarDeg != curDeg){
			break;
		}

		if (i+1 < AA_SIZE(c)){
			curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
		} else {
			curDeg = -1;
		}
	}

	AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	poly->alloc = polyElemsAlloc;
	poly->size = polyElemsSize;
	poly->elems = polyElems;
	poly->unpacked = 0;
	poly->nvar = a->nvar;

	*Pe = poly;
	return c;
}

void subByLeadingTerm_AAZ (AltArrZ_t** a, int nvar)
{
	AltArrZ_t* aa = *a;

	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (AA_SIZE (aa) == 1) {
		freePolynomial_AAZ (aa);
		*a = NULL;
		return;
	}

	if (aa->unpacked) {
		AltArrZ_t* lt = leadingTerm_AAZ (aa, nvar);
		*a = subPolynomials_AAZ_inp (aa, lt, nvar);
		freePolynomial_AAZ (lt);
		return;
	}

	AltArrZ_t* ret = makePolynomial_AAZ (AA_ALLOC(aa), nvar);
	memcpy (ret->elems, aa->elems + 1, sizeof (AAZElem_t)*(AA_SIZE(aa)-1));

	AA_SIZE (ret) = AA_SIZE (aa) - 1;

	mpz_clear (aa->elems->coef);
	free (aa);

	*a = ret;
}


/*****************
 * SMZP Multiplication & Helpers
 *****************/

#define PRODHEAP_LAZY_INIT 1

ProductHeap_AAZ* prodheapInit_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, int nvar) {
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	h->elements = (ProductHeapElem_AAZ*) malloc(sizeof(ProductHeapElem_AAZ)*AA_SIZE(a));

	h->elements[0].chain = prodheapMakeChain_AAZ(0, 0, NULL);
	addExponentVectors(a->elems->degs, b->elems->degs, h->elements->degs);

	h->heapSize = 1;
	h->maxHeapSize = AA_SIZE(a);
	h->nvar = nvar;

#if PRODHEAP_LAZY_INIT
#else
	int asize = a->size;
	AAZElem_t* aElems = a->elems;
	degrees_t bDegs = b->elems->degs;
	ProductHeapChain_AAZ* ch;
	for (int i = 1; i < asize; ++i) {
		ch = prodheapMakeChain_AAZ(i, 0, NULL);
		prodheapInsert_AAZ(h, ch, aElems[i].degs + bDegs);
	}
#endif
	return h;
}

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
void prodheapInsert_AAZ(ProductHeap_AAZ* h, ProductHeapChain_AAZ* chain, register degrees_t degs) {
	register int s = h->heapSize;
	ProductHeapElem_AAZ* elems = h->elements;

	if (s == 0) {
		elems[0].degs = degs;
		elems[0].chain = chain;
		h->heapSize = 1;
		return;
	}

	//first check if we can chain off the root
	if (isEqualExponentVectors(elems[0].degs, degs)) {
		chain->next = elems[0].chain;
		elems[0].chain = chain;
		return;
	}

	//otherwise, we must search the heap to find the new product's insertion point
	//note that since we are looking for chains we cannot use the simple swim method
	//we sort of fake the swimming, looking for a chain to be made or the eventual
	//place where the swim would stop. At this point, we insert the new elem
	//in that spot, and "push" the entire path we took down a level. Assuming
	//that we insert e and it ends up at the root, we push down the 'x' path
	//                                      //
	//     x     --->    e                  //
	//    / \           / \                 //
	//   x   o         x   o
	//                /
	//               x

	register int i = (s-1) >> 1; //i is parent
	register int j = s;       //j is current insertion point
	register long long unsigned int path = 1;
	while (j > 0) {
		if (isEqualExponentVectors(elems[i].degs, degs)) {
			chain->next = elems[i].chain;
			elems[i].chain = chain;
			return;
		} else if (isLessExponentVectors(elems[i].degs, degs)) {
			path <<= 1;
			if (!(j & 1)) {
				//set the trailing bit to 1 to distinguish left/right of path
				path += 1;
			}
			j = i;
			i = (i-1) >> 1;
		} else { //cmp > 0
			break;
		}
	}

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	ProductHeapElem_AAZ temp;
	ProductHeapElem_AAZ elem = {degs, chain};

	//TODO use i index again to swap between elements rather than use a second temp elem.
	while (j <= s) {
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = (j << 1) + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
ProductHeapChain_AAZ* prodheapRemoveMax_AAZ(ProductHeap_AAZ* h) {
	ProductHeapElem_AAZ* elems = h->elements;
	ProductHeapChain_AAZ* maxElem = elems[0].chain;
	register int i = 0;
	register int j = 1;
	register int s = --(h->heapSize);

	//promote largest children
	while (j < s) {
		if (j+1 < s && isLessExponentVectors(elems[j].degs, elems[j+1].degs)) {
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (j << 1) + 1;
	}
	//now place last element into i and swim up to make tree complete
	j = (i-1) >> 1;
	while(i > 0) {
		if (isLessExponentVectors(elems[s].degs, elems[j].degs)) {
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j-1) >> 1;
	}
	elems[i] = elems[s];

	return maxElem;
}

/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure.
 * It is assumed that both a and b have compatible exponent vectors.
 *
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
AltArrZ_t* multiplyPolynomials_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, int nvar) {
	if (a == NULL || a->size == 0 || b == NULL || b->size == 0) {
		return NULL;
	}

	if (a->unpacked || b->unpacked) {
		return multiplyPolynomials_AAZ_unpk(a,b,nvar);
	}

	degrees_t aMax = calculateMaxDegs_AAZ(a);
	degrees_t bMax = calculateMaxDegs_AAZ(b);

	if (!monomialMultFitsPacked(aMax, bMax, nvar)) {
		//force unpk to run by unpacking here.
		// fprintf(stderr, "forcing an unpack before multiply!\n");
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) a);
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) b);
		AltArrZ_t* ret = multiplyPolynomials_AAZ_unpk(a,b,nvar);
		packExponentVectors_AAZ_inp((AltArrZ_t*) a);
		packExponentVectors_AAZ_inp((AltArrZ_t*) b);
		return ret;
	}

	// reorder to obtain smaller as a.
	if (b->size < a->size) {
		const AltArrZ_t* temp = a;
		a = b;
		b = temp;
	}

	ProductHeap_AAZ* h = prodheapInit_AAZ(a,b,nvar);
	// mpz_t ccoef;
	// mpz_init(ccoef);

	//TODO smarter allocation here? dynamic reallocating?
	register unsigned long long int allocC = AA_SIZE(a)*AA_SIZE(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;
	AltArrZ_t* c = makePolynomial_AAZ(allocC, nvar);
// fprintf(stderr, "alloced c at  : %llu\n", allocC);

	//k is c, i is a, j is b.
	register int k = 0;
	// register int i = 0;
	// register int j = 0;

	AAZElem_t* __restrict__ cElems = c->elems;
	AAZElem_t* __restrict__ bElems = b->elems;

	AAZElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1;
	register int firstB = 0;

	ProductHeapChain_AAZ* maxElem = NULL;
	ProductHeapChain_AAZ* nextMaxElem = NULL;
	degrees_t* nextDegs;
	while ( (nextDegs = prodheapPeek_AAZ(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		cElems[k].degs = *nextDegs;
		mpz_init(cElems[k].coef);

		while (nextDegs != NULL && isEqualExponentVectors(cElems[k].degs, *nextDegs)) {
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AAZ* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AAZ(h);
			while (maxElem != NULL) {

				mpz_addmul(cElems[k].coef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				// mpz_mul(ccoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				// mpz_add(cElems[k].coef, cElems[k].coef, ccoef);

#if PRODHEAP_LAZY_INIT
				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA) {
					oldMaxElem = prodheapMakeChain_AAZ((maxElem->a_i)+1, firstB, oldMaxElem);
				}
#endif

				//cache next before freeing or overwriting
				nextMaxElem = maxElem->next;

				//If the extracted term has another product in the stream,
				//update the product and push onto the oldMaxElem chain
				if(maxElem->b != lastB) {
					++(maxElem->b);
					maxElem->next = oldMaxElem;
					oldMaxElem = maxElem;
				} else {
					//we are done with the maxElem ProductHeapChain
					maxElem->next = NULL;
					prodheapFreeChain_AAZ(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;

			nextDegs = prodheapPeek_AAZ(h);
		}

		//Commit new term to the product.
		if (mpz_sgn(cElems[k].coef) != 0) {
			++k;
		} else {
			//reset accumulator variables and do not increment k.
			mpz_clear(cElems[k].coef);
		}

		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			prodheapInsert_AAZ(h, maxElem, aElems[maxElem->a_i].degs + bElems[maxElem->b].degs);
			maxElem = nextMaxElem;
		}

		if (k >= allocC) {
			allocC <<= 1;
			resizePolynomial_AAZ(c, allocC);
			// fprintf(stderr, "\n\nRESIZING\n\n");
			cElems = c->elems;
		}
	}

	// mpz_clear(ccoef);
	prodheapFree_AAZ(h);

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);

	return c;
}

/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure.
 * It is assumed that both a and b have compatible exponent vectors.
 *
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
void multiplyPolynomialsPreAlloc_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, AltArrZ_t** cc) {
	if (cc == NULL) {
		return;
	}

	if (a == NULL || a->size == 0 || b == NULL || b->size == 0) {
		if (*cc == NULL) {
			return;
		}
		AltArrZ_t* c = *cc;
		for (polysize_t i = 0; i < c->size; ++i) {
			mpz_clear(c->elems[i].coef);
		}
		c->size = 0;
		return;
	}

	if (a->nvar != b->nvar) {
		fprintf(stderr, "Polynomials are incompatible in multiplyPolynomialPreAlloc\n" );
		free((int*) -1);
		exit(1);
	}

	int nvar = a->nvar;
	if (a->unpacked || b->unpacked) {
		multiplyPolynomialsPreAlloc_AAZ_unpk(a,b,cc);
		return;
	}

	degrees_t aMax = calculateMaxDegs_AAZ(a);
	degrees_t bMax = calculateMaxDegs_AAZ(b);

	if (!monomialMultFitsPacked(aMax, bMax, nvar)) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) a);
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) b);
		multiplyPolynomialsPreAlloc_AAZ_unpk(a,b,cc);
		packExponentVectors_AAZ_inp((AltArrZ_t*) a);
		packExponentVectors_AAZ_inp((AltArrZ_t*) b);
		return;
	}

	// reorder to obtain smaller as a.
	if (b->size < a->size) {
		const AltArrZ_t* temp = a;
		a = b;
		b = temp;
	}

	ProductHeap_AAZ* h = prodheapInit_AAZ(a,b,nvar);
	// mpz_t ccoef;
	// mpz_init(ccoef);

	register unsigned long long int allocC = AA_SIZE(a)*AA_SIZE(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;

	AltArrZ_t* c = *cc;
	if (c == NULL) {
		c = makePolynomial_AAZ(allocC, nvar);
		*cc = c;
	} else {
		if (c->unpacked) {
			clearAndPackExponentVectors_AAZ_inp(c);
		}
		allocC = c->alloc;
	}

	c->nvar = nvar; //since we're packed we can just do this

	//k is c, i is a, j is b.
	register polysize_t k = 0;
	register polysize_t cSize = c->size;
	// register int i = 0;
	// register int j = 0;

	AAZElem_t* __restrict__ cElems = c->elems;
	AAZElem_t* __restrict__ bElems = b->elems;

	AAZElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1;
	register int firstB = 0;

	ProductHeapChain_AAZ* maxElem = NULL;
	ProductHeapChain_AAZ* nextMaxElem = NULL;
	degrees_t* nextDegs;
	while ( (nextDegs = prodheapPeek_AAZ(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		cElems[k].degs = *nextDegs;
		if (k >= cSize) {
			mpz_init(cElems[k].coef);
		} else {
			mpz_set_ui(cElems[k].coef, 0ul);
		}

		while (nextDegs != NULL && isEqualExponentVectors(cElems[k].degs, *nextDegs)) {
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AAZ* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AAZ(h);
			while (maxElem != NULL) {

				mpz_addmul(cElems[k].coef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				// mpz_mul(ccoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				// mpz_add(cElems[k].coef, cElems[k].coef, ccoef);

#if PRODHEAP_LAZY_INIT
				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA) {
					oldMaxElem = prodheapMakeChain_AAZ((maxElem->a_i)+1, firstB, oldMaxElem);
				}
#endif

				//cache next before freeing or overwriting
				nextMaxElem = maxElem->next;

				//If the extracted term has another product in the stream,
				//update the product and push onto the oldMaxElem chain
				if(maxElem->b != lastB) {
					++(maxElem->b);
					maxElem->next = oldMaxElem;
					oldMaxElem = maxElem;
				} else {
					//we are done with the maxElem ProductHeapChain
					maxElem->next = NULL;
					prodheapFreeChain_AAZ(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;

			nextDegs = prodheapPeek_AAZ(h);
		}

		//Commit new term to the product.
		if (mpz_sgn(cElems[k].coef) != 0) {
			++k;
		} else {
			//reset accumulator variables and do not increment k.
			if ( k >= cSize) {
				mpz_clear(cElems[k].coef);
			} else {
				mpz_set_ui(cElems[k].coef, 0ul);
			}
		}

		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			prodheapInsert_AAZ(h, maxElem, aElems[maxElem->a_i].degs + bElems[maxElem->b].degs);
			maxElem = nextMaxElem;
		}

		if (k >= allocC) {
			allocC <<= 1;
			resizePolynomial_AAZ(c, allocC);
			// fprintf(stderr, "\n\nRESIZING\n\n");
			cElems = c->elems;
		}
	}

	// mpz_clear(ccoef);
	prodheapFree_AAZ(h);

	AA_SIZE(c) = k;
	for ( ; k < cSize; ++k) {
		mpz_clear(cElems[k].coef);
	}
	resizePolynomial_AAZ(c, k);
}


AltArrZ_t* multiplyPolynomials_AAZ_inp(AltArrZ_t* a, const AltArrZ_t* b, int nvar) {
	//TODO actually do this in-place.
	AltArrZ_t* prod = multiplyPolynomials_AAZ(a,b,nvar);
	freePolynomial_AAZ(a);
	return prod;
}

AltArrZ_t* multiplyManyPolynomials_AAZ(AltArrZ_t const*const* polys, int n) {
	if (n == 0 || polys == NULL) {
		return NULL;
	}
	if (n < 2) {
		return deepCopyPolynomial_AAZ(polys[0]);
	}

	polysize_t maxSize = 1;
	for (int i = 0; i < n; ++i) {
		maxSize *= polys[i]->size;
	}

	maxSize = maxSize > 60000000 ? 60000000 : maxSize;
	AltArrZ_t* prod = makePolynomial_AAZ(maxSize, 1);
	multiplyPolynomialsPreAlloc_AAZ(polys[0], polys[1], &prod);
	AltArrZ_t* prod2;
	for (int i = 2; i < n; ++i) {
		prod2 = deepCopyPolynomial_AAZ(prod);
		multiplyPolynomialsPreAlloc_AAZ(polys[i], prod2, &prod);
		freePolynomial_AAZ(prod2);
	}

	return prod;
}

AltArrZ_t* multiplyAllButOnePolynomials_AAZ(AltArrZ_t const*const* polys, int n, int idx) {
	if (n == 0 || polys == NULL) {
		return NULL;
	}
	if (n <= 2) {
		if (idx == 0) {
			return deepCopyPolynomial_AAZ(polys[1]);
		}
		return deepCopyPolynomial_AAZ(polys[0]);
	}

	polysize_t maxSize = 1;
	for (int i = 0; i < n; ++i) {
		maxSize *= polys[i]->size;
	}

	maxSize = maxSize > 60000000 ? 60000000 : maxSize;

//1: Naive iteartive way
	AltArrZ_t* prod = makePolynomial_AAZ(maxSize, polys[0]->nvar);
	AltArrZ_t* prod2;
	if (idx < 2) {
		if (idx == 0) {
			multiplyPolynomialsPreAlloc_AAZ(polys[1], polys[2], &prod);
		} else {
			multiplyPolynomialsPreAlloc_AAZ(polys[0], polys[2], &prod);
		}
		for (int i = 3; i < n; ++i) {
			prod2 = deepCopyPolynomial_AAZ(prod);
			multiplyPolynomialsPreAlloc_AAZ(polys[i], prod2, &prod);
			freePolynomial_AAZ(prod2);
		}
	} else {
		multiplyPolynomialsPreAlloc_AAZ(polys[0], polys[1], &prod);
		for (int i = 2; i < n; ++i) {
			if ( i == idx) {
				continue;
			}
			prod2 = deepCopyPolynomial_AAZ(prod);
			multiplyPolynomialsPreAlloc_AAZ(polys[i], prod2, &prod);
			freePolynomial_AAZ(prod2);
		}
	}

	return prod;
}

void multiplyManyPolynomialsPreAlloc_AAZ( AltArrZ_t const*const*  polys, int n, AltArrZ_t** retProd) {
	if (n < 2 || polys == NULL || retProd == NULL) {
		return;
	}

	polysize_t maxSize = 1;
	for (int i = 0; i < n; ++i) {
		maxSize *= polys[i]->size;
	}

	AltArrZ_t* prod = *retProd;
	multiplyPolynomialsPreAlloc_AAZ(polys[0], polys[1], &prod);
	AltArrZ_t* prod2;
	for (int i = 2; i < n; ++i) {
		prod2 = deepCopyPolynomial_AAZ(prod);
		multiplyPolynomialsPreAlloc_AAZ(polys[i], prod2, &prod);
		freePolynomial_AAZ(prod2);
	}

}



void multiplyPolynomialAtIdxByXn_AAZ_inp (AltArrZ_t* aa, int idx, int n, int nvar)
{
	if (aa == NULL) {
		return;
	}

	if (idx < 0 || idx >= nvar) {
		return;
	}

	if (aa->unpacked) {
		// TODO
		// multiplyPolynomialByXnAtIdx_AAZ_inp_unpk(aa, idx, degsList, nvar);
		fprintf(stderr, "SMQP Error, multiplyPolynomialByXnAtIdx_AAZ_inp not yet implemented with exponent unpacking!\n");
		exit (EXIT_FAILURE);
		// return;
	}

	const int* sizes = getExpOffsetArray (nvar);
	const degrees_t* maxDegs = getMaxExpArray(nvar);

	if (n > maxDegs[idx]) {
		fprintf(stderr, "SMQP Error, multiplyPolynomialByXnAtIdx_AAZ_inp doesn't work when n > maxDegs[idx] !\n");
		exit (EXIT_FAILURE);
	}

	for (int i = 0; i < AA_SIZE(aa); i++) {

		aa->elems[i].degs += ((degrees_t) n << sizes[idx]);
	}

}

/*****************
 * Polynomial exponentiation
 *****************/


AltArrZ_t* exponentiatePoly_AAZ(AltArrZ_t* a, unsigned int n, int nvar) {
	if (n == 0) {
		AltArrZ_t* ret = makePolynomial_AAZ(1, nvar);
		mpz_init(ret->elems->coef);
		mpz_set_ui(ret->elems->coef, 1ul);
		ret->elems->degs = 0;
		ret->size = 1;
		return ret;
	} else if (n == 1) {
		return deepCopyPolynomial_AAZ(a);
	}

	AltArrZ_t* r = NULL;
	AltArrZ_t* b = deepCopyPolynomial_AAZ(a);
	while (n > 1) {
		if (n & 1) {
			r = (r == NULL) ? deepCopyPolynomial_AAZ(b) : multiplyPolynomials_AAZ_inp(r, b, nvar);
		}
		b = multiplyPolynomials_AAZ_inp(b, b, nvar);
		n >>= 1;
	}
	r = (r == NULL) ? deepCopyPolynomial_AAZ(b) : multiplyPolynomials_AAZ_inp(r, b, nvar);

	freePolynomial_AAZ(b);
	return r;
}



/*****************
 * Polynomial division
 *****************/

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
void divisionGetNextTerm_AAZ(ProductHeap_AAZ* h, const AAZElem_t* __restrict__ aElems, const AAZElem_t* __restrict__ bElems, mpz_t* retCoef) {
	if (h->heapSize == 0) {
		return;
	}

	int lastB = h->lastB;

	ProductHeapChain_AAZ* insertChain = NULL;
	ProductHeapChain_AAZ* maxElem, *nextMaxElem;

	// mpz_t prodCoef;
	// mpz_init(prodCoef);
	degrees_t* nextDegs = prodheapPeek_AAZ(h);
	register degrees_t maxDegs = *nextDegs;

	while ( nextDegs != NULL && isEqualExponentVectors(maxDegs, *nextDegs)) {
		maxElem = prodheapRemoveMax_AAZ(h);

		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;

			mpz_addmul(*retCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
			// mpz_mul(prodCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
			// mpz_add(*retCoef, *retCoef, prodCoef);
			if (maxElem->b != lastB) {
				++(maxElem->b);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				maxElem->next = NULL;
				prodheapFreeChain_AAZ(maxElem);
			}

			maxElem = nextMaxElem;
		}

		nextDegs = prodheapPeek_AAZ(h);
	}

	while(insertChain != NULL) {
		maxElem = insertChain->next;
		insertChain->next = NULL;
		prodheapInsert_AAZ(h,insertChain, aElems[insertChain->a_i].degs + bElems[insertChain->b].degs);
		insertChain = maxElem;
	}

	// mpz_clear(prodCoef);
}

void divideBySingleTerm_AAZ(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar) {
	if (b == NULL || b->size == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}
	if (nvar == 0) {
		if(mpz_divisible_p(c->elems->coef, b->elems->coef)) {
			*res_a = makeConstPolynomial_AAZ(1, 0, c->elems->coef);
			mpz_divexact((*res_a)->elems->coef, (*res_a)->elems->coef, b->elems->coef);
			*res_r = NULL;
		} else {
			*res_a = NULL;
			*res_r = makeConstPolynomial_AAZ(1, 0, c->elems->coef);
		}
		return;
	}

	if (c->unpacked || b->unpacked) {
		divideBySingleTerm_AAZ_unpk(c, b, res_a, res_r, nvar);
		return;
	}


	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxSize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	AAZElem_t* curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar) && mpz_divisible_p(k->coef, bElems->coef)) {
			mpz_divexact(curA->coef, k->coef, bElems->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			mpz_set(curR->coef, k->coef);
			curR->degs = k->degs;
			++j;
			++(curR);
			mpz_init(curR->coef);
		}
		++k;
	}

	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxSize;
	r->alloc = maxSize;

	*res_a = a;
	*res_r = r;
}

/**
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials_AAZ(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, register int nvar) {

	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AAZ(0, nvar);
		*res_r = NULL;
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		divideBySingleTerm_AAZ(c, b, res_a, res_r, nvar);
		return;
	}

	if (c->unpacked || b->unpacked) {
		dividePolynomials_AAZ_unpk(c, b, res_a, res_r, nvar);
		return;
	}

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	AAZElem_t* __restrict__ curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);


	register degrees_t beta = b->elems->degs;


	//init a with lt(c)/lt(b);
	DivTest_ptr divTest = getMonomialDivideTestFuncPtr(nvar);
	while (k != lenK && !( (*divTest)(k->degs, beta) && mpz_divisible_p(k->coef, b->elems->coef)) ) {
		mpz_set(curR->coef, k->coef);
		curR->degs = k->degs;
		++j;
		if (j >= maxRSize) {
			maxRSize <<= 1;
			r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
			curR = r->elems + j - 1;
		}
		++(curR);
		mpz_init(curR->coef);
		++k;
	}

	if (k == lenK) {
		//no division to do at all!
		mpz_clear(curA->coef);
		mpz_clear(curR->coef);

		AA_SIZE(a) = i;
		AA_SIZE(r) = j;
		a->alloc = maxASize;
		r->alloc = maxRSize;

		*res_a = a;
		*res_r = r;
		return;
	}

	subtractExponentVectors(k->degs, beta, curA->degs);
	//assuming b is monic!
	// mpz_set(curA->coef, k->coef);
	mpz_divexact(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps;
	register cmpExp_t cmp;
	while(delta != NULL || k != lenK) {

		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AAZ(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AAZ(h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if ((*divTest)(eps, beta) && mpz_divisible_p(curA->coef, b->elems->coef)) {
			subtractExponentVectors(eps, beta, curA->degs);
			mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {

			//swap here so that curA becomes 0.
			mpz_swap(curR->coef, curA->coef);
			curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
				curR = r->elems + j - 1;
			}
			++(curR);
			mpz_init(curR->coef);
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	*res_a = a;
	*res_r = r;
}

void exactDividePolynomials_AAZ (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, register int nvar)
{
	if (b == NULL || AA_SIZE(b) == 0) {
	//division by zero
		fprintf (stderr, "Division by zero! Exiting...");
		exit (EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
	//c is zero
		*res_a = makePolynomial_AAZ (0, nvar);
		return;
	}

	if (c->unpacked || b->unpacked) {
		exactDividePolynomials_AAZ_unpk(c, b, res_a, nvar);
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		AltArrZ_t* res_r;
		divideBySingleTerm_AAZ (c, b, res_a, &res_r, nvar);
	//fprintf (stderr, "r := \n" );
	//printAAZ (res_r);
		freePolynomial_AAZ (res_r);
		return;
	}

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	// register int maxRSize = maxASize;
	register int i = 0;
	// register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ (maxASize, nvar);
	// AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	// AAZElem_t* __restrict__ curR = r->elems;
	mpz_init(curA->coef);
	// mpz_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	DivTest_ptr divTest = getMonomialDivideTestFuncPtr(nvar);
	// while (k != lenK) { // && !( (*divTest)(k->degs, beta) && mpz_divisible_p(k->coef, b->elems->coef))
	// mpz_set(curR->coef, k->coef);
	// curR->degs = k->degs;
	// ++j;
	// if (j >= maxRSize) {
	// maxRSize <<= 1;
	//  r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
	//   curR = r->elems + j - 1;
	// }
	// ++(curR);
	// mpz_init(curR->coef);
	// ++k;
	// }

	if (k == lenK) {
	//no division to do at all!
		mpz_clear (curA->coef);
	// mpz_clear(curR->coef);

		AA_SIZE(a) = i;
	// AA_SIZE(r) = j;
		a->alloc = maxASize;
	// r->alloc = maxRSize;

		*res_a = a;
	// *res_r = r;
		return;
	}

	subtractExponentVectors (k->degs, beta, curA->degs);
	//assuming b is monic!
	// mpz_set(curA->coef, k->coef);
	mpz_divexact (curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ (nvar);
	prodheapResize_AAZ (h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ (h, prodheapMakeChain_AAZ (0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	mpz_init (curA->coef);

	degrees_t* delta = prodheapPeek_AAZ (h);
	register degrees_t eps;
	register cmpExp_t cmp;
	while (delta != NULL || k != lenK) {

		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors (*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ (h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn (curA->coef) == 0) {
		//in this case, the term with degree delta ended up
		//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ (h);
				continue;
			} else {
				mpz_neg (curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ (h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn (curA->coef) == 0) {
				delta = prodheapPeek_AAZ (h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub (curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn (curA->coef) == 0) {
					delta = prodheapPeek_AAZ (h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set (curA->coef,k->coef);
			++k;
		}

		if ((*divTest)(eps, beta) && mpz_divisible_p(curA->coef, b->elems->coef)) {
			subtractExponentVectors(eps, beta, curA->degs);
			mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		}
		// else {

			//swap here so that curA becomes 0.
			//mpz_swap(curR->coef, curA->coef);
			//curR->degs = eps;
			//++j;
			/* if (j >= maxRSize) { */
			/* 	maxRSize <<= 1; */
			/* 	r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t)); */
			/* 	curR = r->elems + j - 1; 	 */
			/* } */
			// ++(curR);
			// mpz_init(curR->coef);
		// }

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	//mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	//AA_SIZE(r) = j;
	a->alloc = maxASize;
	//r->alloc = maxRSize;

	*res_a = a;
	//*res_r = r;
	return;
}

void univariatePseudoDivideBySingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy) {
	if (b == NULL || b->size == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (c->unpacked || b->unpacked) {
		univariatePseudoDivideBySingleTerm_AAZ_unpk(c, b, res_a, res_r, e, lazy);
		return;
	}

	int nvar = 1;
	c = deepCopyPolynomial_AAZ(c);

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxSize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	AAZElem_t* curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	int multSteps = 0;

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar)) {
			if (mpz_divisible_p(k->coef, bElems->coef) == 0) {
				multiplyByInteger_AAZ_inp(c, bElems->coef);
				multiplyByInteger_AAZ_inp(a, bElems->coef);
				++multSteps;
			}

			mpz_divexact(curA->coef, k->coef, bElems->coef);
			// mpz_set(curA->coef, k->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			++(a->size);
			mpz_init(curA->coef);
		} else {
			mpz_set(curR->coef, k->coef);
			curR->degs = k->degs;
			++j;
			++(curR);
			mpz_init(curR->coef);
		}
		++k;
	}

	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxSize;
	r->alloc = maxSize;

	if (!lazy) {
		degrees_t d = c->elems->degs - b->elems->degs + 1 - multSteps;
		multSteps += d;

		mpz_t bPow;
		mpz_init(bPow);
		mpz_set(bPow, b->elems->coef);
		mpz_pow_ui(bPow, bPow, d);
		multiplyByInteger_AAZ_inp(a, bPow);
		multiplyByInteger_AAZ_inp(r, bPow);
		mpz_clear(bPow);

		// for (int j = 0; j < d; ++j) {
		// 	multiplyByInteger_AAZ_inp(a, b->elems->coef);
		// 	multiplyByInteger_AAZ_inp(r, b->elems->coef);
		// }
	}

	*res_a = a;
	*res_r = r;

	if (e != NULL) {
		*e = multSteps;
	}

	freePolynomial_AAZ(c);
}

void univariatePseudoDividePolynomials_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	int nvar = 1;
	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AAZ(0, nvar);
		*res_r = NULL;
		return;
	}

	if (c->unpacked || b->unpacked) {
		univariatePseudoDividePolynomials_AAZ_unpk(c, b, res_a, res_r, e, lazy);
		return;
	}

	if (isLessExponentVectors(c->elems->degs, b->elems->degs)) {
		*res_r = deepCopyPolynomial_AAZ(c);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		univariatePseudoDivideBySingleTerm_AAZ(c, b, res_a, res_r, e, lazy);
		return;
	}

	c = deepCopyPolynomial_AAZ(c);

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	AAZElem_t* __restrict__ curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	int multSteps = 0;
	subtractExponentVectors(k->degs, beta, curA->degs);
	if (mpz_divisible_p(k->coef, b->elems->coef) == 0) {
		//set the quotient coefficient before updating dividend!
		mpz_set(curA->coef, k->coef);
		multiplyByInteger_AAZ_inp(c, b->elems->coef);
		++multSteps;
	} else {
		mpz_divexact(curA->coef, k->coef, b->elems->coef);
	}

	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	a->size = 1;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps;
	register cmpExp_t cmp;
	while(delta != NULL || k != lenK) {
		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AAZ(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AAZ(h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (compareExponentVectors(eps, beta) >= 0) {
			if (mpz_divisible_p(curA->coef, b->elems->coef) == 0) {
				multiplyByInteger_AAZ_inp(c, b->elems->coef);
				multiplyByInteger_AAZ_inp(a, b->elems->coef);
				++multSteps;

			} else {
				mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			}

			subtractExponentVectors(eps, beta, curA->degs);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
			++i;
			++(curA);
			++(a->size);
			mpz_init(curA->coef);
		} else {
			//swap here so that curA becomes 0.
			mpz_swap(curR->coef, curA->coef);
			curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
				curR = r->elems + j - 1;
			}
			++(curR);
			mpz_init(curR->coef);
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	if (!lazy) {
		degrees_t d = c->elems->degs - b->elems->degs + 1 - multSteps;
		multSteps += d;

		mpz_t bPow;
		mpz_init(bPow);
		mpz_set(bPow, b->elems->coef);
		mpz_pow_ui(bPow, bPow, d);
		multiplyByInteger_AAZ_inp(a, bPow);
		multiplyByInteger_AAZ_inp(r, bPow);
		mpz_clear(bPow);

		// for (int j = 0; j < d; ++j) {
			// multiplyByInteger_AAZ_inp(a, b->elems->coef);
			// multiplyByInteger_AAZ_inp(r, b->elems->coef);
		// }
	}

	*res_a = a;
	*res_r = r;

	if (e != NULL) {
		*e = multSteps;
	}

	freePolynomial_AAZ(c);
}


int divideTestSingleTerm_AAZ(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {
	if (b == NULL || b->size == 0) {
		return 0;
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		return 1;
	}

	if (c->unpacked || b->unpacked) {
		return divideTestSingleTerm_AAZ_unpk(c, b, res_a, nvar);
	}

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	mpz_init(curA->coef);

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar) && mpz_divisible_p(k->coef, bElems->coef)) {
			mpz_divexact(curA->coef, k->coef, bElems->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			mpz_clear(curA->coef);
			a->size = i;
			freePolynomial_AAZ(a);
			return 0;
		}
		++k;
	}

	mpz_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxSize;

	if (res_a != NULL) {
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}

	return 1;
}

int divideTest_AAZ(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {

	if (b == NULL || AA_SIZE(b) == 0) {
		return 0;
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		*res_a = makePolynomial_AAZ(0, nvar);
		return 1;
	}

	if (c->unpacked || b->unpacked) {
		return divideTest_AAZ_unpk(c, b, res_a, nvar);
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		return	divideTestSingleTerm_AAZ(c, b, res_a, nvar);
	}

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register degrees_t beta = b->elems->degs;

	if (!monomialDivideTest(k->degs, beta, nvar) || !mpz_divisible_p(k->coef, b->elems->coef)) {
		return 0;
	}

	register int maxASize = AA_SIZE(b) > 10 ? AA_SIZE(b) + 1 : 10;
	register int i = 0;
	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	mpz_init(curA->coef);

	//init a with lt(c)/lt(b);
	subtractExponentVectors(k->degs, beta, curA->degs);
	mpz_divexact(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps;
	register cmpExp_t cmp;
	while(k != lenK || delta != NULL) {

		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AAZ(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AAZ(h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest(eps, beta, nvar) && mpz_divisible_p(curA->coef, b->elems->coef)) {
			subtractExponentVectors(eps, beta, curA->degs);
			mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);

			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {

			//divide test fails
			mpz_clear(curA->coef);
			prodheapFree_AAZ(h);

			a->size = i;
			freePolynomial_AAZ(a);
			return 0;
		}

		delta = prodheapPeek_AAZ(h);
	}

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxASize;

	prodheapFree_AAZ(h);

	if (res_a != NULL) {
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
	return 1;
}

void divideByLeadingTerms_AAZ (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar)
{
	if (b == NULL || b->size == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (c->unpacked || b->unpacked) {
		divideByLeadingTerms_AAZ_unpk (c, b, res_a, res_r, nvar);
		return;
	}

	AltArrZ_t* a;
	AltArrZ_t* r;

	/* while (k != lenK) { */
	if (monomialDivideTest (c->elems->degs, b->elems->degs, nvar)) {
		// a:
		a = makePolynomial_AAZ (1, nvar);
		mpz_init (a->elems->coef);
		mpz_div (a->elems->coef, c->elems->coef, b->elems->coef);
		subtractExponentVectors (c->elems->degs, b->elems->degs, a->elems->degs);
		AA_SIZE (a) = 1;

		// r:
		r = NULL;
	} else {
		// a:
		a = makePolynomial_AAZ (1, nvar);
		mpz_init (a->elems->coef);
		mpz_set_si (a->elems->coef, 1l);
		a->elems->degs = 0l;
		AA_SIZE (a) = 1;

		// r:
		r = makePolynomial_AAZ (1, nvar);
		mpz_init (r->elems->coef);
		mpz_set(r->elems->coef, c->elems->coef);
		r->elems->degs = c->elems->degs;
		AA_SIZE (r) = 1;
	}

	*res_a = a;
	*res_r = r;
}



/*****************
 * Derivative / Integral
 *****************/

AltArrZ_t* derivative_AAZ(const AltArrZ_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx < 0 || idx >= aa->nvar) {
		return NULL;
	}

	if (aa->unpacked) {
		return derivative_AAZ_unpk(aa, idx, k);
	}

	AltArrZ_t* ret = makePolynomial_AAZ(aa->size, aa->nvar);

	const int* sizes = getExpOffsetArray(aa->nvar);
	const degrees_t* masks = getExpMaskArray(aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAZElem_t* __restrict__ elems = aa->elems;
	AAZElem_t* __restrict__ retElems = ret->elems;
	degrees_t deg;
	mpz_t mpzDeg;
	mpz_init(mpzDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i) {
		if (!(elems[i].degs & masks[idx])) {
			continue;
		}
		deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);
		if (deg < k) {
			continue;
		}

		retElems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpz_init(retElems[insertIdx].coef);
		mpz_set(retElems[insertIdx].coef, elems[i].coef);

		mpz_set_ui(mpzDeg, deg);
		for (int j = 0; j < k; ++j) {
			mpz_mul(retElems[insertIdx].coef, retElems[insertIdx].coef, mpzDeg);
			--deg;
			mpz_sub(mpzDeg, mpzDeg, mpzOne);
		}
		retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpz_clear(mpzDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	if (insertIdx == 0) {
		freePolynomial_AAZ(ret);
		return NULL;
	}
	return ret;
}

void derivative_AAZ_inp(AltArrZ_t** aa_p, int idx, int k) {
	if (aa_p == NULL) {
		return;
	}

	if (isZero_AAZ(*aa_p)) {
		freePolynomial_AAZ(*aa_p);
		*aa_p = NULL;
		return;
	}

	if (idx < 0 || idx >= (*aa_p)->nvar) {
		freePolynomial_AAZ(*aa_p);
		*aa_p = NULL;
		return;
	}

	if ((*aa_p)->unpacked) {
		derivative_AAZ_inp_unpk(aa_p, idx, k);
		return;
	}

	AltArrZ_t* aa = *aa_p;
	const int* sizes = getExpOffsetArray(aa->nvar);
	const unsigned long long int* masks = getExpMaskArray(aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAZElem_t* elems = aa->elems;
	degrees_t deg;
	for (int i = 0; i < size; ++i) {
		if (!(elems[i].degs & masks[idx])) {
			continue;
		}
		deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);
		if (deg < k) {
			continue;
		}

		elems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpz_set(elems[insertIdx].coef, elems[i].coef);

		for (int j = 0; j < k; ++j) {
			mpz_mul_ui(elems[insertIdx].coef, elems[insertIdx].coef, (unsigned long) deg);
			--deg;
		}
		elems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}
	for (int i = insertIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
		elems[i].degs = 0l;
	}

	aa->size = insertIdx;
	if (insertIdx == 0) {
		freePolynomial_AAZ(aa);
		*aa_p = NULL;
	}
}



/**
 * Integrate with respect to a variable that does not currently exist in aa.
 * It becomes the main variable (that is, to the left).
 *
 */
AltArr_t* integrateExpand_AAZ(AltArrZ_t* aa, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	AltArr_t* ret = deepCopyPolynomial_AAFromAAZ(aa);
	AltArr_t* intRet = integrateExpand_AA(ret, k);
	freePolynomial_AA(ret);
	return intRet;
}


AltArr_t* integral_AAZ(AltArrZ_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	AltArr_t* ret = deepCopyPolynomial_AAFromAAZ(aa);
	AltArr_t* intRet = integral_AA(ret, idx, k);
	freePolynomial_AA(ret);
	return intRet;
}



/*****************
 * Content, PrimitivePart, etc.
 *****************/

void integralContent_AAZ(const AltArrZ_t* aa, mpz_t ret) {
	if (aa == NULL || aa->size <= 0) {
		mpz_set_ui(ret, 1ul);
		return;
	}

	const AAZElem_t* elems = aa->elems;
	int size = aa->size;
	mpz_set_ui(ret, 1ul);

	mpz_t one;
	mpz_init(one);
	mpz_set_si(one, 1l);

	mpz_t cont;
	mpz_init(cont);
	mpz_abs(cont, elems->coef);

	for (int i = 1; i < size; ++i) {
		if (mpz_cmp(one, cont) == 0) {
			break;
		}
		mpz_gcd(cont, cont, elems[i].coef);
	}

	if (mpz_sgn(elems->coef) < 0 && mpz_sgn(cont) > 0) {
		mpz_neg(cont, cont);
	}

	mpz_set(ret, cont);
	mpz_clear(cont);
	mpz_clear(one);
}

AltArrZ_t* primitivePart_AAZ(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	mpz_t content;
	mpz_init(content);

	integralContent_AAZ(aa, content);

	AltArrZ_t* res = deepCopyPolynomial_AAZ(aa);
	if (mpz_cmp_si(content, 1l) == 0) {
		return res;
	}

	AAZElem_t* elems = res->elems;
	int size = res->size;
	for (int i = 0; i < size; ++i) {
		mpz_divexact(elems[i].coef, elems[i].coef, content);
	}

	mpz_clear(content);
	return res;
}

AltArrZ_t* primitivePartAndContent_AAZ(const AltArrZ_t* aa, mpz_t cont) {
	if (aa == NULL || aa->size == 0) {
		mpz_set_si(cont, 0l);
		return NULL;
	}

	integralContent_AAZ(aa, cont);

	AltArrZ_t* res = deepCopyPolynomial_AAZ(aa);
	if (mpz_cmp_si(cont, 1l) == 0) {
		return res;
	}

	AAZElem_t* elems = res->elems;
	int size = res->size;
	for (int i = 0; i < size; ++i) {
		mpz_divexact(elems[i].coef, elems[i].coef, cont);
	}

	return res;
}

AltArrZ_t* primitivePartAndContent_AAZFromAA(AltArr_t* aa, mpq_t cont) {
	if (aa == NULL || aa->size == 0) {
		mpq_set_si(cont, 0l, 1l);
		return NULL;
	}

	if (aa->unpacked) {
		return primitivePartAndContent_AAZFromAA_unpk(aa, cont);
	}

	integralContent_AA(aa, cont);

	AltArrZ_t* res = makePolynomial_AAZ(aa->size, aa->nvar);
	int size = res->size = aa->size;
	AAElem_t* elems = aa->elems;
	AAZElem_t* reselems = res->elems;

	mpq_t tmp;
	mpq_init(tmp);
	for (int i = 0; i < size; ++i) {
		mpq_div(tmp, elems[i].coef, cont);
		mpz_init_set(reselems[i].coef, mpq_numref(tmp));
		reselems[i].degs = elems[i].degs;
	}

	mpq_clear(tmp);

	return res;
}


void primitivePart_AAZ_inp(AltArrZ_t* aa) {
	if (aa == NULL || aa->size <= 0) {
		return;
	}

	mpz_t content;
	mpz_init(content);

	integralContent_AAZ(aa, content);
	if (mpz_cmp_si(content, 1l) == 0) {
		mpz_clear(content);
		return;
	}

	AAZElem_t* elems = aa->elems;
	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpz_divexact(elems[i].coef, elems[i].coef, content);
	}

	mpz_clear(content);
}

void primitivePartAndContent_AAZ_inp(AltArrZ_t* aa, mpz_t content) {
	if (aa == NULL || aa->size <= 0) {
		return;
	}

	integralContent_AAZ(aa, content);
	if (mpz_cmp_si(content, 1l) == 0) {
		return;
	}

	AAZElem_t* elems = aa->elems;
	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpz_divexact(elems[i].coef, elems[i].coef, content);
	}
}

//Get the content of aa w.r.t. the variables marked as active
//by active[i] > 0;
//the content is one such that the leading coefficient
//of the corresponding primitive part is positive.
AltArrZ_t* contentInVars_AAZ(const AltArrZ_t* aa_in, const int* active) {
	if (isZero_AAZ(aa_in)) {
		return NULL;
	}

	if (aa_in->unpacked) {
		return contentInVars_AAZ_unpk(aa_in, active);
	}

	int nvar = aa_in->nvar;

	//need to move all active variables to the front so lex ordering
	//will make the below algorithm work.
	int nactive = 0;
	int maxActiveIdx = 0;
	for (int i = 0; i < nvar; ++i) {
		if (active[i]) {
			nactive++;
			maxActiveIdx = i;
		}
	}

	if (nactive == nvar) {
		mpz_t cont;
		mpz_init(cont);
		integralContent_AAZ(aa_in, cont);
		AltArrZ_t* ret = makeConstPolynomial_AAZ(1, nvar, cont);
		mpz_clear(cont);
		return ret;
	}

	int varMap[nvar];
	int reverseMap[nvar];
	int activeIdx = 0;
	int inactiveIdx = nactive;
	for (int i = 0; i < nvar; ++i) {
		if (active[i]) {
			varMap[i] = activeIdx;
			reverseMap[activeIdx] = i;
			++activeIdx;
		} else {
			varMap[i] = inactiveIdx;
			reverseMap[inactiveIdx] = i;
			++inactiveIdx;
		}
	}

	const AltArrZ_t* aa;
	if (maxActiveIdx >= nactive) {
		aa = deepCopyPolynomial_AAZ(aa_in);
		reorderVars_AAZ((AltArrZ_t*) aa, varMap, nvar);
	} else {
		aa = aa_in;
	}


	AltArrZ_t* cont = NULL;
	const degrees_t* masks = getExpMaskArray(nvar);

	//Here we build up a mask for the first nactive vars.
	//By previous re-ordering they are guaranteed to be at the front of the exp vec.
	int j;
	degrees_t activeMask = 0ll;
	for (j = 0; j < nactive; ++j) {
		activeMask |= masks[j];
	}

	AltArrZ_t localCoef;
	localCoef.unpacked = aa->unpacked;
	localCoef.nvar = aa->nvar;

	AltArrZ_t* tmpCoef = NULL;
	AAZElem_t* elems = aa->elems;
	degrees_t curDegs = aa->elems->degs & activeMask;
	degrees_t nextDegs = 0ll;
	int startIdx = 0;
	int i, k;
	for (i = 0; i < aa->size; ++i) {
		nextDegs = aa->elems[i].degs & activeMask;
		if (nextDegs != curDegs) {
			//we have a new coef from prevCoef,...,i-1
			localCoef.size = i - startIdx;
			localCoef.alloc = localCoef.size;
			localCoef.elems = aa->elems + startIdx;
			deepCopyPolynomial_AAZ_inp(&localCoef, &tmpCoef);

			//clear out active variables from coef poly
			for (k = 0; k < tmpCoef->size; ++k) {
				tmpCoef->elems[k].degs &= (~activeMask);
			}

			if (cont == NULL) {
				cont = tmpCoef;
				tmpCoef = NULL;
			} else {
			// 	char* syms[] = {"x", "y", "z", "s", "t", "w", "v", "a", "b", "c"};
			// 	fprintf(stderr, "cont: ");
			// 	printPoly_AAZ(stderr, cont, syms, nvar);
			// 	fprintf(stderr, "\ncoef: ");
			// 	printPoly_AAZ(stderr, tmpCoef, syms, nvar);
			// 	fprintf(stderr, "\n");
				// AltArrZ_t* tmp2 = gcd_AAZ_tmp(cont, tmpCoef);
				AltArrZ_t* tmp2 = gcd_AAZ(cont, tmpCoef);
				// fprintf(stderr, "tmp2: ");
				// printPoly_AAZ(stderr, tmp2, syms, nvar);
				// fprintf(stderr, "\n\n");
				// freePolynomial_AAZ(cont);
				cont = tmp2;
				if (isOne_AAZ(cont)) {
					break;
				}
			}

			startIdx = i;
			curDegs = nextDegs;
		}
	}
	//process the last coef now
	localCoef.size = i - startIdx;
	localCoef.alloc = localCoef.size;
	localCoef.elems = aa->elems + startIdx;
	deepCopyPolynomial_AAZ_inp(&localCoef, &tmpCoef);

	//clear out active variables from coef poly
	for (k = 0; k < tmpCoef->size; ++k) {
		tmpCoef->elems[k].degs &= (~activeMask);
	}

	if (cont == NULL) {
		cont = tmpCoef;
		tmpCoef = NULL;
	} else {
		AltArrZ_t* tmp2 = gcd_AAZ(cont, tmpCoef);
		// AltArrZ_t* tmp2 = gcd_AAZ_tmp(cont, tmpCoef);
		freePolynomial_AAZ(cont);
		cont = tmp2;
	}



	freePolynomial_AAZ(tmpCoef);

	if (maxActiveIdx >= nactive) {
		freePolynomial_AAZ((AltArrZ_t*) aa);
		reorderVars_AAZ(cont, reverseMap, nvar);
	}

	return cont;
}


/**
 * Returns GCD g where s*a + t*b = g.
 * If s or t are NULL then those bezout coefficients are not returned.
 */
// AltArrZ_t* extendedEuclidean(AltArrZ_t* a, AltArrZ_t* b, AltArrZ_t** s, AltArrZ_t** t) {
// 	if (a->nvar != 1 || b->nvar != 1) {
// 		fprintf(stderr, "SMQP ERROR: Calling univariate GCD on multivariate arrays\n");
// 		exit(1);
// 	}


// 	if (a == NULL || a->size == 0) {
// 		return deepCopyPolynomial_AAZ(b);
// 	}
// 	if (b == NULL || b->size == 0){
// 		return deepCopyPolynomial_AAZ(a);
// 	}

// 	AltArrZ_t* r0 = NULL;
// 	AltArrZ_t* r1 = NULL;
// 	AltArrZ_t* r2 = NULL;

// 	mpz_t c0;
// 	mpz_t c1;
// 	mpz_init(c0);
// 	mpz_init(c1);

// 	if (isGreaterExponentVectors(a->elems->degs, b->elems->degs)) {
// 		r0 = primitivePartAndContent_AAZ(a, c0);
// 		r1 = primitivePartAndContent_AAZ(b, c1);
// 		// r0 = deepCopyPolynomial_AA(a);
// 		// r1 = deepCopyPolynomial_AA(b);
// 	} else {
// 		r0 = primitivePartAndContent_AAZ(b, c0);
// 		r1 = primitivePartAndContent_AAZ(a, c1);
// 		// r0 = deepCopyPolynomial_AA(b);
// 		// r1 = deepCopyPolynomial_AA(a);
// 	}

// 	AltArrZ_t* quo = NULL;
// 	AltArrZ_t* s0, * s1 = NULL, * s2 = NULL;
// 	AltArrZ_t* t0 = NULL, * t1, * t2 = NULL;
// 	mpz_set_si(c0, 1l);
// 	s0 = makeConstPolynomial_AAZ(1, 1, c0);
// 	t1 = makeConstPolynomial_AAZ(1, 1, c0);

// 	int i = 1;
// 	int e;
// 	while (r1 != NULL && r1->size > 0) {
// 		univariatePseudoDividePolynomials_AAZ(r0, r1, &quo, &r2, &e, 1);

// 		fprintf(stderr, "i: %d\n", ++i);

// 		freePolynomial_AAZ(r0);
// 		r0 = r1;
// 		r1 = r2;
// 		r2 = NULL;

// 		mpz_pow_ui(c0, r0->elems->coef, (unsigned long) e);
// 		// primitivePartAndContent_AAZ_inp(quo, c1);
// 		// mpz_divexact(c0, c0, c1);

// 		gmp_fprintf(stderr, "c0: %Zd,       c1: %Zd\n", c0, c1);
// 		gmp_fprintf(stderr, "e: %d,       c0: %Zd\n", e, c0);

// 		s2 = multiplyPolynomials_AAZ(s1, quo, 1);
// 		multiplyByInteger_AAZ_inp(s0, c0);
// 		s0 = subPolynomials_AAZ_inp(s0, s2, 1);
// 		freePolynomial_AAZ(s2);
// 		s2 = s0;
// 		s0 = s1;
// 		s1 = s2;
// 		s2 = NULL;


// 		t2 = multiplyPolynomials_AAZ(t1, quo, 1);
// 		multiplyByInteger_AAZ_inp(t0, c0);
// 		t0 = subPolynomials_AAZ_inp(t0, t2, 1);
// 		freePolynomial_AAZ(t2);
// 		t2 = t0;
// 		t0 = t1;
// 		t1 = t2;
// 		t2 = NULL;


// 		primitivePartAndContent_AAZ_inp(r1, c0);
// 		// mpz_divexact(c0, c0, c1);
// 		divideByIntegerExact_AAZ_inp(s1, c0);
// 		divideByIntegerExact_AAZ_inp(t1, c0);

// 		fprintf(stderr, "r: \n");
// 		printAAZ(r1);
// 		fprintf(stderr, "q: \n");
// 		printAAZ(quo);
// 		fprintf(stderr, "s0: \n");
// 		printAAZ(s0);
// 		fprintf(stderr, "s1: \n");
// 		printAAZ(s1);
// 		fprintf(stderr, "t0: \n");
// 		printAAZ(t0);
// 		fprintf(stderr, "t1: \n");
// 		printAAZ(t1);

// 		freePolynomial_AAZ(quo);
// 		quo = NULL;

// 		fprintf(stderr, "\n");
//  	}

// 	freePolynomial_AAZ(r1);
// 	freePolynomial_AAZ(r2);
// 	if (r0 != NULL && r0->size > 0 && mpz_sgn(r0->elems->coef) < 0) {
// 		negatePolynomial_AAZ(r0);
// 	}

// 	freePolynomial_AAZ(s1);
// 	freePolynomial_AAZ(t1);


// 	if (isGreaterExponentVectors(a->elems->degs, b->elems->degs)) {
// 		if (s != NULL) {
// 			*s = s0;
// 		} else {
// 			freePolynomial_AAZ(s0);
// 		}
// 		if (t != NULL) {
// 			*t = t0;
// 		} else {
// 			freePolynomial_AAZ(t0);
// 		}
// 	}
// 	else {
// 		if (s != NULL) {
// 			*s = t0;
// 		} else {
// 			freePolynomial_AAZ(t0);
// 		}
// 		if (t != NULL) {
// 			*t = s0;
// 		} else {
// 			freePolynomial_AAZ(s0);
// 		}
// 	}


// 	mpz_clear(c0);
// 	mpz_clear(c1);

// 	return r0;

// }

AltArrZ_t* univariateHeuristicGCD_AAZ(const AltArrZ_t* a, const AltArrZ_t* b) {
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ(a);
	}

	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMQP ERROR: Calling univariateHeuristicGCD on multivariate arrays\n");
		exit(1);
	}

	if (a->unpacked || b->unpacked) {
		return NULL; //just fail for now... this shouldn't happen anyway.
	}

	mpz_t conta;
	mpz_t contb;
	mpz_init(conta);
	mpz_init(contb);
	AltArrZ_t* aa = primitivePartAndContent_AAZ(a, conta);
	AltArrZ_t* bb = primitivePartAndContent_AAZ(b, contb);
	mpz_gcd(conta, conta, contb);
	mpz_clear(contb);

	AltArrZ_t* gcd = univariateHeuristicPrimitiveGCD_AAZ(aa, bb);

	multiplyByInteger_AAZ_inp(gcd, conta);
	mpz_clear(conta);
	freePolynomial_AAZ(aa);
	freePolynomial_AAZ(bb);
	return gcd;
}

void genPoly_univgcdheu_AAZ_inp(AltArrZ_t* gcd, mpz_t eps, mpz_t halfeps, mpz_t gamma) {
	AAZElem_t* elems = gcd->elems;
	int curSize = 0;
	mpz_init(elems[curSize].coef);
	for(int i = 0; mpz_sgn(gamma) != 0; ++i) {
		mpz_mod(elems[curSize].coef, gamma, eps);
		if (mpz_cmp(elems[curSize].coef, halfeps) > 0) {
			mpz_sub(elems[curSize].coef, elems[curSize].coef, eps);
		}

		mpz_sub(gamma, gamma, elems[curSize].coef);
		mpz_fdiv_q(gamma, gamma, eps);

		if (mpz_sgn(elems[curSize].coef) != 0) {
			elems[curSize].degs = i;
			++curSize;
			if (curSize >= gcd->alloc) {
				resizePolynomial_AAZ(gcd, gcd->alloc << 1);
				elems = gcd->elems;
			}
			mpz_init(elems[curSize].coef);
		}

	}
	mpz_clear(elems[curSize].coef);

	gcd->size = curSize;
	curSize >>= 1;
	AAZElem_t tmp;
	//reverse array to make degs in right order; curSize is now half the size of gcd
	for (int i = 0; i < curSize; ++i ){
		tmp = elems[i];
		elems[i] = elems[gcd->size - 1 - i];
		elems[gcd->size - 1 - i] = tmp;
	}
}

AltArrZ_t* univariateHeuristicPrimitiveGCD_AAZ(const AltArrZ_t* a, const AltArrZ_t* b) {
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ(a);
	}

	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMQP ERROR: Calling univariateHeuristicGCD on multivariate arrays\n");
		exit(1);
	}

	if (a->unpacked || b->unpacked) {
		return NULL; //just fail for now... this shouldn't happen anyway.
	}

	degrees_t maxDegs = (isLessExponentVectors(a->elems->degs, b->elems->degs) ? b->elems->degs : a->elems->degs);
	mpz_t infa, infb, eps, halfeps;
	mpz_inits(infa, infb, eps, halfeps, NULL);
	infinityNorm_AAZ(a, infa);
	infinityNorm_AAZ(b, infb);

	if (mpz_cmp(infa, infb) < 0) {
		mpz_set(infa, infb);
	}
	mpz_mul_ui(infa, infa, 2ul);
	mpz_add_ui(infa, infa, 2ul);
	mpz_swap(eps, infa);
	long epsSize = mpz_sizeinbase(eps,2);

	AltArrZ_t* gcd = NULL;
	for (int attempts = 0; attempts < 6; ++attempts) {
		//arbitrary heuristic cutoff for too expensive
		if (maxDegs * epsSize > 5000) {
			gcd = NULL;
			break;
		}
		univarEvaluate_AAZ(a, eps, infa);
		univarEvaluate_AAZ(b, eps, infb);
		mpz_gcd(infa, infa, infb);
		//infa is 'gamma', the integer gcd of the evaluation images.

		gcd = makePolynomial_AAZ(10, 1);
		mpz_fdiv_q_ui(halfeps, eps, 2ul);
		genPoly_univgcdheu_AAZ_inp(gcd, eps, halfeps, infa);
		primitivePart_AAZ_inp(gcd);


		if (divideTest_AAZ(a, gcd, NULL, 1) && divideTest_AAZ(b, gcd, NULL, 1)) {
			break;
		}

		freePolynomial_AAZ(gcd);
		gcd = NULL;

		//update eps;b
		mpz_mul_si(eps, eps, 73794l);
		mpz_fdiv_q_ui(eps, eps, 27011ul);
		epsSize = mpz_sizeinbase(eps, 2);
	}

	mpz_clears(infa, infb, eps, NULL);
	return gcd;
}

AltArrZ_t* univariatePrimitiveGCD_AAZ(const AltArrZ_t* a, const AltArrZ_t* b) {
	AltArrZ_t* gAA = univariateHeuristicPrimitiveGCD_AAZ(a,b);
	if (gAA == NULL) {
		DUZP_t* aDUZP = convertFromAltArrZ_DUZP(a);
		DUZP_t* bDUZP = convertFromAltArrZ_DUZP(b);
		DUZP_t* gDUZP = primitiveGCD_DUZP(aDUZP, bDUZP);
		gAA = convertToAltArrZ_DUZP(gDUZP);
		freePolynomial_DUZP(aDUZP);
		freePolynomial_DUZP(bDUZP);
		freePolynomial_DUZP(gDUZP);
	}
	return gAA;
}

AltArrZ_t* univariateGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b) {
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ(a);
	}

	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMQP ERROR: Calling univariate GCD on multivariate arrays\n");
		exit(1);
	}

	//conversion to DUZP handles unpacking. Return is always packed.
	mpz_t conta;
	mpz_t contb;
	mpz_init(conta);
	mpz_init(contb);
	AltArrZ_t* primA = primitivePartAndContent_AAZ(a, conta);
	AltArrZ_t* primB = primitivePartAndContent_AAZ(b, contb);
	mpz_gcd(conta, conta, contb);

	AltArrZ_t* gAA = univariatePrimitiveGCD_AAZ(primA, primB);

	multiplyByInteger_AAZ_inp(gAA, conta);

	freePolynomial_AAZ(primA);
	freePolynomial_AAZ(primB);
	mpz_clear(conta);
	mpz_clear(contb);

	return gAA;

#if 0
	// A simple method based on PRS for gcd.
	if (a->unpacked || b->unpacked) {
		return univariateGCD_AAZ_unpk(a, b);
	}

	AltArrZ_t* r0 = NULL;
	AltArrZ_t* r1 = NULL;
	AltArrZ_t* r2 = NULL;

	mpz_t c0;
	mpz_t c1;
	mpz_init(c0);
	mpz_init(c1);

	if (isGreaterExponentVectors(a->elems->degs, b->elems->degs)) {
		r0 = primitivePartAndContent_AAZ(a, c0);
		r1 = primitivePartAndContent_AAZ(b, c1);
		// r0 = deepCopyPolynomial_AA(a);
		// r1 = deepCopyPolynomial_AA(b);
	} else {
		r0 = primitivePartAndContent_AAZ(b, c0);
		r1 = primitivePartAndContent_AAZ(a, c1);
		// r0 = deepCopyPolynomial_AA(b);
		// r1 = deepCopyPolynomial_AA(a);
	}

	AltArrZ_t* quo = NULL;
	while (r1 != NULL && r1->size > 0) {
		univariatePseudoDividePolynomials_AAZ(r0, r1, &quo, &r2, NULL, 1);
		freePolynomial_AAZ(quo);
		quo = NULL;

		freePolynomial_AAZ(r0);
		r0 = r1;
		r1 = r2;
		primitivePart_AAZ_inp(r1);
		r2 = NULL;
	}

	freePolynomial_AAZ(r1);
	freePolynomial_AAZ(r2);
	if (r0 != NULL && r0->size > 0 && mpz_sgn(r0->elems->coef) < 0) {
		negatePolynomial_AAZ(r0);
	}

	mpz_clear(c0);
	mpz_clear(c1);

	return r0;
#endif
}

AltArrZ_t* commonFactor_AAZ(const AltArrZ_t* a, AltArrZ_t** factored) {
	if (a == NULL || a->size == 0) {
		return NULL;
	}

	if (a->unpacked) {
		return commonFactor_AAZ_unpk(a, factored);
	}

	AltArrZ_t* ret = makePolynomial_AAZ(1, a->nvar);
	mpz_init(ret->elems->coef);
	mpz_set_ui(ret->elems->coef, 1ul);

	if (a->nvar == 1) {
		ret->elems->degs = a->elems[a->size-1].degs;
	} else {
		degrees_t min = a->elems[a->size-1].degs;
		const degrees_t* masks = getExpMaskArray(a->nvar);
		for (int i = a->size-2; i >= 0; --i) {
			for (int j = 0; j < a->nvar; ++j) {
				if ((a->elems[i].degs & masks[j]) < (min & masks[j]) ) {
					min = min & (~masks[j]); //zero out j;
					min |= (a->elems[i].degs & masks[j]);
				}
			}
			if (isZeroExponentVector(min)) {
				break;
			}
		}
		ret->elems->degs = min;
	}
	ret->size = 1;

	if (factored != NULL) {
		AltArrZ_t* factRet = deepCopyPolynomial_AAZ(a);
		if (! isZeroExponentVector(ret->elems->degs)) {
			for (int i = 0; i < a->size; ++i) {
				factRet->elems[i].degs -= ret->elems->degs;
			}
		}
		*factored = factRet;
	}

	return ret;
}

void dehomogenizePolynomial_AAZ_inp(AltArrZ_t* aa, int varIdx) {
	if (isZero_AAZ(aa)) {
		return;
	}

	shrinkNumVarsAtIdx_AAZ(aa, varIdx);
	canonicalizePolynomial_AAZ(aa);
}

//Add in a variable at index varIdx, pushing the variable at that index
//and the following variables to the right in the exponent vector
void homogenizePolynomial_AAZ_inp(AltArrZ_t* aa, int varIdx) {

	if (isZero_AAZ(aa)) {
		if (aa != NULL) {
			expandNumVars_AAZ(aa, aa->nvar+1);
		}
		return;
	}
	if (isConstant_AAZ(aa)) {
		expandNumVars_AAZ(aa, aa->nvar+1);
		return;
	}


	if (aa->unpacked || aa->nvar+1 > MAX_NVAR_PACKED) {
		homogenizePolynomial_AAZ_inp_unpk(aa, varIdx);
		return;
	}

	int size = aa->size;
	degree_t tdeg = totalDegree_AAZ(aa);

	//expand exponent vectors out of place, in case we have to unpack
	AAZElem_t* tmpElems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*size);
	memcpy(tmpElems, aa->elems, sizeof(AAZElem_t)*size);

	int nvar = aa->nvar;
	const degrees_t* __restrict__ oldMasks = getExpMaskArray(nvar);
	const int* __restrict__ oldSizes = getExpOffsetArray(nvar);
	const int* __restrict__ newSizes = getExpOffsetArray(nvar+1);
	const degrees_t* __restrict__ maxExps = getMaxExpArray(nvar+1);

	degrees_t degs;
	degrees_t curDeg;
	degrees_t curTDeg;
	int j;
	for (int i = 0; i < aa->size; ++i) {
		degs = tmpElems[i].degs;
		tmpElems[i].degs = 0;
		curTDeg = 0;

		//re-pack the variables [0,varIdx) to the same indices
		for (j = 0; j < varIdx; ++j) {
			curDeg = GET_NTH_EXP(degs, oldMasks[j], oldSizes[j]);
			if (curDeg > maxExps[j]) {
				free(tmpElems);
				unpackExponentVectors_AAZ_inp(aa);
				homogenizePolynomial_AAZ_inp_unpk(aa, varIdx);
				return;
			}
			curTDeg += curDeg;
			tmpElems[i].degs |= (curDeg << newSizes[j]);
		}

		//re-pack the variables [varIdx,nvar) to be shifted right
		for (j = varIdx; j < nvar; ++j) {
			curDeg = GET_NTH_EXP(degs, oldMasks[j], oldSizes[j]);
			if (curDeg > maxExps[j+1]) {
				free(tmpElems);
				unpackExponentVectors_AAZ_inp(aa);
				homogenizePolynomial_AAZ_inp_unpk(aa, varIdx);
				return;
			}
			curTDeg += curDeg;
			tmpElems[i].degs |= (curDeg << newSizes[j+1]);
		}

		//set the homogenizing variable's degree for this monomial
		curDeg = tdeg - curTDeg;
		if (curDeg > maxExps[varIdx]) {
			free(tmpElems);
			unpackExponentVectors_AAZ_inp(aa);
			homogenizePolynomial_AAZ_inp_unpk(aa, varIdx);
		}
		tmpElems[i].degs |= (curDeg << newSizes[varIdx]);
	}

	//copy back to original poly now
	memcpy(aa->elems, tmpElems, sizeof(AAZElem_t)*size);
	free(tmpElems);
	aa->nvar = nvar+1;

	mergeSortPolynomial_AAZ(aa);
}

AltArrZ_t* homogeneousGCD_AAZ(const AltArrZ_t* aa, const AltArrZ_t* bb) {

	if (isZero_AAZ(aa) || isZero_AAZ(bb)) {
		return NULL;
	}

	if (!isHomogeneous_AAZ(aa) || !isHomogeneous_AAZ(bb)) {
		return NULL;
	}

	if (aa->nvar != bb->nvar) {
		return NULL;
	}

	int nvar = aa->nvar;
    int varsA[nvar];
    int varsB[nvar];
    for (int k = 0; k < nvar; ++k) {
    	varsA[k] = 0;
    	varsB[k] = 0;
    }
    nonZeroVariables_AAZ(aa, varsA);
    nonZeroVariables_AAZ(bb, varsB);
    int doShrink;
    for (int k = 0; k < nvar; ++k) {
    	if (varsA[k] == 0) {
    		if (varsB[k] == 0) {
    			doShrink = 1;
    		} else {
    			return NULL;
    		}
    	} else if (varsB[k] == 0) {
    		return NULL;
    	}
    }

    AltArrZ_t* a = deepCopyPolynomial_AAZ(aa);
    AltArrZ_t* b = deepCopyPolynomial_AAZ(bb);
    mpz_t ca, cb;
    mpz_init(ca);
    mpz_init(cb);
    primitivePartAndContent_AAZ_inp(a, ca);
    primitivePartAndContent_AAZ_inp(b, cb);
    mpz_gcd(ca, ca, cb);

	degree_t ldegs_A[nvar];
	degree_t ldegs_B[nvar];
	removeLowDegrees_AAZ_inp(a, ldegs_A);
	removeLowDegrees_AAZ_inp(b, ldegs_B);
	degree_t ldegs_G[nvar];
	for (int i = 0; i < nvar; ++i) {
		ldegs_G[i] = MIN(ldegs_A[i], ldegs_B[i]);
	}

    //if the polys have the same variables with positive degree
    //and the same variables of degree 0, remove them.
	int varMap[nvar];
    int trueNvar;
    if (doShrink) {
		trueNvar = tryShrinkVariables_AAZ_inp(a, varMap);
		tryShrinkVariables_AAZ_inp(b, varMap);
    }

    //choose the variable to remove based on the max. of min partial degrees
    //i.e., the variable whose degree could be highest in the GCD.
    int varIdx;
    degree_t maxDeg, curMin;
    partialDegrees_AAZ(a, ldegs_A); //reuse ldegs arrays
    partialDegrees_AAZ(b, ldegs_B);
    curMin = ldegs_A[0] < ldegs_B[0] ? ldegs_A[0] : ldegs_B[0]; //get min of 0th index as cur max
    maxDeg = curMin;
    varIdx = 0;
    for (int i = 1; i < trueNvar; ++i) {
    	curMin = ldegs_A[0] < ldegs_B[0] ? ldegs_A[0] : ldegs_B[0];
    	if (curMin > maxDeg) {
    		maxDeg = curMin;
    		varIdx = i;
    	}
    }

	dehomogenizePolynomial_AAZ_inp(a, varIdx);
	dehomogenizePolynomial_AAZ_inp(b, varIdx);

	//recursively get GCD
	AltArrZ_t* g = gcd_AAZ_tmp(a, b);

	homogenizePolynomial_AAZ_inp(g, varIdx);

	reverseShrinkVariables_AAZ_inp(g, nvar, varMap);
	multiplyByMonomial_AAZ_inp(g, ldegs_G); //add back common lowdegrees
	multiplyByInteger_AAZ_inp(g, ca);

	mpz_clear(ca);
	mpz_clear(cb);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);

	return g;
}


//forward declare "private" method
void evalOneVariable_AAZ_inp(AltArrZ_t* aa, int k, degree_t deg, const mpz_t p);

//a and b should be primitive and the same number of variables
//returns 1 if a and b are co-prime.
//returns 0 if a and b may have a common divisor
//returns -1 if the primality test was too expensive
int coprimalityHeuristic_AAZ(const AltArrZ_t* aa, const AltArrZ_t* bb) {
	if (isZero_AAZ(aa) || isZero_AAZ(bb) || aa->nvar != bb->nvar) {
		return 0;
	}

	int nvar = aa->nvar;

	//only evaluate vars where both a and b have positive degree
	degree_t aDegs[nvar];
	degree_t bDegs[nvar];


	partialDegrees_AAZ(aa, aDegs);
	partialDegrees_AAZ(bb, bDegs);

	int allZero = 0;
	for (int i = 0; i < nvar; ++i) {
		allZero |= aDegs[i];
		allZero |= bDegs[i];
	}

	if (allZero == 0) {
		//no variables in common;
		return 1;
	}

	mpz_t ma, mb, tmp;
	mpz_inits(ma, mb, tmp, NULL);

	int ret = 1;

	//copy for in-place evaluations below
	AltArrZ_t* a = deepCopyPolynomial_AAZ(aa);
	AltArrZ_t* b = deepCopyPolynomial_AAZ(bb);
	for (int i = 0; i < nvar; ++i) {
		if (aDegs[i] == 0 || bDegs[i] == 0) {
			//i'th variable is not in common or deg is 0.
			continue;
		}

		//use Cauchy's bound on roots for eval point
		infinityNorm_AAZ(a, ma);
		mpz_fdiv_q(ma, ma, a->elems->coef);
		infinityNorm_AAZ(b, mb);
		mpz_fdiv_q(mb, mb, b->elems->coef);
		mpz_abs(ma, ma);
		mpz_abs(mb, mb);
		if (mpz_cmp(mb, ma) > 0) {
			mpz_swap(ma, mb);
		}
		mpz_mul_ui(ma, ma, 2ul);
		mpz_add_ui(ma, ma, ((unsigned long) i) + 50ul); //add a little fiddle factor

		//bail out if it's too much. 100000 here is a heuristic value...
		//AB, 2020-11-24:
		//Some structured polynomials with very high degree gcds could
		//make this function work way too much for nothing
		//(e.g. in computeing a square free factorization) so let's avoid that.
		if (mpz_size(ma)*aDegs[i] > 100 || mpz_size(ma)*bDegs[i] > 100) {
			ret = -1;
			break;
		}

		//evaluate a
		evalOneVariable_AAZ_inp(a, i, aDegs[i], ma);
		if (isZero_AAZ(a)) {
			ret = 0;
			break;
		}

		//evaluate b
		evalOneVariable_AAZ_inp(b, i, bDegs[i], ma);
		if (isZero_AAZ(b)) {
			ret = 0;
			break;
		}

		primitivePartAndContent_AAZ_inp(a, tmp);
		primitivePartAndContent_AAZ_inp(b, mb);
		mpz_gcd(tmp, tmp, mb);
		mpz_mul_ui(tmp, tmp, 2ul);
		//if gcd of contents is >= evaluation point,
		//theres likely a non-trivial gcd
		if (mpz_cmp(tmp, ma) >= 0) {
			// gmp_fprintf(stderr, "i: %d, gcd: %Zd, ma: %Zd\n", i, tmp, ma);
			ret = 0;
			break;
		}
	}

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	mpz_clears(ma, mb, tmp, NULL);
	return ret;
}


//if varIdx = nvar - 1, g must be at least an allocated polynomial with nvar set correctly.
//otherwise, g is gamma that we are updating.
void genPoly_gcdheu_AAZ_inp(AltArrZ_t* gamma, int varIdx, mpz_t eps, mpz_t halfeps) {
	int nvar = gamma->nvar;
	if (varIdx + 1 == nvar) {
		mpz_t gammaZ;
		mpz_init_set(gammaZ, gamma->elems->coef);
		mpz_clear(gamma->elems->coef); //avoid double init
		genPoly_univgcdheu_AAZ_inp(gamma, eps, halfeps, gammaZ);
		mpz_clear(gammaZ);
		return;
	}

	// gamma->size = 1;
	// mpz_set_ui(gamma->elems->coef, 1ul);
	// gamma->elems->degs = 0ull;
	// return;

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);

	mpz_t infn;
	mpz_init(infn);
	infinityNorm_AAZ(gamma, infn);
	// gmp_fprintf(stderr, "infn: %Zd, eps: %Zd\n", infn, eps);

	long int exp1, exp2;
	double d1 = mpz_get_d_2exp(&exp1, infn);
	double d2 = mpz_get_d_2exp(&exp2, eps);
	d1 = log(d1)  + log(2) * (double) exp1;
	d2 = log(d2)  + log(2) * (double) exp2;
	int maxDeg = (int) ceil(d1/d2);
	// old way: mpz_get_d can cause double to overflow
	// int maxDeg = (int) ceil(log(mpz_get_d(infn)) / log(mpz_get_d(eps)));

	maxDeg += 2; //+2 for buffer room in array
	mpz_clear(infn);

	AAZElem_t* chunks[maxDeg];
	AAZElem_t* chunk;
	int chunkSizes[maxDeg];
	int totalSize = 0;

	int exp, i, gammaSize;
	int nextGammaSize = gamma->size;
	int insertIdx = 0;
	for (exp = 0; nextGammaSize > 0; ++exp) {
		gammaSize = nextGammaSize;
		nextGammaSize = 0;
		int gamma1 = gammaSize+1;
		// fprintf(stderr, "gammaSize: %d, gammaSize+1: %d\n", gammaSize, gamma1);
		chunk = (AAZElem_t*) malloc(sizeof(AAZElem_t)*(gammaSize+1));
		chunks[exp] = chunk;

		//g_exp = chunk[exp] = symmetric mod gamma.
		//simultaneously do gamma = (gamma - g_i) / eps
		//simultaneously multiply g_exp by x_{varidx}^exp
		insertIdx = 0;
		mpz_init(chunk[insertIdx].coef);
		//interate over original number of gamma coefs
		for (i = 0; i < gamma->size; ++i) {
			if (mpz_sgn(gamma->elems[i].coef) == 0) {
				//since we are updating gamma in-place, lazily,
				//some coefs may become 0 before others.
				continue;
			}

			mpz_mod(chunk[insertIdx].coef, gamma->elems[i].coef, eps);
			if(mpz_cmp(chunk[insertIdx].coef, halfeps) > 0) {
				mpz_sub(chunk[insertIdx].coef, chunk[insertIdx].coef, eps);
			}

			mpz_sub(gamma->elems[i].coef, gamma->elems[i].coef, chunk[insertIdx].coef);
			mpz_divexact(gamma->elems[i].coef, gamma->elems[i].coef, eps);
			if (mpz_sgn(gamma->elems[i].coef) != 0) {
				++nextGammaSize;
			}

			//monomial *= x^exp
			chunk[insertIdx].degs = ((degrees_t) exp << sizes[varIdx]) | gamma->elems[i].degs;

			if (mpz_sgn(chunk[insertIdx].coef) != 0) {
				++insertIdx;
				mpz_init(chunk[insertIdx].coef); //HUH MAYBE this goes past the en dof the array..
			}
		}
		mpz_clear(chunk[insertIdx].coef);
		chunkSizes[exp] = insertIdx;
		totalSize += insertIdx;
	}


	//cleanup and prepare gamma to receive new poly
	for (i = 0; i < gamma->size; ++i) {
		mpz_clear(gamma->elems[i].coef);
	}
	free(gamma->elems);

	//we now have valid chunks[i] for i < exp;
	//the chunks need to be reversed in order for canonical form
	AAZElem_t* finalGamma = (AAZElem_t*) malloc(sizeof(AAZElem_t)*totalSize);
	gamma->elems = finalGamma;
	gamma->alloc = totalSize;
	for (i = exp-1; i >= 0; --i) {
		memcpy(finalGamma, chunks[i], sizeof(AAZElem_t)*chunkSizes[i]);
		finalGamma += chunkSizes[i];
		//since we memcpy from chunks, we don't clear the coefs
		free(chunks[i]);
	}
	gamma->size = totalSize;

}


//g must be non-NULL
int heuPrimGCD_inner_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, AltArrZ_t** g, int varIdx, mpz_t ma, mpz_t mb) {

	// fprintf(stderr, "\n@@@@@@@@@@@@@@@@@@@@@@@@@@\nEntering heuPrimGCD with varIdx %d\n\n", varIdx);
	int nvar = a->nvar;
	int failed = 0;
	int dega, degb;

	// degree_t degsA[nvar];
	// degree_t degsB[nvar];
	// partialDegrees_AAZ(a, degsA);
	// partialDegrees_AAZ(a, degsB);
	// int anyInCommon = 0;
	// for (int i = 0; i < nvar; ++i) {
	// 	if(degsA[i] > 0 &&  degsB[i] > 0) {
	// 		anyInCommon = 1;
	// 	}
	// }
	// if (!anyInCommon) {
	// 	return makeConstIntPolynomial_AAZ(1, nvar, 1l);
	// }

	//compute the classic bound from GCDHEU 1986 paper
	infinityNorm_AAZ(a, ma);
	infinityNorm_AAZ(b, mb);
	if (mpz_cmp(ma, mb) > 0) {
		mpz_swap(ma, mb); //make ma the min.
	}

	mpz_t eps;	//keep eps stack local, don't re-use ma or mb
	mpz_init(eps);
	mpz_mul_ui(eps, ma, 2ul);
	mpz_add_ui(eps, eps, 29ul);

	//compute cauchy bound on a and b
	mpz_fdiv_q(ma, ma, a->elems->coef);
	mpz_fdiv_q(mb, mb, b->elems->coef);
	mpz_abs(ma, ma);
	mpz_abs(mb, mb);
	if (mpz_cmp(ma, mb) > 0) {
		mpz_swap(ma, mb); // make ma the min
	}
	mpz_mul_ui(ma, ma, 2ul);
	mpz_add_ui(ma, ma, 2ul);

	//use the larger of the two bounds
	if (mpz_cmp(eps, ma) < 0) {
		mpz_swap(eps, ma);
	}

	dega = mainDegree_AAZ(a);
	degb = mainDegree_AAZ(b);
	int maxdeg = dega > degb ? dega : degb;

	for (int attempts = 0; attempts < HEURISTIC_GCD_ATTEMPTS; ++attempts) {
		if (mpz_sizeinbase(eps, 2) * maxdeg > HEURISTIC_GCD_CUTOFF_SIZE) {
			failed = 1;
			break;
		} else {
			failed = 0;
		}

		//eval a and b
		AltArrZ_t* aa = deepCopyPolynomial_AAZ(a);
		AltArrZ_t* bb = deepCopyPolynomial_AAZ(b);
		evalOneVariable_AAZ_inp(aa, varIdx, dega, eps);
		evalOneVariable_AAZ_inp(bb, varIdx, degb, eps);
		// fprintf(stderr, "\n@@@@@@@@@@@@@@@@@@@@@@@@@@\nEvaluated one var varIdx %d\n\n", varIdx);

		//if we are at the bottom of the recursion
		if (varIdx+1 == nvar) {
			//notice that if only one of a or b is unpacked, then the gcd
			//has to be packed because it can't have degrees higher than
			//the 'smaller' one.
			if (*g != NULL) {
				freePolynomial_AAZ(*g);
			}
			if (a->unpacked && b->unpacked) {
				*g = makePolynomial_AAZ_unpk(10, nvar);
			} else {
				*g = makePolynomial_AAZ(10, nvar);
			}
			mpz_init((*g)->elems->coef);
			mpz_gcd((*g)->elems->coef, aa->elems->coef, bb->elems->coef);
		} else {
			//recursive call
			//make sure "going down" that g is a "clean slate"
			if (*g != NULL) {
				freePolynomial_AAZ(*g);
				*g = NULL;
			}
			failed = heuPrimGCD_inner_AAZ(aa, bb, g, varIdx + 1, ma, mb);
		}


		freePolynomial_AAZ(aa);
		freePolynomial_AAZ(bb);
		if (failed) {
			// fprintf(stderr, "Failed breaking early\n");
			freePolynomial_AAZ(*g);
			*g = NULL;
			break;
		}

		// fprintf(stderr, "\n@@@@@@@@@@@@@@@@@@@@@@@@@@\nEntering genPoly with varIdx %d\n\n", varIdx);

		//gen poly, ma is half eps
		mpz_fdiv_q_ui(ma, eps, 2ul);
		if ((*g)->unpacked) {
			genPoly_gcdheu_AAZ_inp_unpk(*g, varIdx, eps, ma);
		} else {
			genPoly_gcdheu_AAZ_inp(*g, varIdx, eps, ma);
		}
		// AB, 10/15/2020: Actually don't do prim part here
		// in case some variable has degree 0 in one of the input polys,
		// that will cause g to just be an integer as the recursion unwinds
		// and then get transformed to be 1. Only do primpart at very end.
		// primitivePart_AAZ_inp(*g);

		//divide check
		//TODO one can genPoly with co factors
		AltArrZ_t* cofa = NULL;
		AltArrZ_t* cofb = NULL;
		if (divideTest_AAZ(a, *g, &cofa, nvar) && divideTest_AAZ(b, *g, &cofb, nvar)) {
			freePolynomial_AAZ(cofa);
			freePolynomial_AAZ(cofb);
			failed = 0;
			break;
		} else {
			//free here in case g divides a but not the other.
			freePolynomial_AAZ(cofa);
			freePolynomial_AAZ(cofb);
			failed = 1;
		}

		//else fail, update eps and try again
		mpz_mul_ui(eps, eps, 73794ul);
		mpz_fdiv_q_ui(eps, eps, 27011ul);
	}

	mpz_clear(eps);

	return failed;
}

/**
 * Try to compute the GCD of a and b by a sort of Kronecker substitution
 * and then GCDs over Z.
 * In practice, given a highly effective GCD based on hensel lifting, this
 * is rarely useful.
 */
AltArrZ_t* heuristicPrimitiveGCD_AAZ(const AltArrZ_t* a, const AltArrZ_t* b) {

	if (isZero_AAZ(a) || isZero_AAZ(b)) {
		return NULL;
	}


	if (a->nvar != b->nvar) {
		fprintf(stderr, "BPAS ERROR: input polynomials to heuristicGCD_AAZ must be defined in the same polynomial ring\n");
		exit(1);
	}

	if (a->nvar == 1) {
		return univariatePrimitiveGCD_AAZ(a,b);
	}

	mpz_t ma;
	mpz_t mb;
	mpz_init(ma);
	mpz_init(mb);
	infinityNorm_AAZ(a, ma);
	infinityNorm_AAZ(b, mb);
	if (mpz_cmp(ma, mb) > 0) {
		mpz_swap(ma, mb);
	}
	if (mpz_sizeinbase(ma, 2) > HEURISTIC_GCD_CUTOFF_INFN) {
		// don't even bother trying if the inf norm is too high
		return NULL;
	}

	AltArrZ_t* g = NULL;

	//returns g as NULL if not successful.
	int failed = heuPrimGCD_inner_AAZ(a, b, &g, 0, ma, mb);
	if (failed) {
		//make sure g is NULL;
		freePolynomial_AAZ(g);
		g = NULL;
	} else {
		primitivePart_AAZ_inp(g);
	}

	mpz_clear(ma);
	mpz_clear(mb);
	return g;
}

AltArrZ_t* bivariatePrimitiveGCD_AAZ(const AltArrZ_t* aa, const AltArrZ_t* bb) {
	if (aa->nvar != 2 || bb->nvar != 2) {
		fprintf(stderr, "BPAS ERROR in bivariatePrimitiveGCD_AAZ, input polys are not bivariate.\n");
		exit(1);
	}
	if (isExactlyEqual_AAZ(aa, bb)) {
		return deepCopyPolynomial_AAZ(aa);
	}

	// char* syms[] = {"x", "y"};
	// fprintf(stderr, "a:=\n");
	// printPoly_AAZ(stderr, aa, syms, 2);
	// fprintf(stderr, ";\nb:=\n");
	// printPoly_AAZ(stderr, bb, syms, 2);
	// fprintf(stderr, ";\n\n");

	//It is undefined to call this function when the low degrees of aa or bb are not removed.

	//unpacked is not needed in a bivariate world so force packing
	int packedA = 0, packedB = 0;
	if (aa->unpacked) {
		tryPackExponentVectors_AAZ_inp((AltArrZ_t*) aa);
		packedA = 1;
	}
	if (bb->unpacked) {
		tryPackExponentVectors_AAZ_inp((AltArrZ_t*) bb);
		packedB = 1;
	}

	//Here we consider aa, bb \in Z[x > y]. Thus x has index 0, y index 1 in exponent vector.
	//Therefore, we will be performing evaluations and interpolations in y.
	//The degree of x controls the number of interpolations since each coefficient w.r.t x
	//is a univariate polynomial in y. The degree of y controls the number of points used
	//in each interpolation.

	int nvar = 2;
	degree_t pdegs_A[nvar];
	degree_t pdegs_B[nvar];
	partialDegrees_AAZ(aa, pdegs_A);
	partialDegrees_AAZ(bb, pdegs_B);
	if ((pdegs_A[0] == 0 || pdegs_B[0] == 0) && (pdegs_A[1] == 0 || pdegs_B[1] == 0)) {
		return makeConstIntPolynomial_AAZ(1, nvar, 1l);
	}

	degree_t pdegs_G[nvar];
	pdegs_G[0] = MIN(pdegs_A[0], pdegs_B[0]);
	pdegs_G[1] = MIN(pdegs_A[1], pdegs_B[1]);
	int varMap[2];
	varMap[0] = 1;
	varMap[1] = 0;
	int reorderVars = 0;
	if (pdegs_G[0] < pdegs_G[1]) {
		//if degree in x is smaller, swap. We want to do many smaller interpolations later on.
		reorderVars_AAZ((AltArrZ_t*) aa, varMap, nvar);
		reorderVars_AAZ((AltArrZ_t*) bb, varMap, nvar);

		degree_t tmp;
		tmp = pdegs_A[0];
		pdegs_A[0] = pdegs_A[1];
		pdegs_A[1] = tmp;

		tmp = pdegs_B[0];
		pdegs_B[0] = pdegs_B[1];
		pdegs_B[1] = tmp;

		pdegs_G[0] = MIN(pdegs_A[0], pdegs_B[0]);
		pdegs_G[1] = MIN(pdegs_A[1], pdegs_B[1]);

		reorderVars = 1;
	}

	//get lc w.r.t x and y
	//s is coefAy, t is coefBy, the lcoeffs if y was the main var.
	mpz_t s;
	mpz_init_set(s, aa->elems[aa->size-1].coef);
	int dy = 0;
	for (int i = 0; i < aa->size; ++i) {
		if ((aa->elems[i].degs & EXP_2_V2) > dy) {
			dy = aa->elems[i].degs & EXP_2_V2;
			mpz_set(s, aa->elems[i].coef);
		}
	}
	mpz_t t;
	mpz_init_set(t, bb->elems[bb->size-1].coef);
	dy = 0;
	for (int i = 0; i < bb->size; ++i) {
		if ((bb->elems[i].degs & EXP_2_V2) > dy) {
			dy = bb->elems[i].degs & EXP_2_V2;
			mpz_set(t, bb->elems[i].coef);
		}
	}

	//get lcoeff of gcd
	mpz_t lcg, tcg;
	mpz_init(lcg);
	mpz_gcd(lcg, aa->elems->coef, bb->elems->coef);
	mpz_init(tcg);
	mpz_gcd(tcg, aa->elems[aa->size-1].coef, bb->elems[bb->size-1].coef);

	//get bad primes as product of coefs which should not vanish mod p
	mpz_t badPrimes;
	mpz_init_set(badPrimes, aa->elems->coef);
	mpz_mul(badPrimes, badPrimes, bb->elems->coef);
	mpz_divexact(badPrimes, badPrimes, lcg);
	mpz_mul(badPrimes, badPrimes, s);
	mpz_mul(badPrimes, badPrimes, t);
	mpz_gcd(s, s, t);
	mpz_divexact(badPrimes, badPrimes, s);

	//normalization factor; if 1 use lcoeff normalization, else tcoeff normalizaiton
	int use_lcoeff = 1;
	if (mpz_size(tcg)*2 < mpz_size(lcg)) {
		//*2 because it costs extra to find the trailing coef
		use_lcoeff = 0;
		mpz_mul(badPrimes, badPrimes, tcg);
	}

	//loop variables
	DBZP_t* G = allocPolynomial_DBZP(pdegs_G[0]+1, pdegs_G[1]+1);
	for (int i = 0; i <= pdegs_G[0]; ++i) {
		for (int j = 0; j <= pdegs_G[1]; ++j) {
			mpz_init(G->coefsY[i]->coefs[j]);
		}
		G->coefsY[i]->lt = pdegs_G[1];
	}
	G->ltX = pdegs_G[0];
	AltArrZ_t* retG = NULL;
	int resetG = 0;

	mpz_t m, halfm, mpz_prime, tmp;
	mpz_init_set_ui(m, 1ul);
	mpz_inits(halfm, mpz_prime, tmp, NULL);

	degree_t dCy = -1;
	degree_t dGpy;
	const Prime_ptr* Pptr;
	long N = pdegs_G[1];
	long i, j, pt, isEqual;
	elem_t ppt = 0;

	dbspoly_t* Ap = allocPolynomial_spXY(pdegs_A[0]+1, pdegs_A[1]+1);
	dbspoly_t* Bp = allocPolynomial_spXY(pdegs_B[0]+1, pdegs_B[1]+1);
	duspoly_t* cAp = NULL;
	duspoly_t* cBp = NULL;
	duspoly_t* Cp = NULL;
	duspoly_t* gamma = NULL;

	duspoly_t* Api = NULL;
	duspoly_t* Bpi = NULL;
	duspoly_t* Gpi = NULL;
	dbspoly_t* Gxy = NULL;

	elem_t* Xp = (elem_t*) malloc(sizeof(elem_t)*(N+1));
	duspoly_t** Yp = (duspoly_t**) malloc(sizeof(duspoly_t*)*(N+1));

	for (int primeIdx = 0; primeIdx < n_prime64_ptr; ++ primeIdx) {
		Pptr = prime64_ptr + primeIdx;
		if (mpz_fdiv_ui(badPrimes, Pptr->prime) == 0) {
			continue;
		}

		//Ap = Array of Z_p[y];
		//Bp = Array of Z_p[y];
		convertFromAltArrZ_spXY_inp(aa, &Ap, Pptr);
		convertFromAltArrZ_spXY_inp(bb, &Bp, Pptr);

		// fprintf(stderr, "\nAp, Bp:\n" );
		// printPolynomialOutForm_spXY (Ap, Pptr);
		// printPolynomialOutForm_spXY (Bp, Pptr);

		// Ap = primPart(Ap, 'cAp')
		// Bp = primPart(Bp, 'cBp')
		mainPrimitivePart_spXY_inp(Ap, &cAp, Pptr);
		mainPrimitivePart_spXY_inp(Bp, &cBp, Pptr);

		//Cp = Gcd(cAp, cBp);
		freePolynomial_spX(&Cp);
		Cp = NULL;
		GCDInForm_spX (cAp, cBp, &Cp, Pptr);

		//make arguments about degree of this image and of previous images
		if (dCy < 0) {
			//this is out first time through the loop
			dCy = Cp->lt;
		} else if (Cp->lt > dCy) {
			continue;
		} else if (Cp->lt < dCy) {
			dCy = Cp->lt;
			//dont need to free G in this case actually,
			//we will re-use its space with new round of images
			resetG = 1;
			mpz_set_ui(m, 1ul);
			continue;
		}

		//gamma_v = gcd(lc(Ap), lc(Bp)) //a poly in Z_p[y].
		freePolynomial_spX(&gamma);
		gamma = NULL;
		GCDInForm_spX(Ap->coefsY[Ap->ltX], Bp->coefsY[Bp->ltX], &gamma, Pptr);

		//get bound on number of images for interpolation
		N = MIN(pdegs_A[1] - cAp->lt, pdegs_B[1] - cBp->lt);
		N = MIN(N, pdegs_G[1] - dCy + gamma->lt);

		//simple check first
		if (isEqual_spXY(Ap, Bp)) {
			Gxy = Ap;
			Ap = NULL;
		} else {
			//need N+1 points;
			pt = 2; //random starting point
			for (j = 0; j <= N && pt > 0; ) {
				ppt = smallprimefield_convert_in(pt, Pptr);

				//Api = Eval(Ap, y=ip); //poly in Z_p[x]
				evalCoefficientVariableInForm_spXY(Ap, &Api, ppt, Pptr);
				if (Api->lt < pdegs_A[0]) {
					pt = (pt + 1) % Pptr->prime;
					continue;
				}
				//Bpi = Eval(Bp, y=ip); //poly in Z_p[x]
				evalCoefficientVariableInForm_spXY(Bp, &Bpi, ppt, Pptr);
				if (Bpi->lt < pdegs_B[0]) {
					pt = (pt + 1) % Pptr->prime;
					continue;
				}
				if (isEqual_spX(Api, Bpi)) {
					//a weird evaluation cause them to coincide.
					pt = (pt + 1) % Pptr->prime;
					continue;
				}

				// fprintf(stderr, "\nApi, Bpi:\n");
				// printPolynomialOutForm_spX (Api, Pptr);  // TEST
				// printPolynomialOutForm_spX (Bpi, Pptr);  // TEST

				//Gpi = GCD(Api, Bpi);
				freePolynomial_spX(&Gpi);
				Gpi = NULL;
				GCDInForm_spX(Api, Bpi, &Gpi, Pptr);

				// fprintf(stderr, "\nGpi:\n");
				// printPolynomialOutForm_spX (Gpi, Pptr);  // TEST


				//make arguments about degree of this image and of previous images
				if (Gpi->lt == 0 && dCy == 0) {
					retG = makeConstIntPolynomial_AAZ(1, nvar, 1l);
					break;
				} else if (Gpi->lt > pdegs_G[0]) {
					//this image is not right.
					if (j > 0) {
						//then we got at least one correct image with deg pdegs_G[0], so skip just this image
						pt = (pt + 1) % Pptr->prime;
						continue;
					} else {
						if (pt > 10) {
							//we don't like this prime, so get a new one.
							pt = 0; //this triggers a new prime
							break;
						} else {
							//give this prime another chance, maybe the point was bad.
							pt = (pt + 1) % Pptr->prime;
							continue;
						}
					}
				} else if (Gpi->lt < pdegs_G[0]) {
					pdegs_G[0] = Gpi->lt;
					resetG = 1; //degree in G still dropped, so must reset it;
					if (primeIdx > 0 || j > 0) {
						//all previous images were wrong, restart
						mpz_set_ui(m, 1ul);
						pt = 0; //set pt to 0 to restart images
						break;

					}
				}

				//Store this point and image
				Xp[j] = ppt;

				//Gpi comes back monic so normalize lcoeff.
				ppt = evaluateInForm_spX(gamma, ppt, Pptr);
				scalarMulPolynomialInForm_spX_inp(&Gpi, ppt, Pptr);

				Yp[j] = Gpi;
				Gpi = NULL;

				pt = (pt + 1) % Pptr->prime;
				++j; //we got a successful image
			}
			if (retG != NULL && isOne_AAZ(retG)) {
				break;
			}
			if (pt == 0) {
				//we went through all points in the prime field and no luck,
				//or we found that our previous images were bad, so skip this prime
				for (int k = 0; k < j; ++k) {
					freePolynomial_spX(&Yp[k]);
				}
				continue;
			}

			interpPolyUnivarImagesInForm_spXY(Xp, CONSTCONSTCAST(duspoly_t, Yp), N+1, &Gxy, Pptr);
			for (j = 0; j <= N; ++j) {
				freePolynomial_spX(&Yp[j]);
			}

			mainPrimitivePart_spXY_inp(Gxy, NULL, Pptr);
		}

		dGpy = 0;
		if (dCy > 0) {
			for (j = 0; j <= Gxy->ltX; ++j) {
				if (Gxy->coefsY[j] != NULL) {
					mulPolynomialsInForm_spX_inp(Gxy->coefsY + j, Cp, Pptr);
					dGpy = dGpy < Gxy->coefsY[j]->lt ? Gxy->coefsY[j]->lt : dGpy;
				}
			}
		} else {
			for (j = 0; j <= Gxy->ltX; ++j) {
				if (Gxy->coefsY[j] != NULL) {
					dGpy = dGpy < Gxy->coefsY[j]->lt ? Gxy->coefsY[j]->lt : dGpy;
				}
			}
		}

		//make arguments about degree of this image and of previous images
		if (dGpy > pdegs_G[1]) {
			//a bad prime
			continue;
		} else if (dGpy < pdegs_G[1]) {
			//clear previous results, but continue with this image
			pdegs_G[1] = dGpy;
			//dont need to free G in this case actually,
			//we will re-use its space with new round of images
			if (mpz_cmp_ui(m, 1ul) != 0) {
				//we only need to cotninue if we actually reset.
				//it could be that this is just the first iteration
				resetG = 1;
				mpz_set_ui(m, 1ul);
				continue;
			}
		}

		//normalization of Gxy over Z_p
		if (use_lcoeff) {
			ppt = smallprimefield_convert_in(mpz_fdiv_ui(lcg, Pptr->prime), Pptr);
			j = Gxy->coefsY[Gxy->ltX]->lt;
			ppt = smallprimefield_mul(ppt, smallprimefield_inv(Gxy->coefsY[Gxy->ltX]->elems[j], Pptr), Pptr);
		} else {
			ppt = 0;
			for (i = 0; i <= Gxy->ltX && ppt == 0; ++i) {
				for (j = 0; Gxy->coefsY[i] != NULL && j <= Gxy->coefsY[i]->lt; ++j) {
					if (Gxy->coefsY[i]->elems[j] != 0) {
						ppt = smallprimefield_convert_in(mpz_fdiv_ui(tcg, Pptr->prime), Pptr);
						ppt = smallprimefield_mul(ppt, smallprimefield_inv(Gxy->coefsY[i]->elems[j], Pptr), Pptr);
						break;
					}
				}
			}
		}
		for (j = 0; j <= Gxy->ltX; ++j) {
			scalarMulPolynomialInForm_spX_inp(Gxy->coefsY + j, ppt, Pptr);
		}

		mpz_set_si(mpz_prime, Pptr->prime);
		mpz_gcdext(tmp, s, t, mpz_prime, m);

		mpz_mul(m, m, mpz_prime);
		mpz_sub_ui(halfm, m, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		//reset coefs to 0 but re-use space
		if (resetG) {
			for (i = 0; i <= pdegs_G[0]; ++i) {
				for (j = 0; j <= pdegs_G[1]; ++j) {
					mpz_set_ui(G->coefsY[i]->coefs[j], 0ul);
				}
				for (j = pdegs_G[1] + 1; j <= G->coefsY[i]->lt; ++j) {
					mpz_clear(G->coefsY[i]->coefs[j]);
				}
				G->coefsY[i]->lt = pdegs_G[1];
			}
			for (i = pdegs_G[0] + 1; i <= G->ltX; ++i) {
				freePolynomial_DUZP(G->coefsY[i]);
				G->coefsY[i] = NULL;
			}
			G->ltX = pdegs_G[0];
			resetG = 0;
		}

		//CRT
		isEqual = 1;
		unsigned long long modg_out;
		mpz_t* gCoefs;
		for (i = 0; i <= Gxy->ltX; ++i) {
			if (Gxy->coefsY[i] == NULL) {
				continue;
			}
			pt = Gxy->coefsY[i]->lt;
			gCoefs = G->coefsY[i]->coefs;
			for (j = 0; j <= pt; ++j) {
				//reuse t as temp var
				mpz_set(t, gCoefs[j]);

				modg_out = (unsigned long long) smallprimefield_convert_out(Gxy->coefsY[i]->elems[j], Pptr);
				mpz_sub_ui(gCoefs[j], gCoefs[j], modg_out);
				mpz_mul(gCoefs[j], gCoefs[j], s);
				mpz_mul(gCoefs[j], gCoefs[j], mpz_prime);
				mpz_add_ui(gCoefs[j], gCoefs[j], modg_out);
				mpz_mod(gCoefs[j], gCoefs[j], m);
				//normalize in symmetric range
				if(mpz_cmp(gCoefs[j], halfm) > 0) {
					mpz_sub(gCoefs[j], gCoefs[j], m);
				}

				//short circuit evaluation will skip the cmp?
				isEqual = isEqual && (mpz_cmp(t, gCoefs[j]) == 0);
			}
		}

		// fprintf(stderr, "primeIdx: %d, equal: %d\n", primeIdx, isEqual);
		if (isEqual) {
			//trial divisions
			convertFromDBZP_AAZ_inp(G, &retG);
			primitivePart_AAZ_inp(retG);
			if (divideTest_AAZ(aa, retG, NULL, nvar) && divideTest_AAZ(bb, retG, NULL, nvar)) {
				break;
			} else {
				resetG;
			}
		}
	}

	freePolynomial_spX(&cAp);
	freePolynomial_spX(&cBp);
	freePolynomial_spX(&Cp);
	freePolynomial_spX(&gamma);
	freePolynomial_spX(&Api);
	freePolynomial_spX(&Bpi);

	freePolynomial_spXY(Ap);
	freePolynomial_spXY(Bp);
	freePolynomial_spXY(Gxy);

	freePolynomial_DBZP(G);
	free(Xp);
	free(Yp);

	mpz_clears(badPrimes, lcg, tcg, s, t, m, halfm, mpz_prime, tmp, NULL);

	if (packedA) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) aa);
	}
	if (packedB) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) bb);
	}
	if (reorderVars) {
		reorderVars_AAZ((AltArrZ_t*) aa, varMap, nvar);
		reorderVars_AAZ((AltArrZ_t*) bb, varMap, nvar);
		reorderVars_AAZ((AltArrZ_t*) retG, varMap, nvar);

	}

	// fprintf(stderr, "g:=\n" );
	// printPoly_AAZ(stderr, retG, syms, 2);
	// fprintf(stderr, ";\n\n");
	return retG;
}

#if defined(WITH_MAPLE) && (WITH_MAPLE)
AltArrZ_t* gcdByMaple_AAZ(const AltArrZ_t* P, const AltArrZ_t* Q) {
	char* c_names[P->nvar];
	for (int i = 0; i < P->nvar; ++i) {
		c_names[i] = (char*) malloc(sizeof(char)*16);
		sprintf(c_names[i], "x_%d", i);
	}

	char* pBuff, *qBuff;
	size_t size;
	FILE* stream;

	stream = open_memstream (&pBuff, &size);
	printPoly_AAZ(stream, P, (const char**) c_names, P->nvar);
	// fprintf(stream, ":");
	fflush (stream);
	fclose (stream);

	stream = open_memstream (&qBuff, &size);
	printPoly_AAZ(stream, Q, (const char**) c_names, Q->nvar);
	// fprintf(stream, ":");
	fflush(stream);
	fclose(stream);

	char* gstr = (*gcd_maple_string)(pBuff, qBuff);

	AltArr_t* gQ = generate_altarr_var_defined(gstr, (const char**)c_names, P->nvar);
	AltArrZ_t* gz = deepCopyPolynomial_AAZFromAA(gQ);
	freePolynomial_AA(gQ);
	for (int i = 0; i < P->nvar; ++i) {
		free(c_names[i]);
	}
	return gz;
}
#endif


//A top level GCD algorithm.
//TODO: review order of steps (5,6,7,8,9)
AltArrZ_t* gcd_AAZ_tmp(const AltArrZ_t* aa, const AltArrZ_t* bb) {

	if (isZero_AAZ(aa)) {
		if (isZero_AAZ(bb)) {
			return NULL;
		}
		return deepCopyPolynomial_AAZ(bb);
	}

	if (isZero_AAZ(bb)) {
		return deepCopyPolynomial_AAZ(aa);
	}

	if (aa->nvar != bb->nvar) {
		fprintf(stderr, "BPAS ERROR: input polynomials to gcd_AAZ must be defined in the same polynomial ring\n");
		exit(1);
	}

	int nvar = aa->nvar;

	//Step 1: Remove integer content
	mpz_t ac, bc;
	mpz_init(ac);
	mpz_init(bc);
	AltArrZ_t* a = primitivePartAndContent_AAZ(aa, ac);
	AltArrZ_t* b = primitivePartAndContent_AAZ(bb, bc);
	mpz_gcd(ac, ac, bc);
	if (nvar == 0) {
		AltArrZ_t* g = makeConstPolynomial_AAZ(1, 0, ac);
		mpz_clear(ac);
		mpz_clear(bc);
		return g;
	}

	//Step 2: remove low degrees
	degree_t ldegs_A[nvar];
	degree_t ldegs_B[nvar];
	removeLowDegrees_AAZ_inp(a, ldegs_A);
	removeLowDegrees_AAZ_inp(b, ldegs_B);
	degree_t ldegs_G[nvar];
	for (int i = 0; i < nvar; ++i) {
		ldegs_G[i] = MIN(ldegs_A[i], ldegs_B[i]);
	}

	AltArrZ_t* g = NULL;
	degree_t pdegs_A[nvar];
	degree_t pdegs_B[nvar];
	degree_t varsInCommon[nvar];
	int anyInCommon = 0;
	int allInCommon = 1;

	//Step 3: simple checks
	if (divideTest_AAZ(b, a, NULL, nvar)) {
		g = deepCopyPolynomial_AAZ(a);
	} else if (totalDegree_AAZ(a) < 2) {
		g = makeConstIntPolynomial_AAZ(1, nvar, 1l);
	} else if (divideTest_AAZ(a, b, NULL, nvar)) {
		g = deepCopyPolynomial_AAZ(b);
	} else if (totalDegree_AAZ(b) < 2) {
		g = makeConstIntPolynomial_AAZ(1, nvar, 1l);
	} else {
		//Setp 3.2: check intersection of variables with positive degree
		partialDegrees_AAZ(a, pdegs_A);
		partialDegrees_AAZ(b, pdegs_B);

		for (int i = 0; i < nvar; ++i) {
			if(pdegs_A[i] > 0 &&  pdegs_B[i] > 0) {
				varsInCommon[i] = 1;
				anyInCommon = 1;
			} else if (pdegs_A[i] == 0 && pdegs_B[i] == 0) {
				varsInCommon[i] = 1;
			} else {
				varsInCommon[i] = 0;
				allInCommon = 0;
			}
		}
		if (!anyInCommon) {
			g = makeConstIntPolynomial_AAZ(1, nvar, 1l);
		}
	}


	//Step 4: try to find a substitution x -> x^(1/k) foreach x
	//TODO

	//Step 4.2: are they homogeneous? returns NULL if not homogeneous.
	if (g == NULL && allInCommon) {
		g = homogeneousGCD_AAZ(a, b);
	}

	// char* syms[] = {"x", "y", "z", "s", "t", "u", "v", "w", "a", "b", "c"};
	// fprintf(stderr, "g before loop: ");
	// printPoly_AAZ(stderr, g, syms, nvar);
	// fprintf(stderr, "\n" );

	while ( g == NULL ) {
		//Step 5: call univar gcd if possible
		if (nvar == 1) {
			g = univariatePrimitiveGCD_AAZ(a, b);
			break;
		}
		//Step 6: call bivar gcd if possible
		else if (nvar == 2) {
			g = bivariatePrimitiveGCD_AAZ(a, b);
			break;
		}

		// Step 7: Heuristic co-primality check
		int coprime = coprimalityHeuristic_AAZ(a, b);
		// fprintf(stderr, "coprime heruistic said: %d\n", coprime);
		if (coprime == 1) {
			g = makeConstIntPolynomial_AAZ(1, nvar, 1ul);
			break;
		}

		//Step 8: Heuristic gcd
		g = heuristicPrimitiveGCD_AAZ(a, b);
		if (g != NULL) {
			// fprintf(stderr, "@@@@@@@@@@@@@@@@@@@\nRETURNING a Heuristic GCD\n\n");
			break;
		}
		// fprintf(stderr, "@@@@@@@@@@@@@@@@@@@\nFailed to get a Heuristic GCD\n\n");

		//Step 9: Get content of each polynomial w.r.t variables not in common.
		//        Use the content in place of the input polys, thus reducing
		//        number of varibales in each to their intersecton
		int uniqueVars[nvar];
		AltArrZ_t* tmp;
		while (!allInCommon) {
			//if vars(a) != vars(a) \cup vars(b)
			int allMatch = 1;
			for (int i = 0; i < nvar; ++i) {
				uniqueVars[i] = 0;
				if (pdegs_A[i] > 0 && varsInCommon[i] == 0) {
					uniqueVars[i] = 1;
					allMatch = 0;
				}
			}
			if (!allMatch) {
				tmp = contentInVars_AAZ(a, uniqueVars);
				freePolynomial_AAZ(a);
				a = tmp;
				partialDegrees_AAZ(a, pdegs_A);
			}

			//if vars(b) != vars(a) \cup vars(b)
			allMatch = 1;
			for (int i = 0; i < nvar; ++i) {
				uniqueVars[i] = 0;
				if (pdegs_B[i] > 0 && varsInCommon[i] == 0) {
					uniqueVars[i] = 1;
					allMatch = 0;
				}
			}
			if (!allMatch) {
				tmp = contentInVars_AAZ(b, uniqueVars);
				freePolynomial_AAZ(b);
				b = tmp;
				partialDegrees_AAZ(b, pdegs_B);
			}

			//update vars(a) \cup vars(b)
			anyInCommon = 0;
			allInCommon = 1;
			for (int i = 0; i < nvar; ++i) {
				varsInCommon[i] = 0;
				if (pdegs_A[i] > 0 && pdegs_B[i] > 0) {
					varsInCommon[i] = 1;
					anyInCommon = 1;
				} else if (pdegs_A[i] == 0 && pdegs_B[i] == 0) {
					varsInCommon[i] = 1;
				} else {
					allInCommon = 0;
				}
			}
		}
		// At this point, all vars in a and b are either degree 0 or in common.

		//Step 10: Re-do simple checks step 3
		if (!anyInCommon) {
			g = makeConstIntPolynomial_AAZ(1, nvar, 1l);
			break;
		} else if (divideTest_AAZ(b, a, NULL, nvar)) {
			g = deepCopyPolynomial_AAZ(a);
			break;
		} else if (totalDegree_AAZ(a) < 2) {
			g = makeConstIntPolynomial_AAZ(1, nvar, 1l);
			break;
		} else if (divideTest_AAZ(a, b, NULL, nvar)) {
			g = deepCopyPolynomial_AAZ(b);
			break;
		} else if (totalDegree_AAZ(b) < 2) {
			g = makeConstIntPolynomial_AAZ(1, nvar, 1l);
			break;
		}

		int varMap[nvar];
		int trueNvarA = tryShrinkVariables_AAZ_inp(a, varMap);
		int trueNvarB = tryShrinkVariables_AAZ_inp(b, varMap);
		if (trueNvarA != trueNvarB) {
			//this should never ever happen based on the previous while loop.
			fprintf(stderr, "BPAS AAZ GCD ERROR: Shrink maps not equivalent");
			exit(1);
		}

		//Step 11: call univar gcd if possible
		if (trueNvarA == 1) {
			g = univariatePrimitiveGCD_AAZ(a, b);
			reverseShrinkVariables_AAZ_inp(g, nvar, varMap);
			break;
		}
		//Step 12: call bivar gcd if possible
		else if (trueNvarA == 2) {
			// fprintf(stderr, "CALLING BIVAR\n" );
			g = bivariatePrimitiveGCD_AAZ(a, b);
			reverseShrinkVariables_AAZ_inp(g, nvar, varMap);
			break;
		}
		//Step 13: Hensel (or Maple fall back)
		else {
			// fprintf(stderr, "calling primitive gcd hensel\n" );
#if defined(WITH_NTL) && WITH_NTL
			g = primitiveGCDHensel_AAZ(a, b);
#elif defined(WITH_MAPLE) && WITH_MAPLE
			g = gcdByMaple_AAZ(a, b);
#else
			g = makeConstIntPolynomial_AAZ(1,a->nvar,1l);
#endif
			reverseShrinkVariables_AAZ_inp(g, nvar, varMap);

			break;

		}
	}


	multiplyByMonomial_AAZ_inp(g, ldegs_G); //add back common lowdegrees
	multiplyByInteger_AAZ_inp(g, ac); //add back in integer content
	mpz_clear(ac);
	mpz_clear(bc);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);

	// fprintf(stderr, "g before return: ");
	// printPoly_AAZ(stderr, g, syms, nvar);
	// fprintf(stderr, "\n" );
	return g;

}



/*****************
 * Evaluation
 *****************/

void univarEvaluate_AAZ(const AltArrZ_t* aa, const mpz_t point, mpz_t res) {
	if (aa->size == 0) {
		mpz_set_si(res, 0l);
		return;
	}

	if (aa->nvar != 1) {
		fprintf(stderr, "Poly is not univariate in univarEvaluate_AAZ\n");
		exit(1);
	}

	if (aa->unpacked) {
		return univarEvaluate_AAZ_unpk(aa, point, res);
	}

	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	mpz_set(res, elems->coef);
	degrees_t prevDeg = elems->degs;
	degrees_t nextDeg;
	for (int i = 1; i < size; ++i) {
		nextDeg = elems[i].degs;
		for (degrees_t j = prevDeg; j > nextDeg; --j) {
			mpz_mul(res, res, point);
		}
		mpz_add(res, res, elems[i].coef);
		prevDeg = nextDeg;
	}
	for (degrees_t j = prevDeg; j > 0; --j) {
		mpz_mul(res, res, point);
	}
}

//For the polynomial aa, evaluate the k'th variable,
//whose partial degree is deg, at the point p.
//The evaluation is done in-place.
void evalOneVariable_AAZ_inp(AltArrZ_t* aa, int k, degree_t deg, const mpz_t p) {
	if (isZero_AAZ(aa)) {
		return;
	}

	mpz_t powers[deg];
	mpz_init_set(powers[1], p);
	for (int i = 2; i <= deg; ++i) {
		mpz_init(powers[i]);
		mpz_mul(powers[i], powers[i-1], p);
	}

	AAZElem_t* elems = aa->elems;
	int nvar = aa->nvar;
	int size = aa->size;

	int pdeg;
	for (int i = 0; i < size; ++i) {
		pdeg = partialDegreeTerm_AAZ(aa, i, k);
		if (pdeg > 0) {
			mpz_mul(elems[i].coef, elems[i].coef, powers[pdeg]);
			setPartialDegreeTerm_AAZ(aa, i, k, 0);
		}
	}

	for (int i = 1; i <= deg; ++i) {
		mpz_clear(powers[i]);
	}

	mergeSortPolynomial_AAZ(aa);
	condensePolynomial_AAZ(aa);
}



/*****************
 * Conversion for Arithemtic
 *****************/

/**
 * Convert a multivariate polynomial f to a univariate polynomial
 * under some multivariate Kronecker homomorphism.
 * This homomorphism is Z[x_1,...,x_v] -> Z[x] where each monomial maps as follows:
 * x_1^e1x_2^e2...x_v^ev -> x^{e1*alpha1 + e2*alpha2 + ... + ev},
 * where alpha_i = (maxDegs[i+1] + 1)*alpha_{i+1}, and alpha_v = 1.
 *
 * @param f: the SMZP to convert
 * @param maxDegs: the degrees which define the Kronecker substitution.
 *
 * @return a DUZP image of f under the homomorphism.
 */
DUZP_t* convertToDUZP_KS_AAZ(const AltArrZ_t* f, degree_t* maxDegs) {
	if (isZero_AAZ(f)) {
		return NULL;
	}

	degree_t curPow = 1;

	degree_t deg, pdeg, alpha;
	int retDeg = 0;
	int i, j;
	int size = f->size;
	int nvar = f->nvar;
	AAZElem_t* elems = f->elems;

	const degrees_t* masks = NULL;
	const int* sizes;
	if (!f->unpacked) {
		masks = getExpMaskArray(nvar);
		sizes = getExpOffsetArray(nvar);
	}

	int maxDeg = 1;
	for (j = 0; j < nvar; ++j) {
		maxDeg *= (maxDegs[j] + 1);
	}
	//maxDeg -= 1; //needed to be accruate, but let's just hold maxDeg+1 in maxDeg
	DUZP_t* ret = makePolynomial_DUZP(maxDeg);

	//yes, i=1, because makePoly returns constant term = 0.
	for (i = 1; i < maxDeg; ++i) {
		mpz_init(ret->coefs[i]);
	}

	if (f->unpacked) {
		degree_t* degs = (degree_t*) f->elems->degs;
		for (i = 0; i < size; ++i) {
			deg = 0;
			alpha = 1;

			//yes do in reverse order to maintain lex ordering on inverse mapping.
			for (j = nvar-1; j >= 0; --j) {
				pdeg = degs[i*nvar + j];
				deg += pdeg*alpha;
				alpha *= (maxDegs[j] + 1);
			}

			mpz_set(ret->coefs[deg], elems[i].coef);
			retDeg = retDeg > deg ? retDeg : deg;
		}
	} else {
		for (i = 0; i < size; ++i) {
			deg = 0;
			alpha = 1;

			//yes do in reverse order to maintain lex ordering on inverse mapping.
			for (j = nvar-1; j >= 0; --j) {
				pdeg = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
				deg += pdeg*alpha;
				alpha *= (maxDegs[j] + 1);
			}

			mpz_set(ret->coefs[deg], elems[i].coef);
			retDeg = retDeg > deg ? retDeg : deg;
		}
	}

	for (i = retDeg + 1; i < maxDeg; ++i) {
		mpz_clear(ret->coefs[i]);
	}
	ret->lt = retDeg;

	return ret;
}



//it's assumed that the resulting poly can fit packed here

/**
 * Convert a univariate polynomial fd to a multivariate polynomial
 * by the inversion of some multivariate Kronecker homomorphism.
 * This homomorphism is Z[x_1,...,x_v] -> Z[x] where each monomial maps as follows:
 * x_1^e1x_2^e2...x_v^ev -> x^{e1*alpha1 + e2*alpha2 + ... + ev},
 * where alpha_i = (maxDegs[i+1] + 1)*alpha_{i+1}, and alpha_v = 1.
 *
 * @param fd: the DUZP to convert from
 * @param maxDegs: the degrees which define the Kronecker substitution.
 * @param nvar: the number of variables
 *
 * @return a SMZP pre-image of fd under the homomorphism.
 */
AltArrZ_t* convertFromDUZP_KS_AAZ(const DUZP_t* fd, degree_t* maxDegs, int nvar) {
	if (isZero_DUZP(fd)) {
		return NULL;
	}

	const int* sizes = getExpOffsetArray(nvar);

	int fAlloc = 50;
	AltArrZ_t* f = makePolynomial_AAZ(fAlloc, nvar);

	int curSize = 0;
	int i, j;

	degree_t maxAlpha = 1;
	for (j = 0; j < nvar; ++j) {
		maxAlpha *= (maxDegs[j] + 1);
	}

	degree_t alpha;
	degree_t r;
	degrees_t pdeg, mdeg;
	for (i = fd->lt; i >= 0; --i) {
		if (mpz_sgn(fd->coefs[i]) == 0) {
			continue;
		}

		alpha = maxAlpha;
		mdeg = 0ll;
		r = i;
		for (j = 0; j < nvar && r > 0; ++j) {
			//update alpha for current variable
			alpha /= (maxDegs[j] + 1);

			//get partial degree by reversing Kronecker sub
			pdeg = r / alpha;
			r %= alpha;

			//pack the partial degree
			mdeg |= pdeg << sizes[j];
		}

		if (curSize >= fAlloc) {
			fAlloc <<= 1;
			if (fAlloc > fd->lt) {
				fAlloc = fd->lt + 1;
			}
			f->size = curSize;
			resizePolynomial_AAZ(f, fAlloc);
		}

		f->elems[curSize].degs = mdeg;
		mpz_init_set(f->elems[curSize].coef, fd->coefs[i]);
		++curSize;
	}

	f->size = curSize;

	return f;

}


/**
 * Multiply the polynomials f and g using a Kronecker substitution
 * to univariate integer polynomials. These images are then multiplied
 * as univariate polynomials and the multivariate product
 * reconstructed from image product.
 *
 * @see (Moreno Maza and Xie, 2011)
 *
 * @note opeands must exist in the same polynomial ring (i.e. same number of variables).
 *
 */
AltArrZ_t* multiplyPolynomials_KS_AAZ(const AltArrZ_t* f, const AltArrZ_t* g) {
	if (isZero_AAZ(f) || isZero_AAZ(g)) {
		return NULL;
	}

	short nvar = f->nvar;
	if (nvar != g->nvar) {
		fprintf(stderr, "Error in multiplyPolynomials_KS_AAZ: operands must have same number of variables.\n");
		exit(1);
	}

	if (nvar == 1) {
		DUZP_t* fd = convertFromAltArrZ_DUZP(f);
		DUZP_t* gd = convertFromAltArrZ_DUZP(g);
		DUZP_t* pd = multiplyPolynomials_DUZP(fd, gd);
		AltArrZ_t* prod = convertToAltArrZ_DUZP(pd);
		freePolynomial_DUZP(fd);
		freePolynomial_DUZP(gd);
		freePolynomial_DUZP(pd);
		return prod;
	}

	degree_t fDegs[nvar];
	degree_t gDegs[nvar];
	degree_t prodDegs[nvar];
	partialDegrees_AAZ(f, fDegs);
	partialDegrees_AAZ(g, gDegs);

	int i;
	short k;
	for (k = 0; k < nvar; ++k) {
		prodDegs[k] = fDegs[k] + gDegs[k];
	}

	DUZP_t* fd = convertToDUZP_KS_AAZ(f, prodDegs);
	DUZP_t* gd = convertToDUZP_KS_AAZ(g, prodDegs);

	DUZP_t* pd = multiplyPolynomials_DUZP(fd, gd);

	AltArrZ_t* prod = NULL;

	if (isPackableExponentVector(prodDegs, nvar)) {
		prod = convertFromDUZP_KS_AAZ(pd, prodDegs, nvar);
	} else {
		prod = convertFromDUZP_KS_AAZ_unpk(pd, prodDegs, nvar);
	}

	return prod;
}



/****************
* Normal Form (Multi-Divisor Division (MDD))
*****************/

AltArrZ_t* onlyNormalForm_AAZ (AltArrZ_t* f, AltArrZ_t** G, int s, int nvar)
{
	if (s < 0){
		return NULL;
	}

	if (s == 0){
		return deepCopyPolynomial_AAZ (f);
	}

	if (s == 1){
		if (G[0] != NULL && G[0]->size != 0){
			AltArrZ_t* r1 = NULL;
			AltArrZ_t* q1 = NULL;
			dividePolynomials_AAZ (f, *G, &q1, &r1, nvar);
			freePolynomial_AAZ (q1);

			return r1;
		}

		return deepCopyPolynomial_AAZ (f);
	}

	int i; // i= 0 ... (s-1)
	int isDivided;
	AltArrZ_t* h = deepCopyPolynomial_AAZ (f);  // this algorithm runs until h != NULL
	AltArrZ_t* rem = NULL;
	AltArrZ_t* lt = NULL;
	AltArrZ_t* tmp_q;
	AltArrZ_t* tmp_r;

	while (h != NULL && h->size != 0) {

		i = 0;
		isDivided = 0;

		while (i < s && !isDivided)
			{
				if ( h->size != 0 && h!= NULL &&
					G[i] != NULL && G[i]->size != 0 && monomialDivideTest (h->elems[0].degs, G[i]->elems[0].degs, nvar))  {
					tmp_q = NULL;
					tmp_r = NULL;
					divideByLeadingTerms_AAZ (h, G[i], &tmp_q, &tmp_r, nvar);
					freePolynomial_AAZ (tmp_r);

					if (!isOne_AAZ (tmp_q)) {
						tmp_q = multiplyPolynomials_AAZ_inp (tmp_q, G[i], nvar);
						h = subPolynomials_AAZ_inp (h, tmp_q, nvar);
						freePolynomial_AAZ (tmp_q);
					} else {
						freePolynomial_AAZ (tmp_q);
						h = subPolynomials_AAZ_inp (h, G[i], nvar);
					}
					isDivided = 1;
				} else
					i = i + 1;
			}

		if (!isDivided) {
			if (h != NULL && h->size != 0){
				lt = leadingTerm_AAZ (h, nvar);
			} else {
				lt = NULL;
			}

			if (rem == NULL) {
				rem = lt;
			} else {
				if (lt != NULL) {
					rem = addPolynomials_AAZ_inp (rem, lt, nvar);
					freePolynomial_AAZ (lt);
				}
			}
			if (h != NULL && h->size != 0) {
				/* h = subPolynomials_AA_inp (h, lt, nvar); */
				subByLeadingTerm_AAZ (&h, nvar);
			}
		}
	}

	return rem;
}

void heapMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar)
{
	if (s < 0){
		fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0){
		*r = deepCopyPolynomial_AAZ (f);
		return;
	}

	if (s == 1){  // simple multivariate polynomial division
		if (G[0] != NULL && G[0]->size != 0){
			dividePolynomials_AAZ (f, *G, Q, r, nvar);
		} else {
			*r = deepCopyPolynomial_AAZ (f);
			Q[0] = NULL;
		}
		return;
	}

	int i; // i= 0 ... (s-1)
	int isDivided;
	AltArrZ_t* h = deepCopyPolynomial_AAZ (f);  // this algorithm runs until h != NULL
	AltArrZ_t* rem = NULL;
	AltArrZ_t* lt = NULL;
	AltArrZ_t* tmp_q;
	AltArrZ_t* tmp_r;

	while (h != NULL && h->size != 0) {
		i = 0;
		isDivided = 0;

		while (i < s && !isDivided)
			{
				if ( h!= NULL && h->size != 0 && G[i] != NULL && G[i]->size != 0  && monomialDivideTest (h->elems[0].degs, G[i]->elems[0].degs, nvar))  {
					tmp_q = NULL;
					tmp_r = NULL;

					divideByLeadingTerms_AAZ (h, G[i], &tmp_q, &tmp_r, nvar);
					freePolynomial_AAZ (tmp_r);

					if (Q[i] != NULL) {
						Q[i] = addPolynomials_AAZ_inp (Q[i], tmp_q, nvar);
					} else {
						Q[i] = deepCopyPolynomial_AAZ (tmp_q);
					}

					if (!isOne_AAZ (tmp_q)) {
						tmp_q = multiplyPolynomials_AAZ_inp (tmp_q, G[i], nvar);
						h = subPolynomials_AAZ_inp (h, tmp_q, nvar);
						freePolynomial_AAZ (tmp_q);
					} else {
						freePolynomial_AAZ (tmp_q);
						h = subPolynomials_AAZ_inp (h, G[i], nvar);
					}

					/* h = tmp_r; */
					/* i = 0; */
					isDivided = 1;
				} else
					i = i + 1;
			}

		if (!isDivided) {
			if (h != NULL && h->size != 0){
				lt = leadingTerm_AAZ (h, nvar);
			} else {
				lt = NULL;
			}

			if (rem == NULL) {
				rem = lt;
			} else {
				if (lt != NULL) {
					rem = addPolynomials_AAZ_inp (rem, lt, nvar);
					freePolynomial_AAZ (lt);
				}
			}
			if (h != NULL && h->size != 0) {
				/* h = subPolynomials_AAZ_inp (h, lt, nvar); */
				subByLeadingTerm_AAZ (&h, nvar);
			}
		}
	}

	*r = rem;
	return;
}

// Note this is a test
// See the recursiveTeriangularSetMDD_AAZ in SMQP_Support_Recursive-AA.h
void recTriangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar)
{
	if (s < 0){
		fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0){
		*r = deepCopyPolynomial_AAZ (f);
		return;
	}

	if (s == 1){  // simple multivariate polynomial division
		if (G[0] != NULL && G[0]->size != 0){
		dividePolynomials_AAZ (f, *G, Q, r, nvar);
		} else {
		*r = deepCopyPolynomial_AAZ (f);
		Q[0] = NULL;
		}
		return;
	}

	int index;
	int i;
	int* orderList = recursiveLoop (s);
	int orderSize = (1 << s) - 1;
	AltArrZ_t* h = deepCopyPolynomial_AAZ (f);
	AltArrZ_t* rem = NULL;
	AltArrZ_t* lt = NULL;
	AltArrZ_t* tmp_q = NULL;
	AltArrZ_t* tmp_r = NULL;

	while (h != NULL && h->size != 0){
		for (i = 0; i < orderSize; i++) {
			index = orderList[i];
			if (h->size != 0 && h != NULL && G[index] != NULL  && G[index]->size != 0
				&& monomialDivideTest (h->elems[0].degs, G[index]->elems[0].degs, nvar)){
				tmp_q = NULL;
				tmp_r = NULL;

				dividePolynomials_AAZ (h, G[index], &tmp_q, &tmp_r, nvar);
				if (Q[index] != NULL) {
					Q[index] = addPolynomials_AAZ_inp (Q[index], tmp_q, nvar);
				}
				else {
					Q[index] = tmp_q;
				}
				freePolynomial_AAZ (h);
				h = tmp_r;
			}
		}

		if (h->size != 0 && h != NULL)
			lt = leadingTerm_AAZ (h, nvar);
		else
			lt = NULL;

		if (rem == NULL)
			rem = lt;
		else {
			if (lt != NULL){
				rem = addPolynomials_AAZ_inp (rem, lt, nvar);
			}
		}
		if (h != NULL && h->size != 0){
			/* h = subPolynomials_AAZ_inp (h, lt, nvar); */
			subByLeadingTerm_AAZ (&h, nvar);
		}
	}

	*r = rem;
	free(orderList);
	return;
}

// Look at the recursiveTriangularSetMDD_AAZ in SMQP_Support_Recursive-AA.h
void triangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar)
{
	if (s < 0){
		fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0){
		*r = deepCopyPolynomial_AAZ (f);
		return;
	}

	if (s == 1){  // simple multivariate polynomial division
		if (G[0] != NULL && G[0]->size != 0){
			dividePolynomials_AAZ (f, *G, Q, r, nvar);
		} else {
			*r = deepCopyPolynomial_AAZ (f);
			Q[0] = NULL;
		}
		return;
	}

	int i; // i= 0 ... (s-1)
	int isDivided;
	AltArrZ_t* h = deepCopyPolynomial_AAZ (f);
	AltArrZ_t* rem = NULL;
	AltArrZ_t* lt = NULL;
	AltArrZ_t* tmp_q;
	AltArrZ_t* tmp_r;

	while (h != NULL && h->size != 0) {

		i = 0;
		isDivided = 0;

		while (i < s && !isDivided)
			{
				if ( h->size != 0 && h!= NULL && G[i] != NULL && G[i]->size != 0  && monomialDivideTest (h->elems[0].degs, G[i]->elems[0].degs, nvar))  {
					tmp_q = NULL;
					tmp_r = NULL;
					divideByLeadingTerms_AAZ (h, G[i], &tmp_q, &tmp_r, nvar);
					freePolynomial_AAZ (tmp_r);

					if (Q[i] != NULL) {
						Q[i] = addPolynomials_AAZ_inp (Q[i], tmp_q, nvar);
					} else {
						Q[i] = deepCopyPolynomial_AAZ (tmp_q);
					}

					if (!isOne_AAZ (tmp_q)) {
						tmp_q = multiplyPolynomials_AAZ_inp (tmp_q, G[i], nvar);
						h = subPolynomials_AAZ_inp (h, tmp_q, nvar);
						freePolynomial_AAZ (tmp_q);
					} else {
						freePolynomial_AAZ (tmp_q);
						h = subPolynomials_AAZ_inp (h, G[i], nvar);
					}
					isDivided = 1;
				} else
					i = i + 1;
			}

		if (!isDivided) {
			if (h->size != 0 && h != NULL){
				lt = leadingTerm_AAZ (h, nvar);
			} else {
				lt = NULL;
			}

			if (rem == NULL) {
				rem = lt;
			} else {
				if (lt != NULL) {
					rem = addPolynomials_AAZ_inp (rem, lt, nvar);
					freePolynomial_AAZ (lt);
				}
			}
			if (h->size != 0 && h!= NULL) {
				subByLeadingTerm_AAZ (&h, nvar);
			}
		}
	}

	*r = rem;
	return;
}

void normalForm_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar) {
	heapMDD_AAZ (f, G, Q, r, s, nvar);
}

void multiDivisorDivision_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar, int type)
{
	if (type == 0){
		heapMDD_AAZ (f, G, Q, r, s, nvar);
		return;
	}

	triangularSetMDD_AAZ (f, G, Q, r, s, nvar);
}

int multiDivisorDivisionVerification_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t* r, AltArrZ_t* hPow, int nSet, int tnvar)
{
	int nullity = 0;
	for (int i = 0; i < nSet; ++i){
		if (G[i] == NULL || G[i]->size == 0){
			fprintf(stderr, "divisor[%d] is NULL!\n", i);
			nullity++;
		}
		if (Q[i] == NULL || Q[i]->size == 0){
			fprintf(stderr, "quotient[%d] is NULL!\n", i);
			nullity++;
		}
	}
	if (r == NULL || r->size == 0){
		fprintf(stderr, "remainder is NULL!\n");
		nullity++;
	}
	if (hPow == NULL || hPow->size == 0){
		fprintf(stderr, "hPow is NULL!\n");
		nullity++;
	}
	if (f == NULL || f->size == 0){
		fprintf(stderr, "dividend is NULL!\n");
		nullity++;
		exit(1);
	}


	AltArrZ_t* sum = NULL;
	for (int i = 0; i < nSet; ++i){
		if (G[i] != NULL && G[i]->size != 0 && Q[i] != NULL && Q[i]->size != 0){
			if (sum == NULL || sum->size == 0){
				sum = multiplyPolynomials_AAZ (G[i], Q[i], tnvar);
			} else {
				sum = addPolynomials_AAZ_inp (sum , multiplyPolynomials_AAZ (G[i], Q[i], tnvar), tnvar);
			}
		}
	}
	if (r != NULL && r->size != 0){
		sum = addPolynomials_AAZ_inp (sum, r, tnvar);
	}

	AltArrZ_t* hf = NULL;
	if (hPow != NULL && hPow->size != 0){
		 hf = multiplyPolynomials_AAZ (f, hPow, tnvar);
	} else {
		hf = deepCopyPolynomial_AAZ (f);
	}

	AltArrZ_t* res = subPolynomials_AAZ (sum, hf, tnvar);
	if (res == NULL || res->size == 0){
		return 1;
	} else {
		fprintf(stderr, "No. Zero-Poly(s): %d\n", nullity);
		printf("The difference: \n");
		printAAZ (res);
		printf("\n");
		return 0;
	}
}

