

#include "IntegerPolynomial/SMZP_Support_Unpacked.h"
#include "IntegerPolynomial/SMZP_Support.h"
#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "IntegerPolynomial/DUZP_Support.h"

#include <math.h>

AltArrZ_t* makePolynomial_AAZ_unpk(int allocSize, int nvar) {
	if (allocSize < 1) {
		return NULL;
	}

	if (nvar == 0) {
		return makePolynomial_AAZ(allocSize, nvar);
	}

	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	newAA->size = 0;
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->unpacked = 1;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*allocSize);
	degree_t* degs = (degree_t*) calloc(allocSize*nvar, sizeof(degree_t));
	newAA->elems[0].degs = (degrees_t) degs;
	return newAA;
}

AltArrZ_t* makeConstPolynomial_AAZ_unpk(int allocSize, int nvar, const mpz_t coef) {
	if (allocSize < 1) {
		return NULL;
	}

	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	newAA->size = 1;
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->unpacked = 1;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*allocSize);
	degree_t* degs = (degree_t*) calloc(allocSize*nvar, sizeof(degree_t));
	newAA->elems[0].degs = (degrees_t) degs;
	mpz_init_set(newAA->elems[0].coef, coef);
	return newAA;
}

//calcMaxDegs is an internal function,. You better know what you're doing.
//That is, aa must already be unpacked.
void calculateMaxDegsList_AAZ_unpk(const AltArrZ_t* aa, degree_t* maxList) {
	if (isZero_AAZ(aa)) {
		if (aa != NULL) {
			for (int i = 0; i < aa->nvar; ++i) {
				maxList[i] = 0;
			}
		}
		return;
	}

	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	degree_t* degs_unpk;
	for (int i = 0; i < nvar; ++i) {
		maxList[i] = 0;
	}
	for (register int i = 0; i < size; ++i) {
		degs_unpk = (degree_t*) elems[i].degs;
		for (register int j = 0; j < nvar; ++j) {
			maxList[j] = (degs_unpk[j]) > (maxList[j]) ? degs_unpk[j] : maxList[j];
		}
	}

}

//calcMaxDegs is an internal function,. You better know what you're doing.
//That is, aa must already be unpacked, and it returns to you a pointer
//disguised as a degrees_t that needs to be free'd!
degrees_t calculateMaxDegs_AAZ_unpk(const AltArrZ_t* aa) {
	if (isZero_AAZ(aa)) {
		if (aa != NULL) {
			return (degrees_t) calloc(aa->nvar, sizeof(degree_t));
		}
		return 0ll;
	}

	degree_t* maxList = (degree_t*) calloc(aa->nvar, sizeof(degree_t));
	calculateMaxDegsList_AAZ_unpk(aa, maxList);
	return (degrees_t) maxList;
}

void unpackExponentVectors_AAZ_inp(AltArrZ_t* aa) {
	if (aa == NULL || aa->alloc == 0 || aa->nvar == 0 || aa->unpacked) {
		//Note: if nvar < 3 then packed actually has more (or equal) bits of precision
		//Maybe we shouldn't unpack.. but more complex operations
		//would have to handle half packed half unpacked and ughh
		return;
	}

	int nvar = aa->nvar;
	int size = aa->size;
	int nexp = nvar * aa->alloc;

	degree_t* unpackedDegs = (degree_t*) calloc(nexp, sizeof(degree_t));
	aa->unpacked = 1;
	if (size == 0) {
		aa->elems->degs = (degrees_t) unpackedDegs;
		return;
	}

	AAZElem_t* elems = aa->elems;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	degree_t deg;

	for (int i = 0; i < size; ++i) {
		for (int k = 0; k < nvar; ++k) {
			deg = GET_NTH_EXP(elems[i].degs, masks[k], sizes[k]);
			unpackedDegs[i*nvar + k] = deg;
		}
		elems[i].degs = (degrees_t) &(unpackedDegs[i*nvar]);
	}
}

int monomialDivideTest_AAZ_unpk(AltArrZ_t* a, int idxa, AltArrZ_t* b, int idxb) {
	if (a == NULL) {
		return 1;
	}
	if (b == NULL) {
		return 0;
	}
	if (!a->unpacked && !b->unpacked) {
		return monomialDivideTest(a->elems[idxa].degs, b->elems[idxb].degs, a->nvar);
	}
	if (a->unpacked && b->unpacked) {
		return monomialDivideTest_unpk(a->elems[idxa].degs, b->elems[idxb].degs, a->nvar);
	} else if (a->unpacked) {
		//b is packed
		degree_t tmp[b->nvar];
		unpackExponentVector(b->elems[idxb].degs, tmp, b->nvar);
		return monomialDivideTest_unpk(a->elems[idxa].degs, (degrees_t) tmp, b->nvar);
	} else {
		//a is packed
		degree_t tmp[a->nvar];
		unpackExponentVector(a->elems[idxa].degs, tmp, a->nvar);
		return monomialDivideTest_unpk((degrees_t) tmp, b->elems[idxb].degs, a->nvar);
	}
}

void tryPackExponentVectors_AAZ_inp(AltArrZ_t* aa) {
	if (aa == NULL || aa->alloc == 0 || !aa->unpacked) {
		return;
	}
	degree_t* maxDegs = (degree_t*) calculateMaxDegs_AAZ_unpk(aa);
	if (isPackableExponentVector(maxDegs, aa->nvar)) {
		packExponentVectors_AAZ_inp(aa);
	}
	free(maxDegs);
}

void packExponentVectors_AAZ_inp(AltArrZ_t* aa) {
	if (aa == NULL || aa->alloc == 0 || !aa->unpacked) {
		return;
	}
	int nvar = aa->nvar;
	int size = aa->size;

	AAZElem_t* elems = aa->elems;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	degree_t* freeDegs = (degree_t*) elems[0].degs;
	degree_t* degs_unpk;

	for (int i = 0; i < size; ++i) {
		degs_unpk = (degree_t*) elems[i].degs;
		elems[i].degs = 0;
		for (int k = 0; k < nvar; ++k) {
			if (degs_unpk[k] > maxExps[k]) {
				fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for packExponentVectors_AAZ_inp at index %d; %d > %lld.\n", k, degs_unpk[k], maxExps[k]);
				exit(1);
			}
			elems[i].degs |= (((degrees_t) degs_unpk[k]) << sizes[k]);
		}
	}

	aa->unpacked = 0;

	free(freeDegs);
}

void clearAndPackExponentVectors_AAZ_inp(AltArrZ_t* aa) {
	if (aa == NULL || aa->alloc == 0) {
		return;
	}

	//free the degs array
	if (aa->unpacked) {
		degree_t* degs = (degree_t*) aa->elems->degs;
		free(degs);
	}

	//reset all exp vecs to 0.
	for (int i = 0; i < aa->size; ++i) {
		aa->elems[i].degs = 0ll;
	}

	aa->unpacked = 0;
}

void printDegs_AAZ_unpk(FILE* fp, degrees_t degs, const char** syms, int nvar) {
	if (degs == 0) {
		fprintf(stderr, "isZeroExponentVector_unpk: exponent pointer is NULL!\n");
		return;
	}

	degree_t* degs_unpk = (degree_t*) degs;
	fprintf(fp, "%s^%d", syms[0], degs_unpk[0]);
	for (int k = 1; k < nvar; ++k) {
		fprintf(fp, "*%s^%d", syms[k], degs_unpk[k]);
	}
}

void printPoly_AAZ_unpk(FILE* fp, const AltArrZ_t* aa, const char** syms, int nvar) {
	if (aa == NULL || aa->size == 0) {
		fprintf(stderr, "0\n");
		return;
	}

	if (!aa->unpacked) {
		printPoly_AAZ(fp, aa, syms, nvar);
		return;
	}

	gmp_fprintf(fp, "%Zd", aa->elems[0].coef);
	if (!isZeroExponentVector_unpk(aa->elems[0].degs, nvar)) {
		fprintf(fp, "*");
		printDegs_AAZ_unpk(fp, aa->elems[0].degs, syms, nvar);
	}
	for (int i = 1; i < AA_SIZE(aa)-1; ++i) {
		if (mpz_sgn(aa->elems[i].coef) > 0) {
			gmp_fprintf(fp, " + %Zd", aa->elems[i].coef);
		} else {
			mpz_neg(aa->elems[i].coef, aa->elems[i].coef);
			gmp_fprintf(fp, " - %Zd", aa->elems[i].coef);
			mpz_neg(aa->elems[i].coef, aa->elems[i].coef);
		}
		if (!isZeroExponentVector_unpk(aa->elems[i].degs, nvar)) {
			fprintf(fp, "*");
			printDegs_AAZ_unpk(fp, aa->elems[i].degs, syms, nvar);
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
		if (!isZeroExponentVector_unpk(aa->elems[aa->size-1].degs, nvar)) {
			fprintf(fp, "*");
			printDegs_AAZ_unpk(fp, aa->elems[aa->size-1].degs, syms, nvar);
		}
	}

	//fprintf(fp, "\n");
}

int isExactlyEqual_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b) {
	if (a == NULL) {
		return (b == NULL);
	}
	if (b == NULL) {
		return 0;
	}

	if (a->size != b->size || a->nvar != b->nvar) {
		return 0;
	}

	if (!a->unpacked && !b->unpacked) {
		return isExactlyEqual_AAZ(a,b);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*)a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*)b);
		unpackedB = 1;
	}

	int cmp = 1;
	int nvar = a->nvar;
	for (int i = 0; i < a->size; ++i) {
		if (!isEqualExponentVectors_unpk(a->elems[i].degs, b->elems[i].degs, nvar)) {
			cmp = 0;
			break;
		}
		if (mpz_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
			cmp = 0;
			break;
		}
	}

	if (unpackedA) {
		packExponentVectors_AAZ_inp((AltArrZ_t*)a);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*)b);
	}

	return cmp;
}

void nonZeroVariables_AAZ_unpk(const AltArrZ_t* aa, int* foundVar) {
    if (aa == NULL || aa->nvar == 0 || isConstant_AAZ(aa)) {
        return;
    }

    if (!aa->unpacked) {
    	nonZeroVariables_AAZ(aa, foundVar);
    	return;
    }

    int nvar = aa->nvar;
    for (int i = 0; i < nvar; ++i) {
    	foundVar[i] = 0;
    }

    degree_t* aDegs = (degree_t*) aa->elems->degs;
    foundVar[0] = aDegs[0] > 0;
    int searchIdx = 1;
    for (int i = 0; i < aa->size && searchIdx < nvar; ++i) {
    	if (aDegs[i*nvar + searchIdx] > 0) {
            foundVar[searchIdx] = 1;
            while (searchIdx < nvar && foundVar[searchIdx] == 1) {
                ++searchIdx;
            }
    	}

        for (int j = searchIdx; j < nvar; ++j) {
            if (aDegs[i*nvar + j] > 0) {
                foundVar[j] = 1;
            }
        }
    }
}

degree_t totalDegree_AAZ_unpk(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0 || isConstant_AAZ(aa)) {
		return 0;
	}

	if (!aa->unpacked) {
		return totalDegree_AAZ(aa);
	}

	int nvar = aa->nvar;
	int size = aa->size;
	degree_t* aDegs = (degree_t*) aa->elems->degs;
    degree_t total = 0, totalMax = 0;
    for (int i = 0; i < size; ++i) {
        total = 0;
        for (int j = 0; j < nvar; ++j) {
            total += aDegs[i*nvar + j];
        }
        if (total > totalMax) {
            totalMax = total;
        }
    }

    return totalMax;
}

degree_t partialDegree_AAZ_unpk(const AltArrZ_t* aa, int k) {
	if (aa == NULL || aa->size == 0 || k >= aa->nvar) {
		return -1;
	}

	if (!aa->unpacked) {
		return partialDegree_AAZ(aa, k);
	}

	int nvar = aa->nvar;
	degree_t dMax = 0;
	degree_t* degs = (degree_t*) aa->elems->degs;
	for (int i = 0; i < aa->size; ++i) {
		dMax = degs[i*nvar + k] > dMax ? degs[i*nvar + k] : dMax;
	}

	return dMax;
}

void partialDegrees_AAZ_unpk(const AltArrZ_t* aa, degree_t* degsList) {
	if (aa == NULL || aa->size == 0 || aa->nvar == 0) {
		return;
	}

	if (!aa->unpacked) {
		partialDegrees_AAZ(aa, degsList);
		return;
	}

	calculateMaxDegsList_AAZ_unpk(aa, degsList);
}

void lowDegrees_AAZ_unpk(const AltArrZ_t* aa, degree_t* __restrict__ degs) {
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

	if (!aa->unpacked) {
		lowDegrees_AAZ(aa, degs);
		return;
	}


	register int size = aa->size;
	register int nvar = aa->nvar;
	degree_t* __restrict__ aDegs = (degree_t*) aa->elems->degs;

	degree_t allZero = 0;
	for (int i = 0; i < nvar; ++i) {
		degs[i] = aDegs[(size-1)*nvar + i];
		allZero |= degs[i];
	}

	degree_t deg;
	//note the && allZero > 0 for fast bail-out;
	for (register int i = size-2; i >= 0 && allZero > 0; --i) {
		allZero = 0;
		for (register int j = 0; j < nvar; ++j) {
			deg = aDegs[i*nvar + j];
			degs[j] =  deg < (degs[j]) ? deg : degs[j];
			allZero |= degs[j];
		}
	}
}

void removeLowDegrees_AAZ_inp_unpk(AltArrZ_t* aa, degree_t* degs) {
	if (isZero_AAZ(aa)) {
		if (aa != NULL && degs != NULL) {
			int nvar = aa->nvar;
			for (int i = 0; i < nvar; ++i) {
				degs[i] = 0;
			}
		}
		return;
	}

	if (!aa->unpacked) {
		lowDegrees_AAZ(aa, degs);
		return;
	}

	degree_t stackDegs[aa->nvar];
	if (degs == NULL) {
		//Note that doing this will not impact calling scope.
		degs = stackDegs;
	}

	register int size = aa->size;
	register int nvar = aa->nvar;
	degree_t* aDegs = (degree_t*) aa->elems->degs;

	degree_t allZero = 0;
	for (int i = 0; i < nvar; ++i) {
		degs[i] = aDegs[(size-1)*nvar + i];
		allZero |= degs[i];
	}

	degree_t deg;
	//note the && allZero > 0 for fast bail-out;
	for (register int i = size-2; i >= 0 && allZero > 0; --i) {
		allZero = 0;
		for (register int j = 0; j < nvar; ++j) {
			deg = aDegs[i*nvar + j];
			degs[j] =  deg < (degs[j]) ? deg : degs[j];
			allZero |= degs[j];
		}
	}

	if (allZero > 0) {
		for (register int i = 0; i < size; ++i) {
			subtractExponentVectors_unpk(aa->elems[i].degs, (degrees_t) degs, aa->elems[i].degs, nvar);
		}
	}
}



degree_t mainDegree_AAZ_unpk(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return -1;
	}

	if (isConstant_AAZ(aa)) {
		return 0;
	}

	if (!aa->unpacked) {
		return mainDegree_AAZ(aa);
	}

    int nvar = aa->nvar;
    degree_t* degs = (degree_t*) aa->elems->degs;
    for (int j = 0; j < nvar; ++j) {
    	if (degs[j] != 0) {
    		return degs[j];
    	}
    }
    return 0;
}

int mainVariable_AAZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || aa->nvar == 0) {
		return -1;
	}

	if (!aa->unpacked) {
		return mainVariable_AAZ(aa);
	}

    int nvar = aa->nvar;
    degree_t* degs = (degree_t*) aa->elems->degs;
    for (int i = 0; i < nvar; ++i) {
        if (degs[i]) {
        	return i;
        }
    }
    return -1;
}

void coefficient_AAZ_unpk(AltArrZ_t* aa, const degree_t* degsList, int nvar, mpz_t retCoef) {
	if (isZero_AAZ(aa)) {
		mpz_set_si(retCoef, 0l);
        return;
    }
    if (nvar == 0 && aa->nvar == 0) {
        mpz_set(retCoef, aa->elems->coef);
    	return;
	}

	if (!aa->unpacked) {
		coefficient_AAZ(aa, degsList, nvar, retCoef);
		return;
	}

	degrees_t degs = (degrees_t) degsList;
	mpz_set_ui(retCoef, 0l);
    for (int i = 0; i < aa->size; ++i) {
        if (isEqualExponentVectors_unpk(degs, aa->elems[i].degs, nvar)) {
        	mpz_set(retCoef, aa->elems[i].coef);
            break;
        }
        if (isLessExponentVectors_unpk(aa->elems[i].degs, degs, nvar)) {
            //due to ordering of terms, none will have degs as monomial.
            break;
        }
    }
}

void setCoefficient_AAZ_unpk(AltArrZ_t* aa, const degree_t* degsList, int nvar, const mpz_t coef) {
    if (aa == NULL || nvar != aa->nvar) {
    	return;
    }

    if (aa->size == 0 || nvar == 0) {
    	if (aa->alloc < 1) {
    		resizePolynomial_AAZ(aa, 1);
    	}
        mpz_set(aa->elems->coef, coef);
    	return;
    }

    if (!aa->unpacked) {
    	setCoefficient_AAZ(aa, degsList, nvar, coef);
    	return;
    }

    degrees_t degs = (degrees_t) degsList;
    for (int i = 0; i < aa->size; ++i) {
        if (isEqualExponentVectors_unpk(aa->elems[i].degs, degs, nvar)) {
            mpz_set(aa->elems[i].coef, coef);
            return;
        }
        if (isLessExponentVectors_unpk(aa->elems[i].degs, degs, nvar)) {
            //shift everything to the right 1 and insert at i.
            if (aa->size >= aa->alloc) {
                resizePolynomial_AAZ_unpk(aa, aa->size + 1);
            }

            degree_t* aaDegs = (degree_t*) aa->elems->degs;
            aa->elems[aa->size].degs = (degrees_t) (aaDegs + aa->size*nvar);
            mpz_init(aa->elems[aa->size].coef);

            for (int j = aa->size; j > i; --j) {
                mpz_swap(aa->elems[j].coef, aa->elems[j-1].coef);
                setExponentVector_unpk(aa->elems[j].degs, aa->elems[j-1].degs, nvar);
            }
            setExponentVector_unpk(aa->elems[i].degs, degs, nvar);
            aa->elems[i].degs = degs;
            mpz_set(aa->elems[i].coef, coef);
            ++(aa->size);
            return;
        }
    }

    //if we get here then we need to insert at the end of the array;
    if (aa->size >= aa->alloc) {
        resizePolynomial_AAZ_unpk(aa, aa->size + 1);
    }
    degree_t* aaDegs = (degree_t*) aa->elems->degs;
    aa->elems[aa->size].degs = (degrees_t) (aaDegs + aa->size*nvar);

    mpz_init(aa->elems[aa->size].coef);
    mpz_set(aa->elems[aa->size].coef, coef);
    setExponentVector_unpk(aa->elems[aa->size].degs, degs, nvar);
    aa->elems[aa->size].degs = degs;
    ++(aa->size);
}

int isEqualWithVariableOrdering_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b, const int* xs, int xsSize) {
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

	if (!a->unpacked && !b->unpacked) {
		return isEqualWithVariableOrdering_AAZ(a, b, xs, xsSize);
	}

	int nvar = a->nvar;
	int bnvar = b->nvar;
	int v = xsSize / 2;

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

    int ret = 1;
    degree_t* aDegs = (degree_t*) a->elems->degs;
    degree_t* bDegs = (degree_t*) b->elems->degs;
    degree_t adeg, bdeg;
    for (int i = 0; i < a->size; ++i){
        if (mpz_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
            ret = 0;
            break;
        }
        for (int j = 0; j < v; ++j) {
            adeg = aDegs[i*nvar + xs[2*j]-1];
            bdeg = bDegs[i*bnvar + xs[2*j+1]-1];
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

    if (unpackedA) {
    	packExponentVectors_AAZ_inp(a);
    }
    if (unpackedB) {
    	packExponentVectors_AAZ_inp(b);
    }

    return ret;
}



void expandNumVars_AAZ_unpk(AltArrZ_t* aa, int newNvar) {
	if (aa == NULL || aa->alloc == 0) {
		return;
	}

	if (newNvar <= aa->nvar) {
		return;
	}

	// if nvar is 0 we still need to procees as follows in order to allocate
	// the degs array properly!
	// if (aa->nvar == 0) {
	// 	aa->nvar = newNvar;
	// 	return;
	// }

	if (!aa->unpacked) {
	// if (newNvar < 3 || !aa->unpacked) {
		//why are we in unpack...
		//okay maybe for testing we want to explicitly unpack less than 3 variables.
		// fprintf(stderr, "improperly in expandNumVars_unpk\n" );
		expandNumVars_AAZ(aa, newNvar);
		return;
	}

	degree_t* newDegs_p = (degree_t*) calloc(aa->alloc*newNvar, sizeof(degree_t));
	degree_t* oldDegs_p = (degree_t*) aa->elems[0].degs;
	aa->elems[0].degs = (degrees_t) newDegs_p;

	int oldNvar = aa->nvar;
	for (int i = 0; i < aa->size; ++i) {
		aa->elems[i].degs = (degrees_t) (newDegs_p + i*newNvar);
		memcpy(newDegs_p + i*newNvar, oldDegs_p + i*oldNvar, sizeof(degree_t)*oldNvar);
	}

	free(oldDegs_p);
	aa->nvar = newNvar;
}

void expandNumVarsLeft_AAZ_unpk(AltArrZ_t* aa, int newNvar) {
	if (aa == NULL || aa->alloc == 0) {
		return;
	}
	if (newNvar <= aa->nvar) {
		return;
	}

	// if nvar is 0 we still need to procees as follows in order to allocate
	// the degs array properly!
	// if (aa->nvar == 0) {
	// 	aa->nvar = newNvar;
	// 	return;
	// }
	if (!aa->unpacked) {
	// if (newNvar < 3 || !aa->unpacked) {
		//why are we in unpack...
		//okay maybe for testing we want to explicitly unpack less than 3 variables.
		// fprintf(stderr, "improperly in expandNumVarsLeft_unpk\n" );
		expandNumVarsLeft_AAZ(aa, newNvar);
		return;
	}

	degree_t* newDegs_p = (degree_t*) calloc(aa->alloc*newNvar, sizeof(degree_t));
	degree_t* oldDegs_p = (degree_t*) aa->elems[0].degs;
	aa->elems[0].degs = (degrees_t) newDegs_p;


	int oldNvar = aa->nvar;
	int diff = newNvar - aa->nvar;
	for (int i = 0; i < aa->size; ++i) {
		aa->elems[i].degs = (degrees_t) (newDegs_p + i*newNvar);
		memcpy(newDegs_p + i*newNvar + diff, oldDegs_p + i*oldNvar, sizeof(degree_t)*oldNvar);
	}

	free(oldDegs_p);
	aa->nvar = newNvar;
}

void shrinkNumVarsAtIdx_AAZ_unpk(AltArrZ_t* aa, int idx) {
	if (aa == NULL || aa->size < 1 || aa->nvar < 1) {
		return;
	}
	if (aa->nvar == 1) {
		free((degree_t*) aa->elems->degs);
		aa->elems->degs = 0;
		aa->nvar = 0;
		return;
	}

	if (!aa->unpacked) {
		shrinkNumVarsAtIdx_AAZ(aa, idx);
		return;
	}

	int nvar = aa->nvar;
	int newNvar = aa->nvar - 1;

	degree_t* newDegs = (degree_t*) malloc(sizeof(degree_t)*newNvar*aa->alloc);
	degree_t* curDegs = (degree_t*) aa->elems->degs;
	AAZElem_t* elems = aa->elems;

	for (int i = 0; i < aa->size; ++i ) {
		if (idx > 0) {
			memcpy(newDegs + i*newNvar, curDegs + i*nvar, sizeof(degree_t)*idx);
		}
		memcpy(newDegs + i*newNvar + idx, curDegs + i*nvar + idx + 1, sizeof(degree_t)*(nvar-idx-1));
		elems[i].degs = (degrees_t) (newDegs + i*newNvar);
	}

	aa->nvar = newNvar;

	free(curDegs);

	tryPackExponentVectors_AAZ_inp(aa);
}

void shrinkAndReorderVars_AAZ_unpk(AltArrZ_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->nvar < 1) {
		return;
	}

	if (varmapSize > aa->nvar) {
		return;
	}

	if (!aa->unpacked) {
		shrinkAndReorderVars_AAZ(aa, varMap, varmapSize);
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

	if (newNvar == 0) {
		for (int i = 1; i < aa->size; ++i) {
			mpz_clear(aa->elems[i].coef);
		}
		free((degree_t*) aa->elems->degs);
		aa->elems->degs = 0ll;
		aa->size = 1;
		aa->nvar = 0;
		return;
	}

	int nvar = aa->nvar;
	AAZElem_t* elems = aa->elems;
	degree_t* newDegs = (degree_t*) malloc(sizeof(degree_t)*newNvar*aa->alloc);
	degree_t* oldDegs = (degree_t*) elems->degs;
	for (int i = 0; i < aa->size; ++i) {
		for (int j = 0; j < varmapSize; ++j) {
			if (varMap[j] < 0) {
				continue;
			}
			newDegs[i*newNvar + varMap[j]] = oldDegs[i*nvar + j];
		}
		elems[i].degs = (degrees_t) (newDegs + i*newNvar);
	}

	aa->nvar = newNvar;

	free(oldDegs);

	if (needSort) {
		mergeSortPolynomial_AAZ(aa);
	}

	tryPackExponentVectors_AAZ_inp(aa);
}

//     degs[varmap[i]] = curDegs[i];
void reorderVars_AAZ_unpk(AltArrZ_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (!aa->unpacked) {
		reorderVars_AAZ(aa, varMap, varmapSize);
		return;
	}

	int nvar = aa->nvar;

	degree_t tmpDegs[nvar];
	degree_t zeroDegs[nvar];
	for (int i = 0; i < nvar; ++i) {
		zeroDegs[i] = 0;
	}
	degree_t* degs = (degree_t*) aa->elems->degs;

	for (int i = 0; i < aa->size; ++i) {
		memcpy(tmpDegs, degs + i*nvar, sizeof(degree_t)*nvar);
		memcpy(degs + i*nvar, zeroDegs, sizeof(degree_t)*nvar);

		for (int j = 0; j < varmapSize; ++j) {
			degs[i*nvar + varMap[j]] = tmpDegs[j];
		}
	}

	mergeSortPolynomial_AAZ(aa);
}

void setDegrees_AAZ_inp_unpk(AltArrZ_t* aa, int idx, const degree_t* degsList, int nvar) {
	if (aa == NULL || idx >= aa->alloc) {
		return;
	}

	if (!aa->unpacked) {
		setDegrees_AAZ_inp(aa, idx, degsList, nvar);
		return;
	}

	degree_t* degs = (degree_t*) aa->elems[idx].degs;
	for (int k = 0; k < nvar; ++k) {
		degs[k] = degsList[k];
	}
}

void setExponentTerm_AAZ_inp_unpk(AltArrZ_t* aa, int idx, degree_t deg, int k) {
	if (aa==NULL || idx >= aa->alloc) {
		return;
	}

	if (!aa->unpacked) {
		setExponentTerm_AAZ_inp(aa, idx, deg, k);
		return;
	}

	degree_t* degs = (degree_t*) aa->elems[idx].degs;
	degs[k] = deg;
}

void resizePolynomial_AAZ_unpk(AltArrZ_t* aa, int allocSize) {
	if (aa == NULL) {
		return;
	}
	if (!aa->unpacked) {
		resizePolynomial_AAZ(aa, allocSize);
	}

	if (aa->alloc == allocSize) {
		return;
	}

	if (allocSize < aa->size) {
		for (int i = aa->size; i >= allocSize; --i) {
			mpz_clear(aa->elems[i].coef);
		}
		aa->size = allocSize;
	}

	if (allocSize <= 0) {
		//just don't update allocsize like why are doing this non-sense resize?
		return;
	}

	int oldAlloc = aa->alloc;

	aa->elems = (AAZElem_t*) realloc(aa->elems, sizeof(AAZElem_t)*allocSize);
	aa->alloc = allocSize;

	degree_t* degsHead = (degree_t*) aa->elems[0].degs;
	degsHead = realloc(degsHead, sizeof(degree_t)*aa->nvar*allocSize);
	aa->elems[0].degs = (degrees_t) degsHead;

	//do we what calloc would have done;
	if (allocSize > oldAlloc) {
		memset(degsHead + oldAlloc*aa->nvar, 0, sizeof(degree_t)*aa->nvar*(allocSize-oldAlloc));
	}

	int size = aa->size;
	int nvar = aa->nvar;
	AAZElem_t* elems = aa->elems;
	for (int i = 0; i < size; ++i) {
		elems[i].degs = (degrees_t) (degsHead + i*nvar);
	}
}

AltArrZ_t* deepCopyPolynomial_AAZ_unpk(const AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		AltArrZ_t* ret = deepCopyPolynomial_AAZ(aa);
		unpackExponentVectors_AAZ_inp(ret);
		return ret;
	}

	if (aa->nvar == 0) {
		return makeConstPolynomial_AAZ(aa->size, 0, aa->elems->coef);
	}

	int nvar = aa->nvar;
	int size = aa->size;
	int nexp = nvar * aa->alloc;

	degree_t* unpackedDegs = (degree_t*) malloc(sizeof(degree_t) * nexp);
	degree_t* original_degs = (degree_t*) aa->elems[0].degs;
	memcpy(unpackedDegs, original_degs, sizeof(degree_t)*nexp);

	AltArrZ_t* ret = makePolynomial_AAZ(aa->alloc, nvar);
	ret->size = size;

	AAZElem_t* elems = aa->elems;
	AAZElem_t* retelems = ret->elems;
	for (int i = 0; i < size; ++i) {
		mpz_init(retelems[i].coef);
		mpz_set(retelems[i].coef, elems[i].coef);
		retelems[i].degs = (degrees_t) &(unpackedDegs[i*nvar]);
	}

	ret->unpacked = 1;

	return ret;
}

void deepCopyPolynomial_AAZ_inp_unpk(const AltArrZ_t* aa, AltArrZ_t** bb) {
	if (bb == NULL) {
		return;
	}

	if (isZero_AAZ(aa)) {
		if (!isZero_AAZ(*bb)) {
			freePolynomial_AAZ(*bb);
			*bb = NULL;
		}
		return;
	}

	if (!aa->unpacked) {
		if (*bb == NULL) {
			deepCopyPolynomial_AAZ_inp(aa, bb);
			return;
		}
		if ((*bb)->unpacked) {
			free( (degree_t*) (*bb)->elems[0].degs );
			(*bb)->unpacked = 0;
		}

		deepCopyPolynomial_AAZ_inp(aa, bb);
		return;
	}

	AltArrZ_t* b = *bb;
	if (b == NULL) {
		b = deepCopyPolynomial_AAZ_unpk(aa);
		*bb = b;
		return;
	} else if (b->alloc < aa->alloc) {
		resizePolynomial_AAZ_unpk(b, aa->alloc);
	}

	degree_t* unpackedDegs;
	int nvar = aa->nvar;
	if (!b->unpacked) {
		unpackedDegs = (degree_t*) malloc(sizeof(degree_t) * nvar * b->alloc);
		b->unpacked = 1;
	} else {
		unpackedDegs = (degree_t*) b->elems->degs;
		if (b->nvar != nvar) {
			free(unpackedDegs);
			unpackedDegs = (degree_t*) malloc(sizeof(degree_t) * nvar * b->alloc);
			b->nvar = nvar;
		}
	}
	memcpy(unpackedDegs, (degree_t*) aa->elems->degs, sizeof(degree_t)*aa->nvar*aa->size);

	polysize_t i;
	if (b->size < aa->size) {
		for (i = 0; i < b->size; ++i) {
			mpz_set(b->elems[i].coef, aa->elems[i].coef);
			b->elems[i].degs = (degrees_t) (unpackedDegs + (i*nvar));
		}
		for ( ; i < aa->size; ++i) {
			mpz_init_set(b->elems[i].coef, aa->elems[i].coef);
			b->elems[i].degs = (degrees_t) (unpackedDegs + (i*nvar));
		}
	} else {
		for (i = 0; i < aa->size; ++i) {
			mpz_set(b->elems[i].coef, aa->elems[i].coef);
			b->elems[i].degs = (degrees_t) (unpackedDegs + (i*nvar));
		}
		for ( ; i < b->size; ++i) {
			mpz_clear(b->elems[i].coef);
		}
	}

	b->size = aa->size;
}


AltArrZDegList_t* deepCopyPolynomial_AAZDegListFromAA_unpk(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return deepCopyPolynomial_AAZDegListFromAA(aa);
	}

	int nvar = aa->nvar;

	AltArrZDegList_t* ret = makePolynomial_AAZDL(aa->alloc, nvar);
	AAZElem_DegList_t* rElems = ret->elems;
	AAZElem_t* elems = aa->elems;
	degree_t* degs = (degree_t*) elems->degs;

	for (int i = 0; i < aa->size; ++i) {
		mpz_init(rElems[i].coef);
		mpz_set(rElems[i].coef, elems[i].coef);
		rElems[i].degs = malloc(sizeof(degree_t)*nvar);
		for (int j = 0; j < nvar; ++j) {
			rElems[i].degs[j] = degs[i*nvar + j];
		}
	}

	ret->size = aa->size;

	return ret;
}

AltArr_t* deepCopyPolynomial_AAFromAAZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return deepCopyPolynomial_AAFromAAZ(aa);
	}

	if (aa->nvar == 0) {
		mpq_t c;
		mpq_init(c);
		mpz_set(mpq_numref(c), aa->elems->coef);
		AltArr_t* ret = makeConstPolynomial_AA(aa->size, 0, c);
		mpq_clear(c);
		return ret;
	}

	int nvar = aa->nvar;
	int size = aa->size;
	int nexp = nvar * size;

	degree_t* unpackedDegs = (degree_t*) malloc(sizeof(degree_t) * nexp);
	degree_t* original_degs = (degree_t*) aa->elems[0].degs;
	memcpy(unpackedDegs, original_degs, sizeof(degree_t)*nexp);

	AltArr_t* ret = makePolynomial_AA(aa->alloc, nvar);
	ret->size = size;

	AAZElem_t* elems = aa->elems;
	AAElem_t* retelems = ret->elems;
	for (int i = 0; i < size; ++i) {
		mpq_init(retelems[i].coef);
		mpz_set(mpq_numref(retelems[i].coef), elems[i].coef);
		retelems[i].degs = (degrees_t) &(unpackedDegs[i*nvar]);
	}

	ret->unpacked = 1;

	return ret;

}

AltArrZ_t* deepCopyPolynomial_AAZFromAA_unpk(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return deepCopyPolynomial_AAZFromAA(aa);
	}

	if (aa->nvar == 0) {
		if (mpz_cmp_si(mpq_denref(aa->elems->coef), 1l) != 0) {
			fprintf(stderr, "SMZP ERROR: Failed to convert a rational number polynomial to an integer one.\n");
			exit(1);
		}
		return makeConstPolynomial_AAZ(aa->size, 0, mpq_numref(aa->elems->coef));
	}

	int nvar = aa->nvar;
	int size = aa->size;
	int nexp = nvar * aa->alloc;

	degree_t* unpackedDegs = (degree_t*) malloc(sizeof(degree_t) * nexp);
	degree_t* original_degs = (degree_t*) aa->elems[0].degs;
	memcpy(unpackedDegs, original_degs, sizeof(degree_t)*nexp);

	AltArrZ_t* ret = makePolynomial_AAZ(aa->alloc, nvar);
	ret->size = size;

	AAElem_t* elems = aa->elems;
	AAZElem_t* retelems = ret->elems;
	for (int i = 0; i < size; ++i) {
		if (mpz_cmp_si(mpq_denref(elems[i].coef), 1l) != 0) {
			fprintf(stderr, "SMZP ERROR: Failed to convert a rational number polynomial to an integer one.\n");
			exit(1);
		}
		mpz_init(retelems[i].coef);
		mpz_set(retelems[i].coef, mpq_numref(elems[i].coef));
		retelems[i].degs = (degrees_t) &(unpackedDegs[i*nvar]);
	}

	ret->unpacked = 1;

	return ret;

}


AltArrZ_t* sortPolynomial_AAZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}
	if (!aa->unpacked) {
		return sortPolynomial_AAZ(aa);
	}

	AAZElem_t* elems = aa->elems;
	int size = AA_SIZE(aa);
	int nvar = aa->nvar;

	degree_t swapDegs[nvar];
	for (int i = 1; i < size; ++i) {
		for (int j = i; j > 0 && compareExponentVectors_unpk(elems[j-1].degs, elems[j].degs, nvar) < 0; --j) {
			mpz_swap(elems[j-1].coef, elems[j].coef);
			setExponentVector_unpk((degrees_t) swapDegs, elems[j-1].degs, nvar);
			setExponentVector_unpk(elems[j-1].degs, elems[j].degs, nvar);
			setExponentVector_unpk(elems[j].degs,(degrees_t)  swapDegs, nvar);
			// swapDegs = elems[j-1].degs;
			// elems[j-1].degs = elems[j].degs;
			// elems[j].degs = swapDegs;
		}
	}

	condensePolynomial_AAZ_unpk(aa);

	return aa;
}


static void mergeAAElems_unpk(AAZElem_t* __restrict__ a, AAZElem_t* __restrict__ endA, AAZElem_t* __restrict__ b , AAZElem_t* __restrict__ endB, AAZElem_t* __restrict__ sorted, int nvar) {

	int i = 0;
	while (a < endA && b < endB) {
		if (isGreaterExponentVectors_unpk(a->degs, b->degs, nvar)) {
			sorted[i].coef[0] = a->coef[0];
			// sorted[i].coef->_mp_num = a->coef->_mp_num;
			// sorted[i].coef->_mp_den = a->coef->_mp_den;
			memcpy((degree_t*) sorted[i].degs, (degree_t*) a->degs, nvar*sizeof(degree_t));
			// sorted[i] = *a;
			++a;
		} else {
			sorted[i].coef[0] = b->coef[0];
			// sorted[i].coef = b->coef;
			// sorted[i].coef->_mp_num = b->coef->_mp_num;
			// sorted[i].coef->_mp_den = b->coef->_mp_den;
			memcpy((degree_t*) sorted[i].degs, (degree_t*) b->degs, nvar*sizeof(degree_t));
			// sorted[i] = *b;
			++b;
		}
		++i;
	}

	while (a < endA) {
		sorted[i].coef[0] = a->coef[0];
		// sorted[i].coef = a->coef;
		// sorted[i].coef->_mp_num = a->coef->_mp_num;
		// sorted[i].coef->_mp_den = a->coef->_mp_den;
		memcpy((degree_t*) sorted[i].degs, (degree_t*) a->degs, nvar*sizeof(degree_t));
		// sorted[i] = *a;
		++a;
		++i;
	}

	while (b < endB) {
		sorted[i].coef[0] = b->coef[0];
		// sorted[i].coef = b->coef;
		// sorted[i].coef->_mp_num = b->coef->_mp_num;
		// sorted[i].coef->_mp_den = b->coef->_mp_den;
		memcpy((degree_t*) sorted[i].degs, (degree_t*) b->degs, nvar*sizeof(degree_t));
		// sorted[i] = *b;
		++b;
		++i;
	}
}

void mergeSortPolynomial_AAZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (!aa->unpacked) {
		mergeSortPolynomial_AAZ(aa);
		return;
	}

	if (aa->size < 8) {
		sortPolynomial_AAZ_unpk(aa);
		return;
	}

	int nvar = aa->nvar;
	int size = aa->size;
	AAZElem_t* tempElems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*aa->alloc);
	degree_t* tempDegs = (degree_t*) malloc(sizeof(degree_t)*aa->alloc*nvar);
	for (int i = 0; i < size; ++i) {
		tempElems[i].degs = (degrees_t) &(tempDegs[i*nvar]);
	}
	AAZElem_t* elems = aa->elems;
	int end1 = 0;
	int end2 = 0;
	for (int window = 1; window < size; window <<= 1) {
		//merge elems[i]:elems[i+window] with elems[i+window]:elems[i+2*window]
		for (int i = 0; i < size; i += 2*window) {
			end1 = i + window < size ? i + window : size;
			end2 = i + 2*window < size ? i + 2*window : size;
			mergeAAElems_unpk(elems+i, elems+end1, elems+end1, elems+end2, tempElems+i, nvar);
		}

		AAZElem_t* temp = tempElems;
		tempElems = elems;
		elems = temp;
	}

	aa->elems = elems;
	free((degree_t*) tempElems[0].degs);
	free(tempElems);

	condensePolynomial_AAZ_unpk(aa);
}

void condensePolynomial_AAZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (AA_SIZE(aa) < 1) {
		return;
	}

	if (!aa->unpacked) {
		condensePolynomial_AAZ(aa);
		return;
	}

	int size = AA_SIZE(aa);
	int nvar = aa->nvar;
	AAZElem_t* elems = aa->elems;
	int insertIdx = 0;
	int compareIdx = 1;
	while (compareIdx < size) {
		if (compareExponentVectors_unpk(elems[insertIdx].degs, elems[compareIdx].degs, nvar) == 0) {
			mpz_add(elems[insertIdx].coef, elems[insertIdx].coef, elems[compareIdx].coef);
		} else if (mpz_sgn(elems[insertIdx].coef) == 0) {
			//we are about to move to next degree, but the coef is 0, so we
			//need to remove it. Do so by setting insertIdx's degree to the
			//comare idx so it can get overwritten by next compare and add.
			memcpy((degree_t*) elems[insertIdx].degs, (degree_t*) elems[compareIdx].degs, sizeof(degree_t)*nvar);
			// elems[insertIdx].degs = elems[compareIdx].degs;
			--compareIdx;
		} else if(compareIdx - insertIdx > 1) {
			++insertIdx;
			memcpy((degree_t*) elems[insertIdx].degs, (degree_t*) elems[compareIdx].degs, sizeof(degree_t)*nvar);
			// elems[insertIdx].degs = elems[compareIdx].degs;
			mpz_swap(elems[insertIdx].coef, elems[compareIdx].coef);
		} else {
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

void removeZeroTerms_AAZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || aa->size <= 1) {
		return;
	}

	if (!aa->unpacked) {
		removeZeroTerms_AAZ(aa);
		return;
	}

	int size = AA_SIZE(aa);
	AAZElem_t* elems = aa->elems;
	int curIdx = 0;
	int nvar = aa->nvar;
	for (int i = 0; i < size; ++i) {
		if (mpz_cmp_ui(elems[i].coef, 0ul) != 0) {
			if (i != curIdx) {
				mpz_swap(elems[curIdx].coef, elems[i].coef);
				memcpy((degree_t*) elems[curIdx].degs, (degree_t*) elems[i].degs, sizeof(degree_t)*nvar);
			}
			++curIdx;
		}
	}

	//keep one zero around if aa is not null and indentically 0.
	if (curIdx == 0) {
		degree_t* degs = (degree_t*) elems[0].degs;
		for (int k = 0; k < nvar; ++k) {
			degs[k] = 0;
		}
		++curIdx;
	}

	for (int i = curIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
	}

	aa->size = curIdx;
}


void evalPolyToVal_AAZ_unpk(const AltArrZ_t* aa, const mpz_t* vals, short nvar, mpz_t res) {
	if (aa == NULL || nvar != aa->nvar) {
		mpz_set_ui(res, 0ul);
	}

	if (nvar == 0 || aa->nvar == 0) {
		mpz_set(res, aa->elems->coef);
		return;
	}

	if (!aa->unpacked) {
		return evalPolyToVal_AAZ(aa, vals, nvar, res);
	}

	mpz_t* valList[nvar];
	int valListSize[nvar];

	degree_t* maxDegs = (degree_t*) calculateMaxDegs_AAZ_unpk(aa);
	for (int j = 0; j < nvar; ++j) {
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(maxDegs[j]+1));
		valListSize[j] = maxDegs[j]+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < maxDegs[j]+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	mpz_set_ui(res, 0ul);
	mpz_t acc;
	mpz_init(acc);

	int size = aa->size;
	degree_t* degs_unpk;
	for (int i = 0; i < size; ++i) {
		mpz_set(acc, aa->elems[i].coef);
		degs_unpk = (degree_t*) aa->elems[i].degs;
		for (int j = 0; j < nvar; ++j) {
			mpz_mul(acc, acc, valList[j][degs_unpk[j]]);
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
	free(maxDegs);
}

AltArrZ_t* evaluatePoly_AAZ_unpk(const AltArrZ_t* aa, const int* active, const mpz_t* vals, short nvar) {
	if (aa == NULL) {
		return NULL;
	}

	if (nvar == 0 || aa->nvar == 0) {
		return deepCopyPolynomial_AAZ(aa);
	}

	if(!aa->unpacked) {
		return evaluatePoly_AAZ(aa, active, vals, nvar);
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



	mpz_t* valList[nvar];
	int valListSize[nvar];

	degree_t* maxDegs = (degree_t*) calculateMaxDegs_AAZ_unpk(aa);
	for (int j = 0; j < nvar; ++j) {
		if (!active[j]) {
			valList[j] = NULL;
			valListSize[j] = 0;
			continue;
		}
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(maxDegs[j]+1));
		valListSize[j] = maxDegs[j]+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < maxDegs[j]+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}
	free(maxDegs);

	degree_t* degsArr = (degree_t*) aa->elems[0].degs;

	AltArrZ_t* res = makePolynomial_AAZ_unpk(aa->size, newNvar);
	AAZElem_t* resElems = res->elems;
	degree_t* newDegsArr = (degree_t*) resElems[0].degs;
	res->size = aa->size;

	int size = aa->size;
	int k = 0;
	for (int i = 0; i < size; ++i) {
		mpz_init(resElems[i].coef);
		mpz_set(resElems[i].coef, aa->elems[i].coef);
		resElems[i].degs = (degrees_t) &(newDegsArr[i*newNvar]);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = degsArr[i*nvar + j];
			if (valList[j] == NULL) {
				newDegsArr[i*newNvar + k] = deg;
				// newDegs |= deg << newSizes[k];
				++k;
			} else {
				mpz_mul(resElems[i].coef, resElems[i].coef, valList[j][deg]);
			}
		}
		k = 0;
	}

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}

	canonicalizePolynomial_AAZ_unpk(res);

	tryPackExponentVectors_AAZ_inp(res);

	return res;
}


AltArrZ_t* mainLShiftPolynomial_AAZ_unpk (AltArrZ_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (n < 1){
		return deepCopyPolynomial_AAZ (aa);
    }

    if (!aa->unpacked) {
    	return mainLShiftPolynomial_AAZ(aa, n);
    }

    AltArrZ_t* shifted = deepCopyPolynomial_AAZ_unpk(aa);
    return mainLShiftPolynomial_AAZ_inp_unpk(shifted, n);
}


AltArrZ_t* mainLShiftPolynomial_AAZ_inp_unpk (AltArrZ_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (n < 1){
		return aa;
    }

    if (!aa->unpacked) {
    	return mainLShiftPolynomial_AAZ_inp(aa, n);
    }

    int nvar = aa->nvar;
    degree_t* degs = (degree_t*) aa->elems->degs;
    for (int i = 0; i < aa->size; ++i){
    	degs[i*nvar] += n;
    }

    return aa;
}

AltArrZ_t* leadingTerm_AAZ_unpk (AltArrZ_t* aa, int nvar) {
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	if (!aa->unpacked) {
		return leadingTerm_AAZ(aa, nvar);
	}

	AltArrZ_t* lt = makePolynomial_AAZ_unpk(1, aa->nvar);
	mpz_init(lt->elems->coef);
	mpz_set(lt->elems->coef, aa->elems->coef);
	setExponentVector_unpk(lt->elems->degs, aa->elems->degs, aa->nvar);
	lt->size = 1;
	return lt;
}

int leadingVariable_AAZ_unpk (const AltArrZ_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return -2;
	}

	if (!aa->unpacked) {
		return leadingVariable_AAZ(aa);
	}

	degree_t* degs = (degree_t*) aa->elems->degs;
	for (int i = 0; i < aa->nvar; ++i){
		if (degs[i] > 0) {
			return i;
		}
	}

	return -1;
}

int mainLeadingDegree_AAZ_unpk (AltArrZ_t* aa)
{
    if (aa == NULL || aa->size == 0){
		return 0;
    }

    if (!aa->unpacked) {
    	return mainLeadingDegree_AAZ(aa);
    }

    degree_t* degArr = (degree_t*) aa->elems[0].degs;
    return degArr[0];
}

AltArrZ_t* mainLeadingCoefficient_AAZ_unpk (const AltArrZ_t* aa)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }

    if (!aa->unpacked) {
    	return mainLeadingCoefficient_AAZ(aa);
    }

    AAZElem_t* elems = aa->elems;
    degree_t* degs = (degree_t*) elems[0].degs;
    degree_t mvarDeg = degs[0];
    degree_t curDeg = mvarDeg;

    int nvar = aa->nvar;
   	int size = aa->size;
    int polyElemsAlloc = 10;
    int polyElemsSize = 0;
    AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);
    degree_t* polyDegs = (degree_t*) malloc(sizeof(degree_t) * polyElemsAlloc * nvar);

    for (int i = 0; (mvarDeg == curDeg) && i < size; ++i){
		if (polyElemsSize + 1 > polyElemsAlloc){
		    polyElemsAlloc += 10;
		    polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t) * polyElemsAlloc);
		    polyDegs = (degree_t*) realloc(polyDegs, sizeof(degree_t) * polyElemsAlloc * nvar);
		}

		mpz_init (polyElems[i].coef);
		mpz_set (polyElems[i].coef, elems[i].coef);
		setExponentVector_unpk((degrees_t) (polyDegs + i*nvar), elems[i].degs, nvar);
		polyDegs[i*nvar] = 0;
		++polyElemsSize;

		if (i+1 < size){
		    curDeg = degs[(i+1)*nvar];
		} else {
		    curDeg = -1;
		}
    }

    AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
    poly->alloc = polyElemsAlloc;
    poly->size = polyElemsSize;
    poly->elems = polyElems;
    poly->unpacked = 1;
    poly->nvar = aa->nvar;

    size = poly->size;
    for (int i = 0; i < size; ++i) {
    	polyElems[i].degs = (degrees_t) (polyDegs + i*nvar);
    }

    return poly;
}

AltArrZ_t* mainCoefficientAtIdx_AAZ_unpk (AltArrZ_t* aa, int e)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (e < 0){
		return NULL;
    }
    if (!aa->unpacked) {
    	return mainCoefficientAtIdx_AAZ(aa, e);
    }

    AAZElem_t* elems = aa->elems;
    degree_t* degs = (degree_t*) elems->degs;
    degree_t mvarDeg = e;
    degree_t curDeg = degs[0];
    int begin = -1;

    if (mvarDeg > curDeg){
        return NULL;
    }

    int nvar = aa->nvar;
    int size = aa->size;
    int polyElemsAlloc = 10;
    int polyElemsSize = 0;
    AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);
	degree_t* polyDegs = (degree_t*) malloc(sizeof(degree_t) * nvar * polyElemsAlloc);

    for (int i = 0; i < size; ++i){
		if (mvarDeg == curDeg){
		    if (polyElemsSize + 1 > polyElemsAlloc) {
				polyElemsAlloc += 10;
				polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t)*polyElemsAlloc);
		    	polyDegs = (degree_t*) realloc(polyDegs, sizeof(degree_t) * polyElemsAlloc * nvar);
		    }

		    mpz_init (polyElems[polyElemsSize].coef);
		    mpz_set (polyElems[polyElemsSize].coef, elems[i].coef);
		    setExponentVector_unpk((degrees_t) (polyDegs + polyElemsSize*nvar), elems[i].degs, nvar);
		    polyDegs[polyElemsSize*nvar] = 0;
		    ++polyElemsSize;

		    if (begin == -1) {
				begin = i;
		    }
		}

		if (begin != -1 && mvarDeg != curDeg){
		    break;
		}

		if (i+1 < size){
		    curDeg = degs[(i+1)*nvar];
		} else {
		    curDeg = -1;
		}
	}

    AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
    poly->alloc = polyElemsAlloc;
    poly->size = polyElemsSize;
    poly->elems = polyElems;
    poly->unpacked = 1;
    poly->nvar = aa->nvar;

    size = poly->size;
    for (int i = 0; i < size; ++i) {
    	polyElems[i].degs = (degrees_t) (polyDegs + i*nvar);
    }

    return poly;
}

AltArrZ_t* maxPolynomials_AAZ_unpk (AltArrZ_t* a, AltArrZ_t* b) {
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ(a);
	}

	if (a->nvar != b->nvar) {
		fprintf(stderr, "maxPolynomials_AA inputs are not from the same ring\n");
		exit(1);
	}

	if (!a->unpacked && !b->unpacked) {
		return maxPolynomials_AAZ(a,b);
	}

	int cmp;
	if (a->unpacked && b->unpacked) {
		cmp = compareExponentVectors_unpk(a->elems->degs, b->elems->degs, a->nvar);
	} else if (a->unpacked) {
		//b is packed
		degree_t tmp[b->nvar];
		unpackExponentVector(b->elems->degs, tmp, b->nvar);
		cmp = compareExponentVectors_unpk(a->elems->degs, (degrees_t) tmp, b->nvar);
	} else {
		//a is packed
		degree_t tmp[a->nvar];
		unpackExponentVector(a->elems->degs, tmp, a->nvar);
		cmp = compareExponentVectors_unpk((degrees_t) tmp, b->elems->degs, a->nvar);
	}

	if (cmp > 0){
		return deepCopyPolynomial_AAZ(a);
	} else if (cmp < 0){
		return deepCopyPolynomial_AAZ(b);
	} else {
		if (mpz_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return deepCopyPolynomial_AAZ(a);
		} else {
			return deepCopyPolynomial_AAZ(b);
		}
	}

	return NULL;
}

AltArrZ_t* maxPolynomials_AAZ_inp_unpk (AltArrZ_t* a, AltArrZ_t* b) {
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return a;
	}

	if (!a->unpacked && !b->unpacked) {
		return maxPolynomials_AAZ(a,b);
	}

	int cmp;
	if (a->unpacked && b->unpacked) {
		cmp = compareExponentVectors_unpk(a->elems->degs, b->elems->degs, a->nvar);
	} else if (a->unpacked) {
		//b is packed
		degree_t tmp[b->nvar];
		unpackExponentVector(b->elems->degs, tmp, b->nvar);
		cmp = compareExponentVectors_unpk(a->elems->degs, (degrees_t) tmp, b->nvar);
	} else {
		//a is packed
		degree_t tmp[a->nvar];
		unpackExponentVector(a->elems->degs, tmp, a->nvar);
		cmp = compareExponentVectors_unpk((degrees_t) tmp, b->elems->degs, a->nvar);
	}

	if (cmp > 0){
		return a;
	} else if (cmp < 0){
		return deepCopyPolynomial_AAZ_unpk(b);
	} else {
		if (mpz_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return a;
		} else {
			return deepCopyPolynomial_AAZ_unpk(b);
		}
	}
	return NULL;
}

void addInteger_AAZ_inp_unpk(AltArrZ_t* aa, const mpz_t coef) {
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

	if (!aa->unpacked) {
		addInteger_AAZ_inp(aa, coef);
		return;
	}

	int nvar = aa->nvar;
    if (isZeroExponentVector_unpk(aa->elems[aa->size-1].degs, nvar)) {
        mpz_add(aa->elems[aa->size-1].coef, aa->elems[aa->size-1].coef, coef);
        return;
    }

    if (aa->size >= aa->alloc) {
        resizePolynomial_AAZ_unpk(aa, aa->alloc+10);
    }

    mpz_init(aa->elems[aa->size].coef);
    mpz_set(aa->elems[aa->size].coef, coef);

    degree_t* aDegs = (degree_t*) aa->elems->degs;
    aa->elems[aa->size].degs = (degrees_t) (aDegs + aa->size*nvar);
	aDegs = aDegs + aa->size*nvar;
    for (int i = 0; i < nvar; ++i) {
    	aDegs[i] = 0;
    }
    ++(aa->size);
}

AltArrZ_t* addPolynomials_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return addPolynomials_AAZ(a, b, nvar);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ_unpk(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// mpz_t ccoef;
	// mpz_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			mpz_init(cElems[k].coef);
			mpz_add(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
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
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, bElems[j].coef);
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	if (unpackedA) {
		packExponentVectors_AAZ_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	AA_SIZE(c) = k;
	if (k == 0) {
		freePolynomial_AAZ(c);
		return NULL;
	}
	resizePolynomial_AAZ_unpk(c, k);

	tryPackExponentVectors_AAZ_inp(c);

	return c;

}

AltArrZ_t* subPolynomials_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return addPolynomials_AAZ(a, b, nvar);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ_unpk(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// mpz_t ccoef;
	// mpz_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			mpz_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
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
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	if (unpackedA) {
		packExponentVectors_AAZ_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	AA_SIZE(c) = k;
	if (k == 0) {
		freePolynomial_AAZ(c);
		return NULL;
	}
	resizePolynomial_AAZ_unpk(c, k);

	tryPackExponentVectors_AAZ_inp(c);

	return c;

}

AltArrZ_t* addPolynomials_AAZ_inp_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return addPolynomials_AAZ_inp(a, b, nvar);
	}

	int unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ_unpk(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// mpz_t ccoef;
	// mpz_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			// cElems[k] = aElems[i];
			cElems[k].coef[0] = aElems[i].coef[0];
			// cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
			// cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			// mpz_init(cElems[k].coef);
			mpz_add(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpz_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			// cElems[k] = aElems[i];
			cElems[k].coef[0] = aElems[i].coef[0];
			// cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
			// cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			++k;
			++i;
		}
	}

	while(i < asize) {
		cElems[k].coef[0] = aElems[i].coef[0];
		// cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
		// cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, bElems[j].coef);
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(aDegs);
	free(a);

	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	AA_SIZE(c) = k;
	if (k == 0) {
		freePolynomial_AAZ(c);
		return NULL;
	}
	resizePolynomial_AAZ_unpk(c, k);

	tryPackExponentVectors_AAZ_inp(c);

	return c;
}

AltArrZ_t* subPolynomials_AAZ_inp_unpk(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AAZ_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return subPolynomials_AAZ_inp(a, b, nvar);
	}

	int unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ_unpk(asize + bsize, nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;
	AAZElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// mpz_t ccoef;
	// mpz_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			// cElems[k] = aElems[i];
			cElems[k].coef[0] = aElems[i].coef[0];
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			// mpz_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpz_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			// cElems[k] = aElems[i];
			cElems[k].coef[0] = aElems[i].coef[0];
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			++k;
			++i;
		}
	}

	while(i < asize) {
		cElems[k].coef[0] = aElems[i].coef[0];
		// cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
		// cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(aDegs);
	free(a);

	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	AA_SIZE(c) = k;
	if (k == 0) {
		freePolynomial_AAZ(c);
		return NULL;
	}
	resizePolynomial_AAZ_unpk(c, k);

	tryPackExponentVectors_AAZ_inp(c);

	return c;
}

void subPolynomials_AAZ_inpRHS_unpk(const AltArrZ_t* a, AltArrZ_t** bb) {
	if (bb == NULL) {
		return;
	}

	if (isZero_AAZ(a)) {
		negatePolynomial_AAZ(*bb);
		return;
	}

	if (isZero_AAZ(*bb)) {
		deepCopyPolynomial_AAZ_inp_unpk(a, bb);
	}

	if (a->nvar != (*bb)->nvar) {
		fprintf(stderr, "Polynomials are not compatible in subPolynomials_AAZ_inpRHS\n" );
		exit(1);
	}

	AltArrZ_t* b = *bb;

	if (!a->unpacked && !b->unpacked) {
		return subPolynomials_AAZ_inpRHS(a, bb);
	}

	int unpackedA = 0;
	if (!a->unpacked) {
		//yes, unsafe technically. But we convert back in the end.
		//TODO AB 2020-11-01 yeah.. this isn't thread safe so lets fix this
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
	}

	register int nvar = a->nvar;
	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AAZElem_t* cElems = (AAZElem_t*) malloc(sizeof(AAZElem_t) * (asize + bsize));
	degree_t* cDegs = (degree_t*) malloc(sizeof(degree_t)*(asize+bsize)*nvar);

	AAZElem_t* aElems = a->elems;
	AAZElem_t* bElems = b->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// mpz_t ccoef;
	// mpz_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			cElems[k] = bElems[j];
			mpz_neg(cElems[k].coef, cElems[k].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			// cElems[k] = aElems[i];
			cElems[k] = bElems[j];
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			mpz_sub(cElems[k].coef, aElems[i].coef, cElems[k].coef);
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpz_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, aElems[i].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			++k;
			++i;
		}
	}

	while(j < bsize) {
		cElems[k] = bElems[j];
		mpz_neg(cElems[k].coef, cElems[k].coef);
		// cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
		// cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		++k;
		++j;
	}
	while(i < asize) {
		mpz_init_set(cElems[k].coef, aElems[i].coef);
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++i;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(bElems);
	free(bDegs);
	b->elems = cElems;
	b->alloc = asize + bsize;
	b->size = k;

	if (unpackedA) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) a);
	}

	if (k == 0) {
		freePolynomial_AAZ(b);
		*bb = NULL;
		return;
	}

	resizePolynomial_AAZ_unpk(b, k);
	tryPackExponentVectors_AAZ_inp(b);
}


AltArrZ_t* CFDucosOptZ_subPolynomials_AAZ_inp_unpk (AltArrZ_t* a, AltArrZ_t* b, int nvar, AltArrZ_t** Pe, int e)
{
	if (a == NULL && b == NULL) {
		*Pe = NULL;
		return NULL;
	}

	if (!a->unpacked && !b->unpacked) {
		return CFDucosOptZ_subPolynomials_AAZ_inp(a, b, nvar, Pe, e);
	}

	int unpackedB = 0;
	if (a != NULL && !a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
	}
	if (b != NULL && !b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	register int asize = a == NULL ? 0 : AA_SIZE(a);
	register int bsize = b == NULL ? 0 : AA_SIZE(b);

	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a == NULL ? NULL : a->elems;
	AAZElem_t* bElems = b == NULL ? NULL : b->elems;
	AAZElem_t* cElems = c->elems;

	// int polyElemsAlloc = 10;
	// int polyElemsSize = 0;
	// AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);

	// degree_t* aDegs = a == NULL ? NULL : (degree_t*) a->elems->degs;
	// degree_t* bDegs = b == NULL ? NULL : (degree_t*) b->elems->degs;
	degree_t* cDegs = (degree_t*) c->elems->degs;

    // ratNum_t ccoef;
    // mpq_init(ccoef);

    // register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
	    //a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			setExponentVector_unpk(cElems[k].degs, bElems[j].degs, nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
	  	  // a==b
			// cElems[k] = aElems[i];
			cElems[k].coef[0] = aElems[i].coef[0];
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			setExponentVector_unpk(cElems[k].degs, aElems[i].degs, nvar);

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
		    cElems[k].coef[0] = aElems[i].coef[0];
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			setExponentVector_unpk(cElems[k].degs, aElems[i].degs, nvar);
			// cElems[k] = aElems[i];
		    // mpq_init(cElems[k].coef);
		    // mpq_set(cElems[k].coef, aElems[i].coef);
		    // cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		cElems[k].coef[0] = aElems[i].coef[0];
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		setExponentVector_unpk(cElems[k].degs, aElems[i].degs, nvar);
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

    //free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);
	tryPackExponentVectors_AAZ_inp(c);
	//TODO

	if (c == NULL || c->size == 0 || e < 0){
		*Pe = NULL;
	} else {
		*Pe = mainCoefficientAtIdx_AAZ(c, e);
	}

	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}


	return c;
}

ProductHeap_AAZ* prodheapInit_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b, int nvar) {
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	h->elements = (ProductHeapElem_AAZ*) malloc(sizeof(ProductHeapElem_AAZ)*AA_SIZE(a));
	degree_t* heapDegs = (degree_t*) malloc(sizeof(degree_t)*AA_SIZE(a)*nvar);
	h->elements[0].chain = prodheapMakeChain_AAZ(0, 0, NULL);
	int size = AA_SIZE(a);
	for (int i = 0; i < size; ++i ) {
		h->elements[i].degs = (degrees_t) (heapDegs + i*nvar);
	}
	h->unpackedDegs = heapDegs;
	addExponentVectors_unpk(a->elems->degs, b->elems->degs, h->elements->degs, nvar);
	h->heapSize = 1;
	h->maxHeapSize = AA_SIZE(a);
	h->nvar = nvar;
	return h;
}

static inline void prodheapFree_AAZ_unpk(ProductHeap_AAZ* h) {
	free(h->unpackedDegs);
	prodheapFree_AAZ(h);
}

static inline void prodheapSetElem_unpk(ProductHeapElem_AAZ* dest, ProductHeapElem_AAZ* src, int nvar) {
	dest->chain = src->chain;
	setExponentVector_unpk(dest->degs, src->degs, nvar);
}

void prodheapInsert_AAZ_unpk(ProductHeap_AAZ* h, ProductHeapChain_AAZ* chain, degrees_t degs) {
	register int s = h->heapSize;
	ProductHeapElem_AAZ* elems = h->elements;
	int nvar = h->nvar;

	if (s == 0) {
		setExponentVector_unpk(elems[0].degs, degs, nvar);
		elems[0].chain = chain;
		h->heapSize = 1;
		return;
	}

	//first check if we can chain off the root
	if (isEqualExponentVectors_unpk(elems[0].degs, degs, nvar)) {
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
		if (isEqualExponentVectors_unpk(elems[i].degs, degs, nvar)) {
			chain->next = elems[i].chain;
			elems[i].chain = chain;
			return;
		} else if (isLessExponentVectors_unpk(elems[i].degs, degs, nvar)) {
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
	degree_t tempDegs[nvar];
	ProductHeapElem_AAZ temp = {(degrees_t) tempDegs, NULL};
	setExponentVector_unpk(elems[s].degs, degs, nvar);
	ProductHeapElem_AAZ elem = {elems[s].degs, chain};

	//TODO use i index again to swap between elements rather than use a second temp elem.
	while (j <= s) {
		prodheapSetElem_unpk(&temp, elems + j, nvar);
		// temp = elems[j];
		prodheapSetElem_unpk(elems + j, &elem, nvar);
		// elems[j] = elem;
		prodheapSetElem_unpk(&elem, &temp, nvar);
		// elem = temp;
		j = (j << 1) + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

ProductHeapChain_AAZ* prodheapRemoveMax_AAZ_unpk(ProductHeap_AAZ* h) {
	ProductHeapElem_AAZ* elems = h->elements;
	ProductHeapChain_AAZ* maxElem = elems[0].chain;
	register int i = 0;
	register int j = 1;
	register int s = --(h->heapSize);
	int nvar = h->nvar;

	//promote largest children
	while (j < s) {
		if (j+1 < s && isLessExponentVectors_unpk(elems[j].degs, elems[j+1].degs, nvar)) {
			++j;
		}
		prodheapSetElem_unpk(elems + i, elems + j, nvar);
		// elems[i] = elems[j];
		i = j;
		j = (j << 1) + 1;
	}
	//now place last element into i and swim up to make tree complete
	j = (i-1) >> 1;
	while(i > 0) {
		if (isLessExponentVectors_unpk(elems[s].degs, elems[j].degs, nvar)) {
			break;
		}
		prodheapSetElem_unpk(elems + i, elems + j, nvar);
		// elems[i] = elems[j];
		i = j;
		j = (j-1) >> 1;
	}
	prodheapSetElem_unpk(elems + i, elems + s, nvar);
	// elems[i] = elems[s];

	return maxElem;

}

AltArrZ_t* multiplyPolynomials_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b, int nvar) {
	if (a == NULL || a->size == 0 || b == NULL || b->size == 0) {
		return NULL;
	}

	if (!a->unpacked && !b->unpacked) {
		return multiplyPolynomials_AAZ(a,b,nvar);
	}

	// reorder to obtain smaller as a.
	if (b->size < a->size) {
		const AltArrZ_t* temp = a;
		a = b;
		b = temp;
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) b);
		unpackedB = 1;
	}

	ProductHeap_AAZ* h = prodheapInit_AAZ_unpk(a,b,nvar);

	register unsigned long long int allocC = AA_SIZE(a)*AA_SIZE(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;
	AltArrZ_t* c = makePolynomial_AAZ_unpk(allocC, nvar);

	//k is c, i is a, j is b.
	register int k = 0;
	// register int i = 0;
	// register int j = 0;

	AAZElem_t* __restrict__ cElems = c->elems;
	AAZElem_t* __restrict__ bElems = b->elems;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;

	AAZElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1;
	register int firstB = 0;

	ProductHeapChain_AAZ* maxElem = NULL;
	ProductHeapChain_AAZ* nextMaxElem = NULL;
	degrees_t* nextDegs;
	degree_t tempDegs[nvar];
	for (int i = 0; i < nvar; ++i) {
		tempDegs[i] = 0;
	}
	degree_t* tempDegs_p = tempDegs;

	while ( (nextDegs = prodheapPeek_AAZ(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		setExponentVector_unpk((degrees_t) tempDegs_p, h->elements[0].degs, nvar);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		setExponentVector_unpk(cElems[k].degs, *nextDegs, nvar);
		// cElems[k].degs = *nextDegs;
		mpz_init(cElems[k].coef);

		while (nextDegs != NULL && isEqualExponentVectors_unpk(cElems[k].degs, *nextDegs, nvar)) {
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AAZ* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AAZ_unpk(h);
			while (maxElem != NULL) {
				mpz_addmul(cElems[k].coef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA) {
					oldMaxElem = prodheapMakeChain_AAZ((maxElem->a_i)+1, firstB, oldMaxElem);
				}

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
			//will init cElem[k] again on next loop, so clear here.
			mpz_clear(cElems[k].coef);
		}

		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			addExponentVectors_unpk(aElems[maxElem->a_i].degs, bElems[maxElem->b].degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AAZ_unpk(h, maxElem, (degrees_t) tempDegs_p);
			maxElem = nextMaxElem;
		}

		if (k >= allocC) {
			allocC <<= 1;
			resizePolynomial_AAZ_unpk(c, allocC);
			// fprintf(stderr, "\n\nRESIZING\n\n");
			cElems = c->elems;
		}
	}

	prodheapFree_AAZ_unpk(h);

	AA_SIZE(c) = k;

	resizePolynomial_AAZ_unpk(c, k);

	if (unpackedA) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) a);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) b);
	}

	return c;
}

void multiplyPolynomialsPreAlloc_AAZ_unpk(const AltArrZ_t* a, const AltArrZ_t* b, AltArrZ_t** cc) {
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
		exit(1);
	}

	if (!a->unpacked && !b->unpacked) {
		return multiplyPolynomialsPreAlloc_AAZ(a,b,cc);
	}

	// reorder to obtain smaller as a.
	if (b->size < a->size) {
		const AltArrZ_t* temp = a;
		a = b;
		b = temp;
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) b);
		unpackedB = 1;
	}

	int nvar = a->nvar;
	ProductHeap_AAZ* h = prodheapInit_AAZ_unpk(a,b,nvar);

	register unsigned long long int allocC = AA_SIZE(a)*AA_SIZE(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;
	AltArrZ_t* c = *cc;
	if (c == NULL) {
		c = makePolynomial_AAZ_unpk(allocC, nvar);
		*cc = c;
	} else {
		if (!c->unpacked) {
			c->nvar = nvar;
			unpackExponentVectors_AAZ_inp(c);
		}
		allocC = c->alloc;
	}


	//k is c, i is a, j is b.
	register polysize_t k = 0;
	register polysize_t cSize = c->size;
	// register int i = 0;
	// register int j = 0;

	AAZElem_t* __restrict__ cElems = c->elems;
	AAZElem_t* __restrict__ bElems = b->elems;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;


	AAZElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1;
	register int firstB = 0;

	ProductHeapChain_AAZ* maxElem = NULL;
	ProductHeapChain_AAZ* nextMaxElem = NULL;
	degrees_t* nextDegs;
	degree_t tempDegs[nvar];
	degree_t* tempDegs_p = tempDegs;

	while ( (nextDegs = prodheapPeek_AAZ(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		setExponentVector_unpk((degrees_t) tempDegs_p, h->elements[0].degs, nvar);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		setExponentVector_unpk(cElems[k].degs, *nextDegs, nvar);
		// cElems[k].degs = *nextDegs;
		if (k >= cSize) {
			mpz_init(cElems[k].coef);
		} else {
			mpz_set_ui(cElems[k].coef, 0ul);
		}

		while (nextDegs != NULL && isEqualExponentVectors_unpk(cElems[k].degs, *nextDegs, nvar)) {
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AAZ* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AAZ_unpk(h);
			while (maxElem != NULL) {
				mpz_addmul(cElems[k].coef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA) {
					oldMaxElem = prodheapMakeChain_AAZ((maxElem->a_i)+1, firstB, oldMaxElem);
				}

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
			//will init cElem[k] again on next loop, so clear here.
			if (k >= cSize) {
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
			addExponentVectors_unpk(aElems[maxElem->a_i].degs, bElems[maxElem->b].degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AAZ_unpk(h, maxElem, (degrees_t) tempDegs_p);
			maxElem = nextMaxElem;
		}

		if (k >= allocC) {
			allocC <<= 1;
			resizePolynomial_AAZ_unpk(c, allocC);
			// fprintf(stderr, "\n\nRESIZING\n\n");
			cElems = c->elems;
			cDegs = (degree_t*) c->elems[0].degs;
		}
	}

	prodheapFree_AAZ_unpk(h);

	AA_SIZE(c) = k;
	for ( ; k < cSize; ++k) {
		mpz_clear(cElems[k].coef);
	}
	resizePolynomial_AAZ_unpk(c, k);

	if (unpackedA) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) a);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) b);
	}

}


/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
void divisionGetNextTerm_AAZ_unpk(ProductHeap_AAZ* h, const AAZElem_t* __restrict__ aElems, const AAZElem_t* __restrict__ bElems, mpz_t* retCoef) {
	if (h->heapSize == 0) {
		return;
	}

	int lastB = h->lastB;
	int nvar = h->nvar;

	ProductHeapChain_AAZ* insertChain = NULL;
	ProductHeapChain_AAZ* maxElem, *nextMaxElem;

	mpz_t prodCoef;
	mpz_init(prodCoef);
	degrees_t* nextDegs = prodheapPeek_AAZ(h);

	degree_t maxDegsUnpk[nvar];
	degrees_t maxDegs = (degrees_t) maxDegsUnpk;
	setExponentVector_unpk(maxDegs, *nextDegs, nvar);
	// degrees_t maxDegs = *nextDegs;

	while ( nextDegs != NULL && isEqualExponentVectors_unpk(maxDegs, *nextDegs, nvar)) {
		maxElem = prodheapRemoveMax_AAZ_unpk(h);

		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;
			mpz_mul(prodCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
			mpz_add(*retCoef, *retCoef, prodCoef);
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

	degree_t tempDegs[nvar];
	degree_t* tempDegs_p = tempDegs;
	while(insertChain != NULL) {
		maxElem = insertChain->next;
		insertChain->next = NULL;
		addExponentVectors_unpk(aElems[insertChain->a_i].degs, bElems[insertChain->b].degs, (degrees_t) tempDegs_p, nvar);
		prodheapInsert_AAZ_unpk(h, insertChain, (degrees_t) tempDegs_p);
		insertChain = maxElem;
	}

	mpz_clear(prodCoef);
}

void divideBySingleTerm_AAZ_unpk(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar) {
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

	if (!c->unpacked && !b->unpacked) {
		divideBySingleTerm_AAZ(c, b, res_a, res_r, nvar);
		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*)c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*)b);
		unpackedB = 1;
	}

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxSize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ_unpk(maxSize, nvar);

	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	AAZElem_t* curR = r->elems;
	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* rDegs = (degree_t*) r->elems[0].degs;

	mpz_init(curA->coef);
	mpz_init(curR->coef);

	while (k != lenK) {
		if (monomialDivideTest_unpk(k->degs,bElems->degs,nvar)) {
			mpz_div(curA->coef, k->coef, bElems->coef);
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(k->degs, bElems->degs, curA->degs, nvar);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			mpz_set(curR->coef, k->coef);
			curR->degs = (degrees_t) (rDegs + j*nvar);
			setExponentVector_unpk(curR->degs, k->degs, nvar);
			// curR->degs = k->degs;
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

	if (unpackedC) {
		packExponentVectors_AAZ_inp((AltArrZ_t*)c);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*)b);
	}

	tryPackExponentVectors_AAZ_inp(a);
	tryPackExponentVectors_AAZ_inp(r);

	*res_a = a;
	*res_r = r;
}

 void prodheapResize_AAZ_unpk(ProductHeap_AAZ* h, int newAllocSize) {
	h->elements = (ProductHeapElem_AAZ*) realloc(h->elements, sizeof(ProductHeapElem_AAZ)*newAllocSize);
	h->maxHeapSize  = newAllocSize;

	h->unpackedDegs = (degree_t*) realloc(h->unpackedDegs, sizeof(degree_t)*h->nvar*newAllocSize);
	degree_t* heapDegs = h->unpackedDegs;
	int nvar = h->nvar;
	for (int i = 0; i < newAllocSize; ++i ) {
		h->elements[i].degs = (degrees_t) (heapDegs + i*nvar);
	}
}

/**
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials_AAZ_unpk(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, register int nvar) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		if (res_a != NULL) {
			*res_a = makePolynomial_AAZ(0, nvar);
		}
		if (res_r != NULL) {
			*res_r = NULL;
		}
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		divideBySingleTerm_AAZ(c, b, res_a, res_r, nvar);
		return;
	}

	if (!c->unpacked && !b->unpacked) {
		dividePolynomials_AAZ(c, b, res_a, res_r, nvar);
		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*)c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*)b);
		unpackedB = 1;
	}


	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ_unpk(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	AAZElem_t* __restrict__ curR = r->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems[0].degs;
	degree_t* __restrict__ rDegs = (degree_t*) r->elems[0].degs;
	degree_t tempDegs[nvar];
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpz_init(curA->coef);
	mpz_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	while (k != lenK && !monomialDivideTest_unpk(k->degs, beta, nvar)) {
		mpz_set(curR->coef, k->coef);
		curR->degs = (degrees_t) (rDegs + j*nvar);
		setExponentVector_unpk(curR->degs, k->degs, nvar);
		// curR->degs = k->degs;
		++j;
		if (j >= maxRSize) {
			maxRSize <<= 1;
			r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
			degree_t* oldDegs = rDegs;
			rDegs = realloc(rDegs, sizeof(degree_t)*nvar*maxRSize);
			if (oldDegs != rDegs) {
				for (int idx = 0; idx < j; ++idx) {
					r->elems[idx].degs = (degrees_t) (rDegs + idx*nvar);
				}
			}
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

		if (res_a != NULL) {
			*res_a = a;
		} else {
			freePolynomial_AAZ(a);
		}
		if (res_r != NULL) {
			*res_r = r;
		} else {
			freePolynomial_AAZ(r);
		}
		return;
	}

	curA->degs = (degrees_t) (aDegs + i*nvar);
	subtractExponentVectors_unpk(k->degs, beta, curA->degs, nvar);
	mpz_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps = (degrees_t) tempDegs_p;
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
			cmp = compareExponentVectors_unpk(*delta, k->degs, nvar);
		}

		if (cmp > 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
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
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest_unpk(eps, beta, nvar)) {
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			mpz_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {

			//swap here so that curA becomes 0.
			mpz_swap(curR->coef, curA->coef);
			curR->degs = (degrees_t) (rDegs + j*nvar);
			setExponentVector_unpk(curR->degs, eps, nvar);
			// curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
				degree_t* oldDegs = rDegs;
				rDegs = realloc(rDegs, sizeof(degree_t)*nvar*maxRSize);
				if (oldDegs != rDegs) {
					for (int idx = 0; idx < j; ++idx) {
						r->elems[idx].degs = (degrees_t) (rDegs + idx*nvar);
					}
				}
				curR = r->elems + j - 1;
			}
			++(curR);
			mpz_init(curR->coef);
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ_unpk(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	if (unpackedC) {
		packExponentVectors_AAZ_inp((AltArrZ_t*)c);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*)b);
	}

	if (res_a != NULL) {
		tryPackExponentVectors_AAZ_inp(a);
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
	if (res_r != NULL) {
		tryPackExponentVectors_AAZ_inp(r);
		*res_r = r;
	} else {
		freePolynomial_AAZ(r);
	}
}

void exactDividePolynomials_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, register int nvar) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		if (res_a != NULL) {
			*res_a = NULL;
		}
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		AltArrZ_t* res_r = NULL;
		divideBySingleTerm_AAZ(c, b, res_a, &res_r, nvar);
		freePolynomial_AAZ(res_r);
		return;
	}

	if (!c->unpacked && !b->unpacked) {
		exactDividePolynomials_AAZ(c, b, res_a, nvar);
		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	// register int maxRSize = maxASize;
	register int i = 0;
	// register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxASize, nvar);
	// AltArrZ_t* r = makePolynomial_AAZ_unpk(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	// AAZElem_t* __restrict__ curR = r->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems[0].degs;
	// degree_t* __restrict__ rDegs = (degree_t*) r->elems[0].degs;
	degree_t tempDegs[nvar];
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpz_init(curA->coef);
	// mpz_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	// while (k != lenK && !monomialDivideTest_unpk(k->degs, beta, nvar)) {
	// 	mpz_set(curR->coef, k->coef);
	// 	curR->degs = (degrees_t) (rDegs + j*nvar);
	// 	setExponentVector_unpk(curR->degs, k->degs, nvar);
	// 	// curR->degs = k->degs;
	// 	++j;
	// 	if (j >= maxRSize) {
	// 		maxRSize <<= 1;
	// 		r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
	// 		degree_t* oldDegs = rDegs;
	// 		rDegs = realloc(rDegs, sizeof(degree_t)*nvar*maxRSize);
	// 		if (oldDegs != rDegs) {
	// 			for (int idx = 0; idx < j; ++idx) {
	// 				r->elems[idx].degs = (degrees_t) (rDegs + idx*nvar);
	// 			}
	// 		}
	// 		curR = r->elems + j - 1;
	// 	}
	// 	++(curR);
	// 	mpz_init(curR->coef);
	// 	++k;
	// }
	// if (k == lenK) {
	// 	//no division to do at all!
	// 	mpz_clear(curA->coef);
	// 	mpz_clear(curR->coef);

	// 	AA_SIZE(a) = i;
	// 	AA_SIZE(r) = j;
	// 	a->alloc = maxASize;
	// 	r->alloc = maxRSize;

	// 	*res_a = a;
	// 	*res_r = r;
	// 	return;
	// }

	curA->degs = (degrees_t) (aDegs + i*nvar);
	subtractExponentVectors_unpk(k->degs, beta, curA->degs, nvar);
	mpz_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps = (degrees_t) tempDegs_p;
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
			cmp = compareExponentVectors_unpk(*delta, k->degs, nvar);
		}

		if (cmp > 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
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
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest_unpk(eps, beta, nvar)) {
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			mpz_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			mpz_init(curA->coef);
		}
		// else {

		// 	//swap here so that curA becomes 0.
		// 	mpz_swap(curR->coef, curA->coef);
		// 	curR->degs = (degrees_t) (rDegs + j*nvar);
		// 	setExponentVector_unpk(curR->degs, eps, nvar);
		// 	// curR->degs = eps;
		// 	++j;
		// 	if (j >= maxRSize) {
		// 		maxRSize <<= 1;
		// 		r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
		// 		degree_t* oldDegs = rDegs;
		// 		rDegs = realloc(rDegs, sizeof(degree_t)*nvar*maxRSize);
		// 		if (oldDegs != rDegs) {
		// 			for (int idx = 0; idx < j; ++idx) {
		// 				r->elems[idx].degs = (degrees_t) (rDegs + idx*nvar);
		// 			}
		// 		}
		// 		curR = r->elems + j - 1;
		// 	}
		// 	++(curR);
		// 	mpz_init(curR->coef);
		// }

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ_unpk(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	// mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	// AA_SIZE(r) = j;
	a->alloc = maxASize;
	// r->alloc = maxRSize;

	if (unpackedC) {
		packExponentVectors_AAZ_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	if (res_a != NULL) {
		tryPackExponentVectors_AAZ_inp(a);
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
}

//TODO
void univariatePseudoDivideBySingleTerm_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy) {
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

	if (!c->unpacked && !b->unpacked) {
		univariatePseudoDivideBySingleTerm_AAZ(c, b, res_a, res_r, e, lazy);
		return;
	}

	int nvar = 1;
	c = deepCopyPolynomial_AAZ(c);

	int unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp(c);
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	if (isLessExponentVectors_unpk(c->elems->degs, b->elems->degs, nvar)) {
		if (res_r != NULL) {
			*res_r = c; //already did deep copy
			tryPackExponentVectors_AAZ_inp(*res_r);
		} else {
			freePolynomial_AAZ(c);
		}
		if (res_a != NULL) {
			*res_a = NULL;
		}
		if (e != NULL) {
			*e = 0;
		}
		if (unpackedB) {
			packExponentVectors_AAZ_inp(b);
		}
		return;
	}

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxSize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ_unpk(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	AAZElem_t* curR = r->elems;
	degree_t* aDegs = (degree_t*) a->elems->degs;
	degree_t* rDegs = (degree_t*) r->elems->degs;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	int multSteps = 0;

	while (k != lenK) {
		if (monomialDivideTest_unpk(k->degs,bElems->degs,nvar)) {
			if (mpz_divisible_p(k->coef, bElems->coef) == 0) {
				multiplyByInteger_AAZ_inp(c, bElems->coef);
				multiplyByInteger_AAZ_inp(a, bElems->coef);
				++multSteps;
			}

			mpz_divexact(curA->coef, k->coef, bElems->coef);
			// mpz_set(curA->coef, k->coef);
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(k->degs, bElems->degs, curA->degs, nvar);
			++i;
			++(curA);
			++(a->size);
			mpz_init(curA->coef);
		} else {
			mpz_set(curR->coef, k->coef);
			curR->degs = (degrees_t) (rDegs + j*nvar);
			setExponentVector_unpk(curR->degs, k->degs, nvar);
			// curR->degs = k->degs;
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
		degrees_t d = ((degree_t*) c->elems->degs)[0] - ((degree_t*) b->elems->degs)[0] + 1 - multSteps;
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

	if (res_a != NULL) {
		tryPackExponentVectors_AAZ_inp(a);
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
	if (res_r != NULL) {
		tryPackExponentVectors_AAZ_inp(r);
		*res_r = r;
	} else {
		freePolynomial_AAZ(r);
	}

	if (e != NULL) {
		*e = multSteps;
	}

	freePolynomial_AAZ(c);

	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}


}

void univariatePseudoDividePolynomials_AAZ_unpk(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	int nvar = 1;
	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		if (res_a != NULL) {
			*res_a = makePolynomial_AAZ(0, nvar);
		}
		if (res_r != NULL) {
			*res_r = NULL;
		}
		return;
	}

	if (!c->unpacked && !b->unpacked) {
		univariatePseudoDividePolynomials_AAZ(c, b, res_a, res_r, e, lazy);
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		univariatePseudoDivideBySingleTerm_AAZ(c, b, res_a, res_r, e, lazy);
		return;
	}

	c = deepCopyPolynomial_AAZ(c);

	int unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp(c);
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	if (isLessExponentVectors_unpk(c->elems->degs, b->elems->degs, nvar)) {
		if (res_r != NULL) {
			*res_r = c;
			tryPackExponentVectors_AAZ_inp(*res_r);
		} else {
			freePolynomial_AAZ(c);
		}
		if (res_a != NULL) {
			*res_a = NULL;
		}
		if (e != NULL) {
			*e = 0;
		}
		if (unpackedB) {
			packExponentVectors_AAZ_inp(b);
		}
		return;
	}


	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ_unpk(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	AAZElem_t* __restrict__ curR = r->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems->degs;
	degree_t* __restrict__ rDegs = (degree_t*) r->elems->degs;
	degree_t tempDegs[nvar];
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpz_init(curA->coef);
	mpz_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	int multSteps = 0;
	subtractExponentVectors_unpk(k->degs, beta, curA->degs, nvar);
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
	prodheapResize_AAZ_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	a->size = 1;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps = (degrees_t) tempDegs_p;
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
			cmp = compareExponentVectors_unpk(*delta, k->degs, nvar);
		}

		if (cmp > 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
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
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (compareExponentVectors_unpk(eps, beta, nvar) >= 0) {
			if (mpz_divisible_p(curA->coef, b->elems->coef) == 0) {
				multiplyByInteger_AAZ_inp(c, b->elems->coef);
				multiplyByInteger_AAZ_inp(a, b->elems->coef);
				++multSteps;

			} else {
				mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			}

			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = (degree_t*) realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			++(a->size);
			mpz_init(curA->coef);
		} else {
			//swap here so that curA becomes 0.
			mpz_swap(curR->coef, curA->coef);
			curR->degs = (degrees_t) (rDegs + j*nvar);
			setExponentVector_unpk(curR->degs, eps, nvar);
			// curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
				degree_t* oldDegs = rDegs;
				rDegs = realloc(rDegs, sizeof(degree_t)*nvar*maxRSize);
				if (oldDegs != rDegs) {
					for (int idx = 0; idx < j; ++idx) {
						r->elems[idx].degs = (degrees_t) (rDegs + idx*nvar);
					}
				}
				curR = r->elems + j - 1;
			}
			++(curR);
			mpz_init(curR->coef);
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ_unpk(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	if (!lazy) {
		int cDegs = *((degree_t*) c->elems->degs);
		int bDegs = *((degree_t*) b->elems->degs);
		degrees_t d = cDegs - bDegs + 1 - multSteps;
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

	if (res_a != NULL) {
		tryPackExponentVectors_AAZ_inp(a);
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
	if (res_r != NULL) {
		tryPackExponentVectors_AAZ_inp(r);
		*res_r = r;
	} else {
		freePolynomial_AAZ(r);
	}

	if (e != NULL) {
		*e = multSteps;
	}

	freePolynomial_AAZ(c);

	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

}

int divideTestSingleTerm_AAZ_unpk(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {
	if (b == NULL || b->size == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		return 1;
	}

	if (!c->unpacked && !b->unpacked) {
		return divideTestSingleTerm_AAZ(c, b, res_a, nvar);

	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) b);
		unpackedB = 1;
	}

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxSize, nvar);

	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	degree_t* aDegs = (degree_t*) a->elems[0].degs;

	mpz_init(curA->coef);
	int ret = 1;
	while (k != lenK) {
		if (monomialDivideTest_unpk(k->degs,bElems->degs,nvar) && mpz_divisible_p(k->coef, bElems->coef)) {
			mpz_div(curA->coef, k->coef, bElems->coef);
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(k->degs, bElems->degs, curA->degs, nvar);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			mpz_clear(curA->coef);
			a->size = i;
			freePolynomial_AAZ(a);
			a = NULL;
			ret = 0;
			break;
		}
		++k;
	}

	if (a != NULL) {
		mpz_clear(curA->coef);
	}

	if (ret) {
		AA_SIZE(a) = i;
		if (res_a != NULL) {
			tryPackExponentVectors_AAZ_inp(a);
			*res_a = a;
		} else {
			freePolynomial_AAZ(a);
		}
	}

	if (unpackedC) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) c);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) b);
	}

	return ret;
}


int divideTest_AAZ_unpk(const AltArrZ_t* c, const AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {

	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AAZ(0, nvar);
		return 1;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		return divideTestSingleTerm_AAZ_unpk(c, b, res_a, nvar);
	}

	if (!c->unpacked && !b->unpacked) {
		return divideTest_AAZ(c, b, res_a, nvar);
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		//yes, unsafe cast. It'll be put back after.
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) b);
		unpackedB = 1;
	}


	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int i = 0;

	AltArrZ_t* a = makePolynomial_AAZ_unpk(maxASize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems[0].degs;
	degree_t tempDegs[nvar];
	for (int i = 0; i < nvar; ++i) {
		tempDegs[i] = 0;
	}
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpz_init(curA->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	if (!monomialDivideTest_unpk(k->degs, beta, nvar) || !mpz_divisible_p(k->coef, b->elems->coef)) {
		freePolynomial_AAZ(a);
		if (unpackedC) {
			packExponentVectors_AAZ_inp((AltArrZ_t*) c);
		}
		if (unpackedB) {
			packExponentVectors_AAZ_inp((AltArrZ_t*) b);
		}
		return 0;
	}

	curA->degs = (degrees_t) (aDegs + i*nvar);
	subtractExponentVectors_unpk(k->degs, beta, curA->degs, nvar);
	mpz_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps = (degrees_t) tempDegs_p;
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
			cmp = compareExponentVectors_unpk(*delta, k->degs, nvar);
		}

		if (cmp > 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AAZ_unpk(h, a->elems, b->elems, &(curA->coef));
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
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest_unpk(eps, beta, nvar)  && mpz_divisible_p(curA->coef, b->elems->coef)) {
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAZ_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AAZ_unpk(h, prodheapMakeChain_AAZ(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			prodheapFree_AAZ_unpk(h);
			mpz_clear(curA->coef);

			a->size = i;
			freePolynomial_AAZ(a);
			if (unpackedC) {
				packExponentVectors_AAZ_inp((AltArrZ_t*) c);
			}
			if (unpackedB) {
				packExponentVectors_AAZ_inp((AltArrZ_t*) b);
			}
			return 0;
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ_unpk(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxASize;

	if (unpackedC) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) c);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp((AltArrZ_t*) b);
	}

	if (res_a != NULL) {
		tryPackExponentVectors_AAZ_inp(a);
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
	return 1;
}


void divideByLeadingTerms_AAZ_unpk (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar)
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

	if (!c->unpacked && !b->unpacked) {
		divideByLeadingTerms_AAZ (c, b, res_a, res_r, nvar);
 		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AAZ_inp (c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp (b);
		unpackedB = 1;
	}

	AltArrZ_t* a;
	AltArrZ_t* r;

	/* while (k != lenK) { */
	if (monomialDivideTest_AAZ_unpk (c, 0, b, 0)) {
		// a:
		a = makePolynomial_AAZ_unpk (1, nvar);
		mpz_init (a->elems->coef);
		mpz_div (a->elems->coef, c->elems->coef, b->elems->coef);
		subtractExponentVectors_unpk (c->elems->degs, b->elems->degs, a->elems->degs, nvar);
		AA_SIZE (a) = 1;

		// r:
		r = NULL;

	} else {
		// a:
		a = makePolynomial_AAZ_unpk (1, nvar);
		mpz_init (a->elems->coef);
		mpz_set_si (a->elems->coef, 1l);
		a->elems->degs = (degrees_t) 0;
		AA_SIZE (a) = 1;

		// r:
		r = makePolynomial_AAZ_unpk (1, nvar);
		mpz_init (r->elems->coef);
		mpz_set (r->elems->coef, c->elems->coef);
		/* r->elems->degs = c->elems->degs; */
		setExponentVector_unpk (r->elems->degs, c->elems->degs, nvar);
		AA_SIZE (r) = 1;
	}

	if (unpackedC) {
		packExponentVectors_AAZ_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	if (res_a != NULL) {
		tryPackExponentVectors_AAZ_inp(a);
		*res_a = a;
	} else {
		freePolynomial_AAZ(a);
	}
	if (res_r != NULL) {
		tryPackExponentVectors_AAZ_inp(r);
		*res_r = r;
	} else {
		freePolynomial_AAZ(r);
	}
}



/*****************
 * Derivative / Integral
 *****************/

AltArrZ_t* derivative_AAZ_unpk(const AltArrZ_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx < 0 || idx >= aa->nvar) {
		return NULL;
	}

	if (!aa->unpacked) {
		return derivative_AAZ(aa, idx, k);
	}

	int nvar = aa->nvar;
	AltArrZ_t* ret = makePolynomial_AAZ_unpk(aa->size, aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAZElem_t* __restrict__ elems = aa->elems;
	degree_t* degs_p = (degree_t*) elems[0].degs;
	AAZElem_t* __restrict__ retElems = ret->elems;
	degree_t* retDegs_p = (degree_t*) retElems[0].degs;
	degrees_t deg;
	mpz_t mpzDeg;
	mpz_init(mpzDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i) {
		deg = degs_p[i*nvar + idx];
		// deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);
		if (deg < k) {
			continue;
		}

		retElems[insertIdx].degs = (degrees_t) (retDegs_p + insertIdx*nvar);
		setExponentVector_unpk(retElems[insertIdx].degs, elems[i].degs, nvar);
		retDegs_p[insertIdx*nvar + idx] = 0;
		// retElems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpz_init(retElems[insertIdx].coef);
		mpz_set(retElems[insertIdx].coef, elems[i].coef);

		mpz_set_ui(mpzDeg, deg);
		for (int j = 0; j < k; ++j) {
			mpz_mul(retElems[insertIdx].coef, retElems[insertIdx].coef, mpzDeg);
			--deg;
			mpz_sub(mpzDeg, mpzDeg, mpzOne);
		}
		retDegs_p[insertIdx*nvar + idx] = deg;
		// retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpz_clear(mpzDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	if (insertIdx == 0) {
		freePolynomial_AAZ(ret);
		return NULL;
	}

	tryPackExponentVectors_AAZ_inp(ret);

	return ret;
}

void derivative_AAZ_inp_unpk(AltArrZ_t** aa_p, int idx, int k) {
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

	if (!(*aa_p)->unpacked) {
		derivative_AAZ_inp(aa_p, idx, k);
		return;
	}

	AltArrZ_t* aa = *aa_p;
	int nvar = aa->nvar;

	int insertIdx = 0;
	int size = aa->size;
	AAZElem_t* elems = aa->elems;
	degree_t* degs_p = (degree_t*) elems[0].degs;
	degrees_t deg;
	for (int i = 0; i < size; ++i) {
		deg = degs_p[i*nvar + idx];
		// deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);
		if (deg < k) {
			continue;
		}

		setExponentVector_unpk(elems[insertIdx].degs, elems[i].degs, nvar);
		degs_p[insertIdx*nvar + idx] = 0;
		// retElems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpz_set(elems[insertIdx].coef, elems[i].coef);

		for (int j = 0; j < k; ++j) {
			mpz_mul_ui(elems[insertIdx].coef, elems[insertIdx].coef, deg);
			--deg;
		}
		degs_p[insertIdx*nvar + idx] = deg;
		// retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	for (int i = insertIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
	}

	aa->size = insertIdx;
	if (insertIdx == 0) {
		freePolynomial_AAZ(aa);
		*aa_p = NULL;
		return;
	}
	tryPackExponentVectors_AAZ_inp(aa);
}


AltArrZ_t* primitivePartAndContent_AAZFromAA_unpk(AltArr_t* aa, mpq_t cont) {
	if (aa == NULL || aa->size == 0) {
		mpq_set_si(cont, 0l, 1l);
		return NULL;
	}

	if (!aa->unpacked) {
		return primitivePartAndContent_AAZFromAA(aa, cont);
	}

	integralContent_AA(aa, cont);

	AltArrZ_t* res = makePolynomial_AAZ_unpk(aa->size, aa->nvar);
	int size = res->size = aa->size;
	AAElem_t* elems = aa->elems;
	AAZElem_t* reselems = res->elems;
	degree_t* resdegs = (degree_t*) res->elems->degs;

	int nvar = aa->nvar;
	mpq_t tmp;
	mpq_init(tmp);
	for (int i = 0; i < size; ++i) {
		mpq_div(tmp, elems[i].coef, cont);
		mpz_init_set(reselems[i].coef, mpq_numref(tmp));

		reselems[i].degs = (degrees_t) (resdegs + i*nvar);
		setExponentVector_unpk(reselems[i].degs, elems[i].degs, nvar);
	}

	mpq_clear(tmp);

	return res;

}


AltArrZ_t* contentInVars_AAZ_unpk(const AltArrZ_t* aa_in, const int* active) {
	if (isZero_AAZ(aa_in)) {
		return NULL;
	}

	if (!aa_in->unpacked) {
		return contentInVars_AAZ(aa_in, active);
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

	AltArrZ_t localCoef;
	localCoef.unpacked = aa->unpacked;
	localCoef.nvar = aa->nvar;

	AltArrZ_t* tmpCoef = NULL;
	AAZElem_t* elems = aa->elems;
	degree_t* degs = (degree_t*) elems->degs;
	degree_t* tmpDegs = NULL;

	int j;
	degree_t curDegs[nvar];
	degree_t nextDegs[nvar];
	for (j = 0; j < nactive; ++j) {
		curDegs[j] = degs[j];
	}
	for (j = nactive; j < nvar; ++j) {
		curDegs[j] = 0;
		nextDegs[j] = 0;
	}

	int startIdx = 0;
	int i, k, same;
	for (i = 0; i < aa->size; ++i) {
		same = 1;
		for (j = 0; j < nactive; ++j) {
			nextDegs[j] = degs[i*nvar + j];
			if (nextDegs[j] != curDegs[j]) {
				same = 0;
			}
		}

		if (!same) {
			//we have a new coef from prevCoef,...,i-1
			localCoef.size = i - startIdx;
			localCoef.alloc = localCoef.size;
			localCoef.elems = aa->elems + startIdx;
			deepCopyPolynomial_AAZ_inp_unpk(&localCoef, &tmpCoef);
			tmpDegs = (degree_t*) tmpCoef->elems->degs;

			//simultaneously copy degrees into tmpCoef meanwhile zero-ing the active variables
			for (k = 0; k < tmpCoef->size; ++k) {
				for (j = 0; j < nactive; ++j) {
					tmpDegs[k*nvar + j] = 0;
				}
				for (j = nactive; j < nvar; ++j) {
					tmpDegs[k*nvar + j] = degs[(startIdx+k)*nvar + j];
				}
			}

			if (cont == NULL) {
				cont = tmpCoef;
				tmpCoef = NULL;
			} else {
				AltArrZ_t* tmp2 = gcd_AAZ(cont, tmpCoef);
				// AltArrZ_t* tmp2 = gcd_AAZ_tmp(cont, tmpCoef);
				freePolynomial_AAZ(cont);
				cont = tmp2;
				if (isOne_AAZ(cont)) {
					break;
				}
			}

			startIdx = i;
			memcpy(curDegs, nextDegs, nvar*sizeof(degree_t));
		}
	}
	//commit the last coef
	localCoef.size = i - startIdx;
	localCoef.alloc = localCoef.size;
	localCoef.elems = aa->elems + startIdx;
	deepCopyPolynomial_AAZ_inp_unpk(&localCoef, &tmpCoef);
	tmpDegs = (degree_t*) tmpCoef->elems->degs;

	//simultaneously copy degrees into tmpCoef meanwhile zero-ing the active variables
	for (k = 0; k < tmpCoef->size; ++k) {
		for (j = 0; j < nactive; ++j) {
			tmpDegs[k*nvar + j] = 0;
		}
		for (j = nactive; j < nvar; ++j) {
			tmpDegs[k*nvar + j] = degs[(startIdx+k)*nvar + j];
		}
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

	if (!isOne_AAZ(cont)
		&& mpz_sgn(aa->elems->coef) < 0
		&& mpz_sgn(cont->elems->coef) > 0) {
		negatePolynomial_AAZ(cont);
	}

	if (maxActiveIdx >= nactive) {
		freePolynomial_AAZ((AltArrZ_t*) aa);
		reorderVars_AAZ(cont, reverseMap, nvar);
	}
	return cont;
}




AltArrZ_t* univariateGCD_AAZ_unpk(AltArrZ_t* a, AltArrZ_t* b) {
	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMZP ERROR: Calling univariate GCD on multivariate arrays\n");
		exit(1);
	}

	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return univariateGCD_AAZ(a, b);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AAZ_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AAZ_inp(b);
		unpackedB = 1;
	}

	AltArrZ_t* r0 = NULL;
	AltArrZ_t* r1 = NULL;
	AltArrZ_t* r2 = NULL;

	mpz_t c0;
	mpz_t c1;
	mpz_init(c0);
	mpz_init(c1);

	if (isGreaterExponentVectors_unpk(a->elems->degs, b->elems->degs, 1)) {
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

	if (unpackedA) {
		packExponentVectors_AAZ_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AAZ_inp(b);
	}

	return r0;

}

AltArrZ_t* commonFactor_AAZ_unpk(const AltArrZ_t* a, AltArrZ_t** factored) {
	if (a == NULL || a->size == 0) {
		return NULL;
	}

	if (!a->unpacked) {
		return commonFactor_AAZ(a, factored);
	}

	AltArrZ_t* ret = makePolynomial_AAZ_unpk(1, a->nvar);
	mpz_init(ret->elems->coef);
	mpz_set_ui(ret->elems->coef, 1ul);

	int nvar = a->nvar;
	degree_t* retDegs = (degree_t*) ret->elems->degs;
	degree_t* aDegs = (degree_t*) a->elems->degs;
	setExponentVector_unpk(ret->elems->degs, a->elems[a->size-1].degs, nvar);
	for (int i = a->size-2; i >= 0; --i) {
		for (int j = 0; j < nvar; ++j) {
			if ( aDegs[i*nvar + j] < retDegs[j]) {
				retDegs[j] = aDegs[i*nvar + j];
			}
		}
		if (isZeroExponentVector_unpk(ret->elems->degs, nvar)) {
			break;
		}
	}
	ret->size = 1;

	if (factored != NULL) {
		AltArrZ_t* factRet = deepCopyPolynomial_AAZ_unpk(a);
		if (! isZeroExponentVector_unpk(ret->elems->degs, nvar)) {
			for (int i = 0; i < a->size; ++i) {
				subtractExponentVectors_unpk(factRet->elems[i].degs, ret->elems->degs, factRet->elems[i].degs, nvar);
			}
		}
		*factored = factRet;
	}

	return ret;
}


//Add in a variable at index varIdx, pushing the variable at that index
//and the following variables to the right in the exponent vector
void homogenizePolynomial_AAZ_inp_unpk(AltArrZ_t* aa, int varIdx) {

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


	if (!aa->unpacked) {
		unpackExponentVectors_AAZ_inp(aa);
	}

	int size = aa->size;
	degree_t tdeg = totalDegree_AAZ(aa);

	int nvar = aa->nvar;
	int newNvar = aa->nvar+1;

	//expand exponent vectors out of place
	degree_t* tmpDegs = (degree_t*) malloc(sizeof(degree_t)*newNvar*aa->alloc); //must use alloc! not size
	degree_t* degs = (degree_t*) aa->elems->degs;

	degrees_t curTDeg;
	int j;
	for (int i = 0; i < aa->size; ++i) {
		curTDeg = 0;

		//re-pack the variables [0,varIdx) to the same indices
		for (j = 0; j < varIdx; ++j) {
			tmpDegs[i*newNvar + j] = degs[i*nvar + j];
			curTDeg += degs[i*nvar + j];
		}

		//re-pack the variables [varIdx,nvar) to be shifted right
		for (j = varIdx; j < nvar; ++j) {
			tmpDegs[i*newNvar + j + 1] = degs[i*nvar + j];
			curTDeg += degs[i*nvar + j];
		}

		//set the homogenizing variable's degree for this monomial
		curTDeg = tdeg - curTDeg;
		tmpDegs[i*newNvar + varIdx] = curTDeg;

		//replace old degs pointer to new array
		aa->elems[i].degs = (degrees_t) tmpDegs + i*newNvar;
	}

	free(degs);
	aa->nvar = newNvar;
	mergeSortPolynomial_AAZ_unpk(aa);
}



/*******************
 * Evaluate interpolate
 *******************/

void univarEvaluate_AAZ_unpk(const AltArrZ_t* aa, const mpz_t point, mpz_t res) {
	if (aa->size == 0) {
		mpz_set_si(res, 0l);
		return;
	}

	if (aa->nvar != 1) {
		fprintf(stderr, "Poly is not univariate in univarEvaluate_AA_unpk\n");
		exit(1);
	}

	if (!aa->unpacked) {
		return univarEvaluate_AAZ(aa, point, res);
	}

	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	mpz_set(res, elems->coef);
	degree_t* degs = (degree_t*) elems->degs;
	degree_t prevDeg = degs[0];
	degree_t nextDeg;
	for (int i = 1; i < size; ++i) {
		nextDeg = degs[i];
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

void genPoly_univgcdheu_AAZ_inp_unpk(AltArrZ_t* gcd, mpz_t eps, mpz_t halfeps, mpz_t gamma) {
	int nvar = gcd->nvar;
	int varIdx = nvar-1;
	AAZElem_t* elems = gcd->elems;
	degree_t* degs = (degree_t*) gcd->elems->degs;
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
			degs[curSize*nvar + varIdx] = i;
			elems[curSize].degs = (degrees_t) (degs + curSize*nvar);

			++curSize;
			if (curSize >= gcd->alloc) {
				//setting the size is very important!!
				//otherwise, resize will not re-assign elems[i].degs to new degs array
				gcd->size = curSize;
				resizePolynomial_AAZ_unpk(gcd, gcd->alloc << 1);
				elems = gcd->elems;
				degs = (degree_t*) gcd->elems->degs;
			}
			mpz_init(elems[curSize].coef);
		}

	}
	mpz_clear(elems[curSize].coef);

	gcd->size = curSize;
	curSize >>= 1;

	degree_t tmp[nvar];
	degrees_t tmpDegs = (degrees_t) tmp;

	//reverse array to make degs in right order; curSize is now half the size of gcd
	for (int i = 0; i < curSize; ++i ){
		mpz_swap(elems[i].coef, elems[gcd->size - 1 - i].coef);
		//swap the value inside the degs array, not the elems[i].degs pointers
		setExponentVector_unpk(tmpDegs, elems[i].degs, nvar);
		setExponentVector_unpk(elems[i].degs, elems[gcd->size - 1 - i].degs, nvar);
		setExponentVector_unpk(elems[gcd->size - 1 - i].degs, tmpDegs, nvar);
	}
}

//if varIdx = nvar - 1, g must be at least an allocated polynomial with nvar set correctly.
//otherwise, g is gamma that we are updating.
void genPoly_gcdheu_AAZ_inp_unpk(AltArrZ_t* gamma, int varIdx, mpz_t eps, mpz_t halfeps) {
	int nvar = gamma->nvar;
	if (varIdx + 1 == nvar) {
		mpz_t gammaZ;
		mpz_init_set(gammaZ, gamma->elems->coef);
		mpz_clear(gamma->elems->coef); //avoid double init
		genPoly_univgcdheu_AAZ_inp_unpk(gamma, eps, halfeps, gammaZ);
		mpz_clear(gammaZ);
		return;
	}

	mpz_t infn;
	mpz_init(infn);
	infinityNorm_AAZ(gamma, infn);
	long int exp1, exp2;
	double d1 = mpz_get_d_2exp(&exp1, infn);
	double d2 = mpz_get_d_2exp(&exp2, eps);
	d1 = log(d1)  + log(2) * (double) exp1;
	d2 = log(d1)  + log(2) * (double) exp2;
	int maxDeg = (int) ceil(d1/d2);
	maxDeg += 2; //buffer room
	mpz_clear(infn);

	degree_t* gammaDegs = (degree_t*) gamma->elems->degs;

	AAZElem_t* chunks[maxDeg];
	AAZElem_t* chunk;
	degree_t* degChunks[maxDeg];
	degree_t* degChunk;
	int chunkSizes[maxDeg];
	int totalSize = 0;

	int exp, i, gammaSize;
	int nextGammaSize = gamma->size;
	int insertIdx = 0;
	for (exp = 0; nextGammaSize > 0; ++exp) {
		gammaSize = nextGammaSize;
		nextGammaSize = 0;

		//+1 because we mpz_init one past where we end up.
		chunk = (AAZElem_t*) malloc(sizeof(AAZElem_t)*(gammaSize+1));
		chunks[exp] = chunk;
		degChunk = (degree_t*) calloc(gammaSize*nvar, sizeof(degree_t));
		degChunks[exp] = degChunk;

		//g_exp = chunk[exp] = symmetric mod gamma.
		//simultaneously do gamma = (gamma - g_i) / eps
		//simultaneously multiply g_exp by x_{varidx}^exp
		insertIdx = 0;
		mpz_init(chunk[insertIdx].coef);
		chunk[insertIdx].degs = (degrees_t) (degChunk + insertIdx*nvar);
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
			memcpy(degChunk + insertIdx*nvar, gammaDegs + i*nvar, sizeof(degree_t)*nvar);
			degChunk[insertIdx*nvar + varIdx] = exp;

			if (mpz_sgn(chunk[insertIdx].coef) != 0) {
				++insertIdx;
				mpz_init(chunk[insertIdx].coef);
				chunk[insertIdx].degs = (degrees_t) (degChunk + insertIdx*nvar);
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
	free(gammaDegs);

	//we now have valid chunks[i] for i < exp;
	//the chunks need to be reversed in order for canonical form
	AAZElem_t* finalGamma = (AAZElem_t*) malloc(sizeof(AAZElem_t)*totalSize);
	degree_t* finalGammaDegs = (degree_t*) malloc(sizeof(degree_t)*totalSize*nvar);
	gammaDegs = finalGammaDegs; //hold head of the array in gammaDegs

	gamma->elems = finalGamma;
	gamma->alloc = totalSize;

	for (i = exp-1; i >= 0; --i) {
		//do a memcpy on elems so we get coefs cheaply, will have to redo degs later
		memcpy(finalGamma, chunks[i], sizeof(AAZElem_t)*chunkSizes[i]);
		memcpy(finalGammaDegs, degChunks[i], sizeof(degree_t)*chunkSizes[i]*nvar);
		finalGamma += chunkSizes[i];
		finalGammaDegs += chunkSizes[i]*nvar;
		//since we memcpy from chunks, we don't clear the coefs
		free(chunks[i]);
		free(degChunks[i]);
	}
	gamma->size = totalSize;

	//redoing degs since it's now later
	for (i = 0; i < totalSize; ++i) {
		gamma->elems[i].degs = (degrees_t) (gammaDegs + i*nvar);
	}
}


AltArrZ_t* convertFromDUZP_KS_AAZ_unpk(const DUZP_t* fd, const degree_t* maxDegs, int nvar) {
	if (isZero_DUZP(fd)) {
		return NULL;
	}

	int fAlloc = 50;
	AltArrZ_t* f = makePolynomial_AAZ_unpk(fAlloc, nvar);

	int curSize = 0;
	int i, j;

	degree_t maxAlpha = 1;
	for (j = 0; j < nvar; ++j) {
		maxAlpha *= (maxDegs[j] + 1);
	}

	degree_t alpha;
	degree_t r;
	degree_t* degs = (degree_t*) f->elems->degs;
	for (i = fd->lt; i >= 0; --i) {
		if (mpz_sgn(fd->coefs[i]) == 0) {
			continue;
		}

		if (curSize >= fAlloc) {
			fAlloc <<= 1;
			if (fAlloc > fd->lt) {
				fAlloc = fd->lt + 1;
			}
			f->size = curSize;
			resizePolynomial_AAZ_unpk(f, fAlloc);
			degs = (degree_t*) f->elems->degs;
		}

		alpha = maxAlpha;
		r = i;
		for (j = 0; j < nvar; ++j) {
			//update alpha for current variable
			alpha /= (maxDegs[j] + 1);

			//get partial degree by reversing Kronecker sub
			degs[curSize*nvar + j] = r / alpha;
			r %= alpha;
		}

		f->elems[curSize].degs = (degrees_t) (degs + curSize*nvar);
		mpz_init_set(f->elems[curSize].coef, fd->coefs[i]);
		++curSize;
	}

	f->size = curSize;

	return f;
}