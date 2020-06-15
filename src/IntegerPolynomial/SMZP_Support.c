
#include "IntegerPolynomial/SMZP_Support.h"
#include "IntegerPolynomial/DUZP_Support.h"

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

void printPoly_AAZ(FILE* fp, const AltArrZ_t* aa, const char** syms, int nvar) {
	if (aa == NULL || aa->size == 0) {
		fprintf(stderr, "0");
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
		printDegs_AA(fp, aa->elems[0].degs, syms, nvar, masks, sizes);
	}
	for (int i = 1; i < AA_SIZE(aa)-1; ++i) {
		if (mpz_sgn(aa->elems[i].coef) > 0) {
			fprintf(fp, " + ");
		} else {
			fprintf(fp, " ");
		}
		gmp_fprintf(fp, "%Zd", aa->elems[i].coef);
		if (!isZeroExponentVector(aa->elems[i].degs)) {
			fprintf(fp, "*");
			printDegs_AA(fp, aa->elems[i].degs, syms, nvar, masks, sizes);
		}
	}
	if (aa->size > 1) {
		if (mpz_sgn(aa->elems[aa->size-1].coef) > 0) {
			fprintf(fp, " + ");
		} else {
			fprintf(fp, " ");
		}
		gmp_fprintf(fp, "%Zd", aa->elems[aa->size-1].coef);
		if (!isZeroExponentVector(aa->elems[aa->size-1].degs)) {
			fprintf(fp, "*");
			printDegs_AA(fp, aa->elems[aa->size-1].degs, syms, nvar, masks, sizes);
		}
	}

	//fprintf(fp, "\n");

}

void calculateMaxDegsList_AAZ(const AltArrZ_t* aa, degree_t* maxList) {
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

void nonZeroVariables_AAZ(AltArrZ_t* aa, int* foundVar) {
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

degree_t mainDegree_AAZ(AltArrZ_t* aa) {
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
				//TODO
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
				//TODO
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
	if (aa == NULL || AA_SIZE(aa) == 0) {
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

	if (aa->unpacked || (bb != NULL && (*bb)->unpacked)) {
		deepCopyPolynomial_AAZ_inp_unpk(aa, bb);
		return;
	}

	AltArrZ_t* b = *bb;
	if (b == NULL) {
		b = makePolynomial_AAZ(aa->size, aa->nvar);
		*bb = b;
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
	AAZElem_t* tempElems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*size);
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

AltArrZ_t* multiplyByInteger_AAZ(AltArrZ_t* a, mpz_t mult) {
	if (a == NULL || a->size == 0) {
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
	if (a == NULL || a->size == 0) {
		return NULL;
	}

	AltArrZ_t* ret = deepCopyPolynomial_AAZ(a);
	AAZElem_t* retElems = ret->elems;
	for (int i = 0; i < ret->size; ++i) {
		if (!mpz_divisible_p(retElems[i].coef, div)) {
			fprintf(stderr, "SMZP divideByInteger ERROR!\n");
			exit(1);
		}
		mpz_div(retElems[i].coef, retElems[i].coef, div);
	}

	return ret;
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

    AltArrZ_t* poly = makePolynomial_AAZ (coefSize, nvar);
    poly->size = coefSize;
    poly->elems = coef;
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

int leadingVariable_AAZ (AltArrZ_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return -2;
	}

	if (aa->unpacked) {
		return leadingVariable_AAZ_unpk(aa);
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


AltArrZ_t* mainLeadingCoefficient_AAZ (AltArrZ_t* aa)
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
		//TODO!!!
		//for now just call old method for hard exit
		// checkValidMonomialMult(aMax,bMax,nvar);


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

void divideBySingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar) {
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
void dividePolynomials_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, register int nvar) {

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


int divideTestSingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {
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

int divideTest_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {

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
		freePolynomial_AAZ_unpk(aa);
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

AltArrZ_t* univariateHeuristicGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b) {
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
	a = primitivePartAndContent_AAZ(a, conta);
	b = primitivePartAndContent_AAZ(b, contb);
	mpz_gcd(conta, conta, contb);

	degrees_t maxDegs = (isLessExponentVectors(a->elems->degs, b->elems->degs) ? b->elems->degs : a->elems->degs);
	mpz_t infa, infb, eps, halfeps;
	mpz_inits(infa, infb, eps, halfeps, NULL);
	infinityNorm_AAZ(a, infa);
	infinityNorm_AAZ(b, infb);

	if (mpz_cmp(infa, infb) < 0) {
		mpz_set(infa, infb);
	}
	mpz_mul_si(infa, infa, 2l);
	mpz_add_ui(infa, infa, 2ul);
	mpz_swap(eps, infa);
	mpz_fdiv_q_ui(halfeps, eps, 2ul);
	long epsSize = mpz_sizeinbase(eps, 2);

	AltArrZ_t* gcd = NULL;
	for (int attempts = 0; attempts < 6; ++attempts) {
		if (maxDegs * epsSize > 5000) {
			gcd = NULL;
			break;
		}
		univarEvaluate_AAZ(a, eps, infa);
		univarEvaluate_AAZ(b, eps, infb);
		mpz_gcd(infa, infa, infb);


		gcd = makePolynomial_AAZ(10, 1);
		AAZElem_t* elems = gcd->elems;
		int curSize = 0;
		mpz_init(elems[curSize].coef);
		for(int i = 0; mpz_sgn(infa) != 0; ++i) {
			mpz_mod(elems[curSize].coef, infa, eps);
			if (mpz_cmp(elems[curSize].coef, halfeps) > 0) {
				mpz_sub(elems[curSize].coef, elems[curSize].coef, eps);
			}

			mpz_sub(infa, infa, elems[curSize].coef);
			mpz_fdiv_q(infa, infa, eps);

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

		primitivePart_AAZ_inp(gcd);
		if (divideTest_AAZ(a, gcd, NULL, 1) && divideTest_AAZ(b, gcd, NULL, 1)) {
			break;
		}

		freePolynomial_AAZ(gcd);
		gcd = NULL;

		//update eps;
		mpz_mul_si(eps, eps, 73794l);
		mpz_fdiv_q_ui(eps, eps, 27011ul);
		epsSize = mpz_sizeinbase(eps, 2);
	}

	multiplyByInteger_AAZ_inp(gcd, conta);
	mpz_clears(conta, contb, infa, infb, eps, NULL);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	return gcd;
}

AltArrZ_t* univariateHeuristicPrimitiveGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b) {
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
	mpz_mul_si(infa, infa, 2l);
	mpz_add_ui(infa, infa, 2ul);
	mpz_swap(eps, infa);
	mpz_fdiv_q_ui(halfeps, eps, 2ul);
	long epsSize = mpz_sizeinbase(eps,2);

	AltArrZ_t* gcd = NULL;
	for (int attempts = 0; attempts < 6; ++attempts) {
		if (maxDegs * epsSize > 5000) {
			gcd = NULL;
			break;
		}
		univarEvaluate_AAZ(a, eps, infa);
		univarEvaluate_AAZ(b, eps, infb);
		mpz_gcd(infa, infa, infb);


		gcd = makePolynomial_AAZ(10, 1);
		AAZElem_t* elems = gcd->elems;
		int curSize = 0;
		mpz_init(elems[curSize].coef);
		for(int i = 0; mpz_sgn(infa) != 0; ++i) {
			mpz_mod(elems[curSize].coef, infa, eps);
			if (mpz_cmp(elems[curSize].coef, halfeps) > 0) {
				mpz_sub(elems[curSize].coef, elems[curSize].coef, eps);
			}

			mpz_sub(infa, infa, elems[curSize].coef);
			mpz_fdiv_q(infa, infa, eps);

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

		primitivePart_AAZ_inp(gcd);
		if (divideTest_AAZ(a, gcd, NULL, 1) && divideTest_AAZ(b, gcd, NULL, 1)) {
			break;
		}

		freePolynomial_AAZ(gcd);
		gcd = NULL;

		//update eps;
		mpz_mul_si(eps, eps, 73794l);
		mpz_fdiv_q_ui(eps, eps, 27011ul);
		epsSize = mpz_sizeinbase(eps, 2);
	}

	mpz_clears(infa, infb, eps, NULL);
	return gcd;
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
	AltArrZ_t* primA = deepCopyPolynomial_AAZ(a);
	primitivePartAndContent_AAZ_inp(primA, conta);
	AltArrZ_t* primB = deepCopyPolynomial_AAZ(b);
	primitivePartAndContent_AAZ_inp(primB, contb);
	mpz_gcd(conta, conta, contb);

	AltArrZ_t* gAA = univariateHeuristicPrimitiveGCD_AAZ(primA, primB);
	if (gAA == NULL) {
		DUZP_t* aDUZP = convertFromAltArrZ_DUZP(primA);
		DUZP_t* bDUZP = convertFromAltArrZ_DUZP(primB);
		DUZP_t* gDUZP = primitiveGCD_DUZP(aDUZP, bDUZP);
		gAA = convertToAltArrZ_DUZP(gDUZP);
		freePolynomial_DUZP(aDUZP);
		freePolynomial_DUZP(bDUZP);
		freePolynomial_DUZP(gDUZP);
	}

	multiplyByInteger_AAZ_inp(gAA, conta);

	freePolynomial_AAZ(primA);
	freePolynomial_AAZ(primB);
	mpz_clear(conta);
	mpz_clear(contb);

	return gAA;

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



/*****************
 * Evaluation
 *****************/

void univarEvaluate_AAZ(AltArrZ_t* aa, const mpz_t point, mpz_t res) {
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

