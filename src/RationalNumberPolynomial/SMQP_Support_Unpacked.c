#include "RationalNumberPolynomial/SMQP_Support_Unpacked.h"
#include "RationalNumberPolynomial/SMQP_Support-AA.h"

AltArr_t* makePolynomial_AA_unpk(int allocSize, int nvar) {
	if (allocSize < 1) {
		return NULL;
	}

	AltArr_t* newAA = (AltArr_t*) malloc(sizeof(AltArr_t));
	newAA->size = 0; 
	newAA->alloc = allocSize;
	newAA->nvar = nvar;
	newAA->unpacked = 1;
	newAA->elems = (AAElem_t*) malloc(sizeof(AAElem_t)*allocSize);
	degree_t* degs = (degree_t*) calloc(allocSize*nvar, sizeof(degree_t));
	newAA->elems[0].degs = (degrees_t) degs;
	return newAA;
}

void freePolynomial_AA_unpk(AltArr_t* aa) {
	freePolynomial_AA(aa);
}

degrees_t calculateMaxDegs_AA_unpk(const AltArr_t* aa) {

	AAElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	degree_t* maxList = (degree_t*) calloc(nvar, sizeof(degree_t));
	degree_t* degs_unpk;
	for (register int i = 0; i < size; ++i) {
		degs_unpk = (degree_t*) elems[i].degs;
		for (register int j = 0; j < nvar; ++j) {
			maxList[j] = (degs_unpk[j]) > (maxList[j]) ? degs_unpk[j] : maxList[j]; 
		}
	}

	return (degrees_t) maxList;
}

void unpackExponentVectors_AA_inp(AltArr_t* aa) {
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

	AAElem_t* elems = aa->elems;
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
	for (int i = size; i < aa->alloc; ++i) {
		elems[i].degs = (degrees_t) &(unpackedDegs[i*nvar]);
	}
}

int monomialDivideTest_AA_unpk(AltArr_t* a, int idxa, AltArr_t* b, int idxb) {
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

void tryPackExponentVectors_AA_inp(AltArr_t* aa) {
	if (aa == NULL || aa->alloc == 0 || !aa->unpacked) {
		return;
	}
	degree_t* maxDegs = (degree_t*) calculateMaxDegs_AA_unpk(aa);
	if (isPackableExponentVector(maxDegs, aa->nvar)) {
		packExponentVectors_AA_inp(aa);
	}
	free(maxDegs);
}

void packExponentVectors_AA_inp(AltArr_t* aa) {
	if (aa == NULL || aa->alloc == 0 || !aa->unpacked) {
		return;
	}
	int nvar = aa->nvar;
	int size = aa->size;

	AAElem_t* elems = aa->elems;
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
				fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for packExponentVectors_AA_inp at index %d; %d > %lld.\n", k, degs_unpk[k], maxExps[k]);
				// free((int*) -1);
				exit(1);
			}
			elems[i].degs |= (((degrees_t) degs_unpk[k]) << sizes[k]);
		}
	}

	int alloc = aa->alloc;
	for (int i = size; i < alloc; ++i) {
		elems[i].degs = 0;
	}

	aa->unpacked = 0;

	free(freeDegs);
}

void printDegs_AA_unpk(FILE* fp, degrees_t degs, const char** syms, int nvar) {
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

void printPoly_AA_unpk(FILE* fp, const AltArr_t* aa, const char** syms, int nvar) {
	if (aa == NULL || aa->size == 0) {
		fprintf(stderr, "0\n");
		return;
	}

	if (!aa->unpacked) {
		printPoly_AA(fp, aa, syms, nvar);
		return;
	}
	
	gmp_fprintf(fp, "%Qd", aa->elems[0].coef);
	if (!isZeroExponentVector_unpk(aa->elems[0].degs, nvar)) {
		fprintf(fp, "*");
		printDegs_AA_unpk(fp, aa->elems[0].degs, syms, nvar);
	}
	for (int i = 1; i < AA_SIZE(aa)-1; ++i) {
		if (mpq_sgn(aa->elems[i].coef) > 0) {
			fprintf(fp, " + ");
		} else {
			fprintf(fp, " ");
		}
		gmp_fprintf(fp, "%Qd", aa->elems[i].coef);
		if (!isZeroExponentVector_unpk(aa->elems[i].degs, nvar)) {
			fprintf(fp, "*");
			printDegs_AA_unpk(fp, aa->elems[i].degs, syms, nvar);
		}
	}
	if (aa->size > 1) {
		if (mpq_sgn(aa->elems[aa->size-1].coef) > 0) {
			fprintf(fp, " + ");
		} else {
			fprintf(fp, " ");
		}
		gmp_fprintf(fp, "%Qd", aa->elems[aa->size-1].coef);
		if (!isZeroExponentVector_unpk(aa->elems[aa->size-1].degs, nvar)) {
			fprintf(fp, "*");
			printDegs_AA_unpk(fp, aa->elems[aa->size-1].degs, syms, nvar);
		}
	}
	
	//fprintf(fp, "\n");	
}

int isExactlyEqual_AA_unpk(AltArr_t* a, AltArr_t* b) {
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
		return isExactlyEqual_AA(a,b);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
		unpackedA = 1;		
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);		
		unpackedB = 1;		
	}

	int cmp = 1;
	int nvar = a->nvar;
	for (int i = 0; i < a->size; ++i) {
		if (!isEqualExponentVectors_unpk(a->elems[i].degs, b->elems[i].degs, nvar)) {
			cmp = 0;
			break;
		}
		if (mpq_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
			cmp = 0;
			break;
		}
	}

	if (unpackedA) {
		packExponentVectors_AA_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	return cmp;
}

void nonZeroVariables_AA_unpk(AltArr_t* aa, int* foundVar) {
    if (aa == NULL || aa->nvar == 0 || isConstant_AA(aa)) {
        return;
    }

    if (!aa->unpacked) {
    	nonZeroVariables_AA(aa, foundVar);
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

degree_t totalDegree_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0 || isConstant_AA(aa)) {
		return 0;
	}

	if (!aa->unpacked) {
		return totalDegree_AA(aa);
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

degree_t partialDegree_AA_unpk(const AltArr_t* aa, int k) {
	if (aa == NULL || aa->size == 0 || k >= aa->nvar) {
		return -1;
	}

	if (!aa->unpacked) {
		return partialDegree_AA(aa, k);
	}

	int nvar = aa->nvar;
	degree_t dMax = 0;
	degree_t* degs = (degree_t*) aa->elems->degs;
	for (int i = 0; i < aa->size; ++i) {
		dMax = degs[i*nvar + k] > dMax ? degs[i*nvar + k] : dMax; 
	}

	return dMax;
}

void partialDegrees_AA_unpk(const AltArr_t* aa, degree_t* degsList) {
	if (aa == NULL || aa->size == 0 || aa->nvar == 0) {
		return;
	}

	if (!aa->unpacked) {
		partialDegrees_AA(aa, degsList);
		return;
	}

	degrees_t maxDegs = calculateMaxDegs_AA_unpk(aa);
	memcpy(degsList, (degree_t*) maxDegs, sizeof(degree_t)*aa->nvar);
	free((degree_t*) maxDegs);
}

degree_t mainDegree_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return -1;
	}

	if (isConstant_AA(aa)) {
		return 0;
	}

	if (!aa->unpacked) {
		return mainDegree_AA(aa);
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

int mainVariable_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->nvar == 0) {
		return -1;
	}

	if (!aa->unpacked) {
		return mainVariable_AA(aa);
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

void coefficient_AA_unpk(AltArr_t* aa, const degree_t* degsList, int nvar, mpq_t retCoef) {
	if (isZero_AA(aa)) {
		mpq_set_si(retCoef, 0l, 1l);
        return;
    }
    if (nvar == 0 && aa->nvar == 0) {
        mpq_set(retCoef, aa->elems->coef);
    	return;
	}

	if (!aa->unpacked) {
		return coefficient_AA(aa, degsList, nvar, retCoef);
	}

	degrees_t degs = (degrees_t) degsList;
	mpq_set_ui(retCoef, 0l, 1l);
    for (int i = 0; i < aa->size; ++i) {
        if (isEqualExponentVectors_unpk(degs, aa->elems[i].degs, nvar)) {
        	mpq_set(retCoef, aa->elems[i].coef);
            break;
        }
        if (isLessExponentVectors_unpk(aa->elems[i].degs, degs, nvar)) {
            //due to ordering of terms, none will have degs as monomial.
            break;
        }
    }
}

void setCoefficient_AA_unpk(AltArr_t* aa, const degree_t* degsList, int nvar, const mpq_t coef) {
    if (aa == NULL || nvar != aa->nvar) {
    	return;
    }

    if (aa->size == 0 || nvar == 0) {
    	if (aa->alloc < 1) {
    		resizePolynomial_AA(aa, 1);
    	}
        mpq_set(aa->elems->coef, coef);
    	return;
    }

    if (!aa->unpacked) {
    	setCoefficient_AA(aa, degsList, nvar, coef);
    	return;
    }

    degrees_t degs = (degrees_t) degsList;
    for (int i = 0; i < aa->size; ++i) {
        if (isEqualExponentVectors_unpk(aa->elems[i].degs, degs, nvar)) {
            mpq_set(aa->elems[i].coef, coef);
            return;
        }
        if (isLessExponentVectors_unpk(aa->elems[i].degs, degs, nvar)) {
            //shift everything to the right 1 and insert at i.
            if (aa->size >= aa->alloc) {
                resizePolynomial_AA_unpk(aa, aa->size + 1);
            }
            
            degree_t* aaDegs = (degree_t*) aa->elems->degs;
            aa->elems[aa->size].degs = (degrees_t) (aaDegs + aa->size*nvar); 
            mpq_init(aa->elems[aa->size].coef);

            for (int j = aa->size; j > i; --j) {
                mpq_swap(aa->elems[j].coef, aa->elems[j-1].coef);
                setExponentVector_unpk(aa->elems[j].degs, aa->elems[j-1].degs, nvar);
            }
            setExponentVector_unpk(aa->elems[i].degs, degs, nvar);
            aa->elems[i].degs = degs;
            mpq_set(aa->elems[i].coef, coef);
            ++(aa->size);
            return; 
        }
    }

    //if we get here then we need to insert at the end of the array;
    if (aa->size >= aa->alloc) {
        resizePolynomial_AA_unpk(aa, aa->size + 1);
    }
    degree_t* aaDegs = (degree_t*) aa->elems->degs;
    aa->elems[aa->size].degs = (degrees_t) (aaDegs + aa->size*nvar); 

    mpq_init(aa->elems[aa->size].coef);
    mpq_set(aa->elems[aa->size].coef, coef);
    setExponentVector_unpk(aa->elems[aa->size].degs, degs, nvar);
    aa->elems[aa->size].degs = degs;
    ++(aa->size);
}

/**
 * Append a new term ot the end of the AltArr_t aa. 
 * This operations is potentially unsafe; no re-ordering is performed, 
 * it is only a simple append.
 */
void addTerm_AA_unpk(AltArr_t* aa, const degree_t* degsList, int nvar, const mpq_t coef) {
	if (aa == NULL || aa->nvar != nvar) {
		return;
	}

	if (!aa->unpacked) {
		addTerm_AA(aa, degsList, nvar, coef);
		return;
	}

 	if (aa->size >= aa->alloc) {
    	resizePolynomial_AA_unpk(aa, aa->alloc*2);
    }

	degree_t* aaDegs = (degree_t*) aa->elems->degs;
    aa->elems[aa->size].degs = (degrees_t) (aaDegs + aa->size*nvar); 

    mpq_init(aa->elems[aa->size].coef);
    mpq_set(aa->elems[aa->size].coef, coef);
    setExponentVector_unpk(aa->elems[aa->size].degs, (degrees_t) degsList, nvar);
    ++(aa->size);
}

int isEqualWithVariableOrdering_AA_unpk(AltArr_t* a, AltArr_t* b, const int* xs, int xsSize) {
	if (isZero_AA(a)) {
		if (isZero_AA(b)) {
			return 1;
		} else {
			return 0;
		}
	}
	if (isZero_AA(b)) {
		return 0;
	}

	if (a->size != b->size) {
        return 0;
    }

	if (!a->unpacked && !b->unpacked) {
		return isEqualWithVariableOrdering_AA(a, b, xs, xsSize);
	}

	int nvar = a->nvar;
	int bnvar = b->nvar;
	int v = xsSize / 2;

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

    int ret = 1;
    degree_t* aDegs = (degree_t*) a->elems->degs;
    degree_t* bDegs = (degree_t*) b->elems->degs;
    degree_t adeg, bdeg;
    for (int i = 0; i < a->size; ++i){
        if (mpq_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
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
    	packExponentVectors_AA_inp(a);
    }
    if (unpackedB) {
    	packExponentVectors_AA_inp(b);
    }

    return ret;
}



void expandNumVars_AA_unpk(AltArr_t* aa, int newNvar) {
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
		expandNumVars_AA(aa, newNvar);
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

void expandNumVarsLeft_AA_unpk(AltArr_t* aa, int newNvar) {
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
		expandNumVarsLeft_AA(aa, newNvar);
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

void shrinkNumVarsAtIdx_AA_unpk(AltArr_t* aa, int idx) {
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
		shrinkNumVarsAtIdx_AA(aa, idx);
		return;
	}

	int nvar = aa->nvar;
	int newNvar = aa->nvar - 1;

	degree_t* newDegs = (degree_t*) malloc(sizeof(degree_t)*newNvar*aa->alloc);
	degree_t* curDegs = (degree_t*) aa->elems->degs;
	AAElem_t* elems = aa->elems;

	for (int i = 0; i < aa->size; ++i ) {
		if (idx > 0) {
			memcpy(newDegs + i*newNvar, curDegs + i*nvar, sizeof(degree_t)*idx);
		}
		memcpy(newDegs + i*newNvar + idx, curDegs + i*nvar + idx + 1, sizeof(degree_t)*(nvar-idx-1));
		elems[i].degs = (degrees_t) (newDegs + i*newNvar);
	}

	aa->nvar = newNvar;

	free(curDegs);

	tryPackExponentVectors_AA_inp(aa);
}

void shrinkAndReorderVars_AA_unpk(AltArr_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->nvar < 1) {
		return;
	}

	if (varmapSize > aa->nvar) {
		return;
	} 

	if (!aa->unpacked) {
		shrinkAndReorderVars_AA(aa, varMap, varmapSize);
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
			mpq_clear(aa->elems[i].coef);
		}
		free((degree_t*) aa->elems->degs);
		aa->elems->degs = 0ll;
		aa->size = 1;
		aa->nvar = 0;
		return;
	}

	int nvar = aa->nvar;
	AAElem_t* elems = aa->elems;
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
		mergeSortPolynomial_AA(aa);
	}

}

//     degs[varmap[i]] = curDegs[i];
void reorderVars_AA_unpk(AltArr_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (!aa->unpacked) {
		reorderVars_AA(aa, varMap, varmapSize);
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

	tryPackExponentVectors_AA_inp(aa);

	mergeSortPolynomial_AA(aa);
}

void setDegrees_AA_inp_unpk(AltArr_t* aa, int idx, const degree_t* degsList, int nvar) {
	if (aa == NULL || idx > aa->alloc) {
		return;
	}

	if (!aa->unpacked) {
		return setDegrees_AA_inp(aa, idx, degsList, nvar);
	}

	degree_t* degs = (degree_t*) aa->elems[idx].degs;
	for (int k = 0; k < nvar; ++k) {
		degs[k] = degsList[k];
	}
}

void resizePolynomial_AA_unpk(AltArr_t* aa, int allocSize) {
	if (aa == NULL || !aa->unpacked) {
		return resizePolynomial_AA(aa, allocSize);
	}

	if (allocSize < aa->size) {
		for (int i = aa->size; i >= allocSize; --i) {
			mpq_clear(aa->elems[i].coef);
		}
		aa->size = allocSize;
	}

	if (allocSize <= 0) {
		//just don't update allocsize like why are doing this non-sense resize?
		return;
	}

	aa->elems = (AAElem_t*) realloc(aa->elems, sizeof(AAElem_t)*allocSize);
	aa->alloc = allocSize;
	
	degree_t* degsHead = (degree_t*) aa->elems[0].degs;
	degsHead = realloc(degsHead, sizeof(degree_t)*aa->nvar*allocSize);
	aa->elems[0].degs = (degrees_t) degsHead;
	int size = aa->size;
	int nvar = aa->nvar;
	AAElem_t* elems = aa->elems;
	for (int i = 0; i < size; ++i) {
		elems[i].degs = (degrees_t) (degsHead + i*nvar);
	}
}

AltArrDegList_t* deepCopyPolynomial_AADegListFromAA_unpk(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return deepCopyPolynomial_AADegListFromAA(aa);
	}

	int nvar = aa->nvar;

	AltArrDegList_t* ret = makePolynomial_AADL(aa->size, nvar);
	AAElem_DegList_t* rElems = ret->elems;
	AAElem_t* elems = aa->elems;
	degree_t* degs = (degree_t*) aa->elems->degs;

	for (int i = 0; i < aa->size; ++i) {
		mpq_init(rElems[i].coef);
		mpq_set(rElems[i].coef, elems[i].coef);
		rElems[i].degs = malloc(sizeof(degree_t)*nvar);
		for (int j = 0; j < nvar; ++j) {
			rElems[i].degs[j] = degs[i*nvar + j];
		}
	}

	ret->size = aa->size;

	return ret;
}

AltArr_t* deepCopyPolynomial_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		AltArr_t* ret = deepCopyPolynomial_AA(aa);
		unpackExponentVectors_AA_inp(ret);
		return ret;
	}

	if (aa->nvar == 0) {
		return makeConstPolynomial_AA(aa->size, 0, aa->elems->coef);
	}

	int nvar = aa->nvar;
	int size = aa->size;
	int nexp = nvar * aa->size;

	degree_t* unpackedDegs = (degree_t*) malloc(sizeof(degree_t) * nexp);
	degree_t* original_degs = (degree_t*) aa->elems[0].degs;
	memcpy(unpackedDegs, original_degs, sizeof(degree_t)*nexp);

	AltArr_t* ret = makePolynomial_AA(size, nvar);
	ret->size = size;

	AAElem_t* elems = aa->elems;
	AAElem_t* retelems = ret->elems;
	for (int i = 0; i < size; ++i) {
		mpq_init(retelems[i].coef);
		mpq_set(retelems[i].coef, elems[i].coef);
		retelems[i].degs = (degrees_t) &(unpackedDegs[i*nvar]);
	}

	ret->unpacked = 1;

	return ret;
}


AltArr_t* sortPolynomial_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}
	if (!aa->unpacked) {
		return sortPolynomial_AA(aa);
	}

	AAElem_t* elems = aa->elems;
	int size = AA_SIZE(aa);
	int nvar = aa->nvar;

	degree_t swapDegs[nvar];
	for (int i = 1; i < size; ++i) {
		for (int j = i; j > 0 && compareExponentVectors_unpk(elems[j-1].degs, elems[j].degs, nvar) < 0; --j) {
			mpq_swap(elems[j-1].coef, elems[j].coef);
			setExponentVector_unpk((degrees_t) swapDegs, elems[j-1].degs, nvar);
			setExponentVector_unpk(elems[j-1].degs, elems[j].degs, nvar);
			setExponentVector_unpk(elems[j].degs,(degrees_t)  swapDegs, nvar);
			// swapDegs = elems[j-1].degs;
			// elems[j-1].degs = elems[j].degs;
			// elems[j].degs = swapDegs;
		}
	}

	condensePolyomial_AA_unpk(aa);

	return aa;
}


static void mergeAAElems_unpk(AAElem_t* __restrict__ a, AAElem_t* __restrict__ endA, AAElem_t* __restrict__ b , AAElem_t* __restrict__ endB, AAElem_t* __restrict__ sorted, int nvar) {

	int i = 0;
	while (a < endA && b < endB) {
		if (isGreaterExponentVectors_unpk(a->degs, b->degs, nvar)) {
			sorted[i].coef->_mp_num = a->coef->_mp_num;
			sorted[i].coef->_mp_den = a->coef->_mp_den;
			memcpy((degree_t*) sorted[i].degs, (degree_t*) a->degs, nvar*sizeof(degree_t));
			// sorted[i] = *a; 
			++a;
		} else {
			sorted[i].coef->_mp_num = b->coef->_mp_num;
			sorted[i].coef->_mp_den = b->coef->_mp_den;
			memcpy((degree_t*) sorted[i].degs, (degree_t*) b->degs, nvar*sizeof(degree_t));
			// sorted[i] = *b;
			++b;
		}
		++i;
	}

	while (a < endA) {
		sorted[i].coef->_mp_num = a->coef->_mp_num;
		sorted[i].coef->_mp_den = a->coef->_mp_den;
		memcpy((degree_t*) sorted[i].degs, (degree_t*) a->degs, nvar*sizeof(degree_t));
		// sorted[i] = *a;
		++a;
		++i;
	}

	while (b < endB) {
		sorted[i].coef->_mp_num = b->coef->_mp_num;
		sorted[i].coef->_mp_den = b->coef->_mp_den;
		memcpy((degree_t*) sorted[i].degs, (degree_t*) b->degs, nvar*sizeof(degree_t));
		// sorted[i] = *b;
		++b;
		++i;
	}
}

void mergeSortPolynomial_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (!aa->unpacked) {
		mergeSortPolynomial_AA(aa);
		return;
	}

	if (aa->size < 8) {
		sortPolynomial_AA_unpk(aa);
		return;
	}

	int nvar = aa->nvar;
	int size = aa->size;
	AAElem_t* tempElems = (AAElem_t*) malloc(sizeof(AAElem_t)*aa->alloc);
	degree_t* tempDegs = (degree_t*) malloc(sizeof(degree_t)*aa->alloc*nvar);
	for (int i = 0; i < size; ++i) {
		tempElems[i].degs = (degrees_t) &(tempDegs[i*nvar]);
	}
	AAElem_t* elems = aa->elems;
	int end1 = 0;
	int end2 = 0;
	for (int window = 1; window < size; window <<= 1) {
		//merge elems[i]:elems[i+window] with elems[i+window]:elems[i+2*window]
		for (int i = 0; i < size; i += 2*window) {
			end1 = i + window < size ? i + window : size;
			end2 = i + 2*window < size ? i + 2*window : size;
			mergeAAElems_unpk(elems+i, elems+end1, elems+end1, elems+end2, tempElems+i, nvar);
		}

		AAElem_t* temp = tempElems;
		tempElems = elems;
		elems = temp;
	}

	aa->elems = elems;
	free((degree_t*) tempElems[0].degs);
	free(tempElems);

	condensePolyomial_AA_unpk(aa);
}

void condensePolyomial_AA_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (AA_SIZE(aa) < 1) {
		return;
	}

	if (!aa->unpacked) {
		condensePolyomial_AA(aa);
		return;
	}

	int size = AA_SIZE(aa);
	int nvar = aa->nvar;
	AAElem_t* elems = aa->elems;
	int insertIdx = 0;
	int compareIdx = 1;
	while (compareIdx < size) {
		if (compareExponentVectors_unpk(elems[insertIdx].degs, elems[compareIdx].degs, nvar) == 0) {
			mpq_add(elems[insertIdx].coef, elems[insertIdx].coef, elems[compareIdx].coef);
		} else if (mpq_sgn(elems[insertIdx].coef) == 0) {
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
			mpq_swap(elems[insertIdx].coef, elems[compareIdx].coef);
		} else {
			++insertIdx;
		}
		++compareIdx;		
	}

	if (mpq_sgn(elems[insertIdx].coef) != 0) {
		++insertIdx;
	}

	for (int i = insertIdx; i < size; ++i) {
		mpq_clear(elems[i].coef);
	}
	AA_SIZE(aa) = insertIdx;
}

void evalPolyToVal_AA_unpk(const AltArr_t* aa, ratNum_t* vals, short nvar, ratNum_t res) {
	if (aa == NULL || nvar != aa->nvar) {
		mpq_set_ui(res, 0ul, 1ul);
	}

	if (nvar == 0 || aa->nvar == 0) {
		mpq_set(res, aa->elems->coef);
		return;
	}

	if (!aa->unpacked) {
		return evalPolyToVal_AA(aa, vals, nvar, res);
	}

	mpq_t* valList[nvar];
	int valListSize[nvar];

	degree_t* maxDegs = (degree_t*) calculateMaxDegs_AA_unpk(aa);
	for (int j = 0; j < nvar; ++j) {
		valList[j] = (mpq_t*) malloc(sizeof(mpq_t)*(maxDegs[j]+1));
		valListSize[j] = maxDegs[j]+1;

		mpq_init(valList[j][0]);
		mpq_set_ui(valList[j][0], 1ul, 1ul);
		for (int k = 1; k < maxDegs[j]+1; ++k) {
			mpq_init(valList[j][k]);
			mpq_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	mpq_set_ui(res, 0ul, 1ul);
	mpq_t acc; 
	mpq_init(acc);

	int size = aa->size;
	degree_t* degs_unpk;
	for (int i = 0; i < size; ++i) {
		mpq_set(acc, aa->elems[i].coef);
		degs_unpk = (degree_t*) aa->elems[i].degs;
		for (int j = 0; j < nvar; ++j) {
			mpq_mul(acc, acc, valList[j][degs_unpk[j]]);
		}
		mpq_add(res, res, acc);
	}

	mpq_clear(acc);

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpq_clear(valList[j][k]);
		}
		free(valList[j]);
	}
	free(maxDegs);
}

AltArr_t* evaluatePoly_AA_unpk(AltArr_t* aa, int* active, ratNum_t* vals, short nvar) {
	if (aa == NULL) {
		return NULL;
	}

	if (nvar == 0 || aa->nvar == 0) {
		return deepCopyPolynomial_AA(aa);
	}

	if(!aa->unpacked) {
		return evaluatePoly_AA(aa, active, vals, nvar);
	}

	int i;
	int newNvar = nvar;
	for (i = 0; i < nvar; ++i) {
		if (active[i]) {
			--newNvar;
		}
	}

	if (newNvar == 0) {
		mpq_t res;
		mpq_init(res);
		evalPolyToVal_AA(aa, vals, nvar, res);
		AltArr_t* ret = makeConstPolynomial_AA(1, 0, res);
		mpq_clear(res);
		return ret;
	} 

	mpq_t* valList[nvar];
	int valListSize[nvar];

	degree_t* maxDegs = (degree_t*) calculateMaxDegs_AA_unpk(aa);
	for (int j = 0; j < nvar; ++j) {
		if (!active[j]) {
			valList[j] = NULL;
			valListSize[j] = 0;
			continue;
		}
		valList[j] = (mpq_t*) malloc(sizeof(mpq_t)*(maxDegs[j]+1));
		valListSize[j] = maxDegs[j]+1;

		mpq_init(valList[j][0]);
		mpq_set_ui(valList[j][0], 1ul, 1ul);
		for (int k = 1; k < maxDegs[j]+1; ++k) {
			mpq_init(valList[j][k]);
			mpq_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}
	free(maxDegs);

	degree_t* degsArr = (degree_t*) aa->elems[0].degs;

	AltArr_t* res = makePolynomial_AA_unpk(aa->size, newNvar);
	AAElem_t* resElems = res->elems;
	degree_t* newDegsArr = (degree_t*) resElems[0].degs;
	res->size = aa->size;

	int size = aa->size;
	int k = 0;
	for (int i = 0; i < size; ++i) {
		mpq_init(resElems[i].coef);
		mpq_set(resElems[i].coef, aa->elems[i].coef);
		resElems[i].degs = (degrees_t) &(newDegsArr[i*newNvar]);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = degsArr[i*nvar + j];
			if (valList[j] == NULL) {
				newDegsArr[i*newNvar + k] = deg;
				// newDegs |= deg << newSizes[k];
				++k;
			} else {
				mpq_mul(resElems[i].coef, resElems[i].coef, valList[j][deg]);
			}
		}
		k = 0;
	}

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpq_clear(valList[j][k]);
		}
		free(valList[j]);
	}

	mergeSortPolynomial_AA_unpk(res);

	tryPackExponentVectors_AA_inp(res);

	return res;
}


AltArr_t* mainLShiftPolynomial_AA_unpk (AltArr_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (n < 1){
		return deepCopyPolynomial_AA (aa);
    }
    
    if (!aa->unpacked) {
    	return mainLShiftPolynomial_AA(aa, n);
    }

    AltArr_t* shifted = deepCopyPolynomial_AA_unpk(aa);
    return mainLShiftPolynomial_AA_inp_unpk(shifted, n);
}


AltArr_t* mainLShiftPolynomial_AA_inp_unpk (AltArr_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (n < 1){
		return aa;
    }
    
    if (!aa->unpacked) {
    	return mainLShiftPolynomial_AA_inp(aa, n);
    }

    int nvar = aa->nvar;
    degree_t* degs = (degree_t*) aa->elems->degs;
    for (int i = 0; i < aa->size; ++i){
    	degs[i*nvar] += n;
    }   

    return aa;
}

AltArr_t* leadingTerm_AA_unpk (AltArr_t* aa, int nvar) {
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	if (!aa->unpacked) {
		return leadingTerm_AA(aa, nvar);
	}

	AltArr_t* lt = makePolynomial_AA_unpk(1, aa->nvar);
	mpq_init(lt->elems->coef);
	mpq_set(lt->elems->coef, aa->elems->coef);
	setExponentVector_unpk(lt->elems->degs, aa->elems->degs, aa->nvar);
	lt->size = 1;
	return lt;
}

int leadingVariable_AA_unpk (AltArr_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return -2;
	}

	if (!aa->unpacked) {
		return leadingVariable_AA(aa);
	}
	
	degree_t* degs = (degree_t*) aa->elems->degs;
	for (int i = 0; i < aa->nvar; ++i){
		if (degs[i] > 0) {
			return i;
		}
	}

	return -1;
}

int mainLeadingDegree_AA_unpk (AltArr_t* aa)
{
    if (aa == NULL || aa->size == 0){
		return 0;
    }

    if (!aa->unpacked) {
    	return mainLeadingDegree_AA(aa);
    }
    
    degree_t* degArr = (degree_t*) aa->elems[0].degs;
    return degArr[0];
}

AltArr_t* mainLeadingCoefficient_AA_unpk (AltArr_t* aa)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }

    if (!aa->unpacked) {
    	return mainLeadingCoefficient_AA(aa);
    }
    
    AAElem_t* elems = aa->elems;
    degree_t* degs = (degree_t*) elems[0].degs;
    degree_t mvarDeg = degs[0];
    degree_t curDeg = mvarDeg;
    
    int nvar = aa->nvar;
   	int size = aa->size;
    int polyElemsAlloc = 10;
    int polyElemsSize = 0;
    AAElem_t* polyElems = (AAElem_t*) malloc (sizeof(AAElem_t) * polyElemsAlloc);
    degree_t* polyDegs = (degree_t*) malloc(sizeof(degree_t) * polyElemsAlloc * nvar);
    
    for (int i = 0; (mvarDeg == curDeg) && i < size; ++i){
		if (polyElemsSize + 1 > polyElemsAlloc){
		    polyElemsAlloc += 10;
		    polyElems = (AAElem_t*) realloc (polyElems, sizeof (AAElem_t) * polyElemsAlloc);
		    polyDegs = (degree_t*) realloc(polyDegs, sizeof(degree_t) * polyElemsAlloc * nvar);
		}

		mpq_init (polyElems[i].coef);
		mpq_set (polyElems[i].coef, elems[i].coef);
		setExponentVector_unpk((degrees_t) (polyDegs + i*nvar), elems[i].degs, nvar);
		polyDegs[i*nvar] = 0;
		++polyElemsSize;
			
		if (i+1 < size){
		    curDeg = degs[(i+1)*nvar];
		} else {
		    curDeg = -1;
		}
    }

    AltArr_t* poly = (AltArr_t*) malloc(sizeof(AltArr_t));
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

AltArr_t* mainCoefficientAtIdx_AA_unpk (AltArr_t* aa, int e)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (e < 0){
		return NULL;
    }
    if (!aa->unpacked) {
    	return mainCoefficientAtIdx_AA(aa, e);
    }

    AAElem_t* elems = aa->elems;
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
    AAElem_t* polyElems = (AAElem_t*) malloc (sizeof(AAElem_t) * polyElemsAlloc);
	degree_t* polyDegs = (degree_t*) malloc(sizeof(degree_t) * nvar * polyElemsAlloc);    

    for (int i = 0; i < size; ++i){
		if (mvarDeg == curDeg){
		    if (polyElemsSize + 1 > polyElemsAlloc) {
				polyElemsAlloc += 10;
				polyElems = (AAElem_t*) realloc (polyElems, sizeof (AAElem_t)*polyElemsAlloc);
		    	polyDegs = (degree_t*) realloc(polyDegs, sizeof(degree_t) * polyElemsAlloc * nvar);
		    }
	    
		    mpq_init (polyElems[polyElemsSize].coef);
		    mpq_set (polyElems[polyElemsSize].coef, elems[i].coef);
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
    
    AltArr_t* poly = (AltArr_t*) malloc(sizeof(AltArr_t));
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

AltArr_t* maxPolynomials_AA_unpk (AltArr_t* a, AltArr_t* b) {
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AA(a);
	}

	if (a->nvar != b->nvar) {
		fprintf(stderr, "maxPolynomials_AA inputs are not from the same ring\n");
		exit(1);
	}

	if (!a->unpacked && !b->unpacked) {
		return maxPolynomials_AA(a,b);
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
		return deepCopyPolynomial_AA(a);
	} else if (cmp < 0){
		return deepCopyPolynomial_AA(b);
	} else {
		if (mpq_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return deepCopyPolynomial_AA(a);
		} else {
			return deepCopyPolynomial_AA(b);				
		}
	}

	return NULL;
}

AltArr_t* maxPolynomials_AA_inp_unpk (AltArr_t* a, AltArr_t* b) {
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0){
		return a;
	}

	if (!a->unpacked && !b->unpacked) {
		return maxPolynomials_AA(a,b);
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
		return deepCopyPolynomial_AA_unpk(b);
	} else {
		if (mpq_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return a;
		} else {
			return deepCopyPolynomial_AA_unpk(b);
		}
	}
	return NULL;
}

void addRationalNumber_AA_inp_unpk(AltArr_t* aa, const mpq_t coef) {
	if (aa == NULL) {
		return;
	}

	if (aa->size == 0) {
		if (aa->alloc < 1) {
			resizePolynomial_AA(aa, 1);
		}
		mpq_set(aa->elems->coef, coef);
		return;
	}

	if (!aa->unpacked) {
		addRationalNumber_AA_inp(aa, coef);
		return;
	}

	int nvar = aa->nvar;
    if (isZeroExponentVector_unpk(aa->elems[aa->size-1].degs, nvar)) {
        mpq_add(aa->elems[aa->size-1].coef, aa->elems[aa->size-1].coef, coef);
        return;
    } 

    if (aa->size >= aa->alloc) {
        resizePolynomial_AA_unpk(aa, aa->alloc+10);
    }

    mpq_init(aa->elems[aa->size].coef);
    mpq_set(aa->elems[aa->size].coef, coef);

    degree_t* aDegs = (degree_t*) aa->elems->degs;
    aa->elems[aa->size].degs = (degrees_t) (aDegs + aa->size*nvar);
	aDegs = aDegs + aa->size*nvar;
    for (int i = 0; i < nvar; ++i) {
    	aDegs[i] = 0;
    }
    ++(aa->size);
}

AltArr_t* addPolynomials_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AA_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return addPolynomials_AA(a, b, nvar);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b); 
	
	AltArr_t* c = makePolynomial_AA_unpk(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			mpq_init(cElems[k].coef);
			mpq_add(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
			if (mpq_sgn(cElems[k].coef) == 0) {
				mpq_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, aElems[i].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, aElems[i].coef);
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, bElems[j].coef);
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;

	resizePolynomial_AA_unpk(c, k);

	if (unpackedA) {
		packExponentVectors_AA_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(c);

	return c;

}

AltArr_t* subPolynomials_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AA_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return addPolynomials_AA(a, b, nvar);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b); 
	
	AltArr_t* c = makePolynomial_AA_unpk(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_neg(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			mpq_init(cElems[k].coef);
			mpq_sub(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
			if (mpq_sgn(cElems[k].coef) == 0) {
				mpq_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, aElems[i].coef);
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, aElems[i].coef);
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_neg(cElems[k].coef, bElems[j].coef);
		memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	if (unpackedA) {
		packExponentVectors_AA_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	AA_SIZE(c) = k;
	if (k == 0) {
		freePolynomial_AA(c);
		return NULL;
	}
	resizePolynomial_AA_unpk(c, k);

	tryPackExponentVectors_AA_inp(c);

	return c;

}

AltArr_t* addPolynomials_AA_inp_unpk(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AA_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return addPolynomials_AA_inp(a, b, nvar);
	}

	int unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b); 
	
	AltArr_t* c = makePolynomial_AA_unpk(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			// cElems[k] = aElems[i];
			cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
			cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			// mpq_init(cElems[k].coef);
			mpq_add(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpq_sgn(cElems[k].coef) == 0) {
				mpq_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			// cElems[k] = aElems[i];
			cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
			cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			++k;
			++i;
		}
	}

	while(i < asize) {
		cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
		cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// memcpy(cElems + k, aElems + i, sizeof(AAElem_t)*(asize - i));
		++k;
		++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, bElems[j].coef);
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
	AA_SIZE(c) = k;

	resizePolynomial_AA(c, k);

	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(c);

	return c;
}

AltArr_t* subPolynomials_AA_inp_unpk(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA_unpk(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AA_unpk(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return subPolynomials_AA_inp(a, b, nvar);
	}

	int unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b); 
	
	AltArr_t* c = makePolynomial_AA_unpk(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* bDegs = (degree_t*) b->elems[0].degs;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;
	size_t expVecSize = sizeof(degree_t)*nvar;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_neg(cElems[k].coef, bElems[j].coef);
			memcpy(cDegs + k*nvar, bDegs + j*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);
			// cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors_unpk(aElems[i].degs, bElems[j].degs, nvar)) {
			// a==b
			// cElems[k] = aElems[i];
			cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
			cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			// mpq_init(cElems[k].coef);
			mpq_sub(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpq_sgn(cElems[k].coef) == 0) {
				mpq_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			// cElems[k] = aElems[i];
			cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
			cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
			memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
			cElems[k].degs = (degrees_t) (cDegs + k*nvar);

			++k;
			++i;
		}
	}

	while(i < asize) {
		cElems[k].coef->_mp_num = aElems[i].coef->_mp_num;
		cElems[k].coef->_mp_den = aElems[i].coef->_mp_den;
		memcpy(cDegs + k*nvar, aDegs + i*nvar, expVecSize);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		// memcpy(cElems + k, aElems + i, sizeof(AAElem_t)*(asize - i));
		++k;
		++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_neg(cElems[k].coef, bElems[j].coef);
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
		packExponentVectors_AA_inp(b);
	}

	AA_SIZE(c) = k;
	if (k == 0) {
		freePolynomial_AA(c);
		return NULL;
	}
	resizePolynomial_AA(c, k);

	tryPackExponentVectors_AA_inp(c);

	return c;
}

ProductHeap_AA* prodheapInit_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar) {
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	h->elements = (ProductHeapElem_AA*) malloc(sizeof(ProductHeapElem_AA)*AA_SIZE(a));
	degree_t* heapDegs = (degree_t*) malloc(sizeof(degree_t)*AA_SIZE(a)*nvar);
	h->elements[0].chain = prodheapMakeChain_AA(0, 0, NULL);
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

static inline void prodheapFree_AA_unpk(ProductHeap_AA* h) {
	free(h->unpackedDegs);
	prodheapFree_AA(h);
}

static inline void prodheapSetElem_unpk(ProductHeapElem_AA* dest, ProductHeapElem_AA* src, int nvar) {
	dest->chain = src->chain;
	setExponentVector_unpk(dest->degs, src->degs, nvar);
}

void prodheapInsert_AA_unpk(ProductHeap_AA* h, ProductHeapChain_AA* chain, degrees_t degs) {
	register int s = h->heapSize;
	ProductHeapElem_AA* elems = h->elements;
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
	ProductHeapElem_AA temp = {(degrees_t) tempDegs, NULL};
	setExponentVector_unpk(elems[s].degs, degs, nvar);
	ProductHeapElem_AA elem = {elems[s].degs, chain};

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

ProductHeapChain_AA* prodheapRemoveMax_AA_unpk(ProductHeap_AA* h) {
	ProductHeapElem_AA* elems = h->elements;
	ProductHeapChain_AA* maxElem = elems[0].chain;
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

AltArr_t* multiplyPolynomials_AA_unpk(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL || a->size == 0 || b == NULL || b->size == 0) {
		return NULL;
	}

	if (!a->unpacked && !b->unpacked) {
		return multiplyPolynomials_AA(a,b,nvar);
	}

	// reorder to obtain smaller as a. 
	if (b->size < a->size) {
		AltArr_t* temp = a;
		a = b;
		b = temp;
	}
	
	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	ProductHeap_AA* h = prodheapInit_AA_unpk(a,b,nvar);
	ratNum_t ccoef;
	mpq_init(ccoef);

	register unsigned long long int allocC = AA_SIZE(a)*AA_SIZE(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;
	AltArr_t* c = makePolynomial_AA_unpk(allocC, nvar);

	//k is c, i is a, j is b.
	register int k = 0;
	// register int i = 0;
	// register int j = 0;

	AAElem_t* __restrict__ cElems = c->elems;
	AAElem_t* __restrict__ bElems = b->elems;
	degree_t* cDegs = (degree_t*) c->elems[0].degs;

	AAElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1; 
	register int firstB = 0; 	

	ProductHeapChain_AA* maxElem = NULL;
	ProductHeapChain_AA* nextMaxElem = NULL;
	degrees_t* nextDegs;
	degree_t tempDegs[nvar];
	degree_t* tempDegs_p = tempDegs;

	while ( (nextDegs = prodheapPeek_AA(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		setExponentVector_unpk((degrees_t) tempDegs_p, h->elements[0].degs, nvar);
		cElems[k].degs = (degrees_t) (cDegs + k*nvar);
		setExponentVector_unpk(cElems[k].degs, *nextDegs, nvar);
		// cElems[k].degs = *nextDegs;
		mpq_init(cElems[k].coef);

		while (nextDegs != NULL && isEqualExponentVectors_unpk(cElems[k].degs, *nextDegs, nvar)) {
			//we will extract and accumulate the coefficents 
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AA* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AA_unpk(h);
			while (maxElem != NULL) {
				mpq_mul(ccoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				mpq_add(cElems[k].coef, cElems[k].coef, ccoef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA) {
					oldMaxElem = prodheapMakeChain_AA((maxElem->a_i)+1, firstB, oldMaxElem);
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
					prodheapFreeChain_AA(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;	

			nextDegs = prodheapPeek_AA(h);		
		}

		//Commit new term to the product.
		if (mpq_sgn(cElems[k].coef) != 0) {
			++k;
		} else {
			//reset accumulator variables and do not increment k.
			//will init cElem[k] again on next loop, so clear here.
			mpq_clear(cElems[k].coef);
		}
		
		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			addExponentVectors_unpk(aElems[maxElem->a_i].degs, bElems[maxElem->b].degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AA_unpk(h, maxElem, (degrees_t) tempDegs_p);
			maxElem = nextMaxElem;
		}
	
		if (k >= allocC) {
			allocC <<= 1;
			resizePolynomial_AA_unpk(c, allocC);
			// fprintf(stderr, "\n\nRESIZING\n\n");
			cElems = c->elems;
		}
	}

	mpq_clear(ccoef);
	prodheapFree_AA_unpk(h);

	AA_SIZE(c) = k;

	resizePolynomial_AA_unpk(c, k);

	if (unpackedA) {
		packExponentVectors_AA_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	return c;
}

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
void divisionGetNextTerm_AA_unpk(ProductHeap_AA* h, const AAElem_t* __restrict__ aElems, const AAElem_t* __restrict__ bElems, mpq_t* retCoef) {
	if (h->heapSize == 0) {
		return;
	}

	int lastB = h->lastB;
	int nvar = h->nvar;

	ProductHeapChain_AA* insertChain = NULL;
	ProductHeapChain_AA* maxElem, *nextMaxElem;

	mpq_t prodCoef;
	mpq_init(prodCoef);
	degrees_t* nextDegs = prodheapPeek_AA(h);

	degree_t maxDegsUnpk[nvar];
	degrees_t maxDegs = (degrees_t) maxDegsUnpk;
	setExponentVector_unpk(maxDegs, *nextDegs, nvar);
	// degrees_t maxDegs = *nextDegs;

	while ( nextDegs != NULL && isEqualExponentVectors_unpk(maxDegs, *nextDegs, nvar)) {
		maxElem = prodheapRemoveMax_AA_unpk(h);
				
		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;
			mpq_mul(prodCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
			mpq_add(*retCoef, *retCoef, prodCoef);
			if (maxElem->b != lastB) {
				++(maxElem->b);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				maxElem->next = NULL;
				prodheapFreeChain_AA(maxElem);
			}

			maxElem = nextMaxElem;
		}

		nextDegs = prodheapPeek_AA(h);
	}

	degree_t tempDegs[nvar];
	degree_t* tempDegs_p = tempDegs;
	while(insertChain != NULL) {
		maxElem = insertChain->next;
		insertChain->next = NULL;
		addExponentVectors_unpk(aElems[insertChain->a_i].degs, bElems[insertChain->b].degs, (degrees_t) tempDegs_p, nvar);
		prodheapInsert_AA_unpk(h, insertChain, (degrees_t) tempDegs_p);
		insertChain = maxElem;
	}

	mpq_clear(prodCoef);
}

void divideBySingleTerm_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int nvar) {
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
		divideBySingleTerm_AA(c, b, res_a, res_r, nvar);
		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AA_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	AAElem_t* k = c->elems;
	AAElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0; 

	AltArr_t* a = makePolynomial_AA_unpk(maxSize, nvar);
	AltArr_t* r = makePolynomial_AA_unpk(maxSize, nvar);

	AAElem_t* bElems = b->elems;
	AAElem_t* curA = a->elems;
	AAElem_t* curR = r->elems;
	degree_t* aDegs = (degree_t*) a->elems[0].degs;
	degree_t* rDegs = (degree_t*) r->elems[0].degs;

	mpq_init(curA->coef);
	mpq_init(curR->coef);

	while (k != lenK) {
		if (monomialDivideTest_unpk(k->degs,bElems->degs,nvar)) {
			mpq_div(curA->coef, k->coef, bElems->coef);
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(k->degs, bElems->degs, curA->degs, nvar);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {
			mpq_set(curR->coef, k->coef);
			curR->degs = (degrees_t) (rDegs + j*nvar);
			setExponentVector_unpk(curR->degs, k->degs, nvar);
			// curR->degs = k->degs;
			++j;
			++(curR);
			mpq_init(curR->coef);
		}
		++k;
	}

	mpq_clear(curA->coef);
	mpq_clear(curR->coef);

	AA_SIZE(a) = i; 
	AA_SIZE(r) = j; 

	if (unpackedC) {
		packExponentVectors_AA_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(a);
	tryPackExponentVectors_AA_inp(r);

	*res_a = a;
	*res_r = r;
}

 void prodheapResize_AA_unpk(ProductHeap_AA* h, int newAllocSize) {
	h->elements = (ProductHeapElem_AA*) realloc(h->elements, sizeof(ProductHeapElem_AA)*newAllocSize);
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
void dividePolynomials_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, register int nvar) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AA(0, nvar);
		*res_r = NULL;
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		divideBySingleTerm_AA(c, b, res_a, res_r, nvar);
		return;
	}

	if (!c->unpacked && !b->unpacked) {
		dividePolynomials_AA(c, b, res_a, res_r, nvar);
		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AA_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}


	AAElem_t* __restrict__ k = c->elems;
	AAElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0; 

	AltArr_t* a = makePolynomial_AA_unpk(maxASize, nvar);
	AltArr_t* r = makePolynomial_AA_unpk(maxRSize, nvar);
	AAElem_t* __restrict__ curA = a->elems;
	AAElem_t* __restrict__ curR = r->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems[0].degs;
	degree_t* __restrict__ rDegs = (degree_t*) r->elems[0].degs;
	degree_t tempDegs[nvar];
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpq_init(curA->coef);
	mpq_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	while (k != lenK && !monomialDivideTest_unpk(k->degs, beta, nvar)) {
		mpq_set(curR->coef, k->coef);
		curR->degs = (degrees_t) (rDegs + j*nvar);
		setExponentVector_unpk(curR->degs, k->degs, nvar);
		// curR->degs = k->degs;
		++j;
		if (j >= maxRSize) {
			maxRSize <<= 1;
			r->elems = (AAElem_t*) realloc(r->elems, maxRSize*sizeof(AAElem_t));
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
		mpq_init(curR->coef);
		++k;
	}

	if (k == lenK) {
		//no division to do at all!
		mpq_clear(curA->coef);
		mpq_clear(curR->coef);

		AA_SIZE(a) = i;
		AA_SIZE(r) = j;
		a->alloc = maxASize;
		r->alloc = maxRSize;

		*res_a = a;
		*res_r = r;
		return; 
	}

	curA->degs = (degrees_t) (aDegs + i*nvar);
	subtractExponentVectors_unpk(k->degs, beta, curA->degs, nvar);
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	prodheapResize_AA_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AA_unpk(h, prodheapMakeChain_AA(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	mpq_init(curA->coef);

	degrees_t* delta = prodheapPeek_AA(h);
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
			divisionGetNextTerm_AA_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpq_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AA(h);
				continue;
			} else {
				mpq_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AA_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpq_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AA(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpq_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpq_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AA(h);
					continue;
				}
			}
		} else {
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpq_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest_unpk(eps, beta, nvar)) {
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAElem_t*) realloc(a->elems, maxASize*sizeof(AAElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = (degree_t*) realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AA_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AA_unpk(h, prodheapMakeChain_AA(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {

			//swap here so that curA becomes 0.
			mpq_swap(curR->coef, curA->coef);
			curR->degs = (degrees_t) (rDegs + j*nvar);
			setExponentVector_unpk(curR->degs, eps, nvar);
			// curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAElem_t*) realloc(r->elems, maxRSize*sizeof(AAElem_t));
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
			mpq_init(curR->coef);
		}

		delta = prodheapPeek_AA(h);
	}

	prodheapFree_AA_unpk(h);

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);
	mpq_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	if (unpackedC) {
		packExponentVectors_AA_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(a);
	tryPackExponentVectors_AA_inp(r);

	*res_a = a;
	*res_r = r;
} 

void exactDividePolynomials_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, register int nvar) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = NULL;
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		AltArr_t* res_r = NULL;
		divideBySingleTerm_AA(c, b, res_a, &res_r, nvar);
		freePolynomial_AA(res_r);
		return;
	}

	if (!c->unpacked && !b->unpacked) {
		exactDividePolynomials_AA(c, b, res_a, nvar);
		return;
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AA_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	AAElem_t* __restrict__ k = c->elems;
	AAElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	// register int maxRSize = maxASize;
	register int i = 0;
	// register int j = 0; 

	AltArr_t* a = makePolynomial_AA_unpk(maxASize, nvar);
	// AltArr_t* r = makePolynomial_AA_unpk(maxRSize, nvar);
	AAElem_t* __restrict__ curA = a->elems;
	// AAElem_t* __restrict__ curR = r->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems[0].degs;
	// degree_t* __restrict__ rDegs = (degree_t*) r->elems[0].degs;
	degree_t tempDegs[nvar];
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpq_init(curA->coef);
	// mpq_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	// while (k != lenK && !monomialDivideTest_unpk(k->degs, beta, nvar)) {
	// 	mpq_set(curR->coef, k->coef);
	// 	curR->degs = (degrees_t) (rDegs + j*nvar);
	// 	setExponentVector_unpk(curR->degs, k->degs, nvar);
	// 	// curR->degs = k->degs;
	// 	++j;
	// 	if (j >= maxRSize) {
	// 		maxRSize <<= 1;
	// 		r->elems = (AAElem_t*) realloc(r->elems, maxRSize*sizeof(AAElem_t));
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
	// 	mpq_init(curR->coef);
	// 	++k;
	// }
	// if (k == lenK) {
	// 	//no division to do at all!
	// 	mpq_clear(curA->coef);
	// 	mpq_clear(curR->coef);

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
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	prodheapResize_AA_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AA_unpk(h, prodheapMakeChain_AA(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	mpq_init(curA->coef);

	degrees_t* delta = prodheapPeek_AA(h);
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
			divisionGetNextTerm_AA_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpq_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AA(h);
				continue;
			} else {
				mpq_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AA_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpq_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AA(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpq_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpq_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AA(h);
					continue;
				}
			}
		} else {
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpq_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest_unpk(eps, beta, nvar)) {
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAElem_t*) realloc(a->elems, maxASize*sizeof(AAElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AA_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AA_unpk(h, prodheapMakeChain_AA(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} 
		// else {

		// 	//swap here so that curA becomes 0.
		// 	mpq_swap(curR->coef, curA->coef);
		// 	curR->degs = (degrees_t) (rDegs + j*nvar);
		// 	setExponentVector_unpk(curR->degs, eps, nvar);
		// 	// curR->degs = eps;
		// 	++j;
		// 	if (j >= maxRSize) {
		// 		maxRSize <<= 1;
		// 		r->elems = (AAElem_t*) realloc(r->elems, maxRSize*sizeof(AAElem_t));
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
		// 	mpq_init(curR->coef);
		// }

		delta = prodheapPeek_AA(h);
	}

	prodheapFree_AA_unpk(h);

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);
	// mpq_clear(curR->coef);

	AA_SIZE(a) = i;
	// AA_SIZE(r) = j;
	a->alloc = maxASize;
	// r->alloc = maxRSize;

	if (unpackedC) {
		packExponentVectors_AA_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(a);
	// tryPackExponentVectors_AA_inp(r);

	*res_a = a;
	// *res_r = r;
} 

int divideTestSingleTerm_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, int nvar) {
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
		return divideTestSingleTerm_AA(c, b, res_a, nvar);
		
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AA_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	AAElem_t* k = c->elems;
	AAElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;

	AltArr_t* a = makePolynomial_AA_unpk(maxSize, nvar);

	AAElem_t* bElems = b->elems;
	AAElem_t* curA = a->elems;
	degree_t* aDegs = (degree_t*) a->elems[0].degs;

	mpq_init(curA->coef);

	int ret = 1;
	while (k != lenK) {
		if (monomialDivideTest_unpk(k->degs,bElems->degs,nvar)) {
			mpq_div(curA->coef, k->coef, bElems->coef);
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(k->degs, bElems->degs, curA->degs, nvar);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {
			mpq_clear(curA->coef);
			a->size = i;
			freePolynomial_AA_unpk(a);
			a = NULL;
			ret = 0;
			break;
		}
		++k;
	}

	if (a != NULL) {
		mpq_clear(curA->coef);
	}

	if (ret) {
		AA_SIZE(a) = i; 
		tryPackExponentVectors_AA_inp(a);
		*res_a = a;
	}

	if (unpackedC) {
		packExponentVectors_AA_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	return ret;
}


int divideTest_AA_unpk(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, int nvar) {

	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AA(0, nvar);
		return 1;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		return divideTestSingleTerm_AA_unpk(c, b, res_a, nvar);
	}

	if (!c->unpacked && !b->unpacked) {
		return divideTest_AA(c, b, res_a, nvar);
	}

	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AA_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}


	AAElem_t* __restrict__ k = c->elems;
	AAElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int i = 0;

	AltArr_t* a = makePolynomial_AA_unpk(maxASize, nvar);
	AAElem_t* __restrict__ curA = a->elems;
	degree_t* __restrict__ aDegs = (degree_t*) a->elems[0].degs;
	degree_t tempDegs[nvar];
	degree_t* __restrict__ tempDegs_p = tempDegs;

	mpq_init(curA->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	if (k != lenK && !monomialDivideTest_unpk(k->degs, beta, nvar)) {
		freePolynomial_AA_unpk(a);
		if (unpackedC) {
			packExponentVectors_AA_inp(c);
		}
		if (unpackedB) {
			packExponentVectors_AA_inp(b);
		}
		return 0;
	}

	curA->degs = (degrees_t) (aDegs + i*nvar);
	subtractExponentVectors_unpk(k->degs, beta, curA->degs, nvar);
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	prodheapResize_AA_unpk(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
	prodheapInsert_AA_unpk(h, prodheapMakeChain_AA(0, 1, NULL), (degrees_t) tempDegs_p);
	++i;
	++curA;
	mpq_init(curA->coef);

	degrees_t* delta = prodheapPeek_AA(h);
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
			divisionGetNextTerm_AA_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpq_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AA(h);
				continue;
			} else {
				mpq_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			setExponentVector_unpk(eps, *delta, nvar);
			// eps = *delta;
			divisionGetNextTerm_AA_unpk(h, a->elems, b->elems, &(curA->coef));
			if (mpq_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AA(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpq_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpq_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AA(h);
					continue;
				}
			}
		} else {
			setExponentVector_unpk(eps, k->degs, nvar);
			// eps = k->degs;
			mpq_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest_unpk(eps, beta, nvar)) {
			curA->degs = (degrees_t) (aDegs + i*nvar);
			subtractExponentVectors_unpk(eps, beta, curA->degs, nvar);
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAElem_t*) realloc(a->elems, maxASize*sizeof(AAElem_t));
				degree_t* oldDegs = aDegs;
				aDegs = realloc(aDegs, sizeof(degree_t)*nvar*maxASize);
				if (oldDegs != aDegs) {
					for (int idx = 0; idx <= i; ++idx) {
						a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
					}
				}
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AA_unpk(h, maxASize);
			}

			addExponentVectors_unpk(curA->degs, b2Elem->degs, (degrees_t) tempDegs_p, nvar);
			prodheapInsert_AA_unpk(h, prodheapMakeChain_AA(i, 1, NULL), (degrees_t) tempDegs_p);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {
			a->size = i;
			freePolynomial_AA_unpk(a);
			if (unpackedC) {
				packExponentVectors_AA_inp(c);
			}
			if (unpackedB) {
				packExponentVectors_AA_inp(b);
			}
			return 0;
		}

		delta = prodheapPeek_AA(h);
	}

	prodheapFree_AA_unpk(h);

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxASize;

	if (unpackedC) {
		packExponentVectors_AA_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(a);

	*res_a = a;
	return 1;
}

void divideByLeadingTerms_AA_unpk (AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int nvar)
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
		divideByLeadingTerms_AA (c, b, res_a, res_r, nvar);
 		return;
	}
	
	int unpackedC = 0, unpackedB = 0;
	if (!c->unpacked) {
		unpackExponentVectors_AA_inp(c);
		unpackedC = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}
	
	AltArr_t* a;
	AltArr_t* r;
	
	/* while (k != lenK) { */
	if (monomialDivideTest_AA_unpk (c, 0, b, 0)) {
		// a:
		a = makePolynomial_AA_unpk (1, nvar);
		mpq_init (a->elems->coef);
		mpq_div (a->elems->coef, c->elems->coef, b->elems->coef);
		subtractExponentVectors_unpk (c->elems->degs, b->elems->degs, a->elems->degs, nvar);
		AA_SIZE (a) = 1;
		
		// r:
		r = NULL;
		
	} else {
		// a:
		a = makePolynomial_AA_unpk (1, nvar);
		mpq_init (a->elems->coef);
		mpq_set_si (a->elems->coef, 1l, 1lu);
		a->elems->degs = (degrees_t) 0;
		AA_SIZE (a) = 1;

		// r:
		r = makePolynomial_AA_unpk (1, nvar);
		mpq_init (r->elems->coef);
		mpq_set (r->elems->coef, c->elems->coef);
		/* r->elems->degs = c->elems->degs; */
		setExponentVector_unpk (r->elems->degs, c->elems->degs, nvar);
		AA_SIZE (r) = 1;
	}

	if (unpackedC) {
		packExponentVectors_AA_inp(c);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	tryPackExponentVectors_AA_inp(a);
	tryPackExponentVectors_AA_inp(r);
	
	*res_a = a;
	*res_r = r;
}

 
/*****************
 * Derivative / Integral
 *****************/

AltArr_t* derivative_AA_unpk(AltArr_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx < 0 || idx >= aa->nvar) {
		return NULL;
	}

	if (!aa->unpacked) {
		return derivative_AA(aa, idx, k);
	}

	int nvar = aa->nvar;
	AltArr_t* ret = makePolynomial_AA_unpk(aa->size, aa->nvar);
	
	int insertIdx = 0;
	int size = aa->size;
	AAElem_t* __restrict__ elems = aa->elems;
	degree_t* degs_p = (degree_t*) elems[0].degs;
	AAElem_t* __restrict__ retElems = ret->elems;
	degree_t* retDegs_p = (degree_t*) retElems[0].degs;
	degrees_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
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
		mpq_init(retElems[insertIdx].coef);
		mpq_set(retElems[insertIdx].coef, elems[i].coef);

		mpq_set_ui(mpqDeg, deg, 1ul);
		for (int j = 0; j < k; ++j) {
			mpq_mul(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);			
			--deg;	
			mpz_sub(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}
		retDegs_p[insertIdx*nvar + idx] = deg;
		// retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	tryPackExponentVectors_AA_inp(ret);

	return ret;
}

/**
 * Integrate with respect to a variable that does not currently exist in aa.
 * It becomes the main variable (that is, to the left).
 *
 */
AltArr_t* integrateExpand_AA_unpk(AltArr_t* aa, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return integrateExpand_AA(aa, k);
	}

	AltArr_t* ret = deepCopyPolynomial_AA(aa);

	expandNumVarsLeft_AA_unpk(ret, (aa->nvar)+1);
	
	integrateExpand_AA_inp_unpk(ret, k);

	return ret;
}

/**
 * Note, it is assumed aa has already been expanded.
 */
void integrateExpand_AA_inp_unpk(AltArr_t* aa, int k) {
	if (!aa->unpacked) {
		unpackExponentVectors_AA_inp(aa);
	}

	int nvar = aa->nvar;
	int size = aa->size;
	AAElem_t* elems = aa->elems;
	degree_t* degs = (degree_t*) aa->elems[0].degs;

	mpq_t kFact;
	mpq_init(kFact);
	mpz_fac_ui(mpq_numref(kFact), (unsigned long)k);

	for (int i = 0; i < size; ++i) {
		mpq_div(elems[i].coef, elems[i].coef, kFact);
		degs[i*nvar] = k;
	}

	mpq_clear(kFact);
}

AltArr_t* integral_AA_unpk(AltArr_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return integral_AA(aa, idx, k);
	}

	if (idx >= aa->nvar || idx < 0) {
		return integrateExpand_AA_unpk(aa, k);
	}

	AltArr_t* ret = makePolynomial_AA_unpk(aa->size, aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	int nvar = aa->nvar;
	AAElem_t* __restrict__ elems = aa->elems;
	AAElem_t* __restrict__ retElems = ret->elems;

	degree_t* __restrict__ degs = (degree_t*) aa->elems[0].degs;
	degree_t* __restrict__ retDegs = (degree_t*) ret->elems[0].degs;

	degrees_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i) {
		deg = degs[i*nvar + idx];
		
		retDegs[insertIdx*nvar + idx] = 0;
		mpq_init(retElems[insertIdx].coef);
		mpq_set(retElems[insertIdx].coef, elems[i].coef);

		mpq_set_ui(mpqDeg, deg+1, 1ul);
		for (int j = 0; j < k; ++j) {
			mpq_div(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);			
			++deg;
			mpz_add(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}

		memcpy(retDegs + insertIdx*nvar, degs + i*nvar, sizeof(degree_t)*nvar);
		retDegs[insertIdx*nvar + idx] = deg;
		retElems[insertIdx].degs = (degrees_t) (retDegs + insertIdx*nvar);
		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	return ret;
}

AltArr_t* univariateGCD_AA_unpk(AltArr_t* a, AltArr_t* b) {
	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMQP ERROR: Calling univariate GCD on multivariate arrays\n");
		exit(1);
	}

	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AA(a);
	}

	if (!a->unpacked && !b->unpacked) {
		return univariateGCD_AA(a, b);
	}

	int unpackedA = 0, unpackedB = 0;
	if (!a->unpacked) {
		unpackExponentVectors_AA_inp(a);
		unpackedA = 1;
	}
	if (!b->unpacked) {
		unpackExponentVectors_AA_inp(b);
		unpackedB = 1;
	}

	AltArr_t* r0 = NULL;
	AltArr_t* r1 = NULL;
	AltArr_t* r2 = NULL;

	mpq_t c0;
	mpq_t c1;
	mpq_init(c0);
	mpq_init(c1);

	if (isGreaterExponentVectors_unpk(a->elems->degs, b->elems->degs, 1)) {
		r0 = primitivePartAndContent_AA(a, c0);
		r1 = primitivePartAndContent_AA(b, c1);
		// r0 = deepCopyPolynomial_AA(a);
		// r1 = deepCopyPolynomial_AA(b);
	} else {
		r0 = primitivePartAndContent_AA(b, c0);
		r1 = primitivePartAndContent_AA(a, c1);
		// r0 = deepCopyPolynomial_AA(b);
		// r1 = deepCopyPolynomial_AA(a);
	}

	AltArr_t* quo = NULL;
	while (r1 != NULL && r1->size > 0) {
		dividePolynomials_AA(r0, r1, &quo, &r2, 1);
		freePolynomial_AA(quo);
		quo = NULL;

		freePolynomial_AA(r0);
		r0 = r1;
		r1 = r2;
		primitivePart_AA_inp(r1);
		r2 = NULL;
	}

	freePolynomial_AA(r1);
	freePolynomial_AA(r2);
	if (r0 != NULL && r0->size > 0 && mpq_sgn(r0->elems->coef) < 0) {
		negatePolynomial_AA(r0);
	}

	// gcd if we want to include gcd of integer coefficients
	// if (mpz_cmp_si(mpq_denref(c0), 1l) == 0 && mpz_cmp_si(mpq_denref(c1), 1l) == 0) {
	// 	//contents are both integers
	// 	mpq_t gcd;
	// 	mpq_init(gcd);
	// 	mpz_gcd(mpq_numref(gcd), mpq_numref(c0), mpq_numref(c1));
	// 	if (mpz_cmp_si(mpq_numref(gcd), 1l) != 0) {
	// 		for (int i = 0; i < r0->size; ++i) {
	// 			mpq_mul(r0->elems[i].coef, r0->elems[i].coef, gcd);
	// 		}
	// 	}
	// 	mpq_clear(gcd);
	// }

	mpq_clear(c0);
	mpq_clear(c1);

	if (unpackedA) {
		packExponentVectors_AA_inp(a);
	}
	if (unpackedB) {
		packExponentVectors_AA_inp(b);
	}

	return r0;

}


AltArr_t* commonFactor_AA_unpk(AltArr_t* a, AltArr_t** factored) {
	if (a == NULL || a->size == 0) {
		return NULL;
	}

	if (!a->unpacked) {
		return commonFactor_AA(a, factored);
	}

	AltArr_t* ret = makePolynomial_AA_unpk(1, a->nvar);
	mpq_init(ret->elems->coef);
	mpq_set_ui(ret->elems->coef, 1ul, 1ul);

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
		AltArr_t* factRet = deepCopyPolynomial_AA_unpk(a);
		if (! isZeroExponentVector_unpk(ret->elems->degs, nvar)) {
			for (int i = 0; i < a->size; ++i) {
				subtractExponentVectors_unpk(factRet->elems[i].degs, ret->elems->degs, factRet->elems[i].degs, nvar);
			}
		}
		*factored = factRet;
	}

	return ret;
}

void univarEvaluate_AA_unpk(AltArr_t* aa, const mpq_t point, mpq_t res) {
	
	if (aa == NULL || aa->size == 0) {
		mpq_set_si(res, 0l, 1l);
		return;
	}

	if (aa->nvar != 1) {
		fprintf(stderr, "Poly is not univariate in univarEvaluate_AA_unpk\n");
		exit(1);
	}

	if (!aa->unpacked) {
		univarEvaluate_AA(aa, point, res);
		return;
	}

	AAElem_t* elems = aa->elems;
	register int size = aa->size;
	mpq_set(res, elems->coef);
	degree_t* degs = (degree_t*) elems->degs;
	degree_t prevDeg = degs[0];
	degree_t nextDeg;
	for (int i = 1; i < size; ++i) {
		nextDeg = degs[i];
		for (degree_t j = prevDeg; j > nextDeg; --j) {
			mpq_mul(res, res, point);
		}
		mpq_add(res, res, elems[i].coef);
		prevDeg = nextDeg;
	}
	for (degrees_t j = prevDeg; j > 0; --j) {
		mpq_mul(res, res, point);
	}
}

