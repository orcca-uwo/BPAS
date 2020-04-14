
#include "RationalNumberPolynomial/SMQP_Support-AA.h"
#include "IntegerPolynomial/DUZP_Support.h"


/*****************
 * Exponent Vector functions.
 *****************/

#if SMQP_SUPPORT_DEBUG
/**
 * Print the contents of an exponent vector in a comma-separated list
 * d: the exponent vector
 * nvar: the number of elements in the vector
 */
void printDegs(degrees_t d, int nvar) {
	fprintf(stderr, "%llx", d);
	// fprintf(stderr, "%llu, %llu, %llu", (d & FIRST_EXP) >> FIRST_EXP_OFFSET,
										// (d & SECOND_EXP) >> SECOND_EXP_OFFSET,
										// (d & THIRD_EXP) >> THIRD_EXP_OFFSET);
}
#endif




/*****************
 * Node defintion and helper functions.
 *****************/

/**
 * Free a node, and all it's sucessors.
 */
void freePolynomial(Node* node) {
	Node* nextNode;
	while (node != NULL) {
		nextNode = node->next;
		freeNode(node);
		node = nextNode;
	}
}

/**
 * Deep copy the polynomail whose head in the input Node*.
 */
Node* deepCopyPolynomial(Node* node, int nvar) {
	Node* head = NULL;
	Node* tail = NULL;

	// degrees_t degs;
	// degrees_t nodeDegs;
	while (node != NULL) {
		// nodeDegs = node->degs;
		// degs = (degrees_t) malloc(sizeof(degree_t) * nvar);
		// memcpy(degs, nodeDegs, sizeof(degree_t) * nvar);
		tail = addTerm(tail, node->degs, node->coef);
		if (head == NULL) {
			head = tail;
		}
		node = node->next;
	}

	return head;
}

/**
 * Determine the number of terms in a polynomial given the head node
 * of the linked-list representation of the polynomial.
 * returns the number of terms including and following the supplied node
 */
 polysize_t numberOfTermsNode(Node* node) {
	polysize_t val = 0;
	Node *h = node;
	while(h != NULL){
		h = h->next;
		val++;
	}
	return val;
}


/**
 * Given a linked-list of nodes, whose head is poly,
 * sort the list based on the ordering of compareExponentVectors.
 * Updates the polynomial in place by modifying *poly to be the new head node.
 * Also returns a pointer which is the new head node.
 */
Node* sortPolynomial(Node** poly, int nvar) {
	//TODO not insertion sort.
	Node* sorted = NULL;

	Node* node = *poly;
	Node* next;
	while(node != NULL) {
		next = node->next;

		if (sorted == NULL || compareExponentVectors(node->degs, sorted->degs) > 0) {
			//if inserted node is greated than current head
			node->next = sorted;
			sorted = node;
		} else {
			//else traverse and find insertion point
			Node* cur = sorted;
			while (cur->next != NULL && compareExponentVectors(cur->next->degs, node->degs) > 0) {
				cur = cur->next;
			}
			node->next = cur->next;
			cur->next = node;
		}

		node = next;
	}

	condensePolyomial(sorted, nvar);

	*poly = sorted;
	return sorted;
}

/**
 * Given a polynomial in sorted order but with possible duplicates,
 * condense the list such that like-terms are combined into a single node.
 */
void condensePolyomial(Node* poly, int nvar) {
	Node* prev = poly;
	if (prev == NULL ) {
		return;
	}

	Node* node = prev->next;
	while (node != NULL) {
		if (compareExponentVectors(prev->degs, node->degs) == 0) {
			mpq_add(prev->coef, prev->coef, node->coef);
			prev->next = node->next;
			node->next = NULL;
			freeNode(node);
			node = prev->next;
		} else {
			prev = node;
			node = node->next;
		}
	}
}


/**
 * Add a term to the end of the polynomial linked list given the tail node,
 * trailingTerm, the exponent vector for the new term, d, and the coefficient,
 * ratNum_t. A copy of the input coef is made and stored in the node. But the
 * degrees_t is stored directly.
 *
 * Note that the input coef can be null. This signifies to create a Node with a
 * 0 coef.
 *
 * returns a pointer to the new tail of the list (the node created)
 */
Node* addTerm(Node* trailingTerm, degrees_t d, const ratNum_t coef){
	Node* temp = (Node*) malloc(sizeof(Node));

	temp->degs = d;
	mpq_init(temp->coef);
	if (coef != NULL) {
		mpq_set(temp->coef, coef);
	}
	// temp->coef = coef;
	temp->next = NULL;

	if (trailingTerm != NULL) {
		trailingTerm->next = temp;
	}
	return temp;
}

/**
 * Add a term to the end of linked list with tail trailingTerm. The new term
 * will have degs d and a coef of 0.
 *
 * returns a pointer to the new tail of the list (the node created)
 */
Node* addZeroTerm(Node* trailingTerm, degrees_t d) {
	mpq_t mpqZero;
	mpq_init(mpqZero);
	Node* ret = addTerm(trailingTerm, d, mpqZero);
	mpq_clear(mpqZero);
	return ret;
}

/**
 * Given two terms, as nodes a and b, multiply these two and return the
 * single result as a node.
 */
Node* multiplyTerms(Node* a, Node* b, int nvar) {
	if (a == NULL || b == NULL) {
		return NULL;
	}

	ratNum_t coef;
	mpq_init(coef);
	mpq_mul(coef, a->coef, b->coef);
	// ratNum_t coef = a->coef * b->coef;
	// degrees_t degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
	degrees_t degs;
	addExponentVectors(a->degs, b->degs, degs);
	Node* c = addTerm(NULL, degs, coef);
	mpq_clear(coef);
	return c;
}


/**
 * Negate a polynomial.
 * Given the head of a polynomial, negate all coefficients in place.
 */
void negatePolynomial(Node* a) {
	Node* node = a;
	while (node != NULL) {
		mpq_neg(node->coef, node->coef);
		node = node->next;
	}
}

Node* leadingTerm (Node* a){
  Node* n = addTerm(NULL, a->degs, a->coef);
  return n;
}

/*****************
 * Alternating Array helpers
 *****************/

void printAA(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		fprintf(stderr, "0\n");
		return;
	}
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		gmp_fprintf(stderr, "%Qd*%llx + ", aa->elems[i].coef, aa->elems[i].degs);
	}
	fprintf(stderr, "\n");
}

void printDegs_AA(FILE* fp, degrees_t degs, const char** syms, int nvar, const degrees_t* masks, const int* sizes) {
	if (degs == 0) {
		return;
	}
	degree_t deg = GET_NTH_EXP(degs, masks[0], sizes[0]);
	fprintf(fp, "%s^%d", syms[0], deg);
	for (int k = 1; k < nvar; ++k) {
		deg = GET_NTH_EXP(degs, masks[k], sizes[k]);
		fprintf(fp, "*%s^%d", syms[k], deg);
	}
}

void printPositiveDegs_AA(FILE* fp, degrees_t degs, const char** syms, int nvar, const degrees_t* masks, const int* sizes) {
	if (degs == 0) {
		return;
	}
	degree_t deg = GET_NTH_EXP(degs, masks[0], sizes[0]);
	int first = 1;
	for (int k = 0; k < nvar; ++k) {
		deg = GET_NTH_EXP(degs, masks[k], sizes[k]);
		if (deg > 0) {
			if (first) {
				fprintf(fp, "%s^%d", syms[k], deg);
				first = 0;
			} else {
				fprintf(fp, "*%s^%d", syms[k], deg);
			}
		}
	}
}

void printPolyOptions_AA(FILE* fp, const AltArr_t* aa, const char** syms, int nvar, int positiveDegsOnly) {
	if (aa == NULL || aa->size == 0) {
		fprintf(stderr, "0\n");
		return;
	}
	if (aa->unpacked) {
		printPoly_AA_unpk(fp, aa, syms, nvar);
		return;
	}

	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	gmp_fprintf(fp, "%Qd", aa->elems[0].coef);
	if (!isZeroExponentVector(aa->elems[0].degs)) {
		fprintf(fp, "*");
		if (positiveDegsOnly) {
			printPositiveDegs_AA(fp, aa->elems[0].degs, syms, nvar, masks, sizes);
		} else {
			printDegs_AA(fp, aa->elems[0].degs, syms, nvar, masks, sizes);
		}
	}
	for (int i = 1; i < AA_SIZE(aa)-1; ++i) {
		if (mpq_sgn(aa->elems[i].coef) > 0) {
			fprintf(fp, " + ");
		} else {
			fprintf(fp, " ");
		}
		gmp_fprintf(fp, "%Qd", aa->elems[i].coef);
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
		if (mpq_sgn(aa->elems[aa->size-1].coef) > 0) {
			fprintf(fp, " + ");
		} else {
			fprintf(fp, " ");
		}
		gmp_fprintf(fp, "%Qd", aa->elems[aa->size-1].coef);
		if (!isZeroExponentVector(aa->elems[aa->size-1].degs)) {
			fprintf(fp, "*");
			if (positiveDegsOnly) {
				printDegs_AA(fp, aa->elems[aa->size-1].degs, syms, nvar, masks, sizes);
			} else {
				printDegs_AA(fp, aa->elems[aa->size-1].degs, syms, nvar, masks, sizes);
			}
		}
	}

	//fprintf(fp, "\n");
}

degrees_t calculateMaxDegs_AA(const AltArr_t* aa) {
	if (aa == NULL || aa->size == 0 || aa->nvar == 0) {
		return 0;
	}

	AAElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	degrees_t* __restrict__ maxList = (degrees_t*) calloc(nvar, sizeof(degrees_t));
	degrees_t max = 0ll;
	for (register int i = 0; i < size; ++i) {
		for (register int j = 0; j < nvar; ++j) {
			maxList[j] = (elems[i].degs & masks[j]) > (maxList[j]) ? (elems[i].degs & masks[j]) : maxList[j];

			// if ( (elems[i].degs & masks[j]) > (maxList[j]) ) {
				// maxList[j] = (elems[i].degs & masks[j]);
			// }
			// if ( (elems[i].degs & masks[j]) > (max & masks[j]) ) {
			// 	max = max | (elems[i].degs & masks[j]);
			// }
		}
	}

	for (register int j = 0; j < nvar; ++j) {
		max |= maxList[j];
	}

	return max;
}

int isInOrder_AA(AltArr_t* aa) {
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

int isExactlyEqual_AA(AltArr_t* a, AltArr_t* b) {
	if (a == NULL) {
		return (b == NULL);
	}
	if (b == NULL) {
		return 0;
	}

	if (a->size != b->size || a->nvar != b->nvar) {
		return 0;
	}

	if (a->unpacked || b->unpacked) {
		return isExactlyEqual_AA_unpk(a,b);
	}

	for (int i = 0; i < a->size; ++i) {
		if (a->elems[i].degs != b->elems[i].degs) {
			return 0;
		}
		if (mpq_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
			return 0;
		}
	}

	return 1;
}

void nonZeroVariables_AA(AltArr_t* aa, int* foundVar) {
    if (aa == NULL || aa->nvar == 0 || isConstant_AA(aa)) {
        return;
    }

    if (aa->unpacked) {
    	nonZeroVariables_AA_unpk(aa, foundVar);
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

degree_t totalDegree_AA(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0 || isConstant_AA(aa)) {
		return 0;
	}

	if (aa->unpacked) {
		return totalDegree_AA_unpk(aa);
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

degree_t partialDegree_AA(const AltArr_t* aa, int k) {
	if (aa == NULL || aa->size == 0 || k >= aa->nvar) {
		return -1;
	}

	if (aa->unpacked) {
		return partialDegree_AA_unpk(aa, k);
	}

	int nvar = aa->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
    const int* sizes = getExpOffsetArray(nvar);

	degree_t d = 0; degree_t dMax = 0;
	for (int i = 0; i < aa->size; ++i) {
		d = GET_NTH_EXP(aa->elems[i].degs, masks[k], sizes[k]);
		if (d > dMax) {
			dMax = d;
		}
	}

	return dMax;
}

void partialDegrees_AA(const AltArr_t* aa, degree_t* degsList) {
	if (aa == NULL || aa->size == 0 || aa->nvar == 0) {
		return;
	}

	if (aa->unpacked) {
		partialDegrees_AA_unpk(aa, degsList);
		return;
	}

	degrees_t degs = calculateMaxDegs_AA(aa);
	const degrees_t* masks = getExpMaskArray(aa->nvar);
    const int* sizes = getExpOffsetArray(aa->nvar);
    for (int k = 0; k < aa->nvar; ++k) {
    	degsList[k] = GET_NTH_EXP(degs, masks[k], sizes[k]);
    }
}

void setPartialDegreeTerm_AA(AltArr_t* aa, int idx, int k, degree_t deg) {
	if (aa == NULL || aa->size == 0) {
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
			unpackExponentVectors_AA_inp(aa);
			setPartialDegreeTerm_AA(aa, idx, k, deg);
			return;
		}
		aa->elems[idx].degs &= ~(masks[k]);
		aa->elems[idx].degs |= (((degrees_t) deg) << sizes[k]);
	}
	// mergeSortPolynomial_AA(aa);
}


degree_t mainDegree_AA(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return -1;
	}

	if (isConstant_AA(aa)) {
		return 0;
	}

	if (aa->unpacked) {
		return mainDegree_AA_unpk(aa);
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

int mainVariable_AA(AltArr_t* aa) {
	if (aa == NULL || aa->nvar == 0) {
		return -1;
	}

	if (aa->unpacked) {
		return mainVariable_AA_unpk(aa);
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

void coefficient_AA(AltArr_t* aa, const degree_t* degsList, int nvar, mpq_t retCoef) {
	if (isZero_AA(aa)) {
		mpq_set_si(retCoef, 0l, 1l);
		return;
    }
    if (nvar == 0 && aa->nvar == 0) {
        mpq_set(retCoef, aa->elems->coef);
    	return;
	}

	if (aa->unpacked) {
		return coefficient_AA_unpk(aa, degsList, nvar, retCoef);
	}

    const int* sizes = getExpOffsetArray(nvar);
    const degrees_t* maxExps = getMaxExpArray(nvar);

    degrees_t degs = 0;
    for (int j = 0; j < nvar; ++j) {
        if (degsList[j] > maxExps[j]) {
        	//if we are not unpacked and we are requesting a coefficent for a monomial
        	//which cannot be packed then obviously this coefficient does not exist;
        	mpq_set_ui(retCoef, 0l, 1l);
        	return;
        }
        degs |= ((degrees_t)degsList[j] << sizes[j]);
    }

	mpq_set_ui(retCoef, 0l, 1l);
    for (int i = 0; i < aa->size; ++i) {
        if (isEqualExponentVectors(degs, aa->elems[i].degs)) {
        	mpq_set(retCoef, aa->elems[i].coef);
            break;
        }
        if (isLessExponentVectors(aa->elems[i].degs, degs)) {
            //due to ordering of terms, none will have degs as monomial.
            break;
        }
    }
}

void setCoefficient_AA(AltArr_t* aa, const degree_t* degsList, int nvar, const mpq_t coef) {
    if (aa == NULL || nvar != aa->nvar) {
    	return;
    }

    if (aa->size == 0 || nvar == 0) {
    	if (aa->alloc < 1) {
    		resizePolynomial_AA(aa, 1);
    	}
        mpq_set(aa->elems->coef, coef);
        if (nvar != 0) {
            setDegrees_AA_inp(aa, 0, degsList, nvar);
        }
    	return;
    }

    if (aa->unpacked) {
    	setCoefficient_AA_unpk(aa, degsList, nvar, coef);
    	return;
    }

    const int* sizes = getExpOffsetArray(nvar);
    const degrees_t* maxExps = getMaxExpArray(nvar);

    degrees_t degs = 0;
    for (int j = 0; j < nvar; ++j) {
        if (degsList[j] > maxExps[j]) {
        	unpackExponentVectors_AA_inp(aa);
            setCoefficient_AA_unpk(aa, degsList, nvar, coef);
            return;
        }
        degs |= ((degrees_t) degsList[j] << sizes[j]);
    }

    for (int i = 0; i < aa->size; ++i) {
        if (isEqualExponentVectors(aa->elems[i].degs, degs)) {
            mpq_set(aa->elems[i].coef, coef);
            return;
        }
        if (isLessExponentVectors(aa->elems[i].degs, degs)) {
            //shift everything to the right 1 and insert at i.
            if (aa->size >= aa->alloc) {
                resizePolynomial_AA(aa, aa->size + 1);
            }
            mpq_init(aa->elems[aa->size].coef);
            for (int j = aa->size; j > i; --j) {
                mpq_swap(aa->elems[j].coef, aa->elems[j-1].coef);
                aa->elems[j].degs = aa->elems[j-1].degs;
            }
            aa->elems[i].degs = degs;
            mpq_set(aa->elems[i].coef, coef);
            ++(aa->size);
            return;
        }
    }

    //if we get here then we need to insert at the end of the array;

    if (aa->size >= aa->alloc) {
        resizePolynomial_AA(aa, aa->size + 1);
    }
    mpq_init(aa->elems[aa->size].coef);
    mpq_set(aa->elems[aa->size].coef, coef);
    aa->elems[aa->size].degs = degs;
    ++(aa->size);
}

void addTerm_AA(AltArr_t* aa, const degree_t* degsList, int nvar, const mpq_t coef) {
	if (aa == NULL || aa->nvar != nvar) {
		return;
	}

	if (aa->unpacked) {
		addTerm_AA_unpk(aa, degsList, nvar, coef);
		return;
	}


    degrees_t degs = 0;

    if (nvar > 0) {
	    const int* sizes = getExpOffsetArray(nvar);
	    const degrees_t* maxExps = getMaxExpArray(nvar);

	    for (int j = 0; j < nvar; ++j) {
	        if (degsList[j] > maxExps[j]) {
	        	unpackExponentVectors_AA_inp(aa);
	            addTerm_AA_unpk(aa, degsList, nvar, coef);
	            return;
	        }
	        degs |= ((degrees_t) degsList[j] << sizes[j]);
	    }

	    if (aa->size >= aa->alloc) {
	    	resizePolynomial_AA(aa, aa->alloc*2);
	    }
    }

	mpq_init(aa->elems[aa->size].coef);
	mpq_set(aa->elems[aa->size].coef, coef);
	aa->elems[aa->size].degs = degs;
	++(aa->size);
}

int isEqualWithVariableOrdering_AA(AltArr_t* a, AltArr_t*b, const int* xs, int xsSize) {
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

	if (a->unpacked || b->unpacked) {
		return isEqualWithVariableOrdering_AA_unpk(a, b, xs, xsSize);
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
        if (mpq_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
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

AltArr_t* termAtIdx_AA(AltArr_t* a, int idx) {
	if (a == NULL || idx >= a->size) {
		return NULL;
	}

	AltArr_t localA = *a;
	localA.elems = localA.elems + idx;
	localA.size = 1;
	return deepCopyPolynomial_AA(&localA);
}

int isConstantTermZero_AA(AltArr_t* a) {
	if (a == NULL || a->size == 0) {
		return 1;
	}

	if (a->unpacked) {
		return !isZeroExponentVector_unpk(a->elems[a->size-1].degs, a->nvar);
	} else {
		return !isZeroExponentVector(a->elems[a->size-1].degs);
	}
}

void expandNumVars_AA(AltArr_t* aa_in, int newNvar) {
	if (aa_in == NULL) {
		return;
	}
	if (aa_in->alloc == 0) {
		aa_in->nvar = newNvar;
		return;
	}
	if (newNvar <= aa_in->nvar) {
		return;
	}

	if (aa_in->unpacked) {
		expandNumVars_AA_unpk(aa_in, newNvar);
		return;
	}

	if (aa_in->nvar == 0) {
		aa_in->nvar = newNvar;
		return;
	}

	AltArr_t* aa = deepCopyPolynomial_AA(aa_in);

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
				// TODO
				freePolynomial_AA(aa);
				unpackExponentVectors_AA_inp(aa_in);
				expandNumVars_AA_unpk(aa_in, newNvar);
				return;
				// int diff = newNvar - aa->nvar;
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for expand at index %d; %lld > %lld.\n", j+diff, curDeg, maxExps[j+diff]);
				// exit(1);
			}
			aa->elems[i].degs |= (curDeg << newSizes[j]);
		}
	}

	aa_in->nvar = newNvar;
	AAElem_t* tempElems = aa_in->elems;
	aa_in->elems = aa->elems;
	aa->elems = tempElems;

	freePolynomial_AA(aa);
}

void expandNumVarsLeft_AA(AltArr_t* aa_in, int newNvar) {
	if (aa_in == NULL) {
		return;
	}
	if (aa_in->size == 0|| aa_in->alloc == 0) {
		aa_in->nvar = newNvar;
		return;
	}
	if (newNvar <= aa_in->nvar) {
		return;
	}

	if (aa_in->unpacked) {
		expandNumVarsLeft_AA_unpk(aa_in, newNvar);
		return;
	}

	if (aa_in->nvar == 0) {
		aa_in->nvar = newNvar;
		return;
	}


	AltArr_t* aa = deepCopyPolynomial_AA(aa_in);

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
				//then give up and do so in unpacked world.
				// TODO
				freePolynomial_AA(aa);
				unpackExponentVectors_AA_inp(aa_in);
				expandNumVarsLeft_AA_unpk(aa_in, newNvar);
				return;
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for expand at index %d; %lld > %lld.\n", j+diff, curDeg, maxExps[j+diff]);
				// exit(1);
			}
			aa->elems[i].degs |= (curDeg << newSizes[j+diff]);
		}
	}

	aa_in->nvar = newNvar;
	AAElem_t* tempElems = aa_in->elems;
	aa_in->elems = aa->elems;
	aa->elems = tempElems;

	freePolynomial_AA(aa);
}

void shrinkNumVarsAtIdx_AA(AltArr_t* aa, int idx) {
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
		shrinkNumVarsAtIdx_AA_unpk(aa, idx);
		return;
	}

	int nvar = aa->nvar;
	int newNvar = aa->nvar - 1;

	AAElem_t* elems = aa->elems;
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

void shrinkAndReorderVars_AA(AltArr_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->nvar < 1) {
		return;
	}

	if (varmapSize > aa->nvar) {
		return;
	}

	if (aa->unpacked) {
		shrinkAndReorderVars_AA_unpk(aa, varMap, varmapSize);
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
			mpq_clear(aa->elems[i].coef);
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

	AltArr_t* tmpAA = deepCopyPolynomial_AA(aa);
	AAElem_t* elems = tmpAA->elems;
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
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for shrink and reorder vars. At new index %d.", varMap[j]);
				// exit(1);
				//TODO
				freePolynomial_AA(tmpAA);
				unpackExponentVectors_AA_inp(aa);
				shrinkAndReorderVars_AA_unpk(aa, varMap, varmapSize);
				return;
			}
			newDegs |= (curDeg << newSizes[varMap[j]]);
		}
		elems[i].degs = newDegs;
	}

	AAElem_t* tmpElems = aa->elems;
	aa->elems = elems;
	tmpAA->elems = tmpElems;
	freePolynomial_AA(tmpAA);

	aa->nvar = newNvar;

	if (needSort) {
		mergeSortPolynomial_AA(aa);
	}
}

//     degs[varmap[i]] = curDegs[i];
void reorderVars_AA(AltArr_t* aa_in, int* varMap, int varmapSize) {
	if (aa_in == NULL || aa_in->size == 0) {
		return;
	}

	if (aa_in->unpacked) {
		reorderVars_AA_unpk(aa_in, varMap, varmapSize);
		return;
	}

	int nvar = aa_in->nvar;

	const degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	const int* __restrict__ sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	AltArr_t* aa = deepCopyPolynomial_AA(aa_in);

	degrees_t degs;
	degrees_t curDeg;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < varmapSize; ++j) {
			curDeg = GET_NTH_EXP(degs, masks[j], sizes[j]);
			if (curDeg > maxExps[varMap[j]]) {
				// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for reorder at index %d; %lld > %lld.\n", varMap[j], curDeg, maxExps[varMap[j]]);
				// exit(1);
				//TODO give up and unpack
				freePolynomial_AA(aa);
				unpackExponentVectors_AA_inp(aa_in);
				reorderVars_AA_unpk(aa_in, varMap, varmapSize);
				return;
			}
			aa->elems[i].degs |= (curDeg << sizes[varMap[j]]);
		}
	}

	AAElem_t* tmpElems = aa_in->elems;
	aa_in->elems = aa->elems;
	aa->elems = tmpElems;
	freePolynomial_AA(aa);


	mergeSortPolynomial_AA(aa_in);

}


int tryShrinkVariables_AA_inp(AltArr_t* aa, int* varMap) {
	int trueNvar = computeShrinkMap_AA(aa, varMap);
	if (trueNvar != aa->nvar) {
		shrinkAndReorderVars_AA(aa, varMap, aa->nvar);
	}
	return trueNvar;
}

int reverseShrinkVariables_AA_inp(AltArr_t* aa, int varmapSize, int* varMap) {
	if (aa == NULL || varMap == NULL || aa->nvar == varmapSize) {
		return aa == NULL ? 0 : aa->nvar;
	}

	int revMap[varmapSize];
	reverseShrinkMap_AA(varmapSize, aa->nvar, varMap, revMap);

	expandNumVars_AA(aa, varmapSize);
	reorderVars_AA(aa, revMap, varmapSize);
	return varmapSize;
}

void setDegrees_AA_inp(AltArr_t* aa, int idx, const degree_t* degsList, int nvar) {
	if (aa == NULL || idx >= aa->alloc) {
		return;
	}

	if (aa->unpacked) {
		setDegrees_AA_inp_unpk(aa, idx, degsList, nvar);
		return;
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxDegs = getMaxExpArray(nvar);

	aa->elems[idx].degs = 0;

	for (int k = 0; k < nvar; ++k) {
		if (degsList[k] > maxDegs[k]) {
			unpackExponentVectors_AA_inp(aa);
			setDegrees_AA_inp_unpk(aa, idx, degsList, nvar);
			break;

			// fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for setDegrees_AA at index %d; %d > %lld.\n", k, degsList[k], maxDegs[k]);
			// exit(1);
		}

		aa->elems[idx].degs |= ((degrees_t) degsList[k] << sizes[k]);
	}

}

AltArr_t* deepCopyPolynomial_AAFromNode(Node* a, int nvar) {
	if (a == NULL) {
		return NULL;
	}
	polysize_t asize = numberOfTermsNode(a);
	AltArr_t* newAA = makePolynomial_AA(asize, nvar);
	newAA->size = asize;

	int size = AA_SIZE(newAA);
	AAElem_t* elems = newAA->elems;
	for (int i = 0; i < size && a != NULL; ++i) {
		mpq_init(elems[i].coef);
		mpq_set(elems[i].coef, a->coef);
		elems[i].degs = a->degs;
		a = a->next;
	}

	return newAA;
}


Node* deepCopyPolynomial_NodeFromAA(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		fprintf(stderr, "Cannot make Node from unpacked AltArr_t\n");
		exit(1);
	}

	int asize = AA_SIZE(aa);

	AAElem_t* elems = aa->elems;
	Node* head = addTerm(NULL, elems[0].degs, elems[0].coef);
	Node* tail = head;
	for (int i = 1; i < asize; ++i) {
		tail = addTerm(tail, elems[i].degs, elems[i].coef);
	}

	return head;
}

AltArrDegList_t* deepCopyPolynomial_AADegListFromAA(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return deepCopyPolynomial_AADegListFromAA_unpk(aa);
	}

	int nvar = aa->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	AltArrDegList_t* ret = makePolynomial_AADL(aa->alloc, nvar);
	AAElem_DegList_t* rElems = ret->elems;
	AAElem_t* elems = aa->elems;

	for (int i = 0; i < aa->size; ++i) {
		mpq_init(rElems[i].coef);
		mpq_set(rElems[i].coef, elems[i].coef);
		rElems[i].degs = malloc(sizeof(degree_t)*nvar);
		for (int j = 0; j < nvar; ++j) {
			rElems[i].degs[j] = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
		}
	}

	ret->size = aa->size;

	return ret;
}

AltArr_t* deepCopyPolynomial_AA(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return deepCopyPolynomial_AA_unpk(aa);
	}

	AltArr_t* newAA = makePolynomial_AA(aa->alloc, aa->nvar);
	int size = AA_SIZE(newAA) = aa->size;
	AAElem_t* elems = newAA->elems;
	AAElem_t* oldelems = aa->elems;
	for (int i = 0; i < size; ++i) {
		mpq_init(elems[i].coef);
		mpq_set(elems[i].coef, oldelems[i].coef);
		elems[i].degs = oldelems[i].degs;
	}

	return newAA;
}

AltArr_t* sortPolynomial_AA(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return sortPolynomial_AA_unpk(aa);
	}

	AAElem_t* elems = aa->elems;
	int size = AA_SIZE(aa);

	degrees_t swapDegs;
	for (int i = 1; i < size; ++i) {
		for (int j = i; j > 0 && compareExponentVectors(elems[j-1].degs, elems[j].degs) < 0; --j) {
			mpq_swap(elems[j-1].coef, elems[j].coef);
			swapDegs = elems[j-1].degs;
			elems[j-1].degs = elems[j].degs;
			elems[j].degs = swapDegs;
		}
	}

	condensePolyomial_AA(aa);

	return aa;
}


static void mergeAAElems(AAElem_t* __restrict__ a, AAElem_t* __restrict__ endA, AAElem_t* __restrict__ b , AAElem_t* __restrict__ endB, AAElem_t* __restrict__ sorted) {
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

void mergeSortPolynomial_AA(AltArr_t* aa) {
	if (aa == NULL || aa-> size == 0) {
		return;
	}

	if (aa->unpacked) {
		mergeSortPolynomial_AA_unpk(aa);
		return;
	}

	if (aa->size < 8) {
		sortPolynomial_AA(aa);
		return;
	}

	int size = aa->size;
	AAElem_t* tempElems = (AAElem_t*) malloc(sizeof(AAElem_t)*aa->alloc);
	AAElem_t* elems = aa->elems;
	int end1 = 0;
	int end2 = 0;
	for (int window = 1; window < size; window <<= 1) {
		//merge elems[i]:elems[i+window] with elems[i+window]:elems[i+2*window]
		for (int i = 0; i < size; i += 2*window) {
			end1 = i + window < size ? i + window : size;
			end2 = i + 2*window < size ? i + 2*window : size;
			mergeAAElems(elems+i, elems+end1, elems+end1, elems+end2, tempElems+i);
		}

		AAElem_t* temp = tempElems;
		tempElems = elems;
		elems = temp;
	}

	aa->elems = elems;

	free(tempElems);

	condensePolyomial_AA(aa);
}

void condensePolyomial_AA(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) < 1) {
		return;
	}

	if (aa->unpacked) {
		condensePolyomial_AA_unpk(aa);
		return;
	}

	int size = AA_SIZE(aa);
	AAElem_t* elems = aa->elems;
	int insertIdx = 0;
	int compareIdx = 1;
	while (compareIdx < size) {
		if (compareExponentVectors(elems[insertIdx].degs, elems[compareIdx].degs) == 0) {
			mpq_add(elems[insertIdx].coef, elems[insertIdx].coef, elems[compareIdx].coef);
		} else if (mpq_sgn(elems[insertIdx].coef) == 0) {
			//we are about to move to next degree, but the coef is 0, so we
			//need to remove it. Do so by setting insertIdx's degree to the
			//comare idx so it can get overwritten by next compare and add.
			elems[insertIdx].degs = elems[compareIdx].degs;
			--compareIdx;
		} else if(compareIdx - insertIdx > 1) {
			++insertIdx;
			elems[insertIdx].degs = elems[compareIdx].degs;
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

void negatePolynomial_AA(AltArr_t* aa) {
	int size = AA_SIZE(aa);
	AAElem_t* elems = aa->elems;
	for (int i = 0; i < size; ++i) {
		mpq_neg(elems[i].coef, elems[i].coef);
	}
}

AltArr_t* mainLShiftPolynomial_AA (AltArr_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (n < 1){
		return deepCopyPolynomial_AA (aa);
    }

    if (aa->unpacked) {
    	return mainLShiftPolynomial_AA_unpk(aa, n);
    }

    int nvar = aa->nvar;
    int mvarDegOffset = getMVarExpOffset (nvar);

    degrees_t nDegs = ((degrees_t) n) << mvarDegOffset;
    int fitsPacked = monomialMultFitsPacked(aa->elems->degs, nDegs, nvar);
    if (!fitsPacked) {
    	unpackExponentVectors_AA_inp(aa);
    	AltArr_t* ret = mainLShiftPolynomial_AA_unpk(aa, n);
    	packExponentVectors_AA_inp(aa);
    	return ret;
    }
    AltArr_t* shifted = (AltArr_t*) malloc (sizeof (AltArr_t));
    shifted->nvar = aa->nvar;
    shifted->size = aa->size;
    shifted->alloc = aa->alloc;
    shifted->unpacked = 0;
    shifted->elems = (AAElem_t*) malloc (sizeof (AAElem_t)*(aa->size));
    AAElem_t* selems = shifted->elems;
    AAElem_t* aelems = aa->elems;

    for (int i = 0; i < shifted->size; ++i){
		mpq_init (selems[i].coef);
		mpq_set (selems[i].coef, aelems[i].coef);
		selems[i].degs = aelems[i].degs;
		selems[i].degs += (((degrees_t) n) << mvarDegOffset);
    }
    return shifted;
}

AltArr_t* mainLShiftPolynomial_AA_inp (AltArr_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (n < 1){
		return aa;
    }

    if (aa->unpacked) {
    	return mainLShiftPolynomial_AA_inp_unpk(aa, n);
    }

    int nvar = aa->nvar;
    int mvarDegOffset = getMVarExpOffset (nvar);

    degrees_t nDegs = ((degrees_t) n) << mvarDegOffset;
    int fitsPacked = monomialMultFitsPacked(aa->elems->degs, nDegs, nvar);
    if (!fitsPacked) {
    	unpackExponentVectors_AA_inp(aa);
    	return mainLShiftPolynomial_AA_inp_unpk(aa, n);
    }

    AAElem_t* aelems = aa->elems;

    for (int i = 0; i < aa->size; ++i){
		aelems[i].degs += ((degrees_t) n << mvarDegOffset);
    }
    return aa;
}

void evalPolyToVal_AA(const AltArr_t* aa, ratNum_t* vals, short nvar, ratNum_t res) {
	if (aa == NULL || nvar != aa->nvar) {
		mpq_set_ui(res, 0ul, 1ul);
		return;
	}

	if (nvar == 0 || aa->nvar == 0) {
		mpq_set(res, aa->elems->coef);
		return;
	}

	if (aa->unpacked) {
		return evalPolyToVal_AA_unpk(aa, vals, nvar, res);
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);

	mpq_t* valList[nvar];
	int valListSize[nvar];

	degrees_t maxDegs = calculateMaxDegs_AA(aa);
	for (int j = 0; j < nvar; ++j) {
		degree_t deg = GET_NTH_EXP(maxDegs, masks[j], sizes[j]);
		valList[j] = (mpq_t*) malloc(sizeof(mpq_t)*(deg+1));
		valListSize[j] = deg+1;

		mpq_init(valList[j][0]);
		mpq_set_ui(valList[j][0], 1ul, 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpq_init(valList[j][k]);
			mpq_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	mpq_set_ui(res, 0ul, 1ul);
	mpq_t acc;
	mpq_init(acc);

	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpq_set(acc, aa->elems[i].coef);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = GET_NTH_EXP(aa->elems[i].degs, masks[j], sizes[j]);
			mpq_mul(acc, acc, valList[j][deg]);
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

}

AltArr_t* evaluatePoly_AA(AltArr_t* aa, int* active, ratNum_t* vals, short nvar) {
	if (aa == NULL) {
		return NULL;
	}

	if (nvar == 0 || aa->nvar == 0) {
		return deepCopyPolynomial_AA(aa);
	}

	if (aa->unpacked) {
		return evaluatePoly_AA_unpk(aa, active, vals, nvar);
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

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);

	mpq_t* valList[nvar];
	int valListSize[nvar];

	degrees_t maxDegs = calculateMaxDegs_AA(aa);
	for (int j = 0; j < nvar; ++j) {
		degree_t deg = GET_NTH_EXP(maxDegs, masks[j], sizes[j]);
		if (!active[j]) {
			valList[j] = NULL;
			valListSize[j] = 0;
			continue;
		}
		valList[j] = (mpq_t*) malloc(sizeof(mpq_t)*(deg+1));
		valListSize[j] = deg+1;

		mpq_init(valList[j][0]);
		mpq_set_ui(valList[j][0], 1ul, 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpq_init(valList[j][k]);
			mpq_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	AltArr_t* res = makePolynomial_AA(aa->size, newNvar);
	AAElem_t* resElems = res->elems;
	res->size = aa->size;

	const int* newSizes = getExpOffsetArray(newNvar);
	int size = aa->size;
	int k = 0;
	degrees_t newDegs = 0;
	for (int i = 0; i < size; ++i) {
		mpq_init(resElems[i].coef);
		mpq_set(resElems[i].coef, aa->elems[i].coef);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = GET_NTH_EXP(aa->elems[i].degs, masks[j], sizes[j]);
			if (valList[j] == NULL) {
				newDegs |= ((degrees_t) deg) << newSizes[k];
				++k;
			} else {
				mpq_mul(resElems[i].coef, resElems[i].coef, valList[j][deg]);
			}
		}
		resElems[i].degs = newDegs;
		k = 0;
		newDegs = 0;
	}

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpq_clear(valList[j][k]);
		}
		free(valList[j]);
	}

	mergeSortPolynomial_AA(res);

	return res;
}

AltArr_t* convertFromAAElemToAA (AAElem_t* coef, int coefSize, int nvar, int unpacked)
{
    if (coef == NULL || coefSize == 0) {
      return NULL;
    }

    AltArr_t* poly = makePolynomial_AA(coefSize, nvar);
    poly->size = coefSize;
    poly->elems = coef;
    poly->unpacked = unpacked;
    return poly;
}

AltArr_t* swappingExponents_AA (AltArr_t* aa, int idx1, int idx2)
{
  if (aa == NULL) {
      return NULL;
  }

  // As long as we call reorderVars_AA below, then we do not need
  // a _unpk version of this method. Reorder vars handles it.
  // if (aa->unpacked) {
  // 	return swappingExponents_AA_unpk(aa, idx1, idx);
  // }

  AltArr_t* cPoly = deepCopyPolynomial_AA(aa);
  if (idx1 == idx2){
      return cPoly;
  }
  int varMap[cPoly->nvar];
  for (int i = 0; i < cPoly->nvar; ++i){
    varMap[i] = i;
  }
  varMap[idx1] = idx2;
  varMap[idx2] = idx1;

  reorderVars_AA (cPoly, varMap, cPoly->nvar);
  return cPoly;
}

AltArr_t* leadingTerm_AA (AltArr_t* aa, int nvar)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	if (aa->unpacked) {
		return leadingTerm_AA_unpk(aa, nvar);
	}

	AltArr_t* lt = makeConstPolynomial_AA(1, aa->nvar, aa->elems->coef);
	lt->elems->degs = aa->elems->degs;

	return lt;
}

int leadingVariable_AA (AltArr_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return -2;
	}

	if (aa->unpacked) {
		return leadingVariable_AA_unpk(aa);
	}

	const degrees_t* masks = getExpMaskArray(aa->nvar);
	for (int i = 0; i < aa->nvar; ++i){
		if ((aa->elems[0].degs & masks[i]) != 0) {
			return i;
		}
	}

	return -1;
}


int mainLeadingDegree_AA (AltArr_t* aa)
{
    if (aa == NULL || aa->size == 0){
		return 0;
    }

    if (aa->unpacked) {
    	return mainLeadingDegree_AA_unpk(aa);
    }

    int mvarDegOffset = getMVarExpOffset (aa->nvar);
    degrees_t mvarMask = getMVarExpMask (aa->nvar);

    return GET_NTH_EXP(aa->elems->degs, mvarMask, mvarDegOffset);
}

AltArr_t* mainLeadingCoefficient_AA (AltArr_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	if (aa->unpacked) {
		return mainLeadingCoefficient_AA_unpk(aa);
	}

	int mvarDegOffset = getMVarExpOffset(aa->nvar);
	degrees_t mvarMask = getMVarExpMask(aa->nvar);

	AAElem_t* elems = aa->elems;
	degree_t mvarDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
	degree_t curDeg = mvarDeg;

	int polyElemsAlloc = 10;
	int polyElemsSize = 0;
	AAElem_t* polyElems = (AAElem_t*) malloc (sizeof(AAElem_t) * polyElemsAlloc);

	for (int i = 0; (mvarDeg == curDeg) && i < AA_SIZE (aa); ++i){

		if (polyElemsSize + 1 > polyElemsAlloc){
			polyElemsAlloc += 10;
			polyElems = (AAElem_t*) realloc (polyElems, sizeof (AAElem_t) * polyElemsAlloc);
		}

		mpq_init (polyElems[i].coef);
		mpq_set (polyElems[i].coef, elems[i].coef);
		polyElems[i].degs = (elems[i].degs & (~mvarMask));
		++polyElemsSize;

		if (i+1 < AA_SIZE(aa)){
			curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
		} else {
			curDeg = -1;
		}
	}

	AltArr_t* poly = (AltArr_t*) malloc (sizeof (AltArr_t));
	poly->alloc = polyElemsAlloc;
	poly->size = polyElemsSize;
	poly->elems = polyElems;
	poly->unpacked = 0;
	poly->nvar = aa->nvar;

	return poly;
}

AltArr_t* mainCoefficientAtIdx_AA (AltArr_t* aa, int e)
{
    if (aa == NULL || aa->size == 0){
		return NULL;
    }
    if (e < 0){
		return NULL;
    }

    if (aa->unpacked) {
    	return mainCoefficientAtIdx_AA_unpk(aa, e);
    }

    int mvarDegOffset = getMVarExpOffset (aa->nvar);
    degrees_t mvarMask = getMVarExpMask (aa->nvar);
    AAElem_t* elems = aa->elems;
    degree_t mvarDeg = e;
    degree_t curDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
    int begin = -1;

    if (mvarDeg > curDeg){
        return NULL;
    }

    int polyElemsAlloc = 10;
    int polyElemsSize = 0;
    AAElem_t* polyElems = (AAElem_t*) malloc (sizeof(AAElem_t) * polyElemsAlloc);

    for (int i = 0; i < AA_SIZE (aa); ++i){
	if (mvarDeg == curDeg){
	    if (polyElemsSize + 1 > polyElemsAlloc) {
		polyElemsAlloc += 10;
		polyElems = (AAElem_t*) realloc (polyElems, sizeof (AAElem_t)*polyElemsAlloc);
	    }

	    mpq_init (polyElems[polyElemsSize].coef);
	    mpq_set (polyElems[polyElemsSize].coef, elems[i].coef);
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

    AltArr_t* poly = (AltArr_t*) malloc(sizeof(AltArr_t));
    poly->alloc = polyElemsAlloc;
    poly->size = polyElemsSize;
    poly->elems = polyElems;
    poly->unpacked = 0;
    poly->nvar = aa->nvar;

    return poly;
}

void mainCoefficientListAtIdx_AA (AltArr_t* aa, int idx, AltArr_t*** cList, int *sz)
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
		fprintf (stderr, "mainCoefficientListAtIdx_AA_unpk is not implemented yet!\n");
		exit(1);
	}


	AltArr_t* AA = aa;
	int isExpand = 0;
	int nnvar;

	if (idx > 0) {
		AA = aa; //deepCopyPolynomial_AA (aa);

		nnvar = AA->nvar + 1;
		expandNumVarsLeft_AA (AA, nnvar);

		int varmap [nnvar];
		for (int i = 0; i < nnvar; i++) {
			varmap[i] = i;
		}
		varmap[0] = idx+1;
		varmap[idx+1] = 0;
		reorderVars_AA (AA, varmap, nnvar);

		isExpand = 1;
	} else {
		nnvar = AA->nvar;
	}

    int mvarDegOffset = getMVarExpOffset (nnvar);
    degrees_t mvarMask = getMVarExpMask (nnvar);

	AAElem_t* elems = AA->elems;
	degree_t mvarDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
	degree_t curDeg =  mvarDeg;
	degree_t expDeg;

	AltArr_t** mCoefs = (AltArr_t**) calloc (curDeg+1, sizeof (AltArr_t*));
	int len = 0;

	AltArr_t* poly;
	int polyElemsAlloc;
	int polyElemsSize;
	AAElem_t* polyElems;
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
		polyElems = (AAElem_t*) malloc (polyElemsAlloc * sizeof(AAElem_t)); // TODO
		for (int i = start; (curDeg == expDeg) && i < AA_SIZE(AA); i++) {

			if (polyElemsSize + 1 > polyElemsAlloc){
				polyElemsAlloc += 10;
				polyElems = (AAElem_t*) realloc (polyElems, sizeof (AAElem_t) * polyElemsAlloc);
			}

			/* fprintf (stderr, "[SMQP] i = %d\n", i); // TEST */
			mpq_init (polyElems[polyElemsSize].coef);
			mpq_set (polyElems[polyElemsSize].coef, elems[i].coef);
			polyElems[polyElemsSize].degs = (elems[i].degs & (~mvarMask));
			++polyElemsSize;

			if (i+1 < AA_SIZE(AA)){
				curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
			} else {
				curDeg = -1;
			}

			start++;
		}

		poly = (AltArr_t*) malloc (sizeof (AltArr_t));
		poly->alloc = polyElemsAlloc;
		poly->size = polyElemsSize;
		poly->elems = polyElems;
		poly->unpacked = 0;
		poly->nvar = nnvar;

		if (isExpand) {
		    shrinkNumVarsAtIdx_AA (poly, 0);
		}

		/* printAA (poly); // TEST		 */
		mCoefs[expDeg] = poly;
		len++;
	}

	*cList = mCoefs;
	*sz = len;

}

// TODO: used in recursive MDPD
AltArr_t* maxPolynomials_AA (AltArr_t* a, AltArr_t* b)
{
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AA(a);
	}

	if (a->unpacked || b->unpacked) {
		return maxPolynomials_AA_unpk(a,b);
	}

	int cmp = compareExponentVectors(a->elems[0].degs, b->elems[0].degs);

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

// TODO: used in recursive MDPD
AltArr_t* maxPolynomials_AA_inp (AltArr_t* a, AltArr_t* b)
{
	if (a == NULL || a->size == 0){
		if (b == NULL || b->size == 0){
			return NULL;
		}
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0){
		return a;
	}

	if (a->unpacked || b->unpacked) {
		return maxPolynomials_AA_inp_unpk(a, b);
	}

	int cmp = compareExponentVectors(a->elems[0].degs, b->elems[0].degs);

	if (cmp > 0){
		return a;
	} else if (cmp < 0){
		return deepCopyPolynomial_AA(b);
	} else {
		if (mpq_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
			return a;
		} else {
			return deepCopyPolynomial_AA(b);
		}
	}
	return NULL;
}

/*****************
 * SMQP Addition
 *****************/

/**
 * Add two polynomials given their head nodes, a and b.
 * nvar: number of variables in the polynomials.
 * returns a pointer to the head Node of the sum.
 * Note that Nodes and exponent vectors are reused as much as possible.
 */
Node* addPolynomials(Node *a, Node *b, int nvar) {
	Node* c = NULL;
	Node* trailC = c;
	Node* curA = a;
	Node* curB = b;

	ratNum_t ccoef;
	mpq_init(ccoef);
	while (curA != NULL && curB != NULL) {
		cmpExp_t cmp = compareExponentVectors(curA->degs, curB->degs);
		if (cmp < 0) {
			//a < b
			// degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			// memcpy(degs, curB->degs, sizeof(degree_t)*nvar);
			trailC = addTerm (trailC, curB->degs, curB->coef);
			curB = curB->next;
		} else if (cmp == 0) {
			// a==b
			mpq_add(ccoef, curA->coef, curB->coef);
			// ccoef = curA->coef + curB->coef;
			if (mpq_sgn(ccoef) != 0) {
				// degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
				// memcpy(degs, curA->degs, sizeof(degree_t)*nvar);
				trailC = addTerm (trailC, curA->degs, ccoef);
			}
			curA = curA->next;
			curB = curB->next;
		} else {
			//a > b
			// degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			// memcpy(degs, curA->degs, sizeof(degree_t)*nvar);
			trailC = addTerm (trailC, curA->degs, curA->coef);
			curA = curA->next;
		}
		if (c == NULL) {
			c = trailC;
		}
	}

	if (curA != NULL) {
		Node* tail = deepCopyPolynomial(curA, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
		// trailC = addTerm(trailC, curA->degs, curA->coef);
		// if (c == NULL) {
		// 	c = trailC;
		// }
		// curA = curA->next;
	}
	if (curB != NULL) {
		Node* tail = deepCopyPolynomial(curB, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
		// trailC = addTerm(trailC, curB->degs, curB->coef);
		// if (c == NULL) {
		// 	c = trailC;
		// }
		// curB = curB->next;
	}

	mpq_clear(ccoef);
	return c;
}

void addRationalNumber_AA_inp(AltArr_t* aa, const mpq_t coef) {
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

	if (aa->unpacked) {
		addRationalNumber_AA_inp_unpk(aa, coef);
		return;
	}

    if (isZeroExponentVector(aa->elems[aa->size-1].degs)) {
        mpq_add(aa->elems[aa->size-1].coef, aa->elems[aa->size-1].coef, coef);
        return;
    }

    if (aa->size >= aa->alloc) {
        resizePolynomial_AA(aa, aa->alloc+10);
    }

    mpq_init(aa->elems[aa->size].coef);
    mpq_set(aa->elems[aa->size].coef, coef);
    aa->elems[aa->size].degs = 0;
    ++(aa->size);
}

AltArr_t* addPolynomials_AA(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AA(a);
	}

	if (a->unpacked || b->unpacked) {
		return addPolynomials_AA_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArr_t* c = makePolynomial_AA(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			mpq_init(cElems[k].coef);
			mpq_add(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].degs = aElems[i].degs;
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
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;

	resizePolynomial_AA(c, k);
	return c;
}

AltArr_t* subPolynomials_AA(AltArr_t* a, AltArr_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		AltArr_t* ret = deepCopyPolynomial_AA(b);
		negatePolynomial_AA(ret);
		return ret;
	}
	if (b == NULL || b->size == 0) {
		return deepCopyPolynomial_AA(a);
	}

	if (a->unpacked || b->unpacked) {
		return subPolynomials_AA_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArr_t* c = makePolynomial_AA(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			mpq_init(cElems[k].coef);
			mpq_sub(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].degs = aElems[i].degs;
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
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;

	resizePolynomial_AA(c, k);

	return c;
}

Node* addPolynomials_inp(Node* a, Node* b, int nvar) {
	Node* c = NULL;
	Node* trailC = c;
	Node* curA = a;
	Node* curB = b;

	while (curA != NULL && curB != NULL) {
		cmpExp_t cmp = compareExponentVectors(curA->degs, curB->degs);
		if (cmp < 0) {
			//a < b
			if (trailC != NULL) {
				trailC->next = deepCopyNode(curB, nvar);
				trailC = trailC->next;
			} else {
				trailC = deepCopyNode(curB, nvar);
			}
			mpq_set(trailC->coef, trailC->coef);
			curB = curB->next;
		} else if (cmp == 0) {
			// a==b
			mpq_add(curA->coef, curA->coef, curB->coef);
			if (mpq_sgn(curA->coef) == 0) {
				Node* del = curA;
				curA = curA->next;
				freeNode(del);
			} else {
				if(trailC != NULL) {
					trailC->next = curA;
					trailC = trailC->next;
				} else {
					trailC = curA;
				}
				curA = curA->next;
			}
			curB = curB->next;
		} else {
			if (trailC == NULL) {
				trailC = curA;
			} else {
				trailC->next = curA;
				trailC = curA;
			}
			curA = curA->next;
		}

		if (c == NULL) {
			c = trailC;
		}
	}

	if (curA != NULL) {
		if (trailC != NULL) {
			trailC->next = curA;
		} else {
			c = curA;
		}
	}
	if (curB != NULL) {
		Node* tail = deepCopyPolynomial(curB, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
	}

	return c;
}

#define BAD_IN_PLACE 0

AltArr_t* addPolynomials_AA_inp(AltArr_t* a, AltArr_t* b, int nvar) {
#if BAD_IN_PLACE
	AltArr_t* c = addPolynomials_AA(a,b,nvar);
	freePolynomial_AA(a);
	return c;
#else
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		if (a != NULL && a->size == 0) {
			freePolynomial_AA(a);
		}
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0) {
		return a;
	}

	if (a->unpacked || b->unpacked) {
		return addPolynomials_AA_inp_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArr_t* c = makePolynomial_AA(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = aElems[i];
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
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAElem_t)*(asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;
	resizePolynomial_AA(c, k);
	return c;

#endif
}


Node* subPolynomials_inp(Node* a, Node* b, int nvar) {
	Node* c = NULL;
	Node* trailC = c;
	Node* curA = a;
	Node* curB = b;

	while (curA != NULL && curB != NULL) {
		cmpExp_t cmp = compareExponentVectors(curA->degs, curB->degs);
		if (cmp < 0) {
			//a < b
			if (trailC != NULL) {
				trailC->next = deepCopyNode(curB, nvar);
				trailC = trailC->next;
			} else {
				trailC = deepCopyNode(curB, nvar);
			}
			//negate since from b
			mpq_neg(trailC->coef, trailC->coef);

			curB = curB->next;
		} else if (cmp == 0) {
			// a==b
			mpq_sub(curA->coef, curA->coef, curB->coef);
			if (mpq_sgn(curA->coef) == 0) {
				Node* del = curA;
				curA = curA->next;
				freeNode(del);
			} else {
				if(trailC != NULL) {
					trailC->next = curA;
					trailC = trailC->next;
				} else {
					trailC = curA;
				}
				curA = curA->next;
			}
			curB = curB->next;
		} else {
			if (trailC == NULL) {
				trailC = curA;
			} else {
				trailC->next = curA;
				trailC = curA;
			}
			curA = curA->next;
		}

		if (c == NULL) {
			c = trailC;
		}
	}

	if (curA != NULL) {
		if (trailC != NULL) {
			trailC->next = curA;
		} else {
			c = curA;
		}
	}
	if (curB != NULL) {
		Node* tail = deepCopyPolynomial(curB, nvar);
		negatePolynomial(tail);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
	}

	return c;
}

AltArr_t* subPolynomials_AA_inp(AltArr_t* a, AltArr_t* b, int nvar) {

#if BAD_IN_PLACE
	AltArr_t* diff = subPolynomials_AA(a,b,nvar);
	freePolynomial_AA(a);
	return diff;
#else
	if (a == NULL && b == NULL) {
		return NULL;
	}
	if (a == NULL || a->size == 0) {
		if (a != NULL && a->size == 0) {
			freePolynomial_AA(a);
		}
		AltArr_t* ret = deepCopyPolynomial_AA(b);
		negatePolynomial_AA(ret);
		return ret;
	}
	if (b == NULL || b->size == 0) {
		return a;
	}

	if (a->unpacked || b->unpacked) {
		return subPolynomials_AA_inp_unpk(a, b, nvar);
	}

	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	AltArr_t* c = makePolynomial_AA(asize + bsize, nvar);

	AAElem_t* aElems = a->elems;
	AAElem_t* bElems = b->elems;
	AAElem_t* cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpq_init(cElems[k].coef);
			mpq_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = aElems[i];
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
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAElem_t)*(asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while(j < bsize) {
		mpq_init(cElems[k].coef);
		mpq_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;

	resizePolynomial_AA(c, k);
	return c;

#endif
}

void subByLeadingTerm_AA (AltArr_t** a, int nvar)
{
	AltArr_t* aa = *a;

	if (aa == NULL || aa->size == 0) {
		return;
	}

	if (AA_SIZE (aa) == 1) {
		freePolynomial_AA (aa);
		*a = NULL;
		return;
	}

	if (aa->unpacked) {
		AltArr_t* lt = leadingTerm_AA (aa, nvar);
		*a = subPolynomials_AA_inp (aa, lt, nvar);
		freePolynomial_AA (lt);
		return;
	}

	AltArr_t* ret = makePolynomial_AA (AA_ALLOC(aa), nvar);
	memcpy (ret->elems, aa->elems + 1, sizeof (AAElem_t)*(AA_SIZE(aa)-1));

	AA_SIZE (ret) = AA_SIZE (aa) - 1;

	mpq_clear (aa->elems->coef);
	free (aa);

	*a = ret;
}

/*****************
 * SMQP Multiplication & Helpers
 *****************/

#if SMQP_SUPPORT_DEBUG
/**
 * Print the product degrees_t currently in the heap.
 */
void prodheapPrint(ProductHeap* h, int nvar) {
	for (int i = 0; i < h->heapSize; ++i){
		fprintf(stderr, "( %lx,", (unsigned long int) h->elements[i]);
		printDegs(h->elements[i]->product->degs, nvar);
		productHeapElem* next = h->elements[i]->next;
		while (next != NULL) {
			fprintf(stderr, "; ");
			fprintf(stderr, "%lx,", (unsigned long int) next);
			printDegs(next->product->degs, nvar);
			next = next->next;
		}
		fprintf(stderr, ")");
		fprintf(stderr," ");
	}
	fprintf(stderr, "\n");
}

void prodheapPrint_AA(ProductHeap_AA* h) {
	polysize_t s = h->heapSize;
	fprintf(stderr, "size: %ld", s);
	for (int i = 0; i < h->heapSize; ++i){
		fprintf(stderr, "( %llx,", h->elements[i].degs);
		// printDegs(h->elements[i]->product->degs, nvar);
		ProductHeapChain_AA* next = h->elements[i].chain->next;
		while (next != NULL) {
			fprintf(stderr, "; ");
			fprintf(stderr, "%llx,", h->elements[i].degs);
			// printDegs(next->product->degs, nvar);
			next = next->next;
		}
		fprintf(stderr, ")");
		fprintf(stderr," ");
	}
	fprintf(stderr, "\n");
}

#endif

/**
 * Make an element for the product heap, combining nodes a and b as the
 * element's product.
 */
productHeapElem* prodheapMakeElement(ProductHeap* h, Node* a, Node* b) {
	productHeapElem* elem = (productHeapElem*) malloc(sizeof(productHeapElem));
	elem->a_i = a;
	elem->b = b;
	elem->product = multiplyTerms(a, b, h->nvar);
	elem->next = NULL;
	return elem;
}

/**
 * Cleanup memory initialized for the heap
 */
void prodheapFree(ProductHeap* h) {
	productHeapElem** elems = h->elements;
	for (int i = 0; i < h->heapSize; ++i) {
		prodheapFreeElement(elems[i]);
	}
	free(h->elements);
	free(h);
}

/**
 * Create an empty product heap.
 */
ProductHeap* prodheapCreate(int nvar) {
	ProductHeap* h = (ProductHeap*) malloc(sizeof(ProductHeap));
	h->elements = NULL;
	h->heapSize = 0;
	h->maxHeapSize = 0;
	h->nvar = nvar;
	return h;
}

/**
 * Initialize the product heap with the two polynomials to multiply, a and b.
 *
 * We know the maximum heap size is numTerms(a) as at most one of a_i*b
 * is in the heap at once.
 */
ProductHeap* prodheapInit(Node* a, Node* b, int nvar) {
	ProductHeap* h = prodheapCreate(nvar);

	polysize_t num_a = numberOfTermsNode(a);
	productHeapElem** elems = (productHeapElem**) malloc(sizeof(productHeapElem*) * num_a);
	elems[0] = prodheapMakeElement(h, a, b);

	int heapCurSize = 1;
	int heapMaxSize = num_a;
	h->elements = elems;
	h->heapSize = heapCurSize;
	h->maxHeapSize = heapMaxSize;
	return h;
}

ProductHeap_AA* prodheapInit_AA(AltArr_t* a, AltArr_t* b, int nvar) {
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	h->elements = (ProductHeapElem_AA*) malloc(sizeof(ProductHeapElem_AA)*AA_SIZE(a));
#if SMQP_INT_PRODHEAP
	h->elements[0].chain = prodheapMakeChain_AA(0, 0, NULL);
#else
	h->elements[0].chain = prodheapMakeChain_AA(a->elems, b->elems, NULL);
#endif
	addExponentVectors(a->elems->degs, b->elems->degs, h->elements->degs);
	h->heapSize = 1;
	h->maxHeapSize = AA_SIZE(a);
	h->nvar = nvar;
	return h;
}

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
void prodheapInsert(ProductHeap* h, productHeapElem* elem) {

	int s = h->heapSize;
	int N = h->maxHeapSize;
	productHeapElem** elems = h->elements;
#if SMQP_COUNT_CHAINS
	SMQP_INSERTS++;
#endif

	if (s == 0) {
		if (s >= N) {
			N = (s == 0) ? 1 : s * 2;
			h->elements = (productHeapElem**) realloc(h->elements, N*sizeof(productHeapElem*));
			elems = h->elements;
			h->maxHeapSize = N;
		}
		elems[0] = elem;
		h->heapSize = 1;
#if SMQP_COUNT_CHAINS
		SMQP_CHAINS++;
#endif
		return;
	}


	degrees_t prodDegs = elem->product->degs;

	//first check if we can chain off the root
	if (compareExponentVectors(elems[0]->product->degs, prodDegs) == 0) {
		elem->next = elems[0];
		elems[0] = elem;
#if SMQP_COUNT_CHAINS
		SMQP_CHAINS++;
#endif
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

	int i = (s-1)/2; //i is parent
	int j = s;       //j is current insertion point
	long long unsigned int path = 1;
	cmpExp_t cmp;
	while (j > 0) {
		cmp = compareExponentVectors(elems[i]->product->degs, prodDegs);
		if (cmp == 0) {
#if SMQP_COUNT_CHAINS
			SMQP_CHAINS++;
#endif
			elem->next = elems[i];
			elems[i] = elem;
			return;
		} else if (cmp < 0) {
			path <<= 1;
			if (!(j & 1)) {
				//set the trailing bit to 1 to distinguish left/right of path
				path += 1;
			}
			j = i;
			i = (i-1) / 2;
		} else { //cmp > 0
			break;
		}
	}

	//if we made it this far, then we are adding a new element to the heap.
	//resize if necessary
	if (s >= N) {
		N = (s == 0) ? 1 : s * 2;
		h->elements = (productHeapElem**) realloc(h->elements, N*sizeof(productHeapElem*));
		elems = h->elements;
		h->maxHeapSize = N;
	}

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	productHeapElem* temp;
	while (j <= s) {
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = 2*j + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
#if SMQP_INT_PRODHEAP
void prodheapInsert_AA(ProductHeap_AA* h, ProductHeapChain_AA* chain, register degrees_t degs) {
#else
void prodheapInsert_AA(ProductHeap_AA* h, ProductHeapChain_AA* chain) {
#endif

	register int s = h->heapSize;
	ProductHeapElem_AA* elems = h->elements;
#if SMQP_COUNT_CHAINS
	SMQP_INSERTS++;
#endif

#if SMQP_INT_PRODHEAP
#else
	register degrees_t degs;
	addExponentVectors(chain->a_i->degs, chain->b->degs, degs);
#endif

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
	ProductHeapElem_AA temp;
	ProductHeapElem_AA elem = {degs, chain};

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
productHeapElem* prodheapRemoveMax(ProductHeap* h) {
	if (h->heapSize == 0) {
		return NULL;
	}

	productHeapElem** elems = h->elements;
	productHeapElem* maxElem = elems[0];
	int i = 0, j = 1;
	int s = --(h->heapSize);



	//promote largest children
	while (j < s) {
		if (j+1 < s && compareExponentVectors(elems[j]->product->degs, elems[j+1]->product->degs) < 0) {
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (2*j) + 1;
	}
	//now place last element into i and swim up to make tree complete
	j = (i-1)/2;
	while(i > 0) {
		if (compareExponentVectors(elems[s]->product->degs, elems[j]->product->degs) < 0) {
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j-1)/2;
	}
	elems[i] = elems[s];
	return maxElem;
}

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
ProductHeapChain_AA* prodheapRemoveMax_AA(ProductHeap_AA* h) {
	ProductHeapElem_AA* elems = h->elements;
	ProductHeapChain_AA* maxElem = elems[0].chain;
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
Node* multiplyPolynomials(Node* a, Node* b, int nvar) {
	if (a == NULL || b == NULL) {
		return NULL;
	}

	ProductHeap* h = prodheapInit(a,b,nvar);
	// degrees_t cdegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
	degrees_t cdegs;
	ratNum_t ccoef;
	mpq_init(ccoef);

	Node* c = NULL;
	Node* trailC = NULL;

//Note that in this chained heap version of the algorithm the heap
//resuses its product nodes as elements get inserted back into the heap
//and we create additional nodes as needed for the result.
//This differs from the regular heap implementation (below) where the product nodes
//extracted from the heap are used in the result.
	productHeapElem* maxElem = NULL;
	productHeapElem* nextMaxElem = NULL;
	degrees_t* nextDegs;
	while ( (nextDegs = prodheapPeek(h)) != NULL) {
		// memcpy(cdegs, nextDegs, sizeof(degree_t)*nvar);
		cdegs = *nextDegs;

		while (nextDegs != NULL && compareExponentVectors(cdegs, *nextDegs) == 0) {
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			productHeapElem* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax(h);
			while (maxElem != NULL) {
				mpq_add(ccoef, ccoef, maxElem->product->coef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == b && maxElem->a_i->next != NULL) {
					productHeapElem* nextAHeap = prodheapMakeElement(h, maxElem->a_i->next, b);
					//push it to the head of the oldMaxElem
					nextAHeap->next = oldMaxElem;
					oldMaxElem = nextAHeap;
				}

				//cache next before freeing or overwriting
				nextMaxElem = maxElem->next;

				//If the extracted term has another product in the stream,
				//update the product and push onto the oldMaxElem chain
				maxElem->b = (maxElem->b)->next;
				if(maxElem->b != NULL) {
					// ratNum_t tempCoef;
					// mpq_init(tempCoef);
					// degrees_t d = (degrees_t) malloc(sizeof(degree_t)*nvar);
					addExponentVectors(maxElem->a_i->degs, maxElem->b->degs, maxElem->product->degs);
					mpq_mul(maxElem->product->coef, maxElem->a_i->coef, maxElem->b->coef);
					// mpq_clear(tempCoef);
					maxElem->next = oldMaxElem;
					oldMaxElem = maxElem;
				} else {
					//we are done with the maxElem productHeapElem
					prodheapFreeElement(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;

			nextDegs = prodheapPeek(h);
		}

		//Commit new term to the product.
		if (mpq_sgn(ccoef) != 0) {
			trailC = addTerm(trailC, cdegs, ccoef);
			if (c == NULL) {
				c = trailC;
			}

			//reset accumulator variables
			// cdegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			mpq_clear(ccoef);
			mpq_init(ccoef);
		}

		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			prodheapInsert(h,maxElem);
			maxElem = nextMaxElem;

		}
	}
	// free(cdegs);

	mpq_clear(ccoef);
	prodheapFree(h);

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
AltArr_t* multiplyPolynomials_AA(AltArr_t*  a, AltArr_t*  b, int nvar) {
	if (a == NULL || a->size == 0 || b == NULL || b->size == 0) {
		return NULL;
	}

	if (a->unpacked || b->unpacked) {
		return multiplyPolynomials_AA_unpk(a,b,nvar);
	}

	degrees_t aMax = calculateMaxDegs_AA(a);
	degrees_t bMax = calculateMaxDegs_AA(b);

	if (!monomialMultFitsPacked(aMax, bMax, nvar)) {
		//force unpk to run by unpacking here.
		unpackExponentVectors_AA_inp(a);
		unpackExponentVectors_AA_inp(b);
		AltArr_t* ret = multiplyPolynomials_AA_unpk(a,b,nvar);
		packExponentVectors_AA_inp(a);
		packExponentVectors_AA_inp(b);
		return ret;
	}

	// reorder to obtain smaller as a.
	if (b->size < a->size) {
		AltArr_t* temp = a;
		a = b;
		b = temp;
	}

	ProductHeap_AA* h = prodheapInit_AA(a,b,nvar);
	ratNum_t ccoef;
	mpq_init(ccoef);

	register unsigned long long int allocC = AA_SIZE(a)*AA_SIZE(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;
	AltArr_t* c = makePolynomial_AA(allocC, nvar);

	//k is c, i is a, j is b.
	register int k = 0;
	// register int i = 0;
	// register int j = 0;

	AAElem_t* __restrict__ cElems = c->elems;
	AAElem_t* __restrict__ bElems = b->elems;

#if SMQP_INT_PRODHEAP
	AAElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1;
	register int firstB = 0;
#else
	AAElem_t* __restrict__ lastA = &(a->elems[AA_SIZE(a)-1]);
	AAElem_t* __restrict__ lastB = &(b->elems[AA_SIZE(b)-1]);
	AAElem_t* firstB = b->elems;
#endif

	ProductHeapChain_AA* maxElem = NULL;
	ProductHeapChain_AA* nextMaxElem = NULL;
	degrees_t* nextDegs;
	while ( (nextDegs = prodheapPeek_AA(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		cElems[k].degs = *nextDegs;
		mpq_init(cElems[k].coef);

		while (nextDegs != NULL && isEqualExponentVectors(cElems[k].degs, *nextDegs)) {
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AA* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AA(h);
			while (maxElem != NULL) {
#if SMQP_INT_PRODHEAP
				mpq_mul(ccoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
#else
				mpq_mul(ccoef, maxElem->a_i->coef, maxElem->b->coef);
#endif
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
#if SMQP_INT_PRODHEAP
			prodheapInsert_AA(h, maxElem, aElems[maxElem->a_i].degs + bElems[maxElem->b].degs);
#else
			prodheapInsert_AA(h, maxElem);
#endif
			maxElem = nextMaxElem;
		}
	}

	mpq_clear(ccoef);
	prodheapFree_AA(h);

	AA_SIZE(c) = k;

	resizePolynomial_AA(c, k);

	return c;
}



Node* multiplyPolynomials_inp(Node* a, Node* b, int nvar) {
	//TODO actually do this in-place.

	Node* prod = multiplyPolynomials(a, b, nvar);
	freePolynomial(a);
	return prod;
}

AltArr_t* multiplyPolynomials_AA_inp(AltArr_t* a, AltArr_t* b, int nvar) {
	AltArr_t* prod = multiplyPolynomials_AA(a,b,nvar);
	// fprintf(stderr, "multA:= ");
	// for (int idx = 0; idx < a->size; ++idx) {
	// 	gmp_fprintf(stderr, "(%Qd)*y^%lld*z^%lld", a->elems[idx].coef, GET_SECOND_EXP(a->elems[idx].degs), GET_THIRD_EXP(a->elems[idx].degs));
	// 	if (idx < a->size - 1) {
	// 		fprintf(stderr, " +");
	// 	}
	// }
	// fprintf(stderr, ";\nmultB:=");
	// for (int idx = 0; idx < b->size; ++idx) {
	// 	gmp_fprintf(stderr, "(%Qd)*y^%lld*z^%lld", b->elems[idx].coef, GET_SECOND_EXP(b->elems[idx].degs), GET_THIRD_EXP(b->elems[idx].degs));
	// 	if (idx < b->size - 1) {
	// 		fprintf(stderr, " +");
	// 	}
	// }
	// fprintf(stderr, ";\nmultC:=");
	// for (int idx = 0; idx < prod->size; ++idx) {
	// 	gmp_fprintf(stderr, "(%Qd)*y^%lld*z^%lld", prod->elems[idx].coef, GET_SECOND_EXP(prod->elems[idx].degs), GET_THIRD_EXP(prod->elems[idx].degs));
	// 	if (idx < prod->size - 1) {
	// 		fprintf(stderr, " +");
	// 	}
	// }
	// fprintf(stderr, ";\n\n\n");
	freePolynomial_AA(a);
	return prod;
}


void multiplyPolynomialAtIdxByXn_AA_inp (AltArr_t* aa, int idx, int n, int nvar)
{
	if (aa == NULL) {
		return;
	}

	if (idx < 0 || idx >= nvar) {
		return;
	}

	if (aa->unpacked) {
		// TODO
	    // multiplyPolynomialByXnAtIdx_AA_inp_unpk(aa, idx, degsList, nvar);
    	fprintf(stderr, "SMQP Error, multiplyPolynomialByXnAtIdx_AA_inp not yet implemented with exponent unpacking!\n");
    	exit (EXIT_FAILURE);
		// return;
	}

	const int* sizes = getExpOffsetArray (nvar);
	const degrees_t* maxDegs = getMaxExpArray(nvar);

	if (n > maxDegs[idx]) {
    	fprintf(stderr, "SMQP Error, multiplyPolynomialByXnAtIdx_AA_inp doesn't work when n > maxDegs[idx] !\n");
    	exit (EXIT_FAILURE);
	}

	for (int i = 0; i < AA_SIZE(aa); i++) {

		aa->elems[i].degs += ((degrees_t) n << sizes[idx]);
	}

}

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
Node* exponentiatePoly(Node* a, unsigned int n, int nvar) {

	if (n == 0) {
		// degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
		degrees_t degs = 0;
		ratNum_t coef;
		mpq_init(coef);
		mpq_set_ui(coef, 1ul, 1ul);
		Node* ret = addTerm(NULL, degs, coef);
		mpq_clear(coef);
		return ret;
	} else if (n == 1) {
		return deepCopyPolynomial(a, nvar);
	}

	Node* r = NULL;
	Node* b = a;
	while (n > 1) {
		if (n & 1) {
			r = (r == NULL) ? b : multiplyPolynomials(r,b,nvar);
		}
		b = multiplyPolynomials(b,b,nvar);
		n >>= 1;
	}
	r = (r == NULL) ? b : multiplyPolynomials(r,b,nvar);

	return r;
}

AltArr_t* exponentiatePoly_AA(AltArr_t* a, unsigned int n, int nvar) {

	// NOTE! 11/30/18
	//As long as this function simply calls multiplyPolynomials then
	//we don't need to worry about handling packed/unpacked exponents.
	//If this changes then one must create exponentiatePoly_AA_unpk.

	if (n == 0) {
		AltArr_t* ret = makePolynomial_AA(1, nvar);
		mpq_init(ret->elems->coef);
		mpq_set_ui(ret->elems->coef, 1ul, 1ul);
		ret->elems->degs = 0;
		ret->size = 1;
		return ret;
	} else if (n == 1) {
		return deepCopyPolynomial_AA(a);
	}

	AltArr_t* r = NULL;
	AltArr_t* b = deepCopyPolynomial_AA(a);
	while (n > 1) {
		if (n & 1) {
			r = (r == NULL) ? deepCopyPolynomial_AA(b) : multiplyPolynomials_AA_inp(r, b, nvar);

		}
		b = multiplyPolynomials_AA_inp(b, b, nvar);

		n >>= 1;
	}
	r = (r == NULL) ? deepCopyPolynomial_AA(b) : multiplyPolynomials_AA_inp(r, b, nvar);

	freePolynomial_AA(b);
	return r;
}



/*****************
 * Polynomial division
 *****************/

/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients
 * nvar: number of variables of monomials a and b
 */
// int monomialDivideTest(degrees_t adegs, degrees_t bdegs, int nvar) {
// 	return ( (adegs & FIRST_EXP) < (bdegs & SECOND_EXP)) ||
// 		   ( (adegs & SECOND_EXP) < (bdegs & SECOND_EXP)) ||
// 		   ( (adegs & THIRD_EXP) < (bdegs & THIRD_EXP));
// 	// for(int i = 0; i < nvar; ++i) {
// 	// 	if (adegs[i] < bdegs[i]) {
// 	// 		return 0;
// 	// 	}
// 	// }
// 	// return 1;
// }

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
Node* divisionGetNextTerm(ProductHeap* h) {
	productHeapElem** elems = h->elements;
	int size = h->heapSize;

	if (size == 0) {
		return NULL;
	}

	// degrees_t maxDegs = (degrees_t) malloc(sizeof(degree_t) * nvar);
	// memcpy(maxDegs, elems[0]->product->degs, sizeof(degree_t)*nvar);
	degrees_t maxDegs = elems[0]->product->degs;
	Node* ret = addTerm(NULL, maxDegs, NULL);
	productHeapElem* insertChain = NULL;
	productHeapElem* nextMaxElem;

	while(size > 0 && compareExponentVectors(elems[0]->product->degs, maxDegs) == 0) {
		productHeapElem* maxElem = prodheapRemoveMax(h);
		--size;

		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;
			mpq_add(ret->coef, ret->coef, maxElem->product->coef);
			maxElem->b = (maxElem->b)->next;
			if (maxElem->b != NULL) {
				addExponentVectors(maxElem->a_i->degs, maxElem->b->degs, maxElem->product->degs);
				mpq_mul(maxElem->product->coef, maxElem->a_i->coef, maxElem->b->coef);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				prodheapFreeElement(maxElem);
			}
			maxElem = nextMaxElem;
		}
	}

	while(insertChain != NULL) {
		nextMaxElem = insertChain->next;
		insertChain->next = NULL;
		prodheapInsert(h,insertChain);
		insertChain = nextMaxElem;
	}

	if (mpq_sgn(ret->coef) == 0) {
		freeNode(ret);
		ret = NULL;
	}

	return ret;
}

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
#if SMQP_INT_PRODHEAP
void divisionGetNextTerm_AA(ProductHeap_AA* h, const AAElem_t* __restrict__ aElems, const AAElem_t* __restrict__ bElems, mpq_t* retCoef) {
#else
void divisionGetNextTerm_AA(ProductHeap_AA* h, mpq_t* retCoef) {
#endif
	if (h->heapSize == 0) {
		return;
	}

#if SMQP_INT_PRODHEAP
	int lastB = h->lastB;
#else
	AAElem_t* lastB = h->lastB;
#endif

	ProductHeapChain_AA* insertChain = NULL;
	ProductHeapChain_AA* maxElem, *nextMaxElem;

	mpq_t prodCoef;
	mpq_init(prodCoef);
	degrees_t* nextDegs = prodheapPeek_AA(h);

	register degrees_t maxDegs = *nextDegs;

	while ( nextDegs != NULL && isEqualExponentVectors(maxDegs, *nextDegs)) {
		maxElem = prodheapRemoveMax_AA(h);

		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;
#if SMQP_INT_PRODHEAP
			mpq_mul(prodCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
#else
			mpq_mul(prodCoef, maxElem->a_i->coef, maxElem->b->coef);
#endif
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

	while(insertChain != NULL) {
		maxElem = insertChain->next;
		insertChain->next = NULL;
#if SMQP_INT_PRODHEAP
		prodheapInsert_AA(h,insertChain, aElems[insertChain->a_i].degs + bElems[insertChain->b].degs);
#else
		prodheapInsert_AA(h,insertChain);
#endif
		insertChain = maxElem;
	}

	mpq_clear(prodCoef);
}



/**
 * Given a polynomial, c, and a term, b, determine polynomials a and r
 * such that c = b*a + r.
 * a and r are returned in res_a and res_r, respectively.
 */
void divideBySingleTerm(Node* c, Node* b, Node** res_a, Node** res_r, int nvar) {
	if (b == NULL) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	Node* curC = c, *a = NULL, *r = NULL;
	degrees_t degs;
	ratNum_t coef;
	mpq_init(coef);

	*res_a = NULL;
	*res_r = NULL;

	while (curC != NULL) {
		if (monomialDivideTest(curC->degs,b->degs,nvar)) {
			mpq_div(coef, curC->coef, b->coef);
			// coef = curC->coef / b->coef;
			if (mpq_sgn(coef) != 0) {
				// degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
				subtractExponentVectors(curC->degs, b->degs, degs);
				a = addTerm(a,degs,coef);
				if (*res_a == NULL) {
					*res_a = a;
				}
			}
		} else {
			// degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			// memcpy(degs, curC->degs, sizeof(degree_t)*nvar);
			degs = curC->degs;
			mpq_set(coef, curC->coef);
			// coef = curC->coef;
			r = addTerm(r,degs,coef);
			if (*res_r == NULL) {
				*res_r = r;
			}
		}
		curC = curC->next;
	}
	mpq_clear(coef);
}

void divideBySingleTerm_AA(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int nvar) {
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
		divideBySingleTerm_AA_unpk(c, b, res_a, res_r, nvar);
		return;
	}


	AAElem_t* k = c->elems;
	AAElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0;

	AltArr_t* a = makePolynomial_AA(maxSize, nvar);
	AltArr_t* r = makePolynomial_AA(maxSize, nvar);
	AAElem_t* bElems = b->elems;
	AAElem_t* curA = a->elems;
	AAElem_t* curR = r->elems;
	mpq_init(curA->coef);
	mpq_init(curR->coef);

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar)) {
			mpq_div(curA->coef, k->coef, bElems->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {
			mpq_set(curR->coef, k->coef);
			curR->degs = k->degs;
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

	*res_a = a;
	*res_r = r;
}

/**
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials(Node* c, Node* b, Node** res_a, Node** res_r, int nvar) {

	if (b == NULL) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	// b does not divide c
	if(!monomialDivideTest(c->degs,b->degs,nvar)) {
		*res_a = NULL;
		*res_r = deepCopyPolynomial(c, nvar);
		return;
	}

	// b is a monomial so we can do a simple divide
	if (b->next == NULL) {
		divideBySingleTerm(c, b, res_a, res_r, nvar);
		return;
	}

	//init a with lt(c)/lt(b);
	ratNum_t bcoef;
	mpq_init(bcoef);
	mpq_set(bcoef, b->coef);
	// ratNum_t bcoef = b->coef;
	ratNum_t coef;
	mpq_init(coef);
	mpq_div(coef, c->coef, bcoef);
	degrees_t beta = b->degs;
	// degrees_t eps = (degrees_t) malloc(sizeof(degree_t)*nvar);
	degrees_t eps;
	subtractExponentVectors(c->degs, beta, eps);
	Node* a = addTerm(NULL, eps, coef);
	Node* r = NULL;
	*res_a = a; //set head of resulting a to this leading term of a;
	*res_r = r;
	Node* k = c->next;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap* h = prodheapCreate(nvar);
	prodheapInsert(h, prodheapMakeElement(h, a, b->next));

	//loop variables
	// degrees_t delta = NULL;
	// degrees_t adegs = NULL;
	degrees_t* delta;
	degrees_t adegs;
	Node* multerm = NULL;
	cmpExp_t cmp;
	while(k != NULL || h->heapSize > 0) {
		ratNum_t e;
		mpq_init(e);

		delta = prodheapPeek(h);
		if (k == NULL) {
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
			multerm = divisionGetNextTerm(h);
			if (multerm == NULL) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				continue;
			} else {
				eps = multerm->degs;
				mpq_neg(e, multerm->coef);
				// e = mpq_mul-1*multerm->coef;
			}
		} else if (cmp == 0) {
			multerm = divisionGetNextTerm(h);
			if (multerm == NULL) {
				continue;
			} else {
				eps = multerm->degs;
				mpq_sub(e, k->coef, multerm->coef);
				// e = k->coef - multerm->coef;
				k = k->next;
			}
		} else {
			multerm = NULL;
			eps = k->degs;
			mpq_set(e,k->coef);
			// e = k->coef;
			k = k->next;
		}

		if (mpq_sgn(e) != 0) {
			if (monomialDivideTest(eps, b->degs, nvar)) {
				// adegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
				subtractExponentVectors(eps, beta, adegs);
				mpq_div(e, e, bcoef);
				// e /= bcoef;
				a = addTerm(a, adegs, e);
				prodheapInsert(h, prodheapMakeElement(h, a, b->next));
				freeNode(multerm); //node and its degs
			} else {
				if (multerm != NULL) {
					r = addTerm(r,eps,e);
					free(multerm); //just the node, r holds multerm->degs
				} else {
					//else, eps is equal to k->degs, so make a copy.
					// degrees_t temp = (degrees_t) malloc(sizeof(degree_t)*nvar);
					// memcpy(temp, eps, sizeof(degree_t)*nvar);
					r = addTerm(r,eps,e);
				}
				if (*res_r == NULL) {
					*res_r = r;
				}
			}
		} else {
			freeNode(multerm);
		}
		mpq_clear(e);
	}
}



#define MONOMIAL_DIV_PTR 1

/**
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials_AA(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, register int nvar) {

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

	if (c->unpacked || b->unpacked) {
		dividePolynomials_AA_unpk(c, b, res_a, res_r, nvar);
		return;
	}

	AAElem_t* __restrict__ k = c->elems;
	AAElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0;

	AltArr_t* a = makePolynomial_AA(maxASize, nvar);
	AltArr_t* r = makePolynomial_AA(maxRSize, nvar);
	AAElem_t* __restrict__ curA = a->elems;
	AAElem_t* __restrict__ curR = r->elems;
	mpq_init(curA->coef);
	mpq_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
#if MONOMIAL_DIV_PTR
	DivTest_ptr divTest = getMonomialDivideTestFuncPtr(nvar);
	while (k != lenK && !(*divTest)(k->degs, beta)) {
#else
	while (k != lenK && !monomialDivideTest(k->degs, beta, nvar)) {
#endif
		mpq_set(curR->coef, k->coef);
		curR->degs = k->degs;
		++j;
		if (j >= maxRSize) {
			maxRSize <<= 1;
			r->elems = (AAElem_t*) realloc(r->elems, maxRSize*sizeof(AAElem_t));
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

	subtractExponentVectors(k->degs, beta, curA->degs);
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	prodheapResize_AA(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
#if SMQP_INT_PRODHEAP
	prodheapInsert_AA(h, prodheapMakeChain_AA(0, 1, NULL), curA->degs + b2Elem->degs);
#else
	prodheapInsert_AA(h, prodheapMakeChain_AA(curA, b2Elem, NULL));
#endif
	++i;
	++curA;
	mpq_init(curA->coef);

	degrees_t* delta = prodheapPeek_AA(h);
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
#if SMQP_INT_PRODHEAP
			divisionGetNextTerm_AA(h, a->elems, b->elems, &(curA->coef));
#else
			divisionGetNextTerm_AA(h, &(curA->coef));
#endif
			if (mpq_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AA(h);
				continue;
			} else {
				mpq_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
#if SMQP_INT_PRODHEAP
			divisionGetNextTerm_AA(h, a->elems, b->elems, &(curA->coef));
#else
			divisionGetNextTerm_AA(h, &(curA->coef));
#endif
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
			eps = k->degs;
			mpq_set(curA->coef,k->coef);
			++k;
		}

#if MONOMIAL_DIV_PTR
		if ((*divTest)(eps, beta)) {
#else
		if (monomialDivideTest(eps, beta, nvar)) {
#endif
			subtractExponentVectors(eps, beta, curA->degs);
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAElem_t*) realloc(a->elems, maxASize*sizeof(AAElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AA(h, maxASize);
			}
#if SMQP_INT_PRODHEAP
			prodheapInsert_AA(h, prodheapMakeChain_AA(i, 1, NULL), curA->degs + b2Elem->degs);
#else
			prodheapInsert_AA(h, prodheapMakeChain_AA(curA, b2Elem, NULL));
#endif
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {

			//swap here so that curA becomes 0.
			mpq_swap(curR->coef, curA->coef);
			curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAElem_t*) realloc(r->elems, maxRSize*sizeof(AAElem_t));
				curR = r->elems + j - 1;
			}
			++(curR);
			mpq_init(curR->coef);
		}

		delta = prodheapPeek_AA(h);
	}

	prodheapFree_AA(h);

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);
	mpq_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	*res_a = a;
	*res_r = r;
}

void exactDividePolynomials_AA (AltArr_t* c, AltArr_t* b, AltArr_t** res_a, register int nvar) {
    if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf (stderr, "Division by zero! Exiting...");
		exit (EXIT_FAILURE);
    }

    if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AA (0, nvar);
		return;
    }

	if (c->unpacked || b->unpacked) {
		exactDividePolynomials_AA_unpk(c, b, res_a, nvar);
		return;
	}

    // b is a monomial so we can do a simple divide
    if (AA_SIZE(b) == 1) {
	AltArr_t* res_r;
	divideBySingleTerm_AA (c, b, res_a, &res_r, nvar);
	//fprintf (stderr, "r := \n" );
	//printAAZ (res_r);
	freePolynomial_AA (res_r);
	return;
    }

    AAElem_t* __restrict__ k = c->elems;
    AAElem_t* __restrict__ lenK = k + AA_SIZE(c);
    AAElem_t* __restrict__ b2Elem = b->elems + 1;
    register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
    register int i = 0;

    AltArr_t* a = makePolynomial_AA (maxASize, nvar);
    // AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
    AAElem_t* __restrict__ curA = a->elems;
    // AAZElem_t* __restrict__ curR = r->elems;
    mpq_init(curA->coef);
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
	mpq_clear (curA->coef);
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
    /* mpq_divexact (curA->coef, k->coef, b->elems->coef); */
    mpq_div(curA->coef, k->coef, b->elems->coef);
    ++k;

    //init multiplication between a (quotient) and b (divisor)
    ProductHeap_AA* h = prodheapCreate_AA (nvar);
    prodheapResize_AA (h, maxASize);
    h->lastB = AA_SIZE(b) - 1;
    prodheapInsert_AA (h, prodheapMakeChain_AA (0, 1, NULL), curA->degs + b2Elem->degs);
    ++i;
    ++curA;
    mpq_init (curA->coef);

    degrees_t* delta = prodheapPeek_AA (h);
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
		    divisionGetNextTerm_AA (h, a->elems, b->elems, &(curA->coef));
		    if (mpq_sgn (curA->coef) == 0) {
			//in this case, the term with degree delta ended up
			//having its coffeicient cancelled out (i.e. 0)
			delta = prodheapPeek_AA (h);
			continue;
		    } else {
			mpq_neg (curA->coef, curA->coef);
		    }
		} else if (cmp == 0) {
		    eps = *delta;
		    divisionGetNextTerm_AA (h, a->elems, b->elems, &(curA->coef));
		    if (mpq_sgn (curA->coef) == 0) {
			delta = prodheapPeek_AA (h);
			continue; //the chains cancelled themselves out since the peek
		    } else {
			mpq_sub (curA->coef, k->coef, curA->coef);
			++k;
			if (mpq_sgn (curA->coef) == 0) {
			    delta = prodheapPeek_AA (h);
			    continue;
			}
		    }
		} else {
		    eps = k->degs;
		    mpq_set (curA->coef,k->coef);
		    ++k;
		}

		if ((*divTest)(eps, beta)) {
		    subtractExponentVectors(eps, beta, curA->degs);
		    /* mpq_divexact(curA->coef, curA->coef, b->elems->coef); */
		    mpq_div(curA->coef, curA->coef, b->elems->coef);
		    if (i+1 >= maxASize) {
			maxASize <<= 1;
			a->elems = (AAElem_t*) realloc(a->elems, maxASize*sizeof(AAElem_t));
			curA = a->elems + i;
			//prodheap maximum size should be equal to the size of a
			prodheapResize_AA(h, maxASize);
		    }
		    prodheapInsert_AA(h, prodheapMakeChain_AA(i, 1, NULL), curA->degs + b2Elem->degs);
		    ++i;
		    ++(curA);
		    mpq_init(curA->coef);
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

		delta = prodheapPeek_AA(h);
    }

    prodheapFree_AA(h);

    //clear since we always setup one past where we actually are.
    mpq_clear(curA->coef);
    //mpz_clear(curR->coef);

    AA_SIZE(a) = i;
    //AA_SIZE(r) = j;
    a->alloc = maxASize;
    //r->alloc = maxRSize;

    *res_a = a;
    //*res_r = r;
    return;
}

void multiplyByRational_AA_inp(AltArr_t* aa, const mpq_t z) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	for (int i = 0; i < aa->size; ++i) {
		mpq_mul(aa->elems[i].coef, aa->elems[i].coef, z);
	}
}

void univariatePseudoDividePolynomials_AA(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, int lazy) {

	if (b == NULL || b->size == 0) {
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || c->size == 0) {
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (c->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "[BPAS] Error: calling univar pseudo divide with multivariate polynomials.\n");
		exit(EXIT_FAILURE);
	}


	int packedC = 0;
	int packedB = 0;
	if (c->unpacked) {
		packExponentVectors_AA_inp(c);
		packedC = 1;
	}
	if (b->unpacked) {
		packExponentVectors_AA_inp(b);
		packedB = 1;
	}


	if (b->elems->degs > c->elems->degs) {
		*res_r = deepCopyPolynomial_AA(c);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (packedC) {
			unpackExponentVectors_AA_inp(c);
		}
		if (packedB) {
			unpackExponentVectors_AA_inp(b);
		}
		return;
	}

	//since we're in a field, pseudo division is just division.
	dividePolynomials_AA(c, b, res_a, res_r, 1);

	int steps = 0;
	if (!lazy) {
		degrees_t d = c->elems->degs - b->elems->degs + 1;
		steps = d;

		mpq_t bPow;
		mpq_init(bPow);
		mpz_pow_ui(mpq_numref(bPow), mpq_numref(b->elems->coef), d);
		mpz_pow_ui(mpq_denref(bPow), mpq_denref(b->elems->coef), d);
		mpq_canonicalize(bPow);
		multiplyByRational_AA_inp(*res_a, bPow);
		multiplyByRational_AA_inp(*res_r, bPow);
		mpq_clear(bPow);

		// for (int j = 0; j < d; ++j) {
			// multiplyByRational_AA_inp(*res_a, b->elems->coef);
			// multiplyByRational_AA_inp(*res_r, b->elems->coef);
		// }
	}

	if (e != NULL) {
		*e = steps;
	}

	if (packedC) {
		unpackExponentVectors_AA_inp(c);
	}
	if (packedB) {
		unpackExponentVectors_AA_inp(b);
	}

	return;
}

int divideTestSingleTerm_AA(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, int nvar) {
	if (b == NULL || b->size == 0) {
		return 0;
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		return 1;
	}

	if (c->unpacked || b->unpacked) {
		return divideTestSingleTerm_AA_unpk(c, b, res_a, nvar);
	}

	AAElem_t* k = c->elems;
	AAElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;

	AltArr_t* a = makePolynomial_AA(maxSize, nvar);
	AAElem_t* bElems = b->elems;
	AAElem_t* curA = a->elems;
	mpq_init(curA->coef);

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar)) {
			mpq_div(curA->coef, k->coef, bElems->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {
			a->size = i;
			freePolynomial_AA(a);
			return 0;
		}
		++k;
	}

	mpq_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxSize;

	if (res_a != NULL) {
		*res_a = a;
	} else {
		freePolynomial_AA(a);
	}

	return 1;
}

int divideTest_AA(AltArr_t* c, AltArr_t* b, AltArr_t** res_a, int nvar) {

	if (b == NULL || AA_SIZE(b) == 0) {
		return 0;
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		*res_a = makePolynomial_AA(0, nvar);
		return 1;
	}

	if (c->unpacked || b->unpacked) {
		return divideTest_AA_unpk(c, b, res_a, nvar);
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		return divideTestSingleTerm_AA(c, b, res_a, nvar);
	}

	AAElem_t* __restrict__ k = c->elems;
	AAElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAElem_t* __restrict__ b2Elem = b->elems + 1;
	register degrees_t beta = b->elems->degs;

	if (!monomialDivideTest(k->degs, beta, nvar)) {
		return 0;
	}

	register int maxASize = AA_SIZE(c) + 1;
	register int i = 0;
	AltArr_t* a = makePolynomial_AA(maxASize, nvar);
	AAElem_t* __restrict__ curA = a->elems;
	mpq_init(curA->coef);

	//init a with lt(c)/lt(b);
	subtractExponentVectors(k->degs, beta, curA->degs);
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AA* h = prodheapCreate_AA(nvar);
	prodheapResize_AA(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
#if SMQP_INT_PRODHEAP
	prodheapInsert_AA(h, prodheapMakeChain_AA(0, 1, NULL), curA->degs + b2Elem->degs);
#else
	prodheapInsert_AA(h, prodheapMakeChain_AA(curA, b2Elem, NULL));
#endif
	++i;
	++curA;
	mpq_init(curA->coef);

	degrees_t* delta = prodheapPeek_AA(h);
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
#if SMQP_INT_PRODHEAP
			divisionGetNextTerm_AA(h, a->elems, b->elems, &(curA->coef));
#else
			divisionGetNextTerm_AA(h, &(curA->coef));
#endif
			if (mpq_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AA(h);
				continue;
			} else {
				mpq_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
#if SMQP_INT_PRODHEAP
			divisionGetNextTerm_AA(h, a->elems, b->elems, &(curA->coef));
#else
			divisionGetNextTerm_AA(h, &(curA->coef));
#endif
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
			eps = k->degs;
			mpq_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest(eps, beta, nvar)) {
			subtractExponentVectors(eps, beta, curA->degs);
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAElem_t*) realloc(a->elems, maxASize*sizeof(AAElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AA(h, maxASize);
			}
#if SMQP_INT_PRODHEAP
			prodheapInsert_AA(h, prodheapMakeChain_AA(i, 1, NULL), curA->degs + b2Elem->degs);
#else
			prodheapInsert_AA(h, prodheapMakeChain_AA(curA, b2Elem, NULL));
#endif
			++i;
			++(curA);
			mpq_init(curA->coef);
		} else {

			//divide test fails
			prodheapFree_AA(h);
			a->size = i;
			freePolynomial_AA(a);
			return 0;
		}

		delta = prodheapPeek_AA(h);
	}

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxASize;

	if (res_a != NULL) {
		*res_a = a;
	} else {
		freePolynomial_AA(a);
	}
	return 1;
}


void divideByLeadingTerms_AA (AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int nvar)
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
		divideByLeadingTerms_AA_unpk (c, b, res_a, res_r, nvar);
 		return;
	}

	AltArr_t* a;
	AltArr_t* r;

	/* while (k != lenK) { */
	if (monomialDivideTest (c->elems->degs, b->elems->degs, nvar)) {
		// a:
		a = makePolynomial_AA (1, nvar);
		mpq_init (a->elems->coef);
		mpq_div (a->elems->coef, c->elems->coef, b->elems->coef);
		subtractExponentVectors (c->elems->degs, b->elems->degs, a->elems->degs);
		AA_SIZE (a) = 1;

		// r:
		r = NULL;
	} else {
		// a:
		a = makePolynomial_AA (1, nvar);
		mpq_init (a->elems->coef);
		mpq_set_si (a->elems->coef, 1l, 1lu);
		a->elems->degs = 0l;
		AA_SIZE (a) = 1;

		// r:
		r = makePolynomial_AA (1, nvar);
		mpq_init (r->elems->coef);
		mpq_set(r->elems->coef, c->elems->coef);
		r->elems->degs = c->elems->degs;
		AA_SIZE (r) = 1;
	}

	*res_a = a;
	*res_r = r;
}

void applyModuloSymmetricPrimPart_AA_inp(AltArr_t* p, const mpz_t mod) {
	if (isZero_AA(p)) {
		return;
	}

	primitivePart_AA_inp(p);

	polysize_t newSize = 0;
	mpz_t halfMod;
	mpz_init_set(halfMod, mod);
	mpz_div_ui(halfMod, halfMod, 2ul);

	int unpacked = p->unpacked;
	int nvar = p->nvar;
	polysize_t i;
	polysize_t insertIdx = 0;
	polysize_t size = p->size;
	AAElem_t* elems = p->elems;
	for (i = 0; i < size; ++i) {
		// fprintf(stderr, "moding poly %p at %ld\n",p, i);
		mpz_mod(mpq_numref(elems[insertIdx].coef), mpq_numref(elems[i].coef), mod);
		if(mpz_cmp(mpq_numref(elems[insertIdx].coef), halfMod) > 0) {
			mpz_sub(mpq_numref(elems[insertIdx].coef), mpq_numref(elems[insertIdx].coef), mod);
		}

		if (unpacked) {
			 setExponentVector_unpk(elems[insertIdx].degs, elems[i].degs, nvar);
		} else {
			elems[insertIdx].degs = elems[i].degs;
		}

		if (mpz_sgn(mpq_numref(elems[insertIdx].coef)) != 0) {
			++insertIdx;
		}
	}
	for (i = insertIdx; i < size; ++i) {
		mpq_clear(elems[i].coef);
	}

	p->size = insertIdx;
	mpz_clear(halfMod);
}


/*****************
 * Derivative / Integral
 *****************/

AltArr_t* derivative_AA(AltArr_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx < 0 || idx >= aa->nvar) {
		return NULL;
	}

	if (aa->unpacked) {
		return derivative_AA_unpk(aa, idx, k);
	}

	AltArr_t* ret = makePolynomial_AA(aa->size, aa->nvar);

	const int* sizes = getExpOffsetArray(aa->nvar);
	const degrees_t* masks = getExpMaskArray(aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAElem_t* __restrict__ elems = aa->elems;
	AAElem_t* __restrict__ retElems = ret->elems;
	degrees_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
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
		mpq_init(retElems[insertIdx].coef);
		mpq_set(retElems[insertIdx].coef, elems[i].coef);

		mpq_set_ui(mpqDeg, deg, 1ul);
		for (int j = 0; j < k; ++j) {
			mpq_mul(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);
			--deg;
			mpz_sub(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}
		retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	return ret;
}

/**
 * Integrate with respect to a variable that does not currently exist in aa.
 * It becomes the main variable (that is, to the left).
 *
 */
AltArr_t* integrateExpand_AA(AltArr_t* aa, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return integrateExpand_AA_unpk(aa, k);
	}

	const degrees_t* maxs = getMaxExpArray((aa->nvar)+1);
	degrees_t degsK = k;
	if (degsK > maxs[0]) {
		unpackExponentVectors_AA_inp(aa);
		AltArr_t* ret = integrateExpand_AA_unpk(aa, k);
		packExponentVectors_AA_inp(aa);
		return ret;
	}

	AltArr_t* ret = deepCopyPolynomial_AA(aa);
	expandNumVarsLeft_AA(ret, (aa->nvar)+1);
	if (ret->unpacked) {
		//expanding might unpack ret;
		integrateExpand_AA_inp_unpk(ret, k);
		return ret;
	}

	int size = aa->size;
	AAElem_t* elems = ret->elems;

	mpq_t kFact;
	mpq_init(kFact);
	mpz_fac_ui(mpq_numref(kFact), (unsigned long)k);

	int expOffset = getMVarExpOffset(ret->nvar);

	for (int i = 0; i < size; ++i) {
		mpq_div(elems[i].coef, elems[i].coef, kFact);
		elems[i].degs |= (degsK << expOffset);
	}

	mpq_clear(kFact);

	return ret;
}

AltArr_t* integral_AA(AltArr_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx >= aa->nvar || idx < 0) {
		return integrateExpand_AA(aa, k);
	}

	if (aa->unpacked) {
		return integral_AA_unpk(aa, idx, k);
	}

	AltArr_t* ret = makePolynomial_AA(aa->size, aa->nvar);

	const int* sizes = getExpOffsetArray(aa->nvar);
	const degrees_t* masks = getExpMaskArray(aa->nvar);
	const degrees_t* maxs = getMaxExpArray(aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAElem_t* __restrict__ elems = aa->elems;
	AAElem_t* __restrict__ retElems = ret->elems;
	degrees_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i) {
		deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);

		retElems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpq_init(retElems[insertIdx].coef);
		mpq_set(retElems[insertIdx].coef, elems[i].coef);

		mpq_set_ui(mpqDeg, deg+1, 1ul);
		for (int j = 0; j < k; ++j) {
			mpq_div(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);
			++deg;
			mpz_add(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}

		if (deg > maxs[idx]) {
			//then stop what we're doing and integrate with unpacked;
			// fprintf(stderr, "SMQP exponent overflow in integration! Index: %d, value %llu\n", idx, deg);
			// exit(1);
			// TODO!!
			ret->size = insertIdx;
			freePolynomial_AA(ret);
			unpackExponentVectors_AA_inp(aa);
			ret = integral_AA_unpk(aa, idx, k);
			packExponentVectors_AA_inp(aa);
			return ret;
		}

		retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	return ret;
}



/*****************
 * Content, PrimitivePart, etc.
 *****************/

void integralContent_AA(AltArr_t* aa, mpq_t ret) {
	if (aa == NULL || aa->size <= 0) {
		mpq_set_ui(ret, 0ul, 1ul);
		return;
	}

	AAElem_t* elems = aa->elems;
	int size = aa->size;
	mpq_set_ui(ret, 1ul, 1ul);

	mpq_t one;
	mpq_init(one);
	mpq_set_si(one, 1l, 1l);

	mpq_t cont;
	mpq_init(cont);
	mpq_abs(cont, elems->coef);

	mpq_t tempa;
	mpq_t tempb;
	mpq_init(tempa);
	mpq_init(tempb);

	for (int i = 1; i < size; ++i) {
		// if (mpz_cmp(mpq_numref(one), mpq_numref(cont)) == 0) {
		// 	mpq_clear(cont);
		// 	mpq_clear(one);
		// 	mpq_clear(tempa);
		// 	mpq_clear(tempb);
		// 	fprintf(stderr, "gcd is 1\n");
		// 	return;
		// }

		mpq_set(tempa, cont);
		mpq_set(tempb, elems[i].coef);

		mpz_mul(mpq_numref(tempa), mpq_numref(tempa), mpq_denref(tempb));
		mpz_mul(mpq_numref(tempb), mpq_numref(tempb), mpq_denref(tempa));
		mpz_mul(mpq_denref(tempa), mpq_denref(tempa), mpq_denref(tempb));

		mpz_gcd(mpq_numref(cont), mpq_numref(tempa), mpq_numref(tempb));
		mpz_set(mpq_denref(cont), mpq_denref(tempa));
		mpq_canonicalize(cont);
	}

	if (mpq_sgn(elems->coef) < 0 && mpq_sgn(cont) > 0) {
		mpq_neg(cont, cont);
	}

	mpq_set(ret, cont);
	mpq_clear(cont);
	mpq_clear(one);
	mpq_clear(tempa);
	mpq_clear(tempb);
}

AltArr_t* primitivePart_AA(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	mpq_t content;
	mpq_init(content);

	integralContent_AA(aa, content);

	AltArr_t* res = deepCopyPolynomial_AA(aa);
	AAElem_t* elems = res->elems;
	int size = res->size;
	for (int i = 0; i < size; ++i) {
		mpq_div(elems[i].coef, elems[i].coef, content);
	}

	mpq_clear(content);

	return res;
}


AltArr_t* primitivePartAndContent_AA(AltArr_t* aa, mpq_t cont) {
	if (aa == NULL || aa->size == 0) {
		mpq_set_si(cont, 0l, 1l);
		return NULL;
	}

	integralContent_AA(aa, cont);

	AltArr_t* res = deepCopyPolynomial_AA(aa);
	AAElem_t* elems = res->elems;
	int size = res->size;
	for (int i = 0; i < size; ++i) {
		mpq_div(elems[i].coef, elems[i].coef, cont);
	}

	return res;
}

void primitivePart_AA_inp(AltArr_t* aa) {
	if (aa == NULL || aa->size <= 0) {
		return;
	}

	mpq_t content;
	mpq_init(content);

	integralContent_AA(aa, content);

	AAElem_t* elems = aa->elems;
	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpq_div(elems[i].coef, elems[i].coef, content);
	}

	mpq_clear(content);
}

AltArr_t* univariateGCD_AA(AltArr_t* a, AltArr_t* b) {
	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AA(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AA(a);
	}

	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMQP ERROR: Calling univariate GCD on multivariate arrays\n");
		exit(1);
	}

	//conversion to DUZP handles unpacking. Return is always packed.
	mpq_t conta;
	mpq_t contb;
	mpq_init(conta);
	mpq_init(contb);
	AltArr_t* primA = primitivePartAndContent_AA(a, conta);
	AltArr_t* primB = primitivePartAndContent_AA(b, contb);

	DUZP_t* aDUZP = convertFromAltArr_DUZP(primA);
	DUZP_t* bDUZP = convertFromAltArr_DUZP(primB);
	DUZP_t* gDUZP = primitiveGCD_DUZP(aDUZP, bDUZP);

	mpz_mul(mpq_numref(conta), mpq_numref(conta), mpq_denref(contb));
	mpz_mul(mpq_numref(contb), mpq_numref(contb), mpq_denref(conta));
	mpz_mul(mpq_denref(conta), mpq_denref(conta), mpq_denref(contb));
	mpz_gcd(mpq_numref(conta), mpq_numref(conta), mpq_numref(contb));
	mpq_canonicalize(conta);
	AltArr_t* gAA = convertToAltArr_DUZP(gDUZP);
	multiplyByRational_AA_inp(gAA, conta);

	freePolynomial_DUZP(aDUZP);
	freePolynomial_DUZP(bDUZP);
	freePolynomial_DUZP(gDUZP);
	freePolynomial_AA(primA);
	freePolynomial_AA(primB);
	mpq_clear(conta);
	mpq_clear(contb);

	return gAA;

	if (a->unpacked || b->unpacked) {
		return univariateGCD_AA_unpk(a, b);
	}

	AltArr_t* r0 = NULL;
	AltArr_t* r1 = NULL;
	AltArr_t* r2 = NULL;

	mpq_t c0;
	mpq_t c1;
	mpq_init(c0);
	mpq_init(c1);

	if (isGreaterExponentVectors(a->elems->degs, b->elems->degs)) {
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

	return r0;

}

int integerPolynomialTest_AA(AltArr_t* aa) {
	if (isZero_AA(aa)) {
		return 1;
	}
	if (isConstant_AA(aa)) {
		return (mpz_cmp_si(mpq_denref(aa->elems->coef), 1l) == 0);
	}

	for (int i = 0; i < aa->size; ++i) {
		if (mpz_cmp_si(mpq_denref(aa->elems[i].coef), 1l) != 0) {
			return 0;
		}
	}

	return 1;
}


void integerPolynomialTestCont_AA(AltArr_t* aa, mpz_t mpzG) {
	if (aa == NULL || aa->size == 0) {
		mpz_set_si(mpzG, 0l);
		return;
	}

	mpz_set(mpzG, mpq_numref(aa->elems->coef));
	for (int i = 0; i < aa->size; ++i) {
		if (mpz_cmp_si(mpq_denref(aa->elems[i].coef), 1l) != 0) {
			mpz_set_si(mpzG, 0l);
			return;
		}

		if (mpz_cmp_si(mpzG, 1l) != 0) {
			mpz_gcd(mpzG, mpzG, mpq_numref(aa->elems[i].coef));
		}
	}
}

AltArr_t* commonFactor_AA(AltArr_t* a, AltArr_t** factored) {
	if (a == NULL || a->size == 0) {
		return NULL;
	}

	if (a->unpacked) {
		return commonFactor_AA_unpk(a, factored);
	}

	AltArr_t* ret = makePolynomial_AA(1, a->nvar);
	mpq_init(ret->elems->coef);
	mpq_set_ui(ret->elems->coef, 1ul, 1ul);

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
		AltArr_t* factRet = deepCopyPolynomial_AA(a);
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
 * Interpolation / Evaluation
 *****************/

/**
 * Multiply a, univariate, by a univariate binomial (c*x + d) where c and d
 * are constants.
 */
AltArr_t* multiplyByBinomial_AA_inp(AltArr_t* a, AltArr_t* b) {
	int size = a->size;
	AltArr_t* temp = makePolynomial_AA(size, a->nvar);
	for (int i = 0; i < size; ++i) {
		mpq_init(temp->elems[i].coef);
		mpq_mul(temp->elems[i].coef, a->elems[i].coef, b->elems[1].coef);
		temp->elems[i].degs = a->elems[i].degs;

		mpq_mul(a->elems[i].coef, a->elems[i].coef, b->elems[0].coef);
		addExponentVectors(a->elems[i].degs, b->elems[0].degs, a->elems[i].degs);
	}
	temp->size = size;
	a = addPolynomials_AA_inp(a, temp, a->nvar);

	freePolynomial_AA(temp);
	return a;
}


AltArr_t* univarInterpolate_AA(mpq_t* points, mpq_t* vals, int nPoints) {
	mpq_t denoms[nPoints];
	mpq_t diffs[nPoints];

	// for (int i = 0; i < nPoints; ++i) {
	// 	gmp_fprintf(stderr, "points[%d]: %Qd\n", i, points[i]);
	// }
	// for (int i = 0; i < nPoints; ++i) {
	// 	gmp_fprintf(stderr, "vals[%d]: %Qd\n", i, vals[i]);
	// }

	for (int i = 0; i < nPoints; ++i) {
		mpq_init(denoms[i]);
		mpq_init(diffs[i]);
		mpq_set_ui(denoms[i], 1ul, 1ul);
	}

	register int j = 0;
	for (j = 1; j < nPoints; ++j) {
		mpq_sub(diffs[j], points[0], points[j]);
		mpq_mul(denoms[0], denoms[0], diffs[j]);
	}
	mpq_inv(denoms[0], denoms[0]);
	mpq_mul(denoms[0], denoms[0], vals[0]);
	for (int i = 1; i < nPoints; ++i) {
		for (j = 0; j < i-1; ++j) {
			mpq_sub(diffs[j+1], points[i], points[j]);
			mpq_mul(denoms[i], denoms[i], diffs[j+1]);
		}
		mpq_neg(diffs[j+1], diffs[j+1]);
		mpq_mul(denoms[i], denoms[i], diffs[j+1]);
		for (int j = i+1; j < nPoints; ++j) {
			mpq_sub(diffs[j], points[i], points[j]);
			mpq_mul(denoms[i], denoms[i], diffs[j]);
		}

		mpq_inv(denoms[i], denoms[i]);
		mpq_mul(denoms[i], denoms[i], vals[i]);
	}

	// mpq_set_si(diffs[0], 0l, 1l);
	AltArr_t* ret = NULL; //makeConstPolynomial_AA(1, 1, diffs[0]);
	mpq_set_si(diffs[0], 1l, 1l); //set to 1 for prod below


	//Make all the possible factors of numerators now.
	AltArr_t** q = (AltArr_t**) malloc(sizeof(AltArr_t*)*nPoints);
	for (int i = 0; i < nPoints; ++i) {
		q[i] = makePolynomial_AA(2, 1);
		mpq_init(q[i]->elems[0].coef);
		mpq_init(q[i]->elems[1].coef);
		mpq_set_si(q[i]->elems[0].coef, 1l, 1l);
		mpq_neg(q[i]->elems[1].coef, points[i]);
		q[i]->elems[0].degs = 1;
		q[i]->elems[1].degs = 0;
		q[i]->size = 2;
	}

	for (int i = 0; i < nPoints; ++i) {
		AltArr_t* prod = makeConstPolynomial_AA(1, 1, diffs[0]);
		for (int j = 0; j < nPoints; ++j) {
			if (i == j) {
				continue;
			}
			prod = multiplyByBinomial_AA_inp(prod, q[j]);
		}

		//before adding polynomials, update by mutiplying through by denom
		//note that denom is already inverted
		for (int j = 0; j < prod->size; ++j) {
			mpq_mul(prod->elems[j].coef, prod->elems[j].coef, denoms[i]);
		}

		//add should probably take care of this condensation for us.
		// condensePolyomial_AA(prod);
		// if (prod == NULL || prod->size == 0) {
		// 	freePolynomial_AA(prod);
		// 	continue;
		// }

		ret = addPolynomials_AA_inp(ret, prod, 1);
		freePolynomial_AA(prod);
	}

	for (int i = 0; i < nPoints; ++i) {
		freePolynomial_AA(q[i]);
	}
	free(q);
	for (int i = 0; i < nPoints; ++i) {
		mpq_clear(denoms[i]);
	}
	for (int i = 0; i < nPoints; ++i) {
		mpq_clear(diffs[i]);
	}

	return ret;
}

AltArr_t* univarInterpolateDoubles_AA(double* points, double* vals, int nPoints) {
	mpq_t mpqPoints[nPoints];
	mpq_t mpqVals[nPoints];

	for (int i = 0; i < nPoints; ++i) {
		mpq_init(mpqPoints[i]);
		mpq_init(mpqVals[i]);
		mpq_set_d(mpqPoints[i], points[i]);
		mpq_set_d(mpqVals[i], vals[i]);
	}

	AltArr_t* ret = univarInterpolate_AA(mpqPoints, mpqVals, nPoints);
	for (int i = 0; i < nPoints; ++i) {
		mpq_clear(mpqPoints[i]);
		mpq_clear(mpqVals[i]);
	}

	return ret;
}

void univarEvaluate_AA(AltArr_t* aa, const mpq_t point, mpq_t res) {

	if (aa == NULL || aa->size == 0) {
		mpq_set_si(res, 0l, 1l);
		return;
	}

	if (aa->nvar != 1) {
		fprintf(stderr, "Poly is not univariate in univarEvaluate_AA\n");
		exit(1);
	}

	if (aa->unpacked) {
		univarEvaluate_AA_unpk(aa, point, res);
		return;
	}

	AAElem_t* elems = aa->elems;
	register int size = aa->size;
	mpq_set(res, elems->coef);
	degrees_t prevDeg = elems->degs;
	degrees_t nextDeg;
	for (int i = 1; i < size; ++i) {
		nextDeg = elems[i].degs;
		for (degrees_t j = prevDeg; j > nextDeg; --j) {
			mpq_mul(res, res, point);
		}
		mpq_add(res, res, elems[i].coef);
		prevDeg = nextDeg;
	}
	for (degrees_t j = prevDeg; j > 0; --j) {
		mpq_mul(res, res, point);
	}
}




/****************
* Primitive Factorization w.r.t Rational Coefficients
*****************/

content_t* content (Node* a)
{
  Node* cur = a;
  mpzCoefs_t* n = (mpzCoefs_t*) malloc (sizeof (mpzCoefs_t));
  mpzCoefs_t* head;
  mpzCoefs_t* tail;
  mpzCoefs_t* cur_coef;
  content_t* gcd_res = (content_t*) malloc (sizeof (content_t));

  //make a mpzCoefs_t linked list which is coefs converted from mpq_t to mpz_t:
  mpz_init (n->mpzCoef);
  mpz_set (n->mpzCoef, mpq_numref (cur->coef));

  head = n;
  tail = n;

  cur = cur->next;
  while (cur != NULL)
    {
      n = (mpzCoefs_t*) malloc ( sizeof (mpzCoefs_t));
      mpz_init (n->mpzCoef);
      mpz_set (n->mpzCoef, mpq_numref(cur->coef));

      tail->next = n;
      tail = tail->next;

      cur = cur->next;
    }
  tail->next = NULL;

  // compute the content of mpz_t form of the coefs.
  cur_coef = head;
  mpz_init (gcd_res->cont);

  if (cur_coef->next != NULL)
    {
      mpz_gcd (gcd_res->cont, cur_coef->mpzCoef, cur_coef->next->mpzCoef);
      cur_coef = cur_coef->next;

      while (cur_coef->next != NULL)
	{
	  mpz_gcd (gcd_res->cont, gcd_res->cont, cur_coef->next->mpzCoef);
	  cur_coef = cur_coef->next;
	}
      return gcd_res;
    }
  else
    {
      mpz_set_ui (gcd_res->cont, 1ul);
      return gcd_res;
    }
}

Node* primitivePart (Node* a, content_t** returnContent,  int nvar)
{
  Node* cur = a;
  content_t* cont_tmp = content (a);
  Node* n = (Node*) malloc (sizeof (Node));
  Node* head;
  Node* tail;

  n->degs = cur->degs;

  mpq_init(n->coef);
  ratNum_t set_num;
  mpq_init (set_num);
  mpq_set_num (set_num, cont_tmp->cont);
  mpq_div(n->coef, cur->coef, set_num);

  head = n;
  tail = n;
  cur = cur->next;

  while (cur != NULL)
    {
      n = (Node*) malloc (sizeof (Node));
      n->degs = cur->degs;

      mpq_init(n->coef);
      ratNum_t set_num;
      mpq_init (set_num);
      mpq_set_num (set_num, cont_tmp->cont);
      mpq_div(n->coef, cur->coef, set_num);

      tail->next = n;
      tail = tail->next;
      cur = cur->next;
    }

  tail->next = NULL;

  *returnContent = cont_tmp;

  return head;
}



/****************
* Normal Form (Multi-Divisor Division (MDD))
*****************/

AltArr_t* onlyNormalForm_AA (AltArr_t* f, AltArr_t** G, int s, int nvar)
{
	if (s < 0){
		return NULL;
	}

	if (s == 0){
		return deepCopyPolynomial_AA (f);
	}

	if (s == 1){
		if (G[0] != NULL && G[0]->size != 0){
			AltArr_t* r1 = NULL;
			AltArr_t* q1 = NULL;
			dividePolynomials_AA (f, *G, &q1, &r1, nvar);
			freePolynomial_AA (q1);

			return r1;
		}

		return deepCopyPolynomial_AA(f);
	}

	int i; // i= 0 ... (s-1)
	int isDivided;
	AltArr_t* h = deepCopyPolynomial_AA(f);  // this algorithm runs until h != NULL
	AltArr_t* rem = NULL;
	AltArr_t* lt = NULL;
	AltArr_t* tmp_q;
	AltArr_t* tmp_r;

	while (h != NULL && h->size != 0) {

		i = 0;
		isDivided = 0;

		while (i < s && !isDivided)
			{
				if ( h->size != 0 && h!= NULL &&
					G[i] != NULL && G[i]->size != 0 && monomialDivideTest_AA (h, 0, G[i], 0))  {
					tmp_q = NULL;
					tmp_r = NULL;
					divideByLeadingTerms_AA (h, G[i], &tmp_q, &tmp_r, nvar);
					freePolynomial_AA (tmp_r);

					if (!isOne_AA (tmp_q)) {
						tmp_q = multiplyPolynomials_AA_inp (tmp_q, G[i], nvar);
						h = subPolynomials_AA_inp (h, tmp_q, nvar);
						freePolynomial_AA (tmp_q);
					} else {
						freePolynomial_AA (tmp_q);
						h = subPolynomials_AA_inp (h, G[i], nvar);
					}
					isDivided = 1;
				} else
					i = i + 1;
			}

		if (!isDivided) {
			if (h != NULL && h->size != 0){
				lt = leadingTerm_AA (h, nvar);
			} else {
				lt = NULL;
			}

			if (rem == NULL) {
				rem = lt;
			} else {
				if (lt != NULL) {
					rem = addPolynomials_AA_inp (rem, lt, nvar);
					freePolynomial_AA (lt);
				}
			}
			if (h != NULL && h->size != 0) {
				/* h = subPolynomials_AA_inp (h, lt, nvar); */
				subByLeadingTerm_AA (&h, nvar);
			}
		}
	}

	return rem;
}

void heapMDD_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t** r, int s, int nvar)
{
	if (s < 0){
		fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0){
		*r = deepCopyPolynomial_AA (f);
		return;
	}

	if (s == 1){  // simple multivariate polynomial division
		if (G[0] != NULL && G[0]->size != 0){
			dividePolynomials_AA (f, *G, Q, r, nvar);
		} else {
			*r = deepCopyPolynomial_AA(f);
			Q[0] = NULL;
		}
		return;
	}

	int i; // i= 0 ... (s-1)
	int isDivided;
	AltArr_t* h = deepCopyPolynomial_AA(f);  // this algorithm runs until h != NULL
	AltArr_t* rem = NULL;
	AltArr_t* lt = NULL;
	AltArr_t* tmp_q;
	AltArr_t* tmp_r;

	while (h != NULL && h->size != 0) {

		i = 0;
		isDivided = 0;

		while (i < s && !isDivided)
			{
				if ( h->size != 0 && h!= NULL && G[i] != NULL && G[i]->size != 0  && monomialDivideTest_AA(h, 0, G[i], 0))  {
					tmp_q = NULL;
					tmp_r = NULL;
					divideByLeadingTerms_AA (h, G[i], &tmp_q, &tmp_r, nvar);
					freePolynomial_AA (tmp_r);

					if (Q[i] != NULL) {
						Q[i] = addPolynomials_AA_inp (Q[i], tmp_q, nvar);
					} else {
						Q[i] = deepCopyPolynomial_AA (tmp_q);
					}

					if (!isOne_AA (tmp_q)) {
						tmp_q = multiplyPolynomials_AA_inp (tmp_q, G[i], nvar);
						h = subPolynomials_AA_inp (h, tmp_q, nvar);
						freePolynomial_AA (tmp_q);
					} else {
						freePolynomial_AA (tmp_q);
						h = subPolynomials_AA_inp (h, G[i], nvar);
					}
					/* h = tmp_r; */
					/* i = 0; */
					isDivided = 1;
				} else
					i = i + 1;
			}

		if (!isDivided) {
			if (h != NULL && h->size != 0){
				lt = leadingTerm_AA (h, nvar);
			} else {
				lt = NULL;
			}

			if (rem == NULL) {
				rem = lt;
			} else {
				if (lt != NULL) {
					rem = addPolynomials_AA_inp (rem, lt, nvar);
					freePolynomial_AA (lt);
				}
			}
			if (h != NULL && h->size != 0) {
				/* h = subPolynomials_AA_inp (h, lt, nvar); */
				subByLeadingTerm_AA (&h, nvar);
			}
		}
	}

	*r = rem;
	return;
}

// Note this is a test
// See the recursiveTeriangularSetMDD_AA in SMQP_Support_Recursive-AA.h
void recTriangularSetMDD_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t** r, int s, int nvar)
{
	if (s < 0){
	    fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
	    return;
	}

	if (s == 0){
	    *r = deepCopyPolynomial_AA(f);
	    return;
	}

	if (s == 1){  // simple multivariate polynomial division
	    if (G[0] != NULL && G[0]->size != 0){
		dividePolynomials_AA (f, *G, Q, r, nvar);
	    } else {
		*r = deepCopyPolynomial_AA (f);
		Q[0] = NULL;
	    }
	    return;
	}

	int index;
	int i;
	int* orderList = recursiveLoop (s);
	int orderSize = (1 << s) - 1;
	AltArr_t* h = deepCopyPolynomial_AA (f);
	AltArr_t* rem = NULL;
	AltArr_t* lt = NULL;
	AltArr_t* tmp_q = NULL;
	AltArr_t* tmp_r = NULL;

	while (h != NULL && h->size != 0){
		for (i = 0; i < orderSize; i++) {
			index = orderList[i];
			if (h->size != 0 && h != NULL && G[index] != NULL  && G[index]->size != 0
				&& monomialDivideTest_AA(h, 0, G[index], 0)){
			    tmp_q = NULL;
			    tmp_r = NULL;

			    dividePolynomials_AA (h, G[index], &tmp_q, &tmp_r, nvar);
			    if (Q[index] != NULL) {
					Q[index] = addPolynomials_AA_inp (Q[index], tmp_q, nvar);
			    }
			    else {
					Q[index] = tmp_q;
			    }
			    freePolynomial_AA(h);
			    h = tmp_r;
			}
		}

		if (h->size != 0 && h != NULL)
		    lt = leadingTerm_AA (h, nvar);
		else
		    lt = NULL;

		if (rem == NULL)
		    rem = lt;
		else {
		    if (lt != NULL){
				rem = addPolynomials_AA_inp (rem, lt, nvar);
		    }
		}
		if (h != NULL && h->size != 0){
		    /* h = subPolynomials_AA_inp (h, lt, nvar); */
			subByLeadingTerm_AA (&h, nvar);
		}
	}

	*r = rem;
	free(orderList);
	return;
}

// Look at the recursiveTriangularSetMDD_AA in SMQP_Support_Recursive-AA.h
void triangularSetMDD_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t** r, int s, int nvar)
{
	if (s < 0){
		fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0){
		*r = deepCopyPolynomial_AA (f);
		return;
	}

	if (s == 1){  // simple multivariate polynomial division
		if (G[0] != NULL && G[0]->size != 0){
			dividePolynomials_AA (f, *G, Q, r, nvar);
		} else {
			*r = deepCopyPolynomial_AA(f);
			Q[0] = NULL;
		}
		return;
	}

	int i; // i= 0 ... (s-1)
	int isDivided;
	AltArr_t* h = deepCopyPolynomial_AA(f);
	AltArr_t* rem = NULL;
	AltArr_t* lt = NULL;
	AltArr_t* tmp_q;
	AltArr_t* tmp_r;

	while (h != NULL && h->size != 0) {

		i = 0;
		isDivided = 0;

		while (i < s && !isDivided)
			{
				if ( h->size != 0 && h!= NULL && G[i] != NULL && G[i]->size != 0  && monomialDivideTest_AA(h, 0, G[i], 0))  {
					tmp_q = NULL;
					tmp_r = NULL;
					divideByLeadingTerms_AA (h, G[i], &tmp_q, &tmp_r, nvar);
					freePolynomial_AA (tmp_r);

					if (Q[i] != NULL) {
						Q[i] = addPolynomials_AA_inp (Q[i], tmp_q, nvar);
					} else {
						Q[i] = deepCopyPolynomial_AA (tmp_q);
					}

					if (!isOne_AA (tmp_q)) {
						tmp_q = multiplyPolynomials_AA_inp (tmp_q, G[i], nvar);
						h = subPolynomials_AA_inp (h, tmp_q, nvar);
						freePolynomial_AA (tmp_q);
					} else {
						freePolynomial_AA (tmp_q);
						h = subPolynomials_AA_inp (h, G[i], nvar);
					}
					isDivided = 1;
				} else
					i = i + 1;
			}

		if (!isDivided) {
			if (h != NULL && h->size != 0){
				lt = leadingTerm_AA (h, nvar);
			} else {
				lt = NULL;
			}

			if (rem == NULL) {
				rem = lt;
			} else {
				if (lt != NULL) {
					rem = addPolynomials_AA_inp (rem, lt, nvar);
					freePolynomial_AA (lt);
				}
			}
			if (h != NULL && h->size != 0) {
				subByLeadingTerm_AA (&h, nvar);
			}
		}
	}

	*r = rem;
	return;
}

void primitiveFactorTriangularSetMDD_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t** r, int s, int nvar)
{
	if (s < 0){
		fprintf (stderr, "SMQP Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0){
		*r = deepCopyPolynomial_AA (f);
		return;
	}

	if (s == 1){  // simple multivariate polynomial division
		if (G[0] != NULL && G[0]->size != 0){
			dividePolynomials_AA (f, *G, Q, r, nvar);
		} else {
			*r = deepCopyPolynomial_AA (f);
			Q[0] = NULL;
		}
		return;
	}

	mpq_t cont;
	mpq_init (cont);
	AltArr_t* fpp = primitivePartAndContent_AA (f, cont);
	/* AltArr_t** Gpp = (AltArr_t**) malloc (sizeof (AltArr_t*)*s); */
	/* for (int i = 0;  i < s; i++) { */
	/* 	Gpp[i] = primitivePart_AA (G[i]); */
	/* } */

    triangularSetMDD_AA (fpp, G, Q, r, s, nvar);

	mpq_t invC;
	mpq_init (invC);
	mpq_inv (invC, cont);
	multiplyByRational_AA_inp (*r, invC);
	for (int i = 0; i < s; i++) {
		if (Q[i] != NULL && Q[i]->size != 0) {
			multiplyByRational_AA_inp (Q[i], invC);
		}
	}

	mpq_clear (cont);
	mpq_clear (invC);
	freePolynomial_AA (fpp);
	/* for (int i = 0; i < s; i++) { */
	/* 	freePolynomial_AA (Gpp[i]); */
	/* } */
	/* free (Gpp); */
}


void normalForm_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t** r, int s, int nvar) {
	heapMDD_AA (f, G, Q, r, s, nvar);
}

void multiDivisorDivision_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t** r, int s, int nvar, int type)
{
	if (type == 0){
		heapMDD_AA (f, G, Q, r, s, nvar);
		return;
	} else if (type == 1){
		triangularSetMDD_AA (f, G, Q, r, s, nvar);
		return;
	} else{
		primitiveFactorTriangularSetMDD_AA (f, G, Q, r, s, nvar);
	}
	return;
}

int multiDivisorDivisionVerification_AA (AltArr_t* f, AltArr_t** G, AltArr_t** Q, AltArr_t* r, AltArr_t* hPow, int nSet, int tnvar)
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


    AltArr_t* sum = NULL;
    for (int i = 0; i < nSet; ++i){
        if (G[i] != NULL && G[i]->size != 0 && Q[i] != NULL && Q[i]->size != 0){
            if (sum == NULL || sum->size == 0){
                sum = multiplyPolynomials_AA (G[i], Q[i], tnvar);
            } else {
            	sum = addPolynomials_AA_inp (sum , multiplyPolynomials_AA (G[i], Q[i], tnvar), tnvar);
            }
        }
    }
    if (r != NULL && r->size != 0){
        sum = addPolynomials_AA_inp (sum, r, tnvar);
    }

    AltArr_t* hf = NULL;
    if (hPow != NULL && hPow->size != 0){
         hf = multiplyPolynomials_AA (f, hPow, tnvar);
    } else {
        hf = deepCopyPolynomial_AA (f);
    }

    AltArr_t* res = subPolynomials_AA (sum, hf, tnvar);
    if (res == NULL || res->size == 0){
        return 1;
    } else {
        fprintf(stderr, "No. Zero-Poly(s): %d\n", nullity);
        printf("The difference: \n");
        printAA (res);
        printf("\n");
        return 0;
    }
}
