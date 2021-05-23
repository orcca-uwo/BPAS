

#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "IntegerPolynomial/SMZP_Support_Test.h"
#include <time.h>

#if defined(WITH_BLAD) && WITH_BLAD
#include "BLADInterface/bladinterface.h"
#endif

#if defined(WITH_MAPLE) && WITH_MAPLE
#include "Parser/bpas_parser.h"
#endif

void printAAZ (AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
	    fprintf(stderr, "0\n");
		return;
	}
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		gmp_fprintf(stderr, "%Zd*%llx + ", aa->elems[i].coef, aa->elems[i].degs);
	}
	fprintf (stderr, "\n");
}



AltArrZ_t* buildRandomAltArrZPoly_MvarSparse(int univarSparsity, int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg, time_t seed) {
	AltArrZ_t* ret = buildRandomSeededPoly_AAZ_unpk(nvar, nterms, coefBound, sparsity, includeNeg, seed);

	degree_t* degs = (degree_t*) ret->elems->degs;

	degree_t pDeg = degs[0];
	for (int i = 1; i < nvar; ++i) {
		degree_t deg = partialDegree_AAZ(ret, i);
		pDeg = deg > pDeg ? deg : pDeg;
	}

	int coefSize = nterms/(pDeg+1);
	int curCoef = 0;
	degree_t curDeg = pDeg;
	for (int i = nterms-1; i >= 0; --i) {
		degs[i*nvar] = curDeg;
		++curCoef;
		if (curDeg > 0 && curCoef >= coefSize) {
			--curDeg;
			// curDeg -= 2;
			// curDeg = curDeg < 0 ? 0 : curDeg;
			curCoef = 0;
		}
		// degs[i*nvar] = (rand() % pDeg);
	}

	canonicalizePolynomial_AAZ(ret);
	tryPackExponentVectors_AAZ_inp(ret);


	return ret;
}

/**
 * Build a random recursive polynomial.
 * Nterms is the number of terms wrt the main variable.
 * coefNTerms is the number of terms in the coefficients (when viewed recursively, and nvar > 1);
 * coefBound is the absolute upper limit on the rational numbers in the polynomial.
 * sparsity is the difference between degrees (wrt to main variable) of successive terms - 1. Therefore 2 is dense.
 * sparsity is also passed when constructing coefficient polynomials.
 */
RecArrZ_t* buildRandomRecArrZPoly(int nvar, int nterms, int coefNterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	static int initRand = 0;
	static gmp_randstate_t R_STATE;

	if (!initRand) {
		time_t seed = time(NULL);
		// fprintf(stderr, "rec seed: %lu\n", seed);
		srand(seed);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, seed);

		initRand = 1;
	}

	int mvarDegOffset = getMVarExpOffset(nvar);

	if (nvar == 1) {
		Node* head = NULL, *tail = NULL;
		degrees_t mvarDeg;

		mpz_t mpzNum, mpzDen;
		mpz_init(mpzNum);
		mpz_init(mpzDen);

		for (int i = 0; i < nterms; ++i) {
			Node* coefNode = (Node*) malloc(sizeof(Node));
			coefNode->next = NULL;
			coefNode->degs = 0;
			mpq_init(coefNode->coef);

			mpz_urandomb(mpzNum, R_STATE, coefBound);
			mpz_urandomb(mpzDen, R_STATE, coefBound);

			while(mpz_sgn(mpzNum) == 0) {
				mpz_urandomb(mpzNum, R_STATE, coefBound);
			}
			while(mpz_sgn(mpzDen) == 0) {
				mpz_urandomb(mpzDen, R_STATE, coefBound);
			}
			if (includeNeg && rand() % 2) {
				//50/50 chance of being negative
				mpz_neg(mpzNum, mpzNum);
			}

			mpz_set(mpq_numref(coefNode->coef), mpzNum);
			mpz_set(mpq_denref(coefNode->coef), mpzDen);
			mpq_canonicalize(coefNode->coef);

			mvarDeg = (nterms - 1 - i) * (sparsity - 1);
			coefNode->degs += (mvarDeg << mvarDegOffset);

			if (head == NULL) {
				head = coefNode;
			}
			if (tail != NULL) {
				tail->next = coefNode;
			}
			tail = coefNode;
		}

		//Node* head is now the poly we want to make recursive. Let's turn it into an AA.
		AltArrZ_t* aa = deepCopyPolynomial_AAZFromNode(head, nvar);
		freePolynomial(head);
		return convertToRecursiveArrayZ(aa);
	} else {

		Node* head = NULL, *tail = NULL;
		degrees_t mvarDeg;
		int mvarDegOffset = getMVarExpOffset(nvar);


		for (int i = 0; i < nterms; ++i) {
			Node* coefNode = buildRandomZPoly(nvar-1, coefNterms, coefBound, sparsity, includeNeg);
			mvarDeg = (nterms - 1 - i) * (sparsity - 1);
			if (head == NULL) {
				head = coefNode;
			}
			if (tail != NULL) {
				tail->next = coefNode;
			}

			Node* curNode = coefNode;
			while (curNode != NULL) {
				//expand node by 1 var.
				curNode->degs += (mvarDeg << mvarDegOffset);
				if (curNode->next == NULL) {
					tail = curNode;
				}
				curNode = curNode->next;
			}
		}

		//Node* head is now the poly we want to make recursive. Let's turn it into an AA.
		AltArrZ_t* aa = deepCopyPolynomial_AAZFromNode(head, nvar);
		freePolynomial(head);
		return convertToRecursiveArrayZ(aa);
	}
}

RecArrZ_t* convertToRecursiveArrayZ(AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return convertToRecursiveArrayZ_unpk(aa);
	}

	int mvarDegOffset = getMVarExpOffset(aa->nvar);
	degrees_t mvarMask = getMVarExpMask(aa->nvar);

	AAZElem_t* elems = aa->elems;
	degree_t mvarDeg = GET_NTH_EXP(elems[0].degs, mvarMask, mvarDegOffset);
	RecArrZ_t* poly = (RecArrZ_t*) malloc(sizeof(RecArrZ_t));
	poly->alloc = mvarDeg;
	RecArrElemZ_t* recElems = (RecArrElemZ_t*) malloc(sizeof(RecArrElemZ_t)*(mvarDeg+1));
	poly->elems = recElems;

	int curIdx = 0, lastSize = 0;
	degree_t curDeg;
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		curDeg = GET_NTH_EXP(elems[i].degs, mvarMask, mvarDegOffset);
		if (curDeg != mvarDeg) {
			recElems[curIdx].exp = mvarDeg;
			recElems[curIdx].coefSize = i - lastSize;
			recElems[curIdx].coef = elems + lastSize;

			++curIdx;
			lastSize = i;
			mvarDeg = curDeg;
		}

		//throw away exp of mvar for coef.
		elems[i].degs = elems[i].degs & (~mvarMask);
	}

	//one more to commit after loop
	recElems[curIdx].exp = mvarDeg;
	recElems[curIdx].coefSize = AA_SIZE(aa) - lastSize;
	recElems[curIdx].coef = elems + lastSize;
	++curIdx;

	poly->unpacked = 0;
	poly->size = curIdx;
	poly->origAA = aa;

	return poly;
}

AltArrZ_t* convertFromRecursiveArrayZ(RecArrZ_t* poly, int nvar) {
	if (poly == NULL) {
		return NULL;
	}

	if (poly->unpacked) {
		return convertFromRecursiveArrayZ_unpk(poly, nvar);
	}

	AltArrZ_t* aa = NULL;
	if (poly->origAA != NULL) {
		aa = poly->origAA;
	} else {
		aa = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	}

	int mvarDegOffset = getMVarExpOffset(nvar);

	RecArrElemZ_t* recElems = poly->elems;
	AAZElem_t* elems = poly->elems->coef;
	aa->elems = elems;
	register int curSize = 0;
	for (int i = 0; i < poly->size; ++i) {
		degrees_t curDeg = recElems[i].exp;
		for (int j = 0; j < recElems[i].coefSize; ++j) {
			elems[curSize + j].degs += (curDeg << mvarDegOffset);
		}

		curSize += recElems[i].coefSize;
	}

	freeRecArrayZ(poly);

	aa->size = curSize;
	aa->alloc = curSize;
	aa->nvar = nvar;
	aa->unpacked = 0;

	return aa;
}

RecArrZ_t* deepCopyRecArrayZPolynomial(RecArrZ_t* poly, int nvar) {
	if (poly == NULL || poly->size == 0) {
		return NULL;
	}
	int maxCoefSize = 0;
	RecArrElemZ_t* recElems = poly->elems;
	for (int i = 0; i < poly->size; ++i) {
		maxCoefSize += recElems[i].coefSize;
	}

	AltArrZ_t localAA;
	localAA.unpacked = poly->unpacked;
	localAA.elems = poly->elems->coef;
	localAA.size = localAA.alloc = maxCoefSize;
	localAA.nvar = nvar;
	AltArrZ_t* copyAA = deepCopyPolynomial_AAZ(&localAA);

	AAZElem_t* elems = copyAA->elems;
	free(copyAA); //free the container but not underlying elems.

	// AAZElem_t* origElems = poly->elems->coef;
	// AAZElem_t* elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*maxCoefSize);
	// for (int i = 0; i < maxCoefSize; ++i) {
	// 	mpz_init(elems[i].coef);
	// 	mpz_set(elems[i].coef, origElems[i].coef);
	// 	elems[i].degs = origElems[i].degs;
	// }

	RecArrElemZ_t* copyElems = (RecArrElemZ_t*) malloc(sizeof(RecArrElemZ_t)*poly->size);
	memcpy(copyElems, poly->elems, sizeof(RecArrElemZ_t)*poly->size);
	copyElems->coef = elems;
	maxCoefSize = poly->elems->coefSize;
	for (int i = 1; i < poly->size; ++i) {
		copyElems[i].coef = elems + maxCoefSize;
		maxCoefSize += copyElems[i].coefSize;
	}

	RecArrZ_t* retPoly = (RecArrZ_t*) malloc(sizeof(RecArrZ_t));
	retPoly->size = retPoly->alloc = poly->size;
	retPoly->elems = copyElems;
	retPoly->unpacked = poly->unpacked;

	retPoly->origAA = NULL;

	return retPoly;
}

RecArrZ_t* convertToRecursiveArrayZAtIdx (AltArrZ_t* aa, int idx)
{
  if (aa == NULL){
      return NULL;
  }
  if (idx < 0 || idx >= aa->nvar){
    fprintf (stderr, "SMZP Error: idx is out of range!\n");
    exit (1);
  }

  AltArrZ_t* cPoly = deepCopyPolynomial_AAZ (aa);
  if (idx == 0){
    return convertToRecursiveArrayZ (cPoly);
  }

  int varMap[cPoly->nvar];
  for (int i = 0; i < cPoly->nvar; ++i){
    varMap[i] = i;
  }
  varMap[0] = idx;
  varMap[idx] = 0;

  reorderVars_AAZ (cPoly, varMap, cPoly->nvar);
    return convertToRecursiveArrayZ (cPoly);
}

AltArrZ_t* convertFromRecursiveArrayZAtIdx (RecArrZ_t* recPoly, int idx, int nvar)
{
  if (recPoly == NULL){
      return NULL;
  }
  if (idx < 0 || idx >= nvar){
      fprintf (stderr, "SMZP Error: idx is out of range!\n");
      exit (1);
  }

  AltArrZ_t* cPoly = deepCopyPolynomial_AAZ (convertFromRecursiveArrayZ (recPoly, nvar));
  if (idx == 0){
      return cPoly;
  }

  int varMap[nvar];
  for (int i = 0; i < nvar; ++i){
    varMap[i] = i;
  }
  varMap[0] = idx;
  varMap[idx] = 0;

  reorderVars_AAZ (cPoly, varMap, nvar);
  return cPoly;
}


/***********
 * Pesudo Division
 ***********/
AltArrZ_t* recProdHeapGetNextCoef_AAZ(RecProdHeap_AA* h, const RecArrElemZ_t* __restrict__ aElems, const RecArrElemZ_t* __restrict__ bElems, const int* aUnpacked, int bUnpacked) {
	RecProdHeapElem_AA* elems = h->elements;
	int nvar = h->nvar;
	register int lastB = h->lastB;

	AltArrZ_t* ret = NULL;
	register degree_t maxExp = elems->exp;
	register degree_t nextExp = elems->exp;

	RecProdHeapChain_AA* insertChain = NULL;
	RecProdHeapChain_AA* maxElem, * nextMaxElem;

	AltArrZ_t aCoef, bCoef;
	aCoef.nvar = nvar;
	bCoef.nvar = nvar;
	//since b doesn't update through computation, it is statically either packed or unpacked
	bCoef.unpacked = bUnpacked;

	//It's possible that not all elemtns of same degree are chained. So much loop like this.
	while(nextExp == maxExp) {
		maxElem = recProdHeapRemoveMax_AA(h);

		//go through the chain now;
		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;

			aCoef.unpacked = aUnpacked[maxElem->a_i];
			aCoef.alloc = aCoef.size = aElems[maxElem->a_i].coefSize;
			aCoef.elems = aElems[maxElem->a_i].coef;
			bCoef.alloc = bCoef.size = bElems[maxElem->b].coefSize;
			bCoef.elems = bElems[maxElem->b].coef;

			//Notice here we delay the calculation of the actual product coef.
			//This is because a_i will change over the course of the algorithm.
			if (ret == NULL) {
				ret = multiplyPolynomials_AAZ(&aCoef, &bCoef, nvar);
			} else {
				AltArrZ_t* prod = multiplyPolynomials_AAZ(&aCoef, &bCoef, nvar);
				ret = addPolynomials_AAZ_inp(ret, prod, nvar);
				freePolynomial_AAZ(prod);
			}
			if (maxElem->b != lastB) {
				++(maxElem->b);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				maxElem->next = NULL;
				recProdHeapFreeChain_AA(maxElem);
			}
			maxElem = nextMaxElem;
		}

		nextExp = recProdHeapPeek_AA(h);
	}

	while(insertChain != NULL) {
		nextMaxElem = insertChain->next;
		insertChain->next = NULL;
		recProdHeapInsert_AA(h, insertChain, aElems[insertChain->a_i].exp + bElems[insertChain->b].exp);
		insertChain = nextMaxElem;
	}

	return ret;
}

/**
 * Multiply a recursively-viewed polynomial by a "coefficient" polynomial, in place.
 * That is, the coefficient supplied is a polynomial whose variables do not include
 * the mvar of the supplied RecNode_t a.
 *
 * For optimization purposes, it may be that there is a slot in the coefficient poly's
 * exponent vector for the mvar of a, but it should all be 0.
 *
 * Returns a pointer to the new head of a after multiplication.
 */

void recArrayMultiplyByCoefZ_inp(RecArrElemZ_t* aElems, int aSize, AltArrZ_t* coef, int nvar, int* aUnpacked) {
	AltArrZ_t* tempa = malloc(sizeof(AltArrZ_t));
	tempa->nvar = nvar;
	for (int i = 0; i < aSize; ++i) {
		tempa->unpacked = aUnpacked[i];
		tempa->elems = aElems[i].coef;
		tempa->alloc = tempa->size = aElems[i].coefSize;
		tempa = multiplyPolynomials_AAZ_inp(tempa, coef, nvar);
		aUnpacked[i] = tempa->unpacked;
		aElems[i].coef = tempa->elems;
		aElems[i].coefSize = tempa->size;
	}

	free(tempa);
}


void pesudoDivideOneTerm_RecArrayZ(RecArrZ_t* c, RecArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, AltArrZ_t** hPow, int nvar, int lazy) {

	if (b == NULL) {
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (b->elems->exp > c->elems->exp) {
		RecArrZ_t* recR = deepCopyRecArrayZPolynomial(c, nvar);
		*res_r = convertFromRecursiveArrayZ(recR, nvar);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (hPow != NULL) {
			*hPow = makePolynomial_AAZ(1, nvar);
			mpz_init((*hPow)->elems->coef);
			mpz_set_si((*hPow)->elems->coef, 1l);
			(*hPow)->elems->degs = 0;
			(*hPow)->size = 1;
		}
		return;
	}

	RecArrElemZ_t* __restrict__ k = c->elems;
	RecArrElemZ_t* __restrict__ lenK = k + c->size;

	register int mvarOffset = getMVarExpOffset(nvar);
	register int maxASize = k->coefSize;
	register int maxRSize = k->coefSize;
	register int ai = 0;
	register int rj = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);

	register degree_t bDeg = b->elems->exp;
	register degree_t eps = k->exp - bDeg;
	AltArrZ_t kCoef;
	kCoef.alloc = kCoef.size = k->coefSize;
	kCoef.elems = k->coef;
	kCoef.nvar = nvar;
	kCoef.unpacked = c->unpacked;

	int i = 1;
	AltArrZ_t h;
	h.alloc = h.size = b->elems->coefSize;
	h.elems = b->elems->coef;
	h.nvar = nvar;
	h.unpacked = b->unpacked;

	AltArrZ_t* multerm = deepCopyPolynomial_AAZ(&kCoef);
	AltArrZ_t* hI = deepCopyPolynomial_AAZ(&h);

	//do first division manually
	//update degs in-place by degree of the mainvar, eps.
	memcpy(a->elems + ai, multerm->elems, sizeof(AAZElem_t)*multerm->size);
	degree_t* aDegs;
	if (c->unpacked) {
		aDegs = malloc(sizeof(degree_t)*nvar*maxASize);
		memcpy(aDegs, (degree_t*) multerm->elems->degs, sizeof(degree_t)*nvar*multerm->size);
		for (int idx = ai; idx < ai + multerm->size; ++idx) {
			aDegs[idx*nvar] = eps;
			a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
		}
		a->unpacked = 1;
	} else {
		for (int idx = ai; idx < ai + multerm->size; ++idx) {
			a->elems[idx].degs += ((degrees_t)eps << mvarOffset);
		}
	}

	a->size = multerm->size;
	ai += multerm->size;
	++k;
	free(multerm->elems);
	free(multerm);

	while (k != lenK) {
		kCoef.alloc = kCoef.size = k->coefSize;
		kCoef.elems = k->coef;

		eps = k->exp;
// #if PDIVIDE_DIVISBLE_CHECK
// 		// in this case only, when the divisor is one term, divisibility check
// 		// will always pass
// 		multerm = deepCopyPolynomial_AAZ(&kCoef);
// #else
// 		multerm = multiplyPolynomials_AAZ(&kCoef, hI, nvar);
// #endif
		multerm = deepCopyPolynomial_AAZ(&kCoef);
		++k;

		if (eps >= bDeg) {
// #if PDIVIDE_DIVISBLE_CHECK
// 			// AltArr_t* multermQuo = NULL;
// 			// if (divideTest_AA(multerm, &h, &multermQuo, nvar)) {
// 				// freePolynomial_AA(multerm);
// 				// multerm = multermQuo;
// 			// } else {
// 			// }
// #else
// 				a = multiplyPolynomials_AAZ_inp(a, &h, nvar);

// 				//update h counter
// 				//this is only done when we actually add a new element to the quotient.
// 				++i;
// 				hI = multiplyPolynomials_AAZ_inp(hI, &h, nvar); //in-place wrt first arg.
// #endif

			eps -= bDeg;

			//resize if necessary for insertion to a
			if (ai + multerm->size > maxASize) {
				maxASize += (multerm->size << 1);
				resizePolynomial_AAZ(a, maxASize);
			}

			memcpy(a->elems + ai, multerm->elems, sizeof(AAZElem_t)*multerm->size);
			if(multerm->unpacked || a->unpacked) {
				if (!a->unpacked) {
					unpackExponentVectors_AAZ_inp(a);
				}
				if (!multerm->unpacked) {
					unpackExponentVectors_AAZ_inp(multerm);
				}
				aDegs = (degree_t*) a->elems->degs;
				memcpy(aDegs + (ai*nvar), (degree_t*) multerm->elems->degs, sizeof(degree_t)*nvar*multerm->size);
				for (int idx = ai; idx < ai + multerm->size; ++idx) {
					aDegs[idx*nvar] = eps;
					a->elems[idx].degs = (degrees_t) (aDegs + idx*nvar);
				}

				free((degree_t*) multerm->elems->degs);
			} else {
				//update degs in-place by degree of the mainvar, eps.
				for (int idx = ai; idx < ai + multerm->size; ++idx) {
					a->elems[idx].degs += ((degrees_t)eps << mvarOffset);
				}
			}

			ai += multerm->size;
			a->size = ai;;

		} else {
			//accumulate remainder.
			if (rj + multerm->size > maxRSize) {
				r->size = rj;
				maxRSize += (multerm->size << 1);
				resizePolynomial_AAZ(r, maxRSize);
			}

			//copy multerm coef data directly into r
			memcpy(r->elems + rj, multerm->elems, sizeof(AAZElem_t)*multerm->size);
			if (multerm->unpacked || r->unpacked) {
				if (!r->unpacked) {
					unpackExponentVectors_AAZ_inp(r);
				}
				if (!multerm->unpacked) {
					unpackExponentVectors_AAZ_inp(multerm);
				}

				degree_t* rDegs = (degree_t*) r->elems->degs;
				memcpy(rDegs + rj*nvar, (degree_t*) multerm->elems->degs, sizeof(degree_t)*nvar*multerm->size );
				for (int 	idx = rj; idx < rj + multerm->size; ++idx) {
					rDegs[idx*nvar] = eps;
					r->elems[idx].degs = (degrees_t) (rDegs + idx*nvar);
				}

				free((degree_t*) multerm->elems->degs);
			} else {
				//update degs in-place by degree of the mainvar, eps.
				for (int idx = rj; idx < rj + multerm->size; ++idx) {
					r->elems[idx].degs += ((degrees_t)eps << mvarOffset);
				}
			}

			rj += multerm->size;
			r->size = rj;

		}
		free(multerm->elems);
		free(multerm);
	}


	if (ai > 0) {
		a->size = ai;
		a->alloc = maxASize;
	} else {
		freePolynomial_AAZ(a);
		a = NULL;
	}

	if (rj > 0) {
		r->size = rj;
		r->alloc = maxRSize;
	} else {
		freePolynomial_AAZ(r);
		r = NULL;
	}

	if (!lazy) {
		int d = c->elems->exp - b->elems->exp + 1 - i;
		i += d;
		for (int j = 0; j < d; ++j) {
			if (ai > 0) {
				a = multiplyPolynomials_AAZ_inp(a, &h, nvar);
			}
			if (rj > 0) {
				r = multiplyPolynomials_AAZ_inp(r, &h, nvar);
			}
			hI = multiplyPolynomials_AAZ_inp(hI, &h, nvar);
		}
// #if PDIVIDE_DIVISBLE_CHECK
		//r needs one additional one since a and h have one extra from beginning of loop
		if (rj > 0) {
			r = multiplyPolynomials_AAZ_inp(r, &h, nvar);
		}
// #endif
	}

	*res_a = a;
	*res_r = r;

	//Return number of division steps;
	if (e != NULL) {
		*e = i;
	}

	//Return the initial to power of e in hPow;
	if (hPow != NULL) {
		*hPow = hI;
	}

}

void pesudoDivide_RecArrayZ(RecArrZ_t* c, RecArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, AltArrZ_t** hPow, int nvar, int lazy) {

	if (b == NULL || b->size == 0) {
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || c->size == 0) {
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (b->elems->exp > c->elems->exp) {
		RecArrZ_t* recR = deepCopyRecArrayZPolynomial(c, nvar);
		*res_r = convertFromRecursiveArrayZ(recR, nvar);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (hPow != NULL) {
			*hPow = makePolynomial_AAZ(1, nvar);
			mpz_init((*hPow)->elems->coef);
			mpz_set_si((*hPow)->elems->coef, 1l);
			(*hPow)->elems->degs = 0;
			(*hPow)->size = 1;
		}
		return;
	}

	if (b->size == 1) {
		pesudoDivideOneTerm_RecArrayZ(c, b, res_a, res_r, e, hPow, nvar, lazy);
		return;
	}

	RecArrElemZ_t* __restrict__ k = c->elems;
	RecArrElemZ_t* __restrict__ b2Elem = b->elems + 1;
	RecArrElemZ_t* __restrict__ lenK = k + c->size;

	register int mvarOffset = getMVarExpOffset(nvar);
	register int maxASize = c->size < 10 ? 10 : c->size;
	register int maxRSize = maxASize;
	register int ai = 0;
	register int rj = 0;

	int* aUnpacked = calloc(maxASize, sizeof(int));

	AAZElem_t* r = (AAZElem_t*) malloc(sizeof(AAZElem_t)*maxRSize);
	degree_t* rDegs = (degree_t*) malloc(sizeof(degree_t)*maxRSize*nvar);
	RecArrElemZ_t* aElems = (RecArrElemZ_t*) malloc(sizeof(RecArrElemZ_t)*maxASize);
	RecArrElemZ_t* __restrict__ curA = aElems;

	//manually do first div as we know it comes from first term of c;
	int i = 1;
	AltArrZ_t h;
	h.alloc = h.size = b->elems->coefSize;
	h.elems = b->elems->coef;
	h.nvar = nvar;
	h.unpacked = b->unpacked;

	degree_t bDeg = b->elems->exp, eps = k->exp - bDeg;
	AltArrZ_t kCoef;
	kCoef.alloc = kCoef.size = k->coefSize;
	kCoef.elems = k->coef;
	kCoef.nvar = nvar;
	kCoef.unpacked = c->unpacked;

	AltArrZ_t* multerm = deepCopyPolynomial_AAZ(&kCoef);

	curA->coef = multerm->elems;
	curA->coefSize = multerm->size;
	curA->exp = eps;
	aUnpacked[ai] = multerm->unpacked;
	++curA;
	++ai;
	free(multerm);
	++k;

	AltArrZ_t* multerm2 = NULL;
	degree_t delta = 1;
	degree_t kDeg = (k == lenK) ? -1 : k->exp;

	AltArrZ_t* hI = deepCopyPolynomial_AAZ(&h);

	RecProdHeap_AA* prodHeap = recProdHeapCreate_AA(nvar);
	recProdheapResize_AA(prodHeap, maxASize);
	prodHeap->lastB = b->size - 1;
	recProdHeapInsert_AA(prodHeap, recProdHeapMakeChain_AA(0, 1, NULL), aElems->exp + b2Elem->exp);
	//get the leading term of dividend quotient-product difference
	delta = recProdHeapPeek_AA(prodHeap);
	while (delta > -1 || kDeg > -1) {

		if (delta == -1 && kDeg == -1) {
			break;
		}

		if (delta > kDeg) {
			eps = delta;
			multerm = recProdHeapGetNextCoef_AAZ(prodHeap, aElems, b->elems, aUnpacked, b->unpacked);
			if (multerm == NULL || multerm->size == 0) {
				freePolynomial_AAZ(multerm);
				//the elements in the product with degree delta
				// ended up canceling out coefs

				delta = recProdHeapPeek_AA(prodHeap);
				continue;
			}
			negatePolynomial_AAZ(multerm);
		} else if (delta < kDeg) {
			eps = kDeg;
			kCoef.alloc = kCoef.size = k->coefSize;
			kCoef.elems = k->coef;

			//then ctilde is from dividend
			multerm = multiplyPolynomials_AAZ(&kCoef, hI, nvar);
			++k;
			kDeg = (k == lenK) ? -1 : k->exp;
		} else {
			eps = delta;
			//combine both
			multerm2 = recProdHeapGetNextCoef_AAZ(prodHeap, aElems, b->elems, aUnpacked, b->unpacked);
			if (multerm2 == NULL || multerm2->size == 0) {
				freePolynomial_AAZ(multerm2);

				delta = recProdHeapPeek_AA(prodHeap);
				continue;
			}
			kCoef.alloc = kCoef.size = k->coefSize;
			kCoef.elems = k->coef;

			multerm = multiplyPolynomials_AAZ(&kCoef, hI, nvar);
			multerm = subPolynomials_AAZ_inp(multerm, multerm2, nvar);
			freePolynomial_AAZ(multerm2);

			++k;

			kDeg = (k == lenK) ? -1 : k->exp;

			if (multerm == NULL || multerm->size == 0) {
				//if sub resulted in a zero then we we must get a new multerm
				freePolynomial_AAZ(multerm);

				delta = recProdHeapPeek_AA(prodHeap);
				continue;
			}
		}

		//multerm is now the leading coef (eps the degree) of the current difference
		if (eps >= bDeg) {
	#if PDIVIDE_DIVISBLE_CHECK
			AltArrZ_t* multermQuo = NULL;
			if (divideTest_AAZ(multerm, &h, &multermQuo, nvar)) {
				freePolynomial_AAZ(multerm);
				multerm = multermQuo;
			} else {
#endif
				recArrayMultiplyByCoefZ_inp(aElems, ai, &h, nvar, aUnpacked);

				//update h counter
				//this is only done when we actually add a new element to the quotient.
				++i;
				hI = multiplyPolynomials_AAZ_inp(hI, &h, nvar); //in-place wrt first arg.

				// IF the above rec node mult was not done in place, then
				// we would need to do the below. Since we delay the actual
				// calculation of the coef product from heap elements until just
				// before extraction.
				// recProdHeapUpdateQuoByH(prodHeap, h);
#if PDIVIDE_DIVISBLE_CHECK
			}
#endif
			//resize if necessary for insertion to a
			if (ai >= maxASize) {
				maxASize <<= 1;
				aElems = (RecArrElemZ_t*) realloc(aElems, sizeof(RecArrElemZ_t)*maxASize);
				aUnpacked = (int*) realloc(aUnpacked, sizeof(int)*maxASize);
				curA = aElems + ai;

				recProdheapResize_AA(prodHeap, maxASize);
			}
			curA->coef = multerm->elems;
			curA->coefSize = multerm->size;
			curA->exp = eps - bDeg;
			aUnpacked[ai] = multerm->unpacked;
			free(multerm);

			//insert b_2! Since we constructed a to exactly cancel b_1, useless.
			recProdHeapInsert_AA(prodHeap, recProdHeapMakeChain_AA(ai, 1, NULL), curA->exp + b2Elem->exp);
			++ai;
			++curA;

		} else {
			//accumulate remainder.
			if (rj + multerm->size > maxRSize) {
				maxRSize += (multerm->size << 1);
				r = (AAZElem_t*) realloc(r, sizeof(AAZElem_t)*maxRSize);
				degree_t* oldDegs = rDegs;
				rDegs = (degree_t*) realloc(rDegs, sizeof(degree_t)*maxRSize*nvar);
				if (oldDegs != rDegs) {
					for (int idx = 0; idx < rj; ++idx) {
						r[idx].degs = (degrees_t) (rDegs + idx*nvar);
					}
				}
			}
			//copy multerm coef data directly into r
			memcpy(r + rj, multerm->elems, sizeof(AAZElem_t)*multerm->size);
			unpackExponentVectors_AAZ_inp(multerm);
			memcpy(rDegs + rj*nvar, (degree_t*) multerm->elems->degs, sizeof(degree_t)*multerm->size*nvar);

			//update degs in-place by degree of the mainvar, eps.
			for (int idx = rj; idx < rj + multerm->size; ++idx) {
				r[idx].degs = (degrees_t) (rDegs + idx*nvar);
				rDegs[idx*nvar] = eps;
				// r[idx].degs += ((degrees_t)eps << mvarOffset);
			}

			rj += multerm->size;

			free((degree_t*) multerm->elems->degs); //multerm got unpacked
			free(multerm->elems);
			free(multerm);
		}

		delta = recProdHeapPeek_AA(prodHeap);
	}

	//Condense coefs so they lay in the same array.
	int aSize = 0;
	int allPacked = 1;
	for (int idx = 0; idx < ai; ++idx) {
		aSize += aElems[idx].coefSize;
		if (allPacked && aUnpacked[idx]) {
			allPacked = 0;
		}
	}
	AAZElem_t* aBlock = (AAZElem_t*) malloc(sizeof(AAZElem_t)*aSize);

	if (allPacked) {
		aSize = 0;
		//first loop over each recursive element
		for (int idx = 0; idx < ai; ++idx) {
			memcpy(aBlock+aSize, aElems[idx].coef, sizeof(AAZElem_t)*(aElems[idx].coefSize));

			//loop over each element of coefficient and in-place update degs by mainvar exp.
			for (int idx2 = 0; idx2 < aElems[idx].coefSize; ++idx2) {
				aBlock[aSize + idx2].degs += ((degrees_t) aElems[idx].exp << mvarOffset);
			}

			free(aElems[idx].coef);
			aSize += aElems[idx].coefSize;
		}
		free(aElems);
	} else {
		degree_t* aDegs = (degree_t*) malloc(sizeof(degree_t)*nvar*aSize);
		aSize = 0;

		const degrees_t* masks = getExpMaskArray(nvar);
		const int* sizes = getExpOffsetArray(nvar);

		for (int idx = 0; idx < ai; ++idx) {
			//copy all coef data directly
			memcpy(aBlock+aSize, aElems[idx].coef, sizeof(AAZElem_t)*aElems[idx].coefSize);
			if (aUnpacked[idx]) {
				memcpy(aDegs + aSize*nvar, (degree_t*) aElems[idx].coef->degs, sizeof(degree_t)*aElems[idx].coefSize*nvar);
			} else {
				//we need to unpack this coef. Let's do it manually here.
				//start at var idx 1 since the first is 0 now
				for (int idx2 = 0; idx2 < aElems[idx].coefSize; ++idx2) {
					for (int idx3 = 1; idx3 < nvar; ++idx3) {
						aDegs[(aSize+idx2)*nvar + idx3] = GET_NTH_EXP(aElems[idx].coef[idx2].degs, masks[idx3], sizes[idx3]);
					}
				}
			}

			for (int idx2 = 0; idx2 < aElems[idx].coefSize; ++idx2) {
				aDegs[(aSize+idx2)*nvar] = aElems[idx].exp;
				aBlock[aSize + idx2].degs = (degrees_t) (aDegs + (aSize+idx2)*nvar);
			}
			aSize += aElems[idx].coefSize;
		}
	}

	AltArrZ_t* aa = NULL;
	if (aSize > 0) {
		aa = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
		aa->elems = aBlock;
		aa->size = aSize;
		aa->alloc = aSize;
		aa->nvar = nvar;
		aa->unpacked = !allPacked;
	}
	AltArrZ_t* ra = NULL;
	if (rj > 0) {
		ra = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
		ra->elems = r;
		ra->size = rj;
		ra->alloc = maxRSize;
		ra->nvar = nvar;
		ra->unpacked = 1;
		tryPackExponentVectors_AAZ_inp(ra);
	}

	if (!lazy) {
		int d = c->elems->exp - b->elems->exp + 1 - i;
		i += d;
		for (int j = 0; j < d; ++j) {
			aa = aa == NULL ? NULL : multiplyPolynomials_AAZ_inp(aa, &h, nvar);
			ra = ra == NULL ? NULL : multiplyPolynomials_AAZ_inp(ra, &h, nvar);
			hI = multiplyPolynomials_AAZ_inp(hI, &h, nvar);
		}
	}

	*res_a = aa;
	*res_r = ra;

	//Return number of division steps;
	if (e != NULL) {
		*e = i;
	}

	//Return the initial to power of e in hPow;
	if (hPow != NULL) {
		*hPow = hI;
	}
}


//////////////////////////////////
// Multi-divisor Pseudo-division
//////////////////////////////////

void pesudoDivideAtIdx_AAZ (int idx, AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, AltArrZ_t** hPow, int nvar, int lazy)
{
    if (idx >= nvar || idx < 0){
    	printAAZ (c);
		printAAZ (b);
		fprintf(stderr, "SMQP Error: idx(=%d) is out of range while nvar is %d!", idx, nvar);
		exit (EXIT_FAILURE);
    }

    // convert c and b to RecArrZ_t at idx:
    RecArrZ_t* tmpB = convertToRecursiveArrayZAtIdx (b, idx);
    RecArrZ_t* tmpC = convertToRecursiveArrayZAtIdx (c, idx);

    // compute pesudoDivide_RecArr(.):
    AltArrZ_t* q = NULL;
    AltArrZ_t* r = NULL;
    AltArrZ_t* h = NULL;
    int ee = 0;

    pesudoDivide_RecArrayZ (tmpC, tmpB, &q, &r, &ee, &h, nvar, lazy);
    if (idx > 0){
    	// freeRecArrayZ(tmpC);
    	// freeRecArrayZ(tmpB);
		freeRecArrayZAndCoef (tmpC);
		freeRecArrayZAndCoef (tmpB);
    }

    // convert q, r, and hPow at idx and set the outputs:
    if (idx > 0){
		*res_a = swappingExponents_AAZ (q, 0, idx);
		*res_r = swappingExponents_AAZ (r, 0, idx);
		*hPow = swappingExponents_AAZ (h, 0, idx);
		*e = ee;
		freePolynomial_AAZ (q);
		freePolynomial_AAZ (r);
		freePolynomial_AAZ (h);
    } else {
		*res_a = q;
		*res_r = r;
		*hPow = h;
		*e = ee;
    }
}

void naiveMultiDivisorPseudoDivide_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet)
{
    if (nSet  < 0 || nvar < nSet){
		fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
		exit(EXIT_FAILURE);
    }

    mpz_t one;
    mpz_init (one);
    mpz_set_si (one, 1l);

    if (nSet == 0){
		*rem = deepCopyPolynomial_AAZ (c);
		if (*hPow != NULL){
			freePolynomial_AAZ (*hPow);
		}
		*hPow = makeConstPolynomial_AAZ (1, nvar, one);
		quoSet[0] = NULL;
		return;
    }

    int e;
    int mvar;
    AltArrZ_t* h;
    AltArrZ_t* tmpQ;
    AltArrZ_t* tmpR;
    AltArrZ_t* r = deepCopyPolynomial_AAZ (c);
    AltArrZ_t* totalH = NULL;

    /* for (int i = 0; i < nSet; ++i)	 */
    for (int i = nSet-1; i >= 0; --i) {
	    mvar = leadingVariable_AAZ (B[i]);
	    if ((mvar > -1) && (r != NULL && r->size != 0)) {
		    e = 1;
		    h = NULL;
		    tmpQ = NULL;
		    tmpR = NULL;

		    pesudoDivideAtIdx_AAZ (mvar, r, B[i], &tmpQ, &tmpR, &e, &h, nvar, lazy);

		    if (tmpR != NULL && tmpR->size != 0) {
				r = tmpR; // update dividend
			} else {
				r = NULL;
			}

		    if (totalH == NULL || totalH->size == 0) {
				totalH = deepCopyPolynomial_AAZ (h);
		    } else if (h != NULL && h->size != 0){
				totalH = multiplyPolynomials_AAZ_inp (totalH, h, nvar);
		    }

		    if (h != NULL && h->size != 0 && !isOne_AAZ (h)){
				for (int j = 0; j < nSet; ++j){
					if (quoSet[j] != NULL && quoSet[j]->size != 0){
						quoSet[j] = multiplyPolynomials_AAZ_inp (quoSet[j], h, nvar);
					}
				}
		    }

		    if (tmpQ != NULL && tmpQ->size != 0){
		        if (quoSet[i] == NULL || quoSet[i]->size == 0){
		            quoSet[i] = tmpQ;
		        } else {
		            quoSet[i] = addPolynomials_AAZ_inp (quoSet[i], tmpQ, nvar);
					freePolynomial_AAZ (tmpQ);
		        }
		    }

		    freePolynomial_AAZ (h);
		}
	}

    if (totalH == NULL || totalH->size == 0) {
		totalH = makeConstPolynomial_AAZ (1, nvar, one);
    }

	mpz_clear (one);

    *hPow = totalH;
    *rem = r;

    return;
}

void normalizedTriangularSetPseudoDivide_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t*** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet)
{
    if (nSet < 0){
        fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
        exit (EXIT_FAILURE);
	}

	AltArrZ_t** Q; // *quoSet = Q

    if (c == NULL || c->size == 0) {
		mpz_t one1;
		mpz_init (one1);
		mpz_set_si (one1, 1l);

		*rem = NULL;
		*hPow = makeConstPolynomial_AAZ (1, nvar, one1);

		if (nSet > 0) {
			Q = (AltArrZ_t**) calloc (nSet, sizeof (AltArrZ_t*));

			for (int i = 0; i < nSet; ++i){
			    Q[i] = NULL;
			}
			*quoSet = Q;

		} else {
			*quoSet = NULL;
		}

        return;
    }

    if (nSet == 0) {
    	mpz_t one2;
		mpz_init (one2);
		mpz_set_si (one2, 1l);

		*rem = deepCopyPolynomial_AAZ (c);
		*hPow = makeConstPolynomial_AAZ (1, nvar, one2);
		*quoSet = NULL;

		return;
    }

    int isUnpacked = 0;
    if (c->unpacked) {
    	isUnpacked = 1;
    }
    for (int i = 0; i < nSet; ++i) {
    	if (B[i]->unpacked) {
    		isUnpacked = 1;
    	}
    }
    if (isUnpacked) {
    	fprintf(stderr, "normalizedTriangularSetPseudoDivide_AAZ not yet implemented with exponent unpacking!\n");
    	exit (EXIT_FAILURE);
    }

    AltArrZ_t* tmpQ = NULL;
    AltArrZ_t* tmpR = NULL;
    AltArrZ_t* h = NULL;

	mpz_t one;
	mpz_init (one); // was tmpCoef
	mpz_set_si (one, 1l);

    int e = 1;
  	int mvar;
  	int signMvar; //?
    /* int* expOffset = getExpOffsetArrayZ (nvar); // ? */

    if (nSet == 1) {
        mvar = leadingVariable_AAZ (B[0]);

		if (mvar == -2) {
			*rem = deepCopyPolynomial_AAZ (c);
            quoSet[0] = NULL;
            *hPow = makeConstPolynomial_AAZ (1, nvar, one);
            return;
		}

		if (mvar > -1){
			pesudoDivideAtIdx_AAZ (mvar, c, B[0], &tmpQ, &tmpR, &e, &h, nvar, lazy);
		} else{
			pesudoDivideAtIdx_AAZ (0, c, B[0], &tmpQ, &tmpR, &e, &h, nvar, lazy);
		}

        *rem = tmpR;

		Q = (AltArrZ_t**) calloc (nSet, sizeof (AltArrZ_t*));
		Q[0] = tmpQ;

		*quoSet = Q;
        *hPow = h;
        return;
    }

	Q = (AltArrZ_t**) calloc (nSet, sizeof (AltArrZ_t*));

	/////////////// PART 1 ////////////////
	signMvar = leadingVariable_AAZ (B[nSet-1]);
	mvar = (signMvar < 0) ? 0 : signMvar;

	AltArrZ_t** CoefList1;
	int sz_cl1 = 0; // deg(c, mvar)+1
	mainCoefficientListAtIdx_AAZ (c, mvar, &CoefList1, &sz_cl1);

	if (!sz_cl1) {
		fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AA![Part 1]\n");
		exit (EXIT_FAILURE);
	}

    AltArrZ_t** Q1[sz_cl1];

	AltArrZ_t** R1 = (AltArrZ_t**) calloc (sz_cl1, sizeof (AltArrZ_t*));
	AltArrZ_t** HPow1 = (AltArrZ_t**) calloc (sz_cl1, sizeof (AltArrZ_t*));

	for (int i = 0; i < sz_cl1; i++) {
	    normalizedTriangularSetPseudoDivide_AAZ (CoefList1[i], B,
										   &Q1[i], &R1[i], &HPow1[i], nvar, lazy, nSet-1);
	}

	for (int i = 0; i < sz_cl1; i++) {
		if (CoefList1[i] != NULL && CoefList1[i]->size != 0) {
			freePolynomial_AAZ (CoefList1[i]);
		}
	}
	free (CoefList1);

	AltArrZ_t* H1 = NULL;

	for (int i = 0; i < sz_cl1; i++){
		if (HPow1[i] != NULL && HPow1[i]->size != 0){ // TODO
			if (H1 == NULL || H1->size == 0){
				H1 = deepCopyPolynomial_AAZ (HPow1[i]);
			} else {
				H1 = maxPolynomials_AAZ_inp (H1, HPow1[i]); // TODO: change to lcm()
			}
		}
	}

	AltArrZ_t* newR1 = NULL;
	AltArrZ_t* frac_r1;
	AltArrZ_t* frac_q1;

	for (int i = 0; i < sz_cl1; ++i){

		if (H1 == NULL || H1->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AAZ, H1 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}

		dividePolynomials_AAZ (H1, HPow1[i], &frac_q1, &frac_r1, nvar);

		if (frac_r1 != NULL || frac_r1->size != 0) {
			freePolynomial_AAZ (frac_r1);
		}

		if (frac_q1 == NULL || frac_q1->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AAZ, frac_q1 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}

		if (R1[i] != NULL && R1[i]->size != 0) {

			multiplyPolynomialAtIdxByXn_AAZ_inp (R1[i], mvar, i, nvar);
			R1[i] = multiplyPolynomials_AAZ_inp (R1[i], frac_q1, nvar);
			newR1 = addPolynomials_AAZ_inp (newR1, R1[i], nvar);

		}

		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AAZ_inp (Q1[i][j], mvar, i, nvar);
			Q1[i][j] = multiplyPolynomials_AAZ_inp (Q1[i][j], frac_q1, nvar);
			if (i == 0) {
				Q[j] = deepCopyPolynomial_AAZ (Q1[i][j]);
			} else {
				Q[j] = addPolynomials_AAZ_inp (Q[j], Q1[i][j], nvar);
			}
		}

		freePolynomial_AAZ (frac_q1); frac_q1 = NULL;
	}

	for (int i = 0; i < sz_cl1; ++i){
		if (R1[i] != NULL || R1[i]->size != 0) {
			freePolynomial_AAZ (R1[i]);
		}
		if (HPow1[i] != NULL || HPow1[i]->size != 0) {
			freePolynomial_AAZ (HPow1[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q1[i][j] != NULL || Q1[i][j]->size != 0) {
				freePolynomial_AAZ (Q1[i][j]);
			}
		}
	}

	free (R1);
	// free (HPow1);
	// free (Q1);

	/////////////// PART 2 ////////////////

	if (newR1 == NULL || newR1->size == 0) {
		Q[nSet-1] = NULL;

		*quoSet = Q;
		*rem = NULL;
		*hPow = H1;

		return;
	}

	// Cover the crazy corner cases
	if (!isOne_AAZ (H1) && mainVariable_AAZ (newR1) > 1) {
		for (int i = 0; i < nSet-1; i++) {
			if (partialDegree_AAZ (newR1, mainVariable_AAZ (B[i])) >=
				partialDegree_AAZ (B[i], mainVariable_AAZ (B[i]))) {
				AltArrZ_t* r_star = NULL;
				AltArrZ_t* q_star = NULL;

				dividePolynomials_AAZ (newR1, B[i], &q_star, &r_star, nvar);

				if (q_star != NULL) {
						Q[i] = addPolynomials_AAZ_inp (Q[i], q_star, nvar);
						freePolynomial_AAZ (q_star);
				}
				if (r_star != NULL) {
					freePolynomial_AAZ (newR1);
					newR1 = r_star;
				}
				if (newR1 == NULL || newR1->size == 0) {
					Q[nSet-1] = NULL;

					*quoSet = Q;
					*rem = NULL;
					*hPow = H1;

					return;
				}
			}
		}
	}

	AltArrZ_t* q_tilde = NULL;
	AltArrZ_t* r_tilde = NULL;
	AltArrZ_t* h_tilde = NULL;
	int       e_tilde = 0;

	pesudoDivideAtIdx_AAZ (mvar, newR1, B[nSet-1], &q_tilde, &r_tilde, &e_tilde, &h_tilde, nvar, lazy);

	freePolynomial_AAZ (newR1);

	if (h_tilde == NULL || h_tilde->size == 0) {
		fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AAZ, h_tilde couldn't be zero!\n");
		exit (EXIT_FAILURE);
	}

	if (!isOne_AAZ (h_tilde)) {
		for (int j = 0; j < nSet-1; j++) {
			if (Q[j] != NULL && Q[j]->size != 0) {
				Q[j] = multiplyPolynomials_AAZ_inp (Q[j], h_tilde, nvar);
			}
		}
	}

	Q[nSet-1] = q_tilde;

	if (!isOne_AAZ (h_tilde)) {
		H1 = multiplyPolynomials_AAZ_inp (H1, h_tilde, nvar);
	}

	freePolynomial_AAZ (h_tilde);

	/////////////// PART 3 ////////////////

	if (r_tilde == NULL || r_tilde->size == 0) {
		*quoSet = Q;
		*rem = NULL;
		*hPow = H1;

		return;
	}

	AltArrZ_t** CoefList2;
	int sz_cl2 = 0; // deg(r_tilde, mvar)+1
	mainCoefficientListAtIdx_AAZ (r_tilde, mvar, &CoefList2, &sz_cl2);

	if (!sz_cl2) {
		fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AAZ, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AAZ [Part 3]!\n");
		exit (EXIT_FAILURE);
	}

    AltArrZ_t** Q2[sz_cl2];

	AltArrZ_t** R2 = (AltArrZ_t**) calloc (sz_cl2, sizeof (AltArrZ_t*));
	AltArrZ_t** HPow2 = (AltArrZ_t**) calloc (sz_cl2, sizeof (AltArrZ_t*));

	for (int i = 0; i < sz_cl2; i++) {
	    normalizedTriangularSetPseudoDivide_AAZ (CoefList2[i], B,
										   &Q2[i], &R2[i], &HPow2[i], nvar, lazy, nSet-1);
	}

	for (int i = 0; i < sz_cl2; i++) {
		if (CoefList2[i] != NULL && CoefList2[i]->size != 0) {
			freePolynomial_AAZ (CoefList2[i]);
		}
	}
	free (CoefList2);

	AltArrZ_t* H2 = NULL;

	for (int i = 0; i < sz_cl2; i++){
		if (HPow2[i] != NULL && HPow2[i]->size != 0){ // TODO
			if (H2 == NULL || H2->size == 0){
				H2 = deepCopyPolynomial_AAZ (HPow2[i]);
			} else {
				H2 = maxPolynomials_AAZ_inp (H2, HPow2[i]); // TODO: change to lcm()
			}
		}
	}

	///// Update Quotients:
	if (!isOne_AAZ (H2)) {
		for (int j = 0; j < nSet; j++) {
			if (Q[j] != NULL && Q[j]->size != 0) {
				Q[j] = multiplyPolynomials_AAZ_inp (Q[j], H2, nvar);
			}
		}
	}

	AltArrZ_t* newR2 = NULL;
	AltArrZ_t* frac_r2;
	AltArrZ_t* frac_q2;

	for (int i = 0; i < sz_cl2; ++i){

		if (H2 == NULL || H2->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, H2 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}

		dividePolynomials_AAZ (H2, HPow2[i], &frac_q2, &frac_r2, nvar);

		if (frac_r2 != NULL && frac_r2->size != 0) {
			freePolynomial_AAZ (frac_r2);
		}

		if (frac_q2 == NULL || frac_q2->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, frac_q2 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}

		if (R2[i] != NULL && R2[i]->size != 0) {

			multiplyPolynomialAtIdxByXn_AAZ_inp (R2[i], mvar, i, nvar);
			R2[i] = multiplyPolynomials_AAZ_inp (R2[i], frac_q2, nvar);
			newR2 = addPolynomials_AAZ_inp (newR2, R2[i], nvar);

		}

		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AAZ_inp (Q2[i][j], mvar, i, nvar);
			Q2[i][j] = multiplyPolynomials_AAZ_inp (Q2[i][j], frac_q2, nvar);
			Q[j] = addPolynomials_AAZ_inp (Q[j], Q2[i][j], nvar);
		}

		freePolynomial_AAZ (frac_q2); frac_q2 = NULL;
	}

	for (int i = 0; i < sz_cl2; ++i){
		if (R2[i] != NULL || R2[i]->size != 0) {
			freePolynomial_AAZ (R2[i]);
		}
		if (HPow2[i] != NULL || HPow2[i]->size != 0) {
			freePolynomial_AAZ (HPow2[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q2[i][j] != NULL || Q2[i][j]->size != 0) {
				freePolynomial_AAZ (Q2[i][j]);
			}
		}
	}

	free (R2);
	// free (HPow2);
	// free (Q2);

	*quoSet = Q;
	*rem = newR2;
	*hPow = multiplyPolynomials_AAZ_inp (H1, H2, nvar);

	freePolynomial_AAZ (H2);
}


// TODO:
void multiDivisorPseudoDivide_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t** quoSet, AltArrZ_t** rem, AltArrZ_t** hPow, int nvar, int lazy, int nSet)
{
	if (nSet  < 0 || nvar < nSet){
		fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
		exit(EXIT_FAILURE);
    }


    if (nSet == 0){
		mpz_t one;
		mpz_init (one);
		mpz_set_si (one, 1l);

		*rem = deepCopyPolynomial_AAZ (c);
		if (*hPow != NULL){
			freePolynomial_AAZ (*hPow);
		}
		*hPow = makeConstPolynomial_AAZ (1, nvar, one);
		quoSet[0] = NULL;

		mpz_clear (one);
		return;
    }

	int isNormalized = 0;

	if (isNormalized) {
		AltArrZ_t** Q;

		normalizedTriangularSetPseudoDivide_AAZ (c, B, &Q, rem, hPow, nvar, lazy, nSet);

		for (int i = 0; i < nSet; i++) {
			quoSet[i] = Q[i];
		}

		return;
	}

	naiveMultiDivisorPseudoDivide_AAZ (c, B, quoSet, rem, hPow, nvar, lazy, nSet);
}

void recursiveTriangularSetMDD_AAZ (AltArrZ_t* c, AltArrZ_t** B, AltArrZ_t*** quoSet, AltArrZ_t** rem, int nvar, int nSet)
{
    if (nSet < 0){
        fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
        exit (EXIT_FAILURE);
	}

	AltArrZ_t** Q; // *quoSet = Q

    if (c == NULL || c->size == 0) {
		*rem = NULL;
		if (nSet > 0) {
			Q = (AltArrZ_t**) calloc (nSet, sizeof (AltArrZ_t*));

			for (int i = 0; i < nSet; ++i){
			    Q[i] = NULL;
			}
			*quoSet = Q;

		} else {
			*quoSet = NULL;
		}

        return;
    }

    if (nSet == 0) {
		*rem = deepCopyPolynomial_AAZ (c);
		*quoSet = NULL;

		return;
    }

    int isUnpacked = 0;
    if (c->unpacked) {
    	isUnpacked = 1;
    }
    for (int i = 0; i < nSet; ++i) {
    	if (B[i]->unpacked) {
    		isUnpacked = 1;
    	}
    }
    if (isUnpacked) {
    	fprintf(stderr, "recursiveTriangularSetMDD_AAZ not yet implemented with exponent unpacking!\n");
    	exit (EXIT_FAILURE);
    }

    AltArrZ_t* tmpQ = NULL;
    AltArrZ_t* tmpR = NULL;

	/* mpz_t one; */
	/* mpz_init (one); // was tmpCoef */
	/* mpq_set_si (one); */

  	int mvar;
  	int signMvar;

    if (nSet == 1) {
        mvar = leadingVariable_AAZ (B[0]);

		if (mvar == -2) {
			*rem = deepCopyPolynomial_AAZ (c);
            quoSet[0] = NULL;
            return;
		}

		dividePolynomials_AAZ (c, B[0], &tmpQ, &tmpR, nvar);

        *rem = tmpR;

		Q = (AltArrZ_t**) calloc (nSet, sizeof (AltArrZ_t*));
		Q[0] = tmpQ;
		*quoSet = Q;

        return;
    }

	Q = (AltArrZ_t**) calloc (nSet, sizeof (AltArrZ_t*));

	/* int isLazard = 1; */
	/* if (isLazard) { */
	/* 	/\* lazardTriangularSetMDD_AAZ (c, B, Q, rem, nSet, nvar); *\/ */
	/* 	primitiveFactorTriangularSetMDD_AAZ (c, B, Q, rem, nSet, nvar); */
	/* 	*quoSet = Q; */
	/* 	return; */
	/* } */

	if (nSet < 4 || AA_SIZE (c) < 5) {
		normalForm_AAZ (c, B, Q, rem, nSet, nvar);
		*quoSet = Q;
		return;
	}

	/////////////// PART 1 ////////////////
	signMvar = leadingVariable_AAZ (B[nSet-1]);
	mvar = (signMvar < 0) ? 0 : signMvar;

	AltArrZ_t** CoefList1;
	int sz_cl1 = 0; // deg(c, mvar)+1
	mainCoefficientListAtIdx_AAZ (c, mvar, &CoefList1, &sz_cl1);

	if (!sz_cl1) {
		fprintf (stderr, "SMQP Error: In recursiveTriangularSetMDD_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AA![Part 1]\n");
		exit (EXIT_FAILURE);
	}

    AltArrZ_t** Q1[sz_cl1];

	AltArrZ_t** R1 = (AltArrZ_t**) calloc (sz_cl1, sizeof (AltArrZ_t*));

	for (int i = 0; i < sz_cl1; i++) {
		recursiveTriangularSetMDD_AAZ (CoefList1[i], B, &Q1[i], &R1[i], nvar, nSet-1);
	}

	for (int i = 0; i < sz_cl1; i++) {
		if (CoefList1[i] != NULL && CoefList1[i]->size != 0) {
			freePolynomial_AAZ (CoefList1[i]);
		}
	}
	free (CoefList1);

	AltArrZ_t* newR1 = NULL;

	for (int i = 0; i < sz_cl1; ++i){
		if (R1[i] != NULL && R1[i]->size != 0) {
			multiplyPolynomialAtIdxByXn_AAZ_inp (R1[i], mvar, i, nvar);
			newR1 = addPolynomials_AAZ_inp (newR1, R1[i], nvar);
		}

		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AAZ_inp (Q1[i][j], mvar, i, nvar);
			if (i == 0) {
				Q[j] = deepCopyPolynomial_AAZ (Q1[i][j]);
			} else {
				Q[j] = addPolynomials_AAZ_inp (Q[j], Q1[i][j], nvar);
			}
		}
	}

	for (int i = 0; i < sz_cl1; ++i){
		if (R1[i] != NULL || R1[i]->size != 0) {
			freePolynomial_AAZ (R1[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q1[i][j] != NULL || Q1[i][j]->size != 0) {
				freePolynomial_AAZ (Q1[i][j]);
			}
		}
	}

	free (R1);
	// free (Q1);

	/////////////// PART 2 ////////////////

	if (newR1 == NULL || newR1->size == 0) {
		Q[nSet-1] = NULL;

		*quoSet = Q;
		*rem = NULL;

		return;
	}

	AltArrZ_t* q_tilde = NULL;
	AltArrZ_t* r_tilde = NULL;

	if (monomialDivideTest (newR1->elems[0].degs, B[nSet-1]->elems[0].degs, nvar)) {
		dividePolynomials_AAZ (newR1, B[nSet-1], &q_tilde, &r_tilde, nvar);
		freePolynomial_AAZ (newR1);
	} else {
		q_tilde = NULL;
		r_tilde = newR1;
	}

	Q[nSet-1] = q_tilde;

	/////////////// PART 3 ////////////////

	if (r_tilde == NULL || r_tilde->size == 0) {
		*quoSet = Q;
		*rem = NULL;

		return;
	}

	AltArrZ_t** CoefList2;
	int sz_cl2 = 0; // deg(r_tilde, mvar)+1
	mainCoefficientListAtIdx_AAZ (r_tilde, mvar, &CoefList2, &sz_cl2);

	if (!sz_cl2) {
		fprintf (stderr, "SMQP Error: In recursiveTriangularSetMDD_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AAZ [Part 3]!\n");
		exit (EXIT_FAILURE);
	}

    AltArrZ_t** Q2[sz_cl2];

	AltArrZ_t** R2 = (AltArrZ_t**) calloc (sz_cl2, sizeof (AltArrZ_t*));

	for (int i = 0; i < sz_cl2; i++) {
		recursiveTriangularSetMDD_AAZ (CoefList2[i], B, &Q2[i], &R2[i], nvar, nSet-1);
	}

	for (int i = 0; i < sz_cl2; i++) {
		if (CoefList2[i] != NULL && CoefList2[i]->size != 0) {
			freePolynomial_AAZ (CoefList2[i]);
		}
	}
	free (CoefList2);

	AltArrZ_t* newR2 = NULL;

	for (int i = 0; i < sz_cl2; ++i){

		if (R2[i] != NULL && R2[i]->size != 0) {

			multiplyPolynomialAtIdxByXn_AAZ_inp (R2[i], mvar, i, nvar);
			newR2 = addPolynomials_AAZ_inp (newR2, R2[i], nvar);
		}

		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AAZ_inp (Q2[i][j], mvar, i, nvar);
			Q[j] = addPolynomials_AAZ_inp (Q[j], Q2[i][j], nvar);
		}
	}

	for (int i = 0; i < sz_cl2; ++i){
		if (R2[i] != NULL || R2[i]->size != 0) {
			freePolynomial_AAZ (R2[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q2[i][j] != NULL || Q2[i][j]->size != 0) {
				freePolynomial_AAZ (Q2[i][j]);
			}
		}
	}

	free (R2);
	// free (Q2);

	*quoSet = Q;
	*rem = newR2;
}



//////////////////////////
// Subresultant Chains
//////////////////////////


AltArrZ_t* LazardOptZ (AltArrZ_t* Sd, AltArrZ_t* Sdm, AltArrZ_t* s)
{
    if (Sd == NULL || Sd->size == 0){
		fprintf (stderr, "SMZP Error: In LazardOptZ, Sd is NULL!\n");
		exit (1);
    }
    if (Sdm == NULL || Sdm->size == 0){
		return NULL;
    }
    if (Sd->nvar != Sdm->nvar){
		fprintf (stderr, "SMZP Error: In LazardOptZ, Sd->nvar(=%d) != Sdm->nvar(=%d).\n", Sd->nvar, Sdm->nvar);
		exit (1);
    }
	int nvar = Sd->nvar;
    degree_t n = mainLeadingDegree_AAZ(Sd) - mainLeadingDegree_AAZ(Sdm) - 1;//((Sd->elems[0].degs & mvarMask) >> mvarDegOffset) -
    if (n < 1){
		return deepCopyPolynomial_AAZ (Sdm);
    }

    AltArrZ_t* x = mainLeadingCoefficient_AAZ(Sdm);
    AltArrZ_t* y = s;
    AltArrZ_t* c = deepCopyPolynomial_AAZ (x);
    register degree_t a = 1 << Log2_AAZ (n);
    n = n - a;
    AltArrZ_t *tmpz;
	AltArrZ_t* tmp = makePolynomial_AAZ (c->size * c->size, nvar);

	while (a > 1) {
		tmpz = NULL;
		// tmp = multiplyPolynomials_AAZ_inp (c, c, nvar);
		multiplyPolynomialsPreAlloc_AAZ (c, c, &tmp);
		exactDividePolynomials_AAZ (tmp, y, &tmpz, nvar);
		// freePolynomial_AAZ (tmp);
		c = tmpz;
		a = a >> 1;
		if (n >= a) {
			tmpz = NULL;
			if (c != NULL && c->size != 0) {
				// tmp = multiplyPolynomials_AAZ_inp (c, x, nvar);
				multiplyPolynomialsPreAlloc_AAZ (c, x, &tmp);
				exactDividePolynomials_AAZ (tmp, y, &tmpz, nvar);
			}
			// else {
			// 	tmp = NULL;
			// }
			// exactDividePolynomials_AAZ (tmp, y, &tmpz, nvar);
			// freePolynomial_AAZ (tmp);
			c = tmpz;
			n = n - a;
		}
    }
	freePolynomial_AAZ (tmp);

    tmpz = NULL;
    AltArrZ_t* Se = deepCopyPolynomial_AAZ (Sdm);
    Se = multiplyPolynomials_AAZ_inp (Se, c, nvar);
    freePolynomial_AAZ (c);
    exactDividePolynomials_AAZ (Se, y, &tmpz, nvar);
    freePolynomial_AAZ (Se);
	freePolynomial_AAZ (x); // added 30Jan2020
    return tmpz;
}

AltArrZ_t* DucosOptZ (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd)
{
    int nvar = A->nvar;

    if (nvar != Sdm->nvar || nvar != Se->nvar || nvar != sd->nvar){
	fprintf(stderr, "SMZP ERROR: In DucosOpt, inputs' nvars are different!\n");
	exit(1);
    }

    degree_t d = mainLeadingDegree_AAZ(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t e = mainLeadingDegree_AAZ(Sdm);//(Sdm->elems[0].degs & mvarMask) >> mvarDegOffset;
    AltArrZ_t* cdm = mainLeadingCoefficient_AAZ (Sdm);
    AltArrZ_t* se = mainLeadingCoefficient_AAZ (Se);
    AltArrZ_t* tmpz;
    AltArrZ_t* tmpR;
    AltArrZ_t* tmp = NULL;

    if (d == 0){
	fprintf (stderr, "SMZP Error: In DucosOpt, d == 0.\n");
	exit(1);
    }
    if (e == 0 || d < e){
	fprintf (stderr, "SMZP Error: In DucosOpt, e == 0 || d > e.\n");
	exit(1);
    }
    if (cdm == NULL || cdm->size ==0){
	fprintf (stderr, "SMZP Error: In DucosOpt, cdm = NULL.\n");
	exit(1);
    }
    if (se == NULL || se->size ==0){
	fprintf (stderr, "SMZP Error: In DucosOpt, se = NULL.\n");
	exit(1);
    }
    if (d == e){
	// non-defective cases:
	// Streamlined resultant computation for regular case
	AltArrZ_t* dc = mainLeadingCoefficient_AAZ (A);
	AltArrZ_t* rTemp = NULL;
	AltArrZ_t* supTemp = NULL;

	// compute first of two main terms
	AltArrZ_t* res = deepCopyPolynomial_AAZ (A); // res = a
	res = multiplyPolynomials_AAZ_inp (res, se, nvar); // res *= ec
	rTemp = mainCoefficientAtIdx_AAZ(A, e); // rTemp = a.coef(e)
	supTemp = multiplyPolynomials_AAZ (rTemp, Se, nvar); // supTemp = se * rTerm
	if (supTemp != NULL && supTemp->size != 0){
	    res = subPolynomials_AAZ_inp (res, supTemp, nvar); // res -= supTemp
	}
	dividePolynomials_AAZ (res, dc, &tmp, &tmpR, nvar); // res /= dc
	freePolynomial_AAZ (tmpR);
	res = multiplyPolynomials_AAZ (tmp, se, nvar); // res *= ec
	freePolynomial_AAZ (tmp);

	// compute second of two main terms
	rTemp = mainCoefficientAtIdx_AAZ (Se, e-1); // rTemp = se.coef(e-1)
	freePolynomial_AAZ (supTemp); // supTemp.zero()
	supTemp = rTemp; // supTemp.setCoef(0, rTemp)
	rTemp = deepCopyPolynomial_AAZ (se);
	negatePolynomial_AAZ (rTemp); // rTemp = -ec
	rTemp = mainLShiftPolynomial_AAZ_inp (rTemp, 1);

	supTemp = addPolynomials_AAZ (supTemp, rTemp, nvar); // supTemp.setCoef(1, rTemp)
	supTemp = multiplyPolynomials_AAZ_inp (supTemp, Se, nvar); // supTemp *= se
	res = addPolynomials_AAZ_inp (res, supTemp, nvar); // res += supTemp
	freePolynomial_AAZ (supTemp);
	// divide out dc to obtain resultant
	dividePolynomials_AAZ (res, dc, &tmp, &tmpR, nvar); // res /= dc
	freePolynomial_AAZ (tmpR);
	return tmp; // return res

    }
    // defective cases:
    // H =
    AltArrZ_t* H[d];
    for (int j = 0; j <= e; ++j){
    	H[j] = deepCopyPolynomial_AAZ (se);
    	H[j] = mainLShiftPolynomial_AAZ_inp(H[j], j);
    }

    H[e] = subPolynomials_AAZ_inp (H[e], Se, nvar);

    AltArrZ_t* PIe = NULL;
    for (int j = e+1; j < d; ++j){
	//XH_{j-1} =
	H[j] = deepCopyPolynomial_AAZ (H[j-1]);

	mainLShiftPolynomial_AAZ_inp(H[j], 1);

	//syzygy:
	PIe = mainCoefficientAtIdx_AAZ (H[j], e);

	if (PIe != NULL && PIe->size != 0){
	    tmp = multiplyPolynomials_AAZ_inp (PIe, Sdm, nvar);
	    tmpz = NULL;
	    tmpR = NULL;
	    dividePolynomials_AAZ (tmp, cdm, &tmpz, &tmpR, nvar);
	    freePolynomial_AAZ (tmpR);
	    freePolynomial_AAZ (tmp); //  freePolynomial_AAZ (PIe);
	    H[j] = subPolynomials_AAZ_inp (H[j], tmpz, nvar);
	    freePolynomial_AAZ (tmpz);
	}
    }

    // D =
    AltArrZ_t* PIj = NULL;
    AltArrZ_t* D = NULL;
    for (int j = 0; j < d; ++j){
	PIj = mainCoefficientAtIdx_AAZ (A, j);
	if (PIj != NULL && PIj->size != 0 &&
	    H[j] != NULL && H[j]->size != 0){
	    PIj = multiplyPolynomials_AAZ_inp (PIj, H[j], nvar);

	    if (D == NULL || D->size == 0){
		D = deepCopyPolynomial_AAZ (PIj);
	    } else {
		D = addPolynomials_AAZ_inp (D, PIj, nvar);
	    }
	    freePolynomial_AAZ (PIj);
	}
    }

    tmpz = NULL;
    tmpR = NULL;
    AltArrZ_t* lcA = mainLeadingCoefficient_AAZ (A);
    if (lcA != NULL && lcA->size != 0){
	dividePolynomials_AAZ (D, lcA, &tmpz, &tmpR, nvar);
	freePolynomial_AAZ (tmpR);
	freePolynomial_AAZ (D);
	D = tmpz;
    }
    freePolynomial_AAZ (lcA);

    // res =
    AltArrZ_t* res = NULL;
    AltArrZ_t* lPoly = deepCopyPolynomial_AAZ (H[d-1]);
    lPoly = mainLShiftPolynomial_AAZ_inp(lPoly, 1);

    AltArrZ_t* rPoly = mainCoefficientAtIdx_AAZ (lPoly, e);
    if (D != NULL && D->size != 0){
	lPoly = addPolynomials_AAZ_inp (lPoly, D, nvar);
    }
    lPoly = multiplyPolynomials_AAZ_inp (lPoly, cdm, nvar);

    if (rPoly != NULL && rPoly->size != 0){
	rPoly = multiplyPolynomials_AAZ_inp (rPoly, Sdm, nvar);
	res = subPolynomials_AAZ (lPoly, rPoly, nvar);
	freePolynomial_AAZ (lPoly);
	freePolynomial_AAZ (rPoly);
    } else {
	res = lPoly;
    }

    tmpR = NULL;
    tmpz = NULL;
    if (sd != NULL && sd->size != 0){
	    dividePolynomials_AAZ (res, sd, &tmpz, &tmpR, nvar);
    } else {
	tmpz = res;
    }
    freePolynomial_AAZ (tmpR);

    if ((d-e+1)%2 != 0 && tmpz != NULL && tmpz->size != 0){
	negatePolynomial_AAZ (tmpz);
    }

    freePolynomial_AAZ (cdm);
    freePolynomial_AAZ (se);

    return tmpz;
}

AltArrZ_t* mainTailPolynomial_AAZ (AltArrZ_t* a, AltArrZ_t** lc) {
	if (a == NULL || isZero_AAZ (a)) {
		return NULL;
	}
	degree_t d = mainLeadingDegree_AAZ (a);
	if (!d) {
		*lc = deepCopyPolynomial_AAZ (a);
		return NULL;
	}

	AltArrZ_t* lc_a = mainLeadingCoefficient_AAZ (a);
	AltArrZ_t* lcXd = mainLShiftPolynomial_AAZ (lc_a, d);
	AltArrZ_t* tail = subPolynomials_AAZ (a, lcXd, a->nvar);
	if (!isZero_AAZ (lcXd)) freePolynomial_AAZ (lcXd);

	if (lc == NULL) {
		freePolynomial_AAZ (lc_a);
	} else {
		if (isZero_AAZ (lc_a)) {
			*lc = NULL;
		} else {
			*lc = lc_a;
		}
	}
	if (isZero_AAZ (tail)) {
		return NULL;
	}

	return tail;
}

AltArrZ_t* CFDucosOptZ_new (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd, int flag)
{
    int nvar = A->nvar;
    if (nvar != Sdm->nvar || nvar != Se->nvar || nvar != sd->nvar){
		fprintf (stderr, "SMZP ERROR: In CFDucosOptZ, inputs' nvars are different!\n");
		exit (1);
    }

	degree_t d = mainLeadingDegree_AAZ(A);
    degree_t e = mainLeadingDegree_AAZ(Sdm);
    if (d == 0){
		fprintf (stderr, "SMZP Error: In CFDucosOptZ, d == 0.\n");
		exit (1);
    }
    if (e == 0 || d < e){
		fprintf (stderr, "SMZP Error: In CFDucosOptZ, e == 0 || d > e.\n");
		exit (1);
    }
#ifdef SUBRES_TIME_DEBUG
timer_id id;
timer_time t;
double sum = 0;
id = start_timer ();
#endif
	AltArrZ_t *cdm=NULL, *lc_p=NULL, *se=NULL;
	AltArrZ_t* p = mainTailPolynomial_AAZ (A, &lc_p);
	AltArrZ_t* q = mainTailPolynomial_AAZ (Sdm, &cdm);
	AltArrZ_t* h = mainTailPolynomial_AAZ (Se, &se);
	negatePolynomial_AAZ (h); // h = -h;
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt_new, setup : %lf\n", sum);
	// fprintf (stderr, "\n\n\n\n");
	const char* sym[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST
	// fprintf (stderr, "A := ");
	// printPoly_AAZ (stderr, A, sym, A->nvar);
	// fprintf (stderr, "\np := ");
	// printPoly_AAZ (stderr, p, sym, A->nvar);
	// fprintf (stderr, "\nlc_p := ");
	// printPoly_AAZ (stderr, lc_p, sym, lc_p->nvar);
	// fprintf (stderr, "\n");

	// fprintf (stderr, "Sdm := ");
	// printPoly_AAZ (stderr, Sdm, sym, A->nvar);
	// fprintf (stderr, "\nq := ");
	// printPoly_AAZ (stderr, q, sym, A->nvar);
	// fprintf (stderr, "\ncdm := ");
	// printPoly_AAZ (stderr, cdm, sym, lc_p->nvar);
	// fprintf (stderr, "\n");

	// fprintf (stderr, "Se := ");
	// printPoly_AAZ (stderr, Se, sym, A->nvar);
	// fprintf (stderr, "\nh := ");
	// printPoly_AAZ (stderr, h, sym, A->nvar);
	// fprintf (stderr, "\nse := ");
	// printPoly_AAZ (stderr, se, sym, lc_p->nvar);
	// fprintf (stderr, "\n");
#endif
#ifdef SUBRES_TIME_DEBUG
sum = 0;
id = start_timer ();
#endif
	RecArrZ_t* rP = NULL;
	AltArrZ_t tmpCoefA;
	int activeExp[d+1];
	memset (activeExp, -1, (d+1)*sizeof(int));
	if (!isZero_AAZ (p)) {
		rP = convertToRecursiveArrayZ (p);
		// fprintf (stderr, "rP->size = %d\n", rP->size);
		for (int i = 0; i < rP->size; i++) {
			activeExp[rP->elems[i].exp] = i;
			// fprintf (stderr, "activeExp[%d] = %d;\n", rP->elems[i].exp, i);
		}
		tmpCoefA.nvar = nvar;
		tmpCoefA.unpacked = p->unpacked;
	}
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt_new, conversion : %lf\n", sum);
	// for (int i = 0; i < rP->size; i++) {
	// 	fprintf (stderr, "aExp[%d]=%d\t", i, activeExp[i]);
	// }
	// fprintf (stderr, "\n");
#endif
	AltArrZ_t* a, *quo, *t1, *t2;
	// a = coef(p,e)*h
	if (!isZero_AAZ (h) && activeExp[e] != -1) {
		tmpCoefA.size = tmpCoefA.alloc = rP->elems[activeExp[e]].coefSize;
		tmpCoefA.elems = rP->elems[activeExp[e]].coef;
		a = multiplyPolynomials_AAZ (&tmpCoefA, h, nvar);
	} else {
		a = NULL;
	}
#ifdef SUBRES_TIME_DEBUG
	// fprintf (stderr, "\na := ");
	// printPoly_AAZ (stderr, a, sym, lc_p->nvar);
	// fprintf (stderr, "\n");
sum = 0;
id = start_timer ();
#endif
	for (int i = e+1; i < d; i++) {
		if (!isZero_AAZ (h) && mainLeadingDegree_AAZ (h) == e-1) {
			// h = x * tail(h) - exQuo (h_lc*q, cdm);
			t1 = mainTailPolynomial_AAZ (h, &t2); 
			freePolynomial_AAZ (h);
			h = mainLShiftPolynomial_AAZ (t1, 1); // x*tail(h)
			if (!isZero_AAZ (q) && !isZero_AAZ (t2))
				multiplyPolynomialsPreAlloc_AAZ (t2, q, &t1); // t1 is lc((h))*q
			else {t1 = NULL;}
			quo = NULL;
			exactDividePolynomials_AAZ (t1, cdm, &quo, nvar);
			if (t2 != NULL && t2->size != 0) {freePolynomial_AAZ (t2);}
			else {t2 = NULL;}
			if (t1 != NULL && t2->size != 0) {freePolynomial_AAZ (t1);}
			else {t1 = NULL;}
			h = subPolynomials_AAZ_inp (h, quo, nvar); // h -= quo;
			if(quo != NULL && quo->size != 0) {freePolynomial_AAZ (quo);}
			else {quo = NULL;}
		} else {
			// h = x*h
			if (isZero_AAZ (h)) {
				h = NULL;
			} else {
				h = mainLShiftPolynomial_AAZ_inp (h, 1);
			}
		}
		// a = lcrP[i] * h + a
		if (activeExp[i] != -1) {
			tmpCoefA.size = tmpCoefA.alloc = rP->elems[activeExp[i]].coefSize;
			tmpCoefA.elems = rP->elems[activeExp[i]].coef;
			t1 = multiplyPolynomials_AAZ (&tmpCoefA, h, nvar);
			if (a == NULL) {
				a = t1; t1 = NULL;
			} else {
				a = addPolynomials_AAZ_inp (a, t1, nvar);
				freePolynomial_AAZ (t1); t1 = NULL;
			}
		}
	}
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt_new, loop(e+1,..,d) : %lf\n", sum);
	// fprintf (stderr, "\na := ");
	// printPoly_AAZ (stderr, a, sym, lc_p->nvar);
	// fprintf (stderr, "\nh := ");
	// printPoly_AAZ (stderr, h, sym, lc_p->nvar);
	// fprintf (stderr, "\n");
#endif
#ifdef SUBRES_TIME_DEBUG
sum = 0;
id = start_timer ();
#endif
	// new_p = \sum_{i=0}^{e-1} p_iX^i
	AltArrZ_t* new_p = NULL;
	if (!isZero_AAZ (p)) {
	for (int i = 0; i < e; i++) {
		if (activeExp[i] != -1) {
			tmpCoefA.size = tmpCoefA.alloc = rP->elems[activeExp[i]].coefSize;
			tmpCoefA.elems = rP->elems[activeExp[i]].coef;
			t1 = mainLShiftPolynomial_AAZ (&tmpCoefA, i);
			if (new_p == NULL) {
				new_p = t1; t1 = NULL;
			} else {
				new_p = addPolynomials_AAZ_inp (new_p, t1, nvar);
				freePolynomial_AAZ (t1); // TODO:
			}
		}
	}
	}
	// a += se*new_p
	t1 = multiplyPolynomials_AAZ (se, new_p, nvar);
	if (a != NULL) {
		a = addPolynomials_AAZ_inp (a, t1, nvar);
		freePolynomial_AAZ (t1); t1 = NULL;
	} else {
		a = t1; t1 = NULL;
	}
	exactDividePolynomials_AAZ (a, lc_p, &t1, nvar);
	freePolynomial_AAZ (a);
	a = t1; t1 = NULL;
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt_new, D-poly: %lf\n", sum);
sum = 0;
id = start_timer ();
#endif
	AltArrZ_t *lPoly, *rPoly = NULL;
	if (!isZero_AAZ (h) && mainLeadingDegree_AAZ (h) == e-1) {
		// a = cdm*(x*tail(h)+a) - lc(h)*q
		lPoly = mainTailPolynomial_AAZ (h, &rPoly);
		lPoly = mainLShiftPolynomial_AAZ_inp (lPoly, 1);

		lPoly = addPolynomials_AAZ_inp (lPoly, a, nvar);
		t2 = multiplyPolynomials_AAZ (lPoly, cdm, nvar);
		if (isZero_AAZ(lPoly)) {
			freePolynomial_AAZ (lPoly);
			lPoly = makeConstIntPolynomial_AAZ (1, nvar, 1l);
		}
		multiplyPolynomialsPreAlloc_AAZ (rPoly, q, &lPoly);
		t2 = subPolynomials_AAZ_inp (t2, lPoly, nvar);
		freePolynomial_AAZ (lPoly);
		freePolynomial_AAZ (rPoly);
		lPoly = t2; t2 = NULL;

		// a = cdm(x*tail(h)) + cdm*a - lc(h)*q // Doesn't help!
		// lPoly = mainTailPolynomial_AAZ (h, &rPoly);
		// lPoly = mainLShiftPolynomial_AAZ_inp (lPoly, 1);
		// lPoly = multiplyPolynomials_AAZ_inp (lPoly, cdm, nvar); // cdm(x*tail(h))
		// a = multiplyPolynomials_AAZ_inp (a, cdm, nvar); // cdm*a
		// a = addPolynomials_AAZ_inp (a, lPoly, nvar);
		// multiplyPolynomialsPreAlloc_AAZ (rPoly, q, &lPoly);
		// a = subPolynomials_AAZ_inp (a, lPoly, nvar);
		// freePolynomial_AAZ (lPoly);
		// freePolynomial_AAZ (rPoly);
		// lPoly = a; a = NULL;
	} else {
		// a += c*x*h
		lPoly = mainLShiftPolynomial_AAZ (h, 1);
		lPoly = addPolynomials_AAZ_inp (lPoly, a, nvar);
		lPoly = multiplyPolynomials_AAZ_inp (lPoly, cdm, nvar);
	}
    if (sd != NULL){
		exactDividePolynomials_AAZ (lPoly, sd, &t1, nvar);
		freePolynomial_AAZ (lPoly);
    } else {
		t1 = lPoly;
    }
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt_new, last-part: %lf\n", sum);
#endif
	freePolynomial_AAZ (cdm);
	freePolynomial_AAZ (a);
	freePolynomial_AAZ (h);
	freePolynomial_AAZ (new_p);
	freePolynomial_AAZ (lc_p);
	freePolynomial_AAZ (q);
	freeRecArrayZAndCoef (rP);

	if ((d-e+1)&1 && t1 != NULL) {
		negatePolynomial_AAZ (t1);
	}
	return t1;

}


AltArrZ_t* CFDucosOptZ (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd, int flag)
{
    int nvar = A->nvar;
    if (nvar != Sdm->nvar || nvar != Se->nvar || nvar != sd->nvar){
		fprintf (stderr, "SMZP ERROR: In CFDucosOptZ, inputs' nvars are different!\n");
		exit (1);
    }
    register degree_t d = mainLeadingDegree_AAZ(A);// (A->elems[0].degs & mvarMask) >> mvarDegOffset;
    register degree_t e = mainLeadingDegree_AAZ(Sdm);//(Sdm->elems[0].degs & mvarMask) >> mvarDegOffset;
    if (d == 0){
		fprintf (stderr, "SMZP Error: In CFDucosOptZ, d == 0.\n");
		exit (1);
    }
    if (e == 0 || d < e){
		fprintf (stderr, "SMZP Error: In CFDucosOptZ, e == 0 || d > e.\n");
		exit (1);
    }
    AltArrZ_t* cdm = mainLeadingCoefficient_AAZ (Sdm);
    AltArrZ_t* se = mainLeadingCoefficient_AAZ (Se);
    AltArrZ_t *tmpz=NULL, *tmpR=NULL, *tmp=NULL, *res=NULL;
    AltArrZ_t *PIe=NULL, *PIj=NULL;
    AltArrZ_t *D2=NULL, *D=NULL, *H=NULL;
    if (cdm == NULL || cdm->size ==0){
		fprintf (stderr, "SMZP Error: In CFDucosOptZ, cdm = NULL.\n");
		exit (1);
    }
    if (se == NULL || se->size ==0){
		fprintf (stderr, "SMZP Error: In CFDucosOptZ, se = NULL.\n");
		exit (1);
    }

	if (0) { // not-efficient :)
		// prem (A, -Sdm) / (sd^(d-e)*lc(A))
		tmp = deepCopyPolynomial_AAZ (Sdm);
		negatePolynomial_AAZ (tmp); // -B
		int pdp = 0;
		pesudoDivideAtIdx_AAZ (0, A, tmp, &tmpz, &tmpR, &pdp, &res, nvar, 0); // prem(A,-B)
		freePolynomial_AAZ (tmp); tmp = NULL;
		freePolynomial_AAZ (tmpz);
		freePolynomial_AAZ (res);
		res = exponentiatePoly_AAZ (sd, d-e, nvar);
		tmpz = mainLeadingCoefficient_AAZ (A);
		res = multiplyPolynomials_AAZ_inp (res, tmpz, nvar); // sd^(d-e)*lc(A)
		freePolynomial_AAZ (tmpz); tmpz = NULL;
		dividePolynomials_AAZ (tmpR, res, &tmpz, &tmp, nvar);
		freePolynomial_AAZ (tmp);
		freePolynomial_AAZ (tmpR);
		freePolynomial_AAZ (res);

		freePolynomial_AAZ (cdm);
		freePolynomial_AAZ (se);
		return tmpz;
	}

    if (d == e){
	// for non-defective cases:
	fprintf (stderr, "SMZP Warning: In CFDucosOptZ, d == e.\n");

	/* // Streamlined resultant computation for regular case */
	/* AltArrZ_t* dc = mainLeadingCoefficient_AAZ (A); */
	/* AltArrZ_t* rTemp = NULL; */
	/* AltArrZ_t* supTemp = NULL; */

	/* // compute first of two main terms */
	/* AltArrZ_t* res = deepCopyPolynomial_AAZ (A); // res = a */
	/* res = multiplyPolynomials_AAZ_inp (res, se, nvar); // res *= ec */
	/* rTemp = mainCoefficientAtIdx_AAZ(A, e); // rTemp = a.coef(e) */
	/* supTemp = multiplyPolynomials_AAZ (rTemp, Se, nvar); // supTemp = se * rTerm */
	/* if (supTemp != NULL && supTemp->size != 0){ */
	/*     res = subPolynomials_AAZ_inp (res, supTemp, nvar); // res -= supTemp */
	/* } */
	/* /\* dividePolynomials_AAZ (res, dc, &tmp, &tmpR, nvar); // res /= dc *\/ */
	/* exactDividePolynomials_AAZ (res, dc, &tmp, nvar); // res /= dc */
	/* // freePolynomial_AAZ (tmpR); */
	/* res = multiplyPolynomials_AAZ (tmp, se, nvar); // res *= ec */
	/* freePolynomial_AAZ (tmp); */

	/* // compute second of two main terms */
	/* rTemp = mainCoefficientAtIdx_AAZ (Se, e-1); // rTemp = se.coef(e-1) */
	/* freePolynomial_AAZ (supTemp); // supTemp.zero() */
	/* supTemp = rTemp; // supTemp.setCoef(0, rTemp) */
	/* rTemp = deepCopyPolynomial_AAZ (se); */
	/* negatePolynomial_AAZ (rTemp); // rTemp = -ec */
	/* for (int i = 0; rTemp != NULL && i < rTemp->size; ++i){ */
	/*     rTemp->elems[i].degs += ((degrees_t) 1 << mvarDegOffset); */
	/* } */
	/* supTemp = addPolynomials_AAZ (supTemp, rTemp, nvar); // supTemp.setCoef(1, rTemp) */
	/* supTemp = multiplyPolynomials_AAZ_inp (supTemp, Se, nvar); // supTemp *= se */
	/* res = addPolynomials_AAZ_inp (res, supTemp, nvar); // res += supTemp */
	/* freePolynomial_AAZ (supTemp); */
	/* // divide out dc to obtain resultant */
	/* exactDividePolynomials_AAZ (res, dc, &tmp, nvar); // res /= dc */
	/* return tmp; // return res */
    }

    // H[e] =
    H = mainLShiftPolynomial_AAZ (se, e);
    H = subPolynomials_AAZ_inp (H, Se, nvar);
	AltArrZ_t* tmpProduct = makePolynomial_AAZ (1, nvar);

    // D2 =
    PIe = mainCoefficientAtIdx_AAZ (A, e);
    if (PIe != NULL && H != NULL && PIe->size != 0 && H->size != 0){
		// D2 = multiplyPolynomials_AAZ (H, PIe, nvar);
		// freePolynomial_AAZ (PIe);
		D2 = multiplyPolynomials_AAZ_inp (PIe, H, nvar);
		PIe = NULL;
    }

#ifdef SUBRES_TIME_DEBUG
	timer_id id;
	timer_time t;
  	double sum = 0;
	id = start_timer ();
#endif
    for (int j = e+1; j < d; ++j){
		//  H = XH - pi_e (XH)*C
		H = mainLShiftPolynomial_AAZ_inp (H, 1); // H = XH
		PIe = mainCoefficientAtIdx_AAZ (H, e);	// it must be commented for using CFDucosOptZ_subPolynomials_AAZ_inp
		if (PIe != NULL && PIe->size != 0){
			// tmp = multiplyPolynomials_AAZ_inp (PIe, Sdm, nvar);
			multiplyPolynomialsPreAlloc_AAZ (PIe, Sdm, &tmpProduct);
			tmpz = NULL;
			tmpR = NULL;
			exactDividePolynomials_AAZ (tmpProduct, cdm, &tmpz, nvar);
			// freePolynomial_AAZ (tmp); //  freePolynomial_AAZ (PIe);
			H = subPolynomials_AAZ_inp (H, tmpz, nvar);
			// H = CFDucosOptZ_subPolynomials_AAZ_inp (H, tmpz, nvar, &PIe, e-1);
			freePolynomial_AAZ (tmpz);
		}
		tmp = mainCoefficientAtIdx_AAZ (A, j);
		if (tmp != NULL && H != NULL && tmp->size != 0 && H->size != 0){
			tmp = multiplyPolynomials_AAZ_inp (tmp, H, nvar);
			D2 = addPolynomials_AAZ_inp (D2, tmp, nvar);
			freePolynomial_AAZ (tmp);
		}
    }
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum += (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt, loop(e+1,..,d) : %lf\n", sum);
#endif


    // D =
	AltArrZ_t** cList = NULL;
	int len = 0;
	// fprintf (stderr, "start mclist...\n");
    mainCoefficientListAtIdx_AAZ (A, 0, &cList, &len);
	// fprintf (stderr, "stop... (len=%d, e=%d)\n", len, e);
#ifdef SUBRES_TIME_DEBUG
	sum = 0;
	id = start_timer ();
#endif
	for (int j = 0; j < e; ++j){
		if (j >= len) {
			PIj = NULL;
		} else {
			PIj = cList[j]; // mainCoefficientAtIdx_AAZ (A, j);
		}
		if (PIj != NULL && PIj->size != 0){
			PIj = mainLShiftPolynomial_AAZ_inp (PIj, j);
 			// PIj = multiplyPolynomials_AAZ_inp (PIj, se, nvar);
			 multiplyPolynomialsPreAlloc_AAZ (PIj, se, &tmpProduct);
			 freePolynomial_AAZ (PIj);
			 PIj = tmpProduct;
		}
		if (D == NULL || D->size == 0){
			D = deepCopyPolynomial_AAZ (PIj);
		} else {
			D = addPolynomials_AAZ_inp (D, PIj, nvar);
		}
		freePolynomial_AAZ (PIj);
    }
    // free cList
	if (cList != NULL) {
		free (cList);
	}

#ifdef SUBRES_TIME_DEBUG
	sum = 0;
	id = start_timer ();
#endif
    if (D2 != NULL && D2->size != 0){
		D = addPolynomials_AAZ_inp (D, D2, nvar);
		freePolynomial_AAZ (D2);
    }

    tmpz = NULL;
    tmpR = NULL;
    PIe = mainLeadingCoefficient_AAZ (A);
    if (PIe != NULL && PIe->size != 0){
	dividePolynomials_AAZ (D, PIe, &tmpz, &tmpR, nvar);
	freePolynomial_AAZ (tmpR);
	freePolynomial_AAZ (D);
	D = tmpz;
    }
    freePolynomial_AAZ (PIe);
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum += (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt, making D : %lf\n", sum);
#endif

#ifdef SUBRES_TIME_DEBUG
	sum = 0;
	id = start_timer ();
#endif
    AltArrZ_t* lPoly = mainLShiftPolynomial_AAZ (H, 1);
    AltArrZ_t* rPoly = mainCoefficientAtIdx_AAZ (lPoly, e);
    if (D != NULL && D->size != 0){
		lPoly = addPolynomials_AAZ_inp (lPoly, D, nvar);
    }
    if (cdm != NULL){
		lPoly = multiplyPolynomials_AAZ_inp (lPoly, cdm, nvar);
    }

    if (rPoly != NULL && rPoly->size != 0){
		rPoly = multiplyPolynomials_AAZ_inp (rPoly, Sdm, nvar);
		res = subPolynomials_AAZ_inp (lPoly, rPoly, nvar);
		//freePolynomial_AAZ (lPoly);
		freePolynomial_AAZ (rPoly);
    } else {
		res = lPoly;
    }

    tmpR = NULL;
    tmpz = NULL;
    if (sd != NULL && sd->size != 0){
		dividePolynomials_AAZ (res, sd, &tmpz, &tmpR, nvar);
		freePolynomial_AAZ (tmpR);
		freePolynomial_AAZ (res);
    } else {
		tmpz = res;
    }

	// new implementation:
	// AltArrZ_t *lPoly, *rPoly;
	// if (mainLeadingDegree_AAZ (H) == e-1) {
	// 	rPoly = mainLeadingCoefficient_AAZ (H);
	// 	lPoly = subPolynomials_AAZ (H, rPoly, nvar); // tail = H-lcoef(H)
	// 	lPoly = mainLShiftPolynomial_AAZ_inp (lPoly, 1);
	// 	lPoly = addPolynomials_AAZ_inp (lPoly, D, nvar);
	// 	lPoly = multiplyPolynomials_AAZ_inp (lPoly, cdm, nvar);
	// 	rPoly = multiplyPolynomials_AAZ_inp (rPoly, Sdm, nvar);
	// 	lPoly = subPolynomials_AAZ_inp (lPoly, rPoly, nvar);
	// 	freePolynomial_AAZ (rPoly);
	// } else {
	// 	lPoly = mainLShiftPolynomial_AAZ (H, 1);
	// 	lPoly = addPolynomials_AAZ_inp (lPoly, D, nvar);
	// 	lPoly = multiplyPolynomials_AAZ_inp (lPoly, cdm, nvar);
	// }

	// tmpz = NULL;
    // if (sd != NULL && sd->size != 0){
	// 	exactDividePolynomials_AAZ (lPoly, sd, &tmpz, nvar);
	// 	freePolynomial_AAZ (lPoly);
    // } else {
	// 	tmpz = lPoly;
    // }

	// freePolynomial_AAZ (H);
	// freePolynomial_AAZ (D);

#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum += (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In CFDucosOpt, last-part : %lf\n", sum);
#endif

    if ((d-e+1)%2 != 0 && tmpz != NULL && tmpz->size != 0){
		negatePolynomial_AAZ (tmpz);
    }

    freePolynomial_AAZ (cdm);
    freePolynomial_AAZ (se);

    return tmpz;

}


///Use Monomial Basis
// if the number of variables <= 5
// if the number of variables > 5 and  the max number of terms is in the coef (rec. view in main variable) is greater than 2
// TRDsub_resutant_chain
// kernelopts(opaquemodules=false);
// TRDsub_resultant_chain:= RegularCains:-TRDsub....

void DucosSubresultantChainZ_rev (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len, int type)
{

    if (P == NULL || P->size == 0){
		AltArrsZ_t* SS0 = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
		SS0->poly = NULL;
		SS0->next = NULL;
		*SC = SS0;
		*len = 0;
		return;
    }
    if (Q == NULL || Q->size == 0){
		AltArrsZ_t* SS1 = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
		SS1->poly = deepCopyPolynomial_AAZ (P);
		SS1->next = NULL;
		*SC = SS1;
		*len = 1;
		return;
    }
    if (P->nvar != Q->nvar){
		fprintf (stderr, "SMZP Error: In DucosSubresultantChainZ_rev, P->nvar(%d) != Q->nvar(%d).\n",P->nvar, Q->nvar);
		exit(1);
    }

    degree_t degP = mainLeadingDegree_AAZ(P);
    degree_t degQ = mainLeadingDegree_AAZ(Q);

	// fprintf(stderr, "[DUCOS] nvar = %d, deg(P) = %d, deg(Q) = %d\n", P->nvar, degP, degQ);	// TEST

    // adding corner cases to support algebraic extension
    if (degP == 0 && degQ == 0) {
    	mpz_t one;
    	mpz_init (one);
    	mpz_set_si (one, 1l);
    	AltArrsZ_t* SS2 = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
    	SS2->poly = makeConstPolynomial_AAZ (1, P->nvar, one);
    	SS2->next = NULL;
    	*SC = SS2;
    	*len = 1;
    	return;
    } else if (degP == 0) {
    	AltArrsZ_t* SS3 = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
    	SS3->poly = exponentiatePoly_AAZ (P, degQ, P->nvar);
    	SS3->next = NULL;
    	*SC = SS3;
    	*len = 1;
    	return;
    } else if (degQ == 0) {
    	AltArrsZ_t* SS4 = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
    	SS4->poly = exponentiatePoly_AAZ (Q, degP, Q->nvar);
    	SS4->next = NULL;
    	*SC = SS4;
    	*len = 1;
    	return;
    }

///////////////////////// ////////////// ///////////////////////////

    int lazy = 0; // lazy option in pseudoDivide algorithm
    register int nvar = P->nvar;
    int tmpE = 0;
    AltArrZ_t *tmpz=NULL, *tmpH=NULL, *tmpB=NULL;
    AltArrZ_t *s=NULL, *p=NULL, *q=NULL;
    AltArrZ_t *A=NULL, *B=NULL, *C=NULL;
	AltArrsZ_t *newPoly=NULL, *SS=NULL, *tail=NULL;
	AltArrZ_t *tailm=NULL;
    degree_t delta = 0;
    int size = 0;
    if (degP >= degQ){
	p = P;	q = Q;
    } else {
	p = Q; q = P;
    }
    degP = mainLeadingDegree_AAZ(p);
    degQ = mainLeadingDegree_AAZ(q);

	// int isBiCalled = 0;
// 	AltArrsZ_t* SC2 = NULL;
// 	int len2 = 0;
// 	if (nvar == 2) {
// 		degree_t pdeg_p =  partialDegree_AAZ (p, 1);
// 		degree_t pdeg_q =  partialDegree_AAZ (q, 1);
// 		// && isFF && (MIN (mainLeadingDegree_AAZ (P), mainLeadingDegree_AAZ (Q)) > DBZP_MINDEGY_SUBRES_THRESHOLD || isFF2_val)
// 		if (degP > degQ &&
// 			MIN(degP, degQ) > 10) {
// 		// fprintf (stderr, "[nvar = %d] P := ", p->nvar);
// 		// printPoly_AAZ (stderr, p, sym, p->nvar);
// 		// fprintf (stderr, "\n[nvar= %d] q := ", q->nvar);
// 		// printPoly_AAZ (stderr, q, sym, q->nvar);
// 		// fprintf (stderr, "\n");
// 		// fprintf (stderr, "\n");
// 		// fprintf (stderr, "\n");
// 		if (pdeg_p > pdeg_q && MIN(pdeg_p, pdeg_q) > 8) { // TODO: It's a rough threshold!
// 			// fprintf (stderr, "DBZP: bi-modularSubresultantChain is called!\n");
// 			biModularSubresultantChainZ (p, q, &SC2, &len2);
// 			// isBiCalled =1;
// 			*SC = SC2;
// 			*len = len2;
// 			return;
// 		}
// 		}
// 	}
///////////////////////// ////////////// ///////////////////////////
	s = mainLeadingCoefficient_AAZ (q);
    s = exponentiatePoly_AAZ (s, degP - degQ, nvar);
    A = deepCopyPolynomial_AAZ (q);
	int isFlag = 1;
    // B =
    AltArrZ_t* negQ = deepCopyPolynomial_AAZ (q);
    negatePolynomial_AAZ (negQ);

#ifdef SUBRES_TIME_DEBUG
timer_id id;
timer_time t;
double sum = 0;
id = start_timer ();
#endif
	pesudoDivideAtIdx_AAZ (0, p, negQ, &tmpz, &B, &tmpE, &tmpH, nvar, lazy);
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In DucosSubres, pDiv : %lf\n", sum);
#endif
    freePolynomial_AAZ (tmpz);
    freePolynomial_AAZ (tmpH);
    freePolynomial_AAZ (negQ);

    newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
    newPoly->poly = deepCopyPolynomial_AAZ (p);
    newPoly->next = NULL;
    SS = newPoly;
    tail = newPoly;
    size++;
    newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
    if (B == NULL || B->size == 0){
	newPoly->poly = A;
    } else{
	newPoly->poly = deepCopyPolynomial_AAZ (A);
    }
    newPoly->next = NULL;
    tail->next = newPoly;
    tail = newPoly;
    size++;

    degree_t d=0, e=0;
    while (1){
    	if (B == NULL || B->size == 0){ break; }
    	d = mainLeadingDegree_AAZ(A);
    	e = mainLeadingDegree_AAZ(B);
    	newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
    	newPoly->poly = deepCopyPolynomial_AAZ (B);
    	newPoly->next = NULL;
    	tailm = tail->poly;
    	tail->next = newPoly;
    	tail = newPoly;
    	size++;

    	delta = d - e;
    	if (delta > 1){
#ifdef SUBRES_TIME_DEBUG
sum = 0;
id = start_timer ();
#endif
			C = LazardOptZ (tailm, tail->poly, s); // TODO: add s to the inputs
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In DucosSubres, LazardOptZ : %lf\n", sum);
#endif
			newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			newPoly->poly = deepCopyPolynomial_AAZ (C);
			newPoly->next = NULL;
			tail->next = newPoly;
			tail = newPoly;
			size++;
		} else {
				C = deepCopyPolynomial_AAZ (B);
		}
		if (e == 0){ break; } // B.degree() = 0
		if (type == 0){	
			// tmpB = CFDucosOptZ (A, B, C, s, isFlag);
			tmpB = DucosOptZ(A, B, C, s);
		} else {
#ifdef SUBRES_TIME_DEBUG
sum = 0;
id = start_timer ();
#endif
			tmpB = CFDucosOptZ_new (A, B, C, s, isFlag);
#ifdef SUBRES_TIME_DEBUG
t = elapsed_time (&id);
sum = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
fprintf (stderr, "In DucosSubres, CFDucosOpt : %lf\n", sum);
#endif
		}
		freePolynomial_AAZ (B);
		freePolynomial_AAZ (A);
		B = tmpB;
		A = C; // deepCopyPolynomial_AAZ (C);
		freePolynomial_AAZ (s);
		s = mainLeadingCoefficient_AAZ (A);
		isFlag = 0;
		C = NULL; // freePolynomial_AAZ (C);
	}

    // if subresultant is 0, add it to SS:
    if (tail->poly != NULL && tail->poly->size != 0 && mainLeadingDegree_AAZ(tail->poly) > 0) {
		newPoly = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		newPoly->poly = NULL;
		newPoly->next = NULL;
		tail->next = newPoly;
		tail = newPoly;
		size++;
    }
    *len = size;
    *SC = SS;

    return;
}

void DucosSubresultantChainZ (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len)
{
    AltArrsZ_t* invSC;
    int size = 0;
	SubresultantChainWrapperZ_rev (AUTO, P, Q, &invSC, &size);
    // reverse the subresultantchain:
    if (size > 1){
		AltArrsZ_t* curr = invSC, *prev = NULL, *next = NULL;
		while (curr != NULL){
			next = curr->next; curr->next = prev;
			prev = curr; curr = next;
		}
		*SC = prev;
		*len = size;
		return;
    }
    *SC = invSC;
    *len = size;
    return;
}

void DucosSubresultantChainAtIdxZ (AltArrZ_t* P, AltArrZ_t* Q, int idx, AltArrZ_t** SC_idx, AltArrZ_t** SC_idx1)
{
	if (idx < 0) {
		*SC_idx = NULL;
		*SC_idx1 = NULL;
		return;
	}


	AltArrsZ_t *SS;
	int size = 0;

	SubresultantChainMode mode = ChooseTheBestSubresultantChainModeZ(P, Q);
	int degP = mainLeadingDegree_AAZ(P);
	int degQ = mainLeadingDegree_AAZ(Q);
	if (mode == ModUnivar && degP > degQ) {
		DUZP_t* uP = convertFromAltArrZ_DUZP (P);
		DUZP_t* uQ = convertFromAltArrZ_DUZP (Q);
		mpz_t uBound;
		mpz_init (uBound);
		int chain_size = 0;
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** HGCD-ModUnivar SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d, k = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q), idx);
#endif
		DUZP_t** usubres = modularSubresultantAtDeg_DUZP (uP, uQ, idx, &chain_size, uBound, 0, 1) ; // using heuristic &&  hgcd
		// char sym[1] = {'x'}; // TEST
		// fprintf (stderr, "uP := ");
		// printPoly_DUZP (uP, sym);
		// fprintf (stderr, "uQ := ");
		// printPoly_DUZP (uQ, sym);
		// for (int i = 0; i < chain_size; i++) { // TEST
		// 	printPoly_DUZP (usubres[i], sym);
		// }
		freePolynomial_DUZP (uP);
		freePolynomial_DUZP (uQ);
		mpz_clear (uBound);

		if (chain_size==2) {
			*SC_idx1 = convertToAltArrZ_DUZP (usubres[1]);
			*SC_idx = convertToAltArrZ_DUZP (usubres[0]);
			// char sym[1] = {'x'}; // TEST
			// fprintf (stderr, "SC_idx1 := ");
			// printPoly_DUZP (usubres[1], sym);
			// fprintf (stderr, "SC_idx := ");
			// printPoly_DUZP (usubres[0], sym);
		} else if (chain_size==1) {
			if (usubres[0]->lt == idx) {
				*SC_idx = convertToAltArrZ_DUZP (usubres[0]);
				*SC_idx1 = NULL;
			} else {
				*SC_idx1 =  convertToAltArrZ_DUZP (usubres[0]);
				*SC_idx = NULL;
			}
		} else {
			*SC_idx = NULL;
			*SC_idx1 = NULL;
		}
		for (int i = 0; i < chain_size; i++) {
			if (usubres[i] != NULL) { freePolynomial_DUZP (usubres[i]); }
		}
		if (usubres != NULL) { free (usubres); }

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** SRC... DONE\n");
#endif
		return;
	} else if (mode == ModBivarFFT && degP > degQ) {
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** HGCD-ModBivarFFT SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d, k = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q), idx);
#endif
		// int fft_flag = hgcdBiModularFFTSubresultantChainZ(P, Q, idx, &SS, &size);
		int fft_flag = hgcdBiModularFFT4SubresultantChainZ(P, Q, idx, &SS, &size);
		if (fft_flag) {
			if (size == 2) {
				*SC_idx1 = deepCopyPolynomial_AAZ(SS->poly);
				*SC_idx = deepCopyPolynomial_AAZ(SS->next->poly);
			} else if (size == 1) {
				*SC_idx1 = deepCopyPolynomial_AAZ(SS->poly);
				*SC_idx = NULL;
			} else {
				*SC_idx1 = NULL;
				*SC_idx = NULL;
			}
			freeAltArrsZ (SS);
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** SRC... DONE\n");
#endif
			return;
		}
	}

	// char* ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */
	//  fprintf (stderr, "P := ");
	//  printPoly_AAZ_unpk(stderr, P, ch, P->nvar);
	//  fprintf (stderr, "\n");
	//  fprintf (stderr, "Q := ");
	//  printPoly_AAZ_unpk(stderr, Q, ch, Q->nvar);
	//  fprintf (stderr, "\n");
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** HGCD-DUCOS SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d, k = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q), idx);
#endif
	// DucosSubresultantChainZ_rev (P, Q, &SS, &size, 1);
	SubresultantChainWrapperZ_rev (AUTO, P, Q, &SS, &size);

	// fprintf(stderr, "lt = %d |||| idx = %d ||| size = %d\n", mainLeadingDegree_AAZ (SS->poly), idx, size);
	if (size) {
		AltArrsZ_t* cur = SS;
		AltArrsZ_t* curP = cur;
		while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) > idx) {
			// fprintf (stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly)); // TEST
			curP = cur;
			cur = cur->next;
		}
		if (cur != NULL) {
			// it's not end of the chain:
			if (mainLeadingDegree_AAZ (cur->poly) == idx) {
				if (cur->next == NULL || mainLeadingDegree_AAZ (cur->next->poly) < idx) {
					*SC_idx = deepCopyPolynomial_AAZ (cur->poly);
				} else {
					while (cur->next != NULL && mainLeadingDegree_AAZ (cur->next->poly) == idx) {
						cur = cur->next;
					}
					*SC_idx = deepCopyPolynomial_AAZ (cur->poly);
				}
			} else {
				*SC_idx = NULL;
			}
			*SC_idx1 = deepCopyPolynomial_AAZ (curP->poly);
		} else {
			*SC_idx = NULL;
			*SC_idx1 = deepCopyPolynomial_AAZ (curP->poly);
		}
		freeAltArrsZ (SS);
	} else {
		*SC_idx = NULL;
		*SC_idx1 = NULL;
	}
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** SRC... DONE\n");
#endif
}

// for the defective case it always gives you the upper polynoomials in the chain.
void DucosSubresultantChainAtIdxZ_withPrincipleCoefs (AltArrZ_t* P, AltArrZ_t* Q, int idx, AltArrZ_t** SC_idx, AltArrZ_t** SC_idx1, AltArrsZ_t** principle_coefs, int* pcSize)
{
	if (idx < 0) {
		*SC_idx = NULL;
		*SC_idx1 = NULL;
		return;
	}

	AltArrsZ_t* SS;
	int size = 0;
	SubresultantChainWrapperZ_rev (AUTO, P, Q, &SS, &size);

	// fprintf(stderr, "lt = %d |||| idx = %d ||| size = %d\n", mainLeadingDegree_AAZ (SS->poly), idx, size);
	if (size) {
		AltArrsZ_t *node, *head = NULL, *tail = NULL;
		int pcIdx = 0;
		AltArrsZ_t* cur = SS;
		AltArrsZ_t* curP = cur;
		int curDeg = 0;
		int prevDeg = mainLeadingDegree_AAZ (cur->poly);
		while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) > idx) {
			// fprintf (stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly)); // TEST
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			// fprintf (stderr, "prevDeg-curDeg = %d\n", prevDeg-curDeg);
			if (pcIdx > 0 && curDeg == prevDeg) {
				curP = cur;
				// prevDeg = curDeg;
				cur = cur->next;
				continue;
			}
			if (prevDeg-curDeg > 1) {
				for (int i = 0; i < prevDeg-curDeg-1; i++) {
					// fprintf (stderr, "NULL poly\n");
					node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
					node->poly =  NULL;
					node->next = NULL;
					if (head == NULL) {
						head = node;
						tail = node;
					} else {
						tail->next = node;
						tail = node;
					}
					pcIdx++;
				}
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly =  mainLeadingCoefficient_AAZ (cur->poly);
			node->next = NULL;
			if (head == NULL) {
				head = node;
				tail = node;
			} else {
				tail->next = node;
				tail = node;
			}
			pcIdx++;

			curP = cur;
			prevDeg = curDeg;
			cur = cur->next;
		}
		if (cur != NULL) {
			// it's not end of the chain:
			if (mainLeadingDegree_AAZ (cur->poly) == idx) {
				if (cur->next == NULL || mainLeadingDegree_AAZ (cur->next->poly) < idx) {
					*SC_idx = deepCopyPolynomial_AAZ (cur->poly);
				} else {
					while (cur->next != NULL && mainLeadingDegree_AAZ (cur->next->poly) == idx) {
						cur = cur->next;
					}
					*SC_idx = deepCopyPolynomial_AAZ (cur->poly);
				}
			} else {
				*SC_idx = NULL;
			}
			*SC_idx1 = deepCopyPolynomial_AAZ (curP->poly);
		} else {
			*SC_idx = NULL;
			*SC_idx1 = deepCopyPolynomial_AAZ (curP->poly);
		}
		freeAltArrsZ (SS);
		if (pcIdx > 1) {
			AltArrsZ_t* curr = head;
			AltArrsZ_t* prev = NULL;
			AltArrsZ_t* next = NULL;
				while (curr != NULL){
					next = curr->next;
					curr->next = prev;
					prev = curr;
					curr = next;
				}
			*principle_coefs = prev;
		} else {
			*principle_coefs = head;
		}
		*pcSize = pcIdx;
	} else {
		*SC_idx = NULL;
		*SC_idx1 = NULL;
	}
}

AltArrZ_t* DucosResultantZ (AltArrZ_t* P, AltArrZ_t* Q)
{
	if (isZero_AAZ(P)) {
		return deepCopyPolynomial_AAZ(Q);
	}
	if (isZero_AAZ(Q)) {
		return deepCopyPolynomial_AAZ(P);
	}

	AltArrZ_t *sr_idx = NULL, *sr_idx1 = NULL;
	SubresultantChainAtIdxZ (AUTO, P, Q, 0, &sr_idx, &sr_idx1, NULL);
	freePolynomial_AAZ (sr_idx1);
	return sr_idx;
// 	AltArrZ_t* res = NULL;
//     if (P->nvar == 1 && Q->nvar == 1) {
//         DUZP_t* uP = convertFromAltArrZ_DUZP (P);
//         DUZP_t* uQ = convertFromAltArrZ_DUZP (Q);
//         mpz_t uBound;
//         mpz_init_set_si(uBound, -1);
// #if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
// 	fprintf(stderr, "*** ModUnivar Resultant [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
// #endif
//         DUZP_t* ures = modularResultant_DUZP (uP, uQ, uBound, 0, 1); // using hgcd && heuristic
//         freePolynomial_DUZP (uP);
//         freePolynomial_DUZP (uQ);
//         mpz_clear (uBound);
//         res = convertToAltArrZ_DUZP (ures);
//     } /*else if nvar == 2 {

//     } */ else {
//         AltArrZ_t* idx1 = NULL;
// 		DucosSubresultantChainAtIdxZ (P, Q, 0, &res, &idx1);
// 		freePolynomial_AAZ(idx1);
//     }
  	// return res;
}


AltArrZ_t* DucosGCDZ (AltArrZ_t* P, AltArrZ_t* Q)
{
    if (P == NULL || P->size == 0){
    	return deepCopyPolynomial_AAZ (Q);
    }
    if (Q == NULL || Q->size == 0){
    	return deepCopyPolynomial_AAZ (P);
    }
    if (P->nvar != Q->nvar){
    	fprintf (stderr, "SMZP Error: In DucosGCDZ, P->nvar(=%d) != Q->nvar(=%d).\n",
    		 P->nvar, Q->nvar);
    	exit (1);
    }

    int nvar = P->nvar;

    mpz_t one;
    mpz_init(one);
    mpz_set_si (one, 1l);

    degree_t degP = mainLeadingDegree_AAZ(P);//(P->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t degQ = mainLeadingDegree_AAZ(Q);//(Q->elems[0].degs & mvarMask) >> mvarDegOffset;

    if (degP == 0 || degQ == 0){
    	return makeConstPolynomial_AAZ (1, nvar, one);
    }

    AltArrZ_t* res;
    AltArrsZ_t* SC;
    int size = 0;
    DucosSubresultantChainZ (P, Q, &SC, &size);
    AltArrsZ_t* cur = SC;

    if (cur->poly == NULL || cur->poly->size == 0){
    	freeAltArrsZ (SC);
		return makeConstPolynomial_AAZ (1, nvar, one);
    }

    degree_t nd;
    degree_t nnd;
    if (size > 2){
    	nd = mainLeadingDegree_AAZ(P);//(cur->poly->elems[0].degs & mvarMask) >> mvarDegOffset;
    	nnd = mainLeadingDegree_AAZ(cur->next->poly);//(cur->next->poly->elems[0].degs & mvarMask) >> mvarDegOffset;

    	while ( nd == nnd ){
    		if (cur->next == NULL || cur->next->poly == NULL){
    			break;
    		}
    		cur = cur->next;
    		nnd = mainLeadingDegree_AAZ(cur->poly);//(cur->poly->elems[0].degs & mvarMask) >> mvarDegOffset;
    	}
    	if (nd == 0){
    		freeAltArrsZ (SC);
    		return makeConstPolynomial_AAZ (1, nvar, one);
    	} else {
    		res = deepCopyPolynomial_AAZ (cur->poly);
    		freeAltArrsZ (SC);
    		return res;
    	}
    }

    res = deepCopyPolynomial_AAZ (cur->next->poly);
    freeAltArrsZ (SC);

    return res;
}

AltArrZ_t* lastNonZeroChain_AAZ (AltArrZ_t* P, AltArrZ_t* Q)
{
    if (P == NULL || P->size == 0){
    	return deepCopyPolynomial_AAZ (Q);
    }
    if (Q == NULL || Q->size == 0){
    	return deepCopyPolynomial_AAZ (P);
    }
    if (P->nvar != Q->nvar){
    	fprintf (stderr, "SMZP Error: In lastNonZeroChain, P->nvar(=%d) != Q->nvar(=%d).\n",
    		 P->nvar, Q->nvar);
    	exit (1);
    }

    int nvar = P->nvar;

    mpz_t one;
    mpz_init(one);
    mpz_set_si (one, 1l);

    degree_t degP = mainLeadingDegree_AAZ(P);// (P->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t degQ = mainLeadingDegree_AAZ(Q);//(Q->elems[0].degs & mvarMask) >> mvarDegOffset;

    if (nvar == 1){
		return univariateGCD_AAZ (P, Q);
    }

    if (degP == 0 || degQ == 0){
    	return makeConstPolynomial_AAZ (1, nvar, one);
    }


    AltArrsZ_t* SC;
    int size = 0, idx = 0;
    DucosSubresultantChainZ (P, Q, &SC, &size);
    AltArrsZ_t* cur = SC;
    while (idx < size) {

	if (cur->poly != NULL && cur->poly->size != 0){
	    break;
	}

	cur = cur->next;
	++idx;
    }

    AltArrZ_t* res = deepCopyPolynomial_AAZ (cur->poly);
    freeAltArrsZ (SC);

    return res;
}


/*****************************
 * Subresultant Chain Wrapper
******************************/

SubresultantChainMode ChooseTheBestSubresultantChainModeZ (AltArrZ_t* P, AltArrZ_t* Q)
{
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Try To Find The Best for [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif

	if (P->nvar != Q->nvar){
		fprintf (stderr, "SMZP Error, expected inputs with the same number of variables but get %d != %d\n", P->nvar, Q->nvar);
		exit(1);
    }
	degree_t degP = mainLeadingDegree_AAZ(P);
    degree_t degQ = mainLeadingDegree_AAZ(Q);
	degree_t min_deg = MIN(degP, degQ);
	if (P->nvar>2 ||
		P == NULL || P->size == 0 || degP == 0 ||
		Q == NULL || Q->size == 0 || degQ == 0 ) {
		return DUCOS;
	}

	mpz_t inf_norm_P, inf_norm_Q;

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	mpz_inits (inf_norm_P, inf_norm_Q, NULL);
	infinityNorm_AAZ (P, inf_norm_P);
	infinityNorm_AAZ (Q, inf_norm_Q);
	fprintf(stderr, "size(infinityNorm(P) = %ld\n", (long) mpz_size(inf_norm_P));
	fprintf(stderr, "size(infinityNorm(Q) = %ld\n", (long) mpz_size(inf_norm_Q));
	mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);	
#endif

	if (P->nvar == 1) {
	 	if ((degP < 10)  && (degQ < 10)) {
			// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
			return DUCOS;
		} else if (degP == degQ) {
			return DUCOS;
		}

		mpz_inits (inf_norm_P, inf_norm_Q, NULL);
		infinityNorm_AAZ (P, inf_norm_P);
		infinityNorm_AAZ (Q, inf_norm_Q);
		int max_mpz_size = MAX(mpz_size(inf_norm_P), mpz_size(inf_norm_Q));
		if ( min_deg < 12  &&  max_mpz_size > 10 ) {
			mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);	
			return DUCOS;
		}
		mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);	
		return ModUnivar;
	} else if (P->nvar == 2) {
		// if((degP < 4) && (degQ < 4)) {
		// 	// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
		// 	return DUCOS;
		// } else 
		if (min_deg < 17) {
			degree_t deg_p_y = partialDegree_AAZ (P, 1);
			degree_t deg_q_y = partialDegree_AAZ (Q, 1);
			// fprintf(stderr, "deg(p, y) = %d, deg(q, y)=%d\n", deg_p_y, deg_q_y);
			if (deg_p_y < deg_q_y) {
				// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
				return DUCOS;
			} else if ((degP >= deg_p_y) && (degP - deg_p_y < 15) && (degQ >= deg_q_y) && (degQ - deg_q_y < 15)) {
				// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
				// fprintf(stderr, "deg(p, x) = %ld, deg(p, y) = %ld, deg(q, x) = %ld, deg(q, y) = %ld\n", degP, deg_p_y, degQ, deg_q_y);
				return DUCOS;
			} else if ((deg_p_y - degP < 15 ) && (deg_q_y - degQ < 15 )) {
				// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
				return DUCOS;
			}
		}
		// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
		return ModBivarFFT;
	}
	// mpz_clear(inf_norm_P);  mpz_clear(inf_norm_Q);
	return DUCOS;
}

SubresultantChainMode SubresultantChainWrapperZ_rev (SubresultantChainMode mode, AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len)
{

#if defined(DUCOS_SUBRES_MODE) && (DUCOS_SUBRES_MODE == 1) 
	#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
		fprintf(stderr, "*** DUCOS_SUBRES_MODE=ON [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
	#endif
	mode = DUCOS;
#endif

	SubresultantChainMode new_mode, flaged_mode = DUCOS;
	if (mode == AUTO) {
		new_mode = ChooseTheBestSubresultantChainModeZ(P, Q);
	} else {
		new_mode = mode;
	}

	AltArrsZ_t *SS = NULL;
	int flag, size = 0, no_prime = 0;

	if (new_mode == DUCOS) {
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** DUCOS SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif

		DucosSubresultantChainZ_rev(P, Q, &SS, &size, 1); // type = 1 := new cache friendly Ducos Opt. algorithm

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Done! w/ DUCOS SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif
		*len = size;
		*SC = SS;
		return new_mode;
	} else if (new_mode == OldDUCOS) {
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** OldDUCOS SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif

		DucosSubresultantChainZ_rev(P, Q, &SS, &size, 0); // type = 1 := new cache friendly Ducos Opt. algorithm

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Done! w/ OldDUCOS SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif
		*len = size;
		*SC = SS;
		return new_mode;
	}

	degree_t degP = mainLeadingDegree_AAZ(P);
	degree_t degQ = mainLeadingDegree_AAZ(Q);
	AltArrZ_t *p=NULL, *q=NULL;
	if (degP >= degQ){ p = P; q = Q; }
	else { p = Q; q = P; }
	degP = mainLeadingDegree_AAZ(p);
	degQ = mainLeadingDegree_AAZ(q);

	if (new_mode == ModUnivar) {
		AltArrsZ_t *newPoly=NULL, *tail=NULL;
		int chain_size = 0;

		DUZP_t* uP = convertFromAltArrZ_DUZP (p);
		DUZP_t* uQ = convertFromAltArrZ_DUZP (q);

		mpz_t uBound;
		mpz_init (uBound);
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** ModUnivar SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif
		DUZP_t** usubres = modularSubresultantChain_DUZP (uP, uQ, &chain_size, uBound, 0) ; // deterministic = 0

		freePolynomial_DUZP (uP);
		freePolynomial_DUZP (uQ);
		mpz_clear (uBound);

		if (chain_size) {
			newPoly = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			newPoly->poly = deepCopyPolynomial_AAZ (p);
			newPoly->next = NULL;
			SS = newPoly;
			tail = newPoly;
			size++;
			newPoly = NULL;
			newPoly = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			newPoly->poly = deepCopyPolynomial_AAZ (q);
			newPoly->next = NULL;
			tail->next = newPoly;
			tail = newPoly;
			size++;
			for (int i = chain_size-1; i >= 0; --i) {
				newPoly = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
				newPoly->poly = convertToAltArrZ_DUZP (usubres[i]);
				newPoly->next = NULL;
				tail->next = newPoly;
				tail = newPoly;
				size++;
				freePolynomial_DUZP (usubres[i]);
			}
			free (usubres);
			*len = size;
    		*SC = SS;
		}
		return new_mode;
	} else if (new_mode == ModBivar) {
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** ModBivar SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif
		flag = biModularSubresultantChainZ (p, q, &SS, &size, &no_prime);
		// fprintf(stderr, "mode = ModBivar : No. CRT-iter = %d\n", no_prime);
		if (!flag) {
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** ModBivar-FLAG FAILED! [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif
			// fprintf(stderr, "mod = ModBivar : Thresholds Reached!\n");
			new_mode = flaged_mode;
			SubresultantChainWrapperZ_rev(new_mode, p, q, &SS, &size);
		}
		*len = size;
		*SC = SS;
		return new_mode;
	} else if (new_mode == ModBivarFFT) {
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** ModBivarFFT SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
	fprintf(stderr, "*** [deg(P, y) = %d, deg(Q, y) = %d]\n",  partialDegree_AAZ (P, 1),  partialDegree_AAZ (Q, 1));
#endif
		// flag = biModularFFTSubresultantChainZ (p, q, &SS, &size, &no_prime);
		flag = biModularFFT4SubresultantChainZ (p, q, &SS, &size, &no_prime);
		// fprintf(stderr, "mode = ModBivarFFT : No. CRT-iter = %d\n", no_prime);
		if (!flag) { 
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** ModBivarFFT-FLAG FAILED! [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", P->nvar, mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
#endif
			// fprintf(stderr, "mod = ModBivarFFT : Thresholds Reached!\n");
			new_mode = flaged_mode;
			SubresultantChainWrapperZ_rev(new_mode, p, q, &SS, &size);
		}
		*len = size;
		*SC = SS;
		return new_mode;
	}

	*len = 0;
	*SC = NULL;
	return SIGNAL;
}

void _reverseSRCZ (SubresultantChainMode mode, AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len)
{
    AltArrsZ_t* invSC;
    int size = 0;
	SubresultantChainWrapperZ_rev (mode, P, Q, &invSC, &size);

    // reverse the subresultantchain:
    if (size > 1){
	AltArrsZ_t* curr = invSC;
	AltArrsZ_t* prev = NULL;
	AltArrsZ_t* next = NULL;
	while (curr != NULL){
	    next = curr->next;
	    curr->next = prev;
	    prev = curr;
	    curr = next;
	}
	*SC = prev;
	*len = size;
	return;
    }
    *SC = invSC;
    *len = size;
    return;
}

void SubresultantChainZ (SubresultantChainMode mode, AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len)
{
    AltArrsZ_t* invSC;
    int size = 0;
	SubresultantChainWrapperZ_rev (mode, P, Q, &invSC, &size);

    // reverse the subresultantchain:
    if (size > 1){
	AltArrsZ_t* curr = invSC;
	AltArrsZ_t* prev = NULL;
	AltArrsZ_t* next = NULL;
	while (curr != NULL){
	    next = curr->next;
	    curr->next = prev;
	    prev = curr;
	    curr = next;
	}
	*SC = prev;
	*len = size;
	return;
    }
    *SC = invSC;
    *len = size;
    return;
}

void SubresultantChainAtIdxZ (SubresultantChainMode mode, AltArrZ_t* P, AltArrZ_t* Q, int idx, 
							AltArrZ_t **sr_idx, AltArrZ_t **sr_idx1, specSRC_AAZ **srcType) {
#if defined(LAZY_SRC) && !LAZY_SRC
	fprintf (stderr, "SMZP Warning: calling non-lazy SRC for idx=%d\n", idx);
	DucosSubresultantChainAtIdxZ(P, Q, idx, sr_idx, sr_idx1);	
	return; 
#endif

    // Debug Print
	// fprintf(stderr, "[SMZP] SubresultantChainAtIdxZ mdeg(p)= %d mdeg(q)= %d index= %d\n", mainDegree_AAZ(P), mainDegree_AAZ(Q), idx);

	specSRC_AAZ *srcType_local = NULL;
	if (srcType == NULL || (*srcType) == NULL ) {
		srcType_local = (specSRC_AAZ*) calloc (1, sizeof(specSRC_AAZ));
	} else {
		srcType_local = *srcType;
	} srcType_local->types[0] = 0; srcType_local->types[1] = idx;
	if (idx < 0) { // srcType->sr_idx = NULL; srcType->sr_idx1 = NULL;
		if (srcType != NULL && (*srcType) == NULL) {
			*srcType = srcType_local; 
		}
		return;
	} else if (idx >= MIN_spX (mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q))) {
		*sr_idx = deepCopyPolynomial_AAZ (Q);
		*sr_idx1 = deepCopyPolynomial_AAZ (P);
		if (srcType != NULL && (*srcType) == NULL) {
			*srcType = srcType_local; 
		}
		return;	
	} else if (mode == AUTO) {
		mode = ChooseTheBestSubresultantChainModeZ (P, Q);
	}
	int size=0, pcIdx, info_size;
	int curDeg, prevDeg;
	AltArrsZ_t  *SS, *cur, *curP;
	specAQRArray_spX_t *uspecInfoArray = NULL;
	specAQRArray_spXY_t *bspecInfo = NULL;
	if (mode == DUCOS || mode == OldDUCOS) {
		if (!srcType_local->types[2]) {
			_reverseSRCZ (mode, P, Q, &SS, &size); // AltArrsZ_t * cur = SS;
			srcType_local->types[2] = size; srcType_local->gspecInfo = SS;
		}
		if (srcType_local->types[2]) {
			cur = srcType_local->gspecInfo;
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			if (idx == 0) {
				*sr_idx = deepCopyPolynomial_AAZ (cur->poly); 
				cur = cur->next;
				while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
				if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
				else { *sr_idx1 = NULL; }
			} else {
				curP = cur; prevDeg = curDeg;
				while (cur != NULL && curDeg < idx) {
					if (prevDeg < curDeg) { curP = cur; prevDeg = curDeg; }
					cur = cur->next;
					curDeg = mainLeadingDegree_AAZ (cur->poly);
				}
				if (curDeg == idx) {
					*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
					cur = cur->next; 
					while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
					if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
					else { *sr_idx1 = NULL; }
				} else {
					*sr_idx = deepCopyPolynomial_AAZ (curP->poly);
					if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
					else { *sr_idx1 = NULL; }
				}
			}
		}
	} else if (mode == ModUnivar) {
		if (!srcType_local->types[3] 
			&& !srcType_local->types[2]) {
			_reverseSRCZ (mode, P, Q, &SS, &size); 
			srcType_local->types[2] = size; srcType_local->gspecInfo = SS;
		}
		if (srcType_local->types[2]) {
			cur = srcType_local->gspecInfo;
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			if (idx == 0) {
				*sr_idx = deepCopyPolynomial_AAZ (cur->poly); 
				cur = cur->next;
				while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
				if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
				else { *sr_idx1 = NULL; }
			} else {
				curP = cur; prevDeg = curDeg;
				while (cur != NULL && curDeg < idx) {
					if (prevDeg < curDeg) { curP = cur; prevDeg = curDeg; }
					cur = cur->next;
					curDeg = mainLeadingDegree_AAZ (cur->poly);
				}
				if (curDeg == idx) {
					*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
					cur = cur->next; 
					while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
					if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
					else { *sr_idx1 = NULL; }
				} else {
					*sr_idx = deepCopyPolynomial_AAZ (curP->poly);
					if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
					else { *sr_idx1 = NULL; }
				}
			}
		} else {
			AltArrsZ_t *newPoly=NULL, *tail=NULL;
			int chain_size;
			DUZP_t* uP = convertFromAltArrZ_DUZP (P);
			DUZP_t* uQ = convertFromAltArrZ_DUZP (Q);
			mpz_t uBound;
			mpz_init (uBound);
			uspecInfoArray = (specAQRArray_spX_t *) malloc (sizeof(specAQRArray_spX_t));
			uspecInfoArray->uspecArray = srcType_local->uspecInfo;	uspecInfoArray->size = srcType_local->types[3];
			DUZP_t** usubres = regularGCDUnivariateSpeculativeSRC_DUZP (uP, uQ, idx, 
								&chain_size, NULL, uBound, &uspecInfoArray, &info_size, 0); 
			free(uspecInfoArray);
			freePolynomial_DUZP (uP);
			freePolynomial_DUZP (uQ);
			mpz_clear (uBound);
			if (chain_size == 2) {
				*sr_idx = convertToAltArrZ_DUZP(usubres[0]);
				*sr_idx1 = convertToAltArrZ_DUZP(usubres[1]);
				freePolynomial_DUZP(usubres[0]);
				freePolynomial_DUZP(usubres[1]);
			} else if (chain_size == 1){
				*sr_idx = NULL;
				*sr_idx1 = convertToAltArrZ_DUZP(usubres[0]);
				freePolynomial_DUZP(usubres[0]);
			} 
			if (usubres != NULL) { free (usubres); }
		}
	} else if (mode == ModBivar || mode == ModBivarFFT) {
		// if (srcType_local->types[2]) {
		// 	cur = srcType_local->gspecInfo;
		// 	while (cur != NULL) {
		// 		// fprintf(stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 		if (mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 			if (cur->next != NULL && mainLeadingDegree_AAZ(cur->next->poly) > idx) {
		// 				// we find the sr_idx and idx1
		// 				*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
		// 				*sr_idx1 = deepCopyPolynomial_AAZ (cur->next->poly);
		// 			} else {
		// 				curP = cur;
		// 				while (cur != NULL && mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 				// fprintf(stderr, "[inner] mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 					if (!isZero_AAZ (cur->poly)) {
		// 						curP = cur;
		// 					}
		// 					cur = cur->next;
		// 				}
		// 				if (curP != NULL) {
		// 					*sr_idx = deepCopyPolynomial_AAZ (curP->poly);
		// 				} 
		// 				if (cur != NULL) {
		// 					*sr_idx1 = deepCopyPolynomial_AAZ (cur->poly);
		// 				}
		// 			}
		// 			break;
		// 		}
		// 		// fprintf(stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 		cur = cur->next;
		// 	}
		if (!srcType_local->types[4]
			&& !srcType_local->types[2]) {
			if (MIN_spX(partialDegree_AAZ(P, 1), partialDegree_AAZ(Q, 1)) > BILAZY_CROSSOVER 
					&& srcType_local->types[0]) {
				regularGCDBiModularFFTSubresultantChainZ (P, Q, idx, &SS, &size, NULL, &bspecInfo, &info_size, 0);
				srcType_local->bspecInfo = bspecInfo->bspecArray; srcType_local->types[4] = bspecInfo->size;
				free(bspecInfo);
			} else {
				_reverseSRCZ (mode, P, Q, &SS, &size);
				srcType_local->types[2] = size; srcType_local->gspecInfo = SS;
			}
		}
		// 	cur = srcType_local->gspecInfo; curP = cur;
		// 	while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) > idx) { curP = cur; cur = cur->next; }
		// 	if (cur != NULL) {
		// 		if (mainLeadingDegree_AAZ(cur->poly) == idx) {
		// 			if (cur->next == NULL || mainLeadingDegree_AAZ (cur->next->poly) < idx) {
		// 				if (sr_idx != NULL)
		// 					*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
		// 			} else {
		// 				while (cur->next != NULL && mainLeadingDegree_AAZ (cur->next->poly) == idx) {
		// 					cur = cur->next;
		// 				}
		// 				if (sr_idx != NULL) 
		// 					*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
		// 			}
		// 		} else {
		// 			if (sr_idx != NULL)
		// 				*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
		// 		}
		// 		if (sr_idx1 != NULL)  
		// 			*sr_idx1 =  deepCopyPolynomial_AAZ (curP->poly); 
		// 	} else {
		// 		if (sr_idx != NULL)  
		// 			*sr_idx = NULL;
		// 		if (sr_idx1 != NULL)  
		// 			*sr_idx1 =  deepCopyPolynomial_AAZ (curP->poly); 
		// 	}
		// }

		// }
		if (srcType_local->types[2]) {
			cur = srcType_local->gspecInfo;
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			if (idx == 0) {
				*sr_idx = deepCopyPolynomial_AAZ (cur->poly); 
				cur = cur->next;
				while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
				if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
				else { *sr_idx1 = NULL; }
			} else {
				curP = cur; prevDeg = curDeg;
				while (cur != NULL && curDeg < idx) {
					if (prevDeg < curDeg) { curP = cur; prevDeg = curDeg; }
					cur = cur->next;
					curDeg = mainLeadingDegree_AAZ (cur->poly);
				}
				if (curDeg == idx) {
					*sr_idx = deepCopyPolynomial_AAZ (cur->poly);
					cur = cur->next; 
					while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
					if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
					else { *sr_idx1 = NULL; }
				} else {
					*sr_idx = deepCopyPolynomial_AAZ (curP->poly);
					if (cur != NULL) { *sr_idx1 = deepCopyPolynomial_AAZ (cur->poly); }
					else { *sr_idx1 = NULL; }
				}
			}
		} else {
			bspecInfo = (specAQRArray_spXY_t *) malloc (sizeof(specAQRArray_spXY_t));
			bspecInfo->bspecArray = srcType_local->bspecInfo ; bspecInfo->size = srcType_local->types[4];
			regularGCDBiModularFFTSubresultantChainZ (P, Q, idx, &SS, &size, NULL, &bspecInfo, &info_size, 0);
			free(bspecInfo);
			if (size == 2) {
				*sr_idx = deepCopyPolynomial_AAZ (SS->poly);
				*sr_idx1 = deepCopyPolynomial_AAZ (SS->next->poly);
			} else if (size == 1) {
				*sr_idx = deepCopyPolynomial_AAZ (SS->poly);
			}
			freeAltArrsZ(SS);
		}
	}
	if (srcType != NULL) {
		*srcType = srcType_local; 
	} else {
		freeSpecSRC_AAZ (srcType_local);
	}
}

void SubresultantInitAtIdxZ (SubresultantChainMode mode, AltArrZ_t* P, AltArrZ_t* Q, int idx, 
								PolyPairZ_t **s_idx, PolyPairZ_t **s_idx1, specSRC_AAZ **srcType)
{
#if defined(LAZY_SRC) && !LAZY_SRC
	fprintf (stderr, "SMZP Warning: calling non-lazy SRC-initial for idx=%d\n", idx);
	AltArrZ_t *sc_idx = NULL, *sc_idx1 = NULL;
	AltArrsZ_t * pcoefs = NULL; int pc_sz;
	DucosSubresultantChainAtIdxZ_withPrincipleCoefs (P, Q, idx, &sc_idx, &sc_idx1, &pcoefs, &pc_sz);
	*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (sc_idx), mainLeadingDegree_AAZ (sc_idx));
	*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (sc_idx1), mainLeadingDegree_AAZ (sc_idx1));
	freePolynomial_AAZ (sc_idx);
	freePolynomial_AAZ (sc_idx1);
	freeAltArrsZ(pcoefs);
	return;
#endif
	// fprintf(stderr, "[SMZP] SubresultantChainAtIdxZ mdeg(p)= %d mdeg(q)= %d index= %d\n", mainDegree_AAZ(P), mainDegree_AAZ(Q), idx);
	specSRC_AAZ *srcType_local = NULL;
	if (srcType == NULL || (*srcType) == NULL ) {
		srcType_local = (specSRC_AAZ*) calloc (1, sizeof(specSRC_AAZ));
	} else {
		srcType_local = *srcType;
	}	srcType_local->types[0] = 1; srcType_local->types[1] = idx;
	if (idx < 0) {
		if (srcType != NULL && (*srcType) == NULL) {
			*srcType = srcType_local; 
		}
		return;
	} else if (idx >= MIN_spX (mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q))) {
		*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (Q), mainLeadingDegree_AAZ (Q));
		*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (P), mainLeadingDegree_AAZ (P)); 
		if (srcType != NULL && (*srcType) == NULL) {
			*srcType = srcType_local; 
		}
		return;	
	} else if (mode == AUTO) {
		mode = ChooseTheBestSubresultantChainModeZ (P, Q);
	}

	int size, info_size;
	AltArrsZ_t  *SS, *cur, *curP;
	degree_t curDeg, prevDeg;
	specAQRArray_spX_t *uspecInfoArray = NULL;
	specAQRArray_spXY_t *bspecInfo = NULL; 
	if (mode == DUCOS || mode == OldDUCOS) {
		if (!srcType_local->types[2]) {
			_reverseSRCZ (mode, P, Q, &SS, &size);
			srcType_local->types[2] = size; srcType_local->gspecInfo = SS;
		}
		if (srcType_local->types[2]) {
			cur = srcType_local->gspecInfo;
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			if (idx == 0) {
				*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
				cur = cur->next;
				while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
				if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
				else { *s_idx1 = NULL; }
			} else {
				curP = cur; prevDeg = curDeg;
				while (cur != NULL && curDeg < idx) {
					if (prevDeg < curDeg) { curP = cur; prevDeg = curDeg; }
					cur = cur->next;
					curDeg = mainLeadingDegree_AAZ (cur->poly);
				}
				if (curDeg == idx) {
					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
					cur = cur->next; 
					while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
					if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
					else { *s_idx1 = NULL; }
				} else {
					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly));
					if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
					else { *s_idx1 = NULL; }
				}
			}
		}
		// if (srcType_local->types[2]) {
		// 	cur = srcType_local->gspecInfo;
		// 	while (cur != NULL) {
		// 		if (mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 			if (cur->next != NULL && mainLeadingDegree_AAZ(cur->next->poly) > idx) {
		// 				// we find the sr_idx and idx1
		// 				*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
		// 				*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->next->poly), mainLeadingDegree_AAZ (cur->next->poly));
		// 			} else {
		// 				curP = cur;
		// 				while (cur != NULL && mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 					if (!isZero_AAZ (cur->poly)) {
		// 						curP = cur;
		// 					}
		// 					cur = cur->next;
		// 				}
		// 				if (curP != NULL) {
		// 					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly)); 
		// 				} 
		// 				if (cur != NULL) {
		// 					*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
		// 				}
		// 			}
		// 			break;
		// 		}
		// 		cur = cur->next;
		// 	}
		// }
	} else if (mode == ModUnivar) {
		if (!srcType_local->types[2]) {
			_reverseSRCZ (mode, P, Q, &SS, &size);
			srcType_local->types[2] = size; srcType_local->gspecInfo = SS;
		}
		// if (srcType_local->types[2]) {
		// // 	cur = srcType_local->gspecInfo; curP = cur;
		// // 	while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) > idx) {
		// // 		curP = cur; cur = cur->next;				
		// // 	} 
		// // 	if (cur != NULL) {
		// // 		if (mainLeadingDegree_AAZ(cur->poly) == idx) { 
		// // 			if (cur->next == NULL || mainLeadingDegree_AAZ (cur->next->poly) < idx) {
		// // 				if (s_idx != NULL) { *s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
		// // 			} else {
		// // 				while (cur->next != NULL && mainLeadingDegree_AAZ (cur->next->poly) == idx) { cur = cur->next; } // TODO: this is not the case for ModUnivar (remove after testing)
		// // 				if (s_idx != NULL) { *s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
		// // 			}
		// // 		} else { if (s_idx != NULL) { *s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
		// // 		}
		// // 		if (s_idx1 != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly)); }
		// // 	} else {
		// // 		if (s_idx != NULL) { *s_idx = NULL;	}
		// // 		if (s_idx1 != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly)); }
		// // 	}
		// // } else {
		// // 	if (s_idx != NULL) { *s_idx = NULL; }
		// // 	if (s_idx1 != NULL) { *s_idx1 = NULL; }
		// // }
		// 	cur = srcType_local->gspecInfo;
		// 	while (cur != NULL) {
		// 		// fprintf(stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 		if (mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 			if (cur->next != NULL && mainLeadingDegree_AAZ(cur->next->poly) > idx) {
		// 				// we find the sr_idx and idx1
		// 				*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
		// 				*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->next->poly), mainLeadingDegree_AAZ (cur->next->poly));
		// 			} else {
		// 				curP = cur;
		// 				while (cur != NULL && mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 				// fprintf(stderr, "[inner] mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 					if (!isZero_AAZ (cur->poly)) {
		// 						curP = cur;
		// 					}
		// 					cur = cur->next;
		// 				}
		// 				if (curP != NULL) {
		// 					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly)); 
		// 				} 
		// 				if (cur != NULL) {
		// 					*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
		// 				}
		// 			}
		// 			break;
		// 		}
		// 		// fprintf(stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 		cur = cur->next;
		// 	}
		// }
		if (srcType_local->types[2]) {
			cur = srcType_local->gspecInfo;
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			if (idx == 0) {
				*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
				cur = cur->next;
				while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
				if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
				else { *s_idx1 = NULL; }
			} else {
				curP = cur; prevDeg = curDeg;
				while (cur != NULL && curDeg < idx) {
					if (prevDeg < curDeg) { curP = cur; prevDeg = curDeg; }
					cur = cur->next;
					curDeg = mainLeadingDegree_AAZ (cur->poly);
				}
				if (curDeg == idx) {
					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
					cur = cur->next; 
					while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
					if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
					else { *s_idx1 = NULL; }
				} else {
					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly));
					if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
					else { *s_idx1 = NULL; }
				}
			}
		} else {
			AltArrsZ_t *newPoly=NULL, *tail=NULL;
			int chain_size;
			DUZP_t* uP = convertFromAltArrZ_DUZP (P);
			DUZP_t* uQ = convertFromAltArrZ_DUZP (Q);
			mpz_t uBound;
			mpz_init (uBound);
			uspecInfoArray = (specAQRArray_spX_t *) malloc (sizeof(specAQRArray_spX_t));
			uspecInfoArray->uspecArray = srcType_local->uspecInfo;	uspecInfoArray->size = srcType_local->types[3];
			polysize_t degs[2];
			DUZP_t** usubres = regularGCDUnivariateSpeculativeSRC_DUZP (uP, uQ, idx, 
								&chain_size, degs, uBound, &uspecInfoArray, &info_size, 1); 
			free(uspecInfoArray);
			freePolynomial_DUZP (uP);
			freePolynomial_DUZP (uQ);
			mpz_clear (uBound);
			if (chain_size == 2) {
				*s_idx = createPolyPair_AAZ_inp (convertToAltArrZ_DUZP(usubres[0]), degs[0]);
				*s_idx1 = createPolyPair_AAZ_inp (convertToAltArrZ_DUZP(usubres[1]), degs[1]);
				freePolynomial_DUZP(usubres[0]);
				freePolynomial_DUZP(usubres[1]);
			} else if (chain_size == 1){
				*s_idx = NULL;
				*s_idx1 = createPolyPair_AAZ_inp (convertToAltArrZ_DUZP(usubres[0]), degs[0]);
				freePolynomial_DUZP(usubres[0]);
			} 
			if (usubres != NULL) { free (usubres); }
		}
	} else if (mode == ModBivar || mode == ModBivarFFT) {
		if (!srcType_local->types[4] 
			&& !srcType_local->types[2]) {
			if (MIN_spX(partialDegree_AAZ(P, 1), partialDegree_AAZ(Q, 1)) > BILAZY_CROSSOVER && !srcType_local->types[0]) {
				regularGCDBiModularFFTSubresultantChainZ (P, Q, idx, &SS, &size, NULL, &bspecInfo, &info_size, 1);
				srcType_local->bspecInfo = bspecInfo->bspecArray; srcType_local->types[4] = bspecInfo->size;
				free(bspecInfo);
			} else {
				_reverseSRCZ (mode, P, Q, &SS, &size);
				srcType_local->types[2] = size; srcType_local->gspecInfo = SS;
			}
		}
		// if (srcType_local->types[2]) {
		// 	cur = srcType_local->gspecInfo;
		// 	while (cur != NULL) {
		// 		// fprintf(stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 		if (mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 			if (cur->next != NULL && mainLeadingDegree_AAZ(cur->next->poly) > idx) {
		// 				// we find the sr_idx and idx1
		// 				*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
		// 				*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->next->poly), mainLeadingDegree_AAZ (cur->next->poly));
		// 			} else {
		// 				curP = cur;
		// 				while (cur != NULL && mainLeadingDegree_AAZ(cur->poly) <= idx) {
		// 				// fprintf(stderr, "[inner] mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 					if (!isZero_AAZ (cur->poly)) {
		// 						curP = cur;
		// 					}
		// 					cur = cur->next;
		// 				}
		// 				if (curP != NULL) {
		// 					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly)); 
		// 				} 
		// 				if (cur != NULL) {
		// 					*s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
		// 				}
		// 			}
		// 			break;
		// 		}
		// 		// fprintf(stderr, "mainLeadingDegree_AAZ (cur->poly) = %d\n", mainLeadingDegree_AAZ (cur->poly));
		// 		cur = cur->next;
		// 	}
		// }
		if (srcType_local->types[2]) {
			cur = srcType_local->gspecInfo;
			curDeg = mainLeadingDegree_AAZ (cur->poly);
			if (idx == 0) {
				*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
				cur = cur->next;
				while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
				if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
				else { *s_idx1 = NULL; }
			} else {
				curP = cur; prevDeg = curDeg;
				while (cur != NULL && curDeg < idx) {
					if (prevDeg < curDeg) { curP = cur; prevDeg = curDeg; }
					cur = cur->next;
					curDeg = mainLeadingDegree_AAZ (cur->poly);
				}
				if (curDeg == idx) {
					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly));
					cur = cur->next; 
					while (cur != NULL && mainLeadingDegree_AAZ (cur->poly) <= curDeg) { cur = cur->next; }
					if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
					else { *s_idx1 = NULL; }
				} else {
					*s_idx = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (curP->poly), mainLeadingDegree_AAZ (curP->poly));
					if (cur != NULL) { *s_idx1 = createPolyPair_AAZ_inp (mainLeadingCoefficient_AAZ (cur->poly), mainLeadingDegree_AAZ (cur->poly)); }
					else { *s_idx1 = NULL; }
				}
				// fprintf(stderr, "nvar= %d, idx= %d, mdegs = (%d, %d)\n", P->nvar, idx, ((*s_idx) != NULL ? (*s_idx)->d : -1), ((*s_idx1) != NULL ? (*s_idx1)->d : -1));
			}
		} else {
			bspecInfo = (specAQRArray_spXY_t *) malloc (sizeof(specAQRArray_spXY_t));
			bspecInfo->bspecArray = srcType_local->bspecInfo; bspecInfo->size = srcType_local->types[4];
			polysize_t bdegs[2];
			regularGCDBiModularFFTSubresultantChainZ (P, Q, idx, &SS, &size, bdegs, &bspecInfo, &info_size, 1);
			free(bspecInfo);
			if (size == 2) {
				*s_idx = createPolyPair_AAZ_inp (deepCopyPolynomial_AAZ(SS->poly), bdegs[0]);
				*s_idx1 =  createPolyPair_AAZ_inp (deepCopyPolynomial_AAZ(SS->next->poly), bdegs[1]);
			} else if (size == 1) {
				*s_idx = createPolyPair_AAZ_inp (deepCopyPolynomial_AAZ(SS->poly), bdegs[0]);
			}   
				// fprintf(stderr, "nvar= %d, idx= %d, mdegs = (%d, %d)\n", P->nvar, idx, ((*s_idx) != NULL ? (*s_idx)->d : -1), ((*s_idx1) != NULL ? (*s_idx1)->d : -1));
			freeAltArrsZ(SS);
		}
	}
	if (srcType != NULL) {
		*srcType = srcType_local; 
	} else {
		freeSpecSRC_AAZ (srcType_local);
	}
}


///////////////////////////////
// Extended Subresultant Chain
///////////////////////////////

AltArrZ_t* semiLazardOpt_AAZ  (AltArrZ_t* Sd, AltArrZ_t* Sdm, AltArrZ_t* s)
{
    if (Sd == NULL || Sd->size == 0){
		fprintf (stderr, "SMZP Error: In Semi-LazardOpt , Sd is NULL!\n");
		exit (EXIT_FAILURE);
    }

    if (Sdm == NULL || Sdm->size == 0){
		return NULL;
    }

    if (Sd->nvar != Sdm->nvar){
	fprintf (stderr, "SMZP Error: In Semi-LazardOpt , Sd->nvar(=%d) != Sdm->nvar(=%d).\n",
		 Sd->nvar, Sdm->nvar);
	exit (EXIT_FAILURE);
    }

    register int nvar = Sd->nvar;

    degree_t n = mainLeadingDegree_AAZ(Sd) - mainLeadingDegree_AAZ(Sdm) - 1;

    if (n < 1){
    	mpz_t one;
    	mpz_init (one);
    	mpz_set_si (one, 1l);
		return  makeConstPolynomial_AAZ (1, Sdm->nvar, one); // deepCopyPolynomial_AAZ  (Sdm);
    }

    AltArrZ_t* x = mainLeadingCoefficient_AAZ (Sdm);
    AltArrZ_t* y = s;

    if (x == NULL || x->size == 0){
		fprintf (stderr, "SMZP Error: In Semi-LazardOpt, x is NULL!");
		exit (EXIT_FAILURE);
    }
    if (y == NULL || y->size == 0){
		fprintf (stderr, "SMZP Error: In Semi-LazardOpt, y is NULL!");
		exit (EXIT_FAILURE);
    }

    AltArrZ_t* c = deepCopyPolynomial_AAZ (x);
    register degree_t a = 1 << Log2_AAZ (n);
    n = n - a;

    AltArrZ_t* tmpz;
    AltArrZ_t* tmp;
    /* int orign = n; */

    while (a > 1) {
		tmpz = NULL;
		tmp = multiplyPolynomials_AAZ_inp (c, c, nvar);
		exactDividePolynomials_AAZ (tmp, y, &tmpz, nvar);
		freePolynomial_AAZ (tmp);
		c = tmpz;
		a = a >> 1;

		if (n >= a) {
		    tmpz = NULL;
		    if (c != NULL && c->size != 0) {
				tmp = multiplyPolynomials_AAZ_inp (c, x, nvar);
		    } else {
				tmp = NULL;
		    }
		    exactDividePolynomials_AAZ (tmp, y, &tmpz, nvar);
		    freePolynomial_AAZ (tmp);
		    c = tmpz;
		    n = n - a;
		}
    }

    /* if (!orign){ */
    /* 	return Se; */
    /* } */

    AltArrZ_t* sse = NULL; // semi-Se
	exactDividePolynomials_AAZ (c, y, &sse, nvar);
    freePolynomial_AAZ (c);
    freePolynomial_AAZ (x);

    return sse;
}

AltArrZ_t* exLazardOpt_AAZ (AltArrZ_t* Sd, exgcdsZ_t* VSdm, AltArrZ_t* s, AltArrZ_t** hc, AltArrZ_t** qc)
{

	AltArrZ_t* sse = semiLazardOpt_AAZ (Sd, VSdm->r, s); // semi-Se

	if (sse == NULL || sse->size == 0) {
		*hc = NULL;
		*qc = NULL;
		return NULL;
	}

	if (isOne_AAZ (sse)) {
		freePolynomial_AAZ (sse);

		*hc = deepCopyPolynomial_AAZ (VSdm->a);
		*qc = deepCopyPolynomial_AAZ (VSdm->b);
		return deepCopyPolynomial_AAZ (VSdm->r);
	}

	AltArrZ_t* hhc = multiplyPolynomials_AAZ (sse, VSdm->a, sse->nvar);
	AltArrZ_t* qqc = multiplyPolynomials_AAZ (sse, VSdm->b, sse->nvar);
	sse = multiplyPolynomials_AAZ_inp (sse, VSdm->r, sse->nvar); // Se

	*hc = hhc;
	*qc = qqc;
	return sse;
}


AltArrZ_t* exDucosOpt_AAZ (exgcdsZ_t* VSd, exgcdsZ_t* VSdm, AltArrZ_t* Se, AltArrZ_t* sd, AltArrZ_t** hb, AltArrZ_t** qb)
{
    int nvar = VSd->r->nvar;

    if (nvar != VSdm->r->nvar || nvar != Se->nvar || nvar != sd->nvar){
		fprintf (stderr, "SMZP ERROR: In Extended DucosOpt , inputs' nvars are different!\n");
		exit (EXIT_FAILURE);
    }

    degree_t d = mainLeadingDegree_AAZ(VSd->r);
    degree_t e = mainLeadingDegree_AAZ(VSdm->r);

    if (d == 0){
		fprintf (stderr, "SMZP Error: In Extended DucosOpt , d == 0.\n");
		exit (EXIT_FAILURE);
    }
    if (e == 0 || d < e){
		fprintf (stderr, "SMZP Error: In Extended DucosOpt , e == 0 || d > e.\n");
		exit (EXIT_FAILURE);
    }

	AltArrZ_t* cd  = mainLeadingCoefficient_AAZ (VSd->r);
    AltArrZ_t* cdm = mainLeadingCoefficient_AAZ (VSdm->r);
    AltArrZ_t* se  = mainLeadingCoefficient_AAZ (Se);

  	if (cd == NULL || cd->size ==0){
		fprintf (stderr, "SMZP Error: In Extended DucosOpt , cd = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (cdm == NULL || cdm->size ==0){
		fprintf (stderr, "SMZP Error: In Extended DucosOpt , cdm = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (se == NULL || se->size ==0){
		fprintf (stderr, "SMZP Error: In Extended DucosOpt , se = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (d == e){
		// for non-defective cases:
		fprintf (stderr, "SMZP Error: In Extended DucosOpt , d == e.\n"); // TEST
	}

	AltArrZ_t* ssd  = multiplyPolynomials_AAZ (sd, cd, nvar);
	AltArrZ_t* csdm = multiplyPolynomials_AAZ (cdm, se, nvar);
	AltArrZ_t* mul;
	if (csdm != NULL && VSd->r != NULL) {
		mul = multiplyPolynomials_AAZ (csdm, VSd->r, nvar);
	} else {
		mul = NULL;
	}

	AltArrZ_t* r = NULL;
	AltArrZ_t* q = NULL;
	dividePolynomials_AAZ (mul, VSdm->r, &q, &r, nvar);

	freePolynomial_AAZ (mul);
	freePolynomial_AAZ (se);
	freePolynomial_AAZ (cd);
	freePolynomial_AAZ (cdm);

	// Sem:
	AltArrZ_t* Sem = NULL;
	AltArrZ_t* tmp = NULL;
	if (r != NULL && r->size != 0) {
		dividePolynomials_AAZ (r, ssd, &Sem, &tmp, nvar);
		if (tmp != NULL){
			freePolynomial_AAZ (tmp); tmp = NULL;
		}
		freePolynomial_AAZ (r);
	}

	AltArrZ_t* hhb;
	AltArrZ_t* qqb;

	AltArrZ_t* tmp1 = NULL;
	AltArrZ_t* tmp2 = NULL;


	// hb:
	if (VSd->a != NULL) {
		tmp = multiplyPolynomials_AAZ (csdm, VSd->a, nvar);
	}

	if (VSdm->a != NULL) {
		tmp1 = multiplyPolynomials_AAZ (q, VSdm->a, nvar);
	}

	tmp = subPolynomials_AAZ_inp (tmp, tmp1, nvar);
	freePolynomial_AAZ (tmp1); tmp1 = NULL;

	if(tmp == NULL || tmp->size == 0) {
		hhb = NULL;
	} else {
		dividePolynomials_AAZ (tmp, ssd, &tmp1, &tmp2, nvar);
		hhb = tmp1;
		if (tmp2 != NULL) {
			freePolynomial_AAZ (tmp2);
		}
		freePolynomial_AAZ (tmp);
	}

	// qb:
	tmp = NULL;
	if (VSd->b != NULL) {
		tmp = multiplyPolynomials_AAZ (csdm, VSd->b, nvar);
	}

	tmp1 = NULL;
	if (VSdm->b != NULL) {
		tmp1 = multiplyPolynomials_AAZ (q, VSdm->b, nvar);
	}

	tmp = subPolynomials_AAZ_inp (tmp, tmp1, nvar);
	freePolynomial_AAZ (tmp1); tmp1 = NULL;

	if(tmp == NULL || tmp->size == 0) {
		qqb = NULL;
	} else {
		dividePolynomials_AAZ (tmp, ssd, &tmp1, &tmp2, nvar);
		qqb = tmp1;
		if (tmp2 != NULL){
			freePolynomial_AAZ (tmp2);
		}
		freePolynomial_AAZ (tmp);
	}

	if (q != NULL)
		freePolynomial_AAZ (q);

	freePolynomial_AAZ (csdm);
	freePolynomial_AAZ (ssd);

	if ((d-e)%2) {
		*hb = hhb;
		*qb = qqb;
		return Sem;
	}

	// fprintf(stderr, "[SMZP-Ali] neg... \n");

	if (hhb == NULL || hhb->size == 0) {
		*hb = NULL;
	} else {
		negatePolynomial_AAZ (hhb);
		*hb = hhb;
	}

	if (qqb == NULL || qqb->size == 0) {
		*qb = NULL;
	} else {
		negatePolynomial_AAZ (qqb);
		*qb = qqb;
	}

	if (Sem == NULL || Sem->size == 0){
		return NULL;
	}

	negatePolynomial_AAZ (Sem);
	return Sem;

}

void exDucosSubresultantChain_rev_AAZ (AltArrZ_t* P, AltArrZ_t* Q, exgcdsZ_t** SC, int* len, int type)
{

    if (P == NULL || P->size == 0){
		exgcdsZ_t* SS0 = (exgcdsZ_t*) malloc (sizeof(exgcdsZ_t));
		SS0->r = NULL;
		SS0->a = NULL;
		SS0->b = NULL;
		SS0->next = NULL;
		*SC = SS0;
		*len = 0;
		return;
    }

    if (Q == NULL || Q->size == 0){
	    mpz_t pone;
	    mpz_init (pone);
	    mpz_set_si (pone, 1l);

		exgcdsZ_t* SS1 = (exgcdsZ_t*) malloc (sizeof(exgcdsZ_t));
		SS1->r = deepCopyPolynomial_AAZ (P);
		SS1->a = makeConstPolynomial_AAZ (1, P->nvar, pone);
		SS1->b = NULL;
		SS1->next = NULL;
		*SC = SS1;
		*len = 1;
		mpz_clear (pone);
		return;
    }
    if (P->nvar != Q->nvar){
		fprintf (stderr, "SMZP Error: In Extended DucosSubresultantChain_rev, P->nvar(%d) != Q->nvar(%d).\n",P->nvar, Q->nvar);
		exit (EXIT_FAILURE);
	}


    degree_t degP = mainLeadingDegree_AAZ(P);
    degree_t degQ = mainLeadingDegree_AAZ(Q);


    if (degP == 0 || degQ == 0){
        fprintf (stderr, "SMZP: error, Input polynomials to exteneded Subresultant Chain must have positive degree.\n");
		exit (EXIT_FAILURE);
    }

    int lazy = 0; // lazy option in pseudoDivide algorithm
    register int nvar = P->nvar;

    int tmpE = 0;
    AltArrZ_t* tmpQ = NULL;
    AltArrZ_t* tmpH = NULL;
    AltArrZ_t* tmpB = NULL;
    AltArrZ_t* s = NULL;
    AltArrZ_t* p = NULL;
    AltArrZ_t* q = NULL;
    AltArrZ_t* A = NULL;
    AltArrZ_t* B = NULL;
    AltArrZ_t* C = NULL;

    exgcdsZ_t* SS = NULL;
    exgcdsZ_t* tail = NULL;
    exgcdsZ_t* tailm = NULL;
    degree_t delta = 0;
    int size = 0;

    if (degP >= degQ){
		p = P;
		q = Q;
    } else {
		p = Q;
		q = P;

	    degP = mainLeadingDegree_AAZ(p);
    	degQ = mainLeadingDegree_AAZ(q);
    }

    s = mainLeadingCoefficient_AAZ (q);
    s = exponentiatePoly_AAZ (s, degP - degQ, nvar);
    A = deepCopyPolynomial_AAZ (q);


    // fprintf(stderr, "degP: %d, degQ: %d\n", degP, degQ); // TEST


    // B =
    AltArrZ_t* negQ = deepCopyPolynomial_AAZ (q);
    negatePolynomial_AAZ (negQ);

    pesudoDivideAtIdx_AAZ (0, p, negQ, &tmpQ, &B, &tmpE, &tmpH, nvar, lazy);
    freePolynomial_AAZ (negQ);
    // freePolynomial_AAZ (tmpQ);
    // freePolynomial_AAZ (tmpH);

    mpz_t one;
    mpz_init (one);
    mpz_set_si (one, 1l);

    exgcdsZ_t* newPoly = (exgcdsZ_t*) malloc (sizeof(exgcdsZ_t));
    newPoly->r = deepCopyPolynomial_AAZ (p);
    newPoly->a = makeConstPolynomial_AAZ (1, P->nvar, one);
    newPoly->b = NULL;
    newPoly->next = NULL;

    SS = newPoly;
    tail = newPoly;
    size++;

    newPoly = (exgcdsZ_t*) malloc (sizeof(exgcdsZ_t));
    if (B == NULL || B->size == 0){
		newPoly->r = A;

    } else{
		newPoly->r = deepCopyPolynomial_AAZ (A);
    }
    newPoly->a = NULL;
	newPoly->b = makeConstPolynomial_AAZ (1, A->nvar, one);
    newPoly->next = NULL;

	exgcdsZ_t* VA = newPoly;
    exgcdsZ_t* VB = NULL;

    tail->next = newPoly;
    tail = newPoly;
    size++;

    register degree_t d = 0;
    register degree_t e = 0;

    while (1){

    	if (B == NULL || B->size == 0){
    		break;
    	}

		d = mainLeadingDegree_AAZ(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
		e = mainLeadingDegree_AAZ(B);//(B->elems[0].degs & mvarMask) >> mvarDegOffset;

		newPoly = (exgcdsZ_t*) malloc (sizeof(exgcdsZ_t));
		newPoly->r = deepCopyPolynomial_AAZ (B);
		newPoly->a = tmpH;
		newPoly->b = tmpQ;
		newPoly->next = NULL;

		VB = newPoly;

		tailm = tail;
		tail->next = newPoly;
		tail = newPoly;
		size++;
		// fprintf(stderr, "[SMZP-Ali] size = %d\n", size);

		delta = d - e;
		C = NULL; // added after updating, see (*)
		if (delta > 1){
			tmpH = NULL;
			tmpQ = NULL;

			C = exLazardOpt_AAZ (tailm->r, tail, s, &tmpH, &tmpQ);
			newPoly = (exgcdsZ_t*) malloc (sizeof(exgcdsZ_t));
			newPoly->r = deepCopyPolynomial_AAZ (C);
			newPoly->a = tmpH;
			newPoly->b = tmpQ;
			newPoly->next = NULL;

			tail->next = newPoly;

			tail = newPoly;
			size++;

		} else {
			C = deepCopyPolynomial_AAZ (B);
		}

		if (e == 0){ // B.degree() == 0
			break;
		}

		tmpH = NULL;
		tmpQ = NULL;

		tmpB = exDucosOpt_AAZ (VA, VB, C, s, &tmpH, &tmpQ);

		freePolynomial_AAZ (B);
		freePolynomial_AAZ (A);
		B = tmpB;
		A = C; // deepCopyPolynomial_AAZ (C); // TODO: update DucosSubresultant_rev ALgorithm, too. (*)

		VA = newPoly;

		freePolynomial_AAZ (s);
		s = mainLeadingCoefficient_AAZ (A);

		// freePolynomial_AAZ (C);
	}


    // if subresultant is 0, add it to SS:
    if (tail->r != NULL && tail->r->size != 0 && mainLeadingDegree_AAZ(tail->r) > 0) {
		newPoly = (exgcdsZ_t*) malloc (sizeof (exgcdsZ_t));
		newPoly->r = NULL;
		newPoly->a = NULL;
		newPoly->b = NULL;
		newPoly->next = NULL;

		tail->next = newPoly;
		tail = newPoly;
		size++;
    }

    *len = size;
    *SC = SS;
    return;
}

void exDucosSubresultantChain_AAZ (AltArrZ_t* P, AltArrZ_t* Q, exgcdsZ_t** SC, int* len)
{
    exgcdsZ_t* invSC;
    int size = 0;
    exDucosSubresultantChain_rev_AAZ (P, Q, &invSC, &size, 1);

    // reverse the subresultantchain:
    if (size > 1){
		exgcdsZ_t* curr = invSC;
		exgcdsZ_t* prev = NULL;
		exgcdsZ_t* next = NULL;
		while (curr != NULL){
		    next = curr->next;
		    curr->next = prev;
		    prev = curr;
		    curr = next;
		}
		*SC = prev;
		*len = size;
		return;
	}
    *SC = invSC;
    *len = size;
    return;
}

AltArrZ_t* exDucosResultant_AAZ (AltArrZ_t* P, AltArrZ_t* Q, AltArrZ_t** a, AltArrZ_t** b)
{
    exgcdsZ_t* SC;
    int size = 0;
    exDucosSubresultantChain_AAZ (P, Q, &SC, &size);

    if (size){
    	*a = deepCopyPolynomial_AAZ (SC->a);
    	*b = deepCopyPolynomial_AAZ (SC->b);
    	AltArrZ_t* r = deepCopyPolynomial_AAZ (SC->r);

    	freeExgcds_AAZ (SC);

		return r;
    }

    freeExgcds_AAZ (SC);

    *a = NULL;
    *b = NULL;
    return NULL;
}


///////////////////////////////////
///// Triangular Set Normalization
///////////////////////////////////

AltArrZ_t* normalizePolynomial_AAZ (AltArrZ_t* P, AltArrZ_t** T, AltArrZ_t** A, int s, int nvar)
{
	int v = mainVariable_AAZ (P);

	if (v == -1) {
		*A = NULL;
		return NULL;
	}

	if (isConstant_AAZ(P) > 0) {
		mpz_t aone;
		mpz_init (aone);
		mpz_set_si (aone, 1l);
		*A = makeConstPolynomial_AAZ (1, P->nvar, aone);

		return deepCopyPolynomial_AAZ (P);
	}

	int isRes = -1;

	for (int i = 0; i < s; i++) {
		if (v == mainVariable_AAZ (T[i])) {
			isRes = i;
			break;
		}
	}

	if (isRes==-1) {
		degree_t d = mainDegree_AAZ (P);
		AltArrZ_t* h = mainLeadingCoefficient_AAZ (P);
		AltArrZ_t* hv = deepCopyPolynomial_AAZ (h);
		multiplyPolynomialAtIdxByXn_AAZ_inp (h, v, d, nvar); // h = h*v^d
		AltArrZ_t* t = subPolynomials_AAZ (P, hv, nvar); // t = P - h*v^mdeg

		AltArrZ_t* q = NULL;
		AltArrZ_t* r = normalizePolynomial_AAZ (h, T, &q, s, nvar);

		*A = q;

		if (r == NULL || r->size == 0) {
			r = NULL;
		} else {
			multiplyPolynomialAtIdxByXn_AAZ_inp (r, v, d, nvar); // r = r*v^d
		}

		if (t == NULL || q == NULL ||
			t->size == 0 || q->size == 0) {
			return r;
		}

		q = multiplyPolynomials_AAZ_inp (q, t, nvar);
		AltArrZ_t* nf_r = onlyNormalForm_AAZ (q, T, s, nvar);

		r = addPolynomials_AAZ_inp (r, nf_r, nvar);

		freePolynomial_AAZ (hv);
		freePolynomial_AAZ (nf_r);
		freePolynomial_AAZ (h);
		freePolynomial_AAZ (t);

		return r;
	} else {
		AltArrZ_t* a = NULL;
		AltArrZ_t* b = NULL;
		AltArrZ_t* exr = exDucosResultant_AAZ (P, T[isRes], &a, &b);

		freePolynomial_AAZ (b);

		AltArrZ_t* nf_a   = onlyNormalForm_AAZ (a, T, s, nvar);
		AltArrZ_t* nf_exr = onlyNormalForm_AAZ (exr, T, s, nvar);

		freePolynomial_AAZ (a);
		freePolynomial_AAZ (exr);

		*A = nf_a;
		return nf_exr;
	}
}

////////
// GCD
////////
AltArrZ_t* integralGCD_AAZ_polyOut (AltArrZ_t* P, AltArrZ_t* Q)
{

    if (P->size != 1 || Q->size != 1){
	fprintf (stderr, "SMZP Error, In integralGCD_AAZ the input polynomials must have single term." );
	exit(1);
    }

    mpz_t gcd_z;
    mpz_init (gcd_z);
    mpz_gcd (gcd_z, P->elems[0].coef, Q->elems[0].coef);

    return makeConstPolynomial_AAZ (1, P->nvar, gcd_z);
}

AltArrZ_t* gcd_AAZ_Z (AltArrZ_t* P, mpz_t c)
{
    if (P == NULL) {
	return NULL;
    }

    mpz_t ret;
    mpz_init (ret);

    if (P->size == 0) {
	mpz_set_si (ret, 1l);
	return makeConstPolynomial_AAZ (1, P->nvar, ret);
    }

    integralContent_AAZ (P, ret);

    mpz_t gcd_z;
    mpz_init (gcd_z);
    mpz_gcd (gcd_z, ret, c);

    mpz_clear (ret);
    return makeConstPolynomial_AAZ (1, P->nvar, gcd_z);
}

extern float mainPrimTime;
#include "Utils/C_To_Cpp.h"
AltArrZ_t* gcd_AAZ (AltArrZ_t* P, AltArrZ_t* Q)
{

	// const char* ch[] = {"x", "y", "z", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17"}; // TEST
	return gcd_AAZ_tmp(P, Q);

	int mvarP = leadingVariable_AAZ (P);
    int mvarQ = leadingVariable_AAZ (Q);
    if (mvarP < 0){
		if (mvarP == -2){
		    return deepCopyPolynomial_AAZ (Q);
		}
		if (mvarQ == -2) {
		    return deepCopyPolynomial_AAZ (P);
		}
		return gcd_AAZ_Z (Q, P->elems[0].coef);
    }
    if (mvarQ < 0){
		return gcd_AAZ_Z (P, Q->elems[0].coef);
    }

    if (P->nvar == 1 && Q->nvar == 1){
		return univariateGCD_AAZ (P, Q);
    }

#if defined(WITH_MAPLE) && WITH_MAPLE
    return gcdByMaple_AAZ(P,Q);
#elif defined(WITH_BLAD) && WITH_BLAD
    if (P->size > 8 && Q->size > 8) {
    	char* c_names[P->nvar];
    	for (int i = 0; i < P->nvar; ++i) {
    		c_names[i] = (char*) malloc(sizeof(char)*16);
    		sprintf(c_names[i], "x_%d", i);
    	}
    	AltArrZ_t* gz = NULL;
    	gcdBLAD_AAZ(P, Q, (const char**) c_names, &gz);
    	for (int i = 0; i < P->nvar; ++i) {
    		free(c_names[i]);
    	}
    	return gz;
    }
#endif

    if (mvarP < mvarQ){
	return gcd_AAZ (Q, mainContent_AAZ (P));
    }
    if (mvarP > mvarQ){
	return gcd_AAZ (P, mainContent_AAZ (Q));
    }

    // if univariate?(p1) and univariate?(p2) then //TODO:
    // convert p1 and p2 to SUZP and call the SUZP gcd

    AltArrZ_t* contP = NULL;//(AltArrZ_t*) malloc (sizeof (AltArrZ_t));
    AltArrZ_t* contQ = NULL;//(AltArrZ_t*) malloc (sizeof (AltArrZ_t));
    AltArrZ_t* cont = NULL;
    AltArrZ_t* gcd = NULL;
    int nvar = P->nvar;
    int isShrinked = 0;
    // int ppPmvar, ppQmvar;
    if (nvar != Q->nvar){
	fprintf (stderr, "SMZP Error, In gcd_AA, the nvar of input polynomials must be the same.");
	exit(1);
    }

    AltArrZ_t* ppP = mainPrimitiveFactorization_AAZ (P, &contP);

    AltArrZ_t* ppQ = mainPrimitiveFactorization_AAZ (Q, &contQ);

    cont = gcd_AAZ (contP, contQ);

	if (mvarP > 0) {
		int varmap[nvar];
	    int j = 0;

		for (int i = 0; i < mvarP; i++) {
			varmap[i] = -1;
		}

		for (int i = mvarP; i < nvar; i++) {
			varmap[i] = j;
			j++;
		}

		shrinkAndReorderVars_AAZ (ppP, varmap, nvar);
		shrinkAndReorderVars_AAZ (ppQ, varmap, nvar);
		isShrinked = 1;
	}

    /* for (int i = 0; i < mvarP; ++i){ */
	/* shrinkNumVarsAtIdx_AAZ (ppP, 0); */
	/* shrinkNumVarsAtIdx_AAZ (ppQ, 0); */
	/* isShrinked++; */
    /* } */

	/* fprintf (stderr, "[SMZP-Ali-GCD] ppP := "); */
	/* printPoly_AAZ_unpk(stderr, ppP, ch, ppP->nvar); */
	/* fprintf (stderr, "\n"); */

	/* fprintf (stderr, "[SMZP-Ali-GCD] ppQ := "); */
	/* printPoly_AAZ_unpk(stderr, ppQ, ch, ppP->nvar); */
	/* fprintf (stderr, "\n"); */

    AltArrZ_t* lnzch = lastNonZeroChain_AAZ (ppP, ppQ);

    if (isShrinked && lnzch != NULL && lnzch->size != 0){
	expandNumVarsLeft_AAZ (lnzch, nvar);
    }

    if (leadingVariable_AAZ (lnzch) != mvarP){
	mpz_t one;
	mpz_init (one);
	mpz_set_si (one, 1l);
	gcd = makeConstPolynomial_AAZ (1, nvar, one);
    } else {
	gcd = mainPrimitivePart_AAZ (lnzch, mvarP);
    }

    gcd = multiplyPolynomials_AAZ_inp (gcd, cont, nvar);

    freePolynomial_AAZ (contP);
    freePolynomial_AAZ (contQ);
    freePolynomial_AAZ (ppP);
    freePolynomial_AAZ (ppQ);
    freePolynomial_AAZ (cont);

    return gcd;
}

AltArrZ_t* mainPrimitiveFactorization_AAZ (const AltArrZ_t* P, AltArrZ_t** cont)
{
    int mvarP = leadingVariable_AAZ (P);

    if (mvarP < 0){
    	if (mvarP == -1){
/* if (*cont == NULL || (*cont)->size == 0){ */
    		*cont =  makeConstPolynomial_AAZ (1, P->nvar, P->elems[0].coef);
/* } */
    	}
    	return makeConstIntPolynomial_AAZ (1, P->nvar, 1l);
    }

    AltArrZ_t* sP;
    int isShrinked = 0;
    mpz_t one;
    mpz_init (one);
    mpz_set_si (one, 1l);

    if (mvarP > 0){
		sP = deepCopyPolynomial_AAZ (P);
		/* for (int i = 0; i < mvarP; ++i){ */
		/*     shrinkNumVarsAtIdx_AAZ (sP, 0); */
		/* } */
		/* isShrinked = 1; */

		int varmap[P->nvar];
		int j = 0;

		for (int i = 0; i < mvarP; i++) {
			varmap[i] = -1;
		}

		for (int i = mvarP; i < P->nvar; i++) {
			varmap[i] = j;
			j++;
		}

		shrinkAndReorderVars_AAZ (sP, varmap, P->nvar);
		isShrinked = 1;
    } else {
		sP = deepCopyPolynomial_AAZ (P);
    }

    AltArrZ_t* tmp = NULL;
    AltArrZ_t* tmpi = NULL;
    RecArrZ_t* recP = convertToRecursiveArrayZ (sP);
    RecArrElemZ_t* elems = recP->elems;
    int rSize = recP->size;

    int isOne = 0, idx = 1;
    tmpi = convertFromAAZElemToAAZ
						 (elems[0].coef,  elems[0].coefSize,
						  sP->nvar, recP->unpacked);
    AltArrZ_t* coefGCD = deepCopyPolynomial_AAZ (tmpi);
    free(tmpi);

    while (idx < rSize){
    	if (isOne_AAZ(coefGCD)) {
    		isOne = 1;
    		break;
    	}
		// if (coefGCD != NULL && coefGCD->size == 1 &&
		//     leadingVariable_AAZ (coefGCD) == -1 &&
		//     !mpz_cmp(coefGCD->elems->coef, one)){
		//     isOne = 1;
		//     break;
		// }

		tmpi = convertFromAAZElemToAAZ (elems[idx].coef,
									elems[idx].coefSize,
									sP->nvar, recP->unpacked);

		tmp = gcd_AAZ(coefGCD, tmpi);

		freePolynomial_AAZ (coefGCD);
		free(tmpi); //free only the struct not underyling coefs which are shared.
		coefGCD = tmp;
		++idx;
    }

	if (isShrinked && coefGCD != NULL && coefGCD->size != 0){
	    expandNumVarsLeft_AAZ (coefGCD, P->nvar);
	}
	if (coefGCD != NULL && coefGCD->size != 0){
	    *cont = coefGCD;
	} else {
	    fprintf (stderr, "SMZP Error, in mainPrimitiveFactorization_AAZ, coefGCD is NULL!\n");
	    exit(1);
	}
    mpz_clear (one);

	convertFromRecursiveArrayZ(recP, sP->nvar);
	freePolynomial_AAZ(sP);

    if (!isOne){
		AltArrZ_t* tmpR = NULL;
		AltArrZ_t* tmpQ = NULL;

		dividePolynomials_AAZ (P, coefGCD, &tmpQ, &tmpR, P->nvar);
		freePolynomial_AAZ (tmpR);

		return tmpQ;
    }


    return deepCopyPolynomial_AAZ (P);
}

AltArrZ_t* mainPrimitivePart_AAZ (AltArrZ_t* P, int mvar)
{
    int mvarP = mvar;
    AltArrZ_t* sP;
    int isShrinked = 0;
    mpz_t one;
    mpz_init (one);
    mpz_set_si (one, 1l);

    if (mvarP < 0){
	return makeConstPolynomial_AAZ (1, P->nvar, one);
    }

    if (mvarP > 0){
		sP = deepCopyPolynomial_AAZ (P);
		/* for (int i = 0; i < mvarP; ++i){ */
		/*     shrinkNumVarsAtIdx_AAZ (sP, 0); */
		/* } */

		int varmap[P->nvar];
		int j = 0;

		for (int i = 0; i < mvarP; i++) {
			varmap[i] = -1;
		}

		for (int i = mvarP; i < P->nvar; i++) {
			varmap[i] = j;
			j++;
		}

		shrinkAndReorderVars_AAZ (sP, varmap, P->nvar);
		isShrinked = 1;
    } else {
	sP = deepCopyPolynomial_AAZ (P);
    }

    AltArrZ_t* tmp = NULL;
    AltArrZ_t* tmpi = NULL;
    RecArrZ_t* recP = convertToRecursiveArrayZ (sP);
    RecArrElemZ_t* elems = recP->elems;
    int rSize = recP->size;

    int isOne = 0, idx = 1;
    AltArrZ_t* coefGCD = deepCopyPolynomial_AAZ (convertFromAAZElemToAAZ
						 (elems[0].coef,  elems[0].coefSize,
						  sP->nvar, recP->unpacked));

    while (idx < rSize){
	if (coefGCD != NULL && coefGCD->size == 1 &&
	    leadingVariable_AAZ (coefGCD) == -1 &&
	    !mpz_cmp(coefGCD->elems->coef, one)){
	    isOne = 1;
	    break;
	}

	tmpi = deepCopyPolynomial_AAZ (convertFromAAZElemToAAZ (elems[idx].coef,
								elems[idx].coefSize,
								sP->nvar, recP->unpacked));
	tmp = gcd_AAZ (coefGCD, tmpi);
	freePolynomial_AAZ (coefGCD);
	freePolynomial_AAZ (tmpi);
	coefGCD = tmp;
	++idx;
    }

    if (isShrinked && coefGCD != NULL && coefGCD->size != 0){
	expandNumVarsLeft_AAZ (coefGCD, P->nvar);
    }
    mpz_clear (one);

    if (!isOne){
	AltArrZ_t* tmpR = NULL;
	AltArrZ_t* tmpQ = NULL;

	dividePolynomials_AAZ (P, coefGCD, &tmpQ, &tmpR, P->nvar);
	freePolynomial_AAZ (tmpR);

	return tmpQ;
    }

    return deepCopyPolynomial_AAZ (P);
}

AltArrZ_t* mainContent_AAZ (AltArrZ_t* P)
{

    int mvarP = leadingVariable_AAZ (P);
    AltArrZ_t* sP;
    int isShrinked = 0;
    mpz_t one;
    mpz_init (one);
    mpz_set_si (one, 1l);

    if (mvarP == -2){
	return makeConstPolynomial_AAZ (1, P->nvar, one);
    }
    if (mvarP == -1){
	return deepCopyPolynomial_AAZ (P);
    }

    if (mvarP > 0){
		sP = deepCopyPolynomial_AAZ (P);
		/* for (int i = 0; i < mvarP; ++i){ */
		/*     shrinkNumVarsAtIdx_AAZ (sP, 0); */
		/* } */

		int varmap[P->nvar];
		int j = 0;

		for (int i = 0; i < mvarP; i++) {
			varmap[i] = -1;
		}

		for (int i = mvarP; i < P->nvar; i++) {
			varmap[i] = j;
			j++;
		}

		shrinkAndReorderVars_AAZ (sP, varmap, P->nvar);
		isShrinked = 1;
    } else {
	sP = deepCopyPolynomial_AAZ (P);
    }

    AltArrZ_t* tmp = NULL;
    AltArrZ_t* tmpi = NULL;
    RecArrZ_t* recP = convertToRecursiveArrayZ (sP);
    RecArrElemZ_t* elems = recP->elems;
    int rSize = recP->size;

    int idx = 1;
    AltArrZ_t* coefGCD = deepCopyPolynomial_AAZ (convertFromAAZElemToAAZ
						 (elems[0].coef,  elems[0].coefSize,
						  sP->nvar, recP->unpacked));

    while (idx < rSize){
	if (coefGCD != NULL && coefGCD->size == 1 &&
	    leadingVariable_AAZ (coefGCD) == -1 &&
	    !mpz_cmp(coefGCD->elems->coef, one)){
	    break;
	}

	tmpi = deepCopyPolynomial_AAZ (convertFromAAZElemToAAZ (elems[idx].coef,
								elems[idx].coefSize,
								sP->nvar, recP->unpacked));

	tmp = gcd_AAZ (coefGCD, tmpi);
	freePolynomial_AAZ (coefGCD);
	freePolynomial_AAZ (tmpi);
	coefGCD = tmp;
	++idx;
    }

    if (isShrinked && coefGCD != NULL && coefGCD->size != 0){
	expandNumVarsLeft_AAZ (coefGCD, P->nvar);
    }
    if (coefGCD != NULL && coefGCD->size != 0) {
	return coefGCD;
    }

    return makeConstPolynomial_AAZ (1, P->nvar, one);
}



/******************
 * Square Free
 *****************/

/**
 * Computes the square free part of a polynomial.
 */
AltArrZ_t* squareFreePart_AAZ (AltArrZ_t* aa, int nvar)
{

	/* char* ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */

	/* fprintf (stderr, "[SMZP-Ali-SFP] aa := "); */
	/* printPoly_AAZ_unpk(stderr, aa, ch, nvar); */
	/* fprintf (stderr, "\n"); */

	/* fprintf (stderr, "[SMZP-Ali] input polynomial in squareFreePart_AAZ := "); */
	/* printPoly_AAZ_unpk(stderr, aa, ch, nvar); */
	/* fprintf (stderr, "\n"); */

    if (aa == NULL || aa->size == 0) {
		return NULL;
    }

    mpz_t one;
    mpz_init(one);
    mpz_set_si(one, 1l);

    AltArrZ_t* fact = makeConstPolynomial_AAZ(1, nvar, one);
	if (nvar == 0) {
		return fact;
	}

	int lvar = leadingVariable_AAZ (aa);
	if (lvar < 0) {
		return fact;
	}

	AltArrZ_t* primpart = deepCopyPolynomial_AAZ (aa);
    AltArrZ_t* cont = NULL;
    AltArrZ_t* diff = NULL;
    AltArrZ_t* g = NULL;
    AltArrZ_t* next = NULL;
    AltArrZ_t* swappp = NULL;


	int varmap[nvar];
	int nnvar = nvar;
	int isShrinked = 0;
	int j = 0;

	if (lvar > 0) {

		for (int i = 0; i < lvar; i++) {
			varmap[i] = -1;
		}

		for (int i = lvar; i < nvar; i++) {
			varmap[i] = j;
			j++;
		}

		shrinkAndReorderVars_AAZ (primpart, varmap, nvar);
		nnvar = nvar-lvar;

		for (int i = 0; i < nnvar; i++) {
			varmap[i] = i;
		}

		isShrinked = 1;

	} else {

		for (int i = 0; i < nvar; ++i) {
			varmap[i] = i;
		}
	}

	/* for (int i = 0; i < nnvar; i++) { */
	/* 	fprintf (stderr, "varmap[%d] = %d\n", i, varmap[i]); */
	/* } */

	/* fprintf (stderr, "nnvar = %d\n", nnvar); */
	/* fprintf (stderr, "mdeg(primpart) = %d\n", mainDegree_AAZ(primpart)); */

	/* fprintf (stderr, "[SMZP-Ali] primpart after shrink := "); */
	/* printPoly_AAZ_unpk(stderr, primpart, ch, nnvar); */
	/* fprintf (stderr, "\n"); */

    for (int i = 0; i < nnvar; ++i) {
		varmap[0] = i;
		varmap[i] = 0;
		reorderVars_AAZ (primpart, varmap, nnvar); // reorder primpart w.r.t. i-th variable

		// primpart should have non-zero mvar:
		if (leadingVariable_AAZ(primpart)) {
			continue;
		}

		// swappp is prim. part of poly w.r.t. i-th variable
		// cont is next primpart
		swappp = mainPrimitiveFactorization_AAZ (primpart, &cont);

		freePolynomial_AAZ (primpart);

		// diff = d(primpart)/d(x_i)
		diff = derivative_AAZ (swappp, 0, 1);


		/* fprintf (stderr, "[SMZP-Ali] swappp := "); */
		/* printPoly_AAZ_unpk(stderr, swappp, ch, nnvar); */
		/* fprintf (stderr, "\n"); */


		/* fprintf (stderr, "[SMZP-Ali] diff := "); */
		/* printPoly_AAZ_unpk(stderr, diff, ch, nnvar); */
		/* fprintf (stderr, "\n"); */

		// g = gcd(swappp, D(swappp))
		g = gcd_AAZ (swappp, diff);

		/* fprintf (stderr, "[SMZP-Ali] gcd := "); */
		/* printPoly_AAZ_unpk(stderr, g, ch, nnvar); */
		/* fprintf (stderr, "\n"); */


		freePolynomial_AAZ (diff);

		// next = swappp / g
		exactDividePolynomials_AAZ (swappp, g, &next, nnvar);

		freePolynomial_AAZ (g);
		freePolynomial_AAZ (swappp);

		varmap[0] = i;
		varmap[i] = 0;
		reorderVars_AAZ (next, varmap, nnvar);

		if (next != NULL && next->size) {
			if (mpz_sgn(next->elems->coef) < 0) {
				negatePolynomial_AAZ (next);
			}

			fact = multiplyPolynomials_AAZ_inp (fact, next, nnvar); // update fact
			freePolynomial_AAZ (next);
		}

		if (leadingVariable_AAZ(cont) < 0) {
			break;
		} else {
			primpart = cont;
		}
    }

	if (isShrinked) {
		expandNumVarsLeft_AAZ (fact,nvar);
	}

    return fact;
}



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
void squareFree_AAZ(const AltArrZ_t* aa, mpz_t u, AltArrZ_t*** facts_p, degree_t** exps_p, int* nfacts) {

	if (facts_p == NULL || nfacts == NULL) {
		integralContent_AAZ(aa, u);
		if (nfacts != NULL) {
			*nfacts = 0;
		}
		return;
	}

	if (isZero_AAZ(aa)) {
		mpz_set_ui(u, 0ul);
		if (nfacts != NULL) {
			*nfacts = 0;
		}
		return;
	}

	if (isConstant_AAZ(aa) || aa->nvar == 0) {
		mpz_set(u, aa->elems->coef);
		if (nfacts != NULL) {
			*nfacts = 0;
		}
		return;
	}

	AltArrZ_t** facts = *facts_p;
	int factorIdx = 0;
	int factorAlloc = facts == NULL ? 10 : *nfacts;
	if (facts == NULL) {
		facts = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*factorAlloc);
	} else {
		facts = (AltArrZ_t**) realloc(facts, sizeof(AltArrZ_t*)*factorAlloc);
	}

	degree_t* exps;
	if (exps_p == NULL) {
		exps = (degree_t*) malloc(sizeof(degree_t)*factorAlloc);
	} else {
		exps = *exps_p;
		if (exps == NULL) {
			exps = (degree_t*) malloc(sizeof(degree_t)*factorAlloc);
			*exps_p = exps;
		}
	}

// char* syms[] = {"x", "y", "z", "s", "t", "u", "v", "w"};

	int k = 1;
	int nvar = aa->nvar;

	mpz_t tmpCont;
	mpz_init(tmpCont);
	AltArrZ_t* primPart1 = primitivePartAndContent_AAZ(aa, u);
	int mvar = leadingVariable_AAZ(primPart1);

	int activeForContent[nvar];
	for (int i = 0; i < nvar; ++i) {
		activeForContent[i] = 0;
	}
	activeForContent[mvar] = 1;
	AltArrZ_t* content = contentInVars_AAZ(primPart1, activeForContent);
	AltArrZ_t* primPart = NULL;
	exactDividePolynomials_AAZ(primPart1, content, &primPart, primPart1->nvar);
	freePolynomial_AAZ(primPart1);
	// AltArrZ_t* primPart = mainPrimitiveFactorization_AAZ(primPart1, &content);

	for (int varIdx = 0; varIdx < nvar; ++varIdx) {
		if ( (isConstant_AAZ(primPart) && mpz_cmp_ui(primPart->elems->coef, 1ul) != 0)
		     || partialDegreeTerm_AAZ(primPart, 0, mvar) == 1) {

			if (factorIdx >= factorAlloc) {
				factorAlloc <<= 1;
				facts = (AltArrZ_t**) realloc(facts, sizeof(AltArrZ_t*)*factorAlloc);
				exps = (degree_t*) realloc(exps, sizeof(degree_t)*factorAlloc);
			}
			primitivePartAndContent_AAZ_inp(primPart, tmpCont);
			mpz_mul(u, u, tmpCont);
			if(!isConstant_AAZ(primPart)) {
				facts[factorIdx] = primPart;
				primPart = NULL;
				exps[factorIdx] = 1;
				++factorIdx;
			}

		} else if (!isConstant_AAZ(primPart)) {
			AltArrZ_t* dx = derivative_AAZ(primPart, mvar, 1);
			AltArrZ_t* g = gcd_AAZ(primPart, dx);
			primitivePart_AAZ_inp(g);
			freePolynomial_AAZ(dx);

			AltArrZ_t* next = NULL;
			exactDividePolynomials_AAZ (primPart, g, &next, primPart->nvar);
			// if (next != NULL && mpz_sgn(next->elems->coef) < 0) {
			// 	negatePolynomial_AAZ(next);
			// }

			// if (isConstant_AAZ(g)) {
			// 	AltArrZ_t* remFact = NULL;
			// 	AltArrZ_t* comFact = commonFactor_AAZ(next, &remFact);
			// 	if (!isConstant_AAZ(comFact)) {
			// 		freePolynomial_AAZ(next);
			// 		next = remFact;

			// 		degree_t comDegs[nvar];
			// 		partialDegreesTerm_AAZ(comFact, 0, comDegs);
			// 		for (int j = 0; j < nvar; ++j) {
			// 			if (comDegs[j] != 0) {
			// 				AltArrZ_t* temp = makeConstIntPolynomial_AAZ(1, nvar, 1l);
			// 				setExponentTerm_AAZ_inp(temp, 0, 1, j);

			// 				if (factorIdx >= factorAlloc) {
			// 					factorAlloc <<= 1;
			// 					facts = (AltArrZ_t**) realloc(facts, sizeof(AltArrZ_t*)*factorAlloc);
			// 					exps = (degree_t*) realloc(exps, sizeof(degree_t)*factorAlloc);
			// 				}

			// 				facts[factorIdx] = temp;
			// 				exps[factorIdx] = comDegs[j];
			// 				++factorIdx;
			// 			}
			// 		}
			// 	}
			// 	freePolynomial_AAZ(comFact);
			// }

			k = 1;
			while (partialDegreeTerm_AAZ(g, 0, mvar) > 0) {
				AltArrZ_t* y = gcd_AAZ(next, g);
				primitivePart_AAZ_inp(y);
				if (!isExactlyEqual_AAZ(next, y)) {
					if (factorIdx >= factorAlloc) {
						factorAlloc <<= 1;
						facts = (AltArrZ_t**) realloc(facts, sizeof(AltArrZ_t*)*factorAlloc);
						exps = (degree_t*) realloc(exps, sizeof(degree_t)*factorAlloc);
					}

					facts[factorIdx] = NULL;
					exactDividePolynomials_AAZ(next, y, facts+factorIdx, nvar);

					primitivePartAndContent_AAZ_inp(facts[factorIdx], tmpCont);
					mpz_mul(u, u, tmpCont);

					if (isConstant_AAZ(facts[factorIdx])) {
						freePolynomial_AAZ(facts[factorIdx]);
						facts[factorIdx] = NULL;
					} else {
						exps[factorIdx] = k;
						++factorIdx;
					}
				}
				freePolynomial_AAZ(next);
				next = NULL;
				exactDividePolynomials_AAZ(g, y, &next, nvar);
				freePolynomial_AAZ(g);
				g = next;
				next = y;
				if (next != NULL && mpz_cmp_ui(next->elems->coef, 1l) < 0) {
					negatePolynomial_AAZ(next);
				}
				++k;
			}

			if (isConstant_AAZ(next)) {
				mpz_mul(u, u, next->elems->coef);
				freePolynomial_AAZ(next);
			} else {
				if (factorIdx >= factorAlloc) {
					factorAlloc <<= 1;
					facts = (AltArrZ_t**) realloc(facts, sizeof(AltArrZ_t*)*factorAlloc);
					exps = (degree_t*) realloc(exps, sizeof(degree_t)*factorAlloc);
				}

				primitivePartAndContent_AAZ_inp(next, tmpCont);
				mpz_mul(u, u, tmpCont);
				facts[factorIdx] = next;
				exps[factorIdx] = k;
				++factorIdx;
			}

			freePolynomial_AAZ(g);
		}//if not constant

		freePolynomial_AAZ(primPart);
		primPart = content;
		content = NULL;

		activeForContent[mvar] = 0;
		mvar = leadingVariable_AAZ(primPart);
		activeForContent[mvar] = 1;
		content = contentInVars_AAZ(primPart, activeForContent);
		exactDividePolynomials_AAZ(primPart, content, &primPart1, primPart->nvar);
		// primPart1 = mainPrimitiveFactorization_AAZ(primPart, &content);

		freePolynomial_AAZ(primPart);
		primPart = primPart1;
		primPart1 = NULL;
		// zmvar = leadingVariable_AAZ(primPart);

	}

	freePolynomial_AAZ(primPart);
	freePolynomial_AAZ(content);
	mpz_clear(tmpCont);

	*facts_p = facts;
	if (exps_p == NULL) {
		free(exps);
	}
	if (nfacts != NULL) {
		*nfacts = factorIdx;
	}
}

int biModularSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime )
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In modularBiSubresultantChain_DBZP, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}

	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	polysize_t min_Xdeg = MIN (a_degs[1], b_degs[1]);
	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;
	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	int primeIdx = 0;
	Prime_ptr Pptr[1];

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;
	int isSubresEq;

	// int isNULL = 0;
	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		// fprintf(stderr, "In biModularSubresultantChainZ, main-for-loop primeIdx = %d\n", primeIdx); // TEST
		// Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		*Pptr = prime64_ptr[primeIdx];
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			continue;
		}


		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		biSylvSubResultantInForm_spX (modA, a_degs, modB, b_degs, &modSubres, Pptr);
		// fprintf(stderr, "modSubres->n = %ld\n", modSubres->n); // TEST

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		// ***** CRT using gross comparison ***** //
		// // mpz_set (res_work_prev, subres_work[0]->coefs[0]);
		// // subres_work_prev = subres_work
		// if (subres_work_prev == NULL) {
		// 	subres_work_prev = makeBiSubresultant_DBZP (subres_work->n);
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		subres_work_prev->deg[j][0] = subres_work->deg[j][0];
		// 		subres_work_prev->deg[j][1] = subres_work->deg[j][1];
		// 		subres_work_prev->size[j]   = subres_work->size[j];
		// 		subres_work_prev->coefs[j] = (mpz_t*) malloc (subres_work_prev->size[j]*sizeof(mpz_t));
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_init_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// } else {
		// 	if (subres_work->n != subres_work_prev->n) {
		// 		fprintf (stderr, "subres_work->n != subres_work_prev\n");
		// 		exit (1);
		// 	}
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		if (subres_work_prev->size[j] != subres_work->size[j]) {
		// 			fprintf (stderr, "subres_work_prev->size[%ld] != subres_work->size[%ld]\n", j, j);
		// 			exit (1);
		// 		}
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// }

		// for (polysize_t k = 0; k < modSubres->n; k++) {
		// 	if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
		// 		continue;
		// 	}
		// 	polysize_t lt = modSubres->size[k];
		// 	subres_work->deg[k][0] = modSubres->deg[k][0];
		// 	subres_work->deg[k][1] = modSubres->deg[k][1];
		// 	subres_work->size[k] = modSubres->size[k]; // update lt
		// 	for (polysize_t i = 0; i < lt; ++i) {
		// 		coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
		// 		// optimized-CRT (Garner's algorithm)

		// 		mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
		// 		mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
		// 			mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		}
		// 	}
		// }
		// // subres_work == subres_work_prev ?
		// isSubresEq = isBiSubresultantEquals_DBZP (subres_work, subres_work_prev);

		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)

				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);

			}
		}

		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			if (no_prime != NULL) {*no_prime = primeIdx+1;}
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (a);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free (modA); free (modB);
			// reverse the order s.t. tail->head:
			if (sz<2) {
				*Subres = head;
			} else {
				AltArrsZ_t* curr = head, *prev = NULL, *next = NULL;
				while (curr != NULL){
					next = curr->next;
					curr->next = prev;
					prev = curr;
					curr = next;
				}
				*Subres = prev;
			}
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free (modA); free (modB);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	// fprintf (stderr, "In biModularSubresultantChain_DBZP, all primes failed\n"); // TEST
	*Subres = NULL;
	*chain_size = 0;
	return 0;
}

int biModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime)
{

	// const char* sym[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST
	// fprintf (stderr, "In biModularFFTSubresultantChainZ a := ");
	// printPoly_AAZ (stderr, a, sym, a->nvar);
	// fprintf (stderr, "\nIn biModularFFTSubresultantChainZ b := ");
	// printPoly_AAZ (stderr, b, sym, b->nvar);
	// fprintf(stderr, "\n");

	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In biModularFFTSubresultantChainZ, expected only bivariate polynomials.\n");
		exit (1);
	}

	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);

	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	biSubresPr_t* modSubres = NULL;
	unsigned long uprime, coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389;
	const int max_no_failed_omega = 3;
	int isGoodOmega, isSubresEq, no_failed_omega = 0;
	int K, e;
	polysize_t n = _compute_N_gt_By(a_degs, b_degs, &K, &e);

	usfixn64 p_u64, w, w_inv, n_inv;
    montgomery_triple P;
	// fprintf (stderr, "start testing primes...\n"); //TEST
	// int prime_upperbound = no_prime != 0 ? no_prime : n_prime64_ptr; // TODO: don't need this... just for testing :)
	for (; primeIdx < prime_table_size; ++primeIdx) { // primeIdx < n_prime64_ptr
		// *Pptr = prime64_ptr[primeIdx];
		// @note fourier_primes_u64_table is not in-order to use directly in CRT
		// @note subset_fourier_primes_u64_table is a in-order subset of fourier_primes_u64_table
		// @note if CRT fails for large inputs, check the prime list at first!
		// TODO: update fourier_primes_u64_table s.t. FFT-routines in GFPF won't break!
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			free(Pptr);
			continue;
		}

		p_u64 = (unsigned long long) Pptr->prime;
		init_montgomery_triple (&P, p_u64);
		compute_nth_root_of_unity_for_small_prime (p_u64, &w, n);

		// @note we check it now in bivarModularSubResInForm_withFFT_spX
		// @note it's not efficient when we check it here.

		gmp_inv_mod_p_u64 (&w_inv, w, p_u64);
	    gmp_inv_mod_p_u64 (&n_inv, n, p_u64);

		convertIn_GLOBAL_ptr(&w, &P);
	    convertIn_GLOBAL_ptr(&w_inv, &P);
    	convertIn_GLOBAL_ptr(&n_inv, &P);

		// if good prime and omega:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert b to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		isGoodOmega = bivarModularSubResInForm_withFFT_spX(modA, a_degs, modB, b_degs, n, K, e, w, w_inv, n_inv, &modSubres, Pptr);

		if(!isGoodOmega) {
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				free(Pptr);
				continue;
			} else {
				// free aZ, bZ
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				free(Pptr);
				return 0;
			}
		}

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
			}
		}


		free(Pptr);
		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			if (no_prime != NULL) {*no_prime = primeIdx+1;}
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (a);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			// reverse the order s.t. tail->head:
			if (sz<2) {
				*Subres = head;
			} else {
				AltArrsZ_t* curr = head, *prev = NULL, *next = NULL;
				while (curr != NULL){
					next = curr->next;
					curr->next = prev;
					prev = curr;
					curr = next;
				}
				*Subres = prev;
			}
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	// fprintf (stderr, "In biModularFFTSubresultantChain_DBZP, all primes failed\n"); // TEST
	*Subres = NULL;
	*chain_size = 0;
	return 0;
}

int biModularFFT4SubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime)
{
	// const char* sym[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST
	// fprintf (stderr, "In biModularFFTSubresultantChainZ a := ");
	// printPoly_AAZ (stderr, a, sym, a->nvar);
	// fprintf (stderr, "\nIn biModularFFTSubresultantChainZ b := ");
	// printPoly_AAZ (stderr, b, sym, b->nvar);
	// fprintf(stderr, "\n");

	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In biModularFFTSubresultantChainZ, expected only bivariate polynomials.\n");
		exit (1);
	}

	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);

	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	biSubresPr_t* modSubres = NULL;
	unsigned long uprime, coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389;
	const int max_no_failed_omega = 3;
	int isGoodOmega, isSubresEq, no_failed_omega = 0;

	polysize_t md =  b_degs[0] * a_degs[1] + b_degs[1] * a_degs[0] + 1;
	polysize_t r = Log2_AAZ(md)+2;
	polysize_t n = 1 << r;
	polysize_t n2 = n << 1;
	elem_t n_inv;
	elem_t *Omega, *OmegaInv;
	elem_t* _omega, *_omegaInv;

	// fprintf (stderr, "start testing primes...\n"); //TEST
	// int prime_upperbound = no_prime != 0 ? no_prime : n_prime64_ptr; // TODO: don't need this... just for testing :)
	for (; primeIdx < prime_table_size; ++primeIdx) { // primeIdx < n_prime64_ptr
		// *Pptr = prime64_ptr[primeIdx];
		// @note fourier_primes_u64_table is not in-order to use directly in CRT
		// @note subset_fourier_primes_u64_table is a in-order subset of fourier_primes_u64_table
		// @note if CRT fails for large inputs, check the prime list at first!
		// TODO: update fourier_primes_u64_table s.t. FFT-routines in GFPF won't break!
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			free(Pptr);
			continue;
		}

		Omega = (elem_t*) calloc (n2<<1, sizeof (elem_t));
		OmegaInv = (elem_t*) calloc (n2<<1, sizeof (elem_t));
		RootsTable2_spf(n2, r+1, Omega, OmegaInv, Pptr);
		_omega = Omega; _omegaInv = OmegaInv;
		Omega += n2; OmegaInv += n2;
	    n_inv = smallprimefield_inv(smallprimefield_convert_in(n, Pptr), Pptr);

		// if good prime and omega:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert b to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// fprintf(stderr, "running ...\n");
		// call subresultant mod this prime
		isGoodOmega = bivarModularSubResInForm_withFFT4_spX(modA, a_degs, modB, b_degs, n, r, Omega, OmegaInv, n_inv, &modSubres, Pptr);
		// fprintf(stderr,"done...\n");	

		free(_omega);
		free(_omegaInv);

		if(!isGoodOmega) {
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				free(Pptr);
				continue;
			} else {
				// free aZ, bZ
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				free(Pptr);
				return 0;
			}
		}

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
			}
		}


		free(Pptr);
		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			if (no_prime != NULL) {*no_prime = primeIdx+1;}
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (a);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			// reverse the order s.t. tail->head:
			if (sz<2) {
				*Subres = head;
			} else {
				AltArrsZ_t* curr = head, *prev = NULL, *next = NULL;
				while (curr != NULL){
					next = curr->next;
					curr->next = prev;
					prev = curr;
					curr = next;
				}
				*Subres = prev;
			}
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	// fprintf (stderr, "In biModularFFTSubresultantChain_DBZP, all primes failed\n"); // TEST
	*Subres = NULL;
	*chain_size = 0;
	return 0;
}


// descending
int hgcdBiModularSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size)
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In hgcdBiModularSubresultantChainZ, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}
	AltArrsZ_t *node, *head=NULL, *tail=NULL;
	if (isZero_AAZ (a)) {
		if (isZero_AAZ (b)) {
			*ksubres = NULL;
			*chain_size = 0;
		} else {
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			head = node;
			tail = node;
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = NULL;
			node->next = NULL;
			tail->next = node;
			tail = node; // todo
			*chain_size = 2;
			*ksubres = head;
		}
		return 1;
	} else if (isZero_AAZ (b)) {
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ (a);
		node->next = NULL;
		head = node;
		tail = node;
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = NULL;
		node->next = NULL;
		tail->next = node;
		tail = node; // todo
		*chain_size = 2;
		*ksubres = head;
		return 1;
	} else if (mainLeadingDegree_AAZ (a) < mainLeadingDegree_AAZ (b)) {
		return hgcdBiModularSubresultantChainZ (b, a, kk, ksubres, chain_size);
	}
	int sz = 1;
	polysize_t min_Xdeg = MIN (mainLeadingDegree_AAZ(a), mainLeadingDegree_AAZ(b));
	if (kk < 0 || kk >= min_Xdeg) {
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ(a);
		node->next = NULL;
		head = node;
		tail = node;
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ(b);
		node->next = NULL;
		tail->next = node;
		tail = node;
        *chain_size = 2;
		*ksubres = head;
        return 1;
    }
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	// polysize_t min_Xdeg = MIN (a_degs[1], b_degs[1]);
	// fprintf (stderr, "deg(aZ) = [y=%ld, x=%ld]\n", a_degs[0], a_degs[1]);
	// fprintf (stderr, "deg(bZ) = [y=%ld, x=%ld]\n", b_degs[0], b_degs[1]);
	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;
	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	int primeIdx = 0;
	Prime_ptr Pptr[1];

	AltArrZ_t* tmp;
	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;
	int isSubresEq;

	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		*Pptr = prime64_ptr[primeIdx];
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			continue;
		}

		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		// ascending
		hgcdBiSylvSubResultantInForm_spX (modA, a_degs, modB, b_degs, kk, &modSubres, Pptr);

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		// // mpz_set (res_work_prev, subres_work[0]->coefs[0]);
		// // subres_work_prev = subres_work
		// if (subres_work_prev == NULL) {
		// 	subres_work_prev = makeBiSubresultant_DBZP (subres_work->n);
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		subres_work_prev->deg[j][0] = subres_work->deg[j][0];
		// 		subres_work_prev->deg[j][1] = subres_work->deg[j][1];
		// 		subres_work_prev->size[j]   = subres_work->size[j];
		// 		subres_work_prev->coefs[j] = (mpz_t*) malloc (subres_work_prev->size[j]*sizeof(mpz_t));
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_init_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// } else {
		// 	if (subres_work->n != subres_work_prev->n) {
		// 		fprintf (stderr, "subres_work->n != subres_work_prev->n\n");
		// 		exit (1);
		// 	}
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		if (subres_work_prev->size[j] != subres_work->size[j]) {
		// 			fprintf (stderr, "subres_work_prev->size[%ld] != subres_work->size[%ld]\n", j, j);
		// 			exit (1);
		// 		}
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// }

		// for (polysize_t k = 0; k < modSubres->n; k++) {
		// 	if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
		// 		continue;
		// 	}
		// 	polysize_t lt = modSubres->size[k];
		// 	subres_work->deg[k][0] = modSubres->deg[k][0];
		// 	subres_work->deg[k][1] = modSubres->deg[k][1];
		// 	subres_work->size[k] = modSubres->size[k]; // update lt
		// 	for (polysize_t i = 0; i < lt; ++i) {
		// 		coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
		// 		// optimized-CRT (Garner's algorithm)
		// 		mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
		// 		mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
		// 			mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		}
		// 	}
		// }
		// // subres_work == subres_work_prev?
		// int isSubresEq = isBiSubresultantEquals_DBZP (subres_work, subres_work_prev);

		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
			}
		}

		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			// // // To add inputs polynomials a, b to the list... for testing only!
			// @note must be commented!
			// node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			// node->poly = deepCopyPolynomial_AAZ (b);
			// node->next = NULL;
			// tail->next = node;
			// tail = node;
			// sz++;
			// node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			// node->poly = deepCopyPolynomial_AAZ (a);
			// node->next = NULL;
			// tail->next = node;
			// tail = node;
			// sz++;
			// reverse the order s.t. tail->head:
			// if (sz<2) {
			// 	*ksubres = head;
			// } else {
			// 	AltArrsZ_t* curr = head, *prev = NULL, *next = NULL;
			// 	while (curr != NULL){
			// 		next = curr->next;
			// 		curr->next = prev;
			// 		prev = curr;
			// 		curr = next;
			// 	}
			// 	*ksubres = prev;
			// }
			// // // // //
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			*ksubres = head;
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}

	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	// fprintf (stderr, "In hgcdBiModularSubresultantChain_DBZP, all primes failed\n"); // TEST
	*ksubres = NULL;
	*chain_size = 0;
	return 0;
}

int hgcdBiModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size)
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In hgcdBiModularFFTSubresultantChainZ, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}
	AltArrsZ_t *node, *head=NULL, *tail=NULL;
	if (isZero_AAZ (a)) {
		if (isZero_AAZ (b)) {
			*ksubres = NULL;
			*chain_size = 0;
		} else {
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			head = node;
			tail = node;
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = NULL;
			node->next = NULL;
			tail->next = node;
			tail = node; // todo
			*chain_size = 2;
			*ksubres = head;
		}
		return 1;
	} else if (isZero_AAZ (b)) {
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ (a);
		node->next = NULL;
		head = node;
		tail = node;
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = NULL;
		node->next = NULL;
		tail->next = node;
		tail = node; // todo
		*chain_size = 2;
		*ksubres = head;
		return 1;
	} else if (mainLeadingDegree_AAZ (a) < mainLeadingDegree_AAZ (b)) {
		return hgcdBiModularFFTSubresultantChainZ (b, a, kk, ksubres, chain_size);
	}
	int sz = 1;
	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	polysize_t min_Xdeg = MIN (a_degs[1], b_degs[1]);

	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389;
	int isGoodOmega, isSubresEq, no_failed_omega;
	const int max_no_failed_omega = 5;

	int K, e;
	polysize_t n = _compute_N_gt_By(a_degs, b_degs, &K, &e);

	usfixn64 p_u64, w, w_inv, n_inv;
    montgomery_triple P;
	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < prime_table_size; ++primeIdx) {

		// @note see the biModularFFTSubresultantChainZ for more details...
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			free(Pptr);
			continue;
		}

		p_u64 = (unsigned long long) Pptr->prime;
		init_montgomery_triple (&P, p_u64);
		compute_nth_root_of_unity_for_small_prime (p_u64, &w, n);

		// @note we check it now in bivarModularSubResInForm_withFFT_spX
		// @note see the biModularFFTSubresultantChainZ for more details...

		gmp_inv_mod_p_u64 (&w_inv, w, p_u64);
	    gmp_inv_mod_p_u64 (&n_inv, n, p_u64);

		convertIn_GLOBAL_ptr(&w, &P);
	    convertIn_GLOBAL_ptr(&w_inv, &P);
    	convertIn_GLOBAL_ptr(&n_inv, &P);

		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		isGoodOmega = bivarModularHGCDSubResInForm_withFFT_spX(modA, a_degs, modB, b_degs, kk, n, K, e, w, w_inv, n_inv, &modSubres, Pptr);

		if(!isGoodOmega) {
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				continue;
			} else {
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				free(Pptr);
				return 0;
			}
	   }

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		// fprintf(stderr, "modSubres->n = %ld\n", modSubres->n);
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		// // subres_work_prev ?= subres_work
		// if (subres_work_prev == NULL) {
		// 	subres_work_prev = makeBiSubresultant_DBZP (subres_work->n);
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		subres_work_prev->deg[j][0] = subres_work->deg[j][0];
		// 		subres_work_prev->deg[j][1] = subres_work->deg[j][1];
		// 		subres_work_prev->size[j]   = subres_work->size[j];
		// 		subres_work_prev->coefs[j] = (mpz_t*) malloc (subres_work_prev->size[j]*sizeof(mpz_t));
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_init_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// } else {
		// 	if (subres_work->n != subres_work_prev->n) {
		// 		fprintf (stderr, "subres_work->n != subres_work_prev\n");
		// 		exit (1);
		// 	}
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		if (subres_work_prev->size[j] != subres_work->size[j]) {
		// 			fprintf (stderr, "subres_work_prev->size[%ld] != subres_work->size[%ld]\n", j, j);
		// 			exit (1);
		// 		}
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// }

		// for (polysize_t k = 0; k < modSubres->n; k++) {
		// 	if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
		// 		continue;
		// 	}
		// 	polysize_t lt = modSubres->size[k];
		// 	subres_work->deg[k][0] = modSubres->deg[k][0];
		// 	subres_work->deg[k][1] = modSubres->deg[k][1];
		// 	subres_work->size[k] = modSubres->size[k]; // update lt
		// 	for (polysize_t i = 0; i < lt; ++i) {
		// 		coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
		// 		// optimized-CRT (Garner's algorithm)
		// 		mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
		// 		mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
		// 			mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		}
		// 	}
		// }
		// // subres_work == subres_work_prev ?
		// isSubresEq = isBiSubresultantEquals_DBZP (subres_work, subres_work_prev);


		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
			}
		}

		free(Pptr);
		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			*ksubres = head;
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	// fprintf (stderr, "In hgcdBiModularFFTSubresultantChain_DBZP, all primes failed\n"); // TEST
	*ksubres = NULL;
	*chain_size = 0;
	return 0;
}

int hgcdBiModularFFT4SubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size)
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In hgcdBiModularFFTSubresultantChainZ, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}
	AltArrsZ_t *node, *head=NULL, *tail=NULL;
	if (isZero_AAZ (a)) {
		if (isZero_AAZ (b)) {
			*ksubres = NULL;
			*chain_size = 0;
		} else {
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			head = node;
			tail = node;
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = NULL;
			node->next = NULL;
			tail->next = node;
			tail = node; // todo
			*chain_size = 2;
			*ksubres = head;
		}
		return 1;
	} else if (isZero_AAZ (b)) {
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ (a);
		node->next = NULL;
		head = node;
		tail = node;
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = NULL;
		node->next = NULL;
		tail->next = node;
		tail = node; // todo
		*chain_size = 2;
		*ksubres = head;
		return 1;
	} else if (mainLeadingDegree_AAZ (a) < mainLeadingDegree_AAZ (b)) {
		return hgcdBiModularFFT4SubresultantChainZ (b, a, kk, ksubres, chain_size);
	}
	int sz = 1;
	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	polysize_t min_Xdeg = MIN (a_degs[1], b_degs[1]);

	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389;
	int isGoodOmega, isSubresEq, no_failed_omega;
	const int max_no_failed_omega = 5;

	polysize_t md =  b_degs[0] * a_degs[1] + b_degs[1] * a_degs[0] + 1;
	polysize_t r = Log2_AAZ(md)+2;
	polysize_t n = 1 << r;
	polysize_t n2 = n << 1;
	elem_t n_inv;
	elem_t *Omega, *OmegaInv;
	elem_t* _omega, *_omegaInv;

	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < prime_table_size; ++primeIdx) {

		// @note see the biModularFFTSubresultantChainZ for more details...
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			free(Pptr);
			continue;
		}

		Omega = (elem_t*) calloc (n2<<1, sizeof (elem_t));
		OmegaInv = (elem_t*) calloc (n2<<1, sizeof (elem_t));
		RootsTable2_spf(n2, r+1, Omega, OmegaInv, Pptr);
		_omega = Omega; _omegaInv = OmegaInv;
		Omega += n2; OmegaInv += n2;
	    n_inv = smallprimefield_inv(smallprimefield_convert_in(n, Pptr), Pptr);

		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		isGoodOmega = bivarModularHGCDSubResInForm_withFFT4_spX(modA, a_degs, modB, b_degs, kk, n, r, Omega, OmegaInv, n_inv, &modSubres, Pptr);

		free(_omega);
		free(_omegaInv);

		if(!isGoodOmega) {
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				continue;
			} else {
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				free(Pptr);
				return 0;
			}
	   }

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		// fprintf(stderr, "modSubres->n = %ld\n", modSubres->n);
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
			}
		}

		free(Pptr);
		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			*ksubres = head;
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	// fprintf (stderr, "In hgcdBiModularFFTSubresultantChain_DBZP, all primes failed\n"); // TEST
	*ksubres = NULL;
	*chain_size = 0;
	return 0;
}


int regularGCDBiModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size, 
											polysize_t *mdegs, specAQRArray_spXY_t **bspecInfoArray, int *info_size, int results_mode)
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In regularGCDBiModularFFTSubresultantChainZ, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}

	AltArrsZ_t *node, *head=NULL, *tail=NULL;
	if (isZero_AAZ (a)) {
		if (isZero_AAZ (b)) {
			*ksubres = NULL;
			*chain_size = 0;
		} else {
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			head = node;
			tail = node;
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = NULL;
			node->next = NULL;
			tail->next = node;
			tail = node; // todo
			*chain_size = 2;
			*ksubres = head;
		}
		return 1;
	} else if (isZero_AAZ (b)) {
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ (a);
		node->next = NULL;
		head = node;
		tail = node;
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = NULL;
		node->next = NULL;
		tail->next = node;
		tail = node; // todo
		*chain_size = 2;
		*ksubres = head;
		return 1;
	} else if (mainLeadingDegree_AAZ (a) < mainLeadingDegree_AAZ (b)) {
		return regularGCDBiModularFFTSubresultantChainZ (b, a, kk, ksubres, chain_size, mdegs, bspecInfoArray, info_size, results_mode);
	}
	int sz = 1;
	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	polysize_t min_Xdeg = MIN (a_degs[1], b_degs[1]);

	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr, tmpZ; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, tmpZ, NULL);

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389; // TODO: update prime subset.. 
	int isGoodOmega, isSubresEq, no_failed_omega;
	const int max_no_failed_omega = 5;

	int infoIdx = 0, 
		infoAlloc = 25,
		isInfoNull = 0;
	specAQR_spXY_t ** bspecInfo;
	if (*info_size < 1 || bspecInfoArray == NULL) {
		bspecInfo = (specAQR_spXY_t **) malloc (infoAlloc * sizeof(specAQR_spXY_t*));
		for (int i = 0; i < infoAlloc; ++i) {
			bspecInfo[i] = (specAQR_spXY_t *) malloc (sizeof(specAQR_spXY_t));
		}
		isInfoNull = 1;
	} else {
		infoAlloc = *info_size;
		bspecInfo = (*bspecInfoArray)->bspecArray;
	}

	int K, e;
	polysize_t n = _compute_N_gt_By(a_degs, b_degs, &K, &e);

	usfixn64 p_u64, w, w_inv, n_inv;
    montgomery_triple P;
	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < prime_table_size; ++primeIdx) {
		// fprintf(stderr, "primeIdx = %d\n", primeIdx);
		// @note see the biModularFFTSubresultantChainZ for more details...
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			free(Pptr);
			continue;
		}
		p_u64 = (unsigned long long) Pptr->prime;
		init_montgomery_triple (&P, p_u64);
		if (isInfoNull) {
			compute_nth_root_of_unity_for_small_prime (p_u64, &w, n);
			bspecInfo[infoIdx]->W = w;
		} else {
			if (infoIdx > *info_size) {
				fprintf(stderr, "SMZP Error: exceed the upper-bound for the number of primes in lazy computation of subresultants\n");
				exit(EXIT_FAILURE);
			}

			w = bspecInfo[infoIdx]->W;
		}

		// @note we check it now in bivarModularSubResInForm_withFFT_spX
		// @note see the biModularFFTSubresultantChainZ for more details...

		gmp_inv_mod_p_u64 (&w_inv, w, p_u64);
	    gmp_inv_mod_p_u64 (&n_inv, n, p_u64);

		convertIn_GLOBAL_ptr(&w, &P);
	    convertIn_GLOBAL_ptr(&w_inv, &P);
    	convertIn_GLOBAL_ptr(&n_inv, &P);

		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		if (isInfoNull) {
			bspecInfo[infoIdx]->An = createSpecA_spXY (n, a_tdeg); 
			bspecInfo[infoIdx]->Qn = createSpecQ_spXY (n, a_tdeg); 
			bspecInfo[infoIdx]->Rn = createSpecR_spXY (n, a_tdeg); 
		}
		// call subresultant mod this prime
		// polysize_t mdegs[2] = {0, 0};
		isGoodOmega = regularGCDbivarModularHGCDSubResInForm_withFFT_spX (modA, a_degs, modB, b_degs, kk, n, K, e, w, w_inv, n_inv, 
																		&modSubres, mdegs, results_mode, bspecInfo[infoIdx], Pptr);

		infoIdx++;
		if (infoIdx == infoAlloc) {
			infoAlloc += 25;
			specAQR_spXY_t ** tmpspecInfo = (specAQR_spXY_t**) malloc (infoAlloc*sizeof(specAQR_spXY_t*));
			for (int i = 0; i < infoIdx; ++i) {
				tmpspecInfo[i] = bspecInfo[i];
			}
			for (int i = infoIdx; i < infoAlloc; ++i) {
				tmpspecInfo[i] = (specAQR_spXY_t*) malloc (sizeof(specAQR_spXY_t));
			}
			free (bspecInfo);
			bspecInfo = tmpspecInfo;
 		}

		if(!isGoodOmega) {
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				continue;
			} else {
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				free(Pptr);
				return 0;
			}
	    }

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		// fprintf(stderr, "modSubres->n = %ld\n", modSubres->n);
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		isSubresEq = 1;
		elem_t coef_out;
		for (polysize_t k = 0; k < modSubres->n; ++k) {
			if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
				continue;
			}
			polysize_t lt = modSubres->size[k];
			subres_work->deg[k][0] = modSubres->deg[k][0];
			subres_work->deg[k][1] = modSubres->deg[k][1];
			subres_work->size[k] = modSubres->size[k]; // update lt
			for (polysize_t i = 0; i < lt; ++i) {
				mpz_set(tmpZ, subres_work->coefs[k][i]);
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
				mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
				mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
				mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
					mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
				}

				isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
			}
		}

		free(Pptr);
		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			if ((*info_size) < 0 || bspecInfoArray == NULL) {
				for (int i = 0; i < infoIdx; ++i) {
					freeSpecA_spXY (bspecInfo[i]->An);
					freeSpecQ_spXY (bspecInfo[i]->Qn);
					freeSpecR_spXY (bspecInfo[i]->Rn);
				}
				free(bspecInfo);
			// } else if (*info_size < infoIdx) {
			} else if (isInfoNull) {
				if ((*bspecInfoArray) == NULL) {
					*bspecInfoArray = (specAQRArray_spXY_t *) malloc (sizeof (specAQRArray_spXY_t));
				}
				(*bspecInfoArray)->bspecArray = bspecInfo;
				(*bspecInfoArray)->size = infoIdx;
				(*info_size) = infoIdx;
			}
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			*ksubres = head;
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, tmpZ, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	// fprintf (stderr, "In hgcdBiModularFFTSubresultantChain_DBZP, all primes failed\n"); // TEST
	*ksubres = NULL;
	*chain_size = 0;
	return 0;
}