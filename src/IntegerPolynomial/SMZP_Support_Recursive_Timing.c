

#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "IntegerPolynomial/SMZP_Support_Test.h"
#include <time.h>
#include "Utils/Unix_Timer.h"

#if defined(WITH_BLAD)
#include "BLADInterface/bladinterface.h"
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

	mergeSortPolynomial_AAZ(ret);
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
	long long unsigned int* masks = getExpMaskArray(aa->nvar);
	long long unsigned int mvarMask = masks[0];
	free(masks);

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
	while (delta > -1 || kDeg > -1) {
		
		//get the leading term of dividend quotient-product difference
		delta = recProdHeapPeek_AA(prodHeap);

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

		degrees_t* masks = getExpMaskArray(nvar);
		int* sizes = getExpOffsetArray(nvar);

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

void freeAltArrsZ (AltArrsZ_t* AAs)
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
}

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
	fprintf (stderr, "SMZP Error: In LazardOptZ, Sd->nvar(=%d) != Sdm->nvar(=%d).\n",
		 Sd->nvar, Sdm->nvar);
	exit (1);
    }

    /* fprintf (stderr, "Hello Lazard"); */
    /* fprintf (stderr, "Sd, Sdm = \n"); */
    /* printAAZ (Sd); */
    /* printAAZ (Sdm); */
    /* fprintf (stderr, "\n\n"); */
    
    register int nvar = Sd->nvar;
    
    /* degree_t Sd_deg = (Sd->elems[0].degs & mvarMask) >> mvarDegOffset; */
    /* degree_t Sdm_deg = (Sdm->elems[0].degs & mvarMask) >> mvarDegOffset; */
    degree_t n = mainLeadingDegree_AAZ(Sd) - mainLeadingDegree_AAZ(Sdm) - 1;//((Sd->elems[0].degs & mvarMask) >> mvarDegOffset) -
	// ((Sdm->elems[0].degs & mvarMask) >> mvarDegOffset) - 1;
    if (n < 1){
	return deepCopyPolynomial_AAZ (Sdm);
    }
    
    AltArrZ_t* __restrict__ x = mainLeadingCoefficient_AAZ(Sdm);
    AltArrZ_t* __restrict__ y = s; //mainLeadingCoefficient_AAZ(Sd);	    
    
    if (x == NULL || x->size == 0){
	fprintf (stderr, "SMZP Error: In LazardOptZ, x is NULL!");
	exit (1);
    }
    if (y == NULL || y->size == 0){
	fprintf (stderr, "SMZP Error: In LazardOptZ, y is NULL!");
	exit (1);
    }
    
    AltArrZ_t* c = deepCopyPolynomial_AAZ (x);
    register degree_t a = 1 << Log2n (n);
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
	    if (c != NULL && c->size != 0)
		tmp = multiplyPolynomials_AAZ_inp (c, x, nvar);
	    else
		tmp = NULL;
	    exactDividePolynomials_AAZ (tmp, y, &tmpz, nvar);
	    freePolynomial_AAZ (tmp);
	    c = tmpz;
	    n = n - a;
	}
    }
    
    tmpz = NULL;
    AltArrZ_t* Se = deepCopyPolynomial_AAZ (Sdm);
    Se = multiplyPolynomials_AAZ_inp (Se, c, nvar);
    freePolynomial_AAZ (c);
    /* if (!orign){ */
    /* 	return Se; */
    /* } */
    
    exactDividePolynomials_AAZ (Se, y, &tmpz, nvar);
    freePolynomial_AAZ (Se);    
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
    
    return tmpz;    
}

AltArrZ_t* CFDucosOptZ (AltArrZ_t* A, AltArrZ_t* Sdm, AltArrZ_t* Se, AltArrZ_t* sd)
{

    register int nvar = A->nvar;
    
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

    AltArrZ_t* __restrict__ cdm = mainLeadingCoefficient_AAZ (Sdm);
    AltArrZ_t* __restrict__ se = mainLeadingCoefficient_AAZ (Se);

    AltArrZ_t* tmpz = NULL;
    AltArrZ_t* tmpR = NULL;
    AltArrZ_t* tmp = NULL;
    AltArrZ_t* res = NULL;
    AltArrZ_t* PIe = NULL;
    AltArrZ_t* D2 = NULL;
    AltArrZ_t* PIj = NULL;
    AltArrZ_t* D = NULL;	
    AltArrZ_t* H = NULL;
    
    if (cdm == NULL || cdm->size ==0){
	fprintf (stderr, "SMZP Error: In CFDucosOptZ, cdm = NULL.\n");
	exit (1);
    }
    if (se == NULL || se->size ==0){
	fprintf (stderr, "SMZP Error: In CFDucosOptZ, se = NULL.\n");
	exit (1);
    }
    if (d == e){
	// for non-defective cases:
	/* fprintf (stderr, "SMZP Error: In CFDucosOptZ, d == e.\n"); // TEST */

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
    
    // D2 = 
    PIe = mainCoefficientAtIdx_AAZ (A, e);    
    if (PIe != NULL && H != NULL && PIe->size != 0 && H->size != 0){
	D2 = multiplyPolynomials_AAZ (H, PIe, nvar);
	freePolynomial_AAZ (PIe);
    } 

    // PIe = mainCoefficientAtIdx_AAZ (H, e-1);  // it must be uncommented for using CFDucosOptZ_subPolynomials_AAZ_inp
    for (int j = e+1; j < d; ++j){
	//  H = XH - pi_e (XH)*C
	H = mainLShiftPolynomial_AAZ_inp (H, 1); // H = XH
	PIe = mainCoefficientAtIdx_AAZ (H, e);	// it must be commented for using CFDucosOptZ_subPolynomials_AAZ_inp
	if (PIe != NULL && PIe->size != 0){
	    tmp = multiplyPolynomials_AAZ_inp (PIe, Sdm, nvar);
	    tmpz = NULL;
	    tmpR = NULL;
	    exactDividePolynomials_AAZ (tmp, cdm, &tmpz, nvar);
	    freePolynomial_AAZ (tmp); //  freePolynomial_AAZ (PIe);
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
    
    // D =
    for (int j = 0; j < e; ++j){
	PIj = mainCoefficientAtIdx_AAZ (A, j);
	if (PIj != NULL && PIj->size != 0){
	    PIj = mainLShiftPolynomial_AAZ_inp (PIj, j);
	    PIj = multiplyPolynomials_AAZ_inp (PIj, se, nvar);
	}
	if (D == NULL || D->size == 0){
	    D = deepCopyPolynomial_AAZ (PIj);
	} else {
	    D = addPolynomials_AAZ_inp (D, PIj, nvar);
	}
	freePolynomial_AAZ (PIj);
    }
    
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
    
    if ((d-e+1)%2 != 0 && tmpz != NULL && tmpz->size != 0){
	negatePolynomial_AAZ (tmpz);
    }
    
    return tmpz;    
}

void DucosSubresultantChainZ_rev (AltArrZ_t* P, AltArrZ_t* Q, AltArrsZ_t** SC, int* len, int type)
{

	/* char* ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */
	/* fprintf (stderr, "[SMZP-Ali] input polynomial in DucosSubresultant  "); */
	/* fprintf (stderr, "[SMZP-Ali] P :=  "); */
	/* printPoly_AAZ_unpk(stderr, P, ch, P->nvar); */
	/* fprintf (stderr, "\n"); */
	/* fprintf (stderr, "[SMZP-Ali] Q :=  "); */
	/* printPoly_AAZ_unpk(stderr, Q, ch, Q->nvar); */
	/* fprintf (stderr, "\n"); */

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
    if (mainLeadingDegree_AAZ (P) == 0 || mainLeadingDegree_AAZ (Q) == 0){
        fprintf (stderr, "SMZP: error, Input polynomials to subresultantChain must have positive degree.\n");
	exit(1);
    }
    
    int lazy = 0; // lazy option in pseudoDivide algorithm
    register int nvar = P->nvar;
    
    int tmpE = 0;
    AltArrZ_t* tmpz = NULL;
    AltArrZ_t* tmpH = NULL;
    AltArrZ_t* tmpB = NULL;
    AltArrZ_t* s = NULL;
    AltArrZ_t* p = NULL;
    AltArrZ_t* q = NULL;
    AltArrZ_t* A = NULL;
    AltArrZ_t* B = NULL;
    AltArrZ_t* C = NULL;

    AltArrsZ_t* SS = NULL;
    AltArrsZ_t* tail = NULL;
    AltArrZ_t* tailm = NULL;
    degree_t delta = 0;
    int size = 0;
    
    degree_t degP = mainLeadingDegree_AAZ(P);//(P->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t degQ = mainLeadingDegree_AAZ(Q);//(Q->elems[0].degs & mvarMask) >> mvarDegOffset;
    
    if (degP >= degQ){
	p = P;
	q = Q;
    } else {
	p = Q;
	q = P;
    }

    /* char sym[11] = {'x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'x10'} */
    // TEST
     // fprintf (stderr, "nvar = %d\n", p->nvar); 
     // fprintf (stderr, "In DucosSubResult, P = "); 
     // printAAZ (p); 
     // fprintf (stderr, "In DucosSubResult, Q = "); 
     // printAAZ (q); 
        
    degP = mainLeadingDegree_AAZ(p);//(p->elems[0].degs & mvarMask) >> mvarDegOffset;
    degQ = mainLeadingDegree_AAZ(q);//(q->elems[0].degs & mvarMask) >> mvarDegOffset;
    s = mainLeadingCoefficient_AAZ (q);
    s = exponentiatePoly_AAZ (s, degP - degQ, nvar);
    A = deepCopyPolynomial_AAZ (q);
    
    // B =
    AltArrZ_t* negQ = deepCopyPolynomial_AAZ (q);
    negatePolynomial_AAZ (negQ);

	/* fprintf (stderr, "q := "); */
	/* printAAZ (q); */

	/* fprintf (stderr, "negQ := "); */
	/* printAAZ (negQ); */

    // deepCopyPolynomial_AAZ (p), deepCopyPolynomial_AAZ (negQ)
    pesudoDivideAtIdx_AAZ (0, p, negQ, &tmpz, &B, &tmpE, &tmpH, nvar, lazy);

	/* fprintf (stderr, "p := "); */
	/* printAAZ (p); */
	
	/* fprintf (stderr, "negQ := "); */
	/* printAAZ (negQ); */

	/* fprintf (stderr, "p/negQ := "); */
	/* printAAZ (tmpz); */
	
	/* fprintf (stderr, "rem := "); */
	/* printAAZ (B); */

	/* fprintf (stderr, "h := "); */
	/* printAAZ (tmpH); */

    freePolynomial_AAZ (tmpz);
    freePolynomial_AAZ (tmpH);
    freePolynomial_AAZ (negQ);

	
    AltArrsZ_t* newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
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
    
    register degree_t d = 0;
    register degree_t e = 0;
    
    while (1){

    	if (B == NULL || B->size == 0){
    		break;
    	}

    	d = mainLeadingDegree_AAZ(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
    	e = mainLeadingDegree_AAZ(B);//(B->elems[0].degs & mvarMask) >> mvarDegOffset;

    	newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
    	newPoly->poly = deepCopyPolynomial_AAZ (B);
    	newPoly->next = NULL;

    	tailm = tail->poly;
    	tail->next = newPoly;
    	tail = newPoly;
    	size++;

    	delta = d - e;
    	if (delta > 1){

			C = LazardOptZ (tailm, tail->poly, s); // TODO: add s to the inputs 
			newPoly = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			newPoly->poly = deepCopyPolynomial_AAZ (C);
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

		/* fprintf (stderr, "ChainSize := %d\n", size); */
		if (type == 0){
			tmpB = DucosOptZ (A, B, C, s);
		} else {
			tmpB = CFDucosOptZ (A, B, C, s);
		}
		freePolynomial_AAZ (B);
		freePolynomial_AAZ (A);
		B = tmpB;
		A = deepCopyPolynomial_AAZ (C);
		freePolynomial_AAZ (s);
		s = mainLeadingCoefficient_AAZ (A);	
		freePolynomial_AAZ (C);
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
    DucosSubresultantChainZ_rev (P, Q, &invSC, &size, 1);

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

AltArrZ_t* DucosResultantZ (AltArrZ_t* P, AltArrZ_t* Q)
{
    AltArrsZ_t* SC;
    int size = 0;
    DucosSubresultantChainZ (P, Q, &SC, &size);
    if (size){
		return deepCopyPolynomial_AAZ (SC->poly);
    }
    return NULL;
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
    
    AltArrsZ_t* SC;
    int size = 0;
    DucosSubresultantChainZ (P, Q, &SC, &size);
    AltArrsZ_t* cur = SC;

    if (cur->poly == NULL || cur->poly->size == 0){
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
    		return makeConstPolynomial_AAZ (1, nvar, one);
    	} else {
    		return deepCopyPolynomial_AAZ (cur->poly);
    	}
    }
    return deepCopyPolynomial_AAZ (cur->next->poly);
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

     // fprintf (stderr, "\n\n In lastNonZeroSub, P = "); 
     // printAAZ (P); 
     // fprintf (stderr, "\n In lastNonZeroSub, Q = "); 
     // printAAZ (Q); 
     // fprintf (stderr, "\n computing DucosSubresultantChainZ... \n"); 
    
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

AltArrZ_t* gcd_AAZ (AltArrZ_t* P, AltArrZ_t* Q)
{

	/* char* ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */

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

#if defined(WITH_BLAD)
    if (P->size > 8 && Q->size > 8) {
    	char* c_names[P->nvar];
    	for (int i = 0; i < P->nvar; ++i) {
    		c_names[i] = (char*) malloc(sizeof(char)*8);
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

	timer_id id = start_timer();

    AltArrZ_t* lnzch = lastNonZeroChain_AAZ (ppP, ppQ);
    timer_time elapsed1 = elapsed_time(&id);
    double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
    fprintf(stderr, "lnzch time: %f\n", time);
    
    if (isShrinked && lnzch != NULL && lnzch->size != 0){
	expandNumVarsLeft_AAZ (lnzch, nvar);
    }

    if (leadingVariable_AAZ (lnzch) != mvarP){
	mpz_t one;
	mpz_init (one);
	mpz_set_si (one, 1l);
	gcd = makeConstPolynomial_AAZ (1, nvar, one);
    } else {
	 id = start_timer();

  	gcd = mainPrimitivePart_AAZ (lnzch, mvarP);
    elapsed1 = elapsed_time(&id);
     time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
    fprintf(stderr, "prim part of gcd time: %f\n", time);
    }

    gcd = multiplyPolynomials_AAZ_inp (gcd, cont, nvar);

    freePolynomial_AAZ (contP);
    freePolynomial_AAZ (contQ);
    freePolynomial_AAZ (ppP);
    freePolynomial_AAZ (ppQ);
    freePolynomial_AAZ (cont);
    
    return gcd;
}

AltArrZ_t* mainPrimitiveFactorization_AAZ (AltArrZ_t* P, AltArrZ_t** cont)
{
    int mvarP = leadingVariable_AAZ (P);
    AltArrZ_t* sP;
    int isShrinked = 0;
    mpz_t one;
    mpz_init (one);
    mpz_set_si (one, 1l);

    if (mvarP < 0){
	if (mvarP == -1){
	    /* if (*cont == NULL || (*cont)->size == 0){ */
	    mpz_t zz;
	    mpz_init (zz);	    
	    mpz_set(zz, P->elems[0].coef);
	    *cont =  makeConstPolynomial_AAZ (1, P->nvar, zz);
	    /* } */
	}
	
	return makeConstPolynomial_AAZ (1, P->nvar, one);
    }
    
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
    AltArrZ_t* coefGCD = deepCopyPolynomial_AAZ (convertFromAAZElemToAAZ
						 (elems[0].coef,  elems[0].coefSize,
						  sP->nvar, recP->unpacked));
    
    double time = 0;
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
		timer_id id = start_timer();
		tmp = gcd_AAZ (coefGCD, tmpi);
        timer_time elapsed1 = elapsed_time(&id);
        time += (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
		freePolynomial_AAZ (coefGCD);
		freePolynomial_AAZ (tmpi);
		coefGCD = tmp;
		++idx;
    }
    // fprintf(stderr, "primfact gcd time: %f\n", time);
    
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
	
	 char* ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST 

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

        timer_id id = start_timer();
		swappp = mainPrimitiveFactorization_AAZ (primpart, &cont);
        timer_time elapsed1 = elapsed_time(&id);
        double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
        fprintf(stderr, "primfact time: %f\n", time);

		freePolynomial_AAZ (primpart);
		
		// diff = d(primpart)/d(x_i) 
		id = start_timer();

		diff = derivative_AAZ (swappp, 0, 1); 
        elapsed1 = elapsed_time(&id);
        time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
        fprintf(stderr, "diff time: %f\n", time);

		
		 // fprintf (stderr, "G := "); 
		 // printPoly_AAZ_unpk(stderr, swappp, ch, nnvar); 
		 // fprintf (stderr, "\n"); 
		
		
		 // fprintf (stderr, "F := "); 
		 // printPoly_AAZ_unpk(stderr, diff, ch, nnvar); 
		 // fprintf (stderr, "\n"); 
		
		// g = gcd(swappp, D(swappp)) 
		id = start_timer();
		g = gcd_AAZ (swappp, diff);
        elapsed1 = elapsed_time(&id);
        time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
        fprintf(stderr, "gcd time: %f\n", time);

		
		freePolynomial_AAZ (diff);

		// next = swappp / g 
		id = start_timer();
		exactDividePolynomials_AAZ (swappp, g, &next, nnvar);
        elapsed1 = elapsed_time(&id);
        time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
        fprintf(stderr, "divide time: %f\n", time);
		
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
