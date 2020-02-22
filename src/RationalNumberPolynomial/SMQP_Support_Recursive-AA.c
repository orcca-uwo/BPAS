
#include "RationalNumberPolynomial/SMQP_Support_Recursive-AA.h"
#include "RationalNumberPolynomial/SMQP_Support_Test-AA.h"

/** 
 * Build a random recursive polynomial. 
 * Nterms is the number of terms wrt the main variable. 
 * coefNTerms is the number of terms in the coefficients (when viewed recursively, and nvar > 1);
 * coefBound is the upper limit on the number of bits of rational numbers in the polynomial.
 * sparsity is the difference between degrees (wrt to main variable) of successive terms - 1. Therefore 2 is dense. 
 * sparsity is also passed when constructing coefficient polynomials.
 */
RecArr_t* buildRandomRecArrPoly(int nvar, int nterms, int coefNterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
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
		AltArr_t* aa = deepCopyPolynomial_AAFromNode(head, nvar);
		freePolynomial(head);
		return convertToRecursiveArray(aa);
	} else {

		Node* head = NULL, *tail = NULL;
		degrees_t mvarDeg;

		for (int i = 0; i < nterms; ++i) {
			Node* coefNode = buildRandomPoly(nvar-1, coefNterms, coefBound, sparsity, includeNeg);
			mvarDeg = (nterms - 1 - i) * (sparsity - 1);
			if (head == NULL) {
				head = coefNode;
			}
			if (tail != NULL) {
				tail->next = coefNode;
			}

			Node* curNode = coefNode;
			while (curNode != NULL) {
				//If we want integers: 
				// mpz_set_si(mpq_denref(curNode->coef), 1l);
				//expand node by 1 var.
				curNode->degs += (mvarDeg << mvarDegOffset);
				if (curNode->next == NULL) {
					tail = curNode;
				}
				curNode = curNode->next;
			}
		}

		//Node* head is now the poly we want to make recursive. Let's turn it into an AA.
		AltArr_t* aa = deepCopyPolynomial_AAFromNode(head, nvar);
		freePolynomial(head);
		return convertToRecursiveArray(aa);
	}	
}

RecArr_t* convertToRecursiveArray(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (aa->unpacked) {
		return convertToRecursiveArray_unpk(aa);
	}

	int mvarDegOffset = getMVarExpOffset(aa->nvar);
	degrees_t mvarMask = getMVarExpMask(aa->nvar);

	AAElem_t* elems = aa->elems;
	degree_t mvarDeg =  (elems[0].degs & mvarMask) >> mvarDegOffset;
	RecArr_t* poly = (RecArr_t*) malloc(sizeof(RecArr_t));
	poly->alloc = mvarDeg;
	RecArrElem_t* recElems = (RecArrElem_t*) malloc(sizeof(RecArrElem_t)*(mvarDeg+1)); 
	poly->elems = recElems;

	int curIdx = 0, lastSize = 0;
	degree_t curDeg;
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		curDeg = (elems[i].degs & mvarMask) >> mvarDegOffset;
		
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

RecArr_t* convertToRecursiveArrayAtIdx (AltArr_t* aa, int idx) 
{
  if (aa == NULL){
      return NULL;
  }
  if (idx < 0 || idx >= aa->nvar){
    fprintf(stderr, "SMQP Error: idx is out of range!\n");
    exit (EXIT_FAILURE);
  }

  AltArr_t* cPoly = deepCopyPolynomial_AA(aa);
  if (idx == 0){
    return convertToRecursiveArray(cPoly);
  }

  int varMap[cPoly->nvar];
  for (int i = 0; i < cPoly->nvar; ++i){
    varMap[i] = i;
  }
  varMap[0] = idx;
  varMap[idx] = 0;

  reorderVars_AA (cPoly, varMap, cPoly->nvar);

  return convertToRecursiveArray(cPoly);
}

AltArr_t* convertFromRecursiveArray(RecArr_t* poly, int nvar) {
	if (poly == NULL) {
		return NULL;
	}

	if (poly->unpacked) {
		return convertFromRecursiveArray_unpk(poly, nvar);
	}

	AltArr_t* aa = NULL;
	if (poly->origAA != NULL) {
		aa = poly->origAA;
	} else {
		aa = (AltArr_t*) malloc(sizeof(AltArr_t));
	}


	int mvarDegOffset = getMVarExpOffset(nvar);

	RecArrElem_t* recElems = poly->elems;
	AAElem_t* elems = poly->elems->coef;
	aa->elems = elems;
	register int curSize = 0;
	for (int i = 0; i < poly->size; ++i) {
		degrees_t curDeg = recElems[i].exp;
		for (int j = 0; j < recElems[i].coefSize; ++j) {
			elems[curSize + j].degs += (curDeg << mvarDegOffset);
		}

		curSize += recElems[i].coefSize;
	}

	freeRecArray(poly);

	aa->size = curSize;
	aa->alloc = curSize;
	aa->nvar = nvar;
	aa->unpacked = 0;


	return aa;
}

AltArr_t* convertFromRecursiveArrayAtIdx (RecArr_t* recPoly, int idx, int nvar)
{
  if (recPoly == NULL){
      return NULL;
  } 
  if (idx < 0 || idx >= nvar){
      fprintf(stderr, "SMQP Error: idx is out of range!\n");
      exit (EXIT_FAILURE);
  }

  AltArr_t* cPoly = deepCopyPolynomial_AA (convertFromRecursiveArray (recPoly, nvar));
  if (idx == 0){
      return cPoly;
  }

  int varMap[nvar];
  for (int i = 0; i < nvar; ++i){
    varMap[i] = i;
  }
  varMap[0] = idx;
  varMap[idx] = 0;

  reorderVars_AA (cPoly, varMap, nvar);
  return cPoly;
} 

RecArr_t* deepCopyRecArrayPolynomial(RecArr_t* poly, int nvar) {
	if (poly == NULL || poly->size == 0) {
		return NULL;
	}

	// Actually, we rely on deepCopyPolynomial_AA to handle
	// the deep copy of the unpacked exponents and the like 
	// if (poly->unpacked) {
	// 	return deepCopyRecArrayPolynomial_unpk(poly, nvar);
	// }

	int maxCoefSize = 0;
	RecArrElem_t* recElems = poly->elems;
	for (int i = 0; i < poly->size; ++i) {
		maxCoefSize += recElems[i].coefSize;
	}

	AltArr_t localAA;
	localAA.unpacked = poly->unpacked;
	localAA.elems = poly->elems->coef;
	localAA.size = localAA.alloc = maxCoefSize;
	localAA.nvar = nvar;
	AltArr_t* copyAA = deepCopyPolynomial_AA(&localAA);

	AAElem_t* elems = copyAA->elems;
	free(copyAA); //free the container but not underlying elems.

	// AAElem_t* origElems = poly->elems->coef;
	// AAElem_t* elems = (AAElem_t*) malloc(sizeof(AAElem_t)*maxCoefSize);
	// for (int i = 0; i < maxCoefSize; ++i) {
	// 	mpq_init(elems[i].coef);
	// 	mpq_set(elems[i].coef, origElems[i].coef);
	// 	elems[i].degs = origElems[i].degs;
	// }

	RecArrElem_t* copyElems = (RecArrElem_t*) malloc(sizeof(RecArrElem_t)*poly->size);
	memcpy(copyElems, poly->elems, sizeof(RecArrElem_t)*poly->size);
	copyElems->coef = elems;
	maxCoefSize = poly->elems->coefSize;
	for (int i = 1; i < poly->size; ++i) {
		copyElems[i].coef = elems + maxCoefSize;
		maxCoefSize += copyElems[i].coefSize;
	}

	RecArr_t* retPoly = (RecArr_t*) malloc(sizeof(RecArr_t));
	retPoly->size = retPoly->alloc = poly->size;
	retPoly->elems = copyElems;
	retPoly->unpacked = poly->unpacked;

	retPoly->origAA = NULL;

	return retPoly;

}



/***********
 * Pesudo Division
 ***********/

void recProdHeapInsert_AA(RecProdHeap_AA* h, RecProdHeapChain_AA* chain, register degree_t insertExp) {
	
	int s = h->heapSize;
	// int N = h->maxHeapSize;
	RecProdHeapElem_AA* elems = h->elements;
	if (s == 0) {
		elems[0].exp = insertExp;
		elems[0].chain = chain;
		h->heapSize = 1;
		return;
	}

	//first check if we can chain off the root
	if (insertExp == elems->exp) {
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
	//
	//     x     --->    e
	//    / \           / \ 
	//   x   o         x   o
	//                /
 	//               x

	register int i = (s-1)/2; //i is parent 
	register int j = s;       //j is current insertion point
	register long long unsigned int path = 1;
	// degree_t cmpExp;
	while (j > 0) {
		if (elems[i].exp == insertExp) {
			chain->next = elems[i].chain;
			elems[i].chain = chain;
			return;
		} else if (elems[i].exp < insertExp) {
			path <<= 1;
			if (!(j & 1)) {
				//set trailing bit to 1 to distinguish left/right path
				path += 1;
			}
			j = i;
			i = (i - 1) >> 1;
		} else {
			break; //found the insert point
		}
	}

	// //if we made it this far, then we are adding a new element to the heap. 
	// //resize if necessary
	// if (s >= N) {
	// 	N = (s == 0) ? 1 : s * 2;
	// 	h->elements = (RecProdHeapElem_t**) realloc(h->elements, N*sizeof(RecProdHeapElem_t*));
	// 	elems = h->elements;
	// 	h->maxHeapSize = N;
	// }

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	RecProdHeapElem_AA temp;
	RecProdHeapElem_AA elem = {insertExp, chain};

	while (j <= s) {
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = (j << 1) + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

RecProdHeapChain_AA* recProdHeapRemoveMax_AA(RecProdHeap_AA* h) {
	
	RecProdHeapElem_AA* elems = h->elements;
	RecProdHeapChain_AA* maxElem = elems[0].chain;
	register int i = 0;
	register int j = 1;
	register int s = --(h->heapSize);
	
	//promote largest children
	while (j < s) {
		if (j+1 < s && elems[j].exp < elems[j+1].exp) {
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (j << 1) + 1;
	}
	//now place last element into i and swim up to make tree complete 
	j = (i-1) >> 1;
	while(i > 0) {
		if (elems[s].exp < elems[j].exp) {
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j-1) >> 1;
	}
	elems[i] = elems[s]; 

	return maxElem;	
}

AltArr_t* recProdHeapGetNextCoef_AA(RecProdHeap_AA* h, const RecArrElem_t* __restrict__ aElems, const RecArrElem_t* __restrict__ bElems, const int* aUnpacked, int bUnpacked) {
	RecProdHeapElem_AA* elems = h->elements;
	int nvar = h->nvar;
	register int lastB = h->lastB;

	AltArr_t* ret = NULL;
	register degree_t maxExp = elems->exp;
	register degree_t nextExp = elems->exp;

	RecProdHeapChain_AA* insertChain = NULL;
	RecProdHeapChain_AA* maxElem, * nextMaxElem;
	
	AltArr_t aCoef, bCoef;
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
				ret = multiplyPolynomials_AA(&aCoef, &bCoef, nvar);
			} else {
				AltArr_t* prod = multiplyPolynomials_AA(&aCoef, &bCoef, nvar);
				ret = addPolynomials_AA_inp(ret, prod, nvar);
				freePolynomial_AA(prod);
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
void recArrayMultiplyByCoef_inp(RecArrElem_t* aElems, int aSize, AltArr_t* coef, int nvar, int* aUnpacked) {
	AltArr_t* tempa = malloc(sizeof(AltArr_t));
	tempa->nvar = nvar;
	for (int i = 0; i < aSize; ++i) {
		tempa->unpacked = aUnpacked[i];
		tempa->elems = aElems[i].coef;
		tempa->alloc = tempa->size = aElems[i].coefSize;
		tempa = multiplyPolynomials_AA_inp(tempa, coef, nvar);

		aUnpacked[i] = tempa->unpacked;
		aElems[i].coef = tempa->elems;
		aElems[i].coefSize = tempa->size;
	}

	free(tempa);
}


static void printAA(AltArr_t* aa) {
    if (aa == NULL || aa->size == 0 ){
	printf("\n");
	return;
    }
    for (int i = 0; i < AA_SIZE(aa); ++i) {
	gmp_printf( "%Qd*%llx + ", aa->elems[i].coef, aa->elems[i].degs);
    }
    printf("\n");
}


void pesudoDivideOneTerm_RecArray(RecArr_t* c, RecArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy) {

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
		RecArr_t* recR = deepCopyRecArrayPolynomial(c, nvar); 
		*res_r = convertFromRecursiveArray(recR, nvar);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (hPow != NULL) {
			*hPow = makePolynomial_AA(1, nvar);
			mpq_init((*hPow)->elems->coef);
			mpq_set_si((*hPow)->elems->coef, 1l, 1l);
			(*hPow)->elems->degs = 0;
			(*hPow)->size = 1;
		}
		return;
	}

	RecArrElem_t* __restrict__ k = c->elems;
	RecArrElem_t* __restrict__ lenK = k + c->size;

	register int mvarOffset = getMVarExpOffset(nvar);
	register int maxASize = k->coefSize;
	register int maxRSize = k->coefSize;
	register int ai = 0;
	register int rj = 0;

	AltArr_t* a = makePolynomial_AA(maxASize, nvar);
	AltArr_t* r = makePolynomial_AA(maxRSize, nvar);

	register degree_t bDeg = b->elems->exp;
	register degree_t eps = k->exp - bDeg;
	AltArr_t kCoef;
	kCoef.alloc = kCoef.size = k->coefSize;
	kCoef.elems = k->coef;
	kCoef.nvar = nvar;
	kCoef.unpacked = c->unpacked;

	int i = 1;
	AltArr_t h;
	h.alloc = h.size = b->elems->coefSize;
	h.elems = b->elems->coef;
	h.nvar = nvar;
	h.unpacked = b->unpacked;

	AltArr_t* multerm = deepCopyPolynomial_AA(&kCoef);
	AltArr_t* hI = deepCopyPolynomial_AA(&h);

	//do first division manually 
	//update degs in-place by degree of the mainvar, eps.
	memcpy(a->elems + ai, multerm->elems, sizeof(AAElem_t)*multerm->size);
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
// 		multerm = deepCopyPolynomial_AA(&kCoef);
// #else
// 		multerm = multiplyPolynomials_AA(&kCoef, hI, nvar);
// #endif
		multerm = deepCopyPolynomial_AA(&kCoef);
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
// 				a = multiplyPolynomials_AA_inp(a, &h, nvar);	
	
// 				//update h counter
// 				//this is only done when we actually add a new element to the quotient.
// 				++i;
// 				hI = multiplyPolynomials_AA_inp(hI, &h, nvar); //in-place wrt first arg.
// #endif

			eps -= bDeg;

			//resize if necessary for insertion to a
			if (ai + multerm->size > maxASize) {
				maxASize += (multerm->size << 1);
				resizePolynomial_AA(a, maxASize);
			}

			memcpy(a->elems + ai, multerm->elems, sizeof(AAElem_t)*multerm->size);
			if(multerm->unpacked || a->unpacked) {
				if (!a->unpacked) {
					unpackExponentVectors_AA_inp(a);
				}
				if (!multerm->unpacked) {
					unpackExponentVectors_AA_inp(multerm);
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
			a->size = ai;
		} else {
			//accumulate remainder.
			if (rj + multerm->size > maxRSize) {
				r->size = rj;
				maxRSize += (multerm->size << 1);
				resizePolynomial_AA(r, maxRSize);
			}

			//copy multerm coef data directly into r
			memcpy(r->elems + rj, multerm->elems, sizeof(AAElem_t)*multerm->size);
			if (multerm->unpacked || r->unpacked) {
				if (!r->unpacked) {
					unpackExponentVectors_AA_inp(r);
				}
				if (!multerm->unpacked) {
					unpackExponentVectors_AA_inp(multerm);
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
		freePolynomial_AA(a);
		a = NULL;
	}

	if (rj > 0) {
		r->size = rj;
		r->alloc = maxRSize;
	} else {
		freePolynomial_AA(r);
		r = NULL;
	}

	if (!lazy) {
		int d = c->elems->exp - b->elems->exp + 1 - i;
		i += d;
		for (int j = 0; j < d; ++j) {
			if (ai > 0) {
				a = multiplyPolynomials_AA_inp(a, &h, nvar);
			}
			if (rj > 0) {
				r = multiplyPolynomials_AA_inp(r, &h, nvar);
			}
			hI = multiplyPolynomials_AA_inp(hI, &h, nvar);
		}
// #if PDIVIDE_DIVISBLE_CHECK
// 		//r needs one additional one since a and h have one extra from beginning of loop
		if (rj > 0) {
			r = multiplyPolynomials_AA_inp(r, &h, nvar);	
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

void pesudoDivide_RecArray(RecArr_t* c, RecArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy) {

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
		RecArr_t* recR = deepCopyRecArrayPolynomial(c, nvar); 
		*res_r = convertFromRecursiveArray(recR, nvar);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (hPow != NULL) {
			*hPow = makePolynomial_AA(1, nvar);
			mpq_init((*hPow)->elems->coef);
			mpq_set_si((*hPow)->elems->coef, 1l, 1l);
			(*hPow)->elems->degs = 0;
			(*hPow)->size = 1;
		}
		return;
	}

	if (b->size == 1) {
		pesudoDivideOneTerm_RecArray(c, b, res_a, res_r, e, hPow, nvar, lazy);
		return;
	}

	RecArrElem_t* __restrict__ k = c->elems;
	RecArrElem_t* __restrict__ b2Elem = b->elems + 1;
	RecArrElem_t* __restrict__ lenK = k + c->size;

	register int mvarOffset = getMVarExpOffset(nvar);
	register int maxASize = c->size < 2 ? 2 : c->size;
	register int maxRSize = maxASize;
	register int ai = 0;
	register int rj = 0;

	int* aUnpacked = calloc(maxASize, sizeof(int));

	AAElem_t* r = (AAElem_t*) malloc(sizeof(AAElem_t)*maxRSize);
	degree_t* rDegs = (degree_t*) malloc(sizeof(degree_t)*maxRSize*nvar);
	RecArrElem_t* aElems = (RecArrElem_t*) malloc(sizeof(RecArrElem_t)*maxASize); 
	RecArrElem_t* __restrict__ curA = aElems;

	//manually do first div as we know it comes from first term of c;
	int i = 1;
	AltArr_t h;
	h.alloc = h.size = b->elems->coefSize;
	h.elems = b->elems->coef;
	h.nvar = nvar;
	h.unpacked = b->unpacked;

	degree_t bDeg = b->elems->exp, eps = k->exp - bDeg;
	AltArr_t kCoef;
	kCoef.alloc = kCoef.size = k->coefSize;
	kCoef.elems = k->coef;
	kCoef.nvar = nvar;
	kCoef.unpacked = c->unpacked;

	AltArr_t* multerm = deepCopyPolynomial_AA(&kCoef);

	curA->coef = multerm->elems;
	curA->coefSize = multerm->size;
	curA->exp = eps;
	aUnpacked[ai] = multerm->unpacked;
	++curA;
	++ai;
	free(multerm);
	++k;

	AltArr_t* multerm2 = NULL;
	//use long long so it can fit an entire degree_t but can also be compared to -1.
	degree_t delta = 1;
	degree_t kDeg = (k == lenK) ? -1 : k->exp;

	
	AltArr_t* hI = deepCopyPolynomial_AA(&h);

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
			multerm = recProdHeapGetNextCoef_AA(prodHeap, aElems, b->elems, aUnpacked, b->unpacked);			
			if (multerm == NULL || multerm->size == 0) {
				freePolynomial_AA(multerm);
				//the elements in the product with degree delta
				// ended up canceling out coefs
				delta = recProdHeapPeek_AA(prodHeap);
				continue; 
			}
			negatePolynomial_AA(multerm);
		} else if (delta < kDeg) {
			eps = kDeg;
			kCoef.alloc = kCoef.size = k->coefSize;
			kCoef.elems = k->coef;

			//then ctilde is from dividend
			multerm = multiplyPolynomials_AA(&kCoef, hI, nvar);
			++k;
			kDeg = (k == lenK) ? -1 : k->exp;
		} else {
			eps = delta;
			//combine both
			multerm2 = recProdHeapGetNextCoef_AA(prodHeap, aElems, b->elems, aUnpacked, b->unpacked);			
			if (multerm2 == NULL || multerm2->size == 0) {
				freePolynomial_AA(multerm2);

				delta = recProdHeapPeek_AA(prodHeap);
				continue;
			}			
			kCoef.alloc = kCoef.size = k->coefSize;
			kCoef.elems = k->coef;

		  	multerm = multiplyPolynomials_AA(&kCoef, hI, nvar);
			multerm = subPolynomials_AA_inp(multerm, multerm2, nvar);
			freePolynomial_AA(multerm2);
			
			++k;
			kDeg = (k == lenK) ? -1 : k->exp;
			
			if (multerm == NULL || multerm->size == 0) {
				//if sub resulted in a zero then we we must get a new multerm
				freePolynomial_AA(multerm);

				delta = recProdHeapPeek_AA(prodHeap);
				continue;
			}
		}

		//multerm is now the leading coef (eps the degree) of the current difference
		if (eps >= bDeg) {

#if PDIVIDE_DIVISBLE_CHECK
			AltArr_t* multermQuo = NULL;
			if (divideTest_AA(multerm, &h, &multermQuo, nvar)) {
				freePolynomial_AA(multerm);
				multerm = multermQuo;			
			} else {
#endif
				recArrayMultiplyByCoef_inp(aElems, ai, &h, nvar, aUnpacked);
				
				//update h counter
				//this is only done when we actually add a new element to the quotient.
				++i;
				hI = multiplyPolynomials_AA_inp(hI, &h, nvar); //in-place wrt first arg.

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
				aElems = (RecArrElem_t*) realloc(aElems, sizeof(RecArrElem_t)*maxASize);
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
				r = (AAElem_t*) realloc(r, sizeof(AAElem_t)*maxRSize);
				degree_t* oldDegs = rDegs;
				rDegs = (degree_t*) realloc(rDegs, sizeof(degree_t)*maxRSize*nvar);
				if (oldDegs != rDegs) {
					for (int idx = 0; idx < rj; ++idx) {
						r[idx].degs = (degrees_t) (rDegs + idx*nvar);
					}
				}
			}
			//copy multerm coef data directly into r
			memcpy(r + rj, multerm->elems, sizeof(AAElem_t)*multerm->size);
			unpackExponentVectors_AA_inp(multerm);
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
	AAElem_t* aBlock = (AAElem_t*) malloc(sizeof(AAElem_t)*aSize);

	if (allPacked) {
		aSize = 0;
		//first loop over each recursive element
		for (int idx = 0; idx < ai; ++idx) {
			memcpy(aBlock+aSize, aElems[idx].coef, sizeof(AAElem_t)*(aElems[idx].coefSize));
		
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
			memcpy(aBlock+aSize, aElems[idx].coef, sizeof(AAElem_t)*aElems[idx].coefSize);
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

	AltArr_t* aa = NULL;
	if (aSize > 0) {
		aa = (AltArr_t*) malloc(sizeof(AltArr_t));
		aa->elems = aBlock;
		aa->size = aSize;
		aa->alloc = aSize;
		aa->nvar = nvar;
		aa->unpacked = !allPacked;
	}
	AltArr_t* ra = NULL;
	if (rj > 0) {
		ra = (AltArr_t*) malloc(sizeof(AltArr_t));
		ra->elems = r;
		ra->size = rj;
		ra->alloc = maxRSize;
		ra->nvar = nvar;
		ra->unpacked = 1;
		tryPackExponentVectors_AA_inp(ra);
		//TODO
	}

	if (!lazy) {
		int d = c->elems->exp - b->elems->exp + 1 - i;
		i += d;
		for (int j = 0; j < d; ++j) {
			aa = aa == NULL ? NULL : multiplyPolynomials_AA_inp(aa, &h, nvar);
			ra = ra == NULL ? NULL : multiplyPolynomials_AA_inp(ra, &h, nvar);
			hI = multiplyPolynomials_AA_inp(hI, &h, nvar);
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

void pesudoDivideAtIdx_AA(int idx, AltArr_t* c, AltArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy) 
{
    if (idx >= nvar || idx < 0){
    	printAA (c);
		printAA (b);
		fprintf(stderr, "SMQP Error: idx(=%d) is out of range while nvar is %d!", idx, nvar);
		exit (EXIT_FAILURE);
    }

    // convert c and b to RecArr_t at idx:
    RecArr_t* tmpB = convertToRecursiveArrayAtIdx (b, idx);
    RecArr_t* tmpC = convertToRecursiveArrayAtIdx (c, idx);

    // compute pesudoDivide_RecArr(.):
    AltArr_t* q = NULL;
    AltArr_t* r = NULL;
    AltArr_t* h = NULL;
    int ee = 0;

    pesudoDivide_RecArray (tmpC, tmpB, &q, &r, &ee, &h, nvar, lazy);
    if (idx > 0){
    	// freeRecArray(tmpC);
    	// freeRecArray(tmpB);
		freeRecArrayAndCoef (tmpC);
		freeRecArrayAndCoef (tmpB);
    }

    // convert q, r, and hPow at idx and set the outputs:
    if (idx > 0){
		*res_a = swappingExponents_AA (q, 0, idx);
		*res_r = swappingExponents_AA (r, 0, idx);
		*hPow = swappingExponents_AA (h, 0, idx);
		*e = ee;
		freePolynomial_AA (q);
		freePolynomial_AA (r);
		freePolynomial_AA (h);
    } else {
		*res_a = q;
		*res_r = r;
		*hPow = h;
		*e = ee;
    }
}

void naiveMultiDivisorPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet) 
{
    if (nSet  < 0 || nvar < nSet){
		fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
		exit(EXIT_FAILURE);
    }
    
    mpq_t one;
    mpq_init (one);
    mpq_set_si (one, 1l, 1ul);
	
    if (nSet == 0){
		*rem = deepCopyPolynomial_AA(c);
		if (*hPow != NULL){
			freePolynomial_AA (*hPow);
		}
		*hPow = makeConstPolynomial_AA(1, nvar, one);
		quoSet[0] = NULL;
		return;
    }
	
    int e;
    int mvar;
    AltArr_t* h;
    AltArr_t* tmpQ;
    AltArr_t* tmpR;
    AltArr_t* r = deepCopyPolynomial_AA (c);
    AltArr_t* totalH = NULL;
	
    /* for (int i = 0; i < nSet; ++i)	 */
    for (int i = nSet-1; i >= 0; --i) {
	    mvar = leadingVariable_AA (B[i]);
	    if ((mvar > -1) && (r != NULL && r->size != 0)) {    
		    e = 1;
		    h = NULL;
		    tmpQ = NULL;
		    tmpR = NULL;
		    
		    pesudoDivideAtIdx_AA (mvar, r, B[i], &tmpQ, &tmpR, &e, &h, nvar, lazy);
		    		    
		    if (tmpR != NULL && tmpR->size != 0) {
				r = tmpR; // update dividend
			} else {
				r = NULL;
			}
		    
		    if (totalH == NULL || totalH->size == 0) {
				totalH = deepCopyPolynomial_AA (h);
		    } else if (h != NULL && h->size != 0){
				totalH = multiplyPolynomials_AA_inp (totalH, h, nvar);
		    }

		    if (h != NULL && h->size != 0 && !isOne_AA (h)){
				for (int j = 0; j < nSet; ++j){
					if (quoSet[j] != NULL && quoSet[j]->size != 0){
						quoSet[j] = multiplyPolynomials_AA_inp (quoSet[j], h, nvar);
					}
				}
		    }
		    
		    if (tmpQ != NULL && tmpQ->size != 0){
		        if (quoSet[i] == NULL || quoSet[i]->size == 0){
		            quoSet[i] = tmpQ;
		        } else {
		            quoSet[i] = addPolynomials_AA_inp (quoSet[i], tmpQ, nvar);
					freePolynomial_AA (tmpQ);
		        }
		    }
			
		    freePolynomial_AA (h);
		}	
	}
    
    if (totalH == NULL || totalH->size == 0) {
		totalH = makeConstPolynomial_AA (1, nvar, one);
    }

	mpq_clear (one);
	
    *hPow = totalH;
    *rem = r;
    
    return;
}

void normalizedTriangularSetPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t*** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet)
{
    if (nSet < 0){
        fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
        exit (EXIT_FAILURE);
	}

	AltArr_t** Q; // *quoSet = Q
	
    if (c == NULL || c->size == 0) {
		mpq_t one1;
		mpq_init (one1);
		mpq_set_si (one1, 1l, 1ul);

		*rem = NULL;
		*hPow = makeConstPolynomial_AA(1, nvar, one1);

		if (nSet > 0) {
			Q = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*));
			 
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
    	mpq_t one2;
		mpq_init (one2);
		mpq_set_si (one2, 1l, 1ul);
		
		*rem = deepCopyPolynomial_AA(c);
		*hPow = makeConstPolynomial_AA(1, nvar, one2);
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
    	fprintf(stderr, "normalizedTriangularSetPseudoDivide_AA not yet implemented with exponent unpacking!\n");
    	exit (EXIT_FAILURE);
    }

    AltArr_t* tmpQ = NULL;  //(AltArr_t*) calloc (1, sizeof(AltArr_t));
    AltArr_t* tmpR = NULL; //(AltArr_t*) calloc (1, sizeof(AltArr_t));
    AltArr_t* h = NULL;   //(AltArr_t*) calloc (1, sizeof(AltArr_t));

	mpq_t one;
	mpq_init (one); // was tmpCoef
	mpq_set_si (one, 1l, 1ul);
	
    int e = 1;
  	int mvar;
  	int signMvar; //?
    /* int* expOffset = getExpOffsetArray(nvar); // ? */

    if (nSet == 1) {
        mvar = leadingVariable_AA (B[0]);
		
		if (mvar == -2) {
			*rem = deepCopyPolynomial_AA (c);
            quoSet[0] = NULL;
            *hPow = makeConstPolynomial_AA(1, nvar, one);
            return;
		}
 		
		if (mvar > -1){
			pesudoDivideAtIdx_AA (mvar, c, B[0], &tmpQ, &tmpR, &e, &h, nvar, lazy);
		} else{
			pesudoDivideAtIdx_AA (0, c, B[0], &tmpQ, &tmpR, &e, &h, nvar, lazy);
		}
        
        *rem = tmpR;

		Q = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*));
		Q[0] = tmpQ;

		*quoSet = Q;
        *hPow = h;
        return;
    }

	Q = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*));
	
	/////////////// PART 1 ////////////////
	signMvar = leadingVariable_AA (B[nSet-1]);
	mvar = (signMvar < 0) ? 0 : signMvar;

	AltArr_t** CoefList1;
	int sz_cl1 = 0; // deg(c, mvar)+1
	mainCoefficientListAtIdx_AA (c, mvar, &CoefList1, &sz_cl1);

	if (!sz_cl1) {
		fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AA![Part 1]\n");
		exit (EXIT_FAILURE);
	}
	
    AltArr_t** Q1[sz_cl1];
	/* for (int i = 0; i < sz_cl1; i++) { */
	/* 	Q1[i] = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*)); */
	/* } */
	
	AltArr_t** R1 = (AltArr_t**) calloc (sz_cl1, sizeof (AltArr_t*));
	AltArr_t** HPow1 = (AltArr_t**) calloc (sz_cl1, sizeof (AltArr_t*));
	
	for (int i = 0; i < sz_cl1; i++) {
	    normalizedTriangularSetPseudoDivide_AA (CoefList1[i], B,
										   &Q1[i], &R1[i], &HPow1[i], nvar, lazy, nSet-1);
	}

	/* char * ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */

	/* fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* for (int i = 0; i < sz_cl1; i++) { //TEST */
	/* 	fprintf (stderr, "CoefList1[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, CoefList1[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "R1[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, R1[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "HPow1[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, HPow1[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */
	/* } */
	/* fprintf (stderr, "***** *****\n"); // TEST */

	for (int i = 0; i < sz_cl1; i++) {
		if (CoefList1[i] != NULL && CoefList1[i]->size != 0) {
			freePolynomial_AA (CoefList1[i]);
		}
	}
	free (CoefList1);

	AltArr_t* H1 = NULL;
	
	for (int i = 0; i < sz_cl1; i++){
		if (HPow1[i] != NULL && HPow1[i]->size != 0){ // TODO
			if (H1 == NULL || H1->size == 0){
				H1 = deepCopyPolynomial_AA (HPow1[i]);
			} else {
				H1 = maxPolynomials_AA_inp (H1, HPow1[i]); // TODO: change to lcm()
			}
		}
	}

	/* fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* fprintf (stderr, "H1 :=  "); */
	/* printPoly_AA_unpk(stderr, H1, ch, nvar); */
	/* fprintf (stderr, "\n"); */
	/* fprintf (stderr, "***** *****\n"); // TEST */
	
	AltArr_t* newR1 = NULL;
	AltArr_t* frac_r1;
	AltArr_t* frac_q1;

	for (int i = 0; i < sz_cl1; ++i){
		
		if (H1 == NULL || H1->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, H1 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}
		
		dividePolynomials_AA (H1, HPow1[i], &frac_q1, &frac_r1, nvar);

		if (frac_r1 != NULL || frac_r1->size != 0) {
			freePolynomial_AA (frac_r1);
		}
		
		if (frac_q1 == NULL || frac_q1->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, frac_q1 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}

		if (R1[i] != NULL && R1[i]->size != 0) {
			
			multiplyPolynomialAtIdxByXn_AA_inp (R1[i], mvar, i, nvar);		
			R1[i] = multiplyPolynomials_AA_inp (R1[i], frac_q1, nvar);
			newR1 = addPolynomials_AA_inp (newR1, R1[i], nvar);

		}
		
		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AA_inp (Q1[i][j], mvar, i, nvar);
			Q1[i][j] = multiplyPolynomials_AA_inp (Q1[i][j], frac_q1, nvar);
			if (i == 0) {
				Q[j] = deepCopyPolynomial_AA (Q1[i][j]);
			} else {
				Q[j] = addPolynomials_AA_inp (Q[j], Q1[i][j], nvar);
			}
		}
		
		freePolynomial_AA (frac_q1); frac_q1 = NULL;
	}

	/* fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* fprintf (stderr, "newR1 :=  "); */
	/* printPoly_AA_unpk(stderr, newR1, ch, nvar); */
	/* fprintf (stderr, "\n"); */
	/* fprintf (stderr, "***** *****\n"); // TEST */

	for (int i = 0; i < sz_cl1; ++i){
		if (R1[i] != NULL || R1[i]->size != 0) {
			freePolynomial_AA(R1[i]);
		}
		if (HPow1[i] != NULL || HPow1[i]->size != 0) {
			freePolynomial_AA(HPow1[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q1[i][j] != NULL || Q1[i][j]->size != 0) {
				freePolynomial_AA(Q1[i][j]);
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
	if (!isOne_AA (H1) && mainVariable_AA (newR1) > 1) {
		for (int i = 0; i < nSet-1; i++) {
			if (partialDegree_AA (newR1, mainVariable_AA (B[i])) >=
				partialDegree_AA (B[i], mainVariable_AA(B[i]))) {
				AltArr_t* r_star = NULL;
				AltArr_t* q_star = NULL;
				
				dividePolynomials_AA (newR1, B[i], &q_star, &r_star, nvar);

				if (q_star != NULL) {
						Q[i] = addPolynomials_AA_inp (Q[i], q_star, nvar);
						freePolynomial_AA (q_star);
				}
				if (r_star != NULL) {
					freePolynomial_AA (newR1);
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
	
	AltArr_t* q_tilde = NULL;
	AltArr_t* r_tilde = NULL;
	AltArr_t* h_tilde = NULL;
	int       e_tilde = 0;

	pesudoDivideAtIdx_AA (mvar, newR1, B[nSet-1], &q_tilde, &r_tilde, &e_tilde, &h_tilde, nvar, lazy);
	
	freePolynomial_AA (newR1);
	
	if (h_tilde == NULL || h_tilde->size == 0) {
		fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, h_tilde couldn't be zero!\n");
		exit (EXIT_FAILURE);
	}

	if (!isOne_AA (h_tilde)) {
		for (int j = 0; j < nSet-1; j++) {
			if (Q[j] != NULL && Q[j]->size != 0) {
				Q[j] = multiplyPolynomials_AA_inp (Q[j], h_tilde, nvar);
			}
		}
	}

	Q[nSet-1] = q_tilde;

	if (!isOne_AA (h_tilde)) {
		H1 = multiplyPolynomials_AA_inp (H1, h_tilde, nvar); 
	}
	
	freePolynomial_AA (h_tilde);
	
	/////////////// PART 3 ////////////////
	
	if (r_tilde == NULL || r_tilde->size == 0) {
		*quoSet = Q;
		*rem = NULL;
		*hPow = H1;

		return;
	}
	
	AltArr_t** CoefList2;
	int sz_cl2 = 0; // deg(r_tilde, mvar)+1
	mainCoefficientListAtIdx_AA (r_tilde, mvar, &CoefList2, &sz_cl2);

	if (!sz_cl2) {
		fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AA [Part 3]!\n");
		exit (EXIT_FAILURE);
	}
	
    AltArr_t** Q2[sz_cl2];
	/* for (int i = 0; i < sz_cl2; i++) { */
	/* 	Q2[i] = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*)); */
	/* } */
	
	AltArr_t** R2 = (AltArr_t**) calloc (sz_cl2, sizeof (AltArr_t*));
	AltArr_t** HPow2 = (AltArr_t**) calloc (sz_cl2, sizeof (AltArr_t*));
	
	for (int i = 0; i < sz_cl2; i++) {
	    normalizedTriangularSetPseudoDivide_AA (CoefList2[i], B,
										   &Q2[i], &R2[i], &HPow2[i], nvar, lazy, nSet-1);
	}

	/* for (int i = 0; i < sz_cl2; i++) { //TEST */
	/* 	fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* 	fprintf (stderr, "CoefList2[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, CoefList2[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "R2[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, R2[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "HPow2[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, HPow2[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */
	/* } */
	/* fprintf (stderr, "***** *****\n"); // TEST */

	
	for (int i = 0; i < sz_cl2; i++) {
		if (CoefList2[i] != NULL && CoefList2[i]->size != 0) {
			freePolynomial_AA (CoefList2[i]);
		}
	}
	free (CoefList2);
	
	AltArr_t* H2 = NULL;
	
	for (int i = 0; i < sz_cl2; i++){
		if (HPow2[i] != NULL && HPow2[i]->size != 0){ // TODO
			if (H2 == NULL || H2->size == 0){
				H2 = deepCopyPolynomial_AA (HPow2[i]);
			} else {
				H2 = maxPolynomials_AA_inp (H2, HPow2[i]); // TODO: change to lcm()
			}
		}
	}

	///// Update Quotients:
	if (!isOne_AA (H2)) {
		for (int j = 0; j < nSet; j++) {
			if (Q[j] != NULL && Q[j]->size != 0) {
				Q[j] = multiplyPolynomials_AA_inp (Q[j], H2, nvar);
			}	
		}
	}
	
	AltArr_t* newR2 = NULL;
	AltArr_t* frac_r2;
	AltArr_t* frac_q2;

	for (int i = 0; i < sz_cl2; ++i){
		
		if (H2 == NULL || H2->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, H2 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}
		
		dividePolynomials_AA (H2, HPow2[i], &frac_q2, &frac_r2, nvar);

		if (frac_r2 != NULL && frac_r2->size != 0) {
			freePolynomial_AA (frac_r2);
		}
		
		if (frac_q2 == NULL || frac_q2->size == 0){
			fprintf (stderr, "SMQP Error: In normalizedTriangularSetPseudoDivide_AA, frac_q2 couldn't be zero!\n");
			exit (EXIT_FAILURE);
		}

		if (R2[i] != NULL && R2[i]->size != 0) {
			
			multiplyPolynomialAtIdxByXn_AA_inp (R2[i], mvar, i, nvar);		
			R2[i] = multiplyPolynomials_AA_inp (R2[i], frac_q2, nvar);
			newR2 = addPolynomials_AA_inp (newR2, R2[i], nvar);

		}
		
		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AA_inp (Q2[i][j], mvar, i, nvar);
			Q2[i][j] = multiplyPolynomials_AA_inp (Q2[i][j], frac_q2, nvar);
			Q[j] = addPolynomials_AA_inp (Q[j], Q2[i][j], nvar);
		}
		
		freePolynomial_AA (frac_q2); frac_q2 = NULL;
	}
	
	for (int i = 0; i < sz_cl2; ++i){
		if (R2[i] != NULL || R2[i]->size != 0) {
			freePolynomial_AA(R2[i]);
		}
		if (HPow2[i] != NULL || HPow2[i]->size != 0) {
			freePolynomial_AA(HPow2[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q2[i][j] != NULL || Q2[i][j]->size != 0) {
				freePolynomial_AA(Q2[i][j]);
			}
		}
	}

	free (R2);
	// free (HPow2);
	// free (Q2);
	
	*quoSet = Q;
	*rem = newR2;
	*hPow = multiplyPolynomials_AA_inp (H1, H2, nvar);

	freePolynomial_AA (H2);
}

void multiDivisorPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet)
{
	if (nSet  < 0 || nvar < nSet){
		fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
		exit(EXIT_FAILURE);
    }
    
	
    if (nSet == 0){
		mpq_t one;
		mpq_init (one);
		mpq_set_si (one, 1l, 1ul);
		
		*rem = deepCopyPolynomial_AA(c);
		if (*hPow != NULL){
			freePolynomial_AA (*hPow);
		}
		*hPow = makeConstPolynomial_AA(1, nvar, one);
		quoSet[0] = NULL;

		mpq_clear (one);
		return;
    }

	// TODO: 
	int isNormalized = 0;

	if (isNormalized) {
		AltArr_t** Q;
		
		normalizedTriangularSetPseudoDivide_AA  (c, B, &Q, rem, hPow, nvar, lazy, nSet);

		for (int i = 0; i < nSet; i++) {
			quoSet[i] = Q[i];
		}
		
		return;
	}
	
	naiveMultiDivisorPseudoDivide_AA  (c, B, quoSet, rem, hPow, nvar, lazy, nSet);
}


void recursiveTriangularSetMDD_AA (AltArr_t* c, AltArr_t** B, AltArr_t*** quoSet, AltArr_t** rem, int nvar, int nSet)
{
    if (nSet < 0){
        fprintf(stderr, "SMQP Error: nSet(=%d) is out of range!\n", nSet);
        exit (EXIT_FAILURE);
	}

	AltArr_t** Q; // *quoSet = Q
	
    if (c == NULL || c->size == 0) {
		*rem = NULL;
		if (nSet > 0) {
			Q = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*));
			 
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
		*rem = deepCopyPolynomial_AA(c);
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
    	fprintf(stderr, "recursiveTriangularSetMDD_AA not yet implemented with exponent unpacking!\n");
    	exit (EXIT_FAILURE);
    }

    AltArr_t* tmpQ = NULL;
    AltArr_t* tmpR = NULL;

	/* mpq_t one; */
	/* mpq_init (one); // was tmpCoef */
	/* mpq_set_si (one, 1l, 1ul); */
	
  	int mvar;
  	int signMvar;

    if (nSet == 1) {
        mvar = leadingVariable_AA (B[0]);
		
		if (mvar == -2) {
			*rem = deepCopyPolynomial_AA (c);
            quoSet[0] = NULL;
            return;
		}
 		
		dividePolynomials_AA (c, B[0], &tmpQ, &tmpR, nvar);
        
        *rem = tmpR;

		Q = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*));
		Q[0] = tmpQ;
		*quoSet = Q;
		
        return;
    }

	Q = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*));

	/* int isLazard = 1; */
	/* if (isLazard) { */
	/* 	/\* lazardTriangularSetMDD_AA (c, B, Q, rem, nSet, nvar); *\/ */
	/* 	primitiveFactorTriangularSetMDD_AA (c, B, Q, rem, nSet, nvar); */
	/* 	*quoSet = Q; */
	/* 	return; */
	/* } */
	
	if (nSet < 4 || AA_SIZE (c) < 5) {
		normalForm_AA (c, B, Q, rem, nSet, nvar);
		*quoSet = Q;
		return;
	}

	/////////////// PART 1 ////////////////
	signMvar = leadingVariable_AA (B[nSet-1]);
	mvar = (signMvar < 0) ? 0 : signMvar;

	AltArr_t** CoefList1;
	int sz_cl1 = 0; // deg(c, mvar)+1
	mainCoefficientListAtIdx_AA (c, mvar, &CoefList1, &sz_cl1);

	if (!sz_cl1) {
		fprintf (stderr, "SMQP Error: In recursiveTriangularSetMDD_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AA![Part 1]\n");
		exit (EXIT_FAILURE);
	}
	
    AltArr_t** Q1[sz_cl1];
	/* for (int i = 0; i < sz_cl1; i++) { */
	/* 	Q1[i] = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*)); */
	/* } */
	
	AltArr_t** R1 = (AltArr_t**) calloc (sz_cl1, sizeof (AltArr_t*));
	
	for (int i = 0; i < sz_cl1; i++) {
		recursiveTriangularSetMDD_AA (CoefList1[i], B, &Q1[i], &R1[i], nvar, nSet-1);
	}

	/* char * ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */

	/* fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* for (int i = 0; i < sz_cl1; i++) { //TEST */
	/* 	fprintf (stderr, "CoefList1[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, CoefList1[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "R1[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, R1[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */
	/* } */
	/* fprintf (stderr, "***** *****\n"); // TEST */

	for (int i = 0; i < sz_cl1; i++) {
		if (CoefList1[i] != NULL && CoefList1[i]->size != 0) {
			freePolynomial_AA (CoefList1[i]);
		}
	}
	free (CoefList1);
	
	AltArr_t* newR1 = NULL;

	for (int i = 0; i < sz_cl1; ++i){		
		if (R1[i] != NULL && R1[i]->size != 0) {
			multiplyPolynomialAtIdxByXn_AA_inp (R1[i], mvar, i, nvar);		
			newR1 = addPolynomials_AA_inp (newR1, R1[i], nvar);	
		}
		
		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AA_inp (Q1[i][j], mvar, i, nvar);
			if (i == 0) {
				Q[j] = deepCopyPolynomial_AA (Q1[i][j]);
			} else {
				Q[j] = addPolynomials_AA_inp (Q[j], Q1[i][j], nvar);
			}
		}
	}

	/* fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* fprintf (stderr, "newR1 :=  "); */
	/* printPoly_AA_unpk(stderr, newR1, ch, nvar); */
	/* fprintf (stderr, "\n"); */
	/* fprintf (stderr, "***** *****\n"); // TEST */

	for (int i = 0; i < sz_cl1; ++i){
		if (R1[i] != NULL || R1[i]->size != 0) {
			freePolynomial_AA(R1[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q1[i][j] != NULL || Q1[i][j]->size != 0) {
				freePolynomial_AA(Q1[i][j]);
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
	
	AltArr_t* q_tilde = NULL;
	AltArr_t* r_tilde = NULL;

	if (monomialDivideTest_AA (newR1, 0, B[nSet-1], 0)) {
		dividePolynomials_AA (newR1, B[nSet-1], &q_tilde, &r_tilde, nvar);
		freePolynomial_AA (newR1);
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
	
	AltArr_t** CoefList2;
	int sz_cl2 = 0; // deg(r_tilde, mvar)+1
	mainCoefficientListAtIdx_AA (r_tilde, mvar, &CoefList2, &sz_cl2);

	if (!sz_cl2) {
		fprintf (stderr, "SMQP Error: In recursiveTriangularSetMDD_AA, input polynomial is not zero but it is considered as a zero polynomial in mainCoefficientListAtIdx_AA [Part 3]!\n");
		exit (EXIT_FAILURE);
	}
	
    AltArr_t** Q2[sz_cl2];
	/* for (int i = 0; i < sz_cl2; i++) { */
	/* 	Q2[i] = (AltArr_t**) calloc (nSet, sizeof (AltArr_t*)); */
	/* } */
	
	AltArr_t** R2 = (AltArr_t**) calloc (sz_cl2, sizeof (AltArr_t*));
	
	for (int i = 0; i < sz_cl2; i++) {
		recursiveTriangularSetMDD_AA (CoefList2[i], B, &Q2[i], &R2[i], nvar, nSet-1);
	}

	/* for (int i = 0; i < sz_cl2; i++) { //TEST */
	/* 	fprintf (stderr, "***** nvar = %d && nSet = %d *****\n", nvar, nSet); */
	/* 	fprintf (stderr, "CoefList2[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, CoefList2[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "R2[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, R2[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */

	/* 	fprintf (stderr, "HPow2[%d] :=  ", i); */
	/* 	printPoly_AA_unpk(stderr, HPow2[i], ch, nvar); */
	/* 	fprintf (stderr, "\n"); */
	/* } */
	/* fprintf (stderr, "***** *****\n"); // TEST */

	
	for (int i = 0; i < sz_cl2; i++) {
		if (CoefList2[i] != NULL && CoefList2[i]->size != 0) {
			freePolynomial_AA (CoefList2[i]);
		}
	}
	free (CoefList2);
		
	AltArr_t* newR2 = NULL;

	for (int i = 0; i < sz_cl2; ++i){
		
		if (R2[i] != NULL && R2[i]->size != 0) {
			
			multiplyPolynomialAtIdxByXn_AA_inp (R2[i], mvar, i, nvar);		
			newR2 = addPolynomials_AA_inp (newR2, R2[i], nvar);
		}
		
		for (int j = 0; j < nSet-1; j++) {
			multiplyPolynomialAtIdxByXn_AA_inp (Q2[i][j], mvar, i, nvar);
			Q[j] = addPolynomials_AA_inp (Q[j], Q2[i][j], nvar);
		}
	}
	
	for (int i = 0; i < sz_cl2; ++i){
		if (R2[i] != NULL || R2[i]->size != 0) {
			freePolynomial_AA(R2[i]);
		}
		for (int j = 0; j < nSet-1; ++j){
			if (Q2[i][j] != NULL || Q2[i][j]->size != 0) {
				freePolynomial_AA(Q2[i][j]);
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

AltArr_t* LazardOpt  (AltArr_t* Sd, AltArr_t* Sdm, AltArr_t* s)
{
    if (Sd == NULL || Sd->size == 0){
		fprintf (stderr, "SMQP Error: In LazardOpt , Sd is NULL!\n");
		exit (EXIT_FAILURE);
    }
    
    if (Sdm == NULL || Sdm->size == 0){
		return NULL;
    }
    
    if (Sd->nvar != Sdm->nvar){
		fprintf (stderr, "SMQP Error: In LazardOpt , Sd->nvar(=%d) != Sdm->nvar(=%d).\n",
			 Sd->nvar, Sdm->nvar);
		exit (EXIT_FAILURE);
    }

     // fprintf (stderr, "Hello Lazard"); 
     // fprintf (stderr, "Sd, Sdm = \n"); 
     // printAA  (Sd); 
     // printAA  (Sdm); 
     // fprintf (stderr, "\n\n"); 
    
    register int nvar = Sd->nvar;
    
    /* degree_t Sd_deg = (Sd->elems[0].degs & mvarMask) >> mvarDegOffset; */
    /* degree_t Sdm_deg = (Sdm->elems[0].degs & mvarMask) >> mvarDegOffset; */
    degree_t n = mainLeadingDegree_AA(Sd) - mainLeadingDegree_AA(Sdm) - 1;//((Sd->elems[0].degs & mvarMask) >> mvarDegOffset) -
	// ((Sdm->elems[0].degs & mvarMask) >> mvarDegOffset) - 1;
    if (n < 1){
		return deepCopyPolynomial_AA  (Sdm);
    }
    
    AltArr_t* __restrict__ x = mainLeadingCoefficient_AA (Sdm);
    AltArr_t* __restrict__ y = s; //mainLeadingCoefficient_AA (Sd);	    
    
    if (x == NULL || x->size == 0){
		fprintf (stderr, "SMQP Error: In LazardOpt, x is NULL!");
		exit (EXIT_FAILURE);
    }
    if (y == NULL || y->size == 0){
		fprintf (stderr, "SMQP Error: In LazardOpt, y is NULL!");
		exit (EXIT_FAILURE);
    }
    
    AltArr_t* c = deepCopyPolynomial_AA (x);
    register degree_t a = 1 << Log2_AA (n);
    n = n - a;
    
    AltArr_t* tmpq;
    AltArr_t* tmp;
    /* int orign = n; */
    
    while (a > 1) {
	tmpq = NULL;
	tmp = multiplyPolynomials_AA_inp (c, c, nvar);	
	exactDividePolynomials_AA (tmp, y, &tmpq, nvar);
	freePolynomial_AA (tmp);
	c = tmpq;    
	a = a >> 1;
	
	if (n >= a) {
	    tmpq = NULL;
	    if (c != NULL && c->size != 0)
			tmp = multiplyPolynomials_AA_inp (c, x, nvar);
	    else {
			tmp = NULL;
			freePolynomial_AA (c);
	    }
	    exactDividePolynomials_AA (tmp, y, &tmpq, nvar);
	    freePolynomial_AA (tmp);
	    c = tmpq;
	    n = n - a;
		}
    }
    
    tmpq = NULL;
    AltArr_t* Se = deepCopyPolynomial_AA (Sdm);
    Se = multiplyPolynomials_AA_inp (Se, c, nvar);
    freePolynomial_AA (c);
    freePolynomial_AA (x);
    /* if (!orign){ */
    /* 	return Se; */
    /* } */

    exactDividePolynomials_AA (Se, y, &tmpq, nvar);
    freePolynomial_AA (Se);
    return tmpq;
}

AltArr_t* DucosOpt (AltArr_t* A, AltArr_t* Sdm, AltArr_t* Se, AltArr_t* sd)
{
    int nvar = A->nvar;
    
    if (nvar != Sdm->nvar || nvar != Se->nvar || nvar != sd->nvar){
	fprintf(stderr, "SMQP ERROR: In DucosOpt, inputs' nvars are different!\n");
	exit (EXIT_FAILURE);
    }
    
    degree_t d = mainLeadingDegree_AA(A); //(A->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t e = mainLeadingDegree_AA(Sdm); //(Sdm->elems[0].degs & mvarMask) >> mvarDegOffset;
    AltArr_t* cdm = mainLeadingCoefficient_AA (Sdm);
    AltArr_t* se = mainLeadingCoefficient_AA (Se);
    AltArr_t* tmpq;
    AltArr_t* tmpR;
    AltArr_t* tmp = NULL;

    if (d == 0){
	fprintf (stderr, "SMQP Error: In DucosOpt, d == 0.\n");
	exit (EXIT_FAILURE);
    }
    if (e == 0 || d < e){
	fprintf (stderr, "SMQP Error: In DucosOpt, e == 0 || d > e.\n");
	exit (EXIT_FAILURE);
    }
    if (cdm == NULL || cdm->size ==0){
	fprintf (stderr, "SMQP Error: In DucosOpt, cdm = NULL.\n");
	exit (EXIT_FAILURE);
    }
    if (se == NULL || se->size ==0){
	fprintf (stderr, "SMQP Error: In DucosOpt, se = NULL.\n");
	exit (EXIT_FAILURE);
    }
    if (d == e){
	// non-defective cases:
	// Streamlined resultant computation for regular case
	AltArr_t* dc = mainLeadingCoefficient_AA (A); 
	AltArr_t* rTemp = NULL;
	AltArr_t* supTemp = NULL;
	
	// compute first of two main terms
	AltArr_t* res = deepCopyPolynomial_AA (A); // res = a
	res = multiplyPolynomials_AA_inp (res, se, nvar); // res *= ec
	rTemp = mainCoefficientAtIdx_AA(A, e); // rTemp = a.coef(e)
	supTemp = multiplyPolynomials_AA (rTemp, Se, nvar); // supTemp = se * rTerm
	if (supTemp != NULL && supTemp->size != 0){
	    res = subPolynomials_AA_inp (res, supTemp, nvar); // res -= supTemp
	}
	dividePolynomials_AA (res, dc, &tmp, &tmpR, nvar); // res /= dc
	freePolynomial_AA (tmpR); 
	res = multiplyPolynomials_AA (tmp, se, nvar); // res *= ec
	freePolynomial_AA (tmp);

	// compute second of two main terms
	rTemp = mainCoefficientAtIdx_AA (Se, e-1); // rTemp = se.coef(e-1)
	freePolynomial_AA (supTemp); // supTemp.zero()
	supTemp = rTemp; // supTemp.setCoef(0, rTemp)
	rTemp = deepCopyPolynomial_AA (se); 
	negatePolynomial_AA (rTemp); // rTemp = -ec
	rTemp = mainLShiftPolynomial_AA_inp (rTemp, 1);

	supTemp = addPolynomials_AA (supTemp, rTemp, nvar); // supTemp.setCoef(1, rTemp)
	supTemp = multiplyPolynomials_AA_inp (supTemp, Se, nvar); // supTemp *= se
	res = addPolynomials_AA_inp (res, supTemp, nvar); // res += supTemp
	freePolynomial_AA (supTemp);
	// divide out dc to obtain resultant
	dividePolynomials_AA (res, dc, &tmp, &tmpR, nvar); // res /= dc
	freePolynomial_AA (tmpR); 
	return tmp; // return res 
    }
    // defective cases:
    // H = 
    AltArr_t* H[d];
    for (int j = 0; j <= e; ++j){
		H[j] = deepCopyPolynomial_AA (se);
		H[j] = mainLShiftPolynomial_AA_inp (H[j], j);
    }
    
    H[e] = subPolynomials_AA_inp (H[e], Se, nvar);
    
    AltArr_t* PIe = NULL;
    for (int j = e+1; j < d; ++j){
	//XH_{j-1} = 
	H[j] = deepCopyPolynomial_AA (H[j-1]);	    
	
	mainLShiftPolynomial_AA_inp(H[j], 1);
	
	//syzygy:
	PIe = mainCoefficientAtIdx_AA (H[j], e);
	
	if (PIe != NULL && PIe->size != 0){
	    tmp = multiplyPolynomials_AA_inp (PIe, Sdm, nvar);
	    tmpq = NULL;
	    tmpR = NULL;
	    dividePolynomials_AA (tmp, cdm, &tmpq, &tmpR, nvar);
	    freePolynomial_AA (tmpR);
	    freePolynomial_AA (tmp); //  freePolynomial_AA (PIe);
	    H[j] = subPolynomials_AA_inp (H[j], tmpq, nvar);
	    freePolynomial_AA (tmpq);
	}
    }
    
    // D =
    AltArr_t* PIj = NULL;
    AltArr_t* D = NULL;	
    for (int j = 0; j < d; ++j){
	PIj = mainCoefficientAtIdx_AA (A, j);
	if (PIj != NULL && PIj->size != 0 &&
	    H[j] != NULL && H[j]->size != 0){
	    PIj = multiplyPolynomials_AA_inp (PIj, H[j], nvar);
	    
	    if (D == NULL || D->size == 0){
		D = deepCopyPolynomial_AA (PIj);
	    } else {
		D = addPolynomials_AA_inp (D, PIj, nvar);
	    }
	    freePolynomial_AA (PIj);
	}
    }
    
    tmpq = NULL;
    tmpR = NULL;
    AltArr_t* lcA = mainLeadingCoefficient_AA (A);
    if (lcA != NULL && lcA->size != 0){
	dividePolynomials_AA (D, lcA, &tmpq, &tmpR, nvar);
	freePolynomial_AA (tmpR);
	freePolynomial_AA (D);
	D = tmpq;
    }
    freePolynomial_AA (lcA);
    
    // res = 
    AltArr_t* res = NULL;
    AltArr_t* lPoly = deepCopyPolynomial_AA (H[d-1]);
    lPoly = mainLShiftPolynomial_AA_inp(lPoly, 1);
    
    AltArr_t* rPoly = mainCoefficientAtIdx_AA (lPoly, e);
    if (D != NULL && D->size != 0){
	lPoly = addPolynomials_AA_inp (lPoly, D, nvar);
    }
    lPoly = multiplyPolynomials_AA_inp (lPoly, cdm, nvar);
    
    if (rPoly != NULL && rPoly->size != 0){
	rPoly = multiplyPolynomials_AA_inp (rPoly, Sdm, nvar);
	res = subPolynomials_AA (lPoly, rPoly, nvar);
	freePolynomial_AA (lPoly);
	freePolynomial_AA (rPoly);
    } else {
	res = lPoly;
    }
    
    tmpR = NULL;
    tmpq = NULL;
    if (sd != NULL && sd->size != 0){
	    dividePolynomials_AA (res, sd, &tmpq, &tmpR, nvar);
    } else {
	tmpq = res;
    }
    freePolynomial_AA (tmpR);
    
    if ((d-e+1)%2 != 0 && tmpq != NULL && tmpq->size != 0){
	negatePolynomial_AA (tmpq);
    }

    freePolynomial_AA (cdm);
    freePolynomial_AA (se);
    
    return tmpq;    
}

AltArr_t* CFDucosOpt (AltArr_t* A, AltArr_t* Sdm, AltArr_t* Se, AltArr_t* sd)
{

    register int nvar = A->nvar;
    
    if (nvar != Sdm->nvar || nvar != Se->nvar || nvar != sd->nvar){
	fprintf (stderr, "SMQP ERROR: In CFDucosOpt , inputs' nvars are different!\n");
	exit (EXIT_FAILURE);
    }
    
    register degree_t d = mainLeadingDegree_AA(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
    register degree_t e = mainLeadingDegree_AA(Sdm);//(Sdm->elems[0].degs & mvarMask) >> mvarDegOffset;
    
    if (d == 0){
	fprintf (stderr, "SMQP Error: In CFDucosOpt , d == 0.\n");
	exit (EXIT_FAILURE);
    }
    if (e == 0 || d < e){
	fprintf (stderr, "SMQP Error: In CFDucosOpt , e == 0 || d > e.\n");
	exit (EXIT_FAILURE);
    }

    AltArr_t* __restrict__ cdm = mainLeadingCoefficient_AA (Sdm);
    AltArr_t* __restrict__ se = mainLeadingCoefficient_AA (Se);

    AltArr_t* tmpq = NULL;
    AltArr_t* tmpR = NULL;
    AltArr_t* tmp = NULL;
    AltArr_t* res = NULL;
    AltArr_t* PIe = NULL;
    AltArr_t* D2 = NULL;
    AltArr_t* PIj = NULL;
    AltArr_t* D = NULL;	
    AltArr_t* H = NULL;
    
    if (cdm == NULL || cdm->size ==0){
		fprintf (stderr, "SMQP Error: In CFDucosOpt , cdm = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (se == NULL || se->size ==0){
		fprintf (stderr, "SMQP Error: In CFDucosOpt , se = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (d == e){
	// for non-defective cases:
		fprintf (stderr, "SMQP Error: In CFDucosOpt , d == e.\n"); // TEST

	/* // Streamlined resultant computation for regular case */
	/* AltArr_t* dc = mainLeadingCoefficient_AA (A); */
	/* AltArr_t* rTemp = NULL; */
	/* AltArr_t* supTemp = NULL; */
	
	/* // compute first of two main terms */
	/* AltArr_t* res = deepCopyPolynomial_AA (A); // res = a */
	/* res = multiplyPolynomials_AA_inp (res, se, nvar); // res *= ec */
	/* rTemp = mainCoefficientAtIdx_AA(A, e); // rTemp = a.coef(e) */
	/* supTemp = multiplyPolynomials_AA (rTemp, Se, nvar); // supTemp = se * rTerm */
	/* if (supTemp != NULL && supTemp->size != 0){ */
	/*     res = subPolynomials_AA_inp (res, supTemp, nvar); // res -= supTemp */
	/* } */
	/* /\* dividePolynomials_AA (res, dc, &tmp, &tmpR, nvar); // res /= dc *\/ */
	/* exactDividePolynomials_AA (res, dc, &tmp, nvar); // res /= dc */
	/* // freePolynomial_AA (tmpR); */
	/* res = multiplyPolynomials_AA (tmp, se, nvar); // res *= ec */
	/* freePolynomial_AA (tmp); */

	/* // compute second of two main terms */
	/* rTemp = mainCoefficientAtIdx_AA (Se, e-1); // rTemp = se.coef(e-1) */
	/* freePolynomial_AA (supTemp); // supTemp.zero() */
	/* supTemp = rTemp; // supTemp.setCoef(0, rTemp) */
	/* rTemp = deepCopyPolynomial_AA (se); */
	/* negatePolynomial_AA (rTemp); // rTemp = -ec */
	/* for (int i = 0; rTemp != NULL && i < rTemp->size; ++i){ */
	/*     rTemp->elems[i].degs += ((degrees_t) 1 << mvarDegOffset); */
	/* } */
	/* supTemp = addPolynomials_AA (supTemp, rTemp, nvar); // supTemp.setCoef(1, rTemp) */
	/* supTemp = multiplyPolynomials_AA_inp (supTemp, Se, nvar); // supTemp *= se */
	/* res = addPolynomials_AA_inp (res, supTemp, nvar); // res += supTemp */
	/* freePolynomial_AA (supTemp); */
	/* // divide out dc to obtain resultant */
	/* exactDividePolynomials_AA (res, dc, &tmp, nvar); // res /= dc */
	/* return tmp; // return res */
    }
    
    // H[e] = 
    H = mainLShiftPolynomial_AA (se, e);
    H = subPolynomials_AA_inp (H, Se, nvar);
    
    // D2 = 
    PIe = mainCoefficientAtIdx_AA (A, e);    
    if (PIe != NULL && H != NULL && PIe->size != 0 && H->size != 0){
		D2 = multiplyPolynomials_AA (H, PIe, nvar);
		freePolynomial_AA (PIe);
    } 

    // PIe = mainCoefficientAtIdx_AA (H, e-1);  // it must be uncommented for using CFDucosOptZ_subPolynomials_AAZ_inp (see CFDucosOptZ)
    for (int j = e+1; j < d; ++j){
	//  H = XH - pi_e (XH)*C
	H = mainLShiftPolynomial_AA_inp (H, 1); // H = XH
	PIe = mainCoefficientAtIdx_AA (H, e);	// it must be commented for using CFDucosOptZ_subPolynomials_AAZ_inp (see CFDucosOptZ)
	if (PIe != NULL && PIe->size != 0){
	    tmp = multiplyPolynomials_AA_inp (PIe, Sdm, nvar);
	    tmpq = NULL;
	    exactDividePolynomials_AA (tmp, cdm, &tmpq, nvar);
	    freePolynomial_AA (tmp); //  freePolynomial_AA (PIe);
	    H = subPolynomials_AA_inp (H, tmpq, nvar);
	    // H = CFDucosOpt_subPolynomials_AA_inp (H, tmpq, nvar, &PIe, e-1);
	    freePolynomial_AA (tmpq);
	}
	tmp = mainCoefficientAtIdx_AA (A, j);
	if (tmp != NULL && H != NULL && tmp->size != 0 && H->size != 0){
	    tmp = multiplyPolynomials_AA_inp (tmp, H, nvar);
	    D2 = addPolynomials_AA_inp (D2, tmp, nvar);
	    freePolynomial_AA (tmp);
	}
    }
    
    // D =
    for (int j = 0; j < e; ++j){
		PIj = mainCoefficientAtIdx_AA (A, j);
		if (PIj != NULL && PIj->size != 0){
		    PIj = mainLShiftPolynomial_AA_inp (PIj, j);
		    PIj = multiplyPolynomials_AA_inp (PIj, se, nvar);
		}
		if (D == NULL || D->size == 0){
		    D = deepCopyPolynomial_AA (PIj);
		} else {
		    D = addPolynomials_AA_inp (D, PIj, nvar);
		}
		freePolynomial_AA (PIj);
    }
    
    if (D2 != NULL && D2->size != 0){
		D = addPolynomials_AA_inp (D, D2, nvar);
		freePolynomial_AA (D2);
    }    
    
    tmpq = NULL;
    tmpR = NULL;
    PIe = mainLeadingCoefficient_AA (A);
    if (PIe != NULL && PIe->size != 0){
	dividePolynomials_AA (D, PIe, &tmpq, &tmpR, nvar);
	freePolynomial_AA (tmpR);
	freePolynomial_AA (D);
	D = tmpq;
    }
    freePolynomial_AA (PIe);
    
    AltArr_t* lPoly = mainLShiftPolynomial_AA (H, 1);
    AltArr_t* rPoly = mainCoefficientAtIdx_AA (lPoly, e);
    if (D != NULL && D->size != 0){
	lPoly = addPolynomials_AA_inp (lPoly, D, nvar);
    }
    if (cdm != NULL){
	lPoly = multiplyPolynomials_AA_inp (lPoly, cdm, nvar);
    }
	
    if (rPoly != NULL && rPoly->size != 0){
	rPoly = multiplyPolynomials_AA_inp (rPoly, Sdm, nvar);
	res = subPolynomials_AA_inp (lPoly, rPoly, nvar);
	//freePolynomial_AA (lPoly);
	freePolynomial_AA (rPoly);
    } else {
	res = lPoly;
    }
	
    tmpR = NULL;
    tmpq = NULL;
    if (sd != NULL && sd->size != 0){
	dividePolynomials_AA (res, sd, &tmpq, &tmpR, nvar);
	freePolynomial_AA (tmpR);
	freePolynomial_AA (res);
    } else {
	tmpq = res;
    }
    
    if ((d-e+1)%2 != 0 && tmpq != NULL && tmpq->size != 0){
	negatePolynomial_AA (tmpq);
    }
    
    freePolynomial_AA (cdm);
    freePolynomial_AA (se);

    return tmpq;    
}

void DucosSubresultantChain_rev (AltArr_t* P, AltArr_t* Q, AltArrs_t** SC, int* len, int type)
{
    if (P == NULL || P->size == 0){
		AltArrs_t* SS0 = (AltArrs_t*) malloc (sizeof(AltArrs_t));
		SS0->poly = NULL;
		SS0->next = NULL;
		*SC = SS0;
		*len = 0;
		return;
    }
    if (Q == NULL || Q->size == 0){
		AltArrs_t* SS1 = (AltArrs_t*) malloc (sizeof(AltArrs_t));
		SS1->poly = deepCopyPolynomial_AA (P);
		SS1->next = NULL;
		*SC = SS1;
		*len = 1;
		return;
    }
    if (P->nvar != Q->nvar){
		fprintf (stderr, "SMQP Error: In DucosSubresultantChain_rev, P->nvar(%d) != Q->nvar(%d).\n",P->nvar, Q->nvar);
		exit (EXIT_FAILURE);
    }
    

  	degree_t degP = mainLeadingDegree_AA(P);
    degree_t degQ = mainLeadingDegree_AA(Q);

    // adding corner cases for supprt algebraic extensions
    // TODO: revise for efficiency
    if (degP == 0 && degQ == 0) {
    	mpq_t one;
    	mpq_init (one);
    	mpq_set_si (one, 1l, 1ul);
    	AltArrs_t* SS2 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    	SS2->poly = makeConstPolynomial_AA (1, P->nvar, one);
    	SS2->next = NULL;
    	*SC = SS2;
    	*len = 1;
    	return;
    } else if (degP == 0) {
    	AltArrs_t* SS3 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    	SS3->poly = exponentiatePoly_AA (P, degQ, P->nvar);
    	SS3->next = NULL;
    	*SC = SS3;
    	*len = 1;
    	return;
    } else if (degQ == 0) {
    	AltArrs_t* SS4 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    	SS4->poly = exponentiatePoly_AA (Q, degP, Q->nvar);
    	SS4->next = NULL;
    	*SC = SS4;
    	*len = 1;
    	return;
    }

    int lazy = 0; // lazy option in pseudoDivide algorithm
    int nvar = P->nvar;
    
    int tmpE = 0;
    AltArr_t* tmpq = NULL;
    AltArr_t* tmpH = NULL;
    AltArr_t* tmpB = NULL;
    AltArr_t* s = NULL;
    AltArr_t* p = NULL;
    AltArr_t* q = NULL;
    AltArr_t* A = NULL;
    AltArr_t* B = NULL;
    AltArr_t* C = NULL;

    AltArrs_t* SS = NULL;
    AltArrs_t* tail = NULL;
    AltArr_t* tailm = NULL;
    degree_t delta = 0;
    int size = 0;
      
    if (degP >= degQ){
		p = P;
		q = Q;
    } else {
		p = Q;
		q = P;
    
    	degP = mainLeadingDegree_AA(p);
    	degQ = mainLeadingDegree_AA(q);
    }
    
    s = mainLeadingCoefficient_AA (q);
    s = exponentiatePoly_AA (s, degP - degQ, nvar);
    A = deepCopyPolynomial_AA (q);
    

    // fprintf(stderr, "degP: %d, degQ: %d\n", degP, degQ); // TEST

    
    // B =
    AltArr_t* negQ = deepCopyPolynomial_AA (q);
    negatePolynomial_AA (negQ);
    
    // deepCopyPolynomial_AA(p), deepCopyPolynomial_AA(negQ)
    pesudoDivideAtIdx_AA (0, p, negQ, &tmpq, &B, &tmpE, &tmpH, nvar, lazy);
    freePolynomial_AA (tmpq);
    freePolynomial_AA (tmpH);
    freePolynomial_AA (negQ);


    // gmp_fprintf(stderr, "rem: %Qd*x^%lld\n", B->elems->coef, B->elems->degs); // TEST


    AltArrs_t* newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
    newPoly->poly = deepCopyPolynomial_AA (p);
    newPoly->next = NULL;

    SS = newPoly;
    tail = newPoly;
    size++;
    
    newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
    if (B == NULL || B->size == 0){
		newPoly->poly = A;
    } else{
		newPoly->poly = deepCopyPolynomial_AA (A);
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

		d = mainLeadingDegree_AA(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
		e = mainLeadingDegree_AA(B);//(B->elems[0].degs & mvarMask) >> mvarDegOffset;
		
		newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
		newPoly->poly = deepCopyPolynomial_AA (B);
		newPoly->next = NULL;

		tailm = tail->poly;
		tail->next = newPoly;
		tail = newPoly;
		size++;

		delta = d - e;
		if (delta > 1){

			C = LazardOpt (tailm, tail->poly, s); 
			newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
			newPoly->poly = deepCopyPolynomial_AA (C);
			newPoly->next = NULL;

			tail->next = newPoly;
			tail = newPoly;
			size++;

		} else {
			C = deepCopyPolynomial_AA (B);
		}

		if (e == 0){ // B.degree() == 0
			break;
		}
		
		if (type == 0){
			tmpB = DucosOpt (A, B, C, s);
		} else {
			tmpB = CFDucosOpt (A, B, C, s);
		}

		freePolynomial_AA (B);
		freePolynomial_AA (A);
		B = tmpB;
		A = deepCopyPolynomial_AA (C);
		freePolynomial_AA (s);
		s = mainLeadingCoefficient_AA (A);	
		freePolynomial_AA (C);
	}

	// if (B != NULL && B->size != 0) {
	// 	freePolynomial_AA (B);
	// }

	// if (C != NULL && C->size != 0) {
	// 	freePolynomial_AA (C);
	// }

	// freePolynomial_AA (A);
	// freePolynomial_AA (s);

    // if subresultant is 0, add it to SS:
    if (tail->poly != NULL && tail->poly->size != 0 && mainLeadingDegree_AA(tail->poly) > 0) {	
		newPoly = (AltArrs_t*) malloc (sizeof (AltArrs_t));
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

void DucosSubresultantChain (AltArr_t* P, AltArr_t* Q, AltArrs_t** SC, int* len)
{
    AltArrs_t* invSC;
    int size = 0;
    DucosSubresultantChain_rev (P, Q, &invSC, &size, 1);

    // reverse the subresultantchain:
    if (size > 1){
	AltArrs_t* curr = invSC;
	AltArrs_t* prev = NULL;
	AltArrs_t* next = NULL;
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

AltArr_t* DucosResultant  (AltArr_t* P, AltArr_t* Q)
{
    AltArrs_t* SC;
    int size = 0;
    DucosSubresultantChain (P, Q, &SC, &size);
    if (size){
    	AltArr_t* res = deepCopyPolynomial_AA (SC->poly);
		freeAltArrs (SC);
		return res;
    }

    freeAltArrs (SC);
    return NULL;
}

AltArr_t* DucosGCD  (AltArr_t* P, AltArr_t* Q)
{
    if (P == NULL || P->size == 0){
    	return deepCopyPolynomial_AA (Q);
    }
    if (Q == NULL || Q->size == 0){
    	return deepCopyPolynomial_AA (P);
    }
    if (P->nvar != Q->nvar){
    	fprintf (stderr, "SMQP Error: In DucosGCD , P->nvar(=%d) != Q->nvar(=%d).\n",
    		 P->nvar, Q->nvar);
    	exit (EXIT_FAILURE);
    }

    int nvar = P->nvar;
    
    mpq_t one;
    mpq_init(one);
    mpq_set_si (one, 1l, 1l);
    
    degree_t degP = mainLeadingDegree_AA(P);// (P->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t degQ = mainLeadingDegree_AA(Q);//(Q->elems[0].degs & mvarMask) >> mvarDegOffset;
    
    if (degP == 0 || degQ == 0){
    	return makeConstPolynomial_AA (1, nvar, one);
    }
    
    AltArrs_t* SC;
    int size = 0;
    DucosSubresultantChain (P, Q, &SC, &size);
    AltArrs_t* cur = SC;

    if (cur->poly == NULL || cur->poly->size == 0){
    	freeAltArrs (SC);
		return makeConstPolynomial_AA (1, nvar, one);
    }
    
    AltArr_t* res;
    degree_t nd;
    degree_t nnd;
    if (size > 2){
		nd = mainLeadingDegree_AA(cur->poly);//(cur->poly->elems[0].degs & mvarMask) >> mvarDegOffset;
		nnd = mainLeadingDegree_AA(cur->next->poly);//(cur->next->poly->elems[0].degs & mvarMask) >> mvarDegOffset;
	
		while ( nd == nnd ){
		    if (cur->next == NULL || cur->next->poly == NULL){
			break;
		    }
		    cur = cur->next;
		    nnd = mainLeadingDegree_AA(cur->poly);//(cur->poly->elems[0].degs & mvarMask) >> mvarDegOffset;
		}
		if (nd == 0){
			freeAltArrs (SC);
		    return makeConstPolynomial_AA (1, nvar, one);
		} else {
			res = deepCopyPolynomial_AA (cur->poly);
			freeAltArrs (SC);
		    return res;
		}
    }

    res = deepCopyPolynomial_AA (cur->next->poly);
    freeAltArrs (SC);

    return res;
}

AltArr_t* lastNonZeroChain_AA (AltArr_t* P, AltArr_t* Q)
{
    if (P == NULL || P->size == 0){
    	return deepCopyPolynomial_AA (Q);
    }
    if (Q == NULL || Q->size == 0){
    	return deepCopyPolynomial_AA (P);
    }
    if (P->nvar != Q->nvar){
    	fprintf (stderr, "SMQP Error: In lastNonZeroChain, P->nvar(=%d) != Q->nvar(=%d).\n",
    		 P->nvar, Q->nvar);
    	exit (EXIT_FAILURE);
    }

    int nvar = P->nvar;
    
    mpq_t one;
    mpq_init(one);
    mpq_set_si (one, 1l, 1l);
    
    degree_t degP = mainLeadingDegree_AA(P);//(P->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t degQ = mainLeadingDegree_AA(Q);//(Q->elems[0].degs & mvarMask) >> mvarDegOffset;

    if (nvar == 1){
		return univariateGCD_AA (P, Q);
    }
    
    if (degP == 0 || degQ == 0){
    	return makeConstPolynomial_AA (1, nvar, one);
    }

    
    AltArrs_t* SC;
    int size = 0, idx = 0;
    DucosSubresultantChain (P, Q, &SC, &size);
    
    AltArrs_t* cur = SC;
    
    while (idx < size) {
	if (cur->poly != NULL && cur->poly->size != 0){
	    break;
	}
	cur = cur->next;
	++idx;
    }
    AltArr_t* res = deepCopyPolynomial_AA (cur->poly);
    freeAltArrs  (SC);
    
    return res;
}

///////////////////////////////
// Extended Subresultant Chain 
///////////////////////////////

AltArr_t* semiLazardOpt_AA  (AltArr_t* Sd, AltArr_t* Sdm, AltArr_t* s)
{
    if (Sd == NULL || Sd->size == 0){
	fprintf (stderr, "SMQP Error: In Semi-LazardOpt , Sd is NULL!\n");
	exit (EXIT_FAILURE);
    }
    
    if (Sdm == NULL || Sdm->size == 0){
	return NULL;
    }
    
    if (Sd->nvar != Sdm->nvar){
	fprintf (stderr, "SMQP Error: In Semi-LazardOpt , Sd->nvar(=%d) != Sdm->nvar(=%d).\n",
		 Sd->nvar, Sdm->nvar);
	exit (EXIT_FAILURE);
    }

    register int nvar = Sd->nvar;
    
    degree_t n = mainLeadingDegree_AA(Sd) - mainLeadingDegree_AA(Sdm) - 1;

    if (n < 1){
    	mpq_t one;
    	mpq_init (one);
    	mpq_set_si (one, 1l, 1lu);
		return  makeConstPolynomial_AA (1, Sdm->nvar, one); // deepCopyPolynomial_AA  (Sdm);
    }
    
    AltArr_t* x = mainLeadingCoefficient_AA (Sdm);
    AltArr_t* y = s;     
    
    if (x == NULL || x->size == 0){
		fprintf (stderr, "SMQP Error: In Semi-LazardOpt, x is NULL!");
		exit (EXIT_FAILURE);
    }
    if (y == NULL || y->size == 0){
		fprintf (stderr, "SMQP Error: In Semi-LazardOpt, y is NULL!");
		exit (EXIT_FAILURE);
    }
    
    AltArr_t* c = deepCopyPolynomial_AA (x);
    register degree_t a = 1 << Log2_AA (n);
    n = n - a;
    
    AltArr_t* tmpq;
    AltArr_t* tmp;
    /* int orign = n; */
    
    while (a > 1) {
		tmpq = NULL;
		tmp = multiplyPolynomials_AA_inp (c, c, nvar);	
		exactDividePolynomials_AA (tmp, y, &tmpq, nvar);
		freePolynomial_AA (tmp);
		c = tmpq;    
		a = a >> 1;
	
		if (n >= a) {
		    tmpq = NULL;
		    if (c != NULL && c->size != 0) {
				tmp = multiplyPolynomials_AA_inp (c, x, nvar);
		    } else {
				tmp = NULL;
		    }
		    exactDividePolynomials_AA (tmp, y, &tmpq, nvar);
		    freePolynomial_AA (tmp);
		    c = tmpq;
		    n = n - a;
		}
    }

    /* if (!orign){ */
    /* 	return Se; */
    /* } */
    
    AltArr_t* sse = NULL; // semi-Se
	exactDividePolynomials_AA (c, y, &sse, nvar);   
    freePolynomial_AA (c);
    freePolynomial_AA (x); 

    return sse;
}

AltArr_t* exLazardOpt_AA (AltArr_t* Sd, exgcds_t* VSdm, AltArr_t* s, AltArr_t** hc, AltArr_t** qc) 
{

	AltArr_t* sse = semiLazardOpt_AA (Sd, VSdm->r, s); // semi-Se

	if (sse == NULL || sse->size == 0) {
		*hc = NULL;
		*qc = NULL;
		return NULL;
	}

	if (isOne_AA (sse)) {
		freePolynomial_AA (sse);
		
		*hc = deepCopyPolynomial_AA (VSdm->a);
		*qc = deepCopyPolynomial_AA (VSdm->b);
		return deepCopyPolynomial_AA (VSdm->r);
	}

	AltArr_t* hhc = multiplyPolynomials_AA (sse, VSdm->a, sse->nvar);
	AltArr_t* qqc = multiplyPolynomials_AA (sse, VSdm->b, sse->nvar);
	sse = multiplyPolynomials_AA_inp (sse, VSdm->r, sse->nvar); // Se

	*hc = hhc;
	*qc = qqc;
	return sse;
}


AltArr_t* exDucosOpt_AA (exgcds_t* VSd, exgcds_t* VSdm, AltArr_t* Se, AltArr_t* sd, AltArr_t** hb, AltArr_t** qb)
{
    int nvar = VSd->r->nvar;
    
    if (nvar != VSdm->r->nvar || nvar != Se->nvar || nvar != sd->nvar){
		fprintf (stderr, "SMQP ERROR: In Extended DucosOpt , inputs' nvars are different!\n");
		exit (EXIT_FAILURE);
    }
    
    degree_t d = mainLeadingDegree_AA(VSd->r);
    degree_t e = mainLeadingDegree_AA(VSdm->r);
    
    if (d == 0){
		fprintf (stderr, "SMQP Error: In Extended DucosOpt , d == 0.\n");
		exit (EXIT_FAILURE);
    }
    if (e == 0 || d < e){
		fprintf (stderr, "SMQP Error: In Extended DucosOpt , e == 0 || d > e.\n");
		exit (EXIT_FAILURE);
    }

	AltArr_t* cd  = mainLeadingCoefficient_AA (VSd->r);
    AltArr_t* cdm = mainLeadingCoefficient_AA (VSdm->r);
    AltArr_t* se  = mainLeadingCoefficient_AA (Se);

  	if (cd == NULL || cd->size ==0){
		fprintf (stderr, "SMQP Error: In Extended DucosOpt , cd = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (cdm == NULL || cdm->size ==0){
		fprintf (stderr, "SMQP Error: In Extended DucosOpt , cdm = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (se == NULL || se->size ==0){
		fprintf (stderr, "SMQP Error: In Extended DucosOpt , se = NULL.\n");
		exit (EXIT_FAILURE);
    }
    if (d == e){
		// for non-defective cases:
		fprintf (stderr, "SMQP Error: In Extended DucosOpt , d == e.\n"); // TEST
	}

	AltArr_t* ssd  = multiplyPolynomials_AA (sd, cd, nvar);
	AltArr_t* csdm = multiplyPolynomials_AA (cdm, se, nvar);
	AltArr_t* mul;
	if (csdm != NULL && VSd->r != NULL) {
		mul = multiplyPolynomials_AA (csdm, VSd->r, nvar);
	} else {
		mul = NULL;
	}

	AltArr_t* r = NULL;
	AltArr_t* q = NULL;
	dividePolynomials_AA (mul, VSdm->r, &q, &r, nvar);

	freePolynomial_AA (mul);
	freePolynomial_AA (se);
	freePolynomial_AA (cd);
	freePolynomial_AA (cdm);

	// fprintf(stderr, "[SMQP-Ali] Sem... \n");

	// Sem:
	AltArr_t* Sem = NULL;
	AltArr_t* tmp = NULL;
	if (r != NULL && r->size != 0) {
		dividePolynomials_AA (r, ssd, &Sem, &tmp, nvar);
		if (tmp != NULL){
			freePolynomial_AA (tmp); tmp = NULL;
		}
		freePolynomial_AA (r);
	}  

	AltArr_t* hhb;
	AltArr_t* qqb;

	AltArr_t* tmp1 = NULL; 
	AltArr_t* tmp2 = NULL;

	// char* ch[11] = {"x", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j"}; // TEST 

	// fprintf (stderr, "[SMQP-Ali] Sem := ");
	// printPoly_AA_unpk(stderr, Sem, ch, VSd->r->nvar);
	// fprintf (stderr, "\n");


	// fprintf(stderr, "[SMQP-Ali] hb... \n");

	// hb:
	if (VSd->a != NULL) {
		tmp = multiplyPolynomials_AA (csdm, VSd->a, nvar);
	}

	if (VSdm->a != NULL) {
		tmp1 = multiplyPolynomials_AA (q, VSdm->a, nvar);
	}

	tmp = subPolynomials_AA_inp (tmp, tmp1, nvar);
	freePolynomial_AA (tmp1); tmp1 = NULL;

	if(tmp == NULL || tmp->size == 0) {
		hhb = NULL;
	} else {
		dividePolynomials_AA (tmp, ssd, &tmp1, &tmp2, nvar);
		hhb = tmp1;
		if (tmp2 != NULL) {
			freePolynomial_AA (tmp2);
		}
		freePolynomial_AA (tmp);
	}

	// fprintf (stderr, "[SMQP-Ali] hb := ");
	// printPoly_AA_unpk(stderr, hb, ch, VSd->r->nvar);
	// fprintf (stderr, "\n");

	// fprintf(stderr, "[SMQP-Ali] qb... \n");

	// qb:
	tmp = NULL;
	if (VSd->b != NULL) {
		tmp = multiplyPolynomials_AA (csdm, VSd->b, nvar);
	}

	tmp1 = NULL;
	if (VSdm->b != NULL) {
		tmp1 = multiplyPolynomials_AA (q, VSdm->b, nvar);
	}

	tmp = subPolynomials_AA_inp (tmp, tmp1, nvar);
	freePolynomial_AA (tmp1); tmp1 = NULL;

	if(tmp == NULL || tmp->size == 0) {
		qqb = NULL;
	} else {
		dividePolynomials_AA (tmp, ssd, &tmp1, &tmp2, nvar);
		qqb = tmp1;
		if (tmp2 != NULL){
			freePolynomial_AA (tmp2);
		}
		freePolynomial_AA (tmp);
	}


	// fprintf (stderr, "[SMQP-Ali] qb := ");
	// printPoly_AA_unpk(stderr, qb, ch, VSd->r->nvar);
	// fprintf (stderr, "\n");

	if (q != NULL)
		freePolynomial_AA (q);
	
	freePolynomial_AA (csdm);
	freePolynomial_AA (ssd);

	if ((d-e)%2) {
		*hb = hhb;
		*qb = qqb; 
		return Sem;
	}

	// fprintf(stderr, "[SMQP-Ali] neg... \n");

	if (hhb == NULL || hhb->size == 0) {
		*hb = NULL;
	} else {
		negatePolynomial_AA (hhb);
		*hb = hhb;
	}

	if (qqb == NULL || qqb->size == 0) {
		*qb = NULL;
	} else {
		negatePolynomial_AA (qqb);
		*qb = qqb;
	}

	if (Sem == NULL || Sem->size == 0){
		return NULL;
	}

	negatePolynomial_AA (Sem);
	return Sem;

}

void exDucosSubresultantChain_rev_AA (AltArr_t* P, AltArr_t* Q, exgcds_t** SC, int* len, int type)
{
    // char* ch[11] = {"x", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j"}; // TEST 

	// fprintf (stderr, "[SMQP-Ali] P := ");
	// printPoly_AA_unpk(stderr, P, ch, P->nvar);
	// fprintf (stderr, "\n");

	// fprintf (stderr, "[SMQP-Ali] Q := ");
	// printPoly_AA_unpk(stderr, Q, ch, P->nvar);
	// fprintf (stderr, "\n");

    if (P == NULL || P->size == 0){
		exgcds_t* SS0 = (exgcds_t*) malloc (sizeof(exgcds_t));
		SS0->r = NULL;
		SS0->a = NULL;
		SS0->b = NULL;
		SS0->next = NULL;
		*SC = SS0;
		*len = 0;
		return;
    }
    
    if (Q == NULL || Q->size == 0){
	    mpq_t pone; 
	    mpq_init (pone);
	    mpq_set_si (pone, 1l, 1lu);

		exgcds_t* SS1 = (exgcds_t*) malloc (sizeof(exgcds_t));
		SS1->r = deepCopyPolynomial_AA (P);
		SS1->a = makeConstPolynomial_AA (1, P->nvar, pone);
		SS1->b = NULL;
		SS1->next = NULL;
		*SC = SS1;
		*len = 1;
		mpq_clear (pone);
		return;
    }
    if (P->nvar != Q->nvar){
		fprintf (stderr, "SMQP Error: In Extended DucosSubresultantChain_rev, P->nvar(%d) != Q->nvar(%d).\n",P->nvar, Q->nvar);
		exit (EXIT_FAILURE);
	}

  
    //   // adding corner cases for supprt algebraic extensions
    // // TODO: revise for efficiency
    // if (degP == 0 && degQ == 0) {
    // 	mpq_t one;
    // 	mpq_init (one);
    // 	mpq_set_si (one, 1l, 1ul);
    // 	AltArrs_t* SS2 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    // 	SS2->poly = makeConstPolynomial_AA (1, P->nvar, one);
    // 	SS2->next = NULL;
    // 	*SC = SS2;
    // 	*len = 1;
    // 	return;
    // } else if (degP == 0) {
    // 	AltArrs_t* SS3 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    // 	SS3->poly = exponentiatePoly_AA (P, degQ, P->nvar);
    // 	SS3->next = NULL;
    // 	*SC = SS3;
    // 	*len = 1;
    // 	return;
    // } else if (degQ == 0) {
    // 	AltArrs_t* SS4 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    // 	SS4->poly = exponentiatePoly_AA (Q, degP, Q->nvar);
    // 	SS4->next = NULL;
    // 	*SC = SS4;
    // 	*len = 1;
    // 	return;
    // }


    degree_t degP = mainLeadingDegree_AA(P);
    degree_t degQ = mainLeadingDegree_AA(Q);
    
  
    if (degP == 0 || degQ == 0){
        fprintf (stderr, "SMQP: error, Input polynomials to exteneded Subresultant Chain must have positive degree.\n");
		exit (EXIT_FAILURE);
    }
  
    int lazy = 0; // lazy option in pseudoDivide algorithm
    register int nvar = P->nvar;
    
    int tmpE = 0;
    AltArr_t* tmpQ = NULL;
    AltArr_t* tmpH = NULL;
    AltArr_t* tmpB = NULL;
    AltArr_t* s = NULL;
    AltArr_t* p = NULL;
    AltArr_t* q = NULL;
    AltArr_t* A = NULL;
    AltArr_t* B = NULL;
    AltArr_t* C = NULL;

    exgcds_t* SS = NULL;
    exgcds_t* tail = NULL;
    exgcds_t* tailm = NULL;
    degree_t delta = 0;
    int size = 0;
  
    if (degP >= degQ){
		p = P;
		q = Q;
    } else {
		p = Q;
		q = P;
	    
	    degP = mainLeadingDegree_AA(p);
    	degQ = mainLeadingDegree_AA(q);
    }
    
    s = mainLeadingCoefficient_AA (q);
    s = exponentiatePoly_AA (s, degP - degQ, nvar);
    A = deepCopyPolynomial_AA (q);
   

    // fprintf(stderr, "degP: %d, degQ: %d\n", degP, degQ); // TEST
   

    // B =
    AltArr_t* negQ = deepCopyPolynomial_AA (q);
    negatePolynomial_AA (negQ);
    
    pesudoDivideAtIdx_AA (0, p, negQ, &tmpQ, &B, &tmpE, &tmpH, nvar, lazy);
    freePolynomial_AA (negQ);
    // freePolynomial_AA (tmpQ);
    // freePolynomial_AA (tmpH);
    

	// fprintf (stderr, "[SMQP-Ali] B := ");
	// printPoly_AA_unpk(stderr, B, ch, nvar);
	// fprintf (stderr, "\n");

    // gmp_fprintf(stderr, "rem: %Qd*x^%lld\n", B->elems->coef, B->elems->degs); // TEST

    
    mpq_t one;
    mpq_init (one);
    mpq_set_si (one, 1l, 1lu);

    exgcds_t* newPoly = (exgcds_t*) malloc (sizeof(exgcds_t));
    newPoly->r = deepCopyPolynomial_AA (p);
    newPoly->a = makeConstPolynomial_AA (1, P->nvar, one);
    newPoly->b = NULL;
    newPoly->next = NULL;

    SS = newPoly;
    tail = newPoly;
    size++;
    
    newPoly = (exgcds_t*) malloc (sizeof(exgcds_t));
    if (B == NULL || B->size == 0){
		newPoly->r = A;
		
    } else{
		newPoly->r = deepCopyPolynomial_AA (A);
    }
    newPoly->a = NULL;
	newPoly->b = makeConstPolynomial_AA (1, A->nvar, one);
    newPoly->next = NULL;

	exgcds_t* VA = newPoly;    
    exgcds_t* VB = NULL;

    tail->next = newPoly;
    tail = newPoly;
    size++;
    
    register degree_t d = 0;
    register degree_t e = 0;

    while (1){

    	if (B == NULL || B->size == 0){
    		break;
    	}

		d = mainLeadingDegree_AA(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
		e = mainLeadingDegree_AA(B);//(B->elems[0].degs & mvarMask) >> mvarDegOffset;

		newPoly = (exgcds_t*) malloc (sizeof(exgcds_t));
		newPoly->r = deepCopyPolynomial_AA (B);
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
		
			C = exLazardOpt_AA (tailm->r, tail, s, &tmpH, &tmpQ); 
			newPoly = (exgcds_t*) malloc (sizeof(exgcds_t)); 
			newPoly->r = deepCopyPolynomial_AA (C);
			newPoly->a = tmpH;
			newPoly->b = tmpQ;
			newPoly->next = NULL;

			tail->next = newPoly;

			tail = newPoly;
			size++;

		} else {
			C = deepCopyPolynomial_AA (B);
		}

		// fprintf (stderr, "[SMQP-Ali] C := ");
		// printPoly_AA_unpk(stderr, C, ch, C->nvar);
		// fprintf (stderr, "\n");

		if (e == 0){ // B.degree() == 0
			break;
		}
		
		tmpH = NULL;
		tmpQ = NULL;

		// fprintf (stderr, "[SMQP-Ali] VA := ");
		// printPoly_AA_unpk(stderr, VA->r, ch, VA->r->nvar);
		// fprintf (stderr, "\n");

		// fprintf (stderr, "[SMQP-Ali] VB := ");
		// printPoly_AA_unpk(stderr, VB->r, ch, VB->r->nvar);
		// fprintf (stderr, "\n");

		tmpB = exDucosOpt_AA (VA, VB, C, s, &tmpH, &tmpQ);

		// fprintf (stderr, "[SMQP-Ali] tmpB := ");
		// printPoly_AA_unpk(stderr, tmpB, ch, tmpB->nvar);
		// fprintf (stderr, "\n");

		freePolynomial_AA (B);
		freePolynomial_AA (A);
		B = tmpB;
		A = C; // deepCopyPolynomial_AA (C); // TODO: update DucosSubresultant_rev ALgorithm, too. (*)
		
		VA = newPoly;
		
		freePolynomial_AA (s);
		s = mainLeadingCoefficient_AA (A);	
		
		// freePolynomial_AA (C);
	}


	// if (B != NULL && B->size != 0) {
	// 	freePolynomial_AA (B);
	// }

	// if (C != NULL && C->size != 0) {
	// 	freePolynomial_AA (C);
	// }

	// if (A != NULL && A->size != 0) {
	// 	freePolynomial_AA (A);
	// }

	// freePolynomial_AA (s);

    // if subresultant is 0, add it to SS:
    if (tail->r != NULL && tail->r->size != 0 && mainLeadingDegree_AA(tail->r) > 0) {	
		newPoly = (exgcds_t*) malloc (sizeof (exgcds_t));
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

void exDucosSubresultantChain_AA (AltArr_t* P, AltArr_t* Q, exgcds_t** SC, int* len)
{
    exgcds_t* invSC;
    int size = 0;
    exDucosSubresultantChain_rev_AA (P, Q, &invSC, &size, 1);

    // reverse the subresultantchain:
    if (size > 1){
		exgcds_t* curr = invSC;
		exgcds_t* prev = NULL;
		exgcds_t* next = NULL;
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

AltArr_t* exDucosResultant_AA (AltArr_t* P, AltArr_t* Q, AltArr_t** a, AltArr_t** b)
{
    exgcds_t* SC;
    int size = 0;
    exDucosSubresultantChain_AA (P, Q, &SC, &size);
    
    if (size){
    	*a = deepCopyPolynomial_AA (SC->a);
    	*b = deepCopyPolynomial_AA (SC->b);
    	AltArr_t* r = deepCopyPolynomial_AA (SC->r);
    	
    	freeExgcds_AA (SC);
		
		return r;
    }

    freeExgcds_AA (SC);

    *a = NULL;
    *b = NULL;
    return NULL;
}


////////////////////////////
/// Mod Subresultant Chain 
////////////////////////////

// Note: alpha shouldn't be zero! (see the modDucosSubresultantChain_rev)
void modDucosSubresultantChain_rev_AA (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, AltArrs_t** SC, int* len, int type, int isDeepMod)
{
    if (P == NULL || P->size == 0){
		AltArrs_t* SS0 = (AltArrs_t*) malloc (sizeof(AltArrs_t));
		SS0->poly = NULL;
		SS0->next = NULL;
		*SC = SS0;
		*len = 0;
		return;
    }
    if (Q == NULL || Q->size == 0){

    	AltArr_t* mtmp_qq = NULL;
    	AltArr_t* mtmp_rr = NULL;
    	dividePolynomials_AA (P, alpha, &mtmp_qq, &mtmp_rr, P->nvar);
    	freePolynomial_AA (mtmp_qq);

		AltArrs_t* SS1 = (AltArrs_t*) malloc (sizeof(AltArrs_t));
		SS1->poly = mtmp_rr; // deepCopyPolynomial_AA (P);
		SS1->next = NULL;
		*SC = SS1;
		*len = 1;
		return;
    }

    if (P->nvar != Q->nvar || alpha->nvar != Q->nvar){
		fprintf (stderr, "SMQP Error: In DucosSubresultantChain_rev, P->nvar(%d) != Q->nvar(%d) OR alpha->nvar(%d) != Q->nvar(%d).\n",P->nvar, Q->nvar, alpha->nvar, Q->nvar);
		exit (EXIT_FAILURE);
    }


  	degree_t degP = mainLeadingDegree_AA(P);
    degree_t degQ = mainLeadingDegree_AA(Q);
  
  	// commented to support constant polynomials in algebraic extension
  //   if ( degP == 0 || degQ == 0){
  //       fprintf (stderr, "SMQP Error: Input polynomials to subresultantChain must have positive degree.\n");
		// exit (EXIT_FAILURE);
  //   }

    // adding corner cases for supprt algebraic extensions
    // TODO: revise for efficiency
    if (degP == 0 && degQ == 0) {
    	mpq_t one;
    	mpq_init (one);
    	mpq_set_si (one, 1l, 1ul);
    	AltArrs_t* SS2 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    	SS2->poly = makeConstPolynomial_AA (1, P->nvar, one);
    	mpq_clear(one);
    	SS2->next = NULL;
    	*SC = SS2;
    	*len = 1;
    	return;
    } else if (degP == 0) {
    	AltArrs_t* SS3 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    	SS3->poly = exponentiatePoly_AA (P, degQ, P->nvar);
    	SS3->next = NULL;
    	*SC = SS3;
    	*len = 1;
    	return;
    } else if (degQ == 0) {
    	AltArrs_t* SS4 = (AltArrs_t*) malloc (sizeof (AltArrs_t));
    	SS4->poly = exponentiatePoly_AA (Q, degP, Q->nvar);
    	SS4->next = NULL;
    	*SC = SS4;
    	*len = 1;
    	return;
    }



    if (mainVariable_AA(P) >= mainVariable_AA(alpha) || mainVariable_AA(Q) >= mainVariable_AA(alpha) ) {
    	fprintf(stderr, "SMQP Error: The main variable of P and Q in alpha must be one.\n");
    }
  
    int lazy = 0; // lazy option in pseudoDivide algorithm
    register int nvar = P->nvar;
    
    int tmpE = 0;
    AltArr_t* tmpq = NULL;
    AltArr_t* tmpH = NULL;
    AltArr_t* tmpB = NULL;
    AltArr_t* s = NULL;
    AltArr_t* p = NULL;
    AltArr_t* q = NULL;
    AltArr_t* A = NULL;
    AltArr_t* B = NULL;
    AltArr_t* C = NULL;

    AltArrs_t* SS = NULL;
    AltArrs_t* tail = NULL;
    AltArr_t* tailm = NULL;
    degree_t delta = 0;
    int size = 0;
      
    if (degP >= degQ){
		p = P;
		q = Q;
    } else {
		p = Q;
		q = P;
    
    	degP = mainLeadingDegree_AA(p);
    	degQ = mainLeadingDegree_AA(q);
    }
    
    s = mainLeadingCoefficient_AA (q);
    s = exponentiatePoly_AA (s, degP - degQ, nvar);
    A = deepCopyPolynomial_AA (q);
    

    // fprintf(stderr, "degP: %d, degQ: %d\n", degP, degQ); // TEST

    
    // B =
    AltArr_t* negQ = deepCopyPolynomial_AA (q);
    negatePolynomial_AA (negQ);
    
    // deepCopyPolynomial_AA(p), deepCopyPolynomial_AA(negQ)
    pesudoDivideAtIdx_AA (0, p, negQ, &tmpq, &B, &tmpE, &tmpH, nvar, lazy);
    freePolynomial_AA (tmpq);
    freePolynomial_AA (tmpH);
    freePolynomial_AA (negQ);


    // gmp_fprintf(stderr, "rem: %Qd*x^%lld\n", B->elems->coef, B->elems->degs); // TEST

    // computing mod alpha
    AltArr_t* mtmp_r = NULL;
    AltArr_t* mtmp_q = NULL;
    dividePolynomials_AA (p, alpha, &mtmp_q, &mtmp_r, nvar);
    freePolynomial_AA (mtmp_q);

    AltArrs_t* newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
    newPoly->poly = mtmp_r; // deepCopyPolynomial_AA (p);
    newPoly->next = NULL;
    
    SS = newPoly;
    tail = newPoly;
    size++;
    
    // computing mod alpha 
    mtmp_r = NULL;
    mtmp_q = NULL;
    dividePolynomials_AA (A, alpha, &mtmp_q, &mtmp_r, nvar);
    freePolynomial_AA (mtmp_q);

    newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
  //   if (B == NULL || B->size == 0){
		// newPoly->poly = A;
  //   } else{
		// newPoly->poly = deepCopyPolynomial_AA (A);
  //   }
    newPoly->poly = mtmp_r;
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

		// computing mod alpha
		mtmp_r = NULL;
		mtmp_q = NULL;
		dividePolynomials_AA (B, alpha, &mtmp_q, &mtmp_r, nvar);
		freePolynomial_AA (mtmp_q);
		
		newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
		newPoly->poly = mtmp_r; // deepCopyPolynomial_AA (B);
		newPoly->next = NULL;

		tailm = tail->poly;
		tail->next = newPoly;
		tail = newPoly;
		size++;

    	if (mtmp_r == NULL || mtmp_r->size == 0){
    		break;
    	}

		if (isDeepMod) {
			freePolynomial_AA (B);
			B = deepCopyPolynomial_AA (mtmp_r);
		}


		d = mainLeadingDegree_AA(A);//(A->elems[0].degs & mvarMask) >> mvarDegOffset;
		e = mainLeadingDegree_AA(B);//(B->elems[0].degs & mvarMask) >> mvarDegOffset;

		delta = d - e;
		if (delta > 1){

			C = LazardOpt (tailm, tail->poly, s); 

			// computing mod alpha
			mtmp_r = NULL;
			mtmp_q = NULL;
			dividePolynomials_AA (C, alpha, &mtmp_q, &mtmp_r, nvar);
			freePolynomial_AA (mtmp_q);

			newPoly = (AltArrs_t*) malloc (sizeof(AltArrs_t));
			newPoly->poly = mtmp_r;
			newPoly->next = NULL;

			tail->next = newPoly;
			tail = newPoly;
			size++;

			if (isDeepMod) {
				freePolynomial_AA (C);
				C = deepCopyPolynomial_AA (mtmp_r);
			}

		} else {
			C = deepCopyPolynomial_AA (B);
		}

		if (e == 0){ // B.degree() == 0
			break;
		}
		
		if (type == 0){
			tmpB = DucosOpt (A, B, C, s);
		} else {
			tmpB = CFDucosOpt (A, B, C, s);
		}

		freePolynomial_AA (B);
		freePolynomial_AA (A);
		B = tmpB;
		A = deepCopyPolynomial_AA (C);
		freePolynomial_AA (s);
		s = mainLeadingCoefficient_AA (A);	
		freePolynomial_AA (C);
	}
	

	// if (B != NULL && B->size != 0) {
	// 	freePolynomial_AA (B);
	// }

	// if (C != NULL && C->size != 0) {
	// 	freePolynomial_AA (C);
	// }

	// if (A != NULL && A->size != 0) {
	// 	freePolynomial_AA (A);
	// }

	// freePolynomial_AA (s);


    // if subresultant is 0, add it to SS:
    if (tail->poly != NULL && tail->poly->size != 0 && mainLeadingDegree_AA(tail->poly) > 0) {	
		newPoly = (AltArrs_t*) malloc (sizeof (AltArrs_t));
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

void modDucosSubresultantChain_AA (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, AltArrs_t** SC, int* len, int isDeepMod)
{
    AltArrs_t* invSC;
    int size = 0;

    if (alpha == NULL || alpha->size == 0) {
    	DucosSubresultantChain (P, Q, SC, len);
    	return;
    }

	// TODO: isDeepMod needs to be revised for getting more performance!
	modDucosSubresultantChain_rev_AA (P, Q, alpha, &invSC, &size, 1, 0);

    // reverse the subresultantchain:
    if (size > 1){
	AltArrs_t* curr = invSC;
	AltArrs_t* prev = NULL;
	AltArrs_t* next = NULL;
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

AltArr_t* modDucosResultant_AA  (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, int isDeepMod)
{
    AltArrs_t* SC;
    int size = 0;
    modDucosSubresultantChain_AA (P, Q, alpha, &SC, &size, isDeepMod);
    if (size){
    	AltArr_t* res = deepCopyPolynomial_AA (SC->poly);
    	freeAltArrs (SC);
		return res;
    }

    freeAltArrs (SC);

    return NULL;
}

AltArr_t* modLastNonZeroChain_AA (AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, int isDeepMod)
{
    if (P == NULL || P->size == 0){
    	return deepCopyPolynomial_AA (Q);
    }
    if (Q == NULL || Q->size == 0){
    	return deepCopyPolynomial_AA (P);
    }
    if (P->nvar != Q->nvar){
    	fprintf (stderr, "SMQP Error: In lastNonZeroChain, P->nvar(=%d) != Q->nvar(=%d).\n",
    		 P->nvar, Q->nvar);
    	exit (EXIT_FAILURE);
    }

    int nvar = P->nvar;
    
    mpq_t one;
    mpq_init(one);
    mpq_set_si (one, 1l, 1l);
    
    degree_t degP = mainLeadingDegree_AA(P);//(P->elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t degQ = mainLeadingDegree_AA(Q);//(Q->elems[0].degs & mvarMask) >> mvarDegOffset;

  //   if (nvar == 1){
		// return univariateGCD_AA (P, Q);
  //   }
    
    if (degP == 0 || degQ == 0){
    	return makeConstPolynomial_AA (1, nvar, one);
    }
    
    AltArrs_t* SC;
    int size = 0, idx = 0;
    // TODO: isDeepMod has no impact!
    // modDucosSubresultantChain_AA (P, Q, alpha, &SC, &size, isDeepMod);
    DucosSubresultantChain (P, Q, &SC, &size);

    AltArrs_t* cur = SC;
    
    while (idx < size) {
		if (cur->poly != NULL && cur->poly->size != 0){
			AltArr_t* tmp_r = NULL;
			AltArr_t* tmp_q = NULL;
			dividePolynomials_AA (cur->poly, alpha, &tmp_q, &tmp_r, nvar);
			freePolynomial_AA (tmp_q);
			if (tmp_r != NULL && tmp_r->size != 0) {
				freeAltArrs (SC);
				return tmp_r;
			} else {
				tmp_r = NULL;
			}
		}
		cur = cur->next;
		++idx;
    }

    freeAltArrs  (SC);    
    return NULL;
}

// TODO: temporary function ... needed more conditions for integrating into SMQP.
AltArr_t* modKthSubresultantChain_AA (degree_t k, AltArr_t* P, AltArr_t* Q, AltArr_t* alpha, int isDeepMod) {

	degree_t degP = mainLeadingDegree_AA (P);
	degree_t degQ = mainLeadingDegree_AA (Q);

	if (k > degP || k < -1 ) {
		return NULL;
	} else if (k == degP) {
		return deepCopyPolynomial_AA (P);
	} else if (k == degQ) {
		return deepCopyPolynomial_AA (Q);
	}

	AltArrs_t* SC;
	int size = 0;

	modDucosSubresultantChain_AA (P, Q, alpha, &SC, &size, isDeepMod);

	AltArrs_t* cur = SC;
	AltArr_t* res = NULL;
	int idx = 0;

  	if (!(k+1)) {
	    while (idx < size) {
			if (cur->poly != NULL && cur->poly->size != 0) {
			    break;
			}
			cur = cur->next;
			++idx;
	    }
	    res = deepCopyPolynomial_AA (cur->poly);
	    freeAltArrs  (SC);
	    return res;
  	} else {

	  	while (mainLeadingDegree_AA (cur->poly) < k) {
	  		cur = cur->next;
	  	}

	  	while (mainLeadingDegree_AA (cur->poly) == k) {
	  		res = cur->poly;
	  		cur = cur->next;
	  	}

	  	if (res == NULL || res->size == 0) {
	  		freeAltArrs (SC);
	  		return NULL;
	  	} else {
	  		res = deepCopyPolynomial_AA (res);
	  		freeAltArrs (SC);
	  		return res;
	  	}
  	}
}


///////////////////////////////////
///// Triangular Set Normalization
///////////////////////////////////

AltArr_t* normalizePolynomial_AA (AltArr_t* P, AltArr_t** T, AltArr_t** A, int s, int nvar)
{
	int v = mainVariable_AA (P);

	if (v == -1) {
		*A = NULL;
		return NULL;
	}

	if (isConstant_AA(P) > 0) {
		mpq_t aone;
		mpq_init (aone);
		mpq_set_si (aone, 1l, 1lu);
		*A = makeConstPolynomial_AA (1, P->nvar, aone);

		return deepCopyPolynomial_AA (P);
	}

	int isRes = -1;

	for (int i = 0; i < s; i++) {
		if (v == mainVariable_AA (T[i])) {
			isRes = i;
			break;
		}
	} 
	
	if (isRes==-1) {
		degree_t d = mainDegree_AA (P);
		AltArr_t* h = mainLeadingCoefficient_AA (P);
		AltArr_t* hv = deepCopyPolynomial_AA (h); 
		multiplyPolynomialAtIdxByXn_AA_inp (h, v, d, nvar); // h = h*v^d
		AltArr_t* t = subPolynomials_AA (P, hv, nvar); // t = P - h*v^mdeg 
		
		AltArr_t* q = NULL;
		AltArr_t* r = normalizePolynomial_AA (h, T, &q, s, nvar);

		*A = q;
		
		if (r == NULL || r->size == 0) {
			r = NULL;
		} else {
			multiplyPolynomialAtIdxByXn_AA_inp (r, v, d, nvar); // r = r*v^d 
		}

		if (t == NULL || q == NULL || 
			t->size == 0 || q->size == 0) {
			return r;
		}

		q = multiplyPolynomials_AA_inp (q, t, nvar);
		AltArr_t* nf_r = onlyNormalForm_AA (q, T, s, nvar); 

		r = addPolynomials_AA_inp (r, nf_r, nvar);

		freePolynomial_AA (hv);
		freePolynomial_AA (nf_r);
		freePolynomial_AA (h);
		freePolynomial_AA (t);

		return r;
	} else {
		AltArr_t* a = NULL;
		AltArr_t* b = NULL;
		AltArr_t* exr = exDucosResultant_AA (P, T[isRes], &a, &b);

		freePolynomial_AA (b);

		AltArr_t* nf_a   = onlyNormalForm_AA (a, T, s, nvar);
		AltArr_t* nf_exr = onlyNormalForm_AA (exr, T, s, nvar);

		freePolynomial_AA (a);
		freePolynomial_AA (exr);

		*A = nf_a;
		return nf_exr;
	}
}


void normalizeRegularChainDim0 (AltArr_t** T, AltArr_t** N, int s, int nvar) 
{

	if (s < 1 || T == NULL) {
		*N = NULL;
		return;
	}


	AltArr_t** NN = (AltArr_t**) malloc (sizeof(AltArr_t*)*s);

	// If T = {t0}, then the init(T[0]) must be a constant,
	AltArr_t* hi = mainLeadingCoefficient_AA (T[0]);
	AltArr_t* a0 = NULL;
	AltArr_t* r0 = NULL;
	AltArr_t* q0 = NULL;
	// AltArr_t* t0 = NULL;

	if (hi == NULL || hi->size == 0) {
		NN[0] = NULL;
	} else {
		dividePolynomials_AA (T[0], hi, &q0, &r0, nvar);

		if (r0 != NULL && r0->size != 0) {
			freePolynomial_AA (r0); r0 = NULL;
		} 	
		if (hi != NULL && hi->size != 0) {
			freePolynomial_AA (hi); 
		}

		NN[0] = q0;
		q0 = NULL;
	}

	for (int i = 1; i < s; i++) {
		hi = mainLeadingCoefficient_AA (T[i]);
		r0 = normalizePolynomial_AA (hi, T, &a0, i, nvar);
		
		if (hi != NULL && hi->size != 0)	
			freePolynomial_AA (hi);

		if (a0 == NULL || a0->size == 0) {
			q0 = NULL;
		} else {
			dividePolynomials_AA (a0, r0, &q0, &hi, nvar);

			if (hi != NULL && hi->size != 0)	
				freePolynomial_AA (hi);
			if (r0 != NULL && r0->size != 0)	
				freePolynomial_AA (r0);
			
			freePolynomial_AA (a0);
		}

		if (q0 == NULL || q0->size == 0){
			NN[i] = NULL;
		} else  {
			NN[i] = multiplyPolynomials_AA (T[i], q0, nvar);
			freePolynomial_AA (q0);
		}

	}

}

//////////////
// GCD 
//////////////

AltArr_t* integralGCD_AA_polyOut (AltArr_t* P, AltArr_t* Q)
{
    
    if (P->size != 1 || Q->size != 1){
	fprintf (stderr, "SMQP Error, In integralGCD4_AA the input polynomials must have single term." );
	exit (EXIT_FAILURE);
    }

    mpq_t P0;
    mpq_init (P0);

    mpq_t Q0;
    mpq_init (Q0);

    mpq_t gcd_q;
    mpq_init (gcd_q);
    
    /* mpq_gcd (gcd_q, P->elems[0].coef, Q->elems[0].coef); */
    mpz_mul (mpq_numref(P0), mpq_numref(P->elems[0].coef), mpq_denref(Q->elems[0].coef));
    mpz_mul (mpq_numref(Q0), mpq_numref(Q->elems[0].coef), mpq_denref(P->elems[0].coef));
    mpz_mul (mpq_denref(gcd_q), mpq_denref(Q->elems[0].coef), mpq_denref(P->elems[0].coef));
    mpz_gcd (mpq_numref(gcd_q), mpq_numref(P0), mpq_numref(Q0));
    mpq_canonicalize (gcd_q); 

    mpq_clear (P0);
    mpq_clear (Q0);
    return makeConstPolynomial_AA (1, P->nvar, gcd_q);
}

AltArr_t* gcd_AA_Q (AltArr_t* P, mpq_t c)
{
    if (P == NULL) {
	return NULL;
    }

    mpq_t ret;
    mpq_init (ret);
    integralContent_AA (P, ret);
    
    mpq_t gcd_z;
    mpq_init (gcd_z);
    /* mpq_gcd (gcd_z, ret, c); */

    mpq_t P0;
    mpq_init (P0);

    mpq_t Q0;
    mpq_init (Q0);
    
    /* mpq_gcd (gcd_q, P->elems[0].coef, Q->elems[0].coef); */
    mpz_mul (mpq_numref(P0), mpq_numref(ret), mpq_denref(c));
    mpz_mul (mpq_numref(Q0), mpq_numref(c), mpq_denref(ret));
    mpz_mul (mpq_denref(gcd_z), mpq_denref(ret), mpq_denref(c));
    mpz_gcd (mpq_numref(gcd_z), mpq_numref(P0), mpq_numref(Q0));
    mpq_canonicalize (gcd_z); 
    mpq_clear (P0);
    mpq_clear (Q0);
    mpq_clear (ret);
    
    return makeConstPolynomial_AA (1, P->nvar, gcd_z);
}

AltArr_t* gcd_AA (AltArr_t* P, AltArr_t* Q)
{
    int mvarP = leadingVariable_AA (P);
    int mvarQ = leadingVariable_AA (Q);
    
    if (mvarP < 0){
	if (mvarP == -2){
	    return deepCopyPolynomial_AA (Q);
	}
	if (mvarQ == -2) {
	    return deepCopyPolynomial_AA (P);
	}
	return gcd_AA_Q (Q, P->elems[0].coef);
    }
    if (mvarQ < 0){
	return gcd_AA_Q (P, Q->elems[0].coef);
    }

    if (P->nvar == 1 && Q->nvar == 1){
	return univariateGCD_AA (P, Q);
    }
    
    if (mvarP < mvarQ){
	return gcd_AA (Q, mainContent_AA (P));
    }
    if (mvarP > mvarQ){
	return gcd_AA (P, mainContent_AA (Q));
    }
    
    // if univariate?(p1) and univariate?(p2) then //TODO:
    // convert p1 and p2 to SUZP and call the SUZP gcd
    
    AltArr_t* contP = NULL;//(AltArr_t*) malloc (sizeof (AltArr_t));
    AltArr_t* contQ = NULL;//(AltArr_t*) malloc (sizeof (AltArr_t));
    AltArr_t* cont = NULL;
    AltArr_t* gcd = NULL;
    int nvar = P->nvar;
    int isShrinked = 0;
    /* int ppPmvar, ppQmvar; */
    if (nvar != Q->nvar){
	fprintf (stderr, "SMQP Error, In gcd_AA, the nvar of input polynomials must be the same.");
	exit (EXIT_FAILURE);
    }
    
    AltArr_t* ppP = mainPrimitiveFactorization_AA (P, &contP); 
    AltArr_t* ppQ = mainPrimitiveFactorization_AA (Q, &contQ); 
    
    cont = gcd_AA (contP, contQ);

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
		
		shrinkAndReorderVars_AA (ppP, varmap, nvar);
		shrinkAndReorderVars_AA (ppQ, varmap, nvar);
		isShrinked = 1;
	}

    /* for (int i = 0; i < mvarP; ++i){ */
	/* shrinkNumVarsAtIdx_AA (ppP, 0); */
	/* shrinkNumVarsAtIdx_AA (ppQ, 0); */
	/* isShrinked++; */
    /* } */
    
    AltArr_t* lnzch = lastNonZeroChain_AA (ppP, ppQ);
    
    if (isShrinked && lnzch != NULL && lnzch->size != 0){
	expandNumVarsLeft_AA (lnzch, nvar);
    }

    if (leadingVariable_AA (lnzch) != mvarP){
	mpq_t one;
	mpq_init (one);
	mpq_set_si (one, 1l, 1l);
	gcd = makeConstPolynomial_AA (1, nvar, one);
    } else {
	gcd = mainPrimitivePart_AA (lnzch, mvarP);
    }
        
    gcd = multiplyPolynomials_AA_inp (gcd, cont, nvar);

    freePolynomial_AA (contP);
    freePolynomial_AA (contQ);
    freePolynomial_AA (ppP);
    freePolynomial_AA (ppQ);
    freePolynomial_AA (cont);
    
    return gcd;
}

AltArr_t* mainPrimitiveFactorization_AA (AltArr_t* P, AltArr_t** cont)
{
    int mvarP = leadingVariable_AA (P);
    AltArr_t* sP;
    int isShrinked = 0;
    mpq_t one;
    mpq_init (one);
    mpq_set_si (one, 1l, 1l);
    
    if (mvarP < 0){
	if (mvarP == -1){
	    /* if (*cont == NULL || (*cont)->size == 0){ */
	    mpq_t zz;
	    mpq_init (zz);	    
	    mpq_set(zz, P->elems[0].coef);
	    *cont =  makeConstPolynomial_AA (1, P->nvar, zz);
	    /* } */
	}
	
	return makeConstPolynomial_AA (1, P->nvar, one);
    }
    
    if (mvarP){
	sP = deepCopyPolynomial_AA (P);
	/* for (int i = 0; i < mvarP; ++i){ */
	/*     shrinkNumVarsAtIdx_AA (sP, 0); */
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
		
	shrinkAndReorderVars_AA (sP, varmap, P->nvar);
	isShrinked = 1;
	
    } else {
	sP = deepCopyPolynomial_AA (P);
    }

    AltArr_t* tmp = NULL;
    AltArr_t* tmpi = NULL; 
    RecArr_t* recP = convertToRecursiveArray  (sP);
    RecArrElem_t* elems = recP->elems;
    int rSize = recP->size;
    
    int isOne = 0, idx = 1;
    AltArr_t* coefGCD = deepCopyPolynomial_AA (convertFromAAElemToAA
					       (elems[0].coef,  elems[0].coefSize,
						sP->nvar, recP->unpacked));
    
    while (idx < rSize){
	if (coefGCD != NULL && coefGCD->size == 1 &&
	    leadingVariable_AA (coefGCD) == -1 &&
	    !mpq_cmp(coefGCD->elems->coef, one)){
	    isOne = 1;
	    break;
	}
	
	tmpi = deepCopyPolynomial_AA (convertFromAAElemToAA (elems[idx].coef,
							     elems[idx].coefSize,
							     sP->nvar, recP->unpacked));
	tmp = gcd_AA (coefGCD, tmpi);
	freePolynomial_AA (coefGCD);
	freePolynomial_AA (tmpi);
	coefGCD = tmp;
	++idx;
    }
    
	if (isShrinked && coefGCD != NULL && coefGCD->size != 0){
	    expandNumVarsLeft_AA (coefGCD, P->nvar);
	}
	if (coefGCD != NULL && coefGCD->size != 0){
	    *cont = coefGCD;
	} else {
	    fprintf (stderr, "SMQP Error, in mainPrimitiveFactorization_AA, coefGCD is NULL!\n");
	    exit (EXIT_FAILURE);
	}
    mpq_clear (one);
    
    if (!isOne){
	AltArr_t* tmpR = NULL;
	AltArr_t* tmpQ = NULL;
		
	dividePolynomials_AA (P, coefGCD, &tmpQ, &tmpR, P->nvar);
	freePolynomial_AA (tmpR);
		
	return tmpQ;
    }
    
    return deepCopyPolynomial_AA (P);
}

AltArr_t* mainPrimitivePart_AA (AltArr_t* P, int mvar)
{
    int mvarP = mvar;
    AltArr_t* sP;
    int isShrinked = 0;
    mpq_t one;
    mpq_init (one);
    mpq_set_si (one, 1l, 1l);
    
    if (mvarP < 0){
	return makeConstPolynomial_AA (1, P->nvar, one);
    }
    
    if (mvarP){
	sP = deepCopyPolynomial_AA (P);
	/* for (int i = 0; i < mvarP; ++i){ */
	/*     shrinkNumVarsAtIdx_AA (sP, 0); */
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
		
		shrinkAndReorderVars_AA (sP, varmap, P->nvar);
		isShrinked = 1;
    } else {
	sP = deepCopyPolynomial_AA (P);
    }

    AltArr_t* tmp = NULL;
    AltArr_t* tmpi = NULL; 
    RecArr_t* recP = convertToRecursiveArray  (sP);
    RecArrElem_t* elems = recP->elems;
    int rSize = recP->size;
    
    int isOne = 0, idx = 1;
    AltArr_t* coefGCD = deepCopyPolynomial_AA (convertFromAAElemToAA
						 (elems[0].coef,  elems[0].coefSize,
						  sP->nvar, recP->unpacked));
    
    while (idx < rSize){
	if (coefGCD != NULL && coefGCD->size == 1 &&
	    leadingVariable_AA (coefGCD) == -1 &&
	    !mpq_cmp(coefGCD->elems->coef, one)){
	    isOne = 1;
	    break;
	}
	
	tmpi = deepCopyPolynomial_AA (convertFromAAElemToAA (elems[idx].coef,
								elems[idx].coefSize,
								sP->nvar, recP->unpacked));
	tmp = gcd_AA (coefGCD, tmpi);
	freePolynomial_AA (coefGCD);
	freePolynomial_AA (tmpi);
	coefGCD = tmp;
	++idx;
    }
    
    if (isShrinked && coefGCD != NULL && coefGCD->size != 0){
	expandNumVarsLeft_AA (coefGCD, P->nvar);
    }
    mpq_clear (one);
    
    if (!isOne){
	AltArr_t* tmpR = NULL;
	AltArr_t* tmpQ = NULL;

	dividePolynomials_AA (P, coefGCD, &tmpQ, &tmpR, P->nvar);
	freePolynomial_AA (tmpR);
	
	return tmpQ;
    }
    
    return deepCopyPolynomial_AA (P);
}

AltArr_t* mainContent_AA (AltArr_t* P)
{    	

    int mvarP = leadingVariable_AA (P);
    AltArr_t* sP;
    int isShrinked = 0;
    mpq_t one;
    mpq_init (one);
    mpq_set_si (one, 1l, 1l);
    
    if (mvarP == -2){
	return makeConstPolynomial_AA (1, P->nvar, one);
    }
    if (mvarP == -1){
	return deepCopyPolynomial_AA (P);
    }
        
    if (mvarP){
	sP = deepCopyPolynomial_AA (P);
	/* for (int i = 0; i < mvarP; ++i){ */
	/*     shrinkNumVarsAtIdx_AA (sP, 0); */
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
		
		shrinkAndReorderVars_AA (sP, varmap, P->nvar);
	isShrinked = 1;
	
    } else {
	sP = deepCopyPolynomial_AA (P);
    }

    AltArr_t* tmp = NULL;
    AltArr_t* tmpi = NULL; 
    RecArr_t* recP = convertToRecursiveArray  (sP);
    RecArrElem_t* elems = recP->elems;
    int rSize = recP->size;
        
    int idx = 1;
    AltArr_t* coefGCD = deepCopyPolynomial_AA (convertFromAAElemToAA
						 (elems[0].coef,  elems[0].coefSize,
						  sP->nvar, recP->unpacked));
    
    while (idx < rSize){
		if (coefGCD != NULL && coefGCD->size == 1 &&
		    leadingVariable_AA (coefGCD) == -1 &&
		    !mpq_cmp(coefGCD->elems->coef, one)){
		    break;
		}
		
		tmpi = deepCopyPolynomial_AA (convertFromAAElemToAA (elems[idx].coef,
									elems[idx].coefSize,
									sP->nvar, recP->unpacked));
		
		tmp = gcd_AA (coefGCD, tmpi);
		freePolynomial_AA (coefGCD);
		freePolynomial_AA (tmpi);
		coefGCD = tmp;
		++idx;
    }
    
    if (isShrinked && coefGCD != NULL && coefGCD->size != 0){
	expandNumVarsLeft_AA (coefGCD, P->nvar);
    }
    if (coefGCD != NULL && coefGCD->size != 0) {
	return coefGCD;
    }

    return makeConstPolynomial_AA (1, P->nvar, one);
}         


/******************
 * Square Free 
 *****************/

/**
 * Computes the square free part of a polynomial.
 */
AltArr_t* squareFreePart_AA (AltArr_t* aa, int nvar)
{
	
	/* char* ch[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST */

	/* fprintf (stderr, "[SMZP-Ali] input polynomial in squareFreePart_AA := "); */
	/* printPoly_AA_unpk(stderr, aa, ch, nvar); */
	/* fprintf (stderr, "\n"); */

    if (aa == NULL || aa->size == 0) {
		return NULL;
    }
	
    mpq_t one;
    mpq_init(one);
    mpq_set_si(one, 1l, 1l);
    
    AltArr_t* fact = makeConstPolynomial_AA(1, nvar, one);
	int lvar = leadingVariable_AA (aa);
	if (lvar < 0) {
		return fact;
	}
	
	AltArr_t* primpart = deepCopyPolynomial_AA (aa);
    AltArr_t* cont = NULL;
    AltArr_t* diff = NULL;
    AltArr_t* g = NULL;
    AltArr_t* next = NULL;
    AltArr_t* swappp = NULL;
	
	
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
		
		shrinkAndReorderVars_AA (primpart, varmap, nvar);
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
	/* fprintf (stderr, "mdeg(primpart) = %d\n", mainDegree_AA(primpart)); */

	/* fprintf (stderr, "[SMZP-Ali] primpart after shrink := "); */
	/* printPoly_AA_unpk(stderr, primpart, ch, nnvar); */
	/* fprintf (stderr, "\n"); */
	
    for (int i = 0; i < nnvar; ++i) {
		varmap[0] = i;
		varmap[i] = 0;
		reorderVars_AA (primpart, varmap, nnvar); // reorder primpart w.r.t. i-th variable

		// primpart should have non-zero mvar:
		if (leadingVariable_AA(primpart)) {
			continue;
		}
		
		// swappp is prim. part of poly w.r.t. i-th variable
		// cont is next primpart
		swappp = mainPrimitiveFactorization_AA (primpart, &cont);
	
		freePolynomial_AA (primpart);

		// diff = d(primpart)/d(x_i) 
		diff = derivative_AA (swappp, 0, 1); 
		
		// g = gcd(swappp, D(swappp)) 
		g = gcd_AA (swappp, diff);

		freePolynomial_AA (diff);

		// next = swappp / g 
		exactDividePolynomials_AA (swappp, g, &next, nnvar);
		
		freePolynomial_AA (g);
		freePolynomial_AA (swappp);
		
		varmap[0] = i;
		varmap[i] = 0;
		reorderVars_AA (next, varmap, nnvar);
		
		if (next != NULL && next->size) {
			if (mpq_sgn(next->elems->coef) < 0) {
				negatePolynomial_AA (next);
			}
			
			fact = multiplyPolynomials_AA_inp (fact, next, nnvar); // update fact
			freePolynomial_AA (next);
		}
		
		if (leadingVariable_AA (cont) < 0) {
			break;
		} else {
			primpart = cont;
		}
    }

	if (isShrinked) {
		expandNumVarsLeft_AA (fact,nvar);
	}
	
    return fact;
}
