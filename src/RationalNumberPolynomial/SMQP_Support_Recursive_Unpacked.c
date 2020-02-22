
#include "RationalNumberPolynomial/SMQP_Support_Recursive_Unpacked.h"
#include "RationalNumberPolynomial/SMQP_Support_Recursive-AA.h"
#include "RationalNumberPolynomial/SMQP_Support_Unpacked.h"

RecArr_t* convertToRecursiveArray_unpk(AltArr_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return convertToRecursiveArray(aa);
	}

	AAElem_t* elems = aa->elems;
	degree_t* degs = (degree_t*) aa->elems->degs;
	degree_t mvarDeg =  degs[0];
	RecArr_t* poly = (RecArr_t*) malloc(sizeof(RecArr_t));
	poly->alloc = mvarDeg;
	RecArrElem_t* recElems = (RecArrElem_t*) malloc(sizeof(RecArrElem_t)*(mvarDeg+1)); 
	poly->elems = recElems;

	int nvar = aa->nvar;
	int curIdx = 0, lastSize = 0;
	degree_t curDeg;
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		curDeg = degs[i*nvar];
		
		if (curDeg != mvarDeg) {
			recElems[curIdx].exp = mvarDeg;
			recElems[curIdx].coefSize = i - lastSize;
			recElems[curIdx].coef = elems + lastSize;
			++curIdx;
			lastSize = i;
			mvarDeg = curDeg;
		}

		//throw away exp of mvar for coef.
		degs[i*nvar] = 0;
	}

	//one more to commit after loop
	recElems[curIdx].exp = mvarDeg;
	recElems[curIdx].coefSize = AA_SIZE(aa) - lastSize;
	recElems[curIdx].coef = elems + lastSize;
	++curIdx;

	poly->size = curIdx;
	poly->origAA = aa;
	poly->unpacked = 1;

	return poly;

}

AltArr_t* convertFromRecursiveArray_unpk(RecArr_t* poly, int nvar) {
	if (poly == NULL || poly->size == 0) {
		return NULL;
	}

	if (!poly->unpacked) {
		return convertFromRecursiveArray(poly, nvar);
	}

	AltArr_t* aa = NULL;
	if (poly->origAA != NULL) {
		aa = poly->origAA;
	} else {
		aa = (AltArr_t*) malloc(sizeof(AltArr_t));
	}


	RecArrElem_t* recElems = poly->elems;
	AAElem_t* elems = poly->elems->coef;
	degree_t* degs = (degree_t*) elems->degs;
	aa->elems = elems;
	register int curSize = 0;
	for (int i = 0; i < poly->size; ++i) {
		degrees_t curDeg = recElems[i].exp;
		for (int j = 0; j < recElems[i].coefSize; ++j) {
			degs[(curSize + j)*nvar] = curDeg;
		}

		curSize += recElems[i].coefSize;
	}

	freeRecArray(poly);

	aa->size = curSize;
	aa->alloc = curSize;
	aa->nvar = nvar;
	aa->unpacked = 1;
	return aa;
}

// AltArr_t* recProdHeapGetNextCoef_AA_unpk(RecProdHeap_AA* h, const RecArrElem_t* __restrict__ aElems, const RecArrElem_t* __restrict__ bElems) {
// 	RecProdHeapElem_AA* elems = h->elements;
// 	int nvar = h->nvar;
// 	register int lastB = h->lastB;

// 	AltArr_t* ret = NULL;
// 	register degree_t maxExp = elems->exp;
// 	register degree_t nextExp = elems->exp;

// 	RecProdHeapChain_AA* insertChain = NULL;
// 	RecProdHeapChain_AA* maxElem, * nextMaxElem;
	
// 	AltArr_t aCoef, bCoef;
// 	aCoef.nvar = nvar;
// 	bCoef.nvar = nvar;
// 	aCoef.unpacked = h->unpacked;
// 	bCoef.unpacked = h->unpacked;

// 	//It's possible that not all elemtns of same degree are chained. So much loop like this.
// 	while(nextExp == maxExp) {
// 		maxElem = recProdHeapRemoveMax_AA(h);

// 		//go through the chain now;
// 		while (maxElem != NULL) {
// 			nextMaxElem = maxElem->next;

// 			//Notice here we delay the calculation of the actual product coef.
// 			//This is because a_i will change over the course of the algorithm.
// 			if (ret == NULL) {
// 				aCoef.alloc = aCoef.size = aElems[maxElem->a_i].coefSize;
// 				aCoef.elems = aElems[maxElem->a_i].coef;
// 				bCoef.alloc = bCoef.size = bElems[maxElem->b].coefSize;
// 				bCoef.elems = bElems[maxElem->b].coef;
// 				ret = multiplyPolynomials_AA(&aCoef, &bCoef, nvar);
// 			} else {
// 				aCoef.alloc = aCoef.size = aElems[maxElem->a_i].coefSize;
// 				aCoef.elems = aElems[maxElem->a_i].coef;
// 				bCoef.alloc = bCoef.size = bElems[maxElem->b].coefSize;
// 				bCoef.elems = bElems[maxElem->b].coef;
// 				AltArr_t* prod = multiplyPolynomials_AA(&aCoef, &bCoef, nvar);
// 				ret = addPolynomials_AA_inp(ret, prod, nvar);
// 				freePolynomial_AA(prod);
// 			}
// 			if (maxElem->b != lastB) {
// 				++(maxElem->b);
// 				maxElem->next = insertChain;
// 				insertChain = maxElem;
// 			} else {
// 				maxElem->next = NULL;
// 				recProdHeapFreeChain_AA(maxElem);
// 			}
// 			maxElem = nextMaxElem;
// 		}

// 		nextExp = recProdHeapPeek_AA(h);
// 	}

// 	while(insertChain != NULL) {
// 		nextMaxElem = insertChain->next;
// 		insertChain->next = NULL;
// 		recProdHeapInsert_AA(h, insertChain, aElems[insertChain->a_i].exp + bElems[insertChain->b].exp);
// 		insertChain = nextMaxElem;
// 	}

// 	return ret;
// }

// void recArrayMultiplyByCoef_inp_unpk(RecArrElem_t* aElems, int aSize, AltArr_t* coef, int nvar) {
// 	AltArr_t* tempa = malloc(sizeof(AltArr_t));
// 	tempa->nvar = nvar;
// 	tempa->unpacked = 1;
// 	for (int i = 0; i < aSize; ++i) {
// 		tempa->elems = aElems[i].coef;
// 		tempa->alloc = tempa->size = aElems[i].coefSize;
// 		tempa = multiplyPolynomials_AA_inp(tempa, coef, nvar);

// 		aElems[i].coef = tempa->elems;
// 		aElems[i].coefSize = tempa->size;
// 	}

// 	free(tempa);
// }

// void pesudoDivideOneTerm_RecArray_unpk(RecArr_t* c, RecArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy) {

// 	if (b == NULL) {
// 		fprintf(stderr, "Division by zero! Exiting...\n");
// 		exit(EXIT_FAILURE);
// 	}

// 	if (c == NULL) {
// 		*res_a = NULL;
// 		*res_r = NULL;
// 		return;
// 	}

// 	if (!c->unpacked && !b->unpacked) {
// 		pesudoDivideOneTerm_RecArray(c, b, res_a, res_r, e, hPow, nvar, lazy);
// 		return;
// 	}

// 	int unpackedC = 0, unpackedB = 0;
// 	if (!c->unpacked) {
// 		unpackExponentVectors_RecArray_inp(c, nvar);
// 		unpackedC = 1;
// 	}
// 	if (!b->unpacked) {
// 		unpackExponentVectors_RecArray_inp(b, nvar);
// 		unpackedB = 1;
// 	}

// // TODO ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 	if (b->elems->exp > c->elems->exp) {
// 		RecArr_t* recR = deepCopyRecArrayPolynomial(c, nvar); 
// 		*res_r = convertFromRecursiveArray(recR, nvar);
// 		*res_a = NULL;
// 		if (e != NULL) {
// 			*e = 0;
// 		}
// 		if (hPow != NULL) {
// 			*hPow = makePolynomial_AA(1, nvar);
// 			mpq_init((*hPow)->elems->coef);
// 			mpq_set_si((*hPow)->elems->coef, 1l, 1l);
// 			(*hPow)->elems->degs = 0;
// 			(*hPow)->size = 1;
// 		}
// 		return;
// 	}

// 	RecArrElem_t* __restrict__ k = c->elems;
// 	RecArrElem_t* __restrict__ lenK = k + c->size;

// 	register int maxASize = k->coefSize;
// 	register int maxRSize = k->coefSize;
// 	register int ai = 0;
// 	register int rj = 0;

// 	//so actually, we won't call makePolynomial_AA_unpk
// 	//we do this because then we have an extra exponent array
// 	//we will just steal multerm's exponent vector allocation instead
// 	AltArr_t* a = makePolynomial_AA(maxASize, nvar);
// 	AltArr_t* r = makePolynomial_AA(maxRSize, nvar);

// 	register degree_t bDeg = b->elems->exp;
// 	register degree_t eps = k->exp - bDeg;
// 	AltArr_t kCoef;
// 	kCoef.alloc = kCoef.size = k->coefSize;
// 	kCoef.elems = k->coef;
// 	kCoef.nvar = nvar;
// 	kCoef.unpacked = 1;

// 	int i = 1;
// 	AltArr_t h;
// 	h.alloc = h.size = b->elems->coefSize;
// 	h.elems = b->elems->coef;
// 	h.nvar = nvar;
// 	h.unpacked = 1;

// 	AltArr_t* multerm = deepCopyPolynomial_AA(&kCoef);
// 	AltArr_t* hI = deepCopyPolynomial_AA(&h);

// 	//do first division manually 
// 	//we take all of multerm's allocated memory, including unpacked degs
// 	memcpy(a->elems + ai, multerm->elems, sizeof(AAElem_t)*multerm->size);
// 	//update degs in-place by degree of the mainvar, eps.
// 	degree_t* aDegs = (degree_t*) a->elems->degs;
// 	for (int idx = ai; idx < ai + multerm->size; ++idx) {
// 		aDegs[idx*nvar] = eps;
// 		// a->elems[idx].degs += ((degrees_t)eps << mvarOffset);
// 	}
// 	ai += multerm->size;
// 	++k;
// 	free(multerm); //free only the container struct

// 	while (k != lenK) {
// 		kCoef.alloc = kCoef.size = k->coefSize;
// 		kCoef.elems = k->coef;

// 		eps = k->exp;
// #if PDIVIDE_DIVISBLE_CHECK
// 		// in this case only, when the divisor is one term, divisibility check
// 		// will always pass
// 		multerm = deepCopyPolynomial_AA(&kCoef);
// #else
// 		multerm = multiplyPolynomials_AA(&kCoef, hI, nvar);
// #endif
// 		++k;

// 		if (eps >= bDeg) {
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

// 			eps -= bDeg;

// 			//resize if necessary for insertion to a
// 			if (ai + multerm->size > maxASize) {
// 				a->size = ai;
// 				maxASize += (multerm->size << 1);
// 				resizePolynomial_AA(a, maxASize);
// 			}

// 			memcpy(a->elems + ai, multerm->elems, sizeof(AAElem_t)*multerm->size);
// 			//update degs in-place by degree of the mainvar, eps.
// 			for (int idx = ai; idx < ai + multerm->size; ++idx) {
// 				// a->elems[idx].degs += ((degrees_t)eps << mvarOffset);
// 			}
// 			ai += multerm->size;

// 		} else {
// 			//accumulate remainder.
// 			if (rj + multerm->size > maxRSize) {
// 				r->size = rj;
// 				maxRSize += (multerm->size << 1);
// 				resizePolynomial_AA(r, maxRSize);
// 			}
// 			//copy multerm coef data directly into r
// 			memcpy(r->elems + rj, multerm->elems, sizeof(AAElem_t)*multerm->size);

// 			//update degs in-place by degree of the mainvar, eps.
// 			for (int idx = rj; idx < rj + multerm->size; ++idx) {
// 				// r->elems[idx].degs += ((degrees_t)eps << mvarOffset);
// 			}

// 			rj += multerm->size;

// 		}		
// 		free(multerm);
// 	}
	

// 	if (ai > 0) {
// 		a->size = ai;
// 		a->alloc = maxASize;
// 	} else {
// 		freePolynomial_AA(a);
// 		a = NULL;
// 	}

// 	if (rj > 0) {
// 		r->size = rj;
// 		r->alloc = maxRSize;
// 	} else {
// 		freePolynomial_AA(r);
// 		r = NULL;
// 	}

// 	if (!lazy) {
// 		int d = c->elems->exp - b->elems->exp + 1 - i;
// 		i += d;
// 		for (int j = 0; j < d; ++j) {
// 			if (ai > 0) {
// 				a = multiplyPolynomials_AA_inp(a, &h, nvar);
// 			}
// 			if (rj > 0) {
// 				r = multiplyPolynomials_AA_inp(r, &h, nvar);
// 			}
// 			hI = multiplyPolynomials_AA_inp(hI, &h, nvar);
// 		}
// #if PDIVIDE_DIVISBLE_CHECK
// 		//r needs one additional one since a and h have one extra from beginning of loop
// 		if (rj > 0) {
// 			r = multiplyPolynomials_AA_inp(r, &h, nvar);	
// 		}
// #endif
// 	}

// 	*res_a = a;
// 	*res_r = r;

// 	//Return number of division steps;
// 	if (e != NULL) {
// 		*e = i;
// 	}

// 	//Return the initial to power of e in hPow;
// 	if (hPow != NULL) {
// 		*hPow = hI;
// 	}
// }


// void pesudoDivide_RecArray_unpk(RecArr_t* c, RecArr_t* b, AltArr_t** res_a, AltArr_t** res_r, int* e, AltArr_t** hPow, int nvar, int lazy) {

// 	if (b == NULL || b->size == 0) {
// 		fprintf(stderr, "Division by zero! Exiting...\n");
// 		exit(EXIT_FAILURE);
// 	}

// 	if (c == NULL || c->size == 0) {
// 		*res_a = NULL;
// 		*res_r = NULL;
// 		return;
// 	}

// 	if (b->elems->exp > c->elems->exp) {
// 		RecArr_t* recR = deepCopyRecArrayPolynomial(c, nvar); 
// 		*res_r = convertFromRecursiveArray(recR, nvar);
// 		*res_a = NULL;
// 		if (e != NULL) {
// 			*e = 0;
// 		}
// 		if (hPow != NULL) {
// 			*hPow = makePolynomial_AA(1, nvar);
// 			mpq_init((*hPow)->elems->coef);
// 			mpq_set_si((*hPow)->elems->coef, 1l, 1l);
// 			(*hPow)->elems->degs = 0;
// 			(*hPow)->size = 1;
// 		}
// 		return;
// 	}

// 	if (!c->unpacked && !b->unpacked) {
// 		pesudoDivide_RecArray(c, b, res_a, res_r, e, hPow, nvar, lazy);
// 		return;
// 	}

// 	if (b->size == 1) {
// 		pesudoDivideOneTerm_RecArray_unpk(c, b, res_a, res_r, e, hPow, nvar, lazy);
// 		return;
// 	}

// 	int unpackedC = 0, unpackedB = 0;
// 	if (!c->unpacked) {
// 		unpackExponentVectors_RecArray_inp(c, nvar);
// 		unpackedC = 1;
// 	}
// 	if (!b->unpacked) {
// 		unpackExponentVectors_RecArray_inp(b, nvar);
// 		unpackedB = 1;
// 	}

// 	RecArrElem_t* __restrict__ k = c->elems;
// 	RecArrElem_t* __restrict__ b2Elem = b->elems + 1;
// 	RecArrElem_t* __restrict__ lenK = k + c->size;

// 	register int maxASize = c->size < 2 ? 2 : c->size;
// 	register int maxRSize = maxASize;
// 	register int ai = 0;
// 	register int rj = 0;

// 	AAElem_t* r = (AAElem_t*) malloc(sizeof(AAElem_t)*maxRSize);
// 	degree_t* rDegs = (degree_t*) malloc(sizeof(degree_t)*maxRSize*nvar);
// 	RecArrElem_t* aElems = (RecArrElem_t*) malloc(sizeof(RecArrElem_t)*maxASize); 
// 	RecArrElem_t* __restrict__ curA = aElems;


// 	//manually do first div as we know it comes from first term of c;
// 	int i = 1;
// 	AltArr_t h;
// 	h.alloc = h.size = b->elems->coefSize;
// 	h.elems = b->elems->coef;
// 	h.nvar = nvar;
// 	h.unpacked = 1;

// 	degree_t bDeg = b->elems->exp, eps = k->exp - bDeg;
// 	AltArr_t kCoef;
// 	kCoef.alloc = kCoef.size = k->coefSize;
// 	kCoef.elems = k->coef;
// 	kCoef.nvar = nvar;
// 	kCoef.unpacked = 1;

// 	AltArr_t* multerm = deepCopyPolynomial_AA(&kCoef);
// 	curA->coef = multerm->elems;
// 	curA->coefSize = multerm->size;
// 	curA->exp = eps;
// 	++curA;
// 	++ai;
// 	free(multerm);
// 	++k;

// 	AltArr_t* multerm2 = NULL;
// 	//use long long so it can fit an entire degree_t but can also be compared to -1.
// 	degree_t delta = 1;
// 	degree_t kDeg = (k == lenK) ? -1 : k->exp;

	
// 	AltArr_t* hI = deepCopyPolynomial_AA(&h);

// 	RecProdHeap_AA* prodHeap = recProdHeapCreate_AA(nvar);
// 	prodHeap->unpacked = 1;
// 	recProdheapResize_AA(prodHeap, maxASize);
// 	prodHeap->lastB = b->size - 1;
// 	recProdHeapInsert_AA(prodHeap, recProdHeapMakeChain_AA(0, 1, NULL), aElems->exp + b2Elem->exp);

// 	while (delta > -1 || kDeg > -1) {
		
// 		//get the leading term of dividend quotient-product difference
// 		delta = recProdHeapPeek_AA(prodHeap);
// 		if (delta == -1 && kDeg == -1) {
// 			break;
// 		}
	
// 		if (delta > kDeg) {
// 			eps = delta;
// 			multerm = recProdHeapGetNextCoef_AA(prodHeap, aElems, b->elems);			
// 			if (multerm->size == 0) {
// 				freePolynomial_AA(multerm);
// 				//the elements in the product with degree delta
// 				// ended up canceling out coefs
// 				continue; 
// 			}
// 			negatePolynomial_AA(multerm);
// 		} else if (delta < kDeg) {
// 			eps = kDeg;
// 			kCoef.alloc = kCoef.size = k->coefSize;
// 			kCoef.elems = k->coef;

// 			//then ctilde is from dividend
// 			multerm = multiplyPolynomials_AA(&kCoef, hI, nvar);
// 			++k;
// 			kDeg = (k == lenK) ? -1 : k->exp;
// 		} else {
// 			eps = delta;
// 			//combine both
// 			multerm2 = recProdHeapGetNextCoef_AA(prodHeap, aElems, b->elems);			
// 			if (multerm2->size == 0) {
// 				freePolynomial_AA(multerm2);
// 				continue;
// 			}			
// 			kCoef.alloc = kCoef.size = k->coefSize;
// 			kCoef.elems = k->coef;

// 		  	multerm = multiplyPolynomials_AA(&kCoef, hI, nvar);
// 			multerm = subPolynomials_AA_inp(multerm, multerm2, nvar);
// 			freePolynomial_AA(multerm2);
			
// 			++k;
// 			kDeg = (k == lenK) ? -1 : k->exp;
			
// 			if (multerm->size == 0) {
// 				//if sub resulted in a zero then we we must get a new multerm
// 				freePolynomial_AA(multerm);
// 				continue;
// 			}
// 		}

// 		//multerm is now the leading coef (eps the degree) of the current difference
// 		if (eps >= bDeg) {

// #if PDIVIDE_DIVISBLE_CHECK
// 			AltArr_t* multermQuo = NULL;
// 			if (divideTest_AA(multerm, &h, &multermQuo, nvar)) {
// 				freePolynomial_AA(multerm);
// 				multerm = multermQuo;			
// 			} else {
// #endif
// 				recArrayMultiplyByCoef_inp_unpk(aElems, ai, &h, nvar);
				
// 				//update h counter
// 				//this is only done when we actually add a new element to the quotient.
// 				++i;
// 				hI = multiplyPolynomials_AA_inp(hI, &h, nvar); //in-place wrt first arg.

// 				// IF the above rec node mult was not done in place, then
// 				// we would need to do the below. Since we delay the actual
// 				// calculation of the coef product from heap elements until just 
// 				// before extraction.
// 				// recProdHeapUpdateQuoByH(prodHeap, h);			
// #if PDIVIDE_DIVISBLE_CHECK
// 			}
// #endif
// 			//resize if necessary for insertion to a
// 			if (ai >= maxASize) {
// 				maxASize <<= 1;
// 				aElems = (RecArrElem_t*) realloc(aElems, sizeof(RecArrElem_t)*maxASize);
// 				curA = aElems + ai;

// 				recProdheapResize_AA(prodHeap, maxASize);
// 			}

// 			fprintf(stderr, "adding to a multerm->unpacked: %d\n", multerm->unpacked);
// 			fprintf(stderr, "degs: %llx\n", multerm->elems->degs);

// 			curA->coef = multerm->elems;
// 			curA->coefSize = multerm->size;
// 			curA->exp = eps - bDeg;
// 			free(multerm);

// 			//insert b_2! Since we constructed a to exactly cancel b_1, useless.
// 			recProdHeapInsert_AA(prodHeap, recProdHeapMakeChain_AA(ai, 1, NULL), curA->exp + b2Elem->exp);
// 			++ai;
// 			++curA;


// 		} else {
// 			if (rj + multerm->size > maxRSize) {
// 				//accumulate remainder.
// 				maxRSize += (multerm->size << 1);
// 				r = realloc(r, sizeof(AAElem_t)*maxRSize);
// 				degree_t* oldDegs = rDegs;
// 				rDegs = (degree_t*) realloc(rDegs, sizeof(degree_t)*maxRSize*nvar);
// 				if (oldDegs != rDegs) {
// 					for (int idx = 0; idx < rj; ++idx) {
// 						r[idx].degs = (degrees_t) (rDegs + idx*nvar);
// 					}
// 				}
// 			}
// 			//copy multerm coef data directly into r
// 			memcpy(r + rj, multerm->elems, sizeof(AAElem_t)*multerm->size);
// 			memcpy(rDegs + rj, (degree_t*) multerm->elems->degs, sizeof(degree_t)*multerm->size*nvar);

// 			//update degs in-place by degree of the mainvar, eps.
// 			for (int idx = rj; idx < rj + multerm->size; ++idx) {
// 				r[idx].degs = (degrees_t) (rDegs + idx*nvar);
// 				rDegs[idx*nvar] = eps;
// 			}

// 			rj += multerm->size;

// 			free(multerm);
// 		}
// 	}


// 	//Condense coefs so they lay in the same array.
// 	int aSize = 0;
// 	for (int idx = 0; idx < ai; ++idx) {
// 		aSize += aElems[idx].coefSize;
// 	}
// 	AAElem_t* aBlock = (AAElem_t*) malloc(sizeof(AAElem_t)*aSize);
// 	degree_t* aDegs = (degree_t*) malloc(sizeof(degree_t)*aSize*nvar);
// 	aSize = 0;
// 	//first loop over each recursive element
// 	for (int idx = 0; idx < ai; ++idx) {
// 		memcpy(aBlock+aSize, aElems[idx].coef, sizeof(AAElem_t)*(aElems[idx].coefSize));
// 		memcpy(aDegs + aSize*nvar, (degree_t*) aElems[idx].coef->degs, sizeof(degree_t)*(aElems[idx].coefSize)*nvar);

// 		//loop over each element of coefficient and in-place update degs by mainvar exp.
// 		for (int idx2 = 0; idx2 < aElems[idx].coefSize; ++idx2) {
// 			aDegs[(aSize+idx2)*nvar] = aElems[idx].exp;
// 			aBlock[aSize + idx2].degs = (degrees_t) (aDegs + (aSize+idx2)*nvar);
// 			// aBlock[aSize + idx2].degs += ((degrees_t) aElems[idx].exp << mvarOffset);
// 		}

// 		free((degree_t*) aElems[idx].coef->degs);
// 		free(aElems[idx].coef);
// 		aSize += aElems[idx].coefSize;
// 	}
// 	free(aElems);


// 	AltArr_t* aa = NULL;
// 	if (aSize > 0) {
// 		aa = (AltArr_t*) malloc(sizeof(AltArr_t));
// 		aa->elems = aBlock;
// 		aa->size = aSize;
// 		aa->alloc = aSize;
// 		aa->nvar = nvar;
// 		aa->unpacked = 1;
// 	}
// 	AltArr_t* ra = NULL;
// 	if (rj > 0) {
// 		ra = (AltArr_t*) malloc(sizeof(AltArr_t));
// 		ra->elems = r;
// 		ra->size = rj;
// 		ra->alloc = maxRSize;
// 		ra->nvar = nvar;
// 		ra->unpacked = 1;
// 	}

// 	if (!lazy) {
// 		int d = c->elems->exp - b->elems->exp + 1 - i;
// 		i += d;
// 		for (int j = 0; j < d; ++j) {
// 			aa = aa == NULL ? NULL : multiplyPolynomials_AA_inp(aa, &h, nvar);
// 			ra = ra == NULL ? NULL : multiplyPolynomials_AA_inp(ra, &h, nvar);
// 			hI = multiplyPolynomials_AA_inp(hI, &h, nvar);
// 		}
// 	}

// 	*res_a = aa;
// 	*res_r = ra;

// 	//Return number of division steps;
// 	if (e != NULL) {
// 		*e = i;
// 	}

// 	//Return the initial to power of e in hPow;
// 	if (hPow != NULL) {
// 		*hPow = hI;
// 	}


// 	if (unpackedC) {
// 		packExponentVectors_RecArray_inp(c, nvar);
// 	}
// 	if (unpackedB) {
// 		packExponentVectors_RecArray_inp(b, nvar);
// 	}
// }















