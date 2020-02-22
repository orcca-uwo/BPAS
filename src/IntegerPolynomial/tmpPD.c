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
				char* syms[] = {"w"}
				fprintf(stderr, "h: \n" );
				printPolynomial_AAZ(stderr, h, syms, 1);
				fprintf(stderr, "multerm: \n" );
				printPolynomial_AAZ(stderr, multerm, syms, 1);
				fprintf(stderr, "multermquo: \n" );
				printPolynomial_AAZ(stderr, multermquo, syms, 1);
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