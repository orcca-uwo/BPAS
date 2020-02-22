

#include "IntegerPolynomial/SMZP_Support_Test.h"
#include "RationalNumberPolynomial/SMQP_Support_Test-AA.h"

Node* buildRandomZPoly(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

    const int* sizes = getExpOffsetArray(nvar);
    const degrees_t* masks = getExpMaskArray(nvar);

	degrees_t maxTotalDeg = sparsity * nterms;
	degree_t maxUniDeg = (degree_t) ceil(pow(maxTotalDeg, 1.0 / nvar));
	// maxUniDeg = maxUniDeg < 2 ? 2 : maxUniDeg;

	Node* head;
	Node* tail;

	Node* n = (Node*) malloc(sizeof(Node));
	n->next = NULL;
	degrees_t degs = 0;// (degrees_t) calloc(nvar, sizeof(degree_t));
	if (sparsity == 2) {
		//if dense always include the constant
		degs = 0;
	} else {
		degs = rand() % 2; 		//50/50 chance of first term being a constant;
	}
	degrees_t lastDegs = degs;

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);
	mpz_set_si(mpzDen, 1l);

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}
	
	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
	}
	n->degs = degs;
	mpq_init(n->coef);

	mpz_set(mpq_numref(n->coef), mpzNum);
	mpz_set(mpq_denref(n->coef), mpzDen);
	mpq_canonicalize(n->coef);
	
	head = n;
	tail = n;
	--nterms;

	degree_t step = 0;
	for(int i = 0; i < nterms; ++i) {
		n = (Node*) malloc(sizeof(Node));
		n->next = NULL;

		step = 0;
		while(step == 0) {
			step = rand() % sparsity;
			step = step < sparsity / 2 ? step + (sparsity / 2) : step;
			// step %= sparsity;
		}
		degs = getNextDegrees(lastDegs, step, maxUniDeg,nvar, sizes, masks);

		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while(mpz_sgn(mpzNum) == 0) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}

		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
			// coef_l *= -1;
		}

		n->degs = degs;
		lastDegs = degs;
		mpq_init(n->coef);
		mpz_set(mpq_numref(n->coef), mpzNum);
		mpz_set(mpq_denref(n->coef), mpzDen);

		// den is always 1 here so no need to canonicalize.
		// mpq_canonicalize(n->coef);

		tail->next = n;
		tail = n;
	}

	mpz_clear(mpzNum);
	mpz_clear(mpzDen);

	//Now reverse the poly so it is in decreasing order;
	Node* cur = head;
	Node* prev = head;
	Node* next = head->next;
	cur->next = NULL;
	while(next != NULL) {
		cur = next;
		next = cur->next;
		cur->next = prev;
		prev = cur;
	}
	head = cur;

	return head;
}

AltArrZ_t* buildRandomPoly_AAZ_unpk(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

	degrees_t maxTotalDeg = sparsity * nterms;
	degree_t maxUniDeg = (degree_t) ceil(pow(maxTotalDeg, 1.0 / nvar));
	// maxUniDeg = maxUniDeg < 2 ? 2 : maxUniDeg;

	AltArrZ_t* aa = makePolynomial_AAZ_unpk(nterms, nvar);
	degree_t* aDegs = (degree_t*) aa->elems->degs;

	for (int k = 0; k < nvar; ++k) {
		aDegs[k] = 0;
	}
	if (sparsity != 2) {
		aDegs[nvar-1] = rand() % 2; //50/50 chance of first term being a constant;
	}

	mpz_t mpzNum;
	mpz_init(mpzNum);

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}

	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}

	mpz_init(aa->elems[0].coef);
	mpz_set(aa->elems[0].coef, mpzNum);

	degree_t step = 0;
	for(int i = 1; i < nterms; ++i) {
		step = 0;
		while(step == 0) {
			step = rand() % sparsity;
			step = step < sparsity / 2 ? step + (sparsity / 2) : step;
			// step %= sparsity;
		}
		getNextDegrees_unpk((aDegs + (i-1)*nvar), step, maxUniDeg, nvar, (aDegs + i*nvar));
		aa->elems[i].degs = (degrees_t) (aDegs + i*nvar);

		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while(mpz_sgn(mpzNum) == 0) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}

		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpz_init(aa->elems[i].coef);
		mpz_set(aa->elems[i].coef, mpzNum);
	}

	mpz_clear(mpzNum);

	//Now reverse the poly so it is in decreasing order;
	aa->size = nterms;
	mergeSortPolynomial_AAZ(aa);

	return aa;
}

AltArrZ_t* buildRandomSeededPoly_AAZ_unpk(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg, time_t seed) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

	if (seed > 0) {
		srand(seed);
		gmp_randseed_ui(R_STATE, seed);
	}

	degrees_t maxTotalDeg = sparsity * nterms;
	degree_t maxUniDeg = (degree_t) ceil(pow(maxTotalDeg, 1.0 / nvar));
	// maxUniDeg = maxUniDeg < 2 ? 2 : maxUniDeg;

	AltArrZ_t* aa = makePolynomial_AAZ_unpk(nterms, nvar);
	degree_t* aDegs = (degree_t*) aa->elems->degs;

	for (int k = 0; k < nvar; ++k) {
		aDegs[k] = 0;
	}
	if (sparsity != 2) {
		aDegs[nvar-1] = rand() % 2; //50/50 chance of first term being a constant;
	}

	mpz_t mpzNum;
	mpz_init(mpzNum);

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}

	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}

	mpz_init(aa->elems[0].coef);
	mpz_set(aa->elems[0].coef, mpzNum);

	degree_t step = 0;
	for(int i = 1; i < nterms; ++i) {
		step = 0;
		while(step == 0) {
			step = rand() % sparsity;
			step = step < sparsity / 2 ? step + (sparsity / 2) : step;
			// step %= sparsity;
		}
		getNextDegrees_unpk((aDegs + (i-1)*nvar), step, maxUniDeg, nvar, (aDegs + i*nvar));
		aa->elems[i].degs = (degrees_t) (aDegs + i*nvar);

		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while(mpz_sgn(mpzNum) == 0) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}

		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpz_init(aa->elems[i].coef);
		mpz_set(aa->elems[i].coef, mpzNum);
	}

	mpz_clear(mpzNum);

	//Now reverse the poly so it is in decreasing order;
	aa->size = nterms;
	mergeSortPolynomial_AAZ(aa);

	return aa;
}


AltArrZ_t* buildRandomZPolyFromMax(int nvar, const int* maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) {

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);
/*		t = 1572887766; // freeing an invalid pointer in bivariate factor_test()*/
/*		t = 1572893480; // segfault in factor_prim_sqf_inner in bvariate factor_test() (related to previous?)*/
		// t = 1580604430;
		t = 1580938263;
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

	mpz_t mpzNum;
	mpz_init(mpzNum);

	unsigned long long int maxTerms = 1;
	for (int i = 0; i < nvar; ++i) {
		maxTerms *= (maxDegs[i]+1);
	}
	++maxTerms; //+1 for the contstant term;

	AltArrZ_t* res = makePolynomial_AAZ(maxTerms, nvar);
	AAZElem_t* elems = res->elems;
	const int* sizes = getExpOffsetArray(nvar);
	
	float eps = 1e-6;
	if (sparsity < eps) { 
		//generate dense poly;
		degrees_t degs = 0;
		int curIdx = 0;
		degrees_t curDegs[nvar];
		for (int i = 0; i < nvar; ++i) {
			curDegs[i] = maxDegs[i];
		}

		for (int i = nvar-1; i >= 0; --i) {
			if (i == nvar - 1) {
				for (int j = curDegs[i]; j >= 0; --j) {
					mpz_urandomb(mpzNum, R_STATE, coefBound);
					while(mpz_sgn(mpzNum) == 0) {
						mpz_urandomb(mpzNum, R_STATE, coefBound);
					}
					if (includeNeg && rand() % 2) {
						mpz_neg(mpzNum, mpzNum);
					}

					mpz_init(elems[curIdx].coef);
					mpz_set(elems[curIdx].coef, mpzNum);

					degs = 0;
					for (int k = 0; k < nvar; ++k){
						degs |= (curDegs[k] << sizes[k]);
					}
					elems[curIdx].degs = degs;
					
					--(curDegs[i]);
					++curIdx;
				}
			} else if (curDegs[i] > 0) {
				--(curDegs[i]);
				//after decrementing curDegs[i], set all lower-ordered variables
				//to 0 and start the generation again at i = nvar-1;
				for (int k = i+1; k < nvar; ++k) {
					curDegs[k] = maxDegs[k];
				}

				i = nvar; 
			}
		}
		
		res->size = curIdx;

	} else if (sparsity >= 1.0f) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while(mpz_sgn(mpzNum) == 0) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}
		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpz_init(elems[0].coef);
		mpz_set(elems[0].coef, mpzNum);

		elems[0].degs = 0;
		for (int k = 0; k < nvar; ++k){
			elems[0].degs |= ( ((degrees_t)maxDegs[k]) << sizes[k]);
		}

		res->size = 1;
	} else {
		//we are in the general case
		unsigned long long int targetTerms = ceil(maxTerms * (1.0 - sparsity));
	
		//incremement maxDegs for mod operation below
		int maxDegsL[nvar];
		for (int k = 0; k < nvar; ++k) {
			maxDegsL[k] = maxDegs[k] + 1;
		}

		degrees_t curDeg = 0;

		for (int i = 0; i < targetTerms; ++i) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
			while(mpz_sgn(mpzNum) == 0) {
				mpz_urandomb(mpzNum, R_STATE, coefBound);
			}
			if (includeNeg && rand() % 2) {
				mpz_neg(mpzNum, mpzNum);
			}

			mpz_init(elems[i].coef);
			mpz_set(elems[i].coef, mpzNum);

			elems[i].degs = 0;

			if (i == 0) {
				for (int k = 0; k < nvar; ++k){
					elems[0].degs |= ( ((degrees_t)maxDegs[k]) << sizes[k]);
				}			
			} else {
				for (int k = 0; k < nvar; ++k){
					curDeg = rand() % (maxDegsL[k]);
					elems[i].degs |= ( curDeg << sizes[k]);
				}
			}
		}

		res->size = targetTerms;

		sortPolynomial_AAZ(res);
	}

	mpz_clear(mpzNum);

	resizePolynomial_AAZ(res, res->size);
	return res;
}
