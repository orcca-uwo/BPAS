

#include "RationalNumberPolynomial/SMQP_Support_Test-AA.h"

void getNextDegrees_unpk (degree_t* prev, degree_t step, degree_t maxUniDeg, int nvar, degree_t* nextDegs) {
	mpz_t weightedDegree;
	mpz_init(weightedDegree);
	mpz_t weight;
	mpz_init(weight);

	for (int i = 0; i < nvar; ++i) {
		mpz_set_si(weight, maxUniDeg);
		mpz_pow_ui(weight, weight, (nvar-1- i));
		// weight = (degree_t) pow(maxUniDeg, (nvar-1 - i));
		mpz_addmul_ui(weightedDegree, weight, prev[i]);
		// unsigned int prev_i = ((prev & masks[i]) >> sizes[i]); 
		// weightedDegree += (weight*prev_i);
	}
	mpz_add_ui(weightedDegree, weightedDegree, step);
	// weightedDegree += step;
	
	mpz_t deg;
	mpz_init(deg);
	for (int i = 0; i < nvar; ++i) {
		mpz_set_si(weight, maxUniDeg);
		mpz_pow_ui(weight, weight, nvar -1 -i);
		// weight = (degree_t) pow(maxUniDeg,(nvar-1 - i));

		mpz_fdiv_q(deg, weightedDegree, weight);
		// unsigned long long int deg = 0;
		// if (weight != 0) {
		// 	deg = weightedDegree / weight;
		// } else {
		// 	deg = weightedDegree;
		// }

		nextDegs[i] = mpz_get_si(deg);
		// nextDegs += (deg << sizes[i]);

		mpz_submul(weightedDegree, weight, deg);
		// weightedDegree -= weight * deg;
	}

	mpz_clear(weight);
	mpz_clear(weightedDegree);
	mpz_clear(deg);

}

degrees_t getNextDegrees (degrees_t prev, degree_t step, degree_t maxUniDeg, int nvar, const int* sizes, const degrees_t* masks) {
	long long int weightedDegree = 0;
	degree_t weight;

	for (int i = 0; i < nvar; ++i) {
		weight = (degree_t) pow(maxUniDeg, (nvar-1 - i));
		unsigned int prev_i = ((prev & masks[i]) >> sizes[i]); 
		weightedDegree += (weight*prev_i);
	}
	weightedDegree += step;

	degrees_t nextDegs = 0; //(degrees_t) malloc(sizeof(degree_t)*nvar);

	for (int i = 0; i < nvar; ++i) {
		weight = (degree_t) pow(maxUniDeg,(nvar-1 - i));
		unsigned long long int deg = 0;
		if (weight != 0) {
			deg = weightedDegree / weight;
		} else {
			deg = weightedDegree;
		}

		nextDegs += (deg << sizes[i]);
		weightedDegree -= weight * deg;
	}

	return nextDegs;
}

AltArr_t* buildRandomPoly_AA_unpk(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);//1529105633;//time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

	degrees_t maxTotalDeg = sparsity * nterms;
	degree_t maxUniDeg = (degree_t) ceil(pow(maxTotalDeg, 1.0 / nvar));
	// maxUniDeg = maxUniDeg < 2 ? 2 : maxUniDeg;

	AltArr_t* aa = makePolynomial_AA_unpk(nterms, nvar);
	degree_t* aDegs = (degree_t*) aa->elems->degs;

	for (int k = 0; k < nvar; ++k) {
		aDegs[k] = 0;
	}
	if (sparsity != 2) {
		aDegs[nvar-1] = rand() % 2; //50/50 chance of first term being a constant;
	}

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}
	while(mpz_sgn(mpzDen) == 0) {
		mpz_urandomb(mpzDen, R_STATE, coefBound);
	}

	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}

	mpq_init(aa->elems[0].coef);
	mpz_set(mpq_numref(aa->elems[0].coef), mpzNum);
	mpz_set(mpq_denref(aa->elems[0].coef), mpzDen);
	mpq_canonicalize(aa->elems[0].coef);
	// mpq_set_si(n->coef, coef_l, 1ul);
	// n->coef = coef_ul;

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
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}

		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpq_init(aa->elems[i].coef);
		mpz_set(mpq_numref(aa->elems[i].coef), mpzNum);
		mpz_set(mpq_denref(aa->elems[i].coef), mpzDen);
		mpq_canonicalize(aa->elems[i].coef);
	}

	mpz_clear(mpzNum);
	mpz_clear(mpzDen);

	//Now reverse the poly so it is in decreasing order;
	aa->size = nterms;
	mergeSortPolynomial_AA(aa);

	return aa;
}

AltArr_t* buildRandomSeededPoly_AA_unpk(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg, time_t seed) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);//1529105633;//time(NULL);
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

	AltArr_t* aa = makePolynomial_AA_unpk(nterms, nvar);
	degree_t* aDegs = (degree_t*) aa->elems->degs;

	for (int k = 0; k < nvar; ++k) {
		aDegs[k] = 0;
	}
	if (sparsity != 2) {
		aDegs[nvar-1] = rand() % 2; //50/50 chance of first term being a constant;
	}

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}
	while(mpz_sgn(mpzDen) == 0) {
		mpz_urandomb(mpzDen, R_STATE, coefBound);
	}

	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}

	mpq_init(aa->elems[0].coef);
	mpz_set(mpq_numref(aa->elems[0].coef), mpzNum);
	mpz_set(mpq_denref(aa->elems[0].coef), mpzDen);
	mpq_canonicalize(aa->elems[0].coef);
	// mpq_set_si(n->coef, coef_l, 1ul);
	// n->coef = coef_ul;

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
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}

		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpq_init(aa->elems[i].coef);
		mpz_set(mpq_numref(aa->elems[i].coef), mpzNum);
		mpz_set(mpq_denref(aa->elems[i].coef), mpzDen);
		mpq_canonicalize(aa->elems[i].coef);
	}

	mpz_clear(mpzNum);
	mpz_clear(mpzDen);

	//Now reverse the poly so it is in decreasing order;
	aa->size = nterms;
	mergeSortPolynomial_AA(aa);

	return aa;
}


Node* buildRandomPoly(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);//1529105633;//time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

    const int* sizes = getExpOffsetArray(nvar);
    const unsigned long long int* masks = getExpMaskArray(nvar);

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

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}
	while(mpz_sgn(mpzDen) == 0) {
		mpz_urandomb(mpzDen, R_STATE, coefBound);
	}

	// long int coef_l = rand() % coefBound;
	// if (coef_l == 0) {
		// coef_l = 1l;
	// }
	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}
	n->degs = degs;
	mpq_init(n->coef);

	mpz_set(mpq_numref(n->coef), mpzNum);
	mpz_set(mpq_denref(n->coef), mpzDen);
	mpq_canonicalize(n->coef);
	// mpq_set_si(n->coef, coef_l, 1ul);
	// n->coef = coef_ul;
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
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}

		// coef_l = 0l; 
		// while (coef_l == 0) {
			// coef_l = rand() % coefBound;
		// }
		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
			// coef_l *= -1;
		}

		n->degs = degs;
		lastDegs = degs;
		mpq_init(n->coef);
		mpz_set(mpq_numref(n->coef), mpzNum);
		mpz_set(mpq_denref(n->coef), mpzDen);
		mpq_canonicalize(n->coef);
		// mpq_set_si(n->coef, coef_l, 1l);
		// n->coef = coef_ul;

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

AltArr_t* buildRandomPolyFromMax_seeded(time_t seed, int nvar, const int* maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) {

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = seed;
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "build random poly seed: %lu\n", t);
		initRand = 1;
	}

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	unsigned long long int maxTerms = 1;
	for (int i = 0; i < nvar; ++i) {
		maxTerms *= (maxDegs[i]+1);
	}
	++maxTerms; //+1 for the contstant term;

	AltArr_t* res = makePolynomial_AA(maxTerms, nvar);
	AAElem_t* elems = res->elems;
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
					mpz_urandomb(mpzDen, R_STATE, coefBound);
					while(mpz_sgn(mpzDen) == 0) {
						mpz_urandomb(mpzDen, R_STATE, coefBound);
					}
					if (includeNeg && rand() % 2) {
						mpz_neg(mpzNum, mpzNum);
					}

					mpq_init(elems[curIdx].coef);
					mpz_set(mpq_numref(elems[curIdx].coef), mpzNum);
					mpz_set(mpq_denref(elems[curIdx].coef), mpzDen);
					mpq_canonicalize(elems[curIdx].coef);

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
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}
		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpq_init(elems[0].coef);
		mpz_set(mpq_numref(elems[0].coef), mpzNum);
		mpz_set(mpq_denref(elems[0].coef), mpzDen);
		mpq_canonicalize(elems[0].coef);

		elems[0].degs = 0;
		for (int k = 0; k < nvar; ++k){
			elems[0].degs |= ( ((degrees_t)maxDegs[k]) << sizes[k]);
		}

		res->size = 1;
	} else {
		//we are in the general case
		unsigned long long int targetTerms = ceil(maxTerms * (1 - sparsity));
	
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
			mpz_urandomb(mpzDen, R_STATE, coefBound);
			while(mpz_sgn(mpzDen) == 0) {
				mpz_urandomb(mpzDen, R_STATE, coefBound);
			}
			if (includeNeg && rand() % 2) {
				mpz_neg(mpzNum, mpzNum);
			}

			mpq_init(elems[i].coef);
			mpz_set(mpq_numref(elems[i].coef), mpzNum);
			mpz_set(mpq_denref(elems[i].coef), mpzDen);
			mpq_canonicalize(elems[i].coef);

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

		sortPolynomial_AA(res);
	}

	mpz_clear(mpzNum);
	mpz_clear(mpzDen);

	resizePolynomial_AA(res, res->size);
	return res;

}

AltArr_t* buildRandomPolyFromMax(int nvar, const int* maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) {

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);//1527609229;//time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "build random poly seed: %lu\n", t);
		initRand = 1;
	}

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	unsigned long long int maxTerms = 1;
	for (int i = 0; i < nvar; ++i) {
		maxTerms *= (maxDegs[i]+1);
	}
	++maxTerms; //+1 for the contstant term;

	AltArr_t* res = makePolynomial_AA(maxTerms, nvar);
	AAElem_t* elems = res->elems;
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
					mpz_urandomb(mpzDen, R_STATE, coefBound);
					while(mpz_sgn(mpzDen) == 0) {
						mpz_urandomb(mpzDen, R_STATE, coefBound);
					}
					if (includeNeg && rand() % 2) {
						mpz_neg(mpzNum, mpzNum);
					}

					mpq_init(elems[curIdx].coef);
					mpz_set(mpq_numref(elems[curIdx].coef), mpzNum);
					mpz_set(mpq_denref(elems[curIdx].coef), mpzDen);
					mpq_canonicalize(elems[curIdx].coef);

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
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}
		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpq_init(elems[0].coef);
		mpz_set(mpq_numref(elems[0].coef), mpzNum);
		mpz_set(mpq_denref(elems[0].coef), mpzDen);
		mpq_canonicalize(elems[0].coef);

		elems[0].degs = 0;
		for (int k = 0; k < nvar; ++k){
			elems[0].degs |= ( ((degrees_t)maxDegs[k]) << sizes[k]);
		}

		res->size = 1;
	} else {
		//we are in the general case
		unsigned long long int targetTerms = ceil(maxTerms * (1 - sparsity));
	
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
			mpz_urandomb(mpzDen, R_STATE, coefBound);
			while(mpz_sgn(mpzDen) == 0) {
				mpz_urandomb(mpzDen, R_STATE, coefBound);
			}
			if (includeNeg && rand() % 2) {
				mpz_neg(mpzNum, mpzNum);
			}

			mpq_init(elems[i].coef);
			mpz_set(mpq_numref(elems[i].coef), mpzNum);
			mpz_set(mpq_denref(elems[i].coef), mpzDen);
			mpq_canonicalize(elems[i].coef);

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

		sortPolynomial_AA(res);
	}

	mpz_clear(mpzNum);
	mpz_clear(mpzDen);

	resizePolynomial_AA(res, res->size);
	return res;

}

















/****************
* Random Polynomial Algorithm with Maximum Degree
*****************/

/**
 * @param[in] lastDegs Previous degree
 * @param[in] step The distance between the previous and returned degree
 * @param[in] nvar The number of variables
 * \example if the lastDegs = (2 2 4) and step = 3 then output is (1 2 4), and
 *  the sequence of in/out-puts is -> (0 2 4) -> (0 0 4) -> (0 0 3) -> (0 0 2).
 */
degrees_t getNextMaxDegs(degrees_t lastDegs, degrees_t step, int nvar, const int* sizes, const degrees_t* masks)
{
        degrees_t nextDegs = 0; // (degrees_t) calloc (nvar, sizeof(degree_t));
	degrees_t weight = step;


	// fprintf(stderr, "weight := %llx \n", weight); // test!
	// fprintf(stderr, "lastDegs := %llx \n", lastDegs); // test!
	
	for (int i = 0; i < nvar ; ++i)
	  {
	    unsigned long long int lastDegs_i = ((lastDegs & masks[i]) >> sizes[i]);
			// fprintf(stderr, "lastDegs_i := %llx \n", lastDegs_i); // test!
	    unsigned long long int deg = 0;
	    if (weight > 0)
		{
			if (lastDegs_i >= weight){
				deg = lastDegs_i - weight;
				weight = 0;
			} else {
				weight -= lastDegs_i;
				deg = 0;
			}

		} else {
		        deg = lastDegs_i;
		}
	    nextDegs += (deg << sizes[i]);
	}
	return nextDegs;
}

/**
 * Random Polynomial Algorithm with MaxDegs
 * @param[in] nvar The number of variables
 * @param[in] maxDegs The Maximum degrees of the output polynomial
 * @param[in] sparsity
 * @param[in] coefBound
 * @param[in] includeNeg
 * \brief This algorithm generates a random polynomial where the maximum degree is maxDegs.
 */
Node* randomMaxDegsPoly(int nvar, degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg)
{

	if (nvar < 1)
		fprintf(stderr,"nvar is %d < 1", nvar);

	static int initRand = 0; // initialize seed!
	if (!initRand){
		time_t t = time(NULL);
		srand(t);
		initRand = 1;
	}

	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* masks = getExpMaskArray(nvar);

	// define head and tail:
	Node* head;
	Node* tail;

	Node* n = (Node*) malloc (sizeof(Node));
	degrees_t degs = 0; //(degrees_t) calloc (nvar, sizeof(degree_t));
	degrees_t lastDegs = maxDegs;
	long int coef_l = rand() % coefBound;
	if (coef_l == 0)
		coef_l = 1l;
	if (includeNeg && rand()%2)
		coef_l *= -1;


	// set the first term's degree maxDegs:
	mpq_init (n->coef);
	mpq_set_si(n->coef, coef_l, 1ul);
	n->degs = lastDegs;

	head = n;
	tail = n;

	degrees_t step = 0;
        unsigned int lastDegs_m1 = ((lastDegs & masks[nvar-1]) >> sizes[nvar-1]);
	while (lastDegs_m1 != 0){
		n = (Node*) malloc (sizeof(Node));
		n->next = NULL;

		step = 0;
		do{
		  step = (degrees_t) rand() % sparsity;
		} while(step == 0);
		degs = getNextMaxDegs(lastDegs, step, nvar, sizes, masks);
		// fprintf(stderr, "degs := %llx \n", degs);
		coef_l = 0l;
		do{
			coef_l = rand() % coefBound;
		} while (coef_l == 0);
		if (includeNeg && rand() % 2)
			coef_l *= -1;

		n->degs = degs;
		lastDegs = degs;
		lastDegs_m1 = ((lastDegs & masks[nvar-1]) >> sizes[nvar-1]);
		mpq_init(n->coef);
		mpq_set_si(n->coef, coef_l, 1l);

		tail->next = n;
		tail = n;
	}

	return head;
}

/***************
* Random Triangular Set
***************/

/**
 * Convert the exponents of a polynomial
 * @param[in] inPoly The polynomial
 * @param[in] inNvar The number of variables of inPoly
 * @param[in] outNvar The number of variables of output polynomial.
 */
/*
void convertExponentPoly(Node* inPoly, int inNvar, int outNvar)
{
	if (inNvar == outNvar)
		return;

	degrees_t outDegs = (degrees_t) calloc (outNvar, sizeof(degree_t));
	Node* cur = inPoly;

	while(cur != NULL)
	{
		// Ex: (w,v) ---> (x=0,y=0,z=0,w,v)
		for(int i = 0; i < inNvar; ++i)
			outDegs[outNvar - inNvar + i] = cur->degs[i];

		free(cur->degs);
		cur->degs = outDegs;
		outDegs = (degrees_t) calloc (outNvar, sizeof(degree_t));
		cur = cur->next;
	}
	return;
}
*/

/** Build a Random Triangular Set based on maxDegs
 * @param[in] T A pointer to the triangular set
 * @param[in] outNvar The number of variables of output polynomial
 * @param[in] maxDegs The Maximum degree of polynomials
 * @param[in] sparsity
 * @param[in] coefBound
 * @param[in] includeNeg
 * \brief Output is a Lazard or general triangular set (w.r.t the lazard_flag)
 * such that if lazard_flag != 0, then the output set is Lazard triangular set.
 */
void randomTriangularSet (Node** T, int outNvar,degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg, int lazard_flag)
{
	time_t sed = time(NULL);
	srand(sed);

    const int* sizes = getExpOffsetArray(outNvar);
	const degrees_t* masks = getExpMaskArray(outNvar);

	// generate the members
	for (int i = 0; i < outNvar; i++)
	{
		if (i != outNvar-1){
		  	degrees_t tmpDegs = 0; //(degrees_t) calloc (i+1, sizeof(degree_t));
			    
		    for (int j = outNvar-i-1; j < outNvar; ++j){
		      unsigned long long int maxDegs_o = ((maxDegs & masks[j]) >> sizes[j]);
		      tmpDegs += (maxDegs_o << sizes[j]);
		    }
			T[i] = randomMaxDegsPoly(outNvar, tmpDegs, sparsity, coefBound, includeNeg);
		} else {
			T[i] = randomMaxDegsPoly(outNvar, maxDegs, sparsity, coefBound, includeNeg);
		}

		//Lazard Triangular Set:
		if (lazard_flag != 0){ 
		 	mpq_set_si(T[i]->coef, 1l, 1ul); 
		 	for (int j = outNvar-i; j < outNvar; ++j){
		 		unsigned long long int degs = ((T[i]->degs & masks[j]) >> sizes[j]); 
		 		T[i]->degs -= (degs << sizes[j]);
		 	}
		 } 
	}
	return;
}

void randomTriangularSet_AA (AltArr_t** T, int outNvar,degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg, int lazard_flag)
{
	time_t sed = time(NULL);
	srand(sed);

    const int* sizes = getExpOffsetArray(outNvar);
	const degrees_t* masks = getExpMaskArray(outNvar);

	// generate the members
	for (int i = 0; i < outNvar; i++)
	{
		if (i != outNvar-1){
		  	degrees_t tmpDegs = 0; //(degrees_t) calloc (i+1, sizeof(degree_t));
			    
		    for (int j = outNvar-i-1; j < outNvar; ++j){
		      unsigned long long int maxDegs_o = ((maxDegs & masks[j]) >> sizes[j]);
		      tmpDegs += (maxDegs_o << sizes[j]);
		    }
		    //fprintf(stderr, "tmpDegs := %llx \n", tmpDegs);
			T[i] = deepCopyPolynomial_AAFromNode(randomMaxDegsPoly(outNvar, tmpDegs, sparsity, coefBound, includeNeg), outNvar);
		} else {
			T[i] = deepCopyPolynomial_AAFromNode(randomMaxDegsPoly(outNvar, maxDegs, sparsity, coefBound, includeNeg),outNvar);
		}

		//Lazard Triangular Set:
		if (lazard_flag != 0){ 
		 	mpq_set_si(T[i]->elems[0].coef, 1l, 1ul);
		 	for (int j = outNvar-i; j < outNvar; ++j){
		 		unsigned long long int degs = ((T[i]->elems[0].degs & masks[j]) >> sizes[j]);
		 		T[i]->elems[0].degs -= (degs << sizes[j]);
		 	}
		 } 
	}
	return;
}
