
#include "RationalNumberPolynomial/SMQP_Support-AA.h"
#include "RationalNumberPolynomial/SMQP_SparseInterpolation-AA.h"
#include "LinearAlgebra/SystemSolving-IML.h"
#include <math.h>


/********
 * Prime number helpers
 ********/

//The first 32 primes; the current maximum nvar in SMQP.
static int primeNumbers[] = {
	2,
	3,
	5,
	7,
	11,
	13,
	17,
	19,
	23,
	29,
	31,
	37,
	41,
	43,
	47,
	53,
	59,
	61,
	67,
	71,
	73,
	79,
	83,
	89,
	97,
	101,
	103,
	107,
	109,
	113,
	127,
	131
};

mpq_t* getPrimeMatrix(int T, int n) {
	long size = T*n;
	mpq_t* A = (mpq_t*) malloc(sizeof(mpq_t)*size);
	
	int* primes = &(primeNumbers[0]);

	for (int j = 0; j < n; ++j) {
		mpq_init(A[j]);
		mpq_set_si(A[j], 1l, 1l);
	}
	for (int j = 0; j < n; ++j) {
		mpq_init(A[n + j]);
		mpq_set_si(A[n + j], primes[j], 1l);	
	}

	for (int i = 2; i < T; ++i) {
		for (int j = 0; j < n; ++j) {
			mpq_init(A[i*n + j]);
			mpq_mul(A[i*n + j], A[(i-1)*n + j], A[n + j]);	
		}
	}

	return A;
}



/********
 * MPQ Random Generator
 ********/

MPQRandomGenerator_t* initMPQRandomGenerator (unsigned long seed, size_t numBits, size_t denBits) {
	MPQRandomGenerator_t* mpqRGen = (MPQRandomGenerator_t*) malloc(sizeof(MPQRandomGenerator_t));

	gmp_randinit_default (mpqRGen->R_STATE);
	gmp_randseed_ui(mpqRGen->R_STATE, seed);

	mpqRGen->numBits = numBits;
	mpqRGen->denBits = denBits;

	return mpqRGen;
}

void freeMPQRandomGenerator(MPQRandomGenerator_t* mpqRGen) {
	gmp_randclear(mpqRGen->R_STATE);
	free(mpqRGen);
}

void MPQRandomGeneratorGetQ(MPQRandomGenerator_t* mpqRGen, mpq_t result) {
	mpz_urandomb(mpq_numref(result), mpqRGen->R_STATE, mpqRGen->numBits);
	while (mpz_sgn(mpq_numref(result)) == 0) {
		mpz_urandomb(mpq_numref(result), mpqRGen->R_STATE, mpqRGen->numBits);
	}

	mpz_set_ui(mpq_denref(result), 1l);
}

void MPQRandomGeneratorGetZ(MPQRandomGenerator_t* mpqRGen, mpq_t result) {
	mpz_urandomb(mpq_numref(result), mpqRGen->R_STATE, mpqRGen->numBits);
	while (mpz_sgn(mpq_numref(result)) == 0) {
		mpz_urandomb(mpq_numref(result), mpqRGen->R_STATE, mpqRGen->numBits);
	}

	mpz_set_ui(mpq_denref(result), 1l);
}



/********
 * SMQP BlackBox
 ********/

SMQPBlackBox_t* initSMQPBlackBox(AltArr_t* aa) {
	SMQPBlackBox_t* bb = (SMQPBlackBox_t*) malloc(sizeof(SMQPBlackBox_t));
	bb->aa = deepCopyPolynomial_AA(aa);
	return bb;
}

void freeSMQPBlackBox(SMQPBlackBox_t* bb) {
	freePolynomial_AA(bb->aa);
}

void evalSMQPBlackBox(SMQPBlackBox_t* bb, ratNum_t* points, int nvar, ratNum_t res) {
	evalPolyToVal_AA(bb->aa, points, nvar, res);
}



/********
 * Sparse Interpolation Helpers
 ********/

/** 
 * k repesents the index of offset that should occur for degrees within coefPolys
 * is k-1.
 */
AltArr_t* mergePolysWithSkel(AltArr_t* skel, AltArr_t** coefPolys, int nvar, int k) {
	int totalSize = 0;
	for (int i = 0; i < skel->size; ++i) {
		// fprintf(stderr, "poly size: %d\n", (coefPolys[i])->size );
		totalSize += (coefPolys[i])->size;
	}


	int* offsets = getExpOffsetArray(nvar);

	AltArr_t* aa = makePolynomial_AA(totalSize, nvar);
	int curIdx = 0;
	AltArr_t* coefPoly = NULL;
	for (int i = 0; i < skel->size; ++i) {
		coefPoly = coefPolys[i];
		for (int j = 0; j < coefPoly->size; ++j) {
			aa->elems[curIdx + j].degs = skel->elems[i].degs | (coefPoly->elems[j].degs << offsets[k]);
	
			mpq_init(aa->elems[curIdx+j].coef);
			mpq_set(aa->elems[curIdx+j].coef, coefPoly->elems[j].coef);
		}

		curIdx += coefPoly->size;
	}

	free(offsets);

	aa->size = totalSize;
	return aa;
}

/**
 * Given a monomial as an AAElem_t, evaluate that monomial at the point supplied.
 * k is supplied such that only the first k variables of the monomial need be evaluated.
 *
 * The result is returned in res. 
 */
void evalMonomialToVal_AA(AAElem_t* elem, mpq_t* point, int nvar, int k, mpq_t res) {
	degrees_t* masks = getExpMaskArray(nvar);
	int* sizes = getExpOffsetArray(nvar);

	mpq_set_si(res, 1l, 1l);
	mpq_t varExp;
	mpq_init(varExp);
	degree_t deg;
	for (int i = 0; i < k; ++i) {
		deg = GET_NTH_EXP(elem->degs, masks[i], sizes[i]);
		mpz_pow_ui(mpq_numref(varExp), mpq_numref(point[i]), deg);
		mpz_pow_ui(mpq_denref(varExp), mpq_denref(point[i]), deg);
		mpq_canonicalize(varExp);
		mpq_mul(res, res, varExp);
	}

	mpq_clear(varExp);
	free(sizes);
	free(masks);
}

/**
 * Interpolate the first variable for the probabilistic sparse interpolation.
 */
AltArr_t* stage1Interpolate(int degreeBound, int nvar, SMQPBlackBox_t* bb, mpq_t* initPoint, mpq_t initVal, MPQRandomGenerator_t* mpqRGen) {

	mpq_t temp;
	mpq_init(temp);
	mpq_set(temp, initPoint[0]);

	int nPoints = degreeBound + 1;
	mpq_t* vals = (mpq_t*) malloc(sizeof(mpq_t)*nPoints);
	mpq_t* points = (mpq_t*) malloc(sizeof(mpq_t)*nPoints);

	mpq_init(points[0]);
	mpq_set(points[0], initPoint[0]);
	mpq_init(vals[0]);
	mpq_set(vals[0], initVal);

	for (int i = 1; i < nPoints; ++i) {
		mpq_init(vals[i]);
		mpq_init(points[i]);
	
		MPQRandomGeneratorGetQ(mpqRGen, initPoint[0]);
		mpq_set(points[i], initPoint[0]);
		evalSMQPBlackBox(bb, initPoint, nvar, vals[i]);
	}

	// AltArr_t* aa = univariateInterpolation_AA(nPoints-1, points, vals);
	AltArr_t* aa = univarInterpolate_AA(points, vals, nPoints);
	
	int offset = getMVarExpOffset(nvar);
	for (int i = 0; i < aa->size; ++i) {
		aa->elems[i].degs = (aa->elems[i].degs << offset);
	}

	mpq_set(initPoint[0], temp);
	mpq_clear(temp);

	return aa;
}

/** 
 * Interpolate the first variable for the deterministic sparse interpolation.
 */
AltArr_t* stage1InterpolatePrimes(int degreeBound, int nvar, SMQPBlackBox_t* bb, mpq_t* primeMat, int T) {

	mpq_t temp;
	mpq_init(temp);

	int nPoints = degreeBound + 1;
	mpq_t* vals = (mpq_t*) malloc(sizeof(mpq_t)*nPoints);
	mpq_t* points = (mpq_t*) malloc(sizeof(mpq_t)*nPoints);

	//The points will always be the first nPoints of first column of prime mat
	for (int i = 0; i < nPoints; ++i) {
		mpq_init(vals[i]);
		mpq_init(points[i]);
		mpq_set(points[i], primeMat[i*nvar]);
	}

	AltArr_t* longestPoly = NULL;
	//iterate 0..T-1 as the tail ends of the evalPoint (x, p_2^exp, p_3^exp, ...)
	//this corresponds to using row exp of primeMat.
	for (int exp = 0; exp < T; ++exp) {
		mpq_set(temp, primeMat[exp*nvar]);
		for (int i = 0; i < nPoints; ++i) {	
			mpq_set(primeMat[exp*nvar], points[i]); //points[i] is p_1^i.

			//use first row of primeMat as evalPoint (2^i, p_2^exp, ...);
			evalSMQPBlackBox(bb, primeMat + exp*nvar, nvar, vals[i]);
		}
		mpq_set(primeMat[exp*nvar], temp);

		// AltArr_t* aa = univariateInterpolation_AA(nPoints-1, points, vals);
		AltArr_t* aa = univarInterpolate_AA(points, vals, nPoints);
		if (longestPoly == NULL || aa->size > longestPoly->size) {
			freePolynomial_AA(longestPoly);
			longestPoly = aa;
		} else {
			freePolynomial_AA(aa);
		}
	}
	
	int offset = getMVarExpOffset(nvar);
	for (int i = 0; i < longestPoly->size; ++i) {
		longestPoly->elems[i].degs = (longestPoly->elems[i].degs << offset);
	}

	mpq_set(primeMat[0], temp);
	mpq_clear(temp);

	for (int i = 0; i < nPoints; ++i) {
		mpq_clear(points[i]);
		mpq_clear(vals[i]);
	}
	free(vals);
	free(points);

	return longestPoly;
}

AltArr_t* getPolynomialPJ(SMQPBlackBox_t* bb, AltArr_t* skel, int k, int degreeBound, int nvar, mpq_t* initPoint, mpq_t rj, MPQRandomGenerator_t* mpqRGen) {
	
	int t = skel->size;
	mpq_t* A = getMPQIdentity(t);
	mpq_t b[t];
	for (int i = 0; i < t; ++i) {
		mpq_init(b[i]);
	}

	mpq_t point[nvar];
	for (int i = 0; i < nvar; ++i) {
		mpq_init(point[i]);
	}

	mpq_set(point[k], rj);
	for (int i = k+1; i < nvar; ++i) {
		mpq_set(point[i], initPoint[i]);
	}

	//Set up linear equations using the same initialPoint, but varying
	//values for the variable we are trying to interpolate. 
	for (int i = 0; i < t; ++i) {
		for (int l = 0; l < k; ++l) {
			MPQRandomGeneratorGetQ(mpqRGen, point[l]);
		}
		for (int j = 0; j < t; ++j) {
			evalMonomialToVal_AA(&(skel->elems[j]), point, nvar, k, A[i*t + j]);
		}
		evalSMQPBlackBox(bb, point, nvar, b[i]);
	}

	mpq_t X[t];
	for (int i = 0; i < t; ++i) {
		mpq_init(X[i]);
	}
	long* rowList = NULL;

	long rank = solveMPQSystem(A, b, t, X, &rowList);

	if (rank < t) {
		// fprintf(stderr, "Failed to get PJ poly because rank-deficient! Rank: %ld, T: %d\n", rank, t);
		int itr = 0;
		while(1) {
			for (int m = rank; m < t; ++m) {
				int curRow = rowList[m];
				for (int l = 0; l < k; ++l) {
					MPQRandomGeneratorGetQ(mpqRGen, point[l]);
				}
				for (int j = 0; j < t; ++j) {
					evalMonomialToVal_AA(&(skel->elems[j]), point, nvar, k, A[curRow*t + j]);
				}
				evalSMQPBlackBox(bb, point, nvar, b[curRow]);
			}

			free(rowList);
			rowList = NULL;
			rank = solveMPQSystem(A, b, t, X, &rowList);
			if (rank == t) {
				break;
			}

			++itr;
			if (itr > 1000) {
				fprintf(stderr, "After 1000 iterations, could still not find non-singular system!\n");
				exit(1);
			}
		}
	}

	free(rowList);

	AltArr_t* aa = deepCopyPolynomial_AA(skel);
	for (int i = 0; i < t; ++i) {
		mpq_set(aa->elems[i].coef, X[i]);
		mpq_clear(X[i]);
	}
	for (int i = 0; i < t; ++i) {
		mpq_clear(b[i]);
	}
	for (int i = 0; i < nvar; ++i) {
		mpq_clear(point[i]);
	}
	freeMPQMat(A, t, t);

	return aa;
}

AltArr_t* getPolynomialPJPrimes(SMQPBlackBox_t* bb, AltArr_t* skel, int k, int degreeBound, int nvar, mpq_t* primeMat, mpq_t* initPoint, mpq_t rj) {
	int t = skel->size;

	mpq_t* A = getMPQIdentity(t);
	mpq_t b[t];
	for (int i = 0; i < t; ++i) {
		mpq_init(b[i]);
	}

	mpq_t point[nvar];
	for (int i = 0; i < nvar; ++i) {
		mpq_init(point[i]);
	}
	mpq_set(point[k], rj);
	for (int i = k+1; i < nvar; ++i) {
		mpq_set(point[i], initPoint[i]);
	}

	for (int i = 0; i < t; ++i) {
		//instead of random use first k-1 primes, each to the power of i.
		for (int l = 0; l < k; ++l) {
			mpq_set(point[l], primeMat[i*nvar + l]);
		}
		for (int j = 0; j < t; ++j) {
			evalMonomialToVal_AA(&(skel->elems[j]), point, nvar, k, A[i*t + j]);
		}
		evalSMQPBlackBox(bb, point, nvar, b[i]);
	}

	mpq_t X[t];
	for (int i = 0; i < t; ++i) {
		mpq_init(X[i]);
	}
	long* rowList = NULL;

	long rank = solveMPQSystem(A, b, t, X, &rowList);

	if (rank < t) {
		//this is theoretically impossible due to our construction of the matrix. 
		//But, this can occur due to a bug in IML. The matrix is very nearly singular.
		int itr = 0;
		while(1) {
			free(rowList);
			rowList = NULL;
			rank = solveMPQSystem(A, b, t, X, &rowList);
			if (rank == t) {
				break;
			}

			++itr;
			if (itr > 1000) {
				fprintf(stderr, "Failed to interpolate a polynomial due to a bug in IML!\n");
				exit(1);
			}
		}
	}

	free(rowList);

	AltArr_t* aa = deepCopyPolynomial_AA(skel);
	for (int i = 0; i < t; ++i) {
		mpq_set(aa->elems[i].coef, X[i]);
		mpq_clear(X[i]);
	}
	for (int i = 0; i < t; ++i) {
		mpq_clear(b[i]);
	}
	for (int i = 0; i < nvar; ++i) {
		mpq_clear(point[i]);
	}
	freeMPQMat(A, t, t);

	return aa;
}

AltArr_t* multivariateInterpolate_AA(SMQPBlackBox_t* bb, int* degreeBounds, int nvar, int T, double eps) {
	int maxD = 0;
	for (int i = 0; i < nvar; ++i) {
		maxD = maxD < degreeBounds[i] ? degreeBounds[i] : maxD;
	}

	// TODO, use T= prod_i^nvar (degreeBounds[i] +1 )
	if (T < 0) {
		T = (int) (pow(maxD, nvar) + 0.5);
	}

	mpq_t mpEps;
	mpq_init(mpEps);
	mpq_set_d(mpEps, eps);
	mpz_t mpD;
	mpz_init(mpD);
	mpz_set_si(mpD, maxD);

	mpz_t size;
	mpz_init(size);
	mpz_set_si(size, T); //T
	mpz_mul(size, size, mpD); // d*T
	mpz_set_si(mpD, nvar);  
	mpz_mul(size, size, mpD); // ndT
	mpz_mul(size, size, mpD); // n^2dT

	mpz_mul(size, size, mpq_denref(mpEps));
	mpz_cdiv_q(size, size, mpq_numref(mpEps)); // n^2dT / eps;
	size_t nBits = mpz_sizeinbase(size, 2);

	// gmp_fprintf(stderr, "size: %Zd, nBits: %d\n", size, nBits);
	unsigned long seed = time(NULL);//123456ul;
	MPQRandomGenerator_t* mpqRGen = initMPQRandomGenerator(seed, nBits, nBits);

	mpq_t* interpPoints = (mpq_t*) malloc(sizeof(mpq_t)*(maxD + 1));
	mpq_t* interpVals = (mpq_t*) malloc(sizeof(mpq_t)*(maxD + 1));
	mpq_t* RJs = (mpq_t*) malloc(sizeof(mpq_t)*(maxD + 1));
	for (int i = 0; i < maxD+1; ++i) {
		mpq_init(interpPoints[i]);
	}
	for (int i = 0; i < maxD+1; ++i) {
		mpq_init(interpVals[i]);
	}
	for (int i = 0; i < maxD+1; ++i) {
		mpq_init(RJs[i]);
	}

	//get random initial point
	mpq_t* initPoint = (mpq_t*) malloc(sizeof(mpq_t)*(nvar));
	for (int i = 0; i < nvar; ++i) {
		mpq_init(initPoint[i]);
		MPQRandomGeneratorGetZ(mpqRGen, initPoint[i]);
	}
	mpq_t initVal;
	mpq_init(initVal);
	evalSMQPBlackBox(bb, initPoint, nvar, initVal);

	//get initial skeleton
	AltArr_t* skel = stage1Interpolate(degreeBounds[0], nvar, bb, initPoint, initVal, mpqRGen);
	
	mpq_clear(initVal);
	// char* syms[] = {"x", "y", "z", "s", "t", "u", "v", "w"};
	// printPoly_AA(stderr, skel, syms, nvar);
	// fprintf(stderr, "\nskel\n\n");

	//loop over remaining vars
	for (int k = 1; k < nvar; ++k) {
		//loop over degreeBound of that var
		// fprintf(stderr, "k, %d, degreeBounds %d\n", k, degreeBounds[k]);
		AltArr_t* dPolys [degreeBounds[k]]; 
		for (int j = 0; j < degreeBounds[k]; ++j) {
			//pick rj at random
			while (1) {
				int unique = 1;
				MPQRandomGeneratorGetQ(mpqRGen, RJs[j]);
				for (int idx = 0; idx < j; ++idx) {
					if (mpq_cmp(RJs[idx], RJs[j]) == 0) {
						unique = 0;
						break;
					}
				}
				if (unique && mpq_cmp(initPoint[k], RJs[j]) == 0) {
					continue;
				}
				if (unique) {
					break;
				}
			}
			dPolys[j] = getPolynomialPJ(bb, skel, k, degreeBounds[k], nvar, initPoint, RJs[j], mpqRGen);
			// fprintf(stderr, "k: %d, j: %d, poly:\n", k, j);
			// printPoly_AA(stderr, dPolys[j], syms, nvar);
			// fprintf(stderr, "\n\n");
		}

		AltArr_t* coefPolys[skel->size];
		for (int i = 0; i < skel->size; ++i) {
			for (int j = 0; j < degreeBounds[k]; ++j) {
				mpq_set(interpVals[j], dPolys[j]->elems[i].coef);
				mpq_set(interpPoints[j], RJs[j]);
			}
			mpq_set(interpVals[degreeBounds[k]], skel->elems[i].coef);
			mpq_set(interpPoints[degreeBounds[k]], initPoint[k]);
		

			// coefPolys[i] = univariateInterpolation_AA(degreeBounds[k], interpPoints, interpVals);
			coefPolys[i] = univarInterpolate_AA(interpPoints, interpVals, degreeBounds[k]+1);
			// fprintf(stderr, "Poly for monomial: %llx\n", skel->elems[i].degs);
			// printPoly_AA(stderr, coefPolys[i], syms, nvar);
			// fprintf(stderr, "\n\n");
		}

		AltArr_t* newskel = mergePolysWithSkel(skel, coefPolys, nvar, k);
		for (int i = 0; i < skel->size; ++i) {
			freePolynomial_AA(coefPolys[i]);
		}
		freePolynomial_AA(skel);
		skel = newskel;
		// fprintf(stderr, "New Skel: ");
		// printPoly_AA(stderr, skel, syms, nvar);
		// fprintf(stderr, "\n\n");
	}


	for (int i = 0; i < maxD; ++i) {
		mpq_clear(interpVals[i]);
	}
	for (int i = 0; i < maxD; ++i) {
		mpq_clear(interpPoints[i]);
	}
	for (int i = 0; i < maxD; ++i) {
		mpq_clear(RJs[i]);
	}

	free(interpVals);
	free(interpPoints);

	return skel;
}

/**
 * T is the bound on the number of terms in the polynomials.
 * If T < 0, then use T = max_i(degreeBounds[i])^nvar. 
 */
AltArr_t* multivariateInterpolateDeterministic_AA(SMQPBlackBox_t* bb, int* degreeBounds, int nvar, int T) {

	int maxD = 0;
	for (int i = 0; i < nvar; ++i) {
		maxD = maxD < degreeBounds[i] ? degreeBounds[i] : maxD;
	}

	if (T < 0) {
		T = (int) (pow(maxD, nvar) + 0.5);
	}
	int primeMatRows = T >= (maxD + 1) ? T : (maxD + 1);

	mpq_t* primeMat = getPrimeMatrix(primeMatRows, nvar);
	mpq_t* interpPoints = (mpq_t*) malloc(sizeof(mpq_t)*(maxD + 1));
	mpq_t* interpVals = (mpq_t*) malloc(sizeof(mpq_t)*(maxD + 1));
	mpq_t* RJs = (mpq_t*) malloc(sizeof(mpq_t)*(maxD + 1));
	for (int i = 0; i < maxD+1; ++i) {
		mpq_init(interpPoints[i]);
	}
	for (int i = 0; i < maxD+1; ++i) {
		mpq_init(interpVals[i]);
	}
	for (int i = 0; i < maxD+1; ++i) {
		mpq_init(RJs[i]);
	}

	//set initial point to (1,1,1,1,...);
	mpq_t initPoint[nvar];
	for (int i = 0; i < nvar; ++i) {
		mpq_init(initPoint[i]);
		mpq_set(initPoint[i], primeMat[i]);		
	}
	mpq_t bestPoint[nvar];
	for (int i = 0; i < nvar; ++i) {
		mpq_init(bestPoint[i]);
	}

	//get initial skeleton
	AltArr_t* skel = stage1InterpolatePrimes(degreeBounds[0], nvar, bb, primeMat, T);
	if (skel == NULL || skel->size == 0) {
		return NULL;
	}

	//loop over remaining vars; each iteration of this loop is a "stage"
	for (int k = 1; k < nvar; ++k) {
		
		AltArr_t* bestCoefPolys[skel->size]; 
		for (int i = 0; i < skel->size; ++i) {
			bestCoefPolys[i] = NULL;
		}

		//for each stage, we loop over 0..T as exponents of initPoint 
		//to ensure that our dense interpolation yield maximal monomials. 
		for (int exp = 0; exp <= T; ++exp) {
			//set init point to prime_i^exp;
			for (int i = k+1; i < nvar; ++i) {
				mpq_set(initPoint[i], primeMat[exp*nvar + i]);
			}
			
			//loop over degreeBound of that var
			//in constrast to probabilistic alg, we do not re-use the previous
			//point and must loop degreeBound + 1 times.
			AltArr_t* dPolys [degreeBounds[k]+1]; 
			for (int j = 0; j < degreeBounds[k]+1; ++j) {
				//pick rj to be the k'th prime number raised to the power of j
				mpq_set(RJs[j], primeMat[j*nvar + k]);
				dPolys[j] = getPolynomialPJPrimes(bb, skel, k, degreeBounds[k], nvar, primeMat, initPoint, RJs[j]);
			}

			int skip = 0;
			AltArr_t* coefPolys[skel->size];
			for (int i = 0; i < skel->size; ++i) {
				for (int j = 0; j < degreeBounds[k]+1; ++j) {
					mpq_set(interpVals[j], dPolys[j]->elems[i].coef);
					mpq_set(interpPoints[j], RJs[j]);
				}
				coefPolys[i] = univarInterpolate_AA(interpPoints, interpVals, degreeBounds[k]+1);
				if (coefPolys[i] == NULL || coefPolys[i]->size == 0) {
					skip = 1;
				}
			}

			//Determine if this interpolated polynomial is the best one yet
			//That is, will yield the next skeleton with larger number of monomials
			if (skip) {
				for (int i = 0; i < skel->size; ++i) {
					freePolynomial_AA(coefPolys[i]);
				}
			} else {
				skip = 1;
				for (int i = 0; i < skel->size; ++i) {
					if (bestCoefPolys[i] == NULL) {
						skip = 0;
						break;
					}
					if (bestCoefPolys[i]->size < coefPolys[i]->size) {
						skip = 0;
						break;
					}
				}
				if (!skip) {
					for (int i = 0; i < skel->size; ++i) {
						freePolynomial_AA(bestCoefPolys[i]);
						bestCoefPolys[i] = coefPolys[i];
						for (int i = k+1; i < nvar; ++i) {
							mpq_set(bestPoint[i], initPoint[i]);
						}
					}
				} else {
					for (int i = 0; i < skel->size; ++i) {
						freePolynomial_AA(coefPolys[i]);
					}
				}
			}
		}

		AltArr_t* newskel = mergePolysWithSkel(skel, bestCoefPolys, nvar, k);
		for (int i = 0; i < skel->size; ++i) {
			freePolynomial_AA(bestCoefPolys[i]);
		}
		freePolynomial_AA(skel);
		skel = newskel;
	}


	for (int i = 0; i < nvar; ++i) {
		mpq_clear(initPoint[i]);
	}
	for (int i = 0; i < nvar; ++i) {
		mpq_clear(bestPoint[i]);
	}
	for (int i = 0; i < maxD; ++i) {
		mpq_clear(interpVals[i]);
	}
	for (int i = 0; i < maxD; ++i) {
		mpq_clear(interpPoints[i]);
	}
	for (int i = 0; i < maxD; ++i) {
		mpq_clear(RJs[i]);
	}

	free(interpVals);
	free(interpPoints);
	free(RJs);

	freeMPQMat(primeMat, primeMatRows, nvar);

	return skel;
}


AltArr_t* univariateInterpolation_AA(int degreeBound, mpq_t* points, mpq_t* vals) {
	
	int t = degreeBound + 1;
	mpq_t* A = getMPQIdentity(t);
	mpq_t b[t];
	for (int i = 0; i < t; ++i) {
		mpq_init(b[i]);
	}

	//Set up linear equations using the same initialPoint, but varying
	//values for the variable we are trying to interpolate. 
	for (int i = 0; i < t; ++i) {
		mpq_set_si(A[i*t], 1l, 1l);
		mpq_set(A[i*t + 1], points[i]);
		for (int j = 2; j < t; ++j) {
			mpq_mul(A[i*t + j], A[i*t + 1], A[i*t + j -1]);
		}

		mpq_set(b[i], vals[i]);
	}

	mpq_t X[t];
	for (int i = 0; i < t; ++i) {
		mpq_init(X[i]);
	}
	long* rowList = NULL;

	long rank = solveMPQSystem(A, b, t, X, &rowList);
	free(rowList);

	for (int i = 0; i < t; ++i) {
		mpq_clear(b[i]);
	}
	freeMPQMat(A, t, t);

	AltArr_t* aa;
	if (rank < t) {
		aa = univarInterpolate_AA(points, vals, t);
	} else {
		aa = makePolynomial_AA(t, 1);
		int curIdx = 0;
		for (int i = t-1; i >= 0; --i) {
			if (mpq_sgn(X[i]) == 0) {
				continue;
			}
			mpq_init(aa->elems[curIdx].coef);
			mpq_set(aa->elems[curIdx].coef, X[i]);
			aa->elems[curIdx].degs = i;
			++curIdx;
		}
		aa->size = curIdx;
	}
	for (int i = 0; i < t; ++i) {
		mpq_clear(X[i]);
	}

	return aa;
}


