
#ifndef _RANDOM_HELPERS_H_
#define _RANDOM_HELPERS_H_

#include <stdlib.h> 
#include <time.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

static time_t seedRand() {
	static int initRand = 0;
	static time_t t;
	if (!initRand) {
		t = time(NULL);
		srand(t);
		initRand = 1;
	}
	return t;
}

static void rand_mpz_t(unsigned long int coefBound, int includeNeg, mpz_t mpzVal) {
	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = seedRand();

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		initRand = 1;
	}
	
	mpz_init(mpzVal);
	while(mpz_sgn(mpzVal) == 0) {
		mpz_urandomb(mpzVal, R_STATE, coefBound);
	}
	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzVal, mpzVal);
		// coef_l *= -1;
	}
}

static inline void rand_mpq_t(unsigned long int coefBound, int includeNeg, mpq_t mpqVal) {
	mpz_t mpzNum;
	mpz_t mpzDen;
	rand_mpz_t(coefBound,includeNeg,mpzNum);
	rand_mpz_t(coefBound,includeNeg,mpzDen);
	
	mpq_init(mpqVal);
	
	mpz_set(mpq_numref(mpqVal), mpzNum);
	mpz_set(mpq_denref(mpqVal), mpzDen);
	mpq_canonicalize(mpqVal);
}

/**
  * Generate a random value in a inclusive range
  *
  * @param low: lower bound of the range
  * @param high: upper bound of the range
  * @param numElems: number of values in the range;
  **/
static inline int randValInRange(int low,int high) {
    return rand()%(high-low+1)+low;     
}

static inline void rand_mpz_vec(unsigned long int coefBound, int includeNeg, mpz_t* mpzVec, unsigned int n) {
	for (unsigned int i = 0; i < n; ++i) {
		rand_mpz_t(coefBound, includeNeg, mpzVec[i]);
	}
}


#ifdef __cplusplus
}
#endif

#endif
