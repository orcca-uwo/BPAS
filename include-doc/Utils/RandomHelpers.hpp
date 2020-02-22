#ifndef _RANDOM_HELPERS_H_
#define _RANDOM_HELPERS_H_

//#include "../ring.h"
//#include "../RingPolynomial/upolynomial.h"
//#include "../Ring/RationalNumber.hpp"
#include <assert.h>
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include <vector>

static time_t seedRand() {
	static int initRand = 0;
	static time_t t = time(NULL);
	if (!initRand) {
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

static void rand_mpq_t(unsigned long int coefBound, int includeNeg, mpq_t mpqVal) {
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
  * Generate a random value in a range
  *
  * @param low: lower bound of the range
  * @param high: upper bound of the range
  * @param numElems: number of values in the range;
  **/
static int randValInRange(int low,int high) {
    return rand()%(high-low+1)+low;     
}

/**
  * Generate a fixed number of randomly selected values in a range without duplicates
  *
  * @param low: lower bound of the range
  * @param high: upper bound of the range
  * @param numElems: number of values in the range;
  **/
static std::vector<int> randValsInRange(int low,int high,int numElems) {
	assert (numElems <= (high-low+1));
	int index[numElems]={0},n;
	std::vector<int> out;
	bool check(false);
	
	for (auto i=0; i<numElems; ++i) {
		if (i==0)
			index[i] = randValInRange(low,high);
		else {
			while (check==false) {
				check = true;
				n = randValInRange(low,high);
				for (auto j=0; j<i; ++j) {
					if (index[j] == n)
						check = false;
				}
			}
			index[i] = n;
			check = false;
		}
	}
	for (auto i=0; i<numElems; ++i) {
		out.push_back(index[i]);
	}
	return out;
}

#endif
