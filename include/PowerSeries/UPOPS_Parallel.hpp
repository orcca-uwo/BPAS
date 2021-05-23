#ifndef _UPOPS_PARA_H_
#define _UPOPS_PARA_H_

#include "UnivariatePolynomialOverPowerSeries.h"

#include <vector>

/**
 * Update to precision d an array of UPOPS which are the factors of
 * some Hensel Factorization.
 *
 * This function updates all factors simulaneously using
 * multiple threads, as specified by nthreads.
 * If nthreads is 1, computation is serial and uses the calling thread only.
 * If nthreads <= 0, the function estimates the most suitable number of threads to use.
 *
 * @param d: update to precision d
 * @param facts: the array of factors
 * @param nfacts: the number of factors
 * @param nthreads: the total number of threads to use.
 *
 */
void updateHenselFactsParallel_UPOPS(int d, Upops_t** facts, int nfacts, int nthreads);


/**
 * Update to precision d an array of UPOPS which are the factors of
 * some Hensel Factorization.
 *
 * This function updates all factors simulaneously using
 * multiple threads, as specified by threadsPerFactor, a vector of size nfacts.
 * If nfacts[i] = 0, facts[i] is updated by the same threads as assigned for facts[i+1].
 * threadsPerFactor[nfacts-1] should be positive to ensure the last factor gets updated.
 *
 * @param d: update to precision d
 * @param facts: the array of factors
 * @param nfacts: the number of factors
 * @param threadsPerFactor: the number of threads to use to update each factor; an vector of size nfacts.
 *
 */
void updateHenselFactsParallel_UPOPS(int d, Upops_t** facts, int nfacts, std::vector<int> threadsPerFactor);


/**
 * Update to precision d the upops p resulting from a Weierstrass preparation.
 * This automatically updates the corresponding alpha to precision d.
 *
 * This function utilizes up to nthreads threads to perform this update
 * in parallel.
 *
 * @param d: update p and alpha to pecision d
 * @param p: the p resulting from a Weierstrass preparation
 * @param nthreads: the maximum number of threads to use, including the calling threads
 */
void weierstrassUpdateParallel_UPOPS(int d, Upops_t* p, int nthreads);


#endif
