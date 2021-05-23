

#ifndef _DUSP_PARALLEL_H_
#define _DUSP_PARALLEL_H_


#include <gmpxx.h>
#include "ModularPolynomial/DUSP_Support.h"
#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "Utils/Parallel/ExecutorThreadPool.hpp"


///////////////////
namespace DUSP {
///////////////////


//evaluate and compute univar subresultant images for for k1 <= k < k2
//k1: inclusive starting index
//k2: exclusive ending index
//n: total number of images to be computed, needed for 1-D array access of 2D data
//n2: a bound on the size of the univar subresultant chain
//PPtr: the prime pointer for the current prime field
//gpd: partial degrees of the first input poly
//fpd: partial degrees of the second input poly
//t: evaluation points, used later for interp
//g: the first input poly, stored as a dense 2D array
//f: the second input poly stored as a dense 2D array
//S: An array of size n of duspolysA_t* for storing the results
void evalAndGetUnivarImage_SubRes_spX(int k1, int k2, int n, int n2, const Prime_ptr* Pptr, const polysize_t* gpd, const polysize_t* fpd, const elem_t* g, const elem_t* f, elem_t* t, duspolysA_t** S);

/**
 * Bivariate Subresultant Chain  Z[y<x] but polys sorted w.r.t y,x
 *
 * @param g a bivariate polynomial (elem_t*) Z[y<x]
 * @param gpd 2Dim partial degrees of polynomial g s.t. gpd[0]:deg_y and gpd[1]:deg_x
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f s.t. fpd[0]:deg_y and fpd[1]:deg_x
 * @param subres returned subresultant chain modular pr (biSybresPr_t)
 * @param Pptr small prime pointer
 * @param nthreads the maximum number of threads to use; if -1 try to determine best number of threads based on input size
 */
void biSylvSubResultantInForm_parallel_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd,
	biSubresPr_t** subres, const Prime_ptr* Pptr,
	const std::vector<ExecutorThreadPool::threadID>& workers);

/**
 * FFT-based Bivariate Subresultant Chain  Z[y<x] but polys sorted w.r.t y,x
 *
 * @param g a bivariate polynomial (elem_t*) Z[y<x]
 * @param gpd 2Dim partial degrees of polynomial g s.t. gpd[0]:deg_y and gpd[1]:deg_x
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f s.t. fpd[0]:deg_y and fpd[1]:deg_x
 * @param n a power of two s.t. n > B_y = gpd[1]*fpd[0] + gpd[0]*fpd[1] + 1
 * @param K
 * @param e
 * @param w a good primitive root of unity of order m in  Z_p
 * @param w_inv inverse of w mod p
 * @param n_inv inverse of n mod p
 * @param subres returned subresultant chain modular pr
 * @param Pptr small prime pointer
 */
int bivarModularSubResInForm_withFFT_parallel_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr,
											const std::vector<ExecutorThreadPool::threadID>& workers);
int bivarModularHGCDSubResInForm_withFFT_parallel_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr, 
										   const std::vector<ExecutorThreadPool::threadID>& workers);

int biModularSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime, int maxThreads = -1);
int biModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime, int maxThreads = -1);

// int hgcdBiModularSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size, int maxThreads = -1);
int hgcdBiModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size, int maxThreads = -1);

int ParallelSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, long min_mdeg, long max_pdeg, AltArrsZ_t** subres, int* chain_size, int maxThreads = -1);

///////////////////
} // Namespace DUSP
///////////////////


#endif
