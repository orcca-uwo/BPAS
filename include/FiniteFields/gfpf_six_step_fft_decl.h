/* gfpf_six_step_fft.h
 *****************************************
 * an implementation of big prime field fft which relies
 * entirely on radix-based representation of data. Assuming
 * that prime is of the form p=r^k+1, an input vector X
 * of size N can be viewed as a row-major layout of an Nxk
 * matrix.
 * An element X_i is stored as a vector of size k.
 * X_i= (a_0, a_1, ..., a_{k-1}) in base r
 * ie, X_i = a_0 + a_1 .r + ... + a_{k-1}.r^{k-1}.
 *****************************************
 */

#ifndef GFPF_SIX_STEP_FFT_DECL_H_
#define GFPF_SIX_STEP_FFT_DECL_H_

#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "gfpf_arithmetic_decl.h"
#include "fourier_primes_u64.h"
#include "srgfn_primes.h"

/**************************************/

#define DFT_GRAIN_FACTOR (2)
#define CONVOLUTION_GRAIN_FACTOR (2)
#define FFT_BIG_ELEMENTS_PERMUTATION_BLOCKSIZE (32)
#define PERMUTATION_BLOCKSIZE FFT_BIG_ELEMENTS_PERMUTATION_BLOCKSIZE
#define SWAP_BLOCK_SIZE (4)

/**************************************/

#ifndef VERBOSE
#define VERBOSE 0
#endif
/**************************************/

#ifndef VERIFICATION_ENABLED
#define VERIFICATION_ENABLED 0
#endif

/**************************************/
static int supported_gfpf_DFT_K[4] =
  { 8, 16, 32, 64 };

/**************************************/

static void __inline__
DFT2_big_elements (usfixn64* a0, usfixn64* a1, const int k, const usfixn64 r);

/**************************************/

//older implementation of DFT2_big_elements
void
DFT2_big_elements_v0 (usfixn64* a0, usfixn64* a1, int k, usfixn64 r);

/**************************************/

static void __inline__
swap_big_elements (usfixn64* x, usfixn64* y, const int k);

/**************************************/

//older implementation of swap_big_elements
static void __inline__
swap_big_elements_v0 (usfixn64* a, usfixn64* b, int k);

/**************************************/

void
DFT_4_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		    int compute_inverse);

/**************************************/

void
DFT_8_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		    int compute_inverse);

/**************************************/

void
DFT_16_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse);

/**************************************/

void
DFT_32_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse);

/**************************************/

void
DFT_64_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse);

/**************************************/

void
DFT_128_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse);

/**************************************/

void
DFT_256_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse);

/**************************************/

/** compute L_{m}{n} for large coefficients;
 * equivalently, computing transpose {n}{m}
 */
void
stride_permutation_big_elements (usfixn64* A, usfixn64 *B, int m, int n,
				 int n_permutations, int coefficient_size);

/**************************************/

void
mult_pow_R_gmp (mpz_t x_zz, const int s, const usfixn64 k, const usfixn64 r,
		const mpz_t p_zz);

/**************************************/

/** computing twiddle T{m}{n} == D{n}{m}
 * computing D_{K}{J=K^b}
 * part of twiddle computation is done by multPowR;
 * the rest is done by a direct multiplication method
 * multiplication methods used:
 *  1. FFT-based-multiplication
 *  2. gmp-based multiplication
 */
void
twiddle_big_elements (usfixn64* vector, int K, int J, int n_permutations,
		      usfixn64* precomputed_pow_omega_vec, int k, usfixn64 r,
		      mpz_t p, int fft_based_mult_enabled, int compute_inverse,
		      float *t_mult_values, int* n_twiddle_mult);

/**************************************/

/** precompute powers of omega up to (omega)^{K^b}
 * to be used at level b+1
 * array "usfixn64* precomputed_pow_omega_vec" should be be
 * already allocated and passed as a pointer.
 */
void
precompute_sequential_powers_of_omega_big_elements (
    usfixn64* precomputed_pow_omega_vec, int K, int b, usfixn64 radix,
    mpz_t omega, mpz_t p);

/**************************************/

/** compute DFT_{N=K^e} for p=radix^{K/2}+1;
 * omega is omega at the highest level such that (omega^N)=1 mod p
 */
void
DFT_general_big_elements (usfixn64* vector_in, int K, int e, usfixn64 radix,
			  mpz_t p, usfixn64* omega, int fft_based_mult_enabled,
			  int compute_inverse, int verbose);

/**************************************/

int
convolution_mult (usfixn64* x_data, usfixn64* y_data, int n, srgfn_prime *P,
		  mpz_t p_zz, int fft_based_mult_enabled);

/**************************************/

void
test_DFT2_big_elements (usfixn64* x, int K, int e, usfixn64 radix);

/**************************************/

void
test_mult_pow_R_big_elements (usfixn64* x, int K, int e, usfixn64 radix);

/**************************************/

void
test_swap_big_elements (usfixn64* x, int K, int e, usfixn64 radix);

/**************************************/

#endif
