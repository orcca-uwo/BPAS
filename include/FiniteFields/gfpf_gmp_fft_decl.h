/* big_fft_gmp.h
 *****************************************
 * an implementation of big prime field fft
 * which relies entirely on GMP arithmetic.
 *****************************************
 */

#ifndef GFPF_GMP_FFT_DECL_H_
#define GFPF_GMP_FFT_DECL_H_

#include <math.h>
#include <stdlib.h>
#include <gmp.h>
#include "gfpf_types.h"
#include "gfpf_gmp_tools_decl.h"
#include "cpu_timer.h"

#define BLOCKSIZE 32

/**************************************/
#ifndef PROFILING_GMP_ENABLED
#define PROFILING_GMP_ENABLED 0

static int global_gmp_profiling_ongoing = 0;

static int n_dft2_gmp_called = 0;
static int n_swap_gmp_called = 0;
static int n_DFTK_gmp_called = 0;
//int n_mult_pow_R_called = 0;

static int n_stride_permutation_gmp_called = 0;
static int stride_permutation_gmp_called_dims[256][2];

static int n_twiddle_gmp_called = 0;
static int twiddle_gmp_called_dims[256][2];

static int n_basecase_mult_gmp_called = 0;
#endif

/**************************************/

void
add_gmp (mpz_t X, mpz_t Y, mpz_t P);

/**************************************/

void
sub_gmp (mpz_t X, mpz_t Y, mpz_t P);

/**************************************/

void
mult_gmp (mpz_t X, mpz_t Y, mpz_t P);

/**************************************/

void
mult_aux_gmp (mpz_t X, const mpz_t Y, mpz_t aux, const mpz_t P);

/**************************************/

void
convolution_mult_gmp (mpz_t *X, mpz_t *Y, int n, mpz_t P);

/**************************************/
void
mult_gmp_big_elements (usfixn64 *x, usfixn64 *y, const int k, const usfixn64 r,
		       const mpz_t prime);

/**************************************/

void
swap_gmp (mpz_t a, mpz_t b);

/**************************************/

void
swap_aux_gmp (mpz_t a, mpz_t b, mpz_t tmp);

/**************************************/

void
DFT2_gmp (mpz_t a0, mpz_t a1, mpz_t prime);

/**************************************/

void
DFT2_aux_gmp (mpz_t a0, mpz_t a1, mpz_t aux, const mpz_t prime);

/**************************************/

void
DFT_8_gmp (mpz_t* a, mpz_t omega, mpz_t prime, mpz_t * precompute_pow_omega);

/**************************************/

void
DFT_8_aux_gmp (mpz_t* a, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega);

/**************************************/

void
DFT_16_gmp (mpz_t* a, mpz_t omega, mpz_t prime, mpz_t * omega_pow);

/**************************************/

void
DFT_16_aux_gmp (mpz_t* a, mpz_t aux, mpz_t omega, const mpz_t prime, const mpz_t * omega_pow);

/**************************************/

void
DFT_32_gmp (mpz_t* A, mpz_t omega, mpz_t prime, mpz_t * omega_pow);

/**************************************/

void
DFT_32_aux_gmp (mpz_t* A, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * omega_pow);

/**************************************/

void
DFT_64_gmp (mpz_t* A, mpz_t omega, mpz_t prime, mpz_t * precompute_pow_omega);

/**************************************/

void
DFT_64_aux_gmp (mpz_t* A, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega);

/**************************************/

void
DFT_128_aux_gmp (mpz_t* A, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega);

/**************************************/

void
DFT_256_aux_gmp (mpz_t* A, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega);

/**************************************/

void
stride_permutation_gmp (int n_permutations, mpz_t* A, int m, int n);

/**************************************/

void
stride_permutation_aux_gmp (int n_permutations, mpz_t* A, mpz_t* B, int m, int n);

/**************************************/

//precompute powers of omega up to (omega)^;
//to be used at level b+1
//sfixn64* precomputed_pow_omega_vec should be be already allocated!
//just pass a pointer.
void
precompute_sequential_powers_of_omega_gmp (mpz_t* precomputed_pow_omega_vec,
					   const int K, const int b,
					   const mpz_t omega, const mpz_t p);

/**************************************/

void
mult_gmp_for_profiling_only (mpz_t X, mpz_t Y, mpz_t omega, mpz_t P);

/**************************************/

void
twiddle_gmp (int n_permutations, mpz_t* vector, int K, int J, const mpz_t * omega_w,
	     const mpz_t *omega_base_precomputed_vec, const mpz_t prime,
	     float * t_twiddle_values, int * n_twiddle_mult);

/**************************************/

void
twiddle_aux_gmp (int n_permutations, mpz_t* vector, mpz_t* aux, int K, int J,
                 const mpz_t * omega_w, const mpz_t *omega_base_precomputed_vec,
                 const mpz_t prime, float * t_twiddle_values, int * n_twiddle_mult);

/**************************************/

void
DFT_general_gmp (mpz_t* vector, int K, int e, const mpz_t omega_zz, mpz_t prime,
		 int verbose);

/**************************************/
void
DFT_general_aux_memory_gmp (mpz_t* vector, int K, int e, const mpz_t omega_zz, mpz_t prime,
		 int verbose);

/**************************************/

#endif
