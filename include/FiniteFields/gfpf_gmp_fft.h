#ifndef GFPF_GMP_FFT_H_
#define GFPF_GMP_FFT_H_
/* big_fft_gmp.h
 *****************************************
 * an implementation of big prime field fft
 * which relies entirely on GMP arithmetic.
 *****************************************
 */

#include "gfpf_gmp_fft_decl.h"

#ifndef GMP_SERIAL_ENABLED
#define GMP_SERIAL_ENABLED 0
#endif

void
inc_gmp_profiling_counter (int *counter)
{
  if (global_gmp_profiling_ongoing == 0)
    *counter++;
}

/**************************************/

void
add_gmp (mpz_t X, mpz_t Y, mpz_t P)
{
  mpz_add (X, X, Y);
  mpz_mod (X, X, P);
}

/**************************************/

void
sub_gmp (mpz_t X, mpz_t Y, mpz_t P)
{
  mpz_sub (X, X, Y);
  mpz_mod (X, X, P);
}

/**************************************/

void
mult_gmp (mpz_t X, mpz_t Y, mpz_t P)
{
#if PROFILING_GMP_ENABLED == 1
  inc_gmp_profiling_counter (&n_basecase_mult_gmp_called);
#endif

  mpz_mul (X, X, Y);
  mpz_mod (X, X, P);
}

/**************************************/

void __inline__
mult_aux_gmp (mpz_t X, const mpz_t Y, mpz_t aux, const mpz_t P)
{
#if PROFILING_GMP_ENABLED == 1
  inc_gmp_profiling_counter (&n_basecase_mult_gmp_called);
#endif

  mpz_mul (aux, X, Y);
//  mpz_mod (X, aux, P);
  mpz_tdiv_r(X, aux, P);
}


/**************************************/

void
convolution_mult_gmp (mpz_t *X, mpz_t *Y, int n, mpz_t P)
{
#if PROFILING_GMP_ENABLED == 1
  inc_gmp_profiling_counter (&n_basecase_mult_gmp_called);
#endif

  cpu_timer t_convolution_mult;
  timer_record_start (&t_convolution_mult);

  for (int i = 0; i < n; i++)
    {
      mpz_mul (X[i], X[i], Y[i]);
      mpz_mod (X[i], X[i], P);
    }
  timer_record_stop (&t_convolution_mult);
  timer_get_elapsed_time (&t_convolution_mult, "t_convolution_mult_gmp", 1);
}

/**************************************/

void
mult_gmp_big_elements (usfixn64 *x, usfixn64 *y, const int k, const usfixn64 r,
		       const mpz_t prime)
{
  mpz_t X, Y;
  mpz_inits (X, Y, NULL);
  u64vector_to_radix_based_bigint_gmp (X, x, r, k);
  u64vector_to_radix_based_bigint_gmp (Y, y, r, k);

//  X = (X * Y) % prime;
  mpz_mul (X, X, Y);
  mpz_mod (X, X, prime);

  bigint_to_u64vector_radixbased_gmp (x, X, r, k);
  mpz_clears (X, Y, NULL);
}

/**************************************/

void
swap_gmp (mpz_t a, mpz_t b)
{
  mpz_t tmp;
  mpz_init (tmp);
  mpz_set (tmp, a);
  mpz_set (a, b);
  mpz_set (b, tmp);
}

/**************************************/

void __inline__
swap_aux_gmp (mpz_t a, mpz_t b, mpz_t tmp)
{
//  mpz_t tmp;
//  mpz_init (tmp);
  mpz_set (tmp, a);
  mpz_set (a, b);
  mpz_set (b, tmp);
}

/**************************************/

void
DFT2_gmp (mpz_t a0, mpz_t a1, mpz_t prime)
{
#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif
  mpz_t sub;
  mpz_init_set (sub, a0);
  sub_gmp (sub, a1, prime);
  add_gmp (a0, a1, prime);
  mpz_set (a1, sub);
  mpz_clear (sub);
}

/**************************************/

void __inline__
DFT2_aux_gmp (mpz_t a0, mpz_t a1, mpz_t t_zz, const mpz_t prime)
{
#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  mpz_set (t_zz, a0);
  mpz_add (a0, a0, a1);
  mpz_sub (a1, t_zz, a1);

  //  mpz_mod (a0, a0, prime);
  if (mpz_cmp (a0, prime) >= 0)
  mpz_sub (a0, a0, prime);

  //  mpz_mod (a1, a1, prime);
  if (mpz_cmp_ui(a1, 0) < 0)
  mpz_add (a1, a1, prime);

}

/**************************************/

void
DFT_8_gmp (mpz_t* a, mpz_t omega, mpz_t prime, mpz_t * precompute_pow_omega)
{

#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_DFTK_gmp_called);
#endif

  DFT2_gmp (a[0], a[4], prime); // dft on permutated indexes
  DFT2_gmp (a[2], a[6], prime);
  DFT2_gmp (a[1], a[5], prime);
  DFT2_gmp (a[3], a[7], prime);

  mult_gmp (a[6], precompute_pow_omega[2], prime);   //twiddle
  mult_gmp (a[7], precompute_pow_omega[2], prime);

  DFT2_gmp (a[0], a[2], prime); // dft on permutated indexes
  DFT2_gmp (a[4], a[6], prime);
  DFT2_gmp (a[1], a[3], prime);
  DFT2_gmp (a[5], a[7], prime);

  mult_gmp (a[5], precompute_pow_omega[1], prime); // twiddle
  mult_gmp (a[3], precompute_pow_omega[2], prime);
  mult_gmp (a[7], precompute_pow_omega[3], prime);

  DFT2_gmp (a[0], a[1], prime);   // dft on permutated indexes
  DFT2_gmp (a[4], a[5], prime);
  DFT2_gmp (a[2], a[3], prime);
  DFT2_gmp (a[6], a[7], prime);

  // final permutation
  swap_gmp (a[1], a[4]);
  swap_gmp (a[3], a[6]);

}

/**************************************/

void
DFT_8_aux_gmp (mpz_t* a, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega)
{

#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_DFTK_gmp_called);
#endif

  DFT2_aux_gmp (a[0], a[4], aux, prime); // dft on permutated indexes
  DFT2_aux_gmp (a[1], a[5], aux, prime);
  DFT2_aux_gmp (a[2], a[6], aux, prime);
  DFT2_aux_gmp (a[3], a[7], aux, prime);

  mult_aux_gmp (a[6], precompute_pow_omega[2], aux, prime);   //twiddle
  mult_aux_gmp (a[7], precompute_pow_omega[2], aux, prime);

  DFT2_aux_gmp (a[0], a[2], aux, prime); // dft on permutated indexes
  DFT2_aux_gmp (a[1], a[3], aux, prime);
  DFT2_aux_gmp (a[4], a[6], aux, prime);
  DFT2_aux_gmp (a[5], a[7], aux, prime);

  mult_aux_gmp (a[5], precompute_pow_omega[1], aux, prime); // twiddle
  mult_aux_gmp (a[3], precompute_pow_omega[2], aux, prime);
  mult_aux_gmp (a[7], precompute_pow_omega[3], aux, prime);

  DFT2_aux_gmp (a[0], a[1], aux, prime);   // dft on permutated indexes
  DFT2_aux_gmp (a[4], a[5], aux, prime);
  DFT2_aux_gmp (a[2], a[3], aux, prime);
  DFT2_aux_gmp (a[6], a[7], aux, prime);

  // final permutation
  swap_aux_gmp (a[1], a[4], aux);
  swap_aux_gmp (a[3], a[6], aux);

}

/**************************************/

void
DFT_16_gmp (mpz_t* a, mpz_t omega, mpz_t prime, mpz_t * omega_pow)
{

#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_gmp (a[0], a[8], prime);
  DFT2_gmp (a[4], a[12], prime);
  DFT2_gmp (a[2], a[10], prime);
  DFT2_gmp (a[6], a[14], prime);
  DFT2_gmp (a[1], a[9], prime);
  DFT2_gmp (a[5], a[13], prime);
  DFT2_gmp (a[3], a[11], prime);
  DFT2_gmp (a[7], a[15], prime);

  mult_gmp (a[12], omega_pow[4], prime);
  mult_gmp (a[14], omega_pow[4], prime);
  mult_gmp (a[13], omega_pow[4], prime);
  mult_gmp (a[15], omega_pow[4], prime);

  DFT2_gmp (a[0], a[4], prime);
  DFT2_gmp (a[8], a[12], prime);
  DFT2_gmp (a[2], a[6], prime);
  DFT2_gmp (a[10], a[14], prime);
  DFT2_gmp (a[1], a[5], prime);
  DFT2_gmp (a[9], a[13], prime);
  DFT2_gmp (a[3], a[7], prime);
  DFT2_gmp (a[11], a[15], prime);

  mult_gmp (a[10], omega_pow[2], prime);
  mult_gmp (a[6], omega_pow[4], prime);
  mult_gmp (a[14], omega_pow[6], prime);

  mult_gmp (a[11], omega_pow[2], prime);
  mult_gmp (a[7], omega_pow[4], prime);
  mult_gmp (a[15], omega_pow[6], prime);

  DFT2_gmp (a[0], a[2], prime);
  DFT2_gmp (a[8], a[10], prime);
  DFT2_gmp (a[4], a[6], prime);
  DFT2_gmp (a[12], a[14], prime);
  DFT2_gmp (a[1], a[3], prime);
  DFT2_gmp (a[9], a[11], prime);
  DFT2_gmp (a[5], a[7], prime);
  DFT2_gmp (a[13], a[15], prime);

  mult_gmp (a[9], omega_pow[1], prime);
  mult_gmp (a[5], omega_pow[2], prime);
  mult_gmp (a[13], omega_pow[3], prime);
  mult_gmp (a[3], omega_pow[4], prime);
  mult_gmp (a[11], omega_pow[5], prime);
  mult_gmp (a[7], omega_pow[6], prime);
  mult_gmp (a[15], omega_pow[7], prime);

  DFT2_gmp (a[0], a[1], prime);
  DFT2_gmp (a[8], a[9], prime);
  DFT2_gmp (a[4], a[5], prime);
  DFT2_gmp (a[12], a[13], prime);
  DFT2_gmp (a[2], a[3], prime);
  DFT2_gmp (a[10], a[11], prime);
  DFT2_gmp (a[6], a[7], prime);
  DFT2_gmp (a[14], a[15], prime);

  swap_gmp (a[1], a[8]);
  swap_gmp (a[2], a[4]);
  swap_gmp (a[3], a[12]);
  swap_gmp (a[5], a[10]);
  swap_gmp (a[7], a[14]);
  swap_gmp (a[11], a[13]);

}

/**************************************/

void
DFT_16_aux_gmp (mpz_t* a, mpz_t aux, mpz_t omega, const mpz_t prime, const mpz_t * omega_pow)
{

#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_aux_gmp (a[0], a[8], aux, prime);
  DFT2_aux_gmp (a[4], a[12], aux, prime);
  DFT2_aux_gmp (a[2], a[10], aux, prime);
  DFT2_aux_gmp (a[6], a[14], aux, prime);
  DFT2_aux_gmp (a[1], a[9], aux, prime);
  DFT2_aux_gmp (a[5], a[13], aux, prime);
  DFT2_aux_gmp (a[3], a[11], aux, prime);
  DFT2_aux_gmp (a[7], a[15], aux, prime);

  mult_aux_gmp (a[12], omega_pow[4], aux, prime);
  mult_aux_gmp (a[14], omega_pow[4], aux, prime);
  mult_aux_gmp (a[13], omega_pow[4], aux, prime);
  mult_aux_gmp (a[15], omega_pow[4], aux, prime);

  DFT2_aux_gmp (a[0], a[4], aux, prime);
  DFT2_aux_gmp (a[8], a[12], aux, prime);
  DFT2_aux_gmp (a[2], a[6], aux, prime);
  DFT2_aux_gmp (a[10], a[14], aux, prime);
  DFT2_aux_gmp (a[1], a[5], aux, prime);
  DFT2_aux_gmp (a[9], a[13], aux, prime);
  DFT2_aux_gmp (a[3], a[7], aux, prime);
  DFT2_aux_gmp (a[11], a[15], aux, prime);

  mult_aux_gmp (a[10], omega_pow[2], aux, prime);
  mult_aux_gmp (a[6], omega_pow[4], aux, prime);
  mult_aux_gmp (a[14], omega_pow[6], aux, prime);

  mult_aux_gmp (a[11], omega_pow[2], aux, prime);
  mult_aux_gmp (a[7], omega_pow[4], aux, prime);
  mult_aux_gmp (a[15], omega_pow[6], aux, prime);

  DFT2_aux_gmp (a[0], a[2], aux, prime);
  DFT2_aux_gmp (a[8], a[10], aux, prime);
  DFT2_aux_gmp (a[4], a[6], aux, prime);
  DFT2_aux_gmp (a[12], a[14], aux, prime);
  DFT2_aux_gmp (a[1], a[3], aux, prime);
  DFT2_aux_gmp (a[9], a[11], aux, prime);
  DFT2_aux_gmp (a[5], a[7], aux, prime);
  DFT2_aux_gmp (a[13], a[15], aux, prime);

  mult_aux_gmp (a[9], omega_pow[1], aux, prime);
  mult_aux_gmp (a[5], omega_pow[2], aux, prime);
  mult_aux_gmp (a[13], omega_pow[3], aux, prime);
  mult_aux_gmp (a[3], omega_pow[4], aux, prime);
  mult_aux_gmp (a[11], omega_pow[5], aux, prime);
  mult_aux_gmp (a[7], omega_pow[6], aux, prime);
  mult_aux_gmp (a[15], omega_pow[7], aux, prime);

  DFT2_aux_gmp (a[0], a[1], aux, prime);
  DFT2_aux_gmp (a[8], a[9], aux, prime);
  DFT2_aux_gmp (a[4], a[5], aux, prime);
  DFT2_aux_gmp (a[12], a[13], aux, prime);
  DFT2_aux_gmp (a[2], a[3], aux, prime);
  DFT2_aux_gmp (a[10], a[11], aux, prime);
  DFT2_aux_gmp (a[6], a[7], aux, prime);
  DFT2_aux_gmp (a[14], a[15], aux, prime);

  swap_aux_gmp (a[1], a[8], aux);
  swap_aux_gmp (a[2], a[4], aux);
  swap_aux_gmp (a[3], a[12], aux);
  swap_aux_gmp (a[5], a[10], aux);
  swap_aux_gmp (a[7], a[14], aux);
  swap_aux_gmp (a[11], a[13], aux);

}

/**************************************/

void
DFT_32_gmp (mpz_t* A, mpz_t omega, mpz_t prime, mpz_t * omega_pow)
{

#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_gmp (A[0], A[16], prime);
  DFT2_gmp (A[8], A[24], prime);
  DFT2_gmp (A[4], A[20], prime);
  DFT2_gmp (A[12], A[28], prime);
  DFT2_gmp (A[2], A[18], prime);
  DFT2_gmp (A[10], A[26], prime);
  DFT2_gmp (A[6], A[22], prime);
  DFT2_gmp (A[14], A[30], prime);
  DFT2_gmp (A[1], A[17], prime);
  DFT2_gmp (A[9], A[25], prime);
  DFT2_gmp (A[5], A[21], prime);
  DFT2_gmp (A[13], A[29], prime);
  DFT2_gmp (A[3], A[19], prime);
  DFT2_gmp (A[11], A[27], prime);
  DFT2_gmp (A[7], A[23], prime);
  DFT2_gmp (A[15], A[31], prime);

  mult_gmp (A[24], omega_pow[8], prime);
  mult_gmp (A[28], omega_pow[8], prime);
  mult_gmp (A[26], omega_pow[8], prime);
  mult_gmp (A[30], omega_pow[8], prime);
  mult_gmp (A[25], omega_pow[8], prime);
  mult_gmp (A[29], omega_pow[8], prime);
  mult_gmp (A[27], omega_pow[8], prime);
  mult_gmp (A[31], omega_pow[8], prime);

  DFT2_gmp (A[0], A[8], prime);
  DFT2_gmp (A[16], A[24], prime);
  DFT2_gmp (A[4], A[12], prime);
  DFT2_gmp (A[20], A[28], prime);
  DFT2_gmp (A[2], A[10], prime);
  DFT2_gmp (A[18], A[26], prime);
  DFT2_gmp (A[6], A[14], prime);
  DFT2_gmp (A[22], A[30], prime);
  DFT2_gmp (A[1], A[9], prime);
  DFT2_gmp (A[17], A[25], prime);
  DFT2_gmp (A[5], A[13], prime);
  DFT2_gmp (A[21], A[29], prime);
  DFT2_gmp (A[3], A[11], prime);
  DFT2_gmp (A[19], A[27], prime);
  DFT2_gmp (A[7], A[15], prime);
  DFT2_gmp (A[23], A[31], prime);

  mult_gmp (A[20], omega_pow[4], prime);
  mult_gmp (A[12], omega_pow[8], prime);
  mult_gmp (A[28], omega_pow[12], prime);

  mult_gmp (A[22], omega_pow[4], prime);
  mult_gmp (A[14], omega_pow[8], prime);
  mult_gmp (A[30], omega_pow[12], prime);

  mult_gmp (A[21], omega_pow[4], prime);
  mult_gmp (A[13], omega_pow[8], prime);
  mult_gmp (A[29], omega_pow[12], prime);

  mult_gmp (A[23], omega_pow[4], prime);
  mult_gmp (A[15], omega_pow[8], prime);
  mult_gmp (A[31], omega_pow[12], prime);

  DFT2_gmp (A[0], A[4], prime);
  DFT2_gmp (A[16], A[20], prime);
  DFT2_gmp (A[8], A[12], prime);
  DFT2_gmp (A[24], A[28], prime);
  DFT2_gmp (A[2], A[6], prime);
  DFT2_gmp (A[18], A[22], prime);
  DFT2_gmp (A[10], A[14], prime);
  DFT2_gmp (A[26], A[30], prime);
  DFT2_gmp (A[1], A[5], prime);
  DFT2_gmp (A[17], A[21], prime);
  DFT2_gmp (A[9], A[13], prime);
  DFT2_gmp (A[25], A[29], prime);
  DFT2_gmp (A[3], A[7], prime);
  DFT2_gmp (A[19], A[23], prime);
  DFT2_gmp (A[11], A[15], prime);
  DFT2_gmp (A[27], A[31], prime);

  mult_gmp (A[18], omega_pow[2], prime);
  mult_gmp (A[10], omega_pow[4], prime);
  mult_gmp (A[26], omega_pow[6], prime);
  mult_gmp (A[6], omega_pow[8], prime);
  mult_gmp (A[22], omega_pow[10], prime);
  mult_gmp (A[14], omega_pow[12], prime);
  mult_gmp (A[30], omega_pow[14], prime);

  mult_gmp (A[19], omega_pow[2], prime);
  mult_gmp (A[11], omega_pow[4], prime);
  mult_gmp (A[27], omega_pow[6], prime);
  mult_gmp (A[7], omega_pow[8], prime);
  mult_gmp (A[23], omega_pow[10], prime);
  mult_gmp (A[15], omega_pow[12], prime);
  mult_gmp (A[31], omega_pow[14], prime);

  DFT2_gmp (A[0], A[2], prime);
  DFT2_gmp (A[16], A[18], prime);
  DFT2_gmp (A[8], A[10], prime);
  DFT2_gmp (A[24], A[26], prime);
  DFT2_gmp (A[4], A[6], prime);
  DFT2_gmp (A[20], A[22], prime);
  DFT2_gmp (A[12], A[14], prime);
  DFT2_gmp (A[28], A[30], prime);
  DFT2_gmp (A[1], A[3], prime);
  DFT2_gmp (A[17], A[19], prime);
  DFT2_gmp (A[9], A[11], prime);
  DFT2_gmp (A[25], A[27], prime);
  DFT2_gmp (A[5], A[7], prime);
  DFT2_gmp (A[21], A[23], prime);
  DFT2_gmp (A[13], A[15], prime);
  DFT2_gmp (A[29], A[31], prime);

  mult_gmp (A[17], omega_pow[1], prime);
  mult_gmp (A[9], omega_pow[2], prime);
  mult_gmp (A[25], omega_pow[3], prime);
  mult_gmp (A[5], omega_pow[4], prime);
  mult_gmp (A[21], omega_pow[5], prime);
  mult_gmp (A[13], omega_pow[6], prime);
  mult_gmp (A[29], omega_pow[7], prime);
  mult_gmp (A[3], omega_pow[8], prime);
  mult_gmp (A[19], omega_pow[9], prime);
  mult_gmp (A[11], omega_pow[10], prime);
  mult_gmp (A[27], omega_pow[11], prime);
  mult_gmp (A[7], omega_pow[12], prime);
  mult_gmp (A[23], omega_pow[13], prime);
  mult_gmp (A[15], omega_pow[14], prime);
  mult_gmp (A[31], omega_pow[15], prime);

  DFT2_gmp (A[0], A[1], prime);
  DFT2_gmp (A[16], A[17], prime);
  DFT2_gmp (A[8], A[9], prime);
  DFT2_gmp (A[24], A[25], prime);
  DFT2_gmp (A[4], A[5], prime);
  DFT2_gmp (A[20], A[21], prime);
  DFT2_gmp (A[12], A[13], prime);
  DFT2_gmp (A[28], A[29], prime);
  DFT2_gmp (A[2], A[3], prime);
  DFT2_gmp (A[18], A[19], prime);
  DFT2_gmp (A[10], A[11], prime);
  DFT2_gmp (A[26], A[27], prime);
  DFT2_gmp (A[6], A[7], prime);
  DFT2_gmp (A[22], A[23], prime);
  DFT2_gmp (A[14], A[15], prime);
  DFT2_gmp (A[30], A[31], prime);

  swap_gmp (A[1], A[16]);
  swap_gmp (A[2], A[8]);
  swap_gmp (A[3], A[24]);
  swap_gmp (A[5], A[20]);
  swap_gmp (A[6], A[12]);
  swap_gmp (A[7], A[28]);
  swap_gmp (A[9], A[18]);
  swap_gmp (A[11], A[26]);
  swap_gmp (A[13], A[22]);
  swap_gmp (A[15], A[30]);
  swap_gmp (A[19], A[25]);
  swap_gmp (A[23], A[29]);
}

/**************************************/

void
DFT_32_aux_gmp (mpz_t* A, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * omega_pow)
{

#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_aux_gmp (A[0], A[16], aux, prime);
  DFT2_aux_gmp (A[8], A[24], aux, prime);
  DFT2_aux_gmp (A[4], A[20], aux, prime);
  DFT2_aux_gmp (A[12], A[28], aux, prime);
  DFT2_aux_gmp (A[2], A[18], aux, prime);
  DFT2_aux_gmp (A[10], A[26], aux, prime);
  DFT2_aux_gmp (A[6], A[22], aux, prime);
  DFT2_aux_gmp (A[14], A[30], aux, prime);
  DFT2_aux_gmp (A[1], A[17], aux, prime);
  DFT2_aux_gmp (A[9], A[25], aux, prime);
  DFT2_aux_gmp (A[5], A[21], aux, prime);
  DFT2_aux_gmp (A[13], A[29], aux, prime);
  DFT2_aux_gmp (A[3], A[19], aux, prime);
  DFT2_aux_gmp (A[11], A[27], aux, prime);
  DFT2_aux_gmp (A[7], A[23], aux, prime);
  DFT2_aux_gmp (A[15], A[31], aux, prime);

  mult_aux_gmp (A[24], omega_pow[8], aux, prime);
  mult_aux_gmp (A[28], omega_pow[8], aux, prime);
  mult_aux_gmp (A[26], omega_pow[8], aux, prime);
  mult_aux_gmp (A[30], omega_pow[8], aux, prime);
  mult_aux_gmp (A[25], omega_pow[8], aux, prime);
  mult_aux_gmp (A[29], omega_pow[8], aux, prime);
  mult_aux_gmp (A[27], omega_pow[8], aux, prime);
  mult_aux_gmp (A[31], omega_pow[8], aux, prime);

  DFT2_aux_gmp (A[0], A[8], aux, prime);
  DFT2_aux_gmp (A[16], A[24], aux, prime);
  DFT2_aux_gmp (A[4], A[12], aux, prime);
  DFT2_aux_gmp (A[20], A[28], aux, prime);
  DFT2_aux_gmp (A[2], A[10], aux, prime);
  DFT2_aux_gmp (A[18], A[26], aux, prime);
  DFT2_aux_gmp (A[6], A[14], aux, prime);
  DFT2_aux_gmp (A[22], A[30], aux, prime);
  DFT2_aux_gmp (A[1], A[9], aux, prime);
  DFT2_aux_gmp (A[17], A[25], aux, prime);
  DFT2_aux_gmp (A[5], A[13], aux, prime);
  DFT2_aux_gmp (A[21], A[29], aux, prime);
  DFT2_aux_gmp (A[3], A[11], aux, prime);
  DFT2_aux_gmp (A[19], A[27], aux, prime);
  DFT2_aux_gmp (A[7], A[15], aux, prime);
  DFT2_aux_gmp (A[23], A[31], aux, prime);

  mult_aux_gmp (A[20], omega_pow[4], aux, prime);
  mult_aux_gmp (A[12], omega_pow[8], aux, prime);
  mult_aux_gmp (A[28], omega_pow[12], aux, prime);

  mult_aux_gmp (A[22], omega_pow[4], aux, prime);
  mult_aux_gmp (A[14], omega_pow[8], aux, prime);
  mult_aux_gmp (A[30], omega_pow[12], aux, prime);

  mult_aux_gmp (A[21], omega_pow[4], aux, prime);
  mult_aux_gmp (A[13], omega_pow[8], aux, prime);
  mult_aux_gmp (A[29], omega_pow[12], aux, prime);

  mult_aux_gmp (A[23], omega_pow[4], aux, prime);
  mult_aux_gmp (A[15], omega_pow[8], aux, prime);
  mult_aux_gmp (A[31], omega_pow[12], aux, prime);

  DFT2_aux_gmp (A[0], A[4], aux, prime);
  DFT2_aux_gmp (A[16], A[20], aux, prime);
  DFT2_aux_gmp (A[8], A[12], aux, prime);
  DFT2_aux_gmp (A[24], A[28], aux, prime);
  DFT2_aux_gmp (A[2], A[6], aux, prime);
  DFT2_aux_gmp (A[18], A[22], aux, prime);
  DFT2_aux_gmp (A[10], A[14], aux, prime);
  DFT2_aux_gmp (A[26], A[30], aux, prime);
  DFT2_aux_gmp (A[1], A[5], aux, prime);
  DFT2_aux_gmp (A[17], A[21], aux, prime);
  DFT2_aux_gmp (A[9], A[13], aux, prime);
  DFT2_aux_gmp (A[25], A[29], aux, prime);
  DFT2_aux_gmp (A[3], A[7], aux, prime);
  DFT2_aux_gmp (A[19], A[23], aux, prime);
  DFT2_aux_gmp (A[11], A[15], aux, prime);
  DFT2_aux_gmp (A[27], A[31], aux, prime);

  mult_aux_gmp (A[18], omega_pow[2], aux, prime);
  mult_aux_gmp (A[10], omega_pow[4], aux, prime);
  mult_aux_gmp (A[26], omega_pow[6], aux, prime);
  mult_aux_gmp (A[6], omega_pow[8], aux, prime);
  mult_aux_gmp (A[22], omega_pow[10], aux, prime);
  mult_aux_gmp (A[14], omega_pow[12], aux, prime);
  mult_aux_gmp (A[30], omega_pow[14], aux, prime);

  mult_aux_gmp (A[19], omega_pow[2], aux, prime);
  mult_aux_gmp (A[11], omega_pow[4], aux, prime);
  mult_aux_gmp (A[27], omega_pow[6], aux, prime);
  mult_aux_gmp (A[7], omega_pow[8], aux, prime);
  mult_aux_gmp (A[23], omega_pow[10], aux, prime);
  mult_aux_gmp (A[15], omega_pow[12], aux, prime);
  mult_aux_gmp (A[31], omega_pow[14], aux, prime);

  DFT2_aux_gmp (A[0], A[2], aux, prime);
  DFT2_aux_gmp (A[16], A[18], aux, prime);
  DFT2_aux_gmp (A[8], A[10], aux, prime);
  DFT2_aux_gmp (A[24], A[26], aux, prime);
  DFT2_aux_gmp (A[4], A[6], aux, prime);
  DFT2_aux_gmp (A[20], A[22], aux, prime);
  DFT2_aux_gmp (A[12], A[14], aux, prime);
  DFT2_aux_gmp (A[28], A[30], aux, prime);
  DFT2_aux_gmp (A[1], A[3], aux, prime);
  DFT2_aux_gmp (A[17], A[19], aux, prime);
  DFT2_aux_gmp (A[9], A[11], aux, prime);
  DFT2_aux_gmp (A[25], A[27], aux, prime);
  DFT2_aux_gmp (A[5], A[7], aux, prime);
  DFT2_aux_gmp (A[21], A[23], aux, prime);
  DFT2_aux_gmp (A[13], A[15], aux, prime);
  DFT2_aux_gmp (A[29], A[31], aux, prime);

  mult_aux_gmp (A[17], omega_pow[1], aux, prime);
  mult_aux_gmp (A[9], omega_pow[2], aux, prime);
  mult_aux_gmp (A[25], omega_pow[3], aux, prime);
  mult_aux_gmp (A[5], omega_pow[4], aux, prime);
  mult_aux_gmp (A[21], omega_pow[5], aux, prime);
  mult_aux_gmp (A[13], omega_pow[6], aux, prime);
  mult_aux_gmp (A[29], omega_pow[7], aux, prime);
  mult_aux_gmp (A[3], omega_pow[8], aux, prime);
  mult_aux_gmp (A[19], omega_pow[9], aux, prime);
  mult_aux_gmp (A[11], omega_pow[10], aux, prime);
  mult_aux_gmp (A[27], omega_pow[11], aux, prime);
  mult_aux_gmp (A[7], omega_pow[12], aux, prime);
  mult_aux_gmp (A[23], omega_pow[13], aux, prime);
  mult_aux_gmp (A[15], omega_pow[14], aux, prime);
  mult_aux_gmp (A[31], omega_pow[15], aux, prime);

  DFT2_aux_gmp (A[0], A[1], aux, prime);
  DFT2_aux_gmp (A[16], A[17], aux, prime);
  DFT2_aux_gmp (A[8], A[9], aux, prime);
  DFT2_aux_gmp (A[24], A[25], aux, prime);
  DFT2_aux_gmp (A[4], A[5], aux, prime);
  DFT2_aux_gmp (A[20], A[21], aux, prime);
  DFT2_aux_gmp (A[12], A[13], aux, prime);
  DFT2_aux_gmp (A[28], A[29], aux, prime);
  DFT2_aux_gmp (A[2], A[3], aux, prime);
  DFT2_aux_gmp (A[18], A[19], aux, prime);
  DFT2_aux_gmp (A[10], A[11], aux, prime);
  DFT2_aux_gmp (A[26], A[27], aux, prime);
  DFT2_aux_gmp (A[6], A[7], aux, prime);
  DFT2_aux_gmp (A[22], A[23], aux, prime);
  DFT2_aux_gmp (A[14], A[15], aux, prime);
  DFT2_aux_gmp (A[30], A[31], aux, prime);

  swap_aux_gmp (A[1], A[16], aux);
  swap_aux_gmp (A[2], A[8], aux);
  swap_aux_gmp (A[3], A[24], aux);
  swap_aux_gmp (A[5], A[20], aux);
  swap_aux_gmp (A[6], A[12], aux);
  swap_aux_gmp (A[7], A[28], aux);
  swap_aux_gmp (A[9], A[18], aux);
  swap_aux_gmp (A[11], A[26], aux);
  swap_aux_gmp (A[13], A[22], aux);
  swap_aux_gmp (A[15], A[30], aux);
  swap_aux_gmp (A[19], A[25], aux);
  swap_aux_gmp (A[23], A[29], aux);
}

/**************************************/

void
DFT_64_gmp (mpz_t* A, mpz_t omega, mpz_t prime, mpz_t * precompute_pow_omega)
{
#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_gmp (A[0], A[32], prime);
  DFT2_gmp (A[16], A[48], prime);
  DFT2_gmp (A[8], A[40], prime);
  DFT2_gmp (A[24], A[56], prime);
  DFT2_gmp (A[4], A[36], prime);
  DFT2_gmp (A[20], A[52], prime);
  DFT2_gmp (A[12], A[44], prime);
  DFT2_gmp (A[28], A[60], prime);
  DFT2_gmp (A[2], A[34], prime);
  DFT2_gmp (A[18], A[50], prime);
  DFT2_gmp (A[10], A[42], prime);
  DFT2_gmp (A[26], A[58], prime);
  DFT2_gmp (A[6], A[38], prime);
  DFT2_gmp (A[22], A[54], prime);
  DFT2_gmp (A[14], A[46], prime);
  DFT2_gmp (A[30], A[62], prime);

  DFT2_gmp (A[1], A[33], prime);
  DFT2_gmp (A[17], A[49], prime);
  DFT2_gmp (A[9], A[41], prime);
  DFT2_gmp (A[25], A[57], prime);
  DFT2_gmp (A[5], A[37], prime);
  DFT2_gmp (A[21], A[53], prime);
  DFT2_gmp (A[13], A[45], prime);
  DFT2_gmp (A[29], A[61], prime);
  DFT2_gmp (A[3], A[35], prime);
  DFT2_gmp (A[19], A[51], prime);
  DFT2_gmp (A[11], A[43], prime);
  DFT2_gmp (A[27], A[59], prime);
  DFT2_gmp (A[7], A[39], prime);
  DFT2_gmp (A[23], A[55], prime);
  DFT2_gmp (A[15], A[47], prime);
  DFT2_gmp (A[31], A[63], prime);

  //T_{2}^{4}
  mult_gmp (A[48], precompute_pow_omega[16], prime);
  mult_gmp (A[56], precompute_pow_omega[16], prime);
  mult_gmp (A[52], precompute_pow_omega[16], prime);
  mult_gmp (A[60], precompute_pow_omega[16], prime);
  mult_gmp (A[50], precompute_pow_omega[16], prime);
  mult_gmp (A[58], precompute_pow_omega[16], prime);
  mult_gmp (A[54], precompute_pow_omega[16], prime);
  mult_gmp (A[62], precompute_pow_omega[16], prime);
  mult_gmp (A[49], precompute_pow_omega[16], prime);
  mult_gmp (A[57], precompute_pow_omega[16], prime);
  mult_gmp (A[53], precompute_pow_omega[16], prime);
  mult_gmp (A[61], precompute_pow_omega[16], prime);
  mult_gmp (A[51], precompute_pow_omega[16], prime);
  mult_gmp (A[59], precompute_pow_omega[16], prime);
  mult_gmp (A[55], precompute_pow_omega[16], prime);
  mult_gmp (A[63], precompute_pow_omega[16], prime);

  //I_{32} \tx DFT2
  DFT2_gmp (A[0], A[16], prime);
  DFT2_gmp (A[32], A[48], prime);
  DFT2_gmp (A[8], A[24], prime);
  DFT2_gmp (A[40], A[56], prime);
  DFT2_gmp (A[4], A[20], prime);
  DFT2_gmp (A[36], A[52], prime);
  DFT2_gmp (A[12], A[28], prime);
  DFT2_gmp (A[44], A[60], prime);
  DFT2_gmp (A[2], A[18], prime);
  DFT2_gmp (A[34], A[50], prime);
  DFT2_gmp (A[10], A[26], prime);
  DFT2_gmp (A[42], A[58], prime);
  DFT2_gmp (A[6], A[22], prime);
  DFT2_gmp (A[38], A[54], prime);
  DFT2_gmp (A[14], A[30], prime);
  DFT2_gmp (A[46], A[62], prime);
  DFT2_gmp (A[1], A[17], prime);
  DFT2_gmp (A[33], A[49], prime);
  DFT2_gmp (A[9], A[25], prime);
  DFT2_gmp (A[41], A[57], prime);
  DFT2_gmp (A[5], A[21], prime);
  DFT2_gmp (A[37], A[53], prime);
  DFT2_gmp (A[13], A[29], prime);
  DFT2_gmp (A[45], A[61], prime);
  DFT2_gmp (A[3], A[19], prime);
  DFT2_gmp (A[35], A[51], prime);
  DFT2_gmp (A[11], A[27], prime);
  DFT2_gmp (A[43], A[59], prime);
  DFT2_gmp (A[7], A[23], prime);
  DFT2_gmp (A[39], A[55], prime);
  DFT2_gmp (A[15], A[31], prime);
  DFT2_gmp (A[47], A[63], prime);

  mult_gmp (A[40], precompute_pow_omega[8], prime);
  mult_gmp (A[24], precompute_pow_omega[16], prime);
  mult_gmp (A[56], precompute_pow_omega[24], prime);

  mult_gmp (A[44], precompute_pow_omega[8], prime);
  mult_gmp (A[28], precompute_pow_omega[16], prime);
  mult_gmp (A[60], precompute_pow_omega[24], prime);

  mult_gmp (A[42], precompute_pow_omega[8], prime);
  mult_gmp (A[26], precompute_pow_omega[16], prime);
  mult_gmp (A[58], precompute_pow_omega[24], prime);

  mult_gmp (A[46], precompute_pow_omega[8], prime);
  mult_gmp (A[30], precompute_pow_omega[16], prime);
  mult_gmp (A[62], precompute_pow_omega[24], prime);

  mult_gmp (A[41], precompute_pow_omega[8], prime);
  mult_gmp (A[25], precompute_pow_omega[16], prime);
  mult_gmp (A[57], precompute_pow_omega[24], prime);

  mult_gmp (A[45], precompute_pow_omega[8], prime);
  mult_gmp (A[29], precompute_pow_omega[16], prime);
  mult_gmp (A[61], precompute_pow_omega[24], prime);

  mult_gmp (A[43], precompute_pow_omega[8], prime);
  mult_gmp (A[27], precompute_pow_omega[16], prime);
  mult_gmp (A[59], precompute_pow_omega[24], prime);

  mult_gmp (A[47], precompute_pow_omega[8], prime);
  mult_gmp (A[31], precompute_pow_omega[16], prime);
  mult_gmp (A[63], precompute_pow_omega[24], prime);

  DFT2_gmp (A[0], A[8], prime);
  DFT2_gmp (A[32], A[40], prime);
  DFT2_gmp (A[16], A[24], prime);
  DFT2_gmp (A[48], A[56], prime);
  DFT2_gmp (A[4], A[12], prime);
  DFT2_gmp (A[36], A[44], prime);
  DFT2_gmp (A[20], A[28], prime);
  DFT2_gmp (A[52], A[60], prime);

  DFT2_gmp (A[2], A[10], prime);
  DFT2_gmp (A[34], A[42], prime);
  DFT2_gmp (A[18], A[26], prime);
  DFT2_gmp (A[50], A[58], prime);
  DFT2_gmp (A[6], A[14], prime);
  DFT2_gmp (A[38], A[46], prime);
  DFT2_gmp (A[22], A[30], prime);
  DFT2_gmp (A[54], A[62], prime);

  DFT2_gmp (A[1], A[9], prime);
  DFT2_gmp (A[33], A[41], prime);
  DFT2_gmp (A[17], A[25], prime);
  DFT2_gmp (A[49], A[57], prime);
  DFT2_gmp (A[5], A[13], prime);
  DFT2_gmp (A[37], A[45], prime);
  DFT2_gmp (A[21], A[29], prime);
  DFT2_gmp (A[53], A[61], prime);

  DFT2_gmp (A[3], A[11], prime);
  DFT2_gmp (A[35], A[43], prime);
  DFT2_gmp (A[19], A[27], prime);
  DFT2_gmp (A[51], A[59], prime);
  DFT2_gmp (A[7], A[15], prime);
  DFT2_gmp (A[39], A[47], prime);
  DFT2_gmp (A[23], A[31], prime);
  DFT2_gmp (A[55], A[63], prime);

  mult_gmp (A[36], precompute_pow_omega[4], prime);
  mult_gmp (A[20], precompute_pow_omega[8], prime);
  mult_gmp (A[52], precompute_pow_omega[12], prime);
  mult_gmp (A[12], precompute_pow_omega[16], prime);
  mult_gmp (A[44], precompute_pow_omega[20], prime);
  mult_gmp (A[28], precompute_pow_omega[24], prime);
  mult_gmp (A[60], precompute_pow_omega[28], prime);

  mult_gmp (A[38], precompute_pow_omega[4], prime);
  mult_gmp (A[22], precompute_pow_omega[8], prime);
  mult_gmp (A[54], precompute_pow_omega[12], prime);
  mult_gmp (A[14], precompute_pow_omega[16], prime);
  mult_gmp (A[46], precompute_pow_omega[20], prime);
  mult_gmp (A[30], precompute_pow_omega[24], prime);
  mult_gmp (A[62], precompute_pow_omega[28], prime);

  mult_gmp (A[37], precompute_pow_omega[4], prime);
  mult_gmp (A[21], precompute_pow_omega[8], prime);
  mult_gmp (A[53], precompute_pow_omega[12], prime);
  mult_gmp (A[13], precompute_pow_omega[16], prime);
  mult_gmp (A[45], precompute_pow_omega[20], prime);
  mult_gmp (A[29], precompute_pow_omega[24], prime);
  mult_gmp (A[61], precompute_pow_omega[28], prime);
  mult_gmp (A[39], precompute_pow_omega[4], prime);
  mult_gmp (A[23], precompute_pow_omega[8], prime);
  mult_gmp (A[55], precompute_pow_omega[12], prime);
  mult_gmp (A[15], precompute_pow_omega[16], prime);
  mult_gmp (A[47], precompute_pow_omega[20], prime);
  mult_gmp (A[31], precompute_pow_omega[24], prime);
  mult_gmp (A[63], precompute_pow_omega[28], prime);

  DFT2_gmp (A[0], A[4], prime);
  DFT2_gmp (A[32], A[36], prime);
  DFT2_gmp (A[16], A[20], prime);
  DFT2_gmp (A[48], A[52], prime);
  DFT2_gmp (A[8], A[12], prime);
  DFT2_gmp (A[40], A[44], prime);
  DFT2_gmp (A[24], A[28], prime);
  DFT2_gmp (A[56], A[60], prime);

  DFT2_gmp (A[2], A[6], prime);
  DFT2_gmp (A[34], A[38], prime);
  DFT2_gmp (A[18], A[22], prime);
  DFT2_gmp (A[50], A[54], prime);
  DFT2_gmp (A[10], A[14], prime);
  DFT2_gmp (A[42], A[46], prime);
  DFT2_gmp (A[26], A[30], prime);
  DFT2_gmp (A[58], A[62], prime);

  DFT2_gmp (A[1], A[5], prime);
  DFT2_gmp (A[33], A[37], prime);
  DFT2_gmp (A[17], A[21], prime);
  DFT2_gmp (A[49], A[53], prime);
  DFT2_gmp (A[9], A[13], prime);
  DFT2_gmp (A[41], A[45], prime);
  DFT2_gmp (A[25], A[29], prime);
  DFT2_gmp (A[57], A[61], prime);

  DFT2_gmp (A[3], A[7], prime);
  DFT2_gmp (A[35], A[39], prime);
  DFT2_gmp (A[19], A[23], prime);
  DFT2_gmp (A[51], A[55], prime);
  DFT2_gmp (A[11], A[15], prime);
  DFT2_gmp (A[43], A[47], prime);
  DFT2_gmp (A[27], A[31], prime);
  DFT2_gmp (A[59], A[63], prime);

  mult_gmp (A[34], precompute_pow_omega[2], prime);
  mult_gmp (A[18], precompute_pow_omega[4], prime);
  mult_gmp (A[50], precompute_pow_omega[6], prime);
  mult_gmp (A[10], precompute_pow_omega[8], prime);
  mult_gmp (A[42], precompute_pow_omega[10], prime);
  mult_gmp (A[26], precompute_pow_omega[12], prime);
  mult_gmp (A[58], precompute_pow_omega[14], prime);
  mult_gmp (A[6], precompute_pow_omega[16], prime);
  mult_gmp (A[38], precompute_pow_omega[18], prime);
  mult_gmp (A[22], precompute_pow_omega[20], prime);
  mult_gmp (A[54], precompute_pow_omega[22], prime);
  mult_gmp (A[14], precompute_pow_omega[24], prime);
  mult_gmp (A[46], precompute_pow_omega[26], prime);
  mult_gmp (A[30], precompute_pow_omega[28], prime);
  mult_gmp (A[62], precompute_pow_omega[30], prime);

  mult_gmp (A[35], precompute_pow_omega[2], prime);
  mult_gmp (A[19], precompute_pow_omega[4], prime);
  mult_gmp (A[51], precompute_pow_omega[6], prime);
  mult_gmp (A[11], precompute_pow_omega[8], prime);
  mult_gmp (A[43], precompute_pow_omega[10], prime);
  mult_gmp (A[27], precompute_pow_omega[12], prime);
  mult_gmp (A[59], precompute_pow_omega[14], prime);
  mult_gmp (A[7], precompute_pow_omega[16], prime);
  mult_gmp (A[39], precompute_pow_omega[18], prime);
  mult_gmp (A[23], precompute_pow_omega[20], prime);
  mult_gmp (A[55], precompute_pow_omega[22], prime);
  mult_gmp (A[15], precompute_pow_omega[24], prime);
  mult_gmp (A[47], precompute_pow_omega[26], prime);
  mult_gmp (A[31], precompute_pow_omega[28], prime);
  mult_gmp (A[63], precompute_pow_omega[30], prime);

  DFT2_gmp (A[0], A[2], prime);
  DFT2_gmp (A[32], A[34], prime);
  DFT2_gmp (A[16], A[18], prime);
  DFT2_gmp (A[48], A[50], prime);
  DFT2_gmp (A[8], A[10], prime);
  DFT2_gmp (A[40], A[42], prime);
  DFT2_gmp (A[24], A[26], prime);
  DFT2_gmp (A[56], A[58], prime);

  DFT2_gmp (A[4], A[6], prime);
  DFT2_gmp (A[36], A[38], prime);
  DFT2_gmp (A[20], A[22], prime);
  DFT2_gmp (A[52], A[54], prime);
  DFT2_gmp (A[12], A[14], prime);
  DFT2_gmp (A[44], A[46], prime);
  DFT2_gmp (A[28], A[30], prime);
  DFT2_gmp (A[60], A[62], prime);

  DFT2_gmp (A[1], A[3], prime);
  DFT2_gmp (A[33], A[35], prime);
  DFT2_gmp (A[17], A[19], prime);
  DFT2_gmp (A[49], A[51], prime);
  DFT2_gmp (A[9], A[11], prime);
  DFT2_gmp (A[41], A[43], prime);
  DFT2_gmp (A[25], A[27], prime);
  DFT2_gmp (A[57], A[59], prime);

  DFT2_gmp (A[5], A[7], prime);
  DFT2_gmp (A[37], A[39], prime);
  DFT2_gmp (A[21], A[23], prime);
  DFT2_gmp (A[53], A[55], prime);
  DFT2_gmp (A[13], A[15], prime);
  DFT2_gmp (A[45], A[47], prime);
  DFT2_gmp (A[29], A[31], prime);
  DFT2_gmp (A[61], A[63], prime);

  mult_gmp (A[33], precompute_pow_omega[1], prime);
  mult_gmp (A[17], precompute_pow_omega[2], prime);
  mult_gmp (A[49], precompute_pow_omega[3], prime);
  mult_gmp (A[9], precompute_pow_omega[4], prime);
  mult_gmp (A[41], precompute_pow_omega[5], prime);
  mult_gmp (A[25], precompute_pow_omega[6], prime);
  mult_gmp (A[57], precompute_pow_omega[7], prime);
  mult_gmp (A[5], precompute_pow_omega[8], prime);
  mult_gmp (A[37], precompute_pow_omega[9], prime);
  mult_gmp (A[21], precompute_pow_omega[10], prime);
  mult_gmp (A[53], precompute_pow_omega[11], prime);
  mult_gmp (A[13], precompute_pow_omega[12], prime);
  mult_gmp (A[45], precompute_pow_omega[13], prime);
  mult_gmp (A[29], precompute_pow_omega[14], prime);
  mult_gmp (A[61], precompute_pow_omega[15], prime);
  mult_gmp (A[3], precompute_pow_omega[16], prime);
  mult_gmp (A[35], precompute_pow_omega[17], prime);
  mult_gmp (A[19], precompute_pow_omega[18], prime);
  mult_gmp (A[51], precompute_pow_omega[19], prime);
  mult_gmp (A[11], precompute_pow_omega[20], prime);
  mult_gmp (A[43], precompute_pow_omega[21], prime);
  mult_gmp (A[27], precompute_pow_omega[22], prime);
  mult_gmp (A[59], precompute_pow_omega[23], prime);
  mult_gmp (A[7], precompute_pow_omega[24], prime);
  mult_gmp (A[39], precompute_pow_omega[25], prime);
  mult_gmp (A[23], precompute_pow_omega[26], prime);
  mult_gmp (A[55], precompute_pow_omega[27], prime);
  mult_gmp (A[15], precompute_pow_omega[28], prime);
  mult_gmp (A[47], precompute_pow_omega[29], prime);
  mult_gmp (A[31], precompute_pow_omega[30], prime);
  mult_gmp (A[63], precompute_pow_omega[31], prime);

  DFT2_gmp (A[0], A[1], prime);
  DFT2_gmp (A[32], A[33], prime);
  DFT2_gmp (A[16], A[17], prime);
  DFT2_gmp (A[48], A[49], prime);
  DFT2_gmp (A[8], A[9], prime);
  DFT2_gmp (A[40], A[41], prime);
  DFT2_gmp (A[24], A[25], prime);
  DFT2_gmp (A[56], A[57], prime);

  DFT2_gmp (A[4], A[5], prime);
  DFT2_gmp (A[36], A[37], prime);
  DFT2_gmp (A[20], A[21], prime);
  DFT2_gmp (A[52], A[53], prime);
  DFT2_gmp (A[12], A[13], prime);
  DFT2_gmp (A[44], A[45], prime);
  DFT2_gmp (A[28], A[29], prime);
  DFT2_gmp (A[60], A[61], prime);

  DFT2_gmp (A[2], A[3], prime);
  DFT2_gmp (A[34], A[35], prime);
  DFT2_gmp (A[18], A[19], prime);
  DFT2_gmp (A[50], A[51], prime);
  DFT2_gmp (A[10], A[11], prime);
  DFT2_gmp (A[42], A[43], prime);
  DFT2_gmp (A[26], A[27], prime);
  DFT2_gmp (A[58], A[59], prime);

  DFT2_gmp (A[6], A[7], prime);
  DFT2_gmp (A[38], A[39], prime);
  DFT2_gmp (A[22], A[23], prime);
  DFT2_gmp (A[54], A[55], prime);
  DFT2_gmp (A[14], A[15], prime);
  DFT2_gmp (A[46], A[47], prime);
  DFT2_gmp (A[30], A[31], prime);
  DFT2_gmp (A[62], A[63], prime);

  // Final permutation
  swap_gmp (A[1], A[32]);
  swap_gmp (A[2], A[16]);
  swap_gmp (A[3], A[48]);
  swap_gmp (A[4], A[8]);
  swap_gmp (A[5], A[40]);
  swap_gmp (A[6], A[24]);
  swap_gmp (A[7], A[56]);
  swap_gmp (A[9], A[36]);
  swap_gmp (A[10], A[20]);
  swap_gmp (A[11], A[52]);
  swap_gmp (A[13], A[44]);
  swap_gmp (A[14], A[28]);
  swap_gmp (A[15], A[60]);
  swap_gmp (A[17], A[34]);
  swap_gmp (A[19], A[50]);
  swap_gmp (A[21], A[42]);
  swap_gmp (A[22], A[26]);
  swap_gmp (A[23], A[58]);
  swap_gmp (A[25], A[38]);
  swap_gmp (A[27], A[54]);
  swap_gmp (A[29], A[46]);
  swap_gmp (A[31], A[62]);
  swap_gmp (A[35], A[49]);
  swap_gmp (A[37], A[41]);
  swap_gmp (A[39], A[57]);
  swap_gmp (A[43], A[53]);
  swap_gmp (A[47], A[61]);
  swap_gmp (A[55], A[59]);
}

/**************************************/

void
DFT_64_aux_gmp (mpz_t* A, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega)
{
#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_aux_gmp (A[0], A[32], aux, prime);
  DFT2_aux_gmp (A[16], A[48], aux, prime);
  DFT2_aux_gmp (A[8], A[40], aux, prime);
  DFT2_aux_gmp (A[24], A[56], aux, prime);
  DFT2_aux_gmp (A[4], A[36], aux, prime);
  DFT2_aux_gmp (A[20], A[52], aux, prime);
  DFT2_aux_gmp (A[12], A[44], aux, prime);
  DFT2_aux_gmp (A[28], A[60], aux, prime);
  DFT2_aux_gmp (A[2], A[34], aux, prime);
  DFT2_aux_gmp (A[18], A[50], aux, prime);
  DFT2_aux_gmp (A[10], A[42], aux, prime);
  DFT2_aux_gmp (A[26], A[58], aux, prime);
  DFT2_aux_gmp (A[6], A[38], aux, prime);
  DFT2_aux_gmp (A[22], A[54], aux, prime);
  DFT2_aux_gmp (A[14], A[46], aux, prime);
  DFT2_aux_gmp (A[30], A[62], aux, prime);

  DFT2_aux_gmp (A[1], A[33], aux, prime);
  DFT2_aux_gmp (A[17], A[49], aux, prime);
  DFT2_aux_gmp (A[9], A[41], aux, prime);
  DFT2_aux_gmp (A[25], A[57], aux, prime);
  DFT2_aux_gmp (A[5], A[37], aux, prime);
  DFT2_aux_gmp (A[21], A[53], aux, prime);
  DFT2_aux_gmp (A[13], A[45], aux, prime);
  DFT2_aux_gmp (A[29], A[61], aux, prime);
  DFT2_aux_gmp (A[3], A[35], aux, prime);
  DFT2_aux_gmp (A[19], A[51], aux, prime);
  DFT2_aux_gmp (A[11], A[43], aux, prime);
  DFT2_aux_gmp (A[27], A[59], aux, prime);
  DFT2_aux_gmp (A[7], A[39], aux, prime);
  DFT2_aux_gmp (A[23], A[55], aux, prime);
  DFT2_aux_gmp (A[15], A[47], aux, prime);
  DFT2_aux_gmp (A[31], A[63], aux, prime);

  //T_{2}^{4}
  mult_aux_gmp (A[48], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[56], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[52], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[60], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[50], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[58], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[54], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[62], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[49], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[57], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[53], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[61], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[51], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[59], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[55], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[63], precompute_pow_omega[16], aux, prime);

  //I_{32} \tx DFT2
  DFT2_aux_gmp (A[0], A[16], aux, prime);
  DFT2_aux_gmp (A[32], A[48], aux, prime);
  DFT2_aux_gmp (A[8], A[24], aux, prime);
  DFT2_aux_gmp (A[40], A[56], aux, prime);
  DFT2_aux_gmp (A[4], A[20], aux, prime);
  DFT2_aux_gmp (A[36], A[52], aux, prime);
  DFT2_aux_gmp (A[12], A[28], aux, prime);
  DFT2_aux_gmp (A[44], A[60], aux, prime);
  DFT2_aux_gmp (A[2], A[18], aux, prime);
  DFT2_aux_gmp (A[34], A[50], aux, prime);
  DFT2_aux_gmp (A[10], A[26], aux, prime);
  DFT2_aux_gmp (A[42], A[58], aux, prime);
  DFT2_aux_gmp (A[6], A[22], aux, prime);
  DFT2_aux_gmp (A[38], A[54], aux, prime);
  DFT2_aux_gmp (A[14], A[30], aux, prime);
  DFT2_aux_gmp (A[46], A[62], aux, prime);
  DFT2_aux_gmp (A[1], A[17], aux, prime);
  DFT2_aux_gmp (A[33], A[49], aux, prime);
  DFT2_aux_gmp (A[9], A[25], aux, prime);
  DFT2_aux_gmp (A[41], A[57], aux, prime);
  DFT2_aux_gmp (A[5], A[21], aux, prime);
  DFT2_aux_gmp (A[37], A[53], aux, prime);
  DFT2_aux_gmp (A[13], A[29], aux, prime);
  DFT2_aux_gmp (A[45], A[61], aux, prime);
  DFT2_aux_gmp (A[3], A[19], aux, prime);
  DFT2_aux_gmp (A[35], A[51], aux, prime);
  DFT2_aux_gmp (A[11], A[27], aux, prime);
  DFT2_aux_gmp (A[43], A[59], aux, prime);
  DFT2_aux_gmp (A[7], A[23], aux, prime);
  DFT2_aux_gmp (A[39], A[55], aux, prime);
  DFT2_aux_gmp (A[15], A[31], aux, prime);
  DFT2_aux_gmp (A[47], A[63], aux, prime);

  mult_aux_gmp (A[40], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[24], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[56], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[44], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[28], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[60], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[42], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[26], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[58], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[46], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[30], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[62], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[41], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[25], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[57], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[45], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[29], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[61], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[43], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[27], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[59], precompute_pow_omega[24], aux, prime);

  mult_aux_gmp (A[47], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[31], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[63], precompute_pow_omega[24], aux, prime);

  DFT2_aux_gmp (A[0], A[8], aux, prime);
  DFT2_aux_gmp (A[32], A[40], aux, prime);
  DFT2_aux_gmp (A[16], A[24], aux, prime);
  DFT2_aux_gmp (A[48], A[56], aux, prime);
  DFT2_aux_gmp (A[4], A[12], aux, prime);
  DFT2_aux_gmp (A[36], A[44], aux, prime);
  DFT2_aux_gmp (A[20], A[28], aux, prime);
  DFT2_aux_gmp (A[52], A[60], aux, prime);

  DFT2_aux_gmp (A[2], A[10], aux, prime);
  DFT2_aux_gmp (A[34], A[42], aux, prime);
  DFT2_aux_gmp (A[18], A[26], aux, prime);
  DFT2_aux_gmp (A[50], A[58], aux, prime);
  DFT2_aux_gmp (A[6], A[14], aux, prime);
  DFT2_aux_gmp (A[38], A[46], aux, prime);
  DFT2_aux_gmp (A[22], A[30], aux, prime);
  DFT2_aux_gmp (A[54], A[62], aux, prime);

  DFT2_aux_gmp (A[1], A[9], aux, prime);
  DFT2_aux_gmp (A[33], A[41], aux, prime);
  DFT2_aux_gmp (A[17], A[25], aux, prime);
  DFT2_aux_gmp (A[49], A[57], aux, prime);
  DFT2_aux_gmp (A[5], A[13], aux, prime);
  DFT2_aux_gmp (A[37], A[45], aux, prime);
  DFT2_aux_gmp (A[21], A[29], aux, prime);
  DFT2_aux_gmp (A[53], A[61], aux, prime);

  DFT2_aux_gmp (A[3], A[11], aux, prime);
  DFT2_aux_gmp (A[35], A[43], aux, prime);
  DFT2_aux_gmp (A[19], A[27], aux, prime);
  DFT2_aux_gmp (A[51], A[59], aux, prime);
  DFT2_aux_gmp (A[7], A[15], aux, prime);
  DFT2_aux_gmp (A[39], A[47], aux, prime);
  DFT2_aux_gmp (A[23], A[31], aux, prime);
  DFT2_aux_gmp (A[55], A[63], aux, prime);

  mult_aux_gmp (A[36], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[20], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[52], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[12], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[44], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[28], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[60], precompute_pow_omega[28], aux, prime);

  mult_aux_gmp (A[38], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[22], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[54], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[14], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[46], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[30], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[62], precompute_pow_omega[28], aux, prime);

  mult_aux_gmp (A[37], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[21], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[53], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[13], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[45], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[29], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[61], precompute_pow_omega[28], aux, prime);
  mult_aux_gmp (A[39], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[23], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[55], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[15], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[47], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[31], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[63], precompute_pow_omega[28], aux, prime);

  DFT2_aux_gmp (A[0], A[4], aux, prime);
  DFT2_aux_gmp (A[32], A[36], aux, prime);
  DFT2_aux_gmp (A[16], A[20], aux, prime);
  DFT2_aux_gmp (A[48], A[52], aux, prime);
  DFT2_aux_gmp (A[8], A[12], aux, prime);
  DFT2_aux_gmp (A[40], A[44], aux, prime);
  DFT2_aux_gmp (A[24], A[28], aux, prime);
  DFT2_aux_gmp (A[56], A[60], aux, prime);

  DFT2_aux_gmp (A[2], A[6], aux, prime);
  DFT2_aux_gmp (A[34], A[38], aux, prime);
  DFT2_aux_gmp (A[18], A[22], aux, prime);
  DFT2_aux_gmp (A[50], A[54], aux, prime);
  DFT2_aux_gmp (A[10], A[14], aux, prime);
  DFT2_aux_gmp (A[42], A[46], aux, prime);
  DFT2_aux_gmp (A[26], A[30], aux, prime);
  DFT2_aux_gmp (A[58], A[62], aux, prime);

  DFT2_aux_gmp (A[1], A[5], aux, prime);
  DFT2_aux_gmp (A[33], A[37], aux, prime);
  DFT2_aux_gmp (A[17], A[21], aux, prime);
  DFT2_aux_gmp (A[49], A[53], aux, prime);
  DFT2_aux_gmp (A[9], A[13], aux, prime);
  DFT2_aux_gmp (A[41], A[45], aux, prime);
  DFT2_aux_gmp (A[25], A[29], aux, prime);
  DFT2_aux_gmp (A[57], A[61], aux, prime);

  DFT2_aux_gmp (A[3], A[7], aux, prime);
  DFT2_aux_gmp (A[35], A[39], aux, prime);
  DFT2_aux_gmp (A[19], A[23], aux, prime);
  DFT2_aux_gmp (A[51], A[55], aux, prime);
  DFT2_aux_gmp (A[11], A[15], aux, prime);
  DFT2_aux_gmp (A[43], A[47], aux, prime);
  DFT2_aux_gmp (A[27], A[31], aux, prime);
  DFT2_aux_gmp (A[59], A[63], aux, prime);

  mult_aux_gmp (A[34], precompute_pow_omega[2], aux, prime);
  mult_aux_gmp (A[18], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[50], precompute_pow_omega[6], aux, prime);
  mult_aux_gmp (A[10], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[42], precompute_pow_omega[10], aux, prime);
  mult_aux_gmp (A[26], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[58], precompute_pow_omega[14], aux, prime);
  mult_aux_gmp (A[6], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[38], precompute_pow_omega[18], aux, prime);
  mult_aux_gmp (A[22], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[54], precompute_pow_omega[22], aux, prime);
  mult_aux_gmp (A[14], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[46], precompute_pow_omega[26], aux, prime);
  mult_aux_gmp (A[30], precompute_pow_omega[28], aux, prime);
  mult_aux_gmp (A[62], precompute_pow_omega[30], aux, prime);

  mult_aux_gmp (A[35], precompute_pow_omega[2], aux, prime);
  mult_aux_gmp (A[19], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[51], precompute_pow_omega[6], aux, prime);
  mult_aux_gmp (A[11], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[43], precompute_pow_omega[10], aux, prime);
  mult_aux_gmp (A[27], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[59], precompute_pow_omega[14], aux, prime);
  mult_aux_gmp (A[7], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[39], precompute_pow_omega[18], aux, prime);
  mult_aux_gmp (A[23], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[55], precompute_pow_omega[22], aux, prime);
  mult_aux_gmp (A[15], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[47], precompute_pow_omega[26], aux, prime);
  mult_aux_gmp (A[31], precompute_pow_omega[28], aux, prime);
  mult_aux_gmp (A[63], precompute_pow_omega[30], aux, prime);

  DFT2_aux_gmp (A[0], A[2], aux, prime);
  DFT2_aux_gmp (A[32], A[34], aux, prime);
  DFT2_aux_gmp (A[16], A[18], aux, prime);
  DFT2_aux_gmp (A[48], A[50], aux, prime);
  DFT2_aux_gmp (A[8], A[10], aux, prime);
  DFT2_aux_gmp (A[40], A[42], aux, prime);
  DFT2_aux_gmp (A[24], A[26], aux, prime);
  DFT2_aux_gmp (A[56], A[58], aux, prime);

  DFT2_aux_gmp (A[4], A[6], aux, prime);
  DFT2_aux_gmp (A[36], A[38], aux, prime);
  DFT2_aux_gmp (A[20], A[22], aux, prime);
  DFT2_aux_gmp (A[52], A[54], aux, prime);
  DFT2_aux_gmp (A[12], A[14], aux, prime);
  DFT2_aux_gmp (A[44], A[46], aux, prime);
  DFT2_aux_gmp (A[28], A[30], aux, prime);
  DFT2_aux_gmp (A[60], A[62], aux, prime);

  DFT2_aux_gmp (A[1], A[3], aux, prime);
  DFT2_aux_gmp (A[33], A[35], aux, prime);
  DFT2_aux_gmp (A[17], A[19], aux, prime);
  DFT2_aux_gmp (A[49], A[51], aux, prime);
  DFT2_aux_gmp (A[9], A[11], aux, prime);
  DFT2_aux_gmp (A[41], A[43], aux, prime);
  DFT2_aux_gmp (A[25], A[27], aux, prime);
  DFT2_aux_gmp (A[57], A[59], aux, prime);

  DFT2_aux_gmp (A[5], A[7], aux, prime);
  DFT2_aux_gmp (A[37], A[39], aux, prime);
  DFT2_aux_gmp (A[21], A[23], aux, prime);
  DFT2_aux_gmp (A[53], A[55], aux, prime);
  DFT2_aux_gmp (A[13], A[15], aux, prime);
  DFT2_aux_gmp (A[45], A[47], aux, prime);
  DFT2_aux_gmp (A[29], A[31], aux, prime);
  DFT2_aux_gmp (A[61], A[63], aux, prime);

  mult_aux_gmp (A[33], precompute_pow_omega[1], aux, prime);
  mult_aux_gmp (A[17], precompute_pow_omega[2], aux, prime);
  mult_aux_gmp (A[49], precompute_pow_omega[3], aux, prime);
  mult_aux_gmp (A[9], precompute_pow_omega[4], aux, prime);
  mult_aux_gmp (A[41], precompute_pow_omega[5], aux, prime);
  mult_aux_gmp (A[25], precompute_pow_omega[6], aux, prime);
  mult_aux_gmp (A[57], precompute_pow_omega[7], aux, prime);
  mult_aux_gmp (A[5], precompute_pow_omega[8], aux, prime);
  mult_aux_gmp (A[37], precompute_pow_omega[9], aux, prime);
  mult_aux_gmp (A[21], precompute_pow_omega[10], aux, prime);
  mult_aux_gmp (A[53], precompute_pow_omega[11], aux, prime);
  mult_aux_gmp (A[13], precompute_pow_omega[12], aux, prime);
  mult_aux_gmp (A[45], precompute_pow_omega[13], aux, prime);
  mult_aux_gmp (A[29], precompute_pow_omega[14], aux, prime);
  mult_aux_gmp (A[61], precompute_pow_omega[15], aux, prime);
  mult_aux_gmp (A[3], precompute_pow_omega[16], aux, prime);
  mult_aux_gmp (A[35], precompute_pow_omega[17], aux, prime);
  mult_aux_gmp (A[19], precompute_pow_omega[18], aux, prime);
  mult_aux_gmp (A[51], precompute_pow_omega[19], aux, prime);
  mult_aux_gmp (A[11], precompute_pow_omega[20], aux, prime);
  mult_aux_gmp (A[43], precompute_pow_omega[21], aux, prime);
  mult_aux_gmp (A[27], precompute_pow_omega[22], aux, prime);
  mult_aux_gmp (A[59], precompute_pow_omega[23], aux, prime);
  mult_aux_gmp (A[7], precompute_pow_omega[24], aux, prime);
  mult_aux_gmp (A[39], precompute_pow_omega[25], aux, prime);
  mult_aux_gmp (A[23], precompute_pow_omega[26], aux, prime);
  mult_aux_gmp (A[55], precompute_pow_omega[27], aux, prime);
  mult_aux_gmp (A[15], precompute_pow_omega[28], aux, prime);
  mult_aux_gmp (A[47], precompute_pow_omega[29], aux, prime);
  mult_aux_gmp (A[31], precompute_pow_omega[30], aux, prime);
  mult_aux_gmp (A[63], precompute_pow_omega[31], aux, prime);

  DFT2_aux_gmp (A[0], A[1], aux, prime);
  DFT2_aux_gmp (A[32], A[33], aux, prime);
  DFT2_aux_gmp (A[16], A[17], aux, prime);
  DFT2_aux_gmp (A[48], A[49], aux, prime);
  DFT2_aux_gmp (A[8], A[9], aux, prime);
  DFT2_aux_gmp (A[40], A[41], aux, prime);
  DFT2_aux_gmp (A[24], A[25], aux, prime);
  DFT2_aux_gmp (A[56], A[57], aux, prime);

  DFT2_aux_gmp (A[4], A[5], aux, prime);
  DFT2_aux_gmp (A[36], A[37], aux, prime);
  DFT2_aux_gmp (A[20], A[21], aux, prime);
  DFT2_aux_gmp (A[52], A[53], aux, prime);
  DFT2_aux_gmp (A[12], A[13], aux, prime);
  DFT2_aux_gmp (A[44], A[45], aux, prime);
  DFT2_aux_gmp (A[28], A[29], aux, prime);
  DFT2_aux_gmp (A[60], A[61], aux, prime);

  DFT2_aux_gmp (A[2], A[3], aux, prime);
  DFT2_aux_gmp (A[34], A[35], aux, prime);
  DFT2_aux_gmp (A[18], A[19], aux, prime);
  DFT2_aux_gmp (A[50], A[51], aux, prime);
  DFT2_aux_gmp (A[10], A[11], aux, prime);
  DFT2_aux_gmp (A[42], A[43], aux, prime);
  DFT2_aux_gmp (A[26], A[27], aux, prime);
  DFT2_aux_gmp (A[58], A[59], aux, prime);

  DFT2_aux_gmp (A[6], A[7], aux, prime);
  DFT2_aux_gmp (A[38], A[39], aux, prime);
  DFT2_aux_gmp (A[22], A[23], aux, prime);
  DFT2_aux_gmp (A[54], A[55], aux, prime);
  DFT2_aux_gmp (A[14], A[15], aux, prime);
  DFT2_aux_gmp (A[46], A[47], aux, prime);
  DFT2_aux_gmp (A[30], A[31], aux, prime);
  DFT2_aux_gmp (A[62], A[63], aux, prime);

  // Final permutation
  swap_aux_gmp (A[1], A[32], aux);
  swap_aux_gmp (A[2], A[16], aux);
  swap_aux_gmp (A[3], A[48], aux);
  swap_aux_gmp (A[4], A[8], aux);
  swap_aux_gmp (A[5], A[40], aux);
  swap_aux_gmp (A[6], A[24], aux);
  swap_aux_gmp (A[7], A[56], aux);
  swap_aux_gmp (A[9], A[36], aux);
  swap_aux_gmp (A[10], A[20], aux);
  swap_aux_gmp (A[11], A[52], aux);
  swap_aux_gmp (A[13], A[44], aux);
  swap_aux_gmp (A[14], A[28], aux);
  swap_aux_gmp (A[15], A[60], aux);
  swap_aux_gmp (A[17], A[34], aux);
  swap_aux_gmp (A[19], A[50], aux);
  swap_aux_gmp (A[21], A[42], aux);
  swap_aux_gmp (A[22], A[26], aux);
  swap_aux_gmp (A[23], A[58], aux);
  swap_aux_gmp (A[25], A[38], aux);
  swap_aux_gmp (A[27], A[54], aux);
  swap_aux_gmp (A[29], A[46], aux);
  swap_aux_gmp (A[31], A[62], aux);
  swap_aux_gmp (A[35], A[49], aux);
  swap_aux_gmp (A[37], A[41], aux);
  swap_aux_gmp (A[39], A[57], aux);
  swap_aux_gmp (A[43], A[53], aux);
  swap_aux_gmp (A[47], A[61], aux);
  swap_aux_gmp (A[55], A[59], aux);
}

/**************************************/

void
DFT_128_aux_gmp (mpz_t* x, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega)
{
#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_aux_gmp(x[0], x[64], aux, prime);
  DFT2_aux_gmp(x[1], x[65], aux, prime);
  DFT2_aux_gmp(x[2], x[66], aux, prime);
  DFT2_aux_gmp(x[3], x[67], aux, prime);
  DFT2_aux_gmp(x[4], x[68], aux, prime);
  DFT2_aux_gmp(x[5], x[69], aux, prime);
  DFT2_aux_gmp(x[6], x[70], aux, prime);
  DFT2_aux_gmp(x[7], x[71], aux, prime);
  DFT2_aux_gmp(x[8], x[72], aux, prime);
  DFT2_aux_gmp(x[9], x[73], aux, prime);
  DFT2_aux_gmp(x[10], x[74], aux, prime);
  DFT2_aux_gmp(x[11], x[75], aux, prime);
  DFT2_aux_gmp(x[12], x[76], aux, prime);
  DFT2_aux_gmp(x[13], x[77], aux, prime);
  DFT2_aux_gmp(x[14], x[78], aux, prime);
  DFT2_aux_gmp(x[15], x[79], aux, prime);
  DFT2_aux_gmp(x[16], x[80], aux, prime);
  DFT2_aux_gmp(x[17], x[81], aux, prime);
  DFT2_aux_gmp(x[18], x[82], aux, prime);
  DFT2_aux_gmp(x[19], x[83], aux, prime);
  DFT2_aux_gmp(x[20], x[84], aux, prime);
  DFT2_aux_gmp(x[21], x[85], aux, prime);
  DFT2_aux_gmp(x[22], x[86], aux, prime);
  DFT2_aux_gmp(x[23], x[87], aux, prime);
  DFT2_aux_gmp(x[24], x[88], aux, prime);
  DFT2_aux_gmp(x[25], x[89], aux, prime);
  DFT2_aux_gmp(x[26], x[90], aux, prime);
  DFT2_aux_gmp(x[27], x[91], aux, prime);
  DFT2_aux_gmp(x[28], x[92], aux, prime);
  DFT2_aux_gmp(x[29], x[93], aux, prime);
  DFT2_aux_gmp(x[30], x[94], aux, prime);
  DFT2_aux_gmp(x[31], x[95], aux, prime);
  DFT2_aux_gmp(x[32], x[96], aux, prime);
  DFT2_aux_gmp(x[33], x[97], aux, prime);
  DFT2_aux_gmp(x[34], x[98], aux, prime);
  DFT2_aux_gmp(x[35], x[99], aux, prime);
  DFT2_aux_gmp(x[36], x[100], aux, prime);
  DFT2_aux_gmp(x[37], x[101], aux, prime);
  DFT2_aux_gmp(x[38], x[102], aux, prime);
  DFT2_aux_gmp(x[39], x[103], aux, prime);
  DFT2_aux_gmp(x[40], x[104], aux, prime);
  DFT2_aux_gmp(x[41], x[105], aux, prime);
  DFT2_aux_gmp(x[42], x[106], aux, prime);
  DFT2_aux_gmp(x[43], x[107], aux, prime);
  DFT2_aux_gmp(x[44], x[108], aux, prime);
  DFT2_aux_gmp(x[45], x[109], aux, prime);
  DFT2_aux_gmp(x[46], x[110], aux, prime);
  DFT2_aux_gmp(x[47], x[111], aux, prime);
  DFT2_aux_gmp(x[48], x[112], aux, prime);
  DFT2_aux_gmp(x[49], x[113], aux, prime);
  DFT2_aux_gmp(x[50], x[114], aux, prime);
  DFT2_aux_gmp(x[51], x[115], aux, prime);
  DFT2_aux_gmp(x[52], x[116], aux, prime);
  DFT2_aux_gmp(x[53], x[117], aux, prime);
  DFT2_aux_gmp(x[54], x[118], aux, prime);
  DFT2_aux_gmp(x[55], x[119], aux, prime);
  DFT2_aux_gmp(x[56], x[120], aux, prime);
  DFT2_aux_gmp(x[57], x[121], aux, prime);
  DFT2_aux_gmp(x[58], x[122], aux, prime);
  DFT2_aux_gmp(x[59], x[123], aux, prime);
  DFT2_aux_gmp(x[60], x[124], aux, prime);
  DFT2_aux_gmp(x[61], x[125], aux, prime);
  DFT2_aux_gmp(x[62], x[126], aux, prime);
  DFT2_aux_gmp(x[63], x[127], aux, prime);
  mult_aux_gmp(x[96], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[97], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[98], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[99], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[100], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[101], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[102], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[104], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[105], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[106], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[108], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[112], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[113], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[114], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[116], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[120], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[32], aux,	prime);
  DFT2_aux_gmp(x[0], x[32], aux, prime);
  DFT2_aux_gmp(x[1], x[33], aux, prime);
  DFT2_aux_gmp(x[2], x[34], aux, prime);
  DFT2_aux_gmp(x[3], x[35], aux, prime);
  DFT2_aux_gmp(x[4], x[36], aux, prime);
  DFT2_aux_gmp(x[5], x[37], aux, prime);
  DFT2_aux_gmp(x[6], x[38], aux, prime);
  DFT2_aux_gmp(x[7], x[39], aux, prime);
  DFT2_aux_gmp(x[8], x[40], aux, prime);
  DFT2_aux_gmp(x[9], x[41], aux, prime);
  DFT2_aux_gmp(x[10], x[42], aux, prime);
  DFT2_aux_gmp(x[11], x[43], aux, prime);
  DFT2_aux_gmp(x[12], x[44], aux, prime);
  DFT2_aux_gmp(x[13], x[45], aux, prime);
  DFT2_aux_gmp(x[14], x[46], aux, prime);
  DFT2_aux_gmp(x[15], x[47], aux, prime);
  DFT2_aux_gmp(x[16], x[48], aux, prime);
  DFT2_aux_gmp(x[17], x[49], aux, prime);
  DFT2_aux_gmp(x[18], x[50], aux, prime);
  DFT2_aux_gmp(x[19], x[51], aux, prime);
  DFT2_aux_gmp(x[20], x[52], aux, prime);
  DFT2_aux_gmp(x[21], x[53], aux, prime);
  DFT2_aux_gmp(x[22], x[54], aux, prime);
  DFT2_aux_gmp(x[23], x[55], aux, prime);
  DFT2_aux_gmp(x[24], x[56], aux, prime);
  DFT2_aux_gmp(x[25], x[57], aux, prime);
  DFT2_aux_gmp(x[26], x[58], aux, prime);
  DFT2_aux_gmp(x[27], x[59], aux, prime);
  DFT2_aux_gmp(x[28], x[60], aux, prime);
  DFT2_aux_gmp(x[29], x[61], aux, prime);
  DFT2_aux_gmp(x[30], x[62], aux, prime);
  DFT2_aux_gmp(x[31], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[96], aux, prime);
  DFT2_aux_gmp(x[65], x[97], aux, prime);
  DFT2_aux_gmp(x[66], x[98], aux, prime);
  DFT2_aux_gmp(x[67], x[99], aux, prime);
  DFT2_aux_gmp(x[68], x[100], aux, prime);
  DFT2_aux_gmp(x[69], x[101], aux, prime);
  DFT2_aux_gmp(x[70], x[102], aux, prime);
  DFT2_aux_gmp(x[71], x[103], aux, prime);
  DFT2_aux_gmp(x[72], x[104], aux, prime);
  DFT2_aux_gmp(x[73], x[105], aux, prime);
  DFT2_aux_gmp(x[74], x[106], aux, prime);
  DFT2_aux_gmp(x[75], x[107], aux, prime);
  DFT2_aux_gmp(x[76], x[108], aux, prime);
  DFT2_aux_gmp(x[77], x[109], aux, prime);
  DFT2_aux_gmp(x[78], x[110], aux, prime);
  DFT2_aux_gmp(x[79], x[111], aux, prime);
  DFT2_aux_gmp(x[80], x[112], aux, prime);
  DFT2_aux_gmp(x[81], x[113], aux, prime);
  DFT2_aux_gmp(x[82], x[114], aux, prime);
  DFT2_aux_gmp(x[83], x[115], aux, prime);
  DFT2_aux_gmp(x[84], x[116], aux, prime);
  DFT2_aux_gmp(x[85], x[117], aux, prime);
  DFT2_aux_gmp(x[86], x[118], aux, prime);
  DFT2_aux_gmp(x[87], x[119], aux, prime);
  DFT2_aux_gmp(x[88], x[120], aux, prime);
  DFT2_aux_gmp(x[89], x[121], aux, prime);
  DFT2_aux_gmp(x[90], x[122], aux, prime);
  DFT2_aux_gmp(x[91], x[123], aux, prime);
  DFT2_aux_gmp(x[92], x[124], aux, prime);
  DFT2_aux_gmp(x[93], x[125], aux, prime);
  DFT2_aux_gmp(x[94], x[126], aux, prime);
  DFT2_aux_gmp(x[95], x[127], aux, prime);
  mult_aux_gmp(x[48], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[49], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[50], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[51], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[52], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[53], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[54], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[56], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[57], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[58], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[60], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[80], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[81], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[82], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[83], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[84], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[85], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[86], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[88], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[89], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[90], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[92], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[112], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[113], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[114], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[116], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[120], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[48], aux,	prime);
  DFT2_aux_gmp(x[0], x[16], aux, prime);
  DFT2_aux_gmp(x[1], x[17], aux, prime);
  DFT2_aux_gmp(x[2], x[18], aux, prime);
  DFT2_aux_gmp(x[3], x[19], aux, prime);
  DFT2_aux_gmp(x[4], x[20], aux, prime);
  DFT2_aux_gmp(x[5], x[21], aux, prime);
  DFT2_aux_gmp(x[6], x[22], aux, prime);
  DFT2_aux_gmp(x[7], x[23], aux, prime);
  DFT2_aux_gmp(x[8], x[24], aux, prime);
  DFT2_aux_gmp(x[9], x[25], aux, prime);
  DFT2_aux_gmp(x[10], x[26], aux, prime);
  DFT2_aux_gmp(x[11], x[27], aux, prime);
  DFT2_aux_gmp(x[12], x[28], aux, prime);
  DFT2_aux_gmp(x[13], x[29], aux, prime);
  DFT2_aux_gmp(x[14], x[30], aux, prime);
  DFT2_aux_gmp(x[15], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[48], aux, prime);
  DFT2_aux_gmp(x[33], x[49], aux, prime);
  DFT2_aux_gmp(x[34], x[50], aux, prime);
  DFT2_aux_gmp(x[35], x[51], aux, prime);
  DFT2_aux_gmp(x[36], x[52], aux, prime);
  DFT2_aux_gmp(x[37], x[53], aux, prime);
  DFT2_aux_gmp(x[38], x[54], aux, prime);
  DFT2_aux_gmp(x[39], x[55], aux, prime);
  DFT2_aux_gmp(x[40], x[56], aux, prime);
  DFT2_aux_gmp(x[41], x[57], aux, prime);
  DFT2_aux_gmp(x[42], x[58], aux, prime);
  DFT2_aux_gmp(x[43], x[59], aux, prime);
  DFT2_aux_gmp(x[44], x[60], aux, prime);
  DFT2_aux_gmp(x[45], x[61], aux, prime);
  DFT2_aux_gmp(x[46], x[62], aux, prime);
  DFT2_aux_gmp(x[47], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[80], aux, prime);
  DFT2_aux_gmp(x[65], x[81], aux, prime);
  DFT2_aux_gmp(x[66], x[82], aux, prime);
  DFT2_aux_gmp(x[67], x[83], aux, prime);
  DFT2_aux_gmp(x[68], x[84], aux, prime);
  DFT2_aux_gmp(x[69], x[85], aux, prime);
  DFT2_aux_gmp(x[70], x[86], aux, prime);
  DFT2_aux_gmp(x[71], x[87], aux, prime);
  DFT2_aux_gmp(x[72], x[88], aux, prime);
  DFT2_aux_gmp(x[73], x[89], aux, prime);
  DFT2_aux_gmp(x[74], x[90], aux, prime);
  DFT2_aux_gmp(x[75], x[91], aux, prime);
  DFT2_aux_gmp(x[76], x[92], aux, prime);
  DFT2_aux_gmp(x[77], x[93], aux, prime);
  DFT2_aux_gmp(x[78], x[94], aux, prime);
  DFT2_aux_gmp(x[79], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[112], aux, prime);
  DFT2_aux_gmp(x[97], x[113], aux, prime);
  DFT2_aux_gmp(x[98], x[114], aux, prime);
  DFT2_aux_gmp(x[99], x[115], aux, prime);
  DFT2_aux_gmp(x[100], x[116], aux, prime);
  DFT2_aux_gmp(x[101], x[117], aux, prime);
  DFT2_aux_gmp(x[102], x[118], aux, prime);
  DFT2_aux_gmp(x[103], x[119], aux, prime);
  DFT2_aux_gmp(x[104], x[120], aux, prime);
  DFT2_aux_gmp(x[105], x[121], aux, prime);
  DFT2_aux_gmp(x[106], x[122], aux, prime);
  DFT2_aux_gmp(x[107], x[123], aux, prime);
  DFT2_aux_gmp(x[108], x[124], aux, prime);
  DFT2_aux_gmp(x[109], x[125], aux, prime);
  DFT2_aux_gmp(x[110], x[126], aux, prime);
  DFT2_aux_gmp(x[111], x[127], aux, prime);
  mult_aux_gmp(x[24], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[25], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[26], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[27], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[28], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[29], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[30], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[40], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[41], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[42], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[43], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[44], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[45], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[46], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[56], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[57], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[58], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[60], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[72], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[73], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[74], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[75], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[76], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[77], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[78], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[88], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[89], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[90], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[92], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[104], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[105], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[106], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[108], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[120], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[56], aux,	prime);
  DFT2_aux_gmp(x[0], x[8], aux, prime);
  DFT2_aux_gmp(x[1], x[9], aux, prime);
  DFT2_aux_gmp(x[2], x[10], aux, prime);
  DFT2_aux_gmp(x[3], x[11], aux, prime);
  DFT2_aux_gmp(x[4], x[12], aux, prime);
  DFT2_aux_gmp(x[5], x[13], aux, prime);
  DFT2_aux_gmp(x[6], x[14], aux, prime);
  DFT2_aux_gmp(x[7], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[24], aux, prime);
  DFT2_aux_gmp(x[17], x[25], aux, prime);
  DFT2_aux_gmp(x[18], x[26], aux, prime);
  DFT2_aux_gmp(x[19], x[27], aux, prime);
  DFT2_aux_gmp(x[20], x[28], aux, prime);
  DFT2_aux_gmp(x[21], x[29], aux, prime);
  DFT2_aux_gmp(x[22], x[30], aux, prime);
  DFT2_aux_gmp(x[23], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[40], aux, prime);
  DFT2_aux_gmp(x[33], x[41], aux, prime);
  DFT2_aux_gmp(x[34], x[42], aux, prime);
  DFT2_aux_gmp(x[35], x[43], aux, prime);
  DFT2_aux_gmp(x[36], x[44], aux, prime);
  DFT2_aux_gmp(x[37], x[45], aux, prime);
  DFT2_aux_gmp(x[38], x[46], aux, prime);
  DFT2_aux_gmp(x[39], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[56], aux, prime);
  DFT2_aux_gmp(x[49], x[57], aux, prime);
  DFT2_aux_gmp(x[50], x[58], aux, prime);
  DFT2_aux_gmp(x[51], x[59], aux, prime);
  DFT2_aux_gmp(x[52], x[60], aux, prime);
  DFT2_aux_gmp(x[53], x[61], aux, prime);
  DFT2_aux_gmp(x[54], x[62], aux, prime);
  DFT2_aux_gmp(x[55], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[72], aux, prime);
  DFT2_aux_gmp(x[65], x[73], aux, prime);
  DFT2_aux_gmp(x[66], x[74], aux, prime);
  DFT2_aux_gmp(x[67], x[75], aux, prime);
  DFT2_aux_gmp(x[68], x[76], aux, prime);
  DFT2_aux_gmp(x[69], x[77], aux, prime);
  DFT2_aux_gmp(x[70], x[78], aux, prime);
  DFT2_aux_gmp(x[71], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[88], aux, prime);
  DFT2_aux_gmp(x[81], x[89], aux, prime);
  DFT2_aux_gmp(x[82], x[90], aux, prime);
  DFT2_aux_gmp(x[83], x[91], aux, prime);
  DFT2_aux_gmp(x[84], x[92], aux, prime);
  DFT2_aux_gmp(x[85], x[93], aux, prime);
  DFT2_aux_gmp(x[86], x[94], aux, prime);
  DFT2_aux_gmp(x[87], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[104], aux, prime);
  DFT2_aux_gmp(x[97], x[105], aux, prime);
  DFT2_aux_gmp(x[98], x[106], aux, prime);
  DFT2_aux_gmp(x[99], x[107], aux, prime);
  DFT2_aux_gmp(x[100], x[108], aux, prime);
  DFT2_aux_gmp(x[101], x[109], aux, prime);
  DFT2_aux_gmp(x[102], x[110], aux, prime);
  DFT2_aux_gmp(x[103], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[120], aux, prime);
  DFT2_aux_gmp(x[113], x[121], aux, prime);
  DFT2_aux_gmp(x[114], x[122], aux, prime);
  DFT2_aux_gmp(x[115], x[123], aux, prime);
  DFT2_aux_gmp(x[116], x[124], aux, prime);
  DFT2_aux_gmp(x[117], x[125], aux, prime);
  DFT2_aux_gmp(x[118], x[126], aux, prime);
  DFT2_aux_gmp(x[119], x[127], aux, prime);
  mult_aux_gmp(x[12], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[13], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[14], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[15], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[20], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[21], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[22], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[23], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[28], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[29], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[30], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[36], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[37], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[38], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[39], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[44], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[45], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[46], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[52], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[53], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[54], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[60], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[68], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[69], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[70], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[71], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[76], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[77], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[78], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[84], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[85], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[86], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[92], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[100], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[101], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[102], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[108], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[116], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[60], aux,	prime);
  DFT2_aux_gmp(x[0], x[4], aux, prime);
  DFT2_aux_gmp(x[1], x[5], aux, prime);
  DFT2_aux_gmp(x[2], x[6], aux, prime);
  DFT2_aux_gmp(x[3], x[7], aux, prime);
  DFT2_aux_gmp(x[8], x[12], aux, prime);
  DFT2_aux_gmp(x[9], x[13], aux, prime);
  DFT2_aux_gmp(x[10], x[14], aux, prime);
  DFT2_aux_gmp(x[11], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[20], aux, prime);
  DFT2_aux_gmp(x[17], x[21], aux, prime);
  DFT2_aux_gmp(x[18], x[22], aux, prime);
  DFT2_aux_gmp(x[19], x[23], aux, prime);
  DFT2_aux_gmp(x[24], x[28], aux, prime);
  DFT2_aux_gmp(x[25], x[29], aux, prime);
  DFT2_aux_gmp(x[26], x[30], aux, prime);
  DFT2_aux_gmp(x[27], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[36], aux, prime);
  DFT2_aux_gmp(x[33], x[37], aux, prime);
  DFT2_aux_gmp(x[34], x[38], aux, prime);
  DFT2_aux_gmp(x[35], x[39], aux, prime);
  DFT2_aux_gmp(x[40], x[44], aux, prime);
  DFT2_aux_gmp(x[41], x[45], aux, prime);
  DFT2_aux_gmp(x[42], x[46], aux, prime);
  DFT2_aux_gmp(x[43], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[52], aux, prime);
  DFT2_aux_gmp(x[49], x[53], aux, prime);
  DFT2_aux_gmp(x[50], x[54], aux, prime);
  DFT2_aux_gmp(x[51], x[55], aux, prime);
  DFT2_aux_gmp(x[56], x[60], aux, prime);
  DFT2_aux_gmp(x[57], x[61], aux, prime);
  DFT2_aux_gmp(x[58], x[62], aux, prime);
  DFT2_aux_gmp(x[59], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[68], aux, prime);
  DFT2_aux_gmp(x[65], x[69], aux, prime);
  DFT2_aux_gmp(x[66], x[70], aux, prime);
  DFT2_aux_gmp(x[67], x[71], aux, prime);
  DFT2_aux_gmp(x[72], x[76], aux, prime);
  DFT2_aux_gmp(x[73], x[77], aux, prime);
  DFT2_aux_gmp(x[74], x[78], aux, prime);
  DFT2_aux_gmp(x[75], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[84], aux, prime);
  DFT2_aux_gmp(x[81], x[85], aux, prime);
  DFT2_aux_gmp(x[82], x[86], aux, prime);
  DFT2_aux_gmp(x[83], x[87], aux, prime);
  DFT2_aux_gmp(x[88], x[92], aux, prime);
  DFT2_aux_gmp(x[89], x[93], aux, prime);
  DFT2_aux_gmp(x[90], x[94], aux, prime);
  DFT2_aux_gmp(x[91], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[100], aux, prime);
  DFT2_aux_gmp(x[97], x[101], aux, prime);
  DFT2_aux_gmp(x[98], x[102], aux, prime);
  DFT2_aux_gmp(x[99], x[103], aux, prime);
  DFT2_aux_gmp(x[104], x[108], aux, prime);
  DFT2_aux_gmp(x[105], x[109], aux, prime);
  DFT2_aux_gmp(x[106], x[110], aux, prime);
  DFT2_aux_gmp(x[107], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[116], aux, prime);
  DFT2_aux_gmp(x[113], x[117], aux, prime);
  DFT2_aux_gmp(x[114], x[118], aux, prime);
  DFT2_aux_gmp(x[115], x[119], aux, prime);
  DFT2_aux_gmp(x[120], x[124], aux, prime);
  DFT2_aux_gmp(x[121], x[125], aux, prime);
  DFT2_aux_gmp(x[122], x[126], aux, prime);
  DFT2_aux_gmp(x[123], x[127], aux, prime);
  mult_aux_gmp(x[6], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[7], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[10], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[11], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[14], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[15], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[18], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[19], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[22], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[23], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[26], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[27], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[30], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[34], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[35], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[38], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[39], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[42], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[43], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[46], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[50], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[51], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[54], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[58], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[66], precompute_pow_omega[2], aux,	prime);
  mult_aux_gmp(x[67], precompute_pow_omega[2], aux,	prime);
  mult_aux_gmp(x[70], precompute_pow_omega[34], aux,	prime);
  mult_aux_gmp(x[71], precompute_pow_omega[34], aux,	prime);
  mult_aux_gmp(x[74], precompute_pow_omega[18], aux,	prime);
  mult_aux_gmp(x[75], precompute_pow_omega[18], aux,	prime);
  mult_aux_gmp(x[78], precompute_pow_omega[50], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[50], aux,	prime);
  mult_aux_gmp(x[82], precompute_pow_omega[10], aux,	prime);
  mult_aux_gmp(x[83], precompute_pow_omega[10], aux,	prime);
  mult_aux_gmp(x[86], precompute_pow_omega[42], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[42], aux,	prime);
  mult_aux_gmp(x[90], precompute_pow_omega[26], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[26], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[58], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[58], aux,	prime);
  mult_aux_gmp(x[98], precompute_pow_omega[6], aux,	prime);
  mult_aux_gmp(x[99], precompute_pow_omega[6], aux,	prime);
  mult_aux_gmp(x[102], precompute_pow_omega[38], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[38], aux,	prime);
  mult_aux_gmp(x[106], precompute_pow_omega[22], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[22], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[54], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[54], aux,	prime);
  mult_aux_gmp(x[114], precompute_pow_omega[14], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[14], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[46], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[46], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[30], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[30], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[62], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[62], aux,	prime);
  DFT2_aux_gmp(x[0], x[2], aux, prime);
  DFT2_aux_gmp(x[1], x[3], aux, prime);
  DFT2_aux_gmp(x[4], x[6], aux, prime);
  DFT2_aux_gmp(x[5], x[7], aux, prime);
  DFT2_aux_gmp(x[8], x[10], aux, prime);
  DFT2_aux_gmp(x[9], x[11], aux, prime);
  DFT2_aux_gmp(x[12], x[14], aux, prime);
  DFT2_aux_gmp(x[13], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[18], aux, prime);
  DFT2_aux_gmp(x[17], x[19], aux, prime);
  DFT2_aux_gmp(x[20], x[22], aux, prime);
  DFT2_aux_gmp(x[21], x[23], aux, prime);
  DFT2_aux_gmp(x[24], x[26], aux, prime);
  DFT2_aux_gmp(x[25], x[27], aux, prime);
  DFT2_aux_gmp(x[28], x[30], aux, prime);
  DFT2_aux_gmp(x[29], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[34], aux, prime);
  DFT2_aux_gmp(x[33], x[35], aux, prime);
  DFT2_aux_gmp(x[36], x[38], aux, prime);
  DFT2_aux_gmp(x[37], x[39], aux, prime);
  DFT2_aux_gmp(x[40], x[42], aux, prime);
  DFT2_aux_gmp(x[41], x[43], aux, prime);
  DFT2_aux_gmp(x[44], x[46], aux, prime);
  DFT2_aux_gmp(x[45], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[50], aux, prime);
  DFT2_aux_gmp(x[49], x[51], aux, prime);
  DFT2_aux_gmp(x[52], x[54], aux, prime);
  DFT2_aux_gmp(x[53], x[55], aux, prime);
  DFT2_aux_gmp(x[56], x[58], aux, prime);
  DFT2_aux_gmp(x[57], x[59], aux, prime);
  DFT2_aux_gmp(x[60], x[62], aux, prime);
  DFT2_aux_gmp(x[61], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[66], aux, prime);
  DFT2_aux_gmp(x[65], x[67], aux, prime);
  DFT2_aux_gmp(x[68], x[70], aux, prime);
  DFT2_aux_gmp(x[69], x[71], aux, prime);
  DFT2_aux_gmp(x[72], x[74], aux, prime);
  DFT2_aux_gmp(x[73], x[75], aux, prime);
  DFT2_aux_gmp(x[76], x[78], aux, prime);
  DFT2_aux_gmp(x[77], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[82], aux, prime);
  DFT2_aux_gmp(x[81], x[83], aux, prime);
  DFT2_aux_gmp(x[84], x[86], aux, prime);
  DFT2_aux_gmp(x[85], x[87], aux, prime);
  DFT2_aux_gmp(x[88], x[90], aux, prime);
  DFT2_aux_gmp(x[89], x[91], aux, prime);
  DFT2_aux_gmp(x[92], x[94], aux, prime);
  DFT2_aux_gmp(x[93], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[98], aux, prime);
  DFT2_aux_gmp(x[97], x[99], aux, prime);
  DFT2_aux_gmp(x[100], x[102], aux, prime);
  DFT2_aux_gmp(x[101], x[103], aux, prime);
  DFT2_aux_gmp(x[104], x[106], aux, prime);
  DFT2_aux_gmp(x[105], x[107], aux, prime);
  DFT2_aux_gmp(x[108], x[110], aux, prime);
  DFT2_aux_gmp(x[109], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[114], aux, prime);
  DFT2_aux_gmp(x[113], x[115], aux, prime);
  DFT2_aux_gmp(x[116], x[118], aux, prime);
  DFT2_aux_gmp(x[117], x[119], aux, prime);
  DFT2_aux_gmp(x[120], x[122], aux, prime);
  DFT2_aux_gmp(x[121], x[123], aux, prime);
  DFT2_aux_gmp(x[124], x[126], aux, prime);
  DFT2_aux_gmp(x[125], x[127], aux, prime);
  mult_aux_gmp(x[3], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[5], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[7], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[9], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[11], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[13], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[15], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[17], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[19], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[21], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[23], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[25], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[27], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[29], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[33], precompute_pow_omega[2], aux,	prime);
  mult_aux_gmp(x[35], precompute_pow_omega[34], aux,	prime);
  mult_aux_gmp(x[37], precompute_pow_omega[18], aux,	prime);
  mult_aux_gmp(x[39], precompute_pow_omega[50], aux,	prime);
  mult_aux_gmp(x[41], precompute_pow_omega[10], aux,	prime);
  mult_aux_gmp(x[43], precompute_pow_omega[42], aux,	prime);
  mult_aux_gmp(x[45], precompute_pow_omega[26], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[58], aux,	prime);
  mult_aux_gmp(x[49], precompute_pow_omega[6], aux,	prime);
  mult_aux_gmp(x[51], precompute_pow_omega[38], aux,	prime);
  mult_aux_gmp(x[53], precompute_pow_omega[22], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[54], aux,	prime);
  mult_aux_gmp(x[57], precompute_pow_omega[14], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[46], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[30], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[62], aux,	prime);
  mult_aux_gmp(x[65], precompute_pow_omega[1], aux,	prime);
  mult_aux_gmp(x[67], precompute_pow_omega[33], aux,	prime);
  mult_aux_gmp(x[69], precompute_pow_omega[17], aux,	prime);
  mult_aux_gmp(x[71], precompute_pow_omega[49], aux,	prime);
  mult_aux_gmp(x[73], precompute_pow_omega[9], aux,	prime);
  mult_aux_gmp(x[75], precompute_pow_omega[41], aux,	prime);
  mult_aux_gmp(x[77], precompute_pow_omega[25], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[57], aux,	prime);
  mult_aux_gmp(x[81], precompute_pow_omega[5], aux,	prime);
  mult_aux_gmp(x[83], precompute_pow_omega[37], aux,	prime);
  mult_aux_gmp(x[85], precompute_pow_omega[21], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[53], aux,	prime);
  mult_aux_gmp(x[89], precompute_pow_omega[13], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[45], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[29], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[61], aux,	prime);
  mult_aux_gmp(x[97], precompute_pow_omega[3], aux,	prime);
  mult_aux_gmp(x[99], precompute_pow_omega[35], aux,	prime);
  mult_aux_gmp(x[101], precompute_pow_omega[19], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[51], aux,	prime);
  mult_aux_gmp(x[105], precompute_pow_omega[11], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[43], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[27], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[59], aux,	prime);
  mult_aux_gmp(x[113], precompute_pow_omega[7], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[39], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[23], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[55], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[15], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[47], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[31], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[63], aux,	prime);
  DFT2_aux_gmp(x[0], x[1], aux, prime);
  DFT2_aux_gmp(x[2], x[3], aux, prime);
  DFT2_aux_gmp(x[4], x[5], aux, prime);
  DFT2_aux_gmp(x[6], x[7], aux, prime);
  DFT2_aux_gmp(x[8], x[9], aux, prime);
  DFT2_aux_gmp(x[10], x[11], aux, prime);
  DFT2_aux_gmp(x[12], x[13], aux, prime);
  DFT2_aux_gmp(x[14], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[17], aux, prime);
  DFT2_aux_gmp(x[18], x[19], aux, prime);
  DFT2_aux_gmp(x[20], x[21], aux, prime);
  DFT2_aux_gmp(x[22], x[23], aux, prime);
  DFT2_aux_gmp(x[24], x[25], aux, prime);
  DFT2_aux_gmp(x[26], x[27], aux, prime);
  DFT2_aux_gmp(x[28], x[29], aux, prime);
  DFT2_aux_gmp(x[30], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[33], aux, prime);
  DFT2_aux_gmp(x[34], x[35], aux, prime);
  DFT2_aux_gmp(x[36], x[37], aux, prime);
  DFT2_aux_gmp(x[38], x[39], aux, prime);
  DFT2_aux_gmp(x[40], x[41], aux, prime);
  DFT2_aux_gmp(x[42], x[43], aux, prime);
  DFT2_aux_gmp(x[44], x[45], aux, prime);
  DFT2_aux_gmp(x[46], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[49], aux, prime);
  DFT2_aux_gmp(x[50], x[51], aux, prime);
  DFT2_aux_gmp(x[52], x[53], aux, prime);
  DFT2_aux_gmp(x[54], x[55], aux, prime);
  DFT2_aux_gmp(x[56], x[57], aux, prime);
  DFT2_aux_gmp(x[58], x[59], aux, prime);
  DFT2_aux_gmp(x[60], x[61], aux, prime);
  DFT2_aux_gmp(x[62], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[65], aux, prime);
  DFT2_aux_gmp(x[66], x[67], aux, prime);
  DFT2_aux_gmp(x[68], x[69], aux, prime);
  DFT2_aux_gmp(x[70], x[71], aux, prime);
  DFT2_aux_gmp(x[72], x[73], aux, prime);
  DFT2_aux_gmp(x[74], x[75], aux, prime);
  DFT2_aux_gmp(x[76], x[77], aux, prime);
  DFT2_aux_gmp(x[78], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[81], aux, prime);
  DFT2_aux_gmp(x[82], x[83], aux, prime);
  DFT2_aux_gmp(x[84], x[85], aux, prime);
  DFT2_aux_gmp(x[86], x[87], aux, prime);
  DFT2_aux_gmp(x[88], x[89], aux, prime);
  DFT2_aux_gmp(x[90], x[91], aux, prime);
  DFT2_aux_gmp(x[92], x[93], aux, prime);
  DFT2_aux_gmp(x[94], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[97], aux, prime);
  DFT2_aux_gmp(x[98], x[99], aux, prime);
  DFT2_aux_gmp(x[100], x[101], aux, prime);
  DFT2_aux_gmp(x[102], x[103], aux, prime);
  DFT2_aux_gmp(x[104], x[105], aux, prime);
  DFT2_aux_gmp(x[106], x[107], aux, prime);
  DFT2_aux_gmp(x[108], x[109], aux, prime);
  DFT2_aux_gmp(x[110], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[113], aux, prime);
  DFT2_aux_gmp(x[114], x[115], aux, prime);
  DFT2_aux_gmp(x[116], x[117], aux, prime);
  DFT2_aux_gmp(x[118], x[119], aux, prime);
  DFT2_aux_gmp(x[120], x[121], aux, prime);
  DFT2_aux_gmp(x[122], x[123], aux, prime);
  DFT2_aux_gmp(x[124], x[125], aux, prime);
  DFT2_aux_gmp(x[126], x[127], aux, prime);
  swap_aux_gmp(x[1], x[64], aux);
  swap_aux_gmp(x[2], x[32], aux);
  swap_aux_gmp(x[3], x[96], aux);
  swap_aux_gmp(x[4], x[16], aux);
  swap_aux_gmp(x[5], x[80], aux);
  swap_aux_gmp(x[6], x[48], aux);
  swap_aux_gmp(x[7], x[112], aux);
  swap_aux_gmp(x[9], x[72], aux);
  swap_aux_gmp(x[10], x[40], aux);
  swap_aux_gmp(x[11], x[104], aux);
  swap_aux_gmp(x[12], x[24], aux);
  swap_aux_gmp(x[13], x[88], aux);
  swap_aux_gmp(x[14], x[56], aux);
  swap_aux_gmp(x[15], x[120], aux);
  swap_aux_gmp(x[17], x[68], aux);
  swap_aux_gmp(x[18], x[36], aux);
  swap_aux_gmp(x[19], x[100], aux);
  swap_aux_gmp(x[21], x[84], aux);
  swap_aux_gmp(x[22], x[52], aux);
  swap_aux_gmp(x[23], x[116], aux);
  swap_aux_gmp(x[25], x[76], aux);
  swap_aux_gmp(x[26], x[44], aux);
  swap_aux_gmp(x[27], x[108], aux);
  swap_aux_gmp(x[29], x[92], aux);
  swap_aux_gmp(x[30], x[60], aux);
  swap_aux_gmp(x[31], x[124], aux);
  swap_aux_gmp(x[33], x[66], aux);
  swap_aux_gmp(x[35], x[98], aux);
  swap_aux_gmp(x[37], x[82], aux);
  swap_aux_gmp(x[38], x[50], aux);
  swap_aux_gmp(x[39], x[114], aux);
  swap_aux_gmp(x[41], x[74], aux);
  swap_aux_gmp(x[43], x[106], aux);
  swap_aux_gmp(x[45], x[90], aux);
  swap_aux_gmp(x[46], x[58], aux);
  swap_aux_gmp(x[47], x[122], aux);
  swap_aux_gmp(x[49], x[70], aux);
  swap_aux_gmp(x[51], x[102], aux);
  swap_aux_gmp(x[53], x[86], aux);
  swap_aux_gmp(x[55], x[118], aux);
  swap_aux_gmp(x[57], x[78], aux);
  swap_aux_gmp(x[59], x[110], aux);
  swap_aux_gmp(x[61], x[94], aux);
  swap_aux_gmp(x[63], x[126], aux);
  swap_aux_gmp(x[67], x[97], aux);
  swap_aux_gmp(x[69], x[81], aux);
  swap_aux_gmp(x[71], x[113], aux);
  swap_aux_gmp(x[75], x[105], aux);
  swap_aux_gmp(x[77], x[89], aux);
  swap_aux_gmp(x[79], x[121], aux);
  swap_aux_gmp(x[83], x[101], aux);
  swap_aux_gmp(x[87], x[117], aux);
  swap_aux_gmp(x[91], x[109], aux);
  swap_aux_gmp(x[95], x[125], aux);
  swap_aux_gmp(x[103], x[115], aux);
  swap_aux_gmp(x[111], x[123], aux);
}


/**************************************/

void
DFT_256_aux_gmp (mpz_t* x, mpz_t aux, const mpz_t omega, const mpz_t prime, const mpz_t * precompute_pow_omega)
{
#if PROFILING_GMP_ENABLED ==1
  inc_gmp_profiling_counter (&n_dft2_gmp_called);
#endif

  DFT2_aux_gmp(x[0], x[128], aux, prime);
  DFT2_aux_gmp(x[1], x[129], aux, prime);
  DFT2_aux_gmp(x[2], x[130], aux, prime);
  DFT2_aux_gmp(x[3], x[131], aux, prime);
  DFT2_aux_gmp(x[4], x[132], aux, prime);
  DFT2_aux_gmp(x[5], x[133], aux, prime);
  DFT2_aux_gmp(x[6], x[134], aux, prime);
  DFT2_aux_gmp(x[7], x[135], aux, prime);
  DFT2_aux_gmp(x[8], x[136], aux, prime);
  DFT2_aux_gmp(x[9], x[137], aux, prime);
  DFT2_aux_gmp(x[10], x[138], aux, prime);
  DFT2_aux_gmp(x[11], x[139], aux, prime);
  DFT2_aux_gmp(x[12], x[140], aux, prime);
  DFT2_aux_gmp(x[13], x[141], aux, prime);
  DFT2_aux_gmp(x[14], x[142], aux, prime);
  DFT2_aux_gmp(x[15], x[143], aux, prime);
  DFT2_aux_gmp(x[16], x[144], aux, prime);
  DFT2_aux_gmp(x[17], x[145], aux, prime);
  DFT2_aux_gmp(x[18], x[146], aux, prime);
  DFT2_aux_gmp(x[19], x[147], aux, prime);
  DFT2_aux_gmp(x[20], x[148], aux, prime);
  DFT2_aux_gmp(x[21], x[149], aux, prime);
  DFT2_aux_gmp(x[22], x[150], aux, prime);
  DFT2_aux_gmp(x[23], x[151], aux, prime);
  DFT2_aux_gmp(x[24], x[152], aux, prime);
  DFT2_aux_gmp(x[25], x[153], aux, prime);
  DFT2_aux_gmp(x[26], x[154], aux, prime);
  DFT2_aux_gmp(x[27], x[155], aux, prime);
  DFT2_aux_gmp(x[28], x[156], aux, prime);
  DFT2_aux_gmp(x[29], x[157], aux, prime);
  DFT2_aux_gmp(x[30], x[158], aux, prime);
  DFT2_aux_gmp(x[31], x[159], aux, prime);
  DFT2_aux_gmp(x[32], x[160], aux, prime);
  DFT2_aux_gmp(x[33], x[161], aux, prime);
  DFT2_aux_gmp(x[34], x[162], aux, prime);
  DFT2_aux_gmp(x[35], x[163], aux, prime);
  DFT2_aux_gmp(x[36], x[164], aux, prime);
  DFT2_aux_gmp(x[37], x[165], aux, prime);
  DFT2_aux_gmp(x[38], x[166], aux, prime);
  DFT2_aux_gmp(x[39], x[167], aux, prime);
  DFT2_aux_gmp(x[40], x[168], aux, prime);
  DFT2_aux_gmp(x[41], x[169], aux, prime);
  DFT2_aux_gmp(x[42], x[170], aux, prime);
  DFT2_aux_gmp(x[43], x[171], aux, prime);
  DFT2_aux_gmp(x[44], x[172], aux, prime);
  DFT2_aux_gmp(x[45], x[173], aux, prime);
  DFT2_aux_gmp(x[46], x[174], aux, prime);
  DFT2_aux_gmp(x[47], x[175], aux, prime);
  DFT2_aux_gmp(x[48], x[176], aux, prime);
  DFT2_aux_gmp(x[49], x[177], aux, prime);
  DFT2_aux_gmp(x[50], x[178], aux, prime);
  DFT2_aux_gmp(x[51], x[179], aux, prime);
  DFT2_aux_gmp(x[52], x[180], aux, prime);
  DFT2_aux_gmp(x[53], x[181], aux, prime);
  DFT2_aux_gmp(x[54], x[182], aux, prime);
  DFT2_aux_gmp(x[55], x[183], aux, prime);
  DFT2_aux_gmp(x[56], x[184], aux, prime);
  DFT2_aux_gmp(x[57], x[185], aux, prime);
  DFT2_aux_gmp(x[58], x[186], aux, prime);
  DFT2_aux_gmp(x[59], x[187], aux, prime);
  DFT2_aux_gmp(x[60], x[188], aux, prime);
  DFT2_aux_gmp(x[61], x[189], aux, prime);
  DFT2_aux_gmp(x[62], x[190], aux, prime);
  DFT2_aux_gmp(x[63], x[191], aux, prime);
  DFT2_aux_gmp(x[64], x[192], aux, prime);
  DFT2_aux_gmp(x[65], x[193], aux, prime);
  DFT2_aux_gmp(x[66], x[194], aux, prime);
  DFT2_aux_gmp(x[67], x[195], aux, prime);
  DFT2_aux_gmp(x[68], x[196], aux, prime);
  DFT2_aux_gmp(x[69], x[197], aux, prime);
  DFT2_aux_gmp(x[70], x[198], aux, prime);
  DFT2_aux_gmp(x[71], x[199], aux, prime);
  DFT2_aux_gmp(x[72], x[200], aux, prime);
  DFT2_aux_gmp(x[73], x[201], aux, prime);
  DFT2_aux_gmp(x[74], x[202], aux, prime);
  DFT2_aux_gmp(x[75], x[203], aux, prime);
  DFT2_aux_gmp(x[76], x[204], aux, prime);
  DFT2_aux_gmp(x[77], x[205], aux, prime);
  DFT2_aux_gmp(x[78], x[206], aux, prime);
  DFT2_aux_gmp(x[79], x[207], aux, prime);
  DFT2_aux_gmp(x[80], x[208], aux, prime);
  DFT2_aux_gmp(x[81], x[209], aux, prime);
  DFT2_aux_gmp(x[82], x[210], aux, prime);
  DFT2_aux_gmp(x[83], x[211], aux, prime);
  DFT2_aux_gmp(x[84], x[212], aux, prime);
  DFT2_aux_gmp(x[85], x[213], aux, prime);
  DFT2_aux_gmp(x[86], x[214], aux, prime);
  DFT2_aux_gmp(x[87], x[215], aux, prime);
  DFT2_aux_gmp(x[88], x[216], aux, prime);
  DFT2_aux_gmp(x[89], x[217], aux, prime);
  DFT2_aux_gmp(x[90], x[218], aux, prime);
  DFT2_aux_gmp(x[91], x[219], aux, prime);
  DFT2_aux_gmp(x[92], x[220], aux, prime);
  DFT2_aux_gmp(x[93], x[221], aux, prime);
  DFT2_aux_gmp(x[94], x[222], aux, prime);
  DFT2_aux_gmp(x[95], x[223], aux, prime);
  DFT2_aux_gmp(x[96], x[224], aux, prime);
  DFT2_aux_gmp(x[97], x[225], aux, prime);
  DFT2_aux_gmp(x[98], x[226], aux, prime);
  DFT2_aux_gmp(x[99], x[227], aux, prime);
  DFT2_aux_gmp(x[100], x[228], aux, prime);
  DFT2_aux_gmp(x[101], x[229], aux, prime);
  DFT2_aux_gmp(x[102], x[230], aux, prime);
  DFT2_aux_gmp(x[103], x[231], aux, prime);
  DFT2_aux_gmp(x[104], x[232], aux, prime);
  DFT2_aux_gmp(x[105], x[233], aux, prime);
  DFT2_aux_gmp(x[106], x[234], aux, prime);
  DFT2_aux_gmp(x[107], x[235], aux, prime);
  DFT2_aux_gmp(x[108], x[236], aux, prime);
  DFT2_aux_gmp(x[109], x[237], aux, prime);
  DFT2_aux_gmp(x[110], x[238], aux, prime);
  DFT2_aux_gmp(x[111], x[239], aux, prime);
  DFT2_aux_gmp(x[112], x[240], aux, prime);
  DFT2_aux_gmp(x[113], x[241], aux, prime);
  DFT2_aux_gmp(x[114], x[242], aux, prime);
  DFT2_aux_gmp(x[115], x[243], aux, prime);
  DFT2_aux_gmp(x[116], x[244], aux, prime);
  DFT2_aux_gmp(x[117], x[245], aux, prime);
  DFT2_aux_gmp(x[118], x[246], aux, prime);
  DFT2_aux_gmp(x[119], x[247], aux, prime);
  DFT2_aux_gmp(x[120], x[248], aux, prime);
  DFT2_aux_gmp(x[121], x[249], aux, prime);
  DFT2_aux_gmp(x[122], x[250], aux, prime);
  DFT2_aux_gmp(x[123], x[251], aux, prime);
  DFT2_aux_gmp(x[124], x[252], aux, prime);
  DFT2_aux_gmp(x[125], x[253], aux, prime);
  DFT2_aux_gmp(x[126], x[254], aux, prime);
  DFT2_aux_gmp(x[127], x[255], aux, prime);
  mult_aux_gmp(x[192], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[193], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[194], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[195], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[196], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[197], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[198], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[199], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[200], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[201], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[202], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[203], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[204], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[205], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[206], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[207], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[208], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[209], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[210], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[211], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[212], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[213], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[214], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[215], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[216], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[217], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[218], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[219], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[220], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[221], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[222], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[223], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[224], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[225], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[226], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[227], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[228], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[229], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[230], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[231], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[232], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[233], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[234], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[235], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[236], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[237], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[238], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[239], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[240], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[241], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[242], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[243], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[244], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[245], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[246], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[247], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[248], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[249], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[250], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[251], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[252], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[253], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[254], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[64], aux,	prime);
  DFT2_aux_gmp(x[0], x[64], aux, prime);
  DFT2_aux_gmp(x[1], x[65], aux, prime);
  DFT2_aux_gmp(x[2], x[66], aux, prime);
  DFT2_aux_gmp(x[3], x[67], aux, prime);
  DFT2_aux_gmp(x[4], x[68], aux, prime);
  DFT2_aux_gmp(x[5], x[69], aux, prime);
  DFT2_aux_gmp(x[6], x[70], aux, prime);
  DFT2_aux_gmp(x[7], x[71], aux, prime);
  DFT2_aux_gmp(x[8], x[72], aux, prime);
  DFT2_aux_gmp(x[9], x[73], aux, prime);
  DFT2_aux_gmp(x[10], x[74], aux, prime);
  DFT2_aux_gmp(x[11], x[75], aux, prime);
  DFT2_aux_gmp(x[12], x[76], aux, prime);
  DFT2_aux_gmp(x[13], x[77], aux, prime);
  DFT2_aux_gmp(x[14], x[78], aux, prime);
  DFT2_aux_gmp(x[15], x[79], aux, prime);
  DFT2_aux_gmp(x[16], x[80], aux, prime);
  DFT2_aux_gmp(x[17], x[81], aux, prime);
  DFT2_aux_gmp(x[18], x[82], aux, prime);
  DFT2_aux_gmp(x[19], x[83], aux, prime);
  DFT2_aux_gmp(x[20], x[84], aux, prime);
  DFT2_aux_gmp(x[21], x[85], aux, prime);
  DFT2_aux_gmp(x[22], x[86], aux, prime);
  DFT2_aux_gmp(x[23], x[87], aux, prime);
  DFT2_aux_gmp(x[24], x[88], aux, prime);
  DFT2_aux_gmp(x[25], x[89], aux, prime);
  DFT2_aux_gmp(x[26], x[90], aux, prime);
  DFT2_aux_gmp(x[27], x[91], aux, prime);
  DFT2_aux_gmp(x[28], x[92], aux, prime);
  DFT2_aux_gmp(x[29], x[93], aux, prime);
  DFT2_aux_gmp(x[30], x[94], aux, prime);
  DFT2_aux_gmp(x[31], x[95], aux, prime);
  DFT2_aux_gmp(x[32], x[96], aux, prime);
  DFT2_aux_gmp(x[33], x[97], aux, prime);
  DFT2_aux_gmp(x[34], x[98], aux, prime);
  DFT2_aux_gmp(x[35], x[99], aux, prime);
  DFT2_aux_gmp(x[36], x[100], aux, prime);
  DFT2_aux_gmp(x[37], x[101], aux, prime);
  DFT2_aux_gmp(x[38], x[102], aux, prime);
  DFT2_aux_gmp(x[39], x[103], aux, prime);
  DFT2_aux_gmp(x[40], x[104], aux, prime);
  DFT2_aux_gmp(x[41], x[105], aux, prime);
  DFT2_aux_gmp(x[42], x[106], aux, prime);
  DFT2_aux_gmp(x[43], x[107], aux, prime);
  DFT2_aux_gmp(x[44], x[108], aux, prime);
  DFT2_aux_gmp(x[45], x[109], aux, prime);
  DFT2_aux_gmp(x[46], x[110], aux, prime);
  DFT2_aux_gmp(x[47], x[111], aux, prime);
  DFT2_aux_gmp(x[48], x[112], aux, prime);
  DFT2_aux_gmp(x[49], x[113], aux, prime);
  DFT2_aux_gmp(x[50], x[114], aux, prime);
  DFT2_aux_gmp(x[51], x[115], aux, prime);
  DFT2_aux_gmp(x[52], x[116], aux, prime);
  DFT2_aux_gmp(x[53], x[117], aux, prime);
  DFT2_aux_gmp(x[54], x[118], aux, prime);
  DFT2_aux_gmp(x[55], x[119], aux, prime);
  DFT2_aux_gmp(x[56], x[120], aux, prime);
  DFT2_aux_gmp(x[57], x[121], aux, prime);
  DFT2_aux_gmp(x[58], x[122], aux, prime);
  DFT2_aux_gmp(x[59], x[123], aux, prime);
  DFT2_aux_gmp(x[60], x[124], aux, prime);
  DFT2_aux_gmp(x[61], x[125], aux, prime);
  DFT2_aux_gmp(x[62], x[126], aux, prime);
  DFT2_aux_gmp(x[63], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[192], aux, prime);
  DFT2_aux_gmp(x[129], x[193], aux, prime);
  DFT2_aux_gmp(x[130], x[194], aux, prime);
  DFT2_aux_gmp(x[131], x[195], aux, prime);
  DFT2_aux_gmp(x[132], x[196], aux, prime);
  DFT2_aux_gmp(x[133], x[197], aux, prime);
  DFT2_aux_gmp(x[134], x[198], aux, prime);
  DFT2_aux_gmp(x[135], x[199], aux, prime);
  DFT2_aux_gmp(x[136], x[200], aux, prime);
  DFT2_aux_gmp(x[137], x[201], aux, prime);
  DFT2_aux_gmp(x[138], x[202], aux, prime);
  DFT2_aux_gmp(x[139], x[203], aux, prime);
  DFT2_aux_gmp(x[140], x[204], aux, prime);
  DFT2_aux_gmp(x[141], x[205], aux, prime);
  DFT2_aux_gmp(x[142], x[206], aux, prime);
  DFT2_aux_gmp(x[143], x[207], aux, prime);
  DFT2_aux_gmp(x[144], x[208], aux, prime);
  DFT2_aux_gmp(x[145], x[209], aux, prime);
  DFT2_aux_gmp(x[146], x[210], aux, prime);
  DFT2_aux_gmp(x[147], x[211], aux, prime);
  DFT2_aux_gmp(x[148], x[212], aux, prime);
  DFT2_aux_gmp(x[149], x[213], aux, prime);
  DFT2_aux_gmp(x[150], x[214], aux, prime);
  DFT2_aux_gmp(x[151], x[215], aux, prime);
  DFT2_aux_gmp(x[152], x[216], aux, prime);
  DFT2_aux_gmp(x[153], x[217], aux, prime);
  DFT2_aux_gmp(x[154], x[218], aux, prime);
  DFT2_aux_gmp(x[155], x[219], aux, prime);
  DFT2_aux_gmp(x[156], x[220], aux, prime);
  DFT2_aux_gmp(x[157], x[221], aux, prime);
  DFT2_aux_gmp(x[158], x[222], aux, prime);
  DFT2_aux_gmp(x[159], x[223], aux, prime);
  DFT2_aux_gmp(x[160], x[224], aux, prime);
  DFT2_aux_gmp(x[161], x[225], aux, prime);
  DFT2_aux_gmp(x[162], x[226], aux, prime);
  DFT2_aux_gmp(x[163], x[227], aux, prime);
  DFT2_aux_gmp(x[164], x[228], aux, prime);
  DFT2_aux_gmp(x[165], x[229], aux, prime);
  DFT2_aux_gmp(x[166], x[230], aux, prime);
  DFT2_aux_gmp(x[167], x[231], aux, prime);
  DFT2_aux_gmp(x[168], x[232], aux, prime);
  DFT2_aux_gmp(x[169], x[233], aux, prime);
  DFT2_aux_gmp(x[170], x[234], aux, prime);
  DFT2_aux_gmp(x[171], x[235], aux, prime);
  DFT2_aux_gmp(x[172], x[236], aux, prime);
  DFT2_aux_gmp(x[173], x[237], aux, prime);
  DFT2_aux_gmp(x[174], x[238], aux, prime);
  DFT2_aux_gmp(x[175], x[239], aux, prime);
  DFT2_aux_gmp(x[176], x[240], aux, prime);
  DFT2_aux_gmp(x[177], x[241], aux, prime);
  DFT2_aux_gmp(x[178], x[242], aux, prime);
  DFT2_aux_gmp(x[179], x[243], aux, prime);
  DFT2_aux_gmp(x[180], x[244], aux, prime);
  DFT2_aux_gmp(x[181], x[245], aux, prime);
  DFT2_aux_gmp(x[182], x[246], aux, prime);
  DFT2_aux_gmp(x[183], x[247], aux, prime);
  DFT2_aux_gmp(x[184], x[248], aux, prime);
  DFT2_aux_gmp(x[185], x[249], aux, prime);
  DFT2_aux_gmp(x[186], x[250], aux, prime);
  DFT2_aux_gmp(x[187], x[251], aux, prime);
  DFT2_aux_gmp(x[188], x[252], aux, prime);
  DFT2_aux_gmp(x[189], x[253], aux, prime);
  DFT2_aux_gmp(x[190], x[254], aux, prime);
  DFT2_aux_gmp(x[191], x[255], aux, prime);
  mult_aux_gmp(x[96], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[97], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[98], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[99], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[100], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[101], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[102], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[104], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[105], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[106], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[108], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[112], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[113], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[114], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[116], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[120], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[160], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[161], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[162], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[163], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[164], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[165], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[166], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[167], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[168], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[169], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[170], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[171], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[172], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[173], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[174], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[175], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[176], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[177], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[178], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[179], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[180], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[181], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[182], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[183], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[184], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[185], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[186], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[187], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[188], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[189], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[190], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[191], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[224], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[225], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[226], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[227], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[228], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[229], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[230], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[231], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[232], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[233], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[234], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[235], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[236], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[237], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[238], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[239], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[240], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[241], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[242], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[243], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[244], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[245], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[246], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[247], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[248], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[249], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[250], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[251], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[252], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[253], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[254], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[96], aux,	prime);
  DFT2_aux_gmp(x[0], x[32], aux, prime);
  DFT2_aux_gmp(x[1], x[33], aux, prime);
  DFT2_aux_gmp(x[2], x[34], aux, prime);
  DFT2_aux_gmp(x[3], x[35], aux, prime);
  DFT2_aux_gmp(x[4], x[36], aux, prime);
  DFT2_aux_gmp(x[5], x[37], aux, prime);
  DFT2_aux_gmp(x[6], x[38], aux, prime);
  DFT2_aux_gmp(x[7], x[39], aux, prime);
  DFT2_aux_gmp(x[8], x[40], aux, prime);
  DFT2_aux_gmp(x[9], x[41], aux, prime);
  DFT2_aux_gmp(x[10], x[42], aux, prime);
  DFT2_aux_gmp(x[11], x[43], aux, prime);
  DFT2_aux_gmp(x[12], x[44], aux, prime);
  DFT2_aux_gmp(x[13], x[45], aux, prime);
  DFT2_aux_gmp(x[14], x[46], aux, prime);
  DFT2_aux_gmp(x[15], x[47], aux, prime);
  DFT2_aux_gmp(x[16], x[48], aux, prime);
  DFT2_aux_gmp(x[17], x[49], aux, prime);
  DFT2_aux_gmp(x[18], x[50], aux, prime);
  DFT2_aux_gmp(x[19], x[51], aux, prime);
  DFT2_aux_gmp(x[20], x[52], aux, prime);
  DFT2_aux_gmp(x[21], x[53], aux, prime);
  DFT2_aux_gmp(x[22], x[54], aux, prime);
  DFT2_aux_gmp(x[23], x[55], aux, prime);
  DFT2_aux_gmp(x[24], x[56], aux, prime);
  DFT2_aux_gmp(x[25], x[57], aux, prime);
  DFT2_aux_gmp(x[26], x[58], aux, prime);
  DFT2_aux_gmp(x[27], x[59], aux, prime);
  DFT2_aux_gmp(x[28], x[60], aux, prime);
  DFT2_aux_gmp(x[29], x[61], aux, prime);
  DFT2_aux_gmp(x[30], x[62], aux, prime);
  DFT2_aux_gmp(x[31], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[96], aux, prime);
  DFT2_aux_gmp(x[65], x[97], aux, prime);
  DFT2_aux_gmp(x[66], x[98], aux, prime);
  DFT2_aux_gmp(x[67], x[99], aux, prime);
  DFT2_aux_gmp(x[68], x[100], aux, prime);
  DFT2_aux_gmp(x[69], x[101], aux, prime);
  DFT2_aux_gmp(x[70], x[102], aux, prime);
  DFT2_aux_gmp(x[71], x[103], aux, prime);
  DFT2_aux_gmp(x[72], x[104], aux, prime);
  DFT2_aux_gmp(x[73], x[105], aux, prime);
  DFT2_aux_gmp(x[74], x[106], aux, prime);
  DFT2_aux_gmp(x[75], x[107], aux, prime);
  DFT2_aux_gmp(x[76], x[108], aux, prime);
  DFT2_aux_gmp(x[77], x[109], aux, prime);
  DFT2_aux_gmp(x[78], x[110], aux, prime);
  DFT2_aux_gmp(x[79], x[111], aux, prime);
  DFT2_aux_gmp(x[80], x[112], aux, prime);
  DFT2_aux_gmp(x[81], x[113], aux, prime);
  DFT2_aux_gmp(x[82], x[114], aux, prime);
  DFT2_aux_gmp(x[83], x[115], aux, prime);
  DFT2_aux_gmp(x[84], x[116], aux, prime);
  DFT2_aux_gmp(x[85], x[117], aux, prime);
  DFT2_aux_gmp(x[86], x[118], aux, prime);
  DFT2_aux_gmp(x[87], x[119], aux, prime);
  DFT2_aux_gmp(x[88], x[120], aux, prime);
  DFT2_aux_gmp(x[89], x[121], aux, prime);
  DFT2_aux_gmp(x[90], x[122], aux, prime);
  DFT2_aux_gmp(x[91], x[123], aux, prime);
  DFT2_aux_gmp(x[92], x[124], aux, prime);
  DFT2_aux_gmp(x[93], x[125], aux, prime);
  DFT2_aux_gmp(x[94], x[126], aux, prime);
  DFT2_aux_gmp(x[95], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[160], aux, prime);
  DFT2_aux_gmp(x[129], x[161], aux, prime);
  DFT2_aux_gmp(x[130], x[162], aux, prime);
  DFT2_aux_gmp(x[131], x[163], aux, prime);
  DFT2_aux_gmp(x[132], x[164], aux, prime);
  DFT2_aux_gmp(x[133], x[165], aux, prime);
  DFT2_aux_gmp(x[134], x[166], aux, prime);
  DFT2_aux_gmp(x[135], x[167], aux, prime);
  DFT2_aux_gmp(x[136], x[168], aux, prime);
  DFT2_aux_gmp(x[137], x[169], aux, prime);
  DFT2_aux_gmp(x[138], x[170], aux, prime);
  DFT2_aux_gmp(x[139], x[171], aux, prime);
  DFT2_aux_gmp(x[140], x[172], aux, prime);
  DFT2_aux_gmp(x[141], x[173], aux, prime);
  DFT2_aux_gmp(x[142], x[174], aux, prime);
  DFT2_aux_gmp(x[143], x[175], aux, prime);
  DFT2_aux_gmp(x[144], x[176], aux, prime);
  DFT2_aux_gmp(x[145], x[177], aux, prime);
  DFT2_aux_gmp(x[146], x[178], aux, prime);
  DFT2_aux_gmp(x[147], x[179], aux, prime);
  DFT2_aux_gmp(x[148], x[180], aux, prime);
  DFT2_aux_gmp(x[149], x[181], aux, prime);
  DFT2_aux_gmp(x[150], x[182], aux, prime);
  DFT2_aux_gmp(x[151], x[183], aux, prime);
  DFT2_aux_gmp(x[152], x[184], aux, prime);
  DFT2_aux_gmp(x[153], x[185], aux, prime);
  DFT2_aux_gmp(x[154], x[186], aux, prime);
  DFT2_aux_gmp(x[155], x[187], aux, prime);
  DFT2_aux_gmp(x[156], x[188], aux, prime);
  DFT2_aux_gmp(x[157], x[189], aux, prime);
  DFT2_aux_gmp(x[158], x[190], aux, prime);
  DFT2_aux_gmp(x[159], x[191], aux, prime);
  DFT2_aux_gmp(x[192], x[224], aux, prime);
  DFT2_aux_gmp(x[193], x[225], aux, prime);
  DFT2_aux_gmp(x[194], x[226], aux, prime);
  DFT2_aux_gmp(x[195], x[227], aux, prime);
  DFT2_aux_gmp(x[196], x[228], aux, prime);
  DFT2_aux_gmp(x[197], x[229], aux, prime);
  DFT2_aux_gmp(x[198], x[230], aux, prime);
  DFT2_aux_gmp(x[199], x[231], aux, prime);
  DFT2_aux_gmp(x[200], x[232], aux, prime);
  DFT2_aux_gmp(x[201], x[233], aux, prime);
  DFT2_aux_gmp(x[202], x[234], aux, prime);
  DFT2_aux_gmp(x[203], x[235], aux, prime);
  DFT2_aux_gmp(x[204], x[236], aux, prime);
  DFT2_aux_gmp(x[205], x[237], aux, prime);
  DFT2_aux_gmp(x[206], x[238], aux, prime);
  DFT2_aux_gmp(x[207], x[239], aux, prime);
  DFT2_aux_gmp(x[208], x[240], aux, prime);
  DFT2_aux_gmp(x[209], x[241], aux, prime);
  DFT2_aux_gmp(x[210], x[242], aux, prime);
  DFT2_aux_gmp(x[211], x[243], aux, prime);
  DFT2_aux_gmp(x[212], x[244], aux, prime);
  DFT2_aux_gmp(x[213], x[245], aux, prime);
  DFT2_aux_gmp(x[214], x[246], aux, prime);
  DFT2_aux_gmp(x[215], x[247], aux, prime);
  DFT2_aux_gmp(x[216], x[248], aux, prime);
  DFT2_aux_gmp(x[217], x[249], aux, prime);
  DFT2_aux_gmp(x[218], x[250], aux, prime);
  DFT2_aux_gmp(x[219], x[251], aux, prime);
  DFT2_aux_gmp(x[220], x[252], aux, prime);
  DFT2_aux_gmp(x[221], x[253], aux, prime);
  DFT2_aux_gmp(x[222], x[254], aux, prime);
  DFT2_aux_gmp(x[223], x[255], aux, prime);
  mult_aux_gmp(x[48], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[49], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[50], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[51], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[52], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[53], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[54], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[56], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[57], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[58], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[60], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[80], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[81], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[82], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[83], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[84], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[85], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[86], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[88], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[89], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[90], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[92], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[112], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[113], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[114], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[116], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[120], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[144], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[145], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[146], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[147], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[148], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[149], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[150], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[151], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[152], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[153], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[154], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[155], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[156], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[157], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[158], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[159], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[176], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[177], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[178], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[179], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[180], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[181], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[182], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[183], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[184], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[185], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[186], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[187], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[188], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[189], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[190], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[191], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[208], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[209], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[210], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[211], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[212], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[213], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[214], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[215], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[216], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[217], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[218], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[219], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[220], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[221], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[222], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[223], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[240], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[241], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[242], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[243], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[244], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[245], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[246], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[247], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[248], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[249], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[250], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[251], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[252], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[253], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[254], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[112], aux,	prime);
  DFT2_aux_gmp(x[0], x[16], aux, prime);
  DFT2_aux_gmp(x[1], x[17], aux, prime);
  DFT2_aux_gmp(x[2], x[18], aux, prime);
  DFT2_aux_gmp(x[3], x[19], aux, prime);
  DFT2_aux_gmp(x[4], x[20], aux, prime);
  DFT2_aux_gmp(x[5], x[21], aux, prime);
  DFT2_aux_gmp(x[6], x[22], aux, prime);
  DFT2_aux_gmp(x[7], x[23], aux, prime);
  DFT2_aux_gmp(x[8], x[24], aux, prime);
  DFT2_aux_gmp(x[9], x[25], aux, prime);
  DFT2_aux_gmp(x[10], x[26], aux, prime);
  DFT2_aux_gmp(x[11], x[27], aux, prime);
  DFT2_aux_gmp(x[12], x[28], aux, prime);
  DFT2_aux_gmp(x[13], x[29], aux, prime);
  DFT2_aux_gmp(x[14], x[30], aux, prime);
  DFT2_aux_gmp(x[15], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[48], aux, prime);
  DFT2_aux_gmp(x[33], x[49], aux, prime);
  DFT2_aux_gmp(x[34], x[50], aux, prime);
  DFT2_aux_gmp(x[35], x[51], aux, prime);
  DFT2_aux_gmp(x[36], x[52], aux, prime);
  DFT2_aux_gmp(x[37], x[53], aux, prime);
  DFT2_aux_gmp(x[38], x[54], aux, prime);
  DFT2_aux_gmp(x[39], x[55], aux, prime);
  DFT2_aux_gmp(x[40], x[56], aux, prime);
  DFT2_aux_gmp(x[41], x[57], aux, prime);
  DFT2_aux_gmp(x[42], x[58], aux, prime);
  DFT2_aux_gmp(x[43], x[59], aux, prime);
  DFT2_aux_gmp(x[44], x[60], aux, prime);
  DFT2_aux_gmp(x[45], x[61], aux, prime);
  DFT2_aux_gmp(x[46], x[62], aux, prime);
  DFT2_aux_gmp(x[47], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[80], aux, prime);
  DFT2_aux_gmp(x[65], x[81], aux, prime);
  DFT2_aux_gmp(x[66], x[82], aux, prime);
  DFT2_aux_gmp(x[67], x[83], aux, prime);
  DFT2_aux_gmp(x[68], x[84], aux, prime);
  DFT2_aux_gmp(x[69], x[85], aux, prime);
  DFT2_aux_gmp(x[70], x[86], aux, prime);
  DFT2_aux_gmp(x[71], x[87], aux, prime);
  DFT2_aux_gmp(x[72], x[88], aux, prime);
  DFT2_aux_gmp(x[73], x[89], aux, prime);
  DFT2_aux_gmp(x[74], x[90], aux, prime);
  DFT2_aux_gmp(x[75], x[91], aux, prime);
  DFT2_aux_gmp(x[76], x[92], aux, prime);
  DFT2_aux_gmp(x[77], x[93], aux, prime);
  DFT2_aux_gmp(x[78], x[94], aux, prime);
  DFT2_aux_gmp(x[79], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[112], aux, prime);
  DFT2_aux_gmp(x[97], x[113], aux, prime);
  DFT2_aux_gmp(x[98], x[114], aux, prime);
  DFT2_aux_gmp(x[99], x[115], aux, prime);
  DFT2_aux_gmp(x[100], x[116], aux, prime);
  DFT2_aux_gmp(x[101], x[117], aux, prime);
  DFT2_aux_gmp(x[102], x[118], aux, prime);
  DFT2_aux_gmp(x[103], x[119], aux, prime);
  DFT2_aux_gmp(x[104], x[120], aux, prime);
  DFT2_aux_gmp(x[105], x[121], aux, prime);
  DFT2_aux_gmp(x[106], x[122], aux, prime);
  DFT2_aux_gmp(x[107], x[123], aux, prime);
  DFT2_aux_gmp(x[108], x[124], aux, prime);
  DFT2_aux_gmp(x[109], x[125], aux, prime);
  DFT2_aux_gmp(x[110], x[126], aux, prime);
  DFT2_aux_gmp(x[111], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[144], aux, prime);
  DFT2_aux_gmp(x[129], x[145], aux, prime);
  DFT2_aux_gmp(x[130], x[146], aux, prime);
  DFT2_aux_gmp(x[131], x[147], aux, prime);
  DFT2_aux_gmp(x[132], x[148], aux, prime);
  DFT2_aux_gmp(x[133], x[149], aux, prime);
  DFT2_aux_gmp(x[134], x[150], aux, prime);
  DFT2_aux_gmp(x[135], x[151], aux, prime);
  DFT2_aux_gmp(x[136], x[152], aux, prime);
  DFT2_aux_gmp(x[137], x[153], aux, prime);
  DFT2_aux_gmp(x[138], x[154], aux, prime);
  DFT2_aux_gmp(x[139], x[155], aux, prime);
  DFT2_aux_gmp(x[140], x[156], aux, prime);
  DFT2_aux_gmp(x[141], x[157], aux, prime);
  DFT2_aux_gmp(x[142], x[158], aux, prime);
  DFT2_aux_gmp(x[143], x[159], aux, prime);
  DFT2_aux_gmp(x[160], x[176], aux, prime);
  DFT2_aux_gmp(x[161], x[177], aux, prime);
  DFT2_aux_gmp(x[162], x[178], aux, prime);
  DFT2_aux_gmp(x[163], x[179], aux, prime);
  DFT2_aux_gmp(x[164], x[180], aux, prime);
  DFT2_aux_gmp(x[165], x[181], aux, prime);
  DFT2_aux_gmp(x[166], x[182], aux, prime);
  DFT2_aux_gmp(x[167], x[183], aux, prime);
  DFT2_aux_gmp(x[168], x[184], aux, prime);
  DFT2_aux_gmp(x[169], x[185], aux, prime);
  DFT2_aux_gmp(x[170], x[186], aux, prime);
  DFT2_aux_gmp(x[171], x[187], aux, prime);
  DFT2_aux_gmp(x[172], x[188], aux, prime);
  DFT2_aux_gmp(x[173], x[189], aux, prime);
  DFT2_aux_gmp(x[174], x[190], aux, prime);
  DFT2_aux_gmp(x[175], x[191], aux, prime);
  DFT2_aux_gmp(x[192], x[208], aux, prime);
  DFT2_aux_gmp(x[193], x[209], aux, prime);
  DFT2_aux_gmp(x[194], x[210], aux, prime);
  DFT2_aux_gmp(x[195], x[211], aux, prime);
  DFT2_aux_gmp(x[196], x[212], aux, prime);
  DFT2_aux_gmp(x[197], x[213], aux, prime);
  DFT2_aux_gmp(x[198], x[214], aux, prime);
  DFT2_aux_gmp(x[199], x[215], aux, prime);
  DFT2_aux_gmp(x[200], x[216], aux, prime);
  DFT2_aux_gmp(x[201], x[217], aux, prime);
  DFT2_aux_gmp(x[202], x[218], aux, prime);
  DFT2_aux_gmp(x[203], x[219], aux, prime);
  DFT2_aux_gmp(x[204], x[220], aux, prime);
  DFT2_aux_gmp(x[205], x[221], aux, prime);
  DFT2_aux_gmp(x[206], x[222], aux, prime);
  DFT2_aux_gmp(x[207], x[223], aux, prime);
  DFT2_aux_gmp(x[224], x[240], aux, prime);
  DFT2_aux_gmp(x[225], x[241], aux, prime);
  DFT2_aux_gmp(x[226], x[242], aux, prime);
  DFT2_aux_gmp(x[227], x[243], aux, prime);
  DFT2_aux_gmp(x[228], x[244], aux, prime);
  DFT2_aux_gmp(x[229], x[245], aux, prime);
  DFT2_aux_gmp(x[230], x[246], aux, prime);
  DFT2_aux_gmp(x[231], x[247], aux, prime);
  DFT2_aux_gmp(x[232], x[248], aux, prime);
  DFT2_aux_gmp(x[233], x[249], aux, prime);
  DFT2_aux_gmp(x[234], x[250], aux, prime);
  DFT2_aux_gmp(x[235], x[251], aux, prime);
  DFT2_aux_gmp(x[236], x[252], aux, prime);
  DFT2_aux_gmp(x[237], x[253], aux, prime);
  DFT2_aux_gmp(x[238], x[254], aux, prime);
  DFT2_aux_gmp(x[239], x[255], aux, prime);
  mult_aux_gmp(x[24], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[25], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[26], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[27], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[28], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[29], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[30], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[40], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[41], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[42], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[43], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[44], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[45], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[46], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[56], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[57], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[58], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[60], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[72], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[73], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[74], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[75], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[76], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[77], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[78], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[88], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[89], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[90], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[92], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[104], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[105], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[106], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[108], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[120], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[136], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[137], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[138], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[139], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[140], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[141], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[142], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[143], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[152], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[153], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[154], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[155], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[156], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[157], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[158], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[159], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[168], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[169], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[170], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[171], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[172], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[173], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[174], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[175], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[184], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[185], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[186], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[187], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[188], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[189], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[190], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[191], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[200], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[201], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[202], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[203], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[204], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[205], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[206], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[207], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[216], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[217], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[218], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[219], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[220], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[221], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[222], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[223], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[232], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[233], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[234], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[235], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[236], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[237], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[238], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[239], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[248], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[249], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[250], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[251], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[252], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[253], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[254], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[120], aux,	prime);
  DFT2_aux_gmp(x[0], x[8], aux, prime);
  DFT2_aux_gmp(x[1], x[9], aux, prime);
  DFT2_aux_gmp(x[2], x[10], aux, prime);
  DFT2_aux_gmp(x[3], x[11], aux, prime);
  DFT2_aux_gmp(x[4], x[12], aux, prime);
  DFT2_aux_gmp(x[5], x[13], aux, prime);
  DFT2_aux_gmp(x[6], x[14], aux, prime);
  DFT2_aux_gmp(x[7], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[24], aux, prime);
  DFT2_aux_gmp(x[17], x[25], aux, prime);
  DFT2_aux_gmp(x[18], x[26], aux, prime);
  DFT2_aux_gmp(x[19], x[27], aux, prime);
  DFT2_aux_gmp(x[20], x[28], aux, prime);
  DFT2_aux_gmp(x[21], x[29], aux, prime);
  DFT2_aux_gmp(x[22], x[30], aux, prime);
  DFT2_aux_gmp(x[23], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[40], aux, prime);
  DFT2_aux_gmp(x[33], x[41], aux, prime);
  DFT2_aux_gmp(x[34], x[42], aux, prime);
  DFT2_aux_gmp(x[35], x[43], aux, prime);
  DFT2_aux_gmp(x[36], x[44], aux, prime);
  DFT2_aux_gmp(x[37], x[45], aux, prime);
  DFT2_aux_gmp(x[38], x[46], aux, prime);
  DFT2_aux_gmp(x[39], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[56], aux, prime);
  DFT2_aux_gmp(x[49], x[57], aux, prime);
  DFT2_aux_gmp(x[50], x[58], aux, prime);
  DFT2_aux_gmp(x[51], x[59], aux, prime);
  DFT2_aux_gmp(x[52], x[60], aux, prime);
  DFT2_aux_gmp(x[53], x[61], aux, prime);
  DFT2_aux_gmp(x[54], x[62], aux, prime);
  DFT2_aux_gmp(x[55], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[72], aux, prime);
  DFT2_aux_gmp(x[65], x[73], aux, prime);
  DFT2_aux_gmp(x[66], x[74], aux, prime);
  DFT2_aux_gmp(x[67], x[75], aux, prime);
  DFT2_aux_gmp(x[68], x[76], aux, prime);
  DFT2_aux_gmp(x[69], x[77], aux, prime);
  DFT2_aux_gmp(x[70], x[78], aux, prime);
  DFT2_aux_gmp(x[71], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[88], aux, prime);
  DFT2_aux_gmp(x[81], x[89], aux, prime);
  DFT2_aux_gmp(x[82], x[90], aux, prime);
  DFT2_aux_gmp(x[83], x[91], aux, prime);
  DFT2_aux_gmp(x[84], x[92], aux, prime);
  DFT2_aux_gmp(x[85], x[93], aux, prime);
  DFT2_aux_gmp(x[86], x[94], aux, prime);
  DFT2_aux_gmp(x[87], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[104], aux, prime);
  DFT2_aux_gmp(x[97], x[105], aux, prime);
  DFT2_aux_gmp(x[98], x[106], aux, prime);
  DFT2_aux_gmp(x[99], x[107], aux, prime);
  DFT2_aux_gmp(x[100], x[108], aux, prime);
  DFT2_aux_gmp(x[101], x[109], aux, prime);
  DFT2_aux_gmp(x[102], x[110], aux, prime);
  DFT2_aux_gmp(x[103], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[120], aux, prime);
  DFT2_aux_gmp(x[113], x[121], aux, prime);
  DFT2_aux_gmp(x[114], x[122], aux, prime);
  DFT2_aux_gmp(x[115], x[123], aux, prime);
  DFT2_aux_gmp(x[116], x[124], aux, prime);
  DFT2_aux_gmp(x[117], x[125], aux, prime);
  DFT2_aux_gmp(x[118], x[126], aux, prime);
  DFT2_aux_gmp(x[119], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[136], aux, prime);
  DFT2_aux_gmp(x[129], x[137], aux, prime);
  DFT2_aux_gmp(x[130], x[138], aux, prime);
  DFT2_aux_gmp(x[131], x[139], aux, prime);
  DFT2_aux_gmp(x[132], x[140], aux, prime);
  DFT2_aux_gmp(x[133], x[141], aux, prime);
  DFT2_aux_gmp(x[134], x[142], aux, prime);
  DFT2_aux_gmp(x[135], x[143], aux, prime);
  DFT2_aux_gmp(x[144], x[152], aux, prime);
  DFT2_aux_gmp(x[145], x[153], aux, prime);
  DFT2_aux_gmp(x[146], x[154], aux, prime);
  DFT2_aux_gmp(x[147], x[155], aux, prime);
  DFT2_aux_gmp(x[148], x[156], aux, prime);
  DFT2_aux_gmp(x[149], x[157], aux, prime);
  DFT2_aux_gmp(x[150], x[158], aux, prime);
  DFT2_aux_gmp(x[151], x[159], aux, prime);
  DFT2_aux_gmp(x[160], x[168], aux, prime);
  DFT2_aux_gmp(x[161], x[169], aux, prime);
  DFT2_aux_gmp(x[162], x[170], aux, prime);
  DFT2_aux_gmp(x[163], x[171], aux, prime);
  DFT2_aux_gmp(x[164], x[172], aux, prime);
  DFT2_aux_gmp(x[165], x[173], aux, prime);
  DFT2_aux_gmp(x[166], x[174], aux, prime);
  DFT2_aux_gmp(x[167], x[175], aux, prime);
  DFT2_aux_gmp(x[176], x[184], aux, prime);
  DFT2_aux_gmp(x[177], x[185], aux, prime);
  DFT2_aux_gmp(x[178], x[186], aux, prime);
  DFT2_aux_gmp(x[179], x[187], aux, prime);
  DFT2_aux_gmp(x[180], x[188], aux, prime);
  DFT2_aux_gmp(x[181], x[189], aux, prime);
  DFT2_aux_gmp(x[182], x[190], aux, prime);
  DFT2_aux_gmp(x[183], x[191], aux, prime);
  DFT2_aux_gmp(x[192], x[200], aux, prime);
  DFT2_aux_gmp(x[193], x[201], aux, prime);
  DFT2_aux_gmp(x[194], x[202], aux, prime);
  DFT2_aux_gmp(x[195], x[203], aux, prime);
  DFT2_aux_gmp(x[196], x[204], aux, prime);
  DFT2_aux_gmp(x[197], x[205], aux, prime);
  DFT2_aux_gmp(x[198], x[206], aux, prime);
  DFT2_aux_gmp(x[199], x[207], aux, prime);
  DFT2_aux_gmp(x[208], x[216], aux, prime);
  DFT2_aux_gmp(x[209], x[217], aux, prime);
  DFT2_aux_gmp(x[210], x[218], aux, prime);
  DFT2_aux_gmp(x[211], x[219], aux, prime);
  DFT2_aux_gmp(x[212], x[220], aux, prime);
  DFT2_aux_gmp(x[213], x[221], aux, prime);
  DFT2_aux_gmp(x[214], x[222], aux, prime);
  DFT2_aux_gmp(x[215], x[223], aux, prime);
  DFT2_aux_gmp(x[224], x[232], aux, prime);
  DFT2_aux_gmp(x[225], x[233], aux, prime);
  DFT2_aux_gmp(x[226], x[234], aux, prime);
  DFT2_aux_gmp(x[227], x[235], aux, prime);
  DFT2_aux_gmp(x[228], x[236], aux, prime);
  DFT2_aux_gmp(x[229], x[237], aux, prime);
  DFT2_aux_gmp(x[230], x[238], aux, prime);
  DFT2_aux_gmp(x[231], x[239], aux, prime);
  DFT2_aux_gmp(x[240], x[248], aux, prime);
  DFT2_aux_gmp(x[241], x[249], aux, prime);
  DFT2_aux_gmp(x[242], x[250], aux, prime);
  DFT2_aux_gmp(x[243], x[251], aux, prime);
  DFT2_aux_gmp(x[244], x[252], aux, prime);
  DFT2_aux_gmp(x[245], x[253], aux, prime);
  DFT2_aux_gmp(x[246], x[254], aux, prime);
  DFT2_aux_gmp(x[247], x[255], aux, prime);
  mult_aux_gmp(x[12], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[13], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[14], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[15], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[20], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[21], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[22], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[23], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[28], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[29], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[30], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[36], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[37], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[38], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[39], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[44], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[45], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[46], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[52], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[53], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[54], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[60], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[68], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[69], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[70], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[71], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[76], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[77], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[78], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[84], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[85], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[86], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[92], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[100], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[101], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[102], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[108], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[116], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[124], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[132], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[133], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[134], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[135], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[140], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[141], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[142], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[143], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[148], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[149], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[150], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[151], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[156], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[157], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[158], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[159], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[164], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[165], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[166], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[167], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[172], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[173], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[174], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[175], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[180], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[181], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[182], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[183], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[188], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[189], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[190], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[191], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[196], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[197], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[198], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[199], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[204], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[205], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[206], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[207], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[212], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[213], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[214], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[215], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[220], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[221], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[222], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[223], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[228], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[229], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[230], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[231], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[236], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[237], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[238], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[239], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[244], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[245], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[246], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[247], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[252], precompute_pow_omega[124], aux,	prime);
  mult_aux_gmp(x[253], precompute_pow_omega[124], aux,	prime);
  mult_aux_gmp(x[254], precompute_pow_omega[124], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[124], aux,	prime);
  DFT2_aux_gmp(x[0], x[4], aux, prime);
  DFT2_aux_gmp(x[1], x[5], aux, prime);
  DFT2_aux_gmp(x[2], x[6], aux, prime);
  DFT2_aux_gmp(x[3], x[7], aux, prime);
  DFT2_aux_gmp(x[8], x[12], aux, prime);
  DFT2_aux_gmp(x[9], x[13], aux, prime);
  DFT2_aux_gmp(x[10], x[14], aux, prime);
  DFT2_aux_gmp(x[11], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[20], aux, prime);
  DFT2_aux_gmp(x[17], x[21], aux, prime);
  DFT2_aux_gmp(x[18], x[22], aux, prime);
  DFT2_aux_gmp(x[19], x[23], aux, prime);
  DFT2_aux_gmp(x[24], x[28], aux, prime);
  DFT2_aux_gmp(x[25], x[29], aux, prime);
  DFT2_aux_gmp(x[26], x[30], aux, prime);
  DFT2_aux_gmp(x[27], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[36], aux, prime);
  DFT2_aux_gmp(x[33], x[37], aux, prime);
  DFT2_aux_gmp(x[34], x[38], aux, prime);
  DFT2_aux_gmp(x[35], x[39], aux, prime);
  DFT2_aux_gmp(x[40], x[44], aux, prime);
  DFT2_aux_gmp(x[41], x[45], aux, prime);
  DFT2_aux_gmp(x[42], x[46], aux, prime);
  DFT2_aux_gmp(x[43], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[52], aux, prime);
  DFT2_aux_gmp(x[49], x[53], aux, prime);
  DFT2_aux_gmp(x[50], x[54], aux, prime);
  DFT2_aux_gmp(x[51], x[55], aux, prime);
  DFT2_aux_gmp(x[56], x[60], aux, prime);
  DFT2_aux_gmp(x[57], x[61], aux, prime);
  DFT2_aux_gmp(x[58], x[62], aux, prime);
  DFT2_aux_gmp(x[59], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[68], aux, prime);
  DFT2_aux_gmp(x[65], x[69], aux, prime);
  DFT2_aux_gmp(x[66], x[70], aux, prime);
  DFT2_aux_gmp(x[67], x[71], aux, prime);
  DFT2_aux_gmp(x[72], x[76], aux, prime);
  DFT2_aux_gmp(x[73], x[77], aux, prime);
  DFT2_aux_gmp(x[74], x[78], aux, prime);
  DFT2_aux_gmp(x[75], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[84], aux, prime);
  DFT2_aux_gmp(x[81], x[85], aux, prime);
  DFT2_aux_gmp(x[82], x[86], aux, prime);
  DFT2_aux_gmp(x[83], x[87], aux, prime);
  DFT2_aux_gmp(x[88], x[92], aux, prime);
  DFT2_aux_gmp(x[89], x[93], aux, prime);
  DFT2_aux_gmp(x[90], x[94], aux, prime);
  DFT2_aux_gmp(x[91], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[100], aux, prime);
  DFT2_aux_gmp(x[97], x[101], aux, prime);
  DFT2_aux_gmp(x[98], x[102], aux, prime);
  DFT2_aux_gmp(x[99], x[103], aux, prime);
  DFT2_aux_gmp(x[104], x[108], aux, prime);
  DFT2_aux_gmp(x[105], x[109], aux, prime);
  DFT2_aux_gmp(x[106], x[110], aux, prime);
  DFT2_aux_gmp(x[107], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[116], aux, prime);
  DFT2_aux_gmp(x[113], x[117], aux, prime);
  DFT2_aux_gmp(x[114], x[118], aux, prime);
  DFT2_aux_gmp(x[115], x[119], aux, prime);
  DFT2_aux_gmp(x[120], x[124], aux, prime);
  DFT2_aux_gmp(x[121], x[125], aux, prime);
  DFT2_aux_gmp(x[122], x[126], aux, prime);
  DFT2_aux_gmp(x[123], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[132], aux, prime);
  DFT2_aux_gmp(x[129], x[133], aux, prime);
  DFT2_aux_gmp(x[130], x[134], aux, prime);
  DFT2_aux_gmp(x[131], x[135], aux, prime);
  DFT2_aux_gmp(x[136], x[140], aux, prime);
  DFT2_aux_gmp(x[137], x[141], aux, prime);
  DFT2_aux_gmp(x[138], x[142], aux, prime);
  DFT2_aux_gmp(x[139], x[143], aux, prime);
  DFT2_aux_gmp(x[144], x[148], aux, prime);
  DFT2_aux_gmp(x[145], x[149], aux, prime);
  DFT2_aux_gmp(x[146], x[150], aux, prime);
  DFT2_aux_gmp(x[147], x[151], aux, prime);
  DFT2_aux_gmp(x[152], x[156], aux, prime);
  DFT2_aux_gmp(x[153], x[157], aux, prime);
  DFT2_aux_gmp(x[154], x[158], aux, prime);
  DFT2_aux_gmp(x[155], x[159], aux, prime);
  DFT2_aux_gmp(x[160], x[164], aux, prime);
  DFT2_aux_gmp(x[161], x[165], aux, prime);
  DFT2_aux_gmp(x[162], x[166], aux, prime);
  DFT2_aux_gmp(x[163], x[167], aux, prime);
  DFT2_aux_gmp(x[168], x[172], aux, prime);
  DFT2_aux_gmp(x[169], x[173], aux, prime);
  DFT2_aux_gmp(x[170], x[174], aux, prime);
  DFT2_aux_gmp(x[171], x[175], aux, prime);
  DFT2_aux_gmp(x[176], x[180], aux, prime);
  DFT2_aux_gmp(x[177], x[181], aux, prime);
  DFT2_aux_gmp(x[178], x[182], aux, prime);
  DFT2_aux_gmp(x[179], x[183], aux, prime);
  DFT2_aux_gmp(x[184], x[188], aux, prime);
  DFT2_aux_gmp(x[185], x[189], aux, prime);
  DFT2_aux_gmp(x[186], x[190], aux, prime);
  DFT2_aux_gmp(x[187], x[191], aux, prime);
  DFT2_aux_gmp(x[192], x[196], aux, prime);
  DFT2_aux_gmp(x[193], x[197], aux, prime);
  DFT2_aux_gmp(x[194], x[198], aux, prime);
  DFT2_aux_gmp(x[195], x[199], aux, prime);
  DFT2_aux_gmp(x[200], x[204], aux, prime);
  DFT2_aux_gmp(x[201], x[205], aux, prime);
  DFT2_aux_gmp(x[202], x[206], aux, prime);
  DFT2_aux_gmp(x[203], x[207], aux, prime);
  DFT2_aux_gmp(x[208], x[212], aux, prime);
  DFT2_aux_gmp(x[209], x[213], aux, prime);
  DFT2_aux_gmp(x[210], x[214], aux, prime);
  DFT2_aux_gmp(x[211], x[215], aux, prime);
  DFT2_aux_gmp(x[216], x[220], aux, prime);
  DFT2_aux_gmp(x[217], x[221], aux, prime);
  DFT2_aux_gmp(x[218], x[222], aux, prime);
  DFT2_aux_gmp(x[219], x[223], aux, prime);
  DFT2_aux_gmp(x[224], x[228], aux, prime);
  DFT2_aux_gmp(x[225], x[229], aux, prime);
  DFT2_aux_gmp(x[226], x[230], aux, prime);
  DFT2_aux_gmp(x[227], x[231], aux, prime);
  DFT2_aux_gmp(x[232], x[236], aux, prime);
  DFT2_aux_gmp(x[233], x[237], aux, prime);
  DFT2_aux_gmp(x[234], x[238], aux, prime);
  DFT2_aux_gmp(x[235], x[239], aux, prime);
  DFT2_aux_gmp(x[240], x[244], aux, prime);
  DFT2_aux_gmp(x[241], x[245], aux, prime);
  DFT2_aux_gmp(x[242], x[246], aux, prime);
  DFT2_aux_gmp(x[243], x[247], aux, prime);
  DFT2_aux_gmp(x[248], x[252], aux, prime);
  DFT2_aux_gmp(x[249], x[253], aux, prime);
  DFT2_aux_gmp(x[250], x[254], aux, prime);
  DFT2_aux_gmp(x[251], x[255], aux, prime);
  mult_aux_gmp(x[6], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[7], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[10], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[11], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[14], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[15], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[18], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[19], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[22], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[23], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[26], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[27], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[30], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[34], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[35], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[38], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[39], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[42], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[43], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[46], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[50], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[51], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[54], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[58], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[62], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[66], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[67], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[70], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[71], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[74], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[75], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[78], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[82], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[83], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[86], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[90], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[94], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[98], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[99], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[102], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[106], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[110], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[114], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[118], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[122], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[126], precompute_pow_omega[124], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[124], aux,	prime);
  mult_aux_gmp(x[130], precompute_pow_omega[2], aux,	prime);
  mult_aux_gmp(x[131], precompute_pow_omega[2], aux,	prime);
  mult_aux_gmp(x[134], precompute_pow_omega[66], aux,	prime);
  mult_aux_gmp(x[135], precompute_pow_omega[66], aux,	prime);
  mult_aux_gmp(x[138], precompute_pow_omega[34], aux,	prime);
  mult_aux_gmp(x[139], precompute_pow_omega[34], aux,	prime);
  mult_aux_gmp(x[142], precompute_pow_omega[98], aux,	prime);
  mult_aux_gmp(x[143], precompute_pow_omega[98], aux,	prime);
  mult_aux_gmp(x[146], precompute_pow_omega[18], aux,	prime);
  mult_aux_gmp(x[147], precompute_pow_omega[18], aux,	prime);
  mult_aux_gmp(x[150], precompute_pow_omega[82], aux,	prime);
  mult_aux_gmp(x[151], precompute_pow_omega[82], aux,	prime);
  mult_aux_gmp(x[154], precompute_pow_omega[50], aux,	prime);
  mult_aux_gmp(x[155], precompute_pow_omega[50], aux,	prime);
  mult_aux_gmp(x[158], precompute_pow_omega[114], aux,	prime);
  mult_aux_gmp(x[159], precompute_pow_omega[114], aux,	prime);
  mult_aux_gmp(x[162], precompute_pow_omega[10], aux,	prime);
  mult_aux_gmp(x[163], precompute_pow_omega[10], aux,	prime);
  mult_aux_gmp(x[166], precompute_pow_omega[74], aux,	prime);
  mult_aux_gmp(x[167], precompute_pow_omega[74], aux,	prime);
  mult_aux_gmp(x[170], precompute_pow_omega[42], aux,	prime);
  mult_aux_gmp(x[171], precompute_pow_omega[42], aux,	prime);
  mult_aux_gmp(x[174], precompute_pow_omega[106], aux,	prime);
  mult_aux_gmp(x[175], precompute_pow_omega[106], aux,	prime);
  mult_aux_gmp(x[178], precompute_pow_omega[26], aux,	prime);
  mult_aux_gmp(x[179], precompute_pow_omega[26], aux,	prime);
  mult_aux_gmp(x[182], precompute_pow_omega[90], aux,	prime);
  mult_aux_gmp(x[183], precompute_pow_omega[90], aux,	prime);
  mult_aux_gmp(x[186], precompute_pow_omega[58], aux,	prime);
  mult_aux_gmp(x[187], precompute_pow_omega[58], aux,	prime);
  mult_aux_gmp(x[190], precompute_pow_omega[122], aux,	prime);
  mult_aux_gmp(x[191], precompute_pow_omega[122], aux,	prime);
  mult_aux_gmp(x[194], precompute_pow_omega[6], aux,	prime);
  mult_aux_gmp(x[195], precompute_pow_omega[6], aux,	prime);
  mult_aux_gmp(x[198], precompute_pow_omega[70], aux,	prime);
  mult_aux_gmp(x[199], precompute_pow_omega[70], aux,	prime);
  mult_aux_gmp(x[202], precompute_pow_omega[38], aux,	prime);
  mult_aux_gmp(x[203], precompute_pow_omega[38], aux,	prime);
  mult_aux_gmp(x[206], precompute_pow_omega[102], aux,	prime);
  mult_aux_gmp(x[207], precompute_pow_omega[102], aux,	prime);
  mult_aux_gmp(x[210], precompute_pow_omega[22], aux,	prime);
  mult_aux_gmp(x[211], precompute_pow_omega[22], aux,	prime);
  mult_aux_gmp(x[214], precompute_pow_omega[86], aux,	prime);
  mult_aux_gmp(x[215], precompute_pow_omega[86], aux,	prime);
  mult_aux_gmp(x[218], precompute_pow_omega[54], aux,	prime);
  mult_aux_gmp(x[219], precompute_pow_omega[54], aux,	prime);
  mult_aux_gmp(x[222], precompute_pow_omega[118], aux,	prime);
  mult_aux_gmp(x[223], precompute_pow_omega[118], aux,	prime);
  mult_aux_gmp(x[226], precompute_pow_omega[14], aux,	prime);
  mult_aux_gmp(x[227], precompute_pow_omega[14], aux,	prime);
  mult_aux_gmp(x[230], precompute_pow_omega[78], aux,	prime);
  mult_aux_gmp(x[231], precompute_pow_omega[78], aux,	prime);
  mult_aux_gmp(x[234], precompute_pow_omega[46], aux,	prime);
  mult_aux_gmp(x[235], precompute_pow_omega[46], aux,	prime);
  mult_aux_gmp(x[238], precompute_pow_omega[110], aux,	prime);
  mult_aux_gmp(x[239], precompute_pow_omega[110], aux,	prime);
  mult_aux_gmp(x[242], precompute_pow_omega[30], aux,	prime);
  mult_aux_gmp(x[243], precompute_pow_omega[30], aux,	prime);
  mult_aux_gmp(x[246], precompute_pow_omega[94], aux,	prime);
  mult_aux_gmp(x[247], precompute_pow_omega[94], aux,	prime);
  mult_aux_gmp(x[250], precompute_pow_omega[62], aux,	prime);
  mult_aux_gmp(x[251], precompute_pow_omega[62], aux,	prime);
  mult_aux_gmp(x[254], precompute_pow_omega[126], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[126], aux,	prime);
  DFT2_aux_gmp(x[0], x[2], aux, prime);
  DFT2_aux_gmp(x[1], x[3], aux, prime);
  DFT2_aux_gmp(x[4], x[6], aux, prime);
  DFT2_aux_gmp(x[5], x[7], aux, prime);
  DFT2_aux_gmp(x[8], x[10], aux, prime);
  DFT2_aux_gmp(x[9], x[11], aux, prime);
  DFT2_aux_gmp(x[12], x[14], aux, prime);
  DFT2_aux_gmp(x[13], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[18], aux, prime);
  DFT2_aux_gmp(x[17], x[19], aux, prime);
  DFT2_aux_gmp(x[20], x[22], aux, prime);
  DFT2_aux_gmp(x[21], x[23], aux, prime);
  DFT2_aux_gmp(x[24], x[26], aux, prime);
  DFT2_aux_gmp(x[25], x[27], aux, prime);
  DFT2_aux_gmp(x[28], x[30], aux, prime);
  DFT2_aux_gmp(x[29], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[34], aux, prime);
  DFT2_aux_gmp(x[33], x[35], aux, prime);
  DFT2_aux_gmp(x[36], x[38], aux, prime);
  DFT2_aux_gmp(x[37], x[39], aux, prime);
  DFT2_aux_gmp(x[40], x[42], aux, prime);
  DFT2_aux_gmp(x[41], x[43], aux, prime);
  DFT2_aux_gmp(x[44], x[46], aux, prime);
  DFT2_aux_gmp(x[45], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[50], aux, prime);
  DFT2_aux_gmp(x[49], x[51], aux, prime);
  DFT2_aux_gmp(x[52], x[54], aux, prime);
  DFT2_aux_gmp(x[53], x[55], aux, prime);
  DFT2_aux_gmp(x[56], x[58], aux, prime);
  DFT2_aux_gmp(x[57], x[59], aux, prime);
  DFT2_aux_gmp(x[60], x[62], aux, prime);
  DFT2_aux_gmp(x[61], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[66], aux, prime);
  DFT2_aux_gmp(x[65], x[67], aux, prime);
  DFT2_aux_gmp(x[68], x[70], aux, prime);
  DFT2_aux_gmp(x[69], x[71], aux, prime);
  DFT2_aux_gmp(x[72], x[74], aux, prime);
  DFT2_aux_gmp(x[73], x[75], aux, prime);
  DFT2_aux_gmp(x[76], x[78], aux, prime);
  DFT2_aux_gmp(x[77], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[82], aux, prime);
  DFT2_aux_gmp(x[81], x[83], aux, prime);
  DFT2_aux_gmp(x[84], x[86], aux, prime);
  DFT2_aux_gmp(x[85], x[87], aux, prime);
  DFT2_aux_gmp(x[88], x[90], aux, prime);
  DFT2_aux_gmp(x[89], x[91], aux, prime);
  DFT2_aux_gmp(x[92], x[94], aux, prime);
  DFT2_aux_gmp(x[93], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[98], aux, prime);
  DFT2_aux_gmp(x[97], x[99], aux, prime);
  DFT2_aux_gmp(x[100], x[102], aux, prime);
  DFT2_aux_gmp(x[101], x[103], aux, prime);
  DFT2_aux_gmp(x[104], x[106], aux, prime);
  DFT2_aux_gmp(x[105], x[107], aux, prime);
  DFT2_aux_gmp(x[108], x[110], aux, prime);
  DFT2_aux_gmp(x[109], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[114], aux, prime);
  DFT2_aux_gmp(x[113], x[115], aux, prime);
  DFT2_aux_gmp(x[116], x[118], aux, prime);
  DFT2_aux_gmp(x[117], x[119], aux, prime);
  DFT2_aux_gmp(x[120], x[122], aux, prime);
  DFT2_aux_gmp(x[121], x[123], aux, prime);
  DFT2_aux_gmp(x[124], x[126], aux, prime);
  DFT2_aux_gmp(x[125], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[130], aux, prime);
  DFT2_aux_gmp(x[129], x[131], aux, prime);
  DFT2_aux_gmp(x[132], x[134], aux, prime);
  DFT2_aux_gmp(x[133], x[135], aux, prime);
  DFT2_aux_gmp(x[136], x[138], aux, prime);
  DFT2_aux_gmp(x[137], x[139], aux, prime);
  DFT2_aux_gmp(x[140], x[142], aux, prime);
  DFT2_aux_gmp(x[141], x[143], aux, prime);
  DFT2_aux_gmp(x[144], x[146], aux, prime);
  DFT2_aux_gmp(x[145], x[147], aux, prime);
  DFT2_aux_gmp(x[148], x[150], aux, prime);
  DFT2_aux_gmp(x[149], x[151], aux, prime);
  DFT2_aux_gmp(x[152], x[154], aux, prime);
  DFT2_aux_gmp(x[153], x[155], aux, prime);
  DFT2_aux_gmp(x[156], x[158], aux, prime);
  DFT2_aux_gmp(x[157], x[159], aux, prime);
  DFT2_aux_gmp(x[160], x[162], aux, prime);
  DFT2_aux_gmp(x[161], x[163], aux, prime);
  DFT2_aux_gmp(x[164], x[166], aux, prime);
  DFT2_aux_gmp(x[165], x[167], aux, prime);
  DFT2_aux_gmp(x[168], x[170], aux, prime);
  DFT2_aux_gmp(x[169], x[171], aux, prime);
  DFT2_aux_gmp(x[172], x[174], aux, prime);
  DFT2_aux_gmp(x[173], x[175], aux, prime);
  DFT2_aux_gmp(x[176], x[178], aux, prime);
  DFT2_aux_gmp(x[177], x[179], aux, prime);
  DFT2_aux_gmp(x[180], x[182], aux, prime);
  DFT2_aux_gmp(x[181], x[183], aux, prime);
  DFT2_aux_gmp(x[184], x[186], aux, prime);
  DFT2_aux_gmp(x[185], x[187], aux, prime);
  DFT2_aux_gmp(x[188], x[190], aux, prime);
  DFT2_aux_gmp(x[189], x[191], aux, prime);
  DFT2_aux_gmp(x[192], x[194], aux, prime);
  DFT2_aux_gmp(x[193], x[195], aux, prime);
  DFT2_aux_gmp(x[196], x[198], aux, prime);
  DFT2_aux_gmp(x[197], x[199], aux, prime);
  DFT2_aux_gmp(x[200], x[202], aux, prime);
  DFT2_aux_gmp(x[201], x[203], aux, prime);
  DFT2_aux_gmp(x[204], x[206], aux, prime);
  DFT2_aux_gmp(x[205], x[207], aux, prime);
  DFT2_aux_gmp(x[208], x[210], aux, prime);
  DFT2_aux_gmp(x[209], x[211], aux, prime);
  DFT2_aux_gmp(x[212], x[214], aux, prime);
  DFT2_aux_gmp(x[213], x[215], aux, prime);
  DFT2_aux_gmp(x[216], x[218], aux, prime);
  DFT2_aux_gmp(x[217], x[219], aux, prime);
  DFT2_aux_gmp(x[220], x[222], aux, prime);
  DFT2_aux_gmp(x[221], x[223], aux, prime);
  DFT2_aux_gmp(x[224], x[226], aux, prime);
  DFT2_aux_gmp(x[225], x[227], aux, prime);
  DFT2_aux_gmp(x[228], x[230], aux, prime);
  DFT2_aux_gmp(x[229], x[231], aux, prime);
  DFT2_aux_gmp(x[232], x[234], aux, prime);
  DFT2_aux_gmp(x[233], x[235], aux, prime);
  DFT2_aux_gmp(x[236], x[238], aux, prime);
  DFT2_aux_gmp(x[237], x[239], aux, prime);
  DFT2_aux_gmp(x[240], x[242], aux, prime);
  DFT2_aux_gmp(x[241], x[243], aux, prime);
  DFT2_aux_gmp(x[244], x[246], aux, prime);
  DFT2_aux_gmp(x[245], x[247], aux, prime);
  DFT2_aux_gmp(x[248], x[250], aux, prime);
  DFT2_aux_gmp(x[249], x[251], aux, prime);
  DFT2_aux_gmp(x[252], x[254], aux, prime);
  DFT2_aux_gmp(x[253], x[255], aux, prime);
  mult_aux_gmp(x[3], precompute_pow_omega[64], aux,	prime);
  mult_aux_gmp(x[5], precompute_pow_omega[32], aux,	prime);
  mult_aux_gmp(x[7], precompute_pow_omega[96], aux,	prime);
  mult_aux_gmp(x[9], precompute_pow_omega[16], aux,	prime);
  mult_aux_gmp(x[11], precompute_pow_omega[80], aux,	prime);
  mult_aux_gmp(x[13], precompute_pow_omega[48], aux,	prime);
  mult_aux_gmp(x[15], precompute_pow_omega[112], aux,	prime);
  mult_aux_gmp(x[17], precompute_pow_omega[8], aux,	prime);
  mult_aux_gmp(x[19], precompute_pow_omega[72], aux,	prime);
  mult_aux_gmp(x[21], precompute_pow_omega[40], aux,	prime);
  mult_aux_gmp(x[23], precompute_pow_omega[104], aux,	prime);
  mult_aux_gmp(x[25], precompute_pow_omega[24], aux,	prime);
  mult_aux_gmp(x[27], precompute_pow_omega[88], aux,	prime);
  mult_aux_gmp(x[29], precompute_pow_omega[56], aux,	prime);
  mult_aux_gmp(x[31], precompute_pow_omega[120], aux,	prime);
  mult_aux_gmp(x[33], precompute_pow_omega[4], aux,	prime);
  mult_aux_gmp(x[35], precompute_pow_omega[68], aux,	prime);
  mult_aux_gmp(x[37], precompute_pow_omega[36], aux,	prime);
  mult_aux_gmp(x[39], precompute_pow_omega[100], aux,	prime);
  mult_aux_gmp(x[41], precompute_pow_omega[20], aux,	prime);
  mult_aux_gmp(x[43], precompute_pow_omega[84], aux,	prime);
  mult_aux_gmp(x[45], precompute_pow_omega[52], aux,	prime);
  mult_aux_gmp(x[47], precompute_pow_omega[116], aux,	prime);
  mult_aux_gmp(x[49], precompute_pow_omega[12], aux,	prime);
  mult_aux_gmp(x[51], precompute_pow_omega[76], aux,	prime);
  mult_aux_gmp(x[53], precompute_pow_omega[44], aux,	prime);
  mult_aux_gmp(x[55], precompute_pow_omega[108], aux,	prime);
  mult_aux_gmp(x[57], precompute_pow_omega[28], aux,	prime);
  mult_aux_gmp(x[59], precompute_pow_omega[92], aux,	prime);
  mult_aux_gmp(x[61], precompute_pow_omega[60], aux,	prime);
  mult_aux_gmp(x[63], precompute_pow_omega[124], aux,	prime);
  mult_aux_gmp(x[65], precompute_pow_omega[2], aux,	prime);
  mult_aux_gmp(x[67], precompute_pow_omega[66], aux,	prime);
  mult_aux_gmp(x[69], precompute_pow_omega[34], aux,	prime);
  mult_aux_gmp(x[71], precompute_pow_omega[98], aux,	prime);
  mult_aux_gmp(x[73], precompute_pow_omega[18], aux,	prime);
  mult_aux_gmp(x[75], precompute_pow_omega[82], aux,	prime);
  mult_aux_gmp(x[77], precompute_pow_omega[50], aux,	prime);
  mult_aux_gmp(x[79], precompute_pow_omega[114], aux,	prime);
  mult_aux_gmp(x[81], precompute_pow_omega[10], aux,	prime);
  mult_aux_gmp(x[83], precompute_pow_omega[74], aux,	prime);
  mult_aux_gmp(x[85], precompute_pow_omega[42], aux,	prime);
  mult_aux_gmp(x[87], precompute_pow_omega[106], aux,	prime);
  mult_aux_gmp(x[89], precompute_pow_omega[26], aux,	prime);
  mult_aux_gmp(x[91], precompute_pow_omega[90], aux,	prime);
  mult_aux_gmp(x[93], precompute_pow_omega[58], aux,	prime);
  mult_aux_gmp(x[95], precompute_pow_omega[122], aux,	prime);
  mult_aux_gmp(x[97], precompute_pow_omega[6], aux,	prime);
  mult_aux_gmp(x[99], precompute_pow_omega[70], aux,	prime);
  mult_aux_gmp(x[101], precompute_pow_omega[38], aux,	prime);
  mult_aux_gmp(x[103], precompute_pow_omega[102], aux,	prime);
  mult_aux_gmp(x[105], precompute_pow_omega[22], aux,	prime);
  mult_aux_gmp(x[107], precompute_pow_omega[86], aux,	prime);
  mult_aux_gmp(x[109], precompute_pow_omega[54], aux,	prime);
  mult_aux_gmp(x[111], precompute_pow_omega[118], aux,	prime);
  mult_aux_gmp(x[113], precompute_pow_omega[14], aux,	prime);
  mult_aux_gmp(x[115], precompute_pow_omega[78], aux,	prime);
  mult_aux_gmp(x[117], precompute_pow_omega[46], aux,	prime);
  mult_aux_gmp(x[119], precompute_pow_omega[110], aux,	prime);
  mult_aux_gmp(x[121], precompute_pow_omega[30], aux,	prime);
  mult_aux_gmp(x[123], precompute_pow_omega[94], aux,	prime);
  mult_aux_gmp(x[125], precompute_pow_omega[62], aux,	prime);
  mult_aux_gmp(x[127], precompute_pow_omega[126], aux,	prime);
  mult_aux_gmp(x[129], precompute_pow_omega[1], aux,	prime);
  mult_aux_gmp(x[131], precompute_pow_omega[65], aux,	prime);
  mult_aux_gmp(x[133], precompute_pow_omega[33], aux,	prime);
  mult_aux_gmp(x[135], precompute_pow_omega[97], aux,	prime);
  mult_aux_gmp(x[137], precompute_pow_omega[17], aux,	prime);
  mult_aux_gmp(x[139], precompute_pow_omega[81], aux,	prime);
  mult_aux_gmp(x[141], precompute_pow_omega[49], aux,	prime);
  mult_aux_gmp(x[143], precompute_pow_omega[113], aux,	prime);
  mult_aux_gmp(x[145], precompute_pow_omega[9], aux,	prime);
  mult_aux_gmp(x[147], precompute_pow_omega[73], aux,	prime);
  mult_aux_gmp(x[149], precompute_pow_omega[41], aux,	prime);
  mult_aux_gmp(x[151], precompute_pow_omega[105], aux,	prime);
  mult_aux_gmp(x[153], precompute_pow_omega[25], aux,	prime);
  mult_aux_gmp(x[155], precompute_pow_omega[89], aux,	prime);
  mult_aux_gmp(x[157], precompute_pow_omega[57], aux,	prime);
  mult_aux_gmp(x[159], precompute_pow_omega[121], aux,	prime);
  mult_aux_gmp(x[161], precompute_pow_omega[5], aux,	prime);
  mult_aux_gmp(x[163], precompute_pow_omega[69], aux,	prime);
  mult_aux_gmp(x[165], precompute_pow_omega[37], aux,	prime);
  mult_aux_gmp(x[167], precompute_pow_omega[101], aux,	prime);
  mult_aux_gmp(x[169], precompute_pow_omega[21], aux,	prime);
  mult_aux_gmp(x[171], precompute_pow_omega[85], aux,	prime);
  mult_aux_gmp(x[173], precompute_pow_omega[53], aux,	prime);
  mult_aux_gmp(x[175], precompute_pow_omega[117], aux,	prime);
  mult_aux_gmp(x[177], precompute_pow_omega[13], aux,	prime);
  mult_aux_gmp(x[179], precompute_pow_omega[77], aux,	prime);
  mult_aux_gmp(x[181], precompute_pow_omega[45], aux,	prime);
  mult_aux_gmp(x[183], precompute_pow_omega[109], aux,	prime);
  mult_aux_gmp(x[185], precompute_pow_omega[29], aux,	prime);
  mult_aux_gmp(x[187], precompute_pow_omega[93], aux,	prime);
  mult_aux_gmp(x[189], precompute_pow_omega[61], aux,	prime);
  mult_aux_gmp(x[191], precompute_pow_omega[125], aux,	prime);
  mult_aux_gmp(x[193], precompute_pow_omega[3], aux,	prime);
  mult_aux_gmp(x[195], precompute_pow_omega[67], aux,	prime);
  mult_aux_gmp(x[197], precompute_pow_omega[35], aux,	prime);
  mult_aux_gmp(x[199], precompute_pow_omega[99], aux,	prime);
  mult_aux_gmp(x[201], precompute_pow_omega[19], aux,	prime);
  mult_aux_gmp(x[203], precompute_pow_omega[83], aux,	prime);
  mult_aux_gmp(x[205], precompute_pow_omega[51], aux,	prime);
  mult_aux_gmp(x[207], precompute_pow_omega[115], aux,	prime);
  mult_aux_gmp(x[209], precompute_pow_omega[11], aux,	prime);
  mult_aux_gmp(x[211], precompute_pow_omega[75], aux,	prime);
  mult_aux_gmp(x[213], precompute_pow_omega[43], aux,	prime);
  mult_aux_gmp(x[215], precompute_pow_omega[107], aux,	prime);
  mult_aux_gmp(x[217], precompute_pow_omega[27], aux,	prime);
  mult_aux_gmp(x[219], precompute_pow_omega[91], aux,	prime);
  mult_aux_gmp(x[221], precompute_pow_omega[59], aux,	prime);
  mult_aux_gmp(x[223], precompute_pow_omega[123], aux,	prime);
  mult_aux_gmp(x[225], precompute_pow_omega[7], aux,	prime);
  mult_aux_gmp(x[227], precompute_pow_omega[71], aux,	prime);
  mult_aux_gmp(x[229], precompute_pow_omega[39], aux,	prime);
  mult_aux_gmp(x[231], precompute_pow_omega[103], aux,	prime);
  mult_aux_gmp(x[233], precompute_pow_omega[23], aux,	prime);
  mult_aux_gmp(x[235], precompute_pow_omega[87], aux,	prime);
  mult_aux_gmp(x[237], precompute_pow_omega[55], aux,	prime);
  mult_aux_gmp(x[239], precompute_pow_omega[119], aux,	prime);
  mult_aux_gmp(x[241], precompute_pow_omega[15], aux,	prime);
  mult_aux_gmp(x[243], precompute_pow_omega[79], aux,	prime);
  mult_aux_gmp(x[245], precompute_pow_omega[47], aux,	prime);
  mult_aux_gmp(x[247], precompute_pow_omega[111], aux,	prime);
  mult_aux_gmp(x[249], precompute_pow_omega[31], aux,	prime);
  mult_aux_gmp(x[251], precompute_pow_omega[95], aux,	prime);
  mult_aux_gmp(x[253], precompute_pow_omega[63], aux,	prime);
  mult_aux_gmp(x[255], precompute_pow_omega[127], aux,	prime);
  DFT2_aux_gmp(x[0], x[1], aux, prime);
  DFT2_aux_gmp(x[2], x[3], aux, prime);
  DFT2_aux_gmp(x[4], x[5], aux, prime);
  DFT2_aux_gmp(x[6], x[7], aux, prime);
  DFT2_aux_gmp(x[8], x[9], aux, prime);
  DFT2_aux_gmp(x[10], x[11], aux, prime);
  DFT2_aux_gmp(x[12], x[13], aux, prime);
  DFT2_aux_gmp(x[14], x[15], aux, prime);
  DFT2_aux_gmp(x[16], x[17], aux, prime);
  DFT2_aux_gmp(x[18], x[19], aux, prime);
  DFT2_aux_gmp(x[20], x[21], aux, prime);
  DFT2_aux_gmp(x[22], x[23], aux, prime);
  DFT2_aux_gmp(x[24], x[25], aux, prime);
  DFT2_aux_gmp(x[26], x[27], aux, prime);
  DFT2_aux_gmp(x[28], x[29], aux, prime);
  DFT2_aux_gmp(x[30], x[31], aux, prime);
  DFT2_aux_gmp(x[32], x[33], aux, prime);
  DFT2_aux_gmp(x[34], x[35], aux, prime);
  DFT2_aux_gmp(x[36], x[37], aux, prime);
  DFT2_aux_gmp(x[38], x[39], aux, prime);
  DFT2_aux_gmp(x[40], x[41], aux, prime);
  DFT2_aux_gmp(x[42], x[43], aux, prime);
  DFT2_aux_gmp(x[44], x[45], aux, prime);
  DFT2_aux_gmp(x[46], x[47], aux, prime);
  DFT2_aux_gmp(x[48], x[49], aux, prime);
  DFT2_aux_gmp(x[50], x[51], aux, prime);
  DFT2_aux_gmp(x[52], x[53], aux, prime);
  DFT2_aux_gmp(x[54], x[55], aux, prime);
  DFT2_aux_gmp(x[56], x[57], aux, prime);
  DFT2_aux_gmp(x[58], x[59], aux, prime);
  DFT2_aux_gmp(x[60], x[61], aux, prime);
  DFT2_aux_gmp(x[62], x[63], aux, prime);
  DFT2_aux_gmp(x[64], x[65], aux, prime);
  DFT2_aux_gmp(x[66], x[67], aux, prime);
  DFT2_aux_gmp(x[68], x[69], aux, prime);
  DFT2_aux_gmp(x[70], x[71], aux, prime);
  DFT2_aux_gmp(x[72], x[73], aux, prime);
  DFT2_aux_gmp(x[74], x[75], aux, prime);
  DFT2_aux_gmp(x[76], x[77], aux, prime);
  DFT2_aux_gmp(x[78], x[79], aux, prime);
  DFT2_aux_gmp(x[80], x[81], aux, prime);
  DFT2_aux_gmp(x[82], x[83], aux, prime);
  DFT2_aux_gmp(x[84], x[85], aux, prime);
  DFT2_aux_gmp(x[86], x[87], aux, prime);
  DFT2_aux_gmp(x[88], x[89], aux, prime);
  DFT2_aux_gmp(x[90], x[91], aux, prime);
  DFT2_aux_gmp(x[92], x[93], aux, prime);
  DFT2_aux_gmp(x[94], x[95], aux, prime);
  DFT2_aux_gmp(x[96], x[97], aux, prime);
  DFT2_aux_gmp(x[98], x[99], aux, prime);
  DFT2_aux_gmp(x[100], x[101], aux, prime);
  DFT2_aux_gmp(x[102], x[103], aux, prime);
  DFT2_aux_gmp(x[104], x[105], aux, prime);
  DFT2_aux_gmp(x[106], x[107], aux, prime);
  DFT2_aux_gmp(x[108], x[109], aux, prime);
  DFT2_aux_gmp(x[110], x[111], aux, prime);
  DFT2_aux_gmp(x[112], x[113], aux, prime);
  DFT2_aux_gmp(x[114], x[115], aux, prime);
  DFT2_aux_gmp(x[116], x[117], aux, prime);
  DFT2_aux_gmp(x[118], x[119], aux, prime);
  DFT2_aux_gmp(x[120], x[121], aux, prime);
  DFT2_aux_gmp(x[122], x[123], aux, prime);
  DFT2_aux_gmp(x[124], x[125], aux, prime);
  DFT2_aux_gmp(x[126], x[127], aux, prime);
  DFT2_aux_gmp(x[128], x[129], aux, prime);
  DFT2_aux_gmp(x[130], x[131], aux, prime);
  DFT2_aux_gmp(x[132], x[133], aux, prime);
  DFT2_aux_gmp(x[134], x[135], aux, prime);
  DFT2_aux_gmp(x[136], x[137], aux, prime);
  DFT2_aux_gmp(x[138], x[139], aux, prime);
  DFT2_aux_gmp(x[140], x[141], aux, prime);
  DFT2_aux_gmp(x[142], x[143], aux, prime);
  DFT2_aux_gmp(x[144], x[145], aux, prime);
  DFT2_aux_gmp(x[146], x[147], aux, prime);
  DFT2_aux_gmp(x[148], x[149], aux, prime);
  DFT2_aux_gmp(x[150], x[151], aux, prime);
  DFT2_aux_gmp(x[152], x[153], aux, prime);
  DFT2_aux_gmp(x[154], x[155], aux, prime);
  DFT2_aux_gmp(x[156], x[157], aux, prime);
  DFT2_aux_gmp(x[158], x[159], aux, prime);
  DFT2_aux_gmp(x[160], x[161], aux, prime);
  DFT2_aux_gmp(x[162], x[163], aux, prime);
  DFT2_aux_gmp(x[164], x[165], aux, prime);
  DFT2_aux_gmp(x[166], x[167], aux, prime);
  DFT2_aux_gmp(x[168], x[169], aux, prime);
  DFT2_aux_gmp(x[170], x[171], aux, prime);
  DFT2_aux_gmp(x[172], x[173], aux, prime);
  DFT2_aux_gmp(x[174], x[175], aux, prime);
  DFT2_aux_gmp(x[176], x[177], aux, prime);
  DFT2_aux_gmp(x[178], x[179], aux, prime);
  DFT2_aux_gmp(x[180], x[181], aux, prime);
  DFT2_aux_gmp(x[182], x[183], aux, prime);
  DFT2_aux_gmp(x[184], x[185], aux, prime);
  DFT2_aux_gmp(x[186], x[187], aux, prime);
  DFT2_aux_gmp(x[188], x[189], aux, prime);
  DFT2_aux_gmp(x[190], x[191], aux, prime);
  DFT2_aux_gmp(x[192], x[193], aux, prime);
  DFT2_aux_gmp(x[194], x[195], aux, prime);
  DFT2_aux_gmp(x[196], x[197], aux, prime);
  DFT2_aux_gmp(x[198], x[199], aux, prime);
  DFT2_aux_gmp(x[200], x[201], aux, prime);
  DFT2_aux_gmp(x[202], x[203], aux, prime);
  DFT2_aux_gmp(x[204], x[205], aux, prime);
  DFT2_aux_gmp(x[206], x[207], aux, prime);
  DFT2_aux_gmp(x[208], x[209], aux, prime);
  DFT2_aux_gmp(x[210], x[211], aux, prime);
  DFT2_aux_gmp(x[212], x[213], aux, prime);
  DFT2_aux_gmp(x[214], x[215], aux, prime);
  DFT2_aux_gmp(x[216], x[217], aux, prime);
  DFT2_aux_gmp(x[218], x[219], aux, prime);
  DFT2_aux_gmp(x[220], x[221], aux, prime);
  DFT2_aux_gmp(x[222], x[223], aux, prime);
  DFT2_aux_gmp(x[224], x[225], aux, prime);
  DFT2_aux_gmp(x[226], x[227], aux, prime);
  DFT2_aux_gmp(x[228], x[229], aux, prime);
  DFT2_aux_gmp(x[230], x[231], aux, prime);
  DFT2_aux_gmp(x[232], x[233], aux, prime);
  DFT2_aux_gmp(x[234], x[235], aux, prime);
  DFT2_aux_gmp(x[236], x[237], aux, prime);
  DFT2_aux_gmp(x[238], x[239], aux, prime);
  DFT2_aux_gmp(x[240], x[241], aux, prime);
  DFT2_aux_gmp(x[242], x[243], aux, prime);
  DFT2_aux_gmp(x[244], x[245], aux, prime);
  DFT2_aux_gmp(x[246], x[247], aux, prime);
  DFT2_aux_gmp(x[248], x[249], aux, prime);
  DFT2_aux_gmp(x[250], x[251], aux, prime);
  DFT2_aux_gmp(x[252], x[253], aux, prime);
  DFT2_aux_gmp(x[254], x[255], aux, prime);
  swap_aux_gmp(x[1], x[128], aux);
  swap_aux_gmp(x[2], x[64], aux);
  swap_aux_gmp(x[3], x[192], aux);
  swap_aux_gmp(x[4], x[32], aux);
  swap_aux_gmp(x[5], x[160], aux);
  swap_aux_gmp(x[6], x[96], aux);
  swap_aux_gmp(x[7], x[224], aux);
  swap_aux_gmp(x[8], x[16], aux);
  swap_aux_gmp(x[9], x[144], aux);
  swap_aux_gmp(x[10], x[80], aux);
  swap_aux_gmp(x[11], x[208], aux);
  swap_aux_gmp(x[12], x[48], aux);
  swap_aux_gmp(x[13], x[176], aux);
  swap_aux_gmp(x[14], x[112], aux);
  swap_aux_gmp(x[15], x[240], aux);
  swap_aux_gmp(x[17], x[136], aux);
  swap_aux_gmp(x[18], x[72], aux);
  swap_aux_gmp(x[19], x[200], aux);
  swap_aux_gmp(x[20], x[40], aux);
  swap_aux_gmp(x[21], x[168], aux);
  swap_aux_gmp(x[22], x[104], aux);
  swap_aux_gmp(x[23], x[232], aux);
  swap_aux_gmp(x[25], x[152], aux);
  swap_aux_gmp(x[26], x[88], aux);
  swap_aux_gmp(x[27], x[216], aux);
  swap_aux_gmp(x[28], x[56], aux);
  swap_aux_gmp(x[29], x[184], aux);
  swap_aux_gmp(x[30], x[120], aux);
  swap_aux_gmp(x[31], x[248], aux);
  swap_aux_gmp(x[33], x[132], aux);
  swap_aux_gmp(x[34], x[68], aux);
  swap_aux_gmp(x[35], x[196], aux);
  swap_aux_gmp(x[37], x[164], aux);
  swap_aux_gmp(x[38], x[100], aux);
  swap_aux_gmp(x[39], x[228], aux);
  swap_aux_gmp(x[41], x[148], aux);
  swap_aux_gmp(x[42], x[84], aux);
  swap_aux_gmp(x[43], x[212], aux);
  swap_aux_gmp(x[44], x[52], aux);
  swap_aux_gmp(x[45], x[180], aux);
  swap_aux_gmp(x[46], x[116], aux);
  swap_aux_gmp(x[47], x[244], aux);
  swap_aux_gmp(x[49], x[140], aux);
  swap_aux_gmp(x[50], x[76], aux);
  swap_aux_gmp(x[51], x[204], aux);
  swap_aux_gmp(x[53], x[172], aux);
  swap_aux_gmp(x[54], x[108], aux);
  swap_aux_gmp(x[55], x[236], aux);
  swap_aux_gmp(x[57], x[156], aux);
  swap_aux_gmp(x[58], x[92], aux);
  swap_aux_gmp(x[59], x[220], aux);
  swap_aux_gmp(x[61], x[188], aux);
  swap_aux_gmp(x[62], x[124], aux);
  swap_aux_gmp(x[63], x[252], aux);
  swap_aux_gmp(x[65], x[130], aux);
  swap_aux_gmp(x[67], x[194], aux);
  swap_aux_gmp(x[69], x[162], aux);
  swap_aux_gmp(x[70], x[98], aux);
  swap_aux_gmp(x[71], x[226], aux);
  swap_aux_gmp(x[73], x[146], aux);
  swap_aux_gmp(x[74], x[82], aux);
  swap_aux_gmp(x[75], x[210], aux);
  swap_aux_gmp(x[77], x[178], aux);
  swap_aux_gmp(x[78], x[114], aux);
  swap_aux_gmp(x[79], x[242], aux);
  swap_aux_gmp(x[81], x[138], aux);
  swap_aux_gmp(x[83], x[202], aux);
  swap_aux_gmp(x[85], x[170], aux);
  swap_aux_gmp(x[86], x[106], aux);
  swap_aux_gmp(x[87], x[234], aux);
  swap_aux_gmp(x[89], x[154], aux);
  swap_aux_gmp(x[91], x[218], aux);
  swap_aux_gmp(x[93], x[186], aux);
  swap_aux_gmp(x[94], x[122], aux);
  swap_aux_gmp(x[95], x[250], aux);
  swap_aux_gmp(x[97], x[134], aux);
  swap_aux_gmp(x[99], x[198], aux);
  swap_aux_gmp(x[101], x[166], aux);
  swap_aux_gmp(x[103], x[230], aux);
  swap_aux_gmp(x[105], x[150], aux);
  swap_aux_gmp(x[107], x[214], aux);
  swap_aux_gmp(x[109], x[182], aux);
  swap_aux_gmp(x[110], x[118], aux);
  swap_aux_gmp(x[111], x[246], aux);
  swap_aux_gmp(x[113], x[142], aux);
  swap_aux_gmp(x[115], x[206], aux);
  swap_aux_gmp(x[117], x[174], aux);
  swap_aux_gmp(x[119], x[238], aux);
  swap_aux_gmp(x[121], x[158], aux);
  swap_aux_gmp(x[123], x[222], aux);
  swap_aux_gmp(x[125], x[190], aux);
  swap_aux_gmp(x[127], x[254], aux);
  swap_aux_gmp(x[131], x[193], aux);
  swap_aux_gmp(x[133], x[161], aux);
  swap_aux_gmp(x[135], x[225], aux);
  swap_aux_gmp(x[137], x[145], aux);
  swap_aux_gmp(x[139], x[209], aux);
  swap_aux_gmp(x[141], x[177], aux);
  swap_aux_gmp(x[143], x[241], aux);
  swap_aux_gmp(x[147], x[201], aux);
  swap_aux_gmp(x[149], x[169], aux);
  swap_aux_gmp(x[151], x[233], aux);
  swap_aux_gmp(x[155], x[217], aux);
  swap_aux_gmp(x[157], x[185], aux);
  swap_aux_gmp(x[159], x[249], aux);
  swap_aux_gmp(x[163], x[197], aux);
  swap_aux_gmp(x[167], x[229], aux);
  swap_aux_gmp(x[171], x[213], aux);
  swap_aux_gmp(x[173], x[181], aux);
  swap_aux_gmp(x[175], x[245], aux);
  swap_aux_gmp(x[179], x[205], aux);
  swap_aux_gmp(x[183], x[237], aux);
  swap_aux_gmp(x[187], x[221], aux);
  swap_aux_gmp(x[191], x[253], aux);
  swap_aux_gmp(x[199], x[227], aux);
  swap_aux_gmp(x[203], x[211], aux);
  swap_aux_gmp(x[207], x[243], aux);
  swap_aux_gmp(x[215], x[235], aux);
  swap_aux_gmp(x[223], x[251], aux);
  swap_aux_gmp(x[239], x[247], aux);
}


/**************************************/

void
stride_permutation_gmp (int n_permutations, mpz_t* A, int m, int n)
{
#if PROFILING_GMP_ENABLED ==1
  stride_permutation_gmp_called_dims[n_stride_permutation_gmp_called][0] = m;
  stride_permutation_gmp_called_dims[n_stride_permutation_gmp_called][1] = n;
  inc_gmp_profiling_counter (&n_stride_permutation_gmp_called);
#endif

  int blocksize = 1;
  mpz_t* B = (mpz_t*) malloc (n_permutations*m * n * sizeof(mpz_t));
//  for (int i = 0; i < m * n; i++)
//    mpz_init_se (B[i]);

  int stride = m * n;
  usfixn64 mn=m*n;
  int m_log2=log2(m);


#if GMP_SERIAL_ENABLED == 1
  for (int t = 0; t < n_permutations; t++)
    for (int ij = 0; ij < mn; ij++)
#else
#pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
  cilk_for (int t = 0; t < n_permutations; t++)
      cilk_for (int ij = 0; ij < mn; ij++)
#endif
      {
	int i=ij>>m_log2;
	int j=(ij&(m-1));
	  {
	    //  mpz_init(B[k+l*n]);
	    mpz_init_set (B[i + j * n+t*stride], A[j + i * m + t * stride]);
	  }
      }

#if GMP_SERIAL_ENABLED == 1
  for (int t = 0; t < n_permutations; t++)
    for (int i = 0; i < m * n; i++)
#else
    cilk_for (int t = 0; t < n_permutations; t++)
        cilk_for (int i = 0; i < m * n; i++)
#endif
      {
	mpz_set (A[i + t * stride], B[i+t*stride]);
	mpz_clear (B[i + t * stride]);
      }

//  #pragma cilk_grainsize = ((m*n)/(cilk::current_worker_count))
//  cilk_for (int i = 0; i < m * n*n_permutations; i++)
//    mpz_clear (B[i]);
  free (B);
}

/**************************************/

void
stride_permutation_aux_gmp (int n_permutations, mpz_t* A, mpz_t* B, int m, int n)
{
#if PROFILING_GMP_ENABLED ==1
  stride_permutation_gmp_called_dims[n_stride_permutation_gmp_called][0] = m;
  stride_permutation_gmp_called_dims[n_stride_permutation_gmp_called][1] = n;
  inc_gmp_profiling_counter (&n_stride_permutation_gmp_called);
#endif

  int blocksize = 1;
//  mpz_t* B = (mpz_t*) malloc (n_permutations*m * n * sizeof(mpz_t));
//  for (int i = 0; i < m * n; i++)
//    mpz_init_se (B[i]);

  int stride = m * n;
  usfixn64 mn=m*n;
  int m_log2=log2(m);


#if GMP_SERIAL_ENABLED == 1
  for (int t = 0; t < n_permutations; t++)
    for (int ij = 0; ij < mn; ij++)
#else
#pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
  for (int t = 0; t < n_permutations; t++)
	for (int ij = 0; ij < mn; ij++)
#endif
      {
	int i=ij>>m_log2;
	int j=(ij&(m-1));
	  {
	    //  mpz_init(B[k+l*n]);
	    mpz_set (B[i + j * n+t*stride], A[j + i * m + t * stride]);
	  }
      }

#if GMP_SERIAL_ENABLED == 1
  for (int t = 0; t < n_permutations; t++)
    for (int i = 0; i < mn; i++)
#else
  for (int t = 0; t < n_permutations; t++)
    for (int i = 0; i < mn; i++)
#endif
      {
	mpz_set (A[i + t * stride], B[i+t*stride]);
//	mpz_clear (B[i + t * stride]);
      }

//  #pragma cilk_grainsize = ((m*n)/(cilk::current_worker_count))
//  cilk_for (int i = 0; i < m * n*n_permutations; i++)
//    mpz_clear (B[i]);
//  free (B);
}

/**************************************/

/**************************************/
//precompute powers of omega up to (omega)^{K^b}
//to be used at level b+1
//sfixn64* precomputed_pow_omega_vec should be be already allocated!
//just pass a pointer.
void
precompute_sequential_powers_of_omega_gmp (mpz_t* precomputed_pow_omega_vec,
					   const int K, const int b,
					   const mpz_t omega, const mpz_t p)
{
  int n = pow (K, b);
  mpz_t current_twiddle;
  mpz_init_set_ui (current_twiddle, 1);

  for (int i = 0; i < n; i++)
    {
      if (i != 0)
	{
	  mpz_mul (current_twiddle, current_twiddle, omega);
	  mpz_mod (current_twiddle, current_twiddle, p);
	}
//      mpz_init_set (precomputed_pow_omega_vec[i], current_twiddle);
      mpz_set (precomputed_pow_omega_vec[i], current_twiddle);
    }
  mpz_clear (current_twiddle);
}

/**************************************/

void
mult_gmp_for_profiling_only (mpz_t X, mpz_t Y, mpz_t omega, mpz_t P)
{
  mpz_mul (X, X, Y);
  mpz_mod (X, X, P);
}

/**************************************/

void
twiddle_gmp (int n_permutations, mpz_t* vector, int K, int J,
	     const mpz_t * omega_w, const mpz_t *omega_base_precomputed_vec,
	     const mpz_t prime, float * t_twiddle_values, int *n_twiddle_mult)
{

#if PROFILING_GMP_ENABLED ==1
  twiddle_gmp_called_dims[n_twiddle_gmp_called][0] = K;
  twiddle_gmp_called_dims[n_twiddle_gmp_called][1] = J;
  inc_gmp_profiling_counter (&n_twiddle_gmp_called);
#endif



  mpz_t radix_zz;
  mpz_init (radix_zz);
  mpz_set (radix_zz, omega_base_precomputed_vec[1]);

  int J_log2 = log2 (J);

  cpu_timer t_mult;
  float t_mul_pow_r = 0, t_bigint_mult = 0;

  n_twiddle_mult[0] += n_permutations * K * J;
  timer_record_start (&t_mult);

  usfixn64 KJ=K*J;
  usfixn64 J_mask=J-1;

#if GMP_SERIAL_ENABLED == 1
  for (int iteration = 0; iteration < n_permutations; iteration++)
    for (int ij = 0; ij < KJ; ij++)
#else
  #pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
    cilk_for (int iteration = 0; iteration < n_permutations; iteration++)
      cilk_for (int ij = 0; ij < KJ; ij++)
#endif
	{
	  usfixn64 stride = iteration * KJ;
	  int i,j;
	  i=ij>>J_log2;
	  j=ij&(J_mask);

	  //point-wise multiplication
	  //int idx=j*m+i;
	  int idx = i * j;
  //	int q = idx / J;
	  int q = idx >> J_log2;
	  mpz_t t;
	  mpz_init (t);
	  mpz_pow_ui (t, radix_zz, q);
	  mpz_mul (vector[i * J + j + stride], vector[i * J + j + stride], t);
	  mpz_mod (vector[i * J + j + stride], vector[i * J + j + stride],
		   prime);
	  mpz_clear (t);
	}

  timer_record_stop (&t_mult);
  timer_get_elapsed_time (&t_mult, NULL, 1);
  t_mul_pow_r = t_mult.elapsed_time;

  timer_record_start (&t_mult);

#if GMP_SERIAL_ENABLED == 1
  for (int iteration = 0; iteration < n_permutations; iteration++)
    for (int ij = 0; ij < KJ; ij++)
#else
#pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
  cilk_for (int iteration = 0; iteration < n_permutations; iteration++)
    cilk_for (int ij = 0; ij < KJ; ij++)
#endif
      {
	usfixn64 stride = iteration * KJ;
	int i,j;
	i=ij>>J_log2;
	j=ij&(J_mask);
	//int idx=j*m+i;
	int idx = i * j;
//  	int q = idx / J;
	int h = idx & (J_mask);

//	if (h!=0)
	  {
	  mpz_t t;
	  mpz_init (t);
	    //      mpz_set(t, omega_w);
	    //      mpz_pow_ui(t, omega_w, idx);
	    //      mpz_mod(t, t, prime);
	    //      mpz_powm_ui(t, omega_w, idx, prime);
	    //      mpz_powm_ui(t, omega_w, h, prime);
	    mpz_set (t, omega_w[h]);

	    //      mpz_set(t1, omega_base);
	    //      mpz_powm_ui(t1, omega_base, q, prime);
	    //      mpz_mul(t, t, t1);
	    //t=(omega_w)^(i*j);
	    //      mpz_powm_ui(t, omega_w, (long int) (i * j), prime);
	    // vector[j*m+i]=vector[j*m+i]*(omega_w^(i*j));
	    //mult_gmp(vector[i * J + j], t, prime);
	    mpz_mul (vector[i * J + j + stride], vector[i * J + j + stride],
		     t);
	    mpz_mod (vector[i * J + j + stride], vector[i * J + j + stride],
		     prime);
	    mpz_clear (t);
	  }
      }

  n_twiddle_mult[1] += n_permutations * K * J;

  timer_record_stop (&t_mult);
  timer_get_elapsed_time (&t_mult, NULL, 1);
  t_bigint_mult = t_mult.elapsed_time;

  mpz_clear (radix_zz);

  t_twiddle_values[0] = t_mul_pow_r;
  t_twiddle_values[1] = t_bigint_mult;
}

/**************************************/

void
twiddle_aux_gmp (int n_permutations, mpz_t* vector, mpz_t* aux, int K, int J,
	     const mpz_t * omega_w, const mpz_t *omega_base_precomputed_vec,
	     const mpz_t prime, float * t_twiddle_values, int *n_twiddle_mult)
{

#if PROFILING_GMP_ENABLED ==1
  twiddle_gmp_called_dims[n_twiddle_gmp_called][0] = K;
  twiddle_gmp_called_dims[n_twiddle_gmp_called][1] = J;
  inc_gmp_profiling_counter (&n_twiddle_gmp_called);
#endif



  mpz_t radix_zz;
  mpz_init (radix_zz);
  mpz_set (radix_zz, omega_base_precomputed_vec[1]);

  int J_log2 = log2 (J);

  cpu_timer t_mult;
  float t_mul_pow_r = 0, t_bigint_mult = 0;

  n_twiddle_mult[0] += n_permutations * K * J;
  timer_record_start (&t_mult);

  usfixn64 KJ=K*J;
  usfixn64 J_mask=J-1;

#if GMP_SERIAL_ENABLED == 1
  for (int iteration = 0; iteration < n_permutations; iteration++)
    for (int ij = 0; ij < KJ; ij++)
#else
  #pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
    cilk_for (int iteration = 0; iteration < n_permutations; iteration++)
      cilk_for (int ij = 0; ij < KJ; ij++)
#endif
	{
	  usfixn64 stride = iteration * KJ;
	  int i,j;
	  i=ij>>J_log2;
	  j=ij&(J_mask);

	  //point-wise multiplication
	  //int idx=j*m+i;
	  int idx = i * j;
  //	int q = idx / J;
	  int q = idx >> J_log2;
//	  mpz_t t;
//	  mpz_init (t);
	  mpz_pow_ui (aux[i * J + j + stride], radix_zz, q);
	  mpz_mul (vector[i * J + j + stride], vector[i * J + j + stride], aux[i * J + j + stride]);
	  mpz_mod (vector[i * J + j + stride], vector[i * J + j + stride],
		   prime);
//	  mpz_clear (t);
	}

  timer_record_stop (&t_mult);
  timer_get_elapsed_time (&t_mult, NULL, 1);
  t_mul_pow_r = t_mult.elapsed_time;

  timer_record_start (&t_mult);

#if GMP_SERIAL_ENABLED == 1
  for (int iteration = 0; iteration < n_permutations; iteration++)
    for (int ij = 0; ij < KJ; ij++)
#else
#pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
  cilk_for (int iteration = 0; iteration < n_permutations; iteration++)
    cilk_for (int ij = 0; ij < KJ; ij++)
#endif
      {
	usfixn64 stride = iteration * KJ;
	int i,j;
	i=ij>>J_log2;
	j=ij&(J_mask);
	//int idx=j*m+i;
	int idx = i * j;
//  	int q = idx / J;
	int h = idx & (J_mask);

//	if (h!=0)
	  {
//	  mpz_t t;
//	  mpz_init (t);
	    //      mpz_set(t, omega_w);
	    //      mpz_pow_ui(t, omega_w, idx);
	    //      mpz_mod(t, t, prime);
	    //      mpz_powm_ui(t, omega_w, idx, prime);
	    //      mpz_powm_ui(t, omega_w, h, prime);
	    mpz_set (aux[i * J + j + stride], omega_w[h]);

	    //      mpz_set(t1, omega_base);
	    //      mpz_powm_ui(t1, omega_base, q, prime);
	    //      mpz_mul(t, t, t1);
	    //t=(omega_w)^(i*j);
	    //      mpz_powm_ui(t, omega_w, (long int) (i * j), prime);
	    // vector[j*m+i]=vector[j*m+i]*(omega_w^(i*j));
	    //mult_gmp(vector[i * J + j], t, prime);
	    mpz_mul (vector[i * J + j + stride], vector[i * J + j + stride],
		     aux[i * J + j + stride]);
	    mpz_mod (vector[i * J + j + stride], vector[i * J + j + stride],
		     prime);
//	    mpz_clear (t);
	  }
      }

  n_twiddle_mult[1] += n_permutations * K * J;

  timer_record_stop (&t_mult);
  timer_get_elapsed_time (&t_mult, NULL, 1);
  t_bigint_mult = t_mult.elapsed_time;

  mpz_clear (radix_zz);

  t_twiddle_values[0] = t_mul_pow_r;
  t_twiddle_values[1] = t_bigint_mult;
}

/**************************************/

void
DFT_general_gmp (mpz_t* vector, int K, int e, const mpz_t omega_zz, mpz_t prime,
		 int verbose)
{
  usfixn64 N, J;
  N = pow (K, e);
  J = pow (K, e - 1);

  mpz_t omega_base;
  mpz_init_set (omega_base, omega_zz);
  mpz_pow_ui (omega_base, omega_base, pow (K, e - 1));
  mpz_mod (omega_base, omega_base, prime);

  mpz_t current_omega_zz;
  mpz_init (current_omega_zz);
  mpz_set (current_omega_zz, omega_zz);

  // precomputed powers of omega in the base case;
  mpz_t * precomputed_pow_omega_base_vec = (mpz_t*) malloc (K * sizeof(mpz_t));

  float t_stride_total = 0, t_DFT_K_total = 0, t_twiddle_total = 0,
      t_precompute_total, t_mult_pow_r_total = 0, t_arbitrary_mult_total = 0;

  //n_twiddle_mult: [0]:n_mult_pow_r [1]: n_arbitrary_mult
  int n_twiddle_mult[2] =
    { 0, 0 }, n_DFT_K = 0, n_mult_pow_r_in_DFT_K = 0,
      n_mult_pow_r_per_DFT_K = 0;

  cpu_timer t_precompute, t_twiddle, t_DFT_K, t_permutation;

  /////////////////////////////////////
  /////////////////////////////////////

  timer_record_start (&t_precompute);
  for (int i = 0; i < K; i++)
    {
      mpz_init (precomputed_pow_omega_base_vec[i]);
      mpz_powm_ui (precomputed_pow_omega_base_vec[i], omega_base, i, prime);
    }

  /////////////////////////////////////
  ////precomputation of powers of omega: VERIFIED.
  /////////////////////////////////////
  mpz_t ** precomputed_pow_omega_vec_2d = (mpz_t**) malloc (
      (e + 1) * sizeof(mpz_t*));
  mpz_set (current_omega_zz, omega_zz);
  for (int b = e; b >= 2; b--)
    {
//    printf("in precomputing pow omega b=%d\n", b);
      J = pow (K, b - 1);
      precomputed_pow_omega_vec_2d[b] = (mpz_t*) malloc (J * sizeof(mpz_t));

      for (int j = 0; j < J; j++)
	mpz_init (precomputed_pow_omega_vec_2d[b][j]);

      //precompute sequential powers of omega at level=b and store it in 2D vector.
      precompute_sequential_powers_of_omega_gmp (
	  precomputed_pow_omega_vec_2d[b], K, b - 1, current_omega_zz, prime);

      //as we go to a lower level, w=(w^K) mod p
      mpz_powm_ui (current_omega_zz, current_omega_zz, K, prime);
    }
  timer_record_stop (&t_precompute);
  timer_get_elapsed_time (&t_precompute, NULL, 1);
  t_precompute_total += t_precompute.elapsed_time;
  /////////////////////////////////////

  cpu_timer t_dft_gmp;
  timer_record_start (&t_dft_gmp);

  timer_record_start (&t_permutation);
  // step 1: stride permutations on the right-most side of
  // the DFT, computing L{K,J}
  //step1 is verified, but should be done more rigorously.
  for (int b = e; b >= 2; b--)
    {

      J = pow (K, b - 1);
      int n_permutations = N / (K * J);
//      printf ("[step 1: (%-12s) b=%d, RHS L_{K=%d, J=%llu}]\n", b, K, J);
      stride_permutation_gmp (n_permutations, vector, K, J);
    }

  timer_record_stop (&t_permutation);
  timer_get_elapsed_time (&t_permutation, NULL, 1);
  t_stride_total += t_permutation.elapsed_time;

  /////////////////////////////////////
  /////////////////////////////////////
  ////step 2: base-case DFT's at the lowest level
  /////////////////////////////////////
//  printf ("[step 2: (%-12s) b=%d, right DFT_K {K=%d, J=%llu}]\n", 2, K, J);

  //ToDO: automated verification of base-case DFT's

  usfixn64 stride = (K);
  J = pow (K, e - 1);

  void
  (*DFT_K_gmp_ptr) (mpz_t* a, mpz_t omega, mpz_t prime,
		    mpz_t * precompute_pow_omega);

  if (K == 8)
    {
      DFT_K_gmp_ptr = DFT_8_gmp;
      n_mult_pow_r_per_DFT_K = 5;
    }
  else if (K == 16)
    {
      DFT_K_gmp_ptr = DFT_16_gmp;
      n_mult_pow_r_per_DFT_K = 17;
    }
  else if (K == 32)
    {
      DFT_K_gmp_ptr = DFT_32_gmp;
      n_mult_pow_r_per_DFT_K = 66;
    }
  else if (K == 64)
    {
      DFT_K_gmp_ptr = DFT_64_gmp;
      n_mult_pow_r_per_DFT_K = 195;
    }
//  else if (K == 128)
//      {
//      DFT_K_gmp_ptr = DFT_128_gmp;
//      n_mult_pow_r_per_DFT_K = 195;
//      }

  timer_record_start (&t_DFT_K);


#if GMP_SERIAL_ENABLED==1
  for (int j = 0; j < J; j++)
#else
#pragma cilk_grainsize = (J/(cilk::current_worker_count))
  cilk_for (int j = 0; j < J; j++)
#endif
    {
      DFT_K_gmp_ptr (&vector[j * K], omega_base, prime,
		     precomputed_pow_omega_base_vec);
    }
  timer_record_stop (&t_DFT_K);
  timer_get_elapsed_time (&t_DFT_K, NULL, 1);
  t_DFT_K_total += t_DFT_K.elapsed_time;
  n_DFT_K += J;

  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////
  ////step 3:
  ////- twiddle factor multiplication + permutation in between + DFT_K on the left
  ////- beginning at the level b=2;
  /////////////////////////////////////
  mpz_set (current_omega_zz, omega_zz);
  for (int b = 2; b <= e; b++)
    {
      J = pow (K, b - 1);
      int n_permutations = N / (K * J);

//      printf ("[step 3: (%-12s) b=%d, D_{K=%d, J=%llu}, n_permutations=%d]\n",
//	      b, K, J, n_permutations);
      usfixn64 stride = K * J;
      usfixn64 idx = 0;
      mpz_powm_ui (current_omega_zz, omega_zz, pow (K, e - b), prime);

      float t_twiddle_values[2] =
	{ 0, 0 };

      twiddle_gmp (n_permutations, vector, K, J,
		   precomputed_pow_omega_vec_2d[b],
		   precomputed_pow_omega_base_vec, prime, t_twiddle_values,
		   n_twiddle_mult);
      t_mult_pow_r_total += t_twiddle_values[0];
      t_arbitrary_mult_total += t_twiddle_values[1];

      timer_record_start (&t_permutation);
      stride_permutation_gmp (n_permutations, vector, J, K);
      timer_record_stop (&t_permutation);
      timer_get_elapsed_time (&t_permutation, NULL, 1);
      t_stride_total += t_permutation.elapsed_time;

      /////////////////////////////////
      ////step 4:
      ////computing DFT_K on the left.
      /////////////////////////////////
      J = N / K;
      n_permutations = N / (K * J);
      // printf ("-- step 4: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
      // "left DFT_K", b, K, J, n_permutations);

      timer_record_start (&t_DFT_K);
  #if GMP_SERIAL_ENABLED==1
    for (int j = 0; j < J; j++)
  #else
  #pragma cilk_grainsize = (J/(cilk::current_worker_count))
    cilk_for (int j = 0; j < J; j++)
  #endif
	{
	  DFT_K_gmp_ptr (&vector[j * K], omega_base, prime,
			 precomputed_pow_omega_base_vec);
	}

      timer_record_stop (&t_DFT_K);
      timer_get_elapsed_time (&t_DFT_K, NULL, 1);
      t_DFT_K_total += t_DFT_K.elapsed_time;
      n_DFT_K += J;

      /////////////////////////////////
      //// step 5: left-most permutation L_{K,K^e-1} on the right
      /////////////////////////////////
//    printf("-- step 5: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
//    "left perm", b, K, J, n_permutations);

      J = pow (K, b - 1);
      n_permutations = N / (K * J);

      timer_record_start (&t_permutation);
      stride_permutation_gmp (n_permutations, vector, K, J);
      timer_record_stop (&t_permutation);
      timer_get_elapsed_time (&t_permutation, NULL, 1);
      t_stride_total += t_permutation.elapsed_time;
    }

  timer_record_stop (&t_dft_gmp);
  timer_get_elapsed_time (&t_dft_gmp, NULL, 1);

  float t_mult_pow_r_avg;
  t_mult_pow_r_avg = t_mult_pow_r_total / n_twiddle_mult[0];
  int n_mult_pow_r_in_DFT_K_total = (n_mult_pow_r_per_DFT_K * n_DFT_K);
//  n_twiddle_mult[0] += n_mult_pow_r_in_DFT_K_total;
  float t_mult_pow_r_in_DFT_K_total = (t_mult_pow_r_avg
      * n_mult_pow_r_in_DFT_K_total);

//  t_mult_pow_r_total += t_mult_pow_r_in_DFT_K_total;
  t_twiddle_total = t_mult_pow_r_total + t_arbitrary_mult_total;

  if (verbose)
    {
      char suffix[64];
      char msg[256];

      sprintf (suffix, "%s", "gmp");

      sprintf (msg, "%s_%s", "t_precompute", suffix);
      timer_print_time_percentage (t_precompute_total, msg,
				   t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_stride_total", suffix);
      timer_print_time_percentage (t_stride_total, msg, t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_DFT_K_total", suffix);
      timer_print_time_percentage (t_DFT_K_total, msg, t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_DFT_K_total_avg", suffix);
      timer_print_time (t_DFT_K_total / n_DFT_K, msg);

      printf ("%s", quadline);

      sprintf (msg, "%s_%s", "t_mult_pow_r_in_DFT_K_total", suffix);
      timer_print_time_percentage (t_mult_pow_r_in_DFT_K_total, msg,
				   t_dft_gmp.elapsed_time);

      printf ("%s", quadline);

      sprintf (msg, "%s_%s", "t_mult_pow_r_total", suffix);
      timer_print_time_percentage (t_mult_pow_r_total, msg,
				   t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_mult_pow_r_total_avg", suffix);
      timer_print_time (t_mult_pow_r_total / n_twiddle_mult[0], msg);

      sprintf (msg, "%s_%s", "t_arbitrary_mult_total", suffix);
      timer_print_time_percentage (t_arbitrary_mult_total, msg,
				   t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_arbitrary_mult_total_avg", suffix);
      timer_print_time (t_arbitrary_mult_total / n_twiddle_mult[1], msg);

      sprintf (msg, "%s_%s", "t_twiddle_total", suffix);
      timer_print_time_percentage (t_twiddle_total, msg,
				   t_dft_gmp.elapsed_time);
      printf ("%s", quadline);

      sprintf (msg, "%s_%s", "t_dft_six_step", suffix);
      timer_print_time_percentage (t_dft_gmp.elapsed_time, msg,
				   t_dft_gmp.elapsed_time);
      printf ("%s", longline);
    }

  mpz_clear (current_omega_zz);
  mpz_clear (omega_base);
  for (int i = 0; i < K; i++)
    {
      mpz_clear (precomputed_pow_omega_base_vec[i]);
    }
  free (precomputed_pow_omega_base_vec);

  for (int b = e; b >= 2; b--)
    {
      J = pow (K, b - 1);
      for (int j = 0; j < J; j++)
	mpz_clear (precomputed_pow_omega_vec_2d[b][j]);
      free (precomputed_pow_omega_vec_2d[b]);
    }
  free (precomputed_pow_omega_vec_2d);

}

/**************************************/

void
DFT_general_aux_memory_gmp (mpz_t* vector, int K, int e, const mpz_t omega_zz, mpz_t prime,
		 int verbose)
{
  usfixn64 N, J;
  N = pow (K, e);
  J = pow (K, e - 1);


  int max_coeff_bit_cnt=(K*sizeof(mpz_t)*8);
  mpz_t *aux_zz=(mpz_t*)malloc(N*sizeof(mpz_t));
  for(int i=0;i<N;i++)
    mpz_init(aux_zz[i]);
//    mpz_init2(aux_zz[i],max_coeff_bit_cnt);

  mpz_t omega_base;
  mpz_init_set (omega_base, omega_zz);
  mpz_pow_ui (omega_base, omega_base, pow (K, e - 1));
  mpz_mod (omega_base, omega_base, prime);

  mpz_t current_omega_zz;
  mpz_init (current_omega_zz);
  mpz_set (current_omega_zz, omega_zz);

  // precomputed powers of omega in the base case;
  mpz_t * precomputed_pow_omega_base_vec = (mpz_t*) malloc (K * sizeof(mpz_t));

  float t_stride_total = 0, t_DFT_K_total = 0, t_twiddle_total = 0,
      t_precompute_total, t_mult_pow_r_total = 0, t_arbitrary_mult_total = 0;

  //n_twiddle_mult: [0]:n_mult_pow_r [1]: n_arbitrary_mult
  int n_twiddle_mult[2] =
    { 0, 0 }, n_DFT_K = 0, n_mult_pow_r_in_DFT_K = 0,
      n_mult_pow_r_per_DFT_K = 0;

  cpu_timer t_precompute, t_twiddle, t_DFT_K, t_permutation;

  /////////////////////////////////////
  /////////////////////////////////////

  timer_record_start (&t_precompute);
  for (int i = 0; i < K; i++)
    {
      mpz_init (precomputed_pow_omega_base_vec[i]);
      mpz_powm_ui (precomputed_pow_omega_base_vec[i], omega_base, i, prime);
    }

  /////////////////////////////////////
  ////precomputation of powers of omega: VERIFIED.
  /////////////////////////////////////
  mpz_t ** precomputed_pow_omega_vec_2d = (mpz_t**) malloc (
      (e + 1) * sizeof(mpz_t*));
  mpz_set (current_omega_zz, omega_zz);
  for (int b = e; b >= 2; b--)
    {
//    printf("in precomputing pow omega b=%d\n", b);
      J = pow (K, b - 1);
      precomputed_pow_omega_vec_2d[b] = (mpz_t*) malloc (J * sizeof(mpz_t));

      for (int j = 0; j < J; j++)
	mpz_init (precomputed_pow_omega_vec_2d[b][j]);

      //precompute sequential powers of omega at level=b and store it in 2D vector.
      precompute_sequential_powers_of_omega_gmp (
	  precomputed_pow_omega_vec_2d[b], K, b - 1, current_omega_zz, prime);

      //as we go to a lower level, w=(w^K) mod p
      mpz_powm_ui (current_omega_zz, current_omega_zz, K, prime);
    }
  timer_record_stop (&t_precompute);
  timer_get_elapsed_time (&t_precompute, NULL, 1);
  t_precompute_total += t_precompute.elapsed_time;
  /////////////////////////////////////

  cpu_timer t_dft_gmp;
  timer_record_start (&t_dft_gmp);

  timer_record_start (&t_permutation);
  // step 1: stride permutations on the right-most side of
  // the DFT, computing L{K,J}
  //step1 is verified, but should be done more rigorously.
  for (int b = e; b >= 2; b--)
    {
      J = pow (K, b - 1);
      int n_permutations = N / (K * J);
//      printf ("[step 1: (%-12s) b=%d, RHS L_{K=%d, J=%llu}]\n", b, K, J);
      stride_permutation_aux_gmp(n_permutations, vector, aux_zz, K, J);
    }

  timer_record_stop (&t_permutation);
  timer_get_elapsed_time (&t_permutation, NULL, 1);
  t_stride_total += t_permutation.elapsed_time;

  /////////////////////////////////////
  /////////////////////////////////////
  ////step 2: base-case DFT's at the lowest level
  /////////////////////////////////////
//  printf ("[step 2: (%-12s) b=%d, right DFT_K {K=%d, J=%llu}]\n", 2, K, J);

  //ToDO: automated verification of base-case DFT's

  usfixn64 stride = (K);
  J = pow (K, e - 1);

  void
  (*DFT_K_gmp_ptr) (mpz_t* a, mpz_t a_aux, mpz_t omega, mpz_t prime,
		    mpz_t * precompute_pow_omega);

  if (K == 8)
    {
      DFT_K_gmp_ptr = DFT_8_aux_gmp;
      n_mult_pow_r_per_DFT_K = 5;
    }
  else if (K == 16)
    {
      DFT_K_gmp_ptr = DFT_16_aux_gmp;
      n_mult_pow_r_per_DFT_K = 17;
    }
  else if (K == 32)
    {
      DFT_K_gmp_ptr = DFT_32_aux_gmp;
      n_mult_pow_r_per_DFT_K = 66;
    }
  else if (K == 64)
    {
      DFT_K_gmp_ptr = DFT_64_aux_gmp;
      n_mult_pow_r_per_DFT_K = 195;
    }
  else if (K == 128)
    {
      DFT_K_gmp_ptr = DFT_128_aux_gmp;
      //ToDo: set the correct value for the following! -> DONE.
      n_mult_pow_r_per_DFT_K = 321;
    }
  else if (K == 256)
    {
      DFT_K_gmp_ptr = DFT_256_aux_gmp;
      //ToDo: set the correct value for the following! -> DONE.
      n_mult_pow_r_per_DFT_K = 769;
    }

  timer_record_start (&t_DFT_K);

#if GMP_SERIAL_ENABLED==1
  for (int j = 0; j < J; j++)
#else
#pragma cilk_grainsize = (J/(cilk::current_worker_count))
  cilk_for (int j = 0; j < J; j++)
#endif
    {
      DFT_K_gmp_ptr (&vector[j * K], aux_zz[j], omega_base, prime,
		     precomputed_pow_omega_base_vec);
    }
  timer_record_stop (&t_DFT_K);
  timer_get_elapsed_time (&t_DFT_K, NULL, 1);
  t_DFT_K_total += t_DFT_K.elapsed_time;
  n_DFT_K += J;
  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////
  ////step 3:
  ////- twiddle factor multiplication + permutation in between + DFT_K on the left
  ////- beginning at the level b=2;
  /////////////////////////////////////
  mpz_set (current_omega_zz, omega_zz);
  for (int b = 2; b <= e; b++)
    {
      J = pow (K, b - 1);
      int n_permutations = N / (K * J);

//      printf ("[step 3: (%-12s) b=%d, D_{K=%d, J=%llu}, n_permutations=%d]\n",
//	      b, K, J, n_permutations);
      usfixn64 stride = K * J;
      usfixn64 idx = 0;
      mpz_powm_ui (current_omega_zz, omega_zz, pow (K, e - b), prime);

      float t_twiddle_values[2] =
	{ 0, 0 };

      twiddle_aux_gmp (n_permutations, vector, aux_zz, K, J,
      		   precomputed_pow_omega_vec_2d[b],
      		   precomputed_pow_omega_base_vec, prime, t_twiddle_values,
      		   n_twiddle_mult);t_mult_pow_r_total += t_twiddle_values[0];
      t_arbitrary_mult_total += t_twiddle_values[1];

      timer_record_start (&t_permutation);
      stride_permutation_aux_gmp (n_permutations, vector, aux_zz, J, K);
      timer_record_stop (&t_permutation);
      timer_get_elapsed_time (&t_permutation, NULL, 1);
      t_stride_total += t_permutation.elapsed_time;

      /////////////////////////////////
      ////step 4:
      ////computing DFT_K on the left.
      /////////////////////////////////
      J = N / K;
      n_permutations = N / (K * J);
      // printf ("-- step 4: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
      // "left DFT_K", b, K, J, n_permutations);

      timer_record_start (&t_DFT_K);
  #if GMP_SERIAL_ENABLED==1
    for (int j = 0; j < J; j++)
  #else
  #pragma cilk_grainsize = (J/(cilk::current_worker_count))
    cilk_for (int j = 0; j < J; j++)
  #endif
	{
	  DFT_K_gmp_ptr (&vector[j * K], aux_zz[j], omega_base, prime,
			 precomputed_pow_omega_base_vec);
	}

      timer_record_stop (&t_DFT_K);
      timer_get_elapsed_time (&t_DFT_K, NULL, 1);
      t_DFT_K_total += t_DFT_K.elapsed_time;
      n_DFT_K += J;

      /////////////////////////////////
      //// step 5: left-most permutation L_{K,K^e-1} on the right
      /////////////////////////////////
//    printf("-- step 5: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
//    "left perm", b, K, J, n_permutations);

      J = pow (K, b - 1);
      n_permutations = N / (K * J);

      timer_record_start (&t_permutation);
      stride_permutation_aux_gmp (n_permutations, vector, aux_zz, K, J);
      timer_record_stop (&t_permutation);
      timer_get_elapsed_time (&t_permutation, NULL, 1);
      t_stride_total += t_permutation.elapsed_time;
    }

  timer_record_stop (&t_dft_gmp);
  timer_get_elapsed_time (&t_dft_gmp, NULL, 1);

  int n_mult_pow_r_in_DFT_K_total = (n_mult_pow_r_per_DFT_K * n_DFT_K);
  n_twiddle_mult[0]+=n_mult_pow_r_in_DFT_K_total;
  float t_mult_pow_r_avg;
  t_mult_pow_r_avg = t_mult_pow_r_total / n_twiddle_mult[0];

//  n_twiddle_mult[0] += n_mult_pow_r_in_DFT_K_total;
  float t_mult_pow_r_in_DFT_K_total = (t_mult_pow_r_avg
      * n_mult_pow_r_in_DFT_K_total);

//  t_mult_pow_r_total += t_mult_pow_r_in_DFT_K_total;
  t_twiddle_total = t_mult_pow_r_total + t_arbitrary_mult_total;

  if (verbose)
    {
      char suffix[64];
      char msg[256];

      sprintf (suffix, "%s", "gmp");

      sprintf (msg, "%s_%s", "t_precompute", suffix);
      timer_print_time_percentage (t_precompute_total, msg,
				   t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_stride_total", suffix);
      timer_print_time_percentage (t_stride_total, msg, t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_DFT_K_total", suffix);
      timer_print_time_percentage (t_DFT_K_total, msg, t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_DFT_K_total_avg", suffix);
      timer_print_time (t_DFT_K_total / n_DFT_K, msg);

      printf ("%s", quadline);

      sprintf (msg, "%s_%s", "t_mult_pow_r_in_DFT_K_total", suffix);
      timer_print_time_percentage (t_mult_pow_r_in_DFT_K_total, msg,
				   t_dft_gmp.elapsed_time);

      printf ("%s", quadline);

      sprintf (msg, "%s_%s", "t_mult_pow_r_total", suffix);
      timer_print_time_percentage (t_mult_pow_r_total, msg,
				   t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_mult_pow_r_total_avg", suffix);
      timer_print_time (t_mult_pow_r_total / n_twiddle_mult[0], msg);

      sprintf (msg, "%s_%s", "t_arbitrary_mult_total", suffix);
      timer_print_time_percentage (t_arbitrary_mult_total, msg,
				   t_dft_gmp.elapsed_time);

      sprintf (msg, "%s_%s", "t_arbitrary_mult_total_avg", suffix);
      timer_print_time (t_arbitrary_mult_total / n_twiddle_mult[1], msg);

      sprintf (msg, "%s_%s", "t_twiddle_total", suffix);
      timer_print_time_percentage (t_twiddle_total, msg,
				   t_dft_gmp.elapsed_time);
      printf ("%s", quadline);

      sprintf (msg, "%s_%s", "t_dft_six_step_aux", suffix);
      timer_print_time_percentage (t_dft_gmp.elapsed_time, msg,
				   t_dft_gmp.elapsed_time);
      printf ("%s", longline);
    }

  mpz_clear (current_omega_zz);
  mpz_clear (omega_base);
  for (int i = 0; i < K; i++)
    {
      mpz_clear (precomputed_pow_omega_base_vec[i]);
    }
  free (precomputed_pow_omega_base_vec);

  for (int b = e; b >= 2; b--)
    {
      J = pow (K, b - 1);
      for (int j = 0; j < J; j++)
	mpz_clear (precomputed_pow_omega_vec_2d[b][j]);
      free (precomputed_pow_omega_vec_2d[b]);
    }
  free (precomputed_pow_omega_vec_2d);

  for(usfixn64 i=0; i<N; i++)
    mpz_clear(aux_zz[i]);
  free(aux_zz);


}


/*************************************/
#endif
