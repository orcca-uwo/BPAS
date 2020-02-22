#ifndef GFPF_SIX_STEP_FFT_H_
#define GFPF_SIX_STEP_FFT_H_

/* gfpf_six_step_fft.h
 *****************************************
 * an implementation of big prime field fft which relies
 * entirely on radix-based representation of data.
 * Assuming that prime is of the form p=r^k+1, an input
 * vector X of size N can be viewed as a row-major layout
 * of an Nxk matrix.
 * An element X_i is stored as a vector of size k.
 * X_i= (a_0, a_1, ..., a_{k-1}) in base r
 * ie, X_i = a_0 + a_1 .r + ... + a_{k-1}.r^{k-1}.
 *****************************************
 */

#include "gfpf_six_step_fft_decl.h"

/**************************************/

static void __inline__
DFT2_big_elements (usfixn64* a0, usfixn64* a1, const int k, const usfixn64 r)
{

#if PROFILING_ENABLED == 1
  inc_profiling_counter (&n_dft2_called);
#endif

  //a0=a0+a1;
  addition_subtraction_big_elements (a0, a1, k, r);

}

/**************************************/

//older implementation of DFT2_big_elements
void
DFT2_big_elements_v0 (usfixn64* a0, usfixn64* a1, int k, usfixn64 r)
{

#if PROFILING_ENABLED == 1
  inc_profiling_counter (&n_dft2_called);
#endif

  usfixn64* sub = (usfixn64*) malloc (k * sizeof(usfixn64));
//  usfixn64 * sub = global_tmp_dft2;
  memcpy (sub, a0, k * sizeof(usfixn64));

  //a0=a0+a1;
  addition_big_elements_v0 (a0, a1, k, r);

  //sub=a0-a1;
  subtraction_big_elements_v0 (sub, a1, k, r);
  memcpy (a1, sub, k * sizeof(usfixn64));
  free (sub);
}

/**************************************/

static void __inline__
swap_big_elements (usfixn64* x, usfixn64* y, const int k)
{
#if PROFILING_ENABLED == 1
  inc_profiling_counter (&n_swap_called);
#endif

  for (int b = 0; b < k; b += SWAP_BLOCK_SIZE)
#pragma unroll SWAP_BLOCK_SIZE
    for (int i = b; i < b + SWAP_BLOCK_SIZE; i++)
      {
	__asm__ ("xchg %0, %1;\n\t": "+g"(x[i]),"+g"(y[i])::);
//	usfixn64 t = x[i];
//	x[i] = y[i];
//	y[i] = t;
//	      x[i]^=y[i];
//	      y[i]^=x[i];
//	      x[i]^=y[i];
      }
}

/**************************************/

//older implementation of swap_big_elements
static void __inline__
swap_big_elements_v0 (usfixn64* a, usfixn64* b, int k)
{
#if PROFILING_ENABLED == 1
  inc_profiling_counter (&n_swap_called);
#endif
  int n_bytes = k * sizeof(usfixn64);
  usfixn64* tmp = (usfixn64*) malloc (n_bytes);
//  usfixn64 * tmp=swap_memory;
//  usfixn64* tmp = tmp_swap_big_elements;
  memcpy (tmp, a, n_bytes);
  memcpy (a, b, n_bytes);
  memcpy (b, tmp, n_bytes);
  free (tmp);
}

/**************************************/

void
DFT_4_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		    int compute_inverse)
{
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 4];
      //[0,2][1,3]
      DFT2_big_elements (&a[0 * k], &a[2 * k], k, r);
      DFT2_big_elements (&a[1 * k], &a[3 * k], k, r);

      //x3.r^1
      mult_pow_R (&a[3 * k], 1, k, r, compute_inverse);

      //[0,1][2,3]
      DFT2_big_elements (&a[0 * k], &a[1 * k], k, r);
      DFT2_big_elements (&a[2 * k], &a[3 * k], k, r);

      //[0,2][1,3]
      //mapped to
      //[0,1][2,3]
      swap_big_elements (&a[1 * k], &a[2 * k], k);
    }
}

/**************************************/

void
DFT_8_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		    int compute_inverse)
{
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 8];
      //[0,4][2,6][1,5][3,7]
      DFT2_big_elements (&a[0 * k], &a[4 * k], k, r);
      DFT2_big_elements (&a[1 * k], &a[5 * k], k, r);
      DFT2_big_elements (&a[2 * k], &a[6 * k], k, r);
      DFT2_big_elements (&a[3 * k], &a[7 * k], k, r);

      //x6.r^2 , x7.r^2
      mult_pow_R (&a[6 * k], 2, k, r, compute_inverse);
      mult_pow_R (&a[7 * k], 2, k, r, compute_inverse);

      //[0,2][4,6][1,3][5,7]
      DFT2_big_elements (&a[0 * k], &a[2 * k], k, r);
      DFT2_big_elements (&a[1 * k], &a[3 * k], k, r);
      DFT2_big_elements (&a[4 * k], &a[6 * k], k, r);
      DFT2_big_elements (&a[5 * k], &a[7 * k], k, r);

      //[x5,x3,x7] x [r,r^2, r^3]
      mult_pow_R (&a[5 * k], 1, k, r, compute_inverse);
      mult_pow_R (&a[3 * k], 2, k, r, compute_inverse);
      mult_pow_R (&a[7 * k], 3, k, r, compute_inverse);

      //[0,1][4,5][2,3][6,7]
      DFT2_big_elements (&a[0 * k], &a[1 * k], k, r);
      DFT2_big_elements (&a[2 * k], &a[3 * k], k, r);
      DFT2_big_elements (&a[4 * k], &a[5 * k], k, r);
      DFT2_big_elements (&a[6 * k], &a[7 * k], k, r);

      //[0,4][2,6][1,5][3,7]
      //mapped to
      //[0,1][2,3][4,5][6,7]
      swap_big_elements (&a[1 * k], &a[4 * k], k);
      swap_big_elements (&a[3 * k], &a[6 * k], k);
    }
}

/**************************************/

void
DFT_16_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse)
{
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 16];

      DFT2_big_elements (&a[0 * k], &a[8 * k], k, r);
      DFT2_big_elements (&a[1 * k], &a[9 * k], k, r);
      DFT2_big_elements (&a[2 * k], &a[10 * k], k, r);
      DFT2_big_elements (&a[3 * k], &a[11 * k], k, r);
      DFT2_big_elements (&a[4 * k], &a[12 * k], k, r);
      DFT2_big_elements (&a[5 * k], &a[13 * k], k, r);
      DFT2_big_elements (&a[6 * k], &a[14 * k], k, r);
      DFT2_big_elements (&a[7 * k], &a[15 * k], k, r);

      //twiddle
      mult_pow_R (&a[12 * k], 4, k, r, compute_inverse);
      mult_pow_R (&a[14 * k], 4, k, r, compute_inverse);
      mult_pow_R (&a[13 * k], 4, k, r, compute_inverse);
      mult_pow_R (&a[15 * k], 4, k, r, compute_inverse);

      //dft on permutated indices
      DFT2_big_elements (&a[0 * k], &a[4 * k], k, r);
      DFT2_big_elements (&a[1 * k], &a[5 * k], k, r);
      DFT2_big_elements (&a[2 * k], &a[6 * k], k, r);
      DFT2_big_elements (&a[3 * k], &a[7 * k], k, r);
      DFT2_big_elements (&a[8 * k], &a[12 * k], k, r);
      DFT2_big_elements (&a[9 * k], &a[13 * k], k, r);
      DFT2_big_elements (&a[10 * k], &a[14 * k], k, r);
      DFT2_big_elements (&a[11 * k], &a[15 * k], k, r);

      //twiddle
      mult_pow_R (&a[6 * k], 4, k, r, compute_inverse);
      mult_pow_R (&a[7 * k], 4, k, r, compute_inverse);
      mult_pow_R (&a[10 * k], 2, k, r, compute_inverse);
      mult_pow_R (&a[11 * k], 2, k, r, compute_inverse);
      mult_pow_R (&a[14 * k], 6, k, r, compute_inverse);
      mult_pow_R (&a[15 * k], 6, k, r, compute_inverse);

      //dft on permutated indices.
      DFT2_big_elements (&a[0 * k], &a[2 * k], k, r);
      DFT2_big_elements (&a[1 * k], &a[3 * k], k, r);
      DFT2_big_elements (&a[4 * k], &a[6 * k], k, r);
      DFT2_big_elements (&a[5 * k], &a[7 * k], k, r);
      DFT2_big_elements (&a[8 * k], &a[10 * k], k, r);
      DFT2_big_elements (&a[9 * k], &a[11 * k], k, r);
      DFT2_big_elements (&a[12 * k], &a[14 * k], k, r);
      DFT2_big_elements (&a[13 * k], &a[15 * k], k, r);

      mult_pow_R (&a[3 * k], 4, k, r, compute_inverse);
      mult_pow_R (&a[5 * k], 2, k, r, compute_inverse);
      mult_pow_R (&a[7 * k], 6, k, r, compute_inverse);
      mult_pow_R (&a[9 * k], 1, k, r, compute_inverse);
      mult_pow_R (&a[11 * k], 5, k, r, compute_inverse);
      mult_pow_R (&a[13 * k], 3, k, r, compute_inverse);
      mult_pow_R (&a[15 * k], 7, k, r, compute_inverse);

      DFT2_big_elements (&a[0 * k], &a[1 * k], k, r);
      DFT2_big_elements (&a[2 * k], &a[3 * k], k, r);
      DFT2_big_elements (&a[4 * k], &a[5 * k], k, r);
      DFT2_big_elements (&a[6 * k], &a[7 * k], k, r);
      DFT2_big_elements (&a[8 * k], &a[9 * k], k, r);
      DFT2_big_elements (&a[10 * k], &a[11 * k], k, r);
      DFT2_big_elements (&a[12 * k], &a[13 * k], k, r);
      DFT2_big_elements (&a[14 * k], &a[15 * k], k, r);

      //final permutation.
      swap_big_elements (&a[1 * k], &a[8 * k], k);
      swap_big_elements (&a[2 * k], &a[4 * k], k);
      swap_big_elements (&a[3 * k], &a[12 * k], k);
      swap_big_elements (&a[5 * k], &a[10 * k], k);
      swap_big_elements (&a[7 * k], &a[14 * k], k);
      swap_big_elements (&a[11 * k], &a[13 * k], k);
    }
}

/**************************************/

void
DFT_32_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse)
{

//  int k_log2=log2(k);
  const int k_log2 = 4;
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 32];

      DFT2_big_elements (&a[0 << k_log2], &a[16 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[24 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[20 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[28 << k_log2], k, r);
      DFT2_big_elements (&a[2 << k_log2], &a[18 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[26 << k_log2], k, r);
      DFT2_big_elements (&a[6 << k_log2], &a[22 << k_log2], k, r);
      DFT2_big_elements (&a[14 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[1 << k_log2], &a[17 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[25 << k_log2], k, r);
      DFT2_big_elements (&a[5 << k_log2], &a[21 << k_log2], k, r);
      DFT2_big_elements (&a[13 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[3 << k_log2], &a[19 << k_log2], k, r);
      DFT2_big_elements (&a[11 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[7 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[15 << k_log2], &a[31 << k_log2], k, r);

      mult_pow_R (&a[24 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[28 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[26 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[30 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[25 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[29 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[27 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 8, k, r, compute_inverse);

      DFT2_big_elements (&a[0 << k_log2], &a[8 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[24 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[12 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[28 << k_log2], k, r);
      DFT2_big_elements (&a[2 << k_log2], &a[10 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[26 << k_log2], k, r);
      DFT2_big_elements (&a[6 << k_log2], &a[14 << k_log2], k, r);
      DFT2_big_elements (&a[22 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[1 << k_log2], &a[9 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[25 << k_log2], k, r);
      DFT2_big_elements (&a[5 << k_log2], &a[13 << k_log2], k, r);
      DFT2_big_elements (&a[21 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[3 << k_log2], &a[11 << k_log2], k, r);
      DFT2_big_elements (&a[19 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[7 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[23 << k_log2], &a[31 << k_log2], k, r);

      mult_pow_R (&a[20 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[12 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[28 << k_log2], 12, k, r, compute_inverse);

      mult_pow_R (&a[22 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[14 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[30 << k_log2], 12, k, r, compute_inverse);

      mult_pow_R (&a[21 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[13 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[29 << k_log2], 12, k, r, compute_inverse);

      mult_pow_R (&a[23 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[15 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 12, k, r, compute_inverse);

      DFT2_big_elements (&a[0 << k_log2], &a[4 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[20 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[12 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[28 << k_log2], k, r);
      DFT2_big_elements (&a[2 << k_log2], &a[6 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[22 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[14 << k_log2], k, r);
      DFT2_big_elements (&a[26 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[1 << k_log2], &a[5 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[21 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[13 << k_log2], k, r);
      DFT2_big_elements (&a[25 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[3 << k_log2], &a[7 << k_log2], k, r);
      DFT2_big_elements (&a[19 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[11 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[27 << k_log2], &a[31 << k_log2], k, r);

      mult_pow_R (&a[18 << k_log2], 2, k, r, compute_inverse);
      mult_pow_R (&a[10 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[26 << k_log2], 6, k, r, compute_inverse);
      mult_pow_R (&a[6 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[22 << k_log2], 10, k, r, compute_inverse);
      mult_pow_R (&a[14 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[30 << k_log2], 14, k, r, compute_inverse);

      mult_pow_R (&a[19 << k_log2], 2, k, r, compute_inverse);
      mult_pow_R (&a[11 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[27 << k_log2], 6, k, r, compute_inverse);
      mult_pow_R (&a[7 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[23 << k_log2], 10, k, r, compute_inverse);
      mult_pow_R (&a[15 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 14, k, r, compute_inverse);

      DFT2_big_elements (&a[0 << k_log2], &a[2 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[18 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[10 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[26 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[6 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[22 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[14 << k_log2], k, r);
      DFT2_big_elements (&a[28 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[1 << k_log2], &a[3 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[19 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[11 << k_log2], k, r);
      DFT2_big_elements (&a[25 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[5 << k_log2], &a[7 << k_log2], k, r);
      DFT2_big_elements (&a[21 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[13 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[29 << k_log2], &a[31 << k_log2], k, r);

      mult_pow_R (&a[17 << k_log2], 1, k, r, compute_inverse);
      mult_pow_R (&a[9 << k_log2], 2, k, r, compute_inverse);
      mult_pow_R (&a[25 << k_log2], 3, k, r, compute_inverse);
      mult_pow_R (&a[5 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[21 << k_log2], 5, k, r, compute_inverse);
      mult_pow_R (&a[13 << k_log2], 6, k, r, compute_inverse);
      mult_pow_R (&a[29 << k_log2], 7, k, r, compute_inverse);
      mult_pow_R (&a[3 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[19 << k_log2], 9, k, r, compute_inverse);
      mult_pow_R (&a[11 << k_log2], 10, k, r, compute_inverse);
      mult_pow_R (&a[27 << k_log2], 11, k, r, compute_inverse);
      mult_pow_R (&a[7 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[23 << k_log2], 13, k, r, compute_inverse);
      mult_pow_R (&a[15 << k_log2], 14, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 15, k, r, compute_inverse);

      DFT2_big_elements (&a[0 << k_log2], &a[1 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[17 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[9 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[25 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[5 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[21 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[13 << k_log2], k, r);
      DFT2_big_elements (&a[28 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[2 << k_log2], &a[3 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[19 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[11 << k_log2], k, r);
      DFT2_big_elements (&a[26 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[6 << k_log2], &a[7 << k_log2], k, r);
      DFT2_big_elements (&a[22 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[14 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[30 << k_log2], &a[31 << k_log2], k, r);

      swap_big_elements (&a[1 << k_log2], &a[16 << k_log2], k);
      swap_big_elements (&a[2 << k_log2], &a[8 << k_log2], k);
      swap_big_elements (&a[3 << k_log2], &a[24 << k_log2], k);
      swap_big_elements (&a[5 << k_log2], &a[20 << k_log2], k);
      swap_big_elements (&a[6 << k_log2], &a[12 << k_log2], k);
      swap_big_elements (&a[7 << k_log2], &a[28 << k_log2], k);
      swap_big_elements (&a[9 << k_log2], &a[18 << k_log2], k);
      swap_big_elements (&a[11 << k_log2], &a[26 << k_log2], k);
      swap_big_elements (&a[13 << k_log2], &a[22 << k_log2], k);
      swap_big_elements (&a[15 << k_log2], &a[30 << k_log2], k);
      swap_big_elements (&a[19 << k_log2], &a[25 << k_log2], k);
      swap_big_elements (&a[23 << k_log2], &a[29 << k_log2], k);
    }
}

/**************************************/

void
DFT_64_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse)
{

  const int k_log2 = 5;
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 64];

      DFT2_big_elements (&a[0 << k_log2], &a[32 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[48 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[40 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[56 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[36 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[52 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[44 << k_log2], k, r);
      DFT2_big_elements (&a[28 << k_log2], &a[60 << k_log2], k, r);
      DFT2_big_elements (&a[2 << k_log2], &a[34 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[50 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[42 << k_log2], k, r);
      DFT2_big_elements (&a[26 << k_log2], &a[58 << k_log2], k, r);
      DFT2_big_elements (&a[6 << k_log2], &a[38 << k_log2], k, r);
      DFT2_big_elements (&a[22 << k_log2], &a[54 << k_log2], k, r);
      DFT2_big_elements (&a[14 << k_log2], &a[46 << k_log2], k, r);
      DFT2_big_elements (&a[30 << k_log2], &a[62 << k_log2], k, r);

      DFT2_big_elements (&a[1 << k_log2], &a[33 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[49 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[41 << k_log2], k, r);
      DFT2_big_elements (&a[25 << k_log2], &a[57 << k_log2], k, r);
      DFT2_big_elements (&a[5 << k_log2], &a[37 << k_log2], k, r);
      DFT2_big_elements (&a[21 << k_log2], &a[53 << k_log2], k, r);
      DFT2_big_elements (&a[13 << k_log2], &a[45 << k_log2], k, r);
      DFT2_big_elements (&a[29 << k_log2], &a[61 << k_log2], k, r);
      DFT2_big_elements (&a[3 << k_log2], &a[35 << k_log2], k, r);
      DFT2_big_elements (&a[19 << k_log2], &a[51 << k_log2], k, r);
      DFT2_big_elements (&a[11 << k_log2], &a[43 << k_log2], k, r);
      DFT2_big_elements (&a[27 << k_log2], &a[59 << k_log2], k, r);
      DFT2_big_elements (&a[7 << k_log2], &a[39 << k_log2], k, r);
      DFT2_big_elements (&a[23 << k_log2], &a[55 << k_log2], k, r);
      DFT2_big_elements (&a[15 << k_log2], &a[47 << k_log2], k, r);
      DFT2_big_elements (&a[31 << k_log2], &a[63 << k_log2], k, r);

      //T_{2}^{4}
      mult_pow_R (&a[48 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[56 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[52 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[60 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[50 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[58 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[54 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[62 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[49 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[57 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[53 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[61 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[51 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[59 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[55 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[63 << k_log2], 16, k, r, compute_inverse);

      //I_{32} \tx DFT2
      DFT2_big_elements (&a[0 << k_log2], &a[16 << k_log2], k, r);
      DFT2_big_elements (&a[32 << k_log2], &a[48 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[24 << k_log2], k, r);
      DFT2_big_elements (&a[40 << k_log2], &a[56 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[20 << k_log2], k, r);
      DFT2_big_elements (&a[36 << k_log2], &a[52 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[28 << k_log2], k, r);
      DFT2_big_elements (&a[44 << k_log2], &a[60 << k_log2], k, r);
      DFT2_big_elements (&a[2 << k_log2], &a[18 << k_log2], k, r);
      DFT2_big_elements (&a[34 << k_log2], &a[50 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[26 << k_log2], k, r);
      DFT2_big_elements (&a[42 << k_log2], &a[58 << k_log2], k, r);
      DFT2_big_elements (&a[6 << k_log2], &a[22 << k_log2], k, r);
      DFT2_big_elements (&a[38 << k_log2], &a[54 << k_log2], k, r);
      DFT2_big_elements (&a[14 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[46 << k_log2], &a[62 << k_log2], k, r);
      DFT2_big_elements (&a[1 << k_log2], &a[17 << k_log2], k, r);
      DFT2_big_elements (&a[33 << k_log2], &a[49 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[25 << k_log2], k, r);
      DFT2_big_elements (&a[41 << k_log2], &a[57 << k_log2], k, r);
      DFT2_big_elements (&a[5 << k_log2], &a[21 << k_log2], k, r);
      DFT2_big_elements (&a[37 << k_log2], &a[53 << k_log2], k, r);
      DFT2_big_elements (&a[13 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[45 << k_log2], &a[61 << k_log2], k, r);
      DFT2_big_elements (&a[3 << k_log2], &a[19 << k_log2], k, r);
      DFT2_big_elements (&a[35 << k_log2], &a[51 << k_log2], k, r);
      DFT2_big_elements (&a[11 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[43 << k_log2], &a[59 << k_log2], k, r);
      DFT2_big_elements (&a[7 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[39 << k_log2], &a[55 << k_log2], k, r);
      DFT2_big_elements (&a[15 << k_log2], &a[31 << k_log2], k, r);
      DFT2_big_elements (&a[47 << k_log2], &a[63 << k_log2], k, r);

      // I_{8} \tx T_{4}^{8}

      mult_pow_R (&a[40 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[24 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[56 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[44 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[28 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[60 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[42 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[26 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[58 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[46 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[30 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[62 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[41 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[25 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[57 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[45 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[29 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[61 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[43 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[27 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[59 << k_log2], 24, k, r, compute_inverse);

      mult_pow_R (&a[47 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[63 << k_log2], 24, k, r, compute_inverse);

      // I_{32} \tx DFT2
      DFT2_big_elements (&a[0 << k_log2], &a[8 << k_log2], k, r);
      DFT2_big_elements (&a[32 << k_log2], &a[40 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[24 << k_log2], k, r);
      DFT2_big_elements (&a[48 << k_log2], &a[56 << k_log2], k, r);
      DFT2_big_elements (&a[4 << k_log2], &a[12 << k_log2], k, r);
      DFT2_big_elements (&a[36 << k_log2], &a[44 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[28 << k_log2], k, r);
      DFT2_big_elements (&a[52 << k_log2], &a[60 << k_log2], k, r);

      DFT2_big_elements (&a[2 << k_log2], &a[10 << k_log2], k, r);
      DFT2_big_elements (&a[34 << k_log2], &a[42 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[26 << k_log2], k, r);
      DFT2_big_elements (&a[50 << k_log2], &a[58 << k_log2], k, r);
      DFT2_big_elements (&a[6 << k_log2], &a[14 << k_log2], k, r);
      DFT2_big_elements (&a[38 << k_log2], &a[46 << k_log2], k, r);
      DFT2_big_elements (&a[22 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[54 << k_log2], &a[62 << k_log2], k, r);

      DFT2_big_elements (&a[1 << k_log2], &a[9 << k_log2], k, r);
      DFT2_big_elements (&a[33 << k_log2], &a[41 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[25 << k_log2], k, r);
      DFT2_big_elements (&a[49 << k_log2], &a[57 << k_log2], k, r);
      DFT2_big_elements (&a[5 << k_log2], &a[13 << k_log2], k, r);
      DFT2_big_elements (&a[37 << k_log2], &a[45 << k_log2], k, r);
      DFT2_big_elements (&a[21 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[53 << k_log2], &a[61 << k_log2], k, r);

      DFT2_big_elements (&a[3 << k_log2], &a[11 << k_log2], k, r);
      DFT2_big_elements (&a[35 << k_log2], &a[43 << k_log2], k, r);
      DFT2_big_elements (&a[19 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[51 << k_log2], &a[59 << k_log2], k, r);
      DFT2_big_elements (&a[7 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[39 << k_log2], &a[47 << k_log2], k, r);
      DFT2_big_elements (&a[23 << k_log2], &a[31 << k_log2], k, r);
      DFT2_big_elements (&a[55 << k_log2], &a[63 << k_log2], k, r);

      // T_{8}^{16}

      mult_pow_R (&a[36 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[20 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[52 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[12 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[44 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[28 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[60 << k_log2], 28, k, r, compute_inverse);

      mult_pow_R (&a[38 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[22 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[54 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[14 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[46 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[30 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[62 << k_log2], 28, k, r, compute_inverse);

      mult_pow_R (&a[37 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[21 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[53 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[13 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[45 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[29 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[61 << k_log2], 28, k, r, compute_inverse);

      mult_pow_R (&a[39 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[23 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[55 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[15 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[47 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[63 << k_log2], 28, k, r, compute_inverse);

      // I_{32} \ts DFT2
      DFT2_big_elements (&a[0 << k_log2], &a[4 << k_log2], k, r);
      DFT2_big_elements (&a[32 << k_log2], &a[36 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[20 << k_log2], k, r);
      DFT2_big_elements (&a[48 << k_log2], &a[52 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[12 << k_log2], k, r);
      DFT2_big_elements (&a[40 << k_log2], &a[44 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[28 << k_log2], k, r);
      DFT2_big_elements (&a[56 << k_log2], &a[60 << k_log2], k, r);

      DFT2_big_elements (&a[2 << k_log2], &a[6 << k_log2], k, r);
      DFT2_big_elements (&a[34 << k_log2], &a[38 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[22 << k_log2], k, r);
      DFT2_big_elements (&a[50 << k_log2], &a[54 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[14 << k_log2], k, r);
      DFT2_big_elements (&a[42 << k_log2], &a[46 << k_log2], k, r);
      DFT2_big_elements (&a[26 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[58 << k_log2], &a[62 << k_log2], k, r);

      DFT2_big_elements (&a[1 << k_log2], &a[5 << k_log2], k, r);
      DFT2_big_elements (&a[33 << k_log2], &a[37 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[21 << k_log2], k, r);
      DFT2_big_elements (&a[49 << k_log2], &a[53 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[13 << k_log2], k, r);
      DFT2_big_elements (&a[41 << k_log2], &a[45 << k_log2], k, r);
      DFT2_big_elements (&a[25 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[57 << k_log2], &a[61 << k_log2], k, r);

      DFT2_big_elements (&a[3 << k_log2], &a[7 << k_log2], k, r);
      DFT2_big_elements (&a[35 << k_log2], &a[39 << k_log2], k, r);
      DFT2_big_elements (&a[19 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[51 << k_log2], &a[55 << k_log2], k, r);
      DFT2_big_elements (&a[11 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[43 << k_log2], &a[47 << k_log2], k, r);
      DFT2_big_elements (&a[27 << k_log2], &a[31 << k_log2], k, r);
      DFT2_big_elements (&a[59 << k_log2], &a[63 << k_log2], k, r);

      // I_{2} \tx T_{16}^{32}

      mult_pow_R (&a[34 << k_log2], 2, k, r, compute_inverse);
      mult_pow_R (&a[18 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[50 << k_log2], 6, k, r, compute_inverse);
      mult_pow_R (&a[10 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[42 << k_log2], 10, k, r, compute_inverse);
      mult_pow_R (&a[26 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[58 << k_log2], 14, k, r, compute_inverse);
      mult_pow_R (&a[6 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[38 << k_log2], 18, k, r, compute_inverse);
      mult_pow_R (&a[22 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[54 << k_log2], 22, k, r, compute_inverse);
      mult_pow_R (&a[14 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[46 << k_log2], 26, k, r, compute_inverse);
      mult_pow_R (&a[30 << k_log2], 28, k, r, compute_inverse);
      mult_pow_R (&a[62 << k_log2], 30, k, r, compute_inverse);

      mult_pow_R (&a[35 << k_log2], 2, k, r, compute_inverse);
      mult_pow_R (&a[19 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[51 << k_log2], 6, k, r, compute_inverse);
      mult_pow_R (&a[11 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[43 << k_log2], 10, k, r, compute_inverse);
      mult_pow_R (&a[27 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[59 << k_log2], 14, k, r, compute_inverse);
      mult_pow_R (&a[7 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[39 << k_log2], 18, k, r, compute_inverse);
      mult_pow_R (&a[23 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[55 << k_log2], 22, k, r, compute_inverse);
      mult_pow_R (&a[15 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[47 << k_log2], 26, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 28, k, r, compute_inverse);
      mult_pow_R (&a[63 << k_log2], 30, k, r, compute_inverse);

      // I_{2} I_{16} \tx DFT2
      DFT2_big_elements (&a[0 << k_log2], &a[2 << k_log2], k, r);
      DFT2_big_elements (&a[32 << k_log2], &a[34 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[18 << k_log2], k, r);
      DFT2_big_elements (&a[48 << k_log2], &a[50 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[10 << k_log2], k, r);
      DFT2_big_elements (&a[40 << k_log2], &a[42 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[26 << k_log2], k, r);
      DFT2_big_elements (&a[56 << k_log2], &a[58 << k_log2], k, r);

      DFT2_big_elements (&a[4 << k_log2], &a[6 << k_log2], k, r);
      DFT2_big_elements (&a[36 << k_log2], &a[38 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[22 << k_log2], k, r);
      DFT2_big_elements (&a[52 << k_log2], &a[54 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[14 << k_log2], k, r);
      DFT2_big_elements (&a[44 << k_log2], &a[46 << k_log2], k, r);
      DFT2_big_elements (&a[28 << k_log2], &a[30 << k_log2], k, r);
      DFT2_big_elements (&a[60 << k_log2], &a[62 << k_log2], k, r);

      DFT2_big_elements (&a[1 << k_log2], &a[3 << k_log2], k, r);
      DFT2_big_elements (&a[33 << k_log2], &a[35 << k_log2], k, r);
      DFT2_big_elements (&a[17 << k_log2], &a[19 << k_log2], k, r);
      DFT2_big_elements (&a[49 << k_log2], &a[51 << k_log2], k, r);
      DFT2_big_elements (&a[9 << k_log2], &a[11 << k_log2], k, r);
      DFT2_big_elements (&a[41 << k_log2], &a[43 << k_log2], k, r);
      DFT2_big_elements (&a[25 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[57 << k_log2], &a[59 << k_log2], k, r);

      DFT2_big_elements (&a[5 << k_log2], &a[7 << k_log2], k, r);
      DFT2_big_elements (&a[37 << k_log2], &a[39 << k_log2], k, r);
      DFT2_big_elements (&a[21 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[53 << k_log2], &a[55 << k_log2], k, r);
      DFT2_big_elements (&a[13 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[45 << k_log2], &a[47 << k_log2], k, r);
      DFT2_big_elements (&a[29 << k_log2], &a[31 << k_log2], k, r);
      DFT2_big_elements (&a[61 << k_log2], &a[63 << k_log2], k, r);

      // T_{32}^{64}

      mult_pow_R (&a[33 << k_log2], 1, k, r, compute_inverse);
      mult_pow_R (&a[17 << k_log2], 2, k, r, compute_inverse);
      mult_pow_R (&a[49 << k_log2], 3, k, r, compute_inverse);
      mult_pow_R (&a[9 << k_log2], 4, k, r, compute_inverse);
      mult_pow_R (&a[41 << k_log2], 5, k, r, compute_inverse);
      mult_pow_R (&a[25 << k_log2], 6, k, r, compute_inverse);
      mult_pow_R (&a[57 << k_log2], 7, k, r, compute_inverse);
      mult_pow_R (&a[5 << k_log2], 8, k, r, compute_inverse);
      mult_pow_R (&a[37 << k_log2], 9, k, r, compute_inverse);
      mult_pow_R (&a[21 << k_log2], 10, k, r, compute_inverse);
      mult_pow_R (&a[53 << k_log2], 11, k, r, compute_inverse);
      mult_pow_R (&a[13 << k_log2], 12, k, r, compute_inverse);
      mult_pow_R (&a[45 << k_log2], 13, k, r, compute_inverse);
      mult_pow_R (&a[29 << k_log2], 14, k, r, compute_inverse);
      mult_pow_R (&a[61 << k_log2], 15, k, r, compute_inverse);
      mult_pow_R (&a[3 << k_log2], 16, k, r, compute_inverse);
      mult_pow_R (&a[35 << k_log2], 17, k, r, compute_inverse);
      mult_pow_R (&a[19 << k_log2], 18, k, r, compute_inverse);
      mult_pow_R (&a[51 << k_log2], 19, k, r, compute_inverse);
      mult_pow_R (&a[11 << k_log2], 20, k, r, compute_inverse);
      mult_pow_R (&a[43 << k_log2], 21, k, r, compute_inverse);
      mult_pow_R (&a[27 << k_log2], 22, k, r, compute_inverse);
      mult_pow_R (&a[59 << k_log2], 23, k, r, compute_inverse);
      mult_pow_R (&a[7 << k_log2], 24, k, r, compute_inverse);
      mult_pow_R (&a[39 << k_log2], 25, k, r, compute_inverse);
      mult_pow_R (&a[23 << k_log2], 26, k, r, compute_inverse);
      mult_pow_R (&a[55 << k_log2], 27, k, r, compute_inverse);
      mult_pow_R (&a[15 << k_log2], 28, k, r, compute_inverse);
      mult_pow_R (&a[47 << k_log2], 29, k, r, compute_inverse);
      mult_pow_R (&a[31 << k_log2], 30, k, r, compute_inverse);
      mult_pow_R (&a[63 << k_log2], 31, k, r, compute_inverse);

      // I_{32} \tx DFT2
      DFT2_big_elements (&a[0 << k_log2], &a[1 << k_log2], k, r);
      DFT2_big_elements (&a[32 << k_log2], &a[33 << k_log2], k, r);
      DFT2_big_elements (&a[16 << k_log2], &a[17 << k_log2], k, r);
      DFT2_big_elements (&a[48 << k_log2], &a[49 << k_log2], k, r);
      DFT2_big_elements (&a[8 << k_log2], &a[9 << k_log2], k, r);
      DFT2_big_elements (&a[40 << k_log2], &a[41 << k_log2], k, r);
      DFT2_big_elements (&a[24 << k_log2], &a[25 << k_log2], k, r);
      DFT2_big_elements (&a[56 << k_log2], &a[57 << k_log2], k, r);

      DFT2_big_elements (&a[4 << k_log2], &a[5 << k_log2], k, r);
      DFT2_big_elements (&a[36 << k_log2], &a[37 << k_log2], k, r);
      DFT2_big_elements (&a[20 << k_log2], &a[21 << k_log2], k, r);
      DFT2_big_elements (&a[52 << k_log2], &a[53 << k_log2], k, r);
      DFT2_big_elements (&a[12 << k_log2], &a[13 << k_log2], k, r);
      DFT2_big_elements (&a[44 << k_log2], &a[45 << k_log2], k, r);
      DFT2_big_elements (&a[28 << k_log2], &a[29 << k_log2], k, r);
      DFT2_big_elements (&a[60 << k_log2], &a[61 << k_log2], k, r);

      DFT2_big_elements (&a[2 << k_log2], &a[3 << k_log2], k, r);
      DFT2_big_elements (&a[34 << k_log2], &a[35 << k_log2], k, r);
      DFT2_big_elements (&a[18 << k_log2], &a[19 << k_log2], k, r);
      DFT2_big_elements (&a[50 << k_log2], &a[51 << k_log2], k, r);
      DFT2_big_elements (&a[10 << k_log2], &a[11 << k_log2], k, r);
      DFT2_big_elements (&a[42 << k_log2], &a[43 << k_log2], k, r);
      DFT2_big_elements (&a[26 << k_log2], &a[27 << k_log2], k, r);
      DFT2_big_elements (&a[58 << k_log2], &a[59 << k_log2], k, r);

      DFT2_big_elements (&a[6 << k_log2], &a[7 << k_log2], k, r);
      DFT2_big_elements (&a[38 << k_log2], &a[39 << k_log2], k, r);
      DFT2_big_elements (&a[22 << k_log2], &a[23 << k_log2], k, r);
      DFT2_big_elements (&a[54 << k_log2], &a[55 << k_log2], k, r);
      DFT2_big_elements (&a[14 << k_log2], &a[15 << k_log2], k, r);
      DFT2_big_elements (&a[46 << k_log2], &a[47 << k_log2], k, r);
      DFT2_big_elements (&a[30 << k_log2], &a[31 << k_log2], k, r);
      DFT2_big_elements (&a[62 << k_log2], &a[63 << k_log2], k, r);

      // Final permutation
      swap_big_elements (&a[1 << k_log2], &a[32 << k_log2], k);
      swap_big_elements (&a[2 << k_log2], &a[16 << k_log2], k);
      swap_big_elements (&a[3 << k_log2], &a[48 << k_log2], k);
      swap_big_elements (&a[4 << k_log2], &a[8 << k_log2], k);
      swap_big_elements (&a[5 << k_log2], &a[40 << k_log2], k);
      swap_big_elements (&a[6 << k_log2], &a[24 << k_log2], k);
      swap_big_elements (&a[7 << k_log2], &a[56 << k_log2], k);
      swap_big_elements (&a[9 << k_log2], &a[36 << k_log2], k);
      swap_big_elements (&a[10 << k_log2], &a[20 << k_log2], k);
      swap_big_elements (&a[11 << k_log2], &a[52 << k_log2], k);
      swap_big_elements (&a[13 << k_log2], &a[44 << k_log2], k);
      swap_big_elements (&a[14 << k_log2], &a[28 << k_log2], k);
      swap_big_elements (&a[15 << k_log2], &a[60 << k_log2], k);
      swap_big_elements (&a[17 << k_log2], &a[34 << k_log2], k);
      swap_big_elements (&a[19 << k_log2], &a[50 << k_log2], k);
      swap_big_elements (&a[21 << k_log2], &a[42 << k_log2], k);
      swap_big_elements (&a[22 << k_log2], &a[26 << k_log2], k);
      swap_big_elements (&a[23 << k_log2], &a[58 << k_log2], k);
      swap_big_elements (&a[25 << k_log2], &a[38 << k_log2], k);
      swap_big_elements (&a[27 << k_log2], &a[54 << k_log2], k);
      swap_big_elements (&a[29 << k_log2], &a[46 << k_log2], k);
      swap_big_elements (&a[31 << k_log2], &a[62 << k_log2], k);
      swap_big_elements (&a[35 << k_log2], &a[49 << k_log2], k);
      swap_big_elements (&a[37 << k_log2], &a[41 << k_log2], k);
      swap_big_elements (&a[39 << k_log2], &a[57 << k_log2], k);
      swap_big_elements (&a[43 << k_log2], &a[53 << k_log2], k);
      swap_big_elements (&a[47 << k_log2], &a[61 << k_log2], k);
      swap_big_elements (&a[55 << k_log2], &a[59 << k_log2], k);
    }
}

/**************************************/

void
DFT_128_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse)
{

  const int k_log2 = 5;
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 128];

      DFT2_big_elements(&a[0*k], &a[64*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[65*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[66*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[67*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[68*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[69*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[70*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[72*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[73*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[74*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[76*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[15*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[80*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[81*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[82*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[84*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[23*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[88*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[27*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[29*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[30*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[31*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[96*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[97*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[98*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[100*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[39*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[104*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[43*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[45*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[46*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[47*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[112*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[51*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[53*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[54*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[55*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[57*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[58*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[59*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[60*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[61*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[62*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[63*k], &a[127*k], k,	r);
      mult_pow_R(&a[96*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[97*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[98*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[99*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[100*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[101*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[102*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[104*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[105*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[106*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[108*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[112*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[113*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[114*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[116*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[120*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 32, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[32*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[33*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[34*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[35*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[36*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[37*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[38*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[40*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[41*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[42*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[44*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[15*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[48*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[49*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[50*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[52*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[23*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[56*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[27*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[29*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[30*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[31*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[96*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[97*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[98*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[100*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[71*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[104*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[75*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[77*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[78*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[79*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[112*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[83*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[85*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[86*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[87*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[89*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[90*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[91*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[92*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[93*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[94*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[95*k], &a[127*k], k,	r);
      mult_pow_R(&a[48*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[49*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[50*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[51*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[52*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[53*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[54*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[56*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[57*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[58*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[60*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[80*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[81*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[82*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[83*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[84*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[85*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[86*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[88*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[89*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[90*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[92*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[112*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[113*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[114*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[116*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[120*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 48, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[16*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[17*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[18*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[19*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[20*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[21*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[22*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[24*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[25*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[26*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[28*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[15*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[48*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[49*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[50*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[52*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[39*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[56*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[43*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[45*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[46*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[47*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[80*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[81*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[82*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[84*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[71*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[88*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[75*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[77*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[78*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[79*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[112*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[99*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[101*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[102*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[103*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[105*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[106*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[107*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[108*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[109*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[110*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[111*k], &a[127*k], k,	r);
      mult_pow_R(&a[24*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[25*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[26*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[27*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[28*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[29*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[30*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[40*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[41*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[42*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[43*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[44*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[45*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[46*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[56*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[57*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[58*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[60*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[72*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[73*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[74*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[75*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[76*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[77*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[78*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[88*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[89*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[90*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[92*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[104*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[105*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[106*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[108*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[120*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 56, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[8*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[9*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[10*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[11*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[12*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[13*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[14*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[24*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[25*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[26*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[28*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[23*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[40*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[41*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[42*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[44*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[39*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[56*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[51*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[53*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[54*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[55*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[72*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[73*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[74*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[76*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[71*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[88*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[83*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[85*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[86*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[87*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[104*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[99*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[101*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[102*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[103*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[113*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[114*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[115*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[116*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[117*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[118*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[119*k], &a[127*k], k,	r);
      mult_pow_R(&a[12*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[13*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[14*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[15*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[20*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[21*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[22*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[23*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[28*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[29*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[30*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[36*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[37*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[38*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[39*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[44*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[45*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[46*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[52*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[53*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[54*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[60*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[68*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[69*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[70*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[71*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[76*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[77*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[78*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[84*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[85*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[86*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[92*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[100*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[101*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[102*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[108*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[116*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 60, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[4*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[5*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[6*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[7*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[12*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[13*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[14*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[20*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[21*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[22*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[28*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[27*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[36*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[37*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[38*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[44*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[43*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[52*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[51*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[57*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[58*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[59*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[68*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[69*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[70*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[76*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[75*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[84*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[83*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[89*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[90*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[91*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[100*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[99*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[105*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[106*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[107*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[113*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[114*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[115*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[120*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[121*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[122*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[123*k], &a[127*k], k,	r);
      mult_pow_R(&a[6*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[7*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[10*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[11*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[14*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[15*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[18*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[19*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[22*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[23*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[26*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[27*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[30*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[34*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[35*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[38*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[39*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[42*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[43*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[46*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[50*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[51*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[54*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[58*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[66*k], 2, k, r, compute_inverse);
      mult_pow_R(&a[67*k], 2, k, r, compute_inverse);
      mult_pow_R(&a[70*k], 34, k, r, compute_inverse);
      mult_pow_R(&a[71*k], 34, k, r, compute_inverse);
      mult_pow_R(&a[74*k], 18, k, r, compute_inverse);
      mult_pow_R(&a[75*k], 18, k, r, compute_inverse);
      mult_pow_R(&a[78*k], 50, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 50, k, r, compute_inverse);
      mult_pow_R(&a[82*k], 10, k, r, compute_inverse);
      mult_pow_R(&a[83*k], 10, k, r, compute_inverse);
      mult_pow_R(&a[86*k], 42, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 42, k, r, compute_inverse);
      mult_pow_R(&a[90*k], 26, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 26, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 58, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 58, k, r, compute_inverse);
      mult_pow_R(&a[98*k], 6, k, r, compute_inverse);
      mult_pow_R(&a[99*k], 6, k, r, compute_inverse);
      mult_pow_R(&a[102*k], 38, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 38, k, r, compute_inverse);
      mult_pow_R(&a[106*k], 22, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 22, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 54, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 54, k, r, compute_inverse);
      mult_pow_R(&a[114*k], 14, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 14, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 46, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 46, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 30, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 30, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 62, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 62, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[2*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[3*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[6*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[7*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[10*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[11*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[14*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[18*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[19*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[22*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[26*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[29*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[34*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[35*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[38*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[42*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[45*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[50*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[53*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[57*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[60*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[61*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[66*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[67*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[70*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[74*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[77*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[82*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[85*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[89*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[92*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[93*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[98*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[101*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[105*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[108*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[109*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[113*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[116*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[117*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[120*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[121*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[124*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[125*k], &a[127*k], k,	r);
      mult_pow_R(&a[3*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[5*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[7*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[9*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[11*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[13*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[15*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[17*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[19*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[21*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[23*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[25*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[27*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[29*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[33*k], 2, k, r, compute_inverse);
      mult_pow_R(&a[35*k], 34, k, r, compute_inverse);
      mult_pow_R(&a[37*k], 18, k, r, compute_inverse);
      mult_pow_R(&a[39*k], 50, k, r, compute_inverse);
      mult_pow_R(&a[41*k], 10, k, r, compute_inverse);
      mult_pow_R(&a[43*k], 42, k, r, compute_inverse);
      mult_pow_R(&a[45*k], 26, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 58, k, r, compute_inverse);
      mult_pow_R(&a[49*k], 6, k, r, compute_inverse);
      mult_pow_R(&a[51*k], 38, k, r, compute_inverse);
      mult_pow_R(&a[53*k], 22, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 54, k, r, compute_inverse);
      mult_pow_R(&a[57*k], 14, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 46, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 30, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 62, k, r, compute_inverse);
      mult_pow_R(&a[65*k], 1, k, r, compute_inverse);
      mult_pow_R(&a[67*k], 33, k, r, compute_inverse);
      mult_pow_R(&a[69*k], 17, k, r, compute_inverse);
      mult_pow_R(&a[71*k], 49, k, r, compute_inverse);
      mult_pow_R(&a[73*k], 9, k, r, compute_inverse);
      mult_pow_R(&a[75*k], 41, k, r, compute_inverse);
      mult_pow_R(&a[77*k], 25, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 57, k, r, compute_inverse);
      mult_pow_R(&a[81*k], 5, k, r, compute_inverse);
      mult_pow_R(&a[83*k], 37, k, r, compute_inverse);
      mult_pow_R(&a[85*k], 21, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 53, k, r, compute_inverse);
      mult_pow_R(&a[89*k], 13, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 45, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 29, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 61, k, r, compute_inverse);
      mult_pow_R(&a[97*k], 3, k, r, compute_inverse);
      mult_pow_R(&a[99*k], 35, k, r, compute_inverse);
      mult_pow_R(&a[101*k], 19, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 51, k, r, compute_inverse);
      mult_pow_R(&a[105*k], 11, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 43, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 27, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 59, k, r, compute_inverse);
      mult_pow_R(&a[113*k], 7, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 39, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 23, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 55, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 15, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 47, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 31, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 63, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[1*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[3*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[5*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[7*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[9*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[11*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[13*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[17*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[19*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[21*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[25*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[30*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[33*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[35*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[37*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[41*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[46*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[49*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[54*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[58*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[60*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[62*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[65*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[67*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[69*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[73*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[78*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[81*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[86*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[90*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[92*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[94*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[97*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[102*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[106*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[108*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[110*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[114*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[116*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[118*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[120*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[122*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[124*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[126*k], &a[127*k], k,	r);
      swap_big_elements (&a[1*k], &a[64*k], k);
      swap_big_elements (&a[2*k], &a[32*k], k);
      swap_big_elements (&a[3*k], &a[96*k], k);
      swap_big_elements (&a[4*k], &a[16*k], k);
      swap_big_elements (&a[5*k], &a[80*k], k);
      swap_big_elements (&a[6*k], &a[48*k], k);
      swap_big_elements (&a[7*k], &a[112*k], k);
      swap_big_elements (&a[9*k], &a[72*k], k);
      swap_big_elements (&a[10*k], &a[40*k], k);
      swap_big_elements (&a[11*k], &a[104*k], k);
      swap_big_elements (&a[12*k], &a[24*k], k);
      swap_big_elements (&a[13*k], &a[88*k], k);
      swap_big_elements (&a[14*k], &a[56*k], k);
      swap_big_elements (&a[15*k], &a[120*k], k);
      swap_big_elements (&a[17*k], &a[68*k], k);
      swap_big_elements (&a[18*k], &a[36*k], k);
      swap_big_elements (&a[19*k], &a[100*k], k);
      swap_big_elements (&a[21*k], &a[84*k], k);
      swap_big_elements (&a[22*k], &a[52*k], k);
      swap_big_elements (&a[23*k], &a[116*k], k);
      swap_big_elements (&a[25*k], &a[76*k], k);
      swap_big_elements (&a[26*k], &a[44*k], k);
      swap_big_elements (&a[27*k], &a[108*k], k);
      swap_big_elements (&a[29*k], &a[92*k], k);
      swap_big_elements (&a[30*k], &a[60*k], k);
      swap_big_elements (&a[31*k], &a[124*k], k);
      swap_big_elements (&a[33*k], &a[66*k], k);
      swap_big_elements (&a[35*k], &a[98*k], k);
      swap_big_elements (&a[37*k], &a[82*k], k);
      swap_big_elements (&a[38*k], &a[50*k], k);
      swap_big_elements (&a[39*k], &a[114*k], k);
      swap_big_elements (&a[41*k], &a[74*k], k);
      swap_big_elements (&a[43*k], &a[106*k], k);
      swap_big_elements (&a[45*k], &a[90*k], k);
      swap_big_elements (&a[46*k], &a[58*k], k);
      swap_big_elements (&a[47*k], &a[122*k], k);
      swap_big_elements (&a[49*k], &a[70*k], k);
      swap_big_elements (&a[51*k], &a[102*k], k);
      swap_big_elements (&a[53*k], &a[86*k], k);
      swap_big_elements (&a[55*k], &a[118*k], k);
      swap_big_elements (&a[57*k], &a[78*k], k);
      swap_big_elements (&a[59*k], &a[110*k], k);
      swap_big_elements (&a[61*k], &a[94*k], k);
      swap_big_elements (&a[63*k], &a[126*k], k);
      swap_big_elements (&a[67*k], &a[97*k], k);
      swap_big_elements (&a[69*k], &a[81*k], k);
      swap_big_elements (&a[71*k], &a[113*k], k);
      swap_big_elements (&a[75*k], &a[105*k], k);
      swap_big_elements (&a[77*k], &a[89*k], k);
      swap_big_elements (&a[79*k], &a[121*k], k);
      swap_big_elements (&a[83*k], &a[101*k], k);
      swap_big_elements (&a[87*k], &a[117*k], k);
      swap_big_elements (&a[91*k], &a[109*k], k);
      swap_big_elements (&a[95*k], &a[125*k], k);
      swap_big_elements (&a[103*k], &a[115*k], k);
      swap_big_elements (&a[111*k], &a[123*k], k);

    }
}

/**************************************/

void
DFT_256_big_elements (int n, usfixn64* x, const int k, const usfixn64 r,
		     int compute_inverse)
{
  const int k_log2 = 5;
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  cilk_for (int i = 0; i < n; i++)
    {
#if PROFILING_ENABLED == 1
      inc_profiling_counter (&n_DFTK_called);
#endif

      usfixn64 *a = &x[i * k * 256];

      DFT2_big_elements(&a[0*k], &a[64*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[65*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[66*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[67*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[68*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[69*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[70*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[72*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[73*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[74*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[76*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[15*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[80*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[81*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[82*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[84*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[23*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[88*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[27*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[29*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[30*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[31*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[96*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[97*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[98*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[100*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[39*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[104*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[43*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[45*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[46*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[47*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[112*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[51*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[53*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[54*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[55*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[57*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[58*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[59*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[60*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[61*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[62*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[63*k], &a[127*k], k,	r);
      mult_pow_R(&a[96*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[97*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[98*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[99*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[100*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[101*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[102*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[104*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[105*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[106*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[108*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[112*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[113*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[114*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[116*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[120*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 32, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[32*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[33*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[34*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[35*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[36*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[37*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[38*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[40*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[41*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[42*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[44*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[15*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[48*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[49*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[50*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[52*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[23*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[56*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[27*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[29*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[30*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[31*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[96*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[97*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[98*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[100*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[71*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[104*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[75*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[77*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[78*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[79*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[112*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[83*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[85*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[86*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[87*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[89*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[90*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[91*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[92*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[93*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[94*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[95*k], &a[127*k], k,	r);
      mult_pow_R(&a[48*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[49*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[50*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[51*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[52*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[53*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[54*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[56*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[57*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[58*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[60*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[80*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[81*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[82*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[83*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[84*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[85*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[86*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[88*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[89*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[90*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[92*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[112*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[113*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[114*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[116*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[120*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 48, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[16*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[17*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[18*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[19*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[20*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[21*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[22*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[24*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[25*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[26*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[28*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[15*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[48*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[49*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[50*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[52*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[39*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[56*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[43*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[45*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[46*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[47*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[80*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[81*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[82*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[84*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[71*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[88*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[75*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[77*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[78*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[79*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[112*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[99*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[101*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[102*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[103*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[105*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[106*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[107*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[108*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[109*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[110*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[111*k], &a[127*k], k,	r);
      mult_pow_R(&a[24*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[25*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[26*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[27*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[28*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[29*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[30*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[40*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[41*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[42*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[43*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[44*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[45*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[46*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[56*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[57*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[58*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[60*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[72*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[73*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[74*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[75*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[76*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[77*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[78*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[88*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[89*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[90*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[92*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[104*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[105*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[106*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[108*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[120*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 56, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[8*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[9*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[10*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[11*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[12*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[13*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[14*k], k,	r);
      DFT2_big_elements(&a[7*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[24*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[25*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[26*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[28*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[23*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[40*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[41*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[42*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[44*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[39*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[56*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[51*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[53*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[54*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[55*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[72*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[73*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[74*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[76*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[71*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[88*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[83*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[85*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[86*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[87*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[104*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[99*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[101*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[102*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[103*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[120*k], k,	r);
      DFT2_big_elements(&a[113*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[114*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[115*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[116*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[117*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[118*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[119*k], &a[127*k], k,	r);
      mult_pow_R(&a[12*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[13*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[14*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[15*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[20*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[21*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[22*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[23*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[28*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[29*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[30*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[36*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[37*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[38*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[39*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[44*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[45*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[46*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[52*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[53*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[54*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[60*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[68*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[69*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[70*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[71*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[76*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[77*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[78*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[84*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[85*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[86*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[92*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[100*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[101*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[102*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[108*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[116*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[124*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 60, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[4*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[5*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[6*k], k,	r);
      DFT2_big_elements(&a[3*k], &a[7*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[12*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[13*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[14*k], k,	r);
      DFT2_big_elements(&a[11*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[20*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[21*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[22*k], k,	r);
      DFT2_big_elements(&a[19*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[28*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[27*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[36*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[37*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[38*k], k,	r);
      DFT2_big_elements(&a[35*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[44*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[43*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[52*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[51*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[60*k], k,	r);
      DFT2_big_elements(&a[57*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[58*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[59*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[68*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[69*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[70*k], k,	r);
      DFT2_big_elements(&a[67*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[76*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[75*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[84*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[83*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[92*k], k,	r);
      DFT2_big_elements(&a[89*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[90*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[91*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[100*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[99*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[108*k], k,	r);
      DFT2_big_elements(&a[105*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[106*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[107*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[116*k], k,	r);
      DFT2_big_elements(&a[113*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[114*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[115*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[120*k], &a[124*k], k,	r);
      DFT2_big_elements(&a[121*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[122*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[123*k], &a[127*k], k,	r);
      mult_pow_R(&a[6*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[7*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[10*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[11*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[14*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[15*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[18*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[19*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[22*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[23*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[26*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[27*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[30*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[34*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[35*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[38*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[39*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[42*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[43*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[46*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[50*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[51*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[54*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[58*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[62*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[66*k], 2, k, r, compute_inverse);
      mult_pow_R(&a[67*k], 2, k, r, compute_inverse);
      mult_pow_R(&a[70*k], 34, k, r, compute_inverse);
      mult_pow_R(&a[71*k], 34, k, r, compute_inverse);
      mult_pow_R(&a[74*k], 18, k, r, compute_inverse);
      mult_pow_R(&a[75*k], 18, k, r, compute_inverse);
      mult_pow_R(&a[78*k], 50, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 50, k, r, compute_inverse);
      mult_pow_R(&a[82*k], 10, k, r, compute_inverse);
      mult_pow_R(&a[83*k], 10, k, r, compute_inverse);
      mult_pow_R(&a[86*k], 42, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 42, k, r, compute_inverse);
      mult_pow_R(&a[90*k], 26, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 26, k, r, compute_inverse);
      mult_pow_R(&a[94*k], 58, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 58, k, r, compute_inverse);
      mult_pow_R(&a[98*k], 6, k, r, compute_inverse);
      mult_pow_R(&a[99*k], 6, k, r, compute_inverse);
      mult_pow_R(&a[102*k], 38, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 38, k, r, compute_inverse);
      mult_pow_R(&a[106*k], 22, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 22, k, r, compute_inverse);
      mult_pow_R(&a[110*k], 54, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 54, k, r, compute_inverse);
      mult_pow_R(&a[114*k], 14, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 14, k, r, compute_inverse);
      mult_pow_R(&a[118*k], 46, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 46, k, r, compute_inverse);
      mult_pow_R(&a[122*k], 30, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 30, k, r, compute_inverse);
      mult_pow_R(&a[126*k], 62, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 62, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[2*k], k,	r);
      DFT2_big_elements(&a[1*k], &a[3*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[6*k], k,	r);
      DFT2_big_elements(&a[5*k], &a[7*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[10*k], k,	r);
      DFT2_big_elements(&a[9*k], &a[11*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[14*k], k,	r);
      DFT2_big_elements(&a[13*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[18*k], k,	r);
      DFT2_big_elements(&a[17*k], &a[19*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[22*k], k,	r);
      DFT2_big_elements(&a[21*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[26*k], k,	r);
      DFT2_big_elements(&a[25*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[30*k], k,	r);
      DFT2_big_elements(&a[29*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[34*k], k,	r);
      DFT2_big_elements(&a[33*k], &a[35*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[38*k], k,	r);
      DFT2_big_elements(&a[37*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[42*k], k,	r);
      DFT2_big_elements(&a[41*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[46*k], k,	r);
      DFT2_big_elements(&a[45*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[50*k], k,	r);
      DFT2_big_elements(&a[49*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[54*k], k,	r);
      DFT2_big_elements(&a[53*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[58*k], k,	r);
      DFT2_big_elements(&a[57*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[60*k], &a[62*k], k,	r);
      DFT2_big_elements(&a[61*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[66*k], k,	r);
      DFT2_big_elements(&a[65*k], &a[67*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[70*k], k,	r);
      DFT2_big_elements(&a[69*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[74*k], k,	r);
      DFT2_big_elements(&a[73*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[78*k], k,	r);
      DFT2_big_elements(&a[77*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[82*k], k,	r);
      DFT2_big_elements(&a[81*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[86*k], k,	r);
      DFT2_big_elements(&a[85*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[90*k], k,	r);
      DFT2_big_elements(&a[89*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[92*k], &a[94*k], k,	r);
      DFT2_big_elements(&a[93*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[98*k], k,	r);
      DFT2_big_elements(&a[97*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[102*k], k,	r);
      DFT2_big_elements(&a[101*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[106*k], k,	r);
      DFT2_big_elements(&a[105*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[108*k], &a[110*k], k,	r);
      DFT2_big_elements(&a[109*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[114*k], k,	r);
      DFT2_big_elements(&a[113*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[116*k], &a[118*k], k,	r);
      DFT2_big_elements(&a[117*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[120*k], &a[122*k], k,	r);
      DFT2_big_elements(&a[121*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[124*k], &a[126*k], k,	r);
      DFT2_big_elements(&a[125*k], &a[127*k], k,	r);
      mult_pow_R(&a[3*k], 32, k, r, compute_inverse);
      mult_pow_R(&a[5*k], 16, k, r, compute_inverse);
      mult_pow_R(&a[7*k], 48, k, r, compute_inverse);
      mult_pow_R(&a[9*k], 8, k, r, compute_inverse);
      mult_pow_R(&a[11*k], 40, k, r, compute_inverse);
      mult_pow_R(&a[13*k], 24, k, r, compute_inverse);
      mult_pow_R(&a[15*k], 56, k, r, compute_inverse);
      mult_pow_R(&a[17*k], 4, k, r, compute_inverse);
      mult_pow_R(&a[19*k], 36, k, r, compute_inverse);
      mult_pow_R(&a[21*k], 20, k, r, compute_inverse);
      mult_pow_R(&a[23*k], 52, k, r, compute_inverse);
      mult_pow_R(&a[25*k], 12, k, r, compute_inverse);
      mult_pow_R(&a[27*k], 44, k, r, compute_inverse);
      mult_pow_R(&a[29*k], 28, k, r, compute_inverse);
      mult_pow_R(&a[31*k], 60, k, r, compute_inverse);
      mult_pow_R(&a[33*k], 2, k, r, compute_inverse);
      mult_pow_R(&a[35*k], 34, k, r, compute_inverse);
      mult_pow_R(&a[37*k], 18, k, r, compute_inverse);
      mult_pow_R(&a[39*k], 50, k, r, compute_inverse);
      mult_pow_R(&a[41*k], 10, k, r, compute_inverse);
      mult_pow_R(&a[43*k], 42, k, r, compute_inverse);
      mult_pow_R(&a[45*k], 26, k, r, compute_inverse);
      mult_pow_R(&a[47*k], 58, k, r, compute_inverse);
      mult_pow_R(&a[49*k], 6, k, r, compute_inverse);
      mult_pow_R(&a[51*k], 38, k, r, compute_inverse);
      mult_pow_R(&a[53*k], 22, k, r, compute_inverse);
      mult_pow_R(&a[55*k], 54, k, r, compute_inverse);
      mult_pow_R(&a[57*k], 14, k, r, compute_inverse);
      mult_pow_R(&a[59*k], 46, k, r, compute_inverse);
      mult_pow_R(&a[61*k], 30, k, r, compute_inverse);
      mult_pow_R(&a[63*k], 62, k, r, compute_inverse);
      mult_pow_R(&a[65*k], 1, k, r, compute_inverse);
      mult_pow_R(&a[67*k], 33, k, r, compute_inverse);
      mult_pow_R(&a[69*k], 17, k, r, compute_inverse);
      mult_pow_R(&a[71*k], 49, k, r, compute_inverse);
      mult_pow_R(&a[73*k], 9, k, r, compute_inverse);
      mult_pow_R(&a[75*k], 41, k, r, compute_inverse);
      mult_pow_R(&a[77*k], 25, k, r, compute_inverse);
      mult_pow_R(&a[79*k], 57, k, r, compute_inverse);
      mult_pow_R(&a[81*k], 5, k, r, compute_inverse);
      mult_pow_R(&a[83*k], 37, k, r, compute_inverse);
      mult_pow_R(&a[85*k], 21, k, r, compute_inverse);
      mult_pow_R(&a[87*k], 53, k, r, compute_inverse);
      mult_pow_R(&a[89*k], 13, k, r, compute_inverse);
      mult_pow_R(&a[91*k], 45, k, r, compute_inverse);
      mult_pow_R(&a[93*k], 29, k, r, compute_inverse);
      mult_pow_R(&a[95*k], 61, k, r, compute_inverse);
      mult_pow_R(&a[97*k], 3, k, r, compute_inverse);
      mult_pow_R(&a[99*k], 35, k, r, compute_inverse);
      mult_pow_R(&a[101*k], 19, k, r, compute_inverse);
      mult_pow_R(&a[103*k], 51, k, r, compute_inverse);
      mult_pow_R(&a[105*k], 11, k, r, compute_inverse);
      mult_pow_R(&a[107*k], 43, k, r, compute_inverse);
      mult_pow_R(&a[109*k], 27, k, r, compute_inverse);
      mult_pow_R(&a[111*k], 59, k, r, compute_inverse);
      mult_pow_R(&a[113*k], 7, k, r, compute_inverse);
      mult_pow_R(&a[115*k], 39, k, r, compute_inverse);
      mult_pow_R(&a[117*k], 23, k, r, compute_inverse);
      mult_pow_R(&a[119*k], 55, k, r, compute_inverse);
      mult_pow_R(&a[121*k], 15, k, r, compute_inverse);
      mult_pow_R(&a[123*k], 47, k, r, compute_inverse);
      mult_pow_R(&a[125*k], 31, k, r, compute_inverse);
      mult_pow_R(&a[127*k], 63, k, r, compute_inverse);
      DFT2_big_elements(&a[0*k], &a[1*k], k,	r);
      DFT2_big_elements(&a[2*k], &a[3*k], k,	r);
      DFT2_big_elements(&a[4*k], &a[5*k], k,	r);
      DFT2_big_elements(&a[6*k], &a[7*k], k,	r);
      DFT2_big_elements(&a[8*k], &a[9*k], k,	r);
      DFT2_big_elements(&a[10*k], &a[11*k], k,	r);
      DFT2_big_elements(&a[12*k], &a[13*k], k,	r);
      DFT2_big_elements(&a[14*k], &a[15*k], k,	r);
      DFT2_big_elements(&a[16*k], &a[17*k], k,	r);
      DFT2_big_elements(&a[18*k], &a[19*k], k,	r);
      DFT2_big_elements(&a[20*k], &a[21*k], k,	r);
      DFT2_big_elements(&a[22*k], &a[23*k], k,	r);
      DFT2_big_elements(&a[24*k], &a[25*k], k,	r);
      DFT2_big_elements(&a[26*k], &a[27*k], k,	r);
      DFT2_big_elements(&a[28*k], &a[29*k], k,	r);
      DFT2_big_elements(&a[30*k], &a[31*k], k,	r);
      DFT2_big_elements(&a[32*k], &a[33*k], k,	r);
      DFT2_big_elements(&a[34*k], &a[35*k], k,	r);
      DFT2_big_elements(&a[36*k], &a[37*k], k,	r);
      DFT2_big_elements(&a[38*k], &a[39*k], k,	r);
      DFT2_big_elements(&a[40*k], &a[41*k], k,	r);
      DFT2_big_elements(&a[42*k], &a[43*k], k,	r);
      DFT2_big_elements(&a[44*k], &a[45*k], k,	r);
      DFT2_big_elements(&a[46*k], &a[47*k], k,	r);
      DFT2_big_elements(&a[48*k], &a[49*k], k,	r);
      DFT2_big_elements(&a[50*k], &a[51*k], k,	r);
      DFT2_big_elements(&a[52*k], &a[53*k], k,	r);
      DFT2_big_elements(&a[54*k], &a[55*k], k,	r);
      DFT2_big_elements(&a[56*k], &a[57*k], k,	r);
      DFT2_big_elements(&a[58*k], &a[59*k], k,	r);
      DFT2_big_elements(&a[60*k], &a[61*k], k,	r);
      DFT2_big_elements(&a[62*k], &a[63*k], k,	r);
      DFT2_big_elements(&a[64*k], &a[65*k], k,	r);
      DFT2_big_elements(&a[66*k], &a[67*k], k,	r);
      DFT2_big_elements(&a[68*k], &a[69*k], k,	r);
      DFT2_big_elements(&a[70*k], &a[71*k], k,	r);
      DFT2_big_elements(&a[72*k], &a[73*k], k,	r);
      DFT2_big_elements(&a[74*k], &a[75*k], k,	r);
      DFT2_big_elements(&a[76*k], &a[77*k], k,	r);
      DFT2_big_elements(&a[78*k], &a[79*k], k,	r);
      DFT2_big_elements(&a[80*k], &a[81*k], k,	r);
      DFT2_big_elements(&a[82*k], &a[83*k], k,	r);
      DFT2_big_elements(&a[84*k], &a[85*k], k,	r);
      DFT2_big_elements(&a[86*k], &a[87*k], k,	r);
      DFT2_big_elements(&a[88*k], &a[89*k], k,	r);
      DFT2_big_elements(&a[90*k], &a[91*k], k,	r);
      DFT2_big_elements(&a[92*k], &a[93*k], k,	r);
      DFT2_big_elements(&a[94*k], &a[95*k], k,	r);
      DFT2_big_elements(&a[96*k], &a[97*k], k,	r);
      DFT2_big_elements(&a[98*k], &a[99*k], k,	r);
      DFT2_big_elements(&a[100*k], &a[101*k], k,	r);
      DFT2_big_elements(&a[102*k], &a[103*k], k,	r);
      DFT2_big_elements(&a[104*k], &a[105*k], k,	r);
      DFT2_big_elements(&a[106*k], &a[107*k], k,	r);
      DFT2_big_elements(&a[108*k], &a[109*k], k,	r);
      DFT2_big_elements(&a[110*k], &a[111*k], k,	r);
      DFT2_big_elements(&a[112*k], &a[113*k], k,	r);
      DFT2_big_elements(&a[114*k], &a[115*k], k,	r);
      DFT2_big_elements(&a[116*k], &a[117*k], k,	r);
      DFT2_big_elements(&a[118*k], &a[119*k], k,	r);
      DFT2_big_elements(&a[120*k], &a[121*k], k,	r);
      DFT2_big_elements(&a[122*k], &a[123*k], k,	r);
      DFT2_big_elements(&a[124*k], &a[125*k], k,	r);
      DFT2_big_elements(&a[126*k], &a[127*k], k,	r);
      swap_big_elements (&a[1*k], &a[64*k], k);
      swap_big_elements (&a[2*k], &a[32*k], k);
      swap_big_elements (&a[3*k], &a[96*k], k);
      swap_big_elements (&a[4*k], &a[16*k], k);
      swap_big_elements (&a[5*k], &a[80*k], k);
      swap_big_elements (&a[6*k], &a[48*k], k);
      swap_big_elements (&a[7*k], &a[112*k], k);
      swap_big_elements (&a[9*k], &a[72*k], k);
      swap_big_elements (&a[10*k], &a[40*k], k);
      swap_big_elements (&a[11*k], &a[104*k], k);
      swap_big_elements (&a[12*k], &a[24*k], k);
      swap_big_elements (&a[13*k], &a[88*k], k);
      swap_big_elements (&a[14*k], &a[56*k], k);
      swap_big_elements (&a[15*k], &a[120*k], k);
      swap_big_elements (&a[17*k], &a[68*k], k);
      swap_big_elements (&a[18*k], &a[36*k], k);
      swap_big_elements (&a[19*k], &a[100*k], k);
      swap_big_elements (&a[21*k], &a[84*k], k);
      swap_big_elements (&a[22*k], &a[52*k], k);
      swap_big_elements (&a[23*k], &a[116*k], k);
      swap_big_elements (&a[25*k], &a[76*k], k);
      swap_big_elements (&a[26*k], &a[44*k], k);
      swap_big_elements (&a[27*k], &a[108*k], k);
      swap_big_elements (&a[29*k], &a[92*k], k);
      swap_big_elements (&a[30*k], &a[60*k], k);
      swap_big_elements (&a[31*k], &a[124*k], k);
      swap_big_elements (&a[33*k], &a[66*k], k);
      swap_big_elements (&a[35*k], &a[98*k], k);
      swap_big_elements (&a[37*k], &a[82*k], k);
      swap_big_elements (&a[38*k], &a[50*k], k);
      swap_big_elements (&a[39*k], &a[114*k], k);
      swap_big_elements (&a[41*k], &a[74*k], k);
      swap_big_elements (&a[43*k], &a[106*k], k);
      swap_big_elements (&a[45*k], &a[90*k], k);
      swap_big_elements (&a[46*k], &a[58*k], k);
      swap_big_elements (&a[47*k], &a[122*k], k);
      swap_big_elements (&a[49*k], &a[70*k], k);
      swap_big_elements (&a[51*k], &a[102*k], k);
      swap_big_elements (&a[53*k], &a[86*k], k);
      swap_big_elements (&a[55*k], &a[118*k], k);
      swap_big_elements (&a[57*k], &a[78*k], k);
      swap_big_elements (&a[59*k], &a[110*k], k);
      swap_big_elements (&a[61*k], &a[94*k], k);
      swap_big_elements (&a[63*k], &a[126*k], k);
      swap_big_elements (&a[67*k], &a[97*k], k);
      swap_big_elements (&a[69*k], &a[81*k], k);
      swap_big_elements (&a[71*k], &a[113*k], k);
      swap_big_elements (&a[75*k], &a[105*k], k);
      swap_big_elements (&a[77*k], &a[89*k], k);
      swap_big_elements (&a[79*k], &a[121*k], k);
      swap_big_elements (&a[83*k], &a[101*k], k);
      swap_big_elements (&a[87*k], &a[117*k], k);
      swap_big_elements (&a[91*k], &a[109*k], k);
      swap_big_elements (&a[95*k], &a[125*k], k);
      swap_big_elements (&a[103*k], &a[115*k], k);
      swap_big_elements (&a[111*k], &a[123*k], k);

    }
}

/**************************************/

/** compute L_{m}{n} for large coefficients;
 * equivalently, computing transpose {n}{m}
 */
void
stride_permutation_big_elements (usfixn64* A, usfixn64 *B, int m, int n,
				 int n_permutations, int coefficient_size)
{
#if PROFILING_ENABLED ==1

  if (global_profiling_ongoing == 0)
    {
      stride_permutation_called_dims[n_stride_permutation_called][0] = m;
      stride_permutation_called_dims[n_stride_permutation_called][1] = n;
    }
  inc_profiling_counter (&n_stride_permutation_called);
#endif
//  int blocksize = min_three (m, n, PERMUTATION_BLOCKSIZE);
//  usfixn64* B = (usfixn64*) malloc (
//      m * n * coefficient_size * sizeof(usfixn64));

  int n_bytes_per_coefficient = coefficient_size * sizeof(usfixn64);
  //stride = K*J*coefficient_size
  usfixn64 stride = m * n * (coefficient_size);

  int m_log2 = log2 (m);
  int mn=m*n;

#pragma cilk_grainsize = ((n_permutations)/(cilk::current_worker_count))
  cilk_for (int t = 0; t < n_permutations; t++)
#pragma cilk_grainsize = ((mn)/(cilk::current_worker_count))
    cilk_for (int ij = 0; ij < mn; ij++)
      {
//	   int i=ij/m;
//	   int j=ij%m;
	int i = ij >> m_log2;
	int j = ij & (m - 1);
//	idx_src = j + i * m;
	int idx_dest = i + j * n;
	int idx_src = ij;
	memcpy (&B[idx_dest * coefficient_size + t * stride],
		&A[idx_src * coefficient_size + t * stride],
		n_bytes_per_coefficient);
      }
}

/**************************************/

void
mult_pow_R_gmp (mpz_t x_zz, const int s, const usfixn64 k, const usfixn64 r,
		const mpz_t p_zz)
{

  if (s == 0)
    return;
  mpz_t radix;
  mpz_init_set_ui (radix, r);
  mpz_pow_ui (radix, radix, s);
  mpz_mod (radix, radix, p_zz);

  mpz_mul (x_zz, x_zz, radix);
  mpz_mod (x_zz, x_zz, p_zz);

  mpz_clear (radix);
}

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
		      float *t_mult_values, int* n_twiddle_mult)
{
#if PROFILING_ENABLED ==1

  if (global_profiling_ongoing == 0)
    {
      twiddle_called_dims[n_twiddle_called][0] = K;
      twiddle_called_dims[n_twiddle_called][1] = J;
    }
  inc_profiling_counter (&n_twiddle_called);

#endif

  int J_log2 = log2 (J);
  int J_mask = J - 1;

  //K*J*coefficient_size
  usfixn64 stride = K * J * (K / 2);
  usfixn64 N = K * J * n_permutations;

  /////////////////////////////////////
  cpu_timer t_mult_pow_r, t_mult;

  timer_record_start (&t_mult_pow_r);
  int KJ = K * J;
#pragma cilk_grainsize = ((n_permutations*KJ)/(cilk::current_worker_count))
  cilk_for (int t = 0; t < n_permutations; t++)
    cilk_for (int ij = 0; ij < KJ; ij++)
      {
	int idx, q, m;
	int i = ij >> J_log2;
	int j = ij & (J_mask);
	idx = i * j;
	//	    q = idx / J;
	//	    m = idx % J;

	q = idx >> J_log2;
	m = idx & (J_mask);
	  {
	    idx = (i * J + j) * k + t * stride;
	    //cyclic shift part
	    mult_pow_R (&vector[idx], q, k, r, compute_inverse);
	  }
      }

  timer_record_stop (&t_mult_pow_r);

  n_twiddle_mult[0] += (n_permutations * K * J);

  /////////////////////////////////////

//  int * n_gfpf_mult=(int*)malloc(n_permutations)
  timer_record_start (&t_mult);
  if (fft_based_mult_enabled == 1)
    {
#pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
	cilk_for (int t = 0; t < n_permutations; t++)
	  cilk_for (int ij = 0; ij < KJ; ij++)
	    {
	      int idx, q, m;
	      int i = ij >> J_log2;
	      int j = ij & (J_mask);
	      idx = i * j;
	      //	    q = idx / J;
	      //	    m = idx % J;
	      q = idx >> J_log2;
	      m = idx & (J_mask);

	      idx = (i * J + j) * k + t * stride;
	      //		if (m != 0)
		{
		  GFPFMultiplication (&vector[idx],
				      &precomputed_pow_omega_vec[(m) * k], k,
				      &t_crt_data_global, &t_lhc_data_global);
		  //		    n_twiddle_mult[1]++;
		}
	    }
    }
  else
    {
#pragma cilk_grainsize = (n_permutations/(cilk::current_worker_count))
//      for (int i = 0; i < K; i++)
//	for (int j = 0; j < J; j++)
//	  {
	cilk_for (int t = 0; t < n_permutations; t++)
	  cilk_for (int ij = 0; ij < KJ; ij++)
	    {
	      int idx, q, m;
	      int i = ij >> J_log2;
	      int j = ij & (J_mask);
	      idx = i * j;
	      //	    q = idx / J;
	      //	    m = idx % J;
	      q = idx >> J_log2;
	      m = idx & (J_mask);
		{
		  idx = (i * J + j) * k + t * stride;
//		if (m != 0)
		    {

		      mult_gmp_big_elements (
			  &vector[idx], &precomputed_pow_omega_vec[(m) * k], k,
			  r, p);
//		    n_twiddle_mult[1]++;
		    }
		}
	    }

    }

  timer_record_stop (&t_mult);
  n_twiddle_mult[1] += (n_permutations * K * J);

  /////////////////////////////////////

//  free (memory);
  timer_get_elapsed_time (&t_mult_pow_r, NULL, 1);
  timer_get_elapsed_time (&t_mult, NULL, 1);

  t_mult_values[0] = t_mult_pow_r.elapsed_time;
  t_mult_values[1] = t_mult.elapsed_time;
}

/**************************************/

/** precompute powers of omega up to (omega)^{K^b}
* to be used at level b+1
* array "usfixn64* precomputed_pow_omega_vec" should be be
* already allocated and passed as a pointer.
*/
void
precompute_sequential_powers_of_omega_big_elements (
    usfixn64* precomputed_pow_omega_vec, int K, int b, usfixn64 radix,
    mpz_t omega, mpz_t p)
{
  int coefficient_size = K / 2;
  int n = pow (K, b);

  int n_bytes = coefficient_size * sizeof(usfixn64);
  usfixn64 * tmp = (usfixn64*) malloc (n_bytes);

  mpz_t current_twiddle;
  mpz_init_set_ui (current_twiddle, 1);

//	printf("-- precompute powers of omega-at-level-%d up to omega^{K^%d = %d} \n   ",
//			b + 1, b, n);
  for (int i = 0; i < n; i++)
    {
//		progress_bar(i, n);
      if (i != 0)
	{
	  mpz_mul (current_twiddle, current_twiddle, omega);
	  mpz_mod (current_twiddle, current_twiddle, p);
	}
//		printf("pow_omega = %d \n", i);
      bigint_to_u64vector_radixbased_gmp (tmp, current_twiddle, radix,
					  coefficient_size);
      memcpy (&precomputed_pow_omega_vec[i * coefficient_size], tmp, n_bytes);
    }

  free (tmp);
}

/**************************************/
//m=(ab)/c
void mul_ab_div_c_mpq(float *m, float a, int b, int c)
{

  double result;
  mpq_t a_mpq, b_mpq, c_mpq, m_mpq;

  mpq_inits(a_mpq, b_mpq, c_mpq, m_mpq, NULL);
  mpq_set_d(a_mpq, (double)a);
  mpq_set_ui(b_mpq, (usfixn32)b, 1);
  mpq_set_ui(c_mpq, (usfixn32)c, 1);

  mpq_mul(m_mpq, a_mpq, b_mpq);
  mpq_div(m_mpq, m_mpq, c_mpq);

  *m=(float)mpq_get_d(m_mpq);

  mpq_clears(a_mpq, b_mpq, c_mpq, m_mpq, NULL);
}

/**************************************/

/** compute DFT_{N=K^e} for p=radix^{K/2}+1;
* omega is omega at the highest level such that (omega^N)=1 mod p
*/
void
DFT_general_big_elements (usfixn64* vector_in, int K, int e, usfixn64 radix,
			  mpz_t p, usfixn64* omega, int fft_based_mult_enabled,
			  int compute_inverse, int verbose)
{
  //we assume that the highest level is e and the lowest level is 2.
  usfixn64 N, J;
  N = pow (K, e);
  int coefficient_size = K / 2;
  int k = K / 2;

  //precompute sequential powers of omega
  usfixn64** precomputed_pow_omega_vec_2D = (usfixn64**) malloc (
      (e + 1) * sizeof(usfixn64*));

  usfixn64 *vector = vector_in;
  usfixn64 *aux_vector = (usfixn64*) malloc (N * k * sizeof(usfixn64));

  //convert omega at the highest to mpz representation
  mpz_t omega_zz;
  mpz_init (omega_zz);
  u64vector_to_radix_based_bigint_gmp (omega_zz, omega, radix,
				       coefficient_size);

  cpu_timer t_timer;
  float t_stride_total = 0, t_DFT_K_total = 0, t_twiddle_total = 0,
      t_precompute, t_mult_pow_r_total = 0, t_arbitrary_mult_total = 0;

  //n_twiddle_mult: [0]:n_mult_pow_r [1]: n_arbitrary_mult
  int n_twiddle_mult[2] =
    { 0, 0 }, n_DFT_K = 0, n_mult_pow_r_in_DFT_K = 0,
      n_mult_pow_r_per_DFT_K = 0;

  /////////////////////////////////////
  /////////////////////////////////////
  timer_record_start (&t_timer);
  for (int b = e; b >= 2; b--)
    {
#if VERBOSE>=2
      printf ("[precomputing pow_omega for level b=%d]\n", b);
#endif
      J = pow (K, b - 1);
      precomputed_pow_omega_vec_2D[b] = (usfixn64*) malloc (
	  J * coefficient_size * sizeof(usfixn64));

      //precompute sequential powers of omega at level=b and store it in 2D vector.
      precompute_sequential_powers_of_omega_big_elements (
	  precomputed_pow_omega_vec_2D[b], K, b - 1, radix, omega_zz, p);

      //as we go to a lower level, w=(w^K) mod p
      mpz_powm_ui (omega_zz, omega_zz, K, p);
//	printVectorToFile(precomputed_pow_omega_vec_2D[b], J, coefficient_size, "pow_omega");
    }
  timer_record_stop (&t_timer);
  timer_get_elapsed_time (&t_timer, NULL, 1);
  t_precompute = t_timer.elapsed_time;

  /////////////////////////////////////
  /////////////////////////////////////

  int aux_memory_check_bit = 0;
  J = pow (K, e - 1);
//  usfixn64 stride;

  // step 1: stride permutations on the right-most
  // side of the DFT computing L{K,J}
  cpu_timer t_dft_general;
  timer_record_start (&t_dft_general);
  timer_record_start (&t_timer);
  for (int b = e; b >= 2; b--)
    {
      J = pow (K, b - 1);
      int n_permutations = N / (K * J);
      //printf("-- step 1: (%-12s) b=%d , K=%d , J=%d\n", "right perm", b, K,J);
//      stride = (K * J * coefficient_size);
      stride_permutation_big_elements (vector, aux_vector, K, J, n_permutations,
				       coefficient_size);
      swap_ptr (&vector, &aux_vector, &aux_memory_check_bit);
    }
  timer_record_stop (&t_timer);
  timer_get_elapsed_time (&t_timer, NULL, 1);
  t_stride_total += t_timer.elapsed_time;
  //step1 is verified, can be done more rigorously.

  /////////////////////////////////////
  /////////////////////////////////////

// printf("-- step 2: (%-12s) b=%d , K=%d , J=%d \n", "right DFT_K", 2, K, J);
// step 2: base-case DFT's at the lowest level

  J = pow (K, e - 1);
//  stride = (K * coefficient_size);

  void
  (*DFT_K_ptr) (int, usfixn64*, const int, const usfixn64, int);

  timer_record_start (&t_timer);
  switch (K)
    {
    case 4:
      DFT_K_ptr = &DFT_4_big_elements;
      n_mult_pow_r_per_DFT_K = 1;
      break;
    case 8:
      DFT_K_ptr = &DFT_8_big_elements;
      n_mult_pow_r_per_DFT_K = 5;
      break;
    case 16:
      DFT_K_ptr = &DFT_16_big_elements;
      n_mult_pow_r_per_DFT_K = 17;
      break;
    case 32:
      DFT_K_ptr = &DFT_32_big_elements;
      n_mult_pow_r_per_DFT_K = 66;
      break;
    case 64:
      DFT_K_ptr = &DFT_64_big_elements;
      n_mult_pow_r_per_DFT_K = 195;
      break;
    case 128:
      //ToDo: fix the following values.
      DFT_K_ptr = &DFT_128_big_elements;
      n_mult_pow_r_per_DFT_K = 321;
      break;
    case 256:
      //ToDo: fix the following values.
      DFT_K_ptr = &DFT_256_big_elements;
      n_mult_pow_r_per_DFT_K = 769;
      break;
    }



//  CALLGRIND_START_INSTRUMENTATION;
  DFT_K_ptr (J, vector, coefficient_size, radix, compute_inverse);
  n_DFT_K += J;
  timer_record_stop (&t_timer);
  timer_get_elapsed_time (&t_timer, NULL, 1);
  t_DFT_K_total += t_timer.elapsed_time;
//  CALLGRIND_STOP_INSTRUMENTATION;
//  CALLGRIND_DUMP_STATS;
//  printf("EXIT AFTER DFT_K GFPF\n");
//  exit(0);


//ToDO: automated verification of base-case DFT's: DONE.

  ////////////////////////////////////////////////////////////
  //// following printToFile is used for testing purposes
  //// printVectorToFile(vector, N, coefficient_size, "xData_after_base_cases");
  ////////////////////////////////////////////////////////////

  for (int b = 2; b <= e; b++)
    {
      //step 3: twiddle DFT_GRAIN_FACTOR multiplication + permutation in between + DFT_K on the left
      //beginning at the level b=2;
      J = pow (K, b - 1);
      int n_permutations = N / (K * J);

//      printf("-- step 3: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
//             "twiddle+perm", b, K, J, n_permutations);
//      timer_record_start (&t_timer);
      float t_twiddle[2] =
	{ 0, 0 };

//      CALLGRIND_START_INSTRUMENTATION;
      twiddle_big_elements (vector, K, J, n_permutations,
			    precomputed_pow_omega_vec_2D[b], coefficient_size,
			    radix, p, fft_based_mult_enabled, compute_inverse,
			    t_twiddle, n_twiddle_mult);
//      CALLGRIND_STOP_INSTRUMENTATION;
//      CALLGRIND_DUMP_STATS;
//      printf("EXIT AFTER TWIDDLE GFPF\n");
//      exit(0);
//      timer_record_stop (&t_timer);
//      timer_get_elapsed_time (&t_timer, NULL, 1);
//      t_twiddle_total += t_timer.elapsed_time;

      t_mult_pow_r_total += t_twiddle[0];
      t_arbitrary_mult_total += t_twiddle[1];

      /////////////////////////////////////
      /////////////////////////////////////
      timer_record_start (&t_timer);
      stride_permutation_big_elements (vector, aux_vector, J, K, n_permutations,
				       coefficient_size);
      swap_ptr (&vector, &aux_vector, &aux_memory_check_bit);
      timer_record_stop (&t_timer);
      timer_get_elapsed_time (&t_timer, NULL, 1);
      t_stride_total += t_timer.elapsed_time;

      /////////////////////////////////////
      /////////////////////////////////////

      //step 4:
      //computing DFT_K on the left.
      J = N / K;
      n_permutations = N / (K * J);
      timer_record_start (&t_timer);
      DFT_K_ptr (J, vector, coefficient_size, radix, compute_inverse);
      timer_record_stop (&t_timer);
      timer_get_elapsed_time (&t_timer, NULL, 1);
      t_DFT_K_total += t_timer.elapsed_time;
      n_DFT_K += J;
      /////////////////////////////////////
      /////////////////////////////////////

      // step 5: left-most permutation L_{K,K^e-1} on the right
      J = pow (K, b - 1);
      n_permutations = N / (K * J);
      //printf("-- step 5: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n", "left perm", b, K, J, n_permutations);

      timer_record_start (&t_timer);
      stride_permutation_big_elements (vector, aux_vector, K, J, n_permutations,
				       coefficient_size);
      swap_ptr (&vector, &aux_vector, &aux_memory_check_bit);
      timer_record_stop (&t_timer);
      timer_get_elapsed_time (&t_timer, NULL, 1);
      t_stride_total += t_timer.elapsed_time;

      /////////////////////////////////////
      /////////////////////////////////////
    }

  timer_record_stop (&t_dft_general);
  timer_get_elapsed_time (&t_dft_general, NULL, 1);

  if (aux_memory_check_bit == 1)
    {
      swap_ptr (&vector, &aux_vector, &aux_memory_check_bit);
      memcpy (vector, aux_vector, N * k * sizeof(usfixn64));
    }

  int n_mult_pow_r_in_DFT_K_total = (n_mult_pow_r_per_DFT_K * n_DFT_K);
  n_twiddle_mult[0]+=n_mult_pow_r_in_DFT_K_total;

  float t_mult_pow_r_avg = t_mult_pow_r_total / n_twiddle_mult[0];

//  float t_mult_pow_r_in_DFT_K_total=(t_mult_pow_r_avg*n_mult_pow_r_in_DFT_K_total);
//  float t_mult_pow_r_in_DFT_K_total = (t_mult_pow_r_total
//      * n_mult_pow_r_in_DFT_K_total) / n_twiddle_mult[0];

  printf("n_mult_pow_r_in_DFT_K_total=%d\n", n_mult_pow_r_in_DFT_K_total);
//  printf("n_mult_pow_r_in_DFT_K_total=%d\n", n_mult_pow_r_in_DFT_K_total);
  float t_mult_pow_r_in_DFT_K_total;
//  mul_ab_div_c_mpq(&t_mult_pow_r_in_DFT_K_total, t_mult_pow_r_total, n_mult_pow_r_in_DFT_K_total, n_twiddle_mult[0]);

  t_mult_pow_r_in_DFT_K_total=t_mult_pow_r_avg*(float)n_mult_pow_r_in_DFT_K_total;

//  n_twiddle_mult[0]+=n_mult_pow_r_in_DFT_K_total;
//  t_mult_pow_r_total+=t_mult_pow_r_in_DFT_K_total;
  t_twiddle_total = t_mult_pow_r_total + t_arbitrary_mult_total;

  if (verbose)
    {

      char suffix[64];
      char msg[256];

      if (fft_based_mult_enabled)
	{
	  sprintf (suffix, "%s", "gfpf");
	}
      else
	{
	  sprintf (suffix, "%s", "gmp_based");
	}

      sprintf (msg, "%s_%s", "t_precompute", suffix);
      timer_print_time_percentage (t_precompute, msg,
				   t_dft_general.elapsed_time);

      sprintf (msg, "%s_%s", "t_stride_total", suffix);
      timer_print_time_percentage (t_stride_total, msg,
				   t_dft_general.elapsed_time);

      sprintf (msg, "%s_%s", "t_DFT_K_total", suffix);
      timer_print_time_percentage (t_DFT_K_total, msg,
				   t_dft_general.elapsed_time);

      sprintf(msg, "%s_%s", "t_DFT_K_total_avg", suffix);
      timer_print_time(t_DFT_K_total/n_DFT_K, msg);

      printf("%s", quadline);

      sprintf(msg, "%s_%s", "t_mult_pow_r_in_DFT_K_total", suffix);
      timer_print_time_percentage(t_mult_pow_r_in_DFT_K_total, msg, t_dft_general.elapsed_time);

      printf("%s", quadline);
      sprintf(msg, "%s_%s", "t_mult_pow_r_total", suffix);
      timer_print_time_percentage (t_mult_pow_r_total, msg, t_dft_general.elapsed_time);

      sprintf(msg, "%s_%s", "n_mult_pow_r_total", suffix);
      print_quantity(n_twiddle_mult[0], msg);

      sprintf(msg, "%s_%s", "t_mult_pow_r_total_avg", suffix);
      timer_print_time(t_mult_pow_r_total/n_twiddle_mult[0], msg);

      sprintf(msg, "%s_%s", "t_arbitrary_mult_total", suffix);
      timer_print_time_percentage (t_arbitrary_mult_total, msg, t_dft_general.elapsed_time);

      sprintf(msg, "%s_%s", "n_arbitrary_mult_total", suffix);
      print_quantity(n_twiddle_mult[1], msg);

      sprintf(msg, "%s_%s", "t_arbitrary_mult_total_avg", suffix);
      timer_print_time(t_arbitrary_mult_total/n_twiddle_mult[1], msg);

      sprintf (msg, "%s_%s", "t_twiddle_total", suffix);
      timer_print_time_percentage (t_twiddle_total, msg,
			 t_dft_general.elapsed_time);

      sprintf (msg, "%s_%s", "t_dft_six_step", suffix);
      timer_print_time_percentage (t_dft_general.elapsed_time, msg,
				   t_dft_general.elapsed_time);
      printf ("%s", longline);
    }

  free (aux_vector);
  for (int b = e; b >= 2; b--)
    free (precomputed_pow_omega_vec_2D[b]);

  free (precomputed_pow_omega_vec_2D);
//  free(mult_aux_vector);
}     // end DFT_general

/**************************************/

int
convolution_mult (usfixn64* x_data, usfixn64* y_data, int n, srgfn_prime *P,
		  mpz_t p_zz, int fft_based_mult_enabled)
{

  cpu_timer t_convolution_mult;
  timer_record_start (&t_convolution_mult);

  int k = P->k;
  if (fft_based_mult_enabled == 1)
    {
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
      cilk_for (int i = 0; i < n; i++)
	{
	  GFPFMultiplication (&x_data[i * k], &y_data[i * k], k,
			      &t_crt_data_global, &t_lhc_data_global);
	}
//      free (memory);
    }
  /////////////////////////////////////
  else
    {
      mpz_t *xn, *yn;
      xn = (mpz_t*) malloc (n * sizeof(mpz_t));
      yn = (mpz_t*) malloc (n * sizeof(mpz_t));

#pragma cilk_grainsize = (CONVOLUTION_GRAIN_FACTOR*(cilk::current_worker_count))
      cilk_for (int i = 0; i < n; i++)
	{
	  mpz_init (xn[i]);
	  u64vector_to_radix_based_bigint_gmp (xn[i], &x_data[i * P->k],
					       P->radix, P->k);
	}

#pragma cilk_grainsize = (CONVOLUTION_GRAIN_FACTOR*(cilk::current_worker_count))
      cilk_for (int i = 0; i < n; i++)
	{
	  mpz_init (yn[i]);
	  u64vector_to_radix_based_bigint_gmp (yn[i], &y_data[i * P->k],
					       P->radix, P->k);
	}

#pragma cilk_grainsize = (CONVOLUTION_GRAIN_FACTOR*(cilk::current_worker_count))
      cilk_for (int i = 0; i < n; i++)
	{
	  mpz_mul (xn[i], xn[i], yn[i]);
	  mpz_mod (xn[i], xn[i], p_zz);
	}

#pragma cilk_grainsize = (CONVOLUTION_GRAIN_FACTOR*(cilk::current_worker_count))
      cilk_for (int i = 0; i < n; i++)
	bigint_to_u64vector_radixbased_gmp (&x_data[i * P->k], xn[i], P->radix,
					    P->k);

#pragma cilk_grainsize = (CONVOLUTION_GRAIN_FACTOR*(cilk::current_worker_count))
      cilk_for (int i = 0; i < n; i++)
	{
	  mpz_clear (xn[i]);
	  mpz_clear (yn[i]);
	}
      free (xn);
      free (yn);

    }
  timer_record_stop (&t_convolution_mult);
  timer_get_elapsed_time (&t_convolution_mult, "t_convolution_mult", 1);
  printf ("%s", shortline);

  return EXIT_SUCCESS;
}

/**************************************/

void
test_DFT2_big_elements (usfixn64* x, int K, int e, usfixn64 radix)
{

  int n = (int) (pow (K, e));
  int k = K / 2;
  int n_bytes = n * k * sizeof(usfixn64);
  usfixn64 * y = (usfixn64*) malloc (n_bytes);
  memcpy (y, x, n_bytes);
  printf ("[testing DFT2_big_elements K=%d (k=%d), e=%d, radix=%llu] ...\n", K,
	  k, e, radix);

  cpu_timer t_dft2;

  for (int i = 0; i < n; i += 2)
    DFT2_big_elements (&x[i * k], &x[(i + 1) * k], k, radix);
  memcpy (x, y, n_bytes);
  /////////////////////////////////////
  timer_record_start (&t_dft2);
  for (int i = 0; i < n; i += 2)
    DFT2_big_elements (&x[i * k], &x[(i + 1) * k], k, radix);

  timer_record_stop (&t_dft2);
  timer_get_elapsed_time (&t_dft2, "t_dft2", 1);
  /////////////////////////////////////
  timer_record_start (&t_dft2);
  for (int i = 0; i < n; i += 2)
    DFT2_big_elements_v0 (&y[i * k], &y[(i + 1) * k], k, radix);

  timer_record_stop (&t_dft2);
  timer_get_elapsed_time (&t_dft2, "t_dft2_old", 1);
  /////////////////////////////////////
  int verification_result = EXIT_SUCCESS;
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
	if (x[i * k + j] != y[i * k + j])
	  {
	    printf ("-ERROR: mismatch at [%d, %d]\n", i, j);
	    verification_result = EXIT_FAILURE;
	    break;
	  }
      if (verification_result == EXIT_FAILURE)
	break;
    }
  print_verification_msg (__func__, verification_result);
//  printf("%s", shortline);
/////////////////////////////////////
  free (y);
}

/**************************************/

void
test_mult_pow_R_big_elements (usfixn64* x, int K, int e, usfixn64 radix)
{

  int n = (int) (pow (K, e));
  int k = K / 2;
  int n_bytes = n * k * sizeof(usfixn64);
  usfixn64 * y = (usfixn64*) malloc (n_bytes);
  memcpy (y, x, n_bytes);

//  srand (time (NULL));
  srand (0);
  int * sn_vector = (int*) malloc (n * sizeof(int));
  for (int i = 0; i < n; i++)
    sn_vector[i] = rand () % k;

  cpu_timer t_mult_pow_r;
  memcpy (x, y, n_bytes);

  /////////////////////////////////////
  mpz_t *y_zz = (mpz_t*) malloc (n * sizeof(mpz_t));
  for (int i = 0; i < n; i++)
    {
      mpz_init (y_zz[i]);
    }

  cpu_timer t_convert_in, t_convert_out;

  /////////////////////////////////////
  timer_record_start (&t_convert_in);
  for (int i = 0; i < n; i++)
    u64vector_to_radix_based_bigint_gmp (y_zz[i], &y[i * k], radix, k);
  timer_record_stop (&t_convert_in);
  /////////////////////////////////////

  mpz_t p_zz;
  mpz_init (p_zz);
  compute_srgfn_p_gmp (p_zz, radix, k);

  timer_record_start (&t_mult_pow_r);
  for (int i = 0; i < n; i++)
    mult_pow_R_gmp (y_zz[i], sn_vector[i], k, radix, p_zz);
  timer_record_stop (&t_mult_pow_r);
  timer_get_elapsed_time (&t_mult_pow_r, "t_mult_pow_r_gmp", 1);
  timer_get_elapsed_time (&t_mult_pow_r, "t_mult_pow_r_avg_gmp", n);

  /////////////////////////////////////

  timer_record_start (&t_convert_out);
  for (int i = 0; i < n; i++)
    bigint_to_u64vector_radixbased_gmp (&y[i * k], y_zz[i], radix, k);
  timer_record_stop (&t_convert_out);

  timer_get_elapsed_time (&t_convert_in, NULL, 1);
  timer_get_elapsed_time (&t_convert_out, NULL, 1);

  /////////////////////////////////////

  for (int i = 0; i < n; i++)
    mpz_clear (y_zz[i]);

  mpz_clear (p_zz);
  free (y_zz);

  /////////////////////////////////////
  timer_record_start (&t_mult_pow_r);
  for (int i = 0; i < n; i++)
    mult_pow_R (&x[i * k], sn_vector[i], k, radix, 0);

  timer_record_stop (&t_mult_pow_r);

  /////////////////////////////////////
  float t_gmp_based = 0;
  //    t_gmp_based+=t_convert_in.elapsed_time;
  //    t_gmp_based+=t_convert_out.elapsed_time;
  timer_get_elapsed_time (&t_mult_pow_r, NULL, 1);
  t_gmp_based += t_mult_pow_r.elapsed_time;

  timer_print_time (t_gmp_based, "t_mult_pow_r_gmp_based");
  timer_print_time (t_gmp_based / (1.0 * n), "t_mult_pow_r_avg_gmp_based");

  timer_get_elapsed_time (&t_mult_pow_r, "t_mult_pow_r_gfpf", 1);
  timer_get_elapsed_time (&t_mult_pow_r, "t_mult_pow_r_avg_gfpf", n);

  /////////////////////////////////////
  int verification_result = EXIT_SUCCESS;
#if VERIFICATION_ENABLED
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
	if (x[i * k + j] != y[i * k + j])
	  {
	    printf ("-ERROR: mismatch at [%d, %d]\n", i, j);
	    verification_result = EXIT_FAILURE;
	    break;
	  }
      if (verification_result == EXIT_FAILURE)
	break;
    }
  print_verification_msg (__func__, verification_result);
#endif
  printf ("%s", longline);
/////////////////////////////////////
  free (y);
  free (sn_vector);
}

/**************************************/

void
test_swap_big_elements (usfixn64* x, int K, int e, usfixn64 radix)
{

  int n = (int) (pow (K, e));
  int k = K / 2;
  int n_bytes = n * k * sizeof(usfixn64);
  usfixn64 * y = (usfixn64*) malloc (n_bytes);
  memcpy (y, x, n_bytes);
  printf ("[testing DFT2_big_elements K=%d (k=%d), e=%d, radix=%llu] ...\n", K,
	  k, e, radix);

  srand (time (NULL));

  cpu_timer t_swap;

  for (int i = 0; i < n; i += 2)
    swap_big_elements (&x[i * k], &x[(i + 1) * k], k);

  memcpy (x, y, n_bytes);
  /////////////////////////////////////
  timer_record_start (&t_swap);
  for (int i = 0; i < n; i += 2)
    swap_big_elements (&x[i * k], &x[(i + 1) * k], k);
  timer_record_stop (&t_swap);
  timer_get_elapsed_time (&t_swap, "t_swap", 1);
  /////////////////////////////////////
  timer_record_start (&t_swap);
  for (int i = 0; i < n; i += 2)
    swap_big_elements_v0 (&y[i * k], &y[(i + 1) * k], k);

  timer_record_stop (&t_swap);
  timer_get_elapsed_time (&t_swap, "t_swap_old", 1);
  /////////////////////////////////////
  int verification_result = EXIT_SUCCESS;
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
	if (x[i * k + j] != y[i * k + j])
	  {
	    printf ("-ERROR: mismatch at [%d, %d]\n", i, j);
	    verification_result = EXIT_FAILURE;
	    break;
	  }
      if (verification_result == EXIT_FAILURE)
	break;
    }
  print_verification_msg (__func__, verification_result);
//  printf("%s", shortline);
/////////////////////////////////////
  free (y);

}

/**************************************/

#endif
