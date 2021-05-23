//#ifndef SMALL_PRIME_FIELD_FFT_DECL_H_
//#define SMALL_PRIME_FIELD_FFT_DECL_H_
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <gmp.h>
////#include <unistd.h>
//
//#include <math.h>
//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>
////#include <immintrin.h>
//
//#include "gfpf_types.h"
//#include "gfpf_gmp_tools_decl.h"
//#include "fourier_primes_u64.h"
//#include "cpu_timer.h"
//
/////////////////////////////////////////////////////////////
//#ifndef ENABLE_CILK
//#define ENABLE_CILK 1
//#endif
//
//#if ENABLE_CILK == 0
//#define cilk_for for
//#endif
//
/////////////////////////////////////////////////////////////
//
//#ifndef VERIFICATION_ENABLED
//#define VERIFICATION_ENABLED 1
//#endif
//
/////////////////////////////////////////////////////////////
//
//#ifndef LOOP_UNROLLING_ENABLED
//#define LOOP_UNROLLING_ENABLED 1
//#endif
//
/////////////////////////////////////////////////////////////
//
//#ifndef PROFILING_ENABLED
//#define PROFILING_ENABLED 0
//#endif
//
/////////////////////////////////////////////////////////////
//
//#ifndef TWIDDLE_CACHE_SIZE
//#define TWIDDLE_CACHE_SIZE 32
//#endif
//
/////////////////////////////////////////////////////////////
//
//#define SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE 64
//#define SMALL_FFT_OPTIMAL_K 32
//#define SMALL_DFT_GRAIN_FACTOR 2
//#define SMALL_PERMUTATION_GRAIN_FACTOR 8
//#define SMALL_FFT_TWIDDLE_UNROLL_FACTOR 4
//#define SQURE_PERMUTATION_BLOCK_SIZE 4
//
/////////////////////////////////////////////////////////////
//
//#ifndef SMALL_FFT_STEPS_TIMING_ENABLED
//#define SMALL_FFT_STEPS_TIMING_ENABLED 1
//#endif
//
/////////////////////////////////////////////////////////////
//static int supported_small_prime_DFT_K[7] =
//  { 4, 8, 16, 32, 64, 128, 256 };
//
/////////////////////////////////////////////////////////////
//
//typedef struct montgomery_triple
//{
//  usfixn64 p;
//  usfixn64 p_inv;
//  usfixn64 R2_mod_p;
//} montgomery_triple;
//
/////////////////////////////////////////////////////////////
//
//void
//init_montgomery_triple (montgomery_triple*dest, usfixn64 prime);
//
/////////////////////////////////////////////////////////////
//
//usfixn64
//AddModSpe (usfixn64 a, usfixn64 b, const usfixn64 MY_PRIME);
//
/////////////////////////////////////////////////////////////
//
//usfixn64
//SubModSpe (const usfixn64 a, const usfixn64 b, const usfixn64 MY_PRIME);
//
/////////////////////////////////////////////////////////////
//
//usfixn64 //__attribute__ ((noinline))
//mult_ab_mod_p (usfixn64 a, const usfixn64 b, const montgomery_triple * P);
//
/////////////////////////////////////////////////////////////
//
//void //__attribute__ ((noinline))
//mult_ab_mod_p_ptr (usfixn64 * __restrict__ a, const usfixn64 *__restrict__ b,
//		   const montgomery_triple * P);
//
/////////////////////////////////////////////////////////////
////void //__attribute__ ((noinline))
////mult_ab_mod_p_vector(usfixn64 *a, const usfixn64 *b, int n, const montgomery_triple * P);
/////////////////////////////////////////////////////////////
////void //__attribute__ ((noinline))
////mult_ab_mod_p_vector_by_scalar(usfixn64 *a, const usfixn64 b, int n, const montgomery_triple * P);
/////////////////////////////////////////////////////////////
//
////void //__attribute__ ((noinline))
////mult_ab_mod_p_batch (usfixn64 *a0_in, usfixn64 *a1_in, const usfixn64 *b,
////		     const montgomery_triple * P);
///*************************************/
//
////void //__attribute__ ((noinline))
////mult_ab_mod_p_batch4 (usfixn64 *a0, usfixn64 *a1, usfixn64 *a2, usfixn64 *a3,
////		      const usfixn64 *b, const montgomery_triple * P);
///*************************************/
//
////void //__attribute__ ((noinline))
////mult_ab_mod_p_batch8 (usfixn64 *a0, usfixn64 *a1, usfixn64 *a2, usfixn64 *a3,
////		      usfixn64 *a4, usfixn64 *a5, usfixn64 *a6, usfixn64 *a7,
////		      const usfixn64 *b, const montgomery_triple * P);
///*************************************/
////a*b mod n;
//usfixn64
//convertIn_GLOBAL (const usfixn64 a, const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
////a*b mod n;
//void
//convertIn_GLOBAL_ptr (usfixn64 *a, const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//usfixn64
//convertOut_GLOBAL (const usfixn64 r, const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//convertOut_GLOBAL_ptr (usfixn64 * r, const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//usfixn64
//exponent_mod_p_SPF (usfixn64 a, usfixn64 n, const montgomery_triple * P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_2_double_simple (usfixn64* a, const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_2_double_plain (usfixn64* a0, usfixn64* a1, usfixn64* a2, usfixn64* a3,
//		    const usfixn64 prime);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_2_double_array (usfixn64* a0, usfixn64* a1, usfixn64* a2, usfixn64* a3,
//		    const usfixn64 prime);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_2_double (usfixn64* __restrict__ a0, usfixn64* __restrict__ a1,
//	      usfixn64* __restrict__ a2, usfixn64* __restrict__ a3,
//	      const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
////vectorized_avx2
//void
//DFT_2_double_avx2 (usfixn64* a0, usfixn64* a1, usfixn64* a2, usfixn64* a3,
//		   const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//swap (usfixn64* a, usfixn64* b);
//
/////////////////////////////////////////////////////////////
//
////pow_omega is already allocated and is of size n*sizeof(usfixn64);
//void
//precompute_pow_omega (usfixn64* pow_omega, const usfixn64 omega, int n,
//		      const montgomery_triple * P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_4 (int n, usfixn64* x, const usfixn64 *base_case_omega_vec,
//       const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_8 (int n, usfixn64* x, const usfixn64 *base_case_omega_vec,
//       const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_16 (int n, usfixn64* x, const usfixn64 * base_case_omega_vec,
//	const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_32 (int n, usfixn64* __restrict__ x,
//	const usfixn64 * __restrict__ base_case_omega_vec,
//	const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_64 (int n, usfixn64* __restrict__ x,
//	const usfixn64 * __restrict__ base_case_omega_vec,
//	const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_128 (int n, usfixn64* __restrict__ x,
//	 const usfixn64 * __restrict__ base_case_omega_vec,
//	 const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_256 (int n, usfixn64* __restrict__ x,
//	 const usfixn64 * __restrict__ base_case_omega_vec,
//	 const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
////
////void
////DFT_128 (int n, usfixn64* x, const usfixn64 * base_case_omega_vec,
////	const montgomery_triple* P)
////;
////}
////
///////////////////////////////////////////////////////////////
////
////void
////DFT_256 (int n, usfixn64* x, const usfixn64 * base_case_omega_vec,
////	const montgomery_triple* P)
////;
////}
////
/////////////////////////////////////////////////////////////
//
//void
//L_K_K (usfixn64*x, int K);
//
/////////////////////////////////////////////////////////////
//
////L_J_K = transpose {K}{J}
//void
//stride_permutation_small_elements_L_J_K (int n_permutations, usfixn64*x,
//					 const int J, const int K,
//					 usfixn64* aux_memory);
/////////////////////////////////////////////////////////////
//
////L_K_J = transpose {J}{K}
//void
//stride_permutation_small_elements_L_K_J (int n_permutations, usfixn64*x,
//					 const int K, const int J,
//					 usfixn64* aux_memory);
//
/////////////////////////////////////////////////////////////
//
////float total_memcpy=0;
////VERIFIED.
////compute L_{m}{n} for a vector of spf elements
////equivalently, computing transpose {n}{m}
//void
//stride_permutation_small_elements (usfixn64*x, int n_permutations, int m, int n,
//				   usfixn64* aux_memory);
//
/////////////////////////////////////////////////////////////
//
////VERIFIED.
////precompute powers of omega up to (omega)^{K^b}
////to be used at level b+1
////usfixn64* precomputed_pow_omega_vec should be be already allocated!
////just pass a pointer.
////ToDo: this is a prefix computation for mult -> implement it that way and in parallel.
//void
//precompute_sequential_powers_of_omega_small_elements (
//    usfixn64* precomputed_pow_omega_vec, int K, int b, usfixn64 omega,
//    const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
////VERIFIED.
////computing twiddle T{m}{n} == D{n}{m}
////computing D_{K}{J=K^b}
////part of twiddle computation is done by quotient powers;
////the rest is done by a direct power.
//void
//twiddle_small_elements (usfixn64* vector, const int n_permutations, const int K,
//			const int J, const usfixn64* precomputed_pow_omega_vec,
//			const montgomery_triple* P);
/////////////////////////////////////////////////////////////
//
//usfixn64
//pow_K_b (const int K_log2, const int b);
/////////////////////////////////////////////////////////////
//
//// the precomputation of powers of omega: VERIFIED.
//void
//precompute_pow_omega_for_all_levels (usfixn64**precomputed_pow_omega_vec_2D,
//				     const int K, const int e,
//				     const usfixn64 omega,
//				     const montgomery_triple* P);
//
/////////////////////////////////////////////////////////////
//
////swap pointers x and y,
////switch the value of check-bit from 0 to 1 and vice versa.
//void
//swap_ptr (usfixn64**x, usfixn64**y, int * check_bit);
//
/////////////////////////////////////////////////////////////
//
//void
//DFT_general_small_elements (usfixn64* x, int K, int e, montgomery_triple * P,
//			    usfixn64 omega);
//
/////////////////////////////////////////////////////////////
//
//int
//convolution_mult_small_elements (usfixn64* x_data, usfixn64* y_data, int n,
//				 const montgomery_triple *P);
//
/////////////////////////////////////////////////////////////
//
//#endif
