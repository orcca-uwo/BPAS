#ifndef GFPF_ARITHMETIC_DECL_H_
#define GFPF_ARITHMETIC_DECL_H_

/**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "gfpf_types.h"
#include "gfpf_gmp_tools_decl.h"
#include "gfpf_gmp_fft_decl.h"
#include "small_prime_field_fft_decl.h"
#include "cpu_timer.h"

/**************************************/

#ifndef CONVOLUTION_CACHE_SIZE
#define CONVOLUTION_CACHE_SIZE 4
#endif
#define CYCLICSHIFT_CACHE_SIZE CONVOLUTION_CACHE_SIZE

//#define BEGIN_PROFILE CALLGRIND_START_INSTRUMENTATION;
//#define END_PROFILE CALLGRIND_STOP_INSTRUMENTATION; CALLGRIND_DUMP_STATS;

/**************************************/
extern montgomery_triple * global_P1, *global_P2;
extern usfixn64 *global_conv_omega1, *global_conv_omega2;
extern usfixn64 *global_conv_omega1_inv, *global_conv_omega2_inv;
extern usfixn64 *global_conv_omega1_inv, *global_conv_omega2_inv;
extern usfixn64 *global_conv_theta1, *global_conv_theta2;
extern usfixn64 *global_conv_theta1_inv, *global_conv_theta2_inv;
extern usfixn64 global_K1_inv, global_K2_inv;

/**************************************/

#define ARITHMETIC_CACHE_SIZE 1

#ifndef VERIFICATION_ENABLED
#define VERIFICATION_ENABLED 1
#endif

/**************************************/
#ifndef LOOP_UNROLLING_ENABLED
#define LOOP_UNROLLING_ENABLED
#endif

/**************************************/
#ifndef PROFILING_ENABLED
#define PROFILING_ENABLED 0
#endif
static int global_profiling_ongoing = 0;

/**************************************/

/**
 * the following constants are used for montgomery
 * multiplication, addition, sub.
 *
 * R=(2**64)
 * prime = p
 * prime_inv = (-1/p) mod R
 * prime * prime_inv == -1 mod R
 * R1 = (2**128) mod p
 */

/**************************************/
//#define MAX_TMP_ARRAY_SIZE 64
/**************************************/
 void static
inc_profiling_counter (int *counter)
  {
    if (global_profiling_ongoing == 0)
    *counter++;
  }

 static int n_dft2_called = 0;
 static int n_swap_called = 0;
 static int n_DFTK_called = 0;
 static int n_mult_pow_R_called = 0;

 static int n_stride_permutation_called = 0;
 static int stride_permutation_called_dims[256][2];

 static int n_twiddle_called = 0;
 static int twiddle_called_dims[256][2];



/**************************************/

typedef union
{
  struct
  {
    usfixn64 m1_mont;
    usfixn64 m2_mont;

    usfixn64 p1;
    usfixn64 p2;

    usfixn64 p1_inv_m;
    usfixn64 p1_inv_q;
    usfixn64 p2_inv_m;
    usfixn64 p2_inv_q;

    //__int128 p1p2=m+q.u64;
    usfixn64 p1p2_q;
    usfixn64 p1p2_m;

    //// CRITICAL: THIS IS USED FOR NORMALIZING THE RESULT!
    //A=(p1p2-1)/2=a0+a1.u64
    //B=(p1p2-1)  =b0+b1.u64
    usfixn64 p1p2_a0;
    usfixn64 p1p2_a1;
    usfixn64 p1p2_b0;
    usfixn64 p1p2_b1;
  };
  usfixn64 memory[16];
} crt_u192_data;

/**************************************/

typedef union
{
  struct
  {

    //u128/radix=r_inv_0 + r_inv_1.u64
    usfixn64 r_inv_0;
    usfixn64 r_inv_1;
    usfixn64 radix;

    //u64=mb+qb.radix
    usfixn64 qb;
    usfixn64 mb;
  };
  usfixn64 memory[5];
} lhc_u192_data;

/**************************************/

extern crt_u192_data t_crt_data_global;
extern crt_u192_data *t_crt_data_global_ptr;    //=&t_crt_data_global;

extern lhc_u192_data t_lhc_data_global;
extern lhc_u192_data *t_lhc_data_global_ptr;    //=&t_crt_data_global;

/**************************************/
//CHECKED.
int
verify_mult_u64_u128_gmp (const usfixn64 a, const usfixn64 b0,
			  const usfixn64 b1, const usfixn64 s0,
			  const usfixn64 s1);

/**************************************/

//CHECKED.
int
verify_sub_u128_u128_gmp (const usfixn64 a0, const usfixn64 a1,
			  const usfixn64 b0, const usfixn64 b1,
			  const usfixn64 s0, const usfixn64 s1);

/**************************************/

//CHECKED.
int
verify_add_u128_u128_gmp (const usfixn64 a0, const usfixn64 a1,
			  const usfixn64 b0, const usfixn64 b1,
			  const usfixn64 s0, const usfixn64 s1,
			  const usfixn64 s2);

/**************************************/

//CHECKED.
void
mult_u64_u64 (const usfixn64 *a, const usfixn64 *b, usfixn64 *s0_out,
	      usfixn64 *s1_out);

/**************************************/

//CHECKED.
//compute a*(b0+b1.u64);
//results are correct only if the whole product
//is less than u128.
void
mult_u64_u128 (const usfixn64 *a, const usfixn64 *b0, const usfixn64 *b1,
	       usfixn64* s0, usfixn64 *s1);

/**************************************/

//CHECKED.
//(x0,x1)=(x0,x1)-(y0,y1);
void
sub_u128_u128 (const usfixn64 *y0, const usfixn64 *y1, usfixn64* x0,
	       usfixn64 *x1);

/**************************************/

//CHECKED.
// computes (a_u64)*(b0_u64, b1_u64)
// a=a_u64;
// b=b0_u64 + (b1_u64)*u64
void
u64_mod_u64 (usfixn64 *a, const usfixn64 n);

/**************************************/

//CHECKED.
//inv_p_mod_u128: returns [p_inv_m, p_inv_q];
// p_inv = p_inv_m + p_inv_q.u64  = int(u128/p);
void
inv_p_mod_u128 (const usfixn64 * p, usfixn64 * p_inv_m, usfixn64 * p_inv_q);

/**************************************/

//CHECKED.
// - compute [a1m2p1+a2m1p2] and return [u0,u1,u2]
int
verify_crt_mult_sub_u192_with_reduction_gmp (const usfixn64 a1,
					     const usfixn64 a2,
					     crt_u192_data data,
					     const usfixn64 s0,
					     const usfixn64 s1);

/**************************************/

//CHECKED.
int
verify_div_by_const_R_gmp (const usfixn64 *x0_u64, const usfixn64 *x1_u64,
			   const usfixn64 * q, const usfixn64 *m, usfixn64 R);

/**************************************/

//CHECKED.
void
div_by_const_R (const usfixn64 *x0_u64, const usfixn64 *x1_u64,
		const usfixn64 *r0, const usfixn64 *r1, usfixn64 *q,
		usfixn64 *m, usfixn64 R);

/**************************************/

//CHECKED.
//div_values = [x0,x1][rinv_0,rinv_1,R]
//const usfixn64 *x0_u64, const usfixn64 *x1_u64,
//		const usfixn64 *r0, const usfixn64 *r1, usfixn64 *q,
//		usfixn64 *m, usfixn64 R
void
div_by_const_R_ptr (usfixn64* div_values);

/**************************************/

//CHECKED.
//div_values = [x0,x1][rinv_0,rinv_1,R]
//const usfixn64 *x0_u64, const usfixn64 *x1_u64,
//		const usfixn64 *r0, const usfixn64 *r1, usfixn64 *q,
//		usfixn64 *m, usfixn64 R
void
div_by_const_R_ptr_single_digit (usfixn64* div_values);

/**************************************/

//CHECKED.
//lhc_0+=lhc_r mod r;
void
add_lhc_ptr (usfixn64 *__restrict__ lhc_0, const usfixn64* __restrict__ lhc_1,
	     const usfixn64* __restrict__ lhc_2, usfixn64 * __restrict__ R);

/**************************************/

//CHECKED.
void
crt_mult_sub_u192_with_reduction (usfixn64 *a1, usfixn64 *a2, int k,
				  const crt_u192_data * t_crt_data,
				  usfixn32 * sign_u32);

/**************************************/

//CHECKED.
// void  __attribute__((always_inline))
void
lhc_by_R_u128_ptr (usfixn64 * l_vec, usfixn64 * h_vec, usfixn64 * c_vec,
		   usfixn32 * sign_u32, int k);

/**************************************/

//usfixn64 aux_space[256];
//CHECKED.
void
oneShiftRight (usfixn64 * xs, int k);

/**************************************/

//CHECKED.
void
twoShiftRight (usfixn64 * xs, int k);

/**************************************/

void
normalize_positive_negative_vector (usfixn64 *x, const int k, const usfixn64 r);

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
//inline void  //__attribute__((always_inline))
void
addition_hybrid_big_elements (usfixn64 * __restrict__ x,
			      usfixn64 *__restrict__ y,
			      usfixn64 *__restrict__ z, const int k,
			      const usfixn64 r, usfixn64 *c_out,
			      usfixn64*c_negative_out);

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void
addition_big_element_with_single_digit (usfixn64 * x, usfixn64 digit, int plus,
					const int k, const usfixn64 r);

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void
addition_big_elements_v0 (usfixn64 * x, usfixn64 *y, const int k,
			  const usfixn64 r);

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void
addition_big_elements (usfixn64 * x, usfixn64 *y, const int k,
		       const usfixn64 r);

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void
subtraction_big_elements_v0 (usfixn64 *x, usfixn64 *y, const int k,
			     const usfixn64 r);

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void
subtraction_big_elements (usfixn64 *x, usfixn64 *y, const int k,
			  const usfixn64 r);

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void
subtraction_from_zero_big_elements (usfixn64 *y, const int k, const usfixn64 r);

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void
subtraction_upper_minus_lower_big_elements (usfixn64 *x, const int k,
					    const usfixn64 r, int sn);

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void
addition_subtraction_big_elements (usfixn64 * x, usfixn64 *y, const int k,
				   const usfixn64 r);

/**************************************/

void
shift_right_sn (usfixn64*x, int k, usfixn64 r, int s);

/**************************************/

//CHECKED.
//x=x*(r^s)
void
mult_pow_R (usfixn64 *x, int s, const int k, const usfixn64 r,
	    int compute_inverse);

/**************************************/

//CHECKED.
//x=x*(r^s)
void
mult_pow_R_v0 (usfixn64 *x, int s, const int k, const usfixn64 r,
	       int compute_inverse);

/**************************************/

//CHECKED.
void
conv_k (usfixn64 *x, usfixn64 *y, const int k,
	const usfixn64 *base_case_pow_omega,
	const usfixn64 *base_case_pow_omega_inv,
	const usfixn64 *base_case_pow_theta,
	const usfixn64 *base_case_pow_theta_inv, const usfixn64 * K_inv,
	const montgomery_triple*P);

/**************************************/

void
get_proper_conv_p1p2 (usfixn64 *conv_p1, usfixn64 *conv_p2, srgfn_prime prime);

/**************************************/

//CHECKED.
void
init_gfpf_mult_data (crt_u192_data *t_crt_data, lhc_u192_data *t_lhc_data,
		     usfixn64 r, const usfixn64 p1, const usfixn64 p2);

/**************************************/

//CHECKED.
void
init_fft_based_bigint_mult (usfixn64 conv_p1, usfixn64 conv_p2, int k);

/**************************************/

//CHECKED.
void
clear_gfpf_mult_data (crt_u192_data *t_crt_data, lhc_u192_data *t_lhc_data);

/**************************************/

//CHECKED.
void
clear_fft_based_bigint_mult_data ();

/**************************************/

//CHECKED.
void
GFPFMultiplication (usfixn64 *x, const usfixn64 *y, int k,
		    const crt_u192_data *t_crt_data,
		    const lhc_u192_data *t_lhc_data);

/**************************************/

//CHECKED.
void
test_fft_based_arbitrary_mult (usfixn64 *x, usfixn64*y, int n,
			       srgfn_prime prime);

/**************************************/
#endif
