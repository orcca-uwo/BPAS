#ifndef GFPF_GMP_TOOLS_DECL_H_
#define GFPF_GMP_TOOLS_DECL_H_

#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "gfpf_types.h"

/**************************************/

#define KARATSUBA_GRAIN_FACTOR 2

/**************************************/

int
min_three (int a, int b, int c);

/**************************************/

void
mpz_set_u64 (mpz_t xn, usfixn64 x);

/**************************************/

usfixn64
mpz_get_u64 (mpz_t n);

/**************************************/

void
u32vector_to_bigint_gmp (mpz_t bigint, usfixn32 *vector, int input_vector_size);

/**************************************/

///* bigint = (vector) in base=2^32
// * elements of vector are digits of a big integer
// * in base 2^32
// */
//void
//u64vector_to_bigint_gmp (mpz_t & bigint, usfixn32 *vector,
//			 int input_vector_size)
//;
/**************************************/

void
u32vector_to_radix_based_bigint_gmp (mpz_t bigint, usfixn32 *vector,
				     usfixn32 radix, int input_vector_size);
/**************************************/

void
u64vector_to_radix_based_bigint_gmp (mpz_t bigint, usfixn64 *vector,
				     usfixn64 radix, int input_vector_size);
/**************************************/

void
u32vector_to_bigint_gmp_via_str (mpz_t bigint, usfixn32 *vector,
				 int input_vector_size);
/**************************************/
//void
//u64vector_to_bigint_gmp (mpz_t & bigint, usfixn64 *vector,
//			 int input_vector_size)
//;
/**************************************/

void
bigint_to_u32vector_gmp (usfixn32 *vector, mpz_t bigint,
			 int max_input_vector_size);
/**************************************/

void
bigint_to_u32vector_radixbased_gmp (usfixn32 *vector, const mpz_t bigint,
				    usfixn32 radix, int max_input_vector_size);

/**************************************/

void
bigint_to_u64vector_radixbased_gmp (usfixn64 *vector, const mpz_t bigint,
				    usfixn64 radix, int max_input_vector_size);

/**************************************/

void
bigint_to_u64vector_gmp (usfixn64 *vector, mpz_t bigint,
			 int max_input_vector_size);

/**************************************/

/* input:
 * a linear 2d residue list residue_list [n_instances][n_primes]
 * a linear 2d prime_list prime_list [n_instances][n_primes]
 * number of primes n_primes
 * number of instances n_instances
 * function:
 * residue_list[i][j]= ((0xFFFFFFFF)*n_primes)% prime_list[i][j]
 */
void
cpu_generate_big_0xff_residue_list_vector_u32 (usfixn32* residue_list,
					       usfixn32* prime_list,
					       int n_primes, int n_instances);

/**************************************/

/* input:
 * a linear 2d residue list residue_list [n_instances][n_primes]
 * a linear 2d prime_list prime_list [n_instances][n_primes]
 * number of primes n_primes
 * number of instances n_instances
 * function:
 * residue_list[i][j]= ((0xFFFFFFFF)*n_primes)% prime_list[i][j]
 */
void
cpu_generate_single_0xff_residue_list_vector_u32 (usfixn32* residue_list,
						  usfixn32* prime_list,
						  int n_primes,
						  int n_instances);

/**************************************/

/* this function uses GMP to convert the result of
 * CRT to a bigint. Then, matches the residue of that
 * bigint for each of the primes in the prime_list.
 * Finally, it matches the result with the input residue
 * list of the CRT function.
 *
 * If successful, returns 0;
 * otherwise, returns -1;
 */

int
verify_crt_u32vector_gmp (usfixn32* crt_result, usfixn32* residue_list,
			  usfixn32* prime_list, int n_primes, int verbose);

/**************************************/

void
gmp_print_mpz (mpz_t bignum, const char* name, int base);

/**************************************/

void
gmp_bigint_to_radix_based_vec_u32 (usfixn32 *vec, mpz_t bigint, usfixn32 r,
				   int k);

/**************************************/

//compute p where p=r^k+1;
void
gmp_compute_srgfn (mpz_t p, usfixn32 r, int k);

/**************************************/

//compute p where p=r^k+1;
void
gmp_compute_srgfn_u64 (mpz_t p, usfixn64 r, int k);

/**************************************/

// computes (x_vec+y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_add_modp_u32 (usfixn32* x_vector, usfixn32* y_vector, mpz_t p, int k,
		  usfixn64 r);

/**************************************/

// computes (x_vec-y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_sub_modp_u32 (usfixn32* x_vector, usfixn32* y_vector, mpz_t p, int k,
		  usfixn32 r);

/**************************************/

/**************************************/
// computes (x_vec*y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_mult_modp_u32 (usfixn32* x_vector, usfixn32* y_vector, mpz_t p, int k,
		   usfixn32 r);

/**************************************/

// computes (x_vec*y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_mult_pow_r_modp_u32 (usfixn32* x_vector, usfixn32 sn, mpz_t p, int k,
			 usfixn32 r);

/**************************************/

void
gmp_mult_modp_u32_specific (usfixn32* x_vector, usfixn32* y_vector, usfixn32 r,
			    int k, mpz_t p);

/**************************************/

void
print_big_field_element_u32 (usfixn32* bigFieldElement,
			     usfixn32 coefficientSize);

/**************************************/
//
//int
//compute_srgfn_p_gmp (mpz_t p, mpz_t r, int coefficient_size)
//;
/**************************************/

void
bigint_to_bigprime_field_element_u32_gmp (usfixn32 *bigFieldElement,
					  mpz_t bigint, mpz_t p, mpz_t radix,
					  usfixn32 coefficientSize);

/**************************************/

void
vector_to_bigint_radix_based_u32 (mpz_t bigint, usfixn32* vector, mpz_t radix,
				  int coefficient_size);

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_plain_poly_mult_u32 (usfixn32* host_vector_x, usfixn32* host_vector_y,
			 int n, usfixn32 radix, int coefficient_size);

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_plain_poly_mult_u64 (usfixn64* host_vector_x, usfixn64* host_vector_y,
			 int n, usfixn64 radix, int coefficient_size);

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_karatsuba_poly_mult_u32 (usfixn32* host_vector_x, usfixn32* host_vector_y,
			     int n, usfixn32 radix, int coefficient_size);

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_karatsuba_poly_mult_u64 (usfixn64* host_vector_x, usfixn64* host_vector_y,
			     int n, usfixn64 radix, int coefficient_size);

/**************************************/

void
gmp_inv_mod_p (mpz_t inv, mpz_t x, mpz_t p);

/**************************************/

void
gmp_inv_mod_p_u64 (usfixn64 * inv, const usfixn64 x, const usfixn64 p);

/**************************************/

int
gmp_inv_radix_based_bigint_modp_u64 (usfixn64* result_vector,
				     usfixn64* x_vector, usfixn64 radix, int k,
				     mpz_t p);

/**************************************/

int
gmp_mult_radix_based_bigint_by_scalar_modp_u64 (usfixn64* result_vector,
						const mpz_t scalar_zz,
						usfixn64 radix, int k, mpz_t p);

/**************************************/

void
compute_srgfn_p_gmp (mpz_t p, usfixn64 radix, int coefficient_size);

/**************************************/

int
compute_nth_root_of_unity_for_small_prime (usfixn64 prime, usfixn64* omega,
					   int n);

/**************************************/

int
compute_nth_root_of_unity_for_srgfn (srgfn_prime prime, usfixn64* omega, int n);

/**************************************/

#endif //end of GFPF_GMP_TOOLS_H_
