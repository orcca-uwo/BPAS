#ifndef GFPF_GMP_TOOLS_H_
#define GFPF_GMP_TOOLS_H_
#include "gfpf_gmp_tools_decl.h"
#include "errno.h"
/**************************************/

int
min_three (int a, int b, int c)
{
  int min;
  if (a < b)
    {
      min = a;
    }
  else
    {
      min = b;
    }

  if (min < c)
    {
      return min;
    }
  else
    {
      return c;
    }
}

/**************************************/

void __inline__
mpz_set_u64 (mpz_t xn, usfixn64 x)
{
//  mpz_set_ui (xn, (usfixn32) (x >> 32));  //getting higher part
//  mpz_mul_2exp (xn, xn, 32);  //shifting higher part
//  mpz_add_ui (xn, xn, (usfixn32) (x)); //adding lower and higher
  mpz_set_ui (xn, x);  //getting higher part

}

/**************************************/

usfixn64 __inline__
mpz_get_u64 (mpz_t n)
{

//  mpz_t tmp;
//  mpz_init (tmp);
//  mpz_mod_2exp (tmp, n, 64);
//
//  usfixn32 l, h;
//  l = mpz_get_ui (tmp);
//  mpz_div_2exp (tmp, tmp, 32);
//  h = mpz_get_ui (tmp);
//  mpz_clear (tmp);
//
//  return (((usfixn64) h) << 32) + l;
  return mpz_get_ui (n);
}

/**************************************/

void
u32vector_to_bigint_gmp (mpz_t bigint, usfixn32 *vector, int input_vector_size)
{
  mpz_set_ui (bigint, 0);

//////////////////////////////////////
//  mpz_t pow_base;
//  mpz_t m;
//  mpz_init_set_ui (pow_base, 1);
//  mpz_init_set_ui (m, 0);

//  for (int i = 0; i < input_vector_size; i++)
//    {
//      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_mul_ui (pow_base, pow_base, U32_BASE);
//      mpz_add (bigint, bigint, m);
//    }

//  mpz_clear (pow_base);
//  mpz_clear (m);

//////////////////////////////////////

//Horner method for evaluation
  mpz_set_ui (bigint, vector[input_vector_size - 1]);
  for (int i = input_vector_size - 1; i > 0; i--)
    {
      mpz_mul_ui (bigint, bigint, U32_BASE);
      mpz_add_ui (bigint, bigint, vector[i - 1]);
//      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_mul_ui (pow_base, pow_base, U32_BASE);
//      mpz_add (bigint, bigint, m);
    }

}

/**************************************/

///* bigint = (vector) in base=2^32
// * elements of vector are digits of a big integer
// * in base 2^32
// */
//void
//u64vector_to_bigint_gmp (mpz_t & bigint, usfixn32 *vector,
//			 int input_vector_size)
//{
//  mpz_set_ui (bigint, 0);
//
//  mpz_t pow_base;
//  mpz_t m;
//  mpz_init_set_ui (pow_base, 1);
//  mpz_init_set_ui (m, 0);
//
//  mpz_t mpz_u64_base;
//  mpz_init (mpz_u64_base);
//  mpz_set_u64 (mpz_u64_base, U64_BASE);
//
//  mpz_t tmp;
//  mpz_init (tmp);
//
//  for (int i = 0; i < input_vector_size; i++)
//    {
////      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_set_u64 (tmp, vector[i]);
//      mpz_mul (m, pow_base, tmp);
////      mpz_mul_ui (pow_base, pow_base, U64_BASE);
//      mpz_mul (pow_base, pow_base, mpz_u64_base);
//      mpz_add (bigint, bigint, m);
//    }
//
//  mpz_clear (pow_base);
//  mpz_clear (m);
//
//  mpz_clear (mpz_u64_base);
//  mpz_clear (tmp);
//}
//
/**************************************/

void
u32vector_to_radix_based_bigint_gmp (mpz_t bigint, usfixn32 *vector,
				     usfixn32 radix, int input_vector_size)
{
  mpz_set_ui (bigint, 0);

  //////////////////////////////////////
//  mpz_t pow_base;
//  mpz_t m;
//  mpz_init_set_ui (pow_base, 1);
//  mpz_init_set_ui (m, 0);

//  for (int i = 0; i < input_vector_size; i++)
//    {
//      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_mul_ui (pow_base, pow_base, U32_BASE);
//      mpz_add (bigint, bigint, m);
//    }

//  mpz_clear (pow_base);
//  mpz_clear (m);

//////////////////////////////////////

//Horner method for evaluation
  mpz_set_ui (bigint, vector[input_vector_size - 1]);
  for (int i = input_vector_size - 1; i > 0; i--)
    {
      mpz_mul_ui (bigint, bigint, radix);
      mpz_add_ui (bigint, bigint, vector[i - 1]);
//      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_mul_ui (pow_base, pow_base, U32_BASE);
//      mpz_add (bigint, bigint, m);
    }

}

/**************************************/

void
u64vector_to_radix_based_bigint_gmp (mpz_t bigint, usfixn64 *vector,
				     usfixn64 radix, int input_vector_size)
{
//  mpz_init (bigint);
  mpz_set_ui (bigint, 0);

  mpz_t radix_mpz, tmp;
  mpz_inits (radix_mpz, tmp, NULL);

  //////////////////////////////////////
//  mpz_t pow_base;
//  mpz_t m;
//  mpz_init_set_ui (pow_base, 1);
//  mpz_init_set_ui (m, 0);

//  for (int i = 0; i < input_vector_size; i++)
//    {
//      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_mul_ui (pow_base, pow_base, U32_BASE);
//      mpz_add (bigint, bigint, m);
//    }

//  mpz_clear (pow_base);
//  mpz_clear (m);

//////////////////////////////////////

//Horner method for evaluation
//  mpz_set_ui (bigint, vector[input_vector_size - 1]);
  mpz_set_u64 (bigint, vector[input_vector_size - 1]);
  mpz_set_u64 (radix_mpz, radix);
  for (int i = input_vector_size - 1; i > 0; i--)
    {
//      mpz_mul_ui (bigint, bigint, radix);
//      mpz_add_ui (bigint, bigint, vector[i - 1]);

      mpz_mul (bigint, bigint, radix_mpz);
      mpz_set_u64 (tmp, vector[i - 1]);
      mpz_add (bigint, bigint, tmp);
//      mpz_mul_ui (m, pow_base, vector[i]);
//      mpz_mul_ui (pow_base, pow_base, U32_BASE);
//      mpz_add (bigint, bigint, m);
    }

  mpz_clears (radix_mpz, tmp, NULL);
}

/**************************************/

void
u32vector_to_bigint_gmp_via_str (mpz_t bigint, usfixn32 *vector,
				 int input_vector_size)
{

  /*  rewrite digits of bigint in hex format;
   * then write them to a large string
   */

  char * str = (char*) malloc (sizeof(char) * input_vector_size * 9);

  char hex_int[9];

  mpz_init_set_ui (bigint, 0);
  mpz_t t;
  mpz_t m;

  mpz_init_set_ui (t, 1);
  mpz_init_set_ui (m, 0);

  memset (str, 0x00, sizeof(char) * input_vector_size * 9);

//  for (int i = 0; i < input_vector_size; i++)
  for (int i = input_vector_size - 1; i >= 0; i--)
    {
      sprintf (hex_int, "%x", vector[i]);

//      prefixing the empty space with 0s
      for (int j = 0; j < 8 - strlen (hex_int); j++)
	strcat (str, "0");
      strcat (str, hex_int);

    }
//  base is 16 (reading hex digits)
  mpz_set_str (bigint, str, 16);

  mpz_clear (t);
  mpz_clear (m);
  free (str);

}

/**************************************/
//
//void
//u64vector_to_bigint_gmp (mpz_t & bigint, usfixn64 *vector,
//			 int input_vector_size)
//{
//  mpz_set_ui (bigint, 0);
//  mpz_t t;
//  mpz_t m;
//
//  mpz_init_set_ui (t, 1);
//  mpz_init_set_ui (m, 0);
//
//  mpz_t mpz_u64_base;
//  mpz_init (mpz_u64_base);
//  mpz_set_u64 (mpz_u64_base, U64_BASE);
//
//  mpz_t tmp;
//  mpz_init (tmp);
//
//  for (int i = 0; i < input_vector_size; i++)
//    {
////      mpz_mul_ui (m, t, vector[i]);
//      mpz_set_u64 (tmp, vector[i]);
//      mpz_mul (m, t, tmp);
////      mpz_mul_ui (t, t, U64_BASE);
//      mpz_mul (t, t, mpz_u64_base);
//      mpz_add (bigint, bigint, m);
//    }
//
//  mpz_clear (t);
//  mpz_clear (m);
//
//  mpz_clear (mpz_u64_base);
//  mpz_clear (tmp);
//
//}
/**************************************/

void
bigint_to_u32vector_gmp (usfixn32 *vector, mpz_t bigint,
			 int max_input_vector_size)
{

  mpz_t q, r;
  mpz_init_set_ui (q, 0);
  mpz_init_set_ui (r, 0);

  for (int i = 0; i < max_input_vector_size; i++)
    {
//      mpz_fdiv_qr (q, r, bigint, base);
      mpz_fdiv_qr_ui (q, r, bigint, U32_BASE);
      vector[i] = mpz_get_ui (r);
      mpz_set (bigint, q);
    }

  mpz_clear (q);
  mpz_clear (r);

}

/**************************************/

void
bigint_to_u32vector_radixbased_gmp (usfixn32 *vector, const mpz_t bigint,
				    usfixn32 radix, int max_input_vector_size)
{

  mpz_t q, r;
  mpz_init_set_ui (q, 0);
  mpz_init_set_ui (r, 0);

  mpz_set (q, bigint);

  int all_zero = 1;
  for (int i = 0; i < max_input_vector_size; i++)
    {
//      mpz_fdiv_qr (q, r, bigint, base);
      mpz_fdiv_qr_ui (q, r, q, radix);
      vector[i] = mpz_get_ui (r);
//      mpz_set (bigint, q);
      if (vector[i])
	all_zero = 0;
    }

  if (all_zero == 1)
    {
//      printf ("P-1 OCCURED IN CONVERSION\n");
      vector[max_input_vector_size - 1] = radix;
    }

  mpz_clear (q);
  mpz_clear (r);

//	mpz_t r;
////	  mpz_init_set_ui (q, 0);
//	  mpz_init (r);
//
//	  for (int i = 0; i < max_input_vector_size; i++)
//	    {
//	//      mpz_fdiv_qr (q, r, bigint, base);
//	      mpz_fdiv_qr_ui (bigint, r, bigint, radix);
//	      vector[i] = mpz_get_ui (r);
//	    }
//
////	  mpz_clear (q);
//	  mpz_clear (r);
}

/**************************************/

void
bigint_to_u64vector_radixbased_gmp (usfixn64 *vector, const mpz_t bigint,
				    usfixn64 radix, int max_input_vector_size)
{

  mpz_t q, r, radix_mpz;
  mpz_inits (q, r, radix_mpz, NULL);
  mpz_set_u64 (radix_mpz, radix);
  mpz_set (q, bigint);

  int all_zero = 1;
  for (int i = 0; i < max_input_vector_size; i++)
    {
//      mpz_fdiv_qr (q, r, bigint, base);
//      mpz_fdiv_qr_ui (q, r, q, radix);
      mpz_fdiv_qr (q, r, q, radix_mpz);
//      vector[i] = mpz_get_ui (r);
      vector[i] = mpz_get_u64 (r);
//      mpz_set (bigint, q);
      if (vector[i])
	all_zero = 0;
    }

//  if (all_zero == 1)
//    {
//      printf ("P-1 OCCURED IN CONVERSION\n");
//      vector[max_input_vector_size-1]=radix;
//      for(int i=0;i<max_input_vector_size;i++)
//	{
//	  printf("x[%d]=%lu\n", i, vector[i]);
//	}
//    }

  mpz_clears (q, r, radix_mpz, NULL);

}

/**************************************/

void
bigint_to_u64vector_gmp (usfixn64 *vector, mpz_t bigint,
			 int max_input_vector_size)
{

  mpz_t q, r;

  mpz_init_set_ui (q, 0);
  mpz_init_set_ui (r, 0);

  mpz_t mpz_u64_base;
  mpz_init (mpz_u64_base);
  mpz_set_u64 (mpz_u64_base, (usfixn64) (U64_MAX));
  mpz_add_ui (mpz_u64_base, mpz_u64_base, 1);

  for (int i = 0; i < max_input_vector_size; i++)
    {
//      mpz_fdiv_qr_ui (q, r, bigint, U64_BASE);
      mpz_fdiv_qr (q, r, bigint, mpz_u64_base);
//      vector[i] = mpz_get_ui (r);
      vector[i] = mpz_get_u64 (r);
      mpz_set (bigint, q);
    }

  mpz_clear (q);
  mpz_clear (r);
  mpz_clear (mpz_u64_base);
}

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
					       int n_primes, int n_instances)
{
  mpz_t big_n;
  mpz_init_set_ui (big_n, 0);

  mpz_t residue;
  mpz_init_set_ui (residue, 0);

  char * big_n_str = (char*) malloc (
      ((n_primes * (32 / 4) + 32) * sizeof(char)));
  //  sprintf(big_n_str,"0x");
  for (int i = 0; i < n_primes; i++)
    {
      strcat (big_n_str, "FFFFFFFF");
    }

  mpz_init_set_str (big_n, big_n_str, 16);

  //  for (int i = 0; i < n_primes; i++)
  //    {
  //      residue_list[i] = tmp_res_list[i];
  //      prime_list[i] = tmp_prime_list[i];
  //    }

  for (int j = 0; j < n_instances; j++)
    for (int i = 0; i < n_primes; i++)
      {
	mpz_tdiv_r_ui (residue, big_n, prime_list[j * n_primes + i]);
	residue_list[j * n_primes + i] = (mpz_get_ui (residue));
      }
  mpz_clear (big_n);
  mpz_clear (residue);

  free (big_n_str);
}

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
						  int n_primes, int n_instances)
{
  for (int j = 0; j < n_instances; j++)
    for (int i = 0; i < n_primes; i++)
      {
	residue_list[j * n_primes + i] = (0xFFFFFFFF)
	    % prime_list[j * n_primes + i];
      }

}

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
			  usfixn32* prime_list, int n_primes, int verbose)
{
  mpz_t bigint;
  mpz_init_set_ui (bigint, 0);
  u32vector_to_bigint_gmp (bigint, crt_result, n_primes);
  mpz_t res;
  mpz_init_set_ui (res, 0);
  for (int i = 0; i < n_primes; i++)
    {
      mpz_tdiv_r_ui (res, bigint, prime_list[i]);
      if (mpz_get_ui (res) != residue_list[i])
	{

	  printf ("mismatch in result @ i=%d \n", i);
	  printf ("gmp_res = %lu ,  res = %u \n ", mpz_get_ui (res),
		  residue_list[i]);
	  return -1;
	}
      else
	{
	  if (verbose == 1)
	    printf ("... PASS: res [%d] matches ... \n", i);
	}
    }
  if (verbose == 1)
    {
      printf ("PASS!\n");
//      printf (shortline);
    }
  return 0;
}

/**************************************/

void
gmp_print_mpz (mpz_t bignum, const char* name, int base)
{
  char * num_str = (char*) malloc (4096);
  mpz_get_str (num_str, base, bignum);
  printf ("[%s]=%s\n", name, num_str);

}

/**************************************/

void
gmp_bigint_to_radix_based_vec_u32 (usfixn32 *vec, mpz_t bigint, usfixn32 r,
				   int k)
{
  mpz_t tmp;
  mpz_init (tmp);

  for (int i = 0; i < k; i++)
    {
      mpz_mod_ui (tmp, bigint, r);
      mpz_div_ui (bigint, bigint, r);
      vec[i] = mpz_get_ui (tmp);
//      gmp_print_mpz(bigint, "bigint");
//      printf("vec[%d]=%u\n",i, vec[i]);
    }
  //corner case when bigint=r^k = r*(r^{k-1})
  if (mpz_cmp_ui(bigint, 1) == 0)
    {
      vec[k - 1] = r;
    }

  mpz_clear (tmp);
}

/**************************************/

//compute p where p=r^k+1;
void
gmp_compute_srgfn (mpz_t p, usfixn32 r, int k)
{
  mpz_set_ui (p, r);
  mpz_pow_ui (p, p, k);
  mpz_add_ui (p, p, 1);
}

/**************************************/

//compute p where p=r^k+1;
void
gmp_compute_srgfn_u64 (mpz_t p, usfixn64 r, int k)
{
  mpz_set_u64 (p, r);
  mpz_pow_ui (p, p, k);
  mpz_add_ui (p, p, 1);
}

/**************************************/

// computes (x_vec+y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_add_modp_u32 (usfixn32* x_vector, usfixn32* y_vector, mpz_t p, int k,
		  usfixn64 r)
{
  mpz_t x_bignum, y_bignum;
  mpz_init_set_ui (x_bignum, 0);
  mpz_init_set_ui (y_bignum, 0);

  mpz_t tmp;
  mpz_init_set_ui (tmp, 0);

  u32vector_to_radix_based_bigint_gmp (x_bignum, x_vector, r, k);
  u32vector_to_radix_based_bigint_gmp (y_bignum, y_vector, r, k);

  mpz_add (x_bignum, x_bignum, y_bignum);
  mpz_mod (x_bignum, x_bignum, p);

  gmp_bigint_to_radix_based_vec_u32 (x_vector, x_bignum, r, k);

  mpz_clear (x_bignum);
  mpz_clear (y_bignum);
  mpz_clear (tmp);
}

/**************************************/

// computes (x_vec-y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_sub_modp_u32 (usfixn32* x_vector, usfixn32* y_vector, mpz_t p, int k,
		  usfixn32 r)
{
  mpz_t x_bignum, y_bignum;
  mpz_init_set_ui (x_bignum, 0);
  mpz_init_set_ui (y_bignum, 0);

//  int k = COEFFICIENT_SIZE;
//  usfixn32 r = R;

  mpz_t tmp;
  mpz_init_set_ui (tmp, 0);

  u32vector_to_radix_based_bigint_gmp (x_bignum, x_vector, r, k);
  u32vector_to_radix_based_bigint_gmp (y_bignum, y_vector, r, k);

  mpz_sub (x_bignum, x_bignum, y_bignum);
  mpz_mod (x_bignum, x_bignum, p);

  gmp_bigint_to_radix_based_vec_u32 (x_vector, x_bignum, r, k);

  mpz_clear (x_bignum);
  mpz_clear (y_bignum);
  mpz_clear (tmp);
}

/**************************************/

/**************************************/
// computes (x_vec*y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_mult_modp_u32 (usfixn32* x_vector, usfixn32* y_vector, mpz_t p, int k,
		   usfixn32 r)
{
  mpz_t x_bignum, y_bignum;
  mpz_init_set_ui (x_bignum, 0);
  mpz_init_set_ui (y_bignum, 0);

//  int k = COEFFICIENT_SIZE;
//  usfixn32 r = R;

  mpz_t tmp;
  mpz_init_set_ui (tmp, 0);

  u32vector_to_radix_based_bigint_gmp (x_bignum, x_vector, r, k);
  u32vector_to_radix_based_bigint_gmp (y_bignum, y_vector, r, k);

  mpz_mul (x_bignum, x_bignum, y_bignum);
  mpz_mod (x_bignum, x_bignum, p);

  gmp_bigint_to_radix_based_vec_u32 (x_vector, x_bignum, r, k);

  mpz_clear (x_bignum);
  mpz_clear (y_bignum);
  mpz_clear (tmp);
}

/**************************************/

// computes (x_vec*y_vec)mod(p=r^k+1)
// store the result in x_vector
void
gmp_mult_pow_r_modp_u32 (usfixn32* x_vector, usfixn32 sn, mpz_t p, int k,
			 usfixn32 r)
{
//  int k = COEFFICIENT_SIZE;
//  usfixn32 r = R;
  mpz_t x_bignum, y_bignum;
  mpz_init_set_ui (x_bignum, 0);
  mpz_init_set_ui (y_bignum, 0);

  mpz_t tmp;
  mpz_init_set_ui (tmp, r);
  mpz_powm_ui (y_bignum, tmp, sn, p);

  u32vector_to_radix_based_bigint_gmp (x_bignum, x_vector, r, k);
//  gmp_print_mpz(x_bignum,"x");
//  gmp_print_mpz(y_bignum,"y");

  mpz_mul (x_bignum, x_bignum, y_bignum);
  mpz_mod (x_bignum, x_bignum, p);
//  gmp_print_mpz(x_bignum,"x");

  gmp_bigint_to_radix_based_vec_u32 (x_vector, x_bignum, r, k);

  mpz_clear (x_bignum);
  mpz_clear (y_bignum);
  mpz_clear (tmp);
}

/**************************************/

void
gmp_mult_modp_u32_specific (usfixn32* x_vector, usfixn32* y_vector, usfixn32 r,
			    int k, mpz_t p)
{
  mpz_t x_bignum, y_bignum;
  mpz_init_set_ui (x_bignum, 0);
  mpz_init_set_ui (y_bignum, 0);

  mpz_t tmp;
  mpz_init_set_ui (tmp, 0);

  u32vector_to_radix_based_bigint_gmp (x_bignum, x_vector, r, k);
  u32vector_to_radix_based_bigint_gmp (y_bignum, y_vector, r, k);

//  gmp_print_mpz(x_bignum, "x");
//  gmp_print_mpz(y_bignum, "y");

  mpz_mul (x_bignum, x_bignum, y_bignum);
  mpz_mod (x_bignum, x_bignum, p);

  gmp_bigint_to_radix_based_vec_u32 (x_vector, x_bignum, r, k);

  mpz_clear (x_bignum);
  mpz_clear (y_bignum);
  mpz_clear (tmp);
}

/**************************************/

void
print_big_field_element_u32 (usfixn32* bigFieldElement,
			     usfixn32 coefficientSize)
{
  for (int i = 0; i < coefficientSize; i++)
//    cout << "[" << i << "] = " << bigFieldElement[i] << endl;
    printf ("[%d] = %u\n", i, bigFieldElement[i]);
}

/**************************************/
//
//int
//compute_srgfn_p_gmp (mpz_t p, mpz_t r, int coefficient_size)
//{
//  mpz_t m;
//  mpz_init_set_ui (m, 0);
//
//  mpz_pow_ui (m, r, coefficient_size);
//  mpz_add_ui (p, m, 1);
//  mpz_clear (m);
//  return 0;
//}
/**************************************/

void
bigint_to_bigprime_field_element_u32_gmp (usfixn32 *bigFieldElement,
					  mpz_t bigint, mpz_t p, mpz_t radix,
					  usfixn32 coefficientSize)
{

  mpz_mod (bigint, bigint, p);
//  rem (bigInt, bigInt, prime);
//  ZZ m;
  mpz_t m;
  mpz_init_set_ui (m, 0);

  mpz_t q;
  mpz_init_set (q, bigint);

  for (int i = 0; i < coefficientSize; i++)
    {
//      rem (m, bigInt, radix);
////		cout<<"m " <<m<<endl;
//      bigFieldElement[i] = to_ulong (m);
////		cout<<bigFieldElement[i]<<endl;
//      div (bigInt, bigInt, radix);

//      printf("radix = %lu \n", radix);
//      mpz_mod (m, bigint, radix);
//      bigFieldElement[i] = mpz_get_ui (m);
//      printf ("m = %lu \n", mpz_get_ui (m));
//      mpz_tdiv_q (bigint, bigint, radix);

      mpz_mod (m, q, radix);
      bigFieldElement[i] = mpz_get_ui (m);
//      printf ("m = %lu \n", mpz_get_ui (m));
      mpz_tdiv_q (q, q, radix);
    }

  mpz_clear (m);
  mpz_clear (q);
}

/**************************************/

void
vector_to_bigint_radix_based_u32 (mpz_t bigint, usfixn32* vector, mpz_t radix,
				  int coefficient_size)
{
  mpz_t m, pow_radix;
//  mpz_init_set_ui(bigint,0);
  mpz_set_ui (bigint, 0);
  mpz_init_set_ui (pow_radix, 1);
  mpz_init_set_ui (m, 0);

  for (int i = 0; i < coefficient_size; i++)
    {
      mpz_mul_ui (m, pow_radix, vector[i]);
      mpz_add (bigint, bigint, m);
      mpz_mul (pow_radix, pow_radix, radix);
    }

  mpz_clear (m);
  mpz_clear (pow_radix);
}

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_plain_poly_mult_u32 (usfixn32* host_vector_x, usfixn32* host_vector_y,
			 int n, usfixn32 radix, int coefficient_size)
{
  int data_size = n * coefficient_size * sizeof(usfixn32);
  usfixn32 * host_vector_h = (usfixn32*) malloc (data_size);
  memset (host_vector_h, 0x00, data_size);

  mpz_t xn, yn, hn, tmp, p;
  mpz_init (xn);
  mpz_init (yn);
  mpz_init (hn);
  mpz_init (p);
  mpz_init (tmp);

  gmp_compute_srgfn (p, radix, coefficient_size);

  int K_minus_one = 2 * coefficient_size - 1;

  mpz_t *xn_array, *yn_array, *hn_array;

  xn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  yn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  hn_array = (mpz_t*) malloc ((n) * sizeof(mpz_t));
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (xn_array[i]);
      u32vector_to_radix_based_bigint_gmp (xn_array[i],
					   &host_vector_x[i * coefficient_size],
					   radix, coefficient_size);
    }
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (yn_array[i]);
      u32vector_to_radix_based_bigint_gmp (yn_array[i],
					   &host_vector_y[i * coefficient_size],
					   radix, coefficient_size);
    }
  for (int i = 0; i < n; i++)
    {
      mpz_init_set_ui (hn_array[i], 0);
    }

  for (int i = 0; i < n / 2; i++)
    {
      if ((i & K_minus_one) == 0)
	printf (".");

      for (int j = 0; j < n / 2; j++)
	{
	  mpz_mul (tmp, xn_array[i], yn_array[j]);
	  mpz_add (hn_array[i + j], hn_array[i + j], tmp);
	  mpz_mod (hn_array[i + j], hn_array[i + j], p);
	}
    }

  printf ("\n");

  for (int i = 0; i < n; i++)
    {
      gmp_bigint_to_radix_based_vec_u32 (&host_vector_h[(i) * coefficient_size],
					 hn_array[i], radix, coefficient_size);
    }

  memcpy (host_vector_x, host_vector_h, data_size);
  free (host_vector_h);
  mpz_clear (xn);
  mpz_clear (yn);
  mpz_clear (hn);
  mpz_clear (p);
  mpz_clear (tmp);

  for (int i = 0; i < n / 2; i++)
    mpz_clear (xn_array[i]);
  for (int i = 0; i < n / 2; i++)
    mpz_clear (yn_array[i]);
  for (int i = 0; i < n; i++)
    mpz_clear (hn_array[i]);
  free (xn_array);
  free (yn_array);
  free (hn_array);
}

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_plain_poly_mult_u64 (usfixn64* host_vector_x, usfixn64* host_vector_y,
			 int n, usfixn64 radix, int coefficient_size)
{
  int data_size = n * coefficient_size * sizeof(usfixn64);
  usfixn64 * host_vector_h = (usfixn64*) malloc (data_size);
  memset (host_vector_h, 0x00, data_size);

  mpz_t xn, yn, hn, tmp, p;
  mpz_init (xn);
  mpz_init (yn);
  mpz_init (hn);
  mpz_init (p);
  mpz_init (tmp);

  gmp_compute_srgfn_u64 (p, radix, coefficient_size);

  int K_minus_one = 2 * coefficient_size - 1;

  mpz_t *xn_array, *yn_array, *hn_array;

  xn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  yn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  hn_array = (mpz_t*) malloc ((n) * sizeof(mpz_t));
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (xn_array[i]);
      u64vector_to_radix_based_bigint_gmp (xn_array[i],
					   &host_vector_x[i * coefficient_size],
					   radix, coefficient_size);
    }
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (yn_array[i]);
      u64vector_to_radix_based_bigint_gmp (yn_array[i],
					   &host_vector_y[i * coefficient_size],
					   radix, coefficient_size);
    }
  for (int i = 0; i < n; i++)
    {
      mpz_init_set_ui (hn_array[i], 0);
    }

  for (int i = 0; i < n / 2; i++)
    {
      if ((i & K_minus_one) == 0)
	printf (".");

      for (int j = 0; j < n / 2; j++)
	{
//	  printf ("[i=%d,j=%d]\n", i, j);
//	  gmp_print_mpz (xn_array[i], "xn");
//	  gmp_print_mpz (yn_array[j], "yn");
	  mpz_mul (tmp, xn_array[i], yn_array[j]);
	  mpz_mod (tmp, tmp, p);
//	  gmp_print_mpz (tmp, "tmp");
	  mpz_add (hn_array[i + j], hn_array[i + j], tmp);
	  mpz_mod (hn_array[i + j], hn_array[i + j], p);
//	  printf ("\n%s", shortline);
	}
    }
  printf ("\n");

  for (int i = 0; i < n; i++)
    {
      bigint_to_u64vector_radixbased_gmp (
	  &host_vector_h[(i) * coefficient_size], hn_array[i], radix,
	  coefficient_size);
    }

  memcpy (host_vector_x, host_vector_h, data_size);
  free (host_vector_h);
  mpz_clear (xn);
  mpz_clear (yn);
  mpz_clear (hn);
  mpz_clear (p);
  mpz_clear (tmp);

  for (int i = 0; i < n / 2; i++)
    mpz_clear (xn_array[i]);
  for (int i = 0; i < n / 2; i++)
    mpz_clear (yn_array[i]);
  for (int i = 0; i < n; i++)
    mpz_clear (hn_array[i]);
  free (xn_array);
  free (yn_array);
  free (hn_array);
}

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_karatsuba_poly_mult_u32 (usfixn32* host_vector_x, usfixn32* host_vector_y,
			     int n, usfixn32 radix, int coefficient_size)
{
  int data_size = n * coefficient_size * sizeof(usfixn32);
  usfixn32 * host_vector_h = (usfixn32*) malloc (data_size);
  memset (host_vector_h, 0x00, data_size);

  mpz_t tmp, p;
  mpz_init (p);
  mpz_init (tmp);

  gmp_compute_srgfn (p, radix, coefficient_size);

  int K_minus_one = 2 * coefficient_size - 1;

  mpz_t *xn_array, *yn_array, *hn_array, *c1_array;

  xn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  yn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  c1_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  hn_array = (mpz_t*) malloc ((n) * sizeof(mpz_t));

  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (xn_array[i]);
      u32vector_to_radix_based_bigint_gmp (xn_array[i],
					   &host_vector_x[i * coefficient_size],
					   radix, coefficient_size);
    }
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (yn_array[i]);
      u32vector_to_radix_based_bigint_gmp (yn_array[i],
					   &host_vector_y[i * coefficient_size],
					   radix, coefficient_size);
    }
  for (int i = 0; i < n; i++)
    {
      mpz_init_set_ui (hn_array[i], 0);
    }
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init_set_ui (c1_array[i], 0);
    }
//  #a0b0
//  for i in range(n/2):
//	  for j in range(n/2):
//		  c0[i+j]+=a[i]*b[j]
//
//  #a1b1
//  for i in range(n/2):
//	  for j in range(n/2):
//		  c2[i+j]+=a[i+(n/2)]*b[j+(n/2)]
//
//  #a0-a1
//  for i in range(n/2):
//	  a[i]=a[i]-a[i+n/2]
//
//  #b1-b0
//  for i in range(n/2):
//	  b[i]=b[i+n/2]-b[i]
//
//  #c1=(a0-a1)(b1-b0)=a0b1+a1b0-a0b0-a1b1
//  for i in range(n/2):
//	  for j in range(n/2):
//		  c1[i+j]+=a[i]*b[j]
//
//  ##c1 += a0b0+a1b1
//  for i in range(n):
//	  c1[i]+=c0[i]+c2[i]
//
//  for i in range(n):
//	  c[i]+=c0[i]
//
//  for i in range(n):
//	  c[i+n/2]+=c1[i]
//
//  for i in range(n):
//	  c[i+n]+=c2[i]
//

  printf (". ");
  //x0y0
  for (int i = 0; i < n / 4; i++)
    for (int j = 0; j < n / 4; j++)
      {
	mpz_mul (tmp, xn_array[i], yn_array[j]);
	mpz_add (hn_array[i + j], hn_array[i + j], tmp);
	mpz_mod (hn_array[i + j], hn_array[i + j], p);
      }

  printf (". ");
  //x1y1
  for (int i = 0; i < n / 4; i++)
    for (int j = 0; j < n / 4; j++)
      {
	mpz_mul (tmp, xn_array[i + n / 4], yn_array[j + n / 4]);
	mpz_add (hn_array[i + j + n / 2], hn_array[i + j + n / 2], tmp);
	mpz_mod (hn_array[i + j + n / 2], hn_array[i + j + n / 2], p);
      }

  printf (". ");
  //x0=x0-x1
  //y0=y1-y0
  for (int i = 0; i < n / 4; i++)
    {
      mpz_sub (xn_array[i], xn_array[i], xn_array[i + n / 4]);
      mpz_sub (yn_array[i], yn_array[i + n / 4], yn_array[i]);
      mpz_mod (xn_array[i], xn_array[i], p);

    }

  printf (". ");
  //c1=(x0y0)+c0+c2
  for (int i = 0; i < n / 4; i++)
    for (int j = 0; j < n / 4; j++)
      {
	mpz_mul (tmp, xn_array[i], yn_array[j]);
	mpz_add (c1_array[i + j], c1_array[i + j], tmp);
	mpz_mod (c1_array[i + j], c1_array[i + j], p);
      }

  for (int i = 0; i < n / 2; i++)
    {
      mpz_add (c1_array[i], c1_array[i], hn_array[i]);
      mpz_add (c1_array[i], c1_array[i], hn_array[i + n / 2]);
    }

  printf (". ");
  for (int i = 0; i < n / 2; i++)
    {
      mpz_add (hn_array[i + n / 4], hn_array[i + n / 4], c1_array[i]);
      mpz_mod (hn_array[i + n / 4], hn_array[i + n / 4], p);
    }

  printf ("\n");

  printf (". ");
  for (int i = 0; i < n; i++)
    {
      gmp_bigint_to_radix_based_vec_u32 (&host_vector_h[(i) * coefficient_size],
					 hn_array[i], radix, coefficient_size);
    }

  memcpy (host_vector_x, host_vector_h, data_size);
  free (host_vector_h);

  mpz_clear (p);
  mpz_clear (tmp);

  for (int i = 0; i < n / 2; i++)
    mpz_clear (xn_array[i]);
  for (int i = 0; i < n / 2; i++)
    mpz_clear (yn_array[i]);
  for (int i = 0; i < n / 2; i++)
    mpz_clear (c1_array[i]);
  for (int i = 0; i < n; i++)
    mpz_clear (hn_array[i]);
  free (xn_array);
  free (yn_array);
  free (hn_array);
  free (c1_array);
}

/**************************************/

//A naive implementation of polynomial multiplications
//with cofficients in a big prime field.
void
gmp_karatsuba_poly_mult_u64 (usfixn64* host_vector_x, usfixn64* host_vector_y,
			     int n, usfixn64 radix, int coefficient_size)
{
  int data_size = n * coefficient_size * sizeof(usfixn64);
  usfixn64 * host_vector_h = (usfixn64*) malloc (data_size);
  memset (host_vector_h, 0x00, data_size);

  mpz_t tmp, p;
  mpz_init (p);
  mpz_init (tmp);

  gmp_compute_srgfn_u64 (p, radix, coefficient_size);

  int K_minus_one = 2 * coefficient_size - 1;

  mpz_t *xn_array, *yn_array, *hn_array, *c1_array, *tmp_array;

  xn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  yn_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  c1_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  tmp_array = (mpz_t*) malloc ((n / 2) * sizeof(mpz_t));
  hn_array = (mpz_t*) malloc ((n) * sizeof(mpz_t));

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
cilk_for (int i = 0; i < n / 2; i++)
    {
      mpz_init (xn_array[i]);
      u64vector_to_radix_based_bigint_gmp (xn_array[i],
					   &host_vector_x[i * coefficient_size],
					   radix, coefficient_size);
    }

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init (yn_array[i]);
      u64vector_to_radix_based_bigint_gmp (yn_array[i],
					   &host_vector_y[i * coefficient_size],
					   radix, coefficient_size);
    }

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n; i++)
    {
      mpz_init_set_ui (hn_array[i], 0);
    }

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init_set_ui (c1_array[i], 0);
    }

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    {
      mpz_init_set_ui (tmp_array[i], 0);
    }

//  #a0b0
//  for i in range(n/2):
//	  for j in range(n/2):
//		  c0[i+j]+=a[i]*b[j]
//
//  #a1b1
//  for i in range(n/2):
//	  for j in range(n/2):
//		  c2[i+j]+=a[i+(n/2)]*b[j+(n/2)]
//
//  #a0-a1
//  for i in range(n/2):
//	  a[i]=a[i]-a[i+n/2]
//
//  #b1-b0
//  for i in range(n/2):
//	  b[i]=b[i+n/2]-b[i]
//
//  #c1=(a0-a1)(b1-b0)=a0b1+a1b0-a0b0-a1b1
//  for i in range(n/2):
//	  for j in range(n/2):
//		  c1[i+j]+=a[i]*b[j]
//
//  ##c1 += a0b0+a1b1
//  for i in range(n):
//	  c1[i]+=c0[i]+c2[i]
//
//  for i in range(n):
//	  c[i]+=c0[i]
//
//  for i in range(n):
//	  c[i+n/2]+=c1[i]
//
//  for i in range(n):
//	  c[i+n]+=c2[i]
//

  printf ("*");
  //x0y0
//  for (int i = 0; i < n / 4; i++)
//    for (int j = 0; j < n / 4; j++)
//      {
////	mpz_mul (tmp, xn_array[i], yn_array[j]);
////	mpz_add (hn_array[i + j], hn_array[i + j], tmp);
//	mpz_mul (tmp_array[i], xn_array[i], yn_array[j]);
//	mpz_add (hn_array[i + j], hn_array[i + j], tmp_array[i]);
//	mpz_mod (hn_array[i + j], hn_array[i + j], p);
//      }

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int h = 0; h < (n / 2); h++)
    {
      for (int i = 0; i < n / 4; i++)

	{
	  int j = h - i;
	  if ((j < (n / 4)) && (j >= 0))
	    {
	      //	mpz_mul (tmp, xn_array[i], yn_array[j]);
	      //	mpz_add (hn_array[i + j], hn_array[i + j], tmp);
	      mpz_mul (tmp_array[h], xn_array[i], yn_array[j]);
	      mpz_add (hn_array[h], hn_array[h], tmp_array[h]);
	    }

	}
      mpz_mod (hn_array[h], hn_array[h], p);
    }
  printf ("*");
  //x1y1
//      for (int i = 0; i < n / 4; i++)
//	for (int j = 0; j < n / 4; j++)
//	  {
//	    mpz_mul (tmp, xn_array[i + n / 4], yn_array[j + n / 4]);
//	    mpz_add (hn_array[i + j + n / 2], hn_array[i + j + n / 2], tmp);
////	mpz_mul (tmp_array[i], xn_array[i + n / 4], yn_array[j + n / 4]);
////	mpz_add (hn_array[i + j + n / 2], hn_array[i + j + n / 2],
////		 tmp_array[i]);
//	    mpz_mod (hn_array[i + j + n / 2], hn_array[i + j + n / 2], p);
//	  }

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int h = 0; h < n / 2; h++)
    {
      for (int i = 0; i < n / 4; i++)

	{
	  int j = h - i;
	  if ((j >= 0) && (j < n / 4))
	    {
	      mpz_mul (tmp_array[h], xn_array[i + n / 4], yn_array[j + n / 4]);
	      mpz_add (hn_array[h + n / 2], hn_array[h + n / 2], tmp_array[h]);
	      //	mpz_mul (tmp_array[i], xn_array[i + n / 4], yn_array[j + n / 4]);
	      //	mpz_add (hn_array[i + j + n / 2], hn_array[i + j + n / 2],
	      //		 tmp_array[i]);

	    }
	}
      mpz_mod (hn_array[h + n / 2], hn_array[h + n / 2], p);
    }

  printf ("*");
  //x0=x0-x1
  //y0=y1-y0

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 4; i++)
    {
      mpz_sub (xn_array[i], xn_array[i], xn_array[i + n / 4]);
      mpz_sub (yn_array[i], yn_array[i + n / 4], yn_array[i]);
      mpz_mod (xn_array[i], xn_array[i], p);
      mpz_mod (yn_array[i], yn_array[i], p);

    }

  printf ("*");
  //c1=(x0y0)+c0+c2
//	      for (int i = 0; i < n / 4; i++)
//		for (int j = 0; j < n / 4; j++)
//		  {
//		    mpz_mul (tmp, xn_array[i], yn_array[j]);
//		    mpz_add (c1_array[i + j], c1_array[i + j], tmp);
//		    mpz_mod (c1_array[i + j], c1_array[i + j], p);
//		  }
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int h = 0; h < n / 2; h++)
    {
      for (int i = 0; i < n / 4; i++)

	{
	  int j = h - i;

	  if ((j >= 0) && (j < n / 4))
	    {
//			mpz_mul (tmp, xn_array[i], yn_array[j]);
//			mpz_add (c1_array[i + j], c1_array[i + j], tmp);
	      mpz_mul (tmp_array[h], xn_array[i], yn_array[j]);
	      mpz_add (c1_array[h], c1_array[h], tmp_array[h]);
	    }

	}
      mpz_mod (c1_array[h], c1_array[h], p);
    }
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    {
      mpz_add (c1_array[i], c1_array[i], hn_array[i]);
      mpz_add (c1_array[i], c1_array[i], hn_array[i + n / 2]);
    }

  printf ("*");
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    {
      mpz_add (hn_array[i + n / 4], hn_array[i + n / 4], c1_array[i]);
      mpz_mod (hn_array[i + n / 4], hn_array[i + n / 4], p);
    }

  printf ("*");
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n; i++)
    {
      bigint_to_u64vector_radixbased_gmp (
	  &host_vector_h[(i) * coefficient_size], hn_array[i], radix,
	  coefficient_size);

    }

  memcpy (host_vector_x, host_vector_h, data_size);
  free (host_vector_h);

  mpz_clear (p);
  mpz_clear (tmp);

#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    mpz_clear (xn_array[i]);
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    mpz_clear (yn_array[i]);
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    mpz_clear (c1_array[i]);
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n / 2; i++)
    mpz_clear (tmp_array[i]);
#pragma cilk_grainsize = (KARATSUBA_GRAIN_FACTOR*(cilk::current_worker_count))
  //cilk
  for (int i = 0; i < n; i++)
    mpz_clear (hn_array[i]);
  free (xn_array);
  free (yn_array);
  free (hn_array);
  free (c1_array);
  free (tmp_array);
}

/**************************************/

void
gmp_inv_mod_p (mpz_t inv, mpz_t x, mpz_t p)
{
  mpz_t p2;
  mpz_init (p2);

  mpz_sub_ui (p2, p, 2);
  mpz_powm (inv, x, p2, p);

  mpz_clear (p2);
}

/**************************************/

void
gmp_inv_mod_p_u64 (usfixn64 * inv, const usfixn64 x, const usfixn64 p)
{

  mpz_t p_zz, x_zz, inv_zz, p2_zz;
  mpz_inits (p_zz, x_zz, inv_zz, p2_zz, NULL);
  mpz_set_u64 (p_zz, p);
  mpz_set_u64 (x_zz, x);
  mpz_mod (x_zz, x_zz, p_zz);
  mpz_sub_ui (p2_zz, p_zz, 2);
  mpz_powm (inv_zz, x_zz, p2_zz, p_zz);
  *inv = mpz_get_u64 (inv_zz);

//  printf("inv=%llu, n=%llu, p=%llu\n",*inv, x, p);

  mpz_clears (p_zz, x_zz, inv_zz, p2_zz, NULL);

}

/**************************************/

int
gmp_inv_radix_based_bigint_modp_u64 (usfixn64* result_vector,
				     usfixn64* x_vector, usfixn64 radix, int k,
				     mpz_t p)
{

  mpz_t result_zz, x_zz;
  mpz_inits (result_zz, x_zz, NULL);

  u64vector_to_radix_based_bigint_gmp (x_zz, x_vector, radix, k);
  gmp_inv_mod_p (result_zz, x_zz, p);
  bigint_to_u64vector_radixbased_gmp (result_vector, result_zz, radix, k);

  mpz_clears (result_zz, x_zz, NULL);
  return EXIT_SUCCESS;
}

/**************************************/

int
gmp_mult_radix_based_bigint_by_scalar_modp_u64 (usfixn64* result_vector,
						const mpz_t scalar_zz,
						usfixn64 radix, int k, mpz_t p)
{
  mpz_t result_zz;
  mpz_init (result_zz);
  u64vector_to_radix_based_bigint_gmp (result_zz, &result_vector[0], radix, k);

  mpz_mul (result_zz, result_zz, scalar_zz);
  mpz_mod (result_zz, result_zz, p);

  bigint_to_u64vector_radixbased_gmp (result_vector, result_zz, radix, k);

  mpz_clear(result_zz);
  return EXIT_SUCCESS;
}

/**************************************/

void
compute_srgfn_p_gmp (mpz_t p, usfixn64 radix, int coefficient_size)
{
  mpz_t m, r;

  mpz_init_set_ui (m, 0);
  mpz_init_set_ui (r, radix);

  mpz_pow_ui (m, r, coefficient_size);
  mpz_add_ui (p, m, 1);
  mpz_clears (m, r, NULL);
}

/**************************************/
int
compute_nth_root_of_unity_for_small_prime (usfixn64 prime, usfixn64* omega,
					   int n)
{
  if ((prime - 1) % n != 0)
    {
      printf (
	  "Error: %d does not divide %llu - 1!\nCannot compute %d-th primitive root of unity of %llu.\n",
	  n, prime, n, prime);
      return 0;
    }
  mpz_t p_zz;
  mpz_init (p_zz);
  mpz_set_u64 (p_zz, prime);

  mpz_t p1_zz;
  mpz_init (p1_zz);
  mpz_set_u64 (p1_zz, prime - 1);

  mpz_t p2_zz;
  mpz_init (p2_zz);
  mpz_set_u64 (p2_zz, (prime - 1) / 2);

  mpz_t q_zz;
  mpz_init (q_zz);
  mpz_set_u64 (q_zz, (prime - 1) / n);

  usfixn64 c;
  mpz_t c_zz;
  mpz_init (c_zz);

  int i = 0;

  srand (time (NULL));
  while (i < 20)
    {
      c = rand ();
      mpz_set_u64 (c_zz, c);
      mpz_powm (c_zz, c_zz, p2_zz, p_zz);
      if (mpz_cmp (c_zz, p1_zz) == 0)
	{
	  mpz_set_u64 (c_zz, c);
	  mpz_powm (c_zz, c_zz, q_zz, p_zz);
	  *omega = mpz_get_u64 (c_zz);
	  mpz_clears (p_zz, p1_zz, p2_zz, q_zz, c_zz, NULL);
	  return 0;
	}
      i++;
    }
  printf ("No primitive root found!\nPlease try again.\n");
  mpz_clears (p_zz, p1_zz, p2_zz, q_zz, c_zz, NULL);
  *omega = 0;
  return 0;
}

/**************************************/

int
compute_nth_root_of_unity_for_small_prime_py (usfixn64 prime, usfixn64* omega,
					      int n)
{
  int k = 1;
//      int K = (k << 1);
  usfixn64 R = prime - 1;

  char *bpas_src = (char*) malloc (4096);
  char *omega_path = (char*) malloc (4096);
  sprintf (omega_path, "%s/current_omega.tmp", getenv ("PWD"));
  char * cmd_str = (char*) malloc (MAX_CMD_LEN);
  int length = 0;

  sprintf (bpas_src, "%s/src/FiniteFields/", getenv ("BPAS_HOME"));
  length += sprintf (cmd_str + length, "python ");
  length += sprintf (cmd_str + length, "%s", bpas_src);
  length += sprintf (cmd_str + length, "py_roots_of_unity.py ");
  length += sprintf (cmd_str + length, "%d %llu %d ", n, R, k);
  length += sprintf (cmd_str + length, "%s", omega_path);
  length += sprintf (cmd_str + length, " >/dev/null");

#if VERBOSE>=1
  printf ("[compute_nth_root_of_unity_for_srgfn]...\n");
  printf ("[%s] ...\n", cmd_str);
  printf("%s",shortline);
#endif
  if (system (cmd_str) != EXIT_SUCCESS)
    {
      printf ("ERROR: [%s]", cmd_str);
      exit (EXIT_FAILURE);
    }

  readVectorFromFile (omega, 1, k, "current_omega.tmp");

  int verification_enabled = 1; // VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;
  ////////////////////////////////////////////////////
  //verifying if computed omega is n-th root of unity
  ////////////////////////////////////////////////////
  if (verification_enabled == 1)
    {
      mpz_t omega_zz;
      mpz_init (omega_zz);
      u64vector_to_radix_based_bigint_gmp (omega_zz, omega, R, k);

      mpz_t p;
      mpz_init (p);
      compute_srgfn_p_gmp (p, R, k);

      mpz_pow_ui (omega_zz, omega_zz, n);
      mpz_mod (omega_zz, omega_zz, p);

#if VERBOSE
      printf ("[verifying (omega^N) mod p == 1] ...");
#endif
      if (mpz_get_ui (omega_zz) == 1)
	{
	  verification_result = EXIT_SUCCESS;
#if VERBOSE
	  printf ("PASSED\n");
#endif
	}
      else
	{
	  verification_result = EXIT_SUCCESS;
#if VERBOSE
	  printf ("FAILED\n");
#endif
	}
#if VERBOSE
      printf ("%s", shortline);
#endif
      mpz_clear (omega_zz);
      mpz_clear (p);
    }

  sprintf (cmd_str, "%s", omega_path);
  if (remove (cmd_str) != 0)
    {
      printf ("ERROR: [small][%s][ERROR=%s]\n", cmd_str, strerror (errno));
      exit (EXIT_FAILURE);
    }

  free (cmd_str);
  free (bpas_src);
  free (omega_path);
  return verification_result;
}

/**************************************/
int
compute_nth_root_of_unity_for_srgfn (srgfn_prime prime, usfixn64* omega, int n)
{

  mpz_t p_zz;
  mpz_init (p_zz);
  compute_srgfn_p_gmp(p_zz, prime.radix, prime.k);

//  ToDO: fix the following part.
//  if ((prime - 1) % n != 0)
//    {
//      printf (
//	  "Error: %d does not divide %llu - 1!\nCannot compute %d-th primitive root of unity of %llu.\n",
//	  n, prime, n, prime);
//      return 0;
//    }

  mpz_t p1_zz;
  mpz_init (p1_zz);
  mpz_sub_ui(p1_zz, p_zz, 1);

  mpz_t p2_zz;
  mpz_init (p2_zz);
  mpz_div_ui(p2_zz, (p1_zz) , 2);

  mpz_t q_zz;
  mpz_init (q_zz);
  mpz_div_ui(q_zz, (p1_zz),  n);

  usfixn64 c;
  mpz_t c_zz;
  mpz_init (c_zz);

  int etta_exp=n/(2*prime.k);


  mpz_t t_zz;
  mpz_init(t_zz);


  int result=EXIT_FAILURE;
  srand (time (NULL));

//  gmp_randstate_t rand_state;
//  gmp_randinit_mt(rand_state);

  for(int i=0;i< 1000;i++)
    {
      c = rand ();
      mpz_set_u64 (c_zz, c);
//      mpz_urandomm(c_zz, rand_state, p_zz);
      mpz_powm (c_zz, c_zz, p2_zz, p_zz);
      if (mpz_cmp (c_zz, p1_zz) == 0)
	{
	  mpz_set_u64 (c_zz, c);
	  mpz_powm (c_zz, c_zz, q_zz, p_zz);
//	  *omega = mpz_get_u64 (c_zz);
	  mpz_powm_ui(t_zz, c_zz, etta_exp, p_zz);

	  if (mpz_cmp_ui(t_zz, prime.radix)==0)
	    {
	      bigint_to_u64vector_radixbased_gmp(omega, c_zz, prime.radix, prime.k);
	      result=EXIT_SUCCESS;
	      break;
	    }
	}
    }
  mpz_clears (p_zz, p1_zz, p2_zz, q_zz, c_zz, t_zz, NULL);
//  gmp_randclear(rand_state);

  if (result!=EXIT_SUCCESS)
      {
        printf ("[ERROR: No primitive root found!]");
        memset(omega, 0x00, prime.k*sizeof(usfixn64));
        exit(EXIT_FAILURE);
      }
  return result;
}

/**************************************/

int
compute_nth_root_of_unity_for_srgfn_py (srgfn_prime prime, usfixn64* omega, int n)
{
  int k = prime.k;
//      int K = (k << 1);
  usfixn64 R = prime.radix;

  char *omega_path = (char*) malloc (2048);
  sprintf (omega_path, "%s/current_omega.tmp", getenv ("PWD"));

  char * cmd_str = (char*) malloc (MAX_CMD_LEN);
  int length = 0;
  char *bpas_src = (char*) malloc (2048);
  length += sprintf (cmd_str + length, "python ");
  sprintf (bpas_src, "%s/src/FiniteFields/", getenv ("BPAS_HOME"));
  length += sprintf (cmd_str + length, "%s", bpas_src);
  length += sprintf (cmd_str + length, "py_roots_of_unity.py ");
  length += sprintf (cmd_str + length, "%d %llu %d ", n, R, k);
  length += sprintf (cmd_str + length, "%s", omega_path);
  length += sprintf (cmd_str + length, " >/dev/null");

#if VERBOSE>=1
  printf ("[compute_nth_root_of_unity_for_srgfn]...\n");
  printf ("[%s] ...\n", cmd_str);
  printf("%s",shortline);
#endif
  if (system (cmd_str) == -1)
    {

      printf ("ERROR: [%s]", cmd_str);
      exit (EXIT_FAILURE);
    }

  readVectorFromFile (omega, 1, k, omega_path);
  sprintf (cmd_str, "%s", omega_path);
  if (remove (cmd_str) != 0)
    {
      printf ("ERROR: [large][%s][ERROR=%s]\n", cmd_str, strerror (errno));
      exit (EXIT_FAILURE);
    }
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;
  ////////////////////////////////////////////////////
  //verifying if computed omega is n-th root of unity
  ////////////////////////////////////////////////////
  if (verification_enabled == 1)
    {
      mpz_t omega_zz;
      mpz_init (omega_zz);
      u64vector_to_radix_based_bigint_gmp (omega_zz, omega, R, k);

      mpz_t p;
      mpz_init (p);
      compute_srgfn_p_gmp (p, R, k);

      mpz_pow_ui (omega_zz, omega_zz, n);
      mpz_mod (omega_zz, omega_zz, p);

#if VERBOSE
      printf ("[verifying (omega^N) mod p == 1] ...");
#endif
      if (mpz_get_ui (omega_zz) == 1)
	{
	  verification_result = EXIT_SUCCESS;
#if VERBOSE
	  printf ("PASSED\n");
#endif
	}
      else
	{
	  verification_result = EXIT_SUCCESS;
#if VERBOSE
	  printf ("FAILED\n");
#endif
	}
#if VERBOSE
      printf ("%s", shortline);
#endif
      mpz_clear (omega_zz);
      mpz_clear (p);
    }

  free (cmd_str);
  free (bpas_src);
  free (omega_path);
  return verification_result;
}

/**************************************/
#endif
