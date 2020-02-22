#include <stdio.h>
#include <stdlib.h>

#include "bpas.h"

/**************************************/

typedef struct dft_data
{
  usfixn64 *x_data, *y_data, *x_data_copy, *y_data_copy, *omega, *omega_inv;
  srgfn_prime prime;
  int K;
  int e;
  mpz_t p_zz, n_inv_zz, omega_zz, omega_inv_zz;
  int fft_based_mult_enabled;
  int verification_enabled;
  int poly_mult_enabled;
  int data_size;
} dft_data_t;

/**************************************/

int
init_dft_six_step_data (dft_data_t* data, int K, int e,
			int fft_based_mult_enabled, int verification_enabled,
			int poly_mult_enabled)
{
//  if ((K > 128) || (e > 5) || (fft_based_mult_enabled > 1)
  if ((e > 5) || (fft_based_mult_enabled > 1) || (verification_enabled > 1))
    {
      printf ("[ERROR: incorrect or NULL parameter passed to %s]\n", __func__);
      exit (EXIT_FAILURE);
    }
  data->K = K;
  data->e = e;
  data->fft_based_mult_enabled = fft_based_mult_enabled;
  data->verification_enabled = verification_enabled;
  data->poly_mult_enabled = poly_mult_enabled;

  /////////////////////////////////////
  //// select SRGFN
  /////////////////////////////////////
  srgfn_prime * srgfn_ptr;
  switch (data->K / 2)
    {
    case 2:
      srgfn_ptr = p_list_1;
      break;

    case 4:
      srgfn_ptr = &p_list_2[0];
      break;

    case 8:
      srgfn_ptr = &p_list_3[0];
      break;

    case 16:
      srgfn_ptr = &p_list_4[0];
      break;

    case 32:
      srgfn_ptr = &p_list_5[0];
      break;

    case 64:
      srgfn_ptr = &p_list_6[0];
      break;

    case 128:
      srgfn_ptr = &p_list_7[0];
      break;

    case 256:
      srgfn_ptr = &p_list_8[0];
      break;

    default:
      {
	printf ("[ERROR: proper k is not specified!]\n");
	exit (EXIT_FAILURE);
      }
    }

  data->prime.k = srgfn_ptr->k;
  data->prime.radix = srgfn_ptr->radix;

  printf ("[p=r^k+1={%llu}^{%d}+1]\n", data->prime.radix, data->prime.k);

  /////////////////////////////////////
  //// init equivalent gmp prime
  /////////////////////////////////////
  mpz_init (data->p_zz);
  gmp_compute_srgfn_u64 (data->p_zz, data->prime.radix, data->prime.k);
  /////////////////////////////////////
  /////////////////////////////////////
  //// init memory for test data
  /////////////////////////////////////
  usfixn64 n = pow (data->K, data->e);
  usfixn64 data_size = (n * data->prime.k * sizeof(usfixn64));

  data->x_data = (usfixn64*) malloc (data_size);
  memset (data->x_data, 0x00, (data_size));
  if (poly_mult_enabled)
    {
      data->y_data = (usfixn64*) malloc (data_size);
      memset (data->y_data, 0x00, (data_size));
    }
  data->data_size = data_size;

  ///////////////////////////////////////////////////////
  ///// LOADING INPUT VECTORS x and y
  ///////////////////////////////////////////////////////
  ///// WARNING: in case you are using big FFT-based polynomial
  ///// multiplication, make sure each vector of big elements is
  ///// half full (i.e., upper n/2 elements are 0), otherwise,
  ///// multiplication will FAIL.x
  ///////////////////////////////////////////////////////

  srand (time (NULL));
  for (int i = 0; i < n / 2; i++)
    for (int j = 0; j < data->prime.k; j++)
      {
	data->x_data[i * data->prime.k + j] = data->prime.radix - 1 - i; //rand () % prime->radix;  //- 1 - (i);
	if (poly_mult_enabled)
	  data->y_data[i * data->prime.k + j] = data->prime.radix - 2 - i; //rand () % prime.radix;
      }

  if (data->verification_enabled == 1)
    {
#if VERBOSE >=1
      printf (
	  "[copying x_data and y_data which will be used for verification] ...\n");
      printf("%s",shortline);
#endif
      data->x_data_copy = (usfixn64*) malloc (data_size);
      memcpy (data->x_data_copy, data->x_data, data_size);
      if (poly_mult_enabled)
	{
	  data->y_data_copy = (usfixn64*) malloc (data_size);
	  memcpy (data->y_data_copy, data->y_data, data_size);
	}
    }

  /////////////////////////////////////
  //// precompute root of unity and its inverse
  //// (omega and omega_inv, respectively).
  /////////////////////////////////////

  data->omega = (usfixn64*) malloc (data->prime.k * sizeof(usfixn64));
  data->omega_inv = (usfixn64*) malloc (data->prime.k * sizeof(usfixn64));

  compute_nth_root_of_unity_for_srgfn (data->prime, data->omega, n);
  gmp_inv_radix_based_bigint_modp_u64 (data->omega_inv, data->omega,
				       data->prime.radix, data->prime.k,
				       data->p_zz);

  mpz_init (data->omega_zz);
  mpz_init (data->omega_inv_zz);
  u64vector_to_radix_based_bigint_gmp (data->omega_zz, data->omega,
				       data->prime.radix, data->prime.k);
  u64vector_to_radix_based_bigint_gmp (data->omega_inv_zz, data->omega_inv,
				       data->prime.radix, data->prime.k);
  /////////////////////////////////////
  //// this is critical on some platforms,
  //// can lead to terrible gmp error.
  /////////////////////////////////////
  if (sizeof(unsigned long) != sizeof(unsigned long long))
    {
      printf (
	  "   [ERROR:\n"
	  "   On this machine, unsigned long is not as wide as 'unsigned long long'\n"
	  "   This can fail much of the gmp convert in and convert out functions!\n"
	  "   EXIT FAILURE!\n");
      exit (EXIT_FAILURE);
    }

  /////////////////////////////////////
  //// compute inverse of n (size of input vector)
  /////////////////////////////////////

#if VERBOSE>=1
  printf ("[computing n_inv = n^{-1} for inverse DFT]...");
#endif
  mpz_t n_zz;
  mpz_init (n_zz);
  mpz_init (data->n_inv_zz);

  mpz_set_u64 (n_zz, n);
  gmp_inv_mod_p (data->n_inv_zz, n_zz, data->p_zz);
  mpz_clear (n_zz);
#if VERBOSE>=1
  printf ("DONE\n");
  printf("%s",shortline);
#endif

  ///////////////////////////////////////
  //// init data for fft-based multiplication
  //// of two numbers.
  /////////////////////////////////////
  if (data->fft_based_mult_enabled == 1)
    {
      usfixn64 conv_p1, conv_p2;
      get_proper_conv_p1p2 (&conv_p1, &conv_p2, data->prime);
      init_fft_based_bigint_mult (conv_p1, conv_p2, data->prime.k);
      init_gfpf_mult_data (&t_crt_data_global, &t_lhc_data_global,
			   data->prime.radix, conv_p1, conv_p2);
    }
  /////////////////////////////////////

  return EXIT_SUCCESS;
}

/**************************************/

int
clear_dft_six_step_data (dft_data_t* data)
{
  /////////////////////////////////////
  mpz_clear (data->p_zz);
  /////////////////////////////////////
  free (data->x_data);

  /////////////////////////////////////
  if (data->poly_mult_enabled == 1)
    free (data->y_data);
  if (data->verification_enabled == 1)
    {
      free (data->x_data_copy);
      if (data->poly_mult_enabled == 1)
	free (data->y_data_copy);
    }

  /////////////////////////////////////
  free (data->omega);
  free (data->omega_inv);

  mpz_clear (data->omega_zz);
  mpz_clear (data->omega_inv_zz);

  /////////////////////////////////////
  mpz_clear (data->n_inv_zz);
  /////////////////////////////////////
  if (data->fft_based_mult_enabled == 1)
    {
      clear_fft_based_bigint_mult_data ();
      clear_gfpf_mult_data (&t_crt_data_global, &t_lhc_data_global);
    }
  /////////////////////////////////////
  return EXIT_SUCCESS;
}

/**************************************/

int
verify_dft_based_poly_mult (usfixn64* x_data, usfixn64 *x_data_copy,
			    usfixn64 *y_data_copy, int n, srgfn_prime prime)
{

  int verification_result = EXIT_SUCCESS;
  printf ("%s", shortline);
  printf ("[verifying fft_based_poly_mult] ");
  /////////////////////////////////
  //// backup
  //// either of the following functions can
  //// be replaced with a poly multiplication
  //// function, also, keeping coefficients mod p.
  /////////////////////////////////
  //gmp_plain_poly_mult_u64 (xdata_copy, ydata_copy, vector_size, R,
  //			 current_prime.k);
  /////////////////////////////////
  gmp_karatsuba_poly_mult_u64 (x_data_copy, y_data_copy, n, prime.radix,
			       prime.k);

  int test_failed = 0;
  for (int i = 0; i < n; i++)
    {
//      printVector (&x_data[i * prime.k], 1, prime.k, "x");
//      printVector (&x_data_copy[i * prime.k], 1, prime.k, "xcopy");
      for (int j = 0; j < prime.k; j++)
	if (x_data_copy[i * prime.k + j] != x_data[i * prime.k + j])
	  {
	    printf ("ERROR: MISMATCH at [%d,%d] ...\n", i, j);
	    test_failed = 1;
	    break;
	  }
      if (test_failed == 1)
	break;
    }

  if (test_failed == 0)
    {
      verification_result = EXIT_FAILURE;
      printf (" VERIFIED\n");
    }
  else
    {
      verification_result = EXIT_FAILURE;
      printf (" FAILED\n");
    }
  printf ("%s", shortline);
//  print_verification_msg()
  return verification_result;
}

/**************************************/

int
verify_equal_gfpf_vectors (usfixn64* x_data, usfixn64 *x_data_copy, int n,
			   srgfn_prime prime)
{

  int verification_result = EXIT_SUCCESS;
  printf ("%s", shortline);
  printf ("[%s]\n", __func__);

  int test_failed = 0;
  for (int i = 0; i < n; i++)
    {
//      printVector (&x_data[i * prime.k], 1, prime.k, "x");
//      printVector (&x_data_copy[i * prime.k], 1, prime.k, "xcopy");
      for (int j = 0; j < prime.k; j++)
	if (x_data_copy[i * prime.k + j] != x_data[i * prime.k + j])
	  {
	    printf ("ERROR: MISMATCH at [%d,%d] ...\n", i, j);
	    test_failed = 1;
	    break;
	  }
      if (test_failed == 1)
	break;
    }

  if (test_failed == 0)
    {
      verification_result = EXIT_SUCCESS;
    }
  else
    {
      verification_result = EXIT_FAILURE;
    }
  print_verification_msg (__func__, verification_result);
  printf ("%s", shortline);
  return verification_result;
}

/**************************************/

int
test_dft_six_step_gfpf (int K, int e, int fft_based_mult_enabled)
{

  if (K > 64)
    {
      printf ("[WARNING: K>64 is not fully tested!]\n");
//      return EXIT_FAILURE;
    }
  int n = pow (K, e);
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;

  dft_data_t data;
  int poly_mult_enabled = 0;
  init_dft_six_step_data (&data, K, e, fft_based_mult_enabled,
			  verification_enabled, poly_mult_enabled);

  /////////////////////////////////////
  printf ("%s", longline);
  printf (
      "[DFT_N=K^e, n=%d, K=%d, e=%d, p={%llu}^{%d}+1, fft-based-mult=%d] ...\n",
      n, data.K, data.e, data.prime.radix, data.prime.k,
      data.fft_based_mult_enabled);
  printf ("%s", shortline);
  /////////////////////////////////////

  int compute_inverse = 0;
  int verbose_on = 1, verbose_off = 0;
  cpu_timer t_fft_based_mult;

  /////////////////////////////////////
  /////////////////////////////////////

  usfixn64 * x_data_copy = (usfixn64*) malloc (data.data_size);
  memcpy (x_data_copy, data.x_data, data.data_size);

  //x=DFT_n(x)
  DFT_general_big_elements (data.x_data, data.K, data.e, data.prime.radix,
			    data.p_zz, data.omega, data.fft_based_mult_enabled,
			    compute_inverse, verbose_on);

#if VERBOSE>=1
  printf ("DONE\n");
  printf("%s",shortline);
#endif

  ///////////////////////////////////////
  if (verification_enabled)
    {
      //x=DFT^{-1}_n(x)
      compute_inverse = 1;
      DFT_general_big_elements (data.x_data, data.K, data.e, data.prime.radix,
				data.p_zz, data.omega_inv,
				data.fft_based_mult_enabled, compute_inverse,
				verbose_off);

      /**************************************/
#if VERBOSE>=1
      printf ("[multiplying convolution result by n^{-1}]...");
#endif

      //x={1/n}*x
      for (int i = 0; i < n; i++)
	{
	  gmp_mult_radix_based_bigint_by_scalar_modp_u64 (
	      &(data.x_data[i * data.prime.k]), data.n_inv_zz, data.prime.radix,
	      data.prime.k, data.p_zz);
	}
      verification_result = verify_equal_gfpf_vectors (data.x_data,
						       data.x_data_copy, n,
						       data.prime);
    }
  ///////////////////////////////////////
  clear_dft_six_step_data (&data);
  return verification_result;
}

/**************************************/

int
test_dft_based_poly_mult_gfpf (int K, int e, int fft_based_mult_enabled)
{

  int n = pow (K, e);
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;

  dft_data_t data;
  int poly_mult_enabled = 1;
  init_dft_six_step_data (&data, K, e, fft_based_mult_enabled,
			  verification_enabled, poly_mult_enabled);

  /////////////////////////////////////
  printf ("%s", longline);
  printf (
      "[DFT_N=K^e, n=%d, K=%d, e=%d, p={%llu}^{%d}+1, fft-based-mult=%d] ...\n",
      n, data.K, data.e, data.prime.radix, data.prime.k,
      data.fft_based_mult_enabled);
  printf ("%s", longline);
  /////////////////////////////////////

  cpu_timer t_fft_based_mult;
  timer_record_start (&t_fft_based_mult);

  //x=DFT_n(x)
  int compute_inverse = 0;
  int verbose_on = 1, verbose_off = 0;
  DFT_general_big_elements (data.x_data, data.K, data.e, data.prime.radix,
			    data.p_zz, data.omega, data.fft_based_mult_enabled,
			    compute_inverse, verbose_on);

  //y=DFT_n(y)
  DFT_general_big_elements (data.y_data, data.K, data.e, data.prime.radix,
			    data.p_zz, data.omega, data.fft_based_mult_enabled,
			    compute_inverse, verbose_on);
  //x=DFT_n(x)*DFT_n(y)
  convolution_mult (data.x_data, data.y_data, n, &(data.prime), data.p_zz,
		    data.fft_based_mult_enabled);

  //x=DFT^{-1}_n(x)
  compute_inverse = 1;
  DFT_general_big_elements (data.x_data, data.K, data.e, data.prime.radix,
			    data.p_zz, data.omega_inv,
			    data.fft_based_mult_enabled, compute_inverse,
			    verbose_on);

  /**************************************/
#if VERBOSE>=1
  printf ("[multiplying convolution result by n^{-1}]...");
#endif

  //x={1/n}*x
  for (int i = 0; i < n; i++)
    {
      gmp_mult_radix_based_bigint_by_scalar_modp_u64 (
	  &(data.x_data[i * data.prime.k]), data.n_inv_zz, data.prime.radix,
	  data.prime.k, data.p_zz);
    }
#if VERBOSE>=1
  printf ("DONE\n");
  printf("%s",shortline);
#endif

  ///////////////////////////////////////

  if (verification_enabled)
    {
      verification_result = verify_dft_based_poly_mult (data.x_data,
							data.x_data_copy,
							data.y_data_copy, n,
							data.prime);
    }

  ///////////////////////////////////////
  timer_record_stop (&t_fft_based_mult);
  timer_get_elapsed_time (&t_fft_based_mult, "t_fft_based_mult_gfpf", 1);
  printf ("%s", longline);
  ///////////////////////////////////////

  clear_dft_six_step_data (&data);
  return verification_result;
}

/**************************************/

int
test_dft_six_step_aux_memory_gmp (int K, int e)
{
  int n = pow (K, e);
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;

  dft_data_t data;
  int poly_mult_enabled = 0;

  init_dft_six_step_data (&data, K, e, 0, verification_enabled,
			  poly_mult_enabled);

  /////////////////////////////////////
  printf ("%s", longline);
  printf ("[computing DFT_N=K^e, n=%d, K=%d, e=%d, p={%llu}^{%d}+1] ...\n", n,
	  data.K, data.e, data.prime.radix, data.prime.k);
  printf ("%s", shortline);

  /////////////////////////////////////
  mpz_t * x_zz;
  x_zz = (mpz_t*) malloc (n * sizeof(mpz_t));

  for (int i = 0; i < n; i++)
    {
      mpz_init (x_zz[i]);
      u64vector_to_radix_based_bigint_gmp (x_zz[i],
					   &data.x_data[i * data.prime.k],
					   data.prime.radix, data.prime.k);
    }

  /////////////////////////////////////

  int verbose_on = 1, verbose_off = 0;
  cpu_timer t_fft_based_mult;

  timer_record_start (&t_fft_based_mult);
  //x=DFT_n(x)
  int compute_inverse = 0;

  DFT_general_aux_memory_gmp (x_zz, K, e, data.omega_zz, data.p_zz, verbose_on);

  ///////////////////////////////////////
  timer_record_stop (&t_fft_based_mult);
  ///////////////////////////////////////
  ///////////////////////////////////////
  if (verification_enabled)
    {
      //x=DFT^{-1}_n(x)
      DFT_general_aux_memory_gmp (x_zz, K, e, data.omega_inv_zz, data.p_zz,
				  verbose_off);

      /**************************************/
#if VERBOSE>=1
      printf ("[multiplying convolution result by n^{-1}]...");
#endif

      //x={1/n}*x
      for (int i = 0; i < n; i++)
	{
	  mult_gmp (x_zz[i], data.n_inv_zz, data.p_zz);
	}

#if VERBOSE>=1
      printf ("DONE\n");
      printf("%s",shortline);
#endif

      for (int i = 0; i < n; i++)
	bigint_to_u64vector_radixbased_gmp (&(data.x_data[i * data.prime.k]),
					    x_zz[i], data.prime.radix,
					    data.prime.k);

      verification_result = verify_equal_gfpf_vectors (data.x_data,
						       data.x_data_copy, n,
						       data.prime);
    }

  ///////////////////////////////////////

  for (int i = 0; i < n; i++)
    {
      mpz_clear (x_zz[i]);
    }

  free (x_zz);

  clear_dft_six_step_data (&data);
  return verification_result;
}

/**************************************/

int
test_dft_based_poly_mult_gmp (int K, int e)
{
  int n = pow (K, e);
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;

  dft_data_t data;
  int poly_mult_enabled = 1;
  init_dft_six_step_data (&data, K, e, 0, verification_enabled,
			  poly_mult_enabled);

  /////////////////////////////////////
  /////////////////////////////////////
  mpz_t * x_zz, *y_zz;
  x_zz = (mpz_t*) malloc (n * sizeof(mpz_t));
  y_zz = (mpz_t*) malloc (n * sizeof(mpz_t));
  for (int i = 0; i < n; i++)
    {
      mpz_init (x_zz[i]);
      u64vector_to_radix_based_bigint_gmp (x_zz[i],
					   &data.x_data[i * data.prime.k],
					   data.prime.radix, data.prime.k);
      mpz_init (y_zz[i]);
      u64vector_to_radix_based_bigint_gmp (y_zz[i],
					   &data.y_data[i * data.prime.k],
					   data.prime.radix, data.prime.k);
    }

  /////////////////////////////////////

  int verbose_on = 1, verbose_off = 0;
  cpu_timer t_fft_based_mult;
  timer_record_start (&t_fft_based_mult);
//	for(int i=0;i<n_iterations;i++)
  //x=DFT_n(x)
  int compute_inverse = 0;
  DFT_general_gmp (x_zz, K, e, data.omega_zz, data.p_zz, verbose_on);

  //y=DFT_n(y)
  DFT_general_gmp (y_zz, K, e, data.omega_zz, data.p_zz, verbose_on);

  //x=DFT_n(x)*DFT_n(y)
  convolution_mult_gmp (x_zz, y_zz, n, data.p_zz);

  //x=DFT^{-1}_n(x)
  DFT_general_gmp (x_zz, K, e, data.omega_inv_zz, data.p_zz, verbose_on);

  /**************************************/
#if VERBOSE>=1
  printf ("[multiplying convolution result by n^{-1}]...");
#endif

  //x={1/n}*x
  for (int i = 0; i < n; i++)
    {
      mult_gmp (x_zz[i], data.n_inv_zz, data.p_zz);
    }

#if VERBOSE>=1
  printf ("DONE\n");
  printf("%s",shortline);
#endif

  ///////////////////////////////////////
  if (verification_enabled)
    {
      for (int i = 0; i < n; i++)
	bigint_to_u64vector_radixbased_gmp (&(data.x_data[i * data.prime.k]),
					    x_zz[i], data.prime.radix,
					    data.prime.k);

      verification_result = verification_result = verify_dft_based_poly_mult (
	  data.x_data, data.x_data_copy, data.y_data_copy, n, data.prime);
    }

  ///////////////////////////////////////
  timer_record_stop (&t_fft_based_mult);
  timer_get_elapsed_time (&t_fft_based_mult, "t_fft_based_mult_gmp", 1);
  printf ("%s", longline);
  ///////////////////////////////////////

  for (int i = 0; i < n; i++)
    {
      mpz_clear (x_zz[i]);
      mpz_clear (y_zz[i]);
    }

  free (x_zz);
  free (y_zz);

  clear_dft_six_step_data (&data);
  return verification_result;
}

/**************************************/

int
test_arbitrary_mult_gfpf (int K, int e)
{

  int n = pow (K, e);
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;

  dft_data_t data;
  int poly_mult_enabled = 1;
  init_dft_six_step_data (&data, K, e, 1, verification_enabled,
			  poly_mult_enabled);

  /////////////////////////////////////
//  printf ("%s", shortline);
//  printf (
//      "[computing mult=K^e, n=%d, K=%d, e=%d, p={%llu}^{%d}+1, fft-based-mult=%d] ...\n",
//      n, data.K, data.e, data.prime.radix, data.prime.k,
//      data.fft_based_mult_enabled);
//  printf ("%s", shortline);
  /////////////////////////////////////

  test_fft_based_arbitrary_mult (data.x_data, data.y_data, n, data.prime);

  /////////////////////////////////////
  clear_dft_six_step_data (&data);
  return verification_result;
}

/**************************************/

int
test_mult_pow_r_gfpf (int K, int e)
{

  int n = pow (K, e);
  int verification_enabled = VERIFICATION_ENABLED;
  int verification_result = EXIT_SUCCESS;

  dft_data_t data;
  int poly_mult_enabled = 0;
  init_dft_six_step_data (&data, K, e, 0, verification_enabled,
			  poly_mult_enabled);

  /**************************************/
  test_mult_pow_R_big_elements (data.x_data, data.K, data.e, data.prime.radix);
  /**************************************/
  clear_dft_six_step_data (&data);
  return verification_result;
}

/**************************************/

int
main (int argc, char ** argv)
{
  setbuf (stdout, NULL);
  int K = 32;
  int e = 2;
  int fft_based_mult_enabled = 1;
  int use_gmp = 0;

  if ((argc == 2))
    if ((strcmp (argv[1], "-h") == 0) || (strcmp (argv[1], "--help") == 0))
      {
	printf ("%s", shortline);
	printf (
	    "computing DFT_{K^e}:"
	    "args: [1]:K (default=32)  [2]:e (default=2)  [3]:use_gmp (default=0)\n");
	printf ("%s", shortline);
	exit (EXIT_SUCCESS);
      }
  if (argc > 1)
    {
      K = atoi (argv[1]);
    }

  if (argc > 2)
    {
      e = atoi (argv[2]);
    }

  if (argc > 3)
    {
      use_gmp = atoi (argv[3]);
    }

  if (use_gmp == 1)
    {
      test_dft_six_step_aux_memory_gmp (K, e);
    }
  else
    {
      fft_based_mult_enabled = 1;
      test_dft_six_step_gfpf (K, e, fft_based_mult_enabled);
    }

  //  test_dft_based_poly_mult_gfpf (K, e, fft_based_mult_enabled);
  //  test_dft_based_poly_mult_gmp (K, e);

  return EXIT_SUCCESS;
}

/**************************************/
