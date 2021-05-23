#include "bpas.h"
//#include "cpu_timer.h"

///////////////////////////////////////////////////////////

/* test_dft_six_step_small_prime_field_cpu:
 * compute FFT-based multiplication of two
 * polynomials of degree at most N/2 (N=K^e)
 * over small prime field Z/(p_u64)Z.
 */
int
test_dft_six_step_small_prime_field_cpu (int K, int e, usfixn64 p_u64)
{

  /////////////////////////////////////
  //// check if support for the base-case is provided.
  /////////////////////////////////////
  int support_DFT_K_available = 0;
  for (int i = 0; i < sizeof(supported_small_prime_DFT_K) / sizeof(int); i++)
    {
      if (K == supported_small_prime_DFT_K[i])
	{
	  support_DFT_K_available = 1;
	  break;
	}
    }
  if (support_DFT_K_available == 0)
    {
      printf ("ERROR: DFT_{K=%d} is not supported!\n", K);
      exit (EXIT_FAILURE);
    }

  /////////////////////////////////////
  if (sizeof(unsigned long) != sizeof(unsigned long long))
    {
      printf (
	  "-- ERROR:\n"
	  "   On this machine, unsigned long is not as wide as 'unsigned long long'\n"
	  "   This can fail much of the gmp convert in and convert out functions!\n"
	  "   return -1\n");
      return -1;
    }

  /////////////////////////////////////
  //set input vector size N=K^e
  int N = pow (K, e);

  ///////////////////////////////////////
  ///// init montgomery triple
  ///////////////////////////////////////
  montgomery_triple P;
  init_montgomery_triple (&P, p_u64);

  ///////////////////////////////////////////////////////
  ///// LOADING INPUT VECTORS x and y
  ///////////////////////////////////////////////////////
  ///// WARNING: make sure each vector if half full
  ///// (i.e., N/2 elements are 0), otherwise, FFT-based
  ///// polynomial multiplication will FAIL.
  ///////////////////////////////////////////////////////
  int data_size = (N * sizeof(usfixn64));
  usfixn64* x_data = (usfixn64*) malloc (data_size);
  usfixn64* y_data = (usfixn64*) malloc (data_size);

  memset (x_data, 0x00, (data_size));
  memset (y_data, 0x00, (data_size));

  srand(time(NULL));
  usfixn64 rv=rand();
  for (int i = 0; i < N / 2; i++)
    {
      x_data[i] = rv;//p_u64 - 1 - i; //(i % vector_size);
      y_data[i] = rv-1;//p_u64 - 2 - i; //(i % vector_size);
      rv-=2;
    }
//  printf("GENERATING RANDOM DATA!!!\n");

  /////////////////////////////////////
  printf ("%s", shortline);
  printf (
      "[computing small_prime_field DFT_N=K^e, N=%d, K=%d, e=%d, p=%llu] ...\n",
      N, K, e, p_u64 - 1);
  printf ("%s", shortline);
  /////////////////////////////////////
  usfixn64 omega, omega_inv, n_inv;

  //computing N'th root of unity in small prime field (omega).
  compute_nth_root_of_unity_for_small_prime (p_u64, &omega, N);
  //computing inverse of N'th root of unity (omega_inv).
  gmp_inv_mod_p_u64 (&omega_inv, omega, p_u64);
  //computing inverse of N (N^{-1})
  gmp_inv_mod_p_u64 (&n_inv, N, p_u64);

  usfixn64* x_data_copy, *y_data_copy;
  int verification_enabled = VERIFICATION_ENABLED;
  if (verification_enabled)
    {
#if VERBOSE >=1
      printf (
	  "[copying x_data and y_data which will be used for verification] ...\n");
      printf("%s",shortline);
#endif
      x_data_copy = (usfixn64*) malloc (data_size);
      y_data_copy = (usfixn64*) malloc (data_size);

      memcpy (x_data_copy, x_data, data_size);
      memcpy (y_data_copy, y_data, data_size);
    }

  ///////////////////////////////////////
  ///// convert in input data to Montgomery represetnation.
  ///////////////////////////////////////

  for (int i = 0; i < N; i++)
    {
      convertIn_GLOBAL_ptr (&x_data[i], &P);
    }
  for (int i = 0; i < N; i++)
    {
      convertIn_GLOBAL_ptr (&y_data[i], &P);
    }

  ///////////////////////////////////////
  convertIn_GLOBAL_ptr (&omega, &P);
  convertIn_GLOBAL_ptr (&omega_inv, &P);

  //setting up timers.
  cpu_timer t_fft_based_mult;
  timer_record_start (&t_fft_based_mult);

  //x=DFT_n(x)
  DFT_general_small_elements (x_data, K, e, &P, omega);

//  //y=DFT_n(y)
  DFT_general_small_elements (y_data, K, e, &P, omega);
//
//  //x=DFT_n(x)*DFT_n(y)
  convolution_mult_small_elements (x_data, y_data, N, &P);

  //x=DFT^{-1}_n(x)
  DFT_general_small_elements (x_data, K, e, &P, omega_inv);
  ///////////////////////////////////////
#if VERBOSE>=1
  printf ("[computing n_inv = N^{-1} for inverse DFT]...");
#endif

#if VERBOSE>=1
  printf ("DONE\n");
  printf("%s",shortline);
#endif
  ///////////////////////////////////////
  convertIn_GLOBAL_ptr (&n_inv, &P);
  ///////////////////////////////////////

#if VERBOSE>=1
  printf ("[multiplying convolution result by N^{-1}]...");
#endif

  //x={1/n}*x
  for (int i = 0; i < N; i++)
    {
      x_data[i] = mult_ab_mod_p (x_data[i], n_inv, &P);
    }
  ///////////////////////////////////////
  timer_record_stop (&t_fft_based_mult);
  printf ("%s", shortline);
  timer_get_elapsed_time (&t_fft_based_mult, "t_fft_based_mult_total", 1);
  ///////////////////////////////////////

  for (int i = 0; i < N; i++)
    {
      convertOut_GLOBAL_ptr (&x_data[i], &P);
    }

#if VERBOSE>=1
  printf ("DONE\n");
  printf("%s",shortline);
#endif

  ///////////////////////////////////////
  int verification_result = EXIT_SUCCESS;
  if (verification_enabled)
    {
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
      gmp_karatsuba_poly_mult_u64 (x_data_copy, y_data_copy, N, p_u64 - 1, 1);

      int test_failed = 0;
      printVectorToFile (x_data, N, 1, "x.tmp", 0);
      printVectorToFile (x_data_copy, N, 1, "xcopy.tmp", 0);
      for (int i = 0; i < N; i++)
	{
	  if (x_data_copy[i] != x_data[i])
	    {
	      printf ("ERROR: MISMATCH at [%d] ...\n", i);
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
	  exit (EXIT_FAILURE);
	}
      printf ("%s", shortline);
    }

  ///////////////////////////////////////
  free (x_data);
  free (y_data);
  if (verification_enabled)
    {
      free (x_data_copy);
      free (y_data_copy);
    }
  return verification_result;
}

///////////////////////////////////////////////////////////

int
main (int argc, char ** argv)
{

//  test_cilk(argc, argv);
//  return EXIT_SUCCESS;
  setbuf (stdout, NULL);
  int K = 16;
  int e = 2;
  int fft_based_mult_enabled = 1;

  if ((argc == 2))
    if ((strcmp (argv[1], "-h") == 0) || (strcmp (argv[1], "--help") == 0))
      {
	printf ("%s", shortline);
	printf ("computing DFT_{K^e}:"
		"args: [1]:K  [2]:e \n");
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

  test_dft_six_step_small_prime_field_cpu (K, e, fourier_primes_u64_table[0]);

//  //VERIFIED.
//  for (int K = 16; K <= 64; K *= 2)
//    for (int e = 2; e <= 3; e++)
//	 test_dft_six_step_small_prime_field_cpu(K, e, fourier_primes_u64_table[0]);

  return 0;
}

///////////////////////////////////////////////////////////
