#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

typedef unsigned long long int usfixn64;

//usfixn64 conv_p1 = 9223336852482686977;
//usfixn64 conv_p2 = 9220451733971402753;

usfixn64 conv_p1 = 4611615649683210241;
usfixn64 conv_p2 = 4610208274799656961;

void
compute_srgfn_p_gmp (mpz_t p, usfixn64 radix, int k)
{
  mpz_t m;
  mpz_t r;

  mpz_init_set_ui (m, 0);
  mpz_init_set_ui (r, radix);

  mpz_pow_ui (m, r, k);
  mpz_add_ui (p, m, 1);
  mpz_clear (m);
}

/**************************************/
int
compute_primes (int k_log2)
{
  int k = 1 << k_log2;

  usfixn64 radix;

  int stat;
  int n_rep = 50;

//  for (int u3 = 32; u3 < 64; u3++)
//    for (int u2 = 32; u2 < 64; u2++)
//      for (int u1 = 1; u1 < u2; u1++)
//	for (int u0 = 0; u0 < u1; u0++)
//	  {
//	    mpz_t p_zz;
//	    mpz_init (p_zz);
//	    radix = (1UL << u3) + (1UL << u2) + (1UL << u1) + (1UL << u0);
////	printf("radix=%llu,k=%d\n", radix,k);
//	    compute_srgfn_p_gmp (p_zz, radix, k);
//	    stat = mpz_probab_prime_p (p_zz, n_rep);
//	    if (stat != 0)
//	      printf ("PRIME [%d,%d, %d, %d]\n", u3, u2, u1, u0);
//	    mpz_clear (p_zz);
//	  }

  char file_name[4096];
  sprintf (file_name, "srgfn_prime_k%d.h", k);
  FILE * output_file = fopen (file_name, "w");
  setbuf (stdout, NULL);

  printf ("[generating srgfn's for k=%d]...\n", k);
  printf ("writing result to [%s]\n", file_name);

  char buffer[8192];
  sprintf (buffer, "#ifndef SRGFN_k%d_\n"
	   "#define SRGFN_k%d_\n\n"
	   "static srgfn_prime p_list_%d[10] ={\n",
	   k, k, k_log2);
  fprintf (output_file, "%s", buffer);

  mpz_t p_zz, p1p2_zz, kr2_zz;
  mpz_inits (p_zz, p1p2_zz, kr2_zz, NULL);

  mpz_set_ui(p1p2_zz,conv_p1);
  mpz_mul_ui(p1p2_zz,p1p2_zz, conv_p2);
  mpz_sub_ui(p1p2_zz,p1p2_zz, 1);
  mpz_div_ui(p1p2_zz, p1p2_zz, 2);

  int limit = 10;
  int cnt = 0;
  for (int u2 = 62; u2 >= 32; u2--)
    for (int u1 = u2 - 1; u1 >= 1; u1--)
      for (int u0 = u1 - 1; u0 >= 0; u0--)
	{

	  if (cnt>2*limit)
	    break;
	  radix = (1UL << u2) + (1UL << u1) + (1UL << u0);
//	printf("radix=%llu,k=%d\n", radix,k);
	  compute_srgfn_p_gmp (p_zz, radix, k);
	  stat = mpz_probab_prime_p (p_zz, n_rep);
	  if (stat != 0)
	    {

	      mpz_set_ui(kr2_zz, radix);
	      mpz_pow_ui(kr2_zz, kr2_zz, 2);
	      mpz_mul_ui(kr2_zz, kr2_zz, k);
	      mpz_mul_ui(kr2_zz, kr2_zz, 2);
//	      mpz_div_ui(kr2_zz, kr2_zz, 2);
	      if (mpz_cmp(kr2_zz, p1p2_zz)>=0)
		continue;
//		printf ("PRIME [%d,%d, %d]^%d\n", u2, u1, u0,k);
	      printf (".");

	      if (cnt < limit)
		{
		  if (cnt != 0)
		    fprintf (output_file, ",");
		  fprintf (output_file,
			   "{%d, (1UL<<%d)+(1UL<<%d)+(1UL<<%d)}//idx=%d\n", k,
			   u2, u1, u0, cnt);
		}
	      else
		{
		  fprintf (output_file,
			   "//,{%d, (1UL<<%d)+(1UL<<%d)+(1UL<<%d)}//idx=%d\n",
			   k, u2, u1, u0, cnt);
		}
	      cnt++;
	    }

	}

  mpz_clears(p_zz, p1p2_zz, kr2_zz, NULL);
  fprintf (output_file, "};\n#endif\n");

  fflush (output_file);
  fclose (output_file);
  printf ("\n");
  return EXIT_SUCCESS;
}

/**************************************/

void
generate_srgfn (int argc, char**argv)
{
  int k_log2 = 4;
  if (argc > 1)
    k_log2 = atoi (argv[1]);
  for (k_log2 = 6; k_log2 < 10; k_log2++)
//  for (k_log2 = 1; k_log2 < 7; k_log2++)
    compute_primes (k_log2);

}

/**************************************/

void
find_two_largest_u64_fourier_primes ()
{
  //c.2^k+1;

  int stat;
  mpz_t p_zz;
  mpz_init (p_zz);

  usfixn64 max_p1 = 0, max_p2 = 0;
  int n_rep = 20;

  int lower_k_log2 = 45;

  for (int k = 62; k >= lower_k_log2; k--)
    {
      usfixn64 base = (1UL << k);
      int max_c = (1UL << (62 - k));
      for (int c = max_c; c >= 1; c--)
	{
	  usfixn64 p = c * base + 1;
	  mpz_set_ui (p_zz, p);
	  stat = mpz_probab_prime_p (p_zz, n_rep);
	  if (stat != 0)
	    {
//	      printf (".");
	      if (p > max_p1)
		max_p1 = p;
	    }
	}
    }
  printf ("\n");
  printf ("max_p1=%llu\n", max_p1);

  for (int k = 63; k >= lower_k_log2; k--)
    {
      usfixn64 base = (1UL << k);
      int max_c = (1UL << (63 - k));
      for (int c = max_c; c >= 1; c--)
	{
	  usfixn64 p = c * base + 1;
	  mpz_set_ui (p_zz, p);
	  stat = mpz_probab_prime_p (p_zz, n_rep);
	  if (stat != 0)
	    {
	      if (p >= max_p1)
		continue;
	      //	      printf (".");
	      if (p > max_p2)
		{
		  max_p2 = p;
		}
	    }
	}
    }
  printf ("max_p2=%llu\n", max_p2);

}

/**************************************/
int
main (int argc, char ** argv)
{
//  find_two_largest_u64_fourier_primes ();
  generate_srgfn (argc, argv);
  return 0;
}

