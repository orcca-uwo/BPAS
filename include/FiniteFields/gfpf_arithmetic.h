#ifndef GFPF_ARITHMETIC_H_
#define GFPF_ARITHMETIC_H_

#include "gfpf_arithmetic_decl.h"

crt_u192_data t_crt_data_global;
crt_u192_data *t_crt_data_global_ptr;    //=&t_crt_data_global;

lhc_u192_data t_lhc_data_global;
lhc_u192_data *t_lhc_data_global_ptr;    //=&t_crt_data_global;

/**************************************/
montgomery_triple * global_P1, *global_P2;
usfixn64 *global_conv_omega1, *global_conv_omega2;
usfixn64 *global_conv_omega1_inv, *global_conv_omega2_inv;
// usfixn64 *global_conv_omega1_inv, *global_conv_omega2_inv;
usfixn64 *global_conv_theta1, *global_conv_theta2;
usfixn64 *global_conv_theta1_inv, *global_conv_theta2_inv;
usfixn64 global_K1_inv, global_K2_inv;

/**************************************/

//CHECKED.
//int
//verify_mult_u64_u64_gmp (const usfixn64 a, const usfixn64 b, const usfixn64 s0,
//			 const usfixn64 s1)
//{
//  mpz_t x_zz, b_zz, s0_zz, s1_zz, u64_zz;
//  mpz_inits (x_zz, b_zz, s0_zz, s1_zz, u64_zz, NULL);
//
//  mpz_set_u64 (x_zz, a);
//  mpz_set_u64 (b_zz, b);
//  mpz_mul (x_zz, x_zz, b_zz);
//
//  mpz_set_u64 (u64_zz, (usfixn64) (U64_MASK));
//  mpz_add_ui (u64_zz, u64_zz, 1);
//
//  mpz_tdiv_qr (s1_zz, s0_zz, x_zz, u64_zz);
//
//  int verification_result = EXIT_SUCCESS;
//  if (mpz_get_u64 (s0_zz) != s0)
//    {
//      printf ("mismatch in s0!\n");
//      printf ("s0     =%llu\n", s0);
//      printf ("s0[gmp]=%llu\n", mpz_get_u64 (s0_zz));
//      verification_result = (EXIT_FAILURE);
//    }
//  if (mpz_get_u64 (s1_zz) != s1)
//    {
//      printf ("mismatch in s1!\n");
//      printf ("s1=%llu\n", s1);
//      printf ("s1[gmp]=%llu\n", mpz_get_u64 (s1_zz));
//      verification_result = (EXIT_FAILURE);
//    }
//
//  mpz_clears (x_zz, b_zz, s0_zz, s1_zz, u64_zz, NULL);
//  return verification_result;
//}

/**************************************/

//CHECKED.
int
verify_mult_u64_u128_gmp (const usfixn64 a, const usfixn64 b0,
			  const usfixn64 b1, const usfixn64 s0,
			  const usfixn64 s1)
{
  mpz_t x_zz, t_zz, b_zz, s0_zz, s1_zz, u64_zz;
  mpz_inits (x_zz, t_zz, b_zz, s0_zz, s1_zz, u64_zz, NULL);

  //u64
  mpz_set_u64 (u64_zz, (usfixn64) (U64_MASK));
  mpz_add_ui (u64_zz, u64_zz, 1);

  //b=b0+b1.u64
  mpz_set_u64 (b_zz, b1);
  mpz_mul (b_zz, b_zz, u64_zz);
  mpz_set_u64 (t_zz, b0);
  mpz_add (b_zz, b_zz, t_zz);

  //x=a.b
  mpz_set_u64 (x_zz, a);
  mpz_mul (x_zz, x_zz, b_zz);

  //x -> (s0,s1)
  mpz_tdiv_qr (s1_zz, s0_zz, x_zz, u64_zz);

  int verification_result = EXIT_SUCCESS;
  if (mpz_get_u64 (s0_zz) != s0)
    {
      printf ("mismatch in s0!\n");
      printf ("s0     =%llu\n", s0);
      printf ("s0[gmp]=%llu\n", mpz_get_u64 (s0_zz));
      verification_result = EXIT_FAILURE;
    }
  if (mpz_get_u64 (s1_zz) != s1)
    {
      printf ("mismatch in s1!\n");
      printf ("s1=%llu\n", s1);
      printf ("s1[gmp]=%llu\n", mpz_get_u64 (s1_zz));
      verification_result = EXIT_FAILURE;
    }
  mpz_inits (x_zz, t_zz, b_zz, s0_zz, s1_zz, u64_zz, NULL);
  return verification_result;
}

/**************************************/

//CHECKED.
int
verify_sub_u128_u128_gmp (const usfixn64 a0, const usfixn64 a1,
			  const usfixn64 b0, const usfixn64 b1,
			  const usfixn64 s0, const usfixn64 s1)
{
//u64
  mpz_t u64_zz;
  mpz_init (u64_zz);
  mpz_set_u64 (u64_zz, (usfixn64) (U64_MASK));
  mpz_add_ui (u64_zz, u64_zz, 1);

  //u128=u64^2
  mpz_t u128_zz;
  mpz_init (u128_zz);
  mpz_pow_ui (u128_zz, u64_zz, 2);

  //b=b0+b1.u64
  mpz_t b_zz;
  mpz_init (b_zz);
  mpz_set_u64 (b_zz, b1);
  mpz_mul (b_zz, b_zz, u64_zz);

  //b+=b0
  mpz_t b0_zz;
  mpz_init (b0_zz);
  mpz_set_u64 (b0_zz, b0);
  mpz_add (b_zz, b_zz, b0_zz);

  //a=a1*u64;
  mpz_t a_zz;
  mpz_init (a_zz);
  mpz_set_u64 (a_zz, a1);
  mpz_mul (a_zz, a_zz, u64_zz);

  //a+=a0
  mpz_t a0_zz;
  mpz_init (a0_zz);
  mpz_set_u64 (a0_zz, a0);
  mpz_add (a_zz, a_zz, a0_zz);

  //x=a-b
  mpz_sub (a_zz, a_zz, b_zz); //this number is signed.
  mpz_mod (a_zz, a_zz, u128_zz); //get the unsigned result.

  //x -> (s0,s1)
  mpz_t s0_zz, s1_zz;
  mpz_init (s0_zz);
  mpz_init (s1_zz);

  mpz_tdiv_qr (s1_zz, s0_zz, a_zz, u64_zz);

  int result = 0;
  if (mpz_get_u64 (s0_zz) != s0)
    {
      printf ("mismatch in s0!\n");
      printf ("s0     =%llu\n", s0);
      printf ("s0[gmp]=%llu\n", mpz_get_u64 (s0_zz));
      result = -1;
    }
  if (mpz_get_u64 (s1_zz) != s1)
    {
      printf ("mismatch in s1!\n");
      printf ("s1=%llu\n", s1);
      printf ("s1[gmp]=%llu\n", mpz_get_u64 (s1_zz));
      result = -1;
    }

  mpz_clears (a_zz, a0_zz, b_zz, b0_zz, s0_zz, s1_zz, u64_zz, u128_zz, NULL);
  return result;
}

/**************************************/

//CHECKED.
int
verify_add_u128_u128_gmp (const usfixn64 a0, const usfixn64 a1,
			  const usfixn64 b0, const usfixn64 b1,
			  const usfixn64 s0, const usfixn64 s1,
			  const usfixn64 s2)
{
  //u64
  mpz_t u64_zz, b_zz, x_zz, t_zz, s0_zz, s1_zz, s2_zz;
  mpz_inits (u64_zz, b_zz, x_zz, t_zz, s0_zz, s1_zz, s2_zz, NULL);
  mpz_set_u64 (u64_zz, (usfixn64) (U64_MASK));
  mpz_add_ui (u64_zz, u64_zz, 1);

  //b=b0+b1.u64
  mpz_set_u64 (b_zz, b1);
  mpz_mul (b_zz, b_zz, u64_zz);
  mpz_set_u64 (t_zz, b0);
  mpz_add (b_zz, b_zz, t_zz);

  mpz_set_u64 (x_zz, a1);
  mpz_mul (x_zz, x_zz, u64_zz);
  mpz_set_u64 (t_zz, a0);
  mpz_add (x_zz, x_zz, t_zz);

  //x=a-b
  mpz_add (x_zz, x_zz, b_zz);

  //x -> (s0,s1)
  mpz_tdiv_qr (s2_zz, s0_zz, x_zz, u64_zz);
  mpz_tdiv_qr (s2_zz, s1_zz, s2_zz, u64_zz);

  int verification_result = EXIT_SUCCESS;
  if (mpz_get_ui (s0_zz) != s0)
    {
      printf ("mismatch in s0!\n");
      printf ("s0     =%llu\n", s0);
      printf ("s0[gmp]=%llu\n", mpz_get_u64 (s0_zz));
      verification_result = EXIT_FAILURE;
    }
  if (mpz_get_ui (s1_zz) != s1)
    {
      printf ("mismatch in s1!\n");
      printf ("s1=%llu\n", s1);
      printf ("s1[gmp]=%llu\n", mpz_get_u64 (s1_zz));
      verification_result = EXIT_FAILURE;
    }

  if (mpz_get_ui (s2_zz) != s2)
    {
      printf ("mismatch in s2!\n");
      printf ("s2=%llu\n", s2);
      printf ("s2[gmp]=%llu\n", mpz_get_u64 (s2_zz));
      verification_result = EXIT_FAILURE;
    }

  mpz_clears (u64_zz, b_zz, x_zz, t_zz, s0_zz, s1_zz, s2_zz, NULL);
  return verification_result;
}

/**************************************/

//CHECKED.
void __attribute__ ((noinline))
mult_u64_u64 (const usfixn64 *a, const usfixn64 *b, usfixn64 *s0_out,
	      usfixn64 *s1_out)
{
  usfixn64 s0, s1;
  s0 = 0;
  s1 = 0;

  __asm__ __volatile__(
      "movq  %2, %%rax;\n\t"          // rax = a
      "movq  %3, %%rdx;\n\t"// rax = a
      "mulq  %%rdx;\n\t"// rdx:rax = a * b
      "movq  %%rax, %0;\n\t"// s0 = rax
      "movq  %%rdx, %1;\n\t"// s1 = rdx
      : "=&q" (s0),"=&q"(s1)
      : "q"(*a), "q"(*b)
      : "%rax", "%rdx", "memory");

  *s0_out = s0;
  *s1_out = s1;
//#if VERIFICATION_ENABLED == 1
//  int status = verify_mult_u64_u64_gmp (*a, *b, s0, s1);
//  print_verification_msg ("mult_u64_u64", status);
//#endif
}

/**************************************/

//CHECKED.
//compute a*(b0+b1.u64);
//results are correct only if the whole product
//is less than u128.
void __attribute__ ((noinline))
mult_u64_u128 (const usfixn64 *a, const usfixn64 *b0, const usfixn64 *b1,
	       usfixn64* s0, usfixn64 *s1)
{
////	equivalent version using gcc128
//	__int128 mult = (__int128) a * (__int128) b;
//	s0 = mult & (U64_MASK);
//	s1 = mult >> 64;
  //(s0,s1)=a.b0+a.b1.u64
  __asm__ __volatile__(
      "movq  %2, %%rax;\n\t"          // rax = a
      "movq  %3, %%rdx;\n\t"// rax = a
      "mulq  %%rdx;\n\t"// rdx:rax = a * b
      "movq  %%rax, %0;\n\t"// s0 = rax
      "movq  %%rdx, %1;\n\t"// s1 = rdx

      "movq  %2, %%rax;\n\t"// rax = a
      "movq  %4, %%rdx;\n\t"// rax = a
      "mulq  %%rdx; \n\t"
      "addq  %%rax, %1;\n\t"
      : "+&q" (*s0),"+&q"(*s1)
      : "q"(*a), "q"(*b0), "q"(*b1)
      : "%rax", "%rdx", "memory", "cc");

#if VERIFICATION_ENABLED == 1
  //(s0,s1)=(s0,s1)+a.b1.u64
  int status = verify_mult_u64_u128_gmp (*a, *b0, *b1, *s0, *s1);
  print_verification_msg ("mult_u64_u128", status);
#endif
}

/**************************************/

//CHECKED.
//(x0,x1)=(x0,x1)-(y0,y1);
void __attribute__ ((noinline))
sub_u128_u128 (const usfixn64 *y0, const usfixn64 *y1, usfixn64* x0,
	       usfixn64 *x1)
{
////	equivalent version using gcc128
//	__int128 x= (__int128)(*x1);
//	x<<=64;
//	x+=(__int128)(*x0);
//      __int128 y= (__int128)(*y1);
//	y<<=64;
//	y+=(__int128)(*y0);
//	x-=y;
//	*x0=(usfixn64)(x);
//	x>>=64;
//	*x1=(usfixn64)(x);
//	return;
//	__asm__ (
//			"movq  %2, %%rax;\n\t"          // rax = a
//			"mulq  %3;\n\t"// rdx:rax = a * b
//			"movq  %%rax, %0;\n\t"// s0 = rax
//			"movq  %%rdx, %1;\n\t"// s1 = rdx
//			: "=&rm" (s0),"=&rm"(s1)
//			: "rm"(a), "rm"(b)
//			: "%rax", "%rdx");

  //ab0+ab1u64

  //(s0,s1)=a.b0

  usfixn64 a0, a1;
  a0 = *x0;
  a1 = *x1;

//  printf("x0=%llu, x1=%llu\n",*x0,*x1);
//  printf("y0=%llu, y1=%llu\n",*y0,*y1);
  usfixn64 m0, m1;
//  m0 = 0 - 1 - *y0;
//  m1 = 0 - 1 - *y1;
  m0 = U64_MASK - (*y0);
  m1 = U64_MASK - (*y1);
  __asm__ volatile(
      "addq  %2, %0;\n\t"          // rax = a
      "adcq  %3, %1;\n\t"// rax = a
      "addq  $1, %0;\n\t"// rax = a
      "adcq  $0, %1;\n\t"// rax = a
      : "+&qm" (*x0),"+&qm"(*x1)
      : "q"(m0), "q"(m1)
      : "memory", "cc");

//  printf("res:x0=%llu, x1=%llu\n",*x0,*x1);

//#if VERIFICATION_ENABLED == 1
//  //(s0,s1)=(s0,s1)+a.b1.u64
//  int status = verify_sub_u128_u128_gmp (a0, a1, *y0, *y1, *x0, *x1);
//  print_verification_msg ("sub_u128_u128", status);
//#endif
}

/**************************************/

//CHECKED.
//(x0,x1,x2)+=(y0,y1);
//void __attribute__ ((noinline))
//add_u128_u128 (usfixn64 *x0, usfixn64 *x1, usfixn64 *x2, const usfixn64 *y0,
//	       const usfixn64 *y1)
//{
//  usfixn64 a0, a1, a2;
//  a0 = *x0;
//  a1 = *x1;
//  a2 = *x2;
//  __asm__ volatile (
//      "addq  %3, %0;\n\t"          //r0+=y0
//      "adcq  %4, %1;\n\t"
//      "adcq  $0, %2;\n\t"
//      : "+&q" (a0),"+&q"(a1),"+&q"(a2)
//      : "q"(*y0), "q"(*y1)
//      : "memory");
//
//#if VERIFICATION_ENABLED == 1
//  int status = verify_add_u128_u128_gmp (*x0, *x1, *y0, *y1, a0, a1, a2);
//  print_verification_msg ("add_u128_u128", status);
//#endif
//
//  *x0 = a0;
//  *x1 = a1;
//  *x2 = a2;
//
//}
/**************************************/

//CHECKED.
// computes (a_u64)*(b0_u64, b1_u64)
// a=a_u64;
// b=b0_u64 + (b1_u64)*u64
void
u64_mod_u64 (usfixn64 *a, const usfixn64 n)
{

//	a = a % n;
//	return;
  double ninv = 1 / (double) n;
  usfixn64 q = (usfixn64) ((((double) *a)) * ninv);
  usfixn64 res;
  res = (*a) - q * n;
  *a = res & (U64_MASK);
}

/**************************************/

//CHECKED.
//inv_p_mod_u128: returns [p_inv_m, p_inv_q];
// p_inv = p_inv_m + p_inv_q.u64  = int(u128/p);
void
inv_p_mod_u128 (const usfixn64 * p, usfixn64 * p_inv_m, usfixn64 * p_inv_q)
{
  mpz_t u128_zz, u64_zz, q_zz, m_zz, t_zz;
  mpz_inits (u128_zz, u64_zz, q_zz, m_zz, t_zz, NULL);

  mpz_set_u64 (u64_zz, U64_MASK);
  mpz_add_ui (u64_zz, u64_zz, 1);
  mpz_mul (u128_zz, u64_zz, u64_zz);

  mpz_set_u64 (t_zz, *p);
  mpz_tdiv_q (u128_zz, u128_zz, t_zz);
  mpz_tdiv_qr (q_zz, m_zz, u128_zz, u64_zz);
  *p_inv_m = mpz_get_u64 (m_zz);
  *p_inv_q = mpz_get_u64 (q_zz);

  mpz_clears (u128_zz, u64_zz, q_zz, m_zz, t_zz, NULL);
}

/**************************************/

//CHECKED.
// - compute [a1m2p1+a2m1p2] and return [u0,u1,u2]
int
verify_crt_mult_sub_u192_with_reduction_gmp (const usfixn64 a1,
					     const usfixn64 a2,
					     crt_u192_data data,
					     const usfixn64 s0,
					     const usfixn64 s1)
{
  usfixn64 p1 = data.p1;
  usfixn64 p2 = data.p2;

  usfixn64 m1 = 1;          //data.m1;
  usfixn64 m2 = 1;          //data.m2;

  /*  compute the following using gmp:
   * 	s=((a1*m2)%p1)*p2+((a2*m1)%p2)*p1
   *  rewrite s==s0+s1.u64 as two machine words (s0,s1)
   */

  mpz_t t0, t1, p1p2_zz, u64, t_zz, u0, u1;
  mpz_inits (t0, t1, p1p2_zz, u64, t_zz, u0, u1, NULL);

  //u64=2^64
  mpz_set_u64 (u64, U64_MASK);
  mpz_add_ui (u64, u64, 1);

//	char * tmp = (char*) malloc(1024);

  //p1p2=p1*p2
  mpz_set_u64 (p1p2_zz, p1);
  //	mpz_get_str(tmp, 10, p1p2);
  ////	printf("p1p2=%s\n", tmp);

  mpz_set_u64 (t_zz, p2);
  mpz_mul (p1p2_zz, p1p2_zz, t_zz);

  //	mpz_get_str(tmp, 10, p1p2);
  //	printf("p1p2=%s\n", tmp);

  mpz_set_u64 (t0, a2);
  mpz_set_u64 (t_zz, m2);
  mpz_mul (t0, t0, t_zz);
  mpz_set_u64 (t_zz, p1);
  mpz_mul (t0, t0, t_zz);

  //	mpz_get_str(tmp, 10, t0);
  //	printf("a1m2p2=%s\n", tmp);

  mpz_set_u64 (t1, a1);
  mpz_set_u64 (t_zz, m1);
  mpz_mul (t1, t1, t_zz);
  mpz_set_u64 (t_zz, p2);
  mpz_mul (t1, t1, t_zz);

  //	mpz_get_str(tmp, 10, t1);
  //	printf("a2m1p1=%s\n", tmp);

  mpz_add (t0, t0, t1);
  mpz_mod (t0, t0, p1p2_zz);

  //	mpz_get_str(tmp, 10, t0);
  //	printf("t0+t1=%s\n", tmp);

  mpz_tdiv_qr (u1, u0, t0, u64);

  //	printf("s0=%lu \n", s0);
  //	printf("u0=%lu \n", mpz_get_ui(u0));
  //	printf("s1=%lu \n", s1);
  //	printf("u1=%lu \n", mpz_get_ui(u1));

//	mpz_get_str(tmp, 10, t0);
  //	printf("t0=%s\n", tmp);

  int verification_result = EXIT_SUCCESS;
  if (mpz_get_u64 (u0) != s0)
    {
      printf ("mismatch in s0!\n");
      printf ("s0     =%llu\n", s0);
      printf ("s0[gmp]=%llu\n", mpz_get_u64 (u0));
      verification_result = EXIT_FAILURE;
    }

  if (mpz_get_u64 (u1) != s1)
    {
      printf ("mismatch in s1!\n");
      printf ("s1     =%llu\n", s1);
      printf ("s1[gmp]=%llu\n", mpz_get_u64 (u1));
      verification_result = EXIT_FAILURE;
    }

  mpz_clears (t0, t1, p1p2_zz, u64, t_zz, u0, u1, NULL);
  return verification_result;
}

/**************************************/

//CHECKED.
int
verify_div_by_const_R_gmp (const usfixn64 *x0_u64, const usfixn64 *x1_u64,
			   const usfixn64 * q, const usfixn64 *m, usfixn64 R)
{

  mpz_t x_zz, t_zz, q_zz, m_zz;
  mpz_inits (x_zz, t_zz, q_zz, m_zz, NULL);

  mpz_set_u64 (x_zz, U64_MASK);
  mpz_add_ui (x_zz, x_zz, 1);

  mpz_set_u64 (t_zz, *x1_u64);
  mpz_mul (x_zz, x_zz, t_zz);
  mpz_set_u64 (t_zz, *x0_u64);
  mpz_add (x_zz, x_zz, t_zz);

  mpz_set_u64 (t_zz, R);
  mpz_tdiv_qr (q_zz, m_zz, x_zz, t_zz);

  int verification_result = EXIT_SUCCESS;
  if (mpz_get_u64 (q_zz) != *q)
    {
      printf ("mismatch q!\n");
      printf ("q =%llu\n", *q);
      printf ("q[gmp]=%llu\n", mpz_get_u64 (q_zz));
      verification_result = EXIT_FAILURE;
    }

  if (mpz_get_u64 (m_zz) != *m)
    {
      printf ("mismatch m!\n");
      printf ("m =%llu\n", *m);
      printf ("m[gmp]=%llu\n", mpz_get_u64 (m_zz));
      verification_result = EXIT_FAILURE;
    }

  mpz_clears (x_zz, t_zz, q_zz, m_zz, NULL);

  return verification_result;

}

/**************************************/

//CHECKED.
void __attribute__ ((noinline))
div_by_const_R (const usfixn64 *x0_u64, const usfixn64 *x1_u64,
		const usfixn64 *r0, const usfixn64 *r1, usfixn64 *q,
		usfixn64 *m, usfixn64 R)
{
//
//  __int128 x =(__int128)(*x1_u64);
//  x<<=64;
//  x+=(__int128)(*x0_u64);
//  __int128 rr=(__int128)(*r1);
//    rr<<=64;
//    rr+=(__int128)(*r0);

//    printf("=========================\n");
//    printf("r0=%llu, r1=%llu\n", *r0,*r1);

//	r_inv= (u128/r_in);
//	r1,r0=[r_inv/u64, r_inv%u64];
//	x1,x0=[x/u64, x%u64];
//	v0=x0*r0;
//	v1=x0*r1;
//	v2=x1*r0;

//	__int128 v1, v2, q0;
//	usfixn64 x0 = x0_u64;
//	usfixn64 x1 = x1_u64;
//	v0 = 0;
//	v1 = 0;
//	v2 = 0;

  usfixn64 v0_lo = 0, v0_hi = 0;
  usfixn64 v1_lo = 0, v1_hi = 0;
  usfixn64 v2_lo = 0, v2_hi = 0;

//	v0 = (__int128) x0 * (__int128) r0;
  mult_u64_u64 (x0_u64, r0, &v0_lo, &v0_hi);
//	v1 = (__int128) x0 * (__int128) r1;
//	v2 = (__int128) x1 * (__int128) r0;

  mult_u64_u64 (x0_u64, r1, &v1_lo, &v1_hi);
  mult_u64_u64 (x1_u64, r0, &v2_lo, &v2_hi);

//	v0 >>= 64;
//	v1 += (v0_hi);
//	v2 += v1;
//	v2 >>= 64;

  // (v1_lo,v1_hi)+=(v0_hi) //add with carry
//	__asm__ (
//			"movq  %2, %%rax;\n\t"          // rax = a
//			"addq  %%rax, %0;\n\t"// rdx:rax = a * b
//			"adcq  $0x0, %1;\n\t"
//			: "+&rm"(v1_lo),"+&rm"(v1_hi)
//			: "rm" (v0_hi)
//			: "%rax");
  __asm__ volatile(
      "addq  %2, %0;\n\t"		// rdx:rax = a * b
      "adcq  $0x0, %1;\n\t"
      : "+q"(v1_lo),"+q"(v1_hi)
      : "q" (v0_hi)
      : "memory");

  // (v2_lo, v2_hi) += (v1_lo, v1_hi);
//	__asm__ volatile(
//			"movq  %2, %%rax;\n\t"          // rax = a
//			"addq  %%rax, %0;\n\t"// rdx:rax = a * b
//			"adcq  $0x0, %1;\n\t"
//			"movq  %3, %%rax;\n\t"// rax = a
//			"addq  %%rax, %1;\n\t"// rdx:rax = a * b
//			: "+q"(v2_lo),"+q"(v2_hi)
//			: "q" (v1_lo),"q" (v1_hi)
//			: "%rax");

  __asm__ volatile(
      "addq  %2, %0;\n\t"		// rdx:rax = a * b
      "adcq  %3, %1;\n\t"
      : "+rm"(v2_lo),"+rm"(v2_hi)
      : "rm" (v1_lo),"rm" (v1_hi)
      :"memory");

//	q0 = 0;
//	q0 = (__int128) x1 * (__int128) r1;

  usfixn64 q0_hi = 0, q0_lo = 0;

  mult_u64_u64 (x1_u64, r1, &q0_lo, &q0_hi);

  //at this point, only v2_hi is required.
//	q0 += v2;
//	q0+=v2_hi;

//	__asm__ (
//			"movq  %2, %%rax;\n\t"          // rax = a
//			"addq  %%rax, %0;\n\t"// rdx:rax = a * b
//			"adcq  $0x0, %1;\n\t"
//			: "+&rm"(q0_lo),"+&rm"(q0_hi)
//			: "rm" (v2_hi)
//			: "%rax");
  __asm__ volatile (
      "addq  %2, %0;\n\t"		// rdx:rax = a * b
      "adcq  $0x0, %1;\n\t"
      : "+q"(q0_lo),"+q"(q0_hi)
      : "q" (v2_hi)
      : "memory");

//	printf("(q0_l0,hi)+v2= %lu , %lu\n", q0_lo, q0_hi);

  usfixn64 x0 = *x0_u64;
  usfixn64 x1 = *x1_u64;
//	m0=x-(q0*r_in);
//	__int128 m0;
//	m0 = (__int128)(1L << 64);
//	m0 *= (__int128) x1;
//	m0 += (__int128) x0;
//
//	__int128 q0;
//	q0 = (__int128)(1L << 64);
//	q0 *= (__int128) q0_hi;
//	q0 += (__int128) q0_lo;
//	m0 = m0 - q0 * (__int128)(R);

  usfixn64 m0_lo, m0_hi;
  x0 = q0_lo;
  x1 = q0_hi;
  mult_u64_u128 (&R, &q0_lo, &q0_hi, &m0_lo, &m0_hi);

  x0 = *x0_u64;
  x1 = *x1_u64;
  sub_u128_u128 (&m0_lo, &m0_hi, &x0, &x1);

//	usfixn64 m0_lo, m0_hi;

//	if (m0 >= (__int128)(R))
//	{
//		//	printf("carry\n");
//		m0 -= (__int128) R;
//		q0 += 1;
//	}

  if ((x1 >= 1) || ((x1 == 0) && (x0 >= R)))
    {
//#if VERIFICATION_ENABLED ==1
//      printf ("carry\n");
//#endif

      usfixn64 zero = 0;
      sub_u128_u128 (&R, &zero, &x0, &x1);
      q0_lo++;
    }

//	m = (usfixn64) (m0 & U64_MASK);
//	q = (usfixn64) (q0 & (U64_MASK));

  *m = x0;
  *q = q0_lo;
//
//  if ((q0_hi) > 0)
//    {
//      printf ("WARNING: q >= u64!\n");
//      exit (EXIT_FAILURE);
//      ;
//    }

#if VERIFICATION_ENABLED == 1
  printf ("%s", "[verifying div_by_const_R for]...");
  if (verify_div_by_const_R_gmp (x0_u64, x1_u64, q, m, R) != EXIT_SUCCESS)
    {
      printf ("FAILED!\n");
      exit (EXIT_FAILURE);
    }
  else
    {
      printf ("VERIFIED!\n");
    }
  printf ("%s", shortline);
#endif

}

/**************************************/

//CHECKED.
//div_values = [x0,x1][rinv_0,rinv_1,R]
//const usfixn64 *x0_u64, const usfixn64 *x1_u64,
//		const usfixn64 *r0, const usfixn64 *r1, usfixn64 *q,
//		usfixn64 *m, usfixn64 R
void __inline__
div_by_const_R_ptr (usfixn64* div_values)
{
  usfixn64 q0, q1;

  //q0+=high(r0.x0)
  __asm__ __volatile__(
      "movq %1, %%rax;\n\t"
      "mulq %2;\n\t"
      "movq %%rdx, %0;\n\t"
      :"=g"(q0)
      :"g"(div_values[0]),"g"(t_lhc_data_global_ptr->r_inv_0)
      :"%rax","%rdx");

  //(q0,q1)+=high(x0.r1+x1.r0)
  __asm__ __volatile__(
      "movq %1, %%rax;\n\t"
      "mulq %2;\n\t"
      "movq %%rax, %%rsi;\n\t"
      "movq %%rdx, %%rdi;\n\t"

      "movq %3, %%rax;\n\t"
      "mulq %4;\n\t"

      "addq %%rsi, %%rax;\n\t"
      "adcq %%rdi, %%rdx;\n\t"

      "addq %0, %%rax;\n\t"
      "adcq $0x0, %%rdx;\n\t"

      "movq %%rdx, %0;\n\t"
      :"+g"(q0)
      :"g"(div_values[0]),"g"(t_lhc_data_global_ptr->r_inv_1), "g"(div_values[1]),"g"(t_lhc_data_global_ptr->r_inv_0)
      :"%rax","%rdx", "%rsi", "%rdi");

  //s0+s1.u64] -
  __asm__ __volatile__(
      "movq %2, %%rax;\n\t"
      "mulq %3;\n\t"

      "addq %0, %%rax;\n\t"
      "adcq $0x0, %%rdx;\n\t"

      "movq %%rax, %0;\n\t"
      "movq %%rdx, %1;\n\t"
      :"+g"(q0), "=g"(q1)
      :"g"(div_values[1]),"g"(t_lhc_data_global_ptr->r_inv_1)
      :"%rax","%rdx");

//  printf ("x0=%llu, x1=%llu, q0=%llu, q1=%llu\n", div_values[0], div_values[1], q0 , q1);

  usfixn64 s0, s1;

  __asm__ __volatile__(
      "movq %2, %%rax;\n\t"
      "mulq %4;\n\t"

      "movq %%rax, %%rsi;\n\t"
      "movq %%rdx, %%rdi;\n\t"

      "movq %3, %%rax;\n\t"
      "mulq %4;\n\t"

      "addq %%rdi, %%rax;\n\t"
      "adcq $0x0, %%rdx;\n\t"

      "movq %%rsi, %0;\n\t"
      "movq %%rax, %1;\n\t"

      :"=g"(s0), "=g"(s1)
      :"g"(q0), "g"(q1), "g"(t_lhc_data_global_ptr->radix)
      :"%rax","%rdx","%rsi", "%rdi");

//  s0 = U64_MASK - s0;
//  s1 = U64_MASK - s1;
  __asm__ __volatile__(

      "notq %0;\n\t"
      "notq %1;\n\t"
      "addq %2, %0;\n\t"
      "adcq %3, %1;\n\t"
      "addq $0x1, %0;\n\t"
      "adcq $0x0, %1;\n\t"
      :"+g"(s0), "+g"(s1)
      :"g"(div_values[0]), "g"(div_values[1])
      :);

  if ((s1 >= 1) || ((s1 == 0) && (s0 >= t_lhc_data_global_ptr->radix)))
    {
//#if VERIFICATION_ENABLED ==1
//      printf ("carry\n");
//#endif
      //      usfixn64 zero = 0;
      //      sub_u128_u128 (&R, &zero, &x0, &x1);
      __asm__ __volatile__(
	  "addq  %2, %0;\n\t"          // rax = a
	  "adcq  %3, %1;\n\t"// rax = a
	  "addq  $1, %0;\n\t"// rax = a
	  "adcq  $0, %1;\n\t"// rax = a
	  : "+g" (s0),"+g"(s1)
	  : "g"(U64_MASK-t_lhc_data_global_ptr->radix), "g"(U64_MASK)
	  : );
      q0++;
    }

//  printf ("div_by_const_r; x0=%llu; x1=%llu\n", div_values[0], div_values[1]);
//  mpz_t q_zz, m_zz, result_zz, r_zz;
//  mpz_inits (q_zz, m_zz, result_zz, r_zz, NULL);
//
//  mpz_set_u64 (r_zz, t_lhc_data_global_ptr->radix);
//
//  mpz_set_u64 (result_zz, div_values[1]);
//  mpz_mul_2exp (result_zz, result_zz, 64);
//  mpz_add_ui (result_zz, result_zz, div_values[0]);
//
//  mpz_tdiv_qr (q_zz, m_zz, result_zz, r_zz);
//
//  div_values[0] = mpz_get_u64 (m_zz);
//  div_values[1] = mpz_get_u64 (q_zz);
//
////  printf ("gmp; x0=%llu; x1=%llu\n", div_values[0], div_values[1]);
//  mpz_clears (q_zz, m_zz, result_zz, r_zz, NULL);
//
//  return;

  div_values[0] = s0;
  div_values[1] = q0;
}

/**************************************/

//CHECKED.
//div_values = [x0,x1][rinv_0,rinv_1,R]
//const usfixn64 *x0_u64, const usfixn64 *x1_u64,
//		const usfixn64 *r0, const usfixn64 *r1, usfixn64 *q,
//		usfixn64 *m, usfixn64 R
void __inline__          //__attribute__ ((noinline))
div_by_const_R_ptr_single_digit (usfixn64* div_values)
{

//  r0 = t_crt_data->r_inv_0;
//  r1 = t_crt_data->r_inv_1;
//  R = t_crt_data->radix;

  ///////////////////////////////////////
  //  __int128 x =(__int128)(*x1_u64);
  //  x<<=64;
  //  x+=(__int128)(*x0_u64);
  //  __int128 rr=(__int128)(*r1);
  //    rr<<=64;
  //    rr+=(__int128)(*r0);
  ///////////////////////////////////////
  //	r_inv= (u128/r_in);
  //	r1,r0=[r_inv/u64, r_inv%u64];
  //	x1,x0=[x/u64, x%u64];
  //	v0=x0*r0;
  //	v1=x0*r1;
  //	v2=x1*r0;

  usfixn64 q0_lo;

  __asm__ __volatile__(
      "movq  %2, %%rdi;\n\t"          // rax = a
      "movq  %3, %%rax;\n\t"// rax = a
      "mulq  %%rdi;\n\t"// rdx:rax = a * b
      "movq  %%rdx, %%rsi;\n\t"// s1 = rdx
      "movq  %4, %%rax;\n\t"// rax = a
      "mulq  %%rdi;\n\t"// rdx:rax = a * b
      "addq  %%rsi, %%rax;\n\t"// s0 = rax
      "adcq  $0x0, %%rdx;\n\t"// s1 = rdx
      "movq  %%rdx, %0;\n\t"// rax = a
      "movq  %5, %%rax;\n\t"// rax = a
      "mulq %%rdx;\n\t"
      "notq %%rax;\n\t"// rax = a
      "notq %%rdx;\n\t"// rax = a
      "addq  %2, %%rax;\n\t"// rax = a
      "adcq  $0, %%rdx;\n\t"// rax = a
      "addq  $1, %%rax;\n\t"// rax = a
      "adcq  $0, %%rdx;\n\t"// rax = a
      "movq  %%rax, %2;\n\t"// rax = a
      "movq  %%rdx, %1;\n\t"// rax = a
      : "=g"(q0_lo), "+g"(div_values[1]), "+g"(div_values[0])
      : "g"(t_lhc_data_global_ptr->r_inv_0), "g"(t_lhc_data_global_ptr->r_inv_1),"g"(t_lhc_data_global_ptr->radix)
      : "%rax", "%rdx", "%rdi", "%rsi");


  if ((div_values[1] >= 1)
      || ((div_values[1] == 0)
	  && (div_values[0] >= t_lhc_data_global_ptr->radix)))
    {
#if VERIFICATION_ENABLED ==1
      printf ("carry\n");
#endif
      __asm__ __volatile__(
	  "movq %2, %%rdi;\n\t"
	  "notq  %%rdi;\n\t"          // rax = a
	  "addq  %0, %%rdi;\n\t"// rax = a
	  "adcq  %3, %1;\n\t"// rax = a
	  "addq  $1, %%rdi;\n\t"// rax = a
	  "movq %%rdi, %0;\n\t"// rax = a
	  : "+g" (div_values[0]),"+g"(div_values[1])
	  : "g"(t_lhc_data_global_ptr->radix), "g"(U64_MASK)
	  : "%rdi");
      q0_lo++;
    }

//  if ((q0_hi) > 0)
//    {
//      printf ("WARNING: q >= u64!\n");
//      exit (EXIT_FAILURE);
//    }

//  div_values[0] = div_values[0];
  div_values[1] = q0_lo;
}

/**************************************/

//lhc_0+=lhc_r mod r;
void __inline__
add_lhc_ptr (usfixn64 *__restrict__ lhc_0, const usfixn64* __restrict__ lhc_1,
	     const usfixn64* __restrict__ lhc_2, usfixn64 * __restrict__ R)
{
  usfixn64 c_in = 0, c = 0, s;
#pragma unroll
  for (int i = 0; i < 3; i++)
    {
      __asm__ __volatile__(
//      "xorq  %0, %0;\n\t" //s=0
//      "xorq  %1, %1;\n\t"//c=0
	  "movq  %6, %%rdi;\n\t"//s=i0
	  "movq  %2, %%rsi;\n\t"//s=i0
	  "addq  %3, %%rsi;\n\t"//s=i0+i1
	  "adcq  $0, %1;\n\t"//c+=carry;
	  "addq  %5, %%rsi;\n\t"//s=i0+i1
	  "adcq  $0, %1;\n\t"//c+=carry;
	  "cmpq  %%rdi, %%rsi;\n\t"//c+=carry;
	  "jl L0_%=;\n\t"//c+=carry;
	  "subq %%rdi, %%rsi;\n\t"
	  "inc %1;\n\t"
	  "L0_%=:;\n\t"
	  "addq  %4, %%rsi;\n\t"//s=i0+i1
	  "adcq  $0, %1;\n\t"//c+=carry;
	  "cmpq  %%rdi, %%rsi;\n\t"//c+=carry;
	  "jl L1_%=;\n\t"//c+=carry;
	  "subq %%rdi, %%rsi;\n\t"
	  "inc %1;\n\t"
	  "L1_%=:;\n\t"
	  "movq %%rsi, %0;\n\t"
	  : "=g" (s),"+g"(c):
	  "g"(lhc_0[i]), "g"(lhc_1[i]),"g"(lhc_2[i]), "g"(c_in),"g"(*R)
	  : "%rsi", "%rdi");
      lhc_0[i] = s;
      c_in = c;
      c = 0;
    }

//  mpz_t l_zz, h_zz, c_zz;
//  mpz_inits (l_zz, h_zz, c_zz, NULL);
//  mpz_set_u64 (l_zz, lhc_0[0]);
//  mpz_add_ui (l_zz, l_zz, lhc_1[0]);
//  mpz_add_ui (l_zz, l_zz, lhc_2[0]);
//
//  mpz_set_u64 (h_zz, lhc_0[1]);
//  mpz_add_ui (h_zz, h_zz, lhc_1[1]);
//  mpz_add_ui (h_zz, h_zz, lhc_2[1]);
//
//  mpz_set_u64 (c_zz, lhc_0[2]);
//  mpz_add_ui (c_zz, c_zz, lhc_1[2]);
//  mpz_add_ui (c_zz, c_zz, lhc_2[2]);
//
//
//  mpz_mul_ui(h_zz, h_zz, *R);
//  mpz_mul_ui(c_zz, c_zz, *R);
//  mpz_mul_ui(c_zz, c_zz, *R);
//  mpz_add(l_zz, l_zz, h_zz);
//  mpz_add(l_zz, l_zz, c_zz);
//
//  mpz_tdiv_qr_ui(c_zz, l_zz, l_zz, *R);
//  mpz_tdiv_qr_ui(c_zz, h_zz, c_zz, *R);
//
//  lhc_0[0]=mpz_get_u64(l_zz);
//  lhc_0[1]=mpz_get_u64(h_zz);
//  lhc_0[2]=mpz_get_u64(c_zz);
}

/**************************************/

//CHECKED.
void __inline__
crt_mult_sub_u192_with_reduction (usfixn64 *a1, usfixn64 *a2, int k,
				  const crt_u192_data * t_crt_data,
				  usfixn32 * sign_u32)
{

  //// backup gmp implementation.
  usfixn64 s0, s1;
  usfixn64 t[2];
  usfixn64 m0, m1;
//  usfixn64 data_m1, data_m2;
  const usfixn64 one = 1;
//  data_m1 = t_crt_data_global_ptr->m1_mont;
//  data_m2 = t_crt_data_global_ptr->m2_mont;

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int j = 0; j < k; j++)
    {
      t[0] = a1[j];
      t[1] = a2[j];

//  data_m1 = convertIn_GLOBAL (data_m1, global_P1);
//  data_m2 = convertIn_GLOBAL (data_m2, global_P2);
//  t[0] = convertIn_GLOBAL (t[0], global_P1);
//  t[2] = convertIn_GLOBAL (t[2], global_P2);
//  convertIn_GLOBAL_ptr (&t[0], global_P1);
//  convertIn_GLOBAL_ptr (&t[2], global_P2);

//a1m1%p1 and a2m2%p2
      mult_ab_mod_p_ptr (&t[0], &t_crt_data_global_ptr->m1_mont, global_P1);
      mult_ab_mod_p_ptr (&t[1], &t_crt_data_global_ptr->m2_mont, global_P2);

      s0 = t[0];
      s1 = t[1];
////
//  *s0 = convertOut_GLOBAL (*s0, global_P1);
//  *s1 = convertOut_GLOBAL (*s1, global_P2);
////  convert out = montgomery_mult(a.1)

//      mult_ab_mod_p_ptr (&s0, &one, global_P1);
//      mult_ab_mod_p_ptr (&s1, &one, global_P2);

      convertOut_GLOBAL_ptr (&s0, global_P1);
      convertOut_GLOBAL_ptr (&s1, global_P2);

//      mpz_t s0_zz, s1_zz, t_zz, p1p2_zz, p1p2_half_zz;
//      mpz_inits (s0_zz, s1_zz, t_zz, p1p2_zz,p1p2_half_zz, NULL);
//      mpz_set_u64 (s0_zz, s0);
//      mpz_set_u64 (s1_zz, s1);
//      mpz_set_u64 (t_zz, t_crt_data_global_ptr->p2);
//      mpz_mul (s0_zz, s0_zz, t_zz);
//      mpz_set_u64 (t_zz, t_crt_data_global_ptr->p1);
//      mpz_mul (s1_zz, s1_zz, t_zz);
//      mpz_add(s0_zz, s0_zz, s1_zz);
//
//      mpz_set_u64(p1p2_zz, t_crt_data_global_ptr->p1);
//      mpz_set_u64(t_zz, t_crt_data_global_ptr->p2);
//      mpz_mul(p1p2_zz, p1p2_zz, t_zz);
//
//      mpz_mod(s0_zz, s0_zz, p1p2_zz);
//
//      mpz_set(p1p2_half_zz, p1p2_zz);
//      mpz_sub_ui(p1p2_half_zz, p1p2_half_zz, 1);
//      mpz_div_ui(p1p2_half_zz, p1p2_half_zz, 2);
//
//      sign_vector[j]=0;
//      if (mpz_cmp(s0_zz, p1p2_half_zz)>=0)
//	{
//	  sign_vector[j]=1;
//	  mpz_sub(s0_zz, p1p2_zz, s0_zz);
//	  printf("negative sign at [%d]\n",j);
//	}
//
//      mpz_set_u64(t_zz, 1);
//      mpz_mul_2exp(t_zz, t_zz, 64);
//
//      mpz_tdiv_qr(s1_zz, s0_zz, s0_zz, t_zz);
//      a1[j]=mpz_get_u64(s0_zz);
//      a2[j]=mpz_get_u64(s1_zz);
//
//      mpz_clears (s0_zz, s1_zz, t_zz, p1p2_zz,p1p2_half_zz, NULL);
////      return;
//      continue;

////  (a1m1).p2 = s0.p2 -> store in t[0,1];
////  mult_u64_u64 (s0, &data->p2, &t[0], &t[1]);
//  __asm__ __volatile__(
//      "movq  %2, %%rax;\n\t"          // rax = a
//      "movq  %3, %%rdx;\n\t"// rax = a
//      "mulq  %%rdx;\n\t"// rdx:rax = a * b
//      "movq  %%rax, %0;\n\t"// s0 = rax
//      "movq  %%rdx, %1;\n\t"// s1 = rdx
//      : "=q" (t[0]),"=q"(t[1])
//      : "q"(s0), "q"(data->p2)
//      : "%rax", "%rdx");
//
//  ///////////////////////////
//
//  //(a2m2).p1 = s1.p1 -> store in t[2,3]
////  mult_u64_u64 (s1, &data->p1, &t[2], &t[3]);
//  __asm__ __volatile__(
//      "movq  %2, %%rax;\n\t"          // rax = a
//      "movq  %3, %%rdx;\n\t"// rax = a
//      "mulq  %%rdx;\n\t"// rdx:rax = a * b
//      "movq  %%rax, %0;\n\t"// s0 = rax
//      "movq  %%rdx, %1;\n\t"// s1 = rdx
//      : "=q" (t[2]),"=q"(t[3])
//      : "q"(s1), "q"(data->p1)
//      : "%rax", "%rdx");
//
//  //add (t[0,1],carry)+=(t[2,3])
////  usfixn64 carry = 0;
////  add_u128_u128 (&t[0], &t[1], &carry, &t[2], &t[3]);
//
//  __asm__ __volatile__(
//      "addq %2, %0; \n\t"
//      "adcq %3, %1; \n\t;"
//      : "+q" (t[0]),"+q"(t[1])
//      : "q"(t[2]), "q"(t[3])
//      : "cc");

      //  (a1m1).p2 = s0.p2 -> store in t[0,1];
      //  mult_u64_u64 (s0, &data->p2, &t[0], &t[1]);
      //(a2m2).p1 = s1.p1 -> store in t[2,3]
      //  mult_u64_u64 (s1, &data->p1, &t[2], &t[3]);
      __asm__ __volatile__(
	  "movq  %2, %%rax;\n\t"          // rax = a
//	    "movq  %3, %%rdx;\n\t"// rax = a
	  "mulq  %3;\n\t"// rdx:rax = a * b
	  "movq  %%rax, %0;\n\t"// s0 = rax
	  "movq  %%rdx, %1;\n\t"// s1 = rdx

	  "movq  %4, %%rax;\n\t"// rax = a
//	    "movq  %5, %%rdx;\n\t"// rax = a
	  "mulq  %5;\n\t"// rdx:rax = a * b
//	    "addq  %0, %%rax;\n\t"// s0 = rax
//	    "adcq  %1, %%rdx;\n\t"// s1 = rdx
//	    "movq  %%rax, %0;\n\t"// s0 = rax
//	    "movq  %%rdx, %1;\n\t"// s1 = rdx
	  "addq  %%rax, %0;\n\t"// s0 = rax
	  "adcq  %%rdx, %1;\n\t"// s1 = rdx
	  : "+g" (t[0]),"+g"(t[1])
	  : "g"(s0), "g"(t_crt_data_global_ptr->p2),"g"(s1), "g"(t_crt_data_global_ptr->p1)
	  : "%rax", "%rdx");

      ///////////////////////////

      //add (t[0,1],carry)+=(t[2,3])
      //  usfixn64 carry = 0;
      //  add_u128_u128 (&t[0], &t[1], &carry, &t[2], &t[3]);

//  __asm__ __volatile__(
//      "addq %2, %0; \n\t"
//      "adcq %3, %1; \n\t;"
//      : "+q" (t[0]),"+q"(t[1])
//      : "q"(t[2]), "q"(t[3])
//      : "cc");

      if ((t[1] > t_crt_data_global_ptr->p1p2_q)
	  || ((t[1] == t_crt_data_global_ptr->p1p2_q)
	      && (t[0] > t_crt_data_global_ptr->p1p2_m)))
	{
	  m0 = U64_MASK - t_crt_data_global_ptr->p1p2_m;
	  m1 = U64_MASK - t_crt_data_global_ptr->p1p2_q;
	  __asm__ __volatile__(
	      "addq %2, %0; \n\t"
	      "adcq %3, %1; \n\t"
	      "addq $0x1, %0; \n\t"
	      "adcq $0x0, %1; \n\t"
	      : "+g" (t[0]),"+g"(t[1])
	      : "g"(m0), "g"(m1)
	      : );

	}
      //normalizing (i.e., determining sign)

//      sign_vector[j] = 0;
      char sign = 0;
      if ((t[1] > t_crt_data->p1p2_a1)
	  || (t[1] == t_crt_data->p1p2_a1 && t[0] > t_crt_data->p1p2_a0))
	{
	  usfixn64 tmp = t[0] > t_crt_data->p1p2_b0;
	  t[1] = t_crt_data->p1p2_b1 - t[1];
	  t[0] = t_crt_data->p1p2_b0 - t[0] + 1;
	  //	    if (t > t_crt_data->p1p2_b0)
	    {
	      //		s1[j]--;
	      t[1] -= tmp;
	    }
//	  sign_vector[j] = 1;
	  sign = 1;
	}

      a1[j] = t[0];
      a2[j] = t[1];

//      int idx = (j / 32);
//      usfixn32 flag = ((sign == 1) << (j % 32));
      int idx = (j >> 5);
      usfixn32 flag = ((sign == 1) << (j & 0x1f));
      sign_u32[idx] |= flag;

//	s0 = t[0];
//	s1 = t[1];
//
//#if VERIFICATION_ENABLED == 1
//	t[0] = *a1;
//	t[1] = *a2;
//	convertOut_GLOBAL_ptr (&t[0], global_P1);
//	convertOut_GLOBAL_ptr (&t[1], global_P2);
//	int status = verify_crt_mult_sub_u192_with_reduction_gmp (
//	    t[0], t[1], *t_crt_data_global_ptr, s0, s1);
//	print_verification_msg ("crt_mult_sub_u192_with_reduction", status);
//#endif
//	a1[j] = s0;
//	a2[j] = s1;
    }

}

/**************************************/
/**************************************/

//CHECKED.
// void  __attribute__((always_inline))
void __inline__	//__attribute__((noinline))
lhc_by_R_u128_ptr (usfixn64 * l_vec, usfixn64 * h_vec, usfixn64 * c_vec,
		   usfixn32 * sign_u32, int k)
{

//  mpz_t x_zz, t_zz, r_zz;
//  mpz_inits (x_zz, t_zz, r_zz, NULL);
//  mpz_set_u64 (x_zz, lhc_points[2]);
//  mpz_mul_2exp (x_zz, x_zz, 64);
//  mpz_set_u64 (t_zz, lhc_points[1]);
//  mpz_add (x_zz, x_zz, t_zz);
//  mpz_mul_2exp (x_zz, x_zz, 64);
//  mpz_set_u64 (t_zz, lhc_points[0]);
//  mpz_add (x_zz, x_zz, t_zz);
//
//  mpz_set_u64 (r_zz, t_lhc_data_global.radix);
//  mpz_tdiv_qr (x_zz, t_zz, x_zz, r_zz);
//  lhc_points[0] = mpz_get_u64 (t_zz);
//  mpz_tdiv_qr (x_zz, t_zz, x_zz, r_zz);
//  lhc_points[1] = mpz_get_u64 (t_zz);
//  lhc_points[2] = mpz_get_u64 (x_zz);
//
//  mpz_clears (x_zz, t_zz, r_zz, NULL);
//  return;

  usfixn64 m1, q1;
  usfixn64 div_values[2];
//  usfixn64 *div_ptrs[3];
//  div_ptrs[0] = r0;
//  div_ptrs[1] = r1;
//  div_ptrs[2] = R_ptr;
//  usfixn64 *div_ptrs=lhc_ptrs;

  usfixn64 l0_lo;
  usfixn64 l1_lo, l1_hi;
  usfixn64 h0_lo;
  usfixn64 h1_lo_hi[2];
  usfixn64 c0_lo, c0_hi;
  usfixn64 c1_lo, c1_hi;

  //  ### should be precomputed
  //  [mb,qb]=div_by_const_R(u64);
  //  div_by_const_R(0, 1, r0, r1, qb, mb);

  //  qb = u64_mod_R_q;
  //  mb = u64_mod_R_m;

//	[m0,q0]=div_by_const_R(x0);
//  div_by_const_R (x0, &zero, r0, r1, &q0, &m0, R);

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int j = 0; j < k; j++)
    {
      usfixn64 lhc_points[3];
      lhc_points[0] = l_vec[j];
      lhc_points[1] = h_vec[j];
//      lhc_points[2] = 0;

      div_values[0] = lhc_points[0];
//  div_values[1] = 0;
      div_by_const_R_ptr_single_digit (div_values);
      l0_lo = div_values[0];
      h0_lo = div_values[1];

//	# q1=x1/R;s
//	# m1=x1%R;
//	[m1,q1]=div_by_const_R(x1);
//  div_by_const_R (x1, &zero, r0, r1, &q1, &m1, R);

      div_values[0] = lhc_points[1];  //
//  div_values[1] = 0;
      div_by_const_R_ptr_single_digit (div_values);
      m1 = div_values[0];
      q1 = div_values[1];

      __asm__ __volatile__(
	  "movq %8, %%rdi;\n\t"
	  "movq  %6, %%rax;\n\t"          // rax = a
	  "mulq  %%rdi;\n\t"// rdx:rax = a * b
	  "movq  %%rax, %0;\n\t"// s0 = rax
	  "movq  %%rdx, %1;\n\t"// s1 = rdx

	  "movq  %7, %%rax;\n\t"// rax = a
	  "mulq  %%rdi;\n\t"// rdx:rax = a * b
	  "movq  %%rax, %2;\n\t"// s0 = rax
	  "movq  %%rdx, %3;\n\t"// s1 = rdx

	  "movq %6, %%rax;\n\t"
	  "mulq %9;\n\t"
	  "addq  %%rax, %2;\n\t"//r0+=y0
	  "adcq  %%rdx, %3;\n\t"

	  "movq  %7, %%rax;\n\t"// rax = a
	  "mulq  %9;\n\t"// rdx:rax = a * b
	  "movq  %%rax, %4;\n\t"// s0 = rax
	  "movq  %%rdx, %5;\n\t"// s1 = rdx
	  : "+g" (l1_lo),"+g"(l1_hi), "+g" (h1_lo_hi[0]),"+g"(h1_lo_hi[1]), "+g" (c1_lo),"+g"(c1_hi)
	  : "g"(m1), "g"(q1),"g"(t_lhc_data_global_ptr->mb), "g"(t_lhc_data_global_ptr->qb)
	  : "%rax", "%rdx", "%rdi");

//      usfixn64 lhc_ans[3] = { l0_lo, h0_lo, c0_lo };
      lhc_points[0] = l0_lo;
      lhc_points[1] = h0_lo;
      lhc_points[2] = c0_lo;

      usfixn64 lhc_l1c1[3];
      usfixn64 lhc_h1h2[3];
      lhc_l1c1[0] = l1_lo;
      lhc_l1c1[1] = l1_hi;
      div_by_const_R_ptr (lhc_l1c1);
      lhc_l1c1[2] = c1_hi;

      div_by_const_R_ptr (h1_lo_hi);

      lhc_h1h2[0] = 0;
      lhc_h1h2[1] = h1_lo_hi[0];
      lhc_h1h2[2] = h1_lo_hi[1];

      add_lhc_ptr (lhc_points, lhc_l1c1, lhc_h1h2,
		   &t_lhc_data_global_ptr->radix);

//  memcpy (lhc_points, lhc_ans, 3 * sizeof(usfixn64));
//      lhc_points[0] = lhc_ans[0];
//      lhc_points[1] = lhc_ans[1];
//      lhc_points[2] = lhc_ans[2];

      int idx = (j >> 5);
      usfixn32 sign = sign_u32[idx] & 0x1;
//      if (sign_vector[j]!=0)
      sign_u32[idx] >>= 1;
      if (sign)
	{
	  lhc_points[0] = -lhc_points[0];
	  lhc_points[1] = -lhc_points[1];
	  lhc_points[2] = -lhc_points[2];

	}
      l_vec[j] = lhc_points[0];
      h_vec[j] = lhc_points[1];
      c_vec[j] = lhc_points[2];
    }

}

/**************************************/

//usfixn64 aux_space[256];
//CHECKED.
void __inline__
oneShiftRight (usfixn64 * xs, int k)
{
  usfixn64 tmp;
  tmp = xs[k - 1];
//  int n_bytes=k*sizeof(usfixn64);
//  usfixn64 * t=aux_space;//(usfixn64*)malloc(n_bytes);
//
//  memcpy(&t[1],xs, (k-1)*sizeof(usfixn64));
//  t[0]=-tmp;
//  memcpy(xs, t, n_bytes);
//  memmove (&xs[1], xs, (k - 1) * sizeof(usfixn64));
//  free(t);
//  for (int b = k - 1; b > 0; b -= CONVOLUTION_CACHE_SIZE)
#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = (k) - 1; i > 0; i--)
//    for (int i = b; i > b - CONVOLUTION_CACHE_SIZE; i--)
//      {
    xs[i] = xs[i - 1];
//	if (i == 0)
//	  break;
//      }
  xs[0] = -tmp;
}

/**************************************/

//CHECKED.
void __inline__
twoShiftRight (usfixn64 * xs, int k)
{
  usfixn64 t[2];
  t[1] = -xs[k - 1];
  t[0] = -xs[k - 2];
//  int n_bytes=k*sizeof(usfixn64);
//
//  usfixn64 * t=aux_space;//(usfixn64*)malloc(n_bytes);

//  memmove (&xs[2], xs, (k - 2) * sizeof(usfixn64));

//  t[0]=-tmp0;
//  t[1]=-tmp1;
//  memcpy(xs, t, n_bytes);

////	for (short i = k - 1; i > 1; i--)
//  for (int b = k - 1; b > 0; b -= 4)
  //	for (short i = k - 1; i > 0; i--)
//    for (int i = b; i > b - 4; i--)
//  for (int b = k - 1; b > 0; b -= CONVOLUTION_CACHE_SIZE)
//    //  for (short i = (k) - 1; i > 0; i--)
//    for (int i = b; i > b - CONVOLUTION_CACHE_SIZE; i--)
//      {
#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = k - 1; i > 1; i--)
    xs[i] = xs[i - 2];
//	if (i == 0)
//	  break;
//      }
//  xs[0] = -tmp0;
//  xs[1] = -tmp1;
//  memcpy (&xs[0], t, (2) * sizeof(usfixn64));
  xs[0] = t[0];
  xs[1] = t[1];
}

/**************************************/

void
normalize_positive_negative_vector (usfixn64 *x, const int k, const usfixn64 r)
{

  usfixn64 c_negative = 0;
  int i;
  short post;

  for (i = 0; i < k; i++)
    {
      if (x[i] >= c_negative)
	{
	  x[i] -= c_negative;
	  c_negative = 0;
	}
      else
	{
	  c_negative = 1;
	  x[i] = r - 1;
	}

      if (x[i] > r)
	{
	  x[i] = r + x[i];
	  c_negative = 1;
	}
    }

  usfixn64 t;
  //c_negative = count of (-1x-1)'s to add == (X+c_negative);
  for (i = 0; i < k; i++)
    {
      t = c_negative;
      c_negative = 0;
      usfixn64 sum = x[i] + t;
      if ((sum >= r))
	{
	  //	  printf ("first carry!\n");
	  c_negative++;
	  sum -= r;
	}
      x[i] = sum;
    }

//    if (c_negative)
//      printf("c_negative = %llu\n", c_negative);
  if (c_negative)
    {
      c_negative--;
      post = -1;

      for (i = 0; i < k; i++)
	{
	  if (x[i] != 0)
	    {
	      post = i;
	      break;
	    }
	}

      if (post >= 0)
	{
	  for (i = 0; i < post; i++)
	    {
	      x[i] = r - 1;
	    }
	  x[post]--;
	}
      else
	{
	  x[k - 1] = r;

	  for (i = 0; i < k - 1; i++)
	    {
	      x[i] = 0;
	    }
	}
    }

//  if (c_negative)
//    printf("c_negative = %llu\n", c_negative);

}

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
//inline void  //__attribute__((always_inline))
void __inline__  //__attribute__((noinline))
addition_hybrid_big_elements (usfixn64 * __restrict__ x,
			      usfixn64 *__restrict__ y,
			      usfixn64 *__restrict__ z, const int k,
			      const usfixn64 r, usfixn64 *c_out,
			      usfixn64*c_negative_out)
{

  usfixn64 sum = 0, t = 0, c_negative = 0, c = 0;
  usfixn64 c0 = 0;  //, c1 = 0, c2 = 0;

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = 0; i < k; i++)
    {
      t = c0;
      c0 = 0;
      if (x[i] > r)
	{
	  x[i] = r + x[i];
	  c0++;
	}

      if (x[i] >= t)
	{
	  x[i] -= t;
	}
      else
	{
	  x[i] = r - t;
	  c0++;
	}
      ////////////////////
      if (y[i] > r)
	{
	  y[i] = r + y[i];
	  c0++;
	}
      ////////////////////
      if (z[i] > r)
	{
	  z[i] = r + z[i];
	  c0++;
	}
      ////////////////////
      t = c;
      c = 0;
      sum = x[i] + y[i];  //+t+z[i];

      if ((sum >= r))
	{
	  c++;
	  sum -= r;
	}
      sum += t;
      if ((sum >= r))
	{
	  c++;
	  sum -= r;
	}

      sum += z[i];
      if ((sum >= r))
	{
	  c++;
	  sum -= r;
	}
      x[i] = sum;
    }
  c_negative = c0;  // + c1 + c2;
  *c_negative_out += c_negative;
  *c_out += c;

  //need conditions on values of r and p
  //must have 2*r <p
  //also, kr^2<p1p2/2 -> or p1p2 to be sure.
}

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void __inline__
addition_big_element_with_single_digit (usfixn64 * x, usfixn64 digit, int plus,
					const int k, const usfixn64 r)
{

  usfixn64 c = 0, sum = 0;
  //need conditions on values of r and p
  //must have 2*r <p
  //also, kr^2<p1p2/2 -> or p1p2 to be sure.

  if (digit == 0)
    {
//      printf ("SINGLE DIGIT IS ZERO; return\n");
      return;
    }

  if (plus == 1)
    {
      c = digit;

#pragma unroll CONVOLUTION_CACHE_SIZE
      for (int i = 0; i < k; i++)
	{
	  sum = x[i] + c;
	  c = 0;
	  if ((sum >= r))
	    {
	      c++;
	      sum -= r;
	    }
	  x[i] = sum;
	}
    }
  else
    {
      c = -digit + 2;
#pragma unroll CONVOLUTION_CACHE_SIZE
      for (int i = 0; i < k; i++)
	{
	  sum = x[i] + r - 1 + c;
	  c = 0;
	  if ((sum >= r))
	    {
	      c++;
	      sum -= r;
	    }
	  x[i] = sum;
	}
    }

  while (c != 0)
    {
      c--;

      int i = 0;
      for (i = 0; i < k; i++)
	{
	  if (x[i] != 0)
	    {
	      x[i]--;
	      break;
	    }
	  x[i] = r - 1;
	}

      if (i == k)
	{
	  x[k - 1] = r;
	  memset (x, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }
}

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void
addition_big_elements_v0 (usfixn64 * x, usfixn64 *y, const int k,
			  const usfixn64 r)
{

  usfixn64 c = 0, sum = 0;
  int i = 0;
  short post = 0;

  //need conditions on values of r and p
  //must have 2*r <p
  //also, kr^2<p1p2/2 -> or p1p2 to be sure.
  usfixn64 t;

  for (i = 0; i < k; i++)
    {
      t = c;
      c = 0;
      sum = x[i] + y[i] + t;
      if ((sum >= r))
	{
//	  printf ("first carry!\n");
	  c++;
	  sum -= r;
	}
      x[i] = sum;

    }

  while (c != 0)
    {
      c--;
      post = -1;

      for (i = 0; i < k; i++)
	{
	  if (x[i] != 0)
	    {
	      post = i;
	      break;
	    }
	}

      if (post >= 0)
	{
	  for (i = 0; i < post; i++)
	    {
	      x[i] = r - 1;
	    }
	  x[post]--;
	}
      else
	{
	  x[k - 1] = r;

	  for (i = 0; i < k - 1; i++)
	    {
	      x[i] = 0;
	    }
	}
    }
}

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void
addition_big_elements (usfixn64 * x, usfixn64 *y, const int k, const usfixn64 r)
{

  usfixn64 c = 0, sum = 0;
  int i = 0;
  short post = 0;

  //need conditions on values of r and p
  //must have 2*r <p
  //also, kr^2<p1p2/2 -> or p1p2 to be sure.
  usfixn64 t;

  for (i = 0; i < k; i++)
    {
      t = c;
      c = 0;
      sum = x[i] + y[i] + t;
      if ((sum >= r))
	{
//	  printf ("first carry!\n");
	  c++;
	  sum -= r;
	}
      x[i] = sum;

    }

  while (c != 0)
    {
      c--;
//      post = -1;

      for (i = 0; i < k; i++)
	{
	  if (x[i] != 0)
	    {
	      x[i]--;
	      break;
	    }
	  x[i] = r - 1;
	}

      if (i == k)
	{
	  x[k - 1] = r;
	  memset (x, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }
}

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void
subtraction_big_elements_v0 (usfixn64 *x, usfixn64 *y, const int k,
			     const usfixn64 r)
{

  int c = 0;
  int post = 0;
  usfixn64 sub = 0;
  int i = 0;
  for (i = 0; i < k; i++)
    {
      sub = y[i] + c;

      if (x[i] < sub)
	{
	  c = 1;
	  x[i] = r - sub + x[i];
	}
      else
	{
	  c = 0;
	  x[i] = x[i] - sub;
	}
    }

  if (c > 0)
    {
      post = -1;
      for (i = 0; i < k; i++)
	{
	  if (x[i] < (r - 1))
	    {
	      post = i;
	      break;
	    }
	}

      if (post >= 0)
	{
	  for (i = 0; i < post; i++)
	    {
	      x[i] = 0;
	    }

	  x[post]++;
	}
      else
	{
	  x[k - 1] = r;

	  for (i = 0; i < k - 1; i++)
	    {
	      x[i] = 0;
	    }
	}
    }
}

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void __inline__
subtraction_big_elements (usfixn64 *x, usfixn64 *y, const int k,
			  const usfixn64 r)
{

  int c = 0;
  int post = 0;
  usfixn64 sub = 0;
  int i = 0;
  for (i = 0; i < k; i++)
    {
      sub = y[i] + c;

      if (x[i] < sub)
	{
	  c = 1;
	  x[i] = r - sub + x[i];
	}
      else
	{
	  c = 0;
	  x[i] = x[i] - sub;
	}
    }

  if (c > 0)
    {
      post = -1;
      for (i = 0; i < k; i++)
	{
	  if (x[i] < (r - 1))
	    {
	      x[i]++;
	      break;
	    }
	  x[i] = 0;
	}

      if (i == k)
	{
	  x[k - 1] = r;
	  memset (x, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }
}

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void __inline__
subtraction_from_zero_big_elements (usfixn64 *y, const int k, const usfixn64 r)
{

  usfixn64 c = 0;
  usfixn64 sub = 0;

  for (int i = 0; i < k; i++)
    {
//      sub = y[i] + c;
//      if (0 < sub)
//	{
//	  c = 1;
//	  sub = r - sub;
//	}
//      else
//	{
//	  c = 0;
//	  sub = -sub;
//	}
//      y[i] = sub;

      sub = y[i] + c;
      c = (sub != 0);
      y[i] = (-c & r) - sub;

    }

  if (c > 0)
    {
      int i;
      for (i = 0; i < k; i++)
	{
	  if (y[i] < (r - 1))
	    {
	      y[i]++;
	      break;
	    }
	  y[i] = 0;
	}

      if (i == k)
	{
	  y[k - 1] = r;
	  memset (y, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }
}

/**************************************/

//CHECKED.
//x=x-y
// subtraction for GFPF
void __inline__
subtraction_upper_minus_lower_big_elements (usfixn64 *x, const int k,
					    const usfixn64 r, int sn)
{
  int c = 0;
  int post = 0;
  usfixn64 sub = 0;
  int i = 0;

  for (i = 0; i < sn; i++)
    {
      sub = x[i] + c;

//      if (0 < sub)
//	{
//	  c = 1;
//	  sub = r - sub;
//	}
//      else
//	{
//	  c = 0;
//	  sub = 0 - sub;
//	}
//      x[i] = sub;

      c = (sub != 0);
      x[i] = (-c & r) - sub;
    }

  for (i = sn; i < k; i++)
    {
      sub = 0 + c;
//      if (x[i] < sub)
//	{
//	  c = 1;
//	  sub = r - sub + x[i];
//	}
//      else
//	{
//	  c = 0;
//	  sub = x[i] - sub;
//	}

      c = (x[i] < sub);
      x[i] = x[i] - sub + ((-c) & r);
    }

  if (c > 0)
    {
      post = -1;
      for (i = 0; i < k; i++)
	{
	  if (x[i] < (r - 1))
	    {
	      x[i]++;
	      break;
	    }
	  x[i] = 0;
	}

      if (i == k)
	{
	  x[k - 1] = r;
	  memset (x, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }
}

/**************************************/

//CHECKED.
//x=x+y
//addtion for GFPF
void __inline__
addition_subtraction_big_elements (usfixn64 * x, usfixn64 *y, const int k,
				   const usfixn64 r)
{

  usfixn64 c = 0, sum;
  usfixn64 c_sub = 0, sub;
  //need conditions on values of r and p
  //must have 2*r <p
  //also, kr^2<p1p2/2 -> or p1p2 to be sure.

  for (int i = 0; i < k; i++)
    {
      sum = x[i] + y[i] + c;
      sub = y[i] + c_sub;
      c = 0; //CRITICAL.
      if ((sum >= r))
	{
	  c = 1;
	  sum -= r;
	}
//      if (x[i] < sub)
//	{
//	  c_sub = 1;
//	  sub = r - sub + x[i];
//	}
//      else
//	{
//	  c_sub = 0;
//	  sub = x[i] - sub;
//	}

      c_sub = (x[i] < sub);
      if (c_sub)
	{
	  sub = r - sub + x[i];
	}
      else
	{
	  sub = x[i] - sub;
	}

//      sub = x[i]-sub+c_sub*r;
////      sub+=((-c_sub)&r);

      x[i] = sum;
      y[i] = sub;
    }

  //post-processing addition.
  if (c != 0)
    {
      int i;
//      c--;
      for (i = 0; i < k; i++)
	{
	  if (x[i] != 0)
	    {
	      x[i]--;
	      break;
	    }
	  x[i] = r - 1;
	}

      if (i == k)
	{
	  x[k - 1] = r;
	  memset (x, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }

  //post-processing subtraction.
  if (c_sub != 0)
    {
      int i;
//      c_sub--;
      for (i = 0; i < k; i++)
	{
	  if (y[i] < (r - 1))
	    {
	      y[i]++;
	      break;
	    }
	  y[i] = 0;
	}

      if (i == k)
	{
	  y[k - 1] = r;
	  memset (x, 0x00, (k - 1) * sizeof(usfixn64));
	}
    }
}

/**************************************/

void __inline__
shift_right_sn (usfixn64*x, int k, usfixn64 r, int s)
{
//  usfixn64 a[128];
//  int n_bytes = k * sizeof(usfixn64);
//  memset (a, 0x0, n_bytes);
//  memcpy (a, &x[k - s], s * sizeof(usfixn64));
//  if (x[k - 1] == r)
//    {
//      a[s - 1] -= r;
//      a[s]++;
//    }

  usfixn64 t[CYCLICSHIFT_CACHE_SIZE];

  int sn = s;
  int shift_bytes = CYCLICSHIFT_CACHE_SIZE * sizeof(usfixn64);
  int shift_bytes_complement = (k - CYCLICSHIFT_CACHE_SIZE) * sizeof(usfixn64);

  while (sn >= CYCLICSHIFT_CACHE_SIZE)
    {
      sn -= CYCLICSHIFT_CACHE_SIZE;
      memcpy (t, &x[k - CYCLICSHIFT_CACHE_SIZE], (shift_bytes));
      memmove (&x[CYCLICSHIFT_CACHE_SIZE], x, shift_bytes_complement);
      memcpy (x, t, (CYCLICSHIFT_CACHE_SIZE) * sizeof(usfixn64));
    }
  if (sn != 0)
    {
      memcpy (t, &x[k - sn], (sn * sizeof(usfixn64)));
      memmove (&x[sn], x, ((k - sn) * sizeof(usfixn64)));
      memcpy (x, t, (sn) * sizeof(usfixn64));
    }

//  for (int i = 0; i < s; i++)
//    {
//      t = x[k - 1];
//      memcpy(&x[1], x, k1_bytes);
//      x[0] = t;
//    }
//
//  if (x[s]==)
//  for (int i = 0; i < (k - s); i++)
//      {
//        b[i + s] = x[i];
//      }
//  memcpy (&x[s], &x[0], (k - s) * sizeof(usfixn64));
//  memcpy (&x[0], a, s * sizeof(usfixn64));

  if (x[s] == r)
    {
      x[s]++;
      x[s - 1] -= r;
    }

}

/**************************************/

//CHECKED.
//x=x*(r^s)
void __inline__
mult_pow_R (usfixn64 *x, int s, const int k, const usfixn64 r,
	    int compute_inverse)
{
#if PROFILING_ENABLED == 1
  inc_profiling_counter (&n_mult_pow_R_called);
#endif

//  usfixn64 *a = (usfixn64*) malloc (n_bytes);
//  usfixn64 *b = (usfixn64*) malloc (n_bytes);

//  memset (b, 0x0, n_bytes);
//  memset (c, 0x0, n_bytes);

  int all_zero = 1;
  for (int i = 0; i < k; i++)
    if (x[i])
      {
	all_zero = 0;
	break;
      }
  if (all_zero == 1)
    return;

  if (compute_inverse == 1)
    s = (k << 1) - s;
//	s = s % (2 * k);
  s = s & ((k << 1) - 1);

  if ((s == 0))
    {
      return;
    }

  if (s == k)
    {
      ////x=zero-x
      subtraction_from_zero_big_elements (x, k, r);
      return;
    }

  if ((s > k) && (s < (2 * k)))
    {
      s = s - k;
      //mult_pow_R(x, s, k ,r);
      ////x=zero-x
      subtraction_from_zero_big_elements (x, k, r);
    }

  shift_right_sn (x, k, r, s);
  subtraction_upper_minus_lower_big_elements (x, k, r, s);
}

/**************************************/

//CHECKED.
//x=x*(r^s)
void
mult_pow_R_v0 (usfixn64 *x, int s, const int k, const usfixn64 r,
	       int compute_inverse)
{
#if PROFILING_ENABLED == 1
  inc_profiling_counter (&n_mult_pow_R_called);
#endif
  int n_bytes = k * sizeof(usfixn64);
  usfixn64 *a = (usfixn64*) calloc (sizeof(usfixn64), k);
  usfixn64 *b = (usfixn64*) calloc (sizeof(usfixn64), k);
  usfixn64 *c = (usfixn64*) calloc (sizeof(usfixn64), k);

//  usfixn64 *a = global_tmp_a;
//  usfixn64 *b = global_tmp_b;
//  usfixn64 *c = global_tmp_c;
  memset (a, 0x0, n_bytes);
  memset (b, 0x0, n_bytes);
  memset (c, 0x0, n_bytes);

  if (compute_inverse == 1)
    s = (2 * k) - s;
//	s = s % (2 * k);
  s = s & ((k << 1) - 1);

//	if ((s == 0) || (s > (k << 1)))
  if ((s == 0))
    {
      return;
    }
//	if
//		{
//			printf("FAILED!!");
//					return;
//		}
  else if (s == k)
    {
      subtraction_big_elements (c, x, k, r);
      memcpy (x, c, n_bytes);
      free (a);
      free (b);
      free (c);
      return;
    }
  else if ((s > k) && (s < (2 * k)))
//	else if ((s > k) && (s < (k << 1)))
    {
      s = s - k;
//		mult_pow_R(x, s, k ,r);
      subtraction_big_elements (c, x, k, r);
      memcpy (x, c, n_bytes);
    }

  for (int i = 0; i < (k - s); i++)
    {
      b[i + s] = x[i];
    }
  for (int i = k - s; i < k; i++)
    {
      a[i - (k - s)] = x[i];
    }
  if (x[k - 1] == r)
    {
      a[s - 1] -= r;
      a[s]++;
    }
  subtraction_big_elements (b, a, k, r);
  memcpy (x, b, n_bytes);
  free (a);
  free (b);
  free (c);
  return;
}

/**************************************/

//CHECKED.
void __attribute__ ((noinline))
conv_k (usfixn64 *x, usfixn64 *y, const int k,
	const usfixn64 *base_case_pow_omega,
	const usfixn64 *base_case_pow_omega_inv,
	const usfixn64 *base_case_pow_theta,
	const usfixn64 *base_case_pow_theta_inv, const usfixn64 * K_inv,
	const montgomery_triple*P)
{

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = 0; i < k; i++)
    {
      mult_ab_mod_p_ptr (&x[i], &base_case_pow_theta[i], P);
    }

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = 0; i < k; i++)
    {
      mult_ab_mod_p_ptr (&y[i], &base_case_pow_theta[i], P);
    }

  void
  (*DFT_K_ptr) (int n, usfixn64* x, const usfixn64 *base_case_omega_vec,
		const montgomery_triple* P);
  if (k == 4)
    {
      DFT_K_ptr = DFT_4;
    }
  else if (k == 8)
    {
      DFT_K_ptr = DFT_8;
    }
  else if (k == 16)
    {
      DFT_K_ptr = DFT_16;
    }
  else if (k == 32)
    {
      DFT_K_ptr = DFT_32;
    }
  else if (k == 64)
    {
      DFT_K_ptr = DFT_64;
    }
  else if (k == 128)
    {
      DFT_K_ptr = DFT_128;
    }
  else if (k == 256)
    {
      DFT_K_ptr = DFT_256;
    }

  DFT_K_ptr (1, x, base_case_pow_omega, P);
  DFT_K_ptr (1, y, base_case_pow_omega, P);

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = 0; i < k; ++i)
    {
      mult_ab_mod_p_ptr (&x[i], &y[i], P);
    }

  DFT_K_ptr (1, x, base_case_pow_omega_inv, P);

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int i = 0; i < k; ++i)
    {
//      x[i] = mult_ab_mod_p (x[i], *K_inv, P);
//      x[i] = mult_ab_mod_p (x[i], base_case_pow_theta_inv[i], P);
//      mult_ab_mod_p_ptr (&x[i], K_inv, P);
      mult_ab_mod_p_ptr (&x[i], &base_case_pow_theta_inv[i], P);
    }
}

/**************************************/

void
get_proper_conv_p1p2 (usfixn64 *conv_p1, usfixn64 *conv_p2, srgfn_prime prime)
{
//  *conv_p1 = 9223336852482686977;
//  *conv_p2 = 9220451733971402753;

  *conv_p1 = 4611615649683210241;
  *conv_p2 = 4610208274799656961;
  return;
  /////
  int fourier_primes_u64_table_size = 10000;
  int n = fourier_primes_u64_table_size;

  //find first prime larger than the radix, p1
  for (int i = 0; i < n; i++)
    {
      *conv_p1 = fourier_primes_u64_table[i];
      if ((prime.radix) < 2 * (*conv_p1))
	{
#if VERBOSE
	  printf ("PASS: END OF SEARCH FOR P1!\n");
#endif
	  break;
	}
    }

  mpz_t p1p2_zz, p1_zz, r2_zz, r3_zz, k_zz;
  mpz_inits (p1p2_zz, p1_zz, r2_zz, r3_zz, k_zz, NULL);
  mpz_set_u64 (p1_zz, *conv_p1);

  mpz_set_u64 (k_zz, prime.k);
  mpz_set_u64 (r2_zz, prime.radix);
  mpz_pow_ui (r2_zz, r2_zz, 2);
  mpz_mul (r2_zz, r2_zz, k_zz);

  mpz_set_u64 (r3_zz, prime.radix);
  mpz_pow_ui (r3_zz, r3_zz, 3);

  usfixn64 t = 0;
  for (int i = 0; i < n; i++)
    {
      if (*conv_p1 != fourier_primes_u64_table[i])
	{
	  mpz_set_u64 (p1p2_zz, fourier_primes_u64_table[i]);
	  mpz_mul (p1p2_zz, p1p2_zz, p1_zz);
	  mpz_div_ui (p1p2_zz, p1p2_zz, 2); //why 4?
	  if ((mpz_cmp (r2_zz, p1p2_zz) >= 0))
	    {
//	      printf ("[ERROR: must have kr^2 < (1/2).p1p2]\n");
//	      printf(".");
	      continue;
	    }

	  else if ((mpz_cmp (p1p2_zz, r3_zz) >= 0))
	    {
//	      printf ("[ERROR: must have p1p2.(1/2) < r^3]\n");
//	      printf(".");
	      continue;
	    }
//	  printf("\n");
	  else
	    {
	      t = fourier_primes_u64_table[i];
#if VERBOSE
	      printf ("PASS: END OF SEARCH FOR P2!\n");
#endif
	      break;
	    }
	}
    }

  if (t == 0)
    {
      printf ("ERROR: FAILED IN SEARCH FOR P2!\n");
      exit (EXIT_FAILURE);
    }
  *conv_p2 = t;
#if VERBOSE
  printf ("[conv_p1=%llu; conv_p2=%llu;]\n", *conv_p1, *conv_p2);
#endif
  mpz_clears (p1p2_zz, p1_zz, r2_zz, r3_zz, k_zz, NULL);

}

/**************************************/

//CHECKED.
void
init_gfpf_mult_data (crt_u192_data *t_crt_data, lhc_u192_data *t_lhc_data,
		     usfixn64 r, const usfixn64 p1, const usfixn64 p2)
{

  mpz_t u128_zz, u64_zz, r_zz, r_inv_0_zz, r_inv_1_zz;
  mpz_inits (u128_zz, u64_zz, r_zz, r_inv_0_zz, r_inv_1_zz, NULL);

  //1. defining u64 =(1<<64) and u128=(1<<128)
  mpz_set_u64 (u64_zz, (U64_MASK));
  mpz_add_ui (u64_zz, u64_zz, 1);
  mpz_pow_ui (u128_zz, u64_zz, 2);

  ///////////////////////////

  //2. r_inv=(u128)/r;
  mpz_set_u64 (r_zz, r);
  mpz_tdiv_q (r_inv_1_zz, u128_zz, r_zz);

  ///////////////////////////
  //3. r_inv=r_inv_0+r_inv_1*u64
  mpz_tdiv_qr (r_inv_1_zz, r_inv_0_zz, r_inv_1_zz, u64_zz);
  t_lhc_data->r_inv_0 = mpz_get_u64 (r_inv_0_zz);
  t_lhc_data->r_inv_1 = mpz_get_u64 (r_inv_1_zz);

//  printf ("[r_inv_0=%llu, r_inv_1=%llu, r=%llu]\n", t_crt_data->r_inv_0,
//	  t_crt_data->r_inv_1, r);

  ///////////////////////////

  //4. computing m1 and m2
  //m1 = inv (p2) mod p1;
  //m2 = inv (p1) mod p2;
  usfixn64 m1, m2;
  gmp_inv_mod_p_u64 (&m1, p2, p1);
  gmp_inv_mod_p_u64 (&m2, p1, p2);
  usfixn64 m1_mont = m1;
  convertIn_GLOBAL_ptr (&m1_mont, global_P1);
  usfixn64 m2_mont = m2;
  convertIn_GLOBAL_ptr (&m2_mont, global_P2);

  ///////////////////////////

  usfixn64 u64_mod_p1, u64_mod_p2;
  usfixn64 p1_inv_q, p1_inv_m, p2_inv_q, p2_inv_m;
  //5.
  //u64_mod_p1 = (2^64)%p1
  //u64_mod_p2 = (2^64)%p2
  u64_mod_p1 = (U64_MASK);
  u64_mod_p2 = (U64_MASK);
  u64_mod_u64 (&u64_mod_p1, p1);
  u64_mod_u64 (&u64_mod_p2, p2);
  u64_mod_p1 += 1;
  u64_mod_p2 += 1;
  inv_p_mod_u128 (&p1, &p1_inv_m, &p1_inv_q);
  inv_p_mod_u128 (&p2, &p2_inv_m, &p2_inv_q);

//  t_crt_data->m1 = m1;
//  t_crt_data->m2 = m2;
  t_crt_data->m1_mont = m1_mont;
  t_crt_data->m2_mont = m2_mont;

  t_crt_data->p1 = p1;
  t_crt_data->p2 = p2;

  t_crt_data->p1_inv_m = p1_inv_m;
  t_crt_data->p1_inv_q = p1_inv_q;
  t_crt_data->p2_inv_m = p2_inv_m;
  t_crt_data->p2_inv_q = p2_inv_q;

  //t_crt_data->p1p2 = p1p2;
  //t_crt_data->p1p2_q = p1p2 >> 64;
  //t_crt_data->p1p2_m = p1p2 & (U64_MASK);

  ///////////////////////////
  //6. set p1p2=m+q.u64;
  mult_u64_u64 (&p1, &p2, &t_crt_data->p1p2_m, &t_crt_data->p1p2_q);

  ///////////////////////////
  usfixn64 zero = 0, one = 1;
  //7. set u64 = q.r + m (this requires r_inv_0 and r_inv_1);
  div_by_const_R (&zero, &one, &t_lhc_data->r_inv_0, &t_lhc_data->r_inv_1,
		  &t_lhc_data->qb, &t_lhc_data->mb, r);
  ///////////////////////////
  //8.
  t_lhc_data->radix = r;

  ///////////////////////////
  //9.
  //
  //p1*p2-1=B=b0+b1.u64
  //(p1*p2-1)/2=A=a0+b1.u64
  // if s>A: then: s=B-s; and sign(s)=negative
  // else: s is positive.
  //normalization step for setting result of crt in range[-p1p2/2, p1p2-1/2]
  mpz_t a0_zz, a1_zz, b0_zz, b1_zz;
  mpz_inits (a0_zz, a1_zz, b0_zz, b1_zz, NULL);
  mpz_set_u64 (a0_zz, t_crt_data->p1);
  mpz_set_u64 (b0_zz, t_crt_data->p2);
  mpz_mul (a0_zz, a0_zz, b0_zz);
  mpz_sub_ui (a0_zz, a0_zz, 1);

  mpz_tdiv_qr (b1_zz, b0_zz, a0_zz, u64_zz);
  mpz_div_ui (a0_zz, a0_zz, 2);
  mpz_tdiv_qr (a1_zz, a0_zz, a0_zz, u64_zz);

  usfixn64 a0, a1, b0, b1;
  a0 = mpz_get_u64 (a0_zz);
  a1 = mpz_get_u64 (a1_zz);
  b0 = mpz_get_u64 (b0_zz);
  b1 = mpz_get_u64 (b1_zz);

  t_crt_data->p1p2_a0 = a0;
  t_crt_data->p1p2_a1 = a1;
  t_crt_data->p1p2_b0 = b0;
  t_crt_data->p1p2_b1 = b1;

  mpz_clears (a0_zz, a1_zz, b0_zz, b1_zz, NULL);
  ///////////////////////////
  mpz_clears (u128_zz, u64_zz, r_zz, r_inv_0_zz, r_inv_1_zz, NULL);

  ////CRITICAL.
  t_crt_data_global_ptr = &t_crt_data_global;
  t_lhc_data_global_ptr = &t_lhc_data_global;
}

/**************************************/

//CHECKED.
void
init_fft_based_bigint_mult (usfixn64 conv_p1, usfixn64 conv_p2, int k)
{
  //fourier_primes_u64_table[0], conv_p2=fourier_primes_u64_table[1];
  usfixn64 omega1, omega2;
  usfixn64 theta1, theta2;
  usfixn64 omega1_inv, omega2_inv;
  usfixn64 theta1_inv, theta2_inv;
  usfixn64 K_inv_mod_p1, K_inv_mod_p2;
  montgomery_triple *P1, *P2;
  P1 = (montgomery_triple*) malloc (sizeof(montgomery_triple));
  P2 = (montgomery_triple*) malloc (sizeof(montgomery_triple));

  init_montgomery_triple (P1, conv_p1);
  init_montgomery_triple (P2, conv_p2);
  usfixn64 * conv_p1_base_case_omega_pow, *conv_p2_base_case_omega_pow;
  usfixn64 * conv_p1_base_case_omega_inv_pow, *conv_p2_base_case_omega_inv_pow;

  usfixn64 * conv_p1_base_case_theta_pow, *conv_p2_base_case_theta_pow;
  usfixn64 * conv_p1_base_case_theta_inv_pow, *conv_p2_base_case_theta_inv_pow;

  //computing N'th root of unity in small prime field (omega).
  compute_nth_root_of_unity_for_small_prime (conv_p1, &theta1, 2 * k);
  compute_nth_root_of_unity_for_small_prime (conv_p2, &theta2, 2 * k);

//  printf ("theta1=%llu\n", theta1);
//  printf ("theta2=%llu\n", theta2);
//  printf ("omega1=%llu\n", omega1);
//  printf ("omega2=%llu\n", omega2);

  //computing inverse of N'th root of unity (omega_inv).
  gmp_inv_mod_p_u64 (&theta1_inv, theta1, conv_p1);
  gmp_inv_mod_p_u64 (&theta2_inv, theta2, conv_p2);

  //computing inverse of N (N^{-1})
  gmp_inv_mod_p_u64 (&K_inv_mod_p1, k, conv_p1);
  gmp_inv_mod_p_u64 (&K_inv_mod_p2, k, conv_p2);
  convertIn_GLOBAL_ptr (&theta1, P1);
  convertIn_GLOBAL_ptr (&theta2, P2);

  convertIn_GLOBAL_ptr (&theta1_inv, P1);
  convertIn_GLOBAL_ptr (&theta2_inv, P2);

  //omega= (theta)^2
  omega1 = mult_ab_mod_p (theta1, theta1, P1);
  omega2 = mult_ab_mod_p (theta2, theta2, P2);

  omega1_inv = mult_ab_mod_p (theta1_inv, theta1_inv, P1);
  omega2_inv = mult_ab_mod_p (theta2_inv, theta2_inv, P2);

  conv_p1_base_case_omega_pow = (usfixn64*) malloc (k * sizeof(usfixn64));
  conv_p2_base_case_omega_pow = (usfixn64*) malloc (k * sizeof(usfixn64));
  conv_p1_base_case_omega_inv_pow = (usfixn64*) malloc (k * sizeof(usfixn64));
  conv_p2_base_case_omega_inv_pow = (usfixn64*) malloc (k * sizeof(usfixn64));

  conv_p1_base_case_theta_pow = (usfixn64*) malloc (k * sizeof(usfixn64));
  conv_p2_base_case_theta_pow = (usfixn64*) malloc (k * sizeof(usfixn64));
  conv_p1_base_case_theta_inv_pow = (usfixn64*) malloc (k * sizeof(usfixn64));
  conv_p2_base_case_theta_inv_pow = (usfixn64*) malloc (k * sizeof(usfixn64));

  precompute_pow_omega (conv_p1_base_case_omega_pow, omega1, k, P1);
  precompute_pow_omega (conv_p2_base_case_omega_pow, omega2, k, P2);
  precompute_pow_omega (conv_p1_base_case_omega_inv_pow, omega1_inv, k, P1);
  precompute_pow_omega (conv_p2_base_case_omega_inv_pow, omega2_inv, k, P2);

  precompute_pow_omega (conv_p1_base_case_theta_pow, theta1, k, P1);
  precompute_pow_omega (conv_p2_base_case_theta_pow, theta2, k, P2);
  precompute_pow_omega (conv_p1_base_case_theta_inv_pow, theta1_inv, k, P1);
  precompute_pow_omega (conv_p2_base_case_theta_inv_pow, theta2_inv, k, P2);

  /////////////////////////////////////
  //// setting the global variables and pointers
  //// needed for the convolution.
  /////////////////////////////////////

  global_P1 = P1;
  global_P2 = P2;

  global_conv_theta1 = conv_p1_base_case_theta_pow;
  global_conv_theta2 = conv_p2_base_case_theta_pow;
  global_conv_theta1_inv = conv_p1_base_case_theta_inv_pow;
  global_conv_theta2_inv = conv_p2_base_case_theta_inv_pow;

  global_conv_omega1 = conv_p1_base_case_omega_pow;
  global_conv_omega2 = conv_p2_base_case_omega_pow;
  global_conv_omega1_inv = conv_p1_base_case_omega_inv_pow;
  global_conv_omega2_inv = conv_p2_base_case_omega_inv_pow;

  convertIn_GLOBAL_ptr (&K_inv_mod_p1, P1);
  convertIn_GLOBAL_ptr (&K_inv_mod_p2, P2);
  global_K1_inv = K_inv_mod_p1;
  global_K2_inv = K_inv_mod_p2;

  for (int i = 0; i < k; i++)
    {
      global_conv_theta1_inv[i] = mult_ab_mod_p (global_conv_theta1_inv[i],
						 global_K1_inv, P1);
    }

  for (int i = 0; i < k; i++)
    {
      global_conv_theta2_inv[i] = mult_ab_mod_p (global_conv_theta2_inv[i],
						 global_K2_inv, P2);
    }

}

/**************************************/

//CHECKED.
void
clear_gfpf_mult_data (crt_u192_data *t_crt_data, lhc_u192_data *t_lhc_data)
{

//  t_crt_data->m1 = 0;
//  t_crt_data->m2 = 0;
  t_crt_data->m1_mont = 0;
  t_crt_data->m2_mont = 0;

  t_crt_data->p1 = 0;
  t_crt_data->p2 = 0;

  t_crt_data->p1_inv_m = 0;
  t_crt_data->p1_inv_q = 0;
  t_crt_data->p2_inv_m = 0;
  t_crt_data->p2_inv_q = 0;

  t_crt_data->p1p2_q = 0;
  t_crt_data->p1p2_m = 0;

  t_crt_data->p1p2_a0 = 0;
  t_crt_data->p1p2_a1 = 0;
  t_crt_data->p1p2_b0 = 0;
  t_crt_data->p1p2_b1 = 0;

  t_lhc_data->qb = 0;
  t_lhc_data->mb = 0;

  t_lhc_data->radix = 0;
  t_lhc_data->r_inv_0 = 0;
  t_lhc_data->r_inv_1 = 0;
}

/**************************************/

//CHECKED.
void
clear_fft_based_bigint_mult_data ()
{
//  free space for the following allocations:
//  montgomery_triple * global_P1, *global_P2;
//  usfixn64 *global_conv_omega1 , *global_conv_omega2;
//  usfixn64 *global_conv_omega1_inv , *global_conv_omega2_inv;
//  usfixn64 *global_K1_inv, *global_K2_inv

  free (global_P1);
  free (global_P2);

  free (global_conv_theta1);
  free (global_conv_theta2);
  free (global_conv_theta1_inv);
  free (global_conv_theta2_inv);

  free (global_conv_omega1);
  free (global_conv_omega2);
  free (global_conv_omega1_inv);
  free (global_conv_omega2_inv);
}

/**************************************/

//usfixn64 aux_memory[256];
//CHECKED.
void
GFPFMultiplication (usfixn64 *x, const usfixn64 *y, int k,
		    const crt_u192_data *t_crt_data,
		    const lhc_u192_data *t_lhc_data)
{

  int n_bytes = k * sizeof(usfixn64);

  /////////////////////////////////////
  //// CRITICAL
  /////////////////////////////////////
  //// disabling this check can slow down
  //// the multiplication significantly!
  /////////////////////////////////////
  int all_zero_x = 1;
  int all_zero_y = 1;
  for (int i = 1; i < k; i++)
    if (x[i] != 0)
      {
	all_zero_x = 0;
	break;
      }

  for (int i = 1; i < k; i++)
    if (y[i] != 0)
      {
	all_zero_y = 0;
	break;
      }

  if (all_zero_x == 1)
    {
      if (x[0] == 0)
	{
	  memset (x, 0x00, n_bytes);
	  return;
	}
      if (x[0] == 1)
	{
	  memcpy (x, y, n_bytes);
	  return;
	}
    }

  if (all_zero_y == 1)
    {
      if (y[0] == 0)
	{
	  memset (x, 0x00, n_bytes);
	  return;
	}
      if (y[0] == 1)
	{
	  return;
	}
    }

  int x_equl_y = 1;
  for (int i = 0; i < k; i++)
    if (x[i] != y[i])
      {
	x_equl_y = 0;
	break;
      }

//  if (x_equl_y==1)
//    printf("BOTH EQUAL!\n");

  /////////////////////////////////////
  //// END OF CRITICAL
  /////////////////////////////////////
  usfixn64 memory[1024];
  usfixn64 *x1 = x;
  usfixn64 *x2 = &memory[0];
  usfixn64 *c_vec = &memory[k];
  usfixn64 *y1 = c_vec;
  usfixn64 *y2 = c_vec;

  //// enough space for 8x32 digits.
  usfixn32 sign_u32[8] =
    { 0, 0, 0, 0, 0, 0, 0, 0 };
  /////////
  memcpy (x2, x, n_bytes);
  memcpy (y2, y, n_bytes);

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int j = 0; j < k; j++)
    convertIn_GLOBAL_ptr (&x2[j], global_P2);

#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int j = 0; j < k; j++)
    convertIn_GLOBAL_ptr (&y2[j], global_P2);

  ////////////////////////////
  ///// convolution
  ////////////////////////////

  conv_k (x2, y2, k, global_conv_omega2, global_conv_omega2_inv,
	  global_conv_theta2, global_conv_theta2_inv, &global_K2_inv,
	  global_P2);

  /////////
  //    memcpy (x1, x, n_bytes);
  memcpy (y1, y, n_bytes);

  ////////////////////////////
  //// convert in to montgomery form
  ////////////////////////////
#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int j = 0; j < k; j++)
    convertIn_GLOBAL_ptr (&x1[j], global_P1);

  /////////
#pragma unroll CONVOLUTION_CACHE_SIZE
  for (int j = 0; j < k; j++)
    convertIn_GLOBAL_ptr (&y1[j], global_P1);

  ////////////////////////////
  ///// convolution
  ////////////////////////////
  conv_k (x1, y1, k, global_conv_omega1, global_conv_omega1_inv,
	  global_conv_theta1, global_conv_theta1_inv, &global_K1_inv,
	  global_P1);

  ////////////////////////////
  //// convert out from montgomery form: this part is moved to CRT function.
  ////////////////////////////
  //  for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
  //    //#pragma unroll LOOP_UNROLLING_ENABLED
  //    for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
  ////      x1[j] = convertOut_GLOBAL (x1[j], global_P1);
  //      convertOut_GLOBAL_ptr (&x1[j], global_P1);
  //
  //  //#pragma unroll LOOP_UNROLLING_ENABLED
  //  for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
  //    for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
  ////      x2[j] = convertOut_GLOBAL (x2[j], global_P2);
  //      convertOut_GLOBAL_ptr (&x2[j], global_P2);
  ////////////////////////////

  /*
   * result of convolution is stored in two vectors x1 and x2.
   * need to combine elements of the same index from the two vector
   * using CRT; CRT(x1_i, x2_i)=(s0_i+s1_i.u64) => store them in s0 and s1;
   * then, convert the (s0,s1) results to lhc;
   */

  //  crt_mult_sub_u192_with_reduction (&x1[j], &x2[j]);
  crt_mult_sub_u192_with_reduction (x1, x2, k, t_crt_data, sign_u32);
  /////////////////////////////////////
  //// normalize CRT result.
  /////////////////////////////////////
  //p1*p2-1=B=b0+b1.u64
  //(p1*p2-1)/2=A=a0+b1.u64
  // if s>A: then: s=B-s; and sign(s)=negative
  // else: s is positive.
  //normalization step for setting result of crt in range[-p1p2/2, p1p2-1/2]
  //the normalization step is integrated into lhc_by_R_u128_ptr.
  lhc_by_R_u128_ptr (x, x2, c_vec, sign_u32, k);

  //  oneShiftRight (h_vec, &k);
  oneShiftRight (x2, k);
  twoShiftRight (c_vec, k);

  usfixn64 r = t_lhc_data->radix;

  //  addition_big_elements (l_vec, h_vec, k, r);
  //  addition_big_elements (l_vec, c_vec, k, r);
  usfixn64 c_out = 0, c_negative_out = 0;
  addition_hybrid_big_elements (x, x2, c_vec, k, r, &c_out, &c_negative_out);
  if (c_out >= c_negative_out)
    {
      c_out -= c_negative_out;
      addition_big_element_with_single_digit (x, c_out, 0, k, r);
    }
  else
    {
      c_negative_out -= c_out;
      addition_big_element_with_single_digit (x, c_negative_out, 1, k, r);
    }
}

/**************************************/
/**************************************/

void
test_fft_based_arbitrary_mult (usfixn64 *x, usfixn64*y, int n,
			       srgfn_prime prime)
{
  crt_u192_data * t_crt_data = &t_crt_data_global;
  lhc_u192_data * t_lhc_data = &t_lhc_data_global;

  int k = prime.k;
  usfixn64 radix = prime.radix;

//  It is assumed that the following functions are already called.
//  init_fft_based_bigint_mult (conv_p1, conv_p2, k);
//  init_gfpf_mult_data (t_crt_data, t_lhc_data, prime.radix, conv_p1, conv_p2);

  mpz_t *a_zz, *b_zz, p_zz, x_zz;
  mpz_inits (p_zz, x_zz, NULL);
  a_zz = (mpz_t*) malloc (n * sizeof(mpz_t));
  b_zz = (mpz_t*) malloc (n * sizeof(mpz_t));
  for (int i = 0; i < n; i++)
    {
      mpz_init (a_zz[i]);
      mpz_init (b_zz[i]);
    }

  usfixn64 *x_gmp = (usfixn64*) malloc (n * k * sizeof(usfixn64));
  compute_srgfn_p_gmp (p_zz, prime.radix, prime.k);

  char msg[256];
  sprintf (msg, "fft_based_arbitrary_mult (k=%d)", k);

  cpu_timer t_gmp, t_convert_in, t_convert_out;

  /////////////////////////////////////
  timer_record_start (&t_convert_in);
  for (int i = 0; i < n; i++)
    {  //A=bigint(a)
      u64vector_to_radix_based_bigint_gmp (a_zz[i], &x[i * k], radix, k);
    }

  for (int i = 0; i < n; i++)
    { //B=bigint(b)
      u64vector_to_radix_based_bigint_gmp (b_zz[i], &y[i * k], radix, k);
    }
  timer_record_stop (&t_convert_in);
  /////////////////////////////////////

  timer_record_start (&t_gmp);
  //A.B mod P
  for (int i = 0; i < n; i++)
    {
      mpz_mul (a_zz[i], a_zz[i], b_zz[i]);
      mpz_mod (a_zz[i], a_zz[i], p_zz);
    }
  timer_record_stop (&t_gmp);

  timer_get_elapsed_time (&t_gmp, "t_arbitrary_mult_gmp", 1);
  timer_get_elapsed_time (&t_gmp, "t_arbitrary_mult_avg_gmp", n);
  /////////////////////////////////////

  timer_record_start (&t_convert_out);
  for (int i = 0; i < n; i++)
    bigint_to_u64vector_radixbased_gmp (&x_gmp[i * k], a_zz[i], radix, k);
  timer_record_stop (&t_convert_out);
  /////////////////////////////////////
  timer_get_elapsed_time (&t_convert_in, NULL, 1);
  timer_get_elapsed_time (&t_convert_out, NULL, 1);

  float t_gmp_based = 0;
  t_gmp_based += t_convert_in.elapsed_time;
  t_gmp_based += t_gmp.elapsed_time;
  t_gmp_based += t_convert_out.elapsed_time;
  timer_print_time (t_gmp_based, "t_arbitrary_mult_gmp_based");
  timer_print_time (t_gmp_based / (1.0 * n), "t_arbitrary_mult_avg_gmp_based");
  /////////////////////////////////////

  cpu_timer t_gfpf;
  timer_record_start (&t_gfpf);
  //(a.b) mod (p)
  for (int i = 0; i < n; i++)
    {
      GFPFMultiplication (&x[i * k], &y[i * k], k, t_crt_data, t_lhc_data);
    }

  timer_record_stop (&t_gfpf);
  timer_get_elapsed_time (&t_gfpf, "t_arbitrary_mult_gfpf", 1);
  timer_get_elapsed_time (&t_gfpf, "t_arbitrary_mult_avg_gfpf", n);
  /////////////////////////////////////
  printf ("%s", longline);
  /////////////////////////////////////
  int status = EXIT_SUCCESS;
  int verification_enabled = VERIFICATION_ENABLED;

  if (verification_enabled)
    {
      for (int i = 0; i < n; i++)
	{
	  for (int j = 0; j < k; j++)
	    if (x[i * k + j] != x_gmp[i * k + j])
	      {
		printf ("[ERROR: mismatch at i,j=%d,%d]\n", i, j);
		status = EXIT_FAILURE;
		break;
	      }
	  if (status == EXIT_FAILURE)
	    break;
	}
      print_verification_msg (msg, status);
    }
  mpz_clears (p_zz, x_zz, NULL);

  for (int i = 0; i < n; i++)
    {
      mpz_clear (a_zz[i]);
      mpz_clear (b_zz[i]);
    }
  free (x_gmp);
  free (a_zz);
  free (b_zz);
}

/**************************************/

#endif
