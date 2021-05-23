#include "FiniteFields/small_prime_field_fft.h"

///////////////////////////////////////////////////////////
#define FOR_LOOP for 
///////////////////////////////////////////////////////////
#define MAX(x,y) (x>y?x:y)
///////////////////////////////////////////////////////////
#define MIN(x,y) (x<y?x:y)
///////////////////////////////////////////////////////////

void
init_montgomery_triple (montgomery_triple*dest, usfixn64 prime)
{
  /////////////////////////////////////
  dest->p = prime;
  dest->p_negative = -prime;

  mpz_t x_zz, p_zz, u64_zz;
  mpz_inits (x_zz, p_zz, u64_zz, NULL);
  mpz_set_u64 (p_zz, prime);
  mpz_set_u64 (u64_zz, (1UL << 63));
  mpz_mul_ui (u64_zz, u64_zz, 2);

  /////////////////////////////////////
  usfixn64 R2_mod_p;
  mpz_set (x_zz, u64_zz);
  mpz_powm_ui (x_zz, x_zz, 2, p_zz);
  R2_mod_p = mpz_get_u64 (x_zz);
  dest->R2_mod_p = R2_mod_p;

  /////////////////////////////////////

  mpz_set (x_zz, u64_zz);

  mpz_powm_ui (x_zz, x_zz, prime - 2, p_zz); //r_inv
  mpz_mul (x_zz, x_zz, u64_zz);
  mpz_sub_ui (x_zz, x_zz, 1);
  mpz_div (x_zz, x_zz, p_zz);
  dest->p_inv = mpz_get_u64 (x_zz);
  /////////////////////////////////////

#if VERBOSE
  printf ("[init_montgomery_triple]...\n[p=%llu, p_inv=%llu, R2_mod_p=%llu]\n",
      dest->p, dest->p_inv, dest->R2_mod_p);
  printf ("%s", shortline);
#endif

  mpz_clears (x_zz, p_zz, u64_zz, NULL);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

usfixn64 __attribute__ ((noinline))
mult_ab_mod_p (usfixn64 a, const usfixn64 b, const montgomery_triple * P)
{

  //////////////////////////////////////
  //// plain code for testing purposes.
  //////////////////////////////////////
  //    __int128 r = ((__int128)(a) * (__int128 )(b)) % ((__int128) (P->p));
  //     return (usfixn64)(r);

  /////////////////////////////////////
  //// (A+(Ap'%R).p)/R
  ////////////////////////////////////
  // A = a*b
  // A=a0+a1R;
  // m=(a0*p')%R
  // M=(m*P)=m0+m1.R
  // at the end, we need (A+M)/R
  // (a0+m0)/R + (a1+m1)

  ////////////////////////////////////
  // fast implementation via inline assembly
  usfixn64 m = 0;
  __asm__ volatile (
      "movq %1, %%rax\n\t"  //rax=%1;
      "mulq %2\n\t"// mult rax * %2
      "movq %%rax,%%rsi\n\t"//a0=mult.lo=rsi
      "movq %%rdx,%%rdi\n\t"//a1=mult.hi=rdi
      "mulq %3\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
      "mulq %4\n\t"//now mult by p itself.
      "addq %%rsi,%%rax\n\t"//a0+m0
      "adcq %%rdi,%%rdx\n\t"//a1+m1+carry(a0+m0)
      "subq %4,%%rdx\n\t"//m-=p;
      "movq %%rdx, %%rdi\n\t"//rdi=(m-p)
      "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "andq %4, %%rdi\n\t"// m+=(rdi);
      "addq %%rdi, %%rdx\n\t"
      "movq %%rdx,%0\n\t"
      :"=q" (m)
      :"q"(a),"q"(b),"q" (P->p_inv),"q"(P->p)
      :"%rax", "%rdx", "%rsi", "%rdi");

  //if (m>=P->p)
  //{
  //  printf("WARNING OVERFLOW!\n");
  //}

//  printf("m=%llu\n", m);
  return m;

  /////////////////////////////////////
  //// slow implementation via gcc 128-bit type
  /////////////////////////////////////
  //  __int128 A=__int128(a)*__int128(b);
  //  usfixn64 a0=  usfixn64 (A);
  //  A=A>>64;
  //  usfixn64 a1= usfixn64 (A);
  //
  //  usfixn64 m=a0*P->p_inv;
  //  __int128 M = __int128(m)*__int128(P->p);
  //
  //  usfixn64 m0=usfixn64(M);
  //  M=M>>64;
  //  usfixn64 m1 = usfixn64 (M);
  //
  //  usfixn64 s1=(a1)+(m1);
  //  usfixn64 s0=(m0+a0);
  //  if ((s0<m0)||(s0<a0))
  //    s1++;
  //
  ////  s1>>=64;
  //  m=usfixn64(s1);
  //
  //  if (m>=P->p)
  //    {
  //      m-=P->p;
  //    }
  //
  //  if (m>=P->p)
  //    {
  //      m-=P->p;
  //      printf("WARNING OVERFLOW!\n");
  //    }
  //  return m;
}

///////////////////////////////////////////////////////////

void __inline__
mult_ab_mod_p_ptr (usfixn64 * __restrict__ a, const usfixn64 *__restrict__ b,
		   const montgomery_triple * P)
{

  //////////////////////////////////////
  //// plain code for testing purposes.
  //////////////////////////////////////
//    __int128 r = ((__int128)(a) * (__int128 )(b)) % ((__int128) (P->p));
//
//     return (usfixn64)(r);

  /////////////////////////////////////
  //// (A+(Ap'%R).p)/R
  ////////////////////////////////////
  // A = a*b
  // A=a0+a1R;
  // m=(a0*p')%R
  // M=(m*P)=m0+m1.R
  // at the end, we need (A+M)/R
  // (a0+m0)/R + (a1+m1)

  ////////////////////////////////////
  // fast implementation via inline assembly
//  usfixn64 m = 0;
//  __asm__ __volatile__(
//      "movq %1, %%rax\n\t"  //rax=%1;
//      "mulq %2\n\t"// mult rax * %2
//      "movq %%rax,%%rsi\n\t"//a0=mult.lo=rsi
//      "movq %%rdx,%%rdi\n\t"//a1=mult.hi=rdi
//      "mulq %3\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
//      "mulq %4\n\t"//now mult by p itself.
//      "addq %%rsi,%%rax\n\t"//a0+m0
//      "adcq %%rdi,%%rdx\n\t"//a1+m1+carry(a0+m0)
//      "subq %4,%%rdx\n\t"//m-=p;
//      "movq %%rdx, %%rdi\n\t"//rdi=(m-p)
//      "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "andq %4, %%rdi\n\t"// m+=(rdi);
//      "addq %%rdi, %%rdx\n\t"
//      "movq %%rdx,%0\n\t"
//      :"=gm" (m)
//      :"gm"(*a),"q"(*b),"q" (P->p_inv),"q"(P->p)
//      :"%rax", "%rdx", "%rsi", "%rdi");
//
//  usfixn64 m0, m1;
//  __asm__ __volatile__(
//      "movq %2, %%rax\n\t"  //rax=%1;
//      "mulq %3\n\t"// mult rax * %2
//      "movq %%rax,%0\n\t"//a0=mult.lo=rsi
//      "movq %%rdx,%1\n\t"//a1=mult.hi=rdi
//      "mulq %4\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
//      "mulq %5\n\t"//now mult by p itself.
//      "addq %0,%%rax\n\t"//a0+m0
//      "adcq %1,%%rdx\n\t"//a1+m1+carry(a0+m0)
//      "subq %5,%%rdx\n\t"//m-=p;
//      "movq %%rdx, %1\n\t"//rdi=(m-p)
//      "sar $63, %1\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "andq %5, %1\n\t"// m+=(rdi);
//      "addq %1, %%rdx\n\t"
//      "movq %%rdx,%0\n\t"
//      :"+q" (m0),"+g" (m1)
//      :"g"(*a),"g"(*b),"g" (P->p_inv),"g"(P->p)
//      :"%rax", "%rdx");
//
//  //if (m>=P->p)
//  //{
//  //  printf("WARNING OVERFLOW!\n");
//  //}
//
////  printf("m=%llu\n", m);
//  *a = m0;

//  usfixn64 m0, m1;
  __asm__ __volatile__(
      "movq %0, %%rax\n\t"  //rax=%1;
      "mulq %1\n\t"// mult rax * %2
      "movq %%rax,%0\n\t"//a0=mult.lo=rsi
      "movq %%rdx,%%rsi\n\t"//a1=mult.hi=rdi
      "mulq %2\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
      "mulq %3\n\t"//now mult by p itself.
      "addq %0,%%rax\n\t"//a0+m0
      "adcq %%rsi,%%rdx\n\t"//a1+m1+carry(a0+m0)
      "subq %3,%%rdx\n\t"//m-=p;
      "movq %%rdx, %%rsi\n\t"//rdi=(m-p)
      "sar $63, %%rsi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "andq %3, %%rsi\n\t"// m+=(rdi);
      "addq %%rsi, %%rdx\n\t"
      "movq %%rdx,%0\n\t"
      :"+q" (*a)
      :"g"(*b),"g" (P->p_inv),"g"(P->p)
      :"%rax", "%rdx", "%rsi");

  //if (m>=P->p)
  //{
  //  printf("WARNING OVERFLOW!\n");
  //}

//  printf("m=%llu\n", m);
//  *a = m0;

  /////////////////////////////////////
  //// slow implementation via gcc 128-bit type
  /////////////////////////////////////
  //  __int128 A=__int128(a)*__int128(b);
  //  usfixn64 a0=  usfixn64 (A);
  //  A=A>>64;
  //  usfixn64 a1= usfixn64 (A);
  //
  //  usfixn64 m=a0*P->p_inv;
  //  __int128 M = __int128(m)*__int128(P->p);
  //
  //  usfixn64 m0=usfixn64(M);
  //  M=M>>64;
  //  usfixn64 m1 = usfixn64 (M);
  //
  //  usfixn64 s1=(a1)+(m1);
  //  usfixn64 s0=(m0+a0);
  //  if ((s0<m0)||(s0<a0))
  //    s1++;
  //
  ////  s1>>=64;
  //  m=usfixn64(s1);
  //
  //  if (m>=P->p)
  //    {
  //      m-=P->p;
  //    }
  //
  //  if (m>=P->p)
  //    {
  //      m-=P->p;
  //      printf("WARNING OVERFLOW!\n");
  //    }
  //  return m;
}

///////////////////////////////////////////////////////////

void __inline__
convertIn_GLOBAL_ptr (usfixn64 *a, const montgomery_triple* P)
{
  //reduce (a.(R^2%p))
//  return -1;
//  return a;
//  mult_ab_mod_p_ptr (a, &P->R2_mod_p, P);
//  printf("a=%llu, p=%llu, pinv=%llu, r2mod=%llu, m=%llu\n", a, P->p, P->p_inv, P->R2_mod_p, m);

//  usfixn64 m = 0;
//    __asm__ __volatile__(
//        "movq %1, %%rax\n\t"  //rax=%1;
//        "mulq %2\n\t"// mult rax * %2
//        "movq %%rax,%%rsi\n\t"//a0=mult.lo=rsi
//        "movq %%rdx,%%rdi\n\t"//a1=mult.hi=rdi
//        "mulq %3\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
//        "mulq %4\n\t"//now mult by p itself.
//        "addq %%rsi,%%rax\n\t"//a0+m0
//        "adcq %%rdi,%%rdx\n\t"//a1+m1+carry(a0+m0)
//        "subq %4,%%rdx\n\t"//m-=p;
//        "movq %%rdx, %%rdi\n\t"//rdi=(m-p)
//        "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//        "andq %4, %%rdi\n\t"// m+=(rdi);
//        "addq %%rdi, %%rdx\n\t"
//        "movq %%rdx,%0\n\t"
//        :"=g" (m)
//        :"g"(*a),"g"(P->R2_mod_p),"g" (P->p_inv),"g"(P->p)
//        :"%rax", "%rdx", "%rsi", "%rdi");

  usfixn64 m0, m1;
  __asm__ __volatile__(
      "movq %2, %%rax\n\t"  //rax=%1;
      "mulq %3\n\t"// mult rax * %2
      "movq %%rax,%0\n\t"//a0=mult.lo=rsi
      "movq %%rdx,%1\n\t"//a1=mult.hi=rdi
      "mulq %4\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
      "mulq %5\n\t"//now mult by p itself.
      "addq %0,%%rax\n\t"//a0+m0
      "adcq %1,%%rdx\n\t"//a1+m1+carry(a0+m0)
      "subq %5,%%rdx\n\t"//m-=p;
      "movq %%rdx, %1\n\t"//rdi=(m-p)
      "sar $63, %1\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "andq %5, %1\n\t"// m+=(rdi);
      "addq %1, %%rdx\n\t"
      "movq %%rdx,%0\n\t"
      :"+q" (m0),"+g" (m1)
      :"g"(*a),"g"(P->R2_mod_p),"g" (P->p_inv),"g"(P->p)
      :"%rax", "%rdx");

  //if (m>=P->p)
  //{
  //  printf("WARNING OVERFLOW!\n");
  //}

  //  printf("m=%llu\n", m);
  *a = m0;

}

///////////////////////////////////////////////////////////

void __inline__
convertOut_GLOBAL_ptr (usfixn64 * r, const montgomery_triple* P)
{
  // (a+(ap'%R)p)/R
//  usfixn64 m = (a*P->p_inv);
//    __int128 m128=__int128(m)*__int128(P->p);
//    m128+=(__int128(a));
//    m128=m128>>64;
//    return usfixn64 (m128);
//  return -1;
//  usfixn64 one = 1;
//  mult_ab_mod_p_ptr (r, &one, P);

  usfixn64 m0, m1;
  __asm__ __volatile__(
      "movq $1, %%rax\n\t"  //rax=%1;
      "mulq %2\n\t"// mult rax * %2
      "movq %%rax,%0\n\t"//a0=mult.lo=rsi
      "movq %%rdx,%1\n\t"//a1=mult.hi=rdi
      "mulq %3\n\t"//mult (%2(pinv)*rax(a0)=a0*p_inv =  only need the lower part.
      "mulq %4\n\t"//now mult by p itself.
      "addq %0,%%rax\n\t"//a0+m0
      "adcq %1,%%rdx\n\t"//a1+m1+carry(a0+m0)
      "subq %4,%%rdx\n\t"//m-=p;
      "movq %%rdx, %1\n\t"//rdi=(m-p)
      "sar $63, %1\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "andq %4, %1\n\t"// m+=(rdi);
      "addq %1, %%rdx\n\t"
      "movq %%rdx,%0\n\t"
      :"+q" (m0),"+g" (m1)
      :"g"(*r),"g" (P->p_inv),"g"(P->p)
      :"%rax", "%rdx");

  //if (m>=P->p)
  //{
  //  printf("WARNING OVERFLOW!\n");
  //}

  //  printf("m=%llu\n", m);
  *r = m0;

}

///////////////////////////////////////////////////////////

usfixn64
exponent_mod_p_SPF (usfixn64 a, usfixn64 n, const montgomery_triple * P)
{
  int n_log2 = log2 (n);
  usfixn64 x;
  if ((1 << n_log2) != n)
    {
      x = 1;
      convertIn_GLOBAL_ptr (&x, P);
      for (int i = 0; i < n; i++)
	{
	  x = mult_ab_mod_p (x, a, P);
	}
    }
  else // n is a power of 2!
    {
      x = a;
      for (int i = 0; i < n_log2; i++)
	x = mult_ab_mod_p (x, x, P);
    }

  return x;
}

///////////////////////////////////////////////////////////

void __inline__
DFT_2_double(usfixn64* __restrict__ a0, usfixn64* __restrict__ a1,
		      usfixn64* __restrict__ a2, usfixn64* __restrict__ a3,
		      const montgomery_triple* P)
{
//  t[0] = *a0;
//  t[1] = *a1;
//  t[2] = *a2;
//  t[3] = *a3;

//  s[0] = t[0] + t[1];
//  s[1] = t[0] + P->p - t[1];
//  s[2] = t[2] + t[3];
//  s[3] = t[2] + P->p - t[3];

  ////////////////////////////////////////////////
//  usfixn64 s[4];
//  s[0] = *a0 + *a1;
//  s[1] = *a0 + P->p - *a1;
//  s[2] = *a2 + *a3;
//  s[3] = *a2 + P->p - *a3;
//
//  s[0] -= P->p;
//  s[1] -= P->p;
//  s[2] -= P->p;
//  s[3] -= P->p;
//  s[0] += (-(s[0] >> 63)) & P->p;
//  s[1] += (-(s[1] >> 63)) & P->p;
//  s[2] += (-(s[2] >> 63)) & P->p;
//  s[3] += (-(s[3] >> 63)) & P->p;
//
//  *a0 = s[0];
//  *a1 = s[1];
//  *a2 = s[2];
//  *a3 = s[3];
  ////////////////////////////////////////////////

    usfixn64 t1, t3;
    t1 = *a0 - *a1;
    t3 = *a2 - *a3;

    *a0 += *a1;
    *a2 += *a3;

    *a1=t1;
    *a3=t3;

    *a0+=(P->p_negative);
    *a2+=(P->p_negative);

    *a0 += (-(*a0 >> 63)) & P->p;
    *a1 += (-(*a1 >> 63)) & P->p;
    *a2 += (-(*a2 >> 63)) & P->p;
    *a3 += (-(*a3 >> 63)) & P->p;

//    __asm__ __volatile__ (
//    "movq %0, %%rsp\n\t" //rax=%1;
//    "sar $63, %%rsp\n\t" //rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//    "andq %1, %%rsp\n\t" // m+=(rdi);
//    "addq %%rsp, %0\n\t"
//    :"+g" (*a0):"g"(P->p):"%rsp");
//
//    __asm__ __volatile__ (
//    "movq %0, %%rsp\n\t" //rax=%1;
//    "sar $63, %%rsp\n\t" //rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//    "andq %1, %%rsp\n\t" // m+=(rdi);
//    "addq %%rsp, %0\n\t"
//    :"+g" (*a1):"g"(P->p) :"%rsp");
//
//    __asm__ __volatile__ (
//    "movq %0, %%rsp\n\t" //rax=%1;
//    "sar $63, %%rsp\n\t" //rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//    "andq %1, %%rsp\n\t" // m+=(rdi);
//    "addq %%rsp, %0\n\t"
//    :"+g" (*a2):"g"(P->p):"%rsp");
//
//    __asm__ __volatile__ (
//    "movq %0, %%rsp\n\t" //rax=%1;
//    "sar $63, %%rsp\n\t" //rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//    "andq %1, %%rsp\n\t" // m+=(rdi);
//    "addq %%rsp, %0\n\t"
//    :"+g" (*a3):"g"(P->p):"%rsp");

  ////////////////////////////////////////////////
//  usfixn64 t1, t3;
//  t1 = *a0 + P->p - *a1;
//  t3 = *a2 + P->p - *a3;
//
//  *a0 += *a1;
//  *a2 += *a3;
//
//  *a1=t1;
//  *a3=t3;
//
//  ////////////////////////////////////////////////
//  __asm__ __volatile__ (
//      "movq %4, %%rsi\n\t"  //rax=%1;
//      "movq %4, %%rdi\n\t"//rax=%1;
//
//      "neg %%rsi\n\t"//rax=%1;
//
//      "addq %%rsi, %0\n\t"//rax=%1;
//      "addq %%rsi, %1\n\t"//rax=%1;
//      "addq %%rsi, %2\n\t"//rax=%1;
//      "addq %%rsi, %3\n\t"//rax=%1;
//
//      "movq %0, %%r8\n\t"//rax=%1;
//      "movq %1, %%r9\n\t"//rax=%1;
//      "movq %2, %%r10\n\t"//rax=%1;
//      "movq %3, %%r11\n\t"//rax=%1;
//
//      "sar $63, %%r8\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "sar $63, %%r9\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "sar $63, %%r10\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "sar $63, %%r11\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//
//      "andq %%rdi, %%r8\n\t"// m+=(rdi);
//      "andq %%rdi, %%r9\n\t"// m+=(rdi);
//      "andq %%rdi, %%r10\n\t"// m+=(rdi);
//      "andq %%rdi, %%r11\n\t"// m+=(rdi);
//
//      "addq %%r8, %0\n\t"
//      "addq %%r9, %1\n\t"
//      "addq %%r10, %2\n\t"
//      "addq %%r11, %3\n\t"
//
//      :"+g" (*a0),"+g" (*a1),"+g" (*a2),"+g" (*a3)
//      :"g"(P->p)
//      //:"q"(*a),"q"(*b),"q" (*p_inv),"q"(*p)
//      :"%rsi", "%rdi", "%r8", "%r9", "%r10", "%r11");
  ////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////

void 
DFT_2_x4(usfixn64* __restrict__ a0, usfixn64* __restrict__ a1,
		      usfixn64* __restrict__ a2, usfixn64* __restrict__ a3,
		      usfixn64* __restrict__ a4, usfixn64* __restrict__ a5,
		      usfixn64* __restrict__ a6, usfixn64* __restrict__ a7,
		      const montgomery_triple* P)
{

//  t[0] = *a0;
//  t[1] = *a1;
//  t[2] = *a2;
//  t[3] = *a3;

//  s[0] = t[0] + t[1];
//  s[1] = t[0] + P->p - t[1];
//  s[2] = t[2] + t[3];
//  s[3] = t[2] + P->p - t[3];

  ////////////////////////////////////////////////
//  s[0] = *a0 + *a1;
//  s[1] = *a0 + P->p - *a1;
//  s[2] = *a2 + *a3;
//  s[3] = *a2 + P->p - *a3;

//  s[0] -= P->p;
//  s[1] -= P->p;
//  s[2] -= P->p;
//  s[3] -= P->p;
//
//  const usfixn64 u63=(1UL<<63);
////  const usfixn64 u63_mask=(1<<63)-1;
////  const usfixn64 u64_mask_and_p= U64_MASK & & P->p;
////  s[0] += (-(s[0] >> 63)) & P->p;
//  s[0] += (-(s[0] >> 63)) & P->p;
//  s[1] += (-(s[1] >> 63)) & P->p;
//  s[2] += (-(s[2] >> 63)) & P->p;
//  s[3] += (-(s[3] >> 63)) & P->p;

//    *a0 = s[0];
//    *a1 = s[1];
//    *a2 = s[2];
//    *a3 = s[3];
  ////////////////////////////////////////////////
  usfixn64 t1, t3, t5, t7;
  t1 = *a0 + P->p - *a1;
  t3 = *a2 + P->p - *a3;
  t5 = *a4 + P->p - *a5;
  t7 = *a6 + P->p - *a7;

  *a0 += *a1;
  *a2 += *a3;
  *a4 += *a5;
  *a6 += *a7;

  *a1=t1;
  *a3=t3;
  *a5=t5;
  *a7=t7;

  ////////////////////////////////////////////////
  __asm__ __volatile__ (
      "movq %8, %%rsi\n\t"  //rax=%1;
      "movq %8, %%rdi\n\t"//rax=%1;

      "neg %%rsi\n\t"//rax=%1;

      "addq %%rsi, %0\n\t"//rax=%1;
      "addq %%rsi, %1\n\t"//rax=%1;
      "addq %%rsi, %2\n\t"//rax=%1;
      "addq %%rsi, %3\n\t"//rax=%1;
      "addq %%rsi, %4\n\t"//rax=%1;
      "addq %%rsi, %5\n\t"//rax=%1;
      "addq %%rsi, %6\n\t"//rax=%1;
      "addq %%rsi, %7\n\t"//rax=%1;

      "movq %0, %%r8\n\t"//rax=%1;
      "movq %1, %%r9\n\t"//rax=%1;
      "movq %2, %%r10\n\t"//rax=%1;
      "movq %3, %%r11\n\t"//rax=%1;
      "movq %4, %%r12\n\t"//rax=%1;
      "movq %5, %%r13\n\t"//rax=%1;
      "movq %6, %%r14\n\t"//rax=%1;
      "movq %7, %%r15\n\t"//rax=%1;

      "sar $63, %%r8\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r9\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r10\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r11\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r12\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r13\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r14\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
      "sar $63, %%r15\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)

      "andq %%rdi, %%r8\n\t"// m+=(rdi);
      "andq %%rdi, %%r9\n\t"// m+=(rdi);
      "andq %%rdi, %%r10\n\t"// m+=(rdi);
      "andq %%rdi, %%r11\n\t"// m+=(rdi);
      "andq %%rdi, %%r12\n\t"// m+=(rdi);
      "andq %%rdi, %%r13\n\t"// m+=(rdi);
      "andq %%rdi, %%r14\n\t"// m+=(rdi);
      "andq %%rdi, %%r15\n\t"// m+=(rdi);

      "addq %%r8, %0\n\t"
      "addq %%r9, %1\n\t"
      "addq %%r10, %2\n\t"
      "addq %%r11, %3\n\t"
      "addq %%r12, %4\n\t"
      "addq %%r13, %5\n\t"
      "addq %%r14, %6\n\t"
      "addq %%r15, %7\n\t"

      :"+g" (*a0),"+g" (*a1),"+g" (*a2),"+g" (*a3),
       "+g" (*a4),"+g" (*a5),"+g" (*a6),"+g" (*a7)
      :"g"(P->p)
      //:"q"(*a),"q"(*b),"q" (*p_inv),"q"(*p)
      :"%rsi", "%rdi", "%r8", "%r9", "%r10", "%r11", "%r12", "%r13", "%r14", "%r15");
  ////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////

void __inline__
DFT_2_double_working(usfixn64* __restrict__ a0, usfixn64* __restrict__ a1,
	      usfixn64* __restrict__ a2, usfixn64* __restrict__ a3,
	      const montgomery_triple* P)
{
//
//  t[0] = *a0;
//  t[1] = *a1;
//  t[2] = *a2;
//  t[3] = *a3;

//  s[0] = t[0] + t[1];
//  s[1] = t[0] + P->p - t[1];
//  s[2] = t[2] + t[3];
//  s[3] = t[2] + P->p - t[3];

  usfixn64 t0, t1, t2, t3;

  t1 = *a0 + P->p - *a1;
  t3 = *a2 + P->p - *a3;
//  t2 = *a2 + *a3;
//  t0 = *a0 + *a1;

//  *a0 = t0;
//  *a2 = t2;
  *a0+=*a1;
  *a2+=*a3;

  *a1 = t1;
  *a3 = t3;

  *a0 -= P->p;
  *a1 -= P->p;
  *a2 -= P->p;
  *a3 -= P->p;
//
//  *a0 += (-(*a0 >> 63)) & P->p;
//  *a1 += (-(*a1 >> 63)) & P->p;
//  *a2 += (-(*a2 >> 63)) & P->p;
//  *a3 += (-(*a3 >> 63)) & P->p;


  usfixn64 c0,c1,c2,c3;
  c0=-(*a0 >>63);
  c1=-(*a1 >>63);
  c2=-(*a2 >>63);
  c3=-(*a3 >>63);
  c0&=P->p;
  c1&=P->p;
  c2&=P->p;
  c3&=P->p;

  *a0+=c0;
  *a1+=c1;
  *a2+=c2;
  *a3+=c3;

  //  *a0 += (-(*a0 >> 63)) & P->p;
  //  *a1 += (-(*a1 >> 63)) & P->p;
  //  *a2 += (-(*a2 >> 63)) & P->p;
  //  *a3 += (-(*a3 >> 63)) & P->p;


  //////////////////////////////
//  __asm__ __volatile__ (
//      "movq %0, %%rsi\n\t"  //rax=%1;
//      "subq %4, %%rsi\n\t"//rax=%1;
//      "movq %%rsi, %%rdi\n\t"//rax=%1;
//      "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "andq %4, %%rdi\n\t"// m+=(rdi);
//      "addq %%rdi, %%rsi\n\t"
//      "movq %%rsi, %0\n\t"
//
//      "movq %1, %%rsi\n\t"//rax=%1;
//      "subq %4, %%rsi\n\t"//rax=%1;
//      "movq %%rsi, %%rdi\n\t"//rax=%1;
//      "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "andq %4, %%rdi\n\t"// m+=(rdi);
//      "addq %%rdi, %%rsi\n\t"
//      "movq %%rsi, %1\n\t"
//
//      "movq %2, %%rsi\n\t"//rax=%1;
//      "subq %4, %%rsi\n\t"//rax=%1;
//      "movq %%rsi, %%rdi\n\t"//rax=%1;
//      "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "andq %4, %%rdi\n\t"// m+=(rdi);
//      "addq %%rdi, %%rsi\n\t"
//      "movq %%rsi, %2\n\t"
//
//      "movq %3, %%rsi\n\t"//rax=%1;
//      "subq %4, %%rsi\n\t"//rax=%1;
//      "movq %%rsi, %%rdi\n\t"//rax=%1;
//      "sar $63, %%rdi\n\t"//rdi=(m-p)>>63 (arithmetic shift for extending the MSB to the right)
//      "andq %4, %%rdi\n\t"// m+=(rdi);
//      "addq %%rdi, %%rsi\n\t"
//      "movq %%rsi, %3\n\t"
//      :"+g" (s[0]),"+g" (s[1]),"+g" (s[2]),"+g" (s[3])
//      :"g"(P->p)
//      //:"q"(*a),"q"(*b),"q" (*p_inv),"q"(*p)
//      :"%rsi", "%rdi");
//
//  *a0 = s[0];
//  *a1 = s[1];
//  *a2 = s[2];
//  *a3 = s[3];
}

///////////////////////////////////////////////////////////

//vectorized_avx2
void
DFT_2_double_avx2 (usfixn64* a0, usfixn64* a1, usfixn64* a2, usfixn64* a3,
		   const montgomery_triple* P)
{

//  usfixn64 __attribute__((aligned(32))) s[4];
//  s[0] = *a0 + *a1 - P->p;
//  s[1] = *a0 - *a1;
//  s[2] = *a2 + *a3 - P->p;
//  s[3] = *a2 - *a3;
//
//  __m256i zero_vec, s_vec, sub_vec, p_vec;
////  zero_vec = _mm256_setzero_si256 (); //[0,0,0,0]
//  zero_vec = _mm256_xor_si256 (zero_vec, zero_vec);
////  s_vec = _mm256_setr_epi64x (s[0], s[1], s[2], s[3]); //[s0,s1,s2,s3]
//  s_vec = _mm256_load_si256 ((__m256i *) &s); //[0], s[1], s[2], s[3]); //[s0,s1,s2,s3]
////  s_vec = _mm256_setr_epi64x (s[0], s[1], s[2], s[3]); //[s0,s1,s2,s3]
//  p_vec = _mm256_set1_epi64x (P->p); //[p,p,p,p]
////  sub_vec = s_vec;
//
////  s[0] += (s[0] >> 63) & P->p;
////  s[1] += (s[1] >> 63) & P->p;
////  s[2] += (s[2] >> 63) & P->p;
////  s[3] += (s[3] >> 63) & P->p;
//  sub_vec = _mm256_srli_epi64 (s_vec, 63); //[s0>>63,s1>>63,s2>>63,s3>>63]
//  sub_vec = _mm256_sub_epi64 (zero_vec, sub_vec); //[-s0,-s1,-s2,-s3]
//  sub_vec = _mm256_and_si256 (sub_vec, p_vec); //[(-s0)&p,(-s1)&p,(-s2)&p,(-s3)&p]
//  s_vec = _mm256_add_epi64 (s_vec, sub_vec); // addition
//
//  _mm256_store_si256 ((__m256i *) &s, s_vec);
////  usfixn64 *s_ptr = (usfixn64*) &s_vec;
//  *a0 = s[0];
//  *a1 = s[1];
//  *a2 = s[2];
//  *a3 = s[3];
}

///////////////////////////////////////////////////////////

void __inline__
swap (usfixn64* a, usfixn64* b)
{

//  __asm__ volatile (
//      "xchg %0,%1\n\n":
//      "+g"(*a),"+g"(*b):
//      :);

  __asm__ __volatile__ (
      "xchg %0,%1\n\n":
      "+g"(*a),"+g"(*b):
      :);

  //////////////////////////////////
  //// naive
  //////////////////////////////////
//  usfixn64 tmp;
//  tmp = *a;
//  *a = *b;
//  *b = tmp;
  //////////////////////////////////
}

///////////////////////////////////////////////////////////

//swap pointers x and y,
//switch the value of check-bit from 0 to 1 and vice versa.
void
swap_ptr (usfixn64**x, usfixn64**y, int * check_bit)
{
//  if (*check_bit==0)
//    *check_bit=1;
//  else
//    *check_bit=0;
  (*check_bit) += 1;
  (*check_bit) &= 0x1;

//  printf("check-bit after change = %d\n", *check_bit);
  usfixn64 *tmp = *x;
  *x = *y;
  *y = tmp;
}

///////////////////////////////////////////////////////////

//pow_omega is already allocated and is of size n*sizeof(usfixn64);
void
precompute_pow_omega (usfixn64* pow_omega, const usfixn64 omega, int n,
		      const montgomery_triple * P)
{
//  printf("in precompute pow omega !\n");

  pow_omega[0] = 1;
  convertIn_GLOBAL_ptr (&pow_omega[0], P);
  pow_omega[1] = omega; //we assume omega is already in Montgomery representation.
  for (int i = 2; i < n; i++)
    {
      pow_omega[i] = mult_ab_mod_p (pow_omega[i - 1], omega, P);
    }
}

///////////////////////////////////////////////////////////

void
DFT_4 (int n, usfixn64* x, const usfixn64 *base_case_omega_vec,
       const montgomery_triple* P)
{
//#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n; i++)
    {
      usfixn64 *a = &x[i * 4];

      // dft on permutated indices
      DFT_2_double (&a[0], &a[2], &a[1], &a[3], P);

      //twiddle
      mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[1], P);

      // dft on permutated indices
      DFT_2_double (&a[0], &a[1], &a[2], &a[3], P);

      // final permutation
      swap (&a[1], &a[2]);
    }
}

///////////////////////////////////////////////////////////

void
DFT_8 (int n, usfixn64* x, const usfixn64 *base_case_omega_vec,
       const montgomery_triple* P)
{
//#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n;
      i++
    )
      {
	usfixn64 *a = &x[i * 8];

	//dft on permutated indices
	DFT_2_double (&a[0], &a[4], &a[2], &a[6], P);
	DFT_2_double (&a[1], &a[5], &a[3], &a[7], P);

	//twiddle
	mult_ab_mod_p_ptr (&a[6], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[2], P);

	//dft on permutated indices
	DFT_2_double (&a[0], &a[2], &a[1], &a[3], P);
	DFT_2_double (&a[4], &a[6], &a[5], &a[7], P);

	//twiddle
	mult_ab_mod_p_ptr (&a[5], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[3], P);

	//dft on permutated indices
	DFT_2_double (&a[0], &a[1], &a[2], &a[3], P);
	DFT_2_double (&a[4], &a[5], &a[6], &a[7], P);

	//final permutation
	swap (&a[1], &a[4]);
	swap (&a[3], &a[6]);
      }
  }

///////////////////////////////////////////////////////////

void
DFT_8_x4(int n, usfixn64* x, const usfixn64 *base_case_omega_vec,
       const montgomery_triple* P)
{
//#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  FOR_LOOP (int i = 0; i < n; i++ )
      {
	usfixn64 *a = &x[i * 8];

	//dft on permutated indices
	DFT_2_x4  (&a[0], &a[4], &a[2], &a[6], &a[1], &a[5], &a[3], &a[7], P);

	//twiddle
	mult_ab_mod_p_ptr (&a[6], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[2], P);

	//dft on permutated indices
	DFT_2_x4  (&a[0], &a[2], &a[1], &a[3], &a[4], &a[6], &a[5], &a[7], P);

	//twiddle
	mult_ab_mod_p_ptr (&a[5], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[3], P);

	//dft on permutated indices
	DFT_2_x4  (&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], P);

	//final permutation
	swap (&a[1], &a[4]);
	swap (&a[3], &a[6]);
      }
  }

///////////////////////////////////////////////////////////

void
DFT_16 (int n, usfixn64* x, const usfixn64 * base_case_omega_vec,
	const montgomery_triple* P)
{
//#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n;
      i++
    )
      {
	usfixn64 *a = &x[i * 16];

	//dft on permutated indices
	DFT_2_double (&a[0], &a[8], &a[1], &a[9], P);
	DFT_2_double (&a[2], &a[10], &a[3], &a[11], P);
	DFT_2_double (&a[4], &a[12], &a[5], &a[13], P);
	DFT_2_double (&a[6], &a[14], &a[7], &a[15], P);

	//twiddle
	mult_ab_mod_p_ptr (&a[12], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[4], P);

	//dft on permutated indices
	DFT_2_double (&a[0], &a[4], &a[1], &a[5], P);
	DFT_2_double (&a[2], &a[6], &a[3], &a[7], P);
	DFT_2_double (&a[8], &a[12], &a[9], &a[13], P);
	DFT_2_double (&a[10], &a[14], &a[11], &a[15], P);

	//twiddle
	mult_ab_mod_p_ptr (&a[6], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[10], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[6], P);

	//dft on permutated indices.
	DFT_2_double (&a[0], &a[2], &a[1], &a[3], P);
	DFT_2_double (&a[4], &a[6], &a[5], &a[7], P);
	DFT_2_double (&a[8], &a[10], &a[9], &a[11], P);
	DFT_2_double (&a[12], &a[14], &a[13], &a[15], P);

	mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[5], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[9], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[5], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[3], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[7], P);

	DFT_2_double (&a[0], &a[1], &a[2], &a[3], P);
	DFT_2_double (&a[4], &a[5], &a[6], &a[7], P);
	DFT_2_double (&a[8], &a[9], &a[10], &a[11], P);
	DFT_2_double (&a[12], &a[13], &a[14], &a[15], P);

	//final permutation.
	swap (&a[1], &a[8]);
	swap (&a[2], &a[4]);
	swap (&a[3], &a[12]);
	swap (&a[5], &a[10]);
	swap (&a[7], &a[14]);
	swap (&a[11], &a[13]);
      }
  }

///////////////////////////////////////////////////////////
//DFT_32 working with DFT_2_double as the base-case.
//VERIFIED.
void
DFT_32(int n, usfixn64* __restrict__ x,
	const usfixn64 * __restrict__ base_case_omega_vec,
	const montgomery_triple* P)
{

//#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n;
      i++
    )
      {
	usfixn64 *a = &x[i * 32];
	//  usfixn64 p=P->p, p_inv=P->p_inv;

	DFT_2_double (&a[0], &a[16], &a[1], &a[17], P);
	DFT_2_double (&a[2], &a[18], &a[3], &a[19], P);
	DFT_2_double (&a[4], &a[20], &a[5], &a[21], P);
	DFT_2_double (&a[6], &a[22], &a[7], &a[23], P);
	DFT_2_double (&a[8], &a[24], &a[9], &a[25], P);
	DFT_2_double (&a[10], &a[26], &a[11], &a[27], P);
	DFT_2_double (&a[12], &a[28], &a[13], &a[29], P);
	DFT_2_double (&a[14], &a[30], &a[15], &a[31], P);

	mult_ab_mod_p_ptr (&a[24], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[25], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[26], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[28], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[8], P);

	DFT_2_double (&a[0], &a[8], &a[1], &a[9], P);
	DFT_2_double (&a[2], &a[10], &a[3], &a[11], P);
	DFT_2_double (&a[4], &a[12], &a[5], &a[13], P);
	DFT_2_double (&a[6], &a[14], &a[7], &a[15], P);
	DFT_2_double (&a[16], &a[24], &a[17], &a[25], P);
	DFT_2_double (&a[18], &a[26], &a[19], &a[27], P);
	DFT_2_double (&a[20], &a[28], &a[21], &a[29], P);
	DFT_2_double (&a[22], &a[30], &a[23], &a[31], P);

//      for (int i = 20; i < 24; i++)
	mult_ab_mod_p_ptr (&a[20], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[21], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[22], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[4], P);

//      for (int i = 12; i < 16; i++)
	mult_ab_mod_p_ptr (&a[12], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[8], P);

//      for (int i = 28; i < 32; i++)
	mult_ab_mod_p_ptr (&a[28], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[12], P);

	DFT_2_double (&a[0], &a[4], &a[16], &a[20], P);
	DFT_2_double (&a[8], &a[12], &a[24], &a[28], P);
	DFT_2_double (&a[2], &a[6], &a[18], &a[22], P);
	DFT_2_double (&a[10], &a[14], &a[26], &a[30], P);
	DFT_2_double (&a[1], &a[5], &a[17], &a[21], P);
	DFT_2_double (&a[9], &a[13], &a[25], &a[29], P);
	DFT_2_double (&a[3], &a[7], &a[19], &a[23], P);
	DFT_2_double (&a[11], &a[15], &a[27], &a[31], P);

	mult_ab_mod_p_ptr (&a[18], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[19], &base_case_omega_vec[2], P);

	mult_ab_mod_p_ptr (&a[10], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[4], P);

	mult_ab_mod_p_ptr (&a[26], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[6], P);

	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[6], &base_case_omega_vec[8], P);

	mult_ab_mod_p_ptr (&a[22], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[10], P);

	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[12], P);

	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[14], P);

	DFT_2_double (&a[0], &a[2], &a[1], &a[3], P);
	DFT_2_double (&a[4], &a[6], &a[5], &a[7], P);
	DFT_2_double (&a[8], &a[10], &a[9], &a[11], P);
	DFT_2_double (&a[12], &a[14], &a[13], &a[15], P);
	DFT_2_double (&a[16], &a[18], &a[17], &a[19], P);
	DFT_2_double (&a[20], &a[22], &a[21], &a[23], P);
	DFT_2_double (&a[24], &a[26], &a[25], &a[27], P);
	DFT_2_double (&a[28], &a[30], &a[29], &a[31], P);

	mult_ab_mod_p_ptr (&a[17], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr (&a[9], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[25], &base_case_omega_vec[3], P);
	mult_ab_mod_p_ptr (&a[5], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[21], &base_case_omega_vec[5], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[7], P);
	mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[19], &base_case_omega_vec[9], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[11], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[13], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[15], P);

	DFT_2_double (&a[0], &a[1], &a[2], &a[3], P);
	DFT_2_double (&a[4], &a[5], &a[6], &a[7], P);
	DFT_2_double (&a[8], &a[9], &a[10], &a[11], P);
	DFT_2_double (&a[12], &a[13], &a[14], &a[15], P);
	DFT_2_double (&a[16], &a[17], &a[18], &a[19], P);
	DFT_2_double (&a[20], &a[21], &a[22], &a[23], P);
	DFT_2_double (&a[24], &a[25], &a[26], &a[27], P);
	DFT_2_double (&a[28], &a[29], &a[30], &a[31], P);

	swap (&a[1], &a[16]);
	swap (&a[2], &a[8]);
	swap (&a[3], &a[24]);
	swap (&a[5], &a[20]);
	swap (&a[6], &a[12]);
	swap (&a[7], &a[28]);
	swap (&a[9], &a[18]);
	swap (&a[11], &a[26]);
	swap (&a[13], &a[22]);
	swap (&a[15], &a[30]);
	swap (&a[19], &a[25]);
	swap (&a[23], &a[29]);
      }
  }

///////////////////////////////////////////////////////////

void
DFT_32_x4(int n, usfixn64* __restrict__ x,
	const usfixn64 * __restrict__ base_case_omega_vec,
	const montgomery_triple* P)
{

//#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  FOR_LOOP (int i = 0; i < n; i++)
      {
	usfixn64 *a = &x[i * 32];
	//  usfixn64 p=P->p, p_inv=P->p_inv;

	DFT_2_x4 (&a[0], &a[16], &a[1], &a[17], &a[2], &a[18], &a[3], &a[19], P);
	DFT_2_x4(&a[4], &a[20], &a[5], &a[21], &a[6], &a[22], &a[7], &a[23], P);
	DFT_2_x4(&a[8], &a[24], &a[9], &a[25], &a[10], &a[26], &a[11], &a[27], P);
	DFT_2_x4(&a[12], &a[28], &a[13], &a[29], &a[14], &a[30], &a[15], &a[31], P);

	mult_ab_mod_p_ptr (&a[24], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[25], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[26], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[28], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[8], P);

	DFT_2_x4 (&a[0], &a[8], &a[1], &a[9], &a[2], &a[10], &a[3], &a[11], P);
	DFT_2_x4 (&a[4], &a[12], &a[5], &a[13], &a[6], &a[14], &a[7], &a[15], P);
	DFT_2_x4 (&a[16], &a[24], &a[17], &a[25], &a[18], &a[26], &a[19], &a[27], P);
	DFT_2_x4 (&a[20], &a[28], &a[21], &a[29], &a[22], &a[30], &a[23], &a[31], P);

//      for (int i = 20; i < 24; i++)
	mult_ab_mod_p_ptr (&a[20], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[21], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[22], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[4], P);

//      for (int i = 12; i < 16; i++)
	mult_ab_mod_p_ptr (&a[12], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[8], P);

//      for (int i = 28; i < 32; i++)
	mult_ab_mod_p_ptr (&a[28], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[12], P);

	DFT_2_x4 (&a[0], &a[4], &a[16], &a[20], &a[8], &a[12], &a[24], &a[28], P);
	DFT_2_x4 (&a[2], &a[6], &a[18], &a[22], &a[10], &a[14], &a[26], &a[30], P);
	DFT_2_x4 (&a[1], &a[5], &a[17], &a[21], &a[9], &a[13], &a[25], &a[29], P);
	DFT_2_x4 (&a[3], &a[7], &a[19], &a[23], &a[11], &a[15], &a[27], &a[31], P);

	mult_ab_mod_p_ptr (&a[18], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[19], &base_case_omega_vec[2], P);

	mult_ab_mod_p_ptr (&a[10], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[4], P);

	mult_ab_mod_p_ptr (&a[26], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[6], P);

	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[6], &base_case_omega_vec[8], P);

	mult_ab_mod_p_ptr (&a[22], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[10], P);

	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[12], P);

	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[14], P);

	DFT_2_x4  (&a[0], &a[2], &a[1], &a[3], &a[4], &a[6], &a[5], &a[7], P);
	DFT_2_x4  (&a[8], &a[10], &a[9], &a[11], &a[12], &a[14], &a[13], &a[15], P);
	DFT_2_x4  (&a[16], &a[18], &a[17], &a[19], &a[20], &a[22], &a[21], &a[23], P);
	DFT_2_x4  (&a[24], &a[26], &a[25], &a[27], &a[28], &a[30], &a[29], &a[31], P);

	mult_ab_mod_p_ptr (&a[17], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr (&a[9], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[25], &base_case_omega_vec[3], P);
	mult_ab_mod_p_ptr (&a[5], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[21], &base_case_omega_vec[5], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[7], P);
	mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[19], &base_case_omega_vec[9], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[11], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[13], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[15], P);

	DFT_2_x4  (&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], P);
	DFT_2_x4  (&a[8], &a[9], &a[10], &a[11], &a[12], &a[13], &a[14], &a[15], P);
	DFT_2_x4  (&a[16], &a[17], &a[18], &a[19], &a[20], &a[21], &a[22], &a[23], P);
	DFT_2_x4  (&a[24], &a[25], &a[26], &a[27], &a[28], &a[29], &a[30], &a[31], P);

	swap (&a[1], &a[16]);
	swap (&a[2], &a[8]);
	swap (&a[3], &a[24]);
	swap (&a[5], &a[20]);
	swap (&a[6], &a[12]);
	swap (&a[7], &a[28]);
	swap (&a[9], &a[18]);
	swap (&a[11], &a[26]);
	swap (&a[13], &a[22]);
	swap (&a[15], &a[30]);
	swap (&a[19], &a[25]);
	swap (&a[23], &a[29]);
      }
  }

///////////////////////////////////////////////////////////

void
DFT_64 (int n, usfixn64* __restrict__ x,
	const usfixn64 *__restrict__ base_case_omega_vec,
	const montgomery_triple* P)
{
  //#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
//  __cilkrts_set_param	("nworkers","8");
//  printf("nworkers=%d\n", __cilkrts_get_total_workers());

#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n;
      i++
    )
      {
	usfixn64 *a = &x[i * 64];
	DFT_2_double (&a[0], &a[32], &a[1], &a[33], P);
	DFT_2_double (&a[2], &a[34], &a[3], &a[35], P);
	DFT_2_double (&a[4], &a[36], &a[5], &a[37], P);
	DFT_2_double (&a[6], &a[38], &a[7], &a[39], P);
	DFT_2_double (&a[8], &a[40], &a[9], &a[41], P);
	DFT_2_double (&a[10], &a[42], &a[11], &a[43], P);
	DFT_2_double (&a[12], &a[44], &a[13], &a[45], P);
	DFT_2_double (&a[14], &a[46], &a[15], &a[47], P);
	DFT_2_double (&a[16], &a[48], &a[17], &a[49], P);
	DFT_2_double (&a[18], &a[50], &a[19], &a[51], P);
	DFT_2_double (&a[20], &a[52], &a[21], &a[53], P);
	DFT_2_double (&a[22], &a[54], &a[23], &a[55], P);
	DFT_2_double (&a[24], &a[56], &a[25], &a[57], P);
	DFT_2_double (&a[26], &a[58], &a[27], &a[59], P);
	DFT_2_double (&a[28], &a[60], &a[29], &a[61], P);
	DFT_2_double (&a[30], &a[62], &a[31], &a[63], P);

	mult_ab_mod_p_ptr (&a[48], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[56], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[52], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[60], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[50], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[58], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[54], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[62], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[49], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[57], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[53], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[61], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[51], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[59], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[55], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[63], &base_case_omega_vec[16], P);

	DFT_2_double (&a[0], &a[16], &a[1], &a[17], P);
	DFT_2_double (&a[2], &a[18], &a[3], &a[19], P);
	DFT_2_double (&a[4], &a[20], &a[5], &a[21], P);
	DFT_2_double (&a[6], &a[22], &a[7], &a[23], P);
	DFT_2_double (&a[8], &a[24], &a[9], &a[25], P);
	DFT_2_double (&a[10], &a[26], &a[11], &a[27], P);
	DFT_2_double (&a[12], &a[28], &a[13], &a[29], P);
	DFT_2_double (&a[14], &a[30], &a[15], &a[31], P);
	DFT_2_double (&a[32], &a[48], &a[33], &a[49], P);
	DFT_2_double (&a[34], &a[50], &a[35], &a[51], P);
	DFT_2_double (&a[36], &a[52], &a[37], &a[53], P);
	DFT_2_double (&a[38], &a[54], &a[39], &a[55], P);
	DFT_2_double (&a[40], &a[56], &a[41], &a[57], P);
	DFT_2_double (&a[42], &a[58], &a[43], &a[59], P);
	DFT_2_double (&a[44], &a[60], &a[45], &a[61], P);
	DFT_2_double (&a[46], &a[62], &a[47], &a[63], P);

	mult_ab_mod_p_ptr (&a[24], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[25], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[26], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[28], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[16], P);

	mult_ab_mod_p_ptr (&a[40], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[41], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[42], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[43], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[44], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[45], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[46], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[47], &base_case_omega_vec[8], P);

	mult_ab_mod_p_ptr (&a[56], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[57], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[58], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[59], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[60], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[61], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[62], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[63], &base_case_omega_vec[24], P);

	DFT_2_double (&a[0], &a[8], &a[1], &a[9], P);
	DFT_2_double (&a[2], &a[10], &a[3], &a[11], P);
	DFT_2_double (&a[4], &a[12], &a[5], &a[13], P);
	DFT_2_double (&a[6], &a[14], &a[7], &a[15], P);

	DFT_2_double (&a[16], &a[24], &a[17], &a[25], P);
	DFT_2_double (&a[18], &a[26], &a[19], &a[27], P);
	DFT_2_double (&a[20], &a[28], &a[21], &a[29], P);
	DFT_2_double (&a[22], &a[30], &a[23], &a[31], P);

	DFT_2_double (&a[32], &a[40], &a[33], &a[41], P);
	DFT_2_double (&a[34], &a[42], &a[35], &a[43], P);
	DFT_2_double (&a[36], &a[44], &a[37], &a[45], P);
	DFT_2_double (&a[38], &a[46], &a[39], &a[47], P);

	DFT_2_double (&a[48], &a[56], &a[49], &a[57], P);
	DFT_2_double (&a[50], &a[58], &a[51], &a[59], P);
	DFT_2_double (&a[52], &a[60], &a[53], &a[61], P);
	DFT_2_double (&a[54], &a[62], &a[55], &a[63], P);

	mult_ab_mod_p_ptr (&a[12], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[16], P);

	mult_ab_mod_p_ptr (&a[20], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[21], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[22], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[8], P);

	mult_ab_mod_p_ptr (&a[28], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[24], P);

	mult_ab_mod_p_ptr (&a[36], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[37], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[38], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[39], &base_case_omega_vec[4], P);

	mult_ab_mod_p_ptr (&a[44], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr (&a[45], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr (&a[46], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr (&a[47], &base_case_omega_vec[20], P);

	mult_ab_mod_p_ptr (&a[52], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[53], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[54], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[55], &base_case_omega_vec[12], P);

	mult_ab_mod_p_ptr (&a[60], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr (&a[61], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr (&a[62], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr (&a[63], &base_case_omega_vec[28], P);

	DFT_2_double (&a[0], &a[4], &a[1], &a[5], P);
	DFT_2_double (&a[2], &a[6], &a[3], &a[7], P);
	DFT_2_double (&a[8], &a[12], &a[9], &a[13], P);
	DFT_2_double (&a[10], &a[14], &a[11], &a[15], P);
	DFT_2_double (&a[16], &a[20], &a[17], &a[21], P);
	DFT_2_double (&a[18], &a[22], &a[19], &a[23], P);
	DFT_2_double (&a[24], &a[28], &a[25], &a[29], P);
	DFT_2_double (&a[26], &a[30], &a[27], &a[31], P);
	DFT_2_double (&a[32], &a[36], &a[33], &a[37], P);
	DFT_2_double (&a[34], &a[38], &a[35], &a[39], P);
	DFT_2_double (&a[40], &a[44], &a[41], &a[45], P);
	DFT_2_double (&a[42], &a[46], &a[43], &a[47], P);
	DFT_2_double (&a[48], &a[52], &a[49], &a[53], P);
	DFT_2_double (&a[50], &a[54], &a[51], &a[55], P);
	DFT_2_double (&a[56], &a[60], &a[57], &a[61], P);
	DFT_2_double (&a[58], &a[62], &a[59], &a[63], P);

	mult_ab_mod_p_ptr (&a[6], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[10], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[14], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[18], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[19], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[22], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr (&a[26], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[30], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr (&a[34], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[35], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[38], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr (&a[39], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr (&a[42], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[43], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[46], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr (&a[47], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr (&a[50], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[51], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[54], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr (&a[55], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr (&a[58], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[59], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[62], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr (&a[63], &base_case_omega_vec[30], P);

	DFT_2_double (&a[0], &a[2], &a[1], &a[3], P);
	DFT_2_double (&a[4], &a[6], &a[5], &a[7], P);
	DFT_2_double (&a[8], &a[10], &a[9], &a[11], P);
	DFT_2_double (&a[12], &a[14], &a[13], &a[15], P);
	DFT_2_double (&a[16], &a[18], &a[17], &a[19], P);
	DFT_2_double (&a[20], &a[22], &a[21], &a[23], P);
	DFT_2_double (&a[24], &a[26], &a[25], &a[27], P);
	DFT_2_double (&a[28], &a[30], &a[29], &a[31], P);
	DFT_2_double (&a[32], &a[34], &a[33], &a[35], P);
	DFT_2_double (&a[36], &a[38], &a[37], &a[39], P);
	DFT_2_double (&a[40], &a[42], &a[41], &a[43], P);
	DFT_2_double (&a[44], &a[46], &a[45], &a[47], P);
	DFT_2_double (&a[48], &a[50], &a[49], &a[51], P);
	DFT_2_double (&a[52], &a[54], &a[53], &a[55], P);
	DFT_2_double (&a[56], &a[58], &a[57], &a[59], P);
	DFT_2_double (&a[60], &a[62], &a[61], &a[63], P);

	mult_ab_mod_p_ptr (&a[3], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr (&a[5], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr (&a[7], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr (&a[9], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr (&a[11], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr (&a[13], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr (&a[15], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr (&a[17], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr (&a[19], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr (&a[21], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr (&a[23], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr (&a[25], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr (&a[27], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr (&a[29], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr (&a[31], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr (&a[33], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr (&a[35], &base_case_omega_vec[17], P);
	mult_ab_mod_p_ptr (&a[37], &base_case_omega_vec[9], P);
	mult_ab_mod_p_ptr (&a[39], &base_case_omega_vec[25], P);
	mult_ab_mod_p_ptr (&a[41], &base_case_omega_vec[5], P);
	mult_ab_mod_p_ptr (&a[43], &base_case_omega_vec[21], P);
	mult_ab_mod_p_ptr (&a[45], &base_case_omega_vec[13], P);
	mult_ab_mod_p_ptr (&a[47], &base_case_omega_vec[29], P);
	mult_ab_mod_p_ptr (&a[49], &base_case_omega_vec[3], P);
	mult_ab_mod_p_ptr (&a[51], &base_case_omega_vec[19], P);
	mult_ab_mod_p_ptr (&a[53], &base_case_omega_vec[11], P);
	mult_ab_mod_p_ptr (&a[55], &base_case_omega_vec[27], P);
	mult_ab_mod_p_ptr (&a[57], &base_case_omega_vec[7], P);
	mult_ab_mod_p_ptr (&a[59], &base_case_omega_vec[23], P);
	mult_ab_mod_p_ptr (&a[61], &base_case_omega_vec[15], P);
	mult_ab_mod_p_ptr (&a[63], &base_case_omega_vec[31], P);

	DFT_2_double (&a[0], &a[1], &a[2], &a[3], P);
	DFT_2_double (&a[4], &a[5], &a[6], &a[7], P);
	DFT_2_double (&a[8], &a[9], &a[10], &a[11], P);
	DFT_2_double (&a[12], &a[13], &a[14], &a[15], P);
	DFT_2_double (&a[16], &a[17], &a[18], &a[19], P);
	DFT_2_double (&a[20], &a[21], &a[22], &a[23], P);
	DFT_2_double (&a[24], &a[25], &a[26], &a[27], P);
	DFT_2_double (&a[28], &a[29], &a[30], &a[31], P);
	DFT_2_double (&a[32], &a[33], &a[34], &a[35], P);
	DFT_2_double (&a[36], &a[37], &a[38], &a[39], P);
	DFT_2_double (&a[40], &a[41], &a[42], &a[43], P);
	DFT_2_double (&a[44], &a[45], &a[46], &a[47], P);
	DFT_2_double (&a[48], &a[49], &a[50], &a[51], P);
	DFT_2_double (&a[52], &a[53], &a[54], &a[55], P);
	DFT_2_double (&a[56], &a[57], &a[58], &a[59], P);
	DFT_2_double (&a[60], &a[61], &a[62], &a[63], P);

	swap (&a[1], &a[32]);
	swap (&a[2], &a[16]);
	swap (&a[3], &a[48]);
	swap (&a[4], &a[8]);
	swap (&a[5], &a[40]);
	swap (&a[6], &a[24]);
	swap (&a[7], &a[56]);
	swap (&a[9], &a[36]);
	swap (&a[10], &a[20]);
	swap (&a[11], &a[52]);
	swap (&a[13], &a[44]);
	swap (&a[14], &a[28]);
	swap (&a[15], &a[60]);
	swap (&a[17], &a[34]);
	swap (&a[19], &a[50]);
	swap (&a[21], &a[42]);
	swap (&a[22], &a[26]);
	swap (&a[23], &a[58]);
	swap (&a[25], &a[38]);
	swap (&a[27], &a[54]);
	swap (&a[29], &a[46]);
	swap (&a[31], &a[62]);
	swap (&a[35], &a[49]);
	swap (&a[37], &a[41]);
	swap (&a[39], &a[57]);
	swap (&a[43], &a[53]);
	swap (&a[47], &a[61]);
	swap (&a[55], &a[59]);
      }
  }

///////////////////////////////////////////////////////////

void
DFT_128 (int n, usfixn64* __restrict__ x,
	 const usfixn64 *__restrict__ base_case_omega_vec,
	 const montgomery_triple* P)
{
  //#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n;i++)
      {
	usfixn64 *a = &x[i * 128];
	DFT_2_double(&a[0], &a[64], &a[32], &a[96], P);
	DFT_2_double(&a[1], &a[65], &a[33], &a[97], P);
	DFT_2_double(&a[2], &a[66], &a[34], &a[98], P);
	DFT_2_double(&a[3], &a[67], &a[35], &a[99], P);
	DFT_2_double(&a[4], &a[68], &a[36], &a[100], P);
	DFT_2_double(&a[5], &a[69], &a[37], &a[101], P);
	DFT_2_double(&a[6], &a[70], &a[38], &a[102], P);
	DFT_2_double(&a[7], &a[71], &a[39], &a[103], P);
	DFT_2_double(&a[8], &a[72], &a[40], &a[104], P);
	DFT_2_double(&a[9], &a[73], &a[41], &a[105], P);
	DFT_2_double(&a[10], &a[74], &a[42], &a[106], P);
	DFT_2_double(&a[11], &a[75], &a[43], &a[107], P);
	DFT_2_double(&a[12], &a[76], &a[44], &a[108], P);
	DFT_2_double(&a[13], &a[77], &a[45], &a[109], P);
	DFT_2_double(&a[14], &a[78], &a[46], &a[110], P);
	DFT_2_double(&a[15], &a[79], &a[47], &a[111], P);
	DFT_2_double(&a[16], &a[80], &a[48], &a[112], P);
	DFT_2_double(&a[17], &a[81], &a[49], &a[113], P);
	DFT_2_double(&a[18], &a[82], &a[50], &a[114], P);
	DFT_2_double(&a[19], &a[83], &a[51], &a[115], P);
	DFT_2_double(&a[20], &a[84], &a[52], &a[116], P);
	DFT_2_double(&a[21], &a[85], &a[53], &a[117], P);
	DFT_2_double(&a[22], &a[86], &a[54], &a[118], P);
	DFT_2_double(&a[23], &a[87], &a[55], &a[119], P);
	DFT_2_double(&a[24], &a[88], &a[56], &a[120], P);
	DFT_2_double(&a[25], &a[89], &a[57], &a[121], P);
	DFT_2_double(&a[26], &a[90], &a[58], &a[122], P);
	DFT_2_double(&a[27], &a[91], &a[59], &a[123], P);
	DFT_2_double(&a[28], &a[92], &a[60], &a[124], P);
	DFT_2_double(&a[29], &a[93], &a[61], &a[125], P);
	DFT_2_double(&a[30], &a[94], &a[62], &a[126], P);
	DFT_2_double(&a[31], &a[95], &a[63], &a[127], P);
	mult_ab_mod_p_ptr(&a[96], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[97], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[98], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[99], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[100], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[101], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[102], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[104], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[105], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[106], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[108], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[112], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[113], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[114], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[116], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[120], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[32], P);
	DFT_2_double(&a[0], &a[32], &a[64], &a[96], P);
	DFT_2_double(&a[1], &a[33], &a[65], &a[97], P);
	DFT_2_double(&a[2], &a[34], &a[66], &a[98], P);
	DFT_2_double(&a[3], &a[35], &a[67], &a[99], P);
	DFT_2_double(&a[4], &a[36], &a[68], &a[100], P);
	DFT_2_double(&a[5], &a[37], &a[69], &a[101], P);
	DFT_2_double(&a[6], &a[38], &a[70], &a[102], P);
	DFT_2_double(&a[7], &a[39], &a[71], &a[103], P);
	DFT_2_double(&a[8], &a[40], &a[72], &a[104], P);
	DFT_2_double(&a[9], &a[41], &a[73], &a[105], P);
	DFT_2_double(&a[10], &a[42], &a[74], &a[106], P);
	DFT_2_double(&a[11], &a[43], &a[75], &a[107], P);
	DFT_2_double(&a[12], &a[44], &a[76], &a[108], P);
	DFT_2_double(&a[13], &a[45], &a[77], &a[109], P);
	DFT_2_double(&a[14], &a[46], &a[78], &a[110], P);
	DFT_2_double(&a[15], &a[47], &a[79], &a[111], P);
	DFT_2_double(&a[16], &a[48], &a[80], &a[112], P);
	DFT_2_double(&a[17], &a[49], &a[81], &a[113], P);
	DFT_2_double(&a[18], &a[50], &a[82], &a[114], P);
	DFT_2_double(&a[19], &a[51], &a[83], &a[115], P);
	DFT_2_double(&a[20], &a[52], &a[84], &a[116], P);
	DFT_2_double(&a[21], &a[53], &a[85], &a[117], P);
	DFT_2_double(&a[22], &a[54], &a[86], &a[118], P);
	DFT_2_double(&a[23], &a[55], &a[87], &a[119], P);
	DFT_2_double(&a[24], &a[56], &a[88], &a[120], P);
	DFT_2_double(&a[25], &a[57], &a[89], &a[121], P);
	DFT_2_double(&a[26], &a[58], &a[90], &a[122], P);
	DFT_2_double(&a[27], &a[59], &a[91], &a[123], P);
	DFT_2_double(&a[28], &a[60], &a[92], &a[124], P);
	DFT_2_double(&a[29], &a[61], &a[93], &a[125], P);
	DFT_2_double(&a[30], &a[62], &a[94], &a[126], P);
	DFT_2_double(&a[31], &a[63], &a[95], &a[127], P);
	mult_ab_mod_p_ptr(&a[48], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[49], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[50], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[51], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[52], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[53], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[54], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[56], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[57], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[58], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[60], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[80], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[81], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[82], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[83], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[84], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[85], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[86], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[88], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[89], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[90], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[92], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[112], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[113], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[114], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[116], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[120], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[48], P);
	DFT_2_double(&a[0], &a[16], &a[64], &a[80], P);
	DFT_2_double(&a[1], &a[17], &a[65], &a[81], P);
	DFT_2_double(&a[2], &a[18], &a[66], &a[82], P);
	DFT_2_double(&a[3], &a[19], &a[67], &a[83], P);
	DFT_2_double(&a[4], &a[20], &a[68], &a[84], P);
	DFT_2_double(&a[5], &a[21], &a[69], &a[85], P);
	DFT_2_double(&a[6], &a[22], &a[70], &a[86], P);
	DFT_2_double(&a[7], &a[23], &a[71], &a[87], P);
	DFT_2_double(&a[8], &a[24], &a[72], &a[88], P);
	DFT_2_double(&a[9], &a[25], &a[73], &a[89], P);
	DFT_2_double(&a[10], &a[26], &a[74], &a[90], P);
	DFT_2_double(&a[11], &a[27], &a[75], &a[91], P);
	DFT_2_double(&a[12], &a[28], &a[76], &a[92], P);
	DFT_2_double(&a[13], &a[29], &a[77], &a[93], P);
	DFT_2_double(&a[14], &a[30], &a[78], &a[94], P);
	DFT_2_double(&a[15], &a[31], &a[79], &a[95], P);
	DFT_2_double(&a[32], &a[48], &a[96], &a[112], P);
	DFT_2_double(&a[33], &a[49], &a[97], &a[113], P);
	DFT_2_double(&a[34], &a[50], &a[98], &a[114], P);
	DFT_2_double(&a[35], &a[51], &a[99], &a[115], P);
	DFT_2_double(&a[36], &a[52], &a[100], &a[116], P);
	DFT_2_double(&a[37], &a[53], &a[101], &a[117], P);
	DFT_2_double(&a[38], &a[54], &a[102], &a[118], P);
	DFT_2_double(&a[39], &a[55], &a[103], &a[119], P);
	DFT_2_double(&a[40], &a[56], &a[104], &a[120], P);
	DFT_2_double(&a[41], &a[57], &a[105], &a[121], P);
	DFT_2_double(&a[42], &a[58], &a[106], &a[122], P);
	DFT_2_double(&a[43], &a[59], &a[107], &a[123], P);
	DFT_2_double(&a[44], &a[60], &a[108], &a[124], P);
	DFT_2_double(&a[45], &a[61], &a[109], &a[125], P);
	DFT_2_double(&a[46], &a[62], &a[110], &a[126], P);
	DFT_2_double(&a[47], &a[63], &a[111], &a[127], P);
	mult_ab_mod_p_ptr(&a[24], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[25], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[26], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[27], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[28], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[29], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[30], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[40], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[41], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[42], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[43], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[44], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[45], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[46], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[56], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[57], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[58], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[60], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[72], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[73], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[74], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[75], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[76], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[77], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[78], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[88], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[89], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[90], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[92], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[104], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[105], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[106], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[108], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[120], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[56], P);
	DFT_2_double(&a[0], &a[8], &a[64], &a[72], P);
	DFT_2_double(&a[1], &a[9], &a[65], &a[73], P);
	DFT_2_double(&a[2], &a[10], &a[66], &a[74], P);
	DFT_2_double(&a[3], &a[11], &a[67], &a[75], P);
	DFT_2_double(&a[4], &a[12], &a[68], &a[76], P);
	DFT_2_double(&a[5], &a[13], &a[69], &a[77], P);
	DFT_2_double(&a[6], &a[14], &a[70], &a[78], P);
	DFT_2_double(&a[7], &a[15], &a[71], &a[79], P);
	DFT_2_double(&a[16], &a[24], &a[80], &a[88], P);
	DFT_2_double(&a[17], &a[25], &a[81], &a[89], P);
	DFT_2_double(&a[18], &a[26], &a[82], &a[90], P);
	DFT_2_double(&a[19], &a[27], &a[83], &a[91], P);
	DFT_2_double(&a[20], &a[28], &a[84], &a[92], P);
	DFT_2_double(&a[21], &a[29], &a[85], &a[93], P);
	DFT_2_double(&a[22], &a[30], &a[86], &a[94], P);
	DFT_2_double(&a[23], &a[31], &a[87], &a[95], P);
	DFT_2_double(&a[32], &a[40], &a[96], &a[104], P);
	DFT_2_double(&a[33], &a[41], &a[97], &a[105], P);
	DFT_2_double(&a[34], &a[42], &a[98], &a[106], P);
	DFT_2_double(&a[35], &a[43], &a[99], &a[107], P);
	DFT_2_double(&a[36], &a[44], &a[100], &a[108], P);
	DFT_2_double(&a[37], &a[45], &a[101], &a[109], P);
	DFT_2_double(&a[38], &a[46], &a[102], &a[110], P);
	DFT_2_double(&a[39], &a[47], &a[103], &a[111], P);
	DFT_2_double(&a[48], &a[56], &a[112], &a[120], P);
	DFT_2_double(&a[49], &a[57], &a[113], &a[121], P);
	DFT_2_double(&a[50], &a[58], &a[114], &a[122], P);
	DFT_2_double(&a[51], &a[59], &a[115], &a[123], P);
	DFT_2_double(&a[52], &a[60], &a[116], &a[124], P);
	DFT_2_double(&a[53], &a[61], &a[117], &a[125], P);
	DFT_2_double(&a[54], &a[62], &a[118], &a[126], P);
	DFT_2_double(&a[55], &a[63], &a[119], &a[127], P);
	mult_ab_mod_p_ptr(&a[12], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[13], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[14], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[15], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[20], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[21], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[22], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[23], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[28], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[29], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[30], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[36], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[37], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[38], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[39], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[44], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[45], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[46], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[52], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[53], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[54], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[60], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[68], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[69], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[70], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[71], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[76], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[77], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[78], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[84], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[85], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[86], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[92], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[100], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[101], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[102], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[108], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[116], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[60], P);
	DFT_2_double(&a[0], &a[4], &a[64], &a[68], P);
	DFT_2_double(&a[1], &a[5], &a[65], &a[69], P);
	DFT_2_double(&a[2], &a[6], &a[66], &a[70], P);
	DFT_2_double(&a[3], &a[7], &a[67], &a[71], P);
	DFT_2_double(&a[8], &a[12], &a[72], &a[76], P);
	DFT_2_double(&a[9], &a[13], &a[73], &a[77], P);
	DFT_2_double(&a[10], &a[14], &a[74], &a[78], P);
	DFT_2_double(&a[11], &a[15], &a[75], &a[79], P);
	DFT_2_double(&a[16], &a[20], &a[80], &a[84], P);
	DFT_2_double(&a[17], &a[21], &a[81], &a[85], P);
	DFT_2_double(&a[18], &a[22], &a[82], &a[86], P);
	DFT_2_double(&a[19], &a[23], &a[83], &a[87], P);
	DFT_2_double(&a[24], &a[28], &a[88], &a[92], P);
	DFT_2_double(&a[25], &a[29], &a[89], &a[93], P);
	DFT_2_double(&a[26], &a[30], &a[90], &a[94], P);
	DFT_2_double(&a[27], &a[31], &a[91], &a[95], P);
	DFT_2_double(&a[32], &a[36], &a[96], &a[100], P);
	DFT_2_double(&a[33], &a[37], &a[97], &a[101], P);
	DFT_2_double(&a[34], &a[38], &a[98], &a[102], P);
	DFT_2_double(&a[35], &a[39], &a[99], &a[103], P);
	DFT_2_double(&a[40], &a[44], &a[104], &a[108], P);
	DFT_2_double(&a[41], &a[45], &a[105], &a[109], P);
	DFT_2_double(&a[42], &a[46], &a[106], &a[110], P);
	DFT_2_double(&a[43], &a[47], &a[107], &a[111], P);
	DFT_2_double(&a[48], &a[52], &a[112], &a[116], P);
	DFT_2_double(&a[49], &a[53], &a[113], &a[117], P);
	DFT_2_double(&a[50], &a[54], &a[114], &a[118], P);
	DFT_2_double(&a[51], &a[55], &a[115], &a[119], P);
	DFT_2_double(&a[56], &a[60], &a[120], &a[124], P);
	DFT_2_double(&a[57], &a[61], &a[121], &a[125], P);
	DFT_2_double(&a[58], &a[62], &a[122], &a[126], P);
	DFT_2_double(&a[59], &a[63], &a[123], &a[127], P);
	mult_ab_mod_p_ptr(&a[6], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[7], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[10], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[11], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[14], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[15], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[18], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[19], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[22], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[23], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[26], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[27], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[30], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[34], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[35], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[38], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[39], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[42], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[43], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[46], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[50], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[51], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[54], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[58], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[66], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr(&a[67], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr(&a[70], &base_case_omega_vec[34], P);
	mult_ab_mod_p_ptr(&a[71], &base_case_omega_vec[34], P);
	mult_ab_mod_p_ptr(&a[74], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr(&a[75], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr(&a[78], &base_case_omega_vec[50], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[50], P);
	mult_ab_mod_p_ptr(&a[82], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr(&a[83], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr(&a[86], &base_case_omega_vec[42], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[42], P);
	mult_ab_mod_p_ptr(&a[90], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[58], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[58], P);
	mult_ab_mod_p_ptr(&a[98], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr(&a[99], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr(&a[102], &base_case_omega_vec[38], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[38], P);
	mult_ab_mod_p_ptr(&a[106], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[54], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[54], P);
	mult_ab_mod_p_ptr(&a[114], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[46], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[46], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[62], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[62], P);
	DFT_2_double(&a[0], &a[2], &a[64], &a[66], P);
	DFT_2_double(&a[1], &a[3], &a[65], &a[67], P);
	DFT_2_double(&a[4], &a[6], &a[68], &a[70], P);
	DFT_2_double(&a[5], &a[7], &a[69], &a[71], P);
	DFT_2_double(&a[8], &a[10], &a[72], &a[74], P);
	DFT_2_double(&a[9], &a[11], &a[73], &a[75], P);
	DFT_2_double(&a[12], &a[14], &a[76], &a[78], P);
	DFT_2_double(&a[13], &a[15], &a[77], &a[79], P);
	DFT_2_double(&a[16], &a[18], &a[80], &a[82], P);
	DFT_2_double(&a[17], &a[19], &a[81], &a[83], P);
	DFT_2_double(&a[20], &a[22], &a[84], &a[86], P);
	DFT_2_double(&a[21], &a[23], &a[85], &a[87], P);
	DFT_2_double(&a[24], &a[26], &a[88], &a[90], P);
	DFT_2_double(&a[25], &a[27], &a[89], &a[91], P);
	DFT_2_double(&a[28], &a[30], &a[92], &a[94], P);
	DFT_2_double(&a[29], &a[31], &a[93], &a[95], P);
	DFT_2_double(&a[32], &a[34], &a[96], &a[98], P);
	DFT_2_double(&a[33], &a[35], &a[97], &a[99], P);
	DFT_2_double(&a[36], &a[38], &a[100], &a[102], P);
	DFT_2_double(&a[37], &a[39], &a[101], &a[103], P);
	DFT_2_double(&a[40], &a[42], &a[104], &a[106], P);
	DFT_2_double(&a[41], &a[43], &a[105], &a[107], P);
	DFT_2_double(&a[44], &a[46], &a[108], &a[110], P);
	DFT_2_double(&a[45], &a[47], &a[109], &a[111], P);
	DFT_2_double(&a[48], &a[50], &a[112], &a[114], P);
	DFT_2_double(&a[49], &a[51], &a[113], &a[115], P);
	DFT_2_double(&a[52], &a[54], &a[116], &a[118], P);
	DFT_2_double(&a[53], &a[55], &a[117], &a[119], P);
	DFT_2_double(&a[56], &a[58], &a[120], &a[122], P);
	DFT_2_double(&a[57], &a[59], &a[121], &a[123], P);
	DFT_2_double(&a[60], &a[62], &a[124], &a[126], P);
	DFT_2_double(&a[61], &a[63], &a[125], &a[127], P);
	mult_ab_mod_p_ptr(&a[3], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[5], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[7], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[9], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[11], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[13], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[15], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[17], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[19], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[21], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[23], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[25], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[27], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[29], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[33], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr(&a[35], &base_case_omega_vec[34], P);
	mult_ab_mod_p_ptr(&a[37], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr(&a[39], &base_case_omega_vec[50], P);
	mult_ab_mod_p_ptr(&a[41], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr(&a[43], &base_case_omega_vec[42], P);
	mult_ab_mod_p_ptr(&a[45], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[58], P);
	mult_ab_mod_p_ptr(&a[49], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr(&a[51], &base_case_omega_vec[38], P);
	mult_ab_mod_p_ptr(&a[53], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[54], P);
	mult_ab_mod_p_ptr(&a[57], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[46], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[62], P);
	mult_ab_mod_p_ptr(&a[65], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr(&a[67], &base_case_omega_vec[33], P);
	mult_ab_mod_p_ptr(&a[69], &base_case_omega_vec[17], P);
	mult_ab_mod_p_ptr(&a[71], &base_case_omega_vec[49], P);
	mult_ab_mod_p_ptr(&a[73], &base_case_omega_vec[9], P);
	mult_ab_mod_p_ptr(&a[75], &base_case_omega_vec[41], P);
	mult_ab_mod_p_ptr(&a[77], &base_case_omega_vec[25], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[57], P);
	mult_ab_mod_p_ptr(&a[81], &base_case_omega_vec[5], P);
	mult_ab_mod_p_ptr(&a[83], &base_case_omega_vec[37], P);
	mult_ab_mod_p_ptr(&a[85], &base_case_omega_vec[21], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[53], P);
	mult_ab_mod_p_ptr(&a[89], &base_case_omega_vec[13], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[45], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[29], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[61], P);
	mult_ab_mod_p_ptr(&a[97], &base_case_omega_vec[3], P);
	mult_ab_mod_p_ptr(&a[99], &base_case_omega_vec[35], P);
	mult_ab_mod_p_ptr(&a[101], &base_case_omega_vec[19], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[51], P);
	mult_ab_mod_p_ptr(&a[105], &base_case_omega_vec[11], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[43], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[27], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[59], P);
	mult_ab_mod_p_ptr(&a[113], &base_case_omega_vec[7], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[39], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[23], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[55], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[15], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[47], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[31], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[63], P);
	DFT_2_double(&a[0], &a[1], &a[64], &a[65], P);
	DFT_2_double(&a[2], &a[3], &a[66], &a[67], P);
	DFT_2_double(&a[4], &a[5], &a[68], &a[69], P);
	DFT_2_double(&a[6], &a[7], &a[70], &a[71], P);
	DFT_2_double(&a[8], &a[9], &a[72], &a[73], P);
	DFT_2_double(&a[10], &a[11], &a[74], &a[75], P);
	DFT_2_double(&a[12], &a[13], &a[76], &a[77], P);
	DFT_2_double(&a[14], &a[15], &a[78], &a[79], P);
	DFT_2_double(&a[16], &a[17], &a[80], &a[81], P);
	DFT_2_double(&a[18], &a[19], &a[82], &a[83], P);
	DFT_2_double(&a[20], &a[21], &a[84], &a[85], P);
	DFT_2_double(&a[22], &a[23], &a[86], &a[87], P);
	DFT_2_double(&a[24], &a[25], &a[88], &a[89], P);
	DFT_2_double(&a[26], &a[27], &a[90], &a[91], P);
	DFT_2_double(&a[28], &a[29], &a[92], &a[93], P);
	DFT_2_double(&a[30], &a[31], &a[94], &a[95], P);
	DFT_2_double(&a[32], &a[33], &a[96], &a[97], P);
	DFT_2_double(&a[34], &a[35], &a[98], &a[99], P);
	DFT_2_double(&a[36], &a[37], &a[100], &a[101], P);
	DFT_2_double(&a[38], &a[39], &a[102], &a[103], P);
	DFT_2_double(&a[40], &a[41], &a[104], &a[105], P);
	DFT_2_double(&a[42], &a[43], &a[106], &a[107], P);
	DFT_2_double(&a[44], &a[45], &a[108], &a[109], P);
	DFT_2_double(&a[46], &a[47], &a[110], &a[111], P);
	DFT_2_double(&a[48], &a[49], &a[112], &a[113], P);
	DFT_2_double(&a[50], &a[51], &a[114], &a[115], P);
	DFT_2_double(&a[52], &a[53], &a[116], &a[117], P);
	DFT_2_double(&a[54], &a[55], &a[118], &a[119], P);
	DFT_2_double(&a[56], &a[57], &a[120], &a[121], P);
	DFT_2_double(&a[58], &a[59], &a[122], &a[123], P);
	DFT_2_double(&a[60], &a[61], &a[124], &a[125], P);
	DFT_2_double(&a[62], &a[63], &a[126], &a[127], P);
	swap (&a[1], &a[64]);
	swap (&a[2], &a[32]);
	swap (&a[3], &a[96]);
	swap (&a[4], &a[16]);
	swap (&a[5], &a[80]);
	swap (&a[6], &a[48]);
	swap (&a[7], &a[112]);
	swap (&a[9], &a[72]);
	swap (&a[10], &a[40]);
	swap (&a[11], &a[104]);
	swap (&a[12], &a[24]);
	swap (&a[13], &a[88]);
	swap (&a[14], &a[56]);
	swap (&a[15], &a[120]);
	swap (&a[17], &a[68]);
	swap (&a[18], &a[36]);
	swap (&a[19], &a[100]);
	swap (&a[21], &a[84]);
	swap (&a[22], &a[52]);
	swap (&a[23], &a[116]);
	swap (&a[25], &a[76]);
	swap (&a[26], &a[44]);
	swap (&a[27], &a[108]);
	swap (&a[29], &a[92]);
	swap (&a[30], &a[60]);
	swap (&a[31], &a[124]);
	swap (&a[33], &a[66]);
	swap (&a[35], &a[98]);
	swap (&a[37], &a[82]);
	swap (&a[38], &a[50]);
	swap (&a[39], &a[114]);
	swap (&a[41], &a[74]);
	swap (&a[43], &a[106]);
	swap (&a[45], &a[90]);
	swap (&a[46], &a[58]);
	swap (&a[47], &a[122]);
	swap (&a[49], &a[70]);
	swap (&a[51], &a[102]);
	swap (&a[53], &a[86]);
	swap (&a[55], &a[118]);
	swap (&a[57], &a[78]);
	swap (&a[59], &a[110]);
	swap (&a[61], &a[94]);
	swap (&a[63], &a[126]);
	swap (&a[67], &a[97]);
	swap (&a[69], &a[81]);
	swap (&a[71], &a[113]);
	swap (&a[75], &a[105]);
	swap (&a[77], &a[89]);
	swap (&a[79], &a[121]);
	swap (&a[83], &a[101]);
	swap (&a[87], &a[117]);
	swap (&a[91], &a[109]);
	swap (&a[95], &a[125]);
	swap (&a[103], &a[115]);
	swap (&a[111], &a[123]);
      }
  }

///////////////////////////////////////////////////////////

void
DFT_256 (int n, usfixn64* __restrict__ x,
	 const usfixn64 *__restrict__ base_case_omega_vec,
	 const montgomery_triple* P)
{
  //#pragma cilk_grainsize = (SMALL_DFT_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = MAX(n/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (int i = 0; i < n;i++)
      {
	usfixn64 *a = &x[i * 256];
	DFT_2_double(&a[0], &a[128], &a[64], &a[192], P);
	DFT_2_double(&a[1], &a[129], &a[65], &a[193], P);
	DFT_2_double(&a[2], &a[130], &a[66], &a[194], P);
	DFT_2_double(&a[3], &a[131], &a[67], &a[195], P);
	DFT_2_double(&a[4], &a[132], &a[68], &a[196], P);
	DFT_2_double(&a[5], &a[133], &a[69], &a[197], P);
	DFT_2_double(&a[6], &a[134], &a[70], &a[198], P);
	DFT_2_double(&a[7], &a[135], &a[71], &a[199], P);
	DFT_2_double(&a[8], &a[136], &a[72], &a[200], P);
	DFT_2_double(&a[9], &a[137], &a[73], &a[201], P);
	DFT_2_double(&a[10], &a[138], &a[74], &a[202], P);
	DFT_2_double(&a[11], &a[139], &a[75], &a[203], P);
	DFT_2_double(&a[12], &a[140], &a[76], &a[204], P);
	DFT_2_double(&a[13], &a[141], &a[77], &a[205], P);
	DFT_2_double(&a[14], &a[142], &a[78], &a[206], P);
	DFT_2_double(&a[15], &a[143], &a[79], &a[207], P);
	DFT_2_double(&a[16], &a[144], &a[80], &a[208], P);
	DFT_2_double(&a[17], &a[145], &a[81], &a[209], P);
	DFT_2_double(&a[18], &a[146], &a[82], &a[210], P);
	DFT_2_double(&a[19], &a[147], &a[83], &a[211], P);
	DFT_2_double(&a[20], &a[148], &a[84], &a[212], P);
	DFT_2_double(&a[21], &a[149], &a[85], &a[213], P);
	DFT_2_double(&a[22], &a[150], &a[86], &a[214], P);
	DFT_2_double(&a[23], &a[151], &a[87], &a[215], P);
	DFT_2_double(&a[24], &a[152], &a[88], &a[216], P);
	DFT_2_double(&a[25], &a[153], &a[89], &a[217], P);
	DFT_2_double(&a[26], &a[154], &a[90], &a[218], P);
	DFT_2_double(&a[27], &a[155], &a[91], &a[219], P);
	DFT_2_double(&a[28], &a[156], &a[92], &a[220], P);
	DFT_2_double(&a[29], &a[157], &a[93], &a[221], P);
	DFT_2_double(&a[30], &a[158], &a[94], &a[222], P);
	DFT_2_double(&a[31], &a[159], &a[95], &a[223], P);
	DFT_2_double(&a[32], &a[160], &a[96], &a[224], P);
	DFT_2_double(&a[33], &a[161], &a[97], &a[225], P);
	DFT_2_double(&a[34], &a[162], &a[98], &a[226], P);
	DFT_2_double(&a[35], &a[163], &a[99], &a[227], P);
	DFT_2_double(&a[36], &a[164], &a[100], &a[228], P);
	DFT_2_double(&a[37], &a[165], &a[101], &a[229], P);
	DFT_2_double(&a[38], &a[166], &a[102], &a[230], P);
	DFT_2_double(&a[39], &a[167], &a[103], &a[231], P);
	DFT_2_double(&a[40], &a[168], &a[104], &a[232], P);
	DFT_2_double(&a[41], &a[169], &a[105], &a[233], P);
	DFT_2_double(&a[42], &a[170], &a[106], &a[234], P);
	DFT_2_double(&a[43], &a[171], &a[107], &a[235], P);
	DFT_2_double(&a[44], &a[172], &a[108], &a[236], P);
	DFT_2_double(&a[45], &a[173], &a[109], &a[237], P);
	DFT_2_double(&a[46], &a[174], &a[110], &a[238], P);
	DFT_2_double(&a[47], &a[175], &a[111], &a[239], P);
	DFT_2_double(&a[48], &a[176], &a[112], &a[240], P);
	DFT_2_double(&a[49], &a[177], &a[113], &a[241], P);
	DFT_2_double(&a[50], &a[178], &a[114], &a[242], P);
	DFT_2_double(&a[51], &a[179], &a[115], &a[243], P);
	DFT_2_double(&a[52], &a[180], &a[116], &a[244], P);
	DFT_2_double(&a[53], &a[181], &a[117], &a[245], P);
	DFT_2_double(&a[54], &a[182], &a[118], &a[246], P);
	DFT_2_double(&a[55], &a[183], &a[119], &a[247], P);
	DFT_2_double(&a[56], &a[184], &a[120], &a[248], P);
	DFT_2_double(&a[57], &a[185], &a[121], &a[249], P);
	DFT_2_double(&a[58], &a[186], &a[122], &a[250], P);
	DFT_2_double(&a[59], &a[187], &a[123], &a[251], P);
	DFT_2_double(&a[60], &a[188], &a[124], &a[252], P);
	DFT_2_double(&a[61], &a[189], &a[125], &a[253], P);
	DFT_2_double(&a[62], &a[190], &a[126], &a[254], P);
	DFT_2_double(&a[63], &a[191], &a[127], &a[255], P);
	mult_ab_mod_p_ptr(&a[192], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[193], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[194], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[195], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[196], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[197], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[198], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[199], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[200], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[201], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[202], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[203], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[204], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[205], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[206], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[207], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[208], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[209], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[210], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[211], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[212], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[213], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[214], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[215], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[216], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[217], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[218], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[219], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[220], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[221], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[222], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[223], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[224], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[225], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[226], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[227], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[228], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[229], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[230], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[231], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[232], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[233], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[234], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[235], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[236], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[237], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[238], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[239], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[240], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[241], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[242], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[243], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[244], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[245], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[246], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[247], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[248], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[249], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[250], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[251], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[252], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[253], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[254], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[64], P);
	DFT_2_double(&a[0], &a[64], &a[128], &a[192], P);
	DFT_2_double(&a[1], &a[65], &a[129], &a[193], P);
	DFT_2_double(&a[2], &a[66], &a[130], &a[194], P);
	DFT_2_double(&a[3], &a[67], &a[131], &a[195], P);
	DFT_2_double(&a[4], &a[68], &a[132], &a[196], P);
	DFT_2_double(&a[5], &a[69], &a[133], &a[197], P);
	DFT_2_double(&a[6], &a[70], &a[134], &a[198], P);
	DFT_2_double(&a[7], &a[71], &a[135], &a[199], P);
	DFT_2_double(&a[8], &a[72], &a[136], &a[200], P);
	DFT_2_double(&a[9], &a[73], &a[137], &a[201], P);
	DFT_2_double(&a[10], &a[74], &a[138], &a[202], P);
	DFT_2_double(&a[11], &a[75], &a[139], &a[203], P);
	DFT_2_double(&a[12], &a[76], &a[140], &a[204], P);
	DFT_2_double(&a[13], &a[77], &a[141], &a[205], P);
	DFT_2_double(&a[14], &a[78], &a[142], &a[206], P);
	DFT_2_double(&a[15], &a[79], &a[143], &a[207], P);
	DFT_2_double(&a[16], &a[80], &a[144], &a[208], P);
	DFT_2_double(&a[17], &a[81], &a[145], &a[209], P);
	DFT_2_double(&a[18], &a[82], &a[146], &a[210], P);
	DFT_2_double(&a[19], &a[83], &a[147], &a[211], P);
	DFT_2_double(&a[20], &a[84], &a[148], &a[212], P);
	DFT_2_double(&a[21], &a[85], &a[149], &a[213], P);
	DFT_2_double(&a[22], &a[86], &a[150], &a[214], P);
	DFT_2_double(&a[23], &a[87], &a[151], &a[215], P);
	DFT_2_double(&a[24], &a[88], &a[152], &a[216], P);
	DFT_2_double(&a[25], &a[89], &a[153], &a[217], P);
	DFT_2_double(&a[26], &a[90], &a[154], &a[218], P);
	DFT_2_double(&a[27], &a[91], &a[155], &a[219], P);
	DFT_2_double(&a[28], &a[92], &a[156], &a[220], P);
	DFT_2_double(&a[29], &a[93], &a[157], &a[221], P);
	DFT_2_double(&a[30], &a[94], &a[158], &a[222], P);
	DFT_2_double(&a[31], &a[95], &a[159], &a[223], P);
	DFT_2_double(&a[32], &a[96], &a[160], &a[224], P);
	DFT_2_double(&a[33], &a[97], &a[161], &a[225], P);
	DFT_2_double(&a[34], &a[98], &a[162], &a[226], P);
	DFT_2_double(&a[35], &a[99], &a[163], &a[227], P);
	DFT_2_double(&a[36], &a[100], &a[164], &a[228], P);
	DFT_2_double(&a[37], &a[101], &a[165], &a[229], P);
	DFT_2_double(&a[38], &a[102], &a[166], &a[230], P);
	DFT_2_double(&a[39], &a[103], &a[167], &a[231], P);
	DFT_2_double(&a[40], &a[104], &a[168], &a[232], P);
	DFT_2_double(&a[41], &a[105], &a[169], &a[233], P);
	DFT_2_double(&a[42], &a[106], &a[170], &a[234], P);
	DFT_2_double(&a[43], &a[107], &a[171], &a[235], P);
	DFT_2_double(&a[44], &a[108], &a[172], &a[236], P);
	DFT_2_double(&a[45], &a[109], &a[173], &a[237], P);
	DFT_2_double(&a[46], &a[110], &a[174], &a[238], P);
	DFT_2_double(&a[47], &a[111], &a[175], &a[239], P);
	DFT_2_double(&a[48], &a[112], &a[176], &a[240], P);
	DFT_2_double(&a[49], &a[113], &a[177], &a[241], P);
	DFT_2_double(&a[50], &a[114], &a[178], &a[242], P);
	DFT_2_double(&a[51], &a[115], &a[179], &a[243], P);
	DFT_2_double(&a[52], &a[116], &a[180], &a[244], P);
	DFT_2_double(&a[53], &a[117], &a[181], &a[245], P);
	DFT_2_double(&a[54], &a[118], &a[182], &a[246], P);
	DFT_2_double(&a[55], &a[119], &a[183], &a[247], P);
	DFT_2_double(&a[56], &a[120], &a[184], &a[248], P);
	DFT_2_double(&a[57], &a[121], &a[185], &a[249], P);
	DFT_2_double(&a[58], &a[122], &a[186], &a[250], P);
	DFT_2_double(&a[59], &a[123], &a[187], &a[251], P);
	DFT_2_double(&a[60], &a[124], &a[188], &a[252], P);
	DFT_2_double(&a[61], &a[125], &a[189], &a[253], P);
	DFT_2_double(&a[62], &a[126], &a[190], &a[254], P);
	DFT_2_double(&a[63], &a[127], &a[191], &a[255], P);
	mult_ab_mod_p_ptr(&a[96], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[97], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[98], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[99], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[100], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[101], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[102], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[104], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[105], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[106], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[108], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[112], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[113], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[114], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[116], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[120], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[160], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[161], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[162], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[163], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[164], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[165], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[166], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[167], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[168], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[169], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[170], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[171], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[172], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[173], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[174], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[175], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[176], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[177], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[178], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[179], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[180], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[181], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[182], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[183], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[184], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[185], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[186], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[187], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[188], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[189], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[190], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[191], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[224], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[225], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[226], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[227], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[228], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[229], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[230], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[231], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[232], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[233], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[234], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[235], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[236], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[237], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[238], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[239], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[240], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[241], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[242], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[243], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[244], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[245], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[246], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[247], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[248], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[249], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[250], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[251], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[252], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[253], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[254], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[96], P);
	DFT_2_double(&a[0], &a[32], &a[128], &a[160], P);
	DFT_2_double(&a[1], &a[33], &a[129], &a[161], P);
	DFT_2_double(&a[2], &a[34], &a[130], &a[162], P);
	DFT_2_double(&a[3], &a[35], &a[131], &a[163], P);
	DFT_2_double(&a[4], &a[36], &a[132], &a[164], P);
	DFT_2_double(&a[5], &a[37], &a[133], &a[165], P);
	DFT_2_double(&a[6], &a[38], &a[134], &a[166], P);
	DFT_2_double(&a[7], &a[39], &a[135], &a[167], P);
	DFT_2_double(&a[8], &a[40], &a[136], &a[168], P);
	DFT_2_double(&a[9], &a[41], &a[137], &a[169], P);
	DFT_2_double(&a[10], &a[42], &a[138], &a[170], P);
	DFT_2_double(&a[11], &a[43], &a[139], &a[171], P);
	DFT_2_double(&a[12], &a[44], &a[140], &a[172], P);
	DFT_2_double(&a[13], &a[45], &a[141], &a[173], P);
	DFT_2_double(&a[14], &a[46], &a[142], &a[174], P);
	DFT_2_double(&a[15], &a[47], &a[143], &a[175], P);
	DFT_2_double(&a[16], &a[48], &a[144], &a[176], P);
	DFT_2_double(&a[17], &a[49], &a[145], &a[177], P);
	DFT_2_double(&a[18], &a[50], &a[146], &a[178], P);
	DFT_2_double(&a[19], &a[51], &a[147], &a[179], P);
	DFT_2_double(&a[20], &a[52], &a[148], &a[180], P);
	DFT_2_double(&a[21], &a[53], &a[149], &a[181], P);
	DFT_2_double(&a[22], &a[54], &a[150], &a[182], P);
	DFT_2_double(&a[23], &a[55], &a[151], &a[183], P);
	DFT_2_double(&a[24], &a[56], &a[152], &a[184], P);
	DFT_2_double(&a[25], &a[57], &a[153], &a[185], P);
	DFT_2_double(&a[26], &a[58], &a[154], &a[186], P);
	DFT_2_double(&a[27], &a[59], &a[155], &a[187], P);
	DFT_2_double(&a[28], &a[60], &a[156], &a[188], P);
	DFT_2_double(&a[29], &a[61], &a[157], &a[189], P);
	DFT_2_double(&a[30], &a[62], &a[158], &a[190], P);
	DFT_2_double(&a[31], &a[63], &a[159], &a[191], P);
	DFT_2_double(&a[64], &a[96], &a[192], &a[224], P);
	DFT_2_double(&a[65], &a[97], &a[193], &a[225], P);
	DFT_2_double(&a[66], &a[98], &a[194], &a[226], P);
	DFT_2_double(&a[67], &a[99], &a[195], &a[227], P);
	DFT_2_double(&a[68], &a[100], &a[196], &a[228], P);
	DFT_2_double(&a[69], &a[101], &a[197], &a[229], P);
	DFT_2_double(&a[70], &a[102], &a[198], &a[230], P);
	DFT_2_double(&a[71], &a[103], &a[199], &a[231], P);
	DFT_2_double(&a[72], &a[104], &a[200], &a[232], P);
	DFT_2_double(&a[73], &a[105], &a[201], &a[233], P);
	DFT_2_double(&a[74], &a[106], &a[202], &a[234], P);
	DFT_2_double(&a[75], &a[107], &a[203], &a[235], P);
	DFT_2_double(&a[76], &a[108], &a[204], &a[236], P);
	DFT_2_double(&a[77], &a[109], &a[205], &a[237], P);
	DFT_2_double(&a[78], &a[110], &a[206], &a[238], P);
	DFT_2_double(&a[79], &a[111], &a[207], &a[239], P);
	DFT_2_double(&a[80], &a[112], &a[208], &a[240], P);
	DFT_2_double(&a[81], &a[113], &a[209], &a[241], P);
	DFT_2_double(&a[82], &a[114], &a[210], &a[242], P);
	DFT_2_double(&a[83], &a[115], &a[211], &a[243], P);
	DFT_2_double(&a[84], &a[116], &a[212], &a[244], P);
	DFT_2_double(&a[85], &a[117], &a[213], &a[245], P);
	DFT_2_double(&a[86], &a[118], &a[214], &a[246], P);
	DFT_2_double(&a[87], &a[119], &a[215], &a[247], P);
	DFT_2_double(&a[88], &a[120], &a[216], &a[248], P);
	DFT_2_double(&a[89], &a[121], &a[217], &a[249], P);
	DFT_2_double(&a[90], &a[122], &a[218], &a[250], P);
	DFT_2_double(&a[91], &a[123], &a[219], &a[251], P);
	DFT_2_double(&a[92], &a[124], &a[220], &a[252], P);
	DFT_2_double(&a[93], &a[125], &a[221], &a[253], P);
	DFT_2_double(&a[94], &a[126], &a[222], &a[254], P);
	DFT_2_double(&a[95], &a[127], &a[223], &a[255], P);
	mult_ab_mod_p_ptr(&a[48], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[49], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[50], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[51], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[52], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[53], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[54], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[56], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[57], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[58], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[60], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[80], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[81], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[82], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[83], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[84], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[85], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[86], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[88], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[89], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[90], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[92], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[112], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[113], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[114], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[116], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[120], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[144], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[145], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[146], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[147], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[148], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[149], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[150], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[151], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[152], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[153], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[154], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[155], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[156], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[157], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[158], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[159], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[176], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[177], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[178], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[179], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[180], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[181], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[182], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[183], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[184], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[185], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[186], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[187], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[188], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[189], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[190], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[191], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[208], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[209], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[210], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[211], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[212], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[213], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[214], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[215], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[216], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[217], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[218], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[219], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[220], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[221], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[222], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[223], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[240], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[241], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[242], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[243], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[244], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[245], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[246], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[247], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[248], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[249], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[250], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[251], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[252], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[253], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[254], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[112], P);
	DFT_2_double(&a[0], &a[16], &a[128], &a[144], P);
	DFT_2_double(&a[1], &a[17], &a[129], &a[145], P);
	DFT_2_double(&a[2], &a[18], &a[130], &a[146], P);
	DFT_2_double(&a[3], &a[19], &a[131], &a[147], P);
	DFT_2_double(&a[4], &a[20], &a[132], &a[148], P);
	DFT_2_double(&a[5], &a[21], &a[133], &a[149], P);
	DFT_2_double(&a[6], &a[22], &a[134], &a[150], P);
	DFT_2_double(&a[7], &a[23], &a[135], &a[151], P);
	DFT_2_double(&a[8], &a[24], &a[136], &a[152], P);
	DFT_2_double(&a[9], &a[25], &a[137], &a[153], P);
	DFT_2_double(&a[10], &a[26], &a[138], &a[154], P);
	DFT_2_double(&a[11], &a[27], &a[139], &a[155], P);
	DFT_2_double(&a[12], &a[28], &a[140], &a[156], P);
	DFT_2_double(&a[13], &a[29], &a[141], &a[157], P);
	DFT_2_double(&a[14], &a[30], &a[142], &a[158], P);
	DFT_2_double(&a[15], &a[31], &a[143], &a[159], P);
	DFT_2_double(&a[32], &a[48], &a[160], &a[176], P);
	DFT_2_double(&a[33], &a[49], &a[161], &a[177], P);
	DFT_2_double(&a[34], &a[50], &a[162], &a[178], P);
	DFT_2_double(&a[35], &a[51], &a[163], &a[179], P);
	DFT_2_double(&a[36], &a[52], &a[164], &a[180], P);
	DFT_2_double(&a[37], &a[53], &a[165], &a[181], P);
	DFT_2_double(&a[38], &a[54], &a[166], &a[182], P);
	DFT_2_double(&a[39], &a[55], &a[167], &a[183], P);
	DFT_2_double(&a[40], &a[56], &a[168], &a[184], P);
	DFT_2_double(&a[41], &a[57], &a[169], &a[185], P);
	DFT_2_double(&a[42], &a[58], &a[170], &a[186], P);
	DFT_2_double(&a[43], &a[59], &a[171], &a[187], P);
	DFT_2_double(&a[44], &a[60], &a[172], &a[188], P);
	DFT_2_double(&a[45], &a[61], &a[173], &a[189], P);
	DFT_2_double(&a[46], &a[62], &a[174], &a[190], P);
	DFT_2_double(&a[47], &a[63], &a[175], &a[191], P);
	DFT_2_double(&a[64], &a[80], &a[192], &a[208], P);
	DFT_2_double(&a[65], &a[81], &a[193], &a[209], P);
	DFT_2_double(&a[66], &a[82], &a[194], &a[210], P);
	DFT_2_double(&a[67], &a[83], &a[195], &a[211], P);
	DFT_2_double(&a[68], &a[84], &a[196], &a[212], P);
	DFT_2_double(&a[69], &a[85], &a[197], &a[213], P);
	DFT_2_double(&a[70], &a[86], &a[198], &a[214], P);
	DFT_2_double(&a[71], &a[87], &a[199], &a[215], P);
	DFT_2_double(&a[72], &a[88], &a[200], &a[216], P);
	DFT_2_double(&a[73], &a[89], &a[201], &a[217], P);
	DFT_2_double(&a[74], &a[90], &a[202], &a[218], P);
	DFT_2_double(&a[75], &a[91], &a[203], &a[219], P);
	DFT_2_double(&a[76], &a[92], &a[204], &a[220], P);
	DFT_2_double(&a[77], &a[93], &a[205], &a[221], P);
	DFT_2_double(&a[78], &a[94], &a[206], &a[222], P);
	DFT_2_double(&a[79], &a[95], &a[207], &a[223], P);
	DFT_2_double(&a[96], &a[112], &a[224], &a[240], P);
	DFT_2_double(&a[97], &a[113], &a[225], &a[241], P);
	DFT_2_double(&a[98], &a[114], &a[226], &a[242], P);
	DFT_2_double(&a[99], &a[115], &a[227], &a[243], P);
	DFT_2_double(&a[100], &a[116], &a[228], &a[244], P);
	DFT_2_double(&a[101], &a[117], &a[229], &a[245], P);
	DFT_2_double(&a[102], &a[118], &a[230], &a[246], P);
	DFT_2_double(&a[103], &a[119], &a[231], &a[247], P);
	DFT_2_double(&a[104], &a[120], &a[232], &a[248], P);
	DFT_2_double(&a[105], &a[121], &a[233], &a[249], P);
	DFT_2_double(&a[106], &a[122], &a[234], &a[250], P);
	DFT_2_double(&a[107], &a[123], &a[235], &a[251], P);
	DFT_2_double(&a[108], &a[124], &a[236], &a[252], P);
	DFT_2_double(&a[109], &a[125], &a[237], &a[253], P);
	DFT_2_double(&a[110], &a[126], &a[238], &a[254], P);
	DFT_2_double(&a[111], &a[127], &a[239], &a[255], P);
	mult_ab_mod_p_ptr(&a[24], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[25], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[26], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[27], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[28], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[29], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[30], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[40], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[41], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[42], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[43], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[44], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[45], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[46], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[56], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[57], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[58], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[60], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[72], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[73], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[74], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[75], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[76], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[77], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[78], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[88], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[89], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[90], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[92], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[104], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[105], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[106], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[108], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[120], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[136], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[137], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[138], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[139], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[140], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[141], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[142], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[143], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[152], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[153], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[154], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[155], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[156], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[157], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[158], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[159], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[168], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[169], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[170], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[171], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[172], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[173], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[174], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[175], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[184], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[185], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[186], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[187], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[188], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[189], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[190], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[191], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[200], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[201], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[202], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[203], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[204], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[205], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[206], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[207], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[216], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[217], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[218], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[219], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[220], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[221], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[222], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[223], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[232], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[233], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[234], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[235], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[236], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[237], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[238], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[239], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[248], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[249], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[250], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[251], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[252], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[253], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[254], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[120], P);
	DFT_2_double(&a[0], &a[8], &a[128], &a[136], P);
	DFT_2_double(&a[1], &a[9], &a[129], &a[137], P);
	DFT_2_double(&a[2], &a[10], &a[130], &a[138], P);
	DFT_2_double(&a[3], &a[11], &a[131], &a[139], P);
	DFT_2_double(&a[4], &a[12], &a[132], &a[140], P);
	DFT_2_double(&a[5], &a[13], &a[133], &a[141], P);
	DFT_2_double(&a[6], &a[14], &a[134], &a[142], P);
	DFT_2_double(&a[7], &a[15], &a[135], &a[143], P);
	DFT_2_double(&a[16], &a[24], &a[144], &a[152], P);
	DFT_2_double(&a[17], &a[25], &a[145], &a[153], P);
	DFT_2_double(&a[18], &a[26], &a[146], &a[154], P);
	DFT_2_double(&a[19], &a[27], &a[147], &a[155], P);
	DFT_2_double(&a[20], &a[28], &a[148], &a[156], P);
	DFT_2_double(&a[21], &a[29], &a[149], &a[157], P);
	DFT_2_double(&a[22], &a[30], &a[150], &a[158], P);
	DFT_2_double(&a[23], &a[31], &a[151], &a[159], P);
	DFT_2_double(&a[32], &a[40], &a[160], &a[168], P);
	DFT_2_double(&a[33], &a[41], &a[161], &a[169], P);
	DFT_2_double(&a[34], &a[42], &a[162], &a[170], P);
	DFT_2_double(&a[35], &a[43], &a[163], &a[171], P);
	DFT_2_double(&a[36], &a[44], &a[164], &a[172], P);
	DFT_2_double(&a[37], &a[45], &a[165], &a[173], P);
	DFT_2_double(&a[38], &a[46], &a[166], &a[174], P);
	DFT_2_double(&a[39], &a[47], &a[167], &a[175], P);
	DFT_2_double(&a[48], &a[56], &a[176], &a[184], P);
	DFT_2_double(&a[49], &a[57], &a[177], &a[185], P);
	DFT_2_double(&a[50], &a[58], &a[178], &a[186], P);
	DFT_2_double(&a[51], &a[59], &a[179], &a[187], P);
	DFT_2_double(&a[52], &a[60], &a[180], &a[188], P);
	DFT_2_double(&a[53], &a[61], &a[181], &a[189], P);
	DFT_2_double(&a[54], &a[62], &a[182], &a[190], P);
	DFT_2_double(&a[55], &a[63], &a[183], &a[191], P);
	DFT_2_double(&a[64], &a[72], &a[192], &a[200], P);
	DFT_2_double(&a[65], &a[73], &a[193], &a[201], P);
	DFT_2_double(&a[66], &a[74], &a[194], &a[202], P);
	DFT_2_double(&a[67], &a[75], &a[195], &a[203], P);
	DFT_2_double(&a[68], &a[76], &a[196], &a[204], P);
	DFT_2_double(&a[69], &a[77], &a[197], &a[205], P);
	DFT_2_double(&a[70], &a[78], &a[198], &a[206], P);
	DFT_2_double(&a[71], &a[79], &a[199], &a[207], P);
	DFT_2_double(&a[80], &a[88], &a[208], &a[216], P);
	DFT_2_double(&a[81], &a[89], &a[209], &a[217], P);
	DFT_2_double(&a[82], &a[90], &a[210], &a[218], P);
	DFT_2_double(&a[83], &a[91], &a[211], &a[219], P);
	DFT_2_double(&a[84], &a[92], &a[212], &a[220], P);
	DFT_2_double(&a[85], &a[93], &a[213], &a[221], P);
	DFT_2_double(&a[86], &a[94], &a[214], &a[222], P);
	DFT_2_double(&a[87], &a[95], &a[215], &a[223], P);
	DFT_2_double(&a[96], &a[104], &a[224], &a[232], P);
	DFT_2_double(&a[97], &a[105], &a[225], &a[233], P);
	DFT_2_double(&a[98], &a[106], &a[226], &a[234], P);
	DFT_2_double(&a[99], &a[107], &a[227], &a[235], P);
	DFT_2_double(&a[100], &a[108], &a[228], &a[236], P);
	DFT_2_double(&a[101], &a[109], &a[229], &a[237], P);
	DFT_2_double(&a[102], &a[110], &a[230], &a[238], P);
	DFT_2_double(&a[103], &a[111], &a[231], &a[239], P);
	DFT_2_double(&a[112], &a[120], &a[240], &a[248], P);
	DFT_2_double(&a[113], &a[121], &a[241], &a[249], P);
	DFT_2_double(&a[114], &a[122], &a[242], &a[250], P);
	DFT_2_double(&a[115], &a[123], &a[243], &a[251], P);
	DFT_2_double(&a[116], &a[124], &a[244], &a[252], P);
	DFT_2_double(&a[117], &a[125], &a[245], &a[253], P);
	DFT_2_double(&a[118], &a[126], &a[246], &a[254], P);
	DFT_2_double(&a[119], &a[127], &a[247], &a[255], P);
	mult_ab_mod_p_ptr(&a[12], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[13], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[14], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[15], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[20], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[21], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[22], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[23], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[28], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[29], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[30], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[36], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[37], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[38], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[39], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[44], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[45], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[46], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[52], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[53], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[54], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[60], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[68], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[69], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[70], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[71], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[76], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[77], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[78], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[84], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[85], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[86], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[92], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[100], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[101], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[102], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[108], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[116], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[124], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[132], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[133], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[134], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[135], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[140], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[141], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[142], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[143], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[148], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[149], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[150], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[151], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[156], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[157], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[158], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[159], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[164], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[165], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[166], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[167], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[172], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[173], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[174], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[175], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[180], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[181], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[182], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[183], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[188], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[189], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[190], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[191], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[196], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[197], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[198], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[199], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[204], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[205], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[206], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[207], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[212], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[213], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[214], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[215], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[220], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[221], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[222], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[223], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[228], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[229], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[230], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[231], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[236], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[237], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[238], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[239], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[244], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[245], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[246], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[247], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[252], &base_case_omega_vec[124], P);
	mult_ab_mod_p_ptr(&a[253], &base_case_omega_vec[124], P);
	mult_ab_mod_p_ptr(&a[254], &base_case_omega_vec[124], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[124], P);
	DFT_2_double(&a[0], &a[4], &a[128], &a[132], P);
	DFT_2_double(&a[1], &a[5], &a[129], &a[133], P);
	DFT_2_double(&a[2], &a[6], &a[130], &a[134], P);
	DFT_2_double(&a[3], &a[7], &a[131], &a[135], P);
	DFT_2_double(&a[8], &a[12], &a[136], &a[140], P);
	DFT_2_double(&a[9], &a[13], &a[137], &a[141], P);
	DFT_2_double(&a[10], &a[14], &a[138], &a[142], P);
	DFT_2_double(&a[11], &a[15], &a[139], &a[143], P);
	DFT_2_double(&a[16], &a[20], &a[144], &a[148], P);
	DFT_2_double(&a[17], &a[21], &a[145], &a[149], P);
	DFT_2_double(&a[18], &a[22], &a[146], &a[150], P);
	DFT_2_double(&a[19], &a[23], &a[147], &a[151], P);
	DFT_2_double(&a[24], &a[28], &a[152], &a[156], P);
	DFT_2_double(&a[25], &a[29], &a[153], &a[157], P);
	DFT_2_double(&a[26], &a[30], &a[154], &a[158], P);
	DFT_2_double(&a[27], &a[31], &a[155], &a[159], P);
	DFT_2_double(&a[32], &a[36], &a[160], &a[164], P);
	DFT_2_double(&a[33], &a[37], &a[161], &a[165], P);
	DFT_2_double(&a[34], &a[38], &a[162], &a[166], P);
	DFT_2_double(&a[35], &a[39], &a[163], &a[167], P);
	DFT_2_double(&a[40], &a[44], &a[168], &a[172], P);
	DFT_2_double(&a[41], &a[45], &a[169], &a[173], P);
	DFT_2_double(&a[42], &a[46], &a[170], &a[174], P);
	DFT_2_double(&a[43], &a[47], &a[171], &a[175], P);
	DFT_2_double(&a[48], &a[52], &a[176], &a[180], P);
	DFT_2_double(&a[49], &a[53], &a[177], &a[181], P);
	DFT_2_double(&a[50], &a[54], &a[178], &a[182], P);
	DFT_2_double(&a[51], &a[55], &a[179], &a[183], P);
	DFT_2_double(&a[56], &a[60], &a[184], &a[188], P);
	DFT_2_double(&a[57], &a[61], &a[185], &a[189], P);
	DFT_2_double(&a[58], &a[62], &a[186], &a[190], P);
	DFT_2_double(&a[59], &a[63], &a[187], &a[191], P);
	DFT_2_double(&a[64], &a[68], &a[192], &a[196], P);
	DFT_2_double(&a[65], &a[69], &a[193], &a[197], P);
	DFT_2_double(&a[66], &a[70], &a[194], &a[198], P);
	DFT_2_double(&a[67], &a[71], &a[195], &a[199], P);
	DFT_2_double(&a[72], &a[76], &a[200], &a[204], P);
	DFT_2_double(&a[73], &a[77], &a[201], &a[205], P);
	DFT_2_double(&a[74], &a[78], &a[202], &a[206], P);
	DFT_2_double(&a[75], &a[79], &a[203], &a[207], P);
	DFT_2_double(&a[80], &a[84], &a[208], &a[212], P);
	DFT_2_double(&a[81], &a[85], &a[209], &a[213], P);
	DFT_2_double(&a[82], &a[86], &a[210], &a[214], P);
	DFT_2_double(&a[83], &a[87], &a[211], &a[215], P);
	DFT_2_double(&a[88], &a[92], &a[216], &a[220], P);
	DFT_2_double(&a[89], &a[93], &a[217], &a[221], P);
	DFT_2_double(&a[90], &a[94], &a[218], &a[222], P);
	DFT_2_double(&a[91], &a[95], &a[219], &a[223], P);
	DFT_2_double(&a[96], &a[100], &a[224], &a[228], P);
	DFT_2_double(&a[97], &a[101], &a[225], &a[229], P);
	DFT_2_double(&a[98], &a[102], &a[226], &a[230], P);
	DFT_2_double(&a[99], &a[103], &a[227], &a[231], P);
	DFT_2_double(&a[104], &a[108], &a[232], &a[236], P);
	DFT_2_double(&a[105], &a[109], &a[233], &a[237], P);
	DFT_2_double(&a[106], &a[110], &a[234], &a[238], P);
	DFT_2_double(&a[107], &a[111], &a[235], &a[239], P);
	DFT_2_double(&a[112], &a[116], &a[240], &a[244], P);
	DFT_2_double(&a[113], &a[117], &a[241], &a[245], P);
	DFT_2_double(&a[114], &a[118], &a[242], &a[246], P);
	DFT_2_double(&a[115], &a[119], &a[243], &a[247], P);
	DFT_2_double(&a[120], &a[124], &a[248], &a[252], P);
	DFT_2_double(&a[121], &a[125], &a[249], &a[253], P);
	DFT_2_double(&a[122], &a[126], &a[250], &a[254], P);
	DFT_2_double(&a[123], &a[127], &a[251], &a[255], P);
	mult_ab_mod_p_ptr(&a[6], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[7], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[10], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[11], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[14], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[15], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[18], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[19], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[22], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[23], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[26], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[27], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[30], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[34], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[35], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[38], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[39], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[42], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[43], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[46], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[50], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[51], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[54], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[58], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[62], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[66], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[67], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[70], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[71], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[74], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[75], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[78], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[82], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[83], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[86], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[90], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[94], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[98], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[99], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[102], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[106], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[110], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[114], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[118], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[122], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[126], &base_case_omega_vec[124], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[124], P);
	mult_ab_mod_p_ptr(&a[130], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr(&a[131], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr(&a[134], &base_case_omega_vec[66], P);
	mult_ab_mod_p_ptr(&a[135], &base_case_omega_vec[66], P);
	mult_ab_mod_p_ptr(&a[138], &base_case_omega_vec[34], P);
	mult_ab_mod_p_ptr(&a[139], &base_case_omega_vec[34], P);
	mult_ab_mod_p_ptr(&a[142], &base_case_omega_vec[98], P);
	mult_ab_mod_p_ptr(&a[143], &base_case_omega_vec[98], P);
	mult_ab_mod_p_ptr(&a[146], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr(&a[147], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr(&a[150], &base_case_omega_vec[82], P);
	mult_ab_mod_p_ptr(&a[151], &base_case_omega_vec[82], P);
	mult_ab_mod_p_ptr(&a[154], &base_case_omega_vec[50], P);
	mult_ab_mod_p_ptr(&a[155], &base_case_omega_vec[50], P);
	mult_ab_mod_p_ptr(&a[158], &base_case_omega_vec[114], P);
	mult_ab_mod_p_ptr(&a[159], &base_case_omega_vec[114], P);
	mult_ab_mod_p_ptr(&a[162], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr(&a[163], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr(&a[166], &base_case_omega_vec[74], P);
	mult_ab_mod_p_ptr(&a[167], &base_case_omega_vec[74], P);
	mult_ab_mod_p_ptr(&a[170], &base_case_omega_vec[42], P);
	mult_ab_mod_p_ptr(&a[171], &base_case_omega_vec[42], P);
	mult_ab_mod_p_ptr(&a[174], &base_case_omega_vec[106], P);
	mult_ab_mod_p_ptr(&a[175], &base_case_omega_vec[106], P);
	mult_ab_mod_p_ptr(&a[178], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr(&a[179], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr(&a[182], &base_case_omega_vec[90], P);
	mult_ab_mod_p_ptr(&a[183], &base_case_omega_vec[90], P);
	mult_ab_mod_p_ptr(&a[186], &base_case_omega_vec[58], P);
	mult_ab_mod_p_ptr(&a[187], &base_case_omega_vec[58], P);
	mult_ab_mod_p_ptr(&a[190], &base_case_omega_vec[122], P);
	mult_ab_mod_p_ptr(&a[191], &base_case_omega_vec[122], P);
	mult_ab_mod_p_ptr(&a[194], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr(&a[195], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr(&a[198], &base_case_omega_vec[70], P);
	mult_ab_mod_p_ptr(&a[199], &base_case_omega_vec[70], P);
	mult_ab_mod_p_ptr(&a[202], &base_case_omega_vec[38], P);
	mult_ab_mod_p_ptr(&a[203], &base_case_omega_vec[38], P);
	mult_ab_mod_p_ptr(&a[206], &base_case_omega_vec[102], P);
	mult_ab_mod_p_ptr(&a[207], &base_case_omega_vec[102], P);
	mult_ab_mod_p_ptr(&a[210], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr(&a[211], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr(&a[214], &base_case_omega_vec[86], P);
	mult_ab_mod_p_ptr(&a[215], &base_case_omega_vec[86], P);
	mult_ab_mod_p_ptr(&a[218], &base_case_omega_vec[54], P);
	mult_ab_mod_p_ptr(&a[219], &base_case_omega_vec[54], P);
	mult_ab_mod_p_ptr(&a[222], &base_case_omega_vec[118], P);
	mult_ab_mod_p_ptr(&a[223], &base_case_omega_vec[118], P);
	mult_ab_mod_p_ptr(&a[226], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr(&a[227], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr(&a[230], &base_case_omega_vec[78], P);
	mult_ab_mod_p_ptr(&a[231], &base_case_omega_vec[78], P);
	mult_ab_mod_p_ptr(&a[234], &base_case_omega_vec[46], P);
	mult_ab_mod_p_ptr(&a[235], &base_case_omega_vec[46], P);
	mult_ab_mod_p_ptr(&a[238], &base_case_omega_vec[110], P);
	mult_ab_mod_p_ptr(&a[239], &base_case_omega_vec[110], P);
	mult_ab_mod_p_ptr(&a[242], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr(&a[243], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr(&a[246], &base_case_omega_vec[94], P);
	mult_ab_mod_p_ptr(&a[247], &base_case_omega_vec[94], P);
	mult_ab_mod_p_ptr(&a[250], &base_case_omega_vec[62], P);
	mult_ab_mod_p_ptr(&a[251], &base_case_omega_vec[62], P);
	mult_ab_mod_p_ptr(&a[254], &base_case_omega_vec[126], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[126], P);
	DFT_2_double(&a[0], &a[2], &a[128], &a[130], P);
	DFT_2_double(&a[1], &a[3], &a[129], &a[131], P);
	DFT_2_double(&a[4], &a[6], &a[132], &a[134], P);
	DFT_2_double(&a[5], &a[7], &a[133], &a[135], P);
	DFT_2_double(&a[8], &a[10], &a[136], &a[138], P);
	DFT_2_double(&a[9], &a[11], &a[137], &a[139], P);
	DFT_2_double(&a[12], &a[14], &a[140], &a[142], P);
	DFT_2_double(&a[13], &a[15], &a[141], &a[143], P);
	DFT_2_double(&a[16], &a[18], &a[144], &a[146], P);
	DFT_2_double(&a[17], &a[19], &a[145], &a[147], P);
	DFT_2_double(&a[20], &a[22], &a[148], &a[150], P);
	DFT_2_double(&a[21], &a[23], &a[149], &a[151], P);
	DFT_2_double(&a[24], &a[26], &a[152], &a[154], P);
	DFT_2_double(&a[25], &a[27], &a[153], &a[155], P);
	DFT_2_double(&a[28], &a[30], &a[156], &a[158], P);
	DFT_2_double(&a[29], &a[31], &a[157], &a[159], P);
	DFT_2_double(&a[32], &a[34], &a[160], &a[162], P);
	DFT_2_double(&a[33], &a[35], &a[161], &a[163], P);
	DFT_2_double(&a[36], &a[38], &a[164], &a[166], P);
	DFT_2_double(&a[37], &a[39], &a[165], &a[167], P);
	DFT_2_double(&a[40], &a[42], &a[168], &a[170], P);
	DFT_2_double(&a[41], &a[43], &a[169], &a[171], P);
	DFT_2_double(&a[44], &a[46], &a[172], &a[174], P);
	DFT_2_double(&a[45], &a[47], &a[173], &a[175], P);
	DFT_2_double(&a[48], &a[50], &a[176], &a[178], P);
	DFT_2_double(&a[49], &a[51], &a[177], &a[179], P);
	DFT_2_double(&a[52], &a[54], &a[180], &a[182], P);
	DFT_2_double(&a[53], &a[55], &a[181], &a[183], P);
	DFT_2_double(&a[56], &a[58], &a[184], &a[186], P);
	DFT_2_double(&a[57], &a[59], &a[185], &a[187], P);
	DFT_2_double(&a[60], &a[62], &a[188], &a[190], P);
	DFT_2_double(&a[61], &a[63], &a[189], &a[191], P);
	DFT_2_double(&a[64], &a[66], &a[192], &a[194], P);
	DFT_2_double(&a[65], &a[67], &a[193], &a[195], P);
	DFT_2_double(&a[68], &a[70], &a[196], &a[198], P);
	DFT_2_double(&a[69], &a[71], &a[197], &a[199], P);
	DFT_2_double(&a[72], &a[74], &a[200], &a[202], P);
	DFT_2_double(&a[73], &a[75], &a[201], &a[203], P);
	DFT_2_double(&a[76], &a[78], &a[204], &a[206], P);
	DFT_2_double(&a[77], &a[79], &a[205], &a[207], P);
	DFT_2_double(&a[80], &a[82], &a[208], &a[210], P);
	DFT_2_double(&a[81], &a[83], &a[209], &a[211], P);
	DFT_2_double(&a[84], &a[86], &a[212], &a[214], P);
	DFT_2_double(&a[85], &a[87], &a[213], &a[215], P);
	DFT_2_double(&a[88], &a[90], &a[216], &a[218], P);
	DFT_2_double(&a[89], &a[91], &a[217], &a[219], P);
	DFT_2_double(&a[92], &a[94], &a[220], &a[222], P);
	DFT_2_double(&a[93], &a[95], &a[221], &a[223], P);
	DFT_2_double(&a[96], &a[98], &a[224], &a[226], P);
	DFT_2_double(&a[97], &a[99], &a[225], &a[227], P);
	DFT_2_double(&a[100], &a[102], &a[228], &a[230], P);
	DFT_2_double(&a[101], &a[103], &a[229], &a[231], P);
	DFT_2_double(&a[104], &a[106], &a[232], &a[234], P);
	DFT_2_double(&a[105], &a[107], &a[233], &a[235], P);
	DFT_2_double(&a[108], &a[110], &a[236], &a[238], P);
	DFT_2_double(&a[109], &a[111], &a[237], &a[239], P);
	DFT_2_double(&a[112], &a[114], &a[240], &a[242], P);
	DFT_2_double(&a[113], &a[115], &a[241], &a[243], P);
	DFT_2_double(&a[116], &a[118], &a[244], &a[246], P);
	DFT_2_double(&a[117], &a[119], &a[245], &a[247], P);
	DFT_2_double(&a[120], &a[122], &a[248], &a[250], P);
	DFT_2_double(&a[121], &a[123], &a[249], &a[251], P);
	DFT_2_double(&a[124], &a[126], &a[252], &a[254], P);
	DFT_2_double(&a[125], &a[127], &a[253], &a[255], P);
	mult_ab_mod_p_ptr(&a[3], &base_case_omega_vec[64], P);
	mult_ab_mod_p_ptr(&a[5], &base_case_omega_vec[32], P);
	mult_ab_mod_p_ptr(&a[7], &base_case_omega_vec[96], P);
	mult_ab_mod_p_ptr(&a[9], &base_case_omega_vec[16], P);
	mult_ab_mod_p_ptr(&a[11], &base_case_omega_vec[80], P);
	mult_ab_mod_p_ptr(&a[13], &base_case_omega_vec[48], P);
	mult_ab_mod_p_ptr(&a[15], &base_case_omega_vec[112], P);
	mult_ab_mod_p_ptr(&a[17], &base_case_omega_vec[8], P);
	mult_ab_mod_p_ptr(&a[19], &base_case_omega_vec[72], P);
	mult_ab_mod_p_ptr(&a[21], &base_case_omega_vec[40], P);
	mult_ab_mod_p_ptr(&a[23], &base_case_omega_vec[104], P);
	mult_ab_mod_p_ptr(&a[25], &base_case_omega_vec[24], P);
	mult_ab_mod_p_ptr(&a[27], &base_case_omega_vec[88], P);
	mult_ab_mod_p_ptr(&a[29], &base_case_omega_vec[56], P);
	mult_ab_mod_p_ptr(&a[31], &base_case_omega_vec[120], P);
	mult_ab_mod_p_ptr(&a[33], &base_case_omega_vec[4], P);
	mult_ab_mod_p_ptr(&a[35], &base_case_omega_vec[68], P);
	mult_ab_mod_p_ptr(&a[37], &base_case_omega_vec[36], P);
	mult_ab_mod_p_ptr(&a[39], &base_case_omega_vec[100], P);
	mult_ab_mod_p_ptr(&a[41], &base_case_omega_vec[20], P);
	mult_ab_mod_p_ptr(&a[43], &base_case_omega_vec[84], P);
	mult_ab_mod_p_ptr(&a[45], &base_case_omega_vec[52], P);
	mult_ab_mod_p_ptr(&a[47], &base_case_omega_vec[116], P);
	mult_ab_mod_p_ptr(&a[49], &base_case_omega_vec[12], P);
	mult_ab_mod_p_ptr(&a[51], &base_case_omega_vec[76], P);
	mult_ab_mod_p_ptr(&a[53], &base_case_omega_vec[44], P);
	mult_ab_mod_p_ptr(&a[55], &base_case_omega_vec[108], P);
	mult_ab_mod_p_ptr(&a[57], &base_case_omega_vec[28], P);
	mult_ab_mod_p_ptr(&a[59], &base_case_omega_vec[92], P);
	mult_ab_mod_p_ptr(&a[61], &base_case_omega_vec[60], P);
	mult_ab_mod_p_ptr(&a[63], &base_case_omega_vec[124], P);
	mult_ab_mod_p_ptr(&a[65], &base_case_omega_vec[2], P);
	mult_ab_mod_p_ptr(&a[67], &base_case_omega_vec[66], P);
	mult_ab_mod_p_ptr(&a[69], &base_case_omega_vec[34], P);
	mult_ab_mod_p_ptr(&a[71], &base_case_omega_vec[98], P);
	mult_ab_mod_p_ptr(&a[73], &base_case_omega_vec[18], P);
	mult_ab_mod_p_ptr(&a[75], &base_case_omega_vec[82], P);
	mult_ab_mod_p_ptr(&a[77], &base_case_omega_vec[50], P);
	mult_ab_mod_p_ptr(&a[79], &base_case_omega_vec[114], P);
	mult_ab_mod_p_ptr(&a[81], &base_case_omega_vec[10], P);
	mult_ab_mod_p_ptr(&a[83], &base_case_omega_vec[74], P);
	mult_ab_mod_p_ptr(&a[85], &base_case_omega_vec[42], P);
	mult_ab_mod_p_ptr(&a[87], &base_case_omega_vec[106], P);
	mult_ab_mod_p_ptr(&a[89], &base_case_omega_vec[26], P);
	mult_ab_mod_p_ptr(&a[91], &base_case_omega_vec[90], P);
	mult_ab_mod_p_ptr(&a[93], &base_case_omega_vec[58], P);
	mult_ab_mod_p_ptr(&a[95], &base_case_omega_vec[122], P);
	mult_ab_mod_p_ptr(&a[97], &base_case_omega_vec[6], P);
	mult_ab_mod_p_ptr(&a[99], &base_case_omega_vec[70], P);
	mult_ab_mod_p_ptr(&a[101], &base_case_omega_vec[38], P);
	mult_ab_mod_p_ptr(&a[103], &base_case_omega_vec[102], P);
	mult_ab_mod_p_ptr(&a[105], &base_case_omega_vec[22], P);
	mult_ab_mod_p_ptr(&a[107], &base_case_omega_vec[86], P);
	mult_ab_mod_p_ptr(&a[109], &base_case_omega_vec[54], P);
	mult_ab_mod_p_ptr(&a[111], &base_case_omega_vec[118], P);
	mult_ab_mod_p_ptr(&a[113], &base_case_omega_vec[14], P);
	mult_ab_mod_p_ptr(&a[115], &base_case_omega_vec[78], P);
	mult_ab_mod_p_ptr(&a[117], &base_case_omega_vec[46], P);
	mult_ab_mod_p_ptr(&a[119], &base_case_omega_vec[110], P);
	mult_ab_mod_p_ptr(&a[121], &base_case_omega_vec[30], P);
	mult_ab_mod_p_ptr(&a[123], &base_case_omega_vec[94], P);
	mult_ab_mod_p_ptr(&a[125], &base_case_omega_vec[62], P);
	mult_ab_mod_p_ptr(&a[127], &base_case_omega_vec[126], P);
	mult_ab_mod_p_ptr(&a[129], &base_case_omega_vec[1], P);
	mult_ab_mod_p_ptr(&a[131], &base_case_omega_vec[65], P);
	mult_ab_mod_p_ptr(&a[133], &base_case_omega_vec[33], P);
	mult_ab_mod_p_ptr(&a[135], &base_case_omega_vec[97], P);
	mult_ab_mod_p_ptr(&a[137], &base_case_omega_vec[17], P);
	mult_ab_mod_p_ptr(&a[139], &base_case_omega_vec[81], P);
	mult_ab_mod_p_ptr(&a[141], &base_case_omega_vec[49], P);
	mult_ab_mod_p_ptr(&a[143], &base_case_omega_vec[113], P);
	mult_ab_mod_p_ptr(&a[145], &base_case_omega_vec[9], P);
	mult_ab_mod_p_ptr(&a[147], &base_case_omega_vec[73], P);
	mult_ab_mod_p_ptr(&a[149], &base_case_omega_vec[41], P);
	mult_ab_mod_p_ptr(&a[151], &base_case_omega_vec[105], P);
	mult_ab_mod_p_ptr(&a[153], &base_case_omega_vec[25], P);
	mult_ab_mod_p_ptr(&a[155], &base_case_omega_vec[89], P);
	mult_ab_mod_p_ptr(&a[157], &base_case_omega_vec[57], P);
	mult_ab_mod_p_ptr(&a[159], &base_case_omega_vec[121], P);
	mult_ab_mod_p_ptr(&a[161], &base_case_omega_vec[5], P);
	mult_ab_mod_p_ptr(&a[163], &base_case_omega_vec[69], P);
	mult_ab_mod_p_ptr(&a[165], &base_case_omega_vec[37], P);
	mult_ab_mod_p_ptr(&a[167], &base_case_omega_vec[101], P);
	mult_ab_mod_p_ptr(&a[169], &base_case_omega_vec[21], P);
	mult_ab_mod_p_ptr(&a[171], &base_case_omega_vec[85], P);
	mult_ab_mod_p_ptr(&a[173], &base_case_omega_vec[53], P);
	mult_ab_mod_p_ptr(&a[175], &base_case_omega_vec[117], P);
	mult_ab_mod_p_ptr(&a[177], &base_case_omega_vec[13], P);
	mult_ab_mod_p_ptr(&a[179], &base_case_omega_vec[77], P);
	mult_ab_mod_p_ptr(&a[181], &base_case_omega_vec[45], P);
	mult_ab_mod_p_ptr(&a[183], &base_case_omega_vec[109], P);
	mult_ab_mod_p_ptr(&a[185], &base_case_omega_vec[29], P);
	mult_ab_mod_p_ptr(&a[187], &base_case_omega_vec[93], P);
	mult_ab_mod_p_ptr(&a[189], &base_case_omega_vec[61], P);
	mult_ab_mod_p_ptr(&a[191], &base_case_omega_vec[125], P);
	mult_ab_mod_p_ptr(&a[193], &base_case_omega_vec[3], P);
	mult_ab_mod_p_ptr(&a[195], &base_case_omega_vec[67], P);
	mult_ab_mod_p_ptr(&a[197], &base_case_omega_vec[35], P);
	mult_ab_mod_p_ptr(&a[199], &base_case_omega_vec[99], P);
	mult_ab_mod_p_ptr(&a[201], &base_case_omega_vec[19], P);
	mult_ab_mod_p_ptr(&a[203], &base_case_omega_vec[83], P);
	mult_ab_mod_p_ptr(&a[205], &base_case_omega_vec[51], P);
	mult_ab_mod_p_ptr(&a[207], &base_case_omega_vec[115], P);
	mult_ab_mod_p_ptr(&a[209], &base_case_omega_vec[11], P);
	mult_ab_mod_p_ptr(&a[211], &base_case_omega_vec[75], P);
	mult_ab_mod_p_ptr(&a[213], &base_case_omega_vec[43], P);
	mult_ab_mod_p_ptr(&a[215], &base_case_omega_vec[107], P);
	mult_ab_mod_p_ptr(&a[217], &base_case_omega_vec[27], P);
	mult_ab_mod_p_ptr(&a[219], &base_case_omega_vec[91], P);
	mult_ab_mod_p_ptr(&a[221], &base_case_omega_vec[59], P);
	mult_ab_mod_p_ptr(&a[223], &base_case_omega_vec[123], P);
	mult_ab_mod_p_ptr(&a[225], &base_case_omega_vec[7], P);
	mult_ab_mod_p_ptr(&a[227], &base_case_omega_vec[71], P);
	mult_ab_mod_p_ptr(&a[229], &base_case_omega_vec[39], P);
	mult_ab_mod_p_ptr(&a[231], &base_case_omega_vec[103], P);
	mult_ab_mod_p_ptr(&a[233], &base_case_omega_vec[23], P);
	mult_ab_mod_p_ptr(&a[235], &base_case_omega_vec[87], P);
	mult_ab_mod_p_ptr(&a[237], &base_case_omega_vec[55], P);
	mult_ab_mod_p_ptr(&a[239], &base_case_omega_vec[119], P);
	mult_ab_mod_p_ptr(&a[241], &base_case_omega_vec[15], P);
	mult_ab_mod_p_ptr(&a[243], &base_case_omega_vec[79], P);
	mult_ab_mod_p_ptr(&a[245], &base_case_omega_vec[47], P);
	mult_ab_mod_p_ptr(&a[247], &base_case_omega_vec[111], P);
	mult_ab_mod_p_ptr(&a[249], &base_case_omega_vec[31], P);
	mult_ab_mod_p_ptr(&a[251], &base_case_omega_vec[95], P);
	mult_ab_mod_p_ptr(&a[253], &base_case_omega_vec[63], P);
	mult_ab_mod_p_ptr(&a[255], &base_case_omega_vec[127], P);
	DFT_2_double(&a[0], &a[1], &a[128], &a[129], P);
	DFT_2_double(&a[2], &a[3], &a[130], &a[131], P);
	DFT_2_double(&a[4], &a[5], &a[132], &a[133], P);
	DFT_2_double(&a[6], &a[7], &a[134], &a[135], P);
	DFT_2_double(&a[8], &a[9], &a[136], &a[137], P);
	DFT_2_double(&a[10], &a[11], &a[138], &a[139], P);
	DFT_2_double(&a[12], &a[13], &a[140], &a[141], P);
	DFT_2_double(&a[14], &a[15], &a[142], &a[143], P);
	DFT_2_double(&a[16], &a[17], &a[144], &a[145], P);
	DFT_2_double(&a[18], &a[19], &a[146], &a[147], P);
	DFT_2_double(&a[20], &a[21], &a[148], &a[149], P);
	DFT_2_double(&a[22], &a[23], &a[150], &a[151], P);
	DFT_2_double(&a[24], &a[25], &a[152], &a[153], P);
	DFT_2_double(&a[26], &a[27], &a[154], &a[155], P);
	DFT_2_double(&a[28], &a[29], &a[156], &a[157], P);
	DFT_2_double(&a[30], &a[31], &a[158], &a[159], P);
	DFT_2_double(&a[32], &a[33], &a[160], &a[161], P);
	DFT_2_double(&a[34], &a[35], &a[162], &a[163], P);
	DFT_2_double(&a[36], &a[37], &a[164], &a[165], P);
	DFT_2_double(&a[38], &a[39], &a[166], &a[167], P);
	DFT_2_double(&a[40], &a[41], &a[168], &a[169], P);
	DFT_2_double(&a[42], &a[43], &a[170], &a[171], P);
	DFT_2_double(&a[44], &a[45], &a[172], &a[173], P);
	DFT_2_double(&a[46], &a[47], &a[174], &a[175], P);
	DFT_2_double(&a[48], &a[49], &a[176], &a[177], P);
	DFT_2_double(&a[50], &a[51], &a[178], &a[179], P);
	DFT_2_double(&a[52], &a[53], &a[180], &a[181], P);
	DFT_2_double(&a[54], &a[55], &a[182], &a[183], P);
	DFT_2_double(&a[56], &a[57], &a[184], &a[185], P);
	DFT_2_double(&a[58], &a[59], &a[186], &a[187], P);
	DFT_2_double(&a[60], &a[61], &a[188], &a[189], P);
	DFT_2_double(&a[62], &a[63], &a[190], &a[191], P);
	DFT_2_double(&a[64], &a[65], &a[192], &a[193], P);
	DFT_2_double(&a[66], &a[67], &a[194], &a[195], P);
	DFT_2_double(&a[68], &a[69], &a[196], &a[197], P);
	DFT_2_double(&a[70], &a[71], &a[198], &a[199], P);
	DFT_2_double(&a[72], &a[73], &a[200], &a[201], P);
	DFT_2_double(&a[74], &a[75], &a[202], &a[203], P);
	DFT_2_double(&a[76], &a[77], &a[204], &a[205], P);
	DFT_2_double(&a[78], &a[79], &a[206], &a[207], P);
	DFT_2_double(&a[80], &a[81], &a[208], &a[209], P);
	DFT_2_double(&a[82], &a[83], &a[210], &a[211], P);
	DFT_2_double(&a[84], &a[85], &a[212], &a[213], P);
	DFT_2_double(&a[86], &a[87], &a[214], &a[215], P);
	DFT_2_double(&a[88], &a[89], &a[216], &a[217], P);
	DFT_2_double(&a[90], &a[91], &a[218], &a[219], P);
	DFT_2_double(&a[92], &a[93], &a[220], &a[221], P);
	DFT_2_double(&a[94], &a[95], &a[222], &a[223], P);
	DFT_2_double(&a[96], &a[97], &a[224], &a[225], P);
	DFT_2_double(&a[98], &a[99], &a[226], &a[227], P);
	DFT_2_double(&a[100], &a[101], &a[228], &a[229], P);
	DFT_2_double(&a[102], &a[103], &a[230], &a[231], P);
	DFT_2_double(&a[104], &a[105], &a[232], &a[233], P);
	DFT_2_double(&a[106], &a[107], &a[234], &a[235], P);
	DFT_2_double(&a[108], &a[109], &a[236], &a[237], P);
	DFT_2_double(&a[110], &a[111], &a[238], &a[239], P);
	DFT_2_double(&a[112], &a[113], &a[240], &a[241], P);
	DFT_2_double(&a[114], &a[115], &a[242], &a[243], P);
	DFT_2_double(&a[116], &a[117], &a[244], &a[245], P);
	DFT_2_double(&a[118], &a[119], &a[246], &a[247], P);
	DFT_2_double(&a[120], &a[121], &a[248], &a[249], P);
	DFT_2_double(&a[122], &a[123], &a[250], &a[251], P);
	DFT_2_double(&a[124], &a[125], &a[252], &a[253], P);
	DFT_2_double(&a[126], &a[127], &a[254], &a[255], P);
	swap (&a[1], &a[128]);
	swap (&a[2], &a[64]);
	swap (&a[3], &a[192]);
	swap (&a[4], &a[32]);
	swap (&a[5], &a[160]);
	swap (&a[6], &a[96]);
	swap (&a[7], &a[224]);
	swap (&a[8], &a[16]);
	swap (&a[9], &a[144]);
	swap (&a[10], &a[80]);
	swap (&a[11], &a[208]);
	swap (&a[12], &a[48]);
	swap (&a[13], &a[176]);
	swap (&a[14], &a[112]);
	swap (&a[15], &a[240]);
	swap (&a[17], &a[136]);
	swap (&a[18], &a[72]);
	swap (&a[19], &a[200]);
	swap (&a[20], &a[40]);
	swap (&a[21], &a[168]);
	swap (&a[22], &a[104]);
	swap (&a[23], &a[232]);
	swap (&a[25], &a[152]);
	swap (&a[26], &a[88]);
	swap (&a[27], &a[216]);
	swap (&a[28], &a[56]);
	swap (&a[29], &a[184]);
	swap (&a[30], &a[120]);
	swap (&a[31], &a[248]);
	swap (&a[33], &a[132]);
	swap (&a[34], &a[68]);
	swap (&a[35], &a[196]);
	swap (&a[37], &a[164]);
	swap (&a[38], &a[100]);
	swap (&a[39], &a[228]);
	swap (&a[41], &a[148]);
	swap (&a[42], &a[84]);
	swap (&a[43], &a[212]);
	swap (&a[44], &a[52]);
	swap (&a[45], &a[180]);
	swap (&a[46], &a[116]);
	swap (&a[47], &a[244]);
	swap (&a[49], &a[140]);
	swap (&a[50], &a[76]);
	swap (&a[51], &a[204]);
	swap (&a[53], &a[172]);
	swap (&a[54], &a[108]);
	swap (&a[55], &a[236]);
	swap (&a[57], &a[156]);
	swap (&a[58], &a[92]);
	swap (&a[59], &a[220]);
	swap (&a[61], &a[188]);
	swap (&a[62], &a[124]);
	swap (&a[63], &a[252]);
	swap (&a[65], &a[130]);
	swap (&a[67], &a[194]);
	swap (&a[69], &a[162]);
	swap (&a[70], &a[98]);
	swap (&a[71], &a[226]);
	swap (&a[73], &a[146]);
	swap (&a[74], &a[82]);
	swap (&a[75], &a[210]);
	swap (&a[77], &a[178]);
	swap (&a[78], &a[114]);
	swap (&a[79], &a[242]);
	swap (&a[81], &a[138]);
	swap (&a[83], &a[202]);
	swap (&a[85], &a[170]);
	swap (&a[86], &a[106]);
	swap (&a[87], &a[234]);
	swap (&a[89], &a[154]);
	swap (&a[91], &a[218]);
	swap (&a[93], &a[186]);
	swap (&a[94], &a[122]);
	swap (&a[95], &a[250]);
	swap (&a[97], &a[134]);
	swap (&a[99], &a[198]);
	swap (&a[101], &a[166]);
	swap (&a[103], &a[230]);
	swap (&a[105], &a[150]);
	swap (&a[107], &a[214]);
	swap (&a[109], &a[182]);
	swap (&a[110], &a[118]);
	swap (&a[111], &a[246]);
	swap (&a[113], &a[142]);
	swap (&a[115], &a[206]);
	swap (&a[117], &a[174]);
	swap (&a[119], &a[238]);
	swap (&a[121], &a[158]);
	swap (&a[123], &a[222]);
	swap (&a[125], &a[190]);
	swap (&a[127], &a[254]);
	swap (&a[131], &a[193]);
	swap (&a[133], &a[161]);
	swap (&a[135], &a[225]);
	swap (&a[137], &a[145]);
	swap (&a[139], &a[209]);
	swap (&a[141], &a[177]);
	swap (&a[143], &a[241]);
	swap (&a[147], &a[201]);
	swap (&a[149], &a[169]);
	swap (&a[151], &a[233]);
	swap (&a[155], &a[217]);
	swap (&a[157], &a[185]);
	swap (&a[159], &a[249]);
	swap (&a[163], &a[197]);
	swap (&a[167], &a[229]);
	swap (&a[171], &a[213]);
	swap (&a[173], &a[181]);
	swap (&a[175], &a[245]);
	swap (&a[179], &a[205]);
	swap (&a[183], &a[237]);
	swap (&a[187], &a[221]);
	swap (&a[191], &a[253]);
	swap (&a[199], &a[227]);
	swap (&a[203], &a[211]);
	swap (&a[207], &a[243]);
	swap (&a[215], &a[235]);
	swap (&a[223], &a[251]);
	swap (&a[239], &a[247]);
      }
  }

///////////////////////////////////////////////////////////

void
base_L_K_K (usfixn64*x, int K)
{
  int n_bytes = (K * K) * sizeof(usfixn64);
  usfixn64 *t = (usfixn64*) malloc (n_bytes);

  int block_size =
      (K < SQURE_PERMUTATION_BLOCK_SIZE) ? K : SQURE_PERMUTATION_BLOCK_SIZE;
#pragma unroll
  for (int i = 0; i < K; i += block_size)
#pragma unroll
    for (int j = 0; j < K; j += block_size)
#pragma unroll
      for (int k = i; k < (i + block_size); k++)
#pragma unroll
	for (int l = j; l < (j + block_size); l++)
	  t[l * K + k] = x[k * K + l];

  memcpy (x, t, n_bytes);
  free (t);
}

///////////////////////////////////////////////////////////

//L_J_K = transpose {K}{J}
void
stride_permutation_small_elements_L_J_K (int n_permutations, usfixn64*x,
					 const int J, const int K,
					 usfixn64* aux_memory)
{
  int K_bytes = K * sizeof(usfixn64);
  int n_bytes = (J * K_bytes);
//  usfixn64 * y=(usfixn64*)malloc(n_bytes);
  usfixn64 * y = aux_memory;

//  if ((J&(K-1)))
//    {
//      printf("ERROR! FAILED!\n");
//	printf("K=%d , J=%d\n", K, J);
//    }

  int stride = J * K;
  int JK = J / K;
  int K2 = K * K;

//#pragma cilk_grainsize = (n/(cilk::current_worker_count))
////  FOR_LOOP (int t=0;t<n_permutations;t++)
////#pragma cilk_grainsize = (SMALL_PERMUTATION_GRAIN_FACTOR*(cilk::current_worker_count))
//#pragma cilk_grainsize = (JK/(cilk::current_worker_count))
//  FOR_LOOP(int j=0;j<JK;j++)
//    {
//      for(int i=0;i<K;i++)
//	memcpy(&y[j*K2+i*K+t*stride], &x[i*J+j*K+t*stride], K_bytes);
//      L_K_K(&y[j*K2+t*stride], K);
//    }

#pragma cilk_grainsize = ((n*JK)/(cilk::current_worker_count))
  FOR_LOOP (int tidx = 0; tidx < (n_permutations * JK);
      tidx++
    )
//#pragma cilk_grainsize = (SMALL_PERMUTATION_GRAIN_FACTOR*(cilk::current_worker_count))
//  FOR_LOOP(int j=0;j<JK;j++)
      {
	int j = tidx / n_permutations;
	int t = tidx / JK;
	for (int i = 0; i < K; i++)
	memcpy (&y[j * K2 + i * K + t * stride], &x[i * J + j * K + t * stride],
	K_bytes);
	base_L_K_K (&y[j * K2 + t * stride], K);
      }
//  memcpy(x, y, n_bytes);
//  free(y);
  }

///////////////////////////////////////////////////////////

//L_K_J = transpose {J}{K}
void
stride_permutation_small_elements_L_K_J (int n_permutations, usfixn64*x,
					 const int K, const int J,
					 usfixn64* aux_memory)
{
  int K_bytes = K * sizeof(usfixn64);
  int n_bytes = (J * K_bytes);
//  usfixn64 * y=(usfixn64*)malloc(n_bytes);
  usfixn64 * y = aux_memory;
//  if ((J&(K-1)))
//    {
//      printf("ERROR! FAILED!\n");
//      printf("K=%d , J=%d\n", K, J);
//    }

  int stride = J * K;
  int JK = J / K;
  int K2 = K * K;
//cilk_
  for (int t = 0; t < n_permutations; t++)
#pragma cilk_grainsize = (SMALL_PERMUTATION_GRAIN_FACTOR*(cilk::current_worker_count))
//cilk_
    for (int j = 0; j < JK; j++)
      {
	base_L_K_K (&x[j * K2 + t * stride], K);
	for (int i = 0; i < K; i++)
	  memcpy (&y[i * J + j * K + t * stride],
		  &x[i * K + j * K2 + t * stride], K_bytes);
      }

//  FOR_LOOP (int tidx=0;tidx<(n_permutations*JK);tidx++)
//  #pragma cilk_grainsize = (SMALL_PERMUTATION_GRAIN_FACTOR*(cilk::current_worker_count))
////    FOR_LOOP(int j=0;j<JK;j++)
//      {
//        L_K_K(&x[j*K2+t*stride], K);
//        for(int i=0;i<K;i++)
//  	memcpy(&y[i*J+j*K+t*stride], &x[i*K+j*K2+t*stride], K_bytes);
//      }

//  memcpy(x, y, n_bytes);
//  free(y);
}

///////////////////////////////////////////////////////////
//VERIFIED.
//compute L_{m}{n} for a vector of spf elements
//equivalently, computing transpose {n}{m}
void
stride_permutation_small_elements(usfixn64*x, int n_permutations, int m, int n,
				   usfixn64* aux_memory)
{
  const usfixn64 blocksize = MIN(MIN(m,n),SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE);
  const usfixn64 data_size = (n_permutations * m * n * sizeof(usfixn64));
  const usfixn64 stride = m * n;
  usfixn64* y = aux_memory;  //(usfixn64*) malloc (data_size);

  usfixn64 n_ix_blocks=(n/blocksize);
  usfixn64 n_jx_blocks=(m/blocksize);

  FOR_LOOP (usfixn64 t = 0; t < n_permutations;t++)
      {
	usfixn64 t_x_stride=t*stride;
	//#pragma cilk_grainsize = (SMALL_PERMUTATION_GRAIN_FACTOR*(cilk::current_worker_count))
//	for (int i = 0; i < n; i += blocksize)
//	FOR_LOOP (int ix = 0; ix < (n/blocksize); ix +=1)
	FOR_LOOP (usfixn64 ix = 0; ix < (n_ix_blocks); ix +=1)
	  {
	    usfixn64 i=(ix*blocksize);
	    //#pragma cilk_grainsize = (SMALL_PERMUTATION_GRAIN_FACTOR*(cilk::current_worker_count))
	    //cilk_
//	    for (int j = 0; j < m; j += blocksize)
//	    FOR_LOOP(int jx = 0; jx < m/blocksize; jx +=1)
	    FOR_LOOP(usfixn64 jx = 0; jx < n_jx_blocks; jx +=1)
	      {
		usfixn64 j=(jx*blocksize);
		usfixn64 *x_ptr=&x[t_x_stride];
		usfixn64 *y_ptr=&y[t_x_stride];
		// transpose the block beginning at [i,j]
//#pragma unroll
		for (usfixn64 k = i; k < i + blocksize; k++)
		  {
		    usfixn64 km=(k * m);
//#pragma unroll
		    for (usfixn64 l = j; l < j + blocksize; l++)
		      {
//			y[k + l * n + (t * stride)] = x[l + km + (t * stride)];
//			y[k + l * n + (t_x_stride)] = x[l + km + (t_x_stride)];
			y_ptr[k + l * n] = x_ptr[l + km];
		      }
		  }
	      }
	  }
      }

    //// THERE ARE A NUMBER OF SERIOUS BUGS HERE.
//  memcpy (x, y, data_size);
////  free (y);
//  timer_record_stop(&t_memcpy);
//  timer_get_elapsed_time(&t_memcpy,NULL, 1);
//  total_memcpy+=t_memcpy.elapsed_time;
//  timer_print_time(total_memcpy, "t_memcpy");
    /////////////////////////////////////
    //// alternative method
    /////////////////////////////////////
//  if (m >= n)      //Tranpose {n}{m}
//    {
//      int K = n;
//      int J = m;
//      stride_permutation_small_elements_L_J_K (n_permutations, x, J, K,
//					       aux_memory);
//    }
//  else //m<n;
//    {
//      int K = m;
//      int J = n;
//      stride_permutation_small_elements_L_K_J (n_permutations, x, K, J,
//					       aux_memory);
//    }
    /////////////////////////////////////
  }

///////////////////////////////////////////////////////////

//VERIFIED.
//precompute powers of omega up to (omega)^{K^b}
//to be used at level b+1
//usfixn64* precomputed_pow_omega_vec should be be already allocated!
//just pass a pointer.
//ToDo: this is a prefix computation for mult -> implement it that way and in parallel.
void
precompute_sequential_powers_of_omega_small_elements (
    usfixn64* precomputed_pow_omega_vec, int K, int b, usfixn64 omega,
    const montgomery_triple* P)
{

//  printf("-- [smallprimefield]precompute powers of omega-at-level-%d up to omega^{K^%d = %d} \n   ",b + 1, b, n);
//  for (int i = 1; i < (n/K-1)*(K-1); i++)
//    {
////      if (i != 0)
//	{
//	  current_twiddle_u64 = mult_ab_mod_p (current_twiddle_u64, omega, P);
//	}
//      precomputed_pow_omega_vec[i] = current_twiddle_u64;
////      printf ("[w^{%d}=%lu\n", i, precomputed_pow_omega_vec[i]);
//    }

  int J = pow (K, b - 1);
  usfixn64 current_twiddle_u64;
  //form D_{K}{J}
//  usfixn64 omega_K = convertIn_GLOBAL (1, P);
//  int idx = 0;
//  for (int i = 0; i < K; i++)
//    {
//      current_twiddle_u64 = convertIn_GLOBAL (1, P); //w_K
//      for (int j = 0; j < J; j++)
//	{
//	  if ((j != 0) && (i != 0))
//	    current_twiddle_u64 = mult_ab_mod_p (current_twiddle_u64, omega_K,
//						 P);
//	  precomputed_pow_omega_vec[idx] = current_twiddle_u64;
//	  //      printf ("[w^{%d}=%lu\n", i, precomputed_pow_omega_vec[i]);
//	  idx++;
//	}
//      omega_K = mult_ab_mod_p (omega_K, omega, P);
//    }

  usfixn64 *omega_K = (usfixn64*) malloc (K * sizeof(usfixn64));
  omega_K[0] = 1;
  convertIn_GLOBAL_ptr (&omega_K[0], P);
  for (int i = 1; i < K; i++)
    omega_K[i] = mult_ab_mod_p (omega_K[i - 1], omega, P);

  usfixn64 one = 1;
  convertIn_GLOBAL_ptr (&one, P);
// #pragma cilk_grainsize = ((cilk::current_worker_count))
//cilk_
  for (int i = 0; i < K; i++)
    {
      int idx = i * J;
      for (int j = 0; j < J; j++)
	{
	  precomputed_pow_omega_vec[idx] = one; //w_K
	  if ((j != 0) && (i != 0))
	    precomputed_pow_omega_vec[idx] = mult_ab_mod_p (
		precomputed_pow_omega_vec[idx - 1], omega_K[i], P);
//      printf ("[w^{%d}=%lu\n", i, precomputed_pow_omega_vec[i]);
	  idx++;
	}
    }
  free (omega_K);
}

///////////////////////////////////////////////////////////

//VERIFIED.
//computing twiddle T{m}{n} == D{n}{m}
//computing D_{K}{J=K^b}
//part of twiddle computation is done by quotient powers;
//the rest is done by a direct power.
void
twiddle_small_elements_working (usfixn64* vector, const int n_permutations,
				const int K, const int J,
				const usfixn64* precomputed_pow_omega_vec,
				const montgomery_triple* P)
{
//  printf ("*** twiddle_small_elements *** ");
  int vector_idx;
  int KJ = K * J;
  int block_size;

  if (K > TWIDDLE_CACHE_SIZE)
    {
      block_size = TWIDDLE_CACHE_SIZE;
    }
  else
    {
      block_size = K;
    }
//  printf ("twiddle cache=%d\n", block_size);

//  (TWIDDLE_CACHE_SIZE>K?TWIDDLE_CACHE_SIZE:K);
  int n_blocks = KJ / block_size;
  for (int v_idx = 0; v_idx < n_blocks; v_idx += 1)
    {
      int v = v_idx * block_size;
      // FOR_LOOP (int i = 0; i < n_permutations; // causes slowdown! (it was uncommented -r7629).
      //  printf ("FOR_LOOP... [i=%d]\n", i);
      for (int i = 0; i < n_permutations; i++)
	{
//#pragma unroll SMALL_FFT_TWIDDLE_UNROLL_FACTOR
//	    for (int vector_idx = v; vector_idx < v + block_size; vector_idx++)
	  for (int bi = 0; bi < block_size; bi++)
	    {
	      usfixn32 vector_idx = v + bi;
	      //	int i=(vector_idx)/(J);
	      //	int j=(vector_idx)&(J-1);
	      //	twiddle_idx = i * j;
	      //	vector_idx = (i * J + j);
	      mult_ab_mod_p_ptr (&vector[vector_idx + i * KJ],
				 &precomputed_pow_omega_vec[vector_idx], P);
	    }
	}
    }
}

///////////////////////////////////////////////////////////

//VERIFIED.
//computing twiddle T{m}{n} == D{n}{m}
//computing D_{K}{J=K^b}
//part of twiddle computation is done by quotient powers;
//the rest is done by a direct power.
void
twiddle_small_elements (usfixn64* vector, const int n_permutations, const int K,
			const int J, const usfixn64* precomputed_pow_omega_vec,
			const montgomery_triple* P)
{
//  printf ("*** twiddle_small_elements *** ");
  int vector_idx;
  int KJ = K * J;
//  #pragma cilk_grainsize = ((n_permutations*KJ)/(cilk::current_worker_count))
//  FOR_LOOP (int i = 0; i < n_permutations; i++)
//      {
//	usfixn64* vector_ptr = &vector[i * KJ];
//	FOR_LOOP(int vector_idx = 0; vector_idx < KJ; vector_idx++)
//	  {
//	    mult_ab_mod_p_ptr (&vector_ptr[vector_idx],
//	    &precomputed_pow_omega_vec[vector_idx], P);
//	  }
//      }

  const usfixn64 n_iterations=n_permutations*KJ;
  const usfixn64 KJ_mask=KJ-1;
#pragma cilk_grainsize = MAX(n_iterations/(cilk::current_worker_count), SMALL_FFT_OPTIMAL_PERMUTATION_BLOCKSIZE)
  FOR_LOOP (usfixn64 i = 0; i < n_iterations; i++)
      {
//	usfixn64 vector_idx=(i%KJ);
	usfixn64 vector_idx=(i&KJ_mask);
	mult_ab_mod_p_ptr (&vector[i],	&precomputed_pow_omega_vec[vector_idx], P);
      }

}

///////////////////////////////////////////////////////////

// the precomputation of powers of omega: VERIFIED.
void
precompute_pow_omega_for_all_levels (usfixn64**precomputed_pow_omega_vec_2D,
				     const int K, const int e,
				     const usfixn64 omega,
				     const montgomery_triple* P)
{

#if SMALL_FFT_STEPS_TIMING_ENABLED
  cpu_timer t_precompute;
#endif

  usfixn64 current_omega = omega;
  for (int b = e; b >= 2; b--)
    {

#if SMALL_FFT_STEPS_TIMING_ENABLED
      timer_record_start (&t_precompute);
#endif

#if VERBOSE>=2
      printf ("[precomputing pow_omega for level b=%d]\n", b);
#endif
      int J = pow (K, b);
      precomputed_pow_omega_vec_2D[b] = (usfixn64*) malloc (
	  J * sizeof(usfixn64));
//      precomputed_pow_omega_vec_2D[b] = (usfixn64*)aligned_alloc(K*sizeof(usfixn64),J * sizeof(usfixn64));

      //precompute sequential powers of omega at level=b and store it in 2D vector.
      precompute_sequential_powers_of_omega_small_elements (
	  precomputed_pow_omega_vec_2D[b], K, b, current_omega, P);

//      as we go to a lower level, w=(w^K) mod p
      current_omega = exponent_mod_p_SPF (current_omega, K, P);
//      printVectorToFile (precomputed_pow_omega_vec_2D[b], J, coefficient_size,
//			 "pow_omega");
#if SMALL_FFT_STEPS_TIMING_ENABLED
      timer_record_stop (&t_precompute);
      timer_get_elapsed_time (&t_precompute, NULL, 1);
#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
      timer_print_time(t_precompute.elapsed_time, "t_precompute");
#endif
#endif
    }
}

///////////////////////////////////////////////////////////
int find_equivalent_config(int *K_in, int *e_in)
{
  int K=*K_in;
  int e=*e_in;

  int K_log2=log2(K);
  int n2=(K_log2*e);

  int max_K2=7; //6 or 7 or 8 (for 64, 128, and 256, respectively.)
  int best_K2=K_log2;
  int best_e=e;

  for(int K2=K_log2+1;K2<=max_K2;K2++)
    {
      if (n2%(K2))
	continue;
      else
	{
	  best_K2=K2;
	  best_e=n2/K2;
	}
    }

  // printf("best_K=%d (2 ^ %d), best_e2=%d\n", 1<<best_K2, best_K2, best_e);
  K=(1<<best_K2);
  e=best_e;
  *K_in=K;
  *e_in=e;
  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////

void
DFT_general_small_elements (usfixn64* x, int K, int e, montgomery_triple * P,
			    usfixn64 omega)
{
  find_equivalent_config(&K, &e);
  //we assume that the highest level is e  and the lowest level is 2;
  usfixn64 N, J;

//  usfixn64 p_inv = -1; //INV_CONVOLUTION_PRIME_1;
//  usfixn64 p = 4179340454199820289;//CONVOLUTION_PRIME_1;
//  usfixn64 p_inv = 4179340454199820287;
//  usfixn64 p_R = 1878466934230121386;

//  montgomery_triple P;
//  P.p = p;
//  P.p_inv = p_inv;
//  P.R = p_R;

  int K_log2 = log2 (K);

  N = pow (K, e);
  J = pow (K, e - 1);

//  usfixn64 stride;
//  int stride_log2 = K_log2;

  //precompute sequential powers of omega.
  usfixn64** precomputed_pow_omega_vec_2D = (usfixn64**) malloc (
      (e + 1) * sizeof(usfixn64*));

  usfixn64 * stride_permutation_aux = (usfixn64*) malloc (N * sizeof(usfixn64));
  usfixn64 * vector = x;

  int vector_ptr_check_bit = 0;	//0 means vector is pointing to vector_in;

  //////////////////////////////////////
  //// step 1: stride permutations on the right-most
  //// side of the DFT computing L{K,J}
  //////////////////////////////////////

   cpu_timer t_dft_general;
#if SMALL_FFT_STEPS_TIMING_ENABLED
  float t_permutation_total = 0;
  float t_twiddle_total = 0;
  float t_dft_K_total = 0;
  cpu_timer t_permutation;
  cpu_timer t_twiddle;
  cpu_timer t_dft_K;
  cpu_timer t_precompute;
#endif

//   timer_record_start (&t_dft_general);
#if SMALL_FFT_STEPS_TIMING_ENABLED
  // timer_record_start (&t_precompute);
#endif
  precompute_pow_omega_for_all_levels (precomputed_pow_omega_vec_2D, K, e,
				       omega, P);
#if SMALL_FFT_STEPS_TIMING_ENABLED
  // timer_record_stop (&t_precompute);
  // timer_get_elapsed_time (&t_precompute, NULL, 1);
#endif
  int stride_permutation_switch = 0;
  for (int b = e; b >= 2; b--)
    {
      J = pow (K, b - 1);
      int n_permutations = N / (K * J);
      //printf("-- step 1: (%-12s) b=%d , K=%d , J=%d\n", "right perm", b, K,J);

#if SMALL_FFT_STEPS_TIMING_ENABLED
      timer_record_start (&t_permutation);
#endif
      ////////////////////////////////////////////
      //////STRIDE PERMUTATION
      ////////////////////////////////////////////
      stride_permutation_small_elements (vector, n_permutations, K, J,
					 stride_permutation_aux);
      swap_ptr (&vector, &stride_permutation_aux, &vector_ptr_check_bit);
      ////////////////////////////////////////////
#if SMALL_FFT_STEPS_TIMING_ENABLED
      timer_record_stop (&t_permutation);
      timer_get_elapsed_time (&t_permutation, NULL, 1);
      t_permutation_total += t_permutation.elapsed_time;
#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
      timer_print_time(t_permutation.elapsed_time, "t_permutation");
#endif
#endif
    }
//  timer_print_time (t_permutation.elapsed_time, "t_permutation");

//  printf("-- step 2: (%-12s) b=%d , K=%d , J=%d \n", "right DFT_K", 2, K, J);
  //////////////////////////////////////
  //// step 2: base-case DFT's at the lowest level
  //////////////////////////////////////

  J = pow (K, e - 1);

  usfixn64 * precomputed_pow_omega_for_base_case_vec;    //of size K
  precomputed_pow_omega_for_base_case_vec = (usfixn64*) malloc (
      K * sizeof(usfixn64));

  usfixn64 base_case_omega = exponent_mod_p_SPF (omega, J, P);
  precompute_pow_omega (precomputed_pow_omega_for_base_case_vec,
			base_case_omega, K, P);

  /////////////////////////////////////
  //// right hand side DFT_K
  /////////////////////////////////////
#if SMALL_FFT_STEPS_TIMING_ENABLED
    timer_record_start (&t_dft_K);
#endif

  if (K == 4)
    {
      DFT_4 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }
  else if (K == 8)
    {
      DFT_8 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }
  else if (K == 16)
    {
      DFT_16 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }
  else if (K == 32)
    {
      DFT_32 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }
  else if (K == 64)
    {
      DFT_64 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }
  else if (K == 128)
    {
      DFT_128 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }
  else if (K == 256)
    {
      DFT_256 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
    }

#if SMALL_FFT_STEPS_TIMING_ENABLED
    timer_record_stop (&t_dft_K);
    timer_get_elapsed_time (&t_dft_K, NULL, 1);
    timer_print_time (t_dft_K.elapsed_time, "t_dft_K_right");

  // timer_record_stop (&t_dft_K);
  // timer_get_elapsed_time (&t_dft_K, NULL, 1);
#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
  timer_print_time(t_dft_K.elapsed_time, "t_dft_K");
#endif
  // t_dft_K_total += t_dft_K.elapsed_time;
#endif
  //ToDO: automated verification of base-case DFT's: DONE.

  ////////////////////////////////////////////////////////////
  //// following printToFile is used for testing purposes
  //// printVectorToFile(vector, N, coefficient_size, "xData_after_base_cases");
  ////////////////////////////////////////////////////////////

//  cpu_timer t_twiddle;
  for (int b = 2; b <= e; b++)
    {

      ///////////////////
      ////step 3: twiddle SMALL_DFT_GRAIN_FACTOR multiplication + permutation in between + DFT_K on the left
      ////beginning at the level b=2;
      ///////////////////
      //// T_{J}{K}=D_{K}{J}
      ///////////////////
      J = pow (K, b - 1);

      int n_permutations = N / (K * J);

//      printf ("-- step 3: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
//	      "twiddle+perm", b, K, J, n_permutations);

#if SMALL_FFT_STEPS_TIMING_ENABLED
      timer_record_start (&t_twiddle);
#endif
      ////////////////////////////////////////////
      //////TWIDDLE FACTORS
      ////////////////////////////////////////////
      twiddle_small_elements (&vector[0], n_permutations, K, J,
			      precomputed_pow_omega_vec_2D[b], P);
      ////////////////////////////////////////////
#if SMALL_FFT_STEPS_TIMING_ENABLED
      timer_record_stop (&t_twiddle);
      timer_get_elapsed_time (&t_twiddle, NULL, 1);
//#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
//      timer_print_time(t_twiddle.elapsed_time, "t_twiddle");
//#endif
      t_twiddle_total += t_twiddle.elapsed_time;
#endif

      ///////////////////
      //// L_{J}{K}
      ///////////////////
#if SMALL_FFT_STEPS_TIMING_ENABLED
      // timer_record_start (&t_permutation);
#endif
      ////////////////////////////////////////////
      //////STRIDE PERMUTATION
      ////////////////////////////////////////////
      stride_permutation_small_elements (vector, n_permutations, J, K,
					 stride_permutation_aux);
      swap_ptr (&vector, &stride_permutation_aux, &vector_ptr_check_bit);
      ////////////////////////////////////////////
#if SMALL_FFT_STEPS_TIMING_ENABLED
      // timer_record_stop (&t_permutation);
      // timer_get_elapsed_time (&t_permutation, NULL, 1);
      // t_permutation_total += t_permutation.elapsed_time;
#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
      timer_print_time(t_permutation.elapsed_time, "t_permutation");
#endif
#endif
      ///////////////////
      ////step 4:
      ////computing DFT_K on the left.
      ///////////////////
      J = N / K;
      n_permutations = N / (K * J);
//		for (int j = 0; j < pow(K, e - 1); j++)
//		printf("-- step 4: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n",
//				"left DFT_K", b, K, J, n_permutations);

#if SMALL_FFT_STEPS_TIMING_ENABLED
      // timer_record_start (&t_dft_K);
#endif
      if (K == 4)
	{
	  DFT_4 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
      else if (K == 8)
	{
	  DFT_8 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
      else if (K == 16)
	{
	  DFT_16 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
      else if (K == 32)
	{
	  DFT_32 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
      else if (K == 64)
	{
	  DFT_64 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
      else if (K == 128)
	{
	  DFT_128 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
      else if (K == 256)
	{
	  DFT_256 (J, vector, precomputed_pow_omega_for_base_case_vec, P);
	}
#if SMALL_FFT_STEPS_TIMING_ENABLED
      // timer_record_stop (&t_dft_K);
      // timer_get_elapsed_time (&t_dft_K, NULL, 1);
      // t_dft_K_total += t_dft_K.elapsed_time;
#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
      // timer_print_time(t_dft_K.elapsed_time, "t_dft_K");
#endif
#endif

      ///////////////////
      //// step 5: left-most permutation L_{K,K^e-1} on the right
      ///////////////////

      J = pow (K, b - 1);

      n_permutations = N / (K * J);
      //printf("-- step 5: (%-12s) b=%d , K=%d , J=%d , n_permutations = %d \n", "left perm", b, K, J, n_permutations);

#if SMALL_FFT_STEPS_TIMING_ENABLED
      // timer_record_start (&t_permutation);
#endif
      ////////////////////////////////////////////
      //////STRIDE PERMUTATION
      ////////////////////////////////////////////
      stride_permutation_small_elements (vector, n_permutations, K, J,
					 stride_permutation_aux);
      swap_ptr (&vector, &stride_permutation_aux, &vector_ptr_check_bit);
      ////////////////////////////////////////////
#if SMALL_FFT_STEPS_TIMING_ENABLED
      // timer_record_stop (&t_permutation);
      // timer_get_elapsed_time (&t_permutation, NULL, 1);
      // t_permutation_total += t_permutation.elapsed_time;
#if SMALL_FFT_STEPS_TIMING_ENABLED>=2
      timer_print_time(t_permutation.elapsed_time, "t_permutation");
#endif
#endif
    }

  ////////!!!!CRITCAL_REGION!!!!///////
  if (vector_ptr_check_bit == 1)
    {
      swap_ptr (&vector, &stride_permutation_aux, &vector_ptr_check_bit);
      memcpy (vector, stride_permutation_aux, N * sizeof(usfixn64));
    }
  ////////!!!!END_OF_CRITCAL_REGION!!!!///////
//   timer_record_stop (&t_dft_general);
//   timer_get_elapsed_time (&t_dft_general, "dft_general", 1);
#if SMALL_FFT_STEPS_TIMING_ENABLED
  // timer_print_time_percentage (t_precompute.elapsed_time, "t_precompute_total",
  // 		       t_dft_general.elapsed_time);
  // timer_print_time_percentage (t_permutation_total, "t_permutation_total",
  // 		       t_dft_general.elapsed_time);
  // timer_print_time_percentage (t_twiddle_total, "t_twiddle_total",
  // 		       t_dft_general.elapsed_time);
  // timer_print_time_percentage (t_dft_K_total, "t_dft_K_total",
  // 		       t_dft_general.elapsed_time);
#endif
  // timer_print_time_percentage (t_dft_general.elapsed_time, "t_dft_general",
  // 		       t_dft_general.elapsed_time);

  // printf ("%s", shortline);

#if SMALL_FFT_STEPS_TIMING_ENABLED
  timer_print_time (t_twiddle_total, "t_twiddle");
#endif

  for (int i = 2; i < (e + 1); i++)
    free (precomputed_pow_omega_vec_2D[i]);
  free (precomputed_pow_omega_vec_2D);
  free (precomputed_pow_omega_for_base_case_vec);
  free (stride_permutation_aux);
}

///////////////////////////////////////////////////////////

int
convolution_mult_small_elements (usfixn64* x_data, usfixn64* y_data, int n,
				 const montgomery_triple *P)
{
//#pragma cilk_grainsize = (CONVOLUTION_GRAIN_FACTOR*(cilk::current_worker_count))
#pragma cilk_grainsize = (n/(cilk::current_worker_count))
  FOR_LOOP(int i = 0; i < n; i++)
    {
      x_data[i] = mult_ab_mod_p (x_data[i], y_data[i], P);
    }
#if SMALL_FFT_STEPS_TIMING_ENABLED
  // timer_record_stop (&t_convolution_mult);
  // timer_get_elapsed_time (&t_convolution_mult, "t_convolution_mult", 1);
  // printf ("%s", shortline);
#endif
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////
int
test_cilk (int argc, char ** argv)
{
  int n_log2 = 10;

  if (argc > 1)
    n_log2 = atoi (argv[1]);
  usfixn32 n = (1 << n_log2);
  usfixn32* x = (usfixn32*) malloc (n * sizeof(usfixn32));
  usfixn32* y = (usfixn32*) malloc (n * sizeof(usfixn32));
  usfixn32* z = (usfixn32*) malloc (n * sizeof(usfixn32));

  srand (time (NULL));
  usfixn32 rv = rand ();
  for (int i = 0; i < n; i++)
    {
      x[i] = rv;
      rv--;
      y[i] = 2 * rv;
    }

  cpu_timer t_cilk;
  timer_record_start (&t_cilk);

  FOR_LOOP(int i=0;i < n;i++ )
	{
	  z[i]=x[i]+y[i];
	}

      timer_record_stop (&t_cilk);
      timer_get_elapsed_time (&t_cilk, "t_cilk", 1);

      usfixn32 idx = rand () % n / 2;
      idx += 5;
      printf ("z[%d]=%u\n", idx, z[idx]);

      free (x);
      free (y);
      free (z);
      return EXIT_SUCCESS;
    }
