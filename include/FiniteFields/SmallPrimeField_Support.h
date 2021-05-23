#ifndef _SMALLPRIMEFIELD_SUPPORT_H_
#define _SMALLPRIMEFIELD_SUPPORT_H_

#ifdef __cplusplus
    extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
typedef struct Prime_64{
   long long int prime;
   unsigned long long int prime_inv;
   long long int rsquare;
} Prime_ptr;

typedef struct Prime_32{
   int prime;
   unsigned int prime_inv;
   int rsquare;
} Prime_ptr32;

//for primes that fit in "long long int"
long long int smallprimefield_convert_in(long long int a, const Prime_ptr* Pptr);
long long int smallprimefield_convert_out(long long int a, const Prime_ptr* Pptr);
long long int smallprimefield_add(long long int a, long long int b, const Prime_ptr* Pptr);
long long int smallprimefield_sub( long long int a, long long int b, const Prime_ptr* Pptr);
long long int smallprimefield_mul(long long int a, long long int b, const Prime_ptr* Pptr);
long long int smallprimefield_inv(long long int a, const Prime_ptr* Pptr);
long long int smallprimefield_div(long long int a, long long int b, const Prime_ptr* Pptr);
long long int smallprimefield_PrimitiveRootofUnity(long long int n,const Prime_ptr* Pptr);
long long int smallprimefield_exp(long long int a, long long int e,const Prime_ptr* Pptr);
Prime_ptr* smallprimefield_get_prime_constants(long long int prime);

static long long int smallprimefield_AddMod(long long int a, long long int b, long long int prime){
  a = a + b;
  a -= prime;
  a += (a >> 63) & prime;
  return a;
}

static long long int smallprimefield_SubMod(long long int a, long long int b, long long int prime){
   a = a - b;
   a += (a >> 63) & prime;
   return a;
 }

 //compute a* b / R mod prime
static inline long long int smallprimefield_Mul_INLINE(long long int a, long long int b, long long int pinverse, long long int prime){
  __asm__(
      "mulq %2\n\t"
      "movq %%rax,%%rsi\n\t"
      "movq %%rdx,%%rdi\n\t"
      "imulq %3,%%rax\n\t"
      "mulq %4\n\t"
      "add %%rsi,%%rax\n\t"
      "adc %%rdi,%%rdx\n\t"
      "subq %4,%%rdx\n\t"
      "mov %%rdx,%%rax\n\t"
      "sar $63,%%rax\n\t"
      "andq %4,%%rax\n\t"
      "addq %%rax,%%rdx\n\t"
      : "=&d" (a)
      : "a"(a),"rm"(b),"b"((unsigned long long int)pinverse),"c"(prime)
      :"rsi","rdi");
  return a;
}

static inline long long int smallprimefield_Mul_INLINE2(long long int a,long long int b, long long int pinverse, long long int prime){
	asm("mulq %2\n\t"
		"movq %%rax,%%rsi\n\t"
		"movq %%rdx,%%rdi\n\t"
		"imulq %3,%%rax\n\t"
		"mulq %4\n\t"
		"add %%rsi,%%rax\n\t"
		"adc %%rdi,%%rdx\n\t"
		"subq %4,%%rdx\n\t"
		"mov %%rdx,%%rax\n\t"
		"sar $63,%%rax\n\t"
		"andq %4,%%rax\n\t"
		"addq %%rax,%%rdx\n\t"
		"nop\n\t"
		"nop\n\t"
		"nop\n\t"
		: "=d" (a)
		: "a"(a),"d"(b),"b"((long long int) pinverse),"c"((long long int) prime)
		:"rsi","rdi");
	return a;
}


//for primes that fit in "int"
int smallprimefield_convert_in32(int a, const Prime_ptr32* Pptr);
int smallprimefield_convert_out32(int a, const Prime_ptr32* Pptr);
int smallprimefield_add32(int a, int b, const Prime_ptr32* Pptr);
int smallprimefield_sub32(int a, int b, const Prime_ptr32* Pptr);
int smallprimefield_mul32(int a, int b, const Prime_ptr32* Pptr);
int smallprimefield_inv32(int a, const Prime_ptr32* Pptr);
int smallprimefield_div32(int a, int b, const Prime_ptr32* Pptr);
Prime_ptr32* smallprimefield_get_prime_constants32(int prime);
int smallprimefield_PrimitiveRootofUnity32(int n,const Prime_ptr32* Pptr);
int smallprimefield_exp32(int a, int e,const Prime_ptr32* Pptr);

#ifdef __cplusplus
}
#endif

#endif
