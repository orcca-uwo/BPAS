/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __AS_h
#define __AS_h 


#include "Types.h"
#include "CONSTANTS.h"


extern sfixn BASE;
extern sfixn BASE_1;
extern sfixn BASEHALF;



#ifdef LINUXINTEL64

/**
 * MulHiLoUnsigned: 
 * @h: pointer to a sfixn "a".
 * @l: pointer to a sfixn "b".
 *
 * Computing sfixn multiplication and return the high word and low word.
 *           I.e. *h = HighPart(ab), *l = LowPart(ab)
 * Return value: void.
 **/


// If Intel 32 Linux machine use this...

static inline void
MulHiLoUnsigned(sfixn *h, sfixn *l)
{
  //printf("in h= %ld, l= %ld \n",*h,*l);
  __asm__ ("mulq %3" : "=a" (*l), "=d" (* h) : "%0" (* h), "rm" (* l));
  //printf("out h= %ld, l= %ld \n",*h,*l);
}

static inline void
MulAndAddHiLoUnsigned(sfixn *h, sfixn *l,sfixn b,const sfixn p)
{
  //printf("in h= %ld, l= %ld \n",*h,*l);
  __asm__ ("mulq %3\n\t"
		"add %5,%%rax\n\t"
		"adc %4,%%rdx\n\t"
		: "=d" (* h),"=a"(*l)
		: "a" (b), "r" (p),"rm"(*h),"rm"(*l));
  //printf("out h= %ld, l= %ld \n",*h,*l);
}

static inline void AddHiLoUnsigned(sfixn*h, sfixn *l, sfixn a,sfixn b){
  __asm__ ("add %2,%0\n\t"
		"adc %3,%1\n\t"
		:"=r"(*l),"=r"(*h)
		:"r"(b),"r"(a));
}
static inline sfixn MontMultModSpe(sfixn h,sfixn l, const sfixn inv,const sfixn prime){
	asm("mulq %2\n\t"
		"movq %%rax,%%rsi\n\t"
		"movq %%rdx,%%rdi\n\t"
		"imulq %3,%%rax\n\t"
		"mulq %4\n\t"
		"add %%rsi,%%rax\n\t"
		"adc %%rdi,%%rdx\n\t"
		"subq %4,%%rdx\n\t"
		"mov %%rdx,%%rax\n\t"
		"sar %%cl,%%rax\n\t"
		"andq %4,%%rax\n\t"
		"addq %%rax,%%rdx\n\t"
		: "=d" (h)
		: "a"(h),"rm"(l),"r"(inv),"r"(prime),"c"(BASE_1)
		:"rsi","rdi");
	return h;
}
#else


// generic C code.


// the following signature has been modified by Xin 
//  On Nov.19.2008. to fix the g++ compilation (sfixn<->usfixn) problem.
static void
//MulHiLoUnsigned(usfixn *h, usfixn *l)
MulHiLoUnsigned(sfixn *h, sfixn *l)
{
  ulongfixnum prod;
  prod=(ulongfixnum)(usfixn)(*h) * (ulongfixnum)(usfixn)(*l);
  //  printf("a*b=%ld\n", prod);
  

/* hack to fix old GCC on Sparc 32-bit */
/* (big endian)  prod = [hi 32, lo 32] */
#if SOLARIS64
  *h=  *((usfixn *)&prod);
#else
  *h=  (sfixn)(((ulongfixnum)prod)>>BASE);
#endif
  *l=  (sfixn)(usfixn)prod;

  // printf("*h=%ld\n", *h);
  // printf("*l=%ld\n", *l);
}


#endif

#endif
