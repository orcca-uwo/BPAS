/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __inlineFuncs_h
#define __lnlineFuncs_h 

//#include "CONSTANTS.h"
#include "AS.h"

#include "emmintrin.h"

extern sfixn BASE;
extern sfixn BASE_1;
extern sfixn BASEHALF;

//static inline sfixn
//MontMulMod_OPT2_AS_GENE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr){
//  sfixn q2=pPtr->c_sft, c;
//
//  // b has been pre-leftShifted Base_Rpow bits.
//  MulHiLoUnsigned(&a, &b);
//  // now a keeps M1 quo R.  
//  // q2=c<<{Npow}.
//  MulHiLoUnsigned(&b, &q2);
//  // now b keeps M3 quo R;
//  a-=b;
//  // load c;
//  c=pPtr->c;
//  // be careful the sign of q2.
//  q2=(usfixn)q2>>pPtr->Base_Npow;
//  q2*=c;
//  // now q2 keeps M2 quo R;
//  // load P.
//  c=pPtr->P;
//  a+=q2;
//  // compute (M1+M2-M3) quo R.
//  // the result is in range  -p<a<2p.
//  // we need following 3 lines to reduce the error.
//  a += (a >> BASE_1) & c;
//  a -= c;
//  a += (a >> BASE_1) & c;
//  return a;
//}
#if AMD
static inline
void unrolled8MontMul4root16(sfixn* input1, sfixn root1,sfixn root2, sfixn root3,sfixn root4, MONTP_OPT2_AS_GENE * pPtr){
	unsigned long temp1[8]__attribute__((aligned(16)));
	sfixn prime[2]__attribute__((aligned(16))) = {pPtr->P,pPtr->P};
	asm ("movq (%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"movq %%rdx,(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,(%9)\n\t"

			"movq 8(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"movq %%rdx,8(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,8(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,8(%9)\n\t"

			"movq 128(%%rsi),%%rax\n\t"
			"mulq %3\n\t"
			"movq %%rdx,128(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,128(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,16(%9)\n\t"

			"movq 136(%%rsi),%%rax\n\t"
			"mulq %4\n\t"
			"movq %%rdx,136(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,136(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,24(%9)\n\t"

			"movq 512(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"movq %%rdx,512(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,512(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,32(%9)\n\t"

			"movq 520(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"movq %%rdx,520(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,520(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,40(%9)\n\t"

			"movq 640(%%rsi),%%rax\n\t"
			"mulq %3\n\t"
			"movq %%rdx,640(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,640(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,48(%9)\n\t"

			"movq 648(%%rsi),%%rax\n\t"
			"mulq %4\n\t"
			"movq %%rdx,648(%%rsi)\n\t"
			"mulq %5\n\t"
			"subq %%rdx,648(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,56(%9)\n\t"
//--------------------FIRST MULS DONE------------------------//
//----------------Quotients stored in xmm0-3-----------------//
//--------------------SECOND MULS DONE------------------------//
//----------------Quotients stored in xmm5-7-----------------//

			"movq (%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,(%9)\n\t"
			"movq 8(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,8(%9)\n\t"
			"movq 16(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,16(%9)\n\t"
			"movq 24(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,24(%9)\n\t"
			"movq 32(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,32(%9)\n\t"
			"movq 40(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,40(%9)\n\t"
			"movq 48(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,48(%9)\n\t"
			"movq 56(%9),%%rdx\n\t"
			"imulq %6,%%rdx\n\t"
			"movq %%rdx,56(%9)\n\t"
//------------------LAST MULS DONE----------------------------//
//------------------Rem stored in r8-r15-----------------------//
			"movdqa (%%rsi),%%xmm0\n\t"			
			"movdqa 128(%%rsi),%%xmm1\n\t"			
			"movdqa 512(%%rsi),%%xmm2\n\t"			
			"movdqa 640(%%rsi),%%xmm3\n\t"			

			"paddq (%9),%%xmm0\n\t"			
			"paddq 16(%9),%%xmm1\n\t"			
			"paddq 32(%9),%%xmm2\n\t"			
			"paddq 48(%9),%%xmm3\n\t"			
//-------------------Adjust in the right range----------------//

			//Get the prime//
			"movdqa (%8),%%xmm6\n\t"

			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"
//----------------Adjusted one time--------------------//
//----------------Now subtracting prime----------------//
			"psubq %%xmm6,%%xmm0\n\t"
			"psubq %%xmm6,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm6,%%xmm3\n\t"

//-------------------Adjust in the right range----------------//
			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"

			"movdqa %%xmm0,(%%rsi)\n\t"
			"movdqa %%xmm1,128(%%rsi)\n\t"
			"movdqa %%xmm2,512(%%rsi)\n\t"
			"movdqa %%xmm3,640(%%rsi)\n\t"
			:
			: "S" (input1), "g" (root1),"g"(root2),"g"(root3),"g"(root4),"g"(pPtr->c_sft),"g"(pPtr->c),"c"((unsigned char)pPtr->Base_Npow),"r"(prime),"r"(temp1)
			:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","rax","rdx","memory");
	return;
}
static inline
void unrolled8MontMul2root16(sfixn* input1, sfixn root1,sfixn root2, MONTP_OPT2_AS_GENE * pPtr){
	unsigned long temp1[8]__attribute__((aligned(16)));
	sfixn prime[2]__attribute__((aligned(16)))= {pPtr->P,pPtr->P};
	asm ("movq (%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"movq %%rdx,(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,(%7)\n\t"

			"movq 8(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"movq %%rdx,8(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,8(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,8(%7)\n\t"

			"movq 256(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"movq %%rdx,256(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,256(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,16(%7)\n\t"

			"movq 264(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"movq %%rdx,264(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,264(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,24(%7)\n\t"

			"movq 512(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"movq %%rdx,512(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,512(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,32(%7)\n\t"

			"movq 520(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"movq %%rdx,520(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,520(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,40(%7)\n\t"

			"movq 768(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"movq %%rdx,768(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,768(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,48(%7)\n\t"

			"movq 776(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"movq %%rdx,776(%%rsi)\n\t"
			"mulq %3\n\t"
			"subq %%rdx,776(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,56(%7)\n\t"
//--------------------FIRST MULS DONE------------------------//
//--------------------SECOND MULS DONE------------------------//
			"movq (%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,(%7)\n\t"
			"movq 8(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,8(%7)\n\t"
			"movq 16(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,16(%7)\n\t"
			"movq 24(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,24(%7)\n\t"
			"movq 32(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,32(%7)\n\t"
			"movq 40(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,40(%7)\n\t"
			"movq 48(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,48(%7)\n\t"
			"movq 56(%7),%%rdx\n\t"
			"imulq %4,%%rdx\n\t"
			"movq %%rdx,56(%7)\n\t"
//--------------//----LAST MULS DONE----------------------------//
			"movdqa (%%rsi),%%xmm0\n\t"			
			"movdqa 256(%%rsi),%%xmm1\n\t"			
			"movdqa 512(%%rsi),%%xmm2\n\t"			
			"movdqa 768(%%rsi),%%xmm3\n\t"			

			"paddq (%7),%%xmm0\n\t"			
			"paddq 16(%7),%%xmm1\n\t"			
			"paddq 32(%7),%%xmm2\n\t"			
			"paddq 48(%7),%%xmm3\n\t"			
//-------------------Adjust in the right range----------------//

			//Get the prime//
			"movdqa (%6),%%xmm6\n\t"

			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"
//----------------Adjusted one time--------------------//
//----------------Now subtracting prime----------------//
			"psubq %%xmm6,%%xmm0\n\t"
			"psubq %%xmm6,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm6,%%xmm3\n\t"

//-------------------Adjust in the right range----------------//
			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"

			"movdqa %%xmm0,(%%rsi)\n\t"
			"movdqa %%xmm1,256(%%rsi)\n\t"
			"movdqa %%xmm2,512(%%rsi)\n\t"
			"movdqa %%xmm3,768(%%rsi)\n\t"
			:
			: "S" (input1), "g" (root1),"g"(root2),"g"(pPtr->c_sft),"g"(pPtr->c),"c"((unsigned char)pPtr->Base_Npow),"r"(prime),"r"(temp1)
			:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","rax","rdx","memory");
	return;
}
static inline
void unrolled8MontMul(sfixn* input1, sfixn* input2, MONTP_OPT2_AS_GENE * pPtr){
	unsigned long temp1[8]__attribute__((aligned(16)));
	sfixn prime[2]__attribute__((aligned(16)))= {pPtr->P,pPtr->P};

	asm ("movq (%%rsi),%%rax\n\t"
			"mulq (%%rdi)\n\t"
			"movq %%rdx,(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,(%6)\n\t"

			"movq 8(%%rsi),%%rax\n\t"
			"mulq 8(%%rdi)\n\t"
			"movq %%rdx,8(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,8(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,8(%6)\n\t"

			"movq 16(%%rsi),%%rax\n\t"
			"mulq 16(%%rdi)\n\t"
			"movq %%rdx,16(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,16(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,16(%6)\n\t"

			"movq 24(%%rsi),%%rax\n\t"
			"mulq 24(%%rdi)\n\t"
			"movq %%rdx,24(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,24(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,24(%6)\n\t"

			"movq 32(%%rsi),%%rax\n\t"
			"mulq 32(%%rdi)\n\t"
			"movq %%rdx,32(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,32(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,32(%6)\n\t"

			"movq 40(%%rsi),%%rax\n\t"
			"mulq 40(%%rdi)\n\t"
			"movq %%rdx,40(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,40(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,40(%6)\n\t"

			"movq 48(%%rsi),%%rax\n\t"
			"mulq 48(%%rdi)\n\t"
			"movq %%rdx,48(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,48(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,48(%6)\n\t"

			"movq 56(%%rsi),%%rax\n\t"
			"mulq 56(%%rdi)\n\t"
			"movq %%rdx,56(%%rsi)\n\t"
			"mulq %2\n\t"
			"subq %%rdx,56(%%rsi)\n\t"
			"shr %%cl,%%rax\n\t"
			"movq %%rax,56(%6)\n\t"
//--------------------FI//RST MULS DONE------------------------//
//----------------Quotie//nts stored in xmm0-3-----------------//
//--------------------SE//COND MULS DONE------------------------//
//----------------Quotie//nts stored in xmm5-7-----------------//


			"movq (%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,(%6)\n\t"
			"movq 8(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,8(%6)\n\t"
			"movq 16(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,16(%6)\n\t"
			"movq 24(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,24(%6)\n\t"
			"movq 32(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,32(%6)\n\t"
			"movq 40(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,40(%6)\n\t"
			"movq 48(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,48(%6)\n\t"
			"movq 56(%6),%%rdx\n\t"
			"imulq %3,%%rdx\n\t"
			"movq %%rdx,56(%6)\n\t"
//------------------LASTULS DONE----------------------------//
			"movdqa (%%rsi),%%xmm0\n\t"			
			"movdqa 16(%%rsi),%%xmm1\n\t"			
			"movdqa 32(%%rsi),%%xmm2\n\t"			
			"movdqa 48(%%rsi),%%xmm3\n\t"			

			"paddq (%6),%%xmm0\n\t"			
			"paddq 16(%6),%%xmm1\n\t"			
			"paddq 32(%6),%%xmm2\n\t"			
			"paddq 48(%6),%%xmm3\n\t"			

//-------------------Adjust in the right range----------------//

			//Get the prime//
			"movdqa (%5),%%xmm6\n\t"

			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"
//----------------Adjusted one time--------------------//
//----------------Now subtracting prime----------------//
			"psubq %%xmm6,%%xmm0\n\t"
			"psubq %%xmm6,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm6,%%xmm3\n\t"

//-------------------Adjust in the right range----------------//
			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"

			"movdqa %%xmm0,(%%rsi)\n\t"
			"movdqa %%xmm1,16(%%rsi)\n\t"
			"movdqa %%xmm2,32(%%rsi)\n\t"
			"movdqa %%xmm3,48(%%rsi)\n\t"

			:
			: "S" (input1), "D" (input2),"g"(pPtr->c_sft),"r"(pPtr->c),"c"((unsigned char)pPtr->Base_Npow),"r"(prime),"r"(temp1)
			:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","rax","rdx","memory");
	return;
}
#endif
#if XEON5600
static inline
void unrolled8MontMul4root16(sfixn* input1, sfixn root1,sfixn root2, sfixn root3,sfixn root4, MONTP_OPT2_AS_GENE * pPtr){

	asm ("movq (%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"pinsrq $0,%%rdx,%%xmm0\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r8\n\t"
			"pinsrq $0,%%rdx,%%xmm4\n\t"

			"movq 8(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"pinsrq $1,%%rdx,%%xmm0\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r9\n\t"
			"pinsrq $1,%%rdx,%%xmm4\n\t"

			"movq 128(%%rsi),%%rax\n\t"
			"mulq %3\n\t"
			"pinsrq $0,%%rdx,%%xmm1\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r10\n\t"
			"pinsrq $0,%%rdx,%%xmm5\n\t"

			"movq 136(%%rsi),%%rax\n\t"
			"mulq %4\n\t"
			"pinsrq $1,%%rdx,%%xmm1\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r11\n\t"
			"pinsrq $1,%%rdx,%%xmm5\n\t"

			"movq 512(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"pinsrq $0,%%rdx,%%xmm2\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r12\n\t"
			"pinsrq $0,%%rdx,%%xmm6\n\t"

			"movq 520(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"pinsrq $1,%%rdx,%%xmm2\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r13\n\t"
			"pinsrq $1,%%rdx,%%xmm6\n\t"

			"movq 640(%%rsi),%%rax\n\t"
			"mulq %3\n\t"
			"pinsrq $0,%%rdx,%%xmm3\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r14\n\t"
			"pinsrq $0,%%rdx,%%xmm7\n\t"

			"movq 648(%%rsi),%%rax\n\t"
			"mulq %4\n\t"
			"pinsrq $1,%%rdx,%%xmm3\n\t"
			"mulq %5\n\t"
			"movq %%rax,%%r15\n\t"
			"pinsrq $1,%%rdx,%%xmm7\n\t"
//--------------------FIRST MULS DONE------------------------//
//----------------Quotients stored in xmm0-3-----------------//
//--------------------SECOND MULS DONE------------------------//
//----------------Quotients stored in xmm5-7-----------------//

			"shr %%cl,%%r8\n\t"
			"shr %%cl,%%r9\n\t"
			"shr %%cl,%%r10\n\t"
			"shr %%cl,%%r11\n\t"
			"shr %%cl,%%r12\n\t"
			"shr %%cl,%%r13\n\t"
			"shr %%cl,%%r14\n\t"
			"shr %%cl,%%r15\n\t"
			
			"psubq %%xmm4,%%xmm0\n\t"
			"psubq %%xmm5,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm7,%%xmm3\n\t"

			"imulq %6,%%r8\n\t"
			"imulq %6,%%r10\n\t"
			"imulq %6,%%r12\n\t"
			"imulq %6,%%r14\n\t"
			"imulq %6,%%r9\n\t"
			"imulq %6,%%r11\n\t"
			"imulq %6,%%r13\n\t"
			"imulq %6,%%r15\n\t"
//------------------LAST MULS DONE----------------------------//
//------------------Rem stored in r8-r15-----------------------//
			"pinsrq $0,%%r8,%%xmm4\n\t"
			"pinsrq $0,%%r10,%%xmm5\n\t"
			"pinsrq $0,%%r12,%%xmm6\n\t"
			"pinsrq $0,%%r14,%%xmm7\n\t"
			"pinsrq $1,%%r9,%%xmm4\n\t"
			"pinsrq $1,%%r11,%%xmm5\n\t"
			"pinsrq $1,%%r13,%%xmm6\n\t"
			"pinsrq $1,%%r15,%%xmm7\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"
			"paddq %%xmm6,%%xmm2\n\t"
			"paddq %%xmm7,%%xmm3\n\t"
//-------------------Adjust in the right range----------------//

			//Get the prime//
			"pinsrq $0,%8,%%xmm6\n\t"
			"pshufd $68,%%xmm6,%%xmm6\n\t"

			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"
//----------------Adjusted one time--------------------//
//----------------Now subtracting prime----------------//
			"psubq %%xmm6,%%xmm0\n\t"
			"psubq %%xmm6,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm6,%%xmm3\n\t"

//-------------------Adjust in the right range----------------//
			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"

			"movdqa %%xmm0,(%%rsi)\n\t"
			"movdqa %%xmm1,128(%%rsi)\n\t"
			"movdqa %%xmm2,512(%%rsi)\n\t"
			"movdqa %%xmm3,640(%%rsi)\n\t"
			:
			: "S" (input1), "g" (root1),"g"(root2),"g"(root3),"g"(root4),"g"(pPtr->c_sft),"g"(pPtr->c),"c"((unsigned char)pPtr->Base_Npow),"g"(pPtr->P)
			:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","rax","rdx","r8","r9","r10","r11","r12","r13","r14","r15");
	return;
}
static inline
void unrolled8MontMul2root16(sfixn* input1, sfixn root1,sfixn root2, MONTP_OPT2_AS_GENE * pPtr){

	asm ("movq (%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"pinsrq $0,%%rdx,%%xmm0\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r8\n\t"
			"pinsrq $0,%%rdx,%%xmm4\n\t"

			"movq 8(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"pinsrq $1,%%rdx,%%xmm0\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r9\n\t"
			"pinsrq $1,%%rdx,%%xmm4\n\t"

			"movq 256(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"pinsrq $0,%%rdx,%%xmm1\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r10\n\t"
			"pinsrq $0,%%rdx,%%xmm5\n\t"

			"movq 264(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"pinsrq $1,%%rdx,%%xmm1\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r11\n\t"
			"pinsrq $1,%%rdx,%%xmm5\n\t"

			"movq 512(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"pinsrq $0,%%rdx,%%xmm2\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r12\n\t"
			"pinsrq $0,%%rdx,%%xmm6\n\t"

			"movq 520(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"pinsrq $1,%%rdx,%%xmm2\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r13\n\t"
			"pinsrq $1,%%rdx,%%xmm6\n\t"

			"movq 768(%%rsi),%%rax\n\t"
			"mulq %1\n\t"
			"pinsrq $0,%%rdx,%%xmm3\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r14\n\t"
			"pinsrq $0,%%rdx,%%xmm7\n\t"

			"movq 776(%%rsi),%%rax\n\t"
			"mulq %2\n\t"
			"pinsrq $1,%%rdx,%%xmm3\n\t"
			"mulq %3\n\t"
			"movq %%rax,%%r15\n\t"
			"pinsrq $1,%%rdx,%%xmm7\n\t"
//--------------------FIRST MULS DONE------------------------//
//----------------Quotients stored in xmm0-3-----------------//
//--------------------SECOND MULS DONE------------------------//
//----------------Quotients stored in xmm5-7-----------------//

			"shr %%cl,%%r8\n\t"
			"shr %%cl,%%r9\n\t"
			"shr %%cl,%%r10\n\t"
			"shr %%cl,%%r11\n\t"
			"shr %%cl,%%r12\n\t"
			"shr %%cl,%%r13\n\t"
			"shr %%cl,%%r14\n\t"
			"shr %%cl,%%r15\n\t"
			
			"psubq %%xmm4,%%xmm0\n\t"
			"psubq %%xmm5,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm7,%%xmm3\n\t"

			"imulq %4,%%r8\n\t"
			"imulq %4,%%r10\n\t"
			"imulq %4,%%r12\n\t"
			"imulq %4,%%r14\n\t"
			"imulq %4,%%r9\n\t"
			"imulq %4,%%r11\n\t"
			"imulq %4,%%r13\n\t"
			"imulq %4,%%r15\n\t"
//------------------LAST MULS DONE----------------------------//
//------------------Rem stored in r8-r15-----------------------//
			"pinsrq $0,%%r8,%%xmm4\n\t"
			"pinsrq $0,%%r10,%%xmm5\n\t"
			"pinsrq $0,%%r12,%%xmm6\n\t"
			"pinsrq $0,%%r14,%%xmm7\n\t"
			"pinsrq $1,%%r9,%%xmm4\n\t"
			"pinsrq $1,%%r11,%%xmm5\n\t"
			"pinsrq $1,%%r13,%%xmm6\n\t"
			"pinsrq $1,%%r15,%%xmm7\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"
			"paddq %%xmm6,%%xmm2\n\t"
			"paddq %%xmm7,%%xmm3\n\t"
//-------------------Adjust in the right range----------------//

			//Get the prime//
			"pinsrq $0,%6,%%xmm6\n\t"
			"pshufd $68,%%xmm6,%%xmm6\n\t"

			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"
//----------------Adjusted one time--------------------//
//----------------Now subtracting prime----------------//
			"psubq %%xmm6,%%xmm0\n\t"
			"psubq %%xmm6,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm6,%%xmm3\n\t"

//-------------------Adjust in the right range----------------//
			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"

			"movdqa %%xmm0,(%%rsi)\n\t"
			"movdqa %%xmm1,256(%%rsi)\n\t"
			"movdqa %%xmm2,512(%%rsi)\n\t"
			"movdqa %%xmm3,768(%%rsi)\n\t"
			:
			: "S" (input1), "g" (root1),"g"(root2),"g"(pPtr->c_sft),"g"(pPtr->c),"c"((unsigned char)pPtr->Base_Npow),"g"(pPtr->P)
			:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","rax","rdx","r8","r9","r10","r11","r12","r13","r14","r15");
	return;
}
static inline
void unrolled8MontMul(sfixn* input1, sfixn* input2, MONTP_OPT2_AS_GENE * pPtr){

	asm ("movq (%%rsi),%%rax\n\t"
			"mulq (%%rdi)\n\t"
			"pinsrq $0,%%rdx,%%xmm0\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r8\n\t"
			"pinsrq $0,%%rdx,%%xmm4\n\t"

			"movq 8(%%rsi),%%rax\n\t"
			"mulq 8(%%rdi)\n\t"
			"pinsrq $1,%%rdx,%%xmm0\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r9\n\t"
			"pinsrq $1,%%rdx,%%xmm4\n\t"

			"movq 16(%%rsi),%%rax\n\t"
			"mulq 16(%%rdi)\n\t"
			"pinsrq $0,%%rdx,%%xmm1\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r10\n\t"
			"pinsrq $0,%%rdx,%%xmm5\n\t"

			"movq 24(%%rsi),%%rax\n\t"
			"mulq 24(%%rdi)\n\t"
			"pinsrq $1,%%rdx,%%xmm1\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r11\n\t"
			"pinsrq $1,%%rdx,%%xmm5\n\t"

			"movq 32(%%rsi),%%rax\n\t"
			"mulq 32(%%rdi)\n\t"
			"pinsrq $0,%%rdx,%%xmm2\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r12\n\t"
			"pinsrq $0,%%rdx,%%xmm6\n\t"

			"movq 40(%%rsi),%%rax\n\t"
			"mulq 40(%%rdi)\n\t"
			"pinsrq $1,%%rdx,%%xmm2\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r13\n\t"
			"pinsrq $1,%%rdx,%%xmm6\n\t"

			"movq 48(%%rsi),%%rax\n\t"
			"mulq 48(%%rdi)\n\t"
			"pinsrq $0,%%rdx,%%xmm3\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r14\n\t"
			"pinsrq $0,%%rdx,%%xmm7\n\t"

			"movq 56(%%rsi),%%rax\n\t"
			"mulq 56(%%rdi)\n\t"
			"pinsrq $1,%%rdx,%%xmm3\n\t"
			"mulq %2\n\t"
			"movq %%rax,%%r15\n\t"
			"pinsrq $1,%%rdx,%%xmm7\n\t"
//--------------------FIRST MULS DONE------------------------//
//----------------Quotients stored in xmm0-3-----------------//
//--------------------SECOND MULS DONE------------------------//
//----------------Quotients stored in xmm5-7-----------------//

			"shr %%cl,%%r8\n\t"
			"shr %%cl,%%r9\n\t"
			"shr %%cl,%%r10\n\t"
			"shr %%cl,%%r11\n\t"
			"shr %%cl,%%r12\n\t"
			"shr %%cl,%%r13\n\t"
			"shr %%cl,%%r14\n\t"
			"shr %%cl,%%r15\n\t"
			
			"psubq %%xmm4,%%xmm0\n\t"
			"psubq %%xmm5,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm7,%%xmm3\n\t"

			"imulq %3,%%r8\n\t"
			"imulq %3,%%r10\n\t"
			"imulq %3,%%r12\n\t"
			"imulq %3,%%r14\n\t"
			"imulq %3,%%r9\n\t"
			"imulq %3,%%r11\n\t"
			"imulq %3,%%r13\n\t"
			"imulq %3,%%r15\n\t"
//------------------LAST MULS DONE----------------------------//
//------------------Rem stored in r8-r15-----------------------//
			"pinsrq $0,%%r8,%%xmm4\n\t"
			"pinsrq $0,%%r10,%%xmm5\n\t"
			"pinsrq $0,%%r12,%%xmm6\n\t"
			"pinsrq $0,%%r14,%%xmm7\n\t"
			"pinsrq $1,%%r9,%%xmm4\n\t"
			"pinsrq $1,%%r11,%%xmm5\n\t"
			"pinsrq $1,%%r13,%%xmm6\n\t"
			"pinsrq $1,%%r15,%%xmm7\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"
			"paddq %%xmm6,%%xmm2\n\t"
			"paddq %%xmm7,%%xmm3\n\t"
//-------------------Adjust in the right range----------------//

			//Get the prime//
			"pinsrq $0,%5,%%xmm6\n\t"
			"pshufd $68,%%xmm6,%%xmm6\n\t"

			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"
//----------------Adjusted one time--------------------//
//----------------Now subtracting prime----------------//
			"psubq %%xmm6,%%xmm0\n\t"
			"psubq %%xmm6,%%xmm1\n\t"
			"psubq %%xmm6,%%xmm2\n\t"
			"psubq %%xmm6,%%xmm3\n\t"

//-------------------Adjust in the right range----------------//
			"movdqa %%xmm0,%%xmm4\n\t"
			"movdqa %%xmm1,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm0\n\t"
			"paddq %%xmm5,%%xmm1\n\t"

			"movdqa %%xmm2,%%xmm4\n\t"
			"movdqa %%xmm3,%%xmm5\n\t"

			"psrad $31,%%xmm4\n\t"
			"psrad $31,%%xmm5\n\t"

			"pshufd $245,%%xmm4,%%xmm4\n\t"
			"pshufd $245,%%xmm5,%%xmm5\n\t"

			"pand %%xmm6,%%xmm4\n\t"
			"pand %%xmm6,%%xmm5\n\t"

			"paddq %%xmm4,%%xmm2\n\t"
			"paddq %%xmm5,%%xmm3\n\t"

			"movdqa %%xmm0,(%%rsi)\n\t"
			"movdqa %%xmm1,16(%%rsi)\n\t"
			"movdqa %%xmm2,32(%%rsi)\n\t"
			"movdqa %%xmm3,48(%%rsi)\n\t"
			:
			: "S" (input1), "D" (input2),"g"(pPtr->c_sft),"g"(pPtr->c),"c"((unsigned char)pPtr->Base_Npow),"g"(pPtr->P)
			:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","rax","rdx","r8","r9","r10","r11","r12","r13","r14","r15");
	return;
}
#endif
static inline
void  unrolled8AddSubSSEMod(sfixn* a,sfixn* b,const sfixn* prime){
	asm("movdqa (%0),%%xmm1\n\t"
		"movdqa (%1),%%xmm2\n\t"
		"movdqa 16(%0),%%xmm5\n\t"
		"movdqa 16(%1),%%xmm6\n\t"
		"movdqa %2,%%xmm4\n\t"
		"movdqa %%xmm1,%%xmm0\n\t"
		"movdqa %%xmm5,%%xmm7\n\t"

		"paddq %%xmm2,%%xmm1\n\t"
		"psubq %%xmm2,%%xmm0\n\t"
		"paddq %%xmm6,%%xmm5\n\t"
		"psubq %%xmm6,%%xmm7\n\t"
		"psubq %%xmm4,%%xmm1\n\t"
		"psubq %%xmm4,%%xmm5\n\t"

		"movdqa %%xmm0,%%xmm3\n\t"	
		"movdqa %%xmm1,%%xmm2\n\t"
		"movdqa %%xmm7,%%xmm8\n\t"	
		"movdqa %%xmm5,%%xmm6\n\t"

		"psrad $31,%%xmm3\n\t"
		"psrad $31,%%xmm2\n\t"
		"psrad $31,%%xmm8\n\t"
		"psrad $31,%%xmm6\n\t"

		"pshufd $245,%%xmm3,%%xmm3\n\t"
		"pshufd $245,%%xmm2,%%xmm2\n\t"
		"pshufd $245,%%xmm8,%%xmm8\n\t"
		"pshufd $245,%%xmm6,%%xmm6\n\t"

		"pand %%xmm4,%%xmm3\n\t"
		"pand %%xmm4,%%xmm2\n\t"
		"pand %%xmm4,%%xmm8\n\t"
		"pand %%xmm4,%%xmm6\n\t"

		"paddq %%xmm3,%%xmm0\n\t"
		"paddq %%xmm2,%%xmm1\n\t"
		"paddq %%xmm8,%%xmm7\n\t"
		"paddq %%xmm6,%%xmm5\n\t"

		"movdqa %%xmm0,(%1)\n\t"
		"movdqa %%xmm1,(%0)\n\t"
		"movdqa %%xmm7,16(%1)\n\t"
		"movdqa %%xmm5,16(%0)\n\t"
//------------------------------------------------------------//
		"movdqa 32(%0),%%xmm1\n\t"
		"movdqa 32(%1),%%xmm2\n\t"
		"movdqa 48(%0),%%xmm5\n\t"
		"movdqa 48(%1),%%xmm6\n\t"
		"movdqa %%xmm1,%%xmm0\n\t"
		"movdqa %%xmm5,%%xmm7\n\t"

		"paddq %%xmm2,%%xmm1\n\t"
		"psubq %%xmm2,%%xmm0\n\t"
		"paddq %%xmm6,%%xmm5\n\t"
		"psubq %%xmm6,%%xmm7\n\t"
		"psubq %%xmm4,%%xmm1\n\t"
		"psubq %%xmm4,%%xmm5\n\t"

		"movdqa %%xmm0,%%xmm3\n\t"	
		"movdqa %%xmm1,%%xmm2\n\t"
		"movdqa %%xmm7,%%xmm8\n\t"	
		"movdqa %%xmm5,%%xmm6\n\t"

		"psrad $31,%%xmm3\n\t"
		"psrad $31,%%xmm2\n\t"
		"psrad $31,%%xmm8\n\t"
		"psrad $31,%%xmm6\n\t"

		"pshufd $245,%%xmm3,%%xmm3\n\t"
		"pshufd $245,%%xmm2,%%xmm2\n\t"
		"pshufd $245,%%xmm8,%%xmm8\n\t"
		"pshufd $245,%%xmm6,%%xmm6\n\t"

		"pand %%xmm4,%%xmm3\n\t"
		"pand %%xmm4,%%xmm2\n\t"
		"pand %%xmm4,%%xmm8\n\t"
		"pand %%xmm4,%%xmm6\n\t"

		"paddq %%xmm3,%%xmm0\n\t"
		"paddq %%xmm2,%%xmm1\n\t"
		"paddq %%xmm8,%%xmm7\n\t"
		"paddq %%xmm6,%%xmm5\n\t"

		"movdqa %%xmm0,32(%1)\n\t"
		"movdqa %%xmm1,32(%0)\n\t"
		"movdqa %%xmm7,48(%1)\n\t"
		"movdqa %%xmm5,48(%0)\n\t"
		:
		:"D"(a),"S"(b),"m"(*prime)
		:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","memory");
   return;
}
/**
 * AddMod:
 * @a: a int.
 * @b: b int.
 * @p: Prime number.
 * 
 * Return value: (a+b) mod p. 
 **/
static 
sfixn AddMod(sfixn a, sfixn b, sfixn p){
  sfixn r = a + b;
        r -= p;
        r += (r >> BASE_1) & p;
   return r;
}



/**
 * SubMod:
 * @a: a int.
 * @b: b int.
 * @p: Prime number.
 * 
 * Return value: (a-b) mod p. 
 **/
static 
sfixn SubMod(sfixn a, sfixn b, sfixn p){
   sfixn r = a - b;
   r += (r >> BASE_1) & p;
   return r;}

static inline
void  AddSubSSEModInplace(sfixn* a,sfixn* b, const sfixn* prime){
	#if ASM
	asm("movdqa %0,%%xmm1\n\t"
		"movdqa %1,%%xmm2\n\t"
		"movdqa %2,%%xmm4\n\t"
		"movdqa %%xmm1,%%xmm0\n\t"

		"paddq %%xmm2,%%xmm1\n\t"
		"psubq %%xmm4,%%xmm1\n\t"
		"psubq %%xmm2,%%xmm0\n\t"

		"movdqa %%xmm0,%%xmm3\n\t"	
		"movdqa %%xmm1,%%xmm2\n\t"
		"psrad $31,%%xmm3\n\t"
		"psrad $31,%%xmm2\n\t"

		"pshufd $245,%%xmm3,%%xmm3\n\t"
		"pshufd $245,%%xmm2,%%xmm2\n\t"

		"pand %%xmm4,%%xmm3\n\t"
		"pand %%xmm4,%%xmm2\n\t"

		"paddq %%xmm3,%%xmm0\n\t"
		"paddq %%xmm2,%%xmm1\n\t"

		"movdqa %%xmm0,%1\n\t"
		"movdqa %%xmm1,%0\n\t"
		:
		:"m"(*a),"m"(*b),"m"(*prime)
		:"xmm0","xmm1","xmm2","xmm3","xmm4","memory");
	#else
	sfixn u1 = AddMod(*a,*b,*prime);
	sfixn u2 = SubMod(*a,*b,*prime);
	sfixn u3 = AddMod(*(a+1),*(b+1),*prime);
	sfixn u4 = SubMod(*(a+1),*(b+1),*prime);
	*a = u1;
	*(a+1) = u3;
	*b = u2;
	*(b+1) = u4;
	#endif
   return;
}
static inline
void  AddSubSSEMod(sfixn* r1,sfixn* r2,sfixn* a,sfixn* b, const sfixn* prime){
	#if ASM
	asm("movdqa %2,%%xmm1\n\t"
		"movdqa %3,%%xmm2\n\t"
		"movdqa %4,%%xmm4\n\t"
		"movdqa %%xmm1,%%xmm0\n\t"

		"paddq %%xmm2,%%xmm1\n\t"
		"psubq %%xmm4,%%xmm1\n\t"
		"psubq %%xmm2,%%xmm0\n\t"

		"movdqa %%xmm0,%%xmm3\n\t"	
		"movdqa %%xmm1,%%xmm2\n\t"
		"psrad $31,%%xmm3\n\t"
		"psrad $31,%%xmm2\n\t"

		"pshufd $245,%%xmm3,%%xmm3\n\t"
		"pshufd $245,%%xmm2,%%xmm2\n\t"

		"pand %%xmm4,%%xmm3\n\t"
		"pand %%xmm4,%%xmm2\n\t"

		"paddq %%xmm3,%%xmm0\n\t"
		"paddq %%xmm2,%%xmm1\n\t"

		"movdqa %%xmm0,%1\n\t"
		"movdqa %%xmm1,%0\n\t"
		:"=m"(*r1),"=m"(*r2)
		:"m"(*a),"m"(*b),"m"(*prime)
		:"xmm0","xmm1","xmm2","xmm3","xmm4","memory");
	#else
	*r1 = AddMod(*a,*b,*prime);
	*r2 = SubMod(*a,*b,*prime);
	*(r1+1) = AddMod(*(a+1),*(b+1),*prime);
	*(r2+1) = SubMod(*(a+1),*(b+1),*prime);
	#endif
   return;
}


/**
 * NegMod:
 * @a: a int.
 * @p: Prime number.
 * 
 * Return value: (-a) mod p. 
 **/
static 
sfixn NegMod(sfixn a,  sfixn p){
   sfixn r = - a;
   r += (r >> BASE_1) & p;
   return r;}



// This may make things slower!
/**
 * MulMod:
 * @a: An int.
 * @b: An int.
 * @n: An int.
 * 
 * Return value: (a*b) mod n. 
 **/
static  sfixn MulMod(sfixn a, sfixn b, sfixn n)
{
  sfixn q, res;

  /*
   double ninv=1/(double)n;
   q  = (sfixn) ((((double) a) * ((double) b)) * ninv);
   res = a*b - q*n;
  */

  __float128 ninv=1/(__float128)n;
  q  = (sfixn) ((((__float128) a) * ((__float128) b)) * ninv);
  res = (sfixn) (((longfixnum) a)*((longfixnum) b) - ((longfixnum) q)*((longfixnum) n));

  res += (res >> BASE_1) & n;
  res -= n; 
  res += (res >> BASE_1) & n;
  return res;
}

static  sfixn MulModWithInv(sfixn a, sfixn b, sfixn n,__float128 ninv)
{
  sfixn q, res;

  /*
   double ninv=1/(double)n;
   q  = (sfixn) ((((double) a) * ((double) b)) * ninv);
   res = a*b - q*n;
  */

  q  = (sfixn) ((((__float128) a) * ((__float128) b)) * ninv);
  res = (sfixn) (((longfixnum) a)*((longfixnum) b) - ((longfixnum) q)*((longfixnum) n));

  res += (res >> BASE_1) & n;
  res -= n; 
  res += (res >> BASE_1) & n;
  return res;
}


/**
 * egcd:
 * @x: An int.
 * @y: An int.
 * @ao: A pointer of int.
 * @bo: A pointer of int.
 * Extended Euclidean Gcd of machine integers.
 * Return value: (*bo)x + (*ao)y = *vo. 
 **/
static 
void egcd (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo)  {
  sfixn tmp;
  sfixn A,B,C,D,u,v,q;

  u = y; v = x; 
  A=1; B=0; 
  C=0; D=1;

  do {
    q = u / v;
    tmp = u;
    u = v;
    v = tmp - q*v;
    tmp = A;
    A = B;
    B = tmp - q*B;
    tmp = C;
    C = D;
    D = tmp - q*D;
  } while (v != 0);
  *ao=A; 
  *bo=C; 
  *vo=u;
}


/**
 * inverseMod:
 * @n: A sfixn in Z/pZ.
 * @p: A prime.
 * 
 * computer the invers of n modulo p, where n is reduced with respect to p.
 * Return value: (1/n) mod p 
 **/
static 
sfixn inverseMod(sfixn n, sfixn p){
  sfixn a, b, v;
  egcd(n, p, &a, &b, &v);
  if (b < 0)
    b += p;
  return b % p;
}


/**
 * DivMod:
 * @a: A sfixn number in Z/pZ.
 * @b: A sfixn number in Z/pZ.
 * @p: A prime number.
 * Compute the faction a/b modulo p.
 * Return value: (a/b) mod p. 
 **/
static 
sfixn DivMod(sfixn a, sfixn b, sfixn p){
  return MulMod(a, inverseMod(b, p), p);
}


/**
 * QuoMod:
 * @a: A sfixn number in Z/pZ.
 * @b: A sfixn number in Z/pZ.
 * @p: A prime number.
 * Compute the faction a/b modulo p.
 * Return value: (a/b) mod p. 
 **/
static 
sfixn QuoMod(sfixn a, sfixn b, sfixn p){
  return MulMod(a, inverseMod(b, p), p);
}


/**
 * PowerMod:
 * @a: A sfixn number in Z/nZ.
 * @ee: A sfixn number.
 * @n: A moduli.
 * Compute the power a^ee modulo n.
 * Return value: (a^ee) mod n. 
 **/
static 
sfixn PowerMod(sfixn a, sfixn ee, sfixn n)
{
   sfixn x, y;

   usfixn e;

   if (ee == 0) return 1;

   if (ee < 0)
      e = - ((usfixn) ee);
   else
      e = ee;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = MulMod(x, y, n);
      y = MulMod(y, y, n);
      e = e >> 1;
   }

   if (ee < 0) x = inverseMod(x, n);

   return x;
}





static  uint32
partialBitRev(register uint32 x, int32 n)
{
  register uint32 y = 0x55555555;
  x = (((x >> 1) & y) | ((x & y) << 1));
  y = 0x33333333;
  x = (((x >> 2) & y) | ((x & y) << 2));
  y = 0x0f0f0f0f;
  x = (((x >> 4) & y) | ((x & y) << 4));
  y = 0x00ff00ff;
  x = (((x >> 8) & y) | ((x & y) << 8));
  y =((x >> 16) | (x << 16));
  return y>>(32-n);
}



/**
 * MontMulMod_OPT2_AS_GENE:
 * @a: A fixnum.
 * @b: A fixnum.
 * @pPtr: Information for the prime number p.
 * 
 * A improved Montgamoney Trick only works for Fourier Primes.
 * Please see the JSC paper for details.
 * Return value: (a*b)/R mod p, where R is the next power of 2 of p.
 **/
static inline sfixn
MontMulMod_OPT2_AS_GENE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr){
  sfixn q2=pPtr->c_sft, c;
  // b has been pre-leftShifted Base_Rpow bits.
  MulHiLoUnsigned(&a, &b);
  // now a keeps M1 quo R.  
  // q2=c<<{Npow}.
  MulHiLoUnsigned(&b, &q2);
  // now b keeps M3 quo R;
  a-=b;
  // load c;
  c=pPtr->c;
  // be careful the sign of q2.
  q2=(usfixn)q2>>pPtr->Base_Npow;
  q2*=c;
  // now q2 keeps M2 quo R;
  // load P.
  c=pPtr->P;
  a += (a >> BASE_1) & c;
  a+=q2;
  if(a<c)
	return a;
  else
	return a-c;
  // compute (M1+M2-M3) quo R.
  // the result is in range  -p<a<2p.
  // we need following 3 lines to reduce the error.
  //a += (a >> BASE_1) & c;
  //a -= c;
  //a += (a >> BASE_1) & c;
  //return a;
}

#ifndef MONTMULMOD3
#define MONTMULMOD3
static inline sfixn MontMulModSpe_OPT3_AS_GENE_globalfunc(sfixn a,sfixn b,sfixn inv,sfixn prime){
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
		: "=d" (a)
		: "a"(a),"rm"(b),"rm"(inv),"rm"(prime),"c"(BASE_1)
		:"rsi","rdi");
	return a;
}
#endif
// return the first one.
// r2nd fetch the second one.

static  sfixn
MontMulMod_OPT2_AS_Double_GENE(sfixn * r2nd, sfixn a, sfixn b, sfixn x, sfixn y, MONTP_OPT2_AS_GENE * pPtr){
  sfixn q2=pPtr->c_sft, q22=q2;

  // b has been pre-leftShifted Base_Rpow bits.
  MulHiLoUnsigned(&a, &b);
  MulHiLoUnsigned(&b, &q2);
  a-=b;
  MulHiLoUnsigned(&x, &y);
  MulHiLoUnsigned(&y, &q22);
  x-=y;

  y=pPtr->Base_Npow;
  b=pPtr->c;

  // be careful the sign of q2.
  q2=(usfixn)q2>>y;
  q2*=b;
  a+=q2;

  q22=(usfixn)q22>>y;
  y=pPtr->P;
  q22*=b;
  x+=q22;
  // now q2 keeps M2 quo R;
  // load P.


  // compute (M1+M2-M3) quo R.
  // the result is in range  -p<a<2p.
  // we need following 3 lines to reduce the error.


  a += (a >> BASE_1) & y;
  a -= y;
  a += (a >> BASE_1) & y;


  x += (x >> BASE_1) & y;
  x -= y;
  x += (x >> BASE_1) & y;

  *r2nd=x;

  return a;
}


/**
 * MontMulMod_OPT2_AS_GENE_SPE:
 * @a: A fixnum.
 * @b: A fixnum.
 * @pPtr: Information for the prime number p.
 * 
 * A improved Montgamoney Trick only works for Fourier Primes.
 * Please see the JSC paper for details.
 * This rountine is for special Fourier Primes who are in the shape of N+1,
 * where N is a power of 2.
 * Return value: (a*b)/R mod p, where R is the next power of 2 of p.
 **/
static  sfixn
MontMulMod_OPT2_AS_GENE_SPE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr){
  sfixn cpow=pPtr->c_pow, q3=pPtr->R_Npow, p=pPtr->Base_Rpow; 
  // b has been pre-leftShifted Base_Rpow bits.
  MulHiLoUnsigned(&a, &b);
  // now a keeps M1 quo R.

  b=((usfixn)b)>>p;
  b+=b<<cpow;
  p=pPtr->R_Nmask;
  q3=b>>q3;
  a-=q3;
  q3=pPtr->N2_Rpow;
  b&=p;
  b+=b<<cpow;
  b<<=q3;
  p=pPtr->P;
  a+=b;
  // the result is in range  -p<a<2p.
  // we need following 3 lines to reduce the error.
  a += (a >> BASE_1) & p;
  a -= p;
  a += (a >> BASE_1) & p;
  return a;
}


// return the first one.
// r2nd fetch the second one.
static  sfixn
MontMulMod_OPT2_AS_Double_GENE_SPE(sfixn * r2nd, sfixn a, sfixn b, sfixn x, sfixn y, MONTP_OPT2_AS_GENE * pPtr){

  sfixn cpow=pPtr->c_pow, q3=pPtr->R_Npow,  p=pPtr->Base_Rpow; 
  // b has been pre-leftShifted Base_Rpow bits.
  MulHiLoUnsigned(&a, &b);
  b=((usfixn)b)>>p;
  b+=b<<cpow;

  MulHiLoUnsigned(&x, &y);
  // now a keeps M1 quo R.
  y=((usfixn)y)>>p;
  y+=y<<cpow;

  p=pPtr->R_Nmask;


  a-=b>>q3;
  x-=y>>q3;

  q3=pPtr->N2_Rpow;


  b&=p;
  y&=p;

  b+=b<<cpow;
  b<<=q3;
  a+=b;

  y+=y<<cpow;
  y<<=q3;

  p=pPtr->P;

  x+=y;
  // the result is in range  -p<a<2p.
  // we need following 3 lines to reduce the error.
  a += (a >> BASE_1) & p;
  a -= p;
  a += (a >> BASE_1) & p;

  x += (x >> BASE_1) & p;
  x -= p;
  x += (x >> BASE_1) & p;
  
  *r2nd=x;

  return a;

}



//=========================================================
// A technical step before using Montgomery trick.
//=========================================================
static 
sfixn convertToMondMulR(sfixn r,  MONTP_OPT2_AS_GENE * pPtr){
      sfixn p=pPtr->P, R=(1L<<pPtr->Rpow)%p;
      R=(MulMod(r, R, p))<<pPtr->Base_Rpow;
      return R;
}



/**
 * MultiNumbPolyMul_1:
 * @r: r is a scalar.
 * @f: f is a C-Cube polynomial.
 * @pPtr: Information for the prime p.
 * 
 * Compute the product of r and f modulo p.
 * Return value: r*f mod p.
 **/
static 
void MultiNumbPolyMul_1(sfixn r, preFFTRep * f,  MONTP_OPT2_AS_GENE * pPtr){
   register sfixn i;
   sfixn R;
   R=convertToMondMulR(r, pPtr);
   for(i=0; i<SIZ(f); i++)  DATI(f, i)=MontMulMod_OPT2_AS_GENE( DATI(f, i), R, pPtr);
}


//=====================================================
//  copying data from one dense multivariate polynomial
//  to the one.
//=====================================================

static 
void fromtofftRep(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs){
  int32 i;
  sfixn d;
  int32 tmpRes=0, tmpCoeffs=0;
  //printf("N=%ld, Using fromtofftRep().\n", N);
   if(N==0){
    res[0]=coeffs[0];
    return;}
  d=dgs[N];
  if(N==1){
    for(i=0; i<=d; i++){ res[i]=coeffs[i];}
    return;}
  for(i=0; i<=d; i++){
    tmpCoeffs=i*ccum[N];
    tmpRes=i*rccum[N];
    fromtofftRep(N-1, rccum, res+tmpRes, ccum, dgs, coeffs+tmpCoeffs); }
}



/**
 * zeroCoefp:
 * @coefPtr: A vector.
 * @coefSiz: Size of 'coefPtr'.
 * 
 * Test if vector 'coefPtr' of size 'coefSiz' is a zero vector.
 * Return value: 1 for true, 0 for false. 
 **/
static 
int32 zeroCoefp(sfixn * coefPtr, sfixn coefSiz){
  register int32 i;
  for(i=0;i<coefSiz;i++){ if (coefPtr[i]) return 0;}
  return 1;
}


/**
 * shrinkDeg:
 * @deg: The original degree.
 * @poly: The coefficient vector.
 * @coefSiz: The size of a coefficient.
 * 
 * View the 'poly' as an univariate polynomial has 'deg'+1 coefficients.
 * Each coefficient has size of 'coefSiz'.
 * By ignoring the leading zeros, we compute the actual degree.
 * Return value: The actual degree.
 **/
static 
sfixn shrinkDeg(sfixn deg, sfixn * poly, sfixn coefSiz){
  sfixn * tmpPtr=poly+deg*coefSiz;
  while((zeroCoefp(tmpPtr, coefSiz))&& deg>0) {deg--; tmpPtr-=coefSiz;}
  return deg;
}


/**
 * shrinkDegUni:
 * @deg: The original degree.
 * @coef: The coefficient vector.
 * 
 * By ignoring the leading zeros, we compute the actual degree.
 * Return value: The actual degree.
 **/
static 
sfixn shrinkDegUni(sfixn deg, sfixn * cof){
  while((deg>0) && (!cof[deg])) deg--;
  return deg;
}


static 
void nextMCoefData(preFFTRep * Ptr, sfixn N, sfixn M){
   DAT(Ptr)+=M*(CUMI(Ptr, N));
}



/**
 * setLeadingCoefMultiOne:
 * @N: X_N is the main variable of 'f'.
 * @f: a C-Cube polynomial.
 * 
 * Set f's leading coefficient to be 1.
 * Return value: f.
 **/
static 
preFFTRep *
setLeadingCoefMultiOne(sfixn N, preFFTRep* f){
  sfixn d;
  register sfixn i;
  d=shrinkDeg(BUSZSI(f, N), DAT(f), CUMI(f, N));
  nextMCoefData(f,N,d);
  for(i=0; i<CUMI(f,N); i++) DATI(f, i)=0;
  DATI(f,0)=1;
  DAT(f)=DEFDAT(f);
  return f;
}


static 
void
MultiNumbPolyMulMonicize_1(sfixn N, sfixn r, preFFTRep * f,  MONTP_OPT2_AS_GENE * pPtr){
  MultiNumbPolyMul_1(r, f,  pPtr);
  setLeadingCoefMultiOne(N,f);
}




static 
void subEqDgPoly_inner_1
 (sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int32 selector)
{
  int32 i, offset=0;
  if(N==1){
    if(selector==1){
      for(i=0;i<=dgs[1];i++) data1[i]=SubMod(data1[i],data2[i],p);}
    else{
      for(i=0;i<=dgs[1];i++) data2[i]=SubMod(data1[i],data2[i],p);}
    return;}
  for(i=0; i<=dgs[N]; i++){
    offset=accum[N]*i;
    subEqDgPoly_inner_1(N-1, dgs, accum, data1+offset, data2+offset, p, selector);
  } 
}

static 
void subEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p, int32 selector){
  subEqDgPoly_inner_1(N, BUSZS(Ptr1), CUM(Ptr1), DAT(Ptr1), DAT(Ptr2), p, selector);
}


// return 1 means It IS a Zero.
/**
 * zeroPolyp:
 * @polyPtr: A C-Cube polynomial.
 * 
 * To check if 'polyPtr' is a zero polynomial
 * Return value: 1 for true. 0 for false.
 **/
static 
sfixn
zeroPolyp(preFFTRep * polyPtr){
  register int32 i;
  if((!(N(polyPtr))) && (!(DATI(polyPtr,0)))) return 1; 
  for(i=0;i<(SIZ(polyPtr));i++){ if (DATI(polyPtr, i)) return 0;}
  return 1;
}


// return 1 means It IS a Zero.
/**
 * constantPolyp:
 * @polyPtr: A C-Cube polynomial.
 * 
 * To check if 'polyPtr' is a constant polynomial
 * Return value: 1 for true. 0 for false.
 **/
static 
sfixn
constantPolyp(preFFTRep * polyPtr){
  register int32 i;
  if((!(N(polyPtr))) && ((DATI(polyPtr, 0))==1)) return 1; 
  for(i=1;i<(SIZ(polyPtr));i++){ if (DATI(polyPtr, i)) return 0;}
  return 1;
}



/**
 * PolyCleanData:
 * @prePtr:  A C-Cube polynomial.
 * 
 * make prePtr a zero polynomial.
 * Return value: 
 **/
static 
void PolyCleanData(preFFTRep * prePtr){
  register int32 i;
  for(i=0;i<SIZ(prePtr);i++){
    DATI(prePtr, i)=0;
  }
}


/**
 * negatePoly_1:
 * @prePtr:  A C-Cube polynomial.
 * @p: A prime number.
 *
 * Negate all coefficients of 'prePtr'.
 * Return value: 
 **/
static 
void negatePoly_1(preFFTRep * prePtr, sfixn p){
  register int32 i;
  for(i=0;i<SIZ(prePtr);i++){
    DATI(prePtr, i)=p-DATI(prePtr, i);
  }
}


/**
 * setPolyOne:
 * @Ptr: A C-Cube polynomial.
 * 
 * Set 'Ptr' to a constant polynomial has value 1.
 * Return value: 
 **/
static 
void
setPolyOne(preFFTRep *Ptr){
 PolyCleanData(Ptr);
 DATI(Ptr, 0)=1;
}



static 
void addEqDgPoly_inner
    (sfixn N, sfixn *dgs, sfixn *accum, sfixn *data1, sfixn *data2, sfixn p)
{
  int32 i, offset=0;
  if(N==1){
    for(i=0;i<=dgs[1];i++) data1[i]=AddMod(data1[i],data2[i],p);
    return;}
  for(i=0; i<=dgs[N]; i++){
    offset=accum[N]*i;
    addEqDgPoly_inner(N-1, dgs, accum, data1+offset, data2+offset, p);
  } 
}


/**
 * addEqDgPoly_1:
 * @N: Number of variables in 'Ptr1' and 'Ptr2'.
 * @Ptr1: A C-Cube polynomial.
 * @Ptr2: A C-Cube polynomial.
 * @p: A prime number.
 * 
 * Suppose 'Ptr1' and 'Ptr2' has the same dimension and size.
 * Compute the sum of them.
 * Return value: Ptr1 = Ptr1 + Ptr2;
 **/
static 
void addEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p){
  addEqDgPoly_inner(N, BUSZS(Ptr1), CUM(Ptr1), DAT(Ptr1), DAT(Ptr2), p);
}



// we will use the smaller dgs which are 2nd. but accums should used prepectively.
static 
void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p)
{
  int32 i, offset1=0, offset2=0;
  if(N==1){
    for(i=0;i<=dgs2[1];i++) data1[i]=SubMod(data1[i],data2[i],p);
    return;}
  for(i=0; i<=dgs2[N]; i++){
    offset1=accum1[N]*i;
    offset2=accum2[N]*i;
    subPoly_inner_1(N-1, accum1, dgs2, accum2, data1+offset1, data2+offset2, p);
  } 

}

// this in-place fun suppose Ptr1 is the larger buffer on all dimensions.

/**
 * subPoly_1:
 * @N: Number of variables in 'Ptr1' and 'Ptr2'.
 * @Ptr1: A C-Cube polynomial.
 * @Ptr2: A C-Cube polynomial.
 * @p: A prime number.
 * 
 * Suppose 'Ptr1' has the same dimension and larger(or equal) size on each dimension.
 * Compute the difference of them.
 * Return value: Ptr1 = Ptr1 - Ptr2;
 **/
static 
void subPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p){
  subPoly_inner_1(N, CUM(Ptr1), BUSZS(Ptr2), CUM(Ptr2), DAT(Ptr1), DAT(Ptr2), p);
}


static 
void setData(preFFTRep * Ptr, sfixn * dataPtr){ 
   DAT(Ptr)=dataPtr;
}


static 
void backupData(preFFTRep * Ptr){
   TMPDAT(Ptr)=DAT(Ptr);
}


static 
void restoreData(preFFTRep * Ptr){
   DAT(Ptr)=TMPDAT(Ptr);
}


static 
void resumeData(preFFTRep * Ptr){ 
   DAT(Ptr)=DEFDAT(Ptr);
}

static 
void decreaseOneDim(preFFTRep * Ptr){
  SIZ(Ptr)=CUMI(Ptr, N(Ptr));
  (N(Ptr))--;
}


static 
void increaseOneDim(preFFTRep * Ptr){
  (N(Ptr))++;
  SIZ(Ptr)=(SIZ(Ptr))*( BUSZSI(Ptr, N(Ptr)) +1);
}


static 
void resumeDim(preFFTRep * Ptr){
  N(Ptr)=DFN(Ptr);
  SIZ(Ptr)=DFSIZ(Ptr);
  DAT(Ptr)=DEFDAT(Ptr);
}


static 
void nextCoefData(preFFTRep * Ptr, sfixn N){
   DAT(Ptr)+=CUMI(Ptr, N);
}




//yu: is not used anywhere,  
//wrong output for p=38349804861915137 = 2^30 * 35716039 +1
/** 
 * MontMulMod:
 * @a: A sfixn.
 * @b: A sfixn.
 * @pPtr: Information for prime p.
 * 
 * Compute 'a'*'b' modulo p Based on improved Montgomery trick.
 * Return value a*b mod p: 
 **/
static 
sfixn MontMulMod(sfixn a, sfixn b, MONTP_OPT2_AS_GENE *pPtr){
  sfixn brsft;
  brsft=MontMulMod_OPT2_AS_GENE(b, pPtr->R2BRsft, pPtr)<<(pPtr->Base_Rpow);
  return MontMulMod_OPT2_AS_GENE(a,brsft,pPtr);
}


/**
 * coMulVec_1:
 * @co: A coefficient.
 * @deg: Degree of univeriate polynomial f.
 * @vec: Coefficient vector of the univeriate polynomial f.
 * @pPtr: The Info for a prime p.
 * 
 * compute the product of 'co' and 'f' modulo p.
 * Return value: co*f mod p. 
 **/
static 
void
coMulVec_1(sfixn co, sfixn deg, sfixn *vec, MONTP_OPT2_AS_GENE *pPtr){
    int32 i;
    sfixn tmp;
    tmp=(MontMulMod_OPT2_AS_GENE(co,pPtr->R2BRsft,pPtr))<<pPtr->Base_Rpow;
    deg = shrinkDegUni(deg, vec);
    for(i=0; i<=deg; i++){
      vec[i] = MontMulMod_OPT2_AS_GENE(vec[i],tmp,pPtr);
    }

}



/**
 * coMulAddVec:
 * @co: A coefficient.
 * @deg: Degree of univeriate polynomial f1 and f2.
 * @vec1: Coefficient vector of the univeriate polynomial f1.
 * @vec2: Coefficient vector of the univeriate polynomial f2.
 * @pPtr: The Info for a prime p.
 * 
 * compute the product of 'co' and 'f' modulo p.
 * Return value: co*f mod p. 
 **/
static 
void
coMulAddVec(sfixn co, sfixn deg, sfixn *vec1, sfixn *vec2, MONTP_OPT2_AS_GENE *pPtr){
    int32 i;
    sfixn tmp;
    tmp=(MontMulMod_OPT2_AS_GENE(co,pPtr->R2BRsft,pPtr))<<pPtr->Base_Rpow;
    deg = shrinkDegUni(deg, vec2);
    for(i=0; i<=deg; i++){
      vec1[i] = AddMod(vec1[i], MontMulMod_OPT2_AS_GENE(vec2[i], tmp, pPtr), pPtr->P);
    }

}



/**
 * getLeadingCoefMulti:
 * @N: X_N is the main variable of f.
 * @co: co is an empty buffer.
 * @f: A C-Cube polynomial.
 * Copy the leading coefficient of 'f' into 'co'.
 * Return value: co -- the leading coefficient of 'f'.
 **/
static 
preFFTRep *
getLeadingCoefMulti(sfixn N, preFFTRep* co, preFFTRep* f){
  sfixn d;
  d=shrinkDeg(BUSZSI(f, N), DAT(f), CUMI(f, N));
  backupData(f);
  decreaseOneDim(f); 
  nextMCoefData(f,N,d);
  fromtofftRep(N-1,  CUM(co), DAT(co), CUM(f), BUSZS(f), DAT(f));
  increaseOneDim(f);
  restoreData(f);
  return co;
}


/**
 * setCoefMulti:
 * @N: X_N is the main variable of f.
 * @f: A C-Cube polynomial.
 * @co: A coefficient whoes main variable is X_{N-1}.
 * @j: a index number. 
 * Set f's j-th coefficient to be 'co'.
 * Return value: co.
 **/
static  preFFTRep *
setCoefMulti(sfixn N, preFFTRep* f, preFFTRep* co, sfixn j){
  register sfixn i;
  nextMCoefData(f,N,j);
  
  for(i=0; i<SIZ(co); i++) DATI(f, i)=DATI(co, i);
  DAT(f)=DEFDAT(f);
  //printf("f=\n");
  //printPoly(f);
  //fflush(stdout);
  return co;
}






#endif

/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


