namespace PBPAS{
#if AMD
static inline
void unrolledSpe8MontMul(sfixn* input1, sfixn* input2, MONTP_OPT2_AS_GENE * pPtr){
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
void unrolledSpe8MontMul(sfixn* input1, sfixn* input2, MONTP_OPT2_AS_GENE * pPtr){

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
void  unrolledSpe8AddSubSSEMod(sfixn* a,sfixn* b){
	sfixn prime[2] = {MY_PRIME,MY_PRIME};
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
//static sfixn MontMulModSpe_OPT2_AS_GENE(sfixn a, sfixn b){
//	sfixn q2=C_SFT;
//	MulHiLoUnsigned(&a, &b);
//	MulHiLoUnsigned(&b, &q2);
//	a-=b;
//	//q2=(usfixn)q2>>pPtr->Base_Npow;
//	q2=(usfixn)q2>>(BASE-NPOW);
//	q2*=SEE;
//	a += (a >> BASE_1) & MY_PRIME;
//	a+=q2;
//	if(a<MY_PRIME)
//		return a;
//	else
//		return a-MY_PRIME;
//}
static  sfixn MontMulModSpeDEBUG_OPT3_AS_GENE(sfixn a, sfixn b){
	usfixn q2;
	std::cout<<a << " "<<b<<std::endl;
	MulHiLoUnsigned(&a, &b);
	std::cout<<a << " "<<b<<std::endl;
	q2 = b*INV_PRIME;
	std::cout<<q2 << " "<<INV_PRIME<<std::endl;
	MulAndAddHiLoUnsigned(&a, &b,q2,MY_PRIME);
	std::cout<<a << " "<<b<<std::endl;
	a -= MY_PRIME;
	std::cout<<a << " "<<b<<std::endl;
	a += (a >> BASE_1) & MY_PRIME;
	std::cout<<a << " "<<b<<std::endl;
	return a;
}
sfixn testMontMul(usfixn u,usfixn v){
	usfixn w = MulMod(u,v>>RSFT,MY_PRIME);
	w = MulMod(w,RINV,MY_PRIME);
	return w;
}
inline sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){
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
		: "=d" (a)
		: "a"(a),"rm"(b),"b"((sfixn) INV_PRIME),"c"((sfixn) MY_PRIME)
		:"rsi","rdi");
	return a;
}
inline sfixn MontMulModSpe_OPT3_AS_GENE_INLINE2(sfixn a,sfixn b){
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
		: "a"(a),"d"(b),"b"((sfixn) INV_PRIME),"c"((sfixn) MY_PRIME)
		:"rsi","rdi");
	return a;
}
static sfixn AddModSpe(sfixn a, sfixn b){
	sfixn r = a + b;
	r -= MY_PRIME;
	r += (r >> BASE_1) & MY_PRIME;
	return r;
}
static sfixn SubModSpe(sfixn a, sfixn b){
	sfixn r = a - b;
	r += (r >> BASE_1) & MY_PRIME;
	return r;
}

static inline void  AddSubSpeSSEModInplace(sfixn* a,sfixn* b){
	sfixn u1 = AddModSpe(*a,*b);
	sfixn u2 = SubModSpe(*a,*b);
	sfixn u3 = AddModSpe(*(a+1),*(b+1));
	sfixn u4 = SubModSpe(*(a+1),*(b+1));
	*a = u1;
	*(a+1) = u3;
	*b = u2;
	*(b+1) = u4;
	return;
}
static inline
void  AddSubSpeSSEMod(sfixn* r1,sfixn* r2,sfixn* a,sfixn* b){
	*r1 = AddModSpe(*a,*b);
	*r2 = SubModSpe(*a,*b);
	*(r1+1) = AddModSpe(*(a+1),*(b+1));
	*(r2+1) = SubModSpe(*(a+1),*(b+1));
   return;
}
void inline FFT_2POINT(sfixn *A,sfixn *W){
	sfixn u = A[0];
	sfixn t = A[1];
	A[0] = AddModSpe(u,t);
	A[1] = SubModSpe(u,t);
}
void inline FFT_4POINT(sfixn *A,sfixn *W){
	sfixn *Wp = W + (4<<1)-4;
	sfixn w = A[1];
	A[1] = A[2];
	A[2] = w;
	sfixn u = A[0];
	sfixn t = A[1];
	A[0] = AddModSpe(u,t);
	A[1] = SubModSpe(u,t);
	u = A[2];
	t = A[3];
	A[2] = AddModSpe(u,t);
	A[3] = SubModSpe(u,t);
	A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3],*(Wp-3));

	AddSubSpeSSEModInplace(A,A+2);
}
void inline FFT_8POINT(sfixn *A,sfixn *W){
	sfixn *Wp = W + (8<<1)-4;
	sfixn u = A[0];
	sfixn t = A[1];
	A[0] = AddModSpe(u,t);
	A[1] = SubModSpe(u,t);
	u = A[2];
	t = A[3];
	A[2] = AddModSpe(u,t);
	A[3] = SubModSpe(u,t);
	u = A[4];
	t = A[5];
	A[4] = AddModSpe(u,t);
	A[5] = SubModSpe(u,t);
	u = A[6];
	t = A[7];
	A[6] = AddModSpe(u,t);
	A[7] = SubModSpe(u,t);
	A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3],*(Wp-3));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-3));

	AddSubSpeSSEModInplace(A,A+2);
	AddSubSpeSSEModInplace(A+4,A+6);
	A[5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5],*(Wp-11));
	A[6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6],*(Wp-10));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-9));
	AddSubSpeSSEModInplace(A,A+4);
	AddSubSpeSSEModInplace(A+2,A+6);
}
void inline FFT_16POINT(sfixn *A,sfixn *W){
	sfixn *Wp = W + (16<<1)-4;
	sfixn u = A[0];
	sfixn t = A[1];
	A[0] = AddModSpe(u,t);
	A[1] = SubModSpe(u,t);
	u = A[2];
	t = A[3];
	A[2] = AddModSpe(u,t);
	A[3] = SubModSpe(u,t);
	u = A[4];
	t = A[5];
	A[4] = AddModSpe(u,t);
	A[5] = SubModSpe(u,t);
	u = A[6];
	t = A[7];
	A[6] = AddModSpe(u,t);
	A[7] = SubModSpe(u,t);
	u = A[8];
	t = A[9];
	A[8] = AddModSpe(u,t);
	A[9] = SubModSpe(u,t);
	u = A[10];
	t = A[11];
	A[10] = AddModSpe(u,t);
	A[11] = SubModSpe(u,t);
	u = A[12];
	t = A[13];
	A[12] = AddModSpe(u,t);
	A[13] = SubModSpe(u,t);
	u = A[14];
	t = A[15];
	A[14] = AddModSpe(u,t);
	A[15] = SubModSpe(u,t);
	A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[3],*(Wp-3));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[7],*(Wp-3));
	A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[11],*(Wp-3));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[15],*(Wp-3));

	AddSubSpeSSEModInplace(A,A+2);
	AddSubSpeSSEModInplace(A+4,A+6);
	AddSubSpeSSEModInplace(A+8,A+10);
	AddSubSpeSSEModInplace(A+12,A+14);
	A[5] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[5],*(Wp-11));
	A[6] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[6],*(Wp-10));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[7],*(Wp-9));
	A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[13],*(Wp-11));
	A[14] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[14],*(Wp-10));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[15],*(Wp-9));

	AddSubSpeSSEModInplace(A,A+4);
	AddSubSpeSSEModInplace(A+2,A+6);
	AddSubSpeSSEModInplace(A+8,A+12);
	AddSubSpeSSEModInplace(A+10,A+14);

	A[9] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[9],*(Wp-27));
	A[10] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[10],*(Wp-26));
	A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[11],*(Wp-25));
	A[12] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[12],*(Wp-24));
	A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[13],*(Wp-23));
	A[14] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[14],*(Wp-22));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[15],*(Wp-21));
	AddSubSpeSSEModInplace(A,A+8);
	AddSubSpeSSEModInplace(A+2,A+10);
	AddSubSpeSSEModInplace(A+4,A+12);
	AddSubSpeSSEModInplace(A+6,A+14);
}
//void DFT_iter(int r,sfixn *A,sfixn *W){
//	sfixn *Wp = W+(FFT_THRESHOLD<<1)-4;
//	sfixn* Wt;
//	FFT_16_POINT
//	Wp = Wp - 28;
//	BLOCK_3:16,5
//	Wp = Wt-64 - 128;
//	BLOCK_3:128,8
//	Wp = Wt - 512-1024;
//}
DFT_ITERATIVE
void Shuffle(int n, sfixn *A, sfixn *B){
	int n2 = n>>1;
	for (int i=0; i<n2; i+=8){
		int i2 = i<<1;
		A[i] = A[i2];
		B[i] = A[i2+1];	
		A[i+1] = A[i2+2];
		B[i+1] = A[i2+3];
		A[i+2] = A[i2+4];
		B[i+2] = A[i2+5];
		A[i+3] = A[i2+6];
		B[i+3] = A[i2+7];
		A[i+4] = A[i2+8];
		B[i+4] = A[i2+9];
		A[i+5] = A[i2+10];
		B[i+5] = A[i2+11];
		A[i+6] = A[i2+12];
		B[i+6] = A[i2+13];
		A[i+7] = A[i2+14];
		B[i+7] = A[i2+15];
	}
	sfixn *A2 = A + n2;      
	memcpy(A2, B, n2*(sizeof(sfixn)));   
}

//void ArrayBitReversal(int n, sfixn *A, int *RevBidMap){
//	for (int i=0; i<n; i++){
//		int j = RevBidMap[i];
//		if (i<j){
//			sfixn t = A[i];
//			A[i] = A[j];
//			A[j] = t;
//		}
//	}
//}

sfixn testDFT(int n,int num,sfixn* A,sfixn *W){
	sfixn res = A[0];
	sfixn w = *(W+num);
	sfixn root = w;
	for(sfixn k=1;k<n;k++){
		
		sfixn t = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k],*(W+((num*k)%n)));
		res = AddModSpe(res,t);
	}
	return res;
}

void DFT_rec(int n, int r,sfixn *A,sfixn *W,sfixn *B){
	if (n==1) {    
		return;
	}
	else if(n==FFT_THRESHOLD){
		PBPAS::ArrayBitReversalSpe(A);
		DFT_iter(A,W);
		return;
	}
	else{
		int n2 = n>>1;
		int r1 = r-1;
		//Shuffle(n, A, B);
		sfixn *W2 = W+n; 
		DFT_rec(n2, r1, A, W2, B);
		DFT_rec(n2, r1, A+n2, W2, B);

		for (int k=n2; k<(n2<<1); k+=8){
			A[k+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+0],W[k-n2+0]);
			A[k+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+1],W[k-n2+1]);
			A[k+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+2],W[k-n2+2]);
			A[k+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+3],W[k-n2+3]);
			A[k+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+4],W[k-n2+4]);
			A[k+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+5],W[k-n2+5]);
			A[k+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+6],W[k-n2+6]);
			A[k+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+7],W[k-n2+7]);
			sfixn t = A[k-n2];
			sfixn u = A[k];
			A[k-n2] = AddModSpe(t,u);
			A[k] = SubModSpe(t,u);
			t = A[k+1-n2];
			u = A[k+1];
			A[k+1-n2] = AddModSpe(t,u);
			A[k+1] = SubModSpe(t,u);
			t = A[k+2-n2];
			u = A[k+2];
			A[k+2-n2] = AddModSpe(t,u);
			A[k+2] = SubModSpe(t,u);
			t = A[k+3-n2];
			u = A[k+3];
			A[k+3-n2] = AddModSpe(t,u);
			A[k+3] = SubModSpe(t,u);
			t = A[k+4-n2];
			u = A[k+4];
			A[k+4-n2] = AddModSpe(t,u);
			A[k+4] = SubModSpe(t,u);
			t = A[k+5-n2];
			u = A[k+5];
			A[k+5-n2] = AddModSpe(t,u);
			A[k+5] = SubModSpe(t,u);
			t = A[k+6-n2];
			u = A[k+6];
			A[k+6-n2] = AddModSpe(t,u);
			A[k+6] = SubModSpe(t,u);
			t = A[k+7-n2];
			u = A[k+7];
			A[k+7-n2] = AddModSpe(t,u);
			A[k+7] = SubModSpe(t,u);
		}
	}
}
void Shuffle2(int n, sfixn *A, sfixn *B){
	int n2 = n>>2;
	int n1 = n>>1;
	for (int i=0; i<n2; i+=8){
		int i2 = i<<2;
		A[i] = A[i2];
		B[i+n2] = A[i2+1];
		B[i] = A[i2+2];
		B[i+n1] = A[i2+3];
		A[i+1] = A[i2+4];
		B[i+n2+1] = A[i2+5];
		B[i+1] = A[i2+6];
		B[i+n1+1] = A[i2+7];
		A[i+2] = A[i2+8];
		B[i+n2+2] = A[i2+9];
		B[i+2] = A[i2+10];
		B[i+n1+2] = A[i2+11];
		A[i+3] = A[i2+12];
		B[i+n2+3] = A[i2+13];
		B[i+3] = A[i2+14];
		B[i+n1+3] = A[i2+15];
		A[i+4] = A[i2+16];
		B[i+n2+4] = A[i2+17];
		B[i+4] = A[i2+18];
		B[i+n1+4] = A[i2+19];
		A[i+5] = A[i2+20];
		B[i+n2+5] = A[i2+21];
		B[i+5] = A[i2+22];
		B[i+n1+5] = A[i2+23];
		A[i+6] = A[i2+24];
		B[i+n2+6] = A[i2+25];
		B[i+6] = A[i2+26];
		B[i+n1+6] = A[i2+27];
		A[i+7] = A[i2+28];
		B[i+n2+7] = A[i2+29];
		B[i+7] = A[i2+30];
		B[i+n1+7] = A[i2+31];
	}
	sfixn *A2 = A + n2;
	memcpy(A2, B, 3*n2*(sizeof(sfixn)));   
}
void Shuffle3(int n, sfixn *A, sfixn *B){
	int n3 = n>>3;
	for (int i=0; i<n3; i+=2){
		int i2 = i<<3;
		A[i] = A[i2];
		B[i+3*n3] = A[i2+1];
		B[i+n3] = A[i2+2];
		B[i+5*n3] = A[i2+3];
		B[i] = A[i2+4];
		B[i+(n3<<2)] = A[i2+5];
		B[i+(n3<<1)] = A[i2+6];
		B[i+6*n3] = A[i2+7];

		A[i+1] = A[i2+8];
		B[i+1+3*n3] = A[i2+9];
		B[i+1+n3] = A[i2+10];
		B[i+1+5*n3] = A[i2+11];
		B[i+1] = A[i2+12];
		B[i+1+(n3<<2)] = A[i2+13];
		B[i+1+(n3<<1)] = A[i2+14];
		B[i+1+6*n3] = A[i2+15];
	}
	sfixn *A2 = A + n3;
	memcpy(A2, B, 7*n3*(sizeof(sfixn)));   
}
void DFT_eff(int n, int r,sfixn *A,sfixn *W,sfixn *B){
	#if DEBUG
	sfixn test = testDFT(n,1,A,W);
	#endif
	if (n==1) {
		return;
	}
	else if(n==FFT_THRESHOLD){
		PBPAS::ArrayBitReversalSpe(A);
		DFT_iter(A,W);
		return;
	}
	else if(n>FFT_THRESHOLD){
		int n2 = n>>1;
		int r1 = r-1;
		sfixn twoFFTThreshold = FFT_THRESHOLD<<1;
		if(n>(twoFFTThreshold<<1)){
			sfixn pos = 0;
			sfixn m=n;
			unsigned long mask = MY_MASK1 ;//1317624576693539401;//MY_MASK;
			unsigned long mask2 = MY_MASK2;
			unsigned long mask3 = MY_MASK3;

			Shuffle(n,A,B);
			m>>=1;
			#if COMPOSESHUFFLE
			if((m&mask2)!=0){
				Shuffle(m,A,B);
				Shuffle(m,A+m,B);
				m>>=1;
			}
			else if((m&mask3) != 0){
				Shuffle2(m,A,B);
				Shuffle2(m,A+m,B);
				m>>=2;
			}
			while (pos<n){
				while(m>FFT_THRESHOLD){
					Shuffle3(m,A+pos,B);
					m>>=3;
				}
				pos += twoFFTThreshold<<2;
				m = ((pos^(pos-1))+1)>>1;
				if(m&mask3)
					m>>=2;
				if(m&mask2)
					m>>=1;
			}
			//if((n&mask) == 0){
			//	Shuffle(m,A,B);
			//	Shuffle(m,A+m,B);
			//	m>>=1;
			//}
			#else
			while (pos<n){
				//while(m>64*FFT_THRESHOLD){
				//	Shuffle(m,A+pos,B);
				//	Shuffle(m>>1,A+pos,B);
				//	Shuffle(m>>1,A+pos+(m>>1),B);
				//	m>>=2;
				//}
				while(m>FFT_THRESHOLD){
					//Shuffle2(m,A+pos,B);
					//m>>=2;
					Shuffle(m,A+pos,B);
					m>>=1;
				}
				pos += twoFFTThreshold;//<<1;
				m = ((pos^(pos-1))+1)>>1;
				//m -= ((m&mask)>>1);
			}
			#endif
		}
		else if(n==(twoFFTThreshold<<1)){
			Shuffle(n,A,B);
			sfixn m = n>>1;
			Shuffle(m,A,B);
			Shuffle(m,A+m,B);
		}
		else{
			Shuffle(n,A,B);
		}
		//Shuffle(n, A, B);
		sfixn *W2 = W+n; 
		DFT_rec(n2, r1, A, W2, B);
		DFT_rec(n2, r1, A+n2, W2, B);

		for (int k=n2; k<(n2<<1); k+=8){
			A[k+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+0],W[k-n2+0]);
			A[k+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+1],W[k-n2+1]);
			A[k+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+2],W[k-n2+2]);
			A[k+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+3],W[k-n2+3]);
			A[k+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+4],W[k-n2+4]);
			A[k+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+5],W[k-n2+5]);
			A[k+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+6],W[k-n2+6]);
			A[k+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+7],W[k-n2+7]);
			sfixn t = A[k-n2];
			sfixn u = A[k];
			A[k-n2] = AddModSpe(t,u);
			A[k] = SubModSpe(t,u);
			t = A[k+1-n2];
			u = A[k+1];
			A[k+1-n2] = AddModSpe(t,u);
			A[k+1] = SubModSpe(t,u);
			t = A[k+2-n2];
			u = A[k+2];
			A[k+2-n2] = AddModSpe(t,u);
			A[k+2] = SubModSpe(t,u);
			t = A[k+3-n2];
			u = A[k+3];
			A[k+3-n2] = AddModSpe(t,u);
			A[k+3] = SubModSpe(t,u);
			t = A[k+4-n2];
			u = A[k+4];
			A[k+4-n2] = AddModSpe(t,u);
			A[k+4] = SubModSpe(t,u);
			t = A[k+5-n2];
			u = A[k+5];
			A[k+5-n2] = AddModSpe(t,u);
			A[k+5] = SubModSpe(t,u);
			t = A[k+6-n2];
			u = A[k+6];
			A[k+6-n2] = AddModSpe(t,u);
			A[k+6] = SubModSpe(t,u);
			t = A[k+7-n2];
			u = A[k+7];
			A[k+7-n2] = AddModSpe(t,u);
			A[k+7] = SubModSpe(t,u);
		}
	}
	else{
		switch (n){
			case 2:
				FFT_2POINT(A,W);
				break;
			case 4:
				FFT_4POINT(A,W);
				break;
			case 8:
				PBPAS::ArrayBitReversalSpe8(A);
				FFT_8POINT(A,W);
				break;
			case 16:
				PBPAS::ArrayBitReversalSpe16(A);
				FFT_16POINT(A,W);
				break;
			FFT_CASES
			default:
				return;
		}
	}
	#if DEBUG
	if(test!=A[1]){
		std::cout<<n << " "<<test<< " "<<A[1]<<std::endl;
		exit(1);
	}
	#endif
}

void InvDFTKeepMont_eff(int n, int r, sfixn *A,sfixn *W,  sfixn *B,sfixn invn){
	DFT_eff(n,r,A,W,B);
	for(int i=0; i<n; i++) {
		A[i]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[i], invn);
	}
}
void InvDFT_eff(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn){
	DFT_eff(n,r,A,W,B);
	for(int i=0; i<n; i++) {
		A[i]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[i], invn);
	}
}
}
