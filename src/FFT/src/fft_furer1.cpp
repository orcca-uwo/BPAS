#include "../../../include/FFT/src/fft_furer1.h"
#include "../../../include/FFT/src/arraybitreversal.h"
#include "../../../include/FFT/src/modpn.h"
#include <iostream>
#include <string.h>
#define FFT_THRESHOLD 1024
#define FFT_THRESHOLD_LOG 10
#define COMPLETE 18266600088614731775u
#define MY_MASK0 2635249153387078656u
#define MY_MASK1 5270498306774157312u
#define MY_MASK2 10540996613548314624u
#define MY_ROOT0 6917529027641081792u
#define MY_ROOT1 3195247147627680576u
#define MY_ROOT2 7492325059881135552u
#define MY_ROOT3 9770180226917316736u
#define MY_ROOT4 4611686018427388032u
#define MY_ROOT5 8333967898440789248u
#define MY_ROOT6 4036889986187334272u
#define MY_ROOT7 1759034819151153088u
#define LOG_OF_POW 3
#define MY_POW 8
#define CONST_R2 151320947479648669u
#define CONST_R 108086391056891903u
#define BASE_RPOW 6
#define GENERATOR 3
namespace FURERPBPAS1{
static  sfixn MontMulModSpeDEBUG_OPT3_AS_GENE(sfixn a, sfixn b){
	usfixn q2;
	std::cout<<a << " "<<b<<std::endl;
	MulHiLoUnsigned(&a, &b);
	std::cout<<a << " "<<b<<std::endl;
	q2 = b*FURER_INV_PRIME1;
	std::cout<<q2 << " "<<FURER_INV_PRIME1<<std::endl;
	MulAndAddHiLoUnsigned(&a, &b,q2,FURER_MY_PRIME1);
	std::cout<<a << " "<<b<<std::endl;
	a -= FURER_MY_PRIME1;
	std::cout<<a << " "<<b<<std::endl;
	a += (a >> BASE_1) & FURER_MY_PRIME1;
	std::cout<<a << " "<<b<<std::endl;
	return a;
}
sfixn testMontMul(usfixn u,usfixn v){
	usfixn w = MulMod(u,v>>FURER_RSFT1,FURER_MY_PRIME1);
	w = MulMod(w,FURER_RINV1,FURER_MY_PRIME1);
	return w;
}
//  64-bits code
//
//
//inline sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){
//	asm("mulq %2\n\t"
//		"movq %%rax,%%rsi\n\t"
//		"movq %%rdx,%%rdi\n\t"
//		"addq %4,%%rdi\n\t"
//		"imulq %3,%%rax\n\t"
//		"mulq %5\n\t"
//		"xor %%rcx,%%rcx\n\t"
//		"add %%rsi,%%rax\n\t"
//		"adc %%rdi,%%rdx\n\t"
//		"setc %%cl\n\t"
//		"mov $18446744073709551615,%%rsi\n\t"
//		"addq %%rcx,%%rsi\n\t"
//		"andq %4,%%rsi\n\t"
//		"subq %%rsi,%%rdx\n\t"
//		: "=d" (a)
//		: "a"(a),"d"(b),"r"((unsigned long) FURER_INV_PRIME1),"r"((unsigned long) COMPLETE),"r"((unsigned long) FURER_MY_PRIME1)
//		:"rsi","rdi","rcx");
//	return a;
//}
//static sfixn AddModSpe(sfixn a, sfixn b){
//	a+=COMPLETE;
//	unsigned char c = 0;
//	asm("add %2,%%rax\n\t"
//		"setc %%dl"
//			:"=a"(a), "=d"(c)
//			:"r"(a),"a"(b)
//			:);
//	unsigned long g=((unsigned long)c)-1;
//	g&=COMPLETE;
//	return a-g;
//}
//static sfixn SubModSpe(sfixn a, sfixn b){
//	unsigned char c = 0;
//	asm("sub %2,%%rax\n\t"
//		"setc %%dl"
//			:"=a"(a), "=d"(c)
//			:"r"(b),"a"(a)
//			:);
//	unsigned long g=-((unsigned long)c);
//	g&=COMPLETE;
//	return a-g;
//}
#if SPECIALMONT
sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){
	asm("mulq %2\n\t"
	"movq %%rax,%%rsi\n\t"
	"movq %%rdx,%%rdi\n\t"
	"imulq %3,%%rax\n\t"
	"add %%rax,%%rsi\n\t"
	"adc $0,%%rdi\n\t"
	"shr $6,%%rax\n\t"
	"imulq $5,%%rax\n\t"
	"rolq $61,%%rax\n\t"
	"movq %%rax,%%rdx\n\t"
	"andq %5,%%rax\n\t"
	"andq %6,%%rdx\n\t"
	"add %%rsi,%%rax\n\t"
	"adc %%rdi,%%rdx\n\t"
	"subq %4,%%rdx\n\t"
	"mov %%rdx,%%rax\n\t"
	"sar $63,%%rax\n\t"
	"andq %4,%%rax\n\t"
	"addq %%rax,%%rdx\n\t"
	: "=d" (a)
	: "a"(a),"d"(b),"b"((sfixn) FURER_INV_PRIME1),"c"((sfixn) FURER_MY_PRIME1),"rm"(16140901064495857664u),"rm"(2305843009213693951)
	:"rsi","rdi");
	return a;
}
#else
sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){
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
	: "a"(a),"rm"(b),"b"((sfixn) FURER_INV_PRIME1),"c"((sfixn) FURER_MY_PRIME1)
	:"rsi","rdi");
	return a;
}
#endif
//inline sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){
//	asm("mulq %2\n\t"
//		"movq %%rax,%%rsi\n\t"
//		"movq %%rdx,%%rdi\n\t"
//		"imulq %3,%%rax\n\t"
//		"mulq %4\n\t"
//		"add %%rsi,%%rax\n\t"
//		"adc %%rdi,%%rdx\n\t"
//		"subq %4,%%rdx\n\t"
//		"mov %%rdx,%%rax\n\t"
//		"sar $63,%%rax\n\t"
//		"andq %4,%%rax\n\t"
//		"addq %%rax,%%rdx\n\t"
//		: "=d" (a)
//		: "a"(a),"rm"(b),"b"((sfixn) FURER_INV_PRIME1),"c"((sfixn) FURER_MY_PRIME1)
//		:"rsi","rdi");
//	return a;
//}
static sfixn AddModSpe(sfixn a, sfixn b){
	sfixn r = a + b;
	r -= FURER_MY_PRIME1;
	r += (r >> BASE_1) & FURER_MY_PRIME1;
	return r;
}
static sfixn SubModSpe(sfixn a, sfixn b){
	sfixn r = a - b;
	r += (r >> BASE_1) & FURER_MY_PRIME1;
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
static inline void AddSubSpeSSEMod(sfixn* r1,sfixn* r2,sfixn* a,sfixn* b){
	*r1 = AddModSpe(*a,*b);
	*r2 = SubModSpe(*a,*b);
	*(r1+1) = AddModSpe(*(a+1),*(b+1));
	*(r2+1) = SubModSpe(*(a+1),*(b+1));
   return;
}
unsigned long fastExp(unsigned long gen,unsigned long pow){
	if(pow==0)
		return 1;
	if(pow==1)
		return gen;
	if(pow & 1){
		return MontMulModSpe_OPT3_AS_GENE_INLINE(fastExp(MontMulModSpe_OPT3_AS_GENE_INLINE(gen,gen<<BASE_RPOW),pow>>1),gen<<BASE_RPOW);
	}
	else{
		return fastExp(MontMulModSpe_OPT3_AS_GENE_INLINE(gen,gen<<BASE_RPOW),pow>>1);
	}
}
void GenRootOfUnity(sfixn e, sfixn n, sfixn * rootsPtr){
	register int j;
	sfixn root, rootR, R_2=CONST_R2, BRsft=BASE_RPOW;
	rootR=fastExp(MontMulModSpe_OPT3_AS_GENE_INLINE(GENERATOR,R_2<<BRsft),(FURER_MY_PRIME1-1)>>e);
	rootsPtr[0]=CONST_R<<BRsft;
	rootsPtr[1]=rootR;
	for(j=2; j<n; j++) {
		rootsPtr[j] = MontMulModSpe_OPT3_AS_GENE_INLINE(rootsPtr[j-1], rootR<<BRsft);
		rootsPtr[j-1] <<= BRsft;
	}
	rootsPtr[n-1]<<=BRsft;
}
void RootsTableFurer(int n, int r,sfixn *T){
	GenRootOfUnity(r, n, T);
	sfixn *T_ptr = T;
	int ni = n;
	T_ptr = T_ptr + ni;
	ni = n>>1;
	for (int j=0; j<ni; j++){
		T_ptr[j] = T[(j<<1)];
	}
	int f;
	f = ((r-1)%LOG_OF_POW);
	if (f==0)
		f+=LOG_OF_POW;
	if(n>(FFT_THRESHOLD<<1)){
		f = (r-1-FFT_THRESHOLD_LOG)%LOG_OF_POW;
		if (f==0)
			f+=LOG_OF_POW;
		for (int i=f+1; i<=r-FFT_THRESHOLD_LOG; i+=LOG_OF_POW){
			T_ptr = T_ptr + ni;
			ni = n>>i;
			for (int j=0; j<ni; j++){
				T_ptr[j] = T[(j<<i)];
			}
		}
		f = (FFT_THRESHOLD_LOG%LOG_OF_POW);
		if (f==0)
			f+=LOG_OF_POW;
		f+=r-FFT_THRESHOLD_LOG;
		f--;
	}
	for (int i=f+1; i<r; i+=LOG_OF_POW){
		T_ptr = T_ptr + ni;
		ni = n>>i;
		for (int j=0; j<ni; j++){
			T_ptr[j] = T[(j<<i)];
		}
	}
}
void inline FFT_2POINT(sfixn *A,sfixn *W){
#if MODDEBUG
	printf("2points\n");
#endif
	sfixn u = A[0];
	sfixn t = A[1];
	A[0] = AddModSpe(u,t);
	A[1] = SubModSpe(u,t);
}
void inline FFT_4POINT(sfixn *A,sfixn *W){
#if MODDEBUG
	printf("4points\n");
#endif
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
	A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3],*(Wp-3));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-3));
	A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[11],*(Wp-3));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[15],*(Wp-3));

	AddSubSpeSSEModInplace(A,A+2);
	AddSubSpeSSEModInplace(A+4,A+6);
	AddSubSpeSSEModInplace(A+8,A+10);
	AddSubSpeSSEModInplace(A+12,A+14);
	A[5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5],*(Wp-11));
	A[6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6],*(Wp-10));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-9));
	A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[13],*(Wp-11));
	A[14] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[14],*(Wp-10));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[15],*(Wp-9));

	AddSubSpeSSEModInplace(A,A+4);
	AddSubSpeSSEModInplace(A+2,A+6);
	AddSubSpeSSEModInplace(A+8,A+12);
	AddSubSpeSSEModInplace(A+10,A+14);

	A[9] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[9],*(Wp-27));
	A[10] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[10],*(Wp-26));
	A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[11],*(Wp-25));
	A[12] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[12],*(Wp-24));
	A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[13],*(Wp-23));
	A[14] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[14],*(Wp-22));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[15],*(Wp-21));
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
sfixn testDFT(int n,int index,sfixn* A,sfixn *W){
//0<= num <n
	sfixn res = A[0];
	for(int k=1;k<n;k++){
		sfixn t = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k],*(W+((index*k)%n)));
		res = AddModSpe(res,t);
	}
	return res;
}

inline void small_fft0(sfixn* u0, sfixn* u1, sfixn* u2, sfixn* u3, sfixn* u4, sfixn* u5, sfixn* u6, sfixn* u7){
	sfixn t,u;
	u = AddModSpe(*u0,*u1);
	t = SubModSpe(*u0,*u1);
	*u0 = u;
	*u1 = t;
	u = AddModSpe(*u2,*u3);
	t = SubModSpe(*u2,*u3);
	*u2 = u;
	*u3 = t;
	u = AddModSpe(*u4,*u5);
	t = SubModSpe(*u4,*u5);
	*u4 = u;
	*u5 = t;
	u = AddModSpe(*u6,*u7);
	t = SubModSpe(*u6,*u7);
	*u6 = u;
	*u7 = t;
	*u2 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u2,MY_ROOT0);
	*u6 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u6,MY_ROOT0);
	*u3 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u3,MY_ROOT2);
	*u7 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u7,MY_ROOT2);
	u = AddModSpe(*u0,*u2);
	t = SubModSpe(*u0,*u2);
	*u0 = u;
	*u2 = t;
	u = AddModSpe(*u4,*u6);
	t = SubModSpe(*u4,*u6);
	*u4 = u;
	*u6 = t;
	u = AddModSpe(*u1,*u3);
	t = SubModSpe(*u1,*u3);
	*u1 = u;
	*u3 = t;
	u = AddModSpe(*u5,*u7);
	t = SubModSpe(*u5,*u7);
	*u5 = u;
	*u7 = t;
	*u4 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u4,MY_ROOT0);
	*u5 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u5,MY_ROOT1);
	*u6 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u6,MY_ROOT2);
	*u7 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u7,MY_ROOT3);
	u = AddModSpe(*u0,*u4);
	t = SubModSpe(*u0,*u4);
	*u0 = u;
	*u4 = t;
	u = AddModSpe(*u1,*u5);
	t = SubModSpe(*u1,*u5);
	*u1 = u;
	*u5 = t;
	u = AddModSpe(*u2,*u6);
	t = SubModSpe(*u2,*u6);
	*u2 = u;
	*u6 = t;
	u = AddModSpe(*u3,*u7);
	t = SubModSpe(*u3,*u7);
	*u3 = u;
	*u7 = t;
}
inline void small_butterflies0(sfixn* A, int next, int k){
	for(int j=0*next+k;j<1*next+k;j+=8){
		sfixn t,u;
		u = AddModSpe(*(A+j+0),*(A+j+0+2*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+2*next));
		*(A+j+0) = u;
		*(A+j+0+2*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+2*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+2*next));
		*(A+j+1) = u;
		*(A+j+1+2*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+2*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+2*next));
		*(A+j+2) = u;
		*(A+j+2+2*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+2*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+2*next));
		*(A+j+3) = u;
		*(A+j+3+2*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+2*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+2*next));
		*(A+j+4) = u;
		*(A+j+4+2*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+2*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+2*next));
		*(A+j+5) = u;
		*(A+j+5+2*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+2*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+2*next));
		*(A+j+6) = u;
		*(A+j+6+2*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+2*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+2*next));
		*(A+j+7) = u;
		*(A+j+7+2*next) = t;
	}
	for(int j=1*next+k;j<2*next+k;j+=8){
		sfixn t,u;
		*(A+j+0+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+0+2*next),MY_ROOT2);
		*(A+j+1+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+1+2*next),MY_ROOT2);
		*(A+j+2+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+2+2*next),MY_ROOT2);
		*(A+j+3+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+3+2*next),MY_ROOT2);
		*(A+j+4+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+4+2*next),MY_ROOT2);
		*(A+j+5+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+5+2*next),MY_ROOT2);
		*(A+j+6+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+6+2*next),MY_ROOT2);
		*(A+j+7+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+7+2*next),MY_ROOT2);
		u = AddModSpe(*(A+j+0),*(A+j+0+2*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+2*next));
		*(A+j+0) = u;
		*(A+j+0+2*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+2*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+2*next));
		*(A+j+1) = u;
		*(A+j+1+2*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+2*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+2*next));
		*(A+j+2) = u;
		*(A+j+2+2*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+2*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+2*next));
		*(A+j+3) = u;
		*(A+j+3+2*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+2*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+2*next));
		*(A+j+4) = u;
		*(A+j+4+2*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+2*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+2*next));
		*(A+j+5) = u;
		*(A+j+5+2*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+2*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+2*next));
		*(A+j+6) = u;
		*(A+j+6+2*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+2*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+2*next));
		*(A+j+7) = u;
		*(A+j+7+2*next) = t;
	}
	for(int j=4*next+k;j<5*next+k;j+=8){
		sfixn t,u;
		u = AddModSpe(*(A+j+0),*(A+j+0+2*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+2*next));
		*(A+j+0) = u;
		*(A+j+0+2*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+2*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+2*next));
		*(A+j+1) = u;
		*(A+j+1+2*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+2*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+2*next));
		*(A+j+2) = u;
		*(A+j+2+2*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+2*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+2*next));
		*(A+j+3) = u;
		*(A+j+3+2*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+2*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+2*next));
		*(A+j+4) = u;
		*(A+j+4+2*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+2*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+2*next));
		*(A+j+5) = u;
		*(A+j+5+2*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+2*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+2*next));
		*(A+j+6) = u;
		*(A+j+6+2*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+2*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+2*next));
		*(A+j+7) = u;
		*(A+j+7+2*next) = t;
	}
	for(int j=5*next+k;j<6*next+k;j+=8){
		sfixn t,u;
		*(A+j+0+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+0+2*next),MY_ROOT2);
		*(A+j+1+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+1+2*next),MY_ROOT2);
		*(A+j+2+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+2+2*next),MY_ROOT2);
		*(A+j+3+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+3+2*next),MY_ROOT2);
		*(A+j+4+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+4+2*next),MY_ROOT2);
		*(A+j+5+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+5+2*next),MY_ROOT2);
		*(A+j+6+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+6+2*next),MY_ROOT2);
		*(A+j+7+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+7+2*next),MY_ROOT2);
		u = AddModSpe(*(A+j+0),*(A+j+0+2*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+2*next));
		*(A+j+0) = u;
		*(A+j+0+2*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+2*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+2*next));
		*(A+j+1) = u;
		*(A+j+1+2*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+2*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+2*next));
		*(A+j+2) = u;
		*(A+j+2+2*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+2*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+2*next));
		*(A+j+3) = u;
		*(A+j+3+2*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+2*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+2*next));
		*(A+j+4) = u;
		*(A+j+4+2*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+2*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+2*next));
		*(A+j+5) = u;
		*(A+j+5+2*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+2*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+2*next));
		*(A+j+6) = u;
		*(A+j+6+2*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+2*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+2*next));
		*(A+j+7) = u;
		*(A+j+7+2*next) = t;
	}
	for(int j=0*next+k;j<1*next+k;j+=8){
		sfixn t,u;
		u = AddModSpe(*(A+j+0),*(A+j+0+4*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+4*next));
		*(A+j+0) = u;
		*(A+j+0+4*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+4*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+4*next));
		*(A+j+1) = u;
		*(A+j+1+4*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+4*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+4*next));
		*(A+j+2) = u;
		*(A+j+2+4*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+4*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+4*next));
		*(A+j+3) = u;
		*(A+j+3+4*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+4*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+4*next));
		*(A+j+4) = u;
		*(A+j+4+4*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+4*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+4*next));
		*(A+j+5) = u;
		*(A+j+5+4*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+4*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+4*next));
		*(A+j+6) = u;
		*(A+j+6+4*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+4*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+4*next));
		*(A+j+7) = u;
		*(A+j+7+4*next) = t;
	}
	for(int j=1*next+k;j<2*next+k;j+=8){
		sfixn t,u;
		*(A+j+0+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+0+4*next),MY_ROOT1);
		*(A+j+1+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+1+4*next),MY_ROOT1);
		*(A+j+2+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+2+4*next),MY_ROOT1);
		*(A+j+3+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+3+4*next),MY_ROOT1);
		*(A+j+4+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+4+4*next),MY_ROOT1);
		*(A+j+5+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+5+4*next),MY_ROOT1);
		*(A+j+6+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+6+4*next),MY_ROOT1);
		*(A+j+7+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+7+4*next),MY_ROOT1);
		u = AddModSpe(*(A+j+0),*(A+j+0+4*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+4*next));
		*(A+j+0) = u;
		*(A+j+0+4*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+4*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+4*next));
		*(A+j+1) = u;
		*(A+j+1+4*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+4*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+4*next));
		*(A+j+2) = u;
		*(A+j+2+4*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+4*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+4*next));
		*(A+j+3) = u;
		*(A+j+3+4*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+4*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+4*next));
		*(A+j+4) = u;
		*(A+j+4+4*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+4*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+4*next));
		*(A+j+5) = u;
		*(A+j+5+4*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+4*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+4*next));
		*(A+j+6) = u;
		*(A+j+6+4*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+4*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+4*next));
		*(A+j+7) = u;
		*(A+j+7+4*next) = t;
	}
	for(int j=2*next+k;j<3*next+k;j+=8){
		sfixn t,u;
		*(A+j+0+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+0+4*next),MY_ROOT2);
		*(A+j+1+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+1+4*next),MY_ROOT2);
		*(A+j+2+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+2+4*next),MY_ROOT2);
		*(A+j+3+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+3+4*next),MY_ROOT2);
		*(A+j+4+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+4+4*next),MY_ROOT2);
		*(A+j+5+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+5+4*next),MY_ROOT2);
		*(A+j+6+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+6+4*next),MY_ROOT2);
		*(A+j+7+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+7+4*next),MY_ROOT2);
		u = AddModSpe(*(A+j+0),*(A+j+0+4*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+4*next));
		*(A+j+0) = u;
		*(A+j+0+4*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+4*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+4*next));
		*(A+j+1) = u;
		*(A+j+1+4*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+4*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+4*next));
		*(A+j+2) = u;
		*(A+j+2+4*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+4*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+4*next));
		*(A+j+3) = u;
		*(A+j+3+4*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+4*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+4*next));
		*(A+j+4) = u;
		*(A+j+4+4*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+4*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+4*next));
		*(A+j+5) = u;
		*(A+j+5+4*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+4*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+4*next));
		*(A+j+6) = u;
		*(A+j+6+4*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+4*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+4*next));
		*(A+j+7) = u;
		*(A+j+7+4*next) = t;
	}
	for(int j=3*next+k;j<4*next+k;j+=8){
		sfixn t,u;
		*(A+j+0+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+0+4*next),MY_ROOT3);
		*(A+j+1+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+1+4*next),MY_ROOT3);
		*(A+j+2+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+2+4*next),MY_ROOT3);
		*(A+j+3+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+3+4*next),MY_ROOT3);
		*(A+j+4+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+4+4*next),MY_ROOT3);
		*(A+j+5+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+5+4*next),MY_ROOT3);
		*(A+j+6+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+6+4*next),MY_ROOT3);
		*(A+j+7+4*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+7+4*next),MY_ROOT3);
		u = AddModSpe(*(A+j+0),*(A+j+0+4*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+4*next));
		*(A+j+0) = u;
		*(A+j+0+4*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+4*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+4*next));
		*(A+j+1) = u;
		*(A+j+1+4*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+4*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+4*next));
		*(A+j+2) = u;
		*(A+j+2+4*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+4*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+4*next));
		*(A+j+3) = u;
		*(A+j+3+4*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+4*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+4*next));
		*(A+j+4) = u;
		*(A+j+4+4*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+4*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+4*next));
		*(A+j+5) = u;
		*(A+j+5+4*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+4*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+4*next));
		*(A+j+6) = u;
		*(A+j+6+4*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+4*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+4*next));
		*(A+j+7) = u;
		*(A+j+7+4*next) = t;
	}
}
inline void small_fft1(sfixn* u0, sfixn* u1){
	sfixn t,u;
	u = AddModSpe(*u0,*u1);
	t = SubModSpe(*u0,*u1);
	*u0 = u;
	*u1 = t;
}
inline void small_butterflies1(sfixn* A, int next, int k){
}
inline void small_fft2(sfixn* u0, sfixn* u1, sfixn* u2, sfixn* u3){
	sfixn t,u;
	u = AddModSpe(*u0,*u1);
	t = SubModSpe(*u0,*u1);
	*u0 = u;
	*u1 = t;
	u = AddModSpe(*u2,*u3);
	t = SubModSpe(*u2,*u3);
	*u2 = u;
	*u3 = t;
	*u2 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u2,MY_ROOT0);
	*u3 = MontMulModSpe_OPT3_AS_GENE_INLINE(*u3,MY_ROOT2);
	u = AddModSpe(*u0,*u2);
	t = SubModSpe(*u0,*u2);
	*u0 = u;
	*u2 = t;
	u = AddModSpe(*u1,*u3);
	t = SubModSpe(*u1,*u3);
	*u1 = u;
	*u3 = t;
}
inline void small_butterflies2(sfixn* A, int next, int k){
	for(int j=0*next+k;j<1*next+k;j+=8){
		sfixn t,u;
		u = AddModSpe(*(A+j+0),*(A+j+0+2*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+2*next));
		*(A+j+0) = u;
		*(A+j+0+2*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+2*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+2*next));
		*(A+j+1) = u;
		*(A+j+1+2*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+2*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+2*next));
		*(A+j+2) = u;
		*(A+j+2+2*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+2*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+2*next));
		*(A+j+3) = u;
		*(A+j+3+2*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+2*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+2*next));
		*(A+j+4) = u;
		*(A+j+4+2*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+2*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+2*next));
		*(A+j+5) = u;
		*(A+j+5+2*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+2*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+2*next));
		*(A+j+6) = u;
		*(A+j+6+2*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+2*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+2*next));
		*(A+j+7) = u;
		*(A+j+7+2*next) = t;
	}
	for(int j=1*next+k;j<2*next+k;j+=8){
		sfixn t,u;
		*(A+j+0+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+0+2*next),MY_ROOT2);
		*(A+j+1+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+1+2*next),MY_ROOT2);
		*(A+j+2+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+2+2*next),MY_ROOT2);
		*(A+j+3+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+3+2*next),MY_ROOT2);
		*(A+j+4+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+4+2*next),MY_ROOT2);
		*(A+j+5+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+5+2*next),MY_ROOT2);
		*(A+j+6+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+6+2*next),MY_ROOT2);
		*(A+j+7+2*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+7+2*next),MY_ROOT2);
		u = AddModSpe(*(A+j+0),*(A+j+0+2*next));
		t = SubModSpe(*(A+j+0),*(A+j+0+2*next));
		*(A+j+0) = u;
		*(A+j+0+2*next) = t;
		u = AddModSpe(*(A+j+1),*(A+j+1+2*next));
		t = SubModSpe(*(A+j+1),*(A+j+1+2*next));
		*(A+j+1) = u;
		*(A+j+1+2*next) = t;
		u = AddModSpe(*(A+j+2),*(A+j+2+2*next));
		t = SubModSpe(*(A+j+2),*(A+j+2+2*next));
		*(A+j+2) = u;
		*(A+j+2+2*next) = t;
		u = AddModSpe(*(A+j+3),*(A+j+3+2*next));
		t = SubModSpe(*(A+j+3),*(A+j+3+2*next));
		*(A+j+3) = u;
		*(A+j+3+2*next) = t;
		u = AddModSpe(*(A+j+4),*(A+j+4+2*next));
		t = SubModSpe(*(A+j+4),*(A+j+4+2*next));
		*(A+j+4) = u;
		*(A+j+4+2*next) = t;
		u = AddModSpe(*(A+j+5),*(A+j+5+2*next));
		t = SubModSpe(*(A+j+5),*(A+j+5+2*next));
		*(A+j+5) = u;
		*(A+j+5+2*next) = t;
		u = AddModSpe(*(A+j+6),*(A+j+6+2*next));
		t = SubModSpe(*(A+j+6),*(A+j+6+2*next));
		*(A+j+6) = u;
		*(A+j+6+2*next) = t;
		u = AddModSpe(*(A+j+7),*(A+j+7+2*next));
		t = SubModSpe(*(A+j+7),*(A+j+7+2*next));
		*(A+j+7) = u;
		*(A+j+7+2*next) = t;
	}
}
void DFT_iter(sfixn *A,sfixn *W){
	sfixn *Wp = W + 1600;
	sfixn* Wt;
	for(int k=0;k<1024;k+=MY_POW){
	small_fft0((A+k),(A+1+k),(A+2+k),(A+3+k),(A+4+k),(A+5+k),(A+6+k),(A+7+k));
	}
	unsigned long my_array[8]={0,4,2,6,1,5,3,7};
	Wp = Wp- 64;
	for(int k=0;k<1024;k+=64){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[8+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+1],*(Wp+start));
				start+=add;
				A[8+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+2],*(Wp+start));
				start+=add;
				A[8+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+3],*(Wp+start));
				start+=add;
				A[8+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+4],*(Wp+start));
				start+=add;
				A[8+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+5],*(Wp+start));
				start+=add;
				A[8+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+6],*(Wp+start));
				start+=add;
				A[8+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+8+k]);
				t=SubModSpe(A[0+k],A[0+8+k]);
				A[0+k]=u;
				A[0+8+k]=t;
				u=AddModSpe(A[1+k],A[1+8+k]);
				t=SubModSpe(A[1+k],A[1+8+k]);
				A[1+k]=u;
				A[1+8+k]=t;
				u=AddModSpe(A[2+k],A[2+8+k]);
				t=SubModSpe(A[2+k],A[2+8+k]);
				A[2+k]=u;
				A[2+8+k]=t;
				u=AddModSpe(A[3+k],A[3+8+k]);
				t=SubModSpe(A[3+k],A[3+8+k]);
				A[3+k]=u;
				A[3+8+k]=t;
				u=AddModSpe(A[4+k],A[4+8+k]);
				t=SubModSpe(A[4+k],A[4+8+k]);
				A[4+k]=u;
				A[4+8+k]=t;
				u=AddModSpe(A[5+k],A[5+8+k]);
				t=SubModSpe(A[5+k],A[5+8+k]);
				A[5+k]=u;
				A[5+8+k]=t;
				u=AddModSpe(A[6+k],A[6+8+k]);
				t=SubModSpe(A[6+k],A[6+8+k]);
				A[6+k]=u;
				A[6+8+k]=t;
				u=AddModSpe(A[7+k],A[7+8+k]);
				t=SubModSpe(A[7+k],A[7+8+k]);
				A[7+k]=u;
				A[7+8+k]=t;
			for(int i=8+k;i<8+k;i+=8){
				A[8+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i],*(Wp+start));
				start+=add;
				A[8+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+1],*(Wp+start));
				start+=add;
				A[8+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+2],*(Wp+start));
				start+=add;
				A[8+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+3],*(Wp+start));
				start+=add;
				A[8+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+4],*(Wp+start));
				start+=add;
				A[8+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+5],*(Wp+start));
				start+=add;
				A[8+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+6],*(Wp+start));
				start+=add;
				A[8+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+8]);
				t=SubModSpe(A[0+i],A[0+i+8]);
				A[0+i]=u;
				A[0+i+8]=t;
				u=AddModSpe(A[1+i],A[1+i+8]);
				t=SubModSpe(A[1+i],A[1+i+8]);
				A[1+i]=u;
				A[1+i+8]=t;
				u=AddModSpe(A[2+i],A[2+i+8]);
				t=SubModSpe(A[2+i],A[2+i+8]);
				A[2+i]=u;
				A[2+i+8]=t;
				u=AddModSpe(A[3+i],A[3+i+8]);
				t=SubModSpe(A[3+i],A[3+i+8]);
				A[3+i]=u;
				A[3+i+8]=t;
				u=AddModSpe(A[4+i],A[4+i+8]);
				t=SubModSpe(A[4+i],A[4+i+8]);
				A[4+i]=u;
				A[4+i+8]=t;
				u=AddModSpe(A[5+i],A[5+i+8]);
				t=SubModSpe(A[5+i],A[5+i+8]);
				A[5+i]=u;
				A[5+i+8]=t;
				u=AddModSpe(A[6+i],A[6+i+8]);
				t=SubModSpe(A[6+i],A[6+i+8]);
				A[6+i]=u;
				A[6+i+8]=t;
				u=AddModSpe(A[7+i],A[7+i+8]);
				t=SubModSpe(A[7+i],A[7+i+8]);
				A[7+i]=u;
				A[7+i+8]=t;
			}
		pos = 2;
		for(int j=k+2*8;j<k+64;j+=2*8,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+8],*(Wp+start2));
				start2+=add2;
				A[2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+8],*(Wp+start2));
				start2+=add2;
				A[3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+8],*(Wp+start2));
				start2+=add2;
				A[4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+8],*(Wp+start2));
				start2+=add2;
				A[5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+8],*(Wp+start2));
				start2+=add2;
				A[6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+8],*(Wp+start2));
				start2+=add2;
				A[7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+8]);
				t=SubModSpe(A[0+j],A[0+j+8]);
				A[0+j]=u;
				A[0+j+8]=t;
				u=AddModSpe(A[1+j],A[1+j+8]);
				t=SubModSpe(A[1+j],A[1+j+8]);
				A[1+j]=u;
				A[1+j+8]=t;
				u=AddModSpe(A[2+j],A[2+j+8]);
				t=SubModSpe(A[2+j],A[2+j+8]);
				A[2+j]=u;
				A[2+j+8]=t;
				u=AddModSpe(A[3+j],A[3+j+8]);
				t=SubModSpe(A[3+j],A[3+j+8]);
				A[3+j]=u;
				A[3+j+8]=t;
				u=AddModSpe(A[4+j],A[4+j+8]);
				t=SubModSpe(A[4+j],A[4+j+8]);
				A[4+j]=u;
				A[4+j+8]=t;
				u=AddModSpe(A[5+j],A[5+j+8]);
				t=SubModSpe(A[5+j],A[5+j+8]);
				A[5+j]=u;
				A[5+j+8]=t;
				u=AddModSpe(A[6+j],A[6+j+8]);
				t=SubModSpe(A[6+j],A[6+j+8]);
				A[6+j]=u;
				A[6+j+8]=t;
				u=AddModSpe(A[7+j],A[7+j+8]);
				t=SubModSpe(A[7+j],A[7+j+8]);
				A[7+j]=u;
				A[7+j+8]=t;
			for(int i=8;i<8;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+8],*(Wp+start2));
				start2+=add2;
				A[i+1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+8],*(Wp+start2));
				start2+=add2;
				A[i+2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+8],*(Wp+start2));
				start2+=add2;
				A[i+3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+8],*(Wp+start2));
				start2+=add2;
				A[i+4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+8],*(Wp+start2));
				start2+=add2;
				A[i+5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+8],*(Wp+start2));
				start2+=add2;
				A[i+6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+8],*(Wp+start2));
				start2+=add2;
				A[i+7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+8]);
				t=SubModSpe(A[0+i+j],A[0+j+i+8]);
				A[0+j+i]=u;
				A[0+j+i+8]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+8]);
				t=SubModSpe(A[1+i+j],A[1+j+i+8]);
				A[1+j+i]=u;
				A[1+j+i+8]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+8]);
				t=SubModSpe(A[2+i+j],A[2+j+i+8]);
				A[2+j+i]=u;
				A[2+j+i+8]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+8]);
				t=SubModSpe(A[3+i+j],A[3+j+i+8]);
				A[3+j+i]=u;
				A[3+j+i+8]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+8]);
				t=SubModSpe(A[4+i+j],A[4+j+i+8]);
				A[4+j+i]=u;
				A[4+j+i+8]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+8]);
				t=SubModSpe(A[5+i+j],A[5+j+i+8]);
				A[5+j+i]=u;
				A[5+j+i+8]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+8]);
				t=SubModSpe(A[6+i+j],A[6+j+i+8]);
				A[6+j+i]=u;
				A[6+j+i+8]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+8]);
				t=SubModSpe(A[7+i+j],A[7+j+i+8]);
				A[7+j+i]=u;
				A[7+j+i+8]=t;
			}
		}
		unsigned long next = 8;
		small_butterflies0(A,next,k);
	}
	Wp = Wp- 512;
	for(int k=0;k<1024;k+=512){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[64+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+1],*(Wp+start));
				start+=add;
				A[64+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+2],*(Wp+start));
				start+=add;
				A[64+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+3],*(Wp+start));
				start+=add;
				A[64+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+4],*(Wp+start));
				start+=add;
				A[64+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+5],*(Wp+start));
				start+=add;
				A[64+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+6],*(Wp+start));
				start+=add;
				A[64+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+64+k]);
				t=SubModSpe(A[0+k],A[0+64+k]);
				A[0+k]=u;
				A[0+64+k]=t;
				u=AddModSpe(A[1+k],A[1+64+k]);
				t=SubModSpe(A[1+k],A[1+64+k]);
				A[1+k]=u;
				A[1+64+k]=t;
				u=AddModSpe(A[2+k],A[2+64+k]);
				t=SubModSpe(A[2+k],A[2+64+k]);
				A[2+k]=u;
				A[2+64+k]=t;
				u=AddModSpe(A[3+k],A[3+64+k]);
				t=SubModSpe(A[3+k],A[3+64+k]);
				A[3+k]=u;
				A[3+64+k]=t;
				u=AddModSpe(A[4+k],A[4+64+k]);
				t=SubModSpe(A[4+k],A[4+64+k]);
				A[4+k]=u;
				A[4+64+k]=t;
				u=AddModSpe(A[5+k],A[5+64+k]);
				t=SubModSpe(A[5+k],A[5+64+k]);
				A[5+k]=u;
				A[5+64+k]=t;
				u=AddModSpe(A[6+k],A[6+64+k]);
				t=SubModSpe(A[6+k],A[6+64+k]);
				A[6+k]=u;
				A[6+64+k]=t;
				u=AddModSpe(A[7+k],A[7+64+k]);
				t=SubModSpe(A[7+k],A[7+64+k]);
				A[7+k]=u;
				A[7+64+k]=t;
			for(int i=8+k;i<64+k;i+=8){
				A[64+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i],*(Wp+start));
				start+=add;
				A[64+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+1],*(Wp+start));
				start+=add;
				A[64+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+2],*(Wp+start));
				start+=add;
				A[64+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+3],*(Wp+start));
				start+=add;
				A[64+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+4],*(Wp+start));
				start+=add;
				A[64+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+5],*(Wp+start));
				start+=add;
				A[64+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+6],*(Wp+start));
				start+=add;
				A[64+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+64]);
				t=SubModSpe(A[0+i],A[0+i+64]);
				A[0+i]=u;
				A[0+i+64]=t;
				u=AddModSpe(A[1+i],A[1+i+64]);
				t=SubModSpe(A[1+i],A[1+i+64]);
				A[1+i]=u;
				A[1+i+64]=t;
				u=AddModSpe(A[2+i],A[2+i+64]);
				t=SubModSpe(A[2+i],A[2+i+64]);
				A[2+i]=u;
				A[2+i+64]=t;
				u=AddModSpe(A[3+i],A[3+i+64]);
				t=SubModSpe(A[3+i],A[3+i+64]);
				A[3+i]=u;
				A[3+i+64]=t;
				u=AddModSpe(A[4+i],A[4+i+64]);
				t=SubModSpe(A[4+i],A[4+i+64]);
				A[4+i]=u;
				A[4+i+64]=t;
				u=AddModSpe(A[5+i],A[5+i+64]);
				t=SubModSpe(A[5+i],A[5+i+64]);
				A[5+i]=u;
				A[5+i+64]=t;
				u=AddModSpe(A[6+i],A[6+i+64]);
				t=SubModSpe(A[6+i],A[6+i+64]);
				A[6+i]=u;
				A[6+i+64]=t;
				u=AddModSpe(A[7+i],A[7+i+64]);
				t=SubModSpe(A[7+i],A[7+i+64]);
				A[7+i]=u;
				A[7+i+64]=t;
			}
		pos = 2;
		for(int j=k+2*64;j<k+512;j+=2*64,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+64],*(Wp+start2));
				start2+=add2;
				A[2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+64],*(Wp+start2));
				start2+=add2;
				A[3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+64],*(Wp+start2));
				start2+=add2;
				A[4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+64],*(Wp+start2));
				start2+=add2;
				A[5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+64],*(Wp+start2));
				start2+=add2;
				A[6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+64],*(Wp+start2));
				start2+=add2;
				A[7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+64],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+64]);
				t=SubModSpe(A[0+j],A[0+j+64]);
				A[0+j]=u;
				A[0+j+64]=t;
				u=AddModSpe(A[1+j],A[1+j+64]);
				t=SubModSpe(A[1+j],A[1+j+64]);
				A[1+j]=u;
				A[1+j+64]=t;
				u=AddModSpe(A[2+j],A[2+j+64]);
				t=SubModSpe(A[2+j],A[2+j+64]);
				A[2+j]=u;
				A[2+j+64]=t;
				u=AddModSpe(A[3+j],A[3+j+64]);
				t=SubModSpe(A[3+j],A[3+j+64]);
				A[3+j]=u;
				A[3+j+64]=t;
				u=AddModSpe(A[4+j],A[4+j+64]);
				t=SubModSpe(A[4+j],A[4+j+64]);
				A[4+j]=u;
				A[4+j+64]=t;
				u=AddModSpe(A[5+j],A[5+j+64]);
				t=SubModSpe(A[5+j],A[5+j+64]);
				A[5+j]=u;
				A[5+j+64]=t;
				u=AddModSpe(A[6+j],A[6+j+64]);
				t=SubModSpe(A[6+j],A[6+j+64]);
				A[6+j]=u;
				A[6+j+64]=t;
				u=AddModSpe(A[7+j],A[7+j+64]);
				t=SubModSpe(A[7+j],A[7+j+64]);
				A[7+j]=u;
				A[7+j+64]=t;
			for(int i=8;i<64;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+64],*(Wp+start2));
				start2+=add2;
				A[i+1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+64],*(Wp+start2));
				start2+=add2;
				A[i+2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+64],*(Wp+start2));
				start2+=add2;
				A[i+3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+64],*(Wp+start2));
				start2+=add2;
				A[i+4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+64],*(Wp+start2));
				start2+=add2;
				A[i+5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+64],*(Wp+start2));
				start2+=add2;
				A[i+6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+64],*(Wp+start2));
				start2+=add2;
				A[i+7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+64],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+64]);
				t=SubModSpe(A[0+i+j],A[0+j+i+64]);
				A[0+j+i]=u;
				A[0+j+i+64]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+64]);
				t=SubModSpe(A[1+i+j],A[1+j+i+64]);
				A[1+j+i]=u;
				A[1+j+i+64]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+64]);
				t=SubModSpe(A[2+i+j],A[2+j+i+64]);
				A[2+j+i]=u;
				A[2+j+i+64]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+64]);
				t=SubModSpe(A[3+i+j],A[3+j+i+64]);
				A[3+j+i]=u;
				A[3+j+i+64]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+64]);
				t=SubModSpe(A[4+i+j],A[4+j+i+64]);
				A[4+j+i]=u;
				A[4+j+i+64]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+64]);
				t=SubModSpe(A[5+i+j],A[5+j+i+64]);
				A[5+j+i]=u;
				A[5+j+i+64]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+64]);
				t=SubModSpe(A[6+i+j],A[6+j+i+64]);
				A[6+j+i]=u;
				A[6+j+i+64]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+64]);
				t=SubModSpe(A[7+i+j],A[7+j+i+64]);
				A[7+j+i]=u;
				A[7+j+i+64]=t;
			}
		}
		unsigned long next = 64;
		small_butterflies0(A,next,k);
	}
	Wp = Wp- 1024;
	unsigned long pos = 1;
	unsigned long add = my_array[pos]>>2;
	unsigned long start = add;
	A[512+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+1],*(Wp+start));
	start+=add;
	A[512+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+2],*(Wp+start));
	start+=add;
	A[512+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+3],*(Wp+start));
	start+=add;
	A[512+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+4],*(Wp+start));
	start+=add;
	A[512+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+5],*(Wp+start));
	start+=add;
	A[512+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+6],*(Wp+start));
	start+=add;
	A[512+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+7],*(Wp+start));
	start+=add;
	sfixn t,u;
	u=AddModSpe(A[0],A[0+512]);
	t=SubModSpe(A[0],A[0+512]);
	A[0]=u;
	A[0+512]=t;
	u=AddModSpe(A[1],A[1+512]);
	t=SubModSpe(A[1],A[1+512]);
	A[1]=u;
	A[1+512]=t;
	u=AddModSpe(A[2],A[2+512]);
	t=SubModSpe(A[2],A[2+512]);
	A[2]=u;
	A[2+512]=t;
	u=AddModSpe(A[3],A[3+512]);
	t=SubModSpe(A[3],A[3+512]);
	A[3]=u;
	A[3+512]=t;
	u=AddModSpe(A[4],A[4+512]);
	t=SubModSpe(A[4],A[4+512]);
	A[4]=u;
	A[4+512]=t;
	u=AddModSpe(A[5],A[5+512]);
	t=SubModSpe(A[5],A[5+512]);
	A[5]=u;
	A[5+512]=t;
	u=AddModSpe(A[6],A[6+512]);
	t=SubModSpe(A[6],A[6+512]);
	A[6]=u;
	A[6+512]=t;
	u=AddModSpe(A[7],A[7+512]);
	t=SubModSpe(A[7],A[7+512]);
	A[7]=u;
	A[7+512]=t;
	for(int i=8;i<512;i+=8){
		A[512+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i],*(Wp+start));
		start+=add;
		A[512+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+1],*(Wp+start));
		start+=add;
		A[512+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+2],*(Wp+start));
		start+=add;
		A[512+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+3],*(Wp+start));
		start+=add;
		A[512+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+4],*(Wp+start));
		start+=add;
		A[512+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+5],*(Wp+start));
		start+=add;
		A[512+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+6],*(Wp+start));
		start+=add;
		A[512+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[512+i+7],*(Wp+start));
		start+=add;
		u=AddModSpe(A[0+i],A[0+i+512]);
		t=SubModSpe(A[0+i],A[0+i+512]);
		A[0+i]=u;
		A[0+i+512]=t;
		u=AddModSpe(A[1+i],A[1+i+512]);
		t=SubModSpe(A[1+i],A[1+i+512]);
		A[1+i]=u;
		A[1+i+512]=t;
		u=AddModSpe(A[2+i],A[2+i+512]);
		t=SubModSpe(A[2+i],A[2+i+512]);
		A[2+i]=u;
		A[2+i+512]=t;
		u=AddModSpe(A[3+i],A[3+i+512]);
		t=SubModSpe(A[3+i],A[3+i+512]);
		A[3+i]=u;
		A[3+i+512]=t;
		u=AddModSpe(A[4+i],A[4+i+512]);
		t=SubModSpe(A[4+i],A[4+i+512]);
		A[4+i]=u;
		A[4+i+512]=t;
		u=AddModSpe(A[5+i],A[5+i+512]);
		t=SubModSpe(A[5+i],A[5+i+512]);
		A[5+i]=u;
		A[5+i+512]=t;
		u=AddModSpe(A[6+i],A[6+i+512]);
		t=SubModSpe(A[6+i],A[6+i+512]);
		A[6+i]=u;
		A[6+i+512]=t;
		u=AddModSpe(A[7+i],A[7+i+512]);
		t=SubModSpe(A[7+i],A[7+i+512]);
		A[7+i]=u;
		A[7+i+512]=t;
	}
	pos = 2;
	for(int j=2*512;j<1024;j+=2*512,pos+=2){
		unsigned long add1 = my_array[pos]>>2;
		unsigned long start1 = add1;
		A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
		start1+=add1;
		A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
		start1+=add1;
		A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
		start1+=add1;
		A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
		start1+=add1;
		A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
		start1+=add1;
		A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
		start1+=add1;
		A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
		start1+=add1;
		unsigned long add2 = my_array[pos+1]>>2;
		unsigned long start2 = add2;
		A[1+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+512],*(Wp+start2));
		start2+=add2;
		A[2+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+512],*(Wp+start2));
		start2+=add2;
		A[3+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+512],*(Wp+start2));
		start2+=add2;
		A[4+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+512],*(Wp+start2));
		start2+=add2;
		A[5+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+512],*(Wp+start2));
		start2+=add2;
		A[6+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+512],*(Wp+start2));
		start2+=add2;
		A[7+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+512],*(Wp+start2));
		start2+=add2;
		u=AddModSpe(A[0+j],A[0+j+512]);
		t=SubModSpe(A[0+j],A[0+j+512]);
		A[0+j]=u;
		A[0+j+512]=t;
		u=AddModSpe(A[1+j],A[1+j+512]);
		t=SubModSpe(A[1+j],A[1+j+512]);
		A[1+j]=u;
		A[1+j+512]=t;
		u=AddModSpe(A[2+j],A[2+j+512]);
		t=SubModSpe(A[2+j],A[2+j+512]);
		A[2+j]=u;
		A[2+j+512]=t;
		u=AddModSpe(A[3+j],A[3+j+512]);
		t=SubModSpe(A[3+j],A[3+j+512]);
		A[3+j]=u;
		A[3+j+512]=t;
		u=AddModSpe(A[4+j],A[4+j+512]);
		t=SubModSpe(A[4+j],A[4+j+512]);
		A[4+j]=u;
		A[4+j+512]=t;
		u=AddModSpe(A[5+j],A[5+j+512]);
		t=SubModSpe(A[5+j],A[5+j+512]);
		A[5+j]=u;
		A[5+j+512]=t;
		u=AddModSpe(A[6+j],A[6+j+512]);
		t=SubModSpe(A[6+j],A[6+j+512]);
		A[6+j]=u;
		A[6+j+512]=t;
		u=AddModSpe(A[7+j],A[7+j+512]);
		t=SubModSpe(A[7+j],A[7+j+512]);
		A[7+j]=u;
		A[7+j+512]=t;
		for(int i=8;i<512;i+=8){
			A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
			start1+=add1;
			A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
			start1+=add1;
			A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
			start1+=add1;
			A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
			start1+=add1;
			A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
			start1+=add1;
			A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
			start1+=add1;
			A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
			start1+=add1;
			A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
			start1+=add1;
			A[i+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+512],*(Wp+start2));
			start2+=add2;
			A[i+1+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+512],*(Wp+start2));
			start2+=add2;
			A[i+2+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+512],*(Wp+start2));
			start2+=add2;
			A[i+3+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+512],*(Wp+start2));
			start2+=add2;
			A[i+4+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+512],*(Wp+start2));
			start2+=add2;
			A[i+5+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+512],*(Wp+start2));
			start2+=add2;
			A[i+6+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+512],*(Wp+start2));
			start2+=add2;
			A[i+7+j+512] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+512],*(Wp+start2));
			start2+=add2;
			u=AddModSpe(A[0+i+j],A[0+j+i+512]);
			t=SubModSpe(A[0+i+j],A[0+j+i+512]);
			A[0+j+i]=u;
			A[0+j+i+512]=t;
			u=AddModSpe(A[1+i+j],A[1+j+i+512]);
			t=SubModSpe(A[1+i+j],A[1+j+i+512]);
			A[1+j+i]=u;
			A[1+j+i+512]=t;
			u=AddModSpe(A[2+i+j],A[2+j+i+512]);
			t=SubModSpe(A[2+i+j],A[2+j+i+512]);
			A[2+j+i]=u;
			A[2+j+i+512]=t;
			u=AddModSpe(A[3+i+j],A[3+j+i+512]);
			t=SubModSpe(A[3+i+j],A[3+j+i+512]);
			A[3+j+i]=u;
			A[3+j+i+512]=t;
			u=AddModSpe(A[4+i+j],A[4+j+i+512]);
			t=SubModSpe(A[4+i+j],A[4+j+i+512]);
			A[4+j+i]=u;
			A[4+j+i+512]=t;
			u=AddModSpe(A[5+i+j],A[5+j+i+512]);
			t=SubModSpe(A[5+i+j],A[5+j+i+512]);
			A[5+j+i]=u;
			A[5+j+i+512]=t;
			u=AddModSpe(A[6+i+j],A[6+j+i+512]);
			t=SubModSpe(A[6+i+j],A[6+j+i+512]);
			A[6+j+i]=u;
			A[6+j+i+512]=t;
			u=AddModSpe(A[7+i+j],A[7+j+i+512]);
			t=SubModSpe(A[7+i+j],A[7+j+i+512]);
			A[7+j+i]=u;
			A[7+j+i+512]=t;
		}
	}
	unsigned long next = 512;
	small_butterflies1(A,next,0);
}
void DFT_iter512(sfixn *A,sfixn *W){
	sfixn *Wp = W + 576;
	sfixn* Wt;
	for(int k=0;k<512;k+=MY_POW){
	small_fft0((A+k),(A+1+k),(A+2+k),(A+3+k),(A+4+k),(A+5+k),(A+6+k),(A+7+k));
	}
	unsigned long my_array[8]={0,4,2,6,1,5,3,7};
	Wp = Wp- 64;
	for(int k=0;k<512;k+=64){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[8+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+1],*(Wp+start));
				start+=add;
				A[8+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+2],*(Wp+start));
				start+=add;
				A[8+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+3],*(Wp+start));
				start+=add;
				A[8+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+4],*(Wp+start));
				start+=add;
				A[8+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+5],*(Wp+start));
				start+=add;
				A[8+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+6],*(Wp+start));
				start+=add;
				A[8+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+8+k]);
				t=SubModSpe(A[0+k],A[0+8+k]);
				A[0+k]=u;
				A[0+8+k]=t;
				u=AddModSpe(A[1+k],A[1+8+k]);
				t=SubModSpe(A[1+k],A[1+8+k]);
				A[1+k]=u;
				A[1+8+k]=t;
				u=AddModSpe(A[2+k],A[2+8+k]);
				t=SubModSpe(A[2+k],A[2+8+k]);
				A[2+k]=u;
				A[2+8+k]=t;
				u=AddModSpe(A[3+k],A[3+8+k]);
				t=SubModSpe(A[3+k],A[3+8+k]);
				A[3+k]=u;
				A[3+8+k]=t;
				u=AddModSpe(A[4+k],A[4+8+k]);
				t=SubModSpe(A[4+k],A[4+8+k]);
				A[4+k]=u;
				A[4+8+k]=t;
				u=AddModSpe(A[5+k],A[5+8+k]);
				t=SubModSpe(A[5+k],A[5+8+k]);
				A[5+k]=u;
				A[5+8+k]=t;
				u=AddModSpe(A[6+k],A[6+8+k]);
				t=SubModSpe(A[6+k],A[6+8+k]);
				A[6+k]=u;
				A[6+8+k]=t;
				u=AddModSpe(A[7+k],A[7+8+k]);
				t=SubModSpe(A[7+k],A[7+8+k]);
				A[7+k]=u;
				A[7+8+k]=t;
			for(int i=8+k;i<8+k;i+=8){
				A[8+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i],*(Wp+start));
				start+=add;
				A[8+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+1],*(Wp+start));
				start+=add;
				A[8+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+2],*(Wp+start));
				start+=add;
				A[8+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+3],*(Wp+start));
				start+=add;
				A[8+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+4],*(Wp+start));
				start+=add;
				A[8+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+5],*(Wp+start));
				start+=add;
				A[8+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+6],*(Wp+start));
				start+=add;
				A[8+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+8]);
				t=SubModSpe(A[0+i],A[0+i+8]);
				A[0+i]=u;
				A[0+i+8]=t;
				u=AddModSpe(A[1+i],A[1+i+8]);
				t=SubModSpe(A[1+i],A[1+i+8]);
				A[1+i]=u;
				A[1+i+8]=t;
				u=AddModSpe(A[2+i],A[2+i+8]);
				t=SubModSpe(A[2+i],A[2+i+8]);
				A[2+i]=u;
				A[2+i+8]=t;
				u=AddModSpe(A[3+i],A[3+i+8]);
				t=SubModSpe(A[3+i],A[3+i+8]);
				A[3+i]=u;
				A[3+i+8]=t;
				u=AddModSpe(A[4+i],A[4+i+8]);
				t=SubModSpe(A[4+i],A[4+i+8]);
				A[4+i]=u;
				A[4+i+8]=t;
				u=AddModSpe(A[5+i],A[5+i+8]);
				t=SubModSpe(A[5+i],A[5+i+8]);
				A[5+i]=u;
				A[5+i+8]=t;
				u=AddModSpe(A[6+i],A[6+i+8]);
				t=SubModSpe(A[6+i],A[6+i+8]);
				A[6+i]=u;
				A[6+i+8]=t;
				u=AddModSpe(A[7+i],A[7+i+8]);
				t=SubModSpe(A[7+i],A[7+i+8]);
				A[7+i]=u;
				A[7+i+8]=t;
			}
		pos = 2;
		for(int j=k+2*8;j<k+64;j+=2*8,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+8],*(Wp+start2));
				start2+=add2;
				A[2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+8],*(Wp+start2));
				start2+=add2;
				A[3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+8],*(Wp+start2));
				start2+=add2;
				A[4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+8],*(Wp+start2));
				start2+=add2;
				A[5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+8],*(Wp+start2));
				start2+=add2;
				A[6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+8],*(Wp+start2));
				start2+=add2;
				A[7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+8]);
				t=SubModSpe(A[0+j],A[0+j+8]);
				A[0+j]=u;
				A[0+j+8]=t;
				u=AddModSpe(A[1+j],A[1+j+8]);
				t=SubModSpe(A[1+j],A[1+j+8]);
				A[1+j]=u;
				A[1+j+8]=t;
				u=AddModSpe(A[2+j],A[2+j+8]);
				t=SubModSpe(A[2+j],A[2+j+8]);
				A[2+j]=u;
				A[2+j+8]=t;
				u=AddModSpe(A[3+j],A[3+j+8]);
				t=SubModSpe(A[3+j],A[3+j+8]);
				A[3+j]=u;
				A[3+j+8]=t;
				u=AddModSpe(A[4+j],A[4+j+8]);
				t=SubModSpe(A[4+j],A[4+j+8]);
				A[4+j]=u;
				A[4+j+8]=t;
				u=AddModSpe(A[5+j],A[5+j+8]);
				t=SubModSpe(A[5+j],A[5+j+8]);
				A[5+j]=u;
				A[5+j+8]=t;
				u=AddModSpe(A[6+j],A[6+j+8]);
				t=SubModSpe(A[6+j],A[6+j+8]);
				A[6+j]=u;
				A[6+j+8]=t;
				u=AddModSpe(A[7+j],A[7+j+8]);
				t=SubModSpe(A[7+j],A[7+j+8]);
				A[7+j]=u;
				A[7+j+8]=t;
			for(int i=8;i<8;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+8],*(Wp+start2));
				start2+=add2;
				A[i+1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+8],*(Wp+start2));
				start2+=add2;
				A[i+2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+8],*(Wp+start2));
				start2+=add2;
				A[i+3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+8],*(Wp+start2));
				start2+=add2;
				A[i+4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+8],*(Wp+start2));
				start2+=add2;
				A[i+5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+8],*(Wp+start2));
				start2+=add2;
				A[i+6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+8],*(Wp+start2));
				start2+=add2;
				A[i+7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+8]);
				t=SubModSpe(A[0+i+j],A[0+j+i+8]);
				A[0+j+i]=u;
				A[0+j+i+8]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+8]);
				t=SubModSpe(A[1+i+j],A[1+j+i+8]);
				A[1+j+i]=u;
				A[1+j+i+8]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+8]);
				t=SubModSpe(A[2+i+j],A[2+j+i+8]);
				A[2+j+i]=u;
				A[2+j+i+8]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+8]);
				t=SubModSpe(A[3+i+j],A[3+j+i+8]);
				A[3+j+i]=u;
				A[3+j+i+8]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+8]);
				t=SubModSpe(A[4+i+j],A[4+j+i+8]);
				A[4+j+i]=u;
				A[4+j+i+8]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+8]);
				t=SubModSpe(A[5+i+j],A[5+j+i+8]);
				A[5+j+i]=u;
				A[5+j+i+8]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+8]);
				t=SubModSpe(A[6+i+j],A[6+j+i+8]);
				A[6+j+i]=u;
				A[6+j+i+8]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+8]);
				t=SubModSpe(A[7+i+j],A[7+j+i+8]);
				A[7+j+i]=u;
				A[7+j+i+8]=t;
			}
		}
		unsigned long next = 8;
		small_butterflies0(A,next,k);
	}
	Wp = Wp- 512;
	for(int k=0;k<512;k+=512){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[64+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+1],*(Wp+start));
				start+=add;
				A[64+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+2],*(Wp+start));
				start+=add;
				A[64+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+3],*(Wp+start));
				start+=add;
				A[64+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+4],*(Wp+start));
				start+=add;
				A[64+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+5],*(Wp+start));
				start+=add;
				A[64+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+6],*(Wp+start));
				start+=add;
				A[64+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+64+k]);
				t=SubModSpe(A[0+k],A[0+64+k]);
				A[0+k]=u;
				A[0+64+k]=t;
				u=AddModSpe(A[1+k],A[1+64+k]);
				t=SubModSpe(A[1+k],A[1+64+k]);
				A[1+k]=u;
				A[1+64+k]=t;
				u=AddModSpe(A[2+k],A[2+64+k]);
				t=SubModSpe(A[2+k],A[2+64+k]);
				A[2+k]=u;
				A[2+64+k]=t;
				u=AddModSpe(A[3+k],A[3+64+k]);
				t=SubModSpe(A[3+k],A[3+64+k]);
				A[3+k]=u;
				A[3+64+k]=t;
				u=AddModSpe(A[4+k],A[4+64+k]);
				t=SubModSpe(A[4+k],A[4+64+k]);
				A[4+k]=u;
				A[4+64+k]=t;
				u=AddModSpe(A[5+k],A[5+64+k]);
				t=SubModSpe(A[5+k],A[5+64+k]);
				A[5+k]=u;
				A[5+64+k]=t;
				u=AddModSpe(A[6+k],A[6+64+k]);
				t=SubModSpe(A[6+k],A[6+64+k]);
				A[6+k]=u;
				A[6+64+k]=t;
				u=AddModSpe(A[7+k],A[7+64+k]);
				t=SubModSpe(A[7+k],A[7+64+k]);
				A[7+k]=u;
				A[7+64+k]=t;
			for(int i=8+k;i<64+k;i+=8){
				A[64+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i],*(Wp+start));
				start+=add;
				A[64+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+1],*(Wp+start));
				start+=add;
				A[64+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+2],*(Wp+start));
				start+=add;
				A[64+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+3],*(Wp+start));
				start+=add;
				A[64+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+4],*(Wp+start));
				start+=add;
				A[64+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+5],*(Wp+start));
				start+=add;
				A[64+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+6],*(Wp+start));
				start+=add;
				A[64+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+64]);
				t=SubModSpe(A[0+i],A[0+i+64]);
				A[0+i]=u;
				A[0+i+64]=t;
				u=AddModSpe(A[1+i],A[1+i+64]);
				t=SubModSpe(A[1+i],A[1+i+64]);
				A[1+i]=u;
				A[1+i+64]=t;
				u=AddModSpe(A[2+i],A[2+i+64]);
				t=SubModSpe(A[2+i],A[2+i+64]);
				A[2+i]=u;
				A[2+i+64]=t;
				u=AddModSpe(A[3+i],A[3+i+64]);
				t=SubModSpe(A[3+i],A[3+i+64]);
				A[3+i]=u;
				A[3+i+64]=t;
				u=AddModSpe(A[4+i],A[4+i+64]);
				t=SubModSpe(A[4+i],A[4+i+64]);
				A[4+i]=u;
				A[4+i+64]=t;
				u=AddModSpe(A[5+i],A[5+i+64]);
				t=SubModSpe(A[5+i],A[5+i+64]);
				A[5+i]=u;
				A[5+i+64]=t;
				u=AddModSpe(A[6+i],A[6+i+64]);
				t=SubModSpe(A[6+i],A[6+i+64]);
				A[6+i]=u;
				A[6+i+64]=t;
				u=AddModSpe(A[7+i],A[7+i+64]);
				t=SubModSpe(A[7+i],A[7+i+64]);
				A[7+i]=u;
				A[7+i+64]=t;
			}
		pos = 2;
		for(int j=k+2*64;j<k+512;j+=2*64,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+64],*(Wp+start2));
				start2+=add2;
				A[2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+64],*(Wp+start2));
				start2+=add2;
				A[3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+64],*(Wp+start2));
				start2+=add2;
				A[4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+64],*(Wp+start2));
				start2+=add2;
				A[5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+64],*(Wp+start2));
				start2+=add2;
				A[6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+64],*(Wp+start2));
				start2+=add2;
				A[7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+64],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+64]);
				t=SubModSpe(A[0+j],A[0+j+64]);
				A[0+j]=u;
				A[0+j+64]=t;
				u=AddModSpe(A[1+j],A[1+j+64]);
				t=SubModSpe(A[1+j],A[1+j+64]);
				A[1+j]=u;
				A[1+j+64]=t;
				u=AddModSpe(A[2+j],A[2+j+64]);
				t=SubModSpe(A[2+j],A[2+j+64]);
				A[2+j]=u;
				A[2+j+64]=t;
				u=AddModSpe(A[3+j],A[3+j+64]);
				t=SubModSpe(A[3+j],A[3+j+64]);
				A[3+j]=u;
				A[3+j+64]=t;
				u=AddModSpe(A[4+j],A[4+j+64]);
				t=SubModSpe(A[4+j],A[4+j+64]);
				A[4+j]=u;
				A[4+j+64]=t;
				u=AddModSpe(A[5+j],A[5+j+64]);
				t=SubModSpe(A[5+j],A[5+j+64]);
				A[5+j]=u;
				A[5+j+64]=t;
				u=AddModSpe(A[6+j],A[6+j+64]);
				t=SubModSpe(A[6+j],A[6+j+64]);
				A[6+j]=u;
				A[6+j+64]=t;
				u=AddModSpe(A[7+j],A[7+j+64]);
				t=SubModSpe(A[7+j],A[7+j+64]);
				A[7+j]=u;
				A[7+j+64]=t;
			for(int i=8;i<64;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+64],*(Wp+start2));
				start2+=add2;
				A[i+1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+64],*(Wp+start2));
				start2+=add2;
				A[i+2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+64],*(Wp+start2));
				start2+=add2;
				A[i+3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+64],*(Wp+start2));
				start2+=add2;
				A[i+4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+64],*(Wp+start2));
				start2+=add2;
				A[i+5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+64],*(Wp+start2));
				start2+=add2;
				A[i+6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+64],*(Wp+start2));
				start2+=add2;
				A[i+7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+64],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+64]);
				t=SubModSpe(A[0+i+j],A[0+j+i+64]);
				A[0+j+i]=u;
				A[0+j+i+64]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+64]);
				t=SubModSpe(A[1+i+j],A[1+j+i+64]);
				A[1+j+i]=u;
				A[1+j+i+64]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+64]);
				t=SubModSpe(A[2+i+j],A[2+j+i+64]);
				A[2+j+i]=u;
				A[2+j+i+64]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+64]);
				t=SubModSpe(A[3+i+j],A[3+j+i+64]);
				A[3+j+i]=u;
				A[3+j+i+64]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+64]);
				t=SubModSpe(A[4+i+j],A[4+j+i+64]);
				A[4+j+i]=u;
				A[4+j+i+64]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+64]);
				t=SubModSpe(A[5+i+j],A[5+j+i+64]);
				A[5+j+i]=u;
				A[5+j+i+64]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+64]);
				t=SubModSpe(A[6+i+j],A[6+j+i+64]);
				A[6+j+i]=u;
				A[6+j+i+64]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+64]);
				t=SubModSpe(A[7+i+j],A[7+j+i+64]);
				A[7+j+i]=u;
				A[7+j+i+64]=t;
			}
		}
		unsigned long next = 64;
		small_butterflies0(A,next,k);
	}
}
void DFT_iter256(sfixn *A,sfixn *W){
	sfixn *Wp = W + 320;
	sfixn* Wt;
	for(int k=0;k<256;k+=MY_POW){
	small_fft0((A+k),(A+1+k),(A+2+k),(A+3+k),(A+4+k),(A+5+k),(A+6+k),(A+7+k));
	}
	unsigned long my_array[8]={0,4,2,6,1,5,3,7};
	Wp = Wp- 64;
	for(int k=0;k<256;k+=64){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[8+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+1],*(Wp+start));
				start+=add;
				A[8+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+2],*(Wp+start));
				start+=add;
				A[8+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+3],*(Wp+start));
				start+=add;
				A[8+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+4],*(Wp+start));
				start+=add;
				A[8+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+5],*(Wp+start));
				start+=add;
				A[8+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+6],*(Wp+start));
				start+=add;
				A[8+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+8+k]);
				t=SubModSpe(A[0+k],A[0+8+k]);
				A[0+k]=u;
				A[0+8+k]=t;
				u=AddModSpe(A[1+k],A[1+8+k]);
				t=SubModSpe(A[1+k],A[1+8+k]);
				A[1+k]=u;
				A[1+8+k]=t;
				u=AddModSpe(A[2+k],A[2+8+k]);
				t=SubModSpe(A[2+k],A[2+8+k]);
				A[2+k]=u;
				A[2+8+k]=t;
				u=AddModSpe(A[3+k],A[3+8+k]);
				t=SubModSpe(A[3+k],A[3+8+k]);
				A[3+k]=u;
				A[3+8+k]=t;
				u=AddModSpe(A[4+k],A[4+8+k]);
				t=SubModSpe(A[4+k],A[4+8+k]);
				A[4+k]=u;
				A[4+8+k]=t;
				u=AddModSpe(A[5+k],A[5+8+k]);
				t=SubModSpe(A[5+k],A[5+8+k]);
				A[5+k]=u;
				A[5+8+k]=t;
				u=AddModSpe(A[6+k],A[6+8+k]);
				t=SubModSpe(A[6+k],A[6+8+k]);
				A[6+k]=u;
				A[6+8+k]=t;
				u=AddModSpe(A[7+k],A[7+8+k]);
				t=SubModSpe(A[7+k],A[7+8+k]);
				A[7+k]=u;
				A[7+8+k]=t;
			for(int i=8+k;i<8+k;i+=8){
				A[8+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i],*(Wp+start));
				start+=add;
				A[8+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+1],*(Wp+start));
				start+=add;
				A[8+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+2],*(Wp+start));
				start+=add;
				A[8+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+3],*(Wp+start));
				start+=add;
				A[8+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+4],*(Wp+start));
				start+=add;
				A[8+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+5],*(Wp+start));
				start+=add;
				A[8+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+6],*(Wp+start));
				start+=add;
				A[8+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+8]);
				t=SubModSpe(A[0+i],A[0+i+8]);
				A[0+i]=u;
				A[0+i+8]=t;
				u=AddModSpe(A[1+i],A[1+i+8]);
				t=SubModSpe(A[1+i],A[1+i+8]);
				A[1+i]=u;
				A[1+i+8]=t;
				u=AddModSpe(A[2+i],A[2+i+8]);
				t=SubModSpe(A[2+i],A[2+i+8]);
				A[2+i]=u;
				A[2+i+8]=t;
				u=AddModSpe(A[3+i],A[3+i+8]);
				t=SubModSpe(A[3+i],A[3+i+8]);
				A[3+i]=u;
				A[3+i+8]=t;
				u=AddModSpe(A[4+i],A[4+i+8]);
				t=SubModSpe(A[4+i],A[4+i+8]);
				A[4+i]=u;
				A[4+i+8]=t;
				u=AddModSpe(A[5+i],A[5+i+8]);
				t=SubModSpe(A[5+i],A[5+i+8]);
				A[5+i]=u;
				A[5+i+8]=t;
				u=AddModSpe(A[6+i],A[6+i+8]);
				t=SubModSpe(A[6+i],A[6+i+8]);
				A[6+i]=u;
				A[6+i+8]=t;
				u=AddModSpe(A[7+i],A[7+i+8]);
				t=SubModSpe(A[7+i],A[7+i+8]);
				A[7+i]=u;
				A[7+i+8]=t;
			}
		pos = 2;
		for(int j=k+2*8;j<k+64;j+=2*8,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+8],*(Wp+start2));
				start2+=add2;
				A[2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+8],*(Wp+start2));
				start2+=add2;
				A[3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+8],*(Wp+start2));
				start2+=add2;
				A[4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+8],*(Wp+start2));
				start2+=add2;
				A[5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+8],*(Wp+start2));
				start2+=add2;
				A[6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+8],*(Wp+start2));
				start2+=add2;
				A[7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+8]);
				t=SubModSpe(A[0+j],A[0+j+8]);
				A[0+j]=u;
				A[0+j+8]=t;
				u=AddModSpe(A[1+j],A[1+j+8]);
				t=SubModSpe(A[1+j],A[1+j+8]);
				A[1+j]=u;
				A[1+j+8]=t;
				u=AddModSpe(A[2+j],A[2+j+8]);
				t=SubModSpe(A[2+j],A[2+j+8]);
				A[2+j]=u;
				A[2+j+8]=t;
				u=AddModSpe(A[3+j],A[3+j+8]);
				t=SubModSpe(A[3+j],A[3+j+8]);
				A[3+j]=u;
				A[3+j+8]=t;
				u=AddModSpe(A[4+j],A[4+j+8]);
				t=SubModSpe(A[4+j],A[4+j+8]);
				A[4+j]=u;
				A[4+j+8]=t;
				u=AddModSpe(A[5+j],A[5+j+8]);
				t=SubModSpe(A[5+j],A[5+j+8]);
				A[5+j]=u;
				A[5+j+8]=t;
				u=AddModSpe(A[6+j],A[6+j+8]);
				t=SubModSpe(A[6+j],A[6+j+8]);
				A[6+j]=u;
				A[6+j+8]=t;
				u=AddModSpe(A[7+j],A[7+j+8]);
				t=SubModSpe(A[7+j],A[7+j+8]);
				A[7+j]=u;
				A[7+j+8]=t;
			for(int i=8;i<8;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+8],*(Wp+start2));
				start2+=add2;
				A[i+1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+8],*(Wp+start2));
				start2+=add2;
				A[i+2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+8],*(Wp+start2));
				start2+=add2;
				A[i+3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+8],*(Wp+start2));
				start2+=add2;
				A[i+4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+8],*(Wp+start2));
				start2+=add2;
				A[i+5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+8],*(Wp+start2));
				start2+=add2;
				A[i+6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+8],*(Wp+start2));
				start2+=add2;
				A[i+7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+8]);
				t=SubModSpe(A[0+i+j],A[0+j+i+8]);
				A[0+j+i]=u;
				A[0+j+i+8]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+8]);
				t=SubModSpe(A[1+i+j],A[1+j+i+8]);
				A[1+j+i]=u;
				A[1+j+i+8]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+8]);
				t=SubModSpe(A[2+i+j],A[2+j+i+8]);
				A[2+j+i]=u;
				A[2+j+i+8]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+8]);
				t=SubModSpe(A[3+i+j],A[3+j+i+8]);
				A[3+j+i]=u;
				A[3+j+i+8]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+8]);
				t=SubModSpe(A[4+i+j],A[4+j+i+8]);
				A[4+j+i]=u;
				A[4+j+i+8]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+8]);
				t=SubModSpe(A[5+i+j],A[5+j+i+8]);
				A[5+j+i]=u;
				A[5+j+i+8]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+8]);
				t=SubModSpe(A[6+i+j],A[6+j+i+8]);
				A[6+j+i]=u;
				A[6+j+i+8]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+8]);
				t=SubModSpe(A[7+i+j],A[7+j+i+8]);
				A[7+j+i]=u;
				A[7+j+i+8]=t;
			}
		}
		unsigned long next = 8;
		small_butterflies0(A,next,k);
	}
	Wp = Wp- 256;
	unsigned long pos = 1;
	unsigned long add = my_array[pos]>>1;
	unsigned long start = add;
	A[64+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+1],*(Wp+start));
	start+=add;
	A[64+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+2],*(Wp+start));
	start+=add;
	A[64+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+3],*(Wp+start));
	start+=add;
	A[64+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+4],*(Wp+start));
	start+=add;
	A[64+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+5],*(Wp+start));
	start+=add;
	A[64+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+6],*(Wp+start));
	start+=add;
	A[64+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+7],*(Wp+start));
	start+=add;
	sfixn t,u;
	u=AddModSpe(A[0],A[0+64]);
	t=SubModSpe(A[0],A[0+64]);
	A[0]=u;
	A[0+64]=t;
	u=AddModSpe(A[1],A[1+64]);
	t=SubModSpe(A[1],A[1+64]);
	A[1]=u;
	A[1+64]=t;
	u=AddModSpe(A[2],A[2+64]);
	t=SubModSpe(A[2],A[2+64]);
	A[2]=u;
	A[2+64]=t;
	u=AddModSpe(A[3],A[3+64]);
	t=SubModSpe(A[3],A[3+64]);
	A[3]=u;
	A[3+64]=t;
	u=AddModSpe(A[4],A[4+64]);
	t=SubModSpe(A[4],A[4+64]);
	A[4]=u;
	A[4+64]=t;
	u=AddModSpe(A[5],A[5+64]);
	t=SubModSpe(A[5],A[5+64]);
	A[5]=u;
	A[5+64]=t;
	u=AddModSpe(A[6],A[6+64]);
	t=SubModSpe(A[6],A[6+64]);
	A[6]=u;
	A[6+64]=t;
	u=AddModSpe(A[7],A[7+64]);
	t=SubModSpe(A[7],A[7+64]);
	A[7]=u;
	A[7+64]=t;
	for(int i=8;i<64;i+=8){
		A[64+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i],*(Wp+start));
		start+=add;
		A[64+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+1],*(Wp+start));
		start+=add;
		A[64+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+2],*(Wp+start));
		start+=add;
		A[64+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+3],*(Wp+start));
		start+=add;
		A[64+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+4],*(Wp+start));
		start+=add;
		A[64+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+5],*(Wp+start));
		start+=add;
		A[64+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+6],*(Wp+start));
		start+=add;
		A[64+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+7],*(Wp+start));
		start+=add;
		u=AddModSpe(A[0+i],A[0+i+64]);
		t=SubModSpe(A[0+i],A[0+i+64]);
		A[0+i]=u;
		A[0+i+64]=t;
		u=AddModSpe(A[1+i],A[1+i+64]);
		t=SubModSpe(A[1+i],A[1+i+64]);
		A[1+i]=u;
		A[1+i+64]=t;
		u=AddModSpe(A[2+i],A[2+i+64]);
		t=SubModSpe(A[2+i],A[2+i+64]);
		A[2+i]=u;
		A[2+i+64]=t;
		u=AddModSpe(A[3+i],A[3+i+64]);
		t=SubModSpe(A[3+i],A[3+i+64]);
		A[3+i]=u;
		A[3+i+64]=t;
		u=AddModSpe(A[4+i],A[4+i+64]);
		t=SubModSpe(A[4+i],A[4+i+64]);
		A[4+i]=u;
		A[4+i+64]=t;
		u=AddModSpe(A[5+i],A[5+i+64]);
		t=SubModSpe(A[5+i],A[5+i+64]);
		A[5+i]=u;
		A[5+i+64]=t;
		u=AddModSpe(A[6+i],A[6+i+64]);
		t=SubModSpe(A[6+i],A[6+i+64]);
		A[6+i]=u;
		A[6+i+64]=t;
		u=AddModSpe(A[7+i],A[7+i+64]);
		t=SubModSpe(A[7+i],A[7+i+64]);
		A[7+i]=u;
		A[7+i+64]=t;
	}
	pos = 2;
	for(int j=2*64;j<256;j+=2*64,pos+=2){
		unsigned long add1 = my_array[pos]>>1;
		unsigned long start1 = add1;
		A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
		start1+=add1;
		A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
		start1+=add1;
		A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
		start1+=add1;
		A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
		start1+=add1;
		A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
		start1+=add1;
		A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
		start1+=add1;
		A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
		start1+=add1;
		unsigned long add2 = my_array[pos+1]>>1;
		unsigned long start2 = add2;
		A[1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+64],*(Wp+start2));
		start2+=add2;
		A[2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+64],*(Wp+start2));
		start2+=add2;
		A[3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+64],*(Wp+start2));
		start2+=add2;
		A[4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+64],*(Wp+start2));
		start2+=add2;
		A[5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+64],*(Wp+start2));
		start2+=add2;
		A[6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+64],*(Wp+start2));
		start2+=add2;
		A[7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+64],*(Wp+start2));
		start2+=add2;
		u=AddModSpe(A[0+j],A[0+j+64]);
		t=SubModSpe(A[0+j],A[0+j+64]);
		A[0+j]=u;
		A[0+j+64]=t;
		u=AddModSpe(A[1+j],A[1+j+64]);
		t=SubModSpe(A[1+j],A[1+j+64]);
		A[1+j]=u;
		A[1+j+64]=t;
		u=AddModSpe(A[2+j],A[2+j+64]);
		t=SubModSpe(A[2+j],A[2+j+64]);
		A[2+j]=u;
		A[2+j+64]=t;
		u=AddModSpe(A[3+j],A[3+j+64]);
		t=SubModSpe(A[3+j],A[3+j+64]);
		A[3+j]=u;
		A[3+j+64]=t;
		u=AddModSpe(A[4+j],A[4+j+64]);
		t=SubModSpe(A[4+j],A[4+j+64]);
		A[4+j]=u;
		A[4+j+64]=t;
		u=AddModSpe(A[5+j],A[5+j+64]);
		t=SubModSpe(A[5+j],A[5+j+64]);
		A[5+j]=u;
		A[5+j+64]=t;
		u=AddModSpe(A[6+j],A[6+j+64]);
		t=SubModSpe(A[6+j],A[6+j+64]);
		A[6+j]=u;
		A[6+j+64]=t;
		u=AddModSpe(A[7+j],A[7+j+64]);
		t=SubModSpe(A[7+j],A[7+j+64]);
		A[7+j]=u;
		A[7+j+64]=t;
		for(int i=8;i<64;i+=8){
			A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
			start1+=add1;
			A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
			start1+=add1;
			A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
			start1+=add1;
			A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
			start1+=add1;
			A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
			start1+=add1;
			A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
			start1+=add1;
			A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
			start1+=add1;
			A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
			start1+=add1;
			A[i+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+64],*(Wp+start2));
			start2+=add2;
			A[i+1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+64],*(Wp+start2));
			start2+=add2;
			A[i+2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+64],*(Wp+start2));
			start2+=add2;
			A[i+3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+64],*(Wp+start2));
			start2+=add2;
			A[i+4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+64],*(Wp+start2));
			start2+=add2;
			A[i+5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+64],*(Wp+start2));
			start2+=add2;
			A[i+6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+64],*(Wp+start2));
			start2+=add2;
			A[i+7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+64],*(Wp+start2));
			start2+=add2;
			u=AddModSpe(A[0+i+j],A[0+j+i+64]);
			t=SubModSpe(A[0+i+j],A[0+j+i+64]);
			A[0+j+i]=u;
			A[0+j+i+64]=t;
			u=AddModSpe(A[1+i+j],A[1+j+i+64]);
			t=SubModSpe(A[1+i+j],A[1+j+i+64]);
			A[1+j+i]=u;
			A[1+j+i+64]=t;
			u=AddModSpe(A[2+i+j],A[2+j+i+64]);
			t=SubModSpe(A[2+i+j],A[2+j+i+64]);
			A[2+j+i]=u;
			A[2+j+i+64]=t;
			u=AddModSpe(A[3+i+j],A[3+j+i+64]);
			t=SubModSpe(A[3+i+j],A[3+j+i+64]);
			A[3+j+i]=u;
			A[3+j+i+64]=t;
			u=AddModSpe(A[4+i+j],A[4+j+i+64]);
			t=SubModSpe(A[4+i+j],A[4+j+i+64]);
			A[4+j+i]=u;
			A[4+j+i+64]=t;
			u=AddModSpe(A[5+i+j],A[5+j+i+64]);
			t=SubModSpe(A[5+i+j],A[5+j+i+64]);
			A[5+j+i]=u;
			A[5+j+i+64]=t;
			u=AddModSpe(A[6+i+j],A[6+j+i+64]);
			t=SubModSpe(A[6+i+j],A[6+j+i+64]);
			A[6+j+i]=u;
			A[6+j+i+64]=t;
			u=AddModSpe(A[7+i+j],A[7+j+i+64]);
			t=SubModSpe(A[7+i+j],A[7+j+i+64]);
			A[7+j+i]=u;
			A[7+j+i+64]=t;
		}
	}
	unsigned long next = 64;
	small_butterflies2(A,next,0);
}
void DFT_iter128(sfixn *A,sfixn *W){
	sfixn *Wp = W + 192;
	sfixn* Wt;
	for(int k=0;k<128;k+=MY_POW){
	small_fft0((A+k),(A+1+k),(A+2+k),(A+3+k),(A+4+k),(A+5+k),(A+6+k),(A+7+k));
	}
	unsigned long my_array[8]={0,4,2,6,1,5,3,7};
	Wp = Wp- 64;
	for(int k=0;k<128;k+=64){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[8+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+1],*(Wp+start));
				start+=add;
				A[8+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+2],*(Wp+start));
				start+=add;
				A[8+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+3],*(Wp+start));
				start+=add;
				A[8+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+4],*(Wp+start));
				start+=add;
				A[8+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+5],*(Wp+start));
				start+=add;
				A[8+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+6],*(Wp+start));
				start+=add;
				A[8+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+8+k]);
				t=SubModSpe(A[0+k],A[0+8+k]);
				A[0+k]=u;
				A[0+8+k]=t;
				u=AddModSpe(A[1+k],A[1+8+k]);
				t=SubModSpe(A[1+k],A[1+8+k]);
				A[1+k]=u;
				A[1+8+k]=t;
				u=AddModSpe(A[2+k],A[2+8+k]);
				t=SubModSpe(A[2+k],A[2+8+k]);
				A[2+k]=u;
				A[2+8+k]=t;
				u=AddModSpe(A[3+k],A[3+8+k]);
				t=SubModSpe(A[3+k],A[3+8+k]);
				A[3+k]=u;
				A[3+8+k]=t;
				u=AddModSpe(A[4+k],A[4+8+k]);
				t=SubModSpe(A[4+k],A[4+8+k]);
				A[4+k]=u;
				A[4+8+k]=t;
				u=AddModSpe(A[5+k],A[5+8+k]);
				t=SubModSpe(A[5+k],A[5+8+k]);
				A[5+k]=u;
				A[5+8+k]=t;
				u=AddModSpe(A[6+k],A[6+8+k]);
				t=SubModSpe(A[6+k],A[6+8+k]);
				A[6+k]=u;
				A[6+8+k]=t;
				u=AddModSpe(A[7+k],A[7+8+k]);
				t=SubModSpe(A[7+k],A[7+8+k]);
				A[7+k]=u;
				A[7+8+k]=t;
			for(int i=8+k;i<8+k;i+=8){
				A[8+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i],*(Wp+start));
				start+=add;
				A[8+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+1],*(Wp+start));
				start+=add;
				A[8+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+2],*(Wp+start));
				start+=add;
				A[8+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+3],*(Wp+start));
				start+=add;
				A[8+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+4],*(Wp+start));
				start+=add;
				A[8+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+5],*(Wp+start));
				start+=add;
				A[8+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+6],*(Wp+start));
				start+=add;
				A[8+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+8]);
				t=SubModSpe(A[0+i],A[0+i+8]);
				A[0+i]=u;
				A[0+i+8]=t;
				u=AddModSpe(A[1+i],A[1+i+8]);
				t=SubModSpe(A[1+i],A[1+i+8]);
				A[1+i]=u;
				A[1+i+8]=t;
				u=AddModSpe(A[2+i],A[2+i+8]);
				t=SubModSpe(A[2+i],A[2+i+8]);
				A[2+i]=u;
				A[2+i+8]=t;
				u=AddModSpe(A[3+i],A[3+i+8]);
				t=SubModSpe(A[3+i],A[3+i+8]);
				A[3+i]=u;
				A[3+i+8]=t;
				u=AddModSpe(A[4+i],A[4+i+8]);
				t=SubModSpe(A[4+i],A[4+i+8]);
				A[4+i]=u;
				A[4+i+8]=t;
				u=AddModSpe(A[5+i],A[5+i+8]);
				t=SubModSpe(A[5+i],A[5+i+8]);
				A[5+i]=u;
				A[5+i+8]=t;
				u=AddModSpe(A[6+i],A[6+i+8]);
				t=SubModSpe(A[6+i],A[6+i+8]);
				A[6+i]=u;
				A[6+i+8]=t;
				u=AddModSpe(A[7+i],A[7+i+8]);
				t=SubModSpe(A[7+i],A[7+i+8]);
				A[7+i]=u;
				A[7+i+8]=t;
			}
		pos = 2;
		for(int j=k+2*8;j<k+64;j+=2*8,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+8],*(Wp+start2));
				start2+=add2;
				A[2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+8],*(Wp+start2));
				start2+=add2;
				A[3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+8],*(Wp+start2));
				start2+=add2;
				A[4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+8],*(Wp+start2));
				start2+=add2;
				A[5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+8],*(Wp+start2));
				start2+=add2;
				A[6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+8],*(Wp+start2));
				start2+=add2;
				A[7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+8]);
				t=SubModSpe(A[0+j],A[0+j+8]);
				A[0+j]=u;
				A[0+j+8]=t;
				u=AddModSpe(A[1+j],A[1+j+8]);
				t=SubModSpe(A[1+j],A[1+j+8]);
				A[1+j]=u;
				A[1+j+8]=t;
				u=AddModSpe(A[2+j],A[2+j+8]);
				t=SubModSpe(A[2+j],A[2+j+8]);
				A[2+j]=u;
				A[2+j+8]=t;
				u=AddModSpe(A[3+j],A[3+j+8]);
				t=SubModSpe(A[3+j],A[3+j+8]);
				A[3+j]=u;
				A[3+j+8]=t;
				u=AddModSpe(A[4+j],A[4+j+8]);
				t=SubModSpe(A[4+j],A[4+j+8]);
				A[4+j]=u;
				A[4+j+8]=t;
				u=AddModSpe(A[5+j],A[5+j+8]);
				t=SubModSpe(A[5+j],A[5+j+8]);
				A[5+j]=u;
				A[5+j+8]=t;
				u=AddModSpe(A[6+j],A[6+j+8]);
				t=SubModSpe(A[6+j],A[6+j+8]);
				A[6+j]=u;
				A[6+j+8]=t;
				u=AddModSpe(A[7+j],A[7+j+8]);
				t=SubModSpe(A[7+j],A[7+j+8]);
				A[7+j]=u;
				A[7+j+8]=t;
			for(int i=8;i<8;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+8],*(Wp+start2));
				start2+=add2;
				A[i+1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+8],*(Wp+start2));
				start2+=add2;
				A[i+2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+8],*(Wp+start2));
				start2+=add2;
				A[i+3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+8],*(Wp+start2));
				start2+=add2;
				A[i+4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+8],*(Wp+start2));
				start2+=add2;
				A[i+5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+8],*(Wp+start2));
				start2+=add2;
				A[i+6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+8],*(Wp+start2));
				start2+=add2;
				A[i+7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+8]);
				t=SubModSpe(A[0+i+j],A[0+j+i+8]);
				A[0+j+i]=u;
				A[0+j+i+8]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+8]);
				t=SubModSpe(A[1+i+j],A[1+j+i+8]);
				A[1+j+i]=u;
				A[1+j+i+8]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+8]);
				t=SubModSpe(A[2+i+j],A[2+j+i+8]);
				A[2+j+i]=u;
				A[2+j+i+8]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+8]);
				t=SubModSpe(A[3+i+j],A[3+j+i+8]);
				A[3+j+i]=u;
				A[3+j+i+8]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+8]);
				t=SubModSpe(A[4+i+j],A[4+j+i+8]);
				A[4+j+i]=u;
				A[4+j+i+8]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+8]);
				t=SubModSpe(A[5+i+j],A[5+j+i+8]);
				A[5+j+i]=u;
				A[5+j+i+8]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+8]);
				t=SubModSpe(A[6+i+j],A[6+j+i+8]);
				A[6+j+i]=u;
				A[6+j+i+8]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+8]);
				t=SubModSpe(A[7+i+j],A[7+j+i+8]);
				A[7+j+i]=u;
				A[7+j+i+8]=t;
			}
		}
		unsigned long next = 8;
		small_butterflies0(A,next,k);
	}
	Wp = Wp- 128;
	unsigned long pos = 1;
	unsigned long add = my_array[pos]>>2;
	unsigned long start = add;
	A[64+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+1],*(Wp+start));
	start+=add;
	A[64+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+2],*(Wp+start));
	start+=add;
	A[64+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+3],*(Wp+start));
	start+=add;
	A[64+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+4],*(Wp+start));
	start+=add;
	A[64+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+5],*(Wp+start));
	start+=add;
	A[64+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+6],*(Wp+start));
	start+=add;
	A[64+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+7],*(Wp+start));
	start+=add;
	sfixn t,u;
	u=AddModSpe(A[0],A[0+64]);
	t=SubModSpe(A[0],A[0+64]);
	A[0]=u;
	A[0+64]=t;
	u=AddModSpe(A[1],A[1+64]);
	t=SubModSpe(A[1],A[1+64]);
	A[1]=u;
	A[1+64]=t;
	u=AddModSpe(A[2],A[2+64]);
	t=SubModSpe(A[2],A[2+64]);
	A[2]=u;
	A[2+64]=t;
	u=AddModSpe(A[3],A[3+64]);
	t=SubModSpe(A[3],A[3+64]);
	A[3]=u;
	A[3+64]=t;
	u=AddModSpe(A[4],A[4+64]);
	t=SubModSpe(A[4],A[4+64]);
	A[4]=u;
	A[4+64]=t;
	u=AddModSpe(A[5],A[5+64]);
	t=SubModSpe(A[5],A[5+64]);
	A[5]=u;
	A[5+64]=t;
	u=AddModSpe(A[6],A[6+64]);
	t=SubModSpe(A[6],A[6+64]);
	A[6]=u;
	A[6+64]=t;
	u=AddModSpe(A[7],A[7+64]);
	t=SubModSpe(A[7],A[7+64]);
	A[7]=u;
	A[7+64]=t;
	for(int i=8;i<64;i+=8){
		A[64+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i],*(Wp+start));
		start+=add;
		A[64+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+1],*(Wp+start));
		start+=add;
		A[64+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+2],*(Wp+start));
		start+=add;
		A[64+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+3],*(Wp+start));
		start+=add;
		A[64+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+4],*(Wp+start));
		start+=add;
		A[64+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+5],*(Wp+start));
		start+=add;
		A[64+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+6],*(Wp+start));
		start+=add;
		A[64+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[64+i+7],*(Wp+start));
		start+=add;
		u=AddModSpe(A[0+i],A[0+i+64]);
		t=SubModSpe(A[0+i],A[0+i+64]);
		A[0+i]=u;
		A[0+i+64]=t;
		u=AddModSpe(A[1+i],A[1+i+64]);
		t=SubModSpe(A[1+i],A[1+i+64]);
		A[1+i]=u;
		A[1+i+64]=t;
		u=AddModSpe(A[2+i],A[2+i+64]);
		t=SubModSpe(A[2+i],A[2+i+64]);
		A[2+i]=u;
		A[2+i+64]=t;
		u=AddModSpe(A[3+i],A[3+i+64]);
		t=SubModSpe(A[3+i],A[3+i+64]);
		A[3+i]=u;
		A[3+i+64]=t;
		u=AddModSpe(A[4+i],A[4+i+64]);
		t=SubModSpe(A[4+i],A[4+i+64]);
		A[4+i]=u;
		A[4+i+64]=t;
		u=AddModSpe(A[5+i],A[5+i+64]);
		t=SubModSpe(A[5+i],A[5+i+64]);
		A[5+i]=u;
		A[5+i+64]=t;
		u=AddModSpe(A[6+i],A[6+i+64]);
		t=SubModSpe(A[6+i],A[6+i+64]);
		A[6+i]=u;
		A[6+i+64]=t;
		u=AddModSpe(A[7+i],A[7+i+64]);
		t=SubModSpe(A[7+i],A[7+i+64]);
		A[7+i]=u;
		A[7+i+64]=t;
	}
	pos = 2;
	for(int j=2*64;j<128;j+=2*64,pos+=2){
		unsigned long add1 = my_array[pos]>>2;
		unsigned long start1 = add1;
		A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
		start1+=add1;
		A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
		start1+=add1;
		A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
		start1+=add1;
		A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
		start1+=add1;
		A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
		start1+=add1;
		A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
		start1+=add1;
		A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
		start1+=add1;
		unsigned long add2 = my_array[pos+1]>>2;
		unsigned long start2 = add2;
		A[1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+64],*(Wp+start2));
		start2+=add2;
		A[2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+64],*(Wp+start2));
		start2+=add2;
		A[3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+64],*(Wp+start2));
		start2+=add2;
		A[4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+64],*(Wp+start2));
		start2+=add2;
		A[5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+64],*(Wp+start2));
		start2+=add2;
		A[6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+64],*(Wp+start2));
		start2+=add2;
		A[7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+64],*(Wp+start2));
		start2+=add2;
		u=AddModSpe(A[0+j],A[0+j+64]);
		t=SubModSpe(A[0+j],A[0+j+64]);
		A[0+j]=u;
		A[0+j+64]=t;
		u=AddModSpe(A[1+j],A[1+j+64]);
		t=SubModSpe(A[1+j],A[1+j+64]);
		A[1+j]=u;
		A[1+j+64]=t;
		u=AddModSpe(A[2+j],A[2+j+64]);
		t=SubModSpe(A[2+j],A[2+j+64]);
		A[2+j]=u;
		A[2+j+64]=t;
		u=AddModSpe(A[3+j],A[3+j+64]);
		t=SubModSpe(A[3+j],A[3+j+64]);
		A[3+j]=u;
		A[3+j+64]=t;
		u=AddModSpe(A[4+j],A[4+j+64]);
		t=SubModSpe(A[4+j],A[4+j+64]);
		A[4+j]=u;
		A[4+j+64]=t;
		u=AddModSpe(A[5+j],A[5+j+64]);
		t=SubModSpe(A[5+j],A[5+j+64]);
		A[5+j]=u;
		A[5+j+64]=t;
		u=AddModSpe(A[6+j],A[6+j+64]);
		t=SubModSpe(A[6+j],A[6+j+64]);
		A[6+j]=u;
		A[6+j+64]=t;
		u=AddModSpe(A[7+j],A[7+j+64]);
		t=SubModSpe(A[7+j],A[7+j+64]);
		A[7+j]=u;
		A[7+j+64]=t;
		for(int i=8;i<64;i+=8){
			A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
			start1+=add1;
			A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
			start1+=add1;
			A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
			start1+=add1;
			A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
			start1+=add1;
			A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
			start1+=add1;
			A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
			start1+=add1;
			A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
			start1+=add1;
			A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
			start1+=add1;
			A[i+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+64],*(Wp+start2));
			start2+=add2;
			A[i+1+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+64],*(Wp+start2));
			start2+=add2;
			A[i+2+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+64],*(Wp+start2));
			start2+=add2;
			A[i+3+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+64],*(Wp+start2));
			start2+=add2;
			A[i+4+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+64],*(Wp+start2));
			start2+=add2;
			A[i+5+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+64],*(Wp+start2));
			start2+=add2;
			A[i+6+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+64],*(Wp+start2));
			start2+=add2;
			A[i+7+j+64] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+64],*(Wp+start2));
			start2+=add2;
			u=AddModSpe(A[0+i+j],A[0+j+i+64]);
			t=SubModSpe(A[0+i+j],A[0+j+i+64]);
			A[0+j+i]=u;
			A[0+j+i+64]=t;
			u=AddModSpe(A[1+i+j],A[1+j+i+64]);
			t=SubModSpe(A[1+i+j],A[1+j+i+64]);
			A[1+j+i]=u;
			A[1+j+i+64]=t;
			u=AddModSpe(A[2+i+j],A[2+j+i+64]);
			t=SubModSpe(A[2+i+j],A[2+j+i+64]);
			A[2+j+i]=u;
			A[2+j+i+64]=t;
			u=AddModSpe(A[3+i+j],A[3+j+i+64]);
			t=SubModSpe(A[3+i+j],A[3+j+i+64]);
			A[3+j+i]=u;
			A[3+j+i+64]=t;
			u=AddModSpe(A[4+i+j],A[4+j+i+64]);
			t=SubModSpe(A[4+i+j],A[4+j+i+64]);
			A[4+j+i]=u;
			A[4+j+i+64]=t;
			u=AddModSpe(A[5+i+j],A[5+j+i+64]);
			t=SubModSpe(A[5+i+j],A[5+j+i+64]);
			A[5+j+i]=u;
			A[5+j+i+64]=t;
			u=AddModSpe(A[6+i+j],A[6+j+i+64]);
			t=SubModSpe(A[6+i+j],A[6+j+i+64]);
			A[6+j+i]=u;
			A[6+j+i+64]=t;
			u=AddModSpe(A[7+i+j],A[7+j+i+64]);
			t=SubModSpe(A[7+i+j],A[7+j+i+64]);
			A[7+j+i]=u;
			A[7+j+i+64]=t;
		}
	}
	unsigned long next = 64;
	small_butterflies1(A,next,0);
}
void DFT_iter64(sfixn *A,sfixn *W){
	sfixn *Wp = W + 64;
	sfixn* Wt;
	for(int k=0;k<64;k+=MY_POW){
	small_fft0((A+k),(A+1+k),(A+2+k),(A+3+k),(A+4+k),(A+5+k),(A+6+k),(A+7+k));
	}
	unsigned long my_array[8]={0,4,2,6,1,5,3,7};
	Wp = Wp- 64;
	for(int k=0;k<64;k+=64){
		unsigned long pos = 1;
			unsigned long add = my_array[pos];
			unsigned long start = add;
				A[8+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+1],*(Wp+start));
				start+=add;
				A[8+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+2],*(Wp+start));
				start+=add;
				A[8+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+3],*(Wp+start));
				start+=add;
				A[8+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+4],*(Wp+start));
				start+=add;
				A[8+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+5],*(Wp+start));
				start+=add;
				A[8+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+6],*(Wp+start));
				start+=add;
				A[8+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+k+7],*(Wp+start));
				start+=add;
				sfixn t,u;
				u=AddModSpe(A[0+k],A[0+8+k]);
				t=SubModSpe(A[0+k],A[0+8+k]);
				A[0+k]=u;
				A[0+8+k]=t;
				u=AddModSpe(A[1+k],A[1+8+k]);
				t=SubModSpe(A[1+k],A[1+8+k]);
				A[1+k]=u;
				A[1+8+k]=t;
				u=AddModSpe(A[2+k],A[2+8+k]);
				t=SubModSpe(A[2+k],A[2+8+k]);
				A[2+k]=u;
				A[2+8+k]=t;
				u=AddModSpe(A[3+k],A[3+8+k]);
				t=SubModSpe(A[3+k],A[3+8+k]);
				A[3+k]=u;
				A[3+8+k]=t;
				u=AddModSpe(A[4+k],A[4+8+k]);
				t=SubModSpe(A[4+k],A[4+8+k]);
				A[4+k]=u;
				A[4+8+k]=t;
				u=AddModSpe(A[5+k],A[5+8+k]);
				t=SubModSpe(A[5+k],A[5+8+k]);
				A[5+k]=u;
				A[5+8+k]=t;
				u=AddModSpe(A[6+k],A[6+8+k]);
				t=SubModSpe(A[6+k],A[6+8+k]);
				A[6+k]=u;
				A[6+8+k]=t;
				u=AddModSpe(A[7+k],A[7+8+k]);
				t=SubModSpe(A[7+k],A[7+8+k]);
				A[7+k]=u;
				A[7+8+k]=t;
			for(int i=8+k;i<8+k;i+=8){
				A[8+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i],*(Wp+start));
				start+=add;
				A[8+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+1],*(Wp+start));
				start+=add;
				A[8+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+2],*(Wp+start));
				start+=add;
				A[8+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+3],*(Wp+start));
				start+=add;
				A[8+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+4],*(Wp+start));
				start+=add;
				A[8+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+5],*(Wp+start));
				start+=add;
				A[8+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+6],*(Wp+start));
				start+=add;
				A[8+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+7],*(Wp+start));
				start+=add;
				u=AddModSpe(A[0+i],A[0+i+8]);
				t=SubModSpe(A[0+i],A[0+i+8]);
				A[0+i]=u;
				A[0+i+8]=t;
				u=AddModSpe(A[1+i],A[1+i+8]);
				t=SubModSpe(A[1+i],A[1+i+8]);
				A[1+i]=u;
				A[1+i+8]=t;
				u=AddModSpe(A[2+i],A[2+i+8]);
				t=SubModSpe(A[2+i],A[2+i+8]);
				A[2+i]=u;
				A[2+i+8]=t;
				u=AddModSpe(A[3+i],A[3+i+8]);
				t=SubModSpe(A[3+i],A[3+i+8]);
				A[3+i]=u;
				A[3+i+8]=t;
				u=AddModSpe(A[4+i],A[4+i+8]);
				t=SubModSpe(A[4+i],A[4+i+8]);
				A[4+i]=u;
				A[4+i+8]=t;
				u=AddModSpe(A[5+i],A[5+i+8]);
				t=SubModSpe(A[5+i],A[5+i+8]);
				A[5+i]=u;
				A[5+i+8]=t;
				u=AddModSpe(A[6+i],A[6+i+8]);
				t=SubModSpe(A[6+i],A[6+i+8]);
				A[6+i]=u;
				A[6+i+8]=t;
				u=AddModSpe(A[7+i],A[7+i+8]);
				t=SubModSpe(A[7+i],A[7+i+8]);
				A[7+i]=u;
				A[7+i+8]=t;
			}
		pos = 2;
		for(int j=k+2*8;j<k+64;j+=2*8,pos+=2){
			unsigned long add1 = my_array[pos];
			unsigned long start1 = add1;
				A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
				start1+=add1;
				A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
				start1+=add1;
				A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
				start1+=add1;
				A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
				start1+=add1;
				A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
				start1+=add1;
				A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
				start1+=add1;
				A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
				start1+=add1;
			unsigned long add2 = my_array[pos+1];
			unsigned long start2 = add2;
				A[1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+8],*(Wp+start2));
				start2+=add2;
				A[2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+8],*(Wp+start2));
				start2+=add2;
				A[3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+8],*(Wp+start2));
				start2+=add2;
				A[4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+8],*(Wp+start2));
				start2+=add2;
				A[5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+8],*(Wp+start2));
				start2+=add2;
				A[6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+8],*(Wp+start2));
				start2+=add2;
				A[7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+j],A[0+j+8]);
				t=SubModSpe(A[0+j],A[0+j+8]);
				A[0+j]=u;
				A[0+j+8]=t;
				u=AddModSpe(A[1+j],A[1+j+8]);
				t=SubModSpe(A[1+j],A[1+j+8]);
				A[1+j]=u;
				A[1+j+8]=t;
				u=AddModSpe(A[2+j],A[2+j+8]);
				t=SubModSpe(A[2+j],A[2+j+8]);
				A[2+j]=u;
				A[2+j+8]=t;
				u=AddModSpe(A[3+j],A[3+j+8]);
				t=SubModSpe(A[3+j],A[3+j+8]);
				A[3+j]=u;
				A[3+j+8]=t;
				u=AddModSpe(A[4+j],A[4+j+8]);
				t=SubModSpe(A[4+j],A[4+j+8]);
				A[4+j]=u;
				A[4+j+8]=t;
				u=AddModSpe(A[5+j],A[5+j+8]);
				t=SubModSpe(A[5+j],A[5+j+8]);
				A[5+j]=u;
				A[5+j+8]=t;
				u=AddModSpe(A[6+j],A[6+j+8]);
				t=SubModSpe(A[6+j],A[6+j+8]);
				A[6+j]=u;
				A[6+j+8]=t;
				u=AddModSpe(A[7+j],A[7+j+8]);
				t=SubModSpe(A[7+j],A[7+j+8]);
				A[7+j]=u;
				A[7+j+8]=t;
			for(int i=8;i<8;i+=8){
				A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
				start1+=add1;
				A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
				start1+=add1;
				A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
				start1+=add1;
				A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
				start1+=add1;
				A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
				start1+=add1;
				A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
				start1+=add1;
				A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
				start1+=add1;
				A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
				start1+=add1;
				A[i+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+8],*(Wp+start2));
				start2+=add2;
				A[i+1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+8],*(Wp+start2));
				start2+=add2;
				A[i+2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+8],*(Wp+start2));
				start2+=add2;
				A[i+3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+8],*(Wp+start2));
				start2+=add2;
				A[i+4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+8],*(Wp+start2));
				start2+=add2;
				A[i+5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+8],*(Wp+start2));
				start2+=add2;
				A[i+6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+8],*(Wp+start2));
				start2+=add2;
				A[i+7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+8],*(Wp+start2));
				start2+=add2;
				u=AddModSpe(A[0+i+j],A[0+j+i+8]);
				t=SubModSpe(A[0+i+j],A[0+j+i+8]);
				A[0+j+i]=u;
				A[0+j+i+8]=t;
				u=AddModSpe(A[1+i+j],A[1+j+i+8]);
				t=SubModSpe(A[1+i+j],A[1+j+i+8]);
				A[1+j+i]=u;
				A[1+j+i+8]=t;
				u=AddModSpe(A[2+i+j],A[2+j+i+8]);
				t=SubModSpe(A[2+i+j],A[2+j+i+8]);
				A[2+j+i]=u;
				A[2+j+i+8]=t;
				u=AddModSpe(A[3+i+j],A[3+j+i+8]);
				t=SubModSpe(A[3+i+j],A[3+j+i+8]);
				A[3+j+i]=u;
				A[3+j+i+8]=t;
				u=AddModSpe(A[4+i+j],A[4+j+i+8]);
				t=SubModSpe(A[4+i+j],A[4+j+i+8]);
				A[4+j+i]=u;
				A[4+j+i+8]=t;
				u=AddModSpe(A[5+i+j],A[5+j+i+8]);
				t=SubModSpe(A[5+i+j],A[5+j+i+8]);
				A[5+j+i]=u;
				A[5+j+i+8]=t;
				u=AddModSpe(A[6+i+j],A[6+j+i+8]);
				t=SubModSpe(A[6+i+j],A[6+j+i+8]);
				A[6+j+i]=u;
				A[6+j+i+8]=t;
				u=AddModSpe(A[7+i+j],A[7+j+i+8]);
				t=SubModSpe(A[7+i+j],A[7+j+i+8]);
				A[7+j+i]=u;
				A[7+j+i+8]=t;
			}
		}
		unsigned long next = 8;
		small_butterflies0(A,next,k);
	}
}
void DFT_iter32(sfixn *A,sfixn *W){
	sfixn *Wp = W + 32;
	sfixn* Wt;
	for(int k=0;k<32;k+=MY_POW){
	small_fft0((A+k),(A+1+k),(A+2+k),(A+3+k),(A+4+k),(A+5+k),(A+6+k),(A+7+k));
	}
	unsigned long my_array[8]={0,4,2,6,1,5,3,7};
	Wp = Wp- 32;
	unsigned long pos = 1;
	unsigned long add = my_array[pos]>>1;
	unsigned long start = add;
	A[8+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+1],*(Wp+start));
	start+=add;
	A[8+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+2],*(Wp+start));
	start+=add;
	A[8+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+3],*(Wp+start));
	start+=add;
	A[8+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+4],*(Wp+start));
	start+=add;
	A[8+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+5],*(Wp+start));
	start+=add;
	A[8+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+6],*(Wp+start));
	start+=add;
	A[8+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+7],*(Wp+start));
	start+=add;
	sfixn t,u;
	u=AddModSpe(A[0],A[0+8]);
	t=SubModSpe(A[0],A[0+8]);
	A[0]=u;
	A[0+8]=t;
	u=AddModSpe(A[1],A[1+8]);
	t=SubModSpe(A[1],A[1+8]);
	A[1]=u;
	A[1+8]=t;
	u=AddModSpe(A[2],A[2+8]);
	t=SubModSpe(A[2],A[2+8]);
	A[2]=u;
	A[2+8]=t;
	u=AddModSpe(A[3],A[3+8]);
	t=SubModSpe(A[3],A[3+8]);
	A[3]=u;
	A[3+8]=t;
	u=AddModSpe(A[4],A[4+8]);
	t=SubModSpe(A[4],A[4+8]);
	A[4]=u;
	A[4+8]=t;
	u=AddModSpe(A[5],A[5+8]);
	t=SubModSpe(A[5],A[5+8]);
	A[5]=u;
	A[5+8]=t;
	u=AddModSpe(A[6],A[6+8]);
	t=SubModSpe(A[6],A[6+8]);
	A[6]=u;
	A[6+8]=t;
	u=AddModSpe(A[7],A[7+8]);
	t=SubModSpe(A[7],A[7+8]);
	A[7]=u;
	A[7+8]=t;
	for(int i=8;i<8;i+=8){
		A[8+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i],*(Wp+start));
		start+=add;
		A[8+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+1],*(Wp+start));
		start+=add;
		A[8+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+2],*(Wp+start));
		start+=add;
		A[8+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+3],*(Wp+start));
		start+=add;
		A[8+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+4],*(Wp+start));
		start+=add;
		A[8+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+5],*(Wp+start));
		start+=add;
		A[8+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+6],*(Wp+start));
		start+=add;
		A[8+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[8+i+7],*(Wp+start));
		start+=add;
		u=AddModSpe(A[0+i],A[0+i+8]);
		t=SubModSpe(A[0+i],A[0+i+8]);
		A[0+i]=u;
		A[0+i+8]=t;
		u=AddModSpe(A[1+i],A[1+i+8]);
		t=SubModSpe(A[1+i],A[1+i+8]);
		A[1+i]=u;
		A[1+i+8]=t;
		u=AddModSpe(A[2+i],A[2+i+8]);
		t=SubModSpe(A[2+i],A[2+i+8]);
		A[2+i]=u;
		A[2+i+8]=t;
		u=AddModSpe(A[3+i],A[3+i+8]);
		t=SubModSpe(A[3+i],A[3+i+8]);
		A[3+i]=u;
		A[3+i+8]=t;
		u=AddModSpe(A[4+i],A[4+i+8]);
		t=SubModSpe(A[4+i],A[4+i+8]);
		A[4+i]=u;
		A[4+i+8]=t;
		u=AddModSpe(A[5+i],A[5+i+8]);
		t=SubModSpe(A[5+i],A[5+i+8]);
		A[5+i]=u;
		A[5+i+8]=t;
		u=AddModSpe(A[6+i],A[6+i+8]);
		t=SubModSpe(A[6+i],A[6+i+8]);
		A[6+i]=u;
		A[6+i+8]=t;
		u=AddModSpe(A[7+i],A[7+i+8]);
		t=SubModSpe(A[7+i],A[7+i+8]);
		A[7+i]=u;
		A[7+i+8]=t;
	}
	pos = 2;
	for(int j=2*8;j<32;j+=2*8,pos+=2){
		unsigned long add1 = my_array[pos]>>1;
		unsigned long start1 = add1;
		A[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));
		start1+=add1;
		A[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));
		start1+=add1;
		A[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));
		start1+=add1;
		A[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));
		start1+=add1;
		A[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));
		start1+=add1;
		A[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));
		start1+=add1;
		A[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));
		start1+=add1;
		unsigned long add2 = my_array[pos+1]>>1;
		unsigned long start2 = add2;
		A[1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+8],*(Wp+start2));
		start2+=add2;
		A[2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+8],*(Wp+start2));
		start2+=add2;
		A[3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+8],*(Wp+start2));
		start2+=add2;
		A[4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+8],*(Wp+start2));
		start2+=add2;
		A[5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+8],*(Wp+start2));
		start2+=add2;
		A[6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+8],*(Wp+start2));
		start2+=add2;
		A[7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+8],*(Wp+start2));
		start2+=add2;
		u=AddModSpe(A[0+j],A[0+j+8]);
		t=SubModSpe(A[0+j],A[0+j+8]);
		A[0+j]=u;
		A[0+j+8]=t;
		u=AddModSpe(A[1+j],A[1+j+8]);
		t=SubModSpe(A[1+j],A[1+j+8]);
		A[1+j]=u;
		A[1+j+8]=t;
		u=AddModSpe(A[2+j],A[2+j+8]);
		t=SubModSpe(A[2+j],A[2+j+8]);
		A[2+j]=u;
		A[2+j+8]=t;
		u=AddModSpe(A[3+j],A[3+j+8]);
		t=SubModSpe(A[3+j],A[3+j+8]);
		A[3+j]=u;
		A[3+j+8]=t;
		u=AddModSpe(A[4+j],A[4+j+8]);
		t=SubModSpe(A[4+j],A[4+j+8]);
		A[4+j]=u;
		A[4+j+8]=t;
		u=AddModSpe(A[5+j],A[5+j+8]);
		t=SubModSpe(A[5+j],A[5+j+8]);
		A[5+j]=u;
		A[5+j+8]=t;
		u=AddModSpe(A[6+j],A[6+j+8]);
		t=SubModSpe(A[6+j],A[6+j+8]);
		A[6+j]=u;
		A[6+j+8]=t;
		u=AddModSpe(A[7+j],A[7+j+8]);
		t=SubModSpe(A[7+j],A[7+j+8]);
		A[7+j]=u;
		A[7+j+8]=t;
		for(int i=8;i<8;i+=8){
			A[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));
			start1+=add1;
			A[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));
			start1+=add1;
			A[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));
			start1+=add1;
			A[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));
			start1+=add1;
			A[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));
			start1+=add1;
			A[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));
			start1+=add1;
			A[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));
			start1+=add1;
			A[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));
			start1+=add1;
			A[i+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+8],*(Wp+start2));
			start2+=add2;
			A[i+1+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+8],*(Wp+start2));
			start2+=add2;
			A[i+2+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+8],*(Wp+start2));
			start2+=add2;
			A[i+3+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+8],*(Wp+start2));
			start2+=add2;
			A[i+4+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+8],*(Wp+start2));
			start2+=add2;
			A[i+5+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+8],*(Wp+start2));
			start2+=add2;
			A[i+6+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+8],*(Wp+start2));
			start2+=add2;
			A[i+7+j+8] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+8],*(Wp+start2));
			start2+=add2;
			u=AddModSpe(A[0+i+j],A[0+j+i+8]);
			t=SubModSpe(A[0+i+j],A[0+j+i+8]);
			A[0+j+i]=u;
			A[0+j+i+8]=t;
			u=AddModSpe(A[1+i+j],A[1+j+i+8]);
			t=SubModSpe(A[1+i+j],A[1+j+i+8]);
			A[1+j+i]=u;
			A[1+j+i+8]=t;
			u=AddModSpe(A[2+i+j],A[2+j+i+8]);
			t=SubModSpe(A[2+i+j],A[2+j+i+8]);
			A[2+j+i]=u;
			A[2+j+i+8]=t;
			u=AddModSpe(A[3+i+j],A[3+j+i+8]);
			t=SubModSpe(A[3+i+j],A[3+j+i+8]);
			A[3+j+i]=u;
			A[3+j+i+8]=t;
			u=AddModSpe(A[4+i+j],A[4+j+i+8]);
			t=SubModSpe(A[4+i+j],A[4+j+i+8]);
			A[4+j+i]=u;
			A[4+j+i+8]=t;
			u=AddModSpe(A[5+i+j],A[5+j+i+8]);
			t=SubModSpe(A[5+i+j],A[5+j+i+8]);
			A[5+j+i]=u;
			A[5+j+i+8]=t;
			u=AddModSpe(A[6+i+j],A[6+j+i+8]);
			t=SubModSpe(A[6+i+j],A[6+j+i+8]);
			A[6+j+i]=u;
			A[6+j+i+8]=t;
			u=AddModSpe(A[7+i+j],A[7+j+i+8]);
			t=SubModSpe(A[7+i+j],A[7+j+i+8]);
			A[7+j+i]=u;
			A[7+j+i+8]=t;
		}
	}
	unsigned long next = 8;
	small_butterflies2(A,next,0);
}
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


void DFT_rec_p1(int n,sfixn *A,sfixn *W,sfixn *B){
	if (n==1) {    
		return;
	}
	else if(n==FFT_THRESHOLD){
		PBPAS::ArrayBitReversalSpe(A);
		DFT_iter(A,W);
		return;
	}
	else{
		int next = n>>LOG_OF_POW;
		sfixn *W2 = W+n; 
		for (int k=0;k<n;k+=next){
			DFT_rec_p1(next, A+k, W2, B);
		}
		unsigned long array_of_pos[MY_POW] = {0,4,2,6,1,5,3,7};
		unsigned long pos=1;
		unsigned long add = array_of_pos[pos];
		unsigned long start=add;
		A[next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+1],W[start]);
		start+=add;
		A[next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+2],W[start]);
		start+=add;
		A[next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+3],W[start]);
		start+=add;
		A[next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+4],W[start]);
		start+=add;
		A[next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+5],W[start]);
		start+=add;
		A[next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+6],W[start]);
		start+=add;
		A[next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+7],W[start]);
		start+=add;

//----------------------------------------------------

		sfixn t,u;
		u = AddModSpe(A[0],A[next]);
		t = SubModSpe(A[0],A[next]);
		A[0] = u;
		A[next] = t;
		u = AddModSpe(A[1],A[1+next]);
		t = SubModSpe(A[1],A[1+next]);
		A[1] = u;
		A[1+next] = t;
		u = AddModSpe(A[2],A[2+next]);
		t = SubModSpe(A[2],A[2+next]);
		A[2] = u;
		A[2+next] = t;
		u = AddModSpe(A[3],A[3+next]);
		t = SubModSpe(A[3],A[3+next]);
		A[3] = u;
		A[3+next] = t;
		u = AddModSpe(A[4],A[4+next]);
		t = SubModSpe(A[4],A[4+next]);
		A[4] = u;
		A[4+next] = t;
		u = AddModSpe(A[5],A[5+next]);
		t = SubModSpe(A[5],A[5+next]);
		A[5] = u;
		A[5+next] = t;
		u = AddModSpe(A[6],A[6+next]);
		t = SubModSpe(A[6],A[6+next]);
		A[6] = u;
		A[6+next] = t;
		u = AddModSpe(A[7],A[7+next]);
		t = SubModSpe(A[7],A[7+next]);
		A[7] = u;
		A[7+next] = t;
		for(int j=8;j<next;j+=8){
			A[j+next+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+0],W[start]);
			start+=add;
			A[j+next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+1],W[start]);
			start+=add;
			A[j+next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+2],W[start]);
			start+=add;
			A[j+next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+3],W[start]);
			start+=add;
			A[j+next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+4],W[start]);
			start+=add;
			A[j+next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+5],W[start]);
			start+=add;
			A[j+next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+6],W[start]);
			start+=add;
			A[j+next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+7],W[start]);
			start+=add;

//----------------------------------------------------

			u = AddModSpe(A[j],A[j+next]);
			t = SubModSpe(A[j],A[j+next]);
			A[j] = u;
			A[j+next] = t;
			u = AddModSpe(A[j+1],A[j+1+next]);
			t = SubModSpe(A[j+1],A[j+1+next]);
			A[j+1] = u;
			A[j+1+next] = t;
			u = AddModSpe(A[j+2],A[j+2+next]);
			t = SubModSpe(A[j+2],A[j+2+next]);
			A[j+2] = u;
			A[j+2+next] = t;
			u = AddModSpe(A[j+3],A[j+3+next]);
			t = SubModSpe(A[j+3],A[j+3+next]);
			A[j+3] = u;
			A[j+3+next] = t;
			u = AddModSpe(A[j+4],A[j+4+next]);
			t = SubModSpe(A[j+4],A[j+4+next]);
			A[j+4] = u;
			A[j+4+next] = t;
			u = AddModSpe(A[j+5],A[j+5+next]);
			t = SubModSpe(A[j+5],A[j+5+next]);
			A[j+5] = u;
			A[j+5+next] = t;
			u = AddModSpe(A[j+6],A[j+6+next]);
			t = SubModSpe(A[j+6],A[j+6+next]);
			A[j+6] = u;
			A[j+6+next] = t;
			u = AddModSpe(A[j+7],A[j+7+next]);
			t = SubModSpe(A[j+7],A[j+7+next]);
			A[j+7] = u;
			A[j+7+next] = t;
		}
		pos = 2;
		for(int k=2*next; k<n;k+=2*next,pos+=2){
			unsigned long add1 = array_of_pos[pos];
			unsigned long start1=add1;
			A[k+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+1],W[start1]);
			start1+=add1;
			A[k+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+2],W[start1]);
			start1+=add1;
			A[k+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+3],W[start1]);
			start1+=add1;
			A[k+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+4],W[start1]);
			start1+=add1;
			A[k+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+5],W[start1]);
			start1+=add1;
			A[k+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+6],W[start1]);
			start1+=add1;
			A[k+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+7],W[start1]);
			start1+=add1;
			unsigned long add2 = array_of_pos[pos+1];
			unsigned long start2=add2;
			A[k+next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+1],W[start2]);
			start2+=add2;
			A[k+next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+2],W[start2]);
			start2+=add2;
			A[k+next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+3],W[start2]);
			start2+=add2;
			A[k+next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+4],W[start2]);
			start2+=add2;
			A[k+next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+5],W[start2]);
			start2+=add2;
			A[k+next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+6],W[start2]);
			start2+=add2;
			A[k+next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+7],W[start2]);
			start2+=add2;

//----------------------------------------------------

			u = AddModSpe(A[k],A[k+next]);
			t = SubModSpe(A[k],A[k+next]);
			A[k] = u;
			A[k+next] = t;
			u = AddModSpe(A[k+1],A[k+1+next]);
			t = SubModSpe(A[k+1],A[k+1+next]);
			A[k+1] = u;
			A[k+1+next] = t;
			u = AddModSpe(A[k+2],A[k+2+next]);
			t = SubModSpe(A[k+2],A[k+2+next]);
			A[k+2] = u;
			A[k+2+next] = t;
			u = AddModSpe(A[k+3],A[k+3+next]);
			t = SubModSpe(A[k+3],A[k+3+next]);
			A[k+3] = u;
			A[k+3+next] = t;
			u = AddModSpe(A[k+4],A[k+4+next]);
			t = SubModSpe(A[k+4],A[k+4+next]);
			A[k+4] = u;
			A[k+4+next] = t;
			u = AddModSpe(A[k+5],A[k+5+next]);
			t = SubModSpe(A[k+5],A[k+5+next]);
			A[k+5] = u;
			A[k+5+next] = t;
			u = AddModSpe(A[k+6],A[k+6+next]);
			t = SubModSpe(A[k+6],A[k+6+next]);
			A[k+6] = u;
			A[k+6+next] = t;
			u = AddModSpe(A[k+7],A[k+7+next]);
			t = SubModSpe(A[k+7],A[k+7+next]);
			A[k+7] = u;
			A[k+7+next] = t;
			for(int j=8;j<next;j+=8){
				A[k+j+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+0],W[start1]);
				start1+=add1;
				A[k+j+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+1],W[start1]);
				start1+=add1;
				A[k+j+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+2],W[start1]);
				start1+=add1;
				A[k+j+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+3],W[start1]);
				start1+=add1;
				A[k+j+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+4],W[start1]);
				start1+=add1;
				A[k+j+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+5],W[start1]);
				start1+=add1;
				A[k+j+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+6],W[start1]);
				start1+=add1;
				A[k+j+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+7],W[start1]);
				start1+=add1;
				A[k+j+next+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+0],W[start2]);
				start2+=add2;
				A[k+j+next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+1],W[start2]);
				start2+=add2;
				A[k+j+next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+2],W[start2]);
				start2+=add2;
				A[k+j+next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+3],W[start2]);
				start2+=add2;
				A[k+j+next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+4],W[start2]);
				start2+=add2;
				A[k+j+next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+5],W[start2]);
				start2+=add2;
				A[k+j+next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+6],W[start2]);
				start2+=add2;
				A[k+j+next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+7],W[start2]);
				start2+=add2;

//----------------------------------------------------

				u = AddModSpe(A[k+j],A[k+j+next]);
				t = SubModSpe(A[k+j],A[k+j+next]);
				A[k+j] = u;
				A[k+j+next] = t;
				u = AddModSpe(A[k+j+1],A[k+j+1+next]);
				t = SubModSpe(A[k+j+1],A[k+j+1+next]);
				A[k+j+1] = u;
				A[k+j+1+next] = t;
				u = AddModSpe(A[k+j+2],A[k+j+2+next]);
				t = SubModSpe(A[k+j+2],A[k+j+2+next]);
				A[k+j+2] = u;
				A[k+j+2+next] = t;
				u = AddModSpe(A[k+j+3],A[k+j+3+next]);
				t = SubModSpe(A[k+j+3],A[k+j+3+next]);
				A[k+j+3] = u;
				A[k+j+3+next] = t;
				u = AddModSpe(A[k+j+4],A[k+j+4+next]);
				t = SubModSpe(A[k+j+4],A[k+j+4+next]);
				A[k+j+4] = u;
				A[k+j+4+next] = t;
				u = AddModSpe(A[k+j+5],A[k+j+5+next]);
				t = SubModSpe(A[k+j+5],A[k+j+5+next]);
				A[k+j+5] = u;
				A[k+j+5+next] = t;
				u = AddModSpe(A[k+j+6],A[k+j+6+next]);
				t = SubModSpe(A[k+j+6],A[k+j+6+next]);
				A[k+j+6] = u;
				A[k+j+6+next] = t;
				u = AddModSpe(A[k+j+7],A[k+j+7+next]);
				t = SubModSpe(A[k+j+7],A[k+j+7+next]);
				A[k+j+7] = u;
				A[k+j+7+next] = t;
			}
		}
		small_butterflies0(A,next,0);
	}
}
void DFT_eff_p1(int n, int r,sfixn *A,sfixn *W,sfixn *B){
	#if MODDEBUG
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
		int r1 = r-1;
		sfixn pos = 0;
		sfixn m=n;
		int next;
		unsigned long adjust=0;
		sfixn twoFFTThreshold = FFT_THRESHOLD<<1;
		if(m & MY_MASK0){
			next=m>>3;
		}
		if(m & MY_MASK1){
			next=m>>1;
			adjust=2;
		}
		if(m & MY_MASK2){
			next=m>>2;
			adjust=1;
		}
		if(n>twoFFTThreshold){
			while (pos<n){
				while(m>FFT_THRESHOLD){
					Shuffle(m,A+pos,B);
					m>>=1;
				}
				pos += twoFFTThreshold;
				m = ((pos^(pos-1))+1)>>1;
			}
		}
		else{
			Shuffle(n,A,B);
		}
		sfixn *W2 = W+n; 
		for (int k=0;k<n;k+=next){
			DFT_rec_p1(next, A+k, W2, B);
		}
		unsigned long array_of_pos[MY_POW] = {0,4,2,6,1,5,3,7};
		pos=1;

//--------------------FIRST TWIDDLE FACTORS -------

		unsigned long add = array_of_pos[pos]>>adjust;
		unsigned long start=add;
		A[next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+1],W[start]);
		start+=add;
		A[next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+2],W[start]);
		start+=add;
		A[next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+3],W[start]);
		start+=add;
		A[next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+4],W[start]);
		start+=add;
		A[next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+5],W[start]);
		start+=add;
		A[next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+6],W[start]);
		start+=add;
		A[next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[next+7],W[start]);
		start+=add;

//----------------------------------------------------

		sfixn t,u;
		u = AddModSpe(A[0],A[next]);
		t = SubModSpe(A[0],A[next]);
		A[0] = u;
		A[next] = t;
		u = AddModSpe(A[1],A[1+next]);
		t = SubModSpe(A[1],A[1+next]);
		A[1] = u;
		A[1+next] = t;
		u = AddModSpe(A[2],A[2+next]);
		t = SubModSpe(A[2],A[2+next]);
		A[2] = u;
		A[2+next] = t;
		u = AddModSpe(A[3],A[3+next]);
		t = SubModSpe(A[3],A[3+next]);
		A[3] = u;
		A[3+next] = t;
		u = AddModSpe(A[4],A[4+next]);
		t = SubModSpe(A[4],A[4+next]);
		A[4] = u;
		A[4+next] = t;
		u = AddModSpe(A[5],A[5+next]);
		t = SubModSpe(A[5],A[5+next]);
		A[5] = u;
		A[5+next] = t;
		u = AddModSpe(A[6],A[6+next]);
		t = SubModSpe(A[6],A[6+next]);
		A[6] = u;
		A[6+next] = t;
		u = AddModSpe(A[7],A[7+next]);
		t = SubModSpe(A[7],A[7+next]);
		A[7] = u;
		A[7+next] = t;
		for(int j=8;j<next;j+=8){
			A[j+next+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+0],W[start]);
			start+=add;
			A[j+next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+1],W[start]);
			start+=add;
			A[j+next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+2],W[start]);
			start+=add;
			A[j+next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+3],W[start]);
			start+=add;
			A[j+next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+4],W[start]);
			start+=add;
			A[j+next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+5],W[start]);
			start+=add;
			A[j+next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+6],W[start]);
			start+=add;
			A[j+next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+next+7],W[start]);
			start+=add;

//----------------------------------------------------

			u = AddModSpe(A[j],A[j+next]);
			t = SubModSpe(A[j],A[j+next]);
			A[j] = u;
			A[j+next] = t;
			u = AddModSpe(A[j+1],A[j+1+next]);
			t = SubModSpe(A[j+1],A[j+1+next]);
			A[j+1] = u;
			A[j+1+next] = t;
			u = AddModSpe(A[j+2],A[j+2+next]);
			t = SubModSpe(A[j+2],A[j+2+next]);
			A[j+2] = u;
			A[j+2+next] = t;
			u = AddModSpe(A[j+3],A[j+3+next]);
			t = SubModSpe(A[j+3],A[j+3+next]);
			A[j+3] = u;
			A[j+3+next] = t;
			u = AddModSpe(A[j+4],A[j+4+next]);
			t = SubModSpe(A[j+4],A[j+4+next]);
			A[j+4] = u;
			A[j+4+next] = t;
			u = AddModSpe(A[j+5],A[j+5+next]);
			t = SubModSpe(A[j+5],A[j+5+next]);
			A[j+5] = u;
			A[j+5+next] = t;
			u = AddModSpe(A[j+6],A[j+6+next]);
			t = SubModSpe(A[j+6],A[j+6+next]);
			A[j+6] = u;
			A[j+6+next] = t;
			u = AddModSpe(A[j+7],A[j+7+next]);
			t = SubModSpe(A[j+7],A[j+7+next]);
			A[j+7] = u;
			A[j+7+next] = t;
		}
		pos = 2;
		for(int k=2*next; k<n;k+=2*next,pos+=2){
			unsigned long add1 = array_of_pos[pos]>>adjust;
			unsigned long start1=add1;
			A[k+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+1],W[start1]);
			start1+=add1;
			A[k+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+2],W[start1]);
			start1+=add1;
			A[k+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+3],W[start1]);
			start1+=add1;
			A[k+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+4],W[start1]);
			start1+=add1;
			A[k+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+5],W[start1]);
			start1+=add1;
			A[k+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+6],W[start1]);
			start1+=add1;
			A[k+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+7],W[start1]);
			start1+=add1;
			unsigned long add2 = array_of_pos[pos+1]>>adjust;
			unsigned long start2=add2;
			A[k+next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+1],W[start2]);
			start2+=add2;
			A[k+next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+2],W[start2]);
			start2+=add2;
			A[k+next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+3],W[start2]);
			start2+=add2;
			A[k+next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+4],W[start2]);
			start2+=add2;
			A[k+next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+5],W[start2]);
			start2+=add2;
			A[k+next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+6],W[start2]);
			start2+=add2;
			A[k+next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+next+7],W[start2]);
			start2+=add2;

//----------------------------------------------------

			u = AddModSpe(A[k],A[k+next]);
			t = SubModSpe(A[k],A[k+next]);
			A[k] = u;
			A[k+next] = t;
			u = AddModSpe(A[k+1],A[k+1+next]);
			t = SubModSpe(A[k+1],A[k+1+next]);
			A[k+1] = u;
			A[k+1+next] = t;
			u = AddModSpe(A[k+2],A[k+2+next]);
			t = SubModSpe(A[k+2],A[k+2+next]);
			A[k+2] = u;
			A[k+2+next] = t;
			u = AddModSpe(A[k+3],A[k+3+next]);
			t = SubModSpe(A[k+3],A[k+3+next]);
			A[k+3] = u;
			A[k+3+next] = t;
			u = AddModSpe(A[k+4],A[k+4+next]);
			t = SubModSpe(A[k+4],A[k+4+next]);
			A[k+4] = u;
			A[k+4+next] = t;
			u = AddModSpe(A[k+5],A[k+5+next]);
			t = SubModSpe(A[k+5],A[k+5+next]);
			A[k+5] = u;
			A[k+5+next] = t;
			u = AddModSpe(A[k+6],A[k+6+next]);
			t = SubModSpe(A[k+6],A[k+6+next]);
			A[k+6] = u;
			A[k+6+next] = t;
			u = AddModSpe(A[k+7],A[k+7+next]);
			t = SubModSpe(A[k+7],A[k+7+next]);
			A[k+7] = u;
			A[k+7+next] = t;
			for(int j=8;j<next;j+=8){
				A[k+j+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+0],W[start1]);
				start1+=add1;
				A[k+j+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+1],W[start1]);
				start1+=add1;
				A[k+j+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+2],W[start1]);
				start1+=add1;
				A[k+j+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+3],W[start1]);
				start1+=add1;
				A[k+j+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+4],W[start1]);
				start1+=add1;
				A[k+j+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+5],W[start1]);
				start1+=add1;
				A[k+j+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+6],W[start1]);
				start1+=add1;
				A[k+j+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+7],W[start1]);
				start1+=add1;
				A[k+j+next+0]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+0],W[start2]);
				start2+=add2;
				A[k+j+next+1]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+1],W[start2]);
				start2+=add2;
				A[k+j+next+2]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+2],W[start2]);
				start2+=add2;
				A[k+j+next+3]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+3],W[start2]);
				start2+=add2;
				A[k+j+next+4]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+4],W[start2]);
				start2+=add2;
				A[k+j+next+5]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+5],W[start2]);
				start2+=add2;
				A[k+j+next+6]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+6],W[start2]);
				start2+=add2;
				A[k+j+next+7]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+next+7],W[start2]);
				start2+=add2;

//----------------------------------------------------

				u = AddModSpe(A[k+j],A[k+j+next]);
				t = SubModSpe(A[k+j],A[k+j+next]);
				A[k+j] = u;
				A[k+j+next] = t;
				u = AddModSpe(A[k+j+1],A[k+j+1+next]);
				t = SubModSpe(A[k+j+1],A[k+j+1+next]);
				A[k+j+1] = u;
				A[k+j+1+next] = t;
				u = AddModSpe(A[k+j+2],A[k+j+2+next]);
				t = SubModSpe(A[k+j+2],A[k+j+2+next]);
				A[k+j+2] = u;
				A[k+j+2+next] = t;
				u = AddModSpe(A[k+j+3],A[k+j+3+next]);
				t = SubModSpe(A[k+j+3],A[k+j+3+next]);
				A[k+j+3] = u;
				A[k+j+3+next] = t;
				u = AddModSpe(A[k+j+4],A[k+j+4+next]);
				t = SubModSpe(A[k+j+4],A[k+j+4+next]);
				A[k+j+4] = u;
				A[k+j+4+next] = t;
				u = AddModSpe(A[k+j+5],A[k+j+5+next]);
				t = SubModSpe(A[k+j+5],A[k+j+5+next]);
				A[k+j+5] = u;
				A[k+j+5+next] = t;
				u = AddModSpe(A[k+j+6],A[k+j+6+next]);
				t = SubModSpe(A[k+j+6],A[k+j+6+next]);
				A[k+j+6] = u;
				A[k+j+6+next] = t;
				u = AddModSpe(A[k+j+7],A[k+j+7+next]);
				t = SubModSpe(A[k+j+7],A[k+j+7+next]);
				A[k+j+7] = u;
				A[k+j+7+next] = t;
			}
		}
		if(n& MY_MASK1){
		small_butterflies1(A,next,0);
		}
		if(n& MY_MASK2){
		small_butterflies2(A,next,0);
		}
		if(n& MY_MASK0){
		small_butterflies0(A,next,0);
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
			case 32:
				PBPAS::ArrayBitReversalSpe32(A);
				DFT_iter32(A,W);
				break;
			case 64:
				PBPAS::ArrayBitReversalSpe64(A);
				DFT_iter64(A,W);
				break;
			case 128:
				PBPAS::ArrayBitReversalSpe128(A);
				DFT_iter128(A,W);
				break;
			case 256:
				PBPAS::ArrayBitReversalSpe256(A);
				DFT_iter256(A,W);
				break;
			case 512:
				PBPAS::ArrayBitReversalSpe512(A);
				DFT_iter512(A,W);
				break;
			default:
				return;
		}
	}
	#if MODDEBUG
	if(test!=A[1]){
		std::cout<<n << " "<<test<< " "<<A[1]<<std::endl;
		exit(1);
	}
	#endif
}

void InvDFTKeepMont_eff_p1(int n, int r, sfixn *A,sfixn *W,  sfixn *B,sfixn invn){
	DFT_eff_p1(n,r,A,W,B);
	for(int i=0; i<n; i++) {
		A[i]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[i], invn);
	}
}
void InvDFT_eff_p1(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn){
	DFT_eff_p1(n,r,A,W,B);
	for(int i=0; i<n; i++) {
		A[i]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[i], invn);
	}
}
}
