namespace FURERPBPAS{
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
//		: "a"(a),"d"(b),"r"((unsigned long) INV_PRIME),"r"((unsigned long) COMPLETE),"r"((unsigned long) MY_PRIME)
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
MONTGOMERYMUL
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
//		: "a"(a),"rm"(b),"b"((sfixn) INV_PRIME),"c"((sfixn) MY_PRIME)
//		:"rsi","rdi");
//	return a;
//}
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
	rootR=fastExp(MontMulModSpe_OPT3_AS_GENE_INLINE(GENERATOR,R_2<<BRsft),(MY_PRIME-1)>>e);
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

GENERATE_SMALL_FFT
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


void DFT_rec(int n,sfixn *A,sfixn *W,sfixn *B){
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
			DFT_rec(next, A+k, W2, B);
		}
		unsigned long array_of_pos[MY_POW] = ARRAY_FOR_MULS;
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
		SMALL_FFTR
	}
}
void DFT_eff(int n, int r,sfixn *A,sfixn *W,sfixn *B){
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
		FIRST_LEVEL_SHUFFLE
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
			DFT_rec(next, A+k, W2, B);
		}
		unsigned long array_of_pos[MY_POW] = ARRAY_FOR_MULS;
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
		SMALL_FFTs
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
	#if MODDEBUG
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
