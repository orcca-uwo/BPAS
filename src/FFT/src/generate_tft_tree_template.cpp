//using namespace PBPAS1;

namespace TFT_tree{

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

static inline void  TFT_AddSubSpeSSEModInplace(sfixn* a0,sfixn* a1, sfixn* a2, sfixn* a3){
  sfixn u1 = AddModSpe(*a0,*a2);
  sfixn u2 = SubModSpe(*a0,*a2);
  sfixn u3 = AddModSpe(*a1,*a3);
  sfixn u4 = SubModSpe(*a1,*a3);
  *a0 = u1;
  *a1 = u3;
  *a2 = u2;
  *a3 = u4;
  return;
}

void inline TFT_2POINT(sfixn *A,sfixn *W){
	sfixn u = A[0];
	sfixn t = A[1];
	A[0] = AddModSpe(u,t);
	A[1] = SubModSpe(u,t);
}

void inline TFT_4POINT(sfixn *A,sfixn *W){
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

	TFT_AddSubSpeSSEModInplace(A,A+1,A+2,A+3);
	
	sfixn tmp =A[2];
	A[2]=A[1];
	A[1]=tmp;


}

void inline TFT_8POINT(sfixn *A,sfixn *W){
	
	sfixn *Wp = W + (8<<1)-4;
	sfixn u = A[0];
	sfixn t = A[4];
	A[0] = AddModSpe(u,t);
	A[4] = SubModSpe(u,t);
	u = A[2];
	t = A[6];
	A[2] = AddModSpe(u,t);
	A[6] = SubModSpe(u,t);
	u = A[1];
	t = A[5];
	A[1] = AddModSpe(u,t);
	A[5] = SubModSpe(u,t);
	u = A[3];
	t = A[7];
	A[3] = AddModSpe(u,t);
	A[7] = SubModSpe(u,t);


	A[6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6],*(Wp-3));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-3));

	TFT_AddSubSpeSSEModInplace(A,A+4,A+2,A+6);
	TFT_AddSubSpeSSEModInplace(A+1,A+5,A+3,A+7);

	A[5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5],*(Wp-11));
	A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3],*(Wp-10));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-9));


	TFT_AddSubSpeSSEModInplace(A,A+4,A+1,A+5);
	TFT_AddSubSpeSSEModInplace(A+2,A+6,A+3,A+7);


}


void inline TFT_16POINT(sfixn *A,sfixn *W){
	sfixn *Wp = W + (16<<1)-4;
	sfixn u = A[0];
	sfixn t = A[8];
	A[0] = AddModSpe(u,t);
	A[8] = SubModSpe(u,t);
	u = A[4];
	t = A[12];
	A[4] = AddModSpe(u,t);
	A[12] = SubModSpe(u,t);
	u = A[2];
	t = A[10];
	A[2] = AddModSpe(u,t);
	A[10] = SubModSpe(u,t);
	u = A[6];
	t = A[14];
	A[6] = AddModSpe(u,t);
	A[14] = SubModSpe(u,t);
	u = A[1];
	t = A[9];
	A[1] = AddModSpe(u,t);
	A[9] = SubModSpe(u,t);
	u = A[5];
	t = A[13];
	A[5] = AddModSpe(u,t);
	A[13] = SubModSpe(u,t);
	u = A[3];
	t = A[11];
	A[3] = AddModSpe(u,t);
	A[11] = SubModSpe(u,t);
	u = A[7];
	t = A[15];
	A[7] = AddModSpe(u,t);
	A[15] = SubModSpe(u,t);
	A[12] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[12],*(Wp-3));
	A[14] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[14],*(Wp-3));
	A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[13],*(Wp-3));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[15],*(Wp-3));


	TFT_AddSubSpeSSEModInplace(A,A+8,A+4,A+12);
	TFT_AddSubSpeSSEModInplace(A+2,A+10,A+6,A+14);
	TFT_AddSubSpeSSEModInplace(A+1,A+9,A+5,A+13);
	TFT_AddSubSpeSSEModInplace(A+3,A+11,A+7,A+15);


	A[10] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[10],*(Wp-11));
	A[6] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[6],*(Wp-10));
	A[14] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[14],*(Wp-9));
	A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[11],*(Wp-11));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[7],*(Wp-10));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[15],*(Wp-9));

	TFT_AddSubSpeSSEModInplace(A,A+8,A+2,A+10);
	TFT_AddSubSpeSSEModInplace(A+4,A+12,A+6,A+14);
	TFT_AddSubSpeSSEModInplace(A+1,A+9,A+3,A+11);
	TFT_AddSubSpeSSEModInplace(A+5,A+13,A+7,A+15);

	A[9] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[9],*(Wp-27));
	A[5] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[5],*(Wp-26));
	A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[13],*(Wp-25));
	A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[3],*(Wp-24));
	A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[11],*(Wp-23));
	A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[7],*(Wp-22));
	A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE2(A[15],*(Wp-21));

	TFT_AddSubSpeSSEModInplace(A,A+8,A+1,A+9);
	TFT_AddSubSpeSSEModInplace(A+4,A+12,A+5,A+13);
	TFT_AddSubSpeSSEModInplace(A+2,A+10,A+3,A+11);
	TFT_AddSubSpeSSEModInplace(A+6,A+14,A+7,A+15);
}
void inline TFT_iter32(sfixn *A,sfixn *W){

	sfixn *Wp = W + (32<<1)-4;
	sfixn* Wt;

	
		sfixn u = A[ 0];
		sfixn t = A[ 16];
		A[ 0] = AddModSpe(u,t);
		A[ 16] = SubModSpe(u,t);
		u = A[ 8];
		t = A[ 24];
		A[ 8] = AddModSpe(u,t);
		A[ 24] = SubModSpe(u,t);
		u = A[ 4];
		t = A[ 20];
		A[ 4] = AddModSpe(u,t);
		A[ 20] = SubModSpe(u,t);
		u = A[ 12];
		t = A[ 28];
		A[ 12] = AddModSpe(u,t);
		A[ 28] = SubModSpe(u,t);
		u = A[ 2];
		t = A[ 18];
		A[ 2] = AddModSpe(u,t);
		A[ 18] = SubModSpe(u,t);
		u = A[ 10];
		t = A[ 26];
		A[ 10] = AddModSpe(u,t);
		A[ 26] = SubModSpe(u,t);
		u = A[ 6];
		t = A[ 22];
		A[ 6] = AddModSpe(u,t);
		A[ 22] = SubModSpe(u,t);
		u = A[ 14];
		t = A[ 30];
		A[ 14] = AddModSpe(u,t);
		A[ 30] = SubModSpe(u,t);
		A[ 24] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 24],*(Wp-3));
		A[ 28] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 28],*(Wp-3));
		A[ 26] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 26],*(Wp-3));
		A[ 30] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 30],*(Wp-3));


		TFT_AddSubSpeSSEModInplace(A,A+16,A+8,A+24);
		TFT_AddSubSpeSSEModInplace(A+4,A+20,A+12,A+28);
		TFT_AddSubSpeSSEModInplace(A+2,A+18,A+10,A+26);
		TFT_AddSubSpeSSEModInplace(A+6,A+22,A+14,A+30);


		A[ 20] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 20],*(Wp-11));
		A[ 12] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 12],*(Wp-10));
		A[ 28] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 28],*(Wp-9));
		A[ 22] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 22],*(Wp-11));
		A[ 14] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 14],*(Wp-10));
		A[ 30] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 30],*(Wp-9));

	
		TFT_AddSubSpeSSEModInplace(A,A+16,A+4,A+20);
		TFT_AddSubSpeSSEModInplace(A+8,A+24,A+12,A+28);
		TFT_AddSubSpeSSEModInplace(A+2,A+18,A+6,A+22);
		TFT_AddSubSpeSSEModInplace(A+10,A+26,A+14,A+30);


		A[ 18] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 18],*(Wp-27));
		A[ 10] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 10],*(Wp-26));
		A[ 26] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 26],*(Wp-25));
		A[ 6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 6],*(Wp-24));
		A[ 22] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 22],*(Wp-23));
		A[ 14] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 14],*(Wp-22));
		A[ 30] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[ 30],*(Wp-21));
	

		TFT_AddSubSpeSSEModInplace(A,A+16,A+2,A+18);
		TFT_AddSubSpeSSEModInplace(A+8,A+24,A+10,A+26);
		TFT_AddSubSpeSSEModInplace(A+4,A+20,A+6,A+22);
		TFT_AddSubSpeSSEModInplace(A+12,A+28,A+14,A+30);

	
		//>>>>>>>>k+16>>>>>
	
		u = A[1];
		t = A[17];
		A[1] = AddModSpe(u,t);
		A[17] = SubModSpe(u,t);
		u = A[9];
		t = A[25];
		A[9] = AddModSpe(u,t);
		A[25] = SubModSpe(u,t);
		u = A[5];
		t = A[21];
		A[5] = AddModSpe(u,t);
		A[21] = SubModSpe(u,t);
		u = A[13];
		t = A[29];
		A[13] = AddModSpe(u,t);
		A[29] = SubModSpe(u,t);
		u = A[3];
		t = A[19];
		A[3] = AddModSpe(u,t);
		A[19] = SubModSpe(u,t);
		u = A[11];
		t = A[27];
		A[11] = AddModSpe(u,t);
		A[27] = SubModSpe(u,t);
		u = A[7];
		t = A[23];
		A[7] = AddModSpe(u,t);
		A[23] = SubModSpe(u,t);
		u = A[15];
		t = A[31];
		A[15] = AddModSpe(u,t);
		A[31] = SubModSpe(u,t);
		A[25] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[25],*(Wp-3));
		A[29] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[29],*(Wp-3));
		A[27] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[27],*(Wp-3));
		A[31] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[31],*(Wp-3));


		TFT_AddSubSpeSSEModInplace(A+1,A+17,A+9,A+25);
		TFT_AddSubSpeSSEModInplace(A+5,A+21,A+13,A+29);
		TFT_AddSubSpeSSEModInplace(A+3,A+19,A+11,A+27);
		TFT_AddSubSpeSSEModInplace(A+7,A+23,A+15,A+31);



		A[21] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[21],*(Wp-11));
		A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[13],*(Wp-10));
		A[29] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[29],*(Wp-9));
		A[23] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[23],*(Wp-11));
		A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[15],*(Wp-10));
		A[31] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[31],*(Wp-9));

	

		TFT_AddSubSpeSSEModInplace(A+1,A+17,A+5,A+21);
		TFT_AddSubSpeSSEModInplace(A+9,A+25,A+13,A+29);
		TFT_AddSubSpeSSEModInplace(A+3,A+19,A+7,A+23);
		TFT_AddSubSpeSSEModInplace(A+11,A+27,A+15,A+31);
	


		A[19] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[19],*(Wp-27));
		A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[11],*(Wp-26));
		A[27] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[27],*(Wp-25));
		A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp-24));
		A[23] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[23],*(Wp-23));
		A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[15],*(Wp-22));
		A[31] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[31],*(Wp-21));


		TFT_AddSubSpeSSEModInplace(A+1,A+17,A+3,A+19);
		TFT_AddSubSpeSSEModInplace(A+9,A+25,A+11,A+27);
		TFT_AddSubSpeSSEModInplace(A+5,A+21,A+7,A+23);
		TFT_AddSubSpeSSEModInplace(A+13,A+29,A+15,A+31);
	
	Wp = Wp - 28;
	Wp = Wp-(1L<<5);
	Wt = Wp;
	Wp = Wt;
	A[17] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[17],*(Wp+1));

	TFT_AddSubSpeSSEModInplace(A+0,A+16,A+1,A+17);

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	Wp = Wt;
		A[9] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[9],*(Wp+2+0));
		A[25] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[25],*(Wp+2+1));
	
		TFT_AddSubSpeSSEModInplace(A+8,A+24,A+9,A+25);
	Wp = Wt;
		A[5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5],*(Wp+4+0));
		A[21] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[21],*(Wp+4+1));
	
		TFT_AddSubSpeSSEModInplace(A+4,A+20,A+5,A+21);
	Wp = Wt;
		A[13] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[13],*(Wp+6+0));
		A[29] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[29],*(Wp+6+1));

		TFT_AddSubSpeSSEModInplace(A+12,A+28,A+13,A+29);

	Wp = Wt;
		A[3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3],*(Wp+8+0));
		A[19] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[19],*(Wp+8+1));
	
		TFT_AddSubSpeSSEModInplace(A+2,A+18,A+3,A+19);
	Wp = Wt;
		A[11] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[11],*(Wp+10+0));
		A[27] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[27],*(Wp+10+1));

		TFT_AddSubSpeSSEModInplace(A+10,A+26,A+11,A+27);
	Wp = Wt;

		A[7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7],*(Wp+12+0));
		A[23] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[23],*(Wp+12+1));

		TFT_AddSubSpeSSEModInplace(A+6,A+22,A+7,A+23);
	Wp = Wt;
		A[15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[15],*(Wp+14+0));
		A[31] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[31],*(Wp+14+1));
	
		TFT_AddSubSpeSSEModInplace(A+14,A+30,A+15,A+31);
}



  void TFT_twiddle(sfixn *A, sfixn *KRT,int s, int r, int n, int oldn){

    int a=oldn,b=0;
    if(n<oldn){
      
      for(int i=0;i<(log2 (oldn))-(log2 (n));i++){
	b=b+a;
	a=a/2;
      }
      
    }	 
    
   for(int i=0;i<s;i++){
     for(int j=0;j<r;j++){
       if( i*j > 0){ 
	 // std::cout << "i*r+j  is "<< i*r+j << std::endl;
	 A[i*r+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i*r+j], KRT[i*j+b]);  
	 // std::cout << "i*j+b is "<<i*j+b << std::endl;
       }
     }
   }
     
   //permutation//
   sfixn *invectmp = (sfixn *)calloc(s, sizeof(sfixn)); 
   for(int i=0;i<s;i++){
     invectmp[i]=0;
   }
   sfixn *B = (sfixn *)calloc(s, sizeof(sfixn)); 
    for(int i=0;i<s;i++){
     B[i]=0;
   }
    // std::cout << "arrive here  " << std::endl;   
   for(int i=0;i<r;i++){
     for(int j=0;j<s;j++){

       invectmp[j] = A[j*r+i];

     }

     if(s==4){
       //reverse to tft order
       sfixn temp = invectmp[1];
       invectmp[1]=invectmp[2];
       invectmp[2]=temp;	
     }
     else if(s==8){
       //std::cout << " bit reverse at size 8   " << std::endl;
       PBPAS::ArrayBitReversalSpe8(invectmp);
 
     }
     else if(s==16){
       PBPAS::ArrayBitReversalSpe16(invectmp);
     }
     else{
       //permutate s>16
       // not sure
       Shuffle(s,invectmp,B);
     }


     for(int j=0;j<s;j++){

       //if(j*r+i<n){//only calculate relaxed size
	 A[j*r+i] = invectmp[j];
	 //}
     }
   }

   free(invectmp);
   free(B);
    
  }

  void Shuffle_tft(int s, sfixn *testvec, sfixn *SB1){

	//permutation for tft order
        if(s==2){
	  return;
        }
        else if(s==4){
	  sfixn temp = testvec[1];
	  testvec[1]=testvec[2];
	  testvec[2]=temp;
	}
	else if(s==8){
	  PBPAS::ArrayBitReversalSpe8(testvec);
	}
	else if(s==16){
	  PBPAS::ArrayBitReversalSpe16(testvec);
	}
	else if(s==32){
	  PBPAS::ArrayBitReversalSpe32(testvec);
	}
	else if(s==64){
	  PBPAS::ArrayBitReversalSpe64(testvec);
	}
	else if(s==128){
	  PBPAS::ArrayBitReversalSpe128(testvec);
	}
	else if(s==256){
	  PBPAS::ArrayBitReversalSpe256(testvec);
	}
	else if(s==512){
	  PBPAS::ArrayBitReversalSpe512(testvec);
	}
	else if(s==1024){
	  PBPAS::ArrayBitReversalSpe(testvec);
	}
	else
	  {
	    Shuffle(s,testvec,SB1);
	    
	    for(int m=s>>1; m>2047; m=m>>1){
	      //std::cout <<"Shuffle "<< m<< std::endl;
		for(int i=0;i<s/m;i++){			
		  Shuffle(m,testvec+i*m,SB1);
		  //std::cout <<"Shuffle "<< m << " ,testvec+  "<< i<<" * "<<m  << std::endl;
		}
	    }

	    for(int i=0;i< (s>>10);i++){
		PBPAS::ArrayBitReversalSpe(testvec+ (i<<10));	
	    }

	  }

  }
  void TFT_Basecase(int n, int r,sfixn *A,sfixn *W,sfixn *B){

	switch (n){
			case 2:
				TFT_2POINT(A,W);
				break;
			case 4:
				TFT_4POINT(A,W);
				break;
			case 8:
			        //PBPAS::ArrayBitReversalSpe8(A);
				TFT_8POINT(A,W);
				
				break;
		        case 16:
			        //PBPAS::ArrayBitReversalSpe16(A);
				TFT_16POINT(A,W);
			
				break;
		        case 32:

			        TFT_iter32(A,W);
			
				break;
               
	}


  }



  /* l is the length of the input vector
   invec is the input vector
   invectmp is the workspace vector
   m is the length of the input vector
   oldn is the smallest power of 2 larger than the size
         of the original input vector
   n  is the smallest power of 2 larger than the size
         of the current input vector, that is, invec
   basecase is the base case size from which a SLP code is used
   relax_size is the size after zero-padding (in the sense of the
       relaxed Cooley-Tukey)
   e is the log_2 of oldn
   SB1 is another work space array
   See void test9(int K)
   in Main/tests/FFT/src/fft.cpp
  */

  void TFT_Core(int oldn, int n, int l, int m, int basecase, int relax_size, sfixn *KRT, int e, sfixn *SB1,sfixn *invec, sfixn *invectmp ){




  int s,r,l_s,m_s,l_r,m_r;//for break down TFT
  int new_r_n, new_r_l, new_r_m, new_r_a, new_r_b;//for right part
  int new_l_n, new_l_l, new_l_m, new_l_a, new_l_b, new_l_a2;//for left part
  int bigern;
  bigern = pow(2, ceil(log2 (n)));

 
  //std::memcpy(invec, A, sizeof(A)); //can be faster
     
  if(n<=basecase){//N=17,basecase=32

    if(n<bigern){
      /*
      for(int i=0;i<n;i++){
	//cilk_for(int i=0;i<n;i++){
	invectmp[i]=invec[i];
      }
      
      
      //std::memcpy(invectmp, invec, sizeof(invectmp)/oldn*relax_size); //can be faster

      
      //for(int i=n;i<bigern;i++){
	//invectmp[i]=0;
	//}
      
      
      DFT_eff(bigern, e, invectmp, KRT, SB1);
      
      Shuffle_tft(bigern,invectmp,SB1);
      
      for(int i=0;i<n;i++){
	//cilk_for(int i=0;i<n;i++){
	invec[i]=invectmp[i];
      }
      */
    

      
 
      if(bigern>32){ 


	DFT_eff(bigern, e, invec, KRT, SB1);
	
	Shuffle_tft(bigern,invec,SB1);
      }
      else{

	TFT_Basecase(bigern, e, invec, KRT, SB1);

      }





      
    }
    else{
  
      if(bigern>32){ 

	  DFT_eff(n, e, invec, KRT, SB1);
	
	  Shuffle_tft(n,invec,SB1);
      }
      else{

	TFT_Basecase(n, e, invec, KRT, SB1);

      }


    }
  
    return;
  }


  if (bigern > pow(basecase,2)){
        s =  basecase;
        r =  bigern/basecase;
  }
  else{

    s = pow(2, ceil((log2 (bigern))/2.0));
    r = pow(2, floor( (log2 (bigern))/2.0));
  }
 
  if (s < 1){
    return;
  }
 
  l_s = ceil(l/(r*1.0));
  m_s = ceil(m/(r*1.0));
  l_r = std::min(l,r);
  m_r = std::min(m,r);

  //int total_step= (log2 (oldn));
  int bigern_step=(log2 (bigern));
 ///////right part and twddile part////////////////////////////////////////

    new_r_n = s;
    new_r_l = l_s;
    new_r_m = m_s;
 
    new_r_a = relax_size/l_s/l_r;
    new_r_b = l_r;
   

   

    /*
    std::cout<<" relax size is "<<relax_size<<std::endl;

    std::cout<<" s  is "<<s<<" r is  "<<r<<"  l_s is  "<<l_s<<" m_s  is  "<<m_s<<" l_r  is "<<l_r<<"  m_r is  "<<m_r<<std::endl;

    std::cout<<" new_r_n  is "<<new_r_n<<",  new_r_l  is "<<new_r_l<<",  new_r_m is "<<new_r_m<<",  new_r_a  is "<<new_r_a<<",  new_r_b is "<<new_r_b<<std::endl;
    */

    if( new_r_n <= basecase ){
      //sfixn *invectmp = (sfixn *)calloc(new_r_n, sizeof(sfixn));

      	int right_step= log2 (new_r_n);
	
	int ff=oldn,f=0;

	if(new_r_n<oldn){

	  for(int i=0;i< (e - right_step);i++){
	    f=f+ff;
	    ff=ff>>1;
	  }

	}
	int t_a=oldn,t_b=0;//for twiddle matrix omega

	if(bigern<oldn){      
	  for(int i=0;i<( e - bigern_step );i++){
	    t_b=t_b+t_a;
	    t_a=t_a>>1;
	  }
	}

	/*
	for(int j=0;j<new_r_b;j++){
	  
	  int pi=j*new_r_n;
	  
	  for(int i=0;i<new_r_l;i++){ //i<new_r_l only calculate relaxed size
	    
	    invectmp[i+pi]=invec[i*new_r_b+j];
	  }

	  DFT_eff(new_r_n, right_step, invectmp+pi, KRT+f, SB1);//KRT+n is for level 2 omega

	  for(int i=0;i<new_r_n;i++){//i<new_r_n only calculate relaxed size
	    invectmp[i+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[i+pi] , KRT[i*j+t_b]);
	  }
	  
	  if(s>=4){
	    Shuffle_tft(s,invectmp+pi,SB1);
	  }
	  
	  for(int i=0;i<new_r_l;i++){
	    invec[i*new_r_b +j]=invectmp[i+pi];  	 
	  }
	  }*/
	
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
	//transpose_serial(invec,new_r_b,invectmp,s,0,s,0,new_r_b);
	transpose(invec,new_r_b,invectmp,s,0,s,0,new_r_b);


	# pragma cilk_grainsize = TFT_grainsize;
	cilk_for(int j=0;j<new_r_b;j++){
	  int pi=j*new_r_n;

	  
	  //DFT_eff(new_r_n, right_step, invectmp+pi, KRT+f, SB1);//KRT+n is for level 2 omega


	  if(new_r_n>32){ 
	  
	    DFT_eff(new_r_n, right_step, invectmp+pi, KRT+f, SB1);//KRT+n is for level 2 omega

	  }
	  else{
	    TFT_Basecase(new_r_n, right_step, invectmp+pi, KRT+f, SB1);//KRT+n is for level 2 omega
	  }



	  if(new_r_n == 4){

	    invectmp[0+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[0+pi] , KRT[0*j+t_b]);
	    invectmp[1+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[1+pi] , KRT[2*j+t_b]);
	    invectmp[2+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[2+pi] , KRT[1*j+t_b]);
	    invectmp[3+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[3+pi] , KRT[3*j+t_b]);
	  }	  
	  else if(new_r_n == 8){
	    invectmp[0+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[0+pi] , KRT[0*j+t_b]);
	    invectmp[1+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[1+pi] , KRT[4*j+t_b]);
	    invectmp[2+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[2+pi] , KRT[2*j+t_b]);
	    invectmp[3+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[3+pi] , KRT[6*j+t_b]);
	    invectmp[4+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[4+pi] , KRT[1*j+t_b]);
	    invectmp[5+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[5+pi] , KRT[5*j+t_b]);
	    invectmp[6+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[6+pi] , KRT[3*j+t_b]);
	    invectmp[7+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[7+pi] , KRT[7*j+t_b]);

	  }
	  else if(new_r_n == 16){
	    invectmp[0+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[0+pi] , KRT[0*j+t_b]);
	    invectmp[1+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[1+pi] , KRT[8*j+t_b]);
	    invectmp[2+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[2+pi] , KRT[4*j+t_b]);
	    invectmp[3+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[3+pi] , KRT[12*j+t_b]);
	    invectmp[4+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[4+pi] , KRT[2*j+t_b]);
	    invectmp[5+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[5+pi] , KRT[10*j+t_b]);
	    invectmp[6+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[6+pi] , KRT[6*j+t_b]);
	    invectmp[7+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[7+pi] , KRT[14*j+t_b]);
	    invectmp[8+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[8+pi] , KRT[1*j+t_b]);
	    invectmp[9+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[9+pi] , KRT[9*j+t_b]);
	    invectmp[10+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[10+pi] , KRT[5*j+t_b]);
	    invectmp[11+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[11+pi] , KRT[13*j+t_b]);
	    invectmp[12+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[12+pi] , KRT[3*j+t_b]);
	    invectmp[13+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[13+pi] , KRT[11*j+t_b]);
	    invectmp[14+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[14+pi] , KRT[7*j+t_b]);
	    invectmp[15+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[15+pi] , KRT[15*j+t_b]);


	  }
	  else if(new_r_n == 32){
	    invectmp[0+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[0+pi] , KRT[0*j+t_b]);
	    invectmp[1+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[1+pi] , KRT[16*j+t_b]);
	    invectmp[2+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[2+pi] , KRT[8*j+t_b]);
	    invectmp[3+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[3+pi] , KRT[24*j+t_b]);
	    invectmp[4+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[4+pi] , KRT[4*j+t_b]);
	    invectmp[5+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[5+pi] , KRT[20*j+t_b]);
	    invectmp[6+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[6+pi] , KRT[12*j+t_b]);
	    invectmp[7+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[7+pi] , KRT[28*j+t_b]);
	    invectmp[8+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[8+pi] , KRT[2*j+t_b]);
	    invectmp[9+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[9+pi] , KRT[18*j+t_b]);
	    invectmp[10+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[10+pi] , KRT[10*j+t_b]);
	    invectmp[11+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[11+pi] , KRT[26*j+t_b]);
	    invectmp[12+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[12+pi] , KRT[6*j+t_b]);
	    invectmp[13+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[13+pi] , KRT[22*j+t_b]);
	    invectmp[14+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[14+pi] , KRT[14*j+t_b]);
	    invectmp[15+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[15+pi] , KRT[30*j+t_b]);
	    invectmp[16+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[16+pi] , KRT[1*j+t_b]);
	    invectmp[17+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[17+pi] , KRT[17*j+t_b]);
	    invectmp[18+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[18+pi] , KRT[9*j+t_b]);
	    invectmp[19+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[19+pi] , KRT[25*j+t_b]);
	    invectmp[20+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[20+pi] , KRT[5*j+t_b]);
	    invectmp[21+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[21+pi] , KRT[21*j+t_b]);
	    invectmp[22+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[22+pi] , KRT[13*j+t_b]);
	    invectmp[23+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[23+pi] , KRT[29*j+t_b]);
	    invectmp[24+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[24+pi] , KRT[3*j+t_b]);
	    invectmp[25+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[25+pi] , KRT[19*j+t_b]);
	    invectmp[26+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[26+pi] , KRT[11*j+t_b]);
	    invectmp[27+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[27+pi] , KRT[27*j+t_b]);
	    invectmp[28+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[28+pi] , KRT[7*j+t_b]);
	    invectmp[29+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[29+pi] , KRT[23*j+t_b]);
	    invectmp[30+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[30+pi] , KRT[15*j+t_b]);
	    invectmp[31+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[31+pi] , KRT[31*j+t_b]);

	  }


	  else{
	    for(int i=0;i<new_r_n;i++){//i<new_r_n only calculate relaxed size
	      invectmp[i+pi]= MontMulModSpe_OPT3_AS_GENE_INLINE(invectmp[i+pi] , KRT[i*j+t_b]);
	      //std::cout << " i is   "<<i<<" i*j is "<<i*j << std::endl;
	      
	    }  
	  }

	
	  //std::cout << " s is   "<< s << std::endl;
	  if(s>32){
	    Shuffle_tft(s,invectmp+pi,SB1);
	  }

	}
	
	//transpose_serial(invectmp,s,invec,new_r_b,0,new_r_b,0,s);
	transpose(invectmp,s,invec,new_r_b,0,new_r_b,0,s);	
	
    }
    else{  //new_r_n > basecase
      std::cout << " right part break down  "<< std::endl;
    }
    
    

    

 ///////left part////////////////////////////////////////
  new_l_n = r;
  new_l_l = l_r;
  new_l_m = m_r;
  new_l_a2 = m_s;  //for tft
  new_l_a = n/new_l_n;// temporary for fft
  new_l_b = 1;
  //new_l_a2 = relax_size/new_l_l;

  /*
  std::cout << " new_l_n is  "<< new_l_n  << "  new_l_l is "<< new_l_l<<  " new_l_m is  "<<new_l_m<<
    "  new_l_a2  is "<<new_l_a2<< " new_l_a is   "<< new_l_a<< "new_l_b is  "<<  new_l_b  << std::endl;
  */

    if( new_l_n <= basecase ){
      //sfixn *invectmp_left = (sfixn *)calloc(new_l_n, sizeof(sfixn));
	int c=oldn,d=0;
	int left_step= (log2 (new_l_n));
	if(new_l_n<oldn){

	  for(int i=0;i< (e - left_step);i++){
	    d=d+c;
	    c=c>>1;
	  }

	} 
      //std::cout << " left part cal fft  "<< std::endl; 
	
	for(int j=0;j<new_l_a2;j++){
	  //cilk_for(int j=0;j<new_l_a2;j++){
	  int pii=j*new_l_n;
	  
	  
	/*
	  for(int i=0;i<new_l_n;i++){
	  //invectmp_left[i]=invec[j*new_l_n+i];
	  invectmp[i+pii]=invec[pii+i];
	  }
	*/
	  
	  
	//sfixn *B = (sfixn *)calloc(new_l_n, sizeof(sfixn));
	  
	//std::cout << " d is   "<<d<< std::endl;


	if(new_l_n>32){
	  DFT_eff(new_l_n, left_step, invec+pii, KRT+d,SB1+pii);//KRT+n is for level 2 omega
	  Shuffle_tft(new_l_n, invec+pii,SB1+pii);
	}
	else{

	  TFT_Basecase(new_l_n, left_step, invec+pii, KRT+d,SB1+pii);//KRT+n is for level 2 omega
	}

	/*
	for(int i=0;i<new_l_n;i++){
	  //std::cout << " invec here "<< std::endl;
	  invec[pii+i]=invectmp[i+pii];	  
	  //std::cout << " invec  "<<j*new_l_n+i<<" is "<< invec[j*new_l_n+i] << std::endl;
	}
	*/
      }
      //DFT_eff(n, e, invec, KRT, SB1);
      
      //free(invectmp_left);
    }
    else{  //new_l_n > basecase

      # pragma cilk_grainsize = TFT_grainsize; 
      cilk_for(int i=0;i<new_l_a2;i++){
	//for(int i=0;i<new_l_a2;i++){

	TFT_Core(oldn, new_l_n, new_l_l, new_l_m, basecase, relax_size,KRT,e,SB1+i*new_l_n,invec+i*new_l_n,invectmp+i*new_l_n);


      }
      

    }


 
  }



  /**************************
  /    MontDivMod:
  /    calculate num1/num2 mod p;
  /    all the numbers here are in mont mode;
  /
  /***************************/
  sfixn MontDivMod(sfixn num1,sfixn a){

    sfixn u=a,v=4179340454199820289,p=4179340454199820289;//v is prime number
    sfixn x1=num1,x2=0;

    while(u!=1 && v!=1){
      while(u % 2 == 0 ){
	u=u>>1;
	if(x1 % 2 == 0){
	  x1=x1>>1;
	}
	else{
	  x1=(x1+p)>>1;
	}
      }

      while(v % 2 == 0 ){
	v=v/2;
	if(x2 % 2 == 0){
	  x2=x2/2;
	}
	else{
	  x2=(x2+p)/2;
	}
      }

      if(u>=v){
	u=u-v;
	x1=x1-x2;
      }
      else{
	v=v-u;
	x2=x2-x1;
      }
    }
    if(u==1){
      if(x1<p){
	return x1;
      }
      else{
	return x1+p;
      }
    }
    else{
      if(x2<p){
	return x2;
      }
      else{
	return x2+p;
      }
    }

  }

  /**************************
  /    normaltomont:
  /    sfixn p=4179340454199820289;
  /    sfixn R=inverseMod(RINV,p);
  /    std::cout << "R is "<<R<<std::endl;
  /    sfixn R2= MulMod(R,R<<RSFT,p);
  /
  /***************************/
  sfixn normaltomont(sfixn normalnum, sfixn montnum){ 
    montnum = MontMulModSpe_OPT3_AS_GENE_INLINE(normalnum,R2);
    montnum = montnum<<RSFT;
    //std::cout << "KRT[1] in mont mode is "<< montnum <<std::endl;
    //std::cout << "R2 is "<< R2 <<std::endl;
    return montnum;
  }


  void ITFT_DFT(int L2,int step_L2,sfixn *KRT,int d_L2,sfixn *SB, sfixn *SB2){

    Shuffle_tft(L2, SB,SB2);
    DFT_eff(L2, step_L2, SB, KRT+d_L2, SB2);
      
  }


  void ITFT_Core(int oldn, int L, sfixn zeta,sfixn onedivzeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec, sfixn *SB2){

    long int temp,temp2;	
    sfixn *SB = (sfixn *)calloc(oldn, sizeof(sfixn));
    sfixn p=4179340454199820289;

    /*
    sfixn mont_two, mont_invtwo;
    mont_two = normaltomont(2,mont_two);
    mont_invtwo = normaltomont(2089670227099910145,mont_invtwo);
    std::cout << "2 in mont is "<< mont_two <<std::endl;
    std::cout << "2^-1 in mont is "<< mont_invtwo <<std::endl;
    //2 in mont is 3458764513820540920;
    //2^-1 in mont is -9223372036854775808;

    return;*/

     /*
     sfixn R = inverseMod(RINV,p);
     std::cout << "R is "<<R<<std::endl;
     sfixn R_square = MulMod(R,R,p);
     std::cout << "R2 is "<<R_square<<std::endl;

     sfixn test1 = MontMulModSpe_OPT3_AS_GENE_INLINE(653598720217121106,3787527286618587136);
     std::cout << " tidl inverse is  "<<test1 <<std::endl;
     //R is 432345564227567615
     //R_square is 1684656853714315195

    long int abc=1,def=0;
    def =  MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[5]);//from mont mode to normal mode
    std::cout << "KRT[5] in mont mode is "<<KRT[5] <<std::endl;
    std::cout << "KRT[5] in normal mode is "<< def <<std::endl;
    sfixn montdef;
    montdef = normaltomont(def,montdef);
    std::cout << "KRT[5] in mont mode is "<< montdef <<std::endl;

    sfixn num1 = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[5]);
    sfixn num2 = inverseMod(num1,p);
    sfixn montnum2 = normaltomont(num2,montnum2);
    std::cout << "num 1  is " <<num1<<", num2 is "<<num2<<" inverse is "<<montnum2<<std::endl;

    return;*/

    /*
    for(int i=0; i<8;i++)
      std::cout << i <<" is "<<KRT[i]<<std::endl;

    sfixn num1 = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[5]);
    sfixn num2 = inverseMod(num1,p);
    sfixn montnum2 = normaltomont(num2,montnum2);
    std::cout << "num 1  is " <<num1<<", num2 is "<<num2<<" inverse is "<<montnum2<<std::endl;


    return;*/

    /*
    sfixn num1=0,num2=0,num3=0,montnum3=0,montnum4=0;
    //std::cout << KRT[3]<<" div "<<KRT[2] <<std::endl;
    num1 = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[1]);//from mont mode to normal mode
    std::cout << KRT[1]<<" div "<<KRT[2] <<std::endl;
    num2 = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[2]);//from mont mode to normal mode
    std::cout << num1 <<" div "<< num2 <<std::endl;
    num3 = DivMod(num1,num2,p);
    montnum3 = normaltomont(num3,montnum3);
    std::cout << num1<<" div "<<num2<<" = "<< num3 <<std::endl;
    std::cout << KRT[3]<<" div "<<KRT[2]<<" = "<< montnum3 <<std::endl;

    montnum4 =  MontDivMod(KRT[3],KRT[2]);
    std::cout << KRT[3]<<" div "<<KRT[2] <<" = "<< montnum4<<std::endl;

    return;*/
 
    if(L==2){
      //basecase
      if(n==2){
	temp = invec[0];
	invec[0] = AddMod(temp,MulMod(inverseMod(zeta,p),invec[1],p),p);
	invec[1] = SubMod(temp,MulMod(inverseMod(zeta,p),invec[1],p),p);
      }
      else if(n==1 && f==1 && z==2){
	temp = invec[0];
	invec[0] = SubMod(MulMod(2,temp,p),invec[1],p);
	invec[1] = MulMod(zeta,SubMod(temp,invec[1],p),p);
      }
      else if(n==1 && f==1 && z==1){
	temp = invec[0];
	invec[0] = MulMod(2,temp,p);
	invec[1] = MulMod(zeta,temp,p); 
      }
      else if(n==1 && f==0 && z==2){
	temp = invec[0];
	invec[0] = SubMod(MulMod(2,temp,p),invec[1],p);
      }
      else if(n==1 && f==0 && z==1){
	invec[0] = MulMod(2,invec[0],p); 
      }
      else if(n==0 && z==2){
	invec[0] = DivMod(AddMod(invec[0],invec[1],p),2,p);
      }
      else if(n==0 && z==1){
	invec[0] = DivMod(invec[0],2,p);
      }

    
      return;
    }

    //recursive case
    int l,L1,L2,n1,n2,z1,z2,f_quote,z2_quote,m,m_quote;
    l=log2 (L);
    L1=pow(2, floor(l/2.0));
    L2=pow(2, ceil(l/2.0));
    n2=n%L2;
    n1=floor(n/L2);
    z2=z%L2;
    z1=floor(z/L2);
    if(n2+f>0){
      f_quote=1;
    }
    else{
      f_quote=0;
    }
    if(z1>0){
      z2_quote=L2;
    }
    else{
      z2_quote=z2;
    }
    m=std::min(n2,z2);
    m_quote=std::max(n2,z2);

    /*
    std::cout<<" l is "<<l<<" ,L1 is "<<L1<<"  ,L2 is "<<L2<<" ,n1 is "<<n1<<"  ,n2 is "<<n2<<
      "  ,z1 is "<<z1<<"  ,z2 is "<<z2<<"  ,f_quote is " <<f_quote<<" ,z2_quote is "<<z2_quote<<
      " ,m is "<<m<<" ,m_quote is "<<m_quote<<std::endl;*/
    

    //new_zeta1 is zeta^L1 
    sfixn new_zeta1=zeta;
    for(int i=0;i<L1-1;i++){
      new_zeta1 = MulMod(new_zeta1,zeta,p);
    }
    //std::cout<<" zeta is "<<zeta<<" L1 is " <<L1<<" new_zeta1 is "<<new_zeta1<<std::endl;

    //new_zeta2 is W_L^u * zeta
    int d=0,c=oldn;
    int totalstep= (log2 (oldn));
    if(l<totalstep){
      for(int i=0;i< (totalstep-l);i++){
	d=d+c;
	c=c>>1;
      }
    } 
    //std::cout<<" d   is "<<d<<std::endl;
   
    

    //row transforms
    if(L2== pow(2, ceil(totalstep/2.0)) ){
      
      int step_L2 = log2 (L2);
      int d_L2=0,c_L2=oldn;
      //int totalstep= (log2 (oldn));
      if(step_L2<totalstep){
	for(int i=0;i< (totalstep-step_L2);i++){
	  d_L2=d_L2+c_L2;
	  c_L2=c_L2>>1;
	}
      } 
      
      for(int u=0;u<n1;u++){//u is indexing for row
	
	for(int i=0;i<L2;i++){//for each row, there are L2 items
	  SB[i]=invec[u*L2+i];	  
	  //std::cout<<" input invec "<< u*L2+i<<" is "  << invec[u*L2+i] <<std::endl;
	  }
	ITFT_DFT(L2,step_L2,KRT,d_L2,SB,SB2);
	
	invec[u*L2]=SB[0];
	//std::cout<<"output invec "<<u*L2<<" is "  << invec[u*L2] <<std::endl;
	for(int i=1;i<L2;i++){//for each row, there are L2 items	  
	  invec[u*L2+i]=SB[L2-i];
	  //std::cout<<"output invec "<<u*L2+i <<" is "  << invec[u*L2+i] <<std::endl;
	}
      }
    }
    else{
      for(int u=0;u<n1;u++){//u is indexing for row
	ITFT_Core(oldn,L2,new_zeta1,onedivzeta,L2,L2,0,KRT,invKRT,invec+u*L2,SB);
	//for(int i=0;i<L2;i++){//for each row, there are L2 items
	//invec[u*L2+i]=SB[i];
	//}
      }
      
    }
    //for(int i = 0; i < oldn; ++i   )
    //std::cout <<"after row transform  tft "<<i<<"  is  " <<invec[i] << std::endl;

    /*
    for(int i=0;i<oldn;i++){
      std::cout <<" invec "<< i<<" is "<<invec[i] << std::endl;
    }
    std::cout <<" row  <<<<<<<<<<<<<  " << std::endl;*/

    
    //rightmost column transforms
    for(int u=n2;u<m_quote;u++){// u is indexing for column
      for(int i=0;i<L1;i++){// for each col, there are L1 items
	SB[i]=invec[u+i*L2];
	//std::cout<<" rightmost 1 input invec "<< u+i*L2 <<" is "  << invec[u+i*L2] <<std::endl;
      }
      
      sfixn new_omega = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[d+u]);//w_L^u  in normal mode
      sfixn new_zeta2 =  MulMod(new_omega,zeta,p);//w_L^u * zeta 
      sfixn new_onedivzeta = MulMod(MontMulModSpe_OPT3_AS_GENE_INLINE(1,invKRT[d+u]),1,p);

      ITFT_Core(oldn,L1,new_zeta2,new_onedivzeta,z1+1,n1,f_quote,KRT,invKRT,SB,SB);//not sure for omega
 
      for(int i=0;i<L1;i++){// for each col, there are L1 items
	invec[u+i*L2]=SB[i];
	//std::cout<<"rightmost 1 input invec "<< u+i*L2 <<" is "  << invec[u+i*L2] <<std::endl;
      }

    }


    for(int u=m_quote;u<z2_quote;u++){// u is indexing for column
      for(int i=0;i<L1;i++){// for each col, there are L1 items
	SB[i]=invec[u+i*L2];
	//std::cout<<" rightmost 2 SB_"<<i<<" is "<<SB[i]<<std::endl;
      }
 

      sfixn new_omega = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[d+u]);//w_L^u  in normal mode
      sfixn new_zeta2 =  MulMod(new_omega,zeta,p);//w_L^u * zeta 
      sfixn new_onedivzeta = MulMod(MontMulModSpe_OPT3_AS_GENE_INLINE(1,invKRT[d+u]),1,p);

      ITFT_Core(oldn,L1,new_zeta2,new_onedivzeta,z1,n1,f_quote,KRT,invKRT,SB,SB);
   
      for(int i=0;i<L1;i++){// for each col, there are L1 items
	invec[u+i*L2]=SB[i];
	//std::cout<<"AFTER rightmost 2 SB_"<<i<<" is "<<SB[i]<<std::endl;
      }
      
    }
    // for(int i = 0; i < oldn; ++i   )
    //std::cout <<"after rightmost  tft "<<i<<"  is  " <<invec[i] << std::endl;

    //return;
    //std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "<<std::endl;

    /*
    for(int i=0;i<oldn;i++){
      std::cout <<" invec "<< i<<" is "<<invec[i] << std::endl;
    }
    std::cout <<" right most  <<<<<<<<<<<<<  " << std::endl;*/
   

    //last row transforms
    if(f_quote==1){

      /*
      for(int i=0;i<L2;i++){
	SB[i]=invec[n1*L2+i];
	//std::cout<<" last row SB_"<<i<<" is "<<SB[i]<<std::endl;
	//std::cout<<" last row SB_  [n1*L2+i]  is "<<n1*L2+i <<std::endl;
      }
      //for(int i = 0; i < oldn; ++i   )
      // std::cout <<"before last row  tft "<<i<<"  is  " <<invec[i] << std::endl;

     sfixn *invectmp = (sfixn *)calloc(oldn, sizeof(sfixn));
     for(int i = 0; i < oldn; ++i)
       invectmp[i]=invec[i];
     

      */
      ITFT_Core(oldn,L2,new_zeta1,onedivzeta,z2_quote,n2,f,KRT,invKRT,invec+n1*L2,SB);
     
     /*
     for(int i = 0; i < oldn; ++i)
       invec[i]=invectmp[i];

     //for(int i = 0; i < oldn; ++i)
     //std::cout <<"middle last row  tft "<<i<<"  is  " <<SB[i] << std::endl;
     

      for(int i=0;i<L2;i++){
	invec[n1*L2+i]=SB[i];
	//std::cout<<"AFTER last row SB_"<<i<<" is "<<SB[i]<<std::endl;
	//std::cout<<"AFTER last row SB_  [n1*L2+i]  is "<<n1*L2+i <<std::endl;
      }
      free(invectmp);
     */
    }

 

    
    //leftmost column transforms
    for(int u=0;u<m;u++){
      for(int i=0;i<L1;i++){
	SB[i]=invec[u+i*L2];
	//std::cout<<" leftmost 1 SB_"<<i<<" is "<<SB[i]<<std::endl;
      }
     
      sfixn new_omega = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[d+u]);
      sfixn new_zeta2 = MulMod(new_omega,zeta,p);
      ITFT_Core(oldn,L1,new_zeta2,onedivzeta,z1+1,n1+1,0,KRT,invKRT,SB,SB);//not sure for omega
      for(int i=0;i<L1;i++){// for each col, there are L1 items
	invec[u+i*L2]=SB[i];
	//std::cout<<"AFTER left most 1 SB_"<<i<<" is "<<SB[i]<<std::endl;
      }
    }
 
    for(int u=m;u<n2;u++){
      for(int i=0;i<L1;i++){
	SB[i]=invec[u+i*L2];	
	//std::cout<<"left most 2 SB_"<<i<<" is "<<SB[i]<<std::endl;
      }
     
      sfixn new_omega = MontMulModSpe_OPT3_AS_GENE_INLINE(1,KRT[d+u]);
      sfixn new_zeta2 = MulMod(new_omega,zeta,p);
      ITFT_Core(oldn,L1,new_zeta2,onedivzeta,z1,n1+1,0,KRT,invKRT,SB,SB);//not sure for omega
      for(int i=0;i<L1;i++){// for each col, there are L1 items
	invec[u+i*L2]=SB[i];
	//std::cout<<"AFTER left most 2 SB_"<<i<<" is "<<SB[i]<<std::endl;
      }
      }

    /*
    for(int i=0;i<oldn;i++){
      std::cout <<" invec "<< i<<" is "<<invec[i] << std::endl;
    }
    std::cout <<" last row  <<<<<<<<<<<<<  " << std::endl;*/

    free(SB);


  }

  void ITFT_Wrapper(int oldn, int L, sfixn zeta,sfixn onedivzeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec, sfixn *SB,int basecase, sfixn *invectmp, int relax_size){

    /*
  ITFT_Core(oldn,L,zeta,onedivzeta,z,n,f,KRT,invKRT,invec,SB2);
  for (int i = 0; i < oldn; ++i){
    invec[i] = DivMod(invec[i],oldn,4179340454199820289);
    //std::cout <<" after tft invec "<<i<<"  is  " << invec[i] << std::endl;
    }*/

   
    //sfixn trynum = MontMulModSpe_OPT3_AS_GENE_INLINE(RINV,2701251764010509990);
    //std::cout <<" A inverse   is  " <<  trynum  << std::endl;
    //return;
     
    /*
    for (int i = 0; i < oldn; ++i){
      invec[i] = i+1;
    }
    for (int i = 0; i < oldn; ++i){
      std::cout << invec[i] << std::endl;
    }
    PBPAS::sqtranspose(invec,4,0,4,0,1,4);
    std::cout << "<<<<<<<<<<<<<" << std::endl;
    for (int i = 0; i < oldn; ++i){
      std::cout << invec[i] << std::endl;
    }
    return;*/

  sfixn p=4179340454199820289;
  sfixn mont_one = 1729382256910270460;
  int s,r,l_s,m_s,l_r,m_r;//for break down TFT
  int new_r_n, new_r_l, new_r_m, new_r_a, new_r_b;//for right part
  int new_l_n, new_l_l, new_l_m, new_l_a, new_l_b, new_l_a2;//for left part
  int bigern;
  bigern = pow(2, ceil(log2 (z)));
  sfixn mont_2 = normaltomont(2,mont_2);
  sfixn mont_inv2 = normaltomont(inverseMod(2,p),mont_inv2);


  //sfixn *SB2 = (sfixn *)calloc(3*oldn, sizeof(sfixn));
  //sfixn *SB3 = (sfixn *)calloc(3*oldn, sizeof(sfixn));
 

  //unsigned long start = __cilkview_getticks();
  //ITFT_Core(oldn,L,zeta,onedivzeta,z,n,f,KRT,invKRT,invec,SB); //compute in normal mode 
  Mont_ITFT_Core(oldn,L,0,z,n,f,KRT,invKRT,invec,mont_inv2,SB,p);//compute in montgemory mode

  //unsigned long end = __cilkview_getticks();
  //std::cout << L <<"&" << (end - start) / 1000.f << std::endl;

  //free(SB2);
  //free(SB3);
 
  /*
  std::ofstream myfile;
  //myfile.open ("itft_result_right.txt"); // result in normal mode
  myfile.open ("itft_result_wrong.txt");// result in mont mode
 
  for(int i = 0; i < z; ++i) {
    myfile << invec[i]<<"\n";

    }
  

    myfile.close();*/

  /*
  // only need for verify
  for (int i = 0; i < bigern; ++i){
	invec[i] = DivMod(invec[i],bigern,4179340454199820289);
	//std::cout <<" after tft invec "<<i<<"  is  " << invectmp[i] << std::endl;
	}
  */


  return;

}

  void Mont_ITFT_Core(int oldn, int L, int index_zeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec,  sfixn mont_inv2, sfixn *SB, sfixn p){

    long int temp,temp2;	
    //sfixn *SB3 = (sfixn *)calloc(oldn, sizeof(sfixn));
    //sfixn p=4179340454199820289;
    //sfixn Mont_invtwo = -9223372036854775808;


    if(L==2){
      //basecase
      if(n==2){
	temp = invec[0];
	invec[0] = AddModSpe(temp,MontMulModSpe_OPT3_AS_GENE_INLINE(invKRT[index_zeta],invec[1]));
	invec[1] = SubModSpe(temp,MontMulModSpe_OPT3_AS_GENE_INLINE(invKRT[index_zeta],invec[1]));

	//std::cout <<"        11111111          " << std::endl;

      }
      else if(n==1 && f==1 && z==2){
	temp = invec[0];
	invec[0] = SubModSpe(MontMulModSpe_OPT3_AS_GENE_INLINE(Mont_two,temp),invec[1]);
 	invec[1] = MontMulModSpe_OPT3_AS_GENE_INLINE(KRT[index_zeta],SubModSpe(temp,invec[1]));

	//std::cout <<"       222222222         " << std::endl;
      }
      else if(n==1 && f==1 && z==1){
	temp = invec[0];
	invec[0] = MontMulModSpe_OPT3_AS_GENE_INLINE(Mont_two,temp);
	invec[1] = MontMulModSpe_OPT3_AS_GENE_INLINE(KRT[index_zeta],temp); 

	//std::cout <<"        333333333          " << std::endl;
      }
      else if(n==1 && f==0 && z==2){
	temp = invec[0];
	invec[0] = SubModSpe(MontMulModSpe_OPT3_AS_GENE_INLINE(Mont_two,temp),invec[1]);

	//std::cout <<"        444444444         " << std::endl;
      }
      else if(n==1 && f==0 && z==1){
	invec[0] = MontMulModSpe_OPT3_AS_GENE_INLINE(Mont_two,invec[0]); 

	//std::cout <<"       5555555555         " << std::endl;
      }
      else if(n==0 && z==2){

	invec[0] = MontMulModSpe_OPT3_AS_GENE_INLINE(AddModSpe(invec[0],invec[1]),  mont_inv2 );
	
	//std::cout <<"       6666666666       " << std::endl;
      }
      else if(n==0 && z==1){
	invec[0] = MontMulModSpe_OPT3_AS_GENE_INLINE(invec[0], mont_inv2);
	
	//std::cout <<"       7777777777       " << std::endl;
      }

    
      return;
    }

    //recursive case
    int l,L1,L2,n1,n2,z1,z2,f_quote,z2_quote,m,m_quote;
    l=log2 (L);
    L1=pow(2, floor(l/2.0));
    L2=pow(2, ceil(l/2.0));
    n2=n%L2;
    n1=floor(n/L2);
    z2=z%L2;
    z1=floor(z/L2);
    if(n2+f>0){
      f_quote=1;
    }
    else{
      f_quote=0;
    }
    if(z1>0){
      z2_quote=L2;
    }
    else{
      z2_quote=z2;
    }
    m=std::min(n2,z2);
    m_quote=std::max(n2,z2);

    /*
    std::cout<<" l is "<<l<<" ,L1 is "<<L1<<"  ,L2 is "<<L2<<" ,n1 is "<<n1<<"  ,n2 is "<<n2<<
      "  ,z1 is "<<z1<<"  ,z2 is "<<z2<<"  ,f_quote is " <<f_quote<<" ,z2_quote is "<<z2_quote<<
      " ,m is "<<m<<" ,m_quote is "<<m_quote<<std::endl;*/
    


    int totalstep= (log2 (oldn));
    int new_index_zeta1;
    if(index_zeta==0){
      new_index_zeta1 = 0;
    }
    else
      {
	new_index_zeta1 = (index_zeta*L1) % oldn;	
      }

    int temp_zeta2 = pow(2, (totalstep - l )) ;
 
    
    //row transforms
    if(L2 == pow(2, ceil(totalstep/2.0)) ){

      int step_L2 = log2 (L2);
      int d_L2=0,c_L2=oldn;
      if(step_L2<totalstep){
	for(int i=0;i< (totalstep-step_L2);i++){
	  d_L2=d_L2+c_L2;
	  c_L2=c_L2>>1;
	}
      } 

    
      sfixn *SB3 = (sfixn *)calloc(L2*n1, sizeof(sfixn));

      # pragma cilk_grainsize = TFT_grainsize;
      cilk_for(int u=0;u<n1;u++){//u is indexing for row
	//for(int u=0;u<n1;u++){//u is indexing for row
	
	//sfixn *SB3 = (sfixn *)calloc(L2, sizeof(sfixn));//for parallel
	
	for(int i=0;i<L2;i++){//for each row, there are L2 items
	  SB[u*L2+i]=invec[u*L2+i];	  
	}
	ITFT_DFT(L2,step_L2,KRT,d_L2,SB+u*L2,SB3+u*L2);
	invec[u*L2]=SB[0+u*L2];
	for(int i=1;i<L2;i++){//for each row, there are L2 items	  
	  invec[u*L2+i]=SB[(u+1)*L2-i];
	}
	
	
	//free(SB3);//for parallel
      }
      free(SB3);
  
    }
    else{

      # pragma cilk_grainsize = TFT_grainsize;
      cilk_for(int u=0;u<n1;u++){//u is indexing for row  
	  
	//Mont_ITFT_Core(oldn,L2,new_index_zeta1,L2,L2,0,KRT,invKRT,invec+u*L2,mont_inv2,SB);//for serial
	Mont_ITFT_Core(oldn,L2,new_index_zeta1,L2,L2,0,KRT,invKRT,invec+u*L2,mont_inv2,SB+u*L2,p);
      }
      
  
    }

    /*
    if(L1==L2){  
      //std::cout <<"   L1==L2      " << std::endl;
      PBPAS::sqtranspose(invec,L2,0,L1,0,1,L2);
      //sfixn *SB4 = (sfixn *)calloc(L1, sizeof(sfixn));
      
      //rightmost column transforms
      cilk_for(int u=n2;u<m_quote;u++){// u is indexing for column		   
	
	int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;
	Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1+1,n1,f_quote,KRT,invKRT,invec+u*L1,mont_inv2,SB+u*L1,p);
 
      }
      
      cilk_for(int u=m_quote;u<z2_quote;u++){// u is indexing for column
	
	int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;      
	Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1,n1,f_quote,KRT,invKRT,invec+u*L1,mont_inv2,SB+u*L1,p);
      
      }
      
      
    //last row transforms
      if(f_quote==1){
	
      sfixn *temp_row = (sfixn *)calloc(L2<<2, sizeof(sfixn));
      for(int i=0;i<L2;i++){
	temp_row[i]=invec[n1+i*L1];
      }
      
      Mont_ITFT_Core(oldn,L2,new_index_zeta1,z2_quote,n2,f,KRT,invKRT,temp_row,mont_inv2,SB,p);
      
      for(int i=0;i<L2;i++){
	invec[n1+i*L1]=temp_row[i];
      }
      
      free(temp_row);
      
      }
      
      //leftmost column transforms
      cilk_for(int u=0;u<m;u++){
        
	int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;
	
	Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1+1,n1+1,0,KRT,invKRT,invec+u*L1,mont_inv2,SB+u*L1,p);//not sure for omega
	
    }
 
      cilk_for(int u=m;u<n2;u++){
 
      int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;
 
      Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1,n1+1,0,KRT,invKRT,invec+u*L1,mont_inv2,SB+u*L1,p);//not sure for omega
      
    }
    
    //transpose_serial(SB,L1,invec,L2,0,L2,0,L1);	
    PBPAS::sqtranspose(invec,L1,0,L2,0,1,L1);
    //free(SB4);
    
    }
    else{*/


    //transpose_serial(invec,L2,SB,L1,0,L1,0,L2);//for serial
    transpose(invec,L2,SB,L1,0,L1,0,L2);

    //rightmost column transforms

    # pragma cilk_grainsize = TFT_grainsize;
    cilk_for(int u=n2;u<m_quote;u++){// u is indexing for column		   
	
      int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;
      //Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1+1,n1,f_quote,KRT,invKRT,SB+u*L1,mont_inv2,invec);//for serial
      Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1+1,n1,f_quote,KRT,invKRT,SB+u*L1,mont_inv2,invec+u*L1,p);
    }

    # pragma cilk_grainsize = TFT_grainsize;
    cilk_for(int u=m_quote;u<z2_quote;u++){// u is indexing for column
    
      int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;
      //Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1,n1,f_quote,KRT,invKRT,SB+u*L1,mont_inv2,invec);//for serial
      Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1,n1,f_quote,KRT,invKRT,SB+u*L1,mont_inv2,invec+u*L1,p);
    }


    //last row transforms
    if(f_quote==1){
      
      sfixn *temp_row = (sfixn *)calloc(L2<<2, sizeof(sfixn));
      for(int i=0;i<L2;i++){
	temp_row[i]=SB[n1+i*L1];
      }

      Mont_ITFT_Core(oldn,L2,new_index_zeta1,z2_quote,n2,f,KRT,invKRT,temp_row,mont_inv2,invec,p);

      for(int i=0;i<L2;i++){
	SB[n1+i*L1]=temp_row[i];
      }

      free(temp_row);
      
    }
    
    //leftmost column transforms

    # pragma cilk_grainsize = TFT_grainsize;
    cilk_for(int u=0;u<m;u++){
        
      int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;

      //Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1+1,n1+1,0,KRT,invKRT,SB+u*L1,mont_inv2,invec);//not sure for omega //for serial
      Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1+1,n1+1,0,KRT,invKRT,SB+u*L1,mont_inv2,invec+u*L1,p);//not sure for omega 
    }
 
    # pragma cilk_grainsize = TFT_grainsize;
    cilk_for(int u=m;u<n2;u++){
 
      int index_new_zeta2 = ( temp_zeta2 * u + index_zeta) % oldn;
 
      //Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1,n1+1,0,KRT,invKRT,SB+u*L1,mont_inv2,invec);//not sure for omega//for serial
      Mont_ITFT_Core(oldn,L1,index_new_zeta2,z1,n1+1,0,KRT,invKRT,SB+u*L1,mont_inv2,invec+u*L1,p);//not sure for omega
      }
    
    //transpose_serial(SB,L1,invec,L2,0,L2,0,L1); //for serial	
    transpose(SB,L1,invec,L2,0,L2,0,L1);

    //} 
 
  }


void InvDFT_eff(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn){
	DFT_eff(n,r,A,W,B);
	for(int i=0; i<n; i++) {
		A[i]=MontMulModSpe_OPT3_AS_GENE_INLINE(A[i], invn);
	}
}

}


   /*  //test shuffle function
    for(int i=0;i<n;i++){
      A[i]=i+1;
    }
    sfixn *B = (sfixn *)calloc(n, sizeof(sfixn)); 
    for(int i=0;i<n;i++){
      B[i]=0;
    } 

    Shuffle(n,A,B);

    for(int i=0;i<n;i++){
     
      std::cout << " after shuffle A "<<i<<"  is  "<< A[i] << std::endl;
    }
  
    free(B);
    */



/*
parallel to serial
make serial

AND

transpose to transpose_serial    both in TFT_core and Mont_itft_core


 */
