

#include <string.h>

//foward declares
static void Shuffle2(int n, long long int* A,long long int* B);
static void Shuffle(int n, long long int* A,long long int* B);
static long long int testDFT(int n,int num,long long int* A,long long int *W, long long int INV_PRIME, long long int MY_PRIME);

static void extended_gcd(long long int a, long long int b, long long int* rr, long long int* uu, long long int* vv) {
	long long int r = a;
	long long int u = 1;
	long long int v = 0;
	long long int r_ = b;
	long long int u_ = 0;
	long long int v_ = 1;
	long long int q, tr, tu, tv;
	while (r_!=0) {
		q = r/r_;
		tr = r_;
		tu = u_;
		tv = v_;
		r_ = r - q*r_;
		u_ = u - q*u_;
		v_ = v - q*v_;
		r = tr;
		u = tu;
		v = tv;
	}

	if (rr != NULL) {
		*rr = r;
	}
	if (uu != NULL) {
		*uu = u;
	}
	if (vv != NULL) {
		*vv = v;
	}
}

static inline void  AddSubSpecSSEModInplace(long long int* a,long long int* b, long long int MY_PRIME){
	long long int u1 = smallprimefield_AddMod(*a,*b, MY_PRIME);
	long long int u2 = smallprimefield_SubMod(*a,*b, MY_PRIME);
	long long int u3 = smallprimefield_AddMod(*(a+1),*(b+1), MY_PRIME);
	long long int u4 = smallprimefield_SubMod(*(a+1),*(b+1), MY_PRIME);
	*a = u1;
	*(a+1) = u3;
	*b = u2;
	*(b+1) = u4;
	return;
}
static inline
void  AddSubSpecSSEMod(long long int* r1,long long int* r2,long long int* a,long long int* b, long long int MY_PRIME){
	*r1 = smallprimefield_AddMod(*a,*b, MY_PRIME);
	*r2 = smallprimefield_SubMod(*a,*b, MY_PRIME);
	*(r1+1) = smallprimefield_AddMod(*(a+1),*(b+1), MY_PRIME);
	*(r2+1) = smallprimefield_SubMod(*(a+1),*(b+1), MY_PRIME);
   return;
}
static void inline FFT_2POINT(long long int *A,long long int *W, long long int MY_PRIME){
	long long int u = A[0];
	long long int t = A[1];
	A[0] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[1] = smallprimefield_SubMod(u,t,MY_PRIME);
}
static void inline FFT_4POINT(long long int *A,long long int *W, long long int INV_PRIME, long long int MY_PRIME){
	long long int *Wp = W + (4<<1)-4;
	long long int w = A[1];
	A[1] = A[2];
	A[2] = w;
	long long int u = A[0];
	long long int t = A[1];
	A[0] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[1] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[2];
	t = A[3];
	A[2] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[3] = smallprimefield_SubMod(u,t,MY_PRIME);
	A[3] = smallprimefield_Mul_INLINE(A[3],*(Wp-3),INV_PRIME,MY_PRIME);

	AddSubSpecSSEModInplace(A,A+2,MY_PRIME);
}
static void inline FFT_8POINT(long long int *A,long long int *W, long long int INV_PRIME, long long int MY_PRIME){
	long long int *Wp = W + (8<<1)-4;
	long long int u = A[0];
	long long int t = A[1];
	A[0] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[1] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[2];
	t = A[3];
	A[2] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[3] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[4];
	t = A[5];
	A[4] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[5] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[6];
	t = A[7];
	A[6] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[7] = smallprimefield_SubMod(u,t,MY_PRIME);
	A[3] = smallprimefield_Mul_INLINE(A[3],*(Wp-3),INV_PRIME,MY_PRIME);
	A[7] = smallprimefield_Mul_INLINE(A[7],*(Wp-3),INV_PRIME,MY_PRIME);

	AddSubSpecSSEModInplace(A,A+2,MY_PRIME);
	AddSubSpecSSEModInplace(A+4,A+6,MY_PRIME);
	A[5] = smallprimefield_Mul_INLINE(A[5],*(Wp-11),INV_PRIME,MY_PRIME);
	A[6] = smallprimefield_Mul_INLINE(A[6],*(Wp-10),INV_PRIME,MY_PRIME);
	A[7] = smallprimefield_Mul_INLINE(A[7],*(Wp-9),INV_PRIME,MY_PRIME);
	AddSubSpecSSEModInplace(A,A+4,MY_PRIME);
	AddSubSpecSSEModInplace(A+2,A+6,MY_PRIME);
}
static void inline FFT_16POINT(long long int *A,long long int *W,long long int INV_PRIME, long long int MY_PRIME){
	long long int *Wp = W + (16<<1)-4;
	long long int u = A[0];
	long long int t = A[1];
	A[0] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[1] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[2];
	t = A[3];
	A[2] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[3] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[4];
	t = A[5];
	A[4] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[5] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[6];
	t = A[7];
	A[6] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[7] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[8];
	t = A[9];
	A[8] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[9] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[10];
	t = A[11];
	A[10] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[11] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[12];
	t = A[13];
	A[12] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[13] = smallprimefield_SubMod(u,t,MY_PRIME);
	u = A[14];
	t = A[15];
	A[14] = smallprimefield_AddMod(u,t,MY_PRIME);
	A[15] = smallprimefield_SubMod(u,t,MY_PRIME);
	A[3] = smallprimefield_Mul_INLINE2(A[3],*(Wp-3),INV_PRIME,MY_PRIME);
	A[7] = smallprimefield_Mul_INLINE2(A[7],*(Wp-3),INV_PRIME,MY_PRIME);
	A[11] = smallprimefield_Mul_INLINE2(A[11],*(Wp-3),INV_PRIME,MY_PRIME);
	A[15] = smallprimefield_Mul_INLINE2(A[15],*(Wp-3),INV_PRIME,MY_PRIME);

	AddSubSpecSSEModInplace(A,A+2,MY_PRIME);
	AddSubSpecSSEModInplace(A+4,A+6,MY_PRIME);
	AddSubSpecSSEModInplace(A+8,A+10,MY_PRIME);
	AddSubSpecSSEModInplace(A+12,A+14,MY_PRIME);
	A[5] = smallprimefield_Mul_INLINE2(A[5],*(Wp-11),INV_PRIME,MY_PRIME);
	A[6] = smallprimefield_Mul_INLINE2(A[6],*(Wp-10),INV_PRIME,MY_PRIME);
	A[7] = smallprimefield_Mul_INLINE2(A[7],*(Wp-9),INV_PRIME,MY_PRIME);
	A[13] = smallprimefield_Mul_INLINE2(A[13],*(Wp-11),INV_PRIME,MY_PRIME);
	A[14] = smallprimefield_Mul_INLINE2(A[14],*(Wp-10),INV_PRIME,MY_PRIME);
	A[15] = smallprimefield_Mul_INLINE2(A[15],*(Wp-9),INV_PRIME,MY_PRIME);

	AddSubSpecSSEModInplace(A,A+4,MY_PRIME);
	AddSubSpecSSEModInplace(A+2,A+6,MY_PRIME);
	AddSubSpecSSEModInplace(A+8,A+12,MY_PRIME);
	AddSubSpecSSEModInplace(A+10,A+14,MY_PRIME);

	A[9] = smallprimefield_Mul_INLINE2(A[9],*(Wp-27),INV_PRIME,MY_PRIME);
	A[10] = smallprimefield_Mul_INLINE2(A[10],*(Wp-26),INV_PRIME,MY_PRIME);
	A[11] = smallprimefield_Mul_INLINE2(A[11],*(Wp-25),INV_PRIME,MY_PRIME);
	A[12] = smallprimefield_Mul_INLINE2(A[12],*(Wp-24),INV_PRIME,MY_PRIME);
	A[13] = smallprimefield_Mul_INLINE2(A[13],*(Wp-23),INV_PRIME,MY_PRIME);
	A[14] = smallprimefield_Mul_INLINE2(A[14],*(Wp-22),INV_PRIME,MY_PRIME);
	A[15] = smallprimefield_Mul_INLINE2(A[15],*(Wp-21),INV_PRIME,MY_PRIME);
	AddSubSpecSSEModInplace(A,A+8,MY_PRIME);
	AddSubSpecSSEModInplace(A+2,A+10,MY_PRIME);
	AddSubSpecSSEModInplace(A+4,A+12,MY_PRIME);
	AddSubSpecSSEModInplace(A+6,A+14,MY_PRIME);
}

DFT_ITERATIVE

static void Shuffle(int n, long long int *A, long long int *B){
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
	long long int *A2 = A + n2;
	memcpy(A2, B, n2*(sizeof(long long int)));
}

static long long int testDFT(int n,int num,long long int* A,long long int *W, long long int INV_PRIME, long long int MY_PRIME){
	long long int res = A[0];
	long long int w = *(W+num);
	long long int root = w;
	for(long long int k=1;k<n;k++){

		long long int t = smallprimefield_Mul_INLINE(A[k],*(W+((num*k)%n)),INV_PRIME,MY_PRIME);
		res = smallprimefield_AddMod(res,t,MY_PRIME);
	}
	return res;
}

void DFT_rec(int n, int r,long long int *A,long long int *W,long long int *B, long long int INV_PRIME, long long int MY_PRIME){
	if (n==1) {
		return;
	}
	else if(n==FFT_THRESHOLD){
		ArrayBitReversalSpec(A);
		DFT_iter(A,W,INV_PRIME,MY_PRIME);
		return;
	}
	else{
		int n2 = n>>1;
		int r1 = r-1;
		//Shuffle(n, A, B);
		long long int *W2 = W+n;
		DFT_rec(n2, r1, A, W2, B, INV_PRIME, MY_PRIME);
		DFT_rec(n2, r1, A+n2, W2, B, INV_PRIME, MY_PRIME);

		for (int k=n2; k<(n2<<1); k+=8){
			A[k+0]=smallprimefield_Mul_INLINE(A[k+0],W[k-n2+0],INV_PRIME,MY_PRIME);
			A[k+1]=smallprimefield_Mul_INLINE(A[k+1],W[k-n2+1],INV_PRIME,MY_PRIME);
			A[k+2]=smallprimefield_Mul_INLINE(A[k+2],W[k-n2+2],INV_PRIME,MY_PRIME);
			A[k+3]=smallprimefield_Mul_INLINE(A[k+3],W[k-n2+3],INV_PRIME,MY_PRIME);
			A[k+4]=smallprimefield_Mul_INLINE(A[k+4],W[k-n2+4],INV_PRIME,MY_PRIME);
			A[k+5]=smallprimefield_Mul_INLINE(A[k+5],W[k-n2+5],INV_PRIME,MY_PRIME);
			A[k+6]=smallprimefield_Mul_INLINE(A[k+6],W[k-n2+6],INV_PRIME,MY_PRIME);
			A[k+7]=smallprimefield_Mul_INLINE(A[k+7],W[k-n2+7],INV_PRIME,MY_PRIME);
			long long int t = A[k-n2];
			long long int u = A[k];
			A[k-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+1-n2];
			u = A[k+1];
			A[k+1-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+1] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+2-n2];
			u = A[k+2];
			A[k+2-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+2] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+3-n2];
			u = A[k+3];
			A[k+3-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+3] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+4-n2];
			u = A[k+4];
			A[k+4-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+4] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+5-n2];
			u = A[k+5];
			A[k+5-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+5] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+6-n2];
			u = A[k+6];
			A[k+6-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+6] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+7-n2];
			u = A[k+7];
			A[k+7-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+7] = smallprimefield_SubMod(t,u,MY_PRIME);
		}
	}
}
static void Shuffle2(int n, long long int *A, long long int *B){
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
	long long int *A2 = A + n2;
	memcpy(A2, B, 3*n2*(sizeof(long long int)));
}

void DFT_eff(int n, int r,long long int *A,long long int *W,long long int *B, const Prime_ptr* Pptr){
	#if DEBUG
	long long int test = testDFT(n,1,A,W);
	#endif
	if (n==1) {
		return;
	}

	long long int MY_PRIME = Pptr->prime;
	long long int INV_PRIME = Pptr->prime_inv;

	if(n==FFT_THRESHOLD){
		ArrayBitReversalSpec(A);
		DFT_iter(A,W,INV_PRIME,MY_PRIME);
		return;
	}
	else if(n>FFT_THRESHOLD){
		long n2 = n>>1;
		long r1 = r-1;
		long long int twoFFTThreshold = FFT_THRESHOLD<<1;
		if(n>(twoFFTThreshold<<1)){
			long long int pos = 0;
			long long int m=n;

			Shuffle(n,A,B);
			m>>=1;
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
		}
		else if(n==(twoFFTThreshold<<1)){
			Shuffle(n,A,B);
			long long int m = n>>1;
			Shuffle(m,A,B);
			Shuffle(m,A+m,B);
		}
		else{
			Shuffle(n,A,B);
		}
		//Shuffle(n, A, B);
		long long int *W2 = W+n;
		DFT_rec(n2, r1, A, W2, B, INV_PRIME, MY_PRIME);
		DFT_rec(n2, r1, A+n2, W2, B, INV_PRIME, MY_PRIME);

		for (int k=n2; k<(n2<<1); k+=8){
			A[k+0]=smallprimefield_Mul_INLINE(A[k+0],W[k-n2+0],INV_PRIME,MY_PRIME);
			A[k+1]=smallprimefield_Mul_INLINE(A[k+1],W[k-n2+1],INV_PRIME,MY_PRIME);
			A[k+2]=smallprimefield_Mul_INLINE(A[k+2],W[k-n2+2],INV_PRIME,MY_PRIME);
			A[k+3]=smallprimefield_Mul_INLINE(A[k+3],W[k-n2+3],INV_PRIME,MY_PRIME);
			A[k+4]=smallprimefield_Mul_INLINE(A[k+4],W[k-n2+4],INV_PRIME,MY_PRIME);
			A[k+5]=smallprimefield_Mul_INLINE(A[k+5],W[k-n2+5],INV_PRIME,MY_PRIME);
			A[k+6]=smallprimefield_Mul_INLINE(A[k+6],W[k-n2+6],INV_PRIME,MY_PRIME);
			A[k+7]=smallprimefield_Mul_INLINE(A[k+7],W[k-n2+7],INV_PRIME,MY_PRIME);
			long long int t = A[k-n2];
			long long int u = A[k];
			A[k-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+1-n2];
			u = A[k+1];
			A[k+1-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+1] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+2-n2];
			u = A[k+2];
			A[k+2-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+2] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+3-n2];
			u = A[k+3];
			A[k+3-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+3] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+4-n2];
			u = A[k+4];
			A[k+4-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+4] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+5-n2];
			u = A[k+5];
			A[k+5-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+5] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+6-n2];
			u = A[k+6];
			A[k+6-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+6] = smallprimefield_SubMod(t,u,MY_PRIME);
			t = A[k+7-n2];
			u = A[k+7];
			A[k+7-n2] = smallprimefield_AddMod(t,u,MY_PRIME);
			A[k+7] = smallprimefield_SubMod(t,u,MY_PRIME);
		}
	}
	else{
		switch (n){
			case 2:
				FFT_2POINT(A,W,MY_PRIME);
				break;
			case 4:
				FFT_4POINT(A,W,INV_PRIME,MY_PRIME);
				break;
			case 8:
				ArrayBitReversalSpec8(A);
				FFT_8POINT(A,W,INV_PRIME,MY_PRIME);
				break;
			case 16:
				ArrayBitReversalSpec16(A);
				FFT_16POINT(A,W,INV_PRIME,MY_PRIME);
				break;
			FFT_CASES
			default:
				return;
		}
	}
}

void InvDFTKeepMont_eff(int n, int r, long long int *A,long long int *W,  long long int *B,long long int invn, const Prime_ptr* Pptr){

	DFT_eff(n,r,A,W,B,Pptr);

	long long int MY_PRIME = Pptr->prime;
	long long int INV_PRIME = Pptr->prime_inv;

	for(int i=0; i<n; i++) {
		A[i]=smallprimefield_Mul_INLINE(A[i], invn, INV_PRIME, MY_PRIME);
	}
}
void InvDFT_eff(int n, int r, long long int *A,long long int *W, long long int *B,long long int invn, const Prime_ptr* Pptr){

	DFT_eff(n,r,A,W,B,Pptr);

	long long int MY_PRIME = Pptr->prime;
	long long int INV_PRIME = Pptr->prime_inv;

	for(int i=0; i<n; i++) {
		A[i]=smallprimefield_Mul_INLINE(A[i], invn, INV_PRIME, MY_PRIME);
	}
}

