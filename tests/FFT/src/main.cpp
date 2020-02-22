#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>
#include <cilktools/cilkview.h>

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <bpas.h>
#include "oldfft.h"
//typedef __int64_t   sfixn;

using namespace std;

int d1 = 32;//;//32768;//1024;//4; //1; //2;
int d2 = 32;//16385;//32769;//1025;//4; //1; //2;
int v = 2;
int K = 16*1024*1024;//32768; //65536;//2048;//16;//2; //4;
int H = 1024;
void test6(sfixn d1,sfixn d2,sfixn K,sfixn v,sfixn p){
	printf("%lu\n",p);
	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	printf("%lu\n",pPtr->Rpow);
	sfixn e = logceiling(K+1);
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1
	int es1 = PBPAS::ceil_log2_long(H);
	e = es1+1;
	int K2 = H<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTable(K2, e, KRT, pPtr);
	FILE *fichier;
	fichier = fopen("roots","w");
	for(int i=0;i<2*H-2;i++){
		fprintf(fichier,"%lu\n",*(KRT+K2+i));	
	}
	fclose(fichier);
	free(KRT);

}
void test7(sfixn d1,sfixn d2,sfixn K,sfixn v,sfixn p){
	cout<<"The prime: " <<p<<endl;
	sfixn *Ap = (sfixn *)calloc(K, sizeof(sfixn));
	Ap = EX_RandomUniPolyCoeffsVec(K, p);
	//for(int i =0;i<K;i++)
	//	Ap[i]= i;
	int * RevBidMap = (int *)my_malloc(H*sizeof(sfixn));  
	for(sfixn i=0; i<H; ++i)   
	    RevBidMap[i] = i;
	PBPAS::RevBitInd(H, RevBidMap);

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	printf("%lu\n",pPtr->Rpow);
	sfixn e = logceiling(K);
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1
	sfixn *SB1 = (sfixn *)my_calloc(K, sizeof(sfixn)); 
	int K2 = K<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);
	KRT += K2;
	unsigned long start = __cilkview_getticks();
	PBPAS::DFT_eff(K, e, Ap, KRT, pPtr, H, RevBidMap, SB1,1);
	unsigned long end = __cilkview_getticks();
	cout << "new time (ms): " <<"\t" << (end - start) /1.f << endl;
	//---------------------------------------------------------------

	my_free(SB1);
	my_free(KRT-K2);
	free(Ap);
}

void test4(sfixn d1,sfixn d2,sfixn K,sfixn v,sfixn p){
	cout<<"The prime: " <<p<<endl;
	sfixn *Ap = (sfixn *)calloc(K, sizeof(sfixn));
	Ap = EX_RandomUniPolyCoeffsVec(K, p);
	//for(int i =0;i<K;i++)
	//	Ap[i]= i;
	int * RevBidMap = (int *)my_malloc(H*sizeof(sfixn));  
	for(sfixn i=0; i<H; ++i)   
	    RevBidMap[i] = i;
	PBPAS::RevBitInd(H, RevBidMap);

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	printf("%lu\n",pPtr->Rpow);
	sfixn e = logceiling(K+1);
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1
	sfixn *SB1 = (sfixn *)my_calloc(K, sizeof(sfixn)); 
	int es1 = PBPAS::ceil_log2_long(K);
	e = es1+1;
	int K2 = K<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTable(K2, e, KRT, pPtr);
	KRT += K2;
	unsigned long start = __cilkview_getticks();
	PBPAS::DFT_eff(K, es1, Ap, KRT, pPtr, H, RevBidMap, SB1,1);
	unsigned long end = __cilkview_getticks();
	cout << "new time (ms): " <<"\t" << (end - start) /1.f << endl;
	//---------------------------------------------------------------
	int es2 = logceiling(ls2);
	int dims2 = 1<<es2; 
	//dimension (power of 2) of y //y-size of the 2-D FFT of the result

	my_free(SB1);
	my_free(KRT-K2);
	free(Ap);
}
void test3(sfixn d1,sfixn d2,sfixn L,sfixn v,sfixn p){
	  
	  //sfixn p = 10645405947658241; //54 bits, 2^27 * 79314455 + 1, correct cyclic
	  //8529266726959324, 8529266726959324, 8529266726959324, 8529266726959324, 8, 8, 8, 8, 2116143784101677, 2116143784101677, 2116143784101677, 2116143784101677

	  //sfixn p = 10645410511060993; //54 bits, 2^27 * 79314489 + 1
	  //sfixn p = 24686718519083009; //55 bits, 2^28 * 91965193 + 1
	  //sfixn p = 24686728719630337; //55 bits, 2^28 * 91965231 + 1
	  //sfixn p = 38349804861915137; //56 bits, 2^30 * 35716039 +1

	cout<<"The prime: " <<p<<endl;
	sfixn *Ap = (sfixn *)calloc(K, sizeof(sfixn));
	for(sfixn j=0; j<K; ++j)
		Ap[j] = j;
	  
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(K+1);
	sfixn n1 = 1L<<e;
	sfixn* rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	sfixn* dftAPtr = (sfixn *) my_calloc(n1,sizeof(sfixn));
	EX_Mont_GetNthRoots_OPT2_AS_GENE(e,n1,rootsPtr,pPtr);
	unsigned long long start = __cilkview_getticks();
	oldMont_dft_OPT2_AS_GENE(n1,e,rootsPtr,dftAPtr,K,Ap,pPtr);
	unsigned long long end = __cilkview_getticks();
	cout << "old time (ms): " <<"\t" << (end - start) /1.f << endl;
	free(rootsPtr);
	free(dftAPtr);
	//----------------------------------------------------------------
	free(Ap);

}
void test2(sfixn d1,sfixn d2,sfixn L,sfixn v,sfixn p){

	cout<<"The prime: " <<p<<endl;
	sfixn *Ap = (sfixn *)calloc(K*d1, sizeof(sfixn));
	for(sfixn i=0;i<d1;i++)
		for(sfixn j=0; j<K; ++j)
			Ap[j*d1+i] = v;
	  
	sfixn *Bp = (sfixn *)calloc(K*d2, sizeof(sfixn));
	  
	for(sfixn i=0;i<d2;i++)
		for(sfixn j=0; j<K; ++j)
			Bp[j*d2+i] = v;
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(K+1);
	sfixn n1 = 1L<<e;
	sfixn* rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	sfixn* dftAPtr = (sfixn *) my_calloc(n1,sizeof(sfixn));
	EX_Mont_GetNthRoots_OPT2_AS_GENE(e,n1,rootsPtr,pPtr);
	unsigned long long start = __cilkview_getticks();
	//Mont_dft_OPT2_AS_GENE(n1,e,rootsPtr,dftAPtr,K,Ap,pPtr);
	unsigned long long end = __cilkview_getticks();
	cout << "old time (ms): " <<"\t" << (end - start) /1.f << endl;
	free(rootsPtr);
	free(dftAPtr);
	//----------------------------------------------------------------
	int * RevBidMap = (int *)my_malloc(H*sizeof(sfixn));  
	for(sfixn i=0; i<H; ++i)   
	    RevBidMap[i] = i;
	PBPAS::RevBitInd(H, RevBidMap);

	int es1 = PBPAS::ceil_log2_long(K);
	e = es1+1;
	start = __cilkview_getticks();
	sfixn * c2 = PBPAS::TwoConvolutionModNew(Ap,Bp,d1,d2,K, pPtr,H,RevBidMap,1);
	end = __cilkview_getticks();
	cout << "new time (ms): " <<"\t" << (end - start) /1.f << endl;
	sfixn* c22 = c2;
	cout<<c22[0+K*(d1+d2)]<<endl;
	cout<<c22[0]<<endl;
	cout<<c22[1]<<endl;
	cout<<c22[2]<<endl;
	//---------------------------------------------------------------
	int es2 = logceiling(ls2);
	int dims2 = 1<<es2; 
	//dimension (power of 2) of y //y-size of the 2-D FFT of the result

	my_free(c2);

}
void test5(sfixn d1,sfixn d2,sfixn L,sfixn v,sfixn p){
	  
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1
	cout<<"The prime: " <<p<<endl;
	sfixn *Ap= (sfixn *)malloc(K*d1*sizeof(sfixn));//= EX_RandomUniPolyCoeffsVec(K*d1, p);//
	for(sfixn i=0;i<d1;i++)
		for(sfixn j=0; j<K; ++j)
			Ap[j*d1+i] = v;
	  
	sfixn *Bp = (sfixn *)malloc(K*d2*sizeof(sfixn));
	  
	for(sfixn i=0;i<d2;i++)
		for(sfixn j=0; j<K; ++j)
			Bp[j*d2+i] = v;
	sfixn s = K*ls2; //size of result // real size of the result from convolution
	printf("%lu %lu\n",s,(s<<1)*sizeof(sfixn));
	//sfixn *ncc = (sfixn *)my_calloc(s<<1, sizeof(sfixn));       
	sfixn *ncc = (sfixn*)calloc((s<<1),sizeof(sfixn));
	if(ncc == NULL)
		printf("error\n");
	*ncc = 0;
	free(ncc);
	free(Ap);
	free(Bp);
	return;
}
void test1(sfixn d1,sfixn d2,sfixn L,sfixn v,sfixn p){
	  
	  //sfixn p = 10645405947658241; //54 bits, 2^27 * 79314455 + 1, correct cyclic
	  //8529266726959324, 8529266726959324, 8529266726959324, 8529266726959324, 8, 8, 8, 8, 2116143784101677, 2116143784101677, 2116143784101677, 2116143784101677

	  //sfixn p = 10645410511060993; //54 bits, 2^27 * 79314489 + 1
	  //sfixn p = 24686718519083009; //55 bits, 2^28 * 91965193 + 1
	  //sfixn p = 24686728719630337; //55 bits, 2^28 * 91965231 + 1
	  //sfixn p = 38349804861915137; //56 bits, 2^30 * 35716039 +1

	sfixn *Ap;// = (sfixn *)calloc(K*d1, sizeof(sfixn));
	//for(sfixn i=0;i<d1;i++)
		Ap = EX_RandomUniPolyCoeffsVec(d1*K, p);
		//for(sfixn j=0; j<K; ++j)
		//	Ap[j*d1+i] = v;
	  
	sfixn *Bp ;//= (sfixn *)calloc(K*d2, sizeof(sfixn));
	  
	//for(sfixn i=0;i<d2;i++)
		Bp = EX_RandomUniPolyCoeffsVec(d2*K, p);
		//for(sfixn j=0; j<K; ++j)
		//	Bp[j*d2+i] = v;
	int ls2 = d1+d2-1; // y-size of A*B mod x^K + 1

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(K+1);
	sfixn n1 = 1L<<e;
	sfixn* rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	sfixn* dftAPtr = (sfixn *) my_calloc(n1,sizeof(sfixn));
	EX_Mont_GetNthRoots_OPT2_AS_GENE(e,n1,rootsPtr,pPtr);
	printf("%lu\n",pPtr->Max_Root);
	unsigned long long start = __cilkview_getticks();
	sfixn * c1 = PBPAS::TwoConvolutionMod(Ap,Bp,d1,d2,p,K);
	//Mont_dft_OPT2_AS_GENE(n1,e,rootsPtr,dftAPtr,K,Ap,pPtr);
	unsigned long long end = __cilkview_getticks();
	cout << "old time (ms): " <<"\t" << (end - start) /1.f << endl;
	free(rootsPtr);
	free(dftAPtr);
	sfixn * c12 = c1;
	//----------------------------------------------------------------
	int * RevBidMap = (int *)my_malloc(H*sizeof(sfixn));  
	for(sfixn i=0; i<H; ++i)   
	    RevBidMap[i] = i;
	PBPAS::RevBitInd(H, RevBidMap);

	start = __cilkview_getticks();
	sfixn * c2 = PBPAS::TwoConvolutionModNew(Ap,Bp,d1,d2,K, pPtr,H,RevBidMap,1);
	end = __cilkview_getticks();
	cout << "new time (ms): " <<"\t" << (end - start) /1.f << endl;
	//---------------------------------------------------------------
	//dimension (power of 2) of y //y-size of the 2-D FFT of the result

	sfixn * c22 = c2;
	int flag=1;
	if (flag==1) {
		for (sfixn i=0; i<K*ls2; i++){
			if (c12[i]!=c22[i]){
				std::cout<<"new and old does not match at i, c12, c22: "<<i<<", "<<c12[i]<<", "<<c22[i]<<"\n\n";
				flag=0;
				break;
			}
		}
	}

	cout<<"display"<<endl;
	cout<<c12[0+K*(d1+d2)]<<endl;
	cout<<c22[0+K*(d1+d2)]<<endl;
	cout<<c12[0]<<endl;
	cout<<c22[0]<<endl;
	cout<<c12[1]<<endl;
	cout<<c22[1]<<endl;
	cout<<c12[2]<<endl;
	cout<<c22[2]<<endl;
	if (flag==1) {
		std::cout<<"new and old matches \n\n";
	}

	 my_free(c1);
	my_free(c2);

}

//test mul two 64-bit int mod a 64-bit prime
int main(int argc, char *argv[])
{//------------------------------------------
	cout << "#workers:\t" << __cilkrts_get_nworkers() << endl;
	if(argc > 1) d1 = atoi(argv[1]);
	if(argc > 2) d2 = atoi(argv[2]);
	if(argc > 3) K = atoi(argv[3]);
	if(argc > 4) v = atoi(argv[4]);
	sfixn p = 4179340454199820289;//3799912185593857; //
	//test3(d1,d2,K,v,p);
	//test4(d1,d2,K,v,p);
	test7(d1,d2,K,v,p);

  return 0;
}


/* K=4, d1=d2=2, all 1, wrong
---Negacyclic reslut: 
10133099161583469, 38136497311074514, 103503911730324612, 516024752800238164, 7881299347898365, 0, 4, 8, 7881299347898367, 0, 2, 4, 

7881299347898367, 0, 2, 4, 7881299347898365, 0, 4, 8, 7881299347898367, 0, 2, 4

---Cyclic reslut: 
4, 4, 4, 4, 8, 8, 8, 8, 4, 4, 4, 4, 

before handling R:
3940649673949181, 3940649673949181, 3940649673949181, 3940649673949181, 7881299347898362, 7881299347898362, 7881299347898362, 7881299347898362, 3940649673949181, 3940649673949181, 3940649673949181, 3940649673949181

*/

/* half p-1
---Negacyclic reslut: 
1, 2, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0, 

---Cyclic reslut: 
1, 2, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0,
 */

/* all p-1
---Negacyclic reslut: 
7881299347898367, 0, 2, 4, 7881299347898365, 0, 4, 8, 7881299347898367, 0, 2, 4, 

---Cyclic reslut: 
4, 4, 4, 4, 8, 8, 8, 8, 4, 4, 4, 4
 */

/* all p+1, wrong
---Negacyclic reslut: 
10133099161583469, 38136497311074514, 103503911730324612, 516024752800238164, 7881299347898365, 0, 4, 8, 7881299347898367, 0, 2, 4, 

---Cyclic reslut: 
4, 4, 4, 4, 8, 8, 8, 8, 4, 4, 4, 4,

 //DnC FFT test
The prime: 3799912185593857
int d1 = 16; 
int d2 = 16; 
int v = 1;
int K = 16;

---Ap: 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
---Bp: 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
---Negacyclic out: 
3799912185593843, 3799912185593845, 3799912185593847, 3799912185593849, 3799912185593851, 3799912185593853, 3799912185593855, 0, 2, 4, 6, 8, 10, 12, 14, 16, 3799912185593829, 3799912185593833, 3799912185593837, 3799912185593841, 3799912185593845, 3799912185593849, 3799912185593853, 0, 4, 8, 12, 16, 20, 24, 28, 32, 3799912185593815, 3799912185593821, 3799912185593827, 3799912185593833, 3799912185593839, 3799912185593845, 3799912185593851, 0, 6, 12, 18, 24, 30, 36, 42, 48, 3799912185593801, 3799912185593809, 3799912185593817, 3799912185593825, 3799912185593833, 3799912185593841, 3799912185593849, 0, 8, 16, 24, 32, 40, 48, 56, 64, 3799912185593787, 3799912185593797, 3799912185593807, 3799912185593817, 3799912185593827, 3799912185593837, 3799912185593847, 0, 10, 20, 30, 40, 50, 60, 70, 80, 3799912185593773, 3799912185593785, 3799912185593797, 3799912185593809, 3799912185593821, 3799912185593833, 3799912185593845, 0, 12, 24, 36, 48, 60, 72, 84, 96, 3799912185593759, 3799912185593773, 3799912185593787, 3799912185593801, 3799912185593815, 3799912185593829, 3799912185593843, 0, 14, 28, 42, 56, 70, 84, 98, 112, 3799912185593745, 3799912185593761, 3799912185593777, 3799912185593793, 3799912185593809, 3799912185593825, 3799912185593841, 0, 16, 32, 48, 64, 80, 96, 112, 128, 3799912185593731, 3799912185593749, 3799912185593767, 3799912185593785, 3799912185593803, 3799912185593821, 3799912185593839, 0, 18, 36, 54, 72, 90, 108, 126, 144, 3799912185593717, 3799912185593737, 3799912185593757, 3799912185593777, 3799912185593797, 3799912185593817, 3799912185593837, 0, 20, 40, 60, 80, 100, 120, 140, 160, 3799912185593703, 3799912185593725, 3799912185593747, 3799912185593769, 3799912185593791, 3799912185593813, 3799912185593835, 0, 22, 44, 66, 88, 110, 132, 154, 176, 3799912185593689, 3799912185593713, 3799912185593737, 3799912185593761, 3799912185593785, 3799912185593809, 3799912185593833, 0, 24, 48, 72, 96, 120, 144, 168, 192, 3799912185593675, 3799912185593701, 3799912185593727, 3799912185593753, 3799912185593779, 3799912185593805, 3799912185593831, 0, 26, 52, 78, 104, 130, 156, 182, 208, 3799912185593661, 3799912185593689, 3799912185593717, 3799912185593745, 3799912185593773, 3799912185593801, 3799912185593829, 0, 28, 56, 84, 112, 140, 168, 196, 224, 3799912185593647, 3799912185593677, 3799912185593707, 3799912185593737, 3799912185593767, 3799912185593797, 3799912185593827, 0, 30, 60, 90, 120, 150, 180, 210, 240, 3799912185593633, 3799912185593665, 3799912185593697, 3799912185593729, 3799912185593761, 3799912185593793, 3799912185593825, 0, 32, 64, 96, 128, 160, 192, 224, 256, 3799912185593647, 3799912185593677, 3799912185593707, 3799912185593737, 3799912185593767, 3799912185593797, 3799912185593827, 0, 30, 60, 90, 120, 150, 180, 210, 240, 3799912185593661, 3799912185593689, 3799912185593717, 3799912185593745, 3799912185593773, 3799912185593801, 3799912185593829, 0, 28, 56, 84, 112, 140, 168, 196, 224, 3799912185593675, 3799912185593701, 3799912185593727, 3799912185593753, 3799912185593779, 3799912185593805, 3799912185593831, 0, 26, 52, 78, 104, 130, 156, 182, 208, 3799912185593689, 3799912185593713, 3799912185593737, 3799912185593761, 3799912185593785, 3799912185593809, 3799912185593833, 0, 24, 48, 72, 96, 120, 144, 168, 192, 3799912185593703, 3799912185593725, 3799912185593747, 3799912185593769, 3799912185593791, 3799912185593813, 3799912185593835, 0, 22, 44, 66, 88, 110, 132, 154, 176, 3799912185593717, 3799912185593737, 3799912185593757, 3799912185593777, 3799912185593797, 3799912185593817, 3799912185593837, 0, 20, 40, 60, 80, 100, 120, 140, 160, 3799912185593731, 3799912185593749, 3799912185593767, 3799912185593785, 3799912185593803, 3799912185593821, 3799912185593839, 0, 18, 36, 54, 72, 90, 108, 126, 144, 3799912185593745, 3799912185593761, 3799912185593777, 3799912185593793, 3799912185593809, 3799912185593825, 3799912185593841, 0, 16, 32, 48, 64, 80, 96, 112, 128, 3799912185593759, 3799912185593773, 3799912185593787, 3799912185593801, 3799912185593815, 3799912185593829, 3799912185593843, 0, 14, 28, 42, 56, 70, 84, 98, 112, 3799912185593773, 3799912185593785, 3799912185593797, 3799912185593809, 3799912185593821, 3799912185593833, 3799912185593845, 0, 12, 24, 36, 48, 60, 72, 84, 96, 3799912185593787, 3799912185593797, 3799912185593807, 3799912185593817, 3799912185593827, 3799912185593837, 3799912185593847, 0, 10, 20, 30, 40, 50, 60, 70, 80, 3799912185593801, 3799912185593809, 3799912185593817, 3799912185593825, 3799912185593833, 3799912185593841, 3799912185593849, 0, 8, 16, 24, 32, 40, 48, 56, 64, 3799912185593815, 3799912185593821, 3799912185593827, 3799912185593833, 3799912185593839, 3799912185593845, 3799912185593851, 0, 6, 12, 18, 24, 30, 36, 42, 48, 3799912185593829, 3799912185593833, 3799912185593837, 3799912185593841, 3799912185593845, 3799912185593849, 3799912185593853, 0, 4, 8, 12, 16, 20, 24, 28, 32, 3799912185593843, 3799912185593845, 3799912185593847, 3799912185593849, 3799912185593851, 3799912185593853, 3799912185593855, 0, 2, 4, 6, 8, 10, 12, 14, 16, 

---Cyclic out: 
16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 224, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 


*/
