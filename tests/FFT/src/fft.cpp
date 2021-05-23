
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <bpas.h>
#include "oldfft.h"
#include "../../../include/FFT/src/arraybitreversal.h"
#include "../../../include/FFT/src/fft_furer1.h"
#include "../../../include/FFT/src/tft_tree1.h"
#include "../../../include/FFT/src/tft_tree2.h"
#include "../../../include/global.h"

//typedef __int64_t   sfixn;

using namespace std;

int d1 = 32;//;//32768;//1024;//4; //1; //2;
int d2 = 32;//16385;//32769;//1025;//4; //1; //2;
int v = 2;
int H = 1024;
int K=16*1024*1024;
int fft=1;

int basecase = 16;

void test0(sfixn p,int K){
	sfixn *Ap = (sfixn *)calloc(K, sizeof(sfixn));
	Ap = EX_RandomUniPolyCoeffsVec(K, p);
	  

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	  
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(K+1);
	sfixn n1 = 1L<<e;
	sfixn* rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	sfixn* dftAPtr = (sfixn *) my_calloc(n1,sizeof(sfixn));
	EX_Mont_GetNthRoots_OPT2_AS_GENE(e,n1,rootsPtr,pPtr);
	unsigned long long start;
	startTimer(&start);
 	oldMont_dft_OPT2_AS_GENE(n1,e,rootsPtr,dftAPtr,K,Ap,pPtr);
	float elapsed;
 	stopTimer(&start,&elapsed);
	cout << "Modpn library (sec): " <<"\t" << elapsed << endl;
	free(rootsPtr);
	free(dftAPtr);
	free(Ap);
	//----------------------------------------------------------------

}


void testf(sfixn p,int K, int check){
	printf("fft_furer1\n");
	p = 180143985094819841;
	sfixn *Ap = (sfixn *)calloc(K, sizeof(sfixn));
	Ap = EX_RandomUniPolyCoeffsVec(K, p);
	sfixn *Bp = (sfixn *)calloc(K, sizeof(sfixn));

	/*
	for(int i = 0; i < K; ++i)
	Ap[i] = i+1;*/



	for(int i = 0; i < K; ++i)
		Bp[i] = Ap[i];
	sfixn e = logceiling(K);
	printf("%ld,%ld\n",p,e);
	sfixn *SB1 = (sfixn *)my_calloc(K, sizeof(sfixn)); 
	int K2 = K<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	FURERPBPAS1::RootsTableFurer(K2, e+1, KRT);
	KRT += K2;
    
	
	unsigned long long start;
	startTimer(&start);
 	FURERPBPAS1::DFT_eff_p1(K, e, Ap, KRT, SB1);
	float elapsed;
 	stopTimer(&start,&elapsed);
	
	cout << "1-D FFTs (sec): " <<"\t" << elapsed << endl;

	if (check == 1) {
	  bool isEq = 1;
	  for (int i = 0; i < K; ++i) {
	    sfixn r = FURERPBPAS1::testDFT(K, i, Bp, KRT);
	    if (r != Ap[i]) {
	      isEq = 0;
	      printf("error in index %i!\n", i);
	      break;
	    }
	  }
	  if (isEq) { cout << "pass." << endl; }
	}

	my_free(SB1);
	my_free(KRT-K2);
	//free(invKRT);
	free(Ap);
}


void test1(sfixn p,int K, int check){
	sfixn *Ap = (sfixn *)calloc(K, sizeof(sfixn));
	Ap = EX_RandomUniPolyCoeffsVec(K, p);
	sfixn *Bp = (sfixn *)calloc(K, sizeof(sfixn));

	/*
	for(int i = 0; i < K; ++i)
	Ap[i] = i+1;*/



	for(int i = 0; i < K; ++i)
		Bp[i] = Ap[i];
	//for(int i =0;i<K;i++)
	//	Ap[i]= i;
	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(K);
	sfixn *SB1 = (sfixn *)my_calloc(K, sizeof(sfixn)); 
	int K2 = K<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);
	//FURERPBPAS1::RootsTableFurer(K2, e+1,KRT);
	KRT += K2;
    
	
	unsigned long long start;
	startTimer(&start);
 	PBPAS::DFT_eff(K, e, Ap, KRT, pPtr, H, NULL, SB1,1);
 //	FURERPBPAS1:: DFT_eff_p1(K, e,Ap,KRT,SB1);
	float elapsed;
 	stopTimer(&start,&elapsed);
	
	cout << "1-D FFTs (sec): " <<"\t" << elapsed << endl;
	//cout << (log2 (K)) <<"  " << log10 ((end - start) / 1000.f) << endl;
	//cout << (log2 (K)) <<"  " << (end - start) / 1000.f << endl;
	//cout << (log2 (K)) <<"  "<< (end - start) / 1000.f<<"  "<< log10 ((end - start) / 1000.f) << endl;



	/*
	TFT_tree1::Shuffle_tft(K,Ap,SB1);
	//permutation for tft order

	ofstream myfile;
	myfile.open ("fft_check_result.txt");
	for(int i = 0; i < K; ++i) {
	  myfile << Ap[i]<<"\n";
	}
	myfile.close();*/





	if (check == 1) {
	  bool isEq = 1;
	  for (int i = 0; i < K; ++i) {
	    sfixn r = PBPAS1::testDFT(K, i, Bp, KRT);
	   // sfixn r = FURERPBPAS1::testDFT(K, i, Bp, KRT);
	    
	    if (r != Ap[i]) {
	      isEq = 0;
	      printf("error in index %i!\n", i);
	      break;
	    }
	  }
	  if (isEq) { cout << "pass." << endl; }
	}

	
	// for(int i =0;i<K;i++){
	//   cout << "A "<<i<<"  is "<<Ap[i] << endl;

	//   }


	
	
	// sfixn R=1L<<(pPtr->Rpow); 
	// R=R%pPtr->P;
	// sfixn invn1=inverseMod(K,pPtr->P); 
	// invn1 = MulMod(R,invn1,pPtr->P);
	// invn1 = MulMod(R,invn1,pPtr->P);
	// invn1 <<= pPtr->Base_Rpow;


	
	// PBPAS::DFT_eff(K, e, Ap, KRT, pPtr, H, NULL, SB1,1);
	// PBPAS::InvDFT_eff(K, e, Ap, KRT2,pPtr,H,NULL, SB1,invn1,1);

	
	// for(int i =0;i<K;i++){
	//   cout << "Inverse FFT omega "<<i<<"  is "<<invKRT[i] << endl;

	// }
	
	my_free(SB1);
	my_free(KRT-K2);
//	free(invKRT);
	free(Ap);
}
/*
void Shufflefttest(int K, sfinx *A){
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);

	PBPAS1::Shuffle(K,A,B);
	PBPAS1::Shuffle(K>>1,A,B);
	PBPAS1::Shuffle(K>>1,A+(K>>1),B);
	free(B);
}
*/




void test11(int K){//test tft shuffle order
  sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
  sfixn* SB1 = (sfixn*) malloc(sizeof(sfixn)*K);

 
  for(int k = 0;k<K;k++)
    A[k] = k;

  for(int k = 0;k<K;k++)
    cout << k << "before shuffle is "<<A[k] << endl;



  TFT_tree1::Shuffle_tft(K,A,SB1);
  
  for(int k = 0;k<K;k++)
    cout << "B["<< k <<"]="<<A[k]<<";" << endl;
  
  /*
	ofstream myfile;
	myfile.open ("tft_shuffle_result.txt");
	for(int i = 0; i < K; ++i) {
	  myfile << A[i]<<"\n";
	}
	myfile.close();*/


	free(A);
	free(SB1);

}


void test2(int K){
	sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
	for(int k = 0;k<K;k++)
		A[k] = 2;
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);
	PBPAS1::Shuffle2(K,A,B);
}
void Shuffle3(int n, sfixn *A, sfixn *B){
	asm("nop");
	asm("nop");
	asm("nop");
	int n3 = n>>3;
	int n2 = n>>2;
	int n1 = n>>1;
	for (int i=0; i<n3; i+=1){
		int i2 = i<<3;
		A[i] = A[i2];
		B[i+3*n3] = A[i2+1];
		B[i+1*n3] = A[i2+2];
		B[i+5*n3] = A[i2+3];
		B[i+0*n3] = A[i2+4];
		B[i+4*n3] = A[i2+5];
		B[i+2*n3] = A[i2+6];
		B[i+6*n3] = A[i2+7];
	}
	//for (int i=0; i<n3; i+=128){
	//	for (int j=0;j<128;j++){
	//		int i2 = j<<3;
	//		int i3 = i<<3;
	//		B[i+j] = A[i3+i2+4];
	//	}
	//	for (int j=0;j<128;j++){
	//		int i2 = j<<3;
	//		int i3 = i<<3;
	//		B[i+j+3*n3] = A[i3+i2+1];
	//		B[i+j+1*n3] = A[i3+i2+2];
	//		B[i+j+5*n3] = A[i3+i2+3];
	//		B[i+j+4*n3] = A[i3+i2+5];
	//		B[i+j+2*n3] = A[i3+i2+6];
	//		B[i+j+6*n3] = A[i3+i2+7];
	//		A[i+j] = A[i3+i2];
	//	}
	//}
	sfixn *A2 = A + n3;
	memcpy(A2, B, 7*n3*(sizeof(sfixn)));   
	asm("nop");
}
void Shuffle2(int n, sfixn *A, sfixn *B){
	asm("nop");
	asm("nop");
	int n2 = n>>2;
	int n1 = n>>1;
	for (int i=0; i<n2; i+=1){
		int i2 = i<<2;
		A[i] = A[i2];
		B[i+n2] = A[i2+1];
		B[i] = A[i2+2];
		B[i+n1] = A[i2+3];
	}
	//for (int i=0; i<n2; i+=2048){
	//	for(int j=0;j<2048;j++){
	//		int i3 = i>>2;
	//		int i2 = j<<2;
	//		B[i3+j+n2] = A[i+i2+1];
	//		B[i3+j] = A[i+i2+2];
	//		B[i3+j+n1] = A[i+i2+3];
	//		A[i3+j] = A[i+i2];
	//	}
	//}
	sfixn *A2 = A + n2;
	memcpy(A2, B, 3*n2*(sizeof(sfixn)));   
	asm("nop");
}
void Shuffle1(int n, sfixn *A, sfixn *B){
	asm("nop");
	int n1 = n>>1;
	for (int i=0; i<n1; i+=1){
		int i2 = i<<1;
		A[i] = A[i2];
		B[i] = A[i2+1];
	}
	//for (int i=0; i<n1; i+=2048){
	//	for(int j=0;j<2048;j++){
	//		int i3 = i>>1;
	//		int i2 = j<<1;
	//		B[i3+j] = A[i2+i+1];
	//		A[i3+j] = A[i2+i];
	//	}
	//}
	sfixn *A2 = A + n1;
	memcpy(A2, B, n1*(sizeof(sfixn)));   
	asm("nop");
}
void test7(int K){
	sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
	for(int k = 0;k<K;k++)
		A[k] = 2;
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);
	free(A);
	free(B);
}
void test6(int K){
	sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
	for(int k = 0;k<K;k++)
		A[k] = 2;
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);

	for(int i=0;i<K;i++)
	  {
	   	A[i] = i+1;
	  }

	for(int i=0;i<K;i++)
	  {
	    cout << "before test6 shuffle3 :   A "<<i<<" is "<< A[i] << endl;
	  }

	Shuffle3(K,A,B);

	for(int i=0;i<K;i++)
	  {
	    cout << "result of test6 shuffle3 :   A "<<i<<" is "<< A[i] << endl;
	  }





	free(A);
	free(B);
}
void test5(int K){
	sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
	for(int k = 0;k<K;k++)
		A[k] = 2;
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);

	for(int i=0;i<K;i++)
	  {
	   	A[i] = i+1;
	  }

	for(int i=0;i<K;i++)
	  {
	    cout << "before test5 shuffle1+2 :   A "<<i<<" is "<< A[i] << endl;
	  }



	Shuffle1(K,A,B);
	PBPAS1::Shuffle2(K>>1,A,B);
	PBPAS1::Shuffle2(K>>1,A+(K>>1),B);
	//Shuffle1(K>>1,A+(K>>1),B);
	//Shuffle1(K>>2,A+2*(K>>2),B);
	//Shuffle1(K>>2,A+3*(K>>2),B);

	for(int i=0;i<K;i++)
	  {
	    cout << "result of test5 shuffle1+2 :   A "<<i<<" is "<< A[i] << endl;
	  }




	free(A);
	free(B);
}
void test4(int K){
	sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
	for(int k = 0;k<K;k++)
		A[k] = 2;
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);

	for(int i=0;i<K;i++)
	  {
	   	A[i] = i+1;
	  }

	for(int i=0;i<K;i++)
	  {
	    cout << "before test4 shuffle1 :   A "<<i<<" is "<< A[i] << endl;
	  }



	Shuffle1(K,A,B);
	Shuffle1(K>>1,A,B);
	Shuffle1(K>>1,A+(K>>1),B);
	Shuffle1(K>>2,A,B);
	Shuffle1(K>>2,A+(K>>2),B);
	Shuffle1(K>>2,A+2*(K>>2),B);
	Shuffle1(K>>2,A+3*(K>>2),B);
	//shuffle1 is regular permutaion

	for(int i=0;i<K;i++)
	  {
	    cout << "result of test4 shuffle1 :   A "<<i<<" is "<< A[i] << endl;
	  }


	free(A);
	free(B);
}

void test3(int K){
	sfixn* A = (sfixn*) malloc(sizeof(sfixn)*K);
	for(int k = 0;k<K;k++)
		A[k] = 2;

	for(int i=0;i<K;i++)
	  {
	   	A[i] = i+1;
	  }

	for(int i=0;i<K;i++)
	  {
	    cout << "before test3 shuffle :   A "<<i<<" is "<< A[i] << endl;
	  }
	
	sfixn* B = (sfixn*) malloc(sizeof(sfixn)*K);
	PBPAS1::Shuffle(K,A,B);
	PBPAS1::Shuffle(K>>1,A,B);
	PBPAS1::Shuffle(K>>1,A+(K>>1),B);


	for(int i=0;i<K;i++)
	  {
	    cout << "result of test3 shuffle :   A "<<i<<" is "<< A[i] << endl;
	  }


	free(A);
	free(B);
}

void test8(sfixn p,int K,int basecasefortft){
  int oldn;
  oldn = pow(2, ceil(log2 (K)));

	sfixn *Ap = (sfixn *)calloc(oldn, sizeof(sfixn));
	//Ap = EX_RandomUniPolyCoeffsVec(oldn, p);

	
	for(int i =0;i<K;i++)
		Ap[i]= i+1;
	for(int i =K;i<oldn;i++)
	  Ap[i]= 0;

	
	int s,r,l,m,relax_size,l_s,m_s,l_r,m_r;
        l = K;
        m = K;
     

  if (oldn > pow(basecasefortft,2)){
        s =  basecasefortft;
        r =  oldn/basecasefortft;
  }
  else{

    s = pow(2, ceil((log2 (oldn))/2.0));
    r = pow(2, floor( (log2 (oldn))/2.0));
  }
 
  
  l_s = ceil(l/(r*1.0));
  m_s = ceil(m/(r*1.0));
  l_r = std::min(l,r);
  m_r = std::min(m,r);
  

  relax_size = l_s * l_r;
  sfixn *invectmp = (sfixn *)calloc(oldn, sizeof(sfixn));
  for(int i = 0; i < oldn; ++i)
    {invectmp[i] = 0;}



	//cout <<" oldn is "<<oldn<<" k is "<<K << endl;
	sfixn *Dp = (sfixn *)calloc(oldn, sizeof(sfixn));
	//sfixn *Dp = (sfixn *)calloc(relax_size, sizeof(sfixn));
	for(int i = 0; i < K; ++i){
	  Dp[i] = Ap[i];
	}
	for(int i= K; i<relax_size; ++i){
	  Dp[i]=0;
	}
	
	//Dp[0]=9;Dp[1]=3;Dp[2]=13;Dp[3]=4;Dp[4]=3;Dp[5]=0;Dp[6]=0;Dp[7]=0;

	//Dp[0]=1;Dp[1]=2;Dp[2]=3;Dp[3]=4;Dp[4]=5;Dp[5]=0;Dp[6]=0;Dp[7]=0;

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(oldn);
	sfixn *SB1 = (sfixn *)my_calloc(oldn, sizeof(sfixn)); 
	int K2 = oldn<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);

	sfixn* invKRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::InverseRootsTable(K2,KRT,invKRT);
	KRT += K2;//for omega
	invKRT += K2;// for omega^-1


	//__cilkrts_set_param("nworkers", "4");
	__cilkrts_set_param("nworkers", "12");
	cout<<"Using "<<__cilkrts_get_nworkers()<<" workers"<<endl;

	unsigned long long start;
	startTimer(&start);
	TFT_tree1::TFT_Core_p1(oldn,K,l,m,basecasefortft,relax_size,KRT,e,SB1,Dp,invectmp);
 	float elapsed;
 	stopTimer(&start,&elapsed);
	
	//PBPAS::DFT_eff(oldn, e, Dp, KRT, pPtr, H, NULL, SB1,1);
	//cout << "Forward TFTs (sec): " <<"\t" << (end - start) / 1000.f << endl;
	//cout << (end - start) / 1000.f << endl;
	cout << (log2 (K)) <<"  "<< elapsed <<"  "<< log10 (elapsed) << endl;



	my_free(SB1);
	my_free(KRT-K2);
	free(Ap);
	//free(Bp);	
	//free(Cp);	
	free(Dp);

	free(invectmp);
	//free(invKRT);

}

void test9(sfixn p,int K,int basecasefortft){
  int oldn;
  oldn = pow(2, ceil(log2 (K)));

	sfixn *Ap = (sfixn *)calloc(oldn, sizeof(sfixn));
	//Ap = EX_RandomUniPolyCoeffsVec(oldn, p);

	
	for(int i =0;i<K;i++)
		Ap[i]= i+1;
	for(int i =K;i<oldn;i++)
	  Ap[i]= 0;

	
	int s,r,l,m,relax_size,l_s,m_s,l_r,m_r;
        l = K;
        m = K;
     

  if (oldn > pow(basecasefortft,2)){
        s =  basecasefortft;
        r =  oldn/basecasefortft;
  }
  else{

    s = pow(2, ceil((log2 (oldn))/2.0));
    r = pow(2, floor( (log2 (oldn))/2.0));
  }
 
  
  l_s = ceil(l/(r*1.0));
  m_s = ceil(m/(r*1.0));
  l_r = std::min(l,r);
  m_r = std::min(m,r);
  

  relax_size = l_s * l_r;
  sfixn *invectmp = (sfixn *)calloc(oldn, sizeof(sfixn));
  for(int i = 0; i < oldn; ++i)
    {invectmp[i] = 0;}


  
  // cout << "relaxe size is  " <<  relax_size  <<" oldn is " <<oldn<< endl;




	///verify part///
	/*
	sfixn *Dp = (sfixn *)calloc(oldn, sizeof(sfixn));
	if(K<oldn){
	  for(int i = 0; i < K; ++i)
	    {Dp[i] = i+1;}
	  for(int i = K; i < oldn; ++i)
	    {Dp[i] = 0;}
	}
	else{
	  for(int i = 0; i < K; ++i)
	    {Dp[i] = Ap[i];}
	}
	*/
	//cout <<" oldn is "<<oldn<<" k is "<<K << endl;
	sfixn *Dp = (sfixn *)calloc(oldn, sizeof(sfixn));
	//sfixn *Dp = (sfixn *)calloc(relax_size, sizeof(sfixn));
	for(int i = 0; i < K; ++i){
	  Dp[i] = Ap[i];
	}
	for(int i= K; i<relax_size; ++i){
	  Dp[i]=0;
	}
	
	//Dp[0]=9;Dp[1]=3;Dp[2]=13;Dp[3]=4;Dp[4]=3;Dp[5]=0;Dp[6]=0;Dp[7]=0;

	//Dp[0]=1;Dp[1]=2;Dp[2]=3;Dp[3]=4;Dp[4]=5;Dp[5]=0;Dp[6]=0;Dp[7]=0;

	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(oldn);
	sfixn *SB1 = (sfixn *)my_calloc(oldn, sizeof(sfixn)); 
	int K2 = oldn<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);

	sfixn* invKRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::InverseRootsTable(K2,KRT,invKRT);
	KRT += K2;//for omega
	invKRT += K2;// for omega^-1



	unsigned long long start;
	startTimer(&start);
 	TFT_tree1::TFT_Core_p1(oldn,K,l,m,basecasefortft,relax_size,KRT,e,SB1,Dp,invectmp);
	float elapsed;
 	stopTimer(&start,&elapsed);
	
	
	//PBPAS::DFT_eff(oldn, e, Dp, KRT, pPtr, H, NULL, SB1,1);
	//cout << "Forward TFTs (sec): " <<"\t" << (end - start) / 1000.f << endl;
	cout << elapsed << endl;
	
	/*
	for (int i = 0; i < oldn; ++i) {
	  cout<<"KRT "<< i<<" is " << KRT[i] << endl;
	}
	return;*/



	//Cp[0]=1;Cp[1]=2;Cp[2]=3;Cp[3]=4;Cp[4]=5;Cp[5]=0;Cp[6]=0;Cp[7]=0;
	///verify part///
	sfixn *testvec = (sfixn *)calloc(oldn, sizeof(sfixn));
	/*
	for (int i = 0; i < oldn; ++i) {
	  testvec[i] = TFT_tree1::testDFT(oldn, i, Cp, KRT);
	  }*/

	
	for (int i = 0; i < K; ++i) {
	  testvec[i] = Ap[i];
	}
	for (int i = K; i < oldn; ++i) {
	  testvec[i] = 0;
	}
	PBPAS::DFT_eff(oldn, e, testvec, KRT, pPtr, H, NULL, SB1,1);

	TFT_tree1::Shuffle_tft(oldn,testvec,SB1);
	//permutation for tft order

		


	bool isEqtft = 1;
	for (int i = 0; i < K; ++i) {
	  if (Dp[i] != testvec[i]) {
			isEqtft = 0;
			printf("tft error in index %i!\n", i);
			break;
		}
	}
	if (isEqtft) { cout << "tft pass." << endl; }

	/*
	ofstream myfile;
	myfile.open ("tft_parallel_result.txt");
	for(int i = 0; i < K; ++i) {
	  myfile << Dp[i]<<"\n";
	}
	myfile.close();*/



	my_free(SB1);
	my_free(KRT-K2);
	free(Ap);
	//free(Bp);	
	//free(Cp);	
	free(Dp);
	free(testvec);
	free(invectmp);
	//free(invKRT);

}


void test_itft_datacollect(sfixn p,int K,int basecasefortft){
  int oldn;
  oldn = pow(2, ceil(log2 (K)));
  

  int s,r,l,m,relax_size,l_s,m_s,l_r,m_r;
  l = K;
  m = K;
  
  
  if (oldn > pow(basecasefortft,2)){
    s =  basecasefortft;
    r =  oldn/basecasefortft;
  }
  else{
    
    s = pow(2, ceil((log2 (oldn))/2.0));
    r = pow(2, floor( (log2 (oldn))/2.0));
  }
 
  
  l_s = ceil(l/(r*1.0));
  m_s = ceil(m/(r*1.0));
  l_r = std::min(l,r);
  m_r = std::min(m,r);
  

  relax_size = l_s * l_r;
  sfixn *invectmp = (sfixn *)calloc(oldn, sizeof(sfixn));
  for(int i = 0; i < oldn; ++i)
    {invectmp[i] = 0;}
  
  
  sfixn *Dp = (sfixn *)calloc(oldn, sizeof(sfixn));
  sfixn *Ap = (sfixn *)calloc(oldn, sizeof(sfixn));//for test
  for(int i = 0; i < K; ++i){
    Dp[i]=i+1;  
    Ap[i]=i+1;
  }
  for(int i = K; i < oldn; ++i){
    Dp[i]=0;
    Ap[i]=0;
  }




	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(oldn);
	sfixn *SB1 = (sfixn *)my_calloc(oldn, sizeof(sfixn)); 
	int K2 = oldn<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);

	sfixn* invKRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::InverseRootsTable(K2,KRT,invKRT);
	KRT += K2;//for omega
	invKRT += K2;// for omega^-1

	//PBPAS::DFT_eff(oldn, e, Dp, KRT, pPtr, H, NULL, SB1,1);	
	//TFT_tree1::Shuffle_tft(oldn,Dp,SB1);
	//permutation for tft order


	//TFT_tree1::TFT_Core_p1(oldn,K,l,m,basecasefortft,relax_size,KRT,e,SB1,Dp,invectmp);

	//for(int i = K; i < oldn; ++i){
	//Dp[i]=0;
	//}




	//pre calculate invn for inverseTFT
	sfixn R=1L<<(pPtr->Rpow); 
	R=R%pPtr->P;	
	int invnstep = (log2 (oldn)) + 1;
	sfixn* invn = (sfixn *)my_calloc(invnstep, sizeof(sfixn));
	int original_n = oldn;
	for(int i=0;i< invnstep;i++){
	  invn[i] =inverseMod(original_n,pPtr->P);
	  invn[i] = MulMod(R,invn[i],pPtr->P);
	  invn[i] <<= pPtr->Base_Rpow;
	  original_n = original_n / 2;
	}

	//cout <<" invn  is  " << invn << endl;
	//TFT_tree1::Shuffle_tft(oldn,Dp,SB1);//shuffle the output of TFT
	//TFT_tree1::InvDFT_eff_p1( oldn, log2 (oldn), Dp, invKRT, Dp, invn);

	/*
	///////////////////////////// for collecting running time // same input with tft
	for(int i = 0; i < K; ++i){
	  Dp[i]=i+1;  
	}
	for(int i = K; i < oldn; ++i){
	  Dp[i]=0;
	  }*/
	/////////////////////////////////


	/*
	for(int i = 0; i < K; ++i){
	  Dp[i]=i+1;  
	  }*/

  sfixn mont_2 = TFT_tree1::normaltomont(2,mont_2);
  sfixn mont_inv2 = TFT_tree1::normaltomont(inverseMod(2,p),mont_inv2);
  sfixn primenum = 4179340454199820289;


  //__cilkrts_set_param("nworkers", "4");
  __cilkrts_set_param("nworkers", "12");
	cout<<"Using "<<__cilkrts_get_nworkers()<<" workers"<<endl;

	unsigned long long start;
	startTimer(&start);
 	//TFT_tree1::ITFT_Core_p1(oldn,oldn,1,1,K,K,0,KRT,invKRT,Dp,SB1);
	//TFT_tree1::ITFT_Wrapper_p1(oldn,oldn,1,1,K,K,0,KRT,invKRT,Dp,SB1,basecasefortft,invectmp,relax_size);
	TFT_tree1::Mont_ITFT_Core_p1(oldn,oldn,0,K,K,0,KRT,invKRT,Dp,mont_inv2,SB1,primenum);
	float elapsed;
 	stopTimer(&start,&elapsed);
	
	//cout << "Inverse TFTs (sec): " <<"\t" << (end - start) / 1000.f << endl;	
	//cout << K <<"&" << (end - start) / 1000.f << endl;
	cout << (log2 (K)) <<"  "<< elapsed <<"  "<< log10 (elapsed) << endl;






	//for(int i = 0; i < K; ++i   )
	//cout <<" inverse tft "<<i<<"  is  " <<Dp[i] << endl;


	free(Dp);
	//free(SB1);
	//free(KRT);


}






void test_itft(sfixn p,int K,int basecasefortft){
  int oldn;
  oldn = pow(2, ceil(log2 (K)));
  

  int s,r,l,m,relax_size,l_s,m_s,l_r,m_r;
  l = K;
  m = K;
  
  
  if (oldn > pow(basecasefortft,2)){
    s =  basecasefortft;
    r =  oldn/basecasefortft;
  }
  else{
    
    s = pow(2, ceil((log2 (oldn))/2.0));
    r = pow(2, floor( (log2 (oldn))/2.0));
  }
 
  
  l_s = ceil(l/(r*1.0));
  m_s = ceil(m/(r*1.0));
  l_r = std::min(l,r);
  m_r = std::min(m,r);
  

  relax_size = l_s * l_r;
  sfixn *invectmp = (sfixn *)calloc(oldn, sizeof(sfixn));
  for(int i = 0; i < oldn; ++i)
    {invectmp[i] = 0;}
  
  
  sfixn *Dp = (sfixn *)calloc(oldn, sizeof(sfixn));
  sfixn *Ap = (sfixn *)calloc(oldn, sizeof(sfixn));//for test
  for(int i = 0; i < K; ++i){
    Dp[i]=i+1;  
    Ap[i]=i+1;
  }
  for(int i = K; i < oldn; ++i){
    Dp[i]=0;
    Ap[i]=0;
  }




	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(oldn);
	sfixn *SB1 = (sfixn *)my_calloc(oldn, sizeof(sfixn)); 
	int K2 = oldn<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);

	sfixn* invKRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::InverseRootsTable(K2,KRT,invKRT);
	KRT += K2;//for omega
	invKRT += K2;// for omega^-1

	PBPAS::DFT_eff(oldn, e, Dp, KRT, pPtr, H, NULL, SB1,1);	
	TFT_tree1::Shuffle_tft(oldn,Dp,SB1);
	//permutation for tft order


	//TFT_tree1::TFT_Core_p1(oldn,K,l,m,basecasefortft,relax_size,KRT,e,SB1,Dp,invectmp);

	for(int i = K; i < oldn; ++i){
	  Dp[i]=0;
	}
	//for(int i = 0; i < K; ++i){
	//Dp[i]=i+1;
	//	}
	//Dp[0]=15,Dp[1]=1,Dp[2]=3,Dp[3]=12,Dp[4]=7,Dp[5]=0,Dp[6]=0,Dp[7]=0;
	
	//Dp[0]=1,Dp[1]=2,Dp[2]=3,Dp[3]=4,Dp[4]=5,Dp[5]=6,Dp[6]=7,Dp[7]=8,Dp[8]=9;
	//pre calculate invn for inverseTFT
	sfixn R=1L<<(pPtr->Rpow); 
	R=R%pPtr->P;	
	int invnstep = (log2 (oldn)) + 1;
	sfixn* invn = (sfixn *)my_calloc(invnstep, sizeof(sfixn));
	int original_n = oldn;
	for(int i=0;i< invnstep;i++){
	  invn[i] =inverseMod(original_n,pPtr->P);
	  invn[i] = MulMod(R,invn[i],pPtr->P);
	  invn[i] <<= pPtr->Base_Rpow;
	  original_n = original_n / 2;
	}

	//cout <<" invn  is  " << invn << endl;
	//TFT_tree1::Shuffle_tft(oldn,Dp,SB1);//shuffle the output of TFT
	//TFT_tree1::InvDFT_eff_p1( oldn, log2 (oldn), Dp, invKRT, Dp, invn);

	/*
	///////////////////////////// for collecting running time // same input with tft
	for(int i = 0; i < K; ++i){
	  Dp[i]=i+1;  
	}
	for(int i = K; i < oldn; ++i){
	  Dp[i]=0;
	  }*/
	/////////////////////////////////


	/*
	for(int i = 0; i < K; ++i){
	  Dp[i]=i+1;  
	  }*/


	  sfixn mont_2 = TFT_tree1::normaltomont(2,mont_2);
	  sfixn mont_inv2 = TFT_tree1::normaltomont(inverseMod(2,p),mont_inv2);
	  sfixn primenum = 4179340454199820289;


	unsigned long long start;
	startTimer(&start);
 	//TFT_tree1::ITFT_Core_p1(oldn,oldn,1,1,K,K,0,KRT,invKRT,Dp,SB1);
	TFT_tree1::ITFT_Wrapper_p1(oldn,oldn,1,1,K,K,0,KRT,invKRT,Dp,SB1,basecasefortft,invectmp,relax_size);
	//TFT_tree1::Mont_ITFT_Core_p1(oldn,oldn,0,K,K,0,KRT,invKRT,Dp,mont_inv2,SB1,primenum);
	float elapsed;
 	stopTimer(&start,&elapsed);
	
	cout << "Inverse TFTs (sec): " <<"\t" << elapsed << endl;	
	//cout << K <<"&" << (end - start) / 1000.f << endl;
	//for(int i = 0; i < K; ++i   )
	//cout <<" inverse tft "<<i<<"  is  " <<Dp[i] << endl;



	


	for (int i = 0; i < K; ++i){
	  Dp[i] = DivMod(Dp[i],oldn,4179340454199820289);
	  //std::cout <<" after tft invec "<<i<<"  is  " << invectmp[i] << std::endl;
	}

	ofstream myfile;
	myfile.open ("inverse_tft_serial_result.txt");
	for(int i = 0; i < K; ++i) {
	  myfile << Dp[i]<<"\n";
	}
	myfile.close();

	/*
	for (int i = 0; i < K; ++i){
	  std::cout <<" after itft invec "<<i<<"  is  " << Dp[i] << std::endl;
	  }*/

	
	bool isEqtft = 1;
	for (int i = 0; i < K; ++i) {
	  
	  if (Dp[i] != Ap[i]) {
			isEqtft = 0;
			printf("inverse tft error in index %i!\n", i);
			break;
		}
	}
	if (isEqtft) { cout << "inverse tft pass." << endl; }
	
	
	//forward use DFT_eff
	/*
	for(int i = 0; i < K; ++i){
	  Dp[i]=i+1;  
	  Ap[i]=i+1;
	}
	for(int i = K; i < oldn; ++i){
	  Dp[i]=0;
	  Ap[i]=0;
	}

	PBPAS::DFT_eff(K, e, Dp, KRT, pPtr, H, NULL, SB1,1); //prepare for InvDFT_Eff
	sfixn R=1L<<(pPtr->Rpow); 
	R=R%pPtr->P;	
	sfixn invn=inverseMod(K,pPtr->P);
	invn = MulMod(R,invn,pPtr->P);
	invn <<= pPtr->Base_Rpow;
	  cout <<" invn  is  " << invn << endl;
	PBPAS::InvDFT_eff(K, e, Dp, invKRT, pPtr, H, NULL, Dp,invn,1);

	for(int i = 0; i < K; ++i   )
	  cout <<" inverse fft "<<i<<"  is  " <<Dp[i] << endl;
	*/

	free(Ap);
	free(Dp);
	//free(SB1);
	//free(KRT);


}


void test10(sfixn p,int K,int basecasefortft){
  int oldn;
  oldn = pow(2, ceil(log2 (K)));
  sfixn *invectmp = (sfixn *)calloc(oldn, sizeof(sfixn));
  
	sfixn *Ap = (sfixn *)calloc(oldn, sizeof(sfixn));
	for(int i =0;i<K;i++)
		Ap[i]= i+1;
	for(int i =K;i<oldn;i++)
		Ap[i]= 0;

	int s,r,l,m,relax_size,l_s,m_s,l_r,m_r;
        l = K;
        m = K;
     

  if (K > pow(basecasefortft,2)){
        s =  basecasefortft;
        r =  K/basecasefortft;
  }
  else{

    s = pow(2, ceil((log2 (K))/2.0));
    r = pow(2, floor( (log2 (K))/2.0));
  }
 
  
  l_s = ceil(l/(r*1.0));
  m_s = ceil(m/(r*1.0));
  l_r = std::min(l,r);
  m_r = std::min(m,r);
  

  relax_size = l_s * l_r;

  
  //cout << "relaxe size is  " <<  relax_size  <<" oldn is " <<oldn<< endl;





	sfixn *Bp = (sfixn *)calloc(K, sizeof(sfixn));
	for(int i = 0; i < K; ++i)
		Bp[i] = Ap[i];


	///verify part///
	sfixn *Cp = (sfixn *)calloc(oldn, sizeof(sfixn));
	if(K<oldn){
	  for(int i = 0; i < K; ++i)
	    {Cp[i] = Ap[i];}
	  for(int i = K; i < oldn; ++i)
	    {Cp[i] = 0;}
	}
	else{
	  for(int i = 0; i < K; ++i)
	    {Cp[i] = Ap[i];}
	}


	///verify part///
	sfixn *Dp = (sfixn *)calloc(oldn, sizeof(sfixn));
	if(K<oldn){
	  for(int i = 0; i < K; ++i)
	    {Dp[i] = Ap[i];}
	  for(int i = K; i < oldn; ++i)
	    {Dp[i] = 0;}
	}
	else{
	  for(int i = 0; i < K; ++i)
	    {Dp[i] = Ap[i];}
	}



	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	sfixn e = logceiling(oldn);
	sfixn *SB1 = (sfixn *)my_calloc(oldn, sizeof(sfixn)); 
	int K2 = oldn<<1;
	sfixn* KRT = (sfixn *)my_calloc(K2<<1, sizeof(sfixn));
	PBPAS::RootsTableSpe(K2, e+1, KRT, pPtr,1);
	KRT += K2;

	unsigned long long start;
	startTimer(&start);
 	//TFT_tree1::DFT_eff_p1(K, e, Ap, KRT, SB1);

	TFT_tree1::TFT_Core_p1(oldn,oldn,l,m,basecasefortft,relax_size,KRT,e,SB1,Dp,invectmp);
        //TFT_tree1::InitTFT_Tree_p1(  basecasefortft,Ap   );
	float elapsed;
 	stopTimer(&start,&elapsed);
	

	
	//cout <<(log2 (K)) <<" " << log10  ((end - start) / 1000.f) << endl;
	cout <<(log2 (K)) <<" " << elapsed << endl;	





	my_free(SB1);
	my_free(KRT-K2);
	free(Ap);
	free(Bp);	
	free(Cp);	
	free(Dp);
	free(invectmp);
	//free(testvec);
}


bool isInteger(double a){
    double b=round(a),epsilon=1e-9; //some small range of error
    return (a<=b+epsilon && a>=b-epsilon); 
}
//test mul two 64-bit int mod a 64-bit prime
int main(int argc, char *argv[])
{//------------------------------------------
  //cout << "#workers:\t" << __cilkrts_get_nworkers() << endl;
	if(argc > 2) K = atoi(argv[2]);
	if(argc > 1) fft = atoi(argv[1]);
	if(argc > 3) basecase = atoi(argv[3]);//basecase for tft

	if(basecase<=2 || basecase==K || !isInteger(log2 (basecase)) ){
	  cout << "Please choose correct basecase." << endl;
	  return 0;
	}


	sfixn p = 4179340454199820289;//3799912185593857; //
	//test3(d1,d2,K,v,p);
	//test4(d1,d2,K,v,p);
	if(fft==1)
	  test1(p,K,1);
	else if(fft==2)
		test2(K);
	else if(fft==3)
		test3(K);
	else if(fft==4)
		test4(K);
	else if(fft==5)
		test5(K);
	else if(fft==6)
		test6(K);
	else if(fft==7)
		test7(K);	
	else if(fft==8)
	        test8(p,K,basecase);
	else if(fft==9)
	        test9(p,K,basecase);
	else if(fft==10)
	        test10(p,K,basecase);	
	else if(fft==11)
	        test11(K);	
	else if(fft==14)
	        test_itft(p,K,basecase);
	else if(fft==15)
                test_itft_datacollect(p,K,basecase);
        else if(fft==16)
	         testf(p,K,1);
	else
		test0(p,K);

  return 0;
}


/*
for size=4, omega is  3360066027580426122;
for size=8, omega is  3324705732702508476; omega^-1 is 1324460247237108937;
for size=16, omega is 1807568976454537418;
for size=32, omega is 1618839207552315154;
for size=512, omega is 440865335274882933;
for size=1024, omega is 2755315832028327088;
for size=2048, omega is 3149061058162746603;
 */
/*
InitTFT_Tree(4179340454199820289, 1618839207552315154, 8, Vector([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]))
*/



