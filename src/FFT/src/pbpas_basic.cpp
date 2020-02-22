// -*- C++ -*- 
// @author Yuzhen Xie

#if DEBUG
#include <cilk/cilk_api.h>
#include <cilktools/cilkview.h>
#endif

#include <cilk/cilk.h>

#include "../../../include/FFT/src/modpn.h"

#include "../../../include/FFT/src/basic_routine.h"

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef __64int_t sfixn;
//typedef __int32_t int32;
using namespace std;
namespace PBPAS {
    
    /*
      C = A*B mod p
      p: prime number
      N: number of variables
      C, A, B: dense recursive representation of C, A and B
      dc, da, db: partial degree vectors for C, A and B
      C should be initialized to zero
      default wp (whichprime) =0 to call the fft for a general Frouier prime
     */
    bool DMPMul(sfixn p,  int N,    
		sfixn *C, int *dc, 
		sfixn *A, int *da,
		sfixn *B, int *db,
		int wp=0)
    {//-------------------------------------------------------------
	MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*) my_malloc(sizeof(MONTP_OPT2_AS_GENE));
	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	//ToDo: N=1, large, contract to 2D
	if (N==1){
	    cout<<"ToDo: N=1, large N?, integrate contract to 2D method"<<endl;
	}

	if (N==2){
	    BivarMul(C,dc,A,da,B,db,pPtr,wp);
	    return 1;
	}

	//N>2

	my_free(pPtr);
	return 0;
    }

    /*
      C = A*B mod pPtr->P
      two variables
      C, A, B: dense recursive representation of C, A and B
      dc, da, db: partial degree vectors for C, A and B
      dc[0] = partial degree of variable 1
      pPtr: struct for Montgmery mod arithmetic
      C should be initialized to zero
     */
    void BivarMul(sfixn *C, int *dc, 
		  sfixn *A, int *da,
		  sfixn *B, int *db,
		  MONTP_OPT2_AS_GENE *pPtr,
		  int wp)
    {//---------------------------------------
	static const int H = 1024; //=DFTBASESIZE
	int ls1 = dc[0]+1;
	int es1 = logceiling(ls1);
	int dims1 = 1<<es1;
	int ls2 = dc[1]+1;
	int es2 = logceiling(ls2);
	int dims2 = 1<<es2;
	
	if ( (es2 > (pPtr->Npow)) || (es1 > (pPtr->Npow)) ){
	    std::cout<<"BPAS: exception, the prime number is too small for the FFT problem."<<"es1: "<<es1<<", es2: "<<es2<<", p: "<<pPtr->P<<std::endl;
	    exit(1); //checkFrourierPrime
	}
	//improve, just compute for the largest dimension
	int RT1_size = (dims1<<1)-2; //2*dims1-2	
	sfixn *RT1   = (sfixn *) my_calloc(RT1_size<<1, sizeof(sfixn));
	PBPAS::RootsTable(dims1, es1, RT1, pPtr);
	sfixn *invRT1 = RT1 + RT1_size;
	PBPAS::InverseRootsTable(dims1, RT1, invRT1);

	int RT2_size = (dims2<<1)-2; //2*dims2-2
	sfixn *RT2   = (sfixn *) my_calloc(RT2_size<<1, sizeof(sfixn));
	PBPAS::RootsTable(dims2, es2, RT2, pPtr);
	sfixn *invRT2 = RT2 + RT2_size;
	PBPAS::InverseRootsTable(dims2, RT2, invRT2);
	//--------------------------------
	/*
	sfixn *RT1    = (sfixn *) my_calloc(dims1<<1 + dims2<<1, sizeof(sfixn));
	sfixn *invRT1  = RT1+dims1;
	sfixn *RT2    = invRT1+dims1;
	sfixn *invRT2 = RT2+dims2;
	//--------------------------------------
	cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es1,dims1,RT1,pPtr);	
	PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2,dims2,RT2,pPtr);
	cilk_sync;

	invRT1[0] = RT1[0];
	for(int i=1;i<dims1;++i)
	    invRT1[i] = RT1[dims1-1];

	invRT2[0] = RT2[0];
	for(int i=1;i<dims2;++i)
	    invRT2[i] = RT2[dims2-1];
	*/

	//RevBidMap
	int *RevBidMap = new int[H];
	for(int i=0; i<H; ++i)   
	    RevBidMap[i] = i;
	PBPAS::RevBitInd(H, RevBidMap);

	int s = ls1*ls2;
#if __GNUC__ == 4
	memset(C, 0, s*sizeof(sfixn)); 
#else
	cilk_for(int j = 0; j < s*sizeof(sfixn); j++) C[j] = 0;
#endif
	cilk_spawn 
	    Evaluation2D(C,es1,es2,dims1,dims2,ls1,ls2,A,da[0]+1,da[1]+1,RT1,RT2,pPtr,H,RevBidMap,wp);

	
	sfixn *tmp = (sfixn *)my_calloc(s, sizeof(sfixn));
	Evaluation2D(tmp,es1,es2,dims1,dims2,ls1,ls2,B,db[0]+1,db[1]+1,RT1,RT2,pPtr,H,RevBidMap,wp);
	cilk_sync;
	//cout<<"after Evaluation2D"<<endl;

	//result is in C
	PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(s, C, tmp, pPtr);
	//cout<<"after PairwiseMul"<<endl;
	my_free(tmp);
	//cout<<"before Interpolation2D"<<endl;
	Interpolation2D(C,es1,es2,dims1,dims2,ls1,ls2,invRT1,invRT2,pPtr,H,RevBidMap,wp); 

	delete [] RevBidMap;
	my_free(RT1);
    }
    /* 
       res should be initialized to zero
       s1 = partial degree of v1 of A + 1
       s2 = partial degree of v2 of A + 1
    */
    void Evaluation2D(sfixn *res, int es1, int es2, 
		      int dims1, int dims2,
		      int ls1, int ls2, 
		      sfixn *A, int s1, int s2,
		      sfixn *RT1, sfixn *RT2,
		      MONTP_OPT2_AS_GENE *pPtr,
		      int H, int *RevBidMap, int wp)
    {//-----------------------------------------------------------
      //cout<<"in Evaluation2D"<<endl;
      //cout<<"dims1: "<<dims1<<endl;
      //cout<<"ls1: "<<ls1<<endl;
	if (dims1==ls1){//ls1 is a power of 2, use fft
	  int h = dims1>>1; //dims1/2, size of the shuffling buffer for fft
	    size_t st = s1*sizeof(sfixn);
	int *RevBidMapt = new int[H];
	for(int i=0; i<H; ++i)   
	    RevBidMapt[i] = i;
	if(H<=dims1)
	PBPAS::RevBitInd(H, RevBidMapt);
	else
	PBPAS::RevBitInd(dims1, RevBidMapt);
	    cilk_for(int i=0; i<s2; ++i){
		sfixn *SB = (sfixn *)my_malloc(h*sizeof(sfixn)); 
		//copy s1 coefficients of v1 for v2^i in A to the right vector in res
		sfixn *vec = res+i*dims1; //size of vec is a power of 2
		memcpy(vec, A+i*s1, st);
		//whichprime=0 to call the fft for a general Fourier prime number
		//cout<<"in Evaluation2D, before DFT_eff"<<endl;
		PBPAS::DFT_eff(dims1, es1, vec, RT1, pPtr, H, RevBidMapt, SB, wp);
		my_free(SB);
	    }
	}else{
	    cout<<"ToDo: size of v1 is not a power of 2. Use TFT"<<endl;
	    return;
	}
	//cout<<"in Evaluation2D, after eva v1"<<endl;
	if ( ls1==ls2 )
	    PBPAS::sqtranspose(res, ls2, 0, ls1, 0, 1, ls1);
	else{	
	    sfixn n = ls1*ls2;
	    sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
	    PBPAS::transpose(res, ls1, B, ls2, 0, ls2, 0, ls1);
	    memcpy(res, B, n*sizeof(sfixn));
	    my_free(B);
	}
	//cout<<"in Evaluation2D, before eva v2"<<endl;
	if (dims2==ls2) {//Size of v2 is a power of 2, use FFT
	int *RevBidMapt = new int[H];
	for(int i=0; i<H; ++i)   
	    RevBidMapt[i] = i;
	if(H<=dims1)
	PBPAS::RevBitInd(H, RevBidMapt);
	else
	PBPAS::RevBitInd(dims2, RevBidMapt);
	  int h = dims2>>1; 
	    cilk_for(int i=0; i<s1; ++i){
		sfixn *SB = (sfixn *)my_malloc(h*sizeof(sfixn)); 
		sfixn *vec = res+i*dims2; 
		PBPAS::DFT_eff(dims2, es2, vec, RT2, pPtr, H, RevBidMapt, SB, wp);	
		my_free(SB);
	    }
	}else{
	    cout<<"ToDo: size of v2 is not a power of 2. Use TFT"<<endl;
	    return;
	} 
	//cout<<"in Evaluation2D, after eva v2"<<endl;
    }

    void Interpolation2D(sfixn *res, int es1, int es2,
			 int dims1, int dims2, 
			 int ls1, int ls2,
			 sfixn *invRT1, sfixn *invRT2,
			 MONTP_OPT2_AS_GENE *pPtr, 
			 int H, int *RevBidMap, int wp)
    {//---------------------------------------------------------------	
      //cout<<"in Interpolation2D"<<endl;
	sfixn R=1L<<(pPtr->Rpow); 
	R=R%pPtr->P;
	
	sfixn invn2=inverseMod(dims2,pPtr->P);
	invn2 = MulMod(R,invn2,pPtr->P);
	invn2 <<= pPtr->Base_Rpow;
	
	sfixn invn1=inverseMod(dims1,pPtr->P);
	invn1 = MulMod(R,invn1,pPtr->P);
	invn1 = MulMod(R,invn1,pPtr->P);
	invn1 <<= pPtr->Base_Rpow;

	if (dims2==ls2){//size of v2 is a power of 2, use FFT
	int *RevBidMapt = new int[H];
	for(int i=0; i<H; ++i)   
	    RevBidMapt[i] = i;
	if(H<=dims1)
	PBPAS::RevBitInd(H, RevBidMapt);
	else
	PBPAS::RevBitInd(dims2, RevBidMapt);
	    int h = dims2>>1; 
	    cilk_for(int i=0; i<ls1; ++i){
		sfixn *SB = (sfixn *)my_malloc(h*sizeof(sfixn)); //shuffle buffer
		PBPAS::InvDFT_eff_keepMontgomery(dims2, es2, res+i*dims2, invRT2, pPtr, H, RevBidMapt, SB, invn2, wp);
	    }
	}else{
	    cout<<"ToDo: size of v2 is not a power of 2. Use TFT"<<endl;
	    return;
	}
	
	if ( ls1==ls2 )
	  PBPAS::sqtranspose(res, ls1, 0, ls2, 0, 1, ls2);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
	  PBPAS::transpose(res, ls2, B, ls1, 0, ls1, 0, ls2);	  
	  memcpy(res, B, n*sizeof(sfixn));
	  my_free(B);
	}
	
	if (dims1==ls1){//ls1 is a power of 2, use fft
	int *RevBidMapt = new int[H];
	for(int i=0; i<H; ++i)   
	    RevBidMapt[i] = i;
	if(H<=dims1)
	PBPAS::RevBitInd(H, RevBidMapt);
	else
	PBPAS::RevBitInd(dims1, RevBidMapt);
	    int h = dims1>>1; //dims1/2, size of the shuffling buffer for fft
	    cilk_for(int j=0; j<ls2; ++j){
		sfixn *SB = (sfixn *)my_malloc(h*sizeof(sfixn)); //shuffle buffer
		PBPAS::InvDFT_eff(dims1, es1, res+j*dims1, invRT1, pPtr, H, RevBidMapt, SB, invn1, wp);
		my_free(SB);
	    }
	}else{
	    cout<<"ToDo: size of v1 is not a power of 2. Use TFT"<<endl;
	    return;
	}
    }

	// FFT-based multiplication
	bool ks_mul (int l, int s, sfixn* f, sfixn* g, sfixn p, int basesize) {
		MONTP_OPT2_AS_GENE* pPtr = (MONTP_OPT2_AS_GENE*) my_malloc(sizeof(MONTP_OPT2_AS_GENE));
		int check = 0;
		if (p == 4179340454199820289) { check = 1; }
		else if (p == 3799912185593857) { check = 2; }

		EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
		int* RevBidMap = new int[basesize];
		for (int i = 0; i < basesize; ++i)
			RevBidMap[i] = i;
		if( s>=basesize)
		PBPAS::RevBitInd(basesize, RevBidMap);
		else
		PBPAS::RevBitInd(s, RevBidMap);

		int z = s << 1;
		sfixn* sb1 = (sfixn *) my_calloc(s, sizeof(sfixn));
		sfixn* krt = (sfixn *) my_calloc(z<<1, sizeof(sfixn));
		PBPAS::RootsTableSpe(z, l+1, krt, pPtr, 1);
		PBPAS::DFT_eff(s, l, f, krt+z, pPtr, basesize, RevBidMap, sb1, check);
		PBPAS::DFT_eff(s, l, g, krt+z, pPtr, basesize, RevBidMap, sb1, check);
		#pragma cilk_grainsize = 1024;
		cilk_for (int i = 0; i < s; ++i)
			f[i] = MulMod(f[i], g[i], p);
		sfixn R = (1L << (pPtr->Rpow)) % p;
		sfixn invn = inverseMod(s, p);
		invn = MulMod(R, invn, p);
		invn <<= pPtr->Base_Rpow;
		sfixn* ikrt = (sfixn *) my_calloc(z<<1, sizeof(sfixn));
		PBPAS::InverseRootsTable(z, krt, ikrt);
		PBPAS::InvDFT_eff(s, l, f, ikrt+z, pPtr, basesize, RevBidMap, sb1, invn, check);
		my_free(sb1);
		my_free(krt);
		my_free(ikrt);
		my_free(pPtr);
		delete [] RevBidMap;
		return 1;
	}
} //end of PBPAS
