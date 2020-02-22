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
#include "../../../include/ModularPolynomial/src/general_routine.h"

#include "./example_util_gettime.h"

//typedef __int64_t   sfixn;

using namespace std;

int n = 1024;//256;//1024;////16;
int r = 10;//8;//10;//
int H = 8;//8;//16;//4; //16; //4; 
int v = 1;

//test new recursive DFT verse the MODPN DFT 
int main(int argc, char *argv[])
{//------------------------------------------
  cout << "#workers:\t" << __cilkrts_get_nworkers() << endl;
  //if(argc > 1) d1 = atoi(argv[1]);
  
  //cilkview_data_t start_data, end_data;
  //unsigned long long start = __cilkview_getticks();
  //__cilkview_query(start_data);

  //sfixn p = 180143985094819841; //ifactor(180143985094819841 - 1)
  //sfixn p = 7881299347898369; //53 bits
  //sfixn p = 99194853094755497; //57 bits, Npow=3

  //sfixn p = 17; //test new DFT
  sfixn myp = 257;

  //recommanded ones
  //sfixn p = 3553621580972033; //52 bits, 2^45 * 101 + 1
  //sfixn p = 3799912185593857; //52 bits, 2^47 * 27 + 1
  
  //sfixn p = 10645405947658241; //54 bits, 2^27 * 79314455 + 1, correct cyclic
  //8529266726959324, 8529266726959324, 8529266726959324, 8529266726959324, 8, 8, 8, 8, 2116143784101677, 2116143784101677, 2116143784101677, 2116143784101677

  //sfixn p = 10645410511060993; //54 bits, 2^27 * 79314489 + 1
  //sfixn p = 24686718519083009; //55 bits, 2^28 * 91965193 + 1
  //sfixn p = 24686728719630337; //55 bits, 2^28 * 91965231 + 1
  //sfixn p = 38349804861915137; //56 bits, 2^30 * 35716039 +1

  cout<<"The prime: " <<myp<<endl;
 
  //test RevBidInd
  int * RevBidMap = (int *)my_malloc(H*sizeof(sfixn));
  
  for(int i=0; i<H; ++i)   
      RevBidMap[i] = i;
  /*
  std::cout<<"RevBidMap In: "<<std::endl;
  for(int i=0; i<H; ++i)
    std::cout<<RevBidMap[i]<<", ";
  std::cout<<std::endl;
  */
  PBPAS::RevBitInd(H, RevBidMap);
  /*
  std::cout<<"RevBidMap Out: "<<std::endl;
  for(int i=0; i<H; ++i)
    std::cout<<RevBidMap[i]<<", ";
  std::cout<<std::endl;
  */
  MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
  
  EX_MontP_Init_OPT2_AS_GENE(pPtr, myp);
  

  //test RootsTable
  int size = (n<<1) - 2; //2n-2
  sfixn *T = (sfixn *)my_malloc(size*sizeof(sfixn));
  sfixn *Ti = (sfixn *)my_malloc(size*sizeof(sfixn));

  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(r, n, T, pPtr);

  //for(int i=0; i<n; ++i)
  //  T[i] = i+1;

  //std::cout<<"T In: "<<std::endl;
  //for(int i=0; i<size; ++i)
  //  std::cout<<T[i]<<", ";
  //std::cout<<std::endl;

  PBPAS::RootsTable(n,r,T,pPtr);
  
  /*
  std::cout<<"T Out: "<<std::endl;
  int k=0;
  for(int i=0; i<r; ++i){
    k += n>>i;
    std::cout<<T[k]<<", ";
  }
  std::cout<<std::endl;
  std::cout<<std::endl;
  for(int i=0; i<size; ++i)
    std::cout<<T[i]<<", ";
  std::cout<<std::endl;
  */
  PBPAS::InverseRootsTable(n,T,Ti);
  /*
  std::cout<<"Ti Out: "<<std::endl;
  int ki=0;
  for(int i=0; i<r; ++i){
    ki += n>>i;
    std::cout<<Ti[ki]<<", ";
  }
  std::cout<<std::endl;
  std::cout<<std::endl;
  for(int i=0; i<size; ++i)
    std::cout<<Ti[i]<<", ";
  std::cout<<std::endl;
  */
  
  //test new DFT
  sfixn *A = (sfixn *)my_calloc(n, sizeof(sfixn));
  sfixn *As = (sfixn *)my_calloc(n>>1, sizeof(sfixn));
  sfixn *B = (sfixn *)my_calloc(n, sizeof(sfixn));
  sfixn *C = (sfixn *)my_calloc(n, sizeof(sfixn));

  srand(getSeed());

  for(int i=0; i<n; ++i){ 
    int rv = rand() % myp;
    A[i] = rv;
    B[i] = rv;
    C[i] = rv;
  }
  
  /*
  cout<<"\n---C in: \n";
  for(int i=0; i<n; ++i)
    cout<<C[i]<<", ";
  std::cout<<std::endl;
  */
  /*
  PBPAS::Shuffle(n,C,As);
  cout<<"\n---Shuffled C: \n";
  for(int i=0; i<n; ++i)
    cout<<C[i]<<", ";
  std::cout<<std::endl;
  */

  cout<<"\n---sum A in: \n";
  sfixn sum = 0;
  for(int i=0; i<n; ++i){
    sum += A[i];
    sum = sum %myp;
    //cout<<A[i]<<", ";
  } 
  std::cout<<"sum: "<<sum<<std::endl;
  
  
  //unsigned long long start = __cilkview_getticks();
 
  //long start = example_get_time();

  PBPAS::DFT_eff(n,r,A,T,pPtr,H,RevBidMap,As);

  //long end = example_get_time();
  //float mytime = (end - start) / 1.f;
  //std::cout << "DFT_eff time (ms): "<< mytime << std::endl;

  //unsigned long long end = __cilkview_getticks();
  
  //cout << "DFT_eff time (ms): " <<"\t" << (end - start) /1.f << endl;
  
  cout<<"\n---dft A[0]: "<< A[0]<<endl;
 
  //start = __cilkview_getticks();

  //start = example_get_time();

  EX_Mont_DFT_OPT2_AS_GENE_1(n, r, T, B, pPtr);
  //EX_Mont_INVDFT_OPT2_AS_GENE_1(n, r, T, B, pPtr);

  //end = example_get_time();
  //mytime = (end - start) / 1.f;
  // std::cout << "EX_Mont_DFT time (ms): "<< mytime << std::endl;

  //end = __cilkview_getticks();
  //cout << "EX_Mont_DFT time (sec): " <<"\t" << (end - start) /1000.f << endl;
 
  //test against MODPN
  int * RevBidMap_n = (int *)calloc(n, sizeof(sfixn));
  
  for(int i=0; i<n; ++i)   
      RevBidMap_n[i] = i;
  
  PBPAS::RevBitInd(n, RevBidMap_n);

  PBPAS::ArrayBitReversal(n,B,RevBidMap_n);
  cout<<"\n---EX_Mont_DFT B[0]: "<<B[0]<<endl;
  
  cout<<"\n---B out after RevBidMap_n: \n";
  for(int i=0; i<n; ++i)
    cout<<B[i]<<", ";
  std::cout<<std::endl;
  
  int flag=1;
  for(int i=0; i<n; ++i){
    if (A[i]!=B[i]){
      std::cout<<"DFT A, EX_Mont_DFT B Does not match at i, A, B: "<<i<<", "<<A[i]<<", "<<B[i]<<std::endl;
      flag=0;
      break;
    }
  }
  if (flag==1)
    std::cout<<"***DFT A, EX_Mont_DFT B matches"<<std::endl;

  PBPAS::InvDFT_eff(n,r,A,Ti,pPtr,H,RevBidMap,As);

  cout<<"\n---invdft A[0]: "<<A[0]<<endl;
  cout<<"\n---invdft A out sum: \n";
  sum = 0;
  for(int i=0; i<n; ++i){
    sum += A[i];
    sum = sum %myp;
    //cout<<A[i]<<", ";
  }
  std::cout<<"sum: "<<sum<<std::endl;
  
  /*
  cout<<"\n--- invdft A out: \n";
  for(int i=0; i<n; ++i)
    cout<<A[i]<<", ";
  std::cout<<std::endl;
  */
  flag=1;
  for(int i=0; i<n; ++i){
    if (A[i]!=C[i]){
      std::cout<<"Does not match at i, invA, C: "<<i<<", "<<A[i]<<", "<<C[i]<<std::endl;
      flag=0;
      break;
    }
  }
  if (flag==1)
    std::cout<<"***invDFT A, C matches"<<std::endl;

  my_free(RevBidMap);
  
  my_free(pPtr);
  my_free(T);
  my_free(Ti);
  my_free(A);
  my_free(As);
  my_free(B);
  my_free(C);

  return 0;
}

