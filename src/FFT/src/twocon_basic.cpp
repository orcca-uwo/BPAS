// -*- C++ -*-
// @author Yuzhen Xie

#if DEBUG
#include <cilktools/cilkview.h>
#endif
#include "../../../include/FFT/src/basic_routine.h"
#include "../../../include/FFT/src/modpn.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <iostream>
#include <string.h>

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef __64int_t sfixn;
//typedef __int32_t int32;
using namespace std;
namespace PBPAS {
  /*-------------------------------------------------------------------------
   N=K*M
   p=RecoveringPrime(d1, d2, K, M)
   where K=2^k, 2^+`n>=2+log_2(dK)+2M and p>=2^n+1
   A(x^(K-1),y^(d1-1)) (recursive rep: y^0, y^1, ..., y^(d1-1))
   B(x^(K-1),y^(d2-1))
   because of theta (for negcyclic convolution), K is a power of 2

   A and B are bivariate, encoded in dense recursive representation
   for y > x, so A and B appear as univariate polynomials
   with monomials (y^0,...,y^d1) and coefficients in that are univariate in x
   Assumptions:
   deg(A,x) = deg(B,x) = K-1  (may be thanks to padding either A or B)
   deg(A,y) = d1 - 1
   deg(B,y) = d2 - 1

   The output is A*B mod <x^K-1,p> and A*B mod <x^K+1,p>
  ---------------------------------------------------------------------------*/

  sfixn * TwoConvolutionModNew(sfixn *Ap, sfixn *Bp,
			       sfixn d1, sfixn d2, sfixn K,
			       MONTP_OPT2_AS_GENE *pPtr,
			       int H, int *RevBidMap,sfixn num)
  {//------------------------------------------------------------
    //for A*B mod x^K+1 mod p and A*B mod x^K-1 mod p
    //K=2^k=dims1=ls1
    int es1 = PBPAS::ceil_log2_long(K); //logceiling(K);
    int es2K = es1+1;
    int K2 = K<<1;      //length of theta vector

    int ls2 = d1+d2; // y-size of A*B mod x^K + 1
    int es2 = logceiling(ls2);
    int dims2 = 1<<es2; //dimension (power of 2) of y //y-size of the 2-D FFT of the result

    //std::cout<<"es2, dims2: "<<es2<<"," <<dims2<<std::endl;
    //std::cout<<"K, H: "<<K<<"," <<H<<std::endl;

    if ( (es2K > (pPtr->Npow))
	 || (es2 > (pPtr->Npow))
	 || (es1 > (pPtr->Npow)) ){
      std::cout<<"BPAS: exception, the prime number is too small for the FFT problem."<<std::endl;
      exit(1); //checkFrourierPrime
    }

    // for FFT (omega and gamma), theta^(2K-1), (theta^(-1))^0=1:
    // first (K2+1) slots for the successive powers of theta and 1
    // next 2K-2 slots, roots table omega (fft x-direction)
    // next 2*dims2-2 slots, roots table gamma (fft y-direction)
    // next 2K-2 slots, roots table omega (invfft x-direction)
    // last 2*dims2-2 slots, roots table gamma (invfft y-direction)

    //compute roots for the largest dimension, then find the right pos for others
    //two roots tables and scaling factors for negacyclic convolution
    int K2RT_size = (K2<<1)-2; //2*K2-2
    int dRT_size = (dims2<<1)-2; //2*dims2-2

    sfixn *thetaPtr = NULL;
    sfixn *invthetaPtr = NULL;
    sfixn *dRT = NULL;
    sfixn *KRT = NULL;
    sfixn *invKRT = NULL;
    sfixn *invdRT = NULL;
    sfixn *RT = NULL;

    if (es2>=es2K){
      dRT = (sfixn *)my_calloc((dRT_size<<1), sizeof(sfixn));
      RT = dRT;
      PBPAS::RootsTable(dims2, es2, dRT, pPtr);
      invdRT = dRT + dRT_size;
      PBPAS::InverseRootsTable(dims2, dRT, invdRT);

      int g = es2 - es2K;
      thetaPtr = dRT;
      invthetaPtr = invdRT;

      int i = 0;
      while (i<g){
	int dims2i = dims2>>i;
	thetaPtr += dims2i;
	invthetaPtr += dims2i;

	i++;
      }
      invKRT = invthetaPtr + K2;
      KRT = thetaPtr + K2;

    }else{
      thetaPtr = (sfixn *)my_calloc(K2RT_size<<1, sizeof(sfixn));
      RT = thetaPtr;
      PBPAS::RootsTable(K2, es2K, thetaPtr, pPtr);
      invthetaPtr = thetaPtr + K2RT_size;
      PBPAS::InverseRootsTable(K2, thetaPtr, invthetaPtr);

      KRT = thetaPtr + K2;
      invKRT = invthetaPtr + K2;

      int g = es2K - es2;
      dRT = thetaPtr;
      invdRT = invthetaPtr;

      int i = 0;
      while (i<g){
	int K2i = K2>>i;
	dRT += K2i;
	invdRT += K2i;
	i++;
      }
    }

    //sfixn s = K*ls2; //size of the result

    sfixn s = K*dims2; //size of intermediate computation
    // real size of the result from convolution
    // first s slots for the result of negacyclic convolution mod p
    // next s slots for the result of cyclic convolution mod p
    sfixn *ncc = (sfixn *)my_calloc(s<<1, sizeof(sfixn));
    sfixn *cc = ncc+s;
    sfixn R=1L<<(pPtr->Rpow);
    R=R%pPtr->P;

    sfixn invn=inverseMod(dims2,pPtr->P);
    invn = MulMod(R,invn,pPtr->P);
    invn <<= pPtr->Base_Rpow;
    sfixn invn1=inverseMod(K,pPtr->P);
    invn1 = MulMod(R,invn1,pPtr->P);
    invn1 = MulMod(R,invn1,pPtr->P);
    invn1 <<= pPtr->Base_Rpow;
	#if DEBUG
	unsigned long time = -__cilkview_getticks();
	#endif
    cilk_spawn
    PBPAS::NegacyclicConvolution(ncc, s, es1, es2, K, dims2, Ap, Bp, d1, d2, K2, thetaPtr, KRT, dRT, invKRT, invdRT, pPtr, H, RevBidMap,invn,invn1,num);
    PBPAS::CyclicConvolution(cc, s, es1, es2, K, dims2, Ap, Bp, d1, d2, KRT, dRT, invKRT, invdRT, pPtr, H, RevBidMap,invn,invn1,num);
    cilk_sync;
	#if DEBUG
	time += __cilkview_getticks();
	cout << "cyclic and negacyclic (ms): " <<"\t" << time /1.f << endl;
	#endif

    //std::cout<<" CyclicConvolution \n";
    my_free(RT);

    return ncc;
    //next is to CRT for ncc and cc and then Recovering Product for (ncc, K, ls2, M, N);
  }

  void CyclicConvolution(sfixn *res, sfixn s,
			 sfixn es1, sfixn es2,
			 sfixn K, sfixn dims2,
			 sfixn *A, sfixn *B, sfixn dA, sfixn dB,
			 sfixn *KRT, sfixn *dRT,
			 sfixn *invKRT, sfixn *invdRT,
			 MONTP_OPT2_AS_GENE *pPtr,
			 int H, int *RevBidMap, sfixn invn, sfixn invn1,sfixn num)
  {//---------------------------------------------------------------
    //evaluation
//    cout<<"Cyclic conv"<<endl;
	#if DEBUG
	unsigned long time = -__cilkview_getticks();
	#endif
    cilk_spawn
      AdaptiveEvaluation(res, es1, es2, K, dims2, A, dA, KRT, dRT, pPtr, H, RevBidMap, num);
//	unsigned long end = __cilkview_getticks();
//	cout << "adaptive evaluation 1 (ms): " <<"\t" << (end - start) /1.f << endl;

    sfixn *tmp = (sfixn *)my_calloc(s, sizeof(sfixn));
//	start = __cilkview_getticks();
    AdaptiveEvaluation(tmp, es1, es2, K, dims2, B, dB, KRT, dRT, pPtr, H, RevBidMap, num);

    cilk_sync;
	#if DEBUG
	time += __cilkview_getticks();
	cout << "adaptive evaluation in cyclic (ms): " <<"\t" << time /1.f << endl;
	#endif

    /*
//    std::cout<<"CyclicConvolution, AdaptiveEvaluation, res: \n";
    for (sfixn i=0; i<s; i++){
//      std::cout<<res[i]<<", ";
    }
//    std::cout<<"\n";

//    std::cout<<"CyclicConvolution, AdaptiveEvaluation, tmp: \n";
    for (sfixn i=0; i<s; i++){
//      std::cout<<tmp[i]<<", ";
    }
//    std::cout<<"\n";
    */
    //pairwise multiply
	#if DEBUG
	time = -__cilkview_getticks();
	#endif
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(s, res, tmp, pPtr);
	#if DEBUG
	time += __cilkview_getticks();
	cout << "pairwise mul in cyclic (ms): " <<"\t" << time /1.f << endl;
	#endif
    /*
//    std::cout<<"CyclicConvolution, AdaptiveEvaluation, PairwiseMul res: \n";
    for (sfixn i=0; i<s; i++){
//      std::cout<<res[i]<<", ";
    }
//    std::cout<<"\n";
    */

    my_free(tmp);

    //interpolation
	#if DEBUG
	time = -__cilkview_getticks();
	#endif
    AdaptiveInterpolation(res, K, es1, es2, dims2, invKRT, invdRT, pPtr, H, RevBidMap,invn,invn1,num);
	#if DEBUG
	time += __cilkview_getticks();
	cout << "interpolation in cyclic (ms): " <<"\t" << time /1.f << endl;
	#endif
//    sfixn R=1L<<(pPtr->Rpow);
//    R=R%pPtr->P;
//    for(sfixn j=0;j<s;j++){
//		res[j] = MulMod(res[j],R,pPtr->P);
//		//printf("%lu ",res[j]);
//	}
	//printf("\n");
  }

  // dA = deg(A,y)+1
  void AdaptiveEvaluation(sfixn *res, sfixn es1, sfixn es2,
			  sfixn K, sfixn dims2,
			  sfixn *A, sfixn dA,
			  sfixn *KRT, sfixn *dRT,
			  MONTP_OPT2_AS_GENE *pPtr,
			  int H, int *RevBidMap,sfixn num)
  {//----------------------------------------------------------------
    //K is a power of 2, use FFT to evaluate in x-direction

    sfixn Kh = K>>1; //K/2

    cilk_for(sfixn i=0; i<dA; ++i){
      sfixn *SB1 = (sfixn *)my_malloc(Kh*sizeof(sfixn));
      sfixn m = i*K;
      //copy K coefficients of x for y^i to the right vector in res
      sfixn *vec = res+m;
      memcpy(vec, A+m, K*sizeof(sfixn));
      PBPAS::DFT_eff(K, es1, vec, KRT, pPtr, H, RevBidMap, SB1,num);
      my_free(SB1);
    }
    if ( K==dims2 )
      PBPAS::sqtranspose(res, dims2, 0, K, 0, 1, K);
    else{
      //sfixn n = K*ls2;
      sfixn n = K*dims2;
      sfixn *B = (sfixn *) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, K, B, dims2, 0, dims2, 0, K);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    sfixn Dh = dims2>>1; //dims2/2

    cilk_for(sfixn i=0; i<K; ++i){
      sfixn *SB2 = (sfixn *)my_malloc(Dh*sizeof(sfixn)); //shuffle buffer
      sfixn *resi = res+i*dims2;
      PBPAS::DFT_eff(dims2, es2, resi, dRT, pPtr, H, RevBidMap, SB2,num);
      my_free(SB2);
    }

  }

  void AdaptiveInterpolation(sfixn *res, sfixn K, sfixn es1,
			     sfixn es2, sfixn dims2,
			     sfixn *invKRT, sfixn *invdRT,
			     MONTP_OPT2_AS_GENE *pPtr,
			     int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num)
  {//---------------------------------------------------------------
    //input res' layout is ls2 x K, from evaluation
    //interpolate in y (d direction)
    /*
//    std::cout<<"---AdaptiveInterpolation "<<"\n";

//    std::cout<<"\t contents of invKRT: "<<"\n";
    for (sfixn u=0; u<((K<<1)-2); u++){
//      std::cout<<invKRT[u]<<", ";
    }
//    std::cout<<"\n";

//    std::cout<<"\t contents of invdRT: "<<"\n";
    for (sfixn u=0; u<((dims2<<1)-2); u++){
//      std::cout<<invdRT[u]<<", ";
    }
//    std::cout<<"\n";


//    std::cout<<"---AdaptiveInterpolation: over dims2 for K num "<<"\n";
    */

    //dims2 is a power of 2, use FFT
    sfixn Dh = dims2>>1; //dims2/2

    cilk_for(sfixn j=0; j<K; ++j){
      sfixn *SB2 = (sfixn *)my_malloc(Dh*sizeof(sfixn)); //shuffle buffer
      sfixn *resj = res+j*dims2;

      /*/--------------------------
      sfixn sum = 0;
      for (sfixn u=0; u<dims2; u++){
	sum = sum+resj[u];
	sum = sum % pPtr->P;
      }
//      std::cout<<"j, in sum mod: "<<j<<", "<<sum<<"\n";
      */
      PBPAS::InvDFT_eff_keepMontgomery(dims2, es2, resj, invdRT, pPtr, H, RevBidMap, SB2,invn,num);
//      //std::cout<<"j, out resj[0]: "<<j<<", "<<resj[0]<<"\n";

      my_free(SB2);
    }


    //transform
    if ( K==dims2 )
      PBPAS::sqtranspose(res, K, 0, dims2, 0, 1, dims2);
    else{
      sfixn n = K*dims2;
      sfixn *B = (sfixn * ) my_malloc(n*sizeof(sfixn));
      PBPAS::transpose(res, dims2, B, K, 0, K, 0, dims2);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    //interpolate in x, K is a power of 2
//    //std::cout<<"---AdaptiveInterpolation: over K for dims2 num "<<"\n";

    sfixn Kh = K>>1; //K/2

    cilk_for(sfixn j=0; j<dims2; ++j){
      sfixn *SB1 = (sfixn *)my_malloc(Kh*sizeof(sfixn));

      sfixn *resj = res+j*K;

      /*/--------------------------
      sfixn sum = 0;
      for (sfixn u=0; u<K; u++){
	sum = sum+resj[u];
	sum = sum % pPtr->P;
      }
//      std::cout<<"j, in sum mod: "<<j<<", "<<sum<<"\n";
      */
      PBPAS::InvDFT_eff(K, es1, resj, invKRT, pPtr, H, RevBidMap, SB1,invn1,num);
//      //std::cout<<"j, out resj[0]: "<<j<<", "<<resj[0]<<"\n";

      my_free(SB1);
    }
    //now res' layout is K x ls2
  }

  void NegacyclicConvolution(sfixn *res, sfixn s,
			     sfixn es1, sfixn es2,
			     sfixn K, sfixn dims2,
			     sfixn *A, sfixn *B,
			     sfixn dA, sfixn dB, sfixn K2,
			     sfixn *thetaPtr,
			     sfixn *KRT, sfixn *dRT,
			     sfixn *invKRT, sfixn *invdRT,
			     MONTP_OPT2_AS_GENE *pPtr,
			     int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num)
  {//--------------------------------------------------------------
    //apply weight vector and evaluate
//    cout<<"negacyclic convolution"<<endl;
	#if DEBUG
	unsigned long time = -__cilkview_getticks();
	#endif
    cilk_spawn
      WeightVectorAdaptiveEvaluation(res, es1, es2, K, dims2, A, dA, thetaPtr, KRT, dRT, pPtr, H, RevBidMap,num);
//	time += __cilkview_getticks();
//	cout << "adaptive evaluation 1 :" <<time<<endl;
    sfixn *tmp = (sfixn *)my_calloc(s, sizeof(sfixn));
//	time = -__cilkview_getticks();
    WeightVectorAdaptiveEvaluation(tmp, es1, es2, K, dims2, B, dB, thetaPtr, KRT, dRT, pPtr, H, RevBidMap,num);
    cilk_sync;
	#if DEBUG
	time += __cilkview_getticks();
	cout << "adaptive evaluation negacyclic:" <<time<<endl;
	#endif

    //pairwise multiply
	#if DEBUG
	time = -__cilkview_getticks();
	#endif
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(s, res, tmp, pPtr);
    //PBPAS::EX_Mont_PairwiseMul_OPT2_AS(s, res, tmp, pPtr->P);//MulMod
	#if DEBUG
	time += __cilkview_getticks();
	cout << "pairwise :" <<time<<endl;
	#endif

    my_free(tmp);
	#if DEBUG
	time = -__cilkview_getticks();
	#endif
    AdaptiveInterpolationWeightVector(res, K, es1, es2, dims2, thetaPtr, invKRT, invdRT, pPtr, H, RevBidMap,invn, invn1,num);
	#if DEBUG
	time += __cilkview_getticks();
	cout << "interpolation in negacyclic :" <<time<<endl;
	#endif
//    sfixn R=1L<<(pPtr->Rpow);
//    R=R%pPtr->P;
//    for(sfixn j=0;j<s;j++){
//		res[j] = MulMod(res[j],R,pPtr->P);
//		//printf("%lu ",res[j]);
//	}
	//printf("\n");
  }

  void WeightVectorAdaptiveEvaluation(sfixn *res, sfixn es1,
				      sfixn es2, sfixn K,
				      sfixn dims2,
				      sfixn *A, sfixn dA,
				      sfixn *thetaPtr,
				      sfixn *KRT, sfixn *dRT,
				      MONTP_OPT2_AS_GENE *pPtr,
				      int H, int *RevBidMap,sfixn num)
  {//------------------------------------------------------------
    sfixn Kh = (K>>1); //K/2 shuffle buffer

    //store each weighted vector (coeff of x) of y^i to res and then FFT inplace
    //cilk_for(sfixn i=0; i<dA; i++){
    cilk_for(sfixn i=0; i<dA; i++){
      sfixn *SB1 = (sfixn *)my_malloc(Kh*sizeof(sfixn));
      sfixn m = i*K;
      sfixn *A_vec = A+m;
      sfixn *res_vec = res+m;
      for(sfixn j=0; j<K; j++){
	res_vec[j]=MontMulMod_OPT2_AS_GENE(A_vec[j], thetaPtr[j], pPtr); //weight vector
      }
	sfixn sum =0;
      PBPAS::DFT_eff(K, es1, res_vec, KRT, pPtr, H, RevBidMap, SB1,num);
      my_free(SB1);
    }
    //transform
    if ( K==dims2 )
      PBPAS::sqtranspose(res, dims2, 0, K, 0, 1, K);
    else{
      sfixn n = K*dims2;
      sfixn *B = (sfixn *) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, K, B, dims2, 0, dims2, 0, K);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    //evaluate in y

    //Size of y dims2 is a power of 2, use FFT in-place
    sfixn Dh = dims2>>1; //dims2/2

    //cilk_for(sfixn i=0; i<K; ++i){
    cilk_for(sfixn i=0; i<K; ++i){
      sfixn *SB2 = (sfixn *)my_malloc(Dh*sizeof(sfixn)); //shuffle buffer
      PBPAS::DFT_eff(dims2, es2, res+i*dims2, dRT, pPtr, H, RevBidMap, SB2,num);
      my_free(SB2);
    }
    //now res' layout is ls2 x K
  }

  void AdaptiveInterpolationWeightVector(sfixn *res, sfixn K,
					 sfixn es1, sfixn es2,
					 sfixn dims2,
					 sfixn *thetaPtr,
					 sfixn *invKRT,
					 sfixn *invdRT,
					 MONTP_OPT2_AS_GENE *pPtr,
					 int H, int *RevBidMap,sfixn invn,sfixn invn1,sfixn num)
  {//---------------------------------------------------------------
    //input res' layout is ls2 x K, from evaluation
    //interpolate in y (d direction)

    //dims2 is a power of 2, use FFT
    sfixn Dh = dims2>>1; //dims2/2

    //cilk_for(sfixn j=0; j<K; ++j){
    cilk_for(sfixn j=0; j<K; ++j){
      sfixn *SB2 = (sfixn *)my_malloc(Dh*sizeof(sfixn)); //shuffle buffer
      PBPAS::InvDFT_eff_keepMontgomery(dims2, es2, res+j*dims2, invdRT, pPtr, H, RevBidMap, SB2,invn,num);
      my_free(SB2);
    }

    //transform
    if ( K==dims2 )
      PBPAS::sqtranspose(res, K, 0, dims2, 0, 1, dims2);
    else{
      sfixn n = K*dims2;
      sfixn *B = (sfixn *)my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, dims2, B, K, 0, K, 0, dims2);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    //interpolate in x, K is a power of 2
    //then unweight

    sfixn Kh = K>>1; //K/2

    sfixn K2 = K<<1; //2K
    //cilk_for(sfixn j=0; j<dims2; ++j){
    cilk_for(sfixn j=0; j<dims2; ++j){
      sfixn *SB1 = (sfixn *)my_malloc(Kh*sizeof(sfixn));
      sfixn *res_vec = res+j*K;
      PBPAS::InvDFT_eff(K, es1, res_vec, invKRT, pPtr, H, RevBidMap, SB1,invn1,num);
      for(sfixn j=0; j<K; j++){//_OPT2_AS_GENE
	//res_vec[j]=MontMulMod_OPT2_AS_GENE(res_vec[j], invthetaPtr[j], pPtr);
	res_vec[j]=MontMulMod_OPT2_AS_GENE(res_vec[j], thetaPtr[K2-j], pPtr);

	//the ith-powers of theta's inverse is backward from thetaPtr[K2] to [K2-K+1]
      }
      my_free(SB1);
    }
  }



  /*old~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *old~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *old~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  sfixn * TwoConvolutionMod(sfixn *Ap, sfixn *Bp,
			     sfixn d1, sfixn d2,
			     sfixn p, sfixn K)
  {//---------------------------------------------------------------
    //debug
    /*
    std::ofstream ofs("BivarModP.input", std::ofstream::out);
    for(sfixn j = 0; j < K*d1; ++j)
      ofs << Ap[j] << " ";
    ofs << "Ap\n\n";
    for(sfixn j = 0; j < K*d2; ++j)
      ofs << Bp[j] << " ";
    ofs << "Bp\n\n";
    ofs.close();
    */

    //for A*B mod x^K+1 mod p and A*B mod x^K-1 mod p
    //K=2^k=dims1=ls1
    sfixn es1 = PBPAS::ceil_log2_long(K); //logceiling(K);
    sfixn es2K = es1+1;
    sfixn K2 = K<<1;      //length of theta vector

    sfixn ls2 = d1+d2-1; // y-size of A*B mod x^K + 1
    sfixn es2 = logceiling(ls2);
    sfixn dims2 = 1<<es2; //dimension (power of 2) of y //y-size of the 2-D FFT of the result

    MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
//    //std::cout<<"pPtr->Npow: "<<pPtr->Npow<<std::endl;
    if ( (es2K > (pPtr->Npow))
	 || (es2 > (pPtr->Npow))
	 || (es1 > (pPtr->Npow)) ){
//      std::cout<<"BPAS: exception, the prime number is too small for the FFT problem."<<std::endl;
      exit(1); //checkFrourierPrime
    }

    // for FFT (omega and gamma), theta^(2K-1), (theta^(-1))^0=1:
    // first K slots for the successive powers of omega (x-direction)
    // next dims2 slots for the successive powers of gamma  (y-direction)
    // next (K2+1) slots for the successive powers of theta and 1
    sfixn sr = K+dims2+K2+1;
    sfixn *rootsPtr = (sfixn *) my_calloc(sr, sizeof(sfixn));
    sfixn *rootsPtr_2 = rootsPtr+K;
    sfixn *thetaPtr = rootsPtr_2+dims2;

    //gamma
    cilk_spawn
      PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dims2, rootsPtr_2, pPtr);

    // theta0 theta^(2*K-1)
    PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2K, K2, thetaPtr, pPtr);
    thetaPtr[K2] = thetaPtr[0]; //special Mongomery rep for (theta^(-1))^0=1,

    cilk_sync;

    //omega
    for (sfixn i=0; i<K; i++){ //use cilk_for or not? test to see
      rootsPtr[i] = thetaPtr[i<<1];
    }

    //debug
//    //std::cout<<"rootsPtr"<<std::endl;
    //for (sfixn i=0; i<sr; i++){ //use cilk_for or not? test to see
//    // std::cout<<rootsPtr[i]<<", ";
    //}
//    //std::cout<<std::endl;

    sfixn s = K*ls2; //size of result // real size of the result from convolution
    // first s slots for the result of negacyclic convolution mod p
    // next s slots for the result of cyclic convolution mod p
    sfixn *ncc = (sfixn *)my_calloc(s<<1, sizeof(sfixn));
    sfixn *cc = ncc+s;

    cilk_spawn
//	unsigned long time = -__cilkview_getticks();
      PBPAS::NegacyclicConvolution1(ncc, s, es1, es2, K, dims2, ls2, Ap, Bp, d1, d2, K2, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"negacyclic convolution " <<time <<endl;
//	time = -__cilkview_getticks();
    PBPAS::CyclicConvolution1(cc, s, es1, es2, K, dims2, ls2, Ap, Bp, d1, d2, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"cyclic convolution " <<time <<endl;
    cilk_sync;

    my_free(pPtr);
    my_free(rootsPtr);

    return ncc;
    //next is to CRT for ncc and cc and then Recovering Product for (ncc, K, ls2, M, N);
  }

  void CyclicConvolution1(sfixn *res, sfixn s,
			  sfixn es1, sfixn es2,
			  sfixn K, sfixn dims2, sfixn ls2,
			  sfixn *A, sfixn *B, sfixn dA, sfixn dB,
			  sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr)
  {//--------------------------------------------------------------------------
    //evaluation
//    cout<<"cyclic "<<endl;
    cilk_spawn
//	unsigned long time =-__cilkview_getticks();
      AdaptiveEvaluation1(res, es1, es2, K, dims2, ls2, A, dA, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"adaptive evaluation 1 "<<time <<endl;
    sfixn *tmp = (sfixn *)my_calloc(s, sizeof(sfixn));
//	time =-__cilkview_getticks();
    AdaptiveEvaluation1(tmp, es1, es2, K, dims2, ls2, B, dB, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"adaptive evaluation 1 "<<time <<endl;
    cilk_sync;

    //pairwise multiply
    //PBPAS::EX_Mont_PairwiseMul_OPT2_AS(s, res, tmp, pPtr->P); //output is in res

//	time =-__cilkview_getticks();
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(s, res, tmp, pPtr); //test??
//	time+=__cilkview_getticks();
//	cout<<"pairwise muls  "<<time <<endl;
    //sfixn BRsft=pPtr->Base_Rpow;
    //cilk_for(sfixn i=0; i<s; i++)
    // res[i]=MontMulMod_OPT2_AS_GENE(res[i], tmp[i]<<BRsft, pPtr);

    my_free(tmp);

    //interpolation
//	time =-__cilkview_getticks();
    AdaptiveInterpolation1(res, K, es1, es2, dims2, ls2, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"interpolation  "<<time <<endl;
  }

  // dA = deg(A,y)+1
  void AdaptiveEvaluation1(sfixn *res, sfixn es1, sfixn es2,
			   sfixn K, sfixn dims2, sfixn ls2,
			   sfixn *A, sfixn dA,
			   sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr)
  {//------------------------------------------------------------------
    // evaluate in x
    //size_t sx = K*sizeof(sfixn);

    //#pragma cilk_grainsize = 8192; //need tune
    //K is a power of 2, use FFT to evaluate in x-direction
    cilk_for(sfixn i=0; i<dA; ++i){
      sfixn m = i*K;
      //copy the K coefficients of x for y^i to the right vector in res
      sfixn *vec = res+m;
      memcpy(vec, A+m, K*sizeof(sfixn));
      EX_Mont_DFT_OPT2_AS_GENE_1(K, es1, rootsPtr, vec, pPtr); //inplace FFT
      //EX_Mont_TDFT_OPT2_AS_GENE_1(K, rootsPtr, vec, pPtr); //inplace TFT
    }

    //transform
    if ( K==ls2 )
      PBPAS::sqtranspose(res, ls2, 0, K, 0, 1, K);
    else{
      sfixn n = K*ls2;
      sfixn *B = (sfixn *) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, K, B, ls2, 0, ls2, 0, K);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    //evaluate in y
    sfixn *rootsPtr2 = rootsPtr+K;
    if (dims2==ls2) {//Size of y is a power of 2, use FFT
      cilk_for(sfixn i=0; i<K; ++i){
	EX_Mont_DFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, res+i*ls2, pPtr); //inplace FFT
	//EX_Mont_TDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, res+i*ls2, pPtr); //inplace TFT
      }
    }else{//ls2 is not a power of 2, use TFT
      //size_t s = ls2*sizeof(sfixn);
      cilk_for(sfixn j=0; j<K; ++j){
	sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn)); //ToDo: allocate once
	sfixn * pos = res+j*ls2;
	//memcpy(tmpVec, pos, s);
	memcpy(tmpVec, pos, ls2*sizeof(sfixn));
	EX_Mont_TDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, tmpVec, pPtr);
	//EX_Mont_DFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, tmpVec, pPtr);

	//memcpy(pos, tmpVec, s);
	memcpy(pos, tmpVec, ls2*sizeof(sfixn));
	my_free(tmpVec);
      }
    }
    //now res' layout is ls2 x K
  }

  void AdaptiveInterpolation1(sfixn *res, sfixn K, sfixn es1,
			      sfixn es2, sfixn dims2, sfixn ls2,
			      sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr)
  {//---------------------------------------------------------------------
    //input res' layout is ls2 x K, from evaluation

    //interpolate in y
    sfixn * rootsPtr2 = rootsPtr+K;

    if (dims2==ls2){//ls2 is a power of 2, use FFT
      cilk_for(sfixn j=0; j<K; ++j)
	EX_Mont_INVDFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, res+j*ls2, pPtr);
	//EX_Mont_INVTDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, res+j*ls2, pPtr);

	//EX_Mont_INVDFT_OPT2_AS_GENE_R_1(ls2, es2, rootsPtr2, res+j*ls2, pPtr);

    }else{//use out-of-place TFT
      size_t s = ls2*sizeof(sfixn);
      cilk_for(sfixn j=0; j<K; j++){
	sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn));
	sfixn * pos = res+j*ls2;
	memcpy(tmpVec, pos, s);
	//memcpy(tmpVec, pos, ls2*sizeof(sfixn));
	EX_Mont_INVTDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, tmpVec, pPtr);//taken

	//EX_Mont_INVDFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, tmpVec, pPtr);

	//EX_Mont_INVTDFT_OPT2_AS_GENE_R_1(ls2, rootsPtr2, tmpVec, pPtr);
	memcpy(pos, tmpVec, s);
	//memcpy(pos, tmpVec, ls2*sizeof(sfixn));
	my_free(tmpVec);
      }
    }
    //transform
    if ( K==ls2 )
      PBPAS::sqtranspose(res, K, 0, ls2, 0, 1, ls2);
    else{
      sfixn n = K*ls2;
      sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, ls2, B, K, 0, K, 0, ls2);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    //interpolate in x, K is a power of 2
    cilk_for(sfixn j=0; j<ls2; ++j)
      //EX_Mont_INVDFT_OPT2_AS_GENE_1(K, es1, rootsPtr, res+j*K, pPtr); //inplace inverse FFT
      //EX_Mont_INVTDFT_OPT2_AS_GENE_1(K, rootsPtr, res+j*K, pPtr); //inplace inverse TFT
      EX_Mont_INVDFT_OPT2_AS_GENE_R_1(K, es1, rootsPtr, res+j*K, pPtr);
    //now res' layout is K x ls2
  }

  void NegacyclicConvolution1(sfixn *res, sfixn s,
			      sfixn es1, sfixn es2,
			      sfixn K, sfixn dims2, sfixn ls2,
			      sfixn *A, sfixn *B, sfixn dA, sfixn dB, sfixn K2,
			      sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr)
  {//---------------------------------------------------------------------------

    //apply weight vector and evaluate
//    cout << "negacyclic"<<endl;
    cilk_spawn
//	unsigned long time= -__cilkview_getticks();
      WeightVectorAdaptiveEvaluation1(res, es1, es2, K, K2, dims2, ls2, A, dA, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"adaptive evaluation 1 "<<time <<endl;
    sfixn *tmp = (sfixn *)my_calloc(s, sizeof(sfixn));
//	time= -__cilkview_getticks();
    WeightVectorAdaptiveEvaluation1(tmp, es1, es2, K, K2, dims2, ls2, B, dB, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"adaptive evaluation 2 "<<time <<endl;
    cilk_sync;

    //pairwise multiply
    //PBPAS::EX_Mont_PairwiseMul_OPT2_AS(s, res, tmp, pPtr->P); //use MulMod, output is in res
//	time= -__cilkview_getticks();
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(s, res, tmp, pPtr); //test??
//	time+=__cilkview_getticks();
//	cout<<"componentwise muls "<<time <<endl;

    my_free(tmp);

//	time= -__cilkview_getticks();
    AdaptiveInterpolationWeightVector1(res, K, es1, es2, dims2, ls2, rootsPtr, pPtr);
//	time+=__cilkview_getticks();
//	cout<<"interpolation "<<time <<endl;

  }

  void WeightVectorAdaptiveEvaluation1(sfixn *res, sfixn es1, sfixn es2,
				       sfixn K, sfixn K2,
				       sfixn dims2, sfixn ls2,
				       sfixn *A, sfixn dA,
				       sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr)
  {//---------------------------------------------------------------------------
    sfixn *thetaPtr = rootsPtr+K+dims2;

    //store each weighted vector (coeff of x) of y^i to res and then FFT inplace
    cilk_for(sfixn  i=0; i<dA; i++){

      sfixn m = i*K;
      sfixn *A_vec = A+m;
      sfixn *res_vec = res+m;
      for(sfixn j=0; j<K; j++){
	res_vec[j]=MontMulMod_OPT2_AS_GENE(A_vec[j], thetaPtr[j], pPtr); //weight vector
	//res_vec[j]=MontMulMod(A_vec[j], thetaPtr[j], pPtr);
      }
      EX_Mont_DFT_OPT2_AS_GENE_1(K, es1, rootsPtr, res_vec, pPtr); //inplace FFT
      //EX_Mont_TDFT_OPT2_AS_GENE_1(K, rootsPtr, res_vec, pPtr); //inplace TFT
    }

    //transform
    if ( K==ls2 )
      PBPAS::sqtranspose(res, ls2, 0, K, 0, 1, K);
    else{
      sfixn n = K*ls2;
      sfixn *B = (sfixn *) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, K, B, ls2, 0, ls2, 0, K);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    //evaluate in y
    sfixn *rootsPtr2 = rootsPtr+K;
    if (dims2==ls2) {//Size of y is a power of 2, use FFT
      cilk_for(sfixn i=0; i<K; ++i){
	EX_Mont_DFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, res+i*ls2, pPtr); //inplace FFT
	//EX_Mont_TDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, res+i*ls2, pPtr);
      }
    }else{//ls2 is not a power of 2, use TFT
      size_t s = ls2*sizeof(sfixn);
      cilk_for(sfixn j=0; j<K; ++j){
	sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn));
	sfixn * pos = res+j*ls2;

	memcpy(tmpVec, pos, s);
	//memcpy(tmpVec, pos, ls2*sizeof(sfixn));
	EX_Mont_TDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, tmpVec, pPtr);
	//EX_Mont_DFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, tmpVec, pPtr);
	memcpy(pos, tmpVec, s);
	//memcpy(pos, tmpVec, ls2*sizeof(sfixn));
	my_free(tmpVec);
      }
    }
    //now res' layout is ls2 x K
  }

void AdaptiveInterpolationWeightVector1(sfixn *res, sfixn K, sfixn es1,
					sfixn es2, sfixn dims2, sfixn ls2,
					sfixn *rootsPtr,
					MONTP_OPT2_AS_GENE *pPtr)
{//---------------------------------------------------------------------
    //input res' layout is ls2 x K, from evaluation

    //interpolate in y
    sfixn * rootsPtr2 = rootsPtr+K;

    if (dims2==ls2){//ls2 is a power of 2, use FFT
      cilk_for(sfixn j=0; j<K; ++j)
	EX_Mont_INVDFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, res+j*ls2, pPtr);
	//EX_Mont_INVTDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, res+j*ls2, pPtr);

	//EX_Mont_INVDFT_OPT2_AS_GENE_R_1(ls2, es2, rootsPtr2, res+j*ls2, pPtr);

    }else{//use out-of-place TFT
      size_t s = ls2*sizeof(sfixn);
      cilk_for(sfixn j=0; j<K; j++){
	sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn));
	sfixn * pos = res+j*ls2;
	memcpy(tmpVec, pos, s);
	//memcpy(tmpVec, pos, ls2*sizeof(sfixn));
	EX_Mont_INVTDFT_OPT2_AS_GENE_1(ls2, rootsPtr2, tmpVec, pPtr);//taken

	//EX_Mont_INVDFT_OPT2_AS_GENE_1(ls2, es2, rootsPtr2, tmpVec, pPtr);

	//EX_Mont_INVTDFT_OPT2_AS_GENE_R_1(ls2, rootsPtr2, tmpVec, pPtr);
	memcpy(pos, tmpVec, s);
	//memcpy(pos, tmpVec, ls2*sizeof(sfixn));
	my_free(tmpVec);
      }
    }
    //transform
    if ( K==ls2 )
      PBPAS::sqtranspose(res, K, 0, ls2, 0, 1, ls2);
    else{
      sfixn n = K*ls2;
      sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(res, ls2, B, K, 0, K, 0, ls2);
      memcpy(res, B, n*sizeof(sfixn));
      my_free(B);
    }

    sfixn *thetaPtr = rootsPtr+K+dims2;
    sfixn K2 = K<<1;

    //sfixn R=1L<<(pPtr->Rpow);
    //R=R%pPtr->P;

    //sfixn R2 = MulMod(R,R,pPtr->P);

    //interpolate in x, K is a power of 2
    cilk_for(sfixn j=0; j<ls2; ++j){
      sfixn *res_vec = res+j*K;
      //EX_Mont_INVDFT_OPT2_AS_GENE_1(K, es1, rootsPtr, res_vec, pPtr); //inplace inverse FFT
      //EX_Mont_INVTDFT_OPT2_AS_GENE_1(K, rootsPtr, res+j*K, pPtr); //inplace inverse TFT
      EX_Mont_INVDFT_OPT2_AS_GENE_R_1(K, es1, rootsPtr, res+j*K, pPtr);

      for(sfixn j=0; j<K; j++){//_OPT2_AS_GENE
	res_vec[j]=MontMulMod_OPT2_AS_GENE(res_vec[j], thetaPtr[K2-j], pPtr);

	//sfixn tmp = MontMulMod_OPT2_AS_GENE(res_vec[j], thetaPtr[K2-j], pPtr);
	//res_vec[j]= MontMulMod_OPT2_AS_GENE(tmp, R2<<=pPtr->Base_Rpow, pPtr);
	//the ith-powers of theta's inverse is backward from thetaPtr[K2] to [K2-K+1]
      }
    }
    //now res' layout is K x ls2
  }


} //end of PBPAS
