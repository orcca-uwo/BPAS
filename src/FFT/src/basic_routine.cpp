// -*- C++ -*- 
// basic_routine.cilk
//@auther Yuzhen Xie

#include <cilk/cilk.h>
#include "../../../include/FFT/src/modpn.h"
#include "../../../include/FFT/src/modpn_export.h"

#include "../../../include/FFT/src/general_routine.h"
#include "../../../include/FFT/src/basic_routine.h"

#include "../../../include/FFT/src/ser_general_routine.h"
#include "../../../include/FFT/src/ser_basic_routine.h"

#include "../../../include/FFT/src/include/example_util_gettime.h"
#include <iostream>
//#include <reducer_opadd.h>

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef int sfixn;

// sfixn MulCut=42;
// sfixn DivCut1=146;
// sfixn DivCut2=5;
// sfixn DivCut3=3;
// sfixn DivCutN=2;

extern sfixn MulCut;
extern sfixn DivCut1;
extern sfixn DivCut2;
extern sfixn DivCut3;
extern sfixn DivCutN;

namespace PBPAS {


  /**
   * fftMultiD_test_1:
   * @coeffs1: coefficient vector for 'f1'.
   * @coeffs2: coefficient vector for 'f2'.
   * @N: number of variables.
   * @es: 2^es[i] = dims [i] for i = 1..N.
   *      es[i] is the exponent of dims[i];
   * @dims: the FFT sizes on each dimension, dims[i] is the FFT size
   *        on dimension i, i=1..N.
   * @pPtr: the information of the prime number.
   * es[0] and dims[0] are useless;
   * seperated evaluation by FFT
   * Return value: 
   **/
  void fftMultiD_test_1(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr){
    register int i;
    sfixn m=0, n=1;
    sfixn *rootsPtr, *tmprootsPtr;
	sfixn *tmprootsPtr1;
	
	//std::cout << "entering fftMultiD_test_1" << std::endl;

	//check Frourier prime
	for(i=1; i<=N; i++) {
	  PBPAS::checkFrourierPrime(es[i], pPtr);
	}

    //compute the primitive roots for all dimensions
    for(i=1; i<=N; i++) { n *= dims[i]; m += dims[i]; } //use reducer??

    rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));

	//:( :(
	//if use tmprootsPtr in computing the roots here, got
	// cilk++ compile error:
	//mulfft.cilk:570: error: invalid operand to binary operator
	//tmprootsPtrD.35412
	//mulfft.cilk:570: internal compiler error: verify_stmts failed

    tmprootsPtr1 = rootsPtr;	

    for(i=1; i<=N; i++) {//parallel???
      PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr1, pPtr);
      tmprootsPtr1 += dims[i];
    }
	//------------------------------------
// 	cilk_for(int j=0; j<N; j++) { //segmentation fault
// 	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[j+1], dims[j+1], rootsPtr+j*dims[j], pPtr);
// 	}

    //make a copy of es and dims for poly2's FFT use in parallel 
    sfixn *es_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    sfixn *dims_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));

    for(int j=1; j<=N; j++) { es_cp[j] = es[j]; dims_cp[j] = dims[j]; }

#ifdef TIMEFN
	float pairwisemul_time = 0.0;
	float pareval_time = 0.0;
	float interpol_time = 0.0;
	long start1, end1;
	start1 = example_get_time();
#endif

#ifdef SERIALFT
	tmprootsPtr = PBPAS::MultiEvalByFFT(coeffs1, N, n, es, dims, pPtr, rootsPtr);
	PBPAS::MultiEvalByFFT(coeffs2, N, n, es_cp, dims_cp, pPtr, rootsPtr); 
#else	

    tmprootsPtr = cilk_spawn PBPAS::MultiEvalByFFT(coeffs1, N, n, es, dims, pPtr, rootsPtr);

    PBPAS::MultiEvalByFFT(coeffs2, N, n, es_cp, dims_cp, pPtr, rootsPtr); 

    cilk_sync;
#endif

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pareval_time += (end1 - start1) / 1.f;
#endif

	my_free(es_cp);
    my_free(dims_cp);

#ifdef TIMEFN
	start1 = example_get_time();
#endif
    // Pairwise-Mul, result is in coeffs1
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs2, pPtr->P);

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pairwisemul_time += (end1 - start1) / 1.f;
#endif

#ifdef TIMEFN
	start1 = example_get_time();
#endif
    PBPAS::InterpolByFFT(coeffs1, N, n, es, dims, pPtr, tmprootsPtr);
#ifdef TIMEFN  
	  end1 = example_get_time();
	  interpol_time += (end1 - start1) / 1.f;
#endif

	//std::cout <<"before exit fftMultiD_test_1" << std::endl;
    my_free(rootsPtr);
#ifdef TIMEFN
	std::cout << pareval_time << " " << pairwisemul_time << " "<< interpol_time << " " ;
#endif
  }

  /**-------------------------------------------------------
   *
   **/
  sfixn * MultiEvalByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr){
    register int i;
    //sfixn j, tmp;
	sfixn tmp, d;
	sfixn * tmprootsPtr;

#ifdef TIMEFN
	float tran_time = 0.0;
	float dft_time = 0.0;
	long start, end;
#endif
	
	tmprootsPtr = rootsPtr;

#ifdef TIMEFN	  
	  start = example_get_time();
#endif	
    if(es[1]){ //mulit-DFT evaluation for dimension 1
//       for(j=0; j<n; j += dims[1]){
// 		cilk_spawn EX_Mont_DFT_OPT2_AS_GENE_1 ( dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);      
//       }  
// 	  cilk_sync;
//------------------------------------
	  d = n/dims[1];
	  cilk_for(int k=0; k<d; k++){
		EX_Mont_DFT_OPT2_AS_GENE_1 ( dims[1], es[1], tmprootsPtr, coeffs1+k*dims[1], pPtr);
	  }
    }
#ifdef TIMEFN  
	  end = example_get_time();
	  dft_time += (end - start) / 1.f;
#endif
	
    for(i=2; i<=N; i++){
      tmprootsPtr += dims[1];

#ifdef TIMEFN
	  start = example_get_time();
#endif
      PBPAS::multi_mat_transpose (N, n, i, dims, coeffs1);
#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif


#ifdef TIMEFN	  
	  start = example_get_time();
#endif
      if( es[i] ){
// 		for(j=0; j<n; j += dims[i]){
// 		  cilk_spawn EX_Mont_DFT_OPT2_AS_GENE_1 (dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
// 		}
// 		cilk_sync;
//--------------------------------------------
		d = n/dims[i];
		cilk_for(int k=0; k<d; k++){
		  EX_Mont_DFT_OPT2_AS_GENE_1 (dims[i], es[i], tmprootsPtr, coeffs1+k*dims[i], pPtr);
		}
      }
#ifdef TIMEFN  
	  end = example_get_time();
	  dft_time += (end - start) / 1.f;
#endif	  

      //multi_mat_transpose didn't change dims accordingly, so we do it here.
      tmp=dims[1];
      dims[1]=dims[i];
      dims[i]=tmp;
      tmp=es[1];
      es[1]=es[i];
      es[i]=tmp;
    }
#ifdef TIMEFN  	
    //std::cout << "MultiEvalByFFT Mat Tran took " << tran_time << "
    //ms." << std::endl;
	std::cout << tran_time << " "<<dft_time << " ";
#endif	
	//return tmprootsPtr;
	sfixn *retPtr = tmprootsPtr;
	return retPtr;
  }

  /**---------------------------------------------------
   *InterpolByFFT: inverse DFT
   **/
  void InterpolByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr){
    register int i;
    sfixn tmp, d;

#ifdef TIMEFN
	float tran_time = 0.0;
	float invdft_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    if( es[1] ){
//       for(j=0; j<n; j+=dims[1]){
// 		cilk_spawn EX_Mont_INVDFT_OPT2_AS_GENE_1( dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);
//       }
// 	  cilk_sync;
//----------------------------------
	  d = n/dims[1];
	  cilk_for(int j=0; j<d; j++){
		EX_Mont_INVDFT_OPT2_AS_GENE_1( dims[1], es[1], tmprootsPtr, coeffs1+j*dims[1], pPtr);
	  }		
    }	
#ifdef TIMEFN  
	  end = example_get_time();
	  invdft_time += (end - start) / 1.f;
#endif

    for(i=N; i>=2; i--){
      tmprootsPtr -= dims[i];
	  //tmprootsPtr -= dims[1];
#ifdef TIMEFN
	  start = example_get_time();
#endif
	  PBPAS::multi_mat_transpose (N, n, i, dims, coeffs1);

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
      if( es[i] ){
// 		for(int j=0; j<n; j+=dims[i]){
// 		  cilk_spawn EX_Mont_INVDFT_OPT2_AS_GENE_1(dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
// 		}
// 		cilk_sync;
//-----------------------------------
		d = n/dims[i];
		cilk_for(int j=0; j<d; j++){
		  EX_Mont_INVDFT_OPT2_AS_GENE_1(dims[i], es[i], tmprootsPtr, coeffs1+j*dims[i], pPtr);
		}
      }
#ifdef TIMEFN  
	  end = example_get_time();
	  invdft_time += (end - start) / 1.f;
#endif	  
      //multi_mat_transpose didn't change dims accordingly, so we do
      //it here.
      tmp=dims[1];
      dims[1]=dims[i];
      dims[i]=tmp;
      tmp=es[1];
      es[1]=es[i];
      es[i]=tmp;
    }
#ifdef TIMEFN  
	std::cout << tran_time << " "<< invdft_time << " ";
    //std::cout << "InterpolByFFT Mat Tran took " << tran_time << " ms." << std::endl;
#endif	
  }
  
  //-------------------------------------
  //N=2 
  void MultiplyBy2DTFT_no_transp(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
	if (N !=2 ) {
	  std::cout<<"Number of variables is not 2"<<std::endl;
	  exit(1);
	}
	
#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);	

    cilk_spawn PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));

    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
	cilk_sync;

#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif

	  //---------------------------
	  PBPAS::tft2D_no_transp(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);

#ifdef TIMEFN
	  start = example_get_time();
#endif					
	  PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
	  PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	  std::cout << initkn_time << " "<< reckn_time << " ";
#endif

  }

  //===============================================================
  //N=2
  void tft2D_no_transp(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr){

	register int i;
	sfixn m=0, n=1, maxdim=0;
	sfixn * rootsPtr, * tmprootsPtr;
	sfixn * tmprootsPtr1;

	for(i=1; i<=N; i++) {
	  checkFrourierPrime(es[i], pPtr);
	}

	for(i=1; i<=N; i++) {//N is usually small
	  n *= ls[i]; 
      m += dims[i]; 
      if(dims[i] > maxdim) maxdim=dims[i]; 
    }
	rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
    tmprootsPtr1 = rootsPtr;
	//compute the primitive roots for all dimensions
	for(i=1; i<=N; i++) {
	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr1, pPtr);
      tmprootsPtr1 += dims[i];
    }
    //make a copy of es, dims and ls for poly2's TFT evaluation in parallel     
	sfixn *es_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    sfixn *dims_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    sfixn *ls_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
	for(i=1; i<=N; i++) { 
	  es_cp[i] = es[i];
	  dims_cp[i] = dims[i]; 
	  ls_cp[i] = ls[i];
	}

#ifdef TIMEFN
	float pairwisemul_time = 0.0;
	float pareval_time = 0.0;
	float interpol_time = 0.0;
	long start1, end1;
	start1 = example_get_time();
#endif
	//parallel run two
    cilk_spawn MultiEvalBy2DTFT_no_transp(coeffs1, n, maxdim, es, dims, ls, pPtr, rootsPtr);

	MultiEvalBy2DTFT_no_transp(coeffs2, n, maxdim, es_cp, dims_cp, ls_cp, pPtr, rootsPtr);
	
	cilk_sync;

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pareval_time += (end1 - start1) / 1.f;
#endif

    my_free(es_cp);
    my_free(dims_cp);
	my_free(ls_cp);	

#ifdef TIMEFN
	start1 = example_get_time();
#endif
	// Pairwise-Mul, result is in coeffs1
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs2, pPtr->P);

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pairwisemul_time += (end1 - start1) / 1.f;
#endif	

#ifdef TIMEFN
	start1 = example_get_time();
#endif
	
    PBPAS::InterpolBy2DTFT_no_transp(coeffs1, n, maxdim, es, dims, ls, pPtr, rootsPtr);
#ifdef TIMEFN  
	  end1 = example_get_time();
	  interpol_time += (end1 - start1) / 1.f;
#endif
    my_free(rootsPtr);

#ifdef TIMEFN
	std::cout << pareval_time << " " << pairwisemul_time << " "<< interpol_time << " " ;
#endif

  }

  void InterpolBy2DTFT_no_transp(sfixn *coeffs1, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr) {
	
	sfixn d;
	sfixn * tmprootsPtr = rootsPtr;

#ifdef TIMEFN
	float invtft_time = 0.0;
	long start, end;
	start = example_get_time();
#endif

	if( es[1] ){
	  d = n/ls[1];
	  cilk_for(int j=0; j<d; j++){
		PBPAS::EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(coeffs1+j*ls[1],  maxdim, ls[1], tmprootsPtr, pPtr);
	  }
	}

	if( es[2] ){	
		d = n/ls[2];
		tmprootsPtr += dims[1];
		PBPAS::inv_TDFT_dim2_no_transp(d, ls[2], tmprootsPtr, maxdim, coeffs1, pPtr);
	}
#ifdef TIMEFN  
	  end = example_get_time();
	  invtft_time += (end - start) / 1.f;
	  std::cout << invtft_time << " ";
#endif		

  }

  //================================================================
  void inv_TDFT_dim2_no_transp(sfixn d, sfixn ls2, sfixn * tmprootsPtr, sfixn maxdim, sfixn * coeffs1, MONTP_OPT2_AS_GENE *pPtr) {

	cilk_for(int i=0; i<d; i++){
	  sfixn * tmpVec = (sfixn *)my_calloc(maxdim, sizeof(sfixn));

	  cilk_for (int j=0; j<ls2; j++) {
		tmpVec[j] = coeffs1[i+j*d];
	  }
	  //invtft
	  EX_Mont_INVTDFT_OPT2_AS_GENE_1( ls2, tmprootsPtr, tmpVec, pPtr );

	  cilk_for (int j=0; j<ls2; j++) {
		coeffs1[i+j*d] = tmpVec[j];
	  }

	  my_free(tmpVec);
	}
  }

  //================================================================
  //N=2
  void MultiEvalBy2DTFT_no_transp(sfixn * coeffs1, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr) {
	
	//register int i;
    //sfixn tmp;
	sfixn * tmprootsPtr;
	sfixn d;

	tmprootsPtr = rootsPtr;

#ifdef TIMEFN
	float tft_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN	  
	  start = example_get_time();
#endif

	if( es[1] ) { //mulit-DFT evaluation for dimension 1
	  d = n/ls[1];
	  cilk_for(int j=0; j<d; j++){
		PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1_par(ls[1], tmprootsPtr, maxdim, coeffs1+j*ls[1], pPtr);
	  }
	}

	if( es[2] ){
	  tmprootsPtr += dims[1]; //move pointer to dim-2
	  d = n/ls[2]; //size of dim-1
	  PBPAS::TDFT_dim2_no_transp(d, ls[2], tmprootsPtr, maxdim, coeffs1, pPtr);
	}
	
#ifdef TIMEFN  
	  end = example_get_time();
	  tft_time += (end - start) / 1.f;
	  std::cout << tft_time << " ";
#endif	

  }
  
  //================================================================
  //d = ls[1]  
  void TDFT_dim2_no_transp(sfixn d, sfixn ls2, sfixn * tmprootsPtr, sfixn maxdim, sfixn * coeffs1, MONTP_OPT2_AS_GENE *pPtr) {
	cilk_for(int i=0; i<d; i++){
	  sfixn * tmpVec = (sfixn *)my_calloc(maxdim, sizeof(sfixn));

	  cilk_for (int j=0; j<ls2; j++) {
		tmpVec[j] = coeffs1[i+j*d];
	  }
	  EX_Mont_TDFT_OPT2_AS_GENE_1 ( ls2, tmprootsPtr, tmpVec, pPtr);
	  cilk_for (int j=0; j<ls2; j++) {
		coeffs1[i+j*d] = tmpVec[j];
	  }
	  my_free(tmpVec);
	}

  }

  //================================================================
  // Parallel Multi-dimensional TFT by seperated evaluation 
  //================================================================
  /**
   * tftMultiD_test_1:
   * @coeffs1: coefficient vector for 'f1'.
   * @coeffs2: coefficient vector for 'f2'.
   * @N: number of variables.
   * @es: 2^es[i] = dims [i] for i = 1..n.
   * @dims: the FFT sizes on each dimension, dims[i] is the FFT on dimensional i, i=1..N.
   * @pPtr: the information of the prime number.
   * 
   *Notes: each TDFT uses a temp vector of size maxdim;
   *       each INVTDFT uses a temp vector of size maxdim
   *
   * Return value: 
   **/
  void tftMultiD_test_1(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr){
	register int i;
	sfixn m=0, n=1, maxdim=0;
	sfixn * rootsPtr, * tmprootsPtr;
	sfixn * tmprootsPtr1;
	//std::cout <<"enter tftMultiD_test_1" << std::endl;

	//check Frourier prime
	for(i=1; i<=N; i++) {
	  checkFrourierPrime(es[i], pPtr);
	}

	for(i=1; i<=N; i++) {//N is usually small
	  n *= ls[i]; 
      m += dims[i]; 
      if(dims[i] > maxdim) maxdim=dims[i]; 
    }

	rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
    tmprootsPtr1 = rootsPtr;
	//compute the primitive roots for all dimensions
	for(i=1; i<=N; i++) {
	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr1, pPtr);
      tmprootsPtr1 += dims[i];
    }
	//----------------------------	
	// 	cilk_for(int j=0; j<N; j++) { //free(): invalid next size (fast)
	// 	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[j+1], dims[j+1], rootsPtr+j*dims[j], pPtr);
	// 	}
    
    //make a copy of es, dims and ls for poly2's TFT evaluation in parallel     
	sfixn *es_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    sfixn *dims_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    sfixn *ls_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
	for(i=1; i<=N; i++) { 
	  es_cp[i] = es[i];
	  dims_cp[i] = dims[i]; 
	  ls_cp[i] = ls[i];
	}
	//----------------------------
// 	cilk_for(int j=1; j<=N; j++) { //N is usually small
// 	  es_cp[j] = es[j]; 
// 	  dims_cp[j] = dims[j]; 
// 	  ls_cp[j] = ls[j]; 
// 	}

	//std::cout <<"before MultiEvalByTFT" << std::endl;

#ifdef TIMEFN
	float pairwisemul_time = 0.0;
	float pareval_time = 0.0;
	float interpol_time = 0.0;
	long start1, end1;
	start1 = example_get_time();
#endif

#ifdef SERIALFT
	tmprootsPtr = MultiEvalByTFT(coeffs1, N, n, maxdim, es, dims, ls, pPtr, rootsPtr);
	MultiEvalByTFT(coeffs2, N, n, maxdim, es_cp, dims_cp, ls_cp, pPtr, rootsPtr);
#else
	//parallel run two
    tmprootsPtr = cilk_spawn MultiEvalByTFT(coeffs1, N, n, maxdim, es, dims, ls, pPtr, rootsPtr);

	MultiEvalByTFT(coeffs2, N, n, maxdim, es_cp, dims_cp, ls_cp, pPtr, rootsPtr);	
	cilk_sync;

#endif

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pareval_time += (end1 - start1) / 1.f;
#endif

    my_free(es_cp);
    my_free(dims_cp);
	my_free(ls_cp);

	//std::cout <<"after MultiEvalByTFT" << std::endl;

#ifdef TIMEFN
	start1 = example_get_time();
#endif
	// Pairwise-Mul, result is in coeffs1
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs2, pPtr->P);

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pairwisemul_time += (end1 - start1) / 1.f;
#endif
	//my_free(coeffs2); //???????

#ifdef TIMEFN
	start1 = example_get_time();
#endif
	//is tmprootsPtr1 pointing to the last dimension??? 
    PBPAS::InterpolByTFT(coeffs1, N, n, maxdim, es, dims, ls, pPtr, tmprootsPtr);
#ifdef TIMEFN  
	  end1 = example_get_time();
	  interpol_time += (end1 - start1) / 1.f;
#endif
	//	std::cout <<"before exit tftMultiD_test_1" << std::endl;
    my_free(rootsPtr);

#ifdef TIMEFN
	std::cout << pareval_time << " " << pairwisemul_time << " "<< interpol_time << " " ;
#endif
	//std::cout <<"after free before exit tftMultiD_test_1" << std::endl;
  }


  /**-------------------------------------------------------
   *
   **/
  sfixn * MultiEvalByTFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr) {
	register int i;
    sfixn tmp;
	sfixn * tmprootsPtr;
	sfixn d;

	//std::cout <<"enter MultiEvalByTFT" << std::endl;
	//for(int j = 0; j < n; j++) std::cout<< coeffs1[j]<<", ";
	//std::cout<<"----end of initial data----" <<std::endl;
	
#ifdef TIMEFN
	float tran_time = 0.0;
	float tft_time = 0.0;
	long start, end;
#endif

	tmprootsPtr = rootsPtr;

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
	if( es[1] ) { //mulit-DFT evaluation for dimension 1
//       for(j=0; j<n; j += ls[1]){//each TDFT uses a temp vector of size maxdim
// 		cilk_spawn PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1(ls[1], tmprootsPtr, maxdim, coeffs1+j, pPtr);
// 	  }
// 	  cilk_sync;
//--------------------------------------
	  d = n/ls[1];
	  cilk_for(int j=0; j<d; j++){
		PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1_par(ls[1], tmprootsPtr, maxdim, coeffs1+j*ls[1], pPtr);
	  }
	}
#ifdef TIMEFN  
	  end = example_get_time();
	  tft_time += (end - start) / 1.f;
#endif
		
	//std::cout <<"after ls[1]" << std::endl;

	for(i=2; i<=N; i++){

      tmprootsPtr += dims[1]; //dims[i] changes to dims[1]

#ifdef TIMEFN
	  start = example_get_time();
#endif	

	  //std::cout<<"----in data----" <<std::endl;
	  //for (int j = 1; j <= N; j++) {
	  //	std::cout<< ls[j] <<", ";
	  //}
	  //std::cout<<std::endl;

	  //for(int j = 0; j < n; j++) std::cout<< coeffs1[j]<<", ";
	  //std::cout<<"----end of in data----" <<std::endl;
	  
	  //------------------------------
      PBPAS::multi_mat_transpose (N, n, i, ls, coeffs1);
	  
	  //for(int j = 0; j < n; j++) std::cout<< coeffs1[j]<<", ";
	  //std::cout<<"----end of out data----" <<std::endl;

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif
	  
	  //std::cout <<"after multi_mat_transpose" << std::endl;

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
	  if( es[i] ){ //each TDFT uses a temp vector of size maxdim
		d = n/ls[i];
		cilk_for(int j=0; j<d; j++){
		  PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1_par(ls[i], tmprootsPtr, maxdim, coeffs1+j*ls[i], pPtr);
		}
	  }
#ifdef TIMEFN  
	  end = example_get_time();
	  tft_time += (end - start) / 1.f;
#endif
	  //std::cout <<"after ls[i]: " << i << std::endl;
	  //multi_mat_transpose didn't change dims accordingly, so we do
	  //it here.
	  tmp=dims[1];
      dims[1]=dims[i];
      dims[i]=tmp;

      tmp=es[1];
      es[1]=es[i];
      es[i]=tmp;

      tmp=ls[1];
      ls[1]=ls[i];
      ls[i]=tmp;
	}
#ifdef TIMEFN  
	std::cout << tran_time << " "<<tft_time << " ";
    //std::cout << "MultiEvalByTFTSer Mat Tran took " << tran_time << " ms." << std::endl;
#endif
	//return tmprootsPtr;
	sfixn *tmpPtr = tmprootsPtr;
	return tmpPtr;
  }

  /**-------------------------------------------------------
   *
   **/
  void EX_Mont_TDFT_OPT2_AS_GENE_1_par(sfixn lsi, sfixn * tmprootsPtr, sfixn maxdim, sfixn * coeffs1j, MONTP_OPT2_AS_GENE *pPtr){

	//for parallel run
	sfixn * tmpVec = (sfixn *)my_calloc(maxdim, sizeof(sfixn));

	PBPAS::copyVec_0_to_d (lsi - 1, tmpVec, coeffs1j);
	EX_Mont_TDFT_OPT2_AS_GENE_1 ( lsi, tmprootsPtr, tmpVec, pPtr);
	PBPAS::copyVec_0_to_d (lsi - 1, coeffs1j, tmpVec);

	my_free(tmpVec);
  }


  //------------------------------
  void EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(sfixn * coeffs1j, sfixn maxdim, sfixn lsi, sfixn * tmprootsPtr, MONTP_OPT2_AS_GENE *pPtr) {

	//for parallel run
	sfixn * tmpVec = (sfixn *)my_calloc(maxdim, sizeof(sfixn));

	PBPAS::copyVec_0_to_d(lsi-1, tmpVec, coeffs1j); 
	EX_Mont_INVTDFT_OPT2_AS_GENE_1( lsi, tmprootsPtr, tmpVec, pPtr);
	PBPAS::copyVec_0_to_d(lsi-1, coeffs1j, tmpVec);
	
	my_free(tmpVec);
  }

  
  /**---------------------------------------------------
   *InterpolByTFT: inverse DFT
   **/
  void InterpolByTFT(sfixn *coeffs1, sfixn N, sfixn n, sfixn maxdim, sfixn * es, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr) {
	register int i;
    sfixn d, tmp;

#ifdef TIMEFN
	float tran_time = 0.0;
	float invtft_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
	if( es[1] ){
	  d = n/ls[1];
	  cilk_for(int j=0; j<d; j++){
		PBPAS::EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(coeffs1+j*ls[1],  maxdim, ls[1], tmprootsPtr, pPtr);
	  }
	}
#ifdef TIMEFN  
	  end = example_get_time();
	  invtft_time += (end - start) / 1.f;
#endif

	for(i=N; i>=2; i--){

      tmprootsPtr -= dims[i];

#ifdef TIMEFN
	  start = example_get_time();
#endif
	  
	  PBPAS::multi_mat_transpose (N, n, i, ls, coeffs1);

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
	  if( es[i] ){	
		d = n/ls[i];
		cilk_for(int j=0; j<d; j++){//each INVTDFT uses a temp vector of size maxdim
		  PBPAS::EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(coeffs1+j*ls[i],  maxdim, ls[i], tmprootsPtr, pPtr);
		}
	  }
#ifdef TIMEFN  
	  end = example_get_time();
	  invtft_time += (end - start) / 1.f;
#endif
	  // multi_mat_transpose didn't change dims accordingly, 
      // so we need to do this here.
      tmp=dims[1];
      dims[1]=dims[i];
      dims[i]=tmp;

      tmp=es[1];
      es[1]=es[i];
      es[i]=tmp;

      tmp=ls[1];
      ls[1]=ls[i];
      ls[i]=tmp;
	}
#ifdef TIMEFN  
	std::cout << tran_time << " "<< invtft_time << " ";
    //std::cout << "InterpolByTFT Mat Tran took " << tran_time << " ms." << std::endl;
#endif
  }


  //-------------------------------------
  void MultiplyByFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    KroFFTRep * kPtr = (KroFFTRep *)my_calloc(1, sizeof(KroFFTRep));
    PBPAS::InitKroFFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

	PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif
      
    //----------------------------------------------
    PBPAS::fftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), pPtr);
	
	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );

#ifdef TIMEFN
	  start = example_get_time();
#endif				
    PBPAS::fromtofftRepMultiD(N, CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    //void freeKroFFTRep(KroFFTRep * x);
    PBPAS::freeKroFFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }

  //-------------------------------------
  //ToDo:can save by using rep12 for the evaluation of rep1
  //can save by not TFT on zeros in the evaluation
  void MultiplyByKronecker1DTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));

#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif

	  sfixn deg1 = 1;
	  for (int i=1; i<=N; i++) {
		deg1 *= LSI(kPtr, i);
	  }
	  deg1 = deg1-1;

	  std::cout<<deg1<<" ";
	  
	  sfixn * ES_1D = (sfixn * )my_calloc(2, sizeof(sfixn));
	  sfixn * DIMS_1D=(sfixn * )my_calloc(2, sizeof(sfixn));
	  sfixn * LS_1D=(sfixn *)my_calloc(2, sizeof(sfixn));
	  
	  ES_1D[1] = logceiling(deg1+1);
	  DIMS_1D[1] = 1 << ES_1D[1];
	  LS_1D[1] = deg1+1;

	  //---------------------------
	  //1-dimensional
	  PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), 1, ES_1D, DIMS_1D, LS_1D, pPtr);

	  my_free(ES_1D);
	  my_free(DIMS_1D);
	  my_free(LS_1D);


#ifdef TIMEFN
	  start = example_get_time();
#endif					
    PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }



  //-------------------------------------
  void MultiplyByKroneckerFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
    sfixn kdg1=0, kdg2=0;
    
    KroFFTRep * kPtr = (KroFFTRep *)my_calloc(1, sizeof(KroFFTRep));
    PBPAS::InitKroFFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);
    
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));

	//-----------
	checkFrourierPrime(KE(kPtr), pPtr);

    kdg1= (BUSZSI(rep1, N)+1)*CUMI(kPtr, N)-1;
    kdg2= (BUSZSI(rep2, N)+1)*CUMI(kPtr, N)-1;

    KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));

    //------------------------------------------------------------------------
    PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr),pPtr);

    //------------------------------------------------------------------------
    PBPAS::EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), kdg2, DATSI(kPtr, 1), pPtr);

	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );
				

    PBPAS::fromtofftRepMultiD(N, CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    //void freeKroFFTRep(KroFFTRep * x);
    PBPAS::freeKroFFTRep(kPtr);    
  }

 /**
 * EX_Mont_FFTMul_OPT2_AS_GENE_1: in-place.
 * @n: the FFT size
 * @e: log2 of n
 * @degRes: degree of result
 * @rootsPtr: powers of the primitive n-th root of unity
 * @degA: (Input and Output) degree of A
 * @APtr: coefficient vector of A 
 * @degB: degree of B
 * @BPtr: coefficient vector of B
 * @pPtr: prime number structure
 * 
 * FFT-based multiplication of A by B. In place: result in A.
 * note: for Kronecker.
 * 
 * Return value: the product of two polynomials.
 **/   
  void EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr){
    //register int i;
    sfixn l_AB, e1, n1, sA, sB;

    l_AB = degRes + 1; // compute l_AB for FFT(A,B).
    e1 = logceiling(l_AB);
	checkFrourierPrime(e1, pPtr);
	
    n1 = 1L << e1;
    
    // force this FFT use the n from parameter if n<n1.
    if ( (n<n1) || (degRes==0) ) { n1=n; e1=e;}

	sA = degA+1;
	sB =degB+1;
    cilk_for(int i=sA; i<n1; i++) APtr[i]=0;
    cilk_for(int i=sB; i<n1; i++) BPtr[i]=0;  
    
    cilk_spawn EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, APtr, pPtr);
    EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, BPtr, pPtr);
    cilk_sync;
    
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n1, APtr, BPtr, pPtr);
    
    EX_Mont_INVDFT_OPT2_AS_GENE_R_1(n1, e1, rootsPtr, APtr, pPtr);    
  }

  //-------------------------------------
  //contract from 4D (d1, d2, d3, d4) to 2D:
  //deg1=d1; deg2=d2+(d2+1)*d3+(d2+(d2+1)*d3+1)*d4;
  //deg2=d2 + (d2+1)*d3 + (d2+1)*(d3+1)*d4;
  void MultiplyByTFT_3D2(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
	
#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif
	
#ifdef TIMEFN
	start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
#ifdef TIMEFN  
	end = example_get_time();
	initkn_time += (end - start) / 1.f;
#endif

	sfixn deg1 = BUSZSI(rep12,1);
	sfixn deg2 = BUSZSI(rep12,2) + (BUSZSI(rep12,2)+1)*BUSZSI(rep12,3) + (BUSZSI(rep12,2)+1)*(BUSZSI(rep12,3)+1)*BUSZSI(rep12,4);

	sfixn * ES_2D = (sfixn * )my_calloc(3, sizeof(sfixn));
	sfixn * DIMS_2D=(sfixn * )my_calloc(3, sizeof(sfixn));
	sfixn * LS_2D=(sfixn *)my_calloc(3, sizeof(sfixn));
	
	ES_2D[1] = logceiling(deg1+1);
	ES_2D[2] = logceiling(deg2+1);
	DIMS_2D[1] = 1 << ES_2D[1];
	DIMS_2D[2] = 1 << ES_2D[2];
	LS_2D[1] = deg1+1;
	LS_2D[2] = deg2+1;
	
	//---------------------------
	PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), 2, ES_2D, DIMS_2D, LS_2D, pPtr);
	//PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	
	my_free(ES_2D);
	my_free(DIMS_2D);
	my_free(LS_2D);

#ifdef TIMEFN
	  start = example_get_time();
#endif					
    PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }


  //-------------------------------------
  //contract from 4D (d1, d2, d3, d4) to 2D:
  //deg1=d1+(d1+1)*d2; deg2=d3+(d3+1)*d4;
  void MultiplyByTFT_2D2(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif
	
#ifdef TIMEFN
	start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
#ifdef TIMEFN  
	end = example_get_time();
	initkn_time += (end - start) / 1.f;
#endif

	sfixn deg1 = BUSZSI(rep12, 1) + (BUSZSI(rep12, 1) + 1)*BUSZSI(rep12, 2);
	sfixn deg2 = BUSZSI(rep12, 3) + (BUSZSI(rep12, 3) + 1)*BUSZSI(rep12, 4);

	sfixn * ES_2D = (sfixn * )my_calloc(3, sizeof(sfixn));
	sfixn * DIMS_2D=(sfixn * )my_calloc(3, sizeof(sfixn));
	sfixn * LS_2D=(sfixn *)my_calloc(3, sizeof(sfixn));
	
	ES_2D[1] = logceiling(deg1+1);
	ES_2D[2] = logceiling(deg2+1);
	DIMS_2D[1] = 1 << ES_2D[1];
	DIMS_2D[2] = 1 << ES_2D[2];
	LS_2D[1] = deg1+1;
	LS_2D[2] = deg2+1;
	
	//---------------------------
	PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), 2, ES_2D, DIMS_2D, LS_2D, pPtr);
	//PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	
	my_free(ES_2D);
	my_free(DIMS_2D);
	my_free(LS_2D);

#ifdef TIMEFN
	  start = example_get_time();
#endif					
    PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }

  //bug!!!!! in mixing extension and contraction
  //-------------------------------------
  //contract from ND (d1, d2, ..., dn) to 2D:
  //trial version: dm is splitted into two parts
  void MultiplyByTFT_RBB(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif
	
#ifdef TIMEFN
	start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
#ifdef TIMEFN  
	end = example_get_time();
	initkn_time += (end - start) / 1.f;
#endif

	sfixn * factm = (sfixn * )my_calloc(3, sizeof(sfixn));
	PBPAS::ReductionToBalancedBivar(N, LS(kPtr), factm); 

	sfixn deg1 = 0;
	sfixn deg2 = 0;

	if ( factm[1]==1 || factm[2]==1 ) {
	  my_free(factm);
	  PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	}else{
	  if (factm[0] == 1) {
		deg1 = factm[1]-1;
		
		deg2 = factm[2]-1;
		if (N>1) {
		  sfixn ls = factm[2];
		  deg2 += ls*BUSZSI(rep12, factm[0]+1);
		  for (int i=factm[0]+2; i<=N; i++) {
			ls *= LSI(kPtr, i-1);
			deg2 += ls*BUSZSI(rep12, i);
		  }
		}
		my_free(factm);
	  }else{
		if (factm[0] == N) {
		  deg1 = BUSZSI(rep12, 1);
		  sfixn ls = 1;
		  for (int i=2; i<N; i++) {
			ls *= LSI(kPtr, i-1);
			deg1 += ls*BUSZSI(rep12, i);
		  }
		  ls *= factm[1];
		  deg1 += ls*(factm[1]-1);
		  
		  deg2 = factm[2]-1;
		  my_free(factm);
		}else{
		  deg1 = BUSZSI(rep12, 1);
		  sfixn ls = 1;
		  for (int i=2; i<factm[0]; i++) {
			ls *= LSI(kPtr, i-1);
			deg1 += ls*BUSZSI(rep12, i);
		  }
		  ls *= factm[1];
		  deg1 += ls*(factm[1]-1);
		  
		  deg2 = factm[2]-1;
		  ls = factm[2];
		  deg2 += ls*BUSZSI(rep12, factm[0]+1);
		  for (int i=factm[0]+2; i<=N; i++) {
			ls *= LSI(kPtr, i-1);
			deg2 += ls*BUSZSI(rep12, i);
		  }
		  my_free(factm);
		}
	  }
	  
	  std::cout<<deg1<<" "<<deg2<<std::endl;

	  sfixn * ES_2D = (sfixn * )my_calloc(3, sizeof(sfixn));
	  sfixn * DIMS_2D=(sfixn * )my_calloc(3, sizeof(sfixn));
	  sfixn * LS_2D=(sfixn *)my_calloc(3, sizeof(sfixn));
	  
	  ES_2D[1] = logceiling(deg1+1);
	  ES_2D[2] = logceiling(deg2+1);
	  DIMS_2D[1] = 1 << ES_2D[1];
	  DIMS_2D[2] = 1 << ES_2D[2];
	  LS_2D[1] = deg1+1;
	  LS_2D[2] = deg2+1;
	  
	  //---------------------------
	  PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), 2, ES_2D, DIMS_2D, LS_2D, pPtr);
	  //PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	  
	  my_free(ES_2D);
	  my_free(DIMS_2D);
	  my_free(LS_2D);
	}

#ifdef TIMEFN
	  start = example_get_time();
#endif					
    PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }

  /*
    c = a*b mod p
  */
  void AdaptiveBivarMul(sfixn *c, int *dc, 
			sfixn *a, int *da, 
			sfixn *b, int *db,
			sfixn p)
  {//---------------------------------------------
    //test using fft first
    int ls1 = dc[0]+1;
    int es1 = logceiling(ls1);
    int ls2 = dc[1]+1;
    int es2 = logceiling(ls2);
    
  }

  //-------------------------------------
  // simple contract from ND (d1, d2, ..., dn) to 2D:
  // no extension
  // no copy to init KN construct, and saved TFT in dimension-1
  // avoid TFT work on zero vectors for the evaluation of 1st variable
  // for 1D, copy directly from rep1 to a vector of size dims1 for TFT and
  // put result to rep12, then for dimension-2 TFT 
  void MultiplyByTFT_RBBnoE_saveOnDim1(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

	if (N==2) {
	  bivarMultiplyBy2DTFT(2, rep12, rep1, rep2, pPtr);
	} else {
	   
	  sfixn * ls = (sfixn * )my_calloc(N+1,sizeof(sfixn));
	  for(int j=1; j<=N; j++){
		ls[j] = BUSZS(rep12)[j] +1;
	  }
	  sfixn t = PBPAS::noExtensionBalancedBivar(N, ls);
	  sfixn ls1 = 1;
	  sfixn ls2 = 1;	  

	  for (int i=1; i<=t; i++) {
		ls1 *= ls[i];
	  }

	  for (int i=t+1; i<=N; i++) {
		ls2 *= ls[i];
	  }

	  //std::cout<<ls1-1<<" "<<ls2-1<<" ";
	  
	  sfixn es1 = logceiling(ls1);
	  sfixn es2 = logceiling(ls2);

	  checkFrourierPrime(es1, pPtr);
	  checkFrourierPrime(es2, pPtr);

	  sfixn dim1 = 1<<es1;
	  sfixn dim2 = 1<<es2;
	  //sfixn n = ls1*ls2;
	  
	  sfixn *rootsPtr=(sfixn *) my_calloc(dim1+dim2, sizeof(sfixn));
	  //--------------------------------------
	  cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es1, dim1, rootsPtr, pPtr);
	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dim2, rootsPtr+dim1, pPtr);
	  cilk_sync;

	  cilk_spawn 
		PBPAS::MultiEvalByTFT_RBBnoE_saveOnDim1(N, t, DAT(rep12), CUM(rep12), ls1, ls2, dim1, dim2, DAT(rep1), BUSZS(rep1), CUM(rep1), rootsPtr, pPtr);
	  
	  //for the 2nd tft transform in the evaluation of rep2
	  sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	  PBPAS::MultiEvalByTFT_RBBnoE_saveOnDim1(N, t, tftData2, CUM(rep12), ls1, ls2, dim1, dim2, DAT(rep2), BUSZS(rep2), CUM(rep2), rootsPtr, pPtr);
	  cilk_sync;

	  //PBPAS::EX_Mont_PairwiseMul_OPT2_AS(SIZ(rep12), DAT(rep12), tftData2, pPtr->P); //diff by R?

	  PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(SIZ(rep12), DAT(rep12), tftData2, pPtr);
	  my_free(tftData2);

	  //now data is in ls2xls1 
	  PBPAS::bivarInterpolBy2DTFT(DAT(rep12), dim1, dim2, ls1, ls2, rootsPtr, pPtr);
	}
  }


  //------------------------------------
  //t1: index where partial degrees of the result rpdegs split into 2 groups
  //ls1=(rpdegs[1]+1)*...*(rpdegs[t1-1]+1), dims1
  //ls2=(rpdegs[t1]+1)*...*(rpdegs[N]+1), dims2
  //rccum: accumulated size in the product res
  //ccum: accumulated size in a given poly coeffs
  //dgs: partial degree vector of coeffs
  void MultiEvalByTFT_RBBnoE_saveOnDim1(sfixn N, sfixn t, sfixn *resCoeffs, sfixn *rccum, sfixn ls1, sfixn ls2, sfixn dim1, sfixn dim2, sfixn *coeffs, sfixn *dgs, sfixn *ccum, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	
	if (ls1) {
	  evalDim1(N, t, resCoeffs, rccum, ls1, dim1, coeffs, dgs, ccum, rootsPtr, pPtr);
	}
	
	if (ls2) {

	  if ( ls1==ls2 )
		PBPAS::sqtranspose(resCoeffs, ls2, 0, ls1, 0, 1, ls1);
	  else{	
		sfixn n = ls1*ls2;
		sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
		PBPAS::transpose(resCoeffs, ls1, B, ls2, 0, ls2, 0, ls1);
		memcpy(resCoeffs, B, n*sizeof(sfixn));
		//cilk_for(int j = 0; j < n; j++) resCoeffs[j] = B[j];
		my_free(B);
	  }

	  sfixn * tmprootsPtr=rootsPtr+dim1;
	  size_t s = ls2*sizeof(sfixn);
	  
	  //#pragma cilk_grainsize = 1;
	  cilk_for(int j=0; j<ls1; j++){
		sfixn * tmpVec = (sfixn *)my_calloc(dim2, sizeof(sfixn));
		sfixn * pos = resCoeffs+j*ls2;
		memcpy(tmpVec, pos, s);
		EX_Mont_TDFT_OPT2_AS_GENE_1( ls2, tmprootsPtr, tmpVec, pPtr);
		memcpy(pos, tmpVec, s);
		
		my_free(tmpVec);
	  }
	  //now resCoeffs' layout is ls2 x ls1
	}
  }

  //--------------------------------------
  void evalDim1(sfixn N, sfixn t, sfixn *resCoeffs, sfixn *rccum, sfixn ls1, sfixn dim1, sfixn *coeffs, sfixn *dgs, sfixn *ccum, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	sfixn d = dgs[N];
	if (N==t) {
	  sfixn * tmpVec = (sfixn *)my_calloc(dim1, sizeof(sfixn));

	  //move the data to the right slots according to the degree
	  //pattern of the result 
	  fromtofftRep(N, rccum, tmpVec, ccum, dgs, coeffs);

	  EX_Mont_TDFT_OPT2_AS_GENE_1( ls1, rootsPtr, tmpVec, pPtr);
	  memcpy(resCoeffs, tmpVec, ls1*sizeof(sfixn));
	  my_free(tmpVec);
	  
	}else{
	  //#pragma cilk_grainsize = 1;
	  cilk_for(int i=0; i<=d; i++){
		PBPAS::evalDim1(N-1, t, resCoeffs+i*rccum[N], rccum, ls1, dim1, coeffs+i*ccum[N], dgs, ccum, rootsPtr, pPtr);
	  }
	}
  }

  //-------------------------------------
  //simple contract from ND (d1, d2, ..., dn) to 2D:
  //no extension
  void MultiplyByTFT_RBBnoE(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

	if ( N==2 ) {
	  bivarMultiplyBy2DTFT(N, rep12, rep1, rep2, pPtr);
		//PBPAS::tftMultiD_test_1_2DTran(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);

	}else{
// 	  std::cout<< "---rep1---"<<std::endl;
// 	  for (int i=0; i<SIZ(rep1); i++) {
// 		std::cout<< DAT(rep1)[i] <<", ";
// 	  }
// 	  std::cout<< "---end of rep1---"<<std::endl;

// 	  std::cout<< "---rep2---"<<std::endl;
// 	  for (int i=0; i<SIZ(rep2); i++) {
// 		std::cout<< DAT(rep2)[i] <<", ";
// 	  }
// 	  std::cout<< "---end of rep2---"<<std::endl;
 
#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif
	
#ifdef TIMEFN
	start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    cilk_spawn PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
	cilk_sync;
	
#ifdef TIMEFN  
	end = example_get_time();
	initkn_time += (end - start) / 1.f;
#endif

	  sfixn t = PBPAS::noExtensionBalancedBivar(N, LS(kPtr)); 

	  sfixn ls1 = 1;
	  sfixn ls2 = 1;	  

	  for (int i=1; i<=t; i++) {
		ls1 *= LSI(kPtr, i);
	  }

	  for (int i=t+1; i<=N; i++) {
		ls2 *= LSI(kPtr, i);
	  }

	  std::cout<<ls1-1<<" "<<ls2-1<<" ";

	  sfixn * ES_2D = (sfixn * )my_calloc(3, sizeof(sfixn));
	  sfixn * DIMS_2D=(sfixn * )my_calloc(3, sizeof(sfixn));
	  sfixn * LS_2D=(sfixn *)my_calloc(3, sizeof(sfixn));
	  
	  ES_2D[1] = logceiling(ls1);
	  ES_2D[2] = logceiling(ls2);
	  DIMS_2D[1] = 1 << ES_2D[1];
	  DIMS_2D[2] = 1 << ES_2D[2];
	  LS_2D[1] = ls1;
	  LS_2D[2] = ls2;
	  
	  //---------------------------
	  //2 dimensional
	  PBPAS::tftMultiD_test_1_2DTran(DATSI(kPtr, 0), DATSI(kPtr, 1), 2, ES_2D, DIMS_2D, LS_2D, pPtr);
	  //PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	  
	  my_free(ES_2D);
	  my_free(DIMS_2D);
	  my_free(LS_2D);

#ifdef TIMEFN
	  start = example_get_time();
#endif					
    PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif

	}
  }

  //-------------------------------------
  //contract from ND (d1, d2, ..., dn) to 2D by:
  //first convert to 1D by Kronecker
  //then extend to 2D
  //no temprary storage of data
  //not used, slower than multi-D
  void MultiplyByTFT_RBB_KN1D_to_2D(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    KroTFTRep* kPtr1 = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr1, BUSZS(rep12), N, 2, pPtr);

    cilk_spawn PBPAS::fromtofftRepMultiD(N,  CUM(kPtr1), DATSI(kPtr1, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr1), DATSI(kPtr1, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
	cilk_sync;

	sfixn deg1d_1 = BUSZSI(rep1, 1);
		for (int i=2; i<=N; i++) {
		  deg1d_1 += CUMI(kPtr1, i)*BUSZSI(rep1, i);
	}

	sfixn deg1d_2 = BUSZSI(rep2, 1);
		for (int i=2; i<=N; i++) {
	  deg1d_2 += CUMI(kPtr1, i)*BUSZSI(rep2, i);
	}	  

		std::cout<<deg1d_1<<" "<<deg1d_2<<std::endl;

	  sfixn b = PBPAS::V1to2V_pivot(deg1d_1, deg1d_2);

	  KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
	  PBPAS::InitKroTFTRep_1V2V(kPtr, b, deg1d_1, deg1d_2);
	  //extend 1V to 2V and store data for TFT
	  PBPAS::from1Vto2VTFTRep(b, deg1d_1, DATSI(kPtr1, 0), CUMI(kPtr, 2), DATSI(kPtr, 0));
	  PBPAS::from1Vto2VTFTRep(b, deg1d_2, DATSI(kPtr1, 1), CUMI(kPtr, 2), DATSI(kPtr, 1));
	  PBPAS::freeKroTFTRep(kPtr1);

#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif

	  //---------------------------
	  // 2 variables
	  PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), 2, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	  
#ifdef TIMEFN
	  start = example_get_time();
#endif

	  PBPAS::from2VTFTRepto1V(b, LSI(kPtr, 1), LSI(kPtr, 2), DATSI(kPtr, 0), SIZ(rep12), DAT(rep12), pPtr);

	  //std::cout<<"done from1Vto2VtoTFTRep"<<std::endl;
	  
	  PBPAS::freeKroTFTRep(kPtr);
	  //std::cout<<"done freeKroTFTRep"<<std::endl;
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif

  }


  //-------------------------------------
  //rep1 and rep2 are univariate
  //extend to bivariate and multiply by 2D-TFT
  //no KroTFTRep
  void MultiplyByTFT_1V2V(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif

	sfixn b = PBPAS::V1to2V_pivot(BUSZSI(rep1, 1), BUSZSI(rep2, 1));
	sfixn ls1 = 2*b-1;
	sfixn ls2 = ( BUSZSI(rep1,1) -(BUSZSI(rep1,1))%b )/b + ( BUSZSI(rep2,1) -(BUSZSI(rep2,1))%b )/b +1;

	//std::cout<<"b, ls1, ls2= "<<b<<", "<<ls1<<", "<<ls2<<std::endl;

	sfixn es1 = logceiling(ls1);
	sfixn es2 = logceiling(ls2);
	sfixn dims1 = 1<<es1;
	sfixn dims2 = 1<<es2;

#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif

	  checkFrourierPrime(es1, pPtr);
	  checkFrourierPrime(es2, pPtr);

	  sfixn *rootsPtr=(sfixn *) my_calloc(dims1 +dims2, sizeof(sfixn));

	  cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es1, dims1, rootsPtr, pPtr);
	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dims2, rootsPtr+dims1, pPtr);
	  cilk_sync;

	  sfixn siz=ls1*ls2;

	  sfixn * tftData1=(sfixn * )my_calloc(siz, sizeof(sfixn));
	  sfixn * tftData2=(sfixn * )my_calloc(siz, sizeof(sfixn));

	  cilk_spawn PBPAS::V1to2VEvalBy2DTFT(tftData1, dims1, dims2, ls1, ls2, DAT(rep1), BUSZSI(rep1, 1), b, rootsPtr, pPtr);
	  
	  PBPAS::V1to2VEvalBy2DTFT(tftData2, dims1, dims2, ls1, ls2, DAT(rep2), BUSZSI(rep2, 1), b, rootsPtr, pPtr);
	  cilk_sync;

	  PBPAS::EX_Mont_PairwiseMul_OPT2_AS(siz, tftData1, tftData2, pPtr->P);

	  my_free(tftData2);

	  //input DATSI(kPtr,0) in ls2 x ls1 layout
	  PBPAS::bivarInterpolBy2DTFT(tftData1, dims1, dims2, ls1, ls2, rootsPtr, pPtr);
	  
#ifdef TIMEFN
	  start = example_get_time();
#endif

	  PBPAS::from2VTFTRepto1V(b, ls1, ls2, tftData1, SIZ(rep12), DAT(rep12), pPtr);

	  my_free(tftData1);
	  
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif

  }

  // work directly from input poly coeffs array for dimension-1
  // evaluation, put result into the right position in resCoeffs
  // b: pivot for extending v1 to 2v
  // d: degree of univariate poly coeffs
  void V1to2VEvalBy2DTFT(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn d, sfixn b, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	sfixn d2 = (d-d%b)/b; //d quo b, partial degree of var2 after converting
						  //poly coeffs1 from univar to bivar

	size_t s = b*sizeof(sfixn);
	size_t res = ls1*sizeof(sfixn);

	//#pragma cilk_grainsize = 1;
	cilk_for(int j=0; j<d2; j++){
	  //evaluate 1st var for the num of d2 coefficients
	  //saved by not tft on ls1*ls2/d2 number of vectors of zeros
	  //for parallel run
	  sfixn * tmpVec = (sfixn *)my_calloc(dims1, sizeof(sfixn));
	  memcpy(tmpVec, coeffs+j*b, s);
	  EX_Mont_TDFT_OPT2_AS_GENE_1( ls1, rootsPtr, tmpVec, pPtr);
	  memcpy(resCoeffs+j*ls1, tmpVec, res);

	  my_free(tmpVec);
	}

	sfixn tail = d+1 - d2*b;
	if ( tail > 0 ) {
	  sfixn * tmpVec = (sfixn *)my_calloc(dims1, sizeof(sfixn));
	  //std::cout<<"before memcpy tail"<<std::endl;
	  memcpy(tmpVec, coeffs+d2*b, tail*sizeof(sfixn));
	  //std::cout<<"after memcpy tail"<<std::endl;
	  EX_Mont_TDFT_OPT2_AS_GENE_1( ls1, rootsPtr, tmpVec, pPtr);
	  memcpy(resCoeffs+d2*ls1, tmpVec, res);

	  my_free(tmpVec);
	}

	if ( ls1==ls2 )
	  PBPAS::sqtranspose(resCoeffs, ls2, 0, ls1, 0, 1, ls1);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
	  PBPAS::transpose(resCoeffs, ls1, B, ls2, 0, ls2, 0, ls1);
	  memcpy(resCoeffs, B, n*sizeof(sfixn));
	  //cilk_for(int j = 0; j < n; j++) resCoeffs[j] = B[j];
	  my_free(B);
	}
	
	sfixn * tmprootsPtr=rootsPtr+dims1;
	s = ls2*sizeof(sfixn);

	//#pragma cilk_grainsize = 1;
	cilk_for(int j=0; j<ls1; j++){
	  sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn));
	  memcpy(tmpVec, resCoeffs+j*ls2, s);
	  EX_Mont_TDFT_OPT2_AS_GENE_1( ls2, tmprootsPtr, tmpVec, pPtr);
	  memcpy(resCoeffs+j*ls2, tmpVec, s);
	  
	  my_free(tmpVec);
	}
	//std::cout<<"end of eval"<<std::endl;
	//now resCoeffs' layout is ls2 x ls1
	
  }

  //----------------
  //contract d1,...dn to 2 vars
  //d1 is split into 2 parts
  void MultiplyByTFT_1Vto2V_multiV_2V(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr) {
	
#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif

	  sfixn sigma = 1;
	  for (int i=2; i<=N; i++)
		sigma *= BUSZSI(rep1, i) + BUSZSI(rep2, i) +1;

	  sfixn b = PBPAS::multi_V1to2V_pivot_opt(sigma, BUSZSI(rep1, 1), BUSZSI(rep2, 1));
	  //std::cout <<"sigma, b= "<<sigma<<", "<<b <<std::endl;
	  if (b==0) {
		std::cout <<"Exit 1, b is zero"<<std::endl;
		exit(1);
	  }
	  sfixn q1rep1 = (BUSZSI(rep1, 1) - BUSZSI(rep1, 1) % b)/b;
  	  sfixn q1rep2 = (BUSZSI(rep2, 1) - BUSZSI(rep2, 1) % b)/b;
	  
	  KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
 	  PBPAS::InitKroTFTRep_1V2V_multi(kPtr, b, sigma, q1rep1, q1rep2);

	  //extend 1V to 2V and store multi var data in 2 var format for
	  //2D-TFT
	  cilk_spawn PBPAS::fromMulti_1Vto2VTFTRep(N, b, q1rep1, q1rep2, BUSZS(rep1), DAT(rep1), BUSZS(rep2), DATSI(kPtr, 0));
	  PBPAS::fromMulti_1Vto2VTFTRep(N, b, q1rep2, q1rep1, BUSZS(rep2), DAT(rep2), BUSZS(rep1), DATSI(kPtr, 1));
	  cilk_sync;

	  //sfixn r1rep1 = BUSZSI(rep1, 1) % b;
	  //sfixn r1rep2 = BUSZSI(rep2, 1) % b;

 	  //extend 1V to 2V and store data for TFT
 	  //PBPAS::fromMulti_1Vto2VTFTRep_1(N, b, q1rep1, BUSZS(rep1), DAT(rep1), sigma, q1rep1+q1rep2+1, DATSI(kPtr, 0));
	  //PBPAS::fromMulti_1Vto2VTFTRep_1(N, b, q1rep2, BUSZS(rep2), DAT(rep2), sigma, q1rep1+q1rep2+1, DATSI(kPtr, 1));	  

  #ifdef TIMEFN  
  	  end = example_get_time();
 	  initkn_time += (end - start) / 1.f;
  #endif

  	  //=====================================
  	  // 2 variables
  	  PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), 2, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
  	  //=====================================
	  
	  //size of DATSI(kPtr, 0) can be larger (twice) than size of DAT(rep12)
	  PBPAS::from2VTFTReptoMultiV(b, sigma, q1rep1+q1rep2+1, DATSI(kPtr, 0), BUSZSI(rep1,1)+BUSZSI(rep2,1)+1, DAT(rep12), pPtr);

 	  PBPAS::freeKroTFTRep(kPtr);

#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif
	  

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif	  
	  
  }

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  //work-in-progress
  //----------------------
  //do not use KroTFTRep
  //make use of rep12 for the evaluation of the first poly
  //save from not doing tft on zeros in 1D-tft on the first variable
  //use size of dims[i] for copy out tft in dimension i 
  //use Matteo's 2D matrix transposition for multi-D matrix transposition
  void MultiplyByTFT_2DTran_noKroRep(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

	sfixn * es = (sfixn * )my_calloc(N+1,sizeof(sfixn));
	sfixn * dims = (sfixn * )my_calloc(N+1,sizeof(sfixn));
	sfixn * ls = (sfixn * )my_calloc(N+1,sizeof(sfixn));
	
	for(int j=1; j<=N; j++){
	  ls[j] = BUSZS(rep12)[j] +1;
      es[j] = logceiling( ls[j] );
      dims[j] = 1<< es[j];      
	}
	//check Frourier prime
	for(int i=1; i<=N; i++) {
	  checkFrourierPrime(es[i], pPtr);
	}

	sfixn m=0;
	for(int i=1; i<=N; i++) {//N is usually small
      m += dims[i]; 
    }
	sfixn * rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
    sfixn * tmprootsPtr1 = rootsPtr;
	for(int i=1; i<=N; i++) {
	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr1, pPtr);
      tmprootsPtr1 += dims[i];
    }

#ifdef TIMEFN
	//------------------------------------
	//run two evaluations in serial to measure the performance of 2D-TFT
	long start = example_get_time();	
	PBPAS::MultiEvalByTFT_2DTran_saveOn1stVar(DAT(rep12), N, SIZ(rep12), dims, ls, CUM(rep12), DAT(rep1), SIZ(rep1), BUSZS(rep1), CUM(rep1), rootsPtr, pPtr);
	long end = example_get_time();
	float oneEvaltime = (end - start) / 1.f;
	std::cout<< oneEvaltime << " ";

	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	sfixn * tmprootsPtr;
	tmprootsPtr = PBPAS::MultiEvalByTFT_2DTran_saveOn1stVar(tftData2, N, SIZ(rep12), dims, ls, CUM(rep12), DAT(rep2), SIZ(rep2), BUSZS(rep2), CUM(rep2), rootsPtr, pPtr);

#else	
	//evaluation-----------------------------
	cilk_spawn
	PBPAS::MultiEvalByTFT_2DTran_saveOn1stVar(DAT(rep12), N, SIZ(rep12), dims, ls, CUM(rep12), DAT(rep1), SIZ(rep1), BUSZS(rep1), CUM(rep1), rootsPtr, pPtr);

	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	sfixn * tmprootsPtr;
	tmprootsPtr = PBPAS::MultiEvalByTFT_2DTran_saveOn1stVar(tftData2, N, SIZ(rep12), dims, ls, CUM(rep12), DAT(rep2), SIZ(rep2), BUSZS(rep2), CUM(rep2), rootsPtr, pPtr);

	cilk_sync;
#endif	
	//pairwise mul----------------------------
	PBPAS::EX_Mont_PairwiseMul_OPT2_AS(SIZ(rep12), DAT(rep12), tftData2, pPtr->P);

	my_free(tftData2);

#ifdef TIMEFN
float interpol_time = 0.0;
	long start1, end1;
	start1 = example_get_time();
#endif

	PBPAS::InterpolByTFT_2DTran(DAT(rep12), N, SIZ(rep12), dims, ls, pPtr, tmprootsPtr);

#ifdef TIMEFN  
	  end1 = example_get_time();
	  interpol_time += (end1 - start1) / 1.f;
	  std::cout << interpol_time << " " ;
#endif

  }

  //------------------------------------
  //resCoeffs: array of result size [0..n-1]
  //N: number of vars. input order is x_1, x_2, ... x_N
  //n: size of resCoeffs
  //es: vector of powers [1..N]. 2^es[i]=dims[i]
  //dims: vector of dimension size [1..N]
  //ls: vector of (partial degree +1) for each var in the result
  //[1..N]. none should be zero.
  //rccum: CUM(rep12)
  //coeffs: data array of input poly
  //n1: size of coeffs
  //pdegs: partial degree vector [1..N] of input poly
  //ccum: CUM(rep1)
  //rootsPtr: primitive root (powered) vector of result for all
  //dimensions [1..N]
  sfixn * MultiEvalByTFT_2DTran_saveOn1stVar(sfixn * resCoeffs, sfixn N, sfixn n, sfixn *dims, sfixn * ls, sfixn * rccum, sfixn * coeffs, sfixn n1, sfixn * pdegs, sfixn * ccum, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){

// 	std::cout <<"enter MultiEvalByTFT_2DTran_saveOn1stVar" << std::endl;
// 	for(int j = 0; j < n1; j++) std::cout<< coeffs[j]<<", ";
// 	std::cout<<"----end of initial data----" <<std::endl;

	//for 1st var evaluation, copy out each segment of data from
	//coeffs and copy to the right position in resCoeffs after tft
	//do pdegs[2]*...*pdegs[N]
	//do not do tft for n/(pdegs[2]*...*pdegs[N]) zeros
	if (ls[1]) {
	  PBPAS::eval1stVar(N, rccum, resCoeffs, ccum, pdegs, coeffs, dims[1], ls[1], rootsPtr, pPtr);
	}
	sfixn * tmprootsPtr = rootsPtr;
	for(int i=2; i<=N; i++){

// 	  std::cout<<"----in data----" <<std::endl;
// 	  for (int j = 1; j <= N; j++) {
// 	  	std::cout<< ls[j] <<", ";
// 	  }
// 	  std::cout<<std::endl;

// 	  for(int j = 0; j < n; j++) std::cout<< resCoeffs[j]<<", ";
// 	  std::cout<<"----end of in data----" <<std::endl;

	  //------------------------------
      PBPAS::multi_mat_transpose_2DTran_eval (N, n, i, ls, resCoeffs);
	  
// 	  for(int j = 0; j < n; j++) std::cout<< resCoeffs[j]<<", ";
// 	  std::cout<<"----end of out data----" <<std::endl;	  

	  //point to the start position of ith dim's roots in rootsPtr
      tmprootsPtr += dims[i-1]; 
	  if( ls[i] ){ //each TDFT uses a temp vector of size dims[i]
		sfixn d = n/ls[i];
		//size is small so grainsize relys on cilk_for
		cilk_for(int j=0; j<d; j++){
		  PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1_par(ls[i], tmprootsPtr, dims[i], resCoeffs+j*ls[i], pPtr);
		}
	  }
	}
	//now data is in the order of x_N, x_{N-1}, ..., x_1
	//tmprootsPtr now points to the begining of dims[N]'s roots
	//return tmprootsPtr, cilk++ needs
	sfixn *tmpPtr = tmprootsPtr;
	return tmpPtr;
  }
 
  //------------------------------------
  //evaluate 1st var for the number of non-zero slots.
  //rccum: accumulated size in the product res
  //ccum: accumulated size in a given poly coeffs
  //dgs: partial degree vector of coeffs
  //dim1 and ls1: dimension (power of 2) and real size of
  //              the first var in res
  void eval1stVar(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum, sfixn * dgs, sfixn * coeffs, sfixn dim1, sfixn ls1, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){

	sfixn d = dgs[N];
	if(N==1){
	  sfixn * tmpVec = (sfixn *)my_calloc(dim1, sizeof(sfixn));
	  memcpy(tmpVec, coeffs, (d+1)*sizeof(sfixn));
	  EX_Mont_TDFT_OPT2_AS_GENE_1( ls1, rootsPtr, tmpVec, pPtr);
	  memcpy(res, tmpVec, ls1*sizeof(sfixn));
	  my_free(tmpVec);
	  //return;
	}else{
	  //size is small so grainsize relys on cilk_for
	  cilk_for(int i=0; i<=d; i++){
		PBPAS::eval1stVar(N-1, rccum, res+i*rccum[N], ccum, dgs, coeffs+i*ccum[N], dim1, ls1, rootsPtr, pPtr);
	  }
	}
  }

  //-------------------------------------
  void MultiplyByTFT_2DTran(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

//     cilk_spawn PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
//     PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
// 	cilk_sync;
	
	cilk_spawn PBPAS::fromtofftRep(N,  CUM(kPtr), DAT(rep12), CUM(rep1), BUSZS(rep1), DAT(rep1));
	PBPAS::fromtofftRep(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
	cilk_sync;

#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif

    //---------------------------
    PBPAS::tftMultiD_test_1_2DTran(DAT(rep12), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);

	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );

#ifdef TIMEFN
	  start = example_get_time();
#endif					
	  //PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }
  //================================================================
  // Parallel Multi-dimensional TFT by seperated evaluation 
  //================================================================
  /**
   * tftMultiD_test_1:
   * @coeffs1: coefficient vector for 'f1'.
   * @coeffs2: coefficient vector for 'f2'.
   * @N: number of variables.
   * @es: 2^es[i] = dims [i] for i = 1..n.
   * @dims: the FFT sizes on each dimension, dims[i] is the FFT on dimensional i, i=1..N.
   * @pPtr: the information of the prime number.
   * 
   *Notes: each TDFT uses a temp vector of size maxdim;
   *       each INVTDFT uses a temp vector of size maxdim
   *
   * Return value: 
   * do not use maxdim
   **/
  void tftMultiD_test_1_2DTran(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr){
	register int i;
	sfixn m=0, n=1;
	sfixn * tmprootsPtr;

	//std::cout <<"enter tftMultiD_test_1" << std::endl;

	//check Frourier prime
	for(i=1; i<=N; i++) {
	  checkFrourierPrime(es[i], pPtr);
	}

	for(i=1; i<=N; i++) {//N is usually small
	  n *= ls[i]; 
      m += dims[i]; 
    }

	sfixn * rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
    sfixn * tmprootsPtr1 = rootsPtr;
	//compute the primitive roots for all dimensions
	for(i=1; i<=N; i++) {
	  PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr1, pPtr);
      tmprootsPtr1 += dims[i];
    }

	//std::cout <<"before MultiEvalByTFT" << std::endl;

#ifdef TIMEFN
	float pairwisemul_time = 0.0;
	float pareval_time = 0.0;
	float interpol_time = 0.0;
	long start1, end1;
	start1 = example_get_time();
#endif

#ifdef SERIALFT
	tmprootsPtr = MultiEvalByTFT_2DTran(coeffs1, N, n, dims, ls, pPtr, rootsPtr);
	MultiEvalByTFT_2DTran(coeffs2, N, n, dims, ls, pPtr, rootsPtr);
#else
	//parallel run two
//     cilk_spawn MultiEvalByTFT_2DTran(coeffs1, N, n, dims, ls, pPtr, rootsPtr);

// 	tmprootsPtr = MultiEvalByTFT_2DTran(coeffs2, N, n, dims, ls, pPtr, rootsPtr);	
// 	cilk_sync;

	tmprootsPtr = MultiEvalByTFT_2DTran(coeffs1, N, n, dims, ls, pPtr, rootsPtr);
	MultiEvalByTFT_2DTran(coeffs2, N, n, dims, ls, pPtr, rootsPtr);

#endif

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pareval_time += (end1 - start1) / 1.f;
#endif

	//std::cout <<"after MultiEvalByTFT" << std::endl;

#ifdef TIMEFN
	start1 = example_get_time();
#endif
	// Pairwise-Mul, result is in coeffs1
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs2, pPtr->P);

#ifdef TIMEFN  
	  end1 = example_get_time();
	  pairwisemul_time += (end1 - start1) / 1.f;
#endif
	//my_free(coeffs2); //???????

#ifdef TIMEFN
	start1 = example_get_time();
#endif
	//is tmprootsPtr1 pointing to the last dimension??? 
    PBPAS::InterpolByTFT_2DTran(coeffs1, N, n, dims, ls, pPtr, tmprootsPtr);
#ifdef TIMEFN  
	  end1 = example_get_time();
	  interpol_time += (end1 - start1) / 1.f;
#endif
	//	std::cout <<"before exit tftMultiD_test_1" << std::endl;
    my_free(rootsPtr);

#ifdef TIMEFN
	std::cout << pareval_time << " " << pairwisemul_time << " "<< interpol_time << " " ;
#endif
	//std::cout <<"after free before exit tftMultiD_test_1" << std::endl;
  }

  /**-------------------------------------------------------
   *
   **/
  sfixn * MultiEvalByTFT_2DTran(sfixn * coeffs1, sfixn N, sfixn n, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr) {

// 	std::cout <<"enter MultiEvalByTFT_2DTran" << std::endl;
// 	for(int j = 0; j < n; j++) std::cout<< coeffs1[j]<<", ";
// 	std::cout<<"----end of initial data----" <<std::endl;

	sfixn d;

#ifdef TIMEFN
	float tran_time = 0.0;
	float tft_time = 0.0;
	long start, end;
#endif

	sfixn * tmprootsPtr = rootsPtr;

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
	if( ls[1] ) { 
	  d = n/ls[1];
	  cilk_for(int j=0; j<d; j++){
		PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1_par(ls[1], tmprootsPtr, dims[1], coeffs1+j*ls[1], pPtr);
	  }
	}
#ifdef TIMEFN  
	  end = example_get_time();
	  tft_time += (end - start) / 1.f;
#endif
		
	//std::cout <<"after ls[1]" << std::endl;

	for(int i=2; i<=N; i++){
	  
	  //point to the start position of ith dim's roots in rootsPtr
      tmprootsPtr += dims[i-1]; 

#ifdef TIMEFN
	  start = example_get_time();
#endif	

// 	  std::cout<<"----in data----" <<std::endl;
// 	  for (int j = 1; j <= N; j++) {
// 	  	std::cout<< ls[j] <<", ";
// 	  }
// 	  std::cout<<std::endl;

// 	  for(int j = 0; j < n; j++) std::cout<< coeffs1[j]<<", ";
// 	  std::cout<<"----end of in data----" <<std::endl;
	  
	  //------------------------------
      PBPAS::multi_mat_transpose_2DTran_eval (N, n, i, ls, coeffs1);
	  
// 	  for(int j = 0; j < n; j++) std::cout<< coeffs1[j]<<", ";
// 	  std::cout<<"----end of out data----" <<std::endl;

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif
	  
	  //std::cout <<"after multi_mat_transpose" << std::endl;

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
	  if( ls[i] ){ //each TDFT uses a temp vector of size dims[i]

		d = n/ls[i];
		cilk_for(int j=0; j<d; j++){
		  PBPAS::EX_Mont_TDFT_OPT2_AS_GENE_1_par(ls[i], tmprootsPtr, dims[i], coeffs1+j*ls[i], pPtr);
		}
	  }
#ifdef TIMEFN  
	  end = example_get_time();
	  tft_time += (end - start) / 1.f;
#endif

	}
#ifdef TIMEFN  
	std::cout << tran_time << " "<<tft_time << " ";
    //std::cout << "MultiEvalByTFT Mat Tran took " << tran_time << " ms." << std::endl;
#endif
	//now data is in the order of x_N, x_{N-1}, ..., x_1
	//tmprootsPtr now points to the begining of dims[N]'s roots
	//return tmprootsPtr, cilk++ needs
	//return tmprootsPtr;
	sfixn *tmpPtr = tmprootsPtr;
	return tmpPtr;
  }
  
  /**---------------------------------------------------
   *InterpolByTFT: inverse DFT
   *input data is in the order of x_N, x_{N-1}, ..., x_1
   **/
  void InterpolByTFT_2DTran(sfixn *coeffs1, sfixn N, sfixn n, sfixn * dims, sfixn * ls, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr) {

    sfixn d;

#ifdef TIMEFN
	float tran_time = 0.0;
	float invtft_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
	if( ls[N] ){
	  d = n/ls[N];
	  //size is small so grainsize relys on cilk_for
	  cilk_for(int j=0; j<d; j++){
		PBPAS::EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(coeffs1+j*ls[N], dims[N], ls[N], tmprootsPtr, pPtr);
	  }
	}
#ifdef TIMEFN  
	  end = example_get_time();
	  invtft_time += (end - start) / 1.f;
#endif

	for(int i=N-1; i>=1; i--){

      tmprootsPtr -= dims[i];

#ifdef TIMEFN
	  start = example_get_time();
#endif
	  
	  PBPAS::multi_mat_transpose_2DTran_interp (N, n, i, ls, coeffs1);

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
	  if( ls[i] ){	
		d = n/ls[i];
		//size is small so grainsize relys on cilk_for
		cilk_for(int j=0; j<d; j++){//each INVTDFT uses a temp vector of size maxdim
		  PBPAS::EX_Mont_INVTDFT_OPT2_AS_GENE_1_par(coeffs1+j*ls[i], dims[i], ls[i], tmprootsPtr, pPtr);
		}
	  }
#ifdef TIMEFN  
	  end = example_get_time();
	  invtft_time += (end - start) / 1.f;
#endif

	}
	// multi_mat_transpose_2DTran_interp does not change 
	// the order of dims and ls 
	// now data is for the order of x_1, x_2, ... x_N
#ifdef TIMEFN  
	std::cout << tran_time << " "<< invtft_time << " ";
    //std::cout << "InterpolByTFT Mat Tran took " << tran_time << " ms." << std::endl;
#endif
  }

  //newnewnewnnewnewnenwnewnenwnewnenwnen


  //-------------------------------------
  void MultiplyByTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){

#ifdef TIMEFN
	float initkn_time = 0.0;
	float reckn_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
    PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

//     cilk_spawn PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
//     PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
// 	cilk_sync;
	
	cilk_spawn PBPAS::fromtofftRep(N,  CUM(kPtr), DAT(rep12), CUM(rep1), BUSZS(rep1), DAT(rep1));
	PBPAS::fromtofftRep(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
	cilk_sync;

#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif

    //---------------------------
    PBPAS::tftMultiD_test_1(DAT(rep12), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);

	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );

#ifdef TIMEFN
	  start = example_get_time();
#endif					
	  //PBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    PBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }


  //-------------------------------------
  // for bivariate multiplication
  // N=2
  // no copy to init KN construct, and saved TFT in dimension-1
  // avoid TFT work on zero vectors for the evaluation of 1st variable
  // for 1D, copy directly from rep1 to a vector of size dims1 for TFT and
  // put result to rep12, then for dimension-2 TFT 
  // saved 30% w.r.t original 2D-TFT
  void bivarMultiplyBy2DTFT_for_spawn(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
	sfixn ls1 = BUSZSI(rep12, 1)+1;
	sfixn ls2 = BUSZSI(rep12, 2)+1;
	sfixn es1 = logceiling(ls1);
	sfixn es2 = logceiling(ls2);

	checkFrourierPrime(es1, pPtr);
	checkFrourierPrime(es2, pPtr);

	sfixn dims1 = 1<<es1;
	sfixn dims2 = 1<<es2;
	sfixn n = ls1*ls2;

// 	sfixn maxdim = dims1;
// 	if (dims1 < dims2)
// 	  maxdim = dims2;

	sfixn *rootsPtr=(sfixn *) my_calloc(dims1+dims2, sizeof(sfixn));

	//--------------------------------------
	cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es1, dims1, rootsPtr, pPtr);
	PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dims2, rootsPtr+dims1, pPtr);
	cilk_sync;
	//--------------------------------------

#ifdef TIMEFN
		//------------------------------------
	//run two evaluations in serial to measure the performance of 2D-TFT
	long start = example_get_time();	
	PBPAS::bivarEvalBy2DTFT_for_spawn(DAT(rep12), dims1, dims2, ls1, ls2, DAT(rep1), BUSZS(rep1), rootsPtr, pPtr);  
	long end = example_get_time();
	float oneEvaltime = (end - start) / 1.f;
	std::cout<< oneEvaltime << " ";
	
	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	
	PBPAS::bivarEvalBy2DTFT_for_spawn(tftData2, dims1, dims2, ls1, ls2, DAT(rep2), BUSZS(rep2), rootsPtr, pPtr); 
	//----------------------------------------------
#else
	//------------------------------------
	cilk_spawn PBPAS::bivarEvalBy2DTFT_for_spawn(DAT(rep12), dims1, dims2, ls1, ls2, DAT(rep1), BUSZS(rep1), rootsPtr, pPtr);  

	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	
	PBPAS::bivarEvalBy2DTFT_for_spawn(tftData2, dims1, dims2, ls1, ls2, DAT(rep2), BUSZS(rep2), rootsPtr, pPtr); 
	cilk_sync;
	//----------------------------------------------
#endif

	PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, DAT(rep12), tftData2, pPtr->P);

	my_free(tftData2);

#ifdef TIMEFN
	long start1 = example_get_time();
#endif
	//now data is in ls2xls1 
	PBPAS::bivarInterpolBy2DTFT_for_spawn(DAT(rep12), dims1, dims2, ls1, ls2, rootsPtr, pPtr);
#ifdef TIMEFN
	long end1 = example_get_time();
	float interp1 = (end1 - start1) / 1.f;
	std::cout<< interp1 << " ";
#endif
  }

  //--------------------------------------------
  //bivariate, ls1 and ls2 are not 0
  void bivarInterpolBy2DTFT_for_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	sfixn * tmprootsPtr = rootsPtr +dims1;
	size_t s = ls2*sizeof(sfixn);
	cilk_for(int j=0; j<ls1; j++){
	  cilk_spawn PBPAS::copyDataOutInvTFT(resCoeffs+j*ls2, dims2, s, ls2, tmprootsPtr, pPtr);
	}
	//cilk_sync;

	if ( ls1==ls2 )
	  PBPAS::sqtranspose(resCoeffs, ls1, 0, ls2, 0, 1, ls2);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(ls1*ls2, sizeof(sfixn));
	  PBPAS::transpose(resCoeffs, ls2, B, ls1, 0, ls1, 0, ls2);
#if __GNUC__ == 4
	  memcpy(resCoeffs, B, n*sizeof(sfixn));
#else
	  cilk_for(int j = 0; j < n; j++) resCoeffs[j] = B[j];
#endif
	  my_free(B);
	}	

	s = ls1*sizeof(sfixn);

	cilk_for(int j=0; j<ls2; j++){
	  cilk_spawn PBPAS::copyDataOutInvTFT(resCoeffs+j*ls1, dims1, s, ls1, rootsPtr, pPtr);
	}
	//cilk_sync;
  }

  //--------------------------------------------
  //bivariate, ls1 and ls2 are not 0
  void bivarEvalBy2DTFT_for_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn * pdegs, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	
	sfixn s1 = pdegs[1]+1;
	size_t s = s1*sizeof(sfixn);
	size_t res = ls1*sizeof(sfixn);
	cilk_for(int j=0; j<=pdegs[2]; j++){
	  //evaluate 1st var for the num of pdegs[2]+1 coefficients
	  //do not use ls2. saving from not fft on padding zeros
	  //for parallel run

	  cilk_spawn PBPAS::copyDataOutTFT(coeffs+j*s1, resCoeffs+j*ls1, dims1, s, res, ls1, rootsPtr, pPtr);
	}
	//cilk_sync;
	
	if ( ls1==ls2 )
	  PBPAS::sqtranspose(resCoeffs, ls2, 0, ls1, 0, 1, ls1);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
	  PBPAS::transpose(resCoeffs, ls1, B, ls2, 0, ls2, 0, ls1);
#if __GNUC__ == 4
	  memcpy(resCoeffs, B, n*sizeof(sfixn));
#else
	  cilk_for(int j = 0; j < n; j++) resCoeffs[j] = B[j];
#endif
	  my_free(B);
	}

	sfixn * tmprootsPtr=rootsPtr+dims1;
	s = ls2*sizeof(sfixn);
	cilk_for(int j=0; j<ls1; j++){
	  cilk_spawn PBPAS::copyDataOutTFT(resCoeffs+j*ls2, resCoeffs+j*ls2, dims2, s, s, ls2, tmprootsPtr, pPtr);

	}
	//cilk_sync;
	//now resCoeffs' layout is ls2 x ls1
  }

  //-------------------------------------
  // for bivariate multiplication
  // N=2
  // no copy to init KN construct, and saved TFT in dimension-1
  // avoid TFT work on zero vectors for the evaluation of 1st variable
  // for 1D, copy directly from rep1 to a maximum vector for TFT and
  // put result to rep12, then for dimension-2 TFT 
  // saved 30% w.r.t original 2D-TFT
  void bivarMultiplyBy2DTFT_spawn(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
	sfixn ls1 = BUSZSI(rep12, 1)+1;
	sfixn ls2 = BUSZSI(rep12, 2)+1;
	sfixn es1 = logceiling(ls1);
	sfixn es2 = logceiling(ls2);

	checkFrourierPrime(es1, pPtr);
	checkFrourierPrime(es2, pPtr);

	sfixn dims1 = 1<<es1;
	sfixn dims2 = 1<<es2;
	sfixn n = ls1*ls2;

// 	sfixn maxdim = dims1;
// 	if (dims1 < dims2)
// 	  maxdim = dims2;

	sfixn *rootsPtr=(sfixn *) my_calloc(dims1+dims2, sizeof(sfixn));

	//--------------------------------------
	cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es1, dims1, rootsPtr, pPtr);
	PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dims2, rootsPtr+dims1, pPtr);
	cilk_sync;
	//--------------------------------------

	//------------------------------------
	cilk_spawn PBPAS::bivarEvalBy2DTFT_spawn(DAT(rep12), dims1, dims2, ls1, ls2, DAT(rep1), BUSZS(rep1), rootsPtr, pPtr);  

	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	
	PBPAS::bivarEvalBy2DTFT_spawn(tftData2, dims1, dims2, ls1, ls2, DAT(rep2), BUSZS(rep2), rootsPtr, pPtr); 
	cilk_sync;
	//----------------------------------------------

	PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, DAT(rep12), tftData2, pPtr->P);

	my_free(tftData2);

	//now data is in ls2xls1 
	PBPAS::bivarInterpolBy2DTFT_spawn(DAT(rep12), dims1, dims2, ls1, ls2, rootsPtr, pPtr);
	
  }

  //--------------------------------------------
  //bivariate, ls1 and ls2 are not 0
  void bivarInterpolBy2DTFT_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	sfixn * tmprootsPtr = rootsPtr +dims1;
	size_t s = ls2*sizeof(sfixn);
	for(int j=0; j<ls1; j++){
	  cilk_spawn PBPAS::copyDataOutInvTFT(resCoeffs+j*ls2, dims2, s, ls2, tmprootsPtr, pPtr);
	}
	cilk_sync;

	if ( ls1==ls2 )
	  PBPAS::sqtranspose(resCoeffs, ls1, 0, ls2, 0, 1, ls2);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(ls1*ls2, sizeof(sfixn));
	  PBPAS::transpose(resCoeffs, ls2, B, ls1, 0, ls1, 0, ls2);
#if __GNUC__ == 4
	  memcpy(resCoeffs, B, n*sizeof(sfixn));
#else
	  cilk_for(int j = 0; j < n; j++) resCoeffs[j] = B[j];
#endif
	  my_free(B);
	}	

	s = ls1*sizeof(sfixn);

	for(int j=0; j<ls2; j++){
	  cilk_spawn PBPAS::copyDataOutInvTFT(resCoeffs+j*ls1, dims1, s, ls1, rootsPtr, pPtr);
	}
	cilk_sync;
  }

  void copyDataOutInvTFT(sfixn * data, sfixn dims, size_t size, sfixn ls, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	sfixn * tmpVec = (sfixn *)my_calloc(dims, sizeof(sfixn));
	memcpy(tmpVec, data, size);
	EX_Mont_INVTDFT_OPT2_AS_GENE_1( ls, rootsPtr, tmpVec, pPtr);
	memcpy(data, tmpVec, size);
	my_free(tmpVec);
  }

  //--------------------------------------------
  //bivariate, ls1 and ls2 are not 0
  void bivarEvalBy2DTFT_spawn(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn * pdegs, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	
	sfixn s1 = pdegs[1]+1;
	size_t s = s1*sizeof(sfixn);
	size_t res = ls1*sizeof(sfixn);
	for(int j=0; j<=pdegs[2]; j++){
	  //evaluate 1st var for the num of pdegs[2]+1 coefficients
	  //do not use ls2. saving from not fft on padding zeros
	  //for parallel run

	  cilk_spawn PBPAS::copyDataOutTFT(coeffs+j*s1, resCoeffs+j*ls1, dims1, s, res, ls1, rootsPtr, pPtr);
	}
	cilk_sync;
	
	if ( ls1==ls2 )
	  PBPAS::sqtranspose(resCoeffs, ls2, 0, ls1, 0, 1, ls1);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
	  PBPAS::transpose(resCoeffs, ls1, B, ls2, 0, ls2, 0, ls1);
#if __GNUC__ == 4
	  memcpy(resCoeffs, B, n*sizeof(sfixn));
#else
	  cilk_for(int j = 0; j < n; j++) resCoeffs[j] = B[j];
#endif
	  my_free(B);
	}

	sfixn * tmprootsPtr=rootsPtr+dims1;
	s = ls2*sizeof(sfixn);
	for(int j=0; j<ls1; j++){
	  cilk_spawn PBPAS::copyDataOutTFT(resCoeffs+j*ls2, resCoeffs+j*ls2, dims2, s, s, ls2, tmprootsPtr, pPtr);

	}
	cilk_sync;
	//now resCoeffs' layout is ls2 x ls1
  }

  void copyDataOutTFT(sfixn * indata, sfixn * outdata, sfixn dims, size_t insize, size_t outsize, sfixn ls, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	sfixn * tmpVec = (sfixn *)my_calloc(dims, sizeof(sfixn));
	memcpy(tmpVec, indata, insize);
	EX_Mont_TDFT_OPT2_AS_GENE_1( ls, rootsPtr, tmpVec, pPtr);
	memcpy(outdata, tmpVec, outsize);
	my_free(tmpVec);
  }

  //-------------------------------------
  // use cilk_for for bivariate multiplication
  // N=2
  // no copy to init KN construct, and saved TFT in dimension-1
  // avoid TFT work on zero vectors for the evaluation of 1st variable
  // for 1D, copy directly from rep1 to a maximum vector for TFT and
  // put result to rep12, then for dimension-2 TFT 
  // saved 30% w.r.t original 2D-TFT
  void bivarMultiplyBy2DTFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
	sfixn ls1 = BUSZSI(rep12, 1)+1;
	sfixn ls2 = BUSZSI(rep12, 2)+1;
	sfixn es1 = logceiling(ls1);
	sfixn es2 = logceiling(ls2);

	//checkFrourierPrime(es1, pPtr);
	//checkFrourierPrime(es2, pPtr);
	if ( (es1 > (pPtr->Npow)) || (es2 > (pPtr->Npow)) )
	  exit (1); //checkFrourierPrime

	sfixn dims1 = 1<<es1;
	sfixn dims2 = 1<<es2;
	sfixn n = ls1*ls2;

// 	sfixn maxdim = dims1;
// 	if (dims1 < dims2)
// 	  maxdim = dims2;

	sfixn *rootsPtr=(sfixn *) my_calloc(dims1+dims2, sizeof(sfixn));
	//--------------------------------------
	cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es1, dims1, rootsPtr, pPtr);	
	PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dims2, rootsPtr+dims1, pPtr);
	cilk_sync;
	
	//std::cout<<"bivarMultiplyBy2DTFT pPtr->P = "<<pPtr->P<<std::endl;

	//std::cout<<"powers of root in X direction -------"<<std::endl;
	//for (int i=0; i<dims1; ++i)
	//  std::cout<<rootsPtr[i]<<", ";
	//std::cout<<std::endl;
	//std::cout<<"powers of root in Y direction -------"<<std::endl;
	//for (int i=dims1; i<dims1+dims2; ++i)
	//  std::cout<<rootsPtr[i]<<", ";
	//std::cout<<std::endl;

#ifdef TIMEFN
	//------------------------------------
	//run two evaluations in serial to measure the performance of 2D-TFT
	long start = example_get_time();
	PBPAS::bivarEvalBy2DTFT(DAT(rep12), dims1, dims2, ls1, ls2, DAT(rep1), BUSZS(rep1)[1], BUSZS(rep1)[2], rootsPtr, pPtr);  
	long end = example_get_time();
	float oneEvaltime = (end - start) / 1.f;
	std::cout<< oneEvaltime << " ";

	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));
	
	PBPAS::bivarEvalBy2DTFT(tftData2, dims1, dims2, ls1, ls2, DAT(rep2), BUSZS(rep2)[1], BUSZS(rep2)[2], rootsPtr, pPtr); 

	//----------------------------------------------
#else	

	cilk_spawn 
	  PBPAS::bivarEvalBy2DTFT(DAT(rep12), dims1, dims2, ls1, ls2, DAT(rep1), BUSZS(rep1)[1], BUSZS(rep1)[2], rootsPtr, pPtr);

	//for the 2nd tft transform in the evaluation of rep2
	sfixn * tftData2 = (sfixn * )my_calloc(SIZ(rep12), sizeof(sfixn));	    
	PBPAS::bivarEvalBy2DTFT(tftData2, dims1, dims2, ls1, ls2, DAT(rep2), BUSZS(rep2)[1], BUSZS(rep2)[2], rootsPtr, pPtr); 
	
	cilk_sync;
	
	//----------------------------------------------
#endif

	//PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, DAT(rep12), tftData2, pPtr->P); //diff by R?
	PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n, DAT(rep12), tftData2, pPtr);

	my_free(tftData2);

#ifdef TIMEFN
	long start1 = example_get_time();
#endif
	//now data is in ls2xls1 
	PBPAS::bivarInterpolBy2DTFT(DAT(rep12), dims1, dims2, ls1, ls2, rootsPtr, pPtr);

	my_free(rootsPtr);
#ifdef TIMEFN	
	long end1 = example_get_time();
	float interp1 = (end1 - start1) / 1.;
	std::cout<< interp1 << " ";
#endif
  }

  //--------------------------------------------
  //bivariate, ls1 and ls2 are not 0
  void bivarInterpolBy2DTFT(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
	//input resCoeffs' layout is ls2 x ls1, from evaluation
	sfixn * tmprootsPtr = rootsPtr +dims1;
	size_t s = ls2*sizeof(sfixn);

	//#pragma cilk_grainsize = 1;
	cilk_for(int j=0; j<ls1; j++){
	  sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn));
	  sfixn * pos = resCoeffs+j*ls2;
	  memcpy(tmpVec, pos, s);
	  EX_Mont_INVTDFT_OPT2_AS_GENE_1( ls2, tmprootsPtr, tmpVec, pPtr);
	  memcpy(pos, tmpVec, s);
	  
	  my_free(tmpVec);
	}

	if ( ls1==ls2 )
	  PBPAS::sqtranspose(resCoeffs, ls1, 0, ls2, 0, 1, ls2);
	else{	
	  sfixn n = ls1*ls2;
	  sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
	  PBPAS::transpose(resCoeffs, ls2, B, ls1, 0, ls1, 0, ls2);	  
	  memcpy(resCoeffs, B, n*sizeof(sfixn));
	  my_free(B);
	}

	s = ls1*sizeof(sfixn);

	//#pragma cilk_grainsize = 1;
	cilk_for(int j=0; j<ls2; j++){
	  sfixn * tmpVec = (sfixn *)my_calloc(dims1, sizeof(sfixn));
	  sfixn * pos = resCoeffs+j*ls1;
	  memcpy(tmpVec, pos, s);
	  EX_Mont_INVTDFT_OPT2_AS_GENE_1( ls1, rootsPtr, tmpVec, pPtr);
	  memcpy(pos, tmpVec, s);
	  
	  my_free(tmpVec);
	}
  }

  //--------------------------------------------
  //bivariate, ls1 and ls2 are greater than 0
  void  bivarEvalBy2DTFT(sfixn * resCoeffs, sfixn dims1, sfixn dims2, sfixn ls1, sfixn ls2, sfixn * coeffs, sfixn pd1, sfixn pd2, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr){
  
    sfixn s1 = pd1+1;
    size_t s = s1*sizeof(sfixn);
    size_t res = ls1*sizeof(sfixn);
   
    //#pragma cilk_grainsize = 1;
    cilk_for(sfixn j=0; j<=pd2; j++){
      //evaluate 1st var for the num of pd1+1 coefficients
      //do not use ls2. saving from not fft on padding zeros
      //for parallel run

      sfixn* tmpVec = (sfixn *)my_calloc(dims1, sizeof(sfixn));     	  
      memcpy(tmpVec, coeffs+j*s1, s);
      EX_Mont_TDFT_OPT2_AS_GENE_1(ls1, rootsPtr, tmpVec, pPtr);

      memcpy(resCoeffs+j*ls1, tmpVec, res); 
      my_free(tmpVec);
    }
    
    if ( ls1==ls2 )
      PBPAS::sqtranspose(resCoeffs, ls2, 0, ls1, 0, 1, ls1);
    else{	
      sfixn n = ls1*ls2;
      sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
      PBPAS::transpose(resCoeffs, ls1, B, ls2, 0, ls2, 0, ls1);
      memcpy(resCoeffs, B, n*sizeof(sfixn));
      my_free(B);
    }
      
      sfixn * tmprootsPtr=rootsPtr+dims1;
      s = ls2*sizeof(sfixn);
      
      //#pragma cilk_grainsize = 1;
      cilk_for(int j=0; j<ls1; j++){
	sfixn * tmpVec = (sfixn *)my_calloc(dims2, sizeof(sfixn));
	sfixn * pos = resCoeffs+j*ls2;
	memcpy(tmpVec, pos, s);
	EX_Mont_TDFT_OPT2_AS_GENE_1(ls2, tmprootsPtr, tmpVec, pPtr);

	memcpy(pos, tmpVec, s);
	
	my_free(tmpVec);
      }
      //now resCoeffs' layout is ls2 x ls1
  }

  //================================================================
  // Classical Multiplication.
  // dgs1[0], dgs2[0] are useless.
  // ccumx[i] keeps the base before ith-dimension.
  //================================================================
  void plainMultiDMul(sfixn N, sfixn * ccum, sfixn * res, sfixn * ccum1, sfixn * dgs1, sfixn * ccum2, sfixn * dgs2, sfixn * coeffs1, sfixn * coeffs2, MONTP_OPT2_AS_GENE * pPtr){
	int d;
	
	sfixn p=pPtr->P, R=(1L<<pPtr->Rpow)%p, SFT=pPtr->Base_Rpow;

	d=shrinkDeg(dgs1[N], coeffs1, ccum1[N]);

	//cilk_for(int i=0; i<=d; i++){
	for(int i=0; i<=d; i++){//race with inline MontMulMod_OPT2_AS_GENE:a += (a >> BASE_1) & c;
	  PBPAS::decomposePoly(N, ccum, res+ccum[N]*i,  N-1, dgs1, coeffs1+ccum1[N]*i, ccum1, N, dgs2, coeffs2, ccum2, pPtr, R, SFT);
	}
  } 


  //-------------------------------
  void
  decomposePoly(sfixn N, sfixn * ccum, sfixn * res, sfixn N1, sfixn * dgs1, sfixn * coeffs1, sfixn * ccum1, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr, sfixn R, sfixn SFT){
	int d;
	if(N1==0){
	  if(! coeffs1[0]) return;
	  R=( MulMod(coeffs1[0], R, pPtr->P))<<SFT;    
	  PBPAS::decomposePoly2(N, ccum, res, R, N2, dgs2, coeffs2, ccum2, pPtr); 
	  return;
	}

	d=shrinkDeg(dgs1[N1], coeffs1, ccum1[N1]);
	//cilk_for(int i=0; i<=d; i++){//race with inline MontMulMod_OPT2_AS_GENE:a += (a >> BASE_1) & c;
	for(int i=0; i<=d; i++){
	  PBPAS::decomposePoly(N, ccum, res+ccum[N1]*i, N1-1, dgs1, coeffs1+ccum1[N1]*i, ccum1, N2, dgs2, coeffs2, ccum2, pPtr, R, SFT);
	}	
  }

  //-------------------------------------
  void
  decomposePoly2(sfixn N, sfixn * ccum, sfixn * res, sfixn num, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr){
	int d;
	if(N2==0){
      res[0]=AddMod(res[0],  MontMulMod_OPT2_AS_GENE(coeffs2[0], num, pPtr), pPtr->P);
	  return;
	}
	
	if(N2==1){ 
	  d=shrinkDeg(dgs2[1], coeffs2, ccum2[1]);
	  cilk_for(int i=0; i<=d; i++) 
	  //for(int i=0; i<=d; i++)
        PBPAS::decomposePoly2(N, ccum, res+ccum[1]*i, num, 0, dgs2, coeffs2+ccum2[1]*i, ccum2, pPtr);
	  return;
	}
	
	d=shrinkDeg(dgs2[N2], coeffs2, ccum2[N2]);
	cilk_for(int i=0; i<=d; i++){
	  //for(int i=0; i<=d; i++){
	  PBPAS::decomposePoly2(N, ccum, res+ccum[N2]*i, num, N2-1, dgs2, coeffs2+ccum2[N2]*i, ccum2, pPtr);
	}	
  }


  // inputPolydegs is dgs.
  // useful data is dgs[1..N].
  /**
   * getRevInvTiSet:
   * @dgs: The partial degrees of the triangular set.
   * @N: The number of variables.
   * @tris: (output) the inverses (reverse ordered) of ts.
   * @ts: The triangular set.
   * @pPtr: The information of the prime.
   *
   * This is used for fast division (or normal form) modulo a triangular set.
   * We pre-compute with this functionn the inverse of the polynomial `ts.i`
   * w.r.t. 'Xi' modulo the underlying triangular set.
   *
   * Compute the inverses (reverse ordered) of given triangular set 'ts'.
   *
   * Return value: 
   **/
  void getRevInvTiSet(sfixn *dgs, sfixn N, TriRevInvSet * tris, TriSet * ts,  MONTP_OPT2_AS_GENE * pPtr){
    register int i;
    preFFTRep tRIPtrtmp;

    for(i=1; i<=N; i++) {
	  SERBPAS::InitOneRevInv(i, &tRIPtrtmp, BDS(ts), dgs[i]);

      if(!EXSTI(tris, i))
		//=========================================================
		PBPAS::NewtonRevInverse(i, &tRIPtrtmp, ELEMI(tris, i), ELEMI(ts, i), ts, tris, pPtr );
	  //=========================================================

      freePoly(&tRIPtrtmp); //modpn MPMMTS.h
	}
  }

  //====================================================================
  //====================================================================
  //
  //
  //
  //     Newton Reverse Inverse.
  //
  //
  //====================================================================
  //====================================================================
  // output tRIPt.
  /**
   * NewtonRevInverse:
   * @N: number of variables of 'ts'.
   * @tRIPtrtmp: temporary buffer.
   * @tRIPtr:(output) the inverse of T_N modulo {T_1,...,T_{N-1}}.
   * @tPtr: a triangular set {T_1,...,T_N}.
   * @tris: the inverses of 'tPtr'.
   * @pPtr: the information of prime p.
   *
   * Compute the inverse of T_N modulo 'tPtr'.
   *
   * Return value: 
   **/
  void NewtonRevInverse(sfixn N, preFFTRep * tRIPtrtmp, preFFTRep * tRIPtr, preFFTRep * tPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	register int i; 
	sfixn e,n;
	sfixn * tmpPtr;
	sfixn cutDeg;
	sfixn degF;
	sfixn d1, d2, d3, d4, d5;
	KroFFTRep * krofftrep, * krofftrep2;
	preFFTRep * resPoly, * resPoly2;
	
	backupData(tPtr); //inlineFuncs.h
	
	if(N==1){ 
	  n=(BUSZSI(tRIPtr, 1))+1;
	  e=logceiling(n);
	  tmpPtr=(sfixn *)my_calloc(((BUSZSI(tPtr, 1))+1),(sizeof(sfixn)));
	  PBPAS::reverseUni(BUSZSI(tPtr, 1), tmpPtr, DAT(tPtr));
	  degF=BUSZSI(tPtr, 1);

	  //=========================================================
	  PBPAS::modularInvPM(BUSZSI(tRIPtr, 1), DAT(tRIPtr), degF, tmpPtr, e, n, pPtr);
	  //=========================================================
	  //use dft, how about when DAT(tRIPtr) is small?

	  my_free(tmpPtr);
	  EXSTI(tris, N)=1;
	  return;
	}

	krofftrep=(KroFFTRep *) my_calloc(1,sizeof(KroFFTRep));
	krofftrep2=(KroFFTRep *) my_calloc(1,sizeof(KroFFTRep));
	resPoly=(preFFTRep *) my_calloc(1,sizeof(preFFTRep));
	resPoly2=(preFFTRep *) my_calloc(1,sizeof(preFFTRep));
	n=(BUSZSI(tRIPtr, N))+1;
	e=logceiling(n);
	tmpPtr=(sfixn *)my_calloc(SIZ(tPtr),(sizeof(sfixn)));

	PBPAS::reverseMulti(BUSZSI(tPtr, N), CUMI(tPtr, N), tmpPtr, DAT(tPtr));
	setData(tPtr, tmpPtr); //inline

	DATI(tRIPtr, 0)=1;
	SERBPAS::InitResPoly(resPoly, N, BUSZS(tRIPtr) , BUSZS(tRIPtr));
	SERBPAS::InitResPoly(resPoly2, N, BUSZS(tPtr) , BUSZS(tRIPtr));
	SERBPAS::InitKroFFTRep(krofftrep, BUSZS(resPoly), N, 2, pPtr);
	SERBPAS::InitKroFFTRep(krofftrep2, BUSZS(resPoly2), N, 2, pPtr);
	
	//=====================
	// Newton Iteration  ==
	//=====================
	
	d1=BUSZSI(tRIPtr, N);
	d2=BUSZSI(tPtr, N);
	d3=BUSZSI(resPoly, N);
	d4=BUSZSI(tRIPtrtmp, N);
	d5=BUSZSI(resPoly2, N);
	
	for (i=1; i<=e; i++){ //serial x^e
	  cutDeg=(1<<i)-1;
	  
	  if(d1>cutDeg) 
		BUSZSI(tRIPtr, N)=cutDeg; 
	  else 
		BUSZSI(tRIPtr, N)=d1;

	  if(d2>cutDeg) 
		BUSZSI(tPtr, N)=cutDeg; 
	  else 
		BUSZSI(tPtr, N)=d2;
	  
	  if(i>1) 
		PBPAS::PolyCleanData(resPoly);
	  
	  if(d3>(2*cutDeg)) 
		BUSZSI(resPoly, N)=2*cutDeg; 
	  else 
		BUSZSI(resPoly, N)=d3;
	  
	  SERBPAS::decreaseKroFFTRep(krofftrep, BUSZS(resPoly));
	  
	  //=========================================================
	  PBPAS::squarePoly_FFT(N, krofftrep, resPoly, tRIPtr, pPtr);
	  // =============================================================
	  
	  
	  PBPAS::PolyCleanData(tRIPtrtmp);
	  
	  if(d4> cutDeg) 
		BUSZSI(tRIPtrtmp, N)=cutDeg; 
	  else 
		BUSZSI(tRIPtrtmp, N)=d4;

	  if(d3>cutDeg) 
		BUSZSI(resPoly, N)=cutDeg; 
	  else 
		BUSZSI(resPoly, N)=d3;

	  //=========================================================
	  PBPAS::reduceCoeffs(N, tRIPtrtmp, resPoly, ts, tris, pPtr);
	  //=========================================================

	  if(i>1)  
		PBPAS::KroneckerCleanData(krofftrep2);

	  if(i>1) 
		PBPAS::PolyCleanData(resPoly2);
	  
	  if(d5>2*cutDeg) 
		BUSZSI(resPoly2, N)=2*cutDeg; 
	  else 
		BUSZSI(resPoly2, N)=d5;
	  
	  SERBPAS::decreaseKroFFTRep(krofftrep2, BUSZS(resPoly2));
	  
	  //====================================================================
	  PBPAS::mulPoly_FFT(N, krofftrep2, resPoly2, tRIPtrtmp, tPtr, pPtr);
	  //==================================================================== 
	  //actually kn fft, can not handle 4 65 65 65 65 1
	  //shall use multi dft/fft
	  
	  if(d5>cutDeg) 
		BUSZSI(resPoly2, N)=cutDeg; 
	  else 
		BUSZSI(resPoly2, N)=d5;
	  
	  BUSZSI(tRIPtrtmp, N)=d4;
	  PBPAS::PolyCleanData(tRIPtrtmp);
	  BUSZSI(tRIPtrtmp, N)=cutDeg;
	  
	  PBPAS::reduceCoeffs(N, tRIPtrtmp, resPoly2, ts, tris, pPtr);
	  
	  //============================================
	  PBPAS::addEqDgPoly_1(N, tRIPtr, tRIPtr, pPtr->P); //inline
	  //============================================

	  //=================================================
	  PBPAS::subEqDgPoly_1(N, tRIPtr, tRIPtrtmp, pPtr->P, 1);//inline
	  //=============================================
	  
	}
	
	// actually free tmpPtr -- the reversal of tPtr->data.
	my_free(DAT(tPtr));
	restoreData(tPtr);
		
	BUSZSI(tRIPtr, N)=d1;
	BUSZSI(tPtr, N)=d2;
	BUSZSI(resPoly, N)=d3;
	BUSZSI(tRIPtrtmp, N)=d4;
	BUSZSI(resPoly2, N)=d5;
	
	freePoly(resPoly);  
	freePoly(resPoly2); 
	SERBPAS::freeKroFFTRep(krofftrep);
	SERBPAS::freeKroFFTRep(krofftrep2);
	my_free(krofftrep);
	my_free(krofftrep2);
	my_free(resPoly);
	my_free(resPoly2);
	EXSTI(tris, N)=1; 
	//
	return;
  }
  
  // Input: degG=n-1, degF<l where x^l is the modulus.
  //        degF< l implies F is reduced r.w.t x^l.
  //        r=ceiling(log_2(l)), n=2^r.
  // Output: G with degG=n-1. G is the inverse of F modulo x^r.
  //         =>  G is also the inverse of F modulo x^l.
  // note: good for all Fourier Primes.
  /**
   * modularInvPM: Computes the inverse of F modulo x^r 
   *               This assumes F has a non-zero trailing coefficient.
   *               This works for all Fourier Primes.
   * @degG: degree of G (output) 
   * @GPtr: coefficient vector of G (output) 
   * @degF: degree of F
   * @FPtr: coefficient vector of F
   * @r: 
   * @n: equals 2^r
   * @pPtr: prime number structure
   * 
   * Using the Middle product trick.
   * So the running time is expected to be in 2 M(n) + o(M(n)) machine 
   * word operations.
   * Return value: G, the inverse of F modulo x^r 
   **/   
  sfixn *
  modularInvPM(sfixn degG, // degG=n-1;
			   sfixn * GPtr, sfixn degF, sfixn * FPtr, 
			   sfixn r, sfixn n,
			   MONTP_OPT2_AS_GENE * pPtr)
  {
	int i;
	sfixn nn, halfnn;
	sfixn * rootsPtr=(sfixn *)my_calloc(n, sizeof(sfixn)), 
	  * tmpGPtr=(sfixn *)my_calloc(n, sizeof(sfixn)), 
	  * tmpFPtr=(sfixn *)my_calloc(n, sizeof(sfixn));
	
	GPtr[0]=1;
	
	for(i=1; i<=r; i++){ //serial, dependent on previous GPtr
	  nn=1<<i;
	  halfnn=nn>>1;
	  
	  SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(i, nn, rootsPtr, pPtr);
	  
	  EX_Mont_DFT_OPT2_AS_GENE( nn, i, rootsPtr, tmpGPtr, nn-1, GPtr, pPtr);
	  
	  if(degF>=(nn-1)) 
		EX_Mont_DFT_OPT2_AS_GENE( nn, i, rootsPtr, tmpFPtr, nn-1, FPtr, pPtr);
	  else    
		EX_Mont_DFT_OPT2_AS_GENE( nn, i, rootsPtr, tmpFPtr, degF, FPtr, pPtr);

	  PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(nn, tmpFPtr, tmpGPtr, pPtr);
	  
	  
	  EX_Mont_INVDFT_OPT2_AS_GENE_R_1 (nn, i, rootsPtr, tmpFPtr, pPtr);
	  
	  cilk_for(int j=0; j<halfnn; j++) 
		tmpFPtr[j]=tmpFPtr[j+halfnn];
	  
	  cilk_for(int j=halfnn; j<nn; j++) 
		tmpFPtr[j]=0;
	  
	  
	  EX_Mont_DFT_OPT2_AS_GENE_1 ( nn, i, rootsPtr, tmpFPtr, pPtr );
	  
	  PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(nn, tmpFPtr, tmpGPtr, pPtr);
	  
	  
	  EX_Mont_INVDFT_OPT2_AS_GENE_R_1 (nn, i, rootsPtr, tmpFPtr, pPtr);
	  
	  cilk_for(int j=halfnn; j<nn; j++) 
		GPtr[j] = SubMod(GPtr[j], tmpFPtr[j-halfnn], pPtr->P);
	  
	};
	
	my_free (rootsPtr);
	my_free (tmpGPtr);
	my_free (tmpFPtr);
	return GPtr;
  }
  
  //===========================================================
  // rPtr = p1Ptr ^ 2
  //===========================================================
  /**
   * squarePoly_FFT:
   * @N: number of variables.
   * @kPtr: the Kronecker presentation.
   * @rPtr: (output) C-Cube polynomial to keep the square (p1)^2.
   * @p1Ptr: C-Cube prolynomial p1.
   * @pPtr: Information for the prime number.
   * use KN FFT or FFT; shall change to 2D TFT or plain multiplication
   * 
   * Return value: 
   **/
  void squarePoly_FFT(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep *p1Ptr,  MONTP_OPT2_AS_GENE * pPtr){
    sfixn kdg1;
    int switcher=0;
	
    if((KE(kPtr))>(pPtr->Npow)) 
      { 
		//printf("FFT size larger than the prime %ld can handle-> Using MultiD-FFT solving this problem.", pPtr->P);
         switcher=1;}
	
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(p1Ptr), BUSZS(p1Ptr), DAT(p1Ptr));
	
    if(switcher){
      //
	  //
	  PBPAS::fftMultiD_square_test_1(DATSI(kPtr, 0), N, ES(kPtr), DIMS(kPtr), pPtr); }
    else{ //use KN 1D FFT
      //
	  kdg1= (BUSZSI(p1Ptr, N)+1)*CUMI(kPtr, N)-1;
	  KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));
	  //
	  SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr),pPtr);
	  //
	  PBPAS::EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), pPtr);

    }
	
    PBPAS::fromtofftRepMultiD(N,  CUM(rPtr), DAT(rPtr), CUM(kPtr), BUSZS(rPtr), DATSI(kPtr, 0));
  }
  
  //================================================================
  // Squaring by Multi-dimensional FFT.
  //
  // es[0] and dims[0] are useless;
  // es[i] is the exponent of dims[i];
  // dims[i] is the # on dims-i
  //================================================================
  //parallelized
  void fftMultiD_square_test_1(sfixn * coeffs1, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr)
  {
	register int i;
    sfixn m=0, n=1;
    sfixn *rootsPtr, *tmprootsPtr;
	sfixn *tmprootsPtr1;
	
	//std::cout << "entering fftMultiD_test_1" << std::endl;

	for(i=1; i<=N; i++) { //N is small
	  checkFrourierPrime(es[i], pPtr);
	}

    //compute the primitive roots for all dimensions
    for(i=1; i<=N; i++) { n *= dims[i]; m += dims[i]; } 
    //use reducer?? N is small

    rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));

    tmprootsPtr1 = rootsPtr;	

    for(i=1; i<=N; i++) {//no parallel, dependency
      SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr1, pPtr);
      tmprootsPtr1 += dims[i];
    }

    tmprootsPtr = PBPAS::MultiEvalByFFT(coeffs1, N, n, es, dims, pPtr, rootsPtr);

    // Pairwise-Mul, result is in coeffs1
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs1, pPtr->P);

    PBPAS::InterpolByFFT(coeffs1, N, n, es, dims, pPtr, tmprootsPtr);

	//std::cout <<"before exit fftMultiD_square_test_1" << std::endl;
    my_free(rootsPtr);

  }
  
  // fastest FFT and good for FFT-primes.
  /**
   * EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1:
   * @n: the FFT size
   * @e: log2 of n
   * @degRes: degree of result
   * @rootsPtr: powers of the primitive n-th root of unity
   * @degA: (Input and Output) degree of A
   * @APtr: coefficient vector of A 
   * @pPtr: prime number structure
   * 
   * FFT-based square of A. In place: result in A.
   * 
   * 
   * Return value: the square of a polynomial.
   **/   
  void EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn * rootsPtr, sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr){
	sfixn l_AB, e1, n1, sA;
	
    l_AB = degRes + 1; // compute l_AB for FFT(A,B).
    e1 = logceiling(l_AB);

	checkFrourierPrime(e1, pPtr);
	
    n1 = 1L << e1;
    
    // force this FFT use the n from parameter if n<n1.
    if ( (n<n1) || (degRes==0) ) { n1=n; e1=e;}
	
	sA = degA+1;
    cilk_for(int i=sA; i<n1; i++) APtr[i]=0;
    
    EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, APtr, pPtr);
    
    PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n1, APtr, APtr, pPtr);
    
    EX_Mont_INVDFT_OPT2_AS_GENE_R_1(n1, e1, rootsPtr, APtr, pPtr);
    	
  }
  
  /**
   * reduceCoeffs:
   * @N: number of variables.
   * @toPtr: (output) 'fromPtr' modulo 'ts'.
   * @fromPtr: a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverses of 'ts'.
   * @pPtr: information of the prime.
   * 
   * 'ts' has main variable N-1. 'fromPtr' has main variable N.
   * Reduce the coefficients of 'fromPtr' by 'ts'. 
   *
   * Return value: 
   **/
  void reduceCoeffs(sfixn N, preFFTRep * toPtr, preFFTRep * fromPtr, TriSet * ts, TriRevInvSet * tris,   MONTP_OPT2_AS_GENE * pPtr)
  { 
// register int i;
	
// 	backupData(fromPtr); //inline
// 	backupData(toPtr);

// 	decreaseOneDim(fromPtr); //inline
// 	decreaseOneDim(toPtr);

// 	for(i=0; i<=BUSZSI(toPtr, N); i++){
// 	  //===============================================
// 	  SERBPAS::MultiMod(N-1, toPtr, fromPtr, ts, tris, pPtr);
// 	  //===============================================
// 	  nextCoefData(fromPtr, N); //inline DAT(fromPtr)+=CUMI(fromPtr, N);
// 	  nextCoefData(toPtr, N); 
// 	}

// 	increaseOneDim(fromPtr); //inline
// 	increaseOneDim(toPtr);

// 	restoreData(fromPtr); //inline
// 	restoreData(toPtr);

	sfixn deg = BUSZSI(toPtr, N);
	preFFTRep ** arrayTo = (preFFTRep **)my_malloc((deg+1)*sizeof(preFFTRep *));
	preFFTRep ** arrayFrom = (preFFTRep **)my_malloc((deg+1)*sizeof(preFFTRep *));
	
	cilk_for(int i=0; i<=deg; i++){
	  arrayTo[i] = SERBPAS::CopyOnePoly_sameData(toPtr);
	  arrayFrom[i] = SERBPAS::CopyOnePoly_sameData(fromPtr);
	  //arrayTo[i] = (preFFTRep *)my_malloc(sizeof(preFFTRep));
	  //SERBPAS::copyPolyPointers(arrayTo[i], toPtr);
	  //race between copyVec_1_to_n and CopyOnePoly_sameData 
	  N(arrayTo[i]) = N-1;
	  SIZ(arrayTo[i]) = CUMI(toPtr, N);
	  DAT(arrayTo[i]) = DAT(toPtr)+i*CUMI(toPtr, N);

	  //arrayFrom[i] = (preFFTRep *)my_malloc(sizeof(preFFTRep));
	  //SERBPAS::copyPolyPointers(arrayFrom[i], fromPtr);
	  N(arrayFrom[i]) = N-1;
	  SIZ(arrayFrom[i]) = CUMI(fromPtr, N);
	  DAT(arrayFrom[i]) = DAT(fromPtr)+i*CUMI(fromPtr, N);

	  PBPAS::MultiMod(N-1, arrayTo[i], arrayFrom[i], ts, tris, pPtr);

	  SERBPAS::freePoly_not_data(arrayTo[i]);
	  SERBPAS::freePoly_not_data(arrayFrom[i]);
	}

	//cilk_for(int i=0; i<=deg; i++){
	//SERBPAS::freePoly_not_data(arrayTo[i]);
	//SERBPAS::freePoly_not_data(arrayFrom[i]);
	//}

//use MultiMod_DF
// 	if (zeroCoefp(DAT(fromPtr), SIZ(fromPtr)) ) return;
// 	PBPAS::MultiMod_1_par(N, fromPtr, ts, tris, pPtr);
// 	PBPAS::fromtofftRep(N, CUM(toPtr), DAT(toPtr), CUM(fromPtr),
// 	BUSZS(toPtr), DAT(fromPtr));

  }
  
  /**
   * MultiMod:
   * @N: number of variables.
   * @outPtr: (output) the image of 'inPtr' modulo 'ts'.
   * @inPtr: a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverse of 'ts'.
   * @pPtr: the information of prime number p.
   * 
   * Reduce 'inPtr' w.r.t 'ts' and save the result in 'outPtr'.
   *
   * Return value: 
   **/
  void
  MultiMod(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr)
  { 
	if(zeroPolyp(inPtr)) return;
	
	//----yu
	//why 1 here to use MultiMod_DF
	//but in modpn:-NormalForm passed in opt=0 to MultiMod_OPT to use MultiMod_BULL
	if(1){ 
	  // classical/FFT/TFT, depth-first reduction
	  PBPAS::MultiMod_DF(N, outPtr, inPtr, ts, tris, pPtr);	  
	}else
	  {
		// pure TFT, Bottom-up level-by-level reduction.
		SERBPAS::MultiMod_BULL(N, outPtr, inPtr, ts, tris, pPtr);
		
	  }
	
  }
  
  
  //========================================================================
  //   Multiavriate polynomial Fast Mod. (recursive.)
  //   using FFT.
  // 
  // will destruct inPtr.
  // N is the # of current dimensions.
  // outPtr keeps the output of the reduction.
  // -- The output is bounded by ts->bounds (inclusively).
  // inPtr keeps the input to be reduced.
  // -- The input is required bounded by 2 times ts->bounds (inslucively).
  // ts is the triangular set.
  // tris is the inverses of reverse-coef-ordered ts. 
  //========================================================================
  
  /**
   * MultiMod_DF:
   * @N: number of variables.
   * @outPtr: (output) the image of 'inPtr' modulo 'ts'.
   * @inPtr: a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverse of 'ts'.
   * @pPtr: the information of prime number p.
   * 
   * Reduce 'inPtr' w.r.t 'ts' and save the result in 'outPtr'.
   *
   * Return value: 
   **/
  void
  MultiMod_DF(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr)
  { 
	//std::cout<<"PBPAS::MultiMod_DF"<<std::endl;
	
	if (zeroCoefp(DAT(inPtr), SIZ(inPtr)) ) return;
	//preFFTRep * inPtrOri = SERBPAS::CopyOnePoly_sameData(inPtr);
	PBPAS::MultiMod_1_par(N, inPtr, ts, tris, pPtr);
	PBPAS::fromtofftRep(N, CUM(outPtr), DAT(outPtr), CUM(inPtr),  BUSZS(outPtr), DAT(inPtr));
	
  }
  
  //============================================================
  //
  //
  //
  //  Multiavriate polynomial Fast Mod. (recursive.)
  //
  //
  //=============================================================
  /**
   * MultiMod_1_par:
   * @N: number of variables.
   * @inPtr: a C-Cube polynomial.
   * @ts: a zero-dim trinagular set.
   * @tris: inverses of 'ts'.
   * @pPtr: the information of prime number p.
   * //@data0: DAT(inPtr)
   * Reduce 'inPtr' by 'ts' in-place.
   * I.e. 'inPtr'= 'inPtr' modulo 'ts'.
   * Return value: 
   **/
  void MultiMod_1_par(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	//register int i;
	int deg1=0, deg2=0, deg3=0, d;
	
	sfixn * restoreDat;
	sfixn restoreSiz;
	preFFTRep tmpIn;
	
	//std::cout<<"PBPAS::MultiMod_1_par, N="<< N <<", DAT(inPtr)=" << DAT(inPtr)<<std::endl;

	if(N==1){ 
	  
	  deg1 = shrinkDeg(BUSZSI(inPtr, N), DAT(inPtr), 1); //inline

	  //std::cout<<"PBPAS::MultiMod_1_par, N==1, deg1="<<deg1<<std::endl;
	  if (!deg1) return;
	  if ( ( (BDSI(ts, 1)) > DivCut1 ) ) {
		PBPAS::UniFastMod_1(deg1, DAT(inPtr), BUSZSI((ELEMI(ts, 1)), 1), DAT(ELEMI(ts, 1)), (NLB(tris))[1], DAT(ELEMI(tris, 1)), pPtr);  
	  }
	  else{
		//std::cout<<"PBPAS::MultiMod_1_par, N==1, UniPlainMod_1"<<std::endl;
		PBPAS::UniPlainMod_1(deg1, DAT(inPtr), BUSZSI((ELEMI(ts, 1)), 1), 
							 DAT(ELEMI(ts, 1)), pPtr); 
	  }
	  return;
	}
	//std::cout<<"PBPAS::MultiMod_1_par, BUSZSI(inPtr, N)=" <<BUSZSI(inPtr, N)<<std::endl;
	//std::cout<<"PBPAS::MultiMod_1_par, CUMI(inPtr, N)=" <<CUMI(inPtr, N)<<std::endl;
	deg2 = shrinkDeg(BUSZSI(inPtr, N), DAT(inPtr), CUMI(inPtr, N));

	//std::cout<<"PBPAS::MultiMod_1_par, deg2=" <<deg2<<std::endl;

	restoreDat = DAT(inPtr); //????
	restoreSiz = SIZ(inPtr); 

// 	N(inPtr) = N-1;
// 	SIZ(inPtr) = CUMI(inPtr, N);
		
// 	for(int i=0; i<=deg2; i++){
// 	  PBPAS::MultiMod_1_par(N-1, inPtrOri, inPtr, ts, tris, pPtr);
// 	  DAT(inPtr) += CUMI(inPtr, N); 
// 	}

	preFFTRep ** cpinPtr = (preFFTRep **)my_malloc((deg2+1)*sizeof(preFFTRep *));

	cilk_for(int i=0; i<=deg2; i++){
	  cpinPtr[i] = SERBPAS::CopyOnePoly_sameData(inPtr); 
	  //calloc inside

	  //cpinPtr[i]=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	  //SERBPAS::copyPolyPointers(cpinPtr[i], inPtr);

	  N(cpinPtr[i]) = N-1;
	  SIZ(cpinPtr[i]) = CUMI(inPtr, N);
	  DAT(cpinPtr[i]) = DAT(inPtr)+i*CUMI(inPtr, N);
	  //std::cout<<"PBPAS::MultiMod_1_par, arrayPoly[i] i=" << i<<std::endl;
	  //std::cout<<"PBPAS::MultiMod_1_par, BUSZSI=" << BUSZSI(arrayPoly[i], N) <<std::endl;
	  //std::cout<<"PBPAS::MultiMod_1_par, CUMI=" << CUMI(arrayPoly[i], N) <<std::endl;
	  PBPAS::MultiMod_1_par(N-1, cpinPtr[i], ts, tris, pPtr);
	  my_free(cpinPtr[i]);
	}

	//cilk_for(int i=0; i<=deg2; i++){
	//my_free(cpinPtr[i]);
	//}
	//std::cout<<"PBPAS::MultiMod_1_par, after recursive"<<std::endl;
	  
	N(inPtr) = N;
	DAT(inPtr) = restoreDat; 
	SIZ(inPtr) = restoreSiz;

	d = BDSI(ts, N);

	sfixn * cpBDS = (sfixn * )my_calloc(N+1, sizeof(sfixn));
	PBPAS::copyVec_1_to_n(N, cpBDS, BDS(ts));

	PBPAS::copyVec_1_to_n(N-1, CUTS(inPtr), BDS(ts));
	deg3 = shrinkDeg(CUTSI(inPtr, N), DAT(inPtr), CUMI(inPtr, N));
	PBPAS::copyVec_1_to_n(N-1, CUTS(inPtr), BUSZS(inPtr));
	
	//BDSI(ts, N) = deg;
	cpBDS[N] = deg3;
	//std::cout<<"PBPAS::MultiMod_1_par, before InitOneReducedPoly"<<std::endl;
	SERBPAS::InitOneReducedPoly(&tmpIn, N, cpBDS);  
	my_free(cpBDS);

	PBPAS::fromtofftRep(N, tmpIn.accum, tmpIn.data, CUM(inPtr), BUSZSD(tmpIn), DAT(inPtr));
	//std::cout<<"PBPAS::MultiMod_1_par, after fromtofftRep"<<std::endl;
	
	//BDSI(ts, N) = d;
	
	if( (((N==2) && (d>DivCut2)) || ((N==3) && (d>DivCut3)) || ((N>3) && (d>DivCutN))) ){
	  //printf("using MultiUniFastMod_1\n");
	  //fflush(stdout);
	  //std::cout<<"PBPAS::MultiMod_1_par, call PBPAS::MultiUniFastMod_1"<<std::endl;
	  PBPAS::MultiUniFastMod_1(N, &tmpIn, deg3, inPtr, ts, tris, pPtr);	
	  //std::cout<<"PBPAS::MultiMod_1_par, after PBPAS::MultiUniFastMod_1"<<std::endl;
    }
	else{
	  //printf("using MultiUniPlainMod_1\n");
	  //fflush(stdout);
	  //std::cout<<"PBPAS::MultiMod_1_par, call PBPAS::MultiUniPlainMod_1"<<std::endl;
	  PBPAS::MultiUniPlainMod_1(N, deg3, &tmpIn, BDSI(ts, N)+1, ELEMI(ts, N), ts, tris, pPtr);
	  //std::cout<<"PBPAS::MultiMod_1_par, after PBPAS::MultiUniPlainMod_1"<<std::endl;
	  PBPAS::fromtofftRep(N, CUM(inPtr), DAT(inPtr), tmpIn.accum, BUSZSD(tmpIn), tmpIn.data);
	  
	}
	
	freePoly(&tmpIn);  

  }
  
  //============================================================
  //
  //
  //
  //  Multiavriate (as univariate) polynomial Fast Mod.
  //
  //
  //=============================================================
  // Modulo poly by T_n, then by [T_{n-1}, T_{n-2}, ...,T_1] 
  /**
   * MultiUniFastMod_1:
   * @N: number of variables.
   * @tmpIn: a temp buffer.
   * @n: The next power of 2 of (deg(inPtr, X_N)-deg(T_N, X_N)+1).
   * @inPtr:(output) a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverses of 'ts'.
   * @pPtr: the information of the prime p.
   * 
   * Return value: 
   **/
  void MultiUniFastMod_1(sfixn N, preFFTRep * tmpIn,  sfixn n, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	sfixn m,e,d1,d2;
	sfixn cutDg;
	
	preFFTRep * resPoly, * resPoly2, * out;
	KroFFTRep * krofftrep, * krofftrep2;
	
	//
    m=BUSZSI((ELEMI(ts, N)), N);
	
	if(n<m){ // if(DEBUG) printf("Leaving MultiUniFastMod_1.\n"); 
	  return;}
	
	resPoly=(preFFTRep *)my_calloc(1, sizeof(preFFTRep));
	resPoly2=(preFFTRep *)my_calloc(1, sizeof(preFFTRep));
	out=(preFFTRep *)my_calloc(1, sizeof(preFFTRep));
	
	krofftrep=(KroFFTRep *)my_calloc(1, sizeof(KroFFTRep));
	krofftrep2=(KroFFTRep *)my_calloc(1, sizeof(KroFFTRep));
	e=logceiling(n-m+1);
	
	if(! (EXSTI(tris, N))) {  
	  SERBPAS::initCopyOneRevInv(out, ELEMI(tris, N));
	  SERBPAS::NewtonRevInverse(N, out, ELEMI(tris, N), ELEMI(ts, N), ts, tris, pPtr );
	  freePoly(out);
	}

	cutDg = n-m;
	d1 = BUSZSI(tmpIn, N);
	BUSZSI(tmpIn, N) = cutDg;

	//BUSZS--partial degree vector 1..N
	//d2 = BUSZSI(ELEMI(tris, N), N);
	//BUSZSI(ELEMI(tris, N), N) = cutDg;

	PBPAS::reverseMulti_1(n, CUMI(tmpIn, N), DAT(tmpIn));

	//avoid from writing to tris
	preFFTRep * cpNtris = SERBPAS::CopyOnePoly_sameData(ELEMI(tris, N));
	BUSZSI(cpNtris, N) = cutDg;

	//SERBPAS::InitResPoly(resPoly, N, BUSZS(tmpIn), BUSZS(ELEMI(tris, N)));
	SERBPAS::InitResPoly(resPoly, N, BUSZS(tmpIn), BUSZS(cpNtris));

	//ELEMI(tris, N) is a preFFTRep *. Its BUSZS is used in mulPoly
	//=====================================================================
	//SERBPAS::EX_mulPoly(N, resPoly, tmpIn, ELEMI(tris, N), pPtr);
	PBPAS::EX_mulPoly(N, resPoly, tmpIn, cpNtris, pPtr);
	SERBPAS::freePoly_not_data(cpNtris);
	//=====================================================================
	//std::cout<<"after 1st EX_mulPoly"<<std::endl;

	//std::cout<<"MultiUniFastMod_1, cutDg="<<cutDg<<", CUMI(resPoly, N)="<<CUMI(resPoly, N)<<", DAT(resPoly)="<<DAT(resPoly)<<std::endl;
	PBPAS::reverseMulti_1(cutDg, CUMI(resPoly, N), DAT(resPoly));
	//std::cout<<"after reverseMulti_1"<<std::endl;
	BUSZSI(tmpIn, N) = d1;

	//BUSZSI(ELEMI(tris, N), N) = d2;
	
	PBPAS::PolyCleanData(tmpIn);
	
	BUSZSI(tmpIn, N) = cutDg;
	//std::cout<<"before 1st reduceCoeffs"<<std::endl;
	PBPAS::reduceCoeffs(N, tmpIn, resPoly, ts, tris, pPtr);
	//std::cout<<"after 1st reduceCoeffs"<<std::endl;
	freePoly(resPoly); 
	
	SERBPAS::InitResPoly(resPoly2, N, BUSZS(ELEMI(ts, N)), BUSZS(tmpIn));
	
	//=====================================================================
	PBPAS::EX_mulPoly(N, resPoly2, ELEMI(ts, N), tmpIn,  pPtr);
	//=====================================================================
	//std::cout<<"after 2nd EX_mulPoly"<<std::endl;
	BUSZSI(tmpIn, N)=d1;
	
	PBPAS::PolyCleanData(tmpIn);
	
	PBPAS::reduceCoeffs(N, tmpIn, resPoly2, ts, tris, pPtr);
	//std::cout<<"after 2nd reduceCoeffs"<<std::endl;
	freePoly(resPoly2);
	
	PBPAS::subPoly_1(N, inPtr, tmpIn, pPtr->P);
	//std::cout<<"after subPoly_1"<<std::endl;
	
	my_free(resPoly); 
	my_free(resPoly2);
	my_free(out);
	my_free(krofftrep); 
	my_free(krofftrep2);
	//std::cout<<"before my_free(cpBUSZS)"<<std::endl;
	//my_free(cpBUSZS);
  }
  

  /**
   * MultiUniPlainMod_1:
   * @N: number of variables.
   * @d1: degree of 'FPtr'.
   * @FPtr: a C-Cube polynomial, the divident.
   * @d2: degree of 'GPtr'.
   * @GPtr: a C-Cube polynomial, the divisor.
   * @ts: a triangular set.
   * @tris: modular inverses of 'ts'.
   * @pPtr: the information of prime p.
   * 
   * In formula FPtr = QPtr * GPtr + rPtr,
   * 'FPtr' and 'GPtr' are input. 
   * 'rPtr' is the output remainder. And 'rPtr' uses 'FPtr' as the in-place buffer.
   * Return value: the degree of the 'QPtr'.
   **/
  preFFTRep *
  MultiUniPlainMod_1(sfixn N, sfixn d1, preFFTRep * FPtr, sfixn d2, preFFTRep * GPtr, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr){
	sfixn d, a, b;
	
	preFFTRep co;
	  
	d1=shrinkDeg(BUSZSI(FPtr, N), DAT(FPtr), CUMI(FPtr, N));
	d2=shrinkDeg(BUSZSI(GPtr, N), DAT(GPtr), CUMI(GPtr, N));
	
	d=d1-d2;
	if(d<0) return FPtr;
	if(N==0){
	  a=DATI(FPtr, 0);
	  b=DATI(GPtr, 0);
	  DATI(FPtr, 0)=0;
	  return FPtr;
	}
	
	if(N==1){
	  if((BDSI(ts, 1))>DivCut1){ 
		PBPAS::UniFastMod_1(d1, DAT(FPtr), d2, DAT(GPtr), 
					 (NLB(tris))[1], DAT(ELEMI(tris, 1)) ,pPtr);}
	  else{
		PBPAS::UniPlainMod_1(d1, DAT(FPtr), d2, DAT(GPtr), pPtr);}
	}
	else{
	  SERBPAS::InitOneReducedPoly(&co, N-1, BUSZS(FPtr));

      while(d>=0){
		PBPAS::PolyCleanData(&co);
        PBPAS::getCoefMulti(N, &co, FPtr, d1);
		//std::cout<<"before SERBPAS::fmecgEqDg_1"<<std::endl;
        d1 = PBPAS::fmecgEqDg_1(N, FPtr, d, &co, GPtr, ts, tris, pPtr);
		//std::cout<<"after SERBPAS::fmecgEqDg_1"<<std::endl;	
	
        if(d1==-1) {
		  freePoly(&co); 
		  return FPtr;
		}
		d=d1-d2;
	  }
	  
      freePoly(&co);
	} 
	return FPtr;
	
  }
  
  //========================================================================
  //  f1 = f1 - co * X**e * f2
  // 
  //  
  //========================================================================
  sfixn
  fmecgEqDg_1(sfixn N, preFFTRep * f1, sfixn e, preFFTRep * co, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr){
	//register sfixn i;
	sfixn d1,d2,d;
	//preFFTRep out, res;
	//std::cout<<"enter fmecgEqDg_1"<<std::endl;
	d1=shrinkDeg(BUSZSI(f1, N), DAT(f1), CUMI(f1, N));
	d2=d1-e;
	//std::cout<<"fmecgEqDg_1, d2="<<d2<<std::endl;
	//preFFTRep *out1;	
	//PBPAS::InitOneReducedPoly(out1, N-1, BUSZS(f1));
	//std::cout<<"fmecgEqDg_1, after out1"<<std::endl;

	preFFTRep **out = (preFFTRep **)my_malloc((d2+1)*sizeof(preFFTRep *));
	cilk_for(int i=0; i<=d2; i++){
	  out[i] = (preFFTRep *)my_calloc(1, sizeof(preFFTRep));
	  PBPAS::InitOneReducedPoly(out[i], N-1, BUSZS(f1));
	} 
	//SERBPAS::InitOneReducedPoly(&out, N-1, BUSZS(f1)); 
	//std::cout<<"fmecgEqDg_1, after init out"<<std::endl;

	//backupData(f2);
	//decreaseOneDim(f2);
	
	preFFTRep **res = (preFFTRep **)my_malloc((d2+1)*sizeof(preFFTRep *));
	cilk_for(int i=0; i<=d2; i++){
	  res[i] = (preFFTRep *)my_calloc(1, sizeof(preFFTRep));
	  SERBPAS::InitResPoly(res[i], N-1, BUSZS(co), BUSZS(f2));
	}
	//SERBPAS::InitResPoly(&res, N-1,  BUSZS(co), BUSZS(f2));
	//std::cout<<"fmecgEqDg_1, after init res"<<std::endl;

	//backupData(f1);
	//decreaseOneDim(f1); //SIZ(f1)=CUMI(f1, N(f1)); (N(f1))--;
	//nextMCoefData(f1,N,e); //DAT(f1)+=e*(CUMI(f1, N));

	sfixn * dataf1 = DAT(f1)+e*(CUMI(f1, N));
	//std::cout<<"fmecgEqDg_1, dataf1="<<dataf1<<std::endl;
	preFFTRep ** cpf1 = (preFFTRep **)my_malloc((d2+1)*sizeof(preFFTRep *));
	preFFTRep ** cpf2 = (preFFTRep **)my_malloc((d2+1)*sizeof(preFFTRep *));

	cilk_for(int i=0; i<=d2; i++){
	  //std::cout<<"fmecgEqDg_1, i="<<i<<std::endl;
	  //SERBPAS::PolyCleanData(&out);
	  //PBPAS::PolyCleanData(out[i]);

	  //SERBPAS::PolyCleanData(&res);
	  //PBPAS::PolyCleanData(res[i]);

	  cpf1[i] = SERBPAS::CopyOnePoly_sameData(f1);
	  SIZ(cpf1[i]) = CUMI(f1, N);
	  N(cpf1[i]) = N-1;
	  DAT(cpf1[i]) = dataf1 + i*CUMI(f1, N);

	  cpf2[i] = SERBPAS::CopyOnePoly_sameData(f2);
	  SIZ(cpf2[i]) = CUMI(f2, N);
	  N(cpf2[i]) = N-1;
	  DAT(cpf2[i]) = DAT(f2) + i*CUMI(f2, N);

	  PBPAS::mul_Reduced(res[i], N-1, out[i], co, cpf2[i], ts, tris, pPtr);
	  //SERBPAS::mul_Reduced(&res, N-1, &out, co, f2, ts, tris, pPtr);
	  
	  SERBPAS::subEqDgPoly_1(N-1, cpf1[i], out[i], pPtr->P, 1); //inline
	  //subEqDgPoly_1(N-1, f1, &out, pPtr->P, 1); //inline

	  //nextCoefData(f1, N);   
	  //nextCoefData(f2, N);
	  freePoly(res[i]);
	  freePoly(out[i]);
	}
	
	//cilk_for(int i=0; i<=d2; i++){
	//freePoly(res[i]);
	//freePoly(out[i]);
	//}

	//freePoly(&res);
	//increaseOneDim(f1);
	//increaseOneDim(f2);
	//restoreData(f1);
	//restoreData(f2);
	//freePoly(&out);

	d=shrinkDeg(BUSZSI(f1, N), DAT(f1), CUMI(f1, N));

	if((zeroCoefp(DAT(f1), CUMI(f1,N)))&& (d==0)) return -1; else return d; //inline
  }
  
  //===================================================
  // MulMod
  // (1) resPtr = f1 * f2.
  // (2) out = resPtr modulo TS.
  //  * resPtr is pre-allocated and passed as a parameter.
  //  * resPtr will be destructed during MultiMod operation.
  //===================================================
  /**
   * mul_Reduced:
   * @resPtr: (output) the product: 'f1' * 'f2'.
   * @N: the number of variables.
   * @out: (output) the modular product: 'f1' * 'f2' mod 'ts' .
   * @f1: C-Cube polynomial.
   * @f2: C-Cube polynomial.
   * @ts: the triangular set.
   * @tris: the inverses of 'ts'.
   * @pPtr: the information of the prime.
   * 
   * Return value: 
   **/
  preFFTRep * 
  mul_Reduced(preFFTRep * resPtr, sfixn N, preFFTRep * out, preFFTRep * f1, preFFTRep * f2, TriSet * ts, TriRevInvSet * tris, MONTP_OPT2_AS_GENE * pPtr){
	KroFFTRep kro;
	
	sfixn d1,d2,d3,d4,d5,sz1, sz2;
	d1=shrinkDeg(BUSZSI(f1, N), DAT(f1), CUMI(f1, N));
	d2=shrinkDeg(BUSZSI(f2, N), DAT(f2), CUMI(f2, N));

	sz1=SERBPAS::getDenseSiz(N, f1, BUSZSI(f1, N), DAT(f1), CUMI(f1, N));
	sz2=SERBPAS::getDenseSiz(N, f2, BUSZSI(f2, N), DAT(f2), CUMI(f2, N));

	//d3=BUSZSI(f1, N);
	//d4=BUSZSI(f2, N);
	
	if((sz1>=MulCut)||(sz2>=MulCut)){
	  d5= BUSZSI(resPtr, N);

	  //race here
	  //BUSZSI(f1, N)=d1;
	  //BUSZSI(f2, N)=d2;

	  preFFTRep * cpf1 = SERBPAS::CopyOnePoly_sameData(f1);
	  preFFTRep * cpf2 = SERBPAS::CopyOnePoly_sameData(f2);
	  BUSZSI(cpf1, N)=d1;
	  BUSZSI(cpf2, N)=d2;

	  BUSZSI(resPtr, N)=d1+d2;

	  SERBPAS::InitKroFFTRep(&kro, BUSZS(resPtr), N, 2, pPtr); 
	  //================================================================
	  PBPAS::mulPoly_FFT(N, &kro, resPtr, cpf1, cpf2,  pPtr);
	  //================================================================

	  SERBPAS::freePoly_not_data(cpf1);
	  SERBPAS::freePoly_not_data(cpf2);
	  //race here	  
	  //BUSZSI(f1, N)=d3;
	  //BUSZSI(f2, N)=d4;

	  BUSZSI(resPtr, N)=d5;
	  
	  SERBPAS::freeKroFFTRep(&kro);}
	else{
	  
	  if(N==1){
		PBPAS::EX_Mont_PlainMul_OPT2_AS_GENE(BUSZSI(resPtr, 1), DAT(resPtr), d1, DAT(f1), d2, DAT(f2), pPtr); 
		
	  } else{
		PBPAS::plainMultiDMul(N, CUM(resPtr), DAT(resPtr), CUM(f1), BUSZS(f1), CUM(f2), BUSZS(f2), DAT(f1), DAT(f2), pPtr); 
	  }
	  
	}
	
	PBPAS::MultiMod(N, out, resPtr, ts, tris, pPtr);
	return out; //not necessary --yu
  }
  
  //***************************************************
  //normal form starts
  //***************************************************
  
  //============================================================
  //
  //
  //
  //  Univariate polynomial Fast Mod.
  //
  //
  //=============================================================
  
  // suppose the RevInvPtr is known. or Say the first poly's RevInv in TriSet is computed.
  //destructive APtr into RPtr.
  
  /**
   * UniFastMod_1:
   * @degA: the degree of univariate polynomial 'A'.
   * @APtr: the coefficient vector of polynomial 'A'.
   * @degB: the degree of univariate polynomial 'B'.
   * @BPtr: the coefficient vector of polynomial 'B'.
   * @Lbound: The next power of 2 of (degA-degB) .
   * @BRevInvPtr: the modular inverse of 'B'.
   * @pPtr: the information of prime number p.
   * 
   * Compute the remainder of A divied by B.
   *
   * Return value: 
   **/
  void
  UniFastMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, sfixn Lbound, sfixn * BRevInvPtr,  MONTP_OPT2_AS_GENE * pPtr){
    
	sfixn * FPtr,  *QPtr, * GPtr;
	sfixn degF, power2,power3, n2,n3,l1,l2,l3, tmp, sz;
	sfixn dg, da;
	sfixn degQ=degA-degB;
	
	//std::cout<<"PBPAS::UniFastMod_1"<<std::endl;

	degF=degQ;
	
	if(degF<0) {
	  //std::cout<<"PBPAS::UniFastMod_1, degF<0"<<std::endl;
	  return;}
	
	l1=degF+1;
	
	n2=1; power2=0; l2=(l1<<1)-1;tmp=l2;
	while(tmp){tmp>>=1; n2<<=1; power2++;}
	n3=1; power3=0;  l3=degB+degQ+1, tmp=l3;
	while(tmp){tmp>>=1; n3<<=1; power3++;}
	
	dg = da = degF;
	
	sz = n2;
	if( n3>n2 ) sz=n3;
	
	FPtr=(sfixn * )my_calloc(sz ,sizeof(sfixn));
	GPtr=(sfixn * )my_calloc(sz ,sizeof(sfixn));
	
	
	cilk_for(int i=0; i<=dg; i++) GPtr[i] = BRevInvPtr[i];
	
	FPtr = PBPAS::reverseUni(degB, FPtr, BPtr);
	
	FPtr = PBPAS::reverseUni(degA, FPtr, APtr);
	
	PBPAS::EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
	
	QPtr=(sfixn * ) my_calloc(degQ+1, sizeof(sfixn));
	
	QPtr = PBPAS::reverseUni(degQ, QPtr, FPtr);
	
	PBPAS::cleanVec(n3-1, FPtr);
	
	
	PBPAS::EX_Mont_FFTMul_OPT2_AS_GENE(n3, power3, degQ+degB, FPtr, degQ, QPtr, degB, BPtr, pPtr);
	
	
	cilk_for(int j=0; j<l3; j++) APtr[j] = SubMod(APtr[j], FPtr[j], pPtr->P);
	cilk_for(int j=l3; j<=degA; j++) APtr[j] = 0;
	
	my_free(FPtr);
	my_free(GPtr);
	my_free(QPtr);

	//std::cout<<"end PBPAS::UniFastMod_1"<<std::endl;
  } 

  //============================================
  // Univariate Plain division.
  // type: in-place.
  // output: the remainder.
  //============================================
  
  /** 
   * UniPlainMod_1:
   * @degA: The degree of the univeriate polynomial 'A'.
   * @APtr: The coefficient vector of univeriate polynomial 'A'.
   * @degB: The degree of the univeriate polynomial 'B'.
   * @BPtr: The coefficient vector of univeriate polynomial 'B'.
   * @pPtr: The information of the prime.
   * 
   * Compute the remainder R, where A = QB+R mod p.
   * The content of input 'A' will be replaced by 'R'. 
   * Return value: 
   **/
  void 
  UniPlainMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr ){
	register sfixn i,j,p=pPtr->P;
	
	sfixn tmp;

	//std::cout<<"PBPAS::UniPlainMod_1"<<std::endl;
	
	if((degA-degB)<0) return;
	if(degB>10) {
      for(i=degA-degB; i>=0; i--){ //race for normalform 1 100 3 
        tmp=MontMulMod_OPT2_AS_GENE(APtr[degB+i], pPtr->R2BRsft, pPtr)<<(pPtr->Base_Rpow);
        for(j=0; j<degB; j++) 
		  APtr[i+j]=SubMod(APtr[i+j], MontMulMod_OPT2_AS_GENE(BPtr[j], tmp, pPtr), p);

        APtr[i+degB]=0;
	  }
	}
	else {
	  for(i=degA-degB; i>=0; i--){ 
		for(j=0; j<degB; j++) 
		  APtr[i+j]=SubMod(APtr[i+j], MulMod(BPtr[j], APtr[degB+i], p), p);

	    APtr[i+degB]=0;
	  }
	}
  }
  

  /**
   * EX_Mont_FFTMul_OPT2_AS_GENE_1:
   * @n: the FFT size
   * @e: log2 of n
   * @degRes: degree of result
   * @degA: (Input and Output) degree of A
   * @APtr: coefficient vector of A 
   * @degB: degree of B
   * @BPtr: coefficient vector of B
   * @pPtr: prime number structure
   * 
   * FFT-based multiplication of A by B. In place: result in A.
   * 
   * 
   * Return value: the product of two polynomials.
   **/   
  void EX_Mont_FFTMul_OPT2_AS_GENE_1(sfixn n, sfixn e, sfixn degRes, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr){
	//register int i;
	sfixn l_AB, e1, n1, sA, sB;
	sfixn * rootsPtr;

	//std::cout<<"PBPAS::EX_Mont_FFTMul_OPT2_AS_GENE_1"<<std::endl;

	l_AB=degRes+1; // compute l_AB for TFT(A,B).
	e1=logceiling(l_AB);
	n1=1L<<e1;
	// force this FFT use the n from parameter if n<n1.
	if ( (n<n1) || (degRes==0) ) {n1=n; e1=e;}
	PBPAS::checkFrourierPrime(e1, pPtr);
	
	rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	
	PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(e1, n1, rootsPtr, pPtr);
	
	sA = degA+1;
	sB = degB+1;
	for(int i=sA; i<n1; i++) APtr[i]=0;
	for(int i=sB; i<n1; i++) BPtr[i]=0; 
 
	cilk_spawn EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, APtr, pPtr);
	EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, BPtr, pPtr);
	cilk_sync;

	PBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n1, APtr, BPtr, pPtr);
	
	EX_Mont_INVDFT_OPT2_AS_GENE_R_1(n1, e1, rootsPtr, APtr, pPtr);
	
	my_free(rootsPtr);
	
  }
  
  /**
   * EX_Mont_FFTMul_OPT2_AS_GENE:
   * @n: the FFT size
   * @e: log2 of n
   * @degRes: degree of result
   * @resPtr: coefficient vector of result
   * @degA: degree of A
   * @APtr: coefficient vector of A 
   * @degB: degree of B
   * @BPtr: coefficient vector of B
   * @pPtr: prime number structure
   * 
   * FFT-based multiplication of A by B. Result is in resPtr.
   * 
   * 
   * Return value: the product of two polynomials.
   **/   
  void EX_Mont_FFTMul_OPT2_AS_GENE(sfixn n, sfixn e, sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr){
	//register int i;
	sfixn l_AB, e1, n1, BRsft=pPtr->Base_Rpow;
	sfixn * rootsPtr, * dftAPtr, * dftBPtr;

	//std::cout<<"PBPAS::EX_Mont_FFTMul_OPT2_AS_GENE"<<std::endl;

	l_AB=degRes+1; // compute l_AB for TFT(A,B).
	e1=logceiling(l_AB);
	n1=1L<<e1;
	// force this FFT use the n from parameter if n<n1.
	if ((n<n1) || (degRes==0)) {n1=n; e1=e;}

	PBPAS::checkFrourierPrime(e1, pPtr);
	
	rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	dftAPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	dftBPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));	
	
	if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){	  
	  PBPAS::Mont_GetNthRoots_OPT2_AS_GENE(e1,n1,rootsPtr,pPtr);
	  
	  cilk_spawn MODPN::Mont_dft_OPT2_AS_GENE( n1, e1, rootsPtr, dftAPtr, degA, APtr, pPtr);
	  MODPN::Mont_dft_OPT2_AS_GENE( n1, e1, rootsPtr, dftBPtr, degB, BPtr, pPtr);
	  cilk_sync;

	  for(int i=0; i<n1; i++) dftAPtr[i]=MontMulMod_OPT2_AS_GENE(dftAPtr[i], dftBPtr[i]<<BRsft, pPtr);
	  
	  MODPN::Mont_invdft_OPT2_AS_GENE_R(n1, e1, rootsPtr, dftAPtr, degRes, resPtr, pPtr);
	  
	}
	else{
	  PBPAS::Mont_GetNthRoots_OPT2_AS_GENE_SPE(e1,n1,rootsPtr,pPtr);

	  cilk_spawn MODPN::Mont_dft_OPT2_AS_GENE_SPE( n1, e1, rootsPtr, dftAPtr, degA, APtr, pPtr);
	  MODPN::Mont_dft_OPT2_AS_GENE_SPE( n1, e1, rootsPtr, dftBPtr, degB, BPtr, pPtr);
	  cilk_sync;

	  for(int i=0; i<n1; i++) dftAPtr[i]=MontMulMod_OPT2_AS_GENE_SPE(dftAPtr[i], dftBPtr[i]<<BRsft, pPtr);

	  MODPN::Mont_invdft_OPT2_AS_GENE_SPE_R(n1, e1, rootsPtr, dftAPtr, degRes, resPtr, pPtr);
	  
	}
		
	my_free(dftAPtr);
	my_free(dftBPtr);
	my_free(rootsPtr);
  }
  

  /**
   * EX_mulPoly:
   * @N: Number of variables.
   * @resPtr: (output) The output prod of 'p1Ptr' and 'p2Ptr'.
   * @f1: A C-Cube polynomial.
   * @f2: A C-Cube polynomial.
   * @pPtr: Information for the Fourier prime number p. 
   *
   * Compute the product of f1 and f2  by either Kronecker-based univariate FFT approach or by direct Multi-dimensional FFT or by classical multiplication.
   *
   * Return value:.
   **/
  void EX_mulPoly(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr){
	sfixn sz1, sz2;
	sz1=getDenseSiz(N, f1, BUSZSI(f1, N), DAT(f1), CUMI(f1, N));
	sz2=getDenseSiz(N, f2, BUSZSI(f2, N), DAT(f2), CUMI(f2, N));
	
	if((sz1>=MulCut)||(sz2>=MulCut)){
	  //EX_mulPoly_FFT( N, resPtr, f1, f2,  pPtr);
	  
	  PBPAS::EX_mulPoly_TFTFFT(N, resPtr, f1, f2, pPtr);
	  
	}else{
	  PBPAS::plainMultiDMul(N, CUM(resPtr), DAT(resPtr), CUM(f1), BUSZS(f1), CUM(f2), BUSZS(f2), DAT(f1), DAT(f2), pPtr);
	}
  }
  
  /**
   * EX_mulPoly_TFTFFT:
   * @N: the number of variables.
   * @resPtr: (output) the product of 'f1' and 'f2', the C-Cube polynomial .
   * @f1: a C-Cube polynomial.
   * @f2: a C-Cube polynomial.
   * @pPtr: the information of the prime p.
   * 
   * The product f1 and f2.
   *
   * Return value: 
   **/
  void EX_mulPoly_TFTFFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE * pPtr)
  {
	KroTFTRep kro;

	//magic 8 here ---yu
	if( SERBPAS::forcutoff_Multi_TFT_FFT(N, BUSZS(f1), BUSZS(f2), 8) ==0 ){
	  PBPAS::EX_mulPoly_FFT(N, resPtr, f1, f2, pPtr);
	   
	}
	else{
	  SERBPAS::InitKroTFTRep(&kro, BUSZS(resPtr), N, 2, pPtr); 
	  //================================================================
	  PBPAS::polyMul_TFT(N, &kro, resPtr, f1, f2,  pPtr);
	  //===============================================================
	  SERBPAS::freeKroTFTRep(&kro);
	}
    return;
  }
  
  /**
   * EX_mulPoly_FFT:
   * @N: Number of variables.
   * @resPtr: (output) The output prod of 'p1Ptr' and 'p2Ptr'.
   * @f1: A C-Cube polynomial.
   * @f2: A C-Cube polynomial.
   * @pPtr: Information for the Fourier prime number p. 
   *
   * Compute the product of f1 and f2  by either Kronecker-based 
   * univariate FFT approach or by direct Multi-dimensional FFT.
   *
   * Return value: The product of f1 and f2.
   **/
  void EX_mulPoly_FFT(sfixn N, preFFTRep * resPtr, preFFTRep * f1, preFFTRep * f2,  MONTP_OPT2_AS_GENE *pPtr){
    KroFFTRep kro;
	
    SERBPAS::InitKroFFTRep(&kro, BUSZS(resPtr), N, 2, pPtr);
    //================================================================
    PBPAS::mulPoly_FFT(N, &kro, resPtr, f1, f2,  pPtr);
    //===============================================================
    SERBPAS::freeKroFFTRep(&kro);
    return;
  }
  
  //============================================================
  // rPtr = p1Ptr * p2Ptr
  //============================================================
  
  /**
   * mulPoly_FFT:
   * @N: Number of variables.
   * @kPtr: The Kronecker data structure.
   * @rPtr: (output) The output prod of 'p1Ptr' and 'p2Ptr'.
   * @p1Ptr: A C-Cube polynomial.
   * @p2Ptr: A C-Cube polynomial.
   * @pPtr: Information for the Fourier prime number p. 
   *
   * Compute the product of p1Ptr and p2Ptr  by either Kronecker-based 
   * univariate FFT approach or by direct Multi-dimensional FFT.
   *
   * Return value: 
   **/
  void mulPoly_FFT(sfixn N, KroFFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr){
    sfixn kdg1=0, kdg2=0;
    int switcher=0;	
	
    if((KE(kPtr))>(pPtr->Npow)) 
      { 
		//printf("FFT size larger than the prime %ld can handle!  using MultiDFFT to solve the problem.  -> ", pPtr->P);
		switcher=1;}
	
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(p1Ptr), BUSZS(p1Ptr), DAT(p1Ptr));
    PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(p2Ptr), BUSZS(p2Ptr), DAT(p2Ptr));
	
    //if((switcher) || (KN(kPtr)<4)){
    //magic number 4 here, for debugging when size of fft is very small ---yu
	if( switcher) {
	  //================================================================
	  PBPAS::fftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), pPtr);
	  //SERBPAS::fftMultiD_test(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), pPtr);
	  //================================================================
	}
    else{
      //KN FFT
	  kdg1= (BUSZSI(p1Ptr, N)+1)*CUMI(kPtr, N)-1;
	  kdg2= (BUSZSI(p2Ptr, N)+1)*CUMI(kPtr, N)-1;
	  KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));
	  
	  SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr),pPtr);
	  
	  //================================================================
	  PBPAS::EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), kdg2, DATSI(kPtr, 1), pPtr);
	  //================================================================
	  
    }

    PBPAS::fromtofftRepMultiD(N,  CUM(rPtr), DAT(rPtr), CUM(kPtr), BUSZS(rPtr), DATSI(kPtr, 0));

  }
  
  
  //-----------------------------------------------
  void polyMul_TFT(sfixn N, KroTFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr){
	
	PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(p1Ptr), BUSZS(p1Ptr), DAT(p1Ptr));
	
	PBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(p2Ptr), BUSZS(p2Ptr), DAT(p2Ptr));
	
	//================================================================
	PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	//SERBPAS::tftMultiD_test(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	//================================================================
	
	PBPAS::fromtofftRepMultiD(N,  CUM(rPtr), DAT(rPtr), CUM(kPtr), BUSZS(rPtr), DATSI(kPtr, 0));
  }
  
  //=========================================================================
  // Bottom-up level-by-level reduction.
  // using TFTs.
  //=========================================================================
  void
  MultiMod_BULL(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr)
  { 
	//std::cout<<"PBPAS::MultiMod_BULL"<<std::endl;
	if (zeroCoefp(DAT(inPtr), SIZ(inPtr)) ) 
	  return;
	
	//=========================================
	PBPAS::MultiMod_1_BULL_par(N, inPtr, ts, tris, pPtr);
	//=========================================
	
	PBPAS::fromtofftRep(N, CUM(outPtr), DAT(outPtr), CUM(inPtr),  BUSZS(outPtr), DAT(inPtr));
	
  }
  
  /**
   * MultiMod_1_BULL_par:
   * @N: number of variables.
   * @outPtr: (output) the image of 'inPtr' modulo 'ts'.
   * @inPtr: a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverse of 'ts'.
   * @pPtr: the information of prime number p.
   * 
   * Reduce 'inPtr' w.r.t 'ts' and save the result in 'outPtr'.
   *
   * Return value: 
   **/
  void MultiMod_1_BULL_par(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	
	sfixn *lastNonZeroPtr, *datPtr;// *tmpVec;
	sfixn siz; //deg;
	//preFFTRep * tmpInPtr, *cpShellInPtr;
	
	//std::cout<<"PBPAS::MultiMod_BULL_1"<<std::endl;
	if(N<1) { std::cout<<"Number of variables < 1."<<std::endl; exit(1); }
	
	siz=BUSZSI(inPtr, 1)+1;
	
	lastNonZeroPtr = DAT(inPtr) + shrinkDeg(SIZ(inPtr)-1, DAT(inPtr), 1);  
	
	//for(datPtr=DAT(inPtr); datPtr<=lastNonZeroPtr; datPtr+=siz){
	sfixn deg2 = (lastNonZeroPtr - DAT(inPtr))/siz;
	cilk_for (int k=0; k<=deg2; k++) {
	  sfixn *datPtrk = DAT(inPtr) + k*siz;
	  sfixn degk=shrinkDeg(BUSZSI(inPtr, 1), datPtrk, 1);
	  
	  if (((BDSI(ts, 1))>DivCut1)) {
		PBPAS::UniFastMod_1(degk, datPtrk, BUSZSI((ELEMI(ts, 1)), 1), DAT(ELEMI(ts, 1)), 
					 (NLB(tris))[1], DAT(ELEMI(tris, 1)), pPtr); 
	  }
	  else{
		PBPAS::UniPlainMod_1(degk, datPtrk, BUSZSI((ELEMI(ts, 1)), 1), 
					  DAT(ELEMI(ts, 1)), pPtr); 
	  }
	}
		
	if(N==1)  return;
	
	//tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	//tmpVec=(sfixn *)my_calloc(N+1, sizeof(sfixn));
	
	for(int i=2; i<N-1; i++){  

	  //for(datPtr=DAT(inPtr); datPtr<=lastNonZeroPtr;
	  //datPtr+=CUMI(inPtr, i+1)){
	  sfixn degiplus1 = (lastNonZeroPtr - DAT(inPtr))/CUMI(inPtr, i+1);
	  cilk_for (int k=0; k<=degiplus1; k++) {
		sfixn *datPtrik = DAT(inPtr) + k*CUMI(inPtr, i+1);
		sfixn degik = shrinkDeg(CUTSI(inPtr, i), datPtrik, CUMI(inPtr, i));    
		
		preFFTRep * tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
		sfixn *tmpVec=(sfixn *)my_calloc(N+1, sizeof(sfixn));
		
		for(int j=0; j<i; j++) 
		  tmpVec[j]=BDSI(ts, j);
		
		tmpVec[i]=degik;
        
		SERBPAS::InitOneReducedPoly(tmpInPtr, i, tmpVec);
		PBPAS::fromtofftRep(i, tmpInPtr->accum, tmpInPtr->data, CUM(inPtr), 
							  BUSZS(tmpInPtr), datPtrik);
		
		preFFTRep *cpShellInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
		
		SERBPAS::copyPolyPointers(cpShellInPtr, inPtr);
		
		N(cpShellInPtr)=i;
		SIZ(cpShellInPtr)=CUMI(inPtr, i+1);
		DAT(cpShellInPtr)=datPtrik;
		
		if((((i==2)&&(BDSI(ts, i)>DivCut2)) || ((i==3)&&(BDSI(ts, i)>DivCut3)) || ((i>3)&&(BDSI(ts, i)>DivCutN)))){		  
		  PBPAS::MultiUniFastMod_1_TFT(i, tmpInPtr, degik, cpShellInPtr, ts, tris, pPtr); 
		  //why just use TFT here? ---yu
															  
		}
		else{			
		  PBPAS::MultiUniPlainMod_1(i, degik, tmpInPtr, BDSI(ts, i)+1,  ELEMI(ts, i), ts, tris, pPtr);
		  
		  PBPAS::fromtofftRep(i, CUM(cpShellInPtr), DAT(cpShellInPtr), CUM(tmpInPtr),  BUSZS(tmpInPtr), DAT(tmpInPtr));
		  
		}
		
		freePoly(tmpInPtr);
		my_free(tmpInPtr);
		my_free(tmpVec);
		my_free(cpShellInPtr);
	  }
	  
	}
	
	//my_free(tmpVec);
	//my_free(tmpInPtr);    
	
	if(N>=3){
	  
	  //tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	  //tmpVec=(sfixn *)my_calloc(N, sizeof(sfixn));

	  //for(datPtr=DAT(inPtr); datPtr<=lastNonZeroPtr; datPtr+=CUMI(inPtr, N)){
	  sfixn degNplus1 = (lastNonZeroPtr - DAT(inPtr))/CUMI(inPtr, N);
	  cilk_for (int k=0; k<=degNplus1; k++) {
		sfixn *datPtrNk = DAT(inPtr) + k*CUMI(inPtr, N);
		sfixn degNk = shrinkDeg(CUTSI(inPtr, N-1), datPtrNk, CUMI(inPtr, N-1));
 
		preFFTRep * tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
		sfixn *tmpVec=(sfixn *)my_calloc(N+1, sizeof(sfixn));
     
		for(int j=0; j<N-1; j++) 
		  tmpVec[j]=BDSI(ts, j);
		
		tmpVec[N-1]=degNk;
		SERBPAS::InitOneReducedPoly(tmpInPtr, N-1, tmpVec);
		PBPAS::fromtofftRep(N-1, tmpInPtr->accum, tmpInPtr->data, CUM(inPtr), 
					 BUSZS(tmpInPtr), datPtrNk);

		preFFTRep *cpShellInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
		SERBPAS::copyPolyPointers(cpShellInPtr, inPtr);

		N(cpShellInPtr)=N-1;
		SIZ(cpShellInPtr)=CUMI(inPtr, N);
		DAT(cpShellInPtr)=datPtrNk;
		
		if((((N-1==2)&&(BDSI(ts, N-1)>DivCut2)) || ((N-1==3)&&(BDSI(ts, N-1)>DivCut3)) || ((N-1>3)&&(BDSI(ts, N-1)>DivCutN)))){
		  
		  PBPAS::MultiUniFastMod_1_TFT(N-1, tmpInPtr, degNk, cpShellInPtr, ts, tris, pPtr);
		}
		else{
		  
		  PBPAS::MultiUniPlainMod_1(N-1, degNk, tmpInPtr, BDSI(ts, N-1)+1,  ELEMI(ts, N-1), ts, tris, pPtr);
		  PBPAS::fromtofftRep(N-1, CUM(cpShellInPtr), DAT(cpShellInPtr), CUM(tmpInPtr),  BUSZS(tmpInPtr), DAT(tmpInPtr));
		}
		
		my_free(tmpVec);
		freePoly(tmpInPtr);
		my_free(tmpInPtr);
		my_free(cpShellInPtr);
	  }
	  //my_free(tmpVec);
	  //my_free(tmpInPtr);
	  
	}
	
	preFFTRep * tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	sfixn *tmpVec=(sfixn *)my_calloc(N+1, sizeof(sfixn));
	sfixn degN=shrinkDeg(BUSZSI(inPtr, N), DAT(inPtr), CUMI(inPtr, N));

	for(int j=0; j<N; j++) 
	  tmpVec[j]=BDSI(ts, j);

	tmpVec[N]=degN;
	SERBPAS::InitOneReducedPoly(tmpInPtr, N, tmpVec);
	PBPAS::fromtofftRep(N, tmpInPtr->accum, tmpInPtr->data, CUM(inPtr), 
				 BUSZS(tmpInPtr), DAT(inPtr));
	PBPAS::MultiUniFastMod_1_TFT(N, tmpInPtr, degN, inPtr, ts, tris, pPtr);
	
	freePoly(tmpInPtr);
	my_free(tmpVec);
	my_free(tmpInPtr);
	
  }
  
  //=============================================================
  // Modulo poly by T_n, then by [T_{n-1}, T_{n-2}, ...,T_1] 
  /**
   * MultiUniFastMod_1_TFT:
   * @N: number of variables.
   * @tmpIn: a temp buffer.
   * @n: The next power of 2 of (deg(inPtr, X_N)-deg(T_N, X_N)+1).
   * @inPtr:(output) a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverses of 'ts'.
   * @pPtr: the information of the prime p.
   * 
   * Return value: 
   **/
  void MultiUniFastMod_1_TFT(sfixn N, preFFTRep * tmpIn,  sfixn n, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	sfixn m,e,d1,d2;
	sfixn cutDg;

	preFFTRep * resPoly, * resPoly2, * out;
	KroTFTRep * krotftrep, * krotftrep2;
	
	
	m=BUSZSI((ELEMI(ts, N)), N);
	
	if(n<m){ if(DEBUG2) printf("Leaving MultiUniFastMod_1.\n"); 
	  return;}
	
	resPoly=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	resPoly2=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	out=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	
	krotftrep=(KroTFTRep *)my_calloc(1,sizeof(KroTFTRep));
	krotftrep2=(KroTFTRep *)my_calloc(1,sizeof(KroTFTRep));
	e=logceiling(n-m+1);

	if(! (EXSTI(tris, N))){  
	  SERBPAS::initCopyOneRevInv(out, ELEMI(tris, N));
	  SERBPAS::NewtonRevInverse(N, out, ELEMI(tris, N), ELEMI(ts, N), ts, tris, pPtr );
	  freePoly(out);
	}
	
	cutDg=n-m;
	d1=BUSZSI(tmpIn, N);
	//d2=BUSZSI(ELEMI(tris, N), N);
	BUSZSI(tmpIn, N)=cutDg;
	//BUSZSI(ELEMI(tris, N), N)=cutDg;

	//avoid from writing to tris
	preFFTRep * cpNtris = SERBPAS::CopyOnePoly_sameData(ELEMI(tris, N));
	BUSZSI(cpNtris, N) = cutDg;

	PBPAS::reverseMulti_1(n, CUMI(tmpIn, N), DAT(tmpIn));
	//SERBPAS::InitResPoly(resPoly, N, BUSZS(tmpIn), BUSZS(ELEMI(tris, N)));
	SERBPAS::InitResPoly(resPoly, N, BUSZS(tmpIn), BUSZS(cpNtris));

	SERBPAS::InitKroTFTRep(krotftrep, BUSZS(resPoly), N, 2, pPtr);

	//=====================================================================
	//SERBPAS::polyMul_TFT(N, krotftrep, resPoly, tmpIn, cpNtris, pPtr);
	PBPAS::polyMul_TFT(N, krotftrep, resPoly, tmpIn, cpNtris, pPtr);
	SERBPAS::freePoly_not_data(cpNtris);
	//============================================================

	SERBPAS::freeKroTFTRep(krotftrep);
	
	PBPAS::reverseMulti_1(cutDg, CUMI(resPoly, N), DAT(resPoly));
	
	BUSZSI(tmpIn, N)=d1;
	//BUSZSI(ELEMI(tris, N), N)=d2;
	
	PBPAS::PolyCleanData(tmpIn);
	
	BUSZSI(tmpIn, N)=cutDg;
	
	//SERBPAS::reduceCoeffs(N, tmpIn, resPoly, ts, tris, pPtr);
	PBPAS::reduceCoeffs(N, tmpIn, resPoly, ts, tris, pPtr);

	freePoly(resPoly); 
	
	SERBPAS::InitResPoly(resPoly2, N,  BUSZS(ELEMI(ts, N)), BUSZS(tmpIn));
	SERBPAS::InitKroTFTRep(krotftrep2, BUSZS(resPoly2), N, 2, pPtr); 
	
	//if(N==topN) t=gettime(); 
	//================================================================
	//SERBPAS::polyMul_TFT(N, krotftrep2, resPoly2, ELEMI(ts, N), tmpIn,  pPtr);
	PBPAS::polyMul_TFT(N, krotftrep2, resPoly2, ELEMI(ts, N), tmpIn,  pPtr);
	//================================================================
	//if(N==topN) top2MulTime1+=gettime()-t;
	
	SERBPAS::freeKroTFTRep(krotftrep2);
	BUSZSI(tmpIn, N)=d1;
	
	PBPAS::PolyCleanData(tmpIn);
	
	// PureFast_
	
	//SERBPAS::reduceCoeffs(N, tmpIn, resPoly2, ts, tris, pPtr);
	PBPAS::reduceCoeffs(N, tmpIn, resPoly2, ts, tris, pPtr);

	freePoly(resPoly2);
	
	PBPAS::subPoly_1(N, inPtr, tmpIn, pPtr->P);
	
	my_free(resPoly); 
	my_free(resPoly2);
	my_free(out);
	my_free(krotftrep); 
	my_free(krotftrep2);
  }
  
  // type: exported
  // input: A = (degA, APtr), B = (degB, BPtr) and p.
  // ouput: (degRes, resPtr) = A*B mod p.
  // use: classical univariate polynomial multiplication. 
  // note: good for all machine word Fourier Prime.
  /**
   * EX_Mont_PlainMul_OPT2_AS_GENE:
   * @degRes: degree of result
   * @resPtr: coefficient vector of result 
   * @degA:  degree of arg 1
   * @APtr:  coefficient vector of  arg 1
   * @degB: degree of arg 2
   * @BPtr: coefficient vector of  arg 2
   * @pPtr: info about the prime number 'p'
   * 
   * Claasical univariate multiplication 
   * in prime characterisitc 'p' (machine word)
   * using the improved Montgommery trick.
   * Will run special if p-1 is a power of 2.
   * 
   * Return value: 
   **/
  void 
  EX_Mont_PlainMul_OPT2_AS_GENE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr ){
	sfixn d1=degA, d2=degB;
	d1=shrinkDegUni(d1, APtr); //inline
	d2=shrinkDegUni(d2, BPtr);

	if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){
	  SERBPAS::Mont_PlainMul_OPT2_AS_GENE(degRes, resPtr, d1, APtr, d2, BPtr, pPtr );
	}
	else{
	  SERBPAS::Mont_PlainMul_OPT2_AS_GENE_SPE(degRes, resPtr, d1, APtr, d2, BPtr, pPtr );   
	}
  }
  
  //======================================================================
  // Plain univariate polynomial multiplication over Z/pZ.
  //======================================================================
  
  // type: local
  // note: good for all machine word Fourier Prime.
  /**
   * Mont_PlainMul_OPT2_AS_GENE:
   * @degRes: degree of result
   * @resPtr: coefficient vector of result 
   * @degA:  degree of arg 1
   * @APtr:  coefficient vector of  arg 1
   * @degB: degree of arg 2
   * @BPtr: coefficient vector of  arg 2
   * @pPtr: info about the prime number 'p'
   * 
   * Classical univariate multiplication 
   * in prime characterisitc 'p' (machine word)
   * using the improved Montgommery trick.
   * 
   * Return value: 
   **/
  void
  Mont_PlainMul_OPT2_AS_GENE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr ){
    int i,j;
	
    sfixn p=pPtr->P, R=(1L<<pPtr->Rpow)%p, BRsft= pPtr->Base_Rpow,
	  R2=(MulMod(R,R,p))<<BRsft, tmp; //inline

    for(i=0; i<=degA; i++){
      if(! APtr[i]) continue;
      if(APtr[i]==1){
		for(j=0; j<=degB; j++) 
		  resPtr[i+j] = AddMod(BPtr[j], resPtr[i+j], p); //inline
      }
      else{		
        tmp = (MontMulMod_OPT2_AS_GENE(APtr[i],R2,pPtr))<<BRsft; //inline
		for(j=0; j<=degB; j++){
		  resPtr[i+j] = AddMod( MontMulMod_OPT2_AS_GENE(BPtr[j],tmp,pPtr), resPtr[i+j], p); //inline
		  
		}
      }	  
    }
  }
  
  
  // type: local
  // note: for machine word Fourier Prime in special shape.
  /**
   * Mont_PlainMul_OPT2_AS_GENE_SPE:
   * @degRes: degree of result
   * @resPtr: coefficient vector of result 
   * @degA:  degree of arg 1
   * @APtr:  coefficient vector of  arg 1
   * @degB: degree of arg 2
   * @BPtr: coefficient vector of  arg 2
   * @pPtr: info about the prime number 'p'
   * 
   * Claasical univariate multiplication 
   * in prime characterisitc 'p' (machine word)
   * using the improved Montgommery trick.
   * Assume p = c*N+1. Both c,N, need to be power of 2. 
   * 
   * Return value: 
   **/
  void
  Mont_PlainMul_OPT2_AS_GENE_SPE(sfixn degRes, sfixn * resPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr ){
    int i,j;
    sfixn p=pPtr->P, R=(1L<<pPtr->Rpow)%p, BRsft= pPtr->Base_Rpow,
	  R2=(MulMod(R,R,p))<<BRsft, tmp;

    for(i=0; i<=degA; i++){
      if(! APtr[i]) continue;
      if(APtr[i]==1){
        for(j=0; j<=degB; j++) resPtr[i+j]=AddMod(BPtr[j], resPtr[i+j], p);
      }
      else{
        tmp=(MontMulMod_OPT2_AS_GENE_SPE(APtr[i],R2,pPtr))<<BRsft;
		for(j=0; j<=degB; j++){
		  resPtr[i+j]=AddMod( MontMulMod_OPT2_AS_GENE_SPE(BPtr[j],tmp,pPtr), resPtr[i+j], p);
		}
	  }
	}
  }
  
  // 1 -> DF
  // 0 -> BULL
  /**
   * MultiMod_OPT:
   * @N: number of variables.
   * @outPtr: (output) the image of 'inPtr' modulo 'ts'.
   * @inPtr: a C-Cube polynomial.
   * @ts: a triangular set.
   * @tris: modular inverse of 'ts'.
   * @pPtr: the information of prime number p.
   * @opt: switching to 1 using Depth-first traversal.
   *       switching to 0 using Level-order traversal.
   *
   * Reduce 'inPtr' w.r.t 'ts' and save the result in 'outPtr'.
   *
   * Return value: 
   **/
  void
  MultiMod_OPT(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr, int opt){
	
	if(zeroPolyp(inPtr)) return;
	
	if(opt){
	  PBPAS::MultiMod_DF(N, outPtr, inPtr, ts, tris, pPtr);
	}else
	  {
	    PBPAS::MultiMod_BULL(N, outPtr, inPtr, ts, tris, pPtr);
		
	  }
  }
  
  /**
   * EX_getRevInvTriSet:
   * @N: the number of variables.
   * @ts: a triangualr set.
   * @pPtr: the information the prime.
   * 
   * To compute the modular inverses of 'ts'.
   *
   * Return value: the modular inverses of 'ts'.
   **/
  TriRevInvSet *
  EX_getRevInvTriSet(sfixn N,  TriSet * ts,  MONTP_OPT2_AS_GENE * pPtr){
	int i;
	TriRevInvSet *tris;
	sfixn *dgs;
	dgs=(sfixn *)my_calloc(N+1, sizeof(sfixn));
	for(i=1; i<=N; i++){
	  dgs[i] = BDSI(ts, i) + 1;
	}
	tris=(TriRevInvSet *)my_malloc(sizeof(TriRevInvSet));
	PBPAS::initTriRevInvSet( dgs, N, tris, ts);
	PBPAS::getRevInvTiSet(dgs, N, tris, ts, pPtr);
	my_free(dgs);
	return tris;
	
  }

  /*--------------------------------------------------------------------------------------
   N=K*M
   p=RecoveringModulus(sfixn K, sfixn M)
   where K=2^k, 2^+`n>=2+k+log_2(d)+2M and p>=2^n+1
   A(x^(K-1),y^(d1-1))
   B(x^(K-1),y^(d2-1))
   because of theta (for negcyclic converlution), K is a power of 2

   A and B are bivariate, encoded in dense recursive representation 
   for y > x, so A and B appear as univariate polynomials 
   with monomials (y^0,...,y^d1) and coefficients in that are univariate in x
   Assumptions:
   deg(A,x) = deg(B,x) = K-1  (may be thanks to padding either A or B)
   deg(A,y) = d1 - 1
   deg(B,y) = d2 - 1  

   The output is A*B mod x^K+1 mod p 
  ----------------------------------------------------------------------------------------*/
  sfixn * IntPolyMulMod_2DTFT(sfixn *A, sfixn *B, 
			      sfixn d1, sfixn d2, 
			      sfixn K, sfixn p)
  {//--------------------------------------------------------------------------------------
    MONTP_OPT2_AS_GENE *pPtr = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

    //for A*B mod x^K+1 mod p 
    sfixn es1 = logceiling(K); //K=2^k=dims1=ls1      
    sfixn es2K = es1+1;
    sfixn ls2 = d1+d2-1; // y-size of A*B mod x^K + 1
    sfixn es2 = logceiling(ls2);
    
    if ( (es2K > (pPtr->Npow)) || (es2 > (pPtr->Npow)) || (es1 > (pPtr->Npow)) )
      exit (1); //checkFrourierPrime
    
    sfixn dims2 = 1<<es2; //dimension (power of 2) of y //y-size of the 2-D FFT of the result
    sfixn s = K*ls2;      //size of result // real size of the result
    sfixn s1 = K*d1;  //size of A
    sfixn s2 = K*d2;  //size of B
    sfixn K2 = K<<1;      //length of theta vector
  
    //for FFT (omega and gamma), theta^(2K-1), (theta^(-1))^0=1, AP,BP,res   	
    sfixn *rootsPtr = (sfixn *) my_calloc(K+dims2+K2+1+s1+s2+s, sizeof(sfixn));
    // first K slots for the successive powers of omega (x-direction)
    // next dims2 slots for the successive powers of gamma  (y-direction)
    // next (K2+1) slots for the successive powers of theta and 1 
    // next s1 slots for the weigth vwctor A'
    // next s2 slots for the weigth vector B'
    // next s slots are for intermediate computations of the result
    sfixn *thetaPtr = rootsPtr+K+dims2;
    sfixn *rootsPtr_2 = rootsPtr+K;
    //gamma
    cilk_spawn PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2, dims2, rootsPtr_2, pPtr);
    
    // theta^0, ... theta^(2*K-1)    
    PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es2K, K2, thetaPtr, pPtr);
    
    //omega
    for (int i=0; i<K; i++){ //use cilk_for or not? test to see
      rootsPtr[i] = thetaPtr[i<<1];
    }

    cilk_sync;

    //debuging print
    //std::cout<<"IntPolyMulMod_2DTFT pPtr->P = "<<pPtr->P<<std::endl;

    //std::cout<<"powers of omega -------"<<std::endl;
    //for (int i=0; i<K; ++i)
    //  std::cout<<rootsPtr[i]<<", ";
    //std::cout<<std::endl;

    //std::cout<<"powers of gamma -------"<<std::endl;
    //for (int i=0; i<dims2; ++i)
    //  std::cout<<rootsPtr_2[i]<<", ";
    //std::cout<<std::endl;
    
    //apply WeightVector and evaluate
    sfixn * AP = thetaPtr+K2+1;
    sfixn * BP = AP+s1;   
    sfixn * res = BP+s2;
    
    cilk_spawn PBPAS::WeightVectorEval(thetaPtr, A, AP, res, K, d1, dims2, ls2, rootsPtr, pPtr);
    
    sfixn * res2 = (sfixn *) my_calloc(s, sizeof(sfixn)); ; 
    PBPAS::WeightVectorEval(thetaPtr, B, BP, res2, K, d2, dims2, ls2, rootsPtr, pPtr);
    cilk_sync;

    PBPAS::EX_Mont_PairwiseMul_OPT2_AS(s, res2, res, pPtr->P);

    PBPAS::bivarInterpolBy2DTFT(res2, K, dims2, K, ls2, rootsPtr, pPtr);

    sfixn R=(1L<<pPtr->Rpow)%pPtr->P;
    sfixn BRsft=pPtr->Base_Rpow;
    thetaPtr[K2]=R<<BRsft; //(theta^(-1))^0=1, special Mongomery rep    
 
    //std::cout<<"powers of theta -------"<<std::endl;
    //for (int i=0; i<=K2; ++i)
    //  std::cout<<thetaPtr[i]<<", ";
    //std::cout<<std::endl;

    //UnweightVector
    cilk_for(int i=0; i<ls2; i++){
      int m = i*K;
      for(int j=0; j<K; j++){
	int n = m+j;	
	res2[n]=MontMulMod_OPT2_AS_GENE(res2[n], thetaPtr[K2-j], pPtr); 
	//the ith-powers of theta's inverse is backward from thetaPtr[K2] to [K2-K+1]
      }
    }

    my_free(pPtr);
    my_free(rootsPtr);

    return res2;
  }

  // A: matrix of K*d
  // thetaPtr: theta^0, ... theta^(2*K-1)
  // ls2   = d1+d2-1;
  // es2   = logceiling(ls2);
  // dims2 = 1<<es2;
  void WeightVectorEval(sfixn *thetaPtr, sfixn *A, sfixn *AP, sfixn *res, 		
			sfixn K, sfixn d, sfixn dims2, sfixn ls2, 
			sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr)
  {//----------------------------------------------------------------------
    //change coordiate
    cilk_for(int i=0; i<d; i++){
      int m = i*K;
      for(int j=0; j<K; j++){
	int n = m+j;	
	AP[n]=MontMulMod_OPT2_AS_GENE(A[n], thetaPtr[j], pPtr); 
      }
    }
    PBPAS::bivarEvalBy2DTFT(res, K, dims2, K, ls2, AP, K-1, d-1, rootsPtr, pPtr);  
  }
}//end of PBPAS


