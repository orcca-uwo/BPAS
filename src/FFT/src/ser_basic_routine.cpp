// -*- C++ -*- 
// @author Yuzhen Xie

#include <cilk/cilk.h>

#include "../../../include/FFT/src/modpn.h"
#include "../../../include/FFT/src/modpn_export.h"
//#include "general_routine.h"
//#include "basic_routine.h"
#include "../../../include/FFT/src/ser_general_routine.h"
#include "../../../include/FFT/src/ser_basic_routine.h"
#include "../../../include/FFT/src/include/example_util_gettime.h"
#include <iostream>

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


namespace SERBPAS {

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
    SERBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif
      
    //---------------------------
    SERBPAS::tftMultiD_test(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
    
	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );
#ifdef TIMEFN
	  start = example_get_time();
#endif				
    SERBPAS::fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    
    SERBPAS::freeKroTFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }

  /**
   * tftMultiD_test:
   * @coeffs1: coefficient vector for 'f1'.
   * @coeffs2: coefficient vector for 'f2'.
   * @N: number of variables.
   * @es: 2^es[i] = dims [i] for i = 1..n.
   * @dims: the FFT sizes on each dimension, dims[i] is the FFT on
   *        dimensional i, i=1..N.
   * @pPtr: the information of the prime number.
   * 
   * Return value: 
   **/
  void tftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr)
  { sfixn i,j,k,  m=0, n=1, tmp, maxdim=0;
    sfixn * rootsPtr, * tmprootsPtr, * tmpVec;

#ifdef TIMEFN
	float f1tran_time = 0.0;
	float f2tran_time = 0.0;
	float btran_time = 0.0;
	float f1tft_time = 0.0;
	float f2tft_time = 0.0;
	float pairwisemul_time = 0.0;
	float invtft_time = 0.0;
	long start, end;	
#endif
     
	//check Frourier prime
	for(i=1; i<=N; i++) {
	  SERBPAS::checkFrourierPrime(es[i], pPtr);
	}	

    for(i=1; i<=N; i++) {n*=ls[i]; 
      m+=dims[i]; 
      if(dims[i]>maxdim) maxdim=dims[i];
    } 
    rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
    tmprootsPtr=rootsPtr;
    
    tmpVec=(sfixn *)my_calloc(maxdim, sizeof(sfixn));

    // Multi-DFT
    if( es[1] ){
      //getNthRoots(es[1], dims[1], tmprootsPtr, pPtr);
      SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[1], dims[1], tmprootsPtr,  pPtr);

     
      //tmpVec=(sfixn *)my_calloc(dims[1], sizeof(sfixn));
      for(j=0; j<n; j+=ls[1]){

#ifdef TIMEFN	  
	  start = example_get_time();
#endif		
		copyVec_0_to_d(ls[1]-1, tmpVec, coeffs1+j);

		EX_Mont_TDFT_OPT2_AS_GENE_1 ( ls[1], tmprootsPtr, tmpVec, pPtr);

		copyVec_0_to_d(ls[1]-1, coeffs1+j, tmpVec);
#ifdef TIMEFN  
	  end = example_get_time();
	  f1tft_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
		copyVec_0_to_d(ls[1]-1, tmpVec, coeffs2+j);

		for(k=ls[1]; k<dims[1]; k++) tmpVec[k]=0; 
      
		EX_Mont_TDFT_OPT2_AS_GENE_1 ( ls[1], tmprootsPtr, tmpVec, pPtr);

		copyVec_0_to_d(ls[1]-1, coeffs2+j, tmpVec);

		for(k=ls[1]; k<dims[1]; k++) tmpVec[k]=0; 
#ifdef TIMEFN  
	  end = example_get_time();
	  f2tft_time += (end - start) / 1.f;
#endif
      }
      //my_free(tmpVec);      
    }
    
    for(i=2;i<=N;i++){
      // last dim -> [1]
      
      //if(Interrupted==1) {
      //my_free(rootsPtr);
      //return; 
      //}
      
      tmprootsPtr+=dims[1];

#ifdef TIMEFN
	  start = example_get_time();
#endif	  
      SERBPAS::multi_mat_transpose (N, n, i, ls, coeffs1);
#ifdef TIMEFN  
	  end = example_get_time();
	  f1tran_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
      SERBPAS::multi_mat_transpose (N, n, i, ls, coeffs2);  
#ifdef TIMEFN  
	  end = example_get_time();
	  f2tran_time += (end - start) / 1.f;
#endif		  
	  
      if(es[i]){
		//getNthRoots(es[i], dims[i], tmprootsPtr, pPtr);
		SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr, pPtr);
	    
		//tmpVec=(sfixn *)my_calloc(maxdim, sizeof(sfixn));
		cleanVec(maxdim-1, tmpVec);

		//use only one tmpVec
		for(j=0; j<n; j+=ls[i]){

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
		  copyVec_0_to_d(ls[i]-1, tmpVec, coeffs1+j); 

		  EX_Mont_TDFT_OPT2_AS_GENE_1 (ls[i], tmprootsPtr, tmpVec, pPtr);

		  copyVec_0_to_d(ls[i]-1, coeffs1+j, tmpVec);
#ifdef TIMEFN  
	  end = example_get_time();
	  f1tft_time += (end - start) / 1.f;
#endif


#ifdef TIMEFN	  
	  start = example_get_time();
#endif
		  copyVec_0_to_d(ls[i]-1, tmpVec, coeffs2+j);

		  for(k=ls[i]; k<dims[i]; k++) tmpVec[k]=0;  
		  
		  EX_Mont_TDFT_OPT2_AS_GENE_1 (ls[i], tmprootsPtr, tmpVec, pPtr);
		  
		  copyVec_0_to_d(ls[i]-1, coeffs2+j, tmpVec); 
		  
		  for(k=ls[i]; k<dims[i]; k++) tmpVec[k]=0;
#ifdef TIMEFN  
	  end = example_get_time();
	  f2tft_time += (end - start) / 1.f;
#endif
		}
		//my_free(tmpVec);
      }
      
      //  multi_mat_transpose didn't change dims accordings, so we need to do this.
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
	  start = example_get_time();
#endif    
    // Pairwise-Mul
    //for(i=0; i<n; i++) coeffs1[i]=MulMod(coeffs1[i], coeffs2[i],pPtr->P);
    SERBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs2, pPtr->P);
#ifdef TIMEFN  
	  end = example_get_time();
	  pairwisemul_time += (end - start) / 1.f;
#endif
    
    // Multi-invdft
#ifdef TIMEFN
	  start = example_get_time();
#endif    
    if(es[1]){
      //tmpVec=(sfixn *)my_calloc(dims[1], sizeof(sfixn));
      for(j=0; j<n; j+=ls[1]){
		
		copyVec_0_to_d(ls[1]-1, tmpVec, coeffs1+j);   
		
		for(k=ls[1]; k<dims[1]; k++) tmpVec[k]=0;  
		
		EX_Mont_INVTDFT_OPT2_AS_GENE_1( ls[1], tmprootsPtr, tmpVec, pPtr);
		
		copyVec_0_to_d(ls[1]-1, coeffs1+j, tmpVec);
		
		for(k=ls[1]; k<dims[1]; k++) tmpVec[k]=0;  
      }
      //my_free(tmpVec);
    }
#ifdef TIMEFN  
	  end = example_get_time();
	  invtft_time += (end - start) / 1.f;
#endif
    
    for(i=N;i>=2;i--){
      tmprootsPtr-=dims[i];

#ifdef TIMEFN
	  start = example_get_time();
#endif
	  SERBPAS::multi_mat_transpose (N, n, i, ls, coeffs1);

#ifdef TIMEFN  
	  end = example_get_time();
	  btran_time += (end - start) / 1.f;
#endif
      
#ifdef TIMEFN
	  start = example_get_time();
#endif
      if(es[i]){
	//tmpVec=(sfixn *)my_calloc(maxdim, sizeof(sfixn));
		cleanVec(maxdim-1, tmpVec);
		
		for(j=0; j<n; j+=ls[i]){
		  
		  copyVec_0_to_d(ls[i]-1, tmpVec, coeffs1+j); 
		  
		  EX_Mont_INVTDFT_OPT2_AS_GENE_1(ls[i], tmprootsPtr, tmpVec, pPtr);
		  
		  copyVec_0_to_d(ls[i]-1, coeffs1+j, tmpVec);  
		  
		  for(k=ls[i]; k<dims[i]; k++) tmpVec[k]=0;
		}
		//my_free(tmpVec);
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
    
    my_free(rootsPtr);
    my_free(tmpVec);

#ifdef TIMEFN  
	std::cout << f1tran_time << " "<<f1tft_time << " "<<f2tran_time << " "<<f2tft_time << " "<<btran_time<<" "<<invtft_time << " "<<pairwisemul_time<< " ";
    //std::cout << "Mat Tran took " << tran_time << " ms." << std::endl;
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
    SERBPAS::InitKroFFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);

    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
#ifdef TIMEFN  
	  end = example_get_time();
	  initkn_time += (end - start) / 1.f;
#endif      
    //----------------------------------------------
    SERBPAS::fftMultiD_test(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), pPtr);
      
	
	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );

#ifdef TIMEFN
	  start = example_get_time();
#endif					
    SERBPAS::fromtofftRepMultiD(N, CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    //void freeKroFFTRep(KroFFTRep * x);
    SERBPAS::freeKroFFTRep(kPtr);
#ifdef TIMEFN  
	  end = example_get_time();
	  reckn_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN  
	std::cout << initkn_time << " "<< reckn_time << " ";
#endif
  }

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
   * 
   * Return value: 
   **/
  void fftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr){
    register int i;
    sfixn m=0, n=1;
    sfixn * rootsPtr, *tmprootsPtr;
	//std::cout << "entering fftMultiD_test_1" << std::endl;

	//check Frourier prime
	for(i=1; i<=N; i++) {
	  SERBPAS::checkFrourierPrime(es[i], pPtr);
	}

    //compute the primitive roots for all dimensions
    for(i=1; i<=N; i++) { n *= dims[i]; m += dims[i]; } 
    rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
    tmprootsPtr = rootsPtr;
    for(i=1; i<=N; i++) {//parallel???
      SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr, pPtr);
      tmprootsPtr += dims[i];
    }
 
    
    //make a copy of es and dims for poly2's FFT use in parallel 
    sfixn *es_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    sfixn *dims_cp = (sfixn *) my_calloc(N+1, sizeof(sfixn));
    for(i=1; i<=N; i++) { es_cp[i] = es[i]; dims_cp[i] = dims[i]; }
   
    tmprootsPtr = SERBPAS::MultiEvalByFFT(coeffs1, N, n, es, dims, pPtr, rootsPtr);
    SERBPAS::MultiEvalByFFT(coeffs2, N, n, es_cp, dims_cp, pPtr, rootsPtr);  

    my_free(es_cp);
    my_free(dims_cp);

#ifdef TIMEFN
	float pairwisemul_time = 0.0;
	long start = example_get_time();
#endif    
    // Pairwise-Mul, result is in coeffs1
    SERBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs2, pPtr->P);
#ifdef TIMEFN  
	  long end = example_get_time();
	  pairwisemul_time += (end - start) / 1.f;
#endif   
    //is tmprootsPtr1 pointing to the last dimension??? 
    SERBPAS::InterpolByFFT(coeffs1, N, n, es, dims, pPtr, tmprootsPtr);

	//std::cout <<"before exit fftMultiD_test_1" << std::endl;
    my_free(rootsPtr);

#ifdef TIMEFN
	std::cout << pairwisemul_time << " ";
#endif
	//std::cout <<"after free before exit fftMultiD_test_1" << std::endl;
  }

  /**-------------------------------------------------------
   *
   **/
  sfixn * MultiEvalByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * rootsPtr){
    register int i;
    sfixn j, tmp;
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
      for(j=0; j<n; j += dims[1]){
		EX_Mont_DFT_OPT2_AS_GENE_1 ( dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);      
      }  
    }
#ifdef TIMEFN  
	  end = example_get_time();
	  dft_time += (end - start) / 1.f;
#endif
	
    for(i=2;i<=N;i++){
      tmprootsPtr+=dims[1];

#ifdef TIMEFN
	  start = example_get_time();
#endif	

      SERBPAS::multi_mat_transpose (N, n, i, dims, coeffs1);

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN	  
	  start = example_get_time();
#endif
      if(es[i]){
		for(j=0; j<n; j += dims[i]){
		  EX_Mont_DFT_OPT2_AS_GENE_1 (dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
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
	std::cout << tran_time << " "<<dft_time << " ";
    //std::cout << "MultiEvalByFFT Mat Tran took " << tran_time << " ms." << std::endl;
#endif	

	return tmprootsPtr;
  }

  /**---------------------------------------------------
   *InterpolByFFT: inverse DFT
   **/
  void InterpolByFFT(sfixn * coeffs1, sfixn N, sfixn n, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE *pPtr, sfixn * tmprootsPtr){
    register int i;
    sfixn j, tmp;

#ifdef TIMEFN
	float tran_time = 0.0;
	float invdft_time = 0.0;
	long start, end;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
    if(es[1]){
      for(j=0; j<n; j+=dims[1]){
		EX_Mont_INVDFT_OPT2_AS_GENE_1( dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);
      }
    }
#ifdef TIMEFN  
	  end = example_get_time();
	  invdft_time += (end - start) / 1.f;
#endif

    for(i=N; i>=2; i--){
      tmprootsPtr -= dims[i];

#ifdef TIMEFN
	  start = example_get_time();
#endif	
      
      SERBPAS::multi_mat_transpose (N, n, i, dims, coeffs1);

#ifdef TIMEFN  
	  end = example_get_time();
	  tran_time += (end - start) / 1.f;
#endif

#ifdef TIMEFN
	  start = example_get_time();
#endif
      if(es[i]){
		for(j=0; j<n; j+=dims[i]){
		  EX_Mont_INVDFT_OPT2_AS_GENE_1(dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
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
  void MultiplyByKroneckerFFT(sfixn N, preFFTRep *rep12, preFFTRep * rep1, preFFTRep * rep2,  MONTP_OPT2_AS_GENE *pPtr){
    sfixn kdg1=0, kdg2=0;
    
    KroFFTRep * kPtr = (KroFFTRep *)my_calloc(1, sizeof(KroFFTRep));
    SERBPAS::InitKroFFTRep(kPtr, BUSZS(rep12), N, 2, pPtr);
    
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));

	//-----------
	SERBPAS::checkFrourierPrime(KE(kPtr), pPtr);

    kdg1= (BUSZSI(rep1, N)+1)*CUMI(kPtr, N)-1;
    kdg2= (BUSZSI(rep2, N)+1)*CUMI(kPtr, N)-1;

    KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));

    //------------------------------------------------------------------------
    SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr),pPtr);

    //------------------------------------------------------------------------
    SERBPAS::EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), kdg2, DATSI(kPtr, 1), pPtr);

	//result is in DATSI(kPtr, 0)
	//my_free(DATSI(kPtr, 1));
			
	//delayed allocation
	//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );
				

    SERBPAS::fromtofftRepMultiD(N, CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
    //void freeKroFFTRep(KroFFTRep * x);
    SERBPAS::freeKroFFTRep(kPtr);
    
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
    register int i;
    sfixn l_AB, e1, n1;

    l_AB = degRes + 1; // compute l_AB for FFT(A,B).
    e1 = logceiling(l_AB);
	checkFrourierPrime(e1, pPtr);

    n1 = 1L << e1;
    
    // force this FFT use the n from parameter if n<n1.
    if ( (n<n1) || (degRes==0) ) { n1=n; e1=e;}
    for(i=degA+1; i<n1; i++) APtr[i]=0;
    for(i=degB+1; i<n1; i++) BPtr[i]=0;  
    
    EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, APtr, pPtr);
    EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, BPtr, pPtr);
    
    
    SERBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n1, APtr, BPtr, pPtr);
    
    EX_Mont_INVDFT_OPT2_AS_GENE_R_1(n1, e1, rootsPtr, APtr, pPtr);
    
  }


  //================================================================
  // Classical Multiplication.
  // dgs1[0], dgs2[0] are useless.
  // ccumx[i] keeps the base before ith-dimension.
  //================================================================
  void plainMultiDMul(sfixn N, sfixn * ccum, sfixn * res, sfixn * ccum1, sfixn * dgs1, sfixn * ccum2, sfixn * dgs2, sfixn * coeffs1, sfixn * coeffs2, MONTP_OPT2_AS_GENE * pPtr){
	int i, d;
	
	sfixn p=pPtr->P, R=(1L<<pPtr->Rpow)%p, SFT=pPtr->Base_Rpow;

	d=shrinkDeg(dgs1[N], coeffs1, ccum1[N]);

	for(i=0; i<=d; i++){
	  SERBPAS::decomposePoly(N, ccum, res+ccum[N]*i,  N-1, dgs1, coeffs1+ccum1[N]*i, ccum1, N, dgs2, coeffs2, ccum2, pPtr, R, SFT);
	}
  } 

  //-------------------------------
  void
  decomposePoly(sfixn N, sfixn * ccum, sfixn * res, sfixn N1, sfixn * dgs1, sfixn * coeffs1, sfixn * ccum1, sfixn N2, sfixn * dgs2, sfixn * coeffs2, sfixn * ccum2, MONTP_OPT2_AS_GENE * pPtr, sfixn R, sfixn SFT){
	int i,d;
	if(N1==0){
	  if(! coeffs1[0]) return;
	  R=( MulMod(coeffs1[0], R, pPtr->P))<<SFT;    
	  SERBPAS::decomposePoly2(N, ccum, res, R, N2, dgs2, coeffs2, ccum2, pPtr); 
	  return;
	}

	d=shrinkDeg(dgs1[N1], coeffs1, ccum1[N1]);
	for(i=0; i<=d; i++){
	  SERBPAS::decomposePoly(N, ccum, res+ccum[N1]*i, N1-1, dgs1, coeffs1+ccum1[N1]*i, ccum1, N2, dgs2, coeffs2, ccum2, pPtr, R, SFT);
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
	  for(int i=0; i<=d; i++)
        SERBPAS::decomposePoly2(N, ccum, res+ccum[1]*i, num, 0, dgs2, coeffs2+ccum2[1]*i, ccum2, pPtr);
	  return;
	}
	
	d=shrinkDeg(dgs2[N2], coeffs2, ccum2[N2]);
	for(int i=0; i<=d; i++){
	  SERBPAS::decomposePoly2(N, ccum, res+ccum[N2]*i, num, N2-1, dgs2, coeffs2+ccum2[N2]*i, ccum2, pPtr);
	}	
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
   * Computer the remainder of A is divied by B.
   *
   * Return value: 
   **/
  void
  UniFastMod_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, sfixn Lbound, sfixn * BRevInvPtr,  MONTP_OPT2_AS_GENE * pPtr){
    
	sfixn * FPtr,  *QPtr, * GPtr;
	sfixn degF, power2,power3, n2,n3,l1,l2,l3, tmp, sz;
	sfixn dg, da;
	sfixn degQ=degA-degB;
	
	
	degF=degQ;
	
	if(degF<0) {
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
	
	
	for(int i=0; i<=dg; i++) GPtr[i] = BRevInvPtr[i];
	
	FPtr = SERBPAS::reverseUni(degB, FPtr, BPtr);
	
	FPtr = SERBPAS::reverseUni(degA, FPtr, APtr);
	
	SERBPAS::EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
	
	QPtr=(sfixn * ) my_calloc(degQ+1, sizeof(sfixn));
	
	QPtr = SERBPAS::reverseUni(degQ, QPtr, FPtr);
	
	SERBPAS::cleanVec(n3-1, FPtr);
	
	
	SERBPAS::EX_Mont_FFTMul_OPT2_AS_GENE(n3, power3, degQ+degB, FPtr, degQ, QPtr, degB, BPtr, pPtr);
	
	
	for(int j=0; j<l3; j++) APtr[j] = SubMod(APtr[j], FPtr[j], pPtr->P);
	for(int j=l3; j<=degA; j++) APtr[j] = 0;
	
	my_free(FPtr);
	my_free(GPtr);
	my_free(QPtr);
	
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
	if((degA-degB)<0) return;
	if(degB>10) {
      for(i=degA-degB; i>=0; i--){
        tmp=MontMulMod_OPT2_AS_GENE(APtr[degB+i],pPtr->R2BRsft , pPtr)<<(pPtr->Base_Rpow);
        for(j=0; j<degB; j++) 
		  APtr[i+j]=SubMod(APtr[i+j], MontMulMod_OPT2_AS_GENE(BPtr[j],tmp,pPtr), p);

        APtr[i+degB]=0;
	  }
	}
	else {
	  for(i=degA-degB; i>=0; i--){
		for(j=0; j<degB; j++) 
		  APtr[i+j]=SubMod(APtr[i+j], MulMod(BPtr[j],APtr[degB+i],p), p);

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
	
	l_AB=degRes+1; // compute l_AB for TFT(A,B).
	e1=logceiling(l_AB);
	n1=1L<<e1;
	// force this FFT use the n from parameter if n<n1.
	if ( (n<n1) || (degRes==0) ) {n1=n; e1=e;}
	SERBPAS::checkFrourierPrime(e1, pPtr);
	
	rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	
	SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(e1, n1, rootsPtr, pPtr);
	
	sA = degA+1;
	sB = degB+1;
	for(int i=sA; i<n1; i++) APtr[i]=0;
	for(int i=sB; i<n1; i++) BPtr[i]=0; 
 
	EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, APtr, pPtr);
	EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, BPtr, pPtr);
	//cilk_sync;

	SERBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n1, APtr, BPtr, pPtr);
	
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
	l_AB=degRes+1; // compute l_AB for TFT(A,B).
	e1=logceiling(l_AB);
	n1=1L<<e1;
	// force this FFT use the n from parameter if n<n1.
	if ((n<n1) || (degRes==0)) {n1=n; e1=e;}

	SERBPAS::checkFrourierPrime(e1, pPtr);
	
	rootsPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	dftAPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));
	dftBPtr=(sfixn *)my_calloc(n1, sizeof(sfixn));	
	
	if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){	  
	  SERBPAS::Mont_GetNthRoots_OPT2_AS_GENE(e1,n1,rootsPtr,pPtr);
	  
	  MODPN::Mont_dft_OPT2_AS_GENE( n1, e1, rootsPtr, dftAPtr, degA, APtr, pPtr);
	  MODPN::Mont_dft_OPT2_AS_GENE( n1, e1, rootsPtr, dftBPtr, degB, BPtr, pPtr);
	  //cilk_sync;

	  for(int i=0; i<n1; i++) dftAPtr[i]=MontMulMod_OPT2_AS_GENE(dftAPtr[i], dftBPtr[i]<<BRsft, pPtr);
	  
	  MODPN::Mont_invdft_OPT2_AS_GENE_R(n1, e1, rootsPtr, dftAPtr, degRes, resPtr, pPtr);
	  
	}
	else{
	  SERBPAS::Mont_GetNthRoots_OPT2_AS_GENE_SPE(e1,n1,rootsPtr,pPtr);

	  MODPN::Mont_dft_OPT2_AS_GENE_SPE( n1, e1, rootsPtr, dftAPtr, degA, APtr, pPtr);
	  MODPN::Mont_dft_OPT2_AS_GENE_SPE( n1, e1, rootsPtr, dftBPtr, degB, BPtr, pPtr);
	  //cilk_sync;

	  for(int i=0; i<n1; i++) dftAPtr[i]=MontMulMod_OPT2_AS_GENE_SPE(dftAPtr[i], dftBPtr[i]<<BRsft, pPtr);

	  MODPN::Mont_invdft_OPT2_AS_GENE_SPE_R(n1, e1, rootsPtr, dftAPtr, degRes, resPtr, pPtr);
	  
	}
	
	
	my_free(dftAPtr);
	my_free(dftBPtr);
	my_free(rootsPtr);
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
		SERBPAS::NewtonRevInverse(i, &tRIPtrtmp, ELEMI(tris, i), ELEMI(ts, i), ts, tris, pPtr );
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
	  SERBPAS::reverseUni(BUSZSI(tPtr, 1), tmpPtr, DAT(tPtr));
	  degF=BUSZSI(tPtr, 1);

	  //=========================================================
	  SERBPAS::modularInvPM(BUSZSI(tRIPtr, 1), DAT(tRIPtr), degF, tmpPtr, e, n, pPtr);
	  //=========================================================

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

	SERBPAS::reverseMulti(BUSZSI(tPtr, N), CUMI(tPtr, N), tmpPtr, DAT(tPtr));
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
	
	for (i=1; i<=e; i++){
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
		SERBPAS::PolyCleanData(resPoly);
	  
	  if(d3>(2*cutDeg)) 
		BUSZSI(resPoly, N)=2*cutDeg; 
	  else 
		BUSZSI(resPoly, N)=d3;
	  
	  SERBPAS::decreaseKroFFTRep(krofftrep, BUSZS(resPoly));
	  
	  //=========================================================
	  SERBPAS::squarePoly_FFT(N, krofftrep, resPoly, tRIPtr, pPtr);
	  // =============================================================
	  
	  SERBPAS::PolyCleanData(tRIPtrtmp);
	  
	  if(d4> cutDeg) 
		BUSZSI(tRIPtrtmp, N)=cutDeg; 
	  else 
		BUSZSI(tRIPtrtmp, N)=d4;

	  if(d3>cutDeg) 
		BUSZSI(resPoly, N)=cutDeg; 
	  else 
		BUSZSI(resPoly, N)=d3;

	  //=========================================================
	  SERBPAS::reduceCoeffs(N, tRIPtrtmp, resPoly, ts, tris, pPtr);
	  //=========================================================

	  if(i>1)  
		SERBPAS::KroneckerCleanData(krofftrep2);

	  if(i>1) 
		SERBPAS::PolyCleanData(resPoly2);
	  
	  if(d5>2*cutDeg) 
		BUSZSI(resPoly2, N)=2*cutDeg; 
	  else 
		BUSZSI(resPoly2, N)=d5;
	  
	  SERBPAS::decreaseKroFFTRep(krofftrep2, BUSZS(resPoly2));
	  
	  //====================================================================
	  SERBPAS::mulPoly_FFT(N, krofftrep2, resPoly2, tRIPtrtmp, tPtr, pPtr);
	  //==================================================================== 
	  
	  if(d5>cutDeg) 
		BUSZSI(resPoly2, N)=cutDeg; 
	  else 
		BUSZSI(resPoly2, N)=d5;
	  
	  BUSZSI(tRIPtrtmp, N)=d4;
	  SERBPAS::PolyCleanData(tRIPtrtmp);
	  BUSZSI(tRIPtrtmp, N)=cutDeg;
	  
	  SERBPAS::reduceCoeffs(N, tRIPtrtmp, resPoly2, ts, tris, pPtr);
	  
	  //============================================
	  SERBPAS::addEqDgPoly_1(N, tRIPtr, tRIPtr, pPtr->P); //inline
	  //============================================

	  //=================================================
	  SERBPAS::subEqDgPoly_1(N, tRIPtr, tRIPtrtmp, pPtr->P, 1);//inline
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

	  SERBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(nn, tmpFPtr, tmpGPtr, pPtr);
	  
	  
	  EX_Mont_INVDFT_OPT2_AS_GENE_R_1 (nn, i, rootsPtr, tmpFPtr, pPtr);
	  
	  for(int j=0; j<halfnn; j++) 
		tmpFPtr[j]=tmpFPtr[j+halfnn];
	  
	  for(int j=halfnn; j<nn; j++) 
		tmpFPtr[j]=0;
	  
	  
	  EX_Mont_DFT_OPT2_AS_GENE_1 ( nn, i, rootsPtr, tmpFPtr, pPtr );
	  
	  SERBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(nn, tmpFPtr, tmpGPtr, pPtr);
	  
	  
	  EX_Mont_INVDFT_OPT2_AS_GENE_R_1 (nn, i, rootsPtr, tmpFPtr, pPtr);
	  
	  for(int j=halfnn; j<nn; j++) 
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
	
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(p1Ptr), BUSZS(p1Ptr), DAT(p1Ptr));
	
    if(switcher){
      //
	  //
	  SERBPAS::fftMultiD_square_test(DATSI(kPtr, 0), N, ES(kPtr), DIMS(kPtr), pPtr); }
    else{ //use KN 1D FFT
      //
	  kdg1= (BUSZSI(p1Ptr, N)+1)*CUMI(kPtr, N)-1;
	  KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));
	  //
	  SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr),pPtr);
	  //
	  SERBPAS::EX_KN_Mont_FFTSQUARE_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), pPtr);

    }
	
    SERBPAS::fromtofftRepMultiD(N,  CUM(rPtr), DAT(rPtr), CUM(kPtr), BUSZS(rPtr), DATSI(kPtr, 0));
  }
  
  //================================================================
// Squaring by Multi-dimensional FFT.
//
// es[0] and dims[0] are useless;
// es[i] is the exponent of dims[i];
// dims[i] is the # on dims-i
//================================================================

void fftMultiD_square_test(sfixn * coeffs1, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr)
{ sfixn i,j, m=0, n=1, tmp;
  sfixn * rootsPtr, * tmprootsPtr;

  // printf("\ndims in fftMultiD_test");
  // printVec(N,dims);
  // printf("\n");
  for(i=1; i<=N; i++) {n*=dims[i]; m+=dims[i];} 
  rootsPtr=(sfixn *) my_calloc(m, sizeof(sfixn));
  tmprootsPtr=rootsPtr;

  // Multi-DFT
  if(es[1]){
    //getNthRoots(es[1], dims[1], tmprootsPtr, pPtr);
    SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[1], dims[1], tmprootsPtr, pPtr);
    for(j=0; j<n; j+=dims[1]){
      //dft_1 ( dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);
      EX_Mont_DFT_OPT2_AS_GENE_1 (dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);
    }  
  }


  for(i=2;i<=N;i++){
    // last dim -> [1]
    
    tmprootsPtr+=dims[1];
    SERBPAS::multi_mat_transpose (N, n, i, dims, coeffs1);
    // multi_mat_transpose (N, i, dims, coeffs2);

    if(es[i]){
      //getNthRoots(es[i], dims[i], tmprootsPtr, pPtr);
      SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(es[i], dims[i], tmprootsPtr, pPtr);
      for(j=0; j<n; j+=dims[i]){
        //dft_1 ( dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
        EX_Mont_DFT_OPT2_AS_GENE_1 ( dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
      }}

    //  multi_mat_transpose didn't change dims accordings, so we need to do this.
    tmp=dims[1];
    dims[1]=dims[i];
    dims[i]=tmp;
    tmp=es[1];
    es[1]=es[i];
    es[i]=tmp;
  }

  // Pairwise-Mul
 
  //for(i=0; i<n; i++) coeffs1[i]=MulMod(coeffs1[i], coeffs1[i],pPtr->P);
  SERBPAS::EX_Mont_PairwiseMul_OPT2_AS(n, coeffs1, coeffs1, pPtr->P);
  // Multi-invdft

  if(es[1]){
    for(j=0; j<n; j+=dims[1]){
      //invdft ( dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);
      EX_Mont_INVDFT_OPT2_AS_GENE_1(dims[1], es[1], tmprootsPtr, coeffs1+j, pPtr);
    }}

  for(i=N;i>=2;i--){

    tmprootsPtr-=dims[i];
   
    SERBPAS::multi_mat_transpose (N, n, i, dims, coeffs1);

    if(es[i]){
      for(j=0; j<n; j+=dims[i]){
        //invdft ( dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
		EX_Mont_INVDFT_OPT2_AS_GENE_1(dims[i], es[i], tmprootsPtr, coeffs1+j, pPtr);
      }}
 
  // multi_mat_transpose didn't change dims accordings, so we need to do this.
    tmp=dims[1];
    dims[1]=dims[i];
    dims[i]=tmp;
    tmp=es[1];
    es[1]=es[i];
    es[i]=tmp;
  }

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
    for(int i=sA; i<n1; i++) APtr[i]=0;
	//4 96 96 96 96 1,
	//degA=33554431
	//Segmentation fault
    
    EX_Mont_DFT_OPT2_AS_GENE_1( n1, e1, rootsPtr, APtr, pPtr);
    
    SERBPAS::EX_Mont_PairwiseMul_OPT2_AS_R(n1, APtr, APtr, pPtr);
    
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
  { register int i;
	
	backupData(fromPtr); //inline
	backupData(toPtr);

	decreaseOneDim(fromPtr); //inline
	decreaseOneDim(toPtr);

	for(i=0; i<=BUSZSI(toPtr, N); i++){
	  //===============================================
	  SERBPAS::MultiMod(N-1, toPtr, fromPtr, ts, tris, pPtr);
	  //===============================================
	  nextCoefData(fromPtr, N); //inline DAT(fromPtr)+=CUMI(fromPtr, N);
	  nextCoefData(toPtr, N); 
	}

	increaseOneDim(fromPtr); //inline
	increaseOneDim(toPtr);

	restoreData(fromPtr); //inline
	restoreData(toPtr);
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
	  SERBPAS::MultiMod_DF(N, outPtr, inPtr, ts, tris, pPtr);	  
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
	
	if (zeroCoefp(DAT(inPtr), SIZ(inPtr)) ) return;
	
	SERBPAS::MultiMod_1(N, inPtr, ts, tris, pPtr);
	SERBPAS::fromtofftRep(N, CUM(outPtr), DAT(outPtr), CUM(inPtr),  BUSZS(outPtr), DAT(inPtr));
	
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
   * MultiMod_1:
   * @N: number of variables.
   * @inPtr: a C-Cube polynomial.
   * @ts: a zero-dim trinagular set.
   * @tris: inverses of 'ts'.
   * @pPtr: the information of prime number p.
   * 
   * Reduce 'inPtr' by 'ts' in-place.
   * I.e. 'inPtr'= 'inPtr' modulo 'ts'.
   * Return value: 
   **/
  void MultiMod_1(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	register int i;
	int deg=0, d;
	
	sfixn * restoreDat;
	sfixn restoreSiz;
	preFFTRep  tmpIn;
	
	if(N==1){ 
	  deg=shrinkDeg(BUSZSI(inPtr, N), DAT(inPtr), 1); //inline
	  
	  if (!deg) return;
	  if ( ( (BDSI(ts, 1)) > DivCut1 ) ) {
		SERBPAS::UniFastMod_1(deg, DAT(inPtr), BUSZSI((ELEMI(ts, 1)), 1), DAT(ELEMI(ts, 1)), 
							(NLB(tris))[1], DAT(ELEMI(tris, 1)) ,pPtr);  
	  }
	  else{
		SERBPAS::UniPlainMod_1(deg, DAT(inPtr), BUSZSI((ELEMI(ts, 1)), 1), 
							 DAT(ELEMI(ts, 1)), pPtr); 
	  }
	  return;
	}
	
	deg=shrinkDeg(BUSZSI(inPtr, N), DAT(inPtr), CUMI(inPtr, N));
	restoreDat=DAT(inPtr);
	N(inPtr)=N-1;
	restoreSiz=SIZ(inPtr);  
	SIZ(inPtr)=CUMI(inPtr, N);
	
	for(i=0; i<=deg; i++){
	  SERBPAS::MultiMod_1(N-1, inPtr, ts, tris, pPtr);
	  DAT(inPtr)+=CUMI(inPtr, N); 
	}
	
	N(inPtr)=N;
	DAT(inPtr)=restoreDat;
	SIZ(inPtr)=restoreSiz;
	d=BDSI(ts, N);
	
	SERBPAS::copyVec_1_to_n(N-1, CUTS(inPtr), BDS(ts));
	deg=shrinkDeg(CUTSI(inPtr, N), DAT(inPtr), CUMI(inPtr, N));
	SERBPAS::copyVec_1_to_n(N-1, CUTS(inPtr), BUSZS(inPtr));
	
	BDSI(ts, N)=deg;
	SERBPAS::InitOneReducedPoly(&tmpIn, N, BDS(ts));  
	SERBPAS::fromtofftRep(N, tmpIn.accum, tmpIn.data, CUM(inPtr),  BUSZSD(tmpIn), DAT(inPtr));
	BDSI(ts, N)=d;
	
	if( (((N==2) && (d>DivCut2)) || ((N==3) && (d>DivCut3)) || ((N>3) && (d>DivCutN))) ){
	  //printf("using MultiUniFastMod_1\n");
	  //fflush(stdout);
	  
	  SERBPAS::MultiUniFastMod_1(N, &tmpIn, deg, inPtr, ts, tris, pPtr);
	  
    }
	else{
	  //printf("using MultiUniPlainMod_1\n");
	  //fflush(stdout);
	  
	  SERBPAS::MultiUniPlainMod_1(N, deg, &tmpIn, BDSI(ts, N)+1,  ELEMI(ts, N), ts, tris, pPtr);
	  SERBPAS::fromtofftRep(N, CUM(inPtr), DAT(inPtr), tmpIn.accum,  BUSZSD(tmpIn), tmpIn.data);
	  
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
	
	resPoly=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	resPoly2=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	out=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	
	krofftrep=(KroFFTRep *)my_calloc(1,sizeof(KroFFTRep));
	krofftrep2=(KroFFTRep *)my_calloc(1,sizeof(KroFFTRep));
	e=logceiling(n-m+1);
	
	if(! (EXSTI(tris, N))) {  
	  SERBPAS::initCopyOneRevInv(out, ELEMI(tris, N));
	  SERBPAS::NewtonRevInverse(N, out, ELEMI(tris, N), ELEMI(ts, N), ts, tris, pPtr );
	  freePoly(out);
	}
	cutDg=n-m;
	d1=BUSZSI(tmpIn, N);
	d2=BUSZSI(ELEMI(tris, N), N);
	BUSZSI(tmpIn, N)=cutDg;
	BUSZSI(ELEMI(tris, N), N)=cutDg;
	
	SERBPAS::reverseMulti_1(n, CUMI(tmpIn, N), DAT(tmpIn));
	SERBPAS::InitResPoly(resPoly, N, BUSZS(tmpIn), BUSZS(ELEMI(tris, N)));
	
	//=====================================================================
	SERBPAS::EX_mulPoly(N, resPoly, tmpIn, ELEMI(tris, N), pPtr);
	//=====================================================================
	
	SERBPAS::reverseMulti_1(cutDg, CUMI(resPoly, N), DAT(resPoly));
	
	BUSZSI(tmpIn, N)=d1;
	BUSZSI(ELEMI(tris, N), N)=d2;
	
	SERBPAS::PolyCleanData(tmpIn);
	
	BUSZSI(tmpIn, N)=cutDg;
	
	SERBPAS::reduceCoeffs(N, tmpIn, resPoly, ts, tris, pPtr);
	
	freePoly(resPoly); 
	
	SERBPAS::InitResPoly(resPoly2, N,  BUSZS(ELEMI(ts, N)), BUSZS(tmpIn));
	
	//=====================================================================
	SERBPAS::EX_mulPoly(N, resPoly2, ELEMI(ts, N), tmpIn,  pPtr);
	//=====================================================================
	
	BUSZSI(tmpIn, N)=d1;
	
	SERBPAS::PolyCleanData(tmpIn);
	
	SERBPAS::reduceCoeffs(N, tmpIn, resPoly2, ts, tris, pPtr);
	
	freePoly(resPoly2);
	
	SERBPAS::subPoly_1(N, inPtr, tmpIn, pPtr->P);
	
	my_free(resPoly); 
	my_free(resPoly2);
	my_free(out);
	my_free(krofftrep); 
	my_free(krofftrep2);
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
	  
	  SERBPAS::EX_mulPoly_TFTFFT(N, resPtr, f1, f2, pPtr);
	  
	}else{
	  SERBPAS::plainMultiDMul(N, CUM(resPtr), DAT(resPtr), CUM(f1), BUSZS(f1), CUM(f2), BUSZS(f2), DAT(f1), DAT(f2), pPtr);
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
	//double t;
	
	if( SERBPAS::forcutoff_Multi_TFT_FFT(N, BUSZS(f1), BUSZS(f2), 8) ==0 ){
	  SERBPAS::EX_mulPoly_FFT(N, resPtr, f1, f2, pPtr);
	   
	}
	else{
	  SERBPAS::InitKroTFTRep(&kro, BUSZS(resPtr), N, 2, pPtr); 
	  //================================================================
	  SERBPAS::polyMul_TFT(N, &kro, resPtr, f1, f2,  pPtr);
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
    SERBPAS::mulPoly_FFT(N, &kro, resPtr, f1, f2,  pPtr);
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
	
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(p1Ptr), BUSZS(p1Ptr), DAT(p1Ptr));
    SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(p2Ptr), BUSZS(p2Ptr), DAT(p2Ptr));
	
    if((switcher) || (KN(kPtr)<4)){
      //
	  //================================================================
	  //SERBPAS::fftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), pPtr);
	  SERBPAS::fftMultiD_test(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), pPtr);
	  //================================================================
	}
    else{
      //KN FFT
	  kdg1= (BUSZSI(p1Ptr, N)+1)*CUMI(kPtr, N)-1;
	  kdg2= (BUSZSI(p2Ptr, N)+1)*CUMI(kPtr, N)-1;
	  KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));
	  
	  SERBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr),pPtr);
	  
	  //================================================================
	  SERBPAS::EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), kdg2, DATSI(kPtr, 1), pPtr);
	  //================================================================
	  
    }

    SERBPAS::fromtofftRepMultiD(N,  CUM(rPtr), DAT(rPtr), CUM(kPtr), BUSZS(rPtr), DATSI(kPtr, 0));

  }
  
  
  //-----------------------------------------------
  void polyMul_TFT(sfixn N, KroTFTRep * kPtr, preFFTRep * rPtr, preFFTRep * p1Ptr, preFFTRep * p2Ptr,  MONTP_OPT2_AS_GENE * pPtr){
	
	SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(p1Ptr), BUSZS(p1Ptr), DAT(p1Ptr));
	
	SERBPAS::fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(p2Ptr), BUSZS(p2Ptr), DAT(p2Ptr));
	
	//================================================================
	//SERBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	SERBPAS::tftMultiD_test(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), pPtr);
	//================================================================
	
	SERBPAS::fromtofftRepMultiD(N,  CUM(rPtr), DAT(rPtr), CUM(kPtr), BUSZS(rPtr), DATSI(kPtr, 0));
  }
  
  //=========================================================================
  // Bottom-up level-by-level reduction.
  // using TFTs.
  //=========================================================================
  void
  MultiMod_BULL(sfixn N, preFFTRep * outPtr, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr)
  { 
	
	if (zeroCoefp(DAT(inPtr), SIZ(inPtr)) ) 
	  return;
	
	//=========================================
	SERBPAS::MultiMod_1_BULL(N, inPtr, ts, tris, pPtr);
	//=========================================
	
	SERBPAS::fromtofftRep(N, CUM(outPtr), DAT(outPtr), CUM(inPtr),  BUSZS(outPtr), DAT(inPtr));
	
  }
  
  /**
   * MultiMod_1_BULL:
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
  void MultiMod_1_BULL(sfixn N, preFFTRep * inPtr, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr){
	
	sfixn *lastNonZeroPtr, *datPtr, *tmpVec;
	sfixn siz, deg, i, j;
	preFFTRep * tmpInPtr, *cpShellInPtr;
	
	if(N<1) { std::cout<<"Number of variables < 1."<<std::endl; exit(1); }
	
	siz=BUSZSI(inPtr, 1)+1;
	
	lastNonZeroPtr=DAT(inPtr)+shrinkDeg(SIZ(inPtr)-1, DAT(inPtr), 1);  
	
	for(datPtr=DAT(inPtr); datPtr<=lastNonZeroPtr; datPtr+=siz){
	  deg=shrinkDeg(BUSZSI(inPtr, 1), datPtr, 1);
	  
	  if (((BDSI(ts, 1))>DivCut1)) {
		SERBPAS::UniFastMod_1(deg, datPtr, BUSZSI((ELEMI(ts, 1)), 1), DAT(ELEMI(ts, 1)), 
					 (NLB(tris))[1], DAT(ELEMI(tris, 1)) ,pPtr); 
	  }
	  else{
		SERBPAS::UniPlainMod_1(deg, datPtr, BUSZSI((ELEMI(ts, 1)), 1), 
					  DAT(ELEMI(ts, 1)), pPtr); 
	  }
	}
	
	
	if(N==1)  return;
	
	tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	tmpVec=(sfixn *)my_calloc(N+1, sizeof(sfixn));
	
	for(i=2; i<N-1; i++){  
	  
	  for(datPtr=DAT(inPtr); datPtr<=lastNonZeroPtr; datPtr+=CUMI(inPtr, i+1)){
		deg=shrinkDeg(CUTSI(inPtr, i), datPtr, CUMI(inPtr, i));    
  
		for(j=0; j<i;j++) 
		  tmpVec[j]=BDSI(ts, j);
		
		tmpVec[i]=deg;
        
		SERBPAS::InitOneReducedPoly(tmpInPtr, i, tmpVec);
		SERBPAS::fromtofftRep(i, tmpInPtr->accum, tmpInPtr->data, CUM(inPtr), 
					 BUSZS(tmpInPtr), datPtr);

		cpShellInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));

		SERBPAS::copyPolyPointers(cpShellInPtr, inPtr);

		N(cpShellInPtr)=i;
		SIZ(cpShellInPtr)=CUMI(inPtr, i+1);
		DAT(cpShellInPtr)=datPtr;
		
		if((((i==2)&&(BDSI(ts, i)>DivCut2)) || ((i==3)&&(BDSI(ts, i)>DivCut3)) || ((i>3)&&(BDSI(ts, i)>DivCutN)))){		  
		  SERBPAS::MultiUniFastMod_1_TFT(i, tmpInPtr, deg, cpShellInPtr, ts, tris, pPtr);
		}
		else{			
		  SERBPAS::MultiUniPlainMod_1(i, deg, tmpInPtr, BDSI(ts, i)+1,  ELEMI(ts, i), ts, tris, pPtr);
		  
		  SERBPAS::fromtofftRep(i, CUM(cpShellInPtr), DAT(cpShellInPtr), CUM(tmpInPtr),  BUSZS(tmpInPtr), DAT(tmpInPtr));
		  
		}
		
		freePoly(tmpInPtr);
		my_free(cpShellInPtr);
	  }
	  
	}
	
	my_free(tmpVec);
	my_free(tmpInPtr);    
	
	if(N>=3){
	  
	  tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	  tmpVec=(sfixn *)my_calloc(N, sizeof(sfixn));
	  for(datPtr=DAT(inPtr); datPtr<=lastNonZeroPtr; datPtr+=CUMI(inPtr, N)){
		deg=shrinkDeg(CUTSI(inPtr, N-1), datPtr, CUMI(inPtr, N-1));
      
		for(j=0; j<N-1;j++) 
		  tmpVec[j]=BDSI(ts, j);
		
		tmpVec[N-1]=deg;
		SERBPAS::InitOneReducedPoly(tmpInPtr, N-1, tmpVec);
		SERBPAS::fromtofftRep(N-1, tmpInPtr->accum, tmpInPtr->data, CUM(inPtr), 
					 BUSZS(tmpInPtr), datPtr);

		cpShellInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
		SERBPAS::copyPolyPointers(cpShellInPtr, inPtr);

		N(cpShellInPtr)=N-1;
		SIZ(cpShellInPtr)=CUMI(inPtr, N);
		DAT(cpShellInPtr)=datPtr;
		
		if((((N-1==2)&&(BDSI(ts, N-1)>DivCut2)) || ((N-1==3)&&(BDSI(ts, N-1)>DivCut3)) || ((N-1>3)&&(BDSI(ts, N-1)>DivCutN)))){
		  
		  SERBPAS::MultiUniFastMod_1_TFT(N-1, tmpInPtr, deg, cpShellInPtr, ts, tris, pPtr);
		}
		else{
		  
		  SERBPAS::MultiUniPlainMod_1(N-1, deg, tmpInPtr, BDSI(ts, N-1)+1,  ELEMI(ts, N-1), ts, tris, pPtr);
		  SERBPAS::fromtofftRep(N-1, CUM(cpShellInPtr), DAT(cpShellInPtr), CUM(tmpInPtr),  BUSZS(tmpInPtr), DAT(tmpInPtr));
		}
		
		
		freePoly(tmpInPtr);
		my_free(cpShellInPtr);
	  }
	  my_free(tmpVec);
	  my_free(tmpInPtr);
	  
	}
	
	tmpInPtr=(preFFTRep *)my_malloc(sizeof(preFFTRep));
	tmpVec=(sfixn *)my_calloc(N+1, sizeof(sfixn));
	deg=shrinkDeg(BUSZSI(inPtr, N), DAT(inPtr), CUMI(inPtr, N));

	for(j=0; j<N; j++) 
	  tmpVec[j]=BDSI(ts, j);

	tmpVec[N]=deg;
	SERBPAS::InitOneReducedPoly(tmpInPtr, N, tmpVec);
	SERBPAS::fromtofftRep(N, tmpInPtr->accum, tmpInPtr->data, CUM(inPtr), 
				 BUSZS(tmpInPtr), DAT(inPtr));
	SERBPAS::MultiUniFastMod_1_TFT(N, tmpInPtr, deg, inPtr, ts, tris, pPtr);
	
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
	d2=BUSZSI(ELEMI(tris, N), N);
	BUSZSI(tmpIn, N)=cutDg;
	BUSZSI(ELEMI(tris, N), N)=cutDg;

	SERBPAS::reverseMulti_1(n, CUMI(tmpIn, N), DAT(tmpIn));
	SERBPAS::InitResPoly(resPoly, N, BUSZS(tmpIn), BUSZS(ELEMI(tris, N)));

	SERBPAS::InitKroTFTRep(krotftrep, BUSZS(resPoly), N, 2, pPtr);

	//=====================================================================
    SERBPAS::polyMul_TFT(N, krotftrep, resPoly, tmpIn, ELEMI(tris, N), pPtr);
	//============================================================

	SERBPAS::freeKroTFTRep(krotftrep);
	
	SERBPAS::reverseMulti_1(cutDg, CUMI(resPoly, N), DAT(resPoly));
	
	BUSZSI(tmpIn, N)=d1;
	BUSZSI(ELEMI(tris, N), N)=d2;
	
	SERBPAS::PolyCleanData(tmpIn);
	
	BUSZSI(tmpIn, N)=cutDg;
	
	SERBPAS::reduceCoeffs(N, tmpIn, resPoly, ts, tris, pPtr);
	
	freePoly(resPoly); 
	
	SERBPAS::InitResPoly(resPoly2, N,  BUSZS(ELEMI(ts, N)), BUSZS(tmpIn));
	SERBPAS::InitKroTFTRep(krotftrep2, BUSZS(resPoly2), N, 2, pPtr); 
	
	//if(N==topN) t=gettime(); 
	//================================================================
    SERBPAS::polyMul_TFT(N, krotftrep2, resPoly2, ELEMI(ts, N), tmpIn,  pPtr);
	//================================================================
	//if(N==topN) top2MulTime1+=gettime()-t;
	
	SERBPAS::freeKroTFTRep(krotftrep2);
	BUSZSI(tmpIn, N)=d1;
	
	SERBPAS::PolyCleanData(tmpIn);
	
	// PureFast_
	
	SERBPAS::reduceCoeffs(N, tmpIn, resPoly2, ts, tris, pPtr);
	
	freePoly(resPoly2);
	
	SERBPAS::subPoly_1(N, inPtr, tmpIn, pPtr->P);
	
	my_free(resPoly); 
	my_free(resPoly2);
	my_free(out);
	my_free(krotftrep); 
	my_free(krotftrep2);
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
		SERBPAS::UniFastMod_1(d1, DAT(FPtr), d2, DAT(GPtr), 
					 (NLB(tris))[1], DAT(ELEMI(tris, 1)) ,pPtr);}
	  else{
		SERBPAS::UniPlainMod_1(d1, DAT(FPtr), d2, DAT(GPtr), pPtr);}
	}
	else{
	  SERBPAS::InitOneReducedPoly(&co, N-1, BUSZS(FPtr));
      while(d>=0){
		SERBPAS::PolyCleanData(&co);
        SERBPAS::getCoefMulti(N, &co, FPtr, d1);

        d1 = SERBPAS::fmecgEqDg_1(N, FPtr, d, &co, GPtr, ts, tris, pPtr);
		
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
	register sfixn i;
	sfixn d1,d2,d;
	preFFTRep out, res;

	SERBPAS::InitOneReducedPoly(&out, N-1, BUSZS(f1)); 
	d1=shrinkDeg(BUSZSI(f1, N), DAT(f1), CUMI(f1, N));
	d2=d1-e;
	backupData(f1);
	backupData(f2);
	decreaseOneDim(f1);  
	decreaseOneDim(f2);
	nextMCoefData(f1,N,e); 
 
	SERBPAS::InitResPoly(&res, N-1,  BUSZS(co), BUSZS(f2));
	for(i=0; i<=d2; i++){
	  SERBPAS::PolyCleanData(&out);
	  SERBPAS::PolyCleanData(&res);
	  	  
	  SERBPAS::mul_Reduced(&res, N-1, &out, co, f2, ts, tris, pPtr);
	  
	  SERBPAS::subEqDgPoly_1(N-1, f1, &out, pPtr->P, 1); //inline
	  nextCoefData(f1, N);   
	  nextCoefData(f2, N);
	}
	freePoly(&res);
	increaseOneDim(f1);
	increaseOneDim(f2);
	restoreData(f1);
	restoreData(f2);
	freePoly(&out);

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
	d3=BUSZSI(f1, N);
	d4=BUSZSI(f2, N);
	
	if((sz1>=MulCut)||(sz2>=MulCut)){
	  d5= BUSZSI(resPtr, N);
	  BUSZSI(f1, N)=d1;
	  BUSZSI(f2, N)=d2;
	  BUSZSI(resPtr, N)=d1+d2;

	  SERBPAS::InitKroFFTRep(&kro, BUSZS(resPtr), N, 2, pPtr); 
	  //================================================================
	  SERBPAS::mulPoly_FFT(N, &kro, resPtr, f1, f2,  pPtr);
	  //================================================================
	  BUSZSI(f1, N)=d3;
	  BUSZSI(f2, N)=d4;
	  BUSZSI(resPtr, N)=d5;
	  
	  SERBPAS::freeKroFFTRep(&kro);
	}else{
	  
	  if(N==1){
		SERBPAS::EX_Mont_PlainMul_OPT2_AS_GENE(BUSZSI(resPtr, 1), DAT(resPtr), d1, DAT(f1), d2, DAT(f2), pPtr); 
		
	  } else{
		SERBPAS::plainMultiDMul(N, CUM(resPtr), DAT(resPtr), CUM(f1), BUSZS(f1), CUM(f2), BUSZS(f2), DAT(f1), DAT(f2), pPtr); 
	  }
	  
	}
	
	SERBPAS::MultiMod(N, out, resPtr, ts, tris, pPtr);
	return out;
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
   * Claasical univariate multiplication 
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
	  SERBPAS::MultiMod_DF(N, outPtr, inPtr, ts, tris, pPtr);
	}else
	  {
		SERBPAS::MultiMod_BULL(N, outPtr, inPtr, ts, tris, pPtr);
		
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
	SERBPAS::initTriRevInvSet( dgs, N, tris, ts);
	SERBPAS::getRevInvTiSet(dgs, N, tris, ts, pPtr);
	my_free(dgs);
	return tris;
	
  }
  
    
}
//end of SERBPAS
