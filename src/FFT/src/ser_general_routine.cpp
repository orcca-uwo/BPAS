// -*- C++ -*- 
// ser_general_routine.cilk
//@author Yuzhen Xie

#include <cilk/cilk.h>

#include "../../../include/FFT/src/modpn.h"
#include "../../../include/FFT/src/include/example_util_gettime.h"
#include "../../../include/FFT/src/modpn_export.h"
#include "../../../include/FFT/src/ser_general_routine.h"

#include <iostream>
//#include <reducer_opadd.h>

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef int sfixn;

namespace SERBPAS {

  //---------------------------------------------------
  void checkFrourierPrime(sfixn e, MONTP_OPT2_AS_GENE * pPtr){
	if(e > (pPtr->Npow)){
	  std::cout << "The given Fourier prime number can't hanlde this FFT multiplication in terms of the FFT-size." << pPtr->P << ", "  << e << std::endl;
	  std::cout << "Please choose another Fourier prime with bigger FFT-size." << std::endl;
	  exit (1);
	}
  }

  /**--------------------------------------------------------
   * EX_Mont_GetNthRoots_OPT2_AS_GENE:
   * @e: 
   * @n: equals 2^e and divides p-1
   * @rootsPtr: (output) an array of size n which contains the powers of the primitive root
   * @pPtr: prime number structire for the prime number p
   * 
   *  Use the special code (for p-1 is a power of 2) if applicable
   *  otherwise use the generic code.
   * 
   * Return value: returns `rootsPtr`
   **/
  void EX_Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr){

    if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){
      SERBPAS::Mont_GetNthRoots_OPT2_AS_GENE(e, n, rootsPtr, pPtr);
    }
    else{
      SERBPAS::Mont_GetNthRoots_OPT2_AS_GENE_SPE(e, n, rootsPtr, pPtr); 
    }
  }

  //======================================================================
  // get the n-th root of unity. 
  //======================================================================
  
  // type: local.
  // note: good for all machine word Fourier Prime.
  /**
   * Mont_GetNthRoots_OPT2_AS_GENE:
   * @e: 
   * @n: equals 2^e and divides p-1
   * @rootsPtr: (output) an array of size n which contains the powers of the primitive root
   * @pPtr: prime number structire for the prime number p
   *  
   * Return value: in place of `rootsPtr`
   **/
  void Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr){
    //register int j;
    sfixn root, rootR, R=(1L<<pPtr->Rpow)%pPtr->P, R_2=MulMod(R, R, pPtr->P), BRsft=pPtr->Base_Rpow;
    
    // printf("The input FFT-N:%ld, this p:%ld can handle at most: %ld. !!!\n\n\n\n",n,pPtr->P, (1L<<pPtr->Npow));
    
    //root=EX_GetPrimitiveNthRoot(e, n, pPtr->P);
    root=PowerMod(pPtr->Max_Root, 1L<<((pPtr->Npow)-e), pPtr->P);
    
    rootR=MontMulMod_OPT2_AS_GENE(root, R_2<<BRsft, pPtr);
    rootsPtr[0]=R<<BRsft;
    rootsPtr[1]=rootR;
    //cilk_for(int j=2; j<n; j++) { data dependence 
	for(int j=2; j<n; j++) {
      rootsPtr[j] = MontMulMod_OPT2_AS_GENE(rootsPtr[j-1], rootR<<BRsft, pPtr);
      rootsPtr[j-1] <<= BRsft;
    }
    rootsPtr[n-1]<<=BRsft;
  }

  /**
   * Mont_GetNthRoots_OPT2_AS_GENE_SPE:
   * @e: 
   * @n: equals 2^e and divides p-1
   * @rootsPtr: (output) an array of size n which contains the powers of the primitive root
   * @pPtr: prime number structire for the prime number p
   * 
   *  Assume p is p=N-1, where N is a power of 2.
   *
   * Return value: returns `rootsPtr`
   **/
  void Mont_GetNthRoots_OPT2_AS_GENE_SPE(sfixn e, sfixn n, sfixn * rootsPtr, MONTP_OPT2_AS_GENE * pPtr){
    //register int j;
    sfixn root, rootR, R=(1L<<pPtr->Rpow)%pPtr->P, R_2=MulMod(R, R, pPtr->P), BRsft=pPtr->Base_Rpow;
 
    root=PowerMod(pPtr->Max_Root, 1L<<((pPtr->Npow)-e), pPtr->P);
    rootR=MontMulMod_OPT2_AS_GENE_SPE(root, R_2<<BRsft, pPtr);
    rootsPtr[0]=R<<BRsft;
    rootsPtr[1]=rootR;

    //cilk_for (int j=2; j<n; j++) { //potential of data race, use reducer??
	for (int j=2; j<n; j++) {
      rootsPtr[j] = MontMulMod_OPT2_AS_GENE_SPE(rootsPtr[j-1], rootR<<BRsft, pPtr);
      rootsPtr[j-1]<<=BRsft;
    }
    rootsPtr[n-1]<<=BRsft;
  }
  
  /**
   * EX_Mont_PairwiseMul_OPT2_AS:
   * @n: the FFT size
   * @APtr: coefficient vector of polynomial A
   * @BPtr: coefficient vector of polynomial B
   * @p: prime number
   * 
   * Pairwise multiplicaiton. Exported. No worries about the big R.
   * 
   * Return value: 
   **/    
  void 
  EX_Mont_PairwiseMul_OPT2_AS(sfixn n, sfixn * APtr, sfixn * BPtr, sfixn p){
    //register int i;
	for (int i=0; i<n; i++) APtr[i] = MulMod(APtr[i], BPtr[i], p);
  }

  /**
   * EX_Mont_PairwiseMul_OPT2_AS_R:
   * @n: the FFT size
   * @APtr: coefficient vector of polynomial A
   * @BPtr: coefficient vector of polynomial B
   * @pPtr: prime number structure
   * 
   * Pairwise multiplicaiton. Big R issue *NOT* handled. 
   * This function is for expert usage only.
   * 
   * Return value: 
   **/    
  void 
  EX_Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr){
     
    if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){
      SERBPAS::Mont_PairwiseMul_OPT2_AS_R(n, APtr, BPtr, pPtr);
    }
    else{
      SERBPAS::Mont_PairwiseMul_OPT2_AS_SPE_R(n, APtr, BPtr, pPtr);
    }
  }

  //======================================================================
  // pairwise mul
  //====================================================================== 
  void 
  Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr){
    //register int i;
    sfixn BRsft=pPtr->Base_Rpow;
    for(int i=0; i<n; i++) APtr[i]=MontMulMod_OPT2_AS_GENE(APtr[i], BPtr[i]<<BRsft, pPtr);
  }

  //------------------------------
  void 
  Mont_PairwiseMul_OPT2_AS_SPE_R(sfixn n, sfixn * APtr, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr){
    //register int i;
    sfixn BRsft=pPtr->Base_Rpow;
    for(int i=0; i<n; i++) APtr[i]=MontMulMod_OPT2_AS_GENE_SPE(APtr[i], BPtr[i]<<BRsft, pPtr);
  }
  
  //Matteo's rectangular matrix transpose---------------
  //out-of-place transpose A[i0..i1][j0..j1] into B[j0..j1][i0..i1]
  //then copy back to A
  // n: size of A
  //row major layout
  void 
  transpose(sfixn *A, sfixn lda, sfixn *B, sfixn ldb,
				 sfixn i0, sfixn i1, sfixn j0, sfixn j1) {
    const int THRESHOLD = 8;
    //const int THRESHOLD = 16;
    //const int THRESHOLD = 17;

  tail:
	sfixn di = i1 - i0, dj = j1 - j0;
	if (di >= dj && di > THRESHOLD) {
	  sfixn im = (i0 + i1) / 2;
	  SERBPAS::transpose(A, lda, B, ldb, i0, im, j0, j1);
	  i0 = im; goto tail;
	} else if (dj > THRESHOLD) {
	  sfixn jm = (j0 + j1) / 2;
	  SERBPAS::transpose(A, lda, B, ldb, i0, i1, j0, jm);
	  j0 = jm; goto tail;
	} else {
	  for (int i = i0; i < i1; ++i)
		for (int j = j0; j < j1; ++j) 
		  B[j * ldb + i] = A[i * lda + j];
	}
	//cilk_for(int j = 0; j < n; j++) A[j] = B[j];
  }

  /* Matteo: Traverse the trapezoidal space (i, j) where
   i0 <= i < i1
   j0 + (i - i0) * dj0 <= j < j1 
 */
  //square matrix A in place transposition
  void sqtranspose(sfixn *A, sfixn lda,
				   sfixn i0, sfixn i1,
				   sfixn j0, sfixn dj0, sfixn j1 /*, int dj1 = 0 */){
    const int THRESHOLD = 8;
    //const int THRESHOLD = 16;
    //const int THRESHOLD = 17;

  tail:
	sfixn di = i1 - i0, dj = j1 - j0;
	if (dj >= 2 * di && dj > THRESHOLD) {
	  sfixn dj2 = dj / 2;
	  SERBPAS::sqtranspose(A, lda, i0, i1, j0, dj0, j0 + dj2);
	  j0 += dj2; dj0 = 0; goto tail;
	} else if (di > THRESHOLD) {
	  sfixn di2 = di / 2;
	  SERBPAS::sqtranspose(A, lda, i0, i0 + di2, j0, dj0, j1);
	  i0 += di2; j0 += dj0 * di2; goto tail;
	} else {
	  for (int i = i0; i < i1; ++i) {
		for (int j = j0; j < j1; ++j) {
		  sfixn x = A[j * lda + i];
		  A[j * lda + i] = A[i * lda + j];
		  A[i * lda + j] = x;
		}
		j0 += dj0;
	  }
	}
  }
  
  /**------------------------------------------------------------------
   * multi_mat_transpose:
   * @N: number of dimensions of Matrix M.
   * @n: number of entries in M.
   * @dm: dm-th dimension.
   * @dims: dimensions vector. dims[i] keep the size of dimension i 
   *        (dims[0] is not used). the size of dims[] is N+1.
   * @data: entries of M, stored in a consecutive memory space of type
   *        sfixn aligned by 16 bytes using my_calloc().
   *        
   * We transpose matrix M from dimension 'dm' to dimension 1. 
   * the rest size-2 dimensions remain.
   * dims array  should be changed corresponding to the matrix
   * tranposition, BUT not in this function. Do it after MultiD-FFT/TFT.
   * 
   * Return value: 
   **/
  void 
  multi_mat_transpose(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data){
	//std::cout<<"dims[1], dims[2]: "<< dims[1] <<", "<< dims[2] <<std::endl;
	//for(int j = 0; j < n; j++) std::cout<<data[j]<<", ";
	//std::cout<<"----end of in data----" <<std::endl;

	//#ifdef MATTEOTRAN2
	if (N==2) { //transpose from dim-2 to dim-1
	  sfixn ni = dims[1];
	  sfixn nj = dims[2];
	  sfixn lda = nj;
	  sfixn ldb = ni;
	  if (ni==nj){
		SERBPAS::sqtranspose(data, lda, 0, ni, 0, 1, ni); 
	  }else{
		//std::cout<<"Matteo's rect-transpose" <<std::endl;
		sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
		SERBPAS::transpose(data, ldb, B, lda, 0, nj, 0, ni);

		for(int j = 0; j < n; j++) data[j] = B[j];
		my_free(B);
	  }
	}else{
	  //#endif
	  //std::cout<<"multi-decompose" <<std::endl;
	  register int i;
	  sfixn * accum, * tmpPtr;

	  //std::cout <<"enter multi_mat_transpose" << std::endl;
	  
	  accum = (sfixn *)my_calloc(N+2, sizeof(sfixn));
	  
	  // accum[1] keeps dim_1's interval, which is 1.
	  // accum[2] keeps dim_2's interval, which is #(dim_1);
	  // accum[3] keeps dim_3's inverval, which is #(dim_1)*#(dim_2)
	  // ...
	  accum[1] = 1;
	  for(i = 2; i <= N+1; i++) accum[i] = dims[i-1] * accum[i-1];
	  //n=accum[N]*dims[N];
	  
#ifdef TRACETRAN
	  std::cout << "mat_tran N, n, dm: "<< N <<", "<<n <<", "<<dm <<", "<<std::endl;
	  std::cout << "dims: ";
	  for (i=1; i<=N; i++) {
		std::cout << i <<" "<<dims[i] <<", ";
	  } 
	  std::cout << std::endl;
	  std::cout << "accums: ";
	  for (i=1; i<=N+1; i++) {
		std::cout << i <<" "<<accum[i] <<", ";
	  }
	  std::cout << std::endl;
#endif
	  tmpPtr = (sfixn * ) my_calloc(n, sizeof(sfixn));
#ifdef TRACETRAN
	  float decomposePar_time = 0.0;
	  long start = example_get_time();
#endif    
	  SERBPAS::decompose(N, dm, accum, dims, tmpPtr, data);
#ifdef TRACETRAN 
	  long end = example_get_time();
	  decomposePar_time += (end - start) / 1.f;
	  std::cout << "decomposePar="<<decomposePar_time<< std::endl;
#endif   
	  //for(i = 0; i < n; i++) data[i] = tmpPtr[i]; 
	  for(int j = 0; j < n; j++) data[j] = tmpPtr[j]; 
	  
	  //std::cout <<"before exit multi_mat_transpose" << std::endl;
	  my_free(accum);
	  my_free(tmpPtr);
	  //std::cout <<"after free in multi_mat_transpose" << std::endl;
	  //#ifdef MATTEOTRAN2
	}
	//#endif 
	//for(int j = 0; j < n; j++) std::cout<<data[j]<<", ";
	//std::cout<<"----end of out data----" <<std::endl;
  }
  
  
  /**-----------------------------------------------------------
   * decompose:
   * @N: number of dimensions of a matrix M.
   * @dm: index for the dm_th dimension.
   * @dims: The dimensions vector (dims[0] is not used). 
   *        dims[i] keep the size on dimension i.
   * @data: entries of the matrix in a consecutive array.
   * 
   * We transpose matrix M from dimension 'dm' to dimension 1. 
   * Return value: 
   **/
  void
  decompose(sfixn N, sfixn dm, sfixn * accum, sfixn * dims, sfixn * tmpPtr, sfixn * data){

	if (N == dm){
	  int ac = accum[N+1] / dims[1];
	  sfixn acN = ac/dims[N];
	  //cilk_for(int k=0; k<acN; k++) { //for works
	  //for(int i = 0; i < dims[1]; i++){ 
	  //change to 
	  for(int k=0; k<acN; k++) { 
		for(int i = 0; i < dims[1]; i++){ 
		  for(int j = 0; j < dims[dm]; j++){
			tmpPtr[j + i*ac + k*dims[N]] = data[j*accum[dm] + i + k*dims[1]];
		  }
		}
	  }
	  return;
	}
	
	sfixn acN = accum[N+1]/accum[N];
	for(int k=0; k<acN; k++) {
	  SERBPAS::decompose(N-1, dm, accum, dims, tmpPtr+k*accum[N], data+k*accum[N]);
	}
	
  } 
  

  /**---------------------------------------------------------
   * getDenseSiz:
   * @N: number of variables.
   * @poly: A C-Cube polynomial.
   * @buszsN: BUSZSI('poly', 'N').
   * @data: The data field of 'poly'.
   * @cumN: CUMI('poly', 'N').
   *
   * Data size of the C-Cube polynomial 'poly' minus all leading zeros in it.
   *
   * Return value: SIZ(poly)- number of all nonleading zeros in the least variable. 
   **/
  sfixn getDenseSiz(sfixn N, preFFTRep *poly, sfixn buszsN, sfixn *data, sfixn cumN){
    sfixn deg, siz=0;  
    
    deg = shrinkDeg(buszsN, data, cumN);       
    
    if(N==1) return deg +1;
    
    if(deg>0) return (deg+1)*CUMI(poly, N);
    
    if(deg ==0){
  
      siz = SERBPAS::getDenseSiz(N-1, poly, BUSZSI(poly, N-1), data, CUMI(poly, N-1));
      return siz; 
    }
    
    return -1;
  }

  /**
   * InitKroFFTRep:
   * @kPtr: (output) A data structure to prepare the Kronecker-based FFT multiplication.
   * @resDgs: The degree vector of the result polynomial.
   * @N: The number of variables.
   * @M: The number of polynomials to multiply (poly1* poly2 * ... * polyM.
   * @pPtr: The information of the prime p.
   * 
   * Create 'kPtr' -- a data structure of type KroFFTRep.
   * This data structure will be used for Kronecker-based FFT multiplicaiton.
   *
   * Return value: 
   **/
  void
  InitKroFFTRep(KroFFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr){
	register int j;
	N(kPtr)=N;     
	M(kPtr)=M;
	ES(kPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn)); 
	//
	DIMS(kPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
	//
	CUM(kPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
	//
	CUMI(kPtr, 1)= 1; SIZ(kPtr)=1;

	for(j=1;j<=N;j++){
	  ESI(kPtr, j)=logceiling(resDgs[j]+1);
	  DIMSI(kPtr, j)=1<<(ESI(kPtr, j));
	  SIZ(kPtr)*=DIMSI(kPtr, j); 
	  if(j>=2){CUMI(kPtr, j)=CUMI(kPtr, j-1)*DIMSI(kPtr, j-1);}
	}
//N is small --------------------------------------------
// 	cilk_for(int i=1; i<=N; i++){
// 	  ESI(kPtr, i) = logceiling(resDgs[i]+1);
// 	}
// 	cilk_for(int i=1; i<=N; i++){
// 	  DIMSI(kPtr, i) = 1<<(ESI(kPtr, i));
// 	}
// 	for(j=1;j<=N;j++){
// 	  SIZ(kPtr) *= DIMSI(kPtr, j);
// 	  if(j>=2){
// 		CUMI(kPtr, j) = CUMI(kPtr, j-1)*DIMSI(kPtr, j-1);
// 	  }
// 	}

	DATS(kPtr)=(sfixn **)my_calloc(M,sizeof(sfixn *));
	//
	for(j=0;j<M;j++){
	  DATSI(kPtr, j)=(sfixn * )my_calloc(SIZ(kPtr),sizeof(sfixn));
	  //
	}

	KN(kPtr)=1; KE(kPtr)=0;

	for(j=1;j<=N;j++) {
	  KN(kPtr)*=DIMSI(kPtr, j); 
	  KE(kPtr)+=ESI(kPtr, j);
	}

	KROOTS(kPtr)=NULL;
	kPtr->Defsize=SIZ(kPtr);
  }

  // free a Kronecker data-structure
  void freeKroFFTRep(KroFFTRep * x){
	//int i;
	if(KROOTS(x) != NULL) {my_free(KROOTS(x)); KROOTS(x)=NULL;}
	if(CUM(x) !=NULL) {my_free(CUM(x)); CUM(x)=NULL;}
	if(ES(x) !=NULL) {my_free(ES(x)); ES(x)=NULL;}
	if(DIMS(x) !=NULL) {my_free(DIMS(x)); DIMS(x)=NULL;}
	if(DATS(x) != NULL){
	  for(int i=0; i<M(x); i++) 
		if(DATSI(x, i) !=NULL) {my_free(DATSI(x, i)); DATSI(x, i)=NULL;}
	  my_free(DATS(x));
	  DATS(x)=NULL;}
  }

  //-------------------------------------
  void
  InitKroTFTRep(KroTFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr){
    register int j;
    N(kPtr)=N;     
    M(kPtr)=M;
    ES(kPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn)); 
    //
    DIMS(kPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    //
    LS(kPtr)=(sfixn *)my_calloc(N+1,sizeof(sfixn));
    
    CUM(kPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    //
    CUMI(kPtr, 1)= 1; SIZ(kPtr)=1;
    for(j=1;j<=N;j++){
      ESI(kPtr, j)=logceiling(resDgs[j]+1);
      DIMSI(kPtr, j)=1<<(ESI(kPtr, j));
      LSI(kPtr, j)=resDgs[j]+1;
      SIZ(kPtr)*=LSI(kPtr, j); 
      if(j>=2){CUMI(kPtr, j)=CUMI(kPtr, j-1)*LSI(kPtr, j-1);}
    }
    
    DATS(kPtr)=(sfixn **)my_calloc(M,sizeof(sfixn *));
    //
    for(j=0;j<M;j++){
      DATSI(kPtr, j)=(sfixn * )my_calloc(SIZ(kPtr),sizeof(sfixn));
      //
    }
    KN(kPtr)=1; KE(kPtr)=0;
    for(j=1;j<=N;j++) {KN(kPtr)*=DIMSI(kPtr, j); KE(kPtr)+=ESI(kPtr, j);}
    KROOTS(kPtr)=NULL;
    kPtr->Defsize=SIZ(kPtr);
  }
  
  //----------------------------
  void freeKroTFTRep(KroTFTRep * x){
    int i;
    if(KROOTS(x) != NULL) {my_free(KROOTS(x)); KROOTS(x)=NULL;}
    if(CUM(x) !=NULL) {my_free(CUM(x)); CUM(x)=NULL;}
    if(LS(x) !=NULL) {my_free(LS(x)); LS(x)=NULL;}
    if(ES(x) !=NULL) {my_free(ES(x)); ES(x)=NULL;}
    if(DIMS(x) !=NULL) {my_free(DIMS(x)); DIMS(x)=NULL;}
    if(DATS(x) != NULL){
      for(i=0; i<M(x); i++) 
		if(DATSI(x, i) !=NULL) {my_free(DATSI(x, i)); DATSI(x, i)=NULL;}
      my_free(DATS(x));
      DATS(x)=NULL;}
  }

  //cutoff=8
  // return 0, indicating using Kronecker-FFT.
  // return 1, indicating using multiD-TFT.
  int
  forcutoff_Multi_TFT_FFT(sfixn N, sfixn *dgs1, sfixn *dgs2, int cutoff){
	int i;
	for(i=1; i<=N; i++){
	  if(dgs1[i]<cutoff) return 0;
	  if(dgs2[i]<cutoff) return 0;
	}
	return 1;
  }
  
  /**
   * InitResPoly:
   * @rPtr: (output) A C-Cube polynomial.
   * @N: Number of variables.
   * @p1dgs: A degree vector for some C-Cube polynomial p1.
   * @p2dgs: A degree vector for some C-Cube polynomial p2.
   * 
   * Initialized a 'clean' (all zero-coefficients) C-Cube polynomial whoes
   * data space can exactly keep the product of p1 and p2.
   *
   * Return value: 
   **/
  void
  InitResPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs, sfixn * p2dgs){  
	register int j;
    N(rPtr)=N;
    BUSZS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUTS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    //
    SIZ(rPtr)=1;
    CUM(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    //
    CUMI(rPtr, 1)= 1;
    for(j=1;j<=N;j++){
      BUSZSI(rPtr, j)=p1dgs[j]+p2dgs[j];
      CUTSI (rPtr, j)=   BUSZSI(rPtr, j);
      SIZ(rPtr)=SIZ(rPtr)*(BUSZSI(rPtr, j)+1);
      if(j>=2){
		CUMI(rPtr, j)=CUMI(rPtr, j-1)*(BUSZSI(rPtr, j-1)+1);
      }
    }
    OFST(rPtr)=0;
    DAT(rPtr)=(sfixn * )my_calloc( SIZ(rPtr),sizeof(sfixn));
	//DAT(rPtr)=NULL;
    //can be delayed to after calculating the result?
    DFN(rPtr)=N(rPtr);
    DFSIZ(rPtr)=SIZ(rPtr);
    DEFDAT(rPtr)=DAT(rPtr);
    //BUSZSI(rPtr, 0)=BUSZSI(rPtr, N);
  }
  
  
  // copy elements 0..d.
  /**
   * copyVec_0_to_d:
   * @n:
   * @desV:
   * @srcV: 
   * assignment:
   * desV[0..d] = srcV[0..d]. 
   * Return value: 
   **/
  void copyVec_0_to_d(int d, sfixn * desV, sfixn * srcV){
	//register int i;
	for(int i=0; i<=d; i++) { desV[i] = srcV[i]; }  //cilk_for slows down
  }
  
  // copy elements 1..n.
  /**
   * copyVec_1_to_n:
   * @n:
   * @desV:
   * @srcV: 
   * assignment:
   * desV[1..n] = srcV[1..n]. 
   * Return value: 
   **/
  void copyVec_1_to_n(int n, sfixn * desV, sfixn * srcV){
	//register int i;
	for(int i=1; i<=n; i++) desV[i]=srcV[i];  
  }

  //==========================================================
  // compare two vectors. equal return 1, otherwise return 0.
  //==========================================================
  /**
   * compareVec:
   * @deg: degree of both 'vec1' and 'vec2'.
   * @vec1: a vector.
   * @vec2: a vector.
   * To compare two input vectors 'vec1' and 'vec2'.
   * if they are equal returns 1 otherwise return 0.
   * Return value: 
   **/
  int
  compareVec(sfixn deg, sfixn * vec1, sfixn * vec2){
	register int i;
	for(i=0; i<=deg; i++){
	  if(vec1[i]!=vec2[i]) {
		std::cout << "----diff at index: " << i << std::endl;
		std::cout << "----res1 here: " << vec1[i] << std::endl;
		std::cout << "----res2 here: " << vec2[i] << std::endl;
		return(0);
	  };
	}
	return 1;
  }
  
  // 1 -- yes they are euqal.
  // 0 -- No they are NOT equal.
  /**
   * EX_IsEqualPoly:
   * @Ptr1: a C-Cube polynomial 'P1'.
   * @Ptr2: a C-Cube polynomial 'P2'.
   * 
   * To compare if two polynomials are equal.
   *
   * Return value: if 'P1' is equal to 'P2', then return 1. Otherwise return 0.
   **/
  int
  EX_IsEqualPoly(preFFTRep * Ptr1, preFFTRep * Ptr2){
	if(N(Ptr1) != N(Ptr2)) return 0;
	if (! SERBPAS::compareVec(N(Ptr1), BUSZS(Ptr1), BUSZS(Ptr2))) return 0;
	if (! SERBPAS::compareVec(SIZ(Ptr1)-1, DAT(Ptr1), DAT(Ptr2))) return 0;
	return 1;
  }
  
  
  // suppose y>x
  // m rows, n columns
  // dgs is always the smaller one.  
  void fromtofftRepMultiD(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs){
	//int i;
	int d;
	//int tmpRes=0, tmpCoeffs=0;
	if(N==0){
	  res[0]=coeffs[0];
	  return;
	}
		
	if(N==1){
	  d=shrinkDeg(dgs[1], coeffs, 1);
	  for(int i=0; i<=d; i++){
		res[i]=coeffs[i];
	  }
	  return;
	}
	
	
	d=shrinkDeg(dgs[N], coeffs, ccum[N]);
	//   for(i=0; i<=d; i++){
	//     tmpCoeffs=i*ccum[N];
	//     tmpRes=i*rccum[N];
	//     fromtofftRepMultiD(N-1, rccum, res+tmpRes, ccum, dgs, coeffs+tmpCoeffs);
	//   }
	for(int i=0; i<=d; i++){
	  SERBPAS::fromtofftRepMultiD(N-1, rccum, res+i*rccum[N], ccum, dgs, coeffs+i*ccum[N]);
	}
  }
  
  /**
   * reverseUni:
   * @deg: 'deg'+1 is the size for both 'vec1', and 'vec2'.
   * @vec1: a destination vector.
   * @vec2: a source vector.
   * 
   * copy the data in 'vec2' in a reverse order into 'vec1'.
   * 
   * Return value: 'vec1'
   **/
  sfixn *
  reverseUni(sfixn deg, sfixn * vec1, sfixn * vec2){
	for(int i=0; i<=deg; i++){
	  vec1[i] = vec2[deg-i];
	}
	return vec1;
  }

  //==========================================================
  // clean vectors.
  //==========================================================
  // note: clean the whole vector.
  void cleanVec(sfixn deg, sfixn * cof){
	for(int i=0; i<=deg; i++) cof[i]=0;
  }

  /**
   * InitOneRevInv:
   * @N: number of variables for the output polynomial
   * @tRIPtr: (output) buffer for the inverse (in reversed order)  
   *          of the N-th polynomial in a triangular set 
   *          with partial degrees strictly bounded by `bounds`
   * @bounds: vector 
   * @di: the truncation degree for the inverse 
   * 
   *   Initializes a buffer for a modular inverse computation
   *   of a polynomial modulo a triangular set.
   * 
   * Return value: 
   **/   
  // degree of RevInv is always power of 2 which is slightly more than enough.
  void
  InitOneRevInv(sfixn N, preFFTRep * tRIPtr, sfixn * bounds, sfixn di){
    register int j;

    sfixn eidi2 = di ;
	if(eidi2<=0) { 
	  std::cout<<"The truncation degree for the inverse <= 0."<<std::endl; 
	  exit(1);
	}

    N(tRIPtr)=N;
    BUSZS(tRIPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUTS(tRIPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUM(tRIPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn)); 
    CUMI(tRIPtr, 1)=1;
    for(j=1; j<=N; j++) {
      BUSZSI(tRIPtr, j)=bounds[j];
      CUTSI(tRIPtr, j)=bounds[j];
      if(j>=2){
		CUMI(tRIPtr, j)=CUMI(tRIPtr, j-1)*(BUSZSI(tRIPtr, j-1)+1);}
    }
	
    BUSZSI(tRIPtr, N)=(1<<((logceiling(eidi2))))-1;

    CUTSI(tRIPtr, N)=BUSZSI(tRIPtr, N);
    SIZ(tRIPtr)=CUMI(tRIPtr, N)*(BUSZSI(tRIPtr, N)+1);
    OFST(tRIPtr)=0;
    DAT(tRIPtr)=(sfixn * )my_calloc(SIZ(tRIPtr),sizeof(sfixn));
    DFN(tRIPtr)=N(tRIPtr);
    DFSIZ(tRIPtr)=SIZ(tRIPtr);
    DEFDAT(tRIPtr)=DAT(tRIPtr);
  }
  
  /**
   * reverseMulti:
   * @deg: ('deg'+1)*'sizOfCoef' is the size for both 'outVec', and 'inVec'.
   * @sizOfCoef: An integer number. 
   * @outVec: A destination vector.
   * @inVec: A source vector.
   * 
   * Copy 'deg'+1 data blocks from 'outVec' in a reverse order into 'outVec'.
   * each data blocks have size 'sizOfCoef'.
   * Return value: 'outVec'
   **/
  sfixn *
  reverseMulti(sfixn deg, sfixn sizOfCoef, sfixn * outVec, sfixn * inVec){
	register sfixn i;
	sfixn end=deg*sizOfCoef;
	for(i=0; i<=end; i+=sizOfCoef){
	  for(int j=0; j<sizOfCoef; j++) outVec[end-i+j]=inVec[i+j];
	}
	return outVec;
  }

  //=================================================================
  // reverse poly in f[x1,...,xn][y].
  //=================================================================
  // The Input multivariate polynomial f[x1,...,xn][y] encoded in a linear array.
  // deg is the degree of f in y.
  // sizOfCoef is the size of coefficient of f in y.
  // degree >=0, sizOfCoef >0;
  sfixn *
  reverseMulti_1(sfixn deg, sfixn sizOfCoef, sfixn * vec){
	register int i, j;
	sfixn * tmpCoef=(sfixn * )my_malloc(sizOfCoef*sizeof(sfixn));
	sfixn * tmp1Ptr=vec, * tmp2Ptr=vec+sizOfCoef*deg;
	for(i=0; i<((deg+1)/2); i++){
	  for(j=0; j<sizOfCoef; j++) tmpCoef[j]=tmp1Ptr[j];
	  for(j=0; j<sizOfCoef; j++) tmp1Ptr[j]=tmp2Ptr[j];
	  for(j=0; j<sizOfCoef; j++) tmp2Ptr[j]=tmpCoef[j];
	  tmp1Ptr+=sizOfCoef;
	  tmp2Ptr-=sizOfCoef;
	}
	my_free(tmpCoef);
	return vec;
  }

  /**inline before
   * PolyCleanData:
   * @prePtr:  A C-Cube polynomial.
   * 
   * make prePtr a zero polynomial.
   * Return value: 
   **/
  void PolyCleanData(preFFTRep * prePtr){
	sfixn size = SIZ(prePtr);
	
	for(int i=0; i<size; i++){
	  DATI(prePtr, i)=0;
	}
  }  
  
  // cannot decrease dimension, only es, dims, siz, aacum.
  void
  decreaseKroFFTRep(KroFFTRep * kPtr, sfixn * resDgs){
	register int i,j;
	CUMI(kPtr, 1)= 1; SIZ(kPtr)=1;
	for(j=1;j<=N(kPtr);j++){
	  ESI(kPtr, j)=logceiling(resDgs[j]+1);
	  DIMSI(kPtr, j)=1<<(ESI(kPtr, j));
	  SIZ(kPtr)*=DIMSI(kPtr, j); 
	  if(j>=2) {
		CUMI(kPtr, j)=CUMI(kPtr, j-1)*DIMSI(kPtr, j-1);
	  }
	}
	
	KN(kPtr)=1; KE(kPtr)=0;
	for(j=1;j<=N(kPtr);j++) {
	  KN(kPtr)*=DIMSI(kPtr, j); 
	  KE(kPtr)+=ESI(kPtr, j);
	}
    
	//
	for(i=0;i<M(kPtr);i++){
	  for(j=0;j<kPtr->Defsize;j++){
		(DATSI(kPtr, i))[j]=0;
	  }
	}
	if(KROOTS(kPtr)!=NULL) {
	  my_free(KROOTS(kPtr)); 
	  KROOTS(kPtr)=NULL;
	}
  }

  //------------------------------------------------
 void KroneckerCleanData(KroFFTRep * kPtr){
  for(int i=0; i<M(kPtr); i++){
    for(int j=0; j<kPtr->Defsize; j++){
      (DATSI(kPtr, i))[j]=0;
	}
  }
}
 
  //=====================================================
  //  copying data from one dense multivariate polynomial
  //  in coeffs to the one in res.
  //=====================================================
  void fromtofftRep(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs){
	
	sfixn d;

	if(N==0){
	  res[0]=coeffs[0];
	  return;
	}

	d=dgs[N];
	if(N==1){
	  for(int i=0; i<=d; i++){ 
		res[i]=coeffs[i];
	  }
	  return;
	}
// 	for(i=0; i<=d; i++){
// 	  tmpCoeffs=i*ccum[N];
// 	  tmpRes=i*rccum[N];
// 	  fromtofftRep(N-1, rccum, res+tmpRes, ccum, dgs, coeffs+tmpCoeffs); 
// 	}
	for(int i=0; i<=d; i++){
	  SERBPAS::fromtofftRep(N-1, rccum, res+i*rccum[N], ccum, dgs, coeffs+i*ccum[N]);
	}

  }

  void InitOneReducedPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs){  
	register int j; 
	
    N(rPtr)=N;
	
    BUSZS(rPtr)=(sfixn * )my_calloc(N+1, sizeof(sfixn));
	
    CUTS(rPtr)=(sfixn * )my_calloc(N+1, sizeof(sfixn));
    SIZ(rPtr)=1;
    CUM(rPtr)=(sfixn * )my_calloc(N+1, sizeof(sfixn));
    CUMI(rPtr, 1)= 1;
	
    for(j=1;j<=N;j++){
	  
      BUSZSI(rPtr, j)=p1dgs[j];
      CUTSI (rPtr, j)=p1dgs[j];   
      SIZ(rPtr)=SIZ(rPtr)*(BUSZSI(rPtr,j)+1);
      if(j>=2){
		CUMI(rPtr, j)=CUMI(rPtr, j-1)*(BUSZSI(rPtr,j-1)+1);
      }
    }

    OFST(rPtr)=0;
    DAT(rPtr)=(sfixn * )my_calloc( SIZ(rPtr),sizeof(sfixn));  
    //
    //BUSZSI(rPtr, 0)=BUSZSI(rPtr, N);
    DFN(rPtr)=N(rPtr);
    DFSIZ(rPtr)=SIZ(rPtr);
    DEFDAT(rPtr)=DAT(rPtr);
  }
  
  void
  initCopyOneRevInv(preFFTRep * outPtr, preFFTRep * inPtr){
    register int j;
    N(outPtr)=N(inPtr);
    BUSZS(outPtr)=(sfixn * )my_calloc(N(inPtr)+1,sizeof(sfixn)); 
    CUTS(outPtr)=(sfixn * )my_calloc(N(inPtr)+1,sizeof(sfixn)); 
    CUM(outPtr)=(sfixn * )my_calloc(N(inPtr)+1,sizeof(sfixn)); 
    for(j=1;j<=N(inPtr);j++){ BUSZSI(outPtr, j)=BUSZSI(inPtr, j);
                              CUTSI(outPtr, j)=CUTSI(inPtr, j);
                              CUMI(outPtr, j)=CUMI(inPtr, j);}
    SIZ(outPtr)=SIZ(inPtr);
    OFST(outPtr)=OFST(inPtr);
    DAT(outPtr)=(sfixn * )my_calloc(SIZ(outPtr),sizeof(sfixn));
    DFN(outPtr)=N(outPtr);
    DFSIZ(outPtr)=SIZ(outPtr);
    DEFDAT(outPtr)=DAT(outPtr);
    //BUSZSI(outPtr, 0)=BUSZSI(outPtr, N(outPtr));
  }

  // this in-place fun suppose Ptr1 is the larger buffer on all dimensions.
  /**
   * subPoly_1:
   * @N: Number of variables in 'Ptr1' and 'Ptr2'.
   * @Ptr1: A C-Cube polynomial.
   * @Ptr2: A C-Cube polynomial.
   * @p: A prime number.
   * 
   * Suppose 'Ptr1' has the same dimension and larger(or equal) size on each dimension.
   * Compute the difference of them.
   * Return value: Ptr1 = Ptr1 - Ptr2;
   **/
  void subPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p){
	SERBPAS::subPoly_inner_1(N, CUM(Ptr1), BUSZS(Ptr2), CUM(Ptr2), DAT(Ptr1), DAT(Ptr2), p);
  }
  
  // we will use the smaller dgs which are 2nd. 
  // but accums should used prepectively.
  void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p)
  {
	int i, offset1=0, offset2=0;
	if(N==1){
	  for(i=0;i<=dgs2[1];i++) data1[i]=SubMod(data1[i],data2[i],p);
	  return;}
	for(i=0; i<=dgs2[N]; i++){
	  offset1=accum1[N]*i;
	  offset2=accum2[N]*i;
	  SERBPAS::subPoly_inner_1(N-1, accum1, dgs2, accum2, data1+offset1, data2+offset2, p);
	} 	
  }

  //--------------------------------
  void subEqDgPoly_inner_1
  (sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int selector)
  {
	int i, offset=0;
	if(N==1){
	  if(selector==1){
		for(i=0;i<=dgs[1];i++) data1[i]=SubMod(data1[i],data2[i],p);}
	  else{
		for(i=0;i<=dgs[1];i++) data2[i]=SubMod(data1[i],data2[i],p);}
	  return;}

	for(i=0; i<=dgs[N]; i++){
	  offset=accum[N]*i;
	  subEqDgPoly_inner_1(N-1, dgs, accum, data1+offset, data2+offset, p, selector);
	} 
  }

  //-------------------------------------
  void subEqDgPoly_1
  (sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p, int selector)
  {
	subEqDgPoly_inner_1(N, BUSZS(Ptr1), CUM(Ptr1), DAT(Ptr1), DAT(Ptr2), p, selector);
  }
    //---------------------------------
  void addEqDgPoly_inner
  (sfixn N, sfixn *dgs, sfixn *accum, sfixn *data1, sfixn *data2, sfixn p)
  {
	int i, offset=0;
	if(N==1){
	  for(i=0;i<=dgs[1];i++) data1[i]=AddMod(data1[i],data2[i],p);
	  return;}
	for(i=0; i<=dgs[N]; i++){
	  offset=accum[N]*i;
	  addEqDgPoly_inner(N-1, dgs, accum, data1+offset, data2+offset, p);
	} 
  }

  /**
   * addEqDgPoly_1:
   * @N: Number of variables in 'Ptr1' and 'Ptr2'.
   * @Ptr1: A C-Cube polynomial.
   * @Ptr2: A C-Cube polynomial.
   * @p: A prime number.
   * 
   * Suppose 'Ptr1' and 'Ptr2' has the same dimension and size.
   * Compute the sum of them.
   * Return value: Ptr1 = Ptr1 + Ptr2;
   **/
  void addEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p){
	addEqDgPoly_inner(N, BUSZS(Ptr1), CUM(Ptr1), DAT(Ptr1), DAT(Ptr2), p);
  }
  
  /**
   * copyPolyPointers:
   * @D: Destination C-cube polynomial.
   * @S: Source C-cube polynomial.
   * 
   * Assume 'D' has type of preFFTRep *, but its sub-fields 
   * have not been initialized.  
   * This function copies S's fields to 'D'.
   *
   * Return value: 
   **/
  void 
  copyPolyPointers(preFFTRep *D, preFFTRep *S)
  {
	D->N=S->N; 
	D->defN=S->defN;
	D->bufSizs=S->bufSizs;
	D->cuts=S->cuts; 
	D->accum=S->accum;
	D->size=S->size;  
	D->defSize=S->defSize;
	D->offset=S->offset;
	D->data=S->data;
	D->tmpData=S->tmpData; 
	D->defData=S->defData; 
  }

  /**
   * getCoefMulti:
   * @N: Number of variables.
   * @co: (output) the coefficient of 'f' in degree `i` w.r.t. `X_N`.
   * @f: A C-Cube polynomial in `N` variables
   * @i: An index.
   * Make a copy of the i-th coefficient of 'f' and save the copy in co.
   * Return value: The copy of i-th coefficient of 'f'. 
   **/
  preFFTRep *
  getCoefMulti(sfixn N, preFFTRep* co, preFFTRep* f, sfixn i){
	decreaseOneDim(f);
	backupData(f); 
	nextMCoefData(f,N,i);
	SERBPAS::fromtofftRep(N-1,  CUM(co), DAT(co), CUM(f), BUSZS(f), DAT(f));
	increaseOneDim(f);
	restoreData(f);
	return co;
  }

  /**
   * initTriRevInvSet:
   * @dgs: A sfixn vector which keeps the truncation degrees for the 
   *         modular inverses. Except for normal form computations,
   *        these should be paritial degrees of the input triangualr set 'tPtr' .
   * @N: Number of variables.
   * @tRevInvPtr: (output) The buffers for the inverses (reverse-ordered) 
   *            of the input triangular set 'tPtr'.
   * @tPtr: A triangular set in dimension-zero.
   * 
   * Create the buffers for the inverses (reverse-ordered) of the intput
   *         triangular set 'tPtr'.
   *
   * Return value: 
   **/
  void initTriRevInvSet(sfixn *dgs, sfixn N, TriRevInvSet *tRevInvPtr, TriSet * tPtr){
	register int i;
	N(tRevInvPtr)=N;
	EXST(tRevInvPtr)=(sfixn *)my_calloc(N+1,sizeof(sfixn) );
	NLB(tRevInvPtr)=(sfixn *)my_calloc(N+1,sizeof(sfixn) );
	NSB(tRevInvPtr)=(sfixn *)my_calloc(N+1,sizeof(sfixn) );
	ELEM(tRevInvPtr)=(preFFTRep **)my_calloc((N+1),sizeof(preFFTRep *) );
	for(i=1; i<=N; i++){
	  ELEMI(tRevInvPtr, i)=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	  
	  SERBPAS::InitOneRevInv(i, ELEMI(tRevInvPtr, i), BDS(tPtr), dgs[i]); 
	  
	  (NLB(tRevInvPtr))[i]= BUSZSI(ELEMI(tRevInvPtr, i), i);
	  (NSB(tRevInvPtr))[i]=BDSI(tPtr, i)-1;
	}
  }
  
  //------------------------------
  TriSet * Ex_initRandomMonicDenseTriSet(sfixn N, sfixn * dgs, MONTP_OPT2_AS_GENE * pPtr){
	TriSet *tPtr;
	tPtr=(TriSet *)my_calloc(1,sizeof(TriSet));
	SERBPAS::initRandomDenseTriSet(N, dgs, tPtr, pPtr);
	return tPtr;
  }
  

  /**
   * initRandomDenseTriSet:
   * @N: The number of variables.
   * @dgs: the partial degrees in the output triangular set
   * @tPtr: (output) the triangular set.
   * @pPtr: The information of the prime.
   * 
   * To create a random monic, zero-dimensional triangular set,
   * with prescribed partial degrees.
   * Return value: The newly created random triangular set. 
   **/
  void initRandomDenseTriSet( sfixn N, sfixn * dgs, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr){
	int i;
	N(tPtr)=N;
	ELEM(tPtr)=(preFFTRep **)my_calloc((N+1), sizeof(preFFTRep *) );
	BDS(tPtr)=(sfixn *) my_calloc((N+1), sizeof(sfixn));
	srand(getSeed());
	for(i=1;i<=N;i++)  BDSI(tPtr, i) = dgs[i-1]-1;
	
	if(MODPN::checkDgsOfST(N, tPtr)==-1) {
	  std::cout<<"Invalid number of variables"<<std::endl;
	  exit(1);
	}
	
	for(i=1; i<=N; i++){
	  ELEMI(tPtr, i)=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	  SERBPAS::InitOneMonicRandomPolys(BDS(tPtr), i, ELEMI(tPtr, i), pPtr, i); 
	  for(int k=1; k<=i; k++){
		std::cout<< BUSZSI ( ELEMI(tPtr, i), k ) <<";";
	  }
	  std::cout<<std::endl;
	}	
  }
  

  /**
   * EX_initRandomTriSet:
   * @N: Number of variables.
   * @dgbound: A sfixn which is a strict upper bound for all partial
   *           degrees in the output triangular set
   * @pPtr: Information of the prime number.
   * Create an randome Lazard trianuglar set in dimension zero
   *       with `N` variables and partial degrees strictly bounded by `degbound`
   * Return value:  an randome Lazard trianuglar set in dimension 
   **/
  TriSet *
  EX_initRandomTriSet( sfixn N, sfixn dgbound, MONTP_OPT2_AS_GENE * pPtr){
	TriSet *tPtr;
	tPtr=(TriSet *)my_calloc(1,sizeof(TriSet));
	SERBPAS::initRandomTriSet(N, dgbound, tPtr, pPtr);
	return tPtr;
  }
  
  // Subroutine of the above function
  void initRandomTriSet( sfixn N, sfixn dgbound, TriSet * tPtr,  MONTP_OPT2_AS_GENE * pPtr){
	int i;
	N(tPtr)=N;
	NMLZ(tPtr)=0;
	ELEM(tPtr)=(preFFTRep **)my_calloc((N+1),sizeof(preFFTRep *) );
	BDS(tPtr)=(sfixn *) my_calloc((N+1),sizeof(sfixn));
	srand(getSeed());
	//for(i=1;i<=N;i++) {while(! BDSI(tPtr,i) ) BDSI(tPtr, i)=rand()%dgbound;}
	if(MODPN::checkDgsOfST(N, tPtr)==-1) {
	  std::cout<<"Invalid number of variables"<<std::endl;
	  exit(1);
	}
	for(i=1;i<=N;i++) {while( BDSI(tPtr,i)<=0 ) BDSI(tPtr, i)=rand()%dgbound;}
	for(i=1; i<=N; i++){
	  ELEMI(tPtr, i)=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	  SERBPAS::InitOneMonicRandomPolys(BDS(tPtr), i, ELEMI(tPtr, i), pPtr, i);
	  for(int k=1; k<=i; k++){
		std::cout<< BUSZSI ( ELEMI(tPtr, i), k ) <<";";
	  }
	  std::cout<<std::endl;
	}
	
  }

  // Subroutine for building a monic random poly
  void
  InitOneMonicRandomPolys(sfixn * bounds, sfixn N, preFFTRep * p1Ptr,  MONTP_OPT2_AS_GENE * pPtr, sfixn seed){
	register int j;
    N(p1Ptr)=N;
    BUSZS(p1Ptr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUTS(p1Ptr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUM(p1Ptr)=(sfixn * )my_calloc(N+1,sizeof(sfixn)); 
    CUMI(p1Ptr, 1)= 1;
    SIZ(p1Ptr)=1;
    OFST(p1Ptr)=0; 
    srand(getSeed());
    for(j=1;j<=N;j++){
      BUSZSI(p1Ptr, j)=bounds[j];
      CUTSI(p1Ptr, j)=bounds[j];
      if(j==N) { BUSZSI(p1Ptr, j)=bounds[j]+1;
		CUTSI(p1Ptr, j)=bounds[j]+1;
      }
      SIZ(p1Ptr)*=BUSZSI(p1Ptr, j)+1;
      if(j>=2){
		CUMI(p1Ptr, j)=CUMI(p1Ptr, j-1)*(BUSZSI(p1Ptr, j-1)+1);
      }
	}
    DAT(p1Ptr)=(sfixn * )my_calloc(SIZ(p1Ptr),sizeof(sfixn));
    SERBPAS::randomMonicPoly(p1Ptr, pPtr->P);
    DFN(p1Ptr)=N(p1Ptr);
    DFSIZ(p1Ptr)=SIZ(p1Ptr);
    DEFDAT(p1Ptr)=DAT(p1Ptr);
    //BUSZSI(p1Ptr, 0)=BUSZSI(p1Ptr, N);
  }
  
  //---------------------------
  void randomMonicPoly(preFFTRep * Ptr, sfixn p){
	SERBPAS::randomMonicPoly_inner(N(Ptr), N(Ptr), BUSZS(Ptr), CUM(Ptr), DAT(Ptr), p);  
  }

  // Subroutine for the above function
  void randomMonicPoly_inner(sfixn N, sfixn N1, sfixn * dgs, sfixn * accum, sfixn * data, sfixn p){
	int i, offset=0;
	if(N1==1){
	  srand(getSeed() );
	  for(i=0; i<=dgs[N1]; i++){ data[i]=(rand()%p); }
	  if (! data[dgs[N1]]) data[dgs[N1]]=1;
	  if(N==1) data[dgs[N]]=1;
	  return;
	}
	
	if(N==N1){
	  for(i=0; i<dgs[N1]; i++){
		offset=accum[N1]*i;
		randomMonicPoly_inner(N, N1-1, dgs, accum, data+offset, p);
	  }
	  (data+accum[N]*dgs[N])[0]=1;}
	
	else{
	  for(i=0; i<=dgs[N1]; i++){
		offset=accum[N1]*i;
		SERBPAS::randomMonicPoly_inner(N, N1-1, dgs, accum, data+offset, p);
	  }
	  
	}
  }

  /**
   * EX_randomPoly:
   * @N: Number of variables.
   * @dgs: A partial-degree vector of size N.
   * @p: A prime number. 
   *
   * To initialize a C-Cube polynomial whoes partial degrees are defined in 'dgs'.
   * Then generate random coefficients in Z/pZ for this C-Cube polynomial.
   * Return value: The newly created C-Cube polynomial
   **/
  preFFTRep *
  EX_randomPoly(sfixn N, sfixn * dgs, sfixn p){
	preFFTRep *poly;
	poly = SERBPAS::EX_InitOnePoly(N, dgs);
	SERBPAS::randomPoly(poly, p);
	return poly;
  }
  
  /**
   * EX_InitOnePoly:
   * @N: Number of variables.
   * @dgs: A partial-degree vector of size N.
   * 
   * To initialize a C-Cube polynomial whoes partial 
   * degrees are defined in 'dgs'.
   *
   * Return value: A newly created C-Cube polynomial
   **/
  preFFTRep *
  EX_InitOnePoly(sfixn N, sfixn * dgs){
	preFFTRep *poly= (preFFTRep *) my_malloc(sizeof(preFFTRep));
	SERBPAS::InitOnePoly(poly, N, dgs);
	return poly;
  }
  
  // Subroutine for the above one
  void InitOnePoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs){  
	register int j; 
    N(rPtr)=N;
    BUSZS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUTS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    SIZ(rPtr)=1;
    CUM(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUMI(rPtr, 1)= 1;
    for(j=1;j<=N;j++){
      BUSZSI(rPtr, j)=p1dgs[j];
      CUTSI (rPtr, j)=p1dgs[j];   
      SIZ(rPtr)=SIZ(rPtr)*(BUSZSI(rPtr,j)+1);
      if(j>=2){
		CUMI(rPtr, j)=CUMI(rPtr, j-1)*(BUSZSI(rPtr,j-1)+1);
      }
    }
    OFST(rPtr)=0;
    DAT(rPtr)=(sfixn * )my_calloc( SIZ(rPtr),sizeof(sfixn));  
    //
    //BUSZSI(rPtr, 0)=BUSZSI(rPtr, N);
    DFN(rPtr)=N(rPtr);
    DFSIZ(rPtr)=SIZ(rPtr);
    DEFDAT(rPtr)=DAT(rPtr);
  }
  
  /**
   * randomPoly:
   * @Ptr: A C-Cube polynomial.
   * @p: A prime number.
   * 
   * The input 'Ptr' is an 'clean' C-Cube polynomial, i.e. 
   * 'Ptr' is initialized with zero coefficients and partial degrees.
   * This routine fills random numbers in Z/pZ into 'Ptr's 
   * coefficient vector.
   *
   * Return value: 
   **/
  void randomPoly(preFFTRep * Ptr, sfixn p){
	SERBPAS::randomPoly_inner(N(Ptr), BUSZS(Ptr), CUM(Ptr), DAT(Ptr), p);  
  }
  
  // Subroutine for the above function
  void randomPoly_inner(sfixn N, sfixn * dgs, sfixn * accum, sfixn * data, sfixn p){
	int i, offset=0;
	
	if(N==1){
	  srand(getSeed());
	  for(i=0; i<=dgs[N]; i++){
		data[i]=(rand()%p);
	  }
	  if (! data[dgs[N]]) data[dgs[N]]=1;
	  return;
	}
	
	
	for(i=0; i<=dgs[N]; i++){
	  offset=accum[N]*i;
	  SERBPAS::randomPoly_inner(N-1, dgs, accum, data+offset, p);
	} 	
  }

    /**
   * CopyOnePoly:
   * @rPtr: Destination polynomial.
   * @fPtr: Source polynomial.
   * 
   * Make a deep copy of 'fPtr' and save the copy in 'rPtr'.
   *
   * Return value: 'rPtr'.
   **/
  void CopyOnePoly(preFFTRep * rPtr, preFFTRep * fPtr ){  
	register int j;
    sfixn N=N(rPtr);
    N(rPtr)=N;
    SIZ(rPtr)=1;
    for(j=1;j<=N;j++){
      BUSZSI(rPtr, j)=BUSZSI(fPtr, j);    
      CUMI(rPtr, j)=CUMI(fPtr, j);
      SIZ(rPtr)=SIZ(rPtr)*(BUSZSI(rPtr, j)+1);
    } 
    OFST(rPtr)=0;
    for(j=0; j<SIZ(rPtr);j++) DATI(rPtr, j)=DATI(fPtr, j);
    //BUSZSI(rPtr, 0)=BUSZSI(rPtr, N);
    DFN(rPtr)=N(rPtr);
    DFSIZ(rPtr)=SIZ(rPtr);
    DEFDAT(rPtr)=DAT(rPtr);
  }
  
      /**
   * CopyOnePoly_sameData:
   * @rPtr: Destination polynomial.
   * @fPtr: Source polynomial.
   * 
   * Make a deep copy of 'fPtr' and save the copy in 'rPtr'.
   *
   * Return value: 'rPtr'.
   **/
  preFFTRep * CopyOnePoly_sameData(preFFTRep * fPtr ){  
	register int j;
	preFFTRep * rPtr = (preFFTRep *)my_calloc(1, sizeof(preFFTRep));

    sfixn N = N(fPtr);
    N(rPtr)=N;
    SIZ(rPtr)=SIZ(fPtr);

	BUSZS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUTS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
	CUM(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));

    for(j=1;j<=N;j++){
      BUSZSI(rPtr, j)=BUSZSI(fPtr, j);  
	  CUTSI(rPtr, j)=CUTSI(fPtr, j);
      CUMI(rPtr, j)=CUMI(fPtr, j);
      //SIZ(rPtr)=SIZ(rPtr)*(BUSZSI(rPtr, j)+1);
    } 
    OFST(rPtr) = OFST(fPtr);
	DAT(rPtr) = DAT(fPtr); //point to the same data array
    //for(j=0; j<SIZ(rPtr);j++) DATI(rPtr, j)=DATI(fPtr, j);
    //BUSZSI(rPtr, 0)=BUSZSI(rPtr, N);
    DFN(rPtr)=N(fPtr);
    DFSIZ(rPtr)=SIZ(fPtr);
    DEFDAT(rPtr)=DAT(fPtr);
	return rPtr;
  }
  
  /**
 * freePoly_not_data:
 * @x: A C-Cube polynomial.
 * 
 * Free a C-Cube polynomial.
 *
 * Return value: 
 **/
void freePoly_not_data(preFFTRep * x){

  if(BUSZS(x) !=NULL) {my_free(BUSZS(x)); BUSZS(x)=NULL;}
  if(CUTS(x) !=NULL) {my_free(CUTS(x)); CUTS(x)=NULL;}
  if(CUM(x) !=NULL) {my_free(CUM(x)); CUM(x)=NULL;}
  //if(DAT(x) !=NULL) {my_free(DAT(x)); DAT(x)=NULL;}

}

}//end of SERBPAS

