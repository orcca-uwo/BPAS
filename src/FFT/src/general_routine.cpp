// -*- C++ -*-
// general_routine.cilk
//@author Yuzhen Xie

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include<cmath>
#include "../../../include/FFT/src/modpn.h"
#include "../../../include/FFT/src/include/example_util_gettime.h"
#include "../../../include/FFT/src/general_routine.h"
#include <iostream>
#include <string.h>

//#include <reducer_opadd.h>

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef int sfixn;

//assume max(deg) <= 1000,000
//size of array factor is log2(1000000) <= 20
// factor[1] is the number of factors
#define NUMFACT 21

namespace PBPAS {

  //---------------------------------------------------
  void checkFrourierPrime(sfixn e, MONTP_OPT2_AS_GENE * pPtr){
	if(e > (pPtr->Npow)){
	  std::cout << "The given Fourier prime number can't handle this FFT multiplication in terms of the FFT-size." << pPtr->P << ", "  << e << std::endl;
	  std::cout << "Please choose another Fourier prime with bigger FFT-size." << std::endl;
	  exit (1);
	}
  }

  /**------------------------------------------------------------------
   * multi_mat_transpose for evaluation (forward FFT)
   * @N: number of dimensions of Matrix M. N>2
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
  void multi_mat_transpose_2DTran_eval(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data){

	//std::cout<<"dims[1], dims[2]: "<< dims[1] <<", "<< dims[2] <<std::endl;
	//for(int j = 0; j < n; j++) std::cout<<data[j]<<", ";
	//std::cout<<"----end of in data----" <<std::endl;

	if (dm==N) {
	  sfixn nj = dims[dm];
	  sfixn ni = n/nj;
	  sfixn lda = nj;
	  sfixn ldb = ni;

	  if (ni==nj){
		PBPAS::sqtranspose(data, lda, 0, ni, 0, 1, ni);
	  }else{
		//std::cout<<"Matteo's rect-transpose" <<std::endl;
		sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
		PBPAS::transpose(data, ldb, B, lda, 0, nj, 0, ni);
		memcpy(data, B, n*sizeof(sfixn));
		//cilk_for(int j = 0; j < n; j++) data[j] = B[j];
		my_free(B);
	  }
	}else{
	  sfixn ni = 1;
	  for (int i=1; i<dm; i++){
		ni *= dims[i];
	  }

	  sfixn nj = dims[dm];
	  sfixn lda = nj;
	  sfixn ldb = ni;
	  sfixn s1 = ni*nj;
	  //std::cout<<"n= "<<n<<std::endl;
	  //std::cout<<"s1= "<<s1<<std::endl;
	  sfixn s2 = n/s1;

	  if (ni==nj){
		cilk_for (int i=0; i<s2; i++) {
		  PBPAS::sqtranspose(data+s1*i, lda, 0, ni, 0, 1, ni);
		}
	  }else{
		size_t s = s1*sizeof(sfixn);
		cilk_for (int i=0; i<s2; i++) {
		  sfixn *B = (sfixn * ) my_calloc(s1, sizeof(sfixn));
		  PBPAS::transpose(data+s1*i, ldb, B, lda, 0, nj, 0, ni);
		  memcpy(data+s1*i, B, s);
		  //cilk_for(int j = 0; j < n; j++) data[j] = B[j];
		  my_free(B);
		}
	  }
	}
  }

  //---------------------------
  //data transposition in inverse tft
  void multi_mat_transpose_2DTran_interp(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data){

	//std::cout<<"dims[1], dims[2]: "<< dims[1] <<", "<< dims[2] <<std::endl;
	//for(int j = 0; j < n; j++) std::cout<<data[j]<<", ";
	//std::cout<<"----end of in data----" <<std::endl;

	if (dm==1) {
	  sfixn nj = dims[dm];
	  sfixn ni = n/nj;
	  sfixn lda = nj;
	  sfixn ldb = ni;

	  if (ni==nj){
		PBPAS::sqtranspose(data, lda, 0, ni, 0, 1, ni);
	  }else{
		//std::cout<<"Matteo's rect-transpose" <<std::endl;
		sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
		PBPAS::transpose(data, ldb, B, lda, 0, nj, 0, ni);
		memcpy(data, B, n*sizeof(sfixn));
		//cilk_for(int j = 0; j < n; j++) data[j] = B[j];
		my_free(B);
	  }
	}else{
	  sfixn ni = 1;
	  for (int i=dm+1; i<=N; i++){
		ni *= dims[i];
	  }
	  sfixn nj = dims[dm];
	  sfixn lda = nj;
	  sfixn ldb = ni;
	  sfixn s1 = ni*nj;
	  //std::cout<<"n= "<<n<<std::endl;
	  //std::cout<<"s1= "<<s1<<std::endl;
	  sfixn s2 = n/s1;

	  if (ni==nj){
		cilk_for (int i=0; i<s2; i++) {
		  PBPAS::sqtranspose(data+s1*i, lda, 0, ni, 0, 1, ni);
		}
	  }else{
		size_t s = s1*sizeof(sfixn);
		cilk_for (int i=0; i<s2; i++) {
		  sfixn *B = (sfixn * ) my_calloc(s1, sizeof(sfixn));
		  PBPAS::transpose(data+s1*i, ldb, B, lda, 0, nj, 0, ni);
		  memcpy(data+s1*i, B, s);
		  //cilk_for(int j = 0; j < n; j++) data[j] = B[j];
		  my_free(B);
		}
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
  void multi_mat_transpose(sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data){
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
		PBPAS::sqtranspose(data, lda, 0, ni, 0, 1, ni);
	  }else{
		//std::cout<<"Matteo's rect-transpose" <<std::endl;
		sfixn *B = (sfixn * ) my_calloc(n, sizeof(sfixn));
		PBPAS::transpose(data, ldb, B, lda, 0, nj, 0, ni);
		memcpy(data, B, n*sizeof(sfixn));
		//cilk_for(int j = 0; j < n; j++) data[j] = B[j];
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
	  PBPAS::decompose(N, dm, accum, dims, tmpPtr, data);
#ifdef TRACETRAN
	  long end = example_get_time();
	  decomposePar_time += (end - start) / 1.f;
	  std::cout << "decomposePar="<<decomposePar_time<< std::endl;
#endif
	  //for(i = 0; i < n; i++) data[i] = tmpPtr[i];
	  //cilk_for(int j = 0; j < n; j++) data[j] = tmpPtr[j];
	  memcpy(data, tmpPtr, n*sizeof(sfixn));
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
	  //cilk_for(int k=0; k<acN; k++) {
	  for(int k=0; k<acN; k++) {
		cilk_for(int i = 0; i < dims[1]; i++){
		  for(int j = 0; j < dims[dm]; j++){
			tmpPtr[j + i*ac + k*dims[N]] = data[j*accum[dm] + i + k*dims[1]];
		  }
		}
	  }
	  return;
	}

	sfixn acN = accum[N+1]/accum[N];
	cilk_for(int k=0; k<acN; k++) {
	  PBPAS::decompose(N-1, dm, accum, dims, tmpPtr+k*accum[N], data+k*accum[N]);
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

      siz = PBPAS::getDenseSiz(N-1, poly, BUSZSI(poly, N-1), data, CUMI(poly, N-1));
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
  // use the data field of result for 1st multiplier by TFT
  //thus initialize only one data array for 2nd multiplier
  //M=1, to call freeKroTFTRep
  void
  InitKroTFTRep_multiplier_1(KroTFTRep * kPtr, sfixn * resDgs, sfixn N, sfixn M,   MONTP_OPT2_AS_GENE * pPtr){
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
    //initialize only one data array for 2nd multiplier
	//use the data struct of result poly for 1st poly data
	//result is directly in the data struct of result poly :)
    DATS(kPtr)=(sfixn **)my_calloc(1, sizeof(sfixn *));
    //
	DATSI(kPtr, 0)=(sfixn * )my_calloc(SIZ(kPtr), sizeof(sfixn));

    KN(kPtr)=1; KE(kPtr)=0;
    for(j=1;j<=N;j++) {KN(kPtr)*=DIMSI(kPtr, j); KE(kPtr)+=ESI(kPtr, j);}
    KROOTS(kPtr)=NULL;
    kPtr->Defsize=SIZ(kPtr);
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
	//std::cout<<"enter freeKroTFTRep"<<std::endl;
    int i;
    if(KROOTS(x) != NULL) {my_free(KROOTS(x)); KROOTS(x)=NULL;}
    if(CUM(x) !=NULL) {my_free(CUM(x)); CUM(x)=NULL;}
    if(LS(x) !=NULL) {my_free(LS(x)); LS(x)=NULL;}
    if(ES(x) !=NULL) {my_free(ES(x)); ES(x)=NULL;}
    if(DIMS(x) !=NULL) {my_free(DIMS(x)); DIMS(x)=NULL;}
	//std::cout<<"before DATS"<<std::endl;
	//std::cout<<"M(x) "<<M(x) <<std::endl;
    if(DATS(x) != NULL){
      for(i=0; i<M(x); i++)
		if(DATSI(x, i) !=NULL) {my_free(DATSI(x, i)); DATSI(x, i)=NULL;}
      my_free(DATS(x));
      DATS(x)=NULL;}
	//std::cout<<"after DATS"<<std::endl;
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
   * Initialized a 'clean' (all zero-coefficients) C-Cube polynomial whose
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
      BUSZSI(rPtr, j) = p1dgs[j]+p2dgs[j];
      CUTSI (rPtr, j) = BUSZSI(rPtr, j);
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
	//std::cout<<"InitResPoly, SIZ="<<SIZ(rPtr)<<std::endl;
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
	//cilk_for(int i=0; i<=d; i++) { desV[i] = srcV[i]; }  //cilk_for slows down
	memcpy(desV, srcV, (d+1)*sizeof(sfixn));
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
	//cilk_for(int i=1; i<=n; i++) desV[i]=srcV[i];
	sfixn * desV0 = desV + 1;
	sfixn * srcV0 = srcV + 1;
	memcpy(desV0, srcV0, n*sizeof(sfixn));
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
	if (! PBPAS::compareVec(N(Ptr1), BUSZS(Ptr1), BUSZS(Ptr2))) return 0;
	if (! PBPAS::compareVec(SIZ(Ptr1)-1, DAT(Ptr1), DAT(Ptr2))) return 0;
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
	  cilk_for(int i=0; i<=d; i++){
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
	cilk_for(int i=0; i<=d; i++){
	  PBPAS::fromtofftRepMultiD(N-1, rccum, res+i*rccum[N], ccum, dgs, coeffs+i*ccum[N]);
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
	cilk_for(int i=0; i<=deg; i++){
	  vec1[i] = vec2[deg-i];
	}
	return vec1;
  }

  //==========================================================
  // clean vectors.
  //==========================================================
  // note: clean the whole vector.
  void cleanVec(sfixn deg, sfixn * cof){
	cilk_for(int i=0; i<=deg; i++) cof[i]=0;
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
// 	register sfixn i;
// 	sfixn end=deg*sizOfCoef;
// 	for(i=0; i<=end; i+=sizOfCoef){
// 	  for(int j=0; j<sizOfCoef; j++) outVec[end-i+j]=inVec[i+j];
// 	}

	sfixn end = deg*sizOfCoef;
	cilk_for(int k=0; k<=deg; k++) {
	  sfixn i = k*sizOfCoef;
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
	// register int i, j;
// 	sfixn * tmpCoef=(sfixn * )my_malloc(sizOfCoef*sizeof(sfixn));
// 	sfixn * tmp1Ptr=vec, * tmp2Ptr=vec+sizOfCoef*deg;

// 	for(i=0; i<((deg+1)/2); i++){
// 	  for(j=0; j<sizOfCoef; j++) tmpCoef[j]=tmp1Ptr[j];
// 	  for(j=0; j<sizOfCoef; j++) tmp1Ptr[j]=tmp2Ptr[j];
// 	  for(j=0; j<sizOfCoef; j++) tmp2Ptr[j]=tmpCoef[j];
// 	  tmp1Ptr+=sizOfCoef;
// 	  tmp2Ptr-=sizOfCoef;
// 	}
// 	my_free(tmpCoef);

//-----------------------------------------
//4535ms for ./normalform 2 1000 1000 2 on diplodocus
	sfixn * tmp1Ptr=vec;
	sfixn * tmp2Ptr=vec + sizOfCoef*deg;
	sfixn * tmpCoef=(sfixn * )my_malloc(sizOfCoef*sizeof(sfixn));

	for(int i=0; i<((deg+1)/2); i++){
	  cilk_for(int j=0; j<sizOfCoef; j++) tmpCoef[j]=tmp1Ptr[j];
	  cilk_for(int j=0; j<sizOfCoef; j++) tmp1Ptr[j]=tmp2Ptr[j];
	  cilk_for(int j=0; j<sizOfCoef; j++) tmp2Ptr[j]=tmpCoef[j];

	  tmp1Ptr+=sizOfCoef;
	  tmp2Ptr-=sizOfCoef;
	}
	my_free(tmpCoef);

//-----------------------------------------------------
// 	cilk_for(int i=0; i<((deg+1)/2); i++){
// 	  sfixn * tmp1Ptr=vec + i*sizOfCoef;
// 	  sfixn * tmp2Ptr=vec + sizOfCoef*deg - i*sizOfCoef;
// 	  sfixn * tmpCoef=(sfixn * )my_malloc(sizOfCoef*sizeof(sfixn));
// 	  //increase memory usage,
// 	  //cause memo swap for cilkscreen ./normalform 2 1000 1000 2
// 	  //4581ms for ./normalform 2 1000 1000 2

// 	  cilk_for(int j=0; j<sizOfCoef; j++) tmpCoef[j]=tmp1Ptr[j];
// 	  cilk_for(int j=0; j<sizOfCoef; j++) tmp1Ptr[j]=tmp2Ptr[j];
// 	  cilk_for(int j=0; j<sizOfCoef; j++) tmp2Ptr[j]=tmpCoef[j];

// 	  my_free(tmpCoef);
// 	}

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

	cilk_for(int i=0; i<size; i++){
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
  cilk_for(int i=0; i<M(kPtr); i++){
    cilk_for(int j=0; j<kPtr->Defsize; j++){
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
	  //cilk_for(int i=0; i<=d; i++){
	  //	res[i]=coeffs[i];
	  //}
	  memcpy(res, coeffs, (d+1)*sizeof(sfixn));
	  return;
	}

	if ( __cilkrts_get_nworkers() >=10 ) {
	  //#pragma cilk_grainsize = 10000;
	  for(int i=0; i<=d; i++){
		PBPAS::fromtofftRep(N-1, rccum, res+i*rccum[N], ccum, dgs, coeffs+i*ccum[N]);
	  }
	}else{
	  //#pragma cilk_grainsize = 1;
	  cilk_for(int i=0; i<=d; i++){
		PBPAS::fromtofftRep(N-1, rccum, res+i*rccum[N], ccum, dgs, coeffs+i*ccum[N]);
	  }
	}

  }

  void InitOneReducedPoly(preFFTRep * rPtr, sfixn N, sfixn * p1dgs){
	register int j;
	//std::cout<<"enter PBPAS::InitOneReducedPoly, N="<<N<<std::endl;
	//std::cout<<"enter PBPAS::InitOneReducedPoly, rPtr="<<rPtr<<std::endl;
	//std::cout<<"PBPAS::InitOneReducedPoly, p1dgs[N]="<<p1dgs[N]<<std::endl;
    N(rPtr)=N;
	//std::cout<<"PBPAS::InitOneReducedPoly, before BUSZS(rPtr)"<<std::endl;
    BUSZS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
	//std::cout<<"PBPAS::InitOneReducedPoly, after BUSZS(rPtr)"<<std::endl;
    CUTS(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    SIZ(rPtr)=1;
    CUM(rPtr)=(sfixn * )my_calloc(N+1,sizeof(sfixn));
    CUMI(rPtr, 1)= 1;
	//std::cout<<"PBPAS::InitOneReducedPoly, before for"<<std::endl;
    for(j=1;j<=N;j++){
	  //std::cout<<"PBPAS::InitOneReducedPoly, for j="<<j<<std::endl;
      BUSZSI(rPtr, j)=p1dgs[j];
      CUTSI (rPtr, j)=p1dgs[j];
      SIZ(rPtr)=SIZ(rPtr)*(BUSZSI(rPtr,j)+1);
      if(j>=2){
		CUMI(rPtr, j)=CUMI(rPtr, j-1)*(BUSZSI(rPtr,j-1)+1);
      }
    }
	//std::cout<<"InitOneReducedPoly, after for"<<std::endl;
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
	PBPAS::subPoly_inner_1(N, CUM(Ptr1), BUSZS(Ptr2), CUM(Ptr2), DAT(Ptr1), DAT(Ptr2), p);
  }


  // we will use the smaller dgs which are 2nd.
  // but accums should used prepectively.
  void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p)
  {
	//int i, offset1=0, offset2=0;
	if(N==1){
	  cilk_for(int i=0; i<=dgs2[1]; i++) data1[i] = SubMod(data1[i], data2[i], p);
	  return;
	}

	cilk_for(int i=0; i<=dgs2[N]; i++){
	  //offset1=accum1[N]*i;
	  //offset2=accum2[N]*i;
	  PBPAS::subPoly_inner_1(N-1, accum1, dgs2, accum2, data1+accum1[N]*i, data2+accum2[N]*i, p);
	}
  }

  //--------------------------------
  void subEqDgPoly_inner_1
  (sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int selector)
  {
	//int i, offset=0;
	if(N==1){
	  if(selector==1){
		cilk_for(int i=0; i<=dgs[1]; i++) data1[i]=SubMod(data1[i],data2[i],p);}
	  else{
		cilk_for(int i=0;i<=dgs[1];i++) data2[i]=SubMod(data1[i],data2[i],p);}
	  return;}

	cilk_for(int i=0; i<=dgs[N]; i++){
	  //offset=accum[N]*i;
	  subEqDgPoly_inner_1(N-1, dgs, accum, data1+accum[N]*i, data2+accum[N]*i, p, selector);
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

  //------------------------------------
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
	PBPAS::fromtofftRep(N-1,  CUM(co), DAT(co), CUM(f), BUSZS(f), DAT(f));
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
	  PBPAS::InitOneRevInv(i, ELEMI(tRevInvPtr, i), BDS(tPtr), dgs[i]);
	  (NLB(tRevInvPtr))[i]= BUSZSI(ELEMI(tRevInvPtr, i), i);
	  (NSB(tRevInvPtr))[i]=BDSI(tPtr, i)-1;
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
	PBPAS::initRandomTriSet(N, dgbound, tPtr, pPtr);
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

	for(i=1;i<=N;i++) {while( BDSI(tPtr,i)<=0 ) BDSI(tPtr, i)=rand()%dgbound;}
	for(i=1; i<=N; i++){
	  ELEMI(tPtr, i)=(preFFTRep *)my_calloc(1,sizeof(preFFTRep));
	  PBPAS::InitOneMonicRandomPolys(BDS(tPtr), i, ELEMI(tPtr, i), pPtr, i);
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
    PBPAS::randomMonicPoly(p1Ptr, pPtr->P);
    DFN(p1Ptr)=N(p1Ptr);
    DFSIZ(p1Ptr)=SIZ(p1Ptr);
    DEFDAT(p1Ptr)=DAT(p1Ptr);
    //BUSZSI(p1Ptr, 0)=BUSZSI(p1Ptr, N);
  }

  //---------------------------
  void randomMonicPoly(preFFTRep * Ptr, sfixn p){
	PBPAS::randomMonicPoly_inner(N(Ptr), N(Ptr), BUSZS(Ptr), CUM(Ptr), DAT(Ptr), p);
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
		PBPAS::randomMonicPoly_inner(N, N1-1, dgs, accum, data+offset, p);
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
	poly = PBPAS::EX_InitOnePoly(N, dgs);
	PBPAS::randomPoly(poly, p);
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
	PBPAS::InitOnePoly(poly, N, dgs);
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
	PBPAS::randomPoly_inner(N(Ptr), BUSZS(Ptr), CUM(Ptr), DAT(Ptr), p);
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
	  PBPAS::randomPoly_inner(N-1, dgs, accum, data+offset, p);
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

  //return index t in ls such that the difference of
  //ls[1]*...*ls[i] and ls[i+1]*...*ls[N] is minimum
  sfixn noExtensionBalancedBivar(sfixn N, sfixn* ls) {
	sfixn diff, difft, sl, sr, slt, srt, t, tt;

	if (N>1) {
	  slt=ls[1];
	  srt=1;

	  if (N==2) return 1;

	  for (int i=3; i<=N; i++) srt *=ls[i];

	  if (slt<=srt) {
		sl = slt*ls[2];
		sr = srt;
		t = 2;
	  }else{
		sr = srt*ls[2];
		sl = slt;
		t = 1;
	  }

	  diff = abs(sl-sr);

	  slt = ls[1]*ls[2];
	  for (int i=3; i<=N; i++) {
		srt = srt/ls[i];
		if (slt<=srt) {
		  sl = slt*ls[i];
		  sr = srt;
		  tt = i;
		}else{
		  sr = srt*ls[i];
		  sl = slt;
		  tt = i-1;
		}

		difft = abs(sl-sr);
		if (difft<diff) {
		  diff = difft;
		  t = tt;
		}
		slt *= ls[i];
	  }
	  return t;
	}else{
	  std::cout<<"Number of variables is less than 1"<<std::endl;
	  exit(1);
	}
  }

  /*-------------------------------------
   *assume max(deg) <= 1000,000
   *size of array factor is log2(1000000) <= 20
   *@ls[i]=partialdeg[i]+1
   *@factm[0]=m, index of variable being splitted
   *factm[1]=s1, factm[2]=s2 s.t. ls[m]=s1*s2
   *------------------------------------------*/
  void ReductionToBalancedBivar(sfixn N, sfixn* ls, sfixn *factm) {
	sfixn factor[NUMFACT];
	sfixn diff, difft, m, ml, mr, sl, sr, slt, srt, nf, mnf;

	if (N>=1) {
	  sl = 1;
	  sr = 1;
	  for (int i=2; i<=N; i++) sr *= ls[i];
	  slt = sl;
	  srt = sr;

	  prime_factorization(ls[1], factor);
	  nf = factor[0]; //number of prime factors
	  for (int j=nf; j>0; j--) {
		if (slt <= srt)
		  slt = slt*factor[j];
		else
		  srt = srt*factor[j];
	  }
	  m = 1;
	  ml = slt/sl;
	  mr = srt/sr;
	  diff = abs(slt-srt);

	  for (int i=2; i<=N; i++) {
		sl *= ls[i-1];
		sr /= ls[i];

		slt = sl;
		srt = sr;
		prime_factorization(ls[i], factor);
		nf = factor[0];
		for (int j=1; j<=nf; j++) {
		  if (slt <= srt)
			slt = slt*factor[j];
		  else
			srt = srt*factor[j];
		}
		difft = abs(slt-srt);
		if (difft < diff) {
		  diff = difft;
		  m = i;
		  ml = slt/sl;
		  mr = srt/sr;
		}
	  }
	  std::cout<<slt<<";"<<srt<<std::endl;
	  std::cout<<m<<";"<<ml<<";"<<mr<<std::endl;
	  factm[0]=m;
	  factm[1]=ml;
	  factm[2]=mr;
	}else{
	  std::cout<<"Number of variables is less than 1"<<std::endl;
	  exit(1);
	}
  }

  //size of fact is 20
  void prime_factorization(sfixn x, sfixn fact[]) {
	sfixn i, j, k;			/* counter */
	sfixn c;			/* remaining product to factor */

	c = x;
	j = 1;
	k = 0;

	while ((c % 2) == 0 && j<NUMFACT-1) {
	  //printf("%ld\n",2);
	  fact[j]=2;
	  j++;
	  k++;
	  c = c / 2;
	}

	i = 3;

	while (i <= (sqrt(c)+1) && j<NUMFACT-1) {
	  if ((c % i) == 0) {
		//printf("%ld\n",i);
		fact[j]=i;
		j++;
		k++;
		c = c / i;
	  }
	  else
		  i = i + 2;
	}

	if (c > 1) {
	  fact[j]=c;//printf("%ld\n",c);
	  k++;
	}

	fact[0] = k; //number of factors
	std::cout<<x<<std::endl;
	for(int m=0; m<NUMFACT; m++)
	  std::cout<<fact[m]<<";";

	std::cout<<std::endl;
  }

  //------------------------
  //compute the b for two univariate polynomial
  //with partial degrees d1 and d2
  //such that the degrees of bivariate representation
  //are sorted of balanced.
  sfixn V1to2V_pivot(sfixn d1, sfixn d2) {
	sfixn k = floor(sqrt((d1+d2)/2));

	sfixn b = k+1;
	sfixn r1 = d1 % b;
	sfixn r2 = d2 % b;
	sfixn q1 = (d1-r1)/b;
	sfixn q2 = (d2-r2)/b;
	//add abs after testing
	sfixn diffs1 = 2*b-q1-q2-2;
	if ( diffs1 == 0 ) {
	  //std::cout<<"balanced, k+1="<<b<<std::endl;
	  return b;
	}
	sfixn s1 = 2*(d1-r1+d2-r2)+diffs1+1;

	b = k;
	r1 = d1 % b;
	r2 = d2 % b;
	q1 = (d1-r1)/b;
	q2 = (d2-r2)/b;
	//add abs after testing
	sfixn diffs2 = 2*b-q1-q2-2;
	if ( diffs2 == 0 ) {
	  //std::cout<<"balanced, k="<<b<<std::endl;
	  return b;
	}

	//remove abs after testing
	if ( abs(diffs1) == abs(diffs2) ) {
	  sfixn s2 = 2*(d1-r1+d2-r2)+diffs2+1;
	  if ( s1 <= s2 ) {
		//std::cout<<"s1<=s2, k+1="<<s1<<";"<<s2<<";"<<k+1<<std::endl;
		return k+1;
	  }
	  else {
		//std::cout<<"s1>s2, k+1="<<s1<<";"<<s2<<";"<<k<<std::endl;
		return k;
	  }
	}
	else {
	  if ( abs(diffs1) < abs(diffs2) ) {
		//std::cout<<"diffs1<diffs2, k+1="<<diffs1<<";"<<diffs2<<";"<<k+1<<std::endl;
		return k+1;
	  }else {
		//std::cout<<"diffs1>diffs2, k="<<diffs1<<";"<<diffs2<<";"<<k<<std::endl;
		return k;
	  }
	}
  }

  //---------------------------------------------------------
  //deg1 and deg2 are partial degrees of poly f and g w.r.t 1st var
  //assumption they are far more than the other degrees
  sfixn multi_V1to2V_pivot_opt(sfixn sigma, sfixn deg1, sfixn deg2) {

	if (deg1+deg2 <= sigma) {
	  std::cout << "Degree of var 1 is too small" <<std::endl;
	  exit(1);
	}

	sfixn k = floor(sqrt( (sigma+1)*(sigma+1) + 8*(deg1+deg2)*sigma ));

	sfixn bt[27];
	sfixn size[27];
	sfixn diffDim[27];
	sfixn b, s1, s2, bf, pos1, minpos, posDiff0=-1;

	for (int i=-102; i<=5; i++) {
	  bf = sigma+k+i;
	  if ( (sigma+k+i) >=8 && (bf%4) == 0 ) {
		pos1 = i;
		b = bf/4;
		s1 = 2*b-1;
		s2 = ((deg1 - deg1 % b)/b + (deg2 - deg2 % b)/b + 1)*sigma;
		size[0] = s1*s2;
		minpos = 0;
		bt[0]=b;
		diffDim[0]=s1-s2;
		if ( diffDim[0]==0 )
		  posDiff0 = 0;

		break;
	  }
	}
	sfixn pos = 0;
	for (int i=pos1+1; i<=5; i++) {
	  bf = sigma+k+i;
	  if ( (bf%4) == 0 ) {
		b = bf/4;
		s1 = 2*b-1;
		s2 = ((deg1 - deg1 % b)/b + (deg2 - deg2 % b)/b + 1)*sigma;

		pos++;
		bt[pos]=b;
		diffDim[pos]=s1-s2;
		size[pos]=s1*s2;

		if (diffDim[pos]==0)
		  posDiff0 = pos;

		if (size[pos]<size[minpos])
		  minpos=pos;

	  }
	}
	if ( posDiff0 >=0 ) {
	  if( (float)size[posDiff0]/size[minpos] <= 1.05 )
		return bt[posDiff0];
	  else
		return bt[minpos];
	}else{
	  return bt[minpos];
	}
  }

  //----------------------------------------------------------
  //degs1[1] to degs1[N] partial degrees of poly 1
  //assumption degs1[1] >> degs1[i]
  sfixn multi_V1to2V_pivot(sfixn N, sfixn * degs1, sfixn *degs2) {
	sfixn sigma = 1;
	for (int i=2; i<=N; i++)
	  sigma *= degs1[i]+degs2[i]+1;

	std::cout<<"multi_V1to2V_pivot, sigma= "<<sigma<<std::endl;
	std::cout<<"multi_V1to2V_pivot, sigma quo 2= "<<(sigma-sigma%2)/2<<std::endl;

	sfixn k = floor(sqrt( (sigma+1)*(sigma+1) + 8*(degs1[1]+degs2[1])*sigma ));
	std::cout<<"multi_V1to2V_pivot, deltsquare= "<<(sigma+1)*(sigma+1) + 8*(degs1[1]+degs2[1])*sigma<<std::endl;
	std::cout<<"multi_V1to2V_pivot, k= "<<k<<std::endl;
	sfixn bf, b, bt, s1, s2, s, st, diff_s1s2, diff_s1s2t, pos1, q1, q2;
	for (int i=-102; i<=5; i++) {
	  bf = sigma+k+i;

	  std::cout<<"multi_V1to2V_pivot, bf= "<<bf<<std::endl;
	  if ( (bf%4) == 0 ) {
		pos1 = i;
		b = bf/4;
		std::cout<<"multi_V1to2V_pivot, first b= "<<b<<std::endl;
		s1 = 2*b-1;
		q1 = (degs1[1] - degs1[1] % b)/b;
		q2 = (degs2[1] - degs2[1] % b)/b;
		s2 = (q1 + q2 + 1)*sigma;

		std::cout<<"multi_V1to2V_pivot, s1, s2= "<<s1<<", "<<s2<<std::endl;
		diff_s1s2 = s1 - s2;
		s = s1*s2;
		std::cout<<"multi_V1to2V_pivot, pos1, diff_s1s2, s = "<<pos1<<", "<<diff_s1s2<<", "<<s<<std::endl;
		break;
	  }
	}
	for (int i=pos1+1; i<=5; i++) {
	  bf = sigma+k+i;
	  if ( bf%4 == 0 ) {
		bt = bf/4;
		std::cout<<"multi_V1to2V_pivot, bt= "<<bt<<std::endl;
		s1 = 2*bt-1;
		s2 = ((degs1[1]-degs1[1]%bt)/bt + (degs2[1]-degs2[1]%bt)/bt + 1)*sigma;
		diff_s1s2t = s1 - s2;
		st = s1*s2;
		std::cout<<"multi_V1to2V_pivot, i, diff_s1s2t, st = "<<i<<", "<<diff_s1s2t<<", "<<st<<std::endl;
		if ( diff_s1s2t <= diff_s1s2 ){
		  if ( diff_s1s2t < diff_s1s2 ) {
			b = bt;
			diff_s1s2 = diff_s1s2t;
			s = st;
		  }else{
			if (st < s) {
			  b = bt;
			  diff_s1s2 = diff_s1s2t;
			  s = st;
			}
		  }
		}
	  }
	}
	std::cout<<"multi_V1to2V_pivot, b= "<<b<<std::endl;
	return b;

  }

  //-----------------------
  //extend univar to bivar by pivot b
  //resDgs1=delt1 = 2b-2, resDgs2=delt2=d1 quo b + d2 quo b
  //
  void InitKroTFTRep_1V2V(KroTFTRep * kPtr, sfixn b, sfixn d1, sfixn d2){
    register int j;
    N(kPtr)=2;
    M(kPtr)=2;
    ES(kPtr)=(sfixn * )my_calloc(3,sizeof(sfixn));
    //
    DIMS(kPtr)=(sfixn * )my_calloc(3,sizeof(sfixn));
    //
    LS(kPtr)=(sfixn *)my_calloc(3,sizeof(sfixn));

    CUM(kPtr)=(sfixn * )my_calloc(3,sizeof(sfixn));
    //
    CUMI(kPtr, 1)= 1; SIZ(kPtr)=1;

	sfixn resDgs1 = 2*b-2;
	sfixn resDgs2 = (d1-d1%b)/b + (d2-d2%b)/b;
	std::cout<<resDgs1<<" "<<resDgs2<<" ";
	//std::cout<<"resDgs1 ="<<resDgs1<<std::endl;
	//std::cout<<"resDgs2 ="<<resDgs2<<std::endl;

	ESI(kPtr, 1)=logceiling(resDgs1+1);
	ESI(kPtr, 2)=logceiling(resDgs2+1);

	LSI(kPtr, 1)=resDgs1+1;
	LSI(kPtr, 2)=resDgs2+1;

    for(j=1;j<=2;j++){
      DIMSI(kPtr, j)=1<<(ESI(kPtr, j));
      SIZ(kPtr)*=LSI(kPtr, j);
      if(j>=2){CUMI(kPtr, j)=CUMI(kPtr, j-1)*LSI(kPtr, j-1);}
    }

    DATS(kPtr)=(sfixn **)my_calloc(2, sizeof(sfixn *));
    //
    for(j=0; j<2; j++){
      DATSI(kPtr, j)=(sfixn * )my_calloc(SIZ(kPtr),sizeof(sfixn));
      //
    }

    KN(kPtr)=1; KE(kPtr)=0;
    for(j=1;j<=2;j++) {KN(kPtr)*=DIMSI(kPtr, j); KE(kPtr)+=ESI(kPtr, j);}
    KROOTS(kPtr)=NULL;
    kPtr->Defsize=SIZ(kPtr);
  }

  //-------------------------
  // sigma = (d2f+d2g+1)*...*(dnf+dng+1)
  //q1f = deg(f,v1) quo b;
  void InitKroTFTRep_1V2V_multi(KroTFTRep *kPtr, sfixn b, sfixn sigma, sfixn q1f, sfixn q1g ) {
	register int j;
	N(kPtr)=2;
    M(kPtr)=2;
    ES(kPtr)=(sfixn * )my_calloc(3,sizeof(sfixn));
    //
    DIMS(kPtr)=(sfixn * )my_calloc(3,sizeof(sfixn));
    //
    LS(kPtr)=(sfixn *)my_calloc(3,sizeof(sfixn));

    CUM(kPtr)=(sfixn * )my_calloc(3,sizeof(sfixn));
    //
    CUMI(kPtr, 1)= 1; SIZ(kPtr)=1;

	sfixn resDgs1 = 2*b-2;
	sfixn resDgs2 = (q1f + q1g + 1)*sigma - 1;
	std::cout<<resDgs1<<" "<<resDgs2<<" ";

	ESI(kPtr, 1)=logceiling(resDgs1+1);
	ESI(kPtr, 2)=logceiling(resDgs2+1);

	LSI(kPtr, 1)=resDgs1+1;
	LSI(kPtr, 2)=resDgs2+1;

    for(j=1;j<=2;j++){
      DIMSI(kPtr, j)=1<<(ESI(kPtr, j));
      SIZ(kPtr)*=LSI(kPtr, j);
      if(j>=2){CUMI(kPtr, j)=CUMI(kPtr, j-1)*LSI(kPtr, j-1);}
    }

    DATS(kPtr)=(sfixn **)my_calloc(2, sizeof(sfixn *));
    //
    for(j=0; j<2; j++){
      DATSI(kPtr, j)=(sfixn * )my_calloc(SIZ(kPtr),sizeof(sfixn));
      //
    }

    KN(kPtr)=1; KE(kPtr)=0;
    for(j=1;j<=2;j++) {KN(kPtr)*=DIMSI(kPtr, j); KE(kPtr)+=ESI(kPtr, j);}
    KROOTS(kPtr)=NULL;
    kPtr->Defsize=SIZ(kPtr);

  }

  //--------------------------------------------------
  // extend 1V to 2V and store data to KroTFTRep data
  // for multiplication by TFT
  // d--degree of poly coeffs
  void from1Vto2VTFTRep(sfixn b, sfixn d, sfixn * coeffs, sfixn rccum2, sfixn * data){
	// b is pivot, also the size f dimension 1
	// b-1 is the degree of 1st var in 2V representation
	// d2 is the degree of 2nd var in 2V representation
	sfixn d2 = (d-d%b)/b;

	//std::cout<<"from1Vto2VtoTFTRep, d2="<<d2<<std::endl;

	cilk_for (int i=0; i<d2; i++) {
	  //sfixn pos1 = i*b;
	  //sfixn pos2 = i*rccum2;
	  for (int j=0; j<b; j++)
		data[i*rccum2+j] = coeffs[i*b+j];
	}

	//when i=d2, d2*b+b could be larger than d
	sfixn pos3 = d2*b;
	sfixn pos4 = d2*rccum2;
	for (int j=0; j<b; j++) {
	  if ( (pos3 + j) <= d )
		data[pos4+j] = coeffs[pos3+j];
	}

	//for (int i=0; i<=d; i++)
	//std::cout<<coeffs[i]<<";";

	//std::cout<<"----f coeffs----"<<std::endl;

	//for (int i=0; i<s; i++)
	//std::cout<<data[i]<<";";

	//std::cout<<"----res f coeffs----"<<std::endl;

  }

  //------------------------------
  // q1 = pdegs[1] quo b;
  // num_boxes_B = sigma = (d2+d2'+1)*...(dn+dn'+1);
  // num_blocks_B = (q1+q1'+1)
  void fromMulti_1Vto2VTFTRep_1(sfixn N, sfixn b, sfixn q1, sfixn *pdegs, sfixn *M1, sfixn num_boxes_B, sfixn num_blocks_B, sfixn *B) {
	sfixn box_size_M1 = pdegs[1]+1;
	sfixn block_size_M1 = b;
	sfixn num_blocks_M1 = q1+1;

	sfixn num_boxes_M1 = 1;
	for (int i=2; i<=N; i++)
	  num_boxes_M1 *= pdegs[i]+1;

	sfixn block_size_B = 2*b-1;
	sfixn box_size_B = block_size_B * num_blocks_B;
	std::cout<<"copy forward, num_boxes_M1= "<<num_boxes_M1<<std::endl;
	std::cout<<"copy forward, num_blocks_M1= "<<num_blocks_M1<<std::endl;
	std::cout<<"copy forward, box_size_B= "<<box_size_B<<std::endl;
	std::cout<<"copy forward, box_size_M1= "<<box_size_M1<<std::endl;

	cilk_for (int box=0; box<num_boxes_M1; box++) {
	  for (int block=0; block<q1; block++) {
		for (int slot=0; slot<b; slot++){
		  B[ box*box_size_B + block*block_size_B + slot ] =
			M1[ box*box_size_M1 + block*block_size_M1 + slot ];
		}
	  }

	  for (int slot=0; slot <=pdegs[1]%b; slot++)
		B[ box*box_size_B + q1*block_size_B + slot ] =
			M1[ box*box_size_M1 + q1*block_size_M1 + slot ];
	}

	for (int i=0; i<num_boxes_M1*box_size_M1; i++)
	  std::cout << M1[i]<<" ";
	std::cout <<"---M1---" <<std::endl;

	for (int i=0; i<box_size_B*num_boxes_B; i++)
	  std::cout << B[i]<<" ";
	std::cout <<"---copy M1 to B---" <<std::endl;
  }

  //-----------------------------------------------
  /**
   * @n: the number of variables.
   * @degs_M1: partial degree vector of poly f
   * @degs_M2: partial degree vector of poly g
   * size of degs_M1 and degs_M2 are n+1.
   * degs_M1[0] and degs_M2[0] are not used.
   * the other slots store the partial degrees.
   * copy coefficients in M1 to B while extending x1 to x1 and y
   * s.t. deg(f,x1)=degs_M1 rem b, deg(f,y)=degs_M1 quo b,
   * contract y and the rest variables into y
   * @q1f: = degs_M1[1] quo b
   * @q1g: = degs_M2[1] quo b
   *
   *#ifdef LINUXINTEL64
   *typedef int sfixn;
   **/
 void fromMulti_1Vto2VTFTRep(sfixn n, sfixn b, sfixn q1f, sfixn q1g, sfixn *degs_M1, sfixn *M1, sfixn *degs_M2, sfixn *B) {

   //sfixn degs_B[n];
   //degs_B[0]=0;
   sfixn *degs_B = (sfixn * ) my_calloc(n, sizeof(sfixn));
   degs_B[1]=(2*b-1)*(q1f +q1g +1)-1;

   //sfixn psize_B[n+1];
   //psize_B[0] = 0;
   sfixn * psize_B = (sfixn * ) my_calloc(n+1, sizeof(sfixn));
   psize_B[1] = degs_B[1] +1;
   sfixn xn_size_B = psize_B[1];
   for (int i=2; i<n; i++) {
	 degs_B[i] = degs_M1[i] +degs_M2[i];
	 psize_B[i] = degs_B[i] +1;
	 xn_size_B *= psize_B[i];
   }
   psize_B[n] = degs_M1[n] +degs_M2[n] +1;

   //sfixn psize_M1[n+1];
   //psize_M1[0] = 0;
   sfixn * psize_M1 = (sfixn * ) my_calloc(n+1, sizeof(sfixn));
   sfixn xn_size_M1 = 1;
   for (int i=1; i<n; i++) {
	 psize_M1[i] = degs_M1[i] +1;
	 xn_size_M1 *= psize_M1[i];
   }
   psize_M1[n] = degs_M1[n] +1;

   sfixn block_size_M1 = b;
   sfixn block_size_B = 2*b-1;
   sfixn r1 = degs_M1[1] % b;//=============

// 	std::cout <<"---xn_size_M1= " <<xn_size_M1<<std::endl;
// 	std::cout <<"---xn_size_B= " <<xn_size_B<<std::endl;
// 	for (int i=0; i<xn_size_B*(degs_M1[n]+degs_M2[n]+1); i++)
// 	  std::cout << B[i]<<" ";
// 	std::cout <<"---in B---" <<std::endl;

// 	for (int i=0; i<xn_size_M1*(degs_M1[n]+1); i++)
// 	  std::cout << M1[i]<<" ";
// 	std::cout <<"---M1---" <<std::endl;

	//can not use cilk_for, nextPosition, internal compiler error:
	//in gimplify_addr_expr, at gimplify.c:3904
	cilk_for (int en=0; en<=degs_M1[n]; en++) {
	  //sfixn p_M1[n+1];
	  //sfixn p_B[n+1];
	  sfixn *p_M1 = (sfixn * ) my_calloc(n+1, sizeof(sfixn));
	  sfixn *p_B = (sfixn * ) my_calloc(n+1, sizeof(sfixn));
	  p_M1[0] = en * xn_size_M1;
	  p_B[0] = en * xn_size_B;
	  p_M1[n]=en;
	  p_B[n]=en;
	  //sfixn offset_en_M1 = p_M1[0] + xn_size_M1;
	  //sfixn offset_en_B = p_B[0] + xn_size_B;

	  for (int i=1; i<n; i++) {
		p_M1[i]=0;
		p_B[i]=0;
	  }

	  int notdone = 0; // true in C
	  while (notdone==0) {
		for (int block=0; block<q1f; block++) {
		  for (int slot=0; slot<b; slot++) {
			B[ p_B[0] + block*block_size_B +slot ] = M1[ p_M1[0] + block*block_size_M1 +slot ];
		  }
		}

		for (int slot=0; slot<=r1; slot++) {
		  B[ p_B[0] + q1f*block_size_B +slot ] = M1[ p_M1[0] + q1f*block_size_M1 +slot ];
		}

// 		for (int i=0; i<xn_size_B*(degs_M1[n]+degs_M2[n]+1); i++)
// 		  std::cout << B[i] <<" ";
// 		std::cout <<"---B, en, p_B[0], p_M1[0]=---" <<en<<", "<<p_B[0]<<", "<<p_M1[0]<<std::endl;

		nextPosition(p_M1, degs_M1, n, psize_M1);
		if (p_M1[0]==-1)
		  notdone = -1; //false
		else {
		  nextPosition(p_B, degs_M1, n, psize_B);
// 		  for (int i=0; i<=n; i++)
// 			std::cout<<p_M1[i]<<" ";
// 		  std::cout <<"---p_M1, en=---" <<en<<std::endl;
// 		  for (int i=0; i<=n; i++)
// 			std::cout<<p_B[i]<<" ";
// 		  std::cout <<"---p_B, en=---" <<en<<std::endl;
		}
	  }
	  my_free(p_M1);
	  my_free(p_B);
	}
	my_free(psize_M1);
	my_free(psize_B);
  }

  //---------------------
  void nextPosition(sfixn *exp, sfixn *deg, sfixn n, sfixn *size) {
	//deg [0,0,d2,...,d_{n-1}]
	//exp [P,0,e2,...,en]
	//size [0,s1,...,sn]
	int found = -1;
	int notdone = 0; //true
	int i=2;
	while (notdone==0 && i<n) {
	  //for (int i=2; i<n; i++) {
	  if (exp[i]<deg[i]) {
		exp[i]++;
		if (i==2)
		  exp[0] += size[1];
		else{
		  for (int k=2; k<i; k++) //k<=i-1
			exp[k]=0;

		  sfixn acc = 1;
		  exp[0] = 0;
		  for (int k=1; k<i; k++)
			acc *= size[k];

		  for (int k=i; k<=n; k++) {
			exp[0] = exp[0] + acc*exp[k];
			acc *= size[k];
		  }
		  //std::cout<<"exp[0], i != 2: "<<exp[0]<<std::endl;
		}
		found=1;
		notdone = -1; //done
		//break;
	  }
	  i++;
	}
	if (found==-1) exp[0]=-1;

  }

  //-------------------
  //extend 2V *data to 1V and store data in *coeffs
  //after multiplication by 2D TFT
  //b is pivot
  //d is d1+d2, degree of product
  //s1 is the degree+1 of 1st var of the product in 2V representation
  //s2 is the degree+1 of 2nd var of the product in 2V representation
  void from2VTFTRepto1V_2Traversal(sfixn b, sfixn s1, sfixn s2, sfixn * data, sfixn siz, sfixn * coeffs, MONTP_OPT2_AS_GENE *pPtr){
	//	std::cout<<"from2VTFTRepto1V, siz="<<siz <<std::endl;

	sfixn delt2 = s2-1;
	sfixn delt1 = s1-1;
	sfixn s = s1*s2;
	//std::cout<<"delt2*b="<<delt2*b<<std::endl;
	//std::cout<<"delt2*s1="<<delt2*s1<<std::endl;

	//-----------------------------------
	// for (int v=0; v<delt2; v++) {
// 	  for (int u=0; u<=delt1; u++) {
// 		coeffs[v*b+u] = AddMod(coeffs[v*b+u], data[v*s1+u], pPtr->P);
// 	  }
// 	}
	//----------------------------------
	sfixn m = floor((delt2-1)/2);

	cilk_for (int v=0; v<=m; v++) {
	  for (int u=0; u<=delt1; u++) {
		coeffs[2*v*b+u] = data[2*v*s1+u];
	  }
	}

	sfixn r = (delt2-1) % 2;
	if ( r == 0 ) {
	  cilk_for (int v=0; v<m; v++) {
		for (int u=0; u<=delt1; u++) {
		  coeffs[(2*v+1)*b+u] = AddMod(coeffs[(2*v+1)*b+u], data[(2*v+1)*s1+u], pPtr->P);
		}
	  }
	}else{
	  cilk_for (int v=0; v<=m; v++) {
		for (int u=0; u<=delt1; u++) {
		  coeffs[(2*v+1)*b+u] = AddMod(coeffs[(2*v+1)*b+u], data[(2*v+1)*s1+u], pPtr->P);
		}
	  }
	}

	//for (int i=0; i<siz; i++)
	//  std::cout<<coeffs[i]<<";";
	//std::cout<<"----first coeffs----"<<std::endl;

	sfixn pos = delt2*b;
	sfixn pos2 = delt2*s1;
	//std::cout<<"delt2*b="<<pos<<std::endl;
	//std::cout<<"delt2*s1="<<delt2*s1<<std::endl;
	for (int u=0; u<=delt1; u++) {
	  //std::cout<<"u="<<u<<std::endl;
	  //sfixn pos = pos0+u;
	  //sfixn pos2 = pos02+u;
	  //std::cout<<"before if pos="<<pos<<std::endl;
	  if ( (pos+u<siz) && (pos2+u<s) ) {
		//std::cout<<"pos="<<pos<<std::endl;
		//std::cout<<"delt2*s1+u="<<delt2*s1+u<<std::endl;
		coeffs[pos+u] = AddMod(coeffs[pos+u], data[pos2+u], pPtr->P);
	  }
	}

	//for (int i=0; i<s; i++)
	//  std::cout<<data[i]<<";";

	//std::cout<<"----data----"<<std::endl;

    //for (int i=0; i<siz; i++)
	//  std::cout<<coeffs[i]<<";";

	//std::cout<<"----coeffs----"<<std::endl;
  }

  //-------------------
  //extend 2V *data to 1V and store data in *coeffs
  //after multiplication by 2D TFT
  //b is pivot
  //d is d1+d2, degree of product
  //s1 is the degree+1 of 1st var of the product in 2V representation
  //s2 is the degree+1 of 2nd var of the product in 2V representation
  void from2VTFTRepto1V(sfixn b, sfixn s1, sfixn s2, sfixn * data, sfixn siz, sfixn * coeffs, MONTP_OPT2_AS_GENE *pPtr){
	//	std::cout<<"from2VTFTRepto1V, siz="<<siz <<std::endl;

	//sfixn r = siz % b; //???????
	sfixn q = s2-1;
	sfixn r = siz-q*b-b;
	//std::cout<<"b, r, q " <<b<<", "<<r<<", "<<q<<std::endl;
	//for (int u=0; u<=b-1; u++) coeffs[u] = data[u];
	memcpy(coeffs, data, b*sizeof(sfixn));
	//std::cout<<"before 2nd for"<<std::endl;

	cilk_for(int w=1; w<=q; w++) {
	  sfixn X = w*b;
	  sfixn Y =(w-1)*(2*b-1)+b;
	  sfixn Z = w*(2*b-1);

	  for (int u=0; u<=b-2; u++) {
		coeffs[X + u] = AddMod(data[Y + u], data[Z + u], pPtr->P);
	  }
	  coeffs[X + (b-1)] = data[Z + (b-1)];
	}

	//std::cout<<"before last for"<<std::endl;
// 	for (int u=0; u<=r; u++) {
// 		coeffs[q*b+b+u] = data[q*(2*b-1)+b + u];
// 	}
	if ( r > 0 ) {
	  sfixn *coeffs_tail = coeffs + q*b+b;
	  sfixn *data_tail = data + q*(2*b-1)+b;
	  //std::cout<<"from2VTFTRepto1V, before memcpy tail"<<std::endl;
	  memcpy(coeffs_tail, data_tail, r*sizeof(sfixn));
	  //std::cout<<"from2VTFTRepto1V, after memcpy tail"<<std::endl;
	  //std::cout<<"after last for"<<std::endl;
	}
	//for (int i=0; i<s1*s2; i++)
	//  std::cout<<data[i]<<";";
	//std::cout<<"----data----"<<std::endl;
	//for (int i=0; i<siz; i++)
	//  std::cout<<coeffs[i]<<";";
	//std::cout<<"----coeffs----"<<std::endl;
  }

  //---------------------------------------
  //num_block_B = q1 + q1'+1
  //box_size_M = d1 + d1'+1
  void from2VTFTReptoMultiV(sfixn b, sfixn sigma, sfixn num_block_B, sfixn *B, sfixn box_size_M, sfixn *M, MONTP_OPT2_AS_GENE *pPtr){

	sfixn num_box_M = sigma;
	sfixn block_size_B = 2*b -1;
	sfixn box_size_B = block_size_B * num_block_B;
	sfixn num_box_B = sigma;

	sfixn q = num_block_B - 1;

// 	std::cout<<"copy back, num_box_B= "<<num_box_B<<std::endl;
// 	std::cout<<"copy back, box_size_B= "<<box_size_B<<std::endl;
// 	std::cout<<"copy back, box_size_M= "<<box_size_M<<std::endl;

// 	for (int i=0; i<block_size_B*num_block_B*sigma; i++)
// 	  std::cout<<B[i]<<";";
// 	std::cout<<"---2V result---"<<std::endl;


	cilk_for (int box=0; box<num_box_B; box++) {
	  sfixn p_B = box * box_size_B;
	  sfixn p_M = box * box_size_M;

	  for (int u=0; u<b; u++) {
		M[p_M+u] = B[p_B+u];
	  }
	  for (int w=1; w<=q; w++) {
		for (int u=0; u<=b-2; u++) {
		  M[p_M+w*b+u] = AddMod(B[p_B+(w-1)*(2*b-1)+b+u], B[p_B+w*(2*b-1)+u], pPtr->P);
		}
		M[p_M+w*b + b-1] = B[p_B+w*(2*b-1) + b-1];
	  }

	  sfixn r =box_size_M-1 -q*b -b;
	  for (int u=0; u<=r; u++)
		M[p_M + q*b+b+u] = B[p_B +q*(2*b-1)+b + u];
	}

  }

}//end of PBPAS

