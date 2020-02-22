// -*- C++ -*- 
// general_routine.cilk
//@author Yuzhen Xie

#include "../../../include/FFT/src/general_routine.h"
#include "malloc.h"
namespace PBPAS {

  /**-----------------------------------------------------------
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
  void EX_Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, 
					sfixn * rootsPtr, 
					MONTP_OPT2_AS_GENE * pPtr)
  {//---------------------------------------------------------------
    if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){
      Mont_GetNthRoots_OPT2_AS_GENE(e, n, rootsPtr, pPtr);
    }
    else{
      Mont_GetNthRoots_OPT2_AS_GENE_SPE(e, n, rootsPtr, pPtr); 
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
  void Mont_GetNthRoots_OPT2_AS_GENE(sfixn e, sfixn n, sfixn * rootsPtr, 
				     MONTP_OPT2_AS_GENE * pPtr)
  {//--------------------------------------------------------------------
   //std::cout<<"Mont_GetNthRoots_OPT2_AS_GENE\n"; 
   register int j;
    sfixn root, rootR, R=(1L<<pPtr->Rpow)%pPtr->P, R_2=MulMod(R, R, pPtr->P), BRsft=pPtr->Base_Rpow;
    //std::cout<<"R_2: "<<R_2<<"\n";
    //std::cout<<"BRsft: "<<BRsft<<"\n";
    // printf("The input FFT-N:%ld, this p:%ld can handle at most: %ld. !!!\n\n\n\n",n,pPtr->P, (1L<<pPtr->Npow));
    
    //root=EX_GetPrimitiveNthRoot(e, n, pPtr->P);
    root=PowerMod(pPtr->Max_Root, 1L<<((pPtr->Npow)-e), pPtr->P);
    //std::cout<<"root: "<<root<<"\n";
    rootR=MontMulMod_OPT2_AS_GENE(root, R_2<<BRsft, pPtr);
    rootsPtr[0]=R<<BRsft;
    rootsPtr[1]=rootR;
    //std::cout<<"rootR: "<<rootR<<"\n";
    for(j=2; j<n; j++) { //data dependence 
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
   * @pPtr: prime number structure for the prime number p
   * 
   *  Assume p is p=N-1, where N is a power of 2.
   *
   * Return value: returns `rootsPtr`
   **/
  void Mont_GetNthRoots_OPT2_AS_GENE_SPE(sfixn e, sfixn n, sfixn * rootsPtr, 
					 MONTP_OPT2_AS_GENE * pPtr)
  {//-------------------------------------------------------------------------
    //std::cout<<"Mont_GetNthRoots_OPT2_AS_GENE_SPE\n";
    register int j;
    sfixn root, rootR, R=(1L<<pPtr->Rpow)%pPtr->P, R_2=MulMod(R, R, pPtr->P), BRsft=pPtr->Base_Rpow;
    //std::cout<<"R_2: "<<R_2<<"\n";
    //std::cout<<"BRsft: "<<BRsft<<"\n";
    root=PowerMod(pPtr->Max_Root, 1L<<((pPtr->Npow)-e), pPtr->P);
    //std::cout<<"root: "<<root<<"\n";
    rootR=MontMulMod_OPT2_AS_GENE_SPE(root, R_2<<BRsft, pPtr);
    //std::cout<<"rootR: "<<rootR<<"\n";
    rootsPtr[0]=R<<BRsft;
    rootsPtr[1]=rootR;

    //cilk_for (int j=2; j<n; j++) { //potential of data race, use reducer??
    for (j=2; j<n; j++) {
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
  void EX_Mont_PairwiseMul_OPT2_AS(sfixn n, sfixn * APtr, sfixn * BPtr, sfixn p)
  {//------------------------------------------------------------------------
    //register int i;
    //#pragma cilk_grainsize = 2*n/__cilkrts_get_nworkers(); //8192; //need tune
    cilk_for (sfixn i=0; i<n; i++) 
      APtr[i] = MulMod(APtr[i], BPtr[i], p); //?why not MontMulMod_OPT2_AS_GENE
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
  void EX_Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, 
				     MONTP_OPT2_AS_GENE * pPtr)
  {//---------------------------------------------------------------------
     
    if((pPtr->c_pow==0) || ((pPtr->c_pow+pPtr->Rpow)>=BASE) ){
      PBPAS::Mont_PairwiseMul_OPT2_AS_R(n, APtr, BPtr, pPtr);
    }
    else{
      PBPAS::Mont_PairwiseMul_OPT2_AS_SPE_R(n, APtr, BPtr, pPtr);
    }
  }

  //======================================================================
  // pairwise mul
  //====================================================================== 
  void Mont_PairwiseMul_OPT2_AS_R(sfixn n, sfixn * APtr, sfixn * BPtr, 
				  MONTP_OPT2_AS_GENE * pPtr)
  {//------------------------------------------------------------
    //register int i;
    sfixn BRsft=pPtr->Base_Rpow;
    cilk_for(sfixn i=0; i<n; i++) 
      APtr[i]=MontMulMod_OPT2_AS_GENE(APtr[i], BPtr[i]<<BRsft, pPtr);

  }

  //------------------------------
  void Mont_PairwiseMul_OPT2_AS_SPE_R(sfixn n, sfixn * APtr, sfixn * BPtr, 
				      MONTP_OPT2_AS_GENE * pPtr)
  {//---------------------------------------------------------------
    //register int i;
    sfixn BRsft=pPtr->Base_Rpow;
    cilk_for(int i=0; i<n; i++) 
      APtr[i]=MontMulMod_OPT2_AS_GENE_SPE(APtr[i], BPtr[i]<<BRsft, pPtr);
  }
 
  //Matteo's rectangular matrix transpose---------------
  //out-of-place transpose A[i0..i1][j0..j1] into B[j0..j1][i0..i1]
  //then copy back to A
  // n: size of A
  //row major layout
  void transpose(sfixn *A, sfixn lda, sfixn *B, sfixn ldb,
		 sfixn i0, sfixn i1, sfixn j0, sfixn j1) 
  {//--------------------------------------------------------
    //const int THRESHOLD = 8; //cagnode18, cache size 4096 KB
    const int THRESHOLD = 16; //tune?
    //const int THRESHOLD = 17;

  tail:
	sfixn di = i1 - i0, dj = j1 - j0;
	if (di >= dj && di > THRESHOLD) {
	  sfixn im = (i0 + i1) / 2;
	  cilk_spawn PBPAS::transpose(A, lda, B, ldb, i0, im, j0, j1);
	  i0 = im; goto tail;
	} else if (dj > THRESHOLD) {
	  sfixn jm = (j0 + j1) / 2;
	  cilk_spawn PBPAS::transpose(A, lda, B, ldb, i0, i1, j0, jm);
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
		   sfixn j0, sfixn dj0, sfixn j1 /*, int dj1 = 0 */)
  {//------------------------------------------------------------------
    //const int THRESHOLD = 8; 
    const int THRESHOLD = 16; //cagnode18, cache line size 64bytes
                              //L1 data 32 KB, L2 4096 KB
    //const int THRESHOLD = 17;

  tail:
	sfixn di = i1 - i0, dj = j1 - j0;
	if (dj >= 2 * di && dj > THRESHOLD) {
	  sfixn dj2 = dj / 2;
	  cilk_spawn PBPAS::sqtranspose(A, lda, i0, i1, j0, dj0, j0 + dj2);
	  j0 += dj2; dj0 = 0; goto tail;
	} else if (di > THRESHOLD) {
	  sfixn di2 = di / 2;
	  cilk_spawn PBPAS::sqtranspose(A, lda, i0, i0 + di2, j0, dj0, j1);
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
  
  //Compute fast log base 2 ceiling of long int
  int ceil_log2_long(unsigned long long x)
  {//-------------------------------------------    
    static const unsigned long long t[6] = {
      0xFFFFFFFF00000000ull,
      0x00000000FFFF0000ull,
      0x000000000000FF00ull,
      0x00000000000000F0ull,
      0x000000000000000Cull,
      0x0000000000000002ull
    };
    
    int y = (((x & (x - 1)) == 0) ? 0 : 1);
    int j = 32;
 
    for (int i = 0; i < 6; i++) {
      int k = (((x & t[i]) == 0) ? 0 : j);
      y += k;
      x >>= k;
      j >>= 1;
    }
    
    return y;
  }
 

  //new DFT
  /*
   * n=2^r
   * W should be a RoosTable
   */
  void general_DFT_iter(int n, int r, 
		sfixn *A, 
		sfixn *W, 
		MONTP_OPT2_AS_GENE *pPtr,sfixn *B)
  {//-------------------------------------
    /*
    std::cout<<"--DFT_iter, n, A: "<<n<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */

    //improve when W[s]=1 by rm t
	if(n == 1)
		return;
    sfixn* u = B;
    const sfixn prime_sse[2]__attribute__((aligned(16))) = {pPtr->P,pPtr->P};

///____________________________________________________

    sfixn *Wp = W+(n<<1)-4; //W+n*2-4, point to the two elements at the last row
    for(int k=0;k<n;k+=16){
	sfixn u = A[k];
	sfixn t = A[k+1];
	A[k] = AddMod(u,t,pPtr->P);
	A[k+1] = SubMod(u,t,pPtr->P);
	u = A[k+2];
	t = A[k+3];
	A[k+2] = AddMod(u,t,pPtr->P);
	A[k+3] = SubMod(u,t,pPtr->P);
	u = A[k+4];
	t = A[k+5];
	A[k+4] = AddMod(u,t,pPtr->P);
	A[k+5] = SubMod(u,t,pPtr->P);
	u = A[k+6];
	t = A[k+7];
	A[k+6] = AddMod(u,t,pPtr->P);
	A[k+7] = SubMod(u,t,pPtr->P);
	u = A[k+8];
	t = A[k+9];
	A[k+8] = AddMod(u,t,pPtr->P);
	A[k+9] = SubMod(u,t,pPtr->P);
	u = A[k+10];
	t = A[k+11];
	A[k+10] = AddMod(u,t,pPtr->P);
	A[k+11] = SubMod(u,t,pPtr->P);
	u = A[k+12];
	t = A[k+13];
	A[k+12] = AddMod(u,t,pPtr->P);
	A[k+13] = SubMod(u,t,pPtr->P);
	u = A[k+14];
	t = A[k+15];
	A[k+14] = AddMod(u,t,pPtr->P);
	A[k+15] = SubMod(u,t,pPtr->P);
        A[k+3] = MontMulMod_OPT2_AS_GENE(A[k+3],*(Wp-3),pPtr);
        A[k+7] = MontMulMod_OPT2_AS_GENE(A[k+7],*(Wp-3),pPtr);
        A[k+11] = MontMulMod_OPT2_AS_GENE(A[k+11],*(Wp-3),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(A[k+15],*(Wp-3),pPtr);
        
        AddSubSSEModInplace(A+k,A+k+2,prime_sse);
        AddSubSSEModInplace(A+k+4,A+k+6,prime_sse);
        AddSubSSEModInplace(A+k+8,A+k+10,prime_sse);
        AddSubSSEModInplace(A+k+12,A+k+14,prime_sse);
        A[k+5] = MontMulMod_OPT2_AS_GENE(A[k+5],*(Wp-11),pPtr);
        A[k+6] = MontMulMod_OPT2_AS_GENE(A[k+6],*(Wp-10),pPtr);
        A[k+7] = MontMulMod_OPT2_AS_GENE(A[k+7],*(Wp-9),pPtr);
        A[k+13] = MontMulMod_OPT2_AS_GENE(A[k+13],*(Wp-11),pPtr);
        A[k+14] = MontMulMod_OPT2_AS_GENE(A[k+14],*(Wp-10),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(A[k+15],*(Wp-9),pPtr);
        
        AddSubSSEModInplace(A+k,A+k+4,prime_sse);
        AddSubSSEModInplace(A+k+2,A+k+6,prime_sse);
        AddSubSSEModInplace(A+k+8,A+k+12,prime_sse);
        AddSubSSEModInplace(A+k+10,A+k+14,prime_sse);
        
        A[k+9] = MontMulMod_OPT2_AS_GENE(A[k+9],*(Wp-27),pPtr);
        A[k+10] = MontMulMod_OPT2_AS_GENE(A[k+10],*(Wp-26),pPtr);
        A[k+11] = MontMulMod_OPT2_AS_GENE(A[k+11],*(Wp-25),pPtr);
        A[k+12] = MontMulMod_OPT2_AS_GENE(A[k+12],*(Wp-24),pPtr);
        A[k+13] = MontMulMod_OPT2_AS_GENE(A[k+13],*(Wp-23),pPtr);
        A[k+14] = MontMulMod_OPT2_AS_GENE(A[k+14],*(Wp-22),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(A[k+15],*(Wp-21),pPtr);
        AddSubSSEModInplace(A+k,A+k+8,prime_sse);
        AddSubSSEModInplace(A+k+2,A+k+10,prime_sse);
        AddSubSSEModInplace(A+k+4,A+k+12,prime_sse);
        AddSubSSEModInplace(A+k+6,A+k+14,prime_sse);
    }   //
    Wp = Wp - 60; //Wp-2^(i+1)

    for (int i=5; i<=r; i++){
      //traversing the tree, bottom-up
      int m = 1L<<i; // 2^i
      int m2 = m>>1; // m/2
      //std::cout<<"m: "<<m<<"\n";
      for (int k=0;k<n;k+=m){
		for (int j=0; j<m2; j+=8){
			A[k+j+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+m2], Wp[j], pPtr);
			A[k+j+1+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+1+m2], Wp[j+1], pPtr);
			A[k+j+2+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+2+m2], Wp[j+2], pPtr);
			A[k+j+3+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+3+m2], Wp[j+3], pPtr);
			A[k+j+4+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+4+m2], Wp[j+4], pPtr);
			A[k+j+5+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+5+m2], Wp[j+5], pPtr);
			A[k+j+6+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+6+m2], Wp[j+6], pPtr);
			A[k+j+7+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+7+m2], Wp[j+7], pPtr);
			AddSubSSEModInplace(A+k+j,A+k+j+m2,prime_sse);
			AddSubSSEModInplace(A+k+j+2,A+k+j+2+m2,prime_sse);
			AddSubSSEModInplace(A+k+j+4,A+k+j+4+m2,prime_sse);
			AddSubSSEModInplace(A+k+j+6,A+k+j+6+m2,prime_sse);
		}	
      }
      Wp = Wp - (1L<<(i+1)); //Wp-2^(i+1)
    }
  }
  void small_DFT_iter(int n, int r, 
		sfixn *A, 
		sfixn *W, 
		MONTP_OPT2_AS_GENE *pPtr,sfixn *B)
  {//-------------------------------------
    /*
    std::cout<<"--DFT_iter, n, A: "<<n<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */

    //improve when W[s]=1 by rm t
	if(n == 1)
		return;
    sfixn* u = B;
    const sfixn prime_sse[2]__attribute__((aligned(16))) = {pPtr->P,pPtr->P};

///____________________________________________________

    sfixn *Wp = W+(n<<1)-4; //W+n*2-4, point to the two elements at the last row
    for(int k=0;k<n;k+=2){
	sfixn u = A[k];
	sfixn v = A[k+1];
	A[k] = AddMod(u,v,pPtr->P);
	A[k+1] = SubMod(u,v,pPtr->P);
    }
    Wp = Wp - 4; //Wp-2^(i+1)

    for (int i=2; i<=r; i++){
      //traversing the tree, bottom-up
      int m = 1L<<i; // 2^i
      int m2 = m>>1; // m/2
      //std::cout<<"m: "<<m<<"\n";
      for (int k=0;k<n;k+=m){
		for (int j=0; j<m2; j+=2){
			A[k+j+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+m2], Wp[j], pPtr);
			A[k+j+1+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+1+m2], Wp[j+1], pPtr);
			AddSubSSEModInplace(A+k+j,A+k+j+m2,prime_sse);
		}	
      }
      Wp = Wp - (1L<<(i+1)); //Wp-2^(i+1)
    }
  }
  void DFT_eff(int n, int r, 
	       sfixn *A, 
	       sfixn *W, 
	       MONTP_OPT2_AS_GENE *pPtr,
	       int H, int *RevBidMap,
	       sfixn *B,sfixn whichprime)
  {//-------------------------------------
    /*
    std::cout<<"--DFT_eff, n, A:" <<n<<"\n";
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */
    if (n==1) { //rm for H>1     
      return;
    }
    if(whichprime ==1){
#if FURER
     // printf("FURER DFT 1\n");
    	      FURERPBPAS1::DFT_eff_p1(n,r,A,W,B);
#else
     //       printf("Other DFT 1\n");
    	      PBPAS1::DFT_eff_p1(n,r,A,W,B);
#endif
          }
    else if(whichprime == 2)
#if FURER
    	      FURERPBPAS2::DFT_eff_p2(n,r,A,W,B);
#else
    	      PBPAS2::DFT_eff_p2(n,r,A,W,B);
#endif
    else{
	    if (n>H){
	      DFT_rec(n,r,A,W,pPtr,H,RevBidMap,B);
	    }
	    else{   
	      
	      ArrayBitReversal(n,A,RevBidMap);
	      /*
	      std::cout<<"--DFT_eff, after ArrayBitReversal:" <<n<<"\n";
	      for(int i=0; i<n; ++i)
		std::cout<<A[i]<<", ";
	      std::cout<<std::endl;
	      */
	      if (n>=16){
		      general_DFT_iter(n,r,A,W,pPtr,B);      
		}
		else
			small_DFT_iter(n,r,A,W,pPtr,B);
	    }
   }
  }

  /*
   * n=2^r
   */
  void DFT_rec(int n, int r, 
	       sfixn *A, 
	       sfixn *W, 
	       MONTP_OPT2_AS_GENE *pPtr,
	       int H, int *RevBidMap,
	       sfixn *B)
  {//-------------------------------------

    /*
    std::cout<<"--DFT_rec, n, A: "<<n<<std::endl;
    for(int i=0; i<n; ++i)
	std::cout<<A[i]<<", ";
      std::cout<<std::endl;
    */
    if (n==1) {    
      return;
    }
    else if(n<=H){
      ArrayBitReversal(n,A,RevBidMap);
      assert(H>=16);
      general_DFT_iter(n,r,A,W,pPtr,B);
      return;
    }
    int n2 = n>>1; //n/2
    int r1 = r-1;
    
    //shuffling data
    Shuffle(n, A, B);

    /*
    std::cout<<"--DFT_rec, after shuffle: "<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */

    sfixn *W2 = W+n; 
    DFT_rec(n2, r1, A,    W2, pPtr, H, RevBidMap, B);
    DFT_rec(n2, r1, A+n2, W2, pPtr, H, RevBidMap, B);

    /*
    std::cout<<"--DFT_rec, n2, A: "<<n2<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */

    const sfixn prime_sse[2]__attribute__((aligned(16))) = {pPtr->P,pPtr->P};

    for (int k=n2; k<(n2<<1); k+=8){
	//sfixn u[8] = {A[k],A[k+1],A[k+2],A[k+3],A[k+4],A[k+5],A[k+6],A[k+7]};
	//sfixn v[8] = {A[k],A[k+1],A[k+2],A[k+3],A[k+4],A[k+5],A[k+6],A[k+7]};
	#if ASM
	unrolled8MontMul(A+k,W+k-n2,pPtr);
	unrolled8AddSubSSEMod(A+k-n2,A+k,prime_sse);
	#else
	A[k+0]=MontMulMod_OPT2_AS_GENE(A[k+0],W[k-n2+0],pPtr);
	A[k+1]=MontMulMod_OPT2_AS_GENE(A[k+1],W[k-n2+1],pPtr);
	A[k+2]=MontMulMod_OPT2_AS_GENE(A[k+2],W[k-n2+2],pPtr);
	A[k+3]=MontMulMod_OPT2_AS_GENE(A[k+3],W[k-n2+3],pPtr);
	A[k+4]=MontMulMod_OPT2_AS_GENE(A[k+4],W[k-n2+4],pPtr);
	A[k+5]=MontMulMod_OPT2_AS_GENE(A[k+5],W[k-n2+5],pPtr);
	A[k+6]=MontMulMod_OPT2_AS_GENE(A[k+6],W[k-n2+6],pPtr);
	A[k+7]=MontMulMod_OPT2_AS_GENE(A[k+7],W[k-n2+7],pPtr);
	sfixn t = A[k-n2];
	sfixn u = A[k];
	A[k-n2] = AddMod(t,u,pPtr->P);
	A[k] = SubMod(t,u,pPtr->P);
	t = A[k+1-n2];
	u = A[k+1];
	A[k+1-n2] = AddMod(t,u,pPtr->P);
	A[k+1] = SubMod(t,u,pPtr->P);
	t = A[k+2-n2];
	u = A[k+2];
	A[k+2-n2] = AddMod(t,u,pPtr->P);
	A[k+2] = SubMod(t,u,pPtr->P);
	t = A[k+3-n2];
	u = A[k+3];
	A[k+3-n2] = AddMod(t,u,pPtr->P);
	A[k+3] = SubMod(t,u,pPtr->P);
	t = A[k+4-n2];
	u = A[k+4];
	A[k+4-n2] = AddMod(t,u,pPtr->P);
	A[k+4] = SubMod(t,u,pPtr->P);
	t = A[k+5-n2];
	u = A[k+5];
	A[k+5-n2] = AddMod(t,u,pPtr->P);
	A[k+5] = SubMod(t,u,pPtr->P);
	t = A[k+6-n2];
	u = A[k+6];
	A[k+6-n2] = AddMod(t,u,pPtr->P);
	A[k+6] = SubMod(t,u,pPtr->P);
	t = A[k+7-n2];
	u = A[k+7];
	A[k+7-n2] = AddMod(t,u,pPtr->P);
	A[k+7] = SubMod(t,u,pPtr->P);
	#endif
    }    
    //for (int k=0; k<n2; k++){
    //  //which w? tmpw=rootsPtr[partialBitRev(i,power)];
    //  sfixn t = A[k+n2];
    //  sfixn u = A[k];
    //  //sfixn t = MontMulMod_OPT2_AS_GENE(v, W[k], pPtr);
    //  A[k] = AddMod(u, t, pPtr->P);
    //  A[k+n2] = SubMod(u, t, pPtr->P);
    //}    
  }
  
  /*
   * size of B is n/2
   */
  void Shuffle(int n, sfixn *A, sfixn *B)
  {//-------------------------------------

    /*
    std::cout<<"--Shuffle, n, A in: "<<n<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */
    //test allocate B locally???
    int n2 = n>>1;
    //shuffling data
    
    for (int i=0; i<n2; i++){
      int i2 = i<<1; //2i
      A[i] = A[i2];
      B[i] = A[i2+1];	
    }
    sfixn *A2 = A + n2;      
    memcpy(A2, B, n2*(sizeof(sfixn)));   
  }

  void ArrayBitReversal(int n, sfixn *A, 
			int *RevBidMap)
  {//------------------------------------- 
    /*
    std::cout<<"--ArrayBitReversal, n, A in: "<<n<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */
    for (int i=0; i<n; i++){
      //pre-computed table
      int j = RevBidMap[i];
      if (i<j){ //change position
	sfixn t = A[i];
	A[i] = A[j];
	A[j] = t;
      }
    }
    #if ASM
    for(int i=0;i<n;i+=4){
    	sfixn t = A[i+1];
	A[i+1] = A[i+2];
	A[i+2] = t;
    }
    #endif
  }

  /*
   * n=2^r
   */

  void safe_DFT_iter(int n, int r, 
		sfixn *A, 
		sfixn *W, 
		MONTP_OPT2_AS_GENE *pPtr,sfixn *B)
  {//-------------------------------------
    /*
    std::cout<<"--DFT_iter, n, A: "<<n<<std::endl;
    for(int i=0; i<n; ++i)
      std::cout<<A[i]<<", ";
    std::cout<<std::endl;
    */

    //improve when W[s]=1 by rm t

    sfixn* u = B;

///____________________________________________________

    const sfixn prime_sse[2]__attribute__((aligned(16))) = {pPtr->P,pPtr->P};
    sfixn *Wp = W+(n<<1)-4; //W+n*2-4, point to the two elements at the last row
    for(int k=0;k<n;k+=4){
	AddSubSSEMod(u,u+2,A+k,A+k+2,prime_sse);
	A[k+0] = u[0];
	A[k+1] = u[2];
	A[k+2] = u[1];
        A[k+3] = u[3];
    }
    Wp = Wp - 4; //Wp-2^(i+1)

    for (int i=2; i<=r; i++){
      //traversing the tree, bottom-up
      int m = 1L<<i; // 2^i
      int m2 = m>>1; // m/2
      //std::cout<<"m: "<<m<<"\n";
      for (int k=0;k<n;k+=m){
		for (int j=0; j<m2; j+=2){
			A[k+j+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+m2], Wp[j], pPtr);
			A[k+j+1+m2] = MontMulMod_OPT2_AS_GENE(A[k+j+1+m2], Wp[j+1], pPtr);
			AddSubSSEModInplace(A+k+j,A+k+j+m2,prime_sse);
		}	
      }
      Wp = Wp - (1L<<(i+1)); //Wp-2^(i+1)
    }
  }
  void DFT_iter(int n, int r, 
		sfixn *A, 
		sfixn *W, 
		MONTP_OPT2_AS_GENE *pPtr,
		sfixn *B)
  {//-------------------------------------

    //improve when W[s]=1 by rm t

    sfixn* u = B;

///____________________________________________________

    const sfixn prime_sse[2]__attribute__((aligned(16))) = {pPtr->P,pPtr->P};
    const unsigned char shuffle_array[16]__attribute__((aligned(16))) = {4,5,6,7,4,5,6,7,12,13,14,15,12,13,14,15};
    sfixn *Wp = W+(n<<1)-4; //W+n*2-4, point to the two elements at the last row
    
    for(int k=0;k<n;k+=16){
	#if ASM
	AddSubSSEMod(u,u+2,A+k,A+k+2,prime_sse);
	AddSubSSEMod(u+4,u+6,A+k+4,A+k+6,prime_sse);
	AddSubSSEMod(u+8,u+10,A+k+8,A+k+10,prime_sse);
	AddSubSSEMod(u+12,u+14,A+k+12,A+k+14,prime_sse);
        A[k+3] = MontMulMod_OPT2_AS_GENE(u[3],*(Wp-3),pPtr);
        A[k+7] = MontMulMod_OPT2_AS_GENE(u[7],*(Wp-3),pPtr);
        A[k+11] = MontMulMod_OPT2_AS_GENE(u[11],*(Wp-3),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(u[15],*(Wp-3),pPtr);
	A[k+0] = u[0];
	A[k+1] = u[2];
	A[k+2] = u[1];
	A[k+4] = u[4];
	A[k+5] = u[6];
	A[k+6] = u[5];
	A[k+8] = u[8];
	A[k+9] = u[10];
	A[k+10] = u[9];
	A[k+12] = u[12];
	A[k+13] = u[14];
	A[k+14] = u[13];
	#else
	sfixn u = A[k];
	sfixn t = A[k+1];
	A[k] = AddMod(u,t,pPtr->P);
	A[k+1] = SubMod(u,t,pPtr->P);
	u = A[k+2];
	t = A[k+3];
	A[k+2] = AddMod(u,t,pPtr->P);
	A[k+3] = SubMod(u,t,pPtr->P);
	u = A[k+4];
	t = A[k+5];
	A[k+4] = AddMod(u,t,pPtr->P);
	A[k+5] = SubMod(u,t,pPtr->P);
	u = A[k+6];
	t = A[k+7];
	A[k+6] = AddMod(u,t,pPtr->P);
	A[k+7] = SubMod(u,t,pPtr->P);
	u = A[k+8];
	t = A[k+9];
	A[k+8] = AddMod(u,t,pPtr->P);
	A[k+9] = SubMod(u,t,pPtr->P);
	u = A[k+10];
	t = A[k+11];
	A[k+10] = AddMod(u,t,pPtr->P);
	A[k+11] = SubMod(u,t,pPtr->P);
	u = A[k+12];
	t = A[k+13];
	A[k+12] = AddMod(u,t,pPtr->P);
	A[k+13] = SubMod(u,t,pPtr->P);
	u = A[k+14];
	t = A[k+15];
	A[k+14] = AddMod(u,t,pPtr->P);
	A[k+15] = SubMod(u,t,pPtr->P);
        A[k+3] = MontMulMod_OPT2_AS_GENE(A[k+3],*(Wp-3),pPtr);
        A[k+7] = MontMulMod_OPT2_AS_GENE(A[k+7],*(Wp-3),pPtr);
        A[k+11] = MontMulMod_OPT2_AS_GENE(A[k+11],*(Wp-3),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(A[k+15],*(Wp-3),pPtr);
	#endif
        
        AddSubSSEModInplace(A+k,A+k+2,prime_sse);
        AddSubSSEModInplace(A+k+4,A+k+6,prime_sse);
        AddSubSSEModInplace(A+k+8,A+k+10,prime_sse);
        AddSubSSEModInplace(A+k+12,A+k+14,prime_sse);
        A[k+5] = MontMulMod_OPT2_AS_GENE(A[k+5],*(Wp-11),pPtr);
        A[k+6] = MontMulMod_OPT2_AS_GENE(A[k+6],*(Wp-10),pPtr);
        A[k+7] = MontMulMod_OPT2_AS_GENE(A[k+7],*(Wp-9),pPtr);
        A[k+13] = MontMulMod_OPT2_AS_GENE(A[k+13],*(Wp-11),pPtr);
        A[k+14] = MontMulMod_OPT2_AS_GENE(A[k+14],*(Wp-10),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(A[k+15],*(Wp-9),pPtr);
        
        AddSubSSEModInplace(A+k,A+k+4,prime_sse);
        AddSubSSEModInplace(A+k+2,A+k+6,prime_sse);
        AddSubSSEModInplace(A+k+8,A+k+12,prime_sse);
        AddSubSSEModInplace(A+k+10,A+k+14,prime_sse);
        
        A[k+9] = MontMulMod_OPT2_AS_GENE(A[k+9],*(Wp-27),pPtr);
        A[k+10] = MontMulMod_OPT2_AS_GENE(A[k+10],*(Wp-26),pPtr);
        A[k+11] = MontMulMod_OPT2_AS_GENE(A[k+11],*(Wp-25),pPtr);
        A[k+12] = MontMulMod_OPT2_AS_GENE(A[k+12],*(Wp-24),pPtr);
        A[k+13] = MontMulMod_OPT2_AS_GENE(A[k+13],*(Wp-23),pPtr);
        A[k+14] = MontMulMod_OPT2_AS_GENE(A[k+14],*(Wp-22),pPtr);
        A[k+15] = MontMulMod_OPT2_AS_GENE(A[k+15],*(Wp-21),pPtr);
        AddSubSSEModInplace(A+k,A+k+8,prime_sse);
        AddSubSSEModInplace(A+k+2,A+k+10,prime_sse);
        AddSubSSEModInplace(A+k+4,A+k+12,prime_sse);
        AddSubSSEModInplace(A+k+6,A+k+14,prime_sse);
    }   //
    Wp = Wp - 60; //Wp-2^(i+1)
    sfixn* Wt = Wp;
    for(int k=0;k<n;k+=128){
	Wp = Wt;
	A[k+17] = MontMulMod_OPT2_AS_GENE(A[k+17],*(Wp+1),pPtr);
	A[k+49] = MontMulMod_OPT2_AS_GENE(A[k+49],*(Wp+1),pPtr);
	A[k+81] = MontMulMod_OPT2_AS_GENE(A[k+81],*(Wp+1),pPtr);
	A[k+113] = MontMulMod_OPT2_AS_GENE(A[k+113],*(Wp+1),pPtr);
	AddSubSSEModInplace(A+k,A+k+16,prime_sse);
	AddSubSSEModInplace(A+k+32,A+k+48,prime_sse);
	AddSubSSEModInplace(A+k+64,A+k+80,prime_sse);
	AddSubSSEModInplace(A+k+96,A+k+112,prime_sse);
	Wp = Wp-64;

	A[k+33] = MontMulMod_OPT2_AS_GENE(A[k+33],*(Wp+1),pPtr);
	A[k+97] = MontMulMod_OPT2_AS_GENE(A[k+97],*(Wp+1),pPtr);
	A[k+48] = MontMulMod_OPT2_AS_GENE(A[k+48],*(Wp+16),pPtr);
	A[k+112] = MontMulMod_OPT2_AS_GENE(A[k+112],*(Wp+16),pPtr);
	A[k+49] = MontMulMod_OPT2_AS_GENE(A[k+49],*(Wp+17),pPtr);
	A[k+113] = MontMulMod_OPT2_AS_GENE(A[k+113],*(Wp+17),pPtr);
	AddSubSSEModInplace(A+k,A+k+32,prime_sse);
	AddSubSSEModInplace(A+k+16,A+k+48,prime_sse);
	AddSubSSEModInplace(A+k+64,A+k+96,prime_sse);
	AddSubSSEModInplace(A+k+80,A+k+112,prime_sse);
	Wp = Wp - 128;
	A[k+65] = MontMulMod_OPT2_AS_GENE(A[k+65],*(Wp+1),pPtr);
	A[k+80] = MontMulMod_OPT2_AS_GENE(A[k+80],*(Wp+16),pPtr);
	A[k+81] = MontMulMod_OPT2_AS_GENE(A[k+81],*(Wp+17),pPtr);
	A[k+96] = MontMulMod_OPT2_AS_GENE(A[k+96],*(Wp+32),pPtr);
	A[k+97] = MontMulMod_OPT2_AS_GENE(A[k+97],*(Wp+33),pPtr);
	A[k+112] = MontMulMod_OPT2_AS_GENE(A[k+112],*(Wp+48),pPtr);
	A[k+113] = MontMulMod_OPT2_AS_GENE(A[k+113],*(Wp+49),pPtr);
	AddSubSSEModInplace(A+k,A+k+64,prime_sse);
	AddSubSSEModInplace(A+k+16,A+k+80,prime_sse);
	AddSubSSEModInplace(A+k+32,A+k+96,prime_sse);
	AddSubSSEModInplace(A+k+48,A+k+112,prime_sse);
	for(int j=2;j<16;j+=2){
		Wp = Wt;
		#if ASM
		unrolled8MontMul2root16(A+k+j+16,*(Wp+j),*(Wp+j+1),pPtr);
		#else
		A[k+j+16] = MontMulMod_OPT2_AS_GENE(A[k+j+16],*(Wp+j),pPtr);
		A[k+j+3*16] = MontMulMod_OPT2_AS_GENE(A[k+j+3*16],*(Wp+j),pPtr);
		A[k+j+5*16] = MontMulMod_OPT2_AS_GENE(A[k+j+5*16],*(Wp+j),pPtr);
		A[k+j+7*16] = MontMulMod_OPT2_AS_GENE(A[k+j+7*16],*(Wp+j),pPtr);
		A[k+j+16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+16+1],*(Wp+j+1),pPtr);
		A[k+j+3*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+3*16+1],*(Wp+j+1),pPtr);
		A[k+j+5*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+5*16+1],*(Wp+j+1),pPtr);
		A[k+j+7*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+7*16+1],*(Wp+j+1),pPtr);
		#endif
		AddSubSSEModInplace(A+k+j,A+k+j+16,prime_sse);
		AddSubSSEModInplace(A+k+j+32,A+k+j+48,prime_sse);
		AddSubSSEModInplace(A+k+j+64,A+k+j+80,prime_sse);
		AddSubSSEModInplace(A+k+j+96,A+k+j+112,prime_sse);
		Wp = Wp-64;

		#if ASM
		unrolled8MontMul4root16(A+k+j+32,*(Wp+j),*(Wp+j+1),*(Wp+j+16),*(Wp+j+17),pPtr);
		#else
		A[k+j+2*16] = MontMulMod_OPT2_AS_GENE(A[k+j+2*16],*(Wp+j),pPtr);
		A[k+j+6*16] = MontMulMod_OPT2_AS_GENE(A[k+j+6*16],*(Wp+j),pPtr);
		A[k+j+2*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+2*16+1],*(Wp+j+1),pPtr);
		A[k+j+6*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+6*16+1],*(Wp+j+1),pPtr);
		A[k+j+3*16] = MontMulMod_OPT2_AS_GENE(A[k+j+3*16],*(Wp+j+16),pPtr);
		A[k+j+7*16] = MontMulMod_OPT2_AS_GENE(A[k+j+7*16],*(Wp+j+16),pPtr);
		A[k+j+3*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+3*16+1],*(Wp+j+16+1),pPtr);
		A[k+j+7*16+1] = MontMulMod_OPT2_AS_GENE(A[k+j+7*16+1],*(Wp+j+16+1),pPtr);
		#endif
		AddSubSSEModInplace(A+k+j,A+k+j+32,prime_sse);
		AddSubSSEModInplace(A+k+j+16,A+k+j+48,prime_sse);
		AddSubSSEModInplace(A+k+j+64,A+k+j+96,prime_sse);
		AddSubSSEModInplace(A+k+j+80,A+k+j+112,prime_sse);
		Wp = Wp - 128;
		A[k+j+64] = MontMulMod_OPT2_AS_GENE(A[k+j+64],*(Wp+j),pPtr);
		A[k+j+65] = MontMulMod_OPT2_AS_GENE(A[k+j+65],*(Wp+j+1),pPtr);
		A[k+j+80] = MontMulMod_OPT2_AS_GENE(A[k+j+80],*(Wp+j+16),pPtr);
		A[k+j+81] = MontMulMod_OPT2_AS_GENE(A[k+j+81],*(Wp+j+17),pPtr);
		A[k+j+96] = MontMulMod_OPT2_AS_GENE(A[k+j+96],*(Wp+j+32),pPtr);
		A[k+j+97] = MontMulMod_OPT2_AS_GENE(A[k+j+97],*(Wp+j+33),pPtr);
		A[k+j+112] = MontMulMod_OPT2_AS_GENE(A[k+j+112],*(Wp+j+48),pPtr);
		A[k+j+113] = MontMulMod_OPT2_AS_GENE(A[k+j+113],*(Wp+j+49),pPtr);
		AddSubSSEModInplace(A+k+j,A+k+j+64,prime_sse);
		AddSubSSEModInplace(A+k+j+16,A+k+j+80,prime_sse);
		AddSubSSEModInplace(A+k+j+32,A+k+j+96,prime_sse);
		AddSubSSEModInplace(A+k+j+48,A+k+j+112,prime_sse);
	}
    }
    Wp = Wt-64 - 128 - 256;
    Wt = Wp;
    int i=1;
    const int jump = 16*(1L<<(3));
    for(int k=0;k<n;k+=jump<<3){
	printf("%d ahoiy\n",k);
	for(int j=0;j<128;j+=2){
		printf("%d ahoiy\n",j);
		Wp = Wt;
		A[k+j+128] = MontMulMod_OPT2_AS_GENE(A[k+j+128],*(Wp+j),pPtr);
		A[k+j+3*128] = MontMulMod_OPT2_AS_GENE(A[k+j+3*128],*(Wp+j),pPtr);
		A[k+j+5*128] = MontMulMod_OPT2_AS_GENE(A[k+j+5*128],*(Wp+j),pPtr);
		A[k+j+7*128] = MontMulMod_OPT2_AS_GENE(A[k+j+7*128],*(Wp+j),pPtr);
		A[k+j+128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+128+1],*(Wp+j+1),pPtr);
		A[k+j+3*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+3*128+1],*(Wp+j+1),pPtr);
		A[k+j+5*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+5*128+1],*(Wp+j+1),pPtr);
		A[k+j+7*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+7*128+1],*(Wp+j+1),pPtr);
		AddSubSSEModInplace(A+k+j,A+k+j+128,prime_sse);
		AddSubSSEModInplace(A+k+j+2*128,A+k+j+3*128,prime_sse);
		AddSubSSEModInplace(A+k+j+4*128,A+k+j+5*128,prime_sse);
		AddSubSSEModInplace(A+k+j+6*128,A+k+j+7*128,prime_sse);
		Wp = Wp-(1L<<(3*i+6));

		A[k+j+2*128] = MontMulMod_OPT2_AS_GENE(A[k+j+2*128],*(Wp+j),pPtr);
		A[k+j+6*128] = MontMulMod_OPT2_AS_GENE(A[k+j+6*128],*(Wp+j),pPtr);
		A[k+j+2*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+2*128+1],*(Wp+j+1),pPtr);
		A[k+j+6*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+6*128+1],*(Wp+j+1),pPtr);
		A[k+j+3*128] = MontMulMod_OPT2_AS_GENE(A[k+j+3*128],*(Wp+j+128),pPtr);
		A[k+j+7*128] = MontMulMod_OPT2_AS_GENE(A[k+j+7*128],*(Wp+j+128),pPtr);
		A[k+j+3*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+3*128+1],*(Wp+j+128+1),pPtr);
		A[k+j+7*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+7*128+1],*(Wp+j+128+1),pPtr);
		AddSubSSEModInplace(A+k+j,A+k+j+2*128,prime_sse);
		AddSubSSEModInplace(A+k+j+128,A+k+j+3*128,prime_sse);
		AddSubSSEModInplace(A+k+j+4*128,A+k+j+6*128,prime_sse);
		AddSubSSEModInplace(A+k+j+5*128,A+k+j+7*128,prime_sse);
		Wp = Wp - (1L<<(3*i+7));
		A[k+j+4*128] = MontMulMod_OPT2_AS_GENE(A[k+j+4*128],*(Wp+j),pPtr);
		A[k+j+4*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+4*128+1],*(Wp+j+1),pPtr);
		A[k+j+5*128] = MontMulMod_OPT2_AS_GENE(A[k+j+5*128],*(Wp+j+128),pPtr);
		A[k+j+5*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+5*128+1],*(Wp+j+128+1),pPtr);
		A[k+j+6*128] = MontMulMod_OPT2_AS_GENE(A[k+j+6*128],*(Wp+j+2*128),pPtr);
		A[k+j+6*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+6*128+1],*(Wp+j+2*128+1),pPtr);
		A[k+j+7*128] = MontMulMod_OPT2_AS_GENE(A[k+j+7*128],*(Wp+j+3*128),pPtr);
		A[k+j+7*128+1] = MontMulMod_OPT2_AS_GENE(A[k+j+7*128+1],*(Wp+j+3*128+1),pPtr);
		AddSubSSEModInplace(A+k+j,A+k+j+4*128,prime_sse);
		AddSubSSEModInplace(A+k+j+128,A+k+j+5*128,prime_sse);
		AddSubSSEModInplace(A+k+j+2*128,A+k+j+6*128,prime_sse);
		AddSubSSEModInplace(A+k+j+3*128,A+k+j+7*128,prime_sse);
	}
    }
    //Wp = Wp - (1L<<(3*i+8)); //Wp-2^(i+1)
    Wp = Wt - 512-1024-2048;
    for (int i=11; i<=r; i++){
      //traversing the tree, bottom-up
      int m = 1L<<i; // 2^i
      int m2 = m>>1; // m/2
      //std::cout<<"m: "<<m<<"\n";
      for (int k=0;k<n;k+=m){
        	for (int j=0; j<m2; j+=16){
        		#if ASM
        		unrolled8MontMul(A+k+j+m2,Wp+j,pPtr);
        		unrolled8AddSubSSEMod(A+k+j,A+k+j+m2,prime_sse);
        		unrolled8MontMul(A+k+j+8+m2,Wp+j+8,pPtr);
        		unrolled8AddSubSSEMod(A+k+8+j,A+k+8+j+m2,prime_sse);
        		#else
        		A[k+j+m2+0] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+0],*(Wp+j+0),pPtr);
        		A[k+j+m2+1] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+1],*(Wp+j+1),pPtr);
        		A[k+j+m2+2] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+2],*(Wp+j+2),pPtr);
        		A[k+j+m2+3] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+3],*(Wp+j+3),pPtr);
        		A[k+j+m2+4] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+4],*(Wp+j+4),pPtr);
        		A[k+j+m2+5] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+5],*(Wp+j+5),pPtr);
        		A[k+j+m2+6] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+6],*(Wp+j+6),pPtr);
        		A[k+j+m2+7] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+7],*(Wp+j+7),pPtr);
        		A[k+j+m2+8] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+8],*(Wp+j+8),pPtr);
        		A[k+j+m2+9] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+9],*(Wp+j+9),pPtr);
        		A[k+j+m2+10] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+10],*(Wp+j+10),pPtr);
        		A[k+j+m2+11] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+11],*(Wp+j+11),pPtr);
        		A[k+j+m2+12] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+12],*(Wp+j+12),pPtr);
        		A[k+j+m2+13] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+13],*(Wp+j+13),pPtr);
        		A[k+j+m2+14] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+14],*(Wp+j+14),pPtr);
        		A[k+j+m2+15] = MontMulMod_OPT2_AS_GENE(A[k+j+m2+15],*(Wp+j+15),pPtr);
        		AddSubSSEModInplace(A+k+j+0,A+k+j+m2+0,prime_sse);
        		AddSubSSEModInplace(A+k+j+2,A+k+j+m2+2,prime_sse);
        		AddSubSSEModInplace(A+k+j+4,A+k+j+m2+4,prime_sse);
        		AddSubSSEModInplace(A+k+j+6,A+k+j+m2+6,prime_sse);
        		AddSubSSEModInplace(A+k+j+8,A+k+j+m2+8,prime_sse);
        		AddSubSSEModInplace(A+k+j+10,A+k+j+m2+10,prime_sse);
        		AddSubSSEModInplace(A+k+j+12,A+k+j+m2+12,prime_sse);
        		AddSubSSEModInplace(A+k+j+14,A+k+j+m2+14,prime_sse);
        		#endif
        	}	
      }
      Wp = Wp - (1L<<(i+1)); //Wp-2^(i+1)
    }
  }

  void RevBitIncr( int *n, int bit )
  {//--------------------------------
    do{
      bit >>= 1;
      *n ^= bit;
    } while( (*n & bit) == 0 && bit != 1 );
  }
  
  /*
   * #A = max
   *A In: 
   * 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
   *A Out: 
   * 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
   */
  void RevBitInd(int max, int *A)
  {//-----------------------------
    //int max = 0x100; //256
    //int max = 16;
    //printf( "max: %d\n", max);
    if (max==1) return;
    int i, j;
    for( i = 0, j = 0; i != max; ++i, RevBitIncr(&j, max) ) {
      if( i < j ){
	int t = A[i];
	A[i]  = j;
	A[j]  = t;
	//printf( "%02x <-> %02x\n", i, j );
	//printf( "%d <-> %d\n", i, j );
      }
    }
  }
   
  /*
   *n = 2^r
   *e.g. n = 16, contents of rootsTable:
   * 1, w, w^2, w^3, w^4, w^5, w^6, w^7, w^8, w^9, w^10, w^11, w^12, w^13, w^14, w^15
   * 1, w^2, w^4, w^6, w^8, w^10, w^12, w^14
   * 1, w^4, w^8, w^12
   * 1, w^8
   */
  void RootsTable(int n, int r, 
		  sfixn *T,
		  MONTP_OPT2_AS_GENE *pPtr)
  {//---------------------------------------------------------------
    
    PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(r, n, T, pPtr);
    
    // std::cout<<"EX_Mont_GetNthRoots_OPT2_AS_GENE T Out: "<<std::endl;
    // for(int i=0; i<n; ++i)
    //   std::cout<<T[i]<<", ";
    // std::cout<<std::endl;
    
    sfixn *T_ptr = T;
    int ni = n;
    for (int i=1; i<r; i++){
      //std::cout<<"i: "<<i<<"\n";
      T_ptr = T_ptr + ni;
      ni = n>>i; //n/2^i
      for (int j=0; j<ni; j++){
	T_ptr[j] = T[(j<<i)]; //j*2^i
      }
    }
  }

  void RootsTableSpe(int n, int r, 
		  sfixn *T,MONTP_OPT2_AS_GENE *p,
		  int m)
  {
#if FURER
	if (m==1){
		FURERPBPAS1::RootsTableFurer(n,r,T);
    //printf("FURER Table 1\n");
	}
	else
		RootsTable(n,r,T,p);
#else
  //  printf("else\n");
		RootsTable(n,r,T,p);
#endif
  }
  // Ti is the reversal of T except the first column
  /*
   *n = 2^r
   *e.g. n = 16, contents of inverserootsTable:
   *  m, m^2, m^3, m^4, m^5, m^6, m^7, m^8, m^9, m^10, m^11, m^12, m^13, m^14, m^15
   *  m^2, m^4, m^6, m^8, m^10, m^12, m^14
   *  m^4, m^8, m^12
   *  m^8
   */



  void InverseRootsTable(int n, sfixn *T, sfixn *Ti)
  {//-------------------------------------------------
    //std::cout<<"--InverseRootsTable, n: "<<n<<std::endl;
    Ti[0] = T[0];
    int m=n;
    //int end=n;
    int i=0; 
    sfixn *Ti_p = Ti;
    sfixn *T_p = T;
    while(m>=2){
      Ti_p[0] = T_p[0];
      for (i=1; i<m; i++){
	Ti_p[m-i] = T_p[i];
      }
      Ti_p+=m;
      T_p +=m;
      m=(m>>1);
      //end+=m;
    }
    //std::cout<<end<<"\n";
  }
  
  /*
   *fft and inv of d dim (y)
   */
  void RootsTable2(int n, int r, 
		   sfixn *T, sfixn *Ti, 
		   MONTP_OPT2_AS_GENE *pPtr)
  {//----------------------------------------------
    RootsTable(n,r,T,pPtr);
    InverseRootsTable(n,T,Ti);
  }


  /*
   *K=n/2
   *scaling factors, roots tables fft and inv for K (x-dim)
   */
  void Weight_RootsTable2(int n, int r, int K,
			  sfixn *Th, sfixn *T, sfixn *Ti, 
			  MONTP_OPT2_AS_GENE *pPtr)
  {//--------------------------------------------------------
    PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(r, n, Th, pPtr);
    Th[K<<1] = Th[0];

    for (int i=0; i<K; i++){
      T[i] = Th[i<<1];
    }
    sfixn *T_ptr = T;
    int ni = K;
    for (int i=1; i<r-1; i++){
      //std::cout<<"i: "<<i<<"\n";
      T_ptr = T_ptr + ni;
      ni = K>>i; //K/2^i
      for (int j=0; j<ni; j++){
	T_ptr[j] = T[(j<<i)]; //j*2^i
      }
    }
    
    InverseRootsTable(K,T,Ti);
  }
  void InvDFT_eff_keepMontgomery(int n, int r, 
		  sfixn *A, 
		  sfixn *W, 
		  MONTP_OPT2_AS_GENE *pPtr,
		  int H, int *RevBidMap,
		  sfixn *B,sfixn invn,sfixn whichprime)
  {//-------------------------------------
    if(whichprime ==1)
#if FURER
    	      FURERPBPAS1::InvDFTKeepMont_eff_p1(n,r,A,W,B,invn);
#else
    	      PBPAS1::InvDFTKeepMont_eff_p1(n,r,A,W,B,invn);
#endif
    else if(whichprime == 2)
#if FURER
    	      FURERPBPAS2::InvDFTKeepMont_eff_p2(n,r,A,W,B,invn);
#else
    	      PBPAS2::InvDFTKeepMont_eff_p2(n,r,A,W,B,invn);
#endif
    else{
	    DFT_eff(n,r,A,W,pPtr,H,RevBidMap,B,whichprime);
	    //std::cout<<"InvDFT_eff before inv, A[0]: "<<A[0]<<"\n";
	    /*
	    std::cout<<"\n--- InvDFT_eff before div by n : \n";
	    for(int i=0; i<n; ++i)
	      std::cout<<A[i]<<", ";
	    std::cout<<std::endl;
	    */

	    for(int i=0; i<n; i++) {
			A[i]=MontMulMod_OPT2_AS_GENE(A[i], invn, pPtr);
		}
	}
}
  void InvDFT_eff(int n, int r, 
		  sfixn *A, 
		  sfixn *W, 
		  MONTP_OPT2_AS_GENE *pPtr,
		  int H, int *RevBidMap,
		  sfixn *B,sfixn invn,sfixn whichprime)
  {//-------------------------------------
    if(whichprime ==1)
#if FURER
    	      FURERPBPAS1::InvDFT_eff_p1(n,r,A,W,B,invn);
#else
    	      PBPAS1::InvDFT_eff_p1(n,r,A,W,B,invn);
#endif
    else if(whichprime == 2)
#if FURER
    	      FURERPBPAS2::InvDFT_eff_p2(n,r,A,W,B,invn);
#else
    	      PBPAS2::InvDFT_eff_p2(n,r,A,W,B,invn);
#endif
    else{
	    DFT_eff(n,r,A,W,pPtr,H,RevBidMap,B,whichprime);
	    //std::cout<<"InvDFT_eff before inv, A[0]: "<<A[0]<<"\n";
	    /*
	    std::cout<<"\n--- InvDFT_eff before div by n : \n";
	    for(int i=0; i<n; ++i)
	      std::cout<<A[i]<<", ";
	    std::cout<<std::endl;
	    */

	    for(int i=0; i<n; i++) {
			A[i]=MontMulMod_OPT2_AS_GENE(A[i], invn, pPtr);
		}
	}
  }
} //end of PBPAS
