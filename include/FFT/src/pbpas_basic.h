// basic-routine.h

#ifndef _pbpas_basic_h_
#define _pbpas_basic_h_

#include "general_routine.h"

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
    void DMPMul(sfixn p,  int N,    
		sfixn *C, int *dc, 
		sfixn *A, int *da,
		sfixn *B, int *db,
		int wp);

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
		  int wp);

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
		      int H, int *RevBidMap, int wp);

    void Interpolation2D(sfixn *res, int es1, int es2,
			 int dims1, int dims2, 
			 int ls1, int ls2,
			 sfixn *invRT1, sfixn *invRT2,
			 MONTP_OPT2_AS_GENE *pPtr, 
			 int H, int *RevBidMap, int wp);
}

#endif

