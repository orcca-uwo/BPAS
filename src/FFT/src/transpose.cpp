#include "../../../include/FFT/src/transpose.h"

/* out-of-place transpose A[i0..i1][j0..j1] into B[j0..j1][i0..i1] */
//    int lda = nj;  // how many columns
//    int ldb = ni;  // how many rows
/* A has ni rows and nj columns, row-major layout */
void transpose_serial(sfixn *A, int lda, sfixn *B, int ldb,
	       int i0, int i1, int j0, int j1)
{
 tail:
     int di = i1 - i0, dj = j1 - j0;
     if (di >= dj && di > TRANSPOSETHRESHOLD) {
	  int im = (i0 + i1) / 2;
	  transpose_serial(A, lda, B, ldb, i0, im, j0, j1);
	  i0 = im; goto tail;
     } else if (dj > TRANSPOSETHRESHOLD) {
	  int jm = (j0 + j1) / 2;
	  transpose_serial(A, lda, B, ldb, i0, i1, j0, jm);
	  j0 = jm; goto tail;
     } else {
	  for (int i = i0; i < i1; ++i)
	       for (int j = j0; j < j1; ++j) 
		    B[j * ldb + i] = A[i * lda + j];
     }
}


void transpose(sfixn *A, int lda, sfixn *B, int ldb,
	       int i0, int i1, int j0, int j1)
{
 tail:
     int di = i1 - i0, dj = j1 - j0;
     if (di >= dj && di > TRANSPOSETHRESHOLD) {
	  int im = (i0 + i1) / 2;
	  cilk_spawn transpose(A, lda, B, ldb, i0, im, j0, j1);
	  i0 = im; goto tail;
     } else if (dj > TRANSPOSETHRESHOLD) {
	  int jm = (j0 + j1) / 2;
	  cilk_spawn transpose(A, lda, B, ldb, i0, i1, j0, jm);
	  j0 = jm; goto tail;
     } else {
	  for (int i = i0; i < i1; ++i)
	       for (int j = j0; j < j1; ++j) 
		    B[j * ldb + i] = A[i * lda + j];
     }
}

 
//invec is a matrix 2*4, SB is the output 
/*
  transpose_serial(invec,4,SB,2,0,2,0,4);
*/
 
