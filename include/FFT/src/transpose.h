#ifndef _TFT_TRANSPOSE_H_
#define _TFT_TRANSPOSE_H_

#include "modpn.h"
#define TRANSPOSETHRESHOLD 16


void transpose_serial(sfixn *A, int lda, sfixn *B, int ldb,
	       int i0, int i1, int j0, int j1);

void transpose(sfixn *A, int lda, sfixn *B, int ldb,
	       int i0, int i1, int j0, int j1);

#endif
