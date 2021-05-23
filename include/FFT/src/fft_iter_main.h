


#ifndef _FFT_ITER_MAIN_
#define _FFT_ITER_MAIN_

#ifdef __cplusplus
extern "C" {
#endif

#include "malloc.h"
#include <string.h>

#include "../../FiniteFields/SmallPrimeField_Support.h"


/*
*n = 2^r
*e.g. n = 16, contents of rootsTable:
* 1, w, w^2, w^3, w^4, w^5, w^6, w^7, w^8, w^9, w^10, w^11, w^12, w^13, w^14, w^15
* 1, w^2, w^4, w^6, w^8, w^10, w^12, w^14
* 1, w^4, w^8, w^12
* 1, w^8
*/
void RootsTable_spf(int n, int r,
	long long int *T, const Prime_ptr* Pptr);


// Ti is the reversal of T except the first column
/*
*n = 2^r
*e.g. n = 16, contents of inverserootsTable:
*  m, m^2, m^3, m^4, m^5, m^6, m^7, m^8, m^9, m^10, m^11, m^12, m^13, m^14, m^15
*  m^2, m^4, m^6, m^8, m^10, m^12, m^14
*  m^4, m^8, m^12
*  m^8
*/
void InverseRootsTable(int n, long long int *T, long long int *Ti);


/*
* Compute roots table and inverse roots table for n = 2^r.
*/
static inline void RootsTable2_spf(int n, int r, long long int *T, long long int *Ti, const Prime_ptr* Pptr)
{//----------------------------------------------
	RootsTable_spf(n,r,T,Pptr);
	InverseRootsTable(n,T,Ti);
}


/**
 * Compute an DFT over the data A of length n,
 * with n = 2^r.
 * W is a pointer to second half or roots table (RootsTable_spf(n,r,Omega,Pptr); Omega += n;)
 * H is FFT threshold (1024).
 * B is temporary working space of size n*sizeof(long long int).
 *
 */
void DFT_eff(int n, int r,
	       long long int *A,
	       long long int *W,
	       const Prime_ptr* Pptr,
	       int H,
	       long long int *B);



void InvDFT_eff(int n, int r,
	  long long int *A,
	  long long int *W,
	  const Prime_ptr* Pptr,
	  int H,
	  long long int *B,
	  long long int invn);

#ifdef __cplusplus
}
#endif

#endif