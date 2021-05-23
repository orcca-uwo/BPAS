
#ifndef _DUSP_FFT_Support_H_
#define _DUSP_FFT_Support_H_

#include "../FFT/src/general_routine.h"
#include "../FFT/src/fft_furer1.h"
// #include "../FiniteFields/small_prime_field_fft.h"
// #include "../FiniteFields/gfpf_arithmetic.h"
#include "DUSP_Support.h"
#include "DUSP_Support_Test.h"
#include "math.h"

// Convert a DUSP Polynomial to C array of type sfixn*
static inline void convertDUSPToArray (duspoly_t* a, sfixn** coefs) {
    normalize_spX (&a);
    *coefs = (long*) a->elems;
}
// Convert a C array of type sfixn* to DUSP Polynomial
static inline void convertDUSPFromArray (duspoly_t** a, sfixn* coefs, polysize_t sz, Prime_ptr* Pptr) {
    setCoefsInForm_spX (a,(long long*) coefs, sz);    
}

/***********************************
 * Basic Arithmetic 
 * For details, see ./DUSP_Support.h 
***********************************/

/**
 * FFT-based Multiplication
 * @return the multiplication of a and b in c where
 */	
void fastMulPolynomialInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** c, Prime_ptr* Pptr);

void fastPSInversionInForm_FFT_spX (duspoly_t* a, duspoly_t** ai, polysize_t lgdg, Prime_ptr* Pptr);

void fastDivPolynomialInForm_wPSInv_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** r, duspoly_t** q, Prime_ptr* Pptr);

void mulMat4ToVecInForm_FFT_spX_inp (dusmat4_t* M, duspoly_t** u, duspoly_t** v, Prime_ptr* Pptr);

void mulMat4ToMat4InForm_FFT_spX_inp (dusmat4_t* A, dusmat4_t* B, dusmat4_t** C, Prime_ptr* Pptr);


/**
 * Yap Half-GCD algorithm with GCD and ExtGCD; Both are using Yap Half-GCD.
*/
void YapHalfGCDMatrixInForm_wPSInv_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, Prime_ptr* Pptr);
void GCDMatrixInForm_wHGCD_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, Prime_ptr* Pptr);
void ExtGCDInForm_wHGCD_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, Prime_ptr* Pptr);


/**
 * updated version of MCA HGCD algorithms
 */

void iterHalfGCDMatrixInForm_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, Prime_ptr* Pptr);   
void halfGCDMatrixInForm_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, Prime_ptr* Pptr);

void extHalfGCDMatrixInForm_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, Prime_ptr* Pptr);
void extHalfGCDMatrixInForm_FFT_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, Prime_ptr* Pptr);

void halfGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** ap, duspoly_t** bp, Prime_ptr* Pptr);
void halfGCDInForm_FFT_spX_inp (duspoly_t** a, duspoly_t** b, Prime_ptr* Pptr);

void GCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** g, Prime_ptr* Pptr);
void extGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, Prime_ptr* Pptr);

	
/** 
 * Sylvester Subresultant Chain Algorithm
 */

// void YapHGCD_BaseCaseInFrom_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, Prime_ptr* Pptr);
// void YapHGCD_TotalCaseInFrom_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, Prime_ptr* Pptr);

// void extYapHGCD_BaseCaseInFrom_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, Prime_ptr* Pptr);
// void extYapHGCD_TotalCaseInFrom_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, Prime_ptr* Pptr);
		
// void MCAHalfGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t mca_k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, polysize_t* h, Prime_ptr* Pptr);
// void MCAHGCD_BaseCaseInForm_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t k, polysize_t mca_k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, polysize_t* h, Prime_ptr* Pptr);
// void MCAHGCD_TotalCaseInFrom_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, Prime_ptr* Pptr);

// void quotientsToRemainders_spX (duspoly_t** Quo, duspoly_t*** Rem, polysize_t qlen, Prime_ptr* Pptr);
// void resultantYapHGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** res, Prime_ptr* Pptr);
// void subResultantYapHGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspolys_t** subres, polysize_t* sz, Prime_ptr* Pptr);

// void kSubResultantYapHGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t k, duspoly_t** ksubres, Prime_ptr* Pptr);
// void kSubResultantMCAHGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, polysize_t k, duspoly_t** ksubres, Prime_ptr* Pptr);



#endif
