
#ifndef DUSP_NTL_SUPPORT_H
#define DUSP_NTL_SUPPORT_H

#if defined(WITH_NTL) && WITH_NTL

#include <gmpxx.h>
#include <iostream>
#include <time.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/pair_ZZ_pX_long.h>
#include <NTL/ZZ_pXFactoring.h>

#include "DUSP_Support.h"
#include "DUSP_Support_Test.h"
#include "DUSP_Support_Factoring.h"

#define MUL_NTLvsBPAS_THRESHOLD 10000
#define SQR_NTLvsBPAS_THRESHOLD 10000
#define EXP_NTLvsBPAS_THRESHOLD 10000
#define DIV_NTLvsBPAS_THRESHOLD 5000
#define GCD_NTLvsBPAS_THRESHOLD 5000
#define ECD_NTLvsBPAS_THRESHOLD 5000
#define FAC_NTLvsBPAS_THRESHOLD 1

// Convertors 
static inline void convertToNTLZZ_spX (const elem_t a, NTL::ZZ& f) { f = a; }
void convertToNTLPolyZZX_spX (const duspoly_t* a, NTL::ZZX& f);
void convertToNTLPolyZZpX_spX (const duspoly_t* a, NTL::ZZ_pX& f, const Prime_ptr* Pptr);

static inline void convertFromNTLZZ_spX (elem_t* a, const NTL::ZZ& f) { *a = NTL::conv<long>(f); }
void convertFromNTLZZX_spX (duspoly_t** a, const NTL::ZZX& f);
void convertFromNTLPolyZZpX_spX (duspoly_t** a, const NTL::ZZ_pX& f, const Prime_ptr* Pptr);

// Basic Operations 
	// Mul
void ntlMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr);
void ntlMulPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr);

void ntlSqrPolynomialInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr);

void ntlSSMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr);
void ntlSSMulPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr);

void ntlSSSqrPolynomialInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr);

	// Div
void ntlDivPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, duspoly_t** q, const Prime_ptr* Pptr);
void ntlDivPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr);

	// Rem
void ntlRemPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, const Prime_ptr* Pptr);
void ntlRemPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr);

	// Quo
void ntlQuoPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr);

	// GCD
void ntlGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, const Prime_ptr* Pptr);
void ntlExtGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, const Prime_ptr* Pptr);

// Factorial
void ntlFactoringInForm_spX (const duspoly_t* a, factors_t** facts, const Prime_ptr* Pptr);

// DUSP Factorization with using ntl-mul and ntl-div algorithms
// See the DUSP_Support_Factoring.h for the documentation
void squareFreeFactorizationInFormll1_spX (const duspoly_t* a, factsll_t** F, const Prime_ptr* Pptr);

void exponentiatePolynomialInForm1_spX (const duspoly_t* a, polysize_t n, duspoly_t** an, const Prime_ptr* Pptr);
void modExponentiatePolynomialInForm1_spX (const duspoly_t* a, polysize_t n, const duspoly_t* f, duspoly_t** an, const Prime_ptr* Pptr);

duspoly_t* expXModSubInForm1_spX (polysize_t q, const duspoly_t* f, duspoly_t** h, const Prime_ptr* Pptr);

void distinctDegFactorizationInFormll1_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr);

void modPoweringForEDF1_spX (const duspoly_t* a, polysize_t d, const duspoly_t* f, duspoly_t** mul_pows, const Prime_ptr* Pptr);

int equalDegSplittingInForm1_spX (const duspoly_t* f, polysize_t d, duspoly_t** g, const Prime_ptr* Pptr);
void equalDegFactorizationInFormll1_spX (const duspoly_t* f, polysize_t d, factsll_t** G, const Prime_ptr* Pptr);
void equalDegFactorsInFormll1_spX (factsll_t* D, factsll_t** G, const Prime_ptr* Pptr);

void modFactorizationInFormll1_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr);

int modFactorVerificationInFormll1_spX (const duspoly_t* f, factsll_t* G, const Prime_ptr* Pptr); 
 

#endif

#endif // DUSP_NTL_SUPPORT_H
