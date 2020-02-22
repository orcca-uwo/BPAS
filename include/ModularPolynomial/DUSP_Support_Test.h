
#ifndef _DUSP_SUPPORT_TEST_H_
#define _DUSP_SUPPORT_TEST_H_

#ifdef __cplusplus
    extern "C" {
#endif

#include "./DUSP_Support.h"
#include "./DUSP_Support_Factoring.h"
    	
// Print polynomial
void printPolynomial_spX (const duspoly_t* a, const Prime_ptr* Pptr);
void printPolynomialOutForm_spX (const duspoly_t* a, const Prime_ptr* Pptr);
    
// print Matrix 2*2 of polynomials 
void printMat4_spX (dusmat4_t* A, const Prime_ptr* Pptr);
void printMat4OutForm_spX (dusmat4_t* A, const Prime_ptr* Pptr);

void printFactors_spX (factors_t* f, const Prime_ptr* Pptr);
void printFactsll_spX (factsll_t* f, const Prime_ptr* Pptr);

void printFactorsOutForm_spX (factors_t* f, const Prime_ptr* Pptr);
void printFactsllOutForm_spX (factsll_t* f, const Prime_ptr* Pptr);

void purePrintFactsll_spX (factsll_t* f, const Prime_ptr* Pptr);
void purePrintFactsllOutForm_spX (factsll_t* f, const Prime_ptr* Pptr);

void purePrintFactorsOutForm_spX (factors_t* f, const Prime_ptr* Pptr);

// Random Polynomial
duspoly_t* randomFullyDensePolynomial_spX (polysize_t sz, const Prime_ptr* Pptr);

int DivisionTest (const duspoly_t* a, const duspoly_t* b, const duspoly_t* r, const duspoly_t* q, const Prime_ptr* Pptr);

int ExtEuclideanTest (const duspoly_t* a, const duspoly_t* b, const duspoly_t* s, const duspoly_t* t, const duspoly_t* g, const Prime_ptr* Pptr);
int ExtEuclideanTestInForm (const duspoly_t* a, const duspoly_t* b, const duspoly_t* s, const duspoly_t* t, const duspoly_t* g, const Prime_ptr* Pptr);

int SFFVerificationInForm_spX (const duspoly_t* a, factors_t* f, const Prime_ptr* Pptr);

		
#ifdef __cplusplus
}
#endif

#endif
