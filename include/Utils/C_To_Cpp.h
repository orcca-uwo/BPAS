

#ifndef _C_2_CPP_H_
#define _C_2_CPP_H_


// #ifdef __cplusplus
// extern "C" {
// #endif

#if defined(WITH_MAPLE) && WITH_MAPLE
extern char* (*gcd_maple_string)(const char*, const char*);
#endif


#if defined(WITH_NTL) && WITH_NTL
#include "../IntegerPolynomial/SMZP_Support.h"
extern void (*factor_AAZ)(const AltArrZ_t* f_in, AltArrZ_t*** factors, int** exps, int* size, mpz_t content);
#endif




// #ifdef __cplusplus
// }
// #endif

#endif