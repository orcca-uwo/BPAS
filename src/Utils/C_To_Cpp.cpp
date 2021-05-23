

#if defined(WITH_NTL) && WITH_NTL
#include <gmpxx.h>
#endif

#include "Utils/C_To_Cpp.h"



#if defined(WITH_MAPLE) && WITH_MAPLE
#include "MapleInterface/MapleInterfaceStream.hpp"
char* (*gcd_maple_string)(const char*, const char*) = &MapleInterfaceStream::gcd_string;
#endif

#if defined(WITH_NTL) && WITH_NTL
#include "IntegerPolynomial/SMZP_Support_Factoring.hpp"
void (*factor_AAZ)(const AltArrZ_t*, AltArrZ_t***, int**, int*, mpz_t) = SMZP::Factoring::factor;
#endif
