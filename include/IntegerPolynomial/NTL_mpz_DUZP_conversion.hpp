#ifndef NTL_MPZ_CONVERSION
#define NTL_MPZ_CONVERSION

#include <gmpxx.h>
#include "NTL/ZZ.h"
#include "NTL/ZZX.h"
#include <gmp.h>
#include "../IntegerPolynomial/DUZP_Support.h"

/*
 * Set an NTL integer (ZZ) using an mpz_class of GMP
 *
 */
void NTLZZ_set_mpz_class (NTL::ZZ& g, const mpz_class& f);

/*
 * Set an NTL integer (ZZ) using an mpz_t of GMP
 *
 */
void NTLZZ_set_mpz_t (NTL::ZZ& g, const mpz_t f1);

/*
 * Set an mpz_class of GMP using an integer (ZZ)
 *
 */
void mpz_class_set_NTLZZ (mpz_class& f, const NTL::ZZ& g);

/*
 * Set an mpz_t of GMP using an integer (ZZ)
 *
 */
void mpz_t_set_NTLZZ (mpz_t& f, const NTL::ZZ& g);

/*
 * Set an NTL polynomial (ZZX) using a DUZP_t
 *
 */
void NTLZZX_set_DUZP_t (NTL::ZZX& g, const DUZP_t* f);

/*
 * Set a DUZP_t using an NTL polynomial (ZZX)
 *
 */
void DUZP_t_set_NTLZZX (DUZP_t** f, const NTL::ZZX& g);

#endif
