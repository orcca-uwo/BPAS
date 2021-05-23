
/*****
 * Supporting methods for sparse multivariate rational polynomials written in
 * pure C.
 *
 * Throughout this it is assumed that polynomials are compatiable. They must have
 * the same number of variables and the same variable ordering.
 *
 * Lexicographical ordering is used for term ordering throughout.
 *
 * Functions assume that the number of variables is at least 1. But, in the most
 * basic case, a NULL exponent vector encodes nvar = 0. Similar to how Node* = NULL
 * encodes the zero polynomial.
 *
 * NOTE: functions herein are still experimental, may not be bug-free,
 *       and may change in future releases.
 *
 *****/

#ifndef _SMZP_SUPPORT_FACTORING_H_
#define _SMZP_SUPPORT_FACTORING_H_

//#ifdef __cplusplus
//extern "C" {
//#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>

/*#include "../RationalNumberPolynomial/ExponentPackingDefines.h"*/
/*#include "../RationalNumberPolynomial/SMQP_Support-AA.h"*/
#include "SMZP_Support_Unpacked.h"
#include "IntegerPolynomial/DUZP_Support.h"
#include "IntegerPolynomial/DUZP_Support_Factoring.hpp"

namespace SMZP {

	namespace Factoring {

		int reconstructLeadingCoefficientBivariate(DUZP_t*** lcs, AltArrZ_t* f, const DUZP::Factoring::vec_DUZP_long_t* lc_factors, const mpz_t* F_tilde, DUZP::Factoring::vec_DUZP_long_t* U0_factors, const mpz_t delta_in, const mpz_t* values);

		void chooseRandomEvaluationPoint(mpz_t** values, AltArrZ_t** U0, mpz_t delta, const AltArrZ_t* lc, const AltArrZ_t* f_in, mpz_t bound);

		int chooseSpecialIntegers(mpz_t* d, const mpz_t Omega, const mpz_t* F_tilde, long k, const mpz_t delta);

		void factor(const AltArrZ_t* f_in, AltArrZ_t*** factors, int** exps, int* size, mpz_t content);

		void factorBivariate(const AltArrZ_t* f, AltArrZ_t*** factors, int** exps, int* size, mpz_t content);

	}

}

//#ifdef __cplusplus
//}
//#endif

#endif
