
#ifndef _SMZP_SUPPORT_RECURSIVE_UNPACKED_H_
#define _SMZP_SUPPORT_RECURSIVE_UNPACKED_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "../RationalNumberPolynomial/ExponentPackingDefines.h"
#include "SMZP_Support_Unpacked.h"

#include <gmp.h>

typedef struct RecArrZ RecArrZ_t;

/***********
 * Conversion to and from recursive view.
 ***********/

/**
 * Convert a polynomial to a recursively viewed polynomial using the variable
 * of highest precedence (index 0) as the main variable.
 *
 * Note this conversion is done INPLACE. And therefore the input Node* is 
 * invalidated by this operation.
 *
 * returns teh recursively viewed polynomial.
 */
RecArrZ_t* convertToRecursiveArrayZ_unpk(AltArrZ_t* aa);

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done INPLACE and therefore the original recursively
 * viewed polynomial is then invalidated.
 * 
 * nvar is the number of variables the returning AltArr_t should have. That is, 
 * the number of variables in the RecArr's coefficients + 1.
 *
 * returns the distributed polynomial.
 */
AltArrZ_t* convertFromRecursiveArrayZ_unpk(RecArrZ_t* poly, int nvar);



#ifdef __cplusplus
}
#endif


#endif