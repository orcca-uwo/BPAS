
#ifndef _SMQP_SUPPORT_RECURSIVE_UNPACKED_AA_H_
#define _SMQP_SUPPORT_RECURSIVE_UNPACKED_AA_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ExponentPackingDefines.h"
#include "SMQP_Support_Unpacked.h"

#include <gmp.h>

typedef struct RecArrElem RecArrElem_t;
typedef struct RecArr RecArr_t;
typedef struct RecProdHeap_AA RecProdHeap_AA;
typedef struct RecProdHeapChain_AA RecProdHeapChain_AA;
typedef struct AltArrs AltArrs_t;


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
RecArr_t* convertToRecursiveArray_unpk(AltArr_t* aa);

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
AltArr_t* convertFromRecursiveArray_unpk(RecArr_t* poly, int nvar);


#ifdef __cplusplus
}
#endif


#endif

