

#ifndef _SPARSE_INTERPOLATOR_H_
#define _SPARSE_INTERPOLATOR_H_

#include "RationalNumberPolynomial/SMQP_Support-AA.h"

#ifdef __cplusplus
extern "C" {
#endif



/**
 * MPQ Random number generator.
 */

typedef struct MPQRandomGenerator {
	gmp_randstate_t R_STATE;
	size_t numBits;
	size_t denBits;
} MPQRandomGenerator_t;

MPQRandomGenerator_t* initMPQRandomGenerator(unsigned long seed, size_t numBits, size_t denBits);

void MPQRandomGeneratorGetQ(MPQRandomGenerator_t* mpqRGen, mpq_t result);

void MPQRandomGeneratorGetZ(MPQRandomGenerator_t* mpqRGen, mpq_t result);

void freeMPQRandomGenerator(MPQRandomGenerator_t* mpqRGen);



/**
 * SMQP Black box function evaluation.
 */

typedef struct SMQPBlackBox {
	AltArr_t* aa;
} SMQPBlackBox_t;

SMQPBlackBox_t* initSMQPBlackBox(AltArr_t* aa);

void freeSMQPBlackBox(SMQPBlackBox_t* bb);

void evalSMQPBlackBox(SMQPBlackBox_t* bb, ratNum_t* points, int nvar, ratNum_t res);



/** 
 * Univariate dense interpolation. 
 * 
 * This is very slow. Better to use 
 * lagrange interpolation provided by univarInterpolate_AA (SMQP_Support).
 */
AltArr_t* univariateInterpolation_AA(int degreeBound, mpq_t* points, mpq_t* vals);



/** 
 * Multivariate interpolation.
 */
// AltArr_t* multivariateInterpolate_AA(SMQPBlackBox_t* bb, int* degreeBounds, int nvar);
AltArr_t* multivariateInterpolate_AA(SMQPBlackBox_t* bb, int* degreeBounds, int nvar, int T, double eps);

AltArr_t* multivariateInterpolateDeterministic_AA(SMQPBlackBox_t* bb, int* degreeBounds, int nvar, int T);

#ifdef __cplusplus
}
#endif


#endif //include guard


