

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


/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


