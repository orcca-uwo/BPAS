/**
	Implementation of 'MulKS.h'

	@author Farnam Mansouri
*/

#include "../../../include/IntegerPolynomial/Multiplication/MulKS.h"

void MulKS::mulBigIntSer(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){
	
	/*int digitCount = ((a->getSize() < b->getSize()) ?  ceil(a->getSize() / a->getRepresentationBase()) :
                                                        ceil(b->getSize() / b->getRepresentationBase()))
                                + a->getCoefficientDigits() + b->getCoefficientDigits();

        mpz_class aInteger = a->getBigIntegerRepresentation(digitCount);
        mpz_class bInteger = b->getBigIntegerRepresentation(digitCount);

        mpz_t aMpz = aInteger.get_mpz_t();
        mpz_t bMpz = bInteger.get_mpz_t();
        mp_int aSize = mpz_getSize()(aMpz);
        mp_int bSize = mpz_getSize()(bMpz);
        mp_limb_t * aLimbs = new mp_limb_t[aSize];
        mp_limb_t * bLimbs = new mp_limb_t[bSize];

        for(int i = 0; i < aSize; i++){
                lA[i] = mpz_getlimbn(aMpz, i);
        }

        for(int i = 0; i < bSize; i++){
                lB[i] = mpz_getlimbn(bMpz, i);
        }

        mp_limb_t *cLimbs = new mp_limb_t[aSize + bSize];

        mpn_mul(cLims, aLimbs, aSize, bLimbs,  bSize);

        //mpz_class result(climbs);
        */

	a->fixDegree();
	b->fixDegree();


        setBitCounts(a, b, c);

        int digitCount = MAX_DIGITS_P(a, b) + 1;

        mpz_class aInteger = a->getBigIntegerSigned(digitCount);
        mpz_class bInteger = b->getBigIntegerSigned(digitCount);

        if(a->isNegativeLC()) aInteger = -aInteger;
        if(b->isNegativeLC()) bInteger = -bInteger;

        mpz_class result = aInteger * bInteger;

        int range = RESULT_DEGREE_P(a, b);

        c->convertFromBigInteger(result, startIndex, startIndex + range, digitCount);

}


void MulKS::mulBigIntSerReciprocalUnsigned(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){
        int digitCount = MAX_DIGITS_FAST_P(a, b);

        mpz_class c1Int = generateFirstInteger(a, b, digitCount);
        mpz_class c2Int = generateSecondInteger(a, b, digitCount);

        c->extractCoeffsicients(c1Int, c2Int, startIndex, startIndex + RESULT_DEGREE_P(a, b), digitCount);
}


mpz_class MulKS::generateFirstInteger(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, int digitCount){
        mpz_class aInt = a->getBigIntegerUnsigned(digitCount);
        mpz_class bInt = b->getBigIntegerUnsigned(digitCount);
        mpz_class cInt = aInt * bInt;
        return cInt;
}


mpz_class MulKS::generateSecondInteger(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, int digitCount){
        mpz_class aInt = a->getReverseBigIntegerUnsigned(digitCount);
        mpz_class bInt = b->getReverseBigIntegerUnsigned(digitCount);
        mpz_class cInt = aInt * bInt;
        return cInt;
}
