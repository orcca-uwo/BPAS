/**
	Implementation of 'MulDnC.h'

	@author Farnam Mansouri
*/

#include "../../../include/IntegerPolynomial/Multiplication/MulDnC.h"

MulDnC::MulDnC(){
        // Default Threshold
        BASE_MUL = 512;
}


void MulDnC::mulDnCBigIntPar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){
        if(a->getSize() <= BASE_MUL && b->getSize() <= BASE_MUL){
		MulKS m;
                m.mulBigIntSer(a, b, c, startIndex);
        } else {

                UnivariateIntegerPolynomial a0(a->getSize() / 2, a->getCoefficientDigits(), a->getCoefficients());
                UnivariateIntegerPolynomial a1(a->getSize() - a0.getSize(), a->getCoefficientDigits(), a->getPointer(a0.getSize()));
                UnivariateIntegerPolynomial b0(b->getSize() / 2, b->getCoefficientDigits(), b->getCoefficients());
                UnivariateIntegerPolynomial b1(b->getSize() - b0.getSize(), b->getCoefficientDigits(), b->getPointer(b0.getSize()));

                cilk_spawn mulDnCBigIntPar(&a0, &b0, c, startIndex);
                cilk_spawn mulDnCBigIntPar(&a1, &b1, c,  startIndex + a0.getSize() + b0.getSize());

		/*mulDnRBigIntPar(&a0, &b1, c, startIndex + b0.getSize());
                mulDnRBigIntPar(&a1, &b0, c, startIndex + a0.getSize());*/

		/*UnivariateIntegerPolynomial tmp(RESULT_DEGREE(a1, b0), MAX_DIGITS(a1, b0));

                cilk_spawn mulDnRBigIntPar(&a0, &b1, c, startIndex + b0.getSize());
                mulDnRBigIntPar(&a1, &b0, &tmp, 0);
                cilk_sync;

                int stId = startIndex + a0.getSize();
                cilk_for(int i = 0; i < a1.getSize() + b0.getSize() - 1; ++i)
                        c->setCoefficient(stId + i, c.getCoefficient(stId + i) + tmp.getCoefficient(i));
                tmp.freeHeap();*/

                UnivariateIntegerPolynomial tmp0(RESULT_DEGREE(a0, b1), MAX_DIGITS(a0, b1));
                UnivariateIntegerPolynomial tmp1(RESULT_DEGREE(a1, b0), MAX_DIGITS(a1, b0));

                cilk_spawn mulDnCBigIntPar(&a0, &b1, &tmp0, 0);//&tmp0, startIndex);
                mulDnCBigIntPar(&a1, &b0, &tmp1, 0);//startIndex + a0.getSize() + b1.getSize());
                cilk_sync;

                int stId = startIndex + b0.getSize();
                for(int i = 0; i < a0.getSize() + b1.getSize() - 1; ++i)
                        c->setCoefficient(stId + i, c->getCoefficient(stId + i) + tmp0.getCoefficient(i));

                stId = startIndex + a0.getSize();
                for(int i = 0; i < a1.getSize() + b0.getSize() - 1; ++i)
                        c->setCoefficient(stId + i, c->getCoefficient(stId + i) + tmp0.getCoefficient(i));

                tmp0.freeHeap(); tmp1.freeHeap();

        }
}

UnivariateIntegerPolynomial MulDnC::mulDnC(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
        UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));

	setBitCounts(a, b);

        int digitCount = MAX_DIGITS_P(a, b) + 1;

        //DIVIDING
        UnivariateIntegerPolynomial a0Poly(a->getSize() / 4, a->getCoefficientDigits(), a->getCoefficients());
        UnivariateIntegerPolynomial a1Poly(a->getSize() / 4, a->getCoefficientDigits(),  a->getPointer(a0Poly.getSize()));
        UnivariateIntegerPolynomial a2Poly(a->getSize() / 4, a->getCoefficientDigits(),  a->getPointer(a0Poly.getSize() + a1Poly.getSize()));
        UnivariateIntegerPolynomial a3Poly(a->getSize() - a0Poly.getSize() - a1Poly.getSize() - a2Poly.getSize(),
                        a->getCoefficientDigits(),  a->getPointer(a0Poly.getSize() + a1Poly.getSize() + a2Poly.getSize()));

        UnivariateIntegerPolynomial b0Poly(b->getSize() / 4, b->getCoefficientDigits(), b->getCoefficients());
        UnivariateIntegerPolynomial b1Poly(b->getSize() / 4, b->getCoefficientDigits(), b->getPointer(b0Poly.getSize()));
        UnivariateIntegerPolynomial b2Poly(b->getSize() / 4, b->getCoefficientDigits(), b->getPointer(b0Poly.getSize() + b1Poly.getSize()));
        UnivariateIntegerPolynomial b3Poly(b->getSize() - b0Poly.getSize() - b1Poly.getSize() - b2Poly.getSize(),
                        b->getCoefficientDigits(), b->getPointer(b0Poly.getSize() + b1Poly.getSize() + b2Poly.getSize()));


        //CONVERSION
#if __GNUC__ == 4
        mpz_class a0 = cilk_spawn
                        a0Poly.getBigIntegerSigned(digitCount);
        mpz_class a1 = cilk_spawn
                        a1Poly.getBigIntegerSigned(digitCount);
        mpz_class a2 = cilk_spawn
                        a2Poly.getBigIntegerSigned(digitCount);
        mpz_class a3 = cilk_spawn
                        a3Poly.getBigIntegerSigned(digitCount);
        mpz_class b0 = cilk_spawn
                        b0Poly.getBigIntegerSigned(digitCount);
        mpz_class b1 = cilk_spawn
                        b1Poly.getBigIntegerSigned(digitCount);
        mpz_class b2 = cilk_spawn
                        b2Poly.getBigIntegerSigned(digitCount);
        mpz_class b3 =
                        b3Poly.getBigIntegerSigned(digitCount);
        cilk_sync;
#else
	mpz_class a0, a1, a2, a3, b0, b1, b2, b3;
  	cilk_spawn a0Poly.getBigIntegerSigned(&a0, digitCount);
        cilk_spawn a1Poly.getBigIntegerSigned(&a1, digitCount);
        cilk_spawn a2Poly.getBigIntegerSigned(&a2, digitCount);
        cilk_spawn a3Poly.getBigIntegerSigned(&a3, digitCount);
        cilk_spawn b0Poly.getBigIntegerSigned(&b0, digitCount);
        cilk_spawn b1Poly.getBigIntegerSigned(&b1, digitCount);
        cilk_spawn b2Poly.getBigIntegerSigned(&b2, digitCount);
        b3Poly.getBigIntegerSigned(&b3, digitCount);
        cilk_sync;
#endif

        if(a0Poly.isNegativeLC()) a0 = -a0;
        if(a1Poly.isNegativeLC()) a1 = -a1;
        if(a2Poly.isNegativeLC()) a2 = -a2;
        if(a3Poly.isNegativeLC()) a3 = -a3;
        if(b0Poly.isNegativeLC()) b0 = -b0;
        if(b1Poly.isNegativeLC()) b1 = -b1;
        if(b2Poly.isNegativeLC()) b2 = -b2;
        if(b3Poly.isNegativeLC()) b3 = -b3;

	//EVALUATION & MULTIPLYING
        mpz_t * w = new mpz_t[16];
        for(int i = 0; i < 16; i++)     mpz_init(w[i]);
        cilk_spawn
                mpz_mul(w[0], a0.get_mpz_t(), b0.get_mpz_t());
        cilk_spawn
                mpz_mul(w[1], a0.get_mpz_t(), b1.get_mpz_t());
        cilk_spawn
                mpz_mul(w[2], a0.get_mpz_t(), b2.get_mpz_t());
        cilk_spawn
                mpz_mul(w[3], a0.get_mpz_t(), b3.get_mpz_t());
        cilk_spawn
                mpz_mul(w[4], a1.get_mpz_t(), b0.get_mpz_t());
        cilk_spawn
                mpz_mul(w[5], a1.get_mpz_t(), b1.get_mpz_t());
        cilk_spawn
                mpz_mul(w[6], a1.get_mpz_t(), b2.get_mpz_t());
        cilk_spawn
                mpz_mul(w[7], a1.get_mpz_t(), b3.get_mpz_t());
        cilk_spawn
                mpz_mul(w[8], a2.get_mpz_t(), b0.get_mpz_t());
        cilk_spawn
                mpz_mul(w[9], a2.get_mpz_t(), b1.get_mpz_t());
        cilk_spawn
                mpz_mul(w[10], a2.get_mpz_t(), b2.get_mpz_t());
        cilk_spawn
                mpz_mul(w[11], a2.get_mpz_t(), b3.get_mpz_t());
        cilk_spawn
                mpz_mul(w[12], a3.get_mpz_t(), b0.get_mpz_t());
        cilk_spawn
                mpz_mul(w[13], a3.get_mpz_t(), b1.get_mpz_t());
        cilk_spawn
                mpz_mul(w[14], a3.get_mpz_t(), b2.get_mpz_t());

                mpz_mul(w[15], a3.get_mpz_t(), b3.get_mpz_t());
        cilk_sync;

        mpz_class w0(w[0]); mpz_class w1(w[1]); mpz_class w2(w[2]); mpz_class w3(w[3]);
        mpz_class w4(w[4]); mpz_class w5(w[5]); mpz_class w6(w[6]); mpz_class w7(w[7]);
        mpz_class w8(w[8]); mpz_class w9(w[9]); mpz_class w10(w[10]); mpz_class w11(w[11]);
        mpz_class w12(w[12]); mpz_class w13(w[13]); mpz_class w14(w[14]); mpz_class w15(w[15]);

        //INTERPOLATION
        mpz_class * co = new mpz_class[7];
        co[0] = w0;
        co[1] = w1 + w4;
        co[2] = w2 + w5 + w8;
        co[3] = w3 + w6 + w9 + w12;
        co[4] = w7 + w10 + w13;
        co[5] = w11 + w14;
        co[6] = w15;

	//CONVERSION
        int range = RESULT_DEGREE(a0Poly, b0Poly);
        int rangeHalf = range / 2;
	
	UnivariateIntegerPolynomial cTmp(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));

        int start = range - rangeHalf;

        cilk_spawn
                c.convertFromBigInteger(co[0], 0, range , digitCount);
        cilk_spawn
                c.convertFromBigInteger(co[2], 2 * start, 2 * start + range, digitCount);
        cilk_spawn
                c.convertFromBigInteger(co[4], 4 * start, 4 * start + range, digitCount);
        cilk_spawn
                c.convertFromBigInteger(co[6], 6 * start, 6 * start + range, digitCount);
        cilk_spawn
                cTmp.convertFromBigInteger(co[1], start, start + range, digitCount);
        cilk_spawn
                cTmp.convertFromBigInteger(co[3], 3 * start, 3 * start + range, digitCount);

                cTmp.convertFromBigInteger(co[5], 5 * start, 5 * start + range, digitCount);
        cilk_sync;

	delete[] co; delete[] w;

        int end = start * 5 + range;
        for(int i = start; i < end; ++i)
                c.setCoefficient(i, c.getCoefficient(i) + cTmp.getCoefficient(i));

        cTmp.freeHeap();

        return c;

}


void MulDnC::mulDnCBigIntParStatic(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, UnivariateIntegerPolynomial* tmp, int startIndex){
        if(a->getSize() <= BASE_MUL && b->getSize() <= BASE_MUL){
		MulKS m;
                m.mulBigIntSer(a, b, c, startIndex);
        } else {

                UnivariateIntegerPolynomial a0(a->getSize() / 2, a->getCoefficientDigits(), a->getCoefficients());
                UnivariateIntegerPolynomial a1(a->getSize() - a0.getSize(), a->getCoefficientDigits(), a->getPointer(a0.getSize()));
                UnivariateIntegerPolynomial b0(b->getSize() / 2, b->getCoefficientDigits(), b->getCoefficients());
                UnivariateIntegerPolynomial b1(b->getSize() - b0.getSize(), b->getCoefficientDigits(), b->getPointer(b0.getSize()));

                cilk_spawn mulDnCBigIntParStatic(&a0, &b0, c, tmp, startIndex);
                cilk_spawn mulDnCBigIntParStatic(&a1, &b1, c, tmp, startIndex + a0.getSize() + b0.getSize());

                cilk_sync;

                cilk_spawn mulDnCBigIntParStatic(&a0, &b1, c, tmp, startIndex + b0.getSize());
                mulDnCBigIntParStatic(&a1, &b0, tmp, c, startIndex + a0.getSize());
                cilk_sync;
        }
}


void MulDnC::mulDnCBigIntParStatic(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, UnivariateIntegerPolynomial* tmp0, UnivariateIntegerPolynomial* tmp1, int startIndex){
        if(a->getSize() <= BASE_MUL && b->getSize() <= BASE_MUL){
		MulKS m;		
                m.mulBigIntSer(a, b, c, startIndex);
        } else {

                UnivariateIntegerPolynomial a0(a->getSize() / 2, a->getCoefficientDigits(), a->getCoefficients());
                UnivariateIntegerPolynomial a1(a->getSize() - a0.getSize(), a->getCoefficientDigits(), a->getPointer(a0.getSize()));
                UnivariateIntegerPolynomial b0(b->getSize() / 2, b->getCoefficientDigits(), b->getCoefficients());
                UnivariateIntegerPolynomial b1(b->getSize() - b0.getSize(), b->getCoefficientDigits(), b->getPointer(b0.getSize()));

                cilk_spawn mulDnCBigIntParStatic(&a0, &b0, c, tmp0, tmp1, startIndex);
                cilk_spawn mulDnCBigIntParStatic(&a1, &b1, c, tmp0, tmp1,  startIndex + a0.getSize() + b0.getSize());

                cilk_spawn mulDnCBigIntParStatic(&a0, &b1, tmp0, tmp1, c, 0);
                mulDnCBigIntParStatic(&a1, &b0, tmp1, tmp0, c, 0);
                cilk_sync;
        }
}


void MulDnC::mulIterBigIntPar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){
	MulKS m;
        if(a->getSize() <= BASE_MUL && b->getSize() <= BASE_MUL){
                m.mulBigIntSer(a, b, c, startIndex);
        } else {

                int polyPairs = ( MAX(a->getSize(), b->getSize()) / (2 * BASE_MUL) );

                int aLeftDeg, aRightDeg, bLeftDeg, bRightDeg;
                //aLeftDeg = a->getSize() / (1 << (polyPairs + 1));
                //aRightDeg = (a->getSize() >> polyPairs) - aLeftDeg;
                //bLeftDeg = b->getSize() / (1 << (polyPairs + 1));
                //bRightDeg = (b->getSize() >> polyPairs) - bLeftDeg;
                aLeftDeg = BASE_MUL; aRightDeg = BASE_MUL; bLeftDeg = BASE_MUL; bRightDeg = BASE_MUL;

                int opSize = aLeftDeg + aRightDeg + bLeftDeg + bRightDeg;

                UnivariateIntegerPolynomial a0(aLeftDeg, a->getCoefficientDigits(), a->getCoefficients());
                UnivariateIntegerPolynomial a1(aRightDeg, a->getCoefficientDigits(), a->getCoefficients());
                UnivariateIntegerPolynomial b0(bLeftDeg, b->getCoefficientDigits(), b->getCoefficients());
                UnivariateIntegerPolynomial b1(bRightDeg, b->getCoefficientDigits(), b->getCoefficients());

                UnivariateIntegerPolynomial tmp(opSize, c->getCoefficientDigits());

                int aStartIndex, bStartIndex, cStartIndex;
                for(int i = 0; i < polyPairs; ++i){
                        for(int j = 0; j < polyPairs; ++j){

                                aStartIndex = i * (aLeftDeg + aRightDeg);
                                bStartIndex = j * (bLeftDeg + bRightDeg);
                                cStartIndex = aStartIndex + bStartIndex;

                                a0.setCoefficients(a->getPointer(aStartIndex));
                                a1.setCoefficients(a->getPointer(aStartIndex + aLeftDeg));
                                b0.setCoefficients(b->getPointer(bStartIndex));
                                b1.setCoefficients(b->getPointer(bStartIndex + bLeftDeg));

                                tmp.setToZero();

                                cilk_spawn m.mulBigIntSer(&a0, &b0, c, cStartIndex);
                                cilk_spawn m.mulBigIntSer(&a1, &b1, c, cStartIndex + aLeftDeg + bLeftDeg);
                                cilk_spawn m.mulBigIntSer(&a0, &b1, &tmp, 0);
                                m.mulBigIntSer(&a1, &b0, &tmp, aLeftDeg + bRightDeg);
                                cilk_sync;

                                int s = cStartIndex + bLeftDeg; int e = s + aLeftDeg + bRightDeg;
                                cilk_for(int k = s; k < e; ++k)
                                        c->setCoefficient(k, c->getCoefficient(k) + tmp.getCoefficient(k - s));

                                s = cStartIndex + aLeftDeg; e = s + bLeftDeg + aRightDeg;
                                int tmpIndex = aLeftDeg + bRightDeg - s;
                                cilk_for(int k = s; k < e; ++k)
					c->setCoefficient(k, c->getCoefficient(k) + tmp.getCoefficient(tmpIndex + k));
                        }

                }

                tmp.freeHeap();

        }
}


UnivariateIntegerPolynomial MulDnC::mulDnCBigIntParallelStatic1(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
        int size = RESULT_DEGREE_P(a, b);
        UnivariateIntegerPolynomial c(size, MAX_DIGITS_P(a, b));
        setBitCounts(a, b, &c);
        UnivariateIntegerPolynomial tmp(size, MAX_DIGITS_P(a, b));
        tmp.setRepresentationBase(a->getRepresentationBase());
        mulDnCBigIntParStatic(a, b, &c, &tmp, 0);

        //int a0 = a->getSize() / 2; int b0 = b->getSize() / 2;
        //int a1 = a->getSize() - a0; int b1 = b->getSize() - b0;

        for(int i = 0; i < size; ++i)
                c.setCoefficient(i, c.getCoefficient(i) + tmp.getCoefficient(i));

        tmp.freeHeap();

        return c;
}


UnivariateIntegerPolynomial MulDnC::mulDnCBigIntParallelStatic2(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
        int size = RESULT_DEGREE_P(a, b);
        UnivariateIntegerPolynomial c(size, MAX_DIGITS_P(a, b));
        setBitCounts(a, b, &c);
        UnivariateIntegerPolynomial tmp0(size/2, MAX_DIGITS_P(a, b));
        UnivariateIntegerPolynomial tmp1(size/2, MAX_DIGITS_P(a, b));
        tmp0.setRepresentationBase(a->getRepresentationBase());
        tmp1.setRepresentationBase(a->getRepresentationBase());
        mulDnCBigIntParStatic(a, b, &c, &tmp0, &tmp1, 0);

        int a0 = a->getSize() / 2; int b0 = b->getSize() / 2;
        int a1 = a->getSize() - a0; int b1 = b->getSize() - b0;

        for(int i = 0; i < a0 + b1 - 1; ++i)
                c.setCoefficient(b0 + i, c.getCoefficient(b0 + i) + tmp0.getCoefficient(i));

        for(int i = 0; i < a1 + b0 - 1; ++i)
		c.setCoefficient(a0 + i, c.getCoefficient(a0 + i) + tmp1.getCoefficient(i));

        tmp0.freeHeap(); tmp1.freeHeap();

        return c;
}


UnivariateIntegerPolynomial MulDnC::mulSignedDecomposed(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){

        UnivariateIntegerPolynomial aPlus(a->getSize(), a->getCoefficientDigits());
        UnivariateIntegerPolynomial bPlus(b->getSize(), b->getCoefficientDigits());
        UnivariateIntegerPolynomial aMinus(a->getSize(), a->getCoefficientDigits());
        UnivariateIntegerPolynomial bMinus(b->getSize(), b->getCoefficientDigits());
        aPlus.setToZero();      bPlus.setToZero();
        aMinus.setToZero();     bMinus.setToZero();

        int i, aPlusMax = 0, aMinusMax = 0, bPlusMax = 0, bMinusMax = 0, plusTmp, minusTmp;
        for(i = 0; i < a->getSize(); ++i){
                if(a->getCoefficient(i) >= 0){
                        aPlus.setCoefficient(i, a->getCoefficient(i));
                        plusTmp = mpz_sizeinbase(aPlus.getCoefficient(i).get_mpz_t(), 2);
                        if(plusTmp > aPlusMax)          aPlusMax = plusTmp;
                } else{
                        aMinus.setCoefficient(i, 0 - a->getCoefficient(i));
                        minusTmp = mpz_sizeinbase(aMinus.getCoefficient(i).get_mpz_t(), 2);
                        if(minusTmp > aMinusMax)        aMinusMax = minusTmp;
                }
        }

        for(i = 0; i < b->getSize(); ++i){
                if(b->getCoefficient(i) >= 0){
                        bPlus.setCoefficient(i, b->getCoefficient(i));
                        plusTmp = mpz_sizeinbase(bPlus.getCoefficient(i).get_mpz_t(), 2);
                        if(plusTmp > bPlusMax)          bPlusMax = plusTmp;
                } else{
                        bMinus.setCoefficient(i, 0 - b->getCoefficient(i));
                        minusTmp = mpz_sizeinbase(bMinus.getCoefficient(i).get_mpz_t(), 2);
                        if(minusTmp > bMinusMax)        bMinusMax = minusTmp;
                }
        }

        UnivariateIntegerPolynomial cPlusPlus(RESULT_DEGREE(aPlus, bPlus), MAX_DIGITS(aPlus, bPlus));
        UnivariateIntegerPolynomial cPlusMinus(RESULT_DEGREE(aPlus, bMinus), MAX_DIGITS(aPlus, bMinus));
        UnivariateIntegerPolynomial cMinusPlus(RESULT_DEGREE(aMinus, bPlus), MAX_DIGITS(aMinus, bPlus));
        UnivariateIntegerPolynomial cMinusMinus(RESULT_DEGREE(aMinus, bMinus), MAX_DIGITS(aMinus, bMinus));

        if((aPlusMax > 0 && aPlusMax < 32) || (aMinusMax > 0 && aMinusMax < 32) ||
                (bPlusMax > 0 && bPlusMax < 32) || (bMinusMax > 0 && bMinusMax < 32)){
                aPlus.setBitPackage(1); aMinus.setBitPackage(1);
                bPlus.setBitPackage(1); bMinus.setBitPackage(1);
                cPlusPlus.setBitPackage(1); cPlusMinus.setBitPackage(1);
                cMinusPlus.setBitPackage(1); cMinusMinus.setBitPackage(1);
                aPlus.setRepresentationBase(2); aMinus.setRepresentationBase(2);
                bPlus.setRepresentationBase(2); bMinus.setRepresentationBase(2);
                cPlusPlus.setRepresentationBase(2); cPlusMinus.setRepresentationBase(2);
                cMinusPlus.setRepresentationBase(2); cMinusMinus.setRepresentationBase(2);
        }

        aPlus.setCoefficientDigits(aPlusMax / aPlus.getBitPackage());
        aMinus.setCoefficientDigits(aMinusMax / aMinus.getBitPackage());
        bPlus.setCoefficientDigits(bPlusMax / bPlus.getBitPackage());
        bMinus.setCoefficientDigits(bMinusMax / bMinus.getBitPackage());

        //BASE_MUL = MAX(a->getSize(), b->getSize()) / 4;

	MulKS m;
        if(aPlus.getCoefficientDigits() > 0 && bPlus.getCoefficientDigits() > 0)
                cilk_spawn m.mulBigIntSer(&aPlus, &bPlus, &cPlusPlus, 0);
        if(aPlus.getCoefficientDigits() > 0 && bMinus.getCoefficientDigits() > 0)
                cilk_spawn m.mulBigIntSer(&aPlus, &bMinus, &cPlusMinus, 0);
        if(aMinus.getCoefficientDigits() > 0 && bPlus.getCoefficientDigits() > 0)
                cilk_spawn m.mulBigIntSer(&aMinus, &bPlus, &cMinusPlus, 0);
        if(aMinus.getCoefficientDigits() > 0 && bMinus.getCoefficientDigits() > 0)
                m.mulBigIntSer(&aMinus, &bMinus, &cMinusMinus, 0);
        cilk_sync;

        aPlus.freeHeap();       aMinus.freeHeap();
        bPlus.freeHeap();       bMinus.freeHeap();

        UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));

        for(i = 0; i < c.getSize(); ++i)
                c.setCoefficient(i, cPlusPlus.getCoefficient(i) - cPlusMinus.getCoefficient(i)
                                - cMinusPlus.getCoefficient(i) + cMinusMinus.getCoefficient(i));

        cPlusPlus.freeHeap();   cPlusMinus.freeHeap();
        cMinusPlus.freeHeap();  cMinusMinus.freeHeap();

        return c;

}

