/**
	Implementation of 'MulNaive.cpp'

	@author Farnam Mansouri
*/

#include "../../../include/IntegerPolynomial/Multiplication/MulNaive.h"

MulNaive::MulNaive(){
        // Default Threshold
	BASE_MUL = 512;
}


void MulNaive::mulNaiveSer(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){
        for (int x = 0; x < a->getSize(); ++x)
                for (int y = 0; y < b->getSize(); ++y)
                        c->setCoefficient(startIndex + x + y, c->getCoefficient(startIndex + x + y) + a->getCoefficient(x) * b->getCoefficient(y));
}


void MulNaive::mulDnRNaivePar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){

        if(a->getSize() <= BASE_MUL && b->getSize() <= BASE_MUL){
                mulNaiveSer(a, b, c, startIndex);
        } else {

                UnivariateIntegerPolynomial a0(a->getSize() / 2, a->getCoefficientDigits(), a->getCoefficients());
		mpz_class a1Start = a->getCoefficient(a0.getSize());
                UnivariateIntegerPolynomial a1(a->getSize() - a0.getSize(), a->getCoefficientDigits(), &a1Start);
                UnivariateIntegerPolynomial b0(b->getSize() / 2, b->getCoefficientDigits(), b->getCoefficients());
		mpz_class b1Start = b->getCoefficient(b0.getSize());
                UnivariateIntegerPolynomial b1(b->getSize() - b0.getSize(), b->getCoefficientDigits(), &b1Start);

                cilk_spawn mulDnRNaivePar(&a0, &b0, c, startIndex);
                cilk_spawn mulDnRNaivePar(&a1, &b1, c, startIndex + a0.getSize() + b0.getSize());

                UnivariateIntegerPolynomial tmp0(RESULT_DEGREE(a0, b1), MAX_DIGITS(a0, b1));
                UnivariateIntegerPolynomial tmp1(RESULT_DEGREE(a1, b0), MAX_DIGITS(a1, b0));

                cilk_spawn mulDnRNaivePar(&a0, &b1, &tmp0, 0);
                mulDnRNaivePar(&a1, &b0, &tmp1, 0);
                cilk_sync;

                int stId = startIndex + b0.getSize();
                cilk_for(int i = 0; i < a0.getSize() + b1.getSize() - 1; ++i)
                        c->setCoefficient(stId + i , c->getCoefficient(stId + i) + tmp0.getCoefficient(i));

                stId = startIndex + a0.getSize();
                cilk_for(int i = 0; i < a1.getSize() + b0.getSize() - 1; ++i)
                        c->setCoefficient(stId + i , c->getCoefficient(stId + i) + tmp1.getCoefficient(i));


                tmp0.freeHeap(); tmp1.freeHeap();
        }
}


void MulNaive::mulIterNaivePar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex){

        if(a->getSize() <= BASE_MUL && b->getSize() <= BASE_MUL){
                mulNaiveSer(a, b, c, startIndex);
        } else {

		int aPolyCount = a->getSize() / BASE_MUL;
                int aLastPolySize = a->getSize() % BASE_MUL;
                int bPolyCount = b->getSize() / BASE_MUL;
                int bLastPolySize = b->getSize() % BASE_MUL;

                int aStart, bStart;

		mpz_class aS, bS;

                for(int i = 0; i < aPolyCount; ++i){
                        for(int j = 0; j < bPolyCount; ++j){

                                aStart = i * BASE_MUL;
                                bStart = j * BASE_MUL;

				aS = a->getCoefficient(aStart);
				bS = b->getCoefficient(bStart);

                                UnivariateIntegerPolynomial a0(BASE_MUL, a->getCoefficientDigits(), &aS);
                                UnivariateIntegerPolynomial b0(BASE_MUL, b->getCoefficientDigits(), &bS);

                                int destStart = aStart + bStart;

                                mulNaiveSer(&a0, &b0, c, startIndex + destStart);

                        }
                }

                if(aLastPolySize != 0){
                        aStart = aPolyCount * BASE_MUL;
                        for(int j = 0; j < bPolyCount; ++j){
                                bStart = j * BASE_MUL;

				aS = a->getCoefficient(aStart);
				bS = b->getCoefficient(bStart);

                                UnivariateIntegerPolynomial a0(aLastPolySize, a->getCoefficientDigits(), &aS);
                                UnivariateIntegerPolynomial b0(BASE_MUL, b->getCoefficientDigits(), &bS);

                                int destStart = aStart + bStart;

                                mulNaiveSer(&a0, &b0, c, startIndex + destStart);
                        }
                }

                if(bLastPolySize != 0){
                        bStart = bPolyCount * BASE_MUL;
                        for(int i = 0; i < aPolyCount; ++i){
                                aStart = i * BASE_MUL;

				aS = a->getCoefficient(aStart);
				bS = b->getCoefficient(bStart);

                                UnivariateIntegerPolynomial a0(BASE_MUL, a->getCoefficientDigits(), &aS);
                                UnivariateIntegerPolynomial b0(bLastPolySize, b->getCoefficientDigits(), &bS);

                                int destStart = aStart + bStart;

                                mulNaiveSer(&a0, &b0, c, startIndex + destStart);
                        }
                        if(aLastPolySize != 0){
                                aStart = aPolyCount * BASE_MUL;

				aS = a->getCoefficient(aStart);
				bS = b->getCoefficient(bStart);

                                UnivariateIntegerPolynomial a0(aLastPolySize, a->getCoefficientDigits(), &aS);
                                UnivariateIntegerPolynomial b0(bLastPolySize, b->getCoefficientDigits(), &bS);

                                int destStart = aStart + bStart;

                                mulNaiveSer(&a0, &b0, c, startIndex + destStart);
                        }
                }
        }
}
