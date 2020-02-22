/**
	Implementation of 'MulToom4.h'

	@author Farnam Mansouri
*/

#include "../../../include/IntegerPolynomial/Multiplication/MulToom4.h"


void MulToom4::growArrays(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
	int max = MAX(a->getSize(), b->getSize());
	int r = max % 4;
	if(r != 0) max += 4 - r;
	a->grow(max);
	b->grow(max);	
}


void MulToom4::evaluate(mpz_class * r, mpz_class * c){
        mpz_class tmp0 = c[0] + c[2];
        mpz_class tmp1 = c[1] + c[3];
        mpz_class tmp2 = 8 * c[0] + 2 * c[2];
        mpz_class tmp3 = 4 * c[1] + c[3];

        r[0] = c[0];
        r[1] = tmp0 + tmp1;
        r[2] = tmp0 - tmp1;
        r[3] = tmp2 + tmp3;
        r[4] = tmp2 - tmp3;
        r[5] = c[0] + 2 * c[1] + 4 * c[2] + 8 * c[3];
        r[6] = c[3];

	// SERIAL, NOT OPTIMIZED	
	/*mpz_class w_a0 = a0;
        mpz_class w_a1 = a0 + a1 + a2 + a3;
        mpz_class w_a2 = a0 - a1 + a2 - a3;
        mpz_class w_a3 = 8 * a0 + 4 * a1 + 2 * a2 + a3;
        mpz_class w_a4 = 8 * a0 - 4 * a1 + 2 * a2 - a3;
        mpz_class w_a5 = a0 + 2 * a1 + 4 * a2 + 8 * a3;
        mpz_class w_a6 = a3;

        mpz_class w_b0 = b0;
        mpz_class w_b1 = b0 + b1 + b2 + b3;
        mpz_class w_b2 = b0 - b1 + b2 - b3;
        mpz_class w_b3 = 8 * b0 + 4 * b1 + 2 * b2 + b3;
        mpz_class w_b4 = 8 * b0 - 4 * b1 + 2 * b2 - b3;
        mpz_class w_b5 = b0 + 2 * b1 + 4 * b2 + 8 * b3;
        mpz_class w_b6 = b3;*/
}


void MulToom4::interpolate(mpz_class * w){

	// SERIAL, CONVERTING TO THE IDENTITY MATRIX
	w[5] += w[3]; w[4] += w[3]; w[4] /= 2; w[2] += w[1]; w[2] /= 2; w[3] -= w[4];
        w[1] -= w[2]; w[5] -= 65 * w[2]; w[2] -= w[6]; w[2] -= w[0];
        w[5] += 45 * w[2]; w[4] -= w[6]; w[4] /= 4; w[3] /= 2; w[5] -= 4 * w[3];
        w[3] -= w[1]; w[3] /= 3; w[4] -= 16 * w[0];
        w[4] -= 4 * w[2]; w[4] /= -3; w[2] -= w[4];
        w[5] /= 30; w[1] -= w[5]; w[1] -= w[3]; w[1] /= -3;
        w[3] -= 5 * w[1]; w[5] += w[1];

	// SERIAL, NOT OPTIMIZED, USING INVERSE OF MATRIX
	/*mpz_class c0 = w0;
        mpz_class c1 = (-180 * w0 - 120 * w1 + 40 * w2 + 10 * w3 - 6 * w4 + 4 * w5 - 180 * w6) / 360;
        mpz_class c2 = (- 1800 * w0 - 60 * w1 - 60 * w2 + 15 * w3 + 15 * w4 + 0 * w5 + 90 * w6) / 360;
        mpz_class c3 = (900 * w0 + 540 * w1 - 140 * w2 - 20 * w3 + 0 * w4 - 20 * w5 + 900 * w6) / 360;
        mpz_class c4 = (1440 * w0 + 240 * w1 + 240 * w2 - 15 * w3 - 15 * w4 + 0 * w5 - 450 * w6) / 360;
        mpz_class c5 = (- 720 * w0 - 240 * w1 - 80 * w2 + 10 * w3 + 6 * w4 + 16 * w5 - 720 * w6) / 360;
        mpz_class c6 = w6;*/
}


void MulToom4::mulToom4Par(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){

	setBitCounts(a, b);
	
        int digitCount = MAX_DIGITS_P(a, b) + 1;

	//DIVISION
	UnivariateIntegerPolynomial * aPoly = new UnivariateIntegerPolynomial[4];
        UnivariateIntegerPolynomial * bPoly = new UnivariateIntegerPolynomial[4];

        for(int i = 0; i < 4; ++i){
                aPoly[i].setSize(a->getSize() / 4);
                bPoly[i].setSize(b->getSize() / 4);
                aPoly[i].setCoefficientDigits(a->getCoefficientDigits());
                bPoly[i].setCoefficientDigits(b->getCoefficientDigits());
		aPoly[i].setCoefficients(a->getPointer(i * (a->getSize() / 4)));
                bPoly[i].setCoefficients(b->getPointer(i * (b->getSize() / 4)));
        }
	//CONVERSION
	
	mpz_class *aInt = new mpz_class[4];
        mpz_class *bInt = new mpz_class[4];

        cilk_for(int i = 0; i < 4; ++i) {
                aPoly[i].getBigIntegerSigned(&aInt[i], digitCount);
                bPoly[i].getBigIntegerSigned(&bInt[i], digitCount);                
        }
        // for(int i = 0; i < 3; ++i){
        //         cilk_spawn
        //                 aPoly[i].getBigIntegerSigned(&aInt[i], digitCount);
        //         cilk_spawn
        //                 bPoly[i].getBigIntegerSigned(&bInt[i], digitCount);
        // }
        // cilk_spawn
        //         aPoly[3].getBigIntegerSigned(&aInt[3], digitCount);

        //         bPoly[3].getBigIntegerSigned(&bInt[3], digitCount);
        // cilk_sync;

        for(int i = 0; i < 4; ++i){
                if(aPoly[i].isNegativeLC()) aInt[i] = -aInt[i];
                if(bPoly[i].isNegativeLC()) bInt[i] = -bInt[i];
        }

	//EVALUATION
	mpz_class * w_a = new mpz_class[7];
        mpz_class * w_b = new mpz_class[7];

	setNumberOfWorkers((char*) "2");

        cilk_spawn
                evaluate(w_a, aInt);
                evaluate(w_b, bInt);
        cilk_sync;

	delete[] aInt; delete[] bInt;

	//MULTIPLYING

	setNumberOfWorkers((char*) "7");

        // for (int i = 0; i < 7; ++i) {
        //         std::cerr << "w_a[" << i << "]: " << mpz_sizeinbase(w_a[i].get_mpz_t(), 2) << std::endl;
        // }
        // for (int i = 0; i < 7; ++i) {
        //         std::cerr << "w_b[" << i << "]: " << mpz_sizeinbase(w_b[i].get_mpz_t(), 2) << std::endl;
        // }


	mpz_class *w = new mpz_class[7];


        // mpz_class w0,w1,w2,w3,w4,w5,w6;

        // cilk_spawn 
        //         multiplyGMP(&w0, w_a[0], w_b[0]);
        // cilk_spawn
        //         multiplyGMP(&w1, w_a[1], w_b[1]);
        // cilk_spawn
        //         multiplyGMP(&w2, w_a[2], w_b[2]);
        // cilk_spawn
        //         multiplyGMP(&w3, w_a[3], w_b[3]);
        // cilk_spawn
        //         multiplyGMP(&w4, w_a[4], w_b[4]);
        // cilk_spawn
        //         multiplyGMP(&w5, w_a[5], w_b[5]);
        // multiplyGMP(&w6, w_a[6], w_b[6]);
        // cilk_sync;

        // w[0] = w0;
        // w[1] = w1;
        // w[2] = w2;
        // w[3] = w3;
        // w[4] = w4;
        // w[5] = w5;
        // w[6] = w6;

	// for(int i = 0; i < 6; i++)
 //                multiplyGMPI(w, i, w_a[i], w_b[i]);
 //                // cilk_spawn 
	// multiplyGMPI(w, 6, w_a[6], w_b[6]);
	// cilk_sync;

        cilk_for(int i = 0; i < 7; ++i) {
                multiplyGMP(&w[i], w_a[i], w_b[i]);
        }

        for (int i = 0; i < 7; ++i) {
                std::cerr << "w[" << i << "]: " << mpz_sizeinbase(w[i].get_mpz_t(), 2) << std::endl;
        }

	delete[] w_a; delete[] w_b;

	//INTERPOLATION
	interpolate(w);

	//CONVERSION, MERGE
	int range = RESULT_DEGREE(aPoly[0], bPoly[0]);
        int rangeHalf = range / 2;

	delete[] aPoly; delete[] bPoly;

	UnivariateIntegerPolynomial cTmp(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));

        int start = range - rangeHalf;

	setNumberOfWorkers((char*) "8");

	cilk_for(int i = 0; i < 7; i += 2) {
        // cilk_spawn
                   c->convertFromBigInteger(w[i], i * start, i * start + range, digitCount);
        }
        cilk_for(int i = 1; i <= 5; i += 2) {
        // for(int i = 1; i < 5; i += 2) {
                // cilk_spawn
                        cTmp.convertFromBigInteger(w[i], i * start, i * start + range, digitCount);
        }
        // cTmp.convertFromBigInteger(w[5], 5 * start, 5 * start + range, digitCount);
        // cilk_sync;

	delete[] w;


        int end = start * 5 + range;
        for(int i = start; i < end; ++i)
                c->setCoefficient(i, c->getCoefficient(i) + cTmp.getCoefficient(i));

        cTmp.freeHeap();

}
