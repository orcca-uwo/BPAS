/**
	Implementation of 'MulToom8.h'

	@author Farnam Mansouri
*/

#include "../../../include/IntegerPolynomial/Multiplication/MulToom8.h"


void MulToom8::growArrays(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
        int max = MAX(a->getSize(), b->getSize());
        int r = max % 8;
        if(r != 0) max += 8 - r;
        a->grow(max);
        b->grow(max);
}

/**
 * Auxilary function for intermediate computations.
 * Computing i0 * c0 + i1 * c1 + i2 * c2 + i3 * c3.
 */
mpz_class mult(mpz_class c0, mpz_class c1, mpz_class c2, mpz_class c3, mpz_class i0, mpz_class i1, mpz_class i2, mpz_class i3){
	return  i0 * c0 + i1 * c1 + i2 * c2 + i3 * c3;
}

/**
 * Auxilary function for intermediate computations.
 * Computing i0 * c0 + i1 * c1
 */
mpz_class mult(mpz_class c0, mpz_class c1, mpz_class i0, mpz_class i1){
	return i0 * c0 + i1 * c1;
}

/**
 * Auxilary function for intermediate computations.
 * Computing r1 = tmp1 + tmp2 and r2 = tmp1 - tmp2.
 */
void addSub(mpz_class *r1, mpz_class *r2, mpz_class tmp1, mpz_class tmp2){
        *r1 = tmp1 + tmp2;
        *r2 = tmp1 - tmp2;
}

/**
 * Auxilary function for intermediate computations.
 * Computing r = a + b.
 */
void add(mpz_class *r, mpz_class a, mpz_class b){
	*r = a + b;
}

/**
 * Auxilary function for intermediate computations.
 * Computing r = (a + b) / d
 */
void addDiv(mpz_class *r, mpz_class a, mpz_class b, mpz_class d){
	*r = (a + b) / d;
}

/**
 * Auxilary function for intermediate computations.
 * Computing r = (a + b + c + d) / div.
 */
void addDiv(mpz_class *r, mpz_class a, mpz_class b, mpz_class c, mpz_class d, mpz_class div){
	*r = (a + b + c + d) / div;
}

/**
 * Computing the first partition of evaluation.
 */
void intermEval1(mpz_class *tmp, mpz_class c1, mpz_class c2, mpz_class c3, mpz_class c4){
	tmp[0] = mult(c1, c2, c3, c4, 1, 1, 1, 1);
	tmp[1] = mult(c1, c2, c3, c4, 128, 32, 8, 2);
	tmp[2] = mult(c1, c2, c3, c4, 1, 4, 16, 64);
	tmp[3] = mult(c1, c2, c3, c4, 16384, 1024, 64, 4);
	tmp[4] = mult(c1, c2, c3, c4, 1, 16, 256, 4096);
	tmp[5] = mult(c1, c2, c3, c4, 2097152, 32768, 512, 8);
	tmp[6] = mult(c1, c2, c3, c4, 1, 64, 4096, 262144);
}

/**
 * Computing the second partition of evaluation.
 */
void intermEval2(mpz_class *tmp, mpz_class c1, mpz_class c2, mpz_class c3, mpz_class c4){
	tmp[0] = mult(c1, c2, c3, c4, 1, 1, 1, 1);
	tmp[1] = mult(c1, c2, c3, c4, 64, 16, 4, 1);
	tmp[2] = mult(c1, c2, c3, c4, 2, 8, 32, 128);
	tmp[3] = mult(c1, c2, c3, c4, 4096, 256, 16, 1);
	tmp[4] = mult(c1, c2, c3, c4, 4, 64, 1024, 16384);
	tmp[5] = mult(c1, c2, c3, c4, 262144, 4096, 64, 1);
	tmp[6] = mult(c1, c2, c3, c4, 8, 512, 32768, 2097152);
}


void MulToom8::evaluate(mpz_class * r, mpz_class *c){
	
	r[0 ] = c[0];

	mpz_class * tmp1 = new mpz_class[7];
	mpz_class * tmp2 = new mpz_class[7];

	// cilk_spawn 
		intermEval1(tmp1, c[0], c[2], c[4], c[6]);
		intermEval2(tmp2, c[1], c[3], c[5], c[7]);
	// cilk_sync;

	// cilk_spawn 
	addSub(&r[1], &r[2], tmp1[0], tmp2[0]);
	// cilk_spawn 
	addSub(&r[3], &r[4], tmp1[1], tmp2[1]);
	// cilk_spawn 
	addSub(&r[5], &r[6], tmp1[2], tmp2[2]);
	// cilk_spawn 
	addSub(&r[7], &r[8], tmp1[3], tmp2[3]);
	// cilk_spawn 
	addSub(&r[9], &r[10], tmp1[4], tmp2[4]);
	// cilk_spawn 
	addSub(&r[11], &r[12], tmp1[5], tmp2[5]);
	add(&r[13], tmp1[6], tmp2[6]);
	// cilk_sync;

	r[14] = c[7];

	delete[] tmp1; delete[] tmp2;
	
	/*r[0 ] = c0;
        r[1 ] = c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7;
        r[2 ] = c0 - c1 + c2 - c3 + c4 - c5 + c6 - c7;
        r[3 ] = 128 * c0 + 64 * c1 + 32 * c2 + 16 * c3 + 8 * c4 + 4 * c5 + 2 * c6 + c7;
        r[4 ] = 128 * c0 - 64 * c1 + 32 * c2 - 16 * c3 + 8 * c4 - 4 * c5 + 2 * c6 - c7;
        r[5 ] = c0 + 2 * c1 + 4 * c2 + 8 * c3 + 16 * c4 + 32 * c5 + 64 * c6 + 128 * c7;
        r[6 ] = c0 - 2 * c1 + 4 * c2 - 8 * c3 + 16 * c4 - 32 * c5 + 64 * c6 - 128 * c7;
        r[7 ] = 16384 * c0 + 4096 * c1 + 1024 * c2 + 256 * c3 + 64 * c4 + 16 * c5 + 4 * c6 + c7;
        r[8 ] = 16384 * c0 - 4096 * c1 + 1024 * c2 - 256 * c3 + 64 * c4 - 16 * c5 + 4 * c6 - c7;
        r[9 ] = c0 + 4 * c1 + 16 * c2 + 64 * c3 + 256 * c4 + 1024 * c5 + 4096 * c6 + 16384 * c7;
        r[10] = c0 - 4 * c1 + 16 * c2 - 64 * c3 + 256 * c4 - 1024 * c5 + 4096 * c6 - 16384 * c7;
        r[11] = 2097152 * c0 + 262144 * c1 + 32768 * c2 + 4096 * c3 + 512 * c4 + 64 * c5 + 8 * c6 + c7;
        r[12] = 2097152 * c0 - 262144 * c1 + 32768 * c2 - 4096 * c3 + 512 * c4 - 64 * c5 + 8 * c6 - c7;
        r[13] = c0 + 8 * c1 + 64 * c2 + 512 * c3 + 4096 * c4 + 32768 * c5 + 262144 * c6 + 2097152 * c7;
        r[14] = c7;*/
	
}

/**
 * Computing the intermediate results for interpolation.
 */
void intermFour1(mpz_class *tmp, mpz_class w014, mpz_class w1, mpz_class w2, mpz_class w3){
	tmp[0] = mult(w014, w1, w2, w3, -6056933132250, -13890539520, 10803752960, 101285184);
	tmp[1] = mult(w014, w1, w2, w3, 6056933132250, 13748072448, -10661285888, -96833088);
	tmp[2] = mult(w014, w1, w2, w3, -24227732529000, -17816217600, 5469071360, 47336976);
	tmp[3] = mult(w014, w1, w2, w3, -6056933132250, -13186037760, 10099251200, 80052864);
	tmp[4] = mult(w014, w1, w2, w3, -387643720464000, -111124316160, -86430023680, 405140736);
	tmp[5] = mult(w014, w1, w2, w3, 3028466566125, 5530892544, -3987499264, -14214984);
	tmp[6] = mult(w014, w1, w2, w3, 96910930116000, 36898971648, 12489613312, -119093568);
}

/**
 * Computing the intermediate results for interpolation.
 */
void intermFour2(mpz_class *tmp, mpz_class w0, mpz_class w122, mpz_class w34, mpz_class w56){
	tmp[0] = mult(w0, w122, w34, w56, -8075910843000, -23744512, 371008, 21824);
	tmp[1] = mult(w0, w122, w34, w56, -6153074928000, -359206912, 4625728, 345344);
	tmp[2] = mult(w0, w122, w34, w56, 6460728674400, 93933568, -1415488, -87104);
	tmp[3] = mult(w0, w122, w34, w56, 6080685811200, 1184186368, -5612608, -1332224);
	tmp[4] = mult(w0, w122, w34, w56, 6058412236800, 1519648768, -5936128, -5586944);
	tmp[5] = mult(w0, w122, w34, w56, -6062853888000, -1449459712, 5870848, 4542464);
}

/**
 * Computing the intermediate results for interpolation.
 */
void intermFour3(mpz_class *tmp, mpz_class w4, mpz_class w5, mpz_class w6, mpz_class w7){
	tmp[0] = mult(w4, w5, w6, w7, 84917184, -29724288, 17808384, 147444);
	tmp[1] = mult(w4, w5, w6, w7, -89369280, 29789760, -17873856, -180180);
	tmp[2] = mult(w4, w5, w6, w7, 357477120, 476636160, 285981696, -360360);
	tmp[3] = mult(w4, w5, w6, w7, -68136960, 29463120, -17547216, -47100);
	tmp[4] = mult(w4, w5, w6, w7, -71560896, -186202368, -4452096, 98274);
	tmp[5] = mult(w4, w5, w6, w7, 326640, 97926720, -50263104, -32730);
	tmp[6] = mult(w4, w5, w6, w7, 8257032, -14214984, 8257032, 8127);
}

/**
 * Computing the intermediate results for interpolation.
 */
void intermFour4(mpz_class *tmp, mpz_class w78, mpz_class w910, mpz_class w1112, mpz_class w14){
	tmp[0] = mult(w78, w910, w1112, w14, -5204, -64, 1, 23752678950);
	tmp[1] = mult(w78, w910, w1112, w14, 5396, 256, -1, -96141795750);
	tmp[2] = mult(w78, w910, w1112, w14, -1364, -4, 1, 1479104550);
	tmp[3] = mult(w78, w910, w1112, w14, -5444, -1024, 1, 403795542150);
	tmp[4] = mult(w78, w910, w1112, w14, 5456, 4096, -1, -2018977710750);
	tmp[5] = mult(w78, w910, w1112, w14, 4436, 16, -1, -5920755750);
}

/**
 * Computing the intermediate results for interpolation.
 */
void intermFour5(mpz_class *tmp, mpz_class w8, mpz_class w9, mpz_class w10, mpz_class w11){
	tmp[0] = mult(w8, w9, w10, w11, -10950, -63480, 19800, 5);
	tmp[1] = mult(w8, w9, w10, w11, 36180, -16350, 5430, 5);
	tmp[2] = mult(w8, w9, w10, w11, 169260, -16380, 5460, 65);
	tmp[3] = mult(w8, w9, w10, w11, -136524, 16374, -5454, -17);
	tmp[4] = mult(w8, w9, w10, w11, 76446, 229344, -54624, -17);
	tmp[5] = mult(w8, w9, w10, w11, -338520, -524160, -174720, 65);
	tmp[6] = mult(w8, w9, w10, w11, -2667, +8127, -2667, -1);	
}

/**
 * Computing the intermediate results for interpolation.
 */
void intermFour6(mpz_class *tmp, mpz_class w12, mpz_class w13){
	tmp[0] = mult(w12, w13, -15, -32);
	tmp[1] = mult(w12, w13, 63, 128);
	tmp[2] = mult(w12, w13, 3, 8);
	tmp[3] = mult(w12, w13, -3, 2);
	tmp[4] = mult(w12, w13, 15, -2);
	tmp[5] = mult(w12, w13, -63, 2);
}

void MulToom8::interpolate(mpz_class *w){

	std::cerr << "starting interpolate: " << std::endl;

	mpz_class w014 = w[0] + w[14];
        mpz_class w122 = w[1] + w[2];
        mpz_class w34 = w[3] + w[4];
        mpz_class w56 = w[5] + w[6];
        mpz_class w78 = w[7] + w[8];
        mpz_class w910 = w[9] + w[10];
        mpz_class w1112 = w[11] + w[12];

	mpz_class *tmp = new mpz_class[39];

	setNumberOfWorkers((char*) "6");

	// cilk_spawn 	
		intermFour1(tmp, w014, w[1], w[2], w[3]);
	// cilk_spawn
		intermFour2(&tmp[7], w[0], w122, w34, w56);
	// cilk_spawn
		intermFour3(&tmp[13], w[4], w[5], w[6], w[7]);
	// cilk_spawn
		intermFour4(&tmp[20], w78, w910, w1112, w[14]);
	// cilk_spawn
		intermFour5(&tmp[26], w[8], w[9], w[10], w[11]);
	intermFour6(&tmp[33], w[12], w[13]);
	// cilk_sync;
	
	setNumberOfWorkers((char*) "13");

	mpz_class tmp_w13 = -w[13];

        // cilk_spawn
                addDiv(&w[1 ], tmp[0], tmp[14], tmp[28], tmp[38], 48455465058000);
        // cilk_spawn
                addDiv(&w[2 ], tmp[7], tmp[22], 94662691200);
        // cilk_spawn
                addDiv(&w[3 ], tmp[1], tmp[13], tmp[29], tmp[37], 567976147200);
        // cilk_spawn
                addDiv(&w[4 ], tmp[9], tmp[25], 4441651200);
        // cilk_spawn
                addDiv(&w[5 ], tmp[3], tmp[16], tmp[27], tmp[36], 33312384000);
        // cilk_spawn
                addDiv(&w[6 ], tmp[8], tmp[20], 1045094400);
        // cilk_spawn
                addDiv(&w[7 ], tmp[5], tmp[19], tmp[32], tmp_w13, 4115059200);
        // cilk_spawn
                addDiv(&w[8 ], tmp[10], tmp[21], 1045094400);
        // cilk_spawn
                addDiv(&w[9 ], tmp[2], tmp[18], tmp[26], tmp[35], 33312384000);
        // cilk_spawn
                addDiv(&w[10], tmp[12], tmp[23], 4441651200);
        // cilk_spawn
                addDiv(&w[11], tmp[6], tmp[17], tmp[30], tmp[33], 567976147200);
        // cilk_spawn
                addDiv(&w[12], tmp[11], tmp[24], 94662691200);
        
                addDiv(&w[13], tmp[4], tmp[15], tmp[31], tmp[34], 48455465058000);
	
        // cilk_sync;

	// delete[] tmp;

	/*mpz_class c0  = w0;

        mpz_class c1  = (-6056933132250*w0 -13890539520*w1 +10803752960*w2 +101285184*w3 -89369280*w4 +29789760*w5 -17873856*w6 -180180*w7 +169260*w8 -16380*w9 +5460*w10 +65*w11 -63*w12 +2*w13 -6056933132250*w14)/48455465058000;
        mpz_class c2  = (-8075910843000*w0 -23744512*w1 -23744512*w2 +371008*w3 +371008*w4 +21824*w5 +21824*w6 -1364*w7 -1364*w8 -4*w9 -4*w10 +w11 +w12 +0*w13 +1479104550*w14)/94662691200;
        mpz_class c3  = (6056933132250*w0 +13748072448*w1 -10661285888*w2 -96833088*w3 +84917184*w4 -29724288*w5 +17808384*w6 +147444*w7 -136524*w8 +16374*w9 -5454*w10 -17*w11 +15*w12 -2*w13 +6056933132250*w14)/567976147200;
        mpz_class c4  = (6460728674400*w0 +93933568*w1 +93933568*w2 -1415488*w3 -1415488*w4 -87104*w5 -87104*w6 +4436*w7 +4436*w8 +16*w9 +16*w10 -w11 -w12 +0*w13 -5920755750*w14)/4441651200;
        mpz_class c5  = (-6056933132250*w0 -13186037760*w1 +10099251200*w2 +80052864*w3 -68136960*w4 +29463120*w5 -17547216*w6 -47100*w7 +36180*w8 -16350*w9 +5430*w10 +5*w11 -3*w12 +2*w13 -6056933132250*w14)/33312384000;
        mpz_class c6  = (-6153074928000*w0 -359206912*w1 -359206912*w2 +4625728*w3 +4625728*w4 +345344*w5 +345344*w6 -5204*w7 -5204*w8 -64*w9 -64*w10 +w11 +w12 +0*w13 +23752678950*w14)/1045094400;
        mpz_class c7  = (3028466566125*w0 +5530892544*w1 -3987499264*w2 -14214984*w3 +8257032*w4 -14214984*w5 +8257032*w6 +8127*w7 -2667*w8 +8127*w9 -2667*w10 -w11 +0*w12 -w13 +3028466566125*w14)/4115059200;
        mpz_class c8  = (6080685811200*w0 +1184186368*w1 +1184186368*w2 -5612608*w3 -5612608*w4 -1332224*w5 -1332224*w6 +5396*w7 +5396*w8 +256*w9 +256*w10 -w11 -w12 +0*w13 -96141795750*w14)/1045094400;
        mpz_class c9 = (-24227732529000*w0 -17816217600*w1 +5469071360*w2 +47336976*w3 +326640*w4 +97926720*w5 -50263104*w6 -32730*w7 -10950*w8 -63480*w9 +19800*w10 +5*w11 +3*w12 +8*w13 -24227732529000*w14)/33312384000;
        mpz_class c10 = (-6062853888000*w0 -1449459712*w1 -1449459712*w2 +5870848*w3 +5870848*w4 +4542464*w5 +4542464*w6 -5444*w7 -5444*w8 -1024*w9 -1024*w10 +w11 +w12 +0*w13 +403795542150*w14)/4441651200;
        mpz_class c11 = (96910930116000*w0 +36898971648*w1 +12489613312*w2 -119093568*w3 -71560896*w4 -186202368*w5 -4452096*w6 +98274*w7 +76446*w8 +229344*w9 -54624*w10 -17*w11 -15*w12 -32*w13 +96910930116000*w14)/567976147200;
        mpz_class c12 = (6058412236800*w0 +1519648768*w1 +1519648768*w2 -5936128*w3 -5936128*w4 -5586944*w5 -5586944*w6 +5456*w7 +5456*w8 +4096*w9 +4096*w10 -w11 -w12 +0*w13 -2018977710750*w14)/94662691200;
        mpz_class c13 = (-387643720464000*w0 -111124316160*w1 -86430023680*w2 +405140736*w3 +357477120*w4 +476636160*w5 +285981696*w6 -360360*w7 -338520*w8 -524160*w9 -174720*w10 +65*w11 +63*w12 +128*w13 -387643720464000*w14)/48455465058000;

        mpz_class c14 =  w14;*/
}


void MulToom8::mulToom8Par(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){

	setBitCounts(a, b);

        int digitCount = MAX_DIGITS_P(a, b) + 1;

	//DIVISION
	UnivariateIntegerPolynomial * aPoly = new UnivariateIntegerPolynomial[8];
	UnivariateIntegerPolynomial * bPoly = new UnivariateIntegerPolynomial[8];

	for(int i = 0; i < 8; ++i){
		aPoly[i].setSize(a->getSize() / 8);
		bPoly[i].setSize(b->getSize() / 8);
		aPoly[i].setCoefficientDigits(a->getCoefficientDigits());
		bPoly[i].setCoefficientDigits(b->getCoefficientDigits());
		aPoly[i].setCoefficients(a->getPointer(i * (a->getSize() / 8)));
		bPoly[i].setCoefficients(b->getPointer(i * (b->getSize() / 8)));
	}

	//CONVERSION
	mpz_class *aInt = new mpz_class[8];
	mpz_class *bInt = new mpz_class[8];

	for(int i = 0; i < 7; ++i){
		cilk_spawn 
			aPoly[i].getBigIntegerSigned(&aInt[i], digitCount);
		cilk_spawn
			bPoly[i].getBigIntegerSigned(&bInt[i], digitCount);
	}
	cilk_spawn 
		aPoly[7].getBigIntegerSigned(&aInt[7], digitCount);
		
		bPoly[7].getBigIntegerSigned(&bInt[7], digitCount);
	cilk_sync;

	for(int i = 0; i < 8; ++i){
		if(aPoly[i].isNegativeLC()) aInt[i] = -aInt[i];
		if(bPoly[i].isNegativeLC()) bInt[i] = -bInt[i];
	}

	for (int i = 0; i < 8; ++i) {
		std::cerr << "aInt[" << i << "]: " << mpz_sizeinbase(aInt[i].get_mpz_t(), 2) << std::endl;
	}
	for (int i = 0; i < 8; ++i) {
		std::cerr << "bInt[" << i << "]: " << mpz_sizeinbase(bInt[i].get_mpz_t(), 2) << std::endl;
	}
	//EVALUATION
	mpz_class * w_a = new mpz_class[15];
        mpz_class * w_b = new mpz_class[15];

	// cilk_spawn
                evaluate(w_a, aInt);
                evaluate(w_b, bInt);
        // cilk_sync;

	delete[] aInt; delete[] bInt;


	for (int i = 0; i < 15; ++i) {
		std::cerr << "w_a[" << i << "]: " << mpz_sizeinbase(w_a[i].get_mpz_t(), 2) << std::endl;
	}
for (int i = 0; i < 15; ++i) {
		std::cerr << "w_b[" << i << "]: " << mpz_sizeinbase(w_b[i].get_mpz_t(), 2) << std::endl;
	}

	//MULTIPLYING
	
	setNumberOfWorkers((char*) "15");

	mpz_class *w = new mpz_class[15];
	for(int i = 0; i < 14; i++)
		multiplyGMP(&w[i], w_a[i], w_b[i]);
	multiplyGMP(&w[14], w_a[14], w_b[14]);
	// cilk_sync;

	delete[] w_a; delete[] w_b;

	for (int i = 0; i < 15; ++i) {
		std::cerr << "w[" << i << "]: " << mpz_sizeinbase(w[i].get_mpz_t(), 2) << std::endl;
	}

	//INTERPOLATION
	
	interpolate(w);

	cilk_sync;

	//CONVERSION, MERGE
	int range = RESULT_DEGREE(aPoly[0], bPoly[0]);
        int rangeHalf = range / 2;

	delete[] aPoly; delete[] bPoly;

	UnivariateIntegerPolynomial cTmp(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
	
	int start = range - rangeHalf;

	setNumberOfWorkers((char*) "16");

	for(int i = 0; i < 15; i += 2){
		cilk_spawn 
			c->convertFromBigInteger(w[i], i * start, i * start + range, digitCount);
	}
	for(int i = 1; i < 13; i += 2){
		cilk_spawn
			cTmp.convertFromBigInteger(w[i], i * start, i * start + range, digitCount);
	}
	cTmp.convertFromBigInteger(w[13], 13 * start, 13 * start + range, digitCount);
        cilk_sync;

	delete[] w;

	int end = start * 13 + range;
        for(int i = start; i < end; ++i)
                c->setCoefficient(i, c->getCoefficient(i) + cTmp.getCoefficient(i));

        cTmp.freeHeap();

}
