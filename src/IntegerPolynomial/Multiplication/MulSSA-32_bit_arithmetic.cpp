/**
 * Implementation of 'MulSSA.h'
 *
 * @author Ning Xie and Yuzhen Xie
 */

#include "../../../include/IntegerPolynomial/Multiplication/MulSSA.h"

/*----------------------------------------
Our new 2-convolution multiplication
-----------------------------------------*/
void MulSSA::mul2C(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) {
	cilk_spawn
		a->setCoefficientDigits(a->maxCoefficientSize());
		b->setCoefficientDigits(b->maxCoefficientSize());
	cilk_sync;

	int digitCount = MAX_COEFF_DIGITS(a, b);
	sfixn d = MAX(a->getSize(), b->getSize());
	determineGoodN(digitCount, d);

//std::cout << "N: " << N << ", K: " << K << ", M: " << M << std::endl;

	if (N < 4096)
		mul2C2(a, b, c);
	else
		mul2C3(a, b, c);
}

void MulSSA::mul2C2(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c)
{//---------------------------------------------------------------

//unsigned long long start, end;

  sfixn d1 = a->getSize();
  sfixn d2 = b->getSize();

//start = __cilkview_getticks();
  BivariatePolynomial *biA = cilk_spawn 
    ToBivarTwoMod(a->getCoefficients(), d1, prime1, prime2);
  BivariatePolynomial *biB = ToBivarTwoMod(b->getCoefficients(), d2, prime1, prime2);
  cilk_sync;

//end = __cilkview_getticks();
//std::cout << "ToBivarTwoMod (sec): \t" << (end - start) /1000.f << std::endl;

//start = __cilkview_getticks();
  // cyclic and negacyclic convolution
  sfixn *ncc1 = PBPAS::TwoConvolutionMod(biA->getCoefficients(), biB->getCoefficients(), biA->getCoefficients()+biA->getSize(), biB->getCoefficients()+biB->getSize(), d1, d2, prime1, prime2, biA->getK());
//end = __cilkview_getticks();
//cout << "TwoConvolutionMod (sec): \t" << (end - start) /1000.f << endl;

  biA->freeHeap();
  biB->freeHeap();

//start = __cilkview_getticks();
  CRT_ToGMP_Recovering2(ncc1, d1, d2, c->getCoefficients());
//end = __cilkview_getticks();
//cout << "CRT_ToGMP_Recovering (sec): \t" << (end - start) /1000.f << endl;

  my_free(ncc1);
}

void MulSSA::mul2C3(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) {
//unsigned long long start, end;
	sfixn d1 = a->getSize();
	sfixn d2 = b->getSize();

//start = __cilkview_getticks();
	BivariatePolynomial *biA = cilk_spawn toBivariateMod3(a->getCoefficients(), d1, prime1, prime2, prime3);
	BivariatePolynomial *biB = toBivariateMod3(b->getCoefficients(), d2, prime1, prime2, prime3);
	cilk_sync;
//end = __cilkview_getticks();
//cout << "ConvertIn (sec): \t" << (end - start) /1000.f << endl;

//start = __cilkview_getticks();
	// cyclic and negacyclic convolution
	sfixn *ncc1 = PBPAS::TwoConvolutionMod3(biA->getCoefficients(), biB->getCoefficients(), d1, d2, prime1, prime2, prime3, biA->getK());
//end = __cilkview_getticks();
//cout << "TwoConvolution (sec): \t" << (end - start) /1000.f << endl;

	biA->freeHeap();
	biB->freeHeap();

//start = __cilkview_getticks();
	CRT_ToGMP_Recovering(ncc1, d1, d2, c->getCoefficients());
//end = __cilkview_getticks();
//cout << "ConvertOut (sec): \t" << (end - start) /1000.f << endl;

	delete [] ncc1;
}

/**
 * Convert from CRT representation to a mpz_t object
 *
 * Output:
 * @zp: Big integer
 *
 * Input:
 * @a: An array mod prime1: a[0] + a[1] * 2^M + ...
 * @b: An array mod prime2: b[0] + b[1] * 2^M + ...
 * @c: An array mod prime3: c[0] + c[1] * 2^M + ...
 * @size: Size of each array
 **/
void MulSSA::CRTtoMPZ(mpz_t zp, sfixn* a, sfixn* b, sfixn* c, int size) {
	mpz_t zn;
	mpz_init2(zp, LIMB_BITS);
	mpz_init2(zn, LIMB_BITS);

	int idx = 0;
	unsigned __int128 postainer[2] = {0}, negtainer[2] = {0};
	sfixn shifter = 0;
	for (int i = 0; i < size; ++i) {
		__int128 elem = (__int128) (((long int) a[i] * U23) % prime1) * P2_P3
				+ (__int128) (((long int) b[i] * U31) % prime2) * P3_P1
				+ (__int128) (((long int) c[i] * U12) % prime3) * P1_P2;
		if (elem > HALF_P1_P2_P3) { elem -= P1_P2_P3; }
		else if (elem < N_HALF_P1_P2_P3) { elem += P1_P2_P3; }

		if (elem < 0) {
			elem = -elem;
			unsigned __int128 tmp = elem << shifter;
			negtainer[0] += tmp;
			bool carry = negtainer[0] < tmp;
			if (shifter)
				tmp = elem >> (128 - shifter);
			else { tmp = 0; }
			negtainer[1] += tmp + carry;
		}
		else if (elem > 0) {
			unsigned __int128 tmp = elem << shifter;
			postainer[0] += tmp;
			bool carry = postainer[0] < tmp;
			if (shifter)
				tmp = elem >> (128 - shifter);
			else { tmp = 0; }
			postainer[1] += tmp + carry;
		}
		shifter += M;

		if (shifter >= 128) {
			if (postainer[0] > 0) {
				zp->_mp_d[idx] = (mp_limb_t) postainer[0];
				zp->_mp_d[idx+1] = (mp_limb_t) (postainer[0] >> GMP_LIMB_BITS);
				zp->_mp_size = idx + 2;
			}
			else { zp->_mp_d[idx] = 0; zp->_mp_d[idx+1] = 0; }
			postainer[0] = postainer[1];
			postainer[1] = 0;

			if (negtainer[0] > 0) {
				zn->_mp_d[idx] = (mp_limb_t) negtainer[0];
				zn->_mp_d[idx+1] = (mp_limb_t) (negtainer[0] >> GMP_LIMB_BITS);
				zn->_mp_size = idx + 2;
			}
			else { zn->_mp_d[idx] = 0; zn->_mp_d[idx+1] = 0; }
			negtainer[0] = negtainer[1];
			negtainer[1] = 0;

			shifter -= 128;
			idx += 2;
		}
	}

	if (postainer[0] > 0) {
		zp->_mp_d[idx] = (mp_limb_t) postainer[0];
		zp->_mp_d[idx+1] = (mp_limb_t) (postainer[0] >> GMP_LIMB_BITS);
		zp->_mp_size = idx + 2;
	}
	if (negtainer[0] > 0) {
		zn->_mp_d[idx] = (mp_limb_t) negtainer[0];
		zn->_mp_d[idx+1] = (mp_limb_t) (negtainer[0] >> GMP_LIMB_BITS);
		zn->_mp_size = idx + 2;
	}
	idx += 2;
	if (postainer[1] > 0) {
		zp->_mp_d[idx] = (mp_limb_t) postainer[1];
		zp->_mp_d[idx+1] = (mp_limb_t) (postainer[1] >> GMP_LIMB_BITS);
		zp->_mp_size = idx + 2;
	}
	if (negtainer[1] > 0) {
		zn->_mp_d[idx] = (mp_limb_t) negtainer[1];
		zn->_mp_d[idx+1] = (mp_limb_t) (negtainer[1] >> GMP_LIMB_BITS);
		zn->_mp_size = idx + 2;
	}

	if (zp->_mp_size && !zp->_mp_d[zp->_mp_size-1])
		zp->_mp_size--;
	if (zn->_mp_size && !zn->_mp_d[zn->_mp_size-1])
		zn->_mp_size--;

	if (zn > 0) { mpz_sub (zp, zp, zn); }
	mpz_clear(zn);
}

/**
 * Convert from (a mod p1, b mod p2) to GMP class
 *
 * @size: the size of array a & b
 **/
void MulSSA::CRTtoMPZ(mpz_t zp, sfixn* a, sfixn* b, int size) {
	mpz_t zn;
	mpz_init2(zp, LIMB_BITS);
	mpz_init2(zn, LIMB_BITS);
  
	int idx = 0;
	__int128 Tp = 0, Tn = 0;
	sfixn shifter = 0;
	for (int i = 0; i < size; ++i) {
		long int elem = (long int) ((b[i] - a[i]) % prime2) * U1 % prime2;
		elem = a[i] + elem * prime1;
		if (elem > HALF_P1_P2) { elem -= P1_P2; }
		else if (elem < N_HALF_P1_P2) { elem += P1_P2; }
		if (elem < 0)
			Tn += (__int128) (-elem) << shifter;
		else if (elem)
			Tp += (__int128) elem << shifter;
		shifter += M;

		if (shifter >= GMP_LIMB_BITS) {
			if (Tp > 0) {
				zp->_mp_d[idx] = (mp_limb_t) Tp;
				zp->_mp_size = idx + 1;
			}
			else { zp->_mp_d[idx] = 0; }
			Tp >>= GMP_LIMB_BITS;
			if (Tn > 0) {
				zn->_mp_d[idx] = (mp_limb_t) Tn;
				zn->_mp_size = idx + 1;
			}
			else { zn->_mp_d[idx] = 0; }
			Tn >>= GMP_LIMB_BITS;
			shifter -= GMP_LIMB_BITS;
			idx++;
		}
	}
	if (Tp > 0) {
		zp->_mp_d[idx] = (mp_limb_t) Tp;
		zp->_mp_size = idx + 1;
		Tp >>= GMP_LIMB_BITS;
		if (Tp > 0) {
			zp->_mp_d[idx+1] = (mp_limb_t) Tp;
			zp->_mp_size = idx + 2;
		}
	}
	if (Tn > 0) {
		zn->_mp_d[idx] = (mp_limb_t) Tn;
		zn->_mp_size = idx + 1;
		Tn >>= GMP_LIMB_BITS;
		if (Tn > 0) {
			zn->_mp_d[idx+1] = (mp_limb_t) Tn;
			zn->_mp_size = idx + 2;
		}
	}
	if (zn > 0) { mpz_sub (zp, zp, zn); }
	mpz_clear(zn);
}

/**
 * for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
 * which is a large integer encoded as a polynomial in x, i.e.
 * ncc1_i(x), ncc2_i(x), ncc3_i(x), cc1_i(x), cc2_i(x), cc3(x)
 * (1) for each coefficient of x, apply CRT for the three prime numbers,
 *     get ncc_i and cc_i
 * (2) convert ncc_i and cc_i to GMP, get u_i and v_i
 * (3) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N
 *
 * Output:
 * @c: An array storing big integer coefficients, that is c = a * b
 *
 * Input:
 * @ncc1: An array storing ncc1_i(x), ncc2_i(x), ncc3_i(x), cc1_i(x), cc2_i(x), cc3(x)
 * @d1: Degree of polynomial a
 * @d2: Degree of polynomial b
 **/
void MulSSA::CRT_ToGMP_Recovering(sfixn *ncc1, int d1, int d2, mpz_class *c) {
	int csize = d1 + d2 - 1;
	int fullsize = K * csize;
	sfixn *ncc2 = ncc1 + fullsize;
	sfixn *ncc3 = ncc2 + fullsize;
	sfixn *cc1 = ncc3 + fullsize;
	sfixn *cc2 = cc1 + fullsize;
	sfixn *cc3 = cc2 + fullsize;
  
	#pragma cilk_grainsize = 1024;
	cilk_for (int i = 0; i < csize; ++i) {
		int Ks = i * K;
		mpz_t u, v;
		CRTtoMPZ(u, ncc1+Ks, ncc2+Ks, ncc3+Ks, K);
		CRTtoMPZ(v, cc1+Ks, cc2+Ks, cc3+Ks, K);

		mpz_add (c[i].get_mpz_t(), u, v);
		c[i] >>= 1;
		mpz_sub (v, v, u);
		c[i] += mpz_class (v) << (N - 1);

		mpz_clear(u);
		mpz_clear(v);
	}
}

/*---------------------------------------------------------------------
  for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
  which is a large integer encoded as a polynomial in x, i.e.
  ncc1_i(x), ncc2_i(x), cc1_i(x), cc2_i(x)
  (1) for each coefficient of x, apply CRT for the two prime numbers, 
      get ncc_i and cc_i
  (2) convert ncc_i and cc_i to GMP, get u_i and v_i
  (3) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N
------------------------------------------------------------------------*/
void MulSSA::CRT_ToGMP_Recovering2(sfixn *ncc1, int d1, int d2, mpz_class *c) {
	int csize = d1 + d2 - 1;
	int fullsize = K * csize;
	sfixn *ncc2 = ncc1 + fullsize;
	sfixn *cc1 = ncc2 + fullsize;
	sfixn *cc2 = cc1 + fullsize;
  
	#pragma cilk_grainsize = 1024;
	cilk_for (int i = 0; i < csize; ++i) {
		int Ks = i * K;
		mpz_t u, v;
		CRTtoMPZ(u, ncc1+Ks, ncc2+Ks, K);
		CRTtoMPZ(v, cc1+Ks, cc2+Ks, K);

		mpz_add (c[i].get_mpz_t(), u, v);
		c[i] >>= 1;
		mpz_sub (v, v, u);
		c[i] += mpz_class (v) << (N - 1);

		mpz_clear(u);
		mpz_clear(v);
	}
}

/**
 * Convert a univariate polynomial with large coefficients to 
 * a bivariate polynomial representation with relatively 
 * small coefficients.
 *
 * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
 * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
 *
 * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
 *
 * @param coeff the coefficients of large integers
 * @param d partial degree of the univariate polynomial plus one
 */
BivariatePolynomial* MulSSA::toBivariateMod3(mpz_class *coeff, int d, sfixn p1, sfixn p2, sfixn p3) {
	// 3 times bigger size
	BivariatePolynomial *biA = new BivariatePolynomial(d, K, M, 3);

	#pragma cilk_grainsize = 1024;
	cilk_for (int i = 0; i < d; ++i) {
		if (coeff[i] != 0) {
		        int base = i * K, base1 = biA->getSize() + base, base2 = biA->getSize() + base1;
		        sfixn* a1 = biA->getCoefficients() + base;
		        sfixn* a2 = biA->getCoefficients() + base1;
		        sfixn* a3 = biA->getCoefficients() + base2;
			mpzToPolyMod3(a1, a2, a3, coeff[i].get_mpz_t(), M, p1, p2, p3);
		}
	}
  
	return biA;
}

void MulSSA::mpzToPolyMod3(sfixn* a1, sfixn* a2, sfixn* a3, mpz_t coef, int M, sfixn p1, sfixn p2, sfixn p3) {
	bool isNeg = 0;
	int s = coef->_mp_size; //abs(s)>=1
	if (s < 0) {
		isNeg = 1;
		s = -s;
	}

	mp_ptr limbs = coef->_mp_d;
	mp_limb_t carry = 0;
	mp_size_t carry_size = 0;

	// all the limbs except for the last one are full-size (64 bits)
	int idx = 0;
	for (int i = 0; i < s-1; ++i) {
		mp_limb_t w = limbs[i], q;
		size_t k = 0; // current position in w
		if (!carry_size) {
			while (GMP_LIMB_BITS - k >= M) {
				q = w >> M;
				a1[idx] = (sfixn) ((w) ^ (q << M));
				if (isNeg) {
					a1[idx] = -a1[idx];
					a2[idx] = a1[idx] + p2;
					a3[idx] = a1[idx] + p3;
					a1[idx] = a1[idx] + p1;
				}
				else {
					a2[idx] = a1[idx];
					a3[idx] = a1[idx];
				}
				idx++;
				w = q;
				k += M;
			}
		}
		else {
			size_t l = M - carry_size;
			while (GMP_LIMB_BITS - k >= l) {
				q = w >> l;
				a1[idx] = (sfixn) ((((w) ^ (q << l)) << carry_size) + carry); //w % 2^a + carry
				if (isNeg) {
					a1[idx] = -a1[idx];
					a2[idx] = a1[idx] + p2;
					a3[idx] = a1[idx] + p3;
					a1[idx] = a1[idx] + p1;
				}
				else {
					a2[idx] = a1[idx];
					a3[idx] = a1[idx];
				}
				w = q;
				idx++;
				k += l;

				l = M;
				carry_size = 0;
				carry = 0;
			}
		}
		carry_size = GMP_LIMB_BITS - k;
		carry = w;
	}

	// last limb, its bit size may be smaller than long
	mp_limb_t w = limbs[s-1];
	size_t b = log2(w)+1; // num of bits
	size_t k = 0; // current position in w

	if (!carry_size) {
		while (b - k >= M) {
			mp_limb_t q = w >> M;
			a1[idx] = (sfixn) ((w) ^ (q << M));
			if (isNeg) {
				a1[idx] = -a1[idx];
				a2[idx] = a1[idx] + p2;
				a3[idx] = a1[idx] + p3;
				a1[idx] = a1[idx] + p1;
			}
			else {
				a2[idx] = a1[idx];
				a3[idx] = a1[idx];
			}
			idx++;
			w = q;
			k += M;
		}
		if (w) {
			a1[idx] = (sfixn) w;
			if (isNeg) {
				a1[idx] = -a1[idx];
				a2[idx] = a1[idx] + p2;
				a3[idx] = a1[idx] + p3;
				a1[idx] = a1[idx] + p1;
			}
			else {
				a2[idx] = a1[idx];
				a3[idx] = a1[idx];
			}
		}
	}
	else {
		size_t l = M - carry_size;
		while (b - k >= l) {
			mp_limb_t q = w >> l; //quotient = w / 2^a
			a1[idx] = (sfixn) ((((w) ^ (q << l)) << carry_size) + carry); //w % 2^a + carry
			if (isNeg) {
				a1[idx] = -a1[idx];
				a2[idx] = a1[idx] + p2;
				a3[idx] = a1[idx] + p3;
				a1[idx] = a1[idx] + p1;
			}
			else {
				a2[idx] = a1[idx];
				a3[idx] = a1[idx];
			}
			w = q;
			k += l;
			idx++;

			l = M;
			carry_size = 0;
			carry = 0;
		}
		if (w) {
			a1[idx] = (sfixn) ((w << carry_size) + carry);
			if (isNeg) {
				a1[idx] = -a1[idx];
				a2[idx] = a1[idx] + p2;
				a3[idx] = a1[idx] + p3;
				a1[idx] = a1[idx] + p1;
			}
			else {
				a2[idx] = a1[idx];
				a3[idx] = a1[idx];
			}
		}
	}
}

/**
 * Convert a univariate polynomial with large coefficients to 
 * a bivariate polynomial representation with relatively 
 * small coefficients.
 *
 * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
 * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
 *
 * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
 *
 * @param coeff the coefficients of large integers
 * @param d partial degree of the univariate polynomial plus one
 */
BivariatePolynomial * MulSSA::ToBivarTwoMod(mpz_class *coeff, int d, 
					    int p1, int p2)
{//-------------------------------------------------------------------
  BivariatePolynomial *biA = new BivariatePolynomial(d, K, M, 2);//double size
  sfixn *A2 = biA->getCoefficients() + biA->getSize(); 
  // biA->coefficients has length of K*d*2

  //convert a large int coeff[i] to a univar poly of degree K-1
  //by a base of 2^M
  //store the coefficients of the univar poly to the coeff of y^i
  #pragma cilk_grainsize = 1024;
  cilk_for(int i=0; i<d; ++i){
      sfixn ci = i*K;
      sfixn *Ai = biA->getCoefficients() + ci;
      sfixn *A2i = A2 + ci;
      mpzToPolyTwoMod((coeff[i]).get_mpz_t(), M, p1, p2, Ai, A2i);
  }

  return biA;
}

//---------------------------------------------------
//just for binary representation mpz
//precondition: z!=0, number of bits of z / M <= K
//K = size of X1 
void MulSSA::mpzToPolyTwoMod(mpz_t z, int M, 
			     sfixn prime1, sfixn prime2,
			     sfixn *X1, sfixn *X2) 
{//----------------------------------------------------------  
  int s = z->_mp_size; //abs(s)>=1
  int isNeg = 0;
  if (s<0) {
    isNeg = 1;
    s = -s;
  }

  mp_ptr limbs = z->_mp_d;
  mp_limb_t carry = 0;
  mp_size_t carry_size = 0;

  int e = s-1; 
  int index = 0;

  //all the limbs except for the last one are full-size (64 bits)
  size_t b = GMP_LIMB_BITS;
  
  int i = 0;
  while (i<e){//1st to the 2nd last limb
    long unsigned int w = limbs[i];
    
    i++;  
    size_t j = 0; //current position in w
    if (carry_size == 0){
      while (b-j>=M){
	long unsigned int q = w>>M;
	X1[index] = (sfixn) ((w) ^ (q<<M));	
	if (isNeg==1){
	  X1[index] = -X1[index];
	  X2[index] = X1[index]+prime2;
	  X1[index] = X1[index]+prime1;	  
	}else
	  X2[index] = X1[index];

	index++;
	w = q;
	j += M; 
      }
    }else{
      size_t a = M - carry_size;    
      while (b-j>=a){ 	
	long unsigned int q = (w>>a); //quotient = w / 2^a
	X1[index] = (sfixn) ((((w) ^ (q<<a)) <<carry_size) + carry); //w % 2^a + carry, 
	if (isNeg==1){
	  X1[index] = -X1[index];
	  X2[index] = X1[index]+prime2;
	  X1[index] = X1[index]+prime1;	  
	}else
	  X2[index] = X1[index];
	
	w = q;
	j += a;
	index++;
	
	a = M;
	carry_size = 0;
	carry = 0;
      }  
    }
    carry_size = b - j;
    carry = w;  
  }
  
  //last limb, its bit size may be smaller than long
  long unsigned int w = limbs[e];
  b = log2(w)+1; //num of bits
  
  size_t j = 0; //current position in w    
  if (carry_size == 0){
    while (b-j>=M){
      long unsigned int q = w>>M;  
 
      X1[index] = (sfixn) ((w) ^ (q<<M));
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];

      index++;
      w = q;
      j += M; 
    }
    if (w!=0) {//carry
      X1[index] = (sfixn) w;
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];
    }
 
  }else{
    size_t a = M - carry_size;  
    while (b-j>=a){ 
      long unsigned int q = (w>>a); //quotient = w / 2^a
      X1[index] = (sfixn) ((((w) ^ (q<<a)) <<carry_size) + carry); //w % 2^a + carry, 
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];

      w = q;
      j += a;
      index++;
      
      a = M;
      carry_size = 0;
      carry = 0;
    } 

    if (w!=0){
      X1[index] = (sfixn) ((w <<carry_size) + carry); 
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];
    }
  } 
}

int MulSSA::binarySearch(int n, int start, int end) {
	int mid = (start + end) / 2;
  
	if (n == goodN[mid])
		return mid;
	else if (n < goodN[mid] && mid > 0 && n > goodN[mid - 1]) 
		return mid;
	else if (n > goodN[mid] && n <= goodN[end])
		return binarySearch(n, mid, end);
	else if (n < goodN[mid] && n >= goodN[start]) 
		return binarySearch(n, start, mid);
	else {
		std::cout<<"Cannot find good N, K, M.\n";
		exit(1);
	}	  
}

void MulSSA::determineGoodN(int n, sfixn d) {
	int index = binarySearch(n, 0, size-1); 
	N = goodN[index];
	K = 1 << goodK[index];
	M = goodM[index];	
	LIMB_BITS = (1 + goodK[index] + logceiling(d) + M + N) + 64;
}
