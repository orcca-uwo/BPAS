
/**
 * Implementation of 'MulSSA.h'
 *
 * @author Yuzhen Xie
 */

#if TDEBUG
#include <cilktools/cilkview.h>
#endif
#include "../../../include/IntegerPolynomial/Multiplication/MulSSA.h"

/*----------------------------------------
Our new 2-convolution multiplication
-----------------------------------------*/
void MulSSA::mul2C(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) 
{//----------------------------------------
	cilk_spawn a->setCoefficientDigits(a->maxCoefficientSize());
	b->setCoefficientDigits(b->maxCoefficientSize());
	cilk_sync;

	int digitCount = MAX_COEFF_DIGITS(a, b);
	determineGoodN(digitCount, MAX(a->getSize(),b->getSize()));
//std::cout << "N: " << N << ", K: " << K << ", M: " << M << std::endl;

	//mul2C1(a, b, c);
	mul2C2(a, b, c);
}

void MulSSA::mul2C1(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) 
{//--------------------------------------------  
	int d1 = a->getSize();
	int d2 = b->getSize();

	int csize = d1 + d2 ;
	int es2 = logceiling(csize);
	int dims2 = 1<<es2;

	BivariatePolynomial *biA = cilk_spawn ToBivarMod(a->getCoefficients(), d1, P);
	BivariatePolynomial *biB = ToBivarMod(b->getCoefficients(), d2, P);
	cilk_sync;

	sfixn *ncc = NULL;
	
	//sfixn *
	if (csize==dims2)
	  ncc = PBPAS::TwoConvolutionModNew(biA->getCoefficients(), biB->getCoefficients(), d1, d2, K, MontP, DFTBASESIZE, RevBidMap,0);
	  else
	    ncc = PBPAS::TwoConvolutionMod(biA->getCoefficients(), biB->getCoefficients(), d1, d2, P, K);
	  

	biA->freeHeap();
	biB->freeHeap();
	csize--;
	RecoverProduct(ncc, csize, dims2, c->getCoefficients());

	delete [] ncc;
}

BivariatePolynomial* MulSSA::ToBivarMod(mpz_class* coef, int d, sfixn p) {
	BivariatePolynomial* biA = new BivariatePolynomial(d, K, M);

	#pragma cilk_grainsize = 1024;
	cilk_for (int i = 0; i < d; ++i) {
		sfixn* a = biA->getCoefficients() + i * K;
		mpzToPolyMod(a, coef[i].get_mpz_t(), M, p);
	}

	return biA;
}

void MulSSA::mpzToPolyMod(sfixn* a, mpz_t coef, int M, sfixn p) {
	if (coef->_mp_size == 0) {
		for (int i = 0; i < K; ++i)
			a[i] = 0;
		return;
	}
        bool isNeg = 0;
        int s = coef->_mp_size;
        if (s < 0) {
                isNeg = 1;
                s = -s;
        }

        int idx = 0, carry_size = GMP_LIMB_BITS;
        mp_limb_t carry = coef->_mp_d[idx];
        for (int i = 0; i < K; ++i) {
                int k = M;
                mp_limb_t w = carry;
                if (carry_size < M) {
                        idx++;
                        k = M - carry_size;
                        if (idx < s) {
                                carry = coef->_mp_d[idx];
                                w += carry << carry_size;
                                carry_size = GMP_LIMB_BITS;
                        }
                        else {
                                carry = 0;
                                carry_size = 0;
                                a[i] = (sfixn) w;
                        }
                }

                if (carry_size >= M) {
                        a[i] = (sfixn) (w ^ ((w >> M) << M));
                        carry >>= k;
                        carry_size -= k;
                }

                if (isNeg)
                        a[i] = p - a[i];
        }
}

void MulSSA::RecoverProduct(sfixn* ncc, int csize, int dims2, mpz_class* c) {
  sfixn* cc = ncc;
  if (csize+1!=dims2)
    cc += K * csize;
  else
    cc += K * dims2;
 
        #pragma cilk_grainsize = 1024;
        cilk_for (int i = 0; i < csize; ++i) {
                int Ks = i * K;
                mpz_t u, v;
                ToMPZ(u, ncc+Ks);
                ToMPZ(v, cc+Ks);

                mpz_add (c[i].get_mpz_t(), u, v);
                c[i] >>= 1;
                mpz_sub (v, v, u);
                c[i] += mpz_class (v) << (N - 1);

                mpz_clear(u);
                mpz_clear(v);
        }
}

void MulSSA::ToMPZ(mpz_t zp, sfixn* a) {
	mpz_t zn;
	mpz_init2(zp, LIMB_BITS);
	mpz_init2(zn, LIMB_BITS);

	int idx = 0;
	__int128 Tp = 0, Tn = 0;
	sfixn shifter = 0;
	for (int i = 0; i < K; ++i) {
		if (a[i] > HALF_P)
			Tn += (__int128) (P - a[i]) << shifter;
		else if (a[i] > 0)
			Tp += (__int128) a[i] << shifter;
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

void MulSSA::mul2C2(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) 
{//------------------------------------------------
  //std::cout<<"mul2C2, K: "<<K<<std::endl;
  //if (MontP1==NULL) 
  //  std::cout<<"MontP1 null \n ";
  //else{
  //  std::cout<<"MontP1: "<<MontP1->P<<std::endl;
  //  std::cout<<"mul2C2, before a->getSize \n";
  //}

 
  //	std::cout<<"mul2C2, before a->getSize 2\n";
	int d1 = a->getSize();
	int d2 = b->getSize();
	#if TDEBUG
	unsigned long time;
	#endif
	//assert(K>=DFTBASESIZE);
	//assert(d1+d2>=DFTBASESIZE);
	#if TDEBUG
	std::cout<<"MUL2C2"<<std::endl;
	std::cout<<"d1 " <<d1<<" d2 "<<d2<<std::endl;
	#endif

	int csize = d1 + d2 ;
	int es2 = logceiling(csize);
	int dims2 = 1<<es2;
	
	
	//std::cout<<"mul2C2, before ToBivarMod \n";
//mpz_class *gmpA = a->getCoefficients();
//for (int i = 0; i < d1; ++i)
//std::cout << gmpA[i] << " ";
//std::cout << std::endl;
	
	BivariatePolynomial* biA = cilk_spawn ToBivarMod(a->getCoefficients(), d1, P1, P2);
	BivariatePolynomial* biB = ToBivarMod(b->getCoefficients(), d2, P1, P2);
	cilk_sync;

//std::cout<<"mul2C2, after ToBivarMod \n";
//sfixn *tmpA = biA->getCoefficients();
//for (int i = 0; i < d1; ++i)
//std::cout << tmpA[i] << " ";
//std::cout << std::endl;

	sfixn *nc1 = NULL;
	sfixn *nc2 = NULL;
	if (csize==dims2) {
	#if TDEBUG
	std::cout << "Entering New TwoConvolution"<<std::endl;
	time = -__cilkview_getticks();
	#endif
	  nc1 = cilk_spawn PBPAS::TwoConvolutionModNew(biA->getCoefficients(), biB->getCoefficients(), d1, d2, K, MontP1, DFTBASESIZE, NULL,1);
	  nc2 = PBPAS::TwoConvolutionModNew(biA->getCoefficients()+biA->getSize(), biB->getCoefficients()+biB->getSize(), d1, d2, K, MontP2, DFTBASESIZE, NULL,2);

	  cilk_sync;
	#if TDEBUG
	time += __cilkview_getticks();
	std::cout << "Two convolution time :" <<time<<std::endl;
	#endif
	}else{
	#if TDEBUG
	std::cout << "Entering Old TwoConvolution"<<std::endl;
	unsigned long time = -__cilkview_getticks();
	#endif
	  nc1 = cilk_spawn PBPAS::TwoConvolutionMod(biA->getCoefficients(), biB->getCoefficients(), d1, d2, P1, K);
	  nc2 = PBPAS::TwoConvolutionMod(biA->getCoefficients()+biA->getSize(), biB->getCoefficients()+biB->getSize(), d1, d2, P2, K);
	  cilk_sync;
	#if TDEBUG
	time += __cilkview_getticks();
	std::cout << "Two convolution time :" <<time<<std::endl;
	#endif
	}
	//std::cout<<"mul2C2, after TwoConvolutionModNew \n";

	biA->freeHeap();
	biB->freeHeap();
	//std::cout<<"mul2C2, after biB->freeHeap \n";
	csize--;
	#if TDEBUG
	time = -__cilkview_getticks();
	std::cout<< M <<std::endl;
	#endif
	RecoverProduct(nc1, nc2, csize, dims2, c->getCoefficients());
	#if TDEBUG
	time += __cilkview_getticks();
	std::cout << "RecoverProduct time :" <<time<<std::endl;
	#endif
	my_free(nc1);
	my_free(nc2);
	//#if DEBUG
	//std::cout<<"free done"<<std::endl;
	//#endif
	//std::cout<<"mul2C2, after delete [] nc2 \n";
}

BivariatePolynomial* MulSSA::ToBivarMod(mpz_class* coef, int d, sfixn p1, sfixn p2) 
{//-------------------------------------

	BivariatePolynomial* biA = new BivariatePolynomial(d, K, M, 2);

	#pragma cilk_grainsize = 1024;
	cilk_for (int i = 0; i < d; ++i) {
		int base = i * K, base1 = biA->getSize() + base;
		sfixn* a1 = biA->getCoefficients() + base;
		sfixn* a2 = biA->getCoefficients() + base1;
		mpzToPolyMod(a1, a2, coef[i].get_mpz_t(), M, p1, p2);
	}
	return biA;
}

void MulSSA::mpzToPolyMod(sfixn* a1, sfixn* a2, mpz_t coef, int M, sfixn p1, sfixn p2) {
	if (coef->_mp_size == 0) {
		for (int i = 0; i < K; ++i)
			a2[i] = a1[i] = 0;
			
		return;
	}

	bool isNeg = 0;
	int s = coef->_mp_size;
	if (s < 0) {
		isNeg = 1;
		s = -s;
	}

	int idx = 0, carry_size = GMP_LIMB_BITS;
	mp_limb_t carry = coef->_mp_d[idx];
	for (int i = 0; i < K; ++i) {
		int k = M;
		mp_limb_t w = carry;
		if (carry_size < M) {
			idx++;
			k = M - carry_size;
			if (idx < s) {
				carry = coef->_mp_d[idx];
				w += carry << carry_size;
				carry_size = GMP_LIMB_BITS;
			}
			else {
				carry = 0;
				carry_size = 0;
				a1[i] = (sfixn) w;
			}
		}
		
		if (carry_size >= M) {
			a1[i] = (sfixn) (w ^ ((w >> M) << M));
			carry >>= k;
			carry_size -= k;
		}

		if (isNeg) {
			a1[i] = -a1[i];
			a2[i] = a1[i] + p2;
			a1[i] = a1[i] + p1;
		}
		else { a2[i] = a1[i]; }
	}
}

void MulSSA::RecoverProduct(sfixn* nc1, sfixn* nc2, int csize, int dims2, mpz_class* c) {

	int fullsize = K * dims2;
	if (csize+1!=dims2)
	  fullsize = K * csize; //for old TwoConvolutionMod
	
	sfixn* ncc1 = nc1;
	sfixn* cc1 = nc1 + fullsize;
	sfixn* ncc2 = nc2;
	sfixn* cc2 = nc2 + fullsize;

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

void MulSSA::CRTtoMPZ(mpz_t zp, sfixn* a, sfixn* b, int size) {
	mpz_t zn;
	mpz_init2(zp, LIMB_BITS);
	mpz_init2(zn, LIMB_BITS);

	int idx = 0;
	unsigned __int128 postainer[2] = {0}, negtainer[2] = {0};
	sfixn shifter = 0;
	for (int i = 0; i < size; ++i) {
		//sfixn diff =(a[i] - b[i]);
		//elem = elem * P2 + b[i];
		sfixn diff = (a[i] - b[i]);
		if(diff<0){
			diff+=P1;
		}
		__int128 elem = MontMulModSpe_OPT3_AS_GENE_globalfunc(diff,U2_R1_sft,INV_PRIME1,MY_PRIME1);
		elem = elem * P2 + b[i];
		if (elem > HALF_P1_P2) { elem -= P1_P2; }
		else if (elem < N_HALF_P1_P2) { elem += P1_P2; }

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
		std::cout << "BPAS: error, cannot find good N, K, M." << std::endl;;
		exit(1);
	}
}

void MulSSA::determineGoodN(int n, sfixn d) {
	int index = binarySearch(n, 0, size-1); 
	N = goodN[index];
	K = 1 << goodK[index];
	M = goodM[index];
	//u_i, v_i has bit bound 2dK2^(M+N) 
	//+1 for the extra limb recommended by mpz_t
	LIMB_BITS = (1 + goodK[index] + logceiling(d) + M + N) + 128; //tight bound
}
