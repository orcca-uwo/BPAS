


#include "IntegerPolynomial/DUZP_Hensel.h"
#include "Utils/MacroHelpers.h"

/**
 * Solve the two-term diophantine equation u*sigma + w*tau = c
 * for sigma nad tau. 
 *
 * sigma and tau are returned by pointer.
 *
 * returns 1 iff the equation can be solved (i.e. gcd(u,w) does not divide c )
 */
int UDP_spX(const duspoly_t* u, const duspoly_t* w, const duspoly_t* c, duspoly_t** sigma, duspoly_t** tau, const Prime_ptr* Pptr) {
	//Special case down to just EEA
	if (isOneInForm_spX(c, Pptr)) {
		duspoly_t* g = NULL;
		extGCDInForm_spX(u, w, sigma, tau, &g, Pptr);
		// fprintf(stderr, "u:\n");
		// printPolynomialOutForm_spX(u, Pptr);
		// fprintf(stderr, "w:\n");
		// printPolynomialOutForm_spX(w, Pptr);
		// fprintf(stderr, "sigma:\n");
		// printPolynomialOutForm_spX(*sigma, Pptr);
		// fprintf(stderr, "tau:\n");
		// printPolynomialOutForm_spX(*tau, Pptr);


		if (!isOneInForm_spX(g, Pptr)) {
			printPolynomialOutForm_spX(g, Pptr);
			freePolynomial_spX (&g);
			freePolynomial_spX (&*sigma);
			freePolynomial_spX (&*tau);
			return 0;
		}
		freePolynomial_spX (&g);
		return 1;
	}

	//Otherwise do full solve.
	duspoly_t* g = NULL;
	duspoly_t* s = NULL;
	duspoly_t* t = NULL;
	extGCDInForm_spX(u, w, &s, &t, &g, Pptr);

	//simplified problem where we don't need to divide by g.
	if (isOneInForm_spX(g, Pptr)) {
		//sigma = s * c
		//tau = t * c
		mulPolynomialsInForm_spX (s, c, sigma, Pptr);
		mulPolynomialsInForm_spX (t, c, tau, Pptr);
		
		//free in prep for re-using variable as quotient or product
		freePolynomial_spX (&t);
		freePolynomial_spX (&s);

		//sigma = rem(sigma, w); s = quo(sgima, w)
		plainDivPolynomialsInForm_spX_inp(sigma, w, &s, Pptr);

		//tau = t*c + quo(sigma,w)*u
		mulPolynomialsInForm_spX (s, u, &t, Pptr);
		addPolynomialsInForm_spX_inp (tau, t, Pptr);

		freePolynomial_spX (&s);
		freePolynomial_spX (&t);
		return 1;
	}

	//full solution to problem
	if (!isDividablePolys_spX(c, g, Pptr)) {
		freePolynomial_spX (&g);
		freePolynomial_spX (&s);
		freePolynomial_spX (&t);
		return 0;		
	}

	//s = s*c/g, t = t*c/g
	mulPolynomialsInForm_spX (s, c, sigma, Pptr);
	mulPolynomialsInForm_spX (t, c, tau, Pptr);
	freePolynomial_spX (&s);
	freePolynomial_spX (&t);
	plainExactDivPolynomialsInForm_spX(*sigma, g, &s, Pptr);
	plainExactDivPolynomialsInForm_spX(*tau, g, &t, Pptr);
	freePolynomial_spX (&*tau);

	//sigma = rem(s, w/g), s = quo(s, w/g) 
	duspoly_t* wgQuo;
	plainExactDivPolynomialsInForm_spX(w, g, &wgQuo, Pptr);
	plainDivPolynomialsInForm_spX_inp(&s, wgQuo, sigma, Pptr);
	duspoly_t* temp = NULL;
	temp = s;
	s = *sigma;
	*sigma = temp;

	//tau = t*c/g + quo(s, w/g)*(u/g)
	duspoly_t* ugQuo;
	plainExactDivPolynomialsInForm_spX(u, g, &ugQuo, Pptr);
	mulPolynomialsInForm_spX (s, ugQuo, tau, Pptr);
	addPolynomialsInForm_spX_inp(tau, t, Pptr);

	freePolynomial_spX (&s);
	freePolynomial_spX (&t);
	freePolynomial_spX (&ugQuo);
	freePolynomial_spX (&wgQuo);
	return 1;
}


/**
 * Solve the multi-term univariate diophantine problem for sigma_i in 
 * b[0]*sigma[0] + ... + b[r-1]*sigma[r-1] = c mod p 
 * where b[j] = prod f[i] for i != j.
 * Each sigma_i returned in the pre-allocated array sigmas.
 *
 * c : the rhs of the diophantine equation
 * fs : an array of the factors of the b polynomials in the diophantine equation
 *
 * returns 1 iff the problem could be solved (i.e. if gcd(f[i],f[j]) == 1 for all i,j
 * and p does not divide l.c. of f[i] for all i).
 */
int multiUDP_spX(const duspoly_t* c, duspoly_t const*const* fs, int r, duspoly_t** sigmas, const Prime_ptr* Pptr) {

	if (r == 2) {
		return UDP_spX(fs[1], fs[0], c, sigmas, sigmas + 1, Pptr);
	}

	/* First solve b[0]*s[0] + ... + b[r-1]*s[i-1] = 1 for each s*/

	//compute the products of f, beta_i = prod(f_j), j >= i+1 
	//temporarily use sigmas to store products of f 
	//for j = r-1 and r-2 just use f[r-1], f[r-2]
	mulPolynomialsInForm_spX (fs[r-1], fs[r-2], sigmas + r-2, Pptr);
	for (int i = r-3; i > 0; --i) {
		mulPolynomialsInForm_spX (fs[i], sigmas[i+1], sigmas + i, Pptr);
	}

	duspoly_t* beta = makePolynomial_spX(1);
	duspoly_t* nextBeta = NULL;
    beta->elems[0] = smallprimefield_convert_in (1, Pptr);
	for (int j = 0; j < r-2; ++j) {
		if (!UDP_spX(fs[j], sigmas[j+1], beta, &nextBeta, sigmas + j, Pptr)) {
			fprintf(stderr, "UDP failed! fs are not co-prime!!!!!!!!\n\n\n\n");
			exit(1);
		}
		freePolynomial_spX (&sigmas[j+1]);
		freePolynomial_spX (&beta);
		beta = nextBeta;
		nextBeta = NULL;
	}
	//do last solve manually using f[r-2], f[r-1]
	// fprintf(stderr, "about to do manual solve in multiUDP\n");
	if (!UDP_spX(fs[r-2], fs[r-1], beta, sigmas + r-1, sigmas + r-2, Pptr)) {
		fprintf(stderr, "UDP failed! fs are not co-prime!!!!!!!!\n\n\n\n");
		exit(1);
	}
	// fprintf(stderr, "got UDP solution1: \n");
	// printPolynomialOutForm_spX(sigmas[r-1], Pptr);
	// fprintf(stderr, "got UDP solution2: \n" );
	// printPolynomialOutForm_spX(sigmas[r-2], Pptr);
	freePolynomial_spX (&beta);
	beta = NULL;

	if (!isOneInForm_spX(c, Pptr)) {
		//Multiply through by c now and get rem by fs[i] to get sigma[i].
		for (int j = 0; j < r; ++j) {
			mulPolynomialsInForm_spX (c, sigmas[j], &beta, Pptr);
			freePolynomial_spX (&sigmas[j]);
			// fprintf(stderr, "c*sigmas[%d]\n", j);
			// printPolynomialOutForm_spX(beta, Pptr);
			plainRemPolynomialsInForm_spX_inp(&beta, fs[j], Pptr);
			// fprintf(stderr, "rem beta:[%d]\n", j);
			// printPolynomialOutForm_spX(beta, Pptr);
			sigmas[j] = beta;
			beta = NULL;
		}
	}

	return 1;
}


// a += p * mod
void padicUpdate(DUZP_t* a, duspoly_t* p, mpz_t modulus, const Prime_ptr* Pptr) {
	if (p == NULL || p->alloc == 0) {
		return;
	}

	if (a->lt < p->lt) {
		//this shouldn't happen?
		fprintf(stderr, "nope in padicUpdate\n");

		fprintf(stderr, "a: \n");
		printPoly_DUZP(a, "x");

		fprintf(stderr, "p: \n");
		printPolynomialOutForm_spX(p, Pptr);

		exit(1);
	}

	for (int i = 0; i <= p->lt; ++i) {
		mpz_addmul_ui(a->coefs[i], modulus, (unsigned long long) smallprimefield_convert_out(p->elems[i], Pptr));
	}

}


// a += p * mod
void symmetricPadicUpdate(DUZP_t* a, duspoly_t* p, mpz_t modulus, const Prime_ptr* Pptr, register unsigned long long halfP, mpz_t workZ) {
	if (p == NULL || p->alloc == 0) {
		return;
	}

	if (a->lt <= p->lt) {
		//this shouldn't happen?
		fprintf(stderr, "nope in symmetric padicUpdate\n");

		fprintf(stderr, "a: \n");
		printPoly_DUZP(a, "x");

		fprintf(stderr, "p: \n");
		printPolynomialOutForm_spX(p, Pptr);

		exit(1);
	}


	register unsigned long long pelem;
// 1
	register long long prime = Pptr->prime;
// 2	
	// register unsigned long long prime = Pptr->prime;
	for (int i = 0; i <= p->lt; ++i) {
		pelem = (unsigned long long) smallprimefield_convert_out(p->elems[i], Pptr);
//1
		if (pelem > halfP) {
			mpz_set_si(workZ, (long long) pelem - prime);
			mpz_addmul(a->coefs[i], modulus, workZ);
		} else {
			mpz_addmul_ui(a->coefs[i], modulus, pelem);
		}
//2
		// if (pelem > halfP) {
		// 	mpz_addmul_ui(a->coefs[i], modulus, pelem);
		// 	mpz_submul_ui(a->coefs[i], modulus, prime);
		// } else {
		// 	mpz_addmul_ui(a->coefs[i], modulus, pelem);
		// }

	}
}

// a += p * mod
void symmetricQuadraticPadicUpdate(DUZP_t* a, DUZP_t* p, mpz_t modulus, mpz_t halfMod) {
	if (p == NULL || p->alloc == 0) {
		return;
	}

	if (a->lt <= p->lt) {
		//this shouldn't happen?
		fprintf(stderr, "nope in symmetric quadratic padicUpdate\n");

		fprintf(stderr, "a: \n");
		printPoly_DUZP(a, "x");

		fprintf(stderr, "p: \n");
		printPoly_DUZP(p, "x");

		exit(1);
	}

	for (int i = 0; i <= p->lt; ++i) {
//1
		if (mpz_cmp(p->coefs[i], halfMod) > 0) {
			mpz_sub(p->coefs[i], p->coefs[i], modulus);
		}
		mpz_addmul(a->coefs[i], p->coefs[i], modulus);
//2
		// if (pelem > halfP) {
		// 	mpz_addmul_ui(a->coefs[i], modulus, pelem);
		// 	mpz_submul_ui(a->coefs[i], modulus, prime);
		// } else {
		// 	mpz_addmul_ui(a->coefs[i], modulus, pelem);
		// }

	}

}

// a += p * mod
void quadraticPadicUpdate(DUZP_t* a, const DUZP_t* p, const mpz_t modulus) {
	if (p == NULL || p->alloc == 0) {
		return;
	}

	if (a->lt < p->lt) {
		//this shouldn't happen?
		fprintf(stderr, "nope in quadratic padicUpdate\n");

		fprintf(stderr, "a: \n");
		printPoly_DUZP(a, "x");

		fprintf(stderr, "p: \n");
		printPoly_DUZP(p, "x");

		exit(1);
	}

	for (int i = 0; i <= p->lt; ++i) {
		mpz_addmul(a->coefs[i], p->coefs[i], modulus);
	}

}


int monicPositiveTwoTermPadicLift(const DUZP_t* a, const duspoly_t* u, const duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const Prime_ptr* Pptr) {
	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);


	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	int primePow = 1;


	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;

	extGCDInForm_spX(u, w, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX (&sp_one);

	DUZP_t* uu = convertFromDUSP_DUZP(u, Pptr);
	DUZP_t* ww = convertFromDUSP_DUZP(w, Pptr);
	DUZP_t* error = multiplyPolynomials_DUZP(uu, ww);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);

		//sigma = s * c
		//tau = t * c
		duspoly_t* sigma = multiplyAndModDUSP_DUZP(error, sp_s, Pptr);
		duspoly_t* tau = multiplyAndModDUSP_DUZP(error, sp_t, Pptr);

		//sigma = rem(sigma, w); s = quo(sgima, w)
		duspoly_t* q;
		plainDivPolynomialsInForm_spX_inp(&sigma, w, &q, Pptr);
		//in place or fast div???
		// duspoly_t* new_sigma;
		// fastDivPolynomialInForm_wPSInv_spX (sigma, w, &new_sigma, &q, Pptr);


		//tau = t*c + quo(sigma,w)*u
		duspoly_t* quProd;
		mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
		freePolynomial_spX (&quProd);

		padicUpdate(uu, tau, mod, Pptr);
		padicUpdate(ww, sigma, mod, Pptr);

		//update the error.
		multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}

	mpz_clears(bound, mod, NULL);

	freePolynomial_DUZP(error);
	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);

	if (liftedU != NULL) {
		*liftedU = uu;
	} else {
		freePolynomial_DUZP(uu);
	}

	if (liftedW != NULL) {
		*liftedW = ww;
	} else {
		freePolynomial_DUZP(ww);
	}

	return primePow;
}

/**
 * This one uses symmetric range of finite field in case any coef in a is negative.
 */
int monicTwoTermPadicLift(const DUZP_t* a, const duspoly_t* u, const duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const Prime_ptr* Pptr) {
	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);
	int primePow = 1;

	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;

	extGCDInForm_spX(u, w, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		printPolynomialOutForm_spX(u, Pptr);
		printPolynomialOutForm_spX(w, Pptr);

		exit(1);
	}
	freePolynomial_spX (&sp_one);

	DUZP_t* uu = convertFromDUSPSymmetric_DUZP(u, Pptr);
	DUZP_t* ww = convertFromDUSPSymmetric_DUZP(w, Pptr);
	DUZP_t* error = multiplyPolynomials_DUZP(uu, ww);

	subtractPolynomials_DUZP_inpRHS(a, &error);

	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {

		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);

		//sigma = s * c
		//tau = t * c
		duspoly_t* sigma = multiplyAndModDUSP_DUZP(error, sp_s, Pptr);
		duspoly_t* tau = multiplyAndModDUSP_DUZP(error, sp_t, Pptr);

		//sigma = rem(sigma, w); s = quo(sgima, w)
		duspoly_t* q;
		plainDivPolynomialsInForm_spX_inp(&sigma, w, &q, Pptr);
		//TODO 
		//in place or fast div???
		// duspoly_t* new_sigma;
		// fastDivPolynomialInForm_wPSInv_spX (sigma, w, &new_sigma, &q, Pptr);

		//tau = t*c + quo(sigma,w)*u
		duspoly_t* quProd;
		mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
		freePolynomial_spX (&quProd);

		symmetricPadicUpdate(uu, tau, mod, Pptr, halfP, workZ);
		symmetricPadicUpdate(ww, sigma, mod, Pptr, halfP, workZ);

		//update the error.
		multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}


	mpz_clears(bound, mod, workZ, NULL);

	freePolynomial_DUZP(error);

	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);

	if (liftedU != NULL) {
		*liftedU = uu;
	} else {
		freePolynomial_DUZP(uu);
	}

	if (liftedW != NULL) {
		*liftedW = ww;
	} else {
		freePolynomial_DUZP(ww);
	}

	return primePow;

}


/**
 *
 * u and w are modified throughout this algorithm to save space.
 */ 
int twoTermPadicLift(const DUZP_t* aa, duspoly_t* u, duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const mpz_t gamma, const Prime_ptr* Pptr) {
	if (liftedU == NULL && liftedW == NULL) {
		return 0;;
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(aa, bound);
	mpz_mul_si(bound, bound, 2l);
	mpz_mul(bound, bound, gamma);
	mpz_abs(bound, bound);

	mpz_t alpha;
	mpz_init(alpha);
	mpz_set(alpha, aa->coefs[aa->lt]);


//Two strategies here for leading coefficient correct. 
//First (commented) is based off Algorithms for Computer Algebra. 
//One sets leading cofficients of both u and w so that their product 
//has the leading coefficient of a at every step. Hence deg(error) < deg(a).
//Second is similar to ensure deg(error) < deg(a). Force w to be monic
//and force lc(u) = lc(a). This second technique is used in quadratic lifting
//and multi-term lifting.

//1
	// DUZP_t* a = multiplyByInteger_DUZP(aa, gamma);
	// prime_t gamma_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);
	// if (smallprimefield_convert_out(u->elems[u->lt], Pptr) != 1) {
	// 	gamma_sp = smallprimefield_mul(gamma_sp, smallprimefield_inv(u->elems[u->lt], Pptr), Pptr);
	// }
	// scalarMulPolynomialInForm_spX_inp (&u, gamma_sp, Pptr);

	// prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);	
	// if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
	// 	alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	// }
	// scalarMulPolynomialInForm_spX_inp (&w, alpha_sp, Pptr);
	
//2
	DUZP_t* a = (DUZP_t*) aa;
	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);
	if (smallprimefield_convert_out(u->elems[u->lt], Pptr) != 1) {
		alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(u->elems[u->lt], Pptr), Pptr);
	}
	scalarMulPolynomialInForm_spX_inp (&u, alpha_sp, Pptr);

	if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
		scalarMulPolynomialInForm_spX_inp (&w, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	}
	
	// fprintf(stderr, "u after lc adjust:\n");
	// printPolynomialOutForm_spX(u, Pptr);
	// fprintf(stderr, "w after lc adjust:\n");
	// printPolynomialOutForm_spX(w, Pptr);

	//get bezout coefficients
	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;
	extGCDInForm_spX(u, w, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX (&sp_one);

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);

	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);
	int primePow = 1;

	DUZP_t* uu = convertFromDUSPSymmetric_DUZP(u, Pptr);
	DUZP_t* ww = convertFromDUSPSymmetric_DUZP(w, Pptr);
	
//1
	// mpz_set(uu->coefs[uu->lt], gamma);
	// mpz_set(ww->coefs[ww->lt], alpha);
//2
	mpz_set(uu->coefs[uu->lt], alpha);
	
	DUZP_t* error = multiplyPolynomials_DUZP(uu, ww);
	subtractPolynomials_DUZP_inpRHS(a, &error);
	// fprintf(stderr, "error:\n");
	// printPoly_DUZP(error, "x");

	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {

		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);

		//sigma = s * c
		//tau = t * c
		duspoly_t* sigma = multiplyAndModDUSP_DUZP(error, sp_s, Pptr);
		duspoly_t* tau = multiplyAndModDUSP_DUZP(error, sp_t, Pptr);

		//sigma = rem(sigma, w); s = quo(sgima, w)
		duspoly_t* q;
		plainDivPolynomialsInForm_spX_inp(&sigma, w, &q, Pptr);
		//in place or fast div???
		// duspoly_t* new_sigma;
		// fastDivPolynomialInForm_wPSInv_spX (sigma, w, &new_sigma, &q, Pptr);

		//tau = t*c + quo(sigma,w)*u
		duspoly_t* quProd;
		// mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
		freePolynomial_spX (&quProd);

		symmetricPadicUpdate(uu, tau, mod, Pptr, halfP, workZ);
		symmetricPadicUpdate(ww, sigma, mod, Pptr, halfP, workZ);

		//update the error.
		multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}

	//make lifted polys primitive
	mpz_t cont;
	mpz_init(cont);
	primitivePartAndContent_DUZP_inp(uu, cont);
	
//1
	// freePolynomial_DUZP(a);
//2
	//nothing

	freePolynomial_DUZP(error);
	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);

	if (liftedU != NULL) {
		*liftedU = uu;
	} else {
		freePolynomial_DUZP(uu);
	}

	if (liftedW != NULL) {
		mpz_divexact(cont, gamma, cont);
		if (mpz_cmp_si(cont, 1l) != 0) {
			divideByIntegerExact_DUZP_inp(ww, cont);
		}
		*liftedW = ww;
	} else {
		freePolynomial_DUZP(ww);
	}

	mpz_clears(bound, alpha, mod, cont, workZ, NULL);
	return primePow;
}

/**
 * A two-term padic lift specialized for GCD computation. 
 * aa is one of the polynomials whose gcd is being computed,
 * u is the gcd image, w is the cofactor. 
 * It is assumed that the gcd image is monic in the finite field. 
 * Hence, lc(w) = lc(aa) mod Pptr->prime.
 */
int gcdPadicLift(const DUZP_t* aa, duspoly_t* u, duspoly_t* w, DUZP_t** liftedU, const Prime_ptr* Pptr) {
	if (liftedU == NULL) {
		return 0;
	}

//Two strategies here for leading coefficient correct. 
//First (commented) is based off Algorithms for Computer Algebra. 
//One sets leading cofficients of both u and w so that their product 
//has the leading coefficient of a at every step. Hence deg(error) < deg(a).
//Second is similar to ensure deg(error) < deg(a). Force w to be monic
//and force lc(u) = lc(a). This second technique is used in quadratic lifting
//and multi-term lifting. It won't necessarily get the true lifted cofactor...


	mpz_t alpha;
	mpz_init(alpha);
	mpz_set(alpha, aa->coefs[aa->lt]);

//1
	// DUZP_t* a = multiplyByInteger_DUZP(aa, gamma);
	// prime_t gamma_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);
	// if (smallprimefield_convert_out(u->elems[u->lt], Pptr) != 1) {
	// 	gamma_sp = smallprimefield_mul(gamma_sp, smallprimefield_inv(u->elems[u->lt], Pptr), Pptr);
	// }
	// scalarMulPolynomialInForm_spX_inp (&u, gamma_sp, Pptr);

	// prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);	
	// if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
	// 	alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	// }
	// scalarMulPolynomialInForm_spX_inp (&w, alpha_sp, Pptr);
	
//2


	DUZP_t* a = (DUZP_t*) aa;
	//make w monic and give u the lc(a) mod p.
	if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
		scalarMulPolynomialInForm_spX_inp (&u, w->elems[w->lt], Pptr);
		scalarMulPolynomialInForm_spX_inp (&w, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	}
	
	//get bezout coefficients
	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;
	extGCDInForm_spX(u, w, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX (&sp_one);

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);

	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);

	int primePow = 1;


	DUZP_t* uu = convertFromDUSPSymmetric_DUZP(u, Pptr);
	DUZP_t* ww = convertFromDUSPSymmetric_DUZP(w, Pptr);
	
//1
	// mpz_set(uu->coefs[uu->lt], gamma);
	// mpz_set(ww->coefs[ww->lt], alpha);
//2
	mpz_set(uu->coefs[uu->lt], alpha);

	mpz_abs(alpha, alpha);
	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(aa, bound);
	mpz_mul_si(bound, bound, 2l);
	mpz_mul(bound, bound, alpha);
	
	DUZP_t* error = multiplyPolynomials_DUZP(uu, ww);
	
	subtractPolynomials_DUZP_inpRHS(a, &error);

	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);

		//sigma = s * c
		//tau = t * c
		duspoly_t* sigma = multiplyAndModDUSP_DUZP(error, sp_s, Pptr);
		duspoly_t* tau = multiplyAndModDUSP_DUZP(error, sp_t, Pptr);

		//sigma = rem(sigma, w); s = quo(sgima, w)
		duspoly_t* q;
		plainDivPolynomialsInForm_spX_inp(&sigma, w, &q, Pptr);
		//in place or fast div???
		// duspoly_t* new_sigma;
		// fastDivPolynomialInForm_wPSInv_spX (sigma, w, &new_sigma, &q, Pptr);

		//tau = t*c + quo(sigma,w)*u
		duspoly_t* quProd;
		// mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
		freePolynomial_spX (&quProd);

		symmetricPadicUpdate(uu, tau, mod, Pptr, halfP, workZ);
		symmetricPadicUpdate(ww, sigma, mod, Pptr, halfP, workZ);

		//update the error.
		multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;

	}

	//make lifted polys primitive
	mpz_t cont;
	mpz_init(cont);
	primitivePartAndContent_DUZP_inp(uu, cont);
	*liftedU = uu;
	
//1
	// freePolynomial_DUZP(a);
//2
	//nothing

	freePolynomial_DUZP(error);
	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);
	freePolynomial_DUZP(ww);
	mpz_clears(bound, alpha, mod, cont, workZ, NULL);

	return primePow;
}

//Returns 0 iff step was unsuccessful. Either u and w are not co-prime or the error was zero.
//This method uses a bunch of temporary memory and out-of-place arithmetic. The following
//function is much better.
int quadraticPadicHenselStep_DUZP(const DUZP_t* a, DUZP_t* u, DUZP_t* w, DUZP_t* s, DUZP_t* t, mpz_t mod, DUZP_t** uu, DUZP_t** ww, DUZP_t** ss, DUZP_t** tt) {

	DUZP_t* error = multiplyPolynomials_DUZP(u,w);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	if (isZero_DUZP(error)) {
		freePolynomial_DUZP(error);
		return 0;
	}

	mpz_mul(mod, mod, mod);

	applyModuloSymmetric_DUZP_inp(error, mod); //This mod is useful before se = s*error if |a|_inf >> mod.

	DUZP_t* se = multiplyPolynomials_DUZP(s, error);
	applyModuloSymmetric_DUZP_inp(se, mod);
	DUZP_t* quo, *rem;
	dividePolynomials_DUZP(se, w, &quo, &rem);
	DUZP_t* wNext = addPolynomials_DUZP(w, rem);
	applyModuloSymmetric_DUZP_inp(wNext, mod);
	freePolynomial_DUZP(rem);
	freePolynomial_DUZP(se);

	DUZP_t* te = multiplyPolynomials_DUZP(t, error);
	DUZP_t* qu = multiplyPolynomials_DUZP(quo, u);
	freePolynomial_DUZP(quo);
	DUZP_t* tmp = addPolynomials_DUZP(u, te);
	DUZP_t* uNext = addPolynomials_DUZP(qu, tmp);
	applyModuloSymmetric_DUZP_inp(uNext, mod); 
	freePolynomial_DUZP(te);
	freePolynomial_DUZP(qu);
	freePolynomial_DUZP(tmp);

	//lift bezout coefficients now
	DUZP_t* su = multiplyPolynomials_DUZP(s, uNext);
	DUZP_t* tw = multiplyPolynomials_DUZP(t, wNext);
	DUZP_t* b = addPolynomials_DUZP(su,tw);
	mpz_sub_ui(b->coefs[0], b->coefs[0], 1ul);
	freePolynomial_DUZP(su);
	freePolynomial_DUZP(tw);

	DUZP_t* sb = multiplyPolynomials_DUZP(s, b);
	applyModuloSymmetric_DUZP_inp(sb, mod);
	// applyModulo_DUZP_inp(sb, mod);

	dividePolynomials_DUZP(sb, wNext, &quo, &rem);
	DUZP_t* sNext = subtractPolynomials_DUZP(s, rem);
	applyModuloSymmetric_DUZP_inp(sNext, mod);
	freePolynomial_DUZP(rem);

	DUZP_t* tb = multiplyPolynomials_DUZP(t, b);
	qu = multiplyPolynomials_DUZP(quo, uNext);
	DUZP_t* sum = addPolynomials_DUZP(tb, qu);
	DUZP_t* tNext = subtractPolynomials_DUZP(t, sum);
	applyModuloSymmetric_DUZP_inp(tNext, mod);

	freePolynomial_DUZP(b);
	freePolynomial_DUZP(qu);
	freePolynomial_DUZP(sum);
	freePolynomial_DUZP(quo);

	*uu = uNext;
	*ww = wNext;
	*ss = sNext;
	*tt = tNext;

	return 1;
}

//Returns 0 iff step was unsuccessful. Either u and w are not co-prime or the error was zero.
int quadraticPadicHenselStepOpt_DUZP(const DUZP_t* a, DUZP_t* u, DUZP_t* w, DUZP_t* s, DUZP_t* t, mpz_t mod, DUZP_t** error, DUZP_t** se, DUZP_t** te) {

	multiplyPolynomialsPreAlloc_DUZP(u, w, error);
	subtractPolynomials_DUZP_inpRHS(a, error);

	if (isZero_DUZP(*error)) {
		return 0;
	}

	mpz_mul(mod, mod, mod);

	applyModuloSymmetric_DUZP_inp(*error, mod); //This mod is useful before se = s*error if |a|_inf >> mod.

	multiplyPolynomialsPreAlloc_DUZP(s, *error, se);
	applyModuloSymmetric_DUZP_inp(*se, mod);
	DUZP_t* quo, *rem;
	dividePolynomials_DUZP(*se, w, &quo, &rem);
	
	addPolynomials_DUZP_inp(&w, rem);
	applyModuloSymmetric_DUZP_inp(w, mod);
	freePolynomial_DUZP(rem);

	multiplyPolynomialsPreAlloc_DUZP(quo, u, te);
	addPolynomials_DUZP_inp(&u, *te);
	freePolynomial_DUZP(quo);
	multiplyPolynomialsPreAlloc_DUZP(t, *error, te);
	addPolynomials_DUZP_inp(&u, *te);
	applyModuloSymmetric_DUZP_inp(u, mod); 

	//lift bezout coefficients now
	multiplyPolynomialsPreAlloc_DUZP(t, w, te);
	multiplyPolynomialsPreAlloc_DUZP(s, u, se);
	addPolynomials_DUZP_inp(se,*te);
	mpz_sub_ui((*se)->coefs[0], (*se)->coefs[0], 1ul);
	//b is now se.

	//compute sb and tb using error and te as working space
	multiplyPolynomialsPreAlloc_DUZP(s, *se, error);
	multiplyPolynomialsPreAlloc_DUZP(t, *se, te);
	applyModuloSymmetric_DUZP_inp(*error, mod);

	dividePolynomials_DUZP(*error, w, &quo, &rem);
	subtractPolynomials_DUZP_inp(&s, rem);
	applyModuloSymmetric_DUZP_inp(s, mod);
	freePolynomial_DUZP(rem);

	multiplyPolynomialsPreAlloc_DUZP(quo, u, se);
	freePolynomial_DUZP(quo);
	addPolynomials_DUZP_inp(te, *se);
	subtractPolynomials_DUZP_inp(&t, *te);
	applyModuloSymmetric_DUZP_inp(t, mod);

	return 1;

}

int monicTwoTermQuadraticPadicLift(const DUZP_t* aa, duspoly_t* up, duspoly_t* wp, DUZP_t** liftedU, DUZP_t** liftedW, const Prime_ptr* Pptr) {
	if (liftedU == NULL && liftedW == NULL) {
		return 0;
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(aa, bound);
	mpz_mul_si(bound, bound, 2l);

	//get bezout coefficients
	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;
	extGCDInForm_spX(up, wp, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX (&sp_one);

	mpz_t mod; 
	mpz_init_set_ui(mod, Pptr->prime);
	int primePow = 1;

	DUZP_t* u = convertFromDUSPSymmetric_DUZP(up, Pptr);
	DUZP_t* w = convertFromDUSPSymmetric_DUZP(wp, Pptr);
	DUZP_t* s = convertFromDUSPSymmetric_DUZP(sp_s, Pptr);
	DUZP_t* t = convertFromDUSPSymmetric_DUZP(sp_t, Pptr);	
	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);

	DUZP_t* work1 = makePolynomial_DUZP(aa->lt+1);
	DUZP_t* work2 = makePolynomial_DUZP(aa->lt+1);
	DUZP_t* work3 = makePolynomial_DUZP(aa->lt+1);
	int successfulStep = 1;

	while( successfulStep && mpz_cmp(mod, bound) < 0 ) {
		successfulStep = quadraticPadicHenselStepOpt_DUZP(aa, u, w, s, t, mod, &work1, &work2, &work3);

		primePow = primePow*2;
	}

	freePolynomial_DUZP(s);
	freePolynomial_DUZP(t);
	freePolynomial_DUZP(work1);
	freePolynomial_DUZP(work2);
	freePolynomial_DUZP(work3);


	if (liftedU != NULL) {
		*liftedU = u;
	} else {
		freePolynomial_DUZP(u);
	}

	if (liftedW != NULL) {
		*liftedW = w;
	} else {
		freePolynomial_DUZP(w);
	}

	mpz_clears(bound, mod, NULL);

	return primePow;

}


//TODO this doesn't work at all
// void twoTermQuadraticPadicLift(const DUZP_t* aa, duspoly_t* up, duspoly_t* wp, DUZP_t** liftedU, DUZP_t** liftedW, mpz_t gamma, const Prime_ptr* Pptr) {
// 	if (liftedU == NULL && liftedW == NULL) {
// 		return;
// 	}

// 	mpz_t bound;
// 	mpz_init(bound);
// 	infinityNorm_DUZP(aa, bound);
// 	mpz_mul_si(bound, bound, 2l);
// 	mpz_abs(gamma, gamma);
// 	mpz_mul(bound, bound, gamma);

// 	mpz_t alpha;
// 	mpz_init(alpha);
// 	mpz_set(alpha, aa->coefs[aa->lt]);

// 	DUZP_t* a = multiplyByInteger_DUZP(aa, gamma);
	
// 	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);	
// 	if (smallprimefield_convert_out(up->elems[up->lt], Pptr) != 1) {
// 		alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(up->elems[up->lt], Pptr), Pptr);
// 	}
// 	scalarMulPolynomialInForm_spX_inp (&up, alpha_sp, Pptr);

// 	if (smallprimefield_convert_out(up->elems[up->lt], Pptr) != 1) {
// 		scalarMulPolynomialInForm_spX_inp (&wp, smallprimefield_inv(wp->elems[wp->lt], Pptr), Pptr);
// 	}

// 	//get bezout coefficients
// 	duspoly_t* sp_s, *sp_t;
// 	duspoly_t* sp_one;
// 	extGCDInForm_spX(up, wp, &sp_s, &sp_t, &sp_one, Pptr);

// 	if (!isOneInForm_spX(sp_one, Pptr)) {
// 		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
// 		printPolynomialOutForm_spX(sp_one, Pptr);
// 		exit(1);
// 	}
// 	freePolynomial_spX (&sp_one);


// 	DUZP_t* u = convertFromDUSPSymmetric_DUZP(up, Pptr);
// 	// mpz_set(u->coefs[u->lt], gamma);
// 	DUZP_t* w = convertFromDUSPSymmetric_DUZP(wp, Pptr);
// 	// mpz_set(w->coefs[w->lt], alpha);
// 	DUZP_t* s = convertFromDUSPSymmetric_DUZP(sp_s, Pptr);
// 	DUZP_t* t = convertFromDUSPSymmetric_DUZP(sp_t, Pptr);	
// 	freePolynomial_spX (&sp_s);
// 	freePolynomial_spX (&sp_t);

// 	mpz_t mod; 
// 	mpz_init(mod);
// 	mpz_set_ui(mod, Pptr->prime);

// 	DUZP_t *uu, *ww, *ss, *tt;
// 	int successfulStep = 1;

// 	while( successfulStep && mpz_cmp(mod, bound) < 0 ) {

// 		successfulStep = quadraticPadicHenselStep_DUZP(aa, u, w, s, t, mod, &uu, &ww, &ss, &tt);

// 		freePolynomial_DUZP(u);
// 		freePolynomial_DUZP(w);
// 		freePolynomial_DUZP(s);
// 		freePolynomial_DUZP(t);

// 		u = uu;
// 		w = ww;
// 		s = ss;
// 		t = tt;
// 		uu = ww = ss = tt = NULL;
// 	}


// 	// DUZP_t* error = multiplyPolynomials_DUZP(uu, ww);
// 	// subtractPolynomials_DUZP_inpRHS(a, &error);

// 	// while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
// 	// 	//c = error / p^k
// 	// 	divideByIntegerExact_DUZP_inp(error, mod);

// 	// 	//sigma = s * c
// 	// 	//tau = t * c
// 	// 	duspoly_t* sigma = multiplyAndModDUSP_DUZP(error, sp_s, Pptr);
// 	// 	duspoly_t* tau = multiplyAndModDUSP_DUZP(error, sp_t, Pptr);

// 	// 	//sigma = rem(sigma, w); s = quo(sgima, w)
// 	// 	duspoly_t* q;
// 	// 	plainDivPolynomialsInForm_spX_inp(&sigma, w, &q, Pptr);
// 	// 	//in place or fast div???
// 	// 	// duspoly_t* new_sigma;
// 	// 	// fastDivPolynomialInForm_wPSInv_spX (sigma, w, &new_sigma, &q, Pptr);

// 	// 	//tau = t*c + quo(sigma,w)*u
// 	// 	duspoly_t* quProd;
// 	// 	mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
// 	// 	addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
// 	// 	freePolynomial_spX (&quProd);

// 	// 	symmetricPadicUpdate(uu, tau, mod, Pptr, halfP, workZ);
// 	// 	symmetricPadicUpdate(ww, sigma, mod, Pptr, halfP, workZ);

// 	// 	//update the error.
// 	// 	multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
// 	// 	subtractPolynomials_DUZP_inpRHS(a, &error);

// 	// 	mpz_mul_si(mod, mod, Pptr->prime);
// 	// }

// 	//make lifted polys primitive
// 	mpz_t cont;
// 	mpz_init(cont);
// 	primitivePartAndContent_DUZP_inp(u, cont);
	
// 	if (liftedU != NULL) {
// 		*liftedU = u;
// 	} else {
// 		freePolynomial_DUZP(u);
// 	}

// 	if (liftedW != NULL) {
// 		mpz_divexact(cont, gamma, cont);
// 		if (mpz_cmp_si(cont, 1l) != 0) {
// 			divideByIntegerExact_DUZP_inp(w, cont);
// 		}
// 		*liftedW = w;
// 	} else {
// 		freePolynomial_DUZP(w);
// 	}

// 	mpz_clears(bound, alpha, mod, cont, NULL);

// }

int gcdQuadraticPadicLift(const DUZP_t* a, duspoly_t* up, duspoly_t* wp, DUZP_t** liftedU, const Prime_ptr* Pptr) {
	if (liftedU == NULL) {
		return 0;
	}

	mpz_t alpha;
	mpz_init(alpha);
	mpz_set(alpha, a->coefs[a->lt]);

	//make w monic and give u the lc(a) mod p.
	if (smallprimefield_convert_out(wp->elems[wp->lt], Pptr) != 1) {
		scalarMulPolynomialInForm_spX_inp (&up, wp->elems[wp->lt], Pptr);
		scalarMulPolynomialInForm_spX_inp (&wp, smallprimefield_inv(wp->elems[wp->lt], Pptr), Pptr);
	}

	//get bezout coefficients
	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;
	extGCDInForm_spX(up, wp, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX (&sp_one);


	DUZP_t* u = convertFromDUSPSymmetric_DUZP(up, Pptr);
	// mpz_set(u->coefs[u->lt], gamma);
	DUZP_t* w = convertFromDUSPSymmetric_DUZP(wp, Pptr);
	// mpz_set(w->coefs[w->lt], alpha);
	DUZP_t* s = convertFromDUSPSymmetric_DUZP(sp_s, Pptr);
	DUZP_t* t = convertFromDUSPSymmetric_DUZP(sp_t, Pptr);	
	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	int primePow = 1;

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);
	mpz_abs(alpha, alpha);
	mpz_mul(bound, bound, alpha);


	// DUZP_t *uu, *ww, *ss, *tt;
	DUZP_t* work1 = makePolynomial_DUZP(a->lt+1);
	DUZP_t* work2 = makePolynomial_DUZP(a->lt+1);
	DUZP_t* work3 = makePolynomial_DUZP(a->lt+1);
	int successfulStep = 1;

	while( successfulStep && mpz_cmp(mod, bound) < 0 ) {

		// successfulStep = quadraticPadicHenselStep_DUZP(a, u, w, s, t, mod, &uu, &ww, &ss, &tt);
		successfulStep = quadraticPadicHenselStepOpt_DUZP(a, u, w, s, t, mod, &work1, &work2, &work3);

		primePow *= 2;

		// freePolynomial_DUZP(u);
		// freePolynomial_DUZP(w);
		// freePolynomial_DUZP(s);
		// freePolynomial_DUZP(t);

		// u = uu;
		// w = ww;
		// s = ss;
		// t = tt;
		// uu = ww = ss = tt = NULL;
	}


	primitivePart_DUZP_inp(u);	
	*liftedU = u;

	freePolynomial_DUZP(w);
	freePolynomial_DUZP(work1);
	freePolynomial_DUZP(work2);
	freePolynomial_DUZP(work3);

	mpz_clears(bound, alpha, mod, NULL);

	return primePow;

}

/**
 * u and w are modified throughout this algorithm to save space.
 */ 
int positiveTwoTermPadicLift(const DUZP_t* aa, duspoly_t* u, duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, const mpz_t gamma, const Prime_ptr* Pptr) {
	if (liftedU == NULL && liftedW == NULL) {
		return 0;
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(aa, bound);
	mpz_mul_si(bound, bound, 2l);
	mpz_mul(bound, bound, gamma);

	mpz_t alpha;
	mpz_init(alpha);
	mpz_set(alpha, aa->coefs[aa->lt]);

	//scale polynomials to fix l.c. problem
	DUZP_t* a = multiplyByInteger_DUZP(aa, gamma);
	
	prime_t gamma_sp = smallprimefield_convert_in(mpz_fdiv_ui(gamma, Pptr->prime), Pptr);
	if (smallprimefield_convert_out(u->elems[u->lt], Pptr) != 1) {
		gamma_sp = smallprimefield_mul(gamma_sp, smallprimefield_inv(u->elems[u->lt], Pptr), Pptr);
	}
	scalarMulPolynomialInForm_spX_inp (&u, gamma_sp, Pptr);

	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);	
	if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
		alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	}
	scalarMulPolynomialInForm_spX_inp (&w, alpha_sp, Pptr);

	//get bezout coefficients
	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;
	extGCDInForm_spX(u, w, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX (&sp_one);

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	int primePow = 1;

	DUZP_t* uu = convertFromDUSP_DUZP(u, Pptr);
	DUZP_t* ww = convertFromDUSP_DUZP(w, Pptr);
	mpz_set(uu->coefs[uu->lt], gamma);
	mpz_set(ww->coefs[ww->lt], alpha);

	DUZP_t* error = multiplyPolynomials_DUZP(uu, ww);
	subtractPolynomials_DUZP_inpRHS(a, &error);
	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);

		//sigma = s * c
		//tau = t * c
		duspoly_t* sigma = multiplyAndModDUSP_DUZP(error, sp_s, Pptr);
		duspoly_t* tau = multiplyAndModDUSP_DUZP(error, sp_t, Pptr);

		//sigma = rem(sigma, w); s = quo(sgima, w)
		duspoly_t* q;
		plainDivPolynomialsInForm_spX_inp(&sigma, w, &q, Pptr);
		//in place or fast div???
		// duspoly_t* new_sigma;
		// fastDivPolynomialInForm_wPSInv_spX (sigma, w, &new_sigma, &q, Pptr);

		//tau = t*c + quo(sigma,w)*u
		duspoly_t* quProd;
		// mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		mulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
		freePolynomial_spX (&quProd);

		padicUpdate(uu, tau, mod, Pptr);
		padicUpdate(ww, sigma, mod, Pptr);

		//update the error.
		multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}

	//make lifted polys primitive
	mpz_t cont;
	mpz_init(cont);
	primitivePartAndContent_DUZP_inp(uu, cont);
	
	freePolynomial_DUZP(error);
	freePolynomial_DUZP(a);
	freePolynomial_spX (&sp_s);
	freePolynomial_spX (&sp_t);

	if (liftedU != NULL) {
		*liftedU = uu;
	} else {
		freePolynomial_DUZP(uu);
	}

	if (liftedW != NULL) {
		mpz_divexact(cont, gamma, cont);
		if (mpz_cmp_si(cont, 1l) != 0) {
			divideByIntegerExact_DUZP_inp(ww, cont);
		}
		*liftedW = ww;
	} else {
		freePolynomial_DUZP(ww);
	}

	mpz_clears(bound, alpha, mod, cont, NULL);

	return primePow;
}

/**
 * Given a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime,
 * Where f[i] are monic and pair-wise co-prime,
 * Lift f[i] such that a = lc(a)lf[0]*...lf[nf-1] over Z. 
 */ 
int monicMultiTermPadicLift_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return 0;
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);

	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*)*nf);	

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);
	int primePow = 1;

	for (int i = 0; i < nf; ++i) {
		liftedF[i] = convertFromDUSPSymmetric_DUZP(f[i], Pptr);
	}

	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	duspoly_t* c;
	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);

		c = convertToDUSP_DUZP(error, Pptr);
		int pass = multiUDP_spX(c, f, nf, sigmas, Pptr);
		if (!pass) {
			break;
		}

		for (int i = 0; i < nf; ++i) {
			symmetricPadicUpdate(liftedF[i], sigmas[i], mod, Pptr, halfP, workZ);
			freePolynomial_spX (&sigmas[i]);
		}

		//update the error.
		freePolynomial_DUZP(error);
		error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
		
		subtractPolynomials_DUZP_inpRHS(a, &error);
		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;

	}
	
	freePolynomial_DUZP(error);
	free(sigmas);

	mpz_clears(mod, bound, workZ, NULL);

	return primePow;
}

/**
 * Given a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime,
 * Where f[i] are monic and pair-wise co-prime,
 * Lift f[i] such that a = lc(a)lf[0]*...lf[nf-1] over Z. 
 */ 
int monicMultiTermPadicLiftOpt_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return 0;
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);

	//get bezout coefficients
	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*)*nf);
	duspoly_t* sp_one = makePolynomial_spX(1);
	sp_one->elems[0] = smallprimefield_convert_in(1, Pptr);
	sp_one->lt = 1;
	int pass = multiUDP_spX(sp_one, f, nf, sigmas, Pptr);
	freePolynomial_spX (&sp_one);
	if (!pass) {
		fprintf(stderr, "Input polynomials to multiTermPadicLift were not co-prime!\n");
		exit(1);
	}

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);
	int primePow = 1;

	for (int i = 0; i < nf; ++i) {
		liftedF[i] = convertFromDUSPSymmetric_DUZP(f[i], Pptr);
	}

	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	duspoly_t* c;
	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {

		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);
		c = convertToDUSP_DUZP(error, Pptr);

		duspoly_t* prod;
		for (int i = 0; i < nf; ++i) {
			// mulPolynomialsInForm_spX (c, sigmas[i], &prod, Pptr);
			mulPolynomialsInForm_spX (c, sigmas[i], &prod, Pptr);
			plainRemPolynomialsInForm_spX_inp(&prod, f[i], Pptr);
			symmetricPadicUpdate(liftedF[i], prod, mod, Pptr, halfP, workZ);
			freePolynomial_spX (&prod);
			prod = NULL;
		}

		//update the error.
		multiplyManyPolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, &error);
		subtractPolynomials_DUZP_inpRHS(a, &error);
		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;

	}

	//make lifted polys primitive
	for (int i = 0; i < nf; ++i) {
		freePolynomial_spX (&sigmas[i]);
	}
	freePolynomial_DUZP(error);
	free(sigmas);

	mpz_clears(mod, bound, workZ, NULL);

	return primePow;
}

int multiTermPadicLiftResume(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, duspoly_t const*const* sigmas, short nf, int e, const mpz_t bound, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return 0;
	}

	mpz_t mod; 
	mpz_init_set_ui(mod, Pptr->prime);
	mpz_pow_ui(mod, mod, e);
	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);
	int primePow = e;

	mpz_set(liftedF[0]->coefs[liftedF[0]->lt], a->coefs[a->lt]);

	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	if (isZero_DUZP(error)) {
		// primitivePart_DUZP_inp(liftedF[0]);
		freePolynomial_DUZP(error);
		mpz_clears(mod, workZ, NULL);
		return e;
	}

	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(a->coefs[a->lt], Pptr->prime), Pptr);	
	duspoly_t* tmp_sp;
	scalarMulPolynomialInForm_spX(f[0], alpha_sp, &tmp_sp, Pptr);

	duspoly_t* c;
	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);
		c = convertToDUSP_DUZP(error, Pptr);

		duspoly_t* prod;
		// mulPolynomialsInForm_spX (c, sigmas[0], &prod, Pptr);
		mulPolynomialsInForm_spX (c, sigmas[0], &prod, Pptr);
		plainRemPolynomialsInForm_spX_inp(&prod, tmp_sp, Pptr);
		symmetricPadicUpdate(liftedF[0], prod, mod, Pptr, halfP, workZ);
		freePolynomial_spX (&prod);
		prod = NULL;
		for (int i = 1; i < nf; ++i) {
			// mulPolynomialsInForm_spX (c, sigmas[i], &prod, Pptr);
			mulPolynomialsInForm_spX(c, sigmas[i], &prod, Pptr);
			plainRemPolynomialsInForm_spX_inp(&prod, f[i], Pptr);
			symmetricPadicUpdate(liftedF[i], prod, mod, Pptr, halfP, workZ);
			freePolynomial_spX (&prod);
			prod = NULL;
		}

		//update the error.
		multiplyManyPolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, &error);
		subtractPolynomials_DUZP_inpRHS(a, &error);
		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}

	//make lifted non-monic poly primitive
	// primitivePart_DUZP_inp(liftedF[0]);

	if (!isZero_DUZP(error)) {
		//we got here because bound exceeded.
		//mod out l.c. for return;
		mpz_mod(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
		mpz_t halfMod;
		mpz_init_set(halfMod, mod);
		mpz_div_ui(halfMod, halfMod, 2ul);
		if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod) > 0) {
			mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
		}
		mpz_clear(halfMod);
	}

	freePolynomial_DUZP(error);
	mpz_clears(mod, workZ, NULL);

	return primePow;
}

/**
 * Given a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime,
 * Where f[i] are monic and pair-wise co-prime,
 * Lift f[i] such that a = lc(a)lf[0]*...lf[nf-1] over Z. 
 */ 
int multiTermPadicLiftStart(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, duspoly_t*** retSigmas, short nf, const mpz_t bound, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return 0;
	}

	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(a->coefs[a->lt], Pptr->prime), Pptr);	
	duspoly_t* tmp_sp;
	scalarMulPolynomialInForm_spX(f[0], alpha_sp, &tmp_sp, Pptr);

	liftedF[0] = convertFromDUSPSymmetric_DUZP(tmp_sp, Pptr);
	mpz_set(liftedF[0]->coefs[liftedF[0]->lt], a->coefs[a->lt]);
	for (int i = 1; i < nf; ++i) {
		liftedF[i] = convertFromDUSPSymmetric_DUZP(f[i], Pptr);
	}


	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	if (isZero_DUZP(error)) {
		freePolynomial_DUZP(error);
		freePolynomial_spX (&tmp_sp);
		return 1;
	}


	//get bezout coefficients
	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*)*nf);
	duspoly_t* sp_one = makePolynomial_spX(1);
	sp_one->elems[0] = smallprimefield_convert_in(1, Pptr);
	sp_one->lt = 1;

	const duspoly_t* udpInput[nf];
	udpInput[0] = tmp_sp;
	memcpy(udpInput + 1, f + 1, sizeof(duspoly_t*)*(nf-1));

	int pass = multiUDP_spX(sp_one, udpInput, nf, sigmas, Pptr);
	freePolynomial_spX (&sp_one);
	if (!pass) {
		fprintf(stderr, "Input polynomials to multiTermPadicLift were not co-prime!\n");
		exit(1);
	}

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);
	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);
	int primePow = 1;

	duspoly_t* c;
	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_DUZP_inp(error, mod);
		c = convertToDUSP_DUZP(error, Pptr);

		duspoly_t* prod = NULL;
		mulPolynomialsInForm_spX (c, sigmas[0], &prod, Pptr);

		// TEST
		fprintf(stderr, "[before] prod = \n");
		printPolynomialOutForm_spX (prod, Pptr);

		plainRemPolynomialsInForm_spX_inp(&prod, tmp_sp, Pptr);

		// TEST
		fprintf(stderr, "[after] prod = \n");
		printPolynomialOutForm_spX (prod, Pptr);

		symmetricPadicUpdate(liftedF[0], prod, mod, Pptr, halfP, workZ);
		freePolynomial_spX (&prod);
		prod = NULL;
		for (int i = 1; i < nf; ++i) {
			mulPolynomialsInForm_spX (c, sigmas[i], &prod, Pptr);

			// TEST
			fprintf(stderr, "[2-before] prod = \n");
			printPolynomialOutForm_spX (prod, Pptr);

			// TEST
			fprintf(stderr, "[2-before] f[i] = \n");
			printPolynomialOutForm_spX (f[i], Pptr);

			plainRemPolynomialsInForm_spX_inp(&prod, f[i], Pptr);

			// TEST
			fprintf(stderr, "[2-after] prod = \n");
			printPolynomialOutForm_spX (prod, Pptr);

			symmetricPadicUpdate(liftedF[i], prod, mod, Pptr, halfP, workZ);
			freePolynomial_spX (&prod);
			prod = NULL;
		}

		//update the error.
		multiplyManyPolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, &error);
		subtractPolynomials_DUZP_inpRHS(a, &error);
		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;

		freePolynomial_spX(&c);
	}

	
	//make lifted non-monic poly primitive
	// mpz_t cont;
	// mpz_init(cont);
	// primitivePart_DUZP_inp(liftedF[0]);
	// primitivePartAndContent_DUZP_inp(liftedF[0], cont);
	// gmp_fprintf(stderr, "Content: %Zd\n", cont);
	// mpz_clear(cont);

	if (!isZero_DUZP(error)) {
		//we got here because bound exceeded.
		//mod out l.c. for return;
		mpz_mod(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
		mpz_t halfMod;
		mpz_init_set(halfMod, mod);
		mpz_div_ui(halfMod, halfMod, 2ul);
		if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod) > 0) {
			mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
		}
		mpz_clear(halfMod);
	}

	if (retSigmas == NULL) {
		for (int i = 0; i < nf; ++i) {
			freePolynomial_spX (&sigmas[i]);
		}
		free(sigmas);	
	} else {
		*retSigmas = sigmas;
	}

	freePolynomial_DUZP(error);
	mpz_clears(mod, workZ, NULL);

	//reset input poly to be the monic one
	freePolynomial_spX (&tmp_sp);

	return primePow;
}

/**
 * Given a monic polynomial, a, a = f[0]*f[1]*...*f[nf-1] mod Pptr->prime, where f[i] are monic and pair-wise co-prime,
 * lift f[i] such that a = lc(a)liftedF[0]*...liftedF[nf-1] over Z.
 * This process is recursive in that it lifts by multiplying several factors together 
 * and solving many two-factor lifts.
 * It has problems with these products of factors being not co-prime for a two-factor lift.
 *
 * @param a a monic polynomial over the integers whose factors are to be lifted. 
 * @param f an array of nf polynomials ove ra finite field which are factors of a mouldo the prime.
 * @param[out] liftedF a pre-allocated array of size nf to store the lifted factors.
 * @param nf the number of factors in f to be lifted.
 * @param Pptr the prime of the finite field.
 * 
 */ 
void monicMultiTermPadicLiftOpt_Recursive(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return;
	}

	if (nf == 2) {
		monicTwoTermPadicLift(a, f[0], f[1], &(liftedF[0]), &(liftedF[1]), Pptr);
		return;
	}
	if (nf == 3) {
		monicMultiTermPadicLiftOpt_Iterative(a, f, liftedF, nf, Pptr);
		return;
	}

	//partition f into two sets and make intermediate products. 
	duspoly_t* f1, *f2;
	int f1Size = nf/2;
	int f2Size = nf - (nf/2);
	f1 = plainMulManyPolynomialsInForm_spX( CONSTCONSTCAST(duspoly_t, f), f1Size, Pptr);
	f2 = plainMulManyPolynomialsInForm_spX( CONSTCONSTCAST(duspoly_t, f + f1Size), f2Size, Pptr);

	DUZP_t* liftedF1, *liftedF2;
	monicTwoTermPadicLift(a, f2, f1, &(liftedF2), &(liftedF1), Pptr);

	//recursive calls;
	monicMultiTermPadicLiftOpt_Recursive(liftedF1, f, liftedF, f1Size, Pptr);
	monicMultiTermPadicLiftOpt_Recursive(liftedF2, f + f1Size, liftedF + f1Size, f2Size, Pptr);
}



/**
 * This sometimes segfaults...
 * It's also very slow compared ot the iterative approach...
 */
void multiTermPadicLiftOpt_Recursive(const DUZP_t* a, duspoly_t** f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return;
	}

	if (nf == 2) {
		if (mpz_cmp_ui(a->coefs[a->lt], 1ul)) {
			monicTwoTermPadicLift(a, f[0], f[1], &(liftedF[0]), &(liftedF[1]), Pptr);
			return;
		}
		mpz_t gamma;
		mpz_init(gamma);
		mpz_set_ui(gamma, 1ul);
		//in twoTermPadicLift, W always gets lt(a) as l.c. So swap 0 and 1 in order of params.
		twoTermPadicLift(a, f[1], f[0], &(liftedF[1]), &(liftedF[0]), gamma, Pptr);
		mpz_clear(gamma);
		return;
	}
	if (nf == 3) {
		multiTermPadicLiftOpt_Iterative(a, CONSTCONSTCAST(duspoly_t, f), liftedF, nf, Pptr);
		return;
	}

	//partition f into two sets and make intermediate products. 
	duspoly_t* f1, *f2;
	int f1Size = nf/2;
	int f2Size = nf - (nf/2);
	f1 = plainMulManyPolynomialsInForm_spX( CONSTCONSTCAST(duspoly_t, f), f1Size, Pptr);
	f2 = plainMulManyPolynomialsInForm_spX( CONSTCONSTCAST(duspoly_t, f + f1Size), f2Size, Pptr);

	mpz_t gamma;
	mpz_init(gamma);
	mpz_set_ui(gamma, 1ul);
	DUZP_t* liftedF1, *liftedF2;
	twoTermPadicLift(a, f2, f1, &(liftedF2), &(liftedF1), gamma, Pptr);
	mpz_clear(gamma);

	//recursive calls;
	multiTermPadicLiftOpt_Recursive(liftedF1, f, liftedF, f1Size, Pptr);
	multiTermPadicLiftOpt_Recursive(liftedF2, f + f1Size, liftedF + f1Size, f2Size, Pptr);

}

/**
 * Given a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime,
 * Where f[i] are monic and pair-wise co-prime,
 * Lift f[i] such that a = lc(a)lf[0]*...lf[nf-1] over Z.
 * nf >= 2. 
 */ 
int monicMultiTermQuadraticPadicLiftOpt_Iterative(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, short nf, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return 0;
	}

	if (nf < 2) {
		//not supported
		return 0;
	}

	for (int i = 0; i < nf; ++i) {
		liftedF[i] = convertFromDUSPSymmetric_DUZP(f[i], Pptr);
	}

	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	if (isZero_DUZP(error)) {
		//no lifting at all to do
		freePolynomial_DUZP(error);
		return 1;	
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(a, bound);
	mpz_mul_si(bound, bound, 2l);

	//get bezout coefficients
	duspoly_t** sigmas_sp = (duspoly_t**) malloc(sizeof(duspoly_t*)*nf);
	duspoly_t* sp_one = makePolynomial_spX(1);
	sp_one->elems[0] = smallprimefield_convert_in(1, Pptr);
	sp_one->lt = 1;
	int pass = multiUDP_spX(sp_one, f, nf, sigmas_sp, Pptr);
	freePolynomial_spX (&sp_one);
	if (!pass) {
		fprintf(stderr, "Input polynomials to multiTermPadicLift were not co-prime!\n");
		exit(1);
	}

	DUZP_t** sigmas = (DUZP_t**) sigmas_sp; //all pointers are same size so reuse allocation
	DUZP_t* tmp;
	for(int i = 0; i < nf; ++i) {
		//keep sigmas in positive prime field range under padic update
		tmp = convertFromDUSP_DUZP(sigmas_sp[i], Pptr);
		freePolynomial_spX (&sigmas_sp[i]);
		sigmas[i] = tmp;
	}

	mpz_t mod, mod2, halfMod; 
	mpz_init_set_ui(halfMod, Pptr->prime);
	mpz_div_ui(halfMod, halfMod, 2ul);
	mpz_init_set_ui(mod, Pptr->prime);
	mpz_init2(mod2, 128 +  mp_bits_per_limb); //prime <= 64 bits
	mpz_mul(mod2, mod, mod);
	int primePow = 1;


	DUZP_t* sigmaError = makePolynomial_DUZP(2*a->lt);
	DUZP_t* prod = makePolynomial_DUZP(2 * a->lt); //prealloc for work in loop
	DUZP_t** prods = (DUZP_t**) malloc(sizeof(DUZP_t*)*nf);
	for (int i = 0; i < nf; ++i) {
		prods[i] = makePolynomial_DUZP(2 * a->lt);
	}

	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {
		
		//c = error / p^k
		//keep error in positive prime field range under padic update
		applyModulo_DUZP_inp(error, mod2);
		divideByIntegerExact_DUZP_inp(error, mod);

		for (int i = 0; i < nf; ++i) {
			multiplyPolynomialsPreAlloc_DUZP(error, sigmas[i], &prod);
			remainder_DUZP_inp(&prod, liftedF[i]);
			applyModulo_DUZP_inp(prod, mod);
			symmetricQuadraticPadicUpdate(liftedF[i], prod, mod, halfMod);
		}

		//update the error.
		multiplyManyPolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, &error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		if (isZero_DUZP(error) || mpz_cmp(mod2, bound) > 0) {
			//then next loop will break so don't bother lifting sigmas
			break;
		}

		//update prods for sigma lift
		for (int i = 0; i < nf; ++i) {
			multiplyAllButOnePolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, i, &prods[i]);
		}

		//lift the sigmas now
		multiplyPolynomials_DUZP_inp(sigmas[0], &prods[0]);
		multiplyPolynomials_DUZP_inp(sigmas[1], &prods[1]);
		addPolynomialsPreAlloc_DUZP(prods[0], prods[1], &sigmaError);
		for(int i = 2; i < nf; ++i) {
			multiplyPolynomials_DUZP_inp(sigmas[i], &prods[i]);
			addPolynomials_DUZP_inp(&sigmaError, prods[i]);
		}
		negatePolynomial_DUZP_inp(sigmaError);
		mpz_add_ui(sigmaError->coefs[0], sigmaError->coefs[0], 1ul);

		applyModulo_DUZP_inp(sigmaError, mod2);
		divideByIntegerExact_DUZP_inp(sigmaError, mod);

		for (int i = 0; i < nf; ++i) {
			multiplyPolynomialsPreAlloc_DUZP(sigmaError, sigmas[i], &prod);
			remainder_DUZP_inp(&prod, liftedF[i]);
			applyModulo_DUZP_inp(prod, mod);
			quadraticPadicUpdate(sigmas[i], prod, mod); //keep sigmas in positive representaiton of Z/modZ
		}

		//update mods
		mpz_set(mod, mod2);
		mpz_div_ui(halfMod, mod, 2ul);
		mpz_mul(mod2, mod2, mod2);

		primePow *= 2;
	}

	for (int i = 0; i < nf; ++i) {
		freePolynomial_DUZP(sigmas[i]);
		freePolynomial_DUZP(prods[i]);
	}
	free(prods);
	free(sigmas);
	freePolynomial_DUZP(error);
	freePolynomial_DUZP(prod);
	freePolynomial_DUZP(sigmaError);

	mpz_clears(mod, bound, NULL);

	return primePow;
}

static inline void multiTermQuadraticPadicHenselStepSigmas_DUZP(DUZP_t* sigmaError, DUZP_t** liftedF, DUZP_t** sigmas, DUZP_t** prods, short nf, DUZP_t** prod, mpz_t mod, mpz_t mod2) {

	//update prods for sigma lift
	for (int i = 0; i < nf; ++i) {
		multiplyAllButOnePolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, i, &prods[i]);
	}

	//lift the sigmas now
	multiplyPolynomials_DUZP_inp(sigmas[0], &prods[0]);
	multiplyPolynomials_DUZP_inp(sigmas[1], &prods[1]);
	addPolynomialsPreAlloc_DUZP(prods[0], prods[1], &sigmaError);
	for(int i = 2; i < nf; ++i) {
		multiplyPolynomials_DUZP_inp(sigmas[i], &prods[i]);
		addPolynomials_DUZP_inp(&sigmaError, prods[i]);
	}
	negatePolynomial_DUZP_inp(sigmaError);
	mpz_add_ui(sigmaError->coefs[0], sigmaError->coefs[0], 1ul);

	applyModulo_DUZP_inp(sigmaError, mod2);
	divideByIntegerExact_DUZP_inp(sigmaError, mod);

	multiplyPolynomialsPreAlloc_DUZP(sigmaError, sigmas[0], prod);
	remainderMod_DUZP_inp(prod, liftedF[0], mod);
	quadraticPadicUpdate(sigmas[0], *prod, mod); //keep sigmas in positive representaiton of Z/modZ
	for (int i = 1; i < nf; ++i) {
		multiplyPolynomialsPreAlloc_DUZP(sigmaError, sigmas[i], prod);
		remainder_DUZP_inp(prod, liftedF[i]);
		applyModulo_DUZP_inp(*prod, mod);
		quadraticPadicUpdate(sigmas[i], *prod, mod); //keep sigmas in positive representaiton of Z/modZ
	}
}


static inline void multiTermQuadraticPadicHenselStepFactors_DUZP(const DUZP_t* a, DUZP_t** error, DUZP_t** liftedF, DUZP_t** sigmas, short nf, DUZP_t** prod, mpz_t mod, mpz_t mod2, mpz_t mod4, mpz_t halfMod, mpz_t halfMod2) {

	//c = error / p^k
	//keep error in positive prime field range under padic update
	applyModulo_DUZP_inp(*error, mod2);
	divideByIntegerExact_DUZP_inp(*error, mod);

	multiplyPolynomialsPreAlloc_DUZP(*error, sigmas[0], prod);
	remainderMod_DUZP_inp(prod, liftedF[0], mod);

	symmetricQuadraticPadicUpdate(liftedF[0], *prod, mod, halfMod);
	for (int i = 1; i < nf; ++i) {
		multiplyPolynomialsPreAlloc_DUZP(*error, sigmas[i], prod);
		remainder_DUZP_inp(prod, liftedF[i]);
		applyModulo_DUZP_inp(*prod, mod);
		symmetricQuadraticPadicUpdate(liftedF[i], *prod, mod, halfMod);
	}

	//update l.c. of f[0] for calculation of next error
	mpz_mul(mod4, mod2, mod2);
	mpz_set(halfMod, halfMod2);
	mpz_div_ui(halfMod2, mod4, 2ul);
	mpz_mod(liftedF[0]->coefs[liftedF[0]->lt], a->coefs[a->lt], mod4);
	if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod2) > 0) {
		mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod4);
	}

	//update the error.
	multiplyManyPolynomialsPreAlloc_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf, error);
	subtractPolynomials_DUZP_inpRHS(a, error);
}

/**
 * Given a = lc(a)*f[0]*f[1]*...*f[nf-1] mod Pptr->prime,
 * Where f[i] are monic and pair-wise co-prime,
 * Lift f[i] such that a = lc(a)lf[0]*...lf[nf-1] over Z.
 * nf >= 2. 
 */ 
int multiTermQuadraticPadicLiftStart(const DUZP_t* a, duspoly_t const*const* f, DUZP_t** liftedF, DUZP_t*** retSigmas, short nf, const mpz_t bound, const Prime_ptr* Pptr) {
	if (liftedF == NULL) {
		return 0;
	}

	if (nf < 2) {
		//not supported
		return 0;
	}

	mpz_t mod, mod2, halfMod, halfMod2; 
	mpz_init_set_ui(halfMod, Pptr->prime);
	mpz_div_ui(halfMod, halfMod, 2ul);
	mpz_init_set_ui(mod, Pptr->prime);
	mpz_init2(mod2, 128 +  mp_bits_per_limb); //prime <= 64 bits
	mpz_mul(mod2, mod, mod);
	mpz_init_set(halfMod2, mod2);
	mpz_div_ui(halfMod2, halfMod2, 2ul);
	int primePow = 1;

	//temporarily set f[0] to be multiplied through by lc(a).
	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(a->coefs[a->lt], Pptr->prime), Pptr);	
	duspoly_t* tmp2;
	scalarMulPolynomialInForm_spX(f[0], alpha_sp, &tmp2, Pptr);

	liftedF[0] = convertFromDUSPSymmetric_DUZP(tmp2, Pptr);
	for (int i = 1; i < nf; ++i) {
		liftedF[i] = convertFromDUSPSymmetric_DUZP(f[i], Pptr);
	}

	//set lc(lf[0]) to be lc(a) mod p^2.
	mpz_mod(liftedF[0]->coefs[liftedF[0]->lt], a->coefs[a->lt], mod2);
	if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod2) > 0) {
		mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod2);
	}

	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	if (isZero_DUZP(error)) {
		//no lifting at all to do
		freePolynomial_DUZP(error);

		freePolynomial_spX (&tmp2);

		mpz_clears(mod, mod2, halfMod, halfMod2, NULL);
		return 1;	
	}

	const duspoly_t* udpInput[nf];
	udpInput[0] = tmp2;
	memcpy(udpInput + 1, f + 1, sizeof(duspoly_t*)*(nf-1));

	//get bezout coefficients
	duspoly_t** sigmas_sp = (duspoly_t**) malloc(sizeof(duspoly_t*)*nf);
	duspoly_t* sp_one = makePolynomial_spX(1);
	sp_one->elems[0] = smallprimefield_convert_in(1, Pptr);
	sp_one->lt = 1;
	int pass = multiUDP_spX(sp_one, udpInput, nf, sigmas_sp, Pptr);
	freePolynomial_spX (&sp_one);
	if (!pass) {
		fprintf(stderr, "Input polynomials to multiTermPadicLift were not co-prime!\n");
		exit(1);
	}

	DUZP_t** sigmas = (DUZP_t**) sigmas_sp; //all pointers are same size so reuse allocation
	DUZP_t* tmp;
	for(int i = 0; i < nf; ++i) {
		//keep sigmas in positive prime field range under padic update
		tmp = convertFromDUSP_DUZP(sigmas_sp[i], Pptr);
		freePolynomial_spX (&sigmas_sp[i]);
		sigmas[i] = tmp;
	}
	freePolynomial_spX (&tmp2);

	//prealloc for work in loop
	DUZP_t* sigmaError = makePolynomial_DUZP(2*a->lt);
	DUZP_t* prod = makePolynomial_DUZP(2 * a->lt);
	DUZP_t** prods = (DUZP_t**) malloc(sizeof(DUZP_t*)*nf);
	for (int i = 0; i < nf; ++i) {
		prods[i] = makePolynomial_DUZP(2 * a->lt);
	}

	mpz_t mod4;
	mpz_init(mod4);
	
	while( !isZero_DUZP(error) && mpz_cmp(mod, bound) < 0 ) {

		multiTermQuadraticPadicHenselStepFactors_DUZP(a, &error, liftedF, sigmas, nf, &prod, mod, mod2, mod4, halfMod, halfMod2);

		if (retSigmas == NULL && (isZero_DUZP(error) || mpz_cmp(mod2, bound) > 0)) {
			//then next loop will break so don't bother lifting sigmas
			primePow *= 2;
			break;
		}

		multiTermQuadraticPadicHenselStepSigmas_DUZP(sigmaError, liftedF, sigmas, prods, nf, &prod, mod, mod2);

		//update mods
		//halfMod and halfMod2 updated above at mod4
		mpz_set(mod, mod2);
		mpz_set(mod2, mod4);

		primePow *= 2;
	}

	//clear out lc before return.
	// primitivePart_DUZP_inp(liftedF[0]);
	// mpz_t cont;
	// mpz_init(cont);
	// primitivePartAndContent_DUZP_inp(liftedF[0], cont);
	// gmp_fprintf(stderr, "quadratic content: %Zd\n", cont);
	// mpz_clear(cont);

	// if (!isZero_DUZP(error)) {
	// 	mpz_mod(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
	// 	if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod) > 0) {
	// 		mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
	// 	}
	// }

	if (retSigmas == NULL) {
		for (int i = 0; i < nf; ++i) {
			freePolynomial_DUZP(sigmas[i]);
			freePolynomial_DUZP(prods[i]);
		}
		free(sigmas);
	} else {
		for (int i = 0; i < nf; ++i) {
			freePolynomial_DUZP(prods[i]);
		}
		*retSigmas = sigmas;
	}

	free(prods);
	freePolynomial_DUZP(error);
	freePolynomial_DUZP(prod);
	freePolynomial_DUZP(sigmaError);

	mpz_clears(mod, mod2, mod4, halfMod, halfMod2, NULL);

	return primePow;
}

int multiTermQuadraticPadicLiftResume(const DUZP_t* a, DUZP_t** liftedF, DUZP_t ** sigmas, short nf, int e, const mpz_t targetBound, const Prime_ptr* Pptr) {

	if (liftedF == NULL || sigmas == NULL) {
		return -1;
	}

	if (nf < 2) {
		//not supported
		return -1;
	}


	mpz_t mod, mod2, halfMod, halfMod2; 
	mpz_init_set_ui(mod, Pptr->prime);
	mpz_pow_ui(mod, mod, e);

	mpz_init_set(halfMod, mod);	
	mpz_div_ui(halfMod, halfMod, 2ul);
	mpz_init_set(mod2, mod);
	mpz_mul(mod2, mod2, mod);

	mpz_init_set(halfMod2, mod2);
	mpz_div_ui(halfMod2, halfMod2, 2ul);
	int primePow = e;

	//set lc(lf[0]) to be lc(a) mod p^2.
	mpz_mod(liftedF[0]->coefs[liftedF[0]->lt], a->coefs[a->lt], mod2);
	if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod2) > 0) {
		mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod2);
	}

	DUZP_t* error = multiplyManyPolynomials_DUZP( CONSTCONSTCAST(DUZP_t, liftedF), nf);
	subtractPolynomials_DUZP_inpRHS(a, &error);

	if (isZero_DUZP(error)) {
		//no lifting at all to do
		freePolynomial_DUZP(error);

		mpz_clears(mod, mod2, halfMod, halfMod2, NULL);
		return primePow;	
	}


	//prealloc for work in loop
	DUZP_t* sigmaError = makePolynomial_DUZP(2*a->lt);
	DUZP_t* prod = makePolynomial_DUZP(2 * a->lt);
	DUZP_t** prods = (DUZP_t**) malloc(sizeof(DUZP_t*)*nf);
	for (int i = 0; i < nf; ++i) {
		prods[i] = makePolynomial_DUZP(2 * a->lt);
	}

	mpz_t mod4;
	mpz_init(mod4);

	while( !isZero_DUZP(error) && mpz_cmp(mod, targetBound) < 0 ) {

		multiTermQuadraticPadicHenselStepFactors_DUZP(a, &error, liftedF, sigmas, nf, &prod, mod, mod2, mod4, halfMod, halfMod2);
		
		multiTermQuadraticPadicHenselStepSigmas_DUZP(sigmaError, liftedF, sigmas, prods, nf, &prod, mod, mod2);

		//update mods
		//halfMod and halfMod2 updated above at mod4
		mpz_set(mod, mod2);
		mpz_set(mod2, mod4);

		primePow *= 2;
	}

	//clear out lc before return.
	// primitivePart_DUZP_inp(liftedF[0]);

	// if (!isZero_DUZP(error)) {
	// 	mpz_mod(liftedF[0]->elems[liftedF[0]->lt], liftedF[0]->elems[liftedF[0]->lt], mod);
	// 	if (mpz_cmp(liftedF[0]->coefs[liftedF[0]->lt], halfMod) > 0) {
	// 		mpz_sub(liftedF[0]->coefs[liftedF[0]->lt], liftedF[0]->coefs[liftedF[0]->lt], mod);
	// 	}
	// }

	for (int i = 0; i < nf; ++i) {
		freePolynomial_DUZP(prods[i]);
	}
	free(prods);

	freePolynomial_DUZP(error);
	freePolynomial_DUZP(prod);
	freePolynomial_DUZP(sigmaError);

	mpz_clears(mod, mod2, mod4, halfMod, halfMod2, NULL);

	return primePow;
}




























/***********************************
 *
 * Operations specialized to GCD
 *
 ***********************************/


DUZP_t* monicUnivarHenselGCD_DUZP(const DUZP_t* a, const DUZP_t* b, const Prime_ptr* Pptr) {

	duspoly_t* ap = convertToDUSP_DUZP(a, Pptr);
	duspoly_t* bp = convertToDUSP_DUZP(b, Pptr);
	duspoly_t* g = NULL;

	GCDInForm_spX(ap, bp, &g, Pptr);
	if (isOneInForm_spX(g, Pptr)) {
		DUZP_t* ret = convertFromDUSP_DUZP(g, Pptr);
		freePolynomial_spX (&ap);
		freePolynomial_spX (&bp);
		freePolynomial_spX (&g);
		return ret;
	}

	//get cofactor
	duspoly_t* cofact_a;
	plainExactDivPolynomialsInForm_spX(ap, g, &cofact_a, Pptr);
	duspoly_t* one_sp;
	GCDInForm_spX(g, cofact_a, &one_sp, Pptr);

	DUZP_t* liftedG;
	if (!isOneInForm_spX(one_sp, Pptr)) {
		freePolynomial_spX (&one_sp);
		freePolynomial_spX (&cofact_a);
		plainExactDivPolynomialsInForm_spX(bp, g, &cofact_a, Pptr);
		GCDInForm_spX(g, cofact_a, &one_sp, Pptr);
		if (!isOneInForm_spX(one_sp, Pptr)) {
			freePolynomial_spX (&one_sp);
			freePolynomial_spX (&cofact_a);
			freePolynomial_spX (&ap);
			freePolynomial_spX (&bp);
			return primitiveGCD_DUZP(a,b);
			// fprintf(stderr, "Both cofactors are not relatively prime to GCD!\n");
			// fprintf(stderr, "a: \n");
			// printPoly_DUZP(a, "x");
			// fprintf(stderr, "b: \n");
			// printPoly_DUZP(b, "x");
		}
		freePolynomial_spX (&one_sp);
		monicTwoTermQuadraticPadicLift(b, g, cofact_a, &liftedG, NULL, Pptr);
	} else {
	 	// monicTwoTermPadicLift(a, g, cofact_a, &liftedG, NULL, Pptr);
		monicTwoTermQuadraticPadicLift(a, g, cofact_a, &liftedG, NULL, Pptr);
	}

	freePolynomial_spX (&ap);
	freePolynomial_spX (&bp);
 	freePolynomial_spX (&g);
 	freePolynomial_spX (&cofact_a);
 	return liftedG;
}

DUZP_t* univarHenselGCD_DUZP(const DUZP_t* a, const DUZP_t* b, const Prime_ptr* Pptr) {

	//get contents and prim parts first then convert then get gcd.
	mpz_t contA, contB, contG, gamma;
	mpz_inits(contA, contB, contG, gamma, NULL);

	DUZP_t* primA = primitivePartAndContent_DUZP(a, contA);
	DUZP_t* primB = primitivePartAndContent_DUZP(b, contB);
	mpz_gcd(contG, contA, contB);
	mpz_divexact(contA, contA, contG);
	mpz_gcd(gamma, primA->coefs[primA->lt], primB->coefs[primB->lt]);
	//TODO use better gamma?

	duspoly_t* ap = convertToDUSP_DUZP(primA, Pptr);
	duspoly_t* bp = convertToDUSP_DUZP(primB, Pptr);
	duspoly_t* g = NULL;

	GCDInForm_spX(ap, bp, &g, Pptr);
	if (isOneInForm_spX(g, Pptr)) {
		DUZP_t* ret = convertFromDUSP_DUZP(g, Pptr);
		freePolynomial_spX (&ap);
		freePolynomial_spX (&bp);
		freePolynomial_spX (&g);
		freePolynomial_DUZP(primA);
		freePolynomial_DUZP(primB);
		mpz_clears(contA, contB, contG, gamma, NULL);
		return ret;
	}

	//get cofactor
	duspoly_t* cofact_a;
	plainExactDivPolynomialsInForm_spX(ap, g, &cofact_a, Pptr);

	DUZP_t* liftedG;
	duspoly_t* one_sp;
	GCDInForm_spX(g, cofact_a, &one_sp, Pptr);
	if (!isOneInForm_spX(one_sp, Pptr)) {
		freePolynomial_spX (&one_sp);
		freePolynomial_spX (&cofact_a);
		plainExactDivPolynomialsInForm_spX(bp, g, &cofact_a, Pptr);
		GCDInForm_spX(g, cofact_a, &one_sp, Pptr);
		if (!isOneInForm_spX(one_sp, Pptr)) {
			liftedG = primitiveGCD_DUZP(primA, primB);
			// fprintf(stderr, "a: \n");
			// printPoly_DUZP(a, "x");
			// fprintf(stderr, "b: \n");
			// printPoly_DUZP(b, "x");
		} else {
		 	// twoTermQuadraticPadicLift(b, g, cofact_a, &liftedG, NULL, gamma, Pptr);
		 	// twoTermPadicLift(primB, g, cofact_a, &liftedG, NULL, gamma, Pptr);
		 	gcdPadicLift(primB, g, cofact_a, &liftedG, Pptr);
		 	// gcdQuadraticPadicLift(primB, g, cofact_a, &liftedG, Pptr);
		}
	} else {
	 	// twoTermQuadraticPadicLift(a, g, cofact_a, &liftedG, NULL, gamma, Pptr);
	 	// twoTermPadicLift(primA, g, cofact_a, &liftedG, NULL, gamma, Pptr);
	 	gcdPadicLift(primA, g, cofact_a, &liftedG, Pptr);
	 	// gcdQuadraticPadicLift(primA, g, cofact_a, &liftedG, Pptr);
		
	}

	//now multiply by proper content
	// primitivePart_DUZP_inp(liftedG); //already primitive by padic lift
	multiplyByInteger_DUZP_inp(liftedG, contG);
	mpz_clears(contA, contB, contG, gamma, NULL);

	freePolynomial_spX (&one_sp);
 	freePolynomial_spX (&g);
 	freePolynomial_spX (&cofact_a);

 	freePolynomial_DUZP(primA);
 	freePolynomial_DUZP(primB);
 	return liftedG;
}

DUZP_t* univarQuadraticHenselGCD_DUZP(const DUZP_t* a, const DUZP_t* b, const Prime_ptr* Pptr) {

	//get contents and prim parts first then convert then get gcd.
	mpz_t contA, contB, contG, gamma;
	mpz_inits(contA, contB, contG, gamma, NULL);

	DUZP_t* primA = primitivePartAndContent_DUZP(a, contA);
	DUZP_t* primB = primitivePartAndContent_DUZP(b, contB);
	mpz_gcd(contG, contA, contB);
	mpz_divexact(contA, contA, contG);
	mpz_gcd(gamma, primA->coefs[primA->lt], primB->coefs[primB->lt]);
	//TODO use better gamma?

	duspoly_t* ap = convertToDUSP_DUZP(primA, Pptr);
	duspoly_t* bp = convertToDUSP_DUZP(primB, Pptr);
	duspoly_t* g = NULL;

	GCDInForm_spX(ap, bp, &g, Pptr);
	if (isOneInForm_spX(g, Pptr)) {
		DUZP_t* ret = convertFromDUSP_DUZP(g, Pptr);
		freePolynomial_spX (&ap);
		freePolynomial_spX (&bp);
		freePolynomial_spX (&g);
		return ret;
	}

	//get cofactor
	duspoly_t* cofact_a;
	plainExactDivPolynomialsInForm_spX(ap, g, &cofact_a, Pptr);

	DUZP_t* liftedG;
	duspoly_t* one_sp;
	GCDInForm_spX(g, cofact_a, &one_sp, Pptr);
	if (!isOneInForm_spX(one_sp, Pptr)) {
		freePolynomial_spX (&one_sp);
		freePolynomial_spX (&cofact_a);
		plainExactDivPolynomialsInForm_spX(bp, g, &cofact_a, Pptr);
		GCDInForm_spX(g, cofact_a, &one_sp, Pptr);
		if (!isOneInForm_spX(one_sp, Pptr)) {
			liftedG = primitiveGCD_DUZP(primA, primB);
			// fprintf(stderr, "a: \n");
			// printPoly_DUZP(a, "x");
			// fprintf(stderr, "b: \n");
			// printPoly_DUZP(b, "x");
		} else {
		 	// twoTermQuadraticPadicLift(a, g, cofact_a, &liftedG, NULL, gamma, Pptr);
		 	// twoTermPadicLift(primA, g, cofact_a, &liftedG, NULL, gamma, Pptr);
		 	// gcdPadicLift(primB, g, cofact_a, &liftedG, Pptr);
		 	gcdQuadraticPadicLift(primB, g, cofact_a, &liftedG, Pptr);
		}
	} else {
	 	// twoTermQuadraticPadicLift(a, g, cofact_a, &liftedG, NULL, gamma, Pptr);
	 	// twoTermPadicLift(primA, g, cofact_a, &liftedG, NULL, gamma, Pptr);
	 	gcdQuadraticPadicLift(primA, g, cofact_a, &liftedG, Pptr);
		
	}

	//now multiply by proper content
	// primitivePart_DUZP_inp(liftedG); //already primitive by padic lift
	multiplyByInteger_DUZP_inp(liftedG, contG);
	mpz_clears(contA, contB, contG, gamma, NULL);

	freePolynomial_spX (&one_sp);
 	freePolynomial_spX (&g);
 	freePolynomial_spX (&cofact_a);

 	freePolynomial_DUZP(primA);
 	freePolynomial_DUZP(primB);
 	return liftedG;
}

