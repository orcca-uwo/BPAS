void twoTermPadicLift(const DUZP_t* aa, duspoly_t* u, duspoly_t* w, DUZP_t** liftedU, DUZP_t** liftedW, mpz_t gamma, Prime_ptr* Pptr) {
	if (liftedU == NULL && liftedW == NULL) {
		return;
	}

	mpz_t bound;
	mpz_init(bound);
	infinityNorm_DUZP(aa, bound);
	mpz_mul_si(bound, bound, 2l);
	mpz_mul(bound, bound, gamma);

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
	// scalarMulPolynomialInForm_spX_inp(u, gamma_sp, Pptr);

	// prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);	
	// if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
	// 	alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	// }
	// scalarMulPolynomialInForm_spX_inp(w, alpha_sp, Pptr);
	
//2
	DUZP_t* a = (DUZP_t*) aa;
	prime_t alpha_sp = smallprimefield_convert_in(mpz_fdiv_ui(alpha, Pptr->prime), Pptr);
	if (smallprimefield_convert_out(u->elems[u->lt], Pptr) != 1) {
		alpha_sp = smallprimefield_mul(alpha_sp, smallprimefield_inv(u->elems[u->lt], Pptr), Pptr);
	}
	scalarMulPolynomialInForm_spX_inp(u, alpha_sp, Pptr);

	if (smallprimefield_convert_out(w->elems[w->lt], Pptr) != 1) {
		scalarMulPolynomialInForm_spX_inp(w, smallprimefield_inv(w->elems[w->lt], Pptr), Pptr);
	}
	
	fprintf(stderr, "u after lc adjust:\n");
	printPolynomialOutForm_spX(u, Pptr);
	fprintf(stderr, "w after lc adjust:\n");
	printPolynomialOutForm_spX(w, Pptr);

	//get bezout coefficients
	duspoly_t* sp_s, *sp_t;
	duspoly_t* sp_one;
	extGCDInForm_spX(u, w, &sp_s, &sp_t, &sp_one, Pptr);

	if (!isOneInForm_spX(sp_one, Pptr)) {
		fprintf(stderr, "Co-factors are not relatively prime in monicTwoTermPadicLift!\n");
		printPolynomialOutForm_spX(sp_one, Pptr);
		exit(1);
	}
	freePolynomial_spX(sp_one);

	mpz_t mod; 
	mpz_init(mod);
	mpz_set_ui(mod, Pptr->prime);

	unsigned long long int halfP = (Pptr->prime - 1) >> 1;
	mpz_t workZ;
	mpz_init2(workZ, mp_bits_per_limb + 64);

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
		fprintf(stderr, "\nlift loop\n" );
		gmp_fprintf(stderr, "a->lc: %Zd\nu->lt: %Zd\n", a->coefs[a->lt], uu->coefs[uu->lt]);
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
		plainMulPolynomialsInForm_spX (q, u, &quProd, Pptr);
		addPolynomialsInForm_spX_inp (&tau, quProd, Pptr);
		freePolynomial_spX(quProd);

		symmetricPadicUpdate(uu, tau, mod, Pptr, halfP, workZ);
		symmetricPadicUpdate(ww, sigma, mod, Pptr, halfP, workZ);

		//update the error.
		multiplyPolynomialsPreAlloc_DUZP(uu,ww,&error);
		subtractPolynomials_DUZP_inpRHS(a, &error);

		mpz_mul_si(mod, mod, Pptr->prime);
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
	freePolynomial_spX(sp_s);
	freePolynomial_spX(sp_t);

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

}
