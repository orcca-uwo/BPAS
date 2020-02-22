

#include "IntegerPolynomial/SMZP_Hensel.h"
#include "IntegerPolynomial/DUZP_Hensel.h"
#include "LinearAlgebra/Vandermonde.h"

#include "Utils/MacroHelpers.h"
#include "Utils/RandomHelpers.h"

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
const char* henselsyms[] = {"x", "y", "z", "s", "t", "u", "v", "w", "a", "b", "c", "d", "e", "f"};
#endif

duspoly_t* univarToPrimeField_AAZ(const AltArrZ_t* a, const Prime_ptr* Pptr) {
	if (isZero_AAZ(a)) {
		return zeroPolynomial_spX();
	}

	if (a->nvar > 1) {
		fprintf(stderr, "univarToPrimeField_AAZ ERROR: Input poly is not univariate!\n");
		exit(1);
	}

	if (a->nvar == 0) {
		duspoly_t* ppoly = zeroPolynomial_spX();
		ppoly->elems[0] = smallprimefield_convert_in(mpz_fdiv_ui(a->elems->coef, Pptr->prime), Pptr);
		return ppoly;
	}

	polysize_t maxExp = totalDegree_AAZ(a);
	duspoly_t* ppoly= makePolynomial_spX(maxExp+1);
	ppoly->lt = maxExp;
	elem_t* pelems = ppoly->elems;

	AAZElem_t* aelems = a->elems;	

	register unsigned long uP = (unsigned long) Pptr->prime;
	register polysize_t size = a->size;

	degree_t deg;
	for (polysize_t i = 0; i < size; ++i) {
		deg = partialDegreeTerm_AAZ(a, i, 0);
		pelems[deg] = smallprimefield_convert_in(mpz_fdiv_ui(aelems[i].coef, uP), Pptr);
	}

	return ppoly;
}


void univarToPrimeFieldPreAlloc_AAZ(const AltArrZ_t* a, const Prime_ptr* Pptr, duspoly_t** pp) {
	if (pp == NULL) {
		return;
	}

	if (isZero_AAZ(a)) {
		if ((*pp) == NULL) {
			*pp = zeroPolynomial_spX();
		} else {
			(*pp)->lt = 0;
			(*pp)->elems[0] = smallprimefield_convert_in(0ul, Pptr);
		}
		return;
	}

	if (a->nvar > 1) {
		fprintf(stderr, "univarToPrimeField_AAZ ERROR: Input poly is not univariate!\n");
		exit(1);
	}

	if (a->nvar == 0) {
		if (*pp == NULL) {
			*pp = zeroPolynomial_spX();
		}
		(*pp)->lt = 0;
		(*pp)->elems[0] = smallprimefield_convert_in(mpz_fdiv_ui(a->elems->coef, Pptr->prime), Pptr);
		return;
	}

	polysize_t maxExp = a->elems->degs;
	duspoly_t* ppoly = *pp;

	if (ppoly == NULL) {
		ppoly = makePolynomial_spX(maxExp+1);
	} else if (ppoly->alloc <= maxExp) {
		reallocatePolynomial_spX(&ppoly, maxExp+1);
	}
	ppoly->lt = maxExp;
	elem_t* pelems = ppoly->elems;

	AAZElem_t* aelems = a->elems;	

	register unsigned long uP = (unsigned long) Pptr->prime;
	register polysize_t size = a->size;

	for (polysize_t i = 0; i <= maxExp; ++i) {
		pelems[i] = smallprimefield_convert_in(0, Pptr);
	}

	degree_t deg;
	for (polysize_t i = 0; i < size; ++i) {
		deg = partialDegreeTerm_AAZ(a, i, 0);
		pelems[deg] = smallprimefield_convert_in(mpz_fdiv_ui(aelems[i].coef, uP), Pptr);
	}

	return;
}

AltArrZ_t* deepCopyPolynomial_AAZFromDUSP(const duspoly_t* ppoly, const Prime_ptr* Pptr) {
	if (isZero_spX(ppoly)) {
		return NULL;
	}

	polysize_t size = ppoly->lt + 1;
	elem_t* pelems = ppoly->elems;

	AltArrZ_t* aa = makePolynomial_AAZ(size, 1);
	AAZElem_t* aelems = aa->elems;
	int insertIdx = 0;
	for(polysize_t i = ppoly->lt; i >= 0; --i) {
		if (pelems[i] != 0) {
			aelems[insertIdx].degs = i;
			mpz_init_set_ui(aelems[insertIdx].coef, smallprimefield_convert_out(pelems[i], Pptr));
			++insertIdx;
		}
	}

	aa->size = insertIdx;
	// resizePolynomial_AAZ(aa, insertIdx); //isn't really needed
	return aa;
}

//convert a DUSP to a bivariate polynomial where the main variable is the same var as the DUSP.
AltArrZ_t* convertFromDUSPBivariate_AAZ(const duspoly_t* ppoly, const Prime_ptr* Pptr) {
	if (isZero_spX(ppoly)) {
		return NULL;
	}

	polysize_t size = ppoly->lt + 1;
	elem_t* pelems = ppoly->elems;

	AltArrZ_t* aa = makePolynomial_AAZ(size, 2);
	AAZElem_t* aelems = aa->elems;
	int insertIdx = 0;
	elem_t halfP = Pptr->prime >> 1;
	elem_t tmp;
	for(polysize_t i = ppoly->lt; i >= 0; --i) {
		if (pelems[i] != 0) {
			aelems[insertIdx].degs = (i << EXP_OFFSET_1_V2);
			tmp = smallprimefield_convert_out(pelems[i], Pptr);
			if (tmp > halfP) {
				tmp -= Pptr->prime;
			}
			mpz_init_set_si(aelems[insertIdx].coef, tmp);
			++insertIdx;
		}
	}

	aa->size = insertIdx;
	// resizePolynomial_AAZ(aa, insertIdx); //isn't really needed
	return aa;
}

//it is assumed that n < MAX_EXP_NVAR_VNVAR
//(i.e. the maximum exponent size of the nth variable in an exponent vector of size n)
AltArrZ_t* binomialExpansion(const mpz_t x, unsigned int n, int nvar) {

	if (mpz_sgn(x) == 0) {
		AltArrZ_t* a = makePolynomial_AAZ(1, nvar);
		mpz_init_set_ui(a->elems[0].coef, 1ul);
		a->elems[0].degs = n;
		a->size = 1;
		return a;
	}

	AltArrZ_t* a = makePolynomial_AAZ(n+1u, nvar);
	degrees_t i = n;

	mpz_t nfact;
	mpz_init_set_ui(nfact, 1ul);
	for (unsigned int r = 2; r <= n; ++r) {
		mpz_mul_ui(nfact, nfact, r);
	}

	mpz_t xPow;
	mpz_init_set(xPow, x);

	mpz_t nrfact;
	mpz_init_set(nrfact, nfact);
	mpz_t rfact;
	mpz_init_set_ui(rfact, 1ul);
	mpz_t coef;
	mpz_init(coef);

	AAZElem_t* elems = a->elems;
	mpz_init_set_ui(elems[0].coef, 1ul);
	elems[0].degs = i;
	--i;
	for (unsigned int r = 1; r <= n; ++r) {
		mpz_divexact_ui(nrfact, nrfact, n - r + 1u );

		mpz_mul_ui(rfact, rfact, r);

		elems[r].degs = i;
		--i;
		mpz_init_set(elems[r].coef, xPow);
		mpz_mul(coef, rfact, nrfact);
		mpz_divexact(coef, nfact, coef);
		mpz_mul(elems[r].coef, elems[r].coef, coef);
		if (r % 2 == 1) {
			mpz_neg(elems[r].coef, elems[r].coef);
		}

		mpz_mul(xPow, xPow, x);
	}

	mpz_clears(xPow, coef, nrfact, rfact, nfact, NULL);


	a->size = n+1;
	return a;
}

/**
 * a = a + p*(x2 - x)^k
 */
void idealadicUpdate(AltArrZ_t** a, const duspoly_t* p, const mpz_t x, int k, const Prime_ptr* Pptr, const mpz_t mpPrime) {
	if (a == NULL) {
		 return;
	}

	AltArrZ_t* b = binomialExpansion(x, k, 2);
	AltArrZ_t* pa = convertFromDUSPBivariate_AAZ(p, Pptr);
	AltArrZ_t* up = multiplyPolynomials_AAZ(b, pa, 2);

	*a = addPolynomials_AAZ_inp(*a, up, 2);
	//is the mod really neaded? YES! b will have large coefficeints as powers of evaluation point
	applyModuloSymmetric_AAZ_inp(*a, mpPrime);
	
	freePolynomial_AAZ(up);
	freePolynomial_AAZ(pa);
}

/** 
 * a = a + p*(x_v - x)^k
 */
void multivarIdealadicUpdate(AltArrZ_t** a, AltArrZ_t* p, const mpz_t x, int k, const mpz_t mpPrime) {
	if (a == NULL) {
		return;
	}

	int nvar = (*a)->nvar;
	AltArrZ_t* b = binomialExpansion(x, k, nvar);
	expandNumVars_AAZ(p, nvar); //expand right and full in gap with x_v from b
	AltArrZ_t* up = multiplyPolynomials_AAZ(p, b, nvar);
	shrinkNumVarsAtIdx_AAZ(p, nvar-1);//shrink p back
	*a = addPolynomials_AAZ_inp(*a, up, nvar);
	applyModuloSymmetric_AAZ_inp(*a, mpPrime);

	freePolynomial_AAZ(up);
}


//assumes variable to lift is index 1 in a. i.e. the second variable, i.e. not the main variable
//a, uu, ww, are viewed as being over the finite field Z / Pptr->prime*Z
int monicBivarTwoFactorHenselLift_AAZ(const AltArrZ_t* a, const AltArrZ_t* uu, const AltArrZ_t* ww, const mpz_t x, const Prime_ptr* Pptr, AltArrZ_t** liftedU, AltArrZ_t** liftedW) {

	if (liftedU == NULL && liftedW == NULL) {
		return 0;
	}

	if (uu->nvar != 1 || ww->nvar != 1) {
		fprintf(stderr, "Images were not univariate in monicBivarTwoFactorHenselLift_AAZ\n" );
		return 0;
	}

	mpz_t mpPrime;
	mpz_init_set_si(mpPrime, Pptr->prime);

	AltArrZ_t* u = deepCopyPolynomial_AAZ(uu);
	expandNumVars_AAZ(u, 2);
	AltArrZ_t* w = deepCopyPolynomial_AAZ(ww);
	expandNumVars_AAZ(w, 2);

	AltArrZ_t* error = multiplyPolynomials_AAZ(u,w,2);
	subPolynomials_AAZ_inpRHS(a, &error);	 
	applyModuloSymmetric_AAZ_inp(error, mpPrime);

	duspoly_t* up = univarToPrimeField_AAZ(uu, Pptr);
	duspoly_t* wp = univarToPrimeField_AAZ(ww, Pptr);
	duspoly_t* sigma, *tau;

	degree_t deg = partialDegree_AAZ(a, 1);
	AltArrZ_t* kDeriv;
	AltArrZ_t* c;
	duspoly_t* cp = makePolynomial_spX(a->size + 1);

	mpz_t vals[2];
	mpz_init(vals[0]);
	mpz_set_ui(vals[0], 1ul);
	vals[1][0] = x[0]; //temporarily put x into vals array
	int active[2] = {0, 1};//only evaluate second variable 
	
	for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {
		kDeriv = derivative_AAZ(error, 1, k);

		c = evaluatePoly_AAZ(kDeriv, active, vals, 2);
		freePolynomial_AAZ(kDeriv);
		divideByIntegerExact_AAZ_inp(c, vals[0]); //divide by k!		
		univarToPrimeFieldPreAlloc_AAZ(c, Pptr, &cp);
		freePolynomial_AAZ(c);

		UDP_spX(up, wp, cp, &sigma, &tau, Pptr);

		//ideal-adic update
		//u += sigma * (x2 - x);
		idealadicUpdate(&u, tau, x, k, Pptr, mpPrime);
		//w += tau * (x2 - x);
		idealadicUpdate(&w, sigma, x, k, Pptr, mpPrime);

		freePolynomial_spX (&sigma);
		freePolynomial_spX (&tau);

		multiplyPolynomialsPreAlloc_AAZ(u, w, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);

		mpz_mul_ui(vals[0], vals[0], k+1);
	}	

	mpz_clear(vals[0]);

	freePolynomial_spX (&up);
	freePolynomial_spX (&wp);

	int ret;
	if (!isZero_AAZ(error)) {
		ret = 0;
	} else {
		ret = 1;
	}

	freePolynomial_AAZ(error);

	if (liftedU != NULL) {
		*liftedU = u;
	} else {
		freePolynomial_AAZ(u);
	}

	if (liftedW != NULL) {
		*liftedW = w;
	} else {
		freePolynomial_AAZ(w);
	}

	return ret;
}


//assumes variable to lift is index 1 in a. i.e. the second variable, i.e. not the main variable
//a, uu, ww, are viewed as being over the finite field Z / Pptr->prime*Z
int monicBivarHenselLift_AAZ(const AltArrZ_t* a, AltArrZ_t const*const* fs, unsigned int r, const mpz_t x, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	if (a == NULL || fs == NULL || liftedF == NULL || r < 2) {
		return 0;
	}

	if (r == 2) {
		monicBivarTwoFactorHenselLift_AAZ(a, fs[0], fs[1], x, Pptr, liftedF, liftedF + 1);
		return 0;
	}

	for (unsigned int i = 0; i < r; ++i) {
		if (fs[i]->nvar != 1) {
			fprintf(stderr, "Images were not univariate in monicBivarTwoFactorHenselLift_AAZ\n" );
			return 0;
		}
	}

	mpz_t mpPrime;
	mpz_init_set_si(mpPrime, Pptr->prime);

	duspoly_t** fps = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	for (unsigned int i = 0; i < r; ++i) {
		fps[i] = univarToPrimeField_AAZ(fs[i], Pptr);
		liftedF[i] = deepCopyPolynomial_AAZ(fs[i]);
		expandNumVars_AAZ(liftedF[i], 2);
	}

	AltArrZ_t* error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r);
	subPolynomials_AAZ_inpRHS(a, &error);	 
	applyModuloSymmetric_AAZ_inp(error, mpPrime);

	degree_t deg = partialDegree_AAZ(a, 1);
	AltArrZ_t* kDeriv;
	AltArrZ_t* c;
	duspoly_t* cp = makePolynomial_spX(a->size + 1);

	mpz_t vals[2];
	mpz_init(vals[0]);
	mpz_set_ui(vals[0], 1ul);
	vals[1][0] = x[0]; //temporarily put x into vals array
	int active[2] = {0, 1};//only evaluate second variable 
	
	for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {
		kDeriv = derivative_AAZ(error, 1, k);

		c = evaluatePoly_AAZ(kDeriv, active, vals, 2);
		freePolynomial_AAZ(kDeriv);
		divideByIntegerExact_AAZ_inp(c, vals[0]); //divide by k!		
		univarToPrimeFieldPreAlloc_AAZ(c, Pptr, &cp);
		freePolynomial_AAZ(c);

		multiUDP_spX(cp, CONSTCONSTCAST(duspoly_t, fps), r, sigmas, Pptr);

		for (unsigned int i = 0; i < r; ++i) {
			idealadicUpdate(liftedF + i, sigmas[i], x, k, Pptr, mpPrime);
			freePolynomial_spX (&sigmas[i]);
		}


		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);

		mpz_mul_ui(vals[0], vals[0], k+1);
	}	

	mpz_clear(vals[0]);

	for (unsigned int i = 0; i < r; ++i) {
		freePolynomial_spX (&fps[i]);
	}
	free(fps);		
	free(sigmas);

	int ret = isZero_AAZ(error);
	freePolynomial_AAZ(error);
	return ret;
}

AltArrZ_t* bivarReplaceLC_DUZP(const DUZP_t* p, const DUZP_t* lc, const Prime_ptr* Pptr, const mpz_t mpPrime) {

	if (p == NULL || lc == NULL || Pptr == NULL) {
		return NULL;
	}

	if (lc->lt == 0 && mpz_cmp_ui(lc->coefs[0], 1l) == 0) {
		AltArrZ_t* ret = convertToAltArrZ_DUZP(p);
		expandNumVars_AAZ(ret, 2);
		return ret;
	}

	AltArrZ_t* a = makePolynomial_AAZ(p->lt + lc->lt + 2, 2); //nvar = 2;

	//assuming a is packed;
	degree_t deg = p->lt;
	AAZElem_t* elems = a->elems;
	int aIdx = 0;
	for (int i = lc->lt; i >= 0; --i) {
		if (mpz_sgn(lc->coefs[i]) == 0) {
			continue;
		}
		mpz_init_set(elems[aIdx].coef, lc->coefs[i]);
		// mpz_fdiv_ui(elems[aIdx].coef, Pptr->prime);
		elems[aIdx].degs = ((degrees_t) deg << EXP_OFFSET_1_V2) | ((degrees_t) i);
		++aIdx;
	}

	for (int i = p->lt-1; i >= 0; --i) {
		if (mpz_sgn(p->coefs[i]) == 0) {
			continue;
		}
		mpz_init_set(elems[aIdx].coef, p->coefs[i]);

		// mpz_fdiv_ui(elems[aIdx].coef, Pptr->prime);
		elems[aIdx].degs = ((degrees_t) i << EXP_OFFSET_1_V2);
		++aIdx;
	}

	a->size = aIdx;

	return a;
}

//forward declare
int bivarHenselLiftCoefs_AAZ(const AltArrZ_t* a, AltArrZ_t** Fs, unsigned int r, const mpz_t bound, const Prime_ptr* Pptr);

//assumes variable to lift is index 1 in a. i.e. the second variable, i.e. not the main variable
//a, uu, ww, are viewed as being over the finite field Z / Pptr->prime*Z
int bivarHenselLiftVars_AAZ(const AltArrZ_t* aa, DUZP_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t x, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	if (aa == NULL || lcF == NULL || fs == NULL || liftedF == NULL || r < 2) {
		return 0;
	}

	mpz_t mpPrime;
	mpz_init_set_si(mpPrime, Pptr->prime);

	AltArrZ_t* a = deepCopyPolynomial_AAZ(aa);
	applyModuloSymmetric_AAZ_inp(a, mpPrime);


	duspoly_t** fps = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	for (unsigned int i = 0; i < r; ++i) {
		fps[i] = convertToDUSP_DUZP(fs[i], Pptr);
		liftedF[i] = bivarReplaceLC_DUZP(fs[i], lcF[i], Pptr, mpPrime);
		applyModuloSymmetric_AAZ_inp(liftedF[i], mpPrime);
	}

	AltArrZ_t* error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r);
	subPolynomials_AAZ_inpRHS(a, &error);	 
	applyModuloSymmetric_AAZ_inp(error, mpPrime);

	degree_t deg = partialDegree_AAZ(a, 1);
	AltArrZ_t* kDeriv;
	AltArrZ_t* c;
	duspoly_t* cp = makePolynomial_spX(a->size + 1);

	mpz_t vals[2];
	mpz_init(vals[0]);
	mpz_set_ui(vals[0], 1ul);
	vals[1][0] = x[0]; //temporarily put x into vals array
	int active[2] = {0, 1};//only evaluate second variable 

	for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {
		kDeriv = derivative_AAZ(error, 1, k);

		c = evaluatePoly_AAZ(kDeriv, active, vals, 2);
		freePolynomial_AAZ(kDeriv);
		divideByIntegerExact_AAZ_inp(c, vals[0]); //divide by k!		
		univarToPrimeFieldPreAlloc_AAZ(c, Pptr, &cp);
		freePolynomial_AAZ(c);

		multiUDP_spX(cp, CONSTCONSTCAST(duspoly_t, fps), r, sigmas, Pptr);

		for (unsigned int i = 0; i < r; ++i) {
			idealadicUpdate(liftedF + i, sigmas[i], x, k, Pptr, mpPrime);
			freePolynomial_spX (&sigmas[i]);
		}


		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);

		mpz_mul_ui(vals[0], vals[0], k+1);
	}	

	mpz_clear(vals[0]);

	for (unsigned int i = 0; i < r; ++i) {
		freePolynomial_spX (&fps[i]);
	}
	free(fps);		
	free(sigmas);

	int ret = isZero_AAZ(error);
	freePolynomial_AAZ(error);
	freePolynomial_AAZ(a);
	return ret;
}

int bivarHenselLift_AAZ(const AltArrZ_t* a, DUZP_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t x, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	int ret = bivarHenselLiftVars_AAZ(a, CONSTCONSTCAST(DUZP_t,lcF), fs, r, x, Pptr, liftedF);

	if (ret) {
		mpz_t bound;
		mpz_init(bound);
		infinityNorm_AAZ(a, bound);
		mpz_mul(bound, bound, a->elems->coef);
		degree_t pdegs[a->nvar];
		double degProd = 1;
		partialDegrees_AAZ(a, pdegs);
		for (int i = 0; i < a->nvar; ++i) {
			degProd *= (pdegs[i] + 1.0);
		}
		degProd = sqrt(degProd);
		mpz_mul_si(bound, bound, (long) (degProd + 0.5));
		degree_t deg = totalDegree_AAZ(a);
		mpz_mul_2exp(bound, bound, deg);

		//if bound is already smaller than prime, we must have captured true coefs
		//with just variable lifting over Z_p.
		if (mpz_cmp_ui(bound, Pptr->prime) < 0) {
			mpz_clear(bound);
			return ret;
		}

		AltArrZ_t* lc;
		for (int i = 0; i < r; ++i) {
			// fprintf(stderr, "before replace: \n" );
			// printPoly_AAZ(stderr, liftedF[i], henselsyms, liftedF[i]->nvar);
			// fprintf(stderr, "\n" );
			// printPoly_DUZP(lcF[i], "x");
			lc = convertToAltArrZ_DUZP(lcF[i]);

			replaceLC_AAZ_inp(liftedF + i, lc);
			freePolynomial_AAZ(lc);

			// fprintf(stderr, "\nafter replace: \n" );
			// printPoly_AAZ(stderr, liftedF[i], henselsyms, liftedF[i]->nvar);
			// fprintf(stderr, "\n" );
		}

		// fprintf(stderr, "a: \n" );
		// printPoly_AAZ(stderr, a, henselsyms, a->nvar);
		// fprintf(stderr, "\n\n\n" );


		ret = bivarHenselLiftCoefs_AAZ(a, liftedF, r, bound, Pptr);

		mpz_t tmp;
		mpz_init_set_ui(tmp, Pptr->prime);
		mpz_pow_ui(tmp, tmp, ret);
		ret = mpz_cmp(tmp, bound) < 0;

		mpz_clear(bound);
		mpz_clear(tmp);
	}


	return ret;
}


/**
 * Multiply the polynomial pointed to by a_ptr by the binomial (x-b), in place.
 */
void univarMultByBinomial_AAZ_inp(AltArrZ_t** a_ptr, const mpz_t b) {
	if (a_ptr == NULL) {
		return;
	}

	AltArrZ_t* a = *a_ptr;
	if (a == NULL) {
		a = makePolynomial_AAZ(2, 1);
		a->elems[0].degs = 1ull;
		mpz_init_set_ui(a->elems[0].coef, 1ul);
		a->elems[1].degs = 0ull;
		mpz_init_set(a->elems[1].coef, b);
		a->size = 2;
		*a_ptr = a;
		return;
	}

	if (a->alloc < a->size + 1) {
		resizePolynomial_AAZ(a, a->size*2); //ammortize resize cost hopefully;
	}

	AAZElem_t* elems = a->elems;
	elems[a->size].degs = 0ull;
	mpz_init(elems[a->size].coef);
	mpz_mul(elems[a->size].coef, elems[a->size-1].coef, b);
	mpz_neg(elems[a->size].coef, elems[a->size].coef);
	for (int i = a->size-1; i > 0; --i) {
		mpz_submul(elems[i].coef, elems[i-1].coef, b);
		elems[i].degs += 1ull;
	}
	elems[0].degs += 1ull;
	++(a->size);
}


/**
 * Incrementally interpolate the jth point (x[j],c) using the polynomial pointed to
 * by a_ptr as the existing interpolant. Returns a 1 on successful interpolation.
 * Failure may occur if the points being interpolated cannot be fitted by an integer
 * polynomial. 
 *
 *
 * @param a_ptr a pointer to the existing interpolating polynomial to be updated
 * @param x the list of interpolation nodes
 * @param j the index of the new point to interpolate within x. 
 * @param bj the value corresponding with x[j] forming the point to interpolate.  
 *
 * @return 1 if the interpolation was successful, 
 *         0 if the interpolation did not add degree to the interpolant,
 *        -1 if the interpolation was not successful.
 */
int univarIncrNewtonInterp_AAZ(AltArrZ_t** a_ptr, const mpz_t* x, int j, const mpz_t bj) {

	if (a_ptr == NULL) {
		return 0;
	}

	AltArrZ_t* a = *a_ptr;
	if (j == 0) {
		if (a == NULL) {
			a = makeConstPolynomial_AAZ(1, 1, bj);
			*a_ptr = a;
			return 1;
		}
		if (a->size > 1) {
			for (int i = 1; i < a->size; ++i) {
				mpz_clear(a->elems[i].coef);
			}
			mpz_set(a->elems[0].coef, bj);
		} else if (a->size == 0) {
			mpz_init_set(a->elems[0].coef, bj);
		}
		return 1;
	}



	//phi_j = \prod_{k=1..j-1}(x - x[k])
	//p_3(x) = p_2(x) + (bj - p_2(xi))/()

	//build the newton polynomial for this step while simultaneously evaluating it
	mpz_t s,t;
	mpz_init(s);
	mpz_init_set(t, x[j]);
	mpz_sub(t, t, x[0]);
	AltArrZ_t* phi = makePolynomial_AAZ(j+1,1);
	mpz_init_set_ui(phi->elems[0].coef, 1ul);
	mpz_init_set(phi->elems[1].coef, x[0]);
	mpz_neg(phi->elems[1].coef, phi->elems[1].coef);
	phi->elems[0].degs = 1;
	phi->elems[1].degs = 0;
	phi->size = 2;
	for (int i = 1; i < j; ++i) {
		univarMultByBinomial_AAZ_inp(&phi, x[i]);	

		mpz_sub(s, x[j], x[i]);
		mpz_mul(t, t, s);	
	}

	if (mpz_sgn(t) == 0) {
		//interpolation nodes were not unique!
		mpz_clears(s,t,NULL);
		return 0;
	}

	mpz_t y;
	mpz_init(y);
	univarEvaluate_AAZ(a, x[j], y);
	// gmp_fprintf(stderr, "a evaluated at %Zd to %Zd\n", x[j], y);
	mpz_sub(y, bj, y);
	int ret;
	// gmp_fprintf(stderr, "bg - y = %Zd\n", y);
	if (mpz_sgn(y) == 0) {
		ret = 0;
	} else if (mpz_divisible_p(y, t)) {
		mpz_div(y, y, t);
		multiplyByInteger_AAZ_inp(phi, y);
		int degPrev = a->elems[0].degs;
		a = addPolynomials_AAZ_inp(a, phi, 1);
		*a_ptr = a;

		int degsAfter = a->elems[0].degs;
		if (degPrev == degsAfter) {
			ret = 0;
		} else {
			ret = 1;
		}
	} else {
		//The thing being interpolated is not an integer polynomial.
		fprintf(stderr, "univarIncrNewtonInterp_AAZ failed to interpolate an Integer polynomial.\n");
		// exit(1);
		ret = -1;
	}


	mpz_clears(s, t, y, NULL);

	return ret;
}

/**
 * Incrementally interpolate the jth point (x[j],c) using the polynomial pointed to
 * by a_ptr as the existing interpolant. Returns 1 on successful interpolation.
 * Failure may occur if the interpolation nodes are not unique.
 * The interpolation occurs modulo Pptr->prime.
 *
 *
 * @param a_ptr a pointer to the existing interpolating polynomial to be updated
 * @param x the list of interpolation nodes
 * @param j the index of the new point to interpolate within x. 
 * @param bj the value corresponding with x[j] forming the point to interpolate.  
 * @param Pptr the modulus 
 *
 * @return 1 if the interpolation was successful, 
 *         0 if the interpolation did not add degree to the interpolant,
 *        -1 if the interpolation was not successful.
 */
int univarIncrNewtonInterpModP_AAZ(AltArrZ_t** a_ptr, const mpz_t* x, int j, const mpz_t bj, const Prime_ptr* Pptr) {

	if (a_ptr == NULL) {
		return 0;
	}

	AltArrZ_t* a = *a_ptr;
	if (j == 0) {
		if (a == NULL) {
			a = makeConstPolynomial_AAZ(1, 1, bj);
			*a_ptr = a;
			return 1;
		}
		if (a->size > 1) {
			for (int i = 1; i < a->size; ++i) {
				mpz_clear(a->elems[i].coef);
			}
			mpz_set(a->elems[0].coef, bj);
		} else if (a->size == 0) {
			mpz_init_set(a->elems[0].coef, bj);
		}
		return 1;
	}


	mpz_t p;
	mpz_init_set_ui(p, Pptr->prime);

	//phi_j = \prod_{k=1..j-1}(x - x[k])
	//p_3(x) = p_2(x) + (bj - p_2(xi))/()


	//build the newton polynomial for this step while simultaneously evaluating it
	mpz_t s,t;
	mpz_init(s);
	mpz_init_set(t, x[j]);
	mpz_sub(t, t, x[0]);
	AltArrZ_t* phi = makePolynomial_AAZ(j+1,1);
	mpz_init_set_ui(phi->elems[0].coef, 1ul);
	mpz_init_set(phi->elems[1].coef, x[0]);
	mpz_neg(phi->elems[1].coef, phi->elems[1].coef);
	phi->elems[0].degs = 1;
	phi->elems[1].degs = 0;
	phi->size = 2;
	for (int i = 1; i < j; ++i) {
		univarMultByBinomial_AAZ_inp(&phi, x[i]);	

		mpz_sub(s, x[j], x[i]);
		mpz_mul(t, t, s);	
	}

	if (mpz_sgn(t) == 0) {
		//interpolation nodes were not unique!
		mpz_clears(s,t,NULL);
		return 0;
	}

	mpz_t y;
	mpz_init(y);
	univarEvaluate_AAZ(a, x[j], y);
	// gmp_fprintf(stderr, "a evaluated at %Zd to %Zd\n", x[j], y);
	mpz_sub(y, bj, y);
	int ret;
	// gmp_fprintf(stderr, "bg - y = %Zd\n", y);
	if (mpz_sgn(y) == 0) {
		ret = 0;
	} else if (mpz_divisible_p(y, t)) {
		mpz_div(y, y, t);
		multiplyByInteger_AAZ_inp(phi, y);
		int degPrev = a->elems[0].degs;
		a = addPolynomials_AAZ_inp(a, phi, 1);
		*a_ptr = a;

		int degsAfter = a->elems[0].degs;
		if (degPrev == degsAfter) {
			ret = 0;
		} else {
			ret = 1;
		}
	} else {
		mpz_invert(t, t, p);
		mpz_mul(y, y, t);
		multiplyByInteger_AAZ_inp(phi, y);

		int degPrev = a->elems[0].degs;
		a = addPolynomials_AAZ_inp(a, phi, 1);
		applyModuloSymmetric_AAZ_inp(a, p);
		*a_ptr = a;

		int degsAfter = a->elems[0].degs;
		if (degPrev == degsAfter) {
			ret = 0;
		} else {
			ret = 1;
		}
	}


	mpz_clears(s, t, y, p, NULL);

	return ret;
}

/**
 * Divide the univariate polynomial (a) by a monic lienar term (b). 
 * 
 * @param a the dividend
 * @param b the divisor
 * @param[out] q a pointer to the resulting quotient; may be NULL;
 * @param[out] rem a pointer to the resulting remainder; may be NULL. 
 *
 * @return 1 if the division was exact, 
 *         0 if the division was not exact, 
 *        -1 if an error occurred (such as invalid parameters; 
 */
//TODO NOTE This only works if a is dense. 
int divideUnivarMonicLinear_AAZ(AltArrZ_t* a, AltArrZ_t* b, AltArrZ_t** q, mpz_t* rem) {
	if (q == NULL && rem == NULL) {
		return -1;
	}

	if (a->nvar != b->nvar || b->nvar != 1) {
		return -1;
	}
	if (b->size != 2 || mpz_cmp_ui(b->elems->coef, 1ul) != 0 || b->elems->degs != 1ull) {
		return -1;
	}

	degree_t curDeg; 
	AAZElem_t* aelems = a->elems;
	int asize = a->size;
	mpz_t* div = &(b->elems[1].coef);
	//+1 here because we are using the last term of quo below to compute remainder.
	AltArrZ_t* quo = makePolynomial_AAZ(a->elems->degs+1, 1);
	mpz_init_set(quo->elems->coef, a->elems->coef);
	curDeg = quo->elems->degs = a->elems->degs - 1;
	int quoIdx = 1;
	int prevZero = 0;
	for (int i = 1; i < asize-1; ++i) {
		--curDeg;
		quo->elems[quoIdx].degs = curDeg;

		if (prevZero) {
			mpz_set(quo->elems[quoIdx].coef, aelems[i].coef);
			++quoIdx;
		} else {
			mpz_init_set(quo->elems[quoIdx].coef,  aelems[i].coef);
			mpz_submul(quo->elems[quoIdx].coef, quo->elems[quoIdx-1].coef, *div);
			if (mpz_sgn(quo->elems[quoIdx].coef) == 0) {
				prevZero = 1;
			} else {
				++quoIdx;
			}
		}
	}

	if (!prevZero) {
		mpz_init(quo->elems[quoIdx].coef);
		mpz_set(quo->elems[quoIdx].coef, aelems[a->size-1].coef);
		mpz_submul(quo->elems[quoIdx].coef, quo->elems[quoIdx-1].coef, *div);
	} else {
		mpz_set(quo->elems[quoIdx].coef, aelems[a->size].coef);
	}

	int ret = (mpz_sgn(quo->elems[quoIdx].coef) == 0);
	if (rem != NULL) {
		mpz_set(*rem, quo->elems[quoIdx].coef);
	} 
	mpz_clear(quo->elems[quoIdx].coef);
	
	quo->size = quoIdx;
	if (q == NULL) {
		freePolynomial_AAZ(quo);
	} else {
		*q = quo;
	}

	return ret;
}

/** 
 * Evaluate all variables except the main variable of the input polynomial,
 * returning the result as a univariate DUZP polynomial.
 * 
 * @param aa the polynomial to evaluate
 * @param vals a list of aa->nvar-1 values at which to evaluate aa.
 *
 * @return the resulting evaluated polynomial in a dense representation DUZP.  
 */
DUZP_t* evaluatePolyToDUZP_AAZ(const AltArrZ_t* aa, const mpz_t* vals) {
	if (vals == NULL || isZero_AAZ(aa)) {
		return NULL;
	}

	if (aa->nvar == 0) {
		return makeBigConstPolynomial_DUZP(1, aa->elems[0].coef);
	}
	if (aa->nvar == 1) {
		return convertFromAltArrZ_DUZP(aa);
	}

	short nvar = aa->nvar;

	// int* sizes = getExpOffsetArray(nvar);
	// unsigned long long int* masks = getExpMaskArray(nvar);

	mpz_t* valList[nvar-1];
	int valListSize[nvar-1];

	degree_t maxDegs[nvar];
	partialDegrees_AAZ(aa, maxDegs);
	for (int j = 0; j < nvar-1; ++j) {
		//j+1 because we evaluate all but main variable
		degree_t deg = maxDegs[j+1];
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(deg+1));
		valListSize[j] = deg+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	DUZP_t* res = makePolynomial_DUZP(maxDegs[0]+1);
	res->lt = maxDegs[0];
	int lt_r = res->lt;
	for (int i = 0; i <= lt_r; ++i) {
		mpz_init(res->coefs[i]);
	}

	// degree_t curDeg = maxDegs[0];
	mpz_t curVal;
	mpz_init(curVal);
	degree_t termDegs[nvar];
	for (int i = 0; i < aa->size; ++i) {
		mpz_set(curVal, aa->elems[i].coef);
		partialDegreesTerm_AAZ(aa, i, termDegs);
		for (int j = 0; j < nvar-1; ++j) {
			mpz_mul(curVal, curVal, valList[j][termDegs[j+1]]);
		}
		// curDeg = GET_NTH_EXP(aa->elems[i].degs, masks[0], sizes[0]);
		mpz_add(res->coefs[termDegs[0]], res->coefs[termDegs[0]], curVal);
	}

	for (int j = 0; j < nvar-1; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}
	mpz_clear(curVal);

	return res;
}

/** 
 * Evaluate all variables except the main variable of the input polynomial,
 * returning the result as a univariate DUZP polynomial.
 * 
 * @param aa the polynomial to evaluate
 * @param vals a list of aa->nvar-1 values at which to evaluate aa.
 * @param res_p a pointer to the (possibly pre-allocated) duspoly_t to be returned after evaluation.
 * @param Pptr the Prime_ptr describing the finite field.
 *
 * @return the resulting evaluated polynomial in a dense representation DUZP.  
 */
void evaluatePolyToDUSP_AAZ(const AltArrZ_t* aa, const mpz_t* vals, duspoly_t** res_p, const Prime_ptr* Pptr) {
	if (res_p == NULL) {
		return;
	}
	if (vals == NULL || isZero_AAZ(aa)) {
		*res_p = NULL;
		return;
	}

	if (aa->nvar == 0) {
		if (*res_p == NULL) {
			*res_p = makePolynomial_spX(1);
		}
		(*res_p)->lt = 0;
		(*res_p)->elems[0] = smallprimefield_convert_in(mpz_fdiv_ui(aa->elems[0].coef, Pptr->prime), Pptr);
		return;
	}

	if (aa->nvar == 1) {
		if (*res_p != NULL) {
			free(*res_p);
		}
		*res_p = univarToPrimeField_AAZ(aa, Pptr);
		return;
	}


	short nvar = aa->nvar;

	// int* sizes = getExpOffsetArray(nvar);
	// unsigned long long int* masks = getExpMaskArray(nvar);

	mpz_t* valList[nvar-1];
	int valListSize[nvar-1];

	degree_t maxDegs[nvar];
	partialDegrees_AAZ(aa, maxDegs);
	for (int j = 0; j < nvar-1; ++j) {
		//j+1 because we evaluate all but main variable
		degree_t deg = maxDegs[j+1];
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(deg+1));
		valListSize[j] = deg+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	duspoly_t* res = *res_p;
	if (res == NULL) {
		res = makePolynomial_spX(maxDegs[0]+1);
		*res_p = res;
	} else if (res->alloc < maxDegs[0] + 1) {
    	reallocatePolynomial_spX (res_p, maxDegs[0]+1);
    	res = *res_p;
   	}


   	// fprintf(stderr, "doing work in eval to DUSP\n" );
	res->lt = maxDegs[0];
	int lt_r = res->lt;
	memset(res->elems, 0, sizeof(elem_t)*(lt_r+1));

	// degree_t curDeg = maxDegs[0];
	mpz_t curVal;
	elem_t modVal;
	mpz_init(curVal);
	degree_t termDegs[nvar];

	for (int i = 0; i < aa->size; ++i) {
		mpz_set(curVal, aa->elems[i].coef);
		partialDegreesTerm_AAZ(aa, i, termDegs);
		for (int j = 0; j < nvar-1; ++j) {
			// curDeg = GET_NTH_EXP(aa->elems[i].degs, masks[j+1], sizes[j+1]);
			mpz_mul(curVal, curVal, valList[j][termDegs[j+1]]);
		}
		// curDeg = GET_NTH_EXP(aa->elems[i].degs, masks[0], sizes[0]);
		modVal = smallprimefield_convert_in(mpz_fdiv_ui(curVal, Pptr->prime), Pptr);
		res->elems[termDegs[0]] = smallprimefield_add(res->elems[termDegs[0]], modVal, Pptr);
	}

	for (int j = 0; j < nvar-1; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}
	mpz_clear(curVal);
}



// void interpCoefsWithSkeleton_AAZ(const AltArrZ_t* skeleton )

/**
 * Combine many polynomials which represent the dense representation of
 * a polynomial with polynomial coefficients in a sparse polynomial in one more variable.
 * That is the index of the polynomial in coefs represents the degree of
 * a new variable corresponding to that coefficient. 
 *
 * The polynomials in coefs must all have the same number of variables.
 *
 * @param coefs the dense list of polynomial coefficients
 * @param n the number of coefficients in coefs.
 * 
 * @return a sparse polynomial in coefs[0]->nvar + 1 variables
 *         whose terms are built from coefs. 
 */
AltArrZ_t* combinePolyCoefs(AltArrZ_t const*const* coefs, int n) {
	int finalSize = 0;
	int i;
	int nvar = 0;
	int unpacked = 0;
	for (i = 0; i < n; ++i) {
		finalSize += coefs[i]->size;
		unpacked |= coefs[i]->unpacked;
		nvar |= coefs[i]->nvar;
	}

	if (nvar != coefs[0]->nvar) {
		fprintf(stderr, "ERROR in combinePolyCoefs: Number of variables in coefs not consistent.\n");
		return NULL;
	}

	if (unpacked) {
		AltArrZ_t** tmpCoefs = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*n);
		unpacked = 0;
		for (int i = 0; i < n; ++i) {
			tmpCoefs[i] = deepCopyPolynomial_AAZ(coefs[i]);
			tryPackExponentVectors_AAZ_inp(tmpCoefs[i]);
			unpacked |= tmpCoefs[i]->unpacked;
		}
		if (unpacked) {
			fprintf(stderr, "ERROR in combinePolyCoefs: Unpacked version not yet implemented!\n");
			return NULL;
		}

		AltArrZ_t* ret = combinePolyCoefs(CONSTCONSTCAST(AltArrZ_t,tmpCoefs), n);
		for (int i = 0; i < n; ++i) {
			freePolynomial_AAZ(tmpCoefs[i]);
		}
		free(tmpCoefs);
		return ret;
	}


	const degrees_t* __restrict__ oldMasks = getExpMaskArray(nvar);
	const int* __restrict__ oldSizes = getExpOffsetArray(nvar);
	const int* __restrict__ newSizes = getExpOffsetArray(nvar+1);
	const degrees_t* __restrict__ maxExps = getMaxExpArray(nvar+1);

	AltArrZ_t* aa = makePolynomial_AAZ(finalSize, nvar);
	aa->size = finalSize;
	const AltArrZ_t* curCoef;
	AAZElem_t* curElems = aa->elems;

	int j, k;
	polysize_t t;
	degrees_t curDeg; //use degrees_t for cast-less shift.
	for (i = n-1; i >= 0; --i) {
		curCoef = coefs[i];
		if (curCoef == NULL) {
			continue;
		}

		t = curCoef->size;
		for (j = 0; j < t; ++j) {
			mpz_init_set(curElems[j].coef, curCoef->elems[j].coef);
			curElems[j].degs = 0;
			for (k = 0; k < nvar; ++k) {
				curDeg = GET_NTH_EXP(curCoef->elems[j].degs, oldMasks[k], oldSizes[k]); 
				if (curDeg > maxExps[k+1]) {
					fprintf(stderr, "ERROR in combinePolyCoefs: High degree caused unpack; unpacked version not yet implemented!\n");
					return NULL;
				}
				curElems[j].degs |= curDeg << newSizes[k+1];
			}
			curElems[j].degs |= ((degrees_t) i << newSizes[0]);
		} 
		curElems = curElems + t;
	}

	return aa;
}

/**
 * Combine many univariate polynomials which represent the dense representation of
 * a bivariate polynomial with.
 * That is, the index of the polynomial in coefs represents the degree of
 * a new variable corresponding to that coefficient. 
 *
 * The polynomials in coefs must all have the same number of variables.
 *
 * @param coefs the dense list of polynomial coefficients
 * @param n the number of coefficients in coefs.
 * 
 * @return a sparse polynomial in coefs[0]->nvar + 1 variables
 *         whose terms are built from coefs. 
 */
AltArrZ_t* combineUnivarPolyCoefs(AltArrZ_t const*const* coefs, int n) {
	int finalSize = 0;
	int i;
	int nvar = 0;
	int unpacked = 0;
	for (i = 0; i < n; ++i) {
		if (!isZero_AAZ(coefs[i])) {
			finalSize += coefs[i]->size;
			unpacked |= coefs[i]->unpacked;
			nvar |= coefs[i]->nvar;
		}
	}

	for (int i = 0; i < n; ++i) {
		if (!isZero_AAZ(coefs[i]) && nvar != coefs[i]->nvar) {
			fprintf(stderr, "ERROR in combineUnivarPolyCoefs: Number of variables in coefs not consistent.\n");
			return NULL;
		}
	}

	if (unpacked) {
		//this should never happen but here we are
		AltArrZ_t** tmpCoefs = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*n);
		unpacked = 0;
		for (int i = 0; i < n; ++i) {
			tmpCoefs[i] = deepCopyPolynomial_AAZ(coefs[i]);
			tryPackExponentVectors_AAZ_inp(tmpCoefs[i]);
			unpacked |= tmpCoefs[i]->unpacked;
		}
		if (unpacked) {
			fprintf(stderr, "ERROR in combineUnivarPolyCoefs: Unpacked version not yet implemented!\n");
			return NULL;
		}

		AltArrZ_t* ret = combineUnivarPolyCoefs(CONSTCONSTCAST(AltArrZ_t, tmpCoefs), n);
		for (int i = 0; i < n; ++i) {
			freePolynomial_AAZ(tmpCoefs[i]);
		}
		free(tmpCoefs);
		return ret;
	}

	if (finalSize == 0) {
		return NULL;
	}

	AltArrZ_t* aa = makePolynomial_AAZ(finalSize, nvar+1);
	aa->size = finalSize;
	const AltArrZ_t* curCoef;
	AAZElem_t* curElems = aa->elems;

	int j;
	polysize_t t;
	for (i = n-1; i >= 0; --i) {
		curCoef = coefs[i];
		if (isZero_AAZ(curCoef)) {
			continue;
		}

		t = curCoef->size;
		for (j = 0; j < t; ++j) {
			mpz_init_set(curElems[j].coef, curCoef->elems[j].coef);
			curElems[j].degs = curCoef->elems[j].degs;
			curElems[j].degs |= ((degrees_t) i << EXP_OFFSET_1_V2);
		} 
		curElems = curElems + t;
	}

	return aa;
}



/** 
 * Solves the diophantine equation c = sum_i^r(us[i], sigmas[i]). 
 * where us[i], sigmas[i], c are bivariate polynomials over a finite field.
 * NOTE: assumed sigmas[i] is NULL; 
 */
int multiBDP_AAZ(const AltArrZ_t* c, AltArrZ_t const* const* us, AltArrZ_t** sigmas, unsigned int r, const Prime_ptr* Pptr, const mpz_t mpPrime) {
	if (sigmas == NULL || us == NULL) {
		return 0;
	}

	//get size of prime for rand beta
	elem_t lprime = Pptr->prime;
	int msb = 0; 
    while (lprime != 0) { 
        lprime = lprime / 2; 
        ++msb; 
    } 

	size_t betasAlloc = 10;
	mpz_t* betas = (mpz_t*) malloc(sizeof(mpz_t)*betasAlloc);
	mpz_init(betas[0]);


	rand_mpz_t(msb, 0, betas[0]);
	mpz_mod_ui(betas[0], betas[0], Pptr->prime); //just in case

	duspoly_t* evalC = NULL;
	evaluatePolyToDUSP_AAZ(c, betas, &evalC, Pptr);
	duspoly_t* evalUs[r];
	for (int i = 0; i < r; ++i) {
		if (us[i]->nvar != 2) {
			fprintf(stderr, "In multiBDP_AAZ an input polynomial was not bivariate!\n");
			exit(1);
		}
		evalUs[i] = NULL;
		evaluatePolyToDUSP_AAZ(us[i], betas, evalUs + i, Pptr);
		// gmp_fprintf(stderr, "beta: %Zd\n", betas[0]);
		// printPolynomialOutForm_spX(evalUs[i], Pptr);
	}

	// fprintf(stderr, "evalC: ");
	// printPolynomialOutForm_spX(evalC, Pptr);


	//TODO: probably we can save time if we pre-compute the products of the Us and then evaluate, 
	//passing the products to multiUDP

	duspoly_t** sigmas_sp = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	// fprintf(stderr, "calling multiUDP first time\n");
	int udpPass = multiUDP_spX(evalC, CONSTCONSTCAST(duspoly_t, evalUs), r, sigmas_sp, Pptr);
	if (!udpPass) {
		//TODO try again with a new evaluation point?
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "UDP FAILED in BDP!!\n");
#endif
		return 0;
	} else {
		// fprintf(stderr, "First UDP passed in BDP!!\n");
		// for (int i = 0; i < r; ++i) {
			// fprintf(stderr, "sigmas[%d]: ", i);
			// printPolynomialOutForm_spX(sigmas_sp[i], Pptr);
		// }
	}



	//sigmas interp is a matrix of polynomials.
	//row i is a list of polynomials associated with the the ith true sigma.
	//this list has one polynomial per coefficient of the true sigma and its
	//column indicates the degree of the main variable the coefficient poly
	//is associated with. 
	AltArrZ_t*** sigmas_interp = (AltArrZ_t***) malloc(sizeof(AltArrZ_t**) * r);
	int* sigmas_sizes = (int*) malloc(sizeof(int)*r);

	int curDeg = 0;
	mpz_t curCoef;
	mpz_init(curCoef);
	for (int i = 0; i < r; ++i) {
		if (sigmas_sp[i] == NULL) {
			sigmas_sizes[i] = 1;
			sigmas_interp[i] = (AltArrZ_t**) calloc(sigmas_sizes[i], sizeof(AltArrZ_t**));
			mpz_set_ui(curCoef, 0ul);
			univarIncrNewtonInterpModP_AAZ(sigmas_interp[i], betas, curDeg, curCoef, Pptr);
			// fprintf(stderr, "sigmas_interp[%d][0] = ", i);
			// printPoly_AAZ(stderr, sigmas_interp[i][0], syms, 1);
			// fprintf(stderr, "\n");
		} else {
			sigmas_sizes[i] = sigmas_sp[i]->lt + 1;		
			sigmas_interp[i] = (AltArrZ_t**) calloc(sigmas_sizes[i], sizeof(AltArrZ_t**));
			for (int k = 0; k < sigmas_sizes[i]; ++k) {
				mpz_set_ui(curCoef, smallprimefield_convert_out(sigmas_sp[i]->elems[k], Pptr));
				univarIncrNewtonInterpModP_AAZ(sigmas_interp[i] + k, betas, curDeg, curCoef, Pptr);
				// fprintf(stderr, "sigmas_interp[%d][%d] = ", i, k);
				// printPoly_AAZ(stderr, sigmas_interp[i][k], syms, 1);
				// fprintf(stderr, "\n");
			}
			freePolynomial_spX(sigmas_sp + i);
		}
		sigmas[i] = NULL;
	}
	++curDeg;

	AltArrZ_t* tmpProd = makePolynomial_AAZ(c->size, 2);

	int maxTries = 100;
	while(curDeg < maxTries) {
		if (curDeg >= betasAlloc) {
			betasAlloc <<= 1;
			betas = (mpz_t*) realloc(betas, sizeof(mpz_t)*betasAlloc);
		}

		mpz_init(betas[curDeg]);
		rand_mpz_t(msb, 0, betas[curDeg]);
		mpz_mod_ui(betas[curDeg], betas[curDeg], Pptr->prime);
		evaluatePolyToDUSP_AAZ(c, betas+curDeg, &evalC, Pptr);
		// gmp_fprintf(stderr, "beta_%d = %Zd\n", curDeg, betas[curDeg]);
		for (int i = 0; i < r; ++i){
			evaluatePolyToDUSP_AAZ(us[i], betas+curDeg, evalUs + i, Pptr);
			// fprintf(stderr, "u[%d](y = beta_i) = ", i);
			// printPolynomialOutForm_spX(evalUs[i], Pptr);
		}
		
		udpPass = multiUDP_spX(evalC, CONSTCONSTCAST(duspoly_t, evalUs), r, sigmas_sp, Pptr);
		if (!udpPass) {
			//TODO try again with a new evaluation point??
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr, "UDP FAILED in BDP!! %d\n", curDeg);
#endif
			return 0;
		}
		// for (int i = 0; i < r; ++i) {
			// fprintf(stderr, "sigmas[%d]: ", i);
			// printPolynomialOutForm_spX(sigmas_sp[i], Pptr);
		// }
		//now add new points to interpolation of each coef.
		int changed = 0;
		for (int i = 0; i < r; ++i) {
			if ((sigmas_sp[i] == NULL && sigmas_sizes[i] != 1) || 
				(sigmas_sp[i] != NULL && sigmas_sizes[i] != sigmas_sp[i]->lt + 1)) {
				for (int k = 0; k < r; ++k) {
					for (int k2 = 0; k2 < sigmas_sizes[k]; ++k2) {
						free(sigmas_interp[k][k2]);
					}
					free(sigmas_interp[k]);
					free(sigmas_sp[k]);
				}
				free(sigmas_interp);
				free(sigmas_sp);
				for (int k = 0; k <= curDeg; ++k) {
					mpz_clear(betas[k]);
				}
				free(betas);
				return 0; //FAIL
			}
			for (int k = 0; k < sigmas_sizes[i]; ++k) {
				if (sigmas_sp[i] == NULL) {
					mpz_set_ui(curCoef, 0l);
				} else {
					mpz_set_ui(curCoef, smallprimefield_convert_out(sigmas_sp[i]->elems[k], Pptr));
				}

				//if at least one sigmas coef gets changed then changed will be non-zero.
				changed |= univarIncrNewtonInterpModP_AAZ(sigmas_interp[i] + k, betas, curDeg, curCoef, Pptr);
				// fprintf(stderr, "sigmas_interp[%d][%d] = ", i, k);
				// printPoly_AAZ(stderr, sigmas_interp[i][k], syms, 1);
				// fprintf(stderr, "\n");
			}
		}
		// fprintf(stderr, "\n");


		if (!changed) {
			//combine interpolated coefs;
			for (int i = 0; i < r; ++i) {
				sigmas[i] = combineUnivarPolyCoefs(CONSTCONSTCAST(AltArrZ_t, sigmas_interp[i]), sigmas_sizes[i]);
			}


			//now check if correct;
			AltArrZ_t* sum = NULL;
			AltArrZ_t* b;
			for (int i = 0; i < r; ++i) {
				b = multiplyAllButOnePolynomials_AAZ(us, r, i);
				// fprintf(stderr, "b->nvar: %d\n", b == NULL ? 0 : b->nvar);
				// fprintf(stderr, "sigmas[%d]->nvar: %d\n", i, sigmas[i] == NULL ? 0 : sigmas[i]->nvar);
				// fprintf(stderr, "b: \n");
				// printPoly_AAZ(stderr, b, henselsyms, b->nvar);
				// fprintf(stderr, "\nsigmas[%d]:\n", i);
				// printPoly_AAZ(stderr, sigmas[i], henselsyms, sigmas[i]->nvar);
				multiplyPolynomialsPreAlloc_AAZ(sigmas[i], b, &tmpProd);
				sum = addPolynomials_AAZ_inp(sum, tmpProd, 2);
				freePolynomial_AAZ(b);
			}

			applyModuloSymmetric_AAZ_inp(sum, mpPrime);
			// fprintf(stderr, "\nsum mod:\n" );
			// printPoly_AAZ(stderr, sum, henselsyms, sum->nvar);
			// fprintf(stderr, "\n\n");

			int done = isExactlyEqual_AAZ(c, sum);
			freePolynomial_AAZ(sum);

			if (done) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr, "\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nmultiBDP DONE!\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n" );
#endif
				break;
			} else {
				for (int i = 0; i < r; ++i) {
					freePolynomial_AAZ(sigmas[i]);
					sigmas[i] = NULL;
				}
			}
		}

		++curDeg;
		// fprintf(stderr, "curDeg in multiBDP loop: %ld\n", curDeg);
	}

	//cleanup
	for (int k = 0; k < r; ++k) {
		for (int k2 = 0; k2 < sigmas_sizes[k]; ++k2) {
			free(sigmas_interp[k][k2]);
		}
		free(sigmas_interp[k]);
		free(sigmas_sp[k]);
	}
	free(sigmas_interp);
	free(sigmas_sp);
	for (int k = 0; k < curDeg; ++k) {
		mpz_clear(betas[k]);
	}
	free(betas);
	freePolynomial_AAZ(tmpProd);
	
	return (curDeg < maxTries);
}


//liftedF already have true LC multipled into them.
int bivarHenselLiftCoefs_AAZ(const AltArrZ_t* a, AltArrZ_t** Fs, unsigned int r, const mpz_t bound, const Prime_ptr* Pptr) {

	if (Fs == NULL) {
		return 0;
	}


	AltArrZ_t* error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, Fs), r);
	// fprintf(stderr, "error: \n" );
	// printPoly_AAZ(stderr, error, henselsyms, error->nvar);
	subPolynomials_AAZ_inpRHS(a, &error);	 
	// fprintf(stderr, "\nerror after: \n" );
	// printPoly_AAZ(stderr, error, henselsyms, error == NULL ? 0 : error->nvar);
	// fprintf(stderr, "\n" );

	if (isZero_AAZ(error)) {
		// fprintf(stderr, "No lifting to do in lift coef bivar\n");
		freePolynomial_AAZ(error);
		return 1;
	}

	//alloc for bezout coefficients
	AltArrZ_t** sigmas = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*r);

	mpz_t mod; 
	mpz_init_set_ui(mod, Pptr->prime);
	mpz_t mpPrime;
	mpz_init_set(mpPrime, mod);
	int primePow = 1;

	int success = 0;
	while( !isZero_AAZ(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
		divideByIntegerExact_AAZ_inp(error, mod);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);
		
		success = 0;
		for (int nsteps = 0; nsteps < 5 && !success; ++nsteps) {
			//try multiBDP a few times 
			success = multiBDP_AAZ(error, CONSTCONSTCAST(AltArrZ_t, Fs), sigmas, r, Pptr, mpPrime);
		}

		if (!success) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr, "in bivarLiftCoefs: multi BDP FAILED!!! for primePow=%d\n", primePow);
#endif
			break;
		}

		for (int i = 0; i < r; ++i) {
			multiplyByInteger_AAZ_inp(sigmas[i], mod);
			Fs[i] = addPolynomials_AAZ_inp(Fs[i], sigmas[i], Fs[i]->nvar);
			freePolynomial_AAZ(sigmas[i]);
		}

		//update the error.
		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, Fs), r, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}

	if (!isZero_AAZ(error)) {
		//we got here because bound exceeded.
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		gmp_fprintf(stderr, "Bound exceeded in bivarHenselLiftCoefs_AAZ, current mod: %Zd\n", mod);
#endif
	}

	free(sigmas);	

	freePolynomial_AAZ(error);
	mpz_clears(mod, mpPrime, NULL);

	return primePow;
}



AltArrZ_t* replaceLC_DUZP(const DUZP_t* p, const AltArrZ_t* lc, mpz_t mpPrime) {

	if (p == NULL || lc == NULL)
	 {
		return NULL;
	}

	if (isConstant_AAZ(lc)) {
		AltArrZ_t* ret = convertToAltArrZ_DUZP(p);
		expandNumVars_AAZ(ret, lc->nvar);
		mpz_set(ret->elems->coef, lc->elems->coef);
		return ret;
	}

	int nvar = lc->nvar;

	AltArrZ_t* a = deepCopyPolynomial_AAZ(lc);
	// expandNumVarsLeft_AAZ(a, nvar); //lc already bivar
	resizePolynomial_AAZ(a, p->lt + lc->size + 2);

	degree_t deg = p->lt;
	AAZElem_t* elems = a->elems;
	for (int i = 0; i < lc->size; ++i) {
		setExponentTerm_AAZ_inp(a, i, deg, 0);
		// elems[i].degs |= ((degrees_t) deg << offset);
	}

	degree_t degs[nvar];
	for (int i = 1; i < nvar; ++i) {
		degs[i] = 0;
	}
	int aIdx = lc->size;
	for (int i = p->lt-1; i >= 0; --i) {
		mpz_init_set(elems[aIdx].coef, p->coefs[i]);
		degs[0] = i;
		setDegrees_AAZ_inp(a, aIdx, degs, nvar);
		// elems[aIdx].degs = ((degrees_t) i << offset);
		++aIdx;
	}

	a->size = aIdx;
	applyModuloSymmetric_AAZ_inp(a, mpPrime); //apply modulo will get rid of terms with 0 coefs since 0 mod p is 0.
	return a;
}

void replaceLC_AAZ_inp_unpk(AltArrZ_t** pp, const AltArrZ_t* lc_in) {
	if (pp == NULL || lc_in == NULL) {
		return;
	}

	if (*pp == NULL) {
		*pp = deepCopyPolynomial_AAZ(lc_in);
		return;
	}

	AltArrZ_t* p = *pp;
	if (!p->unpacked) {
		unpackExponentVectors_AAZ_inp(p);
	}

	const AltArrZ_t* lc = lc_in;
	if (!lc_in->unpacked) {
		lc = deepCopyPolynomial_AAZ(lc_in);
		unpackExponentVectors_AAZ_inp((AltArrZ_t*) lc);
	}


	degree_t* degs = (degree_t*) p->elems->degs;
	int nvar = p->nvar;

	int recLCSize = 1;
	degrees_t mdeg = degs[0];

	int i = 0;
	while (++i < p->size && (degs[i*nvar] == mdeg)) {
		++recLCSize;
	}

	int newSize = p->size + lc->size - recLCSize;
	degree_t* lcDegs = (degree_t*) lc->elems->degs;

	if (recLCSize < lc->size) {
		//allocate space, shift terms right, and then replace the LC 

		resizePolynomial_AAZ(p, newSize);
		degs = (degree_t*) p->elems->degs;
		for (int i = p->size; i < newSize; ++i) {
			mpz_init(p->elems[i].coef);
			p->elems[i].degs = (degrees_t) degs + i*nvar;
		}

		memmove(degs+lc->size, degs+recLCSize, sizeof(degree_t)*nvar*(p->size - recLCSize));


		for (int i = p->size-1; i >= recLCSize; ++i) {
			mpz_set(p->elems[i - recLCSize + lc->size].coef, p->elems[i].coef);
		}

		memcpy(degs, lcDegs, sizeof(degree_t)*nvar*lc->size);
		for (int i = 0; i < lc->size; ++i) {
			mpz_set(p->elems[i].coef, lc->elems[i].coef);
			degs[i*nvar] = mdeg;
		}
	}
	else if (lc->size < recLCSize) {
		//shift terms left, free un-needed space, and then replace LC

		memmove(degs+lc->size, degs+recLCSize, sizeof(degree_t)*nvar*(p->size - recLCSize));
		for (int i = lc->size; i < p->size; ++i) {
			mpz_set(p->elems[i].coef, p->elems[i + recLCSize - lc->size].coef);
		}

		memcpy(degs, lcDegs, sizeof(degree_t)*nvar*lc->size);
		for (int i = 0; i < lc->size; ++i) {
			mpz_set(p->elems[i].coef, lc->elems[i].coef);
			degs[i*nvar] = mdeg;
		}

	} else {
		memcpy(degs, lcDegs, sizeof(degree_t)*nvar*lc->size);
		for (int i = 0; i < lc->size; ++i) {
			mpz_set(p->elems[i].coef, lc->elems[i].coef);
			degs[i*nvar] = mdeg;
		}
	}

	p->size = newSize;

	if (!lc_in->unpacked) {
		freePolynomial_AAZ((AltArrZ_t*) lc);
	}
}

void replaceLC_AAZ_inp(AltArrZ_t** pp, const AltArrZ_t* lc) {

	if (pp == NULL || lc == NULL) {
		return;
	}

	if (*pp == NULL) {
		*pp = deepCopyPolynomial_AAZ(lc);
		return;		
	}

	AltArrZ_t* p = *pp;

	if (p->unpacked || lc->unpacked) {
		replaceLC_AAZ_inp_unpk(pp, lc);
		return;
	}

	int nvar = p->nvar;

	int recLCSize = 1;
	degrees_t mvarMask = getMVarExpMask(nvar);
	degrees_t mdeg = p->elems->degs & mvarMask;
	// fprintf(stderr, "mdeg: %llx\n", mdeg);
	int i = 0;
	while (++i < p->size && ( (p->elems[i].degs & mvarMask) == mdeg)) {
		// fprintf(stderr, "degs[%d]: %llx\n", i, p->elems[i].degs);
		++recLCSize;
	}

	// fprintf(stderr, "recLCSize: %d\n",recLCSize );
	int idxLC = 0;
	if (recLCSize < lc->size) {
		resizePolynomial_AAZ(p, (p->size + lc->size - recLCSize));
		//shift terms right
		memmove(p->elems+lc->size, p->elems+recLCSize, sizeof(AAZElem_t)*(p->size-recLCSize));
		for (i = 0; i < recLCSize; ++i) {
			mpz_set(p->elems[i].coef, lc->elems[idxLC].coef);
			p->elems[i].degs = mdeg | lc->elems[idxLC++].degs;
		}
		for ( ; i < lc->size; ++i) {
			mpz_init_set(p->elems[i].coef, lc->elems[idxLC].coef);
			p->elems[i].degs = mdeg | lc->elems[idxLC++].degs;
		}
	} else if (recLCSize > lc->size) {
		for (i = lc->size; i < recLCSize; ++i) {
			mpz_clear(p->elems[i].coef);
		}
		//shift terms left after clearing overlapping terms
		memmove(p->elems+lc->size, p->elems+recLCSize, sizeof(AAZElem_t)*(p->size-recLCSize));
		for (i = 0; i < lc->size; ++i) {
			mpz_set(p->elems[i].coef, lc->elems[idxLC].coef);
			p->elems[i].degs = mdeg | lc->elems[idxLC++].degs;
		}
	} else {
		for (i = 0; i < recLCSize; ++i) {
			mpz_set(p->elems[i].coef, lc->elems[idxLC].coef);
			p->elems[i].degs = mdeg | lc->elems[idxLC++].degs;
		}
	}
	p->size = p->size - recLCSize + lc->size;

}

//Assuming p has x_1...x_j and lc has x_1,x_2,...x_n, j < n, replace lc(p) with lc.
AltArrZ_t* replaceLC_AAZ(const AltArrZ_t* p, const AltArrZ_t* lc) {
	AltArrZ_t* ret = deepCopyPolynomial_AAZ(p);
	expandNumVars_AAZ(ret, lc->nvar); //p is missing at least the least variable of lc.
	replaceLC_AAZ_inp(&ret, lc);
	return ret;
}


int trivarHenselLiftVars_AAZ(const AltArrZ_t* aa, AltArrZ_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t* x, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	if (aa == NULL || lcF == NULL || fs == NULL || liftedF == NULL || r < 2 || x == NULL) {
		return 0;
	}

	mpz_t mpPrime;
	mpz_init_set_si(mpPrime, Pptr->prime);

	short nvar = aa->nvar;
	int active[3] = {0, 0, 1}; //evaluate a and LC's to bivar first

	AltArrZ_t* aMod = deepCopyPolynomial_AAZ(aa);
	applyModuloSymmetric_AAZ_inp(aMod, mpPrime);
	
	AltArrZ_t* a = evaluatePoly_AAZ(aMod, active, x, nvar);
	applyModuloSymmetric_AAZ_inp(a, mpPrime);
	
	duspoly_t** fps = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	AltArrZ_t** tmpLifted = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*) * r);
	for (unsigned int i = 0; i < r; ++i) {
		fps[i] = convertToDUSP_DUZP(fs[i], Pptr);
		AltArrZ_t* tmpLC = evaluatePoly_AAZ(lcF[i], active, x, nvar); 
		tmpLifted[i] = replaceLC_DUZP(fs[i], tmpLC, mpPrime);
		freePolynomial_AAZ(tmpLC);
	}

	AltArrZ_t* error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, tmpLifted), r);
	subPolynomials_AAZ_inpRHS(a, &error);
	applyModuloSymmetric_AAZ_inp(error, mpPrime);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
	fprintf(stderr, "\ngot first error: ");
	printPoly_AAZ(stderr, error, henselsyms, error->nvar);
	fprintf(stderr, "\n\n");
#endif


	degree_t deg = partialDegree_AAZ(a, 1);
	AltArrZ_t* kDeriv;
	AltArrZ_t* c;
	duspoly_t* cp = makePolynomial_spX(a->size + 1);

	mpz_t kfact;
	mpz_init_set_ui(kfact, 1ul);
	active[1] = 1; //evaluate error term to univar

	for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {
		kDeriv = derivative_AAZ(error, 1, k);

		c = evaluatePoly_AAZ(kDeriv, active, x, 2);
		freePolynomial_AAZ(kDeriv);
		divideByIntegerExact_AAZ_inp(c, kfact); //divide by k!		
		univarToPrimeFieldPreAlloc_AAZ(c, Pptr, &cp);
		freePolynomial_AAZ(c);

		multiUDP_spX(cp, CONSTCONSTCAST(duspoly_t, fps), r, sigmas, Pptr);

		for (unsigned int i = 0; i < r; ++i) {
			idealadicUpdate(tmpLifted + i, sigmas[i], x[1], k, Pptr, mpPrime);
			freePolynomial_spX (&sigmas[i]);
		}


		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, tmpLifted), r, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);

		mpz_mul_ui(kfact, kfact, k+1ul);
	}	

	//cleanup and check success before next round of lifting.
	for (unsigned int i = 0; i < r; ++i) {
		freePolynomial_spX (&fps[i]);
	}
	free(fps);	
	
	freePolynomial_AAZ(a); //can just use aMod now
	a = NULL;
	freePolynomial_spX(&cp);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
	fprintf(stderr, "Checking error after first round of lift\n");
	printPoly_AAZ(stderr, error, henselsyms, error->nvar);
	fprintf(stderr, "\n\n" );
#endif

	if (!isZero_AAZ(error)) {
		//lift from univar to bivar failed!
		freePolynomial_AAZ(error);
		freePolynomial_AAZ(aMod);
		mpz_clear(kfact);
		free(sigmas);
		return 0;
	}

	//reset everything now for bivar lift
	mpz_set_ui(kfact, 1ul);
	active[1] = 0;
	active[2] = 1; //evaluate error term to bivar

	for (unsigned int i = 0; i < r; ++i) {
		//replaceLC also expands nvar by 1.
		// fprintf(stderr, "before replace %d\n",i );
		// printPoly_AAZ(stderr, tmpLifted[i], henselsyms, tmpLifted[i]->nvar);
		// fprintf(stderr, "\n" );
		// fprintf(stderr, "lc %d\n",i );
		// printPoly_AAZ(stderr, lcF[i], henselsyms, lcF[i]->nvar);
		// fprintf(stderr, "\n" );
		
		liftedF[i] = replaceLC_AAZ(tmpLifted[i], lcF[i]);
		// fprintf(stderr, "after replace %d\n",i );
		// printPoly_AAZ(stderr, liftedF[i], henselsyms, liftedF[i]->nvar);
		// fprintf(stderr, "\n" );
		applyModuloSymmetric_AAZ_inp(liftedF[i], mpPrime);
	}

	deg = partialDegree_AAZ(aMod, 2);

	AltArrZ_t** sigmas_aaz;
	if (sizeof(AltArrZ_t*) == sizeof(duspoly_t*)) {
		sigmas_aaz = (AltArrZ_t**) sigmas;
	} else {
		free(sigmas);
		sigmas_aaz = (AltArrZ_t**) malloc(sizeof(AltArrZ_t)*r);
	}
	sigmas = NULL;


#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
	fprintf(stderr, "Finished first reset in trivar var lift\n");
#endif


	error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r);

	subPolynomials_AAZ_inpRHS(aMod, &error);	 
	applyModuloSymmetric_AAZ_inp(error, mpPrime);
	int dioSuccess = 0;
	for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {

		//TODO use improved version of getting taylor coef, avoid computing full derivative.
		//see Monogan and Pearce: POLY: New dast struct in Maple 17.
		//derive and evaluate w.r.t to third var
		kDeriv = derivative_AAZ(error, 2, k);
		c = evaluatePoly_AAZ(kDeriv, active, x, 3);

		divideByIntegerExact_AAZ_inp(c, kfact);
		applyModuloSymmetric_AAZ_inp(c, mpPrime);
		freePolynomial_AAZ(kDeriv);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "c for k = %d \n", k);
		printPoly_AAZ(stderr, c, henselsyms, c->nvar);
		fprintf(stderr, "\n" );
		for (int l = 0; l < r; ++l){
			fprintf(stderr, "tmpLifted[%d]: \n", l);
			printPoly_AAZ(stderr, tmpLifted[l], henselsyms, tmpLifted[l]->nvar);
			fprintf(stderr, "\n");
		}
#endif

		dioSuccess = multiBDP_AAZ(c, CONSTCONSTCAST(AltArrZ_t, tmpLifted), sigmas_aaz, r, Pptr, mpPrime);
		freePolynomial_AAZ(c);

		if (!dioSuccess) {
			for (unsigned int i = 0; i < r; ++i) {
				freePolynomial_AAZ(sigmas_aaz[i]);
			}
			break;
		}

		for (unsigned int i = 0; i < r; ++i) {
			multivarIdealadicUpdate(liftedF + i, sigmas_aaz[i], x[2], k, mpPrime);
			freePolynomial_AAZ(sigmas_aaz[i]);
		}

		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r, &error);
		subPolynomials_AAZ_inpRHS(aMod, &error);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "\n\nupdated error: \n" );
		printPoly_AAZ(stderr, error, henselsyms, error == NULL ? 0 : error->nvar);
		fprintf(stderr, "\n\n" );
#endif
		mpz_mul_ui(kfact, kfact, k+1ul);
	}

	for (unsigned int i = 0; i < r; ++i) {
		freePolynomial_AAZ(tmpLifted[i]);
	}
	free(tmpLifted);
	free(sigmas_aaz);

	mpz_clear(kfact);

	int ret = isZero_AAZ(error);
	freePolynomial_AAZ(error);
	freePolynomial_AAZ(aMod);

	return ret;

}

typedef struct RecBivarElem {
	degrees_t degs12;
	int coefSize;
	AAZElem_t* coef;
} RecBivarElem_t;

typedef struct RecBivarAAZ {
	int totalNvar;
	int totalSize;
	int size;
	int unpacked;
	RecBivarElem_t* elems;
} RecBivarAAZ_t;

RecBivarAAZ_t* convertToRecBivar_AAZ_unpk(const AltArrZ_t* aa) {
	if (isZero_AAZ(aa)) {
		return NULL;
	}

	AltArrZ_t* a = deepCopyPolynomial_AAZ(aa);
	RecBivarAAZ_t* rec = (RecBivarAAZ_t*) malloc(sizeof(RecBivarAAZ_t));
	rec->totalNvar = a->nvar;
	rec->unpacked = 1;

	int size = 1;
	int nvar = a->nvar;
	degree_t* degs = (degree_t*) a->elems->degs;
	degree_t curX = degs[0];
	degree_t curY = degs[1];
	for (int i = 1; i < a->size; ++i){
		if (degs[i*nvar] != curX || degs[i*nvar + 1] != curY) {
			++size;
			curX = degs[i*nvar];
			curY = degs[i*nvar + 1];
		}
	}

	RecBivarElem_t* recElems = (RecBivarElem_t*) malloc(sizeof(RecBivarElem_t)*size);
	rec->size = size;
	rec->totalSize = a->size;
	rec->elems = recElems;

	curX = degs[0];
	curY = degs[1];
	degs[0] = 0;
	degs[1] = 1;
	recElems->degs12 = (((degrees_t) curX) << 32) | curY;
	recElems->coef = a->elems;
	int curRec = 0;
	recElems[curRec].coefSize = 1;
	for (int i = 1; i < a->size; ++i){
		if (degs[i*nvar] != curX || degs[i*nvar + 1] != curY) {
			++curRec;
			curX = degs[i*nvar];
			curY = degs[i*nvar + 1];
			recElems[curRec].degs12 = (((degrees_t) curX) << 32) | curY;
			recElems[curRec].coefSize = 1;
			recElems[curRec].coef = a->elems + i;
		} else {
			++(recElems[curRec].coefSize);
		}

		degs[i*nvar] = 0;
		degs[i*nvar + 1] = 0;
	}	

	//yes free a only don't use freePolynomial_AAZ since rec is using its AAZElem_t
	free(a);

	return rec;
}

RecBivarAAZ_t* convertToRecBivar_AAZ(const AltArrZ_t* aa, degrees_t biMask, const degrees_t* masks, const int* offsets) {
	if (isZero_AAZ(aa)) {
		return NULL;
	}

	if (aa->unpacked) {
		return convertToRecBivar_AAZ_unpk(aa);
	}

	AltArrZ_t* a = deepCopyPolynomial_AAZ(aa);
	RecBivarAAZ_t* rec = (RecBivarAAZ_t*) malloc(sizeof(RecBivarAAZ_t));
	rec->totalNvar = a->nvar;
	rec->unpacked = 0;

	int size = 1;
	degrees_t curMon = a->elems->degs & biMask;
	for (int i = 1; i < a->size; ++i){
		if (curMon != (a->elems[i].degs & biMask)) {
			++size;
			curMon = a->elems[i].degs & biMask;
		}
	}

	RecBivarElem_t* recElems = (RecBivarElem_t*) malloc(sizeof(RecBivarElem_t)*size);
	rec->size = size;
	rec->totalSize = a->size;
	rec->elems = recElems;

	curMon = a->elems->degs & biMask;
	a->elems[0].degs &= ~biMask;
	recElems->degs12 = ((degrees_t) GET_NTH_EXP(curMon, masks[0], offsets[0])) << EXP_OFFSET_1_V2 | GET_NTH_EXP(curMon, masks[1], offsets[1]);
	recElems->coef = a->elems;
	int curRec = 0;
	recElems[curRec].coefSize = 1;
	for (int i = 1; i < a->size; ++i){
		if (curMon != (a->elems[i].degs & biMask)) {
			++curRec;
			curMon = a->elems[i].degs & biMask;
			recElems[curRec].degs12 = ((degrees_t) GET_NTH_EXP(curMon, masks[0], offsets[0])) << EXP_OFFSET_1_V2 | GET_NTH_EXP(curMon, masks[1], offsets[1]);
			recElems[curRec].coefSize = 1;
			recElems[curRec].coef = a->elems + i;
		} else {
			++(recElems[curRec].coefSize);
		}

		a->elems[i].degs &= ~biMask; //clear out all exponents not in top 2, since those are cached in degs12.
	}	

	//yes free a only don't use freePolynomial_AAZ since rec is using its AAZElem_t
	free(a);

	return rec;
}


#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
void printRecBivar_AAZ(RecBivarAAZ_t* rec, const char** syms) {
	AltArrZ_t coef;
	coef.unpacked = 0;
	coef.nvar = rec->totalNvar;

	for(int i = 0; i < rec->size; ++i) {
		coef.size = coef.alloc = rec->elems[i].coefSize;
		coef.elems = rec->elems[i].coef;
		fprintf(stderr, "(" );
		printPoly_AAZ(stderr, &coef, syms, coef.nvar);
		fprintf(stderr, ")*%s%s^(%llx)", syms[0], syms[1], rec->elems[i].degs12);
		if (i + 1 < rec->size) {
			fprintf(stderr, " + ");
		}
	}
	fprintf(stderr, "\n");
}
#endif

void freeRecBivar_AAZ(RecBivarAAZ_t* rec) {
	if (rec == NULL) {
		return;
	}

	AAZElem_t* elems = rec->elems->coef;
	for (int i = 0; i < rec->totalSize; ++i) {
		mpz_clear(elems[i].coef);
	}

	if (rec->unpacked) {
		degree_t* degs = (degree_t*) elems->degs;
		free(degs);
	}

	free(elems);
	free(rec->elems);
	free(rec);
}

void countBivarSupportCoefSize_AAZ_unpk(const AltArrZ_t* a, unsigned int* supportSize_p, unsigned int* maxCoefSize_p ) {
	if (isZero_AAZ(a)) {
		if (supportSize_p != NULL) {
			*supportSize_p = 0;
		}
		if (maxCoefSize_p != NULL) {
			*maxCoefSize_p = 0;
		}
		return;
	}

	int size = a->size;
	int supportSize = 1;
	int maxCoefSize = 1;
	int curCoefSize = 1;
	int nvar = a->nvar;

	degree_t* degs = (degree_t*) a->elems->degs;
	degree_t curX = degs[0];
	degree_t curY = degs[1];
	for (int i = 1; i < size; ++i) {
		if (degs[i*nvar] != curX || degs[i*nvar + 1] != curY) {
			++supportSize;
			curX = degs[i*nvar];
			curY = degs[i*nvar + 1];
			maxCoefSize = maxCoefSize < curCoefSize ? curCoefSize : maxCoefSize;
			curCoefSize = 1;
		} else {
			++curCoefSize;
		}
	}
	//one more check for the last coef 
	maxCoefSize = maxCoefSize < curCoefSize ? curCoefSize : maxCoefSize;


	if (supportSize_p != NULL) {
		*supportSize_p = supportSize;	
	}
	if (maxCoefSize_p != NULL) {
		*maxCoefSize_p = maxCoefSize;
	}
}

//viewing the polynomial recursively as a bivariate polynomial, 
//count the number of monomials in the bivariate support and the
//max number of terms in any coeffient (which is a poly in x_3...x_nvar)
void countBivarSupportCoefSize_AAZ(const AltArrZ_t* a, degrees_t biMask, unsigned int* supportSize_p, unsigned int* maxCoefSize_p ) {
	if (isZero_AAZ(a)) {
		if (supportSize_p != NULL) {
			*supportSize_p = 0;
		}
		if (maxCoefSize_p != NULL) {
			*maxCoefSize_p = 0;
		}
		return;
	}

	if (a->unpacked) {
		countBivarSupportCoefSize_AAZ_unpk(a, supportSize_p, maxCoefSize_p);
		return;
	}

	int size = a->size;
	int supportSize = 1;
	int maxCoefSize = 1;
	int curCoefSize = 1;
	degrees_t curMon = a->elems->degs & biMask;
	for (int i = 1; i < size; ++i){
		if (curMon != (a->elems[i].degs & biMask)) {
			++supportSize;
			curMon = a->elems[i].degs & biMask;
			maxCoefSize = maxCoefSize < curCoefSize ? curCoefSize : maxCoefSize;
			curCoefSize = 1;
		} else {
			++curCoefSize;
		}
	}
	//one more check for the last coef 
	maxCoefSize = maxCoefSize < curCoefSize ? curCoefSize : maxCoefSize;


	if (supportSize_p != NULL) {
		*supportSize_p = supportSize;	
	}
	if (maxCoefSize_p != NULL) {
		*maxCoefSize_p = maxCoefSize;
	}
}


void evalRecMonomials(const RecBivarAAZ_t* rec, int recIdx, const elem_t* point, elem_t* evals, int nvar, const degrees_t* masks, const int* offsets, const Prime_ptr* Pptr) {

 	const RecBivarElem_t* recElem = &rec->elems[recIdx];
	const AAZElem_t* elems = recElem->coef;
	degree_t deg;

	if (rec->unpacked) {
		degree_t* degs = (degree_t*) elems->degs;
		for (int i = 0; i < recElem->coefSize; ++i) {
			evals[i] = smallprimefield_convert_in(1l, Pptr);
			//we start at 2 because the coef of RecBivar has zero degree for j={0,1};
			for (int j = 2; j < nvar; ++j) {
				//TODO unpacked
				deg = degs[i*nvar + j];
				//point has index 0 for variable index 2
				evals[i] = smallprimefield_mul(evals[i], smallprimefield_exp(point[j-2], deg, Pptr), Pptr);
			}		
		}
	} else {
		for (int i = 0; i < recElem->coefSize; ++i) {
			evals[i] = smallprimefield_convert_in(1l, Pptr);
			//we start at 2 because the coef of RecBivar has zero degree for j={0,1};
			for (int j = 2; j < nvar; ++j) {
				//TODO unpacked
				deg = GET_NTH_EXP(elems[i].degs, masks[j], offsets[j]);
				//point has index 0 for variable index 2
				evals[i] = smallprimefield_mul(evals[i], smallprimefield_exp(point[j-2], deg, Pptr), Pptr);
			}
		}
	}
}



int uniqueElemTSet(elem_t* set, int setSize) {
	for (int k = 0; k < setSize; ++k) {
		for (int l = k+1; l < setSize; ++l) {
			if (set[k] == set[l]) {
				return 0;
			}
		}
	}
	return 1;
}

void evalToBivarSymmetricSPF_AAZ_unpk(AltArrZ_t* ret, const AltArrZ_t* a, elem_t const*const* pointPowers, int j, const Prime_ptr* Pptr) {
	if (ret == NULL || isZero_AAZ(a)) {
		return;
	}


	//Allocate a little bit extra in place of iterating though a to get true bivaraite support ammount.
	AAZElem_t* resElems = ret->elems;
	int size = a->size;
	int nvar = a->nvar;
	int k = 0;
	int v;
	degree_t* degs = (degree_t*) a->elems->degs;
	degree_t curX = degs[0];
	degree_t curY = degs[1];
	elem_t curCoef = smallprimefield_convert_in(0, Pptr);
	elem_t tmp;
	unsigned long lPrime = Pptr->prime;
	unsigned long halfP = lPrime / 2;
	for (int i = 0; i < size; ++i) {
	if (degs[i*nvar] != curX || degs[i*nvar + 1] != curY) {
			//commit previous result to ret;
			tmp = smallprimefield_convert_out(curCoef, Pptr);
			if (tmp > halfP) {
				tmp -= lPrime;
			}
			mpz_set_si(resElems[k].coef, tmp);
			resElems[k].degs = 0ull;
			resElems[k].degs |= ( ((degrees_t) curX) << EXP_OFFSET_1_V2);
			resElems[k].degs |= ( ((degrees_t) curY) << EXP_OFFSET_2_V2);
			++k;
		
			curCoef = smallprimefield_convert_in(0, Pptr);
			curX = degs[i*nvar];
			curY = degs[i*nvar + 1];
		}

		tmp = smallprimefield_convert_in(1l, Pptr);
		for (v = 2; v < nvar; ++v) {
			tmp = smallprimefield_mul(tmp, pointPowers[v-2][j* degs[i*nvar + v] ], Pptr);
		}
		tmp = smallprimefield_mul(tmp, smallprimefield_convert_in(mpz_fdiv_ui(a->elems[i].coef, lPrime), Pptr), Pptr);
		curCoef = smallprimefield_add(curCoef, tmp, Pptr);
	}


	tmp = smallprimefield_convert_out(curCoef, Pptr);
	if (tmp > halfP) {
		tmp -= lPrime;
	}
	//one more to commit the final term since it wasn't triggered in the loop
	mpz_set_si(resElems[k].coef, tmp);
	resElems[k].degs = 0ull;
	resElems[k].degs |= ( ((degrees_t) curX) << EXP_OFFSET_1_V2);
	resElems[k].degs |= ( ((degrees_t) curX) << EXP_OFFSET_2_V2);
	++k;

	ret->size = k;
} 

//point powers has powers increasing along columns, hence the basic point is in first column.
//eval the poly as if it is being evaluated at (basic point)^j
//ret should be initialized to a polynomial of size of the bivariate support of a 
//with all mpz coefs initialized
void evalToBivarSymmetricSPF_AAZ(AltArrZ_t* ret, const AltArrZ_t* a, elem_t const*const* pointPowers, int j, const degrees_t* masks, const int* offsets, degrees_t biMask, const Prime_ptr* Pptr) {
	if (ret == NULL || isZero_AAZ(a)) {
		return;
	}

	if (a->unpacked) {
		evalToBivarSymmetricSPF_AAZ_unpk(ret, a, pointPowers, j, Pptr);
		return;
	}

	//Allocate a little bit extra in place of iterating though a to get true bivaraite support ammount.
	AAZElem_t* resElems = ret->elems;
	int size = a->size;
	int nvar = a->nvar;
	int k = 0;
	int v;
	degrees_t degs;
	degrees_t curDeg = a->elems->degs & biMask;
	elem_t curCoef = smallprimefield_convert_in(0, Pptr);
	elem_t tmp;
	unsigned long lPrime = Pptr->prime;
	unsigned long halfP = lPrime / 2;
	for (int i = 0; i < size; ++i) {
		degs = a->elems[i].degs;
		if (curDeg != (degs & biMask)) {
			//commit previous result to ret;
			tmp = smallprimefield_convert_out(curCoef, Pptr);
			if (tmp > halfP) {
				tmp -= lPrime;
			}
			mpz_set_si(resElems[k].coef, tmp);
			resElems[k].degs = 0ull;
			resElems[k].degs |= (GET_NTH_EXP(curDeg, masks[0], offsets[0]) << EXP_OFFSET_1_V2);
			resElems[k].degs |= (GET_NTH_EXP(curDeg, masks[1], offsets[1]) << EXP_OFFSET_2_V2);
			++k;
		
			curCoef = smallprimefield_convert_in(0, Pptr);
			curDeg = degs & biMask;
		}

		tmp = smallprimefield_convert_in(1l, Pptr);
		for (v = 2; v < nvar; ++v) {
			tmp = smallprimefield_mul(tmp, pointPowers[v-2][j*GET_NTH_EXP(degs, masks[v], offsets[v])], Pptr);
		}
		tmp = smallprimefield_mul(tmp, smallprimefield_convert_in(mpz_fdiv_ui(a->elems[i].coef, lPrime), Pptr), Pptr);
		curCoef = smallprimefield_add(curCoef, tmp, Pptr);
	}


	tmp = smallprimefield_convert_out(curCoef, Pptr);
	if (tmp > halfP) {
		tmp -= lPrime;
	}
	//one more to commit the final term since it wasn't triggered in the loop
	mpz_set_si(resElems[k].coef, tmp);
	resElems[k].degs = 0ull;
	resElems[k].degs |= (GET_NTH_EXP(degs, masks[0], offsets[0]) << EXP_OFFSET_1_V2);
	resElems[k].degs |= (GET_NTH_EXP(degs, masks[1], offsets[1]) << EXP_OFFSET_2_V2);
	++k;

	ret->size = k;
} 


int checkSigmasSupports(AltArrZ_t** sigmas, int nSigmas, int nImages) {

	//NOTE we assume sigmas are packed here since they are only ever bivariate

	//j increments along columns, i.e. unique sigmas
	//i increments along sets of images
	for (int j = 0; j < nSigmas; ++j) {
		int size = sigmas[j] == NULL ? 0 : sigmas[j]->size;
		if (size == 0) {
			for (int i = 1; i < nImages; ++i) {
				if (!isZero_AAZ(sigmas[nSigmas*i + j])) {
					return 0;
				}
			}	
		} else {
			for (int i = 1; i < nImages; ++i) {
				if (size != sigmas[nSigmas*i + j]->size) {
					return 0;
				}
			}
			degrees_t degs;
			for (int k = 0; k < size; ++k) {
				degs = sigmas[j]->elems[k].degs;
				for (int i = 1; i < nImages; ++i) {
					if (degs != sigmas[nSigmas*i + j]->elems[k].degs) {
						return 0;
					}
				}
			}
		}
	}

	return 1;
}

//sigmas should be pre-generated with a skeleton of the sigma. 
//c->nvar >= 2; Moreover, should have positive degree in main variable
int multiMDP_AAZ(const AltArrZ_t* cc, AltArrZ_t const* const* us, AltArrZ_t** sigmas, unsigned int r, const Prime_ptr* Pptr, const mpz_t mpPrime) {
// static inline void rand_mpz_vec(unsigned long int coefBound, int includeNeg, mpz_t* mpzVec, unsigned int n) {

	if (cc->nvar == 2) {
		return multiBDP_AAZ(cc, us, sigmas, r, Pptr, mpPrime);
	}

	if (cc->nvar < 3) {
		return 0;
	}

	AltArrZ_t* c = deepCopyPolynomial_AAZ(cc);

	short nvar = c->nvar;
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* offsets  = getExpOffsetArray(nvar);
	degrees_t biMask = masks[0] | masks[1];

	unsigned int i, j;

	RecBivarAAZ_t** recSigmas = (RecBivarAAZ_t**) malloc(sizeof(RecBivarAAZ_t*)*r);
	unsigned int tMax = 0u;
	for (i = 0u; i < r; ++i) {
		recSigmas[i] = convertToRecBivar_AAZ(sigmas[i], biMask, masks, offsets);
		
		for (j = 0u; recSigmas[i] != NULL &&  j < recSigmas[i]->size; ++j) {
			tMax = tMax < recSigmas[i]->elems[j].coefSize ? recSigmas[i]->elems[j].coefSize : tMax;
		}
	}

	if (tMax == 0) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "tMax was 0 in multiMDP... should you be calling multiBDP?\n");
#endif
		return 0;
	}

	unsigned int evalSetSize = nvar - 2;

	//compute evaluation points for images as powers of a random vector
	//compute monomial evaluations for later Vandermonde matrices
	//must do this first to ensure out evalPoint is good, i.e., 
	//all monomial evaluations come out to be unique.
	elem_t* evalPoint = malloc(sizeof(elem_t)*evalSetSize);
	elem_t** monomialEvalSets = (elem_t**) malloc(sizeof(elem_t*)*r);
	for (i = 0u; i < r; ++i) {
		if (sigmas[i] == NULL) {
			monomialEvalSets[i] = NULL;
		} else {
			monomialEvalSets[i] = (elem_t*) malloc(sizeof(elem_t)*sigmas[i]->size);
		}
	}

	int success = 0;
	int maxTries = 3;
	for (int ntries = 0; ntries < maxTries && !success; ++ntries) {
		for (i = 0u; i < evalSetSize; ++i) {
			elem_t randVal = randValInRange(1, Pptr->prime-1);
			evalPoint[i] = smallprimefield_convert_in(randVal, Pptr);	
		}

		success = 1;
		for (i = 0u; i < r && success; ++i) {
			if (recSigmas[i] == NULL) {
				continue;
			}

			//here, if you sum the size of each recursive coefficient, it's 
			//just the total size of the expanded sigma.
			elem_t* evalSet = monomialEvalSets[i];

			int curIdx = 0;
			//so now evaluate each monomial of each coefficient polynomial.
			for (j = 0u; j < recSigmas[i]->size; ++j) {
				evalRecMonomials(recSigmas[i], j, evalPoint, evalSet + curIdx, nvar, masks, offsets, Pptr);
				if (!uniqueElemTSet(evalSet+curIdx, recSigmas[i]->elems[j].coefSize)) {
					success = 0;
					break;
				}

				curIdx += recSigmas[i]->elems[j].coefSize;
			}
		}

	}
	if (!success) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "multiMDP failed; could not find unique set of evaluations\n" );
#endif
		for (i = 0u; i < r; ++i) {
			free(monomialEvalSets[i]);
			freeRecBivar_AAZ(recSigmas[i]);
		}
		free(monomialEvalSets);
		free(evalPoint);
		return 0;
	}

	//the set is good, so let's compute powers of this random point 
	//to later be used as a lookup table for evaluations 

	//get max partial degrees across all polys to evaluate;
	degree_t maxDegs[nvar];
	degree_t degs[nvar];
	partialDegrees_AAZ(c, maxDegs);
	for (i = 0u; i < r; ++i) {
		partialDegrees_AAZ(us[i], degs);
		for (j = 2u; j < nvar; ++j) {
			// fprintf(stderr, "degs[%d]: %d, maxDges[%d], %d\n", j, degs[j], j, maxDegs[j]);
			maxDegs[j] = maxDegs[j] < degs[j] ? degs[j] : maxDegs[j];
		}
	}

	//build a table of evalPoint powers for evaluation.
	elem_t** pointPowers = (elem_t**) malloc(sizeof(elem_t*)*evalSetSize);
	degree_t powerSize;
	degree_t deg;
	for (i = 0u; i < evalSetSize; ++i) {
		// +2 since eval point is in x_3...x_nvar
		deg = maxDegs[i+2];
		deg = deg <= 0 ? 1 : deg;
		powerSize = deg <= 1 ? 2 : deg; //make deg at least 2 so that malloc >= 2
		powerSize *= tMax;
		pointPowers[i] = (elem_t*) malloc(sizeof(elem_t)*powerSize+1); //(point^tMax)^deg, +1 for indexing
		pointPowers[i][0] = smallprimefield_convert_in(1ll, Pptr);
		pointPowers[i][1] = evalPoint[i];
		for (int k = 2; k <= deg*tMax; ++k) {
			pointPowers[i][k] = smallprimefield_mul(pointPowers[i][k-1], evalPoint[i], Pptr);
		}
	}

	//setup polynomials to store the resulting evaluations images which are inputs to multiBDP.
	AltArrZ_t** uImages = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*r);
	unsigned int bivarSuppSize;
	countBivarSupportCoefSize_AAZ(c, biMask, &bivarSuppSize, NULL);	
	AltArrZ_t* cImage = makePolynomial_AAZ(bivarSuppSize,2); //nvar=2
	for (i = 0u; i < bivarSuppSize; ++i) {
		mpz_init(cImage->elems[i].coef);
	}
	cImage->size = bivarSuppSize;

	for (i = 0u; i < r; ++i) {
		countBivarSupportCoefSize_AAZ(us[i], biMask, &bivarSuppSize, NULL);	
		uImages[i] = makePolynomial_AAZ(bivarSuppSize,2);
		for (j = 0u; j < bivarSuppSize; ++j) {
			mpz_init(uImages[i]->elems[j].coef);
		}
		uImages[i]->size = bivarSuppSize;
	}
	
	//each row of sigmasImages is one multiBDP solution, we need tMax solutions.
	AltArrZ_t** sigmasImages = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*r*tMax);
	//solve tMax BDPs
	int ntries;
	for (i = 0u; i < tMax; ++i) {
		
		evalToBivarSymmetricSPF_AAZ(cImage, c, CONSTCONSTCAST(elem_t, pointPowers), i+1, masks, offsets, biMask, Pptr);
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "c: \n");
		printPoly_AAZ(stderr, c, henselsyms, c->nvar);
		fprintf(stderr, "\nat z = %lld\n,cImage: \n", smallprimefield_convert_out(pointPowers[0][1], Pptr));
		printPoly_AAZ(stderr, cImage, henselsyms, 2);
#endif
		for (j = 0u; j < r; ++j) {
			evalToBivarSymmetricSPF_AAZ(uImages[j], us[j], CONSTCONSTCAST(elem_t, pointPowers), i+1, masks, offsets, biMask, Pptr);
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr, "\nu[%d]: \n", j);
			printPoly_AAZ(stderr, us[j], henselsyms, c->nvar);
			fprintf(stderr, "\nat z = %lld\n,uimage[%d]: \n", smallprimefield_convert_out(pointPowers[0][1], Pptr), j);
			printPoly_AAZ(stderr, uImages[j], henselsyms, 2);
#endif
		}

		success = 0;
		for (ntries = 0; ntries < 5 && !success; ++ntries) {
			success = multiBDP_AAZ(cImage, CONSTCONSTCAST(AltArrZ_t, uImages), sigmasImages + (i*r), r, Pptr, mpPrime);
			// fprintf(stderr, "\nBDP ntry %d\n", ntries+1);
		}

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "\nBDP success for i = %d: %d\n", i, success);
#endif
		
		if (!success) {
			for (int  k = 0; k < i; ++k) {
				for (int l = 0; l < r; ++l) {
					freePolynomial_AAZ(sigmasImages[k*r + l]);
				}
			}
			break;
		}
	}

	freePolynomial_AAZ(cImage);
	for (i = 0u; i < r; ++i) {
		freePolynomial_AAZ(uImages[i]);
	}
	free(uImages);

	// fprintf(stderr, "finished doing bdp\n" );

	//determine if supports are the same for matchin sigmas
	if (success) {
		success = checkSigmasSupports(sigmasImages, r, tMax); 
	}

	// fprintf(stderr, "sigmsa supports sucess: %d\n", success );

	//we can sparsely interpolate sigmas now 
	//foreach sigma[i], i = 1..r
	//    foreach bivar coef, j = i..w, w <= tMax
	//        setup and solve the Vandermonde system
	if (success) {
		// foreach bivar coef
		int curIdx = 0;
		int curImgIdx = 0;
		int curT = 0;
		int transposed = 1; //when solvng system, exponents increase down a column instead of across a row
		int startExp = 1;
		RecBivarAAZ_t* rec;
		degrees_t bivarDegs;
		elem_t* vandRes;
		elem_t halfP = Pptr->prime >> 1;
		elem_t tmp;
		if (evalSetSize > 1) {
			//reuse pointPowers to store b, each are guaranteed to have at least tMax 
			vandRes = pointPowers[1];
		} else {
			vandRes = (elem_t*) malloc(sizeof(elem_t)*tMax);
		}

		for (i = 0u; i < r; ++i) {
			// fprintf(stderr, "sparse interp for sigma[%d]\n", i);

			if (recSigmas[i] == NULL) {
				continue;
			}

			//if so, all sigmasImages for i are zero by previous checkSigmasSupports
			//hence, sigmas must be identically 0.
			if (isZero_AAZ(sigmasImages[i])) {
				freePolynomial_AAZ(sigmas[i]);
				sigmas[i] = NULL;
				continue;
			}


			rec = recSigmas[i];
			curIdx = 0;

			//rec->elems[j] and sigmasImages[...]->elems[j] won't align if sigmasImages has some which evalute to 0.
			//use curImgIdx for proper index sigmasImages->elems. It is incremented when we find non-0 coef which
			//matches rec->elems[j].
			curImgIdx = 0;


			for (j = 0u; j < rec->size; ++j) {
				bivarDegs = rec->elems[j].degs12;					

				curT = rec->elems[j].coefSize;
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr, "curT: %d\n", curT);
#endif

				//Check if bivarDegs is the same as sigmasImages[...]->elems[j].degs. 
				//Otherwise, coef of bivarDegs in true sigmas is 0.
				//We only need to check one sigma image since all sigmas have same support by checkSigmasSupports above.
				if (bivarDegs != sigmasImages[i]->elems[curImgIdx].degs) {
					// fprintf(stderr, "skipping degs %llx \n", bivarDegs);
					// fprintf(stderr, "sigmas images degs %llx \n", sigmasImages[i]->elems[curImgIdx].degs);
					for (int k = 0; k < curT; ++k) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
						if (curIdx + k >= sigmas[i]->size) {
							fprintf(stderr, "\n\nWOOOAH NOW\n");
							free((int*)-1);
						}
#endif
						// fprintf(stderr, "setting sigmas[%d]->elems[%d] = %ld\n", i, curIdx+k, tmp);
						mpz_set_ui(sigmas[i]->elems[curIdx+k].coef, 0ul);
					}
					curIdx += curT;
					continue;
				}

				for (int k = 0; k < curT; ++k) {
					pointPowers[0][k] = smallprimefield_convert_in(mpz_fdiv_ui(sigmasImages[(k*r) + i]->elems[curImgIdx].coef, Pptr->prime), Pptr);
				}

				//setup the Vandermonde system
				success = solveFFVandermondeSystemInForm(monomialEvalSets[i] + curIdx, pointPowers[0], curT, startExp, transposed, &vandRes, Pptr);

				if (!success) {
					// Well this shouldn't happend since we checked that monomialEvalSets was unique above
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
					fprintf(stderr, "Solving Vandermonde system in multiMDP_AAZ failed!?\n");
#endif
					break;
				}

				for (int k = 0; k < curT; ++k) {
					tmp = smallprimefield_convert_out(vandRes[k], Pptr);
					//symmetric convert out!
					if (tmp > halfP) {
						tmp -= Pptr->prime;
					}

					//update the true sigmas now!
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
					if (curIdx + k >= sigmas[i]->size) {
						fprintf(stderr, "\n\nWOOOAH NOW\n");
						free((int*)-1);
					}
#endif
					mpz_set_si(sigmas[i]->elems[curIdx+k].coef, tmp);
				}

				// fprintf(stderr, "finished updating the sigma\n\n" );

				curIdx += curT;
				++curImgIdx;
			} //j = 0..curT

			canonicalizePolynomial_AAZ(sigmas[i]);

			if (!success) {
				break;
			}
		} // i = 0..r	

		if (evalSetSize <= 1) {
			free(vandRes);
		}

		// TODO
		//probabilistic verification

		//now check if correct;
		AltArrZ_t* tmpProd = makePolynomial_AAZ(c->size, c->nvar);
		AltArrZ_t* sum = NULL;
		AltArrZ_t* b;
		for (int i = 0; i < r; ++i) {
			b = multiplyAllButOnePolynomials_AAZ(us, r, i);
			// fprintf(stderr, "b: \n");
			// printPoly_AAZ(stderr, b, henselsyms, b->nvar);
			// fprintf(stderr, "\nsigmas[%d]:\n", i);
			// printPoly_AAZ(stderr, sigmas[i], henselsyms, sigmas[i] == NULL ? 0 :sigmas[i]->nvar);
			multiplyPolynomialsPreAlloc_AAZ(sigmas[i], b, &tmpProd);
			sum = addPolynomials_AAZ_inp(sum, tmpProd, tmpProd->nvar);
			freePolynomial_AAZ(b);
		}

		applyModuloSymmetric_AAZ_inp(sum, mpPrime);
		// fprintf(stderr, "\nsum mod:\n" );
		// printPoly_AAZ(stderr, sum, henselsyms, sum->nvar);
		// fprintf(stderr, "\n\n");
		// fprintf(stderr, "\nc:\n" );
		// printPoly_AAZ(stderr, c, henselsyms, c->nvar);
		// fprintf(stderr, "\n\n");
		success = isExactlyEqual_AAZ(c, sum);
		freePolynomial_AAZ(sum);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		if (success) {
			fprintf(stderr, "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n    Multi MDP Success!\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
		} else {
			fprintf(stderr, "\n\nMDP Failed!\n\n" );
		}
#endif
		for (j = 0u; j < tMax; ++j) {
			for (unsigned int k = 0; k < r; ++k) {
				freePolynomial_AAZ(sigmasImages[(j*r) + k]);				
			}
		}
	}

	//individual polys either freed as BDP fail or just above at end of
	free(sigmasImages);
	for (i = 0u; i < evalSetSize; ++i) {
		free(pointPowers[i]);
	}
	free(pointPowers);
	for (i = 0u; i < r; ++i) {
		freeRecBivar_AAZ(recSigmas[i]);	
	}
	free(recSigmas);

	for (i = 0u; i < r; ++i) {
		free(monomialEvalSets[i]);
	}
	free(monomialEvalSets);

	return 1;
}


int multivarHenselLiftVars_AAZ(const AltArrZ_t* aa, AltArrZ_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t* x, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	if (aa == NULL || lcF == NULL || fs == NULL || liftedF == NULL || r < 2 || x == NULL) {
		return 0;
	}

	mpz_t mpPrime;
	mpz_init_set_si(mpPrime, Pptr->prime);

	short nvar = aa->nvar;
	int active[nvar];
	//evaluate a and LC's to bivar first
	active[0] = 0;
	active[1] = 0;
	for (int i = 2; i < nvar; ++i) {
		active[i] = 1;
	}

	AltArrZ_t* aMod = deepCopyPolynomial_AAZ(aa);
	applyModuloSymmetric_AAZ_inp(aMod, mpPrime);
	
	AltArrZ_t* a = evaluatePoly_AAZ(aMod, active, x, nvar);
	applyModuloSymmetric_AAZ_inp(a, mpPrime);
	
	duspoly_t** fps = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	duspoly_t** sigmas = (duspoly_t**) malloc(sizeof(duspoly_t*) * r);
	AltArrZ_t** tmpLifted = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*) * r);
	for (unsigned int i = 0; i < r; ++i) {
		fps[i] = convertToDUSP_DUZP(fs[i], Pptr);
		AltArrZ_t* tmpLC = evaluatePoly_AAZ(lcF[i], active, x, nvar); 
		tmpLifted[i] = replaceLC_DUZP(fs[i], tmpLC, mpPrime);
		freePolynomial_AAZ(tmpLC);
	}

	AltArrZ_t* error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, tmpLifted), r);
	subPolynomials_AAZ_inpRHS(a, &error);
	applyModuloSymmetric_AAZ_inp(error, mpPrime);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
	fprintf(stderr, "\ngot first error: ");
	printPoly_AAZ(stderr, error, henselsyms, error->nvar);
	fprintf(stderr, "\n\n");
#endif


	//here we manually do the first lift since the factors we are lifting from can be encoded
	//as DUSP. For every other variable it's basically the same

	degree_t deg = partialDegree_AAZ(a, 1);
	AltArrZ_t* kDeriv;
	AltArrZ_t* c;
	duspoly_t* cp = makePolynomial_spX(a->size + 1);

	mpz_t kfact;
	mpz_init_set_ui(kfact, 1ul);
	active[1] = 1; //evaluate error term to univar

	for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {
		kDeriv = derivative_AAZ(error, 1, k);

		c = evaluatePoly_AAZ(kDeriv, active, x, 2);
		freePolynomial_AAZ(kDeriv);
		divideByIntegerExact_AAZ_inp(c, kfact); //divide by k!		
		univarToPrimeFieldPreAlloc_AAZ(c, Pptr, &cp);
		freePolynomial_AAZ(c);

		multiUDP_spX(cp, CONSTCONSTCAST(duspoly_t, fps), r, sigmas, Pptr);

		for (unsigned int i = 0; i < r; ++i) {
			idealadicUpdate(tmpLifted + i, sigmas[i], x[1], k, Pptr, mpPrime);
			freePolynomial_spX (&sigmas[i]);
		}

		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, tmpLifted), r, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);

		mpz_mul_ui(kfact, kfact, k+1ul);
	}	

	//cleanup and check success before next round of lifting.
	for (unsigned int i = 0; i < r; ++i) {
		freePolynomial_spX (&fps[i]);
	}
	free(fps);	
	
	freePolynomial_AAZ(a); //will re-use variable a for subsequent loops;
	a = NULL;
	freePolynomial_spX(&cp);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
	fprintf(stderr, "Checking error after first round of lift\n");
	printPoly_AAZ(stderr, error, henselsyms, error->nvar);
	fprintf(stderr, "\n\n" );
#endif

	if (!isZero_AAZ(error)) {
		//lift from univar to bivar failed!
		freePolynomial_AAZ(error);
		freePolynomial_AAZ(aMod);
		mpz_clear(kfact);
		free(sigmas);
		return 0;
	}

	//reset everything now and start general lifting loop

	AltArrZ_t** sigmas_aaz;
	if (sizeof(AltArrZ_t*) == sizeof(duspoly_t*)) {
		sigmas_aaz = (AltArrZ_t**) sigmas;
	} else {
		free(sigmas);
		sigmas_aaz = (AltArrZ_t**) malloc(sizeof(AltArrZ_t)*r);
	}
	sigmas = NULL;
	for (int i = 0; i < r; ++i) {
		sigmas_aaz[i] = NULL;
	}

	for (int v = 2; v < nvar; ++v) {
		active[v-1] = 0; //reset last loop's active for error eval
		active[v] = 0;

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "\n\n=================================================================\nstarting lift v=%d\n=================================================================\n\n\n", v);
#endif
		if (v < nvar - 1) {
			a = evaluatePoly_AAZ(aMod, active, x, nvar);
		} else {
			a = aMod;
			aMod = NULL;
		}
		deg = partialDegree_AAZ(a, v);
		if (deg == 0) {
			freePolynomial_AAZ(a);
			a = NULL;
			continue;
		}
		applyModuloSymmetric_AAZ_inp(a, mpPrime);


		if (v > 2) {
			for (unsigned int i = 0; i < r; ++i) {
				sigmas_aaz[i] = deepCopyPolynomial_AAZ(tmpLifted[i]);
			}
		}

		//replaceLC also expands nvar by 1.
		if (v < nvar-1) {
			for (unsigned int i = 0; i < r; ++i) {
				// fprintf(stderr, "before replace %d\n",i );
				// printPoly_AAZ(stderr, tmpLifted[i], henselsyms, tmpLifted[i]->nvar);
				// fprintf(stderr, "\n" );
				// fprintf(stderr, "lc %d\n",i );
				// printPoly_AAZ(stderr, lcF[i], henselsyms, lcF[i]->nvar);
				// fprintf(stderr, "\n" );

				AltArrZ_t* tmpLC = evaluatePoly_AAZ(lcF[i], active, x, nvar); 

				// fprintf(stderr, "tmplc %d\n",i );
				// printPoly_AAZ(stderr, tmpLC, henselsyms, tmpLC->nvar);
				// fprintf(stderr, "\n" );
				liftedF[i] = replaceLC_AAZ(tmpLifted[i], tmpLC);
				freePolynomial_AAZ(tmpLC);
				applyModuloSymmetric_AAZ_inp(liftedF[i], mpPrime);	
			
				// fprintf(stderr, "after replace %d\n",i );
				// printPoly_AAZ(stderr, liftedF[i], henselsyms, liftedF[i]->nvar);
				// fprintf(stderr, "\n" );
			}		
		} else {
			for (unsigned int i = 0; i < r; ++i) {
				liftedF[i] = replaceLC_AAZ(tmpLifted[i], lcF[i]);
				applyModuloSymmetric_AAZ_inp(liftedF[i], mpPrime);
			} 
		}

		mpz_set_ui(kfact, 1ul);
		active[v] = 1; //evaluate coefs of error to one less variables.

		error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r);

 		subPolynomials_AAZ_inpRHS(a, &error);	 
		applyModuloSymmetric_AAZ_inp(error, mpPrime);
		int dioSuccess = 0;
		for (unsigned long k = 1; k <= deg && !isZero_AAZ(error); ++k) {

			//TODO use improved version of getting taylor coef, avoid computing full derivative.
			//see Monogan and Pearce: POLY: New dast struct in Maple 17.
			//derive and evaluate w.r.t to third var
			kDeriv = derivative_AAZ(error, v, k);
			c = evaluatePoly_AAZ(kDeriv, active, x, kDeriv->nvar);
			divideByIntegerExact_AAZ_inp(c, kfact);
			applyModuloSymmetric_AAZ_inp(c, mpPrime);
			freePolynomial_AAZ(kDeriv);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr, "c for k = %d \n", k);
			printPoly_AAZ(stderr, c, henselsyms, c->nvar);
			fprintf(stderr, "\n" );
			for (int l = 0; l < r; ++l){
				fprintf(stderr, "tmpLifted[%d]: \n", l);
				printPoly_AAZ(stderr, tmpLifted[l], henselsyms, tmpLifted[l]->nvar);
				fprintf(stderr, "\n");
				fprintf(stderr, "sigmas_aaz[%d]: \n", l);
				printPoly_AAZ(stderr, sigmas_aaz[l], henselsyms, sigmas_aaz[l] == NULL ? 0 : sigmas_aaz[l]->nvar);
				fprintf(stderr, "\n");
			}
#endif

			dioSuccess = multiMDP_AAZ(c, CONSTCONSTCAST(AltArrZ_t, tmpLifted), sigmas_aaz, r, Pptr, mpPrime);
			freePolynomial_AAZ(c);

			if (!dioSuccess) {
				break;
			}

			for (unsigned int i = 0; i < r; ++i) {
				multivarIdealadicUpdate(liftedF + i, sigmas_aaz[i], x[v], k, mpPrime);
				if (v == 2) {
					//for v > 2 sigmas are used as skeleton in multiMDP
					freePolynomial_AAZ(sigmas_aaz[i]);
					sigmas_aaz[i] = NULL;
				}
			}

			multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, liftedF), r, &error);
			subPolynomials_AAZ_inpRHS(a, &error);
			applyModuloSymmetric_AAZ_inp(error, mpPrime);

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr, "\n\nupdated error: \n" );
			printPoly_AAZ(stderr, error, henselsyms, error == NULL ? 0 : error->nvar);
			fprintf(stderr, "\n\n" );
#endif
			mpz_mul_ui(kfact, kfact, k+1ul);
		}

		//after each variable lift, make newly liftedF tmpLifted for next round
		if (v < nvar-1) {
			for (int i = 0; i < r; ++i) {
				freePolynomial_AAZ(tmpLifted[i]);
				tmpLifted[i] = liftedF[i];
				liftedF[i] = NULL;
			}
		} else {
			for (unsigned int i = 0; i < r; ++i) {
				freePolynomial_AAZ(tmpLifted[i]);
			}
		}

		for (unsigned int i = 0; i < r; ++i) {
			freePolynomial_AAZ(sigmas_aaz[i]);
		}

		freePolynomial_AAZ(a); //this takes care of freeing aMod also.
		a = NULL;
	}

	free(tmpLifted);
	free(sigmas_aaz);

	mpz_clear(kfact);

	int ret = isZero_AAZ(error);
	freePolynomial_AAZ(error);
	return ret;

}



//liftedF already have true LC multipled into them.
int henselLiftCoefs_AAZ(const AltArrZ_t* a, AltArrZ_t** Fs, unsigned int r, const mpz_t bound, const Prime_ptr* Pptr) {

	if (Fs == NULL) {
		return 0;
	}

	AltArrZ_t* error = multiplyManyPolynomials_AAZ( CONSTCONSTCAST(AltArrZ_t, Fs), r);
	subPolynomials_AAZ_inpRHS(a, &error);	 

	if (isZero_AAZ(error)) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "No lifting to do in hensel lift coefs\n");
#endif
		freePolynomial_AAZ(error);
		return 1;
	}

	//alloc for bezout coefficients
	AltArrZ_t** sigmas = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*r);
	for (int i = 0; i < r; ++i) {
		sigmas[i] = deepCopyPolynomial_AAZ(Fs[i]);
	}

	mpz_t mod; 
	mpz_init_set_ui(mod, Pptr->prime);
	mpz_t mpPrime;
	mpz_init_set(mpPrime, mod);
	int primePow = 1;

	int success = 0;
	while( !isZero_AAZ(error) && mpz_cmp(mod, bound) < 0 ) {
		//c = error / p^k
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr, "error before: \n");
		printPoly_AAZ(stderr, error, henselsyms, error->nvar);
#endif
		divideByIntegerExact_AAZ_inp(error, mod);
		applyModuloSymmetric_AAZ_inp(error, mpPrime);
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG	
		fprintf(stderr, "\nerror after: \n");
		printPoly_AAZ(stderr, error, henselsyms, error->nvar);
		fprintf(stderr, "\n");
#endif
		
		success = multiMDP_AAZ(error, CONSTCONSTCAST(AltArrZ_t, Fs), sigmas, r, Pptr, mpPrime);

		if (!success) {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr, "in henselLiftCoefs: multi MDP FAILED!!! for primePow=%d\n", primePow);
#endif
			break;
		}

		for (int i = 0; i < r; ++i) {
			multiplyByInteger_AAZ_inp(sigmas[i], mod);
			// fprintf(stderr, "\nbefore Fs[%d]: \n", i);
			// printPoly_AAZ(stderr, Fs[i], henselsyms, Fs[i]->nvar);
			Fs[i] = addPolynomials_AAZ_inp(Fs[i], sigmas[i], Fs[i]->nvar);
			// fprintf(stderr, "\nafter Fs[%d]: \n", i);
			// printPoly_AAZ(stderr, Fs[i], henselsyms, Fs[i]->nvar);
		}
		// fprintf(stderr, "\n" );

		//update the error.
		multiplyManyPolynomialsPreAlloc_AAZ( CONSTCONSTCAST(AltArrZ_t, Fs), r, &error);
		subPolynomials_AAZ_inpRHS(a, &error);
		mpz_mul_si(mod, mod, Pptr->prime);
		++primePow;
	}

	if (!isZero_AAZ(error)) {
		//we got here because bound exceeded.
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		gmp_fprintf(stderr, "Failure or Bound exceeded in henselLiftCoefs_AAZ, current mod: %Zd\n", mod);
#endif
	}

	for (int i = 0; i < r; ++i) {
		freePolynomial_AAZ(sigmas[i]);
	}
	free(sigmas);	

	freePolynomial_AAZ(error);
	mpz_clears(mod, mpPrime, NULL);

	return primePow;
}

int trivarHenselLift_AAZ(const AltArrZ_t* a, AltArrZ_t const*const* lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t* x, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	int ret = trivarHenselLiftVars_AAZ(a, lcF, fs, r, x, Pptr, liftedF);

	if (ret) {
		mpz_t bound;
		mpz_init(bound);
		infinityNorm_AAZ(a, bound);
		mpz_mul(bound, bound, a->elems->coef);
		degree_t pdegs[a->nvar];
		double degProd = 1;
		partialDegrees_AAZ(a, pdegs);
		for (int i = 0; i < a->nvar; ++i) {
			degProd *= (pdegs[i] + 1.0);
		}
		degProd = sqrt(degProd);
		mpz_mul_si(bound, bound, (long) (degProd + 0.5));
		degree_t deg = totalDegree_AAZ(a);
		mpz_mul_2exp(bound, bound, deg);

		//if bound is already smaller than prime, we must have captured true coefs
		//with just variable lifting over Z_p.
		if (mpz_cmp_ui(bound, Pptr->prime) < 0) {
			mpz_clear(bound);
			return ret;
		}

		AltArrZ_t* lc;
		for (int i = 0; i < r; ++i) {
			lc = deepCopyPolynomial_AAZ(lcF[i]);			
			replaceLC_AAZ_inp(liftedF + i, lc);
			freePolynomial_AAZ(lc);
		}

		ret = henselLiftCoefs_AAZ(a, liftedF, r, bound, Pptr);

		mpz_t tmp;
		mpz_init_set_ui(tmp, Pptr->prime);
		mpz_pow_ui(tmp, tmp, ret);
		ret = mpz_cmp(tmp, bound) < 0;

		mpz_clear(tmp);
		mpz_clear(bound);
	}


	return ret;

}


int henselLift_AAZ(const AltArrZ_t* aa, AltArrZ_t** lcF, DUZP_t const*const* fs, unsigned int r, const mpz_t* x_in, const Prime_ptr* Pptr, AltArrZ_t** liftedF) {

	int ringNvar = aa->nvar;
	int varMap[ringNvar];

	AltArrZ_t* a = deepCopyPolynomial_AAZ(aa);
	int trueNvar = tryShrinkVariables_AAZ_inp(a, varMap);

	if (trueNvar == 1) {
		freePolynomial_AAZ(a);
		return 0;
	}

	mpz_t x[trueNvar];
	if (trueNvar != ringNvar) {
		for (int i = 0; i < r; ++i) {
			shrinkAndReorderVars_AAZ((AltArrZ_t*) lcF[i], varMap, ringNvar);
		}

		int j = 0;
		for (int i = 0; i < ringNvar; ++i) {
			if (varMap[i] >= 0) {
				mpz_init_set(x[j], x_in[i]);
				++j;
			}
		}

	} else {
		a = (AltArrZ_t*) aa;
		memcpy(x, x_in, sizeof(mpz_t)*trueNvar);
	}

	int nvar = a->nvar;
	int ret;
	if (nvar == 2) {
		DUZP_t** dzLC = (DUZP_t**) malloc(sizeof(DUZP_t*)*r);
		for (int i = 0; i < r; ++i) {
			if (partialDegree_AAZ(lcF[i], 0) != 0) {
				fprintf(stderr, "Hensel lifting failed, inconsistent product and true leading coefficients.\n" );
				return 0;
			}
			shrinkNumVarsAtIdx_AAZ(lcF[i], 0);
			dzLC[i] = convertFromAltArrZ_DUZP(lcF[i]);
		}
		
		ret =  bivarHenselLift_AAZ(a, CONSTCONSTCAST(DUZP_t,dzLC), fs, r, x[1], Pptr, liftedF);


		for (int i = 0; i < r; ++i) {
			freePolynomial_DUZP(dzLC[i]);
		}
		free(dzLC);

	} else if (nvar == 3) {
		ret = trivarHenselLift_AAZ(a, CONSTCONSTCAST(AltArrZ_t, lcF), fs, r, x, Pptr, liftedF);
	} else {

		ret = multivarHenselLiftVars_AAZ(a, CONSTCONSTCAST(AltArrZ_t, lcF), fs, r, x, Pptr, liftedF);

		if (ret) {
			mpz_t bound;
			mpz_init(bound);
			infinityNorm_AAZ(a, bound);
			mpz_mul(bound, bound, a->elems->coef);
			degree_t pdegs[a->nvar];
			double degProd = 1;
			partialDegrees_AAZ(a, pdegs);
			for (int i = 0; i < a->nvar; ++i) {
				degProd *= (pdegs[i] + 1.0);
				fprintf(stderr, "pdegs[%d]: %d\n",i, pdegs[i]);
			}
			degProd = sqrt(degProd);
			fprintf(stderr, "sqrt degProd: %f, ceill: %ld\n", degProd, (long) (degProd + 0.5));
			mpz_mul_si(bound, bound, (long) (degProd + 0.5));
			degree_t deg = totalDegree_AAZ(a);
			mpz_mul_2exp(bound, bound, deg);

			//if bound is already smaller than prime, we must have captured true coefs
			//with just variable lifting over Z_p.
			if (mpz_cmp_ui(bound, Pptr->prime) < 0) {
				mpz_clear(bound);
			} else {
				//replace with proper multi-precision LC's
				AltArrZ_t* lc;
				for (int i = 0; i < r; ++i) {
					lc = deepCopyPolynomial_AAZ(lcF[i]);			
					expandNumVarsLeft_AAZ(lc, a->nvar);
					replaceLC_AAZ_inp(liftedF + i, lc);
					freePolynomial_AAZ(lc);
				}

				ret = henselLiftCoefs_AAZ(a, liftedF, r, bound, Pptr);

				mpz_t tmp;
				mpz_init_set_ui(tmp, Pptr->prime);
				mpz_pow_ui(tmp, tmp, ret);
				ret = mpz_cmp(tmp, bound) < 0;

				mpz_clear(tmp);
				mpz_clear(bound);
			}
		}

	}

	if (trueNvar != ringNvar) {
		freePolynomial_AAZ((AltArrZ_t*) a);
		for (int i = 0; i < trueNvar; ++i) {
			mpz_clear(x[i]);
		}	

		int revMap[ringNvar];
		reverseShrinkMap_AAZ(ringNvar, trueNvar, varMap, revMap);

		for (int i = 0; i < r; ++i) {
			expandNumVars_AAZ(liftedF[i], ringNvar);
			reorderVars_AAZ(liftedF[i], revMap, ringNvar);
		}
	}

	return ret;
}

