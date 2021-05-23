

#include "IntegerPolynomial/SMZP_Factoring_Support.h"
#include "IntegerPolynomial/DUZP_Support.h"
#include <time.h>

#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
#define FACT_AAZ_VERBOSE 1
#else
#define FACT_AAZ_VERBOSE 0
#endif

//void random_mpz_init(gmp_randstate_t** R_STATE, int** init) {
void random_mpz_init(gmp_randstate_t** R_STATE) {

	static int RM_initRand = 0;
	static gmp_randstate_t RM_R_STATE;
	if (!RM_initRand) {
		time_t t = time(NULL);

		gmp_randinit_default (RM_R_STATE);
		gmp_randseed_ui(RM_R_STATE, t);

		// fprintf(stderr, "fac seed: %lu\n", t);
		RM_initRand = 1;
	}

	*R_STATE = &RM_R_STATE;
//	*init = &RM_initRand;

}

// void chooseRandomEvaluationPoint_factAAZ(mpz_t* values, AltArrZ_t** U0, mpz_t delta, const AltArrZ_t* lc, const AltArrZ_t* f_in, mpz_t bound) {

// 	if (f_in->nvar < 2) {
// 		fprintf(stderr,"BPAS, error: input to chooseRandomEvaluationPoint must be at least bivariate!");
// 	}


// 	int n = f_in->nvar;
// 	int* active = (int*) malloc(sizeof(int)*n);
// 	active[0] = 0;

// 	mpz_t prod, g;
// 	mpz_init(prod);
// 	mpz_init(g);

// 	gmp_randstate_t* RS;
// //	int* init;
// //	random_mpz_init(&RS,&init);
// 	random_mpz_init(&RS);
// 	AltArrZ_t* tmp;
// //	*lc = mainLeadingCoefficient_AAZ (f_in); // TODO: request constant input for this function call
// //	AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
// //	*lc = mainLeadingCoefficient_AAZ (f);
// //	freePolynomial_AAZ(f);

// 	randomize:

// 	mpz_set_ui(prod, 1u);
// 	int maxTries = 20*n;
// 	int curTry = 0;

// 	for (int i=1; i<n && curTry < maxTries; ++i) { // TODO: include a max tries cutoff (in case the bound is too
// 							  //       small and the number of variables is too large

// 		mpz_urandomm(values[i],*RS,bound);
// 		mpz_gcd(g,prod,values[i]);
// 		if (mpz_cmp_ui(g,1) != 0) {
// 			--i;
// 			++curTry;
// 		}
// 		else {
// 			mpz_mul(prod,prod,values[i]);
// 		}
// 	}

// 	for (int i=1; i<n; ++i) {
// 		active[i] = 1;
// 	}


// 	tmp = evaluatePoly_AAZ(lc,active,values,n);

// //	const char* vs2[] = {"x_1", "x_2"};
// //	fprintf(stderr,"lc(a_2) = ");
// //	printPoly_AAZ(stderr,tmp,vs2,tmp->nvar);
// //	fprintf(stderr,"\n");

// 	if (isZero_AAZ(tmp)) {
// 		if (FACT_AAZ_VERBOSE) {
// 			fprintf(stderr,"  [CREP] re-randomizing (lc zero)...\n");
// 		}
// 		freePolynomial_AAZ(tmp);
// 		goto randomize;
// 	}

// 	*U0 = evaluatePoly_AAZ(f_in,active,values,n);
// 	const char* vs2[] = {"x_1", "x_2"};
// 	if (FACT_AAZ_VERBOSE) {
// 		fprintf(stderr,"  [CREP] U0 = ");
// 		printPoly_AAZ(stderr,*U0,vs2,(*U0)->nvar);
// 		fprintf(stderr,"\n");
// 	}

// 	freePolynomial_AAZ(tmp);
// 	tmp = primitivePartAndContent_AAZ(*U0, delta);
// 	freePolynomial_AAZ(*U0);
// 	*U0 = tmp;
// 	tmp = squareFreePart_AAZ(*U0,(*U0)->nvar);


// 	// TODO: ensure that u is squarefree
// 	if (!isExactlyEqual_AAZ(*U0,tmp)) {
// 		if (FACT_AAZ_VERBOSE) {
// 			fprintf(stderr,"  [CREP] re-randomizing (U0 not sqrfree)...\n");
// //			fprintf(stderr,"U0 nvar = %d, tmp nvar = %d\n",(*U0)->nvar,tmp->nvar);
// 		}
// 		freePolynomial_AAZ(tmp);
// 		goto randomize;
// 	}

// 	mpz_clear(prod);
// 	mpz_clear(g);

// 	free(active);
// //	random_mpz_free(RS,init);
// 	freePolynomial_AAZ(tmp);
// }


void chooseRandomEvaluationPointHensel_AAZ(AltArrZ_t const*const* nonzero, int nnzero, mpz_t bound, mpz_t* values, int nvar) {

	gmp_randstate_t* RS;
	random_mpz_init(&RS);

	mpz_t prod, g;
	mpz_init(prod);
	mpz_init(g);

	AltArrZ_t* tmp;
	int active[nvar];
	active[0] = 0;
	for (int i=1; i<nvar; ++i) {
		active[i] = 1;
	}

	int good = 0;
	while(!good) {
		mpz_set_ui(prod, 1u);
		int maxTries = 20*nvar;
		int curTry = 0;

		for (int i=1; i<nvar && curTry < maxTries; ++i) {
			mpz_urandomm(values[i],*RS,bound);
			mpz_gcd(g,prod,values[i]);
			if (mpz_cmp_ui(g,1) != 0) {
				--i;
				++curTry;
			}
			else {
				mpz_mul(prod,prod,values[i]);
			}
		}

		good = 1;
		for (int i = 0; i < nnzero; ++i) {
			tmp = evaluatePoly_AAZ(nonzero[i],active,values,nvar);
			if (isZero_AAZ(tmp)) {
				freePolynomial_AAZ(tmp);
				good = 0;
				if (FACT_AAZ_VERBOSE) {
					fprintf(stderr,"  [CREP] re-randomizing (nonzero was zero)...\n");
				}
				mpz_mul_2exp(bound, bound, 1); //double bound and try again
				break;
			}
			freePolynomial_AAZ(tmp);
		}
	}

	mpz_clear(prod);
	mpz_clear(g);
}


/**
 * Algorithm N "Nondivisors" from Wang, 1978.
 */
int checkEvaluationPoint_factAAZ(mpz_t* d, const mpz_t Omega, const mpz_t* F_tilde, long k, const mpz_t delta) {

	mpz_t d0;
	mpz_t q;
	mpz_t r;
	mpz_init(d0);
	mpz_init(q);
	mpz_init(r);
	mpz_mul(d0,delta,Omega);
	if (FACT_AAZ_VERBOSE) {
		gmp_fprintf(stderr,"  [CSI] d_0 = %Zd\n",d0);
	}

	int ret = 1;
	for (int i=1; i<=k; ++i) {

		mpz_abs(q,F_tilde[i-1]);
		if (FACT_AAZ_VERBOSE) {
			gmp_fprintf(stderr,"  [CSI] A1: q = %Zd\n",q);
		}

		for (int j=i-1; j>=0; --j) {

			if (j == 0) {
				mpz_set(r,d0);
			}
			else {
				mpz_set(r,d[j-1]);
			}

			if (FACT_AAZ_VERBOSE) {
				gmp_fprintf(stderr,"  [CSI] B1: r = %Zd\n",r);
			}

			while(mpz_cmp_ui(r, 1) != 0) {
				mpz_gcd(r,r,q);
				mpz_divexact(q,q,r);

				if (FACT_AAZ_VERBOSE) {
					gmp_fprintf(stderr,"  [CSI] B2: r = %Zd\n",r);
					gmp_fprintf(stderr,"  [CSI] B2: q = %Zd\n",q);
				}

			}
			if (mpz_cmp_ui(q,1) == 0) {
				// special integers cannot be generated given the input
				ret = 0;
				break;
			}
		}

		mpz_set(d[i-1],q);
		if (FACT_AAZ_VERBOSE) {
			gmp_fprintf(stderr,"  [CSI] B1: d_%d = %Zd\n",i,(*d)[i-1]);
		}
	}


	mpz_clear(d0);
	mpz_clear(q);
	mpz_clear(r);

	return ret; // special integers generated successfully

}


int reconstructLeadingCoefficients_factAAZ(AltArrZ_t** lcs, AltArrZ_t* f, AltArrZ_t const*const* lc_factors, const int* lc_exps, int nlcfs,
	const mpz_t* F_tilde, DUZP_t** U0_factors, int nufs, const mpz_t delta_in, const mpz_t* values)
{
	if (lcs == NULL) {
		return 0;
	}

	int nvar = f->nvar;
	mpz_t one;
	mpz_init_set_ui(one,1);
	for (int i=0; i<nufs; ++i) {
		lcs[i] = makeConstPolynomial_AAZ(1, nvar, one);
	}

	// for (int i = 0; i < nlcfs; ++i) {
	// 	gmp_fprintf(stderr, "F[%d]: %Zd\n", i, F_tilde[i]);
	// }
	mpz_t* lc_terms = (mpz_t*) malloc(sizeof(mpz_t)*nufs);
	for (int i=0; i<nufs; ++i) {
		const DUZP_t* poly = U0_factors[i];
		mpz_init_set(lc_terms[i],poly->coefs[poly->lt]);
		mpz_mul(lc_terms[i],lc_terms[i],delta_in);
		// gmp_fprintf(stderr, "lc_terms[%d]: %Zd\n", i, lc_terms[i]);
	}
	int* alloc_count = (int*) calloc(nlcfs,sizeof(int));
	mpz_t temp;
	mpz_t delta;
	mpz_t d;
	mpz_init(temp);
	mpz_init_set(delta,delta_in);
	mpz_init(d);

	int sign = mpz_sgn(delta);
	if (sign < 0) {
		mpz_neg(delta, delta);
	}


	// "allocate"/"distribute" factors of the leading coefficient of the input polynomial to factors of the univariate image U0 of f
	for (int i=nlcfs-1; i>-1; --i) {
		for (int j=0; j<nufs; ++j) {
			while ( mpz_divisible_p(lc_terms[j],F_tilde[i]) ) {
				mpz_divexact(lc_terms[j],lc_terms[j],F_tilde[i]);
				lcs[j] = multiplyPolynomials_AAZ_inp(lcs[j],lc_factors[i],nvar);
				alloc_count[i]++;
			}
		}
	}

	// check whether the allocation was successful
	for (int i=0; i<nlcfs; ++i) {
		if (alloc_count[i] != lc_exps[i]) {
			// allocation of factors failed
			return 0;
		}
	}

	if (sign < 0) {
		negatePolynomial_DUZP_inp(U0_factors[nufs-1]);
	}

	/// distribute the content of U0  ///

	// compute D_tilde (evaluate elements of lcs)
	mpz_t* D_tilde = lc_terms;
	for (int i=0; i<nufs; ++i) {
		evalPolyToVal_AAZ(lcs[i], values, nvar, D_tilde[i]);
	}

	// distribute the integer content in the best way possible
	//
	// AB, 2020-10-29: This is just what was left in lc_terms after dividing by
	// numerous F_tildes..
	if (mpz_cmp_ui(delta,1) == 0) { // if delta == 1, multiply lcs by the lcs of U0_factors/D_tilde
		for (int i=0; i<nufs; ++i) {
			DUZP_t* poly = U0_factors[i];
			mpz_divexact(temp,poly->coefs[poly->lt],D_tilde[i]);
			multiplyByInteger_AAZ_inp(lcs[i], temp);
		}
	}
	else { // otherwise, use algo on p. 1219
		for (int i=0; i<nufs; ++i) {
			DUZP_t* poly = U0_factors[i];
			mpz_gcd(d,poly->coefs[poly->lt],D_tilde[i]);				// d = GCD(lc(u_i),D~_i)
			mpz_divexact(temp,poly->coefs[poly->lt],d);
			multiplyByInteger_AAZ_inp(lcs[i], temp); // C_i = lc(u_i)/d)*D_i

			mpz_divexact(temp,D_tilde[i],d);
			multiplyByInteger_DUZP_inp(U0_factors[i], temp); // u_i = (D~_i/d)*u_i

			mpz_divexact(delta,delta,temp);
		}

		if (mpz_cmp_ui(delta,1) != 0) {
			for (int i=0; i<nufs; ++i) {
				multiplyByInteger_DUZP_inp(U0_factors[i], delta);
				multiplyByInteger_AAZ_inp(lcs[i], delta);
			}
			mpz_pow_ui(temp,delta,nufs-1);
			multiplyByInteger_AAZ_inp(f,temp);
		}
	}

	free(alloc_count);
	for (int i=0; i<nufs; ++i) {
		mpz_clear(D_tilde[i]);
	}
	free(D_tilde); //lc_terms is alias to D_tilde
	mpz_clear(temp);
	mpz_clear(delta);
	mpz_clear(d);

	return 1;

}
