#include <gmpxx.h>
#include <gmp.h>

#include "IntegerPolynomial/SMZP_Support.h"
#include "IntegerPolynomial/SMZP_Hensel.h"
#include "IntegerPolynomial/DUZP_Support.h"
#include "IntegerPolynomial/DUZP_Support_Factoring.hpp"
#include "IntegerPolynomial/SMZP_Support_Factoring.hpp"
#include "IntegerPolynomial/SMZP_Support_Recursive.h"

using namespace SMZP::Factoring;
using namespace DUZP::Factoring;

// const char* vs[] = {"x_1","x_2","x_3","x_4","x_5","x_6","x_7","x_8","x_9","x_10"};
const char* vs[] = {"x", "y", "z", "s", "t", "u", "v", "w", "x_1","x_2","x_3","x_4","x_5","x_6","x_7","x_8","x_9","x_10"};
const char* v1[] = {"x_1"};
const char* v2[] = {"x_2"};
const char* v3[] = {"x_3"};
const char* v4[] = {"x_4"};
const char* v5[] = {"x_5"};
const char* v6[] = {"x_6"};
const char* v7[] = {"x_7"};
const char* v8[] = {"x_8"};
const char* v9[] = {"x_9"};
const char* v10[] = {"x_10"};
char firv[] = "x_1";
char secv[] = "x_2";
char thiv[] = "x_3";
char fouv[] = "x_4";
char fifv[] = "x_5";
char sixv[] = "x_6";
char sevv[] = "x_7";
char eigv[] = "x_8";
char ninv[] = "x_9";
char tenv[] = "x_10";

namespace SMZP {
	namespace Factoring {
#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		static int verbose = 1;
#else 
		static int verbose = 0;
#endif
	}
}

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

//void random_mpz_free(gmp_randstate_t* R_STATE, int* init) {
//	
//	gmp_randclear(*R_STATE);
//	*init = 0;

//}

void SMZP::Factoring::chooseRandomEvaluationPoint(mpz_t** values, AltArrZ_t** U0, mpz_t delta, const AltArrZ_t* lc, const AltArrZ_t* f_in, mpz_t bound) {

	if (f_in->nvar < 2) {
		fprintf(stderr,"BPAS, error: input to chooseRandomEvaluationPoint must be at least bivariate!");
	}
	
	
	int n = f_in->nvar;
	int* active = (int*) malloc(sizeof(int)*n);
	active[0] = 0;
	
	mpz_t prod, g;
	mpz_init(prod);
	mpz_init(g);
	
	gmp_randstate_t* RS;
//	int* init;
//	random_mpz_init(&RS,&init);
	random_mpz_init(&RS);
	AltArrZ_t* tmp;
//	*lc = mainLeadingCoefficient_AAZ (f_in); // TODO: request constant input for this function call
//	AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
//	*lc = mainLeadingCoefficient_AAZ (f);
//	freePolynomial_AAZ(f);
	
	randomize:

	mpz_set_ui(prod, 1u);
	int maxTries = 20*n;
	int curTry = 0;

	for (int i=1; i<n && curTry < maxTries; ++i) { // TODO: include a max tries cutoff (in case the bound is too
							  //       small and the number of variables is too large

		mpz_urandomm((*values)[i],*RS,bound);
		mpz_gcd(g,prod,(*values)[i]);
		if (mpz_cmp_ui(g,1) != 0) {
			--i;
			++curTry;
		}
		else {
			mpz_mul(prod,prod,(*values)[i]);
		}
	}
	
	for (int i=1; i<n; ++i)
		active[i] = 1;

	
	tmp = evaluatePoly_AAZ(lc,active,*values,n);

//	const char* vs2[] = {"x_1", "x_2"};
//	fprintf(stderr,"lc(a_2) = ");
//	printPoly_AAZ(stderr,tmp,vs2,tmp->nvar);
//	fprintf(stderr,"\n");
	
	if (isZero_AAZ(tmp)) {
		if (verbose)
			fprintf(stderr,"  [CREP] re-randomizing (lc zero)...\n");
		freePolynomial_AAZ(tmp);
		goto randomize;
	}
	
	*U0 = evaluatePoly_AAZ(f_in,active,*values,n);
	const char* vs2[] = {"x_1", "x_2"};
	if (verbose) {
		fprintf(stderr,"  [CREP] U0 = ");
		printPoly_AAZ(stderr,*U0,vs2,(*U0)->nvar);
		fprintf(stderr,"\n");
	}
	
	freePolynomial_AAZ(tmp);
	tmp = primitivePartAndContent_AAZ(*U0, delta);
	freePolynomial_AAZ(*U0);
	*U0 = tmp;
	tmp = squareFreePart_AAZ(*U0,(*U0)->nvar);
	
	
	// TODO: ensure that u is squarefree
	if (!isExactlyEqual_AAZ(*U0,tmp)) {
		if (verbose) {
			fprintf(stderr,"  [CREP] re-randomizing (U0 not sqrfree)...\n");
//			fprintf(stderr,"U0 nvar = %d, tmp nvar = %d\n",(*U0)->nvar,tmp->nvar);
		}
		freePolynomial_AAZ(tmp);
		goto randomize;
	}
	
	mpz_clear(prod);
	mpz_clear(g);
	
	free(active);
//	random_mpz_free(RS,init);
	freePolynomial_AAZ(tmp);
	
}

int SMZP::Factoring::chooseSpecialIntegers(mpz_t** d, const mpz_t Omega, const mpz_t* F_tilde, long k, const mpz_t delta) {

	mpz_t d0;
	mpz_t q;
	mpz_t r;
	mpz_init(d0);
	mpz_init(q);
	mpz_init(r);
	mpz_mul(d0,delta,Omega);
	if (verbose)
		gmp_fprintf(stderr,"  [CSI] d_0 = %Zd\n",d0);
	
	for (int i=1; i<=k; ++i) {
	
		mpz_abs(q,F_tilde[i-1]);
		if (verbose)
			gmp_fprintf(stderr,"  [CSI] A1: q = %Zd\n",q);
		
		for (int j=i-1; j>=0; --j) {
			
			if (j == 0)
				mpz_set(r,d0);
			else
				mpz_set(r,(*d)[j-1]);
			
		if (verbose)	
			gmp_fprintf(stderr,"  [CSI] B1: r = %Zd\n",r);
			
			reduce:
			
				mpz_gcd(r,r,q);
				mpz_div(q,q,r);
			
			if (verbose) {	
				gmp_fprintf(stderr,"  [CSI] B2: r = %Zd\n",r);
				gmp_fprintf(stderr,"  [CSI] B2: q = %Zd\n",q);
			}
			
			if (!(mpz_cmp_ui(r,1) == 0)) goto reduce;
			
			if (mpz_cmp_ui(q,1) == 0) return 0; // special integers cannot be generated given the input
		}
		
		mpz_set((*d)[i-1],q);
		if (verbose)
			gmp_fprintf(stderr,"  [CSI] B1: d_%d = %Zd\n",i,(*d)[i-1]);
	}
	
	
	mpz_clear(d0);
	mpz_clear(q);
	mpz_clear(r);
	
	return 1; // special integers generated successfully

}

int reconstructLeadingCoefficients(AltArrZ_t*** lcs, AltArrZ_t* f, AltArrZ_t** lc_factors, const int* lc_exps, int nlcfs, const mpz_t* F_tilde, vec_DUZP_long_t* U0_factors, const mpz_t delta_in, mpz_t* values) {
	// TODO: for some reason I can't make lc_factors here const?


	int nvar = f->nvar;
//	long nlcfs = lc_factors->size;
	long nufs = U0_factors->size;
	// fprintf(stderr,"  [RLC] nlcs = %d\n",nlcfs);
	// fprintf(stderr,"  [RLC] nufs = %ld\n",nufs);
	*lcs = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*nufs);
	mpz_t one;
	mpz_init_set_ui(one,1);
	for (int i=0; i<nufs; ++i) {
//		(*lcs)[i] = makeConstPolynomial_DUZP(1,1);
		(*lcs)[i] = makeConstPolynomial_AAZ(1, nvar, one);
	}
	mpz_t* lc_terms = (mpz_t*) malloc(sizeof(mpz_t)*nufs);
	for (int i=0; i<nufs; ++i) {
		const DUZP_t* poly = U0_factors->pairs[i].a;
		mpz_init_set(lc_terms[i],poly->coefs[poly->lt]);
		mpz_mul(lc_terms[i],lc_terms[i],delta_in);
	}
	int* alloc_count = (int*) calloc(nlcfs,sizeof(int));
	mpz_t temp;
	mpz_t delta;
	mpz_t d;
	mpz_init(temp);
	mpz_init_set(delta,delta_in);
	mpz_init(d);
	
	// allocate factors of the leading coefficient of the input polynomial to factors of the univariate image U0 of f
	if (verbose)
		fprintf(stderr,"  [RLC] allocating factors...");
	for (int i=nlcfs-1; i>-1; --i) {
		for (int j=0; j<nufs; ++j) {
			while ( mpz_divisible_p(lc_terms[j],F_tilde[i]) ) {
				mpz_divexact(lc_terms[j],lc_terms[j],F_tilde[i]);
//				multiplyPolynomials_DUZP_inp(lc_factors->pairs[i].a,&((*lcs)[j]));
				if (verbose) {
					fprintf(stderr,"\n  [RLC] j = %d, lc_factors[%d] = ",j,i);
					printPoly_AAZ(stderr,lc_factors[i],vs,nvar);
				}
				(*lcs)[j] = multiplyPolynomials_AAZ_inp((*lcs)[j],lc_factors[i],nvar);
				alloc_count[i]++;
			}
		}
	}
//	char fvar[] = "x_1";
//	char svar[] = "x_2";
	if (verbose) {
		fprintf(stderr,"done.\n");
		for (int i=0; i<nufs; ++i) {
			fprintf(stderr,"  [RLC] lcs[%d] = ",i);
			printPoly_AAZ(stderr,(*lcs)[i],vs,nvar);
			fprintf(stderr,"\n");
		}
	}
	
	// check whether the allocation was successful
	if (verbose)
		fprintf(stderr,"  [RLC] checking allocation...");
	for (int i=0; i<nlcfs; ++i) {
		if (alloc_count[i] != lc_exps[i]) { // allocation of factors failed
			if (verbose)
				fprintf(stderr,"done.\n");
			return 0;
		}
	}
	if (verbose)
		fprintf(stderr,"done.\n");
	
	/// distribute the content of U0  ///
	
	// compute D_tilde (evaluate elements of *lcs)
	if (verbose)
		fprintf(stderr,"  [RLC] computing D_tilde...");
	mpz_t* D_tilde = (mpz_t*) malloc(sizeof(mpz_t)*nufs); // TODO: consider re-using the lc_terms array for this
//	mpz_set(temp,values[1]);
	for (int i=0; i<nufs; ++i) {
		mpz_init(D_tilde[i]);
//		evaluate_DUZP((*lcs)[i], temp, D_tilde[i]);
		evalPolyToVal_AAZ((*lcs)[i], values, nvar, D_tilde[i]); // TODO: this should work with values as const!
	}
	if (verbose) {
		fprintf(stderr,"done.\n");
//		for (int i=0; i<nufs; ++i) {
//			gmp_fprintf(stderr,"  [RLCB] lcs[%d](%Zd) = %Zd\n",i,temp,D_tilde[i]);
//		}
	}
	
	// distribute the integer content in the best way possible
	if (verbose)
		fprintf(stderr,"  [RLC] distributing integer content...");
	if (mpz_cmp_ui(delta,1) == 0) { // if delta == 1, multiply *lcs by the lcs of U0_factors/D_tilde
	
		for (int i=0; i<nufs; ++i) {
			const DUZP_t* poly = U0_factors->pairs[i].a;
			mpz_divexact(temp,poly->coefs[poly->lt],D_tilde[i]);
//			if (verbose)
//				gmp_fprintf(stderr,"  [RLC] lc(u_%d)/D~_%d = %Zd\n",i,i,temp);
//			multiplyByInteger_DUZP_inp((*lcs)[i],temp);
			multiplyByInteger_AAZ_inp((*lcs)[i], temp);
		}
	}
	else { // otherwise, use algo on p. 1219
		
		for (int i=0; i<nufs; ++i) {
			const DUZP_t* poly = U0_factors->pairs[i].a;
			mpz_gcd(d,poly->coefs[poly->lt],D_tilde[i]);				// d = GCD(lc(u_i),D~_i)
			mpz_divexact(temp,poly->coefs[poly->lt],d);
//			multiplyByInteger_DUZP_inp((*lcs)[i],temp);					// C_i = lc(u_i)/d)*D_i
			multiplyByInteger_AAZ_inp((*lcs)[i], temp);
			
			mpz_divexact(temp,D_tilde[i],d);
			multiplyByInteger_DUZP_inp(U0_factors->pairs[i].a,temp);	// u_i = (D~_i/d)*u_i
//			multiplyByInteger_AAZ_inp(U0_factors[i], temp);
			
			mpz_divexact(delta,delta,temp);
		}
		
		if (mpz_cmp_ui(delta,1) != 0) {
			// fprintf(stderr,"  [RLC] exceptional case: delta != 1");
			for (int i=0; i<nufs; ++i) {
				multiplyByInteger_DUZP_inp(U0_factors->pairs[i].a,delta);
//				multiplyByInteger_DUZP_inp((*lcs)[i],delta);
//				multiplyByInteger_AAZ_inp(U0_factors[i], delta);
				multiplyByInteger_AAZ_inp((*lcs)[i], delta);
			}
			mpz_pow_ui(temp,delta,nufs-1);
			multiplyByInteger_AAZ_inp(f,temp);
		}
	}
	if (verbose) {
		fprintf(stderr,"done.\n");
		for (int i=0; i<nufs; ++i) {
			fprintf(stderr,"  [RLC] lcs[%d] = ",i);
			printPoly_AAZ(stderr,(*lcs)[i],vs,(*lcs)[i]->nvar);
			fprintf(stderr,"\n");
//			printPoly_DUZP((*lcs)[i],svar);
		}
//		const char* vars[] = {"x_1","x_2"};
		fprintf(stderr,"f = ");
		printPoly_AAZ(stderr,f,vs,f->nvar);
		fprintf(stderr,"\n");
	}
	
	free(alloc_count);
	for (int i=0; i<nufs; ++i) {
		mpz_clear(lc_terms[i]);
		mpz_clear(D_tilde[i]);
	}	
	free(lc_terms);
	free(D_tilde);
	mpz_clear(temp);
	mpz_clear(delta);
	mpz_clear(d);
	
	return 1;

}

int SMZP::Factoring::reconstructLeadingCoefficientBivariate(DUZP_t*** lcs, AltArrZ_t* f, const vec_DUZP_long_t* lc_factors, const mpz_t* F_tilde, vec_DUZP_long_t* U0_factors, const mpz_t delta_in, const mpz_t* values) {


	long nlcfs = lc_factors->size;
	long nufs = U0_factors->size;
	*lcs = (DUZP_t**) malloc(sizeof(DUZP_t*)*nufs);
	for (int i=0; i<nufs; ++i) {
		(*lcs)[i] = makeConstPolynomial_DUZP(1,1);
	}
	mpz_t* lc_terms = (mpz_t*) malloc(sizeof(mpz_t)*nufs);
	for (int i=0; i<nufs; ++i) {
		const DUZP_t* poly = U0_factors->pairs[i].a;
		mpz_init_set(lc_terms[i],poly->coefs[poly->lt]);
		mpz_mul(lc_terms[i],lc_terms[i],delta_in);
	}
	int* alloc_count = (int*) calloc(nlcfs,sizeof(int));
	mpz_t temp;
	mpz_t delta;
	mpz_t d;
	mpz_init(temp);
	mpz_init_set(delta,delta_in);
	mpz_init(d);
	
	// allocate factors of the leading coefficient of the input polynomial to factors of the univariate image U0 of f
	if (verbose)
		fprintf(stderr,"  [RLCB] allocating factors...");
	for (int i=nlcfs-1; i>-1; --i) {
		for (int j=0; j<nufs; ++j) {
			while ( mpz_divisible_p(lc_terms[j],F_tilde[i]) ) {
				mpz_divexact(lc_terms[j],lc_terms[j],F_tilde[i]);
				multiplyPolynomials_DUZP_inp(lc_factors->pairs[i].a,&((*lcs)[j]));
				alloc_count[i]++;
			}
		}
	}
	char fvar[] = "x_1";
	char svar[] = "x_2";
	if (verbose) {
		fprintf(stderr,"done.\n");
		for (int i=0; i<nufs; ++i) {
			fprintf(stderr,"  [RLCB] lcs[%d] = ",i);
			printPoly_DUZP((*lcs)[i],svar);
		}
	}
	
	// check whether the allocation was successful
	if (verbose)
		fprintf(stderr,"  [RLCB] checking allocation...");
	for (int i=0; i<nlcfs; ++i) {
		if (alloc_count[i] != lc_factors->pairs[i].b) { // allocation of factors failed
			if (verbose)
				fprintf(stderr,"done.\n");
			return 0;
		}
	}
	if (verbose)
		fprintf(stderr,"done.\n");
	
	/// distribute the content of U0  ///
	
	// compute D_tilde (evaluate elements of *lcs)
	if (verbose)
		fprintf(stderr,"  [RLCB] computing D_tilde...");
	mpz_t* D_tilde = (mpz_t*) malloc(sizeof(mpz_t)*nufs); // TODO: consider re-using the lc_terms array for this
	mpz_set(temp,values[1]);
	for (int i=0; i<nufs; ++i) {
		mpz_init(D_tilde[i]);
		evaluate_DUZP((*lcs)[i], temp, D_tilde[i]);
	}
	if (verbose) {
		fprintf(stderr,"done.\n");
//		for (int i=0; i<nufs; ++i) {
//			gmp_fprintf(stderr,"  [RLCB] lcs[%d](%Zd) = %Zd\n",i,temp,D_tilde[i]);
//		}
	}
	
	// distribute the integer content in the best way possible
	if (verbose)
		fprintf(stderr,"  [RLCB] distributing integer content...");
	if (mpz_cmp_ui(delta,1) == 0) { // if delta == 1, multiply *lcs by the lcs of U0_factors/D_tilde
	
		for (int i=0; i<nufs; ++i) {
			const DUZP_t* poly = U0_factors->pairs[i].a;
			mpz_divexact(temp,poly->coefs[poly->lt],D_tilde[i]);
//			if (verbose)
//				gmp_fprintf(stderr,"  [RLCB] lc(u_%d)/D~_%d = %Zd\n",i,i,temp);
			multiplyByInteger_DUZP_inp((*lcs)[i],temp);
		}
	}
	else { // otherwise, use algo on p. 1219
		
		for (int i=0; i<nufs; ++i) {
			DUZP_t* poly = U0_factors->pairs[i].a;
			mpz_gcd(d,poly->coefs[poly->lt],D_tilde[i]);	// d = GCD(lc(u_i),D~_i)
			mpz_divexact(temp,poly->coefs[poly->lt],d);
			multiplyByInteger_DUZP_inp((*lcs)[i],temp);		// C_i = lc(u_i)/d)*D_i
			
			mpz_divexact(temp,D_tilde[i],d);
			multiplyByInteger_DUZP_inp(poly,temp);			// u_i = (D~_i/d)*u_i
			
			mpz_divexact(delta,delta,temp);
		}
		
		if (mpz_cmp_ui(delta,1) != 0) {
			// fprintf(stderr,"  [RLCB] exceptional case: delta != 1");
			for (int i=0; i<nufs; ++i) {
				multiplyByInteger_DUZP_inp(U0_factors->pairs[i].a,delta);
				multiplyByInteger_DUZP_inp((*lcs)[i],delta);
			}
			mpz_pow_ui(temp,delta,nufs-1);
			multiplyByInteger_AAZ_inp(f,temp);
		}
	}
	if (verbose) {
		fprintf(stderr,"done.\n");
//		for (int i=0; i<nufs; ++i) {
//			fprintf(stderr,"  [RLCB] lcs[%d] = ",i);
//			printPoly_DUZP((*lcs)[i],svar);
//		}
//		const char* vars[] = {"x_1","x_2"};
//		fprintf(stderr,"f = ");
//		printPoly_AAZ(stderr,f,vars,f->nvar);
	}
	
	free(alloc_count);
	for (int i=0; i<nufs; ++i) {
		mpz_clear(lc_terms[i]);
		mpz_clear(D_tilde[i]);
	}	
	free(lc_terms);
	free(D_tilde);
	mpz_clear(temp);
	mpz_clear(delta);
	mpz_clear(d);
	
	return 1;

}

void factorBivariate_prim_sqf(const AltArrZ_t* f_in, AltArrZ_t*** factors, int* size) {

	degree_t* degs = (degree_t*) malloc(sizeof(degree_t)*2); // holds the partial degrees of the input
	partialDegrees_AAZ(f_in, degs);
	
	if (degs[0] == 0 && degs[1] == 0) { // constant input
		free(degs);
		return;
	}
	
	const char* vars[] = {"x_1", "x_2"};
	const char* var1[] = {"x_1"};
	const char* var2[] = {"x_2"};
	char fvar[] = "x_1";
	char svar[] = "x_2";
	
	int F_tilde_alloc = 0;
	int d_alloc = 0;
	
	mpz_t delta;
	mpz_t Omega;
	mpz_init(delta);
	mpz_init(Omega);
	AltArrZ_t* u;
	AltArrZ_t* lc;
	AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
	lc = mainLeadingCoefficient_AAZ (f);
//	freePolynomial_AAZ(f);	
	integralContent_AAZ(lc, Omega);
	
	DUZP_t* U0;
	DUZP_t** fs;
	AltArrZ_t** fsm;
	AltArrZ_t** liftedF;
	long r;
	
	mpz_t cc;
	mpz_init(cc);
	vec_DUZP_long_t U0_factors;
	vec_DUZP_long_t_init2(&U0_factors,0);
	
	mpz_t* values = (mpz_t*) malloc(sizeof(mpz_t)*2);
	for (int i=0; i<2; ++i)
		mpz_init(values[i]);
	
	int success = 0;
//		int ntries = 5;
	int tryno = 1;
	mpz_t(bound);
	mpz_init_set_ui(bound,200);
	mpz_t(ceiling);
	mpz_init_set_ui(ceiling,1000); // actual ceiling will be twice this value

	if (verbose) {
		fprintf(stderr," [FBPSF] Leading coefficient of input:\n [FBPSF] lc(f) = ");
		printPoly_AAZ(stderr,lc,vars,lc->nvar);
		printf("\n");
	}
	
	degree_t* degs2 = (degree_t*) malloc(sizeof(degree_t)*2); // holds the partial degrees of the leading coefficient
	partialDegrees_AAZ(lc, degs2);
			
	// choose prime for lifting
	// TODO: determine how we want to make this choice
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(9223372036854602819);
	
	AltArrZ_t* mulc;
	DUZP_t* ulc;
	int varMap[2] = {0,1};
	
	vec_DUZP_long_t lc_factors;
	vec_DUZP_long_t_init2(&lc_factors,0);
	
	if (degs2[1] == 0) { // trivial leading coefficient
		if (verbose)
			fprintf(stderr," [FBPSF] lc(f) is constant!\n");
		goto find_substitution;
	}
		
	// factor lc(f)
	mulc = deepCopyPolynomial_AAZ(lc); // not needed if we don't need to preserve lc
	reorderVars_AAZ(mulc,varMap,mulc->nvar);
	mulc->nvar = 1;
	ulc = convertFromAltArrZ_DUZP(mulc);
	freePolynomial_AAZ(mulc);
	if (verbose) {
		fprintf(stderr," [FBPSF] Univariate leading coefficient of input:\n [FBPSF] ulc(f) = ");
		printPoly_DUZP(ulc,svar);
	}
	
	DUZP_t** lcs;
	
	if (verbose)
		fprintf(stderr," [FBPSF] calling factorUnivariate:\n");
		
	vec_DUZP_long_t_clear(&lc_factors);
	factor(&cc,&lc_factors,ulc);
	
	freePolynomial_DUZP(ulc);
	
	if (degs[0] == 0) { // f = lc(f), so return factors of lc(f)
		*factors = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*lc_factors.size);
		for (int i=0; i<lc_factors.size; ++i) {
			(*factors)[i] = convertToAltArrZ_DUZP(lc_factors.pairs[i].a); // TODO: return exponentiated lc factors here!!!
			(*factors)[i]->nvar = 2;
			reorderVars_AAZ((*factors)[i], varMap, 2);
		}
		*size = lc_factors.size;
		if (verbose)
			fprintf(stderr," [FBPSF] freeing from early point:\n");
		goto early_point_free;
	}
		
	mpz_t* F_tilde;
	F_tilde = (mpz_t*) malloc(sizeof(mpz_t)*lc_factors.size);
	F_tilde_alloc = 1;
	
	for (int i=0; i<lc_factors.size; ++i) {
		mpz_init(F_tilde[i]);
	}
	
	if (lc_factors.size == 1 && lc_factors.pairs[0].b == 1) {
		if (verbose)
			fprintf(stderr," [FBPSF] lc(f) is irreducible!\n");
		// choose a random eval pt, skip to factorization of U0, check to see that F~ divides only one lc(u_i), and
		// then repeat if necessary
		goto find_substitution;
	}
	
	mpz_t* d;
	d = (mpz_t*) malloc(sizeof(mpz_t)*lc_factors.size);
	d_alloc = 1;
	
	for (int i=0; i<lc_factors.size; ++i) {
		mpz_init(d[i]);
	}
	
	// TODO: check to see if deg(lc,x_2) = 0 or deg(f,x_1) = 0 and handle these cases specially
	
	find_special_substitution:
	
	while (true) {
		if (verbose)
			fprintf(stderr," [FBPSF] try %d \n",tryno);
	
		chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
		
		if (verbose) {
			for (int i=0; i<2; ++i)
				gmp_fprintf(stderr," [FBPSF] values[%d] = %Zd\n",i,values[i]);
			fprintf(stderr," [FBPSF] Evaluated polynomial:\nU0 = ");
			printPoly_AAZ(stderr,u,var1,u->nvar);
			printf("\n");
			gmp_fprintf(stderr," [FBPSF] delta = %Zd\n",delta);
			gmp_fprintf(stderr," [FBPSF] Omega = %Zd\n",Omega);
		
			for (int i=0; i<lc_factors.size; ++i) {
				fprintf(stderr," [FBPSF] mult = %ld, F[%d] = ",lc_factors.pairs[i].b,i);
				printPoly_DUZP(lc_factors.pairs[i].a,svar);
			}
		}
		
		// compute F_tilde values (evaluation of the factors of the lc)
		for (int i=0; i<lc_factors.size; ++i) {
			evaluate_DUZP(lc_factors.pairs[i].a, values[1], F_tilde[i]);
			if (verbose)
				gmp_fprintf(stderr," [FBPSF] F_tilde[%d] = %Zd\n",i,F_tilde[i]);
		}
		
//			for (int i=0; i<lc_factors.size; ++i)
//				gmp_fprintf(stderr,"[FBPSF] d[%d] = %Zd\n",i,d[i]);	
		
		success = chooseSpecialIntegers(&d, Omega, F_tilde, lc_factors.size, delta);
					
		if (success) {
			if (verbose)
				fprintf(stderr," [FBPSF] sucess after %d tries!\n",tryno);
			break;
		}
		
		freePolynomial_AAZ(u);
		
		if (mpz_cmp(bound,ceiling) <= 0)
			mpz_mul_ui(bound,bound,2);
		tryno++;
		
	}
	
	goto univariate_image_factorization;
	
	find_substitution:

	chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
	
	if (verbose) {
		for (int i=0; i<2; ++i)
			gmp_fprintf(stderr," [FBPSF] values[%d] = %Zd\n",i,values[i]);
		fprintf(stderr," [FBPSF] Evaluated polynomial:\nU0 = ");
		printPoly_AAZ(stderr,u,var1,u->nvar);
		printf("\n");
		gmp_fprintf(stderr," [FBPSF] delta = %Zd\n",delta);
		gmp_fprintf(stderr," [FBPSF] Omega = %Zd\n",Omega);
	}
	
	univariate_image_factorization:
	
	// univariate image factorization
	U0 = convertFromAltArrZ_DUZP(u);
	freePolynomial_AAZ(u);
	
	if (verbose) {
		fprintf(stderr," [FBPSF] U0 = ");
		printPoly_DUZP(U0,fvar);
	}
	
	vec_DUZP_long_t_clear(&U0_factors);
	
	if (verbose) {
		fprintf(stderr," [FBPSF] calling factor:\n");
	}
	factor(&cc, &U0_factors, U0);
			
	r = U0_factors.size;
	
	if (verbose) {
		for (int i=0; i<r; ++i) {
			for (int j=0; j<U0_factors.pairs[i].b; ++j) {
				fprintf(stderr," [FBPSF] b = %ld, a = ",U0_factors.pairs[i].b);
				printPoly_DUZP(U0_factors.pairs[i].a,fvar);
			}
		}
	}
	
	if (r == 1) {
		if (verbose)
			fprintf(stderr," [FBPSF] single factor of univariate image of f:\n");
		*factors = (AltArrZ_t**) malloc(sizeof(AltArr_t*));
		*factors[0] = deepCopyPolynomial_AAZ(f);
		*size = 1;
		if (verbose)
			fprintf(stderr," [FBPSF] freeing from midpoint:\n");
		goto midpoint_free;
	}
	
	if (lc_factors.size == 1 && lc_factors.pairs[0].b == 1) {
		evaluate_DUZP(lc_factors.pairs[0].a, values[1], F_tilde[0]);
		if (verbose)
			gmp_fprintf(stderr," [FBPSF] F_tilde[%d] = %Zd\n",0,F_tilde[0]);
		
		int div_count = 0;
		mpz_t lcv;
		mpz_init(lcv);
		for (int i=0; i<r; ++i) {
			if (verbose) {
				mpz_set(lcv,(U0_factors.pairs[i].a)->coefs[(U0_factors.pairs[i].a)->lt]);
				gmp_fprintf(stderr," [FBPSF] lcv = %Zd\n",lcv);
			}
			if ( mpz_divisible_p((U0_factors.pairs[i].a)->coefs[(U0_factors.pairs[i].a)->lt],F_tilde[0]) )
				++div_count;
		}
		mpz_clear(lcv);
		
		if (div_count != 1) {
			if (verbose)
				fprintf(stderr," [FBPSF] singe lc factor, unique divisibility of U0 factors test failed!\n");
			goto find_substitution;
		}
	}
	
	if (degs2[1] != 0) {
		// distribute factors of the leading coefficient of f
		success = reconstructLeadingCoefficientBivariate(&lcs, f, &lc_factors, F_tilde, &U0_factors, delta, values);
		
		if (!success) {
			if (verbose)
				fprintf(stderr," [FBPSF] lc reconstruction failed: choosing new substitution...\n");
			for (int i=0; i<U0_factors.size; ++i) {
				freePolynomial_DUZP(lcs[i]);
			}
			free(lcs);
			if (mpz_cmp(bound,ceiling) <= 0)
				mpz_mul_ui(bound,bound,2);
			if (lc_factors.size == 1 && lc_factors.pairs[0].b == 1)
				goto find_substitution;
			else
				goto find_special_substitution;
		}
	}
	else { // set each leading coefficient to the content of U0 to allow lifting to work
		fprintf(stderr, "\n\nelse branch!\n" );
		lcs = (DUZP_t**) malloc(sizeof(DUZP_t*)*r);
		for (int i=0; i<r; ++i) {
			lcs[i] = makeConstPolynomial_DUZP(1,1);
			multiplyByInteger_DUZP_inp(lcs[i], delta);
			printPoly_DUZP(lcs[i], "Z");
		}
	}

	// package the univariate factors for lifting
	if (verbose)
		fprintf(stderr," [FBPSF] r = %ld\n",r);
		
	fs = (DUZP_t**) malloc(sizeof(DUZP_t*)*r);
	
//		const char* vars[] = {"x", "y"};
	for (int i=0; i<r; ++i) {
		fs[i] = U0_factors.pairs[i].a;
		if (verbose) {
			fprintf(stderr,"[FBPSF] fs[%d] = ",i);
			printPoly_DUZP(fs[i],fvar);
//			printf("\n");
		}
	}
	
	// do the lift
	liftedF = (AltArrZ_t**) malloc(sizeof(AltArr_t*)*r);

	if (verbose)
		fprintf(stderr," [FBPSF] calling bivarHenselLift_AAZ...");
		
	success = bivarHenselLift_AAZ(f, lcs, fs, r, values[1], Pptr, liftedF);
	
	if (verbose)
		fprintf(stderr,"done.\n");
	
	if (!success) {
		if (verbose)
			fprintf(stderr," [FBPSF] Hensel lifting failed: choosing new substitution...\n");
		for (int i=0; i<r; ++i) {
			freePolynomial_DUZP(lcs[i]);
		}
		free(lcs);
		free(liftedF);
		free(fs);
		if (mpz_cmp(bound,ceiling) <= 0)
			mpz_mul_ui(bound,bound,2);
		if (degs2[1] == 0 || (lc_factors.size == 1 && lc_factors.pairs[0].b == 1))
			goto find_substitution;
		else
			goto find_special_substitution;
	}
	
	for (int i=0; i<r; ++i) {
		primitivePart_AAZ_inp(liftedF[i]);
		if (verbose) {
			printPoly_AAZ(stderr,liftedF[i],vars,2);
			printf("\n");
		}
	}
	
	// set output
	
	*factors = liftedF;
	*size = r;
	
	// cleanup
	if (verbose)
		fprintf(stderr," [FBPSF] freeing from end point:\n");
	
	// specific to this point //
	for (int i=0; i<r; ++i) {
		freePolynomial_DUZP(lcs[i]);
	}
	free(lcs);
	free(fs);
	
	// from midpoint //
	
	midpoint_free:
	
	fprintf(stderr, "in midpoint_free\n" );
	vec_DUZP_long_t_free(&U0_factors);
	if (F_tilde_alloc) {
		for (int i=0; i<lc_factors.size; ++i) {
			mpz_clear(F_tilde[i]);
		}
		free(F_tilde);
	}
	if (d_alloc) {
		for (int i=0; i<lc_factors.size; ++i) {
			mpz_clear(d[i]);
		}
		free(d);
	}
	
	// from early point //
	
	early_point_free:
	
	vec_DUZP_long_t_free(&lc_factors);
	free(degs);
	free(degs2);
	mpz_clear(delta);
	mpz_clear(Omega);
	mpz_clear(cc);
	mpz_clear(bound);
	mpz_clear(ceiling);
	for (int i=0; i<2; ++i)
		mpz_clear(values[i]);
	free(values);
	freePolynomial_AAZ(lc);
	freePolynomial_AAZ(f);
	free(Pptr);
		

	fprintf(stderr, "returning from factor\n" );
	return;
}

void factor_prim_sqf(const AltArrZ_t* f_in, AltArrZ_t*** factors, int* size) {


	
	if (verbose) {
		fprintf(stderr,"\n [FACTORPSF] Entering factor_prim_sqf.\n");
		fprintf(stderr," [FACTORPSF] f_in = ");
		printPoly_AAZ(stderr,f_in,vs,f_in->nvar);
		fprintf(stderr,"\n");
	}
	
	int nvar = f_in->nvar;
	degree_t* degs = (degree_t*) malloc(sizeof(degree_t)*nvar); // holds the partial degrees of the input
	partialDegrees_AAZ(f_in, degs);
//	for (int i=0; i<nvar; ++i)
//		fprintf(stderr," [FACTORPSF] degs[%d] = %d\n",i,degs[i]);
	
	int degSum = 0;
	for (int i=0; i<nvar; ++i)
		degSum += degs[i];
	
	if (!degSum) { // constant input
		if (verbose)
			fprintf(stderr," [FACTORPSF] Constant input, returning nothing.\n");
		free(degs);
		return;
	}
	
	int F_tilde_alloc = 0;
	int d_alloc = 0;
	int lcs_alloc = 0;
	int lcexps_alloc = 0;
	
	mpz_t delta; 	// content of U0 (evaluation of f at non-main vars)
	mpz_t Omega; 	// content of leading coefficient of input f
	mpz_init(delta);
	mpz_init(Omega);
	AltArrZ_t* u;	// holds multivariate version of U0
	AltArrZ_t* lc;	// holds leading coefficient of f
	AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
	lc = mainLeadingCoefficient_AAZ(f);
	integralContent_AAZ(lc, Omega);
	
	DUZP_t* U0;				// evaluation of f at non-main vars
	DUZP_t** fs;			// 
	AltArrZ_t** fsm;		//
	AltArrZ_t** liftedF;	// holds lifted factors of U0
	long r;					// number of factors of U0
//	int* active;			// array to indicate which variables to evaluate
	
	mpz_t cc;								// stores unneeded content
	mpz_init(cc);
	vec_DUZP_long_t U0_factors;				// container for factors of U0
	vec_DUZP_long_t_init2(&U0_factors,0);
	
	mpz_t* values = (mpz_t*) malloc(sizeof(mpz_t)*nvar); // holds evaluation point
	for (int i=0; i<nvar; ++i)
		mpz_init(values[i]);
	
	int success = 0;				// boolean for whether a routine is successful
//		int ntries = 5;
	int tryno = 1;					// counter for number of evaluation points tried
	mpz_t(bound);					// bound for the range of evaluation point selection
	mpz_init_set_ui(bound,200);		// initially set to 200, increasing to a ceiling
	mpz_t(ceiling);					// half of the ceiling value for evaluation points
	mpz_init_set_ui(ceiling,1000);	// actual ceiling will be twice this value


	AltArrZ_t** lcs;	// holds the factors of the leading coefficient 
	int* lcexps = NULL;		// holds the exponents of the factors of the leading coefficient
	int nlcs = 0;			// the number of factors of the leading coefficient
	AltArrZ_t** rlcs;	// holds the reconstructed leading coefficients of the factors of the input
	int* rlcexps;		// holds the exponents of the reconstructed leading coefficients of the factors of the input
	

	if (verbose) {
		fprintf(stderr," [FACTORPSF] Leading coefficient of input:\n [FBPSF] lc(f) = ");
		printPoly_AAZ(stderr,lc,vs,lc->nvar);
		printf("\n");
	}
	
	degree_t* degs2 = (degree_t*) malloc(sizeof(degree_t)*nvar); // holds the partial degrees of the leading coefficient
	partialDegrees_AAZ(lc, degs2);
	
	int varMap[nvar];	// for removing the main variable from the leading coefficient
	varMap[0] = -1;
	for (int i=1; i<nvar; ++i)
		varMap[i] = i-1;
			
	// choose prime for lifting
	// TODO: determine how we want to make this choice
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(9223372036854602819);
	
	degSum = 0;
	for (int i=1; i<nvar; ++i)
		degSum += degs2[i];
	if (degSum == 0) { // trivial leading coefficient
		if (verbose)
			fprintf(stderr," [FACTORPSF] lc(f) is constant!\n");
		goto find_substitution;
	}
		
	// factor lc(f)
	
	AltArrZ_t* rlc;
	DUZP_t* ulc;
	
	rlc = deepCopyPolynomial_AAZ(lc); // not needed if we don't need to preserve lc; TODO: shrink and expand 
	                                  //                                                   not deepcopy
	shrinkNumVarsAtIdx_AAZ(rlc,0);
	
	if (verbose) {
		fprintf(stderr," [FACTORPSF] Reduced leading coefficient of input:\n [FACTORPSF] rlc(f) = ");
		printPoly_AAZ(stderr,rlc,vs,rlc->nvar);
		fprintf(stderr,"\n");
	}

	if (verbose)
		fprintf(stderr," [FACTORPSF] recursively calling factor:\n");

	factor(rlc,&lcs,&lcexps,&nlcs,cc);
	lcs_alloc = 1;
	lcexps_alloc = 1;
	for (int i=0; i<nlcs; ++i) {
		expandNumVarsLeft_AAZ(lcs[i],nvar);
	}

	freePolynomial_AAZ(rlc);
	if (verbose) {
		for (int i=0; i<nlcs; ++i) {
			fprintf(stderr," [FACTORPSF] [");
			printPoly_AAZ(stderr,lcs[i],vs,lcs[i]->nvar);
			fprintf(stderr,",%d]\n",lcexps[i]);
		}
	}

	/// THIS IS NOT CORRECT AS CURRENTLY IMPLEMENTED /// TODO: really?
	if (degs[0] == 0) { // f = lc(f), so return factors of lc(f)
		// fprintf(stderr, "[FACTORPSF]: degs is degree 0!!!!\n" );
		*factors = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*nlcs);
		for (int i=0; i<nlcs; ++i) {
			(*factors)[i] = lcs[i]; // TODO: return exponentiated lc factors here!!!
//			(*factors)[i]->nvar = 2;
//			reorderVars_AAZ((*factors)[i], varMap, 2);
		}
		*size = nlcs;
		lcs_alloc = 0;

		if (verbose)
			fprintf(stderr," [FACTORPSF] freeing from early point:\n");
		goto early_point_free;
	}
		
	mpz_t* F_tilde;
	F_tilde = (mpz_t*) malloc(sizeof(mpz_t)*nlcs);
	F_tilde_alloc = 1;
	
	for (int i=0; i<nlcs; ++i) {
		mpz_init(F_tilde[i]);
	}
	
	if (nlcs == 1 && lcexps[0] == 1) {
		if (verbose)
			fprintf(stderr," [FACTORPSF] lc(f) is irreducible!\n");
		// choose a random eval pt, skip to factorization of U0, check to see that 
		// F~ divides only one lc(u_i), and then repeat if necessary
		goto find_substitution;
	}
	
	mpz_t* d;
	d = (mpz_t*) malloc(sizeof(mpz_t)*nlcs);
	d_alloc = 1;
	
	for (int i=0; i<nlcs; ++i) {
		mpz_init(d[i]);
	}
	
	// TODO: check to see if deg(lc,x_2) = 0 or deg(f,x_1) = 0 and handle these cases specially
	
	find_special_substitution:
	
//	active = (int*) calloc(nvar,sizeof(int);
//	for (int i=1; i<nvar; ++i)
//		active[i] = 1;
	
	while (true) {
		if (verbose)
			fprintf(stderr," [FACTORPSF] try %d \n",tryno);
	
		chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
		
		if (verbose) {
			for (int i=0; i<nvar; ++i)
				gmp_fprintf(stderr," [FACTORPSF] values[%d] = %Zd\n",i,values[i]);
			fprintf(stderr," [FACTORPSF] Evaluated polynomial:\nU0 = ");
			printPoly_AAZ(stderr,u,vs,u->nvar);
			printf("\n");
			gmp_fprintf(stderr," [FACTORPSF] delta = %Zd\n",delta);
			gmp_fprintf(stderr," [FACTORPSF] Omega = %Zd\n",Omega);
		
			for (int i=0; i<nlcs; ++i) {
				fprintf(stderr," [FACTORPSF] mult = %d, F[%d] = ",lcexps[i],i);
				printPoly_AAZ(stderr,lcs[i],vs,lcs[i]->nvar);
				fprintf(stderr,"\n");
			}
		}
		
		// compute F_tilde values (evaluation of the factors of the lc)
		for (int i=0; i<nlcs; ++i) {
//			evaluate_DUZP(lc_factors.pairs[i].a, values[1], F_tilde[i]);
			evalPolyToVal_AAZ(lcs[i],values, nvar, F_tilde[i]);
//			evaluatePoly_AAZ(lcs[i], active, values, lcs[i]->nvar);
			if (verbose)
				gmp_fprintf(stderr," [FACTORPSF] F_tilde[%d] = %Zd\n",i,F_tilde[i]);
		}
		
//			for (int i=0; i<lc_factors.size; ++i)
//				gmp_fprintf(stderr,"[FBPSF] d[%d] = %Zd\n",i,d[i]);	
		
		success = chooseSpecialIntegers(&d, Omega, F_tilde, nlcs, delta);
					
		if (success) {
			if (verbose)
				fprintf(stderr," [FACTORPSF] sucess after %d tries!\n",tryno);
			break;
		}
		
		freePolynomial_AAZ(u);
		
		if (mpz_cmp(bound,ceiling) <= 0)
			mpz_mul_ui(bound,bound,2);
		tryno++;
		
	}
	
//	free(active);
	
	goto univariate_image_factorization;
	
	find_substitution:

	chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
	
	if (verbose) {
		for (int i=0; i<nvar; ++i)
			gmp_fprintf(stderr," [FACTORPSF] values[%d] = %Zd\n",i,values[i]);
		fprintf(stderr," [FACTORPSF] Evaluated polynomial:\n  U0 = ");
		printPoly_AAZ(stderr,u,vs,u->nvar);
		printf("\n");
		gmp_fprintf(stderr," [FACTORPSF] delta = %Zd\n",delta);
		gmp_fprintf(stderr," [FACTORPSF] Omega = %Zd\n",Omega);
	}
	
	univariate_image_factorization:
	
	// univariate image factorization
	U0 = convertFromAltArrZ_DUZP(u);
	freePolynomial_AAZ(u);
	
	if (verbose) {
		fprintf(stderr," [FACTORPSF] U0 = ");
		printPoly_DUZP(U0,firv);
	}
	
	vec_DUZP_long_t_clear(&U0_factors);
	
	if (verbose) {
		fprintf(stderr," [FACTORPSF] calling factor:\n");
	}
	factor(&cc, &U0_factors, U0);
			
	r = U0_factors.size;
	
	if (verbose) {
		for (int i=0; i<r; ++i) {
			for (int j=0; j<U0_factors.pairs[i].b; ++j) {
				fprintf(stderr," [FACTORPSF] b = %ld, a = ",U0_factors.pairs[i].b);
				printPoly_DUZP(U0_factors.pairs[i].a,firv);
			}
		}
	}
	
	if (r == 1) {
		if (verbose)
			fprintf(stderr," [FACTORPSF] single factor of univariate image of f:\n");
		*factors = (AltArrZ_t**) malloc(sizeof(AltArr_t*));
		(*factors)[0] = deepCopyPolynomial_AAZ(f);
//		fprintf(stderr," [FACTORPSF] f = ");
//		printPoly_AAZ(stderr,f,vs,f->nvar);
//		fprintf(stderr,"\n");
//		fprintf(stderr," [FACTORPSF] (*factors)[0] = ");
//		printPoly_AAZ(stderr,(*factors)[0],vs,(*factors)[0]->nvar);
//		fprintf(stderr,"\n");
		*size = 1;
		if (verbose)
			fprintf(stderr," [FACTORPSF] freeing from midpoint:\n");
		goto midpoint_free;
	}
	
	if (nlcs == 1 && lcexps[0] == 1) { // leading coefficient is irreducible
		evalPolyToVal_AAZ(lcs[0],values, nvar, F_tilde[0]);
		if (verbose)
			gmp_fprintf(stderr," [FACTORPSF] F_tilde[%d] = %Zd\n",0,F_tilde[0]);
		
		int div_count = 0;
		for (int i=0; i<r; ++i) {
			if (verbose) {
				mpz_t lcv;
				mpz_init(lcv);
				mpz_set(lcv,(U0_factors.pairs[i].a)->coefs[(U0_factors.pairs[i].a)->lt]);
				gmp_fprintf(stderr," [FACTORPSF] lcv = %Zd\n",lcv);
				mpz_clear(lcv);
			}
			if ( mpz_divisible_p((U0_factors.pairs[i].a)->coefs[(U0_factors.pairs[i].a)->lt],F_tilde[0]) )
				++div_count;
		}
		
		if (div_count != 1) {
			if (verbose)
				fprintf(stderr," [FACTORPSF] singe lc factor, unique divisibility of U0 factors test failed!\n");
			goto find_substitution;
		}
	}
	
	if (degSum != 0) {
		// distribute factors of the leading coefficient of f
		success = reconstructLeadingCoefficients(&rlcs, f, lcs, lcexps, nlcs, F_tilde, &U0_factors, delta, values);
//		success = 1;
		
		if (!success) {
			if (verbose)
				fprintf(stderr," [FACTORPSF] lc reconstruction failed: choosing new substitution...\n");
			for (int i=0; i<r; ++i) {
				freePolynomial_AAZ(rlcs[i]);
			}
			free(rlcs);
			if (mpz_cmp(bound,ceiling) <= 0)
				mpz_mul_ui(bound,bound,2);
			if (nlcs == 1 && lcexps[0] == 1) // leading coefficient is irreducible
				goto find_substitution;
			else
				goto find_special_substitution;
		}
	}
	else { // set each leading coefficient to the content of U0 to allow lifting to work
		rlcs = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*r);
		for (int i=0; i<r; ++i) {
//			rlcs[i] = makeConstPolynomial_AAZ(1, nvar, delta);
//			multiplyByInteger_DUZP_inp(lcs[i], delta);
			rlcs[i] = multiplyByInteger_AAZ(makeConstPolynomial_AAZ(1, nvar, delta), delta);
		}
	}

	// package the univariate factors for lifting
	if (verbose)
		fprintf(stderr," [FACTORPSF] r = %ld\n",r);
		
	fs = (DUZP_t**) malloc(sizeof(DUZP_t*)*r);
	
	for (int i=0; i<r; ++i) {
		fs[i] = U0_factors.pairs[i].a;
		if (verbose) {
			fprintf(stderr," [FACTORPSF] fs[%d] = ",i);
			printPoly_DUZP(fs[i],firv);
//			printf("\n");
		}
	}
	
	// do the lift
	liftedF = (AltArrZ_t**) malloc(sizeof(AltArr_t*)*r);
//	for (int i=0; i<r; ++i)
//		liftedF = NULL;

	if (verbose)
		fprintf(stderr," [FACTORPSF] calling henselLift_AAZ...");
		
	success = henselLift_AAZ(f, rlcs, fs, r, values, Pptr, liftedF);
//	success = bivarHenselLift_AAZ(f, lcs, fs, r, values[1], Pptr, liftedF);
	
	if (verbose)
		fprintf(stderr,"done.\n");
	
	if (!success) {
		if (verbose)
			fprintf(stderr," [FACTORPSF] Hensel lifting failed: choosing new substitution...\n");
//		exit(1);
		for (int i=0; i<r; ++i) {
			freePolynomial_AAZ(rlcs[i]);
//			freePolynomial_DUZP(lcs[i]);
		}
		free(rlcs);
		free(liftedF);
		free(fs);
		if (mpz_cmp(bound,ceiling) <= 0)
			mpz_mul_ui(bound,bound,2);
		if (degSum == 0 || (nlcs == 1 && lcexps[0] == 1))
			goto find_substitution;
		else
			goto find_special_substitution;
	}
	
	for (int i=0; i<r; ++i) {
		primitivePart_AAZ_inp(liftedF[i]);
		if (verbose) {
			printPoly_AAZ(stderr,liftedF[i],vs,liftedF[i]->nvar);
			printf("\n");
		}
	}
	
	// set output
	
	*factors = liftedF;
	*size = r;

//	// TEMPORARILY RETURN THE INPUT //
//	*factors = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*));
//	(*factors)[0] = deepCopyPolynomial_AAZ(f_in);
//	*size = 1;
//	return;
//	//////////////////////////////////

	
	// cleanup
	if (verbose)
		fprintf(stderr," [FACTORPSF] freeing from end point:\n");
	
	// specific to this point //
	for (int i=0; i<r; ++i) {
		freePolynomial_AAZ(rlcs[i]);
//		freePolynomial_DUZP(lcs[i]);
	}
	free(rlcs);
	free(fs);
	
	// from midpoint //
	
	midpoint_free:
	
	fprintf(stderr,"1\n");
	vec_DUZP_long_t_free(&U0_factors);
	if (F_tilde_alloc) {
		fprintf(stderr,"2\n");
		for (int i=0; i<nlcs; ++i) {
			mpz_clear(F_tilde[i]);
		}
		free(F_tilde);
	}
	if (d_alloc) {
		fprintf(stderr,"3\n");
		for (int i=0; i<nlcs; ++i) {
			mpz_clear(d[i]);
		}
		free(d);
	}
	
	// from early point //
	
	early_point_free:
	
	if (lcs_alloc) {
		fprintf(stderr,"4\n");
		for (int i=0; i<nlcs; ++i) {
			freePolynomial_AAZ(lcs[i]);
		}
		free(lcs);
		free(lcexps);
	}
//	fprintf(stderr,"5\n");
	free(degs);
//	fprintf(stderr,"6\n");
	free(degs2);
//	fprintf(stderr,"7\n");
	mpz_clear(delta);
//	fprintf(stderr,"8\n");
	mpz_clear(Omega);
//	fprintf(stderr,"9\n");
	mpz_clear(cc);
//	fprintf(stderr,"10\n");
	mpz_clear(bound);
//	fprintf(stderr,"11\n");
	mpz_clear(ceiling);
//	fprintf(stderr,"12\n");
	for (int i=0; i<nvar; ++i)
		mpz_clear(values[i]);
//	fprintf(stderr,"13\n");
	free(values);
//	fprintf(stderr,"14\n");
	freePolynomial_AAZ(lc);
//	fprintf(stderr,"15\n");
	freePolynomial_AAZ(f);
//	fprintf(stderr,"16\n");
	free(Pptr);
	fprintf(stderr,"17\n");
		
	return;
}

void SMZP::Factoring::factorBivariate(const AltArrZ_t* f_in, AltArrZ_t*** factors, int** exps, int* size, mpz_t content) {

	if (f_in->nvar != 2) {
		fprintf(stderr,"BPAS, error: input to factorBivariate must be bivariate!");
	}
	
	// fprintf(stderr,"[FACTORBV] Entering factorBivariate:\n");
	
	*size = 0;
	const char* vars[] = {"x_1", "x_2"};
	const char* var1[] = {"x_1"};
	const char* var2[] = {"x_2"};
	char fvar[] = "x_1";
	char svar[] = "x_2";
	
	// remove powers of x_1 and x_2 and make the result primitive 
	// AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
	// AltArrZ_t* temp;
	AltArrZ_t* f = NULL;
	AltArrZ_t* pcontent = commonFactor_AAZ(f_in, &f);
	if (verbose) {
		fprintf(stderr,"[FACTORBV] pcontent = ");
		printPoly_AAZ(stderr,pcontent,vars,pcontent == NULL ? 0 : pcontent->nvar);
		fprintf(stderr,"\n");
		fprintf(stderr,"[FACTORBV] f = ");
		printPoly_AAZ(stderr,f,vars,f->nvar);
		fprintf(stderr,"\n");
	}
	
	degree_t* degs = (degree_t*) malloc(sizeof(degree_t)*2); // holds the partial degrees of the common factor
	degree_t* degsList = (degree_t*) malloc(sizeof(degree_t)*2); // for setting degrees of the common factor monomial
	partialDegrees_AAZ(pcontent, degs);
	for (int i=0; i<2; ++i) { // counting the number of monomial factors
		if (degs[i] != 0)
			(*size)++;
	}
	*factors = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*(*size));
	*exps = (int*) malloc(sizeof(int)*(*size));
	int idx = 0;
	for (int i=0; i<2; ++i) { // package the monomial factors
		if (degs[i] != 0) {
			for (int j=0; j<2; ++j)
				degsList[j] = i==j;
			setDegrees_AAZ_inp(pcontent, 0, degsList, 2);
			(*factors)[idx] = deepCopyPolynomial_AAZ(pcontent);
			(*exps)[idx] = degs[i];
			idx++;
		}
	}
	
//	primitivePartContent_AAZ_inp(f, content);
	
	// do some cleanup
	free(degs);
	free(degsList);
	freePolynomial_AAZ(pcontent);
	
	// squarefree factorization of f
	AltArrZ_t** sqffacts = NULL;
	degree_t* sqfexps = NULL;
	int nsqffacts;
	// fprintf(stderr,"[FACTORBV] calling squareFree_AAZ...");
	squareFree_AAZ(f, content, &sqffacts, &sqfexps, &nsqffacts);
	// fprintf(stderr,"done.\n");
	
	if (verbose) {
		for (int i=0; i<nsqffacts; ++i) {
			fprintf(stderr,"[FACTORBV] [");
			printPoly_AAZ(stderr,sqffacts[i],vars,2);
			fprintf(stderr,",%d]\n",sqfexps[i]);
		}
	}
	
	// factorize each squarefree factor
	AltArrZ_t** irrfactors;
	int irrsize = 0;

	int ringNvar = 2;
	degree_t pdegs[ringNvar];
	int shrinkLeft;
	int trueNvar;

	for (int i=0; i<nsqffacts; ++i) {
		trueNvar = 0;
		partialDegrees_AAZ(sqffacts[i], pdegs);
		if (pdegs[0] == 0) {
			shrinkLeft = 1;
			shrinkNumVarsAtIdx_AAZ(sqffacts[i], 0);
			trueNvar = 1;
		} else if(pdegs[1] == 0) {
			shrinkLeft = 0;
			shrinkNumVarsAtIdx_AAZ(sqffacts[i], 1);
			trueNvar = 1;
		}

		if (trueNvar == 1) {
			DUZP_t* duzp = convertFromAltArrZ_DUZP(sqffacts[i]);
			vec_DUZP_t ex;
			vec_DUZP_t_init2(&ex,0);
			
			//univariate factor
			long bnd = MaxBits(duzp) + (NumBits(duzp->lt+1)+1)/2;
			DUZP::Factoring::factor_prim_sqf(&ex, duzp, bnd);
			
			idx = *size;
			*size = *size + ex.size;
			*factors = (AltArrZ_t**) realloc(*factors,sizeof(AltArrZ_t*)*(*size));
			*exps = (int*) realloc(*exps,sizeof(int)*(*size));
		
			for (int j=0; j<ex.size; ++j) {
				(*factors)[idx + j] = convertToAltArrZ_DUZP(ex.polys[j]);
				if (shrinkLeft) {
					expandNumVarsLeft_AAZ((*factors)[idx + j], 2);
				} else {
					expandNumVars_AAZ((*factors)[idx + j], 2);
				}
				(*exps)[idx + j] = sqfexps[i];
			}
			
			freePolynomial_DUZP(duzp);
			vec_DUZP_t_free(&ex);

		} else {
			factorBivariate_prim_sqf(sqffacts[i], &irrfactors, &irrsize);
			idx = *size;
			*size = *size + irrsize;
			*factors = (AltArrZ_t**) realloc(*factors,sizeof(AltArrZ_t*)*(*size));
			*exps = (int*) realloc(*exps,sizeof(int)*(*size));
			
			for (int j=0; j<irrsize; ++j) {
				(*factors)[idx+j] = irrfactors[j];
				(*exps)[idx+j] = sqfexps[i];
			}
		
			free(irrfactors);
		}
	}
	
	
	// cleanup
	freePolynomial_AAZ(f);
	for (int i=0; i<nsqffacts; ++i)
		freePolynomial_AAZ(sqffacts[i]);
	free(sqffacts);
	free(sqfexps);
	

}

void SMZP::Factoring::factor(const AltArrZ_t* f_in, AltArrZ_t*** factors, int** exps, int* size, mpz_t content) {
	
	// fprintf(stderr,"[FACTOR] Entering factor:\n");
	
	if (f_in->nvar == 1) {
		// fprintf(stderr,"[FACTOR] Univariate input, calling factor from DUZP:\n");
		
		DUZP_t* f = convertFromAltArrZ_DUZP(f_in);
		vec_DUZP_long_t u_factors;
		vec_DUZP_long_t_init2(&u_factors,0);
		mpz_t cc;
		mpz_init(cc);
		
		factor(&cc,&u_factors,f);
		
		mpz_set(content,cc);
		mpz_clear(cc);
		*factors = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*u_factors.size);
		*exps = (int*) malloc(sizeof(int)*u_factors.size);
		*size = u_factors.size;
		
		for (int i=0; i<u_factors.size; ++i) {
			(*factors)[i] = convertToAltArrZ_DUZP(u_factors.pairs[i].a);
			(*exps)[i] = (int) u_factors.pairs[i].b;
		}
		
		freePolynomial_DUZP(f);
		vec_DUZP_long_t_free(&u_factors);
		return;
	}
	
	if (f_in->nvar == 2) {
		// TODO: consider eliminating the high level factorBivariate function, and having 
		//       this function call factorBivariate_prim_sqf at the appropriate place.
		// fprintf(stderr,"[FACTOR] Bivariate input, calling factorBivariate:\n");
		factorBivariate(f_in,factors,exps,size,content);
		return;
	}
	
	int nvar = f_in->nvar;
	*size = 0;
//	const char* vars[] = {"x_1", "x_2"};
//	const char* var1[] = {"x_1"};
//	const char* var2[] = {"x_2"};
//	char fvar[] = "x_1";
//	char svar[] = "x_2";
	
	// remove powers of x_i //
	AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
	AltArrZ_t* temp;
	AltArrZ_t* pcontent = commonFactor_AAZ(f, &temp);
	freePolynomial_AAZ(f);
	f = temp;
	if (verbose) {
		fprintf(stderr,"[FACTOR] pcontent = ");
		printPoly_AAZ(stderr,pcontent,vs,pcontent->nvar);
		fprintf(stderr,"\n");
		fprintf(stderr,"[FACTOR] f = ");
		printPoly_AAZ(stderr,f,vs,f->nvar);
		fprintf(stderr,"\n");
	}
	
	degree_t* degs = (degree_t*) malloc(sizeof(degree_t)*nvar); 	// holds the partial degrees of the common factor
	degree_t* degsList = (degree_t*) malloc(sizeof(degree_t)*nvar); // for setting degrees of the common factor monomial
	partialDegrees_AAZ(pcontent, degs);
	for (int i=0; i<nvar; ++i) { // counting the number of monomial factors
		if (degs[i] != 0)
			(*size)++;
//		fprintf(stderr,"degs[%d] = %d\n",i,degs[i]);
	}
	*factors = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*(*size));
	*exps = (int*) malloc(sizeof(int)*(*size));
	int idx = 0;
	for (int i=0; i<nvar; ++i) { // package the monomial factors
		if (degs[i] != 0) {
			for (int j=0; j<nvar; ++j)
				degsList[j] = i==j; // pick out the current (ith) variable
			setDegrees_AAZ_inp(pcontent, 0, degsList, nvar);
			(*factors)[idx] = deepCopyPolynomial_AAZ(pcontent);
			(*exps)[idx] = degs[i];
			idx++;
		}
	}
	
	// do some cleanup
	free(degs);
	free(degsList);
	freePolynomial_AAZ(pcontent);
	
	// squarefree factorization of f
	AltArrZ_t** sqffacts = NULL;
	degree_t* sqfexps = NULL;
	int nsqffacts;
	fprintf(stderr,"[FACTOR] calling squareFree_AAZ...");
	squareFree_AAZ(f, content, &sqffacts, &sqfexps, &nsqffacts);
	fprintf(stderr,"done.\n");
	
	if (verbose) {
		for (int i=0; i<nsqffacts; ++i) {
			fprintf(stderr,"[FACTOR] [");
			printPoly_AAZ(stderr,sqffacts[i],vs,sqffacts[i]->nvar);
			fprintf(stderr,",%d]\n",sqfexps[i]);
		}
	}
	
	// factorize each squarefree factor
	AltArrZ_t** irrfactors;
	int irrsize = 0;
	
	int varMap[nvar];
	int trueNvar;
	for (int i=0; i<nsqffacts; ++i) {
		trueNvar = tryShrinkVariables_AAZ_inp(sqffacts[i], varMap);
		if (trueNvar == 1) {
			DUZP_t* duzp = convertFromAltArrZ_DUZP(sqffacts[i]);
			vec_DUZP_t ex;
			vec_DUZP_t_init2(&ex,0);
			
			//univariate factor
			long bnd = MaxBits(duzp) + (NumBits(duzp->lt+1)+1)/2;
			DUZP::Factoring::factor_prim_sqf(&ex, duzp, bnd);
			
			idx = *size;
			*size = *size + ex.size;
			*factors = (AltArrZ_t**) realloc(*factors,sizeof(AltArrZ_t*)*(*size));
			*exps = (int*) realloc(*exps,sizeof(int)*(*size));
		
			for (int j=0; j<ex.size; ++j) {
				(*factors)[idx + j] = convertToAltArrZ_DUZP(ex.polys[j]);
				reverseShrinkVariables_AAZ_inp((*factors)[idx + j], nvar, varMap);
				(*exps)[idx + j] = sqfexps[i];
			}
			
			freePolynomial_DUZP(duzp);
			vec_DUZP_t_free(&ex);
		} else {

			if (trueNvar == 2) {
				fprintf(stderr, "\n\n\n\n\n\n\n\n\nbivar factor in factor sqf prim\n\n\n\n\n" );
				factorBivariate_prim_sqf(sqffacts[i], &irrfactors, &irrsize);
			} else {
				factor_prim_sqf(sqffacts[i], &irrfactors, &irrsize);
				fprintf(stderr, "factoring %d of %d\n",i, nsqffacts );
			}
			for (int j = 0; j < irrsize; ++j) {
				reverseShrinkVariables_AAZ_inp(irrfactors[j], nvar, varMap);
			}
		
			idx = *size;
			*size = *size + irrsize;
			*factors = (AltArrZ_t**) realloc(*factors,sizeof(AltArrZ_t*)*(*size));
			*exps = (int*) realloc(*exps,sizeof(int)*(*size));
			
			for (int j=0; j<irrsize; ++j) {
				(*factors)[idx+j] = irrfactors[j];
				(*exps)[idx+j] = sqfexps[i];
			}
			
			free(irrfactors);
		}
	}
	
	
	// cleanup
	freePolynomial_AAZ(f);
	for (int i=0; i<nsqffacts; ++i)
		freePolynomial_AAZ(sqffacts[i]);
	free(sqffacts);
	free(sqfexps);
	
	fprintf(stderr, "returning from SMZP::Factor\n" );
}

//void factorBivariate(const AltArrZ_t* f_in, AltArrZ_t*** factors, long* size) {

//	if (f_in->nvar != 2) {
//		fprintf(stderr,"BPAS, error: input to factorBivariate must be bivariate!");
//	}
//	
//	// make the input primitive 
//	AltArrZ_t* contentX1;
//	AltArrZ_t* contentX2;
//	AltArrZ_t* tmp;
//	AltArrZ_t* f = deepCopyPolynomial_AAZ(f_in);
//	int varMap[] = {1,0};
//	reorderVars_AAZ(f, varMap, 2);
//	fprintf(stderr,"mainPrimitiveFactorization_AAZ: X1\n");
//	tmp = mainPrimitiveFactorization_AAZ(f,&contentX1);
//	reorderVars_AAZ(contentX1, varMap, 2);
//	fprintf(stderr,"free and assign\n");
//	freePolynomial_AAZ(f);
//	f = tmp;
//	reorderVars_AAZ(f, varMap, 2);
//	fprintf(stderr,"mainPrimitiveFactorization_AAZ: X2\n");
//	tmp = mainPrimitiveFactorization_AAZ(f,&contentX2);
//	fprintf(stderr,"free and assign\n");
//	freePolynomial_AAZ(f);
//	f = tmp;
////	fprintf(stderr,"mainPrimitivePart_AAZ X2\n");
////	tmp = mainPrimitivePart_AAZ(f,1);
////	fprintf(stderr,"exactDividePolynomials_AAZ\n");
////	exactDividePolynomials_AAZ(f,tmp,&contentX2,2);
////	fprintf(stderr,"free and assign\n");
////	freePolynomial_AAZ(f);
////	f = tmp;
//	const char* vs2[] = {"x_1", "x_2"};
//	fprintf(stderr,"contentX2 = ");
//	printPoly_AAZ(stderr,contentX2,vs2,2);
//	printf("\n");
//	const char* vs[] = {"x_2", "x_1"};
//	fprintf(stderr,"contentX1 = ");
//	printPoly_AAZ(stderr,contentX1,vs,2);
//	printf("\n");
//	fprintf(stderr,"pp(f) = ");
//	printPoly_AAZ(stderr,f,vs,2);
//	printf("\n");
//	

//	// evaluation point
//	mpz_t x;
//	mpz_init_set_ui(x,7);
//	
//	// set values for evaluation (TODO: determine randomly according to a computed bound)
//	int active[2] = {0,1};
//	mpz_t* vals = (mpz_t*) malloc(sizeof(mpz_t)*2);
//	for (int i=0; i<2; ++i) {
//		mpz_init(vals[i]);
//		mpz_set(vals[i],x);
//	}
//	
//	// do the evaluation
//	fprintf(stderr,"calling evaluatePoly_AAZ:\n");
//	AltArrZ_t* ff = evaluatePoly_AAZ(f,active,vals,2);
////	const char* vs[] = {"x_1", "x_2"};
//	printPoly_AAZ(stderr,ff,vs,2);
//	printf("\n");
//	
//	mpz_t content;
//	mpz_init(content);
//	vec_DUZP_long_t u_factors;
//	vec_DUZP_long_t_init(&u_factors);
//	
//	// make the evaluated poly a DUZP
//	fprintf(stderr,"calling convertFromAltArrZ_DUZP:\n");
//	DUZP_t* fff = convertFromAltArrZ_DUZP(ff);
//	char var[3] = "x\0";
//	printPoly_DUZP(fff,var);
//	
//	// univariate factorization
//	fprintf(stderr,"calling factor:\n");
//	factor(&content, &u_factors, fff);
//	
////	long r = 0;
//	char xv[] = {'x'};
//	for (int i=0; i<u_factors.size; ++i) {
//		for (int j=0; j<u_factors.pairs[i].b; ++j) {
//			fprintf(stderr,"b = %ld, a = ",u_factors.pairs[i].b);
//			printPoly_DUZP(u_factors.pairs[i].a,xv);
//		}
////		r += u_factors.pairs[i].b;
//	}

//	// package the univariate factors for lifting
//	long r = u_factors.size;
//	fprintf(stderr,"r = %ld\n",r);
//	AltArrZ_t** fs = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*r);
//	
//	const char* vars[] = {"x", "y"};
//	for (int i=0; i<u_factors.size; ++i) {
//		fs[i] = convertToAltArrZ_DUZP(u_factors.pairs[i].a);
//		printPoly_AAZ(stderr,fs[i],vars,1);
//		printf("\n");
//	}
////	AltArrZ_t* tmp;
////	int k=0;
////	for (int i=0; i<u_factors.size; ++i) {
////		tmp = convertToAltArrZ_DUZP(u_factors.pairs[i].a);
////		fs[k] = tmp;
////		fprintf(stderr,"k = %d\n",k);
////		printPoly_AAZ(stderr,fs[k],vars,1);
////		printf("\n");
////		for (int j=1; j<u_factors.pairs[i].b; ++j) {
////			k++;
////			fs[k] = deepCopyPolynomial_AAZ(tmp);
////			fprintf(stderr,"k = %d\n",k);
////			printPoly_AAZ(stderr,fs[k],vars,1);
////			printf("\n");
////		}
////		k++;
////	}
//		
//	// choose prime for lifting
////	Prime_ptr* Pptr = smallprimefield_get_prime_constants(11);
//	Prime_ptr* Pptr = smallprimefield_get_prime_constants(9223372036854602819);
//	
//	// do the lift
//	int extra = 0;
//	if (!isConstant_AAZ(contentX1))
//		extra++;
//	if (!isConstant_AAZ(contentX2))
//		extra++;
//	AltArrZ_t** liftedF = (AltArrZ_t**) malloc(sizeof(AltArr_t*)*(r+extra));

//	fprintf(stderr,"calling monicBivarHenselLift_AAZ...");
//	monicBivarHenselLift_AAZ(f, fs, r, x, Pptr, liftedF);
//	fprintf(stderr,"done.\n");
//	
//	for (int i=0; i<r; ++i) {
//		printPoly_AAZ(stderr,liftedF[i],vars,2);
//		printf("\n");
//	}
//	
//	int index = r;
//	if (!isConstant_AAZ(contentX1)) {
//		liftedF[index] = contentX1;
//		index++;
//	}
//	if (!isConstant_AAZ(contentX2)) {
//		liftedF[index] = contentX2;
//	}
//	*size = r+extra;
//	
//	// cleanup
//	fprintf(stderr,"freeing mpz_t:\n");
//	mpz_clear(x);
//	mpz_clear(content);
//	for (int i=0; i<2; ++i)
//		mpz_clear(vals[i]);
//	free(vals);
//	fprintf(stderr,"freeing AAZ_t:\n");
//	for (int i=0; i<r; ++i)
//		freePolynomial_AAZ(fs[i]);
//	free(fs);
//	free(Pptr);
//	fprintf(stderr,"freeing DUZP_long_t:\n");
//	vec_DUZP_long_t_free(&u_factors);
///*	for (int i=0; i<r; ++i)*/
///*		freePolynomial_AAZ(liftedF[i]);*/
///*	free(liftedF);*/
//	fprintf(stderr,"returning factors:\n");
//	*factors = liftedF;
//	

//}
