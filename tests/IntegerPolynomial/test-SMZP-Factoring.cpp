#include <gmpxx.h>

#include "IntegerPolynomial/SMZP_Support_Factoring.hpp"
#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "IntegerPolynomial/SMZP_Support_Test.h"
#include "IntegerPolynomial/SMZP_Support.h"

#include "IntegerPolynomial/DUZP_Support_Factoring.hpp"
#include "IntegerPolynomial/NTL_mpz_DUZP_conversion.hpp"

#include "Parser/bpas_parser.h"

polysize_t maxDeg = 15;
int coefBoundBits = 5ul;
float sparsity = 0.1;
int includeNeg = 1;

char xvar[1] = {'x'};
const char* vars[] = {"x_1","x_2","x_3","x_4","x_5","x_6","x_7","x_8","x_9","x_10"};

extern void init_vars(char*** vars);
extern void factor_test();
extern void bv_factor_test();
extern void factor_test_old();
extern void evaluation_point_test();
extern void special_integers_test();
extern void lc_reconstruction_test();
/*extern void factor_test(int deg, int n, int monic);*/
/*extern void factor_test(int* sizes, int n, int monic);*/

using namespace SMZP::Factoring;
using namespace DUZP::Factoring;

int main() {

	factor_test();
//	bv_factor_test();
//	lc_reconstruction_test();
//	special_integers_test();
//	evaluation_point_test();
//	factor_test_old();
	
	return 0;
}

AltArrZ_t* randomBivariatePoly(int* maxDegs, int nf, int monic) {

	int constantTerm = 0;
	int nvar = 2;
//	const char* vars[] = {"x_1", "x_2"};
	degrees_t degs = maxDegs[0]+1; // one more than the first of max degrees
	degree_t deg[] = {(degree_t) degs,0}; // depends on choice of degs
	degree_t deg2[] = {0,0};
	mpz_t one;
	mpz_init_set_ui(one,1);
	
	AltArrZ_t* a;
	AltArrZ_t* f = makeConstPolynomial_AAZ(1, nvar, one);
	AltArrZ_t* g;
	
	for (int i=0; i<nf; ++i) {
		a = buildRandomZPolyFromMax(nvar, maxDegs, coefBoundBits, sparsity, includeNeg);
		if (monic) {
			setCoefficient_AAZ(a, deg, nvar, one);
			sortPolynomial_AAZ(a);
		}
		if (constantTerm) {
			setCoefficient_AAZ(a, deg2, nvar, one);
			sortPolynomial_AAZ(a);
		}
		printPoly_AAZ(stderr,a,vars,2);
		printf("\n");
		g = multiplyPolynomials_AAZ(a,f,nvar);
		freePolynomial_AAZ(f);
		f = g;
		freePolynomial_AAZ(a);
	}
//	printPoly_AAZ(stderr,f,vars,2);
//	printf("\n");
	
	mpz_clear(one);
	return f;

}

AltArrZ_t* randomBivariatePoly(int* maxDegs, int* exps, int nf, int monic) {

	int constantTerm = 0;
	int nvar = 2;
//	const char* vars[] = {"x_1", "x_2"};
	degrees_t degs = maxDegs[0]+1; // one more than the first of max degrees
	degree_t deg[] = {(degree_t) degs,0}; // depends on choice of degs
	degree_t deg2[] = {0,0};
	mpz_t one;
	mpz_init_set_ui(one,1);
	
	AltArrZ_t* a;
	AltArrZ_t* f = makeConstPolynomial_AAZ(1, nvar, one);
	AltArrZ_t* g;
	AltArrZ_t* tmp;
	
	for (int i=0; i<nf; ++i) {
		a = buildRandomZPolyFromMax(nvar, maxDegs, coefBoundBits, sparsity, includeNeg);
		if (monic) {
			setCoefficient_AAZ(a, deg, nvar, one);
			sortPolynomial_AAZ(a);
		}
		if (constantTerm) {
			setCoefficient_AAZ(a, deg2, nvar, one);
			sortPolynomial_AAZ(a);
		}
		tmp = deepCopyPolynomial_AAZ(a);
		for (int j=1; j<exps[i]; ++j) {
			fprintf(stderr,"multiplying by a, j=%d\n",j);
			g = multiplyPolynomials_AAZ(a,tmp,a->nvar);
			freePolynomial_AAZ(a);
			a = g;
		}
		freePolynomial_AAZ(tmp);
		printPoly_AAZ(stderr,a,vars,2);
		printf("\n");
		g = multiplyPolynomials_AAZ(a,f,nvar);
		freePolynomial_AAZ(f);
		f = g;
		freePolynomial_AAZ(a);
	}
//	printPoly_AAZ(stderr,f,vars,2);
//	printf("\n");
	
	mpz_clear(one);
	return f;

}

AltArrZ_t* randomPoly(int nvar, int* maxDegs, int* exps, int nf, int monic) {

	int constantTerm = 0;
//	const char* vars[] = {"x_1", "x_2"};
	degrees_t degs = maxDegs[0]+1; // one more than the first of max degrees
	degree_t* deg = (degree_t*) calloc(nvar,sizeof(degree_t));
	deg[0] = (degree_t) degs; // depends on choice of degs
	degree_t* deg2 = (degree_t*) calloc(nvar,sizeof(degree_t));
	mpz_t one;
	mpz_init_set_ui(one,1);
	
	AltArrZ_t* a;
	AltArrZ_t* f = makeConstPolynomial_AAZ(1, nvar, one);
	AltArrZ_t* g;
	AltArrZ_t* tmp;
	
	AltArrZ_t* lc;
	
	for (int i=0; i<nf; ++i) {
		a = buildRandomZPolyFromMax(nvar, maxDegs, coefBoundBits, sparsity, includeNeg);
		if (monic) {
			setCoefficient_AAZ(a, deg, nvar, one);
			sortPolynomial_AAZ(a);
		}
		if (constantTerm) {
			setCoefficient_AAZ(a, deg2, nvar, one);
			sortPolynomial_AAZ(a);
		}
		tmp = deepCopyPolynomial_AAZ(a);
		fprintf(stderr,"exps[%d] = %d\n",i,exps[i]);
		for (int j=1; j<exps[i]; ++j) {
			fprintf(stderr,"multiplying by a, j=%d\n",j);
			g = multiplyPolynomials_AAZ(a,tmp,a->nvar);
			freePolynomial_AAZ(a);
			a = g;
		}
		fprintf(stderr,"irreducible factor: a = ");
		printPoly_AAZ(stderr,tmp,vars,tmp->nvar);
		printf("\n");
		lc = mainLeadingCoefficient_AAZ (tmp);
		shrinkNumVarsAtIdx_AAZ(lc,0);
		fprintf(stderr,"lc(a) = ");
		printPoly_AAZ(stderr,lc,&vars[1],lc->nvar);
		printf("\n");
		freePolynomial_AAZ(tmp);
		freePolynomial_AAZ(lc);
		g = multiplyPolynomials_AAZ(a,f,nvar);
		freePolynomial_AAZ(f);
		f = g;
		freePolynomial_AAZ(a);
	}
//	printPoly_AAZ(stderr,f,vars,2);
//	printf("\n");
	
	mpz_clear(one);
	return f;

}

void lc_reconstruction_test() {
	
	mpz_t* values = (mpz_t*) malloc(sizeof(mpz_t)*2);
	for (int i=0; i<2; ++i)
		mpz_init(values[i]);
	int maxDegs[2] = {2,2};
	
	AltArrZ_t* f = randomBivariatePoly(maxDegs,2,0);
	AltArrZ_t* tmp = squareFreePart_AAZ(f,f->nvar);
	freePolynomial_AAZ(f);
	f = tmp;
	tmp = primitivePart_AAZ(f);
	freePolynomial_AAZ(f);
	f = tmp;
//	mpz_t one;
//	mpz_init_set_ui(one,1);
//	AAZElem_t x1;
//	mpz_set(x1.coef,one);
//	int* degsx1[] {1,0};
//	x1.degs = degsx1;
//	AltArrZ_t* x1 = makeConstPolynomial_AAZ(1, 2, one);
	
	
//	const char* vars[] = {"x_1", "x_2"};
	const char* var1[] = {"x_1"};
	const char* var2[] = {"x_2"};
	char fvar[] = "x_1";
	char svar[] = "x_2";
	fprintf(stderr,"Input polynomial:\nf = ");
	printPoly_AAZ(stderr,f,vars,2);
	printf("\n");
	
	mpz_t(bound);
	mpz_init_set_ui(bound,200);
	
	mpz_t delta;
	mpz_init(delta);
	mpz_t Omega;
	mpz_init(Omega);
	
	AltArrZ_t* u;
	AltArrZ_t* lc;
	AltArrZ_t* f2 = deepCopyPolynomial_AAZ(f);
	lc = mainLeadingCoefficient_AAZ (f2);
	freePolynomial_AAZ(f2);
	integralContent_AAZ(lc, Omega);
	
	fprintf(stderr,"Leading coefficient of input:\nlc(f) = ");
	printPoly_AAZ(stderr,lc,vars,lc->nvar);
	printf("\n");
	
	AltArrZ_t* mulc = deepCopyPolynomial_AAZ(lc);
	int varMap[2] = {0,1};
	reorderVars_AAZ(mulc,varMap,mulc->nvar);
	mulc->nvar = 1;
	DUZP_t* ulc = convertFromAltArrZ_DUZP(mulc);
	freePolynomial_AAZ(mulc);
	fprintf(stderr,"Univariate eading coefficient of input:\nulc(f) = ");
	printPoly_DUZP(ulc,svar);
	
	mpz_t cc;
	mpz_init(cc);
	vec_DUZP_long_t lc_factors;
	vec_DUZP_long_t_init2(&lc_factors,0);
	
	fprintf(stderr,"calling BPAS factor:\n");
	
	factor(&cc,&lc_factors,ulc);
	
	mpz_t* F_tilde;
	mpz_t* d;
	
	int ntries = 10;
	int success;
	
	for (int i=0; i<ntries; ++i) {
		fprintf(stderr,"try %d \n",i+1);
	
		chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
		
		for (int i=0; i<2; ++i)
			gmp_fprintf(stderr,"values[%d] = %Zd\n",i,values[i]);
		fprintf(stderr,"Evaluated polynomial:\nU0 = ");
		printPoly_AAZ(stderr,u,var1,u->nvar);
		printf("\n");
		gmp_fprintf(stderr,"delta = %Zd\n",delta);
		gmp_fprintf(stderr,"Omega = %Zd\n",Omega);
		
		for (int i=0; i<lc_factors.size; ++i) {
			fprintf(stderr,"mult = %ld, F[%d] = ",lc_factors.pairs[i].b,i);
			printPoly_DUZP(lc_factors.pairs[i].a,xvar);
		}
		
//		integralContent_AAZ(u, delta);
//		integralContent_AAZ(lc, Omega);
//		gmp_fprintf(stderr,"delta = %Zd\n",delta);
//		gmp_fprintf(stderr,"Omega = %Zd\n",Omega);
		
		F_tilde = (mpz_t*) malloc(sizeof(mpz_t)*lc_factors.size);
		d = (mpz_t*) malloc(sizeof(mpz_t)*lc_factors.size);
		
		for (int i=0; i<lc_factors.size; ++i) {
			mpz_init(F_tilde[i]);
			mpz_init(d[i]);
			evaluate_DUZP(lc_factors.pairs[i].a, values[1], F_tilde[i]);
			gmp_fprintf(stderr,"F_tilde[%d] = %Zd\n",i,F_tilde[i]);
		}
		
		for (int i=0; i<lc_factors.size; ++i)
			gmp_fprintf(stderr,"d[%d] = %Zd\n",i,d[i]);	
		
		success = chooseSpecialIntegers(&d, Omega, F_tilde, lc_factors.size, delta);
		
		
		if (success) {
			fprintf(stderr,"sucess after %d tries!\n",i+1);
			break;
		}
		
		for (int i=0; i<lc_factors.size; ++i) {
			mpz_clear(F_tilde[i]);
			mpz_clear(d[i]);
		}
		free(F_tilde);
		free(d);
		freePolynomial_AAZ(u);
		
		mpz_mul_ui(bound,bound,2);
		
		if (i = ntries)
			return;
		
	}
	
	for (int i=0; i<lc_factors.size; ++i)
		gmp_fprintf(stderr,"d[%d] = %Zd\n",i,d[i]);
	
	DUZP_t* U0 = convertFromAltArrZ_DUZP(u);
	printPoly_DUZP(U0,fvar);
	vec_DUZP_long_t U0_factors;
	vec_DUZP_long_t_init2(&U0_factors,0);
	
	gmp_fprintf(stderr,"delta = %Zd\n",delta);
	gmp_fprintf(stderr,"Omega = %Zd\n",Omega);
	
	if (success) {	
		fprintf(stderr,"Reconstructing leading coefficients of the factors of U0\n");
		
		fprintf(stderr,"calling BPAS factor:\n");
		
		factor(&cc,&U0_factors,U0);
		for (int i=0; i<U0_factors.size; ++i) {
			fprintf(stderr,"U0[%d] = ",i);
			printPoly_DUZP(U0_factors.pairs[i].a,svar);
		}
		
		DUZP_t** lcs;
		
		success = reconstructLeadingCoefficientBivariate(&lcs, f, &lc_factors, F_tilde, &U0_factors, delta, values);
		
		fprintf(stderr,"success = %d\n",success);
		for (int i=0; i<U0_factors.size; ++i) {
			printPoly_DUZP(lcs[i],svar);
		}
		
		for (int i=0; i<U0_factors.size; ++i)
			freePolynomial_DUZP(lcs[i]);
		free(lcs);
	}
	
	vec_DUZP_long_t_free(&U0_factors);
	vec_DUZP_long_t_free(&lc_factors);
	mpz_clear(delta);
	mpz_clear(cc);
	mpz_clear(bound);
	mpz_clear(Omega);
	freePolynomial_AAZ(f);
//	freePolynomial_AAZ(u);
	freePolynomial_AAZ(lc);
	for (int i=0; i<2; ++i)
		mpz_clear(values[i]);
	for (int i=0; i<lc_factors.size; ++i) {
		mpz_clear(F_tilde[i]);
		mpz_clear(d[i]);
	}
	free(values);
	free(F_tilde);
	free(d);
	freePolynomial_AAZ(u);
//	for (int i=0; i<lc_factors.size; ++i) {
//		mpz_clear(F_tilde[i]);
//		mpz_clear(d[i]);
//	}

}

void special_integers_test() {
	
	mpz_t* values = (mpz_t*) malloc(sizeof(mpz_t)*2);
	for (int i=0; i<2; ++i)
		mpz_init(values[i]);
	int maxDegs[2] = {2,2};
	
	AltArrZ_t* f = randomBivariatePoly(maxDegs,2,0);
	AltArrZ_t* tmp = squareFreePart_AAZ(f,f->nvar);
	freePolynomial_AAZ(f);
	f = tmp;
	tmp = primitivePart_AAZ(f);
	freePolynomial_AAZ(f);
	f = tmp;
//	mpz_t one;
//	mpz_init_set_ui(one,1);
//	AAZElem_t x1;
//	mpz_set(x1.coef,one);
//	int* degsx1[] {1,0};
//	x1.degs = degsx1;
//	AltArrZ_t* x1 = makeConstPolynomial_AAZ(1, 2, one);
	
	
//	const char* vars[] = {"x_1", "x_2"};
	const char* var1[] = {"x_1"};
	const char* var2[] = {"x_2"};
	char svar[] = "x_2";
	fprintf(stderr,"Input polynomial:\nf = ");
	printPoly_AAZ(stderr,f,vars,2);
	printf("\n");
	
	mpz_t(bound);
	mpz_init_set_ui(bound,200);
	
	mpz_t delta;
	mpz_init(delta);
	mpz_t Omega;
	mpz_init(Omega);
	
	AltArrZ_t* u;
	AltArrZ_t* lc;
	AltArrZ_t* f2 = deepCopyPolynomial_AAZ(f);
	lc = mainLeadingCoefficient_AAZ (f2);
	freePolynomial_AAZ(f2);
	integralContent_AAZ(lc, Omega);
	
	fprintf(stderr,"Leading coefficient of input:\nlc(f) = ");
	printPoly_AAZ(stderr,lc,vars,lc->nvar);
	printf("\n");
	
	AltArrZ_t* mulc = deepCopyPolynomial_AAZ(lc);
	int varMap[2] = {0,1};
	reorderVars_AAZ(mulc,varMap,mulc->nvar);
	mulc->nvar = 1;
	DUZP_t* ulc = convertFromAltArrZ_DUZP(mulc);
	freePolynomial_AAZ(mulc);
	fprintf(stderr,"Univariate eading coefficient of input:\nulc(f) = ");
	printPoly_DUZP(ulc,svar);
	
	mpz_t cc;
	mpz_init(cc);
	vec_DUZP_long_t uni_factors;
	vec_DUZP_long_t_init2(&uni_factors,0);
	
	fprintf(stderr,"calling BPAS factor:\n");
	
	factor(&cc,&uni_factors,ulc);
	
	mpz_t* F_tilde;
	mpz_t* d;
	
	int ntries = 10;
	int success;
	
	for (int i=0; i<ntries; ++i) {
		fprintf(stderr,"try %d \n",i+1);
	
		chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
		
		for (int i=0; i<2; ++i)
			gmp_fprintf(stderr,"values[%d] = %Zd\n",i,values[i]);
		fprintf(stderr,"Evaluated polynomial:\nU0 = ");
		printPoly_AAZ(stderr,u,var1,u->nvar);
		printf("\n");
		gmp_fprintf(stderr,"delta = %Zd\n",delta);
		gmp_fprintf(stderr,"Omega = %Zd\n",Omega);
		
		for (int i=0; i<uni_factors.size; ++i) {
			fprintf(stderr,"mult = %ld, F[%d] = ",uni_factors.pairs[i].b,i);
			printPoly_DUZP(uni_factors.pairs[i].a,xvar);
		}
		
//		integralContent_AAZ(u, delta);
//		integralContent_AAZ(lc, Omega);
//		gmp_fprintf(stderr,"delta = %Zd\n",delta);
//		gmp_fprintf(stderr,"Omega = %Zd\n",Omega);
		
		F_tilde = (mpz_t*) malloc(sizeof(mpz_t)*uni_factors.size);
		d = (mpz_t*) malloc(sizeof(mpz_t)*uni_factors.size);
		
		for (int i=0; i<uni_factors.size; ++i) {
			mpz_init(F_tilde[i]);
			mpz_init(d[i]);
			evaluate_DUZP(uni_factors.pairs[i].a, values[1], F_tilde[i]);
			gmp_fprintf(stderr,"F_tilde[%d] = %Zd\n",i,F_tilde[i]);
		}
		
		for (int i=0; i<uni_factors.size; ++i)
			gmp_fprintf(stderr,"d[%d] = %Zd\n",i,d[i]);	
		
		success = chooseSpecialIntegers(&d, Omega, F_tilde, uni_factors.size, delta);
		
		for (int i=0; i<uni_factors.size; ++i) {
			mpz_clear(F_tilde[i]);
			mpz_clear(d[i]);
		}
		free(F_tilde);
		free(d);
		freePolynomial_AAZ(u);
		
		if (success) {
			fprintf(stderr,"sucess after %d tries!\n",i+1);
			break;
		}
		
		mpz_mul_ui(bound,bound,2);
		
	}
	
	for (int i=0; i<uni_factors.size; ++i)
		gmp_fprintf(stderr,"d[%d] = %Zd\n",i,d[i]);
	
	vec_DUZP_long_t_free(&uni_factors);
	mpz_clear(delta);
	mpz_clear(cc);
	mpz_clear(bound);
	mpz_clear(Omega);
	freePolynomial_AAZ(f);
//	freePolynomial_AAZ(u);
	freePolynomial_AAZ(lc);
	for (int i=0; i<2; ++i)
		mpz_clear(values[i]);
	free(values);
//	for (int i=0; i<uni_factors.size; ++i) {
//		mpz_clear(F_tilde[i]);
//		mpz_clear(d[i]);
//	}

}

void evaluation_point_test() {
	
	mpz_t* values = (mpz_t*) malloc(sizeof(mpz_t)*2);
	for (int i=0; i<2; ++i)
		mpz_init(values[i]);
	int maxDegs[2] = {3,3};
	
	AltArrZ_t* f = randomBivariatePoly(maxDegs,2,0);
	AltArrZ_t* tmp = squareFreePart_AAZ(f,f->nvar);
	freePolynomial_AAZ(f);
	f = tmp;
	
//	const char* vars[] = {"x_1", "x_2"};
	fprintf(stderr,"Input polynomial:\nf = ");
	printPoly_AAZ(stderr,f,vars,2);
	printf("\n");
	
	mpz_t bound;
	mpz_init_set_ui(bound,44);
	mpz_t delta;
	mpz_init(delta);
	
	AltArrZ_t* u;
	AltArrZ_t* lc;
	AltArrZ_t* f2 = deepCopyPolynomial_AAZ(f);
	lc = mainLeadingCoefficient_AAZ (f2);
	freePolynomial_AAZ(f2);
	
	chooseRandomEvaluationPoint(&values, &u, delta, lc, f, bound);
	
	for (int i=0; i<2; ++i)
		gmp_fprintf(stderr,"values[%d] = %Zd\n",i,values[i]);
	fprintf(stderr,"Evaluated polynomial:\nu = ");
	printPoly_AAZ(stderr,u,vars,2);
	printf("\n");
	fprintf(stderr,"Leading coefficient of input:\nlc(f) = ");
	printPoly_AAZ(stderr,lc,vars,2);
	printf("\n");
	gmp_fprintf(stderr,"delta = %Zd\n",delta);
	
	mpz_clear(bound);
	mpz_clear(delta);
	freePolynomial_AAZ(f);
	freePolynomial_AAZ(u);
	freePolynomial_AAZ(lc);
	for (int i=0; i<2; ++i)
		mpz_clear(values[i]);
	free(values);

}

void factor_test() {

	int nf = 2;
	int nvar = 3;
	int maxDegs[] = {2
	,1,0};	 // size same as nvar
	int exponents[] = {1,1,1,1}; // size same as nf
//	int maxDegs[] = {3,2};
//	int exponents[] = {1,1};
	
//	printPoly_AAZ(stderr,f,vars,2);
//	printf("\n");

	AltArrZ_t* f = randomPoly(nvar, maxDegs, exponents, nf, 0);
	f->elems->degs = 2ll << EXP_OFFSET_1_V3;
	mpz_t cc;
	mpz_init_set_ui(cc, 456ul);
	mpz_set_ui(f->elems->coef, 1234ll);
	AltArrZ_t* g = makeConstPolynomial_AAZ(2, nvar, cc);
	g->elems->degs = 2ll;// << EXP_OFFSET_1_V3;
	g->elems[1].degs = 0ll;
	mpz_init_set_ui(g->elems[1].coef, 1ll);
	g->size = 2;
	f = multiplyPolynomials_AAZ(f, g, nvar);
	
	// nvar = 3;
	// altarr_pack* pack = generate_altarr_pack("43904*x_1^21*x_2 + 2192064*x_1^20*x_2 + 50062656*x_1^19*x_2 + 692621568*x_1^18*x_2 + 6476737632*x_1^17*x_2 + 43203657552*x_1^16*x_2 + 211092107856*x_1^15*x_2 + 760716429696*x_1^14*x_2 + 1991664554792*x_1^13*x_2 + 3571252805844*x_1^12*x_2 + 3474257197596*x_1^11*x_2 - 1489426063056*x_1^10*x_2 - 11716652926950*x_1^9*x_2 - 21147449029017*x_1^8*x_2 - 21083369832141*x_1^7*x_2 - 11290310317152*x_1^6*x_2 - 3005016572014*x_1^5*x_2 + 623182743543*x_1^4*x_2 + 1409552636247*x_1^3*x_2 + 613898076330*x_1^2*x_2 + 73718023050*x_1*x_2");
	// altarr_pack* pack = generate_altarr_pack("43904*x_1^21*x_2*x_3+2192064*x_1^20*x_2*x_3+50062656*x_1^19*x_2*x_3+692621568*x_1^18*x_2*x_3+6476737632*x_1^17*x_2*x_3+43203657552*x_1^16*x_2*x_3+211092107856*x_1^15*x_2*x_3+760716429696*x_1^14*x_2*x_3+1991664554792*x_1^13*x_2*x_3+3571252805844*x_1^12*x_2*x_3+3474257197596*x_1^11*x_2*x_3-1489426063056*x_1^10*x_2*x_3-11716652926950*x_1^9*x_2*x_3-21147449029017*x_1^8*x_2*x_3-21083369832141*x_1^7*x_2*x_3-11290310317152*x_1^6*x_2*x_3-3005016572014*x_1^5*x_2*x_3+623182743543*x_1^4*x_2*x_3+1409552636247*x_1^3*x_2*x_3+613898076330*x_1^2*x_2*x_3+73718023050*x_1*x_2*x_3");
    // freePolynomial_AAZ(f);
	// f = deepCopyPolynomial_AAZFromAA(pack->altarr_t_data);
    
    // freePolynomial_AA(pack->altarr_t_data);
    // free(pack->vars);
    // free(pack);

	fprintf(stderr,"Input polynomial:\n");
	printPoly_AAZ(stderr,f,vars,nvar);
	printf("\n");
	
	AltArrZ_t** factors;
	int* exps;
	int size;
	mpz_t c;
	mpz_init(c);

	factor(f, &factors, &exps, &size, c);
	
	fprintf(stderr,"Factors from factor:\n");
	gmp_fprintf(stderr,"[%Zd]\n",c);
	for (int i=0; i<size; ++i) {
		fprintf(stderr,"[");
		printPoly_AAZ(stderr,factors[i],vars,nvar);
		fprintf(stderr,",%d]\n",exps[i]);
	}
	
	AltArrZ_t* prod = makeConstPolynomial_AAZ(1, nvar, c);
	for (int i=0; i<size; ++i) {
		for (int j=0; j<exps[i]; ++j)
			prod = multiplyPolynomials_AAZ(factors[i],prod,nvar);
	}
	fprintf(stderr,"Product:\n");
	printPoly_AAZ(stderr,prod,vars,nvar);
	printf("\n");
	
	AltArrZ_t* diff = subPolynomials_AAZ(prod, f, nvar);
	fprintf(stderr,"Difference:\n");
	printPoly_AAZ(stderr,diff,vars,nvar);
	printf("\n");
	
	
	freePolynomial_AAZ(f);
	for (int i=0; i<size; ++i)
		freePolynomial_AAZ(factors[i]);
	free(factors);
	free(exps);
	mpz_clear(c);

}

void bv_factor_test() {

	int nf = 2;
	int nvar = 2;
	int maxDegs[2] = {5,9};
	int exponents[2] = {2,1};
//	int maxDegs[2] = {3,2};
//	int exponents[2] = {1,1};
//	const char* vars[] = {"x_1", "x_2"};
	
//	printPoly_AAZ(stderr,f,vars,2);
//	printf("\n");

	AltArrZ_t* f = randomBivariatePoly(maxDegs, exponents, nf, 0);
	
	fprintf(stderr,"Input polynomial:\n");
	printPoly_AAZ(stderr,f,vars,2);
	printf("\n");
	
	AltArrZ_t** factors;
	int* exps;
	int size;
	mpz_t c;
	mpz_init(c);

	factorBivariate(f, &factors, &exps, &size, c);
	
	fprintf(stderr,"Factors from factorBivariate:\n");
	gmp_fprintf(stderr,"[%Zd]\n",c);
	for (int i=0; i<size; ++i) {
		fprintf(stderr,"[");
		printPoly_AAZ(stderr,factors[i],vars,2);
		fprintf(stderr,",%d]\n",exps[i]);
	}
	
	AltArrZ_t* prod = makeConstPolynomial_AAZ(1, nvar, c);
	for (int i=0; i<size; ++i) {
		for (int j=0; j<exps[i]; ++j)
			prod = multiplyPolynomials_AAZ(factors[i],prod,nvar);
	}
	fprintf(stderr,"Product:\n");
	printPoly_AAZ(stderr,prod,vars,2);
	printf("\n");
	
	AltArrZ_t* diff = subPolynomials_AAZ(prod, f, nvar);
	fprintf(stderr,"Difference:\n");
	printPoly_AAZ(stderr,diff,vars,2);
	printf("\n");
	
	
	freePolynomial_AAZ(f);
	for (int i=0; i<size; ++i)
		freePolynomial_AAZ(factors[i]);
	free(factors);
	free(exps);
	mpz_clear(c);

}

void factor_test_old() {

	int nf = 2;
	int nvar = 2;
	int maxDegs[2] = {3,3};
//	const char* vars[] = {"x_1", "x_2"};
	degrees_t degs = 4; // one more than the first of max degrees
	degree_t deg[] = {4,0}; // depends on choice of degs
	mpz_t one;
	mpz_init_set_ui(one,1);
	
//	printPoly_AAZ(stderr,f,vars,2);
//	printf("\n");

	AltArrZ_t* f = randomBivariatePoly(maxDegs, nf, 1);
	fprintf(stderr,"Input polynomial:\n");
	printPoly_AAZ(stderr,f,vars,2);
	printf("\n");
	
	AltArrZ_t** factors;
	int* exps;
	int size;
	mpz_t c;
	mpz_init(c);

	factorBivariate(f, &factors, &exps, &size, c);
	
	fprintf(stderr,"Factors from factorBivariate:\n");
	for (int i=0; i<size; ++i) {
		printPoly_AAZ(stderr,factors[i],vars,2);
		printf("\n");
	}
	AltArrZ_t* prod = makeConstPolynomial_AAZ(1, nvar, one);
	for (int i=0; i<size; ++i) {
		prod = multiplyPolynomials_AAZ(factors[i],prod,nvar);
	}
	fprintf(stderr,"Product:\n");
	printPoly_AAZ(stderr,prod,vars,2);
	printf("\n");
	
	AltArrZ_t* diff = subPolynomials_AAZ(prod, f, nvar);
	fprintf(stderr,"Difference:\n");
	printPoly_AAZ(stderr,diff,vars,2);
	printf("\n");
	
	
	freePolynomial_AAZ(f);
	for (int i=0; i<size; ++i)
		freePolynomial_AAZ(factors[i]);
	mpz_clear(one);
	mpz_clear(c);
	

}

//void factor_test(int deg, int n, int monic) {

//	DUZP_t* f;
//	if (n > 1) {
//		DUZP_t** F = makePolynomials(deg,n,monic);
//		
//		f = makeConstPolynomial_DUZP(1,1);
//		
//		for (int i=0; i<n; ++i) {
//			multiplyPolynomials_DUZP_inp(F[i],&f);
//			freePolynomial_DUZP(F[i]);
//		}
//		free(F);
//			
//	}
//	else if (n == 1) {
//		f = buildRandomPoly_DUZP(deg,coefBoundBits,sparsity,includeNeg);
//		if (monic) {
//			mpz_t one;
//			mpz_init(one);
//			mpz_set_ui(one,1);
//			mpz_set(f->coefs[f->lt],one);
//			mpz_clear(one);
//		}
//	}
//	else {
//		fprintf(stderr,"Invalid number of factors!");
//		exit(1);
//	}
//	
//	printPoly_DUZP(f,xvar);
//	
//	NTL::ZZX ff;
//	NTLZZX_set_DUZP_t(ff,f);
//	NTL::ZZ c;
//	NTL::vec_pair_ZZX_long factors;
//	
////	factor_data_t fd;
//	factor(c, factors, ff, 1, 0);
////	factor(c, factors, ff, 1, 0, fd);
//	
//	DUZP_t* g = makeConstPolynomial_DUZP(1,0);
//	for (int i=1; i<=factors.length(); ++i) {
//		DUZP_t_set_NTLZZX(&g,factors(i).a);
//		printPoly_DUZP(g,xvar);
//	}
//	
//	freePolynomial_DUZP(g);
//	
//	mpz_t cc;
//	mpz_init(cc);
//	vec_DUZP_long_t factors2;
//	vec_DUZP_long_t_init2(&factors2,0);
//	
//	fprintf(stderr,"calling BPAS factor:\n");
//	
//	factor(&cc,&factors2,f);
//	
//	g = makeConstPolynomial_DUZP(1,1);
//	DUZP_t* h;
//	for (int i=0; i<factors2.size; ++i) {
//		h = makeConstPolynomial_DUZP(1,1);
//		fprintf(stderr,"exp = %ld\n",factors2.pairs[i].b);
//		for (int j=0; j<factors2.pairs[i].b; ++j)
//			multiplyPolynomials_DUZP_inp(factors2.pairs[i].a,&h);
//		multiplyPolynomials_DUZP_inp(h,&g);
//		printPoly_DUZP(h,xvar);
//		freePolynomial_DUZP(h);
//	}
//	multiplyByInteger_DUZP_inp(g,cc);
//	
//	fprintf(stderr,"g = ");
//	printPoly_DUZP(g,xvar);
//	fprintf(stderr,"f = ");
//	printPoly_DUZP(f,xvar);
//	g = subtractPolynomials_DUZP(f,g);
//	fprintf(stderr,"f - g = ");
//	printPoly_DUZP(g,xvar);
//	
//	
////	print_factor_data(fd);
//	
//	mpz_clear(cc);
//	freePolynomial_DUZP(f);
//	freePolynomial_DUZP(g);
//	vec_DUZP_long_t_free(&factors2);

//}

//void factor_test(int* sizes, int n, int monic) {

//	DUZP_t** F = makePolynomials(sizes,n,monic);
//	
//	DUZP_t* f = makeConstPolynomial_DUZP(1,1);
//	
//	for (int i=0; i<n; ++i)
//		multiplyPolynomials_DUZP_inp(F[i],&f);

////	DUZP_t* f = buildRandomPoly_DUZP(16,coefBoundBits,sparsity,includeNeg);
//	
//	printPoly_DUZP(f,xvar);
//	
//	NTL::ZZX ff;
//	NTLZZX_set_DUZP_t(ff,f);
//	NTL::ZZ c;
//	NTL::vec_pair_ZZX_long factors;
//	
////	factor_data_t fd;
//	factor(c, factors, ff, 1, 0);
//	
//	DUZP_t* g = makeConstPolynomial_DUZP(1,0);
//	for (int i=1; i<=factors.length(); ++i) {
//		DUZP_t_set_NTLZZX(&g,factors(i).a);
//		printPoly_DUZP(g,xvar);
//	}
//	
//	freePolynomial_DUZP(g);
//	
//	mpz_t cc;
//	mpz_init(cc);
//	vec_DUZP_long_t factors2;
//	vec_DUZP_long_t_init2(&factors2,0);
//	
//	factor(&cc,&factors2,f);
//	
//	g = makeConstPolynomial_DUZP(1,1);
//	DUZP_t* h;
//	for (int i=0; i<factors2.size; ++i) {
//		h = makeConstPolynomial_DUZP(1,1);
//		fprintf(stderr,"exp = %ld\n",factors2.pairs[i].b);
//		for (int j=0; j<factors2.pairs[i].b; ++j)
//			multiplyPolynomials_DUZP_inp(factors2.pairs[i].a,&h);
//		multiplyPolynomials_DUZP_inp(h,&g);
//		printPoly_DUZP(h,xvar);
//		freePolynomial_DUZP(h);
//	}
//	multiplyByInteger_DUZP_inp(g,cc);
//	
//	fprintf(stderr,"g = ");
//	printPoly_DUZP(g,xvar);
//	fprintf(stderr,"f = ");
//	printPoly_DUZP(f,xvar);
//	g = subtractPolynomials_DUZP(f,g);
//	fprintf(stderr,"f - g = ");
//	printPoly_DUZP(g,xvar);
//	
////	print_factor_data(fd);
//	
//	freePolynomial_DUZP(f);
//	freePolynomial_DUZP(g);

//}
