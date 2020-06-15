

#include <gmpxx.h>

#include "PowerSeries/UPOPS_Hensel.hpp"
#include "PowerSeries/UPOPS_Weierstrass.h"

#if defined(WITH_NTL) && WITH_NTL
#include "IntegerPolynomial/DUZP_Support.h"
#include "IntegerPolynomial/DUZP_Support_Factoring.hpp"
using namespace DUZP::Factoring;

void evalToZeroAndFactor_UPOPS(Upops_t* f, mpq_t** c, int** k, int* n) {
	if (f == NULL) {
		if (n != NULL) {
			*n = 0;
		}
		return;
	}
	if (f->deg < 0) {
		if (n != NULL) {
			*n = 0;
		}
		return;
	}

	AltArr_t* q = makePolynomial_AA(f->deg+1, 1);
	int curQ = 0;
	for (int i = f->deg; i >= 0; --i) {
		if (f->data[i] != NULL && f->data[i]->deg >= 0 && !isZero_AA(f->data[i]->polys[0])) {
			mpq_init(q->elems[curQ].coef);
			mpq_set(q->elems[curQ].coef, f->data[i]->polys[0]->elems[0].coef);
			q->elems[curQ].degs = i;
			++curQ;
		}
	}
	q->size = curQ;

	primitivePart_AA_inp(q);
	DUZP_t* p = convertFromAltArr_DUZP(q);
	freePolynomial_AA(q);

	vec_DUZP_long_t u_factors;
	vec_DUZP_long_t_init2(&u_factors,0);
	mpz_t cc;
	mpz_init(cc);

	factor(&cc,&u_factors,p);

	mpz_clear(cc);

	if (n != NULL) {
		*n = u_factors.size;
	}

	if (k != NULL) {
		int* exps = *k;
		if (exps == NULL) {
			exps = (int*) malloc(sizeof(int)*u_factors.size);
			*k = exps;
		}
		for (int i=0; i<u_factors.size; ++i) {
			exps[i] = (int) u_factors.pairs[i].b;
		}
	}

	if (c != NULL) {
		mpq_t* facts = *c;
		DUZP_t* fact;
		if (facts == NULL) {
			facts = (mpq_t*) malloc(sizeof(mpq_t)*u_factors.size);
			for (int i=0; i<u_factors.size; ++i) {
				mpq_init(facts[i]);
			}
			*c = facts;
		}
		for (int i=0; i<u_factors.size; ++i) {
			fact = u_factors.pairs[i].a;
			if (fact->lt != 1) {
				fprintf(stderr, "evalToZeroAndFactor_UPOPS ERROR Input UPOPS did not factor to linear factors.\n");
				exit(1);
			}
			mpz_neg(mpq_numref(facts[i]), fact->coefs[0]);
			mpz_set(mpq_denref(facts[i]), fact->coefs[1]);
			mpq_canonicalize(facts[i]);
		}
	}

	vec_DUZP_long_t_free(&u_factors);
	freePolynomial_DUZP(p);
}
#endif


/**
 * c_r : the rational number
 * k : the integer
 * Computes lower triangular matrix corresponding to c_r (it is a modified version of the pascal triangle)
 *
 */
mpq_list_t* getPascalTriLowerMat(const mpq_t c_r, int k) {
    mpq_t* data = (mpq_t*) malloc(sizeof(mpq_t)*(k+1)*(k+1));

    mpq_list_t* mpql = (mpq_list_t*) malloc(sizeof(mpq_list_t));
    mpql->size = (k+1)*(k+1);
    mpql->data = data;
    mpql->refCount = 1;

    mpq_t firstEntry;
    mpq_init(firstEntry);
    mpq_t middleEntry;
    mpq_init(middleEntry);

    unsigned long k1 = k+1;
    mpq_t fracUpTo;
    mpq_init(fracUpTo);
	for (unsigned long t = 0; t < k1*k1; ++t) {
	    mpq_init(data[t]);
	}
	for (unsigned i = 0; i < k1; ++i) {
	    for (unsigned long j = 0; j <= i; ++j) {
	        if (j == 0) {
	            mpz_set_ui(mpq_numref(firstEntry), 1ul);
	            mpz_set_ui(mpq_denref(firstEntry), 1ul);
	            mpz_set_ui(mpq_numref(data[i*k1 + j]), 1ul);

	        } else {

	            mpz_set_ui(mpq_numref(middleEntry), i-j+1);
	            mpz_set_ui(mpq_denref(middleEntry), j);
	            mpq_canonicalize(middleEntry);
	            mpq_mul(firstEntry, firstEntry, middleEntry);

	            mpz_pow_ui(mpq_numref(fracUpTo), mpq_numref(c_r), j);
	            mpz_pow_ui(mpq_denref(fracUpTo), mpq_denref(c_r), j);
	            mpq_canonicalize(fracUpTo);

	            mpq_mul(data[i*k1 + j], firstEntry, fracUpTo);
	            mpq_canonicalize(data[i*k1 + j]);
	        }

	    }
	}

	mpq_clear(firstEntry);
	mpq_clear(middleEntry);
	mpq_clear(fracUpTo);
	return mpql;
}




/**
 * d : the degree
 * l : the l-th entry of W (which is W[l])
 * S : the array of rational numbers (row matrix)
 * V : the upops
 * Get the homogPart of degree d for W[l]
 */
Poly_ptr linearCombGen(int d, long long l, mpq_list_t* S, Upops_t*  V) {


	int upopDeg = V->deg;
	int k = upopDeg;
	Poly_ptr ret = NULL;
	int j;
	int k1 = k+1;

	for (int i = l; i <= k; ++i) {
		j = i - l;
		Poly_ptr homog = homogPart_PS(d, V->data[i]);
		AltArr_t* aa = deepCopyPolynomial_AA(homog);
		multiplyByRational_AA_inp(aa, S->data[i*k1+j]);
		ret = addPolynomials_AA_inp(ret, aa, aa == NULL ? 0 : aa->nvar);
		freePolynomial_AA(aa);
	}

	return ret;
}


/**
*
*/
Poly_ptr linearCombGenVoid(int d, void* l, void* param1, void* param2) {
        return linearCombGen(d, (long long) l, (mpq_list_t*) param1, (Upops_t*) param2);
}

/**
 * c_r : the rational number
 * f : the input upop
 * W_out : the output containing W[i]s
 */
Upops_t* taylorShift_UPOPS(Upops_t* f, const mpq_t c_r) {

	int upopDeg = f->deg;
	int k = upopDeg;
	Poly_ptr homog;
	int k1 = k+1;
	int j;

	mpq_list_t* S = getPascalTriLowerMat(c_r, k);

	PowerSeries_t** W = (PowerSeries_t**) malloc(sizeof(PowerSeries_t*) * (k + 1));

	for (int m = 0; m <= k; ++m) {
		W[m] = allocatePowerSeries_PS(f->data[m]->alloc);
	}

	for (long long n = 0; n <= k; ++n) {
		W[n]->genOrder = 3;
		W[n]->genParam1 = (void*) n;
		W[n]->genParam2 = S;
		reserveMPQList_PS(S);
		reserve_UPOPS(f);
		W[n]->genParam3 = f;
		W[n]->paramType1 = PLAIN_DATA;
		W[n]->paramType2 = MPQ_LIST;
		W[n]->paramType3 = UPOPS;
		W[n]->gen.tertiaryGen = &(linearCombGenVoid);
		W[n]->polys[0] = linearCombGen(0, n, S, f);
		W[n]->deg = 0;
	}

	destroyMPQList_PS(S);

	Upops_t* convertToUpop = allocateUnivariatePolynomialOverPowerSeries_UPOPS(0);
	convertToUpop->alloc = k+1;
	convertToUpop->data = W;
	convertToUpop->deg = upopDeg;

	return convertToUpop;
}

#if defined(WITH_NTL) && WITH_NTL
void HenselFactorization_UPOPS(Upops_t* f, Upops_t*** facts_out, int* n) {
	if (facts_out == NULL || n == NULL) {
		return;
	}

	mpq_t* C = NULL;
	evalToZeroAndFactor_UPOPS(f, &C, NULL, n);

	int nfacts = *n;

	Upops_t** facts = (Upops_t**) malloc(sizeof(Upops_t*)*nfacts);


	Upops_t* fstar = f;
	reserve_UPOPS(fstar);
	Upops_t* g;
	Upops_t* p;
	Upops_t* alpha;

	for (int i = 0; i < nfacts; ++i) {
		g = taylorShift_UPOPS(fstar, C[i]);
		weierstrassPreparation_UPOPS(g, &p, &alpha);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(g);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(fstar);

		mpq_neg(C[i], C[i]);
		facts[i] = taylorShift_UPOPS(p, C[i]);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(p);

		if (i+1 < nfacts) {
			fstar = taylorShift_UPOPS(alpha, C[i]);
		}
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha);
	}


	for (int i = 0; i < nfacts; ++i) {
		mpq_clear(C[i]);
	}
	free(C);

	*facts_out = facts;
}
#endif

/**
 * f : the input upop
 * C : the given roots
 * facts : the array of upops
 */
void _FactorizationViaHensel_UPOPS(Upops_t* f, mpq_t* C,  Upops_t** facts, int n) {

	Upops_t* W_out = taylorShift_UPOPS(f, C[0]);
	Upops_t* p_out = NULL;
	Upops_t* alpha_out = NULL;
	weierstrassPreparation_UPOPS(W_out, &p_out, &alpha_out);

	mpq_t negateVal;
	mpq_init(negateVal);
	mpq_neg(negateVal, C[0]);
	Upops_t* W_p_out = taylorShift_UPOPS(p_out, negateVal);
	facts[0] = W_p_out;
	Upops_t* W_alpha_out = taylorShift_UPOPS(alpha_out, negateVal);
	Upops_t* fstar = W_alpha_out;
	if (n > 1) {
		_FactorizationViaHensel_UPOPS(fstar, C+1, facts+1, n-1);
	}
}



// /**
//  * f : the input upop
//  * C : the given roots
//  * n : the number of factors
//  * facts_out : the output
//  */
// void HenselFactorization_UPOPS(Upops_t* f, mpq_t* C, int n, Upops_t*** facts_out) {
//   if(facts_out == NULL) { return; }

//    Upops_t** facts = (Upops_t**) malloc(sizeof(Upops_t*)*n);
//    _FactorizationViaHensel_UPOPS(f, C, facts, n);

//    *facts_out = facts;

// }

void HenselFactorization_UPOPS(Upops_t* f, mpq_t* C, int nfacts, Upops_t*** facts_out) {
	if(facts_out == NULL) { return; }

	Upops_t** facts = (Upops_t**) malloc(sizeof(Upops_t*)*nfacts);

	Upops_t* fstar = f;
	reserve_UPOPS(fstar);
	Upops_t* g;
	Upops_t* p;
	Upops_t* alpha;

	for (int i = 0; i < nfacts; ++i) {
		g = taylorShift_UPOPS(fstar, C[i]);
		weierstrassPreparation_UPOPS(g, &p, &alpha);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(g);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(fstar);

		mpq_neg(C[i], C[i]);
		facts[i] = taylorShift_UPOPS(p, C[i]);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(p);

		if (i+1 < nfacts) {
			fstar = taylorShift_UPOPS(alpha, C[i]);
		}
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha);
	}

	*facts_out = facts;
}
