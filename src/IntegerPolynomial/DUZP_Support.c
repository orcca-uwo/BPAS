
#include "IntegerPolynomial/DUZP_Support.h"
#include <math.h>
#include <time.h>


duspoly_t* convertToPrimeField_DUZP(const DUZP_t* p, const Prime_ptr* pptr) {
	if (p == NULL || p->alloc == 0) {
		return NULL;
	}

	duspoly_t* pp = makePolynomial_spX(p->lt+1);
	polysize_t lt = pp->lt = p->lt;
	mpz_t* coefs = p->coefs;
	elem_t* primeElems = pp->elems;

	unsigned long uPrime = (unsigned long) pptr->prime;

	for (polysize_t i = 0; i <= lt; ++i) {
		primeElems[i] = smallprimefield_convert_in(mpz_fdiv_ui(coefs[i], uPrime), pptr);
	}

	pp->lt = lt;

	return pp;
}

void convertToPrimeField_DUZP_inp(const DUZP_t* p, const Prime_ptr* pptr, duspoly_t** pp) {
	if (pp == NULL) {
		return;
	}

	if (p == NULL || p->alloc == 0) {
		freePolynomial_spX (&*pp);
		*pp = NULL;
		return;
	}

	if (*pp == NULL) {
		*pp = makePolynomial_spX(p->lt+1);
	} else if ((*pp)->alloc < p->lt+1) {
    	reallocatePolynomial_spX(pp, p->lt+1);
	}
	polysize_t lt = (*pp)->lt = p->lt;
	mpz_t* coefs = p->coefs;
	elem_t* primeElems = (*pp)->elems;

	unsigned long uPrime = (unsigned long) pptr->prime;

	for (polysize_t i = 0; i <= lt; ++i) {
		primeElems[i] = smallprimefield_convert_in(mpz_fdiv_ui(coefs[i], uPrime), pptr);
	}

	(*pp)->lt = lt;
}

void applyModulo_DUZP_inp(DUZP_t* p, const mpz_t mod) {
	if (p == NULL || p->alloc == 0) {
		return;
	}

	polysize_t lt = p->lt;
	int foundLt = 0;
	for (polysize_t i = lt; i >= 0; --i) {
		mpz_mod(p->coefs[i], p->coefs[i], mod);
		if (!foundLt && mpz_sgn(p->coefs[i]) == 0) {
			mpz_clear(p->coefs[i]);
			--lt;
		} else {
			foundLt = 1;
		}
	}
	if (lt < 0) {
		p->lt = 0l;
		mpz_init(p->coefs[0]);
	} else {
		p->lt = lt;
	}
}

void applyModuloSymmetric_DUZP_inp(DUZP_t* p, const mpz_t mod) {
	if (p == NULL || p->alloc == 0) {
		return;
	}

	polysize_t lt = p->lt;
	mpz_t halfMod;
	mpz_init_set(halfMod, mod);
	mpz_div_ui(halfMod, halfMod, 2ul);
	int foundLt = 0;
	for (polysize_t i = lt; i >= 0; --i) {
		// fprintf(stderr, "moding poly %p at %ld\n",p, i);
		mpz_mod(p->coefs[i], p->coefs[i], mod);
		if(mpz_cmp(p->coefs[i], halfMod) > 0) {
			mpz_sub(p->coefs[i], p->coefs[i], mod);
		}

		if (!foundLt && mpz_sgn(p->coefs[i]) == 0) {
			mpz_clear(p->coefs[i]);
			--lt;
		} else {
			foundLt = 1;
		}
	}
	if (lt < 0) {
		p->lt = 0l;
		mpz_init(p->coefs[0]);
	} else {
		p->lt = lt;
	}

	mpz_clear(halfMod);
}

void evaluate_DUZP(const DUZP_t* p, const mpz_t c, mpz_t val) {
	if (c == NULL || val == NULL) {
		return;
	}
	if (isZero_DUZP(p)) {
		mpz_set_ui(val, 0ul);
		return;
	}
	int n = p->lt;
	mpz_set(val,p->coefs[n]);
	for (int i = n-1; i >= 0; --i) {
		mpz_mul(val,val,c);
		mpz_add(val,val,p->coefs[i]);
	}
}

DUZP_t* deepCopyPolynomial_DUZP(const DUZP_t* p) {
	if (p == NULL) {
		return NULL;
	}

	DUZP_t* ret = makePolynomial_DUZP(p->lt + 1);
	polysize_t lt = ret->lt = p->lt;
	mpz_t* retcoefs = ret->coefs;
	mpz_t* coefs = p->coefs;
	mpz_set(retcoefs[0], coefs[0]);
	for (polysize_t i = 1; i <= lt; ++i) {
		mpz_init_set(retcoefs[i], coefs[i]);
	}

	return ret;
}

void deepCopyPolynomial_DUZP_inp(const DUZP_t* a, DUZP_t** bb) {
	if (bb == NULL) {
		return;
	}

	if (a == NULL) {
		if (*bb != NULL) {
			freePolynomial_DUZP(*bb);
			*bb = NULL;
		}
		return;
	}

	DUZP_t* b = *bb;
	if (b == NULL) {
		b = makePolynomial_DUZP(a->lt + 1);
		*bb = b;
	} else if (b->alloc <= a->lt) {
		resizePolynomial_DUZP(b, a->lt + 1);
	}

	int i;
	for (i = b->lt + 1; i <= a->lt; ++i) {
		mpz_init(b->coefs[i]);
	}
	for (i = a->lt + 1; i <= b->lt; ++i) {
		mpz_clear(b->coefs[i]);
	}
	for (i = 0; i <= a->lt; ++i) {
		mpz_set(b->coefs[i], a->coefs[i]);
	}
	b->lt = a->lt;
}

DUZP_t* deepCopyPolynomial_DUZPFromspX(const duspoly_t* p, const Prime_ptr* pptr) {
	if (p == NULL || p->alloc == 0) {
		return NULL;
	}

	DUZP_t* pp = makePolynomial_DUZP(p->alloc);
	polysize_t lt = pp->lt = p->lt;
	mpz_t* coefs = pp->coefs;
	elem_t* primeElems = p->elems;
	for (polysize_t i = 0; i <= lt; ++i) {
		mpz_init(coefs[i]);
		mpz_set_si(coefs[i], smallprimefield_convert_out(primeElems[i], pptr));
	}

	return pp;
}

DUZP_t* deepCopyPolynomial_DUZPFromspXSymmetric(const duspoly_t* p, const Prime_ptr* pptr) {
	if (p == NULL || p->alloc == 0) {
		return NULL;
	}

	DUZP_t* pp = makePolynomial_DUZP(p->lt + 1);
	polysize_t lt = pp->lt = p->lt;
	mpz_t* coefs = pp->coefs;
	elem_t* primeElems = p->elems;
	long long int halfP = (pptr->prime - 1) >> 1;
	long long int tmp;
	tmp = smallprimefield_convert_out(primeElems[0], pptr);
	if (tmp > halfP) {
		tmp -= pptr->prime;
	}
	mpz_set_si(coefs[0], tmp);
	for (polysize_t i = 1; i <= lt; ++i) {
		tmp = smallprimefield_convert_out(primeElems[i], pptr);
		if (tmp > halfP) {
			tmp -= pptr->prime;
		}
		mpz_init_set_si(coefs[i], tmp);
	}

	return pp;
}

int isEqual_DUZP(const DUZP_t* a, const DUZP_t* b) {
	if (a == b) {
		return 1;
	}
	if (a == NULL || b == NULL) {
		//from previous check, both a and b cannot be null
		return 0;
	}

	polysize_t lt_a = a->lt, lt_b = b->lt;

	if (lt_a != lt_b) {
		return 0;
	}

	mpz_t* acoefs = a->coefs, * bcoefs = b->coefs;
	for (lt_b = 0; lt_b <= lt_a; ++lt_b) {
		if (mpz_cmp(acoefs[lt_b], bcoefs[lt_b]) != 0) {
			return 0;
		}
	}

	return 1;
}


void printPoly_DUZP(const DUZP_t* poly, char* sym) {
	if (poly == NULL || poly->coefs == NULL || (poly->lt == 0 && mpz_sgn(poly->coefs[0]) == 0)) {
		fprintf(stderr, "0\n");
		return;
	}

	polysize_t lt = poly->lt;
	mpz_t* coefs = poly->coefs;
	int first = 1;
	mpz_t tmp;
	mpz_init(tmp);
	for (polysize_t i = 0; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) == 0) {
			continue;
		}
		if (first) {
			gmp_fprintf(stderr, "%Zd*%s^%ld", coefs[i], sym, i);
			first = 0;
			continue;
		}
		if (mpz_sgn(coefs[i]) > 0) {
			fprintf(stderr, " + ");
			gmp_fprintf(stderr, "%Zd*%s^%ld", coefs[i], sym, i);
		} else {
			fprintf(stderr, " - ");
			mpz_abs(tmp, coefs[i]);
			gmp_fprintf(stderr, "%Zd*%s^%ld", tmp, sym, i);
		}
	}
	fprintf(stderr, "\n");

	mpz_clear(tmp);
}

void printPoly_DUZP_maple(DUZP_t* poly, char* sym, char* maplevar) {
	if (poly == NULL) {
		fprintf(stderr, "%s := 0;\n", maplevar);
		return;
	}

	fprintf(stderr, "%s := ", maplevar);
	polysize_t lt = poly->lt;
	mpz_t* coefs = poly->coefs;
	int first = 1;
	mpz_t tmp;
	mpz_init(tmp);
	for (polysize_t i = 0; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) == 0) {
			continue;
		}
		if (first) {
			gmp_fprintf(stderr, "%Zd*%s^%ld", coefs[i], sym, i);
			first = 0;
			continue;
		}
		if (mpz_sgn(coefs[i]) > 0) {
			fprintf(stderr, " + ");
			gmp_fprintf(stderr, "%Zd*%s^%ld", coefs[i], sym, i);
		} else {
			fprintf(stderr, " - ");
			mpz_abs(tmp, coefs[i]);
			gmp_fprintf(stderr, "%Zd*%s^%ld", tmp, sym, i);
		}
	}
	fprintf(stderr, ";\n");

	mpz_clear(tmp);
}

void printPolyToFile_DUZP(FILE* fp, DUZP_t* poly, char* sym) {
	if (poly == NULL || (poly->lt == 0 && mpz_sgn(poly->coefs[0]) == 0)) {
		fprintf(fp, "0\n");
		return;
	}

	polysize_t lt = poly->lt;
	mpz_t* coefs = poly->coefs;
	int first = 1;
	mpz_t tmp;
	mpz_init(tmp);
	for (polysize_t i = 0; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) == 0) {
			continue;
		}
		if (first) {
			gmp_fprintf(fp, "%Zd*%s^%ld", coefs[i], sym, i);
			first = 0;
			continue;
		}
		if (mpz_sgn(coefs[i]) > 0) {
			fprintf(fp, " + ");
			gmp_fprintf(fp, "%Zd*%s^%ld", coefs[i], sym, i);
		} else {
			fprintf(fp, " - ");
			mpz_abs(tmp, coefs[i]);
			gmp_fprintf(fp, "%Zd*%s^%ld", tmp, sym, i);
		}
	}
	fprintf(fp, "\n");

	mpz_clear(tmp);
}

DUZP_t* buildRandomPoly_DUZP(polysize_t maxDeg, int coefBoundBits, float sparsity, int includeNeg) {
	polysize_t nterms = ceil((maxDeg + 1) * (1.0 - sparsity));

	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL); //1566936954;
		// t = 1579622689; // fail subresultant
		// t = 1580176795; // fail modularsubresultant-hgcd
/*		time_t t = 173;*/
/*		t = 1568063634;*/
/*		t = 1568402466; // hanging at modular factorization stage, searching for 2 primes*/
/*		t = 1568403655; // same*/
/*		t = 1568403775; // hanging at modular factorization stage, searching for 1 prime*/
/*		t = 1568403844; // slow at modular factorization stage, but finishest, for 1 prime*/
/*		t = 1568403947; // same*/
/*		t = 1568403986; // same*/
/*		t = 1568403414; // not divisible yo*/
/*		t = 1569260375; // irreducible deg 20 monic; not divisible yo*/
/*		t = 1570045809; // not divisible yo*/
/*		t = 1570469451; // vH algo failure for recombine_test() sizes[4] = {2,3,4,5,6} LM bound*/
/*		t = 1570469634;*/
		fprintf(stderr,"t = %ld\n",t);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

	DUZP_t* p = makePolynomial_DUZP(maxDeg+1);
	mpz_t* coefs = p->coefs;
	p->lt = maxDeg;

	if (sparsity > 0.9999) {
		mpz_urandomb(coefs[0], R_STATE, coefBoundBits);
		if (includeNeg && rand() % 2) {
			mpz_neg(coefs[0], coefs[0]);
		}

		for (polysize_t i = 1; i < nterms; ++i) {
			mpz_init(coefs[i]);
			mpz_urandomb(coefs[i], R_STATE, coefBoundBits);
			if (includeNeg && rand() % 2) {
				mpz_neg(coefs[i], coefs[i]);
			}
		}
		while(mpz_sgn(coefs[maxDeg]) == 0) {
			mpz_urandomb(coefs[maxDeg], R_STATE, coefBoundBits);
		}
		return p;
	}

	for (polysize_t i = 1; i < maxDeg; ++i) {
		mpz_init(coefs[i]);
	}
	mpz_init(coefs[maxDeg]);

	while(mpz_sgn(coefs[maxDeg]) == 0) {
		mpz_urandomb(coefs[maxDeg], R_STATE, coefBoundBits);
	}
	if (includeNeg && rand() % 2) {
		mpz_neg(coefs[maxDeg], coefs[maxDeg]);
	}
	--nterms;

	polysize_t idx;
	for (polysize_t i = 0; i < nterms; ++i) {
		idx = rand() % maxDeg;
		while(mpz_sgn(coefs[idx]) == 0) {
			mpz_urandomb(coefs[idx], R_STATE, coefBoundBits);
		}
		if (includeNeg && rand() % 2) {
			mpz_neg(coefs[idx], coefs[idx]);
		}
	}
	return p;
}

void content_DUZP(const DUZP_t* p, mpz_t c) {
	if (p == NULL || p->alloc == 0){
		mpz_set_ui(c, 0ul);
		return;
	}

	mpz_t* coefs = p->coefs;
	polysize_t lt = p->lt;

	polysize_t i = 0;
	for (; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) != 0) {
			mpz_set(c, coefs[i]);
			break;
		}
	}

	for (; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) != 0) {
			mpz_gcd(c, c, coefs[i]);
			if (mpz_cmp_si(c, 1l) == 0) {
				if (mpz_sgn(coefs[lt]) < 0) {
					mpz_neg(c,c);
				}
				return;
			}
		}
	}

	if (mpz_sgn(coefs[lt]) < 0 && mpz_sgn(c) > 0) {
		mpz_neg(c, c);
	}

	return;
}

DUZP_t* primitivePart_DUZP(const DUZP_t* p) {

	DUZP_t* pp = deepCopyPolynomial_DUZP(p);
	mpz_t cont;
	mpz_init(cont);
	content_DUZP(pp, cont);
	if (mpz_cmp_si(cont, 1l) == 0) {
		mpz_clear(cont);
		return pp;
	}

	divideByIntegerExact_DUZP_inp(pp, cont);
	mpz_clear(cont);
	return pp;
}

void primitivePart_DUZP_inp(DUZP_t* p) {
	mpz_t cont;
	mpz_init(cont);
	content_DUZP(p, cont);
	if (mpz_cmp_si(cont, 1l) == 0) {
		mpz_clear(cont);
		return;
	}

	divideByIntegerExact_DUZP_inp(p, cont);
	mpz_clear(cont);
}

void primitivePartAndContent_DUZP_inp(DUZP_t* p, mpz_t cont) {
	content_DUZP(p, cont);
	if (mpz_cmp_si(cont, 1l) == 0) {
		return;
	}

	divideByIntegerExact_DUZP_inp(p, cont);
}


DUZP_t* primitivePartAndContent_DUZP(const DUZP_t* p, mpz_t cont) {
	DUZP_t* pp = deepCopyPolynomial_DUZP(p);
	content_DUZP(pp, cont);
	if (mpz_cmp_si(cont, 1l) == 0) {
		return pp;
	}

	divideByIntegerExact_DUZP_inp(pp, cont);
	return pp;
}





/************************
 * Arithmetic Functions
 ************************/

DUZP_t* subtractDUSP_DUZP(DUZP_t* a, duspoly_t* b, const Prime_ptr* Pptr) {
	if (isZero_DUZP(a)) {
		DUZP_t* ret = convertFromDUSP_DUZP(b, Pptr);
		negatePolynomial_DUZP_inp(ret);
		return ret;
	}

	if (isZero_spX(b)) {
		return deepCopyPolynomial_DUZP(a);
	}

	polysize_t diffSize = MAX_spX(a->lt, b->lt) + 1;
	DUZP_t* d = makePolynomial_DUZP(diffSize);

	polysize_t minSize = MIN_spX(a->lt, b->lt) + 1;

	mpz_sub_ui(d->coefs[0], a->coefs[0], (unsigned long long) smallprimefield_convert_out(b->elems[0], Pptr));
	for (polysize_t i = 1; i < minSize; ++i) {
		mpz_init(d->coefs[i]);
		mpz_sub_ui(d->coefs[i], a->coefs[i], (unsigned long long) smallprimefield_convert_out(b->elems[i], Pptr));
	}

	//if deg a  < deg b
	if (minSize == a->lt + 1) {
		for (polysize_t i = minSize; i <= b->lt; ++i) {
			mpz_init(d->coefs[i]);
			mpz_set_si(d->coefs[i], smallprimefield_convert_out(b->elems[i], Pptr) * -1);
		}
	} else {
		for (polysize_t i = minSize; i <= a->lt; ++i) {
			mpz_init_set(d->coefs[i], a->coefs[i]);
		}
	}

	d->lt = diffSize - 1;

	return d;
}

DUZP_t* addPolynomials_DUZP(DUZP_t* a, DUZP_t* b) {
	if (isZero_DUZP(a)) {
		return deepCopyPolynomial_DUZP(b);
	}

	if (isZero_DUZP(b)) {
		return deepCopyPolynomial_DUZP(a);
	}

	polysize_t diffSize = MAX_spX(a->lt, b->lt) + 1;
	DUZP_t* d = makePolynomial_DUZP(diffSize);

	polysize_t minSize = MIN_spX(a->lt, b->lt) + 1;

	mpz_add(d->coefs[0], a->coefs[0], b->coefs[0]);
	for (polysize_t i = 1; i < minSize; ++i) {
		mpz_init(d->coefs[i]);
		mpz_add(d->coefs[i], a->coefs[i], b->coefs[i]);
	}

	//if deg a  < deg b
	if (minSize == a->lt + 1) {
		for (polysize_t i = minSize; i <= b->lt; ++i) {
			mpz_init_set(d->coefs[i], b->coefs[i]);
		}
	} else {
		for (polysize_t i = minSize; i <= a->lt; ++i) {
			mpz_init_set(d->coefs[i], a->coefs[i]);
		}
	}

	d->lt = diffSize - 1;
	for (int i = d->lt; i >= 1; --i) {
    	if (mpz_sgn(d->coefs[i]) == 0) {
    		mpz_clear(d->coefs[i]);
    		--(d->lt);
    	} else {
    		break;
    	}
    }


	return d;
}

void addPolynomials_DUZP_inp(DUZP_t** aa, DUZP_t* b) {
	if (aa == NULL) {
		return;
	}
	if (isZero_DUZP(b)) {
		return;
	}

	DUZP_t* a = *aa;

	if (isZero_DUZP(a)) {
		if (a == NULL) {
			a = makePolynomial_DUZP(b->lt + 1);
		} else if (a->alloc <= b->lt) {
			resizePolynomial_DUZP(a, b->lt + 1);
		}
		polysize_t i;
		for (i = 0; i <= a->lt; ++i) {
			mpz_set(a->coefs[i], b->coefs[i]);
		}
		for ( ; i <= b->lt; ++i) {
			mpz_init_set(a->coefs[i], b->coefs[i]);
		}
		a->lt = b->lt;
		return;
	}

	polysize_t a_lt = a->lt;
	polysize_t b_lt = b->lt;
	polysize_t i;
	if (a_lt < b_lt) {
		if (a->alloc <= b_lt) {
			resizePolynomial_DUZP(a, b->lt + 1);
		}
		for (i = 0; i <= a_lt; ++i) {
			mpz_add(a->coefs[i], a->coefs[i], b->coefs[i]);
		}
		for (; i <= b_lt; ++i) {
			mpz_init_set(a->coefs[i], b->coefs[i]);
		}
		a->lt = b->lt;
	} else {
		for (i = 0; i <= b_lt; ++i) {
			mpz_add(a->coefs[i], a->coefs[i], b->coefs[i]);
		}
	}

	for (i = a->lt; i >= 1; --i) {
		if (mpz_sgn(a->coefs[i]) == 0) {
			mpz_clear(a->coefs[i]);
			--(a->lt);
		} else {
			break;
		}
	}

}

void addPolynomialsPreAlloc_DUZP(DUZP_t* a, DUZP_t* b, DUZP_t** sump) {
	if (sump == NULL) {
		return;
	}

	if (isZero_DUZP(a) && isZero_DUZP(b)) {
		freePolynomial_DUZP(*sump);
		*sump = NULL;
	}

	DUZP_t* sum = *sump;

	if (isZero_DUZP(a)) {
		if (sum == NULL) {
			sum = makePolynomial_DUZP(b->lt + 1);
			*sump = sum;
		} else if (sum->alloc <= b->lt) {
			resizePolynomial_DUZP(sum, b->lt + 1);
		}
		polysize_t i;
		polysize_t slt = sum->lt;
		sum->lt = b->lt;
		if (slt > b->lt) {
			for (i = 0; i <= b->lt; ++i) {
				mpz_set(sum->coefs[i], b->coefs[i]);
			}
			for (; i <= slt; ++i) {
				mpz_clear(sum->coefs[i]);
			}
		} else {
			for (i = 0; i <= slt; ++i) {
				mpz_set(sum->coefs[i], b->coefs[i]);
			}
			for (; i <= b->lt; ++i) {
				mpz_init_set(sum->coefs[i], b->coefs[i]);
			}
		}
		return;
	}

	if (isZero_DUZP(b)) {
		if (sum == NULL) {
			sum = makePolynomial_DUZP(a->lt + 1);
			*sump = sum;
		} else if (sum->alloc <= a->lt) {
			resizePolynomial_DUZP(sum, a->lt + 1);
		}
		polysize_t i;
		polysize_t slt = sum->lt;
		sum->lt = a->lt;
		if (slt > a->lt) {
			for (i = 0; i <= a->lt; ++i) {
				mpz_set(sum->coefs[i], a->coefs[i]);
			}
			for (; i <= slt; ++i) {
				mpz_clear(sum->coefs[i]);
			}
		} else {
			for (i = 0; i <= slt; ++i) {
				mpz_set(sum->coefs[i], a->coefs[i]);
			}
			for (; i <= a->lt; ++i) {
				mpz_init_set(sum->coefs[i], a->coefs[i]);
			}
		}
		return;
	}


	polysize_t sumSize = MAX_spX(a->lt, b->lt) + 1;
	if (sum == NULL) {
		sum = makePolynomial_DUZP(sumSize);
		*sump = sum;
	} else if (sum->alloc <= sumSize) {
		resizePolynomial_DUZP(sum, sumSize);
	}

	polysize_t minSize = MIN_spX(a->lt, b->lt) + 1;
	polysize_t slt = sum->lt;

	if (slt < minSize) {
		polysize_t i;
		for (i = 0; i <= sum->lt; ++i) {
			mpz_add(sum->coefs[i], a->coefs[i], b->coefs[i]);
		}
		for (; i < minSize; ++i) {
			mpz_init(sum->coefs[i]);
			mpz_add(sum->coefs[i], a->coefs[i], b->coefs[i]);
		}
		//if deg a  < deg b
		if (minSize == a->lt + 1) {
			for (i = minSize; i <= b->lt; ++i) {
				mpz_init_set(sum->coefs[i], b->coefs[i]);
			}
			sum->lt = b->lt;
		} else {
			for (i = minSize; i <= a->lt; ++i) {
				mpz_init_set(sum->coefs[i], a->coefs[i]);
			}
			sum->lt = a->lt;
		}
	} else {
		polysize_t i;
		for (i = 0; i < minSize; ++i) {
			mpz_add(sum->coefs[i], a->coefs[i], b->coefs[i]);
		}

		if (minSize == a->lt + 1) {
			if (slt < b->lt) {
				for ( ; i <= slt; ++i) {
					mpz_set(sum->coefs[i], b->coefs[i]);
				}
				for ( ; i <= b->lt; ++i) {
					mpz_init_set(sum->coefs[i], b->coefs[i]);
				}
			} else {
				for ( ; i <= b->lt; ++i) {
					mpz_set(sum->coefs[i], b->coefs[i]);
				}
			}
		} else {
			if (slt < a->lt) {
				for ( ; i <= slt; ++i) {
					mpz_set(sum->coefs[i], a->coefs[i]);
				}
				for ( ; i <= a->lt; ++i) {
					mpz_init_set(sum->coefs[i], a->coefs[i]);
				}
			} else {
				for ( ; i <= slt; ++i) {
					mpz_set(sum->coefs[i], a->coefs[i]);
				}
			}
		}

		//resuing min var to be MAX
		minSize = MAX_spX(a->lt, b->lt);
		for ( i = minSize + 1 ; i <= slt; ++i) {
			mpz_clear(sum->coefs[i]);
		}
		sum->lt = minSize;
	}
}


DUZP_t* subtractPolynomials_DUZP(DUZP_t* a, DUZP_t* b) {
	if (isZero_DUZP(a)) {
		DUZP_t* ret = deepCopyPolynomial_DUZP(b);
		negatePolynomial_DUZP_inp(ret);
		return ret;
	}

	if (isZero_DUZP(b)) {
		return deepCopyPolynomial_DUZP(a);
	}

	polysize_t diffSize = MAX_spX(a->lt, b->lt) + 1;
	DUZP_t* d = makePolynomial_DUZP(diffSize);

	polysize_t minSize = MIN_spX(a->lt, b->lt) + 1;

	mpz_sub(d->coefs[0], a->coefs[0], b->coefs[0]);
	for (polysize_t i = 1; i < minSize; ++i) {
		mpz_init(d->coefs[i]);
		mpz_sub(d->coefs[i], a->coefs[i], b->coefs[i]);
	}

	//if deg a  < deg b
	if (minSize == a->lt + 1) {
		for (polysize_t i = minSize; i <= b->lt; ++i) {
			mpz_init(d->coefs[i]);
			mpz_neg(d->coefs[i], b->coefs[i]);
		}
	} else {
		for (polysize_t i = minSize; i <= a->lt; ++i) {
			mpz_init_set(d->coefs[i], a->coefs[i]);
		}
	}

	d->lt = diffSize - 1;
	for (int i = d->lt; i >= 1; --i) {
    	if (mpz_sgn(d->coefs[i]) == 0) {
    		mpz_clear(d->coefs[i]);
    		--(d->lt);
    	} else {
    		break;
    	}
    }

	return d;
}

void subtractPolynomials_DUZP_inp(DUZP_t** aa, DUZP_t* b) {
	if (aa == NULL) {
		return;
	}
	DUZP_t* a = *aa;

	if (isZero_DUZP(a)) {
		if (a->alloc <= b->lt) {
			resizePolynomial_DUZP(a, b->lt + 1);
		}
		polysize_t i;
		for (i = 0; i <= a->lt; ++i) {
			mpz_neg(a->coefs[i], b->coefs[i]);
		}
		for ( ; i <= b->lt; ++i) {
			mpz_init_set(a->coefs[i], b->coefs[i]);
			mpz_neg(a->coefs[i], a->coefs[i]);
		}
		return;
	}

	if (isZero_DUZP(b)) {
		return;
	}

	polysize_t a_lt = a->lt;
	polysize_t b_lt = b->lt;
	polysize_t i;
	if (a_lt < b_lt) {
		if (a->alloc <= b_lt) {
			resizePolynomial_DUZP(a, b->lt + 1);
		}
		for (i = 0; i <= a_lt; ++i) {
			mpz_sub(a->coefs[i], a->coefs[i], b->coefs[i]);
		}
		for (; i <= b_lt; ++i) {
			mpz_init_set(a->coefs[i], b->coefs[i]);
			mpz_neg(a->coefs[i], a->coefs[i]);
		}
		a->lt = b->lt;
	} else {
		for (i = 0; i <= b_lt; ++i) {
			mpz_sub(a->coefs[i], a->coefs[i], b->coefs[i]);
		}
	}
	for (int i = a->lt; i >= 1; --i) {
    	if (mpz_sgn(a->coefs[i]) == 0) {
    		mpz_clear(a->coefs[i]);
    		--(a->lt);
    	} else {
    		break;
    	}
    }

}

void subtractPolynomials_DUZP_inpRHS(const DUZP_t* a, DUZP_t** b) {
	if (b == NULL) {
		return;
	}

	if (isZero_DUZP(*b)) {
		// freePolynomial_DUZP(*b);
		*b = deepCopyPolynomial_DUZP(a);
		return;
	}

	polysize_t minLT = MIN_spX(a->lt,(*b)->lt);

	mpz_t* belems = (*b)->coefs;
	mpz_t* aelems = a->coefs;
	for (polysize_t i = 0; i <= minLT; ++i) {
		mpz_sub(belems[i], a->coefs[i], belems[i]);
	}

	if (a->lt > minLT) {
		if (a->lt + 1 > (*b)->alloc) {
			resizePolynomial_DUZP(*b, a->lt + 1);
			belems = (*b)->coefs;
		}
		for (polysize_t i = minLT + 1; i <= a->lt; ++i) {
			mpz_init_set(belems[i], aelems[i]);
		}
		(*b)->lt = a->lt;
	} else if ((*b)->lt > minLT) {
		polysize_t bLT = (*b)->lt;
		for (polysize_t i = minLT + 1; i <= bLT; ++i) {
			mpz_neg(belems[i], belems[i]);
		}
	}

	for (int i = (*b)->lt; i >= 0; --i) {
		if(mpz_sgn(belems[i]) == 0) {
			mpz_clear(belems[i]);
		} else {
			(*b)->lt = i;
			break;
		}

		if (i == 0) {
			(*b)->lt = -1;
			freePolynomial_DUZP(*b);
			*b = NULL;
		}
	}
}

duspoly_t* multiplyAndModDUSP_DUZP(DUZP_t* a, duspoly_t* b, const Prime_ptr* Pptr) {

	DUZP_t* bb = convertFromDUSP_DUZP(b, Pptr);
	DUZP_t* c = multiplyPolynomials_DUZP(a, bb);
	freePolynomial_DUZP(bb);

	duspoly_t* ret = convertToPrimeField_DUZP(c, Pptr);
	freePolynomial_DUZP(c);
	return ret;
}

void multiplyAndModDUSP_DUZP_inp(DUZP_t* a, duspoly_t** b, const Prime_ptr* Pptr) {

	DUZP_t* bb = convertFromDUSP_DUZP(*b, Pptr);
	DUZP_t* c = multiplyPolynomials_DUZP(a, bb);

	convertToPrimeField_DUZP_inp(c, Pptr, b);
}



DUZP_t* multiplyPolynomials_DUZP(const DUZP_t* a, const DUZP_t* b) {
    // poly * 0 = 0
    if (a == NULL || b == NULL) {
		return NULL;
    }

    if (a->lt == 0) {
    	return multiplyByInteger_DUZP(b, a->coefs[0]);
    }
    if (b->lt == 0) {
    	return multiplyByInteger_DUZP(a, b->coefs[0]);
    }

    polysize_t lt_a = a->lt;
    polysize_t lt_b = b->lt;
    polysize_t lt_c = lt_a + lt_b;

    DUZP_t* c = makePolynomial_DUZP(lt_c+1);
    c->lt = lt_c;

    register polysize_t i, j;
    register polysize_t jmin = 0;
    register polysize_t jmax = 0;

    mpz_t* a_coefs = a->coefs;
    mpz_t* b_coefs = b->coefs;
    mpz_t* c_coefs = c->coefs;

    // poly * poly =
    mpz_addmul(c_coefs[0], a_coefs[0], b_coefs[0]);
    for (i = 1; i <= lt_c; i++) {
		mpz_init(c_coefs[i]);
		jmin = MIN_spX(i, lt_a);
		jmax = MAX_spX(i - lt_b, 0);
		for (j = jmax; j <= jmin; j++) {
		    // c[i] += a[j]*b[i-j];
		    mpz_addmul(c_coefs[i], a_coefs[j], b_coefs[i-j]);
		}
    }

    return c;
}

// DUZP_t* multiplyAllButOnePolynomials_DUZP(DUZP_t** polys, int n, int excludeIdx) {
// 	if (n == 0 || polys == NULL) {
// 		return NULL;
// 	}
// 	if (n < 2) {
// 		return deepCopyPolynomial_DUZP(polys[0]);
// 	}

// 	if (excludeIdx == 0) {
// 		DUZP_t* prod = multiplyByInteger_DUZP(polys[1], polys[2]);
// 		DUZP_t* tmpProd;
// 		for (int i = 2; i < n; ++i) {
// 			tmpProd = multiplyPolynomials_DUZP(prod, polys[i]);
// 			freePolynomial_DUZP(prod);
// 			prod = tmpProd;
// 		}
// 		return prod;
// 	} else if (excludeIdx == 1) {
// 		DUZP_t* prod = multiplyPolynomials_DUZP(polys[0], polys[2]);
// 		DUZP_t* tmpProd;
// 		for (int i = 2; i < n; ++i) {
// 			tmpProd = multiplyPolynomials_DUZP(prod, polys[i]);
// 			freePolynomial_DUZP(prod);
// 			prod = tmpProd;
// 		}
// 		return prod;
// 	} else {
// 		DUZP_t* prod = multiplyPolynomials_DUZP(polys[0], polys[1]);
// 		DUZP_t* tmpProd;
// 		for (int i = 2; i < n; ++i) {
// 			if (i == excludeIdx) {
// 				continue;
// 			}
// 			tmpProd = multiplyPolynomials_DUZP(prod, polys[i]);
// 			freePolynomial_DUZP(prod);
// 			prod = tmpProd;
// 		}
// 		return prod;
// 	}
// }

DUZP_t* multiplyAllButOnePolynomials_DUZP(DUZP_t const*const* polys, int n, int idx) {
	if (n == 0 || polys == NULL) {
		return NULL;
	}
	if (n <= 2) {
		if (idx == 0) {
			return deepCopyPolynomial_DUZP(polys[1]);
		}
		return deepCopyPolynomial_DUZP(polys[0]);
	}

//1: Naive iteartive way
	polysize_t totalDeg = 0;
	for (int i = 0; i < n; ++i) {
		totalDeg += polys[i]->lt;
	}

	DUZP_t* prod = makePolynomial_DUZP(totalDeg + 1);
	if (idx < 2) {
		if (idx == 0) {
			multiplyPolynomialsPreAlloc_DUZP(polys[1], polys[2], &prod);
		} else {
			multiplyPolynomialsPreAlloc_DUZP(polys[0], polys[2], &prod);
		}
		for (int i = 3; i < n; ++i) {
			multiplyPolynomials_DUZP_inp(polys[i], &prod);
		}
	} else {
		multiplyPolynomialsPreAlloc_DUZP(polys[0], polys[1], &prod);
		for (int i = 2; i < n; ++i) {
			if ( i == idx) {
				continue;
			}
			multiplyPolynomials_DUZP_inp(polys[i], &prod);
		}
	}

	// DUZP_t* tmpProd, *tmp;
	// if (idx == 0 ){
	// 	tmpProd = multiplyPolynomials_DUZP(polys[1], polys[2]);
	// 	for (int i = 3; i < n; ++i) {
	// 		tmp = multiplyPolynomials_DUZP(tmpProd, polys[i]);
	// 		tmpProd = tmp;
	// 	}
	// } else if (idx == 1) {
	// 	tmpProd = multiplyPolynomials_DUZP(polys[0], polys[2]);
	// 	for (int i = 3; i < n; ++i) {
	// 		tmp = multiplyPolynomials_DUZP(tmpProd, polys[i]);
	// 		tmpProd = tmp;
	// 	}
	// } else {
	// 	tmpProd = multiplyPolynomials_DUZP(polys[0], polys[1]);
	// 	for (int i = 2; i < n; ++i) {
	// 		if (i ==  idx) {
	// 			continue;
	// 		}
	// 		tmp = multiplyPolynomials_DUZP(tmpProd, polys[i]);
	// 		tmpProd = tmp;
	// 	}
	// }

	// if(!isEqual_DUZP(tmpProd, prod)) {
	// 	fprintf(stderr, "\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" );
	// 	fprintf(stderr, "prod:= \n");
	// 	printPoly_DUZP(prod, "x");
	// 	fprintf(stderr, "tmpProd:= \n");
	// 	printPoly_DUZP(tmpProd, "x");
	// }

	return prod;
}

void multiplyAllButOnePolynomialsPreAlloc_DUZP(DUZP_t const*const* polys, int n, int idx, DUZP_t** pprod) {
	if (n == 0 || polys == NULL || pprod == NULL) {
		return;
	}
	DUZP_t* prod = *pprod;
	if (n <= 2) {
		int k = idx == 0 ? 1 : 0;
		if (prod->alloc <= polys[k]->lt) {
			resizePolynomial_DUZP(prod, polys[k]->lt);
		}
		polysize_t i;
		polysize_t plt = prod->lt;
		prod->lt = polys[k]->lt;
		if (plt > polys[k]->lt) {
			for (i = 0; i <= polys[k]->lt; ++i) {
				mpz_set(prod->coefs[i], polys[k]->coefs[i]);
			}
			for (; i <= plt; ++i) {
				mpz_clear(prod->coefs[i]);
			}
		} else {
			for (i = 0; i <= plt; ++i) {
				mpz_set(prod->coefs[i], polys[k]->coefs[i]);
			}
			for (; i <= polys[k]->lt; ++i) {
				mpz_init_set(prod->coefs[i], polys[k]->coefs[i]);
			}
		}
		return;
	}

//1: Naive iteartive way
	polysize_t totalDeg = 0;
	for (int i = 0; i < n; ++i) {
		totalDeg += polys[i]->lt;
	}
	if (prod->alloc <= totalDeg + 1) {
		resizePolynomial_DUZP(prod, totalDeg + 1);
	}

	if (idx < 2) {
		if (idx == 0) {
			multiplyPolynomialsPreAlloc_DUZP(polys[1], polys[2], &prod);
		} else if (idx == 1) {
			multiplyPolynomialsPreAlloc_DUZP(polys[0], polys[2], &prod);
		}
		for (int i = 3; i < n; ++i) {
			multiplyPolynomials_DUZP_inp(polys[i], &prod);
		}
	} else {
		multiplyPolynomialsPreAlloc_DUZP(polys[0], polys[1], &prod);
		for (int i = 2; i < n; ++i) {
			if ( i == idx) {
				continue;
			}
			multiplyPolynomials_DUZP_inp(polys[i], &prod);
		}
	}
}

DUZP_t* multiplyManyPolynomials_DUZP(DUZP_t const*const* polys, int n) {
	if (n == 0 || polys == NULL) {
		return NULL;
	}
	if (n < 2) {
		return deepCopyPolynomial_DUZP(polys[0]);
	}

//1: Naive iteartive way
	polysize_t totalDeg = 0;
	for (int i = 0; i < n; ++i) {
		totalDeg += polys[i]->lt;
	}

	DUZP_t* prod = makePolynomial_DUZP(totalDeg + 1);
	multiplyPolynomialsPreAlloc_DUZP(polys[0], polys[1], &prod);
	for (int i = 2; i < n; ++i) {
		multiplyPolynomials_DUZP_inp(polys[i], &prod);
	}
//2: Do a tree;
	// //do a serial reduce for elements of index larger than a power of 2
	// int max = 2;
	// while ( 2*max <= n) { max <<= 1; }
	// DUZP_t* prod, *tmpProd;
	// if (max+1 < n) {
	// 	//at least one mult to do.
	// 	// fprintf(stderr, "doing serial:\n");
	// 	// fprintf(stderr, "\t%d * %d\n", max, max+1);
	// 	prod = multiplyPolynomials_DUZP(polys[max], polys[max+1]);
	// 	for (int i = max+2; i < n; ++i) {
	// 		// fprintf(stderr, "\t%d * %d\n", i-1, i);
	// 		tmpProd = multiplyPolynomials_DUZP(prod, polys[i]);
	// 		freePolynomial_DUZP(prod);
	// 		prod = tmpProd;
	// 	}
	// } else if (max < n) {
	// 	//n = max + 1;
	// 	prod = polisEqualys[max];
	// }

	// // fprintf(stderr, "doing first level:\n");
	// int halfmax = max/2;
	// DUZP_t** tmpPolys = (DUZP_t**) malloc(sizeof(DUZP_t*)*halfmax);
	// for (int i = 0; i < halfmax; ++i) {
	// 	// fprintf(stderr, "\t%d * %d\n", 2*i, 2*i + 1);
	// 	tmpPolys[i] = multiplyPolynomials_DUZP(polys[2*i], polys[2*i + 1]);
	// }

	// //add in the previous serial reduciton of last few elements.
	// if (max < n) {
	// 	// fprintf(stderr, "adding in serial prod to index:%d\n", halfmax-1);
	// 	tmpProd = multiplyPolynomials_DUZP(tmpPolys[halfmax-1], prod);
	// 	freePolynomial_DUZP(tmpPolys[halfmax-1]);
	// 	tmpPolys[halfmax-1] = tmpProd;
	// 	if (max+1 < n) {
	// 		freePolynomial_DUZP(prod);
	// 	}
	// }

	// int step = 1; //actually step is 2 w.r.t original list

	// // fprintf(stderr, "\ndoing generic steps:\n");
	// DUZP_t* tmp;
	// while (step < halfmax) {
	// 	for (int j = 0; j < halfmax; j += (2*step)) {
	// 		// fprintf(stderr, "\t%d * %d\n", j, j+step);
	// 		tmp = multiplyPolynomials_DUZP(tmpPolys[j], tmpPolys[j+step]);
	// 		freePolynomial_DUZP(tmpPolys[j]);
	// 		tmpPolys[j] = tmp;
	// 	}
	// 	step <<= 1;
	// 	// fprintf(stderr, "\n");
	// }
	// for (int i = 1; i < halfmax; ++i) {
	// 	freePolynomial_DUZP(tmpPolys[i]);
	// }

	// prod = tmpPolys[0];

	// free(tmpPolys);

	return prod;
}

void multiplyManyPolynomialsPreAlloc_DUZP(DUZP_t const*const* polys, int n, DUZP_t** prod) {
	if (n == 0 || polys == NULL || prod == NULL) {
		return;
	}


	DUZP_t* c = *prod;

	if (n < 2) {
		const DUZP_t* a = polys[0];
		polysize_t lt = a->lt;
		if (c == NULL) {
			c = makePolynomial_DUZP(a->lt + 1);
			*prod = c;
		} else if (c->alloc < lt + 1) {
			resizePolynomial_DUZP(c, lt + 1);
		}

		polysize_t lt_c = c->lt;
		if (lt_c <= lt) {
			for (polysize_t i = 0; i <= lt_c; ++i) {
				mpz_set(c->coefs[i], a->coefs[i]);
			}
			for (polysize_t i = lt_c + 1; i <= lt; ++i) {
				mpz_init(c->coefs[i]);
				mpz_set(c->coefs[i], a->coefs[i]);
			}
		} else {
			for (polysize_t i = 0; i <= lt; ++i) {
				mpz_set(c->coefs[i], a->coefs[i]);
			}
			for (polysize_t i = lt + 1; i <= lt_c; ++i) {
				mpz_clear(c->coefs[i]);
			}
		}
		return;
	}

	polysize_t maxDeg = 0;
	for (int i = 0; i < n; ++i) {
		maxDeg += polys[i]->lt;
	}

	if (c == NULL) {
		c = makePolynomial_DUZP(maxDeg + 1);
		*prod = c;
	} else if (c->alloc < maxDeg + 1) {
		resizePolynomial_DUZP(c, maxDeg + 1);
	}

	register polysize_t lt_c = c->lt;
	for (polysize_t i = lt_c + 1; i <= maxDeg; ++i) {
		mpz_init(c->coefs[i]);
	}

	lt_c = polys[0]->lt;
	for (polysize_t i = 0; i <= lt_c; ++i) {
		mpz_set(c->coefs[i], polys[0]->coefs[i]);
	}

    register polysize_t i, j;
    register polysize_t jmin = 0;
    register polysize_t jmax = 0;
    register polysize_t lt_a, lt_prod = lt_c;
    register int k = 1;
    const DUZP_t* a;

	//c = c * polys[k];
    while (lt_c == 0 && k < n) {
		a = polys[k];
		lt_a = a->lt;

		for (j = 0; j <= lt_a; ++j) {
			mpz_mul(c->coefs[j], c->coefs[0], a->coefs[j]);
		}

		lt_c = lt_a;
		++k;
    }

	for (; k < n; ++k) {
		//c = c * polys[k];
		a = polys[k];
		lt_a = a->lt;
		if (lt_a == 0) {
			for (j = 0; j <= lt_c; ++j) {
				mpz_mul(c->coefs[j], c->coefs[j], a->coefs[0]);
			}
		} else {

		    lt_prod += lt_a;
			for (i = lt_prod; i >= 0; --i) {
				jmin = MIN_spX(i, lt_c);
				jmax = MAX_spX(i - lt_a, 0);
				mpz_mul(c->coefs[i], c->coefs[jmin], a->coefs[i-jmin]);
				for (j = jmin-1; j >= jmax; --j) {
				    // c[i] += a[j]*b[i-j];
				    mpz_addmul(c->coefs[i], c->coefs[j], a->coefs[i-j]);
				}
		    }
		    lt_c = lt_prod;

		}

	}

	c->lt = lt_c;
}

void multiplyPolynomialsPreAlloc_DUZP(const DUZP_t* a, const DUZP_t* b, DUZP_t** c) {
	if (c == NULL) {
		return;
	}

    // poly * 0 = 0
    if (a == NULL || b == NULL) {
    	DUZP_t* cc = *c;
    	if (!isZero_DUZP(cc)) {
    		mpz_set_ui(cc->coefs[0], 0ul);
    		for (polysize_t i = 1; i <= cc->lt; ++i) {
    			mpz_clear(cc->coefs[i]);
    		}
    		cc->lt = 0;
    	}
    	// freePolynomial_DUZP(*c);
    	// *c = NULL;
    	return;
    }

    if (a->lt == 0) {
    	// return multiplyByInteger_DUZP(b, a->coefs[0]);
    	DUZP_t* cc = *c;
    	polysize_t bSize = b->lt + 1;
    	if (cc == NULL) {
    		cc = makePolynomial_DUZP(bSize);
    		*c = cc;
    	}
    	if (cc->alloc < bSize) {
			resizePolynomial_DUZP(cc, bSize);
    	}
    	polysize_t cSize = cc->lt + 1;
    	cc->lt = b->lt;

		mpz_t* ccoefs = cc->coefs;
		mpz_t* bcoefs = b->coefs;

		if (cSize < bSize) {
			for (polysize_t i = 0; i < cSize; ++i) {
				mpz_mul(ccoefs[i], bcoefs[i], a->coefs[0]);
			}
			for (polysize_t i = cSize; i < bSize; ++i) {
				mpz_init(ccoefs[i]);
				mpz_mul(ccoefs[i], bcoefs[i], a->coefs[0]);
			}
		} else {
			for (polysize_t i = 0; i < bSize; ++i) {
				mpz_mul(ccoefs[i], bcoefs[i], a->coefs[0]);
			}
			for (polysize_t i = bSize; i < cSize; ++i) {
				mpz_clear(ccoefs[i]);
			}
		}
		return;
    }

    if (b->lt == 0) {
    	// return multiplyByInteger_DUZP(a, b->coefs[0]);
    	DUZP_t* cc = *c;
    	polysize_t aSize = a->lt + 1;
    	if (cc == NULL) {
    		cc = makePolynomial_DUZP(aSize);
    		*c = cc;
    	}
    	polysize_t cSize = cc->lt + 1;
    	cc->lt = a->lt;
    	if (cc->alloc < aSize) {
			resizePolynomial_DUZP(cc, aSize);
    	}
		mpz_t* ccoefs = cc->coefs;
		mpz_t* acoefs = a->coefs;
		if (cSize < aSize) {
			for (polysize_t i = 0; i < cSize; ++i) {
				mpz_mul(ccoefs[i], acoefs[i], b->coefs[0]);
			}
			for (polysize_t i = cSize; i < aSize; ++i) {
				mpz_init(ccoefs[i]);
				mpz_mul(ccoefs[i], acoefs[i], b->coefs[0]);
			}
		} else {
			for (polysize_t i = 0; i < aSize; ++i) {
				mpz_mul(ccoefs[i], acoefs[i], b->coefs[0]);
			}
			for (polysize_t i = aSize; i < cSize; ++i) {
				mpz_clear(ccoefs[i]);
			}
		}
		return;
    }

	DUZP_t* cc = *c;

    polysize_t lt_a = a->lt;
    polysize_t lt_b = b->lt;
    polysize_t lt_c = lt_a + lt_b;

    if (cc == NULL) {
    	// fprintf(stderr, "making poly\n" );
    	cc = makePolynomial_DUZP(lt_c+1);
    	*c = cc;
    }

    if (cc->alloc < lt_c + 1) {
    	// fprintf(stderr, "resize poly\n" );
		resizePolynomial_DUZP(cc, lt_c+1);
    }

    polysize_t lt_c_old = cc->lt;

    mpz_t* a_coefs = a->coefs;
    mpz_t* b_coefs = b->coefs;
    mpz_t* c_coefs = cc->coefs;

    cc->lt = lt_c;

    register polysize_t i, j;
    register polysize_t jmin = 0;
    register polysize_t jmax = 0;


    if (lt_c_old > lt_c) {
    	for (i = 0; i <= lt_c; ++i) {
	    	// fprintf(stderr, "pralloc %p set2 0 %ld\n", cc, i);
    		mpz_set_ui(c_coefs[i], 0ul);
			jmin = MIN_spX(i, lt_a);
			jmax = MAX_spX(i - lt_b, 0);
			for (j = jmax; j <= jmin; j++) {
			    // c[i] += a[j]*b[i-j];
			    mpz_addmul(c_coefs[i], a_coefs[j], b_coefs[i-j]);
			}
	    }
    	for (int i = lt_c + 1; i < lt_c_old; ++i) {
	    	// fprintf(stderr, "pralloc %p clear %ld\n", cc, i);
			mpz_clear(c_coefs[i]);
    	}
    } else {
	    for (i = 0; i <= lt_c_old; ++i) {
	    	// fprintf(stderr, "pralloc %p set 0 %ld\n", cc, i);
    		mpz_set_ui(c_coefs[i], 0ul);
			jmin = MIN_spX(i, lt_a);
			jmax = MAX_spX(i - lt_b, 0);
			for (j = jmax; j <= jmin; j++) {
			    // c[i] += a[j]*b[i-j];
			    mpz_addmul(c_coefs[i], a_coefs[j], b_coefs[i-j]);
			}
	    }
	    for (i = lt_c_old + 1; i <= lt_c; ++i) {
	    	// fprintf(stderr, "pralloc %p init %ld\n", cc, i);
			mpz_init_set_ui(c_coefs[i], 0ul);
			jmin = MIN_spX(i, lt_a);
			jmax = MAX_spX(i - lt_b, 0);
			for (j = jmax; j <= jmin; j++) {
			    // c[i] += a[j]*b[i-j];
			    mpz_addmul(c_coefs[i], a_coefs[j], b_coefs[i-j]);
			}
	    }
    }
}

void multiplyPolynomials_DUZP_inp(const DUZP_t* a, DUZP_t** bb) {
	if (bb == NULL) {
		return;
	}
	if (isZero_DUZP(a) || isZero_DUZP(*bb)) {
		freePolynomial_DUZP(*bb);
		*bb = NULL;
		return;
	}

	if (a->lt == 0) {
		multiplyByInteger_DUZP_inp(*bb, a->coefs[0]);
		return;
	} else if ((*bb)->lt == 0) {
		DUZP_t* ret = multiplyByInteger_DUZP(a, (*bb)->coefs[0]);
		freePolynomial_DUZP(*bb);
		*bb = ret;
		return;
	}


	DUZP_t* b = *bb;
	polysize_t lt_a = a->lt;
    polysize_t lt_b = b->lt;
    polysize_t lt_c = lt_a + lt_b;

    if (b->alloc < lt_c + 1) {
    	resizePolynomial_DUZP(b, lt_c+1);
    }

    register polysize_t i, j;
    register polysize_t jmin = 0;
    register polysize_t jmax = 0;

    mpz_t* __restrict__ a_coefs = a->coefs;
    mpz_t* __restrict__ b_coefs = b->coefs;

    for (i = lt_b + 1; i <= lt_c; ++i) {
    	mpz_init(b_coefs[i]);
    }

    for (i = lt_c; i >= 0; --i) {
		jmin = MIN_spX(i, lt_b);
		jmax = MAX_spX(i - lt_a, 0);
		mpz_mul(b_coefs[i], b_coefs[jmin], a_coefs[i-jmin]);
		for (j = jmin-1; j >= jmax; --j) {
		    // c[i] += a[j]*b[i-j];
		    mpz_addmul(b_coefs[i], b_coefs[j], a_coefs[i-j]);
		}
    }

    b->lt = lt_c;
}


void multiplyByBinomial_DUZP_inp(DUZP_t** a_ptr, const mpz_t b) {
	if (a_ptr == NULL) {
		return;
	}

	DUZP_t* a = *a_ptr;
	if (a == NULL) {
		a = makePolynomial_DUZP(2);

		mpz_init_set_ui(a->coefs[0], 1ul);
		mpz_init_set(a->coefs[1], b);
		a->lt = 1;
		*a_ptr = a;
		return;
	}

	if (a->alloc < a->lt + 1) {
		resizePolynomial_DUZP(a, a->lt*2); //ammortize resize cost hopefully;
	}


	mpz_t tmp;
	mpz_init(tmp);
	mpz_mul(tmp, a->coefs[0], b);
	mpz_neg(tmp, tmp);

	for (int i = 0; i < a->lt; ++i) {
		mpz_submul(a->coefs[i], a->coefs[i+1], b);
	}

	//shift right
	memmove(a->coefs+1, a->coefs, sizeof(mpz_t)*(a->lt+1));
	//strict copy of struct contents; don't clear tmp.
	a->coefs[0][0] = tmp[0];

	++(a->lt);
}


/**
 * Does b divide a ?
 */
int divideTest_DUZP(const DUZP_t* a, const DUZP_t* b, DUZP_t** q) {
	if (a == NULL) {
		return 1;
	}
	if (b == NULL) {
		fprintf(stderr, "Tried to divide by 0 in divideTest_DUZP!\n");
		return 0;
	}

	polysize_t lt_a = a->lt;
	polysize_t lt_b = b->lt;

	if (lt_b > lt_a) {
		return 0;
	}

	polysize_t diff_deg = lt_a - lt_b;

	mpz_t* lcb = &(b->coefs[lt_b]);

	DUZP_t* rr = deepCopyPolynomial_DUZP(a);
	mpz_t* rr_coefs = rr->coefs;

	DUZP_t* quo = makePolynomial_DUZP(diff_deg + 1);
	mpz_clear(quo->coefs[0]);

	polysize_t i = 0, j = 0;
    for (i = diff_deg; i >= 0; --i) {
    	if (mpz_sgn(rr_coefs[lt_b + i]) != 0) {
    		if (!mpz_divisible_p(rr_coefs[lt_b + i], *lcb)) {
			    rr->lt = lt_a;
			    freePolynomial_DUZP(rr);
			    for (j = i+1; j <= diff_deg; ++j) {
			    	mpz_clear(quo->coefs[j]);
			    }
			    quo->lt = -1;
			    freePolynomial_DUZP(quo);
    			return 0;
    		}

    		mpz_init(quo->coefs[i]);
    		mpz_divexact(quo->coefs[i], rr_coefs[lt_b + i], *lcb);

  			mpz_set_si(rr_coefs[lt_b + i], 0l);
    	    --(rr->lt);
    	    for (j = 0; j < lt_b; ++j) {
    	    	mpz_submul(rr_coefs[i+j], b->coefs[j], quo->coefs[i]);
    	    }
    	} else {
    		mpz_init(quo->coefs[i]); //quo[i] = 0;
    	}
    }

    for (i = rr->lt; i >= 0; --i) {
    	if (mpz_sgn(rr_coefs[i]) != 0) {
		    rr->lt = lt_a;
		    freePolynomial_DUZP(rr);
		    quo->lt = diff_deg;
		   	freePolynomial_DUZP(quo);
    		return 0;
    	}
    }

    rr->lt = lt_a;
    freePolynomial_DUZP(rr);
	quo->lt = diff_deg;

	if (q != NULL) {
		*q = quo;
	} else {
		freePolynomial_DUZP(quo);
	}
    return 1;
}

int divideByMonicLinear_DUZP(const DUZP_t* a, const mpz_t b, DUZP_t** q, mpz_t* rem) {
	if (q == NULL && rem == NULL) {
		return -1;
	}

	if (a == NULL) {
		if (q != NULL) {
			freePolynomial_DUZP(*q);
			*q = NULL;
		}
		if (rem != NULL) {
			mpz_set_ui(*rem, 0ul);
		}
		return 1;
	}

	DUZP_t* quo = makePolynomial_DUZP(a->lt);
	quo->lt = a->lt - 1;
	mpz_init_set(quo->coefs[quo->lt], a->coefs[a->lt]);
	for (int i = quo->lt-1; i >= 0; --i) {
		mpz_init_set(quo->coefs[i], a->coefs[i+1]);
		mpz_addmul(quo->coefs[i], quo->coefs[i+1], b);
	}

	int ret;
	if (rem != NULL) {
		mpz_set(*rem, a->coefs[0]);
		mpz_addmul(*rem, quo->coefs[0], b);
		ret = (mpz_sgn(*rem) == 0);
	} else {
		mpz_t temp;
		mpz_init_set(temp, a->coefs[0]);
		mpz_addmul(temp, quo->coefs[0], b);
		ret = (mpz_sgn(temp) == 0);
		mpz_clear(temp);
	}

	if (q != NULL) {
		*q = quo;
	} else {
		freePolynomial_DUZP(quo);
	}

	return ret;

}

int dividePolynomials_DUZP(const DUZP_t* a, const DUZP_t* b, DUZP_t** q, DUZP_t** r) {

    // a / 0 = Error
    if (isZero_DUZP (b)) {
    	fprintf (stderr, "DUZP Error: Division by zero!\n");
    	exit (1);
    }

    // 0 / b = 0
    if (isZero_DUZP (a)) {
    	if (r != NULL) {
    		*r = NULL;
    	}
    	if (q != NULL) {
    		*q = NULL;
    	}
    	return 1;
    }

    // a/a = 1
    if (isEqual_DUZP (a, b)) {
		if (r != NULL) {
			*r = NULL;
		}
		if (q != NULL) {
			*q = makeConstPolynomial_DUZP(1, 1l);
		}
		return 1;
    }

    polysize_t deg_a = a->lt;
    polysize_t deg_b = b->lt;

    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
    	if (r != NULL) {
    		*r = deepCopyPolynomial_DUZP (a);
    	}
    	if (q != NULL) {
			*q = NULL;
    	}
    	return 1;
    }

    polysize_t diff_deg = deg_a - deg_b;

    DUZP_t* rr = deepCopyPolynomial_DUZP(a); // rr = a;
    DUZP_t* qq = makePolynomial_DUZP(diff_deg+1);

    for (int i = 1; i <= diff_deg; ++i) {
    	mpz_init(qq->coefs[i]);
    }
    qq->lt = diff_deg;

    register polysize_t i, j;
    mpz_t* lc_b = b->coefs + b->lt;
    mpz_t* rCoefs = rr->coefs;
    mpz_t* qCoefs = qq->coefs;

    // a = b*qq + rr
    for (i = diff_deg; i >= 0; --i) {
    	if (mpz_sgn(rCoefs[deg_b + i])) {
    		if (mpz_divisible_p(rCoefs[deg_b + i], *lc_b)) {
    			mpz_divexact(qCoefs[i], rCoefs[deg_b + i], *lc_b);

    			// mpz_set_ui(rCoefs[deg_b + i], 0ul);

    			mpz_clear(rCoefs[deg_b + i]);
				--(rr->lt);

    			for (j = deg_b-1; j >= 0; --j) {
    				mpz_submul(rCoefs[i+j], b->coefs[j], qCoefs[i]);
    			}
	    	} else {
	    		//bail out?
	    		freePolynomial_DUZP(rr);
	    		freePolynomial_DUZP(qq);
	    		if (r != NULL) {
	    			*r = NULL;
	    		}
	    		if (q != NULL) {
	    			*q = NULL;
	    		}
	    		return 0;
	    	}
		} else {
			mpz_clear(rCoefs[deg_b + i]);
			--(rr->lt);
		}
    }


  //   // a = b*qq + rr
  //   for (i = diff_deg; i >= 0; i--) {
		// lc_rr_in = idxCoeffInForm_spX (rr, deg_b + i);
		// /* lc_rr_in = smallprimefield_convert_in (lc_rr_in, Pptr); */

		// if (lc_rr_in) {
		// 	lc_rr_in = smallprimefield_mul (lc_rr_in, lc_b_in, Pptr);

		// 	qq->elems[i] =  lc_rr_in; // smallprimefield_convert_out (lc_rr_in, Pptr);

		// 	rr->elems[deg_b+i] = 0;
		// 	rr->lt -= 1;
		// 	rr->alloc -= 1;
		// 	for (j = deg_b-1;j >= 0; j--) {
		// 		/* b_in = smallprimefield_convert_in (b->elems[j], Pptr); */
		// 		b_in = smallprimefield_mul (b->elems[j], lc_rr_in, Pptr);
		// 		/* b_in = smallprimefield_convert_out (b_in, Pptr); */

		// 		rr->elems[i+j] = smallprimefield_sub (rr->elems[i+j], b_in, Pptr);
		// 	}
		// }
  //   }

    for (i = rr->lt; i >= 1; --i) {
    	if (mpz_sgn(rCoefs[i]) == 0) {
    		mpz_clear(rCoefs[i]);
    		--(rr->lt);
    	} else {
    		break;
    	}
    }

    if (q != NULL) {
    	*q = qq;
    }
	if (r != NULL) {
    	*r = rr;
	}

	return 1;
}

//get the remainder of an exact division
//return 1 iff the division was exact
//If it returns 0 then a is left in a modified state.
int remainder_DUZP_inp(DUZP_t** aa, DUZP_t* b) {
	if (aa == NULL) {
		return 0;
	}

    // a / 0 =
    if (isZero_DUZP (b)) {
    	fprintf (stderr, "DUZP Error: Division by zero!\n");
    	exit (1);
    }

    // 0 / b = 0
    if (isZero_DUZP (*aa)) {
    	return 1;
    }

    // a / const =
    if (b->lt == 0) {
    	mpz_t cont;
    	mpz_init(cont);
    	content_DUZP(*aa, cont);
    	if (mpz_divisible_p(cont, b->coefs[0])) {
			freePolynomial_DUZP(*aa);
			*aa = NULL;
			mpz_clear(cont);
    		return 1;
    	}
		mpz_clear(cont);
		return 0;
    }

    // a/a = 1
    if (isEqual_DUZP (*aa, b)) {
		freePolynomial_DUZP(*aa);
		*aa = NULL;
		return 1;
    }


    DUZP_t* r = *aa;
    polysize_t deg_a = r->lt;
    polysize_t deg_b = b->lt;

    if (deg_a < deg_b) {
    	//rem is itself
    	return 1;
    }

    polysize_t diff_deg = deg_a - deg_b;
    polysize_t i, j;
    mpz_t* lc_b = b->coefs + b->lt;
    mpz_t* rCoefs = r->coefs;

    // a = b*qq + rr
    for (i = diff_deg; i >= 0; --i) {
    	if (mpz_sgn(rCoefs[deg_b + i])) {
    		if (mpz_divisible_p(rCoefs[deg_b + i], *lc_b)) {
    			mpz_divexact(rCoefs[deg_b + i], rCoefs[deg_b + i], *lc_b);

    			for (j = deg_b-1; j >= 0; --j) {
    				mpz_submul(rCoefs[i+j], b->coefs[j], rCoefs[deg_b + i]);
    			}

    			mpz_clear(rCoefs[deg_b + i]);
				--(r->lt);

	    	} else {
	    		//	?
	    		fprintf(stderr, "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@\nNahhhh in divide\n\n" );
	    		return 0;
	    	}
		} else {
			mpz_clear(rCoefs[deg_b + i]);
			--(r->lt);
		}
    }

    for (i = r->lt; i >= 1; --i) {
    	if (mpz_sgn(rCoefs[i]) == 0) {
    		mpz_clear(rCoefs[i]);
    		--(r->lt);
    	} else {
    		break;
    	}
    }
    return 1;
}


//get the Rem(aa,b) modulo mod;
//return 1 iff the division was exact
//If it returns 0 then a is left in a modified state.
int remainderMod_DUZP_inp(DUZP_t** aa, const DUZP_t* b, const mpz_t mod) {
	if (aa == NULL) {
		return 0;
	}

    // a / 0 =
    if (isZero_DUZP (b)) {
    	fprintf (stderr, "DUZP Error: Division by zero!\n");
    	exit (1);
    }

    // 0 / b = 0
    if (isZero_DUZP (*aa)) {
    	return 1;
    }

    // a / const =
    if (b->lt == 0) {
    	mpz_t lcinv;
    	mpz_init(lcinv);
    	if(!mpz_invert(lcinv, b->coefs[0], mod)) {
    		//lc(b) is a zero divisor
    		mpz_clear(lcinv);
    		return 0;
    	}

    	multiplyByInteger_DUZP_inp(*aa, lcinv);
    	mpz_clear(lcinv);
    	applyModulo_DUZP_inp(*aa, mod);
		return 1;
    }

    // a/a = 1
    if (isEqual_DUZP (*aa, b)) {
		freePolynomial_DUZP(*aa);
		*aa = NULL;
		return 1;
    }


    DUZP_t* r = *aa;
    polysize_t deg_a = r->lt;
    polysize_t deg_b = b->lt;

    if (deg_a < deg_b) {
    	//rem is itself
    	applyModulo_DUZP_inp(*aa, mod);
    	return 1;
    }

	mpz_t lcinv;
	mpz_init(lcinv);
	if(!mpz_invert(lcinv, b->coefs[b->lt], mod)) {
		//lc(b) is a zero divisor
		mpz_clear(lcinv);
		return 0;
	}

    polysize_t diff_deg = deg_a - deg_b;
    polysize_t i, j;
    // mpz_t* lc_b = b->coefs + b->lt;
    mpz_t* rCoefs = r->coefs;

    // a = b*qq + rr
    for (i = diff_deg; i >= 0; --i) {
    	if (mpz_sgn(rCoefs[deg_b + i])) {
    		mpz_mul(rCoefs[deg_b + i], rCoefs[deg_b + i], lcinv);

			for (j = deg_b-1; j >= 0; --j) {
				mpz_submul(rCoefs[i+j], b->coefs[j], rCoefs[deg_b + i]);
			}

			mpz_clear(rCoefs[deg_b + i]);
			--(r->lt);
		} else {
			mpz_clear(rCoefs[deg_b + i]);
			--(r->lt);
		}
    }
	mpz_clear(lcinv);

    for (i = r->lt; i >= 1; --i) {
    	if (mpz_sgn(rCoefs[i]) == 0) {
    		mpz_clear(rCoefs[i]);
    		--(r->lt);
    	} else {
    		break;
    	}
    }

    applyModulo_DUZP_inp(r, mod);
    return 1;
}


/**
 * Get the GCD of two DUZP primitive polynomials.
 */
DUZP_t* primitiveGCD_DUZP(const DUZP_t* a, const DUZP_t* b) {
	if (a == NULL || b == NULL) {
		return NULL;
	}

	mpz_t lcg;
	mpz_init(lcg);
	mpz_gcd(lcg, a->coefs[a->lt], b->coefs[b->lt]);

	polysize_t d = MIN_spX(a->lt, b->lt);

	DUZP_t* g_work = makePolynomial_DUZP(d+1); //gcd currrently being built up by CRT
	mpz_t* g_work_coefs = g_work->coefs;
	for (polysize_t i = 1; i <= d; ++i) {
		mpz_init(g_work_coefs[i]);
	}
	g_work->lt = d;

	mpz_t m; //product of primes;
	mpz_init_set_ui(m, 1l);
	mpz_t halfm;
	mpz_init(halfm);
	mpz_t newm;
	mpz_init(newm);
	mpz_t lcgm; //product of lcg and m to test for bad primes
	mpz_init_set(lcgm, lcg);

	mpz_t s, t, mpz_prime; //for CRT
	mpz_inits(s, t, mpz_prime, NULL);

	int primeIdx = 0;
	Prime_ptr Pptr[1];

	duspoly_t* moda = makePolynomial_spX(a->lt + 1);
	duspoly_t* modb = makePolynomial_spX(b->lt + 1);
	duspoly_t* modg = NULL;

	//TODO instead of fibs maybe just store two g_work_coefs and
	//do divide test when it stops changing.
	int fibs[] = {0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584};
	int fibIdx = 0, nFibs = 17;

	for(; primeIdx < n_prime64_ptr; ++primeIdx) {
		*Pptr = prime64_ptr[primeIdx];
		if (mpz_divisible_ui_p(lcgm, (unsigned long long) Pptr->prime)) {
			//bad prime
			continue;
		}

		//obtain modular images for this prime
		convertToPrimeField_DUZP_inp(a, Pptr, &moda);
		convertToPrimeField_DUZP_inp(b, Pptr, &modb);

//1
		// plainGCDInForm_spX(moda, modb, &modg, Pptr);
//2
		GCDInForm_spX (moda, modb, &modg, Pptr);

		if (modg->lt == 0) {
			//then gcd is 1
			DUZP_t* ret = makePolynomial_DUZP(1);
			mpz_init(ret->coefs[0]);
			mpz_set_si(ret->coefs[0], 1l);
			ret->lt = 0;
			return ret;
		} else if (modg->lt < g_work->lt) {
			//then previous images were bad, start over
			polysize_t lt = modg->lt;
			polysize_t i = 0;
			elem_t* primeElems = modg->elems;
			long long int halfP = (Pptr->prime - 1) >> 1;
			long long int tmp;
			long long int lcg_pf = smallprimefield_convert_in(mpz_fdiv_ui(lcg, Pptr->prime), Pptr);
			for (; i <= lt; ++i) {
				tmp = smallprimefield_mul(primeElems[i], lcg_pf, Pptr);
				tmp = smallprimefield_convert_out(tmp, Pptr);
				if (tmp > halfP) {
					tmp -= Pptr->prime;
				}
				mpz_set_si(g_work_coefs[i], tmp);
			}
			lt = g_work->lt;
			for (; i <= lt; ++i) {
				mpz_clear(g_work_coefs[i]);
			}
			g_work->lt = modg->lt;

			mpz_set_si(m, Pptr->prime);
			mpz_mul(lcgm, m, lcg);

			//TODO fix this so that GCD puts result into already allocated g
			freePolynomial_spX (&modg);
			modg = NULL;
		} else if (modg->lt > g_work->lt) {
			//current image is bad, ignore it and move on
			//TODO fix this so that GCD puts result into already allocated g
			freePolynomial_spX (&modg);
			modg = NULL;
			continue;
		} else {

			long long int lcg_pf = smallprimefield_convert_in(mpz_fdiv_ui(lcg, Pptr->prime), Pptr);
			scalarMulPolynomialInForm_spX_inp (&modg, lcg_pf, Pptr);

			mpz_set_si(mpz_prime, Pptr->prime);
			mpz_gcdext(lcgm, s, t, mpz_prime, m);

			mpz_mul(newm, m, mpz_prime);
			// mpz_sub_ui(halfm, mpz_prime, 1ul);
			mpz_sub_ui(halfm, newm, 1ul);
			mpz_fdiv_q_2exp(halfm, halfm, 1ul);

			//CRT
			polysize_t lt = modg->lt;
			unsigned long long modg_out;
			for (int i = 0; i <= lt; ++i) {
				//only use lcgm as a temporary holder since it gets overwritten at end of loop
				//s*prime + t*m = 1; t is s for current prime

				//naive application of bezout coefficients
				// mpz_mul(g_work_coefs[i], g_work_coefs[i], s);
				// mpz_mod(g_work_coefs[i], g_work_coefs[i], m);
				// mpz_mul(g_work_coefs[i], g_work_coefs[i], mpz_prime);
				// mpz_mul_si(lcgm, t, smallprimefield_convert_out(modg->elems[i], Pptr));
				// mpz_mod(lcgm, lcgm, mpz_prime);
				// mpz_mul(lcgm, lcgm, m);
				// mpz_add(g_work_coefs[i], g_work_coefs[i], lcgm);

				//optimized application of bezout coefs by rearranging
				// s*prime + t*m = 1 after multiplying through by modg[i].
				// we have the two extended gcd equations: s*prime + t*m = 1, s'*m + t'*prime = 1;
				// therefore s = t', t = s'. Using this fact to rearrange as follows:
				// 1) modg[i]*s'*m + modg[i]*t'*prime = modg[i]
				// 2) modg[i] - modg[i]*s*prime = modg[i]*s'*m.
				// Notice RHS of eq. is a term of CRT summation for modg[i] (except for mod prime).
				// 3) we want (g[i]*s mod m)*prime + (modg[i]*s' mod prime)*m. Replace mods with u*m and w*prime, respectively.
				// 4) (g[i]*s + u*m)*prime + (modg[i]*s' + w*prime)*m
				// 5) g[i]*s*prime + modg[i]*s'*m - (u+w)*(prime*m)
				// 6) (5) & (2) -> g[i]*s*prime + modg[i] - modg[i]*s*prime + (u*w)*(prime*m);
				// 7) modg[i] + s*prime*(g[i] - modg[i]) + (u*w)*(prime*m);
				// Finally, we implement (7) mod (prime*m)
				// This reduces cost from 4 mults, 2 mods, 1 add to 2 mults, 1 mod, 2 adds.

				modg_out = (unsigned long long) smallprimefield_convert_out(modg->elems[i], Pptr);
				mpz_sub_ui(g_work_coefs[i], g_work_coefs[i], modg_out);
				mpz_mul(g_work_coefs[i], g_work_coefs[i], s);
				mpz_mul(g_work_coefs[i], g_work_coefs[i], mpz_prime);
				mpz_add_ui(g_work_coefs[i], g_work_coefs[i], modg_out);
				mpz_mod(g_work_coefs[i], g_work_coefs[i], newm);
				//normalize in symmetric range
				if(mpz_cmp(g_work_coefs[i], halfm) > 0) {
					mpz_sub(g_work_coefs[i], g_work_coefs[i], newm);
				}
			}

			mpz_swap(m, newm);
			mpz_mul(lcgm, m, lcg);

			freePolynomial_spX (&modg);
		}

		int doDivideTest = fibIdx < nFibs ? (primeIdx == fibs[fibIdx++]) : 1;
		if (doDivideTest) {
			DUZP_t* primG = primitivePart_DUZP(g_work);
			if (divideTest_DUZP(a, primG, NULL) && divideTest_DUZP(b, primG, NULL)) {
				mpz_clears(m, halfm, newm, lcgm, s, t, mpz_prime, NULL);
				freePolynomial_spX (&moda);
				freePolynomial_spX (&modb);
				freePolynomial_DUZP(g_work);

				return primG;
			}
		}

	}

	mpz_clears(lcg, lcgm, s, t, mpz_prime, newm, halfm, m, NULL);


	//all primes failed....
	DUZP_t* ret = makePolynomial_DUZP(1);
	mpz_set_si(ret->coefs[0], 1l);
	ret->lt = 0;

	mpz_clears(m, halfm, newm, lcgm, s, t, mpz_prime, NULL);
	freePolynomial_spX (&moda);
	freePolynomial_spX (&modb);
	freePolynomial_DUZP(g_work);

	return ret;
}

void differentiate_DUZP_inp(DUZP_t** a) {

	if (isZero_DUZP(*a)) {
		return;
	}

	mpz_t t1;
	mpz_init(t1);

	for (int i=0; i<(*a)->lt; ++i) {
		mpz_set_ui(t1,i+1);
		mpz_mul(t1,t1,(*a)->coefs[i+1]);
		mpz_set((*a)->coefs[i],t1);
	}

	if ((*a)->lt != 0) {
		mpz_clear((*a)->coefs[(*a)->lt]);
		(*a)->lt--;
	}
	else
		mpz_set_ui((*a)->coefs[(*a)->lt],0);

	mpz_clear(t1);


}

DUZP_t* GCD_DUZP(DUZP_t* a, DUZP_t* b) {

	mpz_t conta;
	mpz_t contb;
	mpz_init(conta);
	mpz_init(contb);
	DUZP_t* ap = primitivePartAndContent_DUZP(a, conta);
	DUZP_t* bp = primitivePartAndContent_DUZP(b, contb);

	DUZP_t* gp = primitiveGCD_DUZP(ap,bp);

	mpz_gcd(conta, conta, contb);
	multiplyByInteger_DUZP_inp(gp, conta);
	mpz_clear(conta);
	mpz_clear(contb);
	freePolynomial_DUZP(ap);
	freePolynomial_DUZP(bp);
	return gp;
}

DUZP_t* convertFromAltArr_DUZP(const AltArr_t* aa) {
	if (aa == NULL || aa->alloc == 0 || aa->size == 0) {
		return NULL;
	}
	if (aa->nvar != 1) {
		fprintf(stderr, "BPAS ERROR: Trying to convert multivariate AltArr_t to DUZP\n");
		exit(1);
	}

	DUZP_t* retPoly = NULL;
	if (aa->unpacked) {
		degree_t* degs = (degree_t*) aa->elems->degs;
		degree_t curDeg = degs[0];
		retPoly = makePolynomial_DUZP(curDeg + 1);
		retPoly->lt = curDeg;

		for (int i = 1; i <= curDeg; ++i) {
			mpz_init(retPoly->coefs[i]);
		}

		int aaSize = aa->size;
		AAElem_t* aaElems = aa->elems;
		for (int i = 0; i < aaSize; ++i) {
			curDeg = degs[i];
			if (mpz_cmp_si(mpq_denref(aaElems[i].coef), 1l) != 0) {
				fprintf(stderr, "BPAS ERROR: Trying to convert non-integer AltArr_t to DUZP\n");
				exit(1);
			}
			mpz_set(retPoly->coefs[curDeg], mpq_numref(aaElems[i].coef));
		}
	} else {
		AAElem_t* aaElems = aa->elems;
		degrees_t curDeg = aaElems->degs;
		retPoly = makePolynomial_DUZP(curDeg + 1);
		retPoly->lt = curDeg;

		for (int i = 1; i <= curDeg; ++i) {
			mpz_init(retPoly->coefs[i]);
		}

		int aaSize = aa->size;
		for (int i = 0; i < aaSize; ++i) {
			curDeg = aaElems[i].degs;
			if (mpz_cmp_si(mpq_denref(aaElems[i].coef), 1l) != 0) {
				fprintf(stderr, "BPAS ERROR: Trying to convert non-integer AltArr_t to DUZP\n");
				exit(1);
			}

			mpz_set(retPoly->coefs[curDeg], mpq_numref(aaElems[i].coef));
		}
	}

	return retPoly;
}

AltArr_t* convertToAltArr_DUZP(const DUZP_t* poly) {
	if (poly == NULL) {
		return NULL;
	}

	int polySize = 0;
	for (int i = 0; i <= poly->lt; ++i){
		if (mpz_sgn(poly->coefs[i]) != 0) {
			++polySize;
		}
	}

	AltArr_t* retPoly = makePolynomial_AA(polySize, 1);
	polySize = 0;
	for (int i = poly->lt; i >= 0; --i){
		if (mpz_sgn(poly->coefs[i]) != 0) {
			mpq_init(retPoly->elems[polySize].coef);
			mpz_set(mpq_numref(retPoly->elems[polySize].coef), poly->coefs[i]);
			retPoly->elems[polySize].degs = i;
			++polySize;
		}
	}
	retPoly->size = polySize;

	return retPoly;
}

DUZP_t* convertFromDUSP_DUZP(const duspoly_t* a, const Prime_ptr* Pptr) {
	if (isZero_spX(a)) {
		return NULL;
	}

	polysize_t size = a->lt + 1;
	elem_t* pelems = a->elems;

	DUZP_t* ret = makePolynomial_DUZP(size);
	mpz_set_si(ret->coefs[0], smallprimefield_convert_out(pelems[0], Pptr));
	for(polysize_t i = 1; i < size; ++i) {
		mpz_init(ret->coefs[i]);
		mpz_set_si(ret->coefs[i], smallprimefield_convert_out(pelems[i], Pptr));
	}

	ret->lt = a->lt;
	return ret;
}

DUZP_t* convertFromDUSP_DUZP_inRange (const duspoly_t* a, const Prime_ptr* Pptr) {
	if (isZero_spX(a)) {
		return NULL;
	}

	polysize_t size = a->lt + 1;
	elem_t* pelems = a->elems;

	DUZP_t* ret = makePolynomial_DUZP(size);
	mpz_set_si(ret->coefs[0], smallprimefield_convert_out(pelems[0], Pptr));
	if (mpz_cmp_si (ret->coefs[0], Pptr->prime>>1) > 0) {
		mpz_sub_ui (ret->coefs[0], ret->coefs[0], (unsigned long) Pptr->prime);
	}

	for(polysize_t i = 1; i < size; ++i) {
		mpz_init(ret->coefs[i]);
		mpz_set_si(ret->coefs[i], smallprimefield_convert_out(pelems[i], Pptr));
		if (mpz_cmp_si (ret->coefs[i], Pptr->prime>>2) > 0) {
			mpz_sub_ui (ret->coefs[i], ret->coefs[i], (unsigned long) Pptr->prime);
		}
	}

	ret->lt = a->lt;
	return ret;
}

duspoly_t* convertToDUSP_DUZP(const DUZP_t* a, const Prime_ptr* Pptr) {
	if (isZero_DUZP(a)) {
		return NULL;
	}

	polysize_t size = a->lt + 1;
	duspoly_t* ret = makePolynomial_spX(size);

	register unsigned long P = (unsigned long) Pptr->prime;
	register elem_t tmp;

	for (polysize_t i = 0; i < size; ++i) {
		tmp = mpz_fdiv_ui(a->coefs[i], P);
/*		fprintf(stderr, "tmp: %lld\n", tmp);*/
/*		fprintf(stderr, "prime: %lld\n", Pptr->prime);*/
		ret->elems[i] = smallprimefield_convert_in(tmp, Pptr);
	}

	ret->lt = size - 1;
	return ret;
}


DUZP_t* convertFromAltArrZ_DUZP(const AltArrZ_t* aa) {
	if (aa == NULL || aa->alloc == 0 || aa->size == 0) {
		return NULL;
	}
	if (aa->nvar != 1) {
		fprintf(stderr, "BPAS ERROR: Trying to convert multivariate AltArr_t to DUZP\n");
		exit(1);
	}

	DUZP_t* retPoly = NULL;
	if (aa->unpacked) {
		degree_t* degs = (degree_t*) aa->elems->degs;
		degree_t curDeg = degs[0];
		retPoly = makePolynomial_DUZP(curDeg + 1);
		retPoly->lt = curDeg;

		for (int i = 0; i <= curDeg; ++i) {
			mpz_init(retPoly->coefs[i]);
		}

		int aaSize = aa->size;
		AAZElem_t* aaElems = aa->elems;
		for (int i = 0; i < aaSize; ++i) {
			curDeg = degs[i];
			mpz_set(retPoly->coefs[curDeg], aaElems[i].coef);
		}
	} else {
		AAZElem_t* aaElems = aa->elems;
		degrees_t curDeg = aaElems->degs;
		retPoly = makePolynomial_DUZP(curDeg + 1);
		retPoly->lt = curDeg;

		for (int i = 0; i <= curDeg; ++i) {
			mpz_init(retPoly->coefs[i]);
		}

		int aaSize = aa->size;
		for (int i = 0; i < aaSize; ++i) {
			curDeg = aaElems[i].degs;
			mpz_set(retPoly->coefs[curDeg], aaElems[i].coef);
		}
	}

	return retPoly;
}

AltArrZ_t* convertToAltArrZ_DUZP(const DUZP_t* poly) {
	if (poly == NULL) {
		return NULL;
	}

	int polySize = 0;
	for (int i = 0; i <= poly->lt; ++i){
		if (mpz_sgn(poly->coefs[i]) != 0) {
			++polySize;
		}
	}

	AltArrZ_t* retPoly = makePolynomial_AAZ(polySize, 1);
	polySize = 0;
	for (int i = poly->lt; i >= 0; --i){
		if (mpz_sgn(poly->coefs[i]) != 0) {
			mpz_init(retPoly->elems[polySize].coef);
			mpz_set(retPoly->elems[polySize].coef, poly->coefs[i]);
			retPoly->elems[polySize].degs = i;
			++polySize;
		}
	}
	if (polySize) {
		retPoly->size = polySize;
	}

	return retPoly;
}

DUZP_t* modularResultant_DUZP (const DUZP_t* a, const DUZP_t* b, mpz_t uBound, int deterministic, int hgcd)
{
	// corner cases
	if (isZero_DUZP (a)) {
		if (isZero_DUZP (b) || b->lt) {
			return NULL;
		} else {
			return makeConstPolynomial_DUZP (1, 1l);
		}
	} else if (isZero_DUZP (b)) {
		if (a->lt) {
			return NULL;
		} else {
			return makeConstPolynomial_DUZP (1, 1l);
		}
	}

	if (!a->lt && !b->lt) {
		return makeConstPolynomial_DUZP (1, 1l);
	}

	mpz_t bound;
	mpz_init (bound);
	if (deterministic) {
		// coef bound
		if (uBound == NULL || mpz_sgn (uBound) <= 0) {
			mpz_t* coefs = a->coefs;
			mpz_t maxA, maxB;
			mpz_t n1, m1;
			// compute Norm_max (f) and Norm_max (g)
			mpz_init (maxA);
			mpz_init (maxB);
			for (polysize_t i = 0; i <= a->lt; i++) {
				if (mpz_cmp (coefs[i], maxA) > 0)
					mpz_set (maxA, coefs[i]);
			}
			coefs = b->coefs;
			for (polysize_t i = 0; i <= b->lt; i++) {
				if (mpz_cmp (coefs[i], maxB) > 0)
					mpz_set (maxB, coefs[i]);
			}
			// compute resultant-bound
			mpz_pow_ui (maxA, maxA, (unsigned long) b->lt); // maxA^deg(b)
			mpz_pow_ui (maxB, maxB, (unsigned long) a->lt); // maxB^deg(a)
			mpz_init_set_si (n1, a->lt + 1); // n1 = deg(a) + 1
			mpz_init_set_si (m1, b->lt + 1); // m1 = deg(b) + 1
			mpz_pow_ui (n1, n1, (b->lt)>>1); // n1 = n1^(deg(b)/2)
			mpz_pow_ui (m1, m1, (a->lt)>>1); // m1 = m1^(deg(a)/2)
			mpz_set (bound, n1);
			mpz_mul (bound, bound, m1);
			mpz_mul (bound, bound, maxA);
			mpz_mul (bound, bound, maxB);
			// mpz_mul_si (bound, bound, 2); // TODO: see page 112 MCA

			if (uBound == NULL) {
				mpz_init_set (uBound, bound);
			} else {
				mpz_clear (uBound);
				mpz_init_set (uBound, bound);
			}
			mpz_clears (maxA, maxB, n1, m1, NULL);
		} else {
			mpz_set (bound, uBound);
		}
	}


	DUZP_t* res_work;
	mpz_t  res_work_coef, res_work_prev;
	mpz_inits (res_work_coef, res_work_prev, NULL);

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, NULL);

	int primeIdx = 0;
	Prime_ptr Pptr[1];

	duspoly_t* modA = makePolynomial_spX (a->lt + 1);
	duspoly_t* modB = makePolynomial_spX (b->lt + 1);
	duspoly_t* modRes= NULL;
	unsigned long res_out = 0;

	int fibs[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584}; // using Fib for determinisitic=0
	int fibIdx=0, nFibs=17;
	int isDone = 0;

	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		*Pptr = prime64_ptr[primeIdx];
		if (mpz_divisible_ui_p (a->coefs[a->lt], (unsigned long) Pptr->prime) ||
		    mpz_divisible_ui_p (b->coefs[b->lt], (unsigned long) Pptr->prime) ) {
				// bad prime
				// fprintf (stderr, "bad prime: %lld\n", (&prime64_ptr[primeIdx])->prime);
				continue;
			}

		// if good prime:
		// convert modular images
		convertToPrimeField_DUZP_inp (a, Pptr, &modA);
		convertToPrimeField_DUZP_inp (b, Pptr, &modB);
		// call resultant mod this prime
		if (!hgcd) {
			sylvResultantInForm_spX (modA, modB, &modRes, Pptr);
		} else {
			hgcdResultantInForm_spX (modA, modB, &modRes, Pptr);
		}
		// @note modRes != NULL for all modA, modB \in DUZP

		mpz_set_si (mpz_pr, Pptr->prime);
		mpz_gcdext (g,s,t, mpz_pr, m);

		mpz_mul (newm, m, mpz_pr);
		mpz_sub_ui (halfm, newm, 1ul);
		mpz_fdiv_q_2exp (halfm, halfm, 1ul);

		mpz_set (res_work_prev, res_work_coef);
		res_out = (unsigned long) smallprimefield_convert_out (modRes->elems[0], Pptr);

		// MCA-CRT (not-efficient!!)
		// mpz_gcdext (g, s, t, divm, mpz_pr);
		// mpz_mul (vs, mpz_coef, s);
		// mpz_mod (cc, vs, mpz_pr);
		// mpz_mul (ccdivm, cc, divm);
		// mpz_add (tsum, tsum, ccdivm);

		// optimized-CRT (Garner's algorithm)
		mpz_sub_ui (res_work_coef, res_work_coef, res_out);
		mpz_mul (res_work_coef, res_work_coef, s);
		mpz_mul (res_work_coef, res_work_coef, mpz_pr);
		mpz_add_ui (res_work_coef, res_work_coef, res_out);
		mpz_mod (res_work_coef, res_work_coef, newm);

		if (mpz_cmp (res_work_coef, halfm) > 0) {
			mpz_sub (res_work_coef, res_work_coef, newm);
		}

		freePolynomial_spX (&modRes);

		if (!deterministic) {
			int doDivideTest = fibIdx < nFibs ? (primeIdx == fibs[fibIdx++]) : 1;
			// fprintf (stderr, "doDivideTest = %d\n", doDivideTest); // TEST
			// gmp_printf ("For prime_Idx=%d => res=%Zd\n", primeIdx, res_work_coef); // TEST
			if (doDivideTest || !mpz_cmp(res_work_coef, res_work_prev)) {
				isDone = 1;
			}
		} else {
			if (mpz_cmp (m, bound) > 0) {
				isDone = 1;
			}
		}
		if (isDone) {
			if (!mpz_sgn (res_work_coef)) {
				res_work = NULL;
			} else {
				res_work = makePolynomial_DUZP (1);
				mpz_init_set (res_work->coefs[0], res_work_coef);
				res_work->lt = 0;
			}

			freePolynomial_spX (&modA);
			freePolynomial_spX (&modB);
			mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_coef, res_work_prev, NULL);
			return res_work;
		}
		mpz_set (m, newm);
	}

	freePolynomial_spX (&modA);
	freePolynomial_spX (&modB);
	mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_coef, res_work_prev, NULL);

	// all primes failed!!
	// fprintf (stderr, "all primes failed\n"); // TEST
	return NULL;

}

DUZP_t** modularSubresultantChain_DUZP (const DUZP_t* a, const DUZP_t* b, int* chain_size, mpz_t uBound, int deterministic)
{
	// corner cases
	DUZP_t** subres;
	if (isZero_DUZP (a)) {
		if (isZero_DUZP (b) || b->lt) {
			return NULL;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	} else if (isZero_DUZP (b)) {
		if (a->lt) {
			return NULL;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	}

	if (!a->lt) {
		if (!b->lt) {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = modularResultant_DUZP (a, b, uBound, deterministic, 0);
			*chain_size = 1;
			return subres;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	} else if (!b->lt) {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
	}

	mpz_t bound;
	mpz_init (bound);
	if (deterministic) {
		if (uBound == NULL || mpz_sgn (uBound) <= 0) {
			mpz_t* coefs = a->coefs;
			mpz_t maxA, maxB;
			mpz_t n1, m1;
			// compute Norm_max (f) and Norm_max (g)
			mpz_init (maxA);
			mpz_init (maxB);
			for (polysize_t i = 0; i <= a->lt; i++) {
				if (mpz_cmp (coefs[i], maxA) > 0)
					mpz_set (maxA, coefs[i]);
			}
			coefs = b->coefs;
			for (polysize_t i = 0; i <= b->lt; i++) {
				if (mpz_cmp (coefs[i], maxB) > 0)
					mpz_set (maxB, coefs[i]);
			}
			// compute resultant-bound
			mpz_pow_ui (maxA, maxA, (unsigned long) b->lt); // maxA^deg(b)
			mpz_pow_ui (maxB, maxB, (unsigned long) a->lt); // maxB^deg(a)
			mpz_init_set_si (n1, a->lt + 1); // n1 = deg(a) + 1
			mpz_init_set_si (m1, b->lt + 1); // m1 = deg(b) + 1
			mpz_pow_ui (n1, n1, (b->lt)>>1); // n1 = n1^(deg(b)/2)
			mpz_pow_ui (m1, m1, (a->lt)>>1); // m1 = m1^(deg(a)/2)
			mpz_set (bound, n1);
			mpz_mul (bound, bound, m1);
			mpz_mul (bound, bound, maxA);
			mpz_mul (bound, bound, maxB);
			// mpz_mul_si (bound, bound, 2); // TODO: see page 112 MCA

			if (uBound == NULL) {
				mpz_init_set (uBound, bound);
			} else {
				mpz_clear (uBound);
				mpz_init_set (uBound, bound);
			}
			mpz_clears (maxA, maxB, n1, m1, NULL);
		} else {
			mpz_set (bound, uBound);
		}
	}

	polysize_t min_deg = MIN_spX (a->lt, b->lt);
	DUZP_t** subres_work = (DUZP_t**) malloc (min_deg * sizeof(DUZP_t*));
	for (int k = 0; k < min_deg; k++) {
		subres_work[k] = makePolynomial_DUZP (k+1);
		for (int i = 0; i <= k; i++) {
			mpz_init (subres_work[k]->coefs[i]);
		}
	}
	// fprintf (stderr, "create subres_work\n"); //TEST
	mpz_t res_work_prev;
	mpz_init (res_work_prev);
	DUZP_t** subres_work_prev = (DUZP_t**) malloc (min_deg * sizeof(DUZP_t*));
	for (int k = 0; k < min_deg; k++) {
		subres_work_prev[k] = makePolynomial_DUZP (k+1);
		for (int i = 0; i <= k; i++) {
			mpz_init (subres_work_prev[k]->coefs[i]);
		}
	}


	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, NULL);

	int primeIdx = 0;
	Prime_ptr Pptr[1];

	duspoly_t* modA = makePolynomial_spX (a->lt + 1);
	duspoly_t* modB = makePolynomial_spX (b->lt + 1);
	duspolysA_t* modSubres= NULL;
	unsigned long coef_out = 0;
	polysize_t mod_chain_size = 0;
	polysize_t mod_chain_size_prev = 0;

	int fibs[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584}; // using Fib for determinisitic=0
	int fibIdx=0, nFibs=17;
	int isNULL = 0;
	int isDone = 0;
	mpz_t tc;

	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		*Pptr = prime64_ptr[primeIdx];

		mpz_init(tc);
		mpz_mod_ui(tc, a->coefs[a->lt],  (unsigned long) Pptr->prime);
		if(!smallprimefield_convert_in(mpz_get_si (tc), Pptr)) {
			mpz_clear(tc);
			continue;
		}
		mpz_mod_ui(tc, b->coefs[b->lt],  (unsigned long) Pptr->prime);
		if(!smallprimefield_convert_in(mpz_get_si (tc), Pptr)) {
			mpz_clear(tc);
			continue;
		}
		mpz_clear(tc);

		// if good prime:
		// convert modular images
		convertToPrimeField_DUZP_inp (a, Pptr, &modA);
		convertToPrimeField_DUZP_inp (b, Pptr, &modB);
		// call subresultant mod this prime
		subresultantChainInForm_spX (modA, modB, &modSubres, &mod_chain_size, Pptr);

		// fprintf (stderr, "subresultant called for pr_idx=%d\n", primeIdx); //TEST
		// for (int kk = 0; kk < mod_chain_size; kk++) { // TEST
		// 	fprintf (stderr, "modsubres[%d]=", kk);
		// 	printPolynomialOutForm_spX (modSubres->polys[kk], Pptr);
		// } // TEST

		if (mod_chain_size_prev && mod_chain_size_prev != mod_chain_size) {
			// bad prime
			fprintf (stderr, "bad prime (mod_chain_size is different: %ld != %ld): %lld\n", mod_chain_size, mod_chain_size_prev, (&prime64_ptr[primeIdx])->prime);
			exit(1);
		} else if (!mod_chain_size_prev) {
			mod_chain_size_prev = mod_chain_size;
		}

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);

		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// mpz_set (res_work_prev, subres_work[0]->coefs[0]);
		// subres_work_prev = subres_work
		for (int k = 0; k < min_deg; k++) {
			for (int i = 0; i <= subres_work[k]->lt; i++) {
				mpz_set (subres_work_prev[k]->coefs[i], subres_work[k]->coefs[i]);
			}
			subres_work_prev[k]->lt = subres_work[k]->lt;
		}

		for (int k = 0; k < mod_chain_size; k++) {
			if (isZero_spX (modSubres->polys[k])) {
				continue;
			}
			polysize_t lt = modSubres->polys[k]->lt;
			subres_work[lt]->lt = lt; // update lt
			for (int i = 0; i <= lt; ++i) {
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->polys[k]->elems[i], Pptr);
				// optimized-CRT (Garner's algorithm)
				mpz_sub_ui (subres_work[lt]->coefs[i], subres_work[lt]->coefs[i], coef_out);
				mpz_mul (subres_work[lt]->coefs[i], subres_work[lt]->coefs[i], s);
				mpz_mul (subres_work[lt]->coefs[i], subres_work[lt]->coefs[i], mpz_pr);
				mpz_add_ui (subres_work[lt]->coefs[i], subres_work[lt]->coefs[i], coef_out);
				mpz_mod (subres_work[lt]->coefs[i], subres_work[lt]->coefs[i], newm);

				if (mpz_cmp (subres_work[lt]->coefs[i], halfm) > 0) {
					mpz_sub(subres_work[lt]->coefs[i], subres_work[lt]->coefs[i], newm);
				}
			}
		}

		freePolysA_spX (modSubres);
		modSubres=NULL;

		if (!deterministic) {
			int doDivideTest = fibIdx < nFibs ? (primeIdx == fibs[fibIdx++]) : 1;
			// fprintf (stderr, "doDivideTest = %d\n", doDivideTest); // TEST
			// gmp_printf ("For prime_Idx=%d => res=%Zd\n", primeIdx, subres_work[0]->coefs[0]); // TEST
			// subres_work == subres_work_prev?
			int isSubresEq = 1;
			for (int k = 0; k < min_deg; k++){
				if (subres_work_prev[k]->lt != subres_work[k]->lt) {
					isSubresEq = 0;
					break;
				}
				for (int i = 0; i <= subres_work[k]->lt; i++) {
					if (mpz_cmp(subres_work_prev[k]->coefs[i], subres_work[k]->coefs[i])) {
						isSubresEq = 0;
						break;
					}
				}
				if (!isSubresEq) {
					break;
				}
			}

			if (doDivideTest || isSubresEq) { // !mpz_cmp(subres_work[0]->coefs[0], res_work_prev)
				isDone = 1;
			}
		} else {
			if (mpz_cmp (m, bound) > 0) {
				isDone = 1;
			}
		}
		if (isDone) {
			if (!mpz_sgn (subres_work[0]->coefs[0])) {
				freePolynomial_DUZP (subres_work[0]);
				subres_work[0] = NULL;
			}
			for (int k = 1; k < min_deg; k++) {
				if (!subres_work[k]->lt) {
					subres_work[k]->lt = k;
					freePolynomial_DUZP (subres_work[k]);
					subres_work[k] = NULL;
				}
			}

			freePolynomial_spX (&modA);
			freePolynomial_spX (&modB);
			mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_prev, NULL);
			*chain_size = min_deg;
			return subres_work;
		}

		mpz_set (m, newm);
	}

	freePolynomial_spX (&modA);
	freePolynomial_spX (&modB);
	mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_prev, NULL);
	*chain_size = min_deg;

	// all primes failed!! TODO:...
	// fprintf (stderr, "all primes failed\n"); // TEST
	return subres_work;
}

DUZP_t** modularSubresultantAtDeg_DUZP (const DUZP_t* a, const DUZP_t* b, int k, int* chain_size, mpz_t uBound, int deterministic, int hgcd)
{

	// char sym[1] = {'x'}; // TEST
	// fprintf (stderr, "In modularSubresultantAtDeg_DUZP: a := ");
	// printPoly_DUZP (a, sym); // TEST
	// fprintf (stderr, "In modularSubresultantAtDeg_DUZP: b := ");
	// printPoly_DUZP (b, sym); // TEST

	// corner cases
	DUZP_t** subres;
	if (isZero_DUZP (a)) {
		if (isZero_DUZP (b) || b->lt) {
			return NULL;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	} else if (isZero_DUZP (b)) {
		if (a->lt) {
			return NULL;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	}

	polysize_t min_deg = MIN_spX (a->lt, b->lt);
	if (k < 0 || k >= min_deg) {
		subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
		subres[0] = deepCopyPolynomial_DUZP (b);
		subres[1] = deepCopyPolynomial_DUZP (a);
		*chain_size = 2;
		return subres;
	}

	if (!k && min_deg < 2) {
		subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
		subres[0] = modularResultant_DUZP (a, b, uBound, deterministic, hgcd);
		subres[1] = (a->lt >= b->lt) ? deepCopyPolynomial_DUZP (b) : deepCopyPolynomial_DUZP (a);
		*chain_size = 2;
		return subres;
	}

	if (!hgcd) {
		int all_chain_size, all_idx=0;
		DUZP_t** all_subres = modularSubresultantChain_DUZP (a, b, &all_chain_size, uBound, deterministic);

		// char sym[1] = {'x'}; // TEST
		// for (int i = 0; i < all_chain_size; i++) { // TEST
		// 	fprintf (stderr, "all_subres[%d] := ", i);
		// 	printPoly_DUZP (all_subres[i], sym);
		// } // TEST

		subres = (DUZP_t**) malloc (sizeof (DUZP_t*));
		while (all_idx < all_chain_size) {
			if (isZero_DUZP (all_subres[all_idx]) || all_subres[all_idx]->lt < k) {
				all_idx++;
			} else {
				break;
			}
		}

		// fprintf (stderr, "all_chain_size = %d ||| all_idx = %d\n", all_chain_size, all_idx); // TEST

		if (all_idx < all_chain_size-1) {
			subres[0] = deepCopyPolynomial_DUZP (all_subres[all_idx+1]);
			subres[1] = deepCopyPolynomial_DUZP (all_subres[all_idx]);
			*chain_size = 2;

			for (int i = 0; i < all_chain_size; i++) {
				freePolynomial_DUZP (all_subres[i]);
			}
			free (all_subres);

			return subres;
		} else {
			fprintf (stderr, "DUZP Warning: In modularSubresultantAtDeg, all_idx exceeds valid range.\n");
		}

		for (int i = 0; i < all_chain_size; i++) {
			freePolynomial_DUZP (all_subres[i]);
		}
		free (all_subres);
		free (subres);
		*chain_size = 0;
		return NULL;
	}

	mpz_t bound;
	mpz_init(bound);
	if (deterministic) {
		if (uBound == NULL || mpz_sgn (uBound) <= 0) {
			mpz_t* coefs = a->coefs;
			mpz_t maxA, maxB;
			mpz_t n1, m1;
			// compute Norm_max (f) and Norm_max (g)
			mpz_init (maxA);
			mpz_init (maxB);
			for (polysize_t i = 0; i <= a->lt; i++) {
				if (mpz_cmp (coefs[i], maxA) > 0)
					mpz_set (maxA, coefs[i]);
			}
			coefs = b->coefs;
			for (polysize_t i = 0; i <= b->lt; i++) {
				if (mpz_cmp (coefs[i], maxB) > 0)
					mpz_set (maxB, coefs[i]);
			}
			// compute resultant-bound
			mpz_pow_ui (maxA, maxA, (unsigned long) b->lt); // maxA^deg(b)
			mpz_pow_ui (maxB, maxB, (unsigned long) a->lt); // maxB^deg(a)
			mpz_init_set_si (n1, a->lt + 1); // n1 = deg(a) + 1
			mpz_init_set_si (m1, b->lt + 1); // m1 = deg(b) + 1
			mpz_pow_ui (n1, n1, (b->lt)>>1); // n1 = n1^(deg(b)/2)
			mpz_pow_ui (m1, m1, (a->lt)>>1); // m1 = m1^(deg(a)/2)
			mpz_set (bound, n1);
			mpz_mul (bound, bound, m1);
			mpz_mul (bound, bound, maxA);
			mpz_mul (bound, bound, maxB);
			// mpz_mul_si (bound, bound, 2); // TODO: see page 112 MCA

			if (uBound == NULL) {
				mpz_init_set (uBound, bound);
			} else {
				mpz_clear (uBound);
				mpz_init_set (uBound, bound);
			}
			mpz_clears (maxA, maxB, n1, m1, NULL);
		} else {
			mpz_set (bound, uBound);
		}
	}

	// polysize_t min_deg = MIN_spX (a->lt, b->lt):
	DUZP_t** subres_work = (DUZP_t**) malloc (2 * sizeof(DUZP_t*));
	for (int i = 0; i < 2; i++) {
		subres_work[i] = makePolynomial_DUZP (k+10);
		for (int j = 0; j < k+10; j++) {
			mpz_init (subres_work[i]->coefs[j]);
		}
	}
	DUZP_t** subres_work_prev = (DUZP_t**) malloc (2* sizeof (DUZP_t*));
	for (int i = 0; i < 2; i++) {
		subres_work_prev[i] = makePolynomial_DUZP (k+10);
		for (int j = 0; j < k+10; j++) {
			mpz_init (subres_work_prev[i]->coefs[j]);
		}
	}

	// fprintf (stderr, "create subres_work\n"); //TEST
	mpz_t res_work_prev, res1_work_prev;
	mpz_inits (res_work_prev, res1_work_prev, NULL);

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	// fprintf (stderr, "start first init set...\n"); // TEST
	mpz_inits (g, s, t, mpz_pr, NULL);
	// fprintf (stderr, "finish first init set...\n"); // TEST

	int primeIdx = 0;
	Prime_ptr Pptr[1];

	mpz_t tc;
	mpz_t *tmp_coef;
	duspoly_t* modA = makePolynomial_spX (a->lt + 1);
	duspoly_t* modB = makePolynomial_spX (b->lt + 1);
	duspolysA_t* modSubres= NULL;
	unsigned long coef_out = 0;
	polysize_t mod_chain_size = 0;
	polysize_t mod_chain_size_prev = 0;
	polysize_t max_allocs;

	int fibs[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584}; // using Fib for determinisitic=0
	int fibIdx=0, nFibs=17;
	int isNULL = 0;
	int isDone = 0;

	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		*Pptr = prime64_ptr[primeIdx];

		mpz_init(tc);
		mpz_mod_ui(tc, a->coefs[a->lt],  (unsigned long) Pptr->prime);
		if(!smallprimefield_convert_in(mpz_get_si (tc), Pptr)) {
			mpz_clear(tc);
			continue;
		}
		mpz_mod_ui(tc, b->coefs[b->lt],  (unsigned long) Pptr->prime);
		if(!smallprimefield_convert_in(mpz_get_si (tc), Pptr)) {
			mpz_clear(tc);
			continue;
		}
		mpz_clear(tc);

		// if good prime:
		// convert modular images
		convertToPrimeField_DUZP_inp (a, Pptr, &modA);
		convertToPrimeField_DUZP_inp (b, Pptr, &modB);
		// call subresultant mod this prime
		hgcdSubResultantInFormA_spX (modA, modB, k, &modSubres, &mod_chain_size, Pptr);

		// fprintf (stderr, "In modularSubresultantAtDeg_DUZP: subresultant called for pr_idx=%d (mod_chain_size = %ld)\n", primeIdx, mod_chain_size); //TEST
		// for (int kk = 0; kk < 2; kk++) { // TEST
		// 	fprintf (stderr, "modularSubresultantAtDeg_DUZP: modsubres[%d]=", kk);
		// 	printPolynomialOutForm_spX (modSubres->polys[kk], Pptr);
		// } // TEST

		if (mod_chain_size_prev && mod_chain_size_prev != mod_chain_size) {
			// bad prime
			fprintf (stderr, "bad prime (mod_chain_size is different: %ld != %ld): %lld\n", mod_chain_size, mod_chain_size_prev, (&prime64_ptr[primeIdx])->prime);
			exit(1);
			// TODO: return 0;
		} else if (!mod_chain_size_prev) {
			mod_chain_size_prev = mod_chain_size;
		}

		max_allocs = -1;
		for(int i = 0; i < 2; i++) {
			if(modSubres->polys[mod_chain_size-1-i] != NULL &&
			subres_work[i]->alloc < modSubres->polys[mod_chain_size-1-i]->alloc &&
			max_allocs < modSubres->polys[mod_chain_size-1-i]->alloc) {
				max_allocs = modSubres->polys[mod_chain_size-1-i]->alloc;
			}
		}
		if (max_allocs != -1) {
			// TODO: not a good implementation! this case rarely happens though..
			// fprintf (stderr, "max_allocs = %ld\n", max_allocs);
			for (int i = 0; i < 2; i++) {
				// fprintf (stderr, "subres_work[%d]->alloc = %ld\n", i, subres_work[i]->alloc);
				tmp_coef = (mpz_t*) malloc (max_allocs * sizeof(mpz_t));
				for (int j = 0; j < subres_work[i]->alloc; j++) {
					mpz_init(tmp_coef[j]);
					mpz_set(tmp_coef[j], subres_work[i]->coefs[j]);
					mpz_clear(subres_work[i]->coefs[j]);
				}
				for (int j = subres_work[i]->alloc; j < max_allocs; j++) {
					mpz_init (tmp_coef[j]);
				}
				free(subres_work[i]->coefs);
				subres_work[i]->coefs = tmp_coef;
				tmp_coef = (mpz_t*) malloc (max_allocs * sizeof(mpz_t));
				for (int j = 0; j < subres_work_prev[i]->alloc; j++) {
					mpz_init(tmp_coef[j]);
					mpz_set(tmp_coef[j], subres_work_prev[i]->coefs[j]);
					mpz_clear(subres_work_prev[i]->coefs[j]);
				}
				for (int j = subres_work_prev[i]->alloc; j < max_allocs; j++) {
					mpz_init (tmp_coef[j]);
				}
				free(subres_work_prev[i]->coefs);
				subres_work_prev[i]->coefs = tmp_coef;
				subres_work[i]->alloc = max_allocs;
				subres_work_prev[i]->alloc = max_allocs;
			}
		}

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);

		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// mpz_set (res_work_prev, subres_work[0]->coefs[0]);
		// mpz_set (res1_work_prev, subres_work[1]->coefs[0]);
		for (int kk = 0; kk < 2; kk++) {
			for (int i = 0; i <= subres_work[kk]->lt; i++) {
				mpz_set (subres_work_prev[kk]->coefs[i], subres_work[kk]->coefs[i]);
				// gmp_fprintf (stderr, "subres_work_prev[%d]->coefs[%d] = %Zd\t subres_work[%d]->coefs[%d] = %Zd\n", kk, i, subres_work_prev[kk]->coefs[i], kk, i, subres_work[kk]->coefs[i]); // TEST
			}
		}

		for (int kk = 0; kk < mod_chain_size; kk++) {

			// fprintf (stderr, "[kk=%d] modSubres->poly[kk] := ", kk); // TEST
			// printPolynomialOutForm_spX (modSubres->polys[kk], Pptr); // TEST

			if (isZero_spX (modSubres->polys[kk])) {
				continue;
			}
			polysize_t lt = modSubres->polys[kk]->lt;
			polysize_t kk_idx = mod_chain_size-1-kk;
			subres_work[kk_idx]->lt = lt; // update lt
			for (int i = 0; i <= lt; ++i) {
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->polys[kk]->elems[i], Pptr);
				// optimized-CRT (Garner's algorithm)
				// gmp_fprintf (stderr, "subres_work[%d]->coefs[%d] = %Zd\t coef_out = %lu\n", kk_idx, i, subres_work[kk_idx]->coefs[i], coef_out); // TEST
				mpz_sub_ui (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], coef_out);
				mpz_mul (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], s);
				mpz_mul (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], mpz_pr);
				mpz_add_ui (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], coef_out);
				mpz_mod (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], newm);

				if (mpz_cmp (subres_work[kk_idx]->coefs[i], halfm) > 0) {
					mpz_sub(subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], newm);
				}
			}
		}

		freePolysA_spX (modSubres);
		modSubres=NULL;

		if (!deterministic) {
			int doDivideTest = fibIdx < nFibs ? (primeIdx == fibs[fibIdx++]) : 1;
			// fprintf (stderr, "doDivideTest = %d\n", doDivideTest); // TEST
			// gmp_printf ("For prime_Idx=%d => res=%Zd\n", primeIdx, subres_work[0]->coefs[0]); // TEST
			int EqualityTest = 1;
			for (int kk = 0; kk < 2; kk++) {
				for (int i = 0; i <= subres_work[kk]->lt; i++) {
					if (mpz_cmp(subres_work[kk]->coefs[i], subres_work_prev[kk]->coefs[i])) {
						EqualityTest = 0;
						break;
					}
				}
				if (!EqualityTest) {
					break;
				}
			}

			if (doDivideTest || EqualityTest) {
				isDone = 1;
			}
		} else {
			if (mpz_cmp (m, bound) > 0) {
				isDone = 1;
			}
		}

		// fprintf(stderr, "isDone := %d\n",isDone); // TEST

		if (isDone) {
			if (!subres_work[0]->lt && !mpz_sgn (subres_work[0]->coefs[0])) {
				freePolynomial_DUZP (subres_work[0]);
				subres_work[0] = NULL;
			}
			if (!subres_work[1]->lt && !mpz_sgn (subres_work[1]->coefs[0])) {
				freePolynomial_DUZP (subres_work[1]);
				subres_work[1] = NULL;
			}

			freePolynomial_spX (&modA);
			freePolynomial_spX (&modB);
			mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_prev, res1_work_prev, NULL);
			*chain_size = 2;
			return subres_work;
		}

		mpz_set (m, newm);
	}

	freePolynomial_spX (&modA);
	freePolynomial_spX (&modB);
	mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_prev, res1_work_prev, NULL);
	*chain_size = 2;

	// all primes failed!! TODO:...
	// fprintf (stderr, "all primes failed\n"); // TEST
	return subres_work;
}

DUZP_t** regularGCDUnivariateSpeculativeSRC_DUZP (const DUZP_t* a, const DUZP_t* b, int k, int* chain_size, 
											polysize_t *degs, mpz_t uBound, specAQRArray_spX_t **uspecInfoArray, 
											int *info_size, int results_mode)
{
	// corner cases
	DUZP_t** subres;
	if (isZero_DUZP (a)) {
		if (isZero_DUZP (b) || b->lt) {
			return NULL;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	} else if (isZero_DUZP (b)) {
		if (a->lt) {
			return NULL;
		} else {
			subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
			subres[0] = makeConstPolynomial_DUZP (1, 1l);
			*chain_size = 1;
			return 	subres;
		}
	}

	polysize_t min_deg = MIN_spX (a->lt, b->lt);
	if (k < 0 || k >= min_deg) {
		subres = (DUZP_t**) malloc (sizeof(DUZP_t*));
		subres[0] = deepCopyPolynomial_DUZP (b);
		subres[1] = deepCopyPolynomial_DUZP (a);
		*chain_size = 2;
		return subres;
	}

	mpz_t bound;
	mpz_init(bound);

	// polysize_t min_deg = MIN_spX (a->lt, b->lt):
	DUZP_t** subres_work = (DUZP_t**) malloc (2 * sizeof(DUZP_t*));
	for (int i = 0; i < 2; i++) {
		subres_work[i] = makePolynomial_DUZP (k+10);
		for (int j = 0; j < k+10; j++) {
			mpz_init (subres_work[i]->coefs[j]);
		}
	}
	DUZP_t** subres_work_prev = (DUZP_t**) malloc (2* sizeof (DUZP_t*));
	for (int i = 0; i < 2; i++) {
		subres_work_prev[i] = makePolynomial_DUZP (k+10);
		for (int j = 0; j < k+10; j++) {
			mpz_init (subres_work_prev[i]->coefs[j]);
		}
	}

	// fprintf (stderr, "create subres_work\n"); //TEST
	mpz_t res_work_prev, res1_work_prev;
	mpz_inits (res_work_prev, res1_work_prev, NULL);

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	// fprintf (stderr, "start first init set...\n"); // TEST
	mpz_inits (g, s, t, mpz_pr, NULL);
	// fprintf (stderr, "finish first init set...\n"); // TEST

	int primeIdx = 0;
	Prime_ptr Pptr[1];
	
	mpz_t tc;
	mpz_t *tmp_coef;
	duspoly_t* modA = makePolynomial_spX (a->lt + 1);
	duspoly_t* modB = makePolynomial_spX (b->lt + 1);
	duspolysA_t* modSubres= NULL;
	unsigned long coef_out = 0;
	polysize_t mod_chain_size = 0, 
			mod_chain_size_prev = 0,
			max_allocs, td;
	int isNULL = 0;

	int infoIdx = 0,
		infoAlloc = 25, 
		isInfoNull = 0;
	specAQR_spX_t ** uspecInfo = NULL;
	if (*info_size < 1 || uspecInfoArray == NULL) {
		uspecInfo = (specAQR_spX_t**) malloc (infoAlloc * sizeof(specAQR_spX_t*));
		for (int i = 0; i < infoAlloc; ++i) {
			uspecInfo[i] = (specAQR_spX_t *) malloc (sizeof(specAQR_spX_t));
		}
		isInfoNull = 1;
	} else {
		infoAlloc = *info_size;
		uspecInfo = (*uspecInfoArray)->uspecArray; 
	}

	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		// fprintf(stderr, "primeIdx = %d\n", primeIdx);
		*Pptr = prime64_ptr[primeIdx];

		mpz_init(tc);
		mpz_mod_ui(tc, a->coefs[a->lt],  (unsigned long) Pptr->prime);
		if(!smallprimefield_convert_in(mpz_get_si (tc), Pptr)) {
			mpz_clear(tc);
			continue;
		}
		mpz_mod_ui(tc, b->coefs[b->lt],  (unsigned long) Pptr->prime);
		if(!smallprimefield_convert_in(mpz_get_si (tc), Pptr)) {
			mpz_clear(tc);
			continue;
		}
		mpz_clear(tc);

		// if good prime:
		// convert modular images
		convertToPrimeField_DUZP_inp (a, Pptr, &modA);
		convertToPrimeField_DUZP_inp (b, Pptr, &modB);
		// call subresultant mod this prime

		if (isInfoNull) {
			// fprintf(stderr, "Initialize A, Q, R for primeIdx = %d\n", primeIdx);
			uspecInfo[infoIdx]->A = createSpecA_spX (modA->alloc);
			// cerr << "A: size := " << A->size << "  alloc := " << A->alloc << endl;
			uspecInfo[infoIdx]->Q = createSpecQ_spX (modA->alloc);
			// cerr << "Q: size := " << Q->size << "  alloc := " << Q->alloc << endl;
			uspecInfo[infoIdx]->R = createSpecR_spX (modA->alloc);
			// cerr << "R: size := " << R->size << "  alloc := " << R->alloc << endl;
		}
		// fprintf(stderr, "primeIdx = %d, infoIdx = %d, isInfoNull = %d\n", primeIdx, infoIdx, isInfoNull);
		// hgcdSubResultantInFormA_spX (modA, modB, k, &modSubres, &mod_chain_size, Pptr);
		// fprintf(stderr, "[old] modSubres->size = %d mod_chain_size = %d\n", modSubres->size, mod_chain_size);
		// fprintf(stderr, "[old] modSubres->poly[0] = ");
		// printPolynomial_spX (modSubres->polys[0], Pptr);
		// fprintf(stderr, "[old] modSubres->poly[1] = ");
		// printPolynomial_spX (modSubres->polys[1], Pptr);
		// freePolysA_spX (modSubres); modSubres = NULL;
		regularGCDSpecSRC_spX (modA, modB, k, &modSubres, degs, results_mode, 
						uspecInfo[infoIdx]->A, uspecInfo[infoIdx]->Q, uspecInfo[infoIdx]->R, Pptr);
		// fprintf(stderr, "[new] modSubres->size = %d\n", modSubres->size);
		// fprintf(stderr, "[new] modSubres->poly[0] = ");
		// printPolynomial_spX (modSubres->polys[0], Pptr);
		// fprintf(stderr, "[new] modSubres->poly[1] = ");
		// printPolynomial_spX (modSubres->polys[1], Pptr);
		// fprintf(stderr, "regularGCDSpecSRC_spX... DONE  size=%d\n", modSubres->size);
		mod_chain_size = modSubres->size;
		infoIdx++;		

		if (infoIdx == infoAlloc) {
			// fprintf(stderr, "reallocate specInfo... \n");
			infoAlloc += 25; 
			specAQR_spX_t ** tmpSpecInfo = (specAQR_spX_t **) malloc (infoAlloc * sizeof(specAQR_spX_t *));
			for (int i = 0; i < infoIdx; ++i) {
				tmpSpecInfo[i] = uspecInfo[i];
			}
			for (int i = infoIdx; i < infoAlloc; ++i) {
				tmpSpecInfo[i] = (specAQR_spX_t*) malloc (sizeof(specAQR_spX_t));
			}
			free(uspecInfo);
			uspecInfo = tmpSpecInfo;
		}

		// fprintf (stderr, "In modularSubresultantAtDeg_DUZP: subresultant called for pr_idx=%d (mod_chain_size = %ld)\n", primeIdx, mod_chain_size); //TEST
		// for (int kk = 0; kk < 2; kk++) { // TEST
		// 	fprintf (stderr, "modularSubresultantAtDeg_DUZP: modsubres[%d]=", kk);
		// 	printPolynomialOutForm_spX (modSubres->polys[kk], Pptr);
		// } // TEST

		if (mod_chain_size_prev && mod_chain_size_prev != mod_chain_size) {
			// bad prime
			fprintf (stderr, "bad prime (mod_chain_size is different: %ld != %ld): %lld\n", mod_chain_size, mod_chain_size_prev, (&prime64_ptr[primeIdx])->prime);
			exit(1);
			// TODO: return 0;
		} else if (!mod_chain_size_prev) {
			mod_chain_size_prev = mod_chain_size;
		}

		max_allocs = -1;
		for(int i = 0; i < 2; i++) {
			if(modSubres->polys[mod_chain_size-1-i] != NULL &&
			subres_work[i]->alloc < modSubres->polys[mod_chain_size-1-i]->alloc &&
			max_allocs < modSubres->polys[mod_chain_size-1-i]->alloc) {
				max_allocs = modSubres->polys[mod_chain_size-1-i]->alloc;
			}
		}
		if (max_allocs != -1) {
			// TODO: not a good implementation! this case rarely happens though..
			// fprintf (stderr, "max_allocs = %ld\n", max_allocs);
			for (int i = 0; i < 2; i++) {
				// fprintf (stderr, "subres_work[%d]->alloc = %ld\n", i, subres_work[i]->alloc);
				tmp_coef = (mpz_t*) malloc (max_allocs * sizeof(mpz_t));
				for (int j = 0; j < subres_work[i]->alloc; j++) {
					mpz_init(tmp_coef[j]);
					mpz_set(tmp_coef[j], subres_work[i]->coefs[j]);
					mpz_clear(subres_work[i]->coefs[j]);
				}
				for (int j = subres_work[i]->alloc; j < max_allocs; j++) {
					mpz_init (tmp_coef[j]);
				}
				free(subres_work[i]->coefs);
				subres_work[i]->coefs = tmp_coef;
				tmp_coef = (mpz_t*) malloc (max_allocs * sizeof(mpz_t));
				for (int j = 0; j < subres_work_prev[i]->alloc; j++) {
					mpz_init(tmp_coef[j]);
					mpz_set(tmp_coef[j], subres_work_prev[i]->coefs[j]);
					mpz_clear(subres_work_prev[i]->coefs[j]);
				}
				for (int j = subres_work_prev[i]->alloc; j < max_allocs; j++) {
					mpz_init (tmp_coef[j]);
				}
				free(subres_work_prev[i]->coefs);
				subres_work_prev[i]->coefs = tmp_coef;
				subres_work[i]->alloc = max_allocs;
				subres_work_prev[i]->alloc = max_allocs;
			}
		}

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);

		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		for (int kk = 0; kk < 2; kk++) {
			for (int i = 0; i <= subres_work[kk]->lt; i++) {
				mpz_set (subres_work_prev[kk]->coefs[i], subres_work[kk]->coefs[i]);
			}
		}

		for (int kk = 0; kk < mod_chain_size; kk++) {
			if (isZero_spX (modSubres->polys[kk])) {
				continue;
			}
			polysize_t lt = modSubres->polys[kk]->lt;
			polysize_t kk_idx = mod_chain_size-1-kk;
			subres_work[kk_idx]->lt = lt; // update lt
			for (int i = 0; i <= lt; ++i) {
				coef_out = (unsigned long) smallprimefield_convert_out (modSubres->polys[kk]->elems[i], Pptr);
				// gmp_fprintf (stderr, "subres_work[%d]->coefs[%d] = %Zd\t coef_out = %lu\n", kk_idx, i, subres_work[kk_idx]->coefs[i], coef_out); // TEST
				mpz_sub_ui (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], coef_out);
				mpz_mul (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], s);
				mpz_mul (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], mpz_pr);
				mpz_add_ui (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], coef_out);
				mpz_mod (subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], newm);

				if (mpz_cmp (subres_work[kk_idx]->coefs[i], halfm) > 0) {
					mpz_sub(subres_work[kk_idx]->coefs[i], subres_work[kk_idx]->coefs[i], newm);
				}
			}
		}

		freePolysA_spX (modSubres);
		modSubres=NULL;

		int EqualityTest = 1;
		for (int kk = 0; kk < 2; kk++) {
			for (int i = 0; i <= subres_work[kk]->lt; i++) {
				if (mpz_cmp(subres_work[kk]->coefs[i], subres_work_prev[kk]->coefs[i])) {
					EqualityTest = 0;
					break;
				}
			}
			if (!EqualityTest) {
				break;
			}
		}

		if (EqualityTest) {

			if (*info_size < 0 || uspecInfoArray == NULL) {
				for (int i = 0; i < infoIdx; ++i) {
					freeSpecA_spX (uspecInfo[i]->A);
					freeSpecQ_spX (uspecInfo[i]->Q);
					freeSpecR_spX (uspecInfo[i]->R);
				}
				free(uspecInfo);
			} else if (*info_size < infoIdx) {
				if (*uspecInfoArray == NULL) {
					(*uspecInfoArray) = (specAQRArray_spX_t *) malloc (sizeof(specAQRArray_spX_t));
				}
				(*uspecInfoArray)->uspecArray = uspecInfo;
				(*uspecInfoArray)->size = infoIdx;
				(*info_size) = infoIdx;
			}

			if (!subres_work[0]->lt && !mpz_sgn (subres_work[0]->coefs[0])) {
				freePolynomial_DUZP (subres_work[0]);
				subres_work[0] = NULL;
			}
			if (!subres_work[1]->lt && !mpz_sgn (subres_work[1]->coefs[0])) {
				freePolynomial_DUZP (subres_work[1]);
				subres_work[1] = NULL;
			}

			freePolynomial_spX (&modA);
			freePolynomial_spX (&modB);
			mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_prev, res1_work_prev, NULL);
			*chain_size = 2;
			td = degs[0]; degs[0] = degs[1]; degs[1] = td;
			return subres_work;
		}

		mpz_set (m, newm);
	}

	freePolynomial_spX (&modA);
	freePolynomial_spX (&modB);
	mpz_clears (bound, m, halfm, newm, g, s, t, mpz_pr, res_work_prev, res1_work_prev, NULL);
	*chain_size = 2;
	td = degs[0]; degs[0] = degs[1]; degs[1] = td;
	// all primes failed!! TODO:...
	// fprintf (stderr, "all primes failed\n"); // TEST
	return subres_work;
}



duspolysA_t* _transposePolysA_and_ConvertOut_spX (duspolysA_t* aa, polysize_t max_deg)
{
    if (aa == NULL || aa->size == 0) {
        return NULL;
    }
    // aa: Mat(size, max_deg+1)
    int size = aa->size;
    polysize_t t;
    // aa_T : Mat(max_deg+1, size)
    duspolysA_t* aa_T = makePolysAwithMaxAlloc_spX (max_deg+1, size);
	Prime_ptr Pptr[1];
    for (int j=0; j<=max_deg; j++) { // traverse monomials
        for (int i=0; i<size; i++) { // traverse polynomials
			if (aa->polys[i] != NULL && aa->polys[i]->lt >= j) {
				// TODO: check the correctness:
				aa_T->polys[j]->elems[i] = smallprimefield_convert_out (aa->polys[i]->elems[j], &prime64_ptr[i]);
			} else {
				aa_T->polys[j]->elems[i] = 0;
			}

        }
        aa_T->polys[j]->lt = size-1;

        if (!aa_T->polys[j]->elems[size-1]) {
            t = size - 2;
            while (!aa_T->polys[j]->elems[t]) {
                t--;
            }
            if (t < 0) {
                freePolynomial_spX (&aa_T->polys[j]);
            } else {
                aa_T->polys[j]->lt = t;
                // TODO: need reallocation?
            }
        }
    }

    return aa_T;
}
