#include "IntegerPolynomial/DUZP_Support_Factoring.hpp"
#include "IntegerPolynomial/NTL_mpz_DUZP_conversion.hpp"
#include "IntegerPolynomial/DUZP_Support.h"
#include "ModularPolynomial/DUSP_NTL_Support.h"
#include "NTL/vec_vec_long.h"
#include "NTL/ZZ.h"
#include "NTL/ZZXFactoring.h"
#include "NTL/mat_ZZ.h"
#include "NTL/LLL.h"
#include "Utils/MacroHelpers.h"
#include <math.h>

using namespace DUZP::Factoring;

#define DO_NTL_UNIVAR_FACT 1

char x[1] = {'x'};

/**
 * A structure for storing information about
 * the modular factorization process:
 *
 * Pptr: Prime_ptr for the best prime
 *       (smallest number of modular factors).
 * num_primes: The number of primes found that
 *             permit modular factorization.
 * possibleDegrees: bit vector of possible
 *                  degrees of the rational
 *                  factors (used for
 *                  irreducibility detection).
 * p: vector of the primes that permit modular
 *    factorization.
 * pattern: the degree patterns obtained after
 *          distinct degree factorization for
 *          the corresponding primes, stored as
 *          a vector of vectors of frequencies
 *          of degrees (given by the index).
 * s: an NTL class for generating a sequence of
 *    prime numbers.
 */
typedef struct factoring_info {
	Prime_ptr* Pptr;
	long num_primes;
	mpz_t possibleDegrees;
	NTL::vec_long p;
	NTL::vec_vec_long pattern;
	NTL::PrimeSeq s;
} factoring_info_t;

/**
 * Initialize a factoring_info_t structure.
 */
void factoring_info_t_init(factoring_info_t* fi) {
	mpz_init(fi->possibleDegrees);
	fi->Pptr = NULL;
}

/**
 * Free a factoring_info_t structure.
 */
void factoring_info_t_free(factoring_info_t* fi) {
	mpz_clear(fi->possibleDegrees);
	if (fi->Pptr != NULL)
		free(fi->Pptr);
}

/**
 * Print an NTL multiprecision integer matrix.
 */
void print_NTL_mat_ZZ(NTL::mat_ZZ& M) {

	for (int i=1; i<=M.NumRows(); ++i) {
		for (int j=1; j<=M.NumCols(); ++j)
			std::cerr << M(i,j) << "\t";
		std::cerr << std::endl;
	}

}

/**
 * Print a vector of mpz_t.
 **/
void DUZP::Factoring::print_vec_mpz(const vec_mpz_t a) {
	for (int i=0; i<a.size; ++i)
		gmp_fprintf(stderr,"vec[%d]=%Zd\n",i,a.array[i]);
}

/**
 * Print a vector of DUZP_t.
 **/
void DUZP::Factoring::print_vec_DUZP(const vec_DUZP_t a) {
	for (int i=0; i<a.size; ++i)
		printPoly_DUZP(a.polys[i],x);
}

/**
 * Print a vector of duspoly_t.
 **/
void DUZP::Factoring::print_vec_duspoly(const vec_duspoly_t a) {
	for (int i=0; i<a.size; ++i)
		printPolynomialOutForm_spX(a.polys[i],a.Pptr);
}

/**
 * Print a vector of DUZP, long pairs.
 **/
void DUZP::Factoring::print_vec_DUZP_long(const vec_DUZP_long_t a) {
	for (int i=0; i<a.size; ++i) {
		fprintf(stderr,"b = %ld, a = ",a.pairs[i].b);
		printPoly_DUZP(a.pairs[i].a,x);
	}
}

/**
 * print an array of mpz_t.
 **/
void DUZP::Factoring::print_mpz_array(const mpz_t* a, long len) {
	for (int i=0; i<len; ++i)
		gmp_fprintf(stderr,"arr[%d]=%Zd\n",i,a[i]);
}

/**
 * Print an array of long_t.
 **/
void DUZP::Factoring::print_long_t_array(const long_t* a, long len) {
	for (int i=0; i<len; ++i)
		gmp_fprintf(stderr,"arr[%d]=%Ld\n",i,a[i]);
}

/**
 * Print a vector of long (currently takes an NTL
 * type, but will eventually be a BPAS vec_long.
 **/
void print_vec_long(const NTL::vec_long& a) {
	for (int i=1; i<=a.length(); ++i)
		gmp_fprintf(stderr,"vec[%d]=%ld\n",i,a(i));
}

/**
 * Print a vector of vectors of mpz_t.
 **/
void DUZP::Factoring::print_vec_vec(const vec_vec_mpz_t* vec) {
	for (int i=0; i<vec->length; ++i)
		for (int j=0; j<vec->depth; ++j)
			gmp_fprintf(stderr,"vec[%d][%d] = %Zd\n",i,j,vec->array[i][j]);
}

/**
 * Set an mpz_t using a signed long long int.
 **/
void mpz_set_sll(mpz_t n, long long sll) {

    mpz_set_si(n,(int)(sll >> 32));
    mpz_mul_2exp(n,n,32);
    mpz_add_ui(n,n,(unsigned int)sll);
}

/**
 * Multiply an mpz_t by a signed long long int.
 **/
void mpz_mul_sll(mpz_t r, mpz_t n, long long sll) {

	mpz_t tmp;
	mpz_init(tmp);
	mpz_set_sll(tmp,sll);
	mpz_mul(r,n,tmp);
	mpz_clear(tmp);
}

/**
 * Set an mpz_t using an unsigned long long int.
 **/
void mpz_set_ull(mpz_t n, unsigned long long ull) {

    mpz_set_ui(n,(unsigned int)(ull >> 32));
    mpz_mul_2exp(n,n,32);
    mpz_add_ui(n,n,(unsigned int)ull);
}

/**
 * Get an unsigned long long int from an mpz_t.
 **/
unsigned long long mpz_get_ull(mpz_t n) {

    unsigned int lo, hi;
    mpz_t tmp;

    mpz_init(tmp);
    mpz_mod_2exp(tmp,n,64);   /* tmp = (lower 64 bits of n) */

    lo = mpz_get_ui(tmp);     /* lo = tmp & 0xffffffff */
    mpz_div_2exp(tmp,tmp,32); /* tmp >>= 32 */
    hi = mpz_get_ui(tmp);     /* hi = tmp & 0xffffffff */

    mpz_clear(tmp);

    return (((unsigned long long)hi) << 32) + lo;
}

/**
 * Get a signed long long int from an mpz_t.
 **/
long long mpz_get_sll(mpz_t n) {

	int negate = 0;
	if (mpz_cmp_ui(n,0) < 0) {
		negate = 1;
		mpz_neg(n,n);
	}

	long long out = (long long)mpz_get_ull(n);
	if (negate)
		 out = -out;
    return out;
}

/**
 * Push a DUZP_t to a vector of DUZP_t.
 **/
void pushFactor(vec_DUZP_t* factors, DUZP_t* f) {

	fprintf(stderr,"[PushFactor] adding new factor:\n");
	if (factors->alloc == factors->size) {
		factors->alloc = factors->alloc*2;
		factors->polys = (DUZP_t**) realloc(factors->polys,sizeof(DUZP_t*)*factors->alloc);
	}
	factors->polys[factors->size] = deepCopyPolynomial_DUZP(f);
	factors->size++;
}

/**
 * Initialize a vector of DUZP_t (allocation of 1, size 0).
 **/
void DUZP::Factoring::vec_DUZP_t_init(vec_DUZP_t* a) {

	a->polys = (DUZP_t**) malloc(sizeof(DUZP_t*));
	a->polys[0] = makeConstPolynomial_DUZP(1,0);
	a->alloc = 1;
	a->size = 0;

}

/**
 * Initialize a vector of DUZP_t (allocation alloc; size 0).
 **/
void DUZP::Factoring::vec_DUZP_t_init2(vec_DUZP_t* a, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->polys = (DUZP_t**) malloc(sizeof(DUZP_t*)*alloc);
	for (int i=0; i<alloc; ++i)
		a->polys[i] = makeConstPolynomial_DUZP(1,0);
	a->alloc = alloc;
	a->size = 0;

}

/**
 * Initialize a vector of DUZP_t (allocation alloc; size alloc).
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_DUZP_t_init_set(vec_DUZP_t* a, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->polys = (DUZP_t**) malloc(sizeof(DUZP_t*)*alloc);
	for (int i=0; i<alloc; ++i)
		a->polys[i] = makeConstPolynomial_DUZP(1,0);
	a->alloc = alloc;
	a->size = alloc;

}

/**
 * Set the size of a vector of DUZP_t; if alloc < size, then
 * the vector will be realloc'd to allow the size to be set.
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_DUZP_t_set(vec_DUZP_t* a, long size) {

	if (size < 0) {
		fprintf(stderr,"BPAS, error: vector size must be non-negative.\n");
		exit(1);
	}

	if (size > a->alloc) {
		a->polys = (DUZP_t**) realloc(a->polys,sizeof(DUZP_t*)*size);
		for (int i=a->alloc; i<size; ++i)
			a->polys[i] = makeConstPolynomial_DUZP(1,0);
		a->alloc = size;
	}
	a->size = size;

}

/**
 * Clear a vector of DUZP_t. The size will be set to
 * zero and the alloc to 1.
 **/
void DUZP::Factoring::vec_DUZP_t_clear(vec_DUZP_t* a) {

	if (a->alloc == 0) {
		fprintf(stderr,"BPAS, error: vector must be initialized to clear it.\n");
		exit(1);
	}

	for (int i=0; i<a->alloc; ++i)
		freePolynomial_DUZP(a->polys[i]);
	a->polys = (DUZP_t**) realloc(a->polys,sizeof(DUZP_t*));
	a->polys[0] = makeConstPolynomial_DUZP(1,0);
	a->alloc = 1;
	a->size = 0;

}

/**
 * Push a DUZP_t onto the end of a vector of DUZP_t.
 *
 * This function automatically expands the vector (by doubling)
 * if the allocated space is full.
 **/
void DUZP::Factoring::vec_DUZP_t_push(vec_DUZP_t* a, const DUZP_t* b) {

	if (a->size == a->alloc) {
		a->polys = (DUZP_t**) realloc(a->polys,sizeof(DUZP_t*)*2*a->alloc);
		for (int i=a->alloc; i<2*a->alloc; ++i) {
			a->polys[i] = makeConstPolynomial_DUZP(1,0);
		}
		a->alloc *= 2;
	}
	DUZP_t* c = deepCopyPolynomial_DUZP(b);
	freePolynomial_DUZP(a->polys[a->size]);
	a->polys[a->size] = c;
	a->size = a->size+1;

}

/**
 * Pop a DUZP_t from the end of a vector of DUZP_t.
 *
 * This function automatically shrinks the vector when
 * the allocated space is 3 times or more than the
 * size of the vector.
 **/
DUZP_t* DUZP::Factoring::vec_DUZP_t_pop(vec_DUZP_t* a) {

	if (a->size == 0) {
		fprintf(stderr,"BPAS, error: cannot pop an empty vector.\n");
		exit(1);
	}
	DUZP_t* ret = makeConstPolynomial_DUZP(1,0);
	DUZP_t* tmp;
	tmp = a->polys[a->size-1];
	a->polys[a->size-1] = ret;
	ret = tmp;
	a->size--;
	if (a->alloc >= 3*a->size) {
		long newalloc = a->alloc/2;
		for (int i=newalloc; i<a->alloc; ++i)
			freePolynomial_DUZP(a->polys[i]);
		a->polys = (DUZP_t**) realloc(a->polys,sizeof(DUZP_t*)*newalloc);
		a->alloc = a->alloc/2;
	}
	return ret;
}

/**
 * Free the allocated space of a vector of DUZP_t.
 **/
void DUZP::Factoring::vec_DUZP_t_free(vec_DUZP_t* a) {

	for (int i=0; i<a->alloc; ++i)
		freePolynomial_DUZP(a->polys[i]);
	free(a->polys);

}

/**
 * Initialize a vector of DUZP_t/long pairs
 * (allocation of 1, size 0).
 **/
void DUZP::Factoring::vec_DUZP_long_t_init(vec_DUZP_long_t* a) {

	a->pairs = (DUZP_long_t*) malloc(sizeof(DUZP_long_t));
	a->pairs[0].a = makeConstPolynomial_DUZP(1,0);
	a->pairs[0].b = 0;
	a->alloc = 1;
	a->size = 0;

}


/**
 * Initialize a vector of DUZP_t/long pairs
 * (allocation alloc; size 0).
 **/
void DUZP::Factoring::vec_DUZP_long_t_init2(vec_DUZP_long_t* a, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->pairs = (DUZP_long_t*) malloc(sizeof(DUZP_long_t)*alloc);
	for (int i=0; i<alloc; ++i) {
		a->pairs[i].a = makeConstPolynomial_DUZP(1,0);
		a->pairs[i].b = 0;
	}
	a->alloc = alloc;
	a->size = 0;

}


/**
 * Initialize a vector of DUZP_t/long pairs
 * (allocation alloc; size alloc).
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_DUZP_long_t_init_set(vec_DUZP_long_t* a, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->pairs = (DUZP_long_t*) malloc(sizeof(DUZP_long_t)*alloc);
	for (int i=0; i<alloc; ++i) {
		a->pairs[i].a = makeConstPolynomial_DUZP(1,0);
		a->pairs[i].b = 0;
	}
	a->alloc = alloc;
	a->size = alloc;

}

/**
 * Set the size of a vector of DUZP_t/long pairs. If
 * alloc < size, then the vector will be realloc'd to allow
 * the size to be set.
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_DUZP_long_t_set(vec_DUZP_long_t* a, long size) {

	if (size < 0) {
		fprintf(stderr,"BPAS, error: vector size must be non-negative.\n");
		exit(1);
	}

	if (size > a->alloc) {
		a->pairs = (DUZP_long_t*) realloc(a->pairs,sizeof(DUZP_long_t)*size);
		for (int i=a->alloc; i<size; ++i) {
			a->pairs[i].a = makeConstPolynomial_DUZP(1,0);
			a->pairs[i].b = 0;
		}
		a->alloc = size;
	}
	a->size = size;

}

/**
 * Clear a vector of DUZP_t/long pairs. The size will be set
 * to zero and the alloc to 1.
 **/
void DUZP::Factoring::vec_DUZP_long_t_clear(vec_DUZP_long_t* a) {

	if (a->alloc == 0) {
		fprintf(stderr,"BPAS, error: vector must be initialized to clear it.\n");
		exit(1);
	}

	for (int i=0; i<a->alloc; ++i)
		freePolynomial_DUZP(a->pairs[i].a);
	a->pairs = (DUZP_long_t*) realloc(a->pairs,sizeof(DUZP_long_t));
	a->pairs[0].a = makeConstPolynomial_DUZP(1,0);
	a->alloc = 1;
	a->size = 0;

}


/**
 * Push a DUZP_t/long pair onto the end of a vector of
 * DUZP_t/long pairs.
 *
 * This function automatically expands the vector (by doubling)
 * if the allocated space is full.
 **/
void DUZP::Factoring::vec_DUZP_long_t_push(vec_DUZP_long_t* a, const DUZP_long_t b) {

	if (a->size == a->alloc) {
		a->pairs = (DUZP_long_t*) realloc(a->pairs,sizeof(DUZP_long_t)*2*a->alloc);
		for (int i=a->alloc; i<2*a->alloc; ++i) {
			a->pairs[i].a = makeConstPolynomial_DUZP(1,0);
			a->pairs[i].b = 0;
		}
		a->alloc *= 2;
	}
	DUZP_t* c = deepCopyPolynomial_DUZP(b.a);
	freePolynomial_DUZP(a->pairs[a->size].a);
	a->pairs[a->size].a = c;
	a->pairs[a->size].b = b.b;
	a->size = a->size+1;

}

/**
 * Pop a DUZP_t/long pair from the end of a vector of
 * DUZP_t/long pairs.
 *
 * This function automatically shrinks the vector when
 * the allocated space is 3 times or more than the
 * size of the vector.
 **/
DUZP_long_t DUZP::Factoring::vec_DUZP_long_t_pop(vec_DUZP_long_t* a) {

	if (a->size == 0) {
		fprintf(stderr,"BPAS, error: cannot pop an empty vector.\n");
		exit(1);
	}
	DUZP_long_t ret;
	DUZP_long_t tmp;
	tmp.a = a->pairs[a->size-1].a;
	tmp.b = a->pairs[a->size-1].b;
	a->pairs[a->size-1].a = makeConstPolynomial_DUZP(1,0);
	a->pairs[a->size-1].b = 0;
	ret.a = tmp.a;
	ret.b = tmp.b;
	a->size--;
	if (a->alloc >= 3*a->size) {
		long newalloc = a->alloc/2;
		for (int i=newalloc; i<a->alloc; ++i)
			freePolynomial_DUZP(a->pairs[i].a);
		a->pairs = (DUZP_long_t*) realloc(a->pairs,sizeof(DUZP_long_t)*newalloc);
		a->alloc = a->alloc/2;
	}
	return ret;
}

/**
 * Free the allocated space of a vector of DUZP_t/long pairs.
 **/
void DUZP::Factoring::vec_DUZP_long_t_free(vec_DUZP_long_t* a) {

	for (int i=0; i<a->alloc; ++i)
		freePolynomial_DUZP(a->pairs[i].a);
	free(a->pairs);

}

/**
 * Initialize a vector of duspoly_t (allocation of 1, size 0).
 **/
void DUZP::Factoring::vec_duspoly_t_init(vec_duspoly_t* a, Prime_ptr* Pptr) {

	a->polys = (duspoly_t**) malloc(sizeof(duspoly_t*));
	a->polys[0] = makePolynomial_spX(1);
	a->Pptr = Pptr;
	a->alloc = 1;
	a->size = 0;

}


/**
 * Initialize a vector of duspoly_t (allocation alloc; size 0).
 **/
void DUZP::Factoring::vec_duspoly_t_init2(vec_duspoly_t* a, Prime_ptr* Pptr, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->polys = (duspoly_t**) malloc(sizeof(duspoly_t*)*alloc);
	for (int i=0; i<alloc; ++i)
		a->polys[i] = makePolynomial_spX(1);
	a->Pptr = Pptr;
	a->alloc = alloc;
	a->size = 0;

}


/**
 * Initialize a vector of duspoly_t (allocation alloc; size alloc).
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_duspoly_t_init_set(vec_duspoly_t* a, Prime_ptr* Pptr, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->polys = (duspoly_t**) malloc(sizeof(duspoly_t*)*alloc);
	for (int i=0; i<alloc; ++i)
		a->polys[i] = makePolynomial_spX(1);
	a->Pptr = Pptr;
	a->alloc = alloc;
	a->size = alloc;

}


/**
 * Set the size of a vector of duspoly_t; if alloc < size, then
 * the vector will be realloc'd to allow the size to be set.
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_duspoly_t_set(vec_duspoly_t* a, long size) {

	if (size < 0) {
		fprintf(stderr,"BPAS, error: vector size must be non-negative.\n");
		exit(1);
	}

	if (size > a->alloc) {
		a->polys = (duspoly_t**) realloc(a->polys,sizeof(duspoly_t*)*size);
		for (int i=a->alloc; i<size; ++i)
			a->polys[i] = makePolynomial_spX(1);
		a->alloc = size;
	}
	a->size = size;

}

/**
 * Push a DUZP_t onto the end of a vector of duspoly_t.
 *
 * This function automatically expands the vector (by doubling)
 * if the allocated space is full.
 **/
void DUZP::Factoring::vec_duspoly_t_push(vec_duspoly_t* a, const duspoly_t* b) {

	if (a->size == a->alloc) {
		a->polys = (duspoly_t**) realloc(a->polys,sizeof(duspoly_t*)*2*a->alloc);
		for (int i=a->alloc; i<2*a->alloc; ++i) {
			a->polys[i] = makePolynomial_spX(1);
		}
		a->alloc *= 2;
	}
	duspoly_t* c = deepCopyPolynomial_spX(b);
	freePolynomial_spX(&(a->polys[a->size]));
	a->polys[a->size] = c;
	a->size = a->size+1;

}

/**
 * Pop a duspoly_t from the end of a vector of duspoly_t.
 *
 * This function automatically shrinks the vector when
 * the allocated space is 3 times or more than the
 * size of the vector.
 **/
duspoly_t* DUZP::Factoring::vec_duspoly_t_pop(vec_duspoly_t* a) {

	if (a->size == 0) {
		fprintf(stderr,"BPAS, error: cannot pop an empty vector.\n");
		exit(1);
	}
	duspoly_t* ret = makePolynomial_spX(1);
	duspoly_t* tmp;
	tmp = a->polys[a->size-1];
	a->polys[a->size-1] = ret;
	ret = tmp;
	a->size--;
	if (a->alloc >= 3*a->size) {
		long newalloc = a->alloc/2;
		for (int i=newalloc; i<a->alloc; ++i)
			freePolynomial_spX(&(a->polys[i]));
		a->polys = (duspoly_t**) realloc(a->polys,sizeof(duspoly_t*)*newalloc);
		a->alloc = a->alloc/2;
	}
	return ret;
}

/**
 * Free the allocated space of a vector of duspoly_t.
 **/
void DUZP::Factoring::vec_duspoly_t_free(vec_duspoly_t* a) {

	for (int i=0; i<a->alloc; ++i)
		freePolynomial_spX(&(a->polys[i]));
	free(a->Pptr);
	free(a->polys);

}

/**
 * Initialize a vector of mpz_t (allocation 1; size 0).
 **/
void DUZP::Factoring::vec_mpz_t_init(vec_mpz_t* a) {

	a->array = (mpz_t*) malloc(sizeof(mpz_t));
	mpz_init(a->array[0]);
	a->alloc = 1;
	a->size = 0;

}

/**
 * Initialize a vector of mpz_t (allocation alloc; size 0).
 **/
void DUZP::Factoring::vec_mpz_t_init2(vec_mpz_t* a, long alloc) {

	if (alloc <= 0)
		alloc = 1;
	a->array = (mpz_t*) malloc(sizeof(mpz_t)*alloc);
	for (int i=0; i<alloc; ++i)
		mpz_init(a->array[i]);
	a->alloc = alloc;
	a->size = 0;

}

/**
 * Set the size of a vector of mpz_t; if alloc < size, then
 * the vector will be realloc'd to allow the size to be set.
 *
 * This function is intended for use only when the user is
 * manually maintaining the size of the vector for efficiency
 * or control purposes. For automatic size handling the push
 * and pop functions should be used instead.
 **/
void DUZP::Factoring::vec_mpz_t_set(vec_mpz_t* a, long size) {

	if (size < 0) {
		fprintf(stderr,"BPAS, error: vector size must be non-negative.\n");
		exit(1);
	}

	if (size > a->alloc) {
		a->array = (mpz_t*) realloc(a->array,sizeof(mpz_t)*size);
		for (int i=a->alloc; i<size; ++i)
			mpz_init(a->array[i]);
		a->alloc = size;
	}
	a->size = size;

}

/**
 * Push a mpz_t onto the end of a vector of mpz_t.
 *
 * This function automatically expands the vector (by doubling)
 * if the allocated space is full.
 **/
void DUZP::Factoring::vec_mpz_t_push(vec_mpz_t* a, const mpz_t b) {

	if (a->size == a->alloc) {
		a->array = (mpz_t*) realloc(a->array,sizeof(mpz_t)*2*a->alloc);
		for (int i=a->alloc; i<2*a->alloc; ++i) {
			mpz_init(a->array[i]);
		}
		a->alloc *= 2;
	}

	mpz_set(a->array[a->size],b);
	a->size = a->size+1;

}

/**
 * Push a long int onto a vector of mpz_t.
 *
 * This function automatically expands the vector (by doubling)
 * if the allocated space is full.
 **/
void DUZP::Factoring::vec_mpz_t_push_l(vec_mpz_t* a, long b) {

	mpz_t c;
	mpz_init(c);
	mpz_set_si(c,b);
	if (a->size == a->alloc) {
		a->array = (mpz_t*) realloc(a->array,sizeof(mpz_t)*2*a->alloc);
		for (int i=a->alloc; i<2*a->alloc; ++i)
			mpz_init(a->array[i]);
		a->alloc *= 2;
	}
	mpz_set(a->array[a->size],c);
	a->size = a->size+1;
	mpz_clear(c);

}

/**
 * Pop a mpz_t from the end of a vector of DUZP_t.
 *
 * This function automatically shrinks the vector when
 * the allocated space is 3 times or more than the
 * size of the vector.
 **/
void DUZP::Factoring::vec_mpz_t_pop(mpz_t* r, vec_mpz_t* a) {

	if (a->size == 0) {
		fprintf(stderr,"BPAS, error: cannot pop an empty vector.\n");
		exit(1);
	}
	mpz_clear(*r);
	mpz_init(*r);
	mpz_set(*r,a->array[a->size-1]);
	mpz_set_ui(a->array[a->size-1],0);
	a->size--;
	if (a->alloc >= 3*a->size) {
		long newalloc = a->alloc/2;
		for (int i=newalloc; i<a->alloc; ++i)
			mpz_clear(a->array[i]);
		a->array = (mpz_t*) realloc(a->array,sizeof(mpz_t)*newalloc);
		a->alloc = a->alloc/2;
	}
}

/**
 * Free the allocated space of a vector of mpz_t.
 **/
void DUZP::Factoring::vec_mpz_t_free(vec_mpz_t* a) {

	for (int i=0; i<a->alloc; ++i)
		mpz_clear(a->array[i]);
	free(a->array);

}

/**
 * Initialize a vector of vectors of mpz_t.
 *
 * The outer vector is set to have the supplied length, with
 * each inner vector in it being initialized to size 1.
 **/
void DUZP::Factoring::vec_vec_mpz_t_init(vec_vec_mpz_t* a, long length) {

	a->array = (mpz_t**) malloc(sizeof(mpz_t*)*length);
	a->length = length;
	for (int i=0; i<length; ++i) {
		a->array[i] = (mpz_t*) malloc(sizeof(mpz_t));
		mpz_init(a->array[i][0]);
	}
	a->depth = 1;
}

/**
 * Set the size of a vector of vectors of mpz_t based on a
 * supplied length of the outer vector and depth of the inner
 * vector.
 *
 * This routine will automatically extend the current length
 * and depth of the vector of vectors to acommodate the required
 * new length and depth by reallocating the underlying arrays.
 * Accordingly there is no information loss when the length and
 * depth are non-decreasing, with information only being lost
 * when either the length or depth decreases. Note that each
 * inner vector in the outer vector has the same depth.
 **/
void DUZP::Factoring::vec_vec_mpz_t_set(vec_vec_mpz_t* a, long length, long depth) {

	long l_diff = length - a->length;
	long d_diff = depth - a->depth;

	if (l_diff < 0) { // new length is smaller
//		fprintf(stderr,"decreasing length...\n");
		for (int i=length; i<a->length; ++i) {
			for (int j=0; j<a->depth; ++j) {
				mpz_clear(a->array[i][j]);
			}
			free(a->array[i]);
		}
		a->array = (mpz_t**) realloc(a->array,sizeof(mpz_t*)*length);
	}
	else if (l_diff > 0) { // new length is larger
//		fprintf(stderr,"increasing length...\n");
		a->array = (mpz_t**) realloc(a->array,sizeof(mpz_t*)*length);
		for (int i=a->length; i<length; ++i) {
			a->array[i] = (mpz_t*) malloc(sizeof(mpz_t)*depth);
			for (int j=0; j<depth; ++j) {
				mpz_init(a->array[i][j]);
			}
		}
	}
	if (d_diff < 0) { // new depth is smaller
//		fprintf(stderr,"decreasing depth...\n");
		long bnd = MIN(length,a->length);
		for (int i=0; i<bnd; ++i) {
			for (int j=depth; j<a->depth; ++j) {
				mpz_clear(a->array[i][j]);
			}
			a->array[i] = (mpz_t*) realloc(a->array[i],sizeof(mpz_t)*depth);
		}
	}
	else if (d_diff > 0) { // new depth is larger
//		fprintf(stderr,"increasing depth...\n");
		long bnd = MIN(length,a->length);
		for (int i=0; i<bnd; ++i) {
			a->array[i] = (mpz_t*) realloc(a->array[i],sizeof(mpz_t)*depth);
			for (int j=a->depth; j<depth; ++j) {
				mpz_init(a->array[i][j]);
			}
		}
	}
	a->length = length;
	a->depth = depth;
}

/**
 * Free the allocated space of a vector of vectors of mpz_t.
 **/
void DUZP::Factoring::vec_vec_mpz_t_free(vec_vec_mpz_t* a) {

	for (int i=0; i<a->length; ++i) {
		for (int j=0; j<a->depth; ++j)
			mpz_clear(a->array[i][j]);
		free(a->array[i]);
	}
	a->depth = 0;
	free(a->array);
	a->length = 0;
}

/**
 * Possible degree computation: rational version.
 *
 * Computes all the possible degrees of products of
 * p-adic factors in polys in the following way:
 *
 *   - degs[i] encodes as a bit vector the possible;
 *   - degrees of products of size m in the set
 *     {polys[i],...,polys[len-1]}.
 */
void DUZP::Factoring::ComputePossibleDegreesRational(mpz_t* degs, const vec_DUZP_t* factors, long m) {


	long len = factors->size;
	DUZP_t** polys = factors->polys;
	mpz_t tmp,old,one;
	mpz_init(one);
	mpz_init(old);
	mpz_init(tmp);
	mpz_set_ui(one,1);

	if (m == 0)
		goto cleanup;

	if (m < 1 || m > len) {
		fprintf(stderr,"BPAS, error: invalid input to computePossibleDegrees\n");
		fprintf(stderr,"m = %ld, len = %ld\n",m,len);
		exit(1);
	}

	mpz_set(degs[len-1],one);
	mpz_mul_2exp(degs[len-1],degs[len-1],polys[len-1]->lt);

	for (int i = len-2; i >= 0; --i) {
		mpz_set(tmp,one);
		mpz_mul_2exp(tmp,tmp,polys[i]->lt);
		mpz_ior(degs[i],tmp,degs[i+1]);
	}

	for (int i = 2; i <= m; ++i) {
		mpz_set(old,degs[len-i]);
		mpz_mul_2exp(degs[len-i],degs[len-i+1],polys[len-i]->lt);

		for (int j = len-i-1; j >= 0; --j) {
			mpz_mul_2exp(tmp,old,polys[j]->lt);
			mpz_set(old,degs[j]);
			mpz_ior(degs[j],degs[j+1],tmp);
		}
	}

	cleanup:

	mpz_clear(one);
	mpz_clear(old);
	mpz_clear(tmp);
}

/**
 * Possible degree computation: rational version.
 *
 * Computes all the possible degrees of products of
 * modular factors in polys in the following way:
 *
 *   - degs[i] encodes as a bit vector the possible;
 *   - degrees of products of size m in the set
 *     {polys[i],...,polys[len-1]}.
 */
void DUZP::Factoring::ComputePossibleDegreesModular(mpz_t* degs, const factors_t* factors, long m, const Prime_ptr* pp) {


	long len = factors->alloc;
	duspoly_t** polys = factors->polys;
	mpz_t tmp,old,one;
	mpz_init(one);
	mpz_init(old);
	mpz_init(tmp);
	mpz_set_ui(one,1);

	if (m == 0)
		return;

	if (m < 1 || m > len)
		fprintf(stderr,"BPAS, error: invalid input to computePossibleDegrees\n");

	mpz_set(degs[len-1],one);
	mpz_mul_2exp(degs[len-1],degs[len-1],polys[len-1]->lt);

	for (int i = len-2; i >= 0; --i) {
		mpz_set(tmp,one);
		mpz_mul_2exp(tmp,tmp,polys[i]->lt);
		mpz_ior(degs[i],tmp,degs[i+1]);
	}

	for (int i = 2; i <= m; ++i) {
		mpz_set(old,degs[len-i]);
		mpz_mul_2exp(degs[len-i],degs[len-i+1],polys[len-i]->lt);

		for (int j = len-i-1; j >= 0; --j) {
			mpz_mul_2exp(tmp,old,polys[j]->lt);
			mpz_set(old,degs[j]);
			mpz_ior(degs[j],degs[j+1],tmp);
		}
	}
}

/**
 * Remove m p-adic factors from a vector according
 * to an index array.
 **/
void RemoveFactorsRational(vec_DUZP_t* padicFactors, const long* idxs, long m) {

	DUZP_t** polys = padicFactors->polys;
	DUZP_t* tmp;
	long n = padicFactors->alloc;

	int i=0;
	for (int j=0; j < n; j++) {
		if (i < m && j == idxs[i])
			i++;
		else { // swap polys[j] and polys[j-i]
			tmp = polys[j];
			polys[j] = polys[j-i];
			polys[j-i] = tmp;
		}
	}

	DUZP_t* r;
	for (i=n-m; i<n; ++i) {
		r = vec_DUZP_t_pop(padicFactors);
		freePolynomial_DUZP(r);
	}
//	DUZP_t** re_alloc = NULL;
//	for (i=n-m; i<n; ++i) {
//		freePolynomial_DUZP(padicFactors->polys[i]);
//	}
//	re_alloc = (DUZP_t**) realloc(padicFactors->polys,sizeof(DUZP_t*)*(n-m));
//	if (re_alloc != NULL) {
//		padicFactors->polys = re_alloc;
//	}
//	padicFactors->alloc = n-m;
//	padicFactors->size = n-m;

}

/**
 * Remove m modular factors from a vector according
 * to an index array.
 **/
void DUZP::Factoring::RemoveFactorsModular(factors_t* modularFactors, const long* idxs, long m) {

	duspoly_t** polys = modularFactors->polys;
	duspoly_t* tmp;
	long n = modularFactors->alloc;

	int i=0;
	for (int j=0; j < n; j++) {
		if (i < m && j == idxs[i])
			i++;
		else { // swap polys[j] and polys[j-i]
			tmp = polys[j];
			polys[j] = polys[j-i];
			polys[j-i] = tmp;
		}
	}

	duspoly_t** re_alloc = NULL;
	for (i=n-m; i<n; ++i) {
		freePolynomial_spX(&(modularFactors->polys[i]));
	}
	re_alloc = (duspoly_t**) realloc(modularFactors->polys,sizeof(duspoly_t*)*(n-m));
	if (re_alloc != NULL) {
		modularFactors->polys = re_alloc;
	}
	modularFactors->alloc = n-m;
}

/**
 * Convert a bit vector v to an array x of integers
 * of size n.
 **/
void DUZP::Factoring::UnpackBitVector(int* x, const mpz_t v, long n) {

	for (int i=0; i<n; ++i)
		x[i] = mpz_tstbit(v,i);
}

/**
 * Test whether the constant term of the current product
 * of m p-adic factors (determined by idxs) mod P divides
 * the constant term of the polynomial being factored.
 * The leading coefficient of the polynomial being factored
 * is used to construct the corresponding constant term from
 * the monic p-adic factors.
 **/
int TestConstantTermRational(mpz_t* prod, long ProdLen, const vec_DUZP_t* padicFactors, const mpz_t P, const long* idxs, long m, const mpz_t constTerm, const mpz_t leadingCoeff) {
	// Note that P = p^e is the power of the prime to which the modular factors are lifted //
	mpz_t t1, t2;
	mpz_init(t1);
	mpz_init(t2);

	if (ProdLen == 0) {
		mpz_mul(prod[0],leadingCoeff,padicFactors->polys[idxs[0]]->coefs[0]);
		ProdLen++;
	}

	for (int i=ProdLen; i < m; i++)
		mpz_mul(prod[i],prod[i-1],padicFactors->polys[idxs[i]]->coefs[0]);

	ProdLen = m-1;

	mpz_mod(t1,prod[m-1],P);
	mpz_set(t2,P);
	mpz_fdiv_q_2exp(t2,t2,1);
	if (mpz_cmp(t1,t2) > 0) {
		mpz_set(t2,P);
		mpz_sub(t1,t1,t2);
	}

	int ret = mpz_divisible_p(constTerm,t1);
	mpz_clear(t1);
	mpz_clear(t2);

	return ret;
}

/**
 * Computes a potential factor g of the polynomial being
 * factored from m of the p-adic factors (polys) mod P
 *  according to the index array idxs.
 **/
void DUZP::Factoring::ComputePotentialFactorRational(DUZP_t** g, DUZP_t** polys, const long* idxs, long m, const mpz_t P) {

	if (m != 0) {
		DUZP_t* res = deepCopyPolynomial_DUZP(polys[idxs[0]]);

		for (int i=1; i<m; ++i) {
			multiplyPolynomials_DUZP_inp(polys[idxs[i]], &res);
			applyModuloSymmetric_DUZP_inp(res,P);
		}

		if (*g != NULL)
			freePolynomial_DUZP(*g);
		*g = res;
	}
}

/**
 * Computes the complement g of a potential factor h of
 * the polynomial f being factored, i.e., f = gh, from
 * m of the n p-adic factors (polys) mod P according to
 * the index array idxs.
 **/
void ComputeComplementaryPotentialFactorRational(DUZP_t** g, DUZP_t** polys, long n, const long* idxs, long m, const mpz_t P) {

	if (m != 0) {
		DUZP_t* res = makeConstPolynomial_DUZP(1,1);

		int i=0;
		for (int j=0; j<n; ++j) {
			if (i < m && j == idxs[i]) {
				++i;
			}
			else {
				multiplyPolynomials_DUZP_inp(polys[j], &res);
				applyModuloSymmetric_DUZP_inp(res,P);
			}
		}

		if (*g != NULL)
			freePolynomial_DUZP(*g);
		*g = res;
	}
}

/**
 * Returns the maximum number of bits needed to store the
 * coefficients of a polynomial f.
 **/
long DUZP::Factoring::MaxBits(const DUZP_t* f)
{
	long max,tmp;
	max = 0;

	for (int i=0; i<=f->lt; i++) {
		tmp = mpz_sizeinbase(f->coefs[i],2);
		if (tmp > max)
			max = tmp;
	}

	return max;
}

/**
 * Computes the log2 of the Landau-Mignotte bound
 *      ∥h∥∞ ≤ (n + 1)^0.5 * 2^k * ∥f∥∞
 * where n = deg(f) and k = deg(h), h a rational factor of f.
 **/
long DUZP::Factoring::LandauMignotteBoundInBits(const DUZP_t* f, long k) {

	long log_norm = MaxBits(f);
	long n = f->lt;
	long ret = (long) ceil(log2(n+1));
	ret = 1 + ((ret - 1) / 2); // ceiling division ret/2 avoiding overflow
	ret += k + log_norm;
	return ret;
}

/**
 * Exponential complexity search for the rational factors of
 * the polynomial f_in from the p-adic factors of f_in mod P,
 * P = p^e. The algorithm will search for rational factors
 * formed from m p-adic factors, based on whether the product
 * of p-adic factors has coefficients less than a supplied
 * bound.
 **/
void DUZP::Factoring::NaiveFactorSearch(DUZP_t** f_in, vec_DUZP_t* rationalFactors, vec_DUZP_t* padicFactors, mpz_t P, long_t p, int e, long m, long_t bound) {

//	fprintf(stderr,"    [NFS] Entering NaiveFactorSearch:\n");

	DUZP_t* f = deepCopyPolynomial_DUZP(*f_in);

	DUZP_t** polys = padicFactors->polys; // the modular factors left to consider
	long n = padicFactors->size; // number of modular factors left to consider

	long* idxs = (long*) malloc(sizeof(long)*m); // list of indices in polys for potential factor
	long* degs = (long*) malloc(sizeof(long)*m); // degrees of partial products of polys in idxs up to a given index
	mpz_t* prod = (mpz_t*) malloc(sizeof(mpz_t)*m); // partial products of constant terms of polys in idxs up to a given index
	for (int i=0; i<m; ++i) {
		mpz_init(prod[i]);
	}

	long pLSize = n;
	long prodLength = 0; // number of terms used for the current valid product of constant terms

	// possible degrees of factors of size m
	mpz_t* possDeg = (mpz_t*) malloc(sizeof(mpz_t)*pLSize);
	for (int i=0; i<pLSize; ++i)
		mpz_init(possDeg[i]);
	ComputePossibleDegreesRational(possDeg,padicFactors,m);

	// unpacked bit vector of possible degrees
	int* upkdPossDeg = (int*) malloc(sizeof(int)*f->lt);

	// constant term of f
	mpz_t constTerm;
	mpz_init(constTerm);
	mpz_mul(constTerm,f->coefs[0],f->coefs[f->lt]);

	mpz_t leadingCoeff; // leading coefficient of f mod P
	mpz_init(leadingCoeff);

	// reduce mpz_t lc mod P
	mpz_mod(leadingCoeff,f->coefs[f->lt],P);
//	gmp_fprintf(stderr,"    [NFS] leadingCoeff = %Zd\n",leadingCoeff);

	DUZP_t* g = NULL; // potential factor
	DUZP_t* h = NULL; // quotient of f by potential factor

	// bookkeeping variables
	long i = 0;
//	long count = 0;
	long state = 0;

	idxs[0] = 0;

//	fprintf(stderr,"    [NFS] Starting while loop:\n");

	while (idxs[0] <= n-m) {

		UnpackBitVector(upkdPossDeg,possDeg[idxs[0]],f->lt); // f->lt could perhaps be smaller here

		degs[0] = polys[idxs[0]]->lt;

		i = 1;
		state = 0;
		prodLength = 0;

		for (;;) {

//			fprintf(stderr,"    [NFS] I = {");
//			for (int k=0; k<i; ++k) {
//				fprintf(stderr,"%ld",idxs[k]);
//				if (k != i-1)
//					fprintf(stderr,", ");
//			}
//			fprintf(stderr,"}\n");

			if (i < prodLength)
				prodLength = i;

			if (i == m) {

//				fprintf(stderr,"    [NFS] i=m condition met:\n");
//				fprintf(stderr,"f:\n");
//				printPoly_DUZP(f,"x");

				state = 1;  // default continuation state

//				fprintf(stderr,"    [NFS] checking possible degrees:\n");
				if (!upkdPossDeg[degs[m-1]]) { // product is not possible
					i--;
//					count++;
					continue;
				}

//				fprintf(stderr,"    [NFS] performing constant term test:\n");
				if (!TestConstantTermRational(prod, prodLength, padicFactors, P, idxs, m, constTerm, leadingCoeff)) {
					i--;
//					count += 100; // accounting for when to change primes, I think
					continue;
				}

//				count += 1000; // accounting for when to change primes, I think
				if (2*degs[m-1] <= f->lt) {
//					fprintf(stderr,"    [NFS] starting regular potential factor computation:\n");

					ComputePotentialFactorRational(&g,polys,idxs,m,P);
					multiplyByInteger_DUZP_inp(g, leadingCoeff);
					applyModuloSymmetric_DUZP_inp(g,P);

//					fprintf(stderr,"    [NFS] checking max bits:\n");
					if (MaxBits(g) > bound) {
						i--;
						continue;
					}

//					fprintf(stderr,"    [NFS] computing primitive part:\n");

//					g = primitivePart_DUZP(g);

					DUZP_t* gg;
					gg = primitivePart_DUZP(g);
					freePolynomial_DUZP(g);
					g = gg;

//					fprintf(stderr,"    [NFS] performing division test:\n");
					if (!divideTest_DUZP(f, g, &h)) {
						i--;
						continue;
					}

					// we found a factor!
					vec_DUZP_t_push(rationalFactors,g);
					#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
						fprintf(stderr,"    [NFS] degree %ld factor found\n",g->lt);
					#endif

					freePolynomial_DUZP(f);
					f = h;
					h = NULL;

					mpz_mod(leadingCoeff,f->coefs[f->lt],P);
					mpz_mul(constTerm,f->coefs[0],f->coefs[f->lt]);

				}
				else {
//					fprintf(stderr,"    [NFS] starting complementary potential factor computation:\n");

					ComputeComplementaryPotentialFactorRational(&g,polys,n,idxs,m,P);
					multiplyByInteger_DUZP_inp(g, leadingCoeff);
					applyModuloSymmetric_DUZP_inp(g,P);

//					fprintf(stderr,"    [NFS] checking max bits:\n");
					if (MaxBits(g) > bound) {
						i--;
						continue;
					}

//					g = primitivePart_DUZP(g);

					DUZP_t* gg;
					gg = primitivePart_DUZP(g);
					freePolynomial_DUZP(g);
					g = gg;

//					fprintf(stderr,"    [NFS] performing division test:\n");
					if (!divideTest_DUZP(f, g, &h)) {
						i--;
						continue;
					}

					// we found a factor!
					vec_DUZP_t_push(rationalFactors,h);
					#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
						fprintf(stderr,"    [NFS] degree %ld factor found\n",h->lt);
					#endif
//					fprintf(stderr,"    [NFS] h =\n");
//					printPoly_DUZP(h,x);

					freePolynomial_DUZP(f);
					f = g;
					g = NULL;

					mpz_mod(leadingCoeff,f->coefs[f->lt],P);
					mpz_mul(constTerm,f->coefs[0],f->coefs[f->lt]);

				}

//				fprintf(stderr,"    [NFS] before removing factors\n");
				RemoveFactorsRational(padicFactors,idxs,m);
//				fprintf(stderr,"    [NFS] after removing factors\n");
				n = padicFactors->size;
//				count = 0;

				if (2*m > n && m != n)
					goto done;
				else
					break;
			}
			else if (state == 0) {
				idxs[i] = idxs[i-1] + 1;
				degs[i] = degs[i-1] + polys[idxs[i]]->lt;
				i++;
			}
			else { // state == 1
				idxs[i]++;
				if (i == 0) break;

				if (idxs[i] > n-m+i)
					i--;
				else {
					degs[i] = degs[i-1] + polys[idxs[i]]->lt;
					i++;
					state = 0;
				}
			}
		}

	}

	done:

//	fprintf(stderr,"    [NFS] leaving NaiveFactorSearch.\n");
	freePolynomial_DUZP(*f_in);
	*f_in = f;
//	*f_in = deepCopyPolynomial_DUZP(f);
	freePolynomial_DUZP(g);
	freePolynomial_DUZP(h);
	free(idxs);
	free(degs);
	for (int i=0; i<m; ++i) {
		mpz_clear(prod[i]);
	}
	free(prod);
	for (int i=0; i<pLSize; ++i)
		mpz_clear(possDeg[i]);
	free(possDeg);
	free(upkdPossDeg);
	mpz_clear(constTerm);
	mpz_clear(leadingCoeff);
}

/**
 * Compute a bound bnd on the size of the roots of the
 * polynomial f.
 **/
void DUZP::Factoring::computeRootBound(mpz_t* bnd, const DUZP_t* f) {

	int NTL_ROOT_BOUND = 1; // whether to use NTL's root bound algorithm (better than Cauchy bound)
	int MULTIPLY_BY_LC = 1; // whether to multiply the root by the leading coefficient of f (needed for factorization algorithm

//	mpz_init(*bnd);

	DUZP_t* g;
	int N;

	N = f->lt;
	g = deepCopyPolynomial_DUZP(f);
	mpz_t* coefs = g->coefs;

	if (NTL_ROOT_BOUND) {

		if (mpz_sgn(f->coefs[0]) == 0) {
			gmp_fprintf(stderr,"BPAS, error: constant term must be non-zero to compute root bound.\n");
			exit(1);
		}

		mpz_t lower,upper,mean,check;

		if (mpz_sgn(coefs[N]) == -1) mpz_neg(coefs[N], coefs[N]);

		for (int i=0; i<N; ++i) {
			if (mpz_sgn(coefs[i]) == 1) mpz_neg(coefs[i], coefs[i]);
		}


		mpz_init(lower);
		mpz_init(upper);
		mpz_init(mean);

		mpz_set_si(lower,0);
		mpz_set_si(upper,1);

//		AltArrZ_t* G,H;

		mpz_t val;
		mpz_init(val);

		evaluate_DUZP(g,upper,val);

		while (mpz_sgn(val) == -1) {
			mpz_mul_ui(upper,upper,2);
			evaluate_DUZP(g,upper,val);
		}

		mpz_init(check);
		mpz_sub(check,upper,lower);
		mpz_sub_ui(check,check,1);

		// lower < root <= upper

		while (mpz_sgn(check) == 1) { // check if (upper - lower) > 1
			mpz_add(mean,upper,lower);
			mpz_fdiv_q_ui(mean,mean,2);

			evaluate_DUZP(g,mean,val);
			if (mpz_sgn(val) == -1)
				mpz_set(lower,mean);
			else
				mpz_set(upper,mean);

			mpz_sub(check,upper,lower);
			mpz_sub_ui(check,check,1);
		}


		if (MULTIPLY_BY_LC)
			mpz_mul(*bnd,upper,coefs[N]);
		else
			mpz_set(*bnd,upper);
		mpz_clear(lower);
		mpz_clear(upper);
		mpz_clear(mean);
		mpz_clear(check);
		mpz_clear(val);

	}
	else { // Compute Cauchy bound

		mpf_t val2,den,max;
		mp_bitcnt_t cnt = 64;
		mpf_set_default_prec (cnt);
		mpf_init(val2);
		mpf_init(den);
		mpf_init(max);

		mpf_set_z(den,coefs[N]);

		for (int i=0; i<N; ++i) {

			mpf_set_z(val2,coefs[i]);
			mpf_div(val2,val2,den);
			mpf_abs(val2,val2);

			if (mpf_cmp(max,val2) < 0)
				mpf_set(max,val2);

		}


		mpf_add_ui(max,max,1);
		mpf_ceil(max,max);

		mpz_set_f(*bnd,max);
		if (MULTIPLY_BY_LC)
			mpz_mul(*bnd,*bnd,coefs[N]);

		mpf_clear(val2);
		mpf_clear(den);
		mpf_clear(max);

	}

	freePolynomial_DUZP(g);
}

/**
 * Computes a cut-off value determining when the
 * factorization algorithm should switch from
 * sparse to dense mode (the latter using a
 * random linear transformation A of traces of the
 * polynomial being factored. After this cut-off,
 * the value return determines the number of
 * random mixtures of the traces computed
 * (i.e., the dimension of the image of A).
 **/
long DUZP::Factoring::d1_value(long deltaInBits, long r, long s) {

	double tmp;
	tmp = (double) r;
	tmp *= (double) s;
	tmp /= (double) deltaInBits;
	return (long) ( 0.30*tmp ) + 1;
}

/**
 * Compute the excess precision needed for the current
 * iteration of the van Hoeij algorithm given as the
 * prime p to the power of delta, where delta is given
 * as a number of bits. Then delta and p^delta are
 * returned.
 **/
void DUZP::Factoring::Computepdelta(long* delta, mpz_t* pdelta, long_t p, long deltaInBits) {

	mpz_t res;
	mpz_init(res);
	long i;

	i = *delta;
	mpz_set(res,*pdelta);

	while (mpz_sizeinbase(res,2) <= deltaInBits) {
		i++;
		mpz_mul_sll(res,res,p);
	}

	*delta = i;
	mpz_set(*pdelta,res);
	mpz_clear(res);
}

/**
 * Compute the vector of values b for p^b such that p^{b_i}
 * is a bound for the ith trace of a rational factor of the
 * polynomial f being factored. This computation requires a
 * bound on the roots of f, with n = deg(f).
 **/
void Computepb(NTL::vec_long& b, vec_mpz_t* pb, long_t p, long d, const mpz_t rootBound, long n) {


	mpz_t tmp1,tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);
	long i;

	mpz_set(tmp1,rootBound);
	mpz_pow_ui(tmp1,tmp1,d);
	mpz_mul_ui(tmp1,tmp1,2*n); // tmp1 = 2*n*rootBound^d;

	if (d == 1) {
		i = 0;
		mpz_set_ui(tmp2,1);
	}
	else {
		i = b(d-1);
      	mpz_set(tmp2,pb->array[d-2]);
	}

	while (mpz_cmp(tmp2,tmp1) <= 0) {
		i++;
		mpz_mul_sll(tmp2,tmp2,p);
	}

	b.SetLength(d);
	b(d) = i;

	vec_mpz_t_push(pb,tmp2);

	mpz_clear(tmp1);
	mpz_clear(tmp2);

}

/**
 * Compute an effective value b_eff for the bound p^{b_eff} on
 * traces of a rational factor of the polynomial being factored
 * when a random linear transformation of the traces is being
 * used. This computation requires a bound on the roots of f,
 * with n = deg(f).
 **/
void DUZP::Factoring::Computepb_eff(long* b_eff, mpz_t* pb_eff, long_t p, long d, const mpz_t rootBound, long n, long ran_bits) {

	mpz_clear(*pb_eff);
	mpz_init(*pb_eff);
	mpz_t tmp1,tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);
	long i;

	if (mpz_cmp_ui(rootBound,1) == 0) {
		mpz_set_ui(tmp1,n);
		mpz_mul_ui(tmp1,tmp1,d);
		mpz_mul_2exp(tmp1,tmp1,ran_bits+1);
	}
	else {
		mpz_set(tmp1,rootBound);
		mpz_pow_ui(tmp1,tmp1,d);
		mpz_mul_ui(tmp1,tmp1,n); // tmp1 = n*rootBound^d;
		mpz_mul_2exp(tmp1,tmp1,ran_bits+2);
	}

	i = 0;
	mpz_set_ui(tmp2,1);

	while (mpz_cmp(tmp2,tmp1) <= 0) {
		i++;
		mpz_mul_sll(tmp2,tmp2,p);
	}

	*b_eff = i;
	mpz_set(*pb_eff,tmp2);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
}

/**
 * Computes the dth trace of of the kth p-adic factor f mod P,
 * with Tr storing all of the computed traces of all of the
 * p-adic factors mod P.
 *
 * The routine ensures that the input f is monic, and deg(f)>0.
 * The prime power P must be > 1.
 * Tr->length must be >= d and Tr->array[i-1], for i = 1..d-1,
 * should be the Tr_i(f) mod P (in van Hoeij's notation).
 * The quantity Tr_d(f) mod P is computed, and stored in
 * Tr->array[d-1].
 */
void DUZP::Factoring::ComputeTrace(vec_vec_mpz_t* Tr, long k, const DUZP_t* f, long d, const mpz_t P) {

	long n = f->lt;

	// check arguments

	if (n <= 0 || mpz_cmp_ui(f->coefs[n],1) != 0) {
		fprintf(stderr,"BPAS, error: in ComputeTrace, degree of f is zero or f is not monic.\n");
		exit(1);
	}
	if (d <= 0) {
		fprintf(stderr,"BPAS, error: in ComputeTrace, d is not one or greater.\n");
		exit(1);
	}
	if (Tr->depth < d) {
		fprintf(stderr,"BPAS, error: in ComputeTrace, vec_vec_mpz_t has depth less than d.\n");
		exit(1);
	}
	if (mpz_cmp_ui(P,1) <= 0) {
		fprintf(stderr,"BPAS, error: in ComputeTrace, P too small.\n");
		exit(1);
	}

	mpz_t t1;
	mpz_t t2;
	mpz_init(t1);
	mpz_init(t2);

	// treat d > deg(f) separately

	if (d > n) {

		mpz_set_ui(t1,0);

		for (int i=1; i<=n; i++) {
			mpz_mul(t2,Tr->array[k][i+d-n-2],f->coefs[i-1]);
			mpz_add(t1,t1,t2);
		}
	}
	else {

		mpz_mul_ui(t1,f->coefs[n-d], d);

		for (int i=1; i<d; i++) {
			mpz_mul(t2,Tr->array[k][i-1],f->coefs[n-d+i]);
			mpz_add(t1,t1,t2);
		}
	}

	mpz_mod(t1,t1,P);
	mpz_sub(t1,P,t1);
	mpz_set(Tr->array[k][d-1],t1);

	mpz_clear(t1);
	mpz_clear(t2);
}

/**
 * Computes the two-sided cut (see van Hoeij's paper) of d
 * traces of the kth p-adic factor mod P of the polynomial
 * being factored. This computation requires the vector of
 * bounds p^{b_i} on the traces, the excess precision p^delta
 * above the largest bound p^{b_i} to which the p-adic factors
 * have been computed, as well as the leading coefficient lc
 * of the polynomial being factored.
 *
 * The routine requires that the number of traces d is > 0.
 * The depth of the vectors cutTraces and traces must be >= d.
 * The length of the vector pb must have length >= d
 */
void DUZP::Factoring::CutTraces(vec_vec_mpz_t* cutTraces, const vec_vec_mpz_t* traces, long k, long d, const vec_mpz_t pb, const mpz_t pdelta, const mpz_t P, const mpz_t lc) {

	if (d <= 0) {
		fprintf(stderr,"BPAS, error: in CutTraces, number of traces must be positive.\n");
		exit(1);
	}
	if (cutTraces->depth < d) {
		fprintf(stderr,"BPAS, error: in CutTraces, cutTraces vector must have depth >= d.\n");
		exit(1);
	}
	if (traces->depth < d) {
		fprintf(stderr,"BPAS, error: in CutTraces, traces vector must have depth >= d.\n");
		exit(1);
	}
	if (pb.size < d) {
		fprintf(stderr,"BPAS, error: in CutTraces, pb vector must have length >= d.\n");
		exit(1);
	}
	if (mpz_cmp_ui(P,1) <= 0) {
		fprintf(stderr,"BPAS, error: in CutTraces, prime power is too small.\n");
		exit(1);
	}

	mpz_t lcpow;
	mpz_t lcred;
	mpz_init_set_ui(lcpow,1);
	mpz_init(lcred);
	mpz_mod(lcred,lc,P);

	mpz_t pdelta_2;
	mpz_init(pdelta_2);
	mpz_fdiv_q_2exp(pdelta_2,pdelta,1);

	mpz_t t1;
	mpz_t t2;
	mpz_init(t1);
	mpz_init(t2);

	for (int i=0; i<d; i++) {
		mpz_mul(lcpow,lcpow,lcred);
		mpz_mod(lcpow,lcpow,P);
		mpz_mod(t1,traces->array[k][i],P);
		mpz_mul(t1,lcpow,t1);
		mpz_mod(t1,t1,P);

		mpz_fdiv_q_2exp(t2,pb.array[i],1);
		mpz_add(t1,t1,t2);
		mpz_fdiv_q(t1,t1,pb.array[i]); // TODO: ensure that floor division is what we want here
		mpz_mod(t1,t1,pdelta);
		if (mpz_cmp(t1,pdelta_2) > 0)
			mpz_sub(t1,t1,pdelta);

		mpz_set(cutTraces->array[k][i],t1);
	}

	mpz_clear(lcpow);
	mpz_clear(lcred);
	mpz_clear(pdelta_2);
	mpz_clear(t1);
	mpz_clear(t2);
}

/**
 * Computes the d1 two-sided cuts (see van Hoeij's paper) of
 * the d traces computed for the kth p-adic factor mod P of
 * the polynomial being factored in the case where the d
 * traces are transformed by a random linear transformation A.
 * This computation requires the vector of
 * bounds p^{b_i} on the traces, the excess precision p^delta
 * (above the largest bound p^{b_i}) to which the p-adic
 * factors have been computed, as well as the leading
 * coefficient lc of the polynomial being factored.
 *
 * The routine requires that the number of traces d is > 0.
 * The depth of the vectors cutTraces and traces must be >= d.
 * The length of the vector pb must have length >= d
 */
void DenseCutTraces(vec_vec_mpz_t* cutTraces, const vec_vec_mpz_t* traces, long k, long d, long d1, const mpz_t pb_eff, const mpz_t pdelta, const mpz_t P, const mpz_t lc, const vec_vec_mpz_t* A) {

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"    [DCT] DenseCutTraces:\n");
	#endif

	mpz_t pdelta_2;
	mpz_init(pdelta_2);
	mpz_div_2exp(pdelta_2,pdelta,1);

	mpz_t pb_eff_2;
	mpz_init(pb_eff_2);
	mpz_div_2exp(pb_eff_2,pb_eff,1);

	mpz_t acc;
	mpz_t t1;
	mpz_t t2;
	mpz_init(acc);
	mpz_init(t1);
	mpz_init(t2);

	mpz_t lcpow;
	mpz_t lcred;
	mpz_init(lcpow);
	mpz_init(lcred);
	mpz_mod(lcred,lc,P);

	for (int i=0; i<d1; ++i) {
		mpz_set_ui(lcpow,1);
		mpz_set_ui(acc,0);

		for (int j=0; j<d; ++j) {
			mpz_mul(lcpow,lcpow,lcred);
			mpz_mod(lcpow,lcpow,P);
			mpz_mul(t1,lcpow,traces->array[k][j]);
			mpz_mod(t1,t1,P);
			mpz_mod(t2,A->array[i][j],P);
			mpz_mul(t1,t1,t2);
			mpz_mod(t1,t1,P);
			mpz_add(acc,acc,t1);
			mpz_mod(acc,acc,P);
		}

		mpz_set(t1,acc);
		mpz_add(t1, t1, pb_eff_2);
		mpz_div(t1, t1, pb_eff);
		mpz_mod(t1, t1, pdelta);
		if (mpz_cmp(t1,pdelta_2) > 0) {
			mpz_sub(t1, t1, pdelta);
		}

		mpz_set(cutTraces->array[k][i],t1);
	}

	mpz_clear(acc);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(lcpow);
	mpz_clear(lcred);
}

/**
 * Computes the basis matrix for the knapsack problem defined
 * by van Hoeij (2002). This is formed from the lattice
 * basis matrix B_L for the lattice of (powers of) rational
 * factors of the input polynomial, the two-sided cuts of the
 * (possibly transformed) traces cutTraces of the p-adic
 * factors, the excess precision pdelta that determines the
 * precision of the two-sided cuts, the number d of traces
 * used and the number r of p-adic factors.
 *
 * The routine computes and returns a `balancing constant' C
 * to prevent significant size differences between different
 * entries in the basis matrix, and returns the reduction
 * matrix M.
 */
void BuildReductionMatrix(NTL::mat_ZZ& M, long* C, long r, long d, const mpz_t pdelta, const vec_vec_mpz_t* cutTraces, const NTL::mat_ZZ& B_L) {

	long s = B_L.NumRows();

	double dtmp;
	dtmp = (double) d;
	dtmp *= (double) r;
	dtmp /= 2.0;
	dtmp = sqrt(dtmp);
	dtmp += 1;
	*C = (long) dtmp; // sqrt(d*r/2)+1

	M.SetDims(s+d, r+d);
	clear(M);

	mpz_t t1;
	mpz_t t2;
	mpz_t pdelta2;
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(pdelta2);
	mpz_class t3;
	NTL::ZZ t4;

	for (int i = 1; i <= s; ++i)
		for (int j = 1; j <= r; ++j)
			NTL::mul(M(i, j), B_L(i, j), *C);

	mpz_div_2exp(pdelta2,pdelta,1);

	long maxbits = 0;

	for (int i = 1; i <= s; ++i)
		for (int j = 1; j <= d; ++j) {
			mpz_set_ui(t1,0);
			for (int k = 1; k <= r; ++k) {
				mpz_class_set_NTLZZ(t3,B_L(i,k)); // TODO: this is inefficient, need to do conversion only once
				mpz_mul(t2,t3.get_mpz_t(),cutTraces->array[k-1][j-1]);
				mpz_add(t1,t1,t2);
			}

			mpz_mod(t1,t1,pdelta);
			if (mpz_cmp(t1,pdelta2)>0)
				mpz_sub(t1,t1,pdelta);

			maxbits = MAX(maxbits, mpz_sizeinbase(t1,2));

			t3 = mpz_class(t1);
			NTLZZ_set_mpz_class(t4,t3);
			M(i, j+r) = t4;
		}


	t3 = mpz_class(pdelta);
	NTLZZ_set_mpz_class(t4,t3);
	for (int i = 1; i <= d; ++i)
		M(i+s, i+r) = t4;

	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(pdelta2);
}

/**
 * Computes the projection of the reduced basis matrix
 * computed by LLL onto the first r coordinates, generating
 * an updated basis for the the lattice of (powers of)
 * rational factors of the input polynomial, after
 * determining the number of rows that are `short' in
 * the sense of LLL reduction, which is determined from
 * the vector D of sizes of the Gram-Schmidt basis vectors
 * returned by the LLL_plus algorithm.
 *
 * This computation requires the `balancing constant' C,
 * the number r of p-adic factors, the number d of
 * (transformed) traces computed in order to compute the
 * size bound on all vectors in the solution set of the
 * knapsack problem must have. This determines the number
 * of rows in the projected basis.
 */
void ProjectBasis(NTL::mat_ZZ& B1, NTL::vec_ZZ& D, NTL::mat_ZZ& M, long C, long r, long d) {

	long k = M.NumRows();
	mpz_class d1(d);
	mpz_class r1(r);
	mpz_class C1(C);
	mpz_class bound;
	bound = 4*C1*C1*r1 + d1*r1*r1;
	NTL::ZZ bnd;
	NTLZZ_set_mpz_class(bnd,bound);

	while (k >= 1 && 4*D[k] > bnd*D[k-1]) k--;

	NTL::mat_ZZ B2;

	B2.SetDims(k, r);

	for (int i = 1; i <= k; ++i)
		for (int j = 1; j <= r; ++j)
			div(B2(i, j), M(i, j), C);

	M.kill(); // save space
	D.kill();

	NTL::ZZ det2;
	long rnk;

	rnk = image(det2, B2);

	B1.SetDims(rnk, r);
	for (int i = 1; i <= rnk; ++i)
		for (int j = 1; j <= r; ++j)
			B1(i, j) = B2(i + k - rnk, j);
}

/**
 * Computes an integer d together with an n x m matrix R
 * so that R/d is the reduced row echelon form of the
 * input n x m integer matrix M, which is assumed to have
 * linearly independent rows.
 *
 * This routine is probabilistic in the following sense: the
 * result is always correct, but with a negligible
 * probability (specifically, if NTL::GenPrime returns a
 * composite, and the modular NTL::gauss routine can't
 * invert a non-zero element).
 */
void ReducedRowEchelonForm(NTL::ZZ& d_out, NTL::mat_ZZ& R_out, const NTL::mat_ZZ& M) {

	long n = M.NumRows();
	long m = M.NumCols();

	if (n == 0 || m == 0) fprintf(stderr,"BPAS, error: input matrix has a zero dimension.");

	NTL::zz_pBak bak;
	bak.save();

	NTL::mat_ZZ S;
	S.SetDims(n, n);
	NTL::mat_ZZ S_inv;
	NTL::mat_ZZ R;
	NTL::mat_zz_p MM;

	for (;;) {
		long p = NTL::GenPrime_long(NTL_SP_NBITS);
		NTL::zz_p::init(p);

		conv(MM, M);

		long r = gauss(MM);
		if (r < n) continue;

		// compute pos(1..n), so that pos(i) is the index
		// of the i-th pivot column

		NTL::vec_long pos;
		pos.SetLength(n);

		long i, j;
		for (i=j=1; i<=n; ++i) {
			while (MM(i, j) == 0) j++;
			pos(i) = j;
			j++;
		}

		// compute the n x n sub-matrix consisting of the
		// pivot columns of M


		for (i=1; i <= n; ++i)
			for (j=1; j<=n; ++j)
				S(i, j) = M(i, pos(j));

		NTL::ZZ d;

		inv(d, S_inv, S);
		if (d == 0) continue;

		mul(R, S_inv, M);

		// now check that R is of the right form, which it will be
		// if we were not unlucky

		long OK = 1;

		for (i=1; i<=n && OK; ++i) {
			for (j=1; j<pos(i) && OK; ++j)
				if (R(i, j) != 0) OK = 0;

			if (R(i, pos(i)) != d) OK = 0;

			for (j=1; j<i && OK; ++j)
				if (R(j, pos(i)) != 0) OK = 0;
		}

		if (!OK) continue;

		d_out = d;
		R_out = R;
		break;
	}

	MM.kill();
	S.kill();
	S_inv.kill();
	R.kill();
}

//void ReducedRowEchelonForm(NTL::ZZ& d_out, NTL::mat_ZZ& R_out, const NTL::mat_ZZ& M) { // backup

//	long n = M.NumRows();
//	long m = M.NumCols();

//	if (n == 0 || m == 0) fprintf(stderr,"BPAS, error: input matrix has a zero dimension.");

//	NTL::zz_pBak bak;
//	bak.save();

//	for (;;) {
//		long p = NTL::GenPrime_long(NTL_SP_NBITS);
//		NTL::zz_p::init(p);

//		NTL::mat_zz_p MM;
//		conv(MM, M);

//		long r = gauss(MM);
//		if (r < n) continue;

//		// compute pos(1..n), so that pos(i) is the index
//		// of the i-th pivot column

//		NTL::vec_long pos;
//		pos.SetLength(n);

//		long i, j;
//		for (i=j=1; i<=n; ++i) {
//			while (MM(i, j) == 0) j++;
//			pos(i) = j;
//			j++;
//		}

//		// compute the n x n sub-matrix consisting of the
//		// pivot columns of M

//		NTL::mat_ZZ S;
//		S.SetDims(n, n);

//		for (i=1; i <= n; ++i)
//			for (j=1; j<=n; ++j)
//				S(i, j) = M(i, pos(j));

//		NTL::mat_ZZ S_inv;
//		NTL::ZZ d;

//		inv(d, S_inv, S);
//		if (d == 0) continue;

//		NTL::mat_ZZ R;
//		mul(R, S_inv, M);

//		// now check that R is of the right form, which it will be
//		// if we were not unlucky

//		long OK = 1;

//		for (i=1; i<=n && OK; ++i) {
//			for (j=1; j<pos(i) && OK; ++j)
//				if (R(i, j) != 0) OK = 0;

//			if (R(i, pos(i)) != d) OK = 0;

//			for (j=1; j<i && OK; ++j)
//				if (R(j, pos(i)) != 0) OK = 0;
//		}

//		if (!OK) continue;

//		d_out = d;
//		R_out = R;
//		break;
//	}
//}

/**
 * Computes a potential rational factor of the polynomial
 * being factored from the p-adic factors mod P according
 * to the index vector I.
 */
void ComputePotentialRationalFactor(DUZP_t** factor, const vec_DUZP_t* padicFactors, const mpz_t P, const NTL::vec_long& I) {

	DUZP_t* ret = makeConstPolynomial_DUZP(1,1);

	for (int i=0; i<I.length(); ++i) {
		multiplyPolynomials_DUZP_inp(padicFactors->polys[I[i]],&ret);
		applyModuloSymmetric_DUZP_inp(ret,P);
	}

	*factor = ret;
}

/**
 * Checks that the success conditions A and B for van
 * Hoeij's algorithm are met. If both conditions are
 * met then the vector of rational factors of the
 * polynomial f being factored is returned.
 *
 * Condition A is that each column of B_L contains
 * exactly one 1 with all other entries being zero,
 * required if the basis vectors are to determine a
 * rational factor.
 *
 * Condition B is that each row provides the indices
 * of the p-adic factors mod P that when multiplied
 * together and reduced mod P produce a rational
 * factor of f.
 */
long ConditionsAreMet(vec_DUZP_t* factors, const NTL::mat_ZZ& B_L, const vec_DUZP_t* padicFactors, const mpz_t P, const DUZP_t* f, long bound) {

//	gmp_fprintf(stderr,"P = %Zd\n",P);

	// double tt0, tt1;
	NTL::ZZ det;
	NTL::mat_ZZ R;
	long s, r;
	long i,j,cnt;

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"    [CAM] checking A (s = %ld): reduced row echelon form\n",B_L.NumRows());
	#endif

	ReducedRowEchelonForm(det, R, B_L);

	// check if condition A holds

	s = B_L.NumRows();
	r = B_L.NumCols();

	for (j = 0; j<r; ++j) {
		cnt = 0;
		for (i = 0; i<s; ++i) {
			if (R[i][j] == 0) continue;
			if (R[i][j] != det) {
				fprintf(stderr,"    [CAM] non-unity value in basis matrix: Failed :(\n");
				return 0;
			}
			cnt++;
		}

		if (cnt != 1) {
			fprintf(stderr,"    [CAM] multiple ones in basis matrix column: Failed :(\n");
			return 0;
		}
	}

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"    [CAM] Passed A check, checking B...\n");
	#endif

	// extract relevant information from R

	NTL::vec_vec_long I_vec;
	I_vec.SetLength(s);

	NTL::vec_long deg_vec;
	deg_vec.SetLength(s); // TODO: could use this information to allocate space
						  //       for the products and compute in place

	for (i = 0; i<s; ++i) {
		long dg = 0;

		for (j = 0; j<r; ++j) {
			if (R[i][j] != 0) append(I_vec[i], j);
			dg += padicFactors->polys[j]->lt;
		}

		deg_vec[i] = dg;
	}


	R.kill(); // save space

	// check if any candidate factor is the product of too few
	// modular factors

	for (i = 0; i<s; ++i)
		if (I_vec[i].length() <= CARDINALITY_THRESHOLD) {
			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"    [CAM] potential factor produced by too few modular factors: Failed :(\n");
			#endif
			return 0;
		}

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"    [CAM] 1..");
	#endif


	// sort deg_vec, I_vec in order of increasing degree

	for (i=0; i<s-1; ++i)
		for (int j=0; j < s-1-i; ++j)
			if (deg_vec[j] > deg_vec[j+1]) {
				_ntl_swap(deg_vec[j], deg_vec[j+1]);
				swap(I_vec[j], I_vec[j+1]);
			}

	// perform constant term tests

	mpz_t ct;
	mpz_init(ct);
	mpz_mul(ct,f->coefs[f->lt],f->coefs[0]);

	mpz_t half_P;
	mpz_init(half_P);
	mpz_fdiv_q_2exp(half_P,P,1);

	mpz_t lc;
	mpz_t prod;
	mpz_init(lc);
	mpz_init(prod);
	mpz_mod(lc,f->coefs[f->lt],P);

	mpz_t t1;
	mpz_init(t1);

	int exit = 0;

	for (i=0; i<s; ++i) {
		NTL::vec_long& I = I_vec[i];
		mpz_set(prod,lc);
		for (j=0; j<I.length(); ++j) {
			mpz_mul(prod,prod,(padicFactors->polys[I[j]])->coefs[0]);
			mpz_mod(prod,prod,P);
		}

		mpz_mod(t1,prod,P);
		if (mpz_cmp(t1,half_P)>0)
			mpz_sub(t1,t1,P);

		if (!mpz_divisible_p(ct, t1)) {
			exit = 1;
			break;
		}
	}

	mpz_clear(ct);
	mpz_clear(half_P);
	mpz_clear(prod);
	mpz_clear(t1);

	if (exit) {
		fprintf(stderr,"\n    [CAM] constant term division test failed: Failed :(\n");
		mpz_clear(lc);
		return 0;
	}

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"2..");
	#endif

	// multiply out polynomials and perform size tests

	vec_DUZP_t fac;
	vec_DUZP_t_init(&fac);
	DUZP_t* g;

	for (i=0; i<s-1; ++i) {
		NTL::vec_long& I = I_vec[i];
		ComputePotentialRationalFactor(&g,padicFactors,P,I);
		multiplyByInteger_DUZP_inp(g,lc);
		applyModuloSymmetric_DUZP_inp(g,P);

		if (MaxBits(g) > bound) {
			fprintf(stderr,"\n    [CAM] potential factor violates bound: Failed :(\n");
			freePolynomial_DUZP(g);
			mpz_clear(lc);
			return 0;
		}

		primitivePart_DUZP_inp(g);
		vec_DUZP_t_push(&fac,g);
		freePolynomial_DUZP(g);
	}

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"3..");
	#endif

////	// finally...trial division

	DUZP_t* f1 = deepCopyPolynomial_DUZP(f);
	DUZP_t* h;

	for (i=0; i<s-1; ++i) {
		if (!divideTest_DUZP(f1, fac.polys[i], &h)) {
			fprintf(stderr,"\n    [CAM] division test failed: Failed :(\n");
			mpz_clear(lc);
			return 0;
		}

		freePolynomial_DUZP(f1);
		f1 = h;
	}

	// got them!

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"done!\n");
	#endif

	for (i=0; i<fac.size; ++i)
		vec_DUZP_t_push(factors,fac.polys[i]);
	vec_DUZP_t_push(factors, f1);

	freePolynomial_DUZP(f1);
	vec_DUZP_t_free(&fac);

	mpz_clear(lc);

	return 1;
}

/**
 * Continues the lifting process for a set of p-adic
 * factors mod P1 = p^e1 of f up to a new lifting
 * bound. The new lifting exponent e1 and corresponding
 * power P1 are returned by over-writing the supplied
 * values.
 *
 * Currently, linear lifting is being used, requiring
 * that the original modular factors and their sigmas
 * are also passed to the routine. The p-adic sigmas,
 * required for the quadratic lifting are not currently
 * being used.
 */
void AdditionalLifting(mpz_t* P1, int* e1, vec_DUZP_t* padicFactors, DUZP_t** sigmas, long p, int new_bound, const DUZP_t* f, long doubling, const factors_t* modularFactors, duspoly_t** sigmas2) {

	int new_e1;

	if (doubling)
		new_e1 = MAX(2*(*e1), new_bound); // at least double e1
	else
		new_e1 = new_bound;

	fprintf(stderr,"    [AL] additional hensel lifting from %d to %d...\n",*e1,new_e1);

	mpz_t new_P1;
	mpz_init(new_P1);

	mpz_ui_pow_ui(new_P1, p, new_e1);

	// make f monic mod P1 before lifting

	Prime_ptr* Pptr = smallprimefield_get_prime_constants(p);
	int nf = padicFactors->size;

//	DUZP_t** liftedPolys = padicFactors->polys;

	mpz_t(targetBound);
	mpz_init_set_ui(targetBound,p);
	mpz_pow_ui(targetBound,targetBound,new_e1);
	fprintf(stderr,"    [AL] modularFactors:\n");
	for (int i=0; i<nf; ++i)
		printPolynomialOutForm_spX(modularFactors->polys[i],Pptr);
	fprintf(stderr,"    [AL] padicFactors->polys:\n");
	for (int i=0; i<nf; ++i)
		printPoly_DUZP(padicFactors->polys[i],x);
	fprintf(stderr,"    [AL] sigmas:\n");
	for (int i=0; i<nf; ++i)
		printPolynomialOutForm_spX(sigmas2[i],Pptr);
	gmp_fprintf(stderr,"    [AL] targetBound = %Zd\n",targetBound);
	fprintf(stderr,"    [AL] e1 before lift: %d\n",*e1);
	fprintf(stderr,"    [AL] f = ");
	printPoly_DUZP(f,x);


	fprintf(stderr,"    [AL] resuming lift...");
	*e1 = multiTermPadicLiftResume(f, CONSTCONSTCAST(duspoly_t,modularFactors->polys), padicFactors->polys, sigmas2, nf, *e1, targetBound, Pptr);
//	*e1 = multiTermQuadraticPadicLiftResume(f, liftedPolys, sigmas, nf, *e1, targetBound, Pptr);
	fprintf(stderr,"done.\n");

//	for (int i=0; i<nf; ++i)
//		freePolynomial_DUZP(padicFactors->polys[i]);
//	free(padicFactors->polys);
//	padicFactors->polys = liftedPolys;

	mpz_clear(targetBound);
	mpz_set(*P1,new_P1);
	mpz_clear(new_P1);
	free(Pptr);
	*e1 = new_e1;

}

/**
 * Computes the recombination step for univariate
 * factorization following the algorithm of van Hoeij
 * (2002) and the implementation strategy of the NTL
 * library.
 *
 * The routine takes as input the polynomial f_in to
 * be factored, the factorization of f_in modularFactors
 * computed modulo p, the vector of p-adic factors
 * lifted modulo P_in = p^{e_in}, the lifting bound, as
 * well as the linear (sigmas2) and quadratic (sigmas)
 * coefficients in the multi-term diophantine problem
 * solved in the lifting process.
 *
 * Currently linear lifting is being used, so sigmas
 * is not being used.
 */
void DUZP::Factoring::Recombine(vec_DUZP_t* rationalFactors, const DUZP_t* f_in, const vec_DUZP_t* padicFactors_in, DUZP_t** sigmas, const mpz_t P_in, long_t p, int e_in, long_t bound, const factors_t* modularFactors, duspoly_t** sigmas2) {

	if (rationalFactors == NULL)
		return;

	vec_DUZP_t_clear(rationalFactors);

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"   [RECOMB] Recombine:\n");
	#endif

	long r = padicFactors_in->size;

	mpz_t P;
	mpz_init(P);
	mpz_set(P,P_in);
	int e = e_in;

	DUZP_t* f = deepCopyPolynomial_DUZP(f_in);

	vec_DUZP_t padicFactors;
	vec_DUZP_t_init2(&padicFactors,padicFactors_in->size);
	for (int i=0; i<padicFactors_in->size; ++i)
		vec_DUZP_t_push(&padicFactors,padicFactors_in->polys[i]);
//	vec_DUZP_t_init_set(&padicFactors,padicFactors_in->size);
//	for (int i=0; i<padicFactors_in->size; ++i) {
//		padicFactors.polys[i] = deepCopyPolynomial_DUZP(padicFactors_in->polys[i]);
//	}
//	fprintf(stderr,"   [RECOMB] init: padicFactors.size = %ld, padicFactors.alloc = %ld:\n",padicFactors.size,padicFactors.alloc);

	vec_DUZP_t_clear(rationalFactors);
//	vec_DUZP_t_free(rationalFactors);
//	vec_DUZP_t_init(rationalFactors);

//	static int initRand = 0;
//	static gmp_randstate_t R_STATE;
//	if (!initRand) {
//		time_t t = time(NULL);
//
//		gmp_randinit_default (R_STATE);
//		gmp_randseed_ui(R_STATE, t);

//		// fprintf(stderr, "seed: %lu\n", t);
//		initRand = 1;
//	}
//	mpz_t fiftyfifty;
//	mpz_init(fiftyfifty);

	long m = 1;
	// TODO: at some point, enable NaiveFactorSearch as THE recombination method when neither NTL
	//       nor Maple are installed.
	while (2*m <= padicFactors.size && (m <= CARDINALITY_THRESHOLD || padicFactors.size <= SIZE_THRESHOLD)) {

		if (padicFactors.size == 0)
			break;

		NaiveFactorSearch(&f, rationalFactors, &padicFactors, P, p, e, m, bound);
//		fprintf(stderr,"   [RECOMB] after NFS: padicFactors.size = %ld, padicFactors.alloc = %ld:\n",padicFactors.size,padicFactors.alloc);
		m++;

	}

	if (2*m > padicFactors.size) {

		if (f->lt > 0)
			vec_DUZP_t_push(rationalFactors,f);

	}
	else {

		#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr,"   [RECOMB] Beginning van Hoeij Algorithm...\n");
		#endif

		long n = f->lt;

		r = padicFactors.size;

		long d = 0;				// number of traces to consider
		long deltaInBits = 0;	// excess precision in bits max(a_i - b_i)

		long delta = 0;		// excess precision max(a_i - b_i)
		mpz_t pdelta;
		mpz_init(pdelta);
		mpz_set_ui(pdelta,1);

		NTL::vec_long b;	// b_i
		vec_mpz_t pb;		// p^(b_i)
		vec_mpz_t_init2(&pb,0);

		vec_vec_mpz_t traces;		// vector of vectors of traces
		vec_vec_mpz_t cutTraces;	// vector of vectors of two-sided cut traces
		vec_vec_mpz_t_init(&traces,r);
		vec_vec_mpz_t_init(&cutTraces,r);

		mpz_t RB;
		mpz_init(RB);
		computeRootBound(&RB,f);

		NTL::mat_ZZ B_L;			// factor exponent lattice basis
		ident(B_L, r);

		// potentially needed variables
		long dense = 0;
		long ran_bits = 32;
		long s = r;

		while (1) {

			// potentially needed variables
			long d_last;
			long d_inc;
			long d_index;

			d_last = d;

//			fprintf(stderr,"d_last = %ld, d_inc = %ld\n",d_last,d_inc);

			// NTL's Optimization Magix! //

			// set d_inc:

			if (!dense) {
				d_inc = 1 + d/8;
			}
			else {
				d_inc = 1 + d/4;
			}

			d_inc = MIN(d_inc, n-1-d);

//			fprintf(stderr,"d_inc = %ld\n",d_inc);

			d += d_inc;

			// set bit_delta:

			if (deltaInBits == 0) { // set initial value...don't make it any smaller than 2*r

				deltaInBits = 2*r;
			}
			else {

				long extra_bits;

				if (!dense) {
					extra_bits = 1 + deltaInBits/8;
				}
				else if (d_inc != 0) {

					if (d1_value(deltaInBits, r, s) > 1)
						extra_bits = 1 + deltaInBits/16;
					else
						extra_bits = 0;

				}
				else
					extra_bits = 1 + deltaInBits/8;

				deltaInBits += extra_bits;
			}

			///////////////////////////////

//			fprintf(stderr,"deltaInBits = %ld\n",deltaInBits);

			if (d > d1_value(deltaInBits, r, s))
				dense = 1;

//			fprintf(stderr,"dense = %ld\n",dense);

			Computepdelta(&delta, &pdelta, p, deltaInBits);

//			gmp_fprintf(stderr,"pdelta = %Zd\n",pdelta);

			long d1;
			long b_eff;
			mpz_t pb_eff;
			mpz_init(pb_eff);

			if (!dense) {
				for (d_index = d_last+1; d_index <= d; d_index++) {
//					fprintf(stderr,"d_index = %ld, p = %Ld\n",d_index,p);
					Computepb(b, &pb, p, d_index, RB, n);
				}

//				gmp_fprintf(stderr,"d = %ld, RB = %Zd\n",d,RB);

				d1 = d;
				b_eff = b(d);
//				b_eff = b[d-1];
				mpz_set(pb_eff,pb.array[d-1]);
//				gmp_fprintf(stderr,"pb_eff = %Zd\n",pb_eff);
			}
			else {
				d1 = d1_value(deltaInBits, r, s);
//				fprintf(stderr,"d1 = %ld\n",d1);
				Computepb_eff(&b_eff, &pb_eff, p, d, RB, n, ran_bits);
//				gmp_fprintf(stderr,"pb_eff = %Zd\n",pb_eff);
			}

//			gmp_fprintf(stderr,"b_eff = %ld, pb_eff = %Zd\n",b_eff,pb_eff);

//			gmp_fprintf(stderr,"d = %ld, s = %ld, delta = %ld, b_eff = %ld\n",d,s,delta,b_eff);

//			if (dense)
//				fprintf(stderr,"d1 = %ld\n",d1);

//			fprintf(stderr,"b_eff + delta > e = %d\n",(b_eff + delta > e));
			if (b_eff + delta > e) {
				long doubling;

				// TODO: decide whether to implement this (equivalent to at least one quadratic lift).
				doubling = 1;

//				fprintf(stderr,"    [AL] padicFactors:\n");
//				for (int i=0; i<padicFactors.size; ++i)
//					printPoly_DUZP(padicFactors.polys[i],x);

				AdditionalLifting(&P,&e,&padicFactors,sigmas,p, b_eff + delta,f,doubling,modularFactors,sigmas2);
//				fprintf(stderr,"   [RECOMB] after: AL padicFactors.size = %ld, padicFactors.alloc = %ld:\n",padicFactors.size,padicFactors.alloc);

				#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
					fprintf(stderr,"   [RECOMB] recomputing traces:\n");
				#endif

				vec_vec_mpz_t_set(&traces,r,d_last);

//				fprintf(stderr,"r = %ld, d_last = %ld\n",r,d_last);
//				fprintf(stderr,"padicFactors:\n");
//				for (int i=0; i<r; ++i)
//					printPoly_DUZP(padicFactors.polys[i],x);

				for (int i=0; i<r; i++) {

					for (d_index = 1; d_index <= d_last; d_index++) {
//						fprintf(stderr,"computing traces[%d][%ld]...",i,d_index-1);
						ComputeTrace(&traces,i,padicFactors.polys[i],d_index,P);
//						fprintf(stderr,"done.\n");
					}
				}
//				fprintf(stderr,"   [RECOMB] after CT: padicFactors.size = %ld, padicFactors.alloc = %ld:\n",padicFactors.size,padicFactors.alloc);
			}

			vec_vec_mpz_t A;
			vec_vec_mpz_t_init(&A,0);
			NTL::ZZ zzt;

			if (dense) {
				#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
					fprintf(stderr,"   [RECOMB] computing random matrix A\n");
				#endif
				vec_vec_mpz_t_set(&A,d1,d);
				for (int i=0; i<d1; ++i) {
					for (int j=0; j<d; ++j) {
					    RandomBits(zzt, ran_bits); // NTL (Not So) Super Secure Random Bit magix!
					    mpz_t_set_NTLZZ(A.array[i][j],zzt);
//						mpz_urandomb(A.array[i][j],R_STATE,ran_bits);
//						mpz_urandomb(fiftyfifty,R_STATE,1);
////						gmp_fprintf(stderr,"\nfiftyfifty = %Zd\n",fiftyfifty);
//						if (mpz_cmp_ui(fiftyfifty,1) == 0)
//							mpz_neg(A.array[i][j],A.array[i][j]);

//						gmp_fprintf(stderr,"%Zd ",A.array[i][j]);
					}
//					fprintf(stderr,"\n");
				}
			}

			vec_vec_mpz_t_set(&traces,r,d);
			vec_vec_mpz_t_set(&cutTraces,r,d1);

			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"   [RECOMB] computing traces\n");
			#endif
			for (int i=0; i<r; i++) {

				for (d_index=d_last+1; d_index<=d; d_index++) {
					ComputeTrace(&traces,i,padicFactors.polys[i],d_index,P);
				}

//				#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
//					fprintf(stderr,"   [RECOMB] computing two-sided cut traces:\n");
//				#endif
				if (!dense) {
//					fprintf(stderr," sparse case\n");
					CutTraces(&cutTraces,&traces,i,d,pb,pdelta,P,f->coefs[f->lt]);
				}
				else {
//					fprintf(stderr," dense case\n");
					DenseCutTraces(&cutTraces,&traces,i,d,d1,pb_eff,pdelta,P,f->coefs[f->lt],&A);
				}
			}

			mpz_clear(pb_eff);
			vec_vec_mpz_t_free(&A);

			NTL::mat_ZZ M;
			long C;

//			fprintf(stderr,"building reduction matrix\n");
			BuildReductionMatrix(M, &C, r, d1, pdelta, &cutTraces, B_L);

			// TODO: decide whether to implement this NTL optimization
//			if (SkipSparse) {
//				if (!dense) {
//					fprintf(stderr,"skipping LLL\n");
//						continue;
//				}
//			}

			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"   [RECOMB] computing lattice basis reduction\n");
			#endif
			NTL::vec_ZZ D;
			long rnk = LLL_plus(D, M);

			if (rnk != s + d1) {
				fprintf(stderr,"BPAS, error: non-sensical rank from LLL in Recombine\n");
				exit(1);
			}

			NTL::mat_ZZ B1;

			ProjectBasis(B1, D, M, C, r, d1);

			if (B1.NumRows() >= s) continue;
			// no progress...try again

			// otherwise, update B_L and test if we are done

			swap(B1, B_L);
			B1.kill();
			s = B_L.NumRows();

			if (s == 0) {
				fprintf(stderr,"BPAS, error: lattice basis cannot have dimension 0!\n");
				exit(1);
			}

			if (s == 1) {
				vec_DUZP_t_push(rationalFactors, f);
				break;
			}

			if (s > r / (CARDINALITY_THRESHOLD + 1)) continue;
			// dimension too high...we can't be done

			if (ConditionsAreMet(rationalFactors,B_L,&padicFactors,P,f,bound)) break;

		}

		mpz_clear(pdelta);
		mpz_clear(RB);
		vec_mpz_t_free(&pb);
		vec_vec_mpz_t_free(&traces);
		vec_vec_mpz_t_free(&cutTraces);

	}

//	fprintf(stderr,"   [RECOMB] free: padicFactors.size = %ld, padicFactors.alloc = %ld:\n",padicFactors.size,padicFactors.alloc);
	freePolynomial_DUZP(f);
	mpz_clear(P);
	vec_DUZP_t_free(&padicFactors);
//	mpz_clear(fiftyfifty);

}

/**
 * Computes the degree pattern of a set of modular
 * factors computed using distinct degree
 * factorization (where the exponents indicate the
 * degree of the ultimate modular factors and the
 * polynomials are products of modular factors of
 * the given degree).
 *
 * The degree pattern is encoded so that pat[i]
 * gives the number of modular factors of degree i.
 */
void RecordDegreePattern(NTL::vec_long& pat, factors_t* modularFactors) {

	long n = pat.length();

	for (int i=0; i<n; ++i)
		pat[i] = 0;

	long k = modularFactors->alloc;
	long d, m;

	for (int i=0; i<k; ++i) {
		d = modularFactors->exps[i];
		m = (modularFactors->polys[i]->lt)/d;

		pat[d] += m;
	}
}

/**
 * Computes the degree pattern of a set of modular
 * factors computed using distinct degree
 * factorization (where the exponents indicate the
 * degree of the ultimate modular factors and the
 * polynomials are products of modular factors of
 * the given degree).
 *
 * The degree pattern is encoded so that pat[i]
 * gives the number of modular factors of degree i.
 */
void RecordDegreePattern(NTL::vec_long& pat, factsll_t* modularFactors) {

	long n = pat.length();

	for (int i=0; i<n; ++i)
		pat[i] = 0;

	long d, m;
	factsll_t* currFactor = modularFactors;

	while (currFactor != NULL) {
		d = currFactor->exp;
		m = (currFactor->poly->lt)/d;

		pat[d] += m;

		currFactor = currFactor->next;
	}
}

/**
 * Computes the number of modular factors given
 * the degree pattern obtained from distinct
 * degree factorization.
 */
long NumberOfFactors(const NTL::vec_long& pat) // n = len(pat)
{
	long n = pat.length();
	long res = 0;

	for (int i=0; i<n; ++i)
		res += pat[i];

	return res;
}

/**
 * Computes the possible degree of rational factors
 * obtained from modular factors given the degree
 * pattern pat obtained from distinct degree
 * factorization. This is used to detect irreducible
 * polynomials early in the factorization process.
 */
void ComputePossibleDegree(mpz_t* degs, const NTL::vec_long& pat) {

	long n = pat.length();

	mpz_set_si(*degs,1);
	mpz_t tmp;
	mpz_init(tmp);

	for (int deg=1; deg<n; ++deg) {
		for (int i=0; i<pat[deg]; ++i) {
			mpz_mul_2exp(tmp,*degs,deg);
			mpz_ior(*degs,*degs,tmp);
		}
	}

	mpz_clear(tmp);
}

/**
 * Computes the factorization of the input polynomial f
 * modulo a number of small primes (small in the sense
 * of starting from 3 and working upwards) and computes
 * the factorization of modulo the prime that has the
 * smallest number of factors of the primes tried.
 *
 * The input f is partially factored using distinct
 * degree factoriztion if the prime does not divide
 * lc(f) and f is squarefree modulo that prime. Of the
 * primes for which distinct degree factorization is
 * performed (numbering BPAS_SPFACT_INIT_NUM_PRIMES),
 * the prime with the smallest number of modular factors
 * is the factored fully using same degree factorization.
 * The result of this computation is returned.
 *
 * Degree pattern computations are used for early
 * detection of irreducible input.
 */
long SmallPrimeFactorize(factors_t** modularFactors, factoring_info_t* factoringInfo, const DUZP_t* f) {

	if (modularFactors == NULL) {
		fprintf(stderr,"BPAS, error: NULL pointer passed to SmallPrimeFactorization, cannot return result\n");
		exit(1);
	}

	// some sanity checking...

	if (BPAS_SPFACT_INIT_NUM_PRIMES < 1 || BPAS_SPFACT_INIT_NUM_PRIMES > 10000) {
		fprintf(stderr,"BPAS, error: bad value for BPAS_SPFACT_INIT_NUM_PRIMES: %d\n",BPAS_SPFACT_INIT_NUM_PRIMES);
		exit(1);
	}

//	if (BPAS_SPFACT_MAX_NUM_PRIMES < BPAS_SPFACT_INIT_NUM_PRIMES || BPAS_SPFACT_MAX_NUM_PRIMES > 10000) {
//		fprintf(stderr,"BPAS, error: bad value for BPAS_SPFACT_MAX_NUM_PRIMES: %d\n",BPAS_SPFACT_MAX_NUM_PRIMES);
//		exit(1);
//	}

	factoringInfo->p.SetLength(BPAS_SPFACT_INIT_NUM_PRIMES);
	factoringInfo->pattern.SetLength(BPAS_SPFACT_INIT_NUM_PRIMES);

	// set bits 0..n of factoringInfo->PossibleDegrees
	mpz_setbit(factoringInfo->possibleDegrees,f->lt+1);
	mpz_sub_ui(factoringInfo->possibleDegrees,factoringInfo->possibleDegrees,1);

	mpz_t lc;
	mpz_init_set(lc,f->coefs[f->lt]);

	Prime_ptr* Pptr;

	long num_primes = 0;
	factoringInfo->num_primes = num_primes;
	factoringInfo->Pptr = smallprimefield_get_prime_constants(2);

	factsll_t* thisFactorsll; // = initFactsll_spX(1);
	factsll_t* bestFactorsll = initFactsll_spX(1);

	factors_t* outputFactors;

	long long int modlc;
	long p = factoringInfo->s.next();
	long degf = f->lt;
	long bestp_index;

	long min_k = degf+1;
	long irred = 0;

	mpz_t pd;
	mpz_init(pd);

	duspoly_t* ff;
	duspoly_t* ffp;
	duspoly_t* g;

	for (; num_primes<BPAS_SPFACT_INIT_NUM_PRIMES ;) {

		p = factoringInfo->s.next();

		if (!p) {
			fprintf(stderr,"BPAS, error: we've run out of small primes! :o\n");
			exit(1);
		}
		if (mpz_divisible_ui_p(lc,p)) {
			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"   [SPF] skipping %ld\t(lc test)\n",p);
			#endif
			continue;
		}


		Pptr = smallprimefield_get_prime_constants(p);

		ff = convertToDUSP_DUZP(f,Pptr);

		// TODO: use the inp version of this function
		monicPolynomialInForm_spX (ff,&ffp,&modlc,Pptr); // make ff monic
		freePolynomial_spX(&ff);
		ff = ffp;

		derivativePolyInForm_spX (ff, &ffp, Pptr);
		plainGCDInForm_spX (ff, ffp, &g, Pptr); // check ff is squarefree

		if (!isOneInForm_spX (g,Pptr)) {
			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"   [SPF] skipping %ld\t(sqf test)\n",p);
			#endif
			free(Pptr);
			freePolynomial_spX(&ff);
			freePolynomial_spX(&ffp);
			freePolynomial_spX(&g);
			continue;
		}

		#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr,"   [SPF] distinct degree factorization mod %ld...",p);
		#endif

		distinctDegFactorizationInFormll1_spX (ff, &thisFactorsll, Pptr);

		#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr,"done.\n");
		#endif

		factoringInfo->p[num_primes] = p;

		NTL::vec_long& pattern = factoringInfo->pattern[num_primes];
		pattern.SetLength(degf+1);

		RecordDegreePattern(pattern, thisFactorsll);
		long k = NumberOfFactors(pattern);

		// fprintf(stderr,"   [SPF] degree sequence: ");
		// for (int i=0; i<=degf; ++i) {
		// 	if (pattern[i])
		// 		fprintf(stderr,"%ld*%d ",pattern[i],i);
		// }
		// fprintf(stderr,"\n");

		if (k == 1) {
			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"   [SPF] only one modular factor; irreducible polynomial found!\n");
			#endif
			irred = 1;
			free(Pptr);
			freePolynomial_spX(&ff);
			freePolynomial_spX(&ffp);
			freePolynomial_spX(&g);
			freeFactsll_spX(&thisFactorsll);
			break;
		}

//       update admissibility info

		ComputePossibleDegree(&pd,pattern);
		mpz_and(factoringInfo->possibleDegrees,factoringInfo->possibleDegrees,pd);

		if (mpz_popcount(factoringInfo->possibleDegrees) == 2) {
			#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
				fprintf(stderr,"   [SPF] only one nontrivial possible degree; irreducible polynomial found!\n");
			#endif
			irred = 1;
			free(Pptr);
			freePolynomial_spX(&ff);
			freePolynomial_spX(&ffp);
			freePolynomial_spX(&g);
			freeFactsll_spX(&thisFactorsll);
			break;
		}

		if (k < min_k) {

			min_k = k;
			freeFactsll_spX(&bestFactorsll);
			bestFactorsll = thisFactorsll;
			free(factoringInfo->Pptr);
			factoringInfo->Pptr = Pptr;
			bestp_index = num_primes;
		}
		else {
			freeFactsll_spX(&thisFactorsll);
			free(Pptr);
		}

		num_primes++;

		freePolynomial_spX(&ff);
		freePolynomial_spX(&ffp);
		freePolynomial_spX(&g);

	}


	if (!irred) {
		// complete factorization

		// TODO: decide whether we need to do this accounting from NTL: remove best prime from LocalInfo
//			swap(LocalInfo.pattern[bestp_index], LocalInfo.pattern[NumPrimes-1]);
//			LocalInfo.p[bestp_index] = LocalInfo.p[NumPrimes-1];
//			NumPrimes--;

		#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr,"   [SPF] completing equal degree factorization mod %Ld...",factoringInfo->Pptr->prime);
		#endif

		equalDegFactorsInFormll_spX (bestFactorsll, &thisFactorsll, factoringInfo->Pptr);

		#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
			fprintf(stderr,"done.\n");
		#endif

		outputFactors = convertToFactorsList_spX (thisFactorsll, factoringInfo->Pptr);
		freeFactsll_spX(&thisFactorsll);
		*modularFactors = outputFactors;
	}

	factoringInfo->num_primes = num_primes;

	mpz_clear(pd);
	mpz_clear(lc);
	freeFactsll_spX(&bestFactorsll);
	return irred;

}


/**
 * The core factorization routine for a primitive and
 * squarefree input polynomial f_in(X) that has a
 * positive leading coefficient.
 *
 * This routine checks for two special cases:
 *   - X is a factor
 *   - X-1 is a factor
 * and removes these initially before proceeding.
 * It then computes a small prime factorization of the
 * input using the irreducibility detection of the
 * SmallPrimeFactorization routine.
 * If this leads to multiple modular factors, these
 * factors are lifted according to the supplied bound,
 * provided the supplied bound meets or exceeds the
 * Landau-Mignotte bound. Currently this process is
 * using linear Hensel lifting.
 * After the lifting process, the p-adic factors mod
 * a power of the small prime are recombined using
 * the van Hoeij method implemented as Recombine.
 */
void factor_prim_sqf_inner(vec_DUZP_t* factors, const DUZP_t* f_in, long bnd)
{
// input is primitive and square-free, with positive leading coefficient
	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"  [FPSI] factor_prim_sqf_inner:\n");
	#endif

	if (factors == NULL)
		return;

	vec_DUZP_t_clear(factors);

	if (f_in->lt <= 1) {
		vec_DUZP_t_push(factors,f_in);
//		fprintf(stderr,"factor_prim_sqf_inner: trivial case 1.\n");
		return;
	}

	// remove a factor of X, if necessary

	DUZP_t* f;
	DUZP_t* r;
	DUZP_t* ex = makePolynomial_DUZP(2);
	DUZP_t* tmp;
	mpz_init(ex->coefs[1]);
	mpz_set_ui(ex->coefs[1],1);
	ex->lt = 1;
	DUZP_t* exm1 = deepCopyPolynomial_DUZP(ex);
	mpz_set_ui(exm1->coefs[0],-1);
	long xfac;

	if (mpz_cmp_ui(f_in->coefs[0],0) == 0) {
		DUZP_t* r;

		dividePolynomials_DUZP(f_in, ex, &f, &r);
		freePolynomial_DUZP(r);
		xfac = 1;
	}
	else {
		f = deepCopyPolynomial_DUZP(f_in);
		xfac = 0;
	}

	// return a factor of X-1 if necessary

	long x1fac = 0;

	mpz_t c1;
	mpz_init(c1);

	mpz_set_ui(c1,0);
	for (int i=0; i<f->lt; ++i)
		mpz_add(c1,c1,f->coefs[i]);

	if (c1 == 0) {
		x1fac = 1;
		mpz_set_ui(exm1->coefs[0],-1);
		dividePolynomials_DUZP(f, exm1, &tmp, &r);
		freePolynomial_DUZP(f);
		freePolynomial_DUZP(r);
		f = tmp;

	}

	mpz_set_ui(c1,0);
	for (int i=0; i<f->lt; ++i)
		mpz_add(c1,c1,f->coefs[i]);

	if (f->lt <= 1) {
		vec_DUZP_t_set(factors,0);
		if (f->lt > 0) {
			vec_DUZP_t_push(factors,f);
		}
		if (xfac) {
			vec_DUZP_t_push(factors,ex);
		}

		if (x1fac) {
			vec_DUZP_t_push(factors,exm1);
		}

//		fprintf(stderr,"factor_prim_sqf_inner: trivial case 2.\n");

		mpz_clear(c1);
		freePolynomial_DUZP(f);
		freePolynomial_DUZP(ex);
		freePolynomial_DUZP(exm1);
		return;
	}

	// TODO: reverse f if this makes lead coefficient smaller

	mpz_t t1;
	mpz_t t2;
	mpz_init(t1);
	mpz_init(t2);

	mpz_abs(t1,f->coefs[f->lt]);
	mpz_abs(t2,f->coefs[0]);

	// TODO: decide if we want to implement reversal strategy
//	if (mpz_cmp(t1,t2) > 0) {
////		inplace_rev(f);
//		rev = 1;
//	}
//	else
//		rev = 0;

	// obtain factorization modulo small primes

	factoring_info_t factoringInfo;
	factoring_info_t_init(&factoringInfo);
	factors_t* modularFactors;

	int irred = SmallPrimeFactorize(&modularFactors, &factoringInfo, f);

	if (irred) {
		// f was found to be irreducible

//		if (rev)
//			inplace_rev(f);

		vec_DUZP_t_push(factors,f);

		if (xfac) {
			vec_DUZP_t_push(factors,ex);
		}

		if (x1fac) {
			vec_DUZP_t_push(factors,exm1);
		}

		mpz_clear(c1);
		mpz_clear(t1);
		mpz_clear(t2);
		freePolynomial_DUZP(f);
		freePolynomial_DUZP(ex);
		freePolynomial_DUZP(exm1);
		return;
	}

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"  [FPSFI] number of modular factors = %ld\n",modularFactors->alloc);
	#endif

	// prepare for Hensel lifting

	// first, calculate bit bound

	long bnd1;
	long n = f->lt;
	long i = 0;
	long e;
	mpz_t P;
	mpz_init(P);
	long p = factoringInfo.Pptr->prime;

	bnd1 = MaxBits(f) + (NumBits(n+1)+1)/2;
//	fprintf(stderr,"bnd1 = %ld\n",bnd1);

	if (!bnd || bnd1 < bnd)
		bnd = bnd1;
	//	fprintf(stderr,"bnd = %ld\n",bnd);

	// TODO: use PossibleDegrees information here
//	i = n/2;
//	while (!bit(LocalInfo.PossibleDegrees, i))
//		i--;

	// NTL Coefficient-Setting Magic! //

	long lc_bnd = NumBits(f->coefs[f->lt]);
//	fprintf(stderr,"lc_bnd = %ld\n",lc_bnd);

	long coeff_bnd = bnd + lc_bnd + i;

	long lift_bnd;

	lift_bnd = coeff_bnd + 15;
//	fprintf(stderr,"coeff_bnd + 15 = %ld\n",lift_bnd);
	// +15 helps avoid trial divisions...can be any number >= 0

	lift_bnd = MAX(lift_bnd, bnd + lc_bnd + 2*NumBits(n) + DUZP_OVERLIFT);
//	fprintf(stderr,"max(lift_bnd, bnd + lc_bnd + 2*NumBits(n) + ZZX_OVERLIFT) = %ld\n",lift_bnd);
	// facilitates "n-1" and "n-2" tests

	lift_bnd = MAX(lift_bnd, lc_bnd + NumBits(c1));
//	fprintf(stderr,"max(lift_bnd, lc_bnd + NumBits(c1)) = %ld\n",lift_bnd);
	// facilitates f(1) test

	lift_bnd += 2;
	// +2 needed to get inequalities right
//	fprintf(stderr,"lift_bnd = %ld\n",lift_bnd);

	////////////////////////////////////


	double dt = (double) lift_bnd;
	dt /= log2((double) p);
	e = (long) floor(dt);
	mpz_ui_pow_ui(P,p,e);

	while (NumBits(P) <= lift_bnd) {
		mpz_mul_ui(P, P, p);
		e++;
	}
//	fprintf(stderr,"lifting bound = %ld bits\n",lift_bnd);

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"  [FPSFI] Hensel lifting to exponent %ld\n",e);
	#endif

	// third, compute f1 so that it is monic and equal to f mod P

	DUZP_t* f1;

	if (mpz_cmp_si(f->coefs[f->lt],1) == 0) {
		f1 = deepCopyPolynomial_DUZP(f);
	}
	else if (mpz_cmp_si(f->coefs[f->lt],-1) == 0) {
		f1 = deepCopyPolynomial_DUZP(f);
		negatePolynomial_DUZP_inp(f1);
	}
	else {
		mpz_mod(t1,f->coefs[f->lt],P);
		if (mpz_sgn(P) < 0) {
			fprintf(stderr,"BPAS, error: what what?! negative prime power?!\n");
			exit(1);
		}
		mpz_invert(t1,t1,P);
		f1 = makePolynomial_DUZP(n+1);
		for (i = 1; i <=n; ++i)
			mpz_init(f1->coefs[i]);
		f1->lt = n;
		for (i = 0; i <= n; ++i) {
			mpz_mul(t2,f->coefs[i],t1);
			mpz_fdiv_r(f1->coefs[i],t2,P);
		}
	}


	// Do Hensel lift

	int nf = modularFactors->alloc;
	vec_DUZP_t padicFactors;
	vec_DUZP_t_init_set(&padicFactors,nf);
	for (int i=0; i<nf; ++i) {
		freePolynomial_DUZP(padicFactors.polys[i]);
		padicFactors.polys[i] = NULL;
	}
	DUZP_t** sigmas = NULL;
	long long int p1 = factoringInfo.Pptr->prime;
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(p);
	mpz_t targetBound;
	mpz_init_set(targetBound,P);

	duspoly_t** sigmas2 = NULL;

//	int e1 = multiTermQuadraticPadicLiftStart(f1, CONSTCONSTCAST(duspoly_t,modularFactors->polys), padicFactors.polys, &sigmas, nf, targetBound, Pptr);
	int e1 = multiTermPadicLiftStart(f1, CONSTCONSTCAST(duspoly_t,modularFactors->polys), padicFactors.polys, &sigmas2, nf, targetBound, Pptr);

	// search for true factors

	// TODO: NTL has a second option here when the number of factors is small
	Recombine(factors,f,&padicFactors,sigmas,P,p1,e1,coeff_bnd,modularFactors,sigmas2);

//	if (rev) {
//		for (i = 0; i < r; i++) {
//			inplace_rev(factors[i]);
//			if (sign(LeadCoeff(factors[i])) < 0)
//				negate(factors[i], factors[i]);
//		}
//	}

	if (xfac) {
		vec_DUZP_t_push(factors,ex);
	}

	if (x1fac) {
		vec_DUZP_t_push(factors,exm1);
	}

	// that's it!!

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"  [FPSFI] factoring completed. degree sequence:\n  ");
		for (int i=0; i<factors->size; ++i)
			fprintf(stderr,"%ld ",factors->polys[i]->lt);
		fprintf(stderr,"\n");
//		fprintf(stderr,"  [FPSFI] factors->size = %ld\n",factors->size);
//		for (int i=0; i<factors->size; ++i)
//			printPoly_DUZP(factors->polys[i],x);
//		fprintf(stderr,"\n");
	#endif

	mpz_clear(c1);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(P);
	mpz_clear(targetBound);
	freePolynomial_DUZP(f1);
	freePolynomial_DUZP(f);
	freePolynomial_DUZP(ex);
	freePolynomial_DUZP(exm1);
	freeFactors_spX(&modularFactors);
	vec_DUZP_t_free(&padicFactors);
	if (sigmas2 != NULL) {
		for (int i=0; i<nf; ++i) {
			freePolynomial_spX(&sigmas2[i]);
		}
		free(sigmas2);
	}
	free(Pptr);
	factoring_info_t_free(&factoringInfo);
}

/**
 * A wrapper for factor_prim_sqf_inner, which computes
 * the factorization of a primitive and squarefree
 * polynomial.
 *
 * In the NTL approach, this is used for another level
 * of optimization, doing a deflation test and controls
 * behaviour depending on whether their power hack is
 * used. If we decide to add these, they would be
 * implemented here; if not, this routine can be
 * merged with factor_prim_sqf_inner.
 */
void DUZP::Factoring::factor_prim_sqf(vec_DUZP_t* factors, const DUZP_t* ff, long bnd)
{
	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr," [FPS] factor_prim_sqf:\n");
	#endif
//	fprintf(stderr,"factor_prim_sqf input: ff = ");
//	printPoly_DUZP(ff,x);

// input is primitive and square-free, with positive leading
// coefficient

	if (factors == NULL)
		return;

	vec_DUZP_t_clear(factors);

	if (isZero_DUZP(ff)) {
		fprintf(stderr,"BPAS, error: input f to factor_prim_sqf cannot be zero.\n");
		exit(1);
	}

	if (ff->lt <= 0) {
		return;
	}


#if DO_NTL_UNIVAR_FACT
	NTL::ZZX ntlFF;
	DUZP_t* tmpFact = NULL;
	DUZP_long_t pair;
	NTLZZX_set_DUZP_t(ntlFF,ff);
	NTL::ZZ cont;
	NTL::vec_pair_ZZX_long NTLfactors;
	factor(cont, NTLfactors, ntlFF, 0 /*verbose*/, bnd);
	for (int i=1; i<=NTLfactors.length(); ++i) {
		DUZP_t_set_NTLZZX(&tmpFact,NTLfactors(i).a);
		pair.a = tmpFact;
		vec_DUZP_t_push(factors,tmpFact);
		freePolynomial_DUZP(tmpFact);
		tmpFact = NULL;
	}
#else
	// TODO: determine if we want to implement NTL's power hack and deflation check (see NTL SFFactor code)
	factor_prim_sqf_inner(factors,ff,bnd);
//	fprintf(stderr," [FPS] factors->size = %ld\n",factors->size);
	return;
#endif
}

/**
 * Computes the squarefree decomposition of the input
 * polynomial ff using Yun's squarefree factorization
 * algorithm.
 */
void DUZP::Factoring::SquareFreeDecompose(vec_DUZP_long_t* u, const DUZP_t* ff) {

	if (ff->lt <= 0)
		return;

	if (u == NULL)
		return;

	if (u->alloc == 0)
		vec_DUZP_long_t_init2(u,0);
	else
		vec_DUZP_long_t_set(u,0);

	// input is primitive
	DUZP_t* f = deepCopyPolynomial_DUZP(ff);

	DUZP_t* d;
	DUZP_t* v;
	DUZP_t* w;
	DUZP_t* s = makeConstPolynomial_DUZP(1,0);
	DUZP_t* t1;
	DUZP_t* r;
	long i;
	DUZP_long_t pair;

	t1 = differentiate_DUZP(f);
	d = primitiveGCD_DUZP(f,t1);

	if (d->lt == 0) {
		pair.a = f;
		pair.b = 1L;
		vec_DUZP_long_t_push(u,pair);
		freePolynomial_DUZP(f);
		freePolynomial_DUZP(t1);
		freePolynomial_DUZP(d);
		freePolynomial_DUZP(s);
		return;
	}

	dividePolynomials_DUZP(f,d,&v,&r);
	freePolynomial_DUZP(r);
	dividePolynomials_DUZP(t1,d,&w,&r);
	i = 0;

	for (;;) {
		i = i + 1;

		freePolynomial_DUZP(t1);
		t1 = differentiate_DUZP(v);
		freePolynomial_DUZP(s);
		s = subtractPolynomials_DUZP(w,t1);

		if (isZero_DUZP(s)) {
			if (v->lt != 0) {
				pair.a = v;
				pair.b = i;
				vec_DUZP_long_t_push(u,pair);
			}
			freePolynomial_DUZP(f);
			freePolynomial_DUZP(t1);
			freePolynomial_DUZP(d);
			freePolynomial_DUZP(v);
			freePolynomial_DUZP(w);
			freePolynomial_DUZP(r);
			freePolynomial_DUZP(s);
			return;
		}

		freePolynomial_DUZP(d);
		d = GCD_DUZP(v,s); // d = gcd(v,s)
		freePolynomial_DUZP(r);
		freePolynomial_DUZP(t1);
		dividePolynomials_DUZP(v,d,&t1,&r);
		freePolynomial_DUZP(v);
		v = deepCopyPolynomial_DUZP(t1); // v = v/d
		freePolynomial_DUZP(r);
		freePolynomial_DUZP(w);
		dividePolynomials_DUZP(s,d,&w,&r); // w = s/d

		if (d->lt != 0) {
			pair.a = d;
			pair.b = i;
			vec_DUZP_long_t_push(u,pair);
		}
	}

}

/**
 * Exported univariate factorization routine with an
 * optional lifting bound bnd (if bnd is not supplied,
 * the Landau-Mignotte bound will be used).
 *
 * The routine returns the content c along with the
 * vector of factors of ff over Z.
 */
void DUZP::Factoring::factorWithBound(mpz_t* c, vec_DUZP_long_t* factors, const DUZP_t* ff, long bnd) {

	if (factors == NULL)
		return;

	vec_DUZP_long_t_clear(factors);

	if (ff->lt <= 0) {
		mpz_set(*c,ff->coefs[0]);
		return;
	}

	DUZP_t* f = deepCopyPolynomial_DUZP(ff);
	primitivePartAndContent_DUZP_inp(f, *c);
	if (mpz_cmp_ui(f->coefs[f->lt],0) < 0) {
		mpz_neg(*c,*c);
		negatePolynomial_DUZP_inp(f);
	}

#if DO_NTL_UNIVAR_FACT
	NTL::ZZX ntlFF;
	DUZP_t* tmpFact = NULL;
	DUZP_long_t pair;
	NTLZZX_set_DUZP_t(ntlFF,f);
	NTL::ZZ cont;
	NTL::vec_pair_ZZX_long NTLfactors;
	factor(cont, NTLfactors, ntlFF, 0 /*verbose*/, bnd);
	for (int i=1; i<=NTLfactors.length(); ++i) {
		DUZP_t_set_NTLZZX(&tmpFact,NTLfactors(i).a);
		pair.a = tmpFact;
		pair.b = NTLfactors(i).b;
		vec_DUZP_long_t_push(factors, pair);
		freePolynomial_DUZP(tmpFact);
		tmpFact = NULL;
	}
	freePolynomial_DUZP(f);
#else

	// TODO: return here if the input is linear
	// TODO: other irreducibility heuristics/tests to add? Eisenstein's criterion.

	long bnd1 = MaxBits(f) + (NumBits(f->lt+1)+1)/2;
	if (!bnd || bnd > bnd1)
		bnd = bnd1;

	vec_DUZP_long_t sfd;
	vec_DUZP_long_t_init2(&sfd,0);
	DUZP_long_t pair;

	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"[FACTOR] square-free decomposition...");
	#endif
	SquareFreeDecompose(&sfd, f);
	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"done.\n");
	fprintf(stderr,"[FACTOR] square-free decomposition output:\n");
	print_vec_DUZP_long(sfd);
	#endif

	vec_DUZP_t ex;
	vec_DUZP_t_init2(&ex,0);

	for (int i=0; i<sfd.size; ++i) {

		factor_prim_sqf(&ex, sfd.pairs[i].a, bnd);
//		fprintf(stderr,"[FACTOR] ex.>size = %ld\n",ex.size);

		for (int j=0; j<ex.size; ++j) {
			pair.a = ex.polys[j];
			pair.b = sfd.pairs[i].b;
			vec_DUZP_long_t_push(factors, pair);
		}
	}

//	for (int i=0; i<factors->size; ++i) {
//		fprintf(stderr,"e = %ld, f = ",factors->pairs[i].b);
//		printPoly_DUZP(factors->pairs[i].a,x);
//	}
	#if defined(SMZP_FACTORING_DEBUG) && SMZP_FACTORING_DEBUG
		fprintf(stderr,"[FACTOR] factorization complete.\n");
	#endif

	freePolynomial_DUZP(f);
	vec_DUZP_t_free(&ex);
	vec_DUZP_long_t_free(&sfd);
#endif
}
