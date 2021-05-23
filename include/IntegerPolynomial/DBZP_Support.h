/**
 * Supporting methods for dense bivariate polynomials over integers (DBZP)
 *  written in C.
 *
 */

#ifndef _DBZP_SUPPORT_
#define _DBZP_SUPPORT_

#ifdef __cplusplus
    extern "C" {
#endif

#include "DUZP_Support.h"

/**
 * Dense Bivariate Integer Polynomial C struct Z[x > y].
 * Dense array of DUZP_t polynomials over Z[y], coefsY.
 * partial leading term based on X, ltX and coefsY allocation size, allocX.
 */
typedef struct DBZP {
	polysize_t ltX;
	polysize_t allocX;
	DUZP_t** coefsY;
} DBZP_t;


/**
 * Allocate a DBZP with alloc number of coefficients.
 * i.e., alloc = partial maximum degree - 1.
 * returns the newly allocated DBZP.
 */
static inline DBZP_t* makePolynomial_DBZP(int allocX) {
	if (allocX <= 0) {
		allocX = 1;
	}
	DBZP_t* p = (DBZP_t*) malloc(sizeof(DBZP_t));
	p->coefsY = (DUZP_t**) malloc (allocX*sizeof(DUZP_t*));
	for (polysize_t i = 0; i < allocX; i++) {
		p->coefsY[i] = NULL;
	}
	p->allocX = allocX;
	p->ltX = 0;
	return p;
}

/**
 * Allocate a DBZP with allocX number of coefficients and each
 * coefficient of allocation allocY.
 * Like allocPolynomial_DUZP, no coefficient is yet initialized.
 * @return the newly allocated DBZP
 */
static inline DBZP_t* allocPolynomial_DBZP(int allocX, int allocY) {
	if (allocX <= 0) {
		allocX = 1;
	}
	if (allocY <= 0) {
		allocY = 1;
	}

	DBZP_t* p = (DBZP_t*) malloc(sizeof(DBZP_t));
	p->coefsY = (DUZP_t**) malloc(allocX*sizeof(DUZP_t*));
	for (int i = 0; i < allocX; ++i) {
		p->coefsY[i] = allocPolynomial_DUZP(allocY);
	}
	p->allocX = allocX;
	p->ltX = 0;
	return p;
}

/**
 * Make a constant DBZP_t polynomial
 *
 * @param allocX allocation size
 * @param z coef (signed long)
 * @return
 */
static inline DBZP_t* makeConstPolynomial_DBZP(int allocX, signed long z) {
	DBZP_t* ret = makePolynomial_DBZP(allocX);
	ret->coefsY[0] = makeConstPolynomial_DUZP (1, z);
	return ret;
}

/**
 * Make a constant DBZP_t polynomial
 *
 * @param allocX allocation size
 * @param z coef (mpz_t)
 * @return
 */
static inline DBZP_t* makeBigConstPolynomial_DBZP(int allocX, const mpz_t z) {
	DBZP_t* ret = makePolynomial_DBZP(allocX);
	ret->coefsY[0] = makeBigConstPolynomial_DUZP (1, z);
	return ret;
}

/**
 * Free a DBZP.
 * @note
 */
static inline void freePolynomial_DBZP(DBZP_t* p) {
	if (p != NULL) {
		polysize_t lt = p->ltX;
		for (polysize_t i = 0; i <= lt; i++) {
			freePolynomial_DUZP (p->coefsY[i]);
		}
		free(p->coefsY);
		free(p);
	}
}

/**
 * Resize a DBZP polynomial
 * @note this works in-place
 *
 * @param p DBZP_t polynomial
 * @param allocSize new allocation size
 */
static inline void resizePolynomial_DBZP(DBZP_t* p, polysize_t allocSize) {
	if (p == NULL) {
		return;
	}

	if (allocSize == p->allocX) {
		return;
	}

	if (allocSize <= p->ltX) {
		for (polysize_t i = p->ltX; i >= allocSize; --i) {
			freePolynomial_DUZP (p->coefsY[i]);
		}
		p->ltX = allocSize-1;
	}
	p->coefsY = (DUZP_t**) realloc (p->coefsY, sizeof(DUZP_t*)*allocSize);
	for (polysize_t i = p->ltX+1; i < allocSize; i++) {
		p->coefsY[i] = NULL;
	}
	p->allocX = allocSize;
}

/**
 * Evaluate coefsY of a DBZP at a constant c
 *
 * @param p DBZP_t polynomial
 * @param c a evaluating point
 * @return a bivariate polynomial (DBZP_t) over Z[x].
 */
static inline DBZP_t* partialEvaluate_DBZP(const DBZP_t* p, const mpz_t c) {
	if (p == NULL) {
		return NULL;
	}

	polysize_t ltX = p->ltX;
	DBZP_t* px = makePolynomial_DBZP (ltX+1);
	px->ltX = ltX;
	mpz_t vals;
	for (polysize_t i = 0; i <= ltX; i++) {
		mpz_init (vals);
		evaluate_DUZP (p->coefsY[i], c, vals);
		px->coefsY[i] = makeBigConstPolynomial_DUZP (1, vals);
		mpz_clear (vals);
	}

	polysize_t newlt = ltX;
	while (newlt >= 0 && px->coefsY[newlt] == NULL) {
		--newlt;
	}
	if (newlt < 0) {
		freePolynomial_DBZP (px);
		return NULL;
	}
	resizePolynomial_DBZP (px, newlt+1);
	return px;
}

/**
 * Evaluate a DBZP \in Z[x,y] at (x,y) = (c,d)
 *
 * @param p DBZP_t polynomial
 * @param c a X evaluating point
 * @param d a Y evaluating point
 * @param var result of evaluation
 */
static inline void evaluate_DBZP(const DBZP_t* p, const mpz_t c, const mpz_t d, mpz_t var) {
	if (p == NULL) {
		return;
	}

	polysize_t ltX = p->ltX;
	DBZP_t* px = makePolynomial_DBZP (ltX+1);
	px->ltX = ltX;
	mpz_t vals;
	// evaluate at y=d
	for (polysize_t i = 0; i <= ltX; i++) {
		mpz_init (vals);
		evaluate_DUZP (p->coefsY[i], d, vals);
		px->coefsY[i] = makeBigConstPolynomial_DUZP (1, vals);
		mpz_clear (vals);
	}

	polysize_t newlt = ltX;
	while (newlt >= 0 && px->coefsY[newlt] == NULL) {
		--newlt;
	}
	if (newlt < 0) {
		freePolynomial_DBZP (px);
		return;
	}

	// evaluate at x=c
	mpz_t result;
	mpz_init_set (result, px->coefsY[newlt]->coefs[0]);
	for (polysize_t i = newlt-1; i >= 0; i++) {
		mpz_mul (result, result, c);
		if (px->coefsY[i] != NULL) {
			mpz_add (result, result, px->coefsY[i]->coefs[0]);
		}
	}

	mpz_init_set (var, result);
}

/**
 * Deep copy a DBZP.
 * @return the new DBZP.
 */
static inline DBZP_t* deepCopyPolynomial_DBZP(const DBZP_t* p) {
	if (p == NULL) {
		return NULL;
	}

	polysize_t ltX = p->ltX;
	DBZP_t* dp = makePolynomial_DBZP (ltX+1);
	for (polysize_t i = 0; i < ltX; i++) {
		dp->coefsY[i] = deepCopyPolynomial_DUZP (p->coefsY[i]);
	}
	dp->ltX = ltX;
	return dp;
}

/**
 * Deep copy the DBZP a into the DBZP pointed to by bb.
 * If the DBZP pointed to by bb is not NULL, re-use any allocation
 * available to store the copy of a.
 *
 * @param a the DBZP to copy
 * @param[in,out] bb a pointer to a (possibly) pre-allocated DBZP to store the copy.
 */
void deepCopyPolynomial_DBZP_inp(const DBZP_t* a, DBZP_t** bb);

static inline int isEqual_DBZP(const DBZP_t* p, const DBZP_t* q) {
	if (p == NULL) {
		if (q == NULL) {
			return 1;
		} else {
			return 0;
		}
	}
	if (q == NULL) {
		return 0;
	}

	if (p->ltX != q->ltX) {
		return 0;
	}

	for (int i = 0; i <= p->ltX; ++i) {
		if (!isEqual_DUZP(p->coefsY[i], q->coefsY[i])) {
			return 0;
		}
	}

	return 1;
}

/**
 * (distributed) Bivariate Subresultant Chain Data-Type over Z
 * @note each subresultant-chain is stored in one C-array
 * @note using distributed-data-type led to a better performance for modularSubres over Z[x>y].
 */
typedef struct {
    polysize_t n; // Number of polynomials in a subresultant-chain, at least 2
    polysize_t** deg; // Partial degree of each polynomial (bi-var: deg[2])
    polysize_t* size; // size of each polynomial array
    mpz_t** coefs; // coefficients of polynomials
} biSubresZ_t;

static inline biSubresZ_t* makeBiSubresultant_DBZP (polysize_t n) {
    if (n < 0) {
        return NULL;
    }
    biSubresZ_t* res = (biSubresZ_t*) malloc (sizeof(biSubresZ_t));
    res->n = n;
    if (n > 0) {
		res->deg = (polysize_t**) malloc (n*sizeof(polysize_t*));
		res->coefs = (mpz_t**) malloc (n*sizeof(mpz_t*));
        res->size = (polysize_t*) calloc (n, sizeof (polysize_t));
        for (polysize_t i = 0; i < n; i++) {
            res->deg[i] = (polysize_t*) calloc (2, sizeof(polysize_t));
        }
    } else {
        res->n = 0;
        res->deg = NULL;
        res->size = NULL;
        res->coefs = NULL;
    }
	return res;
}

static inline void freeBiSubresultant_DBZP (biSubresZ_t* s) {
    if (s == NULL) {
        return;
    }
    if (s->n) {
        for (polysize_t i = 0; i < s->n; i++) {
            free (s->deg[i]);
            if (s->size[i]) {
				for (polysize_t j = 0; j < s->size[i]; j++) {
					mpz_clear(s->coefs[i][j]);
				}
				free (s->coefs[i]);
			}
        }
        free (s->deg);
        free (s->size);
		free (s->coefs);
    }
    free (s);
}

static inline biSubresZ_t* deepCopyBiSubresultant_DBZP (biSubresZ_t* a) {
	if (a == NULL) {
		return NULL;
	}

	biSubresZ_t* ac = makeBiSubresultant_DBZP (a->n);
	if (!a->n) {
		return ac;
	}
	for (polysize_t i = 0; i < a->n; i++) {
		ac->deg[i][0] = a->deg[i][0];
		ac->deg[i][1] = a->deg[i][1];
		ac->size[i]   = a->size[i];
		if (a->size[i]) {
			for (polysize_t j = 0; j < a->size[i]; j++) {
				mpz_set (ac->coefs[i][j], a->coefs[i][j]);
			}
		}
	}
	return ac;
}

static inline int isBiSubresultantEquals_DBZP (biSubresZ_t* a, biSubresZ_t* b) {
	if (a == NULL && b == NULL) {
		return 1;
	} else if (a == NULL) {
		return 0;
	} else if (b == NULL) {
		return 0;
	}
	if (a->n != b->n) {
		// fprintf (stderr, "[DBZP] a->n != b->n\n");
		return 0;
	}
	for (polysize_t i = 0; i < a->n; i++) {
		if (a->deg[i][0] != b->deg[i][0] || a->deg[i][1] != b->deg[i][1]) {
			// fprintf (stderr, "[DBZP] a->deg[i][0] != b->deg[i][0] || a->deg[i][1] != b->deg[i][1]\n");
			return 0;
		} else if (a->size[i] != b->size[i]) {
			// fprintf (stderr, "[DBZP] a->size[i] != b->size[i]\n");
			return 0;
		}
		for (polysize_t j = 0; j < a->size[i]; j++) {
			if (mpz_cmp (a->coefs[i][j], b->coefs[i][j]) ) {
				// fprintf (stderr, "[DBZP] a->coefs[i][j] != b->coefs[i][j]\n");
				return 0;
			}
		}
	}
	return 1;
}

static inline AltArrZ_t* convertToAltArrZ_DBZP (const mpz_t* a, polysize_t* apd)
{
	// // apd[0]: max-deg w.r.t y
	// // apd[1]; max-deg w.r.t x
	// // poly is sorted w.r.t x > y

	int size = 0;
	int n = apd[0] + 1;
	for (polysize_t i = 0; i <= apd[1]; i++) {
		for (polysize_t j = 0; j <= apd[0]; j++) {
			size += (mpz_sgn(a[n*i+j]) != 0);
		}
	}
	if (size == 0) {
		return NULL;
	}

	AltArrZ_t* aa = makePolynomial_AAZ(size, 2);

	int curSize = 0;
	//i goes to 0, inlcusive.
	for (degrees_t i = apd[1]+1; i -- > 0; ) {
		for (degrees_t j = apd[0]+1; j -- > 0; ) {
			if (mpz_sgn(a[n*i+j])) {
				mpz_init_set(aa->elems[curSize].coef, a[n*i + j]);
				aa->elems[curSize].degs = ((i) << EXP_OFFSET_1_V2) | (j);
				++curSize;
			}
		}
	}
	aa->size = curSize;
	return aa;

	// AltArrZ_t* aZZ = NULL;
	// degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t));
	// // int n = apd[0] + 1;
	// for (polysize_t i = 0; i <= apd[1]; i++) {
	// 	for (polysize_t j = 0; j <= apd[0]; j++) {
	// 		if (mpz_sgn(a[n*i+j])) {
	// 			// fprintf (stderr, "deg[0] = %ld , deg[1] = %ld\n", j, i);
	// 			// gmp_fprintf (stderr, " a[n*i+j] = %Zd\n", a[n*i+j]);
	// 			// a true monomial with [j][i]
	// 			degs[0] = i; degs[1] = j;
	// 			if (aZZ == NULL) {
	// 				aZZ = makeConstPolynomial_AAZ(1, 2, a[n*i+j]);
	// 				setDegrees_AAZ_inp(aZZ, 0, degs, 2);
	// 			} else {
	// 				setCoefficient_AAZ (aZZ, degs, 2, a[n*i+j]);
	// 			}
	// 		}
	// 	}
	// }
	// free (degs);
	// if (aZZ != NULL) {
	// 	mergeSortPolynomial_AAZ (aZZ);
	// }
	// return aZZ;
}

// TODO remake this like the above convertToAltArrZ_DBZP
static inline AltArrZ_t* convertToAltArrZ_DBSP (const elem_t* a, polysize_t* apd, Prime_ptr* Pptr)
{
	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	AltArrZ_t* aZZ = NULL;

	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t));
	mpz_t tmp;
	mpz_init (tmp);
	int n = apd[0] + 1;
	for (polysize_t i = 0; i <= apd[1]; i++) {
		for (polysize_t j = 0; j <= apd[0]; j++) {
			if (a[n*i+j]) {
				degs[0] = i; degs[1] = j;
				// gmp_fprintf (stderr, " a[n*i+j] = %Zd\n", a[n*i+j]);
				mpz_set_si (tmp, smallprimefield_convert_out (a[n*i+j], Pptr));
				if (aZZ == NULL) {
					aZZ = makeConstPolynomial_AAZ(1, 2, tmp);
					setDegrees_AAZ_inp(aZZ, 0, degs, 2);
				} else {
					setCoefficient_AAZ (aZZ, degs, 2, tmp);
				}
			}
		}
	}
	mpz_clear (tmp);
	free (degs);
	if (aZZ != NULL) {
		mergeSortPolynomial_AAZ (aZZ);
	}
	return aZZ;
}

// don't need to init mpz_t* a
static inline mpz_t* convertFromAltArrZ_DBZP (AltArrZ_t* aZ, polysize_t* apd)
{
	if (isZero_AAZ (aZ) || aZ->nvar != 2) {
		return NULL;
	}

	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	if (apd == NULL) {
		fprintf (stderr, "DUZP Error: apd must be initialized\n");
		exit (1);
	}
	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t));
	partialDegrees_AAZ (aZ, degs);
	apd[0] = degs[1];
	apd[1] = degs[0];

	int allocSz = (degs[1]+1) * (degs[0]+1);
	// fprintf(stderr, "allocSz = %d\n", allocSz);
	mpz_t* a = (mpz_t*) malloc (allocSz *sizeof(mpz_t));
	for (int k = 0; k < allocSz; k++) {
		mpz_init_set_d (a[k], 0);
	}
	int n = apd[1]+1;
	for (int k = 0; k < aZ->size; k++) {
		partialDegreesTerm_AAZ (aZ, k, degs);
		// fprintf (stderr, "a[%d] set!\n", n*degs[1]+degs[0]);
		mpz_set (a[n*degs[1]+degs[0]], aZ->elems[k].coef);
		degs[0] = degs[1] = 0;
	}
	// for (int i = 0; i < allocSz; i++) {
	// 	gmp_fprintf (stderr, "a[%d] = %Zd\n", i, a[i]);
	// }
	free (degs);
	return a;
}

static inline elem_t* convertFromAltArrZ_DBSP (AltArrZ_t* aZ, polysize_t* apd, Prime_ptr* Pptr)
{
	if (isZero_AAZ (aZ) || aZ->nvar != 2) {
		return NULL;
	}

	// apd[0]: max-deg w.r.t y
	// apd[1]; max-deg w.r.t x
	// poly is sorted w.r.t x > y
	if (apd == NULL) {
		fprintf (stderr, "DUZP Error: apd must be initialized\n");
		exit (1);
	}
	degree_t* degs = (degree_t*) calloc (2, sizeof(degree_t));
	partialDegrees_AAZ (aZ, degs);
	apd[0] = degs[1];
	apd[1] = degs[0];

	int allocSz = (degs[1]+1) * (degs[0]+1);
	// fprintf(stderr, "allocSz = %d\n", allocSz);
	elem_t* a = (elem_t*) calloc (allocSz, sizeof(elem_t));
	mpz_t tmp;
	mpz_init (tmp);
	int n = apd[1]+1;
	for (int k = 0; k < aZ->size; k++) {
		partialDegreesTerm_AAZ (aZ, k, degs);
		// fprintf (stderr, "a[%d] set!\n", n*degs[1]+degs[0]);
		mpz_mod_ui (tmp, aZ->elems[k].coef, (unsigned long) Pptr->prime);
		a[n*degs[1]+degs[0]] = mpz_get_si (aZ->elems[k].coef);
		a[n*degs[1]+degs[0]] = smallprimefield_convert_in (a[n*degs[1]+degs[0]], Pptr);
		degs[0] = degs[1] = 0;
	}
	// for (int i = 0; i < allocSz; i++) {
	// 	fprintf (stderr, "a[%d] = %Zd\n", i, a[i]);
	// }
	free (degs);
	return a;
}


static inline void printDBZP_AAZ (mpz_t* a, polysize_t* apd)
{
  if (apd[0]==0 || apd[1]==0) {
    fprintf (stderr, "0;\n");
  }

  polysize_t ad = apd[1] + 1;
  for (int i = 0; i <= apd[0]; i++) {
    for (int j = 0; j <= apd[1]; j++) {
      gmp_fprintf (stderr, "%Zd\t",a[ad*i+j]);
    }
    fprintf (stderr, "\n");
  }
}

static inline void printDBSP_AAZ (elem_t* a, polysize_t* apd, Prime_ptr* Pptr)
{
  if (apd[0]==0 || apd[1]==0) {
    fprintf (stderr, "0;\n");
  }

  polysize_t ad = apd[1] + 1;
  for (int i = 0; i <= apd[0]; i++) {
    for (int j = 0; j <= apd[1]; j++) {
      fprintf (stderr, "%lld\t", smallprimefield_convert_out (a[ad*i+j], Pptr));
    }
    fprintf (stderr, "\n");
  }
}

/**
 * Modular Bivariate Subresultant Chain working with a prime-set
 * @note see fourier_primes_u64_table
 *
 * @param g a distributed 2D polynomial over Z
 * @param gpd partial degrees gpd[0], gpd[1]
 * @param f a distributed 2D polynomial over Z s.t. deg (g, x) > deg (f, x)
 * @param fpd partial degrees fpd[0], fpd[1]
 * @param m returned-middle-computation in CRT
 * @param startIdx starting from prime64_ptr[startIdx] in CRT
 * @param primeSz size of 'primes' list: {startIdx+0, startIdx+1, ..., startIdx+primeSz-1}
 * @return biSubresZ_t mod \prod_{i} primes[i] with CRT properties
 */
biSubresZ_t* modularSetBiSubresultant_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd, mpz_t m, int startIdx, int primeSz);

/**
 * Modular Bivariate Subresultant Chain over Z
 * @note see fourier_primes_u64_table
 *
 * @param g a distributed 2D polynomial over Z
 * @param gpd partial degrees gpd[0], gpd[1]
 * @param f a distributed 2D polynomial over Z s.t. deg (g, x) > deg (f, x)
 * @param fpd partial degrees fpd[0], fpd[1]
 * @return biSubresZ_t mod \prod_{i} primes[i] with CRT properties
 */
biSubresZ_t* modularBiSubresultant_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd);



#ifdef __cplusplus
}
#endif

#endif 