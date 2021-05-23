/**
 * Supporting methods for dense bivariate polynomials over small-prime fields (DBSP)
 *  written in C.
 */

#ifndef _DBSP_SUPPORT_H_
#define _DBSP_SUPPORT_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "DUSP_Support.h"

/******************************************
 * Bivariate Polynomials over Small Primes
 *****************************************/
/**
 * Dense Bivariate Integer Polynomial C struct Z_p[x > y].
 * Dense array of duspoly_t polynomials over Z_p[y], coefsY.
 * partial leading term based on X, ltX and coefsY allocation size, allocX.
 */
typedef struct {
    polysize_t ltX;
    polysize_t allocX;
    duspoly_t** coefsY;
} dbspoly_t;

/**
 * Allocate a dbspoly with alloc number of coefficients.
 * i.e., alloc = partial maximum degree - 1.
 * returns the newly allocated dbpoly.
 */
static inline dbspoly_t* makePolynomial_spXY (int allocX) {
	if (allocX <= 0) {
		allocX = 1;
	}
	dbspoly_t* p = (dbspoly_t*) malloc(sizeof(dbspoly_t));
	p->coefsY = (duspoly_t**) malloc (allocX*sizeof(duspoly_t*));
	for (polysize_t i = 0; i < allocX; i++) {
		p->coefsY[i] = NULL;
	}
	p->allocX = allocX;
	p->ltX = 0;
	return p;
}

static inline dbspoly_t* allocPolynomial_spXY(polysize_t allocX, polysize_t allocY) {
    if (allocX <= 0) {
        allocX = 1;
    }
    if (allocY <= 0) {
        allocY = 1;
    }
    dbspoly_t* ret = (dbspoly_t*) malloc(sizeof(dbspoly_t));
    ret->coefsY = (duspoly_t**) malloc(sizeof(duspoly_t*)*(allocX));
    for (polysize_t i = 0; i < allocX; ++i) {
        ret->coefsY[i] = makePolynomial_spX(allocY);
    }
    ret->allocX = allocX;
    ret->ltX = 0;
    return ret;
}

/**
 * Make a constant dbspoly_t polynomial
 *
 * @param allocX allocation size
 * @param z NOT-in-form coef (elem_t)
 * @return
 */
static inline dbspoly_t* constPolynomialInForm_spXY (int allocX, elem_t z, const Prime_ptr* Pptr) {
	dbspoly_t* ret = makePolynomial_spXY (allocX);
	ret->coefsY[0] = constPolynomialInForm_spX (z, Pptr);
	return ret;
}


/**
 * Make a constant dbspoly_t polynomial
 *
 * @param allocX allocation size
 * @param z InForm coef (elem_t)
 * @return
 */
static inline dbspoly_t* constPolynomial_spXY (int allocX, elem_t z, const Prime_ptr* Pptr) {
	dbspoly_t* ret = makePolynomial_spXY (allocX);
	ret->coefsY[0] = constPolynomial_spX (z, Pptr);
	return ret;
}

/**
 * Free a dbspoly_t.
 * @note
 */
static inline void freePolynomial_spXY (dbspoly_t* p) {
	if (p != NULL) {
		polysize_t lt = p->ltX;
		for (polysize_t i = 0; i <= lt; i++) {
			freePolynomial_spX (&(p->coefsY[i]));
		}
		free(p->coefsY);
		free(p);
	}
}

/**
 * Resize a dbspoly polynomial
 * @note this works in-place
 *
 * @param p dbspoly polynomial
 * @param allocSize new allocation size
 */
static inline void resizePolynomial_spXY (dbspoly_t* p, polysize_t allocSize) {
	if (p == NULL) {
		return;
	}

	if (allocSize == p->allocX) {
		return;
	}

	if (allocSize <= p->ltX) {
		for (polysize_t i = p->ltX; i >= allocSize; --i) {
			freePolynomial_spX (&(p->coefsY[i]));
		}
		p->ltX = allocSize-1;
	}
	p->coefsY = (duspoly_t**) realloc (p->coefsY, sizeof(duspoly_t*)*allocSize);
	for (polysize_t i = p->ltX+1; i < allocSize; i++) {
		p->coefsY[i] = NULL;
	}
	p->allocX = allocSize;
}

static inline int isEqual_spXY(const dbspoly_t* p, const dbspoly_t* q) {
    if (p == NULL) {
        return (q == NULL);
    }
    if (q == NULL) {
        return 0;
    }

    if (p->ltX != q->ltX) {
        return 0;
    }

    for (int i = 0; i <= p->ltX; ++i) {
        if (!isEqual_spX(p->coefsY[i], q->coefsY[i])) {
            return 0;
        }
    }

    return 1;
}

/**
 * Convert a bivariate SMZP polynomial to a dbspoly_t over
 * the finite field described by Pptr.
 *
 * @note it is assumed that q is pre-allocated with enough space to hold
 *       all coefficients of p.
 *
 * @param p the bivariate SMZP to convert
 * @param q[out] the pre-allocated dbsp to be filled with p's data
 * @param Pptr the small prime pointer
 */
void convertFromAltArrZ_spXY_inp(const AltArrZ_t* p, dbspoly_t** q, const Prime_ptr* Pptr);

/**
 * Compute the primitive part and content of p, with resepct to the main variable.
 * For Z_p[x >y] its content is a polynomial in Z_[y], the gcd of all of the
 * univariate polynomial coefficients.
 *
 * This function is in place in the sense that p is itself transformed into its
 * primitive part. The content is returned as a new polynomial pointed to by cont,
 * if cont is non-NULL. If cont already points to a non-NULL duspoly_t*, that polynomial will be free'd.
 *
 * @param p the polynomial to convert to its primitive part
 * @param[out] cont the content of the input polynomial p
 * @param Pptr the small prime pointer
 */
static inline void mainPrimitivePart_spXY_inp(dbspoly_t* p, duspoly_t** cont, const Prime_ptr* Pptr) {
    if (p == NULL) {
        if (cont == NULL) {
            *cont = NULL;
        }
        return;
    }

    if (p->ltX == 0) {
        if (cont != NULL) {
            *cont = p->coefsY[0];
        }
        p->coefsY[0] = constPolynomial_spX (1, Pptr);
        return;
    }

    duspoly_t* g = deepCopyPolynomial_spX(p->coefsY[0]);
    duspoly_t* nextG = NULL;
    for (int i = 1; i <= p->ltX && !isOne_spX(g); ++i) {
        if (isZero_spX(p->coefsY[i])) {
            continue;
        }
        GCDInForm_spX(g, p->coefsY[i], &nextG, Pptr);
        freePolynomial_spX(&g);
        g = nextG;
        nextG = NULL;
    }

    if (!isOne_spX(g)) {
        for (int i = 0; i <= p->ltX; ++i) {
            plainExactDivPolynomialsInForm_spX (p->coefsY[i], g, &nextG, Pptr);
            freePolynomial_spX(p->coefsY + i);
            p->coefsY[i] = nextG;
        }
    }

    if (cont != NULL) {
        if (*cont != NULL) {
            freePolynomial_spX(cont);
        }
        *cont = g;
    } else {
        freePolynomial_spX(&g);
    }


}

/**
 * Evaluate the polynomial p at y=pt, yielding a univariate polynomial in x.
 * The result is stored in the polynomial pointed to by eval. If that polynomial is NULL,
 * a new polynomial is created and stored in eval.
 *
 * @param p the polynomial to evaluate
 * @param eval a pointer to a (possible NULL) polynomial to store the evaluation.
 * @param pt a small prime field element in montgomery form
 * @param Pptr the small prime pointer
 *
 */
static inline void evalCoefficientVariableInForm_spXY(const dbspoly_t* p, duspoly_t** eval, elem_t pt, const Prime_ptr* Pptr) {
    if (eval == NULL) {
        return;
    }

    if (p == NULL) {
        if (*eval != NULL) {
            freePolynomial_spX(eval);
            *eval = NULL;
        }
        return;
    }

    duspoly_t* e = *eval;
    if (e == NULL) {
        e = makePolynomial_spX(p->ltX + 1);
        *eval = e;
    }
    if (e->alloc <= p->ltX) {
        reallocatePolynomial_spX(&e, p->ltX + 1);
    }

    e->lt = 0;
    for (int i = 0; i <= p->ltX; ++i) {
        e->elems[i] = evaluateInForm_spX(p->coefsY[i], pt, Pptr);
        e->lt = (e->elems[i] != 0) ? i : e->lt;
    }
}


/**
 * Interpolate a bivariate polynomial Z_p[x>y] from several evaluation images.
 * images[i] is the evaluation of the polynomial to interpolate at the point y=t[i].
 *
 * Both the points and the images should be in montgomery form.
 *
 * If ret points to a non-NULL dbspoly_t*, that polynomial's allocation is re-used
 * to store the interpolated polynomial.
 *
 * @param t the evaluation points
 * @param images the univariate images from which to interpolate coefficients in y
 * @param n the number of points and images
 * @param[out] ret a pointer to the interpolated polynomial.
 * @param Pptr the small prime pointer.
 *
 */
void interpPolyUnivarImagesInForm_spXY(const elem_t* t, duspoly_t const*const* images, int n, dbspoly_t** ret, const Prime_ptr* Pptr);

/**
 * (distributed) Bivariate Subresultant Chain Data-Type for small prime
 * @note each subresultant-chain is stored in one C-array
 * @note using distributed-data-type led to a better performance for modularSubres over Z[x>y].
 */
typedef struct {
    polysize_t n; // Number of polynomials in a subresultant-chain, at least 2
    polysize_t** deg; // Partial degree of each polynomial (bi-var: deg[2])
    polysize_t* size; // size of each polynomial array
    elem_t** coefs; // coefficients of polynomials
} biSubresPr_t;

static inline biSubresPr_t* makeBiSubresultantInForm_spX (polysize_t n) {
    if (n < 0) {
        return NULL;
    }
    biSubresPr_t* res = (biSubresPr_t*) malloc (sizeof(biSubresPr_t));
    res->n = n;
    if (n > 0) {
        res->deg = (polysize_t**) malloc (n*sizeof(polysize_t*));
        res->coefs = (elem_t**) malloc (n*sizeof(elem_t*));
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

static inline void freeBiSubresultantInForm_spX (biSubresPr_t* s) {
    if (s == NULL) {
        return;
    }
    if (s->n) {
        for (polysize_t i = 0; i < s->n; i++) {
            free (s->deg[i]);
            if (s->size[i])
                free (s->coefs[i]);
        }
        free (s->deg);
        free (s->size);
        free (s->coefs);
    }
    free (s);
}

/**
 * Compute Partial Subresultant Chains w.r.t y in Z[x>y]
 * @note This makes use of lagrange algorithm for interpolation
 * @note There are independent for loops (cilk_for hopes!)
 *
 * @param pSubres partial subresultant chains s.t. every poly in this sequence obtained after specializing x in Z[x>y].
 * @param XYsubres returned bivariate subresultant chains over Z[x>y] modulus pr
 * @param t evaluation points
 * @param n size of the sequence, that is, the degree+1 of the targeted subresultant
 * @param alpha block size
 * @param Pptr small prime pointer
 */
// void _interpolatePartialSubResInForm_spXY (duspolysA_t** pSubres, biSubresPr_t* XYsubres, elem_t* t, polysize_t n, polysize_t alpha, const Prime_ptr* Pptr);

/**
 * Bivariate Subresultant Chain
 *
 * @param g a bivariate polynomial (dbspoly_t)
 * @param gpd 2Dim partial degrees of polynomial g
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f
 * @param subres returned subresultant chain modular pr (BiSybresPr_t)
 * @param Pptr small prime pointer
 */
void biSylvSubResultantInForm_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, biSubresPr_t** subres, const Prime_ptr* Pptr);

/**
 * Speculative Bivariate Subresultant Chain  Z[y<x] but polys sorted w.r.t y,x
 *
 * @param g a bivariate polynomial (elem_t*) Z[y<x]
 * @param gpd 2Dim partial degrees of polynomial g s.t. gpd[0]:deg_y and gpd[1]:deg_x
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f s.t. fpd[0]:deg_y and fpd[1]:deg_x
 * @param kk computes k-th subresultant chain (set it to 0, computes resultant and the previous subresultant chain)
 * @param ksubres returned subresultant chain modular pr (biSybresPr_t)
 * @param Pptr small prime pointer
 */
void hgcdBiSylvSubResultantInForm_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t kk, biSubresPr_t** ksubres, const Prime_ptr* Pptr);

////////////////////////////////
// FFT-based Subresultant Chains
////////////////////////////////

/**
 * Compute N = 2^k > B_y for a k \in \N
 * @note N satifies the DFT algorithm in GFPF_Support.h
 */
static inline polysize_t _compute_N_gt_By (polysize_t* gpd, polysize_t* fpd, int* K_out, int* e_out)
{
    polysize_t B = gpd[0] * fpd[1] + gpd[1] * fpd[0] + 1;
    int K = 2;
    int e = 1;
    while (B > 0) {
        B = B/K;
        e++;
    }
    *K_out = K;
    *e_out = e;
    return pow(K, e); // max_deg + 1 <= N
}

/**
 * convert polynomial elem_t* g with partial degrees deg(g, y) = pdg[0], deg(g, x) = pdg[1] to
 * a polynomial usfixn64* g_new with partial degrees deg(g_new, y) = N, deg(g_new, x) = xpd
 * TODO: It's not an efficient implementation
 */
static inline usfixn64* convert_bivarElemPoly_to_usfixn64Poly (elem_t* g, polysize_t* pdg, polysize_t n)
{
    if (pdg[0] >= n) {
        fprintf(stderr, "Error: expected n be greater than the deg(g, y) in convert_bivarElemPoly_to_usfixn64Poly \n");
        exit(EXIT_FAILURE);
    }
    usfixn64* g_new = (usfixn64*) calloc((pdg[1]+1)*n, sizeof(usfixn64));
    for(polysize_t i = 0; i <= pdg[0]; i++) { // y
        for(polysize_t j = 0; j <= pdg[1]; j++) { // x
            g_new[j * (n) + i] = g[i * (pdg[1]+1) + j];
        }
    }
    return g_new;
}


static inline elem_t* convert_bivarElemPoly_to_elemPoly (elem_t* g, polysize_t* pdg, polysize_t n)
{
    if (pdg[0] >= n) {
        fprintf(stderr, "Error: expected n be greater than the deg(g, y) in convert_bivarElemPoly_to_usfixn64Poly \n");
        exit(EXIT_FAILURE);
    }
    elem_t* g_new = (elem_t*) calloc((pdg[1]+1)*n, sizeof(elem_t));
    for(polysize_t i = 0; i <= pdg[0]; i++) { // y
        for(polysize_t j = 0; j <= pdg[1]; j++) { // x
            g_new[j * (n) + i] = g[i * (pdg[1]+1) + j];
        }
    }
    return g_new;
}

/**
 * Matteo's Transpose (cache-friendly 2D Transpose)
 * out-of-place transpose g[i0..i1-1][j0..j1-1] into f[j0..j1-1][i0..i1-1]
 * g(n*m): gd = m and f(m*n) : fd = n
 * @note integrated in transpose2D_spX
 */
static inline void MatteoTranspose2D_spX (usfixn64 * g, polysize_t gd,
                                    usfixn64* f, polysize_t fd,
                                    polysize_t i0, polysize_t i1,
                                    polysize_t j0, polysize_t j1)
{
    const int THRESHOLD = 16;
    polysize_t di, dj, im, jm;
    polysize_t i, j;

    tail:
    di = i1 - i0; dj = j1 - j0;
    if (di >= dj && di > THRESHOLD) {
        im = (i0 + i1)/2;
        MatteoTranspose2D_spX(g, gd, f, fd,
                            i0, im, j0, j1);
        i0 = im;
        goto tail;
    } else if (dj > THRESHOLD) {
        jm = (j0 + j1)/2;
        MatteoTranspose2D_spX(g, gd, f, fd,
                            i0, i1, j0, jm);
        j0 = jm;
        goto tail;
    } else {
        for (i = i0; i < i1; i++) {
            for (j = j0; j < j1; j++) {
                f[j * fd + i] = g[i * gd + j];
            }
        }
    }
}

static inline void MatteoTranspose2D_elem_spX (elem_t* g, polysize_t gd,
                                    elem_t* f, polysize_t fd,
                                    polysize_t i0, polysize_t i1,
                                    polysize_t j0, polysize_t j1)
{
    const int THRESHOLD = 16;
    polysize_t di, dj, im, jm;
    polysize_t i, j;

    tail:
    di = i1 - i0; dj = j1 - j0;
    if (di >= dj && di > THRESHOLD) {
        im = (i0 + i1)/2;
        MatteoTranspose2D_elem_spX(g, gd, f, fd,
                            i0, im, j0, j1);
        i0 = im;
        goto tail;
    } else if (dj > THRESHOLD) {
        jm = (j0 + j1)/2;
        MatteoTranspose2D_elem_spX(g, gd, f, fd,
                            i0, i1, j0, jm);
        j0 = jm;
        goto tail;
    } else {
        for (i = i0; i < i1; i++) {
            for (j = j0; j < j1; j++) {
                f[j * fd + i] = g[i * gd + j];
            }
        }
    }
}

/**
 * An Efficient Transpose2D n * m => m * n
 */
static inline usfixn64* transpose2D_spX (usfixn64* g, polysize_t n, polysize_t m)
{
    usfixn64* g_out = (usfixn64*) calloc(n*m, sizeof(usfixn64));

    // matteo 2D transpose
    MatteoTranspose2D_spX (g, m, g_out, n, 0, n, 0, m);

    // // naive 2D transpose
    // for(polysize_t i = 0; i < n; i++) {
    //     for(polysize_t j = 0; j < m; j++) {
    //         g_out[j * n + i] = g[i * m + j]; // m * n  = n * m
    //     }
    // }

    free(g);
    return g_out;
}

static inline elem_t* transpose2D_elem_spX (elem_t* g, polysize_t n, polysize_t m)
{
    elem_t* g_out = (elem_t*) calloc(n*m, sizeof(elem_t));

    // matteo 2D transpose
    MatteoTranspose2D_elem_spX (g, m, g_out, n, 0, n, 0, m);

    // // naive 2D transpose
    // for(polysize_t i = 0; i < n; i++) {
    //     for(polysize_t j = 0; j < m; j++) {
    //         g_out[j * n + i] = g[i * m + j]; // m * n  = n * m
    //     }
    // }

    free(g);
    return g_out;
}

/**
 * Row-wise DFT to apply on Bivariate Polynomial Z[x > y] (in-place)
 * @param g a bivariate polynomial Z[x > y] of size (xpd+1)*n
 * @param xpd degree of g w.r.t x
 * @param n = K^e computed N by _compute_N_gt_By
 * @param K ...
 * @param e ...
 * @param w a good primitive root of unity of order m in  Z_p
 * @param Pptr small prime pointer
 */
static inline usfixn64* bivarRowWiseDFT_spX (usfixn64* g, polysize_t gpd_x, polysize_t n,
                        int K, int e, usfixn64 w, montgomery_triple P)
{
    for(polysize_t i = 0; i <= gpd_x; ++i) {
        DFT_general_small_elements (&(g[i*n]), K, e, &P, w);
    }
    return g;
}

static inline elem_t* bivarRowWiseDFT4_spX (elem_t* g, polysize_t gpd_x, polysize_t n, polysize_t r, 
                                            elem_t* Omega, const Prime_ptr* Pptr)
{
    elem_t *SB1;
    for (polysize_t i = 0; i <= gpd_x; ++i) {
        SB1 = (elem_t*) calloc (n, sizeof(elem_t)); // TODO: double-check
        DFT_eff(n, r, &(g[i*n]), Omega, Pptr, (int) FFT_H_SIZE, SB1);
        free(SB1);
    }
    return g; 
}    

/**
 * Row-wise DFT^{-1} to apply on Bivariate Polynomial Z[x > y] (in-place)
 * @note Assume w_inv and n_inv are already convertIn_GLOBAL_ptr
 */
static inline usfixn64* bivarRowWiseDFTInv_spX (usfixn64* g, polysize_t gpd_x, polysize_t n,
                        int K, int e, usfixn64 w_inv, usfixn64 n_inv, montgomery_triple P)
{
    // polysize_t sz = (gpd_x+1)*n;
    //x=DFT_n^{-1}(x)
    for(polysize_t i = 0; i <= gpd_x; i++) {
        DFT_general_small_elements (&(g[i*n]), K, e, &P, w_inv);
        for(polysize_t j = 0; j < n; j++) {
            g[i*n+j] = mult_ab_mod_p (g[i*n+j], n_inv, &P);
        }
    }
    // for(polysize_t j = 0; j < sz; j++) {
    //     g[j] = mult_ab_mod_p (g[j], n_inv, &P); }
    return g;
}

static inline elem_t* bivarRowWiseDFT4Inv_spX (elem_t* g, polysize_t gpd_x, polysize_t n,
                        polysize_t r, elem_t* OmegaInv, elem_t n_inv, const Prime_ptr* Pptr)
{
    //x=DFT_n^{-1}(x)
    elem_t *SB1;
    for(polysize_t i = 0; i <= gpd_x; i++) {
        SB1 = (elem_t*) calloc (n, sizeof(elem_t)); // TODO: double-check?
        InvDFT_eff (n, r, &(g[i*n]), OmegaInv, Pptr, (int) FFT_H_SIZE, SB1, n_inv);
        free(SB1);
    }
    return g;
}


// not functional anymore...
// typedef struct {
//     usfixn64** polys;
//     polysize_t chain_size;
//     polysize_t* alloc_sizes;
// } bivarSubresFFT_t;

/**
 * compute DFT-inverse of generated subresultant chains by
 * _modularSubres_from_usfixn64 and _compute_DFTInv_modularSubres_inputs,
 * @param SS an array of _modularSubres_from_usfixn64 arrays
 * @param subres_size the chain size from _modularSubres_from_usfixn64
 * @param max_chain_alloc_sizes the maximum allocation size for every member of the chain from _modularSubres_from_usfixn64
 * @param n the K^e computed N by _compute_N_gt_By
 * @param K
 * @param e
 * @param w a good primitive root of unity of order m in  Z_p
 * @param Pptr small prime pointer
 *
 * @note subres_size and max_chain_alloc_sizes must be consistent in entire SS
 */
usfixn64** _compute_DFTInv_modularSubres (usfixn64** SS, polysize_t n, int K, int e, usfixn64 w_inv, usfixn64 n_inv,
                                        polysize_t subres_size, polysize_t alpha, polysize_t* max_chain_alloc_sizes,
                                        const Prime_ptr* Pptr);

/**
 * FFT-based Bivariate Subresultant Chain  Z[y<x] but polys sorted w.r.t y,x
 *
 * @param g a bivariate polynomial (elem_t*) Z[y<x]
 * @param gpd 2Dim partial degrees of polynomial g s.t. gpd[0]:deg_y and gpd[1]:deg_x
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f s.t. fpd[0]:deg_y and fpd[1]:deg_x
 * @param k_sub to return k-th and k+1-th subresultants in the chain based on Half-GCD
 * @param n a power of two s.t. n > B_y = gpd[1]*fpd[0] + gpd[0]*fpd[1] + 1
 * @param K
 * @param e
 * @param w a good primitive root of unity of order m in  Z_p
 * @param w_inv inverse of w mod p
 * @param n_inv inverse of n mod p
 * @param subres returned subresultant chain modular pr
 * @param Pptr small prime pointer
 */
int bivarModularSubResInForm_withFFT_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr);
int bivarModularSubResInForm_withFFT4_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd,
                                           polysize_t n, polysize_t r, elem_t* Omega, elem_t* OmegaInv, elem_t n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr);


int bivarModularHGCDSubResInForm_withFFT_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr);
int bivarModularHGCDSubResInForm_withFFT4_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, polysize_t r, elem_t* Omega, elem_t* OmegaInv, elem_t n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr);

typedef struct {
    specA_spX_t ** A;
    polysize_t size; 
} specA_spXY_t;

typedef struct {
    specQ_spX_t ** Q;
    polysize_t size; 
} specQ_spXY_t;

typedef struct {
    specR_spX_t ** R;
    polysize_t size; 
} specR_spXY_t;

typedef struct {
    specA_spXY_t* An;
    specQ_spXY_t* Qn;
    specR_spXY_t* Rn;
    usfixn64 W; 
} specAQR_spXY_t;

typedef struct {
    specAQR_spXY_t ** bspecArray;
    polysize_t size;   
} specAQRArray_spXY_t;

static inline void freeSpecA_spXY(specA_spXY_t *An) {
    if (An == NULL) { 
        return; 
    }
    for (polysize_t i = 0; i < An->size; ++i) {
        freeSpecA_spX (An->A[i]);
    }
    free(An->A);
    free(An);
}

static inline void freeSpecQ_spXY(specQ_spXY_t *Qn) {
    if (Qn == NULL) { 
        return; 
    }
    for (polysize_t i = 0; i < Qn->size; ++i) {
        freeSpecQ_spX (Qn->Q[i]);
    }
    free(Qn->Q);
    free(Qn);
}

static inline void freeSpecR_spXY(specR_spXY_t *Rn) {
    if (Rn == NULL) { 
        return; 
    }
    for (polysize_t i = 0; i < Rn->size; ++i) {
        freeSpecR_spX (Rn->R[i]);
    }
    free(Rn->R);
    free(Rn);
}


static inline specA_spXY_t* createSpecA_spXY (polysize_t n, polysize_t alloc) {
    specA_spXY_t* An = (specA_spXY_t*) malloc (sizeof(specA_spXY_t));
    An->A = (specA_spX_t **) malloc (n * sizeof(specA_spX_t));
    for (polysize_t i = 0; i < n; ++i) {
        An->A[i] = createSpecA_spX (alloc);
    }
    An->size = n;
    return An;
}

static inline specQ_spXY_t* createSpecQ_spXY (polysize_t n, polysize_t alloc) {
    specQ_spXY_t* Qn = (specQ_spXY_t*) malloc (sizeof(specQ_spXY_t));
    Qn->Q = (specQ_spX_t **) malloc (n * sizeof(specQ_spX_t));
    for (polysize_t i = 0; i < n; ++i) {
        Qn->Q[i] = createSpecQ_spX (alloc);
    }
    Qn->size = n;
    return Qn;
}

static inline specR_spXY_t* createSpecR_spXY (polysize_t n, polysize_t alloc) {
    specR_spXY_t* Rn = (specR_spXY_t*) malloc (sizeof(specR_spXY_t));
    Rn->R = (specR_spX_t **) malloc (n * sizeof(specR_spX_t));
    for (polysize_t i = 0; i < n; ++i) {
        Rn->R[i] = createSpecR_spX (alloc);
    }
    Rn->size = n;
    return Rn;
}

/**
 * Computing Subresultant Chain lazily using Half-GCD over Z_{p}
 * 
 * @param a DUSP polynomial 
 * @param b DUSP polynomial 
 * @param k subresultant index 
 * @param results return a dense list of speculative subresultants 
 * @param degs [2] return the degree of speculative subresultants
 * @param results_mode 0: compute the entire S_{k+1}, S_{k} in results    
 *                     1: compute the leading coeff s_{k+1}, s_{k} in results 
 * @param A store lazy information (leading coefficients) 
 * @param Q store lazy information (quotients)
 * @param R store lazy information (remainders)
 * @param Pptr small prime pointer
 * 
 * @note A, Q must be initialized before the first call 
 */
void regularGCDSpecSRC_spX (duspoly_t *a, duspoly_t *b, polysize_t k,
                            duspolysA_t **results, polysize_t *degs, int results_mode, 
                            specA_spX_t *A, specQ_spX_t *Q, specR_spX_t *R, const Prime_ptr *Pptr);

int regularGCDbivarModularHGCDSubResInForm_withFFT_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, polysize_t *mdegs, int results_mode,  
                                           specAQR_spXY_t* specInfo, const Prime_ptr* Pptr);
int regularGCDbivarModularHGCDSubResInForm_withFFT4_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, polysize_t *mdegs, int results_mode,  
                                           specAQR_spXY_t* specInfo, const Prime_ptr* Pptr);


#ifdef __cplusplus
}
#endif

#endif

