#include "ModularPolynomial/DBSP_Support.h"

void convertFromAltArrZ_spXY_inp(const AltArrZ_t* p, dbspoly_t** q_in, const Prime_ptr* Pptr) {
    if (q_in == NULL) {
        return;
    }

    dbspoly_t* q = *q_in;

    if (isZero_AAZ(p)) {
        if (q == NULL) {
            *q_in = constPolynomial_spXY (1, 1, Pptr);
        } else {
            q->ltX = 0;
            if (q->coefsY[0] != NULL) {
                q->coefsY[0]->lt = 0;
                q->coefsY[0]->elems[0] = 0;
            } else {
                q->coefsY[0] = constPolynomial_spX(1, Pptr);
            }
        }
        return;
    }

    if (p->nvar != 2 || p->unpacked) {
        fprintf(stderr, "DUSP ERROR: Trying to convert a non-bivariate SMZP to dbspoly_t\n");
        exit(1);
    }

    int size = p->size;
    AAZElem_t* elems = p->elems;

    degree_t degs[2];
    partialDegrees_AAZ(p, degs);
    degree_t degX = degs[0];
    if (q == NULL) {
        q = allocPolynomial_spXY(degs[0]+1, degs[1]+1);
        *q_in = q;
    }
    if (q->allocX <= degX) {
        resizePolynomial_spXY(q, degX+1);
        for (int i = q->ltX+1; i <= degX; ++i) {
            q->coefsY[i] = makePolynomial_spX(degs[1]+1);
        }
    }

    duspoly_t** polys = q->coefsY;
    for (int i = 0; i <= degX; ++i) {
        if (polys[i] == NULL) {
            polys[i] = makePolynomial_spX(degs[1]+1);
        }
        polys[i]->lt = 0;
        memset(polys[i]->elems, 0, sizeof(elem_t)*polys[i]->alloc);
    }
    q->ltX = degX;

    degree_t degY;
    for (int i = 0; i < size; ++i) {
        degX = GET_NTH_EXP(elems[i].degs, EXP_1_V2, EXP_OFFSET_1_V2);
        degY = GET_NTH_EXP(elems[i].degs, EXP_2_V2, EXP_OFFSET_2_V2);
        if (polys[degX]->alloc <= degY) {
            //this should do the setting of new allocation to 0
            reallocatePolynomial_spX(polys + degX, degY+1);
        }
        polys[degX]->elems[degY] = smallprimefield_convert_in(mpz_fdiv_ui(elems[i].coef, Pptr->prime), Pptr);
        polys[degX]->lt = polys[degX]->lt < degY ? degY : polys[degX]->lt;
    }
}

void _lagrangeBasisInForm_spX (const elem_t* t, const polysize_t n, elem_t* a, const Prime_ptr* Pptr)
{
    a[0] = smallprimefield_sub (0, t[0], Pptr); // neg(t[0], pr)
    if (n==1) {
        return;
    }
    a[1] = smallprimefield_convert_in (1, Pptr);
    for (polysize_t i = 2; i <= n; i++) {
        for (polysize_t j = i; j >= 0; --j) {
            a[j] = smallprimefield_mul (a[j], smallprimefield_sub(0,t[i-1], Pptr), Pptr);
            if (j) { a[j] = smallprimefield_add (a[j], a[j-1], Pptr); }
        }
    }
}

// q = g / (x-u) mod p
// n is the size of q
// q/a is the result
elem_t _evalDivPolynomials_spX (const elem_t* a, elem_t u, elem_t* q, polysize_t n, const Prime_ptr* Pptr)
{
    elem_t u_inv = smallprimefield_inv (u, Pptr);
    elem_t eval_a = smallprimefield_sub (0, smallprimefield_mul (a[0], u_inv, Pptr), Pptr); q[0] = eval_a;
    elem_t tmp = u;
    for (polysize_t i = 1; i < n; i++) {
        q[i] = smallprimefield_mul (smallprimefield_sub (q[i-1], a[i], Pptr), u_inv, Pptr);
        eval_a = smallprimefield_add (eval_a, smallprimefield_mul (q[i], tmp, Pptr), Pptr);
        tmp = smallprimefield_mul (tmp, u, Pptr);
    }
    return eval_a;
}


/**
 * Interpolate a bivariate polynomial Z_p[x>y] from several evaluation images.
 * images[i] is the evaluation of the polynomial to interpolate at the point y=t[i].
 *
 * Both the points and the images should be in montgomery form.
 */
void interpPolyUnivarImagesInForm_spXY(const elem_t* t, duspoly_t const*const* images, int n, dbspoly_t** ret, const Prime_ptr* Pptr) {
    if (t == NULL || images == NULL || ret == NULL) {
        return;
    }

    //all images have same degree
    int degX = images[0]->lt;

    //Compute q, the lagrange basis polynomials stored densely as a matrix
    //q[i] is the ith lagrange basis poly l_i.
    //These q are re-used for all interpolations.
    elem_t* g = (elem_t*) calloc(n+1, sizeof(elem_t));
    _lagrangeBasisInForm_spX(t, n, g, Pptr);
    elem_t* q = (elem_t*) malloc(sizeof(elem_t*)*n*n);
    elem_t e;
    for (int i = 0; i < n; ++i) {
        e = _evalDivPolynomials_spX(g, t[i], q + i*n, n, Pptr);
        e = smallprimefield_inv(e, Pptr);
        for (int j = 0; j < n; ++j) {
            q[i*n + j] = smallprimefield_mul(q[i*n + j], e, Pptr);
        }
    }
    free(g);


    dbspoly_t* r = *ret;
    if (r == NULL) {
        r = allocPolynomial_spXY(degX+1, n);
        *ret = r;
    } else {
        for (int i = 0; i <= degX; ++i) {
            if (r->coefsY[i] == NULL) {
                r->coefsY[i] = makePolynomial_spX(n);
            } else if (r->coefsY[i]->alloc < n) {
                reallocatePolynomial_spX(r->coefsY + i, n);
            }
        }
    }

    elem_t curVal;
    elem_t* curInterp;
    int i, j;
    //interpolate one coefficient in Y at a time
    for (int k = 0; k <= degX; ++k) {
        curInterp = r->coefsY[k]->elems;
        memset(curInterp, 0, sizeof(elem_t)*n);
        //for each image coef of degree i
        for (i = 0; i < n; ++i) {
            //notice k is second index; we are accessing images by column
            //but, since we only access each val one it doesn't make sense
            //to transpose it. Especailly since we must iterate through entire
            //curInterp before accessing the next image value.
            curVal = images[i]->elems[k];
            //for each coef in the interpolation
            for (j = 0; j < n; ++j) {
                curInterp[j] = smallprimefield_add(curInterp[j],
                                    smallprimefield_mul(curVal, q[i*n + j], Pptr), Pptr);
            }

        }
        r->coefsY[k]->lt = n-1;
        normalize_spX(r->coefsY + k); //in case deg is not n-1
    }
    free(q);

    r->ltX = degX;
}

/**
 * Compute Partial Subresultant Chains w.r.t y in Z[y>x]
 * @note This makes use of lagrange algorithm for interpolation
 * @note There are independent for-loops (cilk_for hopes!)
 *
 * @param pSubres partial subresultant chains s.t. every poly in this sequence obtained after evaluation at y
 * @param XYsubres returned bivariate subresultant chains over Z[y>x] modulus pr
 * @param t evaluation points
 * @param n size of the sequence, that is, the degree+1 of the targeted resultant
 * @param alpha block size (deg_y(g) + deg_y (f))
 * @param Pptr small prime pointer
 */
static void _interpolatePartialSubResInForm_spXY (duspolysA_t** pSubres, biSubresPr_t* XYsubres, elem_t* t, polysize_t n, polysize_t alpha, const Prime_ptr* Pptr)
{
    // To store all images:
    elem_t* g = (elem_t*) calloc (XYsubres->n * (n+1), sizeof (elem_t));
    elem_t* q = (elem_t*) calloc (XYsubres->n * n, sizeof (elem_t));
    // for each row of the subresultant chain
    // TODO: cilk_for, cilk_grainsize optimization needed!
    for (polysize_t k = 0; k < XYsubres->n; k++) {
        if (!XYsubres->size[k]) {
            continue;
        }
        // compute delta, offset
        polysize_t delta=n-k*alpha, offset=n*k;
        _lagrangeBasisInForm_spX (t, delta, &g[offset+k], Pptr);
        // for each specialized subresultant chain at row k
        for (polysize_t i = 0; i < delta; i++) {
            if (isZero_spX (pSubres[i]->polys[k])) {
                continue;
            }
            elem_t e = _evalDivPolynomials_spX (&g[offset+k], t[i], &q[offset], delta, Pptr);
            e = smallprimefield_inv (e, Pptr);
            for (polysize_t j = 0; j < delta; j++) { // compute lagrange-basis element
                q[offset+j] = smallprimefield_mul (q[offset+j], e, Pptr);
            }
            // for the s-th coef of subresultant chain in row k, add the contribution v[s] * q[j] to it,
            // the goal is computing \sum_{s=0}{pSubres->polys[i]->lt???} v[s] * q[s]
            // TODO: cilk_for, cilk_grainsize optimization needed!
            duspoly_t* tmpZ = pSubres[i]->polys[k];
            elem_t cc;
            for (polysize_t s = 0; s <= tmpZ->lt; s++) {
                // TODO: check the uppor bound ??
                    if (tmpZ->elems[s]) {
                        cc = tmpZ->elems[s];
                        for (polysize_t j = 0; j < delta; j++) {
                            XYsubres->coefs[k][s*delta+j] = smallprimefield_add (XYsubres->coefs[k][s*delta+j] ,
                                                                smallprimefield_mul (cc, q[offset+j], Pptr), Pptr);
                        }
                    }
            }
        }
    }
    free (q);
    free (g);
}

/**
 * Bivariate Subresultant Chain  Z[y<x] but polys sorted w.r.t y,x
 *
 * @param g a bivariate polynomial (elem_t*) Z[y<x]
 * @param gpd 2Dim partial degrees of polynomial g s.t. gpd[0]:deg_y and gpd[1]:deg_x
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f s.t. fpd[0]:deg_y and fpd[1]:deg_x
 * @param subres returned subresultant chain modular pr (biSybresPr_t)
 * @param Pptr small prime pointer
 */
void biSylvSubResultantInForm_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, biSubresPr_t** subres, const Prime_ptr* Pptr)
{
	// corner cases
	// set-up:
	// n2, n3 are total degrees of g, f
	int n, n2=0, n3=0;
	n=gpd[0]+1; // deg_y(g) + 1
	for (polysize_t i = 0; i <= gpd[1]; i++) {
		for (polysize_t j = 0; j <= gpd[0]; ++j) {
			if (g[n*i+j] != 0 && (i+j) > n2) {
				n2 = i+j; }
		}
	}
	n=fpd[0]+1; // deg_y(f) + 1
	for (polysize_t i = 0; i <= fpd[1]; i++) {
		for (polysize_t j = 0; j <= fpd[0]; ++j) {
			if (f[n*i+j] != 0 && (i+j) > n3) {
				n3 = i+j; }
		}
	}
	n2 = n2*n3 + 1; // to find the min number of eval. points
	n = fpd[1]*gpd[0] + fpd[0]*gpd[1] + 1;
    // n: size of the resultant based on the Sylv Matrix (see MCA)
    // n: number of evaluation points (t)
	if (n < 2) { n = 2; }
    // fprintf (stderr, "[biSbResPr] n = %d\n", n); // TEST
	n2 = (fpd[1] > 2) ? fpd[1] : 2;
    // fprintf (stderr, "n2' = %d , n = %d\n", n2, n); // TEST
    // n2: number of subresultants to be computed based on x in Z[x>y]
	int as = gpd[1] + 1, bs = fpd[1] + 1;
    // degree in x of f(t,x) and g(t,x) mod (t-t0)
    // in evaluation of y at t.
	elem_t* a = (elem_t*) calloc (as, sizeof(elem_t)); // to store evaluated polys in Z(t_i)[x] s.t. i \in {0, ..., n-1}
	elem_t* b = (elem_t*) calloc (bs, sizeof(elem_t)); // to store evaluated polys in Z(t_i)[x] s.t. i \in {0, ..., n-1}
	elem_t* t = (elem_t*) calloc (n, sizeof(elem_t)); // to store evaluation points
	duspolysA_t** S = (duspolysA_t**) malloc (n*sizeof(duspolysA_t*)); // copies of subres over pr
    polysize_t true_sz = 0;
    polysize_t at=gpd[0]+1, bt=fpd[0]+1;
    elem_t e, tt;
    // evaluate t and compute subres of the second variable (x)
	// TODO: cilk_for, cilk_grainsize optimization needed!
    elem_t tmp;
	for (polysize_t k = 0; k < n; k++) {
		// work on evaluation point l*n + k + 1
		polysize_t ad=0, bd=0, l=0;
		while (ad < gpd[1] || bd < fpd[1]) {
			t[k] = smallprimefield_convert_in (l*n + k+1, Pptr);
			// a:
            for (polysize_t i = 0; i <= gpd[1]; i++) {
                e = g[(gpd[1]+1)*0+i];
                tmp = t[k];
                for (polysize_t j = 1; j <= gpd[0]; j++) {
                    // e += a[(gpd[1]+1)*j + i]*t[k]^j
                    tt = smallprimefield_mul (g[(gpd[1]+1)*j+i], tmp, Pptr);
                    e = smallprimefield_add (e, tt, Pptr);
                    tmp = smallprimefield_mul(tmp, t[k], Pptr);
                }
                a[i] = e;
                if (a[i]) { ad = i; }
                // a[k*as+i] = e;
                // if (a[k*as+i]) { ad = i; }
            }
			// b:
            for (polysize_t i = 0; i <= fpd[1]; i++) {
                e = f[(fpd[1]+1)*0+i];
                tmp = t[k];
                for (polysize_t j = 1; j <= fpd[0]; j++) {
                    // e += a[(gpd[1]+1)*j + i]*t[k]^j
                    tt = smallprimefield_mul (f[(fpd[1]+1)*j+i], tmp, Pptr);
                    e = smallprimefield_add (e, tt, Pptr);
                    tmp = smallprimefield_mul(tmp, t[k], Pptr);
                }
                b[i] = e;
                if (b[i]) { bd = i; }
                // b[k*bs+i] = e;
                // if (b[k*bs+i]) { bd = i; }
            }
			l++;
		}
		// Do sylvSubres to compute S[k] over Z_pr[y]
		duspoly_t* a_dusp = makePolynomial_spX (ad+1);
		duspoly_t* b_dusp = makePolynomial_spX (bd+1);
        a_dusp->elems = memcpy (a_dusp->elems, a, a_dusp->alloc * sizeof(elem_t));
        b_dusp->elems = memcpy (b_dusp->elems, b, b_dusp->alloc * sizeof(elem_t));
		a_dusp->lt = ad;
		b_dusp->lt = bd;
		polysize_t sz = 0;
        duspolysA_t* SS = NULL;

        // fprintf (stderr, "a_dusp := ");
        // printPolynomialOutForm_spX (a_dusp, Pptr);
        // fprintf (stderr, "b_dusp := ");
        // printPolynomialOutForm_spX (b_dusp, Pptr);

        // sylvSubResultantInFormA_spX (a_dusp, b_dusp, &SS, &sz, Pptr);
        subresultantChainInForm_spX (a_dusp, b_dusp, &SS, &sz, Pptr);
        if (sz > n2) {
            fprintf (stderr, "DUSP ERROR: size of returned Subres is larger than expectation!\n");
            exit (1);
        } else {
            S[k] = makePolysA_spX (n2);
            for (polysize_t i = 0; i < sz; i++) {
                if (isZero_spX (SS->polys[i])) {
                    S[k]->polys[i] = NULL;
                } else {
                    S[k]->polys[i] = deepCopyPolynomial_spX (SS->polys[i]);
                }
                // @note Don't use the following way
                // if (isZero_spX (SS->polys[sz-1-i])) {
                //     S[k]->polys[i] = NULL;
                // } else {
                //     S[k]->polys[i] = deepCopyPolynomial_spX (SS->polys[sz-1-i]);
                // }
                // @note Don't use the following way
                // if (!isZero_spX (SS->polys[i])) {
                //     if (S[k]->polys[SS->polys[i]->lt] == NULL)
                //         S[k]->polys[SS->polys[i]->lt] = deepCopyPolynomial_spX (SS->polys[i]);
                // }
            }
        }
        // for (polysize_t tt = 0; tt < n2; tt++) { // TEST
        // fprintf (stderr, "[k=%ld] Subres[%ld] := ", k, tt);
        // printPolynomialOutForm_spX (S[k]->polys[tt], Pptr); } // TEST
        freePolysA_spX (SS);
        freePolynomial_spX (&a_dusp);
        freePolynomial_spX (&b_dusp);
	}
    // for (polysize_t i = 0; i < n; i++) {
    //     fprintf (stderr, "t[%ld] = %lld\n", i, smallprimefield_convert_out (t[i], Pptr)); }
	// Final Step
	// lagrange interpolation (with a 2D distributed-data-type)
	polysize_t k = gpd[0] + fpd[0];
	biSubresPr_t* s = makeBiSubresultantInForm_spX (n2);
	for (polysize_t i = 0; i < n2; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = (n-i*k-1); // deg_y // it's always >= 0
		s->deg[i][1] = i;  // deg_x
		s->size[i] = (n-i*k-1) > -1 ? (s->deg[i][0]+1)*(s->deg[i][1]+1) : 0;
		s->coefs[i] = (elem_t*) calloc (s->size[i], sizeof(elem_t));
	}
    // for (int k = 0; k < n; k++) {
    //     for (polysize_t tt = 0; tt < n2; tt++) { // TEST
    //         fprintf (stderr, "[k=%ld] Subres[%ld] := ", k, tt);
    //         printPolynomialOutForm_spX (S[k]->polys[tt], Pptr); } // TEST
    //         fprintf (stderr, "\n\n"); } // TEST
	_interpolatePartialSubResInForm_spXY (S, s, t, n, k, Pptr);
    // for (polysize_t i = 0; i < s->n; i++) {
    //     for (polysize_t j = 0; j < s->size[i]; j++) {
    //         fprintf (stderr, "s[%ld][%ld] = %lld\n", i, j, s->coefs[i][j]); } }
	free (a); free (b); free (t);
	for (polysize_t i = 0; i < n; i++) {
		freePolysA_spX (S[i]);
	}
	free (S);
    *subres = s;
}

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
void hgcdBiSylvSubResultantInForm_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t kk, biSubresPr_t** ksubres, const Prime_ptr* Pptr)
{
	// corner cases
	// set-up:
	// n2, n3 are total degrees of g, f
	int n, n2=0, n3=0;
	n=gpd[0]+1; // deg_y(g) + 1
	for (polysize_t i = 0; i <= gpd[1]; i++) {
		for (polysize_t j = 0; j <= gpd[0]; ++j) {
			if (g[n*i+j] != 0 && (i+j) > n2) {
				n2 = i+j; }
		}
	}
	n=fpd[0]+1; // deg_y(f) + 1
	for (polysize_t i = 0; i <= fpd[1]; i++) {
		for (polysize_t j = 0; j <= fpd[0]; ++j) {
			if (f[n*i+j] != 0 && (i+j) > n3) {
				n3 = i+j; }
		}
	}
	n2 = n2*n3 + 1; // to find the min number of eval. points
	n = fpd[1]*gpd[0] + fpd[0]*gpd[1] + 1;
    // n: size of the resultant based on the Sylv Matrix (see MCA)
    // n: number of evaluation points (t)
	if (n < 2) { n = 2; }
    // fprintf (stderr, "[biSbResPr] n = %d\n", n); // TEST
	n2 = (fpd[1] > 2) ? fpd[1] : 2;
    // fprintf (stderr, "n2' = %d , n = %d\n", n2, n);
    // n2: number of subresultants to be computed based on x in Z[x>y]
	int as = gpd[1] + 1, bs = fpd[1] + 1;
    // degree in x of f(t,x) and g(t,x) mod (t-t0)
    // in evaluation of y at t.
	elem_t* a = (elem_t*) calloc (as*n, sizeof(elem_t)); // to store evaluated polys in Z(t_i)[x] s.t. i \in {0, ..., n-1}
	elem_t* b = (elem_t*) calloc (bs*n, sizeof(elem_t)); // to store evaluated polys in Z(t_i)[x] s.t. i \in {0, ..., n-1}
	elem_t* t = (elem_t*) calloc (n, sizeof(elem_t)); // to store evaluation points
	duspolysA_t** S = (duspolysA_t**) malloc (n*sizeof(duspolysA_t*)); // copies of subres over pr
    polysize_t true_sz = 0;
    polysize_t at=gpd[0]+1, bt=fpd[0]+1;
    elem_t e, tt;
    // evaluate t and compute subres of the second variable (x)
	// TODO: cilk_for, cilk_grainsize optimization needed!
	for (polysize_t k = 0; k < n; k++) {
		// work on evaluation point l*n + k + 1
		polysize_t ad=0, bd=0, l=0;
		while (ad < gpd[1] || bd < fpd[1]) {
			t[k] = smallprimefield_convert_in (l*n + k+1, Pptr);
			// a:
            for (polysize_t i = 0; i <= gpd[1]; i++) {
                e = 0;
                for (polysize_t j = 0; j <= gpd[0]; j++) {
                    // e += a[(gpd[1]+1)*j + i]*t[k]^j
                    tt = smallprimefield_mul (g[(gpd[1]+1)*j+i], smallprimefield_exp (t[k], j, Pptr), Pptr);
                    e = smallprimefield_add (e, tt, Pptr);
                }
                a[k*as+i] = e;
                if (a[k*as+i]) { ad = i; }
            }
			// b:
            for (polysize_t i = 0; i <= fpd[1]; i++) {
                e = 0;
                for (polysize_t j = 0; j <= fpd[0]; j++) {
                    // e += a[(gpd[1]+1)*j + i]*t[k]^j
                    tt = smallprimefield_mul (f[(fpd[1]+1)*j+i], smallprimefield_exp (t[k], j, Pptr), Pptr);
                    e = smallprimefield_add (e, tt, Pptr);
                }
                b[k*bs+i] = e;
                if (b[k*bs+i]) { bd = i; }
            }
			l++;
		}
		// Do sylvSubres to compute S[k] over Z_pr[y]
		duspoly_t* a_dusp = makePolynomial_spX (ad+1);
		duspoly_t* b_dusp = makePolynomial_spX (bd+1);
		// a_dusp->elems = (elem_t*) malloc (a_dusp->alloc * sizeof(elem_t)); // &a[k*as];
        a_dusp->elems = memcpy (a_dusp->elems, &a[k*as], a_dusp->alloc * sizeof(elem_t));
		// b_dusp->elems = (elem_t*) malloc (b_dusp->alloc * sizeof(elem_t)); // &b[k*bs];
        b_dusp->elems = memcpy (b_dusp->elems, &b[k*bs], b_dusp->alloc * sizeof(elem_t));
		a_dusp->lt = ad;
		b_dusp->lt = bd;
		polysize_t sz = 0;
        duspolysA_t* SS = NULL;

        // fprintf (stderr, "a_dusp := ");
        // printPolynomialOutForm_spX (a_dusp, Pptr);
        // fprintf (stderr, "b_dusp := ");
        // printPolynomialOutForm_spX (b_dusp, Pptr);

        // sylvSubResultantInFormA_spX (a_dusp, b_dusp, &SS, &sz, Pptr);
        // subresultantChainInForm_spX (a_dusp, b_dusp, &SS, &sz, Pptr);
        hgcdSubResultantInFormA_spX (a_dusp, b_dusp, kk, &SS, &sz, Pptr);
        if (sz > n2) {
            fprintf (stderr, "DUSP ERROR: size of returned Subres is larger than expectation!\n");
            exit (1);
        } else {
            S[k] = makePolysA_spX (n2);
            for (polysize_t i = 0; i < sz; i++) {
                if (isZero_spX (SS->polys[i])) {
                    S[k]->polys[i] = NULL;
                } else {
                    S[k]->polys[i] = deepCopyPolynomial_spX (SS->polys[i]);
                }
                // @note Don't use following way
                // if (isZero_spX (SS->polys[sz-1-i])) {
                //     S[k]->polys[i] = NULL;
                // } else {
                //     S[k]->polys[i] = deepCopyPolynomial_spX (SS->polys[sz-1-i]);
                // }
                // @note Don't use following way
                // if (!isZero_spX (SS->polys[i])) {
                //     if (S[k]->polys[SS->polys[i]->lt] == NULL)
                //         S[k]->polys[SS->polys[i]->lt] = deepCopyPolynomial_spX (SS->polys[i]);
                // }
            }
        }
        // for (polysize_t tt = 0; tt < n2; tt++) { // TEST
        // fprintf (stderr, "[k=%ld] Subres[%ld] := ", k, tt);
        // printPolynomialOutForm_spX (S[k]->polys[tt], Pptr); } // TEST
        freePolysA_spX (SS);
        freePolynomial_spX (&a_dusp);
        freePolynomial_spX (&b_dusp);
	}
    // for (polysize_t i = 0; i < n; i++) {
    //     fprintf (stderr, "t[%ld] = %lld\n", i, smallprimefield_convert_out (t[i], Pptr)); }
	// Final Step
	// lagrange interpolation (with a 2D distributed-data-type)
	polysize_t k = gpd[0] + fpd[0];
	biSubresPr_t* s = makeBiSubresultantInForm_spX (n2);
	for (polysize_t i = 0; i < n2; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = (n-i*k-1); // deg_y // it's always >= 0
		s->deg[i][1] = i;  // deg_x
		s->size[i] = (n-i*k-1) > -1 ? (s->deg[i][0]+1)*(s->deg[i][1]+1) : 0;
		s->coefs[i] = (elem_t*) calloc (s->size[i], sizeof(elem_t));
	}
    // for (int k = 0; k < n; k++) {
    //     for (polysize_t tt = 0; tt < n2; tt++) { // TEST
    //         fprintf (stderr, "[k=%ld] Subres[%ld] := ", k, tt);
    //         printPolynomialOutForm_spX (S[k]->polys[tt], Pptr); } // TEST
    //         fprintf (stderr, "\n\n"); } // TEST
	_interpolatePartialSubResInForm_spXY (S, s, t, n, k, Pptr);
    // for (polysize_t i = 0; i < s->n; i++) {
    //     for (polysize_t j = 0; j < s->size[i]; j++) {
    //         fprintf (stderr, "s[%ld][%ld] = %lld\n", i, j, s->coefs[i][j]); } }
	free (a); free (b); free (t);
	for (polysize_t i = 0; i < n; i++) {
		freePolysA_spX (S[i]);
	}
	free (S);
    *ksubres = s;
}

////////////////////////////////
// FFT-based Subresultant Chains
////////////////////////////////

/**
 * static function to compute modularSubres from two usfixn64 arrays
 * representing modular univariate polynomials in Montgomery form
 */
static usfixn64* _modularSubres_from_usfixn64 (usfixn64* g, polysize_t gpdx, usfixn64* f, polysize_t fpdx,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, const Prime_ptr* Pptr)
{
        duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
		duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { a_dusp->elems[i] = (elem_t) g[i]; }
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { b_dusp->elems[i] = (elem_t) f[i]; }
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);
        // fprintf(stderr, "TEST: a_dusp = "); // TEST
        // printPolynomial_spX(a_dusp, Pptr); // TEST
        // fprintf(stderr, "TEST: b_dusp = ");  // TEST
        // printPolynomial_spX(b_dusp, Pptr); // TEST

        polysize_t sz = 0;
        duspolysA_t* S = NULL;
        subresultantChainInForm_spX (a_dusp, b_dusp, &S, &sz, Pptr);

        // fprintf(stderr, "TEST: subres_sz = %ld\n", sz); // TEST
        // for (polysize_t i = 0; i < sz; i++) { // TEST
        //     fprintf(stderr, "TEST: subres[%ld] = ", i); // TEST
        //     printPolynomial_spX(S->polys[i], Pptr); // TEST
        // } // TEST
        if (S->polys[0] == NULL) {
            S->polys[0] = zeroPolynomial_spX();
        }
        polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
            if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            lcl_chain_alloc_sizes[i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        usfixn64* subres_array = (usfixn64*) calloc(sz * total_alloc_size, sizeof(usfixn64));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), lcl_chain_alloc_sizes[i] * sizeof(usfixn64));
            start += lcl_chain_alloc_sizes[i];
        }
        *chain_alloc_sizes = lcl_chain_alloc_sizes;
        *subres_size = sz;

        // ///////////// ADD a, b to the chain ////////////
        // @note Must be commented... for DEBUGGING only!

        // polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc ((sz+2)*sizeof(polysize_t));
        // polysize_t total_alloc_size = 0;
        // for(polysize_t i = 0; i < sz; i++) {
        //     lcl_chain_alloc_sizes[i+2] = S->polys[i]->alloc;
        //     total_alloc_size += lcl_chain_alloc_sizes[i+2];
        // }
        // if (gpdx >= fpdx){
        //     lcl_chain_alloc_sizes[0] = a_dusp->alloc;
        //     lcl_chain_alloc_sizes[1] = b_dusp->alloc;
        //     total_alloc_size += lcl_chain_alloc_sizes[0] + lcl_chain_alloc_sizes[1];
        // } else {
        //     lcl_chain_alloc_sizes[0] = b_dusp->alloc;
        //     lcl_chain_alloc_sizes[1] = a_dusp->alloc;
        //     total_alloc_size += lcl_chain_alloc_sizes[0] + lcl_chain_alloc_sizes[1];
        // }
        // usfixn64* subres_array = (usfixn64*) calloc((sz+2) * total_alloc_size, sizeof(usfixn64));
        // polysize_t start = 0;
        // if (gpdx >= fpdx) {
        //     memcpy(&(subres_array[start]), (a_dusp->elems), lcl_chain_alloc_sizes[0] * sizeof(usfixn64));
        //     start += lcl_chain_alloc_sizes[0];
        //     memcpy(&(subres_array[start]), (b_dusp->elems), lcl_chain_alloc_sizes[1] * sizeof(usfixn64));
        //     start += lcl_chain_alloc_sizes[1];
        // } else {
        //     memcpy(&(subres_array[start]), (b_dusp->elems), lcl_chain_alloc_sizes[0] * sizeof(usfixn64));
        //     start += lcl_chain_alloc_sizes[0];
        //     memcpy(&(subres_array[start]), (a_dusp->elems), lcl_chain_alloc_sizes[1] * sizeof(usfixn64));
        //     start += lcl_chain_alloc_sizes[1];
        // }
        // for(polysize_t i = 2; i < sz+2; i++) {
        //     memcpy(&(subres_array[start]), (S->polys[i-2]->elems), lcl_chain_alloc_sizes[i] * sizeof(usfixn64));
        //     start += lcl_chain_alloc_sizes[i];
        // }
        // *chain_alloc_sizes = lcl_chain_alloc_sizes;
        // *subres_size = sz+2;
        // ///////////// /////////////////// ////////////

        freePolynomial_spX(&a_dusp); // TODO: we can do better
        freePolynomial_spX(&b_dusp);
        freePolysA_spX(S);
        return subres_array;
}

static elem_t* _modularSubres_from_elems (elem_t* g, polysize_t gpdx, elem_t* f, polysize_t fpdx,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, const Prime_ptr* Pptr)
{
        duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
		duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { a_dusp->elems[i] = g[i]; }
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { b_dusp->elems[i] = f[i]; }
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);
        // fprintf(stderr, "TEST: a_dusp = "); // TEST
        // printPolynomial_spX(a_dusp, Pptr); // TEST
        // fprintf(stderr, "TEST: b_dusp = ");  // TEST
        // printPolynomial_spX(b_dusp, Pptr); // TEST

        polysize_t sz = 0;
        duspolysA_t* S = NULL;
        subresultantChainInForm_spX (a_dusp, b_dusp, &S, &sz, Pptr);

        // fprintf(stderr, "TEST: subres_sz = %ld\n", sz); // TEST
        // for (polysize_t i = 0; i < sz; i++) { // TEST
        //     fprintf(stderr, "TEST: subres[%ld] = ", i); // TEST
        //     printPolynomial_spX(S->polys[i], Pptr); // TEST
        // } // TEST
        if (S->polys[0] == NULL) {
            S->polys[0] = zeroPolynomial_spX();
        }
        polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
            if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            lcl_chain_alloc_sizes[i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        elem_t* subres_array = (elem_t*) calloc(sz * total_alloc_size, sizeof(elem_t));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), lcl_chain_alloc_sizes[i] * sizeof(elem_t));
            start += lcl_chain_alloc_sizes[i];
        }
        *chain_alloc_sizes = lcl_chain_alloc_sizes;
        *subres_size = sz;

        freePolynomial_spX(&a_dusp); 
        freePolynomial_spX(&b_dusp);
        freePolysA_spX(S);
        return subres_array;
}

/**
 * see _compute_DFTInv_modularSubres
 */
static usfixn64** _compute_DFTInv_modularSubres_inputs (usfixn64** SS, polysize_t n, polysize_t subres_size,
                                                        polysize_t alpha, polysize_t* max_chain_alloc_sizes)
{
    usfixn64* subres;
    polysize_t start, loc;
    usfixn64** result = (usfixn64**) malloc(subres_size * sizeof(usfixn64*));
    for(int j = 0; j < subres_size; j++) {
        result[j] = (usfixn64*) calloc(n * max_chain_alloc_sizes[j], sizeof(usfixn64));
    }
    polysize_t delta;
    start = 0;
    for(polysize_t k = 0; k < subres_size; k++) {
        // delta = n - k*alpha; // TODO: to use this optimization we need truncated-FFT 
        loc = max_chain_alloc_sizes[k];
        for(polysize_t i = 0; i < n; i++) {
            memcpy(&(result[k][i*loc]), &(SS[i][start]), loc * sizeof(usfixn64)); // n * loc
        }
        start += loc;
    }

    return result;
}

static elem_t** _compute_DFT4Inv_modularSubres_inputs (elem_t** SS, polysize_t n, polysize_t subres_size,
                                                        polysize_t alpha, polysize_t* max_chain_alloc_sizes)
{
    elem_t* subres;
    polysize_t start, loc;
    elem_t** result = (elem_t**) malloc(subres_size * sizeof(elem_t*));
    for(int j = 0; j < subres_size; j++) {
        result[j] = (elem_t*) calloc(n * max_chain_alloc_sizes[j], sizeof(elem_t));
    }
    polysize_t delta;
    start = 0;
    for(polysize_t k = 0; k < subres_size; k++) {
        // delta = n - k*alpha; // TODO: to use this optimization we need truncated-FFT 
        loc = max_chain_alloc_sizes[k];
        for(polysize_t i = 0; i < n; i++) {
            memcpy(&(result[k][i*loc]), &(SS[i][start]), loc * sizeof(elem_t)); // n * loc
        }
        start += loc;
    }

    return result;
}


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
                                        const Prime_ptr* Pptr)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    // evaluate_subres_n[i] = n * max_chain_alloc_sizes[i] for i = 0, ..., subres_size-1
    usfixn64** evaluated_subres_n = _compute_DFTInv_modularSubres_inputs(SS, n, subres_size, alpha, max_chain_alloc_sizes);
    for(polysize_t i = 0; i < subres_size; i++) {
       evaluated_subres_n[i] = transpose2D_spX(evaluated_subres_n[i], n, max_chain_alloc_sizes[i]); // max_chain_alloc_sizes[i] * n
       evaluated_subres_n[i] = bivarRowWiseDFTInv_spX(evaluated_subres_n[i], max_chain_alloc_sizes[i]-1, n, K, e, w_inv, n_inv, P);
    }
    return evaluated_subres_n;
}

elem_t** _compute_DFT4Inv_modularSubres (elem_t** SS, polysize_t n, polysize_t r, elem_t* OmegaInv, elem_t n_inv,
                                        polysize_t subres_size, polysize_t alpha, polysize_t* max_chain_alloc_sizes,
                                        const Prime_ptr* Pptr)
{
    // evaluate_subres_n[i] = n * max_chain_alloc_sizes[i] for i = 0, ..., subres_size-1
    elem_t** evaluated_subres_n = _compute_DFT4Inv_modularSubres_inputs(SS, n, subres_size, alpha, max_chain_alloc_sizes);
    for(polysize_t i = 0; i < subres_size; i++) {
       evaluated_subres_n[i] = transpose2D_elem_spX(evaluated_subres_n[i], n, max_chain_alloc_sizes[i]); // max_chain_alloc_sizes[i] * n
       evaluated_subres_n[i] = bivarRowWiseDFT4Inv_spX(evaluated_subres_n[i], max_chain_alloc_sizes[i]-1, n, r, OmegaInv, n_inv, Pptr);
    }
    return evaluated_subres_n;
}

/**
 * FFT-based Bivariate Subresultant Chain  Z[y<x] but polys sorted w.r.t y,x
 *
 * @param g a bivariate polynomial (elem_t*) Z[y<x]
 * @param gpd 2Dim partial degrees of polynomial g s.t. gpd[0]:deg_y and gpd[1]:deg_x
 * @param f a bivariate polynomial (dbspoly_t) s.t. deg(g, x) > deg(f, x)
 * @param fpd 2Dim partial degrees of polynomial f s.t. fpd[0]:deg_y and fpd[1]:deg_x
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
                                           biSubresPr_t** subres, const Prime_ptr* Pptr)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    usfixn64* g_eval = convert_bivarElemPoly_to_usfixn64Poly(g, gpd, n); // (gpd[1]+1) * n
    usfixn64* f_eval = convert_bivarElemPoly_to_usfixn64Poly(f, fpd, n); // (fpd[1]+1) * n
    g_eval = bivarRowWiseDFT_spX(g_eval, gpd[1], n, K, e, w, P); // row-wise FFT(g_eval)
    f_eval = bivarRowWiseDFT_spX(f_eval, fpd[1], n, K, e, w, P); // row-wise FFT(f_eval)
    g_eval = transpose2D_spX(g_eval, gpd[1]+1, n); // n * (gpd[1]+1)
    f_eval = transpose2D_spX(f_eval, fpd[1]+1, n); // n * (fpd[1]+1)
    usfixn64** subres_arrays = (usfixn64**) malloc(n * sizeof(usfixn64*));
    polysize_t *chain_alloc_sizes = NULL, *max_chain_alloc_sizes = NULL;
    polysize_t subres_size, test_subres_size = 0;
    // cilk-for
    for(polysize_t i = 0; i < n; i++) {
        subres_size = 0;
        subres_arrays[i] = _modularSubres_from_usfixn64(&(g_eval[i * (gpd[1]+1)]), gpd[1],
                                                        &(f_eval[i * (fpd[1]+1)]), fpd[1],
                                                        &subres_size, &chain_alloc_sizes, Pptr);
        // check the goodness of omega OR initialize test_subres_size and max_chain_alloc_sizes
        if (test_subres_size != subres_size) {
            if (test_subres_size == 0 && i == 0) {
                test_subres_size = subres_size;
                max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
                for(polysize_t k = 0; k < test_subres_size; k++) {
                    max_chain_alloc_sizes[k] = chain_alloc_sizes[k];
                }
            }
            else {
                free(g_eval);
                free(f_eval);
                for (polysize_t j = 0; j < i; j++) {
                    free(subres_arrays[j]);
                }
                free(subres_arrays);
                if (max_chain_alloc_sizes != NULL)
                    free(max_chain_alloc_sizes);
                if (chain_alloc_sizes != NULL)
                    free(chain_alloc_sizes);
                return 0;
            }
        }
        // check the goodness of omega
        for(polysize_t k = 0; k < test_subres_size &&
            max_chain_alloc_sizes[k] != chain_alloc_sizes[k]; k++) {
            free(g_eval);
            free(f_eval);
            for (polysize_t j = 0; j < i; j++) {
                free(subres_arrays[j]);
            }
            free(subres_arrays);
            if (max_chain_alloc_sizes != NULL)
                free(max_chain_alloc_sizes);
            if (chain_alloc_sizes != NULL)
                free(chain_alloc_sizes);
            return 0;
        }
        if (chain_alloc_sizes != NULL) {
            free(chain_alloc_sizes);
            chain_alloc_sizes = NULL;
        } 
    }
    polysize_t alpha = gpd[0] + fpd[0];
    usfixn64** polys = _compute_DFTInv_modularSubres(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
	biSubresPr_t* s = makeBiSubresultantInForm_spX (test_subres_size);
	for (polysize_t i = 0; i < test_subres_size; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = n-1; // (n-i*alpha-1); // deg_y
		s->deg[i][1] = max_chain_alloc_sizes[i]-1; // deg_x
		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1); // (n-i*alpha-1) > -1 ? (s->deg[i][0]+1)*(s->deg[i][1]+1) : 0;
		s->coefs[i] = (elem_t*) malloc (s->size[i] * sizeof(elem_t));
        memcpy(s->coefs[i], polys[i], s->size[i]*sizeof(elem_t));
	}
    *subres = s;
    free(g_eval); free(f_eval);
    for (polysize_t j = 0; j < n; j++) {
        free(subres_arrays[j]);
    }
    free(subres_arrays);
    for (polysize_t j = 0; j < test_subres_size; j++) {
        free(polys[j]);
    }
    free(polys);
    return 1;
}

int bivarModularSubResInForm_withFFT4_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd,
                                           polysize_t n, polysize_t r, elem_t* Omega, elem_t* OmegaInv, elem_t n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr)
{
    elem_t* g_eval = convert_bivarElemPoly_to_elemPoly(g, gpd, n); // (gpd[1]+1) * n
    elem_t* f_eval = convert_bivarElemPoly_to_elemPoly(f, fpd, n); // (fpd[1]+1) * n
    g_eval = bivarRowWiseDFT4_spX(g_eval, gpd[1], n, r, Omega, Pptr); // row-wise FFT(g_eval)
    f_eval = bivarRowWiseDFT4_spX(f_eval, fpd[1], n, r, Omega, Pptr); // row-wise FFT(f_eval)
    g_eval = transpose2D_elem_spX(g_eval, gpd[1]+1, n); // n * (gpd[1]+1)
    f_eval = transpose2D_elem_spX(f_eval, fpd[1]+1, n); // n * (fpd[1]+1)
    elem_t** subres_arrays = (elem_t**) malloc(n * sizeof(elem_t*));
    polysize_t *chain_alloc_sizes = NULL, *max_chain_alloc_sizes = NULL;
    polysize_t subres_size, test_subres_size = 0;
    // cilk-for
    for(polysize_t i = 0; i < n; i++) {
        subres_size = 0;
        subres_arrays[i] = _modularSubres_from_elems(&(g_eval[i * (gpd[1]+1)]), gpd[1],
                                                        &(f_eval[i * (fpd[1]+1)]), fpd[1],
                                                        &subres_size, &chain_alloc_sizes, Pptr);
        // check the goodness of omega OR initialize test_subres_size and max_chain_alloc_sizes
        if (test_subres_size != subres_size) {
            if (test_subres_size == 0 && i == 0) {
                test_subres_size = subres_size;
                max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
                for(polysize_t k = 0; k < test_subres_size; k++) {
                    max_chain_alloc_sizes[k] = chain_alloc_sizes[k];
                }
            }
            else {
                free(g_eval);
                free(f_eval);
                for (polysize_t j = 0; j < i; j++) {
                    free(subres_arrays[j]);
                }
                free(subres_arrays);
                if (max_chain_alloc_sizes != NULL)
                    free(max_chain_alloc_sizes);
                if (chain_alloc_sizes != NULL)
                    free(chain_alloc_sizes);
                return 0;
            }
        }
        // check the goodness of omega
        for(polysize_t k = 0; k < test_subres_size &&
            max_chain_alloc_sizes[k] != chain_alloc_sizes[k]; k++) {
            free(g_eval);
            free(f_eval);
            for (polysize_t j = 0; j < i; j++) {
                free(subres_arrays[j]);
            }
            free(subres_arrays);
            if (max_chain_alloc_sizes != NULL)
                free(max_chain_alloc_sizes);
            if (chain_alloc_sizes != NULL)
                free(chain_alloc_sizes);
            return 0;
        }
        if (chain_alloc_sizes != NULL) {
            free(chain_alloc_sizes);
            chain_alloc_sizes = NULL;
        }
    }
    polysize_t alpha = gpd[0] + fpd[0];
    elem_t** polys = _compute_DFT4Inv_modularSubres(subres_arrays, n, r, OmegaInv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
	biSubresPr_t* s = makeBiSubresultantInForm_spX (test_subres_size);
	for (polysize_t i = 0; i < test_subres_size; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = n-1; // (n-i*alpha-1); // deg_y
		s->deg[i][1] = max_chain_alloc_sizes[i]-1; // deg_x
		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1); // (n-i*alpha-1) > -1 ? (s->deg[i][0]+1)*(s->deg[i][1]+1) : 0;
		s->coefs[i] = (elem_t*) malloc (s->size[i] * sizeof(elem_t));
        memcpy(s->coefs[i], polys[i], s->size[i]*sizeof(elem_t));
	}
    *subres = s;
    free(g_eval); free(f_eval);
    for (polysize_t j = 0; j < n; j++) {
        free(subres_arrays[j]);
    }
    free(subres_arrays);
    for (polysize_t j = 0; j < test_subres_size; j++) {
        free(polys[j]);
    }
    free(polys);
    return 1;
}

static usfixn64* _modularHGCDSubres_from_usfixn64 (usfixn64* g, polysize_t gpdx, usfixn64* f, polysize_t fpdx, polysize_t sub_k,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, const Prime_ptr* Pptr)
{
        duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
		duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { a_dusp->elems[i] = (elem_t) g[i]; }
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { b_dusp->elems[i] = (elem_t) f[i]; }
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);
        // fprintf(stderr, "TEST: a_dusp = "); // TEST
        // printPolynomial_spX(a_dusp, Pptr); // TEST
        // fprintf(stderr, "TEST: b_dusp = ");  // TEST
        // printPolynomial_spX(b_dusp, Pptr); // TEST

        polysize_t sz = 0;
        duspolysA_t* S = NULL;
        hgcdSubResultantInFormA_spX (a_dusp, b_dusp, sub_k, &S, &sz, Pptr);
        polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
         if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            lcl_chain_alloc_sizes[i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        usfixn64* subres_array = (usfixn64*) calloc(sz * total_alloc_size, sizeof(usfixn64));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), lcl_chain_alloc_sizes[i] * sizeof(usfixn64));
            start += lcl_chain_alloc_sizes[i];
        }
        *chain_alloc_sizes = lcl_chain_alloc_sizes;
        *subres_size = sz;

        freePolynomial_spX(&a_dusp); // TODO: we can do better
        freePolynomial_spX(&b_dusp); // TODO: we can do better
        freePolysA_spX(S);
        return subres_array;
}

static elem_t* _modularHGCDSubres_from_elems (elem_t* g, polysize_t gpdx, elem_t* f, polysize_t fpdx, polysize_t sub_k,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, const Prime_ptr* Pptr)
{ 
        duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
		duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { a_dusp->elems[i] = (elem_t) g[i]; }
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { b_dusp->elems[i] = (elem_t) f[i]; }
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);

        polysize_t sz = 0;
        duspolysA_t* S = NULL;
        hgcdSubResultantInFormA_spX (a_dusp, b_dusp, sub_k, &S, &sz, Pptr);
        polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
         if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            lcl_chain_alloc_sizes[i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        elem_t* subres_array = (elem_t*) calloc(sz * total_alloc_size, sizeof(elem_t));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), lcl_chain_alloc_sizes[i] * sizeof(elem_t));
            start += lcl_chain_alloc_sizes[i];
        }
        *chain_alloc_sizes = lcl_chain_alloc_sizes;
        *subres_size = sz;

        freePolynomial_spX(&a_dusp);
        freePolynomial_spX(&b_dusp);
        freePolysA_spX(S);
        return subres_array;
}


int bivarModularHGCDSubResInForm_withFFT_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    usfixn64* g_eval = convert_bivarElemPoly_to_usfixn64Poly(g, gpd, n); // (gpd[1]+1) * n
    usfixn64* f_eval = convert_bivarElemPoly_to_usfixn64Poly(f, fpd, n); // (fpd[1]+1) * n
    g_eval = bivarRowWiseDFT_spX(g_eval, gpd[1], n, K, e, w, P); // row-wise FFT(g_eval)
    f_eval = bivarRowWiseDFT_spX(f_eval, fpd[1], n, K, e, w, P); // row-wise FFT(f_eval)
    g_eval = transpose2D_spX(g_eval, gpd[1]+1, n); // n * (gpd[1]+1)
    f_eval = transpose2D_spX(f_eval, fpd[1]+1, n); // n * (fpd[1]+1)
    usfixn64** subres_arrays = (usfixn64**) malloc(n * sizeof(usfixn64*));
    polysize_t *chain_alloc_sizes = NULL, *max_chain_alloc_sizes = NULL;
    polysize_t subres_size, test_subres_size = 0;
    for(polysize_t i = 0; i < n; i++) {
        subres_size = 0;
        subres_arrays[i] = _modularHGCDSubres_from_usfixn64(&(g_eval[i * (gpd[1]+1)]), gpd[1],
                                                        &(f_eval[i * (fpd[1]+1)]), fpd[1], k_sub,
                                                        &subres_size, &chain_alloc_sizes, Pptr);
        // check the goodness of omega OR initialize test_subres_size and max_chain_alloc_sizes
        if (test_subres_size != subres_size) {
            if (test_subres_size == 0 && i == 0) {
                test_subres_size = subres_size;
                max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
                for(polysize_t k = 0; k < test_subres_size; k++) {
                    max_chain_alloc_sizes[k] = chain_alloc_sizes[k];
                }
            }
            else {
                free(g_eval); free(f_eval);
                for (polysize_t j = 0; j < i; j++) { free(subres_arrays[j]);
                } free(subres_arrays);
                if (max_chain_alloc_sizes != NULL) {free(max_chain_alloc_sizes);}
                if (chain_alloc_sizes != NULL) {free(chain_alloc_sizes);}
                return 0;
            }
        }
        // check the goodness of omega
        for(polysize_t k = 0; k < test_subres_size &&
            max_chain_alloc_sizes[k] != chain_alloc_sizes[k]; k++) {
            free(g_eval); free(f_eval);
            for (polysize_t j = 0; j < i; j++) { free(subres_arrays[j]);
            } free(subres_arrays);
            if (max_chain_alloc_sizes != NULL) {free(max_chain_alloc_sizes);}
            if (chain_alloc_sizes != NULL) {free(chain_alloc_sizes);}
            return 0;
        }
        if (chain_alloc_sizes != NULL) {
            free(chain_alloc_sizes);
            chain_alloc_sizes = NULL;
        }
    }
    polysize_t alpha = gpd[0] + fpd[0];
    // bivarSubresFFT_t* subres_out = (bivarSubresFFT_t*) malloc(sizeof(bivarSubresFFT_t));
    // subres_out->polys = _compute_DFTInv_modularSubres(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
    // subres_out->chain_size = test_subres_size;
    // subres_out->alloc_sizes = max_chain_alloc_sizes;
    // *subres = subres_out;

    usfixn64** polys = _compute_DFTInv_modularSubres(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
    // fprintf(stderr, "TEST: subres_size = %ld\n", test_subres_size); // TEST
	biSubresPr_t* s = makeBiSubresultantInForm_spX (test_subres_size);
	for (polysize_t i = 0; i < test_subres_size; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = n-1; // deg_y
		s->deg[i][1] = max_chain_alloc_sizes[i]-1; // deg_x
		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1);
		s->coefs[i] = (elem_t*) malloc (s->size[i] * sizeof(elem_t));
        memcpy(s->coefs[i], polys[i], s->size[i]*sizeof(elem_t));
        // for(polysize_t j = 0; j < s->size[i]; j++) {
        //     s->coefs[i][j] = (elem_t) polys[i][j];
        // }
	}
    *subres = s;
    free(g_eval); free(f_eval);
    for (polysize_t j = 0; j < n; j++) { free(subres_arrays[j]);
    }
    free(subres_arrays);
    for (polysize_t j = 0; j < test_subres_size; j++) {
        free(polys[j]);
    }
    free(polys);
    return 1;
}

/**
 * static function to compute modular hgcd-based Subres from two usfixn64 arrays
 * representing modular univariate polynomials in Montgomery form
 */
static usfixn64* _RegularGCDModularSubres_from_usfixn64_spXY (usfixn64* g, polysize_t gpdx, usfixn64* f, polysize_t fpdx, polysize_t sub_k,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, polysize_t *mdegs, int results_mode,  
                                            specAQR_spXY_t* specInfo, polysize_t idx, const Prime_ptr* Pptr)
{
        duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
		duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { a_dusp->elems[i] = (elem_t) g[i]; }
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { b_dusp->elems[i] = (elem_t) f[i]; }
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);
        // fprintf(stderr, "TEST: a_dusp = "); // TEST
        // printPolynomial_spX(a_dusp, Pptr); // TEST
        // fprintf(stderr, "TEST: b_dusp = ");  // TEST
        // printPolynomial_spX(b_dusp, Pptr); // TEST

        duspolysA_t* S = NULL;
        regularGCDSpecSRC_spX (a_dusp, b_dusp, sub_k, &S, mdegs, results_mode, 
                            specInfo->An->A[idx], specInfo->Qn->Q[idx], specInfo->Rn->R[idx], Pptr);
        if (S == NULL) {
            S = makePolysA_spX (1);
            S->polys[0] = zeroPolynomial_spX ();
        }
        polysize_t sz = S->size;
        polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
         if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            lcl_chain_alloc_sizes[i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        usfixn64* subres_array = (usfixn64*) calloc(sz * total_alloc_size, sizeof(usfixn64));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), lcl_chain_alloc_sizes[i] * sizeof(usfixn64));
            start += lcl_chain_alloc_sizes[i];
        }
        *chain_alloc_sizes = lcl_chain_alloc_sizes;
        *subres_size = sz;

        freePolynomial_spX(&a_dusp); 
        freePolynomial_spX(&b_dusp);
        freePolysA_spX(S);
        return subres_array;
}


static polysize_t _findRemainderIndex_spX(polysize_t *n, polysize_t n_size, polysize_t k)
{
    for (polysize_t h = 0; h < n_size -1; ++h) {
        if (n[h] > k && n[h+1] <= k) {
            return h;
        }
    }
    return n_size-2;
}

void regularGCDSpecSRC_spX (duspoly_t *a, duspoly_t *b, polysize_t k,
                            duspolysA_t **results, polysize_t *degs, int results_mode, 
                            specA_spX_t *A, specQ_spX_t *Q, specR_spX_t *R, const Prime_ptr *Pptr)
{
    duspolysA_t* sres = NULL;
    if (a == NULL || b == NULL || A == NULL || Q == NULL || R == NULL) { 
        return; 
    } else if (a->lt < b->lt) {
        regularGCDSpecSRC_spX (b, a, k, results, degs, results_mode, A, Q, R, Pptr);
        return;
    } else if (k < 0 || k > b->lt) {
        sres = makePolysA_spX(2);
        if (results_mode == 0) {
            sres->polys[0] = deepCopyPolynomial_spX (a);
            sres->polys[1] = deepCopyPolynomial_spX (b);
        } else {
            if (isZero_spX (a)) {
                sres->polys[0] = zeroPolynomial_spX ();
            } else {
                sres->polys[0] = constPolynomial_spX (a->elems[0], Pptr);
            }
            if (isZero_spX (b)) {
                sres->polys[1] = zeroPolynomial_spX ();
            } else {
                sres->polys[1] = constPolynomial_spX (b->elems[0], Pptr);
            }
        }
        if (results != NULL) {
            *results = sres;
        } else {
            freePolysA_spX (sres);
        }
        if (degs != NULL) {
            degs[0] = degPolynomial_spX (a);
            degs[1] = degPolynomial_spX (b);
        }
        return; 
    }

    int are_A_and_Q_already_computed;
    if (A->size == 0 && Q->size == 0) {
        are_A_and_Q_already_computed = 0;
    } else if (A->size != 0 && Q->size != 0) {
        are_A_and_Q_already_computed = 1;
    } else {
        fprintf(stderr, "DUSP Error: expected A and Q with compatible sizes, but received A->size = %ld, Q->size = %ld\n", A->size, Q->size);
        exit(EXIT_FAILURE);
    }
    if (k == 0 || !are_A_and_Q_already_computed) {
        sres = makePolysA_spX(2);
        elem_t* lcrem = NULL; // TODO: _HGCD... should use specA_spX_t data-struct instead
        _HGCDInFormAtMaxK_spX (a, b, 0, &(R->r[1]), &(R->r[0]), &lcrem, Q, Pptr);
        int idx = 2, min_deg = MIN_spX(a->lt, b->lt); // lcrem size!
        A->n[0] = a->lt; A->lc[0] = a->elems[a->lt]; 
        A->n[1] = b->lt; A->lc[1] = b->elems[b->lt]; 
        for(int i = min_deg-1; i >= 0; --i) { 
            if (lcrem[i]) {
                A->lc[idx] = lcrem[i]; A->n[idx] = i;
                idx++; 
            }
        }
        if (lcrem != NULL) {
            free(lcrem);
        }
        if (isZero_spX (R->r[0])) {
            A->lc[idx] = 0; A->n[idx] = 0;
            idx++;
        }
        A->size = idx; 
        R->size = 2;
    }

    polysize_t h;
    if (k == 0) {
        h = A->size - 2;
    } else {
        h = _findRemainderIndex_spX(A->n, A->size, k); // n_{h} > k and n_{h+1} <= k 
    }

    unsigned long t_h = 0, t_h1 = 0;
    elem_t alpha = smallprimefield_convert_in(1, Pptr), 
        alpha1 = smallprimefield_convert_in(1, Pptr),
        m1 = smallprimefield_convert_in (-1, Pptr);
    for (polysize_t i = 1; i < h; ++i) {
        t_h += (A->n[i-1] - A->n[h])*(A->n[i]-A->n[h]);
        t_h1 += (A->n[i-1] - A->n[h+1])*(A->n[i]-A->n[h+1]);
        alpha = smallprimefield_mul(alpha, 
                    smallprimefield_exp (A->lc[i], A->n[i-1] - A->n[i+1], Pptr), 
                    Pptr); 
    }
    t_h1 += (A->n[h-1] - A->n[h+1])*(A->n[h]-A->n[h+1]);
    alpha1 = smallprimefield_mul(alpha, 
                smallprimefield_exp(A->lc[h], A->n[h-1] - A->n[h+1], Pptr), 
                Pptr); 

    alpha = t_h & 1 ? smallprimefield_mul (alpha, m1, Pptr) : alpha; 
    alpha1 = t_h1 & 1 ? smallprimefield_mul (alpha1, m1, Pptr) : alpha1;
    if (results_mode == 1) {
        // compute s_{k+1}, s_{k}
        sres = makePolysA_spX(2);
        sres->polys[0] = constPolynomial_spX (smallprimefield_mul(alpha, A->lc[h], Pptr), Pptr);
        sres->polys[1] = constPolynomial_spX (smallprimefield_mul(alpha1, A->lc[h+1], Pptr), Pptr);
        if (results != NULL) {
            *results = sres;
        } else {
            freePolysA_spX (sres);
        }
        if (degs != NULL) {
            degs[0] = A->n[h];
            degs[1] = A->n[h+1];
        }
        return;
    } 

    if (R->size < 2) {
        fprintf(stderr, "DUSP Error: bad input (R->size = 0)\n");
        exit(EXIT_FAILURE);
    }
    polysize_t i = 2;
    duspoly_t *r_k1, *r_k, *tmp = NULL;
    // compute r_k, r_k1 s.t. deg(r_k1) = n_{h} > k deg(r_k) = n_{h+1} <= k
    while(degPolynomial_spX(R->r[i-2]) != A->n[h+1] && degPolynomial_spX(R->r[i-1]) != A->n[h]) {
    // fprintf (stderr, "create rem_{k}, rem_{k+1} for i = %ld\n", i);
        if (R->size <= i) {
            // we have to compute upper remainders (R->r[i]) 
            // using list of quotients and the last two remainders (r_k, r_k1)
            // R[i] = R[i-1] Q[size-1] + R[i-2] for 2 <= i < R->size
            // fprintf(stderr, "[new] q[%ld-1] := ", i); printPolynomialOutForm_spX (Q->q[i-1], Pptr);
            // fprintf(stderr, "[new] r[%ld-1] := ", i); printPolynomialOutForm_spX (R->r[i-1], Pptr);
            // fprintf(stderr, "[new] r[%ld-2] := ", i); printPolynomialOutForm_spX (R->r[i-2], Pptr);
            mulPolynomialsInForm_spX (R->r[i-1], Q->q[Q->size-i+1], &tmp, Pptr);
            addPolynomialsInForm_spX_inp (&tmp, R->r[i-2], Pptr);
            R->r[i] = tmp; tmp = NULL;
            // freePolynomial_spX (&tmp); 
            R->size++;
        }
        i++;            
    }
    r_k = R->r[i-2];
    r_k1 = R->r[i-1];
    
    // compute S_{k+1}, S_{k}
    sres = makePolysA_spX(2);
    scalarMulPolynomialInForm_spX (r_k1, alpha, &(sres->polys[0]), Pptr);
    scalarMulPolynomialInForm_spX (r_k, alpha1, &(sres->polys[1]), Pptr);
    if (degs != NULL) {
        degs[0] = (sres->polys[0])->lt;
        degs[1] = (sres->polys[1])->lt;
    }
    if (results != NULL) {
        *results = sres;
    } else {
        freePolysA_spX (sres);
    }
}

int regularGCDbivarModularHGCDSubResInForm_withFFT_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, polysize_t *mdegs, int results_mode,  
                                           specAQR_spXY_t* specInfo, const Prime_ptr* Pptr)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    usfixn64* g_eval = convert_bivarElemPoly_to_usfixn64Poly(g, gpd, n); 
    usfixn64* f_eval = convert_bivarElemPoly_to_usfixn64Poly(f, fpd, n); 
    g_eval = bivarRowWiseDFT_spX(g_eval, gpd[1], n, K, e, w, P); 
    f_eval = bivarRowWiseDFT_spX(f_eval, fpd[1], n, K, e, w, P); 
    g_eval = transpose2D_spX(g_eval, gpd[1]+1, n); 
    f_eval = transpose2D_spX(f_eval, fpd[1]+1, n); 
    usfixn64** subres_arrays = (usfixn64**) malloc(n * sizeof(usfixn64*));
    polysize_t *chain_alloc_sizes = NULL, *max_chain_alloc_sizes = NULL;
    polysize_t subres_size, test_subres_size = 0;
    for(polysize_t i = 0; i < n; i++) {
        subres_size = 0;
        subres_arrays[i] = _RegularGCDModularSubres_from_usfixn64_spXY (&(g_eval[i * (gpd[1]+1)]), gpd[1], 
                                                                    &(f_eval[i * (fpd[1]+1)]), fpd[1], k_sub,
                                                                    &subres_size, &chain_alloc_sizes, mdegs, results_mode,
                                                                    specInfo, i, Pptr);
        if (test_subres_size != subres_size) {
            if (test_subres_size == 0 && i == 0) {
                test_subres_size = subres_size;
                max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
                for(polysize_t k = 0; k < test_subres_size; k++) {
                    max_chain_alloc_sizes[k] = chain_alloc_sizes[k];
                }
            }
            else {
                free(g_eval); free(f_eval);
                for (polysize_t j = 0; j < i; j++) { free(subres_arrays[j]);
                } free(subres_arrays);
                if (max_chain_alloc_sizes != NULL) {free(max_chain_alloc_sizes);}
                if (chain_alloc_sizes != NULL) {free(chain_alloc_sizes);}
                return 0;
            }
        }
        // check the goodness of omega
        for(polysize_t k = 0; k < test_subres_size &&
            max_chain_alloc_sizes[k] != chain_alloc_sizes[k]; k++) {
            free(g_eval); free(f_eval);
            for (polysize_t j = 0; j < i; j++) { free(subres_arrays[j]);
            } free(subres_arrays);
            if (max_chain_alloc_sizes != NULL) {free(max_chain_alloc_sizes);}
            if (chain_alloc_sizes != NULL) {free(chain_alloc_sizes);}
            return 0;
        }
        if (chain_alloc_sizes != NULL) {
            free(chain_alloc_sizes);
            chain_alloc_sizes = NULL;
        }
    }
    polysize_t alpha = gpd[0] + fpd[0];
    usfixn64** polys = _compute_DFTInv_modularSubres(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
	biSubresPr_t* s = makeBiSubresultantInForm_spX (test_subres_size);
	for (polysize_t i = 0; i < test_subres_size; i++) { 
		s->deg[i][0] = n-1; // deg_y
		s->deg[i][1] = max_chain_alloc_sizes[i]-1; // deg_x
		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1);
		s->coefs[i] = (elem_t*) malloc (s->size[i] * sizeof(elem_t));
        memcpy(s->coefs[i], polys[i], s->size[i]*sizeof(elem_t));
	}
    *subres = s;
    free(g_eval); free(f_eval);
    for (polysize_t j = 0; j < n; j++) { 
        free(subres_arrays[j]);
    }
    free(subres_arrays);
    for (polysize_t j = 0; j < test_subres_size; j++) {
        free(polys[j]);
    }
    free(polys);
    return 1;
}

int bivarModularHGCDSubResInForm_withFFT4_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, polysize_t r, elem_t* Omega, elem_t* OmegaInv, elem_t n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr)
{
    elem_t* g_eval = convert_bivarElemPoly_to_elemPoly(g, gpd, n); // (gpd[1]+1) * n
    elem_t* f_eval = convert_bivarElemPoly_to_elemPoly(f, fpd, n); // (fpd[1]+1) * n

    g_eval = bivarRowWiseDFT4_spX(g_eval, gpd[1], n, r, Omega, Pptr); // row-wise FFT(g_eval)
    f_eval = bivarRowWiseDFT4_spX(f_eval, fpd[1], n, r, Omega, Pptr); // row-wise FFT(f_eval)
    g_eval = transpose2D_elem_spX(g_eval, gpd[1]+1, n); // n * (gpd[1]+1)
    f_eval = transpose2D_elem_spX(f_eval, fpd[1]+1, n); // n * (fpd[1]+1)
    elem_t** subres_arrays = (elem_t**) malloc(n * sizeof(elem_t*));
    polysize_t *chain_alloc_sizes = NULL, *max_chain_alloc_sizes = NULL;
    polysize_t subres_size, test_subres_size = 0;
    for(polysize_t i = 0; i < n; i++) {
        subres_size = 0;
        subres_arrays[i] = _modularHGCDSubres_from_elems(&(g_eval[i * (gpd[1]+1)]), gpd[1],
                                                        &(f_eval[i * (fpd[1]+1)]), fpd[1], k_sub,
                                                        &subres_size, &chain_alloc_sizes, Pptr);
        // check the goodness of omega OR initialize test_subres_size and max_chain_alloc_sizes
        if (test_subres_size != subres_size) {
            if (test_subres_size == 0 && i == 0) {
                test_subres_size = subres_size;
                max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
                for(polysize_t k = 0; k < test_subres_size; k++) {
                    max_chain_alloc_sizes[k] = chain_alloc_sizes[k];
                }
            }
            else {
                free(g_eval); free(f_eval);
                for (polysize_t j = 0; j < i; j++) { free(subres_arrays[j]);
                } free(subres_arrays);
                if (max_chain_alloc_sizes != NULL) {free(max_chain_alloc_sizes);}
                if (chain_alloc_sizes != NULL) {free(chain_alloc_sizes);}
                return 0;
            }
        }
        // check the goodness of omega
        for(polysize_t k = 0; k < test_subres_size &&
            max_chain_alloc_sizes[k] != chain_alloc_sizes[k]; k++) {
            free(g_eval); free(f_eval);
            for (polysize_t j = 0; j < i; j++) { free(subres_arrays[j]);
            } free(subres_arrays);
            if (max_chain_alloc_sizes != NULL) {free(max_chain_alloc_sizes);}
            if (chain_alloc_sizes != NULL) {free(chain_alloc_sizes);}
            return 0;
        }
        if (chain_alloc_sizes != NULL) {
            free(chain_alloc_sizes);
            chain_alloc_sizes = NULL;
        }
    }
    polysize_t alpha = gpd[0] + fpd[0];
    elem_t** polys = _compute_DFT4Inv_modularSubres(subres_arrays, n, r, OmegaInv, n_inv, 
                        test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
	biSubresPr_t* s = makeBiSubresultantInForm_spX (test_subres_size);
	for (polysize_t i = 0; i < test_subres_size; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = n-1; // deg_y
		s->deg[i][1] = max_chain_alloc_sizes[i]-1; // deg_x
		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1);
		s->coefs[i] = (elem_t*) malloc (s->size[i] * sizeof(elem_t));
        memcpy(s->coefs[i], polys[i], s->size[i]*sizeof(elem_t));
	}
    *subres = s;
    free(g_eval); free(f_eval);
    for (polysize_t j = 0; j < n; j++) { free(subres_arrays[j]);
    }
    free(subres_arrays);
    for (polysize_t j = 0; j < test_subres_size; j++) {
        free(polys[j]);
    }
    free(polys);
    return 1;
}
