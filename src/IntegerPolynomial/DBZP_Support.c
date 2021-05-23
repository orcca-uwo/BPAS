#include "IntegerPolynomial/DBZP_Support.h"

void deepCopyPolynomial_DBZP_inp(const DBZP_t* a, DBZP_t** bb) {
	if (bb == NULL) {
		return;
	}

	if (a == NULL) {
		if (*bb != NULL) {
			freePolynomial_DBZP(*bb);
			*bb = NULL;
		}
		return;
	}

	DBZP_t* b = *bb;
	if (b == NULL) {
		b = makePolynomial_DBZP(a->ltX + 1);
		*bb = b;
	} else if (b->allocX <= a->ltX) {
		resizePolynomial_DBZP(b, a->ltX + 1);
	}

	int i;
	for (i = 0; i <= a->ltX; ++i) {
		//this will handle the case that b->coefsY[i] is NULL;
		deepCopyPolynomial_DUZP_inp(a->coefsY[i], b->coefsY + i);
	}
	for (i = a->ltX + 1; i <= b->ltX; ++i) {
		freePolynomial_DUZP(b->coefsY[i]);
		b->coefsY[i] = NULL;
	}

	b->ltX = a->ltX;
}



// /**
//  * Subresultant Chian Modular a prime
//  *
//  * @param g a bivariate polynomial (DBZP)
//  * @param gpd 2Dim partial degrees of polynomial g
//  * @param f a bivariate polynomial (DBZP) s.t. deg(g, x) > deg(f, x)
//  * @param fpd 2Dim partial degrees of polynomial f
//  * @param subres returned subresultant chain modular pr
//  * @param Pptr small prime pointer
//  */
// void subresultantModPr_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd, DBZP_t* subres, const Prime_ptr* Pptr)
// {
// 	// corner cases

// 	// set-up:
// 	// n2, n3 are total degrees of g, f
// 	int n, n2=0, n3=0;
// 	n=gpd[0]+1; // #polys + 1
// 	for (polysize_t i = 0; i <= gpd[1]; i++) {
// 		for (polysize_t j = 0; j <= gpd[0]; ++j) {
// 			if (g[n*i+j] != 0 && (i+j) > n3) {
// 				n2 = i+j;
// 			}
// 		}
// 	}
// 	n=fpd[0]+1; // #polys + 1
// 	for (polysize_t i = 0; i <= fpd[1]; i++) {
// 		for (polysize_t j = 0; j <= fpd[0]; ++j) {
// 			if (f[n*i+j] != 0 && (i+j) > n3) {
// 				n3 = i+j;
// 			}
// 		}
// 	}
// 	n2 = n2*n3 + 1; // double-check!
// 	n = fpd[1]*gpd[0] + fpd[0]*gpd[1] + 1; // size of the resultant based on the Sylv Matrix
// 	if (n2 < n) {
// 		n = n2;
// 	}
// 	if (n < 2) {
// 		n = 2;
// 	}

// 	n2 = (fpd[1] > 2) ? fpd[1] : 2; // number of subresultants to be computed
// 	int as = gpd[1] + 1, bs = fpd[1] + 1; // degree in x of f(t,x) and g(t,x) mod (t-t0)
// 	elem_t* a = (elem_t*) calloc (as*n, sizeof(elem_t));
// 	elem_t* b = (elem_t*) calloc (bs*n, sizeof(elem_t));
// 	elem_t* t = (elem_t*) calloc (n, sizeof(elem_t)); // evaluation points

// 	duspolysA_t** S = (duspolysA_t**) malloc (n*sizeof(duspolysA_t*)); // copies of subres over pr
// 	// TODO: can implement better...
// 	for (polysize_t i = 0; i < n; i++) {
// 		S[i] = makePolysA_spX (n2);
// 	}

// 	mpz_t m, e;
// 	mpz_init_set_si (m, Pptr->prime);
// 	// evaluate t and compute subres of main variable (x)
// 	// TODO: cilk_for, cilk_grainsize optimization needed!
// 	for (polysize_t k = 0; k < n; k++) {
// 		// work on evaluation point l*n + k + 1
// 		t[k] = k+1;
// 		polysize_t ad=0, bd=0, l=0;
// 		polysize_t at=gpd[0]+1, bt=fpd[0]+1;
// 		while (ad < gpd[1] || bd < fpd[1]) {
// 			t[k] = smallprimefield_add (t[k], smallprimefield_mul (l, n, Pptr), Pptr);
// 			// a:
// 			for (polysize_t i = 0; i <= gpd[1]; i++) {
// 				mpz_init_set_si (e, smallprimefield_convert_out (g[at*i+gpd[0]], Pptr));
// 				for (polysize_t j = gpd[0]-1; j > -1; --j) {
// 					mpz_mul_si (e, e, smallprimefield_convert_out (t[k], Pptr)); // e *= t[k]
// 					mpz_add_ui (e, e, smallprimefield_convert_out (g[at*i+j], Pptr)); // e += g[at*i+j]
// 				}
// 				mpz_mod (e, e, m); // e = e mod m(pr)
// 				a[k*as+i] = smallprimefield_convert_in (mpz_get_si (e) , Pptr);
// 				if (a[k*as+i]) { ad = i; } // uppdate a->size
// 				mpz_clear (e);
// 			}
// 			// b:
// 			for (polysize_t i = 0; i <= fpd[1]; i++) {
// 				mpz_init_set_si (e, smallprimefield_convert_out (f[bt*i+fpd[0]], Pptr));
// 				for (polysize_t j = fpd[0]-1; j > -1; --j) {
// 					mpz_mul_si (e, e, smallprimefield_convert_out (t[k], Pptr)); // e *= t[k]
// 					mpz_add_ui (e, e, smallprimefield_convert_out (f[bt*i+j], Pptr)); // e += f[bt*i+j]
// 				}
// 				mpz_mod (e, e, m); // e = e mod m(pr)
// 				b[k*bs+i] = smallprimefield_convert_in (mpz_get_si (e) , Pptr);
// 				if (b[k*bs+i]) { bd = i; } // uppdate b->size
// 				mpz_clear (e);
// 			}
// 			l++;
// 		}
// 		// Do sylvSubres to compute S[k] over Z_pr[y]
// 		duspoly_t* a_dusp = makePolynomial_spX (ad+1);
// 		duspoly_t* b_dusp = makePolynomial_spX (bd+1);
// 		a_dusp->elems = a[k*as];
// 		b_dusp->elems = b[k*as];
// 		a_dusp->lt = ad;
// 		b_dusp->lt = bd;
// 		int sz = 0;
// 		subresultantChainInForm_spX (a_dusp, b_dusp, &S[k], &sz, Pptr);
// 		if (a_dusp != NULL) { a_dusp->lt = 0; free (a_dusp); }
// 		if (b_dusp != NULL) { b_dusp->lt = 0; free (b_dusp); }
// 	}

// 	// Final Step
// 	// lagrange interpolation (with a 2D distributed-data-type)
// 	polysize_t k = gpd[0] + fpd[0];
// 	biSubresPr_t* s = makeBiSubresultantInForm_spX (n2);
// 	for (polysize_t i = 0; i < n2; i++) { // set-up returned-interpolated-subresultant-chain
// 		s->deg[i][0] = n-i*k-1;
// 		s->deg[i][1] = i;
// 		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1);
// 		s->coefs[i] = (elem_t*) calloc (s->size[i], sizeof(elem_t));
// 	}
// 	_interpolatePartialSubResInForm_spXY (S, &s, t, n, k, Pptr);

// 	free (a);
// 	free (b);
// 	free (t);
// 	for (polysize_t i = 0; i < n; i++) {
// 		freePolysA_spX (S[i]);
// 	}
// 	free (S);
// 	*subres = convertFromBiSubRes_DUZP_inp (s);
// 	return;
// }

biSubresZ_t* modularSetBiSubresultant_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd, mpz_t m, int startIdx, int primeSz)
{
	// first modular subres
	Prime_ptr  Pptr[1];
	*Pptr = prime64_ptr[startIdx];
	biSubresPr_t* S = NULL;
	mpz_t tc;
	mpz_init (tc);
	// Convert g to Z_p[0]
	polysize_t g_tdeg = (gpd[0]+1)*(gpd[1]+1);
	elem_t* gg = (elem_t*) malloc (g_tdeg*sizeof(elem_t));
	for (polysize_t i = 0; i < g_tdeg; i++) {
		mpz_mod_ui (tc, g[i], (unsigned long) Pptr->prime);
		gg[i] = smallprimefield_convert_in (mpz_get_si (tc), Pptr);
	}
	// Convert f to Z_p[0]
	polysize_t f_tdeg = (fpd[0]+1)*(fpd[1]+1);
	elem_t* ff = (elem_t*) malloc (f_tdeg*sizeof(elem_t));
	for (polysize_t i = 0; i < f_tdeg; i++) {
		mpz_mod_ui (tc, f[i], (unsigned long) Pptr->prime);
		ff[i] = smallprimefield_convert_in (mpz_get_si (tc), Pptr);
	}
	mpz_clear (tc);

	biSylvSubResultantInForm_spX (gg, gpd, ff, fpd, &S, Pptr);
	if (!S->n) {
		free (gg);
		free (ff);
		return NULL;
	}

	biSubresZ_t* a = makeBiSubresultant_DBZP (S->n);
	mpz_t halfm;
	mpz_init_set_si (halfm, (Pptr->prime-1)>>1);
	// construct a w.r.t pr[0]
	for (polysize_t j = 0; j < S->n; j++) {
		a->deg[j][0] = S->deg[j][0];
		a->deg[j][1] = S->deg[j][1];
		a->size[j]   = S->size[j];
		a->coefs[j] = (mpz_t*) malloc (a->size[j]*sizeof(mpz_t));
		for (polysize_t l = 0; l < a->size[j]; l++) {
			mpz_init_set_si (a->coefs[j][l], smallprimefield_convert_out (S->coefs[j][l], Pptr));
			if (mpz_cmp (a->coefs[j][l], halfm) > 0) {
				mpz_sub_ui (a->coefs[j][l], a->coefs[j][l], (unsigned long) Pptr->prime);
			}
		}
	}
	mpz_set_si (m, Pptr->prime);
	freeBiSubresultantInForm_spX (S); S=NULL;

	for (int i = 1; i < primeSz; i++) {
		// next modular subres
		*Pptr =  prime64_ptr[startIdx+i];
		mpz_init (tc);
		// Convert g to Z_p[i]
		for (polysize_t idx = 0; idx < g_tdeg; idx++) {
			mpz_mod_ui (tc, g[idx], (unsigned long) Pptr->prime);
			gg[idx] = smallprimefield_convert_in (mpz_get_si (tc), Pptr);
		}
		// Convert f to Z_p[i]
		for (polysize_t idx = 0; idx < f_tdeg; idx++) {
			mpz_mod_ui (tc, f[idx], (unsigned long) Pptr->prime);
			ff[idx] = smallprimefield_convert_in (mpz_get_si (tc), Pptr);
		}

		biSylvSubResultantInForm_spX (gg, gpd, ff, fpd, &S, Pptr);
		mpz_t t, x, e, o;
		mpz_inits (t, e, o);
		mpz_init_set_si (x, Pptr->prime);
		mpz_gcdext (t, e, o, x, m);
		mpz_mul (t, m, x);
		mpz_fdiv_q_2exp (halfm, t, 1ul);
		// CRT
		for (polysize_t j = 0; j < a->n; j++) {
			for (polysize_t l = 0; l < a->size[j]; l++) {
				mpz_ui_sub (tc, smallprimefield_convert_out (S->coefs[j][l], Pptr), a->coefs[j][l]);
				mpz_mul (tc, tc, o);
				mpz_mod (tc, tc, x);
				mpz_mul (tc, tc, m);
				mpz_add (a->coefs[j][l], a->coefs[j][l], tc);
				if (mpz_cmp (a->coefs[j][l], halfm) > 0) {
					mpz_sub (a->coefs[j][l], a->coefs[j][l], t);
				}
			}
		}
		mpz_set (m, t);
		mpz_clears (t, x, e, o, tc, NULL);
		freeBiSubresultantInForm_spX (S); S=NULL;
	}

	return a;
}

biSubresZ_t* modularBiSubresultant_DBZP (mpz_t* g, polysize_t* gpd, mpz_t* f, polysize_t* fpd)
{
	int primeIdx=0, k=1;
	Prime_ptr Pptr[1] = {prime64_ptr[primeIdx]}; // TODO: check bad-prime?

	//  the first biSubResultant call
	mpz_t m;
	mpz_init (m);
	biSubresZ_t* s = modularSetBiSubresultant_DBZP (g, gpd , f, fpd, m, primeIdx, k);

	biSubresZ_t* a = makeBiSubresultant_DBZP (s->n);
	// construct a w.r.t pr[0]
	for (polysize_t j = 0; j < s->n; j++) {
		a->deg[j][0] = s->deg[j][0];
		a->deg[j][1] = s->deg[j][1];
		a->size[j]   = s->size[j];
		a->coefs[j] = (mpz_t*) malloc (a->size[j]*sizeof(mpz_t));
		for (polysize_t l = 0; l < a->size[j]; l++) {
			mpz_init (a->coefs[j][l]);
		}
	}

	while (1) { // check stability at the end!
		primeIdx += k;
		mpz_t x, e, o, t, halft;
		mpz_inits (x, e, o, t, halft, NULL);

		biSubresZ_t* h = modularSetBiSubresultant_DBZP (g, gpd, f, gpd, x, primeIdx, k);
		mpz_gcdext (t, e, o, x, m);

		mpz_mul (t, x, m);
		mpz_tdiv_q_2exp (halft, t, 1); // double-check!
		// CRT
		for (polysize_t i = 0; i < a->n; i++) {
			for (polysize_t j = 0; j < a->size[i]; j++) {
				mpz_sub (a->coefs[i][j], h->coefs[i][j], s->coefs[i][j]);
				mpz_mul (a->coefs[i][j], a->coefs[i][j], o);
				mpz_mod (a->coefs[i][j], a->coefs[i][j], x);
				mpz_mul (a->coefs[i][j], a->coefs[i][j], m);
				mpz_add (a->coefs[i][j], a->coefs[i][j], s->coefs[i][j]);
				if (mpz_cmp(a->coefs[i][j], halft) > 0) {
					mpz_sub (a->coefs[i][j], a->coefs[i][j], t);
				}
			}
		}

		// Check if it's stable or not (s == a):
		int isDone = 0;
		if (s->n == a->n) {
			int isStable = 1;
			for (polysize_t j = 0; j < a->n; j++) {
				if (a->size[j] == s->size[j]) {
					for (polysize_t l = 0; l < a->size[j]; l++) {
						if (mpz_cmp (a->coefs[j][l], s->coefs[j][l])) { isStable = 0; break; }
					}
				} else { isStable = 0; break; }
				if (!isStable) { break; }
			}
			if (isStable) { isDone = 1; }
		}

		mpz_clears (x, e, o, halft, NULL);
		freeBiSubresultant_DBZP (h); h=NULL;
		if (isDone) {
			mpz_clear (t);
			break;
		}

		freeBiSubresultant_DBZP (s); s=NULL;
		s = deepCopyBiSubresultant_DBZP (a);
		mpz_set (m, t);
	}

	mpz_clear (m);
	freeBiSubresultant_DBZP (a);
	return s;
}