#include "ModularPolynomial/DUSP_Support_Factoring.h"


void computePolyToPowerPInv_spX (const duspoly_t* a, duspoly_t** ap, const Prime_ptr* Pptr) 
{

	if (ap == NULL) {
		return;
	}

	if (isZero_spX (a)) {
		freePolynomial_spX (&*ap);
		return;
	}

	polysize_t k = degPolynomial_spX (a)/Pptr->prime;
	duspoly_t* aa = makePolynomial_spX (k+1);
	aa->lt = k; // bound

	for (polysize_t i = 0; i < k+1; i++) {
		aa->elems[i] = a->elems[i*(Pptr->prime)];
	}

	normalize_spX (&aa);

	*ap = aa;
}


void computePolyToPowerPInv_spX_inp (duspoly_t** a, const Prime_ptr* Pptr) 
{
	if (a == NULL) {
		return;
	}

	if (isZero_spX (*a)) {
		return;
	}

	polysize_t k = degPolynomial_spX (*a)/Pptr->prime;
	duspoly_t* aa = *a;

	polysize_t i;
	for (i = 1; i < k+1; i++) {
		aa->elems[i] = aa->elems[i*(Pptr->prime)];
	}

	for (i = k+1; i < aa->alloc; i++){
		aa->elems[i] = 0;
	}

	normalize_spX (&aa);
}

void squareFreeFactorizationInForm_spX (const duspoly_t* a, factors_t** f, const Prime_ptr* Pptr) 
{
	if (f == NULL) {
		return;
	}

	if (isZero_spX (a)){
		*f = NULL;
		return;
	}

	elem_t lc_a = smallprimefield_convert_in (1, Pptr);

	if (degPolynomial_spX (a) == 0) {
		*f = initFactors_spX (1);
		(*f)->polys[0] = deepCopyPolynomial_spX (a);
		(*f)->exps[0] = 1;
		return;
	} else if (degPolynomial_spX (a) == 1) {
		
		*f = initFactors_spX (2);
		
		if (a->elems[0] != lc_a) {
			monicPolynomialInForm_spX (a, &(*f)->polys[1], &lc_a, Pptr);
			(*f)->polys[0] = constPolynomialInForm_spX (lc_a, Pptr);
			(*f)->exps[0] = 1;
		} else {
			(*f)->polys[0] = constPolynomialInForm_spX (lc_a, Pptr);
			(*f)->exps[0] = 1;
			(*f)->polys[1] = deepCopyPolynomial_spX (a);
			(*f)->exps[1] = 1;
		}

		return;
	}

	duspoly_t* aa = NULL;
	int isDCp = 1;

	if (a->elems[0] != lc_a) {
		monicPolynomialInForm_spX (a, &aa, &lc_a, Pptr);		
	} else {
		aa = deepCopyPolynomial_spX (a);
	}

	factors_t* ff = initFactors_spX (degPolynomial_spX (aa)+1);
	ff->polys[0] = constPolynomialInForm_spX (lc_a, Pptr);
	ff->exps[0] = 1;

	polysize_t i = 1; 
	polysize_t j = 1; 
	
	duspoly_t* g = NULL;
	duspoly_t* ap= NULL;
	duspoly_t* r = NULL;
	duspoly_t* q = NULL;

	duspoly_t* wg= NULL;
	duspoly_t* wq= NULL;

	duspoly_t* tmp=NULL;

	do {

		derivativePolyInForm_spX (aa, &ap, Pptr);
		GCDInForm_spX (aa, ap, &g, Pptr); //TODO: with HGCD-based GCD, use in-place 
		plainDivPolynomialsInForm_spX (aa, g, &r, &q, Pptr);

		freePolynomial_spX (&r);
		freePolynomial_spX (&ap);

		if (degPolynomial_spX (q) > 0) {
			do {

				GCDInForm_spX (q, g, &wg, Pptr);
				plainDivPolynomialsInForm_spX (q, wg, &r, &wq, Pptr);

				freePolynomial_spX (&r);

				if (degPolynomial_spX (wq) > 0) {
				
					if (ff->polys[i*j] == NULL) {
						ff->polys[i*j] = wq; wq = NULL;
						ff->exps[i*j] = i*j;

					} else {
						mulPolynomialsInForm_spX (ff->polys[i*j], wq, &tmp, Pptr);
						freePolynomial_spX (&ff->polys[i*j]);
						freePolynomial_spX (&wq);
						ff->polys[i*j] = tmp; tmp = NULL;
						ff->exps[i*j] = i*j;
					}
				} else {
					freePolynomial_spX (&wq);
				}

				if (degPolynomial_spX (wg) > 0) {
					
					plainDivPolynomialsInForm_spX (g, wg, &r, &tmp, Pptr);

					freePolynomial_spX (&g);
					freePolynomial_spX (&r);
					freePolynomial_spX (&q);

					g = tmp; tmp = NULL;
					q = wg;	wg = NULL;

					i++;
				} else {
					freePolynomial_spX (&wg);
					break;
				}

			} while (1);

			freePolynomial_spX (&q);
		}

		if (isDCp) {
			freePolynomial_spX (&aa);
		}

		if (degPolynomial_spX (g) > 0) {

			j = j*Pptr->prime;
			computePolyToPowerPInv_spX (g, &aa, Pptr);

			freePolynomial_spX (&g);

		} else {

			freePolynomial_spX (&g);
			*f = ff;
			return;
		}
	} while (1);
}

void squareFreeFactorizationInFormll_spX (const duspoly_t* a, factsll_t** F, const Prime_ptr* Pptr) 
{
	if (F == NULL) {
		return;
	}

	if (isZero_spX (a)){
		*F = NULL;
		return;
	}

	elem_t lc_a = smallprimefield_convert_in (1, Pptr);

	if (degPolynomial_spX (a) == 0) {
		*F = initFactsll_spX (0);
		(*F)->poly = deepCopyPolynomial_spX (a);
		return;
	} else if (degPolynomial_spX (a) == 1) {
		
		factsll_t* node1 = initFactsll_spX (1);

		if (leadingCoeffInForm_spX (a) != lc_a) {
			fprintf(stderr, "DUSP Error, In squareFreeFactorizationInFormll, input must be monic!\n");
			exit(1);

		} else {
			node1->poly = deepCopyPolynomial_spX (a);
			*F = node1;
		}

		return;
	}

	duspoly_t* aa;

	if (leadingCoeffInForm_spX (a) != lc_a) {
			fprintf(stderr, "DUSP Error, In squareFreeFactorizationInFormll, input must be monic!\n");
			exit(1);
	} else {
		aa = deepCopyPolynomial_spX (a);
		// isDCp 0;
	}

	factsll_t* node;
	factsll_t* head = NULL;
	factsll_t* tail;

	polysize_t i = 1; 
	polysize_t j = 1; 
	
	duspoly_t* g = NULL;
	duspoly_t* ap = NULL;
	duspoly_t* r = NULL;
	duspoly_t* q = NULL;

	duspoly_t* wg = NULL;
	duspoly_t* wq = NULL;

	duspoly_t* tmp = NULL;

	do {

		derivativePolyInForm_spX (aa, &ap, Pptr);
		GCDInForm_spX (aa, ap, &g, Pptr); 
		plainDivPolynomialsInForm_spX (aa, g, &r, &q, Pptr);

		freePolynomial_spX (&r);
		freePolynomial_spX (&ap);

		if (degPolynomial_spX (q) > 0) {
			do {

				GCDInForm_spX (q, g, &wg, Pptr);
				plainDivPolynomialsInForm_spX (q, wg, &r, &wq, Pptr);

				freePolynomial_spX (&r);

				if (degPolynomial_spX (wq) > 0) {
						
					node = initFactsll_spX (i*j);
					node->poly = wq; wq = NULL;

					if (head == NULL) {
						head = node;
						tail = node;
					} else {
						tail->next = node;
						tail = node;
					}

				} else {
					freePolynomial_spX (&wq);
				}

				if (degPolynomial_spX (wg) > 0) {
					
					plainDivPolynomialsInForm_spX (g, wg, &r, &tmp, Pptr);

					freePolynomial_spX (&g);
					freePolynomial_spX (&r);
					freePolynomial_spX (&q);

					g = tmp; tmp = NULL;
					q = wg;	wg = NULL;

					i++;
				} else {
					freePolynomial_spX (&wg);
					break;
				}

			} while (1);

			freePolynomial_spX (&q);
		}

		if (degPolynomial_spX (g) > 0) {

			j = j*Pptr->prime;
			freePolynomial_spX (&aa);

			computePolyToPowerPInv_spX (g, &aa, Pptr);

			freePolynomial_spX (&g);

		} else {

			freePolynomial_spX (&aa);
			freePolynomial_spX (&g);

			*F = head;
			return;
		}
	} while (1);
}

duspoly_t* expXModSubInForm_spX (polysize_t q, const duspoly_t* f, duspoly_t** h, const Prime_ptr* Pptr) 
{
	if (q == 0 || isZero_spX (f)) {
		return NULL;
	}

	duspoly_t* xpoly = XpolynomialInForm_spX ();
	duspoly_t* xpow = NULL;
	duspoly_t* res  = NULL;

	if (isZero_spX (*h)) {
		xpow = makePowerOfXPolynomialInForm_spX (q, Pptr);
		plainRemPolynomialsInForm_spX_inp (&xpow, f, Pptr);

		if (h != NULL) {
			*h = xpow;
		}
	} else {
		modExponentiatePolynomialInForm_spX (*h, q, f, &xpow, Pptr);

		freePolynomial_spX (h);
		if (h != NULL) {
			*h = xpow;
		}
	}

	subPolynomialsInForm_spX (xpow, xpoly, &res, Pptr);

	freePolynomial_spX (&xpoly);
	
	// xpow is h now!
	if (h == NULL) {
		freePolynomial_spX (&xpow);
	}
	
	return res;
}

void distinctDegFactorizationInFormll_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr)
{
	if (G == NULL){
		return;
	}

	if (isZero_spX (f)) {
		*G = NULL;
		return;
	}

	if (leadingCoeffInForm_spX (f) != smallprimefield_convert_in (1, Pptr)) {
		fprintf(stderr, "DUSP Error, In distinctDegFactorizationInFormll, input must be monic!\n");
		exit(1);
	}

	polysize_t s = 1;
	
	duspoly_t* t = NULL;
	duspoly_t* h = NULL;
	duspoly_t* f1 = NULL;

	duspoly_t* hmx = expXModSubInForm_spX (Pptr->prime, f, &h, Pptr);
	factsll_t* nodeG = initFactsll_spX (s);

	GCDInForm_spX (hmx, f, &(nodeG->poly), Pptr);
	plainDivPolynomialsInForm_spX (f, nodeG->poly, &t, &f1, Pptr);
	
	freePolynomial_spX (&hmx);
	freePolynomial_spX (&t);

	factsll_t* head = nodeG;
	factsll_t* tail = nodeG;

	duspoly_t* fs = NULL;

	while (f1 != NULL && f1->lt > 0) {
		s++;
		hmx = expXModSubInForm_spX (Pptr->prime, f, &h, Pptr);
		nodeG = initFactsll_spX (s);

		GCDInForm_spX (hmx, f1, &(nodeG->poly), Pptr);
		plainDivPolynomialsInForm_spX (f1, nodeG->poly, &t, &fs, Pptr);

		freePolynomial_spX (&hmx);
		freePolynomial_spX (&t);

		freePolynomial_spX (&f1);
		f1 = fs;

		tail->next = nodeG;
		tail = nodeG;
	}


	freePolynomial_spX (&f1);
	freePolynomial_spX (&h);

	*G = head;
}


void modPoweringForEDF_spX (const duspoly_t* a, polysize_t d, const duspoly_t* f, duspoly_t** mul_pows, const Prime_ptr* Pptr)
{

	if (mul_pows == NULL) {
		return;
	}

	if (d == 0) {
		*mul_pows = constPolynomial_spX (1, Pptr);
		return;
	}


	polysize_t p = Pptr->prime;
	polysize_t n = (p-1)>>1;

	if (n < 0) {
		fprintf(stderr, "DUSP Error: Bad prime in modPoweringForEDF.\n");
		exit(1);
	}

	if (d == 1) {
		modExponentiatePolynomialInForm_spX (a, n, f, mul_pows, Pptr);
		return;
	}

	// compute a**q, a**q**2, a**q**3, ..., a**q**(d-1) 
	// and multiply all together in muls

	duspoly_t* muls = NULL;
	duspoly_t* aqi  = NULL;
	duspoly_t* aqim = NULL;

	// compute aqim = a**p mod f 
	modExponentiatePolynomialInForm_spX (a, p, f, &muls, Pptr);
	
	if (d > 2) {
		aqim = deepCopyPolynomial_spX (muls);
	}

	mulPolynomialsInForm_spX_inp (&muls, a, Pptr);
	plainRemPolynomialsInForm_spX_inp (&muls, f, Pptr); // muls = (a * a**p) (mod f)


	for (polysize_t i = 2; i < d; i++) {
		// compute aqi = aqim**p = a**p**(i) mod f
		modExponentiatePolynomialInForm_spX (aqim, p, f, &aqi, Pptr);

		mulPolynomialsInForm_spX_inp (&muls, aqi, Pptr);
		plainRemPolynomialsInForm_spX_inp (&muls, f, Pptr);

		freePolynomial_spX (&aqim);
		aqim = aqi;
		aqi = NULL;
	}

	freePolynomial_spX (&aqim);

	modExponentiatePolynomialInForm_spX (muls, n, f, &aqi, Pptr);
	freePolynomial_spX (&muls);
	
	*mul_pows = aqi;
}

int equalDegSplittingInForm_spX (const duspoly_t* f, polysize_t d, duspoly_t** g, const Prime_ptr* Pptr)
{
	if (g == NULL) {
		return 0;
	}

	if (isZero_spX (f) || d == 0) {
		*g = NULL;
		return 1;
	}

	duspoly_t* a = randomPolynomialInForm_spX (degPolynomial_spX (f), Pptr);

	if (a == NULL ||  a->lt < 1) {
		freePolynomial_spX (&a);
		return 0;
	}

	duspoly_t* g1 = NULL;
	GCDInForm_spX (f, a, &g1, Pptr);

	if (g1 != NULL && g1->lt > 0) {
		*g = g1;

		freePolynomial_spX (&a);
		return 1;
	} else {
		freePolynomial_spX (&g1);
	}

	duspoly_t* b  = NULL;
	duspoly_t* bb = NULL;

	modPoweringForEDF_spX (a, d, f, &b, Pptr);

	duspoly_t* one = constPolynomialInForm_spX (smallprimefield_convert_in (1, Pptr), Pptr);

	subPolynomialsInForm_spX (b, one, &bb, Pptr);
	
	freePolynomial_spX (&a);
	freePolynomial_spX (&b);
	freePolynomial_spX (&one);

	GCDInForm_spX (bb, f, &g1, Pptr);
	
	freePolynomial_spX (&bb); 

	if (g1 != NULL && g1->lt > 0 && !isEqual_spX (f, g1)) {
		*g = g1;

		return 1;
	} else {
		freePolynomial_spX (&g1);
	}

	return 0;

}

void equalDegFactorizationInFormll_spX (const duspoly_t* f, polysize_t d, factsll_t** G, const Prime_ptr* Pptr)
{
	if (G == NULL) {
		return;
	}

	if (isZero_spX (f)) {
		return;
	}

	if (degPolynomial_spX (f) <= d) {
		factsll_t* node = initFactsll_spX (0);
		node->poly = deepCopyPolynomial_spX (f);

		if (*G == NULL) {
			*G = node;
		} else {
			factsll_t* cur = *G;
			while (cur->next != NULL) {
				cur = cur->next;
			}
			cur->next = node;
		}

		return;
	}

	duspoly_t* g = NULL;
	int isFact = 0;
	
	while (!isFact) {
		isFact = equalDegSplittingInForm_spX (f, d, &g, Pptr);
	}

	duspoly_t* r = NULL;
	duspoly_t* q = NULL;

	plainDivPolynomialsInForm_spX (f, g, &r, &q, Pptr);
	freePolynomial_spX (&r);

	equalDegFactorizationInFormll_spX (g, d, G, Pptr);
	equalDegFactorizationInFormll_spX (q, d, G, Pptr);

	freePolynomial_spX (&g);
	freePolynomial_spX (&q);
}

void equalDegFactorsInFormll_spX (factsll_t* D, factsll_t** G, const Prime_ptr* Pptr)
{
	if (G == NULL) {
		return;
	}

	if (isZeroFactsll_spX (D)) {
		*G = NULL;
		return;
	}

	factsll_t* Dcur = D;

	factsll_t* EDFo = NULL;
	factsll_t* EDFcur = NULL;

	factsll_t* node = NULL;
	factsll_t* head = NULL;
	factsll_t* tail = NULL;

	while (Dcur != NULL) {

		equalDegFactorizationInFormll_spX (Dcur->poly, Dcur->exp, &EDFo, Pptr);

		if (!isZeroFactsll_spX (EDFo)) {
			EDFcur = EDFo;

			while (EDFcur != NULL) {
				
				if (!isConstant_spX (EDFcur->poly)) {
					node = initFactsll_spX (1);
					node->poly = EDFcur->poly;

					if (head == NULL) {
						head = node;
						tail = node;
					} else {
						tail->next = node;
						tail = node; 
					}
				} else {
					freePolynomial_spX (&(EDFcur->poly));
				}

				EDFcur = EDFcur->next;								
			}

			free (EDFo);
			EDFo = NULL;
		}

		Dcur = Dcur->next;
	}

	*G = head;
}

void modFactorizationInFormll_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr)
{

	if (G == NULL) {
		return;
	}	
	
#if TIMER_SUBPROGRAMS
    timer_id id;
    timer_time t;
    double sum1 = -1;

    timer_id id2;
    timer_time t2;
    double sum2 = -1;
#endif 

	if (isZero_spX (f)) {
		*G = NULL;
		return;
	}

	if (leadingCoeffInForm_spX (f) != smallprimefield_convert_in (1, Pptr)) {
		fprintf(stderr, "DUSP Error, In modFactorizationInFormll, input must be monic!\n");
		exit(1);
	}


	polysize_t i = 0;
	polysize_t e;

	int isDiv = 1;

	factsll_t* cur;

	duspoly_t* q = NULL;
	duspoly_t* r = NULL;

	duspoly_t* hi  = NULL;
	duspoly_t* hmx = NULL;
	
	duspoly_t* fi  = deepCopyPolynomial_spX (f);

	factsll_t* Gg = NULL;
	duspoly_t* g  = NULL;

	factsll_t* head = NULL;
	factsll_t* node;
	factsll_t* tail = NULL;


#if TIMER_SUBPROGRAMS
	    id2 = start_timer ();
#endif 

	while (fi != NULL && fi->lt > 0) {

#if DEBUG_PRINT_LINES
		fprintf(stderr, "f[%lu] := \n", i);
		printPolynomial_spX (fi, Pptr);

#endif 

		i += 1;

#if TIMER_SUBPROGRAMS
	    id = start_timer ();
#endif 

		hmx = expXModSubInForm_spX (Pptr->prime, f, &hi, Pptr);		

#if TIMER_SUBPROGRAMS
	    t = elapsed_time (&id);
        sum1 = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [expXModSub = %f ]\n", sum1);
#endif 
#if TIMER_SUBPROGRAMS
	    id = start_timer ();
#endif 

		GCDInForm_spX (hmx, fi, &g, Pptr);

#if TIMER_SUBPROGRAMS
	    t = elapsed_time (&id);
        sum1 = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [GCD = %f ]\n", sum1);
#endif 
		freePolynomial_spX (&hmx);

		if (g == NULL || g->lt < 1) {
			
			freePolynomial_spX (&g);

			continue;
		}

#if  DEBUG_PRINT_LINES
		fprintf(stderr, "gcd(x^{q^%lu}-x, f_{%lu-1}) = ", i, i);
		printPolynomial_spX (g, Pptr);

		if (Gg != NULL) {
			fprintf(stderr, "Gg must be NULL but, \n");
			printFactsllOutForm_spX (Gg, Pptr);
			fprintf(stderr, "\n");
		}
#endif
#if TIMER_SUBPROGRAMS
	    id = start_timer ();
#endif 

		
		equalDegFactorizationInFormll_spX (g, i, &Gg, Pptr);


#if TIMER_SUBPROGRAMS
	    t = elapsed_time (&id);
        sum1 = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [EDF = %f ]\n", sum1);
#endif 
#if  DEBUG_PRINT_LINES
		fprintf(stderr, "Gg[%lu] = \n",i);
		printFactsllOutForm_spX (Gg, Pptr);
		fprintf(stderr, "\n");
		int j = 0; // TEST
#endif

		freePolynomial_spX (&g);

		cur = Gg;
		while (cur != NULL) {
			r = deepCopyPolynomial_spX (cur->poly);

#if  DEBUG_PRINT_LINES
			fprintf(stderr, "cur->poly[%d] := ",j);
			printPolynomialOutForm_spX (cur->poly, Pptr);
			j++;
#endif

#if TIMER_SUBPROGRAMS
	id = start_timer ();
#endif
			
			isDiv = isDividablePolysQuoInForm_spX (fi, r, &q, Pptr);

#if TIMER_SUBPROGRAMS
	    t = elapsed_time (&id);
        sum1 = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [Div = %f ]\n", sum1);
#endif 

			freePolynomial_spX (&fi);
			fi = q; q = NULL;
			e = 1;

#if TIMER_SUBPROGRAMS
	id = start_timer ();
#endif		
		
			isDiv = isDividablePolysQuoInForm_spX (fi, r, &q, Pptr);

#if TIMER_SUBPROGRAMS
	    t = elapsed_time (&id);
        sum1 = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [Div = %f ]\n", sum1);
#endif

			while (isDiv) {
				e += 1;

				freePolynomial_spX (&fi);
				fi = q; q = NULL;
				
#if TIMER_SUBPROGRAMS
	id = start_timer ();
#endif		

				isDiv = isDividablePolysQuoInForm_spX (fi, r, &q, Pptr);
			
#if TIMER_SUBPROGRAMS
	    t = elapsed_time (&id);
        sum1 = (t.tv_sec + ((double)t.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [Div = %f ]\n", sum1);
#endif

			}

			freePolynomial_spX (&q);

			node = initFactsll_spX (e);
			node->poly = r; r = NULL;

			if (head == NULL) {
				head = node;
				tail = node;
			} else {
				tail->next = node;
				tail = node;
			}

			cur = cur->next;
		}

		freeFactsll_spX (&Gg);
		Gg = NULL;
	}

#if TIMER_SUBPROGRAMS
	    t2 = elapsed_time (&id2);
        sum2 = (t2.tv_sec + ((double)t2.tv_usec / CLOCKS_PER_SEC));
	    fprintf(stderr, "[time] [whileLoop = %f ]\n", sum2);
#endif 	

	freePolynomial_spX (&hi);
	freePolynomial_spX (&fi);


	*G = head;
}

int modFactorVerificationInFormll_spX (const duspoly_t* f, factsll_t* G, const Prime_ptr* Pptr) 
{
	if (isZero_spX (f)) {
		if (isZeroFactsll_spX (G)) 
			return 1;
		else 
			return 0;
	} 

	factsll_t* cur = G;
	duspoly_t* ff = constPolynomial_spX (1, Pptr);
	duspoly_t* tmp = NULL;
	duspoly_t* tmp2 = NULL;

	while (cur != NULL) {
		if (cur->poly != NULL) {
			if (cur->exp == 0) {
				fprintf(stderr, "DUSP Error, In modFactorVerificationInFormll: G couldn't have polynomial with idx=0 as a factor!\n");
				exit (1);
			}

			exponentiatePolynomialInForm_spX (cur->poly, cur->exp, &tmp, Pptr);
			mulPolynomialsInForm_spX (ff, tmp, &tmp2, Pptr);
			
			freePolynomial_spX (&tmp);
			freePolynomial_spX (&ff);
			
			ff = tmp2;
			tmp2 = NULL;
		} else {
			fprintf(stderr, "DUSP Error, In modFactorVerificationInFormll: G couldn't have zero polynomial as a factor!\n");
			exit (1);
		}

		cur = cur->next;
	}

	if (isEqual_spX (f, ff)) {
		return 1;
	} else {
		fprintf(stderr, "f := ");
		printPolynomialOutForm_spX (f, Pptr);

		fprintf(stderr, "prod(facts) := ");
		printPolynomialOutForm_spX (ff, Pptr);

		fprintf(stderr, "G = \n");
		purePrintFactsllOutForm_spX (G, Pptr);

		fprintf(stderr, "DUSP Error, In modFactorVerificationInFormll is FAILED!\n");
		exit(1);

	}
}

factors_t* convertToFactorsList_spX (factsll_t* F, const Prime_ptr* Pptr) 
{

	if (isZeroFactsll_spX (F)) {
		return NULL;
	}

	polysize_t sz = lenFactsll_spX (F);
	factors_t* Fl = initFactors_spX (sz);

	factsll_t* cur = F;
	for (polysize_t i = 0; i < sz; i++) {
		Fl->polys[i] = deepCopyPolynomial_spX (cur->poly);
		Fl->exps[i] = cur->exp;

		cur = cur->next;
	}

	return Fl;
}

