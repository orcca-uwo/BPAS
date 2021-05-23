

#include "ModularPolynomial/DUSP_Parallel.hpp"

// static int totalUnivarImages = 0;
// static float totalInterpTime = 0;
// static float totalUnivarTime= 0;
// static float totalCRTTime = 0;
static int THREAD_ELEMENTS_BOUNDARY = 128;
static int THREAD_UNIVAR_MIN_SIZE = 128	;

extern "C" {
	void _lagrangeBasisInForm_spX (const elem_t* t, const polysize_t n, elem_t* a, const Prime_ptr* Pptr);
	elem_t _evalDivPolynomials_spX (elem_t* a, elem_t u, elem_t* q, polysize_t n, const Prime_ptr* Pptr);
}

//evaluate and compute univar subresultant images for for k1 <= k < k2
//k1: inclusive starting index
//k2: exclusive ending index
//n: total number of images to be computed, needed for 1-D array access of 2D data
//n2: a bound on the number of subresultants in the chain
//PPtr: the prime pointer for the current prime field
//gpd: partial degrees of the first input poly
//fpd: partial degrees of the second input poly
//t: evaluation points, used later for interp
//g: the first input poly, stored as a dense 2D array
//f: the second input poly stored as a dense 2D array
//S: An array of size n of duspolysA_t* for storing the results
void DUSP::evalAndGetUnivarImage_SubRes_spX(int k1, int k2, int n, int n2, const Prime_ptr* Pptr, const polysize_t* gpd, const polysize_t* fpd, const elem_t* g, const elem_t* f, elem_t* t, duspolysA_t** S) {
	// work on evaluation point l*n + k + 1
	polysize_t ad, bd, l;
    elem_t e, tt;

    duspoly_t* a_dusp = makePolynomial_spX(gpd[1]+1);
    duspoly_t* b_dusp = makePolynomial_spX(fpd[1]+1);
    elem_t* a = a_dusp->elems;
    elem_t* b = b_dusp->elems;

    elem_t tmp;
	for (polysize_t k = k1; k < k2; k++) {
		ad = 0; bd = 0; l = 0;
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
	        }
			l++;
		}
		// Do sylvSubres to compute S[k] over Z_pr[y]
		// duspoly_t* a_dusp = makePolynomial_spX (ad+1);
		// duspoly_t* b_dusp = makePolynomial_spX (bd+1);
		// // a_dusp->elems = (elem_t*) malloc (a_dusp->alloc * sizeof(elem_t)); // &a[k*as];
	 //    a_dusp->elems = memcpy (a_dusp->elems, &a[k*as], a_dusp->alloc * sizeof(elem_t));
		// // b_dusp->elems = (elem_t*) malloc (b_dusp->alloc * sizeof(elem_t)); // &b[k*bs];
	 //    b_dusp->elems = memcpy (b_dusp->elems, &b[k*bs], b_dusp->alloc * sizeof(elem_t));
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
	    // ++totalUnivarImages;
	    if (sz > n2) {
	        fprintf (stderr, "DUSP ERROR: size of returned Subres is larger than expectation!\n");
	        exit (1);
	    } else {
	        S[k] = makePolysA_spX (n2);
	        for (polysize_t i = 0; i < sz; i++) {
	            if (isZero_spX (SS->polys[i])) {
	                S[k]->polys[i] = NULL;
	            } else {
	                S[k]->polys[i] = SS->polys[i];
	                SS->polys[i] = NULL;
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
	    // for (polysize_t tt = 0; tt < n2; tt++) {  
	    // fprintf (stderr, "[k=%ld] Subres[%ld] := ", k, tt);
	    // printPolynomialOutForm_spX (S[k]->polys[tt], Pptr); }  

	    // Don't free all the polys in the duspolysA_t since we copied into S[k].
	    // freePolysA_spX (SS);
	    free(SS->polys);
	    free(SS);
	}
    freePolynomial_spX (&a_dusp);
    freePolynomial_spX (&b_dusp);
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
void interpolatePartialSubResInForm_parallel_spXY (elem_t* g, elem_t* q, polysize_t k1, polysize_t k2,
	duspolysA_t** pSubres, biSubresPr_t* XYsubres, elem_t* t, polysize_t n, polysize_t alpha, const Prime_ptr* Pptr,
	const std::vector<ExecutorThreadPool::threadID>& workers)
{
    // for each row of the subresultant chain
    // TODO: cilk_for, cilk_grainsize optimization needed!
    for (polysize_t k = k1; k < k2; k++) {
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
            for (polysize_t s = 0; s <= tmpZ->lt; s++) { // TODO: check the uppor bound ??
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
void DUSP::biSylvSubResultantInForm_parallel_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, biSubresPr_t** subres, const Prime_ptr* Pptr,
	const std::vector<ExecutorThreadPool::threadID>& workers)
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
    // fprintf (stderr, "[biSbResPr] n = %d\n", n);  
	n2 = (fpd[1] > 2) ? fpd[1] : 2;
    // fprintf (stderr, "n2' = %d , n = %d\n", n2, n);
    // n2: number of subresultants to be computed based on x in Z[x>y]
	int as = gpd[1] + 1, bs = fpd[1] + 1;
    // degree in x of f(t,x) and g(t,x) mod (t-t0)
    // in evaluation of y at t.
	// elem_t* a = (elem_t*) calloc (as*n, sizeof(elem_t)); // to store evaluated polys in Z(t_i)[x] s.t. i \in {0, ..., n-1}
	// elem_t* b = (elem_t*) calloc (bs*n, sizeof(elem_t)); // to store evaluated polys in Z(t_i)[x] s.t. i \in {0, ..., n-1}
	elem_t* t = (elem_t*) calloc (n, sizeof(elem_t)); // to store evaluation points
	duspolysA_t** S = (duspolysA_t**) malloc (n*sizeof(duspolysA_t*)); // copies of subres over pr
    polysize_t true_sz = 0;
    polysize_t at=gpd[0]+1, bt=fpd[0]+1;
    // evaluate t and compute subres of the second variable (x)
	// TODO: cilk_for, cilk_grainsize optimization needed!

	// for (polysize_t k = 0; k < n; k++) {
	// 	// evalAndGetUnivarImage_SubRes_spX(int k1, int k2, int n, int n2, Prime_ptr* Pptr, polysize_t* gpd, polysize_t* fpd, elem_t* t, elem_t* g, elem_t* f, duspolysA_t** S);
	// }

    // unsigned long long start;
	// float univarTime, interpTime;
	// _startTimer(&start);
    int nworkers = workers.size();
    int k1 = 0, k2 = 0;
    for (int k = 0; k < nworkers; k++) {
    	//get number of images to compute, rounding k2 down
    	// to nearest multiple of MIN_IMAGES_PER_THREAD
    	k2 = (n / (nworkers+1)); //+1 becuase current thread does work to
    	k2 = (k2 / THREAD_ELEMENTS_BOUNDARY) * THREAD_ELEMENTS_BOUNDARY;

    	//add offset
    	k2 += k1;

    	std::function<void()> task = std::bind(DUSP::evalAndGetUnivarImage_SubRes_spX,
			k1, k2, n, n2, Pptr, gpd, fpd, g, f, t, S);
    	ExecutorThreadPool::getThreadPool().executeTask(workers[k], task);
    	// fprintf(stderr, "TEST: worker %d evaling [%d, %d)\n", k, k1, k2);
   		k1 = k2;
    }
    //do the last one on the current thread
    //if nthreads does not evenly divide n, do the extra images on the current thread
    k2 = n;
	// fprintf(stderr, "TEST: launch local worker: %d\n", workerIdx);
	// fprintf(stderr, "TEST: worker %d evaling [%d, %d)\n", nworkers, k1, k2);
	evalAndGetUnivarImage_SubRes_spX(k1, k2, n, n2, Pptr, gpd, fpd, g, f, t, S);

	//sync threads
	ExecutorThreadPool::getThreadPool().waitForThreads(workers);
	// _stopTimer(&start, &univarTime);
	// totalUnivarTime += univarTime;

    // for (polysize_t i = 0; i < n; i++) {
    //     fprintf (stderr, "t[%ld] = %lld\n", i, smallprimefield_convert_out (t[i], Pptr)); }
	// Final Step
	// lagrange interpolation (with a 2D distributed-data-type)
	polysize_t degsum = gpd[0] + fpd[0];
	biSubresPr_t* s = makeBiSubresultantInForm_spX (n2);
	for (polysize_t i = 0; i < n2; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = (n-i*degsum-1); // deg_y // it's always >= 0
		s->deg[i][1] = i;  // deg_x
		s->size[i] = (n-i*degsum-1) > -1 ? (s->deg[i][0]+1)*(s->deg[i][1]+1) : 0;
		s->coefs[i] = (elem_t*) calloc (s->size[i], sizeof(elem_t));
	}
    // for (int k = 0; k < n; k++) {
    //     for (polysize_t tt = 0; tt < n2; tt++) {  
    //         fprintf (stderr, "[k=%ld] Subres[%ld] := ", k, tt);
    //         printPolynomialOutForm_spX (S[k]->polys[tt], Pptr); }  
    //         fprintf (stderr, "\n\n"); }  



    // To store all images:
    elem_t* gwork = (elem_t*) calloc (n2 * (n+1), sizeof (elem_t));
    elem_t* qwork = (elem_t*) calloc (n2 * n, sizeof (elem_t));
    //n is deg in y we want to interpolate
    //sum_i=0^n2 (n - i*degsum)
    // polysize_t delta=n-chainIndex*degsum, offset=n*chainIndex;

	// _startTimer(&start);
    polysize_t amtWork = (n*n2) - (degsum*(n2-1)*n2/2);
    amtWork /= (nworkers+1); //estimate of work per thread;
    k1 = 0, k2 = 0;
    polysize_t curWork = 0;
    int workerIdx = 0;
    for (int k = 0; k < n2; ++k) {
    	curWork += n-k*degsum;
    	if (curWork > amtWork) {
    		k2 = k + 1; // exclusive upper bound
	    	std::function<void()> task = std::bind(interpolatePartialSubResInForm_parallel_spXY,
				gwork, qwork, k1, k2, S, s, t, n, degsum, Pptr, workers);
	    	ExecutorThreadPool::getThreadPool().executeTask(workers[workerIdx], task);
	    	// fprintf(stderr, "worker %d doing [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
	    	k1 = k2;
	    	++workerIdx;
			curWork = 0;
    	}
    }
    k2 = n2;
	// fprintf(stderr, "worker %d doing [%d,%d), curWork: %d\n", nworkers, k1, k2, curWork);
	interpolatePartialSubResInForm_parallel_spXY(gwork, qwork, k1, k2, S, s, t, n, degsum, Pptr, workers);

	//sync threads
	ExecutorThreadPool::getThreadPool().waitForThreads(workers);
	// _stopTimer(&start, &interpTime);
	// totalInterpTime += interpTime;

    free (qwork);
    free (gwork);

	free (t);
	for (polysize_t i = 0; i < n; i++) {
		freePolysA_spX (S[i]);
	}
	free (S);
    *subres = s;
}

void biModularSubresultantCRT_parallel(int k1, int k2, biSubresZ_t* subres_work, const biSubresPr_t* modSubres, int* subResEq, mpz_t tmpZ, const mpz_t mpz_pr, const mpz_t s, const mpz_t halfm, const mpz_t newm, const Prime_ptr* Pptr) {
	int isSubresEq = 1;
	elem_t coef_out;
	for (polysize_t k = k1; k < k2; ++k) {
		if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
			continue;
		}
		polysize_t lt = modSubres->size[k];
		subres_work->deg[k][0] = modSubres->deg[k][0];
		subres_work->deg[k][1] = modSubres->deg[k][1];
		subres_work->size[k] = modSubres->size[k]; // update lt
		for (polysize_t i = 0; i < lt; ++i) {
			mpz_set(tmpZ, subres_work->coefs[k][i]);
			coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
			// optimized-CRT (Garner's algorithm)

			mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
			mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
			mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
			mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
			mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
			if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
				mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
			}

			isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);

		}
	}
	*subResEq = isSubresEq;
}


int DUSP::biModularSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime, int maxThreads)
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In modularBiSubresultantChain_DBZP, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}

	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	biSubresZ_t* subres_work = NULL;
	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, NULL);


	int primeIdx = 0;
	Prime_ptr Pptr[1];

	AltArrZ_t* tmp;
	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);
	elem_t* modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t* modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime = (unsigned long) Pptr->prime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;


	polysize_t chainSize = (b_degs[1] > 2) ? b_degs[1] : 2;

	ExecutorThreadPool& pool = ExecutorThreadPool::getThreadPool();

	int nthreads; // No. Threads
	if (maxThreads < 0) {
		int nUnivarImages = b_degs[1]*a_degs[0] + b_degs[0]*a_degs[1] + 1;
		// fprintf(stderr, "TEST: Univar Images Needed: %d\n", nUnivarImages );

		nthreads = ((nUnivarImages) / THREAD_UNIVAR_MIN_SIZE) - 1; //cur thread counts as one
		nthreads = (chainSize-1) < nthreads ? (chainSize-1) : nthreads;
		// fprintf(stderr, "TEST: nthreads: %d\n", nthreads );
		// fprintf(stderr, "TEST: chain size: %d\n", chainSize );

		infinityNorm_AAZ(a, s);
		infinityNorm_AAZ(b, t);
		if (mpz_cmp(s, t) < 0) {
			mpz_swap(s,t); //make s max
		}
		int nprimes = mpz_size(s) * a_degs[1] * 3/2; //an estimate for nthreads
		if (nprimes < 125) {
			//do nothing
		} else if (nprimes < 150) {
			if (nthreads < 1) {
				nthreads = 1;
			}
		} else if (nprimes < 175) {
			if (nthreads < 2) {
				nthreads = 2;
			}
		} else if (nprimes < 250) {
			if (nthreads < 3) {
				nthreads = 3;
			}
		} else if (nprimes < 300) {
			if (nthreads < 4) {
				nthreads = 4;
			}
		} else {
			nthreads = ExecutorThreadPool::maxThreads;
		}
		nthreads = nthreads < 0 ? 0 : nthreads;
	} else {
		nthreads = maxThreads;
	}

	std::vector<ExecutorThreadPool::threadID> workers;
	nthreads = pool.obtainThreads(nthreads, workers);
	// fprintf(stderr, "TEST final nthreads: %d\n", nthreads);
	// float elapsed;
	// unsigned long long start;
	// float CRTTime = 0.0f;

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Parallel- nworkers=%d\n", (int) workers.size());
#endif

	int subResEq[nthreads+1];
	mpz_t tmpZ[nthreads+1];
	for (int i = 0; i < nthreads+1; ++i) {
		mpz_init(tmpZ[i]);
		subResEq[i] = 0;
	}

	// int isNULL = 0;
	for (; primeIdx < n_prime64_ptr; ++primeIdx) {
		// fprintf(stderr, "In biModularSubresultantChainZ, main-for-loop primeIdx = %d\n", primeIdx);  
		// Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		*Pptr = prime64_ptr[primeIdx];
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) || 
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			continue;
		}

		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		biSylvSubResultantInForm_parallel_spX (modA, a_degs, modB, b_degs, &modSubres, Pptr, workers);

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		int isSubresEq = 1;
		// _startTimer(&start);
		if (primeIdx < 256) {
			biModularSubresultantCRT_parallel(0, modSubres->n, subres_work, modSubres, &subResEq[0], 
											tmpZ[0], mpz_pr, s, halfm, newm, Pptr);
			isSubresEq = subResEq[0];
		} else {
			int amtWork = 0;
			for (polysize_t k = 0; k < modSubres->n; ++k) {
				amtWork += modSubres->size[k];
			}
			amtWork /= (nthreads+1);
			int k1 = 0, k2 = 0;
			int workerIdx = 0;
			int curWork = 0;
			for (int k = 0; k < modSubres->n; ++k) {
				curWork += modSubres->size[k];
				if (curWork > amtWork) {
					k2 = k;
					if (!subResEq[workerIdx]) {
						std::function<void()> task = std::bind(biModularSubresultantCRT_parallel, k1, k2,
							subres_work, modSubres, &subResEq[workerIdx], tmpZ[workerIdx], mpz_pr, s, halfm, newm, Pptr);
						ExecutorThreadPool::getThreadPool().executeTask(workers[workerIdx], task);
						// fprintf(stderr, "TEST: worker %d CRT [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
					} else {
						// fprintf(stderr, "TEST: worker %d CRT [%d,%d), skipping work\n", workerIdx, k1, k2, curWork);
					}
					k1 = k;
					++workerIdx;
					curWork = 0;
				}
			}
			k2 = modSubres->n;
			if (!subResEq[workerIdx]) {
				// fprintf(stderr, "TEST: main t %d CRT [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
				biModularSubresultantCRT_parallel(k1, k2,
					subres_work, modSubres, &subResEq[workerIdx], tmpZ[workerIdx], mpz_pr, s, halfm, newm, Pptr);
				++workerIdx;
			} else {
				// fprintf(stderr, "TEST: worker %d CRT [%d,%d), skipping work\n", workerIdx, k1, k2, curWork);
			}

			//sync threads
			ExecutorThreadPool::getThreadPool().waitForThreads(workers);
			// _stopTimer(&start, &CRTTime);
			// totalCRTTime += CRTTime;

			//Notice here we use workerIdx and not nthreads+1
			//Not all threads may have done work in the CRT loop
			for (int i = 0; i < workerIdx; ++i) {
				isSubresEq = isSubresEq && subResEq[i];
			}
		}

		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			if (no_prime != NULL) {*no_prime = primeIdx+1;}
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (a);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			freeBiSubresultant_DBZP (subres_work);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free (modA); free (modB);
			for (int i = 0; i < nthreads+1; ++i) {
				mpz_clear(tmpZ[i]);
			}

			// reverse the order s.t. tail->head:
			if (sz<2) {
				*Subres = head;
			} else {
				AltArrsZ_t* curr = head, *prev = NULL, *next = NULL;
				while (curr != NULL){
					next = curr->next;
					curr->next = prev;
					prev = curr;
					curr = next;
				}
				*Subres = prev;
			}
			*chain_size = sz;

			pool.returnThreads(workers);

			// fprintf(stderr, "TEST: total elapsed: %g\n", totalElapsed );
			// fprintf(stderr, "TEST: total univar: %g\n", totalUnivarTime );
			// fprintf(stderr, "TEST: total interp: %g\n", totalInterpTime);
			// fprintf(stderr, "TEST: total CRT: %g\n", totalCRTTime);
			// fprintf(stderr, "TEST: TOTAL NVAR IMAGES: %d\n", totalUnivarImages );
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
	for (int i = 0; i < nthreads+1; ++i) {
		mpz_clear(tmpZ[i]);
	}
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free (modA); free (modB);
	freeBiSubresultant_DBZP (subres_work);
	// fprintf (stderr, "In biModularSubresultantChain_DBZP, all primes failed\n");  
	*Subres = NULL;
	*chain_size = 0;

	pool.returnThreads(workers);

	return 0;
}

void modularSubres_from_usfixn64_parallel (int thr_k, int k1, int k2, usfixn64* g, polysize_t gpdx, usfixn64* f, polysize_t fpdx,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, const Prime_ptr* Pptr, 
											usfixn64** subres_array_out)
{

	duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
	duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
	polysize_t** chain_alloc_size_array = NULL;
	polysize_t* subres_size_array = NULL;
	if((k2-k1)) {
		chain_alloc_size_array = (polysize_t**) malloc((k2-k1)*sizeof(polysize_t*));
		subres_size_array = (polysize_t*) malloc((k2-k1)*sizeof(polysize_t));
	}
	// fprintf(stderr, "TEST: threadIdx = %d, k1 = %d, k2= %d\n", thr_k, k1, k2);
	for(int k = k1, idx = 0; k < k2; k++, idx++) {
		// fprintf(stderr, "TEST: k = %d\n", k);
		usfixn64* g_k = &(g[k * (gpdx+1)]);
		usfixn64* f_k = &(f[k * (fpdx+1)]);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { 
			a_dusp->elems[i] = (elem_t) g_k[i]; 
		} 
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { 
			b_dusp->elems[i] = (elem_t) f_k[i]; 
		}
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);
        // fprintf(stderr, "Par-TEST: a_dusp = ");  
        // printPolynomial_spX(a_dusp, Pptr);  
        // fprintf(stderr, "Par-TEST: b_dusp = ");   
        // printPolynomial_spX(b_dusp, Pptr);  
        polysize_t sz = 0;
        duspolysA_t* S = NULL;
        subresultantChainInForm_spX (a_dusp, b_dusp, &S, &sz, Pptr);
        // fprintf(stderr, "TEST: brown subres_sz = %ld\n", sz);  
        // for (polysize_t i = 0; i < sz; i++) {  
        //     fprintf(stderr, "TEST: brown subres[%ld] = ", i);  
        //     printPolynomial_spX(S->polys[i], Pptr);  
        // }  

        // brown corner case:
        if (S->polys[0] == NULL) {
            S->polys[0] = zeroPolynomial_spX();
        }
        // polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        chain_alloc_size_array[idx] =  (polysize_t*) malloc (sz*sizeof(polysize_t));
		polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
            if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            chain_alloc_size_array[idx][i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        usfixn64* subres_array = (usfixn64*) calloc(sz * total_alloc_size, sizeof(usfixn64));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), chain_alloc_size_array[idx][i] * sizeof(usfixn64));
            start += chain_alloc_size_array[idx][i];
        }
        // *chain_alloc_sizes = lcl_chain_alloc_sizes;
        // *subres_size = sz;
		subres_size_array[idx] = sz;
		freePolysA_spX(S);
		
		subres_array_out[k] = subres_array;
        
	}
	
	freePolynomial_spX(&a_dusp); // TODO: we can do better
	freePolynomial_spX(&b_dusp);
	
	if ((k2-k1)) {
		
		polysize_t test_subres_size_array = subres_size_array[0]; 
		for(int j = 1; j < k2-k1 && (test_subres_size_array != subres_size_array[j]); j++) {
			test_subres_size_array = -1;
			break;
		}

		subres_size[thr_k] = test_subres_size_array;
		chain_alloc_sizes[thr_k] = chain_alloc_size_array[0];


		free(subres_size_array);
		for(int i = 1; i < k2-k1; i++) {
			free(chain_alloc_size_array[i]);
		}
		free(chain_alloc_size_array);
	} else {
		subres_size[thr_k] = -1;
		chain_alloc_sizes[thr_k] = NULL;
	}
}

void bivarRowWiseDFT_k1_to_k2_spX (int k1, int k2, usfixn64* g, polysize_t n,
                        int K, int e, usfixn64 w, montgomery_triple P)
{
	for(int i = k1; i < k2; i++) {
		DFT_general_small_elements (&(g[i*n]), K, e, &P, w);
	}
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
usfixn64* bivarRowWiseDFT_spX_parallel (usfixn64* g, polysize_t gpd_x, polysize_t n,
                        int K, int e, usfixn64 w, montgomery_triple P,
						const std::vector<ExecutorThreadPool::threadID>& workers)
{

    //x=DFT_n(x)
	if (gpd_x < 32 && n < 257) {
		for(polysize_t i = 0; i <= gpd_x; i++) {
			DFT_general_small_elements (&(g[i*n]), K, e, &P, w);
		}		
	} else {
		int nworkers = workers.size();
		int k1 = 0, k2 = 0;
		for(int k = 0; k < nworkers; k++) {
			k2 = (gpd_x+1) / (nworkers+1);
			k2 = (k2/THREAD_ELEMENTS_BOUNDARY) * THREAD_ELEMENTS_BOUNDARY;
			k2 += k1;
			std::function<void()> task = std::bind(bivarRowWiseDFT_k1_to_k2_spX, k1, k2, g, n, K, e, w, P);
			ExecutorThreadPool::getThreadPool().executeTask(workers[k], task);
			k1 = k2;
		}

		k2 = gpd_x+1;
		bivarRowWiseDFT_k1_to_k2_spX (k1, k2, g, n, K, e, w, P);

		ExecutorThreadPool::getThreadPool().waitForThreads(workers);

	}
    return g;
}

void bivarRowWiseDFTInv_k1_to_k2_spX (int k1, int k2, usfixn64* g, polysize_t n,
                        int K, int e, usfixn64 w_inv, usfixn64 n_inv, montgomery_triple P)
{
	for(int i = k1; i < k2; i++) {
		DFT_general_small_elements (&(g[i*n]), K, e, &P, w_inv);
		for(polysize_t j = 0; j < n; j++) {
			g[i*n+j] = mult_ab_mod_p (g[i*n+j], n_inv, &P); 
		}
	}
}


/**
 * Row-wise DFT^{-1} to apply on Bivariate Polynomial Z[x > y] (in-place)
 * @note Assume w_inv and n_inv are already convertIn_GLOBAL_ptr
 */
usfixn64* bivarRowWiseDFTInv_spX_parallel (usfixn64* g, polysize_t gpd_x, polysize_t n,
                        int K, int e, usfixn64 w_inv, usfixn64 n_inv, montgomery_triple P,
						const std::vector<ExecutorThreadPool::threadID>& workers)
{
	if (gpd_x < 32 && n < 257) {
		//x=DFT_n^{-1}(x)
		for(polysize_t i = 0; i <= gpd_x; i++) {
			DFT_general_small_elements (&(g[i*n]), K, e, &P, w_inv);
			for(polysize_t j = 0; j < n; j++) {
				g[i*n+j] = mult_ab_mod_p (g[i*n+j], n_inv, &P); 
			}
		}
	} else {
		int nworkers = workers.size();
		int k1 = 0, k2 = 0;
		for (int k = 0; k < nworkers; k++) {
			k2 = (gpd_x+1) / (nworkers+1);
			k2 = (k2/THREAD_ELEMENTS_BOUNDARY) * THREAD_ELEMENTS_BOUNDARY;
			k2 += k1;
			std::function<void()> task = std::bind(bivarRowWiseDFTInv_k1_to_k2_spX, k1, k2, g, n, K, e, w_inv, n_inv, P);
			ExecutorThreadPool::getThreadPool().executeTask(workers[k], task);
			k1 = k2;
		}

		k2 = gpd_x+1;
		bivarRowWiseDFTInv_k1_to_k2_spX(k1, k2, g, n, K, e, w_inv, n_inv, P);

		ExecutorThreadPool::getThreadPool().waitForThreads(workers);
	}
    return g;
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
usfixn64** compute_DFTInv_modularSubres_parallel (usfixn64** SS, polysize_t n, int K, int e, usfixn64 w_inv, usfixn64 n_inv,
                                        polysize_t subres_size, polysize_t alpha, polysize_t* max_chain_alloc_sizes,
                                        const Prime_ptr* Pptr, const std::vector<ExecutorThreadPool::threadID>& workers)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    // evaluate_subres_n[i] = n * max_chain_alloc_sizes[i] for i = 0, ..., subres_size-1
    usfixn64** evaluated_subres_n = _compute_DFTInv_modularSubres_inputs(SS, n, subres_size, alpha, max_chain_alloc_sizes);
    for(polysize_t i = 0; i < subres_size; i++) {
       evaluated_subres_n[i] = transpose2D_spX(evaluated_subres_n[i], n, max_chain_alloc_sizes[i]); // max_chain_alloc_sizes[i] * n
       evaluated_subres_n[i] = bivarRowWiseDFTInv_spX_parallel(evaluated_subres_n[i], max_chain_alloc_sizes[i]-1, n, K, e, w_inv, n_inv, P, workers);
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
int DUSP::bivarModularSubResInForm_withFFT_parallel_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr,
										   const std::vector<ExecutorThreadPool::threadID>& workers)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    usfixn64* g_eval = convert_bivarElemPoly_to_usfixn64Poly(g, gpd, n); // (gpd[1]+1) * n
    usfixn64* f_eval = convert_bivarElemPoly_to_usfixn64Poly(f, fpd, n); // (fpd[1]+1) * n
    // fprintf(stderr, "\nTEST: (g, f) convert_bivarElemPoly_to_usfixn64Poly\n");  
    g_eval = bivarRowWiseDFT_spX_parallel(g_eval, gpd[1], n, K, e, w, P, workers); // row-wise FFT(g_eval)
    f_eval = bivarRowWiseDFT_spX_parallel(f_eval, fpd[1], n, K, e, w, P, workers); // row-wise FFT(f_eval)
    // fprintf(stderr, "TEST: (g, f) bivarRowWiseDFT_spX\n");  
    g_eval = transpose2D_spX(g_eval, gpd[1]+1, n); // n * (gpd[1]+1)
    f_eval = transpose2D_spX(f_eval, fpd[1]+1, n); // n * (fpd[1]+1)
    // fprintf(stderr, "TEST: (g, f) transpose2D_spX\n");  
    usfixn64** subres_arrays = (usfixn64**) malloc(n * sizeof(usfixn64*));

    // cilk-for
	int nworkers = workers.size();
	// fprintf(stderr, "TEST: n = %ld, nworkers = %d\n", n, nworkers);
    polysize_t** chain_alloc_sizes = (polysize_t**) malloc((nworkers+1) * sizeof(polysize_t*)); 
	polysize_t* max_chain_alloc_sizes = NULL;
    polysize_t* subres_size = (polysize_t*) calloc((nworkers+1), sizeof(polysize_t)); 
	polysize_t test_subres_size = 0;

	int k1 = 0, k2 = 0;
	for (int k = 0; k < nworkers; k++) {
		//get number of images to compute, rounding k2 down
    	// to nearest multiple of MIN_IMAGES_PER_THREAD
    	k2 = (n / (nworkers+1)); //+1 becuase current thread does work to
    	k2 = (k2 / THREAD_ELEMENTS_BOUNDARY) * THREAD_ELEMENTS_BOUNDARY;

		// add offset
		k2 += k1;

    	std::function<void()> task = std::bind(modularSubres_from_usfixn64_parallel,
			k, k1, k2, g_eval, gpd[1], f_eval, fpd[1], subres_size, chain_alloc_sizes, Pptr, subres_arrays);
    	ExecutorThreadPool::getThreadPool().executeTask(workers[k], task);

		k1 = k2;
	}
    //do the last one on the current thread
    //if nthreads does not evenly divide n, do the extra images on the current thread
	k2 = n;

	// fprintf(stderr, "TEST: launch local worker: %d\n", workerIdx);
	// fprintf(stderr, "TEST: worker %d evaling [%d, %d)\n", nworkers, k1, k2);

	modularSubres_from_usfixn64_parallel (nworkers, k1, k2, g_eval, gpd[1], f_eval, fpd[1], 
										subres_size, chain_alloc_sizes, Pptr, subres_arrays);

	//sync threads
	ExecutorThreadPool::getThreadPool().waitForThreads(workers);

	for (int i = 0; i <= nworkers && (chain_alloc_sizes[i] == NULL); i++) {
		fprintf(stderr, "DUSP Parallel Error: works are not devided equally\n");
		exit(1);
	}

	test_subres_size = subres_size[0];
	int isGoodOmega = 1;
	if(test_subres_size == -1) {
		isGoodOmega = 0;
	} else {
		for (int i = 1; i <= nworkers && (test_subres_size != subres_size[i]); i++) {
			isGoodOmega = 0;
			break;
		}
	}
	
	if(!isGoodOmega) {
		// fprintf(stderr, "TEST: (g, f) test_subres_size (%ld) != subres_size (%ld)\n", test_subres_size, subres_size);
		free(g_eval);
		free(f_eval);
		for (polysize_t j = 0; j < n; j++) {
			free(subres_arrays[j]);
		}
		free(subres_arrays);
		if (chain_alloc_sizes != NULL) {
			for(int j = 0; j <= nworkers && (chain_alloc_sizes[j] != NULL); j++) {
				free(chain_alloc_sizes[j]);
			}
			free(chain_alloc_sizes);
		}
		free(subres_size);

		// pool.returnThreads(workers);

		return 0;		
	}

	// TODO: 
	max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
	for(polysize_t j = 0; j < test_subres_size; j++) {
		max_chain_alloc_sizes[j] = chain_alloc_sizes[0][j];
	}
	for(int j = 0; j <= nworkers; j++) {
		free(chain_alloc_sizes[j]);
	}
	free(chain_alloc_sizes);


    polysize_t alpha = gpd[0] + fpd[0];
    usfixn64** polys = compute_DFTInv_modularSubres_parallel(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr, workers);
    // fprintf(stderr, "TEST: subres_size = %ld\n", test_subres_size);  
	biSubresPr_t* s = makeBiSubresultantInForm_spX (test_subres_size);
	for (polysize_t i = 0; i < test_subres_size; i++) { // set-up returned-interpolated-subresultant-chain
		s->deg[i][0] = n-1; // (n-i*alpha-1); // deg_y
		s->deg[i][1] = max_chain_alloc_sizes[i]-1; // deg_x
		s->size[i] = (s->deg[i][0]+1)*(s->deg[i][1]+1); // (n-i*alpha-1) > -1 ? (s->deg[i][0]+1)*(s->deg[i][1]+1) : 0;
		s->coefs[i] = (elem_t*) malloc (s->size[i] * sizeof(elem_t));
        memcpy(s->coefs[i], polys[i], s->size[i]*sizeof(elem_t));
        // for(polysize_t j = 0; j < s->size[i]; j++) {  ing ...
        //     s->coefs[i][j] = (elem_t) polys[i][j];
        // }
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

int DUSP::biModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, AltArrsZ_t** Subres, int* chain_size, int* no_prime, int maxThreads)
{

	// const char* sym[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};  
	// fprintf (stderr, "In biModularFFTSubresultantChainZ a := ");
	// printPoly_AAZ (stderr, a, sym, a->nvar);
	// fprintf (stderr, "\nIn biModularFFTSubresultantChainZ b := ");
	// printPoly_AAZ (stderr, b, sym, b->nvar);
	// fprintf(stderr, "\n");

	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In biModularFFTSubresultantChainZ, expected only bivariate polynomials.\n");
		exit (1);
	}

	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	biSubresZ_t* subres_work = NULL;
	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, NULL);

	AltArrZ_t* tmp;
	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	biSubresPr_t* modSubres = NULL;
	// biSubresZ_t* subres_work_prev = NULL;
	unsigned long uprime, coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389; 
	usfixn64 p_u64, w, w_inv, n_inv;
    montgomery_triple P;
	const int max_no_failed_omega = 3;
	int isGoodOmega, isSubresEq, no_failed_omega = 0;
	int K, e; 
	polysize_t n = _compute_N_gt_By(a_degs, b_degs, &K, &e); 
	polysize_t chainSize = (b_degs[1] > 2) ? b_degs[1] : 2;
 
	ExecutorThreadPool& pool = ExecutorThreadPool::getThreadPool(); // start the ThreadPool

	int nthreads; // No. Threads
	if (maxThreads < 0) {
		int nUnivarImages = b_degs[1]*a_degs[0] + b_degs[0]*a_degs[1] + 1;
		// fprintf(stderr, "TEST: Univar Images Needed: %d\n", nUnivarImages );

		nthreads = ((nUnivarImages) / THREAD_UNIVAR_MIN_SIZE) - 1; //cur thread counts as one
		nthreads = (chainSize-1) < nthreads ? (chainSize-1) : nthreads;
		// fprintf(stderr, "TEST: nthreads: %d\n", nthreads );
		// fprintf(stderr, "TEST: chain size: %d\n", chainSize );

		infinityNorm_AAZ(a, s);
		infinityNorm_AAZ(b, t);
		if (mpz_cmp(s, t) < 0) {
			mpz_swap(s,t); //make s max
		}
		int nprimes = mpz_size(s) * a_degs[1] * 3/2; //an estimate for nthreads
		if (nprimes < 125) {
			//do nothing
		} else if (nprimes < 150) {
			if (nthreads < 1) {
				nthreads = 1;
			}
		} else if (nprimes < 175) {
			if (nthreads < 2) {
				nthreads = 2;
			}
		} else if (nprimes < 250) {
			if (nthreads < 3) {
				nthreads = 3;
			}
		} else if (nprimes < 300) {
			if (nthreads < 4) {
				nthreads = 4;
			}
		} else {
			nthreads = ExecutorThreadPool::maxThreads;
		}
		nthreads = nthreads < 0 ? 0 : nthreads;
	} else {
		nthreads = maxThreads;
	}

	std::vector<ExecutorThreadPool::threadID> workers;
	nthreads = pool.obtainThreads(nthreads, workers);
	// fprintf(stderr, "TEST: final nthreads: %d\n", nthreads);
	// float elapsed;
	// unsigned long long start;
	// float CRTTime = 0.0f;

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Parallel- bivarPS-FFT nworkers=%d\n", (int) workers.size());
#endif


	int subResEq[nthreads+1];
	mpz_t tmpZ[nthreads+1];
	for (int i = 0; i < nthreads+1; ++i) {
		mpz_init(tmpZ[i]);
		subResEq[i] = 0;
	}

	// fprintf (stderr, "start testing primes...\n"); //TEST
	// int prime_upperbound = no_prime != 0 ? no_prime : n_prime64_ptr; // TODO: don't need this... just for testing :)
	for (; primeIdx < prime_table_size; ++primeIdx) { // primeIdx < n_prime64_ptr
		// *Pptr = prime64_ptr[primeIdx];
		// @note fourier_primes_u64_table is not in-order to use directly in CRT 
		// @note subset_fourier_primes_u64_table is a in-order subset of fourier_primes_u64_table 
		// @note if CRT fails for large inputs, check the prime list at first!
		// TODO: update fourier_primes_u64_table s.t. FFT-routines in GFPF won't break! 
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) || 
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			
			free(Pptr);
			continue;
		}

		p_u64 = (unsigned long long) Pptr->prime;
		init_montgomery_triple (&P, p_u64);
		compute_nth_root_of_unity_for_small_prime (p_u64, &w, n);
		
		// @note we check it now in bivarModularSubResInForm_withFFT_spX
		// @note it's not efficient when we check it here (in top-level). 

		gmp_inv_mod_p_u64 (&w_inv, w, p_u64);
	    gmp_inv_mod_p_u64 (&n_inv, n, p_u64);

		convertIn_GLOBAL_ptr(&w, &P);
	    convertIn_GLOBAL_ptr(&w_inv, &P);
    	convertIn_GLOBAL_ptr(&n_inv, &P);

		// if good prime and omega:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], (unsigned long) Pptr->prime), Pptr);
		}
		// Convert b to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], (unsigned long) Pptr->prime), Pptr);
		}

		// call subresultant mod this prime
		isGoodOmega = bivarModularSubResInForm_withFFT_parallel_spX(modA, a_degs, modB, b_degs, n, K, e, w, w_inv, n_inv, &modSubres, Pptr, workers);

		if(!isGoodOmega) { 
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				free(Pptr);
				continue; 
			} else {
				// free aZ, bZ
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				for (int i = 0; i < nthreads+1; ++i) {
					mpz_clear(tmpZ[i]);
				}
				free(Pptr);
				return 0;
			}
		} 

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		int isSubresEq = 1;
		// _startTimer(&start);
		if (primeIdx < 50) {
			biModularSubresultantCRT_parallel(0, modSubres->n, subres_work, modSubres, &subResEq[0], 
											tmpZ[0], mpz_pr, s, halfm, newm, Pptr);
			isSubresEq = subResEq[0];
		} else {
			int amtWork = 0;
			for (polysize_t k = 0; k < modSubres->n; ++k) {
				amtWork += modSubres->size[k];
			}
			amtWork /= (nthreads+1);
			int k1 = 0, k2 = 0;
			int workerIdx = 0;
			int curWork = 0;
			for (int k = 0; k < modSubres->n; ++k) {
				curWork += modSubres->size[k];
				if (curWork > amtWork) {
					k2 = k;
					if (!subResEq[workerIdx]) {
						std::function<void()> task = std::bind(biModularSubresultantCRT_parallel, k1, k2,
							subres_work, modSubres, &subResEq[workerIdx], tmpZ[workerIdx], mpz_pr, s, halfm, newm, Pptr);
						ExecutorThreadPool::getThreadPool().executeTask(workers[workerIdx], task);
						// fprintf(stderr, "TEST: worker %d CRT [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
					} else {
						// fprintf(stderr, "TEST: worker %d CRT [%d,%d), skipping work\n", workerIdx, k1, k2, curWork);
					}
					k1 = k;
					++workerIdx;
					curWork = 0;
				}
			}
			k2 = modSubres->n;
			if (!subResEq[workerIdx]) {
				// fprintf(stderr, "TEST: main t %d CRT [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
				biModularSubresultantCRT_parallel(k1, k2,
					subres_work, modSubres, &subResEq[workerIdx], tmpZ[workerIdx], mpz_pr, s, halfm, newm, Pptr);
				++workerIdx;
			} else {
				// fprintf(stderr, "TEST: worker %d CRT [%d,%d), skipping work\n", workerIdx, k1, k2, curWork);
			}
			//sync threads
			ExecutorThreadPool::getThreadPool().waitForThreads(workers);
			// _stopTimer(&start, &CRTTime);
			// totalCRTTime += CRTTime;

			//Notice here we use workerIdx and not nthreads+1
			//Not all threads may have done work in the CRT loop
			for (int i = 0; i < workerIdx; ++i) {
				isSubresEq = isSubresEq && subResEq[i];
			}
		}

		free(Pptr);
		if (modSubres != NULL) { 
			freeBiSubresultantInForm_spX (modSubres); 
		}
		modSubres = NULL;  


		// *********** Serial CRT ************//
		// // subres_work_prev ?= subres_work
		// if (subres_work_prev == NULL) {
		// 	subres_work_prev = makeBiSubresultant_DBZP (subres_work->n);
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		subres_work_prev->deg[j][0] = subres_work->deg[j][0];
		// 		subres_work_prev->deg[j][1] = subres_work->deg[j][1];
		// 		subres_work_prev->size[j]   = subres_work->size[j];
		// 		subres_work_prev->coefs[j] = (mpz_t*) malloc (subres_work_prev->size[j]*sizeof(mpz_t));
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_init_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// } else {
		// 	if (subres_work->n != subres_work_prev->n) {
		// 		fprintf (stderr, "subres_work->n != subres_work_prev\n");
		// 		exit (1);
		// 	}
		// 	for (polysize_t j = 0; j < subres_work->n; j++) {
		// 		if (subres_work_prev->size[j] != subres_work->size[j]) {
		// 			fprintf (stderr, "subres_work_prev->size[%ld] != subres_work->size[%ld]\n", j, j);
		// 			exit (1);
		// 		}
		// 		for (polysize_t l = 0; l < subres_work_prev->size[j]; l++) {
		// 			mpz_set (subres_work_prev->coefs[j][l], subres_work->coefs[j][l]);
		// 		}
		// 	}
		// }

		// for (polysize_t k = 0; k < modSubres->n; k++) {
		// 	if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
		// 		continue;
		// 	}
		// 	polysize_t lt = modSubres->size[k];
		// 	subres_work->deg[k][0] = modSubres->deg[k][0];
		// 	subres_work->deg[k][1] = modSubres->deg[k][1];
		// 	subres_work->size[k] = modSubres->size[k]; // update lt
		// 	for (polysize_t i = 0; i < lt; ++i) {
		// 		coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
		// 		// optimized-CRT (Garner's algorithm)
		// 		mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
		// 		mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
		// 			mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		}
		// 	}
		// }

		// free(Pptr);
		// if (modSubres != NULL) { freeBiSubresultantInForm_spX (modSubres); }
		// modSubres = NULL; // TEST
		// // subres_work == subres_work_prev ?
		// isSubresEq = isBiSubresultantEquals_DBZP (subres_work, subres_work_prev);
		////////////////////////////////

		if (isSubresEq) {
			if (no_prime != NULL) {*no_prime = primeIdx+1;} 
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (a);
			node->next = NULL;
			tail->next = node;
			tail = node;
			sz++;
			freeBiSubresultant_DBZP (subres_work);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			for (int i = 0; i < nthreads+1; ++i) {
				mpz_clear(tmpZ[i]);
			}
			// reverse the order s.t. tail->head:
			if (sz<2) {
				*Subres = head;
			} else {
				AltArrsZ_t* curr = head, *prev = NULL, *next = NULL;
				while (curr != NULL){
					next = curr->next;
					curr->next = prev;
					prev = curr;
					curr = next;
				}
				*Subres = prev;
			}
			*chain_size = sz;

			pool.returnThreads(workers);
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
	for (int i = 0; i < nthreads+1; ++i) {
		mpz_clear(tmpZ[i]);
	}
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	freeBiSubresultant_DBZP (subres_work);
	// fprintf (stderr, "In biModularFFTSubresultantChain_DBZP, all primes failed\n");  
	*Subres = NULL;
	*chain_size = 0;
	
	pool.returnThreads(workers);
	
	return 0;
}

void modularHGCDSubres_from_usfixn64_parallel (int thr_k, int k1, int k2, usfixn64* g, polysize_t gpdx, usfixn64* f, polysize_t fpdx, polysize_t sub_k,
                                            polysize_t* subres_size, polysize_t** chain_alloc_sizes, const Prime_ptr* Pptr, 
											usfixn64** subres_array_out)
{

	duspoly_t* a_dusp = makePolynomial_spX (gpdx+1);
	duspoly_t* b_dusp = makePolynomial_spX (fpdx+1);
	polysize_t** chain_alloc_size_array = NULL;
	polysize_t* subres_size_array = NULL;
	if((k2-k1)) {
		chain_alloc_size_array = (polysize_t**) malloc((k2-k1)*sizeof(polysize_t*));
		subres_size_array = (polysize_t*) malloc((k2-k1)*sizeof(polysize_t));
	}
	// fprintf(stderr, "TEST: threadIdx = %d, k1 = %d, k2= %d\n", thr_k, k1, k2);
	for(int k = k1, idx = 0; k < k2; k++, idx++) {
		// fprintf(stderr, "TEST: k = %d\n", k);
		usfixn64* g_k = &(g[k * (gpdx+1)]);
		usfixn64* f_k = &(f[k * (fpdx+1)]);
        for(polysize_t i = 0; i < a_dusp->alloc; i++) { 
			a_dusp->elems[i] = (elem_t) g_k[i]; 
		} 
        for(polysize_t i = 0; i < b_dusp->alloc; i++) { 
			b_dusp->elems[i] = (elem_t) f_k[i]; 
		}
		a_dusp->lt = gpdx;
		b_dusp->lt = fpdx;
        normalize_spX(&a_dusp);
        normalize_spX(&b_dusp);
        // fprintf(stderr, "Par-TEST: a_dusp = ");  
        // printPolynomial_spX(a_dusp, Pptr);  
        // fprintf(stderr, "Par-TEST: b_dusp = ");   
        // printPolynomial_spX(b_dusp, Pptr);  
        polysize_t sz = 0;
        duspolysA_t* S = NULL;
        // subresultantChainInForm_spX (a_dusp, b_dusp, &S, &sz, Pptr);
		hgcdSubResultantInFormA_spX(a_dusp, b_dusp, sub_k, &S, &sz, Pptr);
        // fprintf(stderr, "TEST: brown subres_sz = %ld\n", sz);  
        // for (polysize_t i = 0; i < sz; i++) {  
        //     fprintf(stderr, "TEST: brown subres[%ld] = ", i);  
        //     printPolynomial_spX(S->polys[i], Pptr);  
        // }  

        // brown corner case:
        if (S->polys[0] == NULL) {
            S->polys[0] = zeroPolynomial_spX();
        }
        // polysize_t* lcl_chain_alloc_sizes = (polysize_t*) malloc (sz*sizeof(polysize_t));
        chain_alloc_size_array[idx] =  (polysize_t*) malloc (sz*sizeof(polysize_t));
		polysize_t total_alloc_size = 0;
        for(polysize_t i = 0; i < sz; i++) {
            if (S->polys[i] == NULL) { S->polys[i] = zeroPolynomial_spX(); }
            chain_alloc_size_array[idx][i] = S->polys[i]->alloc;
            total_alloc_size += S->polys[i]->alloc;
        }
        usfixn64* subres_array = (usfixn64*) calloc(sz * total_alloc_size, sizeof(usfixn64));
        polysize_t start = 0;
        for(polysize_t i = 0; i < sz; i++) {
            memcpy(&(subres_array[start]), (S->polys[i]->elems), chain_alloc_size_array[idx][i] * sizeof(usfixn64));
            start += chain_alloc_size_array[idx][i];
        }
        // *chain_alloc_sizes = lcl_chain_alloc_sizes;
        // *subres_size = sz;
		subres_size_array[idx] = sz;
		freePolysA_spX(S);
		
		subres_array_out[k] = subres_array;
        
	}
	
	freePolynomial_spX(&a_dusp); // TODO: we can do better
	freePolynomial_spX(&b_dusp);
	
	if ((k2-k1)) {
		
		polysize_t test_subres_size_array = subres_size_array[0]; 
		for(int j = 1; j < k2-k1 && (test_subres_size_array != subres_size_array[j]); j++) {
			test_subres_size_array = -1;
			break;
		}

		subres_size[thr_k] = test_subres_size_array;
		chain_alloc_sizes[thr_k] = chain_alloc_size_array[0];


		free(subres_size_array);
		for(int i = 1; i < k2-k1; i++) {
			free(chain_alloc_size_array[i]);
		}
		free(chain_alloc_size_array);
	} else {
		subres_size[thr_k] = -1;
		chain_alloc_sizes[thr_k] = NULL;
	}
}

int DUSP::bivarModularHGCDSubResInForm_withFFT_parallel_spX (elem_t* g, polysize_t* gpd, elem_t* f, polysize_t* fpd, polysize_t k_sub,
                                           polysize_t n, int K, int e, usfixn64 w, usfixn64 w_inv, usfixn64 n_inv,
                                           biSubresPr_t** subres, const Prime_ptr* Pptr, 
										   const std::vector<ExecutorThreadPool::threadID>& workers)
{
    montgomery_triple P;
    init_montgomery_triple (&P, (unsigned long long) Pptr->prime);
    usfixn64* g_eval = convert_bivarElemPoly_to_usfixn64Poly(g, gpd, n); // (gpd[1]+1) * n
    usfixn64* f_eval = convert_bivarElemPoly_to_usfixn64Poly(f, fpd, n); // (fpd[1]+1) * n
    g_eval = bivarRowWiseDFT_spX_parallel(g_eval, gpd[1], n, K, e, w, P, workers); // row-wise FFT(g_eval)
    f_eval = bivarRowWiseDFT_spX_parallel(f_eval, fpd[1], n, K, e, w, P, workers); // row-wise FFT(f_eval)
    g_eval = transpose2D_spX(g_eval, gpd[1]+1, n); // n * (gpd[1]+1)
    f_eval = transpose2D_spX(f_eval, fpd[1]+1, n); // n * (fpd[1]+1)
	usfixn64** subres_arrays = (usfixn64**) malloc(n * sizeof(usfixn64*));

    // cilk-for
	int nworkers = workers.size();
	// fprintf(stderr, "TEST: n = %ld, nworkers = %d\n", n, nworkers);
    polysize_t** chain_alloc_sizes = (polysize_t**) malloc((nworkers+1) * sizeof(polysize_t*)); 
	polysize_t* max_chain_alloc_sizes = NULL;
    polysize_t* subres_size = (polysize_t*) calloc((nworkers+1), sizeof(polysize_t)); 
	polysize_t test_subres_size = 0;

	int k1 = 0, k2 = 0;
	for (int k = 0; k < nworkers; k++) {
		//get number of images to compute, rounding k2 down
    	// to nearest multiple of MIN_IMAGES_PER_THREAD
    	k2 = (n / (nworkers+1)); //+1 becuase current thread does work to
    	k2 = (k2 / THREAD_ELEMENTS_BOUNDARY) * THREAD_ELEMENTS_BOUNDARY;

		// add offset
		k2 += k1;

    	std::function<void()> task = std::bind(modularHGCDSubres_from_usfixn64_parallel,
			k, k1, k2, g_eval, gpd[1], f_eval, fpd[1], k_sub, subres_size, chain_alloc_sizes, Pptr, subres_arrays);
    	ExecutorThreadPool::getThreadPool().executeTask(workers[k], task);

		k1 = k2;
	}
    //do the last one on the current thread
    //if nthreads does not evenly divide n, do the extra images on the current thread
	k2 = n;

	// fprintf(stderr, "TEST: launch local worker: %d\n", workerIdx);
	// fprintf(stderr, "TEST: worker %d evaling [%d, %d)\n", nworkers, k1, k2);

	modularHGCDSubres_from_usfixn64_parallel (nworkers, k1, k2, g_eval, gpd[1], f_eval, fpd[1], k_sub,
										subres_size, chain_alloc_sizes, Pptr, subres_arrays);

	//sync threads
	ExecutorThreadPool::getThreadPool().waitForThreads(workers);

	for (int i = 0; i <= nworkers && (chain_alloc_sizes[i] == NULL); i++) {
		fprintf(stderr, "DUSP_Parallel Error: works are not devided equally\n");
		exit(1);
	}

	test_subres_size = subres_size[0];
	int isGoodOmega = 1;
	if(test_subres_size == -1) {
		isGoodOmega = 0;
	} else {
		for (int i = 1; i <= nworkers && (test_subres_size != subres_size[i]); i++) {
			isGoodOmega = 0;
			break;
		}
	}
	
	if(!isGoodOmega) {
		// fprintf(stderr, "TEST: (g, f) test_subres_size (%ld) != subres_size (%ld)\n", test_subres_size, subres_size);
		free(g_eval);
		free(f_eval);
		for (polysize_t j = 0; j < n; j++) {
			free(subres_arrays[j]);
		}
		free(subres_arrays);
		if (chain_alloc_sizes != NULL) {
			for(int j = 0; j <= nworkers && (chain_alloc_sizes[j] != NULL); j++) {
				free(chain_alloc_sizes[j]);
			}
			free(chain_alloc_sizes);
		}
		free(subres_size);

		// pool.returnThreads(workers);

		return 0;		
	}

	// TODO: 
	max_chain_alloc_sizes = (polysize_t*) malloc(test_subres_size * sizeof(polysize_t));
	for(polysize_t j = 0; j < test_subres_size; j++) {
		max_chain_alloc_sizes[j] = chain_alloc_sizes[0][j];
	}
	for(int j = 0; j <= nworkers; j++) {
		free(chain_alloc_sizes[j]);
	}
	free(chain_alloc_sizes);

    polysize_t alpha = gpd[0] + fpd[0];
    // bivarSubresFFT_t* subres_out = (bivarSubresFFT_t*) malloc(sizeof(bivarSubresFFT_t));
    // subres_out->polys = _compute_DFTInv_modularSubres(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr);
    // subres_out->chain_size = test_subres_size;
    // subres_out->alloc_sizes = max_chain_alloc_sizes;
    // *subres = subres_out;

    usfixn64** polys = compute_DFTInv_modularSubres_parallel(subres_arrays, n, K, e, w_inv, n_inv, test_subres_size, alpha, max_chain_alloc_sizes, Pptr, workers);
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

int DUSP::hgcdBiModularFFTSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, int kk, AltArrsZ_t** ksubres, int* chain_size, int maxThreads)
{
	if (a->nvar != 2 || b->nvar != 2) {
		fprintf (stderr, "DBZP Error, In hgcdBiModularFFTSubresultantChainZ, a->nvar, b->nvar must be 2.\n");
		exit (1);
	}
	AltArrsZ_t *node, *head=NULL, *tail=NULL;
	if (isZero_AAZ (a)) {
		if (isZero_AAZ (b)) {
			*ksubres = NULL;
			*chain_size = 0;
		} else {
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = deepCopyPolynomial_AAZ (b);
			node->next = NULL;
			head = node;
			tail = node;
			node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
			node->poly = NULL;
			node->next = NULL;
			tail->next = node;
			tail = node; // todo
			*chain_size = 2;
			*ksubres = head;
		}
		return 1;
	} else if (isZero_AAZ (b)) {
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = deepCopyPolynomial_AAZ (a);
		node->next = NULL;
		head = node;
		tail = node;
		node = (AltArrsZ_t*) malloc (sizeof (AltArrsZ_t));
		node->poly = NULL;
		node->next = NULL;
		tail->next = node;
		tail = node; // todo
		*chain_size = 2;
		*ksubres = head;
		return 1;
	} else if (mainLeadingDegree_AAZ (a) < mainLeadingDegree_AAZ (b)) {
		return DUSP::hgcdBiModularFFTSubresultantChainZ (b, a, kk, ksubres, chain_size, maxThreads);
	}
	int sz = 1;
	// convert a, b
	mpz_t *aZ, *bZ;
	polysize_t a_degs[2], b_degs[2];
	aZ = convertFromAltArrZ_DBZP (a, a_degs);
	bZ = convertFromAltArrZ_DBZP (b, b_degs);
	polysize_t min_Xdeg = MIN (a_degs[1], b_degs[1]);

	biSubresZ_t* subres_work = NULL;
	biSubresZ_t* subres_work_prev = NULL;

	mpz_t m; // product of primes
	mpz_init_set_ui (m, 1l);
	mpz_t halfm, newm;
	mpz_inits (halfm, newm, NULL);
	mpz_t g, s, t, mpz_pr; // bezout coefs and (cast) prime
	mpz_inits (g, s, t, mpz_pr, NULL);

	polysize_t a_tdeg = (a_degs[0]+1)*(a_degs[1]+1);
	polysize_t b_tdeg = (b_degs[0]+1)*(b_degs[1]+1);

	AltArrZ_t* tmp;
	elem_t *modA = (elem_t*) malloc (a_tdeg*sizeof(elem_t));
	elem_t *modB = (elem_t*) malloc (b_tdeg*sizeof(elem_t));
	unsigned long uprime;
	biSubresPr_t* modSubres = NULL;
	unsigned long coef_out = 0;

	int primeIdx = 0;
	const int prime_table_size = 389;
	usfixn64 p_u64, w, w_inv, n_inv;
    montgomery_triple P;
	int isGoodOmega, isSubresEq, no_failed_omega;
	const int max_no_failed_omega = 5;

	int K, e;
	polysize_t n = _compute_N_gt_By(a_degs, b_degs, &K, &e);
	polysize_t chainSize = (b_degs[1] > 2) ? b_degs[1] : 2;

	ExecutorThreadPool& pool = ExecutorThreadPool::getThreadPool(); // start the ThreadPool

	int nthreads; // No. Threads
	if (maxThreads < 0) {
		int nUnivarImages = b_degs[1]*a_degs[0] + b_degs[0]*a_degs[1] + 1;
		// fprintf(stderr, "TEST: Univar Images Needed: %d\n", nUnivarImages );

		nthreads = ((nUnivarImages) / THREAD_UNIVAR_MIN_SIZE) - 1; //cur thread counts as one
		nthreads = (chainSize-1) < nthreads ? (chainSize-1) : nthreads;
		// fprintf(stderr, "TEST: nthreads: %d\n", nthreads );
		// fprintf(stderr, "TEST: chain size: %d\n", chainSize );

		infinityNorm_AAZ(a, s);
		infinityNorm_AAZ(b, t);
		if (mpz_cmp(s, t) < 0) {
			mpz_swap(s,t); //make s max
		}
		int nprimes = mpz_size(s) * a_degs[1] * 3/2; //an estimate for nthreads
		if (nprimes < 125) {
			//do nothing
		} else if (nprimes < 150) {
			if (nthreads < 1) {
				nthreads = 1;
			}
		} else if (nprimes < 175) {
			if (nthreads < 2) {
				nthreads = 2;
			}
		} else if (nprimes < 250) {
			if (nthreads < 3) {
				nthreads = 3;
			}
		} else if (nprimes < 300) {
			if (nthreads < 4) {
				nthreads = 4;
			}
		} else {
			nthreads = ExecutorThreadPool::maxThreads;
		}
		nthreads = nthreads < 0 ? 0 : nthreads;
	} else {
		nthreads = maxThreads;
	}

	std::vector<ExecutorThreadPool::threadID> workers;
	nthreads = pool.obtainThreads(nthreads, workers);
	// fprintf(stderr, "TEST: final nthreads: %d\n", nthreads);
	// float elapsed;
	// unsigned long long start;
	// float CRTTime = 0.0f;

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Parallel- bivarPS-FFT nworkers=%d\n", (int) workers.size());
#endif


	int subResEq[nthreads+1];
	mpz_t tmpZ[nthreads+1];
	for (int i = 0; i < nthreads+1; ++i) {
		mpz_init(tmpZ[i]);
		subResEq[i] = 0;
	}

	// fprintf (stderr, "start testing primes...\n"); //TEST
	for (; primeIdx < prime_table_size; ++primeIdx) {

		// @note see the biModularFFTSubresultantChainZ for more details...
		Prime_ptr *Pptr = smallprimefield_get_prime_constants(subset_fourier_primes_u64_table[primeIdx]);
		uprime = (unsigned long) Pptr->prime;

		if(!smallprimefield_convert_in (mpz_fdiv_ui(aZ[0], uprime), Pptr) ||
			!smallprimefield_convert_in (mpz_fdiv_ui(bZ[0], uprime), Pptr)) {
			continue;
		}

		p_u64 = (unsigned long long) Pptr->prime;
		init_montgomery_triple (&P, p_u64);
		compute_nth_root_of_unity_for_small_prime (p_u64, &w, n);

		// @note we check it now in bivarModularSubResInForm_withFFT_spX
		// @note see the biModularFFTSubresultantChainZ for more details...

		gmp_inv_mod_p_u64 (&w_inv, w, p_u64);
	    gmp_inv_mod_p_u64 (&n_inv, n, p_u64);

		convertIn_GLOBAL_ptr(&w, &P);
	    convertIn_GLOBAL_ptr(&w_inv, &P);
    	convertIn_GLOBAL_ptr(&n_inv, &P);

		// if good prime:
		// convert modular images
		// Convert a to Z_p[0]
		for (polysize_t i = 0; i < a_tdeg; i++) {
			modA[i] = smallprimefield_convert_in (mpz_fdiv_ui(aZ[i], uprime), Pptr);
		}
		// Convert f to Z_p[0]
		for (polysize_t i = 0; i < b_tdeg; i++) {
			modB[i] = smallprimefield_convert_in (mpz_fdiv_ui(bZ[i], uprime), Pptr);
		}

		// call subresultant mod this prime
		isGoodOmega = bivarModularHGCDSubResInForm_withFFT_parallel_spX(modA, a_degs, modB, b_degs, kk, n, K, e, w, w_inv, n_inv, &modSubres, Pptr, workers);

		if(!isGoodOmega) {
			if (no_failed_omega < max_no_failed_omega) {
				no_failed_omega++;
				continue;
			} else {
				mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
				for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
				free(aZ);
				for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
				free(bZ);
				free(modA); free(modB);
				free(Pptr);
				for (int i = 0; i < nthreads+1; ++i) {
					mpz_clear(tmpZ[i]);
				}
				return 0;
			}
	   }

		mpz_set_si(mpz_pr, Pptr->prime);
		mpz_gcdext(g, s, t, mpz_pr, m);
		mpz_mul(newm, m, mpz_pr);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		// construct subres_work at the first round:
		// fprintf(stderr, "modSubres->n = %ld\n", modSubres->n);
		if (subres_work == NULL) {
			subres_work = makeBiSubresultant_DBZP (modSubres->n);
			for (polysize_t j = 0; j < modSubres->n; j++) {
				subres_work->deg[j][0] = modSubres->deg[j][0];
				subres_work->deg[j][1] = modSubres->deg[j][1];
				subres_work->size[j]   = modSubres->size[j];
				subres_work->coefs[j] = (mpz_t*) malloc (subres_work->size[j]*sizeof(mpz_t));
				for (polysize_t l = 0; l < subres_work->size[j]; l++) {
					mpz_init (subres_work->coefs[j][l]);
				}
			}
		}

		// isSubresEq = 1;
		// elem_t coef_out;
		// for (polysize_t k = 0; k < modSubres->n; ++k) {
		// 	if (modSubres->size[k] == 0 || modSubres->coefs == NULL) {
		// 		continue;
		// 	}
		// 	polysize_t lt = modSubres->size[k];
		// 	subres_work->deg[k][0] = modSubres->deg[k][0];
		// 	subres_work->deg[k][1] = modSubres->deg[k][1];
		// 	subres_work->size[k] = modSubres->size[k]; // update lt
		// 	for (polysize_t i = 0; i < lt; ++i) {
		// 		mpz_set(tmpZ, subres_work->coefs[k][i]);
		// 		coef_out = (unsigned long) smallprimefield_convert_out (modSubres->coefs[k][i], Pptr);
		// 		// optimized-CRT (Garner's algorithm)
		// 		mpz_sub_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], s);
		// 		mpz_mul (subres_work->coefs[k][i], subres_work->coefs[k][i], mpz_pr);
		// 		mpz_add_ui (subres_work->coefs[k][i], subres_work->coefs[k][i], coef_out);
		// 		mpz_mod (subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		if (mpz_cmp (subres_work->coefs[k][i], halfm) > 0) {
		// 			mpz_sub(subres_work->coefs[k][i], subres_work->coefs[k][i], newm);
		// 		}

		// 		isSubresEq = isSubresEq && (mpz_cmp(tmpZ, subres_work->coefs[k][i]) == 0);
		// 	}
		// }

		int isSubresEq = 1;
		// _startTimer(&start);
		if (primeIdx < 50) {
			biModularSubresultantCRT_parallel(0, modSubres->n, subres_work, modSubres, &subResEq[0], 
											tmpZ[0], mpz_pr, s, halfm, newm, Pptr);
			isSubresEq = subResEq[0];
		} else {
			int amtWork = 0;
			for (polysize_t k = 0; k < modSubres->n; ++k) {
				amtWork += modSubres->size[k];
			}
			amtWork /= (nthreads+1);
			int k1 = 0, k2 = 0;
			int workerIdx = 0;
			int curWork = 0;
			for (int k = 0; k < modSubres->n; ++k) {
				curWork += modSubres->size[k];
				if (curWork > amtWork) {
					k2 = k;
					if (!subResEq[workerIdx]) {
						std::function<void()> task = std::bind(biModularSubresultantCRT_parallel, k1, k2,
							subres_work, modSubres, &subResEq[workerIdx], tmpZ[workerIdx], mpz_pr, s, halfm, newm, Pptr);
						ExecutorThreadPool::getThreadPool().executeTask(workers[workerIdx], task);
						// fprintf(stderr, "TEST: worker %d CRT [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
					} else {
						// fprintf(stderr, "TEST: worker %d CRT [%d,%d), skipping work\n", workerIdx, k1, k2, curWork);
					}
					k1 = k;
					++workerIdx;
					curWork = 0;
				}
			}
			k2 = modSubres->n;
			if (!subResEq[workerIdx]) {
				// fprintf(stderr, "TEST: main t %d CRT [%d,%d), curWork: %d\n", workerIdx, k1, k2, curWork);
				biModularSubresultantCRT_parallel(k1, k2,
					subres_work, modSubres, &subResEq[workerIdx], tmpZ[workerIdx], mpz_pr, s, halfm, newm, Pptr);
				++workerIdx;
			} else {
				// fprintf(stderr, "TEST: worker %d CRT [%d,%d), skipping work\n", workerIdx, k1, k2, curWork);
			}
			//sync threads
			ExecutorThreadPool::getThreadPool().waitForThreads(workers);
			// _stopTimer(&start, &CRTTime);
			// totalCRTTime += CRTTime;

			//Notice here we use workerIdx and not nthreads+1
			//Not all threads may have done work in the CRT loop
			for (int i = 0; i < workerIdx; ++i) {
				isSubresEq = isSubresEq && subResEq[i];
			}
		}

		free(Pptr);
		if (modSubres != NULL) {
			freeBiSubresultantInForm_spX (modSubres);
		}
		modSubres = NULL;

		if (isSubresEq) {
			AltArrsZ_t *node, *head=NULL, *tail=NULL;
			int sz = 1;
			if (!subres_work->coefs[0] || !subres_work->size[0]) {
				tmp = NULL;
			} else {
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[0], subres_work->deg[0]);
			}
			node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
			node->poly = tmp; tmp = NULL;
			node->next = NULL;
			head = node;
			tail = node;
			for (polysize_t i = 1; i < subres_work->n; i++) {
				if (!subres_work->size[i]) {
					continue;
				}
				tmp = convertToAltArrZ_DBZP (subres_work->coefs[i], subres_work->deg[i]);
				if (tmp == NULL || tmp->size == 0) {
					freePolynomial_AAZ (tmp);
					tmp = NULL;
					continue;
				}
				node = (AltArrsZ_t*) malloc (sizeof(AltArrsZ_t));
				node->poly = tmp; tmp = NULL;
				node->next = NULL;
				tail->next = node;
				tail = node;
				sz++;
			}
			freeBiSubresultant_DBZP (subres_work);
			freeBiSubresultant_DBZP (subres_work_prev);
			mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
			for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
			free(aZ);
			for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
			free(bZ);
			free(modA); free(modB);
			for (int i = 0; i < nthreads+1; ++i) {
				mpz_clear(tmpZ[i]);
			}
			*ksubres = head;
			*chain_size = sz;
			return 1;
		}
		mpz_set (m, newm);
	}
	mpz_clears (m, halfm, newm, g, s, t, mpz_pr, NULL);
	freeBiSubresultant_DBZP (subres_work);
	freeBiSubresultant_DBZP (subres_work_prev);
	for(polysize_t i = 0; i < (a_degs[1]+1)*(a_degs[0]+1); i++) { mpz_clear(aZ[i]);}
	free(aZ);
	for(polysize_t i = 0; i < (b_degs[1]+1)*(b_degs[0]+1); i++) { mpz_clear(bZ[i]);}
	free(bZ);
	free(modA); free(modB);
	for (int i = 0; i < nthreads+1; ++i) {
		mpz_clear(tmpZ[i]);
	}
	// fprintf (stderr, "In hgcdBiModularFFTSubresultantChain_DBZP, all primes failed\n"); // TEST
	*ksubres = NULL;
	*chain_size = 0;
	return 0;
}

int DUSP::ParallelSubresultantChainZ (AltArrZ_t* a, AltArrZ_t* b, long min_mdeg, long max_pdeg, AltArrsZ_t** subres, int* chain_size, int maxthread)
{
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Parallel SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", a->nvar, mainLeadingDegree_AAZ(a), mainLeadingDegree_AAZ(b));
#endif
	int no_prime = 0;
	AltArrsZ_t* subresRev = NULL;
	int size = 0, flag = 0;
	if (min_mdeg < 12 && max_pdeg < 300) {
		flag = DUSP::biModularSubresultantChainZ(a, b, &subresRev, &size, &no_prime, maxthread);
	} else if (max_pdeg < 150) {
		flag = DUSP::biModularSubresultantChainZ(a, b, &subresRev, &size, &no_prime, maxthread);
	} else {
		flag = DUSP::biModularFFTSubresultantChainZ(a, b, &subresRev, &size, &no_prime, maxthread); 
	}

	// flag = DUSP::biModularSubresultantChainZ(a, b, &subresRev, &size, &no_prime, maxthread);

	// const char* sym[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};  
	// fprintf (stderr, "In biModularFFTSubresultantChainZ a := ");
	// printPoly_AAZ (stderr, a, sym, a->nvar);
	// fprintf (stderr, "\nIn biModularFFTSubresultantChainZ b := ");
	// printPoly_AAZ (stderr, b, sym, b->nvar);
	// fprintf(stderr, "\n");

	if (!flag) { 
#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Parallel- Subres. FAILED!\n");
#endif
		SubresultantChainWrapperZ_rev(DUCOS, a, b, &subresRev, &size);
	}
	*chain_size = size;

    if (size > 1){
		AltArrsZ_t* curr = subresRev;
		AltArrsZ_t* prev = NULL;
		AltArrsZ_t* next = NULL;
		while (curr != NULL){
			next = curr->next;
			curr->next = prev;
			prev = curr;
			curr = next;
		}
			*subres = prev;
    } else {
		*subres = subresRev;
	}

#if defined(PRINT_WHICH_SUBRES) && (PRINT_WHICH_SUBRES == 1)
	fprintf(stderr, "*** Done! Parallel SubresultantChain [nvar = %d, mdeg(P) = %d, mdeg(Q) = %d]\n", a->nvar, mainLeadingDegree_AAZ(a), mainLeadingDegree_AAZ(b));
#endif

	return flag; 
}


