#include "ModularPolynomial/DUSP_FFT_Support.h"


void fastMulPolynomialInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** c, Prime_ptr* Pptr)
{  
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "fMul-FFT\n");
#endif
    
    if (degPolynomial_spX (a) < FFTMUL_CROSSOVER || degPolynomial_spX (b) < FFTMUL_CROSSOVER) {
        // plainMulPolynomialsInForm_spX (a, b, c, Pptr);
        mulPolynomialsInForm_spX (a, b, c, Pptr);
        return;
    }

    polysize_t max_deg = MAX_spX (degPolynomial_spX (a), degPolynomial_spX (b));

    if (Pptr->prime < 1337006139375617UL || Pptr->prime > 9214646312576745473UL) {
        // we likely not using fouried_primes_u64, so,
        KaratsubaMulPolynomialsInForm_spX (a, b, c, Pptr);
        return;
    } else {
        // fprintf (stderr, "6th-step FFT...\n"); // TEST
        // TODO: k is not efficient 
        int K=4, e=1;
        polysize_t l = max_deg+1; // technically, it's max alloc size
        while (l > 0) { l = l/K; e++; } 
        // fprintf (stderr, "DUSP FFT Params: max_deg = %ld\n", max_deg); // TEST
        // fprintf (stderr, "DUSP FFT Params: e = %d\n", e); // TEST
        polysize_t N = pow (K, e); // max_deg + 1 <= N
        polysize_t N2 = N>>1; // polynomials have alloc_size at most N2
        // fprintf (stderr, "DUSP FFT Params: K^e = %ld\n", N); // TEST
        // fprintf (stderr, "DUSP FFT Params: K^e/2 = %ld\n", N2); // TEST
        if (max_deg >= N2) {
            fprintf (stderr, "max_deg = %ld > K^e/2 = %d^%d", max_deg, K, e);
            exit(1);
        }
        
        usfixn64 p_u64 = (unsigned long long) Pptr->prime; //  fourier_primes_u64_table[?];
        // init montgomery triple:
        montgomery_triple P;
        init_montgomery_triple (&P, p_u64);

        usfixn64* x_data = (usfixn64*) calloc (N, sizeof(usfixn64));
        usfixn64* y_data = (usfixn64*) calloc (N, sizeof(usfixn64));
        for (long i = 0; i <= degPolynomial_spX (a); ++i) {
            x_data[i] = (usfixn64) a->elems[i]; }
        for (long i = 0; i <= degPolynomial_spX (b); ++i) {
            y_data[i] = (usfixn64) b->elems[i]; }

        // init dft parameters 
        usfixn64 omega, omega_inv, n_inv;
        //computing N'th root of unity in small prime field (omega).
        compute_nth_root_of_unity_for_small_prime (p_u64, &omega, N);
        //computing inverse of N'th root of unity (omega_inv).
        gmp_inv_mod_p_u64 (&omega_inv, omega, p_u64);
        //computing inverse of N (N^{-1})
        gmp_inv_mod_p_u64 (&n_inv, N, p_u64);
        // convert to montgomery representation
        convertIn_GLOBAL_ptr(&omega, &P);
        convertIn_GLOBAL_ptr(&omega_inv, &P);
        //x=DFT_n(x)
        DFT_general_small_elements (x_data, K, e, &P, omega);
        //y=DFT_n(y)
        DFT_general_small_elements (y_data, K, e, &P, omega);
        //x=DFT_n(x)*DFT_n(y)
        convolution_mult_small_elements (x_data, y_data, N, &P);
        //x=DFT^{-1}_n(x)
        DFT_general_small_elements (x_data, K, e, &P, omega_inv);
        convertIn_GLOBAL_ptr(&n_inv, &P);
        for (int i = 0; i < N; i++) {
            x_data[i] = mult_ab_mod_p (x_data[i], n_inv, &P);
        }

        // make output
        setCoefsInForm_spX (c, (elem_t*) x_data, N);
        normalize_spX (c);

        free (x_data);
        free (y_data);
        return;
    }

    polysize_t r = Log2_spX (max_deg)+2;
    polysize_t n = 1<<r;
    polysize_t n2 = n<<1;
    // std::cout << "r = " << r << "\tn = " << n << "\tn2 = " << n2 << std::endl;
    sfixn* f = (sfixn*) calloc (n, sizeof (sfixn));
    sfixn* g = (sfixn*) calloc (n, sizeof (sfixn));
    for (long i = 0; i <= degPolynomial_spX (a); ++i) {
	f[i] = (long) a->elems[i]; }
    for (long i = 0; i <= degPolynomial_spX (b); ++i) {
	g[i] = (long) b->elems[i]; }
    // generate Omega and OmegaInv:
    MONTP_OPT2_AS_GENE* pPtr = (MONTP_OPT2_AS_GENE*) my_malloc(sizeof (MONTP_OPT2_AS_GENE));
    EX_MontP_Init_OPT2_AS_GENE (pPtr, Pptr->prime); // 
    sfixn R=1L<<(pPtr->Rpow); 
    R=R%pPtr->P;
    int H = 1024, j = 0;

    sfixn *SB1, *Omega, *OmegaInv, *pwm_fg;
    Omega = (sfixn*) calloc (n2<<1, sizeof (sfixn));
    OmegaInv = (sfixn*) calloc (n2<<1, sizeof (sfixn));
    sfixn invn = inverseMod (n, pPtr->P);
    invn = MulMod (R, invn, pPtr->P);
    invn <<= pPtr->Base_Rpow;
    PBPAS::RootsTable2 (n2, r+1, Omega, OmegaInv, pPtr);
    Omega += n2;     
    OmegaInv += n2;
    SB1 = (sfixn*) calloc (n, sizeof (sfixn)); 
    PBPAS::DFT_eff(n, r, f, Omega, pPtr, H, NULL, SB1, 1);
    // //FURERPBPAS1::DFT_eff_p1(K,e,Ap,W,SB1);
    PBPAS::DFT_eff(n, r, g, Omega, pPtr, H, NULL, SB1, 1);
    // //FURERPBPAS1::DFT_eff_p1(K,e,Ap,W,SB1);
    pwm_fg = (sfixn*) calloc (n, sizeof (sfixn));
    for (j = 0; j < n; j++) {
      pwm_fg[j] = (long) smallprimefield_mul ((long long) f[j], (long long) g[j], Pptr);}

    free (f);
    free (g);
    PBPAS::InvDFT_eff(n, r, pwm_fg, OmegaInv, pPtr, H, NULL, SB1, invn, 1);
    // //FURERPBPAS1::InvDFT_eff_p1(K,e,Ap,W,SB1);

    setCoefsInForm_spX (c, (long long*) pwm_fg, n);
    normalize_spX (c);
    // Omega and OmegaInv shouldn't be freed!
    free (SB1);
}

void plainExponentiatePolynomialInForm_FFT_spX (duspoly_t* a, polysize_t n, duspoly_t** an, Prime_ptr* Pptr)
{
    if (n == 0) {
        elem_t one = smallprimefield_convert_in (1, Pptr);
        *an = constPolynomialInForm_spX (one, Pptr);
        return;
    } else if (n == 1) {
        *an = deepCopyPolynomial_spX (a);
        return;
    }

    duspoly_t* mul = NULL;
    duspoly_t* tmp = NULL;
    duspoly_t* b = deepCopyPolynomial_spX (a);

    long nn =  (long) n;

    while (nn > 1) {
        if (nn & 1) {            
            if (isZero_spX (mul)) {
                mul = deepCopyPolynomial_spX (b);
            } else {
                // plainMulPolynomialsInForm_spX (mul, b, &tmp, Pptr);
                fastMulPolynomialInForm_FFT_spX (mul, b, &tmp, Pptr);
                freePolynomial_spX (&mul);
                mul = tmp; tmp = NULL;
            }
        }
        
        // plainSqrPolynomialInForm_spX (b, &tmp, Pptr);
        // plainMulPolynomialsInForm_spX (b, b, &tmp, Pptr);
        fastMulPolynomialInForm_FFT_spX (b, b, &tmp, Pptr);
        freePolynomial_spX (&b);
        b = tmp; tmp = NULL;

        nn >>= 1;
    }

    if (isZero_spX (mul)) {
        mul = deepCopyPolynomial_spX (b);
    } else {
        // plainMulPolynomialsInForm_spX (mul, b, &tmp, Pptr);
        fastMulPolynomialInForm_FFT_spX (mul, b, &tmp, Pptr);
        freePolynomial_spX (&mul);
        mul = tmp;
    }


    freePolynomial_spX (&b);
    *an = mul;
}


// l should be power of 2
void fastPSInversionInForm_FFT_spX (duspoly_t* a, duspoly_t** ai, polysize_t lgdg, Prime_ptr* Pptr)
{
    if (a->elems[0] == 0) {
    	fprintf (stderr, "DUSP Error, input polynomial has a wrong format in power series inversion!\n");
    	exit(1);
    }

    if (Pptr->prime != 4179340454199820289) {
    	fprintf (stderr, "DUSP Error,  fastPSInversionInForm_FFT doesn't work with prime = %lld!\n", Pptr->prime);
    	exit(1);
    }
    
    duspoly_t* aa = NULL; // TODO: 
    elem_t e0 = 0;
    int isOne = 0;
    e0 = smallprimefield_inv (a->elems[0], Pptr); // e^{-1}
    
    if (e0 == smallprimefield_convert_in (1, Pptr)) {
    	isOne = 1;
	   aa = deepCopyPolynomial_spX (a);
    } else {
	   scalarMulPolynomialInForm_spX (a, e0, &aa, Pptr);
    }

    // fprintf (stderr, "a = "); // TEST
    // printPolynomial_spX (aa, Pptr);
     
    polysize_t r = lgdg;    // r
    polysize_t n = 1<<lgdg; // n
    polysize_t min = MIN_spX (n, degPolynomial_spX (aa));
    // if (n <= degPolynomial_spX (aa)) {
    //   fprintf (stderr, "DUSP Error, In power series inversion, 2^ldgdg <= deg(a) !\n");
    // 	exit(1);
    // }
    
    sfixn* f = (sfixn*) calloc (n, sizeof (sfixn));
    sfixn* g = (sfixn*) calloc (n, sizeof (sfixn));

    for (polysize_t i = 0; i <= min; ++i) {
	   f[i] = (long) aa->elems[i];
    }

    g[0] = (long) smallprimefield_convert_in (1, Pptr);
    
    // fprintf (stderr, "reallocate(f, %ld) = \n", n); // TEST
    // for (polysize_t i = 0; i <= n; ++i) {
    // 	fprintf (stderr, "f[%ld] = %ld\n", i, f[i]); 
    // }

    // fprintf (stderr, "\nreallocate(g, %ld) = \n", n); // TEST

    // for (polysize_t i = 0; i <= n; ++i) {
    // 	fprintf (stderr, "g[%ld] = %ld\n", i, g[i]); 
    // }

    
    // generate Omega and OmegaInv:
    // TODO: Don't need pPtr 
    MONTP_OPT2_AS_GENE* pPtr = (MONTP_OPT2_AS_GENE*) my_malloc(sizeof (MONTP_OPT2_AS_GENE));
    EX_MontP_Init_OPT2_AS_GENE (pPtr, Pptr->prime); // 

    sfixn R=1L<<(pPtr->Rpow); 
    R=R%pPtr->P;

    int H = 1024; 
    
    sfixn* SB1;
    sfixn* Omega;
    sfixn* OmegaInv;

    sfixn* dft_f;
    sfixn* dft_g;
    sfixn* pwm_fg;

    sfixn* h;
    sfixn* pwm_gh;

    sfixn* t;
    
    sfixn invn;
    polysize_t pow = 1;
    polysize_t pow1;
    polysize_t pow2;
    polysize_t j;

    for (polysize_t i = 1; i <= r; i++) {

    	pow1 = pow;    // pow1 = 2^{i-1}
    	pow <<= 1;     // pow  = 2^i
    	pow2 = pow<<1; // pow2 = 2^{i+1}
    	
    	Omega = (sfixn*) calloc (pow2<<1, sizeof (sfixn));
    	OmegaInv = (sfixn*) calloc (pow2<<1, sizeof (sfixn));

    	invn = inverseMod (pow, pPtr->P);
    	invn = MulMod (R, invn, pPtr->P);
    	invn <<= pPtr->Base_Rpow;
        
    	PBPAS::RootsTable2 (pow2, i+1, Omega, OmegaInv, pPtr);
    	Omega += pow2;     
    	OmegaInv += pow2;
    	
    	// leftSplitPolynomial of pow:
    	dft_f = (sfixn*) malloc (sizeof (sfixn)*pow);
    	for (j = 0; j < pow; ++j) {
    	    dft_f[j] = f[j];
    	}
    	
    	dft_g = (sfixn*) calloc (pow, sizeof (sfixn));
    	for (j = 0; j < pow; ++j) {
    	    dft_g[j] = g[j];
    	}

    	SB1 = (sfixn*) my_calloc (pow, sizeof (sfixn)); 
    	PBPAS::DFT_eff(pow, i, dft_f, Omega, pPtr, H, NULL, SB1, 1);
    	// //FURERPBPAS1::DFT_eff_p1(K,e,Ap,W,SB1);
    	
    	PBPAS::DFT_eff(pow, i, dft_g, Omega, pPtr, H, NULL, SB1, 1);
    	// //FURERPBPAS1::DFT_eff_p1(K,e,Ap,W,SB1);

    	pwm_fg = (sfixn*) calloc (pow, sizeof (sfixn));

    	for (j = 0; j < pow; j++) {
    	    pwm_fg[j] = (long) smallprimefield_mul ((long long) dft_f[j], (long long) dft_g[j], Pptr);
    	}

    	free (dft_f); dft_f = NULL;

    	PBPAS::InvDFT_eff(pow, i, pwm_fg, OmegaInv, pPtr, H, NULL, SB1, invn, 1);
    	// //FURERPBPAS1::InvDFT_eff_p1(K,e,Ap,W,SB1);
    	
    	h = (sfixn*) calloc (pow, sizeof (sfixn)); 
    	
    	
    	for (j = 0; j < pow1; ++j) {
    	    h[j] = pwm_fg[j + pow1];
    	}

    	free (pwm_fg); pwm_fg = NULL;

    	PBPAS::DFT_eff(pow, i, h, Omega, pPtr, H, NULL, SB1, 1);
    	//FURERPBPAS1::DFT_eff_p1(K,e,Ap,W,SB1);
    	
    	pwm_gh = (sfixn*) calloc (pow, sizeof (sfixn));

    	for (j = 0; j < pow; ++j) { // TODO: ...
    	    pwm_gh[j] = (long) smallprimefield_mul ((long long) dft_g[j], (long long) h[j], Pptr);
    	}

    	free (h); h = NULL;
    	free (dft_g); dft_g = NULL;

    	PBPAS::InvDFT_eff(pow, i, pwm_gh, OmegaInv, pPtr, H, NULL, SB1, invn, 1);

    	free (SB1); SB1 = NULL;
    	//FURERPBPAS1::InvDFT_eff_p1(K,e,Ap,W,SB1);
    	
    	t = (sfixn*) calloc (pow, sizeof (sfixn));
    	
    	for (j = 0; j < pow1; ++j) {
    	    t[j+pow1] = pwm_gh[j];
    	}

    	free (pwm_gh); pwm_gh = NULL;
    	
    	for (j = pow1; j < pow; ++j) {
    	    g[j] = (long) smallprimefield_sub ((long long) g[j], (long long) t[j], Pptr);
    	}

    	free (t); t = NULL;
    }

    // free (f);
    freePolynomial_spX (&aa); // TODO: 
    
    setCoefsInForm_spX (ai, (long long*) g, n);
    normalize_spX (ai);

    if (isZero_spX (*ai) || isOne) {
	   return;
    }

    duspoly_t* aii = NULL;
    scalarMulPolynomialInForm_spX (*ai, e0, &aii, Pptr);

    freePolynomial_spX (&*ai);
    *ai = aii;

}

void fastDivPolynomialInForm_wPSInv_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** r, duspoly_t** q, Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "fDiv-FFT\n");
#endif

  // a / 0 = Error
    if (isZero_spX (b)) {
    	fprintf (stderr, "DUSP Error: Division by zero!\n");
    	exit (1);
    }
    
    // 0 / b = 0
    if (isZero_spX (a)) {
    	*r = NULL;
    	*q = NULL;
    	return;
    }

    // a/a = 1
    if (isEqual_spX (a, b)) {
    	*r = NULL;
    	*q = constPolynomial_spX (1, Pptr);
    	return;
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
    	*r = deepCopyPolynomial_spX (a);
    	*q = NULL;
    	return;
    }

    if (deg_a < PLAINDIV_CROSSOVER || deg_b < PLAINDIV_CROSSOVER) {
        plainDivPolynomialsInForm_spX (a, b, r, q, Pptr);
        return;
    }

    if (deg_a < 1024  && deg_b < 1024 && deg_a - deg_b < 11 ) {
        plainDivPolynomialsInForm_spX (a, b, r, q, Pptr);
        return;
    }


    duspoly_t* ar = NULL;
    duspoly_t* br = NULL;
    duspoly_t* br_inv = NULL;
    duspoly_t* qr = NULL;
    duspoly_t* rr = NULL;

    polysize_t ldeg = logceiling (deg_a - deg_b + 1); // ceil (Log2_spX (deg_a - deg_b + 1));
	
    reversePolynomial_spX (a, &ar, deg_a, Pptr);
    reversePolynomial_spX (b, &br, deg_b, Pptr);
    
    fastPSInversionInForm_FFT_spX (br, &br_inv, ldeg, Pptr);

    freePolynomial_spX (&br);
    reallocatePolynomial_spX (&ar, deg_a - deg_b + 1);
    reallocatePolynomial_spX (&br_inv, deg_a - deg_b + 1);
   
    fastMulPolynomialInForm_FFT_spX (ar, br_inv, &qr, Pptr);
    reallocatePolynomial_spX (&qr, deg_a - deg_b + 1);
    
    freePolynomial_spX (&ar);
    freePolynomial_spX (&br_inv);

    reversePolynomial_spX (qr, q, deg_a - deg_b, Pptr);

    freePolynomial_spX (&qr);
    
    fastMulPolynomialInForm_FFT_spX (b,*q, &rr, Pptr);
    subPolynomialsInForm_spX (a, rr, r, Pptr);
    
    freePolynomial_spX (&rr);
}

void mulMat4ToVecInForm_FFT_spX_inp (dusmat4_t* M, duspoly_t** u, duspoly_t** v, Prime_ptr* Pptr)
{
    if (M == NULL) {
        freePolynomial_spX (&*u);
        freePolynomial_spX (&*v);
        return;
    }

    if (isIdentityMat4InForm_spX (M, Pptr)) {
        return;
    }

    if (isZero_spX (*u) && isZero_spX (*v)) {
        *u = NULL;
        *v = NULL;
        return;
    }

    duspoly_t* M0 = NULL;
    duspoly_t* M1 = NULL;
    duspoly_t* t0 = NULL;
    duspoly_t* t1 = NULL;
    
    fastMulPolynomialInForm_FFT_spX (M->polys[0], *u, &M0, Pptr);
    fastMulPolynomialInForm_FFT_spX (M->polys[1], *v, &M1, Pptr);
    
    addPolynomialsInForm_spX (M0, M1, &t0, Pptr);    
    
    freePolynomial_spX (&M0);
    freePolynomial_spX (&M1);
    
    fastMulPolynomialInForm_FFT_spX (M->polys[2], *u, &M0, Pptr);
    fastMulPolynomialInForm_FFT_spX (M->polys[3], *v, &M1, Pptr);
    
    addPolynomialsInForm_spX (M0, M1, &t1, Pptr);

    freePolynomial_spX (&M0);
    freePolynomial_spX (&M1);

    freePolynomial_spX (&*u); 
    freePolynomial_spX (&*v); 

    *u = t0;
    *v = t1;
}

// A*B = C
void mulMat4ToMat4InForm_FFT_spX_inp (dusmat4_t* A, dusmat4_t* B, dusmat4_t** C, Prime_ptr* Pptr)
{

    if (A == NULL || B == NULL) {
	return;
    }

    duspoly_t* C_0= B->polys[0];
    duspoly_t* C_1= B->polys[1];
    duspoly_t* C_2= B->polys[2];
    duspoly_t* C_3= B->polys[3];
   
    mulMat4ToVecInForm_FFT_spX_inp (A, &C_0, &C_2, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (A, &C_1, &C_3, Pptr);

    // freeMat4_spX (&B);
  
    dusmat4_t* CC = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    CC->polys[0] = C_0;
    CC->polys[1] = C_1;
    CC->polys[2] = C_2;
    CC->polys[3] = C_3;
    
    *C = CC;
}

void YapHalfGCDMatrixInForm_wPSInv_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, Prime_ptr* Pptr)
{
    
    dusmat4_t* E = NULL;
    identityMat4InForm_spX (&E, Pptr);
    
    if (isZero_spX (b) || isZero_spX (a)) {
    	*M = E;
    	return;
    }

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    duspoly_t* A = NULL;
    duspoly_t* B = NULL;
     ;
    if (deg_a > deg_b) {
    	A = deepCopyPolynomial_spX (a);
    	B = deepCopyPolynomial_spX (b);
    } else if (deg_a <= deg_b) {
	   // TODO : Using Iterative HGCD
	   *M = E;
	   return;
    }
    
    /* 
       else {
       fprintf (stderr, "DUSP Error: In HalfGCD, degrees of polynomials are not appropriate!\n");
       fprintf (stderr, "deg(a) = deg(b) = %lu\n", deg_a);
       exit(1);
       } 
    */
    
    polysize_t m = (deg_a + 1)>>1;

    /* fprintf (stderr, "\n**** YapHalfGCDMatrixInForm_wPSInv_FFT_spX (m = %lu) **** \n", m);  // TEST */
    /* fprintf (stderr, "a = ");  // TEST */
    /* printPolynomialOutForm_spX (a, Pptr); */

    /* fprintf (stderr, "b = ");  // TEST */
    /* printPolynomialOutForm_spX (b, Pptr); */
    
    if (deg_b < m) {
	*M = E;
	return;
    }

    duspoly_t* A0 = NULL;
    duspoly_t* B0 = NULL;
    dusmat4_t* RR = NULL;
    
    rightShiftPolynomial_spX (A, &A0, m);
    rightShiftPolynomial_spX (B, &B0, m);

    /* fprintf (stderr, "    rightShiftPolynomial_spX (A, &A0, m); = ");  // TEST */
    /* printPolynomialOutForm_spX (A0, Pptr); */

    /* fprintf (stderr, "    rightShiftPolynomial_spX (B, &B0, m); = ");  // TEST */
    /* printPolynomialOutForm_spX (B0, Pptr); */
    
    YapHalfGCDMatrixInForm_wPSInv_FFT_spX (A0, B0, &RR, Pptr);  // first recursive call

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);
    
    /* fprintf (stderr, "    YapHalfGCDMatrixInForm_wPSInv_FFT_spX (A0, B0, &RR, Pptr);  // first recursive call\n");  // TEST */
    /* printMat4OutForm_spX (RR, Pptr); */


    mulMat4ToVecInForm_FFT_spX_inp (RR, &A, &B, Pptr); // A = Ap = r_{j-1} , B = Bp = r_{j}

    /* fprintf (stderr, "    RR*A0 = ");  // TEST */
    /* printPolynomialOutForm_spX (A, Pptr); */

    /* fprintf (stderr, "    RR*B0 = ");  // TEST */
    /* printPolynomialOutForm_spX (B, Pptr); */
    /* fprintf (stderr, "\n"); */
    
    /* if (degPolynomial_spX (A) < degPolynomial_spX (B)) {  // TEST */
    /* 	fprintf (stderr, "deg(r_{j-1}) > deg(r_j)\n"); */
    /* 	fprintf (stderr, "r_{j-1} = ");  // TEST */
    /* 	printPolynomialOutForm_spX (A, Pptr); */
	
    /* 	fprintf (stderr, "r_j = ");  // TEST */
    /* 	printPolynomialOutForm_spX (B, Pptr); */

    /* 	exit(1); */
    /* } */
    
    if (degPolynomial_spX (B) < m) { // isZero_spX (B)

	freePolynomial_spX (&A); 
	freeMat4_spX (&E);
	
	*M = RR;
	return;
    }

    duspoly_t* q = NULL;
    duspoly_t* r = NULL;
    duspoly_t* D = NULL;
    
    /* plainDivPolynomialsInForm_spX (A, B, &r, &q, Pptr); // TODO : Fast Monic Div */
    // fastDivPolynomialInForm_wPSInv_spX (A, B, &r, &q, Pptr);
    fastDivPolynomialInForm_wPSInv_FFT_spX (A, B, &r, &q, Pptr);
    
    
    /* fprintf (stderr, "    plainDiv (r) = ");  // TEST */
    /* printPolynomialOutForm_spX (r,  Pptr); */

    /* fprintf (stderr, "    plainDiv (q) = ");  // TEST */
    /* printPolynomialOutForm_spX (q, Pptr); */
    
    freePolynomial_spX (&A);
    
    if (isZero_spX (r)) { 
	
    	negPolynomialInForm_spX (q, Pptr);
    	
    	dusmat4_t* TM = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    	dusmat4_t* MM = NULL;
    	TM->polys[0] = NULL;
    	TM->polys[1] = constPolynomial_spX (1, Pptr);
    	TM->polys[2] = constPolynomial_spX (1, Pptr);
    	TM->polys[3] = q;

    	/* fprintf (stderr, "     TM = \n");  // TEST */
    	/* printMat4OutForm_spX (TM, Pptr); */
    	
    	
    	mulMat4ToMat4InForm_FFT_spX_inp (TM, RR, &MM, Pptr);

    	freeMat4_spX (&TM);
    	freeMat4_spX (&E);
    	freePolynomial_spX (&B);
    	
    	/* fprintf (stderr, "     RR*TM = MM = \n");  // TEST */
    	/* printMat4OutForm_spX (MM, Pptr); */
    	
    	*M = MM;
    	return;
    }

    elem_t lc = leadingCoeffInForm_spX (r); // IT SHOULD BE UNCOMMENT

    /* fprintf (stderr, "lc(r) = %lld\n", lc); // TEST */
    
    lc = smallprimefield_inv (lc, Pptr); // IT SHOULD BE UNCOMMENT

    /* fprintf (stderr, "inv (lc(r)) = %lld\n", lc); // TEST */

    scalarMulPolynomialInForm_spX (r, lc, &D, Pptr);

    /* fprintf (stderr, "    r*inv(lc) = ");  // TEST */
    /* printPolynomialOutForm_spX (D,  Pptr); */
    
    freePolynomial_spX (&r);

    duspoly_t* tq = NULL;
    
    scalarMulPolynomialInForm_spX (q, smallprimefield_sub (0, lc, Pptr), &tq, Pptr);

    /* fprintf (stderr, "    -q*inv(lc) = ");  // TEST */
    /* printPolynomialOutForm_spX (tq,  Pptr); */
        
    freePolynomial_spX (&q);
    
    
    dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = tq;
    
    /* fprintf (stderr, "     Q = \n");  // TEST */
    /* printMat4OutForm_spX (Q, Pptr); */
    
    polysize_t k = 2*m - degPolynomial_spX (B);
    
    /* fprintf (stderr, "k = %lu\n", k); // TEST */
    
    duspoly_t* C0 = NULL;
    duspoly_t* D0 = NULL;
    dusmat4_t* S  = NULL;
    
    rightShiftPolynomial_spX (B, &C0, k);
    rightShiftPolynomial_spX (D , &D0, k);

    freePolynomial_spX (&B);
    freePolynomial_spX (&D);

    YapHalfGCDMatrixInForm_wPSInv_FFT_spX (C0, D0, &S, Pptr); // second recursive call

    /* fprintf (stderr, "     S = \n");  // TEST */
    /* printMat4OutForm_spX (S, Pptr); */
    
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);
    
    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;
    mulMat4ToMat4InForm_FFT_spX_inp (Q, RR, &QRR, Pptr);

    /* fprintf (stderr, "     Q*RR = \n");  // TEST */
    /* printMat4OutForm_spX (QRR, Pptr); */

    mulMat4ToMat4InForm_FFT_spX_inp (S, QRR, &SQRR, Pptr);
    
    freeMat4_spX (&Q);
    freeMat4_spX (&S);
    
    *M = SQRR;
    
    /* fprintf (stderr, "     S*Q*RR = \n");  // TEST */
    /* printMat4OutForm_spX (SQRR, Pptr); */
    
    return;
}

void GCDMatrixInForm_wHGCD_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, Prime_ptr* Pptr)
{

  if (degPolynomial_spX (a) == degPolynomial_spX (b)) {
    GCDMatrixInForm_wHGCD_spX (a, b, M, Pptr);
    return;
  }
  
    // fprintf (stderr, "\n**** GCDMatrixInForm_wHGCD_FFT_spX **** \n"); // TEST
    // fprintf (stderr, "a = ");  // TEST
    // printPolynomial_spX (a, Pptr);

    // fprintf (stderr, "b = ");  // TEST
    // printPolynomial_spX (b, Pptr);
    
    dusmat4_t* RR = NULL;
    YapHalfGCDMatrixInForm_wPSInv_FFT_spX (a, b, &RR, Pptr);
    
    // TODO : implement non in-place mulMat4ToVecInForm 
    duspoly_t* a0 = deepCopyPolynomial_spX (a);
    duspoly_t* b0 = deepCopyPolynomial_spX (b);
    
    mulMat4ToVecInForm_FFT_spX_inp (RR, &a0, &b0, Pptr);

    // fprintf (stderr, "RR = \n"); // TEST
    // printMat4_spX (RR, Pptr);

    // fprintf (stderr, "a0 = ");  // TEST
    // printPolynomial_spX (a0, Pptr);
    
    // fprintf (stderr, "b0 = ");  // TEST
    // printPolynomial_spX (b0, Pptr);
    
    if (isZero_spX (b0)) {
	*M = RR;
	return;
    }

    duspoly_t* r   = NULL;
    duspoly_t* q   = NULL;
    duspoly_t* lcr = NULL;
    duspoly_t* lcq = NULL;
    dusmat4_t* Q   = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    
    fastDivPolynomialInForm_wPSInv_FFT_spX (a0, b0, &r, &q, Pptr); 
    negPolynomialInForm_spX (q, Pptr); // q = -q;
    
    freePolynomial_spX (&a0);

    /* fprintf (stderr, "r = ");  // TEST */
    /* printPolynomial_spX (r, Pptr); */

    /* fprintf (stderr, "q = ");  // TEST */
    /* printPolynomial_spX (q, Pptr); */
    
    if (isZero_spX (r)) {
	
	freePolynomial_spX (&b0);
	
	Q->polys[0] = NULL;
	Q->polys[1] = constPolynomial_spX (1, Pptr);
	Q->polys[2] = constPolynomial_spX (1, Pptr);
	Q->polys[3] = q; // q = -q

	/* fprintf (stderr, "Q = \n"); // TEST */
	/* printMat4_spX (Q, Pptr); */

	dusmat4_t* QRR = NULL;
	
	mulMat4ToMat4InForm_FFT_spX_inp (Q, RR, &QRR, Pptr); 
	
	freeMat4_spX (&Q);
	
	*M = QRR;
	return;
    }
    
    
    elem_t lc = leadingCoeffInForm_spX (r);
    lc = smallprimefield_inv (lc, Pptr);
    
    scalarMulPolynomialInForm_spX (r, lc, &lcr, Pptr);
    scalarMulPolynomialInForm_spX (q, lc, &lcq, Pptr);
    
    freePolynomial_spX (&r);
    freePolynomial_spX (&q);
    
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = lcq;

    /* fprintf (stderr, "Q with lcq = \n"); // TEST */
    /* printMat4_spX (Q, Pptr); */
    
    /* fprintf (stderr, "b0 = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    /* fprintf (stderr, "lcr = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    
    dusmat4_t* S = NULL;
    GCDMatrixInForm_wHGCD_FFT_spX (b0, lcr, &S, Pptr); 

    /* fprintf (stderr, "S = \n"); // TEST */
    /* printMat4_spX (S, Pptr); */

    freePolynomial_spX (&b0);
    freePolynomial_spX (&lcr);

    dusmat4_t* RRQ = NULL;
    dusmat4_t* SRRQ = NULL;
    
    mulMat4ToMat4InForm_FFT_spX_inp (Q, RR, &RRQ, Pptr);
    mulMat4ToMat4InForm_FFT_spX_inp (S, RRQ, &SRRQ, Pptr);

    // fprintf (stderr, "M = \n"); // TEST
    // printMat4_spX (SQRR, Pptr);

    freeMat4_spX (&Q);
    freeMat4_spX (&S);
    
    *M = SRRQ;     
}
    
void ExtGCDInForm_wHGCD_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eGCD-FFT\n");
#endif

    if (isZero_spX (a) || isZero_spX (b)) {
	plainExtGCDInForm_spX (a, b, u, v, g, Pptr);
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    if (deg_a > deg_b) {
	
	dusmat4_t* M    = NULL;
	duspoly_t* M00a = NULL;
	duspoly_t* M01b = NULL;
	duspoly_t* gcd  = NULL;
	
	GCDMatrixInForm_wHGCD_FFT_spX (a, b, &M, Pptr);
	
	fastMulPolynomialInForm_FFT_spX (M->polys[0], a, &M00a, Pptr);
	fastMulPolynomialInForm_FFT_spX (M->polys[1], b, &M01b, Pptr);

	addPolynomialsInForm_spX (M00a, M01b, &gcd, Pptr);

	freePolynomial_spX (&M00a);
	freePolynomial_spX (&M01b);
		
	elem_t lc_gcd = 0;
	duspoly_t* mgcd = NULL;
	monicPolynomialInForm_spX (gcd, &mgcd, &lc_gcd, Pptr);
	freePolynomial_spX (&gcd);
	
	duspoly_t* tu = NULL;
	duspoly_t* tv = NULL;
	lc_gcd = smallprimefield_inv (lc_gcd, Pptr);

	scalarMulPolynomialInForm_spX (M->polys[0], lc_gcd, &tu, Pptr);
	scalarMulPolynomialInForm_spX (M->polys[1], lc_gcd, &tv, Pptr);
    
	*g = mgcd;
	*u = tu;
	*v = tv;
	
	freeMat4_spX (&M);
    
	return;

    } else if (deg_a == deg_b) {

	duspoly_t* r = NULL;
	duspoly_t* q = NULL;
	  
	fastDivPolynomialInForm_wPSInv_FFT_spX (a, b, &r, &q, Pptr); // r = a mod b

	ExtGCDInForm_wHGCD_FFT_spX (b, r, v, u, g, Pptr);

	freePolynomial_spX (&r);
	freePolynomial_spX (&q);
	   
	return;
	
    } else {
	   ExtGCDInForm_wHGCD_FFT_spX (b, a, v, u, g, Pptr);
    }
}


void iterHalfGCDMatrixInForm_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, Prime_ptr* Pptr) 
{
    
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "iHGCD-FFT\n");
#endif

    dusmat4_t* E = NULL;
    identityMat4InForm_spX (&E, Pptr);

    if (isZero_spX (a) || isZero_spX (b)) {
        *M = E;
        return;
    }

    long m = degPolynomial_spX (a) - d;
    // if (m < 0) { // commented for getting the last two remainders: (lc(.)*gcd(f,g), 0)
    //     m = 0;
    // }

    if (degPolynomial_spX (b) <= m) {
        *M = E;
        return;
    }

    duspoly_t* A = deepCopyPolynomial_spX (a);
    duspoly_t* B = deepCopyPolynomial_spX (b);

    duspoly_t* tmp_q = NULL;
    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL; 

   // (normal) extended division algorithm
    while (degPolynomial_spX (B) > m) { 
    
#if DEBUG_PRINT_LINES
        fprintf(stderr, "deg(B) = %lu, \t m = %l \n", degPolynomial_spX (B), m);
        fprintf(stderr, "before div:\n");
        fprintf(stderr, "A = ");
        printPolynomial_spX (A, Pptr);
        fprintf(stderr, "B = ");
        printPolynomial_spX (B, Pptr);
        fprintf(stderr, "E = \n");
        printMat4OutForm_spX (E, Pptr);
#endif

        plainDivPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        // ////////////////////
        // fastDivPolynomialInForm_wPSInv_FFT_spX (A, B, &tmp1, &tmp_q, Pptr);
        // freePolynomial_spX (&A);
        // A = tmp1; 
        // tmp1 = NULL;
        // ///////////////////
        swap_spX (&A, &B);

        // plainMulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
        fastMulPolynomialInForm_FFT_spX (E->polys[2], tmp_q, &tmp1, Pptr);
        subPolynomialsInForm_spX (E->polys[0], tmp1, &tmp2, Pptr);
        freePolynomial_spX (&tmp1);
        freePolynomial_spX (&E->polys[0]);

        E->polys[0] = E->polys[2];
        E->polys[2] = tmp2;
        tmp2 = NULL;

        // plainMulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
        fastMulPolynomialInForm_FFT_spX (E->polys[3], tmp_q, &tmp1, Pptr);
        subPolynomialsInForm_spX (E->polys[1], tmp1, &tmp2, Pptr);
        freePolynomial_spX (&tmp1);
        freePolynomial_spX (&E->polys[1]);

        E->polys[1] = E->polys[3];
        E->polys[3] = tmp2;
        tmp2 = NULL;

        freePolynomial_spX (&tmp_q);

        if (isZero_spX (B)) {
            break;
        }
    }

    freePolynomial_spX (&A);
    freePolynomial_spX (&B);

    *M = E;
    return;

}


void iterHalfGCDMatrixInForm_FFT_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, Prime_ptr* Pptr) 
{
    
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "iHGCD\n");
#endif

    dusmat4_t* E;
    identityMat4InForm_spX (&E, Pptr);

    if (isZero_spX (*a) || isZero_spX (*b)) {
        *M = E;
        return;
    }

    long m = degPolynomial_spX (*a) - d;
    // if (m < 0) { // commented for getting the last two remainders: (lc(.)*gcd(f,g), 0)
    //     m = 0;
    // }

    if (degPolynomial_spX (*b) <= m) {
        *M = E;
        return;
    }

    duspoly_t* A = *a; // deepCopyPolynomial_spX (a);
    duspoly_t* B = *b; // deepCopyPolynomial_spX (b);

    duspoly_t* tmp_q = NULL;
    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL; 

    // (normal) extended division algorithm
    while (degPolynomial_spX (B) > m) { 
    
#if DEBUG_PRINT_LINES
        fprintf(stderr, "deg(B) = %lu, \t m = %l \n", degPolynomial_spX (B), m);
        fprintf(stderr, "before div:\n");
        fprintf(stderr, "A = ");
        printPolynomial_spX (A, Pptr);
        fprintf(stderr, "B = ");
        printPolynomial_spX (B, Pptr);
        fprintf(stderr, "E = \n");
        printMat4OutForm_spX (E, Pptr);
#endif

        // using fast div doesn't help here!
        plainDivPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        swap_spX (&A, &B);

        if (tmp_q == NULL || tmp_q->alloc == 0) {
            swap_spX (&E->polys[0], &E->polys[2]);
            swap_spX (&E->polys[1], &E->polys[3]);            
        } else {

            //plainMulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
            fastMulPolynomialInForm_FFT_spX (E->polys[2], tmp_q, &tmp1, Pptr);
            if (tmp1 == NULL || tmp1->alloc == 0) {
                swap_spX (&E->polys[0], &E->polys[2]);
            } else {
                subPolynomialsInForm_spX (E->polys[0], tmp1, &tmp2, Pptr);
                freePolynomial_spX (&tmp1);
                freePolynomial_spX (&E->polys[0]);

                E->polys[0] = E->polys[2];
                E->polys[2] = tmp2;
                tmp2 = NULL;
            }

            // plainMulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
            fastMulPolynomialInForm_FFT_spX (E->polys[3], tmp_q, &tmp1, Pptr);
            if (tmp1 == NULL || tmp1->alloc == 0) {
                swap_spX (&E->polys[1], &E->polys[3]);            
            } else {                
                subPolynomialsInForm_spX (E->polys[1], tmp1, &tmp2, Pptr);
                freePolynomial_spX (&tmp1);
                freePolynomial_spX (&E->polys[1]);

                E->polys[1] = E->polys[3];
                E->polys[3] = tmp2;
                tmp2 = NULL;
            }
        }

        freePolynomial_spX (&tmp_q);

        if (isZero_spX (B)) {
            break;
        }
    }

    *a = deepCopyPolynomial_spX (A);
    *b = deepCopyPolynomial_spX (B);

    freePolynomial_spX (&A);
    freePolynomial_spX (&B);

    *M = E;
    return;

}

void halfGCDMatrixInForm_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDM-FFT\n");
#endif

    dusmat4_t* E = NULL;
    identityMat4InForm_spX (&E, Pptr);


    if (isZero_spX (a) || isZero_spX (b)) {
        *M = E;
        return;
    }

    if (degPolynomial_spX (a) - degPolynomial_spX (b) >= d ||
        degPolynomial_spX (a) - degPolynomial_spX (b) <  0 ) {
        *M = E;
        return;
    }

    if (degPolynomial_spX (a) < 600 && degPolynomial_spX (b) < 600) {
        halfGCDMatrixInForm_spX (a, b, M, d, Pptr);
        return;
    }

    long m = degPolynomial_spX (a) - 2*(d-1);
    
    // corner case:
    if (m < 0) {
        m = 0;
    }

    duspoly_t* A0 = NULL;
    duspoly_t* B0 = NULL; 
      
    rightShiftPolynomial_spX (a, &A0, m);
    rightShiftPolynomial_spX (b, &B0, m);

    if (d < HALFGCD_CROSSOVER) {
        freeMat4_spX (&E);

        iterHalfGCDMatrixInForm_FFT_spX_inp (&A0, &B0, M, d, Pptr);

        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);

        return;
    }

    dusmat4_t* RR = NULL;
    long d1 = (d+1) >> 1;

    // corner cases:
    if (d1 >= d) {
        d1 = d-1;
    } else if (d1 < 1) {
        d1 = 1;
    }

    halfGCDMatrixInForm_FFT_spX (A0, B0, &RR, d1, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (RR, &A0, &B0, Pptr);

    if (isZero_spX (B0)) {
        
        freePolynomial_spX (&A0);
        freeMat4_spX (&E);

        *M = RR;
        return;
    }


    long d2 = degPolynomial_spX (B0) + m + d - degPolynomial_spX (a);
    
    if (d2 <= 0) {
        
        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);
        freeMat4_spX (&E);

        *M = RR;
        return;
    }
    // else do the second part of HGCD:

    duspoly_t* q = NULL;
    dusmat4_t* SS = NULL;

    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL;
    
    // plainDivPolynomialsInForm_spX_inp (&A0, B0, &q, Pptr);
    ////////////////////
    fastDivPolynomialInForm_wPSInv_FFT_spX (A0, B0, &tmp1, &q, Pptr);
    freePolynomial_spX (&A0);
    A0 = tmp1;
    tmp1 = NULL;
    ///////////////////
    swap_spX (&A0, &B0);

    halfGCDMatrixInForm_FFT_spX (A0, B0, &SS, d2, Pptr);

    //plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    fastMulPolynomialInForm_FFT_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    fastMulPolynomialInForm_FFT_spX (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[1], &RR->polys[3]);
    swap_spX (&RR->polys[3], &tmp2); 

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_FFT_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS); 
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 
 
    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);

}


void extHalfGCDMatrixInForm_FFT_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eHGCDM-FFT\n");
#endif

    dusmat4_t* E = NULL;
    identityMat4_spX (&E, Pptr);

    if (isZero_spX (a) || isZero_spX (b)) {
        *M = E;
        return;
    }

    if (degPolynomial_spX (a) - degPolynomial_spX (b) >= d ||
        degPolynomial_spX (a) - degPolynomial_spX (b) <  0 ) {
        *M = E;
        return;
    }

    if (d < HALFGCD_CROSSOVER) {
        freeMat4_spX (&E);
        iterHalfGCDMatrixInForm_FFT_spX (a, b, M, d, Pptr);
        return;
    }

    long d1 = (d+1) >> 1;
    if (d1 == d) {
        d1 = d-1;
    } else if (d1 < 1) {
        d1 = 1;
    }

    dusmat4_t* RR = NULL;
    
    duspoly_t* A0 = deepCopyPolynomial_spX (a);
    duspoly_t* B0 = deepCopyPolynomial_spX (b);

    halfGCDMatrixInForm_FFT_spX (a, b, &RR, d1, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (RR, &A0, &B0, Pptr);

    if (isZero_spX (B0)) {
        freeMat4_spX (&E);
        freePolynomial_spX (&A0);
     
        *M = RR;
        return;
    }

    long d2 = degPolynomial_spX (B0) - degPolynomial_spX (a) + d;
    if (d2 <= 0) {
        freeMat4_spX (&E);
        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);

        *M = RR;
        return;
    }

    duspoly_t* q = NULL;
    dusmat4_t* SS = NULL;

    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL;

    // plainDivPolynomialsInForm_spX_inp (&A0, B0, &q, Pptr);
    //////////////////
    fastDivPolynomialInForm_wPSInv_FFT_spX (A0, B0, &tmp1, &q, Pptr);
    freePolynomial_spX (&A0);
    A0 = tmp1;
    tmp1 = NULL;
    //////////////////
    swap_spX (&A0, &B0);

    extHalfGCDMatrixInForm_FFT_spX (A0, B0, &SS, d2, Pptr);

    // plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    fastMulPolynomialInForm_FFT_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    fastMulPolynomialInForm_FFT_spX  (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[3], &RR->polys[1]);
    swap_spX (&RR->polys[3], &tmp2);

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_FFT_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS);
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);
}

void extHalfGCDMatrixInForm_FFT_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eHGCDMin-FFT\n");
#endif

    dusmat4_t* E = NULL;
    identityMat4_spX (&E, Pptr);

    duspoly_t* A = *a;
    duspoly_t* B = *b;

    if (isZero_spX (A) || isZero_spX (B)) {
        *M = E;
        return;
    }

    if (degPolynomial_spX (A) - degPolynomial_spX (B) >= d ||
        degPolynomial_spX (A) - degPolynomial_spX (B) <  0 ) {
        *M = E;
        return;
    }

    if (d < HALFGCD_CROSSOVER) {
        freeMat4_spX (&E);
        iterHalfGCDMatrixInForm_FFT_spX (A, B, M, d, Pptr);

        mulMat4ToVecInForm_FFT_spX_inp (*M, a, b, Pptr);
        
        return;
    }

    long d1 = (d+1) >> 1;
    if (d1 == d) {
        d1 = d-1;
    } else if (d1 < 1) {
        d1 = 1;
    }

    polysize_t degA = degPolynomial_spX (A);
    dusmat4_t* RR = NULL;
    
    // duspoly_t* A0 = A; // deepCopyPolynomial_spX (a);
    // duspoly_t* B0 = B; // deepCopyPolynomial_spX (b);

    halfGCDMatrixInForm_FFT_spX (A, B, &RR, d1, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (RR, &A, &B, Pptr);

    if (isZero_spX (B)) {
        freeMat4_spX (&E);

        *a = A; // TODO: delete after test!
        *b = NULL; 
    
        *M = RR;
        return;
    }

    long d2 = degPolynomial_spX (B) - degA + d;
    if (d2 <= 0) {
        freeMat4_spX (&E);

        *a = A;
        *b = B;

        *M = RR;
        return;
    }

    duspoly_t* q = NULL;
    dusmat4_t* SS = NULL;

    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL;

    //plainDivPolynomialsInForm_spX_inp (&A, B, &q, Pptr);
    //////////////////
    fastDivPolynomialInForm_wPSInv_FFT_spX (A, B, &tmp1, &q, Pptr);
    freePolynomial_spX (&A);
    A = tmp1;
    tmp1 = NULL;
    /////////////////
    swap_spX (&A, &B);

    extHalfGCDMatrixInForm_FFT_spX_inp (&A, &B, &SS, d2, Pptr);

    // plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    fastMulPolynomialInForm_FFT_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    fastMulPolynomialInForm_FFT_spX (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[3], &RR->polys[1]);
    swap_spX (&RR->polys[3], &tmp2);

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_FFT_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS);
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 

    *a = A; // TODO: better assignmnet... then delete these lines!
    *b = B;

    // freePolynomial_spX (&A0);
    // freePolynomial_spX (&B0);
}

void halfGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** ap, duspoly_t** bp, Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCD-FFT\n");
#endif

    if (isZero_spX (b)) {
        *ap = deepCopyPolynomial_spX (a);
        *bp = deepCopyPolynomial_spX (b);

        return;
    }

    duspoly_t* A; 
    duspoly_t* B;
    if (degPolynomial_spX (a) >= degPolynomial_spX (b)) {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    } else {
        A = deepCopyPolynomial_spX (b);
        B = deepCopyPolynomial_spX (a);
    }

    polysize_t degA = degPolynomial_spX (A);
    long d = (degPolynomial_spX (A)+1) >> 1;

    if (degPolynomial_spX (A) - degPolynomial_spX (B) >= d) {
        *ap = A;
        *bp = B; 

        return;
    }

    long d1 = (d+1) >> 1;
    if (d1 == d) {
        d1 = d-1;
    } else if (d1 < 1) {
        d1 = 1;
    }

    dusmat4_t* M;
    halfGCDMatrixInForm_FFT_spX (A, B, &M, d1, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (M, &A, &B, Pptr);
    freeMat4_spX (&M);

    if (isZero_spX (B)) {
        *ap = A;
        *bp = NULL;

        return;
    }

    long d2 = degPolynomial_spX (B) - degA + d;

    if (d2 <= 0) {
        *ap = A;
        *bp = B;
        
        return;
    }

    duspoly_t* q = NULL;
    duspoly_t* tmp1 = NULL;

    // plainDivPolynomials_spX_inp (&A, B, &q, Pptr);
    ///////////////////
    fastDivPolynomialInForm_wPSInv_FFT_spX (A, B, &tmp1, &q, Pptr);
    freePolynomial_spX (&A);
    A = tmp1;
    ///////////////////
    swap_spX (&A, &B);

    freePolynomial_spX (&q);

    halfGCDMatrixInForm_FFT_spX (A, B, &M, d2, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (M, &A, &B, Pptr);
    freeMat4_spX (&M);

    *ap = A;
    *bp = B; 
}

void halfGCDInForm_FFT_spX_inp (duspoly_t** a, duspoly_t** b, Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDin-FFT\n");
#endif

    if (isZero_spX (*b)) {
        return;
    }

    duspoly_t* A; 
    duspoly_t* B;


    if (degPolynomial_spX (*a) >= degPolynomial_spX (*b)) {
        A = *a;
        B = *b;
    } else {
        A = *b;
        B = *a;
    }

    polysize_t degA = degPolynomial_spX (A);
    long d = (degA+1) >> 1;

    if (degA - degPolynomial_spX (B) >= d) {
        return;
    }

    long d1 = (d+1) >> 1;
    if (d1 == d) {
        d1 = d-1;
    } else if (d1 < 1) {
        d1 = 1;
    }

    dusmat4_t* M;
    halfGCDMatrixInForm_FFT_spX (A, B, &M, d1, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (M, &A, &B, Pptr);
    freeMat4_spX (&M);

    if (isZero_spX (B)) {
        *a = A;
        *b = B;

        return;
    }

    long d2 = degPolynomial_spX (B) - degA + d;

    if (d2 <= 0) {
        *a = A;
        *b = B;

        return;
    }

    duspoly_t* q = NULL;
    duspoly_t* tmp1 = NULL;

    // plainDivPolynomials_spX_inp (&A, B, &q, Pptr);
    ///////////////////
    fastDivPolynomialInForm_wPSInv_FFT_spX (A, B, &tmp1, &q, Pptr);
    freePolynomial_spX (&A);
    A = tmp1;
    ///////////////////
    swap_spX (&A, &B);

    freePolynomial_spX (&q);

    halfGCDMatrixInForm_FFT_spX (A, B, &M, d2, Pptr);
    mulMat4ToVecInForm_FFT_spX_inp (M, &A, &B, Pptr);
    freeMat4_spX (&M);

    *a = A;
    *b = B; 
}


void GCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** g, Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "GCD-FFT\n");
#endif

    if (degPolynomial_spX (a) < PLAINGCD_CROSSOVER ||
        degPolynomial_spX (b) < PLAINGCD_CROSSOVER) {
        plainGCDInForm_spX (a, b, g, Pptr);
        return;
    } 

    duspoly_t* A = NULL;
    duspoly_t* B = NULL;
    
    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL;

    if (degPolynomial_spX (a) == degPolynomial_spX (b)) {
        //plainRemPolynomialsInForm_spX (a, b, &B, Pptr);
        ////////////////////
        fastDivPolynomialInForm_wPSInv_FFT_spX (a, b, &B, &tmp1, Pptr);
        freePolynomial_spX (&tmp1);
        ////////////////////
        A = deepCopyPolynomial_spX (b);
    } else if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    } else {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    }

    while (!isZero_spX (B) && degPolynomial_spX (A) > PLAINGCD_CROSSOVER) {
        halfGCDInForm_FFT_spX_inp (&A, &B, Pptr);

        if (!isZero_spX (B)) {

            //plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
            ////////////////////////
            fastDivPolynomialInForm_wPSInv_FFT_spX (A, B, &tmp1, &tmp2, Pptr);
            freePolynomial_spX (&tmp2);
            freePolynomial_spX (&A);
            A = tmp1;
            ////////////////////////
            swap_spX (&A, &B);
        }
    }

    plainGCDInForm_spX (A, B, g, Pptr);

    freePolynomial_spX (&A);
    freePolynomial_spX (&B);
}


void extGCDInForm_FFT_spX (duspoly_t* a, duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, Prime_ptr* Pptr)
{


#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eGCD-FFT\n");
#endif

    if (degPolynomial_spX (a) < PLAINGCD_CROSSOVER ||
        degPolynomial_spX (b) < PLAINGCD_CROSSOVER) {
        plainExtGCDInForm_spX (a, b, u, v, g, Pptr);
        return;
    }

    duspoly_t* A = NULL;
    duspoly_t* B = NULL;
    duspoly_t* q = NULL;

    int isEqualDeg = 0;
    int isSwap = 0;

    if (degPolynomial_spX (a) == degPolynomial_spX (b)) {
        // plainDivPolynomialsInForm_spX (a, b, &B, &q, Pptr);
        fastDivPolynomialInForm_wPSInv_FFT_spX (a, b, &B, &q, Pptr);
        A = deepCopyPolynomial_spX (b);
        
        isEqualDeg = 1;

    } else if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
        A = deepCopyPolynomial_spX (b);
        B = deepCopyPolynomial_spX (a);

        isSwap = 1;

    } else {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    }
 
    long d = degPolynomial_spX (A) + 1;
    
    dusmat4_t* M = NULL;
    
    duspoly_t* uu = NULL;
    duspoly_t* vv = NULL;
    duspoly_t* gg = NULL;

    duspoly_t* tmp = NULL;

    // fprintf(stderr, "d = $l\n", d);

    extHalfGCDMatrixInForm_FFT_spX_inp (&A, &B, &M, d, Pptr);

    freePolynomial_spX (&B);

    if (isEqualDeg) {
        uu = M->polys[1];
        
        //plainMulPolynomialsInForm_spX (q, M->polys[1], &tmp, Pptr);
        fastMulPolynomialInForm_FFT_spX (q, M->polys[1], &tmp, Pptr);
        subPolynomialsInForm_spX (M->polys[0], tmp, &vv, Pptr);
        freePolynomial_spX (&tmp);
        freePolynomial_spX (&q);
    } else if (isSwap) {
        uu = M->polys[1];
        vv = M->polys[0];
    } else {
        uu = M->polys[0];
        vv = M->polys[1];
    }

    freePolynomial_spX (&M->polys[2]);
    freePolynomial_spX (&M->polys[3]);

    elem_t lc_A;

    monicPolynomialInForm_spX (A, &gg, &lc_A, Pptr);
    
    freePolynomial_spX (&A);

    if (lc_A != smallprimefield_convert_in (1, Pptr)) {
        lc_A = smallprimefield_inv (lc_A, Pptr);
        
        scalarMulPolynomialInForm_spX_inp (&uu, lc_A, Pptr);
        scalarMulPolynomialInForm_spX_inp (&vv, lc_A, Pptr);
    }

    *g = gg;
    *u = uu;
    *v = vv;

}
