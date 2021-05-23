 

#include "ModularPolynomial/DUSP_NTL_Support.h" 

#if defined(WITH_NTL) && WITH_NTL

using namespace NTL;
using namespace std;

////////////////
/// Conversions 
////////////////
void convertToNTLPolyZZX_spX (const duspoly_t* a, ZZX& f)
{
	if (isZero_spX (a)) {
		f = 0;
		return;
	}

	elem_t* elems = a->elems;

	for (polysize_t i = 0; i <= POLY_LT(a); i++) {
		SetCoeff (f, i, elems[i]);
	}

	f.normalize();
}

void convertToNTLPolyZZpX_spX (const duspoly_t* a, ZZ_pX& f, const Prime_ptr* Pptr)
{
    ZZ p = NTL::ZZ(Pptr->prime);
    ZZ_p::init(p);
    
    if (isZero_spX (a)) {
      SetCoeff(f, 0, ZZ_p(0));
      return;
    }
	
    for (long i = 0; i <= POLY_LT(a); i++) {
    	ZZ_p elm = ZZ_p (a->elems[i]);
    	SetCoeff (f, i, elm);
    }

    f.normalize();
}


void convertFromNTLPolyZZX_spX (duspoly_t** a, const ZZX& f)
{
	if (deg(f) == -1) {
		*a = NULL;
		return;
	}

	duspoly_t* tmp = makePolynomial_spX (deg(f)+1);
	tmp->lt = deg(f);

	elem_t coef;

	for (polysize_t i = 0; i <= deg(f); i++) {
		convertFromNTLZZ_spX (&coef, coeff(f, i));
		tmp->elems[i] = coef;
	}

	normalize_spX (&tmp);
	*a = tmp;
}

void convertFromNTLPolyZZpX_spX (duspoly_t** a, const ZZ_pX& f, const Prime_ptr* Pptr)
{
	if (deg(f) == -1) {
		*a = NULL;
		return;
	}

	duspoly_t* tmp = makePolynomial_spX (deg(f)+1);
	tmp->lt = deg(f);

	for (polysize_t i = 0; i <= deg(f); i++) {
		tmp->elems[i] = conv<long>(coeff (f, i));
	}

	normalize_spX (&tmp);
	*a = tmp;

}

/////////////
// Mul / Sqr
/////////////

void ntlMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (a) < MUL_NTLvsBPAS_THRESHOLD || 
		degPolynomial_spX (b) < MUL_NTLvsBPAS_THRESHOLD) {
		mulPolynomialsInForm_spX (a, b, c, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX fg;
	mul (fg, f, g);

	convertFromNTLPolyZZpX_spX (c, fg, Pptr);

	convertPolynomialToMontgomery_spX_inp (c, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlMulPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (*a) < MUL_NTLvsBPAS_THRESHOLD || 
		degPolynomial_spX (b) < MUL_NTLvsBPAS_THRESHOLD) {
		mulPolynomialsInForm_spX_inp (a, b, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (*a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX fg;
	mul (fg, f, g);

	freePolynomial_spX (a);
	convertFromNTLPolyZZpX_spX (a, fg, Pptr);

	convertPolynomialToMontgomery_spX_inp (a, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlSqrPolynomialInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (a) < SQR_NTLvsBPAS_THRESHOLD) {
		sqrPolynomialInForm_spX (a, a2, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);

	ZZ_pX f;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);

	ZZ_pX f2;
	sqr (f2, f);

	convertFromNTLPolyZZpX_spX (a2, f2, Pptr);

	convertPolynomialToMontgomery_spX_inp (a2, Pptr);

	freePolynomial_spX (&ao);	
}


void ntlSSMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (a) < MUL_NTLvsBPAS_THRESHOLD || 
		degPolynomial_spX (b) < MUL_NTLvsBPAS_THRESHOLD) {
		mulPolynomialsInForm_spX (a, b, c, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZX f, g;
	convertToNTLPolyZZX_spX (ao, f);
	convertToNTLPolyZZX_spX (bo, g);

	ZZX fg;
	SSMul (fg, f, g);

	convertFromNTLPolyZZX_spX (c, fg);

	convertPolynomialToMontgomery_spX_inp (c, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlSSMulPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (*a) < MUL_NTLvsBPAS_THRESHOLD || 
		degPolynomial_spX (b) < MUL_NTLvsBPAS_THRESHOLD) {
		mulPolynomialsInForm_spX_inp (a, b, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (*a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZX f, g;
	convertToNTLPolyZZX_spX (ao, f);
	convertToNTLPolyZZX_spX (bo, g);

	ZZX fg;
	SSMul (fg, f, g);

	freePolynomial_spX (a);
	convertFromNTLPolyZZX_spX (a, fg);

	convertPolynomialToMontgomery_spX_inp (a, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlSSSqrPolynomialInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (a) < SQR_NTLvsBPAS_THRESHOLD) {
		sqrPolynomialInForm_spX (a, a2, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);

	ZZX f;
	convertToNTLPolyZZX_spX (ao, f);

	ZZX f2;
	SSSqr (f2, f);

	convertFromNTLPolyZZX_spX (a2, f2);

	convertPolynomialToMontgomery_spX_inp (a2, Pptr);

	freePolynomial_spX (&ao);	
}


///////////////////
// Div / Rem / Quo
///////////////////
void ntlDivPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, duspoly_t** q, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (a) < DIV_NTLvsBPAS_THRESHOLD ||
		degPolynomial_spX (b) < DIV_NTLvsBPAS_THRESHOLD) {
		plainDivPolynomialsInForm_spX (a, b, r, q, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX rr, qq;
	DivRem (qq, rr, f, g);

	convertFromNTLPolyZZpX_spX (r, rr, Pptr);
	convertFromNTLPolyZZpX_spX (q, qq, Pptr);

	convertPolynomialToMontgomery_spX_inp (r, Pptr);
	convertPolynomialToMontgomery_spX_inp (q, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlDivPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (*a) < DIV_NTLvsBPAS_THRESHOLD ||
		degPolynomial_spX (b) < DIV_NTLvsBPAS_THRESHOLD) {
		plainDivPolynomialsInForm_spX_inp (a, b, q, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (*a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX rr, qq;
	DivRem (qq, rr, f, g);

	freePolynomial_spX (a);
	convertFromNTLPolyZZpX_spX (a, qq, Pptr);
	convertFromNTLPolyZZpX_spX (q, qq, Pptr);

	convertPolynomialToMontgomery_spX_inp (a, Pptr);
	convertPolynomialToMontgomery_spX_inp (q, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlRemPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (a) < DIV_NTLvsBPAS_THRESHOLD ||
		degPolynomial_spX (b) < DIV_NTLvsBPAS_THRESHOLD) {
		plainRemPolynomialsInForm_spX (a, b, r, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX rr;
	rem (rr, f, g);

	convertFromNTLPolyZZpX_spX (r, rr, Pptr);

	convertPolynomialToMontgomery_spX_inp (r, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);

}

void ntlRemPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
	if (degPolynomial_spX (*a) < DIV_NTLvsBPAS_THRESHOLD ||
		degPolynomial_spX (b) < DIV_NTLvsBPAS_THRESHOLD) {
		plainRemPolynomialsInForm_spX_inp (a, b, Pptr);
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (*a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX rr;
	rem (rr, f, g);

	freePolynomial_spX (a);
	convertFromNTLPolyZZpX_spX (a, rr, Pptr);	

	convertPolynomialToMontgomery_spX_inp (a, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlQuoPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr)
{

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, g;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, g, Pptr);

	ZZ_pX qq;
	div (qq, f, g);

	convertFromNTLPolyZZpX_spX (q, qq, Pptr);

	convertPolynomialToMontgomery_spX_inp (q, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}


///////////////
// GCD / ExtGCD 
///////////////

void ntlGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, const Prime_ptr* Pptr)
{

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, h;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, h, Pptr);

	ZZ_pX gg;
	GCD (gg, f, h);

	convertFromNTLPolyZZpX_spX (g, gg, Pptr);

	convertPolynomialToMontgomery_spX_inp (g, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}

void ntlExtGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, const Prime_ptr* Pptr)
{

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
	duspoly_t* bo = convertPolynomialFromMontgomery_spX (b, Pptr);

	ZZ_pX f, h;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);
	convertToNTLPolyZZpX_spX (bo, h, Pptr);

	ZZ_pX gg, uu, vv;
	XGCD (gg, uu, vv, f, h);

	convertFromNTLPolyZZpX_spX (g, gg, Pptr);
	convertFromNTLPolyZZpX_spX (u, uu, Pptr);
	convertFromNTLPolyZZpX_spX (v, vv, Pptr);

	convertPolynomialToMontgomery_spX_inp (g, Pptr);
	convertPolynomialToMontgomery_spX_inp (u, Pptr);
	convertPolynomialToMontgomery_spX_inp (v, Pptr);

	freePolynomial_spX (&ao);
	freePolynomial_spX (&bo);
}


/////////////
// Factoring
/////////////

void ntlFactoringInForm_spX (const duspoly_t* a, factors_t** facts, const Prime_ptr* Pptr)
{
	if (isZero_spX (a)) {
		*facts = NULL;
		return;
	}

	duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);

	if (leadingCoeffInForm_spX (ao) != 1) {
		fprintf(stderr, "DUSP Error: In ntlFactoringInForm_spX, the input polynomial must be monic!\n");
		exit (1);
	}

	vec_pair_ZZ_pX_long factPair;
	ZZ_pX f;
	convertToNTLPolyZZpX_spX (ao, f, Pptr);

	berlekamp (factPair, f);

	factors_t* ffacts = initFactors_spX (factPair.length());
	for (int i = 0; i < factPair.length(); i++) {
		convertFromNTLPolyZZpX_spX (&(ffacts->polys[i]), factPair[i].a, Pptr);
		convertPolynomialToMontgomery_spX_inp (&(ffacts->polys[i]), Pptr);
		ffacts->exps[i] = factPair[i].b;
	}

	*facts = ffacts;
	freePolynomial_spX (&ao);
}

void squareFreeFactorizationInFormll1_spX (const duspoly_t* a, factsll_t** F, const Prime_ptr* Pptr) 
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
		ntlGCDInForm_spX (aa, ap, &g, Pptr); 
		plainDivPolynomialsInForm_spX (aa, g, &r, &q, Pptr);

		freePolynomial_spX (&r);
		freePolynomial_spX (&ap);

		if (degPolynomial_spX (q) > 0) {
			do {

				ntlGCDInForm_spX (q, g, &wg, Pptr);
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


void exponentiatePolynomialInForm1_spX (const duspoly_t* a, polysize_t n, duspoly_t** an, const Prime_ptr* Pptr)
{
    if (an == NULL){
        return;
    }

    if (n == 0) {
        elem_t one = smallprimefield_convert_in (1, Pptr);
        *an = constPolynomial_spX (one, Pptr);
        return;
    } else if (n == 1) {
        *an = deepCopyPolynomial_spX (a);
        return;
    }

    if (n < EXP_NTLvsBPAS_THRESHOLD) {
	    duspoly_t* mul = NULL;
	    duspoly_t* tmp;
	    duspoly_t* b = deepCopyPolynomial_spX (a);

	    long nn =  (long) n;

	    while (nn > 1) {
	        if (nn & 1) {            
	            if (isZero_spX (mul)) {
	                mul = deepCopyPolynomial_spX (b);
	            } else {
	                ntlSSMulPolynomialsInForm_spX_inp (&mul, b, Pptr);
	            }
	        }
	        
	        ntlSSSqrPolynomialInForm_spX (b, &tmp, Pptr);

	        freePolynomial_spX (&b);
	        b = tmp; tmp = NULL;

	        nn >>= 1;
	    }

	    if (isZero_spX (mul)) {
	        mul = deepCopyPolynomial_spX (b);
	    } else {
	        ntlSSMulPolynomialsInForm_spX_inp (&mul, b, Pptr);
	    }


	    freePolynomial_spX (&b);
	    *an = mul;
    } else {
		duspoly_t* ao = convertPolynomialFromMontgomery_spX (a, Pptr);
		
		ZZ_pX f;
		convertToNTLPolyZZpX_spX (ao, f, Pptr);
		
		ZZ_pX fn;
		power (fn, f, n);

		convertFromNTLPolyZZpX_spX (an, fn, Pptr);

		convertPolynomialToMontgomery_spX_inp (an, Pptr);

		freePolynomial_spX (&ao);
    }
}


void modExponentiatePolynomialInForm1_spX (const duspoly_t* a, polysize_t n, const duspoly_t* f, duspoly_t** an, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "ModExpPoly\n");
#endif

    if (an == NULL){
        return;
    }

    if (n == 0) {
        elem_t one = smallprimefield_convert_in (1, Pptr);
        *an = constPolynomial_spX (one, Pptr);
        return;
    } else if (n == 1) {
        ntlRemPolynomialsInForm_spX (a, f, an, Pptr);
        // *an = deepCopyPolynomial_spX (a);
        return;
    }

    polysize_t deg_f = degPolynomial_spX (f);
    if (deg_f < 1) {
        fprintf(stderr, "DUSP Error: Bad Input for modulus polynomial, f, in modExponentiatePolynomialInForm.\n");
        exit(1);
    }

    duspoly_t* mul = NULL;
    duspoly_t* tmp;
    duspoly_t* b = deepCopyPolynomial_spX (a);

    long nn =  (long) n;

    while (nn > 1) {
        if (nn & 1) {            
            if (isZero_spX (mul)) {
                mul = deepCopyPolynomial_spX (b);
            } else {
                // mulPolynomialsInForm_spX_inp (&mul, b, Pptr);
                ntlSSMulPolynomialsInForm_spX_inp (&mul, b, Pptr);
            }
        }
        
        // sqrPolynomialInForm_spX (b, &tmp, Pptr);
        ntlSqrPolynomialInForm_spX (b, &tmp, Pptr);

        if (POLY_LT(tmp) > 2*deg_f) {
            // plainRemPolynomialsInForm_spX_inp (&tmp, f, Pptr);
            ntlRemPolynomialsInForm_spX_inp (&tmp, f, Pptr);
        }

        freePolynomial_spX (&b);

        b = tmp; tmp = NULL;

        nn >>= 1;
    }

    if (isZero_spX (mul)) {
        mul = deepCopyPolynomial_spX (b);
    } else {
        // mulPolynomialsInForm_spX_inp (&mul, b, Pptr);
        ntlSSMulPolynomialsInForm_spX_inp (&mul, b, Pptr);
    }

    // plainRemPolynomialsInForm_spX_inp (&mul, f, Pptr);
	ntlRemPolynomialsInForm_spX_inp (&mul, f, Pptr);

    freePolynomial_spX (&b);
    *an = mul;
}



duspoly_t* expXModSubInForm1_spX (polysize_t q, const duspoly_t* f, duspoly_t** h, const Prime_ptr* Pptr) 
{
	if (q == 0 || isZero_spX (f)) {
		return NULL;
	}

	duspoly_t* xpoly = Xpolynomial_spX ();
	duspoly_t* xpow = NULL;
	duspoly_t* res  = NULL;

	if (isZero_spX (*h)) {
		xpow = makePowerOfXpolynomial_spX (q, Pptr);
		ntlRemPolynomialsInForm_spX_inp (&xpow, f, Pptr);
		
		if (h != NULL) {
			*h = xpow;
		}
	} else {
		modExponentiatePolynomialInForm1_spX (*h, q, f, &xpow, Pptr);
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

void distinctDegFactorizationInFormll1_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr)
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

	duspoly_t* hmx = expXModSubInForm1_spX (Pptr->prime, f, &h, Pptr);
	factsll_t* nodeG = initFactsll_spX (s);

	ntlGCDInForm_spX (hmx, f, &(nodeG->poly), Pptr);
	plainDivPolynomialsInForm_spX (f, nodeG->poly, &t, &f1, Pptr);
	
	freePolynomial_spX (&hmx);
	freePolynomial_spX (&t);

	factsll_t* head = nodeG;
	factsll_t* tail = nodeG;

	duspoly_t* fs = NULL;

	while (f1 != NULL && f1->lt > 0) {
		s++;
		hmx = expXModSubInForm1_spX (Pptr->prime, f, &h, Pptr);
		nodeG = initFactsll_spX (s);

		ntlGCDInForm_spX (hmx, f1, &(nodeG->poly), Pptr);
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

void modPoweringForEDF1_spX (const duspoly_t* a, polysize_t d, const duspoly_t* f, duspoly_t** mul_pows, const Prime_ptr* Pptr)
{

	if (mul_pows == NULL) {
		return;
	}

	if (d == 0) {
		*mul_pows = constPolynomialInForm_spX (1, Pptr);
		return;
	}


	polysize_t p = Pptr->prime;
	polysize_t n = (p-1)>>1;

	if (n < 0) {
		fprintf(stderr, "DUSP Error: Bad prime in modPoweringForEDF.\n");
		exit(1);
	}

	if (d == 1) {
		modExponentiatePolynomialInForm1_spX (a, n, f, mul_pows, Pptr);
		return;
	}

	// compute a**q, a**q**2, a**q**3, ..., a**q**(d-1) 
	// and multiply all together in muls

	duspoly_t* muls = NULL;
	duspoly_t* aqi  = NULL;
	duspoly_t* aqim = NULL;

	// compute aqim = a**p mod f 
	modExponentiatePolynomialInForm1_spX (a, p, f, &muls, Pptr);
	
	if (d > 2) {
		aqim = deepCopyPolynomial_spX (muls);
	}

	ntlMulPolynomialsInForm_spX_inp (&muls, a, Pptr);
	ntlRemPolynomialsInForm_spX_inp (&muls, f, Pptr); // muls = (a * a**p) (mod f)


	for (polysize_t i = 2; i < d; i++) {
		// compute aqi = aqim**p = a**p**(i) mod f
		modExponentiatePolynomialInForm1_spX (aqim, p, f, &aqi, Pptr);

		ntlMulPolynomialsInForm_spX_inp (&muls, aqi, Pptr);
		ntlRemPolynomialsInForm_spX_inp (&muls, f, Pptr);

		freePolynomial_spX (&aqim);
		aqim = aqi;
		aqi = NULL;
	}

	freePolynomial_spX (&aqim);

	modExponentiatePolynomialInForm1_spX (muls, n, f, &aqi, Pptr);
	freePolynomial_spX (&muls);
	
	*mul_pows = aqi;
}


int equalDegSplittingInForm1_spX (const duspoly_t* f, polysize_t d, duspoly_t** g, const Prime_ptr* Pptr)
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
	ntlGCDInForm_spX (f, a, &g1, Pptr);

	if (g1 != NULL && g1->lt > 0) {
		*g = g1;

		freePolynomial_spX (&a);
		return 1;
	} else {
		freePolynomial_spX (&g1);
	}

	duspoly_t* b  = NULL;
	duspoly_t* bb = NULL;

	modPoweringForEDF1_spX (a, d, f, &b, Pptr);

	duspoly_t* one = constPolynomial_spX (smallprimefield_convert_in (1, Pptr), Pptr);

	subPolynomialsInForm_spX (b, one, &bb, Pptr);
	
	freePolynomial_spX (&a);
	freePolynomial_spX (&b);
	freePolynomial_spX (&one);

	ntlGCDInForm_spX (bb, f, &g1, Pptr);
	
	freePolynomial_spX (&bb); 

	if (g1 != NULL && g1->lt > 0 && !isEqual_spX (f, g1)) {
		*g = g1;

		return 1;
	} else {
		freePolynomial_spX (&g1);
	}

	return 0;

}

void equalDegFactorizationInFormll1_spX (const duspoly_t* f, polysize_t d, factsll_t** G, const Prime_ptr* Pptr)
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
		isFact = equalDegSplittingInForm1_spX (f, d, &g, Pptr);
	}

	duspoly_t* r = NULL;
	duspoly_t* q = NULL;

	plainDivPolynomialsInForm_spX (f, g, &r, &q, Pptr);
	freePolynomial_spX (&r);

	equalDegFactorizationInFormll1_spX (g, d, G, Pptr);
	equalDegFactorizationInFormll1_spX (q, d, G, Pptr);

	freePolynomial_spX (&g);
	freePolynomial_spX (&q);
}

void equalDegFactorsInFormll1_spX (factsll_t* D, factsll_t** G, const Prime_ptr* Pptr)
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

		equalDegFactorizationInFormll1_spX (Dcur->poly, Dcur->exp, &EDFo, Pptr);

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

void modFactorizationInFormll1_spX (const duspoly_t* f, factsll_t** G, const Prime_ptr* Pptr)
{

	if (G == NULL) {
		return;
	}	

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

	while (fi != NULL && fi->lt > 0) {

		// fprintf(stderr, " i = %ld\n", i); // TEST

		i += 1;

		hmx = expXModSubInForm1_spX (Pptr->prime, f, &hi, Pptr);

		ntlGCDInForm_spX (hmx, fi, &g, Pptr);

		freePolynomial_spX (&hmx);

		if (g == NULL || g->lt < 1) {
			
			freePolynomial_spX (&g);

			continue;
		}

		// fprintf(stderr, "EDF Starts...\n");
		equalDegFactorizationInFormll1_spX (g, i, &Gg, Pptr);
		// fprintf(stderr, "EDF Stops...\n");
		
		freePolynomial_spX (&g);

		cur = Gg;
		while (cur != NULL) {
			r = deepCopyPolynomial_spX (cur->poly);
			
			isDiv = isDividablePolysQuoInForm_spX (fi, r, &q, Pptr);

			freePolynomial_spX (&fi);
			fi = q; q = NULL;
			e = 1;	
		
			isDiv = isDividablePolysQuoInForm_spX (fi, r, &q, Pptr);

			while (isDiv) {
				e += 1;

				freePolynomial_spX (&fi);
				fi = q; q = NULL;	

				isDiv = isDividablePolysQuoInForm_spX (fi, r, &q, Pptr);
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

	freePolynomial_spX (&hi);
	freePolynomial_spX (&fi);

	*G = head;
}

int modFactorVerificationInFormll1_spX (const duspoly_t* f, factsll_t* G, const Prime_ptr* Pptr) 
{
	if (isZero_spX (f)) {
		if (isZeroFactsll_spX (G)) 
			return 1;
		else 
			return 0;
	} 

	factsll_t* cur = G;
	duspoly_t* ff = constPolynomialInForm_spX (1, Pptr);
	duspoly_t* tmp;
	duspoly_t* tmp2;

	while (cur != NULL) {
		if (cur->poly != NULL) {
			if (cur->exp == 0) {
				fprintf(stderr, "DUSP Error, In modFactorVerificationInFormll: G couldn't have polynomial with idx=0 as a factor!\n");
				exit (1);
			}

			exponentiatePolynomialInForm1_spX (cur->poly, cur->exp, &tmp, Pptr);
			// plainMulPolynomialsInForm_spX (ff, tmp, &tmp2, Pptr);
			// mulPolynomialsInForm_spX (ff, tmp, &tmp2, Pptr);
			ntlMulPolynomialsInForm_spX (ff, tmp, &tmp2, Pptr);

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

#endif
