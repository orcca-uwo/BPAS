#include "ModularPolynomial/DUSP_Support.h"
#include "ModularPolynomial/DUSP_Support_Test.h"


void setCoefs_spX (duspoly_t** a, elem_t* coefs, polysize_t sz, const Prime_ptr* Pptr)
{       
    // if (!isZero_spX (*a)) {
	// 	freePolynomial_spX (&*a);
    // }
    if (a == NULL) {
        return;
    }

    duspoly_t* poly = makePolynomial_spX (sz);
    POLY_LT(poly) = sz - 1;
    for (polysize_t i = 0; i < sz; i++) {
        poly->elems[i] = smallprimefield_convert_in (coefs[i], Pptr);
    }

    *a = poly;

    normalize_spX (a);
}


void setCoefsInForm_spX (duspoly_t** a, elem_t* coefs, polysize_t sz)
{       
    if (a == NULL) {
        return;
    }

    duspoly_t* poly = makePolynomial_spX (sz);
    POLY_LT(poly) = sz - 1;
    for (polysize_t i = 0; i < sz; i++) {
        poly->elems[i] = coefs[i];
    }

    *a = poly;

    normalize_spX (a);
}


elem_t* getCoefs_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
		return 0;
    }
    
    elem_t* cp_elems = (elem_t*) malloc (sizeof(elem_t)*POLY_LT(a));
    for (polysize_t i = 0; i < POLY_ALLOC(a); ++i) {
		cp_elems[i] =  smallprimefield_convert_out (a->elems[i], Pptr);
		// smallprimefield_convert_out (a->elems[i], Pptr);
    }
    
    return cp_elems;    
}

elem_t* getCoefsInForm_spX (const duspoly_t* a)
{
    if (isZero_spX (a)) {
		return 0;
    }
    
    elem_t* cp_elems = (elem_t*) malloc (sizeof(elem_t)*POLY_LT(a));
    memcpy (cp_elems, a->elems, sizeof(elem_t)*POLY_LT(a));
  //   for (polysize_t i = 0; i < POLY_ALLOC(a); ++i) {
		// cp_elems[i] = a->elems[i];
  //   }
    
    return cp_elems;    
}


elem_t* getCoefsInForm_spX_inp (duspoly_t* a)
{
    if (isZero_spX (a)) {
		return 0;
    }

    return a->elems;    
}


/**********************
 * Convert Polynomial
 *********************/

duspoly_t* convertPolynomialToMontgomery_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
		return NULL;
    }
    
    duspoly_t* poly = makePolynomial_spX (POLY_ALLOC (a));
    POLY_LT(poly) = POLY_LT(a);
    
    const elem_t* elems = a->elems;
    polysize_t i;
    
    for (i = 0; i <= POLY_LT(a); i++) {
		poly->elems[i] = smallprimefield_convert_in (elems[i], Pptr);
	    // smallprimefield_convert_in (elems[i], Pptr);
    }
    
    normalize_spX (&poly);
    return poly;
}

void convertPolynomialToMontgomery_spX_inp (duspoly_t** a, const Prime_ptr* Pptr)
{
    if (isZero_spX (*a)) {
        freePolynomial_spX (a);
        return;
    }
    
    elem_t* elems = (*a)->elems;
    polysize_t i;
    
    for (i = 0; i <= POLY_LT(*a); i++) {
        (*a)->elems[i] = smallprimefield_convert_in (elems[i], Pptr);
        // smallprimefield_convert_in (elems[i], Pptr);
    }
    
    normalize_spX (a);
}


duspoly_t* convertPolynomialFromMontgomery_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
		return NULL;
    }
    
    duspoly_t* poly = makePolynomial_spX (POLY_ALLOC (a));
    POLY_LT(poly) = POLY_LT(a);
    
    const elem_t* elems = a->elems;
    polysize_t i;
    
    for (i = 0; i <= POLY_LT(a); i++) {
		poly->elems[i] = smallprimefield_convert_out (elems[i], Pptr); // smallprimefield_convert_out (elems[i], Pptr);
    }
    
    normalize_spX (&poly);
    return poly;
}


void convertPolynomialFromMontgomery_spX_inp (duspoly_t** a, const Prime_ptr* Pptr)
{
    if (isZero_spX (*a)) {
        freePolynomial_spX (a);
        return;
    }
    
    elem_t* elems = (*a)->elems;
    polysize_t i;
    
    for (i = 0; i <= POLY_LT(*a); i++) {
        (*a)->elems[i] = smallprimefield_convert_out (elems[i], Pptr);
    }
    
    normalize_spX (a);
}

/*************************
 * Shift/Split Polynomial
 *************************/

void rightShiftPolynomial_spX (const duspoly_t* a , duspoly_t** sa, polysize_t n)
{
    if (sa == NULL) {
        return;
    }

    if (isZero_spX (a)) {
		*sa = NULL;
		return;
    }

    if (n == 0) {
		*sa = deepCopyPolynomial_spX (a);
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    
    if (deg_a < n) {
		*sa = NULL;
		return;
    }
    
    duspoly_t* rsa = makePolynomial_spX (POLY_ALLOC (a) - n);
    POLY_LT (rsa) = deg_a - n;
    
    memmove (rsa->elems, a->elems + n, (POLY_LT(rsa)+1)*sizeof(elem_t));
    
    normalize_spX (&rsa);
    *sa = rsa;
}

void rightShiftPolynomial_spX_inp (duspoly_t** a, polysize_t n)
{    
    duspoly_t* A = *a;

    if (isZero_spX (A)) {
        return;
    }

    if (n == 0) {
        return;
    }
    
    polysize_t deg_a = degPolynomial_spX (A);
    
    if (deg_a < n) {
        freePolynomial_spX (&*a);
        return;
    }
        
    elem_t* Elems = (elem_t*) malloc (sizeof(elem_t)*(POLY_ALLOC(A)-n));

    memmove (Elems, A->elems + n, sizeof(elem_t)*(POLY_ALLOC(A)-n));
    
    free(A->elems);
    A->elems = Elems;

    normalize_spX (&A);
}

void leftShiftPolynomial_spX (const duspoly_t* a , duspoly_t** sa, polysize_t n)
{
    if (sa == NULL) {
        return;
    }

    if (isZero_spX (a)) {
		*sa = NULL;
		return;
    }
    
    if (n == 0) {
		*sa = deepCopyPolynomial_spX (a);
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    
    
    duspoly_t* lsa = makePolynomial_spX (deg_a + n + 1);
    POLY_LT (lsa) = deg_a + n;    
    memmove (lsa->elems + n, a->elems, (POLY_LT(a)+1)*sizeof(elem_t));


    normalize_spX (&lsa);
    *sa = lsa;
}

void leftShiftAddPolynomials_spX_inp (duspoly_t** a, const duspoly_t* b, polysize_t n, const Prime_ptr* Pptr)
{
    if (isZero_spX (b)) {
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (*a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t dmax  = MAX_spX(deg_a, deg_b + n);
    polysize_t i = 0;
    
    reallocatePolynomial_spX (a, dmax+1);

    elem_t* a_elems = (*a)->elems;
    
    for (i = 0; i <= deg_b; ++i) {
		a_elems[i+n] =  smallprimefield_add (a_elems[i+n], b->elems[i], Pptr);
    }
    
    normalize_spX (a);
}


void leftShiftSubPolynomials_spX_inp (duspoly_t** a, const duspoly_t* b, polysize_t n, const Prime_ptr* Pptr)
{    
    if (isZero_spX (b)) {
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (*a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t dmax  = MAX_spX(deg_a, deg_b + n);
    polysize_t i = 0;
    
    reallocatePolynomial_spX (a, dmax+1);
    
    elem_t* a_elems = (*a)->elems;
    
    for (i = 0; i <= deg_b; ++i) {
		a_elems[i+n] = smallprimefield_sub (a_elems[i+n], b->elems[i], Pptr);
    }
    
    normalize_spX (a);
}


void leftSplitPolynomial_spX (const duspoly_t* a, duspoly_t** sa, polysize_t n)
{
    if (sa == NULL) {
        return;
    }

    if (!isZero_spX (*sa)) {
        freePolynomial_spX (&*sa);
    }

    if (isZero_spX (a)) {
		return;
    }

    if (n < 1) {
		*sa = NULL;
		return;
    }

    duspoly_t* lsp = makePolynomial_spX (n);
    POLY_LT (lsp) = n-1;

    for (polysize_t i = 0; i < n && i <= POLY_LT(a) ; ++i) {
		(lsp)->elems[i] = a->elems[i];
    }
    
    normalize_spX (&lsp);
    *sa = lsp;
}


/***********************************
 * Basic Polynoimal Func/Arithmetic 
 **********************************/

void addPolynomials_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL) {
        return;
    }

    if (!isZero_spX (*c)) {
        freePolynomial_spX (&*c);
    }

    // poly + 0 = poly
    if (isZero_spX (a)) {
		*c = deepCopyPolynomial_spX(b);
		return;
    } else if (isZero_spX (b)){
		*c = deepCopyPolynomial_spX(a);
		return;
    }

    polysize_t deg_a = degPolynomial_spX(a);
    polysize_t deg_b = degPolynomial_spX(b);
    polysize_t minab = MIN_spX(deg_a, deg_b);
    polysize_t maxab = MAX_spX(deg_a, deg_b);
    register polysize_t i;

    duspoly_t* cc = makePolynomial_spX (maxab + 1);
    POLY_LT(cc) = maxab;
    
    const elem_t* a_elems = a->elems;
    const elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;

    // poly + elem = 
    if (maxab == 0) {
		cc->elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
		*c = cc;
		return;
    } else if (minab == 0) {
		if (deg_a == 0) {
			for (i = 0; i <= deg_b; i++) {
				c_elems[i] = b_elems[i];
			}
			c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		} else {
			for (i = 0; i <= deg_a; i++) {
				c_elems[i] = a_elems[i];
			}
			c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		}
    }

    // poly + poly = 
    for (i = 0; i <= minab; i++) {
		c_elems[i] = smallprimefield_add(a_elems[i], b_elems[i], Pptr);
	
    }
    if (deg_a > minab && !isEqual_spX(cc, a)) {
		for (i = minab+1; i <= deg_a; i++) {
			c_elems[i] = a_elems[i];
		}
    } else if (deg_b > minab && !isEqual_spX(cc, b)) {
		for (i = minab+1; i <= deg_b; i++) {
			c_elems[i] = b_elems[i];
		}
    } else {
		normalize_spX (&cc);
    }

    if (isZero_spX (cc)) {
        freePolynomial_spX (&cc);
		*c = NULL;
		return;
    }

    *c = cc;
}

void addPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{    
    if (c == NULL) {
        return;
    }

    if (!isZero_spX (*c)) {
        freePolynomial_spX (&*c);
    }

    // poly + 0 = poly
    if (isZero_spX (a)) {
		*c = deepCopyPolynomial_spX(b);
		return;
    } else if (isZero_spX (b)){
		*c = deepCopyPolynomial_spX(a);
		return;
    }

    polysize_t deg_a = degPolynomial_spX(a);
    polysize_t deg_b = degPolynomial_spX(b);
    polysize_t minab = MIN_spX(deg_a, deg_b);
    polysize_t maxab = MAX_spX(deg_a, deg_b);
    register polysize_t i;

    duspoly_t* cc = makePolynomial_spX (maxab+1);
    POLY_LT(cc) = maxab; // upper bound 
    
    const elem_t* a_elems = a->elems;
    const elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;

    // poly + elem = 
    if (maxab == 0) {
		cc->elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
		*c = cc;
		return;
    } else if (minab == 0) {
		if (deg_a == 0) {
			for (i = 0; i <= deg_b; i++) {
				c_elems[i] = b_elems[i];
			}
			c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		} else {
			for (i = 0; i <= deg_a; i++) {
				c_elems[i] = a_elems[i];
			}
			c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		}
    }

    // poly + poly = 
    for (i = 0; i <= minab; i++) {
		c_elems[i] = smallprimefield_add(a_elems[i], b_elems[i], Pptr);
	
    }
    if (deg_a > minab && !isEqual_spX(cc, a)) {
		for (i = minab+1; i <= deg_a; i++) {
			c_elems[i] = a_elems[i];
		}
    } else if (deg_b > minab && !isEqual_spX(cc, b)) {
		for (i = minab+1; i <= deg_b; i++) {
			c_elems[i] = b_elems[i];
		}
    } else {
		normalize_spX (&cc);
    }

    if (isZero_spX (cc)) {
        freePolynomial_spX (&cc);
		*c = NULL;
		return;
    }

    *c = cc;
}

// Addition in-place w.r.t the first input
void addPolynomials_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
    if (isZero_spX (b)) {
		return;
    }
    
    duspoly_t* c = NULL;
    duspoly_t* ap = *a;
    
    addPolynomials_spX (ap, b, &c, Pptr);
    
    freePolynomial_spX (&ap);
    *a = c;
}

// Addition in-place w.r.t the first input
void addPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
    // poly + 0 = poly
    if (isZero_spX (*a)) {
        freePolynomial_spX (a);
        *a = deepCopyPolynomial_spX(b);
        return;
    } else if (isZero_spX (b)){
        return;
    }

    polysize_t deg_a = degPolynomial_spX(*a);
    polysize_t deg_b = degPolynomial_spX(b);
    polysize_t minab = MIN_spX(deg_a, deg_b);
    polysize_t maxab = MAX_spX(deg_a, deg_b);
    register polysize_t i;


    duspoly_t* cc = *a;
    reallocatePolynomial_spX (&cc, maxab + 1);
    cc->lt = maxab;

    elem_t* a_elems = cc->elems;
    const elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;

    // poly + elem = 
    if (maxab == 0) {
        c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
        *a = cc;
        return;
    } else if (minab == 0) {
        if (deg_a == 0) {
            for (i = 1; i <= deg_b; i++) {
                c_elems[i] = b_elems[i];
            }
            c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
            return;
        } else {
            c_elems[0] = smallprimefield_add(a_elems[0], b_elems[0], Pptr);
            return;
        }
    }

    // poly + poly = 
    if (deg_b < maxab) {
        maxab = deg_b;
    }

    for (i = 0; i <= maxab; i++) {
        c_elems[i] = smallprimefield_add(a_elems[i], b_elems[i], Pptr);
    }

    normalize_spX (a);

   if (isZero_spX (cc)) {
        freePolynomial_spX (a);
        return;
    }

}


void subPolynomials_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL) {
        return;
    }

    polysize_t deg_a = degPolynomial_spX(a);
    polysize_t deg_b = degPolynomial_spX(b);
    polysize_t minab = MIN_spX(deg_a, deg_b);
    polysize_t maxab = MAX_spX(deg_a, deg_b);
    register polysize_t i;
    
    // poly - 0 = poly && 0 - poly = -poly
    if (isZero_spX (a)) {
		if (isZero_spX (b)) {
			*c = NULL;
		}
		duspoly_t* tmp = deepCopyPolynomial_spX (b);
		negPolynomial_spX (tmp, Pptr);
		*c = tmp;
		return;
    } else if (isZero_spX (b)){
		printPolynomial_spX (a, Pptr);
		*c = deepCopyPolynomial_spX(a);
		return; 
    }
    
    duspoly_t* cc = makePolynomial_spX (maxab+1);
    POLY_LT(cc) = maxab;

    const elem_t* a_elems = a->elems;
    const elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;
    
    // poly - elem =  && elem - poly = 
    if (maxab == 0) {
		c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
		*c = cc;
		return;
    } else if (minab == 0) {
		if (deg_a == 0) {
			for (i = 0; i <= deg_b; i++) {
				c_elems[i] = smallprimefield_sub (0, b_elems[i], Pptr);
			}
			c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		} else {
			for (i = 0; i <= deg_a; i++) {
				c_elems[i] = a_elems[i];
			}
			c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		}
    }

    // poly - poly = 
    for (i = 0; i <= minab; i++) {
		c_elems[i] = smallprimefield_sub(a_elems[i], b_elems[i], Pptr);
    }
    if (deg_a > minab && !isEqual_spX(cc, a)) {
		for (i = minab+1; i <= deg_a; i++) {
			c_elems[i] = a_elems[i];
		}
    } else if (deg_b > minab && !isEqual_spX(cc, b)) {
		for (i = minab+1; i <= deg_b; i++) {
			c_elems[i] = smallprimefield_sub (0, b_elems[i], Pptr);
		}
    } else {
		normalize_spX (&cc);
    }
    

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
		*c = NULL;
		return;
    }
    
    *c = cc;
    return;
}


void subPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL){
        return;
    }

    polysize_t deg_a = degPolynomial_spX(a);
    polysize_t deg_b = degPolynomial_spX(b);
    polysize_t minab = MIN_spX(deg_a, deg_b);
    polysize_t maxab = MAX_spX(deg_a, deg_b);
    register polysize_t i;
    
    // poly - 0 = poly && 0 - poly = -poly
    if (isZero_spX (a)) {
		if (isZero_spX (b)) {
			*c = NULL;
		}
		duspoly_t* tmp = deepCopyPolynomial_spX (b);
		negPolynomial_spX (tmp, Pptr);
		*c = tmp;
		return;
    } else if (isZero_spX (b)){
		*c = deepCopyPolynomial_spX(a);
		return; 
    }
    
    duspoly_t* cc = makePolynomial_spX (maxab+1);
    POLY_LT(cc) = maxab;

    const elem_t* a_elems = a->elems;
    const elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;
    
    // poly - elem =  && elem - poly = 
    if (maxab == 0) {
		c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
		*c = cc;
		return;
    } else if (minab == 0) {
		if (deg_a == 0) {
			for (i = 0; i <= deg_b; i++) {
				c_elems[i] = smallprimefield_sub (0, b_elems[i], Pptr);
			}
			c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		} else {
			for (i = 0; i <= deg_a; i++) {
				c_elems[i] = a_elems[i];
			}
			c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
			*c = cc;
			return;
		}
    }

    // poly - poly = 
    for (i = 0; i <= minab; i++) {
		c_elems[i] = smallprimefield_sub(a_elems[i], b_elems[i], Pptr);
    }
    if (deg_a > minab && !isEqual_spX(cc, a)) {
		for (i = minab+1; i <= deg_a; i++) {
			c_elems[i] = a_elems[i];
		}
    } else if (deg_b > minab && !isEqual_spX(cc, b)) {
		for (i = minab+1; i <= deg_b; i++) {
			c_elems[i] = smallprimefield_sub (0, b_elems[i], Pptr);
		}
    } else {
		normalize_spX (&cc);
    }
    

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
		*c = NULL;
		return;
    }
    
    *c = cc;
    return;
}


// Subtraction in-place w.r.t the first input
void subPolynomials_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
    if (isZero_spX (b)) {
		return;
    }
    
    duspoly_t* c = NULL;
    duspoly_t* ap = *a;
    
    subPolynomials_spX (ap, b, &c, Pptr);
    
    freePolynomial_spX (&ap);
    *a = c;
}


// Subtraction in-place w.r.t the first input
void subPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
    if (a == NULL) {
        return;
    }

    polysize_t deg_a = degPolynomial_spX(*a);
    polysize_t deg_b = degPolynomial_spX(b);
    polysize_t minab = MIN_spX(deg_a, deg_b);
    polysize_t maxab = MAX_spX(deg_a, deg_b);
    register polysize_t i;
    
    // poly - 0 = poly && 0 - poly = -poly
    if (isZero_spX (*a)) {
        if (isZero_spX (b)) {
            freePolynomial_spX (a);
            return;
        }
        *a = deepCopyPolynomial_spX (b);
        negPolynomial_spX (*a, Pptr);
        return;
    } else if (isZero_spX (b)){
        return; 
    }
    
    duspoly_t* cc = *a;
    reallocatePolynomial_spX (&cc, maxab+1);
    cc->lt = maxab;

    const elem_t* a_elems = cc->elems;
    const elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;
    
    // poly - elem =  && elem - poly = 
    if (maxab == 0) {
        c_elems[0] = smallprimefield_sub(a_elems[0], b_elems[0], Pptr);
        return;
    } else if (minab == 0) {
        if (deg_a == 0) {
            for (i = 1; i <= deg_b; i++) {
                c_elems[i] = smallprimefield_sub (0, b_elems[i], Pptr);
            }
            c_elems[0] = smallprimefield_sub (a_elems[0], b_elems[0], Pptr);

            setLT_spX (&cc, 1);
            return;
        } else {
            c_elems[0] = smallprimefield_sub (a_elems[0], b_elems[0], Pptr);
            return;
        }
    }

    if (deg_b < maxab) {
        maxab = deg_b;
    }

    // poly - poly = 
    for (i = 0; i <= maxab; i++) {
        c_elems[i] = smallprimefield_sub(a_elems[i], b_elems[i], Pptr);
    }

    normalize_spX (a);
    // setLT_spX (&cc, 1);

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
        freePolynomial_spX (a);
        return;
    }
    
    return;
}


void negPolynomial_spX (duspoly_t* a, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
		return;
    }
    
    elem_t* elems = a->elems;
    polysize_t i;
    for (i = 0; i <= POLY_LT(a); i++) {
		elems[i] = smallprimefield_sub (0, elems[i], Pptr);
		// - elems[i] + (((-elems[i]) >> PRIME_BIT_LENGTH) & pr);
    }    
    
    normalize_spX (&a);
}
 
 
void negPolynomialInForm_spX (duspoly_t* a, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
		return;
    }
    
    elem_t* elems = a->elems;
    polysize_t i;
    for (i = 0; i <= POLY_LT(a); i++) {
        elems[i] = smallprimefield_sub (0, elems[i], Pptr);
	    //- elems[i] + (((-elems[i]) >> PRIME_BIT_LENGTH) & pr);
    }    

    normalize_spX (&a);
}


void scalarMulPolynomial_spX  (const duspoly_t* a, elem_t b, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL){
        return;
    }

    if (isZero_spX(a)) {
		*c = zeroPolynomial_spX();
		return; 
    }
    

    elem_t* elems = a->elems;
    elem_t elem_in = 0;
    elem_t b_in = smallprimefield_convert_in (b, Pptr);
    
    polysize_t i;
    
    duspoly_t* cc = makePolynomial_spX (POLY_ALLOC (a));
    POLY_LT(cc) = POLY_LT(a); 

    for (i = 0; i <= POLY_LT (cc); i++){
		elem_in = smallprimefield_convert_in (elems[i], Pptr);
		cc->elems[i] = smallprimefield_convert_out (smallprimefield_mul
													(elem_in, b_in, Pptr), Pptr);
    }
    
    normalize_spX (&cc);
    *c = cc;
}

void scalarMulPolynomialInForm_spX  (const duspoly_t* a, elem_t b, duspoly_t** c, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "sMul\n");
#endif

    if (c == NULL){
        return;
    }

    if (isZero_spX(a)) {
		*c = zeroPolynomial_spX();
		return; 
    }
    

    elem_t* elems = a->elems;
    
    polysize_t i;
    
    duspoly_t* cc = makePolynomial_spX (POLY_ALLOC (a));
    POLY_LT(cc) = POLY_LT(a);
    
    for (i = 0; i <= POLY_LT (cc); i++){
		cc->elems[i] = smallprimefield_mul (elems[i], b, Pptr);
    }
    
    normalize_spX (&cc);
    *c = cc;
}

void scalarMulPolynomialInForm_spX_inp (duspoly_t** a, elem_t b, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "sMulin\n");
#endif

    if (isZero_spX (*a)) {
        return;
    }

    elem_t* elems = (*a)->elems;
    for (int i = 0; i <= POLY_LT(*a); ++i) {
        elems[i] = smallprimefield_mul (elems[i], b, Pptr);
    }

	normalize_spX (a);
}


void monicPolynomial_spX (const duspoly_t* a, duspoly_t** ma, elem_t* lc, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
	   if (ma != NULL){
        	*ma = constPolynomialInForm_spX (1, Pptr);
       }
       if (lc != NULL){
    		*lc = 1;
       }
		return;
    }

    elem_t lc_a = leadingCoeffInForm_spX (a);
    
    if (lc != NULL) {
        *lc = lc_a;
    }
    
    if (ma != NULL) {
        lc_a = smallprimefield_convert_in (lc_a, Pptr);
    
        if (lc_a == 0) {
    		fprintf (stderr, "DUSP Error, leading coefficient of non-zero polynomial cannot be zero!\n");
    		exit(1);
        } else {    
    		lc_a = smallprimefield_inv (lc_a, Pptr);
    		lc_a = smallprimefield_convert_out (lc_a, Pptr);

    		scalarMulPolynomial_spX (a, lc_a, ma, Pptr);
        }
        
        normalize_spX (ma);
    }
}

void monicPolynomialInForm_spX (const duspoly_t* a, duspoly_t** ma, elem_t* lc, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
        if (ma != NULL){
    		*ma = constPolynomial_spX (1, Pptr);
        }
	    if (lc != NULL){
        	*lc = smallprimefield_convert_in (1, Pptr);
        }
		return;
    }

    elem_t lc_a = leadingCoeffInForm_spX (a);
    if (lc != NULL){
        *lc = lc_a;
    }
    
    if (ma != NULL){
        if (lc_a == 0) {
    		fprintf (stderr, "DUSP Error, leading coefficient of non-zero polynomial cannot be zero!\n");
    		exit(1);
        } else {    
    		lc_a = smallprimefield_inv (lc_a, Pptr);      
    		scalarMulPolynomialInForm_spX (a, lc_a, ma, Pptr);
        }
        
        normalize_spX (ma);
    }
}

void monicPolynomialInForm_spX_inp (duspoly_t** a, elem_t* lc, const Prime_ptr* Pptr)
{

    if (isZero_spX (*a)) {
        *a = constPolynomial_spX (1, Pptr);
        if (lc != NULL){
            *lc = smallprimefield_convert_in (1, Pptr);
        }
        return;
    }

    elem_t lc_a = leadingCoeffInForm_spX (*a);

    if (lc != NULL){
        *lc = lc_a;
    }
    
    if (lc_a == 0) {
        fprintf (stderr, "DUSP Error, leading coefficient of non-zero polynomial cannot be zero!\n");
        exit(1);
    } else {    
        lc_a = smallprimefield_inv (lc_a, Pptr);      
        scalarMulPolynomialInForm_spX_inp (a, lc_a, Pptr);
    }
    
    normalize_spX (a);
}

void plainMulPolynomials_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL) {
        return;
    }

    // poly * 0 = 0
    if (isZero_spX (a) || isZero_spX (b)){
		*c = zeroPolynomial_spX ();
		return;
    }

    // TODO : isEqual ?
    
    // poly * elem = 
    if (isConstant_spX (a)) {
		scalarMulPolynomial_spX (b, a->elems[0], c, Pptr);
        return; 
    }
    if (isConstant_spX (b)) {
		scalarMulPolynomial_spX (a, b->elems[0], c, Pptr);
		return; 
    }
				  
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_c = deg_a + deg_b;
    polysize_t i, j;
    polysize_t jmin = 0;
    polysize_t jmax = 0;
    
    duspoly_t* cc = makePolynomial_spX (deg_c + 1);
    POLY_LT (cc) = deg_c;

    elem_t part_sum = 0;
    elem_t a_in = 0;
    elem_t b_in = 0;
    
    elem_t* a_elems = a->elems;
    elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;

    // poly * poly = 
    for (i = 0; i <= deg_c; i++) {
		part_sum = 0;
		jmin = MIN_spX(i, deg_a);
		jmax = MAX_spX(i-deg_b, 0);
		for (j = jmax; j <= jmin; j++) {
			// part_sum += a[j]*b[i-j];
			a_in = smallprimefield_convert_in (a_elems[j], Pptr);
			b_in = smallprimefield_convert_in (b_elems[i-j], Pptr);    
			part_sum = smallprimefield_add ( part_sum, smallprimefield_mul (a_in, b_in, Pptr), Pptr);
		}
		c_elems[i] = smallprimefield_convert_out (part_sum, Pptr); // part_sum % pr; 
    }

    normalize_spX (&cc);

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
		*c = NULL;
		return;
    }

    *c = cc;
}


void plainMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pMul\n");
#endif

    if (c == NULL) {
        return;
    }

    // poly * 0 = 0
    if (isZero_spX (a) || isZero_spX (b)){
		*c = zeroPolynomial_spX ();
		return;
    }

    // poly * elem = 
    if (isConstant_spX (a)) {
		scalarMulPolynomialInForm_spX (b, a->elems[0], c, Pptr);
        return; 
    }
    if (isConstant_spX (b)) {
		scalarMulPolynomialInForm_spX (a, b->elems[0], c, Pptr);
		return; 
    }
				  
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_c = deg_a + deg_b;
    polysize_t i, j;
    polysize_t jmin = 0;
    polysize_t jmax = 0;

    
    duspoly_t* cc = makePolynomial_spX (deg_c + 1);
    POLY_LT (cc) = deg_c;
    elem_t part_sum = 0;
    
    elem_t* a_elems = a->elems;
    elem_t* b_elems = b->elems;
    elem_t* c_elems = cc->elems;

    // poly * poly = 
    for (i = 0; i <= deg_c; i++) {
		part_sum = 0;
		jmin = MIN_spX(i, deg_a);
		jmax = MAX_spX(i-deg_b, 0);
		for (j = jmax; j <= jmin; j++) {
			// part_sum += a[j]*b[i-j];
			part_sum = smallprimefield_add (part_sum, smallprimefield_mul (a_elems[j], b_elems[i-j], Pptr), Pptr);
		}
		c_elems[i] = part_sum; 
    }

    normalize_spX (&cc);

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
		*c = NULL;
		return;
    }

    *c = cc;
}

duspoly_t* plainMulManyPolynomialsInForm_spX(duspoly_t const*const* polys, int n, const Prime_ptr* Pptr) {
    if (n == 0 || polys == NULL) {
        return NULL;
    }
    if (n < 2) {
        return deepCopyPolynomial_spX(polys[0]);
    }

    duspoly_t* prod;
    plainMulPolynomialsInForm_spX(polys[0], polys[1], &prod, Pptr);
    duspoly_t* tmpProd;
    for (int i = 2; i < n; ++i) {
        plainMulPolynomialsInForm_spX(prod, polys[i], &tmpProd, Pptr);
        freePolynomial_spX(&prod);
        prod = tmpProd;
    }

    return prod;
}

void multiplyByBinomialInForm_spX_inp(duspoly_t** a_ptr, elem_t b, const Prime_ptr* Pptr) {
    if (a_ptr == NULL) {
        return;
    }

    duspoly_t* a = *a_ptr;
    if (a == NULL) {
        a = makePolynomial_spX(2);
        a->elems[0] = smallprimefield_convert_in(1, Pptr);
        a->elems[1] = b;
        a->lt = 1;
        *a_ptr = a;
        return;
    }

    if (a->alloc < a->lt + 1) {
        reallocatePolynomial_spX(a_ptr, a->lt*2); //ammortize resize cost hopefully;
    }

    elem_t tmp = smallprimefield_mul(a->elems[0], b, Pptr);
    tmp = smallprimefield_mul(smallprimefield_convert_in(-1, Pptr), tmp, Pptr);

    for (int i = 0; i < a->lt; ++i) {
        //a[i] -= a[i+1]*b
        a->elems[i] = smallprimefield_sub(a->elems[i],smallprimefield_mul(a->elems[i+1], b, Pptr), Pptr);
    }
    
    //shift right 1
    memmove(a->elems+1, a->elems, sizeof(elem_t)*(a->lt+1));
    a->elems[0] = tmp;

    ++(a->lt);
}

void plainSqrPolynomial_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr)
{
    if (a2 == NULL) {
        return;
    }

    // 0*0 = 0
    if (isZero_spX (a)){
		*a2 = zeroPolynomial_spX ();
		return; 
    }
    
    // term*term = 
    if (isConstant_spX (a)) {
		*a2 = constPolynomialInForm_spX ((a->elems[0] * a->elems[0]), Pptr);
        return; 
    }
				  
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_c = deg_a << 1;
    polysize_t i, j;
    polysize_t jmin = 0;
    polysize_t jmax = 0;
    polysize_t m = 0;
    
    duspoly_t* cc = makePolynomial_spX (deg_c + 1);
    POLY_LT (cc) = deg_c;

    elem_t part_sum = 0;
    elem_t a_in_l = 0;
    elem_t a_in_r = 0;
    
    elem_t* a_elems = a->elems;
    elem_t* c_elems = cc->elems;
    
    for (i = 0; i <= deg_c; i++) {
		part_sum = 0;
		jmin = MAX_spX(i-deg_a, 0);
		jmax = MIN_spX(i, deg_a);
		m = jmax - jmin + 1;
		jmax = jmin + (m >> 1) - 1;
		for (j = jmin; j <= jmax; j++) {
			a_in_l = smallprimefield_convert_in (a_elems[j], Pptr);
			a_in_r = smallprimefield_convert_in (a_elems[i-j], Pptr);
			part_sum = smallprimefield_add (part_sum, smallprimefield_mul (a_in_l,  a_in_r, Pptr), Pptr);
		}
		part_sum = smallprimefield_add (part_sum, part_sum, Pptr);
		if (m & 1) {
			a_in_l = smallprimefield_convert_in (a_elems[jmax+1], Pptr);
			part_sum = smallprimefield_add (part_sum, smallprimefield_mul (a_in_l, a_in_l, Pptr), Pptr);
		}
		c_elems[i] = smallprimefield_convert_out (part_sum, Pptr);// part_sum % pr; 
    }

    normalize_spX (&cc);

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
		*a2 = NULL;
		return;
    }
    
    *a2 = cc;
}


void plainSqrPolynomialInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pSqr\n");
#endif

    if (a2 == NULL) {
        return;
    }
 
    // 0*0 = 0
    if (isZero_spX (a)){
		*a2 = zeroPolynomial_spX ();
		return; 
    }
    
    // term*term = 
    if (isConstant_spX (a)) {
		*a2 = constPolynomialInForm_spX (smallprimefield_mul
										 (a->elems[0], a->elems[0], Pptr), Pptr);
        return; 
    }
				  
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_c = deg_a << 1;
    polysize_t i, j;
    polysize_t jmin = 0;
    polysize_t jmax = 0;
    polysize_t m = 0;
    
    duspoly_t* cc = makePolynomial_spX (deg_c + 1);
    POLY_LT (cc) = deg_c;

    elem_t part_sum = 0;
    
    elem_t* a_elems = a->elems;
    elem_t* c_elems = cc->elems;
    
    for (i = 0; i <= deg_c; i++) {

		part_sum = 0;
		jmin = MAX_spX(i-deg_a, 0); 
		jmax = MIN_spX(i, deg_a);   
		m = jmax - jmin + 1;    
		jmax = jmin + (m >> 1) - 1; 

		for (j = jmin; j <= jmax; j++) { 

			part_sum = smallprimefield_add (part_sum, smallprimefield_mul (a_elems[j],  a_elems[i-j], Pptr), Pptr);
		}
		part_sum = smallprimefield_add (part_sum, part_sum, Pptr); // part_sum = 2(a[0]*a[4] + a[1]*a[3])
		if (m & 1) { 

			part_sum = smallprimefield_add (part_sum, smallprimefield_mul (a_elems[jmax+1], a_elems[jmax+1], Pptr), Pptr); 
		}
		c_elems[i] = part_sum; // part_sum % Pptr->prime; // replace to MontRed in BigPrimeField
    }

    normalize_spX (&cc);

    if (cc == NULL || POLY_ALLOC(cc) == 0) {
		*a2 = NULL;
		return;
    }
    
    *a2 = cc;
}

void plainDivPolynomials_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, duspoly_t** q, const Prime_ptr* Pptr)
{
    
    // a / 0 = Error
    if (isZero_spX (b)) {
    	fprintf (stderr, "DUSP Error: Division by zero!\n");
    	exit (1);
    }
    
    // 0 / b = 0
    if (isZero_spX (a)) {
        if (r != NULL){
        	*r = NULL;
        } 
        if (q != NULL){
        	*q = NULL;
        }
    	return;
    }

    // a/a = 1
    if (isEqual_spX (a, b)) {
        if (r != NULL){
		  *r = NULL;
        }
        if (q != NULL){
		  *q = constPolynomialInForm_spX (1, Pptr);
        }
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
        if (r != NULL){
        	*r = deepCopyPolynomial_spX (a);
        }
        if (q != NULL){
        	*q = NULL;
        }
    	return;
    }

    polysize_t diff_deg = deg_a - deg_b;
    
    duspoly_t* rr = deepCopyPolynomial_spX (a); // rr = a;
    duspoly_t* qq = makePolynomial_spX (diff_deg+1);
    POLY_LT(qq) = diff_deg;
    
    polysize_t i, j;
    elem_t lc_b_in = leadingCoeffInForm_spX (b);
    lc_b_in = smallprimefield_convert_in (lc_b_in, Pptr);
    lc_b_in = smallprimefield_inv (lc_b_in, Pptr);
    
    elem_t lc_rr_in = 0;
    elem_t b_in = 0;    
    
    // a = b*qq + rr
    for (i = diff_deg; i >= 0; i--) {	
		lc_rr_in = idxCoeffInForm_spX (rr, deg_b + i);
		lc_rr_in = smallprimefield_convert_in (lc_rr_in, Pptr);

		if (lc_rr_in) {
			lc_rr_in = smallprimefield_mul (lc_rr_in, lc_b_in, Pptr);

			qq->elems[i] = smallprimefield_convert_out (lc_rr_in, Pptr); 
	    
			rr->elems[deg_b+i] = 0;
			rr->lt -= 1;
			for (j = deg_b-1;j >= 0; j--) {		   
				b_in = smallprimefield_convert_in (b->elems[j], Pptr);
				b_in = smallprimefield_mul (b_in, lc_rr_in, Pptr);
				b_in = smallprimefield_convert_out (b_in, Pptr);
		
				rr->elems[i+j] = smallprimefield_sub (rr->elems[i+j], b_in, Pptr);
			}
		}	
    }


    if (q != NULL) {
        normalize_spX (&qq);
        if (isZero_spX (qq)) {
    		freePolynomial_spX (&qq);
    		*q = constPolynomialInForm_spX (1, Pptr);
        } else {
    		*q = qq;
        }        
    } else {
        freePolynomial_spX (&qq);
    }
    
    if (r != NULL){
        normalize_spX (&rr);
        *r = rr;
    } else {
        freePolynomial_spX (&rr);
    }
}

void plainDivPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, duspoly_t** q, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pDiv\n");
#endif


    // a / 0 = Error
    if (isZero_spX (b)) {
    	fprintf (stderr, "DUSP Error: Division by zero!\n");
    	exit (1);
    }
    
    // 0 / b = 0
    if (isZero_spX (a)) {
        if (r != NULL){
            *r = NULL;
        } 
        if (q != NULL){
            *q = NULL;
        }
        return;
    }

    // a/a = 1
    if (isEqual_spX (a, b)) {
        if (r != NULL){
		  *r = NULL;
        }
        if (q != NULL){
		  *q = constPolynomial_spX (1, Pptr); // TODO : 
        }
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
        if (r != NULL){
    	   *r = deepCopyPolynomial_spX (a);
        }
        if (q != NULL){
    	   *q = NULL;
        }
    	return;
    }

    polysize_t diff_deg = deg_a - deg_b;
    
    duspoly_t* rr = deepCopyPolynomial_spX (a); // rr = a;
    duspoly_t* qq = makePolynomial_spX (diff_deg+1);
    POLY_LT(qq) = diff_deg;
    
    polysize_t i, j;
    elem_t lc_b_in = leadingCoeffInForm_spX (b);
    /* lc_b_in = smallprimefield_convert_in (lc_b_in, Pptr); */
    lc_b_in = smallprimefield_inv (lc_b_in, Pptr);
    
    elem_t lc_rr_in = 0;
    elem_t b_in = 0;    
    
    // a = b*qq + rr
    for (i = diff_deg; i >= 0; i--) {	
		lc_rr_in = idxCoeffInForm_spX (rr, deg_b + i);
		/* lc_rr_in = smallprimefield_convert_in (lc_rr_in, Pptr); */

		if (lc_rr_in) {
			lc_rr_in = smallprimefield_mul (lc_rr_in, lc_b_in, Pptr);

			qq->elems[i] = lc_rr_in; // smallprimefield_convert_out (lc_rr_in, Pptr); 
	    
			rr->elems[deg_b+i] = 0;
			rr->lt -= 1;
			for (j = deg_b-1;j >= 0; j--) {		   
				/* b_in = smallprimefield_convert_in (b->elems[j], Pptr); */
				b_in = smallprimefield_mul (b->elems[j], lc_rr_in, Pptr);
				/* b_in = smallprimefield_convert_out (b_in, Pptr); */
		
				rr->elems[i+j] = smallprimefield_sub (rr->elems[i+j], b_in, Pptr);
			}
		}	
    }

     if (q != NULL) {
        normalize_spX (&qq);
        if (isZero_spX (qq)) {
            freePolynomial_spX (&qq);
            *q = constPolynomial_spX (1, Pptr);
        } else {
            *q = qq;
        }        
    } else {
        freePolynomial_spX (&qq);
    }
    
    if (r != NULL){
        normalize_spX (&rr);
        *r = rr;
    } else {
        freePolynomial_spX (&rr);
    }
}

void plainDivPolynomials_spX_inp ( duspoly_t** a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr)
{
    duspoly_t* r = NULL;
    plainDivPolynomials_spX (*a, b, &r, q, Pptr);

    freePolynomial_spX (a);
    *a = r; 
}


void plainDivPolynomialsInForm_spX_inp ( duspoly_t** a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pDivin\n");
#endif

    duspoly_t* r = NULL;
    plainDivPolynomialsInForm_spX (*a, b, &r, q, Pptr);

    freePolynomial_spX (a);
    *a = r; 
}


void plainRemPolynomials_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, const Prime_ptr* Pptr)
{
    
    if (r == NULL) {
        return;
    }

    // a / 0 = 
    if (isZero_spX (b)) {
    	fprintf (stderr, "DUSP Error: Division by zero!\n");
    	exit (1);
    }
    
    // 0 / b = 0
    if (isZero_spX (a)) {
    	*r = NULL;
    	return;
    }

    // a / const = 
    if (isConstant_spX (b)) {
		*r = NULL;
		return;
    }
    
    // a/a = 1
    if (isEqual_spX (a, b)) {
		*r = NULL;
		return;
    }

    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
    	*r = deepCopyPolynomial_spX (a);
    	return;
    }
    
    polysize_t diff_deg = deg_a - deg_b;
    
    duspoly_t* rr = deepCopyPolynomial_spX (a); // rr = a;
    
    polysize_t i, j;
    elem_t lc_b_in = leadingCoeffInForm_spX (b);
    lc_b_in = smallprimefield_convert_in (lc_b_in, Pptr);
    lc_b_in = smallprimefield_inv (lc_b_in, Pptr);
    
    elem_t lc_rr_in = 0;
    elem_t b_in = 0;    
    
    // a = b*qq + rr
    for (i = diff_deg; i >= 0; i--) {	
		lc_rr_in = idxCoeffInForm_spX (rr, deg_b + i);
		lc_rr_in = smallprimefield_convert_in (lc_rr_in, Pptr);
	
		if (lc_rr_in) {
			lc_rr_in = smallprimefield_mul (lc_rr_in, lc_b_in, Pptr);

			rr->elems[deg_b+i] = 0;
			rr->lt -= 1;
			rr->alloc -= 1;
			for (j = deg_b-1;j >= 0; j--) {
				b_in = smallprimefield_convert_in (b->elems[j], Pptr);
				b_in = smallprimefield_mul (b_in, lc_rr_in, Pptr);
				b_in = smallprimefield_convert_out (b_in, Pptr);
		
				rr->elems[i+j] = smallprimefield_sub (rr->elems[i+j], b_in, Pptr);
			}
		}	
    }
    
    normalize_spX (&rr);
    
    if (rr == NULL || POLY_ALLOC(rr) == 0) {
		freePolynomial_spX (&rr);
		*r = NULL;
		return;
    }
    
    *r = rr;
}

void plainRemPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, const Prime_ptr* Pptr)
{
        
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pRem\n");
#endif

    if (r == NULL){
        return;
    }

    // a / 0 = 
    if (isZero_spX (b)) {
    	fprintf (stderr, "DUSP Error: Division by zero!\n");
    	exit (1);
    }
    
    // 0 / b = 0
    if (isZero_spX (a)) {
    	*r = NULL;
    	return;
    }

    // a / const = 
    if (isConstant_spX (b)) {
		*r = NULL;
		return;
    }
    
    // a/a = 1
    if (isEqual_spX (a, b)) {
		*r = NULL;
		return;
    }

    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
    	*r = deepCopyPolynomial_spX (a);
    	return;
    }
    
    polysize_t diff_deg = deg_a - deg_b;
    
    duspoly_t* rr = deepCopyPolynomial_spX (a); // rr = a;
    
    polysize_t i, j;
    elem_t lc_b_in = leadingCoeffInForm_spX (b);
    /* lc_b_in = smallprimefield_convert_in (lc_b_in, Pptr); */
    lc_b_in = smallprimefield_inv (lc_b_in, Pptr);
    
    elem_t lc_rr_in = 0;
    elem_t b_in = 0;    
    
    // a = b*qq + rr
    for (i = diff_deg; i >= 0; i--) {	
		lc_rr_in = idxCoeffInForm_spX (rr, deg_b + i);
		/* lc_rr_in = smallprimefield_convert_in (lc_rr_in, Pptr); */
	
		if (lc_rr_in) {
			lc_rr_in = smallprimefield_mul (lc_rr_in, lc_b_in, Pptr);

			rr->elems[deg_b+i] = 0;
			rr->lt -= 1;
			rr->alloc -= 1;
			for (j = deg_b-1;j >= 0; j--) {
				/* b_in = smallprimefield_convert_in (b->elems[j], Pptr); */
				b_in = smallprimefield_mul (b->elems[j], lc_rr_in, Pptr);
				/* b_in = smallprimefield_convert_out (b_in, Pptr); */
				rr->elems[i+j] = smallprimefield_sub (rr->elems[i+j], b_in, Pptr);
			}
		}	
    }
    
    normalize_spX (&rr);
    
    if (rr == NULL || POLY_ALLOC(rr) == 0) {
		freePolynomial_spX (&rr);
		*r = NULL;
		return;
    }
    
    *r = rr;
}


void plainRemPolynomials_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
    if (a == NULL) {
        return;
    }
 
    duspoly_t* r = NULL;
    plainRemPolynomials_spX (*a, b, &r, Pptr);
    freePolynomial_spX (a);

    *a = r;
}

void plainRemPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{	
    if (a == NULL) {
        return;
    }

    duspoly_t* r = NULL;
    plainRemPolynomialsInForm_spX (*a, b, &r, Pptr);
    freePolynomial_spX (a);

    *a = r;
}

int divideByMonicLinearInForm_spX (const duspoly_t* a, elem_t b, duspoly_t** q, elem_t* rem, const Prime_ptr* Pptr) {
    if (q == NULL && rem == NULL) {
        return -1;
    }
    
    if (a == NULL) {
        if (q != NULL) {
            freePolynomial_spX(q);
            *q = NULL;
        }
        if (rem != NULL) {
            *rem = smallprimefield_convert_in(0, Pptr);
        }
        return smallprimefield_convert_in(1, Pptr);
    }

    if (a->lt == 0) {
        if (q != NULL) {
            *q = NULL;
        }
        if (rem != NULL) {
            *rem = a->elems[0];
        }
        return (smallprimefield_convert_out(a->elems[0], Pptr) == 0);
    }

    duspoly_t* quo = makePolynomial_spX(a->lt);
    quo->lt = a->lt - 1;
    quo->elems[quo->lt] = a->elems[a->lt];
    for (int i = quo->lt-1; i >= 0; --i) {
        quo->elems[i] = smallprimefield_add(a->elems[i+1], smallprimefield_mul(quo->elems[i+1], b, Pptr), Pptr);
    }    

    int ret;
    if (rem != NULL) {
        *rem = smallprimefield_add(a->elems[0], smallprimefield_mul(quo->elems[0], b, Pptr), Pptr);
        ret = (smallprimefield_convert_out(*rem, Pptr) == 0);
    } else {
        elem_t tmp = smallprimefield_add(a->elems[0], smallprimefield_mul(quo->elems[0], b, Pptr), Pptr);
        ret = (smallprimefield_convert_out(tmp, Pptr) == 0);
    }

    if (q != NULL) {
        *q = quo;
    } else {
        freePolynomial_spX(&quo);
    }

    return ret;

}

void plainGCD_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, const Prime_ptr* Pptr)
{
    
    if (g == NULL){
        return;
    }

    elem_t etmp = 0;

    if (isZero_spX (a)) {
		monicPolynomial_spX (b, g, &etmp, Pptr);
		return;
	
    } else if (isZero_spX (b)) {
		monicPolynomial_spX (a, g, &etmp, Pptr);
		return;
	
    } else {
	
		duspoly_t* r0 = NULL;
		duspoly_t* r1 = NULL;
		duspoly_t* t  = NULL;
	
		plainRemPolynomials_spX (a, b, &t, Pptr);
		if (isZero_spX (t)) {
			monicPolynomial_spX (b, g, &etmp, Pptr);
			return;
		} else {
			r0 = deepCopyPolynomial_spX (b);
			r1 = t;
		}
	
		do {	    
			plainRemPolynomials_spX_inp (&r0, r1, Pptr); // r0 = r0 mod r1
			t  = r0;
			r0 = r1; 
			r1 = t;	    
		} while (!isZero_spX(r1));
	
		freePolynomial_spX (&r1);

		if (isZero_spX (r0)) {
            freePolynomial_spX (&r0);
			*g = NULL;
			return;
		}
	
		monicPolynomial_spX (r0, g, &etmp, Pptr);
        freePolynomial_spX (&r0);
    }
}


void plainGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pGCD\n");
#endif

    if (g == NULL){
        return;
    }

    elem_t etmp = 0;

    if (isZero_spX (a)) {
		monicPolynomialInForm_spX (b, g, &etmp, Pptr);
		return;
	
    } else if (isZero_spX (b)) {
		monicPolynomialInForm_spX (a, g, &etmp, Pptr);
		return;
	
    } else if (isEqual_spX (a, b)) {
        monicPolynomialInForm_spX (a, g, &etmp, Pptr);
        return;
    } else {
	   
       const duspoly_t* A;
       const duspoly_t* B;

       if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
            A = b;
            B = a;
       } else {
            A = a;
            B = b;
       }

		duspoly_t* r0 = NULL;
		duspoly_t* r1 = NULL;
		duspoly_t* t  = NULL;
	
		plainRemPolynomialsInForm_spX (A, B, &t, Pptr);
		if (isZero_spX (t)) {
            freePolynomial_spX (&t);
			monicPolynomialInForm_spX (B, g, &etmp, Pptr);
			return;
		} else {
			r0 = deepCopyPolynomial_spX (B);
			r1 = t;
		}
	
		do {	    
			plainRemPolynomialsInForm_spX_inp (&r0, r1, Pptr); // r0 = r0 mod r1
			swap_spX (&r0, &r1);    
		} while (!isZero_spX(r1));
	
        freePolynomial_spX (&r1);

		if (isZero_spX (r0)) {
            freePolynomial_spX (&r0);
			*g = NULL;
			return;
		}
	
		monicPolynomialInForm_spX (r0, g, &etmp, Pptr);
        freePolynomial_spX (&r0);
    }
}

void _plainGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, elem_t** lcrems, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pGCD\n");
#endif

    if (g == NULL){
        return;
    }

    elem_t etmp = 0;

    if (isZero_spX (a)) {
		monicPolynomialInForm_spX (b, g, &etmp, Pptr);
		return;
	
    } else if (isZero_spX (b)) {
		monicPolynomialInForm_spX (a, g, &etmp, Pptr);
		return;
	
    } else if (isEqual_spX (a, b)) {
        monicPolynomialInForm_spX (a, g, &etmp, Pptr);
        return;
    } else {
	   
       const duspoly_t* A;
       const duspoly_t* B;

       polysize_t min_deg = MIN_spX (a->lt, b->lt);
       elem_t* lcRems;
       if (*lcrems == NULL) {
           lcRems = (elem_t*) calloc (min_deg, sizeof(elem_t));
       } else {
           lcRems = *lcrems;
       }

       if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
            A = b;
            B = a;
       } else {
            A = a;
            B = b;
       }

		duspoly_t* r0 = NULL;
		duspoly_t* r1 = NULL;
		duspoly_t* t  = NULL;
	
		plainRemPolynomialsInForm_spX (A, B, &t, Pptr);
    
        if (!isZero_spX(t)){
            lcRems[t->lt] = t->elems[t->lt];
        }
		if (isZero_spX (t)) {
            freePolynomial_spX (&t);
			monicPolynomialInForm_spX (B, g, &etmp, Pptr);
			*lcrems = lcRems;
            return;
		} else {
			r0 = deepCopyPolynomial_spX (B);
			r1 = t;
		}
	
		do {	    
			plainRemPolynomialsInForm_spX_inp (&r0, r1, Pptr); // r0 = r0 mod r1
            if (!isZero_spX(r0)){
            lcRems[r0->lt] = r0->elems[r0->lt];
            }
			swap_spX (&r0, &r1);    
		} while (!isZero_spX(r1));
	
        freePolynomial_spX (&r1);

		if (isZero_spX (r0)) {
            freePolynomial_spX (&r0);
			*g = NULL;
            *lcrems = lcRems;
			return;
		}
	
		monicPolynomialInForm_spX (r0, g, &etmp, Pptr);
        freePolynomial_spX (&r0);
        *lcrems = lcRems;
    }
}

void _plainGCDInFormAtMaxK_spX (const duspoly_t* a, const duspoly_t* b, polysize_t k, duspoly_t** rk1, duspoly_t** rk, elem_t** lcrems, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pGCD\n");
#endif

    if (rk1 == NULL || rk == NULL){
        return;
    }

    elem_t etmp = 0;

    if (isZero_spX (a)) {
        *rk1 = deepCopyPolynomial_spX (b);
        *rk  = NULL;
		return;
    } else if (isZero_spX (b)) {
        *rk1 = deepCopyPolynomial_spX (a);
        *rk  = NULL;
		return;
    } else if (k < 0 || k > MIN_spX (a->lt, b->lt)) {
        return;
    } else if (isEqual_spX (a, b)) {
        if (k == MIN_spX(a->lt, b->lt)) {
            *rk1 = deepCopyPolynomial_spX (a);
            *rk  = deepCopyPolynomial_spX (b);
        } else {
            *rk1 = deepCopyPolynomial_spX (a);
            *rk = NULL;
        }
        return;
    } else {
       const duspoly_t* A;
       const duspoly_t* B;
       polysize_t min_deg = MIN_spX (a->lt, b->lt);
       elem_t* lcRems;
       if (*lcrems == NULL) {
           lcRems = (elem_t*) calloc (min_deg, sizeof(elem_t));
       } else {
           lcRems = *lcrems;
       }

       if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
            A = b;
            B = a;
       } else {
            A = a;
            B = b;
       }

		duspoly_t* r0 = NULL;
		duspoly_t* r1 = NULL;
		duspoly_t* t  = NULL;
	
		plainRemPolynomialsInForm_spX (A, B, &t, Pptr);
        
        // fprintf (stderr, "[k=%ld] t = ", k);  // TEST
        // printPolynomialOutForm_spX (t, Pptr); // TEST
 
        if (!isZero_spX(t)){
            lcRems[t->lt] = t->elems[t->lt];
            if (t->lt <= k) {
                *rk1 = deepCopyPolynomial_spX (B);
                *rk  = t;
                *lcrems = lcRems;
                return;
            }
        }
		if (isZero_spX (t)) {
            freePolynomial_spX (&t);
            *rk1 = deepCopyPolynomial_spX (B);
            *rk = NULL;
			*lcrems = lcRems;
            return;
		} else {
			r0 = deepCopyPolynomial_spX (B);
			r1 = t;
		}
	
		do {	    
			plainRemPolynomialsInForm_spX_inp (&r0, r1, Pptr); // r0 = r0 mod r1
            if (!isZero_spX(r0)){
            lcRems[r0->lt] = r0->elems[r0->lt];
            }
            if (isZero_spX (r0) || r0->lt <= k) {
                *rk1 = r1;
                *rk = r0;
                *lcrems = lcRems;
                return;
            }
			swap_spX (&r0, &r1);    
		} while (!isZero_spX(r1));
	
        freePolynomial_spX (&r1);

		if (isZero_spX (r0)) {
            freePolynomial_spX (&r0);
			*rk1 = NULL;
            *rk = NULL;
            *lcrems = lcRems;
			return;
		}

        *rk1 = r0;
        *rk = r1;
        *lcrems = lcRems;
    }
}


// a*u + b*v = gcd (a,b)
void plainExtGCD_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, const Prime_ptr* Pptr)
{
    
    elem_t etmp = 0;

    // 0 + b * (1/lc(g)) = monic(g)
    if (isZero_spX (a)) {
	   
       if (u != NULL){
		   *u = NULL;
       }

       if (v != NULL){
    		monicPolynomial_spX (b, g, &etmp, Pptr);
    		etmp = smallprimefield_convert_in (etmp, Pptr);
    		etmp = smallprimefield_inv (etmp, Pptr);
    		etmp = smallprimefield_convert_out (etmp, Pptr);
    		*v = constPolynomialInForm_spX (etmp, Pptr);
       }

		return;
	
    } else if (isZero_spX (b)) {
	   if (u != NULL){
    		monicPolynomial_spX (a, g, &etmp, Pptr);
    		etmp = smallprimefield_convert_in (etmp, Pptr);
    		etmp = smallprimefield_inv (etmp, Pptr);
    		etmp = smallprimefield_convert_out (etmp, Pptr);
    		*u = constPolynomialInForm_spX (etmp, Pptr);
       }
       if (v != NULL){
		  *v = NULL;
       }
		return;

    } 

    duspoly_t* q = NULL;
    duspoly_t* tmp = NULL;
    
    duspoly_t* r0 = deepCopyPolynomial_spX (a);
    duspoly_t* r1 = deepCopyPolynomial_spX (b);
    
    duspoly_t* s0 = constPolynomialInForm_spX (1, Pptr);
    duspoly_t* s1 = NULL; //zeroPolynomial_spX ();
    duspoly_t* s2 = NULL;
    
    duspoly_t* t0 = NULL;  // zeroPolynomial_spX ();
    duspoly_t* t1 = constPolynomialInForm_spX (1, Pptr);
    duspoly_t* t2 = NULL;
    
    do {
	
		plainDivPolynomials_spX_inp (&r0, r1, &q, Pptr);
		swap_spX (&r0, &r1);

		plainMulPolynomials_spX (s1, q, &tmp, Pptr);
		subPolynomials_spX (s0, tmp, &s2, Pptr);
		freePolynomial_spX (&s0);
		freePolynomial_spX (&tmp);
		s0 = s1;
		s1 = s2;

		plainMulPolynomials_spX (t1, q, &tmp, Pptr);
		subPolynomials_spX (t0, tmp, &t2, Pptr);
		freePolynomial_spX (&t0);
		freePolynomial_spX (&tmp);
		t0 = t1;
		t1 = t2;

        freePolynomial_spX (&q);
	
    } while (!isZero_spX (r1));

    freePolynomial_spX (&r1);
    freePolynomial_spX (&s2);
    freePolynomial_spX (&t2);
    
    if (isZero_spX (r0) || leadingCoeffInForm_spX (r0) == 1) {
     if (u != NULL){
		*u = s0;
     }
     if (v != NULL){
        *v = t0;
     }	   
     if (g != NULL){
        *g = r0;
     }
		return;
    }
    
    // monic gcd:
    duspoly_t* g_t = NULL;
    duspoly_t* u_t = NULL;
    duspoly_t* v_t = NULL;
    
    monicPolynomial_spX (r0, &g_t, &etmp, Pptr);
    freePolynomial_spX (&r0);
    
    if (g != NULL){
        *g = g_t;
    } else {
        freePolynomial_spX (&g_t);
    }
    
    etmp = smallprimefield_convert_in (etmp, Pptr);
    etmp = smallprimefield_inv (etmp, Pptr);
    etmp = smallprimefield_convert_out (etmp, Pptr);	
    
    if (etmp != 1) {
	   if (u != NULL){
        	scalarMulPolynomial_spX (s0, etmp , &u_t, Pptr);
        	freePolynomial_spX (&s0);
        	*u = u_t;        
       } else {
            freePolynomial_spX (&s0);
       }
       if (v != NULL){
    		scalarMulPolynomial_spX (t0, etmp , &v_t, Pptr);
    		freePolynomial_spX (&t0);
    		*v = v_t;
       } else {
            freePolynomial_spX (&t0);
       }
    } else {
		if (u != NULL){
            *u = s0;
        } else {
            freePolynomial_spX (&s0);
        }

        if (v != NULL){
		  *v = t0;
        } else {
            freePolynomial_spX (&t0);
        }
    }
}

// a*u + b*v = gcd (a,b)
void plainExtGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "pExGCD\n");
#endif
    elem_t etmp = 0;

    // 0 + b * (1/lc(g)) = monic(g)
    if (isZero_spX (a)) {
	
        if (u != NULL){
		  *u = NULL;
        }
        if (v != NULL){
    		monicPolynomialInForm_spX (b, g, &etmp, Pptr);
    		etmp = smallprimefield_inv (etmp, Pptr);
    		*v = constPolynomialInForm_spX (etmp, Pptr);            
        }
		return;
	
    } else if (isZero_spX (b)) {
	   if (u != NULL){
    		monicPolynomialInForm_spX (a, g, &etmp, Pptr);
    		etmp = smallprimefield_inv (etmp, Pptr);
    		*u = constPolynomialInForm_spX (etmp, Pptr);
       }
       if (v != NULL){
		  *v = NULL;
       }
		return;

    } 

    duspoly_t* q = NULL;
    duspoly_t* tmp = NULL;

    duspoly_t* r0;
    duspoly_t* r1;
    int isSwap = 0;
    if (degPolynomial_spX (a) > degPolynomial_spX (b)) {
        r0 = deepCopyPolynomial_spX (a);
        r1 = deepCopyPolynomial_spX (b);
    } else {
        r0 = deepCopyPolynomial_spX (b);
        r1 = deepCopyPolynomial_spX (a);
        isSwap = 1;
    }
    
    duspoly_t* s0 = constPolynomial_spX (1, Pptr);
    duspoly_t* s1 = NULL; //zeroPolynomial_spX ();
    duspoly_t* s2 = NULL;
    
    duspoly_t* t0 = NULL;  // zeroPolynomial_spX ();
    duspoly_t* t1 = constPolynomial_spX (1, Pptr);
    duspoly_t* t2 = NULL;
    
    do {
	    plainDivPolynomialsInForm_spX_inp (&r0, r1, &q, Pptr);
        // fprintf (stderr,"In pExtGCD, r_i := "); // TEST 
        // printPolynomialOutForm_spX (r0, Pptr); // TEST
        swap_spX (&r0, &r1);
		// plainMulPolynomialsInForm_spX (s1, q, &tmp, Pptr);
		mulPolynomialsInForm_spX (s1, q, &tmp, Pptr);
		subPolynomialsInForm_spX (s0, tmp, &s2, Pptr);
		freePolynomial_spX (&s0);
		freePolynomial_spX (&tmp);
		s0 = s1;
		s1 = s2;
		// plainMulPolynomialsInForm_spX (t1, q, &tmp, Pptr);
		mulPolynomialsInForm_spX (t1, q, &tmp, Pptr);
		subPolynomialsInForm_spX (t0, tmp, &t2, Pptr);
		freePolynomial_spX (&t0);
		freePolynomial_spX (&tmp);
		t0 = t1;
		t1 = t2;
        freePolynomial_spX (&q);
    } while (!isZero_spX (r1));

    freePolynomial_spX (&r1); 
    
    freePolynomial_spX (&s2);
    freePolynomial_spX (&t2);

    if (isZero_spX (r0) || leadingCoeffInForm_spX (r0) == smallprimefield_convert_in(1, Pptr)) {

        if (isSwap) {
            swap_spX (&s0, &t0);
        }

		if (u != NULL) {
          *u = s0;
        } else {
            freePolynomial_spX (&s0);
        }

        if (v != NULL) {
		  *v = t0;
        } else {
            freePolynomial_spX (&t0);
        }

		
        if (g != NULL) {
          *g = r0;
        } else {
            freePolynomial_spX (&r0);
        }

		return;
    }
    
    // monic gcd:
    duspoly_t* g_t = NULL;
    duspoly_t* u_t = NULL;
    duspoly_t* v_t = NULL;
    
    monicPolynomialInForm_spX (r0, &g_t, &etmp, Pptr);
    freePolynomial_spX (&r0);

    if (g != NULL){
        *g = g_t;
    } else {
        freePolynomial_spX (&g_t);
    }
	
    if (isSwap) {
        swap_spX (&s0, &t0);
    }

    if (etmp != smallprimefield_convert_in(1, Pptr)) {
        if (u != NULL){
    		etmp = smallprimefield_inv (etmp, Pptr); // TODO : goes into if
    		scalarMulPolynomialInForm_spX (s0, etmp , &u_t, Pptr);
    		freePolynomial_spX (&s0);
    		*u = u_t;            
        } else {
            freePolynomial_spX (&s0);
        }
	
        if (v != NULL){
    		scalarMulPolynomialInForm_spX (t0, etmp , &v_t, Pptr);
    		freePolynomial_spX (&t0);
    		*v = v_t;
        } else {
            freePolynomial_spX (&t0);
        }
    } else {
        if (u != NULL){
		  *u = s0;
        } else {
            freePolynomial_spX (&s0);
        }

        if (v != NULL){
		  *v = t0;
        } else {
            freePolynomial_spX (&t0);
        }
    }
}


void plainMulModPoly_spX (duspoly_t* a, duspoly_t* b, duspoly_t* mod, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL){
        return;
    }

    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }

    duspoly_t* ab = NULL;
    
    plainMulPolynomials_spX (a, b, &ab, Pptr);
    plainRemPolynomials_spX_inp (&ab, mod, Pptr); 

    *c = ab;
}

void plainMulModPolyInForm_spX (duspoly_t* a, duspoly_t* b, duspoly_t* mod, duspoly_t** c, const Prime_ptr* Pptr)
{
    if (c == NULL){
        return;
    }
    
    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }

    duspoly_t* ab = NULL;
    
    plainMulPolynomialsInForm_spX (a, b, &ab, Pptr);
    plainRemPolynomialsInForm_spX_inp (&ab, mod, Pptr); 

    *c = ab;
}


void plainSqrModPoly_spX (duspoly_t* a, duspoly_t* mod, duspoly_t** a2, const Prime_ptr* Pptr)
{
    if (a2 == NULL){
        return;
    }
    
    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }
    
    duspoly_t* aa = NULL;

    plainSqrPolynomial_spX (a, &aa, Pptr);
    plainRemPolynomials_spX_inp (&aa, mod, Pptr);    
    
    *a2 = aa;
}

void plainSqrModPolyInForm_spX (duspoly_t* a, duspoly_t* mod, duspoly_t** a2, const Prime_ptr* Pptr)
{
    if (a2 == NULL){
        return;
    }

    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }
    
    duspoly_t* aa = NULL;
    
    plainSqrPolynomialInForm_spX (a, &aa, Pptr);
    plainRemPolynomialsInForm_spX_inp (&aa, mod, Pptr);    
    
    *a2 = aa;
}

void plainInvModPoly_spX (duspoly_t* a, duspoly_t* mod, duspoly_t** am, const Prime_ptr* Pptr)
{
    if (am == NULL){
        return;
    }

    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }
    
    duspoly_t* g = NULL;
    duspoly_t* r = NULL;
    duspoly_t* s = NULL;

    plainExtGCDInForm_spX (a, mod, &r, &s, &g, Pptr);
    
    freePolynomial_spX (&s);

    if (!isOne_spX (g)) {
		fprintf (stderr, "DUSP Error: Input polynomials are not invertible!\n");
		fprintf (stderr, "gcd of input polynomials is ");
		printPolynomial_spX (g, Pptr);
		exit (1);
    }

    freePolynomial_spX (&g);
    *am = r;
}


void plainInvModPolyInForm_spX (duspoly_t* a, duspoly_t* mod, duspoly_t** am, const Prime_ptr* Pptr)
{
    if (am == NULL){
        return;
    }

    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }
    
    duspoly_t* g = NULL;
    duspoly_t* r = NULL;
    duspoly_t* s = NULL;

    plainExtGCDInForm_spX (a, mod, &r, &s, &g, Pptr);
    
    freePolynomial_spX (&s);

    if (!isOneInForm_spX (g, Pptr)) {
		fprintf (stderr, "DUSP Error: Input polynomials are not invertible!\n");
		fprintf (stderr, "gcd of input polynomials is ");
		printPolynomial_spX (g, Pptr);
		exit (1);
    }

    freePolynomial_spX (&g);
    *am = r;
}


int isInvertibleModPoly_spX (duspoly_t* a, duspoly_t* mod, const Prime_ptr* Pptr)
{
    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }
    
    duspoly_t* g = NULL;

    plainGCD_spX (a, mod, &g, Pptr);
    
    int isOne = isOne_spX (g);
    freePolynomial_spX (&g);
    
    return isOne;
}

int isInvertibleModPolyInForm_spX (duspoly_t* a, duspoly_t* mod, const Prime_ptr* Pptr)
{
    if (isZero_spX (mod)) {
		fprintf (stderr, "DUSP Error: Polynomial Modulus cannot be zero!\n");
		exit (1);
    }
    
    duspoly_t* g = NULL;

    plainGCDInForm_spX (a, mod, &g, Pptr);
    
    int isOne = isOneInForm_spX (g, Pptr);
    freePolynomial_spX (&g);
    
    return isOne;
}

void derivativePolyInForm_spX (duspoly_t* a, duspoly_t** ap, const Prime_ptr* Pptr) 
{
    if (ap == NULL){
        return;
    }

    if (isZero_spX (a) || POLY_LT(a) == 0) {
        freePolynomial_spX (&*ap);
        return;
    }

    duspoly_t* diff = makePolynomial_spX (POLY_ALLOC (a) -1);
    POLY_LT (diff) = POLY_LT (a) - 1; // bound

    polysize_t i;
    for (i = 0; i <= POLY_LT(diff); i++){
        diff->elems[i] = smallprimefield_mul (a->elems[i+1], 
            smallprimefield_convert_in (i+1, Pptr), Pptr);
    }

    normalize_spX (&diff);
    *ap = diff;
}

void derivativePolyInForm_spX_inp (duspoly_t** a, const Prime_ptr* Pptr) 
{    
    if (isZero_spX (*a) || POLY_LT(*a) == 0) {
        freePolynomial_spX (&*a);
        return;
    }

    duspoly_t* diff = *a;
    polysize_t deg_a = degPolynomial_spX (*a);
    
    polysize_t i;
    for (i = 0; i <= deg_a-1; i++){
        diff->elems[i] = smallprimefield_mul (diff->elems[i+1], 
            smallprimefield_convert_in (i+1, Pptr), Pptr);
    }
    diff->elems[deg_a] = 0;  

    normalize_spX (a);
    // *a = diff;
}


/*******************
 * Interpolation 
 ******************/

void evalDivPolynomials_spX (const duspoly_t* a, const elem_t u, duspoly_t** q, elem_t* au, const Prime_ptr* Pptr)
{
    if (isZero_spX (a) || q == NULL) {
        return;
    }
    if (isConstant_spX (a)) {
        *q = deepCopyPolynomial_spX (a);
        *au = a->elems[0];
        return;
    }

    polysize_t deg_a = degPolynomial_spX (a);
    elem_t u_inv = smallprimefield_inv (u, Pptr);
    elem_t aa = smallprimefield_sub (0, smallprimefield_mul (a->elems[0], u_inv, Pptr), Pptr);
    elem_t tmp = u;
    duspoly_t* qq = makePolynomial_spX (deg_a+1);
    qq->lt = deg_a;
    qq->elems[0] = aa;
    for (polysize_t i = 1; i <= deg_a; i++) {
        qq->elems[i] = smallprimefield_mul (smallprimefield_sub (qq->elems[i-1], a->elems[i], Pptr), u_inv, Pptr);
        aa = smallprimefield_add (aa, smallprimefield_mul (qq->elems[i], tmp, Pptr), Pptr);
        tmp = smallprimefield_mul (tmp, u, Pptr);
    }
    normalize_spX (&qq);
    *q = qq;
    *au = aa;
}

void lagrangeBasisInForm_spX (const elem_t* t, const polysize_t n, duspoly_t** a, const Prime_ptr* Pptr)
{
    if (a == NULL || n < 1) {
        return;
    }

    duspoly_t* ret = makePolynomial_spX (n+1);
    ret->elems[0] = smallprimefield_sub (0, t[0], Pptr); // neg(t[0], pr)
    if (n==1) {
        *a = ret; 
        return;
    }
    ret->elems[1] = smallprimefield_convert_in (1, Pptr); 
    for (polysize_t i = 2; i <= n; i++) {
        ret->elems[i] = 0; // TODO:
        for (polysize_t j = i; j >= 0; --j) {
            ret->elems[j] = smallprimefield_mul (ret->elems[j], smallprimefield_sub(0,t[i-1], Pptr), Pptr);
            if (j) {
                ret->elems[j] = smallprimefield_add (ret->elems[j], ret->elems[j-1], Pptr);
            }
        }
    }
    ret->lt = n;
    normalize_spX (&ret);
    *a = ret;
}   


/*****************************
 ** Structured Rand Polynomial
 ****************************/

duspoly_t* randomPolynomialInForm_spX (polysize_t n, const Prime_ptr* Pptr) 
{
    if (n == 0) {
        return NULL;
    }

    static int initTest = 0;
    while(!initTest) {
        srand (time(0));
        initTest = 1;
    }

    duspoly_t* poly = makePolynomial_spX (n);
    poly->lt = n-1;

    for (polysize_t i = 0; i < n; i++) {
        poly->elems[i] = smallprimefield_convert_in (rand(), Pptr);
    }

    normalize_spX (&poly);

    if (isZero_spX (poly)) {
        freePolynomial_spX (&poly);
        poly = randomPolynomialInForm_spX (n, Pptr);
    }

    return poly;
}


duspoly_t* binPolynomial_spX (polysize_t n, const Prime_ptr* Pptr)
{

    if (n == 0) {
        return NULL;
    }

    polysize_t sz = n+1;

    duspoly_t* poly = makePolynomial_spX (sz);
    POLY_LT (poly) = sz-1;

    poly->elems[0] = 1;
    poly->elems[sz-1] = 1;

    return poly;
}


duspoly_t* binPolynomialInForm_spX (polysize_t n, const Prime_ptr* Pptr)
{

    if (n == 0) {
        return NULL;
    }

    polysize_t sz = n+1;

    duspoly_t* poly = makePolynomial_spX (sz);
    POLY_LT (poly) = sz-1;

    elem_t one = smallprimefield_convert_in (1, Pptr);

    poly->elems[0] = one;
    poly->elems[sz-1] = one;

    return poly;
}


/***************************
 ** Karatsuba Based Functions
 ***************************/

void KaratsubaMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "KarMul\n");
#endif

    if (c == NULL){
        return;
    }

    if (isZero_spX (a) || isZero_spX (b)) {
        *c = NULL;
        return;
    }
    
    polysize_t da = degPolynomial_spX (a);
    polysize_t db = degPolynomial_spX (b);

    if (MAX_spX (da,db) < KARATSUBA_CROSSOVER) {
        plainMulPolynomialsInForm_spX (a, b, c, Pptr);
        return;
    }

    const duspoly_t* A;
    const duspoly_t* B;

    if (da < db) {
        A = b;
        B = a;
        da = degPolynomial_spX (A);
        db = degPolynomial_spX (B);
    } else {
        A = a;
        B = b;
    }

    polysize_t k = (da + 1) >> 1;
    polysize_t nn = k << 1;
    // a = a0 + a1*X^k :
    duspoly_t* a0 = NULL;
    duspoly_t* a1 = NULL;
    // b= b0 + b1*X^k :
    duspoly_t* b0 = NULL;
    duspoly_t* b1 = NULL;
    // TODO: could be in-place:
    leftSplitPolynomial_spX  (A, &a0, k); 
    rightShiftPolynomial_spX (A, &a1, k);
    leftSplitPolynomial_spX  (B, &b0, k);
    rightShiftPolynomial_spX (B, &b1, k);

    duspoly_t* h0 = NULL;
    duspoly_t* h1 = NULL;
    duspoly_t* h2 = NULL;
    duspoly_t* h11 = NULL;
    duspoly_t* h21 = NULL;
    duspoly_t* h22 = NULL;
    duspoly_t* a0a1 = NULL;
    duspoly_t* b0b1 = NULL;
    duspoly_t* h0h1 = NULL;
    duspoly_t* h0h1h2 = NULL;

    if (isZero_spX (b1)) {
        KaratsubaMulPolynomialsInForm_spX (a0, b0, &h0, Pptr); // h0 = a0*b0
        addPolynomialsInForm_spX (a0, a1, &a0a1, Pptr); // a0a1 = a0 + a1;
        KaratsubaMulPolynomialsInForm_spX (a0a1, B, &h2, Pptr); // h2 = (a0 + a1)*(b)
        subPolynomialsInForm_spX_inp (&h2, h0, Pptr);
        leftShiftPolynomial_spX (h2, &h21, k); // h21 = h2<<k
        addPolynomialsInForm_spX_inp (&h21, h0, Pptr);
        freePolynomial_spX (&a0a1);
        freePolynomial_spX (&h0);
        freePolynomial_spX (&h2);
        freePolynomial_spX (&a0);
        freePolynomial_spX (&b0);
        freePolynomial_spX (&a1);
        freePolynomial_spX (&b1);
        *c = h21;
        return;
    }

    // normal case: deg(a)/2 < deg(b) :
    KaratsubaMulPolynomialsInForm_spX (a0, b0, &h0, Pptr); // h0 = a0*b0
    KaratsubaMulPolynomialsInForm_spX (a1, b1, &h1, Pptr); // h1 = a1*b1
    addPolynomialsInForm_spX (a0, a1, &a0a1, Pptr); // a0a1 = a0 + a1;
    addPolynomialsInForm_spX (b0, b1, &b0b1, Pptr); // b0b1 = b0 + b1;
    freePolynomial_spX (&a0);
    freePolynomial_spX (&b0);
    freePolynomial_spX (&a1);
    freePolynomial_spX (&b1);
    KaratsubaMulPolynomialsInForm_spX (a0a1, b0b1, &h2, Pptr); // h2 = (a0 + a1)*(b0 + b1)
    freePolynomial_spX (&a0a1); 
    freePolynomial_spX (&b0b1);
    subPolynomialsInForm_spX (h2, h0, &h21, Pptr); // h21 = h2 - h0
    subPolynomialsInForm_spX (h21, h1, &h22, Pptr); // h22 = h21 - h1 = h2 - h0 - h1
    freePolynomial_spX (&h2);
    freePolynomial_spX (&h21);
    leftShiftPolynomial_spX (h1, &h11, nn); // h11 = h1<<n
    freePolynomial_spX (&h1); 
    leftShiftPolynomial_spX (h22, &h21, k); // h21 = h22<<k
    freePolynomial_spX (&h22);
    addPolynomialsInForm_spX (h0, h11, &h0h1, Pptr); // h0h1 = h0 + h11
    addPolynomialsInForm_spX (h0h1, h21, &h0h1h2, Pptr); // h0h1h2 = h0h1 + h21 = h0 + h11 + h21 = h0 + h1<<n + h22<<k
    freePolynomial_spX (&h0); 
    freePolynomial_spX (&h11);
    freePolynomial_spX (&h21);
    freePolynomial_spX (&h0h1);
    *c = h0h1h2;    
}

void KaratsubaMulPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr) 
{
	duspoly_t* cc = NULL;
	KaratsubaMulPolynomialsInForm_spX (*a, b, &cc, Pptr);
	freePolynomial_spX (a);
	*a = cc;
}

void KaratsubaSqrPolynomialsInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr) 
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "KarSqr\n");
#endif

    if (a2 == NULL){
        return;
    }
	duspoly_t* cc = NULL;
	KaratsubaMulPolynomialsInForm_spX (a, a, &cc, Pptr);
	*a2 = cc;
}

/////////////////////
// FFT-based Mul
/////////////////////
void fftMulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr) 
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "6-th_Mul-FFT\n");
#endif

    if (degPolynomial_spX (a) < PLAINMUL_CROSSOVER && 
        degPolynomial_spX (b) < PLAINMUL_CROSSOVER) {
        KaratsubaMulPolynomialsInForm_spX (a, b, c, Pptr);
        return;
    } else if (Pptr->prime < 1337006139375617UL || Pptr->prime > 9214646312576745473UL) {
        // they are min & max of fourier_primes_u64_table[10000]
        KaratsubaMulPolynomialsInForm_spX (a, b, c, Pptr);
        return;
    }
    
   // poly * 0 = 0
    if (isZero_spX (a) || isZero_spX (b)){
	   *c = zeroPolynomial_spX ();
	   return;
    }

    // poly * elem = 
    if (isConstant_spX (a)) {
    	scalarMulPolynomialInForm_spX (b, a->elems[0], c, Pptr);
        return; 
    } else if (isConstant_spX (b)) {
    	scalarMulPolynomialInForm_spX (a, b->elems[0], c, Pptr);
    	return; 
    }

    polysize_t max_deg = MAX_spX (degPolynomial_spX (a), degPolynomial_spX (b));

    // TODO: k is not efficient 
    int K = 4; // try with 8, ... 
    int e = 1;
    polysize_t l = max_deg+1; // technically, it's max alloc size
    while (l > 0) {
        l = l/K;
        e++;
    } 
    polysize_t N = pow (K, e); // max_deg + 1 <= N
    polysize_t N2 = N>>1; // polynomials have alloc_size at most N2
    // fprintf (stderr, "DUSP FFT Params: max_deg = %ld\n", max_deg); // TEST
    // fprintf (stderr, "DUSP FFT Params: e = %d\n", e); // TEST
    // fprintf (stderr, "DUSP FFT Params: K^e = %ld\n", N); // TEST
    // fprintf (stderr, "DUSP FFT Params: K^e/2 = %ld\n", N2); // TEST
    if (max_deg >= N2) {
        fprintf (stderr, "max_deg = %ld > K^e/2 = %d^%d", max_deg, K, e);
        exit(1);
    }
    usfixn64 p_u64 = (unsigned long long) Pptr->prime;
    // init montgomery triple:
    montgomery_triple P;
    init_montgomery_triple (&P, p_u64);
    usfixn64* x_data = (usfixn64*) calloc (N, sizeof(usfixn64));
    usfixn64* y_data = (usfixn64*) calloc (N, sizeof(usfixn64));
    for (long i = 0; i <= degPolynomial_spX (a); ++i) { x_data[i] = (usfixn64) a->elems[i]; }
    for (long i = 0; i <= degPolynomial_spX (b); ++i) { y_data[i] = (usfixn64) b->elems[i]; }
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
    for (int i = 0; i < N; i++) { x_data[i] = mult_ab_mod_p (x_data[i], n_inv, &P); }
    // make output
    setCoefsInForm_spX (c, (elem_t*) x_data, N);
    normalize_spX (c);
    free (x_data);
    free (y_data);
}

void mulPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr) 
{
    if (a == NULL ||  b == NULL) {
        plainMulPolynomialsInForm_spX (a, b, c, Pptr);
        return;
    }
    fftMulPolynomialsInForm_spX (a, b, c, Pptr);
}

void mulPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr) 
{
    if (*a == NULL || b == NULL || Pptr->prime < 1337006139375617UL || Pptr->prime > 9214646312576745473UL) {
        // they are min & max of fourier_primes_u64_table[10000]
    	KaratsubaMulPolynomialsInForm_spX_inp (a, b, Pptr);	
        return;
    } else {
        duspoly_t* c = NULL;
        fftMulPolynomialsInForm_spX (*a, b, &c, Pptr);
        freePolynomial_spX (a);
        *a = c;
    }
}


duspoly_t* mulManyPolynomialsInForm_spX(duspoly_t const*const* polys, int n, const Prime_ptr* Pptr) {
    if (n == 0 || polys == NULL) {
        return NULL;
    }
    if (n < 2) {
        return deepCopyPolynomial_spX(polys[0]);
    }
    duspoly_t* prod;
    mulPolynomialsInForm_spX(polys[0], polys[1], &prod, Pptr);
    duspoly_t* tmpProd;
    for (int i = 2; i < n; ++i) {
        mulPolynomialsInForm_spX(prod, polys[i], &tmpProd, Pptr);
        freePolynomial_spX(&prod);
        prod = tmpProd;
    }
    return prod;
}


void sqrPolynomialInForm_spX (const duspoly_t* a, duspoly_t** a2, const Prime_ptr* Pptr)
{
    if (a2 == NULL) {
        return;
    }
	if (degPolynomial_spX (a) < 1000) {
		plainSqrPolynomialInForm_spX (a, a2, Pptr);
		return;
	}
    fftMulPolynomialsInForm_spX (a, a, a2, Pptr);
}

void exponentiatePolynomialInForm_spX (const duspoly_t* a, polysize_t n, duspoly_t** an, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "ExpPoly\n");
#endif
    if (an == NULL){
        return;
    }
    if (n == 0) {
        elem_t one = smallprimefield_convert_in (1, Pptr);
        *an = constPolynomialInForm_spX (one, Pptr);
        return;
    } else if (n == 1) {
        *an = deepCopyPolynomial_spX (a);
        return;
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
                mulPolynomialsInForm_spX_inp (&mul, b, Pptr);
            }
        }
        
        sqrPolynomialInForm_spX (b, &tmp, Pptr);
        
        freePolynomial_spX (&b);
        b = tmp; tmp = NULL;

        nn >>= 1;
    }

    if (isZero_spX (mul)) {
        mul = deepCopyPolynomial_spX (b);
    } else {
        mulPolynomialsInForm_spX_inp (&mul, b, Pptr);
    }


    freePolynomial_spX (&b);
    *an = mul;
}

void modExponentiatePolynomialInForm_spX (const duspoly_t* a, polysize_t n, const duspoly_t* f, duspoly_t** an, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "ModExpPoly\n");
#endif

    if (an == NULL){
        return;
    }

    if (n == 0) {
        elem_t one = smallprimefield_convert_in (1, Pptr);
        *an = constPolynomialInForm_spX (one, Pptr);
        return;
    } else if (n == 1) {
        plainRemPolynomialsInForm_spX (a, f, an, Pptr);
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
                mulPolynomialsInForm_spX_inp (&mul, b, Pptr);
            }
        }
        
        sqrPolynomialInForm_spX (b, &tmp, Pptr);
        
        if (POLY_LT(tmp) > 2*deg_f) {
            plainRemPolynomialsInForm_spX_inp (&tmp, f, Pptr);
        }

        freePolynomial_spX (&b);

        b = tmp; tmp = NULL;

        nn >>= 1;
    }

    if (isZero_spX (mul)) {
        mul = deepCopyPolynomial_spX (b);
    } else {
        mulPolynomialsInForm_spX_inp (&mul, b, Pptr);
    }

    plainRemPolynomialsInForm_spX_inp (&mul, f, Pptr);

    freePolynomial_spX (&b);
    *an = mul;
}

/********************************
 ** Power Series Based Functions 
 *******************************/

void plainPSInversionInForm_spX (const duspoly_t* a, duspoly_t** ai, polysize_t l, const Prime_ptr* Pptr)
{
    if (ai == NULL) {
        return;
    }
    if (isZero_spX (a)) {
		return;
    }
    if (a->elems[0] == 0) {
		fprintf (stderr, "DUSP Error, input polynomial has a wrong format in power series inversion!\n");
		exit(1);
    }
   
    duspoly_t* aa = NULL;
    elem_t e0 = 0;
    int isOne = 0;
    e0 = smallprimefield_inv (a->elems[0], Pptr); // e^{-1}
    if (e0 == smallprimefield_convert_in (1, Pptr)) {
		isOne = 1;
		aa = deepCopyPolynomial_spX (a);
    } else {
		scalarMulPolynomialInForm_spX (a, e0, &aa, Pptr);
    }
    negPolynomialInForm_spX (aa, Pptr); // aa = -f 
        
    duspoly_t* g0  = constPolynomial_spX (1, Pptr);
    duspoly_t *g1=NULL, *g2=NULL, *gs=NULL, *fgs=NULL;    
    /* polysize_t r = cail (Log2_spX (degPolynomial_spX (aa) + 1)); */
    polysize_t pow = 1;

    for (polysize_t i = 1; i <= l; ++i) {
		addPolynomialsInForm_spX (g0, g0, &g2, Pptr);  // g2 = g0 + g0
		sqrPolynomialInForm_spX (g0, &gs, Pptr);  // gs = g0 * g0
		mulPolynomialsInForm_spX (aa, gs, &fgs, Pptr); // fgs = -f*g0*g0
		addPolynomialsInForm_spX (g2, fgs, &g1, Pptr);      // g1 = g2 + fgs
		freePolynomial_spX (&g0); 
		freePolynomial_spX (&g2); 
		freePolynomial_spX (&gs); 
		freePolynomial_spX (&fgs);
		pow <<= 1; // pow = 2^(i+1);
		leftSplitPolynomial_spX (g1, &g0, pow); // g0 = g1 mod x^2^(i+1)
		freePolynomial_spX (&g1); 
    }
    
    if (isZero_spX (g0)) {
		*ai = NULL;
		return;
    }
    if (isOne) {
		*ai = g0;
		return;
    }
    scalarMulPolynomialInForm_spX (g0, e0, &g2, Pptr);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&g0);
    *ai = g2;
    
}

void reversePolynomial_spX (const duspoly_t* a, duspoly_t** ra, polysize_t d, const Prime_ptr* Pptr) 
{
    if (ra == NULL){
        return;
    }    
    if (isZero_spX (a)) {
		*ra = NULL;
		return;
    }
    polysize_t deg_a = degPolynomial_spX (a);
    if (deg_a > d) {
		*ra = deepCopyPolynomial_spX (a);
		return;
    }
    duspoly_t* rra = makePolynomial_spX (d+1);
    POLY_LT (rra) = d;
    polysize_t i;
    // reverse data in rra:
    for (i = 0; i <= deg_a; ++i) {
		rra->elems[d - deg_a + i] = a->elems[deg_a-i];
    }
    normalize_spX (&rra);
    *ra = rra;
}


void fastDivPolynomialInForm_wPSInv_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, duspoly_t** q, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "fastDiv\n");
#endif

    // a / 0 = Error
    if (isZero_spX (b)) {
    	fprintf (stderr, "DUSP Error: Division by zero!\n");
    	exit (1);
    }
    
    // 0 / b = 0
    if (isZero_spX (a)) {
        if (r != NULL){
    	   *r = NULL; }
        if (q != NULL){
    	   *q = NULL; }
    	return;
    }

    // a/a = 1
    if (isEqual_spX (a, b)) {
		if (r != NULL){
            *r = NULL; }
        if (q != NULL){
     		*q = constPolynomial_spX (1, Pptr); }
		return;
    }
    
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);

    // deg(a) < deg(b) => a = 0*b + a
    if (deg_a < deg_b) {
        if (r != NULL){
    	   *r = deepCopyPolynomial_spX (a); }
        if (q != NULL){
    	   *q = NULL; }
    	return;
    }
    if (deg_a - deg_b < 11 && MAX_spX (deg_a, deg_b) < 1000) {
        plainDivPolynomialsInForm_spX (a, b, r, q, Pptr);
        return;
    }
    if (deg_a < PLAINDIV_CROSSOVER || deg_b < PLAINDIV_CROSSOVER) {
		plainDivPolynomialsInForm_spX (a, b, r, q, Pptr);
		return;
    }

    duspoly_t *ar = NULL, *br=NULL, *br_inv=NULL, *qr=NULL, *qt=NULL, *rr=NULL;
    polysize_t ldeg = LogCeiling_spX (deg_a - deg_b + 1); // ceil (Log2_spX (deg_a - deg_b + 1));
    reversePolynomial_spX (a, &ar, deg_a, Pptr);
    reversePolynomial_spX (b, &br, deg_b, Pptr);
    plainPSInversionInForm_spX (br, &br_inv, ldeg, Pptr); // 
    freePolynomial_spX (&br);
    reallocatePolynomial_spX (&ar, deg_a - deg_b + 1);
    reallocatePolynomial_spX (&br_inv, deg_a - deg_b + 1);
    mulPolynomialsInForm_spX (ar, br_inv, &qr, Pptr);
    // qr = ar * br_inv rem x^{deg_a - deg_b + 1};
    reallocatePolynomial_spX (&qr, deg_a - deg_b + 1);
    freePolynomial_spX (&ar);
    freePolynomial_spX (&br_inv);
    reversePolynomial_spX (qr, &qt, deg_a - deg_b, Pptr);
    freePolynomial_spX (&qr);
    mulPolynomialsInForm_spX (b, qt, &rr, Pptr);
    if (q != NULL) {
        *q = qt;
    } else {
        freePolynomial_spX (&qt);}
    if (r != NULL) {
        subPolynomialsInForm_spX (a, rr, r, Pptr);}
    freePolynomial_spX (&rr);
}

int isDividablePolys_spX (const duspoly_t* a, const duspoly_t* b, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
        if (isZero_spX (b)) {
            return 1;
        }
        return 0;
    }
    
    duspoly_t* r = NULL;
    duspoly_t* q = NULL;

    // plainRemPolynomialsInForm_spX (a, b, &r, Pptr);
    fastDivPolynomialInForm_wPSInv_spX (a, b, &r, &q, Pptr);

    freePolynomial_spX (&q);

    if (isZero_spX (r)) {
        return 1;
    }
    
    freePolynomial_spX (&r);
    return 0;
}

int isDividablePolysQuoInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
        if (isZero_spX (b)) {
            return 1;
        }
        return 0;
    }
    
    duspoly_t* r = NULL;
    
    // plainRemPolynomialsInForm_spX (a, b, &r, Pptr);
    // plainDivPolynomialsInForm_spX (a, b, &r, q, Pptr);
    fastDivPolynomialInForm_wPSInv_spX (a, b, &r, q, Pptr);

    if (isZero_spX (r)) {
        freePolynomial_spX (&r);
        return 1;
    }
    
    freePolynomial_spX (&r);
    
    return 0;
}

void divPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, duspoly_t** q, const Prime_ptr* Pptr)
{
    // fastDivPolynomialInForm_wPSInv_spX (a, b, r, q, Pptr);    
    plainDivPolynomialsInForm_spX (a, b, r, q, Pptr);
}

void divPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, duspoly_t** q, const Prime_ptr* Pptr)
{
plainDivPolynomialsInForm_spX_inp (a, b, q, Pptr);
return;
    if (degPolynomial_spX (*a) < PLAINDIV_CROSSOVER || 
        degPolynomial_spX (b) < PLAINDIV_CROSSOVER) {
        plainDivPolynomialsInForm_spX_inp (a, b, q, Pptr);
        return;
    } else if (Pptr->prime < 1337006139375617UL || Pptr->prime > 9214646312576745473UL) {
        // KaratsubaMulPolynomialsInForm_spX (a, b, c, Pptr);
        plainDivPolynomialsInForm_spX_inp (a, b, q, Pptr);
        return;
    } else {
        duspoly_t* r = NULL;
        fastDivPolynomialInForm_wPSInv_spX (*a, b, q, &r, Pptr);
        freePolynomial_spX (a);
        *a = r;
        return;
    }
}

void remPolynomialsInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** r, const Prime_ptr* Pptr)
{
plainRemPolynomialsInForm_spX (a, b, r, Pptr);
return;
    if (degPolynomial_spX (a) < PLAINDIV_CROSSOVER || 
        degPolynomial_spX (b) < PLAINDIV_CROSSOVER) {
        plainRemPolynomialsInForm_spX (a, b, r, Pptr);
        return;
    } else if (Pptr->prime < 1337006139375617UL || Pptr->prime > 9214646312576745473UL) {
        // KaratsubaMulPolynomialsInForm_spX (a, b, c, Pptr);
        plainRemPolynomialsInForm_spX (a, b, r, Pptr);
        return;
    } else {
        duspoly_t* q = NULL;
        fastDivPolynomialInForm_wPSInv_spX (a, b, &q, r, Pptr);
        freePolynomial_spX (&q);
        return;
    }    
}

void remPolynomialsInForm_spX_inp (duspoly_t** a, const duspoly_t* b, const Prime_ptr* Pptr)
{
plainRemPolynomialsInForm_spX_inp (a, b, Pptr);
return;
    if (degPolynomial_spX (*a) < PLAINDIV_CROSSOVER || 
        degPolynomial_spX (b) < PLAINDIV_CROSSOVER) {
        plainRemPolynomialsInForm_spX_inp (a, b, Pptr);
        return;
    } else if (Pptr->prime < 1337006139375617UL || Pptr->prime > 9214646312576745473UL) {
        // KaratsubaMulPolynomialsInForm_spX (a, b, c, Pptr);
        plainRemPolynomialsInForm_spX_inp (a, b, Pptr);
        return;
    } else {
        duspoly_t *q = NULL, *r = NULL;
        fastDivPolynomialInForm_wPSInv_spX (*a, b, &q, &r, Pptr);
        freePolynomial_spX (&q);
        freePolynomial_spX (a);
        *a = r;
        return;
    }    
}

/****************************
 ** Half GCD Based Functions
 ****************************/

// get polynomial based on matrix indices, i=0..1, j=0..1
duspoly_t* getPolyFromMat4InForm_spX (dusmat4_t* A, polysize_t i, polysize_t j)
{
    if (A == NULL) {
		return NULL;
    }
    
    if (i == 0 && j == 0) {
		return deepCopyPolynomial_spX (A->polys[0]);
    } else if (i == 0 && j == 1) {
		return deepCopyPolynomial_spX (A->polys[1]);
    } else if (i == 1 && j == 0) {
		return deepCopyPolynomial_spX (A->polys[2]);
    } else if (i == 1 && j == 1) {
		return deepCopyPolynomial_spX (A->polys[3]);
    }

    return NULL;
}

// get polynomial based on matrix indices, i=0..1, j=0..1
duspoly_t* getPolyFromMat4InForm_spX_inp (dusmat4_t* A, polysize_t i,polysize_t j)
{
    if (A == NULL) {
		return NULL;
    }
    
    if (i == 0 && j == 0) {
		return (A->polys[0]);
    } else if (i == 0 && j == 1) {
		return (A->polys[1]);
    } else if (i == 1 && j == 0) {
		return (A->polys[2]);
    } else if (i == 1 && j == 1) {
		return (A->polys[3]);
    }
    
    return NULL;
}

// set polynomial based on matrix indices, i=0..1, j=0..1
void setPolyToMat4InForm_spX (dusmat4_t** A, duspoly_t* a, polysize_t i, polysize_t j)
{
    if (i > 1 || j > 1) {
		fprintf (stderr, "DUSP Error, cannot set polynomial to the Matrix!\n");
		exit (1);
    }
    
    if (*A == NULL) {
		*A = (dusmat4_t*) malloc (sizeof (dusmat4_t));
		(*A)->polys[0] = NULL;
		(*A)->polys[1] = NULL;
		(*A)->polys[2] = NULL;
		(*A)->polys[3] = NULL;
    }

    // for non normalized input:
    /* normalize_spX (&a); */
    
    if (i == 0 && j == 0) {
		freePolynomial_spX (&(*A)->polys[0]);
		(*A)->polys[0] = deepCopyPolynomial_spX (a);
    } else if (i == 0 && j == 1) {
		freePolynomial_spX (&(*A)->polys[1]);
		(*A)->polys[1] = deepCopyPolynomial_spX (a);
    } else if (i == 1 && j == 0) {
		freePolynomial_spX (&(*A)->polys[2]);
		(*A)->polys[2] = deepCopyPolynomial_spX (a);
    } else {
		freePolynomial_spX (&(*A)->polys[3]);
		(*A)->polys[3] = deepCopyPolynomial_spX (a);
    }
}

// Mulplication of matrix, M, and vector, [U; V]
// [ U ] =  M * [ U ] = [ M00*U + M01*V ]
// [ V ]        [ V ]   [ M10*U + M11*V ]
void mulMat4ToVec_spX_inp (dusmat4_t* M, duspoly_t** u, duspoly_t** v, const Prime_ptr* Pptr)
{
    if (M == NULL) {
		return;
    }
    
    if (isZero_spX (*u) && isZero_spX (*v)) {
		return;
    }

    duspoly_t* M0 = NULL;
    duspoly_t* M1 = NULL;
    duspoly_t* t0 = NULL;
    duspoly_t* t1 = NULL;
    
    mulPolynomialsInForm_spX (M->polys[0], *u, &M0, Pptr);
    mulPolynomialsInForm_spX (M->polys[1], *v, &M1, Pptr);
    
    
    addPolynomials_spX (M0, M1, &t0, Pptr);    
    
    freePolynomial_spX (&M0);
    freePolynomial_spX (&M1);
    
    mulPolynomialsInForm_spX (M->polys[2], *u, &M0, Pptr);
    mulPolynomialsInForm_spX (M->polys[3], *v, &M1, Pptr);

    addPolynomials_spX (M0, M1, &t1, Pptr);

    freePolynomial_spX (&M0);
    freePolynomial_spX (&M1);

    freePolynomial_spX (&*u); 
    freePolynomial_spX (&*v); 

    *u = t0;
    *v = t1;
}

void mulMat4ToVecInForm_spX_inp (dusmat4_t* M, duspoly_t** u, duspoly_t** v, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "Mat*Vec\n");
#endif

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

    // int mulType = 0;
    duspoly_t* M0 = NULL;
    duspoly_t* M1 = NULL;
    duspoly_t* t0 = NULL;
    duspoly_t* t1 = NULL;

    // plainMulPolynomialsInForm_spX (M->polys[0], *u, &M0, Pptr);
    // plainMulPolynomialsInForm_spX (M->polys[1], *v, &M1, Pptr);
    mulPolynomialsInForm_spX (M->polys[0], *u, &M0, Pptr);
    mulPolynomialsInForm_spX (M->polys[1], *v, &M1, Pptr);
    
    
    addPolynomialsInForm_spX (M0, M1, &t0, Pptr);    
    
    freePolynomial_spX (&M0);
    freePolynomial_spX (&M1);
    
    // plainMulPolynomialsInForm_spX (M->polys[2], *u, &M0, Pptr);
    // plainMulPolynomialsInForm_spX (M->polys[3], *v, &M1, Pptr);
	mulPolynomialsInForm_spX (M->polys[2], *u, &M0, Pptr);
    mulPolynomialsInForm_spX (M->polys[3], *v, &M1, Pptr);

    addPolynomialsInForm_spX (M0, M1, &t1, Pptr);

    freePolynomial_spX (&M0);
    freePolynomial_spX (&M1);

    freePolynomial_spX (&*u); 
    freePolynomial_spX (&*v); 

    *u = t0;
    *v = t1;
}

// A*B = B
void mulMat4ToMat4_spX_inp_inp (dusmat4_t* A, dusmat4_t** B, const Prime_ptr* Pptr)
{
    if (A == NULL || *B == NULL) {
		return;
    }
    
    mulMat4ToVec_spX_inp (A, &(*B)->polys[0], &(*B)->polys[2], Pptr);
    mulMat4ToVec_spX_inp (A, &(*B)->polys[1], &(*B)->polys[3], Pptr);
    
    freeMat4_spX (&A);
    A = NULL;
}

// A*B = B
void mulMat4ToMat4InForm_spX_inp_inp (dusmat4_t* A, dusmat4_t** B, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "Mat*Mat\n");
#endif

    if (A == NULL || *B == NULL) {
		return;
    }
    
    mulMat4ToVecInForm_spX_inp (A, &(*B)->polys[0], &(*B)->polys[2], Pptr);
    mulMat4ToVecInForm_spX_inp (A, &(*B)->polys[1], &(*B)->polys[3], Pptr);
    
    freeMat4_spX (&A);
    A = NULL;
}

// A*B = C
void mulMat4ToMat4InForm_spX_inp (dusmat4_t* A, dusmat4_t* B, dusmat4_t** C, const Prime_ptr* Pptr)
{
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "Mat*Mat\n");
#endif

    if (A == NULL || B == NULL) {
		return;
    }

    duspoly_t* C_0= B->polys[0];
    duspoly_t* C_1= B->polys[1];
    duspoly_t* C_2= B->polys[2];
    duspoly_t* C_3= B->polys[3];
   
    mulMat4ToVecInForm_spX_inp (A, &C_0, &C_2, Pptr);
    mulMat4ToVecInForm_spX_inp (A, &C_1, &C_3, Pptr);

    free (B);
  
    dusmat4_t* CC = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    CC->polys[0] = C_0;
    CC->polys[1] = C_1;
    CC->polys[2] = C_2;
    CC->polys[3] = C_3;
    
    *C = CC;
}


void iterHalfGCDMatrixInForm_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, const Prime_ptr* Pptr) 
{
    
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "iHGCD\n");
#endif

    dusmat4_t* E;
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

        // plainDivPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        divPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        swap_spX (&A, &B);

        // plainMulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
        mulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
        subPolynomialsInForm_spX (E->polys[0], tmp1, &tmp2, Pptr);
        freePolynomial_spX (&tmp1);
        freePolynomial_spX (&E->polys[0]);

        E->polys[0] = E->polys[2];
        E->polys[2] = tmp2;
        tmp2 = NULL;

        // plainMulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
        mulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
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

void iterHalfGCDMatrixInForm_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, const Prime_ptr* Pptr) 
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
        fprintf(stderr, "&A = %p , &B = %p \n", A, B); // TEST
        // fprintf(stderr, "deg(B) = %lu, \t m = %d \n", degPolynomial_spX (B), m);
        fprintf(stderr, "before div:\n");
        fprintf(stderr, "A = ");
        printPolynomial_spX (A, Pptr);
        fprintf(stderr, "B = ");
        printPolynomial_spX (B, Pptr);
        fprintf(stderr, "E = \n");
        printMat4OutForm_spX (E, Pptr);
#endif

        // plainDivPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        divPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        swap_spX (&A, &B);

        if (tmp_q == NULL || tmp_q->alloc == 0) {
            swap_spX (&E->polys[0], &E->polys[2]);
            swap_spX (&E->polys[1], &E->polys[3]);            
        } else {

            // plainMulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
            mulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
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
            mulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
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

    *a = A; //deepCopyPolynomial_spX (A);
    *b = B; //deepCopyPolynomial_spX (B);

    // freePolynomial_spX (&A);
    // freePolynomial_spX (&B);

    *M = E;
    return;

}

void _iterHalfGCDMatrixInForm_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, elem_t** lcrems, const Prime_ptr* Pptr) 
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
        fprintf(stderr, "&A = %p , &B = %p \n", A, B); // TEST
        // fprintf(stderr, "deg(B) = %lu, \t m = %d \n", degPolynomial_spX (B), m);
        fprintf(stderr, "before div:\n");
        fprintf(stderr, "A = ");
        printPolynomial_spX (A, Pptr);
        fprintf(stderr, "B = ");
        printPolynomial_spX (B, Pptr);
        fprintf(stderr, "E = \n");
        printMat4OutForm_spX (E, Pptr);
#endif



        // plainDivPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        divPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        if (!isZero_spX(A)) {
            (*lcrems)[A->lt] = A->elems[A->lt];
        }
        swap_spX (&A, &B);

        if (tmp_q == NULL || tmp_q->alloc == 0) {
            swap_spX (&E->polys[0], &E->polys[2]);
            swap_spX (&E->polys[1], &E->polys[3]);            
        } else {

            // plainMulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
            mulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
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
            mulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
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

    *a = A; //deepCopyPolynomial_spX (A);
    *b = B; //deepCopyPolynomial_spX (B);

    // freePolynomial_spX (&A);
    // freePolynomial_spX (&B);

    *M = E;
    return;

}

void _iterHalfGCDMatrixInFormAtMaxK_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, polysize_t k, elem_t** lcrems, const Prime_ptr* Pptr) 
{
    
#if COUNT_SUBPROGRAMS
    fprintf(stderr, "iHGCD\n");
#endif

    dusmat4_t* E;
    identityMat4InForm_spX (&E, Pptr);

    if (isZero_spX (*a) || isZero_spX (*b) || k < 0 || k > MIN_spX((*a)->lt, (*b)->lt)) {
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
        fprintf(stderr, "&A = %p , &B = %p \n", A, B); // TEST
        // fprintf(stderr, "deg(B) = %lu, \t m = %d \n", degPolynomial_spX (B), m);
        fprintf(stderr, "before div:\n");
        fprintf(stderr, "A = ");
        printPolynomial_spX (A, Pptr);
        fprintf(stderr, "B = ");
        printPolynomial_spX (B, Pptr);
        fprintf(stderr, "E = \n");
        printMat4OutForm_spX (E, Pptr);
#endif

        // plainDivPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);
        divPolynomialsInForm_spX_inp (&A, B, &tmp_q, Pptr);

        if (!isZero_spX(A)) {
            (*lcrems)[A->lt] = A->elems[A->lt];
        }
        if (isZero_spX (A) || A->lt <= k) {
            *a = B;
            *b = A;
            return;
        }

        swap_spX (&A, &B);

        if (tmp_q == NULL || tmp_q->alloc == 0) {
            swap_spX (&E->polys[0], &E->polys[2]);
            swap_spX (&E->polys[1], &E->polys[3]);            
        } else {

            // plainMulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
            mulPolynomialsInForm_spX (E->polys[2], tmp_q, &tmp1, Pptr);
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
            mulPolynomialsInForm_spX (E->polys[3], tmp_q, &tmp1, Pptr);
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

    *a = A; //deepCopyPolynomial_spX (A);
    *b = B; //deepCopyPolynomial_spX (B);

    // freePolynomial_spX (&A);
    // freePolynomial_spX (&B);

    *M = E;
    return;

}


void halfGCDMatrixInForm_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDM\n");
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

        iterHalfGCDMatrixInForm_spX_inp (&A0, &B0, M, d, Pptr);

        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);

        return;
    }

    dusmat4_t* RR = NULL;
    long d1 = (d+1) >> 1;

    // corner cases:
    if (d1 >= d) {
        d1 = d-1;
    }  else if (d1 < 1) {
        d1 = 1;
    }

    halfGCDMatrixInForm_spX (A0, B0, &RR, d1, Pptr);
    mulMat4ToVecInForm_spX_inp (RR, &A0, &B0, Pptr);

    if (isZero_spX (B0)) {
        
        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);
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
    divPolynomialsInForm_spX_inp (&A0, B0, &q, Pptr);
    swap_spX (&A0, &B0);

    halfGCDMatrixInForm_spX (A0, B0, &SS, d2, Pptr);

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);

    // plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[1], &RR->polys[3]);
    swap_spX (&RR->polys[3], &tmp2); 

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS); 
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 
}

void _halfGCDMatrixInForm_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, elem_t** lcrems, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDM\n");
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

        _iterHalfGCDMatrixInForm_spX_inp (&A0, &B0, M, d, lcrems, Pptr);

        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);

        return;
    }

    dusmat4_t* RR = NULL;
    long d1 = (d+1) >> 1;

    // corner cases:
    if (d1 >= d) {
        d1 = d-1;
    }  else if (d1 < 1) {
        d1 = 1;
    }

    _halfGCDMatrixInForm_spX (A0, B0, &RR, d1, lcrems, Pptr);
    mulMat4ToVecInForm_spX_inp (RR, &A0, &B0, Pptr);

    if (isZero_spX (B0)) {
        
        freePolynomial_spX (&A0);
        freePolynomial_spX (&B0);
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
    divPolynomialsInForm_spX_inp (&A0, B0, &q, Pptr);
    if (!isZero_spX(A0)) {
        (*lcrems)[A0->lt] = A0->elems[A0->lt];
    }
    swap_spX (&A0, &B0);

    _halfGCDMatrixInForm_spX (A0, B0, &SS, d2, lcrems, Pptr);

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);

    // plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[1], &RR->polys[3]);
    swap_spX (&RR->polys[3], &tmp2); 

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS); 
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 
}



void extHalfGCDMatrixInForm_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, long d, const Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eHGCDM\n");
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
        iterHalfGCDMatrixInForm_spX (a, b, M, d, Pptr);
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

    halfGCDMatrixInForm_spX (a, b, &RR, d1, Pptr);
    mulMat4ToVecInForm_spX_inp (RR, &A0, &B0, Pptr);

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
    divPolynomialsInForm_spX_inp (&A0, B0, &q, Pptr);
    swap_spX (&A0, &B0);

    extHalfGCDMatrixInForm_spX (A0, B0, &SS, d2, Pptr);

    // plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[3], &RR->polys[1]);
    swap_spX (&RR->polys[3], &tmp2);

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS);
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);
}

void extHalfGCDMatrixInForm_spX_inp (duspoly_t** a, duspoly_t** b, dusmat4_t** M, long d, const Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eHGCDMin\n");
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
        iterHalfGCDMatrixInForm_spX_inp (a, b, M, d, Pptr);
        
        // mulMat4ToVecInForm_spX_inp (*M, a, b, Pptr);
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

    halfGCDMatrixInForm_spX (A, B, &RR, d1, Pptr);
    mulMat4ToVecInForm_spX_inp (RR, &A, &B, Pptr);

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

    // plainDivPolynomialsInForm_spX_inp (&A, B, &q, Pptr);
    divPolynomialsInForm_spX_inp (&A, B, &q, Pptr);
    swap_spX (&A, &B);

    extHalfGCDMatrixInForm_spX_inp (&A, &B, &SS, d2, Pptr);

    // plainMulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[2], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[0], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[0], &RR->polys[2]);
    swap_spX (&RR->polys[2], &tmp2);

    freePolynomial_spX (&tmp2);

    // plainMulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    mulPolynomialsInForm_spX (RR->polys[3], q, &tmp1, Pptr);
    subPolynomialsInForm_spX (RR->polys[1], tmp1, &tmp2, Pptr);
    freePolynomial_spX (&tmp1);

    swap_spX (&RR->polys[3], &RR->polys[1]);
    swap_spX (&RR->polys[3], &tmp2);

    freePolynomial_spX (&tmp2);
    freePolynomial_spX (&q);

    mulMat4ToMat4InForm_spX_inp (SS, RR, M, Pptr);

    freeMat4_spX (&E);
    freeMat4_spX (&SS);
    // freeMat4_spX (&RR); // commented: see mulMat4ToMat4InForm_spX_inp 

    *a = A; // TODO: better assignmnet... then delete these lines!
    *b = B;

    // freePolynomial_spX (&A0);
    // freePolynomial_spX (&B0);
}



void halfGCDInForm_spX (duspoly_t* a, duspoly_t* b, duspoly_t** ap, duspoly_t** bp, const Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCD\n");
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

    dusmat4_t* M = NULL;
    halfGCDMatrixInForm_spX (A, B, &M, d1, Pptr);
    mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
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

    // duspoly_t* q;
    // plainDivPolynomials_spX_inp (&A, B, &q, Pptr);
    // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
    remPolynomialsInForm_spX_inp (&A, B, Pptr);
    swap_spX (&A, &B);

    halfGCDMatrixInForm_spX (A, B, &M, d2, Pptr);
    mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
    freeMat4_spX (&M);

    *ap = A;
    *bp = B; 
}

void halfGCDInForm_spX_inp (duspoly_t** a, duspoly_t** b, const Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDin\n");
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

    dusmat4_t* M = NULL;
    halfGCDMatrixInForm_spX (A, B, &M, d1, Pptr);
    mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
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

    // duspoly_t* q;
    // plainDivPolynomials_spX_inp (&A, B, &q, Pptr);
    // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
    remPolynomialsInForm_spX_inp (&A, B, Pptr);
    swap_spX (&A, &B);
    halfGCDMatrixInForm_spX (A, B, &M, d2, Pptr);
    mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
    freeMat4_spX (&M);

    *a = A;
    *b = B; 
}

void _halfGCDInForm_spX_inp (duspoly_t** a, duspoly_t** b, elem_t** lcrems, const Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDin\n");
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

    dusmat4_t* M = NULL;
    // _halfGCDMatrixInForm_spX (A, B, &M, d1, lcrems, Pptr);
    // mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
      _iterHalfGCDMatrixInForm_spX_inp (&A, &B, &M, d1, lcrems, Pptr); // TODO: fastResultant
    freeMat4_spX (&M);
    if (!isZero_spX (B)) {
        // (*lcrems)[B->lt] = B->elems[B->lt];
    } else {
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

    // duspoly_t* q;
    // plainDivPolynomials_spX_inp (&A, B, &q, Pptr);
    // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
    remPolynomialsInForm_spX_inp (&A, B, Pptr);
    if (!isZero_spX(A)) {
        (*lcrems)[A->lt] = A->elems[A->lt];
    }
    swap_spX (&A, &B);
    // _halfGCDMatrixInForm_spX (A, B, &M, d2, lcrems, Pptr);
    // mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
      _iterHalfGCDMatrixInForm_spX_inp (&A, &B, &M, d2, lcrems, Pptr); // TODO: fastResultant
    // if (!isZero_spX) {
    //     // TODO: lcrems...
    // }
    freeMat4_spX (&M);

    *a = A;
    *b = B; 
}


void _halfGCDInFormAtMaxK_spX_inp (duspoly_t** a, duspoly_t** b, polysize_t k, elem_t** lcrems, const Prime_ptr* Pptr) 
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "HGCDin\n");
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

    if (k < 0 || k >= MIN_spX(A->lt, B->lt)) {
        return;
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
    
    dusmat4_t* M = NULL;
    // _halfGCDMatrixInForm_spX (A, B, &M, d1, lcrems, Pptr);
    // mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
      _iterHalfGCDMatrixInFormAtMaxK_spX_inp (&A, &B, &M, d1, k, lcrems, Pptr); // TODO: fastResultant
    //   fprintf (stderr, "M1 := "); // TEST
    //   printMat4OutForm_spX (M, Pptr); // TEST
    freeMat4_spX (&M);

    if (B->lt <= k) {
        *a = A;
        *b = B;
        return;
    }

    if (!isZero_spX (B)) {
        // (*lcrems)[B->lt] = B->elems[B->lt];
    } else {
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

    // duspoly_t* q;
    // plainDivPolynomials_spX_inp (&A, B, &q, Pptr);
    // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
    remPolynomialsInForm_spX_inp (&A, B, Pptr);
    if (!isZero_spX(A)) {
        (*lcrems)[A->lt] = A->elems[A->lt];
    }
    swap_spX (&A, &B);

    if (B->lt <= k) {
        *a = A;
        *b = B;
        return;
    }
    // _halfGCDMatrixInForm_spX (A, B, &M, d2, lcrems, Pptr);
    // mulMat4ToVecInForm_spX_inp (M, &A, &B, Pptr);
      _iterHalfGCDMatrixInFormAtMaxK_spX_inp (&A, &B, &M, d2, k, lcrems, Pptr); // TODO: fastResultant
    //   fprintf (stderr, "M2 := "); // TEST
    //   printMat4OutForm_spX (M, Pptr); // TEST
    freeMat4_spX (&M);
    *a = A;
    *b = B; 
}


void GCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "GCD\n");
#endif

    if (g == NULL){
        return;
    }

    if (degPolynomial_spX (a) < PLAINGCD_CROSSOVER ||
        degPolynomial_spX (b) < PLAINGCD_CROSSOVER) {
        plainGCDInForm_spX (a, b, g, Pptr);
        return;
    } 

    duspoly_t* A;
    duspoly_t* B;

    if (degPolynomial_spX (a) == degPolynomial_spX (b)) {
        B = NULL;
        // plainRemPolynomialsInForm_spX (a, b, &B, Pptr);
        remPolynomialsInForm_spX (a, b, &B, Pptr);
        A = deepCopyPolynomial_spX (b);
    } else if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    } else {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    }

    while (!isZero_spX (B) && degPolynomial_spX (A) > PLAINGCD_CROSSOVER) {
        halfGCDInForm_spX_inp (&A, &B, Pptr);

        if (!isZero_spX (B)) {
            // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
            remPolynomialsInForm_spX_inp (&A, B, Pptr);
            swap_spX (&A, &B);
        }
    }

    plainGCDInForm_spX (A, B, g, Pptr);

    freePolynomial_spX (&A);
    freePolynomial_spX (&B);
}

void _GCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** g, elem_t** lcrems, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "GCD\n");
#endif

    if (g == NULL){
        return;
    }

    if (degPolynomial_spX (a) < PLAINGCD_CROSSOVER ||
        degPolynomial_spX (b) < PLAINGCD_CROSSOVER) {
        _plainGCDInForm_spX (a, b, g, lcrems, Pptr);
        return;
    } 

    duspoly_t* A;
    duspoly_t* B;
    polysize_t min_deg = MIN_spX (a->lt, b->lt);
    elem_t* lcRems = (elem_t*) calloc (min_deg, sizeof(elem_t));

    if (degPolynomial_spX (a) == degPolynomial_spX (b)) {
        B = NULL;
        // plainRemPolynomialsInForm_spX (a, b, &B, Pptr);
        remPolynomialsInForm_spX (a, b, &B, Pptr);
        if (!isZero_spX (B)){
            lcRems[B->lt] = B->elems[B->lt];
        }
        A = deepCopyPolynomial_spX (b);
    } else if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    } else {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    }

    while (!isZero_spX (B) && degPolynomial_spX (A) > PLAINGCD_CROSSOVER) {
        _halfGCDInForm_spX_inp (&A, &B, &lcRems, Pptr);
        if (!isZero_spX (B)) {
            // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
            remPolynomialsInForm_spX_inp (&A, B, Pptr);
            if (!isZero_spX (A)) {
                lcRems[A->lt] = A->elems[A->lt];
            }
            swap_spX (&A, &B);
        }
    }

    _plainGCDInForm_spX (A, B, g, &lcRems, Pptr);
    *lcrems = lcRems;
    freePolynomial_spX (&A);
    freePolynomial_spX (&B);
}


void _HGCDInFormAtMaxK_spX (const duspoly_t* a, const duspoly_t* b, polysize_t k, duspoly_t** rk1, duspoly_t** rk, elem_t** lcrems, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "GCD\n");
#endif

    if (rk1 == NULL || rk == NULL){
        return;
    }

    // fprintf (stderr, "call _HGCDInFormAtMaxk... deg(a) = %ld and deg(b) = %ld\n", a->lt, b->lt);

    if (degPolynomial_spX (a) < PLAINGCD_CROSSOVER ||
        degPolynomial_spX (b) < PLAINGCD_CROSSOVER) {
        _plainGCDInFormAtMaxK_spX (a, b, k, rk1, rk, lcrems, Pptr);

        // for (int i = 0; i < MIN_spX (a->lt, b->lt);i++) {
        //     fprintf (stderr, "lcrems[%d]=%lld\n", i, (*lcrems)[i]);
        // } // TEST
        // printPolynomialOutForm_spX (*rk1, Pptr); // TEST
        // printPolynomialOutForm_spX (*rk, Pptr); // TEST

        return;
    } 

    if (isZero_spX (a) || isZero_spX (b)) {
        return;
    }

    polysize_t min_deg = MIN_spX (a->lt, b->lt);
    if (k < 0 || k > min_deg) {
        return;
    }

    if (k == min_deg) {
        if (a->lt >= b->lt) {
            *rk1 = deepCopyPolynomial_spX (a);
            *rk = deepCopyPolynomial_spX (b);
        } else {
            *rk1 = deepCopyPolynomial_spX (b);
            *rk = deepCopyPolynomial_spX (a);
        }
    }

    duspoly_t* A;
    duspoly_t* B;
    elem_t* lcRems = (elem_t*) calloc (min_deg, sizeof(elem_t));

    if (degPolynomial_spX (a) == degPolynomial_spX (b)) {
        B = NULL;
        // plainRemPolynomialsInForm_spX (a, b, &B, Pptr);
        remPolynomialsInForm_spX (a, b, &B, Pptr);
        if (!isZero_spX (B)){
            lcRems[B->lt] = B->elems[B->lt];
            if (B->lt <= k) {
                *rk1 = deepCopyPolynomial_spX (b);
                *rk = B;
                *lcrems = lcRems;
                return;
            }
        }
        A = deepCopyPolynomial_spX (b);
    } else if (degPolynomial_spX (a) < degPolynomial_spX (b)) {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    } else {
        A = deepCopyPolynomial_spX (a);
        B = deepCopyPolynomial_spX (b);
    }

    while (!isZero_spX (B) && degPolynomial_spX (A) > PLAINGCD_CROSSOVER) {
        _halfGCDInFormAtMaxK_spX_inp (&A, &B, k, &lcRems, Pptr);
        if (!isZero_spX (B)) {
            if (B->lt <= k) {
                *rk1 = A;
                *rk = B;
                *lcrems = lcRems;
                return;
            }
            // plainRemPolynomialsInForm_spX_inp (&A, B, Pptr);
            remPolynomialsInForm_spX_inp (&A, B, Pptr);
            if (!isZero_spX (A)) {
                lcRems[A->lt] = A->elems[A->lt];
                if (A->lt <= k) {
                    *rk1 = B;
                    *rk = A;
                    *lcrems = lcRems;
                    return;
                }
            }
            swap_spX (&A, &B);
        } else {
            *rk1 = A;
            *rk = NULL;
            *lcrems = lcRems;
            return;
        }
    }

    _plainGCDInFormAtMaxK_spX (A, B, k, rk1, rk, &lcRems, Pptr);
    *lcrems = lcRems;
    freePolynomial_spX (&A);
    freePolynomial_spX (&B);
}


void extGCDInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "eGCD\n");
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
        divPolynomialsInForm_spX (a, b, &B, &q, Pptr);
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

    extHalfGCDMatrixInForm_spX_inp (&A, &B, &M, d, Pptr);

    freePolynomial_spX (&B);

    if (isEqualDeg) {
        uu = M->polys[1];
        
        // plainMulPolynomialsInForm_spX (q, M->polys[1], &tmp, Pptr);
        mulPolynomialsInForm_spX (q, M->polys[1], &tmp, Pptr);
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

    if (g != NULL){
        *g = gg;
    } else {
        freePolynomial_spX (&gg);
    }

    if (u != NULL){
        *u = uu;
    } else {
        freePolynomial_spX (&uu);
    }

    if (v != NULL){
        *v = vv;
    } else {
        freePolynomial_spX (&vv);
    }

    free (M);
}


void YapHalfGCDMatrixInForm_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "YHGCDM\n");
#endif

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
    
    if (deg_a > deg_b) {
		A = deepCopyPolynomial_spX (a);
		B = deepCopyPolynomial_spX (b);
    } else if (deg_a <= deg_b) {
		/* fprintf (stderr, "DUSP Alert: In HalfGCD, deg(a) <= deg(b)!\n"); // TEST */
		*M = E; //TODO:
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
    
    if (deg_b < m) {
		*M = E;
		return;
    }

    duspoly_t* A0 = NULL;
    duspoly_t* B0 = NULL;
    dusmat4_t* RR = NULL;
    
    rightShiftPolynomial_spX (A, &A0, m);
    rightShiftPolynomial_spX (B, &B0, m);
    
    YapHalfGCDMatrixInForm_spX (A0, B0, &RR, Pptr);  // first recursive call

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);
    
   // // TEST:
   //  fprintf(stderr, "0RR = \n");
   //  printMat4_spX (RR, Pptr);

    mulMat4ToVecInForm_spX_inp (RR, &A, &B, Pptr); // A = Ap = r_{j-1} , B = Bp = r_{j}


    // // TEST:
    // fprintf(stderr, "1RR = \n");
    // printMat4_spX (RR, Pptr);
    
    if (degPolynomial_spX (B) < m) { // isZero_spX (B)

		freePolynomial_spX (&A); 
		freeMat4_spX (&E);
	
		*M = RR;
		return;
    }

    duspoly_t* q = NULL;
    duspoly_t* r = NULL;
    duspoly_t* D = NULL;
    
    plainDivPolynomialsInForm_spX (A, B, &r, &q, Pptr); // TODO : Fast Monic Div
    
    freePolynomial_spX (&A);
    
    if (isZero_spX (r)) { 
	
		negPolynomialInForm_spX (q, Pptr);
	
		dusmat4_t* TM = (dusmat4_t*) malloc (sizeof (dusmat4_t));
		dusmat4_t* MM = NULL;
		TM->polys[0] = NULL;
		TM->polys[1] = constPolynomial_spX (1, Pptr);
		TM->polys[2] = constPolynomial_spX (1, Pptr);
		TM->polys[3] = q;	
	
      // // TEST:
      //   fprintf(stderr, "TM = \n");
      //   printMat4_spX (TM, Pptr);
        

      // // TEST:
      //   fprintf(stderr, "RR = \n");
      //   printMat4_spX (RR, Pptr);
        

		// mulMat4ToMat4InForm_spX_inp_inp (RR, &TM, Pptr);
		mulMat4ToMat4InForm_spX_inp (TM, RR, &MM, Pptr);

      // // TEST:
      //   fprintf(stderr, "MM = \n");
      //   printMat4_spX (MM, Pptr);
        
        
		freeMat4_spX (&TM);
		freeMat4_spX (&E);
		freePolynomial_spX (&B);
		
		*M = MM;
		return;
    }

    elem_t lc = leadingCoeffInForm_spX (r); // IT SHOULD BE UNCOMMENT
    
    lc = smallprimefield_inv (lc, Pptr); // IT SHOULD BE UNCOMMENT

    scalarMulPolynomialInForm_spX (r, lc, &D, Pptr);

    freePolynomial_spX (&r);

    duspoly_t* tq = NULL;
    
    scalarMulPolynomialInForm_spX (q, smallprimefield_sub (0, lc, Pptr), &tq, Pptr);

    freePolynomial_spX (&q);    
    
    dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = tq;
        
    polysize_t k = 2*m - degPolynomial_spX (B);    
    duspoly_t* C0 = NULL;
    duspoly_t* D0 = NULL;
    dusmat4_t* S  = NULL;
    
    rightShiftPolynomial_spX (B, &C0, k);
    rightShiftPolynomial_spX (D , &D0, k);

    freePolynomial_spX (&B);
    freePolynomial_spX (&D);

    YapHalfGCDMatrixInForm_spX (C0, D0, &S, Pptr); // second recursive call
    
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);
    
    /* mulMat4ToMat4InForm_spX_inp_inp (RR, &Q, Pptr); */
    /* mulMat4ToMat4InForm_spX_inp_inp (Q, &S, Pptr); */

    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;


    // // TEST:
    // fprintf(stderr, "S = \n");
    // printMat4_spX (S, Pptr);

    // // TEST:
    // fprintf(stderr, "RR = \n");
    // printMat4_spX (RR, Pptr);

    // // TEST:
    // fprintf(stderr, "Q = \n");
    // printMat4_spX (Q, Pptr);


    mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr);
    mulMat4ToMat4InForm_spX_inp (S, QRR, &SQRR, Pptr);

    // // TEST:
    // fprintf(stderr, "QRR = \n");
    // printMat4_spX (QRR, Pptr);

  
    freeMat4_spX (&Q);
    freeMat4_spX (&S);
    
    // // TEST:
    // fprintf(stderr, "SQRR = \n");
    // printMat4_spX (SQRR, Pptr);

    *M = SQRR;

    return;
}


void YapHalfGCDMatrixInForm_wPSInv_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, const Prime_ptr* Pptr)
{
    

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "YHGCDMwPS\n");
#endif

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
    
    if (deg_a > deg_b) {
		A = deepCopyPolynomial_spX (a);
		B = deepCopyPolynomial_spX (b);
    } else if (deg_a <= deg_b) {
		/* fprintf (stderr, "DUSP Alert: In HalfGCD, deg(a) <= deg(b)!\n"); // TEST */
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
    
    if (deg_b < m) {
		*M = E;
		return;
    }

    duspoly_t* A0 = NULL;
    duspoly_t* B0 = NULL;
    dusmat4_t* RR = NULL;
    
    rightShiftPolynomial_spX (A, &A0, m);
    rightShiftPolynomial_spX (B, &B0, m);
    
    YapHalfGCDMatrixInForm_wPSInv_spX (A0, B0, &RR, Pptr);  // first recursive call

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);

    mulMat4ToVecInForm_spX_inp (RR, &A, &B, Pptr); // A = Ap = r_{j-1} , B = Bp = r_{j}
    
    if (degPolynomial_spX (B) < m) { // isZero_spX (B)

		freePolynomial_spX (&A); 
		freeMat4_spX (&E);
	
		*M = RR;
		return;
    }

    duspoly_t* q = NULL;
    duspoly_t* r = NULL;
    duspoly_t* D = NULL;
    
    fastDivPolynomialInForm_wPSInv_spX (A, B, &r, &q, Pptr);
    
    freePolynomial_spX (&A);
    
    if (isZero_spX (r)) { 
	
		negPolynomialInForm_spX (q, Pptr);
	
		dusmat4_t* TM = (dusmat4_t*) malloc (sizeof (dusmat4_t));
		dusmat4_t* MM = NULL;
		TM->polys[0] = NULL;
		TM->polys[1] = constPolynomial_spX (1, Pptr);
		TM->polys[2] = constPolynomial_spX (1, Pptr);
		TM->polys[3] = q;	
	
		// mulMat4ToMat4InForm_spX_inp_inp (RR, &TM, Pptr);
		mulMat4ToMat4InForm_spX_inp (TM, RR, &MM, Pptr);

		freeMat4_spX (&TM);
		freeMat4_spX (&E);
		freePolynomial_spX (&B);
	
		*M = MM;
		return;
    }

    elem_t lc = leadingCoeffInForm_spX (r); // IT SHOULD BE UNCOMMENT
    
    lc = smallprimefield_inv (lc, Pptr); // IT SHOULD BE UNCOMMENT

    scalarMulPolynomialInForm_spX (r, lc, &D, Pptr);
    
    freePolynomial_spX (&r);

    duspoly_t* tq = NULL;
    
    scalarMulPolynomialInForm_spX (q, smallprimefield_sub (0, lc, Pptr), &tq, Pptr);
        
    freePolynomial_spX (&q);
    
    
    dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = tq;
        
    polysize_t k = 2*m - degPolynomial_spX (B);
    
    duspoly_t* C0 = NULL;
    duspoly_t* D0 = NULL;
    dusmat4_t* S  = NULL;
    
    rightShiftPolynomial_spX (B, &C0, k);
    rightShiftPolynomial_spX (D , &D0, k);

    freePolynomial_spX (&B);
    freePolynomial_spX (&D);

    YapHalfGCDMatrixInForm_wPSInv_spX (C0, D0, &S, Pptr); // second recursive call
	
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);

    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;
    mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr);

    mulMat4ToMat4InForm_spX_inp (S, QRR, &SQRR, Pptr);
    
    freeMat4_spX (&Q);
    freeMat4_spX (&S);
    
    *M = SQRR;
    
    return;
}


void GCDMatrixInForm_wHGCD_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "YGCDM\n");
#endif
 
    dusmat4_t* RR = NULL;
    YapHalfGCDMatrixInForm_spX (a, b, &RR, Pptr);
    
    // TODO : implement non in-place mulMat4ToVecInForm 
    duspoly_t* a0 = deepCopyPolynomial_spX (a);
    duspoly_t* b0 = deepCopyPolynomial_spX (b);
    
    mulMat4ToVecInForm_spX_inp (RR, &a0, &b0, Pptr);
    
    if (isZero_spX (b0)) {
		*M = RR;
		return;
    }

    duspoly_t* r   = NULL;
    duspoly_t* q   = NULL;
    duspoly_t* lcr = NULL;
    duspoly_t* lcq = NULL;
    dusmat4_t* Q   = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    
    plainDivPolynomialsInForm_spX (a0, b0, &r, &q, Pptr); // TODO: Fast Monic Divid
    negPolynomialInForm_spX (q, Pptr); // q = -q;
    
    freePolynomial_spX (&a0);

    if (isZero_spX (r)) {
	
		freePolynomial_spX (&b0);
	
		Q->polys[0] = NULL;
		Q->polys[1] = constPolynomial_spX (1, Pptr);
		Q->polys[2] = constPolynomial_spX (1, Pptr);
		Q->polys[3] = q; // q = -q

		dusmat4_t* QRR = NULL;
	
		// mulMat4ToMat4InForm_spX_inp_inp (Q, &RR, Pptr);
		mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr); 
	
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
    
    dusmat4_t* S = NULL;
    GCDMatrixInForm_wHGCD_spX (b0, lcr, &S, Pptr); 

    freePolynomial_spX (&b0);
    freePolynomial_spX (&lcr);

    dusmat4_t* RRQ = NULL;
    dusmat4_t* SRRQ = NULL;
    
    mulMat4ToMat4InForm_spX_inp (Q, RR, &RRQ, Pptr);
    mulMat4ToMat4InForm_spX_inp (S, RRQ, &SRRQ, Pptr);

    freeMat4_spX (&Q);
    freeMat4_spX (&S);
    
    *M = SRRQ;     
}

void GCDInForm_wHGCD_spX (duspoly_t* a, duspoly_t* b, duspoly_t** g, const Prime_ptr* Pptr)
{


#if COUNT_SUBPROGRAMS
    fprintf(stderr, "YGCD\n");
#endif

	if (isZero_spX (a) || isZero_spX (b)) {
		plainGCDInForm_spX (a, b, g, Pptr);
		return;
	}

	if (degPolynomial_spX (a) < PLAINGCD_CROSSOVER || degPolynomial_spX(b) < PLAINGCD_CROSSOVER) {
		plainGCDInForm_spX (a, b, g, Pptr);
		return;
	}

	dusmat4_t* M = NULL;
	duspoly_t* G = deepCopyPolynomial_spX (a);
	duspoly_t* T = deepCopyPolynomial_spX (b);
	GCDMatrixInForm_wHGCD_spX (a, b, &M, Pptr);

	mulMat4ToVecInForm_spX_inp (M, &G, &T, Pptr);

	freePolynomial_spX (&T);
	freeMat4_spX (&M);

	elem_t lc = 0;
	monicPolynomialInForm_spX (G, g, &lc, Pptr);

	freePolynomial_spX (&G);
}

void ExtGCDInForm_wHGCD_spX (duspoly_t* a, duspoly_t* b, duspoly_t** u, duspoly_t** v, duspoly_t** g, const Prime_ptr* Pptr)
{

#if COUNT_SUBPROGRAMS
    fprintf(stderr, "YeGCD\n");
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
	
		GCDMatrixInForm_wHGCD_spX (a, b, &M, Pptr);
	
		mulPolynomialsInForm_spX (M->polys[0], a, &M00a, Pptr);
		mulPolynomialsInForm_spX (M->polys[1], b, &M01b,  Pptr);

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

		plainRemPolynomials_spX (a, b, &r, Pptr); // r = a mod b
		ExtGCDInForm_wHGCD_spX (b, r, v, u, g, Pptr);

		freePolynomial_spX (&r);
	
		return;
	
    } else {
		ExtGCDInForm_wHGCD_spX (b, a, v, u, g, Pptr);
    }
}


void halfGCDMatrix_NTL_spX (duspoly_t* u, duspoly_t* v, dusmat4_t** M, polysize_t m, const Prime_ptr* Pptr)
{

    polysize_t deg_u = degPolynomial_spX (u);
    polysize_t deg_v = degPolynomial_spX (v);
    
    if (isZero_spX (v) || m <= deg_u - deg_v) {
		identityMat4InForm_spX (M, Pptr);
		return;
    }

    long n = deg_u - 2*(m-1);
    
    // n should be non-negative 
    if (n < 0) {
		n = 0;
    }

    dusmat4_t* M_out = NULL;
    
    duspoly_t* u0 = NULL;
    duspoly_t* v0 = NULL;
    /* duspoly_t* rem = NULL; */
    duspoly_t* quo = NULL;

    duspoly_t* t = NULL;
    duspoly_t* tmp1 = NULL;
    duspoly_t* tmp2 = NULL;
    
    rightShiftPolynomial_spX (u, &u0, n);
    rightShiftPolynomial_spX (v, &v0, n);
    
    /* fprintf (stderr, "rightShiftPolynomial_spX (u, &u0, n) = "); */
    /* printPolynomial_spX (u0, Pptr); */
    
    /* fprintf (stderr, "rightShiftPolynomial_spX (v, &v0, n) = "); */
    /* printPolynomial_spX (v0, Pptr);  */
    
    
    // iterative halfGCDMatrix: 
    if (m <= HALFGCD_CROSSOVER) {
		identityMat4InForm_spX (&M_out, Pptr);
		polysize_t deg_u0 = degPolynomial_spX (u0);
        polysize_t deg_v0 = degPolynomial_spX (v0);
	
		long bd = deg_u0 - m;
	
		if (deg_v0 <= bd) {
			*M = M_out;
			return;
		}
	
		while (deg_v0 > bd) {
			plainDivPolynomials_spX_inp (&u0, v0, &quo, Pptr);
			t = u0;
			u0 = v0;
			v0 = t;
	    
			deg_v0 = degPolynomial_spX (v0); // next step
	    
			// set M[0,0] and M[1,0]:
			mulPolynomialsInForm_spX (quo, M_out->polys[2], &tmp1, Pptr);
			subPolynomials_spX (M_out->polys[0], tmp1, &tmp2, Pptr);
	    
			freePolynomial_spX (&tmp1);
			freePolynomial_spX (&M_out->polys[0]);
	    
			M_out->polys[0] = M_out->polys[2];
			M_out->polys[2] = tmp2; // deepCopyPolynomial_spX (tmp2);
	    
			// freePolynomial_spX (&tmp2);

			// set M[0,1] and M[1,1]:
			mulPolynomialsInForm_spX (quo, M_out->polys[3], &tmp1, Pptr);
			subPolynomials_spX (M_out->polys[1], tmp1, &tmp2, Pptr);

			freePolynomial_spX (&tmp1);
			freePolynomial_spX (&M_out->polys[1]);

			M_out->polys[1] = M_out->polys[3];
			M_out->polys[3] = tmp2; //deepCopyPolynomial_spX (tmp2);

			// freePolynomial_spX (&tmp2);
		}
	
		*M = M_out;
		return;
    }

    dusmat4_t* M1 = NULL;
    long m1 = (m+1)>>1; 
    
    if (m1 < 1) {
		m1 = 1;
    }
    
    if (m1 >= m) {
		m1 = m-1;
    }
    
    // recursively call halfGCDMatrix for new threshold m1
    halfGCDMatrix_NTL_spX (u0, v0, &M1, m1, Pptr); 
    mulMat4ToVec_spX_inp (M1, &u0, &v0, Pptr);
    
    dusmat4_t* M2 = NULL;
    long m2 = degPolynomial_spX (v0) - deg_u + n + m; // next threshold 
    
    if (isZero_spX (v0) || m2 <= 0) {
		*M = M1;
		return;
    }
    
    plainDivPolynomials_spX_inp (&u0, v0, &quo, Pptr);
    t  = u0;
    u0 = v0;
    v0 = t; // deepCopyPolynomial_spX (rem);

    // freePolynomial_spX (&rem); // TODO : maybe

    // recusively call halfGCDMatrix for threshold m2:
    halfGCDMatrix_NTL_spX (u0, v0, &M2, m2, Pptr); 
    
    // set M[0,0] and M[1,0]: 
    mulPolynomialsInForm_spX (quo, M1->polys[2], &tmp1, Pptr);
    subPolynomials_spX (M1->polys[0], tmp1, &tmp2, Pptr);

    freePolynomial_spX (&tmp1);
    freePolynomial_spX (&M1->polys[0]);
    
    M1->polys[0] = M1->polys[2];
    M1->polys[2] = tmp2; // deepCopyPolynomial_spX (tmp2);

    // freePolynomial_spX (&tmp2);

    // set M[0,1] and M[1,1]:
    mulPolynomialsInForm_spX (quo, M1->polys[3], &tmp1, Pptr);
    subPolynomials_spX (M1->polys[1], tmp1, &tmp2, Pptr);

    freePolynomial_spX (&tmp1);
    freePolynomial_spX (&M1->polys[1]);
    
    M1->polys[1] = M1->polys[3];
    M1->polys[3] = tmp2; // deepCopyPolynomial_spX (tmp2);

    // freePolynomial_spX (&tmp2);

    mulMat4ToMat4_spX_inp_inp (M1, &M2, Pptr);

    *M = M2;
}


void halfGCD_spX (duspoly_t* u, duspoly_t* v, duspoly_t** up, duspoly_t** vp, const Prime_ptr* Pptr)
{
    // for non-normalized inputs:
    /* normalize_spX (&u); */
    /* normalize_spX (&v); */

    polysize_t deg_u = degPolynomial_spX (u);
    polysize_t deg_v = degPolynomial_spX (v);
    polysize_t m = (deg_u+1)>>1; // ceil(deg_u/2)

    if (isZero_spX (v) || m <= deg_u - deg_v) {
		*up = deepCopyPolynomial_spX (u);
		*vp = deepCopyPolynomial_spX (v);
		return;
    }

    long m1 = (m+1)>>1;

    if (m1 < 1) {
		m1 = 1;
    }

    if (m <= m1) {
		m1 = m-1;
    }
    
    dusmat4_t* M1 = NULL;
    /* duspoly_t* rem = NULL; */
    duspoly_t* u0 = deepCopyPolynomial_spX (u);
    duspoly_t* v0 = deepCopyPolynomial_spX (v);
    duspoly_t* t = NULL;
    
    halfGCDMatrix_NTL_spX (u, v, &M1, m1, Pptr);
    mulMat4ToVec_spX_inp (M1, &u0, &v0, Pptr);
    
    freeMat4_spX (&M1);
    
    long m2 = degPolynomial_spX (v0) - deg_u + m;

    if (isZero_spX (v0) || m2 <= 0) {
		*up = u0;
		*vp = v0;
		return;
    }
        
    plainRemPolynomials_spX_inp (&u0, v0, Pptr);
    t = u0;
    u0 = v0;
    v0 = t;
    
    halfGCDMatrix_NTL_spX (u0, v0, &M1, m2, Pptr);
    mulMat4ToVec_spX_inp (M1, &u0, &v0, Pptr);

    freeMat4_spX (&M1);

    *up = u0;
    *vp = v0;
}

void halfGCD_spX_inp (duspoly_t** u, duspoly_t** v, const Prime_ptr* Pptr)
{
    // for non-normalized inputs:
    /* normalize_spX (&u); */
    /* normalize_spX (&v); */

    polysize_t deg_u = degPolynomial_spX (*u);
    polysize_t deg_v = degPolynomial_spX (*v);
    polysize_t m = (deg_u+1)>>1; // ceil(deg_u/2)
    
    if (isZero_spX (*v) || m <= deg_u - deg_v) {
		return;
    }
    
    long m1 = (m+1)>>1;
    
    if (m1 < 1) {
		m1 = 1;
    }

    if (m <= m1) {
		m1 = m-1;
    }
    
    dusmat4_t* M1 = NULL;
    /* duspoly_t* rem = NULL; */
    duspoly_t* u0 = *u;
    duspoly_t* v0 = *v;
    duspoly_t* t = NULL;
    
    halfGCDMatrix_NTL_spX (u0, v0, &M1, m1, Pptr);
    mulMat4ToVec_spX_inp (M1, &u0, &v0, Pptr);
    
    freeMat4_spX (&M1);
    
    long m2 = degPolynomial_spX (v0) - deg_u + m;

    if (isZero_spX (v0) || m2 <= 0) {
		return;
    }
    
    plainRemPolynomials_spX_inp (&u0, v0, Pptr);
    t = u0;
    u0 = v0;
    v0 = t;
    
    halfGCDMatrix_NTL_spX (u0, v0, &M1, m2, Pptr);
    mulMat4ToVec_spX_inp (M1, &u0, &v0, Pptr);

    freeMat4_spX (&M1);
}


/****************************
 ** Resultant and GCD 
 ****************************/

void sylvResultantInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** res, const Prime_ptr* Pptr)
{
	if (isZero_spX (a)) {
		if (isZero_spX (b) || b->lt) {
            return;
		} else {
			*res = constPolynomial_spX (1, Pptr);	
            return;
		}
	} else if (isZero_spX (b)) {
		if (a->lt) {
			return;
		} else {
			*res = constPolynomial_spX (1, Pptr);
            return;
		}
	}

	if (!a->lt && !b->lt) {
        *res = constPolynomial_spX (1, Pptr);
        return;
	}

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_ab = deg_a * deg_b;
    duspoly_t* res_t = NULL;
    
    if (deg_b == 0) {
		*res = constPolynomialInForm_spX (smallprimefield_exp (leadingCoeffInForm_spX (b), deg_a, Pptr), Pptr);
		return;
    }

    if (deg_a < deg_b) {
        sylvResultantInForm_spX (b, a, &res_t, Pptr);
        if (deg_ab & 1) {
            scalarMulPolynomialInForm_spX_inp (&res_t, smallprimefield_convert_in (-1, Pptr), Pptr);    
        }
        *res = res_t;
		return;
    }

    duspoly_t* r = NULL;
    // plainRemPolynomialsInForm_spX (a, b, &r, Pptr); // TODO: fastDiv
    remPolynomialsInForm_spX (a, b, &r, Pptr); 

    if (isZero_spX (r)) {
		*res = NULL;
		return;
    }

    elem_t coef = smallprimefield_exp (leadingCoeffInForm_spX (b), deg_a - degPolynomial_spX (r), Pptr);
    // fprintf (stderr, "[TEST] lc_b := %lld, deg(a) = %ld, deg(r) = %ld\n", leadingCoeffInForm_spX (b), deg_a, degPolynomial_spX (r));  // TEST 
    // fprintf (stderr, "[TEST] lc_h^(m-p) := %lld\n", coef);  // TEST 

    if (deg_ab & 1) {    
        coef = smallprimefield_mul (coef, smallprimefield_convert_in (-1, Pptr), Pptr);
    }

    sylvResultantInForm_spX (b, r, &res_t, Pptr);
    scalarMulPolynomialInForm_spX_inp (&res_t, coef, Pptr);
	freePolynomial_spX (&r);

    *res = res_t;
    return;    
}

void hgcdResultantInForm_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** res, const Prime_ptr* Pptr)
{
	if (isZero_spX (a)) {
		if (isZero_spX (b) || b->lt) {
            return;
		} else {
			*res = constPolynomial_spX (1, Pptr);	
            return;
		}
	} else if (isZero_spX (b)) {
		if (a->lt) {
			return;
		} else {
			*res = constPolynomial_spX (1, Pptr);
            return;
		}
	}

	if (!a->lt && !b->lt) {
        *res = constPolynomial_spX (1, Pptr);
        return;
	}

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_ab = deg_a * deg_b;
    duspoly_t* res_t = NULL;
    // elem_t mone;

	duspoly_t* tmp = NULL;

    if (deg_a < deg_b) {
        hgcdResultantInForm_spX (b, a, &res_t, Pptr);
        if (deg_ab & 1) {
            scalarMulPolynomialInForm_spX_inp (&res_t, smallprimefield_convert_in (-1, Pptr), Pptr);    
        }
        *res = res_t;
		return;
    }

    elem_t* lcrems = NULL;
    polysize_t rems_size = MIN_spX (a->lt, b->lt);
    _GCDInForm_spX (a, b, &res_t, &lcrems, Pptr);

    if (res_t != NULL && res_t->lt > 0) {
        if (lcrems != NULL) {
            free (lcrems);
        }
        freePolynomial_spX (&res_t);
        *res = NULL;
        return;
    }
    freePolynomial_spX (&res_t);

    if (lcrems == NULL) {
        fprintf (stderr, "DUSP Error: lcrems couldn't be NULL in hgcdResultantInForm.\n");
        exit (1);
    }

    elem_t* alpha = (elem_t*) calloc (rems_size+2, sizeof(elem_t));
    int* n = (int*) calloc (rems_size+2, sizeof(int));
    int idx = 2;
    alpha[0] = a->elems[a->lt]; n[0] = a->lt;
    alpha[1] = b->elems[b->lt]; n[1] = b->lt;
    for (int i = rems_size-1; i >= 0; i--) {
        if (lcrems[i]) {
            alpha[idx] = lcrems[i]; n[idx] = i;
            idx++;
        }
    }
    free (lcrems);
    // fprintf (stderr, "alpha and n size is %ld\n", rems_size+2);  // TEST
    rems_size = idx; // update rems_size

    // for (int i = 0; i < rems_size; i++) { // TEST
    //     fprintf (stderr, "n[%d] = %d\n", i, n[i]);
    // } // TEST

    if (n[idx-1]) {
        fprintf (stderr, "degree must be zero! (n[%ld-1]=%d)\n", rems_size, n[rems_size-1]);
        exit (1);
    }

    unsigned long tha = 0;
    elem_t res_coef = smallprimefield_exp (alpha[rems_size-1], n[rems_size-2], Pptr); // 
    // fprintf (stderr, "[init] res_coef = %lld \n", res_coef); // TEST
    for (int i = 1; i < rems_size-1; i++) {
        tha += n[i-1]*n[i];
        res_coef = smallprimefield_mul (res_coef, smallprimefield_exp (alpha[i], n[i-1]-n[i+1], Pptr), Pptr); // 
    }
    if (tha & 1) {
        res_coef = smallprimefield_mul (res_coef, smallprimefield_convert_in (-1, Pptr), Pptr);
    }
    free (alpha);
    free (n);
    *res = makePolynomial_spX (1);
    (*res)->elems[0] = res_coef;
    // fprintf (stderr, "[return] res_coef = %lld\n", res_coef); // TEST
}



void _premPolynomialsInFrom_spX_inp (elem_t* a, polysize_t* ad, elem_t* b, polysize_t bd, const Prime_ptr* Pptr)
{
    if (*ad == 0) {
        return;
    } else if (bd == 0) {
        fprintf (stderr, "divisor couldn't be zero in prem_spX\n");
        exit (1);
    }

    elem_t* ae = a;
    elem_t* be = b;
    polysize_t a_deg=*ad, b_deg = bd; 
    int delta = (a_deg - b_deg + 1) & 1;
    while (a_deg >= b_deg) {
        for (polysize_t i = a_deg-1, j = b_deg-1; i > -1; --i, --j) {
            if (j > -1) {
                ae[i] = smallprimefield_sub(smallprimefield_mul(ae[i], be[b_deg], Pptr), 
                                            smallprimefield_mul(be[j], ae[a_deg], Pptr), Pptr);
            } else {
                ae[i] = smallprimefield_mul (ae[i], be[b_deg], Pptr);
            }
        }
        --a_deg;
    }
    if (a_deg < 0) { a_deg = 0; }
    if (delta) {
        for (polysize_t i = 0; i < b_deg; i++) {
            if (ae[i]) {
                ae[i] = smallprimefield_sub (0,ae[i], Pptr);
                a_deg = i;
            }
        }
    } else {
        for (int i = b_deg-1; i > -1; --i) {
            if (ae[i]) {
                a_deg = i; 
                break;
            }
        }
    }
    a = ae;
    *ad = a_deg;
}

void brownSubResultantInForm_spX (duspoly_t* a, duspoly_t* b, duspolysA_t** subresA, polysize_t* sz, const Prime_ptr* Pptr)
{
    if (subresA == NULL){
        return;
    }

   	if (isZero_spX (a)) {
		if (isZero_spX (b) || b->lt) {
            *sz = 0;
            return;
		} else {
			*subresA = makePolysA_spX (1);
            (*subresA)->polys[0] = constPolynomial_spX (1, Pptr);
            *sz = 1;
            return;
		}
	} else if (isZero_spX (b)) {
		if (a->lt) {
            *sz = 0;
			return;
		} else {
			*subresA = makePolysA_spX (1);
            (*subresA)->polys[0] = constPolynomial_spX (1, Pptr);
            *sz = 1;
            return;
		}
	}

	if (!a->lt && !b->lt) {
        *subresA = makePolysA_spX (1);
        (*subresA)->polys[0] = constPolynomial_spX (1, Pptr);
        *sz = 1;
        return;
	}

    polysize_t ad=a->lt, bd=b->lt;
    if (ad < bd) {  
        brownSubResultantInForm_spX (b, a, subresA, sz, Pptr);
    }

    elem_t* ae = (elem_t*) malloc ((ad+1)*sizeof(elem_t));
    elem_t* be = (elem_t*) malloc ((bd+1)*sizeof(elem_t));
    memcpy (ae, a->elems, (ad+1)*sizeof(elem_t));
    memcpy (be, b->elems, (bd+1)*sizeof(elem_t));
    duspolysA_t* sres = makePolysA_spX (MIN_spX (ad, bd));
    *sz = MIN_spX (ad, bd);
    elem_t* aa = ae, *bb = be;
    elem_t s = smallprimefield_exp (be[bd], ad-bd, Pptr);
    _premPolynomialsInFrom_spX_inp (ae, &ad, be, bd, Pptr);
    elem_t* t = ae; // ae = prem (ae, bd)
    ae = be; 
    be = t; // swap(ae, be)
    int k = ad;
    ad = bd;
    bd = k; // swap (ad, bd)

    while (bd > 0 || be[0]) {
        int delta = ad-bd;
        if (ad > 0) {
            sres->polys[ad-1] = makePolynomial_spX (bd+1);
            sres->polys[ad-1]->lt = bd;
            memcpy (sres->polys[ad-1]->elems, be, (bd+1)*sizeof(elem_t));
        }
        if (delta > 1) { // lazard-opt in multivariate imp.
            sres->polys[bd] = makePolynomial_spX (bd+1);
            sres->polys[bd]->lt = bd;
            elem_t e = smallprimefield_exp (smallprimefield_div (be[bd], s, Pptr), delta-1, Pptr);
            for (polysize_t i = 0; i <= bd; i++) {
                sres->polys[bd]->elems[i] = smallprimefield_mul (be[i], e, Pptr);
            }
        }
        if (!bd) {break;} // terminate while-loop
        elem_t a = ae[ad];
        _premPolynomialsInFrom_spX_inp (ae, &ad, be, bd, Pptr);
        t = ae;
        ae = be;
        be = t; // swap
        k = ad;
        ad = bd;
        bd = k; // swap
        a = smallprimefield_inv (smallprimefield_mul (a, smallprimefield_exp(s, delta, Pptr), Pptr), Pptr);
        for (polysize_t i = 0; i <= bd; i++) {
            be[i] = smallprimefield_mul (be[i], a, Pptr);
        }
        if (delta > 1) { // double-checked later
            if (ae != NULL || ad > 0) {
                free (ae); ae = NULL;
            }
            if (sres->polys[ad] != NULL) {
                ae = (elem_t*) calloc (sres->polys[ad]->lt+1, sizeof(elem_t));
                memcpy (ae, sres->polys[ad]->elems, (sres->polys[ad]->lt+1)*sizeof(elem_t));
                ad = sres->polys[ad]->lt;
            }
        }
        s = ae[ad];
    }

    if (ae != NULL) {
        free (ae);
    }
    if (be != NULL) {
        free (be);
    }
    *subresA = sres;
}

void sylvSubResultantInForm_spX (duspoly_t* a, duspoly_t* b, duspolys_t** subres, polysize_t* sz, const Prime_ptr* Pptr)
{
    if (subres == NULL){
        return;
    }

    duspolys_t* tt = (duspolys_t*) malloc (sizeof(duspolys_t));

    if (isZero_spX (a)){
        tt->poly = NULL;
		tt->next = NULL;
		*subres = tt;
		*sz = 0;

		return;
    }
    
    if (isZero_spX (b)){
        tt->poly = deepCopyPolynomial_spX (a);
		tt->next = NULL;
		*subres = tt;
		*sz = 1;

		return;
    }

    duspoly_t* A;
    duspoly_t* B;

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_ab = deg_a * deg_b;

    polysize_t size = 0;
    elem_t mone = 0;
    int isSwap = 0;
        
    if (deg_a < deg_b) {
        if (deg_ab & 1) {
            mone = smallprimefield_convert_in (-1, Pptr);
        }
	
		A = b;
		B = a;

		deg_a = deg_b;
		deg_b = degPolynomial_spX (B);
	
		isSwap = 1;

		// fprintf (stderr, "DUSP Error, SubResultant Algorithm doesnt work if deg(f) < deg(g)!\n"); // TEST
		// exit (1);  // TEST
	
    } else {
		A = a;
		B = b;
    }

    duspolys_t* subrr = NULL;
    duspolys_t* tail = NULL;

    tt->poly = deepCopyPolynomial_spX (A);
    tt->next = NULL;
    
    subrr = tt;
    tail = tt;
    size++;

    tt = (duspolys_t*) malloc (sizeof (duspolys_t));
    tt->poly = deepCopyPolynomial_spX (B);
    tt->next = NULL;
    
    tail->next = tt;
    tail = tt;
    size++;

    duspoly_t* Sm = A;
    duspoly_t* h = B;
    duspoly_t* Sn = NULL;
    duspoly_t* Snt = NULL;

    elem_t lc_h;
    elem_t cnst;
    elem_t none;
    elem_t total = smallprimefield_convert_in (1, Pptr);
    polysize_t deg_h;

    if (h->lt == 0) {

		lc_h = leadingCoeffInForm_spX (h);	
		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm), Pptr);
	
		if (isSwap && mone) {
            cnst = smallprimefield_mul (cnst, mone, Pptr);
		}
	
		tt = (duspolys_t*) malloc (sizeof (duspolys_t));
		tt->poly = constPolynomialInForm_spX (cnst, Pptr);
		tt->next = NULL;
	
		tail->next = tt;
		tail = tt;
		size++;
	
		*subres = subrr;
		*sz = size;
	
		return;
    }

    while (1) {

		lc_h = leadingCoeffInForm_spX (h);
		deg_h = degPolynomial_spX (h);
	
		// plainRemPolynomialsInForm_spX (Sm, h, &Sn, Pptr); // TODO: fastDIv
		remPolynomialsInForm_spX (Sm, h, &Sn, Pptr);

		if (isZero_spX (Sn)) {	

			tt = (duspolys_t*) malloc (sizeof (duspolys_t));
			tt->poly = NULL;
			tt->next = NULL;
	    
			tail->next = tt;
			tail = tt;
			size++;
	    
			*subres = subrr;
			*sz = size;
    
			return;
		}

		none = smallprimefield_exp (smallprimefield_convert_in (-1, Pptr),
									deg_h * degPolynomial_spX (Sm), Pptr);

		if (isSwap && mone) {
			none = smallprimefield_mul (none, mone, Pptr);
		}

		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm) -
									degPolynomial_spX (Sn), Pptr);

		// fprintf (stderr, "[TEST] lc_h := %lld, deg(Sm) = %ld, deg(Sn) = %ld\n", lc_h, degPolynomial_spX (Sm), degPolynomial_spX (Sn));  // TEST 
		// fprintf (stderr, "[TEST] lc_h^(m-p) := %lld\n", cnst);  // TEST 

		cnst = smallprimefield_mul (cnst, none, Pptr);
		total = smallprimefield_mul (total, cnst, Pptr);
	
		scalarMulPolynomialInForm_spX (Sn, total, &Snt, Pptr);

		/* freePolynomial_spX (&Sn); Sn = NULL; */
	
		tt = (duspolys_t*) malloc (sizeof (duspolys_t));
		tt->poly = Snt;
		tt->next = NULL;
	
		tail->next = tt;
		tail = tt;
		size++;
	
		Sm = h;
		h = Sn;
		Snt = NULL;
		Sn = NULL;
	
		if (isZero_spX (h) || h->lt == 0) {
			*subres = subrr;
			*sz = size;
			return;
		}
    }
}


void sylvSubResultantInFormA_spX (duspoly_t* a, duspoly_t* b, duspolysA_t** subresA, polysize_t* sz, const Prime_ptr* Pptr)
{
    if (subresA == NULL){
        return;
    }

   	if (isZero_spX (a)) {
		if (isZero_spX (b) || b->lt) {
            *sz = 0;
            return;
		} else {
			*subresA = makePolysA_spX (1);
            (*subresA)->polys[0] = constPolynomial_spX (1, Pptr);
            *sz = 1;
            return;
		}
	} else if (isZero_spX (b)) {
		if (a->lt) {
            *sz = 0;
			return;
		} else {
			*subresA = makePolysA_spX (1);
            (*subresA)->polys[0] = constPolynomial_spX (1, Pptr);
            *sz = 1;
            return;
		}
	}

	if (!a->lt && !b->lt) {
        *subresA = makePolysA_spX (1);
        (*subresA)->polys[0] = constPolynomial_spX (1, Pptr);
        *sz = 1;
        return;
	}

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_ab = deg_a * deg_b;
    duspolysA_t* sres = makePolysA_spX (MIN_spX (deg_a, deg_b));

    duspoly_t* A;
    duspoly_t* B;

    polysize_t size = 0;
    elem_t mone = 0;
    int isSwap = 0;
        
    if (deg_a < deg_b) {
        if (deg_ab & 1) {
            mone = smallprimefield_convert_in (-1, Pptr); }
		A = b;
		B = a;
		deg_a = deg_b;
		deg_b = degPolynomial_spX (B);
		isSwap = 1;
    } else {
		A = a;
		B = b;
    }

    // sres->polys[size] = deepCopyPolynomial_spX (A);
    // size++;
    // sres->polys[size] = deepCopyPolynomial_spX (B);
    // size++;

    duspoly_t *Sm = A, *h = B, *Sn = NULL, *Snt = NULL;
    elem_t lc_h, cnst, none;
    elem_t total = smallprimefield_convert_in (1, Pptr);
    polysize_t deg_h;

    if (h->lt == 0) {
		lc_h = leadingCoeffInForm_spX (h);	
		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm), Pptr);
		if (isSwap && mone) {
            cnst = smallprimefield_mul (cnst, mone, Pptr);
		}
        sres->polys[size] = constPolynomialInForm_spX (cnst, Pptr);
		size++;
		*subresA = sres;
        *sz = size;
		return;
    }

    while (1) {

		lc_h = leadingCoeffInForm_spX (h);
		deg_h = degPolynomial_spX (h);
		// plainRemPolynomialsInForm_spX (Sm, h, &Sn, Pptr); // TODO: fastDIv
		remPolynomialsInForm_spX (Sm, h, &Sn, Pptr);

		if (isZero_spX (Sn)) {	
            sres->polys[size] = NULL;
			size++;	    
			*subresA = sres;
			*sz = size;
			return;
		}
		none = smallprimefield_exp (smallprimefield_convert_in (-1, Pptr), deg_h * degPolynomial_spX (Sm), Pptr);
		if (isSwap && mone) {
            none = smallprimefield_mul (none, mone, Pptr);
		}

		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm) -
									degPolynomial_spX (Sn), Pptr);

		// fprintf (stderr, "[TEST] lc_h := %lld, deg(Sm) = %ld, deg(Sn) = %ld\n", lc_h, degPolynomial_spX (Sm), degPolynomial_spX (Sn));  // TEST 
		// fprintf (stderr, "[TEST] lc_h^(m-p) := %lld\n", cnst);  // TEST 

		cnst = smallprimefield_mul (cnst, none, Pptr);
		total = smallprimefield_mul (total, cnst, Pptr);
	
		scalarMulPolynomialInForm_spX (Sn, total, &Snt, Pptr);

		/* freePolynomial_spX (&Sn); Sn = NULL; */
	    sres->polys[size] = Snt;
    	size++;
	
		Sm = h;
		h = Sn;
		Snt = NULL;
		Sn = NULL;
	
		if (isZero_spX (h) || h->lt == 0) {
			*subresA = sres;
			*sz = size;
			return;
		}
    }
}

void _sylvSubResultantInFormA_spX (duspoly_t* a, duspoly_t* b, duspolysA_t** subresA, const polysize_t first_idx, polysize_t* sz, const Prime_ptr* Pptr)
{
    if (subresA == NULL){
        return;
    }

   	if (isZero_spX (a)) {
		if (isZero_spX (b) || b->lt) {
            *sz = 0;
            return;
		} else {
            (*subresA)->polys[first_idx] = constPolynomial_spX (1, Pptr);
            *sz = 1;
            return;
		}
	} else if (isZero_spX (b)) {
		if (a->lt) {
            *sz = 0;
			return;
		} else {
            (*subresA)->polys[first_idx] = constPolynomial_spX (1, Pptr);
            *sz = 1;
            return;
		}
	}

	if (!a->lt && !b->lt) {
        (*subresA)->polys[first_idx] = constPolynomial_spX (1, Pptr);
        *sz = 1;
        return;
	}

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t deg_ab = deg_a * deg_b;
    // duspolysA_t** sres_polys = (*subresA)->polys;

    duspoly_t* A;
    duspoly_t* B;

    polysize_t size = 0;
    elem_t mone = 0;
    int isSwap = 0;
        
    if (deg_a < deg_b) {
        if (deg_ab & 1) {
            mone = smallprimefield_convert_in (-1, Pptr); 
        }
		A = b;
		B = a;
		deg_a = deg_b;
		deg_b = degPolynomial_spX (B);
		isSwap = 1;
    } else {
		A = a;
		B = b;
    }

    // sres->polys[size] = deepCopyPolynomial_spX (A);
    // size++;
    // sres->polys[size] = deepCopyPolynomial_spX (B);
    // size++;

    duspoly_t* Sm = A;
    duspoly_t* h = B;
    duspoly_t* Sn = NULL;
    duspoly_t* Snt = NULL;

    elem_t lc_h;
    elem_t cnst;
    elem_t none;
    elem_t total = smallprimefield_convert_in (1, Pptr);
    polysize_t deg_h;

    if (h->lt == 0) {

		lc_h = leadingCoeffInForm_spX (h);	
		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm), Pptr);
	
		if (isSwap && mone) {
            cnst = smallprimefield_mul (cnst, mone, Pptr);
		}

        (*subresA)->polys[first_idx+size] = constPolynomialInForm_spX (cnst, Pptr);
		size++;
        *sz = size;
		return;
    }

    while (1) {

		lc_h = leadingCoeffInForm_spX (h);
		deg_h = degPolynomial_spX (h);
		// plainRemPolynomialsInForm_spX (Sm, h, &Sn, Pptr); // TODO: fastDIv
		remPolynomialsInForm_spX (Sm, h, &Sn, Pptr); // TODO: fastDIv

		if (isZero_spX (Sn)) {	
			// tt = (duspolys_t*) malloc (sizeof (duspolys_t));
			// tt->poly = NULL;
			// tt->next = NULL;
	    
			// tail->next = tt;
			// tail = tt;
            (*subresA)->polys[first_idx+size] = NULL;
			size++;
	    	*sz = size;
			return;
		}

		none = smallprimefield_exp (smallprimefield_convert_in (-1, Pptr),
									deg_h * degPolynomial_spX (Sm), Pptr);

		if (isSwap && mone) {
            none = smallprimefield_mul (none, mone, Pptr);
		}

		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm) -
									degPolynomial_spX (Sn), Pptr);

		// fprintf (stderr, "[TEST] lc_h := %lld, deg(Sm) = %ld, deg(Sn) = %ld\n", lc_h, degPolynomial_spX (Sm), degPolynomial_spX (Sn));  // TEST 
		// fprintf (stderr, "[TEST] lc_h^(m-p) := %lld\n", cnst);  // TEST 

		cnst = smallprimefield_mul (cnst, none, Pptr);
		total = smallprimefield_mul (total, cnst, Pptr);
	
		scalarMulPolynomialInForm_spX (Sn, total, &Snt, Pptr);

		/* freePolynomial_spX (&Sn); Sn = NULL; */
	    (*subresA)->polys[first_idx + size] = Snt;
    	size++;
	
		Sm = h;
		h = Sn;
		Snt = NULL;
		Sn = NULL;
	
		if (isZero_spX (h) || h->lt == 0) {
			*sz = size;
			return;
		}
    }
}

void hgcdSubResultantInFormA_spX (duspoly_t* a, duspoly_t* b, polysize_t k, duspolysA_t** ksubresA, polysize_t* sz, const Prime_ptr* Pptr)
{
    if (ksubresA == NULL){
        return;
    }

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    polysize_t min_deg = MIN_spX (deg_a, deg_b);
    duspoly_t* tmp = NULL;
    if (k < 0 || k >= min_deg) {
        *ksubresA = makePolysA_spX (2);
        (*ksubresA)->polys[0] = (deg_a >= deg_b) ? deepCopyPolynomial_spX (a) : deepCopyPolynomial_spX (b);
        (*ksubresA)->polys[1] = (deg_a >= deg_b) ? deepCopyPolynomial_spX (b) : deepCopyPolynomial_spX (a);
        *sz = 2;
        return;
    } else if (!k && min_deg < 2) {
        *ksubresA = makePolysA_spX (2);
        (*ksubresA)->polys[0] = (deg_a >= deg_b) ? deepCopyPolynomial_spX (b) : deepCopyPolynomial_spX (a);
         hgcdResultantInForm_spX (a, b, &tmp, Pptr);
        (*ksubresA)->polys[1] = tmp;
        *sz = 2;
        return;
    }

    polysize_t deg_ab = deg_a * deg_b;
    duspolysA_t* sres = makePolysA_spX (2);
    duspoly_t *A, *B;

    polysize_t size = 0;
    elem_t mone = 0;
    int isSwap = 0;
        
    if (deg_a < deg_b) {
        if (deg_ab & 1) {
            mone = smallprimefield_convert_in (-1, Pptr); 
        }
		A = b;
		B = a;
		deg_a = deg_b;
		deg_b = degPolynomial_spX (B);
		isSwap = 1;
    } else {
		A = a;
		B = b;
    }

    elem_t* lcrems = NULL;
    polysize_t rems_size = MIN_spX (a->lt, b->lt);
    _HGCDInFormAtMaxK_spX (A, B, k, &(sres->polys[0]), &(sres->polys[1]), &lcrems, Pptr);

    if (isZero_spX ((sres->polys[0])) && isZero_spX ((sres->polys[1]))) {
        if (lcrems != NULL)
            free (lcrems);
        *sz = 0;
        return;
    }

    elem_t* alpha = (elem_t*) calloc (rems_size+2, sizeof(elem_t));
    int* n = (int*) calloc (rems_size+2, sizeof(int));
    int idx = 2;
    alpha[0] = A->elems[A->lt]; n[0] = A->lt;
    alpha[1] = B->elems[B->lt]; n[1] = B->lt;
    for (int i = rems_size-1; i >= 0; i--) {
        if (lcrems[i]) {
            alpha[idx] = lcrems[i]; n[idx] = i;
            idx++;
        }
    }
    free (lcrems);
    rems_size = idx; // update rems_size

    // for (int i = 0; i < rems_size; i++) { // TEST
    //     fprintf (stderr, "[original] alpha[%d] = %lld \n", i, alpha[i]);
    // } // TEST

    if (isZero_spX (sres->polys[1])) {
        alpha[rems_size] = 0; n[rems_size] = 0;
        rems_size++;
    }

    if (rems_size == 2) {
        *ksubresA = sres;
        *sz = 2;
        free (alpha);
        free (n);
        return;
    }

    unsigned long tha = 0, tha1 = 0;
    elem_t res_coef, res_coef_k1;
    elem_t m1 = smallprimefield_convert_in (-1, Pptr);
    if (rems_size < 5) {
        tha1 = (n[0]-n[2])*(n[1]-n[2]);
        res_coef_k1 = smallprimefield_exp (alpha[1], n[0]-n[2], Pptr);
        res_coef = res_coef_k1;
        if (tha1 & 1) {
            res_coef_k1 = smallprimefield_mul (res_coef_k1, m1, Pptr);
        } 
        if  (isSwap && mone) {
            res_coef_k1 = smallprimefield_mul (res_coef_k1, m1, Pptr);
        }
        if (rems_size == 3) {
            scalarMulPolynomialInForm_spX_inp (&(sres->polys[1]), res_coef_k1, Pptr);
            *ksubresA = sres;
            *sz = 2;
            free (alpha);
            free (n);
            return;
        } else {
            scalarMulPolynomialInForm_spX_inp (&(sres->polys[0]), res_coef_k1, Pptr);
        }

        tha = (n[0]-n[3])*(n[1]-n[3]) + (n[1]-n[3])*(n[2]-n[3]);
        // res_coef = smallprimefield_exp (alpha[1], n[0]-n[2], Pptr); // already computed for res_coef_k1
        res_coef = smallprimefield_mul (res_coef, smallprimefield_exp (alpha[2], n[1]-n[3], Pptr), Pptr);
        if (tha & 1) {
            res_coef = smallprimefield_mul (res_coef, m1, Pptr);
        } 
        if  (isSwap && mone) {
            res_coef = smallprimefield_mul (res_coef, m1, Pptr);
        }
        scalarMulPolynomialInForm_spX_inp (&(sres->polys[1]), res_coef, Pptr);
        *ksubresA = sres;
        *sz = 2;
        free (alpha);
        free (n);
        return;
    }

    // for (int i = 0; i < rems_size; i++) { // TEST
    //     fprintf (stderr, "n[%d] = %d ||| alpha[%d] = %lld \n", i, n[i], i, alpha[i]);
    // } // TEST

    // rems_size >= 5: 
    res_coef_k1 = smallprimefield_convert_in (1, Pptr);
    int goal_idx = rems_size - 2; // compute subres[k+1]
    for (polysize_t j = 1; j <= goal_idx-1; j++) {
        tha1 += (n[j-1]-n[goal_idx])*(n[j]-n[goal_idx]);
        tha += (n[j-1]-n[goal_idx+1])*(n[j]-n[goal_idx+1]);
        res_coef_k1 = smallprimefield_mul (res_coef_k1, smallprimefield_exp (alpha[j], n[j-1]-n[j+1], Pptr), Pptr);
    }

    tha += (n[rems_size-3]-n[goal_idx+1])*(n[rems_size-2]-n[goal_idx+1]);    
    res_coef = res_coef_k1; // goal_idx = rems_size - 1;
    res_coef = smallprimefield_mul (res_coef, smallprimefield_exp (alpha[rems_size-2], n[rems_size-3]-n[rems_size-1], Pptr), Pptr);

    if (tha & 1) {
        res_coef = smallprimefield_mul (res_coef, m1, Pptr);
    } 
    if (tha1 & 1) {
        res_coef_k1 = smallprimefield_mul (res_coef_k1, m1, Pptr);
    }
    if  (isSwap && mone) {
        res_coef = smallprimefield_mul (res_coef, m1, Pptr);
        res_coef_k1 = smallprimefield_mul (res_coef_k1, m1, Pptr);
    }
    scalarMulPolynomialInForm_spX_inp (&(sres->polys[0]), res_coef_k1, Pptr);
    scalarMulPolynomialInForm_spX_inp (&(sres->polys[1]), res_coef, Pptr);

    *ksubresA = sres;
    *sz = 2;
    free (alpha);
    free (n);
    return;
}

void YapHGCD_BaseCaseInFrom_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, const Prime_ptr* Pptr)
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
    
    if (deg_a > deg_b) {
		A = deepCopyPolynomial_spX (a);
		B = deepCopyPolynomial_spX (b);
    } else if (deg_a <= deg_b) {
		// TODO : Using Iterative HGCD
		/* fprintf (stderr, "DUSP Alert: In HalfGCD, deg(a) <= deg(b)!\n"); // TEST */
		*M = E;
		return;
    }
    
    polysize_t m = (deg_a + 1)>>1;

    /* fprintf (stderr, "\n**** YapHalfGCDMatrixInForm_spX (m = %lu) **** \n", m);  // TEST */
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
    
    YapHGCD_BaseCaseInFrom_spX (A0, B0, &RR, Quo, LC, DD, len, Pptr);  // first recursive call

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);
    
    /* fprintf (stderr, "    YapHalfGCDMatrixInForm_spX (A0, B0, &RR, Pptr);  // first recursive call\n");  // TEST */
    /* printMat4OutForm_spX (RR, Pptr); */

    mulMat4ToVecInForm_spX_inp (RR, &A, &B, Pptr); // A = Ap = r_{j-1} , B = Bp = r_{j}

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
		freePolynomial_spX (&B);
		
		freeMat4_spX (&E);
		
		*M = RR;
		return;
    }
	
	(*LC)[(*len)] = leadingCoeffInForm_spX (B);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (A) + degPolynomial_spX (B);

    duspoly_t* q = NULL;
    duspoly_t* r = NULL;
    duspoly_t* D = NULL;
	
    plainDivPolynomialsInForm_spX (A, B, &r, &q, Pptr); // TODO : Fast Monic Div

	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);	
	(*len) = (*len) + 1;

	
	/* fprintf (stderr, "A := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (A, Pptr), Pptr); */

	/* fprintf (stderr, "B := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (B, Pptr), Pptr); */

	/* fprintf (stderr, "q := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (q, Pptr), Pptr); */

	/* fprintf (stderr, "r := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (r, Pptr), Pptr); */

    /* duspolys_t* tt = (duspolys_t*) malloc (sizeof (duspolys_t)); */
    /* tt->poly = deepCopyPolynomial_spX (q); */
    /* tt->next = NULL; */
    
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
	
	
		// mulMat4ToMat4InForm_spX_inp_inp (RR, &TM, Pptr);
		mulMat4ToMat4InForm_spX_inp (TM, RR, &MM, Pptr);

		freeMat4_spX (&RR);
		freeMat4_spX (&E);
		freePolynomial_spX (&B);
	
		/* fprintf (stderr, "     RR*TM = MM = \n");  // TEST */
		/* printMat4OutForm_spX (MM, Pptr); */
	
		*M = MM;
		return;
    }

    /* elem_t lc = leadingCoeffInForm_spX (r); // IT SHOULD BE UNCOMMENT */

    /* elem_t lc = smallprimefield_convert_in (1, Pptr); // TEST */
    
    /* fprintf (stderr, "lc(r) = %lld\n", lc); // TEST */
    
    /* lc = smallprimefield_inv (lc, Pptr); // IT SHOULD BE UNCOMMENT */

    /* fprintf (stderr, "inv (lc(r)) = %lld\n", lc); // TEST */

    /* scalarMulPolynomialInForm_spX (r, lc, &D, Pptr); */

    /* fprintf (stderr, "    r*inv(lc) = ");  // TEST */
    /* printPolynomialOutForm_spX (D,  Pptr); */
    
    /* freePolynomial_spX (&r); */
	D = r; 

    duspoly_t* tq = NULL;
    
    // lc = smallprimefield_sub (0, lc, Pptr); // lc = - lc;
    /* scalarMulPolynomialInForm_spX (q, smallprimefield_sub (0, lc, Pptr), &tq, Pptr); */
	negPolynomialInForm_spX (q, Pptr);

					   
    /* fprintf (stderr, "    -q*inv(lc) = ");  // TEST */
    /* printPolynomialOutForm_spX (tq,  Pptr); */
        
    /* freePolynomial_spX (&q);     */
	tq = q;
	
    dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomial_spX (1, Pptr); // constPolynomialInForm_spX (lc, Pptr);
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

    YapHGCD_BaseCaseInFrom_spX (C0, D0, &S, Quo, LC, DD, len, Pptr); // second recursive call
	
    /* fprintf (stderr, "     S = \n");  // TEST */
    /* printMat4OutForm_spX (S, Pptr); */
    
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);
    
    /* mulMat4ToMat4InForm_spX_inp_inp (RR, &Q, Pptr); */
    /* mulMat4ToMat4InForm_spX_inp_inp (Q, &S, Pptr); */
	
    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;
    mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr);

    /* fprintf (stderr, "     Q*RR = \n");  // TEST */
    /* printMat4OutForm_spX (QRR, Pptr); */

    mulMat4ToMat4InForm_spX_inp (S, QRR, &SQRR, Pptr);
    
    freeMat4_spX (&RR);
    freeMat4_spX (&QRR);
    
    *M = SQRR;
    
    /* fprintf (stderr, "     S*Q*RR = \n");  // TEST */
    /* printMat4OutForm_spX (SQRR, Pptr); */
    
    /* *M = S; */
    return;
}


void YapHGCD_TotalCaseInFrom_spX (duspoly_t* a, duspoly_t* b, dusmat4_t** M,  duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, const Prime_ptr* Pptr)
{

    /* fprintf (stderr, "\n**** GCDMatrixInForm_wHGCD_spX **** \n");  // TEST */
    /* fprintf (stderr, "a = ");  // TEST */
    /* printPolynomial_spX (a, Pptr); */

    /* fprintf (stderr, "b = ");  // TEST */
    /* printPolynomial_spX (b, Pptr); */
    
    dusmat4_t* RR = NULL;
    YapHGCD_BaseCaseInFrom_spX (a, b, &RR, Quo, LC, DD, len, Pptr);
    
    // TODO : implement non in-place mulMat4ToVecInForm 
    duspoly_t* a0 = deepCopyPolynomial_spX (a);
    duspoly_t* b0 = deepCopyPolynomial_spX (b);
    
    mulMat4ToVecInForm_spX_inp (RR, &a0, &b0, Pptr);

    /* fprintf (stderr, "RR = \n"); // TEST */
    /* printMat4_spX (RR, Pptr); */

    /* fprintf (stderr, "a0 = ");  // TEST */
    /* printPolynomial_spX (a0, Pptr); */

    /* fprintf (stderr, "b0 = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    if (isZero_spX (b0)) {
		*M = RR;
		return;
    }

	(*LC)[(*len)] = leadingCoeffInForm_spX (b0);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (a0) + degPolynomial_spX (b0);
	
    duspoly_t* r   = NULL;
    duspoly_t* q   = NULL;
    dusmat4_t* Q   = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    
    plainDivPolynomialsInForm_spX (a0, b0, &r, &q, Pptr); // TODO: Fast Monic Divid

	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);
	(*len) = (*len) + 1;	

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
	
		// mulMat4ToMat4InForm_spX_inp_inp (Q, &RR, Pptr);
		mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr); 
	
		freeMat4_spX (&RR);
	
		*M = QRR;
		return;
    }
    
    /* elem_t lc = leadingCoeffInForm_spX (r); */
    /* lc = smallprimefield_inv (lc, Pptr); */
    
    /* scalarMulPolynomialInForm_spX (r, lc, &lcr, Pptr); */
    /* scalarMulPolynomialInForm_spX (q, lc, &lcq, Pptr); */

    /* freePolynomial_spX (&r); */
    /* freePolynomial_spX (&q); */
    
    /* Q->polys[0] = NULL; */
    /* Q->polys[1] = constPolynomial_spX (1, Pptr); */
    /* Q->polys[2] = constPolynomialInForm_spX (lc, Pptr); */
    /* Q->polys[3] = lcq; */
		
	Q->polys[0] = NULL;
	Q->polys[1] = constPolynomial_spX (1, Pptr);
	Q->polys[2] = constPolynomial_spX (1, Pptr);
	Q->polys[3] = q; // q = -q

    /* fprintf (stderr, "Q with lcq = \n"); // TEST */
    /* printMat4_spX (Q, Pptr); */
    
    /* fprintf (stderr, "b0 = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    /* fprintf (stderr, "lcr = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    
    dusmat4_t* S = NULL;
    YapHGCD_TotalCaseInFrom_spX (b0, r, &S, Quo, LC, DD, len, Pptr); 

    /* fprintf (stderr, "S = \n"); // TEST */
    /* printMat4_spX (S, Pptr); */

    freePolynomial_spX (&b0);
    freePolynomial_spX (&r);

    dusmat4_t* RRQ = NULL;
    dusmat4_t* SRRQ = NULL;
    
    mulMat4ToMat4InForm_spX_inp (Q, RR, &RRQ, Pptr);
    mulMat4ToMat4InForm_spX_inp (S, RRQ, &SRRQ, Pptr);

    /* fprintf (stderr, "M = \n"); // TEST */
    /* printMat4_spX (SQRR, Pptr); */

    freeMat4_spX (&RR);
    freeMat4_spX (&RRQ);
    
    *M = SRRQ;     
}



void extYapHGCD_BaseCaseInFrom_spX (duspoly_t* a, duspoly_t* b, polysize_t k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, const Prime_ptr* Pptr)
{
	
    dusmat4_t* E = NULL;
    identityMat4InForm_spX (&E, Pptr);
	
    if (isZero_spX (b) || isZero_spX (a)) {
		*M = E;
		return;
    }

	if ( (*len) > 0 && (*DD)[(*len)-1] < k) {
		*M = E;
		return;
	}
	
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
    duspoly_t* A = NULL;
    duspoly_t* B = NULL;
    
    if (deg_a > deg_b) {
		A = deepCopyPolynomial_spX (a);
		B = deepCopyPolynomial_spX (b);
    } else if (deg_a <= deg_b) {
		// TODO : Using Iterative HGCD
		/* fprintf (stderr, "DUSP Alert: In HalfGCD, deg(a) <= deg(b)!\n"); // TEST */
		*M = E;
		return;
    }
    
    polysize_t m = (deg_a + 1)>>1;

    /* fprintf (stderr, "\n**** YapHalfGCDMatrixInForm_spX (m = %lu) **** \n", m);  // TEST */
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
    
    extYapHGCD_BaseCaseInFrom_spX (A0, B0, k, &RR, Quo, LC, DD, len, Pptr);  // first recursive call

    freePolynomial_spX (&A0);
    freePolynomial_spX (&B0);
    
    /* fprintf (stderr, "    YapHalfGCDMatrixInForm_spX (A0, B0, &RR, Pptr);  // first recursive call\n");  // TEST */
    /* printMat4OutForm_spX (RR, Pptr); */

    mulMat4ToVecInForm_spX_inp (RR, &A, &B, Pptr); // A = Ap = r_{j-1} , B = Bp = r_{j}

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
		freePolynomial_spX (&B);
		
		freeMat4_spX (&E);
		
		*M = RR;
		return;
    }
	
	(*LC)[(*len)] = leadingCoeffInForm_spX (B);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (A) + degPolynomial_spX (B);

	if ((*DD)[(*len)] < k) {
		freePolynomial_spX (&A);
		freePolynomial_spX (&B);

		freeMat4_spX (&E);

		(*len) = (*len) + 1;
		*M = RR;
		
		return;
	}

    duspoly_t* q = NULL;
    duspoly_t* r = NULL;
    duspoly_t* D = NULL;
	
    plainDivPolynomialsInForm_spX (A, B, &r, &q, Pptr); // TODO : Fast Monic Div

	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);	
	(*len) = (*len) + 1;
	
	/* fprintf (stderr, "A := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (A, Pptr), Pptr); */

	/* fprintf (stderr, "B := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (B, Pptr), Pptr); */

	/* fprintf (stderr, "q := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (q, Pptr), Pptr); */

	/* fprintf (stderr, "r := " ); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (r, Pptr), Pptr); */

    /* duspolys_t* tt = (duspolys_t*) malloc (sizeof (duspolys_t)); */
    /* tt->poly = deepCopyPolynomial_spX (q); */
    /* tt->next = NULL; */
    
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
		
		// mulMat4ToMat4InForm_spX_inp_inp (RR, &TM, Pptr);
		mulMat4ToMat4InForm_spX_inp (TM, RR, &MM, Pptr);
		
		freeMat4_spX (&RR);
		freeMat4_spX (&E);
		freePolynomial_spX (&B);
		
		/* fprintf (stderr, "     RR*TM = MM = \n");  // TEST */
		/* printMat4OutForm_spX (MM, Pptr); */
		
		*M = MM;
		return;
    }

    /* elem_t lc = leadingCoeffInForm_spX (r); // IT SHOULD BE UNCOMMENT */

    /* elem_t lc = smallprimefield_convert_in (1, Pptr); // TEST */
    
    /* fprintf (stderr, "lc(r) = %lld\n", lc); // TEST */
    
    /* lc = smallprimefield_inv (lc, Pptr); // IT SHOULD BE UNCOMMENT */

    /* fprintf (stderr, "inv (lc(r)) = %lld\n", lc); // TEST */

    /* scalarMulPolynomialInForm_spX (r, lc, &D, Pptr); */

    /* fprintf (stderr, "    r*inv(lc) = ");  // TEST */
    /* printPolynomialOutForm_spX (D,  Pptr); */
    
    /* freePolynomial_spX (&r); */
	D = r; 

    duspoly_t* tq = NULL;
    
    // lc = smallprimefield_sub (0, lc, Pptr); // lc = - lc;
    /* scalarMulPolynomialInForm_spX (q, smallprimefield_sub (0, lc, Pptr), &tq, Pptr); */
	negPolynomialInForm_spX (q, Pptr);

					   
    /* fprintf (stderr, "    -q*inv(lc) = ");  // TEST */
    /* printPolynomialOutForm_spX (tq,  Pptr); */
        
    /* freePolynomial_spX (&q);     */
	tq = q;
	
    dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomial_spX (1, Pptr); // constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = tq;
    
    /* fprintf (stderr, "     Q = \n");  // TEST */
    /* printMat4OutForm_spX (Q, Pptr); */
    
    polysize_t kk = 2*m - degPolynomial_spX (B);
    
    /* fprintf (stderr, "kk = %lu\n", kk); // TEST */
    
    duspoly_t* C0 = NULL;
    duspoly_t* D0 = NULL;
    dusmat4_t* S  = NULL;
    
    rightShiftPolynomial_spX (B, &C0, kk);
    rightShiftPolynomial_spX (D , &D0, kk);

    freePolynomial_spX (&B);
    freePolynomial_spX (&D);

    extYapHGCD_BaseCaseInFrom_spX (C0, D0, kk, &S, Quo, LC, DD, len, Pptr); // second recursive call
	
    /* fprintf (stderr, "     S = \n");  // TEST */
    /* printMat4OutForm_spX (S, Pptr); */
    
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);
    
    /* mulMat4ToMat4InForm_spX_inp_inp (RR, &Q, Pptr); */
    /* mulMat4ToMat4InForm_spX_inp_inp (Q, &S, Pptr); */
	
    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;
    mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr);

    /* fprintf (stderr, "     Q*RR = \n");  // TEST */
    /* printMat4OutForm_spX (QRR, Pptr); */

    mulMat4ToMat4InForm_spX_inp (S, QRR, &SQRR, Pptr);
    
    freeMat4_spX (&RR);
    freeMat4_spX (&QRR);
    
    *M = SQRR;
    
    /* fprintf (stderr, "     S*Q*RR = \n");  // TEST */
    /* printMat4OutForm_spX (SQRR, Pptr); */
    
    /* *M = S; */
    return;
}


void extYapHGCD_TotalCaseInFrom_spX (duspoly_t* a, duspoly_t* b, polysize_t k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, const Prime_ptr* Pptr)
{

    /* fprintf (stderr, "\n**** GCDMatrixInForm_wHGCD_spX **** \n");  // TEST */
    /* fprintf (stderr, "a = ");  // TEST */
    /* printPolynomial_spX (a, Pptr); */

    /* fprintf (stderr, "b = ");  // TEST */
    /* printPolynomial_spX (b, Pptr); */

	if ( (*len) > 0 && (*DD)[(*len)-1] < k) {
		
		dusmat4_t* E = NULL;
		identityMat4InForm_spX (&E, Pptr);
		
		*M = E;
		return;
	}

    dusmat4_t* RR = NULL;
    extYapHGCD_BaseCaseInFrom_spX (a, b, k, &RR, Quo, LC, DD, len, Pptr);
    
    // TODO : implement non in-place mulMat4ToVecInForm 
    duspoly_t* a0 = deepCopyPolynomial_spX (a);
    duspoly_t* b0 = deepCopyPolynomial_spX (b);
    
    mulMat4ToVecInForm_spX_inp (RR, &a0, &b0, Pptr);

    /* fprintf (stderr, "RR = \n"); // TEST */
    /* printMat4_spX (RR, Pptr); */

    /* fprintf (stderr, "a0 = ");  // TEST */
    /* printPolynomial_spX (a0, Pptr); */

    /* fprintf (stderr, "b0 = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    if (isZero_spX (b0)) {
		*M = RR;

		return;
    }

	(*LC)[(*len)] = leadingCoeffInForm_spX (b0);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (a0) + degPolynomial_spX (b0);

	if ((*DD)[(*len)] < k) {		
		(*len) = (*len) + 1;
		*M = RR;

		return;
	}

    duspoly_t* r   = NULL;
    duspoly_t* q   = NULL;
    dusmat4_t* Q   = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    
    plainDivPolynomialsInForm_spX (a0, b0, &r, &q, Pptr); // TODO: Fast Monic Divid

	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);	
	(*len) = (*len) + 1;	
	
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
	
		// mulMat4ToMat4InForm_spX_inp_inp (Q, &RR, Pptr);
		mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr); 
	
		freeMat4_spX (&RR);
	
		*M = QRR;
		return;
    }
    
    /* elem_t lc = leadingCoeffInForm_spX (r); */
    /* lc = smallprimefield_inv (lc, Pptr); */
    
    /* scalarMulPolynomialInForm_spX (r, lc, &lcr, Pptr); */
    /* scalarMulPolynomialInForm_spX (q, lc, &lcq, Pptr); */

    /* freePolynomial_spX (&r); */
    /* freePolynomial_spX (&q); */
    
    /* Q->polys[0] = NULL; */
    /* Q->polys[1] = constPolynomial_spX (1, Pptr); */
    /* Q->polys[2] = constPolynomialInForm_spX (lc, Pptr); */
    /* Q->polys[3] = lcq; */

		
	Q->polys[0] = NULL;
	Q->polys[1] = constPolynomial_spX (1, Pptr);
	Q->polys[2] = constPolynomial_spX (1, Pptr);
	Q->polys[3] = q; // q = -q

    /* fprintf (stderr, "Q with lcq = \n"); // TEST */
    /* printMat4_spX (Q, Pptr); */
    
    /* fprintf (stderr, "b0 = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    /* fprintf (stderr, "lcr = ");  // TEST */
    /* printPolynomial_spX (b0, Pptr); */
    
    
    dusmat4_t* S = NULL;
    extYapHGCD_TotalCaseInFrom_spX (b0, r, k, &S, Quo, LC, DD, len, Pptr); 

    /* fprintf (stderr, "S = \n"); // TEST */
    /* printMat4_spX (S, Pptr); */

    freePolynomial_spX (&b0);
    freePolynomial_spX (&r);

    dusmat4_t* RRQ = NULL;
    dusmat4_t* SRRQ = NULL;
    
    mulMat4ToMat4InForm_spX_inp (Q, RR, &RRQ, Pptr);
    mulMat4ToMat4InForm_spX_inp (S, RRQ, &SRRQ, Pptr);

    /* fprintf (stderr, "M = \n"); // TEST */
    /* printMat4_spX (SQRR, Pptr); */

    freeMat4_spX (&RR);
    freeMat4_spX (&RRQ);
    
    *M = SRRQ;     
}


void resultantYapHGCDInForm_spX (duspoly_t* a, duspoly_t* b, duspoly_t** res, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
		*res = NULL;
		return;
    }

    if (isZero_spX (b)) {
		*res = deepCopyPolynomial_spX (a);
		return;
    }

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    elem_t res_t = smallprimefield_convert_in (1, Pptr);
    elem_t mone;

	duspoly_t* tmp = NULL;
	
	// TODO: uncommented after testing : 
	/* if (deg_b < PLAINRES_CROSSOVER || deg_a <PLAINRES_CROSSOVER) { */
	/* 	sylvResultantInForm_spX (a, b, res, Pptr); */
	/* 	return; */
	/* } */
	
    if (deg_b == 0) {
        mone = smallprimefield_exp (leadingCoeffInForm_spX (b), deg_a, Pptr);
		*res = constPolynomialInForm_spX (mone, Pptr);
		return;
    }
    
    if (deg_a < deg_b) {

		// TODO: OutForm ...
		mone = smallprimefield_convert_in (-1, Pptr);
		mone = smallprimefield_exp (mone, deg_a*deg_b, Pptr);
		
		sylvResultantInForm_spX (b, a, &tmp, Pptr);
		
		scalarMulPolynomialInForm_spX (tmp, mone, res, Pptr);
		
		freePolynomial_spX (&tmp);
	
		return;
    }

    duspoly_t* r = NULL;
	duspoly_t* A = NULL;
	duspoly_t* B = deepCopyPolynomial_spX (b); 

	
	elem_t coef;
	
	if (deg_a == deg_b) {
		
		plainRemPolynomialsInForm_spX (a, b, &r, Pptr); // ُTODO: fastDiv
		
		if (isZero_spX (r)) {
			*res = NULL;
			return;
		}

		A = B;
		B = r;

		// TODO: OutForm ...  
		mone = smallprimefield_convert_in (-1, Pptr);
		mone = smallprimefield_exp (mone, deg_a*deg_b, Pptr);
		
		res_t = smallprimefield_exp (leadingCoeffInForm_spX (b),
									 deg_a - degPolynomial_spX (r), Pptr);
		
		res_t = smallprimefield_mul (res_t, mone, Pptr);
		
	}
	
	if (A == NULL) {
		A = deepCopyPolynomial_spX (a);
	}

	duspoly_t** Quo = (duspoly_t**) malloc (sizeof (duspoly_t*)*(degPolynomial_spX (B) + 2));
	elem_t* LC = (elem_t*) calloc (degPolynomial_spX (B) + 2, sizeof (elem_t));
	polysize_t* DD = (polysize_t*) calloc (degPolynomial_spX (B) + 2, sizeof (polysize_t));
	polysize_t len = 0;
	
	LC[len] = leadingCoeffInForm_spX (A);
	DD[len] = degPolynomial_spX (A);
	len++;

	dusmat4_t* MM = NULL;
	
	YapHGCD_TotalCaseInFrom_spX	(A, B, &MM, &Quo, &LC, &DD, &len, Pptr);
    mulMat4ToVecInForm_spX_inp (MM, &A, &B, Pptr);

	/* for (polysize_t i = 0; i < len; i++) { */
	/* 	freePolynomial_spX (&Rem[i]); */
	/* } */

	/* free (Rem); */
	
	/* // TEST */
	/* for (int i = 0; i < len; i++) { */
	/* 	fprintf (stderr, "LC[%d] = %lld \t DD[%d] = %ld\n", i, smallprimefield_convert_out (LC[i], Pptr), i, DD[i]); */

	/* } */
	/* fprintf (stderr, "\n"); */

	
	/* fprintf (stderr, "A_out :="); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (A, Pptr), Pptr); */

	/* fprintf (stderr, "B_out :="); */
	/* printPolynomial_spX (convertPolynomialFromMontgomery_spX (B, Pptr), Pptr); */
	
	freeMat4_spX (&MM);

	// TODO: we probably dont need this part because of the HGCD structure!
	if (!isZero_spX (B)) {

		LC[len] = leadingCoeffInForm_spX (B);
		DD[len] = degPolynomial_spX (B);
		len++;
		
		r = NULL;
		plainRemPolynomialsInForm_spX (A, B, &r, Pptr); // ُTODO: fastDiv
		
		freePolynomial_spX (&A);
		
		A = B;
		B = r;		
	}

	if (degPolynomial_spX (A) > 0 && isZero_spX (B)) {
		*res = NULL;

		freePolynomial_spX (&A);
		
		return;
	}

	for (polysize_t i = 0; i < len-2; i++) {
		
		// TODO: OutForm ...  
		if (DD[i] & DD[i+1] & 1) {
			mone = smallprimefield_convert_in (-1, Pptr);
		} else {
			mone = smallprimefield_convert_in (1, Pptr);
		}
		 
		coef = smallprimefield_exp (LC[i+1], DD[i]-DD[i+2], Pptr);
		coef =  smallprimefield_mul (coef, mone, Pptr);
		res_t = smallprimefield_mul (res_t, coef, Pptr);
		
	}

	if (degPolynomial_spX (A) == 0) {
		
		coef = smallprimefield_exp (LC[len-1], DD[len-2], Pptr);
		res_t = smallprimefield_mul (res_t, coef, Pptr);

		*res = constPolynomialInForm_spX (res_t, Pptr);
		
		return;
	} 

	// TODO: OutForm ...  
	if ( DD[len-2] & DD[len-1] & 1) {
		mone = smallprimefield_convert_in (-1, Pptr);
	} else {
		mone = smallprimefield_convert_in (1, Pptr);
	}
	
	coef = smallprimefield_exp (LC[len-1], DD[len-2] - degPolynomial_spX (B), Pptr);
	coef = smallprimefield_mul (coef, mone, Pptr);
	res_t = smallprimefield_mul (res_t, coef, Pptr);
	
    r = constPolynomialInForm_spX (res_t, Pptr);

	sylvResultantInForm_spX (A, B, &tmp, Pptr);

	if (!isZero_spX (tmp)) {		
		plainMulPolynomialsInForm_spX (r, tmp, res, Pptr);

		freePolynomial_spX (&tmp);

	} else {
		*res = r;
	}
	 
	freePolynomial_spX (&A);
	freePolynomial_spX (&B);
	freePolynomial_spX (&r);

	free (LC);
	free (DD);

}


void quotientsToRemainders_spX (duspoly_t** Quo, duspoly_t*** Rem, polysize_t qlen, const Prime_ptr* Pptr)
{
	if (qlen < 0) {
		fprintf (stderr, "DUSP Error: quotient len is out of range!\n");
		exit (1);
	}

	duspoly_t** Rem_t = (duspoly_t**) calloc (qlen+1, sizeof (duspoly_t*));

	if (qlen == 0) {
		*Rem = Rem_t;
		return;
	}

	if (!isZero_spX (Quo[0])) {
		fprintf (stderr, "DUSP Error: Quo[0] is not zero!\n");
		exit (1);
	}

	if (isZero_spX ((*Rem)[0]) || isZero_spX((*Rem)[1])) {
		fprintf (stderr, "DUSP Error: Rem[0] and Rem[1] have to non-zero!\n");
		exit (1);
	}

	Rem_t[0] = (*Rem)[0];
	Rem_t[1] = (*Rem)[1];

	duspoly_t* rq = NULL;
	
	for (polysize_t i = 1; i < qlen ; i++) {

		if ( Quo[i] == NULL ) {
			Rem_t[i+1] = Rem_t[i-1];
			continue;
		}
		
		// rq = Rem_t[i]*Quo[i]
		plainMulPolynomialsInForm_spX (Rem_t[i], Quo[i], &rq, Pptr);
		
		// Rem_t[i+1] = Rem_t[i-1] - rq
		subPolynomialsInForm_spX (Rem_t[i-1], rq, &Rem_t[i+1], Pptr);

		freePolynomial_spX (&rq);
		rq = NULL;
	}

	*Rem = Rem_t;

}

void subResultantYapHGCDInForm_spX (duspoly_t* a, duspoly_t* b, duspolys_t** subres, polysize_t* sz, const Prime_ptr* Pptr)
{

    duspolys_t* tt = (duspolys_t*) malloc (sizeof(duspolys_t));

    if (isZero_spX (a)){
        tt->poly = NULL;
		tt->next = NULL;
		*subres = tt;
		*sz = 0;
		
		return;
    }
    
    if (isZero_spX (b)){
        tt->poly = deepCopyPolynomial_spX (a);
		tt->next = NULL;
		*subres = tt;
		*sz = 1;
		
		return;
    }
	
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    /* elem_t res_t = smallprimefield_convert_in (1, Pptr); */
    	
	duspoly_t* A;
    duspoly_t* B;
	
    polysize_t size = 0;
    elem_t mone = 0;
    int isSwap = 0;
    
    if (deg_a < deg_b) {
		mone = smallprimefield_convert_in (-1, Pptr);
		mone = smallprimefield_exp (mone, deg_a*deg_b, Pptr);
	
		A = b;
		B = a;

		deg_a = deg_b;
		deg_b = degPolynomial_spX (B);
	
		isSwap = 1;

    } else {
		A = a;
		B = b;
    }

	duspolys_t* subrr = NULL;
    duspolys_t* tail = NULL;

    tt->poly = deepCopyPolynomial_spX (A);
    tt->next = NULL;
    
    subrr = tt;
    tail = tt;
    size++;

    tt = (duspolys_t*) malloc (sizeof (duspolys_t));
    tt->poly = deepCopyPolynomial_spX (B);
    tt->next = NULL;
    
    tail->next = tt;
    tail = tt;
    size++;

	/* duspoly_t** subR = (duspoly_t**) malloc (sizeof(duspoly_t*)*(deg_b+1)); */
	/* polysize_t subRsz = 2; */
	/* subR[0] = deepCopyPolynomial_spX (a); */
	/* subR[1] = deepCopyPolynomial_spX (b); */

	duspoly_t* Sm = A;
    duspoly_t* h = B;
    duspoly_t* Sn = NULL;
    duspoly_t* Snt = NULL;

    elem_t lc_h;
    elem_t cnst;
    elem_t none;
    elem_t total = smallprimefield_convert_in (1, Pptr);
    polysize_t deg_h;

	if (deg_a == deg_b) {

		h = NULL;
		plainRemPolynomials_spX (A, B, &h, Pptr);
		
		tt = (duspolys_t*) malloc (sizeof (duspolys_t));
		tt->poly = deepCopyPolynomial_spX (h);
		tt->next = NULL;
	
		tail->next = tt;
		tail = tt;
		size++;

		if (isZero_spX (h)) {
			*subres = subrr;
			*sz = size;
			return;
		}

		fprintf (stderr, "DUSP Support: Subresultant algorithm doesnt support polynomials with the same degrees!\n");
		exit (1);
		
		/* Sm = B; */
		/* A = B; */
		/* B = h; */
		
		/* deg_a = degPolynomial_spX (A); */
		/* deg_b = degPolynomial_spX (B); */
		
	}

	if (h->lt == 0) {

				
		lc_h = leadingCoeffInForm_spX (h);	
		cnst = smallprimefield_exp (lc_h, deg_a, Pptr);
		
		if (isSwap) {
			cnst = smallprimefield_mul (cnst, mone, Pptr);
		}

		tt = (duspolys_t*) malloc (sizeof (duspolys_t));
		tt->poly = constPolynomialInForm_spX (cnst, Pptr);
		tt->next = NULL;
	
		tail->next = tt;
		tail = tt;
		size++;
	
		*subres = subrr;
		*sz = size;
	
		return;
    }

	duspoly_t** Quo = (duspoly_t**) calloc ((deg_b + 2), sizeof (duspoly_t*));
	elem_t* LC = (elem_t*) calloc (deg_b + 2, sizeof (elem_t));
	polysize_t* DD = (polysize_t*) calloc (deg_b + 2, sizeof (polysize_t));
	polysize_t len = 0;

	Quo[len] = NULL;
	LC[len] = leadingCoeffInForm_spX (a);
	DD[len] = deg_a;
	len++;
	
	dusmat4_t* MM = NULL;
	
	YapHGCD_TotalCaseInFrom_spX	(A, B, &MM, &Quo, &LC, &DD, &len, Pptr);
	
	duspoly_t** Rem_t = (duspoly_t**) calloc (2, sizeof (duspoly_t*));
	Rem_t[0] = A;
	Rem_t[1] = B;
	
	quotientsToRemainders_spX (Quo, &Rem_t, len, Pptr);

	/* // TEST */
	/* for (int i = 0; i < len+1; i++) { */
	/* 	fprintf (stderr, "Rem_t[%d] = ", i); */
	/* 	printPolynomial_spX (convertPolynomialFromMontgomery_spX (Rem_t[i], Pptr), Pptr); */
	/* } */
	
	/* // TEST */
	/* for (int i = 0; i < len; i++) { */
	/* 	fprintf (stderr, "LC[%d] = %lld \t DD[%d] = %ld\n", i, smallprimefield_convert_out (LC[i], Pptr), i, DD[i]); */
	/* } */
	
	freeMat4_spX (&MM);
	
	for (int i = 0; i < len; i++) {
		if (Quo[i] != NULL) {
			freePolynomial_spX (&Quo[i]);
		}
	}
	free (Quo);
	
	/* res_t = smallprimefield_exp (LC[1], DD[0] - DD[1], Pptr); */
	
	for (polysize_t i = 1; i < len; i++) {
		
		/* // TODO: OutForm ... */
		/* /\* if ((DD[i]-DD[i+1])*(DD[0]-DD[i+1]+i+1) % 2) { *\/ */
		/* if ( DD[i] & DD[i+1] & 1) { */
		/* 	mone = smallprimefield_convert_in (-1, Pptr); */
		/* } else { */
		/* 	mone = smallprimefield_convert_in (1, Pptr); */
		/* } */
		
		/* coef = smallprimefield_mul (LC[i], LC[i+1], Pptr); */
		/* coef = smallprimefield_exp (coef, DD[i]-DD[i+1], Pptr); */
		/* coef = smallprimefield_mul (coef, mone, Pptr); */
		
		/* res_t = smallprimefield_mul (res_t, coef, Pptr); // S[DD[i+1]] */

		/* /\* subR[subRsz] = constPolynomialInForm_spX (res_t, Pptr); *\/ */
		/* /\* subRsz++; *\/ */
		
		/* fprintf (stderr, "S[%ld] = %lld; \n", DD[i+1], smallprimefield_convert_out (res_t, Pptr)); */
		
		lc_h = leadingCoeffInForm_spX (h);
		deg_h = degPolynomial_spX (h);

		Sn = Rem_t[i+1];
		none = smallprimefield_exp (smallprimefield_convert_in (-1, Pptr), // TODO
									deg_h * degPolynomial_spX (Sm), Pptr);

		if (isSwap) {
			none = smallprimefield_mul (none, mone, Pptr);
		}

		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm) -
									degPolynomial_spX (Sn), Pptr);

		cnst = smallprimefield_mul (cnst, none, Pptr);
		total = smallprimefield_mul (total, cnst, Pptr);
		
		scalarMulPolynomialInForm_spX (Sn, total, &Snt, Pptr);

		tt = (duspolys_t*) malloc (sizeof (duspolys_t));
		tt->poly = Snt;
		tt->next = NULL;
		
		tail->next = tt;
		tail = tt;
		size++;

		Sm = h;
		h = Sn;
		Snt = NULL;
		Sn = NULL;
		
		if (isZero_spX (h) || h->lt == 0) {
			break;
		}
		
	}

	*subres = subrr;
	*sz = size;

	for (int i = 2; i < len; i++) {
		if (Rem_t[i] != NULL) {
			freePolynomial_spX (&Rem_t[i]);
		}
	}
	free (Rem_t);

	free (LC);
	free (DD);

}

void kSubResultantYapHGCDInForm_spX (duspoly_t* a, duspoly_t* b, polysize_t k, duspoly_t** ksubres, const Prime_ptr* Pptr)
{

	if (k < 0) {
		fprintf (stderr, "DUSP Error: k is out of renge in subResutlant algorithm!\n");
		exit (1);
	}

	if (k == 0) {
		resultantYapHGCDInForm_spX (a, b, ksubres, Pptr);
		return;
	}

	if (k == degPolynomial_spX (b)) {
		*ksubres = deepCopyPolynomial_spX (b);
		return;
	}

	if (k == degPolynomial_spX (a)) {
		*ksubres = deepCopyPolynomial_spX (a);
		return;
	}

	if (k > degPolynomial_spX (b)) {
		*ksubres = NULL;
		return;
	}
	
    if (isZero_spX (a)){
		if (degPolynomial_spX (b) == k) {
			*ksubres = deepCopyPolynomial_spX (b);
		} else {
			*ksubres = NULL;
		}	

		return;
    }
    
    if (isZero_spX (b)){
		if (degPolynomial_spX (a) == k) {
			*ksubres = deepCopyPolynomial_spX (a);
		} else {
			*ksubres = NULL;
		}	
		
		return;
    }	
	
	fprintf(stderr, "DUSP Error: kSubResultantYapHGCDInForm_spX... COMMENTED!\n");
	exit (1);
 ////////////////////////////// COMMENTED /////////////////////////////////
 //    polysize_t deg_a = degPolynomial_spX (a);
 //    polysize_t deg_b = degPolynomial_spX (b);
    
	// // TODO: threshold for sylvSubResultant...
	
	// duspoly_t* A;
 //    duspoly_t* B;
	
 //    elem_t mone = 0;
 //    int isSwap = 0;
    
 //    if (deg_a < deg_b) {
	// 	mone = smallprimefield_convert_in (-1, Pptr);
	// 	mone = smallprimefield_exp (mone, deg_a*deg_b, Pptr);
	
	// 	A = b;
	// 	B = a;

	// 	deg_a = deg_b;
	// 	deg_b = degPolynomial_spX (B);
	
	// 	isSwap = 1;

 //    } else {
	// 	A = a;
	// 	B = b;
 //    }

	// duspoly_t* Sm = A;
 //    duspoly_t* h = B;
 //    duspoly_t* Sn = NULL;
 //    duspoly_t* Snt = NULL;

 //    elem_t lc_h;
 //    elem_t cnst;
 //    elem_t none;
 //    elem_t total = smallprimefield_convert_in (1, Pptr);
 //    polysize_t deg_h;

	// if (deg_a == deg_b) {

	// 	h = NULL;
	// 	plainRemPolynomials_spX (A, B, &h, Pptr);

	// 	if (degPolynomial_spX (h) == k) {
			
	// 		if (deg_a & deg_b & 1) {
	// 			none = smallprimefield_convert_in (-1, Pptr);
	// 		} else {
	// 			none = smallprimefield_convert_in (1, Pptr);
	// 		}
			
	// 		cnst = smallprimefield_exp (leadingCoeffInForm_spX (B), deg_a - k, Pptr);
			
	// 		if (isSwap) {
	// 		    none = smallprimefield_mul (none, mone, Pptr);
	// 		}

	// 		cnst = smallprimefield_mul (cnst, none, Pptr);
			
	// 		scalarMulPolynomialInForm_spX (h, cnst, &Sn, Pptr);
			
	// 		freePolynomial_spX (&h);
			
	// 		if (degPolynomial_spX (Sn) == k) {
	// 			*ksubres = Sn;
	// 			return;
	// 		} else {
	// 			*ksubres = NULL;
	// 			return;
	// 		}
			
	// 	}
		
	// 	if (degPolynomial_spX (h) < k) {
	// 		*ksubres = NULL;
	// 		return;
	// 	}

	// 	fprintf (stderr, "DUSP Support: Subresultant algorithm doesnt support polynomials with the same degrees!\n");
	// 	exit (1);
		
	// 	/* Sm = B; */
	// 	/* A = B; */
	// 	/* B = h; */
		
	// 	/* deg_a = degPolynomial_spX (A); */
	// 	/* deg_b = degPolynomial_spX (B); */
		
	// }

	// duspoly_t** Quo = (duspoly_t**) calloc ((deg_b + 2), sizeof (duspoly_t*));
	// elem_t* LC = (elem_t*) calloc (deg_b + 2, sizeof (elem_t));
	// polysize_t* DD = (polysize_t*) calloc (deg_b + 2, sizeof (polysize_t));
	// polysize_t len = 0;
	// polysize_t isHalf = 1;
	
	// Quo[len] = NULL;
	// LC[len] = leadingCoeffInForm_spX (a);
	// DD[len] = deg_a;
	// len++;
	
	// dusmat4_t* MM = NULL;

	// /* YapHGCD_TotalCaseInFrom_spX	(A, B, &MM, &Quo, &LC, &DD, &len, Pptr); */	
	// extYapHGCD_TotalCaseInFrom_spX	(A, B, k, &MM, &Quo, &LC, &DD, &len, Pptr);

	// fprintf (stderr, "A := ");
	// printPolynomial_spX (convertPolynomialFromMontgomery_spX (A, Pptr), Pptr);

	// fprintf (stderr, "B := ");
	// printPolynomial_spX (convertPolynomialFromMontgomery_spX (B, Pptr), Pptr);

	// fprintf (stderr, "MM := \n");
	// printMat4OutForm_spX (MM, Pptr);

	// duspoly_t* AA = deepCopyPolynomial_spX (A);
	// duspoly_t* BB = deepCopyPolynomial_spX (B);

	// mulMat4ToVecInForm_spX_inp (MM, &AA, &BB, Pptr);

	// fprintf (stderr, "MM*A := ");
	// printPolynomial_spX (convertPolynomialFromMontgomery_spX (AA, Pptr), Pptr);

	// fprintf (stderr, "MM*B := ");
	// printPolynomial_spX (convertPolynomialFromMontgomery_spX (BB, Pptr), Pptr);

	
	// polysize_t qlen = 0;
	// for (int i = 0; i < len && DD[i] > k; i++) {
	// 	if (DD[i+1] == k) {
	// 		qlen = i+1;
	// 		break;
	// 	}
	// }

	// fprintf (stderr, "qlen = %ld\n", qlen);
	// fprintf (stderr, "k = %ld\n", k);
	
	// duspoly_t** Rem_t;

	// if (degPolynomial_spX (AA) == k) {
	// 	Rem_t = (duspoly_t**) calloc (qlen+1, sizeof (duspoly_t*));

	// 	for (int i = 0; i < qlen;i++) {
	// 		Rem_t[i] = constPolynomialInForm_spX (LC[i], Pptr);
	// 	}
	// 	Rem_t[qlen] = AA;
		
	// 	freePolynomial_spX (&BB);
	// } else if (degPolynomial_spX (BB) == k) {
	// 	Rem_t = (duspoly_t**) calloc (qlen+1, sizeof (duspoly_t*));

	// 	for (int i = 0; i < qlen;i++) {
	// 		Rem_t[i] = constPolynomialInForm_spX (LC[i], Pptr);
	// 	}
	// 	Rem_t[qlen] = BB;
		
	// 	freePolynomial_spX (&AA);		
	// } else {
	// 	Rem_t = (duspoly_t**) calloc (2, sizeof (duspoly_t*));

	// 	Rem_t[0] = A;
	// 	Rem_t[1] = B;
	
	// 	quotientsToRemainders_spX (Quo, &Rem_t, qlen, Pptr);

	// 	freePolynomial_spX (&AA);
	// 	freePolynomial_spX (&BB);
	// 	isHalf = 0;
	// }
	
	// // TEST
	// for (int i = 0; i < len; i++) {
	// 	fprintf (stderr, "Quo_t[%d] = ", i);
	// 	printPolynomial_spX (convertPolynomialFromMontgomery_spX (Quo[i], Pptr), Pptr);
	// }

	// // TEST
	// for (int i = 0; i < qlen+1; i++) {
	// 	fprintf (stderr, "Rem_t[%d] = ", i);
	// 	printPolynomial_spX (convertPolynomialFromMontgomery_spX (Rem_t[i], Pptr), Pptr);
	// }
	
	// // TEST
	// for (int i = 0; i < len; i++) {
	// 	fprintf (stderr, "LC[%d] = %lld \t DD[%d] = %ld\n", i, smallprimefield_convert_out (LC[i], Pptr), i, DD[i]);
	// }
	
	// freeMat4_spX (&MM);
	
	// for (int i = 0; i < len; i++) {
	// 	if (Quo[i] != NULL) {
	// 		freePolynomial_spX (&Quo[i]);
	// 	}
	// }
	// free (Quo);

	// /* if (!isHalf) { */
	// /* 	fprintf (stderr, "isHalf is 0?!\n"); */
	// /* 	exit(1); */
	// /* } */
	
	// /* res_t = smallprimefield_exp (LC[1], DD[0] - DD[1], Pptr); */
	
	// for (polysize_t i = 1; i < qlen; i++) {
		
	// 	/* // TODO: OutForm ... */
	// 	/* /\* if ((DD[i]-DD[i+1])*(DD[0]-DD[i+1]+i+1) % 2) { *\/ */
	// 	/* if ( DD[i] & DD[i+1] & 1) { */
	// 	/* 	mone = smallprimefield_convert_in (-1, Pptr); */
	// 	/* } else { */
	// 	/* 	mone = smallprimefield_convert_in (1, Pptr); */
	// 	/* } */
		
	// 	/* coef = smallprimefield_mul (LC[i], LC[i+1], Pptr); */
	// 	/* coef = smallprimefield_exp (coef, DD[i]-DD[i+1], Pptr); */
	// 	/* coef = smallprimefield_mul (coef, mone, Pptr); */
		
	// 	/* res_t = smallprimefield_mul (res_t, coef, Pptr); // S[DD[i+1]] */

	// 	/* /\* subR[subRsz] = constPolynomialInForm_spX (res_t, Pptr); *\/ */
	// 	/* /\* subRsz++; *\/ */
		
	// 	/* fprintf (stderr, "S[%ld] = %lld; \n", DD[i+1], smallprimefield_convert_out (res_t, Pptr)); */
		
	// 	lc_h = LC[i]; // leadingCoeffInForm_spX (h);
	// 	deg_h = DD[i]; // degPolynomial_spX (h);

	// 	Sn = Rem_t[i+1];
		
	// 	none = smallprimefield_exp (smallprimefield_convert_in (-1, Pptr),
	// 								deg_h * DD[i-1], Pptr); // DD[i-1] = degPolynomial_spX (Sm)

	// 	if (isSwap) {
	// 		none = smallprimefield_mul (none, mone, Pptr);
	// 	}

	// 	cnst = smallprimefield_exp (lc_h, DD[i-1] - DD[i+1], Pptr);// DD[i-1] = degPolynomial_spX (Sm)
	// 	                                                           // DD[i+1] = degPolynomial_spX (Sn)

	// 	cnst = smallprimefield_mul (cnst, none, Pptr);
	// 	total = smallprimefield_mul (total, cnst, Pptr);
		
	// 	scalarMulPolynomialInForm_spX (Sn, total, &Snt, Pptr);

	// 	if (degPolynomial_spX (Snt) == k) {
	// 		*ksubres = Snt;

	// 		break;
	// 	} /* else if (degPolynomial_spX (Snt) < k) { */
	// 	/* 	*ksubres = NULL; */

	// 	/* 	break; */
	// 	/* } */
		
	// 	Sm = h;
	// 	h = Sn;
	// 	Snt = NULL;
	// 	Sn = NULL;
		
		
	// }

	// if (isHalf) {
	// 	freePolynomial_spX (&Rem_t[0]);
	// 	freePolynomial_spX (&Rem_t[1]);
	// }
	
	// for (int i = 2; i < qlen; i++) {
	// 	if (Rem_t[i] != NULL) {
	// 		freePolynomial_spX (&Rem_t[i]);
	// 	}
	// }
	
	// free (Rem_t);

	// free (LC);
	// free (DD);

}

void MCAHalfGCDInForm_spX (duspoly_t* a, duspoly_t* b, polysize_t mca_k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, polysize_t* h, const Prime_ptr* Pptr)
{
	
    dusmat4_t* E = NULL;
    identityMat4InForm_spX (&E, Pptr);
	
    if (isZero_spX (b) || isZero_spX (a)) {
		
		*h = 0;
		*M = E;
		
		return;
    }
	
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
	
	if (mca_k < deg_a - deg_b) { // was < in MCA
		
		*h = 0;
		*M = E;
		
		return;
	}

	if (mca_k == 0 || deg_a == deg_b) {
		*h = 1;
		
		elem_t clc = smallprimefield_div (leadingCoeffInForm_spX (a), leadingCoeffInForm_spX (b), Pptr);
		clc = smallprimefield_inv (clc, Pptr);
		
		duspoly_t* r_t = NULL;
		duspoly_t* q_t = NULL;
		
		(*LC)[(*len)] = leadingCoeffInForm_spX (b);
		(*DD)[(*len)] = degPolynomial_spX (b);
		
		plainDivPolynomialsInForm_spX (a, b, &r_t, &q_t, Pptr);
		
		(*Quo)[(*len)] = q_t;
		(*len) = (*len) + 1;

		freePolynomial_spX (&r_t);
		
		duspoly_t* t1 = E->polys[0];
		duspoly_t* t2 = E->polys[3];
		E->polys[0] = NULL;
		E->polys[1] = t1;
		E->polys[2] = t2;
		E->polys[3] = constPolynomialInForm_spX (clc, Pptr);

		*M = E;

		return;
	}
	
	polysize_t d = (mca_k%2) ? (mca_k>>1)+1 : (mca_k>>1);
	polysize_t m = deg_a - 2*d + 2;

	
	/* fprintf (stderr, "d = %ld\n",d); */
	/* fprintf (stderr, "m = %ld\n",m); */
	
    duspoly_t* A0 = NULL;
    duspoly_t* B0 = NULL;
    dusmat4_t* RR = NULL;

    polysize_t h1 = 0;
 
	rightShiftPolynomial_spX (a, &A0, m);
    rightShiftPolynomial_spX (b, &B0, m); // deg_b - ((2*m)-2-deg_a+deg_b)
	
    MCAHalfGCDInForm_spX (A0, B0, d-1, &RR, Quo, LC, DD, len, &h1, Pptr);  // first recursive call

	freePolynomial_spX (&A0); A0 = NULL;
	freePolynomial_spX (&B0); B0 = NULL;
	
	rightShiftPolynomial_spX (a, &A0, deg_a - 2*mca_k);
    rightShiftPolynomial_spX (b, &B0, deg_a - 2*mca_k); 
	
	mulMat4ToVecInForm_spX_inp (RR, &A0, &B0, Pptr);

	polysize_t delta = degPolynomial_spX (RR->polys[3]);

	polysize_t ds = mca_k - delta - degPolynomial_spX (A0) + degPolynomial_spX (B0); // degPolynomial_spX (B0) - degPolynomial_spX (a) + m + d;
		
	/* fprintf (stderr, "ds = %ld\n", ds); */
	/* fprintf (stderr, "delta = %ld\n", delta); */
	
	if (isZero_spX (B0) || ds < 0) { 

		freePolynomial_spX (&A0);
		freePolynomial_spX (&B0);
		
		freeMat4_spX (&E);
		
		*h = h1;
		*M = RR;
		return;
	}
	
	(*LC)[(*len)] = leadingCoeffInForm_spX (B0);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (A0) + degPolynomial_spX (B0);
	
    duspoly_t* q = NULL;
    duspoly_t* D = NULL;
	
    plainDivPolynomialsInForm_spX (A0, B0, &D, &q, Pptr); // TODO : Fast Monic Div
	
	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);
	(*len) = (*len) + 1;
	
    freePolynomial_spX (&A0);

	negPolynomialInForm_spX (q, Pptr);

	dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomial_spX (1, Pptr); // constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = q;
	
    if (isZero_spX (D)) {
			
		dusmat4_t* MM = NULL;

		mulMat4ToMat4InForm_spX_inp (Q, RR, &MM, Pptr);

		freeMat4_spX (&RR);
		freeMat4_spX (&E);
		freePolynomial_spX (&B0);
		
		*h = h1;
		*M = MM;
	
		return;
    }
	
    duspoly_t* C0 = NULL;
    duspoly_t* D0 = NULL;
    dusmat4_t* S  = NULL;
	polysize_t h2 = 0;

	/* fprintf (stderr, "deg(rj)-2*ds = %ld\n", degPolynomial_spX (B0) - 2*ds); */
	
    rightShiftPolynomial_spX (B0, &C0, degPolynomial_spX (B0) - 2*ds);
    rightShiftPolynomial_spX (D ,&D0, degPolynomial_spX (B0) - 2*ds);

    freePolynomial_spX (&B0);
    freePolynomial_spX (&D);
	
    MCAHalfGCDInForm_spX (C0, D0, ds, &S, Quo, LC, DD, len, &h2, Pptr);
	    
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);
	
    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;
	
	mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr);
	
    mulMat4ToMat4InForm_spX_inp (S, QRR, &SQRR, Pptr);
    
    freeMat4_spX (&RR);
    freeMat4_spX (&QRR);
	
	*h = h1 + h2 + 1;
    *M = SQRR;

    return;
}

void MCAHGCD_BaseCaseInForm_spX (duspoly_t* a, duspoly_t* b, polysize_t k, polysize_t mca_k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, polysize_t* h, const Prime_ptr* Pptr)
{
	
    dusmat4_t* E = NULL;
    identityMat4InForm_spX (&E, Pptr);
	
    if (isZero_spX (b) || isZero_spX (a)) {
		
		*h = 0;
		*M = E;
		
		return;
    }
	
	if ( k > -1 && (*len) > 0 && (*DD)[(*len)-1] < k) {
		*M = E;
		return;
	}

    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
	
	if (mca_k < deg_a - deg_b) { // was < in MCA
		
		*h = 0;
		*M = E;
		
		return;
	}

	if (mca_k == 0 || deg_a == deg_b) {
		*h = 1;
		
		elem_t clc = smallprimefield_div (leadingCoeffInForm_spX (a), leadingCoeffInForm_spX (b), Pptr);
		clc = smallprimefield_inv (clc, Pptr);
		
		duspoly_t* r_t = NULL;
		duspoly_t* q_t = NULL;
		
		(*LC)[(*len)] = leadingCoeffInForm_spX (b);
		(*DD)[(*len)] = degPolynomial_spX (b);
		
		plainDivPolynomialsInForm_spX (a, b, &r_t, &q_t, Pptr);
		
		(*Quo)[(*len)] = q_t;
		(*len) = (*len) + 1;

		freePolynomial_spX (&r_t);
		
		duspoly_t* t1 = E->polys[0];
		duspoly_t* t2 = E->polys[3];
		E->polys[0] = NULL;
		E->polys[1] = t1;
		E->polys[2] = t2;
		E->polys[3] = constPolynomialInForm_spX (clc, Pptr);

		*M = E;

		return;
	}
	
	polysize_t d = (mca_k%2) ? (mca_k>>1)+1 : (mca_k>>1);
	
	fprintf (stderr, "d = %ld\n",d);

	polysize_t m = deg_a - 2*d + 2;
	
	fprintf (stderr, "m = %ld\n",m);
	
    duspoly_t* A0 = NULL;
    duspoly_t* B0 = NULL;
    dusmat4_t* RR = NULL;

    polysize_t h1 = 0;
 
	rightShiftPolynomial_spX (a, &A0, m);
    rightShiftPolynomial_spX (b, &B0, m); // deg_b - ((2*m)-2-deg_a+deg_b)

	fprintf (stderr, "A0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (A0, Pptr), Pptr);

	fprintf (stderr, "B0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (B0, Pptr), Pptr);
	
    MCAHGCD_BaseCaseInForm_spX (A0, B0, k, d-1, &RR, Quo, LC, DD, len, &h1, Pptr);  // first recursive call

	freePolynomial_spX (&A0); A0 = NULL;
	freePolynomial_spX (&B0); B0 = NULL;
	
	rightShiftPolynomial_spX (a, &A0, deg_a - 2*mca_k);
    rightShiftPolynomial_spX (b, &B0, deg_a - 2*mca_k); 

	
	fprintf (stderr, "[NEW] A0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (A0, Pptr), Pptr);

	fprintf (stderr, "[NEW] B0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (B0, Pptr), Pptr);
	
	mulMat4ToVecInForm_spX_inp (RR, &A0, &B0, Pptr);

	polysize_t delta = degPolynomial_spX (RR->polys[3]);

	polysize_t ds = mca_k - delta - degPolynomial_spX (A0) + degPolynomial_spX (B0); // degPolynomial_spX (B0) - degPolynomial_spX (a) + m + d;
	// OR mca_k - degPolynomial_spX (RR->polys[3]) - degPolynomial_spX (A) + degPolynomial_spX (B);

	fprintf (stderr, "RR = \n");
	printMat4OutForm_spX (RR, Pptr);
	
	fprintf (stderr, "ds = %ld\n", ds);
	fprintf (stderr, "delta = %ld\n", delta);
	
	if (isZero_spX (B0) || ds < 0) { 
		// OR if mca_k < degPolynomial_spX (RR->polys[3]) + degPolynomial_spX (A) - degPolynomial_spX (B)		
		freePolynomial_spX (&A0);
		freePolynomial_spX (&B0);
		
		freeMat4_spX (&E);
		
		*h = h1;
		*M = RR;
		return;
	}
	
	(*LC)[(*len)] = leadingCoeffInForm_spX (B0);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (A0) + degPolynomial_spX (B0);

    duspoly_t* q = NULL;
    duspoly_t* D = NULL;
	
    plainDivPolynomialsInForm_spX (A0, B0, &D, &q, Pptr); // TODO : Fast Monic Div

	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);
	(*len) = (*len) + 1;

    freePolynomial_spX (&A0);

	negPolynomialInForm_spX (q, Pptr);

	dusmat4_t* Q = (dusmat4_t*) malloc (sizeof (dusmat4_t));
    Q->polys[0] = NULL;
    Q->polys[1] = constPolynomial_spX (1, Pptr);
    Q->polys[2] = constPolynomial_spX (1, Pptr); // constPolynomialInForm_spX (lc, Pptr);
    Q->polys[3] = q;

    if (isZero_spX (D)) {
			
		dusmat4_t* MM = NULL;

		mulMat4ToMat4InForm_spX_inp (Q, RR, &MM, Pptr);

		freeMat4_spX (&RR);
		freeMat4_spX (&E);
		freePolynomial_spX (&B0);
	
		*h = h1;
		*M = MM;
	
		return;
    }
	
	// Uncommented if you want to use exact MCA algorithm:
	// polysize_t ds = mca_k - degPolynomial_spX (RR->polys[3]) - degPolynomial_spX (A) + degPolynomial_spX (B); 

    duspoly_t* C0 = NULL;
    duspoly_t* D0 = NULL;
    dusmat4_t* S  = NULL;
	polysize_t h2 = 0;

	fprintf (stderr, "deg(rj)-2*ds = %ld\n", degPolynomial_spX (B0) - 2*ds);
	
	fprintf (stderr, "[next] B0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (B0, Pptr), Pptr);
	
	fprintf (stderr, "[next] D = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (D, Pptr), Pptr);
	
    rightShiftPolynomial_spX (B0, &C0, degPolynomial_spX (B0) - 2*ds);
    rightShiftPolynomial_spX (D ,&D0, degPolynomial_spX (B0) - 2*ds);

	fprintf (stderr, "[next] C0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (C0, Pptr), Pptr);
	
	fprintf (stderr, "[next] D0 = ");
	printPolynomial_spX (convertPolynomialFromMontgomery_spX (D0, Pptr), Pptr);

    freePolynomial_spX (&B0);
    freePolynomial_spX (&D);
	
    MCAHGCD_BaseCaseInForm_spX (C0, D0, k, ds, &S, Quo, LC, DD, len, &h2, Pptr); // second recursive call
	    
    freePolynomial_spX (&C0);
    freePolynomial_spX (&D0);
	
	fprintf (stderr, "S = \n");
	printMat4OutForm_spX (S, Pptr);

    dusmat4_t* QRR  = NULL;
    dusmat4_t* SQRR = NULL;
	
	mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr);
	
    mulMat4ToMat4InForm_spX_inp (S, QRR, &SQRR, Pptr);
    
    freeMat4_spX (&RR);
    freeMat4_spX (&QRR);
	
	*h = h1 + h2 + 1;
    *M = SQRR;

	fprintf (stderr, "return = \n");
	printMat4OutForm_spX (SQRR, Pptr);

    return;
}

void MCAHGCD_TotalCaseInFrom_spX (duspoly_t* a, duspoly_t* b, polysize_t k, dusmat4_t** M, duspoly_t*** Quo, elem_t** LC, polysize_t** DD, polysize_t* len, const Prime_ptr* Pptr)
{
	if ((*len) > 0 && (*DD)[(*len)-1] < k) {
		
		dusmat4_t* E = NULL;
		identityMat4InForm_spX (&E, Pptr);
		
		*M = E;
		return;
	}
	
	polysize_t mca_k = degPolynomial_spX (a);
	polysize_t h = 0;
    dusmat4_t* RR = NULL;
    MCAHGCD_BaseCaseInForm_spX (a, b, k, mca_k, &RR, Quo, LC, DD, len, &h, Pptr);
    
    // TODO : implement non in-place mulMat4ToVecInForm 
    duspoly_t* a0 = deepCopyPolynomial_spX (a);
    duspoly_t* b0 = deepCopyPolynomial_spX (b);
    
    mulMat4ToVecInForm_spX_inp (RR, &a0, &b0, Pptr);
    
    if (isZero_spX (b0)) {
		*M = RR;
		return;
    }
	
	(*LC)[(*len)] = leadingCoeffInForm_spX (b0);	
	(*DD)[(*len)] = (*DD)[(*len)-1] - degPolynomial_spX (a0) + degPolynomial_spX (b0);
	
    duspoly_t* r   = NULL;
    duspoly_t* q   = NULL;
    dusmat4_t* Q   = (dusmat4_t*) malloc (sizeof(dusmat4_t));
    
    plainDivPolynomialsInForm_spX (a0, b0, &r, &q, Pptr); // TODO: Fast Monic Divid

	(*Quo)[(*len)] = deepCopyPolynomial_spX (q);	
	(*len) = (*len) + 1;	

	negPolynomialInForm_spX (q, Pptr); // q = -q;
    
    freePolynomial_spX (&a0);
	
	Q->polys[0] = NULL;
	Q->polys[1] = constPolynomial_spX (1, Pptr);
	Q->polys[2] = constPolynomial_spX (1, Pptr);
	Q->polys[3] = q; // q = -q
	
    if (isZero_spX (r)) {
	
		freePolynomial_spX (&b0);
		
		dusmat4_t* QRR = NULL;
		// mulMat4ToMat4InForm_spX_inp_inp (Q, &RR, Pptr);
		mulMat4ToMat4InForm_spX_inp (Q, RR, &QRR, Pptr); 
	
		freeMat4_spX (&RR);
	
		*M = QRR;
		return;
    }
	
    dusmat4_t* S = NULL;
	mca_k = degPolynomial_spX (b0);
	
    MCAHGCD_BaseCaseInForm_spX (b0, r, k, mca_k, &S, Quo, LC, DD, len, &h, Pptr);

    freePolynomial_spX (&b0);
    freePolynomial_spX (&r);

    dusmat4_t* RRQ = NULL;
    dusmat4_t* SRRQ = NULL;
    
    mulMat4ToMat4InForm_spX_inp (Q, RR, &RRQ, Pptr);
    mulMat4ToMat4InForm_spX_inp (S, RRQ, &SRRQ, Pptr);

    freeMat4_spX (&RR);
    freeMat4_spX (&RRQ);
    
    *M = SRRQ;     
}

void kSubResultantMCAHGCDInForm_spX (duspoly_t* a, duspoly_t* b, polysize_t k, duspoly_t** ksubres, const Prime_ptr* Pptr)
{

	if (k < 0) {
		fprintf (stderr, "DUSP Error: k is out of renge in subResutlant algorithm!\n");
		exit (1);
	}

	if (k == 0) {
		resultantYapHGCDInForm_spX (a, b, ksubres, Pptr);
		return;
	}

	if (k == degPolynomial_spX (b)) {
		*ksubres = deepCopyPolynomial_spX (b);
		return;
	}

	if (k == degPolynomial_spX (a)) {
		*ksubres = deepCopyPolynomial_spX (a);
		return;
	}

	if (k > degPolynomial_spX (b)) {
		*ksubres = NULL;
		return;
	}
	
    if (isZero_spX (a)){
		if (degPolynomial_spX (b) == k) {
			*ksubres = deepCopyPolynomial_spX (b);
		} else {
			*ksubres = NULL;
		}	

		return;
    }
    
    if (isZero_spX (b)){
		if (degPolynomial_spX (a) == k) {
			*ksubres = deepCopyPolynomial_spX (a);
		} else {
			*ksubres = NULL;
		}	
		
		return;
    }	
	
    polysize_t deg_a = degPolynomial_spX (a);
    polysize_t deg_b = degPolynomial_spX (b);
    
	// TODO: threshold for sylvSubResultant...
	
	duspoly_t* A;
    duspoly_t* B;
	
    elem_t mone = 0;
    int isSwap = 0;
    
    if (deg_a < deg_b) {
		mone = smallprimefield_convert_in (-1, Pptr);
		mone = smallprimefield_exp (mone, deg_a*deg_b, Pptr);
	
		A = b;
		B = a;

		deg_a = deg_b;
		deg_b = degPolynomial_spX (B);
	
		isSwap = 1;

    } else {
		A = a;
		B = b;
    }

	duspoly_t* Sm = A;
    duspoly_t* h = B;
    duspoly_t* Sn = NULL;
    duspoly_t* Snt = NULL;

    elem_t lc_h;
    elem_t cnst;
    elem_t none;
    elem_t total = smallprimefield_convert_in (1, Pptr);
    polysize_t deg_h;

	if (deg_a == deg_b) {

		h = NULL;
		plainRemPolynomials_spX (A, B, &h, Pptr);

		if (degPolynomial_spX (h) == k) {
			
			if (deg_a & deg_b & 1) {
				none = smallprimefield_convert_in (-1, Pptr);
			} else {
				none = smallprimefield_convert_in (1, Pptr);
			}
			
			cnst = smallprimefield_exp (leadingCoeffInForm_spX (B), deg_a - k, Pptr);
			
			if (isSwap) {
			    none = smallprimefield_mul (none, mone, Pptr);
			}

			cnst = smallprimefield_mul (cnst, none, Pptr);
			
			scalarMulPolynomialInForm_spX (h, cnst, &Sn, Pptr);
			
			freePolynomial_spX (&h);
			
			if (degPolynomial_spX (Sn) == k) {
				*ksubres = Sn;
				return;
			} else {
				*ksubres = NULL;
				return;
			}
			
		}
		
		if (degPolynomial_spX (h) < k) {
			*ksubres = NULL;
			return;
		}

		fprintf (stderr, "DUSP Support: Subresultant algorithm doesnt support polynomials with the same degrees!\n");
		exit (1);
		
		/* Sm = B; */
		/* A = B; */
		/* B = h; */
		
		/* deg_a = degPolynomial_spX (A); */
		/* deg_b = degPolynomial_spX (B); */
		
	}

	duspoly_t** Quo = (duspoly_t**) calloc ((deg_b + 2), sizeof (duspoly_t*));
	elem_t* LC = (elem_t*) calloc (deg_b + 2, sizeof (elem_t));
	polysize_t* DD = (polysize_t*) calloc (deg_b + 2, sizeof (polysize_t));
	polysize_t len = 0;
	polysize_t hh = 0;
	
	Quo[len] = NULL;
	LC[len] = leadingCoeffInForm_spX (a);
	DD[len] = deg_a;
	len++;
	
	dusmat4_t* MM = NULL;
	polysize_t mca_k = degPolynomial_spX (A);
	/* YapHGCD_TotalCaseInFrom_spX	(A, B, &MM, &Quo, &LC, &DD, &len, Pptr); */	
	/* extYapHGCD_TotalCaseInFrom_spX	(A, B, k, &MM, &Quo, &LC, &DD, &len, Pptr); */
	MCAHGCD_BaseCaseInForm_spX (A, B, k, mca_k, &MM, &Quo, &LC, &DD, &len, &hh, Pptr);
		
	polysize_t qlen = 0;
	for (int i = 0; i < len && DD[i] > k; i++) {
		if (DD[i+1] == k) {
			qlen = i+1;
			break;
		}
	}

	/* fprintf (stderr, "qlen = %d\n", qlen); */
	/* fprintf (stderr, "k = %ld\n", k); */
	
	duspoly_t** Rem_t = (duspoly_t**) calloc (2, sizeof (duspoly_t*));
	Rem_t[0] = A;
	Rem_t[1] = B;
	
	quotientsToRemainders_spX (Quo, &Rem_t, qlen, Pptr);

	/* // TEST */
	/* for (int i = 0; i < qlen+1; i++) { */
	/* 	fprintf (stderr, "Quo_t[%d] = ", i); */
	/* 	printPolynomial_spX (convertPolynomialFromMontgomery_spX (Quo[i], Pptr), Pptr); */
	/* } */

	/* // TEST */
	/* for (int i = 0; i < qlen+1; i++) { */
	/* 	fprintf (stderr, "Rem_t[%d] = ", i); */
	/* 	printPolynomial_spX (convertPolynomialFromMontgomery_spX (Rem_t[i], Pptr), Pptr); */
	/* } */
	
	/* // TEST */
	/* for (int i = 0; i < len; i++) { */
	/* 	fprintf (stderr, "LC[%d] = %lld \t DD[%d] = %ld\n", i, smallprimefield_convert_out (LC[i], Pptr), i, DD[i]); */
	/* } */
	
	freeMat4_spX (&MM);
	
	for (int i = 0; i < len; i++) {
		if (Quo[i] != NULL) {
			freePolynomial_spX (&Quo[i]);
		}
	}
	free (Quo);
	
	/* res_t = smallprimefield_exp (LC[1], DD[0] - DD[1], Pptr); */
	
	for (polysize_t i = 1; i < qlen; i++) {
		
		/* // TODO: OutForm ... */
		/* /\* if ((DD[i]-DD[i+1])*(DD[0]-DD[i+1]+i+1) % 2) { *\/ */
		/* if ( DD[i] & DD[i+1] & 1) { */
		/* 	mone = smallprimefield_convert_in (-1, Pptr); */
		/* } else { */
		/* 	mone = smallprimefield_convert_in (1, Pptr); */
		/* } */
		
		/* coef = smallprimefield_mul (LC[i], LC[i+1], Pptr); */
		/* coef = smallprimefield_exp (coef, DD[i]-DD[i+1], Pptr); */
		/* coef = smallprimefield_mul (coef, mone, Pptr); */
		
		/* res_t = smallprimefield_mul (res_t, coef, Pptr); // S[DD[i+1]] */

		/* /\* subR[subRsz] = constPolynomialInForm_spX (res_t, Pptr); *\/ */
		/* /\* subRsz++; *\/ */
		
		/* fprintf (stderr, "S[%ld] = %lld; \n", DD[i+1], smallprimefield_convert_out (res_t, Pptr)); */
		
		lc_h = leadingCoeffInForm_spX (h);
		deg_h = degPolynomial_spX (h);

		Sn = Rem_t[i+1];
		none = smallprimefield_exp (smallprimefield_convert_in (-1, Pptr),
									deg_h * degPolynomial_spX (Sm), Pptr);

		if (isSwap) {
			none = smallprimefield_mul (none, mone, Pptr);
		}

		cnst = smallprimefield_exp (lc_h, degPolynomial_spX (Sm) -
									degPolynomial_spX (Sn), Pptr);

		cnst = smallprimefield_mul (cnst, none, Pptr);
		total = smallprimefield_mul (total, cnst, Pptr);
		
		scalarMulPolynomialInForm_spX (Sn, total, &Snt, Pptr);

		if (degPolynomial_spX (Snt) == k) {
			*ksubres = Snt;

			break;
		} else if (degPolynomial_spX (Snt) < k) {
			*ksubres = NULL;

			break;
		}
		
		Sm = h;
		h = Sn;
		Snt = NULL;
		Sn = NULL;
		
		
	}

	for (int i = 2; i < qlen; i++) {
		if (Rem_t[i] != NULL) {
			freePolynomial_spX (&Rem_t[i]);
		}
	}
	
	free (Rem_t);

	free (LC);
	free (DD);

}




/****************************
 ** Polynomial Factorization 
 ****************************/

// See DUSP_Support_Factoring.c and DUSP_FFT_Support_Factoring.cpp


/******************************************
 * Bivariate Polynomials over Small Primes
 *****************************************/

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
elem_t _evalDivPolynomials_spX (elem_t* a, elem_t u, elem_t* q, polysize_t n, const Prime_ptr* Pptr)
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
void _interpolatePartialSubResInForm_spXY (duspolysA_t** pSubres, biSubresPr_t* XYsubres, elem_t* t, polysize_t n, polysize_t alpha, const Prime_ptr* Pptr)
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
        // fprintf (stderr, "g is computed for k = %ld and delta = %ld\n", k, delta);
        // for (polysize_t sss=0; sss <= delta; sss++) {
        //     fprintf (stderr, "g[%ld] = %lld\n", offset+k+sss, smallprimefield_convert_out (g[offset+k+sss], Pptr)); }
        // for each specialized subresultant chain at row k
        for (polysize_t i = 0; i < delta; i++) {
            // fprintf (stderr, "pSubres[%ld]->polys[%ld] \n", i, k);
            // fprintf (stderr, "pSubres[i] = %p\n", pSubres[i]);
            // fprintf (stderr, "delta = %ld\n", delta);
            // fprintf (stderr, "pSubres[i]->polys = %p\n", pSubres[i]->polys);
            // fprintf (stderr, "pSubres[i]->size = %ld\n", pSubres[i]->size);
            // printPolynomialOutForm_spX (pSubres[14]->polys[1], Pptr);  // TEST 
            if (isZero_spX (pSubres[i]->polys[k])) {
                continue;
            }
            elem_t e = _evalDivPolynomials_spX (&g[offset+k], t[i], &q[offset], delta, Pptr);
            e = smallprimefield_inv (e, Pptr);
            // fprintf (stderr, "q is computed for k = %ld and i = %ld\n", k, i);
            // fprintf (stderr, "e = %lld\n", smallprimefield_convert_out (e, Pptr));
            // for (polysize_t sss=0; sss <= delta; sss++) {
            //     fprintf (stderr, "q[%ld] = %lld\n", offset+sss, smallprimefield_convert_out (q[offset+sss], Pptr));
            // }
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
                            XYsubres->coefs[k][s*delta+j] = smallprimefield_add (XYsubres->coefs[k][s*delta+j] , smallprimefield_mul (cc, q[offset+j], Pptr), Pptr);
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
		a_dusp->elems = (elem_t*) malloc (a_dusp->alloc * sizeof(elem_t)); // &a[k*as];
        a_dusp->elems = memcpy (a_dusp->elems, &a[k*as], a_dusp->alloc * sizeof(elem_t));
		b_dusp->elems = (elem_t*) malloc (b_dusp->alloc * sizeof(elem_t)); // &b[k*bs];
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
        brownSubResultantInForm_spX (a_dusp, b_dusp, &SS, &sz, Pptr);
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
    *subres = s;
}



