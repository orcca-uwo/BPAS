#include "ModularPolynomial/DUSP_Support_Test.h"

// print polynomial
// TODO delete prp
void printPolynomial_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{

    if (a == NULL || POLY_ALLOC(a) == 0) {
    	fprintf(stderr, "  ;\n");
    	return;
    }

    elem_t* elems = a->elems;

    fprintf(stderr, " ");
    if (POLY_LT(a) == 0) {
    	fprintf(stderr, "%lld in Z_{%lld}[x];\n", elems[0], Pptr->prime);
    	return;
    }

    if (elems[0] != 0) {
        fprintf(stderr, "%lld + ", elems[0]);
    }

    for (polysize_t i = 1; i < POLY_LT(a); ++i) {
	    if (elems[i] != 0) {
           fprintf(stderr, "%lld*x^%ld + ", elems[i], i);
        }
    }

    fprintf(stderr, "%lld*x^%ld in Z_{%lld}[x];\n", elems[POLY_LT(a)], POLY_LT(a), Pptr->prime);
}

void purePrintPolySymOutForm_spX (const duspoly_t* a, const Prime_ptr* Pptr, const char* sym)
{

    if (a == NULL || POLY_ALLOC(a) == 0) {
        fprintf(stderr, " ;");
        return;
    }

    elem_t* elems = a->elems;

    if (POLY_LT(a) == 0) {
        fprintf(stderr, "%lld", smallprimefield_convert_out(elems[0], Pptr));
        return;
    }

    if (elems[0] != 0) {
        fprintf(stderr, "%lld + ", smallprimefield_convert_out(elems[0], Pptr));
    }

    for (polysize_t i = 1; i < POLY_LT(a); ++i) {
        if (elems[i] != 0) {
           fprintf(stderr, "%lld*%s^%ld + ", smallprimefield_convert_out(elems[i], Pptr), sym, i);
        }
    }

    fprintf(stderr, "%lld*%s^%ld", smallprimefield_convert_out(elems[POLY_LT(a)], Pptr), sym, POLY_LT(a));
}


void purePrintPolynomial_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{

    if (a == NULL || POLY_ALLOC(a) == 0) {
        fprintf(stderr, " ");
        return;
    }

    elem_t* elems = a->elems;

    fprintf(stderr, " ");
    if (POLY_LT(a) == 0) {
        fprintf(stderr, "%lld", elems[0]);
        return;
    }

    if (elems[0] != 0) {
        fprintf(stderr, "%lld + ", elems[0]);
    }

    for (polysize_t i = 1; i < POLY_LT(a); ++i) {
        if (elems[i] != 0) {
            fprintf(stderr, "%lld*x^%ld + ", elems[i], i);
        }
    }

    fprintf(stderr, "%lld*x^%ld ", elems[POLY_LT(a)], POLY_LT(a));
}

void printPolynomialOutForm_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{
    duspoly_t* a_out = convertPolynomialFromMontgomery_spX (a, Pptr);
    printPolynomial_spX (a_out, Pptr);
    freePolynomial_spX (&a_out);
}

void purePrintPolynomialOutForm_spX (const duspoly_t* a, const Prime_ptr* Pptr)
{
    duspoly_t* a_out = convertPolynomialFromMontgomery_spX (a, Pptr);
    purePrintPolynomial_spX (a_out, Pptr);
    freePolynomial_spX (&a_out);
}

void printPolynomialOutForm_spXY (const dbspoly_t* a, const Prime_ptr* Pptr) {
    if (a == NULL) {
        return;
    }


    for (int i = 0; i < a->ltX; ++i) {
        fprintf(stderr, "(");
        purePrintPolySymOutForm_spX (a->coefsY[i], Pptr, "y");
        fprintf(stderr, ")*x^%d + ", i );
    }
    fprintf(stderr, "(");
    purePrintPolySymOutForm_spX (a->coefsY[a->ltX], Pptr, "y");
    fprintf(stderr, ")*x^%ld in Z_{%lld}[x,y]\n", a->ltX, Pptr->prime);

}

// print Matrix 2*2 of polynomials
// TODO delete prp
void printMat4_spX (dusmat4_t* A, const Prime_ptr* Pptr)
{
    // fprintf(stderr, "printPolynomialInForm_spX... \n");
    if (A == NULL) {
	fprintf(stderr, "  ;\n");
	return;
    }

    fprintf (stderr, "A[0,0] := ");
    printPolynomial_spX (A->polys[0], Pptr);

    fprintf (stderr, "A[0,1] := ");
    printPolynomial_spX (A->polys[1], Pptr);

    fprintf (stderr, "A[1,0] := ");
    printPolynomial_spX (A->polys[2], Pptr);

    fprintf (stderr, "A[1,1] := ");
    printPolynomial_spX (A->polys[3], Pptr);

}

void printMat4OutForm_spX (dusmat4_t* A, const Prime_ptr* Pptr)
{
    // fprintf(stderr, "printPolynomialInForm_spX... \n");
    if (A == NULL) {
	fprintf(stderr, "  ;\n");
	return;
    }

    fprintf (stderr, "A[0,0] := ");
    printPolynomialOutForm_spX (A->polys[0], Pptr);

    fprintf (stderr, "A[0,1] := ");
    printPolynomialOutForm_spX (A->polys[1], Pptr);

    fprintf (stderr, "A[1,0] := ");
    printPolynomialOutForm_spX (A->polys[2], Pptr);

    fprintf (stderr, "A[1,1] := ");
    printPolynomialOutForm_spX (A->polys[3], Pptr);

}

void printFactors_spX (factors_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactors_spX (f)) {
        fprintf(stderr, " ;\n");
    }

    for (polysize_t i = 0; i < f->alloc; i++) {
        fprintf(stderr, "Factor[%lu] := ", i);
        printPolynomial_spX (f->polys[i], Pptr);
    }
}

void printFactorsOutForm_spX (factors_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactors_spX (f)) {
        fprintf(stderr, " ;\n");
    }

    for (polysize_t i = 0; i < f->alloc; i++) {
        fprintf(stderr, "Factor[%lu] := ", i);
        printPolynomialOutForm_spX (f->polys[i], Pptr);
    }
}


void printFactsll_spX (factsll_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactsll_spX (f)) {
        fprintf(stderr, " ;\n");
        return;
    }

    factsll_t* cur = f;
    while (cur != NULL) {
        fprintf(stderr, "Factll[%lu] := ", cur->exp);
        printPolynomial_spX (cur->poly, Pptr);
        cur = cur->next;
    }
}

void purePrintFactsll_spX (factsll_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactsll_spX (f)) {
        fprintf(stderr, " ;\n");
        return;
    }

    fprintf(stderr, "(");
    factsll_t* cur = f;
    while (cur != NULL) {
        // fprintf(stderr, "Factll[%lu] := ", cur->exp);
        purePrintPolynomial_spX (cur->poly, Pptr);
        fprintf(stderr, ")^%lu ", cur->exp);

        if (cur->next != NULL) {
            fprintf(stderr, "* (");
        }

        cur = cur->next;
    }
    fprintf(stderr, ";\n");
}

void printFactsllOutForm_spX (factsll_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactsll_spX (f)) {
        fprintf(stderr, " ;\n");
        return;
    }

    factsll_t* cur = f;
    while (cur != NULL) {
        fprintf(stderr, "Factll[%lu] := ", cur->exp);
        printPolynomialOutForm_spX (cur->poly, Pptr);
        cur = cur->next;
    }
}


void purePrintFactsllOutForm_spX (factsll_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactsll_spX (f)) {
        fprintf(stderr, " ;\n");
        return;
    }

    fprintf(stderr, "(");
    factsll_t* cur = f;
    while (cur != NULL) {
        // fprintf(stderr, "Factll[%lu] := ", cur->exp);
        purePrintPolynomialOutForm_spX (cur->poly, Pptr);
        fprintf(stderr, ")^%lu ", cur->exp);

        if (cur->next != NULL) {
            fprintf(stderr, "* (");
        }

        cur = cur->next;
    }
    fprintf(stderr, ";\n");
}

void purePrintFactorsOutForm_spX (factors_t* f, const Prime_ptr* Pptr)
{
    if (isZeroFactors_spX (f)) {
        fprintf(stderr, " ;\n");
        return;
    }

    fprintf(stderr, "(");
    for(polysize_t i = 0; i < f->alloc; i++) {
        // fprintf(stderr, "Factll[%lu] := ", cur->exp);
        purePrintPolynomialOutForm_spX (f->polys[i], Pptr);
        fprintf(stderr, ")^%lu ", f->exps[i]);

        if (i < f->alloc -1 ) {
            fprintf(stderr, "* (");
        }

    }
    fprintf(stderr, ";\n");
}

/**
 * Random Polynomial
 *
 * @param sz size of polynomial (deg = sz-1)
 * @param sparsity the maximum distance between two successive terms
 * @param Pptr small prime pointer
 */
duspoly_t* randomFullyDensePolynomial_spX (polysize_t sz, const Prime_ptr* Pptr)
{
    if (sz == 0) {
	return NULL;
    }
    polysize_t s = 2; // (sparsity < 2) ? 2: sparsity;

    duspoly_t* poly = makePolynomial_spX (sz);
    POLY_LT(poly) = sz-1;

    elem_t* elems = poly->elems;

    static int initTest = 0;

    while (!initTest) {
        time_t t = time(NULL);
        srand(t);

        initTest = 1;
    }

    //polysize_t step = ((polysize_t) time() % s) > 0 ? ((polysize_t) time() % s) : 1;
    polysize_t step  = s-1;

    for (polysize_t i = 0; i < sz; i += step) {
	elems[i] = (elem_t) rand() % Pptr->prime;
	// smallprimefield_covert_in ((elem_t) rand(), pr, R);
	// step = ((polysize_t) rand() % s) > 0 ? ((polysize_t) time() % s) : 1;
    }

    while (elems[sz-1] == 0) {
	elems[sz-1] = (elem_t) rand() % Pptr->prime;
    }

    normalize_spX (&poly);

    return poly;
}

int DivisionTest (const duspoly_t* a, const duspoly_t* b, const duspoly_t* r, const duspoly_t* q, const Prime_ptr* Pptr) {

    duspoly_t* right = NULL;
    duspoly_t* res = NULL;

    if (degPolynomial_spX (r) >= degPolynomial_spX (b)) {
      return 0;
    }

    plainMulPolynomials_spX (b, q, &right, Pptr);
    addPolynomials_spX (right, r, &res, Pptr);
    freePolynomial_spX (&right);

    return isEqual_spX (res, a);
}

int ExtEuclideanTest (const duspoly_t* a, const duspoly_t* b, const duspoly_t* s, const duspoly_t* t, const duspoly_t* g, const Prime_ptr* Pptr) {

    duspoly_t* sa = NULL;
    duspoly_t* tb = NULL;
    duspoly_t* res = NULL;

    // g = s*a + t*b
    plainMulPolynomials_spX (s, a, &sa, Pptr);
    plainMulPolynomials_spX (t, b, &tb, Pptr);
    addPolynomials_spX (sa, tb, &res, Pptr);
    freePolynomial_spX (&sa);
    freePolynomial_spX (&tb);

    return isEqual_spX (res, g);
}

int ExtEuclideanTestInForm (const duspoly_t* a, const duspoly_t* b, const duspoly_t* s, const duspoly_t* t, const duspoly_t* g, const Prime_ptr* Pptr) {

    duspoly_t* sa = NULL;
    duspoly_t* tb = NULL;
    duspoly_t* res = NULL;

    // g = s*a + t*b
    plainMulPolynomialsInForm_spX (s, a, &sa, Pptr);
    plainMulPolynomialsInForm_spX (t, b, &tb, Pptr);
    addPolynomialsInForm_spX (sa, tb, &res, Pptr);
    freePolynomial_spX (&sa);
    freePolynomial_spX (&tb);

    /* fprintf (stderr, "a*s + b*t = "); */
    /* printPolynomial_spX (res, Pptr); */


    return isEqual_spX (res, g);
}

int SFFVerificationInForm_spX (const duspoly_t* a, factors_t* f, const Prime_ptr* Pptr)
{
    if (isZero_spX (a)) {
        if (isZeroFactors_spX (f)) {
            return 1;
        } else {
            return 0;
        }
    }

    duspoly_t* aa = f->polys[0];
    duspoly_t* t1;
    duspoly_t* t2;

    for (polysize_t i = 1; i < f->alloc; i++) {
        if (!isZero_spX (f->polys[i])) {

            // fprintf(stderr, "i = %lu\n", i);

            exponentiatePolynomialInForm_spX (f->polys[i], i, &t1, Pptr);
            plainMulPolynomialsInForm_spX (aa, t1, &t2, Pptr);

            freePolynomial_spX (&t1);
            freePolynomial_spX (&aa);

            aa = t2; t2 = NULL;
        }
    }

    if (!isEqual_spX (a, aa)) {

        fprintf(stderr, "a = ");
        printPolynomial_spX (a, Pptr);

        fprintf(stderr, "f[0]*...*f[n]) = \n");
        printPolynomial_spX (aa, Pptr);


        fprintf(stderr, "DUSP Error, squareFreeFactorizationInForm is not correct!\n");
        exit(1);
    }

    return 1;
}
