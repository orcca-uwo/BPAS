
#include "PowerSeries/PowerSeries.h"
#include "PowerSeries/UnivariatePolynomialOverPowerSeries.h"

#if PARALLEL_POWER_SERIES
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#endif

/**
 * Create a power series struct with some default allocation size.
 * @param alloc: the size of the array of polys to allocate.
 * @return the newly created power series
 */
PowerSeries_t* allocatePowerSeries_PS(int alloc)  {
    PowerSeries_t* ps = (PowerSeries_t*) malloc(sizeof(PowerSeries_t));
    ps->refCount = 1;

    if (alloc > 0) {
        ps->polys = (Poly_ptr*) malloc(sizeof(Poly_ptr)*alloc);
        ps->alloc = alloc;
    } else {
        ps->polys = NULL;
        ps->alloc = 0;
    }

    ps->deg = -1;

    ps->genOrder = -1;
    ps->gen.nullaryGen = NULL;

    ps->genParam1 = NULL;
    ps->genParam2 = NULL;
    ps->genParam3 = NULL;

    ps->paramType1 = PLAIN_DATA;
    ps->paramType2 = PLAIN_DATA;
    ps->paramType3 = PLAIN_DATA;

    return ps;
}

void resizePowerSeries_PS(PowerSeries_t* f, int newAlloc) {
    if (f == NULL) {
        return;
    }

    if (newAlloc < f->deg) {
        int deg = f->deg;
        for (int i = newAlloc + 1; i <= deg; ++i) {
            freePolynomial_AA(f->polys[i]);
        }
        f->polys = realloc(f->polys, newAlloc*sizeof(Poly_ptr));
    }

    if (newAlloc > f->alloc) {
        f->polys = realloc(f->polys, newAlloc*sizeof(Poly_ptr));
    }

    if (f->polys == NULL) {
        fprintf(stderr, "BPAS Error: Could not allocate memory in resizePowerSeries_PS.\n");
        exit(1);
    }

    f->alloc = newAlloc;
}



/**
 * Destroy the power series.
 * Actually, decerements the reference count and destroys conditionally.
 * @param ps : a power series
 */
void destroyPowerSeries_PS(PowerSeries_t* ps) {

    --(ps->refCount);

    if (ps->refCount <= 0) {
        if(ps->genParam1 != NULL) {
            if (ps->paramType1 == POWER_SERIES) {
                destroyPowerSeries_PS((PowerSeries_t*) ps->genParam1);
            } else if (ps->paramType1 == UPOPS) {
                destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) ps->genParam1);
            } else if (ps->paramType1 == MPQ_LIST) {
                destroyMPQList_PS((mpq_list_t*) ps->genParam1);
            }
        }
        if(ps->genParam2 != NULL) {
            if (ps->paramType2 == POWER_SERIES) {
                destroyPowerSeries_PS((PowerSeries_t*) ps->genParam2);
            } else if (ps->paramType2 == UPOPS) {
                destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) ps->genParam2);
            } else if (ps->paramType2 == MPQ_LIST) {
                destroyMPQList_PS((mpq_list_t*) ps->genParam2);
            }
        }
        if(ps->genParam3 != NULL) {
            if (ps->paramType3 == POWER_SERIES) {
                destroyPowerSeries_PS((PowerSeries_t*) ps->genParam3);
            } else if (ps->paramType3 == UPOPS) {
                destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) ps->genParam3);
            } else if (ps->paramType3 == MPQ_LIST) {
                destroyMPQList_PS((mpq_list_t*) ps->genParam3);
            }
        }

        for (int i = 0; i <= ps->deg; ++i) {
             freePolynomial_AA(ps->polys[i]);
        }
        free(ps->polys);
        free(ps);
    }
}


/**
 * Increment the refrence count of the input power series
 * @param ps : a power series
 */
void reserve_PS(PowerSeries_t* ps) {
    ++(ps->refCount);
}


/**
 * Create the zero power series.
 * @return a pointer to the zero power series.
 */
PowerSeries_t* zeroPowerSeries_PS() {
    PowerSeries_t* zero = allocatePowerSeries_PS(1);
    zero->deg = -1;
    zero->polys[0] = NULL;
    zero->genOrder = -1;
    zero->nvar = -1;
    return zero;
}


/**
 * Create the one power series.
 * @return a pointer to the one power series.
 */
PowerSeries_t* onePowerSeries_PS(int nvar) {
    PowerSeries_t* One = allocatePowerSeries_PS(1);
    One->deg = 0;
    One->polys[0] = makeConstIntPolynomial_AA(1,nvar,1, 1);
    One->nvar = nvar;
    One->genOrder = -1;

    return One;
}


/**
 * Create a constant power series from the constnat coef.
 * @param coef : a constant
 @return a pointer to the constant power series.
 */
PowerSeries_t* constPowerSeries_PS(const mpq_t coef, int nvar) {

    PowerSeries_t* constant = allocatePowerSeries_PS(1);
    constant->deg = 0;
    constant->polys[0] = makeConstPolynomial_AA(1,nvar,coef);
    constant->nvar = nvar;
    constant->genOrder = -1;

    return constant;
}


/**
 * d : the degree
 * nvar : the number of variables
 * compute the homogeneous part of a geometric series
 */
Poly_ptr homogPart_geo_series_PS(int d, long long nvar) {
    int partialDegs[nvar];
    for (int i = 0; i < nvar; ++i) {
        partialDegs[i] = 0;
    }
    Poly_ptr monomial = makeConstIntPolynomial_AA(1,nvar, 1, 1);

    Poly_ptr sum = NULL;
    for (int i = 0; i < nvar; ++i) {
        partialDegs[i] = 1;
        setDegrees_AA_inp(monomial, 0, partialDegs, nvar);
        partialDegs[i] = 0;
        sum = addPolynomials_AA_inp(sum, monomial, nvar);
    }

    Poly_ptr ret = exponentiatePoly_AA(sum, d, nvar);

    freePolynomial_AA(monomial);
    freePolynomial_AA(sum);

    return ret;
}


/**
 * d : the degree
 * nvar : the number of variables
 */
Poly_ptr homogPartVoid_geo_series_PS(int d, void* nvar) {
    return homogPart_geo_series_PS(d, (long long) nvar);
}


/**
 * Create the geometric series of nvar number of variables as a power series.
 * @param nvar : the number of variables
 * @return the geometric series a power series.
 */
PowerSeries_t* geometricSeries_PS(long long nvar) {
    PowerSeries_t* geoSeries = allocatePowerSeries_PS(1);
    geoSeries->deg = 0;
    geoSeries->nvar = nvar;
    geoSeries->polys[0] = makeConstIntPolynomial_AA(1, nvar, 1, 1);
    geoSeries->genOrder = 1;
    geoSeries->genParam1 = (void*) nvar;
    geoSeries->paramType1 = PLAIN_DATA;
    geoSeries->gen.unaryGen = &(homogPartVoid_geo_series_PS);
    return geoSeries;
}


/**
 * d : the degree
 * nvar : the number of variables
 * compute the homogeneous part of the sum of all monomials
 */
Poly_ptr homogPart_sumAll_series_PS(int d, long long nvar) {
    int partialDegs[nvar];
    for (int i = 0; i < nvar; ++i) {
        partialDegs[i] = 0;
    }
    Poly_ptr monomial = makeConstIntPolynomial_AA(1,nvar, 1, 1);

    Poly_ptr sum = NULL;
    for (int i = 0; i < nvar; ++i) {
        partialDegs[i] = 1;
        setDegrees_AA_inp(monomial, 0, partialDegs, nvar);
        partialDegs[i] = 0;
        sum = addPolynomials_AA_inp(sum, monomial, nvar);
    }

    Poly_ptr ret = exponentiatePoly_AA(sum, d, nvar);
    for (int i = 0; i < ret->size; ++i) {
        mpz_set_ui(mpq_numref(ret->elems[i].coef), 1ul);
        mpz_set_ui(mpq_denref(ret->elems[i].coef), 1ul);
    }

    freePolynomial_AA(monomial);
    freePolynomial_AA(sum);

    return ret;
}


/**
 * d : the degree
 * nvar : the number of variables
 */
Poly_ptr homogPartVoid_sumAll_series_PS(int d, void* nvar) {
    return homogPart_geo_series_PS(d, (long long) nvar);
}


/**
 * Create the geometric series of nvar number of variables as a power series.
 * @param nvar : the number of variables
 * @return the geometric series a power series.
 */
PowerSeries_t* sumOfAllMonomials_PS(long long nvar) {
    PowerSeries_t* sumSeries = allocatePowerSeries_PS(1);
    sumSeries->deg = 0;
    sumSeries->nvar = nvar;
    sumSeries->polys[0] = makeConstIntPolynomial_AA(1, nvar, 1, 1);
    sumSeries->genOrder = 1;
    sumSeries->genParam1 = (void*) nvar;
    sumSeries->paramType1 = PLAIN_DATA;
    sumSeries->gen.unaryGen = &(homogPartVoid_sumAll_series_PS);
    return sumSeries;
}



/**
 * Print a power series to the file pointer fp
 * using the symbols sym as symbols of the homogeneous polynomial.
 * @param fp: the file pointer to print to; may be stdout or stderr or something else.
 * @param ps : a power series
 * @param sym : a list of variables
 */
void print_PS(FILE* fp, PowerSeries_t* ps, const char** sym) {
    int First = 1;
    for(int i = 0; i <= ps->deg; ++i) {
        if (!isZero_AA(ps->polys[i])) {
            if (!First) {
                if (mpq_sgn(ps->polys[i]->elems->coef) < 0) {
                    mpq_neg(ps->polys[i]->elems->coef, ps->polys[i]->elems->coef);
                    fprintf(fp, " - ");
                    printPolyClean_AA(fp, ps->polys[i], sym, ps->polys[i]->nvar);
                    mpq_neg(ps->polys[i]->elems->coef, ps->polys[i]->elems->coef);
                } else {
                    fprintf(fp, " + ");
                    printPolyClean_AA(fp, ps->polys[i], sym,  ps->polys[i]->nvar);
                }
             } else {
                printPolyClean_AA(fp, ps->polys[i], sym,  ps->polys[i]->nvar);
                First = 0;
             }
        }
   }
}


/**********************
 * Main functional interface
 **********************/

/**
 * Given a power series f and an integer d, return the homogeneous
 * part of f of degree d.
 *
 * The is the main functional interface for power series.
 * Calling this function will generate terms in the power series as needed.
 * @see polynomialPart_PS.
 *
 * @return the homogeneous part of f of degree d.
 */
Poly_ptr homogPart_PS(int d, PowerSeries_t* f) {

    if (d <= f->deg) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
        fprintf(stderr, "Returning from the homogPart and degree: %d\n", d);
        const char* syms[] = {"z", "x", "y"};
        fprintf(stderr, "f: %p, f->polys[%d]: ", f, d );
        printPoly_AA(stderr, f->polys[d], syms, f->polys[d] == NULL ? 0 :f->polys[d]->nvar);
        fprintf(stderr, "\n" );
#endif
        return f->polys[d];

    } else {
        if (f->genOrder == -1) {
            // fprintf(stderr, "Returning NULL from homogPart_PS with d=%d, f=%p\n", d, f);
            return NULL;
        }

        if (d + 1 > f->alloc) {

            int newAlloc = (2*f->alloc < d+1) ? d + 1 : 2*f->alloc;
            f->polys = (Poly_ptr*) realloc(f->polys, sizeof(Poly_ptr)*newAlloc);
            f->alloc = newAlloc;
            if (f == NULL) {
                fprintf(stderr, "Reallocation Error!\n");
                exit(1);
            }
        }
        for (int i = f->deg+1; i <= d; ++i ) {
            switch (f->genOrder) {
                case -1 : {
                    return NULL;
                }
                case 0: {
                    f->polys[i] = (*f->gen.nullaryGen)(i);
                    break;
                }
                case 1: {
                    f->polys[i] = (*f->gen.unaryGen)(i, f->genParam1);
                    break;
                }
                case 2: {
                    f->polys[i] = (*f->gen.binaryGen)(i, f->genParam1, f->genParam2);
                    break;
                }
                case 3: {
                    f->polys[i] = (*f->gen.tertiaryGen)(i, f->genParam1, f->genParam2, f->genParam3);
                    break;
                }
            }

            f->deg = i;
        }
        return f->polys[d];
    }

}

/**
 * Computes the polynomial part of a power series up to degree d
 * @param d: the degree of the resulting polynomial part
 * @param f: a power series
 *
 * @return the polynomial part of f of degree d.
 */
Poly_ptr polynomialPart_PS(int d,   PowerSeries_t* f) {
    Poly_ptr p = NULL;
    int i;
    const char* syms[] = {"z", "x", "y"};
    for (i = 0; i <= d; i++) {
        // fprintf(stderr, "adding polys in polynomialPart_PS %d\n P: ", i);
        // printPolyClean_AA(stderr, p, syms, p == NULL ? 0 : p->nvar);

        // fprintf(stderr, "\nTmp: ");
        Poly_ptr tmp = homogPart_PS(i, f);
        // printPolyClean_AA(stderr, tmp, syms, tmp == NULL ? 0 : tmp->nvar);
        // fprintf(stderr, "\n");

        p = addPolynomials_AA_inp(p, tmp, tmp == NULL ? 0 : tmp->nvar);
        // fprintf(stderr, "\nresulting poly in polynomialPart_PS %d\n ", i);
        // printPolyClean_AA(stderr, p, syms, p == NULL ? 0 : p->nvar);
        // fprintf(stderr, "\n");
    }

    return p;
}



/**********************
 * Conversion Helpers
 **********************/

/**
 * Given a polynomial p, it returns an array of homogeneous polynomials
 * whose sum is equal to p.
 * In the returned array, the polynomial of index i has degree i.
 * @return an array of homogeneous polynomials
 */
Poly_ptr* convertPolyToArrayOfHomogeneous_PS(const Poly_ptr p) {

    int i;
    int array_size = totalDegree_AA(p);

    Poly_ptr* arrayOfHomogeneous = (Poly_ptr*) malloc(sizeof(Poly_ptr)*((array_size)+1));
    for (i = 0; i <= totalDegree_AA(p); ++i) {
        arrayOfHomogeneous[i] = homogeneousPart_AA(p,i);
    }

    return arrayOfHomogeneous;
}


/**
 * Eliminates terms of degree atleast d+1 from polynomial p.
 * @param p: the polynomial to truncate
 * @param d: the (maximum) degree of the resulting truncated polynomial.
 * @return the truncated polynomial
 */
Poly_ptr truncatePolynomial_PS(const Poly_ptr p, int d) {
    int j;
    Poly_ptr polytrunc = NULL;
    int array_size = p->size;
    Poly_ptr term;
    for (j = 0; j<=array_size; j++) {
        term = termAtIdx_AA(p, j);
        if (totalDegree_AA(term) <= d) {
            polytrunc = addPolynomials_AA_inp(polytrunc, term, p->nvar);
        }
        freePolynomial_AA(term);
    }
    return polytrunc;
}

/**
 * Converts a polynomial to a power series.
 *
 * @param p: the polynomial to convert
 * @return the power series encoding the input polynomial.
 */
PowerSeries_t* convertPolyToPowerSeries_PS(const Poly_ptr p) {
    if (isZero_AA(p)) {
        return zeroPowerSeries_PS();
    }
    degree_t totalDeg = totalDegree_AA(p);
    PowerSeries_t* ps = allocatePowerSeries_PS(0);
    ps->nvar = p == NULL ? 0 : p->nvar;
    ps->deg = totalDeg;
    ps->polys = convertPolyToArrayOfHomogeneous_PS(p);
    ps->alloc = totalDeg +1;
    ps->genOrder = -1;

    return ps;
}




/**********************
 * Addition and Subtraction
 **********************/


/**
 * An internal function.
 * Computes homogeneous part of the sum of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPart_sum_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g) {
    const char* syms[] = {"z", "x", "y"};
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "Entering homogPart_sum_PS: %d with %p, %p\n",d, f, g);
    fprintf(stderr, "f: ");
    print_PS(stderr, f, syms);
    fprintf(stderr, "\ng: ");
    print_PS(stderr, g, syms);
    fprintf(stderr, "\n");

#endif
    Poly_ptr lefthomog = homogPart_PS(d,f);
    // fprintf(stderr, "got left homog\n");
    // printPoly_AA(stderr, lefthomog, syms, lefthomog == NULL ? 0 : lefthomog->nvar);
    // fprintf(stderr, "\n" );

    Poly_ptr righthomog = homogPart_PS(d, g);
    // fprintf(stderr, "got right homog\n");
    Poly_ptr result = addPolynomials_AA(lefthomog, righthomog, lefthomog == NULL ? 0 : lefthomog->nvar);


#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "Leaving  homogPart_sum_PS: %d with %p, %p\n",d, f, g);
#endif
    return result;
}


/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series sum.
 * @param d: the requested degree to generate
 * @param param1: the left-hand side of the addition.
 * @param param2: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPartVoid_sum_PS(int deg, void* param1, void* param2) {
        return homogPart_sum_PS(deg, (PowerSeries_t*) param1, (PowerSeries_t*) param2);
}


/**
 * Given two power series f and g, it returns the sum f+g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeries_t* addPowerSeries_PS(PowerSeries_t* f,  PowerSeries_t* g) {
    if (!isZero_PS(f) && !isZero_PS(g) && f->nvar != g->nvar) {
        fprintf(stderr, "BPAS Error: f and g must have same number of variables in addPowerSeries_PS\n");
        exit(1);
    }


    int i;
    int mindeg;
    int nvar;
    if (g->deg == -1 && f->deg == -1) {
      //  mindeg = 0;
        return zeroPowerSeries_PS();
    } else if (g->deg == -1) {
        mindeg = f->deg;
        nvar = f->nvar;
      //  reserve_PS(f); return f;Example 16 in weierstrass crashes
    } else if (f->deg == -1){
        mindeg = g->deg;
        nvar = g->nvar;
    //  reserve_PS(g); return g; Example 16 in weierstrass crashes
    } else {
        //Half-lazily:
        // mindeg = MIN(f->deg, g->deg);
        mindeg = 0;
        nvar = f->nvar;
    }


    PowerSeries_t* sum = allocatePowerSeries_PS(mindeg+1);
    sum->nvar = nvar;
    sum->genParam1 = f;
    sum->genParam2 = g;
    sum->paramType1 = POWER_SERIES;
    sum->paramType2 = POWER_SERIES;

    reserve_PS(f);
    reserve_PS(g);

    sum->genOrder = 2;
    sum->gen.binaryGen = &(homogPartVoid_sum_PS);
    sum->deg = mindeg;
    for (i = 0; i <= mindeg; i++){
        sum->polys[i] = homogPart_sum_PS(i, f, g);
    }

    return sum;
}


/**
 * An internal function.
 * Computes homogeneous part of the sum of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPart_sub_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g) {

    Poly_ptr lefthomog = homogPart_PS(d,f);
    Poly_ptr righthomog = homogPart_PS(d, g);
    Poly_ptr result;
    if (isZero_AA(lefthomog) && isZero_AA(righthomog)){
         result = NULL;
    } else {
         result = subPolynomials_AA(lefthomog, righthomog, lefthomog == NULL ? 0 : lefthomog->nvar);
    }


    return result;
}


/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series difference.
 * @param d: the requested degree to generate
 * @param param1: the left-hand side of the difference.
 * @param param2: the right-hand side of the difference
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPartVoid_sub_PS(int deg, void* param1, void* param2) {
    return homogPart_sub_PS(deg, (PowerSeries_t*) param1, (PowerSeries_t*) param2);
}


/**
 * Given two power series f and g, it returns the difference f-g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeries_t* subPowerSeries_PS(PowerSeries_t* f,  PowerSeries_t* g) {
    if (!isZero_PS(f) && !isZero_PS(g) && f->nvar != g->nvar) {
        fprintf(stderr, "BPAS Error: f and g must have same number of variables in subPowerSeries_PS\n");
        free((int*)-1);
        exit(1);
    }

    int i;
    int mindeg;
    int nvar;
    if (g->deg == -1 && f->deg == -1) {
        return zeroPowerSeries_PS();
    //mindeg = 0;
    } else if (g->deg == -1) {
        //reserve_PS(f);
        //return f;
        mindeg = f->deg;//negatePowerSeries_PS
        nvar = f->nvar;
    } else if (f->deg == -1){
        //reserve_PS(g);return negatePowerSeries_PS(g);
        mindeg = g->deg;
        nvar = g->nvar;
    } else {
        // mindeg = MIN(f->deg, g->deg);
        mindeg = 0;
        nvar = f->nvar;
    }
    PowerSeries_t* sub = allocatePowerSeries_PS(mindeg+1);
    sub->nvar = nvar;
    sub->genParam1 = f;
    sub->genParam2 = g;
    sub->paramType1 = POWER_SERIES;
    sub->paramType2 = POWER_SERIES;

    reserve_PS(f);
    reserve_PS(g);
    sub->genOrder = 2;
    sub->gen.binaryGen = &(homogPartVoid_sub_PS);
    sub->deg = mindeg;
    for (i = 0; i <= mindeg; i++){
        sub->polys[i] = homogPart_sub_PS(i, f, g);
    }

    return sub;
}


/**
 * An internal function.
 * Generate the homogeneous part of degree d of the
 * negation of power series f.
 * d : the degree of the homogeneous part to generate
 * f : the power series to negate
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPart_negate_PS(int d, PowerSeries_t* f)  {
    Poly_ptr poly = deepCopyPolynomial_AA(homogPart_PS(d, f));
    negatePolynomial_AA(poly);
    return poly;
}


/**
 * An internal function.
 * A void generator wrapper for the negation of a power series.
 *
 * @param d: the degree of the homogeneous part to generate.
 * @param param1: the power series to negate.
 */
Poly_ptr homogPartVoid_negate_PS(int d, void* param1) {
    return homogPart_negate_PS(d, (PowerSeries_t*) param1);
}


/**
 * Negate the power series f, returning a lazily constructed power series.
 * @param f: the power series to negate.
 * @return the negation of the input power series.
 */
PowerSeries_t* negatePowerSeries_PS(PowerSeries_t* f) {

    if (isZero_PS(f)) {
        reserve_PS(f);
        return f;
    }

    PowerSeries_t* ps = allocatePowerSeries_PS(f->deg+1);
    ps->deg = 0;
    ps->nvar = f->nvar;
    ps->genOrder = 1;
    ps->genParam1 = f;
    reserve_PS(f);
    ps->paramType1 = POWER_SERIES;

    ps->gen.unaryGen = &(homogPartVoid_negate_PS);
    for (int i = 0; i <= ps->deg; i++){

        if (f->polys[i] == NULL) {
            ps->polys[i] = NULL;
        } else {
            ps->polys[i] = homogPart_negate_PS(i, f);
        }

    }
    return ps;

}



/**********************
 * Multiplication and Division
 **********************/

/**
 * An internal function.
 * Computes homogeneous part of the product of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand operand of the multiplication
 * @param g: the right-hand operand of the multiplication
 * @return the homogeneous part of degree d of the product.
 */
Poly_ptr homogPart_prod_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "Entering homogPart_prod_PS: %d with %p, %p\n",d, f, g);
#endif
    Poly_ptr prod;
    Poly_ptr sum = NULL;

    // const char* syms[] = {"z", "x", "y"};
    Poly_ptr tmp1, tmp2;
    for (int i = 0; i <= d; i++) {
        tmp1 = homogPart_PS(d-i, f);
        if (isZero_AA(tmp1)) {
            continue;
        }
        tmp2 = homogPart_PS(i, g);
        if (isZero_AA(tmp2)) {
            continue;
        }
        // fprintf(stderr, "tmp1: ");
        // printPoly_AA(stderr, tmp1, syms, tmp1->nvar);
        // fprintf(stderr, "\ntmp2: ");
        // printPoly_AA(stderr, tmp2, syms, tmp2->nvar);
        // fprintf(stderr, "\n");
        prod = multiplyPolynomials_AA(tmp1, tmp2, tmp1->nvar);


        sum = addPolynomials_AA_inp(sum, prod, prod->nvar); //prod is guaranteed to not be 0 by above checks.
        freePolynomial_AA(prod);
    }
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "\n" );
    AltArr_t localSum;
    localSum.nvar = -12;
    if (sum != NULL) {
        localSum = *sum;
    }
    printPoly_AA(stderr, sum, syms, sum == NULL? 0 : sum->nvar);
    fprintf(stderr, "\n");
    fprintf(stderr, "Leaving homogPart_prod_PS: %d with %p, %p\n",d, f, g);
#endif
    return sum;
}


/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series product.
 * @param d: the requested degree to generate
 * @param param1: the left-hand operand of the multiplication
 * @param param2: the right-hand operand of the multiplication
 * @return the homogeneous part of degree d of the product.
 */
Poly_ptr homogPartVoid_prod_PS(int d, void* param1, void* param2) {
    return homogPart_prod_PS(d, (PowerSeries_t*) param1, (PowerSeries_t*) param2);
}

/**
 * Given two power series f and g, it returns the product, constructed lazily.
 * @param f: the multiplier power series
 * @param g: the multiplicand power series
 * @return a pointer to the resulting power series product.
 */
PowerSeries_t* multiplyPowerSeries_PS( PowerSeries_t* f,  PowerSeries_t* g) {
    if (!isZero_PS(f) && !isZero_PS(g) && f->nvar != g->nvar) {
        fprintf(stderr, "BPAS Error: f and g must have same number of variables in multiplyPowerSeries_PS\n");
        free((int*)-1);
        exit(1);
    }

    int i;
    int mindeg;


    if (g->deg == -1 && f->deg == -1) {
        //mindeg = 0;
        return zeroPowerSeries_PS();
    } else if (g->deg == -1) {
       // mindeg = 0;//f->deg;
        return zeroPowerSeries_PS();
    } else if (f->deg == -1){
       // mindeg = 0;//g->deg;
        return zeroPowerSeries_PS();
    } else {
        // mindeg = MIN(f->deg, g->deg);
        mindeg = 0;
    }

    PowerSeries_t* prod = allocatePowerSeries_PS(mindeg+1);
    prod->nvar = f->nvar;
    prod->genParam1 = f;
    prod->genParam2 = g;
    prod->paramType1 = POWER_SERIES;
    prod->paramType2 = POWER_SERIES;

    reserve_PS(f);
    reserve_PS(g);
    prod->genOrder = 2;
    prod->gen.binaryGen = &(homogPartVoid_prod_PS);
    prod->deg = mindeg;
    for (i = 0; i <= (mindeg); i++){
        prod->polys[i] = homogPart_prod_PS(i, f, g);
    }

    return prod;
}



#if PARALLEL_POWER_SERIES
/**
 * Given two power series f and g, it returns the product, constructed lazily.
 * @param f: the multiplier power series
 * @param g: the multiplicand power series
 * @return a pointer to the resulting power series product.
 */
PowerSeries_t* multiplyPowerSeries_PS( PowerSeries_t* f,  PowerSeries_t* g) {

    int i;
    int mindeg;
     /*Later I will remove extra if statements*/
    if (g->deg == -1 && f->deg == -1) {
        //mindeg = 0;
        return zeroPowerSeries_PS();
    } else if (g->deg == -1) {
       // mindeg = 0;//f->deg;
        return zeroPowerSeries_PS();
    } else if (f->deg == -1){
       // mindeg = 0;//g->deg;
        return zeroPowerSeries_PS();
    } else {
        mindeg = MIN(f->deg, g->deg);

    }

    PowerSeries_t* prod = allocatePowerSeries_PS(mindeg+1);
    prod->genParam1 = f;
    prod->genParam2 = g;
    prod->paramType1 = POWER_SERIES;
    prod->paramType2 = POWER_SERIES;

    reserve_PS(f);
    reserve_PS(g);
    prod->genOrder = 2;
    prod->gen.binaryGen = &(homogPartVoid_prod_PS);
    prod->deg = mindeg;
    for (i = 0; i <= (mindeg); i++){
        prod->polys[i] = homogPart_prod_PS(i, f, g);
    }

    return prod;
}
#endif











/**
 * An internal function.
 * Computes homogeneous part of the quotient of two power series of degree d.
 * @param deg: the requested degree to generate
 * @param f: the dividend
 * @param h: the divisor, which must be a unit.
 * @param quo, the quotient being generated.
 * @return the homogeneous part of degree d of the quotient.
 */
Poly_ptr homogPart_quo_PS(int deg,  PowerSeries_t* f,  PowerSeries_t* h, PowerSeries_t* quo) {

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
   fprintf(stderr, "requested degree quo: %d\n", deg);
#endif

    Poly_ptr gd;
    int i;
    Poly_ptr h0 = homogPart_PS(0, h);
    if (h0 == 0){
        fprintf(stderr, "ERROR: In homogPart_quo_PS, the divisor was zero!\n");
        exit(1);
    }

    if (deg <= 1) {

        Poly_ptr q0 = deepCopyPolynomial_AA(homogPart_PS(0,f));
        divideByRational_AA_inp(q0, h0->elems->coef);
        if (deg == 0){
            return q0;
        }

        Poly_ptr h1 = homogPart_PS(1,h);
        Poly_ptr f1 = homogPart_PS(1,f);
        Poly_ptr multipoly = multiplyPolynomials_AA(q0, h1, q0 == NULL ? 0 : q0->nvar);
        freePolynomial_AA(q0);
        multipoly = subPolynomials_AA_inp( multipoly, f1,  f1 == NULL ? 0 : f1->nvar);
        negatePolynomial_AA(multipoly);
        divideByRational_AA_inp(multipoly, h0->elems->coef);
        return multipoly;
    }

    Poly_ptr s = deepCopyPolynomial_AA(homogPart_PS(deg,f));

    for (i = 1; i <= deg; ++i) {
        Poly_ptr homog = homogPart_PS(i,h);
        //Poly_ptr homogquo = homogPart_quo_PS(deg-i, f, h);

        Poly_ptr homogquo = homogPart_PS(deg-i, quo);
        Poly_ptr p = multiplyPolynomials_AA(homog, homogquo, homog == NULL ? 0 : homog->nvar);
        //freePolynomial_AA(homogquo);
        s = subPolynomials_AA_inp(s, p, s == NULL ? 0 : s->nvar);
        freePolynomial_AA(p);
    }

    divideByRational_AA_inp(s, h0->elems->coef);
    return s;
}


/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series quotient.
 * @param d: the requested degree to generate
 * @param param1: the dividend
 * @param param2: the divisor
 * @param param3: the quotient itself
 * @return the homogeneous part of degree d of the quotient.
 */
Poly_ptr homogPartVoid_quo_PS(int d, void* param1, void* param2, void* param3) {
  return homogPart_quo_PS(d, (PowerSeries_t*) param1, (PowerSeries_t*) param2, (PowerSeries_t*) param3);
}


/**
 * Given two power series f and h, it returns the quotient f/h, constructed lazily.
 * @param f: the dividend power series
 * @param h: the divisor power series
 * @return a pointer to the resulting power series quotient.
 */
PowerSeries_t* dividePowerSeries_PS(PowerSeries_t* f,  PowerSeries_t* h)  {
    if (!isUnit_PS(h)) {
        fprintf(stderr, "BPAS Error: divisor is not inverible in dividePowerSeries_PS\n" );
        exit(1);
    }

    if (isZero_PS(f)) {
        reserve_PS(f);
        return f;
    }

    if (!isZero_PS(f) && f->nvar != h->nvar) {
        fprintf(stderr, "BPAS Error: f and h must have same number of variables in dividePowerSeries_PS\n");
        exit(1);
    }


    PowerSeries_t* quo = allocatePowerSeries_PS(1);
    quo->nvar = f->nvar;
    quo->genParam1 = f;
    quo->genParam2 = h;
    quo->genParam3 = quo;

    quo->paramType1 = POWER_SERIES;
    quo->paramType2 = POWER_SERIES;
    quo->paramType3 = PLAIN_DATA; //Yes plain_data. Avoids chicken-and-egg with circular reference.

    reserve_PS(f);
    reserve_PS(h);
    quo->genOrder = 3;
    quo->gen.tertiaryGen = &(homogPartVoid_quo_PS);
    quo->deg = 0;
    //quo->polys[0] = homogPart_quo_PS(0, f, h);WAS
    quo->polys[0] = homogPart_quo_PS(0, f, h, quo);
    return quo;
}


