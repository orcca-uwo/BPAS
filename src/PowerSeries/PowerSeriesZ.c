
#include "PowerSeries/PowerSeriesZ.h"

/**
 * Create a power series struct with some default allocation size.
 * @param alloc: the size of the array of polys to allocate.
 * @return the newly created power series
 */
PowerSeriesZ_t* allocatePowerSeries_PSZ(int alloc)  {
    PowerSeriesZ_t* ps = (PowerSeriesZ_t*) malloc(sizeof(PowerSeriesZ_t));
    ps->refCount = 1;

    if (alloc > 0) {
        ps->polys = (PolyZ_ptr*) malloc(sizeof(PolyZ_ptr)*alloc);
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



/**
 * Destroy the power series.
 * Actually, decerements the reference count and destroys conditionally.
 * @param ps : a power series
 */
void destroyPowerSeries_PSZ(PowerSeriesZ_t* ps) {

    --(ps->refCount);

    if (ps->refCount <= 0) {
        if(ps->genParam1 != NULL) {
            if (ps->paramType1 == POWER_SERIES) {
                destroyPowerSeries_PSZ((PowerSeriesZ_t*) ps->genParam1);
            // } else if (ps->paramType1 == UPOPS) {
                // destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) ps->genParam1);
            // } else if (ps->paramType1 == MPQ_LIST) {
                // destroyMPQList_PSZ((mpq_list_t*) ps->genParam1);
            }
        }
        if(ps->genParam2 != NULL) {
            if (ps->paramType2 == POWER_SERIES) {
                destroyPowerSeries_PSZ((PowerSeriesZ_t*) ps->genParam2);
            // } else if (ps->paramType2 == UPOPS) {
                // destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) ps->genParam2);
            // } else if (ps->paramType2 == MPQ_LIST) {
                // destroyMPQList_PSZ((mpq_list_t*) ps->genParam2);
            }
        }
        if(ps->genParam3 != NULL) {
            if (ps->paramType3 == POWER_SERIES) {
                destroyPowerSeries_PSZ((PowerSeriesZ_t*) ps->genParam3);
            // } else if (ps->paramType3 == UPOPS) {
                // destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) ps->genParam3);
            // } else if (ps->paramType3 == MPQ_LIST) {
                // destroyMPQList_PSZ((mpq_list_t*) ps->genParam3);
            }
        }

        for (int i = 0; i <= ps->deg; ++i) {
             freePolynomial_AAZ(ps->polys[i]);
        }
        free(ps->polys);
        free(ps);
    }
}


/**
 * Increment the refrence count of the input power series
 * @param ps : a power series
 */
void reserve_PSZ(PowerSeriesZ_t* ps) {
    ++(ps->refCount);
}


/**
 * Create the zero power series.
 * @return a pointer to the zero power series.
 */
PowerSeriesZ_t* zeroPowerSeries_PSZ() {
    PowerSeriesZ_t* zero = allocatePowerSeries_PSZ(1);
    zero->deg = -1;
    zero->polys[0] = NULL;
    zero->genOrder = -1;

    return zero;
}


/**
 * Create the one power series.
 * @return a pointer to the one power series.
 */
PowerSeriesZ_t* onePowerSeries_PSZ(int nvar) {
    PowerSeriesZ_t* One = allocatePowerSeries_PSZ(1);
    One->deg = 0;
    One->polys[0] = makeConstIntPolynomial_AAZ(1,nvar,1);
    One->genOrder = -1;

    return One;
}


/**
 * Create a constant power series from the constnat coef.
 * @param coef : a constant
 @return a pointer to the constant power series.
 */
PowerSeriesZ_t* constPowerSeries_PSZ(const mpz_t coef, int nvar) {

    PowerSeriesZ_t* constant = allocatePowerSeries_PSZ(1);
    constant->deg = 0;
    constant->polys[0] = makeConstPolynomial_AAZ(1,nvar,coef);
    constant->genOrder = -1;

    return constant;
}


/**
 * d : the degree
 * nvar : the number of variables
 * compute the homogeneous part of a geometric series
 */
PolyZ_ptr homogPart_geo_series_PSZ(int d, long long nvar) {
   int partialDegs[nvar];
   for (int i = 0; i <= nvar-1; ++i) {
        partialDegs[i] = 0;
    }
   PolyZ_ptr monomial = makeConstIntPolynomial_AAZ(1,nvar, 1);

   PolyZ_ptr sum = NULL;
   for (int i = 0; i <= nvar-1; ++i) {
       partialDegs[i] = 1;
       setDegrees_AAZ_inp(monomial, 0, partialDegs, nvar);
       partialDegs[i] = 0;
       sum = addPolynomials_AAZ_inp(sum, monomial, nvar);
   }


   PolyZ_ptr ret = makeConstIntPolynomial_AAZ(1, nvar, 1);
   for (int j = 1; j <= d; ++j) {
       ret = multiplyPolynomials_AAZ_inp(ret,sum ,nvar);
   }

   freePolynomial_AAZ(monomial);
   freePolynomial_AAZ(sum);

   return ret;
}


/**
 * d : the degree
 * nvar : the number of variables
 */
PolyZ_ptr homogPartVoid_geo_series_PSZ(int d, void* nvar) {
    return homogPart_geo_series_PSZ(d, (long long) nvar);
}


/**
 * Create the geometric series of nvar number of variables as a power series.
 * @param nvar : the number of variables
 * @return the geometric series a power series.
 */
PowerSeriesZ_t* geometricSeries_PSZ(long long nvar) {
    PowerSeriesZ_t* geoSeries = allocatePowerSeries_PSZ(1);
    geoSeries->deg = 0;
    geoSeries->polys[0] = makeConstIntPolynomial_AAZ(1, nvar, 1);
    geoSeries->genOrder = 1;
    geoSeries->genParam1 = (void*) nvar;
    geoSeries->paramType1 = PLAIN_DATA;
    geoSeries->gen.unaryGen = &(homogPartVoid_geo_series_PSZ);
    return geoSeries;
}



/**
 * Print a power series to the file pointer fp
 * using the symbols sym as symbols of the homogeneous polynomial.
 * @param fp: the file pointer to print to; may be stdout or stderr or something else.
 * @param ps : a power series
 * @param sym : a list of variables
 */
void print_PSZ(FILE* fp, PowerSeriesZ_t* ps, const char** sym) {
    int First = 1;
    for(int i = 0; i <= ps->deg; ++i) {
        if (!isZero_AAZ(ps->polys[i])) {
            if (!First) {
                if (mpz_sgn(ps->polys[i]->elems->coef) < 0) {
                    mpz_neg(ps->polys[i]->elems->coef, ps->polys[i]->elems->coef);
                    fprintf(fp, " - ");
                    printPolyClean_AAZ(fp, ps->polys[i], sym, ps->polys[i]->nvar);
                    mpz_neg(ps->polys[i]->elems->coef, ps->polys[i]->elems->coef);
                } else {
                    fprintf(fp, " + ");
                    printPolyClean_AAZ(fp, ps->polys[i], sym,  ps->polys[i]->nvar);
                }
             } else {
                printPolyClean_AAZ(fp, ps->polys[i], sym,  ps->polys[i]->nvar);
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
 * @see polynomialPart_PSZ.
 *
 * @return the homogeneous part of f of degree d.
 */
PolyZ_ptr homogPart_PSZ(int d, PowerSeriesZ_t* f) {

    if (d <= f->deg) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
        fprintf(stderr, "Returning from the homogPart and degree: %d\n", d);
        const char* syms[] = {"z", "x", "y"};
        fprintf(stderr, "f: %p, f->polys[%d]: ", f, d );
        printPoly_AAZ(stderr, f->polys[d], syms, f->polys[d] == NULL ? 0 :f->polys[d]->nvar);
        fprintf(stderr, "\n" );
#endif
        return f->polys[d];

    } else {
        if (f->genOrder == -1) {
            // fprintf(stderr, "Returning NULL from homogPart_PSZ with d=%d, f=%p\n", d, f);
            return NULL;
        }

        if (d + 1 > f->alloc) {

            int newAlloc = (2*f->alloc < d+1) ? d + 1 : 2*f->alloc;
            f->polys = (PolyZ_ptr*) realloc(f->polys, sizeof(PolyZ_ptr)*newAlloc);
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
PolyZ_ptr polynomialPart_PSZ(int d,   PowerSeriesZ_t* f) {
    PolyZ_ptr p = NULL;
    int i;
    const char* syms[] = {"z", "x", "y"};
    for (i = 0; i <= d; i++) {
        // fprintf(stderr, "adding polys in polynomialPart_PSZ %d\n P: ", i);
        // printPolyClean_AA(stderr, p, syms, p == NULL ? 0 : p->nvar);

        // fprintf(stderr, "\nTmp: ");
        PolyZ_ptr tmp = homogPart_PSZ(i, f);
        // printPolyClean_AA(stderr, tmp, syms, tmp == NULL ? 0 : tmp->nvar);
        // fprintf(stderr, "\n");

        p = addPolynomials_AAZ_inp(p, tmp, tmp == NULL ? 0 : tmp->nvar);
        // fprintf(stderr, "\nresulting poly in polynomialPart_PSZ %d\n ", i);
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
PolyZ_ptr* convertPolyToArrayOfHomogeneous_PSZ(const PolyZ_ptr p) {

    int i;
    int array_size = totalDegree_AAZ(p);

    PolyZ_ptr* arrayOfHomogeneous = (PolyZ_ptr*) malloc(sizeof(PolyZ_ptr)*((array_size)+1));
    for (i = 0; i <= totalDegree_AAZ(p); ++i) {
        arrayOfHomogeneous[i] = homogeneousPart_AAZ(p,i);
    }

    return arrayOfHomogeneous;
}


/**
 * Eliminates terms of degree atleast d+1 from polynomial p.
 * @param p: the polynomial to truncate
 * @param d: the (maximum) degree of the resulting truncated polynomial.
 * @return the truncated polynomial
 */
PolyZ_ptr truncatePolynomial_PSZ(const PolyZ_ptr p, int d) {
    int j;
    PolyZ_ptr polytrunc = NULL;
    int array_size = p->size;
    PolyZ_ptr term;
    for (j = 0; j<=array_size; j++) {
        term = termAtIdx_AAZ(p, j);
        if (totalDegree_AAZ(term) <= d) {
            polytrunc = addPolynomials_AAZ_inp(polytrunc, term, p->nvar);
        }
        freePolynomial_AAZ(term);
    }
    return polytrunc;
}

/**
 * Converts a polynomial to a power series.
 *
 * @param p: the polynomial to convert
 * @return the power series encoding the input polynomial.
 */
PowerSeriesZ_t* convertPolyToPowerSeries_PSZ(const PolyZ_ptr p) {
    degree_t totalDeg = totalDegree_AAZ(p);
    PowerSeriesZ_t* ps = allocatePowerSeries_PSZ(0);
    ps->deg = totalDeg;
    ps->polys = convertPolyToArrayOfHomogeneous_PSZ(p);
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
PolyZ_ptr homogPart_sum_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g) {
    const char* syms[] = {"z", "x", "y"};
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "Entering homogPart_sum_PSZ: %d with %p, %p\n",d, f, g);
    fprintf(stderr, "f: ");
    print_PSZ(stderr, f, syms);
    fprintf(stderr, "\ng: ");
    print_PSZ(stderr, g, syms);
    fprintf(stderr, "\n");

#endif
    PolyZ_ptr lefthomog = homogPart_PSZ(d,f);
    // fprintf(stderr, "got left homog\n");
    // printPoly_AA(stderr, lefthomog, syms, lefthomog == NULL ? 0 : lefthomog->nvar);
    // fprintf(stderr, "\n" );

    PolyZ_ptr righthomog = homogPart_PSZ(d, g);
    // fprintf(stderr, "got right homog\n");
    PolyZ_ptr result = addPolynomials_AAZ(lefthomog, righthomog, lefthomog == NULL ? 0 : lefthomog->nvar);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "Leaving  homogPart_sum_PSZ: %d with %p, %p\n",d, f, g);
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
PolyZ_ptr homogPartVoid_sum_PSZ(int deg, void* param1, void* param2) {
        return homogPart_sum_PSZ(deg, (PowerSeriesZ_t*) param1, (PowerSeriesZ_t*) param2);
}


/**
 * Given two power series f and g, it returns the sum f+g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeriesZ_t* addPowerSeries_PSZ(PowerSeriesZ_t* f,  PowerSeriesZ_t* g) {

    int i;
    int mindeg;

    if (g->deg == -1 && f->deg == -1) {
      //  mindeg = 0;
        return zeroPowerSeries_PSZ();
    } else if (g->deg == -1) {
        mindeg = f->deg;
      //  reserve_PSZ(f); return f;Example 16 in weierstrass crashes
    } else if (f->deg == -1){
        mindeg = g->deg;
    //  reserve_PSZ(g); return g; Example 16 in weierstrass crashes
    } else {
        mindeg = MIN(f->deg, g->deg);
    }
    PowerSeriesZ_t* sum = allocatePowerSeries_PSZ(mindeg+1);
    sum->genParam1 = f;
    sum->genParam2 = g;
    sum->paramType1 = POWER_SERIES;
    sum->paramType2 = POWER_SERIES;

    reserve_PSZ(f);
    reserve_PSZ(g);

    sum->genOrder = 2;
    sum->gen.binaryGen = &(homogPartVoid_sum_PSZ);
    sum->deg = mindeg;
    for (i = 0; i <= mindeg; i++){
        sum->polys[i] = homogPart_sum_PSZ(i, f, g);
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
PolyZ_ptr homogPart_sub_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g) {

    PolyZ_ptr lefthomog = homogPart_PSZ(d,f);
    PolyZ_ptr righthomog = homogPart_PSZ(d, g);
    PolyZ_ptr result;
    if (isZero_AAZ(lefthomog) && isZero_AAZ(righthomog)){
         result = NULL;
    } else {
         result = subPolynomials_AAZ(lefthomog, righthomog, lefthomog == NULL ? 0 : lefthomog->nvar);
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
PolyZ_ptr homogPartVoid_sub_PSZ(int deg, void* param1, void* param2) {
    return homogPart_sub_PSZ(deg, (PowerSeriesZ_t*) param1, (PowerSeriesZ_t*) param2);
}


/**
 * Given two power series f and g, it returns the difference f-g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeriesZ_t* subPowerSeries_PSZ(PowerSeriesZ_t* f,  PowerSeriesZ_t* g) {

    int i;
    int mindeg;
    if (g->deg == -1 && f->deg == -1) {
        return zeroPowerSeries_PSZ();
    //mindeg = 0;
    } else if (g->deg == -1) {
        //reserve_PSZ(f);
        //return f;
        mindeg = f->deg;//negatePowerSeries_PSZ
    } else if (f->deg == -1){
        //reserve_PSZ(g);return negatePowerSeries_PSZ(g);
        mindeg = g->deg;
    } else {
        mindeg = MIN(f->deg, g->deg);
    }
    PowerSeriesZ_t* sub = allocatePowerSeries_PSZ(mindeg+1);
    sub->genParam1 = f;
    sub->genParam2 = g;
    sub->paramType1 = POWER_SERIES;
    sub->paramType2 = POWER_SERIES;

    reserve_PSZ(f);
    reserve_PSZ(g);
    sub->genOrder = 2;
    sub->gen.binaryGen = &(homogPartVoid_sub_PSZ);
    sub->deg = mindeg;
    for (i = 0; i <= mindeg; i++){
        sub->polys[i] = homogPart_sub_PSZ(i, f, g);
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
PolyZ_ptr homogPart_negate_PSZ(int d, PowerSeriesZ_t* f)  {
    PolyZ_ptr poly = deepCopyPolynomial_AAZ(homogPart_PSZ(d, f));
    negatePolynomial_AAZ(poly);
    return poly;
}


/**
 * An internal function.
 * A void generator wrapper for the negation of a power series.
 *
 * @param d: the degree of the homogeneous part to generate.
 * @param param1: the power series to negate.
 */
PolyZ_ptr homogPartVoid_negate_PSZ(int d, void* param1) {
    return homogPart_negate_PSZ(d, (PowerSeriesZ_t*) param1);
}


/**
 * Negate the power series f, returning a lazily constructed power series.
 * @param f: the power series to negate.
 * @return the negation of the input power series.
 */
PowerSeriesZ_t* negatePowerSeries_PSZ(PowerSeriesZ_t* f) {

    PowerSeriesZ_t* ps = allocatePowerSeries_PSZ(f->deg+1);
    ps->deg = f->deg;
    ps->genOrder = 1;
    ps->genParam1 = f;
    reserve_PSZ(f);
    ps->paramType1 = POWER_SERIES;

    ps->gen.unaryGen = &(homogPartVoid_negate_PSZ);
    for (int i = 0; i <= f->deg; i++){

        if (f->polys[i] == NULL) {
            ps->polys[i] = NULL;
        } else {
            ps->polys[i] = homogPart_negate_PSZ(i, f);
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
PolyZ_ptr homogPart_prod_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "Entering homogPart_prod_PSZ: %d with %p, %p\n",d, f, g);
#endif
    PolyZ_ptr prod;
    PolyZ_ptr sum = NULL;

    // const char* syms[] = {"z", "x", "y"};
    PolyZ_ptr tmp1, tmp2;
    for (int i = 0; i <= d; i++) {
        tmp1 = homogPart_PSZ(d-i, f);
        if (isZero_AAZ(tmp1)) {
            continue;
        }
        tmp2 = homogPart_PSZ(i, g);
        if (isZero_AAZ(tmp2)) {
            continue;
        }
        // fprintf(stderr, "tmp1: ");
        // printPoly_AA(stderr, tmp1, syms, tmp1->nvar);
        // fprintf(stderr, "\ntmp2: ");
        // printPoly_AA(stderr, tmp2, syms, tmp2->nvar);
        // fprintf(stderr, "\n");
        prod = multiplyPolynomials_AAZ(tmp1, tmp2, tmp1->nvar);
        sum = addPolynomials_AAZ_inp(sum, prod, prod->nvar); //prod is guaranteed to not be 0 by above checks.
        freePolynomial_AAZ(prod);
    }
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "\n" );
    AltArr_t localSum;
    localSum.nvar = -12;
    if (sum != NULL) {
        localSum = *sum;
    }
    printPoly_AAZ(stderr, sum, syms, sum == NULL? 0 : sum->nvar);
    fprintf(stderr, "\n");
    fprintf(stderr, "Leaving homogPart_prod_PSZ: %d with %p, %p\n",d, f, g);
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
PolyZ_ptr homogPartVoid_prod_PSZ(int d, void* param1, void* param2) {
    return homogPart_prod_PSZ(d, (PowerSeriesZ_t*) param1, (PowerSeriesZ_t*) param2);
}

/**
 * Given two power series f and g, it returns the product, constructed lazily.
 * @param f: the multiplier power series
 * @param g: the multiplicand power series
 * @return a pointer to the resulting power series product.
 */
PowerSeriesZ_t* multiplyPowerSeries_PSZ( PowerSeriesZ_t* f,  PowerSeriesZ_t* g) {

    int i;
    int mindeg;
     /*Later I will remove extra if statements*/
    if (g->deg == -1 && f->deg == -1) {
        //mindeg = 0;
        return zeroPowerSeries_PSZ();
    } else if (g->deg == -1) {
       // mindeg = 0;//f->deg;
        return zeroPowerSeries_PSZ();
    } else if (f->deg == -1){
       // mindeg = 0;//g->deg;
        return zeroPowerSeries_PSZ();
    } else {
        mindeg = MIN(f->deg, g->deg);

    }

    PowerSeriesZ_t* prod = allocatePowerSeries_PSZ(mindeg+1);
    prod->genParam1 = f;
    prod->genParam2 = g;
    prod->paramType1 = POWER_SERIES;
    prod->paramType2 = POWER_SERIES;

    reserve_PSZ(f);
    reserve_PSZ(g);
    prod->genOrder = 2;
    prod->gen.binaryGen = &(homogPartVoid_prod_PSZ);
    prod->deg = mindeg;
    for (i = 0; i <= (mindeg); i++){
        prod->polys[i] = homogPart_prod_PSZ(i, f, g);
    }

    return prod;
}


/**
 * An internal function.
 * Computes homogeneous part of the quotient of two power series of degree d.
 * @param deg: the requested degree to generate
 * @param f: the dividend
 * @param h: the divisor, which must be a unit.
 * @param quo, the quotient being generated.
 * @return the homogeneous part of degree d of the quotient.
 */
PolyZ_ptr homogPart_quo_PSZ(int deg,  PowerSeriesZ_t* f,  PowerSeriesZ_t* h, PowerSeriesZ_t* quo) {

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
   fprintf(stderr, "requested degree quo: %d\n", deg);
#endif

    PolyZ_ptr gd;
    int i;
    PolyZ_ptr h0 = homogPart_PSZ(0, h);
    if (h0 == 0){
        fprintf(stderr, "ERROR: In homogPart_quo_PSZ, the divisor was zero!\n");
        exit(1);
    }

    if (deg <= 1) {

        PolyZ_ptr f0 = homogPart_PSZ(0,f);
        PolyZ_ptr q0 = divideByInteger_AAZ(f0, h0->elems->coef);

        if (deg == 0){
            return q0;
        }

        PolyZ_ptr h1 = homogPart_PSZ(1,h);
        PolyZ_ptr f1 = homogPart_PSZ(1,f);
        PolyZ_ptr multipoly = multiplyPolynomials_AAZ(q0, h1, q0 == NULL ? 0 : q0->nvar);
        freePolynomial_AAZ(q0);
        multipoly = subPolynomials_AAZ_inp( multipoly, f1,  f1 == NULL ? 0 : f1->nvar);
        negatePolynomial_AAZ(multipoly);
        divideByIntegerExact_AAZ_inp(multipoly, h0->elems->coef);
        // PolyZ_ptr homogpoly = homogeneousPart_AAZ(multipoly, 1);
        // freePolynomial_AAZ(multipoly);
        return multipoly;
    }

    PolyZ_ptr s = deepCopyPolynomial_AAZ(homogPart_PSZ(deg,f));

    for (i = 1; i <= deg; ++i) {
        PolyZ_ptr homog = homogPart_PSZ(i,h);
        //PolyZ_ptr homogquo = homogPart_quo_PSZ(deg-i, f, h);

        PolyZ_ptr homogquo = homogPart_PSZ(deg-i, quo);
        PolyZ_ptr p = multiplyPolynomials_AAZ(homog, homogquo, homog == NULL ? 0 : homog->nvar);
        //freePolynomial_AA(homogquo);
        s = subPolynomials_AAZ_inp(s, p, s == NULL ? 0 : s->nvar);
        freePolynomial_AAZ(p);
    }

    PolyZ_ptr qlast = NULL;
    PolyZ_ptr rlast = NULL;

    divideByIntegerExact_AAZ_inp(s, h0->elems->coef);
    // gd = homogeneousPart_AAZ(s, deg);
    // freePolynomial_AAZ(s);
    // return gd;
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
PolyZ_ptr homogPartVoid_quo_PSZ(int d, void* param1, void* param2, void* param3) {
  return homogPart_quo_PSZ(d, (PowerSeriesZ_t*) param1, (PowerSeriesZ_t*) param2, (PowerSeriesZ_t*) param3);
}


/**
 * Given two power series f and g, it returns the quotient f/g, constructed lazily.
 * @param f: the dividend power series
 * @param g: the divisor power series
 * @return a pointer to the resulting power series quotient.
 */
PowerSeriesZ_t* dividePowerSeries_PSZ(PowerSeriesZ_t* f,  PowerSeriesZ_t* h)  {

    PowerSeriesZ_t* quo = allocatePowerSeries_PSZ(1);
    quo->genParam1 = f;
    quo->genParam2 = h;
    quo->genParam3 = quo;

    quo->paramType1 = POWER_SERIES;
    quo->paramType2 = POWER_SERIES;
    quo->paramType3 = PLAIN_DATA; //Yes plain_data. Avoids chicken-and-egg with circular reference.

    reserve_PSZ(f);
    reserve_PSZ(h);
    quo->genOrder = 3;
    quo->gen.tertiaryGen = &(homogPartVoid_quo_PSZ);
    quo->deg = 0;
    //quo->polys[0] = homogPart_quo_PSZ(0, f, h);WAS
    quo->polys[0] = homogPart_quo_PSZ(0, f, h, quo);
    return quo;
}


