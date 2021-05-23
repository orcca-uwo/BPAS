#ifndef _POWERSERIESZ_H_
#define _POWERSERIESZ_H_

#include "IntegerPolynomial/SMZP_Support.h"
#include "PowerSeries.h"

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Typedef for polynomial.
 * All polynomials are pointers to structs.
 */
typedef AltArrZ_t* PolyZ_ptr;

/**
 * Typedefs for power series generators.
 */
typedef PolyZ_ptr (*homog_partZ_gen) (int);
typedef PolyZ_ptr (*homog_partZ_gen_unary) (int, void*);
typedef PolyZ_ptr (*homog_partZ_gen_binary) (int, void*, void*);
typedef PolyZ_ptr (*homog_partZ_gen_tertiary) (int, void*, void*, void*);

/**
 * A union type for all possible homogeneous part generators.
 */
typedef union HomogPartZGenerator {
    homog_partZ_gen nullaryGen;
    homog_partZ_gen_unary unaryGen;
    homog_partZ_gen_binary binaryGen;
    homog_partZ_gen_tertiary tertiaryGen;
} HomogPartZGenerator_u;



/**
 * The lazy multivariate power series struct for power series over Q.
 * It has an array of homogeneous polynomials whose index is their
 * degree.
 * Each polynomial in polys should have the same number of variables
 * and exist in the same polynomial ring.
 *
 * deg: the highest known degree of homogeneous polynomials in polys
 * alloc: the allocation size of polys
 * polys: the array of homogeneous polys
 * gen: a function pointer to the polynomial generator
 * genOrder: the order of the generator (how many void* params it has)
 * refCount* the reference count of this power series.
 */
typedef struct PowerSeriesZ {
    int deg;
    int alloc;
    PolyZ_ptr* polys;

    HomogPartZGenerator_u gen;
    int genOrder;
    void* genParam1;
    void* genParam2;
    void* genParam3;
    GenParamType_e paramType1;
    GenParamType_e paramType2;
    GenParamType_e paramType3;

    int refCount;

} PowerSeriesZ_t;

/**
 * Create a power series struct with some default allocation size.
 * @param alloc: the size of the array of polys to allocate.
 * @return the newly created power series
 */
PowerSeriesZ_t* allocatePowerSeries_PSZ(int alloc);


/**
 * Destroy the power series.
 * Actually, decerements the reference count and destroys conditionally.
 * @param ps : a power series
 */
void destroyPowerSeries_PSZ(PowerSeriesZ_t* ps);


/**
 * Increment the refrence count of the input power series
 * @param ps : a power series
 */
void reserve_PSZ(PowerSeriesZ_t* ps);


/**
 * Given a power series f, check whether or not f is a unit.
 * @return 1 iff f is a unit.
 */
static inline int isUnit_PSZ(PowerSeriesZ_t* f) {
    if (f == NULL || f-> deg < 0 || isZero_AAZ(f->polys[0])) {
        return 0;
    }

    return isOne_AAZ(f->polys[0]) || isNegativeOne_AAZ(f->polys[0]);
}

/**
 * Given a power series f, it check whether or not f is truly 0.
 * That is, is it 0 and is it idendtically 0, not simply
 * 0 because it was lazily constructed.
 *
 * @return 1 iff f is 0.
 */
static inline int isZero_PSZ(PowerSeriesZ_t* f) {
    if (f == NULL || f->deg == -1) {
        return 1;
    }
    return 0;
}


/**
 * Create the zero power series.
 * @note: 0 is 0 in all underlying polynomial rings, thus no nvar parameter needed.
 * @return a pointer to the zero power series.
 */
PowerSeriesZ_t* zeroPowerSeries_PSZ();


/**
 * Create the one power series.
 * @param nvar: the number of variables for the underlying power series ring.
 * @return a pointer to the one power series.
 */
PowerSeriesZ_t* onePowerSeries_PSZ(int nvar);


/**
 * Create a constant power series from the constnat coef.
 * @param coef : a constant
 * @param nvar : the number of variables for the underyling polynomial ring.
 * @return a pointer to the constant power series.
 */
PowerSeriesZ_t* constPowerSeries_PSZ(const mpz_t coef, int nvar);


/**
 * Create the geometric series of nvar number of variables as a power series.
 * @param nvar : the number of variables
 * @return the geometric series a power series.
 */
PowerSeriesZ_t* geometricSeries_PSZ(long long nvar);

/**
 * Print a power series to the file pointer fp
 * using the symbols sym as symbols of the homogeneous polynomial.
 * @param fp: the file pointer to print to; may be stdout or stderr or something else.
 * @param ps : a power series
 * @param sym : a list of variables
 */
void print_PSZ(FILE* fp, PowerSeriesZ_t* ps, const char** sym);


/**********************
 * Main functional interface
 **********************/

/**
 * Given a power series f and an integer d, return the homogeneous
 * part of f of degree d.
 *
 * Note that, for performance reasons, this returns a pointer to the
 * same polynomial as is stored in the power series. It should not be free'd.
 * Modifications should be made only to a copy of the returned polynmoial.
 * @see deepCopyPolynomial_AA
 *
 * The is the main functional interface for power series.
 * Calling this function will generate terms in the power series as needed.
 * @see polynomialPart_PSZ.
 *
 * @return the homogeneous part of f of degree d.
 */
PolyZ_ptr homogPart_PSZ(int d, PowerSeriesZ_t* f);

/**
 * Computes the polynomial part of a power series up to degree d
 * @param d: the degree of the resulting polynomial part
 * @param f: a power series
 *
 * @return the polynomial part of f of degree d.
 */
PolyZ_ptr polynomialPart_PSZ(int d,   PowerSeriesZ_t* f);

/**
 * Update the given power series up and including to the requested degree d.
 * @param f : the power series
 * @param d : the degree
 */
static inline void updateToDeg_PSZ(int d, PowerSeriesZ_t* f) {
    homogPart_PSZ(d, f);
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
PolyZ_ptr* convertPolyToArrayOfHomogeneous_PSZ(const PolyZ_ptr p);

/**
 * Eliminates terms of degree atleast d+1 from polynomial p.
 * @param p: the polynomial to truncate
 * @param d: the (maximum) degree of the resulting truncated polynomial.
 * @return the truncated polynomial
 */
PolyZ_ptr truncatePolynomial_PSZ(const PolyZ_ptr p, int d);

/**
 * Converts a polynomial to a power series.
 *
 * @param p: the polynomial to convert
 * @return the power series encoding the input polynomial.
 */
PowerSeriesZ_t* convertPolyToPowerSeries_PSZ(const PolyZ_ptr p);


/**
 * Return an array of homogeneous polynomials, each polynomial the
 * homogeneous part of f from degree 0 up to and including d.
 * f : the power series
 * d : the integer
 */
static inline PolyZ_ptr* copyUpTo_PSZ(PowerSeriesZ_t* f, int d) {
    PolyZ_ptr* arrayOfHomogPart = (PolyZ_ptr*) malloc(sizeof(PolyZ_ptr)*(d+1));
    for (int i = 0; i <= d; ++i) {
        arrayOfHomogPart[i] =  deepCopyPolynomial_AAZ(homogPart_PSZ(i, f));
    }
    return arrayOfHomogPart;
}


/**********************
 * Addition and Subtraction
 **********************/


/**
 * Given two power series f and g, it returns the sum f+g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeriesZ_t* addPowerSeries_PSZ(PowerSeriesZ_t* f,  PowerSeriesZ_t* g);

/**
 * Given two power series f and g, it returns the difference f-g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeriesZ_t* subPowerSeries_PSZ(PowerSeriesZ_t* f,  PowerSeriesZ_t* g);

/**
 * Negate the power series f, returning a lazily constructed power series.
 * @param f: the power series to negate.
 * @return the negation of the input power series.
 */
PowerSeriesZ_t* negatePowerSeries_PSZ(PowerSeriesZ_t* f);


/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series sum.
 * @param d: the requested degree to generate
 * @param param1: the left-hand side of the addition.
 * @param param2: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
PolyZ_ptr homogPartVoid_sum_PSZ(int deg, void* param1, void* param2);

/**
 * An internal function.
 * Computes homogeneous part of the sum of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
PolyZ_ptr homogPart_sum_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g);

/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series difference.
 * @param d: the requested degree to generate
 * @param param1: the left-hand side of the difference.
 * @param param2: the right-hand side of the difference
 * @return the homogeneous part of degree d of the sum.
 */
PolyZ_ptr homogPartVoid_sub_PSZ(int deg, void* param1, void* param2);

/**
 * An internal function.
 * Computes homogeneous part of the sum of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
PolyZ_ptr homogPart_sub_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g);

/**
 * An internal function.
 * A void generator wrapper for the negation of a power series.
 *
 * @param d: the degree of the homogeneous part to generate.
 * @param param1: the power series to negate.
 */
PolyZ_ptr homogPartVoid_negate_ps(int d, void* param1);

/**
 * An internal function.
 * Generate the homogeneous part of degree d of the
 * negation of power series f.
 * d : the degree of the homogeneous part to generate
 * f : the power series to negate
 * @return the homogeneous part of degree d
 */
PolyZ_ptr homogPart_negate_PSZ(int d, PowerSeriesZ_t* f);



/**********************
 * Multiplication and Division
 **********************/

/**
 * Given two power series f and g, it returns the product, constructed lazily.
 * @param f: the multiplier power series
 * @param g: the multiplicand power series
 * @return a pointer to the resulting power series product.
 */
PowerSeriesZ_t* multiplyPowerSeries_PSZ( PowerSeriesZ_t* f,  PowerSeriesZ_t* g);

/**
 * Given two power series f and g, it returns the quotient f/g, constructed lazily.
 * @param f: the dividend power series
 * @param g: the divisor power series
 * @return a pointer to the resulting power series quotient.
 */
PowerSeriesZ_t* dividePowerSeries_PSZ(PowerSeriesZ_t* f,  PowerSeriesZ_t* h);

/**
 * Given power series f, return 1/f, constructed lazily.
 * @param f: the divisor power series
 * @return a pointer to the resulting power series quotient.
 */
static inline PowerSeriesZ_t* inversePowerSeries_PSZ(PowerSeriesZ_t* f) {
    if (f == NULL) {
        return NULL;
    }
    int nvar = f->polys[f->deg]->nvar;
    PowerSeriesZ_t* one = onePowerSeries_PSZ(nvar);
    PowerSeriesZ_t* ret = dividePowerSeries_PSZ(one, f);
    destroyPowerSeries_PSZ(one);
    return ret;
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
PolyZ_ptr homogPartVoid_prod_PSZ(int d, void* param1, void* param2);


/**
 * An internal function.
 * Computes homogeneous part of the product of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand operand of the multiplication
 * @param g: the right-hand operand of the multiplication
 * @return the homogeneous part of degree d of the product.
 */
PolyZ_ptr homogPart_prod_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g);

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
PolyZ_ptr homogPartVoid_quo_PSZ(int d, void* param1, void* param2, void* param3);

/**
 * An internal function.
 * Computes homogeneous part of the quotient of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the dividend
 * @param g: the divisor, which must be a unit.
 * @param quo, the quotient being generated.
 * @return the homogeneous part of degree d of the quotient.
 */
PolyZ_ptr homogPart_quo_PSZ(int d,  PowerSeriesZ_t* f,  PowerSeriesZ_t* g, PowerSeriesZ_t* quo);




#ifdef __cplusplus
}
#endif


#endif
