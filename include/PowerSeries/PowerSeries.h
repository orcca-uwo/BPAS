#ifndef _POWERSERIES_H_
#define _POWERSERIES_H_

#include "RationalNumberPolynomial/SMQP_Support-AA.h"


#ifdef __cplusplus
extern "C" {
#endif

/**
 * Typedef for polynomial.
 * All polynomials are pointers to structs.
 */
typedef AltArr_t* Poly_ptr;

/**
 * Typedefs for power series generators.
 */
typedef Poly_ptr (*homog_part_gen) (int);
typedef Poly_ptr (*homog_part_gen_unary) (int, void*);
typedef Poly_ptr (*homog_part_gen_binary) (int, void*, void*);
typedef Poly_ptr (*homog_part_gen_tertiary) (int, void*, void*, void*);

/**
 * Computes the min value of x and y.
 */
#define MIN(x, y) ((x) < (y) ? (x) : (y))

/**
 * A union type for all possible homogeneous part generators.
 */
typedef union HomogPartGenerator {
    homog_part_gen nullaryGen;
    homog_part_gen_unary unaryGen;
    homog_part_gen_binary binaryGen;
    homog_part_gen_tertiary tertiaryGen;
} HomogPartGenerator_u;

/**
 * An enumeration type for the possible parameters to a generator.
 */
typedef enum GenParamType_e {
    PLAIN_DATA = 0,
    POWER_SERIES = 1,
    UPOPS = 2,
    WEAK_UPOPS = 3, // A weak reference to a upops, this this param was not reserved
    MPQ_LIST = 4
} GenParamType_e;

/**
 * A helper struct for a ref-counted array of mpq_t values.
 */
typedef struct mpq_list {
    int size;
    int refCount;
    mpq_t* data;
} mpq_list_t;

/**
 * Decrement the refcount of a mpq_list_t and conditionally
 * destroy it.
 * @param mpql: the mpq_list_t to destroy
 */
static inline void destroyMPQList_PS(mpq_list_t* mpql) {
    if (mpql != NULL) {
        --(mpql->refCount);
        if (mpql->refCount <= 0) {
            for (int i = 0; i < mpql->size; ++i) {
                mpq_clear(mpql->data[i]);
            }
            free(mpql->data);
            free(mpql);
        }
    }
}

/**
 * Increment the ref count of a mpq_list_t.
 * @param mpql: the mpq_list_t to reserve.
 */
static inline void reserveMPQList_PS(mpq_list_t* mpql) {
    if (mpql != NULL) {
        ++(mpql->refCount);
    }
}

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
typedef struct PowerSeries {
    int deg;
    int alloc;
    Poly_ptr* polys;

    HomogPartGenerator_u gen;
    int genOrder;
    void* genParam1;
    void* genParam2;
    void* genParam3;
    GenParamType_e paramType1;
    GenParamType_e paramType2;
    GenParamType_e paramType3;

    int refCount;

} PowerSeries_t;

/**
 * Create a power series struct with some default allocation size.
 * @param alloc: the size of the array of polys to allocate.
 * @return the newly created power series
 */
PowerSeries_t* allocatePowerSeries_PS(int alloc);


/**
 * Destroy the power series.
 * Actually, decerements the reference count and destroys conditionally.
 * @param ps : a power series
 */
void destroyPowerSeries_PS(PowerSeries_t* ps);


/**
 * Increment the refrence count of the input power series
 * @param ps : a power series
 */
void reserve_PS(PowerSeries_t* ps);


/**
 * Given a power series f, check whether or not f is a unit.
 * @return 1 iff f is a unit.
 */
static inline int isUnit_PS(PowerSeries_t* f) {
    if (f == NULL || f-> deg < 0 || isZero_AA(f->polys[0])) {
        return 0;
    }

    return 1;
}

/**
 * Given a power series f, it check whether or not f is truly 0.
 * That is, is it 0 and is it idendtically 0, not simply
 * 0 because it was lazily constructed.
 *
 * @return 1 iff f is 0.
 */
static inline int isZero_PS(PowerSeries_t* f) {
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
PowerSeries_t* zeroPowerSeries_PS();


/**
 * Create the one power series.
 * @param nvar: the number of variables for the underlying power series ring.
 * @return a pointer to the one power series.
 */
PowerSeries_t* onePowerSeries_PS(int nvar);


/**
 * Create a constant power series from the constnat coef.
 * @param coef : a constant
 * @param nvar : the number of variables for the underyling polynomial ring.
 * @return a pointer to the constant power series.
 */
PowerSeries_t* constPowerSeries_PS(const mpq_t coef, int nvar);


/**
 * Create the geometric series of nvar number of variables as a power series.
 * @param nvar : the number of variables
 * @return the geometric series a power series.
 */
PowerSeries_t* geometricSeries_PS(long long nvar);

/**
 * Print a power series to the file pointer fp
 * using the symbols sym as symbols of the homogeneous polynomial.
 * @param fp: the file pointer to print to; may be stdout or stderr or something else.
 * @param ps : a power series
 * @param sym : a list of variables
 */
void print_PS(FILE* fp, PowerSeries_t* ps, const char** sym);


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
Poly_ptr homogPart_PS(int d, PowerSeries_t* f);

/**
 * Computes the polynomial part of a power series up to degree d
 * @param d: the degree of the resulting polynomial part
 * @param f: a power series
 *
 * @return the polynomial part of f of degree d.
 */
Poly_ptr polynomialPart_PS(int d,   PowerSeries_t* f);

/**
 * Update the given power series up and including to the requested degree d.
 * @param f : the power series
 * @param d : the degree
 */
static inline void updateToDeg_PS(int d, PowerSeries_t* f) {
    homogPart_PS(d, f);
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
Poly_ptr* convertPolyToArrayOfHomogeneous_PS(const Poly_ptr p);

/**
 * Eliminates terms of degree atleast d+1 from polynomial p.
 * @param p: the polynomial to truncate
 * @param d: the (maximum) degree of the resulting truncated polynomial.
 * @return the truncated polynomial
 */
Poly_ptr truncatePolynomial_PS(const Poly_ptr p, int d);

/**
 * Converts a polynomial to a power series.
 *
 * @param p: the polynomial to convert
 * @return the power series encoding the input polynomial.
 */
PowerSeries_t* convertPolyToPowerSeries_PS(const Poly_ptr p);


/**
 * Return an array of homogeneous polynomials, each polynomial the
 * homogeneous part of f from degree 0 up to and including d.
 * f : the power series
 * d : the integer
 */
static inline Poly_ptr* copyUpTo_PS(PowerSeries_t* f, int d) {
    Poly_ptr* arrayOfHomogPart = (Poly_ptr*) malloc(sizeof(Poly_ptr)*(d+1));
    for (int i = 0; i <= d; ++i) {
        arrayOfHomogPart[i] =  deepCopyPolynomial_AA(homogPart_PS(i, f));
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
PowerSeries_t* addPowerSeries_PS(PowerSeries_t* f,  PowerSeries_t* g);

/**
 * Given two power series f and g, it returns the difference f-g, constructed lazily.
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return a pointer to the resulting power series quotient.
 */
PowerSeries_t* subPowerSeries_PS(PowerSeries_t* f,  PowerSeries_t* g);

/**
 * Negate the power series f, returning a lazily constructed power series.
 * @param f: the power series to negate.
 * @return the negation of the input power series.
 */
PowerSeries_t* negatePowerSeries_PS(PowerSeries_t* f);


/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series sum.
 * @param d: the requested degree to generate
 * @param param1: the left-hand side of the addition.
 * @param param2: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPartVoid_sum_PS(int deg, void* param1, void* param2);

/**
 * An internal function.
 * Computes homogeneous part of the sum of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPart_sum_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g);

/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series difference.
 * @param d: the requested degree to generate
 * @param param1: the left-hand side of the difference.
 * @param param2: the right-hand side of the difference
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPartVoid_sub_PS(int deg, void* param1, void* param2);

/**
 * An internal function.
 * Computes homogeneous part of the sum of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand side of the addition
 * @param g: the right-hand side of the addition
 * @return the homogeneous part of degree d of the sum.
 */
Poly_ptr homogPart_sub_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g);

/**
 * An internal function.
 * A void generator wrapper for the negation of a power series.
 *
 * @param d: the degree of the homogeneous part to generate.
 * @param param1: the power series to negate.
 */
Poly_ptr homogPartVoid_negate_ps(int d, void* param1);

/**
 * An internal function.
 * Generate the homogeneous part of degree d of the
 * negation of power series f.
 * d : the degree of the homogeneous part to generate
 * f : the power series to negate
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPart_negate_PS(int d, PowerSeries_t* f);



/**********************
 * Multiplication and Division
 **********************/

/**
 * Given two power series f and g, it returns the product, constructed lazily.
 * @param f: the multiplier power series
 * @param g: the multiplicand power series
 * @return a pointer to the resulting power series product.
 */
PowerSeries_t* multiplyPowerSeries_PS( PowerSeries_t* f,  PowerSeries_t* g);

/**
 * Given two power series f and g, it returns the quotient f/g, constructed lazily.
 * @param f: the dividend power series
 * @param g: the divisor power series
 * @return a pointer to the resulting power series quotient.
 */
PowerSeries_t* dividePowerSeries_PS(PowerSeries_t* f,  PowerSeries_t* h);

/**
 * Given power series f, return 1/f, constructed lazily.
 * @param f: the divisor power series
 * @return a pointer to the resulting power series quotient.
 */
static inline PowerSeries_t* inversePowerSeries_PS(PowerSeries_t* f) {
    if (f == NULL) {
        return NULL;
    }
    int nvar = f->polys[f->deg]->nvar;
    PowerSeries_t* one = onePowerSeries_PS(nvar);
    PowerSeries_t* ret = dividePowerSeries_PS(one, f);
    destroyPowerSeries_PS(one);
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
Poly_ptr homogPartVoid_prod_PS(int d, void* param1, void* param2);


/**
 * An internal function.
 * Computes homogeneous part of the product of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the left-hand operand of the multiplication
 * @param g: the right-hand operand of the multiplication
 * @return the homogeneous part of degree d of the product.
 */
Poly_ptr homogPart_prod_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g);

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
Poly_ptr homogPartVoid_quo_PS(int d, void* param1, void* param2, void* param3);

/**
 * An internal function.
 * Computes homogeneous part of the quotient of two power series of degree d.
 * @param d: the requested degree to generate
 * @param f: the dividend
 * @param g: the divisor, which must be a unit.
 * @param quo, the quotient being generated.
 * @return the homogeneous part of degree d of the quotient.
 */
Poly_ptr homogPart_quo_PS(int d,  PowerSeries_t* f,  PowerSeries_t* g, PowerSeries_t* quo);




#ifdef __cplusplus
}
#endif


#endif
