#ifndef _UPOPSZ_Z_H_
#define _UPOPSZ_Z_H_

#include "PowerSeriesZ.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The univariate polynomial over power series struct.
 * The polynomial is represented densely: the of power series coefficients
 * have their index equal to the degree of their corresponding monomial.
 *
 * NOTE: For ease of conversion and construction, it is assumed that
 * the polynomials representing the power series coefficients include the
 * polynomial variable of the upops.
 * That is, a upops which looks to an element of Q[[X_1,...,X_n]][Y] is actually represented
 * as Q[[Y,X_1,...,X_n]][Y], with each underlying polynomial having degree 0 in Y.
 *
 * @see convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ
 * @see
 *
 * deg : the degree of the univariate polynomial over power series (Upops)
 * alloc : the size of the array of power series associated with the Upops
 * data : the array of power series associated with the Upops
 * genParam1 and genParam2 are the ancestors of a Upops
 * construct a Upops
 */
typedef struct UPOPS_Z {
    int deg;
    int alloc;
    PowerSeriesZ_t** data;

    PowerSeriesZ_t** weierstrassFData;
    int fDataSize;

    int refCount;
} UpopsZ_t;

/**
 * Create a upops with a default allocation size.
 * @param alloc : the size of the array of power series associated with the Upops
 */
UpopsZ_t* allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(int alloc);

/**
 * Destroy the input upops. Actually, decrement its reference count and destroy conditionally.
 * @param upops : a given univariate polynomial over power series
 */
void destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(UpopsZ_t* upops);

/**
 * Increment the reference count of a upops.
 * @param upops : a given univariate polynomial over power series
 */
void reserve_UPOPSZ(UpopsZ_t* upops);

/**
 * Construct the Upops encoding zero.
 * @return 0 as a upops
 */
UpopsZ_t* zero_UPOPSZ();

/**
 * Construct one Upops encoding one.
 * @param nvar: the number of variables of the underlying polynomial ring.
 *              Recall that the upops polynomial variable is included as the first variable therein.
 * @return 1 as a upops
 */
UpopsZ_t* one_UPOPSZ(int nvar);

/**
 * Update the power series coefficients of the upops f
 * so that they are all at least of degree d.
 *
 * @param f: the upops to update
 * @param d: the requested degree.
 */
void  updateToDeg_UPOPSZ(int d, UpopsZ_t* f);

/**
 * upops : a Upops
 * sym : a list of variables
 * print a Upops
 */
void print_UPOPSZ(FILE* fp, UpopsZ_t* upops, const char** sym);


/***********************
 * Conversion Helpers
 ***********************/

/**
 * Convert a polynomial to a univariate polynomial with power series coefficients.
 * The variable of index 0 in p is used as the univariate polynomial variable.
 * Then, the polynomials in the power series coefficients remain
 * to be defined in their original polynomial ring, but have degree 0 in their main variable.
 *
 * @param p : a polynomial to convert.
 * @return p converted to a upops.
 */
UpopsZ_t* convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(const PolyZ_ptr p);

/**
 * Convert an array of power series to a univariate polynomial over power series (upops).
 * Where the index of the power series implies its associated monomial's degree.
 * @note the power series should be defined to exist in Q[[Y,X_1,...,X_n]]
 * but have degree 0 in Y, the eventual polynomial variable of the upops.
 *
 * The power series in the input array are copied (reserved) and so ownership
 * is *not* transferd to the resulting upops.
 *
 * @see convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ
 *
 * @param ps : an array of power series
 * @param size : the size of the array
 * @return a upops whose coefficients are that of the array ps.
 */
UpopsZ_t* convertArrayOfPSToUPOPS_UPOPSZ(PowerSeriesZ_t** ps, int size);

/**
 * Get the polynomial part of a upops, truncating its coefficients
 * at degree d and converting it to a multivarate polynomial.
 *
 * f : a upops
 * d : the degree
 * nvar : the number of variables
 * compute polynomial part of the given power series
 */
PolyZ_ptr polynomialPart_UPOPSZ(int d, UpopsZ_t* f);



/***********************
 * UPOPS Arithmetic
 ***********************/

/**
 * Add two upops, returning the sum constructed lazily.
 * f : a Upops
 * g : a Upops
 * @return the sum f + g
 */
UpopsZ_t*  addUnivariatePolyOverPowerSeries_UPOPSZ(UpopsZ_t* f, UpopsZ_t* g);


/**
 * Subtract two upops, returning the difference constructed lazily.
 * f : a Upops
 * g : a Upops
 * @return the difference f - g
 */
UpopsZ_t*  subUnivariatePolyOverPowerSeries_UPOPSZ(UpopsZ_t* f, UpopsZ_t* g);


/**
 * Multiply two upops, returning the product constructed lazily.
 * f : a Upops
 * g : a Upops
 * @return the product f * g
 */
UpopsZ_t*  multiplyUnivariatePolyOverPowerSeries_UPOPSZ(UpopsZ_t* f, UpopsZ_t* g);




#ifdef __cplusplus
}
#endif

#endif

