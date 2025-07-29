// Include guard
#ifndef _LODEoUPS_H_
#define _LODEoUPS_H_

// Includes from BPAS references
#include <gmp.h>
#include "Parser/bpas_parser.h"
#include "Utils/Unix_Timer.h"
#include "PowerSeries/PowerSeries.h"

// My includes
#include <sys/time.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <assert.h>
#include <time.h>
#include <string.h>

// So this can work in C++
#ifdef __cplusplus
extern "C" {
#endif


/**
 * Macro to return the nth coefficient of a univariate power series as a rational number
 * Only works for use with LODEoUPS structure called lode
 * @param ps the power series, of type PowerSeries_t *
 * @param n integer n to return the nth coefficient 
 * */ 
#define GETCOEF_MPQ(ps, n) (*get_univariate_power_series_coefficient_mpq(ps, n, &(lode->zero_mpq)))

// Macros to get the power series coefficients a[j]_n as an mpq_t rational
#define GC_A(j,n) (GETCOEF_MPQ(lode->coefs[j], n))
#define GC_Y(n) (GETCOEF_MPQ(lode->sol, n)) // Get coefficient n of solution power series, that is c[n]


// Used for recursive construction of ODE solving for the base case for serial recursion
// This should be set to optimize the cache performance
#define SERIAL_SIZE_BASE 8192




/**
 * @struct 
 * @brief 
 * Structure for a Linear ODE with power series coefficients
 * 
 * Fields:
 * @param sol The power series solution of the ODE
 * @param ord: The order of the ODE
 * @param inits: an array of initial conditions as rational numbers mpq_t
 * @param force: the right hand side (forcing function / nonhomogeneous term) of the linear ODE


*/
typedef struct LinearODEOverUnivariatePowerSeries {

    PowerSeries_t * sol;

    int ord;

    PowerSeries_t ** coefs;

    PowerSeries_t * forcing;

    mpq_t * inits; // mpq_t *

    int update_deg;

    mpq_t ** rf;

    mpq_t zero_mpq;


} LODEoUPS_t;


/**
 * Constructor for LODEoUPS structure
 * @param ord the order of the ODE
 * @param coefs an array of univariate power series which represent the coefficients of the ode. These will be deep copied
 * @param forcing the forcing function of the ODE (the LHS of the equation)
 * @param inits the initial conditions for the ODE, these should probably be of polynomial type
*/
LODEoUPS_t * construct_LODEoUPS(int ord, PowerSeries_t ** coefs, PowerSeries_t * forcing, ratNum_t * inits);

/**
 * @param ord : the order of the DE
 * Allocates necessary memory for a LODEoUPS
 */
LODEoUPS_t* allocate_LODEoUPS(int ord); 

/**
 * Function to deallocate and clear memory for LODEoUPS structure
 * @param lode the LODEoUPS structure
 */
void destroy_LODEoUPS(LODEoUPS_t * lode);

/**
 * Function to return the coefficient of a power series
 * Had to make this instead of a macro for the case where power series do not have a generator that returns 0
 * @param ps the power series to return the coefficient from
 * @param n the degree of the power series coefficient to return
 * @param zero_mpq an mpq_t which has the value 0 (0/1). The user has to manually manage (clear) the memory used by this mpq
 */
mpq_t * get_univariate_power_series_coefficient_mpq(PowerSeries_t * ps, int n, mpq_t * zero_mpq);

/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series solution for a Linear ODE with Univariate Power Series coefficients
 * @param d: the requested degree to generate
 * @param param1: the LODEoUPS structure
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPartVoid_LODEoUPS_sol_PS(int d, void* param1);

/**
 * An internal function.
 * Computes the homogeneous part of the
 * power series solution for a Linear ODE with Univariate Power Series coefficients
 * @param d: the requested degree to generate
 * @param param1: the LODEoUPS structure
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPart_LODEoUPS_sol_PS(int d, LODEoUPS_t * lode);

/**
 * Solves a linear ODE up to the given degree using the generator for the solution power series
 * @param deg the precision for the power series solution to be calculated up to
 * @param lode the linear ODE structure to solve
 */
void solve_LODE_upto_degree_using_sol_generator(int deg, LODEoUPS_t * lode);


/**
 * Helper function to make a univariate monomial of degree d with a specific rational coefficient
 * @param deg the degree of the univariate monomial
 * @param coef the rational coefficient of the univariate monomial
 */
Poly_ptr make_univariate_monomial(int deg, mpq_t coef);


/**
 * Internal function
 * This function will update the LODE coefficients and rising factorial table to 
 * The depth required to calculate a desired degree of the power series solution
 * This helps so these operations are only done once, and not many times in the generator
 * @param deg the desired degree of the power series solution
 * @param lode the LODE structure
 */
void update_LODE_auxiliaries(int deg, LODEoUPS_t * lode);


/**
 * Internal function
 * Turns the initial conditions given to the ODE which are y(0), y'(0), y''(0), etc
 * Into coefficients in the power series solution, that is, divide the initial condition by a factorial
 * @param lode the LODE structure for which we are correcting the intial conditions
 */
void turn_LODE_initial_conditions_to_coefficients(LODEoUPS_t * lode);











/**
 * Solves a linear ODE up to the given degree, using a method that treats the degree as a fixed precision
 * The method assumes that no prior coefficients for the power series solution have been computed
 * @param deg the precision for the power series solution to be calculated up to
 * @param lode the linear ODE structure to solve
 */
void solve_LODE_upto_degree_truncated_construction(int deg, LODEoUPS_t * lode);


/**
 * An internal function
 * Recursively calculates the "T" table which at position T[i][i] holds the coefficient  c[d+i] of the power series
 * Solution to an ODE
 * This is the rectangular recursive call, for when the section of the T table being calculated is a rectangle
 * @param left the leftmost index in the T table being calculated in this recursive call
 * @param right the rightmost index in the T table being calculated in this recursive call
 * @param low the lowest vertical index in the T table being calculated in this recursive call
 * @param high the highest vertical index in the T table being calculated in this recursive call
 * @param d the order of the differential equation
 * @param w the w vector, which contains the constant terms 1/rf[d][i] * 1/a[d][0]
 * @param g the g table
 * @param c an mpq_t vector which contains the coefficients of the power series solution to the ODE
 * @param T the T table being calculated

 */
void calculate_T_table_rectangular(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T);


/**
 * An internal function
 * Recursively calculates the "T" table which at position T[i][i] holds the coefficient  c[d+i] of the power series
 * Solution to an ODE
 * This is the triangular recursive call, for when the section of the T table being calculated is triangular
 * @param left the leftmost index in the T table being calculated in this recursive call
 * @param right the rightmost index in the T table being calculated in this recursive call
 * @param low the lowest vertical index in the T table being calculated in this recursive call
 * @param high the highest vertical index in the T table being calculated in this recursive call
 * @param d the order of the differential equation
 * @param w the w vector, which contains the constant terms 1/rf[d][i] * 1/a[d][0]
 * @param g the g table
 * @param c an mpq_t vector which contains the coefficients of the power series solution to the ODE
 * @param T the T table being calculated

 */
void calculate_T_table_triangular(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T);












/**
* Create the table (2-D array) of rising factorials using simple iterative construction
* @param d the order of the DE for which the rising factorial table is being created
* @param k the precision (minus order) of the coefficient of the power series solution being calculated
*/
mpq_t ** rising_factorial_table_mpq_t (int d, int k);


/**
* Deallocate the memory allocated for the rising factorial table
* @param d the order of the DE for which the rising factorial table is being created
* @param k the precision (minus order) of the coefficient of the power series solution being calculated
*/
void destroy_rising_factorial_table_mpq_t(mpq_t ** rf, int d, int k);


/**
* Create the table (2-D array) of rising factorials using recursive construction
* @param d the order of the DE for which the rising factorial table is being created
* @param k the precision (minus order) of the coefficient of the power series solution being calculated
*/
mpq_t ** rising_factorial_table_recursive_mpq_t (int dd, int kk);

/*
Recursive function to call to create rising factorial table recursively dividing only along witdh (treat height as constant)
*/
void rising_factorial_mpq_t (int left_i, int right_i, mpq_t ** rf, int height_i);







/**
 Function to generate the mclaurin series for a univariate exponential function
 Of the form A*e^(B*x) --> front_factor * e ^ (inner_factor * x)
 @param front_factor the front factor of the exponential function so front_factor * e^x
 @param inner_factor the inner factor of the exponential function so e^(inner_factor * x)
*/
PowerSeries_t * exponential_univariate(mpq_t* front_factor, mpq_t* inner_factor);

/**
 * An internal function.
 * A void generator wrapper for the generator of a power series of the exponential function A * e ^ Bx
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPartVoid_expUnivariate_PS(int d, void* param1, void* param2);

/**
 * An internal function.
 * Computes the homogeneous part of a power series of the exponential funciton A * e ^ Bx
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPart_expUnivariate_PS(int d, mpq_t* front_factor , mpq_t* inner_factor); 





/**
 Function to generate the mclaurin series for a univariate exponential function with x to the exponent of some natural number
 Of the form A*e^(B*x^n) --> front_factor * e ^ (inner_factor * x^n)
 @param front_factor the front factor of the exponential function so front_factor * e^(x^n)
 @param inner_factor the inner factor of the exponential function so e^(inner_factor * x^n)
 @param nat the exponent of x within the exponential function, so e^(inner_factor * x^n)
*/
PowerSeries_t * exponential_univariate_natural_power(mpq_t* front_factor, mpq_t* inner_factor, int nat);

/**
 * An internal function.
 * A void generator wrapper for the generator of a power series of the exponential function
 * With x to the exponent of some natural number
 *  A * e ^ (B*x^n)
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 * @param nat the exponent of x within the exponential function
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPartVoid_expUnivariate_natural_power_PS(int d, void* param1, void* param2, void* param3);

/**
 * An internal function.
 * Computes the homogeneous part of a power series of the exponential funciton
 * With x to the exponent of some natural number
 *  A * e ^ (B*x^n)
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 *  @param nat the exponent of x within the exponential function
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPart_expUnivariate_natural_power_PS(int d, mpq_t* front_factor , mpq_t* inner_factor, int nat);










// So this can work in C++
#ifdef __cplusplus
}
#endif


// Include guard
#endif