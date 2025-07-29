// Include guard
#ifndef _LODEoUPS_PARALLEL_HPP_
#define _LODEoUPS_PARALLEL_HPP_

// BPAS Includes
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#include "LinearODEOverUnivariatePowerSeries.h"

// C++ includes
#include <iostream>
#include <chrono>
#include <gmpxx.h>
using namespace std;


// Base Case size for parallel work in calculating the T table
#define PARALLEL_WORK_BASE 4096 


/**
 * Solves a linear ODE up to the given degree, using a method that treats the degree as a fixed precision
 * The method assumes that no prior coefficients for the power series solution have been computed
 * @param deg the precision for the power series solution to be calculated up to
 * @param lode the linear ODE structure to solve
 */
void solve_LODE_upto_degree_truncated_parallel(int deg, LODEoUPS_t * lode);


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
void calculate_T_table_triangular_parallel(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T);


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
void calculate_T_table_rectangular_parallel(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T);







// Include guard
#endif