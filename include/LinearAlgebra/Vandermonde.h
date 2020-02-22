


#ifndef _VANDERMONDE_Z_
#define _VANDERMONDE_Z_

#include "../ModularPolynomial/DUSP_Support.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Solve a system of equations described by a (generalized) Vandermonde matrix.
 * The Vandermonde matrix is defined by the array of distinct values t such that
 * the ith row of the matrix is powers of t[i]. Such powers increase by 1 
 * per column starting from either 0 or 1, depending on the value of startExp.
 *
 * A non-zero value of startExp makes the first column t[i]^1, otherwise the
 * first column is t[i]^0. 
 *
 * If transposed is a non-zero value then the matrix is constructed as above 
 * but then transposed before solving.
 * 
 * @note: Currently only transposed = 1 are valid options.
 *
 * @param t the distinct values defining the matrix.
 * @param b the right hand side vector.
 * @param size the size of the arrays t and b.
 * @param startExp a boolean describing to start the first column with powers of 0 or 1.
 * @param transposed a boolean to transpose the matrix before solving.
 * @param[out] x a point to the solution of the system; if x points to NULL then space is allocated for the result.
 * @param Pptr the Prime_ptr describing the finite field.
 *
 * @return a non-zero value iff the system was solved successfully. 
 *
 */
int solveFFVandermondeSystem(const elem_t* t, const elem_t* b, int size, int startExp, int transposed, elem_t** x, const Prime_ptr* Pptr);

/**
 * Solve a system of equations described by a (generalized) Vandermonde matrix.
 * The Vandermonde matrix is defined by the array of distinct values t such that
 * the ith row of the matrix is powers of t[i]. Such powers increase by 1 
 * per column starting from either 0 or 1, depending on the value of startExp.
 *
 * A non-zero value of startExp makes the first column t[i]^1, otherwise the
 * first column is t[i]^0. 
 *
 * If transposed is a non-zero value then the matrix is constructed as above 
 * but then transposed before solving.
 * 
 * @note: Currently only transposed = 1 are valid options.
 * @note: All finite field elements are assumed to be converted into montgomery form
 *        and the result is returned in montgomery form.
 *
 * @param t the distinct values defining the matrix.
 * @param b the right hand side vector.
 * @param size the size of the arrays t and b.
 * @param startExp a boolean describing to start the first column with powers of 0 or 1.
 * @param transposed a boolean to transpose the matrix before solving.
 * @param[out] x a point to the solution of the system; if x points to NULL then space is allocated for the result.
 * @param Pptr the Prime_ptr describing the finite field.
 *
 * @return a non-zero value iff the system was solved successfully. 
 *
 */
int solveFFVandermondeSystemInForm(const elem_t* t, const elem_t* b, int size, int startExp, int transposed, elem_t** x, const Prime_ptr* Pptr);

#ifdef __cplusplus
}
#endif


#endif
