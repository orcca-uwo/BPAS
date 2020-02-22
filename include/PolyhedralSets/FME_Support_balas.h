/*
 * FME_Support_balas.h
 * ****************************
 * implementation of funcitons to find the
 * redundancy test cone (W^0) of an input polyhedron.
 */

#ifndef BALAS_H
#define BALAS_H

#include "FME_Support_inequality.h"
#include "FME_Support_linAlg.h"
#include "FME_Support_fme.h"
#include "FME_Support_unrolledll.h"

/*
 * find coefficient matrix of eliminating variables.
 * input: the input inequality_t system
 * result: one-dimensional array containing coefficients
 * elimVarsNum: number of varialbes to be eliminated
 */

void whatIsA(inequalityULL_t input, mpq_t * result, int elimVarsNum);

//********************************************************************
/*
 * find remained variables coefficient matrix
 * input: the input inequality_t system
 * result: one-dimensional array containing coefficients
 * elimVarsNum: number of varialbes to be eliminated
 * ineqNum: number of inequalities in the input system
 * varNum: dimension of the polyhedron (number of varialbes in the input system)
 */

void whatIsB(inequalityULL_t input, mpq_t * result, int elimVarsNum, int ineqNum, int varNum);

//***************************************************************************************
/*
 * find right hand side vector
 * input: the input inequality_t system
 * result: one-dimensional array containing coefficients
 * ineqNum: number of inequalities in the input system.
 */

void whatIsd(inequalityULL_t input, mpq_t * result, int ineqNum);

//********************************************************************

/*
 * find B0 matrix from B.
 * The B0 matrix is used to find minimal representation of the projected polyhedron.
 * B is m * q and it is full columns rank
 * B: m*q, full-column rank matrix of remaining variables
 * B0: computed B0 matrix
 * p: number of eliminated variables
 */

void whatIsB0(mpq_t * B, mpq_t * B0, int m, int q, int p);

//********************************************************************

/*
 * find a variable indexed by v from an equation and substitute it in an
 * inequality_t system.
 * s1: the equation used for evaluating
 * s2: system of inequalities
 * v: index of a variable to be substitute
 * ineqNum: number of inequalities in the input system
 * varNum: number of variables in the input system
 */

void eval(inequality_t s1, inequalityULL_t s2, int v, int ineqNum, int varNum);

//***********************************************************************
/*
 * check if coefficient of a variable is zero in all inequalities
 * returns 1 if this is true, zero otherwise.
 * l: system of linear inequalities
 * c: index of the variable to check
 */
int zeroCol(inequalityULL_t l, int c);

//********************************************************************
int isUnit(node * current, int varNum, int j);

//********************************************************************
/*
 * swap two variables in a system
 * l: inequality_t system
 * c1: index of the first variable
 * c2: index of the second variable
 */

void swapVar(inequalityULL_t * l, int c1, int c2);

//********************************************************************

//
/*find and return the redundancy test cone, W0.
 * W0 = {omega (B0 ^ -1)A = 0, omega (B0 ^ -1) >= 0, omega (B0 ^ -1) d >=0 }
 * A0: coefficient vector of variables to eliminate (find by "whatIsA" function )
 * A0: coefficient vector of variables to remain (find by "whatIsB" function )
 * d0: right-hand-side vector (find by "whatIsd" function )
 * m: number of inequalities in the input system
 * q: number of remaining variables
 * p: number of eliminating variables
 * r: is equal to q (for future improvements)
 */
inequalityULL_t whatIsW0(mpq_t * A0, mpq_t * B0, mpq_t * d0, int m, int q, int p, int r);

//***************************************************************************************
/*
 * find coefficient matrix of an inequality_t system
 */
void coefMatrix(inequalityULL_t input, mpq_t * mat, int varNum, int ineqNum);

//********************************************************************
/*
 * put everything together to find redundancy test cone from input system
 * data: input polyhedron
 * elimNumber: number of variables to eliminate
 * varNum: number of variables in the input system
 * ineqNum: number of inequalities in the input system
 * W0ll: the resulting cone
 */

void balasW0(inequalityULL_t data, int elimNumber, int varNum, int ineqNum, inequalityULL_t * W0ll);

//*******************************************************************************************
/*
 * check if a vector is an extreme ray of a cone, using algebraic test
 * returns 1 if it is, 0 if it is not.
 * A: coefficient matrix of the cone
 * v: the vector to check
 * varNum: number of variables
 * ineqNum: number of inequalities
 * vcheck: an integer computed in another function
 */
int checkBoundary(mpq_t * A, mpq_t * v, int varNum, int ineqNum, int vcheck);

#endif
