/*
 * Implementation of mathematical operations on inequalities
 */

#ifndef INEQUALITYOPERATIONS_H
#define INEQUALITYOPERATIONS_H

#include "FME_Support_inequality.h"
#include "FME_Support_linAlg.h"
#include "FME_Support_unrolledll.h"


int ineqSgn(inequality_t a, int varIndex);

void combineTwoIneq(inequality_t ineq1, inequality_t ineq2, inequality_t * result, int varIndex);

void findFMEMultiplier(inequality_t ineq1, inequality_t ineq2,
		int varIndex , mpz_t result1 , mpz_t result2);

void findCommenMult(mpz_t rational1, mpz_t rational2 , mpz_t mult1 , mpz_t mult2);

void multScalarToIneq(inequality_t * a, mpz_t mult, inequality_t * result);

void addTwoIneq(inequality_t ineq1, inequality_t ineq2,
		inequality_t * result);
		
void simplifyIneq(inequality_t * ineq1 , inequality_t * result);


#endif
