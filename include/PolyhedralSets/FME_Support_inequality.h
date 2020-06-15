/*
 * Implementation of inequality strucutre and its related functions
 */

#ifndef INEQUALITY_H
#define INEQUALITY_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define MAX_LEN 500

/*
 * This data structure shows an inequality by its dimension,
 * variable coefficients and the right-hand-side constant.
 */

typedef struct
{
	int dimension;
	mpz_t * coeff;
	mpz_t constant;
} inequality_t;



int allocInequality(int varNum, inequality_t * newIneq);

void setInequality(inequality_t * newIneq , mpz_t * coeffData , mpz_t constantData);

void copyInequality(inequality_t * source , inequality_t * dest);

void printInequality(inequality_t * input , char var);

void getFromFile(mpz_t * data , char * fileName , int varNum , int ineqNum);

void rationalToInt(mpq_t * source , mpz_t * dest , int len);

void getInputArray(char * fileName , inequality_t * data , int varNum , int ineqNum);

void freeInequality(inequality_t * a);




#endif
