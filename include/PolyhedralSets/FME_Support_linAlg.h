
/*
 * This file provides interface to the FLINT library for
 * linear algebra operations.
 */

#ifndef LINALG_H
#define LINALG_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

/*
 * This data structure shows matrix using a list of its coefficients,
 * number of rows and number of columns.
 */

typedef struct
{
	mpq_t * entries;
	int rowNum;
	int colNum;
} matrix_t;


void allocMatrix(int r , int c , matrix_t * m);

void setMatrix(matrix_t * m ,  mpq_t * data);

void printMatrix(matrix_t * m);

void matrixToFlint(matrix_t source , fmpq_mat_t * dest);

void flintToMatrix(fmpq_mat_t source , matrix_t * dest);

void matrixMatrixMult(matrix_t * mult1 , matrix_t * mult2 , matrix_t * out);

int matrixRank(matrix_t * input);

void matrixTranspose(matrix_t * source , matrix_t * dest);

void matrixInverse(matrix_t * mat , matrix_t * inverse);

void subMat(matrix_t * mat , int * index , int size , matrix_t * subMat);

void freeMatrix(matrix_t * m);

#endif
