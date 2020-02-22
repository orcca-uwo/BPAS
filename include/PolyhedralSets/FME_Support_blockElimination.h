/*
 * FME_Support_blockElimination.h
 * *******************************
 * implementation of block elimination to get the
 * projection of the redundancy test cone, W0
 */
#ifndef BLOCKELIMINATION_H
#define BLOCKELIMINAITON_H

#include "FME_Support_extreme.h"
#include "FME_Support_linAlg.h"
#include "FME_Support_balas.h"
#include "FME_Support_unrolledll.h"

/**
 * find p' <= p
 * mat: coefficient matrix of input polyhedron
 * p: number of variables to eliminate
 * m: number of inequalities in the input system
 */
int pPrime(mpq_t * mat, int p, int m);
/***********************************************

/**
 * decompose the input matrix to get coefficient vector of the projection cone
 * mat: input matrix
 * out1: coefficient matrix of projection cone
 * out2: coefficient matrix of other variables
 * pprime: number less that p, computed in pPrime function
 * q: number of variables that remain
 * m: number of inequalities
 */
void projectionCone(mpq_t * mat, mpq_t * out1, mpq_t * out2, int pprime, int q, int m);

/***************************************************************************************

/**
 * use the structure of the cone and prepare date for finding extreme rays
 * mat: coefficient of the input cone
 * extrIn: input data prepared for finding extreme rays
 * m: number of inequalities
 * qprime: integer number related to pprime and q
 */

void makeExtrInput(mpq_t * mat, mpq_t * extrIn, int m, int qprime);

/*******************************************************************
/**
 * put everythings together to find and return projection of W0 cone, using block elimination
 * w0: the cone W0
 * p: number of variables to eliminate
 * q: number of variables that remain
 * m: number of inequalities
 * size: store number of inequalities in the projection
 */
mpq_t * blockElimination(inequalityULL_t w0, int p, int q, int m, int * size);
	
#endif

