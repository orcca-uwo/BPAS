/*
 * Implementation of inequality redundancy test, using Balas' algorithm.
 */

#ifndef REDUNDANCYTEST_H
#define REDUNDANCYTEST_H 
#include "FME_Support_inequality.h"
#include "FME_Support_linAlg.h"
#include "FME_Support_unrolledll.h"
#include "FME_Support_extreme.h"

void getCoeffMatrix(inequality_t * l , matrix_t * mat , int ineqNum);

void getConstVector(inequality_t * l , matrix_t * vec , int ineqNum);

void remVarMat(const matrix_t inputMat , const matrix_t inputVec , matrix_t * outMat , int varNum);

void elimVarMat(const matrix_t inputMat , const matrix_t inputVec , matrix_t * outMat , int varNum);

void blocElimination(const matrix_t inputMat , matrix_t * extrMat , matrix_t * outMat);

void initialTestCone(inequality_t * inputSys , int ineqNum , matrix_t * initTestCone);

void testCone(const matrix_t * initialCone , int step , matrix_t * output);

int extremeRayTest(matrix_t * coeffMat , matrix_t * vec);

int redundancyTest(matrix_t * testCone , inequality_t * toBeTestedIneq , int idx);
#endif

