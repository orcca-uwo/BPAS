/**
 * FME_Support_kohler.h
 * ***********************
 * implementation of checks for redundancies
 */


#ifndef KOHLER_H
#define KOHLER_H

#include "FME_Support_inequality.h"
#include "FME_Support_linAlg.h"
#include "FME_Support_balas.h"
#include "FME_Support_unrolledll.h"


/**
 * checking if an inequality_t pass kohler's check w.r.t. a system
 * returns 0 if does not pass, 1 otherwise
 * a: inequality_t to check
 * d: inequality_t system
 * varNum: number of variables
 * ineqNum: number of inequalities
 */

int kohlerCheck(inequality_t a, inequalityULL_t * d, int v , int varNum , int ineqNum);
//*******************************************************************************

/**
 * balas check function returns 1 if inequality_t is not redundant, 0 otherwise
 * a: inequality_t to check
 * pw0: projection of redundancy test cone
 * varNum: number of variables
 * ineqNum: number of inequalities
 * sizepw0: number of inequalities in pw0
 * elimNum: number of eliminated variables
 */

int balasCheck(inequality_t a, mpq_t * pw0, int varNum, int ineqNum , int sizepw0 , int elimNum);


#endif
