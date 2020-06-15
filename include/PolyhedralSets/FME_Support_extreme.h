/*
 * This file provides an interface to the CDD
 * library.
 */
#ifndef EXTREME_H
#define EXTREME_H

#include "setoper.h"
#include "FME_Support_inequality.h"
#include "cdd.h"
#include "cddmp.h"
#include "FME_Support_linAlg.h"
#include <gmp.h>

/*
 * compute extreme rays of the polyhedron defined by 'coeffMat',
 * place the result in the 'extRays' matrix,
 * using CDD library.
 */

void findExtremeRays(matrix_t coeffMat , matrix_t * extrRays);

#endif
