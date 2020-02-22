/*
 * Implementation of minimal polyhedral projection using FME
 */
#ifndef PROJECTION_H
#define PRIJECTION_H

#include "FME_Support_fme.h"

/*
 * This function reads coefficients of a full-dimensional and pointed
 * polyhedron with 'varNum' varialbes and 'ineqNum' inequalities from
 * the 'file' and finds its minimal triangular form using FME.
 */

void project(char * file, int varNum, int ineqNum);

#endif
