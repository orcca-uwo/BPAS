#ifndef PROJECTION_DNC_H
#define PRIJECTION_DNC_H

#include "FME_Support_fme_dnc.h"

/*
 * This function reads coefficients of a full-dimensional and pointed
 * polyhedron with 'varNum' varialbes and 'ineqNum' inequalities from
 * the 'file' and finds its minimal triangular form using FME.
 */

void project_dnc(char * fileName, int varNum, int ineqNum, int thr);

#endif
