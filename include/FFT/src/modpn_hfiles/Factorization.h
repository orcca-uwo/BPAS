/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __Factorization_h
#define __Factorization_h 

#include "Types.h"
#include "generalFuncs.h"
#include "MPMMTS.h"
#include "FINTERP.h"




sfixn *
LcmPolyPair(sfixn *dgLcmAddr, sfixn d1, sfixn *f1, sfixn d2, sfixn *f2, MONTP_OPT2_AS_GENE * pPtr );

sfixn *
SquareFreeFact(sfixn *degR, sfixn degF, sfixn *FPtr, sfixn p);

#endif
