/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __IsInvertible
#define __IsInvertible

#include "Types.h"
#include "generalFuncs.h"
#include "FINTERP.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include "LinkedList.h"
#include "IteratedResultant.h"
#include "RegularGcd.h"
#include <math.h>

LinkedQueue *
isInvertible_zeroDim(preFFTRep *poly, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);

#endif
