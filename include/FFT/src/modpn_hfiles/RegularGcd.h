/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __RegularGcd
#define __RegularGcd

#include "Types.h"
#include "generalFuncs.h"
#include "FINTERP.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include "LinkedList.h"
#include "IteratedResultant.h"
#include "IsInvertible.h"
#include <math.h>

LinkedQueue *EX_RegularizeList_1(LinkedQueue *RegQueue, LinkedQueue *ToCheckQueue, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *
EX_RegularGcd(preFFTRep *f1, preFFTRep *f2, TriSet *ts,  SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *
EX_RegularGcd_Wrapped(preFFTRep *f1, preFFTRep *f2, TriSet *ts,  sfixn M, MONTP_OPT2_AS_GENE *pPtr);

#endif
