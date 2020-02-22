/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __Iter_h
#define __Iter_h 

#include "Types.h"
#include "generalFuncs.h"
#include "FINTERP.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include <math.h>


sfixn
iteratedResultant_zerodim(preFFTRep *poly, TriSet *ts, MONTP_OPT2_AS_GENE *pPtr);


sfixn *iteratedResultant_onedim(sfixn *resDgAddr, preFFTRep *poly, 
                         TriSet *ts, sfixn bound, sfixn freeVarNo, 
                         MONTP_OPT2_AS_GENE * pPtr);

preFFTRep *EX_Resultant_Multi(preFFTRep *f1, preFFTRep *f2, sfixn N, MONTP_OPT2_AS_GENE *pPtr);

SCUBE *EX_SubResultantChain(preFFTRep *f1, preFFTRep *f2, sfixn N, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *EX_ResultantFromChain(SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr);

SCUBE *EX_SCUBE_Init(sfixn M, sfixn *bounds, sfixn w);

void EX_SCUBE_Free(SCUBE *scube);

void EX_SCUBE_Print(SCUBE *scube);

preFFTRep *
EX_Resultant_Multi_Wrapper(preFFTRep *f1, preFFTRep *f2, sfixn N, MONTP_OPT2_AS_GENE *pPtr);


#endif
