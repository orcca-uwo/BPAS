/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __SubResultantSeq_H
# define __SubResultantSeq_H

#include "Types.h"
#include "generalFuncs.h"
#include "MPMMTS.h"
#include "FDIV.h"

sfixn *
SubResultantSeq(sfixn dg1, sfixn *f1, sfixn dg2, sfixn *f2, MONTP_OPT2_AS_GENE *pPtr);

sfixn *
SubResultantSeq_1(sfixn w, sfixn Ssz, sfixn *S, sfixn dg1, sfixn *f1, sfixn dg2, sfixn *f2, MONTP_OPT2_AS_GENE *pPtr);


sfixn *
SubResultantSeq_1_new(sfixn w, sfixn Ssz, sfixn *S, sfixn dg1, sfixn *f1, sfixn dg2, sfixn *f2, MONTP_OPT2_AS_GENE *pPtr);

void
printSRS(sfixn w, sfixn *SRS);

sfixn EX_Resultant_Uni(sfixn dgP, sfixn *P, sfixn dgQ, sfixn *Q, MONTP_OPT2_AS_GENE *pPtr);


#endif 
