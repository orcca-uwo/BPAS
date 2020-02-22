/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __FDIV_h
#define __FDIV_h

#include "Types.h"
#include "generalFuncs.h"
#include "FMUL.h"


sfixn *
modularInvPM(sfixn, sfixn *, sfixn, sfixn *, sfixn, sfixn, MONTP_OPT2_AS_GENE *);

void 
fmedg_1(sfixn degRes, sfixn * resPtr, sfixn e, sfixn r, sfixn degp2, sfixn * p2Ptr, MONTP_OPT2_AS_GENE *);


void 
fastDiv(sfixn *, sfixn, sfixn *, sfixn, sfixn *, sfixn, sfixn *, MONTP_OPT2_AS_GENE *);

void 
plainDiv(sfixn * RPtr,sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

double
fastDiv_bench(sfixn * RPtr, sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );


void 
plainDivMonic_1(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE *);

void 
plainDiv_1(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr,  MONTP_OPT2_AS_GENE *);

void 
fastDiv_1(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, sfixn * BRevInvPtr, MONTP_OPT2_AS_GENE * pPtr );


void 
plainRemMonic_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void
plainDivNew(sfixn * RPtr,sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void 
plainRem(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void 
fastQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void 
plainQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void 
fastRem(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );


void 
UniQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void 
UniRem(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );


sfixn *
EX_UniRem(sfixn *degRAddr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );


sfixn *
EX_UniQuo(sfixn *degQAddr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void
PlainPseudoRemainder(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);


void
PlainPseudoQuotient(sfixn *degQAddr, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

void 
UniPseuQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );


sfixn *
EX_PQuo_Uni(sfixn *degQAddr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr);

#endif
