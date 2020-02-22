/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __GCD_h
#define __GCD_h 

#include "Types.h"
#include "generalFuncs.h"
#include "FMUL.h"
#include "FDIV.h"




sfixn
gcd_Uni_1(sfixn * change, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr,MONTP_OPT2_AS_GENE * pPtr, sfixn );

sfixn
ExGcd_Uni_1(sfixn * change,sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr,MONTP_OPT2_AS_GENE * pPtr, sfixn * degS, sfixn * degT, sfixn * S1Ptr, sfixn * S2Ptr,sfixn * T1Ptr, sfixn * T2Ptr, sfixn);

void 
plainDiv(sfixn * RPtr, sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void
ExGcd_Uni(sfixn * uPtr, sfixn *dC, sfixn * vPtr, sfixn *dD, sfixn * gcdPtr, sfixn *dG, 
	  sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);

void
Gcd_Uni(sfixn * gcdPtr, sfixn *dG, 
	sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);


sfixn *
EX_GCD_UNI(sfixn *dGAddr, 
	   sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);


int32
ExGcd_Uni_RFR(sfixn d, sfixn * vPtr, sfixn *dD, sfixn * gcdPtr, sfixn *dG,  sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);

void
normalize_1(sfixn deg, sfixn * cof, MONTP_OPT2_AS_GENE * pPtr);


#endif
