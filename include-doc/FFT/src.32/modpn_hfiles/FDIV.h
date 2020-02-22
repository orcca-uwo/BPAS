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
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


