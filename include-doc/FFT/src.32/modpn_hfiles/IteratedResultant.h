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


