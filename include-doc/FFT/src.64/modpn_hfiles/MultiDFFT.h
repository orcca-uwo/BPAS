/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __MultiDFFT_h
#define __MultiDFFT_h 

#include "Types.h"
#include "FMUL.h"
#include "matrix.h"
#include "generalFuncs.h"
#include <stdlib.h>            /* rand() and srand() functions               */
#include <assert.h>            /* assert()                                   */



// dgs1[0], dgs2[0] are useless.
// ccumx[i] keeps the base before ith-dimension.
void plainMultiDMul(sfixn N, sfixn * ccum, sfixn * res, sfixn * ccum1, sfixn * dgs1, sfixn * ccum2, sfixn * dgs2, sfixn * coeffs1, sfixn * coeffs2, MONTP_OPT2_AS_GENE * pPtr);


void fftMultiD_square_test(sfixn * coeffs1, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);

void fftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);



void fromtofftRepMultiD(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs);

void fftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, MONTP_OPT2_AS_GENE * pPtr);

void tftMultiD_test(sfixn * coeffs1, sfixn * coeffs2, sfixn N, sfixn * es, sfixn * dims, sfixn * ls,  MONTP_OPT2_AS_GENE * pPtr);




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


