/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __matrix_h
#define __matrix_h 

#include "Types.h"
#include "UniHensel.h"

#include <stdlib.h>            /* rand() and srand() functions               */
#include <assert.h>            /* assert()                                   */



void multi_mat_transpose (sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data);

void mat_transpose (sfixn, sfixn, sfixn * );

POLYVECTOR_SLG *randomPolyVec_SLG(int32 M, int32 GN, sfixn p, int32 N);

void freeVec_SLG(POLYVECTOR_SLG *vec_slg);

void  printJMatrix(POLYMATRIX * mat);

void freeJMatrix(POLYMATRIX * mat);

POLYMATRIX *createJMatrix(POLYVECTOR_SLG *polyVec_SLG, TriSet *ts, TriRevInvSet * tris, sfixn N, MONTP_OPT2_AS_GENE * pPtr);


POLYMATRIX *initJMatrix(sfixn N, int32 m, int32 n, TriSet *ts);


POLYMATRIX *mulJMatrix(sfixn N, POLYMATRIX * mat1, POLYMATRIX * mat2, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
randomJMatrix(POLYMATRIX * mat, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
INVJMatrix(sfixn N, POLYMATRIX * mat, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
createJMatrix_ForLifting(POLYVECTOR_SLG *polyVec_SLG, TriSet *ts, TriRevInvSet * tris, sfixn N, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
SLG2PolyVecotr(POLYVECTOR_SLG *polyVec_SLG, TriSet *ts, TriRevInvSet * tris, sfixn N, MONTP_OPT2_AS_GENE * pPtr);



POLYMATRIX *
createJMatrix_PolyForlifting(TriSet *ts, TriRevInvSet * tris,  sfixn N, MONTP_OPT2_AS_GENE * pPtr);


POLYMATRIX *
scalarMulJMatrix_1(sfixn r, POLYMATRIX * mat, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
subJMatrix_1(sfixn N, POLYMATRIX * mat1, POLYMATRIX * mat2, MONTP_OPT2_AS_GENE * pPtr) ;


TriSet *
AddHs2TriSet(sfixn N, TriSet *ts, POLYMATRIX *Hs, MONTP_OPT2_AS_GENE * pPtr);



POLYMATRIX *
increaseMatrix_ForLifting(sfixn N, POLYMATRIX *mat, TriSet *ts);


void *BlockDecompose(void * PTR);



POLYVECTOR_SLG * example_1_PolyVec_SLG();


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


