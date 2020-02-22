#ifndef __FINTERP_h
#define __FINTERP_h 

#include "Types.h"
#include "generalFuncs.h"
#include "FMUL.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include <math.h>


subProdTree *
subProdTreeCre( sfixn itemNo, sfixn itemSz, sfixn *items, sfixn p);

void
subProdTreeFree(subProdTree * tree);
void
printSubProdTree(subProdTree * tree);
sfixn *
FastEvaluation(sfixn n, sfixn degf, sfixn *fPtr, subProdTree * tree, sfixn p);
sfixn *
SlowEvaluation(sfixn degf, sfixn *fPtr, sfixn nopts, sfixn *pts, sfixn p);
sfixn *
linearCombineModulus(sfixn *Cs, subProdTree * tree, sfixn p);
sfixn *
fastInterp(sfixn* polyDg, sfixn n, sfixn *Us,  subProdTree * tree, sfixn *Vs, sfixn p);
sfixn *
direvative(sfixn deg, sfixn *coef, sfixn p);
void
mulNodes(sfixn dr, sfixn *srcAddr, sfixn ds, sfixn *Addr1, sfixn *Addr2, sfixn nodeSz, MONTP_OPT2_AS_GENE * pPtr);
sfixn *
createPts(sfixn start, sfixn m, sfixn *bounds, sfixn p);
subProdTree**
createArrOfSPTrees(sfixn m, sfixn *bounds, sfixn *pts, sfixn p);
void
freeArrOfSPTTrees(sfixn m, subProdTree **trees);
void
printArrOfSPTTrees(sfixn m, subProdTree **trees);

void
FastEvaluation_1(sfixn resultsz, sfixn *result, sfixn degf, sfixn *fPtr, subProdTree * tree, sfixn p);

void
freePTS_TREE(PTS_TREE *pt);


void fastEvalMulti_test(sfixn M, sfixn *dims, sfixn *E, preFFTRep* poly, subProdTree** trees, MONTP_OPT2_AS_GENE * pPtr);




preFFTRep*
fastInterpMulti_test(sfixn N, sfixn M, sfixn *dims, sfixn *EEE, sfixn **UsPtr, subProdTree** trees, MONTP_OPT2_AS_GENE * pPtr);

sfixn **
convertpts2ptsPtr(sfixn m, sfixn *bounds, sfixn *pts);

PTS_TREE*
createGoodPtsForf1f2(sfixn N, sfixn d1, preFFTRep *f1, sfixn d2, preFFTRep *f2, sfixn m, sfixn *bounds, sfixn *dims1, sfixn *dims2, MONTP_OPT2_AS_GENE *pPtr);


int
createGoodRootsForf1f2(sfixn N, sfixn M, sfixn d1, preFFTRep *f1, sfixn d2, preFFTRep *f2, sfixn *es, sfixn *dims, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

void
getSubResultantChains(sfixn N, sfixn w, sfixn Ssz, sfixn *S, sfixn* Edimsq,                      sfixn *E1, sfixn *Edims2, sfixn*E2, MONTP_OPT2_AS_GENE * pPtr);

void
printAllSRS(sfixn no, sfixn w, sfixn *AllSRS);


sfixn *
get_ithSlice_fromSubResultantChains(sfixn ith, sfixn N, sfixn w, sfixn slicesz, sfixn *slicedims, sfixn Ssz, sfixn *S);


sfixn *
get_ithDthSubSlice_fromSubResultantChains(sfixn ith, sfixn dth, sfixn N_1, sfixn w, sfixn subslicesz, sfixn *subslicedims, sfixn Ssz, sfixn *S);


int
tracingNextCandidateSlice(int start, sfixn w, sfixn Ssz, sfixn *S);


void
set_ithSliceZero_fromSubResultantChains(sfixn ith, sfixn N, sfixn w, sfixn slicesz, sfixn *slicedims, sfixn Ssz, sfixn *S);


preFFTRep*
interpIthSlice(sfixn ith, sfixn N, sfixn m, sfixn w, sfixn slicesz, sfixn *slicedims, sfixn Ssz, sfixn *S,  PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE * pPtr);


preFFTRep*
interpIthDthSlice(sfixn ith, sfixn dth, sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicedims, sfixn Ssz, sfixn *S,  PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep*
interpNextCandidateSliceLC(int*nextiAddr, int start, sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicedims, sfixn Ssz, sfixn *S,  PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep*
interpNextCandidateSliceLT(int*nextiAddr, int start, sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicedims, sfixn Ssz, sfixn *S,  PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE * pPtr);


void fastDftMulti_test(sfixn M, sfixn *es, sfixn *dims, sfixn *E, preFFTRep* poly,  sfixn *rootsPtr, MONTP_OPT2_AS_GENE * pPtr);

preFFTRep*
fastInvDftMulti_test(sfixn N, sfixn M, sfixn *es, sfixn *dims, sfixn *EEE, sfixn *rootsPtr,MONTP_OPT2_AS_GENE * pPtr);



preFFTRep*
interpIthSliceDFT(sfixn ith, sfixn N, sfixn m, sfixn w, sfixn slicesz, sfixn *slicees, sfixn *slicedims, sfixn Ssz, sfixn *S,  sfixn *rootsPtr, MONTP_OPT2_AS_GENE * pPtr);


preFFTRep*
interpIthDthSliceDFT(sfixn ith, sfixn dth, sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims, sfixn Ssz, sfixn *S,  sfixn *rootsPtr, MONTP_OPT2_AS_GENE * pPtr);


preFFTRep*
interpNextCandidateSliceLCDFT(int *nextiAddr, int start, sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims, sfixn Ssz, sfixn *S, sfixn *rootsPtr, MONTP_OPT2_AS_GENE * pPtr);


sfixn
SlowEvaluation1pt(sfixn degf, sfixn *fPtr, sfixn pt, sfixn p);


void
getQuotients(sfixn N, sfixn dd, sfixn Qsz, sfixn *Q, sfixn* Edims1,                      sfixn *E1, sfixn *Edims2, sfixn*E2, MONTP_OPT2_AS_GENE * pPtr, int opt );

void permuteSlice1toN(sfixn N, sfixn slicesz, sfixn *slicedims, sfixn *slice);

preFFTRep *EX_QuoMulti(preFFTRep *f1, preFFTRep *f2, sfixn N, MONTP_OPT2_AS_GENE * pPtr, int opt);



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


