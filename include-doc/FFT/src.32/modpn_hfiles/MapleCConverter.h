#ifndef __MapleCConverter_h
#define __MapleCConverter_h 

#include "Types.h"
#include "generalFuncs.h"
#include "UniHensel.h"
#include "MPMMTS.h"
#include "FINTERP.h"
#include "LinkedList.h"
#include "IteratedResultant.h"
#include "IsInvertible.h"
//#include "/home/xli96/maple10/extern/include/maplec.h"



preFFTRep * createOneWrapperPoly(sfixn  N, sfixn *dgs, sfixn *data);

SLG * Maple2C_DAG2DAG(int GN, sfixn *MDA);

void create_pdeg_coef_Vec (sfixn *pdegVec, sfixn*coefVec, preFFTRep * poly);

preFFTRep *inverse_create_pdeg_coef_Vec(sfixn N, sfixn *dgs, sfixn *pdegVec, sfixn*coefVec);


//ALGEB M_DECL MyIdentity( MKernelVector kv, ALGEB *args );

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
TestMDag2CDag(sfixn GN,  sfixn *MDA);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall 
#else
void
#endif
TestRecden2C(int N, sfixn  sz, sfixn *dgs,  sfixn *MDA);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
TestC2Recden ( sfixn pdVdeg, sfixn *pdegVec, sfixn cVdeg, sfixn *coefVec);



void getSMPfromC (sfixn size, sfixn *buffer);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
MulPolyTFTFFTCN(sfixn N, sfixn *rdgs, sfixn rBsz, sfixn *resBuffer, 
		sfixn *p1dgs, sfixn p1sz, sfixn *p1Buffer, 
                sfixn *p2dgs, sfixn p2sz, sfixn *p2Buffer, 
                sfixn dVsz, sfixn *pdegVec,
                sfixn cVsz, sfixn *coefVec,
                sfixn p);
int
estimatePartialDegVecSize(sfixn *dgs, sfixn n);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
MulPolyTFTFFTCNC(sfixn N,  sfixn *dgs1,
		 sfixn p1dgssz, sfixn *p1dgs, sfixn p1sz, sfixn *p1Buffer, 
		     sfixn *dgs2,
		 sfixn p2dgssz, sfixn *p2dgs, sfixn p2sz, sfixn *p2Buffer, 
		     sfixn *rdgs,
                sfixn dVsz, sfixn *pdegVec,
                sfixn cVsz, sfixn *coefVec,
		 sfixn p);



#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
 TFTFFTUNIC(sfixn dr, sfixn *resPtr,
		  sfixn d1, sfixn *v1Ptr, sfixn d2, sfixn *v2Ptr,
                  sfixn p);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
FASTDIVC(sfixn degR, sfixn *RPtr, sfixn degQ, sfixn *QPtr, sfixn degA, sfixn *APtr, sfixn degB, sfixn *BPtr, sfixn p);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
PLAINDIVC(sfixn degR, sfixn *RPtr,sfixn degQ, sfixn *QPtr, sfixn degA, sfixn *APtr, sfixn degB, sfixn *BPtr, sfixn p);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
PLAINGCDUNIC(sfixn ud, sfixn *uPtr, sfixn vd, sfixn *vPtr, sfixn gd, sfixn *gcdPtr, sfixn dA, sfixn *APtr, sfixn dB, sfixn *BPtr, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
FASTGCDUNIC(sfixn ud, sfixn *uPtr, sfixn vd, sfixn *vPtr, sfixn gd, sfixn *gcdPtr, sfixn dA, sfixn *APtr, sfixn dB, sfixn *BPtr, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall 
#else
void
#endif
subProdTreeCreWrapC(sfixn h, sfixn levels, sfixn *W, sfixn *NoNodes, 
                   sfixn *Bases, sfixn totSZ, sfixn *data, 
		    sfixn itemNo, sfixn itemSz, sfixn p);

void
subProdTreeFreeWrapC(subProdTree * tree);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
FastEvalWrapC(sfixn n, sfixn *EvalPts, sfixn degf, sfixn *fPtr, sfixn h, sfixn *W, sfixn *NoNodes, sfixn *Bases, sfixn *data, sfixn p);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
FastInterpWrapC(sfixn n, sfixn *InterpedPts, sfixn *EvaluatingPts, sfixn *EvaluatedPts, sfixn h, sfixn *W, sfixn *NoNodes, sfixn *Bases, sfixn *data, sfixn p);


// Creating the pts_tree, and return the result data back to Maple.




#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
PTRTREESCRECN(sfixn *Nmp, sfixn *ptPHWDSZ, sfixn *bounds, sfixn *dims1, sfixn *dims2,
	    sfixn *pts_s,
            sfixn *h_s,
            sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
            sfixn *data_s,
	    sfixn *p1dgs, sfixn p1sz, sfixn *p1Buffer, 
	      sfixn *p2dgs, sfixn p2sz, sfixn *p2Buffer);


// Creating the pts_tree by using the data passed from Maple.
PTS_TREE *
createWrapperPTS_TREE(sfixn N, sfixn m, sfixn *bounds, sfixn pts_sSz, sfixn *pts_s,
                          sfixn h_sSz, sfixn *h_s,
                          sfixn WNB_sSz, sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
		      sfixn data_sSz, sfixn *data_s, sfixn p);


subProdTree *
createWrapperTree(sfixn h, sfixn *W, sfixn *NoNodes, sfixn *Bases, sfixn *data);

void freeWrapperPTS_TREE(PTS_TREE *pts_tree);

#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
FastEvalMultiWrapCN(sfixn *Nmp, sfixn *ptPHWDSZ, 
                   sfixn *bounds, 
	           sfixn *pts_s,
                   sfixn *h_s,
                   sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
                   sfixn *data_s,
                   sfixn *dims, sfixn Esz, sfixn *E,
		    sfixn *fdgs, sfixn fsz, sfixn *fBuffer);




#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
FastInterpMultiWrapCN(sfixn *Nmp, sfixn *ptPHWDSZ,
                      sfixn *bounds, 
	              sfixn *pts_s, sfixn *h_s,
                      sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
                      sfixn *data_s,
                      sfixn *dims, sfixn Esz, sfixn *E,
                      sfixn dVsz, sfixn *pdegVec,
                      sfixn cVsz, sfixn *coefVec);



#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
SubResultantChains(sfixn N, sfixn w, sfixn Ssz, sfixn *S, sfixn* Edims1,                      sfixn E1sz, sfixn *E1, sfixn *Edims2, sfixn E2sz, sfixn*E2, sfixn p);







#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
InterpIthDthMultiWrapCN(sfixn *Nmp, sfixn *ptPHWDSZ, sfixn ith, sfixn dth, sfixn *bounds, 
	           sfixn *pts_s, sfixn *h_s,
                   sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
                   sfixn *data_s,
		   sfixn w, sfixn subslicesz, sfixn *subslicedims, 
                   sfixn Ssz, sfixn *S,
		      //sfixn *fdgs, sfixn fsz, sfixn *fBuffer,
                   sfixn dVsz, sfixn *pdegVec,
                   sfixn cVsz, sfixn *coefVec);



#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
InterpIthMultiWrapCN(sfixn *Nmp, sfixn *ptPHWDSZ, sfixn ith, sfixn *bounds, 
	           sfixn *pts_s, sfixn *h_s,
                   sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
                   sfixn *data_s,
		   sfixn w, sfixn slicesz, sfixn *slicedims, 
                   sfixn Ssz, sfixn *S,
		      //sfixn *fdgs, sfixn fsz, sfixn *fBuffer,
                   sfixn dVsz, sfixn *pdegVec,
		     sfixn cVsz, sfixn *coefVec);



#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
InterpNextLTMultiWrapCN(sfixn *Nmp, sfixn *ptPHWDSZ,sfixn start, sfixn *bounds, 
	           sfixn *pts_s, sfixn *h_s,
                   sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
                   sfixn *data_s,
		   sfixn w, sfixn subslicesz, sfixn *subslicedims, 
                   sfixn Ssz, sfixn *S,
		      //sfixn *fdgs, sfixn fsz, sfixn *fBuffer,
                   sfixn dVsz, sfixn *pdegVec,
			sfixn cVsz, sfixn *coefVec);


#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
InterpNextLCMultiWrapCN(sfixn *Nmp, sfixn *ptPHWDSZ,sfixn start, sfixn *bounds, 
	           sfixn *pts_s, sfixn *h_s,
                   sfixn *W_s, sfixn *NoNodes_s, sfixn *Bases_s,
                   sfixn *data_s,
		   sfixn w, sfixn subslicesz, sfixn *subslicedims, 
                   sfixn Ssz, sfixn *S,
		      //sfixn *fdgs, sfixn fsz, sfixn *fBuffer,
                   sfixn dVsz, sfixn *pdegVec,
			sfixn cVsz, sfixn *coefVec);





#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
DftMultiWrapCN(sfixn *Nmp,
               sfixn *es, sfixn *dims, sfixn Esz, sfixn *E,
	       sfixn *fdgs, sfixn fsz, sfixn *fBuffer, sfixn *rootsPtr);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
InvDftMultiWrapCN(sfixn *Nmp,
		  sfixn *es, sfixn *dims, sfixn Esz, sfixn *E,
                      sfixn dVsz, sfixn *pdegVec,
		  sfixn cVsz, sfixn *coefVec, sfixn *rootsPtr);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
InvDftIthDthMultiWrapCN(sfixn *Nmp, sfixn ith, sfixn dth, 
 		   sfixn w, sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims, 
                   sfixn Ssz, sfixn *S,
                   sfixn dVsz, sfixn *pdegVec,
			sfixn cVsz, sfixn *coefVec, sfixn *rootsPtr);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
InvDftIthMultiWrapCN(sfixn *Nmp, sfixn ith, 
		     sfixn w, sfixn slicesz, sfixn *slicees, sfixn *slicedims, 
                   sfixn Ssz, sfixn *S,
                   sfixn dVsz, sfixn *pdegVec,
		     sfixn cVsz, sfixn *coefVec, sfixn *rootsPtr);

#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
InvDftNextLCMultiWrapCN(sfixn *Nmp, sfixn start, 
			sfixn w, sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims, 
                   sfixn Ssz, sfixn *S,
                   sfixn dVsz, sfixn *pdegVec,
			sfixn cVsz, sfixn *coefVec, sfixn *rootsPtr);

#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
createGoodRootsCN(sfixn N, sfixn M, sfixn *f1dgs, sfixn *f1Buffer, 
                  sfixn *f2dgs, sfixn *f2Buffer,
		  sfixn *es, sfixn *dims, sfixn *rootsPtr, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall 
#else
void
#endif
PLAINRFRUNIC(sfixn d, sfixn vd, sfixn *vPtr, sfixn gd, sfixn *gcdPtr, sfixn dA, sfixn *APtr, sfixn dB, sfixn *BPtr, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) int __stdcall
#else
int
#endif
NewtonLiftUniCN(sfixn *outPDGVECS, sfixn *outCOEFVECS, sfixn Y, sfixn y0, sfixn N,
                sfixn *GNS, sfixn *MDAS, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
MultiModCN(sfixn dVsz, sfixn *pdegVec, sfixn cVsz, sfixn *coefVec, sfixn N, sfixn *fdgs, sfixn *fBuffer, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p, sfixn opt);


#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
NormalizeCN(sfixn *outPDGVECS, sfixn *outCOEFVECS, sfixn N, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
GetQuotientCN(sfixn N, sfixn dd, sfixn Ssz, sfixn *S, sfixn* Edims1,                      sfixn E1sz, sfixn *E1, sfixn *Edims2, sfixn E2sz, sfixn*E2, sfixn p, sfixn opt);

#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
ReduceCoeffCN(sfixn dVsz, sfixn *pdegVec, sfixn cVsz, sfixn *coefVec, sfixn N, sfixn *fdgs, sfixn *fBuffer, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p);

#ifdef WINDOWS
__declspec(dllexport) void __stdcall
#else
void
#endif
QuotientModTriSetCN(sfixn dVsz, sfixn *pdegVec, sfixn cVsz, sfixn *coefVec, sfixn N, sfixn *fdgs1, sfixn *fBuffer1, sfixn *fdgs2, sfixn *fBuffer2, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p);



#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
FastInterpRFRWrapC(sfixn *np, sfixn *InterpedPtsNum, sfixn *InterpedPtsDen, sfixn *EvaluatingPts, sfixn *EvaluatedPts, sfixn h, sfixn *W, sfixn *NoNodes, sfixn *Bases, sfixn *data);


#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
IterResOneDimCN(sfixn *outVec, sfixn M, sfixn *fdgs, sfixn *fBuffer, sfixn N, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn bound, sfixn freeVarNo, sfixn p);

#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
IterResZeroDimCN(sfixn M, sfixn *fdgs, sfixn *fBuffer, sfixn N, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS,  sfixn p);

#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif 
isInvertableCN(sfixn N, sfixn *fdgs, sfixn *fBuffer, sfixn *TS_DGS,
               sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p);


#ifdef WINDOWS
__declspec(dllexport) sfixn __stdcall
#else
sfixn
#endif
IsInvertibleChainCN(sfixn *Ns, sfixn *outPolyPDGVECS, sfixn *outPolyCOEFVECS, sfixn *outTsPDGVECS, sfixn *outTsCOEFVECS, sfixn N, sfixn *fdgs, sfixn *fBuffer, sfixn *TS_DGS, sfixn *inDGS, sfixn *inSIZS, sfixn *inCOEFS, sfixn p);

TriSet *
createWrapperDenseTriSet(sfixn N, sfixn * dgs);



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


