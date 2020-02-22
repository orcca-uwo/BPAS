// -*- C++ -*- 
// DDMP.h
// distributed dense multivariate polynomial

#ifndef __DDMP_h
#define __DDMP_h

#include "modpn.h"

#include "ser_general_routine.h"
#include "ser_basic_routine.h"
#include "general_routine.h"
#include "basic_routine.h"

#include "PrimeField.h"

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef int sfixn;


//DDMP over Z/pZ, p--prime
template <class Field>
class DDMP 
{
private:
  Field* mField;
  preFFTRep *mPoly;
  //static const sfixn DENSECUTOFF = 42; //dense size
  //static const sfixn PARTDEGCUTOFF = 8; //partial degree
  //static const sfixn KNFFTCUTOFF = 4; //Kronecker univariate fft size
                       // why so small??????
  //data struct suitable for FFT/TFT, encoding:
  //N(poly)--number of variables, 
  //DAT(poly)--coefficient vector
  //BUSZS(poly)--partial degree vector
  //CUM(poly)--displacement vector
  //...see Types.h in modpn.h

public:
  DDMP(Field* F): mField (F), mPoly (NULL) {}
  //DDMP(sfixn N, sfixn *pdgs, sfixn *coeffs);

  //why do not need to free mPoly? it uses calloc!
  //~DDMP() { }
  ~DDMP() { 
    if (BUSZS(mPoly) != NULL) {
      my_free(BUSZS(mPoly)); 
      BUSZS(mPoly) = NULL;
    }    
    if (CUTS(mPoly) != NULL) {
      my_free(CUTS(mPoly)); 
      CUTS(mPoly) = NULL;
    }
    if (CUM(mPoly) != NULL) {
      my_free(CUM(mPoly)); 
      CUM(mPoly) = NULL;
    }
    if (DAT(mPoly) != NULL) {
      my_free(DAT(mPoly)); 
      DAT(mPoly) = NULL;
    }
    if (mPoly != NULL) {
      my_free(mPoly); 
      mPoly = NULL;
    }
  }


  preFFTRep* FFTRep() { return mPoly; }

  sfixn NumVars() { return N(mPoly); }

  sfixn* PartialDegrees() { return BUSZS(mPoly); }

  sfixn* coeffs() { return DAT(mPoly); }

  sfixn coeffsSize() { return SIZ(mPoly); }

  Field* CoeffField() { return mField; }

  /**---------------------------------------------------
   * getDenseSize
   * Data size of the C-Cube polynomial 'poly' minus all leading zeros
   * in it.
   **/
  sfixn DenseSize(){
    sfixn N = N(mPoly);
    return PBPAS::getDenseSiz(N, mPoly, BUSZSI(mPoly, N), DAT(mPoly), CUMI(mPoly, N));  
  } 
 
  void setFFTRep(preFFTRep* rep) { mPoly=rep; }

  /**
   * Create data struct for this poly and
   * Set the number of variables, partial degrees and coefficients
   * @param N Number of variables
   * @param pdgs Partial degree vector
   * @param coeffs Coefficient vector (in a prime field)
   **/
  void setData(int N, sfixn* pdgs, sfixn* coeffs){
    register int j;

    mPoly = (preFFTRep *)my_calloc(1, sizeof(preFFTRep));
    N(mPoly) = N;
    //BUSZS(mPoly) = (sfixn * )my_calloc(N+1, sizeof(sfixn));
    //CUTS(mPoly) = (sfixn * )my_calloc(N+1, sizeof(sfixn));
    SIZ(mPoly) = 1;
    CUM(mPoly) = (sfixn * )my_calloc(N+1, sizeof(sfixn));
    CUMI(mPoly, 1) = 1;

    BUSZS(mPoly) = pdgs; //point to
    CUTS(mPoly) = pdgs;

    for(j=1; j<=N; j++){
      //BUSZSI(mPoly, j) = pdgs[j-1];
      //CUTSI (mPoly, j) = pdgs[j-1];   
      SIZ(mPoly) = SIZ(mPoly) * (BUSZSI(mPoly, j) + 1);
      if(j >= 2){
	CUMI(mPoly, j) = CUMI(mPoly, j-1) * (BUSZSI(mPoly, j-1) + 1);
      }
    }
    OFST(mPoly) = 0;
    //DAT(mPoly) = (sfixn * )my_calloc( SIZ(mPoly),sizeof(sfixn));  
    DAT(mPoly) = coeffs; //point to (shallow copy), since coeffs is large
    //real coefficient size, not fft size
  
    DFN(mPoly) = N(mPoly);
    DFSIZ(mPoly) = SIZ(mPoly);
    DEFDAT(mPoly) = DAT(mPoly);    
  }


  /**---------------------------------------------------
   * MultiplyByCutoff: multiply f1 and this 
   * @f1: DDMP
   * @f12: <= f1*this
   * choice of methods determined by cutoffs: 
   *       FFT, TFT, Kronecker, classical
   *      
   * assume f1 and this have the same number of vars
   * Return value: the product of f1 and this
   **/
  // DDMP<Field> MultiplyByCutoff(DDMP<Field>& f1){
  // 	sfixn N, sz1, sz2;
  //   preFFTRep *rep1, *rep2, *rep12;
  //   sfixn *dgs1, *dgs2;
    
  //   N = N(mPoly); //num of vars
  //   rep1 = f1.FFTRep();
  //   rep2 = mPoly;
    
  //   dgs1 = f1.PartialDegrees();
  //   dgs2 = BUSZS(mPoly);

  // 	rep12 = (preFFTRep *) my_calloc(1, sizeof(preFFTRep));
  //   PBPAS::InitResPoly(rep12, N, dgs1, dgs2); 
	
  // 	Field & pf = f1.CoeffField();
  //   MONTP_OPT2_AS_GENE mp = pf.getMontPrime();
	
  // 	sz1= PBPAS::getDenseSiz(N, rep1, BUSZSI(rep1, N), DAT(rep1), CUMI(rep1, N));
  // 	sz2= PBPAS::getDenseSiz(N, rep2, BUSZSI(rep2, N), DAT(rep2), CUMI(rep2, N));
	
  // 	if (( sz1 >= DENSECUTOFF ) || ( sz2 >= DENSECUTOFF )) 
  // 	  {
  // 		if (PBPAS::forcutoff_Multi_TFT_FFT(N, BUSZS(rep1), BUSZS(rep2), PARTDEGCUTOFF) ==0 )
  // 		  {
  // 			//some partial degs are smaller than 8
  // 			KroFFTRep * kPtr = (KroFFTRep *)my_calloc(1, sizeof(KroFFTRep));
  // 			InitKroFFTRep(kPtr, BUSZS(rep12), N, 2, &mp);
			
  // 			fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
  // 			fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
			
  // 			if ( ((KE(kPtr))>(mp.Npow)) || (KN(kPtr)<KNFFTCUTOFF) ) 
  // 			  {
  // 				//FFT size larger than the prime %ld can handle!  
  // 				//using MultiDFFT to solve the problem.
  // 				std::cout << "use FFT" << std::endl;
  // 				//----------------------------------------------
  // 				PBPAS::fftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), &mp);
  // 				//std::cout << "after fftMultiD_test" << std::endl;
  // 				my_free(DATSI(kPtr, 1));
  // 				DATSI(kPtr, 1)=NULL;
				
  // 				//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );
				
  // 				fromtofftRepMultiD(N, CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
  // 				//void freeKroFFTRep(KroFFTRep * x);
  // 				freeKroFFTRep(kPtr);
  // 			  } else 
  // 			  {
  // 				std::cout << "use KNFFT" << std::endl;
  // 				sfixn kdg1= (BUSZSI(rep1, N)+1)*CUMI(kPtr, N)-1;
  // 				sfixn kdg2= (BUSZSI(rep2, N)+1)*CUMI(kPtr, N)-1;
				
  // 				KROOTS(kPtr)=(sfixn *)my_calloc(KN(kPtr), sizeof(sfixn));
				
  // 				//-----------------------------
  // 				PBPAS::EX_Mont_GetNthRoots_OPT2_AS_GENE(KE(kPtr), KN(kPtr), KROOTS(kPtr), &mp);
  // 				//std::cout << "after GetNthRoots" << std::endl;
  // 				//-----------------------------
  // 				PBPAS::EX_KN_Mont_FFTMul_OPT2_AS_GENE_1(KN(kPtr), KE(kPtr), 0, KROOTS(kPtr), kdg1, DATSI(kPtr, 0), kdg2, DATSI(kPtr, 1), &mp);

  // 				//std::cout << "after KN_Mont_FFTMul" << std::endl;
  // 				my_free(DATSI(kPtr, 1));
  // 				DATSI(kPtr, 1)=NULL;
  // 				//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );			
  // 				fromtofftRepMultiD(N, CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
  // 				//std::cout << "after fromtofftRepMultiD" << std::endl;
  // 				//void freeKroFFTRep(KroFFTRep * x);
  // 				freeKroFFTRep(kPtr);
  // 			  }
  // 		  } else 
  // 		  {
  // 			std::cout << "use TFT" << std::endl;
  // 			KroTFTRep* kPtr = (KroTFTRep *)my_calloc(1, sizeof(KroTFTRep));
  // 			PBPAS::InitKroTFTRep(kPtr, BUSZS(rep12), N, 2, &mp);
  
  // 			fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 0), CUM(rep1), BUSZS(rep1), DAT(rep1));
  // 			fromtofftRepMultiD(N,  CUM(kPtr), DATSI(kPtr, 1), CUM(rep2), BUSZS(rep2), DAT(rep2));
  
  // 			//---------------------------
  // 			PBPAS::tftMultiD_test_1(DATSI(kPtr, 0), DATSI(kPtr, 1), N, ES(kPtr), DIMS(kPtr), LS(kPtr), &mp);
  
  // 			//result is in DATSI(kPtr, 0)
  // 			my_free(DATSI(kPtr, 1));
  // 			DATSI(kPtr, 1)=NULL;
  // 			//delayed allocation
  // 			//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );
  
  
  // 			fromtofftRepMultiD(N,  CUM(rep12), DAT(rep12), CUM(kPtr), BUSZS(rep12), DATSI(kPtr, 0));
  
  // 			PBPAS::freeKroTFTRep(kPtr);
  // 		  }
  // 	  }else 
  // 	  {
  // 		std::cout << "use classical" << std::endl;
  // 		//DAT(rep12)=(sfixn * )my_calloc( SIZ(rep12),sizeof(sfixn) );
  // 		PBPAS::plainMultiDMul(N, CUM(rep12), DAT(rep12), CUM(rep1), BUSZS(rep1), CUM(rep2), BUSZS(rep2), DAT(rep1), DAT(rep2), &mp);
  // 	  }
  
  // 	DDMP<Field> f12(mField);    
  // 	f12.setFFTRep(rep12);
  
  // 	return (f12); //return by copy
  // 	//my_free(rep12);
  // }
  

  /**
   * Multiply 2 polynomials with coefficients in a prime field
   *
   * @param f1 Pointer to the coefficient vector of polynomial f1 
   * @param f12 Pointer to the place holder of the product of f1 and this
   * @param opt Choice of method, default is simple contraction to 2D 
   */
  void MultiplyBy(DDMP<Field> *f1, DDMP<Field> *f12, sfixn opt=-11)
  {//----------------------------------------------------------------
    //default to use MultiplyByTFT_RBBnoE_saveOnDim1 to handle both bivar and multivar
 
    sfixn N;
    preFFTRep *rep1, *rep2, *rep12;
    sfixn *dgs1, *dgs2;
    
    N = N(mPoly); //num of vars
    rep1 = f1->FFTRep();
    rep2 = mPoly;
    
    dgs1 = f1->PartialDegrees();
    dgs2 = BUSZS(mPoly);
    
    rep12 = f12->FFTRep();
    //rep12 = (preFFTRep *) my_calloc(1, sizeof(preFFTRep));    
    //PBPAS::InitResPoly(rep12, N, dgs1, dgs2); 
    
    Field* pf= f1->CoeffField();
    MONTP_OPT2_AS_GENE* mp = pf->getMontPrime();

    switch(opt){
  
    case -7: { //use multi-TFT
      PBPAS::MultiplyByTFT_3D2(N, rep12, rep1, rep2, mp);
      break;
    }
    case -6: { //use multi-TFT
      PBPAS::MultiplyByTFT_2D2(N, rep12, rep1, rep2, mp);
      break;
    }
    case -4: { //use multi-TFT
      PBPAS::MultiplyByTFT_RBB(N, rep12, rep1, rep2, mp);
      break;
    } 
    case -3: { //nd to 1d by KN then to 2d, much more work
      PBPAS::MultiplyByTFT_RBB_KN1D_to_2D(N, rep12, rep1, rep2, mp);
      break;
    } 
    case -2: { //extend var1 to 2 vars and contract the 2nd with
      //var2..varn to 1 var
      PBPAS::MultiplyByTFT_1Vto2V_multiV_2V(N, rep12, rep1, rep2, mp);
      break;
    }   
    case -11: { //
      //simple contraction nd to 2d, no KroRep, saved in 1st dim eval 
      PBPAS::MultiplyByTFT_RBBnoE_saveOnDim1(N, rep12, rep1, rep2, mp);
      break;
    }  
    case -1: { //
      //simple contraction nd to 2d, very good
      PBPAS::MultiplyByTFT_RBBnoE(N, rep12, rep1, rep2, mp);
      break;
    }  
    case 0: { //2 processor ~= KN 1D-TFT 1 processor
      //std::cout<<"DDMP, SIZ="<<SIZ(rep12)<<std::endl;
      //rep1 and rep2 are univariate
      //evaluate 1st var for the num of d2 coefficients
      //saved by not tft on ls1*ls2/d2 number of vectors of zeros
      PBPAS::MultiplyByTFT_1V2V(N, rep12, rep1, rep2, mp);
      break;
    }  
    case 1: { //use multi-TFT
      PBPAS::MultiplyByTFT(N, rep12, rep1, rep2, mp);
      break;
    } 
    case 1101: { 
      //use multi-TFT, 2D transposition, 
      //no use of KroRep, saved on 1st var eval
      PBPAS::MultiplyByTFT_2DTran_noKroRep(N, rep12, rep1, rep2, mp);
      break;
    }
    case 110: { //use multi-TFT, 2D transposition, 
      PBPAS::MultiplyByTFT_2DTran(N, rep12, rep1, rep2, mp);
      break;
    }	  
    case 11: { //just for N=2, by TFT, no InitKroTFTRep, only cilk_for
      PBPAS::bivarMultiplyBy2DTFT(N, rep12, rep1, rep2, mp);
      break;
    } 
    case 111: { //just for N=2, by TFT, no InitKroTFTRep, 
      //spawn inside for
      PBPAS::bivarMultiplyBy2DTFT_spawn(N, rep12, rep1, rep2, mp);
      break;
    } 
    case 112: { //just for N=2, by TFT, no InitKroTFTRep, spawn inside
      //cilk_for 
      PBPAS::bivarMultiplyBy2DTFT_for_spawn(N, rep12, rep1, rep2, mp);
      break;
    } 	  
    case 12: { //for bivar mul, use 2D-TFT but no transposition
      PBPAS::MultiplyBy2DTFT_no_transp(N, rep12, rep1, rep2, mp);
      break;
    } 	  
    case 2: {//use multi-FFT
      PBPAS::MultiplyByFFT(N, rep12, rep1, rep2, mp);
      break;
    }
    case 3: {
      //convert to univariate representation by Kronecker substitution
      //then multiply by two evaluation and one interpolation via FFT
      PBPAS::MultiplyByKroneckerFFT(N, rep12, rep1, rep2, mp);
      break;
    }
    case 31: {
      //convert to univariate representation by Kronecker substitution
      //then multiply by two evaluation and one interpolation via TFT
      PBPAS::MultiplyByKronecker1DTFT(N, rep12, rep1, rep2, mp);
      break;
    }  
    case 4: {
      PBPAS::plainMultiDMul(N, CUM(rep12), DAT(rep12), CUM(rep1), BUSZS(rep1), CUM(rep2), BUSZS(rep2), DAT(rep1), DAT(rep2), mp);
    }
    case 5: {
      SERBPAS::MultiplyByTFT(N, rep12, rep1, rep2, mp);
      break;
    }
    case 6: {//use multi-FFT
      SERBPAS::MultiplyByFFT(N, rep12, rep1, rep2, mp);
      break;
    }
    case 7: {
      //convert to univariate representation by Kronecker substitution
      //then multiply by two evaluation and one interpolation via FFT
      SERBPAS::MultiplyByKroneckerFFT(N, rep12, rep1, rep2, mp);
      break;
    }
    case 8: {
      SERBPAS::plainMultiDMul(N, CUM(rep12), DAT(rep12), CUM(rep1), BUSZS(rep1), CUM(rep2), BUSZS(rep2), DAT(rep1), DAT(rep2), mp);
    }
  
    }
    
    //DDMP<Field> f12(mField);    
    //f12->setFFTRep(rep12);
    //return (f12); //return by copy
  }
};

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


