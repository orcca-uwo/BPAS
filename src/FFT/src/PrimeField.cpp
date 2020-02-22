// -*- C++ -*- 
//PrimeField.cilk
//@author Yuzhen Xie

#include <cilk/cilk.h>
#include "../../../include/FFT/src/modpn.h"
#include "../../../include/FFT/src/PrimeField.h"

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef int sfixn;

/**
 * PrimeField: constructor 
 * Create a prime struct for modulo computation by Montgomery trick
 * 
 * @param p a prime number
 */
PrimeField::PrimeField(sfixn p){
  mPrime = p;
  mMontPrime = (MONTP_OPT2_AS_GENE*)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
  EX_MontP_Init_OPT2_AS_GENE(mMontPrime, mPrime);
}

/**---------------------------------------------------------------
 * destructor
 **/
PrimeField:: ~PrimeField(){
  if (mMontPrime != NULL) {
    my_free(mMontPrime); 
    mMontPrime = NULL;
  }
}




