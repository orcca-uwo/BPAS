// -*- C++ -*- 
//PrimeField.h
#ifndef __PrimeField_h
#define __PrimeField_h

#include "modpn.h"

//sfixn is defined in modpn such that
//#ifdef LINUXINTEL64
//typedef int sfixn;

class PrimeField {
public:
  PrimeField(sfixn p);
  ~PrimeField();

  sfixn prime() { return mPrime; }
  MONTP_OPT2_AS_GENE* getMontPrime(){ return mMontPrime; }
  
private:
  //prime coefficient field
  sfixn mPrime;
  //prime struct for modulo computation by Montgomery trick
  MONTP_OPT2_AS_GENE* mMontPrime;

};

#endif
