// -*- C++ -*- 
// modpn-export.h
#ifndef __modpn_export_h
#define __modpn_export_h

#include "modpn.h"

namespace MODPN {

  //exported for c++ use ----yu
  
  sfixn * 
  Mont_dft_OPT2_AS_GENE ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr,  sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr );
  
  //exported for c++ use ----yu
  
  sfixn * 
  Mont_dft_OPT2_AS_GENE_SPE ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr,  sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr );
  
  //exported for c++ use ----yu
  
  sfixn * 
  Mont_invdft_OPT2_AS_GENE_SPE_R ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, sfixn degRes, sfixn * ResPtr, MONTP_OPT2_AS_GENE * pPtr);
  
  //exported for c++ use ----yu
  
  sfixn * 
  Mont_invdft_OPT2_AS_GENE_R ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, sfixn degRes, sfixn * ResPtr, MONTP_OPT2_AS_GENE * pPtr);
  
  int checkDgsOfST(sfixn N, TriSet * ts);

}//end of MODPN

#endif
