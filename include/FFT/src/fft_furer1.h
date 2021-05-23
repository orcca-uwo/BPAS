#include "modpn.h"
#ifndef FURER_FFTSPE1
#define FURER_FFTSPE1
#define FURER_MY_PRIME1 180143985094819841u
#define FURER_INV_PRIME1 288230376151711743u
#define FURER_C_SFT1 180143985094819840
#define FURER_SEE1 5
#define FURER_RINV1 0u
#define FURER_RSFT1 6
#define FURER_NPOW1 55
namespace FURERPBPAS1{
void Shuffle2(int n, sfixn* A,sfixn* B);
void Shuffle(int n, sfixn* A,sfixn* B);
void RootsTableFurer(int n, int r,sfixn *T);
sfixn testDFT(int n,int index,sfixn* A,sfixn *W);
void DFT_eff_p1(int n, int r,sfixn *A,sfixn *W,sfixn *B);

void InvDFTKeepMont_eff_p1(int n, int r, sfixn *A,sfixn *W,  sfixn *B,sfixn invn);

void InvDFT_eff_p1(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn);

}
#endif
