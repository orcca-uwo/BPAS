#include "modpn.h"
#ifndef FURER_FFTSPE2
#define FURER_FFTSPE2
#define FURER_MY_PRIME2 2485986994308513793u
#define FURER_INV_PRIME2 4611686018427387903u
#define FURER_C_SFT2 2485986994308513792
#define FURER_SEE2 69
#define FURER_RINV2 0u
#define FURER_RSFT2 2
#define FURER_NPOW2 55
namespace FURERPBPAS2{
void Shuffle2(int n, sfixn* A,sfixn* B);
void Shuffle(int n, sfixn* A,sfixn* B);
void RootsTableFurer(int n, int r,sfixn *T);
sfixn testDFT(int n,int index,sfixn* A,sfixn *W);
void DFT_eff_p2(int n, int r,sfixn *A,sfixn *W,sfixn *B);

void InvDFTKeepMont_eff_p2(int n, int r, sfixn *A,sfixn *W,  sfixn *B,sfixn invn);

void InvDFT_eff_p2(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn);

}
#endif
