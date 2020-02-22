#include "modpn.h"
#ifndef FFTSPE2
#define FFTSPE2
#define MY_PRIME2 2485986994308513793
#define INV_PRIME2 2485986994308513791
#define C_SFT2 2485986994308513792
#define SEE2 69
#define RINV2 1340102364119433216
#define RSFT2 2
#define NPOW2 55
namespace PBPAS2{
void Shuffle2(int n, sfixn* A,sfixn* B);
void Shuffle(int n, sfixn* A,sfixn* B);
sfixn testDFT(int n,int index,sfixn* A,sfixn *W);
void DFT_eff_p2(int n, int r,sfixn *A,sfixn *W,sfixn *B);

void InvDFTKeepMont_eff_p2(int n, int r, sfixn *A,sfixn *W,  sfixn *B,sfixn invn);

void InvDFT_eff_p2(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn);

}
#endif
