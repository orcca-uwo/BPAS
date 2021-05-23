#include "modpn.h"
#ifndef FFTSPE1
#define FFTSPE1
#define MY_PRIME1 4179340454199820289
#define INV_PRIME1 4611686018427387903
#define C_SFT1 4179340454199820288
#define SEE1 29
#define RINV1 0
#define RSFT1 2
#define NPOW1 57
namespace PBPAS1{
void Shuffle2(int n, sfixn* A,sfixn* B);
void Shuffle(int n, sfixn* A,sfixn* B);
sfixn testDFT(int n,int index,sfixn* A,sfixn *W);
void DFT_eff_p1(int n, int r,sfixn *A,sfixn *W,sfixn *B);

void InvDFTKeepMont_eff_p1(int n, int r, sfixn *A,sfixn *W,  sfixn *B,sfixn invn);

void InvDFT_eff_p1(int n, int r, sfixn *A,sfixn *W, sfixn *B,sfixn invn);

}
#endif
