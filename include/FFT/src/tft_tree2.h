#include <math.h>
#include <algorithm>
#include <iostream>
#include "modpn.h"
#include "transpose.h"
#include "fft_iter2.h"
#ifndef TFTSPE2
#define TFTSPE2
#define MY_PRIME2 2485986994308513793
#define INV_PRIME2 2485986994308513791
#define C_SFT2 2485986994308513792
#define SEE2 69
#define RINV2 1340102364119433216
#define RSFT2 2
#define R2 2559286960657440491
#define Mont_two 3458764513820540920
#define TFT_grainsize 2048
#define NPOW2 55
using namespace PBPAS2;
namespace TFT_tree2{
void Shuffle2(int n, sfixn* A,sfixn* B);
void Shuffle(int n, sfixn* A,sfixn* B);
sfixn testDFT(int n,int index,sfixn* A,sfixn *W);
void Shuffle(int n, sfixn *A, sfixn *B);

static inline void  TFT_AddSubSpeSSEModInplace(sfixn* a0,sfixn* a1, sfixn* a2, sfixn* a3);

void inline TFT_4POINT_p2(sfixn *A,sfixn *W);

void inline TFT_8POINT_p2(sfixn *A,sfixn *W);

void inline TFT_16POINT_p2(sfixn *A,sfixn *W);

void inline TFT_iter32_p2(sfixn *A,sfixn *W);

  void TFT_twiddle(sfixn *A, sfixn *KRT,int s, int r, int n, int oldn);

  void Shuffle_tft(int s, sfixn *testvec, sfixn *SB1);

  void TFT_Basecase_p2(int n, int r,sfixn *A,sfixn *W,sfixn *B);

  void TFT_Core_p2(int oldn, int n, int l, int m, int basecase, int relax_size, sfixn *KRT, int e, sfixn *SB1,sfixn *invec, sfixn *invectmp );

  sfixn MontDivMod(sfixn num1,sfixn a);

  sfixn normaltomont(sfixn normalnum, sfixn montnum);
 
  void ITFT_DFT_p2(int L2,int step_L2,sfixn *KRT,int d_L2,sfixn *SB, sfixn *SB2);

  void ITFT_Core_p2(int oldn, int L, sfixn zeta,sfixn onedivzeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec, sfixn *SB2);

  void ITFT_Wrapper_p2(int oldn, int L, sfixn zeta,sfixn onedivzeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec, sfixn *SB,int basecase, sfixn *invectmp, int relax_size);

  void Mont_ITFT_Core_p2(int oldn, int L, int index_zeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec,  sfixn mont_inv2, sfixn *SB, sfixn p);

}
#endif
