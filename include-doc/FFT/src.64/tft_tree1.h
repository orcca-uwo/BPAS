#include <math.h>
#include <algorithm>
#include <iostream>
#include "modpn.h"
#include "transpose.h"
#include "fft_iter1.h"
#ifndef TFTSPE1
#define TFTSPE1
#define MY_PRIME1 4179340454199820289
#define INV_PRIME1 4179340454199820287
#define C_SFT1 4179340454199820288
#define SEE1 29
#define RINV1 3787527286618587136
#define RSFT1 2
#define R2 2559286960657440491
#define Mont_two 3458764513820540920
#define TFT_grainsize 2048
#define NPOW1 57
using namespace PBPAS1;
namespace TFT_tree1{
void Shuffle2(int n, sfixn* A,sfixn* B);
void Shuffle(int n, sfixn* A,sfixn* B);
sfixn testDFT(int n,int index,sfixn* A,sfixn *W);
void Shuffle(int n, sfixn *A, sfixn *B);

static inline void  TFT_AddSubSpeSSEModInplace(sfixn* a0,sfixn* a1, sfixn* a2, sfixn* a3);

void inline TFT_4POINT_p1(sfixn *A,sfixn *W);

void inline TFT_8POINT_p1(sfixn *A,sfixn *W);

void inline TFT_16POINT_p1(sfixn *A,sfixn *W);

void inline TFT_iter32_p1(sfixn *A,sfixn *W);

  void TFT_twiddle(sfixn *A, sfixn *KRT,int s, int r, int n, int oldn);

  void Shuffle_tft(int s, sfixn *testvec, sfixn *SB1);

  void TFT_Basecase_p1(int n, int r,sfixn *A,sfixn *W,sfixn *B);

  void TFT_Core_p1(int oldn, int n, int l, int m, int basecase, int relax_size, sfixn *KRT, int e, sfixn *SB1,sfixn *invec, sfixn *invectmp );

  sfixn MontDivMod(sfixn num1,sfixn a);

  sfixn normaltomont(sfixn normalnum, sfixn montnum);
 
  void ITFT_DFT_p1(int L2,int step_L2,sfixn *KRT,int d_L2,sfixn *SB, sfixn *SB2);

  void ITFT_Core_p1(int oldn, int L, sfixn zeta,sfixn onedivzeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec, sfixn *SB2);

  void ITFT_Wrapper_p1(int oldn, int L, sfixn zeta,sfixn onedivzeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec, sfixn *SB,int basecase, sfixn *invectmp, int relax_size);

  void Mont_ITFT_Core_p1(int oldn, int L, int index_zeta, int z, int n, int f, sfixn *KRT,sfixn *invKRT, sfixn *invec,  sfixn mont_inv2, sfixn *SB, sfixn p);

}
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


