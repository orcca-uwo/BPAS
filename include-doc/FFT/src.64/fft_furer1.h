#include "modpn.h"
#ifndef FFTSPE1
#define FFTSPE1
#define MY_PRIME1 180143985094819841u
#define INV_PRIME1 180143985094819839u
#define C_SFT1 180143985094819840
#define SEE1 5
#define RINV1 112589990684262400u
#define RSFT1 6
#define NPOW1 55
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


