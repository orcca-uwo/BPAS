// basic-routine.h

#ifndef _pbpas_basic_h_
#define _pbpas_basic_h_

#include "general_routine.h"

namespace PBPAS {
  /*
      C = A*B mod p
      p: prime number
      N: number of variables
      C, A, B: dense recursive representation of C, A and B
      dc, da, db: partial degree vectors for C, A and B
      C should be initialized to zero
      default wp (whichprime) =0 to call the fft for a general Frouier prime
     */
    void DMPMul(sfixn p,  int N,    
		sfixn *C, int *dc, 
		sfixn *A, int *da,
		sfixn *B, int *db,
		int wp);

/*
      C = A*B mod pPtr->P
      two variables
      C, A, B: dense recursive representation of C, A and B
      dc, da, db: partial degree vectors for C, A and B
      dc[0] = partial degree of variable 1
      pPtr: struct for Montgmery mod arithmetic
      C should be initialized to zero
     */
    void BivarMul(sfixn *C, int *dc, 
		  sfixn *A, int *da,
		  sfixn *B, int *db,
		  MONTP_OPT2_AS_GENE *pPtr,
		  int wp);

/* 
       res should be initialized to zero
       s1 = partial degree of v1 of A + 1
       s2 = partial degree of v2 of A + 1
    */
    void Evaluation2D(sfixn *res, int es1, int es2, 
		      int dims1, int dims2,
		      int ls1, int ls2, 
		      sfixn *A, int s1, int s2,
		      sfixn *RT1, sfixn *RT2,
		      MONTP_OPT2_AS_GENE *pPtr,
		      int H, int *RevBidMap, int wp);

    void Interpolation2D(sfixn *res, int es1, int es2,
			 int dims1, int dims2, 
			 int ls1, int ls2,
			 sfixn *invRT1, sfixn *invRT2,
			 MONTP_OPT2_AS_GENE *pPtr, 
			 int H, int *RevBidMap, int wp);
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


