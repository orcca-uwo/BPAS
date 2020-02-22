#ifndef __AS_h
#define __AS_h 


#include "Types.h"
#include "CONSTANTS.h"


extern sfixn BASE;
extern sfixn BASE_1;
extern sfixn BASEHALF;



#ifdef LINUXINTEL32

/**
 * MulHiLoUnsigned: 
 * @h: pointer to a sfixn "a".
 * @l: pointer to a sfixn "b".
 *
 * Computing sfixn multiplication and return the high word and low word.
 *           I.e. *h = HighPart(ab), *l = LowPart(ab)
 * Return value: void.
 **/


// If Intel 32 Linux machine use this...

static inline void
MulHiLoUnsigned(sfixn *h, sfixn *l)
{
  __asm__ ("mull %3" : "=a" (*l), "=d" (* h) : "%0" (* h), "rm" (* l));
}



#else


// generic C code.


// the following signature has been modified by Xin 
//  On Nov.19.2008. to fix the g++ compilation (sfixn<->usfixn) problem.
static void
//MulHiLoUnsigned(usfixn *h, usfixn *l)
MulHiLoUnsigned(sfixn *h, sfixn *l)
{
  ulongfixnum prod;
  prod=(ulongfixnum)(usfixn)(*h) * (ulongfixnum)(usfixn)(*l);
  //  printf("a*b=%ld\n", prod);
  

/* hack to fix old GCC on Sparc 32-bit */
/* (big endian)  prod = [hi 32, lo 32] */
#if SOLARIS64
  *h=  *((usfixn *)&prod);
#else
  *h=  (sfixn)(((ulongfixnum)prod)>>BASE);
#endif
  *l=  (sfixn)(usfixn)prod;

  // printf("*h=%ld\n", *h);
  // printf("*l=%ld\n", *l);
}


#endif

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


