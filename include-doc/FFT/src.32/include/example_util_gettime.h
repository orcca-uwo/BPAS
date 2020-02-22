// -*- C++ -*-

#ifndef CILK_EXAMPLE_UTIL_GETTIME_H_INCLUDED
#define CILK_EXAMPLE_UTIL_GETTIME_H_INCLUDED

/* 
 * This simple utility function hides Linux/Windows differences so that we can provide
 * examples with source code common to the different platforms.  Feel free
 * to use this unsupported file at your own risk.
 */

#ifdef _WIN32
#include <Windows.h>
#else
#include <stdlib.h>
#include <sys/time.h>
#endif

/* example_get_time
   get the time in milliseconds.  This means different things in Windows vs.
   Unix.  In Windows, it's a call to GetTickCount() which is the uptime of
   the system.  In Unix, it is implemented with gettimeofday() and sets the
   counter to zero the first time example_get_time is called.

   returns: the number of milliseconds since the start time.
 */
extern "C++"
inline
int example_get_time () {
#ifdef _WIN32
    // Windows implementation.
    return (int) GetTickCount();
#else
    static timeval *start = NULL;
    struct timezone tzp = { 0, 0 };
    if (NULL == start) {
        // define the current time as 0.
        start = (timeval*) malloc(sizeof(timeval));
        gettimeofday(start, &tzp);
        return 0;
    } else {
        // subtract the start time from the current time.
        timeval end;
        long ms = 0;
        gettimeofday(&end, &tzp);
        ms = (end.tv_sec - start->tv_sec) * 1000;
        ms += ((int)end.tv_usec - (int)start->tv_usec) / 1000;
        return (int) ms;
    }
#endif
}

/* example_random
   hash the value, n, that is passed using a 32-bit LCG:
   (a * n + b) % c
   The values of a, b, and c, are the same as those used in glibc.  This
   function exists to provide random numbers that are consistent across
   compilers and between 32 and 64 bit platforms.  It will also generate
   consistent results in parallel, provided that each node in the call-
   graph provides a consistent n.

   returns: the hashed value of n in the low 32 bits.
 */
extern "C++"
inline
unsigned int example_random (unsigned int n) {
  // 32-bit LCG used by glibc.
  // We use this to keep random numbers consistent between 32 and 64 bit
  // platforms and between compilers.
  return (((unsigned int)1103515245 * n) + (unsigned int)12345)
    % (unsigned int)0xFFFFFFFF;
}

#endif // CILK_EXAMPLE_UTIL_GETTIME_H_INCLUDED
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


