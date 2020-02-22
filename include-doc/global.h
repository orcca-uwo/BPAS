#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <math.h>
#include <unistd.h>

#include <algorithm>
#include <iterator>
#include <cstring>
#include <cctype>
#include <cmath>
#include <map>
#include <list>
#include <ctime>


#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <gmpxx.h>
#include "polynomial.h"		// Polynomial template

// For execution timing
#if !defined(SERIAL) && defined(CILKVIEW_TIMING)
	#include <cilktools/cilkview.h>
#else
	#include <time.h>
#endif

#ifdef SERIAL
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#endif

/**
 * Reserved global variable
 * Interval's width (input) used in real root isolation of regular chains
 **/
extern mpq_class bpas_root_width;

template <class Ring>
class UnivariateTerm {
        public:
                Ring coef;
                int exp;

                /**
                 * Tern constructor
                 *
                 * @param
                 **/
                UnivariateTerm<Ring> () : coef(), exp(0) {}
                UnivariateTerm<Ring> (Ring c, int e) : coef(c), exp(e) {}
                /**
                 * Copy constructor
                 *
                 * @param b: A univariate term
                 **/
                UnivariateTerm<Ring> (const UnivariateTerm<Ring>& b) : coef(b.coef), exp(b.exp) {}

                /**
                 * Overload operator =
                 *
                 * @param b: A univariate term
                 **/
                inline UnivariateTerm<Ring>& operator= (UnivariateTerm<Ring> b) {
                        if (this != &b) {
                                coef = b.coef;
                                exp = b.exp;
                        }
                        return *this;
                }

                inline friend std::ostream& operator<< (std::ostream &out, UnivariateTerm<Ring>& b) {
                        out << b.exp << ": " <<  b.coef;
                        return out;
                }
};

template <class Ring>
class MultivariateTerm {
        public:
                Ring coef;      // Coefficient
                int v;          // Number of variables
                int* degs;      // Degrees, in the ascending order of
                                // the weight of variates
                                // Such that degs[0] < .. < degs[v-1]

                /**
                 * Constructor
                 *
                 * @param
                 **/
                MultivariateTerm<Ring> () : coef(), v(0) { }
                /**
                 * Copy constructor
                 *
                 * @param b: A multivariate term
                 **/
                MultivariateTerm (const MultivariateTerm<Ring>& b) : coef(b.coef), v(b.v) {
                        degs = new int[v];
                        std::copy(b.degs, b.degs+v, degs);
                }
                /**
                 * Destructor
                 *
                 * @param
                 **/
                ~MultivariateTerm<Ring> () { if (v) { delete [] degs; } }

                /**
                 * Overload operator =
                 *
                 * @param b: A multivariate term
                 **/
                inline MultivariateTerm<Ring>& operator= (MultivariateTerm<Ring> b) {
                        if (this != &b) {
                                coef = b.coef;
                                v = b.v;
                                degs = new int[v];
                                std::copy(b.degs, b.degs+v, degs);
                        }
                        return *this;
                }
};

inline void startTimer(unsigned long long *start){
#if !defined(SERIAL) && defined(CILKVIEW_TIMING)
        *start = __cilkview_getticks();
#else
        *start = (unsigned long long)clock();
#endif
}

inline void stopTimer(unsigned long long *start, float *elapsed){
#if !defined(SERIAL) && defined(CILKVIEW_TIMING)
	*elapsed = (__cilkview_getticks() - *start) / 1000.f;
#else
	*elapsed = (float)((unsigned long long)clock() - *start) / CLOCKS_PER_SEC;
#endif
}

inline void stopTimerAddElapsed(unsigned long long *start, float *elapsed){
#if !defined(SERIAL) && defined(CILKVIEW_TIMING)
	*elapsed += (__cilkview_getticks() - *start) / 1000.f;
#else
	*elapsed += (float)((unsigned long long)clock() - *start) / CLOCKS_PER_SEC;
#endif
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


