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

#include "Utils/Unix_Timer.h"

inline void startTimer(unsigned long long *start){
        _startTimer(start);
// #if !defined(SERIAL) && defined(CILKVIEW_TIMING)
//         *start = __cilkview_getticks();
// #else
//         *start = (unsigned long long)clock();
// #endif
}

inline void stopTimer(unsigned long long *start, float *elapsed){
        _stopTimer(start, elapsed);
// #if !defined(SERIAL) && defined(CILKVIEW_TIMING)
// 	*elapsed = (__cilkview_getticks() - *start) / 1000.f;
// #else
// 	*elapsed = (float)((unsigned long long)clock() - *start) / CLOCKS_PER_SEC;
// #endif
}

inline void stopTimerAddElapsed(unsigned long long *start, float *elapsed){
        _stopTimerAddElapsed(start, elapsed);
// #if !defined(SERIAL) && defined(CILKVIEW_TIMING)
// 	*elapsed += (__cilkview_getticks() - *start) / 1000.f;
// #else
// 	*elapsed += (float)((unsigned long long)clock() - *start) / CLOCKS_PER_SEC;
// #endif
}


#endif
