

#include <gmpxx.h>
#include "PowerSeries/PowerSeries_Parallel.hpp"
#include "Utils/Parallel/ExecutorThreadPool.hpp"

#include <cmath>
#include <sstream>

#include "Utils/Unix_Timer.h"

#define WORK_PER_THREAD 5000
#define WORK_PER_THREAD_DIV 5000

//If number of homog parts to compute is <=
//nthreads * COUNT_PER_THREAD_THREAD, distribute
//homog parts to compute equally. Otherwise, use work estimate.
#define COUNT_PER_THREAD_THRESH 4

#if PS_PARALLEL_TIME
float g_PStimer[24];
#endif



typedef ExecutorThreadPool::threadID threadID;


double _estimateWork(int n, int nvar) {
	double dn = (double) n;
	switch (nvar) {
		case 1: {
			return 1;
		}
		case 2: {
			return dn + 1;
		}
		case 3: {
			return 0.5*dn*dn + 1.5*n + 1;
		}
		case 4: {
			return 0.1666*dn*dn*dn + dn*dn + 1.8333*dn + 1;
		}
		case 5: {
			return 0.041666*dn*dn*dn*dn + 0.41666*dn*dn*dn + 1.45833*dn*dn + 2.0833*dn + 1;
		}
		default : {
			double term = 1;
			double denom = 1;
			for (int k = 2; k < nvar; ++k) { denom *= k; }
			denom = 1 / denom;
			term = 1;
			for (int k = 1; k < nvar; ++k) {
				term *= (double) (n + k);
			}
			return term*denom;
		}
	}
}


double _estimateTotalWork(int n, int nvar) {
	//For any homogeneous part of degree i, the maximum number of terms is (n+nvar-1) choose (nvar-1)
	double dn = (double) n;
	switch (nvar) {
		case 1: {
			return dn + 1;
		}
		case 2: {
			return 0.5*(dn+1)*(dn+2);
		}
		case 3: {
			double dn1 = dn+1;
			return 0.166*(dn1*dn1*dn1) + 0.5*(dn1*dn1) + 0.333*dn + 0.333;
		}
		case 4: {
			double dn1 = dn+1;
			return 0.04166*(dn1*dn1*dn1*dn1) + 0.25*(dn1*dn1*dn1) + 0.4583*(dn1*dn1) + 0.25*dn + 0.25;
		}
		case 5: {
			double dn1 = dn+1;
			return 0.00833*(dn1*dn1*dn1*dn1*dn1) + 0.08333*(dn1*dn1*dn1*dn1) + 0.29166*(dn1*dn1*dn1) +
				   0.41666*(dn1*dn1) + 0.2*n + 0.2;
		}
		default : {
			double sum = 0;
			double term = 1;
			double denom = 1;
			for (int k = 2; k < nvar; ++k) { denom *= k; }
			denom = 1 / denom;

			for (int i = 0; i <= n; ++i) {
				term = 1;
				for (int k = 1; k < nvar; ++k) {
					term *= (double) (i + k);
				}
				sum += (term*denom);
			}
			return sum;
		}
	}
}

void _updateAncestors_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {

    if (f->genOrder == 1) {
        if (f->gen.unaryGen == homogPartVoid_negate_PS) {
            updateToDeg_PS(d, (PowerSeries_t*) f->genParam1);
        }
    }

    if (f->genOrder == 2) {
        if (f->gen.binaryGen == homogPartVoid_sum_PS ||
        	f->gen.binaryGen == homogPartVoid_sub_PS ||
            f->gen.binaryGen == homogPartVoid_prod_PS) {

			std::function<void()> task1 = std::bind(updateToDeg_PS, d, (PowerSeries_t*) f->genParam1);
			std::function<void()> task2 = std::bind(updateToDeg_PS, d, (PowerSeries_t*) f->genParam2);
			if (tids.size() > 0) {
				ExecutorThreadPool::getThreadPool().executeTask(tids[0], task1);
			} else {
				task1();
			}
			task2();
			ExecutorThreadPool::getThreadPool().waitForThreads(tids);
        }
    }
}



#if PS_PARALLEL_TIME
void _updateDegRange_PS(int begin, int end, PowerSeries_t* f, float* timer) {
#else
void _updateDegRange_PS(int begin, int end, PowerSeries_t* f) {
#endif
	if (f == NULL || end < begin || f->genOrder == -1) {
		return;
	}

#if PS_PARALLEL_TIME
    unsigned long long start;
    _startTimer(&start);
#endif

	switch (f->genOrder) {
		case 1 : {
			for (int i = begin; i <= end; ++i) {
                f->polys[i] = (*f->gen.unaryGen)(i, f->genParam1);
			}
			break;
		}
		case 2 : {
			for (int i = begin; i <= end; ++i) {
                f->polys[i] = (*f->gen.binaryGen)(i, f->genParam1, f->genParam2);
			}
			break;
		}
		case 3 : {
			for (int i = begin; i <= end; ++i) {
                f->polys[i] = (*f->gen.tertiaryGen)(i, f->genParam1, f->genParam2, f->genParam3);
			}
			break;
		}
	}

#if PS_PARALLEL_TIME
    _stopTimer(&start, timer);
#endif

}


#if 0
void _updateDegRange_quo_PS(int begin, int end, std::mutex& cv_m, std::condition_variable& cv, PowerSeries_t* quo, float* timer) {

	//assert(quo->deg < begin);
	//memory leaks otherwise

	unsigned long long start;
	_startTimer(&start);

	PowerSeries_t* f = (PowerSeries_t*) quo->genParam1;
	PowerSeries_t* h = (PowerSeries_t*) quo->genParam2;

	std::unique_lock<std::mutex> lk(cv_m, std::defer_lock);

	if (begin == 0) {
		Poly_ptr q0 = deepCopyPolynomial_AA(homogPart_PS(0,f));
        divideByRational_AA_inp(q0, h->polys[0]->elems->coef);
        if (begin == 0){
            quo->polys[0] = q0;
        }
        begin = 1;
        lk.lock();
	    quo->deg = 0;
	    lk.unlock();
    }

    if (begin == 1) {
        Poly_ptr h1 = homogPart_PS(1,h);
        Poly_ptr f1 = homogPart_PS(1,f);
        Poly_ptr multipoly = multiplyPolynomials_AA(quo->polys[0], h1, h1 == NULL ? 0 : h1->nvar);
        multipoly = subPolynomials_AA_inp( multipoly, f1,  f1 == NULL ? 0 : f1->nvar);
        negatePolynomial_AA(multipoly);
        divideByRational_AA_inp(multipoly, h->polys[0]->elems->coef);
        quo->polys[1] = multipoly;
        lk.lock();
	    quo->deg = 1;
	    lk.unlock();
	}

	int i;
	// std::stringstream ss;
	// ss << "thread: " << std::this_thread::get_id() << " doing [" << begin << ", " << end << "]" <<  std::endl;
	// std::cerr << ss.str();
	// ss.str("");
	for (int k = begin; k <= end; ++k) {
		// ss << "thread: " << std::this_thread::get_id() << " doing k: " << k << std::endl;
		// std::cerr << ss.str();
		// ss.str("");
		Poly_ptr s = deepCopyPolynomial_AA(homogPart_PS(k,f));

	    for (i = 0; i < k; ++i) {
			// ss << "thread: " << std::this_thread::get_id() << " locking. i = " << i << " q->deg: " << quo->deg << std::endl;
			// std::cerr << ss.str();
			// ss.str("");

        	lk.lock();
	        if (quo->deg < i) {
        		cv.wait(lk, [=]{return quo->deg >= i;});
	        }
    		lk.unlock();
			// ss << "thread: " << std::this_thread::get_id() << " unlocking " << std::endl;
			// std::cerr << ss.str();
			// ss.str("");

	        //If we pass here, then quo->polys[i] is valid and, transitively,
	        //so is quo->polys[j] for j <= i.

	        Poly_ptr homog = homogPart_PS(k-i,h);
	        Poly_ptr p = multiplyPolynomials_AA(homog, quo->polys[i], homog == NULL ? 0 : homog->nvar);
	        s = subPolynomials_AA_inp(s, p, s == NULL ? 0 : s->nvar);
	        freePolynomial_AA(p);
	    }

	    divideByRational_AA_inp(s, h->polys[0]->elems->coef);

	    quo->polys[k] = s;
	    lk.lock();
	    quo->deg = k;
	    lk.unlock();
	    cv.notify_all();
	}

	_stopTimer(&start, timer);

}
#endif


#if PS_PARALLEL_TIME
void _updateDegRange_quo1_PS(int begin, int end, int maxDeg, PowerSeries_t* quo, float* timer) {
#else
void _updateDegRange_quo1_PS(int begin, int end, int maxDeg, PowerSeries_t* quo) {
#endif
	//assert(quo->deg < begin && begin > 1);
	//memory leaks otherwise

#if PS_PARALLEL_TIME
	unsigned long long start;
	_startTimer(&start);
#endif

	PowerSeries_t* f = (PowerSeries_t*) quo->genParam1;
	PowerSeries_t* h = (PowerSeries_t*) quo->genParam2;

	int nvar = quo->nvar;
	int i;
	Poly_ptr p, s;
	for (int k = begin; k <= end; ++k) {
		s = deepCopyPolynomial_AA(homogPart_PS(k,f));
		for (i = 0; i < maxDeg; ++i) {
			p = multiplyPolynomials_AA(homogPart_PS(k-i,h), quo->polys[i], nvar);
			s = subPolynomials_AA_inp(s, p, nvar);
			freePolynomial_AA(p); //TODO use prealloc
		}
		quo->polys[k] = s;
	}


#if PS_PARALLEL_TIME
	_stopTimer(&start, timer);
#endif

}

void _updateDegRange_quo2_PS(int begin, int end, PowerSeries_t* quo) {
	//At this point, quo[k] = f[k] - sum_{j=0}^{begin-1} (q[j]h[k-j]), for k=begin..end

	PowerSeries_t* f = (PowerSeries_t*) quo->genParam1;
	PowerSeries_t* h = (PowerSeries_t*) quo->genParam2;

	int nvar = quo->nvar;
	int i;
	Poly_ptr s, p;
	mpq_t h0;
	mpq_init(h0);
	mpq_set(h0, h->polys[0]->elems->coef);
	for (int k = begin; k <= end; ++k) {
		s = quo->polys[k];
		for (i = begin; i < k; ++i) {
			p = multiplyPolynomials_AA(homogPart_PS(k-i, h), quo->polys[i], nvar);
			s = subPolynomials_AA_inp(s, p, nvar);
			freePolynomial_AA(p); //TODO use prealloc
		}
    	divideByRational_AA_inp(s, h0);
    	quo->polys[k] = s;
	}
	mpq_clear(h0);

	quo->deg = end;
}



void _updateToDegParaBaseCase_PS(int d, PowerSeries_t* f, double totWork, std::vector<threadID> tids) {


	_updateAncestors_PS(d, f, tids);

    if (d + 1 > f->alloc) {
        int newAlloc = (2*f->alloc < d+1) ? d + 1 : 2*f->alloc;
    	resizePowerSeries_PS(f, newAlloc);
    }



#if PS_PARALLEL_TIME
    float timers[12];
#endif

    int nthreads = 	tids.size();
    ExecutorThreadPool& pool = ExecutorThreadPool::getThreadPool();

    int begin = f->deg + 1; int end;
    int workerIdx = 0;
	std::function<void()> task;
    if ((d - f->deg) <= (COUNT_PER_THREAD_THRESH*(nthreads+1))) {
    	//nthreads/2 in numerator does rounding.
    	//if not exact division, last (main) thread gets fewer
    	int targN = (d - f->deg + ((nthreads+1)/2)) / (nthreads+1);
    	// fprintf(stderr, "targN: %d\n", targN);
    	int curN = 0;
    	for (int i = f->deg + 1; i <= d; ++i) {
    		++curN;
    		if (curN >= targN) {
    			end = i;

#if PS_PARALLEL_TIME
    			task = std::bind(_updateDegRange_PS, begin, end, f, timers+workerIdx);
#else
    			task = std::bind(_updateDegRange_PS, begin, end, f);
#endif
    			if (workerIdx < nthreads) {
	    			// fprintf(stderr, "worker %d, update count deg range [%d,%d]\n",workerIdx, begin, end );
    				pool.executeTask(tids[workerIdx], task);
    				++workerIdx;
    			} else {
	    			// fprintf(stderr, "worker main, update count deg range [%d,%d]\n", begin, end );
    				task();
    				++workerIdx;
    			}
    			begin = i+1;
    			curN = 0;
    		}
    	}
    } else {
    	int targWork = totWork / (nthreads + 1);
	    double curWork = 0;
	    for (int i = f->deg + 1; i <= d; ++i) {
	    	curWork += _estimateWork(i, f->nvar);
	    	if (curWork >= targWork) {
	    		end = i;
#if PS_PARALLEL_TIME
    			task = std::bind(_updateDegRange_PS, begin, end, f, timers+workerIdx);
#else
    			task = std::bind(_updateDegRange_PS, begin, end, f);
#endif
	    		if (workerIdx < nthreads) {
	    			// fprintf(stderr, "worker %d, update deg range [%d,%d]\n",workerIdx, begin, end );
	    			pool.executeTask(tids[workerIdx], task);
	    			++workerIdx;
	    		} else {
	    			// fprintf(stderr, "worker main, update deg range [%d,%d]\n", begin, end );
	    			task();
	    		}
	    		begin = i+1;
	    		curWork = 0;
	    	}
	    }
    }

    //if last deg range was not committed inside the loop
    //(i.e. totWork or (d - f->deg) is not divisible by nthreads + 1)
    if (end < d) {
    	end = d;
		// fprintf(stderr, "worker main, update deg range [%d,%d]\n", begin, end );
#if PS_PARALLEL_TIME
		_updateDegRange_PS(begin, end, f, timers + workerIdx);
#else
		_updateDegRange_PS(begin, end, f);
#endif
    	++workerIdx;
    }

    pool.waitForThreads(tids);
#if PS_PARALLEL_TIME
    for (int i = 0; i < workerIdx; ++i) {
    	fprintf(stderr, "time[%d] = %8.4f\n", i, timers[i]);
    }
#endif

    f->deg = d;
}



//Done the pre-processing and determined that f really can be
//udpated in parallel. This is then the DnC function.
void _updateToDegParallel_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {

    if (d + 1 > f->alloc) {
        int newAlloc = (2*f->alloc < d+1) ? d + 1 : 2*f->alloc;
    	resizePowerSeries_PS(f, newAlloc);
    }

    int nvar = f->nvar;
    int nthreads = tids.size() + 1;
    double totWork = _estimateTotalWork(d, nvar) - _estimateTotalWork(d, f->deg);

    double targWork = WORK_PER_THREAD * nthreads;
    double curWork = 0.0;
    int curDeg = f->deg;
    int k = f->deg + 1;
    for ( ; k <= d; ++k) {
    	curWork += _estimateWork(k, f->nvar);
    	//accumulate the minimum amount of work AND
    	//if num of homog parts is less than threshold, ensure each thread
    	//gets a whole number
    	if (curWork > targWork) {
    		if (COUNT_PER_THREAD_THRESH*nthreads < (k - f->deg) || ((k - f->deg) % nthreads) == 0) {
	    		// fprintf(stderr, "calling base case for [%d,%d]\n", f->deg+1, k);
	    		_updateToDegParaBaseCase_PS(k, f, curWork, tids);
	    		break;
    		}
    	}
    }

    if (k > d) {
		// fprintf(stderr, "calling base case for [%d,%d]\n", f->deg+1, d);
    	//then we broke the loop without actually computing anything
    	_updateToDegParaBaseCase_PS(d, f, curWork, tids);
    } else if (k < d) {
    	//otherwise, there's still more work to do. So, recurse.
	    _updateToDegParallel_PS(d, f, tids);
    }
}


void _updateToDegParaDivBaseCase_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {
	// fprintf(stderr, "division base case [%d,%d]\n", f->deg+1, d);

	_updateAncestors_PS(d, f, tids);
	updateToDeg_PS(d, f);
}

void _updateToDegParallelDivide_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {

    if (d + 1 > f->alloc) {
        int newAlloc = (2*f->alloc < d+1) ? d + 1 : 2*f->alloc;
    	resizePowerSeries_PS(f, newAlloc);
    }

	int curDeg = f->deg;
	int nvar = f->nvar;
	double work = _estimateTotalWork(d, nvar) - _estimateTotalWork(f->deg, nvar);

	if (work < WORK_PER_THREAD_DIV || d-curDeg < tids.size()) {
		_updateToDegParaDivBaseCase_PS(d, f, tids);
		return;
	}

	int k = d;
	int nthreads = tids.size();
	double targWork = WORK_PER_THREAD_DIV*nthreads;
	double curWork = 0.0;
	for ( ; k > curDeg; --k) {
		curWork += _estimateWork(k, nvar);
		// if (curWork >= targWork && (d-k) >= (nthreads+1)) {
		if (curWork >= targWork && k < d && (d-k) % (nthreads+1) == 0) {
			break;
		}
	}

	//if we couldn't find a suitable rectangle, then proceed serially
	if (k == curDeg) {
		_updateToDegParaDivBaseCase_PS(d, f, tids);
		return;
	}


	//otherwise, our "rectangle" is [k,d]
	//1: recruse to compute [curDeg+1,k-1]
	//2: compute the rectangle of partial sums
	//3: compute the right triangle of the rectangle

	// Recursive call
	// fprintf(stderr, "division recursing for d=%d\n", k);
	_updateToDegParallelDivide_PS(k, f, tids);

	//Comupute Rectangle
	//First make sure ancestors are updated
	_updateAncestors_PS(d, f, tids);


#if PS_PARALLEL_TIME
	float time[12];
#endif
	int begin = k+1;
	std::function<void()> task;

	int perThread = (d-k) / (nthreads+1);
	int end = begin + perThread - 1;
	for (int i = 0; i < nthreads; ++i) {
		fprintf(stderr, "worker[%d] computing rectangle for [%d,%d], maxDeg=%d\n", i, begin, end, k+1);
#if PS_PARALLEL_TIME
		task = std::bind(_updateDegRange_quo1_PS, begin, end, k+1, f, time + i);
#else
		task = std::bind(_updateDegRange_quo1_PS, begin, end, k+1, f);
#endif
		ExecutorThreadPool::getThreadPool().executeTask(tids[i], task);
		// task();
		begin = end + 1;
		end += perThread;
	}
	fprintf(stderr, "worker main computing rectangle for [%d,%d]\n", begin, end );
#if PS_PARALLEL_TIME
	_updateDegRange_quo1_PS(begin, end, k+1, f, time + nthreads);
#else
	_updateDegRange_quo1_PS(begin, end, k+1, f);
#endif
	ExecutorThreadPool::getThreadPool().waitForThreads(tids);

#if PS_PARALLEL_TIME
	for (int i = 0; i < nthreads+1; ++i) {
		fprintf(stderr, "time[%d]: %10f\n", i, time[i]);
	}
	fprintf(stderr, "\n");
#endif

	//Compute right triangle of bottom trapezoid
	//At this point, quo[i] = f[i] - sum_{j=0}^{k-1} (q[j]h[i-j]), for i=k..d
	_updateDegRange_quo2_PS(k+1, d, f);
}


void updateToDegParallel_PS(int d, PowerSeries_t* f, int nthreads) {
    if (d <= f->deg) {
        return;
    }

    if (f->genOrder == -1) {
        return;
    }

	int inPara = 0;
    if (f->genOrder == 1) {
        inPara = 1;
    }
    if (f->genOrder == 2) {
        if (f->gen.binaryGen == homogPartVoid_sum_PS ||
        	f->gen.binaryGen == homogPartVoid_sub_PS ||
            f->gen.binaryGen == homogPartVoid_prod_PS) {
            inPara = 1;
        }
    } else if (f->genOrder == 3) {
        if (f->gen.tertiaryGen == homogPartVoid_quo_PS) {
            inPara = 2;
        }
    }

    if (inPara == 0) {
        updateToDeg_PS(d, f);
        return;
    }


    std::vector<threadID> tids;
    if (nthreads <= 0) {
	    int nvar = f->nvar;
	    double totWork = _estimateTotalWork(d, nvar) - _estimateTotalWork(d, f->deg);
	    nthreads = totWork / WORK_PER_THREAD;
    }

    if (nthreads > 1) {
    	ExecutorThreadPool::getThreadPool().obtainThreads(nthreads-1, tids);
    }

    if (inPara == 1) {
    	_updateToDegParallel_PS(d, f, tids);
    } else if (inPara == 2) {
    	_updateToDegParallelDivide_PS(d, f, tids);
    }


    ExecutorThreadPool::getThreadPool().returnThreads(tids);
}


void updateToDegParallel_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {
    if (d <= f->deg) {
        return;
    }

    if (f->genOrder == -1) {
        return;
    }

	int inPara = 0;
    if (f->genOrder == 1) {
        inPara = 1;
    }
    if (f->genOrder == 2) {
        if (f->gen.binaryGen == homogPartVoid_sum_PS ||
        	f->gen.binaryGen == homogPartVoid_sub_PS ||
            f->gen.binaryGen == homogPartVoid_prod_PS) {
            inPara = 1;
        }
    } else if (f->genOrder == 3) {
        if (f->gen.tertiaryGen == homogPartVoid_quo_PS) {
            inPara = 2;
        }
    }

    if (inPara == 0) {
        updateToDeg_PS(d, f);
        return;
    }

    if (inPara == 1) {
    	_updateToDegParallel_PS(d, f, tids);
    } else if (inPara == 2) {
    	_updateToDegParallelDivide_PS(d, f, tids);
    }
}




#if PS_PARALLEL_TIME
void _homogParaMultLoop(PowerSeries_t* g, PowerSeries_t* h, int d, int begin, int end, Poly_ptr* res, float* timer) {
#else
void _homogParaMultLoop(PowerSeries_t* g, PowerSeries_t* h, int d, int begin, int end, Poly_ptr* res) {
#endif

	//s = \sum_{i=begin}^end g->polys[d-i] * h->polys[i], i =begin...end

#if PS_PARALLEL_TIME
	unsigned long long start;
	_startTimer(&start);
#endif

	int nvar = g->nvar;
	Poly_ptr s = NULL;
	Poly_ptr  m;
	Poly_ptr tmp1, tmp2;
	for (int i = begin; i <= end; ++i) {
        tmp1 = homogPart_PS(d-i, g);
        if (isZero_AA(tmp1)) {
            continue;
        }
        tmp2 = homogPart_PS(i, h);
        if (isZero_AA(tmp2)) {
            continue;
        }
		m = multiplyPolynomials_AA(tmp1, tmp2, nvar);
		s = addPolynomials_AA_inp(s, m, nvar);
		freePolynomial_AA(m);
	}
	*res = s;

#if PS_PARALLEL_TIME
	_stopTimerAddElapsed(&start, timer);
#endif

}


Poly_ptr _homogParaMultiply_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {

	PowerSeries_t* g = (PowerSeries_t*) f->genParam1;
	PowerSeries_t* h = (PowerSeries_t*) f->genParam2;

	//Update ancestors
	//Play a little trick here for Weierstrass with checking if constant terms are 0
    if (isZero_AA(h->polys[0])) {
    	homogPartParallel_PS(d-1, g, tids);
 	} else {
 		homogPartParallel_PS(d, g, tids);
 	}
	if (isZero_AA(g->polys[0])) {
		homogPartParallel_PS(d-1, h, tids);
	} else {
		homogPartParallel_PS(d, h, tids);
	}

#if PS_PARALLEL_TIME
	unsigned long long start;
	_startTimer(&start);
#endif
	int nthreads = tids.size();
	Poly_ptr res[nthreads + 1];

	int incr = d / (nthreads+1);
	int begin = 0;
	int end = begin + incr;
	std::function<void()> task;
	for (int i = 0; i < nthreads; ++i) {
		// fprintf(stderr, "homog mult worker[%d] doing [%d,%d]\n", i, begin, end );
#if PS_PARALLEL_TIME
		task = std::bind(_homogParaMultLoop, g, h, d, begin, end, res + i, g_PStimer + i);
#else
		task = std::bind(_homogParaMultLoop, g, h, d, begin, end, res + i);
#endif
		ExecutorThreadPool::getThreadPool().executeTask(tids[i], task);
		begin = end + 1;
		end += incr;
	}
	end = d;
	// fprintf(stderr, "homog mult worker main doing [%d,%d]\n", begin, end );
#if PS_PARALLEL_TIME
	_homogParaMultLoop(g, h, d, begin, end, res+nthreads, g_PStimer+nthreads);
#else
	_homogParaMultLoop(g, h, d, begin, end, res+nthreads);
#endif

	ExecutorThreadPool::getThreadPool().waitForThreads(tids);

	int nvar = f->nvar;
	for (int i = 0; i < nthreads; ++i) {
		res[nthreads] = addPolynomials_AA_inp(res[nthreads], res[i], nvar);
		freePolynomial_AA(res[i]);
	}
	f->polys[d] = res[nthreads];

	f->deg = d;

	return f->polys[d];
}

inline Poly_ptr _homogPartParaAddSub_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {

	PowerSeries_t* g = (PowerSeries_t*) f->genParam1;
	PowerSeries_t* h = (PowerSeries_t*) f->genParam2;
	homogPartParallel_PS(d, g, tids);
	homogPartParallel_PS(d, h, tids);
	return homogPart_PS(d, f); //Doing a single homog part for add/sub in parallel makes no sense.
}

Poly_ptr homogPartParallel_PS(int d, PowerSeries_t* f, std::vector<threadID> tids) {
	if (f == NULL) {
		return NULL;
	}

	if (d <= f->deg) {
		return f->polys[d];
	} else if (f->genOrder == -1) {
		return NULL;
	}

	//an arbitrary cutoff to go serial
	//d>=15 implies at least 8 coefs for 2 threads, 4 coefs for 4 threads, etc.
	if (tids.size() == 0 || d < 7) {
		return homogPart_PS(d, f);
	}

	//an arbitrary threshold to call updateToDegParallel
	if (d - f->deg > tids.size()) {
		updateToDegParallel_PS(d, f, tids);
		return f->polys[d];
	}

    if (d + 1 > f->alloc) {
        int newAlloc = (2*f->alloc < d+1) ? d + 1 : 2*f->alloc;
    	resizePowerSeries_PS(f, newAlloc);
    }

	int inPara = 0;
    if (f->genOrder == 2) {
        if (f->gen.binaryGen == homogPartVoid_sum_PS
        	|| f->gen.binaryGen == homogPartVoid_sub_PS) {
        	for (int i = f->deg + 1; i <= d; ++i) {
        		_homogPartParaAddSub_PS(i, f, tids);
        	}
        	return f->polys[d];
        }
        else if (f->gen.binaryGen == homogPartVoid_prod_PS) {
        	for (int i = f->deg + 1; i <= d; ++i) {
        		_homogParaMultiply_PS(i, f, tids);
        	}
        	return f->polys[d];
        }
    }


    //As a last effort, just do in serial
    return homogPart_PS(d, f);


}


