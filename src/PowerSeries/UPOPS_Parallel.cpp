

#include <gmpxx.h>
#include <math.h>
#include "PowerSeries/PowerSeries_Parallel.hpp"
#include "PowerSeries/UPOPS_Parallel.hpp"
#include "PowerSeries/UPOPS_Weierstrass.h"
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#include "Utils/Parallel/AsyncGeneratorPool.hpp"
#include "Utils/Unix_Timer.h"

typedef ExecutorThreadPool::threadID threadID;


#if PS_PARALLEL_TIME

float g_phase1Time = 0;
float g_phase1Other = 0;
float g_timer[24];
float g_alphaTimer[24];
extern float g_PStimer[24];
float g_waitTime[24];

#endif

#if PS_PARALLEL_TIME
void _weierstrassLemmaLoopBody(PowerSeries_t* G, PowerSeries_t* H, int r, int begin, int end, Poly_ptr* res, float* timer) {
#else
void _weierstrassLemmaLoopBody(PowerSeries_t* G, PowerSeries_t* H, int r, int begin, int end, Poly_ptr* res) {
#endif

#if PS_PARALLEL_TIME
	unsigned long long start;
	_startTimer(&start);
#endif

	int nvar = H->nvar;
	Poly_ptr s = NULL;
	Poly_ptr  m;
	for (int i = begin; i <= end; ++i) {
		m = multiplyPolynomials_AA(homogPart_PS(r-i, G), homogPart_PS(i, H), nvar);
		s = addPolynomials_AA_inp(s, m, nvar);
		freePolynomial_AA(m);
	}

	*res = s;

#if PS_PARALLEL_TIME
	_stopTimerAddElapsed(&start, timer);
#endif
}



Poly_ptr weierstrassLemmaParallel_UPOPS(PowerSeries_t* F, PowerSeries_t* G, PowerSeries_t* H, int r, std::vector<threadID> tids) {
	if (F == NULL || G == NULL || H == NULL) {
		return NULL;
	}

	unsigned long long start;
	if (r < 8) {
	_startTimer(&start);
		Poly_ptr ret = lemmaForWeierstrass_UPOPS(F, G, H, r);
// _stopTimerAddElapsed(&start, &g_phase1Other);
		return ret;
	}

	if (r == 200) {
		// const char* syms[] = {"x", "y"};
		// for (int i = 0; i < 200; ++i) {
		// 	if (isZero_AA(G->polys[i]) || isZero_AA(H->polys[r-i])) {
		// 		fprintf(stderr, "G[%d] = 0\n", i);
		// 		continue;
		// 	}
		// 	size_t max = mpz_sizeinbase(mpq_numref(G->polys[i]->elems->coef),2);
		// 	max += mpz_sizeinbase(mpq_denref(G->polys[i]->elems->coef),2);
		// 	max += mpz_sizeinbase(mpq_numref(G->polys[r-i]->elems->coef),2);
		// 	max += mpz_sizeinbase(mpq_denref(G->polys[r-i]->elems->coef),2);
		// 	fprintf(stderr, "G[%d]+H[%d] = %lu\n", i,r-i, max);
		// }
	}

#if PS_PARALLEL_TIME
	_startTimer(&start);
#endif

	int nthreads = tids.size();
	Poly_ptr res[nthreads + 1];

	int incr = r / (nthreads+1);
	int begin = 1;
	int end = begin + incr;
	std::function<void()> task;
	for (int i = 0; i < nthreads; ++i) {
		// fprintf(stderr, "worker[%d] doing [%d,%d]\n", i, begin, end );
#if PS_PARALLEL_TIME
		task = std::bind(_weierstrassLemmaLoopBody, G, H, r, begin, end, res + i, g_timer + i);
#else
		task = std::bind(_weierstrassLemmaLoopBody, G, H, r, begin, end, res + i);
#endif

		ExecutorThreadPool::getThreadPool().executeTask(tids[i], task);
		begin = end + 1;
		end += incr;
	}
	end = r-1;
// _stopTimerAddElapsed(&start, &g_phase1Other);

	// fprintf(stderr, "worker main doing [%d,%d]\n", begin, end );
#if PS_PARALLEL_TIME
	_weierstrassLemmaLoopBody(G, H, r, begin, end, res+nthreads, g_timer+nthreads);
#else
	_weierstrassLemmaLoopBody(G, H, r, begin, end, res+nthreads);
#endif

	ExecutorThreadPool::getThreadPool().waitForThreads(tids);

	int nvar = F->nvar;
	for (int i = 0; i < nthreads; ++i) {
		res[nthreads] = addPolynomials_AA_inp(res[nthreads], res[i], nvar);
		freePolynomial_AA(res[i]);
	}

	if(r > F->deg) {
		if (!isZero_AA(res[nthreads])) {
			negatePolynomial_AA(res[nthreads]);
		}
	} else {
		res[nthreads] = subPolynomials_AA_inp(res[nthreads], F->polys[r], nvar);
		negatePolynomial_AA(res[nthreads]);
	}

	divideByRational_AA_inp(res[nthreads], H->polys[0]->elems[0].coef);

	return res[nthreads];
}



int _computeAlphaWork(int d, int m) {
	//count work by number of power series multiplications
	int totalWork = 0;
	for (int i = m; i >= 0; --i) {
		if ((m-i) >= d) {
			totalWork += d;
		} else {
			//alpha[m] is no work, just a copy;
			//alpha[m-j] does j multiplications
			totalWork += (m-i);
		}
	}
	return totalWork;
}


/**
 * Return true iff distrib was sucessful and [begin[i],end[i]]
 * is the index range for worker[i] to update
 */
int _tryDistribAlphaWork(int totalWork, int monic, int d, int m, int nthreads, std::vector<int>& begin, std::vector<int>& end) {

	begin.resize(nthreads, -1);
	end.resize(nthreads, -1);
	int targWork = (totalWork + ((nthreads)/2)) / (nthreads);
	targWork = targWork <= 0 ? 1 : targWork;
	// fprintf(stderr, "In try distrib; nthreads: %d, totalWork: %d, targWork: %d\n", nthreads, totalWork, targWork );
	int curWork = 0;

	int workerIdx = 0;

	int beginIdx = m;
	int endIdx = m;

	for (int i = m; i >= 0; --i) {
		if ((m-i) > d) {
			curWork += d;
		} else {
			if (monic) {
				curWork += (m-i-1) < 0 ? 0 : m-i-1;
			} else {
				curWork += (m-i) < 0 ? 0 : m-i;
			}
		}

		if (curWork >= targWork) {
			if (workerIdx >= nthreads) {
				// fprintf(stderr, "bad distrib because workerIDx = %d, nthreads = %d\n", workerIdx, nthreads);
				return 0;
			}

			endIdx = i;

			//reverse because we are iterating m --> 0
			end[workerIdx] = beginIdx;
			begin[workerIdx] = i;
			// fprintf(stderr, "worker[%d] doing [%d,%d]\n", workerIdx, begin[workerIdx], end[workerIdx]);

			curWork = 0;
			beginIdx = i-1;
			++workerIdx;
		}
	}

	if (endIdx > 0 || workerIdx != nthreads) {
		// end[workerIdx] = beginIdx;
		// begin[workerIdx] = 0;
		// fprintf(stderr, "worker[%d] doing [%d,%d]\n", workerIdx, begin[workerIdx], end[workerIdx]);
		// return 1;

		return 0;
	} else {
		return 1;
	}

}

/**
 * Return true iff distrib was sucessful and [begin[i],end[i]]
 * represents the range for worker[i*threadsPer]
 */
int _computeAlphaWork(int monic, int d, int m, int nthreads, std::vector<int>& begin, std::vector<int>& end, int& threadsPer) {
	//count work by number of power series multiplications
	int totalWork = 0;
	for (int i = m; i >= 0; --i) {
		if ((m-i) > d) {
			totalWork += d;
			// fprintf(stderr, "work[%d] (m-i=%d) = %d\n", i, (m-i), d);
		} else {
			//alpha[m] is no work, just a copy;
			//alpha[m-j] does j multiplications
			if (monic) {
				totalWork += (m-i-1) < 0 ? 0 : m-i-1;
				// fprintf(stderr, "work[%d] monic (m-i=%d) = %d\n", i, (m-i), m-i-1);
			} else {
				totalWork += (m-i) < 0 ? 0 : m-i;
				// fprintf(stderr, "work[%d] (m-i=%d) = %d\n", i, (m-i), m-i);
			}
		}
	}

	threadsPer = 1;
	int success = _tryDistribAlphaWork(totalWork, monic, d, m, nthreads, begin, end);
	if (!success && nthreads % 2 == 0) {
		success = _tryDistribAlphaWork(totalWork, monic, d, m, nthreads/2, begin, end);
		if (success) {
			threadsPer = 2;
		} else {
			//could try again.. but unlikely when mixed with the hensel pipeline
			//unless using many many threads.
		}
	}

	return success;
}


#if PS_PARALLEL_TIME
void _updateAlpha_UPOPS(int d, Upops_t* alpha, int begin, int end, std::vector<threadID> tids, float* timer) {
	unsigned long long start;
	_startTimer(&start);
#else
void _updateAlpha_UPOPS(int d, Upops_t* alpha, int begin, int end, std::vector<threadID> tids) {
#endif

	for (int i = begin; i <= end; ++i) {
		if (tids.size() > 0) {
			homogPartParallel_PS(d, alpha->data[i], tids);
		} else {
			updateToDeg_PS(d, alpha->data[i]);
		}
	}

#if PS_PARALLEL_TIME
	_stopTimerAddElapsed(&start, timer);
#endif

}

void weierstrassUpdateParallel_UPOPS(Upops_t* p, std::vector<threadID> tids) {
	if (p == NULL || p->data[0]->paramType3 != UPOPS) {
		return;
	}


	Upops_t* alpha = (Upops_t*) p->data[0]->genParam3;


	int curdP = -1;
	int curdA = -1;
	for (int i = 0; i < p->deg; ++i) {
		curdP = curdP == -1 ? p->data[i]->deg : MIN(curdP, p->data[i]->deg);
	}
	for (int i = 0; i < alpha->deg; ++i) {
		curdA = curdA == -1 ? alpha->data[i]->deg : MIN(curdA, alpha->data[i]->deg);
	}

	//we look to increment each by one herein.
	curdP += 1;
	curdA += 1;

	int d = p->deg;
	int m = alpha->deg;
	Poly_ptr polyP;


	PowerSeries_t** F = p->weierstrassFData;

	for (int lP = 0; lP < d; ++lP) {
		// fprintf(stderr, "curdP %d, alloc: %d\n\n", curdP, p->data[lP]->alloc);
		if (curdP + 1 > p->data[lP]->alloc) {
			int newAlloc = (2*(p->data[lP])->alloc < curdP+1) ? curdP + 1 : 2*p->data[lP]->alloc;
			p->data[lP]->polys = (Poly_ptr*) realloc(p->data[lP]->polys, sizeof(Poly_ptr)*newAlloc);
			p->data[lP]->alloc = newAlloc;
		}
	}
	//in serial, alpha's allocation automatically managed by PS.
	//in parallel, since multi-reader one-writer, make sure allocaiton
	//does not change during updates.
	for (int lA = 0; lA <= m; ++lA) {
		if (curdA + 1 > alpha->data[lA]->alloc) {
			int newAlloc = (2*(alpha->data[lA])->alloc < curdA+1) ? curdA + 1 : 2*alpha->data[lA]->alloc;
			alpha->data[lA]->polys = (Poly_ptr*) realloc(alpha->data[lA]->polys, sizeof(Poly_ptr)*newAlloc);
			alpha->data[lA]->alloc = newAlloc;
		}
	}


	for (int lP = 0; lP <= d - 1; ++lP) {
		if (F[lP]->deg == -1){
			polyP = NULL;
		}else{

#if PS_PARALLEL_TIME
			unsigned long long start;
			_startTimer(&start);
#endif
			// fprintf(stderr, "updating F[%d] to %d with tids.size=%d\n", lP, curdP, tids.size());
			homogPartParallel_PS(curdP, F[lP], tids);
			// updateToDeg_PS(curdP,F[lP]); //update F based on new terms of p->data[lP-1]
			// fprintf(stderr, "\n", lP, curdP);
#if PS_PARALLEL_TIME
			_stopTimerAddElapsed(&start, &g_phase1Other);
			_startTimer(&start);
#endif

			polyP = weierstrassLemmaParallel_UPOPS(F[lP], p->data[lP], alpha->data[0], curdP, tids);
			// polyP = lemmaForWeierstrass_UPOPS(F[lP], p->data[lP], alpha->data[0], curdP);
#if PS_PARALLEL_TIME
			_stopTimerAddElapsed(&start, &g_phase1Time);
#endif
		}

		p->data[lP]->polys[curdP] = polyP;
		//manually set new degree since this updated PS is needed for next iteration of the loop.
		p->data[lP]->deg = curdP;
	}

	int nthreads = tids.size();
	std::vector<int> begin, end;
	int threadsPer;
	int goodDistrib;
	if (nthreads > 0) {
		goodDistrib = _computeAlphaWork(isOne_PS(alpha->data[alpha->deg]), d, m, nthreads+1, begin, end, threadsPer);
	} else {
		goodDistrib = 0;
		threadsPer = 1;
	}

	if (!goodDistrib) {
#if PS_PARALLEL_TIME
		unsigned long long start;
		_startTimer(&start);
#endif
		//update here because alpha[m] because it's just a copy
		// fprintf(stderr, "where alpha->deg=%d, not a good distrib! doing each alpha wiht %d\n", m, tids.size()+1 );
		updateToDeg_PS(curdA, alpha->data[m]);
		for (int i = m-1; i >= 0; --i) {
			homogPartParallel_PS(curdA, alpha->data[i],tids);
		}
#if PS_PARALLEL_TIME
		_stopTimerAddElapsed(&start, g_alphaTimer);
#endif
	} else {
		// int begin = 0;
		for (size_t i = 0; i < begin.size()-1; ++i) {
			if (begin[i] >= 0) {
				// fprintf(stderr, "where alpha->deg=%d, worker[%d] updating [%d,%d]\n", m, i*threadsPer, begin[i], end[i]);
				std::vector<threadID> childTids;
				for (int k = 1; k < threadsPer; ++k) {
					childTids.push_back(tids[i*threadsPer+k]);
				}
#if PS_PARALLEL_TIME
				std::function<void()> f = std::bind(_updateAlpha_UPOPS, curdA, alpha, begin[i], end[i], childTids, g_alphaTimer + i);
#else
				std::function<void()> f = std::bind(_updateAlpha_UPOPS, curdA, alpha, begin[i], end[i], childTids);
#endif
				ExecutorThreadPool::getThreadPool().executeTask(tids[i*threadsPer], f);
			}
		}
		if (begin[begin.size()-1] >= 0) {
			// fprintf(stderr, "where alpha->deg=%d, worker main updating [%d,%d]\n", m, begin[begin.size()-1], end[begin.size()-1]);
			std::vector<threadID> childTids;
			for (int k = 0; k < threadsPer-1; ++k) {
				childTids.push_back(tids[(begin.size()-1)*threadsPer+k]);
			}
#if PS_PARALLEL_TIME
			_updateAlpha_UPOPS(curdA, alpha, begin[begin.size()-1], end[begin.size()-1], childTids, g_alphaTimer + nthreads);
#else
			_updateAlpha_UPOPS(curdA, alpha, begin[begin.size()-1], end[begin.size()-1], childTids);
#endif
		}
		ExecutorThreadPool::getThreadPool().waitForThreads(tids);
	}
}


void weierstrassUpdateParallel_UPOPS(int d, Upops_t* p, int nthreads) {
	if (p == NULL) {
		return;
	} else if (p->data[0]->paramType3 != UPOPS) {
		updateToDeg_UPOPS(d, p);
		return;
	}

	std::vector<threadID> tids;
	ExecutorThreadPool::getThreadPool().obtainThreads(nthreads-1, tids);

	int curPrec = -1;
	for (int i = 0; i <= p->deg; ++i) {
		int prec = p->data[i] == NULL ? -1 : p->data[i]->deg;
		if (prec >= 0 && (curPrec < 0 || prec < curPrec)) {
			curPrec = prec;
		}
	}

	for (int i = curPrec+1; i <= d; ++i) {
		weierstrassUpdateParallel_UPOPS(p, tids);
	}


#if PS_PARALLEL_TIME
	for (int i = 0; i < nthreads; ++i) {
		fprintf(stderr, "g_PStimer[%d]: %10f\n",i,  g_PStimer[i]);
	}
	for (int i = 0; i < nthreads; ++i) {
		fprintf(stderr, "g_lemmaTimer[%d]: %10f\n",i,  g_timer[i]);
	}
	for (int i = 0; i < nthreads; ++i) {
		fprintf(stderr, "g_alphaTimer[%d]: %10f\n",i,  g_alphaTimer[i]);
	}

	fprintf(stderr, "g_phase1Other: %10f\n", g_phase1Other);
	fprintf(stderr, "g_phase1Time: %10f\n", g_phase1Time);
#endif

	ExecutorThreadPool::getThreadPool().returnThreads(tids);
}



void tshiftUpdateParallel_UPOPS(int d, Upops_t* f, std::vector<threadID> tids) {
	if (f == NULL || f->deg < 0) {
		return;
	}

	if (f->data[0]->genOrder != 3 || f->data[0]->paramType2 != MPQ_LIST) {
		//then this f is not the result of a taylor shift.
		return;
	}

	Upops_t* ancestor = (Upops_t*) f->data[0]->genParam3;
	int curD = ancestor->data[0]->deg;
	if (ancestor->weierstrassFData != NULL) {
		for (int i = curD+1; i <= d; ++i) {
			weierstrassUpdateParallel_UPOPS(ancestor, tids);
		}
	}

	//do the shift itself in serial
	updateToDeg_UPOPS(d, f);
}

//Middle of the pipeline
#if PS_PARALLEL_TIME
void _henselFactorUpdateGen_UPOPS(int d, Upops_t** facts, int nfacts, int nthreads, float* time, float* waitTime, AsyncGenerator<int>& consume, AsyncGenerator<int>& produce) {
#else
void _henselFactorUpdateGen_UPOPS(int d, Upops_t** facts, int nfacts, int nthreads, AsyncGenerator<int>& consume, AsyncGenerator<int>& produce) {
#endif
	if (facts == NULL || facts[0]->data[facts[0]->deg]->deg >= d) {
		return;
	}

#if PS_PARALLEL_TIME
	unsigned long long start;
	_startTimer(&start);
#endif

	std::vector<threadID> tids;
	ExecutorThreadPool::getThreadPool().obtainThreads(nthreads-1, tids);

	int curPrecision = facts[0]->data[facts[0]->deg]->deg;
	int updateDeg;

#if PS_PARALLEL_TIME
	unsigned long long start2;
	_startTimer(&start2);
	while (consume.getNextObject(updateDeg)) {
		_stopTimerAddElapsed(&start2,waitTime);
#else
	while (consume.getNextObject(updateDeg)) {
#endif

		//iterate from current precision to next available one at a time
		//This allows for consumers later in the pipeline to start work earlier.
		//And, consume may not give a nice sequence, but rather a single interger = d.
		for (int i = curPrecision + 1; i <= updateDeg; ++i) {
			// gmp_fprintf(stderr, "Fact->deg: %d, tc: %Qd, updating to %d\n", fact->deg, fact->data[0]->polys[0]->elems->coef, i);

			for (int k = 0; k < nfacts; ++k) {
				Upops_t* fact = facts[k];
				// gmp_fprintf(stderr, "Fact[%d]->deg: %d, tc: %Qd, updating to %d with tids=%d\n", k, fact->deg, fact->data[0]->polys[0]->elems->coef, i, tids.size());
				if (tids.size() == 0) {
					updateToDeg_UPOPS(i, facts[k]);
				} else {
					tshiftUpdateParallel_UPOPS(i, facts[k], tids);
				}
			}
			// updateToDeg_UPOPS(i, fact);
			// if (i % 8 == 0) {
				// Upops_t* fact = facts[0];
				// gmp_fprintf(stderr, "Fact[%d]->deg: %d, tc: %Qd, producing %d\n", 0, fact->deg, fact->data[0]->polys[0]->elems->coef, i);
				produce.generateObject(i);
			// }
		}
		curPrecision = updateDeg;

#if PS_PARALLEL_TIME
		_startTimer(&start2);
#endif
	}

	if (updateDeg < d) {
		for (int i = 0; i < nfacts; ++i) {
			updateToDeg_UPOPS(d, facts[i]);
		}
	}
	produce.generateObject(d);
	produce.setComplete();

	ExecutorThreadPool::getThreadPool().returnThreads(tids);

#if PS_PARALLEL_TIME
	_stopTimer(&start, time);
#endif

}


//Tail end of the pipeline
#if PS_PARALLEL_TIME
void _henselFactorUpdateGenTail_UPOPS(int d, Upops_t** facts, int nfacts, int nthreads, float* time, float* waitTime, AsyncGenerator<int>& consume) {
#else
void _henselFactorUpdateGenTail_UPOPS(int d, Upops_t** facts, int nfacts, int nthreads, AsyncGenerator<int>& consume) {
#endif
	if (facts == NULL || facts[0]->data[facts[0]->deg]->deg >= d) {
		return;
	}

#if PS_PARALLEL_TIME
	unsigned long long start;
	_startTimer(&start);
#endif

	std::vector<threadID> tids;
	ExecutorThreadPool::getThreadPool().obtainThreads(nthreads-1, tids);


	int curPrecision = facts[0]->data[facts[0]->deg]->deg;
	int updateDeg;

#if PS_PARALLEL_TIME
	unsigned long long start2;
	_startTimer(&start2);
	while (consume.getNextObject(updateDeg)) {
		_stopTimerAddElapsed(&start2,waitTime);
#else
	while (consume.getNextObject(updateDeg)) {
#endif
		//iterate from current precision to next available one at a time
		//This allows for consumers later in the pipeline to start work earlier.
		//And, consume may not give a nice sequence, but rather a single interger = d.
		for (int i = curPrecision + 1; i <= updateDeg; ++i) {

			for (int k = 0; k < nfacts; ++k) {
			Upops_t* fact = facts[k];
			// gmp_fprintf(stderr, "Tail Fact[%d]->deg: %d, tc: %Qd, updating to %d\n", k,fact->deg, fact->data[0]->polys[0]->elems->coef, i);
				if (tids.size() == 0) {
					updateToDeg_UPOPS(i, facts[k]);
				} else {
					tshiftUpdateParallel_UPOPS(i, facts[k], tids);
				}
			}
			// updateToDeg_UPOPS(i, fact);
		}
		curPrecision = updateDeg;

#if PS_PARALLEL_TIME
		_startTimer(&start2);
#endif
	}

	if (updateDeg < d) {
		for (int i = 0; i < nfacts; ++i) {
			updateToDeg_UPOPS(d, facts[i]);
		}
	}

	ExecutorThreadPool::getThreadPool().returnThreads(tids);

#if PS_PARALLEL_TIME
	_stopTimer(&start, time);
#endif
}



void updateHenselFactsParallel_UPOPS(int d, Upops_t** facts, int nfacts, std::vector<int> threadsPerFactor) {
	if (facts == NULL || nfacts == 0) {
		return;
	}

    //If no threads available, just update the last factor
    //to the requested precision. This will trigger a cascade
    //of updates.
	if (threadsPerFactor.size() == 0) {
    	updateToDeg_UPOPS(d, facts[nfacts-1]);
    	return;
    }

    if (threadsPerFactor.size() != nfacts) {
    	int nthreads = 0;
    	for (auto i : threadsPerFactor) {
    		nthreads += i;
    	}
    	updateHenselFactsParallel_UPOPS(d, facts, nfacts, nthreads);
    	return;
    }


    //Create a fake "head" of the pipeline so that we can re-use the
    //"middle" pipeline function for the first factor.
    //This fake head just produces one element: d.
    //We could make it produce currentPecision...d, if neeed.
    std::vector<int> v;
    v.push_back(d);
	AsyncGeneratorPool<int> head(v);

	int pipelineSize = 0;
	int factIdx = -1;
	for (int i = 0; i < nfacts; ++i) {
		if (factIdx < 0 && threadsPerFactor[i] > 0) {
			factIdx = i;
		}
		if (threadsPerFactor[i] > 0) {
			++pipelineSize;
		}
	}

	//current thread used as a pipeline stage, so request 1 less
	std::vector<threadID> pipelineTids;
	ExecutorThreadPool::getThreadPool().obtainThreads(pipelineSize-1, pipelineTids);
#if PS_PARALLEL_TIME
    float timers[pipelineTids.size()+1];
#endif

	//if we did not get multiple threads, just do it serially using the main thread.
	if (pipelineTids.size() == 0) {
#if PS_PARALLEL_TIME
		unsigned long long start;
		for (int i = 0; i < nfacts; ++i) {
			_startTimer(&start);
			updateToDeg_UPOPS(d, facts[i]);
			_stopTimerAddElapsed(&start, timers+i);
		}
		for (int i = 0; i < nfacts; ++i) {
			fprintf(stderr, "Time[%d]: %10f\n", i, timers[i]);
		}
#else
		for (int i = 0; i < nfacts; ++i) {
			updateToDeg_UPOPS(d, facts[i]);
		}
#endif
		return;
	}
//
	//Setup factors to store the resulting data without needing
	//many reallocs.
	for (int i = 0; i < nfacts-1; ++i) {
		resizeCoefficients_UPOPS(facts[i], d+2);
	}


// If we want to test just the parallelism within the first stage
#if 0
	// updateToDeg_UPOPS(d, facts[0]);
	AsyncGeneratorPool<int> producer;
	_henselFactorUpdateGen_UPOPS(d, facts, 1,
	//_henselFactorUpdateGen_UPOPS(d, facts+1, 1,
					threadsPerFactor[0], timers + 0, g_waitTime + 0, std::ref(head), producer);

	ExecutorThreadPool::getThreadPool().waitForThreads(pipelineTids);
	fprintf(stderr, "Phase 1 other: %10f\n", g_phase1Other );
	fprintf(stderr, "Phase 1 time: %10f\n", g_phase1Time );
	int weierstrassThreads = threadsPerFactor[0];
	for (int i = 0; i < weierstrassThreads; ++i) {
		fprintf(stderr, "PSTime[%d]: %10f\n", i, g_PStimer[i]);
	}
	for (int i = 0; i < weierstrassThreads; ++i) {
		fprintf(stderr, "Time[%d]: %10f\n", i, g_timer[i]);
	}
	for (int i = 0; i < 4; ++i) {
		fprintf(stderr, "AlphaTime[%d]: %10f\n", i, g_alphaTimer[i]);
	}
	fprintf(stderr, "Time[0]: %10f\n", timers[0] );

	ExecutorThreadPool::getThreadPool().returnThreads(pipelineTids);
	return;
#endif

	std::vector<AsyncGeneratorPool<int>*> gens;
	gens.reserve(nfacts);
	int prevFact = -1;
	AsyncGeneratorPool<int>* gen;
	for (size_t i = 0; i < pipelineTids.size(); ++i) {
		for (int j = prevFact + 1; j < nfacts; ++j) {
			if (threadsPerFactor[j] > 0) {
				factIdx = j;
				break;
			}
		}
		//fprintf(stderr, "stage %lu has prevFact %d, factIdx %d\n", i, prevFact, factIdx );
		//fprintf(stderr, "stage %lu has facts [%d, %d); nfacts = %d, has threads: %d\n", i, prevFact+1, (factIdx - prevFact) + prevFact + 1, factIdx - prevFact, threadsPerFactor[factIdx]);

		if (prevFact == -1) {
			gen = new AsyncGeneratorPool<int>(pipelineTids[i], _henselFactorUpdateGen_UPOPS,
					d, facts, factIdx - prevFact, threadsPerFactor[factIdx],
#if PS_PARALLEL_TIME
					timers + 0, g_waitTime + 0,
#endif
					std::ref(head));
		} else {
			gen = new AsyncGeneratorPool<int>(pipelineTids[i], _henselFactorUpdateGen_UPOPS,
					d, facts + prevFact + 1, factIdx - prevFact, threadsPerFactor[factIdx],
#if PS_PARALLEL_TIME
					timers + i, g_waitTime + i,
#endif
					std::ref(*gens[i-1]));
		}

		//Uncomment this to see how the stages work if run in serial
		// ExecutorThreadPool::getThreadPool().waitForThreads(pipelineTids);
		gens.push_back(gen);

		prevFact = factIdx;
	}
	int finalThreads = 0;
	for (int j = factIdx + 1; j < nfacts; ++j) {
		// fprintf(stderr, "threads per factor[%d]: %d\n", j, threadsPerFactor[j]);
		finalThreads += threadsPerFactor[j];
	}
	//For the tail, pass the very last factor since nthreads may be less than
	//the number of factors. This ensures the group of factors of index
	// tids.size() to nfacts-1 all get updated.
	// fprintf(stderr, "About to call update for facts[%d] of deg: %d\n",nfacts-1, facts[nfacts-1]->deg );


	//fprintf(stderr, "stage last has facts [%d, %d]; nfacts = %d, has threads: %d\n", prevFact+1, prevFact + 1 + (nfacts-1-prevFact), nfacts-1-prevFact, finalThreads);
	_henselFactorUpdateGenTail_UPOPS(d, facts + prevFact + 1, nfacts - 1 - prevFact, finalThreads,
#if PS_PARALLEL_TIME
		timers + pipelineTids.size(), g_waitTime + pipelineTids.size(),
#endif
		*gens[gens.size()-1]);


	//This should be unnecessary by serial call to above update tail, but we'll do it anyways
	ExecutorThreadPool::getThreadPool().waitForThreads(pipelineTids);

#if PS_PARALLEL_TIME
	for (size_t i = 0; i < pipelineTids.size()+1; ++i) {
		fprintf(stderr, "Time[%lu]: %10f\n", i, timers[i]);
	}
	for (size_t i = 0; i < pipelineTids.size()+1; ++i) {
		fprintf(stderr, "Wait time[%lu]: %10f\n", i, g_waitTime[i]);
	}
#endif

	ExecutorThreadPool::getThreadPool().returnThreads(pipelineTids);

}


//Estimate the cost in a Weierstrass preparation of produced p,alpha of degree d,m
int _weierstrassWorkInHensel(int d, int m) {
	int cost = d; //lemma cost
    int t = d < m ? d : m;

    //phase 1
    for (int i = 0; i <= t && i < d; ++i) {
    	cost += i;
    }
    for (int i = t+1; i <= d-1; ++i) {
    	cost += m;
    }

    //phase 2
    for (int i = 1; i <= t; ++i) {
    	cost += (i-1);
    }
    for (int i = t+1; i <= m; ++i) {
    	cost += d;
    }

    return cost;
}

int _henselWorkPerFactor(Upops_t** facts, int nfacts, std::vector<int>& perFactor) {
	perFactor.clear();

	int tdeg = 0;
	for (int i = 0; i < nfacts; ++i) {
		tdeg += facts[i]->deg;
	}

	int total = 0;
	for (int i = 0; i < nfacts-1; ++i) {
		tdeg -= facts[i]->deg;
		perFactor.push_back(_weierstrassWorkInHensel(facts[i]->deg, tdeg));
		total += perFactor[i];
	}

	//last factor is nothing
	perFactor.push_back(0);


	return total;
}

void updateHenselFactsParallel_UPOPS(int d, Upops_t** facts, int nfacts, int nthreads) {
	if (facts == NULL || nfacts == 0) {
		return;
	}

	// fprintf(stderr, "hensel nthreads in: %d\n", nthreads );

	std::vector<int> workPerFactor;
	int totalWork = _henselWorkPerFactor(facts, nfacts, workPerFactor);

	std::vector<int> threadsPerFactor(nfacts);
	float targetRatio = 1.0f / nthreads;
	float curRatio = 0.0f;
	int assignedThreads = 0;
	for (int i = 0; i < nfacts; ++i) {
	// for (int i = workPerFactor.size() - 1; i >= 0; --i) {
		curRatio += (float) workPerFactor[i] / totalWork;

		// fprintf(stderr, "workPerFactor[%d] : %d\n",  i, workPerFactor[i]);
		// fprintf(stderr, "curRatio[%d] : %10f\n",  i, curRatio);
		if (curRatio >= targetRatio) {
			threadsPerFactor[i] = (int) roundf(curRatio*nthreads);
			assignedThreads += threadsPerFactor[i];
			curRatio = 0.0f;
		}
	}
	for (int i = 0; i < nfacts && assignedThreads < nthreads; ++i) {
		threadsPerFactor[i] += 1;
		++assignedThreads;
	}
	if (threadsPerFactor[nfacts-1] == 0) {
		for (int i = nfacts-2; i >= 0; --i) {
			if (threadsPerFactor[i] > 0) {
				threadsPerFactor[i] -= 1;
				threadsPerFactor[nfacts-1] = 1;
				break;
			}
		}
	}


#if PS_PARALLEL_TIME
	for (int i = 0; i < nfacts; ++i) {
		fprintf(stderr, "fact[%d] threads = %d\n", i, threadsPerFactor[i] );
	}
#endif

	updateHenselFactsParallel_UPOPS(d, facts, nfacts, threadsPerFactor);

}
