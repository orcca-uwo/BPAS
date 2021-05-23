#ifndef _POWERSERIES_PARA_H_
#define _POWERSERIES_PARA_H_

#include "PowerSeries.h"
#include "Utils/Parallel/ExecutorThreadPool.hpp"

/**
 * Update f by computing its homogeneous parts of degree up to and including d.
 *
 * This update occurs concurrently using the specified number of threads.
 * If nthreads is <= 0 then a dynmaic number of threads is chosen
 * determined by the value of d.
 *
 * @param d: the degree of the homogeneous part of f to return
 * @param f: the power series to update or obtain a homogeneous part of
 * @param nthreads: the number of threads to use in the computation.
 *
 * @return: the homogeneous part of f of degree d
 */
void updateToDegParallel_PS(int d, PowerSeries_t* f, int nthreads = 1);

/**
 * Update f by computing its homogeneous parts of degree up to and including d.
 *
 * This update occurs concurrently using the specified vector of
 * threadIDs.
 *
 * @param d: the degree of the homogeneous part of f to return
 * @param f: the power series to update or obtain a homogeneous part of
 * @param tids: the list of ExecutorThreadPool threadIDs to use for concurrent execution.
 *
 * @return: the homogeneous part of f of degree d
 */
Poly_ptr homogPartParallel_PS(int d, PowerSeries_t* f, std::vector<ExecutorThreadPool::threadID> tids);


#endif
