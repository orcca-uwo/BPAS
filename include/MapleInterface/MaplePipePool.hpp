
#ifndef _MAPLE_PIPE_POOL_HPP_
#define _MAPLE_PIPE_POOL_HPP_


#include <queue>
#include <functional>
#include <mutex>
#include <thread>
#include <condition_variable>
#include "MapleInterfacePipe.hpp"

class MaplePipePool {

private:

	unsigned int nPipes;
	MapleInterfacePipe* allPipes;

	std::deque<MapleInterfacePipe*> pipePool;

	std::mutex m_mutex;
	std::condition_variable m_cv;

	friend MapleInterfacePipe;

public:

	/**
	 * Create a pool of FunctionMaplePipes of size n.
	 * n: size of thread pool.
	 */
	MaplePipePool(unsigned int n = std::thread::hardware_concurrency());

	~MaplePipePool();

	/**
	 * Get a MapleInterfacePipe from the pool.
	 * If the pool is empty, this method will block
	 * until the pool is re-populated. 
	 */
	MapleInterfacePipe* getMaplePipe();

	void putbackPipe(MapleInterfacePipe* t);

};


#endif
