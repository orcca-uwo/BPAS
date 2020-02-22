
#include "MapleInterface/MaplePipePool.hpp"

MaplePipePool::MaplePipePool(unsigned int n) :
	nPipes(n), 
	pipePool(),
	m_mutex() 
{
	allPipes = new MapleInterfacePipe[n];
	for (unsigned int i = 0; i < n; ++i) {
		pipePool.push_back(allPipes + i);
	}

}

MaplePipePool::~MaplePipePool() {
	synchronized_nonrecursive(m_mutex) {
		while (pipePool.size() != nPipes) {
			m_cv.wait(lk);
		}
		delete[] allPipes;
		while (!pipePool.empty()) {
			pipePool.pop_front();
		}
	}
}

MapleInterfacePipe* MaplePipePool::getMaplePipe() {
	MapleInterfacePipe* ret = NULL;

	synchronized_nonrecursive(m_mutex) {
		while (pipePool.empty()) {
			m_cv.wait(lk);		
		}

		ret = pipePool.front();
		pipePool.pop_front();
	}
	return ret;
}

void MaplePipePool::putbackPipe(MapleInterfacePipe* t) {
	synchronized_nonrecursive(m_mutex) {
		pipePool.push_front(t);
	}
	m_cv.notify_one();
}
