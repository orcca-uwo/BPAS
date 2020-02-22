
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#include "Utils/Parallel/FunctionExecutorThread.hpp"
#include <iostream>
#include "../../../include/Utils/Parallel/Synchronized.hpp"


ExecutorThreadPool& ExecutorThreadPool::getThreadPool() {
	static ExecutorThreadPool pool(ExecutorThreadPool::maxThreads);
	return pool;
}

int ExecutorThreadPool::maxThreads = std::thread::hardware_concurrency();

ExecutorThreadPool::ExecutorThreadPool(int n) :
	nThreads(n), 
	threadPool(),
	taskPool(),
	nRetiredThreads(0),
	nPriorityThreads(0),
	maxPriorityThreads(n),
	priorityThreads(n),
	m_mutex(), 
	m_cv()
{
	poolThreads = new FunctionExecutorThread[n];
	for (int i = 0; i < n; ++i) {
		poolThreads[i].setPool((ExecutorThreadPool*)this);
		threadPool.push_back(poolThreads + i);
	}

}


ExecutorThreadPool::~ExecutorThreadPool() {
	waitForAllThreads();

	//don't delete threads in threadPool directly since it is a 
	//mixture of both poolThreads and priorityThreads.
	delete[] poolThreads;
	for (size_t i = 0; i < nPriorityThreads; ++i) {
		delete priorityThreads[i];
	}
	priorityThreads.clear();
}

void ExecutorThreadPool::addTask(std::function<void()>& f) {
	synchronized_nonrecursive(m_mutex) {
		taskPool.push_back(std::move(f));
	}
	tryPullTask();
}

void ExecutorThreadPool::addPriorityTask(std::function<void()>& f) {
	bool tryPull = 0;
	synchronized_nonrecursive(m_mutex) {
		if (!threadPool.empty()) {
			taskPool.push_front(std::move(f));
			tryPull = 1;
		} else if (nPriorityThreads < maxPriorityThreads) {
			FunctionExecutorThread* worker = new FunctionExecutorThread();
			worker->setPool((ExecutorThreadPool*)this);

			++nPriorityThreads;
			priorityThreads.push_back(worker);
			
			worker->sendRequest(std::move(f));
		} else {
			taskPool.push_front(std::move(f));
		}
	}

	if (tryPull) {
		tryPullTask();
	}	
}



void ExecutorThreadPool::putbackThread(FunctionExecutorThread* t) {
	bool putback = 0;
	synchronized_nonrecursive(m_mutex) {
		if (nRetiredThreads < nPriorityThreads) {
			//if we launched a priority thread, don't put this thread back in the pool because
			//we want the total number of threads in flight to be nThreads (or, at work, nThreads+1). 
			//And, really, it doesn't matter if the threads in the pool were created for priority tasks
			//or just the standard thread created at the beginning. They're all FunctionExecutorThreads. 
			//Notice we also don't notify the cv in this case because we're waiting for a useable thread
			//to come back.
			++nRetiredThreads;
		} else if (!taskPool.empty()) {
			std::function<void()> task = std::move(taskPool.front());
			taskPool.pop_front();
			t->sendRequest(task);
		} else {
			threadPool.push_front(t);
			putback = 1;
		}
	}
	if (putback) {
		//use the condition variable ot notify thread waiting on all tasks to complete.
		m_cv.notify_all();
	}
}


void ExecutorThreadPool::tryPullTask() {
	synchronized_nonrecursive(m_mutex) {
		if (!taskPool.empty() && !threadPool.empty()) {
			FunctionExecutorThread* worker = threadPool.front();
			threadPool.pop_front();
			std::function<void()> task = std::move(taskPool.front());
			taskPool.pop_front();
			worker->sendRequest(task);
		}
	}
}

bool ExecutorThreadPool::allThreadsBusy() {
	bool busy = 0;
	synchronized_nonrecursive(m_mutex) {
		busy = threadPool.empty();
	}
	return busy;
}

void ExecutorThreadPool::waitForAllThreads() {
	//the below computation is tricky.
	//Since we actually also launch priority threads on occassion, 
	//the threads in the pool are (potentially) a mixture of regular
	//and priority threads. However, due to the way putBackThreads works, 
	//the pool should only ever be at most nThreads. 
	
	bool idle = 0;
	bool empty = 0;
	while(!empty || !idle) {
		synchronized_nonrecursive(m_mutex) {
			empty = taskPool.empty();
			idle = (threadPool.size() == nThreads);
			if (empty && !idle) {
				m_cv.wait(lk);
			}
		}
	
		if (!empty) {
			tryPullTask();
		}
	}
}