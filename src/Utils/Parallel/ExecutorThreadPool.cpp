
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#include "Utils/Parallel/FunctionExecutorThread.hpp"
#include <iostream>
#include "../../../include/Utils/Parallel/Synchronized.hpp"


ExecutorThreadPool& ExecutorThreadPool::getThreadPool() {
	static ExecutorThreadPool pool(ExecutorThreadPool::maxThreads);
	return pool;
}

int ExecutorThreadPool::maxThreads = std::thread::hardware_concurrency() - 1;

const ExecutorThreadPool::threadID ExecutorThreadPool::notAThread = std::thread::id();

ExecutorThreadPool::ExecutorThreadPool(int n) :
	nThreads(n),
	threadPool(),
	taskPool(),
	nRetiredThreads(0),
	nPriorityThreads(0),
	maxPriorityThreads(n),
	priorityThreads(),
	m_mutex(),
	m_cv()
{
	poolThreads = new FunctionExecutorThread[n];
	std::function<void(FunctionExecutorThread*)> cb = std::bind((&ExecutorThreadPool::putbackThread), this, std::placeholders::_1);
	// fprintf(stderr, "Pool construction cb: %p\n", cb.target<void(*)(FunctionExecutorThread*)>());
	// const char* name = cb.target_type().name();
	// fprintf(stderr, "cb target type name: %s\n", name );
	for (int i = 0; i < n; ++i) {
		poolThreads[i].setCallback(cb);
		poolThreads[i].start();
		threadPool.push_back(poolThreads + i);
	}

}


ExecutorThreadPool::~ExecutorThreadPool() {
	// fprintf(stderr, "Killing thread pool, waiting for threads\n");
	waitForAllThreads();

	// fprintf(stderr, "thread pool size at end: %lu\n", threadPool.size());
	// fprintf(stderr, "num priority threads   : %lu\n", priorityThreads.size());

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
		taskPool.push_back(f);
	}
	tryPullTask();
}

void ExecutorThreadPool::addPriorityTask(std::function<void()>& f) {
	bool tryPull = 0;
	synchronized_nonrecursive(m_mutex) {
		if (!threadPool.empty()) {
			taskPool.push_front(f);
			tryPull = 1;
		} else if (nPriorityThreads < maxPriorityThreads) {
			FunctionExecutorThread* worker = new FunctionExecutorThread();
			std::function<void(FunctionExecutorThread*)> cb = std::bind(&ExecutorThreadPool::putbackThread, this, std::placeholders::_1);
			worker->setCallback(cb);
			worker->start();
			++nPriorityThreads;
			priorityThreads.push_back(worker);

			worker->sendRequest(f);
		} else {
			taskPool.push_front(f);
		}
	}

	if (tryPull) {
		tryPullTask();
	}
}



void ExecutorThreadPool::putbackThread(FunctionExecutorThread* t) {
	bool putback = 0;
//	fprintf(stderr, "TRYING PUTBACK %p\n", t);
	synchronized_nonrecursive(m_mutex) {
		if (nRetiredThreads < nPriorityThreads) {
			//if we launched a priority thread, don't put this thread back in the pool because
			//we want the total number of threads in flight to be nThreads (or, at work, nThreads+1).
			//And, really, it doesn't matter if the threads in the pool were created for priority tasks
			//or just the standard thread created at the beginning. They're all FunctionExecutorThreads.
			//Notice we also don't notify the cv in this case because we're waiting for a useable thread
			//to come back.
			++nRetiredThreads;
//			fprintf(stderr, "RETIRED A THREAD\n");
		} else if (!taskPool.empty()) {
//			fprintf(stderr, "TAKING TASK INSTEAD %p\n", t);
			std::function<void()> task = std::move(taskPool.front());
			taskPool.pop_front();
			t->sendRequest(task);
		} else {
			threadPool.push_front(t);
//			fprintf(stderr, "PUTBACK THREAD %d\n", threadPool.size());
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
//			fprintf(stderr, "TASK TO THREAD %p\n", worker);
			std::function<void()> task = std::move(taskPool.front());
			taskPool.pop_front();
			worker->sendRequest(std::move(task));
			// fprintf(stderr, "REQUEST SENT\n" );
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

//	bool idle = 0;
//	bool empty = 0;
//	while(!empty || !idle) {
//		synchronized_nonrecursive(m_mutex) {
//			empty = taskPool.empty();
//			idle = (threadPool.size() == nThreads);
//			bool allbusy = threadPool.empty();
//			fprintf(stderr, "Waiting for all threads. Pool size: %d, nPriority: %d, task size: %d, empty: %d, idle: %d, allbusy: %d\n", threadPool.size(), nPriorityThreads, taskPool.size(), empty, idle, allbusy);
//
//
//			if ((empty && !idle) || allbusy) {
//				m_cv.wait(lk);
//			}
//		}
//
//		if (!empty) {
//			tryPullTask();
//		}
//	}
//
	synchronized_nonrecursive(m_mutex) {
		while(threadPool.size() != nThreads || taskPool.size() > 0) {
//			fprintf(stderr, "Waiting for all threads. Pool size: %d, nPriority: %d, task size: %d\n", threadPool.size(), nPriorityThreads, taskPool.size());
			m_cv.wait(lk);
		}
	}



}

/**
 * Obtain, or reserve a certain number of theads from the pool,
 * returning those reserved threads unique IDs as a vector.
 *
 * These obtained threads are reserved uniquely for the caller of this function
 * and no longer participate in the executor's normal pool of threads.
 *
 * Clients must always return the threads obtained in this way.
 * These threads are temporarily "owned" by the client until they are returned to the pool
 * with returnThreads().
 *
 * @return the number of threads actually obtained, which may be less than requested.
 *
 * @see ExectureThreadPool::returnThreads
 */
int ExecutorThreadPool::obtainThreads(int numThreads, std::vector<ExecutorThreadPool::threadID>& reservedThreads) {

/**
 * To remove the threads from the queue of poolThreads without them actually
 * being busy with a task has interesting effects. In this case we do *not* want
 * the threads to return to the pool after processing a task. Rather, they should
 * be manually returned when done with whatever the client wanted to do by
 * explicitly calling returnThreads().
 *
 * To establish this, the side effect is that the FunctionExecutorThreads
 * temporarily have no "owning pool". Otherwise, they would go back to the
 * owning pool after processing a single request.
 * Clients should always return threads obtained in this way.
 */
	reservedThreads.clear();
	int numObtained = 0;
	synchronized_nonrecursive(m_mutex) {
		if (!threadPool.empty()) {
			int nAvail = threadPool.size();
			numObtained = nAvail < numThreads ? nAvail : numThreads;
			for (int i = 0; i < numObtained; ++i) {
				FunctionExecutorThread* worker = threadPool.front();
				threadPool.pop_front();
				worker->setCallback(std::function<void(FunctionExecutorThread*)>()); //temporarily doesn't exist in a pool, until returned.
				reservedThreads.emplace_back(worker->get_id());
			}
		}

	}

	return numObtained;

}

void ExecutorThreadPool::returnThreads(std::vector<ExecutorThreadPool::threadID>& reservedThreads) {
	std::vector<FunctionExecutorThread*> threadsToReturn;
	threadsToReturn.reserve(reservedThreads.size());

	synchronized_nonrecursive(m_mutex) {
		while(!reservedThreads.empty()) {
			threadID id = reservedThreads.back();
			reservedThreads.pop_back();
			if (id == notAThread) {
				continue;
			}
			for (size_t i = 0; i < nThreads; ++i) {
				if (poolThreads[i].get_id() == id) {
					threadsToReturn.push_back(poolThreads + i);
				}
			}
			for (size_t i = 0; i < nPriorityThreads; ++i) {
				if (priorityThreads[i]->get_id() == id) {
					threadsToReturn.push_back(priorityThreads[i]);
				}
			}
		}
	}

	//_1 is a placeholder for first param of cb
	std::function<void(FunctionExecutorThread*)> cb = std::bind((&ExecutorThreadPool::putbackThread), this, std::placeholders::_1);

	//Make sure you put back the threads outside the synchronized block.
	for (FunctionExecutorThread* ft : threadsToReturn) {
		ft->waitForThread();
		ft->setCallback(cb);
		putbackThread(ft);
	}
}

void ExecutorThreadPool::executeTask(const threadID id, std::function<void()>& f) {
	if (id == notAThread) {
		f();
		return;
	}

	for (size_t i = 0; i < nThreads; ++i) {
		if (poolThreads[i].get_id() == id) {
			poolThreads[i].sendRequest(f);
			return;
		}
	}

	FunctionExecutorThread* targetThread = NULL;
	synchronized_nonrecursive(m_mutex) {
		for (size_t i = 0; i < nPriorityThreads; ++i) {
			if (priorityThreads[i]->get_id() == id) {
				targetThread = priorityThreads[i];
				break;
			}
		}
	}

	if (targetThread != NULL){
		targetThread->sendRequest(f);
	}
}

void ExecutorThreadPool::waitForThreads(const std::vector<ExecutorThreadPool::threadID>& reservedThreads) {
	for (threadID id : reservedThreads) {
		if (id == notAThread) {
			continue;
		}
		bool found = false;
		for (size_t i = 0; i < nThreads; ++i) {
			if (poolThreads[i].get_id() == id) {
				poolThreads[i].waitForThread();
				found = true;
				break;
			}
		}
		if (!found) {
			FunctionExecutorThread* targetThread = NULL;
			//since priority threads are dynamically created, we need to synchronize
			//before iterating over them. Iterators could be invalidated by a
			//push_back.
			synchronized_nonrecursive(m_mutex) {
				for (size_t i = 0; i < nPriorityThreads; ++i) {
					if (priorityThreads[i]->get_id() == id) {
						targetThread = priorityThreads[i];
						break;
					}
				}
			}
			if (targetThread != NULL) {
				targetThread->waitForThread();
			}
		}
	}
}
