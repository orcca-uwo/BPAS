
#ifndef _EXECUTOR_THREAD_POOL_HPP_
#define _EXECUTOR_THREAD_POOL_HPP_


class FunctionExecutorThread;

#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <thread>

class ExecutorThreadPool {

private:

	size_t nThreads;
	FunctionExecutorThread* poolThreads;

	std::deque<FunctionExecutorThread*> threadPool;
	std::deque<std::function<void()>> taskPool;

	size_t nRetiredThreads;
	size_t nPriorityThreads;
	size_t maxPriorityThreads;
	std::vector<FunctionExecutorThread*> priorityThreads;

	std::mutex m_mutex;
	std::condition_variable m_cv;

	void putbackThread(FunctionExecutorThread* t);

	void tryPullTask();

	friend FunctionExecutorThread;

	/**
	 * Create a pool of FunctionExecutorThreads of size n.
	 * n: size of thread pool.
	 */
	ExecutorThreadPool(int n = std::thread::hardware_concurrency());

	~ExecutorThreadPool();

public:

	/**
	 * Obtain a reference to the singleton ExecutorThreadPool object.
	 * @return a reference to the singleton ExecutorThreadPool.
	 */
	static ExecutorThreadPool& getThreadPool();
	
	static int maxThreads;

	/**
	 * Add a task to the executor. 
	 */
	void addTask(std::function<void()>& f);

	/**
	 * Add a task to th executor pool which should be executed before other tasks.
	 * If the thread pool is empty, a temporary new thread will be created to service this task.
	 * This new priority thread will replace a standard thread once the
	 * standard thread returns to the thread pool. 
	 *
	 */
	void addPriorityTask(std::function<void()>& f);


	/**
	 * A query to see if all threads are busy and the pool is empty.
	 *
	 * @return true iff the pool is empty.
	 */
	bool allThreadsBusy();

	/**
	 * A blocking function call to wait for all tasks and threads to finish.
	 */
	void waitForAllThreads();
};


#endif
