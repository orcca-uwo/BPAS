
#ifndef _EXECUTOR_THREAD_POOL_HPP_
#define _EXECUTOR_THREAD_POOL_HPP_


#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "FunctionExecutorThread.hpp"


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
	ExecutorThreadPool(int n = std::thread::hardware_concurrency() - 1);

	~ExecutorThreadPool();

public:

	typedef FunctionExecutorThread::id threadID;

	static const threadID notAThread;

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
	 * @param numThreads : number of threads wishing to be obtained from the pool
	 * @param[out] reservedThreads : the IDs of the threads which were obtained
	 * @return the number of threads actually obtained, which may be less than requested.
	 *
	 * @see ExectureThreadPool::returnThreads
	 */
	int obtainThreads(int numThreads, std::vector<threadID>& reservedThreads);

	/**
	 * Obtain, or reserve a thead from the pool,
	 * returning the reserved thread's unique ID.
	 * If all threads are busy, then a special ID is returned
	 * (ExectutorThreadPool::notAThread) that will still work for
	 * the executeTask, returnThread methods, but calling executeTake
	 * will not actually execute the task in parallel.
	 *
	 * The above says that client codes may follow the pattern of:
	 * 1) obtainThread, 2) executeTask, 3) returnThread;
	 * without regard for if the task is actually executed in parallel
	 * or not.
	 *
	 * The obtained thread is reserved uniquely for the caller of this function
	 * and no longer participates in the executor's normal pool of threads.
	 *
	 * Clients must always return the thread obtained in this way.
	 * The thread is temporarily "owned" by the client until it is returned to the pool
	 * with returnThread().
	 *
	 * @param[out] reservedTID : the ID of the thread which was obtained, or notAThread.
	 *
	 * @see ExectureThreadPool::returnThread
	 */
	void obtainThread(threadID& reservedID) {
		std::vector<threadID> tids;
		int n = this->obtainThreads(1, tids);
		if (n > 0) {
			reservedID = tids[0];
		} else {
			reservedID = notAThread;
		}
	}

	/**
	 * Given a vector of thread IDs where were returned from obtainThreads()
	 * return those threads back to the pool, hence releasing ownership
	 * of them.
	 *
	 * @param reservedThreads: the IDs of the threads to return to the pool.
	 */
	void returnThreads(std::vector<threadID>& reservedThreads);

	/**
	 * Given a thread ID which was obtained from obtainThread(),
	 * return that thread back to the pool, hence releasing ownership
	 * of it.
	 *
	 * @param reservedID: the ID of the thread to return to the pool.
	 */
	void returnThread(const threadID reservedID) {
		if (reservedID == notAThread) {
			return;
		}
		std::vector<threadID> tids;
		tids.push_back(reservedID);
		this->returnThreads(tids);
	}

	/**
	 * Given a thread ID, one obtained from obtainThreads(),
	 * ask the the Executor to execute the task f using
	 * the thread with ID id.
	 *
	 * @param id : the ID of the thread to execute the task
	 * @param f : the task to execute
	 */
	void executeTask(const threadID id, std::function<void()>& f);

	/**
	 * Given a particular thread ID obtained from obtainThreads()
	 * wait (block the current thread) until that thread has finished
	 * processing its tasks and is idle.
	 * @param id the thread id to wait for
	 */
	void waitForThread(threadID id) {
		if (id == notAThread) {
			return;
		}
		std::vector<threadID> ids;
		ids.emplace_back(id);
		waitForThreads(ids);
	}

	/**
	 * Given a list of thread IDs obtained from obtainThreads()
	 * wait (block the current thread) until those threads have
	 * finished processing all their tasks and are idle.
	 * @param reservedTreads the thread IDs to wait for
	 */
	void waitForThreads(const std::vector<threadID>& reservedThreads);


};


#endif
