
#ifndef _FUNCTION_EXECUTOR_THREAD_HPP_
#define _FUNCTION_EXECUTOR_THREAD_HPP_

#include "AsyncGenerator.hpp"
#include "OneWayEventThread.hpp"
#include <exception>

class ExecutorThreadPool;

class FunctionExecutorThread : public OneWayEventThread<std::function<void()>> {

private:

	std::function<void(FunctionExecutorThread*)> callback;

	bool isIdle;

#if ONEWAY_EVENT_THREAD_PARALLEL
	std::mutex m_mutex;
	std::condition_variable m_cv;
#endif

public:

	FunctionExecutorThread() :
		OneWayEventThread<std::function<void()>>(),
		callback(),
		isIdle(true)
	{

	}

	/**
	 * Create a FunctionExecutorThread with a callback which
	 * gets called after each request has finished being processed.
	 * Such a callback gets called with "this" as a parameter.
	 *
	 * @param taskDoneCallBack : the callback to call once for every request processed.
	 *
	 */
	FunctionExecutorThread(std::function<void(FunctionExecutorThread*)> taskDoneCallBack) :
		OneWayEventThread<std::function<void()>>(),
		callback(taskDoneCallBack),
		isIdle(true)
	{

	}

	virtual ~FunctionExecutorThread() {}

	/**
	 * Set the callback which gets called after each request has finished processeing.
	 * Such a callback gets called with "this" as a parameter.
	 *
	 * @param taskDoneCallBack : the callback to call once for every request processed.
	 *
	 */
	void setCallback(const std::function<void(FunctionExecutorThread*)>& taskDoneCallBack) {
#if ONEWAY_EVENT_THREAD_PARALLEL
		synchronized_nonrecursive(m_mutex) {
			callback = taskDoneCallBack;
		}
#else
		callback = taskDoneCallBack;
#endif
	}


	/**
	 * Send a function to the thread to be executed.
	 *
	 * @param f the function to execute on the thread
	 */
	void sendRequest(const std::function<void()>& f) {
#if ONEWAY_EVENT_THREAD_PARALLEL
		synchronized_nonrecursive(m_mutex) {
			isIdle = false;
			OneWayEventThread<std::function<void()>>::sendRequest(f);
		}
#else
		isIdle = false;
		OneWayEventThread<std::function<void()>>::sendRequest(f);
		isIdle = true;
#endif
	}

	/**
	 * Wait until the FunctionExecutorThread is idle.
	 * A blocking function call, that waits until the FunctionExecutorThread
	 * is idle before returning. Can be used a synchronization point.
	 */
	void waitForThread() {
#if ONEWAY_EVENT_THREAD_PARALLEL
		synchronized_nonrecursive(m_mutex) {
			while (!isIdle) {
				m_cv.wait(lk);
			}
		}
#endif
	}

protected:

	virtual void processRequest(const std::function<void()>& f);

#if ONEWAY_EVENT_THREAD_PARALLEL
	virtual void eventLoop() override;
#endif
};



#endif
