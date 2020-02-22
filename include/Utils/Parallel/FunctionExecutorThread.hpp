
#ifndef _FUNCTION_EXECUTOR_THREAD_HPP_
#define _FUNCTION_EXECUTOR_THREAD_HPP_

#include "AsyncGenerator.hpp"
#include "OneWayEventThread.hpp"
#include <exception>

class ExecutorThreadPool;

class FunctionExecutorThread : public OneWayEventThread<std::function<void()>> {

private:

	ExecutorThreadPool* owningPool;

public:
	
	FunctionExecutorThread() : 
		owningPool(NULL),
		OneWayEventThread<std::function<void()>>()  {}

	FunctionExecutorThread(ExecutorThreadPool* pool) : 
		owningPool(pool),
		OneWayEventThread<std::function<void()>>() {}

	void setPool(ExecutorThreadPool* pool) {
		owningPool = pool;
	}

	virtual ~FunctionExecutorThread() {}
	
protected: 

	virtual void processRequest(const std::function<void()>& f);
};



#endif
