
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#include "Utils/Parallel/FunctionExecutorThread.hpp"


void FunctionExecutorThread::processRequest(const std::function<void()>& f) {
	try {
		f();
	} catch (std::exception& e) {
		std::cerr << "Caught exception: " << e.what() << " while processing request." << std::endl;
	}

	if (owningPool != NULL) {
		owningPool->putbackThread(this);
	}
}
