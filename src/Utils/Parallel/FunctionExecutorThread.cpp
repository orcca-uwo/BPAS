
#include "Utils/Parallel/ExecutorThreadPool.hpp"
#include "Utils/Parallel/FunctionExecutorThread.hpp"


void FunctionExecutorThread::processRequest(const std::function<void()>& f) {
	try {
		f();
	} catch (std::exception& e) {
		std::cerr << "Caught exception: " << e.what() << " while processing request." << std::endl;
	}
}

#if ONEWAY_EVENT_THREAD_PARALLEL
void FunctionExecutorThread::eventLoop() {
	std::function<void()> reqObj;
	// fprintf(stderr, "starting event loop: %p\n", this);
	while (requestQueue.getNextObject(reqObj)) {
		// std::cerr << std::this_thread::get_id() << " is about to process request" << std::endl;

		// fprintf(stderr, "processing starting %p.\n",  this);
		processRequest(reqObj);
		// fprintf(stderr, "processing done %p.\n", this );

		std::function<void(FunctionExecutorThread*)> localCB;
		synchronized_nonrecursive(m_mutex) {
			localCB = callback;
		}
		// fprintf(stderr, "about to call callback %p\n", this );
		if (localCB) {
			// fprintf(stderr, "really about to call callback %p\n", this );
			localCB((FunctionExecutorThread*)this);
		}
		// fprintf(stderr, "called back %p\n", this );

		int notify = 0;
		synchronized_nonrecursive(m_mutex) {
			if (requestQueue.streamEmpty()) {
				isIdle = true;
				notify = 1;
			}
		}
		if (notify) {
			m_cv.notify_all();
		}
	}

}
#endif