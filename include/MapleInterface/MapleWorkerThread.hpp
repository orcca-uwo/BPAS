
#ifndef _MAPLE_WORKER_THREAD_HPP_
#define _MAPLE_WORKER_THREAD_HPP_

#include <gmpxx.h>

#include "MapleInterface.h"
#include "../Utils/Parallel/EventThread.hpp"

class MapleWorkerThread : public EventThread<std::vector<std::string>,std::vector<std::string>> {
	
public:

	MapleWorkerThread() : EventThread<std::vector<std::string>,std::vector<std::string>>() {
		std::vector<std::string> request;
		request.emplace_back("startKernel");
		sendRequestAndWait(request);
	}
	
	virtual ~MapleWorkerThread() {
		threadCleanup();
	}

protected:

	virtual bool processTask(const std::vector<std::string>& request, std::vector<std::string>& response);

	virtual void threadCleanup();

private:

	/**
	 * Factor a polynomial. 
	 * request should be {"factors", <the_object_to_factor>}.
	 * response will be returned as {<numeric_factor>, <fact_1>, <exp_1>, ..., <fact_n>, <exp_n>}.
	 * 
	 * returns true iff a response was generated. 
	 */
	bool factors_mapleWorker(const std::vector<std::string>& request, std::vector<std::string>& response);

	/**
     * Get the GCD of two polynomials.
     * request should be {"gcd", <poly1>, <poly2>}
     * response will be returned as {<the_gcd>}
     *
	 * returns true iff a response was generated. 
	 */
	bool gcd_mapleWorker(const std::vector<std::string>& request, std::vector<std::string>& response);

	/**
     * Validate a triangularize result.
     * request should be {"triangularizeValiadte", <isLazard>, <F>, <rc>, <vars>, <our_results>}
     * response will be returned as {<boolean>}
     *
	 * returns true iff a response was generated. 
	 */
	bool triangularizeValidate_mapleWorker(const std::vector<std::string>& request, std::vector<std::string>& resposnse);


};


#endif