


#include "MapleInterface/MapleInterfaceStream.hpp"

#include "IntegerPolynomial/mzpolynomial.hpp"


MapleInterfaceStream MapleInterfaceStream::_theMISInstance;

void MapleInterfaceStream::sendCommand(const std::vector<std::string>& request, std::vector<std::string>& response) {
#if (defined(SERIAL) && SERIAL) || !MAPLE_INTERFACE_POOL
	response = m_thread.sendRequestAndWait(request);
#else
	MapleInterfacePipe* pipe = m_pool.getMaplePipe();
	pipe->sendCommand(request, response);
	m_pool.putbackPipe(pipe);
#endif		
}



// #include "Utils/Unix_Timer.h"
// #include <atomic>

/**
 * Hacky way to give C code access to a C++ method.
 */
char* (*gcd_maple_string)(const char*, const char*) = &MapleInterfaceStream::gcd_string;

float g_mapleWaitTime = 0;

// std::atomic<long> mapleGCDTimeSec;
// std::atomic<long> mapleGCDTimeUSec;

// std::atomic<long> mapleFactorTimeSec;
// std::atomic<long> mapleFactorTimeUSec;

char* MapleInterfaceStream::gcd_via_string(const char* a, const char* b) {

	std::vector<std::string> request, response;
	request.emplace_back("gcd");
	request.emplace_back(a);
	request.emplace_back(b);

	// timer_id id = start_timer();
	// std::cerr << std::this_thread::get_id() << " is about to wait for C GCD" << std::endl;
	unsigned long long start;
	startTimer(&start);
	sendCommand(request, response);
	// response = m_thread.sendRequestAndWait(request);
	stopTimerAddElapsed(&start, &g_mapleWaitTime);
	// std::cerr << std::this_thread::get_id() << " finished waiting for C GCD" << std::endl;

    // timer_time elapsed1 = elapsed_time(&id);
	// mapleGCDTimeSec += elapsed1.tv_sec;
 //    mapleGCDTimeUSec += elapsed1.tv_usec;
 //    double time = (mapleGCDTimeSec + ((double)mapleGCDTimeUSec / 1000000));
 //    fprintf(stderr, "\nmapleGCDTime: %f\n", time);

	char* ret = (char*) malloc(sizeof(char)*response[0].length()+1);
	strcpy(ret, response[0].c_str());
	return ret;
}
