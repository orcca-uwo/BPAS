

#ifndef _MAPLE_INTERFACE_PIPE_HPP_
#define _MAPLE_INTERFACE_PIPE_HPP_

#include "MapleInterface.hpp"

#include "../DataStructures/Factors.hpp"
#include "../Utils/fdstream.hpp"

#include "../Utils/Parallel/Synchronized.hpp"
#include <mutex>

class MapleInterfaceStream;

/**
 * The MapleInterfacePipe is an instance of a pipe to 
 * long-running Maple kernel. It is designed to  
 * serialize requests to a maple kernel.
 *
 * NOTE: This MapleInterfacePipe should either be the
 * first thing allocated within a main method or 
 * allocated via a static instance. Due to the fork
 * required for the pipe it should be done as soon
 * as possible to minimize usage requirements of fork
 * as well as cleaning up the child process.
 *
 * MapleInterfaceStream does exactly this.
 * 
 */
class MapleInterfacePipe : MapleInterface {
	
	friend MapleInterfaceStream;

	int c2p[2];
	int p2c[2];

	boost::fdostream* outstr;
	boost::fdistream* instr;

	bool child;
	int childPID;

	std::recursive_mutex m_mutex;

	void readLoop();

protected:

	virtual void sendCommand(const std::vector<std::string>& request, std::vector<std::string>& response);

public:

	MapleInterfacePipe();

	~MapleInterfacePipe();

};




#endif