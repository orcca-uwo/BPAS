

#ifndef _MAPLE_INTERFACE_STREAM_HPP_
#define _MAPLE_INTERFACE_STREAM_HPP_

#include "MapleInterface.hpp"
#include "MapleInterfacePipe.hpp"
#include "MaplePipePool.hpp"

extern float g_mapleWaitTime;

#define MAPLE_INTERFACE_POOL 1


/**
 * The MapleInterfaceStream is a singleton 
 * designed to serialize requests to the maple kernel.
 * It acts as the single interface to the OpenMaple
 * sub-system (since Maple is not thread safe).
 *
 * Users can get the singleton instance via the 
 * static instance method.
 */
class MapleInterfaceStream : public MapleInterface {
	
#if (defined(SERIAL) && SERIAL) || !MAPLE_INTERFACE_POOL
	MapleWorkerThread m_thread;
#else
	MaplePipePool m_pool;
#endif

	
	MapleInterfaceStream() :
#if (defined(SERIAL) && SERIAL) || !MAPLE_INTERFACE_POOL
		m_thread()
#else
		m_pool()
#endif
	{

	}

	~MapleInterfaceStream() {

	}

	char* gcd_via_string(const char* a, const char* b);

protected:

	virtual void sendCommand(const std::vector<std::string>& request, std::vector<std::string>& response);

public:
	
	static MapleInterfaceStream _theMISInstance;

	/**
	 * Get the MapleInterfaceStream singleton instance.
	 * returns the instance.
	 */
	static MapleInterfaceStream& instance() {
		return MapleInterfaceStream::_theMISInstance;
	}

	/**
	 * A static helper method for getting the gcd of two
	 * polynomial represented as strings.
	 * This is used by C functions to call an externally
	 * defined C++ function.
	 */
	static char* gcd_string(const char* a, const char* b) {
		MapleInterfaceStream& mis = instance();
		return mis.gcd_via_string(a, b);
	}


};




#endif
