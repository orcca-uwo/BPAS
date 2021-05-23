
#ifndef _ONEWAY_EVENT_THREAD_HPP_
#define _ONEWAY_EVENT_THREAD_HPP_

#if defined(SERIAL) && SERIAL
#define ONEWAY_EVENT_THREAD_PARALLEL 0
#else
#define ONEWAY_EVENT_THREAD_PARALLEL 1
#endif

#if ONEWAY_EVENT_THREAD_PARALLEL
#include "Synchronized.hpp"
#include "AsyncObjectStream.hpp"
#include <thread>
#include <condition_variable>
#endif

#include "EventThreadException.hpp"
#include <vector>
#include <map>
#include <iostream>


/**
 * The OneWayEventThread class represents a long-running, event-loop based
 * thread for processing requests. It supports
 * many requesters and exaclty one responder--the event thread.
 *
 * Once a OneWayEventThread is created, one should call the start() method
 * to create and start the processing thread. The processing thread
 * runs continually until this object is destroyed.
 *
 * This class is abstract. One should subclass it and specialize
 * the processRequest method to actually perform the work required
 * to turn the Request object into a Response object.
 *
 * The class is templated by the object that it should
 * receive as a request.
 *
 */
template<class Request>
class OneWayEventThread {

protected:

#if ONEWAY_EVENT_THREAD_PARALLEL
	AsyncObjectStream<Request> requestQueue;

	std::thread m_worker;

	virtual void eventLoop() {
		Request reqObj;
		while (requestQueue.getNextObject(reqObj)) {
			// std::cerr << std::this_thread::get_id() << " is about to process request" << std::endl;
			processRequest(reqObj);
		}

	}

#else
	static long uniqueIDCount;

	long myID;

#endif

protected:

	/**
	 * Clean up the inter-thread communication resources
	 * and any resources used by the thread.
	 *
	 * Dervived classes should call this super-class method
	 * after they have cleaned up resources specific to their
	 * dervied implementation.
	 */
	virtual void threadCleanup() {
#if ONEWAY_EVENT_THREAD_PARALLEL
		requestQueue.resultsFinished();
		if (m_worker.joinable()) {
			m_worker.join();
		}
#endif
	}

	virtual void processRequest(const Request& reqObj) = 0;

public:

#if ONEWAY_EVENT_THREAD_PARALLEL
	typedef std::thread::id id;
#else
	typedef long id;
#endif

	/**
	 * Create a new event thread. Spins up the
	 * worker thread and the collector thread.
	 */
	OneWayEventThread()
#if ONEWAY_EVENT_THREAD_PARALLEL
	:	requestQueue(),
		m_worker()
	{
	
	}
#else
	{
		myID = ++uniqueIDCount;
	}
#endif

	virtual ~OneWayEventThread() {
		threadCleanup();
	}

	void start() {
#if ONEWAY_EVENT_THREAD_PARALLEL
		m_worker = std::thread(&OneWayEventThread::eventLoop, this);
#endif
	}


	/* Event threads are not copyable */
	OneWayEventThread(const OneWayEventThread& other) = delete;
	OneWayEventThread& operator=(const OneWayEventThread&) = delete;

	/**
	 * Get this OneWayEventThread's unique ID.
	 * @return the unique id
	 */
	id get_id() {
#if ONEWAY_EVENT_THREAD_PARALLEL
		return m_worker.get_id();
#else
		return myID;
#endif
	}

	/**
	 * Implements an asynchronous request send.
	 *
	 * reqObj: the Request to send.
	 */
	virtual void sendRequest(const Request& reqObj) {
#if ONEWAY_EVENT_THREAD_PARALLEL
		Request sendCopy = reqObj;
		requestQueue.addResult(sendCopy);
#else
		processRequest(reqObj);
#endif
	}


};

#endif
