
#ifndef _EVENT_THREAD_HPP_
#define _EVENT_THREAD_HPP_

#if defined(SERIAL) && SERIAL
#define EVENT_THREAD_PARALLEL 0
#else 
#define EVENT_THREAD_PARALLEL 1
#endif

#if EVENT_THREAD_PARALLEL
#include "Synchronized.hpp"
#include "AsyncObjectStream.hpp"
#include <thread>
#include <condition_variable>
#endif

#include "EventThreadException.hpp"

#include "../Unix_Timer.h"

#include <vector>
#include <map>
#include <iostream>


extern float g_MapleComputeTime;

/**
 * The EventThread class represents a long-running, event-loop based
 * thread for unified access to some sub-module. It supports
 * many requesters and exaclty one responder--the event thread.
 *
 * The design is such that it automatically sets up a 
 * processing thread on creation. 
 *
 * This class is abstract. One should subclass it and specialize 
 * the processTask method to actually perform the work required
 * to turn the Request object into a Response object. 
 * If no response is required, one simply returns false.
 * See the documentation of the processTask method.
 *
 * The class can easily be used with the intent that
 * the processing thread will be a single point of 
 * access to some sub-system, simply make a singleton
 * of the sub-class 
 *
 * The class is templated by the object that it should
 * receive as a request and return as a response. 
 * 
 */
template<class Request, class Response = bool> 
class EventThread {
	
protected:

#if EVENT_THREAD_PARALLEL
	AsyncObjectStream<Request> requestQueue;
	AsyncObjectStream<std::pair<Request, Response>> responseQueue;

	std::thread m_worker;
	std::thread m_collector;

	/** Synchronization primitives for the response map **/
	std::mutex m_mutex;
	std::condition_variable m_cv;
#endif

	std::multimap<Request,Response> m_map;

#if EVENT_THREAD_PARALLEL
	void eventLoop() {
		Request reqObj;
		while (requestQueue.getNextObject(reqObj)) {
			// std::cerr << std::this_thread::get_id() << " is about to process request" << std::endl;
			Response respObj;
			unsigned long long start;
			_startTimer(&start);
			bool respond = processTask(reqObj, respObj);
			_stopTimerAddElapsed(&start, &g_MapleComputeTime);
			if (respond) {
				responseQueue.addResult(std::pair<Request,Response>(std::move(reqObj),std::move(respObj)));
			}
		}

	}

	void collectLoop() {
		std::pair<Request,Response> respObj;
		while (responseQueue.getNextObject(respObj)) {
			synchronized_nonrecursive(m_mutex) {
				m_map.insert(std::move(respObj));
			}
			m_cv.notify_all();
		} 
	}
#endif

	/**
	 * Clean up the inter-thread communication resources
	 * and any resources used by the thread. 
	 * 
	 * Dervived classes should call this super-class method
	 * after they have cleaned up resources specific to their
	 * dervied implementation. 
	 */
	virtual void threadCleanup() {
#if EVENT_THREAD_PARALLEL
		requestQueue.resultsFinished();
		responseQueue.resultsFinished();
		if (m_worker.joinable()) {
			m_worker.join();
		}
		if (m_collector.joinable()) {
			m_collector.join();
		}
#endif
	}

	virtual bool processTask(const Request& reqObj, Response& respObj) = 0;

public:

	/**
	 * Create a new event thread. Spins up the
	 * worker thread and the collector thread.
	 */
	EventThread() 
#if EVENT_THREAD_PARALLEL
	:	requestQueue(),
		responseQueue(),
		m_worker(),
		m_collector()
	{
		m_worker = std::thread(&EventThread::eventLoop, this);
		m_collector = std::thread(&EventThread::collectLoop, this);
	}
#else
	{}
#endif

	virtual ~EventThread() {
		threadCleanup();
	}

	/* Event threads are not copyable */
	EventThread(const EventThread& other) = delete;
	EventThread& operator=(const EventThread&) = delete;

	/**
	 * Submit a request to the event thread and wait here
	 * for the response. 
	 * This is provided for convenience and works exactly
	 * as if calling sendRequest() followed by getResponse().
	 *
	 * retObj: the request to make.
	 *
	 * returns the response of the request reqObj.
	 */
	Response sendRequestAndWait(const Request& reqObj) {
		sendRequest(reqObj);
		return getResponse(reqObj);
	}

	/**
	 * Implements an asynchronous request send.
	 * If a response is expected some time in the future,
	 * then one can call getResponse() or getResponse(Request)
	 * to get its response some time in the future.
	 *
	 * reqObj: the Request to send.
	 */ 
	void sendRequest(const Request& reqObj) {
#if EVENT_THREAD_PARALLEL
		Request sendCopy = reqObj;
		requestQueue.addResult(sendCopy);
#else 
		Response respObj;
		bool respond = processTask(reqObj, respObj);
		if (respond) {
			std::pair<Request,Response> respPair(reqObj,std::move(respObj));
			m_map.insert(std::move(respPair));
		} 
#endif
	}

	/**
	 * Get the response that matches the Request reqObj.
	 * If the response is not yet ready, it will wait
	 * until it is ready.
	 *
	 * Returns the response.
	 */
	Response getResponse(const Request& reqObj) {
		Response resp;
#if EVENT_THREAD_PARALLEL
		bool found = 0;
		synchronized_nonrecursive(m_mutex) {
			while(!found) {
				auto search = m_map.find(reqObj);
				if (search != m_map.end()) {
					resp = std::move(search->second);
					m_map.erase(search);
					found = 1;
				} else {
					//lk defined in synchronized_nonrecursive macro
					m_cv.wait(lk);
				}
			}
		}
#else 
		auto search = m_map.find(reqObj);
		if (search != m_map.end()) {
			resp = std::move(search->second);
			m_map.erase(search);
		}
#endif
		return resp;
	}

};



#endif

