
#ifndef _ASYNC_OBJ_STREAM_H_
#define _ASYNC_OBJ_STREAM_H_

#if defined(SERIAL) && SERIAL
	#error "AsyncObjectStream cannot be built serially.\nAre you sure this header should be included?"
#endif

#include <thread>
#include <future>
#include <mutex>
#include <queue>
#include <iostream>
#include <condition_variable>

#include "Synchronized.hpp"

/** 
 * A class to implement a (possibly) multi-threaded producer-consumer stream.
 * This supports one producer and one consumer.
 * Moreover, to avoid unnecessary locking and synchronization, some methods
 * should only be called by the consumer and some methods should only be called 
 * by the producer. See individual method documentation.
 *
 */
template <class Object> 
class AsyncObjectStream {

	std::queue<Object> retObjs;
	bool finished;

	std::mutex m_mutex;
	std::condition_variable m_cv;



public:

	/** 
	 * Construct a new AsyncObjectStream
	 */
	AsyncObjectStream() : 
		finished(0),
		m_mutex(),
		m_cv()
	{}

	~AsyncObjectStream() {
		resultsFinished();
	}

	/** 
	 * This method is only to be called by the consumer.
	 * 
	 * Gets the next object on the stream. 
	 * This method will block to wait for a result if one is pending.
	 * 
	 * returns true if the next object was successfully retrieved.
	 * returns false if the stream has finished and no new objects are available.
	 */
	bool getNextObject(Object& res) {

		if (isFinished()) {
			return false;
		}

		bool ret;
		synchronized_nonrecursive(m_mutex) {
			while(retObjs.empty() && !finished) {
				m_cv.wait(lk);
			}
			if (retObjs.empty() && finished) {
				ret = false;
			} else {
				ret = true;
				res = retObjs.front();
				retObjs.pop();	
			}
		}
		return ret;
	}

	/**
	 * This method is only to be called by the producer.
	 *
	 * Adding a new result will move from the Object res
	 * allowing for efficient data transfer.
	 */
	void addResult(Object& res) {	
		if (finished) {
			throw std::runtime_error("AsyncObjectStream declared finished by tried to add new result.");
		}
		synchronized_nonrecursive(m_mutex) {
			retObjs.push(std::move(res));
		}
		m_cv.notify_all();
	}

	void addResult(Object&& res) {	
		if (finished) {
			throw std::runtime_error("AsyncObjectStream declared finished by tried to add new result.");
		}
		synchronized_nonrecursive(m_mutex) {
			retObjs.push(std::move(res));
		}
		m_cv.notify_all();
	}

	/**
	 * This method is only to be called by the producer.
	 *
	 * Declare that the producer has completed producing
	 * and the previous call to addResult was the last result.
	 * Adding a new result will move from the Object res
	 * allowing for efficient data transfer.
	 */
	void resultsFinished() {
		synchronized_nonrecursive(m_mutex) {
			finished = true;
		}
		m_cv.notify_all();
	}

	bool isFinished() {
		bool noMoreObjects;
		synchronized_nonrecursive(m_mutex) {
			noMoreObjects = finished && retObjs.empty();
		}
		if (noMoreObjects) {
			return true;
		}
		return false;
	}

};




#endif