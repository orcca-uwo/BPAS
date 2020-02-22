
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

#include "Synchronized.hpp"

/** 
 * A class to implement a (possibly) multi-threaded producer-consumer stream.
 * This supports one producer and one consumer.
 * Moreover, to avoid unnecessary locking and synchronization, some methods
 * should only be called by the consumer and some threads should only be called 
 * by the producer. See individual method documentation.
 *
 */
template <class Object> 
class AsyncObjectStream {

	std::queue<Object> retObjs;
	bool finished;

	std::promise<Object> promise;
	bool promiseSet;
	
	/**
	 * In general, we must synchronize on the access to the future.
	 * Since the producer thread creates new futures and the consumer thread
	 * accesses that future it must be synchronized. 
	 */ 
	std::future<Object>  future;

	std::recursive_mutex m_mutex;

	void prepareNewFuture() {
		synchronized(m_mutex) {
			if (finished && retObjs.empty()) {
				//do nothing, leave future invalid
				//to mark stream is closed 
			} else {
				bool futureNotTaken = future.valid();
				if (!futureNotTaken) {
					promise = std::promise<Object>();
					future = promise.get_future();
					promiseSet = 0;
				}
			}
		}
	}

	/**
	 * This method is only to be called by the producer.
	 */
	void tryPrepareNextResult() {

		synchronized(m_mutex) {
			if(finished && retObjs.empty()) {
				//nothing to do...
			} else if(!promiseSet && !retObjs.empty()) {
				promise.set_value(std::move(retObjs.front()));
				retObjs.pop();
				promiseSet = 1;
			}
			//Do nothing? Either producer is busy and the queue is empty or we are waiting 
			//for consumer to take a new value.


			//Consumer will prepare new future and new result 
			//  else {
			// 	//promise is set.

			// 	bool futureNotTaken = future.valid();
			// 	//if futureNotTaken then we must just wait for 
			// 	//the consumer to take the value already put
			// 	//into the promise. Otherwise, both the promise
			// 	//was set and the consumer has retrieved the object
			// 	//therefore, we must create both a new promise
			// 	//and a new future, and add a new object from the buffer
			// 	//queue into the new promise.
			// 	if (!futureNotTaken) {
			// 		prepareNewFuture();
			// 		promise = std::promise<Object>();
			// 		future = promise.get_future();
			// 		promise.set_value(std::move(retObjs.pop()));
			// 		promiseSet = 1;
					
			// 	}
			// }

		}
	}

public:

	/** 
	 * Construct a new AsyncObjectStream
	 */
	AsyncObjectStream() : 
		finished(0),
		promiseSet(0)
	{
		future = promise.get_future();
	}

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
		bool noMoreObjects;
		bool valid;
		synchronized(m_mutex) {
			valid = future.valid();
			noMoreObjects = finished && (!valid || !promiseSet) && retObjs.empty();
		}
		if (noMoreObjects) {
			return false;
		}
		if (!valid) {
			return false;
		}

		//WARNING! Can *not* wait within a synchronized block.
		//But, if we get this far (i.e. valid == true) then we know
		//that the future object is valid and simply waiting for
		//the promise to publish a result. Hence, there is no chance
		//that the future object itself could change or become invalid
		//since only the consumer thread calls future.get(); 
		future.wait();
		bool caught;
		synchronized(m_mutex) {
			try {
				res = future.get();
				caught = false;
			} catch (std::logic_error e) {
				caught = true;
			}
		}
		if (caught) {
			return false;
		}

		
		prepareNewFuture();
		tryPrepareNextResult();
		return true;
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
		synchronized(m_mutex) {
			retObjs.push(std::move(res));
		}
		tryPrepareNextResult();
	}

	void addResult(Object&& res) {	
		if (finished) {
			throw std::runtime_error("AsyncObjectStream declared finished by tried to add new result.");
		}
		synchronized(m_mutex) {
			retObjs.push(std::move(res));
		}
		tryPrepareNextResult();
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
		synchronized(m_mutex) {
			finished = true;
			if (retObjs.empty() && !promiseSet) {
				try{
					promise.set_exception(
					std::make_exception_ptr(
						std::logic_error("Cancel waiting; we are done.")));
					promiseSet = 1;
				} catch (std::future_error& fe) {
					std::cerr << "AsyncObjectStream encountered an error cleaning its threads" << std::endl;
					exit(1);
				}
			}
		}
	}

};




#endif