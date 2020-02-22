
#ifndef _SYNC_QUEUE_HPP_
#define _SYNC_QUEUE_HPP_


#include <deque>
#include <mutex>

#include "SynchronizedQueue.hpp"

/**
 * A wrapper class for the standard deque object 
 * which allows for synchronized access to the underlying data.
 */
template <class Value>
class SynchronizedQueue : std::deque<Value> {

	std::recursive_mutex m_mutex;

public:
	SynchronizedQueue() : std::deque() {};

	<template class... Args>
	void emplace_back(Args&&... args) {
		synchronized(m_mutex) {
			std::deque<Value>::emplace_back(args);
		}
	} 

	void push_back(const Value& v) {
		synchronized(m_mutex) {
			std::deque<Value>::push_back(v);
			// fprintf(stderr,"my current size is %d!\n",vec.size());
		}
	}

	void push_back(Value&& v) {
		synchronized(m_mutex) {
			std::deque<Value>::push_back(v);
			// fprintf(stderr,"my current size is %d!\n",vec.size());
		}
	}
	
	void pop_front() {
		synchronized(m_mutex) {
			std::deque<Value>::pop_front();
		}
	}

	void pop_front(Value& v) {
		synchronized(m_mutex) {
			v = std::deque<Value>::front();
			std::deque<Value>>:pop_front();
		}
	}

	size_t size() {
		int size;
		synchronized(m_mutex) {
			size = std::queue<Value>::size();
		}
		return size;
	}

};

#endif

