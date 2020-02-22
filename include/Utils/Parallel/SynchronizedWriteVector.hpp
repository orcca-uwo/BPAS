
#ifndef _SYNC_WRITE_VECTOR_HPP_
#define _SYNC_WRITE_VECTOR_HPP_


#include <vector>
#include <mutex>

#include "Synchronized.hpp"

template <class Value>
class SynchronizedWriteVector {
	
protected:
	std::vector<Value> vec;
	std::recursive_mutex m_mutex;

public:
	SynchronizedWriteVector() {};

	SynchronizedWriteVector(int n) : vec(n) {};

	size_t size() const {
		return vec.size();
	}

	bool empty() const {
		return vec.empty();
	}

	void push_back(const Value& v) {
		synchronized(m_mutex) {
			vec.push_back(v);
			// fprintf(stderr,"my current size is %d!\n",vec.size());
		}
	}

	void push_back(Value&& v) {
		synchronized(m_mutex) {
			vec.push_back(v);
			// fprintf(stderr,"my current size is %d!\n",vec.size());
		}
	}
	
	typename std::vector<Value>::const_iterator begin() {
		return vec.begin();
	}

	typename std::vector<Value>::const_iterator end() {
		return vec.end();
	}

	void insert(typename std::vector<Value>::const_iterator insertPos, typename std::vector<Value>::const_iterator start, typename std::vector<Value>::const_iterator end) {
		synchronized(m_mutex) {
			vec.insert(insertPos,start,end);
		}
	}
	
	void clear() {
		synchronized(m_mutex) {
			vec.clear();
		}
	}
	
	void reserve(size_t size) {
		synchronized(m_mutex) {
			vec.reserve(size);
		}
	}

	const Value& operator[](size_t pos) const {
		return vec[pos];
	}
	
	const std::vector<Value>& vector() {
		return vec;
	}
	
	void moveVectorIn(std::vector<Value>&& in) {
		synchronized(m_mutex) {
			vec.clear();
			vec = in;
		}
	}

	std::vector<Value> moveVectorOut() {
		return std::move(vec);
	}

};

#endif

