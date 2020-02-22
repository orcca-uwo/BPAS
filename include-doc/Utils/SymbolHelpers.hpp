#ifndef _SYMBOL_HELPERS_H_
#define _SYMBOL_HELPERS_H_

#include <vector>
#include <algorithm>
#include "../Symbol/Symbol.hpp"

/// Misc Helper functions ///

static bool contains(std::vector<Symbol>& v,const Symbol& a) {
	int aIndex = std::find(v.begin(),v.end(),a) - v.begin();
	if (aIndex == v.size())
		return false;
	else
		return true;
}

static void switchVariables(const Symbol& a, const Symbol& b, std::vector<Symbol>& v) {
	int aIndex = std::find(v.begin(),v.end(),a) - v.begin();
	int bIndex = std::find(v.begin(),v.end(),b) - v.begin();
	Symbol temp;
	temp = std::move(v[aIndex]);
	v[aIndex] = std::move(v[bIndex]);
	v[bIndex] = std::move(temp);
}

static void moveSymbolToBegin(std::vector<Symbol>& v, const Symbol& a) {
	int aIndex = std::find(v.begin(),v.end(),a) - v.begin();
	if (aIndex == v.size()) {
		std::cerr << "BPAS: error, input vector does not contain the input symbol" << std::endl;
		exit(1); 
	}
	v.erase(v.begin()+aIndex);
	v.insert(v.begin(),a);
}

static void moveInitialSymbolToEnd(std::vector<Symbol>& v) {
	if (!v.empty()) {
		Symbol s(v[0]);
		v.erase(v.begin());
		v.push_back(s);
	}
}

static void printVariables(std::vector<Symbol>&& vv) {
	std::vector<Symbol> v = vv;
	for (int i=0;i<v.size();++i)
		std::cout << "v[" << i << "] = " << v[i] << std::endl;
}

static void printVariables(const std::vector<Symbol>& vv) {
	std::vector<Symbol> v = vv;
	for (int i=0;i<v.size();++i)
		std::cout << "v[" << i << "] = " << v[i] << std::endl;
}

static void printVariables(std::vector<Symbol>&& vv, std::string name) {
	std::vector<Symbol> v = vv;
	for (int i=0;i<v.size();++i)
		std::cout << name << "[" << i << "] = " << v[i] << std::endl;
}

static void printVariables(const std::vector<Symbol>& vv, std::string name) {
	std::vector<Symbol> v = vv;
	for (int i=0;i<v.size();++i)
		std::cout << name << "[" << i << "] = " << v[i] << std::endl;
}

/// Set operations on vectors of Symbols ///

template <typename T>
static bool isSubset(const std::vector<T>& a, const std::vector<T>& b) {
	std::vector<Symbol> v1(a),v2(b);
	
	std::sort(v1.begin(), v1.end());
	std::sort(v2.begin(), v2.end());
	return std::includes(v2.begin(), v2.end(), v1.begin(), v1.end());
}

template <typename T>
static bool isAMemberOf(const T& a, const std::vector<T>& b) {
	return (std::find(b.begin(),b.end(),a) != b.end());
}

// TODO: THESE FUNCTIONS MAY BE INEFFICIENT
static std::vector<Symbol> setUnion(const std::vector<Symbol>& a, const std::vector<Symbol>& b) {
	std::vector<Symbol> v1(a),v2(b),v3;
	
	std::sort(v1.begin(), v1.end());
	std::sort(v2.begin(), v2.end());
	std::set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
	
	return v3;
}

static std::vector<Symbol> setIntersection(const std::vector<Symbol> &a, const std::vector<Symbol> &b) {
    std::vector<Symbol> v1(a),v2(b),v3;
    
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
    
    return v3;
  }

static std::vector<Symbol> setDifference(const std::vector<Symbol> &a, const std::vector<Symbol> &b) {
    std::vector<Symbol> v1(a),v2(b),v3;
    
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
    
    return v3;
  }
  
/*
*
* Find all Symbols in v1 that do not appear in v2 and return the intersection v1 ∩ v2 with elements in the same order
* as they appear in v1.
*
* @param v1: vector from which to remove elements
* @param v2: vector of elements to (potentially) remove
*
*/
// TODO: Make this more efficient if possible
static std::vector<Symbol> orderPreservingSetIntersection(const std::vector<Symbol> &v1, const std::vector<Symbol> &v2) {
	std::vector<Symbol> v3,ret;
	
	v3 = setIntersection(v1,v2);
	// reserve space for ret?
	std::vector<Symbol>::iterator it;
	for (int i=0; i<v1.size(); ++i) {
		if (std::find(v3.begin(),v3.end(),v1[i]) != v3.end())
			ret.push_back(v1[i]);
	}
	return ret;
}
  
/*
*
* Remove all Symbols in v1 that appear in v2 and return the elements in the same order
* as they appear in v1.
*
* @param v1: vector from which to remove elements
* @param v2: vector of elements to (potentially) remove
*
*/
// TODO: Make this more efficient if possible
static std::vector<Symbol> orderPreservingSetDifference(const std::vector<Symbol> &v1, const std::vector<Symbol> &v2) {
	std::vector<Symbol> v3,ret;
	
	v3 = setDifference(v1,v2);
	// reserve space for ret?
	std::vector<Symbol>::iterator it;
	for (int i=0; i<v1.size(); ++i) {
		if (std::find(v3.begin(),v3.end(),v1[i]) != v3.end())
			ret.push_back(v1[i]);
	}
	return ret;
}
  
/*
*
* Collect all Symbols in v1 that appear in v2, preserving the order of elements in v1 followed by the
* elements of v2 \ v1 in the same order they appear in v2.
*
* @param v1: vector from which to remove elements
* @param v2: vector of elements to (potentially) remove
*
*/
// TODO: Make this more efficient if possible
static std::vector<Symbol> orderPreservingSetUnion(const std::vector<Symbol> &v1, const std::vector<Symbol> &v2) {
	std::vector<Symbol> ret(v1);
	printVariables(ret);
	
	// reserve space for ret?
	std::vector<Symbol>::iterator it;
	for (int i=0; i<v2.size(); ++i) {
		if (std::find(v1.begin(),v1.end(),v2[i]) == v1.end()) {
			ret.push_back(v2[i]);
		}
	}
	return ret;
}

#endif
