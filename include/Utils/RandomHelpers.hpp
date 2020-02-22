#ifndef _RANDOM_HELPERS_HPP_
#define _RANDOM_HELPERS_HPP_

//#include "../ring.h"
//#include "../RingPolynomial/upolynomial.h"
//#include "../Ring/RationalNumber.hpp"
#include <vector>
#include <assert.h>

#include "RandomHelpers.h"


/**
  * Generate a fixed number of randomly selected values in a range without duplicates
  *
  * @param low: lower bound of the range
  * @param high: upper bound of the range
  * @param numElems: number of values in the range;
  **/
static std::vector<int> randValsInRange(int low,int high,int numElems) {
	assert (numElems <= (high-low+1));
	int index[numElems]={0},n;
	std::vector<int> out;
	bool check(false);
	
	for (auto i=0; i<numElems; ++i) {
		if (i==0)
			index[i] = randValInRange(low,high);
		else {
			while (check==false) {
				check = true;
				n = randValInRange(low,high);
				for (auto j=0; j<i; ++j) {
					if (index[j] == n)
						check = false;
				}
			}
			index[i] = n;
			check = false;
		}
	}
	for (auto i=0; i<numElems; ++i) {
		out.push_back(index[i]);
	}
	return out;
}

#endif
