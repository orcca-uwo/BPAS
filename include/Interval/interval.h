#ifndef _INTERVAL_H_
#define _INTERVAL_H_

#include "../DyadicRationalNumber/globals.h"
#include "../Symbol/Symbol.hpp"
#include <iostream>

/**
 * Data Structure for interval [a, b].
 */
class Interval {
	public:
		lfixq left;
		lfixq right;
};

/**
 * Interval addition
 * @res = @a + @b
 **/
inline void intervalAddition(Interval* res, Interval* a, Interval* b) {
	res->left = a->left + b->left;
	res->right = a->right + b->right;
}

/**
 * Interval multiplication
 * @res = @a * @b
 **/
void intervalMultiplication(Interval* res, Interval* a, Interval* b);

/** Interval lists for real roots of multivariate polynomials.
 * i.e. []*..*[], .., []*..*[]
 **/
class Intervals {
	private:
		int var;			// Number of variates
		Symbol* names;		// Variable names
		std::vector<Interval> roots;	// Roots, using one dimension
						// to present two dimensions
						// []*..*[], .., []*..*[]
						// Each root is in the ascending
						// order of the weight of variates

	public:
		/**
		 * Construct intervals
		 *
		 * @param
		 **/
		Intervals() : var(1) {
			names = new Symbol[2];
			names[0] = "1";
			names[1] = "_1";
		}
		/**
		 * Construct intervals with number of variates
		 *
		 * @param v: Number of variates
		 **/
		Intervals(int v) : var(v) {
			names = new Symbol[var+1];
			names[0] = "1";
			for (int i = 1; i <= var; ++i) {
				std::ostringstream convert;
				convert << var - i + 1;
				names[i] = "_";
				names[i] += convert.str();
			}
		}
		/**
		 * Copy construct
		 *
		 * @param b: Intervals
		 **/
		Intervals (const Intervals& b) : var(b.var), roots(b.roots) {
			names = new Symbol[var+1];
			std::copy(b.names, b.names+var+1, names);
		}
		/**
		 * Destroy intervals
		 *
		 * @param
		 **/
		~Intervals() {
			delete [] names;
			roots.clear();
		}
		/**
		 * Overload operator =
		 *
		 * @param b: Intervals
		 **/
		inline Intervals& operator= (Intervals b) {
			if (this != &b) {
				delete [] names;
				var = b.var;
				roots = b.roots;
				names = new Symbol[var+1];
				std::copy(b.names, b.names+var+1, names);
			}
			return *this;
		}
		/**
		 * Get number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() {
			return var;
		}

		/**
		 * Get the number of roots
		 *
		 * @param
		 **/
		inline int numberOfIntervals() {
			return roots.size() / var;
		}

		/**
		 * Push an interval to the end of the list
		 *
		 * @param pI: An interval
		 **/
		inline void pushInterval(Interval pI) {
			roots.push_back(pI);
		}
		/**
		 * Push an interval to the end of the list
		 *
		 * @param l: Interval's left
		 * @param r: Interval's right
		 **/
		inline void pushInterval(lfixq l, lfixq r) {
			Interval elem;
			elem.left = l;
			elem.right = r;
			roots.push_back(elem);
		}
		inline void pushInterval(double l, double r) {
			Interval elem;
			elem.left = l;
			elem.right = r;
			roots.push_back(elem);
		}
		/**
		 * Pop the last interval
		 *
		 * @param
		 **/
		inline void popInterval() {
			if (roots.size())
				roots.pop_back();
		}
		/**
		 * Clear all the intervals
		 *
		 * @param
		 **/
		inline void clear() {
			roots.clear();
		}

		/**
		 * Get a root of a variate
		 *
		 * @param k: Offset of roots
		 * @param l: Offset of variables in an ascending order
		 **/
		inline Interval* interval(int k, int l=0) {
			return &roots[k*var+l];
		}
		/**
		 * Get a root of variates
		 *
		 * @param k: Offset of roots
		 **/
		inline Intervals root(int k) {
			Intervals pIs(var);
			for (int i = 0; i < var; ++i)
				pIs.roots.push_back(roots[k*var+i]);
			return pIs;
		}
		/**
		 * Set variable names
		 *
		 * @param xs: Variable names
		 **/
		inline void setVariableNames (const std::vector<Symbol>& xs) {
			int ns = xs.size();
			if (ns >= var) {
				names[0] = "9";
				for (int i = var, j = 0; i > 0; --i, ++j)
					names[i] = xs[j];
			}
			else
				std::cout << "BPAS: warning, not enough variables to set in Intervals." << std::endl;
		}
		/**
		 * Copy the k-th root
		 *
		 * @param pIs: A list of roots
		 * @param k: Offset of roots
		 **/
		void copyFrom(Intervals& pIs, int k);
		/**
		 * Scale to 2^k
		 *
		 * @param k: The power of 2
		 **/
		void scale(int k);
		/**
		 * Scale to that over 2
		 *
		 * @param
		 **/
		void transformLeft();
		/**
		 * Scale to that plus 1 and then over 2
		 *
		 * @param
		 **/
		void transformRight();
		/**
		 * Negate each interval and re-sort
		 *
		 * @param
		 **/
		void negate();
		/**
		 * Merge another intervals into the end of the current one
		 *
		 * @param pIs: A list of intervals
		 **/
		void concatenate(Intervals& pIs);

		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param a: Intervals
		 **/
		friend std::ostream& operator<< (std::ostream &out, Intervals& a);
		bool isExactEqual(Intervals& pIs);
};


/**
 * Swap two intervals
 **/
inline void swapInterval(Interval* a, Interval* b) {
	Interval tmp;
	tmp.left = a->left;
	tmp.right = a->right;
	a->left = b->left;
	a->right = b->right;
	b->left = tmp.left;
	b->right = tmp.right;
}

#endif
