#ifndef _QREGULARCHAIN_H_
#define _QREGULARCHAIN_H_

/**
 * Data Structure for regular chains
 * using univariate & multivariate rational polynomial
 **/

#include "../DyadicRationalNumber/globals.h"
#include "../RationalNumberPolynomial/urpolynomial.h"
#include "../RationalNumberPolynomial/mrpolynomial.h"

class RationalRegularChain {
	private:
		int var;					// Number of variables
		DenseUnivariateRationalPolynomial p;		// Univariate
		SparseMultivariateRationalPolynomial* mps;	// Multivariates
		Symbol* names;				// Variable names, names[0] < names[1] < ...

		/* Real Root Isolation */
		int rcRealRootIsolation(Intervals*, lfixq, int);
		void refineMultivariateInterval(Intervals*, lfixq, int);
		void isolatePositiveMultivariateRealRoots(Intervals*, Intervals*, int, int, lfixq, int);
		void isolateMultivariateRealRoots(Intervals*, Intervals*, int, lfixq, int ts);
	
	public:
		/**
		 * Default constructor
		 *
		 * @param
		 **/ 
		RationalRegularChain () : var(0) {}
		/**
		 * Constructor with number of variables and variable names
		 *
		 * @param v: Number of variables
		 * @param xs: Variable names
		 **/
		RationalRegularChain (int v, Symbol* xs) {
			if (v > 0) {
				var = v;
				if (v > 1)
					mps = new SparseMultivariateRationalPolynomial[var-1];
				if (xs != NULL) {
					names = new Symbol[var];
					for (int i = 0, j = var-1; i < v; ++i, --j)
						names[j] = xs[i];
				}
				else {
					std::cout << "BPAS: error, RationalRegularChain must be defined with variables." << std::endl;
					exit(1);
				}
			}
			else { var = 0; }
		}
		/**
		 * Copy constructor
		 *
		 * @param rc: A regular chain
		 **/ 
		RationalRegularChain (const RationalRegularChain& rc) : var(rc.var), p(rc.p) {
			mps = new SparseMultivariateRationalPolynomial[var-1];
			std::copy(rc.mps, rc.mps+var-1, mps);
			names = new Symbol[var];
			std::copy(rc.names, rc.names+var, names);
		}
		/**
		 * Destructor
		 *
		 * @param
		 **/
		~RationalRegularChain () {
			if (var > 1) { delete [] mps; }
			if (names != NULL) { delete [] names; }
		}
		
//		/**
//		 * Copy an object derived from abstract BPASTriangularSet class to type of current object
//		 *
//		 * @param ts: triangular set to copy
//		 **/
//		inline void copy(const BPASTriangularSet<lfixq,SparseMultivariateRationalPolynomial>& ts) override {
//			if (dynamic_cast<const RationalRegularChain*>(&ts))
//				*this = *dynamic_cast<const RationalRegularChain*>(&ts);
//			else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RationalRegularChain."));
//		}
		
		/**
		 * Number of variables
		 *
		 * @param
		 **/
		int numberOfVariables() const {
			return var;
		}
		/**
		 * Variables
		 *
		 * @param
		 **/
		std::vector<Symbol> variables() const {
			std::vector<Symbol> xs;
			for (int i = var-1; i > -1; --i)
				xs.push_back(names[i]);
			return xs;
		}

		/**
		 * Overload operator =
		 *
		 * @param rc: A regular chain
		 **/
		inline RationalRegularChain& operator= (RationalRegularChain rc) {
			if (this != &rc) {
				var = rc.var;
				p = rc.p;
				delete [] mps;
				mps = new SparseMultivariateRationalPolynomial[var-1];
				std::copy(rc.mps, rc.mps+var-1, mps);
				delete [] names;
				names = new Symbol[var];
				std::copy(rc.names, rc.names+var, names);
			}
			return *this;
		}
		
//		/**
//		 * Assignment operator =
//		 *
//		 * @param a: A BPASTriangularSet
//		 **/
//		inline BPASTriangularSet<lfixq,SparseMultivariateRationalPolynomial>& operator= (const BPASTriangularSet<lfixq,SparseMultivariateRationalPolynomial>& a) override {
//			if (dynamic_cast<const RationalRegularChain*>(&a))
//				*this = dynamic_cast<const RationalRegularChain&>(a);
//			else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RationalRegularChain."));
//			return *this;
//		}

//		//TODO: Need to implement this function as a move, not a copy.
//		/**
//		 * Overload move operator =
//		 *
//		 * @param rc: A regular chain
//		 **/
//		inline RationalRegularChain& operator= (RationalRegularChain&& rc) {
//			if (this != &rc) {
//				var = rc.var;
//				p = rc.p;
//				delete [] mps;
//				mps = new SparseMultivariateRationalPolynomial[var-1];
//				std::copy(rc.mps, rc.mps+var-1, mps);
//				delete [] names;
//				names = new Symbol[var];
//				std::copy(rc.names, rc.names+var, names);
//			}
//			return *this;
//		}
		
		/**
//		 * Move assignment operator =
//		 *
//		 * @param a: A BPASTriangularSet
//		 **/
//		inline BPASTriangularSet<lfixq,SparseMultivariateRationalPolynomial>& operator= (const BPASTriangularSet<lfixq,SparseMultivariateRationalPolynomial>&& a) override {
//			if (dynamic_cast<const RationalRegularChain*>(&a))
//				*this = dynamic_cast<const RationalRegularChain&&>(a);
//			else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RationalRegularChain."));
//			return *this;
//		}

		/**
		 * Overload operator +
		 *
		 * @param up: A univariate polynomial
		 **/
		inline RationalRegularChain operator+ (DenseUnivariateRationalPolynomial up) {
			RationalRegularChain r (*this);
			return (r += up);
		}
		/**
		 * Overload operator +=
		 *
		 * @param up: A univariate polynomial
		 **/
	        RationalRegularChain& operator+= (DenseUnivariateRationalPolynomial up);
		/**
		 * Overload operator +
		 *
		 * @param mp: A multivariate polynomial
		 **/
		inline RationalRegularChain operator+ (SparseMultivariateRationalPolynomial mp) {
			RationalRegularChain r (*this);
			return (r += mp);
		}
		/**
		 * Overload operator +=
		 *
		 * @param mp: A multivariate polynomial
		 **/
	        RationalRegularChain& operator+= (SparseMultivariateRationalPolynomial mp);

		/**
		 * Select a polynomial given the leading variable
		 * 
		 * @param x: The leading variable name
		 **/
        SparseMultivariateRationalPolynomial select(const Symbol& x) const;

		/**
		 * Regular sub-chain with respect to variable s:
		 * Returns the rc consisting of polynomials with 
		 * main variable strictly less than s
		 *
		 * @param s: Symbol of the main variable of specified element of the rc
		 **/
		 void lower(const Symbol& x, RationalRegularChain& ts) const;

		/**
//		 * Returns the rc consisting of polynomials with
//		 * main variable strictly greater than s
//		 *
//		 * @param s: Symbol of the main variable of specified element of the rc
//		 **/
//		 // TODO: NOT IMPLEMENTED YET
//		inline void upper(const Symbol& x, BPASTriangularSet& ts) const {
//			ts = RationalRegularChain();
//		}

		/**
		 * Returns the normal form of p with respect to rc.
		 *
		 * @param p: Polynomial to be normalized modulo the rc
		 **/
		SparseMultivariateRationalPolynomial normalize (const SparseMultivariateRationalPolynomial& p);
		
		/**
		 * Real root isolation
		 *
		 * @param width: The interval's width
		 * @param ts: Taylor Shift option
		 **/
		inline Intervals realRootIsolate (mpq_class width, int ts=-1) {
			Intervals pIs;
			bpas_root_width = width;
			rcRealRootIsolation(&pIs, width, ts);
			pIs.setVariableNames(variables());
			return pIs;
		}

		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param rc: A regular chain
		 **/ 
        friend std::ostream& operator<< (std::ostream &out, RationalRegularChain& rc);
};

/**
 * Real root isolation on a vector of RationalRegularChain
 *
 * @param chains, defined as std::vector<RationalRegularChain> chains (<size>, RationalRegularChain(<# of variables>, <variable names>));
 * @param width: Interval's width
 * @param ts: Taylor shift option
 **/
Intervals RealRootIsolation (std::vector<RationalRegularChain> chains, lfixq width, int ts=-1);

#endif
