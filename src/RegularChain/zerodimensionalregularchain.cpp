
#include "RegularChain/regularchain.hpp"
#include "RegularChain/zerodimensionalregularchain.hpp"
#include "RationalNumberPolynomial/mrpolynomial.h"
#include "Utils/SymbolHelpers.hpp"

int ZDisInvertibleDepth = 0;
int ZDregularGCDDepth = 0;
int ZDregularizeDepth = 0;
int ZDregularizeInitialDepth = 0;
int ZDlnzsrDepth = 0;
int ZDlnzsrInnerDepth = 0;
int ZDdepth = 0;

/// Constructors ///

/**
 * Default constructor: creates an empty triangular set of variable size
 * with empty list of transcendentals
 *
 * @param
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain () : RegularChain<Field,RecursivePoly>() {}

/**
 * Construct an empty triangular set of variable size with
 * empty list of variables and list of transcendental variables given by ps
 *
 * @param ps: The variable names
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const std::vector<Symbol>& ps) : RegularChain<Field,RecursivePoly>() {
	trcVars = ps;
}

template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const std::vector<RecursivePoly> polys) : RegularChain<Field,RecursivePoly>(polys) {

}

/**
 * Construct a variable triangular set containing p, such that the variables of p are treated as algebraic,
 * with empty list of transcendentals
 *
 * @param p: The polynomial to add
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const RecursivePoly& p) : RegularChain<Field,RecursivePoly>(p) {
	if (p.numberOfVariables() != 1) {
		std::cerr << "BPAS: error, cannot construct a ZeroDimensionalRegularChain from a polynomial with more than one variable." << std::endl;
		exit(1);
	}
}

/**
 * Construct a variable triangular set containing p, such that the variables in ts are
 * treated as transcendental, while any remaining variables of p are treated as algebraic
 *
 * @param p: The polynomial to add
 * @param ts: The transcendental variable names
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const RecursivePoly& p, const std::vector<Symbol>& ts) : RegularChain<Field,RecursivePoly>(p,ts) {
	std::vector<Symbol> vs;
	vs.push_back(p.leadingVariable());
	vs.insert(vs.end(),ts.begin(),ts.end());
	vs = orderPreservingSetIntersection(vs,p.variables());
	// TODO: consider sorting before the non-equality test
	if (p.variables() != vs) {
		std::cerr << "BPAS: error, cannot construct a ZeroDimensionalRegularChain from a polynomial with more than one (non-transcendental) variable." << std::endl;
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//			std::cerr << "p.variables:" << std::endl;
//			printVariables(p.variables());
//			std::cerr << "p main variable plus transcendentals in p:" << std::endl;
//			printVariables(vs);
		#endif
		exit(1);
	}
}

/**
 * Copy constructor
 *
 * @param a: A zero dimensional regular chain
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) : RegularChain<Field,RecursivePoly>(a) {
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "copy constructor in ZDRC called..." << std::endl;
	#endif
}

/**
 * Copy constructor
 *
 * @param a: A regular chain
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const RegularChain<Field,RecursivePoly>& a, int options) : RegularChain<Field,RecursivePoly>(a) {
	unsigned long long int timerStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&timerStart);
	#endif
	if ((options & ASSUME_MAKESCHAIN) != ASSUME_MAKESCHAIN) {
		int firstZeroIndex(-1);
		for (int i=vars.size()-1; i>-1; --i) {
			if (set[i].isZero())
				firstZeroIndex = i;
		}
		if (!this->isEmpty()) {
			if (firstZeroIndex == vars.size()-1) {
				std::cerr << "BPAS: error, cannot cast from a positive dimensional regular chain to a zero dimensional regular chain." << std::endl;
				exit(1);
			}
			int secondNonzeroIndex(-1);
			if (firstZeroIndex != -1) {
				for (int i=firstZeroIndex-1; i>-1; --i) {
					if (!set[i].isZero())
						secondNonzeroIndex = i;
				}
			}
			if (secondNonzeroIndex != -1) {
				std::cerr << "BPAS: error, cannot cast from a positive dimensional regular chain to a zero dimensional regular chain." << std::endl;
				exit(1);
			}
		}
		if (firstZeroIndex != -1) {
			if (mode == TS_FIXED) {
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "BPAS: warning, treating a positive dimensional regular chain as zero dimensional." << std::endl;
				#endif
			}
			else {
				vars.erase(vars.begin(),vars.begin()+firstZeroIndex+1);
				set.erase(set.begin(),set.begin()+firstZeroIndex+1);
				if (vars != algVars) {
					std::cerr << "BPAS: error, vars must equal algVars in ZDRC(RC,opts) copy constructor." << std::endl;
					exit(1);
				}
			}
		}
	}
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&timerStart,&zdrcCopyTime);
	#endif
}

/**
 * Move constructor
 *
 * @param a: An r-value reference zero dimensional regular chain
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a) : RegularChain<Field,RecursivePoly>(std::move(a)) {
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "move constructor in ZDRC called..." << std::endl;
	#endif
}

/**
 * Move constructor
 *
 * @param a: An r-value reference regular chain
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (RegularChain<Field,RecursivePoly>&& a, int options) : RegularChain<Field,RecursivePoly>(std::move(a)) {
	if ((options & ASSUME_MAKESCHAIN) != ASSUME_MAKESCHAIN) {
		int firstZeroIndex(-1);
		for (int i=vars.size()-1; i>-1; --i) {
			if (set[i].isZero())
				firstZeroIndex = i;
		}
		if (!this->isEmpty()) {
			if (firstZeroIndex == vars.size()-1) {
				std::cerr << "BPAS: error, cannot cast from a positive dimensional regular chain to a zero dimensional regular chain." << std::endl;
				exit(1);
			}
			int secondNonzeroIndex(-1);
			if (firstZeroIndex != -1) {
				for (int i=firstZeroIndex-1; i>-1; --i) {
					if (!set[i].isZero())
						secondNonzeroIndex = i;
				}
			}
			if (secondNonzeroIndex != -1) {
				std::cerr << "BPAS: error, cannot cast from a positive dimensional regular chain to a zero dimensional regular chain." << std::endl;
				exit(1);
			}
		}
		if (firstZeroIndex != -1) {
			if (mode == TS_FIXED) {
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "BPAS: warning, treating a positive dimensional regular chain as zero dimensional." << std::endl;
				#endif
			}
			else {
				vars.erase(vars.begin(),vars.begin()+firstZeroIndex+1);
				set.erase(set.begin(),set.begin()+firstZeroIndex+1);
				if (vars != algVars) {
					std::cerr << "BPAS: error, vars must equal algVars in ZDRC(RC,opts) move constructor." << std::endl;
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//						printVariables(vars,"vars");
//						printVariables(algVars,"algVars");
//						ZeroDimensionalRegularChain<Field,RecursivePoly> T(*this);
						std::cerr << "*this = " << *this << std::endl;
					#endif
					exit(1);
				}
			}
		}
	}
}

/**
 * Computational constructor: creates a triangular set given all the data
 *
 * @param vs: variables of the triangular set
 * @param avs: algebraic variables of the triangular set
 * @param tvs: transcendental variables of the triangular set
 * @param polys: polynomials of the triangular set
 * @param tsm: whether the triangular set is variable or fixed
 * @param c: characteristic of the triangular set
 * @param normal: whether the triangular set is strongly normalized
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>::ZeroDimensionalRegularChain (const std::vector<Symbol>&& vs, const std::vector<Symbol>&& avs, const std::vector<Symbol>&& tvs, const std::vector<RecursivePoly>&& ts, TriangularSetMode tsm, const mpz_class& c) : RegularChain<Field,RecursivePoly>(std::move(vs),std::move(avs),std::move(tvs),std::move(ts),tsm,c) {}



/// Member Functions ///
		
/**
 * Construct a zero-dimensional regular chain from the current object and an input polynomial
 *
 * @param p: A polynomial
 **/
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::constructChain(const RecursivePoly& p, int options) {
	Symbol v(p.leadingVariable());
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "entering constructChain(p,opts) in ZDRC." << std::endl;
		std::cerr << "opts = " << options << std::endl;
	#endif
	
	if ((options & ASSUME_REGULAR) != ASSUME_REGULAR) {
		// perform a regularity test of p.initial() modulo Sat(T)
		if (!this->isRegular(p.initial())) {
				std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "in constructChain(p,opts) of ZDRC." << std::endl;
				#endif
				exit(1);
		}
	}
	if ((options & ASSUME_ZERODIMENSIONAL) != ASSUME_ZERODIMENSIONAL) {
		std::vector<Symbol> vs,vs2;
		vs2.push_back(v);
	
		if (isAMemberOf(v,algVars)) {
			std::cerr << "BPAS: error, the leading variable of p is an algebraic variable of the current object, so cannot extend it." << std::endl;
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "p = " << p << std::endl;
				std::cerr << "v = " << v << std::endl;
//				printVariables(vars,"vars");
			#endif
			exit(1);
		}
		
		vs = setDifference(p.variables(),this->allVariables());
//			vs = setDifference(p.variables(),vars);
//			vs = setDifference(vs,trcVars);
	
		if (mode == TS_VARIABLE) {
			if (vs != vs2) {
				std::cerr << "BPAS: error, cannot add p to current object because resulting object would not be zero-dimensional." << std::endl;
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//					std::cerr << "p.variables:" << std::endl;
//					printVariables(p.variables());
//					std::cerr << "vars:" << std::endl;
//					printVariables(vars);
//					std::cerr << "trcVars:" << std::endl;
//					printVariables(trcVars);
				#endif
				exit(1);
			}
		}
		else {
			if (!vs.empty()) {
				std::cerr << "BPAS: error, cannot add a polynomial with new variables to a fixed mode triangular set." << std::endl;
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//					std::cerr << "p.variables:" << std::endl;
//					printVariables(p.variables());
//					std::cerr << "vars:" << std::endl;
//					printVariables(vars);
//					std::cerr << "trcVars:" << std::endl;
//					printVariables(trcVars);
				#endif
				exit(1);
				
			}
		}
	}
	
	// make current object temporarily one dimensional to make space for new polynomial
	if (!isAMemberOf(v,vars)) {
		vars.insert(vars.begin(),v);
		set.insert(set.begin(),RecursivePoly());
	}
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "calling constructChain(p,opts) from ZeroDimensionalRegularChain class" << std::endl;
	#endif
	RegularChain<Field,RecursivePoly>::constructChain(p,options);
}
		
/**
 * Construct a zero-dimensional regular chain from the current object and an input upper regular chain
 *
 * @param T: An upper regular chain
 **/
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::constructChain(const RegularChain<Field,RecursivePoly>& T, int options) {
	// TODO: check that the chain is upper
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "entering constructChain(p,opts) in ZDRC." << std::endl;
		std::cerr << "opts = " << options << std::endl;
	#endif
	std::vector<RecursivePoly> polys(T.polynomials());
	for (int i=polys.size()-1; i>-1; --i) {
		if (!polys[i].isZero()) {
			constructChain(polys[i],options);
		}
	}
}

/**
 * Construct a set of regular chains from an input regular chain and polynomial
 *
 * @param T: A regular chain
 * @param p: A polynomial
 **/
template <class Field, class RecursivePoly>
std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> ZeroDimensionalRegularChain<Field,RecursivePoly>::constructChains(const RecursivePoly& p, int options) const {
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "entering ZDRC constructChains(p,opts): " << p << std::endl;
		std::cerr << "this = " << *this << std::endl;
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);
	std::vector<RegularChain<Field,RecursivePoly>> rcResults;
	std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> results;
	Symbol v(p.leadingVariable());
	
	
	if ((options & ASSUME_REGULAR) != ASSUME_REGULAR) {
		// perform a regularity test of p.initial() modulo Sat(T)
		if (!this->isRegular(p.initial())) {
				std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "in constructChains(p,opts) of ZDRC." << std::endl;
				#endif
				exit(1);
		}
	}
	
	if ((options & ASSUME_ZERODIMENSIONAL) != ASSUME_ZERODIMENSIONAL) {
		std::vector<Symbol> vs,vs2;
		vs2.push_back(v);
	
		if (isAMemberOf(v,vars)) {
			std::cerr << "BPAS: error, the leading variable of p is an algebraic variable of the current object, so cannot extend it." << std::endl;
			exit(1);
		}
	
		vs = setDifference(p.variables(),vars);
		vs = setDifference(vs,trcVars);
		if (vs != vs2) {
			std::cerr << "BPAS: error, cannot add p to current object because resulting object would not be zero-dimensional." << std::endl;
			exit(1);
		}
	}
	
	options = options | ASSUME_MAKESCHAIN;
	
	// make root chain T temporarily one dimensional to make space for new polynomial
	std::vector<Symbol> vs;
	std::vector<RecursivePoly> polys;
	vs.insert(vs.begin(),v);
	vs.insert(vs.end(),vars.begin(),vars.end());
	polys.insert(polys.begin(),RecursivePoly());
	polys.insert(polys.end(),set.begin(),set.end());
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//		ZeroDimensionalRegularChain<Field,RecursivePoly> T(vs,algVars,trcVars,polys,TS_VARIABLE,this->isStronglyNormalized());
//		std::cerr << "T = " << T << std::endl;
	
		std::cerr << "calling constructChains(p,opts) from ZeroDimensionalRegularChain class" << std::endl;
	#endif
	RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&ZeroDimensionalRegularChain<Field,RecursivePoly>::constructChainsFromPoly), this, p, true, this->numberOfAlgebraicVariables(), options);
	RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, i);
		rcResults.push_back(currConstructChainsTask);
	}
//	rcResults = this->RegularChain<Field,RecursivePoly>::constructChains(p,true,this->numberOfAlgebraicVariables(),options);
	for (int i=0; i<rcResults.size(); ++i) {
		results.push_back(ZeroDimensionalRegularChain(std::move(rcResults[i]),options));
	}
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "leaving ZDRC constructChains(p,opts): " << p << std::endl;
	#endif
	return results;
}
		
/**
 * Construct a set of regular chains from the current object and an input regular chain above the current object
 *
 * @param T: A regular chain with no overlapping algebraic variables to those of the current object
 **/
template <class Field, class RecursivePoly>
std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> ZeroDimensionalRegularChain<Field,RecursivePoly>::constructChains(const RegularChain<Field,RecursivePoly>& T, int options) const {
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "entering RC constructChains(T,opts)" << std::endl;
	#endif
	std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> results,tasks,moreTasks;
	options = options | ASSUME_MAKESCHAIN;
	std::vector<RecursivePoly> polys(T.polynomials());
	results.push_back(*this);
	for (int i=polys.size()-1; i>-1; --i) {
		if (!polys[i].isZero()) {
			for (int j=0; j<results.size(); ++j) {
				moreTasks = results[j].constructChains(polys[i],options);
				tasks.insert(tasks.end(),moreTasks.begin(),moreTasks.end());
			}
			results = std::move(tasks);
			tasks.clear();
		}
	}
	return results;
}

/**
 * Assignment operator =
 *
 * @param a: A ZeroDimensionalRegularChain
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) {
	RegularChain<Field,RecursivePoly>::operator=(a);
	return *this;
}
		
/**
 * Assignment operator =
 *
 * @param a: A BPASTriangularSet
 **/
template <class Field, class RecursivePoly>
BPASTriangularSet<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (const BPASTriangularSet<Field,RecursivePoly>& a) {
	if (dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to ZeroDimensionalRegularChain."));
	return *this;
}

/**
 * Assignment operator =
 *
 * @param a: A BPASRegularChain
 **/
template <class Field, class RecursivePoly>
BPASRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (const BPASRegularChain<Field,RecursivePoly>& a) {
	if (dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASRegularChain to ZeroDimensionalRegularChain."));
	return *this;
}
		
/**
 * Assignment operator =
 *
 * @param a: A BPASZeroDimensionalRegularChain
 **/
template <class Field, class RecursivePoly>
BPASZeroDimensionalRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (const BPASZeroDimensionalRegularChain<Field,RecursivePoly>& a) {
	if (dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASZeroDimensionalRegularChain to ZeroDimensionalRegularChain."));
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A ZeroDimensionalRegularChain
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a) {
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "Calling RegularChain move assignment operator..." << std::endl;
	#endif
	RegularChain<Field,RecursivePoly>::operator=(std::move(a));
	regularChainOptions = a.regularChainOptions;
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A BPASTriangularSet
 **/
template <class Field, class RecursivePoly>
BPASTriangularSet<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (BPASTriangularSet<Field,RecursivePoly>&& a) {
	if (dynamic_cast<ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<ZeroDimensionalRegularChain<Field,RecursivePoly>&&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to ZeroDimensionalRegularChain."));
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A BPASRegularChain
 **/
template <class Field, class RecursivePoly>
BPASRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (BPASRegularChain<Field,RecursivePoly>&& a) {
	if (dynamic_cast<ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<ZeroDimensionalRegularChain<Field,RecursivePoly>&&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASRegularChain to ZeroDimensionalRegularChain."));
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A BPASZeroDimensionalRegularChain
 **/
template <class Field, class RecursivePoly>
BPASZeroDimensionalRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator= (BPASZeroDimensionalRegularChain<Field,RecursivePoly>&& a) {
	if (dynamic_cast<ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<ZeroDimensionalRegularChain<Field,RecursivePoly>&&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASZeroDimensionalRegularChain to ZeroDimensionalRegularChain."));
	return *this;
}
		
/**
 * Overload operator +
 * Adds a polynomial to a regular chain and returns a new regular chain
 *
 * @param p: A sparse multivariate polynomial
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly> ZeroDimensionalRegularChain<Field,RecursivePoly>::operator+ (const RecursivePoly& p) const {
	ZeroDimensionalRegularChain<Field,RecursivePoly> r(*this);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//		std::cerr << "trcVars in operator+:" << std::endl;
//		printVariables(this->trcVars);
//		std::cerr << "r trcVars in operator+:" << std::endl;
//		printVariables(r.transcendentalVariables());
	#endif
	return (r += p);
}

/**
 * Overload operator +=
 * Adds a polynomial to a regular chain
 *
 * @param p: A recursively viewed polynomial
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator+= (const RecursivePoly& p) {
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "entering ZDRC +=(p)" << std::endl;
		std::cerr << "p = " << p << std::endl;
	#endif
	
	// TODO: check that the condition on the variables of p relative
	//       to the algebraic and transcendental variables of the 
	//       current object is satisfied.
	
	if (!this->isRegular(p.initial())) {
		std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "in operator+=(p) of ZDRC" << std::endl;
		#endif
		exit(1);
	}

	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "calling constructChain(p) from +=(p)" << std::endl;
		std::cerr << "trcVars in operator+=:" << std::endl;
		printVariables(trcVars);
	#endif
	this->constructChain(p,ASSUME_REGULAR);
	return *this;
}
		
/**
 * Add operator +
 * Adds a polynomial to a regular chain and returns a new regular chain
 *
 * @param p: A sparse multivariate polynomial
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly> ZeroDimensionalRegularChain<Field,RecursivePoly>::operator+ (const RegularChain<Field,RecursivePoly>& T) const {
	RegularChain<Field,RecursivePoly> r(*this);
	return (r += T);
}

/**
 * Add assignment operator +=
 * Adds a polynomial to a regular chain
 *
 * @param p: A recursively viewed polynomial
 **/
template <class Field, class RecursivePoly>
ZeroDimensionalRegularChain<Field,RecursivePoly>& ZeroDimensionalRegularChain<Field,RecursivePoly>::operator+= (const RegularChain<Field,RecursivePoly>& T) {
	
	std::vector<RecursivePoly> polys(T.polynomials());
	
	for (int i=polys.size()-1; i>-1; --i) {
		if (!polys[i].isZero()) {
			*this += polys[i];
		}
	}
	// TODO: set the properties of *this according to the prior properties of *this and T
	// TODO: how to handle property conflict between *this and T? This would be a problem if *this requires SQUAREFREE but T does not, e.g.
	//       this is maybe not such an issue if this function is only called internally, but if it is called externally there could be issues.
	
	return *this;
}

/**
 * Overload comparison operator ==
 *
 *
 * @param a: A regular chain
 **/
template <class Field, class RecursivePoly>
bool ZeroDimensionalRegularChain<Field,RecursivePoly>::operator== (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) const {
	return RegularChain<Field,RecursivePoly>::operator==(a);
}

/**
 * Overload comparison operator !=
 *
 *
 * @param a: A regular chain
 **/
template <class Field, class RecursivePoly>
bool ZeroDimensionalRegularChain<Field,RecursivePoly>::operator!= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) const {
	return !(*this==a);
}

template <class Field, class RecursivePoly>	
void ZeroDimensionalRegularChain<Field,RecursivePoly>::lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const {
	RegularChain<Field,RecursivePoly> ts2;
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "moving input" << std::endl;
	#endif
	ts2 = std::move(ts);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "calling lower in RC" << std::endl;
		std::cerr << "current mode is: " << this->mode << std::endl;
	#endif
	RegularChain<Field,RecursivePoly>::lower(s,ts2);
	if (mode == TS_VARIABLE && !ts2.isEmpty()) {
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "ts2 = " << ts2 << std::endl;
			std::cerr << "ts2.isEmpty = " << ts2.isEmpty() << std::endl;
			std::cerr << "calling lowerSlice(s)" << std::endl;
		#endif
		ts2.lowerSlice(s);
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			ts2.display();
		#endif
	}
	ts = ZeroDimensionalRegularChain<Field,RecursivePoly>(std::move(ts2));
}

template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const {
	RegularChain<Field,RecursivePoly> ts2;
	ts2 = std::move(ts);
	RegularChain<Field,RecursivePoly>::upper(s,ts2);
	ts = std::move(ts2);
}

template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::cutChain(const ZeroDimensionalRegularChain<Field,RecursivePoly>& T, const Symbol& v, ZeroDimensionalRegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const {
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "entering cutChain(T,v,Tlv,Tv,Tgv):" << std::endl;
		std::cerr << "calling T.lower(v,Tlv):" << std::endl;
	#endif
	T.lower(v,Tlv);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "calling T.upper(v,Tgv):" << std::endl;
	#endif
	T.upper(v,Tgv);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		printVariables(Tgv.variables(),"Tgv");
		printVariables(Tlv.variables(),"Tlv");
		std::cerr << "calling T.select(v):" << std::endl;
		std::cerr << "this: " << *this << std::endl;
	#endif
	Tv = T.select(v);
}

template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::cutChain(const Symbol& v, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const {
		Tv = select(v);
		upper(v,Tgv);
}

template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::cutChain(const Symbol& v, ZeroDimensionalRegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv) const {
		lower(v,Tlv);
		Tv = select(v);
}
		
/**
 * Regularize the input polynomial with respect to the current regular chain,
 * i.e., compute a splitting of pairs (f_i,T_i) such that h_i is invertible
 * with respect to T_i.
 *
 * @param p: input polynomial
 * @param regularizeInitial: flag to call regularizeInitial or not (default true)
 **/
template <class Field, class RecursivePoly>
std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> ZeroDimensionalRegularChain<Field,RecursivePoly>::intersect(const RecursivePoly& p) const {

	typedef ZeroDimensionalRegularChain<Field,RecursivePoly> RC_INT_OBJ;

	std::vector<RC_INT_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_intersect), this, p);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntResult);
	RC_GEN_CONSUMER_LOOP(intResults, currIntResult, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(intResults, currIntResult, i);
		results.push_back(currIntResult);
	}

	return results;
}
		
/**
 * Compute the intersection of the varieties of the input polynomial and the 
 * current regular chain, i.e., compute a set of components (T_i) of T such 
 * that all zeros of T_i are common to p and all common zeros of p and T are 
 * included among the T_i.
 *
 * @param p: input polynomial
 **/
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::_intersect(const RecursivePoly& p, AsyncGenerator<ZeroDimensionalRegularChain<Field,RecursivePoly>>& results) const {
#else 
template <class Field, class RecursivePoly>
std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> ZeroDimensionalRegularChain<Field,RecursivePoly>::_intersect(const RecursivePoly& p) const {
#endif
	long long unsigned int rcProfilingStart;
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D int] entering intersect: " << p << std::endl;
	#endif
	
	typedef ZeroDimensionalRegularChain<Field,RecursivePoly> RC_INT_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_INT_OBJ, results);
	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegResult);
//	std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> results;
//	std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularComponents;
//	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> currRegularComponent;
	bool reduceMore(false);
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D int] computing squareFreePart:" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	RecursivePoly q(p.squareFreePart());
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
	#endif

	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D int] calling regularize: " << q << std::endl;
	#endif
	RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularize), this, q);
//	regularComponents = regularize(q);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D int] returning to intersect from regularize: " << q << std::endl;
	#endif
	RC_GEN_CONSUMER_LOOP(regResults, currRegResult, i) {
//	for (int i=0; i<regularComponents.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regResults, currRegResult, i);
//		currRegularComponent = regularComponents[i];
//		if (currRegularComponent.poly.isZero()) {
		if (currRegResult.poly.isZero()) {
			RC_GEN_PRODUCER_ACCUMULATE(results, currRegResult.chain);
//			RC_GEN_PRODUCER_ACCUMULATE(results, currRegularComponent.chain);
//			results.push_back(currRegularComponent.chain);
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D int] chain = " << currRegResult.chain << " added to results." << std::endl;
			#endif
		}
	}

//	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//		std::cerr << "[0D int] computing gcd-free (squarefree) factorization of p" << std::endl;
//	#endif
//	RecursivePoly q(p),c;
//	std::vector<RecursivePoly> lf,tasks;
//	lf = this->GCDFreeFactorization(q);
//	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//		std::cerr << "[0D int] GCD-free factors:" << std::endl;
//		for (int i=0; i<lf.size(); ++i) {
//			std::cerr << "lf[" << i << "] = " << lf[i] << std::endl;
//		}
//	#endif
//	for (int i=0; i<lf.size(); ++i) {
//		if (reduceMore) {
//			q = this->reduce(lf[i],c);
//			q *= c;
//		}
//		else {
//			q = lf[i];
//		}
//		if (!this->isConstantPolynomial(q)) {
//			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//				std::cerr << "[0D int] adding " << lf[i] << " to tasks." << std::endl;
//			#endif
//			tasks.push_back(q);
//		}
//	}
//	
//	for (int i=0; i<tasks.size(); ++i) {
//		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//			std::cerr << "[0D int] calling regularize: " << tasks[i] << std::endl;
//		#endif
//		regularComponents = regularize(tasks[i]);
//		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//			std::cerr << "[0D int] returning to intersect from regularize: " << tasks[i] << std::endl;
//		#endif
//		for (int j=0; j<regularComponents.size(); ++j) {
//			currRegularComponent = regularComponents[j];
//			if (currRegularComponent.poly.isZero()) {
//				results.push_back(currRegularComponent.chain);
//				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//					std::cerr << "[0D int] chain = " << currRegularComponent.chain << " added to results." << std::endl;
//				#endif
//			}
//		}
//	}
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//		std::cerr << "[0D int] results.size = " << results.size() << std::endl;
		std::cerr << "[0D int] leaving intersect: " << p << std::endl;
	#endif
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}
		
/**
 * Regularize the input polynomial with respect to the current regular chain,
 * i.e., compute a splitting of pairs (f_i,T_i) such that h_i is invertible
 * with respect to T_i.
 *
 * @param p: input polynomial
 * @param regularizeInitial: flag to call regularizeInitial or not (default true)
 **/
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::regularize(const RecursivePoly& f) const {

	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REG_OBJ;

	std::vector<RC_REG_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularize), this, f);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegResult);
	RC_GEN_CONSUMER_LOOP(regResults, currRegResult, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regResults, currRegResult, i);
		results.push_back(currRegResult);
	}

	return results;
}
		
/**
 * Regularize the input polynomial with respect to the current regular chain,
 * i.e., compute a splitting of pairs (f_i,T_i) such that h_i is invertible
 * with respect to T_i.
 *
 * @param p: input polynomial
 * @param regularizeInitial: flag to call regularizeInitial or not (default true)
 **/
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::_regularize(const RecursivePoly& f, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const {
#else 
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::_regularize(const RecursivePoly& f) const {
#endif
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D reg] entering regularize: " << f << std::endl;
	#endif
	RecursivePoly c;
	
	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REG_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REG_OBJ, results);
//	std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> results;
	// std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> isInv;
	// BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> currIsInv;
	RecursivePoly h(this->reduceMinimal(f));
//	RecursivePoly h(this->reduce(f));
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D reg]\treduceMinimal(f) = " << h << std::endl;
	#endif
	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> t(h,*this);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D reg] regularizing with respect to input polynomial..." << std::endl;
	#endif
	if (this->isConstantPolynomial(h) || this->isEmpty()) {
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D reg] polynomial is constant or chain is empty...regular result obtained." << std::endl;
		#endif
		RC_GEN_PRODUCER_ACCUMULATE(results, t);
		RC_GEN_PRODUCER_COMPLETE(results);
//		results.push_back(t);
//		return results;
	}
	RecursivePoly zero;

	typedef BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_INV_OBJ;

	RC_GEN_CONSUMER_INIT(RC_INV_OBJ, isInv, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_isInvertible), this, h);
	// isInv = this->isInvertible(h);

	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INV_OBJ, currIsInv);
	RC_GEN_CONSUMER_LOOP(isInv, currIsInv, i) {
	// for (int i=0; i<isInv.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(isInv, currIsInv, i);
		// currIsInv = std::move(isInv[i]);

		if (currIsInv.isTrue) {
			t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currIsInv.chain.reduceMinimal(h),std::move(currIsInv.chain));
//			t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currIsInv.chain.reduce(h),std::move(currIsInv.chain));
			RC_GEN_PRODUCER_ACCUMULATE(results, t);
//			results.push_back(t);
		}
		else {
			t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(zero,std::move(currIsInv.chain));
			RC_GEN_PRODUCER_ACCUMULATE(results, t);
//			results.push_back(t);
		}
	}
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D reg] leaving regularize: " << f << std::endl;
	#endif
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}
		
/**
 * Regularize the initial of the input polynomial with respect to the current regular chain,
 * i.e., compute a splitting of pairs (f_i,T_i) such that h_i has invertible
 * initial with respect to T_i.
 *
 * @param p: input polynomial
 **/
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::regularizeInitial(const RecursivePoly& f) const {

	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REGINIT_OBJ;

	std::vector<RC_REGINIT_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_REGINIT_OBJ, regularComponents, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularizeInitial), this, f);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGINIT_OBJ, currRegularComponent);
	// regularComponents = regularizeInitial(f);
	RC_GEN_CONSUMER_LOOP(regularComponents, currRegularComponent, i) {
	// for(int i = 0; i < regularComponents.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regularComponents, currRegularComponent, i);
		results.push_back(currRegularComponent);
	}

	return results;
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularizeInitial(const RecursivePoly& f, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::_regularizeInitial(const RecursivePoly& f) const {
#endif
	++ZDregularizeInitialDepth;
	++ZDdepth;
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] entering regularize: " << f << std::endl;
	#endif
	RecursivePoly c;

	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REGINIT_OBJ;
	// std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> results;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REGINIT_OBJ, results);

	typedef BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_INV_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INV_OBJ, currIsInv);
	// std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> isInv;
	// BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> currIsInv;

	RecursivePoly h(this->reduceMinimal(f));
//	RecursivePoly h(this->reduce(f));
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "]\treduceMinimal(f) = " << h << std::endl;
	#endif
	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> t(std::move(h),*this);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] regularizing with respect to initial of input polynomial..." << std::endl;
	#endif
	std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> tasks;
	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> currTask;
	tasks.push_back(t);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] starting while loop in regularize..." << std::endl;
	#endif
	while (!tasks.empty()) {
		currTask = tasks.back();
		tasks.pop_back();
		if (currTask.chain.isConstantPolynomial(currTask.poly)) {
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] polynomial is constant! regular result obtained." << std::endl;
			#endif
			RC_GEN_PRODUCER_ACCUMULATE(results, currTask);
			// results.push_back(currTask);
		}
		else {
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] polynomial is not constant, continuing in regularize..." << std::endl;
			#endif
			c = currTask.poly.initial();
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] currTask.poly = " << currTask.poly << std::endl;
				std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] currTask.poly.initial = " << currTask.poly.initial() << std::endl;
			#endif
			if (!currTask.chain.isConstantPolynomial(c)) {

				RC_GEN_CONSUMER_INIT(RC_INV_OBJ, isInv, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_isInvertible), &(currTask.chain), c);
				// isInv = currTask.chain.isInvertible(c);
				RC_GEN_CONSUMER_LOOP(isInv, currIsInv, i) {
				// for (int i=0; i<isInv.size(); ++i) {
					RC_GEN_CONSUMER_GET_LOOPELEM(isInv, currIsInv, i);
					// currIsInv = isInv[i];
					if (currIsInv.isTrue) {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] invertible polynomial found! regular result obtained." << std::endl;
						#endif
						t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currIsInv.chain.reduceMinimal(currTask.poly),std::move(currIsInv.chain));
//						t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currIsInv.chain.reduce(currTask.poly),std::move(currIsInv.chain));
						RC_GEN_PRODUCER_ACCUMULATE(results, t);
						// results.push_back(t);
					}
					else {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] invertible polynomial not found, adding new task..." << std::endl;
						#endif
						h = RecursivePoly(currIsInv.chain.reduceMinimal(currTask.poly.tail()));
//						h = RecursivePoly(currIsInv.chain.reduce(currTask.poly.tail()));
						t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(std::move(h),std::move(currIsInv.chain));
						tasks.push_back(t);
					}
				}
			}
			else {
				if (c.isZero()) {
					// currIsInv = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(false,std::move(currTask.chain));
					// isInv.clear();
					// isInv.push_back(currIsInv);
					t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currTask.chain.reduceMinimal(currTask.poly.tail()),std::move(currTask.chain));
//					t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currTask.chain.reduce(currTask.poly.tail()),std::move(currTask.chain));
					tasks.push_back(t);
				}
				else {
					// currIsInv = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(true,std::move(currTask.chain));
					// isInv.clear();
					// isInv.push_back(currIsInv);
					t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currTask.chain.reduceMinimal(currTask.poly),std::move(currTask.chain));
//					t = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(currTask.chain.reduce(currTask.poly),std::move(currTask.chain));
					RC_GEN_PRODUCER_ACCUMULATE(results, t);
				}
			}
//			isInv = currTask.chain.isInvertible(currTask.poly.initial());
//			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//				std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] isInv.size = " << isInv.size() << std::endl;
//			#endif
		}
	}
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regInit][" << ZDdepth << "][" << ZDregularizeInitialDepth << "] leaving regularize: " << f << std::endl;
	#endif
	--ZDregularizeInitialDepth;
	--ZDdepth;
	RC_GEN_PRODUCER_COMPLETE(results);
}


/**
 * Determine whether a polynomial is invertible with respect to the current regular chain,
 * i.e. compute a splitting of pairs (b_i,T_i) such that b_i is true if f is invertible
 * modulo T_i and b_i is false if f is zero modulo T_i.
 *
 * @param p: input polynomial
 **/
template <class Field, class RecursivePoly>
std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::isInvertible(const RecursivePoly& f) const {

	typedef BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_INV_OBJ;

	std::vector<RC_INV_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_INV_OBJ, isInv, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_isInvertible), this, f);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INV_OBJ, currIsInv);
	RC_GEN_CONSUMER_LOOP(isInv, currIsInv, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(isInv, currIsInv, i);
		results.push_back(currIsInv);
	}

	return results;
}


#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::_isInvertible(const RecursivePoly& f, AsyncGenerator<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::_isInvertible(const RecursivePoly& f) const {
#endif
	long long unsigned int rcProfilingStart;
	++ZDisInvertibleDepth;
	++ZDdepth;
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] entering isInvertible: " << f << std::endl;
	#endif
	Symbol v;
	RecursivePoly h,g,q,Tv,temp;
	std::vector<RecursivePoly> polys;
	ZeroDimensionalRegularChain<Field,RecursivePoly> T,L;
	RegularChain<Field,RecursivePoly> U;
	SubResultantChain<RecursivePoly,RecursivePoly> src;
	// std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> results, moreResults;
	
	std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> constructedChains, moreConstructedChains;

	typedef BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_INV_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_INV_OBJ, results);

	// std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> subResultantComponents;
	BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> currResult;
	// PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> currRegularComponent,currSubResultantComponent;

	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REGINIT_OBJ;
	typedef RC_REGINIT_OBJ RC_LNZSR_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_LNZSR_OBJ, currSubResultantComponent);
	
	// std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularComponents;
	RC_GEN_CONSUMER_INIT(RC_REGINIT_OBJ, regularComponents, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularizeInitial), this, f);
	// regularComponents = regularizeInitial(f);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] returning to isInvertible..." << std::endl;
//		std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] regularComponents.size() = " << regularComponents.size() << std::endl;
	#endif
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGINIT_OBJ, currRegularComponent);
	RC_GEN_CONSUMER_LOOP(regularComponents, currRegularComponent, i) {
	// for (int i=0; i<regularComponents.size(); ++i) {
//		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//			std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] i = " << i << std::endl;
//		#endif
		RC_GEN_CONSUMER_GET_LOOPELEM(regularComponents, currRegularComponent, i);
		// currRegularComponent = regularComponents[i];
		h = std::move(currRegularComponent.poly);
		h = h.mainPrimitivePart();
		T = std::move(currRegularComponent.chain);

		// std::cerr << "[0D isInv] h : " << h << std::endl;
		// std::cerr << "[0D isInv] prime: " << T.isSaturatedIdealPrime() << std::endl;
		// std::cerr << "[0D isInv] rc: " << T << std::endl;

		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] h = " << h << std::endl;
			std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] isConstantPolynomial(h) = " << T.isConstantPolynomial(h) << std::endl;
			std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] T = " << T << std::endl;
		#endif
		if (h.isZero()) {
//		if (h.isZero() || T.isEmpty()) {
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] false on T" << std::endl;
			#endif
			currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(false,T);
			RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
			// results.push_back(currResult);
		}
		else if (T.isConstantPolynomial(h)) {
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] true on T" << std::endl;
			#endif
			currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(true,T);
			RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
			// results.push_back(currResult);
		} else if (T.isSaturatedIdealPrime()) {
			if (this->isInSaturatedIdealMinimal(h)) {
				currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(false,T);
			} else {
				currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(true,T);
			}
			RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
		}
		else {
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] processing nonconstant polynomial" << std::endl;
			#endif
			v = h.leadingVariable();
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] v = " << v << std::endl;
			#endif
			cutChain(T,v,L,Tv,U);

			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] T_<v = " << L << std::endl;
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] T_v = " << Tv << std::endl;
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] T_>v = " << U << std::endl;
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] h = " << h << std::endl;
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] computing subresultant chain of T_v and h:" << std::endl;
			#endif
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			src = SubResultantChain<RecursivePoly,RecursivePoly>(Tv,h,v);
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&subresultantChainTime);
			#endif
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] src = " << src << std::endl;
			#endif
			bool assumeReg = true;
			RC_GEN_CONSUMER_INIT(RC_LNZSR_OBJ, subResultantComponents, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::lastNonZeroSubResultant_inner), &(L), Tv, h, v, src, assumeReg);
			// subResultantComponents = L.lastNonZeroSubResultant_inner(Tv,h,v,src,assumeReg);
//			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//				std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] subResultantComponents.size() = " << subResultantComponents.size() << std::endl;
//			#endif

			int constructOptions = ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE | ASSUME_ZERODIMENSIONAL;
			if (this->isSquareFree()) {
				 // since Tv was square free, g must be also.
				constructOptions |= ASSUME_SQUAREFREE;
			}

			RC_GEN_CONSUMER_LOOP(subResultantComponents, currSubResultantComponent, j) {
			// for (int j=0; j<subResultantComponents.size(); ++j) {
//				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//					std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "]\tj = " << j << std::endl;
//				#endif
				RC_GEN_CONSUMER_GET_LOOPELEM(subResultantComponents, currSubResultantComponent, j);
				// currSubResultantComponent = subResultantComponents[j];
				g = std::move(currSubResultantComponent.poly);
				L = std::move(currSubResultantComponent.chain);
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "]\tg = " << g << std::endl;
					std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "]\tL = " << L << std::endl;
				#endif
				if (g.degree(v) <= 0) {
					//deg(g, v) = 0 or g = 0	
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] T variables before merge:" << std::endl;
						printVariables(T.variables());
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] true on U+Tv" << std::endl;
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] Tv = " << Tv << std::endl;
					#endif
					T = std::move(L);
//					T.constructChain(Tv,ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE);

					constructedChains = T.constructChains(Tv, constructOptions);

					for (int i=0; i<constructedChains.size(); ++i) {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] L+Tv = " << constructedChains[i] << std::endl;
						#endif
//						constructedChains[i].constructChain(U,ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE);

						moreConstructedChains = constructedChains[i].constructChains(U, constructOptions);

						for (int k=0; k<moreConstructedChains.size(); ++k) {
							#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
								std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] L+Tv+U = " << moreConstructedChains[k] << std::endl;
							#endif
							currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(true,moreConstructedChains[k]);
							RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
							// results.push_back(currResult);
						}
					}
				} else if (g.degree(v) == Tv.degree(v)) {
					T = std::move(L);
					constructedChains = T.constructChains(Tv, constructOptions);
					for (int i = 0; i < constructedChains.size(); ++i) {
						moreConstructedChains = constructedChains[i].constructChains(U, constructOptions);
						for (int k = 0; k < moreConstructedChains.size(); ++k) {
							currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(false, moreConstructedChains[k]);
						}
					}
				} else {
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] before pseudoDivision:" << std::endl;
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] Tv = " << Tv << std::endl;
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] g: " << g << std::endl;
						printVariables(g.ringVariables(), "gRingVars");
						printVariables(Tv.ringVariables(), "TvRingVars");
					#endif
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&rcProfilingStart);
					#endif
					Tv.pseudoDivide(g,&q,NULL,1);
					g.squareFreePart();
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
					#endif
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] pquo(Tv,g,v) = " << q << std::endl;
						std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] false on T+g" << std::endl;
					#endif
					T = L;
//					T.constructChain(g,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE);
					
					constructedChains = T.constructChains(g,constructOptions);

					for (int i=0; i<constructedChains.size(); ++i) {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] T+g = " << constructedChains[i] << std::endl;
						#endif
//						constructedChains[i].constructChain(U,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE);

						moreConstructedChains = constructedChains[i].constructChains(U,constructOptions);

						for (int k=0; k<moreConstructedChains.size(); ++k) {
							#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
								std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] L+g+U = " << moreConstructedChains[k] << std::endl;
							#endif
							currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(false,moreConstructedChains[k]);
							RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
							// results.push_back(currResult);
						}
					}

					//by check above for deg(g) == deg(Tv), this should never fail
					if (q.degree(v) > 0) {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] true on T+q" << std::endl;
						#endif
						T = std::move(L);
	//					T.constructChain(q.mainPrimitivePart(),ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE);

						constructedChains = T.constructChains(q.mainPrimitivePart(),ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE | ASSUME_ZERODIMENSIONAL);
						
						for (int i=0; i<constructedChains.size(); ++i) {
							#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
								std::cerr << "[0D isInv] L+q = " << constructedChains[i] << std::endl;
							#endif
	//						constructedChains[i].constructChain(U,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE);

							moreConstructedChains = constructedChains[i].constructChains(U,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE | ASSUME_PRIMITIVE | ASSUME_ZERODIMENSIONAL);
							
							for (int k=0; k<moreConstructedChains.size(); ++k) {
								#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
									std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] L+q+U = " << moreConstructedChains[k] << std::endl;
								#endif
								if (this->isSquareFree()) {
									currResult = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(true,moreConstructedChains[k]);
									RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
									// results.push_back(currResult);
								}
								else {
									RC_GEN_CONSUMER_INIT(RC_INV_OBJ, moreResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_isInvertible), &(moreConstructedChains[k]), moreConstructedChains[k].reduceMinimal(f));
	//								RC_GEN_CONSUMER_INIT(RC_INV_OBJ, moreResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_isInvertible), &(moreConstructedChains[j]), moreConstructedChains[j].reduce(f));
									RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INV_OBJ, moreCurrResults);
									RC_GEN_CONSUMER_LOOP(moreResults, moreCurrResults, kk) {
										RC_GEN_CONSUMER_GET_LOOPELEM(moreResults, moreCurrResults, kk);
										RC_GEN_PRODUCER_ACCUMULATE(results, moreCurrResults);
									}
									// moreResults = T.isInvertible(T.reduce(f));
									// results.insert(results.end(),moreResults.begin(),moreResults.end());
								}
							}
						}
					}
				}
			}
		}
	}
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D isInv][" << ZDdepth << "][" << ZDisInvertibleDepth << "] leaving isInvertible: " << f << std::endl;
	#endif
	--ZDisInvertibleDepth;
	--ZDdepth;
	RC_GEN_PRODUCER_COMPLETE(results);
	// return results;
}
		
/**
 * Compute the last nonzero subresultant of the subresultant chain src of f and g modulo the current regular chain,
 * i.e., compute a splitting of pairs (R_i,T_i) such that R_i is the last nonzero subresultant modulo T_i.
 *
 * @param f: input polynomial
 * @param g: input polynomial
 * @param v: common main variable of f and g
 * @param src: subresultant chain of f and g
 * @param assumeRegular: boolean flag specifying whether the function is called for GCD computation
 **/
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::lastNonZeroSubResultant(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, bool assumeRegular, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const {
#else 
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::lastNonZeroSubResultant(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, bool assumeRegular) const {
#endif
	long long unsigned int rcProfilingStart;
	++ZDlnzsrDepth;
	++ZDdepth;
	std::vector<Symbol> vs,vs2,vs3;
//	vs.push_back(v);
//	std::vector<Symbol> fgvs,Tvs,vs;
//	fgvs = setUnion(f.variables(),g.variables());
//	Tvs = setUnion(this->vars,vs);
	
	vs = this->vars;
	if (contains(vs,v)) {
		std::cerr << "BPAS: error, ZeroDimensionalRegularChain cannot contain the variable input to lastNonZeroSubResultant." << std::endl;
		exit(1);
	}
	vs.insert(vs.begin(),v);
//	if (!isSubset(fgvs,Tvs)) {
	if (!isSubset(setUnion(f.variables(),g.variables()),vs)) {
		std::cerr << "BPAS: error, polynomials input to lastNonZeroSubResultant have variables other than v and the those of the ZeroDimensionalRegularChain." << std::endl;
		exit(1);
	}
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] entering lastNonZeroSubResultant: " << std::endl;
//		std::cerr << "vars+v: " << std::endl;
//		printVariables(vs);
	#endif
	RecursivePoly F(f),G(g),one;
	vs2 = orderPreservingSetIntersection(vs,f.ringVariables());
	vs3 = orderPreservingSetDifference(f.ringVariables(),vs2);
	vs2.reserve(vs2.size()+vs3.size());
	vs2.insert(vs2.begin(),vs3.begin(),vs3.end());
	F.setRingVariables(vs2);
	vs2 = orderPreservingSetIntersection(vs,g.ringVariables());
	vs3 = orderPreservingSetDifference(g.ringVariables(),vs2);
	vs2.reserve(vs2.size()+vs3.size());
	vs2.insert(vs2.begin(),vs3.begin(),vs3.end());
	G.setRingVariables(vs2);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//		std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] f vars reordered: " << std::endl;
//		printVariables(F.ringVariables());
//		std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] g vars reordered: " << std::endl;
//		printVariables(G.ringVariables());
		std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] F.initial = " << F.initial() << std::endl;
	#endif
	one.one();
	
	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_LNZSR_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_LNZSR_OBJ, results);
//	std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> results;//,moreResults,regularF,regularG;
	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> currResult;//,currRegularF,currRegularG;

	SubResultantChain<RecursivePoly,RecursivePoly> src;
	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_LNZSR_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_LNZSR_OBJ, currLNZSRIn);	

	typedef RC_LNZSR_OBJ RC_REGINIT_OBJ;
	RC_GEN_CONSUMER_INIT(RC_REGINIT_OBJ, regularF, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularizeInitial), this, F);
	// regularF = regularizeInitial(F);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
	std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] returning to lastNonZeroSubResultant..." << std::endl;
	#endif
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGINIT_OBJ,currRegularF);
	RC_GEN_CONSUMER_LOOP(regularF, currRegularF, i) {
	// for (int i=0; i<regularF.size(); ++i) {
//		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//			std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] i = " << i << std::endl;
//		#endif
		RC_GEN_CONSUMER_GET_LOOPELEM(regularF, currRegularF, i);
		// currRegularF = std::move(regularF[i]);
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] after regularizing f:" << std::endl;
			std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] F = " << currRegularF.poly << " modulo " << currRegularF.chain << std::endl;
		#endif

		//taking pointer to currRegularF is safe because (even with generators) since this inner loop finishes
		//before we get a value of currRegularF in the outer loop.
		RC_GEN_CONSUMER_INIT(RC_REGINIT_OBJ, regularG, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularizeInitial), &(currRegularF.chain), G);
		// regularG = currRegularF.chain.regularizeInitial(G);
		RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGINIT_OBJ, currRegularG);
		RC_GEN_CONSUMER_LOOP(regularG, currRegularG, j) {
		// for (int j=0; j<regularG.size(); ++j) {
//			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//				std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] i = " << i << std::endl;
//			#endif
			RC_GEN_CONSUMER_GET_LOOPELEM(regularG, currRegularG, j);
			// currRegularG = std::move(regularG[j]);
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			currRegularF.poly = currRegularF.poly.mainPrimitivePart();
			currRegularG.poly = currRegularG.poly.mainPrimitivePart();
//			currRegularF.poly = currRegularF.poly.primitivePart(currRegularF.poly.leadingVariable());
//			currRegularG.poly = currRegularG.poly.primitivePart(currRegularG.poly.leadingVariable());
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
			#endif
			#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
				std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] after regularizing g:" << std::endl;
				std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] G = " << currRegularG.poly << " modulo " << currRegularG.chain << std::endl;
			#endif
			if (currRegularF.poly.isZero()) {
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] f is zero!" << std::endl;
				#endif
				RC_GEN_PRODUCER_ACCUMULATE(results, currRegularG);
//				results.push_back(currRegularG);
			}
			else if (currRegularG.poly.isZero()) {
				#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
					std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] g is zero" << std::endl;
				#endif
				currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(std::move(currRegularF.poly),std::move(currRegularG.chain));
				RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
//				results.push_back(currResult);
			}
			else if ((currRegularF.poly.isConstant() == 0) && currRegularF.poly.leadingVariable() == v) {
				if ((currRegularG.poly.isConstant() == 0) && currRegularG.poly.leadingVariable() == v) {
//					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//						std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] f and g both have v as leading variable! (i = " << i << ")" << std::endl;
//					#endif
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&rcProfilingStart);
					#endif
					src = SubResultantChain<RecursivePoly,RecursivePoly>(currRegularF.poly,currRegularG.poly,v);
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&rcProfilingStart,&subresultantChainTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_LNZSR_OBJ, subResultantComponents, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::lastNonZeroSubResultant_inner), &(currRegularG.chain), currRegularF.poly, currRegularG.poly, v, src, assumeRegular);
					// moreResults = currRegularG.chain.lastNonZeroSubResultant_inner(currRegularF.poly,currRegularG.poly,v,src,assumeRegular);
					RC_GEN_CONSUMER_LOOP(subResultantComponents, currLNZSRIn, k) {
						RC_GEN_CONSUMER_GET_LOOPELEM(subResultantComponents, currLNZSRIn, k);
						RC_GEN_PRODUCER_ACCUMULATE(results, currLNZSRIn);
//						results.push_back(currLNZSRIn);
					}
					// results.reserve(results.size()+moreResults.size());
					// results.insert(results.end(),moreResults.begin(),moreResults.end());
				}
				else {
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] g is the lnzsr!" << std::endl;
					#endif
					RC_GEN_PRODUCER_ACCUMULATE(results, currRegularG);
//					results.push_back(currRegularG);
				}
			}
			else {
				if ((currRegularG.poly.isConstant() == 0) && currRegularG.poly.leadingVariable() == v) {
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] f is the lnzsr!" << std::endl;
					#endif
					currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(std::move(currRegularF.poly),std::move(currRegularG.chain));
					RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
//					results.push_back(currResult);
				}
				else {
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] lnzsr is 1!" << std::endl;
					#endif
					currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(one,std::move(currRegularG.chain));
					RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
//					results.push_back(currResult);
				}
			}
		}
	}
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D lnzsr][" << ZDdepth << "][" << ZDlnzsrDepth << "] leaving lastNonZeroSubResultant: " << std::endl;
	#endif
	--ZDlnzsrDepth;
	--ZDdepth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}
		
/**
 * Compute the last nonzero subresultant of the subresultant chain src of f and g modulo the current regular chain,
 * i.e., compute a splitting of pairs (R_i,T_i) such that R_i is the last nonzero subresultant modulo T_i.
 *
 * @param f: input polynomial
 * @param g: input polynomial
 * @param v: common main variable of f and g
 * @param src: subresultant chain of f and g
 * @param assumeRegular: boolean flag specifying whether the function is called for GCD computation
 **/

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::lastNonZeroSubResultant_inner(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool assumeRegular, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::lastNonZeroSubResultant_inner(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool assumeRegular) const {
#endif
	long long unsigned int rcProfilingStart;
	if (v != f.leadingVariable() || v != g.leadingVariable()) {
		std::cerr << "BPAS: error, leading variable of f and g must be v in input to lastNonZeroSubResultant." << std::endl;
		std::cerr << "f: " << f << std::endl << "g: " << g << std::endl;
		exit(1);
	}
	++ZDlnzsrInnerDepth;
	++ZDdepth;
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] entering lastNonZeroSubResultant_inner: " << f << std::endl;
	#endif
	int k(0);
	RecursivePoly s,c;
	RecursivePoly uP;
//	SparseUnivariatePolynomial<RecursivePoly> uP;
	struct Task {
		RecursivePoly a;
		RecursivePoly b;
		ZeroDimensionalRegularChain<Field,RecursivePoly> T;
		int k;
	};

	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_LNZSR_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_LNZSR_OBJ, results);
	// std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> results;
	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> currResult;

	typedef BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_INV_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INV_OBJ, currIsInv);
	// std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> isInv;
	// BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>> currIsInv;

//	while (src.subResultantOfIndex(k).isZero()) {
//		++k;
//	}
	std::vector<Task> tasks;
	Task t = {f,g,*this,k},currTask;
	tasks.push_back(t);
	while (!tasks.empty()) {
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] lnzsr tasks.size = " << tasks.size() << std::endl;
		#endif
		currTask = tasks.back();
		tasks.pop_back();
		k = currTask.k;
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] k = " << k << std::endl;
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] v = " << v << std::endl;
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] T = " << currTask.T << std::endl;
		#endif
		uP = src.subResultantOfIndex(k);
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] src[" << k << "] = " << uP << std::endl;
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] uP: " << uP << std::endl;
		#endif

		RecursivePoly supLead = uP.leadingCoefficient();
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] uP.leadingCoefficient: " << supLead << std::endl;
		#endif

		// TODO: remove either s or uP as both are not needed
		s = uP;
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] s = " << s << std::endl;
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] s.leadingVariable = " << s.leadingVariable() << std::endl;
		#endif

		c = s.leadingCoefficientInVariable(v);
		#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
			std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] c = " << c << std::endl;
		#endif
		if (c.isConstant() == 0) {
			// isInv = currTask.T.isInvertible(c);
			RC_GEN_CONSUMER_INIT(RC_INV_OBJ, isInv, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_isInvertible), &(currTask.T), c);
			RC_GEN_CONSUMER_LOOP(isInv, currIsInv, i) {
			// for (int i=0; i<isInv.size(); ++i) {
				RC_GEN_CONSUMER_GET_LOOPELEM(isInv, currIsInv, i);
				// currIsInv = std::move(isInv[i]);
				if (currIsInv.isTrue) {
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] lnzsr found!" << std::endl;
						std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] T = " << currIsInv.chain << std::endl;
						std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] S before reduceMinimal: " << s << std::endl;
						std::cerr << "s = " << currIsInv.chain.reduceMinimal(s) << std::endl;
//						std::cerr << "s = " << currIsInv.chain.reduce(s) << std::endl;
					#endif
					RecursivePoly temp(currIsInv.chain.reduceMinimal(s));
//					RecursivePoly temp(currIsInv.chain.reduce(s));
					#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
						std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] \t\t\treduceMinimal(s) = " << temp << std::endl;
					#endif
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&rcProfilingStart);
					#endif
					currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(temp.mainPrimitivePart(),std::move(currIsInv.chain));
	//				currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(temp.primitivePart(temp.leadingVariable()),currIsInv.chain);
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
					#endif
					RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
					// results.push_back(currResult);
				}
				else { // c equiv 0 mod(T)
					if (currTask.k == currTask.b.leadingVariableDegree().get_ui()) {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] lnzsr is second last subresultant" << std::endl;
						#endif
						RecursivePoly temp(currIsInv.chain.reduceMinimal(currTask.b));
//						RecursivePoly temp(currIsInv.chain.reduce(currTask.b));
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] \t\t\treduceMinimal(b) = " << temp << std::endl;
						#endif
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&rcProfilingStart);
						#endif
						currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(temp.mainPrimitivePart(),std::move(currIsInv.chain));
	//					currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(temp.primitivePart(temp.leadingVariable()),currIsInv.chain);
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
						#endif
						RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
						// results.push_back(currResult);
					}
					else {
						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
							std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] lnzsr not found, adding new task..." << std::endl;
						#endif
//						RecursivePoly temp(currIsInv.chain.reduce(currTask.a));
//						RecursivePoly temp2(currIsInv.chain.reduce(currTask.b));
//						#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
//							std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] \treduce(a) = " << temp << std::endl;
//							std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] \treduce(b) = " << temp2 << std::endl;
//						#endif
						t = {std::move(currTask.a),std::move(currTask.b),std::move(currIsInv.chain),currTask.k+1};
//						t = {temp,temp2,std::move(currIsInv.chain),currTask.k+1};
						tasks.push_back(t);
					}
				}
			}
		}
		else {
			if (c.isZero()) {
				// currIsInv = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(false,std::move(currTask.T));
				// isInv.clear();
				// isInv.push_back(currIsInv);
				if (currTask.k == currTask.b.leadingVariableDegree().get_ui()) {
					RecursivePoly temp(currTask.T.reduceMinimal(currTask.b));
//					RecursivePoly temp(currTask.T.reduce(currTask.b));
					currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(temp.mainPrimitivePart(),std::move(currTask.T));
					// RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
					RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
				}
				else {
//					RecursivePoly temp(currTask.T.reduce(currTask.a));
//					RecursivePoly temp2(currTask.T.reduce(currTask.b));
					t = {std::move(currTask.a),std::move(currTask.b),std::move(currTask.T),currTask.k+1};
//					t = {temp,temp2,std::move(currTask.T),currTask.k+1};
					tasks.push_back(t);
				}
			}
			else {
				// I think this reduce call is superfluous
				RecursivePoly temp(currTask.T.reduceMinimal(s));
//				RecursivePoly temp(currTask.T.reduce(s));
				currResult = PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>(temp.mainPrimitivePart(),std::move(currTask.T));
				RC_GEN_PRODUCER_ACCUMULATE(results, currResult);
				// results.push_back(currResult);
				// currIsInv = BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>(true,std::move(currTask.T));
				// isInv.clear();
				// isInv.push_back(currIsInv);
			}
		}
//		isInv = currTask.T.isInvertible(c);
	}
	
//	RecursivePoly h(f);
//	ZeroDimensionalRegularChain<Field,RecursivePoly> rc(*this);
//	PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> pcp (h,rc);
//	results.push_back(pcp);
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D lnzsrIn][" << ZDdepth << "][" << ZDlnzsrInnerDepth << "] leaving lastNonZeroSubResultant_inner: " << f << std::endl;
	#endif
	--ZDlnzsrInnerDepth;
	--ZDdepth;
	RC_GEN_PRODUCER_COMPLETE(results);
	// return results;		
}
		
/**
 * Compute the gcd of polynomials f and g with respect to v modulo the current object
 * regular chain
 *
 * @param f: input polynomial
 * @param g: input polynomial
 * @param v: common main variable of f and g
 **/
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::regularGCD(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v) {

	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;

	std::vector<RC_REGGCD_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_REGGCD_OBJ, regGCDResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularGCD), this, f, g, v);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGGCD_OBJ, currRegGCDResult);
	RC_GEN_CONSUMER_LOOP(regGCDResults, currRegGCDResult, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regGCDResults, currRegGCDResult, i);
		results.push_back(currRegGCDResult);
	}

	return results;
}
		
/**
 * Compute the gcd of polynomials f and g with respect to v modulo the current object
 * regular chain
 *
 * @param f: input polynomial
 * @param g: input polynomial
 * @param v: common main variable of f and g
 **/
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::_regularGCD(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) {
#else 
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> ZeroDimensionalRegularChain<Field,RecursivePoly>::_regularGCD(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v) {
#endif
	if (v != f.leadingVariable() || v != g.leadingVariable()) {
		std::cerr << "BPAS: error, leading variable of f and g must be v in input to regularGCD." << std::endl;
		exit(1);
	}
	++ZDregularGCDDepth;
	++ZDdepth;
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regGCD][" << ZDdepth << "][" << ZDregularGCDDepth << "] entering regularGCD: " << std::endl;
	#endif
	
	typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REGGCD_OBJ, results);
	typedef RC_REGGCD_OBJ RC_LNZSR_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_LNZSR_OBJ, currLNZSRResult);
//	std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> result;

	bool assumeReg = true;
	RC_GEN_CONSUMER_INIT(RC_LNZSR_OBJ, LNZSRResults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::lastNonZeroSubResultant), this, f, g, v, assumeReg);
//	result = lastNonZeroSubResultant(f,g,v,true);
	
	RC_GEN_CONSUMER_LOOP(LNZSRResults, currLNZSRResult, i) {
//	for (auto res : result) {
//		RC_GEN_PRODUCER_ACCUMULATE(results, res);
		RC_GEN_CONSUMER_GET_LOOPELEM(LNZSRResults, currLNZSRResult, i);
		RC_GEN_PRODUCER_ACCUMULATE(results, currLNZSRResult);
	}
	
	#ifdef ZERODIMENSIONALREGULARCHAIN_DEBUG
		std::cerr << "[0D regGCD][" << ZDdepth << "][" << ZDregularGCDDepth << "] leaving regularGCD: " << std::endl;
	#endif
	--ZDregularGCDDepth;
	--ZDdepth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return result;
}

/**
 * Generate a random zero dimensional regular chain
 *
 * @param nVars: number of variables = number of algebraic variables
 * @param nTrcVars: number of transcendental variables
 * @param nTerms: maximum number of terms in the polynomials
 * @param coefBound: maximum coefficient size
 * @param pSparsity: sparsity of the polynomials
 * @param includeNeg: whether to include negative coefficients
 **/
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::randomZeroDimensionalRegularChain(int nVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg) {
	this->RegularChain<Field,RecursivePoly>::randomRegularChain(nVars,nVars,nTrcVars,nTerms,coefBound,pSparsity,includeNeg);
}

/**
 * Generate a random zero dimensional regular chain
 *
 * @param nVars: number of variables = number of algebraic variables
 * @param nTrcVars: number of transcendental variables
 * @param maxDegs: maximum degrees among the full list of variables
 * @param coefBound: maximum coefficient size
 * @param pSparsity: sparsity of the polynomials
 * @param includeNeg: whether to include negative coefficients
 **/
template <class Field, class RecursivePoly>
void ZeroDimensionalRegularChain<Field,RecursivePoly>::randomZeroDimensionalRegularChain(int nVars, int nTrcVars, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg) {
	this->RegularChain<Field,RecursivePoly>::randomRegularChain(nVars,nVars,nTrcVars,maxDegs,coefBound,pSparsity,includeNeg);
}

/// Possible instantiations of the RegularChain class ///
template class ZeroDimensionalRegularChain<RationalNumber,SparseMultivariateRationalPolynomial>;
