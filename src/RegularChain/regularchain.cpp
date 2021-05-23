#include "RegularChain/regularchain.hpp"
#include "RegularChain/zerodimensionalregularchain.hpp"
#include "SubResultantChain/subresultantchain.hpp"
#include "RationalNumberPolynomial/mrpolynomial.h"
#include "Utils/RandomHelpers.hpp"
#include "Utils/SymbolHelpers.hpp"
#include "Utils/Parallel/SynchronizedWriteVector.hpp"
#include "Utils/Parallel/TaskScheduler.hpp"
#include "Utils/Unix_Timer.h"
#include "FiniteFields/prime32_constants.h"
#include "global.h"


#define WITH_FACTORIZATION 1

#define RC_TRIANGULARIZE_TREEDATA 0

#define RC_TRIANGULARIZE_PROGRESS 0


extern float SMQP_SRC_TIME;

int intersectDepth = 0;
int regularGCDDepth = 0;
int intersectFreeDepth = 0;
int intersectAlgebraicDepth = 0;
int regularizeDepth = 0;
int regularizeSingleDepth = 0;
int extendDepth = 0;
int cleanChainDepth = 0;
int depth = 0;

int totalIntersects = 0;
int totalSplits = 0;

float reduceMinimalTime(0);
float primitivePartTime(0);
float squareFreePartTime(0);
float subresultantChainTime(0);
float zerodimensionalregularchainTime(0);
float pseudoDivideTime(0);
float normalFormTime(0);
float removeRedundantChainsTime(0);
float factorTime(0);

float triangularizeTime(0);
float intersectTime(0);
float regularGCDTime(0);
float intersectFreeTime(0);
float intersectAlgebraicTime(0);
float regularizeTime(0);
float extendTime(0);
float cleanChainTime(0);
float constructChainTime(0);
float constructChainsTime(0);
float GCDFreeFactorizationTime(0);
float ZDIntersectTime(0);
float ZDRegularizeTime(0);
float isInSaturatedIdealTime(0);
float cleanSetTime(0);

float tsCopyTime(0);
float rcCopyTime(0);
float zdrcCopyTime(0);

// Forward declares
template <class Field, class RecursivePoly>
std::vector<Symbol> findOrderedRing(const std::vector<RecursivePoly>& F);


template <class T>
using SyncVector_ptr = std::shared_ptr<SynchronizedWriteVector<T>>;


/// Some Protected Functions ///


template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::updateRegularChainStates() {
	TriangularSet<Field,RecursivePoly>::updateTriangularSetStates();
	// Primality test: univariate poly irreducible, multivariate polys degree one

	// Squarefree test?
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::updateRegularChainStates(const RecursivePoly& p) {
	TriangularSet<Field,RecursivePoly>::updateTriangularSetStates(p);
	// Primality test:

	// Squarefree test?
}

/// Constructors ///

/**
 * Default constructor: creates an empty triangular set of variable size
 * with empty list of transcendentals
 *
 * @param
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain () : TriangularSet<Field,RecursivePoly>(), squareFreeLevel(0), saturatedIdealPrimeLevel(0) {}

/**
 * Construct an empty triangular set of fixed size in the decreasingly ordered
 * variables given by xs with empty list of transcendentals
 *
 * @param xs: The variable names
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const std::vector<Symbol>& xs) : TriangularSet<Field,RecursivePoly>(xs), squareFreeLevel(0), saturatedIdealPrimeLevel(0) {}

/**
 * Construct an empty triangular set of fixed size in the decreasingly ordered
 * variables given by xs and list of transcendentals given by ts
 *
 * @param xs: The variable names
 * @param ts: The transcendental variable names
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const std::vector<Symbol>& xs, const std::vector<Symbol>& ts) : TriangularSet<Field,RecursivePoly>(xs,ts), squareFreeLevel(0), saturatedIdealPrimeLevel(0) {}

/**
 * Construct a variable triangular set containing p, such that the variables of p are treated as algebraic,
 * with empty list of transcendentals
 *
 * @param p: The polynomial to add
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const RecursivePoly& p) : TriangularSet<Field,RecursivePoly>(p) {
	long long unsigned int rcProfilingStart;
	int index(variableIndex(p.leadingVariable()));
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "index = " << index << std::endl;
		std::cerr << "computing primitive part of " << set[index] << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	set[index] = set[index].mainPrimitivePart();
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
	#endif
	if ((regularChainOptions & MAINTAIN_SQUAREFREE) == MAINTAIN_SQUAREFREE) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "computing squarefree part of " << set[index] << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		set[index] = set[index].squareFreePart();
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
		squareFreeLevel = 1;
	} else if (p.degree() == 1) {
		squareFreeLevel = 1;
	} else {
		squareFreeLevel = 0;
	}

	//TODO: Maintain prime?
	if (p.degree() == 1) {
		saturatedIdealPrimeLevel = 1;
	} else {
		saturatedIdealPrimeLevel = 0;
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
RegularChain<Field,RecursivePoly>::RegularChain (const RecursivePoly& p, const std::vector<Symbol>& ts) : TriangularSet<Field,RecursivePoly>(p,ts) {
	long long unsigned int rcProfilingStart;
	int index(variableIndex(p.leadingVariable()));
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "index = " << index << std::endl;
		std::cerr << "computing primitive part of " << set[index] << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	set[index] = set[index].mainPrimitivePart();
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
	#endif
	if ((regularChainOptions & MAINTAIN_SQUAREFREE) == MAINTAIN_SQUAREFREE) {
		#ifdef REGULARCHAIN_DEBUG
		std::cerr << "computing squarefree part of " << set[index] << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		set[index] = set[index].squareFreePart();
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
		squareFreeLevel = 1;
	} else if (p.degree() == 1) {
		squareFreeLevel = 1;
	}
	else {
		squareFreeLevel = 0;
	}

	if (p.degree() == 1) {
		saturatedIdealPrimeLevel = 1;
	} else {
		saturatedIdealPrimeLevel = 0;
	}
}

template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const std::vector<RecursivePoly> polys) : TriangularSet<Field,RecursivePoly>() {
	std::vector<Symbol> vs = findOrderedRing<Field,RecursivePoly>(polys);
	*this = RegularChain<Field,RecursivePoly>(vs);
	algVars = vars;
	set = polys;

	squareFreeLevel = 0;
	saturatedIdealPrimeLevel = 0;
}

/**
 * Copy constructor
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) : TriangularSet<Field,RecursivePoly>(a), squareFreeLevel(a.squareFreeLevel), saturatedIdealPrimeLevel(a.saturatedIdealPrimeLevel), regularChainOptions(a.options()) {
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "copy constructor in RC for ZDRC called..." << std::endl;
	#endif
}

/**
 * Copy constructor
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const RegularChain<Field,RecursivePoly>& a) : TriangularSet<Field,RecursivePoly>(a), squareFreeLevel(a.squareFreeLevel), saturatedIdealPrimeLevel(a.saturatedIdealPrimeLevel), regularChainOptions(a.options()) {
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "copy constructor in RC called..." << std::endl;
	#endif
}

/**
 * Copy constructor
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (const TriangularSet<Field,RecursivePoly>& a) : TriangularSet<Field,RecursivePoly>(a) {
	//TODO: Implement proper checking to see if this is a valid operation.
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "copy constructor in RC for TS called...NO VALIDITY CHECK IMPLEMENTED" << std::endl;
	#endif
	//TODO: compute squareFreePart of a if MAINTAIN_SQUAREFREE is imposed
	squareFreeLevel = 0;
	saturatedIdealPrimeLevel = 0;
}

/**
 * Move constructor
 *
 * @param a: An r-value reference triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a) : TriangularSet<Field,RecursivePoly>(std::move(a)), squareFreeLevel(a.squareFreeLevel), saturatedIdealPrimeLevel(a.saturatedIdealPrimeLevel), regularChainOptions(a.options()) {
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "move constructor in RC for ZDRC called..." << std::endl;
	#endif
}

/**
 * Move constructor
 *
 * @param a: An r-value reference triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (RegularChain<Field,RecursivePoly>&& a) : TriangularSet<Field,RecursivePoly>(std::move(a)), squareFreeLevel(a.squareFreeLevel), saturatedIdealPrimeLevel(a.saturatedIdealPrimeLevel), regularChainOptions(a.options()) {
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "move constructor in RC called..." << std::endl;
	#endif
}

/**
 * Move constructor
 *
 * @param a: An r-value reference triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>::RegularChain (TriangularSet<Field,RecursivePoly>&& a) : TriangularSet<Field,RecursivePoly>(std::move(a)) {
	//TODO: Implement proper checking to see if this is a valid operation.
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "move constructor in RC for TS called...NO VALIDITY CHECK IMPLEMENTED" << std::endl;
	#endif
	//TODO: compute squareFreePart of a if MAINTAIN_SQUAREFREE is imposed
	squareFreeLevel = 0;
	saturatedIdealPrimeLevel = 0;
}

/**
 * Computational constructor: creates a regular chain given all the data
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
RegularChain<Field,RecursivePoly>::RegularChain (const std::vector<Symbol>&& vs, const std::vector<Symbol>&& avs, const std::vector<Symbol>&& tvs, const std::vector<RecursivePoly>&& ts, TriangularSetMode tsm, const mpz_class& c) : TriangularSet<Field,RecursivePoly>(std::move(vs),std::move(avs),std::move(tvs),std::move(ts),tsm,c) {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "WARNING! NOT FULLY IMPLEMENTED TO MAINTAIN ATTRIBUTES" << std::endl;
	#endif
	//TODO: compute thigns based on MAINTAIN options.
	squareFreeLevel = 0;
	saturatedIdealPrimeLevel = 0;
}

/**
* Find an ordered subring of ringVars containing all and only the variables of p1Vars and p2Vars.
**/
template <class Field, class RecursivePoly>
std::vector<Symbol> findOrderedSubRing(const std::vector<Symbol>& p1Vars, const std::vector<Symbol>& p2Vars, const std::vector<Symbol>& ringVars) {
	std::vector<Symbol> out;
	out = setUnion(p1Vars,p2Vars);
	out = orderPreservingSetIntersection(ringVars,out);
	return out;
}

/**
* Find the ordered subring of ringVars containing all and only the variables of pVars.
**/
template <class Field, class RecursivePoly>
std::vector<Symbol> findOrderedSubRing(const std::vector<Symbol>& pVars, const std::vector<Symbol>& ringVars) {
	return orderPreservingSetIntersection(ringVars,pVars);
}

template <class Field, class RecursivePoly>
RecursivePoly RegularChain<Field,RecursivePoly>::moduloPolysWithConstantInitial(const RecursivePoly& p_in) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[modPolyConstInit] entering moduloPolysWithConstantInitial" << std::endl;
	#endif

	if (this->isConstantPolynomial(p_in) || this->isEmpty())
		return p_in;

	// printVariables(p_in.ringVariables(),"p_in");
	RecursivePoly p(p_in), poly;
	std::vector<RecursivePoly> polys(this->polynomials());
	std::vector<Symbol> subRingVars;

	//TODO: The below is likely slow, instead sprem(p, RC) is RC is strongly normalized.

	for (size_t i=0; i<polys.size(); ++i) {
		if (this->isConstantPolynomial(p))
			break;
		if (polys[i].isZero())
			continue;
		poly = polys[i];
		if (poly.initial().isConstant() != 0) {
			// TODO: Maple has a sophisticated way of computing the univariate sprem
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[modPolyConstInit] p ringVariables:" << std::endl;
				printVariables(p.ringVariables(),"pvar");
				std::cerr << "[modPolyConstInit] poly ringVariables:" << std::endl;
				printVariables(poly.ringVariables(),"polyvar");
				std::cerr << "[modPolyConstInit] *this variables:" << std::endl;
				printVariables(this->allVariables(),"Tvar");
			#endif
			if (mode == TS_FIXED) {
				subRingVars = findOrderedSubRing<Field,RecursivePoly>(p.variables(),poly.variables(),this->allVariables());
				p.setRingVariables(subRingVars);
				poly.setRingVariables(subRingVars);
				p = p.pseudoDivide(poly, NULL, NULL, 1);
				// TODO: when modular computation is implemented, modularize coefficients here
				p = p.primitivePart(this->vars);
			} else {
				p = p.pseudoDivide(poly, NULL, NULL, 1);
				p = p.primitivePart();
			}
		}
	}

	return p;

}

template <class Field, class RecursivePoly>
RecursivePoly RegularChain<Field,RecursivePoly>::makeCanonical(const RecursivePoly& p_in) const {

	RecursivePoly p(p_in);
	// TODO: would modularizeCoefficients here

	if (this->transcendentalVariables().empty()) { // only executed in characteristic zero
		if (p.leadingCoefficient() < 0)
			p = -p;
	}
	// TODO: handle positive characteristic case here (if initial of p is an integer, replace p by p/init(p) mod prime)

	return p;

}

template <class Field, class RecursivePoly>
std::vector<RecursivePoly> RegularChain<Field,RecursivePoly>::factorPolynomial(const RecursivePoly& p_in) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[FP] entering factorPolynomial: " << p_in << std::endl;
	#endif

	RecursivePoly p = reduceMinimalPrimitivePart(p_in);

	long long unsigned int rcProfilingStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	Factors<SparseMultivariateRationalPolynomial> facts = p.factor();
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&factorTime);
	#endif
	std::vector<RecursivePoly> factored;
	for (size_t i=0; i<facts.size(); ++i) {
		factored.push_back(facts[i].first);
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[FP] done factoring, got nfacts: " << facts.size() << std::endl;
	#endif

	// //Don't do this check once we trust factorization
	// RecursivePoly check,temp;
	// check = facts.ringElement();
	// for (size_t j=0; j<facts.size(); ++j) {
	// 	temp = facts.factor(j).first^facts.factor(j).second;
	// 	check *= temp;
	// }
	// if (check != p) {
	// 	std::cerr << "Regular Chain error, factorization is wrong." << std::endl;
	// 	std::cerr << "check = " << check << std::endl;
	// 	std::cerr << "p = " << p << std::endl;
	// 	std::cerr << "ringElement = " << facts.ringElement() << std::endl;
	// 	for (size_t m=0; m<facts.size(); ++m) {
	// 		std::cerr << "factor[" << m <<  "] = [" << facts.factor(m).first << "," << facts.factor(m).second << "]" << std::endl;
	// 	}
	// 	printVariables(p.ringVariables(), "p ringVars");
	// 	printVariables(check.ringVariables(), "check ringVars");
	// 	exit(1);
	// }
	// #ifdef REGULARCHAIN_DEBUG
	// 	std::cerr << "[FP] done factoring check." << std::endl;
	// #endif
	return factored;

	// RegularChain<Field,RecursivePoly> lrc;
	// this->lower(p.mainVariable(),lrc);

	// std::vector<Symbol> pvars;
	// RecursivePoly content;

	// if (lrc.isZeroDimensionalMathematically()) {
	// 	p = reduceMinimalPrimitivePart(p);
	// 	pvars = p.variables();

	// 	if (pvars.size() == 1) // p is univariate (even in terms of params)
	// 		facts = p.factor();
	// 	else
	// 		facts = p.squareFree();

	// 	for (size_t i=0; i<facts.size(); ++i)
	// 		factored.push_back(facts[i].first);
	// }
	// else {
	// 	p = p.mainPrimitivePart();
	// 	p = removeZero(p,this->allVariables());

	// 	RecursivePoly f;
	// 	std::vector<RecursivePoly> tasks;
	// 	tasks.push_back(p);

	// 	while (tasks.size() != 0) {
	// 		f = tasks.back();
	// 		tasks.pop_back();

	// 		f = reduceMinimalPrimitivePart(f);
	// 		pvars = f.variables();

	// 		if (pvars.size() == 1) // f is univariate (even in terms of params)
	// 			facts = f.factor();
	// 		else
	// 			facts = f.squareFree();

	// 		if (facts.size() == 1)
	// 			factored.push_back(facts[0].first);
	// 		else {
	// 			for (size_t i=0; i<facts.size(); ++i)
	// 				tasks.push_back(facts[i].first);
	// 		}

	// 	}
	// }

	// for (size_t i=0; i<factored.size(); ++i) {
	// 	std::cerr << "output from factorPolynomial:" << std::endl;
	// 	std::cerr << "factored[" << i << "] = " << factored[i] << std::endl;
	// }

	// return factored;

// NOTE: This pseudo-code below is from a time when multivariate
//       factorization was impractical.
//	if (this->isZeroDimensionalMathematically()) {
//
	//  	reduce_primpart(p)
	//  	if p is univariate (even in terms of params)
	//  		factor into irreducibles
	//  	else
	//  		do univariate squarefree factorization wrt mvar(p)
//	}
	//  else
	//  	mainPrimitivePart(p)
	//  	this->removeZero(p)
	//  	assert(!(p.isZero()))
	//  	tasks = [p]
	//  	factored = []
	//  	while tasks != [] do
	//  		f = tasks.pop()
	//  		this->moduloPolysWithConstantInitial(f)
	//  		this->removeZero(f)
	//			mainPrimitivePart(f)
	//			if p is univariate
	//  			factors = factor into irreducibles
	//  		else
	//  			factors = do univariate squarefree factorization wrt mvar(p)
	//  		if factors.size == 1
	//  			factored.push(factors[0])
	//  		else
	//  			pile.insert(factors)
	//
}


/// Member Functions ///

/**
 * Construct a regular chain from the current object and an input polynomial
 *
 * @param p: A polynomial
 **/
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::constructChain(const RecursivePoly& p, int options) {
	long long unsigned int rcProfilingStart,constructChainStart;
	bool isReg;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&constructChainStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering constructChain(p,opts) in RC:" << std::endl;
		std::cerr << "opts = " << options << std::endl;
	#endif

	if ((options & ASSUME_REGULAR) != ASSUME_REGULAR) {

		// perform a regularity test of p.initial() modulo Sat(T)
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&constructChainStart,&constructChainTime);
		#endif
		isReg = isRegular(p.initial());
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&constructChainStart);
		#endif

		if (!isReg) {
			std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "in constructChain(p,opts) of RC" << std::endl;
				std::cerr << "*this = " << *this << std::endl;
				std::cerr << "p = " << p << std::endl;
				std::cerr << "p.initial() = " << p.initial() << std::endl;
			#endif
			exit(1);
		}

	}

	RecursivePoly q(p);
//	q.setRingVariables(allVariables());

	if ((options & ASSUME_PRIMITIVE) != ASSUME_PRIMITIVE) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "calling mainPrimitivePart on " << q << " in constructChain(p,opts) of RC class." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		q = q.mainPrimitivePart();
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
	}

	if ((options & ASSUME_SQUAREFREE) != ASSUME_SQUAREFREE) {
		// if we are not maintaining squarefree rcs (at chain level) and we are not
		// assuming the input is squarefree (at the call level), then compute the
		// squarefree part over the base field to reduce expression swell.
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "calling squareFreePart on " << q << " in constructChain(p,opts) of RC class." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		q = q.squareFreePart();
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "returning from squareFreePart call to SMQP" << std::endl;
		#endif
		if ((regularChainOptions & MAINTAIN_SQUAREFREE) == MAINTAIN_SQUAREFREE) {
			std::cerr << "BPAS: warning, if we cannot assume that input is squarefree then we cannot maintain squarefree regular chains without potential splitting; change input options or call ConstructChains instead." << std::endl;
//			exit(1);
		}
	}

	if (((options & ASSUME_REDUCED) != ASSUME_REDUCED) && (this->size() > 0)) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "calling reduce on " << q << " in constructChain(p,opts) of RC class." << std::endl;
			printVariables(vars,"vars");
			printVariables(algVars,"algVars");
		#endif

		q = this->reduceMinimal(q);
//		q = this->reduce(q);

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "calling mainPrimitivePart again on " << q << " in constructChain(p,opts) of RC class." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif

		q = q.mainPrimitivePart();

		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
		#ifdef REGULARCHAIN_DEBUG
		std::cerr << "calling squareFreePart on " << q << " in constructChain(p,opts) of RC class." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif

		q = q.squareFreePart();

		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "returning from squareFreePart call to SMQP" << std::endl;
		#endif
	}

	if ((options & ASSUME_SQUAREFREE)) {
		if (isSquareFree()) {
			++squareFreeLevel;
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "calling operator+=(q) in TS." << std::endl;
	#endif
		TriangularSet<Field,RecursivePoly>::operator+=(q);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "leaving constructChain(p,opts) in RC." << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&constructChainStart,&constructChainTime);
	#endif
}

/**
 * Construct a regular chain from the current object and an input upper regular chain
 *
 * @param T: An upper regular chain
 **/
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::constructChain(const RegularChain<Field,RecursivePoly>& T, int options) {
	// TODO: check that the chain is upper
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering constructChain(p,opts) in RC." << std::endl;
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
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::constructChainsFromPoly(const RecursivePoly& p, bool lazardDecompose, int heightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::constructChainsFromPoly(const RecursivePoly& p, bool lazardDecompose, int heightBound, int options) const {
#endif
	long long unsigned int rcProfilingStart,constructChainsStart;
	bool isReg;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&constructChainsStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering RC constructChainsFromPoly" << std::endl;
		std::cerr << "heightBound = " << heightBound << std::endl;
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_CONSTRUCTCH_OBJ, results);
	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_SFP_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_SFP_OBJ, currSFPTask);
	std::vector<RecursivePoly> factorTasks;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "regularity check" << std::endl;
	#endif
	if ((options & ASSUME_REGULAR) != ASSUME_REGULAR) {

		// perform a regularity test of p.initial() modulo Sat(T)
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&constructChainsStart,&constructChainsTime);
		#endif
		isReg = isRegular(p.initial());
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&constructChainsStart);
		#endif

		if (!isReg) {
			std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "in constructChainsFromPoly of RC." << std::endl;
				std::cerr << "*this = " << *this << std::endl;
				std::cerr << "p = " << p << std::endl;
				std::cerr << "p.initial() = " << p.initial() << std::endl;
			#endif
			exit(1);
		}
	}

	RecursivePoly q(p);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "allVariables():" << std::endl;
		printVariables(allVariables(),"allvars");
	#endif
	if (mode == TS_FIXED) {
		q.setRingVariables(allVariables());
	} else {
		q.setRingVariables(orderPreservingSetUnion(q.ringVariables(), allVariables()));
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "primitive part" << std::endl;
	#endif
	if ((options & ASSUME_PRIMITIVE) != ASSUME_PRIMITIVE) {
		#ifdef REGULARCHAIN_DEBUG
		std::cerr << "calling primitivePart on " << q << " in constructChainsFromPoly of RC class." << std::endl;
//		std::cerr << "q.ringVariables:" << std::endl;
//		printVariables(q.ringVariables(),"qrv");
		#endif

		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		q = q.mainPrimitivePart();
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "primitivePart computed." << std::endl;
			std::cerr << "q = " << q << std::endl;
//			std::cerr << "content = " << q.content() << std::endl;
		#endif
	}

	// TODO: handle the case of maintaining normalized regular chains (see lines 481-508 of regularchains.mm)
	//       The idea here is compute the inverse of the initial of q, which can split the computation,
	//         and then for each branch form a normalized q by reconstructing the head of q and adding
	//         the result of removeZero on the product of tail(q) and the inverse of init(q). The integer
	//         content is removed before adding the normalized q to the component.

//	if ((options & ASSUME_REDUCED) != ASSUME_REDUCED) {
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling reduce on " << q << " in constructChainsFromPoly of RC class." << std::endl;
////			printVariables(vars,"vars");
////			printVariables(algVars,"algVars");
//		#endif
//
//		q = this->reduce(q);
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling primitivePart again on " << q << " in constructChainsFromPoly of RC class." << std::endl;
//		#endif
//
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&rcProfilingStart);
//		#endif
//		q = q.mainPrimitivePart();
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
//		#endif
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "returning to constructChain from primitivePart" << std::endl;
//		#endif
//	}

	// TODO: When normalization is included, another layer of loop needs to be added here.
	// TODO: In Maple, Lazard decomposition and normalizing chains are incompatible; is
	//       this true here?

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "squarefree part" << std::endl;
	#endif
	if (((options & ASSUME_SQUAREFREE) != ASSUME_SQUAREFREE) && ((regularChainOptions & MAINTAIN_SQUAREFREE) == MAINTAIN_SQUAREFREE)) {

		// TODO: maybe?
		// assert(this->isSquareFree && "Construct polys from chain should be square free to start");

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "calling squareFreePart(q,v,hb) in constructChain, v = " << q.leadingVariable() << std::endl;
//			printVariables(q.variables(),"q");
			std::cerr << "calling squareFreePart on " << q << " in constructChainsFromPoly of RC class." << std::endl;
		#endif

		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		RC_GEN_CONSUMER_INIT(RC_SFP_OBJ, squareFreePartTasks, (&RegularChain<Field, RecursivePoly>::_squareFreePartPolynomial), this, q, q.leadingVariable(), lazardDecompose, heightBound, options);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "returning to constructChainsFromPoly from squareFreePart call" << std::endl;
		#endif

		RC_GEN_CONSUMER_LOOP(squareFreePartTasks, currSFPTask, i) {

			RC_GEN_CONSUMER_GET_LOOPELEM(squareFreePartTasks, currSFPTask, i);

			if ((regularChainOptions & CONSTRUCT_FACTORIZE) == CONSTRUCT_FACTORIZE) {

				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "calling factorPolynomial on " << currSFPTask.poly << " wrt " << currSFPTask.chain << " in constructChainsFromPoly of RC class." << std::endl;
				#endif
				factorTasks = currSFPTask.chain.factorPolynomial(currSFPTask.poly);

				for (auto poly : factorTasks) {

					RegularChain<Field,RecursivePoly> rc(currSFPTask.chain);

					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "calling operator+=(q) in TS after makeCanonical." << std::endl;
					#endif
					poly = makeCanonical(poly);
					rc.TriangularSet<Field,RecursivePoly>::operator+=(poly);
					// TODO: when TS::op+=(p) is changed to TS::op+=(p,opts,bool), use regularChainOptions and eliminate next line (obvs); similarly elsewhere in this routine
					//if this was already square free, then no changes needed. It it wasn't, then the result also is not.

					if (currSFPTask.chain.isSquareFree()) {
						rc.squareFreeLevel = currSFPTask.chain.squareFreeLevel + 1;
					} else if (rc.variableIndex(p.leadingVariable()) > rc.variableIndex(rc.algVars[0])) {
						//else if  p is below the current algvars
						rc.squareFreeLevel = 0;
					}

					if ((currSFPTask.chain.isSaturatedIdealPrime() && poly.mainDegree() == 1) || currSFPTask.chain.isEmpty()) {
						rc.saturatedIdealPrimeLevel = currSFPTask.chain.saturatedIdealPrimeLevel + 1;
					} else if (rc.variableIndex(p.leadingVariable()) > rc.variableIndex(rc.algVars[0])) {
						rc.saturatedIdealPrimeLevel = 0;
					}

//					std::cerr << "constructChainsFromPoly output 1:" << p << std::endl;
//					std::cerr << "rc = " << rc << std::endl;
					RC_GEN_PRODUCER_ACCUMULATE(results, rc);
				}
			}
			else {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "calling operator+=(q) in TS after makeCanonical and reduceMinimalPrimitiveSquareFreePart." << std::endl;
				#endif
				currSFPTask.poly = makeCanonical(currSFPTask.chain.reduceMinimalPrimitiveSquareFreePart(currSFPTask.poly));

				bool wasSquareFree = currSFPTask.chain.isSquareFree();
				bool wasSatPrime = currSFPTask.chain.isSaturatedIdealPrime();

				currSFPTask.chain.TriangularSet<Field,RecursivePoly>::operator+=(currSFPTask.poly);

//				std::cerr << "constructChainsFromPoly output 2:" << p << std::endl;
//				std::cerr << "currSFPTask.chain = " << currSFPTask.chain << std::endl;
				if (wasSatPrime && currSFPTask.poly.mainDegree() == 1) {
					currSFPTask.chain.saturatedIdealPrimeLevel += 1;
				} else if (currSFPTask.chain.variableIndex(p.leadingVariable()) > currSFPTask.chain.variableIndex(currSFPTask.chain.algVars[0])) {
					currSFPTask.chain.saturatedIdealPrimeLevel = 0;
				}
				if (wasSquareFree) {
					currSFPTask.chain.squareFreeLevel += 1;
				} else if (currSFPTask.chain.variableIndex(p.leadingVariable()) > currSFPTask.chain.variableIndex(currSFPTask.chain.algVars[0])) {
					currSFPTask.chain.squareFreeLevel = 0;
				}

				RC_GEN_PRODUCER_ACCUMULATE(results, currSFPTask.chain);
			}
		}
	}
	else {
//		RegularChain<Field,RecursivePoly> T(*this);
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling operator+=(q) in TS." << std::endl;
//		#endif
//		T.TriangularSet<Field,RecursivePoly>::operator+=(q);
//		RC_GEN_PRODUCER_ACCUMULATE(results, T);

		// TODO: add new options ASSUME_SQUAREFREE_MODULO_CHAIN and ASSUME_SQUAREFREE_POLYNOMIAL
		//       to replace ASSUME_SQUAREFREE and thereby clarify the semantics, and then
		//       re-factor/re-write this routine to account properly for these cases
		// TODO: Jan23/2020: Let's say that ASSUME_SQUAREFREE always means wrt the RC since
		//       this is the brutal part.

		if (((options & CONSTRUCT_FACTORIZE) == CONSTRUCT_FACTORIZE) || ((regularChainOptions & CONSTRUCT_FACTORIZE) == CONSTRUCT_FACTORIZE)) {

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "calling factorPolynomial on " << q << " wrt " << *this << " in constructChainsFromPoly of RC class." << std::endl;
			#endif
			factorTasks = factorPolynomial(q);

			for (auto poly : factorTasks) {

				poly = makeCanonical(reduceMinimal(poly));
				RegularChain<Field,RecursivePoly> rc(*this);

				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "calling operator+=(q) in TS after makeCanonical." << std::endl;
				#endif
				rc.TriangularSet<Field,RecursivePoly>::operator+=(poly);
				if (this->isSquareFree() && (poly.mainDegree() == 1 || options & ASSUME_SQUAREFREE)) {
					rc.squareFreeLevel = this->squareFreeLevel + 1;
				} else if (rc.variableIndex(p.leadingVariable()) > rc.variableIndex(rc.algVars[0])) {
					//p is not being added on top of this
					rc.squareFreeLevel = 0;
				}
				if ((this->isSaturatedIdealPrime() && poly.mainDegree() == 1) || this->isEmpty()) {
					rc.saturatedIdealPrimeLevel = this->saturatedIdealPrimeLevel + 1;
				} else if (rc.variableIndex(p.leadingVariable()) > rc.variableIndex(rc.algVars[0])) {
					rc.saturatedIdealPrimeLevel = 0;
				}

//				std::cerr << "constructChainsFromPoly output 3:" << p << std::endl;
//				std::cerr << "rc = " << rc << std::endl;
				RC_GEN_PRODUCER_ACCUMULATE(results, rc);
			}
		}
		else {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "calling operator+=(q) in TS after makeCanonical and reduceMinimalPrimitiveSquareFreePart." << std::endl;
			#endif

			q = makeCanonical(reduceMinimalPrimitiveSquareFreePart(q));
			RegularChain<Field,RecursivePoly> rc(*this);
			rc.TriangularSet<Field,RecursivePoly>::operator+=(q);
			if (this->isSquareFree() && (q.mainDegree() == 1 || options & ASSUME_SQUAREFREE)) {
				rc.squareFreeLevel = this->squareFreeLevel + 1;
			} else if (rc.variableIndex(p.leadingVariable()) > rc.variableIndex(algVars[0])) {
				//p is not being added on top of this
				rc.squareFreeLevel = 0;
			}
			if (this->isSaturatedIdealPrime() && q.mainDegree() == 1) {
				rc.saturatedIdealPrimeLevel = this->saturatedIdealPrimeLevel + 1;
			} else if (rc.variableIndex(p.leadingVariable()) > rc.variableIndex(rc.algVars[0])) {
				rc.saturatedIdealPrimeLevel = 0;
			}

//			std::cerr << "constructChainsFromPoly output 4:" << p << std::endl;
//			std::cerr << "rc = " << rc << std::endl;
			RC_GEN_PRODUCER_ACCUMULATE(results, rc);
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "leaving RC constructChainsFromPoly in RC" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&constructChainsStart,&constructChainsTime);
	#endif
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

///**
// * Construct a set of regular chains from an input regular chain and polynomial
// *
// * @param T: A regular chain
// * @param p: A polynomial
// **/
//#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
//template <class Field, class RecursivePoly>
//void RegularChain<Field,RecursivePoly>::constructChainsFromPoly(const RecursivePoly& p, bool lazardDecompose, int heightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
//#else
//template <class Field, class RecursivePoly>
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::constructChainsFromPoly(const RecursivePoly& p, bool lazardDecompose, int heightBound, int options) const {
//#endif
//	long long unsigned int rcProfilingStart,constructChainsStart;
//	bool isReg;
//	#ifdef REGULARCHAIN_PROFILING
//		startTimer(&constructChainsStart);
//	#endif
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "entering RC constructChainsFromPoly" << std::endl;
//		std::cerr << "heightBound = " << heightBound << std::endl;
//	#endif
//	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
//	RC_GEN_PRODUCER_DECLARE_RESULT(RC_CONSTRUCTCH_OBJ, results);
//	typedef RC_CONSTRUCTCH_OBJ RC_SFP_OBJ;
//	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_SFP_OBJ, currSquareFreePartTask);
//	std::vector<RegularChain<Field,RecursivePoly>> moreResults;
////	std::vector<RegularChain<Field,RecursivePoly>> results,moreResults;
////	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularComponents;

////	if ((regularChainOptions & ASSUME_REGULAR) != ASSUME_REGULAR) {
////		regularComponents = regularize(p);
////
////		for (size_t i=0; i<regularComponents.size(); ++i) {
////			if (!regularComponents[i].poly.isZero()) {
////				moreResults = regularComponents[i].chain.constructChains(regularComponents[i].poly, options | ASSUME_REGULAR);
////				results.insert(results.end(),moreResults.begin(),moreResults.end());
////			}
////		}
////		return results;
////	}
//	// TODO: Find a way to avoid the code repeated between here and constructChain

//	if ((options & ASSUME_REGULAR) != ASSUME_REGULAR) {
//
//		// perform a regularity test of p.initial() modulo Sat(T)
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&constructChainsStart,&constructChainsTime);
//		#endif
//		isReg = isRegular(p.initial());
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&constructChainsStart);
//		#endif
//
//		if (!isReg) {
//			std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
//			#ifdef REGULARCHAIN_DEBUG
//				std::cerr << "in constructChainsFromPoly of RC." << std::endl;
//				std::cerr << "*this = " << *this << std::endl;
//				std::cerr << "p = " << p << std::endl;
//				std::cerr << "p.initial() = " << p.initial() << std::endl;
//			#endif
//			exit(1);
//		}
//	}
//
//	RecursivePoly q(p);
////	q.setRingVariables(allVariables());
//
//	if ((options & ASSUME_PRIMITIVE) != ASSUME_PRIMITIVE) {
//		#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "calling primitivePart on " << q << " in constructChainsFromPoly of RC class." << std::endl;
////		std::cerr << "q.ringVariables:" << std::endl;
////		printVariables(q.ringVariables(),"qrv");
//		#endif
//
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&rcProfilingStart);
//		#endif
//		q = q.mainPrimitivePart();
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
//		#endif
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "primitivePart computed." << std::endl;
//			std::cerr << "q = " << q << std::endl;
////			std::cerr << "content = " << q.content() << std::endl;
//		#endif
//	}
//
//	// TODO: handle the case of maintaining normalized regular chains (see lines 481-508 of regularchains.mm)
//	//       The idea here is compute the inverse of the initial of q, which can split the computation,
//	//         and then for each branch form a normalized q by reconstructing the head of q and adding
//	//         the result of removeZero on the product of tail(q) and the inverse of init(q). The integer
//	//         content is removed before adding the normalized q to the component.
//
//	if ((options & ASSUME_REDUCED) != ASSUME_REDUCED) {
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling reduce on " << q << " in constructChainsFromPoly of RC class." << std::endl;
////			printVariables(vars,"vars");
////			printVariables(algVars,"algVars");
//		#endif
//
//		q = this->reduce(q);
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling primitivePart again on " << q << " in constructChainsFromPoly of RC class." << std::endl;
//		#endif
//
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&rcProfilingStart);
//		#endif
//		q = q.mainPrimitivePart();
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
//		#endif
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "returning to constructChain from primitivePart" << std::endl;
//		#endif
//	}
//
//	// TODO: When normalization is included, this needs to be "loopified" to compute the
//	//       squarefree part for each component.
//	// TODO: In Maple, Lazard decomposition and normalizing chains are incompatible; is
//	//       this true here?
//
//	if ((options & ASSUME_SQUAREFREE) != ASSUME_SQUAREFREE) {
////		if (((options & ASSUME_REDUCED) != ASSUME_REDUCED) || ((options & ASSUME_PRIMITIVE) != ASSUME_PRIMITIVE)) {
////			#ifdef REGULARCHAIN_DEBUG
////				std::cerr << "calling primitivePart for squareFreePart on " << q << " in constructChainsFromPoly of RC class." << std::endl;
////			#endif
////
////			#ifdef REGULARCHAIN_PROFILING
////				startTimer(&rcProfilingStart);
////			#endif
////			q = q.mainPrimitivePart();
////			#ifdef REGULARCHAIN_PROFILING
////				stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
////			#endif
////
////			#ifdef REGULARCHAIN_DEBUG
////				std::cerr << "returning to constructChain from primitivePart" << std::endl;
////			#endif
////		}
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling squareFreePart(q,v,hb) in constructChain, v = " << q.leadingVariable() << std::endl;
////			printVariables(q.variables(),"q");
//			std::cerr << "calling squareFreePart on " << q << " in constructChainsFromPoly of RC class." << std::endl;
//		#endif
//
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&rcProfilingStart);
//		#endif
//		RC_GEN_CONSUMER_INIT(RC_SFP_OBJ, squareFreePartTasks, (&RegularChain<Field, RecursivePoly>::_squareFreePart), this, q, q.leadingVariable(), lazardDecompose, heightBound, options);
////		moreResults = _squareFreePart(q,q.leadingVariable(),lazardDecompose,heightBound,options);
////		results = squareFreePart(q,q.leadingVariable(),lazardDecompose,heightBound,options);
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
//		#endif

//		RC_GEN_CONSUMER_LOOP(squareFreePartTasks, currSquareFreePartTask, i) {
////		for (auto mr : moreResults) {
//			RC_GEN_CONSUMER_GET_LOOPELEM(squareFreePartTasks, currSquareFreePartTask, i);
//			RC_GEN_PRODUCER_ACCUMULATE(results, currSquareFreePartTask);
//		}
//
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "returning to constructChainsFromPoly from squareFreePart call" << std::endl;
//		#endif
//	}
//	else {
//		RegularChain<Field,RecursivePoly> T(*this);
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling operator+=(q) in TS." << std::endl;
//		#endif
//		T.TriangularSet<Field,RecursivePoly>::operator+=(q);
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "leaving constructChainsFromPoly in RC." << std::endl;
//		#endif
//		RC_GEN_PRODUCER_ACCUMULATE(results, T);
////		results.emplace_back(T);
//	}
//
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "leaving RC constructChainsFromPoly" << std::endl;
//	#endif
//	#ifdef REGULARCHAIN_PROFILING
//		stopTimerAddElapsed(&constructChainsStart,&constructChainsTime);
//	#endif
//	RC_GEN_PRODUCER_COMPLETE(results);
////	return results;
//}

/**
 * Construct a set of regular chains from the current object and an input regular chain above the current object
 *
 * @param T: A regular chain with no overlapping algebraic variables to those of the current object
 **/
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::constructChainsFromChain(const RegularChain<Field,RecursivePoly>& T, bool lazardDecompose, int heightBound, int options) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering RC constructChainsFromChain" << std::endl;
		std::cerr << "heightBound = " << heightBound << std::endl;
	#endif

	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);
	std::vector<RegularChain<Field,RecursivePoly>> results,tasks,moreTasks;
	std::vector<RecursivePoly> polys = T.polynomials(),polys2;
	std::vector<Symbol> mainVars(this->mainVariables());

	if (T.isEmpty()) {
		results.push_back(*this);
		return results;
	}

	// TODO: implement a mainPolynomials() or nonZeroPolynomials() function in TriangularSet
	for (size_t i=0; i<polys.size(); ++i) {
		if (!polys[i].isZero()) {
			polys2.emplace_back(std::move(polys[i]));
		}
	}

	if ((options & ASSUME_MAKESCHAIN) != ASSUME_MAKESCHAIN) {
		int ind = variableIndex(polys2.back().leadingVariable());
		for (size_t j=0; j<mainVars.size(); ++j) {
			if (ind <= variableIndex(mainVars[j])) {
				std::cerr << "BPAS: error, input regular chain must be entirely above the current object" << std::endl;
				exit(1);
			}
		}
	}

//	std::cerr << "calling constructChainsFromPoly in constructChainsFromChain for the first time..." << std::endl;
	RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), this, polys2.back(), lazardDecompose, heightBound, options);
	RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, i);
		results.push_back(currConstructChainsTask);
	}
//	results = constructChains(polys2.back(),lazardDecompose,heightBound,options);
	if (polys2.size() > 1) {
		for (int i=polys2.size()-2; i>-1; --i) {
			for (size_t j=0; j<results.size(); ++j) {
//				std::cerr << "calling constructChainsFromPoly in constructChainsFromChain" << std::endl;
				RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &(results[j]), polys2[i], lazardDecompose, heightBound, options);
				RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, k) {
					RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, k);
//					std::cerr << "pushing result from constructChainsFromChain" << std::endl;
					tasks.push_back(currConstructChainsTask);
				}
//				moreTasks = results[j].constructChains(polys2[i],lazardDecompose,heightBound,options);
//				tasks.insert(tasks.end(),moreTasks.begin(),moreTasks.end());
			}
			results = std::move(tasks);
			tasks.clear();
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "leaving RC constructChainsFromChain" << std::endl;
	#endif
	return results;
}

/**
 * Construct a set of regular chains from an input triangular set T by calling triangularize on the set of polynomials of T
 *
 * @param T: A triangular set
 **/
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::constructChains(const TriangularSet<Field,RecursivePoly>& T) {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering RC constructChains(T) NOT YET TESTED!" << std::endl;
	#endif
	std::vector<RegularChain<Field,RecursivePoly>> results;
//	std::vector<RecursivePoly> polys = T.polynomials();

	// TODO: Consider adding a query to T for properties to use when initializing rc, which will require defining regularChainOptions for triangular set class.
	RegularChain<Field,RecursivePoly> rc;
	results = rc.triangularize(T.polynomials());

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "leaving RC constructChains(T)" << std::endl;
	#endif
	return results;
}


/**
 * Assignment operator =
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) {
	squareFreeLevel = a.squareFreeLevel;
	saturatedIdealPrimeLevel = a.saturatedIdealPrimeLevel;
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	TriangularSet<Field,RecursivePoly>::operator=(a);
	return *this;
}


/**
 * Assignment operator =
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (const RegularChain<Field,RecursivePoly>& a) {
	squareFreeLevel = a.squareFreeLevel;
	saturatedIdealPrimeLevel = a.saturatedIdealPrimeLevel;
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	TriangularSet<Field,RecursivePoly>::operator=(a);
	return *this;
}

template <class Field, class RecursivePoly>
RegularChain<Field, RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (const TriangularSet<Field,RecursivePoly>& a) {
	TriangularSet<Field,RecursivePoly>::operator=(a);
	return *this;
}

template <class Field, class RecursivePoly>
RegularChain<Field, RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (TriangularSet<Field,RecursivePoly>&& a) {
	TriangularSet<Field,RecursivePoly>::operator=(a);
	return *this;
}

/**
 * Assignment operator =
 *
 * @param a: A BPASTriangularSet
 **/
template <class Field, class RecursivePoly>
BPASTriangularSet<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (const BPASTriangularSet<Field,RecursivePoly>& a) {
	if (dynamic_cast<const RegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<const RegularChain<Field,RecursivePoly>&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RegularChain."));
	return *this;
}

/**
 * Assignment operator =
 *
 * @param a: A BPASRegularChain
 **/
template <class Field, class RecursivePoly>
BPASRegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (const BPASRegularChain<Field,RecursivePoly>& a) {
	if (dynamic_cast<const RegularChain<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<const RegularChain<Field,RecursivePoly>&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASRegularChain to RegularChain."));
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a) {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "RegularChain move assignment operator" << std::endl;
	#endif
	squareFreeLevel = a.squareFreeLevel;
	saturatedIdealPrimeLevel = a.saturatedIdealPrimeLevel;
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	TriangularSet<Field,RecursivePoly>::operator=(std::move(a));
	return  *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (RegularChain<Field,RecursivePoly>&& a) {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "calling move constructor from triangularSet..." << std::endl;
	#endif
	squareFreeLevel = a.squareFreeLevel;
	saturatedIdealPrimeLevel = a.saturatedIdealPrimeLevel;
	regularChainOptions = a.regularChainOptions;
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RegChain_UniqueID = a.RegChain_UniqueID;
#endif
	TriangularSet<Field,RecursivePoly>::operator=(std::move(a));
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "completed move in operator=." << std::endl;
	#endif
	return  *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A BPASTriangularSet
 **/
template <class Field, class RecursivePoly>
BPASTriangularSet<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (BPASTriangularSet<Field,RecursivePoly>&& a) {
	if (dynamic_cast<RegularChain<Field,RecursivePoly>*>(&a)) {
		*this = dynamic_cast<RegularChain<Field,RecursivePoly>&&>(a);
	}
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RegularChain."));
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A BPASRegularChain
 **/
template <class Field, class RecursivePoly>
BPASRegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator= (BPASRegularChain<Field,RecursivePoly>&& a) {
	if (dynamic_cast<RegularChain<Field,RecursivePoly>*>(&a)) {
		*this = dynamic_cast<RegularChain<Field,RecursivePoly>&&>(a);
	}
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RegularChain."));
	return *this;
}

/**
 * Overload operator +
 * Adds a polynomial to a regular chain and returns a new regular chain
 *
 * @param p: A sparse multivariate polynomial
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly> RegularChain<Field,RecursivePoly>::operator+ (const RecursivePoly& p) const {
	RegularChain<Field,RecursivePoly> r(*this);
	return (r += p);
}

/**
 * Overload operator +=
 * Adds a polynomial to a regular chain
 *
 * @param p: A recursively viewed polynomial
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator+= (const RecursivePoly& p) {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering RC +=(p)" << std::endl;
		std::cerr << "p = " << p << std::endl;
	#endif

	if (!isRegular(p.initial())) {
		std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "in operator+=(p) of RC" << std::endl;
		#endif
		exit(1);
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "calling constructChain(p) from +=(p)" << std::endl;
	#endif
	constructChain(p,ASSUME_REGULAR);

	return *this;
}

/**
 * Add operator +
 * Adds a polynomial to a regular chain and returns a new regular chain
 *
 * @param p: A sparse multivariate polynomial
 **/
template <class Field, class RecursivePoly>
RegularChain<Field,RecursivePoly> RegularChain<Field,RecursivePoly>::operator+ (const RegularChain<Field,RecursivePoly>& T) const {
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
RegularChain<Field,RecursivePoly>& RegularChain<Field,RecursivePoly>::operator+= (const RegularChain<Field,RecursivePoly>& T) {
	// this routine assumes that the main variables of T are greater than all of the main variables of the current object
	// this assumption may need to be verified later on, perhaps with a flag that avoids the need to do verifications (added to regularChainOptions)

	std::vector<RecursivePoly> polys(T.polynomials());
	// TODO: consider changing this to modify the current object instead
	RegularChain<Field,RecursivePoly> out(T);

	for (int i=polys.size()-1; i>-1; --i) {
		if (!polys[i].isZero()) {
			out += polys[i];
		}
	}

	*this = out;
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
bool RegularChain<Field,RecursivePoly>::operator== (const RegularChain<Field,RecursivePoly>& a) const {
	return TriangularSet<Field,RecursivePoly>::operator==(a);
}

/**
 * Overload comparison operator !=
 *
 *
 * @param a: A regular chain
 **/
template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::operator!= (const RegularChain<Field,RecursivePoly>& a) const {
	return !(*this==a);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const {
	TriangularSet<Field,RecursivePoly> ts2;
	ts2 = std::move(ts);
	TriangularSet<Field,RecursivePoly>::lower(s,ts2);
	ts = std::move(ts2);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::lower(const Symbol& s, RegularChain<Field,RecursivePoly>& rc) const {
	TriangularSet<Field,RecursivePoly> tmp;
	TriangularSet<Field,RecursivePoly>::lower(s, tmp);
	rc = std::move(tmp);
	rc.squareFreeLevel = MIN(rc.algVars.size(), this->squareFreeLevel);
	rc.saturatedIdealPrimeLevel = MIN(rc.algVars.size(), this->saturatedIdealPrimeLevel);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const {
	TriangularSet<Field,RecursivePoly> ts2;
	ts2 = ts;
	TriangularSet<Field,RecursivePoly>::upper(s,ts2);
	RegularChain<Field,RecursivePoly> tmpRC = RegularChain<Field,RecursivePoly>(std::move(ts2));
	ts = std::move(tmpRC);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::upper(const Symbol& s, RegularChain<Field,RecursivePoly>& ts) const {
	TriangularSet<Field, RecursivePoly> ts2;
	TriangularSet<Field,RecursivePoly>:: upper(s, ts2);
	ts = std::move(ts2);
	if (this->isSquareFree()) {
		ts.squareFreeLevel = ts.algVars.size();
	}
	if (this->isSaturatedIdealPrime()) {
		ts.saturatedIdealPrimeLevel = ts.algVars.size();
	}
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::lowerSlice(const Symbol& s) {
	if (!this->isEmpty()) {
		int index(variableIndex(s)+1);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "index = " << index << std::endl;
//			printVariables(vars,"vars");
			std::cerr << "s = " << s << std::endl;
		#endif
		set.erase(set.begin(),set.begin()+index);
		vars.erase(vars.begin(),vars.begin()+index);
	}
}

template <class Field, class RecursivePoly>
int variableIndex(const Symbol& s,const std::vector<Symbol>& vars) {
	int a = std::find(vars.begin(),vars.end(),s) - vars.begin();
	return a;
}

template <class Field, class RecursivePoly>
bool rankLessThan(const RecursivePoly& a, const RecursivePoly& b, const std::vector<Symbol>& vars) {
	// This algorithm assumes that a and b are the ranks of RecursivePoly elements, thus have a single term //
	//  and that both a and b are in the same ordered ring.                                                 //
	bool out;
	int ai,bi;
	ai = variableIndex<Field,RecursivePoly>(a.leadingVariable(),vars);
	bi = variableIndex<Field,RecursivePoly>(b.leadingVariable(),vars);
	if (ai > bi)
		return true;
	else if (ai == bi) {
		int aDeg,bDeg;
		aDeg = a.leadingVariableDegree().get_si();
		bDeg = b.leadingVariableDegree().get_si();
		if (aDeg < bDeg)
			return true;
		else
			return false;
	}
	else
		return false;
}

template <class Field, class RecursivePoly>
int strictlyLessThanRitt(const RecursivePoly& a, const RecursivePoly& b, const std::vector<Symbol>& vars) {
	if (a.isConstant()) {
		if (b.isConstant()) {
			return false;
		} else {
			return true;
		}
	}

	if (b.isConstant()) {
		return false;
	}

	int ai = variableIndex<Field, RecursivePoly>(a.leadingVariable(), vars);
	int bi = variableIndex<Field, RecursivePoly>(b.leadingVariable(), vars);
	if (ai < bi) {
		return false;
	}
	if (ai > bi) {
		return true;
	}

	//otherwise, ai==bi, and they have same mvar.
	ai = a.leadingVariableDegree().get_si();
	bi = b.leadingVariableDegree().get_si();
	if (ai < bi) {
		return true;
	}
	if (ai > bi) {
		return false;
	}

	//otherwise, they have equal degrees too
	SparseMultivariateRationalPolynomial aa = a.initial();
	SparseMultivariateRationalPolynomial bb = b.initial();
	if (strictlyLessThanRitt<Field, RecursivePoly>(aa, bb, vars)) {
		return true;
	}
	if (strictlyLessThanRitt<Field, RecursivePoly>(bb, aa, vars)) {
		return false;
	}

	//else, initials and ranks are same and we must recurse on tails
	aa = a.tail();
	bb = b.tail();
	return strictlyLessThanRitt<Field, RecursivePoly>(aa, bb, vars);
}

template <class Field, class RecursivePoly>
bool rankGreaterThan(const RecursivePoly& a, const RecursivePoly& b, const std::vector<Symbol>& vars) {
	// This algorithm assumes that a and b are the ranks of RecursivePoly elements, thus have a single term //
	//  and that both a and b are in the same ordered ring.                                                 //
	bool out;
	int ai,bi;
	ai = variableIndex<Field,RecursivePoly>(a.leadingVariable(),vars);
	bi = variableIndex<Field,RecursivePoly>(b.leadingVariable(),vars);
	if (ai < bi)
		return true;
	else if (ai == bi) {
		int aDeg,bDeg;
		aDeg = a.leadingVariableDegree().get_si();
		bDeg = b.leadingVariableDegree().get_si();
		if (aDeg > bDeg)
			return true;
		else
			return false;
	}
	else
		return false;
}

/**
 * Get the encoded options of the regular chain
 *
 **/
template <class Field, class RecursivePoly>
int RegularChain<Field,RecursivePoly>::options() const {
	return regularChainOptions;
}

/**
 * Set the encoded options of the regular chain
 *
 * @param opts: bitwise or of RegularChainOption values
 **/
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::setOptions(int opts) {
	if (opts == 0 || opts == MAINTAIN_SQUAREFREE) {
		regularChainOptions = opts;
	}
	else {
		std::cerr << "BPAS: error, invalid options for RegularChain class." << std::endl;
		exit(1);
	}
}

///**
// * Find out if the regular chain is known to be prime
// *
// * @param
// **/
//template <class Field, class RecursivePoly>
//bool RegularChain<Field,RecursivePoly>::isKnownToBeIrreducible() const {
//	return isIrreducible;
//}


/**
 * Find out if the input polynomial is in the saturated ideal of the current regular chain
 *
 * @param p: a polynomial
 **/
template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isInSaturatedIdeal(const RecursivePoly& p) const {
	RecursivePoly redp;
	return isInSaturatedIdeal(p,redp);
}

/**
 * Find out if the input polynomial is in the saturated ideal of the current regular chain
 * and return the reduced input polynomial
 *
 * @param p: a polynomial
 **/
template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isInSaturatedIdeal(const RecursivePoly& p, RecursivePoly& redp) const {

	// TODO: make this more efficient by not performing the full reduction unless necessary (make routine iteratively check polynomials
	//       in the chain and take into account whether p.leadingVariable is algebraic or not)
	RecursivePoly q;
	redp = this->reduce(p);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "T = " << *this << std::endl;
		std::cerr << "p = " << p << std::endl;
		std::cerr << "T.reduce(p) = " << redp << std::endl;
	#endif
	if (redp.isZero()) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "returning true" << std::endl;
		#endif
		return true;
	}
	else {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "returning false" << std::endl;
		#endif
		return false;
	}
}


/**
 * Find out if the input polynomial is in the saturated ideal of the current regular chain
 *
 * @param p: a polynomial
 **/
template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isInSaturatedIdealMinimal(const RecursivePoly& p) const {
	unsigned long long int inSatMinStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&inSatMinStart);
	#endif

	bool ret = isInSaturatedIdealMinimal_inner(p);

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&inSatMinStart, &isInSaturatedIdealTime);
	#endif

	return ret;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isInSaturatedIdealMinimal_inner(const RecursivePoly& p) const {

//	std::cerr << "there!" << *this << std::endl;
	// std::cerr << "Checking in sat minimal: " << std::endl;
	// std::cerr << "p  : " << p << std::endl;
	// std::cerr << "rc: " << *this << std::endl;

	if (p.isZero())
		return true;
	if (this->isConstantPolynomial(p) || this->isEmpty())
		return false;



	RegularChain<Field,RecursivePoly> Tlv;
	RecursivePoly q,r,Tv;
	Symbol v(p.leadingVariable());

//	std::cerr << "about to cut the chains!!" << std::endl;
	this->cutChain(v,Tlv,Tv);
//	std::cerr << "Tv variables:" << std::endl;
//	printVariables(Tv.ringVariables(),"Tvvars");

	if (!Tv.isZero()) {
		#if REGULARCHAIN_DEBUG
			std::cerr << "[isInSaturatedIdealMinimal] calling pdiv: " << std::endl;
			std::cerr << "p : " << p << std::endl;
			std::cerr << "Tv: " << Tv << std::endl;

		#endif
		q = p.primitivePart().pseudoDivide(Tv).primitivePart();
		//Jan29/2020: We need this right now, because sometimes Tv does not have all variables
		//of the ambient space and pdiv could re-order falsely the vars of the pseudo-remainder.
		q.setRingVariables(this->allVariables());
		if (q.isZero()) {
			return true;
		} else {
			//Jan27/2020: call isInMinimal here
			return Tlv.isInSaturatedIdealMinimal_inner(q);
		}
	}

	if (Tlv.isEmpty())
		return false;

	q = p.initial();
	r = p.tail();
	while (true) {
		if (!Tlv.isInSaturatedIdealMinimal_inner(q))
			return false;
		if (r.isZero())
			return true;
		if (this->isConstantPolynomial(r))
			return false;
		if (r.leadingVariable() == v) {
			q = r.initial();
			r = r.tail();
		}
		else {
			q = r;
			r.zero();
		}
	}

//	std::cerr << "end of isInSaturatedIdealMinimal!" << std::endl;

}

//template <class Field, class RecursivePoly>
//std::vector<int> sortByRank(const std::vector<RecursivePoly>& F, const std::vector<Symbol>& vars) {
//	std::vector<int> sortedIndices;
//	int i,j,x;
//	for (i=0; i<F.size(); ++i)
//		sortedIndices.push_back(i);
//
//	i = 1;
//	while (i < F.size()) {
//		x = sortedIndices[i];
//		j = i-1;
//		while (j >= 0 && rankGreaterThan<Field,RecursivePoly>(F[sortedIndices[j]],F[x],vars)) {
//		    sortedIndices[j+1] = sortedIndices[j];
//		    j = j-1;
//		}
//		sortedIndices[j+1] = x;
//		i = i+1;
//	}
//	return sortedIndices;
//}



template <class Field, class RecursivePoly>
std::vector<Symbol> findOrderedRing(const std::vector<RecursivePoly>& F) {
//	std::vector<Symbol> out,vars;
	if (F.empty())
		return std::vector<Symbol>();
	if (F.size() == 1)
		return F[0].ringVariables();
	else {
		RecursivePoly sum;
		for (size_t i=0; i<F.size(); ++i)
			sum += F[i];
		return sum.ringVariables();
	}
}

//template <class Field, class RecursivePoly>
//std::vector<Symbol> findOrderedRing(const std::vector<RecursivePoly>& F) {
//	if (F.empty())
//		return std::vector<Symbol>();
//	if (F.size() == 1)
//		return F[0].ringVariables();
//	else {
//		std::vector<Symbol> nzvs;
//		for (size_t i=0; i<F.size(); ++i) {
//			nzvs = setUnion(nzvs,F[i].variables()); // needs to replicate ordering of summation
//		}
//		return nzvs;
////		RecursivePoly sum;
////		for (size_t i=0; i<F.size(); ++i)
////			sum += F[i];
////		return sum.ringVariables();
//	}
//}


template <class Field, class RecursivePoly>
std::vector<RecursivePoly> sortDecreasingByRank(const std::vector<RecursivePoly>& F) {
	std::vector<RecursivePoly> sortedF(F);
	std::vector<Symbol> ringVars;
	RecursivePoly x;
	int i,j;

	// TODO: Do I need to check consistency with the ordering of the rc, or should cleanSet do this?
	//  and should behaviour depend on the mode of the triangular set?
	ringVars = findOrderedRing<Field,RecursivePoly>(F);

	i = 1;
	while (i < F.size()) {
		x = sortedF[i];
		j = i-1;
		while (j >= 0 && rankLessThan<Field,RecursivePoly>(sortedF[j],x,ringVars)) {
		    sortedF[j+1] = sortedF[j];
		    j = j-1;
		}
		sortedF[j+1] = x;
		i = i+1;
	}
	return sortedF;
}

template <class Field, class RecursivePoly>
std::vector<RecursivePoly> sortIncreasingByRank(const std::vector<RecursivePoly>& F) {
	std::vector<RecursivePoly> sortedF(F);
	std::vector<Symbol> ringVars;
	RecursivePoly x;
	int i,j;

	// TODO: Do I need to check consistency with the ordering of the rc, or should cleanSet do this?
	//  and should behaviour depend on the mode of the triangular set?
	ringVars = findOrderedRing<Field,RecursivePoly>(F);

	i = 1;
	while (i < F.size()) {
		x = sortedF[i];
		j = i-1;
		while (j >= 0 && rankGreaterThan<Field,RecursivePoly>(sortedF[j],x,ringVars)) {
		    sortedF[j+1] = sortedF[j];
		    j = j-1;
		}
		sortedF[j+1] = x;
		i = i+1;
	}
	return sortedF;
}

template <class Field, class RecursivePoly>
void sortDecreasingStrictRitt(std::vector<RecursivePoly>& F, std::vector<Symbol> ringVars) {
	// std::vector<Symbol> ringVars;
	RecursivePoly x;
	int i,j;

	// ringVars = findOrderedRing<Field,RecursivePoly>(F);


	i = 1;
	while (i < F.size()) {
		x = F[i];
		j = i-1;
		while (j >= 0 && strictlyLessThanRitt<Field,RecursivePoly>(F[j],x,ringVars)) {
		    F[j+1] = F[j];
		    j = j-1;
		}
		F[j+1] = x;
		i = i+1;
	}
}

//this is called TRDreduce_primpart in RegularChains.
template <class Field, class RecursivePoly>
RecursivePoly RegularChain<Field,RecursivePoly>::reduceMinimalPrimitivePart(const RecursivePoly& p_in) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RM] entering reduceMinimalPrimitivePart: " << p_in << std::endl;
		std::cerr << "[RM] *this = " << *this << std::endl;
	#endif

	RecursivePoly p;
	p = reduceMinimal(p_in);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RM] leaving reduceMinimalPrimitivePart input: " << p_in << std::endl;
		std::cerr << "[RM] leaving reduceMinimalPrimitivePart output: " << p << std::endl;
	#endif

	return p.mainPrimitivePart();

}

//this is called TRDreduce_primpart_sqf in RegularChains.
template <class Field, class RecursivePoly>
RecursivePoly RegularChain<Field,RecursivePoly>::reduceMinimalPrimitiveSquareFreePart(const RecursivePoly& p_in) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RM] entering reduceMinimalPrimitiveSquareFreePart: " << p_in << std::endl;
		std::cerr << "[RM] *this = " << *this << std::endl;
	#endif

	RecursivePoly p;
	p = reduceMinimalPrimitivePart(p_in);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RM] leaving reduceMinimalPrimitiveSquareFreePart input: " << p_in << std::endl;
		std::cerr << "[RM] leaving reduceMinimalPrimitiveSquareFreePart output: " << p << std::endl;
	#endif

	return p.squareFreePart();

}

//this is called TRDreduce in RegularChains.
template <class Field, class RecursivePoly>
RecursivePoly RegularChain<Field,RecursivePoly>::reduceMinimal(const RecursivePoly& p) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RM] entering reduceMinimal." << std::endl;
		std::cerr << "[RM] p = " << p << std::endl;
		std::cerr << "[RM] *this = " << *this << std::endl;
	#endif
	if (this->isConstantPolynomial(p) || this->isEmpty()) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[RM] leaving reduceMinimal, trivial case: " << p << std::endl;
		#endif
		return p;
	}

	// TODO: when doing modular computation, a version of TRDmodulo is needed here
	unsigned long long int reduceMinStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&reduceMinStart);
	#endif
	SMQP ret = removeZero(moduloPolysWithConstantInitial(p)); // TODO: compute moduloPolysWithConstantInitial here first
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&reduceMinStart, &reduceMinimalTime);
	#endif
	return  ret;

}

template <class Field, class RecursivePoly>
RecursivePoly RegularChain<Field,RecursivePoly>::removeZero(const RecursivePoly& p) const {
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RZ] entering removeZero: " << p << std::endl;
		printVariables(this->variables(),"thisvars");
	#endif
	if (this->isConstantPolynomial(p) || this->isEmpty())
		return p;

	std::vector<Symbol> ringVars, subRingVars;
//	printVariables(ringVars,"ringVars");
//	std::vector<RecursivePoly> polys(this->polynomials());
//	polys.push_back(p);
//	ringVars = findOrderedRing<Field,RecursivePoly>(polys);
	if (mode == TS_FIXED) {
		ringVars = this->allVariables();
	} else {
		ringVars = orderPreservingSetUnion(p.ringVariables(), this->allVariables());
	}

	RegularChain<Field,RecursivePoly> Tlv;
	RecursivePoly q(p),r,s,Tv;
	Symbol v(p.leadingVariable());
	std::vector<Symbol> vars = {v};

	this->cutChain(v,Tlv,Tv);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "p = " << p << std::endl;
		std::cerr << "T = " << *this << std::endl;
		std::cerr << "v = " << v << std::endl;
		std::cerr << "Tv = " << Tv << std::endl;
		std::cerr << "Tlv = " << Tlv << std::endl;
	#endif

	if (!Tv.isZero()) {
		if (Tv.leadingVariableDegree() == 1 && Tv.tail().isZero()) {
			std::vector<Field> zero = {Field()};
			q = q.evaluate(vars,zero);
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[RZ] making recursive call, position 1." << std::endl;
			#endif
			r = Tlv.removeZero(q);
			assert(this->isInSaturatedIdealMinimal(r-p) && "First assert condition failed in removeZero");
//			subRingVars = findOrderedSubRing<Field,RecursivePoly>(r.variables(),ringVars);
			r.setRingVariables(ringVars);
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[RZ] leaving removeZero position 1: " << r << std::endl;
			#endif
			return r;
		}
		// subRingVars = findOrderedSubRing<Field,RecursivePoly>(Tv.variables(),q.variables(),ringVars);
		#ifdef REGULARCHAIN_DEBUG
			printVariables(ringVars,"ringVars");
		#endif
		Tv.setRingVariables(ringVars);
		q.setRingVariables(ringVars);
		s = q.pseudoDivide(Tv);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[RZ] q.pseudoDivide(Tv) = " << s << std::endl;
			std::cerr << "[RZ] q = " << q << std::endl;
			std::cerr << "[RZ] Tv = " << Tv << std::endl;
			printVariables(q.ringVariables(),"q");
			printVariables(Tv.ringVariables(),"Tv");
			printVariables(s.ringVariables(),"s");
		#endif
//		std::cerr << "hello!" << s << std::endl;
		// printVariables(s.ringVariables(),"svars");

		if (Tlv.isInSaturatedIdealMinimal(s)) {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[RZ] s in Sat(T)" << std::endl;
			#endif
			r.zero();
			//Jan 29/2020 debugging
			// bool isInSatFull = this->isInSaturatedIdeal(r-p);
			// bool isInSatMin = this->isInSaturatedIdealMinimal(r-p);
			// std::cerr << "is in full : " << isInSatFull << std::endl;
			// std::cerr << "is in min : " << isInSatMin << std::endl;

			assert(this->isInSaturatedIdealMinimal(r-p) && "Second assert condition failed in removeZero");
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[RZ] leaving removeZero position 2: " << r << std::endl;
			#endif
			return r;
		}
	}

	if (Tlv.isEmpty()) {
//		r = q;
//		std::cerr << "q = " << q << std::endl;
//		printVariables(q.ringVariables(),"qvars");
//		printVariables(p.ringVariables(),"pvars");
		assert(this->isInSaturatedIdealMinimal(q-p) && "Third assert condition failed in removeZero");
//		subRingVars = findOrderedSubRing<Field,RecursivePoly>(q.variables(),ringVars);
		q.setRingVariables(ringVars);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[RZ] leaving removeZero position 3: " << q << std::endl;
		#endif
		return q;
	}

	s.zero();
	while (q.degree(v) > 0) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[RZ] making recursive call, position 2." << std::endl;
		#endif
		r = Tlv.removeZero(q.initial()) * q.rank();
		subRingVars = findOrderedSubRing<Field,RecursivePoly>(r.variables(),s.variables(),ringVars);
		r.setRingVariables(subRingVars);
		s.setRingVariables(subRingVars);
		s += r;
		q = q.tail();
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RZ] making recursive call, position 3." << std::endl;
	#endif
	r = Tlv.removeZero(q);
	subRingVars = findOrderedSubRing<Field,RecursivePoly>(r.variables(),s.variables(),ringVars);
	r.setRingVariables(subRingVars);
	s.setRingVariables(subRingVars);
	r += s;
	assert(this->isInSaturatedIdealMinimal(r-p) && "Fourth assert condition failed in removeZero");
//	subRingVars = findOrderedSubRing<Field,RecursivePoly>(r.variables(),ringVars);
	r.setRingVariables(ringVars);
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[RZ] leaving removeZero at end: " << r << std::endl;
	#endif
	return r;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::cleanSet(std::vector<RecursivePoly>& polys) const {
	unsigned long long int cleanSetStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&cleanSetStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "entering cleanSet." << std::endl;
	#endif
	std::vector<RecursivePoly> cleanPolys;
	RecursivePoly p;
	for (size_t i=0; i<polys.size(); ++i) {
		p = polys[i];
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "calling removeZero on p: " << p << std::endl;
//		#endif
		// printVariables(p.ringVariables(), "p b4 remZero ");
		p = this->removeZero(p);
		// printVariables(p.ringVariables(), "p b4 modPolyConst ");
		p = this->moduloPolysWithConstantInitial(p); //Jan16/2020: Maple does removeZero then modPolys, but maybe do opposite? as what we do in reduceMinimal
		// printVariables(p.ringVariables(), "p b4 primPart ");
		p = p.primitivePart();
		// printVariables(p.ringVariables(), "p at end ");
		// we might need to modularizeCoefficients here if we have a modular option, but this should be handled by the polynomial ring
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "cleanSet: p after reduction = " << p << std::endl;
		#endif
		if (!p.isZero()) {
			if (this->isConstantPolynomial(p)) {
				polys.clear();
				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&cleanSetStart,&cleanSetTime);
				#endif
				return true;
			}
			else {
				cleanPolys.emplace_back(std::move(p));
			}
		}
	}
	polys = cleanPolys;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "leaving cleanSet." << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&cleanSetStart,&cleanSetTime);
	#endif
	return false;
}

//template <class Field, class RecursivePoly>
//bool RegularChain<Field,RecursivePoly>::cleanSet(const std::vector<RecursivePoly>& inPolys, std::vector<RecursivePoly>& outPolys) const {
//	std::vector<RecursivePoly> cleanPolys;
//	RecursivePoly p;
//	for (size_t i=0; i<inPolys.size(); ++i) {
//		p = inPolys[i];
////		p = this->reduce(p);
//		p = this->removeZero(p);
//		p = p.primitivePart();
//		// we might need to modularizeCoefficients here if we have a modular option, but this should be handled by the polynomial ring
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "cleanSet: p after reduction = " << p << std::endl;
//		#endif
//		if (!p.isZero()) {
//			if (this->isConstantPolynomial(p)) {
//				outPolys.clear();
//				return true;
//			}
//			else {
//				outPolys.push_back(p);
//			}
//		}
//	}
//	return false;
//}


#if RC_TRIANGULARIZE_TREEDATA
int currentLevel;
unsigned long long plotTriTime;
SynchronizedWriteVector<std::tuple<float, int, RegularChain<RationalNumber, SparseMultivariateRationalPolynomial>, int>> timeStamps;
#endif

template <class Field, class RecursivePoly>
void intersectOne(int j, const RegularChain<Field, RecursivePoly>& T, const RecursivePoly& p, int lazardDecompose, int heightBound,
	SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>& nextResults) {
	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;

#ifdef REGULARCHAIN_PROFILING
	unsigned long long start = 0;
	float elapsedTime = 0;
	_startTimer(&start);
#endif

	RecursivePoly q = T.removeZero(p);

	if (q.isZero()) {
		nextResults.push_back(T);
#if RC_TRIANGULARIZE_TREEDATA
		float curTime = 0;
		stopTimerAddElapsed(&plotTriTime, &curTime);
		timeStamps.push_back(std::tuple<float, int, RC_INT_OBJ, int>(curTime, currentLevel, T, j));
#endif
		// std::cerr << "[tri] [" << j << "] qIsZero = " << newResults[j] << endl;
	}
	else if (!T.isConstantPolynomial(q)) { // if p is not inconsistent with T

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[tri] T = " << T << std::endl;
			printVariables(T.variables());
			std::cerr << "[tri] triangularize: calling intersect..." << q << std::endl;
		#endif

		#ifdef REGULARCHAIN_PROFILING
			// stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
		#endif
		RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(T), q, lazardDecompose, heightBound);
		#ifdef REGULARCHAIN_PROFILING
			// startTimer(&triangularizeStart);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[tri] returning to triangularize from intersect..." << q << std::endl;
		#endif

		RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);
		RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, k) {
			RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, k);
			#ifdef REGULARCHAIN_DEBUG
				// std::cerr << "[tri] [" << j << "] currIntTask[" << k << "] = " << currIntTask << endl;
			#endif
#if RC_TRIANGULARIZE_TREEDATA
			float curTime = 0;
			stopTimerAddElapsed(&plotTriTime, &curTime);
			timeStamps.push_back(std::tuple<float, int, RC_INT_OBJ, int>(curTime, currentLevel, currIntTask, j));
#endif
			nextResults.push_back(std::move(currIntTask));
		}
	}
#ifdef REGULARCHAIN_PROFILING
	_stopTimerAddElapsed(&start, &elapsedTime);
#endif
}

#if defined(SERIAL) && SERIAL
#define TRIANGULARIZE_SERIAL 1
#else
#define TRIANGULARIZE_SERIAL 0
#endif

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::triangularize(const std::vector<RecursivePoly>& polys, bool lazardDecompose, int type) {
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::triangularize(const std::vector<RecursivePoly>& polys) {
	long long unsigned int rcProfilingStart,triangularizeStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&triangularizeStart);

		primitivePartTime = 0;
		squareFreePartTime = 0;
		subresultantChainTime = 0;
		zerodimensionalregularchainTime = 0;
		pseudoDivideTime = 0;
		normalFormTime = 0;
		removeRedundantChainsTime = 0;
		factorTime = 0;
		triangularizeTime = 0;
		intersectTime = 0;
		regularGCDTime = 0;
		intersectFreeTime = 0;
		intersectAlgebraicTime = 0;
		regularizeTime = 0;
		extendTime = 0;
		cleanChainTime = 0;
		constructChainTime = 0;
		constructChainsTime = 0;
		GCDFreeFactorizationTime = 0;
		ZDIntersectTime = 0;
		ZDRegularizeTime = 0;
		isInSaturatedIdealTime = 0;
		cleanSetTime = 0;
		tsCopyTime = 0;
		rcCopyTime = 0;
		zdrcCopyTime = 0;
	#endif
	intersectDepth = 0;
	regularGCDDepth = 0;
	intersectFreeDepth = 0;
	intersectAlgebraicDepth = 0;
	regularizeDepth = 0;
	extendDepth = 0;
	cleanChainDepth = 0;
	depth = 0;
	totalIntersects = 0;
	totalSplits = 0;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[tri] entering triangularize..." << std::endl;
		std::cerr << "[tri] regularChainOptions = " << regularChainOptions << std::endl;
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;

	typedef RC_INT_OBJ RC_TRI_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_TRI_OBJ, currTriTask);
	std::vector<RegularChain<Field,RecursivePoly>> results;
//	std::vector<RegularChain<Field,RecursivePoly>> results,newResults,nextResults;
	SynchronizedWriteVector<RegularChain<Field,RecursivePoly>> nextResults;

	if (vars.empty()) {
		vars = findOrderedRing<Field,RecursivePoly>(polys);
		set.reserve(vars.size());
		for (size_t i=0; i<vars.size(); ++i) {
			set.emplace_back(); // initialize elements of set with empty polynomials
		}
		mode = TS_FIXED;
	}

	if (polys.empty()) {
		results.emplace_back(*this);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
		#endif
		return results;
	}

	struct Task {
		std::vector<RecursivePoly> polys;
		RegularChain<Field,RecursivePoly> chain;
	};
	Task currTask,newTask;
	std::vector<Task> tasks,newTasks;
	std::vector<RecursivePoly> F(polys),reducedPolys;
	RecursivePoly p;
	bool inconsistent;
	int size;

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
	#endif
	inconsistent = cleanSet(F);
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&triangularizeStart);
	#endif

	#ifdef REGULARCHAIN_DEBUG
		for (size_t i=0; i<F.size(); ++i)
			std::cerr << "[tri] F[" << i << "] = " << F[i] << std::endl;
		std::cerr << "[tri] inconsistent = " << inconsistent << std::endl;
	#endif

	if (inconsistent) {
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
		#endif
		return results;
	}

	// F = sortDecreasingByRank<Field,RecursivePoly>(F);
	sortDecreasingStrictRitt<Field, RecursivePoly>(F, this->allVariables());

	#ifdef REGULARCHAIN_DEBUG
		for (size_t i=0; i<F.size(); ++i)
			std::cerr << "[tri] F[" << i << "] = " << F[i] << std::endl;
	#endif

	int heightBound;
	if (lazardDecompose) {
		heightBound = this->numberOfVariables();
	} else {
		heightBound = F.size();
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[tri] heightBound = " << heightBound << std::endl;
	#endif

	if (type == 0) { // Bubble Variety of Triangularize
		//old bubble
		// std::vector<RegularChain<Field,RecursivePoly>> G;
		// G.push_back(*this);
		// RC_TRI_CONSUMER_INIT(RC_TRI_OBJ, triTasks, (&RegularChain<Field, RecursivePoly>::_triangularize), this, F, G, lazardDecompose, heightBound);
		// RC_TRI_CONSUMER_LOOP(triTasks, currTriTask, i) {
		// 	RC_TRI_CONSUMER_GET_LOOPELEM(triTasks, currTriTask, i);
		// 	results.push_back(currTriTask);
		// }

		//non-recursive, task-based, cleans list of polynomials waiting for each new component
		results = this->_triangularizeByTasks(F, lazardDecompose, heightBound);


//		results = _triangularize(F,G,lazardDecompose,heightBound);


//		currTask = {F,*this};
//		tasks.push_back(currTask);
//		while (tasks.size() > 0) {
//			currTask = tasks.back();
//			tasks.pop_back();
//			if (currTask.polys.empty()) {
//				results.emplace_back(std::move(currTask.chain));
//			}
//			else {
//				currTask.polys = sortDecreasingByRank<Field,RecursivePoly>(currTask.polys);
//				p = currTask.polys.back();
//				currTask.polys.pop_back();
//				#ifdef REGULARCHAIN_DEBUG
//					std::cerr << "[tri] currTaks.chain = " << currTask.chain << std::endl;
//					printVariables(currTask.chain.variables());
//					std::cerr << "[tri] triangularize: calling intersect..." << p << std::endl;
//				#endif
//
//				#ifdef REGULARCHAIN_PROFILING
//					stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
//				#endif
//				RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(currTask.chain), p, lazardDecompose, heightBound);
////				newResults = currTask.chain.intersect(p,lazardDecompose,heightBound);
//				#ifdef REGULARCHAIN_PROFILING
//					startTimer(&triangularizeStart);
//				#endif
//
//				#ifdef REGULARCHAIN_DEBUG
//					std::cerr << "[tri] returning to triangularize from intersect..." << p << std::endl;
//				#endif
//
//				// TODO: Minimize copies by reducing container use (currTask,newTask,currTaskPolys,reducedPolys)
//
//				RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, j) {
////				for (size_t j=0; j<newResults.size(); ++j) {
//					RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
//					reducedPolys = currTask.polys;
//
//					#ifdef REGULARCHAIN_DEBUG
//						std::cerr << "[tri] newResults[" << j << "] = " << currIntTask << endl;
//						std::cerr << "[tri] reducedPolys before cleanSet:" << std::endl;
//						for (size_t i=0; i<reducedPolys.size(); ++i) {
//							std::cerr << "[tri] reducedPolys[" << i << "] = " << reducedPolys[i] << std::endl;
//						}
//					#endif
//
//					#ifdef REGULARCHAIN_PROFILING
//						stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
//					#endif
//					inconsistent = currIntTask.cleanSet(reducedPolys);
////					inconsistent = newResults[j].cleanSet(reducedPolys);
//					#ifdef REGULARCHAIN_PROFILING
//						startTimer(&triangularizeStart);
//					#endif
//
//					#ifdef REGULARCHAIN_DEBUG
//						std::cerr << "[tri] inconsistent = " << inconsistent << std::endl;
//					#endif
//					if (!inconsistent) {
//						newTask = {reducedPolys,currIntTask};
////						newTask = {reducedPolys,newResults[j]};
//						tasks.push_back(newTask);
//					}
//					#ifdef REGULARCHAIN_DEBUG
//						std::cerr << "[tri] reducedPolys after cleanSet:" << std::endl;
//						for (size_t i=0; i<reducedPolys.size(); ++i) {
//							std::cerr << "[tri] reducedPolys[" << i << "] = " << reducedPolys[i] << std::endl;
//						}
//						std::cerr << "[tri] currTaks.chain = " << currTask.chain << std::endl;
//					#endif
//				}
//			}
//		}
	}
	else { // Level Variety of Triangularize
		results.push_back(*this);
#if RC_TRIANGULARIZE_TREEDATA
		startTimer(&plotTriTime);
		fprintf(stdout, "(%f, %d, %d, %d, %d)\n", 0.0f, 0, 0, -1, this->dimension());
		currentLevel = 1;
#endif
		for (size_t i=F.size(); i-- > 0; ) {
		// while (F.size() > 0) {
			// bool inconsistent = currIntTask.cleanSet(curPolys);
			// sortDecreasingStrictRitt<Field,RecursivePoly>(curTask.first, this->allVariables());

			++totalIntersects;
			totalSplits += results.size();
			p = F[i];
			//inconsistent system

			if (results.size() == 0) {
				results.clear();
				break;
			}

			if (results.size() == 1) {
				// fprintf(stderr, "Manual call with 1 component\n");
				intersectOne(0, results[0], p, lazardDecompose, heightBound, nextResults);
			} else {

#if defined(TRIANGULARIZE_SERIAL) && TRIANGULARIZE_SERIAL
#else
			std::vector<std::thread> thrArr;
#endif
			size_t j;
			for (j=0; j<results.size()-1; ++j) {
#if defined(TRIANGULARIZE_SERIAL) && TRIANGULARIZE_SERIAL
				intersectOne(j, results[j], p, lazardDecompose, heightBound, nextResults);
#else
				std::thread t1(intersectOne<Field,RecursivePoly>, j, std::ref(results[j]), std::ref(p), (int)lazardDecompose, heightBound, std::ref(nextResults));
				thrArr.push_back(std::move(t1));
#endif
			}
			intersectOne(j, results[j], p, lazardDecompose, heightBound, nextResults);
#if defined(TRIANGULARIZE_SERIAL) && TRIANGULARIZE_SERIAL
#else
			for (std::thread& t : thrArr) {
				t.join();
			}
#endif
			} //is newresults > 1

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[tri] results before removeRedundantChain results:" << std::endl;
				for (size_t k=0; k<nextResults.size(); ++k) {
					std::cerr << "[tri] results[" << k << "] = " << nextResults[k] << std::endl;
//					printVariables(nextResults[k].variables(),"nRVars");
				}
			#endif

#if RC_TRIANGULARIZE_TREEDATA
			int parentSize = results.size();
			int beforeSize = nextResults.size();
#endif
			results.clear();
			// RegularChain<Field,RecursivePoly>::removeRedundantChains(nextResults.vector(), results);
			results = nextResults.vector();

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[tri] results after removeRedundantChains call:" << std::endl;
				for (size_t k=0; k<results.size(); ++k) {
					std::cerr << "[tri] results[" << k << "] = " << results[k] << std::endl;
//					printVariables(nextResults[k].variables(),"nRVars");
				}
			#endif

#if RC_TRIANGULARIZE_TREEDATA

			if (beforeSize < results.size()) {
				// fprintf(stderr, "Remove redundant split\n");

				// std::cerr << "Before remove:\n";
				// for (auto t : nextResults.vector()) {
				// 	std::cerr << t << std::endl;
				// }

				// std::cerr << "After remove:\n";
				// for (auto t : results) {
				// 	std::cerr << t << std::endl;
				// }

				exit(1);

			}

			float curTime = 0;
			stopTimerAddElapsed(&plotTriTime, &curTime);
			fprintf(stdout, "level %d timestamp %f\n", currentLevel, curTime);
			if (results.size() == 0) {
				fprintf(stdout, "(%f, %d, %d, %d, %d)\n", curTime, currentLevel, 0, -1, -1);
			}

			for (int j = 0; j < results.size(); ++j) {
				for (size_t i = 0; i < timeStamps.size(); ++i) {
					if (std::get<2>(timeStamps[i]) == results[j]) {
						fprintf(stdout, "(%f, %d, %d, %d, %d)\n", std::get<0>(timeStamps[i]), std::get<1>(timeStamps[i]), j, std::get<3>(timeStamps[i]), results[j].dimension());
						break;
					}
				}
			}

			timeStamps.clear();
			++currentLevel;
#endif

			nextResults.clear();
		}
	}
//	else if (type == 1) {
////		newResults.push_back(*this);
//		results.emplace_back(*this);
//		for (int i=F.size()-1; i>=0; --i) {
//			p = F[i];
//			fprintf(stderr, "Starting loop\n");
//			for (size_t j=0; j<results.size(); ++j) {
////			for (size_t j=0; j<newResults.size(); ++j) {
////				const RegularChain<Field,RecursivePoly>& T = newResults[j];
//				const RegularChain<Field,RecursivePoly>& T = results[j];
//
//				q = T.removeZero(p);
//
//				if (q.isZero()) {
////					nextResults.push_back(newResults[j]);
//					nextResults.emplace_back(std::move(results[j]));
//				}
//				else if (!T.isConstantPolynomial(q)) { // if p is not inconsistent with T
//
//					#ifdef REGULARCHAIN_DEBUG
//						std::cerr << "[tri] T = " << T << std::endl;
//						printVariables(T.variables());
//						std::cerr << "[tri] triangularize: calling intersect..." << q << std::endl;
//					#endif
//
//					#ifdef REGULARCHAIN_PROFILING
//						stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
//					#endif
//					RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(T), q, lazardDecompose, heightBound);
////					newResults = T.intersect(q,lazardDecompose,heightBound);
//					#ifdef REGULARCHAIN_PROFILING
//						startTimer(&triangularizeStart);
//					#endif
//
//					#ifdef REGULARCHAIN_DEBUG
//						std::cerr << "[tri] returning to triangularize from intersect..." << q << std::endl;
//					#endif
//
//					RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, k) {
//						RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, k);
////					for (size_t k=0; k<newResults.size(); ++k) {
////						currIntTask = newResults[k];
//						#ifdef REGULARCHAIN_DEBUG
//							std::cerr << "[tri] currIntTask[" << k << "] = " << currIntTask << endl;
//						#endif
//						nextResults.push_back(std::move(currIntTask));
////						nextResults.emplace_back(std::move(currIntTask));
//					}
////					nextResults.insert(nextResults.end(),newResults.begin(),newResults.end());
//
//				}
//			}
//			fprintf(stderr, "finished loop\n");
//
//			#ifdef REGULARCHAIN_DEBUG
////				std::cerr << "[tri] results before removeRedundantChain results:" << std::endl;
////				for (size_t k=0; k<nextResults.size(); ++k) {
////					std::cerr << "[tri] results[" << k << "] = " << nextResults[k] << std::endl;
////				}
//			#endif
////			results.clear();
////			for (size_t k=0; k<nextResults.size(); ++k) {
////				results.push_back(nextResults[k]);
////			}
////			nextResults.move_vector(results);
////			results = RegularChain<Field,RecursivePoly>::removeRedundantChains(results);
//			results = RegularChain<Field,RecursivePoly>::removeRedundantChains(nextResults);
//			#ifdef REGULARCHAIN_DEBUG
//				std::cerr << "[tri] results after removeRedundantChains call:" << std::endl;
//				for (size_t k=0; k<results.size(); ++k) {
//					std::cerr << "[tri] results[" << k << "] = " << results[k] << std::endl;
//				}
//			#endif
//
//			newResults.clear();
////			newResults.reserve(nn);
////			fprintf(stderr, "packing new results");
////			for (size_t k=0; k<results.size(); ++k) {
////				newResults.push_back(results[k]);
////			}
////			fprintf(stderr,"finished packing new results");
////			nextResults.clear();
////			fprintf(stderr, "clearing next results\n");
////			nextResults.reserve(nn);
//
//		}
//	}

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&triangularizeStart,&triangularizeTime);
	#endif

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[tri] leaving triangularize" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		std::fstream fs;
		fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
		fs << "prime time: " << primitivePartTime << "\t";
		fs << squareFreePartTime << "\t";
		fs << factorTime << "\t";
		fs << pseudoDivideTime << "\t";
		fs << normalFormTime << "\t";
		fs << subresultantChainTime << "\t";
		fs << removeRedundantChainsTime << "\t";
		fs << zerodimensionalregularchainTime << "\t";
		fs << "reduce time: " << reduceMinimalTime << "\t";
		fs << "sat ideal time: " << isInSaturatedIdealTime << "\t";
		fs << "SMQPSRC time: " << SMQP_SRC_TIME << "\t";
		fs.close();
		fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
		fs << intersectTime << "\t";
		fs << regularizeTime << "\t";
		fs << regularGCDTime << "\t";
		fs << intersectFreeTime << "\t";
		fs << intersectAlgebraicTime << "\t";
		fs << extendTime << "\t";
		fs << cleanChainTime << "\t";
		fs << cleanSetTime << "\t";
		fs << GCDFreeFactorizationTime << "\t";
		fs << triangularizeTime << "\t";
		fs << constructChainTime << "\t";
		fs << constructChainsTime << "\t";
		fs << (intersectTime + regularizeTime + regularGCDTime + intersectFreeTime + intersectAlgebraicTime + extendTime + cleanChainTime + cleanSetTime + GCDFreeFactorizationTime + triangularizeTime + constructChainTime + constructChainsTime) << "\t";
		fs.close();
		fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
		fs << tsCopyTime << "\t";
		fs << rcCopyTime << "\t";
		fs << zdrcCopyTime << "\t";
		fs << (intersectTime + regularizeTime + regularGCDTime + intersectFreeTime + intersectAlgebraicTime + extendTime + cleanChainTime + cleanSetTime + GCDFreeFactorizationTime + triangularizeTime + constructChainTime + constructChainsTime) << "\t";
		fs.close();

	#endif

	// fprintf(stderr, "\n\nTotal Intersects: %d\n\n", totalIntersects);
	// fprintf(stderr, "\n\nTotal Splits: %d\n\n", totalSplits);

	if (type == 0) {
		std::vector<RegularChain<Field,RecursivePoly>> results2;
		#if REGULARCHAIN_DEBUG
		for (size_t i=0; i<results.size(); ++i) {
			std::cerr << "results[" << i << "] = " << results[i] << std::endl;
		}
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);

			removeRedundantChains(results,results2);
//			std::vector<RegularChain<Field,RecursivePoly>> results2(removeRedundantChains(results));

			stopTimerAddElapsed(&rcProfilingStart,&removeRedundantChainsTime);

			return results2;
		#else
			RegularChain<Field,RecursivePoly>::removeRedundantChains(results,results2);

			// #if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
			std::stringstream ss;
			for (auto r : results2) {
				int degProd = 1;
				for (auto p : r.polynomials()) {
					if (!p.isZero()) {
						degProd *= p.mainDegree();
					}
				}
				ss << " (" << r.dimension() << "," << degProd << ") ";
			}
			fprintf(stderr, "Component Info:%s\n", ss.str().c_str());
			// #endif

			return results2;
//			return RegularChain<Field,RecursivePoly>::removeRedundantChains(results);
		#endif
	}
	else {
		#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
		for (auto r : results) {
			std::cerr << "Final Dimension: " << r.dimension() << std::endl;
		}
		#endif
		return results;
	}
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_triangularize(const std::vector<RecursivePoly>& polys, const std::vector<RegularChain<Field,RecursivePoly>>& chains, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::_triangularize(const std::vector<RecursivePoly>& polys, std::vector<RegularChain<Field,RecursivePoly>>& chains, bool lazardDecompose, int heightBound) const {
#endif

	typedef RegularChain<Field,RecursivePoly> RC_TRI_OBJ;
	RC_TRI_PRODUCER_DECLARE_RESULT(RC_TRI_OBJ, results);
	RC_TRI_CONSUMER_DECLARE_LOOPELEM(RC_TRI_OBJ, currTriTask);
	typedef RC_TRI_OBJ RC_INT_OBJ;
	RC_TRI_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);

	std::vector<RecursivePoly> F(polys);
	std::vector<RegularChain<Field,RecursivePoly>> tasks,newTasks;
//	std::vector<RegularChain<Field,RecursivePoly>> tasks,newTasks,results;
	RecursivePoly p,q;
	if (polys.size() == 0) {
		for (auto chain : chains) {
			RC_TRI_PRODUCER_ACCUMULATE(results, chain);
		}
		RC_TRI_PRODUCER_COMPLETE(results);
//		return chains;
	}

	// F = sortDecreasingByRank<Field,RecursivePoly>(F);
	sortDecreasingStrictRitt<Field,RecursivePoly>(F, this->allVariables());

	p = F[0];
	F.erase(F.begin());
	RC_TRI_CONSUMER_INIT(RC_TRI_OBJ, triTasks, (&RegularChain<Field, RecursivePoly>::_triangularize), this, F, chains, lazardDecompose, heightBound);
//	tasks = _triangularize(F,chains,lazardDecompose,heightBound);

	RC_TRI_CONSUMER_LOOP(triTasks, currTriTask, i) {
//	for (size_t i=0; i<tasks.size(); ++i) {
		RC_TRI_CONSUMER_GET_LOOPELEM(triTasks, currTriTask, i);
		q = currTriTask.removeZero(p);
//		q = tasks[i].removeZero(p);
		RC_TRI_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(currTriTask), q, lazardDecompose, heightBound);
//		newTasks = currTriTask.intersect(q,lazardDecompose,heightBound);
//		newTasks = tasks[i].intersect(q,lazardDecompose,heightBound);
		RC_TRI_CONSUMER_LOOP(intTasks, currIntTask, j) {
			RC_TRI_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
//		for (size_t j=0; j<newTasks.size(); ++j) {
			RC_TRI_PRODUCER_ACCUMULATE(results, currIntTask);
//			RC_GEN_PRODUCER_ACCUMULATE(results, newTasks[j]);
//			results.push_back(newTasks[j]);
		}
	}
	RC_TRI_PRODUCER_COMPLETE(results);
//	return results;
}

#if RC_TRIANGULARIZE_TASKTREEDATA
#include <atomic>
static std::atomic<int> RC_UNIQUE_ID;
unsigned long long plotTriTime;
#endif

template <class Field, class RecursivePoly>
void triangularizeTask(const RegularChain<Field, RecursivePoly>& rc, std::vector<RecursivePoly>& polys, bool lazardDecompose, int heightBound, TaskScheduler* tasks, SyncVector_ptr<RegularChain<Field,RecursivePoly>> results) {


#if defined(RC_TRIANGULARIZE_PROGRESS) && RC_TRIANGULARIZE_PROGRESS
	fprintf(stderr, "[RegChain] Starting Task: [%lu, %d], Current Results: [", polys.size(), rc.dimension());
	bool first = true;
	for (int i = 0; i < results->size(); ++i) {
		if (first) {
			fprintf(stderr, "%d", (*results)[i].dimension());
			first = false;
		} else {
			fprintf(stderr, ", %d", (*results)[i].dimension());
		}
	}
	fprintf(stderr, "]\n");
#endif

	sortDecreasingStrictRitt<Field,RecursivePoly>(polys, rc.allVariables());
	RecursivePoly p = polys.back();
	polys.pop_back();
	std::vector<RecursivePoly> cleanedPolys;

	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);

	int noComponents = 1;
	RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(rc), p, lazardDecompose, heightBound);
	RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, i);

		noComponents = 0;
		cleanedPolys = polys;
		bool inconsistent = currIntTask.cleanSet(cleanedPolys);

#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
		currIntTask.RegChain_UniqueID = RC_UNIQUE_ID.fetch_add(1);
		float curTime = 0.0f;
		_stopTimerAddElapsed(&plotTriTime, &curTime);
#endif

		if (inconsistent) {
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
			fprintf(stdout, "(%f, %d, %d, %d, %d)\n", curTime, rc.RegChain_UniqueID, currIntTask.RegChain_UniqueID, -1, 0);
#endif
			continue;
		}
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
		fprintf(stdout, "(%f, %d, %d, %d, %ld)\n", curTime, rc.RegChain_UniqueID, currIntTask.RegChain_UniqueID, currIntTask.dimension(), cleanedPolys.size());
#endif

		if (cleanedPolys.size() == 0) {
#if defined(RC_TRIANGULARIZE_PROGRESS) && RC_TRIANGULARIZE_PROGRESS
			fprintf(stderr, "[RegChain] Task Finished: [%lu, %d] -> [%d]\n", polys.size(), rc.dimension(), currIntTask.dimension());
#endif
			results->push_back(currIntTask);
		} else {
#if defined(RC_TRIANGULARIZE_PROGRESS) && RC_TRIANGULARIZE_PROGRESS
				fprintf(stderr, "[RegChain] Scheduling Task: [%lu, %d]\n", cleanedPolys.size(), currIntTask.dimension());
#endif

			if (RC_GEN_COMSUMER_IS_LAST_ITER(intTasks, i)) {
				triangularizeTask(currIntTask, cleanedPolys, lazardDecompose, heightBound, tasks, results);
			} else {
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
				tasks->addPriorityTask((&triangularizeTask<Field,RecursivePoly>), currIntTask, cleanedPolys, lazardDecompose, heightBound, tasks, results);
#else
				tasks->addTask((&triangularizeTask<Field,RecursivePoly>), currIntTask, cleanedPolys, lazardDecompose, heightBound, tasks, results);
#endif
			}
		}
	}
#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	if (noComponents) {
		int tmpID = RC_UNIQUE_ID.fetch_add(1);
		float curTime = 0.0f;
		_stopTimerAddElapsed(&plotTriTime, &curTime);
		fprintf(stdout, "(%f, %d, %d, %d, %d)\n", curTime, rc.RegChain_UniqueID, tmpID, -1, 0);

	}
#endif



}


template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::_triangularizeByTasks(std::vector<RecursivePoly>& polys, bool lazardDecompose, int heightBound) const {

	// typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;

#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
	RC_UNIQUE_ID = 0;
	this->RegChain_UniqueID = RC_UNIQUE_ID.fetch_add(1);
	fprintf(stdout, "(%f, %d, %d, %d, %ld)\n", 0.0f, -1, this->RegChain_UniqueID, this->dimension(), polys.size());
	_startTimer(&plotTriTime);
#endif

	TaskScheduler* tasks = new TaskScheduler;
	SyncVector_ptr<RegularChain<Field,RecursivePoly>> results(new SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>);

	//this recursively schedules resulting tasks via TaskScheduler and its own call stack
	triangularizeTask<Field,RecursivePoly>(*this, polys, lazardDecompose, heightBound, tasks, results);

	//a key synchronization point if tasks split.
	tasks->waitForAllTasks();

	delete tasks;

	// while (tasks.size() > 0) {
	// 	curTask = tasks.front();
	// 	tasks.pop_front();
	// 	std::cerr << "\n\nCurrent Task eqs: \n";
	// 	for (auto curP : curTask.first) {
	// 		std::cerr << curP << std::endl;
	// 		printVariables(curP.ringVariables(), "curP RingVars ");
	// 	}
	// 	std::cerr << "Current rc : \n" << curTask.second << std::endl;
	// 	if (curTask.first.size() == 0) {
	// 		results.push_back(curTask.second);
	// 	} else {

	// 		// std::cerr << "before sort: \n";
	// 		// for (auto curP : curTask.first) {
	// 		// 	std::cerr << curP << std::endl;
	// 		// }
	// 		sortDecreasingStrictRitt<Field,RecursivePoly>(curTask.first, curTask.second.allVariables());
	// 		// std::cerr << "after sort: \n";
	// 		// for (auto curP : curTask.first) {
	// 		// 	std::cerr << curP << std::endl;
	// 		// }

	// 		p = curTask.first.back();
	// 		curTask.first.pop_back();

	// 		std::cerr << "Minimal poly to intersect: \n" << p << std::endl;


	// 		RC_TRI_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(curTask.second), p, lazardDecompose, heightBound);
	// 		RC_TRI_CONSUMER_LOOP(intTasks, currIntTask, i) {
	// 			RC_TRI_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, i);
	// 			curPolys = curTask.first;

	// 			std::cerr << "Before clean set: " << std::endl;
	// 			for (auto tmpP : curPolys) {
	// 				std::cerr << tmpP << std::endl;
	// 				printVariables(tmpP.ringVariables(), "tmpP RingVars ");
	// 			}
	// 			std::cerr << "RC: " << currIntTask << std::endl;
	// 			bool inconsistent = currIntTask.cleanSet(curPolys);
	// 			std::cerr << "After clean set: " << std::endl;
	// 			for (auto tmpP : curPolys) {
	// 				std::cerr << tmpP << std::endl;
	// 				printVariables(tmpP.ringVariables(), "tmpP RingVars ");
	// 			}
	// 			std::cerr << "RC: " << currIntTask << std::endl << std::endl;
	// 			if (inconsistent) {
	// 				std::cerr << "Removing the following component for inconsistency:" << std::endl;
	// 				std::cerr << "RC: " << currIntTask << std::endl;
	// 				continue;
	// 			}
	// 			tasks.emplace_back(std::move(curPolys), std::move(currIntTask));
	// 		}

	// 	}
	// }
	return results->moveVectorOut();
}


//template <class Field, class RecursivePoly>
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::triangularizeRecursive(const std::vector<RecursivePoly>& polys) const {
//	// !!! Current version of the algorithm assumes that *this is an empty regular chain !!! //
//	std::cerr << "entering triangularize..." << std::endl;
//	std::vector<RegularChain<Field,RecursivePoly>> results;
//
//	if (polys.empty()) {
//		results.push_back(*this);
//		return results;
//	}
//
//	std::vector<RecursivePoly> F;
//
//	// TODO: write bool cleanSet(RecursivePoly&) to check the polys and the current object for consistency, and reduce polys modulo *this
//	//        - then if inconsistent return empty rc, if consistent proceed with the reduced polys
//
//	F = sortIncreasingByRank<Field,RecursivePoly>(polys);
//
//	for (size_t i=0; i<F.size(); ++i)
//		std::cerr << "F[" << i << "] = " << F[i] << std::endl;
//
//	results = triangularize_inner(F);
//
//	// TODO: write void removeRedundantRegularChains(std::vector<RegularChain<Field,RecursivePoly>>& results)
//	return results;
//}

//template <class Field, class RecursivePoly>
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::triangularize_inner(std::vector<RecursivePoly>& F) const {
//	// !!! Current version of the algorithm assumes that *this is an empty regular chain !!! //
//	std::cerr << "entering triangularize_inner..." << std::endl;
//	RegularChain<Field,RecursivePoly> temp(*this);
//	std::cerr << "*this = " << temp << std::endl;
//	std::vector<RegularChain<Field,RecursivePoly>> results;
//
//	if (F.empty()) {
//		std::cerr << "Empty polynomial set, returning current object..." << std::endl;
//		results.push_back(*this);
//		std::cerr << "results.size = " << results.size() << std::endl;
//		return results;
//	}
//
//	std::vector<RegularChain<Field,RecursivePoly>> newResults,TT;
//
//	// TODO: write bool cleanSet(RecursivePoly&) to check the polys and the current object for consistency, and reduce polys modulo *this
//	//        - then if inconsistent return empty rc, if consistent proceed with the reduced polys
//
//	RecursivePoly p(F.back());
//	F.pop_back();
//
//	TT = triangularize_inner(F);
//	std::cerr << "returning to triangularize_inner..." << p << std::endl;
//	temp = *this;
//	std::cerr << "*this = " << temp << std::endl;
//	std::cerr << "TT.size = " << TT.size() << std::endl;
//
//	for (size_t i=0; i<TT.size(); ++i) {
//		std::cerr << "triangularize_inner: calling intersect..." << p << std::endl;
//		std::cerr << "TT[" << i << "] = " << TT[i] << std::endl;
//		newResults = TT;
////		for (size_t j=0; j<newResults.size(); ++j)
////			cout << "newResults[" << j << "] = " << newResults[j] << endl;
//		newResults = TT[i].intersect(p);
////		for (size_t j=0; j<newResults.size(); ++j)
////			cout << "newResults[" << j << "] = " << newResults[j] << endl;
//
//		std::cerr << "returning to triangularize_inner from intersect..." << p << std::endl;
//		results.insert(results.end(),newResults.begin(),newResults.end());
//	}
//
//	return results;
//}


template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::squareFreePart(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int inputHeightBound, int options) const {

	typedef RegularChain<Field,RecursivePoly> RC_SFP_OBJ;

	std::vector<RC_SFP_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_SFP_OBJ, squareFreePartResults, (&RegularChain<Field, RecursivePoly>::_squareFreePart), this, p, v, lazardDecompose, inputHeightBound, options);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_SFP_OBJ, currSquareFreePartResult);
	RC_GEN_CONSUMER_LOOP(squareFreePartResults, currSquareFreePartResult, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(squareFreePartResults, currSquareFreePartResult, i);
		results.push_back(currSquareFreePartResult);
	}

	return results;
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_squareFreePart(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int inputHeightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::_squareFreePart(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int inputHeightBound, int options) const {
#endif
	long long unsigned int rcProfilingStart;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqF][v1][" << depth << "] entering _squareFreePart: " << p << ", " << v << std::endl;
		std::cerr << "[sqF][v1][" << depth << "] inputHeightBound = " << inputHeightBound << std::endl;
	#endif

	// The algorithm assumes that init(p) is regular mod Sat(*this), but uncomment to perform a regularity check.
//	if (!isRegular(p.initial())) {
//		std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
//		std::cerr << "in constructChain(p,opts) of RC" << std::endl;
//		exit(1);
//	}

	int heightBound;
	if (lazardDecompose && inputHeightBound == 0)
		heightBound = this->numberOfAlgebraicVariables();
	else
		heightBound = inputHeightBound;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqF][v1][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif

	typedef RegularChain<Field,RecursivePoly> RC_SFP_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_SFP_OBJ, results);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_SFP_OBJ, currSquareFreePartResult);
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);
	std::vector<RegularChain<Field,RecursivePoly>> moreResults;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqF][v1][" << depth << "] computing squareFreePart of p:" << std::endl;
	#endif

	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	RecursivePoly q = p.squareFreePart();
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
	#endif

	if (p.isOne()) {
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqF][v1][" << depth << "] squareFreePart of p computed:" << std::endl;
	#endif
	options = options | ASSUME_SQUAREFREE;

	if (q.mainDegree() == 1 || this->isEmpty()) {
		if (this->numberOfAlgebraicVariables() < heightBound) {
//			RegularChain<Field,RecursivePoly> T(*this);
//			T.constructChain(q,options);

			RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), this, q, lazardDecompose, heightBound, options);

			RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, i) {
				RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, i);

				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqF][v1][" << depth << "] T = " << currConstructChainsTask << " added to results." << std::endl;
				#endif
				RC_GEN_PRODUCER_ACCUMULATE(results, std::move(currConstructChainsTask));
	//			results.emplace_back(T);
			}
		}
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[sqF][v1][" << depth << "] leaving _squareFreePart: " << p << ", " << v << std::endl;
		#endif
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		SubResultantChain<RecursivePoly,RecursivePoly> src(q.derivative(v),q,v,this->allVariables());
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&subresultantChainTime);
			startTimer(&rcProfilingStart);
		#endif

		RC_GEN_CONSUMER_INIT(RC_SFP_OBJ, squareFreePartResults, (&RegularChain<Field, RecursivePoly>::_squareFreePartInner), this, q, v, src, lazardDecompose, heightBound, options);
//		moreResults = _squareFreePartInner(q,v,src,lazardDecompose,heightBound,options);
//		results = _squareFreePartInner(q,v,src,lazardDecompose,heightBound,options);
		RC_GEN_CONSUMER_LOOP(squareFreePartResults, currSquareFreePartResult, i) {
//		for (auto rc : moreResults) {
			RC_GEN_CONSUMER_GET_LOOPELEM(squareFreePartResults, currSquareFreePartResult, i);
			RC_GEN_PRODUCER_ACCUMULATE(results, currSquareFreePartResult);
		}
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[sqF][v1][" << depth << "] " << results.size() << " squareFreePart chains computed " << v << std::endl;
			std::cerr << "[sqF][v1][" << depth << "] leaving _squareFreePart: " << p << ", " << v << std::endl;
		#endif
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_squareFreePartInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::_squareFreePartInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options) const {
#endif
	long long unsigned int rcProfilingStart;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqF][v2][" << depth << "] entering _squareFreePartInner: " << p << ", " << v << std::endl;
		std::cerr << "[sqF][v2][" << depth << "] heightBound " << heightBound << std::endl;
		for (size_t i=0; i<src.size(); ++i) {
			std::cerr << "[sqF][v2][" << depth << "] S[" << i << "] = " << src.subResultantOfIndex(i) << std::endl;
		}
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);
	typedef RegularChain<Field,RecursivePoly> RC_SFP_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_SFP_OBJ, results);
	std::vector<RegularChain<Field,RecursivePoly>> tasks,moreTasks;
	typedef RC_SFP_OBJ RC_INT_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);
	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;
	typedef RC_REGGCD_OBJ RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGGCD_OBJ, currRegGCDTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask2);
//	std::vector<RegularChain<Field,RecursivePoly>> results,tasks,moreTasks;
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> gcdComponents,regularComponents,moreRegularComponents;
	RegularChain<Field,RecursivePoly> currChain,T;
	RecursivePoly r(src.resultant());
	r.setRingVariables(this->allVariables());
	options = options | ASSUME_REGULAR | ASSUME_SQUAREFREE;

	// TODO: Maple organizes this computation differently (see code for TRDsquare_free_poly_src)

	// TODO: do we need a check to verify that this->numberOfAlgebraicVariables() < heightBound before the routine starts?
	// TODO: note that Maple doesn't call regularize unless absolutely necessary (see TRDregular code in solve.mm
	RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), this, r, lazardDecompose, heightBound-1);
//	regularComponents = regularize(r,lazardDecompose,heightBound-1);
	RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, i) {
//	for (size_t i=0; i<regularComponents.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, i);
		const RecursivePoly& f = currRegTask.poly;
		const RegularChain<Field,RecursivePoly>& C = currRegTask.chain;
//		const RecursivePoly& f = regularComponents[i].poly;
//		const RegularChain<Field,RecursivePoly>& C = regularComponents[i].chain;

		if (!C.isInSaturatedIdealMinimal(f)) {
			// C U p
//			currRegTask.chain.constructChain(p,options);

			// TODO: check whether we need to decrease the height bound; any changes to options needed?
			RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &currRegTask.chain, p, lazardDecompose, heightBound, options);

			RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, j) {
				RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, j);
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqF][v2][" << depth << "] chain = " << currConstructChainsTask << " added to results." << std::endl;
				#endif
				RC_GEN_PRODUCER_ACCUMULATE(results, std::move(currConstructChainsTask));
//				RC_GEN_PRODUCER_ACCUMULATE(results, std::move(currRegTask.chain));
	//			results.emplace_back(std::move(regularComponents[i].chain));
			}
		}
		else {
			if (C.dimension() == this->dimension()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqF][v2][" << depth << "] no dimension drop, adding task:" << std::endl;
				#endif
				tasks.emplace_back(std::move(currRegTask.chain));
//				tasks.emplace_back(std::move(regularComponents[i].chain));
			}
			else {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqF][v2][" << depth << "] dimension drop, calling regularize" << std::endl;
				#endif
				// TODO: Maple calls TRDregular_only here (see code in solve.mm)
				RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks2, (&RegularChain<Field, RecursivePoly>::_regularize), this, p.initial(), lazardDecompose, heightBound-1);
//				moreRegularComponents = C.regularize(p.initial(),lazardDecompose,heightBound-1);
				RC_GEN_CONSUMER_LOOP(regTasks2, currRegTask2, j) {
//				for (size_t j=0; j<moreRegularComponents.size(); ++j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(regTasks2, currRegTask2, j);
					if (!currRegTask2.poly.isZero()) {
//					if (!moreRegularComponents[j].poly.isZero()) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[sqF][v2][" << depth << "] adding task:" << std::endl;
						#endif
						tasks.emplace_back(std::move(currRegTask2.chain));
//						tasks.emplace_back(std::move(moreRegularComponents[j].chain));
					}
				}
			}
		}
	}

	while (!tasks.empty()) {
		currChain = tasks.back();
		tasks.pop_back();
		// TODO: This seems problematic since only one of the two inputs is known to have regular initial modulo the chain.

		RC_GEN_CONSUMER_INIT(RC_REGGCD_OBJ, regGCDTasks, (&RegularChain<Field, RecursivePoly>::_regularGCD), &(currChain), src.firstPolynomial(), src.secondPolynomial(), v, src, lazardDecompose, heightBound-1);
//		regularComponents = currChain.regularGCD(src.firstPolynomial(),src.secondPolynomial(),v,src,lazardDecompose,heightBound-1);

		RC_GEN_CONSUMER_LOOP(regGCDTasks, currRegGCDTask, i) {
//		for (size_t i=0; i<regularComponents.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(regGCDTasks, currRegGCDTask, i);
			const RecursivePoly& g = currRegGCDTask.poly;
			const RegularChain<Field,RecursivePoly>& D = currRegGCDTask.chain;
//			const RecursivePoly& g = regularComponents[i].poly;
//			const RegularChain<Field,RecursivePoly>& D = regularComponents[i].chain;

			if (D.dimension() == currChain.dimension()) {
				RecursivePoly q;
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqF][v2][" << depth << "] calling pseudoDivide from SMQP:" << std::endl;
				#endif
				#ifdef REGULARCHAIN_PROFILING
					startTimer(&rcProfilingStart);
				#endif
				p.pseudoDivide(g,&q,NULL,1);
				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
				#endif

				if (this->isConstantPolynomial(q) || q.leadingVariable() != v) {
					std::cerr << "BPAS: error, constant or degree zero in v pseudoquotient in squareFreePart." << std::endl;
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "q = " << q << std::endl;
					#endif
					exit(1);
				}

				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqF][v2][" << depth << "] result added in while loop:" << std::endl;
				#endif
				T = D;
//				T.constructChain(q,options);

				// TODO: check whether we need to decrease the height bound; any changes to options needed?
				RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &T, q, lazardDecompose, heightBound, options);

				RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, i) {
					RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, i);
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[sqF][v2][" << depth << "] T = " << currConstructChainsTask << " added to results." << std::endl;
					#endif

					RC_GEN_PRODUCER_ACCUMULATE(results, std::move(currConstructChainsTask));
	//				results.emplace_back(std::move(T));
				}

				RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), g.initial(), lazardDecompose, heightBound-1);
//				moreTasks = D.intersect(g.initial(),lazardDecompose,heightBound-1);

				RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, j) {
//				for (size_t j=0; j<moreTasks.size(); ++j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
					// TODO: Maple calls TRDregular_only here (see code in solve.mm)
					RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(currIntTask), p.initial(), lazardDecompose, heightBound-1);
//					moreRegularComponents = currIntTask.regularize(p.initial(),lazardDecompose,heightBound-1);
//					moreRegularComponents = moreTasks[j].regularize(p.initial(),lazardDecompose,heightBound-1);

					RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, k) {
//					for (size_t k=0; k<moreRegularComponents.size(); ++k) {
						RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, k);
						if (!currRegTask.poly.isZero()) {
//						if (!moreRegularComponents[k].poly.isZero()) {
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[sqF][v2][" << depth << "] while loop, equal dimension, adding task:" << std::endl;
							#endif
							tasks.emplace_back(std::move(currRegTask.chain));
//							tasks.emplace_back(std::move(moreRegularComponents[k].chain));
						}
					}
				}
			}
			else {
//				moreRegularComponents = D.regularize(p.initial());
				if (D.numberOfAlgebraicVariables() < heightBound) {
					// TODO: Maple calls TRDregular_only here (see code in solve.mm)
					RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(D), p.initial(), lazardDecompose, heightBound-1);
//					moreRegularComponents = D.regularize(p.initial(),lazardDecompose,heightBound-1);

					RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, j) {
//					for (size_t j=0; j<moreRegularComponents.size(); ++j) {
						RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, j);
						if (!currRegTask.poly.isZero()) {
//						if (!moreRegularComponents[j].poly.isZero()) {
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[sqF][v2][" << depth << "] while loop, dimension drop, adding task:" << std::endl;
							#endif
							tasks.emplace_back(std::move(currRegTask.chain));
//							tasks.emplace_back(std::move(moreRegularComponents[j].chain));
						}
					}
				}
			}
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqF][v2][" << depth << "] leaving _squareFreePartInner: " << p << ", " << v << std::endl;
	#endif
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_squareFreePartPolynomial(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int inputHeightBound, int options, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::_squareFreePartPolynomial(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int inputHeightBound, int options) const {
#endif
	long long unsigned int rcProfilingStart;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqFPP][" << depth << "] entering _squareFreePartPolynomial: " << p << ", " << v << std::endl;
		std::cerr << "[sqFPP][" << depth << "] inputHeightBound = " << inputHeightBound << std::endl;
	#endif

	// The algorithm assumes that init(p) is regular mod Sat(*this), but uncomment to perform a regularity check.
//	if (!isRegular(p.initial())) {
//		std::cerr << "BPAS: error, initial of p is not regular with respect to the current regular chain T so T+p is not a regular chain." << std::endl;
//		std::cerr << "in constructChain(p,opts) of RC" << std::endl;
//		exit(1);
//	}

	int heightBound;
	if (lazardDecompose && inputHeightBound == 0)
		heightBound = this->numberOfAlgebraicVariables();
	else
		heightBound = inputHeightBound;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqFPP][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_SFP_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_SFP_OBJ, results);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_SFP_OBJ, currSquareFreePartResult);
	std::vector<RegularChain<Field,RecursivePoly>> moreResults;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqFPP][" << depth << "] computing squareFreePart of p:" << std::endl;
	#endif

	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	RecursivePoly q = p.squareFreePart();
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
	#endif

	if (p.isOne()) {
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqFPP][" << depth << "] squareFreePart of p computed:" << q << std::endl;
	#endif
	options = options | ASSUME_SQUAREFREE;

	if (q.mainDegree() == 1 || this->isEmpty()) {
		if (this->numberOfAlgebraicVariables() < heightBound) {
//			RegularChain<Field,RecursivePoly> T(*this);
//			T.constructChain(q,options);
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[sqFPP][" << depth << "] T = " << *this << " added to results." << std::endl;
			#endif
			RC_GEN_PRODUCER_ACCUMULATE(results, RC_SFP_OBJ(q,*this));
//			results.emplace_back(T);
		}
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[sqFPP][" << depth << "] leaving _squareFreePartPolynomial: " << p << ", " << v << std::endl;
		#endif
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[sqFPP][" << depth << "] computing SRC between q: " << q << std::endl;
			std::cerr << "[sqFPP][" << depth << "] computing SRC and dq/dv: " << q.derivative(v) << std::endl;
		#endif

		SubResultantChain<RecursivePoly,RecursivePoly> src(q.derivative(v),q,v,this->allVariables());
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&subresultantChainTime);
			startTimer(&rcProfilingStart);
		#endif

		RC_GEN_CONSUMER_INIT(RC_SFP_OBJ, squareFreePartResults, (&RegularChain<Field, RecursivePoly>::_squareFreePartPolynomialInner), this, q, v, src, lazardDecompose, heightBound, options);
//		moreResults = _squareFreePartInner(q,v,src,lazardDecompose,heightBound,options);
//		results = _squareFreePartInner(q,v,src,lazardDecompose,heightBound,options);
		RC_GEN_CONSUMER_LOOP(squareFreePartResults, currSquareFreePartResult, i) {
//		for (auto rc : moreResults) {
			RC_GEN_CONSUMER_GET_LOOPELEM(squareFreePartResults, currSquareFreePartResult, i);
			RC_GEN_PRODUCER_ACCUMULATE(results, currSquareFreePartResult);
		}
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[sqFPP][" << depth << "] " << results.size() << " squareFreePart chains computed " << v << std::endl;
			std::cerr << "[sqFPP][" << depth << "] leaving _squareFreePartPolynomial: " << p << ", " << v << std::endl;
		#endif
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_squareFreePartPolynomialInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::_squareFreePartPolynomialInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options) const {
#endif
	long long unsigned int rcProfilingStart;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqFPPI][" << depth << "] entering _squareFreePartPolynomialInner: " << p << ", " << v << std::endl;
		std::cerr << "[sqFPPI][" << depth << "] heightBound " << heightBound << std::endl;
		for (size_t i=0; i<src.size(); ++i) {
			std::cerr << "[sqFPPI][" << depth << "] S[" << i << "] = " << src.subResultantOfIndex(i) << std::endl;
		}
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;
	std::vector<RegularChain<Field,RecursivePoly>> tasks,moreTasks;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);
	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_SFP_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_SFP_OBJ, results);
	typedef RC_SFP_OBJ RC_REGGCD_OBJ;
	typedef RC_REGGCD_OBJ RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGGCD_OBJ, currRegGCDTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask2);
//	std::vector<RegularChain<Field,RecursivePoly>> results,tasks,moreTasks;
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> gcdComponents,regularComponents,moreRegularComponents;
	RegularChain<Field,RecursivePoly> currChain,T;
	RecursivePoly r(src.resultant());
	r.setRingVariables(this->allVariables());
	options = options | ASSUME_REGULAR | ASSUME_SQUAREFREE;

	// TODO: Maple organizes this computation differently (see code for TRDsquare_free_poly_src)

	// TODO: do we need a check to verify that this->numberOfAlgebraicVariables() < heightBound before the routine starts?
	// TODO: note that Maple doesn't call regularize unless absolutely necessary (see TRDregular code in solve.mm
	RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), this, r, lazardDecompose, heightBound-1);
//	regularComponents = regularize(r,lazardDecompose,heightBound-1);
	RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, i) {
//	for (size_t i=0; i<regularComponents.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, i);
		const RecursivePoly& f = currRegTask.poly;
		const RegularChain<Field,RecursivePoly>& C = currRegTask.chain;
//		const RecursivePoly& f = regularComponents[i].poly;
//		const RegularChain<Field,RecursivePoly>& C = regularComponents[i].chain;

		if (!C.isInSaturatedIdealMinimal(f)) {
			// C U p
//			currRegTask.chain.constructChain(p,options);
////			regularComponents[i].chain.constructChain(p,options);

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[sqFPPI][" << depth << "] chain = " << currRegTask.chain << " added to results." << std::endl;
			#endif
			RC_GEN_PRODUCER_ACCUMULATE(results, RC_SFP_OBJ(p,currRegTask.chain));
//			RC_GEN_PRODUCER_ACCUMULATE(results, std::move(regularComponents[i].chain));
//			results.emplace_back(std::move(regularComponents[i].chain));
		}
		else {
			if (C.dimension() == this->dimension()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqFPPI][" << depth << "] no dimension drop, adding task:" << std::endl;
				#endif
				tasks.emplace_back(std::move(currRegTask.chain));
//				tasks.emplace_back(std::move(regularComponents[i].chain));
			}
			else {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqFPPI][" << depth << "] dimension drop, calling regularize" << std::endl;
				#endif
				// TODO: Maple calls TRDregular_only here (see code in solve.mm)
				RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks2, (&RegularChain<Field, RecursivePoly>::_regularize), this, p.initial(), lazardDecompose, heightBound-1);
//				moreRegularComponents = C.regularize(p.initial(),lazardDecompose,heightBound-1);
				RC_GEN_CONSUMER_LOOP(regTasks2, currRegTask2, j) {
//				for (size_t j=0; j<moreRegularComponents.size(); ++j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(regTasks2, currRegTask2, j);
					if (!currRegTask2.poly.isZero()) {
//					if (!moreRegularComponents[j].poly.isZero()) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[sqFPPI][" << depth << "] adding task:" << std::endl;
						#endif
						tasks.emplace_back(std::move(currRegTask2.chain));
//						tasks.emplace_back(std::move(moreRegularComponents[j].chain));
					}
				}
			}
		}
	}

	while (!tasks.empty()) {
		currChain = tasks.back();
		tasks.pop_back();
		// TODO: This seems problematic since only one of the two inputs is known to have regular initial modulo the chain.

		RC_GEN_CONSUMER_INIT(RC_REGGCD_OBJ, regGCDTasks, (&RegularChain<Field, RecursivePoly>::_regularGCD), &(currChain), src.firstPolynomial(), src.secondPolynomial(), v, src, lazardDecompose, heightBound-1);
//		regularComponents = currChain.regularGCD(src.firstPolynomial(),src.secondPolynomial(),v,src,lazardDecompose,heightBound-1);

		RC_GEN_CONSUMER_LOOP(regGCDTasks, currRegGCDTask, i) {
//		for (size_t i=0; i<regularComponents.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(regGCDTasks, currRegGCDTask, i);
			const RecursivePoly& g = currRegGCDTask.poly;
			const RegularChain<Field,RecursivePoly>& D = currRegGCDTask.chain;
//			const RecursivePoly& g = regularComponents[i].poly;
//			const RegularChain<Field,RecursivePoly>& D = regularComponents[i].chain;

			if (D.dimension() == currChain.dimension()) {
				RecursivePoly q;
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqFPPI][" << depth << "] calling pseudoDivide from SMQP:" << std::endl;
					std::cerr << "[sqFPPI][" << depth << "] p = " << p << std::endl;
					std::cerr << "[sqFPPI][" << depth << "] g = " << g << std::endl;

				#endif
				#ifdef REGULARCHAIN_PROFILING
					startTimer(&rcProfilingStart);
				#endif
				p.pseudoDivide(g,&q,NULL,1);
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqFPPI][" << depth << "] obtained pseudoquotient:" << q << std::endl;
				#endif
				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
				#endif

				if (this->isConstantPolynomial(q) || q.leadingVariable() != v) {
					std::cerr << "BPAS: error, constant or degree zero in v pseudoquotient in squareFreePartPolynomialInner." << std::endl;
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "q = " << q << std::endl;
					#endif
					exit(1);
				}

				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqFPPI][" << depth << "] result added in while loop:" << std::endl;
				#endif
//				T = D;
//				T.constructChain(q,options);
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[sqFPPI][" << depth << "] T = " << D << " added to results." << std::endl;
				#endif

				RC_GEN_PRODUCER_ACCUMULATE(results, RC_SFP_OBJ(q,D));
//				RC_GEN_PRODUCER_ACCUMULATE(results, std::move(T));
//				results.emplace_back(std::move(T));

				RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), g.initial(), lazardDecompose, heightBound-1);
//				moreTasks = D.intersect(g.initial(),lazardDecompose,heightBound-1);

				RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, j) {
//				for (size_t j=0; j<moreTasks.size(); ++j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
					// TODO: Maple calls TRDregular_only here (see code in solve.mm)
					RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(currIntTask), p.initial(), lazardDecompose, heightBound-1);
//					moreRegularComponents = currIntTask.regularize(p.initial(),lazardDecompose,heightBound-1);
//					moreRegularComponents = moreTasks[j].regularize(p.initial(),lazardDecompose,heightBound-1);

					RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, k) {
//					for (size_t k=0; k<moreRegularComponents.size(); ++k) {
						RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, k);
						if (!currRegTask.poly.isZero()) {
//						if (!moreRegularComponents[k].poly.isZero()) {
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[sqFPPI][" << depth << "] while loop, equal dimension, adding task:" << std::endl;
							#endif
							tasks.emplace_back(std::move(currRegTask.chain));
//							tasks.emplace_back(std::move(moreRegularComponents[k].chain));
						}
					}
				}
			}
			else {
//				moreRegularComponents = D.regularize(p.initial());
				if (D.numberOfAlgebraicVariables() < heightBound) {
					// TODO: Maple calls TRDregular_only here (see code in solve.mm)
					RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(D), p.initial(), lazardDecompose, heightBound-1);
//					moreRegularComponents = D.regularize(p.initial(),lazardDecompose,heightBound-1);

					RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, j) {
//					for (size_t j=0; j<moreRegularComponents.size(); ++j) {
						RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, j);
						if (!currRegTask.poly.isZero()) {
//						if (!moreRegularComponents[j].poly.isZero()) {
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[sqFPPI][" << depth << "] while loop, dimension drop, adding task:" << std::endl;
							#endif
							tasks.emplace_back(std::move(currRegTask.chain));
//							tasks.emplace_back(std::move(moreRegularComponents[j].chain));
						}
					}
				}
			}
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[sqFPPI][" << depth << "] leaving _squareFreePartPolynomialInner: " << p << ", " << v << std::endl;
	#endif
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

//template <class Field, class RecursivePoly>
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::squareFreePart() const {
//	++depth;
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[sqF][v3][" << depth << "] entering squareFreePart(): " << std::endl;
//	#endif
//	std::vector<RegularChain<Field,RecursivePoly>> results;
//
//	for (size_t i=0; i<set.size(); ++i) {
//		if (!set[i].isZero())
//			set[i] = set[i].squareFreePart(set[i].leadingVariable());
//	}
//
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[sqF][v3][" << depth << "] leaving squareFreePart(): " << std::endl;
//	#endif
//	--depth;
//	return results;
//}

template <class Field, class RecursivePoly>
std::vector<RecursivePoly> RegularChain<Field,RecursivePoly>::GCDFreeFactorization(const RecursivePoly& p, int type) const {
	long long unsigned int rcProfilingStart;
	long long unsigned int GCDFreeFactorizationStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&GCDFreeFactorizationStart);
	#endif
	// TODO: Change to compute an irreducible factorization when this is available.
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[GCDFreeFact] entering GCDFreeFactorization: " << p << std::endl;
	#endif
	std::vector<RecursivePoly> results;
	RecursivePoly q;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[GCDFreeFact] computing integer primitivePart: " << p << std::endl;
	#endif
	q = p.primitivePart();
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
	#endif
	if (type == 0) {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		results.emplace_back(q.squareFreePart());
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
	}
	else if (type == 1) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] computing primitive factorization..." << std::endl;
		#endif
		RecursivePoly c;
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		q = q.mainPrimitivePart(c);
//		c = c.primitivePart();
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
//		results.emplace_back(c);
//		results.emplace_back(q);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] computing squarefree part of c and q" << std::endl;
			std::cerr << "[GCDFreeFact] c = " << c << std::endl;
			std::cerr << "[GCDFreeFact] q = " << q << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		results.emplace_back(c.squareFreePart());
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] c = " << c << std::endl;
		#endif
		results.emplace_back(q.squareFreePart());
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] q = " << q << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
		#endif
	}
	else if (type == 2) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] computing irreducible factorization..." << std::endl;
		#endif
		Factors<RecursivePoly> facts;
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] calling factor..." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
			// fprintf(stderr, "FACTORING!\n" );
			facts = q.factor();
			// fprintf(stderr, "\n\n\n\n\n\n\nDONE FACTORING!\n" );
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&factorTime);
		#endif
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] packaging results..." << std::endl;
		#endif
		results.emplace_back(facts.ringElement());
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&rcProfilingStart);
//		#endif
//		results.push_back(facts.ringElement().primitivePart());
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
//		#endif
		for (size_t i=0; i<facts.size(); ++i) {
			results.emplace_back(facts.factor(i).first);
		}
		// RecursivePoly check,temp;
		// check = facts.ringElement();
		// for (size_t j=0; j<facts.size(); ++j) {
		// 	temp = facts.factor(j).first^facts.factor(j).second;
		// 	check *= temp;
		// }
		// if (check != q) {
		// 	std::cerr << "BLAD error, factorization is wrong." << std::endl;
		// 	std::cerr << "check = " << check << std::endl;
		// 	std::cerr << "q = " << q << std::endl;
		// 	std::cerr << "ringElement = " << facts.ringElement() << std::endl;
		// 	for (size_t m=0; m<facts.size(); ++m) {
		// 		std::cerr << "factor[" << m <<  "] = [" << facts.factor(m).first << "," << facts.factor(m).second << "]" << std::endl;
		// 	}
		// 	printVariables(q.ringVariables(), "q ringVars");
		// 	printVariables(check.ringVariables(), "check ringVars");
		// 	fprintf(stderr, "\n\n\nBAD FACTORIZATION\n" );
		// 	exit(1);
		// }
//		else {
//			std::cerr << "[GCDFreeFact] check = " << check << std::endl;
//			std::cerr << "[GCDFreeFact] q = " << q << std::endl;
//		}
	}
	else {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[GCDFreeFact] computing squarefree factorization..." << std::endl;
		#endif
		Factors<RecursivePoly> facts;
		std::vector<Symbol> vs;
		vs.emplace_back(p.leadingVariable());
		facts = q.squareFree(vs);
		results.emplace_back(facts.ringElement());
//		#ifdef REGULARCHAIN_PROFILING
//			startTimer(&rcProfilingStart);
//		#endif
//		results.push_back(facts.ringElement().primitivePart());
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
//		#endif
		for (size_t i=0; i<facts.size(); ++i) {
			results.emplace_back(facts.factor(i).first);
		}
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[GCDFreeFact] leaving GCDFreeFactorization." << std::endl;
		for (size_t i=0; i<results.size(); ++i)
			std::cerr << "[GCDFreeFact] results[" << i << "] = " << results[i] << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&GCDFreeFactorizationStart,&GCDFreeFactorizationTime);
	#endif
	return results;
}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::intersect(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound) const {

	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;

	std::vector<RC_INT_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intResults, (&RegularChain<Field, RecursivePoly>::_intersect), this, p, lazardDecompose, inputHeightBound);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntResult);
	RC_GEN_CONSUMER_LOOP(intResults, currIntResult, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(intResults, currIntResult, i);
		results.push_back(currIntResult);
	}

	return results;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isIntersectionTrivial(const RecursivePoly& p, RecursivePoly& pReduced) const {

	if (p.isZero() || this->isConstantPolynomial(p)) {
		pReduced = p;
		return true;
	}
	RecursivePoly z = this->pseudoDivide(p);
	if (z.isZero()) {
		pReduced = z;
		return true;
	}
	return false;
}


template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::intersectTrivial(const RecursivePoly& p) const {
	long long unsigned int rcProfilingStart;
	long long unsigned int intersectStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&intersectStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[intTriv] entering intersectTrivial: " << p << ", " << *this << std::endl;
	#endif
	if (!this->isConstantPolynomial(p)) {
		std::cerr << "BPAS: error, cannot call intersectTrivial with non-constant input." << std::endl;
		exit(1);
	}
	std::vector<RegularChain<Field,RecursivePoly>> results;

	if (p.isZero()) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[intTriv] intersectTrivial: input is zero, returning current object..." << std::endl;
		#endif
		results.emplace_back(*this);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectStart,&intersectTime);
		#endif
		return results;
	}
	else {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[intTriv] intersectTrivial: constant polynomial, returning nothing..." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectStart,&intersectTime);
		#endif
		return results;
	}
}


#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_intersect(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::_intersect(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound) const {
#endif
	long long unsigned int rcProfilingStart;
	long long unsigned int intersectStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&intersectStart);
	#endif
	++intersectDepth;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] entering intersect: " << p << ", " << *this << std::endl;
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] inputHeightBound = " << inputHeightBound << std::endl;
	#endif
	bool doAsMapleDoes(true);
	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_INT_OBJ, results);
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);
//	std::vector<RegularChain<Field,RecursivePoly>> results;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] checking zero: " << std::endl;
	#endif
	if (p.isZero()) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] intersect: input is zero, returning current object..." << std::endl;
		#endif
		RC_GEN_PRODUCER_ACCUMULATE(results, *this);
//		results.emplace_back(*this);
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] leaving intersect: " << p << ", " << *this << std::endl;
			for (size_t i=0; i<results.size(); ++i) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
			}
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectStart,&intersectTime);
		#endif
		--intersectDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] checking constant: " << std::endl;
	#endif
	if (this->isConstantPolynomial(p)) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] intersect: constant polynomial, returning nothing..." << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectStart,&intersectTime);
		#endif
		--intersectDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] checking in Sat(T): " << std::endl;
	#endif
	// TODO: consider reducing before this and avoiding this test
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	RecursivePoly z = this->pseudoDivide(p);
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
	#endif
	if (z.isZero()) { // TODO: use isInSaturatedIdealMinimal here (or a modular method)
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] intersect: pseudo remainder is zero, returning current object..." << std::endl;
		#endif
		RC_GEN_PRODUCER_ACCUMULATE(results, *this);
//		results.emplace_back(*this);
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
//			std::cerr << "results.size = " << results.size() << std::endl;
//			std::cerr << "[int][" << intersectDepth << "][" << depth << "] leaving intersect:" << p << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] leaving intersect: " << p << ", " << *this << std::endl;
			for (size_t i=0; i<results.size(); ++i) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
			}
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectStart,&intersectTime);
		#endif
		--intersectDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] checking if we can compute in dimension zero:" << std::endl;
	#endif
	RecursivePoly q,c;
	if (canComputeInDimensionZero(p)) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] can compute in dimension zero!!" << std::endl;
		#endif
		// Note: Maple removes the zero polys from the chain, but we use TS_FIXED chains, that can be treated as zero dimensional.
		// TODO: Keep in mind:
				//#The following three lines are critical for example gamewo5 mod pnum
				//if TRDis_invertible_dim0(clean_p, rc, R) then
				//	return ([ ]);
				//end if:
		q = this->reduceMinimal(p);
//		q = this->reduce(p,c);
//		q *= c; // TODO: consider using this primitive factorization to split the computation.
		// Either of the following two lines (in place of the above two) is currently (Aug 28) equivalent in terms of timing.
//		q = this->reduce(p);
//		q = p;
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] q = T.reduceMinimal(p) = " << q << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] T = " << *this << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] T variables:" << std::endl;
			printVariables(this->variables());
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] p = " << p << std::endl;
		#endif
		ZeroDimensionalRegularChain<Field,RecursivePoly> rc(*this,ASSUME_MAKESCHAIN);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] calling intersect in ZDRC..." << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] *this = " << rc << std::endl;
		#endif
		typedef ZeroDimensionalRegularChain<Field,RecursivePoly> RC_ZDINT_OBJ;
		RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_ZDINT_OBJ, currZDIntTask);
//		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> zdresults;
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		RC_GEN_CONSUMER_INIT(RC_ZDINT_OBJ, zdresults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_intersect), &(rc), q);
//		zdresults = rc.intersect(q);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&zerodimensionalregularchainTime);
		#endif
		#ifdef REGULARCHAIN_DEBUG
//			if (zdresults.empty())
//				std::cerr << "[int][" << intersectDepth << "][" << depth << "] no results from intersect in ZDRC" << std::endl;
		#endif
		RC_GEN_CONSUMER_LOOP(zdresults, currZDIntTask, i) {
//		for (size_t i=0; i<zdresults.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(zdresults, currZDIntTask, i);
			#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "[int][" << intersectDepth << "][" << depth << "] zdresults[" << i << "] = " << zdresults[i] << std::endl;
			#endif
			RC_GEN_PRODUCER_ACCUMULATE(results, std::move(currZDIntTask));
//			RC_GEN_PRODUCER_ACCUMULATE(results, std::move(zdresults[i]));
//			results.emplace_back(RegularChain<Field,RecursivePoly>(std::move(zdresults[i])));
//			#ifdef REGULARCHAIN_DEBUG
//				if (this->variables() != T.variables()) {
//					printVariables(rc.variables(),"rc");
//					printVariables(zdresults[i].chain.variables(),"zdresults");
//					printVariables(T.variables(),"T");
//					printVariables(this->variables(),"this");
//					exit(1);
//				}
//			#endif
		}
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
//			std::cerr << "[int][" << intersectDepth << "][" << depth << "] leaving intersect (inDim0): " << p << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] leaving intersect (inDim0): " << p << ", " << *this << std::endl;
			for (size_t i=0; i<results.size(); ++i) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
			}
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectStart,&intersectTime);
		#endif
		--intersectDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}

	// NB: Maple (merge gcd free) factors all the initials of the regular chain at this point (and then removes constants)

	// TODO: Consider splitting into projection and extension subroutines when factorization is available

	// TODO: We should be computing the gcd-free part of r wrt to the list of intials of T at each step of the projection step,
	//       and continuing if r becomes 0 and replacing r with the gcd-free part otherwise.

	// Setting Height Bound for RegularChain objects //

	int heightBound,newHeightBound;
	if (lazardDecompose)
		heightBound = this->numberOfVariables();
	else
		heightBound = inputHeightBound;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif

	// Projection Part: //

	// GCD-Free-Factorize Input //
	struct MetaTask {
		std::vector<RecursivePoly> P;
		std::vector<SubResultantChain<RecursivePoly,RecursivePoly>> S;
	};
	std::vector<RecursivePoly> P;
	std::vector<SubResultantChain<RecursivePoly,RecursivePoly>> S;
	std::vector<MetaTask> projectionTasks,extensionTasks;
	MetaTask mt;
	std::vector<RecursivePoly> lf;
	RecursivePoly r;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] computing gcd-free (squarefree) factorization of p:" << std::endl;
	#endif

//	r = p; // TODO: why??
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&intersectStart,&intersectTime);
	#endif
	lf = GCDFreeFactorization(p);
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&intersectStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] GCD-free factors:" << std::endl;
		for (size_t i=0; i<lf.size(); ++i) {
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] lf[" << i << "] = " << lf[i] << std::endl;
		}
	#endif
	for (size_t i=0; i<lf.size(); ++i) {
		if (!doAsMapleDoes) {
			q = this->reduceMinimal(lf[i]);
//			q = this->reduce(lf[i],c);
//			q *= c;
		}
		else {
			q = lf[i];
		}
		if (!this->isConstantPolynomial(q)) {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] adding " << lf[i] << " to projection tasks." << std::endl;
			#endif
			P.clear();
			P.emplace_back(q);
//			#ifdef REGULARCHAIN_PROFILING
//				startTimer(&rcProfilingStart);
//			#endif
//			P.push_back(lf[i].squareFreePart()); // computing squareFreePart in GCDFreeFactorization instead
//			#ifdef REGULARCHAIN_PROFILING
//				stopTimerAddElapsed(&rcProfilingStart,&squareFreePartTime);
//			#endif
			mt = {P,S};
			projectionTasks.emplace_back(std::move(mt));
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] Starting projection part of intersect..." << p << std::endl;
	#endif
	std::vector<SubResultantChain<RecursivePoly,RecursivePoly>> srcs;
	SubResultantChain<RecursivePoly,RecursivePoly> src;
	std::vector<Symbol> vs(this->variables()),mainVars(this->mainVariables()),resultantsMainVars;
	RecursivePoly f;
	Symbol v;
	int resultantIndex;

	while (projectionTasks.size()) {
		mt = projectionTasks.back();
		projectionTasks.pop_back();
		P = mt.P;
		r = P.back();
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] current r from projection tasks: " << r << std::endl;
		#endif
		S = mt.S;
		// TODO: Determine whether it is better to recompute this or to maintain this list in each task.
		resultantsMainVars.clear();
		for (size_t i=0; i<P.size(); ++i) {
			resultantsMainVars.emplace_back(P[i].leadingVariable());
		}
		#ifdef REGULARCHAIN_DEBUG
			printVariables(resultantsMainVars,"resultantsMainVars");
			printVariables(mainVars,"mainVars");
		#endif
		if (isAMemberOf(r.leadingVariable(),mainVars)) {
			v = r.leadingVariable();
			#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] computing src of " << r << " and " << select(v) << std::endl;
			#endif
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			// TODO: This statement may require changing if it becomes possible to intersect a polynomial that has variables not appearing in T
			src = SubResultantChain<RecursivePoly,RecursivePoly>(r,select(v),v,this->allVariables());
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&subresultantChainTime);
			#endif
			#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] src = " << src << std::endl;
			#endif
			r = src.resultant();
			r.setRingVariables(this->allVariables());
			S.emplace_back(std::move(src));
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] resultant = " << r << std::endl;
			#endif
			// GCD-Free-Factorize the Resultant //
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] calling GCDFreeFactorization of resultant:" << std::endl;
			#endif

			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&intersectStart,&intersectTime);
			#endif

			lf = GCDFreeFactorization(r);

			#ifdef REGULARCHAIN_PROFILING
				startTimer(&intersectStart);
			#endif

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] GCD-free factors:" << std::endl;
				for (size_t i=0; i<lf.size(); ++i) {
					std::cerr << "[int][" << intersectDepth << "][" << depth << "]lf[" << i << "] = " << lf[i] << std::endl;
				}
				std::cerr << "[int][" << intersectDepth << "][" << depth << "]lf.size = " << lf.size() << std::endl;
			#endif
			for (size_t i=0; i<lf.size(); ++i) {
				q = std::move(lf[i]);
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[int][" << intersectDepth << "][" << depth << "]this->isConstantPolynomial(lf[i]) = " << this->isConstantPolynomial(q) << std::endl;
				#endif
				if (!doAsMapleDoes) {
					q = this->reduceMinimal(q);
//					q = this->reduce(q,c);
//					q *= c;
				}
				if (q.isZero()) { // If this condition is met, the iterated resultant is zero, indicating non-trivial intersection.
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] adding " << q << " to extension tasks." << std::endl;
					#endif
					P.emplace_back(std::move(q));
					mt = {P,S};
					extensionTasks.push_back(mt);
					P.pop_back();
				}
				else if (!this->isConstantPolynomial(q)) { // If this condition is not met, the component is empty, so the task is dropped.
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] adding " << q << " to projection tasks." << std::endl;
					#endif
//					P = mt.P;
					P.emplace_back(std::move(q));
					mt = {P,S};
					projectionTasks.push_back(mt);
					P.pop_back();
				}
			}
		}
		else {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] adding " << r << " to extension tasks." << std::endl;
			#endif
			extensionTasks.push_back(mt);
		}
	}
	projectionTasks.clear();



	// Extension Part: //
	typedef RegularChain<Field,RecursivePoly> RC_INTFREE_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INTFREE_OBJ, currIntFreeTask);
	typedef RC_INTFREE_OBJ RC_INTALG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INTALG_OBJ, currIntAlgTask);
	typedef RC_INTFREE_OBJ RC_CLEANCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CLEANCH_OBJ, currCleanChainTask);

//	std::vector<RegularChain<Field,RecursivePoly>> resultTasks,tasks,moreTasks,intAlgTasks,intTasks;
	std::vector<RegularChain<Field,RecursivePoly>> resultTasks,tasks,moreTasks;
	size_t i,j,k,l,n(this->numberOfVariables());
	bool vIsAResultantsMainVariable,vIsARegularChainMainVariable;
	RegularChain<Field,RecursivePoly> rc;
	#ifdef REGULARCHAIN_DEBUG
	std::cerr << "[int][" << intersectDepth << "][" << depth << "] Starting extension part of intersect..." << p << std::endl;

	std::cerr << "[int][" << intersectDepth << "][" << depth << "] extensionTasks.size = " << extensionTasks.size() << std::endl;
	#endif
	for (l=0; l<extensionTasks.size(); ++l) {
		mt = extensionTasks[l];
		P = mt.P;
		S = mt.S;
		#ifdef REGULARCHAIN_DEBUG
			for (size_t ll=0; ll<P.size(); ++ll) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] P[" << ll << "] = " << P[ll] << std::endl;
			}
			for (size_t ll=0; ll<S.size(); ++ll) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] S[" << ll << "] = " << S[ll] << std::endl;
			}
		#endif
	}

	for (l=0; l<extensionTasks.size(); ++l) {
		i = 1;
		mt = extensionTasks[l];
		P = mt.P;
		S = mt.S;
		resultantsMainVars.clear();
		for (size_t i=0; i<P.size(); ++i) {
			resultantsMainVars.emplace_back(P[i].leadingVariable());
		}
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] initializing empty chain:" << std::endl;
			for (size_t ll=0; ll<P.size(); ++ll) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] P[" << ll << "] = " << P[ll] << std::endl;
			}
			for (size_t ll=0; ll<S.size(); ++ll) {
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] S[" << ll << "] = " << S[ll] << std::endl;
			}
		#endif
		resultTasks.clear();
		resultTasks.emplace_back(RegularChain<Field,RecursivePoly>(vars));
		#ifdef REGULARCHAIN_DEBUG
			printVariables(vars);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] Start of while loop" << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] resultTasks.size = " << resultTasks.size() << std::endl;
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] n = " << n << std::endl;
		#endif

		while (i <= n) {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] i = " << i << std::endl;
			#endif

			v = vs[n-i];
			this->upper(v,rc);
			newHeightBound = heightBound - rc.numberOfAlgebraicVariables();

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] newHeightBound = " << newHeightBound << std::endl;
			#endif

			vIsAResultantsMainVariable = isAMemberOf(v,resultantsMainVars);
			vIsARegularChainMainVariable = isAMemberOf(v,mainVars);

			for (j=0; j<resultTasks.size(); ++j) {
				const RegularChain<Field,RecursivePoly>& C = resultTasks[j];
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[int][" << intersectDepth << "][" << depth << "] j = " << j << std::endl;
					std::cerr << "[int][" << intersectDepth << "][" << depth << "] C = " << resultTasks[j] << std::endl;
//				#endif
//
//				v = vs[n-i];
//				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[int][" << intersectDepth << "][" << depth << "] x[" << i << "] = " << v << std::endl;
				#endif

//				vIsAResultantsMainVariable = isAMemberOf(v,resultantsMainVars);
//				vIsARegularChainMainVariable = isAMemberOf(v,mainVars);
				if (vIsAResultantsMainVariable) {
					resultantIndex = std::find(resultantsMainVars.begin(),resultantsMainVars.end(),v) - resultantsMainVars.begin();
				}

				if (!vIsAResultantsMainVariable && !vIsARegularChainMainVariable) {
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] x[" << i << "] is in neither the resultants list nor the regular chain" << std::endl;
					#endif

					if (i < n) {
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&intersectStart,&intersectTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), this, C, vs[n-i-1], lazardDecompose, newHeightBound);
//						moreTasks = cleanChain(C,vs[n-i-1],lazardDecompose,newHeightBound);
//						moreTasks = cleanChain(C,vs[n-i-1]);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&intersectStart);
						#endif

						#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] returning from cleanChain." << std::endl;
//							std::cerr << "[int][" << intersectDepth << "][" << depth << "] moreTasks.size = " << moreTasks.size() << std::endl;
//						for (size_t z=0; z<moreTasks.size(); ++z)
//							std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << moreTasks[z] << " added to results." << std::endl;
						#endif

						RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, l) {
							RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, l);
							#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << currCleanChainTask << " added to results." << std::endl;
							#endif
							tasks.push_back(currCleanChainTask);
						}
//						tasks.reserve(tasks.size()+moreTasks.size());
//						tasks.insert(tasks.end(),moreTasks.begin(),moreTasks.end());
					}
					else {
						tasks.emplace_back(std::move(resultTasks[j]));
					}
				}
				else if (!vIsAResultantsMainVariable) {
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] x[" << i << "] is only in the regular chain" << std::endl;
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] adding " << select(v) << " to the current loop iteration chain" << std::endl;
					#endif
	//				currChain += select(v);

					if (i < n) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] constructing C U Tv." << std::endl;
						#endif

						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&intersectStart,&intersectTime);
						#endif
//						resultTasks[j].constructChain(select(v),ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_SQUAREFREE); // TODO: should be able to assume primitive here?
						RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &resultTasks[j], select(v), lazardDecompose, newHeightBound, ASSUME_MAKESCHAIN | ASSUME_REGULAR);

						RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, kk) {
							RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, kk);
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&intersectStart);
							#endif

							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << resultTasks[j] << " constructed." << std::endl;
							#endif

							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&intersectStart,&intersectTime);
							#endif
							RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), this, currConstructChainsTask, vs[n-i-1], lazardDecompose, newHeightBound);
	//						moreTasks = cleanChain(resultTasks[j],vs[n-i-1],lazardDecompose,newHeightBound);
	//						moreTasks = cleanChain(resultTasks[j],vs[n-i-1]);
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&intersectStart);
							#endif

							#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] returning from cleanChain." << std::endl;
	//							std::cerr << "[int][" << intersectDepth << "][" << depth << "] moreTasks.size = " << moreTasks.size() << std::endl;
	//						for (size_t z=0; z<moreTasks.size(); ++z)
	//							std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << moreTasks[z] << " added to results." << std::endl;
							#endif

							RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, l) {
								RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, l);
								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << currCleanChainTask << " added to results." << std::endl;
								#endif
								tasks.push_back(currCleanChainTask);
							}
	//						tasks.reserve(tasks.size()+moreTasks.size());
	//						tasks.insert(tasks.end(),moreTasks.begin(),moreTasks.end());
						}
					}
					else {
						tasks.emplace_back(C+select(v));
					}
				}
				else if (!vIsARegularChainMainVariable) {
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] x[" << i << "] is only in the resultants list" << std::endl;
					#endif
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectStart,&intersectTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_INTFREE_OBJ, intFreeTasks, (&RegularChain<Field, RecursivePoly>::intersectFree), &(C), P[resultantIndex], v, lazardDecompose, newHeightBound);
//					intTasks = C.intersectFree(P[resultantIndex],v,lazardDecompose,newHeightBound);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectStart);
					#endif
					RC_GEN_CONSUMER_LOOP(intFreeTasks, currIntFreeTask, k) {
//					for (k=0; k<intTasks.size(); ++k) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] k = " << k << std::endl;
						#endif
						RC_GEN_CONSUMER_GET_LOOPELEM(intFreeTasks, currIntFreeTask, k);
//						const RegularChain<Field,RecursivePoly>& D = intTasks[k];

						if (i < n) {
							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&intersectStart,&intersectTime);
							#endif
							RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), this, currIntFreeTask, vs[n-i-1], lazardDecompose, newHeightBound);
//							moreTasks = cleanChain(currIntFreeTask,vs[n-i-1],lazardDecompose,newHeightBound);
//							moreTasks = cleanChain(D,vs[n-i-1],lazardDecompose,newHeightBound);
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&intersectStart);
							#endif

							#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] returning from cleanChain." << std::endl;
//							std::cerr << "[int][" << intersectDepth << "][" << depth << "] moreTasks.size = " << moreTasks.size() << std::endl;
//							for (size_t z=0; z<moreTasks.size(); ++z)
//								std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << moreTasks[z] << " added to results." << std::endl;
							#endif

							RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, l) {
								RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, l);
								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << currCleanChainTask << " added to results." << std::endl;
								#endif
								tasks.push_back(currCleanChainTask);
							}
//							tasks.reserve(tasks.size()+moreTasks.size());
//							tasks.insert(tasks.end(),moreTasks.begin(),moreTasks.end());
						}
						else {
							tasks.emplace_back(std::move(currIntFreeTask));
//							tasks.emplace_back(std::move(intTasks[k]));
						}
					}
				}
				else {
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[int][" << intersectDepth << "][" << depth << "] x[" << i << "] is in both the resultants list and the regular chain" << std::endl;
					#endif
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectStart,&intersectTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_INTALG_OBJ, intAlgTasks, (&RegularChain<Field, RecursivePoly>::intersectAlgebraic), &(C), P[resultantIndex], *this, v, S[resultantIndex], lazardDecompose, newHeightBound);
//					intAlgTasks = C.intersectAlgebraic(P[resultantIndex],*this,v,S[resultantIndex],lazardDecompose,newHeightBound);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectStart);
					#endif
					RC_GEN_CONSUMER_LOOP(intAlgTasks, currIntAlgTask, k) {
//					for (k=0; k<intAlgTasks.size(); ++k) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] k = " << k << std::endl;
						#endif
						RC_GEN_CONSUMER_GET_LOOPELEM(intAlgTasks, currIntAlgTask, k);
//						const RegularChain<Field,RecursivePoly>& D = intAlgTasks[k];

						if (i < n) {
							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&intersectStart,&intersectTime);
							#endif
							RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), this, currIntAlgTask, vs[n-i-1], lazardDecompose, newHeightBound);
//							moreTasks = cleanChain(currIntAlgTask,vs[n-i-1],lazardDecompose,newHeightBound);
//							moreTasks = cleanChain(D,vs[n-i-1],lazardDecompose,newHeightBound);
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&intersectStart);
							#endif

							#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[int][" << intersectDepth << "][" << depth << "] returning from cleanChain." << std::endl;
//							std::cerr << "[int][" << intersectDepth << "][" << depth << "] moreTasks.size = " << moreTasks.size() << std::endl;
//							for (size_t z=0; z<moreTasks.size(); ++z)
//								std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << moreTasks[z] << " added to results." << std::endl;
							#endif

							RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, l) {
								RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, l);
								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[int][" << intersectDepth << "][" << depth << "] chain = " << currCleanChainTask << " added to results." << std::endl;
								#endif
								tasks.push_back(currCleanChainTask);
							}
//							tasks.reserve(tasks.size()+moreTasks.size());
//							tasks.insert(tasks.end(),moreTasks.begin(),moreTasks.end());
						}
						else {
//							tasks.emplace_back(std::move(intAlgTasks[k]));
							tasks.emplace_back(std::move(currIntAlgTask));
						}
					}
				}
			}
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[int][" << intersectDepth << "][" << depth << "] end of for loop iteration, current resultTasks:" << std::endl;
				for (auto t:tasks) {
					std::cerr << t << std::endl;
				}
			#endif
			resultTasks = tasks;
			tasks.clear();
			++i;
		}
//		// TODO: Determine how best to reserve adequate space in results, so that we don't need to reallocate after each iteration.
//		results.reserve(results.size()+resultTasks.size());
//		results.insert(results.end(),resultTasks.begin(),resultTasks.end());
		for (auto rt : resultTasks) {
			RC_GEN_PRODUCER_ACCUMULATE(results, rt);
		}
	}

	#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
		std::cerr << "[int][" << intersectDepth << "][" << depth << "] leaving intersect: " << p << ", " << *this << std::endl;
		for (size_t i=0; i<results.size(); ++i) {
			std::cerr << "[int][" << intersectDepth << "][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
		}
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&intersectStart,&intersectTime);
	#endif
	--intersectDepth;
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::cleanChain(const RegularChain<Field,RecursivePoly>& C, const Symbol& v, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::cleanChain(const RegularChain<Field,RecursivePoly>& C, const Symbol& v, bool lazardDecompose, int heightBound) const {
#endif
	++cleanChainDepth;
	++depth;
	unsigned long long int cleanChainStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&cleanChainStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] entering cleanChain..." << v << std::endl;
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_CLEANCH_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_CLEANCH_OBJ, results);
//	std::vector<RegularChain<Field,RecursivePoly>> results;
	RegularChain<Field,RecursivePoly> Tlv;
	this->lower(v,Tlv);
	if (C.isEmpty() || !isAMemberOf(v,this->mainVariables()) || (C.dimension() == Tlv.dimension())) {
		RC_GEN_PRODUCER_ACCUMULATE(results, C);
//		results.emplace_back(C);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] leaving cleanChain without changes: " << v << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&cleanChainStart,&cleanChainTime);
		#endif
		--cleanChainDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] dimension drop!" << std::endl;
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] Tlv.dimension = " << Tlv.dimension() << std::endl;
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] C = " << C << std::endl;
//		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] C.variables:" << std::endl;
//		printVariables(rc.variables(),"vars");
//		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] C.algebraicVariables:" << std::endl;
//		printVariables(rc.mainVariables(),"algVars");
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] C.dimension = " << C.dimension() << std::endl;
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] T = " << *this << std::endl;
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] T.dimension = " << this->dimension() << std::endl;
	#endif

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask);
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularComponents;
	RecursivePoly pRed;

	if (C.isRegularizationTrivial(select(v).initial(),pRed)) {
//		std::cerr << "Trivial regularization in cleanChain!" << std::endl;
//		std::cerr << "p.initial = " << select(v).initial() << std::endl;
//		std::cerr << "pRed = " << pRed << std::endl;
//		std::cerr << "this = " << C << std::endl;
		regularComponents = C.regularizeTrivial(select(v).initial(),pRed);

		for (size_t i=0; i<regularComponents.size(); ++i) {
			if (!regularComponents[i].poly.isZero()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] chain = " << regularComponents[i].chain << " added to results." << std::endl;
				#endif
				RC_GEN_PRODUCER_ACCUMULATE(results, regularComponents[i].chain);
			}
		}
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&cleanChainStart,&cleanChainTime);
		#endif
		RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(C), select(v).initial(), lazardDecompose, heightBound);
	//	regularComponents = C.regularize(select(v).initial(),lazardDecompose,heightBound);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&cleanChainStart);
		#endif

		#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] returning from regularize." << std::endl;
		#endif

		RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, i) {
	//	for (size_t i=0; i<regularComponents.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, i);
			if (!currRegTask.poly.isZero()) {
	//		if (!regularComponents[i].poly.isZero()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] chain = " << currRegTask.chain << " added to results." << std::endl;
				#endif

				RC_GEN_PRODUCER_ACCUMULATE(results, currRegTask.chain);
	//			RC_GEN_PRODUCER_ACCUMULATE(results, regularComponents[i].chain);
	//			results.emplace_back(std::move(regularComponents[i].chain));
			}
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[cc][" << cleanChainDepth << "][" << depth << "] leaving cleanChain: " << v << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&cleanChainStart,&cleanChainTime);
	#endif
	--cleanChainDepth;
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::intersectFree(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::intersectFree(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound) const {
#endif
	++intersectFreeDepth;
	++depth;
	unsigned long long int intersectFreeStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&intersectFreeStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] entering intersectFree: " << p << std::endl;
		std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif
	typedef RegularChain<Field,RecursivePoly> RC_INTFREE_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_INTFREE_OBJ, results);
	typedef RC_INTFREE_OBJ RC_INT_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask2);
	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask);
	typedef RC_INTFREE_OBJ RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);

	std::vector<RegularChain<Field,RecursivePoly>> subtasks,subSubtasks,newResults,moreNewResults;
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> tasks;
	PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> currTask;
	RegularChain<Field,RecursivePoly> currSubtask;
	RecursivePoly r,pTail(p.tail());
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] p.tail = " << pTail << std::endl;
	#endif

	if (p.leadingVariable() != v) {
		std::cerr << "BPAS: error, leading variable of p must be v in intersectFree(p,v)" << std::endl;
		exit(1);
	}

	if (isRegularizationTrivial(p.initial(),r)) {
//		std::cerr << "Trivial regularization in intersectFree!" << std::endl;
		tasks = regularizeTrivial(p.initial(),r);

		for (size_t i=0; i<tasks.size(); ++i) {
			#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
				std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] i = " << i << std::endl;
			#endif
			const RegularChain<Field,RecursivePoly>& D = tasks[i].chain;
			if (tasks[i].poly.isZero()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] f = 0" << std::endl;
				#endif

				if (D.isIntersectionTrivial(pTail,r)) {
	//				std::cerr << "First trivial intersection in intersectFree!" << std::endl;
					newResults = D.intersectTrivial(r);
					for (size_t z=0; z<newResults.size(); ++z) {
						RC_GEN_PRODUCER_ACCUMULATE(results, newResults[z]);
					}
				}
				else {
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), pTail, lazardDecompose, heightBound);
		//			subtasks = D.intersect(pTail,lazardDecompose,heightBound);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectFreeStart);
					#endif

					#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from intersect." << std::endl;
		//			for (size_t z=0; z<subtasks.size(); ++z)
		//				std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subtasks[z] << " added to results." << std::endl;
					#endif

		//			results.reserve(results.size()+subtasks.size());
		//			results.insert(results.end(),subtasks.begin(),subtasks.end());
					RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, z) {
		//			for (auto st : subtasks) {
						RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, z);
						RC_GEN_PRODUCER_ACCUMULATE(results, currIntTask);
					}
				}
			}
			else {
				if (D.numberOfAlgebraicVariables() < heightBound) {
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] constructing D U p." << std::endl;
					#endif

					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &(D), p, lazardDecompose, heightBound, ASSUME_MAKESCHAIN | ASSUME_REGULAR);
					RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, j) {
						RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, j);
						subtasks.push_back(currConstructChainsTask);
					}
	//				subtasks = D.constructChains(p,lazardDecompose,heightBound,ASSUME_MAKESCHAIN | ASSUME_REGULAR);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectFreeStart);
					#endif

					#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from constructChains." << std::endl;
					for (size_t z=0; z<subtasks.size(); ++z)
						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subtasks[z] << " added to results." << std::endl;
					#endif

					for (auto st : subtasks) {
						RC_GEN_PRODUCER_ACCUMULATE(results, st);
					}
	//				results.insert(results.end(),subtasks.begin(),subtasks.end());

					if (D.isIntersectionTrivial(p.initial(),r)) {
	//					std::cerr << "Second trivial intersection in intersectFree!" << std::endl;
						newResults = D.intersectTrivial(r);

						for (size_t j=0; j<newResults.size(); ++j) {
							const RegularChain<Field,RecursivePoly>& E = newResults[j];
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] E = " << newResults[j] << std::endl;
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] p.tail = " << pTail << std::endl;
							#endif

							if (E.isIntersectionTrivial(pTail,r)) {
	//							std::cerr << "Third trivial intersection in intersectFree!" << std::endl;
								moreNewResults = E.intersectTrivial(r);
								for (size_t k=0; k<newResults.size(); ++k) {
									RC_GEN_PRODUCER_ACCUMULATE(results, moreNewResults[k]);
								}
							}
							else {
								#ifdef REGULARCHAIN_PROFILING
									stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
								#endif
								RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks2, (&RegularChain<Field, RecursivePoly>::_intersect), &(E), pTail, lazardDecompose, heightBound);
			//					subSubtasks = E.intersect(pTail,lazardDecompose,heightBound);
								#ifdef REGULARCHAIN_PROFILING
									startTimer(&intersectFreeStart);
								#endif

								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from intersect." << std::endl;
			//					for (size_t z=0; z<subSubtasks.size(); ++z)
			//						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subSubtasks[z] << " added to results." << std::endl;
								#endif

								RC_GEN_CONSUMER_LOOP(intTasks2, currIntTask2, k) {
			//					for (auto sst : subSubtasks) {
									RC_GEN_CONSUMER_GET_LOOPELEM(intTasks2, currIntTask2, k);
									RC_GEN_PRODUCER_ACCUMULATE(results, currIntTask2);
								}
			//					results.reserve(results.size()+subSubtasks.size());
			//					results.insert(results.end(),subSubtasks.begin(),subSubtasks.end());
							}
						}
					}
					else {
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), p.initial(), lazardDecompose, heightBound);
		//				subtasks = D.intersect(p.initial(),lazardDecompose,heightBound);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&intersectFreeStart);
						#endif

						RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, j) {
		//				for (size_t j=0; j<subtasks.size(); ++j) {
							RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
							const RegularChain<Field,RecursivePoly>& E = currIntTask;
		//					const RegularChain<Field,RecursivePoly>& E = subtasks[j];
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] E = " << currIntTask << std::endl;
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] p.tail = " << pTail << std::endl;
							#endif

							if (E.isIntersectionTrivial(pTail,r)) {
	//							std::cerr << "Third trivial intersection in intersectFree!" << std::endl;
								newResults = E.intersectTrivial(r);
								for (size_t k=0; k<newResults.size(); ++k) {
									RC_GEN_PRODUCER_ACCUMULATE(results, newResults[k]);
								}
							}
							else {
								#ifdef REGULARCHAIN_PROFILING
									stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
								#endif
								RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks2, (&RegularChain<Field, RecursivePoly>::_intersect), &(E), pTail, lazardDecompose, heightBound);
			//					subSubtasks = E.intersect(pTail,lazardDecompose,heightBound);
								#ifdef REGULARCHAIN_PROFILING
									startTimer(&intersectFreeStart);
								#endif

								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from intersect." << std::endl;
			//					for (size_t z=0; z<subSubtasks.size(); ++z)
			//						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subSubtasks[z] << " added to results." << std::endl;
								#endif

								RC_GEN_CONSUMER_LOOP(intTasks2, currIntTask2, k) {
			//					for (auto sst : subSubtasks) {
									RC_GEN_CONSUMER_GET_LOOPELEM(intTasks2, currIntTask2, k);
									RC_GEN_PRODUCER_ACCUMULATE(results, currIntTask2);
								}
			//					results.reserve(results.size()+subSubtasks.size());
			//					results.insert(results.end(),subSubtasks.begin(),subSubtasks.end());
							}
						}
					}
				}
			}
		}
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
		#endif
		RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), this, p.initial(), lazardDecompose, heightBound);
	//	tasks = regularize(p.initial(),lazardDecompose,heightBound);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&intersectFreeStart);
		#endif

		RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, i) {
	//	for (size_t i=0; i<tasks.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, i);
			#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
				std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] i = " << i << std::endl;
			#endif
			const RegularChain<Field,RecursivePoly>& D = currRegTask.chain;
	//		const RegularChain<Field,RecursivePoly>& D = tasks[i].chain;
			if (currRegTask.poly.isZero()) {
	//		if (tasks[i].poly.isZero()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] f = 0" << std::endl;
				#endif

				if (D.isIntersectionTrivial(pTail,r)) {
	//				std::cerr << "First trivial intersection in intersectFree!" << std::endl;
					newResults = D.intersectTrivial(r);
					for (size_t z=0; z<newResults.size(); ++z) {
						RC_GEN_PRODUCER_ACCUMULATE(results, newResults[z]);
					}
				}
				else {
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), pTail, lazardDecompose, heightBound);
		//			subtasks = D.intersect(pTail,lazardDecompose,heightBound);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectFreeStart);
					#endif

					#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from intersect." << std::endl;
		//			for (size_t z=0; z<subtasks.size(); ++z)
		//				std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subtasks[z] << " added to results." << std::endl;
					#endif

		//			results.reserve(results.size()+subtasks.size());
		//			results.insert(results.end(),subtasks.begin(),subtasks.end());
					RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, z) {
		//			for (auto st : subtasks) {
						RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, z);
						RC_GEN_PRODUCER_ACCUMULATE(results, currIntTask);
					}
				}
			}
			else {
				if (D.numberOfAlgebraicVariables() < heightBound) {
					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] constructing D U p." << std::endl;
					#endif

					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &(D), p, lazardDecompose, heightBound, ASSUME_MAKESCHAIN | ASSUME_REGULAR);
					RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, j) {
						RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, j);
						subtasks.push_back(currConstructChainsTask);
					}
	//				subtasks = D.constructChains(p,lazardDecompose,heightBound,ASSUME_MAKESCHAIN | ASSUME_REGULAR);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectFreeStart);
					#endif

					#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from constructChains." << std::endl;
					for (size_t z=0; z<subtasks.size(); ++z)
						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subtasks[z] << " added to results." << std::endl;
					#endif

					for (auto st : subtasks) {
						RC_GEN_PRODUCER_ACCUMULATE(results, st);
					}
	//				results.insert(results.end(),subtasks.begin(),subtasks.end());

					if (D.isIntersectionTrivial(p.initial(),r)) {
	//					std::cerr << "Second trivial intersection in intersectFree!" << std::endl;
						newResults = D.intersectTrivial(r);

						for (size_t j=0; j<newResults.size(); ++j) {
							const RegularChain<Field,RecursivePoly>& E = newResults[j];
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] E = " << newResults[j] << std::endl;
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] p.tail = " << pTail << std::endl;
							#endif

							if (E.isIntersectionTrivial(pTail,r)) {
	//							std::cerr << "Third trivial intersection in intersectFree!" << std::endl;
								moreNewResults = E.intersectTrivial(r);
								for (size_t k=0; k<newResults.size(); ++k) {
									RC_GEN_PRODUCER_ACCUMULATE(results, moreNewResults[k]);
								}
							}
							else {
								#ifdef REGULARCHAIN_PROFILING
									stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
								#endif
								RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks2, (&RegularChain<Field, RecursivePoly>::_intersect), &(E), pTail, lazardDecompose, heightBound);
			//					subSubtasks = E.intersect(pTail,lazardDecompose,heightBound);
								#ifdef REGULARCHAIN_PROFILING
									startTimer(&intersectFreeStart);
								#endif

								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from intersect." << std::endl;
			//					for (size_t z=0; z<subSubtasks.size(); ++z)
			//						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subSubtasks[z] << " added to results." << std::endl;
								#endif

								RC_GEN_CONSUMER_LOOP(intTasks2, currIntTask2, k) {
			//					for (auto sst : subSubtasks) {
									RC_GEN_CONSUMER_GET_LOOPELEM(intTasks2, currIntTask2, k);
									RC_GEN_PRODUCER_ACCUMULATE(results, currIntTask2);
								}
			//					results.reserve(results.size()+subSubtasks.size());
			//					results.insert(results.end(),subSubtasks.begin(),subSubtasks.end());
							}
						}
					}
					else {
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), p.initial(), lazardDecompose, heightBound);
		//				subtasks = D.intersect(p.initial(),lazardDecompose,heightBound);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&intersectFreeStart);
						#endif

						RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, j) {
		//				for (size_t j=0; j<subtasks.size(); ++j) {
							RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
							const RegularChain<Field,RecursivePoly>& E = currIntTask;
		//					const RegularChain<Field,RecursivePoly>& E = subtasks[j];
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] E = " << currIntTask << std::endl;
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] p.tail = " << pTail << std::endl;
							#endif

							if (E.isIntersectionTrivial(pTail,r)) {
	//							std::cerr << "Third trivial intersection in intersectFree!" << std::endl;
								newResults = E.intersectTrivial(r);
								for (size_t k=0; k<newResults.size(); ++k) {
									RC_GEN_PRODUCER_ACCUMULATE(results, newResults[k]);
								}
							}
							else {
								#ifdef REGULARCHAIN_PROFILING
									stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
								#endif
								RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks2, (&RegularChain<Field, RecursivePoly>::_intersect), &(E), pTail, lazardDecompose, heightBound);
			//					subSubtasks = E.intersect(pTail,lazardDecompose,heightBound);
								#ifdef REGULARCHAIN_PROFILING
									startTimer(&intersectFreeStart);
								#endif

								#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] returning from intersect." << std::endl;
			//					for (size_t z=0; z<subSubtasks.size(); ++z)
			//						std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] chain = " << subSubtasks[z] << " added to results." << std::endl;
								#endif

								RC_GEN_CONSUMER_LOOP(intTasks2, currIntTask2, k) {
			//					for (auto sst : subSubtasks) {
									RC_GEN_CONSUMER_GET_LOOPELEM(intTasks2, currIntTask2, k);
									RC_GEN_PRODUCER_ACCUMULATE(results, currIntTask2);
								}
			//					results.reserve(results.size()+subSubtasks.size());
			//					results.insert(results.end(),subSubtasks.begin(),subSubtasks.end());
							}
						}
					}
				}
			}
		}
	}
	#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
		std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] leaving intersectFree: " << p << std::endl;
		for (size_t i=0; i<results.size(); ++i)
			std::cerr << "[iF][" << intersectFreeDepth << "][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&intersectFreeStart,&intersectFreeTime);
	#endif
	--intersectFreeDepth;
	--depth;
//	return results;
	RC_GEN_PRODUCER_COMPLETE(results);
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::intersectAlgebraic(const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& T, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::intersectAlgebraic(const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& T, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound) const {
#endif
	++intersectAlgebraicDepth;
	++depth;
	unsigned long long int intersectAlgebraicStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&intersectAlgebraicStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] entering intersectAlgebraic:" << p << std::endl;
		std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif

	typedef RegularChain<Field,RecursivePoly> RC_INTALG_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_INTALG_OBJ, results);
	typedef RC_INTALG_OBJ RC_CLEANCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CLEANCH_OBJ, currCleanChainTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INTALG_OBJ, currIntAlgTask);
	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGGCD_OBJ, currRegGCDTask);
	typedef RC_INTALG_OBJ RC_INT_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);
	typedef RC_INTALG_OBJ RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);

//	std::vector<RegularChain<Field,RecursivePoly>> results,subtasks,subSubtasks,subSubSubtasks;
	std::vector<RegularChain<Field,RecursivePoly>> subtasks,subSubtasks,subSubSubtasks,newResults;
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> tasks;
	RecursivePoly r;

	if ((p.leadingVariable() != v) || (!isAMemberOf(v,T.mainVariables()))) {
		std::cerr << "BPAS: error, leading variable of p must be v in intersectAlgebraic(p,v)" << std::endl;
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] p.leadingVariable = " << p.leadingVariable() << ", v = " << v << std::endl;
			std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] isAMemberOf(v,T.mainVariables()) = " << isAMemberOf(v,T.mainVariables()) << std::endl;
		#endif
		exit(1);
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] p = " << p << std::endl;
		std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] T = " << T << std::endl;
		std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] *this = " << *this << std::endl;
	#endif

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
	#endif
	RC_GEN_CONSUMER_INIT(RC_REGGCD_OBJ, intRegGCDTasks, (&RegularChain<Field, RecursivePoly>::_regularGCD), this, p, T.select(v), v, src, lazardDecompose, heightBound-1);
//	tasks = regularGCD(p,T.select(v),v,src,lazardDecompose,heightBound-1);
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&intersectAlgebraicStart);
	#endif

	RC_GEN_CONSUMER_LOOP(intRegGCDTasks, currRegGCDTask, i) {
//	for (size_t i=0; i<tasks.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(intRegGCDTasks, currRegGCDTask, i);
		const RecursivePoly& g = currRegGCDTask.poly;
		const RegularChain<Field,RecursivePoly>& D = currRegGCDTask.chain;
//		const RecursivePoly& g = tasks[i].poly;
//		const RegularChain<Field,RecursivePoly>& D = tasks[i].chain;
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] g = " << currRegGCDTask.poly << std::endl;
			std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] D = " << currRegGCDTask.chain << std::endl;
			std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] *this = " << *this << std::endl;
		#endif

		if (D.dimension() < this->dimension()) {

			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
			#endif
			RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), &(T), D, v, lazardDecompose, heightBound-1);
//			subtasks = T.cleanChain(D,v,lazardDecompose,heightBound-1);
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&intersectAlgebraicStart);
			#endif

			RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, j) {
//			for (size_t j=0; j<subtasks.size(); ++j) {
				RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, j);
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] result from cleanChain call:" << std::endl;
					std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] E = " << currCleanChainTask << std::endl;
					std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] *this = " << *this << std::endl;
				#endif

				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
				#endif
				RC_GEN_CONSUMER_INIT(RC_INTALG_OBJ, intAlgTasks, (&RegularChain<Field, RecursivePoly>::intersectAlgebraic), &(currCleanChainTask), p, T, v, src, lazardDecompose, heightBound);
//				subSubtasks = currCleanChainTask.intersectAlgebraic(p,T,v,src,lazardDecompose,heightBound);
				#ifdef REGULARCHAIN_PROFILING
					startTimer(&intersectAlgebraicStart);
				#endif

				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] returning from recursive iA call." << std::endl;
					for (size_t z=0; z<subSubtasks.size(); ++z)
						std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] chain = " << subSubtasks[z] << " added to results." << std::endl;
				#endif

//				for (auto sst : subSubtasks) {
				RC_GEN_CONSUMER_LOOP(intAlgTasks, currIntAlgTask, k) {
					RC_GEN_CONSUMER_GET_LOOPELEM(intAlgTasks, currIntAlgTask, k);
					RC_GEN_PRODUCER_ACCUMULATE(results, currIntAlgTask);
				}
//				results.reserve(results.size()+subSubtasks.size());
//				results.insert(results.end(),subSubtasks.begin(),subSubtasks.end());
			}
		}
		else {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] constructing D U g." << std::endl;
			#endif

			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
			#endif
			RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &(D), g, lazardDecompose, heightBound, ASSUME_MAKESCHAIN | ASSUME_REGULAR);
			RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, j) {
				RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, j);
				subtasks.push_back(currConstructChainsTask);
			}
//			subtasks = D.constructChains(g,lazardDecompose,heightBound,ASSUME_MAKESCHAIN | ASSUME_REGULAR);
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&intersectAlgebraicStart);
			#endif

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] returning from constructChains." << std::endl;
				for (size_t z=0; z<subtasks.size(); ++z)
					std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] chain = " << subtasks[z] << " added to results." << std::endl;
			#endif

			for (auto st : subtasks) {
				RC_GEN_PRODUCER_ACCUMULATE(results, st);
			}
//			results.insert(results.end(),subtasks.begin(),subtasks.end());

			if (D.numberOfAlgebraicVariables() < heightBound-1) {

				if (D.isIntersectionTrivial(g.initial(),r)) {
//					std::cerr << "Trivial intersection in intersectAlgebraic!" << std::endl;
					newResults = D.intersectTrivial(r);

					for (size_t j=0; j<newResults.size(); ++j) {
						const RegularChain<Field,RecursivePoly>& E = newResults[j];
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] result from intersect call:" << std::endl;
							std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] E = " << newResults[j] << std::endl;
							std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] *this = " << *this << std::endl;
						#endif

						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), &(T), E, v, lazardDecompose, heightBound-1);
	//					subSubtasks = T.cleanChain(E,v,lazardDecompose,heightBound-1);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&intersectAlgebraicStart);
						#endif

						RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, k) {
	//					for (size_t k=0; k<subSubtasks.size(); ++k) {
							RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, k);
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] result from cleanChain call:" << std::endl;
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] F = " << currCleanChainTask << std::endl;
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] T = " << T << std::endl;
							#endif

							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
							#endif
							RC_GEN_CONSUMER_INIT(RC_INTALG_OBJ, intAlgTasks, (&RegularChain<Field, RecursivePoly>::intersectAlgebraic), &(currCleanChainTask), p, T, v, src, lazardDecompose, heightBound);
	//						subSubSubtasks = currCleanChainTask.intersectAlgebraic(p,T,v,src,lazardDecompose,heightBound);
	//						subSubSubtasks = subSubtasks[k].intersectAlgebraic(p,T,v,src,lazardDecompose,heightBound);
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&intersectAlgebraicStart);
							#endif

							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] returning from recursive iA call." << std::endl;
								for (size_t z=0; z<subSubSubtasks.size(); ++z)
									std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] chain = " << subSubSubtasks[z] << " added to results." << std::endl;
							#endif

	//						for (auto ssst : subSubSubtasks) {
							RC_GEN_CONSUMER_LOOP(intAlgTasks, currIntAlgTask, l) {
								RC_GEN_CONSUMER_GET_LOOPELEM(intAlgTasks, currIntAlgTask, l);
								RC_GEN_PRODUCER_ACCUMULATE(results, currIntAlgTask);
							}
	//						results.reserve(results.size()+subSubSubtasks.size());
	//						results.insert(results.end(),subSubSubtasks.begin(),subSubSubtasks.end());
						}
					}
				}
				else {
					#ifdef REGULARCHAIN_PROFILING
						stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
					#endif
					RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), g.initial(), lazardDecompose, heightBound-1);
	//				subtasks = D.intersect(g.initial(),lazardDecompose,heightBound-1);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&intersectAlgebraicStart);
					#endif

					RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, j) {
	//				for (size_t j=0; j<subtasks.size(); ++j) {
						RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, j);
						const RegularChain<Field,RecursivePoly>& E = currIntTask;
	//					const RegularChain<Field,RecursivePoly>& E = subtasks[j];
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] result from intersect call:" << std::endl;
							std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] E = " << currIntTask << std::endl;
							std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] *this = " << *this << std::endl;
						#endif

						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_CLEANCH_OBJ, cleanChainTasks, (&RegularChain<Field, RecursivePoly>::cleanChain), &(T), E, v, lazardDecompose, heightBound-1);
	//					subSubtasks = T.cleanChain(E,v,lazardDecompose,heightBound-1);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&intersectAlgebraicStart);
						#endif

						RC_GEN_CONSUMER_LOOP(cleanChainTasks, currCleanChainTask, k) {
	//					for (size_t k=0; k<subSubtasks.size(); ++k) {
							RC_GEN_CONSUMER_GET_LOOPELEM(cleanChainTasks, currCleanChainTask, k);
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] result from cleanChain call:" << std::endl;
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] F = " << currCleanChainTask << std::endl;
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] T = " << T << std::endl;
							#endif

							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
							#endif
							RC_GEN_CONSUMER_INIT(RC_INTALG_OBJ, intAlgTasks, (&RegularChain<Field, RecursivePoly>::intersectAlgebraic), &(currCleanChainTask), p, T, v, src, lazardDecompose, heightBound);
	//						subSubSubtasks = currCleanChainTask.intersectAlgebraic(p,T,v,src,lazardDecompose,heightBound);
	//						subSubSubtasks = subSubtasks[k].intersectAlgebraic(p,T,v,src,lazardDecompose,heightBound);
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&intersectAlgebraicStart);
							#endif

							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] returning from recursive iA call." << std::endl;
								for (size_t z=0; z<subSubSubtasks.size(); ++z)
									std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] chain = " << subSubSubtasks[z] << " added to results." << std::endl;
							#endif

	//						for (auto ssst : subSubSubtasks) {
							RC_GEN_CONSUMER_LOOP(intAlgTasks, currIntAlgTask, l) {
								RC_GEN_CONSUMER_GET_LOOPELEM(intAlgTasks, currIntAlgTask, l);
								RC_GEN_PRODUCER_ACCUMULATE(results, currIntAlgTask);
							}
	//						results.reserve(results.size()+subSubSubtasks.size());
	//						results.insert(results.end(),subSubSubtasks.begin(),subSubSubtasks.end());
						}
					}
				}
			}
		}
	}
	#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
		std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] leaving intersectAlgebraic " << p << std::endl;
		for (size_t i=0; i<results.size(); ++i)
			std::cerr << "[iA][" << intersectAlgebraicDepth << "][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&intersectAlgebraicStart,&intersectAlgebraicTime);
	#endif
	--intersectAlgebraicDepth;
	--depth;
//	return results;
	RC_GEN_PRODUCER_COMPLETE(results);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::cutChain(const RegularChain<Field,RecursivePoly>& T, const Symbol& v, RegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const {
		T.lower(v,Tlv);
		T.upper(v,Tgv);
		Tv = T.select(v);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::cutChain(const Symbol& v, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const {
		Tv = select(v);
		upper(v,Tgv);
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::cutChain(const Symbol& v, RegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv) const {
		lower(v,Tlv);
		Tv = select(v);
}

// template <class Field, class RecursivePoly>
// bool RegularChain<Field,RecursivePoly>::isRegular(const RecursivePoly& p) const {
// 	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results;
// 	#ifdef REGULARCHAIN_DEBUG
// 		std::cerr << "performing regularity test in isRegular." << std::endl;
// 	#endif
// 	results = regularize(p);
// 	if (results.size() != 1 || results[0].poly.isZero())
// 		return false;
// 	else
// 		return true;
// }


template <class Field, class RecursivePoly>
bool RegularChain<Field, RecursivePoly>::isRegular(const RecursivePoly& p) const {

	if (p.isZero()) {
		return false;
	}
	if (this->isEmpty() || p.isConstant() != 0) {
		return true;
	}

	// std::cerr << "[isRegular] p : " << p << std::endl;
	// std::cerr << "[isRegular] prime: " << this->isSaturatedIdealPrime() << std::endl;
	// std::cerr << "[isRegular] rc: " << *this << std::endl;

	if (this->isSaturatedIdealPrime()) {
		return !isInSaturatedIdealMinimal(p);
	}

	// exit(27);
	return !isIteratedResultantZeroModular(p);

	bool lazyResultant = 1;
	SubResultantChain<RecursivePoly,RecursivePoly> src;
	std::vector<RecursivePoly> lf = GCDFreeFactorization(p);
	for(size_t i = 0; i < lf.size(); ++i) {
		if (canComputeInDimensionZero(lf[i])) {
			lf[i] = this->reduce(lf[i]);
			if (lf[i].isZero()) {
				return false;
			}
			if (lf[i].isConstant() == 0) {
				lf[i] = lf[i].primitivePart();
			}
		}
		if (lf[i].isConstant() == 0) {
			std::vector<Symbol> vars = this->mainVariables();
			for (size_t v = 0; v < vars.size(); ++v) {
				if (lf[i].degree(v) > 0) {
					src = SubResultantChain<RecursivePoly,RecursivePoly>(lf[i],select(v),v);
					lf[i] = src.resultant(lazyResultant).primitivePart();
					if (lf[i].isZero()) {
						return false;
					}
					lf[i].setRingVariables(this->allVariables());
				}
			}
		}
	}
	return true;
}


/**
 * NOTE this has not be thoroughly tested. It's rarely ever called.
 */
template <class Field, class RecursivePoly>
bool RegularChain<Field, RecursivePoly>::isIteratedResultantZeroModular(const RecursivePoly& p) const {

	RecursivePoly pp = p.primitivePart();
	Prime_ptr pr;
	bool good = false;
	for (int i = 0 ; ; ++i) {
		pr = prime32_ptr[i];
		if (Integer(pp.leadingCoefficient()).divisible(pr.prime)) {
			continue;
		}
		for (int j; j < this->set.size(); ++j) {
			if (Integer(set[j].leadingCoefficient()).divisible(pr.prime)) {
				continue;
			}
		}
		break;
	}

	Integer prI(pr.prime);
	bool lazyResultant = 1;
	SubResultantChain<RecursivePoly,RecursivePoly> src;
	std::vector<RecursivePoly> tempFacts;
	std::deque<RecursivePoly> factTasks;

	Symbol v;
	RecursivePoly r;
	RecursivePoly curP, Tv;

	curP = std::move(pp);
	curP.setRingVariables(this->allVariables());
	factTasks.emplace_back(std::move(curP));

	while (factTasks.size() > 0) {
		curP = std::move(factTasks.front());
		factTasks.pop_front();
		std::vector<RecursivePoly> lf = GCDFreeFactorization(curP);
		for (size_t i = 0; i < lf.size(); ++i) {
			v = lf[i].leadingVariable();
			Tv = select(v);
			if (Tv.isZero()) {
				continue;
			}

			src = SubResultantChain<RecursivePoly,RecursivePoly>(lf[i], select(v), v);
			r = src.resultant(lazyResultant);
			r %= prI;

			if (r.isZero()) {
				return 1;
			}
			if (r.isConstant()) {
				continue;
			}

			tempFacts = GCDFreeFactorization(r);
			factTasks.insert(factTasks.end(), std::make_move_iterator(tempFacts.begin()), std::make_move_iterator(tempFacts.end()));
		}
	}

	return 0;
}

//template <class Field, class RecursivePoly>
//std::vector<BoolChainPair<RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::areFactorsRegular(const std::vector<RecursivePoly>& factors) const {
//	struct Task {
//		std::vector<RecursivePoly> lf;
//		RegularChain<Field,RecursivePoly> T;
//	};
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularComponents;
//	std::vector<BoolChainPair<RegularChain<Field,RecursivePoly>>> results;
//	std::vector<Task> tasks;
//	Task task,t;
//	RecursivePoly f;
//	task = {factors,*this};
//	tasks.push_back(t);
//	while (!tasks.empty()) {
//		task = tasks.back();
//		tasks.pop_back();
//		if (task.lf.empty()) {
//			results.push_back(BoolChainPair<RegularChain<Field,RecursivePoly>>(true,task.T));
//		}
//		else {
//			f = task.lf.back();
//			task.lf.pop_back();
//			regularComponents = task.T.regularize(f);
//			for (size_t i=0; i<regularComponents.size(); ++i) {
//				if (regularComponents[i].poly.isZero()) {
//					results.push_back(BoolChainPair<RegularChain<Field,RecursivePoly>>(false,regularComponents[i].chain));
//				}
//				else {
//					if (task.lf.empty()) {
//						results.push_back(BoolChainPair<RegularChain<Field,RecursivePoly>>(true,regularComponents[i].chain));
//					}
//					else {
//						Task t = {task.lf,regularComponents[i].chain};
//						tasks.push_back(t);
//					}
//				}
//			}
//		}
//	}
//	return results;
//}

// ORIGINAL VERSION //
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::regularizeList(const std::vector<RecursivePoly>& knownRegularIn, const std::vector<RecursivePoly>& unknownIfRegularIn, std::vector<RegularChain<Field,RecursivePoly>>& singularComponents, std::vector<RegularChain<Field,RecursivePoly>>& regularComponents, bool lazardDecompose, int heightBound) const {
	unsigned long long int regularizeStart = 0;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] entering regularizeList:" << std::endl;
		std::cerr << "[regList] heightBound = " << heightBound << std::endl;
	#endif

	singularComponents.clear();
	regularComponents.clear();
	if (unknownIfRegularIn.empty()) {
		regularComponents.emplace_back(*this);
		return;
	}

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGSINGLE_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGSINGLE_OBJ, currRegSingleTask);
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularizeResults;
	std::vector<RegularChain<Field,RecursivePoly>> singularComponents2,regularComponents2;
	std::vector<RecursivePoly> knownRegular,unknownIfRegular,knownRegular2,unknownIfRegular2;
	RecursivePoly p;
	knownRegular = knownRegularIn;
	unknownIfRegular = unknownIfRegularIn;
	p = unknownIfRegular.back();
	unknownIfRegular.pop_back();

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] calling regularizeSingle:" << std::endl;
	#endif

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
	RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), this, p, lazardDecompose, heightBound);
//	regularizeResults = regularizeSingle(p,lazardDecompose,heightBound);
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] returning to regularizeList from regularizeSingle:" << std::endl;
	#endif
	RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, i) {
//	for (size_t i=0; i<regularizeResults.size(); ++i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, i);
		if (currRegSingleTask.poly.isZero()) {
//		if (regularizeResults[i].poly.isZero()) {
			singularComponents.emplace_back(std::move(currRegSingleTask.chain));
//			singularComponents.emplace_back(std::move(regularizeResults[i].chain));
			continue;
		}
		else if (lazardDecompose && currRegSingleTask.chain.dimension() < this->dimension()) {
//		else if (lazardDecompose && regularizeResults[i].chain.dimension() < this->dimension()) {
			knownRegular2.emplace_back(p);
			unknownIfRegular2 = knownRegular;
			unknownIfRegular2.insert(unknownIfRegular2.begin(),unknownIfRegular.begin(),unknownIfRegular.end());
		}
		else {
			knownRegular2 = knownRegular;
			knownRegular2.emplace_back(p);
			unknownIfRegular2 = unknownIfRegular;
		}
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regList] recursive call to regularizeList:" << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		currRegSingleTask.chain.regularizeList(knownRegular2,unknownIfRegular2,singularComponents2,regularComponents2,lazardDecompose, heightBound);
//		regularizeResults[i].chain.regularizeList(knownRegular2,unknownIfRegular2,singularComponents2,regularComponents2,lazardDecompose, heightBound);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&regularizeStart);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regList] returning to regularizeList from recursive call:" << std::endl;
		#endif
		singularComponents.insert(singularComponents.begin(),singularComponents2.begin(),singularComponents2.end());
		regularComponents.insert(regularComponents.begin(),regularComponents2.begin(),regularComponents2.end());
		knownRegular2.clear();
		unknownIfRegular2.clear();
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] leaving regularizeList:" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
}

#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::regularizeListStream(const RecursivePoly& P, const std::vector<RecursivePoly>& knownRegularIn, const std::vector<RecursivePoly>& unknownIfRegularIn, bool lazardDecompose, int heightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::regularizeListStream(const RecursivePoly& P, const std::vector<RecursivePoly>& knownRegularIn, const std::vector<RecursivePoly>& unknownIfRegularIn, bool lazardDecompose, int heightBound) const {
#endif
	unsigned long long int regularizeStart = 0;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] entering regularizeListStream:" << std::endl;
		std::cerr << "[regList] heightBound = " << heightBound << std::endl;
	#endif

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGLIST_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REGLIST_OBJ, results);
	RC_REGLIST_OBJ result;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGLIST_OBJ, currRegListTask);
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results,results2;
//	singularComponents.clear();
//	regularComponents.clear();
	if (unknownIfRegularIn.empty()) {
		result = RC_REGLIST_OBJ(P,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(P,*this));
//		return results;
		RC_GEN_PRODUCER_COMPLETE(results);
	}

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGSINGLE_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGSINGLE_OBJ, currRegSingleTask);
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularizeResults;
	std::vector<RecursivePoly> knownRegular,unknownIfRegular,knownRegular2,unknownIfRegular2;
	RecursivePoly p;
	knownRegular = knownRegularIn;
	unknownIfRegular = unknownIfRegularIn;
	p = unknownIfRegular.back();
	unknownIfRegular.pop_back();
	RecursivePoly zero,pRed;
	zero.zero();

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] calling regularizeSingle:" << std::endl;
	#endif

	if (isRegularizationTrivial(p,pRed)) {
//		std::cerr << "Trivial regularization in regularizeListStream!" << std::endl;
		regularizeResults = regularizeTrivial(p,pRed);

		for (size_t i=0; i<regularizeResults.size(); ++i) {
			RC_GEN_PRODUCER_ACCUMULATE(results, regularizeResults[i]);
		}
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), this, p, lazardDecompose, heightBound);
	//	regularizeResults = regularizeSingle(p,lazardDecompose,heightBound);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&regularizeStart);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regList] returning to regularizeListStream from regularizeSingle:" << std::endl;
		#endif
		RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, i) {
	//	for (size_t i=0; i<regularizeResults.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, i);
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[regList] regularizeSingleTask[] = " << currRegSingleTask << std::endl;
			#endif
			if (currRegSingleTask.poly.isZero()) {
	//		if (regularizeResults[i].poly.isZero()) {
				result = RC_REGLIST_OBJ(zero,std::move(currRegSingleTask.chain));
				RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//			results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,std::move(currRegSingleTask.chain)));
				continue;
			}
			else if (lazardDecompose && currRegSingleTask.chain.dimension() < this->dimension()) {
	//		else if (lazardDecompose && regularizeResults[i].chain.dimension() < this->dimension()) {
				knownRegular2.emplace_back(p);
				unknownIfRegular2 = knownRegular;
				unknownIfRegular2.insert(unknownIfRegular2.begin(),unknownIfRegular.begin(),unknownIfRegular.end());
			}
			else {
				knownRegular2 = knownRegular;
				knownRegular2.emplace_back(p);
				unknownIfRegular2 = unknownIfRegular;
			}
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[regList] recursive call to regularizeListStream:" << std::endl;
			#endif
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&regularizeStart,&regularizeTime);
			#endif
			RC_GEN_CONSUMER_INIT(RC_REGLIST_OBJ, regListTasks, (&RegularChain<Field, RecursivePoly>::regularizeListStream), &(currRegSingleTask.chain), P, knownRegular2, unknownIfRegular2, lazardDecompose, heightBound);
	//		results2 = currRegSingleTask.chain.regularizeList(p, knownRegular2,unknownIfRegular2,lazardDecompose, heightBound);
	//		regularizeResults[i].chain.regularizeList(knownRegular2,unknownIfRegular2,singularComponents2,regularComponents2,lazardDecompose, heightBound);
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&regularizeStart);
			#endif

			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[regList] returning to regularizeListStream from recursive call:" << std::endl;
			#endif
			RC_GEN_CONSUMER_LOOP(regListTasks, currRegListTask, j) {
	//		for (auto pcp : results2) {
				RC_GEN_CONSUMER_GET_LOOPELEM(regListTasks, currRegListTask, j);
				RC_GEN_PRODUCER_ACCUMULATE(results, currRegListTask);
	//			results.emplace_back(std::move(pcp));
			}
	//		singularComponents.insert(singularComponents.begin(),singularComponents2.begin(),singularComponents2.end());
	//		regularComponents.insert(regularComponents.begin(),regularComponents2.begin(),regularComponents2.end());
			knownRegular2.clear();
			unknownIfRegular2.clear();
		}
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regList] leaving regularizeListStream:" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
//	return results;
	RC_GEN_PRODUCER_COMPLETE(results);
}

template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::regularize(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound) const {

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REG_OBJ;

	std::vector<RC_REG_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regResults, (&RegularChain<Field, RecursivePoly>::_regularize), this, p, lazardDecompose, inputHeightBound);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegResult);
	RC_GEN_CONSUMER_LOOP(regResults, currRegResult, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regResults, currRegResult, i);
		results.push_back(currRegResult);
	}

	return results;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isRegularizationTrivial(const RecursivePoly& p, RecursivePoly& pReduced) const {

	if (this->isEmpty() || p.isZero() || this->isConstantPolynomial(p)) {
		pReduced = p;
		return true;
	}
	if (this->isInSaturatedIdealMinimal(p)) {
		pReduced.zero();
		return true;
	}
	return false;
}


template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::regularizeTrivial(const RecursivePoly& p, const RecursivePoly& pReduced) const {
	unsigned long long int regularizeStart = 0;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] entering regularizeTrivial: " << p << std::endl;
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] T = " << *this << std::endl;
	#endif

	if (!(this->isConstantPolynomial(p) || this->isEmpty() || pReduced.isZero())) {
		std::cerr << "BPAS: error, input to regularizeTrivial must have a constant polynomial, an empty chain or p must be in Sat(T)." << std::endl;
		exit(1);
	}

	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results;

	if (p.isZero() || pReduced.isZero()) {
		RecursivePoly zero;
		zero.zero();
		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularizeTrivial (zero or in Sat(T)): " << p << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
	}
	else {
		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,*this));
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (const or empty): " << p << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
	}
	return results;
}

// STREAMLINED VERSION //
#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_regularize(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::_regularize(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound) const {
#endif
	++regularizeDepth;
	++depth;
	unsigned long long int regularizeStart = 0;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] entering regularize: " << p << std::endl;
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] inputHeightBound = " << inputHeightBound << std::endl;
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] T = " << *this << std::endl;
	#endif

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REG_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REG_OBJ, results);
	RC_REG_OBJ result;
	typedef RC_REG_OBJ RC_REGLISTSTREAM_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGLISTSTREAM_OBJ, currRegListStreamTask);
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> tasks,moreTasks,evenMoreTasks;
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results,results2,tasks,moreTasks,evenMoreTasks;
	RegularChain<Field,RecursivePoly> T;
	RecursivePoly zero;
	zero.zero();

//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] checking zero:" << std::endl;
//	#endif
	if (p.isZero()) {
		result = RC_REG_OBJ(zero,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (zero): " << p << std::endl;
			for (size_t i=0; i<results.size(); ++i)
				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] checking const or empty:" << std::endl;
//	#endif
	if (this->isConstantPolynomial(p) || this->isEmpty()) {
		result = RC_REG_OBJ(p,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,*this));
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (const or empty): " << p << std::endl;
//			for (size_t i=0; i<results.size(); ++i)
//				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}

	// std::cerr << "[regularize] p : " << p << std::endl;
	// std::cerr << "[regularize] prime: " << this->isSaturatedIdealPrime() << std::endl;
	// std::cerr << "[regularize] rc: " << *this << std::endl;

	if (this->isSaturatedIdealPrime()) {
		if (this->isInSaturatedIdealMinimal(p)) {
			result = RC_REG_OBJ(zero, *this);
		} else {
			result = RC_REG_OBJ(p, *this);
		}
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
		--regularizeDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
	}

	RecursivePoly redP = removeZero(p);
	if (redP.isConstant() != 0) {
		result = RC_REG_OBJ(redP,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (remove zero was constant(T)): " << p << std::endl;
			for (size_t i=0; i<results.size(); ++i)
				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
	}

	// if (isRegular(redP)) {
	// 	result = RC_REG_OBJ(redP,*this);
	// 	RC_GEN_PRODUCER_ACCUMULATE(results, result);
	// 	#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
	// 		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (remove zero was constant(T)): " << p << std::endl;
	// 		for (size_t i=0; i<results.size(); ++i)
	// 			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
	// 	#endif
	// 	#ifdef REGULARCHAIN_PROFILING
	// 		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	// 	#endif
	// 	--regularizeDepth;
	// 	--depth;
	// 	RC_GEN_PRODUCER_COMPLETE(results);
	// }

//Jan27/2020: by calling remove zero, avoid calling isInSaturatedIdeal.
// 	// TODO: this test should be avoided when the polynomial is already reduced
// 	#ifdef REGULARCHAIN_DEBUG
// 		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] checking in Sat(T):" << std::endl;
// 	#endif
// 	if (isInSaturatedIdeal(p)) {
// 		result = RC_REG_OBJ(zero,*this);
// 		RC_GEN_PRODUCER_ACCUMULATE(results, result);
// //		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
// 		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
// 			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (inSat(T)): " << p << std::endl;
// 			for (size_t i=0; i<results.size(); ++i)
// 				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
// 		#endif
// 		#ifdef REGULARCHAIN_PROFILING
// 			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
// 		#endif
// 		--regularizeDepth;
// 		--depth;
// 		RC_GEN_PRODUCER_COMPLETE(results);
// //		return results;
// 	}

	int heightBound;
	std::vector<RecursivePoly> knownRegular,unknownIfRegular;
	std::vector<RegularChain<Field,RecursivePoly>> singularComponents,regularComponents;

	if (lazardDecompose)
		heightBound = this->numberOfVariables();
	else
		heightBound = inputHeightBound;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
	unknownIfRegular = GCDFreeFactorization(redP);
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	if (unknownIfRegular[0].isConstant() != 0)
		unknownIfRegular.erase(unknownIfRegular.begin());

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] calling regularizeList:" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
	RC_GEN_CONSUMER_INIT(RC_REGLISTSTREAM_OBJ, regListStreamTasks, (&RegularChain<Field, RecursivePoly>::regularizeListStream), this, redP, knownRegular, unknownIfRegular, lazardDecompose, heightBound);
//	results2 = regularizeListStream(p, knownRegular,unknownIfRegular,lazardDecompose,heightBound);
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] returning to regularize from regularizeList:" << std::endl;
	#endif
	// TODO: Define a move constructor for PolyChainPair and use std::move on the elements of the vector
//	for (size_t i=0; i<singularComponents.size(); ++i) {
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,singularComponents[i]));
//	}
//	for (size_t i=0; i<regularComponents.size(); ++i) {
////		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,regularComponents[i]));
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(regularComponents[i].reduce(p),regularComponents[i]));
//	}
	RC_GEN_CONSUMER_LOOP(regListStreamTasks, currRegListStreamTask, i) {
//	for (auto rc : results2) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regListStreamTasks, currRegListStreamTask, i);
		RC_GEN_PRODUCER_ACCUMULATE(results, currRegListStreamTask);
//		results.emplace_back(rc);
	}

	#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize: " << p << std::endl;
		for (size_t i=0; i<results.size(); ++i)
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
	--regularizeDepth;
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

// ORIGINAL VERSION //
//template <class Field, class RecursivePoly>
//std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::regularize(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound) const {
//	++regularizeDepth;
//	++depth;
//	unsigned long long int regularizeStart = 0;
//	#ifdef REGULARCHAIN_PROFILING
//		startTimer(&regularizeStart);
//	#endif
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] entering regularize: " << p << std::endl;
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] inputHeightBound = " << inputHeightBound << std::endl;
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] T = " << *this << std::endl;
//	#endif
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results,tasks,moreTasks,evenMoreTasks;
//	RegularChain<Field,RecursivePoly> T;
//	RecursivePoly zero;
//	zero.zero();
//
////	#ifdef REGULARCHAIN_DEBUG
////		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] checking zero:" << std::endl;
////	#endif
//	if (p.isZero()) {
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (zero): " << p << std::endl;
//			for (size_t i=0; i<results.size(); ++i)
//				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
//		#endif
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
//		#endif
//		--regularizeDepth;
//		--depth;
//		return results;
//	}
////	#ifdef REGULARCHAIN_DEBUG
////		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] checking const or empty:" << std::endl;
////	#endif
//	if (this->isConstantPolynomial(p) || this->isEmpty()) {
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,*this));
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (const or empty): " << p << std::endl;
//			for (size_t i=0; i<results.size(); ++i)
//				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
//		#endif
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
//		#endif
//		--regularizeDepth;
//		--depth;
//		return results;
//	}
//	// TODO: this test should be avoided when the polynomial is already reduced
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] checking in Sat(T):" << std::endl;
//	#endif
//	if (isInSaturatedIdeal(p)) {
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (inSat(T)): " << p << std::endl;
//			for (size_t i=0; i<results.size(); ++i)
//				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
//		#endif
//		#ifdef REGULARCHAIN_PROFILING
//			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
//		#endif
//		--regularizeDepth;
//		--depth;
//		return results;
//	}
//
//	int heightBound;
//	std::vector<RecursivePoly> knownRegular,unknownIfRegular;
//	std::vector<RegularChain<Field,RecursivePoly>> singularComponents,regularComponents;
//
//	if (lazardDecompose)
//		heightBound = this->numberOfVariables();
//	else
//		heightBound = inputHeightBound;
//
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
//	#endif

//	#ifdef REGULARCHAIN_PROFILING
//		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
//	#endif
//	unknownIfRegular = GCDFreeFactorization(p);
//	#ifdef REGULARCHAIN_PROFILING
//		startTimer(&regularizeStart);
//	#endif
//
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] calling regularizeList:" << std::endl;
//	#endif
//	#ifdef REGULARCHAIN_PROFILING
//		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
//	#endif
//	regularizeList(knownRegular,unknownIfRegular,singularComponents,regularComponents,lazardDecompose,heightBound);
//	#ifdef REGULARCHAIN_PROFILING
//		startTimer(&regularizeStart);
//	#endif
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] returning to regularize from regularizeList:" << std::endl;
//	#endif
//	// TODO: Define a move constructor for PolyChainPair and use std::move on the elements of the vector
//	for (size_t i=0; i<singularComponents.size(); ++i) {
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,singularComponents[i]));
//	}
//	for (size_t i=0; i<regularComponents.size(); ++i) {
////		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,regularComponents[i]));
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(regularComponents[i].reduce(p),regularComponents[i]));
//	}
//
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize: " << p << std::endl;
//		for (size_t i=0; i<results.size(); ++i)
//			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
//	#endif
//	#ifdef REGULARCHAIN_PROFILING
//		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
//	#endif
//	--regularizeDepth;
//	--depth;
//	return results;
//}


#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::regularizeSingle(const RecursivePoly& p_in, bool lazardDecompose, int heightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::regularizeSingle(const RecursivePoly& p_in, bool lazardDecompose, int heightBound) const {
#endif
	long long unsigned int rcProfilingStart,regularizeStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularizeStart);
	#endif
	++regularizeSingleDepth;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] entering regularizeSingle: " << p_in << std::endl;
		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif
	bool reduceOutput(false);
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);
	std::vector<RegularChain<Field,RecursivePoly>> constructedChains;
	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGGCD_OBJ, currRegGCDTask);
	typedef RC_REGGCD_OBJ RC_REGSINGLE_OBJ;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REGSINGLE_OBJ, results);
	RC_REGSINGLE_OBJ result;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGSINGLE_OBJ, currRegSingleTask);
	typedef RegularChain<Field,RecursivePoly> RC_INT_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_INT_OBJ, currIntTask);
	typedef RC_REGGCD_OBJ RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask2);

	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> tasks,moreTasks,evenMoreTasks,newRegResults,newRegResults2;
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results,tasks,moreTasks,evenMoreTasks;
	std::vector<RegularChain<Field,RecursivePoly>> factorizationTasks,newIntResults;
	RegularChain<Field,RecursivePoly> T;
	RecursivePoly zero,pRed;
	zero.zero();

//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] checking zero:" << std::endl;
//	#endif
	if (p_in.isZero()) {
		result = RC_REGSINGLE_OBJ(zero,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] leaving regularizeSingle (zero): " << p_in << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeSingleDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}
//	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] checking const or empty:" << std::endl;
//	#endif
	if (this->isConstantPolynomial(p_in) || this->isEmpty()) {
		result = RC_REGSINGLE_OBJ(p_in,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
//		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,*this));
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] leaving regularizeSingle (const or empty): " << p_in << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeSingleDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}

	// std::cerr << "[regSingle] p : " << p_in << std::endl;
	// std::cerr << "[regSingle] prime: " << this->isSaturatedIdealPrime() << std::endl;
	// std::cerr << "[regSingle] rc: " << *this << std::endl;

	if (this->isSaturatedIdealPrime()) {
		if (this->isInSaturatedIdealMinimal(p_in)) {
			result = RC_REG_OBJ(zero, *this);
		} else {
			result = RC_REG_OBJ(p_in, *this);
		}
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
		--regularizeSingleDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
	}

	RecursivePoly p = removeZero(p_in);
	if (p.isConstant() != 0) {
		result = RC_REG_OBJ(p,*this);
		RC_GEN_PRODUCER_ACCUMULATE(results, result);
		#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
			std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] leaving regularize (remove zero was constant(T)): " << p << std::endl;
			for (size_t i=0; i<results.size(); ++i)
				std::cerr << "[reg][" << regularizeDepth << "][" << depth << "] results[" << i << "] = (" << results[i].poly << "," << results[i].chain << ")" << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
	}


// Jan27/2020: Call removeZero in favour of this.
// 	// TODO: this test should be avoided when the polynomial is already reduced
// 	#ifdef REGULARCHAIN_DEBUG
// 		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] checking in Sat(T):" << std::endl;
// 	#endif
// 	if (isInSaturatedIdeal(p)) {
// 		result = RC_REGSINGLE_OBJ(zero,*this);
// 		RC_GEN_PRODUCER_ACCUMULATE(results, result);
// //		results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,*this));
// 		#ifdef REGULARCHAIN_DEBUG
// 			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] leaving regularizeSingle (inSat(T)): " << p << std::endl;
// 		#endif
// 		#ifdef REGULARCHAIN_PROFILING
// 			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
// 		#endif
// 		--regularizeSingleDepth;
// 		--depth;
// 		RC_GEN_PRODUCER_COMPLETE(results);
// //		return results;
// 	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] checking can compute in dimension zero:" << std::endl;
	#endif
	if(canComputeInDimensionZero(p)) {
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] can compute in dimension zero!!" << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] T = " << *this << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] T variables:" << std::endl;
			printVariables(this->variables());
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] p = " << p << std::endl;
		#endif
		ZeroDimensionalRegularChain<Field,RecursivePoly> rc(*this,ASSUME_MAKESCHAIN);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] calling regularize in ZDRC..." << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] *this = " << rc << std::endl;
		#endif

		typedef PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>> RC_ZDREG_OBJ;
		RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_ZDREG_OBJ, currZDRegTask);
//		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> zdresults;

		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		RC_GEN_CONSUMER_INIT(RC_ZDREG_OBJ, zdresults, (&ZeroDimensionalRegularChain<Field, RecursivePoly>::_regularize), &(rc), p);
//		zdresults = rc.regularize(p);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&regularizeStart);
		#endif

		#ifdef REGULARCHAIN_DEBUG
//		if (zdresults.empty())
//			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] no results from regularize in ZDRC" << std::endl;
		#endif
		RC_GEN_CONSUMER_LOOP(zdresults, currZDRegTask, i) {
//		for (size_t i=0; i<zdresults.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(zdresults, currZDRegTask, i);
			#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] zdresults[" << i << "] = {" << zdresults[i].poly << "," << zdresults[i].chain << "}" << std::endl;
			#endif
			RegularChain<Field,RecursivePoly> T(std::move(currZDRegTask.chain));
			result = RC_REGSINGLE_OBJ(std::move(currZDRegTask.poly),T);
//			RegularChain<Field,RecursivePoly> T(std::move(zdresults[i].chain));
//			result = RC_REGSINGLE_OBJ(std::move(zdresults[i].poly),T);
			RC_GEN_PRODUCER_ACCUMULATE(results, result);
//			results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(std::move(zdresults[i].poly),T));
		}
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] leaving regularizeSingle (inDim0): " << p << std::endl;
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		--regularizeSingleDepth;
		--depth;
		RC_GEN_PRODUCER_COMPLETE(results);
//		return results;
	}

	int newHeightBound;
	RegularChain<Field,RecursivePoly> E;
	RecursivePoly q;
	Symbol v(p.leadingVariable());
	std::vector<RecursivePoly> lf;
	#ifdef REGULARCHAIN_DEBUG
//		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] this->mainVariables:" << std::endl;
//		printVariables(this->mainVariables());
		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] p = " << p << std::endl;
		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] v = " << v << std::endl;
	#endif
	if (!isAMemberOf(v,this->mainVariables())) {

		if (this->isRegularizationTrivial(p.initial(),pRed)) {
//			std::cerr << "First trivial regularization in regularizeSingle!" << std::endl;
//			std::cerr << "p.initial = " << p.initial() << std::endl;
//			std::cerr << "pRed = " << pRed << std::endl;
//			std::cerr << "this = " << *this << std::endl;
			newRegResults = this->regularizeTrivial(p.initial(),pRed);

			for (size_t i=0; i<newRegResults.size(); ++i) {
				const RecursivePoly& f = newRegResults[i].poly;
				const RegularChain<Field,RecursivePoly>& C = newRegResults[i].chain;

				if (f.isZero()) {

					if (C.isRegularizationTrivial(p.tail(),pRed)) {
//						std::cerr << "Second trivial regularization in regularizeSingle!" << std::endl;
						newRegResults2 = C.regularizeTrivial(p.tail(),pRed);
						for (size_t j=0; j<newRegResults2.size(); ++j) {
							RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults2[j]);
						}
					}
					else {
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&regularizeStart,&regularizeTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks2, (&RegularChain<Field, RecursivePoly>::_regularize), &(C), p.tail(), lazardDecompose, heightBound);
		//				moreTasks = C.regularize(p.tail(),lazardDecompose,heightBound);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&regularizeStart);
						#endif

						RC_GEN_CONSUMER_LOOP(regTasks2, currRegTask2, j) {
							RC_GEN_CONSUMER_GET_LOOPELEM(regTasks2, currRegTask2, j);
							RC_GEN_PRODUCER_ACCUMULATE(results, currRegTask2);
						}
					}
				}
				else {
					if (reduceOutput) {
						result = RC_REGSINGLE_OBJ(C.reduceMinimal(p),C);
//						result = RC_REGSINGLE_OBJ(C.reduce(p),C);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(C.reduce(p),C));
					}
					else {
						result = RC_REGSINGLE_OBJ(p,C);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,C));
					}
				}
			}
		}
		else {
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&regularizeStart);
			#endif
			RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), this, p.initial(), lazardDecompose, heightBound);
	//		tasks = regularize(p.initial(),lazardDecompose,heightBound);
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&regularizeStart);
			#endif

			RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, i) {
	//		for (size_t i=0; i<tasks.size(); ++i) {
				RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, i);
				const RecursivePoly& f = currRegTask.poly;
				const RegularChain<Field,RecursivePoly>& C = currRegTask.chain;
	//			const RecursivePoly& f = tasks[i].poly;
	//			const RegularChain<Field,RecursivePoly>& C = tasks[i].chain;

				if (f.isZero()) {

					if (C.isRegularizationTrivial(p.tail(),pRed)) {
//						std::cerr << "Second trivial regularization in regularizeSingle!" << std::endl;
						newRegResults2 = C.regularizeTrivial(p.tail(),pRed);
						for (size_t j=0; j<newRegResults2.size(); ++j) {
							RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults2[j]);
						}
					}
					else {
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&regularizeStart,&regularizeTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks2, (&RegularChain<Field, RecursivePoly>::_regularize), &(C), p.tail(), lazardDecompose, heightBound);
		//				moreTasks = C.regularize(p.tail(),lazardDecompose,heightBound);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&regularizeStart);
						#endif

						RC_GEN_CONSUMER_LOOP(regTasks2, currRegTask2, j) {
							RC_GEN_CONSUMER_GET_LOOPELEM(regTasks2, currRegTask2, j);
							RC_GEN_PRODUCER_ACCUMULATE(results, currRegTask2);
						}
					}
				}
				else {
					if (reduceOutput) {
						result = RC_REGSINGLE_OBJ(C.reduceMinimal(p),C);
//						result = RC_REGSINGLE_OBJ(C.reduce(p),C);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(C.reduce(p),C));
					}
					else {
						result = RC_REGSINGLE_OBJ(p,C);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,C));
					}
				}
			}
		}
	}
	else {
		std::vector<RegularChain<Field,RecursivePoly>> extendedChains;
		RegularChain<Field,RecursivePoly> Tlv,Tgv;
		RecursivePoly r,Tv;
		cutChain(*this,v,Tlv,Tv,Tgv);
		newHeightBound = heightBound - Tgv.numberOfAlgebraicVariables() - 1;

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] newHeightBound = " << newHeightBound << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tlv = " << Tlv << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tgv = " << Tgv << std::endl;

			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] constructing SubResultantChain in regularizeSingle:" << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] p = " << p << std::endl;
			printVariables(p.ringVariables(),"p");
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
			printVariables(Tv.ringVariables(),"Tv");
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] v = " << v << std::endl;
			if (p.leadingVariable() != Tv.leadingVariable()) {
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] error before src computation." << std::endl;
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] p = " << p << std::endl;
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] v = " << v << std::endl;
				printVariables(this->mainVariables(),"this");
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] this = " << *this << std::endl;
				exit(1);
			}
		#endif
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "setting q = p:" << std::endl;
//		#endif
		q = p;
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "reducing ringVariables of q" << std::endl;
//		#endif
		q.setRingVariables(q.variables());
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "constructing src" << std::endl;
//		#endif
		SubResultantChain<RecursivePoly,RecursivePoly> src(q,Tv,v,this->allVariables());
//		#ifdef REGULARCHAIN_DEBUG
//			std::cerr << "getting resultant" << std::endl;
//		#endif
		r = src.resultant();
		r.setRingVariables(this->allVariables());
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&subresultantChainTime);
		#endif

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] quantities after src computation." << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] p = " << p << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] v = " << v << std::endl;
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] src = " << std::endl;
			for (size_t ll=0; ll<src.size(); ++ll) {
				std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] S[" << ll << "] = " << src.subResultantOfIndex(ll) << std::endl;
			}
//			printVariables(this->mainVariables(),"this");
			std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] this = " << *this << std::endl;
		#endif

		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&regularizeStart,&regularizeTime);
		#endif
		RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(Tlv), r, lazardDecompose, newHeightBound);
//		tasks = Tlv.regularize(r,lazardDecompose,newHeightBound);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&regularizeStart);
		#endif

		RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, i) {
//		for (size_t i=0; i<tasks.size(); ++i) {
			RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, i);
			const RecursivePoly& f = currRegTask.poly;
			const RegularChain<Field,RecursivePoly>& C = currRegTask.chain;
//			const RecursivePoly& f = tasks[i].poly;
//			const RegularChain<Field,RecursivePoly>& C = tasks[i].chain;`

			if (C.dimension() < Tlv.dimension()) {

				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&regularizeStart,&regularizeTime);
				#endif
				extendedChains = C.extend(*this,v,lazardDecompose,heightBound);
//				extendedChains = C.extend(*this,v);
				#ifdef REGULARCHAIN_PROFILING
					startTimer(&regularizeStart);
				#endif

				for (size_t j=0; j<extendedChains.size(); ++j) {

					if (extendedChains[j].isRegularizationTrivial(p,pRed)) {
//						std::cerr << "Third trivial regularization in regularizeSingle!" << std::endl;
						newRegResults = extendedChains[j].regularizeTrivial(p,pRed);
						for (size_t k=0; k<newRegResults.size(); ++k) {
							RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults[k]);
						}
					}
					else {
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&regularizeStart,&regularizeTime);
						#endif
						RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), &(extendedChains[j]), p, lazardDecompose, heightBound);
	//					moreTasks = extendedChains[j].regularizeSingle(p,lazardDecompose,heightBound);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&regularizeStart);
						#endif

						RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, k) {
							RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, k);
							RC_GEN_PRODUCER_ACCUMULATE(results, currRegSingleTask);
						}
					}
				}
			}
			else if (!f.isZero()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] constructing C U T>=v " << std::endl;
				#endif

				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&regularizeStart,&regularizeTime);
				#endif
//				currRegTask.chain.constructChain(Tv,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);
//				currRegTask.chain.constructChain(Tgv,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);
				// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct

				int constructOptions = ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE;
				if (this->isSquareFree()) {
					constructOptions |= ASSUME_SQUAREFREE;
				}
				RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &currRegTask.chain, Tv, lazardDecompose, newHeightBound, constructOptions);

				RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, k1) {
					RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, k1);

					// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
					constructedChains = currConstructChainsTask.constructChainsFromChain(Tgv,lazardDecompose,newHeightBound, constructOptions);

					for (size_t k2=0; k2<constructedChains.size(); k2++) {
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&regularizeStart);
						#endif

						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] chain = " << currRegTask.chain << " added to results." << std::endl;
							std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] after chain construction (!f.isZero())." << std::endl;
							std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
						#endif
						if (reduceOutput) {
							result = RC_REGSINGLE_OBJ(constructedChains[k2].reduceMinimal(p),constructedChains[k2]);
//							result = RC_REGSINGLE_OBJ(constructedChains[k2].reduce(p),constructedChains[k2]);
							RC_GEN_PRODUCER_ACCUMULATE(results, result);
		//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(tasks[i].chain.reduce(p),tasks[i].chain));
						}
						else {
							result = RC_REGSINGLE_OBJ(p,constructedChains[k2]);
		//					result = RC_REGSINGLE_OBJ(p,tasks[i].chain);
							RC_GEN_PRODUCER_ACCUMULATE(results, result);
		//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,tasks[i].chain));
						}
					}
				}
			}
			else {
				#ifdef REGULARCHAIN_DEBUG
					if (p.leadingVariable() != Tv.leadingVariable()) {
						std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] error before regularGCD computation." << std::endl;
						std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] p = " << p << std::endl;
						std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
						std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] v = " << v << std::endl;
//						printVariables(this->mainVariables(),"this");
						std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] this = " << *this << std::endl;
						exit(1);
					}
				#endif

				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&regularizeStart,&regularizeTime);
				#endif
				RC_GEN_CONSUMER_INIT(RC_REGGCD_OBJ, regGCDTasks, (&RegularChain<Field, RecursivePoly>::_regularGCD), &(C), p, Tv, v, src, lazardDecompose, heightBound);
//				moreTasks = C.regularGCD(p,Tv,v,src,lazardDecompose,heightBound);
				#ifdef REGULARCHAIN_PROFILING
					startTimer(&regularizeStart);
				#endif

				RC_GEN_CONSUMER_LOOP(regGCDTasks, currRegGCDTask, j) {
//				for (size_t j=0; j<moreTasks.size(); ++j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(regGCDTasks, currRegGCDTask, j);
					const RecursivePoly& g = currRegGCDTask.poly.mainPrimitivePart();
					const RegularChain<Field,RecursivePoly>& D = currRegGCDTask.chain;

					if (D.dimension() < C.dimension()) {

						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&regularizeStart,&regularizeTime);
						#endif
						extendedChains = D.extend(*this,v,lazardDecompose,heightBound);
//						extendedChains = D.extend(*this,v);
						#ifdef REGULARCHAIN_PROFILING
							startTimer(&regularizeStart);
						#endif

						for (size_t k=0; k<extendedChains.size(); ++k) {

							if (extendedChains[k].isRegularizationTrivial(p,pRed)) {
//								std::cerr << "Fourth trivial regularization in regularizeSingle!" << std::endl;
								newRegResults = extendedChains[k].regularizeTrivial(p,pRed);
								for (size_t l=0; l<newRegResults.size(); ++l) {
									RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults[l]);
								}
							}
							else {
								#ifdef REGULARCHAIN_PROFILING
									stopTimerAddElapsed(&regularizeStart,&regularizeTime);
								#endif
								RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), &(extendedChains[k]), p, lazardDecompose, heightBound);
	//							evenMoreTasks = extendedChains[k].regularizeSingle(p,lazardDecompose,heightBound);
								#ifdef REGULARCHAIN_PROFILING
									startTimer(&regularizeStart);
								#endif

								RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, l) {
									RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, l);
									RC_GEN_PRODUCER_ACCUMULATE(results, currRegSingleTask);
								}
							}
						}
					}
					else {
						if (g.mainDegree() == Tv.mainDegree()) {
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] constructing D U T>=v " << std::endl;
							#endif

							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&regularizeStart,&regularizeTime);
							#endif
//							currRegGCDTask.chain.constructChain(Tv,ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);
//							currRegGCDTask.chain.constructChain(Tgv,ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);

							// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
							int constructOptions = ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE;
							if (this->isSquareFree()) {
								//that is, Tv or Tgv came from a square free RC and dimension has not dropped,
								constructOptions |= ASSUME_SQUAREFREE;
							}
							RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &currRegGCDTask.chain, Tv, lazardDecompose, newHeightBound, constructOptions);

							RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, k1) {
								RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, k1);

								// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
								constructedChains = currConstructChainsTask.constructChainsFromChain(Tgv,lazardDecompose,newHeightBound, constructOptions);

								for (size_t k2=0; k2<constructedChains.size(); k2++) {
									#ifdef REGULARCHAIN_PROFILING
										startTimer(&regularizeStart);
									#endif

									#ifdef REGULARCHAIN_DEBUG
										std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] chain = " << constructedChains[k2] << " added to results." << std::endl;
										std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] after chain construction (g.mainDegree() == Tv.mainDegree())." << std::endl;
										std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;

										std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] potential source of problems, since I'm changing the value of a constant reference here:" << std::endl;
									#endif
									result = RC_REGSINGLE_OBJ(zero,constructedChains[k2]);
									RC_GEN_PRODUCER_ACCUMULATE(results, result);
		//							results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,D));
								}
							}
						}
						else {
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] constructing D U g U T>v " << std::endl;
							#endif

							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&regularizeStart,&regularizeTime);
							#endif
							#if defined(WITH_FACTORIZATION)
								lf = GCDFreeFactorization(g,2);
							#else
								lf = GCDFreeFactorization(g,0);
							#endif
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&regularizeStart);
							#endif

							int constructOptions = ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE;
							if (this->isSquareFree()) {
								//that is, Tv or Tgv came from a square free RC and dimension has not dropped,
								//moreover, since Tv was square free, g must be also.
								constructOptions |= ASSUME_SQUAREFREE;
							}
							for (size_t m=0; m<lf.size(); ++m) {
								if (lf[m].isConstant() == 0) {

									#ifdef REGULARCHAIN_PROFILING
										stopTimerAddElapsed(&regularizeStart,&regularizeTime);
									#endif
									E = D;
//									E.constructChain(lf[m],ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);
//									E.constructChain(Tgv,ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);

									// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
									RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &E, lf[m], lazardDecompose, newHeightBound, constructOptions);

									RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, k1) {
										RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, k1);

										// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
										constructedChains = currConstructChainsTask.constructChainsFromChain(Tgv,lazardDecompose,newHeightBound,constructOptions);

										for (size_t k2=0; k2<constructedChains.size(); k2++) {
											#ifdef REGULARCHAIN_PROFILING
												startTimer(&regularizeStart);
											#endif

											#ifdef REGULARCHAIN_DEBUG
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] E = " << constructedChains[k2] << " added to results." << std::endl;
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] g factor = " << lf[m] << std::endl;
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
											#endif
											result = RC_REGSINGLE_OBJ(zero,constructedChains[k2]);
											RC_GEN_PRODUCER_ACCUMULATE(results, result);
		//									results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(zero,E));
										}
									}
								}
							}
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&rcProfilingStart);
							#endif

							RecursivePoly prem = Tv.pseudoDivide(g,&q,NULL,1);
							q = q.mainPrimitivePart();
							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
							#endif
							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] r = " << prem << std::endl;
								std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] q = " << q << std::endl;
							#endif

							#ifdef REGULARCHAIN_DEBUG
								std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] constructing D U q U T>v " << std::endl;
							#endif

							#ifdef REGULARCHAIN_PROFILING
								stopTimerAddElapsed(&regularizeStart,&regularizeTime);
							#endif
							#if defined(WITH_FACTORIZATION)
								lf = GCDFreeFactorization(q,2);
							#else
								lf = GCDFreeFactorization(q,0);
							#endif
							#ifdef REGULARCHAIN_PROFILING
								startTimer(&regularizeStart);
							#endif

//							lf = GCDFreeFactorization(q);
//							char input;
//							std::cout << "Press enter to continue: ";
//							input = getchar();
//							std::cout << std::endl;
//							std::cerr << "q.leadingVariable() = " << q.leadingVariable() << std::endl;
							for (size_t m=0; m<lf.size(); ++m) {
								if (lf[m].isConstant() == 0 && lf[m].degree(q.leadingVariable()) > 0) {

									#ifdef REGULARCHAIN_PROFILING
										stopTimerAddElapsed(&regularizeStart,&regularizeTime);
									#endif
									E = D;
									// TODO: ONLY PASS ASSUME_SQUAREFREE as an option if
									// ((regularChainOptions & MAINTAIN_SQUAREFREE) == MAINTAIN_SQUAREFREE)?
//									E.constructChain(lf[m],ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);
//									E.constructChain(Tgv,ASSUME_MAKESCHAIN | ASSUME_REGULAR  | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);

									// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
									RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &E, lf[m], lazardDecompose, newHeightBound, constructOptions);

									RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, k1) {
										RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, k1);

										// TODO: Check whether the height bound is correct here; re-verify the assumptions are correct
										constructedChains = currConstructChainsTask.constructChainsFromChain(Tgv,lazardDecompose,newHeightBound, constructOptions);

										for (size_t k2=0; k2<constructedChains.size(); k2++) {
											#ifdef REGULARCHAIN_PROFILING
												startTimer(&regularizeStart);
											#endif

											#ifdef REGULARCHAIN_DEBUG
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] E = " << constructedChains[k2] << " added to tasks for regularize." << std::endl;
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] q factor = " << lf[m] << std::endl;
			//									std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] Tv = " << Tv << std::endl;
											#endif
											factorizationTasks.emplace_back(std::move(constructedChains[k2]));
										}
									}
								}
							}

							for (size_t m=0; m<factorizationTasks.size(); ++m) {
//							RegularChain<Field,RecursivePoly> F(E);
								// TODO: uncomment when isSquareFree is maintained properly
								const RegularChain<Field,RecursivePoly>& F = factorizationTasks[m];
								if (!F.isSquareFree()) {
//								if ((regularChainOptions & MAINTAIN_SQUAREFREE) != MAINTAIN_SQUAREFREE) {

									if (F.isRegularizationTrivial(p,pRed)) {
//										std::cerr << "Fifth trivial regularization in regularizeSingle!" << std::endl;
										newRegResults = F.regularizeTrivial(p,pRed);
										for (size_t n=0; n<newRegResults.size(); ++n) {
											#ifdef REGULARCHAIN_DEBUG
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] " << newRegResults[n] << " added to results." << std::endl;
											#endif
											RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults[n]);
										}
									}
									else {
										#ifdef REGULARCHAIN_PROFILING
											stopTimerAddElapsed(&regularizeStart,&regularizeTime);
										#endif
										RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), &(F), p, lazardDecompose, heightBound);
	//									evenMoreTasks = F.regularizeSingle(p,lazardDecompose,heightBound);
										#ifdef REGULARCHAIN_PROFILING
											startTimer(&regularizeStart);
										#endif

										RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, n) {
											RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, n);
											#ifdef REGULARCHAIN_DEBUG
												std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] " << currRegSingleTask << " added to results." << std::endl;
											#endif
											RC_GEN_PRODUCER_ACCUMULATE(results, currRegSingleTask);
										}
									}
								}
								else {
									if (reduceOutput) {
										result = RC_REGSINGLE_OBJ(F.reduceMinimal(p),F);
										#ifdef REGULARCHAIN_DEBUG
											std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] F = " << result << " added to results." << std::endl;
										#endif
										RC_GEN_PRODUCER_ACCUMULATE(results, result);
//										results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(F.reduce(p),F));
									}
									else {
										result = RC_REGSINGLE_OBJ(p,F);
										#ifdef REGULARCHAIN_DEBUG
											std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] F = " << result << " added to results." << std::endl;
										#endif
										RC_GEN_PRODUCER_ACCUMULATE(results, result);
//										results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(p,F));
									}
								}
							}

							//Jan23/2020: Should be a logical OR here.
							// if (lazardDecompose && D.numberOfAlgebraicVariables() < newHeightBound) {
							if (lazardDecompose || D.numberOfAlgebraicVariables() < newHeightBound) {
								std::vector<RegularChain<Field,RecursivePoly>> intersectedChains;

								#ifdef REGULARCHAIN_PROFILING
									stopTimerAddElapsed(&regularizeStart,&regularizeTime);
								#endif
								#if defined(WITH_FACTORIZATION)
									lf = GCDFreeFactorization(g.initial(),2);
								#else
									lf = GCDFreeFactorization(g.initial(),0);
								#endif
	//							lf = GCDFreeFactorization(g.initial());
								#ifdef REGULARCHAIN_PROFILING
									startTimer(&regularizeStart);
								#endif

								for (size_t m=0; m<lf.size(); ++m) {
									if (lf[m].isConstant() == 0) { // && lf[m].degree(g.leadingVariable()) > 0) {

										if (D.isIntersectionTrivial(lf[m],pRed)) {
//											std::cerr << "Trivial intersection in regularizeSingle!" << std::endl;
											newIntResults = D.intersectTrivial(pRed);

											for (size_t k=0; k<newIntResults.size(); ++k) {

												#ifdef REGULARCHAIN_PROFILING
													stopTimerAddElapsed(&regularizeStart,&regularizeTime);
												#endif
												extendedChains = newIntResults[k].extend(*this,v,lazardDecompose,heightBound);
												#ifdef REGULARCHAIN_PROFILING
													startTimer(&regularizeStart);
												#endif

												for (size_t l=0; l<extendedChains.size(); ++l) {

													if (extendedChains[l].isRegularizationTrivial(p,pRed)) {
//														std::cerr << "Sixth trivial regularization in regularizeSingle!" << std::endl;
														newRegResults = extendedChains[l].regularizeTrivial(p,pRed);
														for (size_t m=0; m<newRegResults.size(); ++m) {
															RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults[m]);
														}
													}
													else {
														#ifdef REGULARCHAIN_PROFILING
															stopTimerAddElapsed(&regularizeStart,&regularizeTime);
														#endif
														RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), &(extendedChains[l]), p, lazardDecompose, heightBound);
		//												evenMoreTasks = extendedChains[l].regularizeSingle(p,lazardDecompose,heightBound);
														#ifdef REGULARCHAIN_PROFILING
															startTimer(&regularizeStart);
														#endif

														RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, m) {
															RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, m);
															RC_GEN_PRODUCER_ACCUMULATE(results, currRegSingleTask);
														}
													}
												}
											}
										}
										else {
											#ifdef REGULARCHAIN_PROFILING
												stopTimerAddElapsed(&regularizeStart,&regularizeTime);
											#endif
											RC_GEN_CONSUMER_INIT(RC_INT_OBJ, intTasks, (&RegularChain<Field, RecursivePoly>::_intersect), &(D), lf[m], lazardDecompose, newHeightBound);
	//										intersectedChains = D.intersect(lf[m],lazardDecompose,newHeightBound);
											#ifdef REGULARCHAIN_PROFILING
												startTimer(&regularizeStart);
											#endif

											RC_GEN_CONSUMER_LOOP(intTasks, currIntTask, k) {
	//										for (size_t k=0; k<intersectedChains.size(); ++k) {
												RC_GEN_CONSUMER_GET_LOOPELEM(intTasks, currIntTask, k);

												#ifdef REGULARCHAIN_PROFILING
													stopTimerAddElapsed(&regularizeStart,&regularizeTime);
												#endif
												extendedChains = currIntTask.extend(*this,v,lazardDecompose,heightBound);
	//											extendedChains = intersectedChains[k].extend(*this,v,lazardDecompose,heightBound);
												#ifdef REGULARCHAIN_PROFILING
													startTimer(&regularizeStart);
												#endif

												for (size_t l=0; l<extendedChains.size(); ++l) {

													if (extendedChains[l].isRegularizationTrivial(p,pRed)) {
//														std::cerr << "Sixth trivial regularization in regularizeSingle!" << std::endl;
														newRegResults2 = extendedChains[l].regularizeTrivial(p,pRed);
														for (size_t m=0; m<newRegResults2.size(); ++m) {
															RC_GEN_PRODUCER_ACCUMULATE(results, newRegResults2[m]);
														}
													}
													else {
														#ifdef REGULARCHAIN_PROFILING
															stopTimerAddElapsed(&regularizeStart,&regularizeTime);
														#endif
														RC_GEN_CONSUMER_INIT(RC_REGSINGLE_OBJ, regSingleTasks, (&RegularChain<Field, RecursivePoly>::regularizeSingle), &(extendedChains[l]), p, lazardDecompose, heightBound);
		//												evenMoreTasks = extendedChains[l].regularizeSingle(p,lazardDecompose,heightBound);
														#ifdef REGULARCHAIN_PROFILING
															startTimer(&regularizeStart);
														#endif

														RC_GEN_CONSUMER_LOOP(regSingleTasks, currRegSingleTask, m) {
															RC_GEN_CONSUMER_GET_LOOPELEM(regSingleTasks, currRegSingleTask, m);
															RC_GEN_PRODUCER_ACCUMULATE(results, currRegSingleTask);
														}
													}
												}
											}
										}
									}
								}
							}
//							intersectedChains = D.intersect(g.initial());
//
//							for (size_t k=0; k<intersectedChains.size(); ++k) {
//								extendedChains = intersectedChains[k].extend(*this,v);
//
//								for (size_t l=0; l<extendedChains.size(); ++l) {
//									evenMoreTasks = extendedChains[l].regularizeSingle(p);
//
//									results.reserve(results.size()+evenMoreTasks.size());
//									results.insert(results.end(),evenMoreTasks.begin(),evenMoreTasks.end());
//								}
//							}
						}
					}
				}
			}
		}
	}
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[regSingle][" << regularizeSingleDepth << "][" << depth << "] leaving regularizeSingle (end of algorithm): " << p << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularizeStart,&regularizeTime);
	#endif
	--regularizeSingleDepth;
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

//#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
//template <class Field, class RecursivePoly>
//void RegularChain<Field,RecursivePoly>::extend(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const {
//#else
template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::extend(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound) const {
//#endif
	++depth;
	unsigned long long int extendStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&extendStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] entering extend(p,v):" << std::endl;
		std::cerr << "[ext][" << depth << "] heightBound = " << heightBound << std::endl;
		std::cerr << "[ext][" << depth << "] p = " << p << std::endl;
		std::cerr << "[ext][" << depth << "] C = " << *this << std::endl;
	#endif
	bool reduceOutput(false);
//	RegularChain<Field,RecursivePoly> E;
//	std::vector<RecursivePoly> lf;
	std::vector<RegularChain<Field,RecursivePoly>> results;
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> tasks;
	typedef RegularChain<Field,RecursivePoly> RC_CONSTRUCTCH_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_CONSTRUCTCH_OBJ, currConstructChainsTask);

	// NB: Maple requires p to be expanded
	// TODO: consider implementing check that the main variable of p is greater than or equal to v (otherwise v is not needed as an input)

	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&extendStart,&extendTime);
	#endif
	tasks = regularize(p.initial(),lazardDecompose,heightBound-1);
//	tasks = regularize(p.initial());
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&extendStart);
	#endif

	for (size_t i=0; i<tasks.size(); ++i) {
		if (!tasks[i].poly.isZero()) {
			if (reduceOutput) {
				tasks[i].poly = tasks[i].chain.reduceMinimal(p);
//				tasks[i].poly = tasks[i].chain.reduce(p);
			}
			else
				tasks[i].poly = p;
			if (!tasks[i].poly.isZero()) {
				#ifdef REGULARCHAIN_DEBUG
					std::cerr << "[ext][" << depth << "] constructing E U p:" << std::endl;
				#endif
//				lf = GCDFreeFactorization(tasks[i].poly);
//				for (size_t m=0; m<lf.size(); ++m) {
//					if (lf[m].isConstant() == 0 && lf[m].degree(p.leadingVariable()) > 0) {
//						char input;
//						input = getchar();
//						E = tasks[i].chain;
//						E.constructChain(lf[m],ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);
//						#ifdef REGULARCHAIN_DEBUG
//							std::cerr << "[ext][" << depth << "] chain = " << tasks[i].chain << " added to results." << std::endl;
//						#endif
//						results.emplace_back(E);
//					}
//				}

				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&extendStart,&extendTime);
				#endif
//				tasks[i].chain.constructChain(p,ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE | ASSUME_SQUAREFREE);

				// TODO: check if we need to decrease the height bound; re-verify assumptions
				// Jan23/2020: Since we called regularize, it can do many weird things and we should not assume squarefreeness.
				RC_GEN_CONSUMER_INIT(RC_CONSTRUCTCH_OBJ, constructChainsTasks, (&RegularChain<Field, RecursivePoly>::constructChainsFromPoly), &tasks[i].chain, p, lazardDecompose, heightBound, ASSUME_MAKESCHAIN | ASSUME_REGULAR | ASSUME_PRIMITIVE);

				RC_GEN_CONSUMER_LOOP(constructChainsTasks, currConstructChainsTask, j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(constructChainsTasks, currConstructChainsTask, j);
					#ifdef REGULARCHAIN_PROFILING
						startTimer(&extendStart);
					#endif

					#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[ext][" << depth << "] chain = " << currConstructChainsTask << " added to results." << std::endl;
					#endif
					results.emplace_back(std::move(currConstructChainsTask));
//					results.emplace_back(std::move(tasks[i].chain));
				}
			}
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] leaving extend(p,v):" << std::endl;
		for (size_t i=0; i<results.size(); ++i) {
			std::cerr << "[ext][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
		}
	#endif
	--depth;
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&extendStart,&extendTime);
	#endif
	return results;
}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::extend(const std::vector<RecursivePoly>& T, const Symbol& v, bool lazardDecompose, int heightBound) const {
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::extend(const std::vector<RecursivePoly>& T, const Symbol& v) const {
	++depth;
	unsigned long long int extendStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&extendStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] entering extend(Tgev.polys,v):" << std::endl;
		std::cerr << "[ext][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif
	std::vector<RegularChain<Field,RecursivePoly>> results,newTasks;

	if (T.empty()) {
		results.emplace_back(*this);
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[ext][" << depth << "] leaving extend(Tgev.polys,v) without calling extend(p,v)." << std::endl;
		#endif
		--depth;
		return results;
	}

	// TODO: consider implementing check that the elements of Tv are ordered increasingly by main variable and the main variables are strictly greater
	//       than any main variables appearing in *this (only if this function could be called from outside would that be needed).

	struct Task {
		std::vector<RecursivePoly> polys;
		RegularChain<Field,RecursivePoly> chain;
	};
	Task task = {T,*this};
	std::vector<Task> tasks;
	int newHeightBound;
	RecursivePoly poly;
	std::vector<RecursivePoly> polys;
	RegularChain<Field,RecursivePoly> chain;

	while (!tasks.empty()) {
		task = tasks.back();
		tasks.pop_back();
		if (task.polys.empty()) {
			results.emplace_back(task.chain);
		}
		polys = task.polys;
		chain = task.chain;
		poly = polys.back();
		polys.pop_back();
		newHeightBound = heightBound - polys.size();

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[ext][" << depth << "] newHeightBound = " << newHeightBound << std::endl;
		#endif

		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&extendStart,&extendTime);
		#endif
		newTasks = chain.extend(poly,v,lazardDecompose,newHeightBound);
//		newTasks = chain.extend(poly,v);
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&extendStart);
		#endif

		tasks.reserve(tasks.size()+newTasks.size());
		for (size_t i=0; i<newTasks.size(); ++i) {
			task = {polys,newTasks[i]};
			tasks.push_back(task);
		}
	}

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] leaving extend(Tgev.polys,v)." << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&extendStart,&extendTime);
	#endif
	--depth;
	return results;
}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::extend(const RegularChain<Field,RecursivePoly>& T, const Symbol& v, bool lazardDecompose, int heightBound) const {
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::extend(const RegularChain<Field,RecursivePoly>& T, const Symbol& v) const {
	++depth;
	unsigned long long int extendStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&extendStart);
	#endif
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] entering extend(T,v):" << std::endl;
		std::cerr << "[ext][" << depth << "] heightBound = " << heightBound << std::endl;
		std::cerr << "[ext][" << depth << "] v = " << v << std::endl;
		std::cerr << "[ext][" << depth << "] C = " << *this << std::endl;
		std::cerr << "[ext][" << depth << "] T = " << T << std::endl;
	#endif

	std::vector<RegularChain<Field,RecursivePoly>> results;
	RegularChain<Field,RecursivePoly> Tgev;
	RecursivePoly Tv;
	T.cutChain(v,Tv,Tgev);
	Tgev.TriangularSet<Field,RecursivePoly>::operator+=(Tv);

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] T>=v = " << Tgev << std::endl;
	#endif

	if (Tgev.isEmpty()) {
		results.emplace_back(*this);

		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[ext][" << depth << "] leaving extend(T,v) without calling extend(polys,v)." << std::endl;
			for (size_t i=0; i<results.size(); ++i) {
				std::cerr << "[ext][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
			}
		#endif
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&extendStart,&extendTime);
		#endif
		--depth;
		return results;
	}

	std::vector<RecursivePoly> polys(Tgev.polynomials()),polys2;
	for (size_t i=0; i<polys.size(); ++i) {
		if (!polys[i].isZero())
			polys2.emplace_back(std::move(polys[i]));
	}
	std::reverse(polys2.begin(),polys2.end());

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[ext][" << depth << "] leaving extend(T,v)," << std::endl;
		for (size_t i=0; i<results.size(); ++i) {
			std::cerr << "[ext][" << depth << "] results[" << i << "] = " << results[i] << std::endl;
		}
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&extendStart,&extendTime);
	#endif
	--depth;
	return extend(polys2,v,lazardDecompose,heightBound);
//	return extend(polys2,v);
}

template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& S, bool lazardDecompose, int inputHeightBound) const {

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;

	std::vector<RC_REGGCD_OBJ> results;
	RC_GEN_CONSUMER_INIT(RC_REGGCD_OBJ, regGCD, (&RegularChain<Field, RecursivePoly>::_regularGCD), this, p, q ,v, S, lazardDecompose, inputHeightBound);
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REGGCD_OBJ, currRegGCD);
	RC_GEN_CONSUMER_LOOP(regGCD, currRegGCD, i) {
		RC_GEN_CONSUMER_GET_LOOPELEM(regGCD, currRegGCD, i);
		results.push_back(currRegGCD);
	}

	return results;
}


#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::_regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& S, bool lazardDecompose, int inputHeightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const {
#else
template <class Field, class RecursivePoly>
std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> RegularChain<Field,RecursivePoly>::_regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& S, bool lazardDecompose, int inputHeightBound) const {
#endif
	unsigned long long int rcProfilingStart,regularGCDStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&regularGCDStart);
	#endif
	++regularGCDDepth;
	++depth;
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] entering regularGCD: " << p << ", " << q << std::endl;
		printVariables(p.ringVariables(), "pVars");
		printVariables(p.ringVariables(), "qVars");
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] inputHeightBound = " << inputHeightBound << std::endl;
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] S = " << std::endl;
		for (size_t k=0; k<S.size(); ++k)
			std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] S[" << k << "] = " << S.subResultantOfIndex(k) << std::endl;
	#endif

	typedef PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>> RC_REGGCD_OBJ;
	RC_REGGCD_OBJ result;
	RC_GEN_PRODUCER_DECLARE_RESULT(RC_REGGCD_OBJ, results);
	typedef RC_REGGCD_OBJ RC_REG_OBJ;
	RC_GEN_CONSUMER_DECLARE_LOOPELEM(RC_REG_OBJ, currRegTask);
//	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results;
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularComponents;
	RecursivePoly fRed;
	struct Task {
		RegularChain<Field,RecursivePoly> chain;
		int i;
	};
	std::vector<Task> tasks;
	Task currTask,newTask;
	RecursivePoly f,tpp,tc,minusOne;
	minusOne.negativeOne();
	int i,d,heightBound;
	if (lazardDecompose) {
		heightBound = this->numberOfVariables();
//		heightBound = this->numberOfAlgebraicVariables();
	}
	else
		heightBound = inputHeightBound;

	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] heightBound = " << heightBound << std::endl;
	#endif

	d = S.secondPolynomial().degree(v).get_si(); // degree of the second poly input to the src
	#ifdef REGULARCHAIN_DEBUG
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] S_" << S.size()-2 << " = " << S.secondPolynomial() << std::endl;
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] deg(S_" << S.size()-2 << ") = d = " << d << std::endl;
	#endif

	currTask = {*this,1};
	tasks.push_back(currTask);
	while (!tasks.empty()) {
		currTask = tasks.back();
		tasks.pop_back();
		i = currTask.i;
		#ifdef REGULARCHAIN_DEBUG
			std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] i = " << i << std::endl;
		#endif
		if (i>d || (i == q.degree(v).get_si())) {
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] case: i>d or i=deg(q)" << std::endl;
			#endif
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			f = q.mainPrimitivePart();
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
			#endif
			result = RC_REGGCD_OBJ(currTask.chain.reduceMinimal(f),currTask.chain);
//			result = RC_REGGCD_OBJ(currTask.chain.reduce(f),currTask.chain);
			RC_GEN_PRODUCER_ACCUMULATE(results, result);
//			results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(currTask.chain.reduce(f),currTask.chain));
		}
		else {
			//Feb2/2021:Reworking "principal subres coef" vs "initial of subres".
			//          We should use initial unless we know the RC is square free
			//          Moreover, SRC object itself handles i >= d case.
			if (currTask.chain.isSquareFree()) {
				f = S.principalSubResultantCoefficientOfIndex(i);
			} else {
			// std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth <<"] getting initial of index " << i << std::endl;
				f = S.subResultantInitialOfIndex(i);
			}
			// std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth <<"] for initial of index " << i << " " << f << std::endl;

			f.setRingVariables(this->allVariables());

			// this next line was "if (i<d) {" following the Maple code
//			if (i<=d) {
			// if (i<d) {
			// 	// f = S.principalSubResultantCoefficientOfIndex(i);
			// 	f = S.subResultantInitialOfIndex(i);
			// 	//Jan29/2020:We need to set vars right now because subresultants do not have their ambient space necessarily the same.
			// 	f.setRingVariables(this->allVariables());
			// 	// f = S.subResultantOfIndex(i).leadingCoefficientInVariable(v);
			// 	// f = RecursivePoly(S.subResultantOfIndex(i)).leadingCoefficientInVariable(v);
			// 	// S.principlalubResultantCoefficientOfIndex(i)
			// }
			// else {
			// 	f = p.initial(); // initial should work here
			// }
			#ifdef REGULARCHAIN_DEBUG
				std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] C = " << currTask.chain << std::endl;
				std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] f = " << f << std::endl;
				printVariables(f.ringVariables(), "fRingVars");
			#endif

			if (currTask.chain.isRegularizationTrivial(f,fRed)) {
//				std::cerr << "Trivial regularization in regularGCD!" << std::endl;
				regularComponents = currTask.chain.regularizeTrivial(f,fRed);

				for (size_t j=0; j<regularComponents.size(); ++j) {
					if (regularComponents[j].chain.dimension() < currTask.chain.dimension()) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] dimension drop" << std::endl;
						#endif
						result = RC_REGGCD_OBJ(RecursivePoly(),regularComponents[j].chain);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
					}
					else if (regularComponents[j].poly.isZero()) { //Jan27/2020
					// else if (regularComponents[j].chain.reduceMinimal(regularComponents[j].poly).isZero()) { // this reduce is not needed if regularize reduces output
//					else if (regularComponents[j].chain.reduce(regularComponents[j].poly).isZero()) { // this reduce is not needed if regularize reduces output
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] s_" << i << " is zero mod D" << std::endl;
						#endif
						newTask = {regularComponents[j].chain,i+1};
						tasks.push_back(newTask);
					}
					else if (i==d) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] case: i==d" << std::endl;
						#endif

						//Jan27/2020: Don't reduce output, but mainPrimPart good.
						result = RC_REGGCD_OBJ(p.mainPrimitivePart(),regularComponents[j].chain);
						// result = RC_REGGCD_OBJ(regularComponents[j].chain.reduceMinimal(p.mainPrimitivePart()),regularComponents[j].chain);
//						result = RC_REGGCD_OBJ(regularComponents[j].chain.reduce(p.mainPrimitivePart()),regularComponents[j].chain);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
					}
					else {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] case: i<d, i = " << i << std::endl;
						#endif
						f = S.subResultantOfIndex(i); //Jan17/2020: reduction causes coefficients to explode.
						#ifdef REGULARCHAIN_DEBUG
						std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth <<"] in reg trivial, i<d, got sub res " << f << std::endl;
						#endif

						f.setRingVariables(this->allVariables());
						// f = regularComponents[j].chain.reduceMinimal(S.subResultantOfIndex(i));
//						f = regularComponents[j].chain.reduce(S.subResultantOfIndex(i));


						#ifdef REGULARCHAIN_PROFILING
							startTimer(&rcProfilingStart);
						#endif
						// Replace with a check that the initial is not equal to the content, and if not return mPP(f)?
						tpp = f.mainPrimitivePart(tc);
						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
						#endif
						if ((f.initial().leadingCoefficient() * tc.leadingCoefficient()) < 0)
							tc *= minusOne;
						if (tc != f.initial()) {
							result = RC_REGGCD_OBJ(tpp,regularComponents[j].chain);
							RC_GEN_PRODUCER_ACCUMULATE(results, result);
						}
						else {
							result = RC_REGGCD_OBJ(f,regularComponents[j].chain);
							RC_GEN_PRODUCER_ACCUMULATE(results, result);
						}
					}
				}
			}
			else {
				#ifdef REGULARCHAIN_PROFILING
					stopTimerAddElapsed(&regularGCDStart,&regularGCDTime);
				#endif
				RC_GEN_CONSUMER_INIT(RC_REG_OBJ, regTasks, (&RegularChain<Field, RecursivePoly>::_regularize), &(currTask.chain), f, lazardDecompose, heightBound);
	//			regularComponents = currTask.chain.regularize(f,lazardDecompose,heightBound);
				#ifdef REGULARCHAIN_PROFILING
					startTimer(&regularGCDStart);
				#endif

				RC_GEN_CONSUMER_LOOP(regTasks, currRegTask, j) {
	//			for (size_t j=0; j<regularComponents.size(); ++j) {
					RC_GEN_CONSUMER_GET_LOOPELEM(regTasks, currRegTask, j);
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] we have a regular polynomial: " << f << std::endl;
						#endif


					if (currRegTask.chain.dimension() < currTask.chain.dimension()) {
	//				if (regularComponents[j].chain.dimension() < currTask.chain.dimension()) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] dimension drop" << std::endl;
						#endif
						result = RC_REGGCD_OBJ(RecursivePoly(),currRegTask.chain);
	//					result = RC_REGGCD_OBJ(RecursivePoly(),regularComponents[j].chain);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(RecursivePoly(),regularComponents[j].chain));
					}
					else if (currRegTask.chain.reduceMinimal(currRegTask.poly).isZero()) { // this reduce is not needed if regularize reduces output
//					else if (currRegTask.chain.reduce(currRegTask.poly).isZero()) { // this reduce is not needed if regularize reduces output
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] s_" << i << " is zero mod D" << std::endl;
						#endif
						newTask = {currRegTask.chain,i+1};
	//					newTask = {regularComponents[j].chain,i+1};
						tasks.push_back(newTask);
					}
					else if (i==d) {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] case: i==d" << std::endl;
						#endif

						result = RC_REGGCD_OBJ(currRegTask.chain.reduceMinimal(p.mainPrimitivePart()),currRegTask.chain);
//						result = RC_REGGCD_OBJ(currRegTask.chain.reduce(p.mainPrimitivePart()),currRegTask.chain);
						RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//					results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(regularComponents[j].chain.reduce(p.mainPrimitivePart()),regularComponents[j].chain));
					}
					else {
						#ifdef REGULARCHAIN_DEBUG
							std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] case: i<d, i = " << i << std::endl;
						#endif
						f = S.subResultantOfIndex(i); //Jan17/2020: reduction causes coefficients to explode.
						f.setRingVariables(this->allVariables());
						// f = currRegTask.chain.reduceMinimal(S.subResultantOfIndex(i));
						// f = currRegTask.chain.reduce(S.subResultantOfIndex(i));



						#ifdef REGULARCHAIN_PROFILING
							startTimer(&rcProfilingStart);
						#endif
						// Replace with a check that the initial is not equal to the content, and if not return mPP(f)?
						tpp = f.mainPrimitivePart(tc);

						#ifdef REGULARCHAIN_PROFILING
							stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
						#endif
						if ((f.initial().leadingCoefficient() * tc.leadingCoefficient()) < 0)
							tc *= minusOne;
	//					std::cerr << "\n\ntc = " << tc << std::endl;
	//					std::cerr << "f.initial() = " << f.initial() << std::endl;
	//					std::vector<Symbol> vtemp = {Symbol("x")};
	//					tc = f.content(vtemp);
	//					std::cerr << "\nf.content(Symbol(\"x\")) = " << f.content(vtemp) << "\n\n\n\n\n\n\n" << std::endl;
						if (tc != f.initial()) {
							result = RC_REGGCD_OBJ(tpp,currRegTask.chain);
	//						result = RC_REGGCD_OBJ(tpp,regularComponents[j].chain);
							RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//						results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(tpp,regularComponents[j].chain));
						}
						else {
							result = RC_REGGCD_OBJ(f,currRegTask.chain);
	//						result = RC_REGGCD_OBJ(f,regularComponents[j].chain);
							RC_GEN_PRODUCER_ACCUMULATE(results, result);
	//						results.emplace_back(PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>(f,regularComponents[j].chain));
						}
					}
				}
			}
		}
	}
	// NB: Maple code returns two lists, one of the results where there is a dimension drop and one of the regular results
	#if defined(REGULARCHAIN_DEBUG) && !defined(RC_WITH_GENERATORS) && !RC_WITH_GENERATORS
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] results: v = " << v << std::endl;
		for (auto t : results)
			std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] g = " << t.poly << ", T = " << t.chain << std::endl;
		std::cerr << "[rGCD][" << regularGCDDepth << "][" << depth << "] leaving regularGCD: " << p << ", " << q << std::endl;
	#endif
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&regularGCDStart,&regularGCDTime);
	#endif
	--regularGCDDepth;
	--depth;
	RC_GEN_PRODUCER_COMPLETE(results);
//	return results;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::isInRadicalSaturatedIdeal(const RecursivePoly& p) const {
	std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> results;
	results = regularize(p);
	for (size_t i=0; i<results.size(); ++i) {
		if (!results[i].poly.isZero()) {
			return false;
		}
	}
	return true;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::compareCertifiedNoSplit(const RegularChain<Field,RecursivePoly>& rc1, const RegularChain<Field,RecursivePoly>& rc2) {

	if (rc1.dimension() > rc2.dimension())
		return false;

	if (!isSubset(rc2.mainVariables(),rc1.mainVariables()))
		return false;

	std::vector<RecursivePoly> polys1,polys2,initials;
	polys1 = rc1.polynomials();
	polys2 = rc2.polynomials();
	RecursivePoly prodInitials,initial;
	prodInitials.one();
	std::vector<Symbol> ls(rc1.allVariables());
	prodInitials.setRingVariables(ls);
	for (size_t i=0; i<polys1.size(); ++i) {
		if (!polys1[i].isZero()) {
			initial = polys1[i].initial();
			initial.setRingVariables(ls);
			prodInitials *= initial;
		}
	}
//	for (size_t i=0; i<polys1.size(); ++i) {
//		if (!polys1[i].isZero()) {
//			prodInitials *= polys1[i].initial();
//		}
//	}
	for (size_t i=0; i<polys2.size(); ++i) {
		if (!polys2[i].isZero()) {
			initials.emplace_back(polys2[i].initial());
		}
	}

	// Test whether all polynomials of rc2 are in √(sat(rc1))
	for (size_t i=0; i<polys2.size(); ++i) {
		if (!rc1.isInRadicalSaturatedIdeal(polys2[i]))
			return false;
	}

	std::vector<RegularChain<Field,RecursivePoly>> lrc1,lrc2;
	RegularChain<Field,RecursivePoly> T;
	polys2 = polys1;
	polys1.clear();
	for (size_t i=0; i<polys2.size(); ++i) {
		if (!polys2[i].isZero())
			polys1.emplace_back(std::move(polys2[i]));
	}

	// Test whether any intersection of V(rc1) with initials of rc2 is in V(h_{rc1}) (so not in W(rc1))
	lrc1 = T.triangularize(polys1); // decomposition of V(rc1)
	for (size_t i=0; i<lrc1.size(); ++i) {
		for (size_t j=0; j<initials.size(); ++j) {
			lrc2 = lrc1[i].intersect(initials[j]);
			for (size_t k=0; k<lrc2.size(); ++k) {
				if (!lrc2[k].isInRadicalSaturatedIdeal(prodInitials))
					return false;
			}
		}
	}
	return true;
}

template <class Field, class RecursivePoly>
bool RegularChain<Field,RecursivePoly>::compareHeuristicNoSplit(const RegularChain<Field,RecursivePoly>& rc1, const RegularChain<Field,RecursivePoly>& rc2) {

	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "test dim(rc1) > dim(rc2)" << std::endl;
	#endif
	if (rc1.dimension() > rc2.dimension())
		return false;

	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "test ¬(mvar(rc2) ⊆ mvar(rc1))" << std::endl;
	#endif
	if (!isSubset(rc2.mainVariables(),rc1.mainVariables()))
		return false;

	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "test ∀ p ∈ rc2 (p ∈ Sat(rc1))" << std::endl;
	#endif
	// Test whether the polynomials of rc2 are zero in W(rc1)
	std::vector<RecursivePoly> polys;
	polys = rc2.polynomials();
	for (size_t i=0; i<polys.size(); ++i) {
		if (!polys[i].isZero()) {
			if (!rc1.isInSaturatedIdealMinimal(polys[i]))
				return false;
		}
	}

	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "test for p = ∏(init(rc2)), V(p) ∩ W(rc1) = ∅" << std::endl;
	#endif
	// Test whether the initials of the polynomials of rc2 have trivial intersection with W(rc1)
	RecursivePoly prodInitials,initial;
	prodInitials.one();
	std::vector<Symbol> ls(rc1.allVariables());
	prodInitials.setRingVariables(ls);
	for (size_t i=0; i<polys.size(); ++i) {
		if (!polys[i].isZero()) {
			initial = polys[i].initial();
			initial.setRingVariables(ls);
			#ifdef REGULARCHAIN_DEBUGII
				std::cerr << "prodInitials = " << prodInitials << std::endl;
				printVariables(prodInitials.variables(),"prodInitials");
				std::cerr << "initial = " << initial << std::endl;
				printVariables(initial.variables(),"initial");
			#endif
			prodInitials *= initial;
		}
	}
	std::vector<RegularChain<Field,RecursivePoly>> lrc;
	lrc = rc1.intersect(prodInitials,true);
	if (lrc.empty())
		return true;
	else
		return false;
}

template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::removeRedundantChains(const std::vector<RegularChain<Field,RecursivePoly>>& lrc, std::vector<RegularChain<Field,RecursivePoly>>& results) {
	int n = lrc.size();
	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "Entering removeRedundantChains: n = " << n << std::endl;
		for (size_t i=0;i<lrc.size();++i) {
			std::cerr << "lrc[" << i << "] = " << lrc[i] << std::endl;
		}
	#endif
	results.clear();
	std::vector<RegularChain<Field,RecursivePoly>> lrcA,lrcB;
	if (n < 2) {
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "Leaving removeRedundantChains in base case: n = " << n << std::endl;
		#endif
		results = lrc;
		return;
	}
	else if (n == 2) {
		results = mergeIrredundantLists(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+1),std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+1,lrc.end()));
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "Returning to removeRedundantChains from recursive call: n = " << n << std::endl;
			for (size_t i=0;i<results.size();++i) {
				std::cerr << "results[" << i << "] = " << results[i] << std::endl;
			}
		#endif
		return;
	}
	else {
		int a = floor(n/2);
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "a = " << a << std::endl;
		#endif
		std::vector<RegularChain<Field,RecursivePoly>> results1,results2;
		lrcA = std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+a);
		lrcB = std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+a,lrc.end());

#if defined(RC_RRC_BY_CILK) && RC_RRC_BY_CILK
		cilk_spawn
		removeRedundantChains(lrcA,results1);
		removeRedundantChains(lrcB,results2);
		cilk_sync;
#elif !(defined(SERIAL) && SERIAL) && (defined(RC_RRC_PARALLEL) && RC_RRC_PARALLEL)
		std::vector<ExecutorThreadPool::threadID> workers;
		int nthreads = ExecutorThreadPool::getThreadPool().obtainThreads(1, workers);
		if (nthreads > 0) {
			std::function<void()> task = std::bind((&RegularChain<Field, RecursivePoly>::removeRedundantChains), std::ref(lrcA), std::ref(results1));
			ExecutorThreadPool::getThreadPool().executeTask(workers[0], task);
			removeRedundantChains(std::ref(lrcB),std::ref(results2));
			ExecutorThreadPool::getThreadPool().waitForThreads(workers);
			ExecutorThreadPool::getThreadPool().returnThreads(workers);
		} else {
			// cilk_spawn
			removeRedundantChains(lrcA,results1);
			removeRedundantChains(lrcB,results2);
			// cilk_sync;
		}
#else
		removeRedundantChains(lrcA,results1);
		removeRedundantChains(lrcB,results2);
#endif

		results = mergeIrredundantLists(results1,results2);
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "Returning to removeRedundantChains from recursive call: n = " << n << std::endl;
			for (size_t i=0;i<results.size();++i) {
				std::cerr << "results[" << i << "] = " << results[i] << std::endl;
			}
		#endif
		return;
	}
}

//template <class Field, class RecursivePoly>
//std::vector<RegularChain<Field,RecursivePoly>> RegularChain<Field,RecursivePoly>::removeRedundantChains(const std::vector<RegularChain<Field,RecursivePoly>>& lrc) {
//	int n = lrc.size();
//	#ifdef REGULARCHAIN_DEBUGII
//		std::cerr << "Entering removeRedundantChains: n = " << n << std::endl;
//		for (size_t i=0;i<lrc.size();++i) {
//			std::cerr << "lrc[" << i << "] = " << lrc[i] << std::endl;
//		}
//		std::vector<RegularChain<Field,RecursivePoly>> ret;
//	#endif
//	if (n < 2) {
//		#ifdef REGULARCHAIN_DEBUGII
//			std::cerr << "Leaving removeRedundantChains in base case: n = " << n << std::endl;
//		#endif
//		return lrc;
//	}
//	else if (n == 2) {
////		return mergeIrredundantLists(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+1),std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+1,lrc.end()));
//		#ifdef REGULARCHAIN_DEBUGII
//			ret = mergeIrredundantLists(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+1),std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+1,lrc.end()));
//			std::cerr << "Returning to removeRedundantChains from recursive call: n = " << n << std::endl;
//			for (size_t i=0;i<ret.size();++i) {
//				std::cerr << "ret[" << i << "] = " << ret[i] << std::endl;
//			}
//			return ret;
//		#else
//			return mergeIrredundantLists(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+1),std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+1,lrc.end()));
//		#endif
//	}
//	else {
//		int a = floor(n/2);
//		#ifdef REGULARCHAIN_DEBUGII
//			std::cerr << "a = " << a << std::endl;
//		#endif
//		#ifdef REGULARCHAIN_DEBUGII
//			ret = mergeIrredundantLists(removeRedundantChains(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+a)),removeRedundantChains(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+a,lrc.end())));
//			std::cerr << "Returning to removeRedundantChains from recursive call: n = " << n << std::endl;
//			for (size_t i=0;i<ret.size();++i) {
//				std::cerr << "ret[" << i << "] = " << ret[i] << std::endl;
//			}
//			return ret;
//		#else
//			return mergeIrredundantLists(removeRedundantChains(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin(),lrc.begin()+a)),removeRedundantChains(std::vector<RegularChain<Field,RecursivePoly>>(lrc.begin()+a,lrc.end())));
//		#endif
//	}
//}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> oneSideCompare(const RegularChain<Field,RecursivePoly>& rc1, const RegularChain<Field,RecursivePoly>& rc2) {
	std::vector<RegularChain<Field,RecursivePoly>> results;
	bool inclusion = RegularChain<Field,RecursivePoly>::compareHeuristicNoSplit(rc1,rc2);
//	bool inclusion = RegularChain<Field,RecursivePoly>::compareCertifiedNoSplit(rc1,rc2);
	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "inclusion? : " << inclusion << std::endl;
	#endif
//	if (!RegularChain<Field,RecursivePoly>::compareCertifiedNoSplit(rc1,rc2)) {
//	if (!RegularChain<Field,RecursivePoly>::compareHeuristicNoSplit(rc1,rc2)) {
	if (!inclusion) {
		results.emplace_back(rc1);
	}
	return results;
}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> oneSideCompare(const std::vector<RegularChain<Field,RecursivePoly>>& lrc, const RegularChain<Field,RecursivePoly>& rc){
	// SynchronizedWriteVector<RegularChain<Field,RecursivePoly>> results;
	std::vector<RegularChain<Field,RecursivePoly>> results;
	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "entering oneSideCompare(lrc,rc)" << std::endl;
		RegularChain<Field,RecursivePoly> rct;
	#endif

	results.reserve(lrc.size());
	for (size_t i=0; i<lrc.size(); ++i) {
	// cilk_for (size_t i=0; i<lrc.size(); ++i) {
		std::vector<RegularChain<Field,RecursivePoly>> interResults;
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "comparing" << std::endl;
			rct = lrc[i];
			std::cerr << rct << std::endl;
			std::cerr << "and" << std::endl;
			rct = rc;
			std::cerr << rct << std::endl;
		#endif
		interResults = oneSideCompare(lrc[i],rc);
		// for (size_t j=0; j<interResults.size(); ++j) {
		// 	results.push_back(interResults[j]);
		// }
		results.insert(results.end(),interResults.begin(),interResults.end());
	}
	return results;
	// return results.moveVectorOut();
}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> oneSideCompare(const std::vector<RegularChain<Field,RecursivePoly>>& lrc1, const std::vector<RegularChain<Field, RecursivePoly>>& lrc2){
	SynchronizedWriteVector<RegularChain<Field,RecursivePoly>> results;
	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "entering oneSideCompare(lrc1,lrc2)" << std::endl;
		RegularChain<Field,RecursivePoly> rc;
	#endif

	results.reserve(lrc1.size()+lrc2.size());
#if defined(RC_RRC_BY_CILK) && RC_RRC_BY_CILK
	cilk_for (size_t i=0; i<lrc1.size(); ++i) {
#else
	for (size_t i=0; i<lrc1.size(); ++i) {
#endif
		std::vector<RegularChain<Field,RecursivePoly>> rc1results;
		rc1results.clear();
		rc1results.emplace_back(lrc1[i]);
		#ifdef REGULARCHAIN_DEBUGII
			for (size_t k=0; k<rc1results.size(); ++k) {
				rc = rc1results[k];
				std::cerr << "rc1results[" << k << "] = " << rc << std::endl;
			}
		#endif
		for (size_t j=0; j<lrc2.size(); ++j) {
			#ifdef REGULARCHAIN_DEBUGII
				std::cerr << "comparing rc1results and" << std::endl;
				rc = lrc2[j];
				std::cerr << "lrc2[" << j << "] = " << rc << std::endl;
			#endif
			rc1results = oneSideCompare(rc1results,lrc2[j]);
		}
		for (size_t j=0; j<rc1results.size(); ++j) {
			results.push_back(rc1results[j]);
		}
//		results.insert(results.end(),rc1results.begin(),rc1results.end());
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "end outer for loop in oneSideCompare(lrc1,lrc2)" << std::endl;
			for (size_t k=0; k<results.size(); ++k) {
				rc = results[k];
				std::cerr << "results[" << k << "] = " << rc << std::endl;
			}
		#endif
	}
	return results.moveVectorOut();
//	return results;
}

template <class Field, class RecursivePoly>
void oneSideCompare_Task(const std::vector<RegularChain<Field,RecursivePoly>>& lrc1, const std::vector<RegularChain<Field, RecursivePoly>>& lrc2, int k1, int k2,
	SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>& results)
{
	#ifdef REGULARCHAIN_DEBUGII
		std::cerr << "entering oneSideCompare(lrc1,lrc2)" << std::endl;
		RegularChain<Field,RecursivePoly> rc;
	#endif

	for (size_t i=k1; i<k2; ++i) {
		std::vector<RegularChain<Field,RecursivePoly>> rc1results;
		rc1results.clear();
		rc1results.emplace_back(lrc1[i]);
		#ifdef REGULARCHAIN_DEBUGII
			for (size_t k=0; k<rc1results.size(); ++k) {
				rc = rc1results[k];
				std::cerr << "rc1results[" << k << "] = " << rc << std::endl;
			}
		#endif
		for (size_t j=0; j<lrc2.size(); ++j) {
			#ifdef REGULARCHAIN_DEBUGII
				std::cerr << "comparing rc1results and" << std::endl;
				rc = lrc2[j];
				std::cerr << "lrc2[" << j << "] = " << rc << std::endl;
			#endif
			rc1results = oneSideCompare(rc1results,lrc2[j]);
		}
		for (size_t j=0; j<rc1results.size(); ++j) {
			results.push_back(rc1results[j]);
		}
//		results.insert(results.end(),rc1results.begin(),rc1results.end());
		#ifdef REGULARCHAIN_DEBUGII
			std::cerr << "end outer for loop in oneSideCompare(lrc1,lrc2)" << std::endl;
			for (size_t k=0; k<results.size(); ++k) {
				rc = results[k];
				std::cerr << "results[" << k << "] = " << rc << std::endl;
			}
		#endif
	}
}

template <class Field, class RecursivePoly>
std::vector<RegularChain<Field,RecursivePoly>> mergeIrredundantLists(const std::vector<RegularChain<Field,RecursivePoly>>& lrc1, const std::vector<RegularChain<Field,RecursivePoly>>& lrc2) {

	#ifdef REGULARCHAIN_DEBUGII
		RegularChain<Field,RecursivePoly> rc;
		std::cerr << "mergeIrredundantLists:" << std::endl;
		for (size_t i=0; i<lrc1.size(); ++i) {
			rc = lrc1[i];
			std::cerr << "lrc1[" << i << "] = " << lrc1[i] << std::endl;
		}
		for (size_t i=0; i<lrc2.size(); ++i) {
			rc = lrc2[i];
			std::cerr << "lrc2[" << i << "] = " << lrc2[i] << std::endl;
		}
	#endif

#if !(defined(RC_RRC_BY_CILK) && RC_RRC_BY_CILK) && !(defined(SERIAL) && SERIAL) && (defined(RC_RRC_PARALLEL) && RC_RRC_PARALLEL)
	SynchronizedWriteVector<RegularChain<Field,RecursivePoly>> results1, results2;

	//do comparison of lrc1, lrc2
	results1.reserve(lrc1.size() + lrc2.size());
	std::vector<ExecutorThreadPool::threadID> workers;
	int size = lrc1.size();
	int nthreads = ExecutorThreadPool::getThreadPool().obtainThreads(size-1, workers);
	int k1 = 0, k2 = 0;
	for (int i = 0; i < nthreads; ++i) {
    	k2 = (size / (nthreads+1)); //+1 becuase current thread does work to
    	k2 += k1; //add offset
    	std::function<void()> task = std::bind(oneSideCompare_Task<Field,RecursivePoly>, std::ref(lrc1), std::ref(lrc2), k1, k2, std::ref(results1));
		ExecutorThreadPool::getThreadPool().executeTask(workers[i], task);
		k1 = k2;
	}
	k2 = size;
	oneSideCompare_Task(lrc1, lrc2, k1, k2, results1);
	ExecutorThreadPool::getThreadPool().waitForThreads(workers);

	//do comparison or lrc2, results1
	size = lrc2.size();
	nthreads = (size-1) < nthreads ? size-1 : nthreads;
	results2.reserve(lrc2.size() + results1.size());
	k1 = 0; k2 = 0;
	for (int i = 0; i < nthreads; ++i) {
    	k2 = (size / (nthreads+1)); //+1 becuase current thread does work to
    	k2 += k1; //add offset
    	std::function<void()> task = std::bind(oneSideCompare_Task<Field,RecursivePoly>, std::ref(lrc2), std::ref(results1.vector()), k1, k2, std::ref(results2));
		ExecutorThreadPool::getThreadPool().executeTask(workers[i], task);
		k1 = k2;
	}
	k2 = size;
	oneSideCompare_Task(lrc2, results1.vector(), k1, k2, results2);
	ExecutorThreadPool::getThreadPool().waitForThreads(workers);

	//done with threads.
	ExecutorThreadPool::getThreadPool().returnThreads(workers);

	std::vector<RegularChain<Field,RecursivePoly>> results = results2.moveVectorOut();
	results.insert(results.end(), results1.begin(), results1.end());
	return results;

#else
	std::vector<RegularChain<Field,RecursivePoly>> result1,result2;

	result1 = oneSideCompare(lrc1,lrc2); 	// determine components of lrc1 not contained in components of lrc2
	#ifdef REGULARCHAIN_DEBUGII
		for (size_t i=0; i<result1.size(); ++i) {
			rc = result1[i];
			std::cerr << "result1[" << i << "] = " << result1[i] << std::endl;
		}
	#endif
	result2 = oneSideCompare(lrc2,result1);	// determine components of lrc2 not contained in remaining components of lrc1
	#ifdef REGULARCHAIN_DEBUGII
		for (size_t i=0; i<result2.size(); ++i) {
			rc = result2[i];
			std::cerr << "result2[" << i << "] = " << result2[i] << std::endl;
		}
	#endif

	result2.reserve(result1.size()+result2.size());
	result2.insert(result2.end(),result1.begin(),result1.end());
	return result2;
#endif
}


/**
 * Generate a random regular chain
 *
 * @param nVars: number of variables
 * @param nAlgVars: number of algebraic variables
 * @param nTrcVars: number of transcendental variables
 * @param nTerms: maximum number of terms in the polynomials
 * @param coefBound: maximum coefficient size
 * @param pSparsity: sparsity of the polynomials
 * @param includeNeg: whether to include negative coefficients
 **/
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::randomRegularChain(int nVars, int nAlgVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg) {

//	bool check = true;
//	bool notRegular = true;
//	std::vector<Symbol> names;
//	RecursivePoly ex(Symbol::randomElement()),exn,temp,one;
//	one.one();
//	mode = TS_FIXED;
//	set.clear();
//	set.reserve(nVars);
//	for (size_t i=0; i<nVars; ++i) {
//		set.emplace_back(); // initialize elements of set with empty polynomials
//	}
//	vars.clear();
//	algVars.clear();
//	trcVars.clear();
//	RecursivePoly p;
//	std::vector<Symbol> vs,vs2,variables,algebraicVariables,transcendentalVariables;
//	std::vector<int> indices;
//	int nvs(nVars+nTrcVars);
//	vs = Symbol::randomElements(nvs);
//	variables.insert(variables.end(),vs.begin(),vs.begin()+nVars);
//	transcendentalVariables.insert(transcendentalVariables.end(),vs.begin()+nVars,vs.begin()+nvs);
//	// Generate random indicies for the algebraic variables
//	indices = randValsInRange(0,nVars-1,nAlgVars);
//	sort(indices.begin(),indices.end());
////	for (auto i=0; i<nAlgVars; ++i) {
////		algebraicVariables.push_back(variables[indices[i]]);
////	}
//	vars = variables;
//	printVariables(vars,"vars");
//	for (size_t i=0; i<indices.size(); ++i)
//		std::cerr << "indices[" << i << "] = " << indices[i] << std::endl;
//	algebraicVariables.push_back(variables[indices.back()]);
//	algVars = algebraicVariables;
//	printVariables(algVars,"algVars");
//	trcVars = transcendentalVariables;
//	printVariables(trcVars,"trcVars");
//	p = TriangularSet<Field,RecursivePoly>::randomTriangularSetPolynomial(variables,indices.back(),transcendentalVariables,nTerms,coefBound,pSparsity,includeNeg);
//	names.clear();
//	names.push_back(variables[indices.back()]);
//	ex.setRingVariables(names);
//	exn = ex^(p.leadingVariableDegree().get_si()+2);
//	temp = one;
//	temp *= exn;
//	p += temp;
//	std::cerr << "p = " << p << std::endl;
//	set[indices.back()] = p;
//	updateTriangularSetStates(p);
//	this->display();
//	std::vector<RecursivePoly> polys;
//	polys = this->polynomials();
//	for (size_t i=0; i<polys.size(); ++i) {
//		std::cerr << "polys[" << i << "] = " << polys[i] << std::endl;
//	}
////	updateTriangularSetStates(p);
//	for (auto i=indices.size()-2; i>=0; --i) {
//		// create a random polynomial in only variables up to the current index, and check to see if it is regular
//		// wrt to the emerging regular chain.
//		while (notRegular) {
//			p = TriangularSet<Field,RecursivePoly>::randomTriangularSetPolynomial(variables,indices[i],transcendentalVariables,nTerms,coefBound,pSparsity,includeNeg);
//			names.clear();
//			names.push_back(p.leadingVariable());
//			ex.setRingVariables(names);
//			exn = ex^(p.leadingVariableDegree().get_si()+2);
//			temp = one;
//			temp *= exn;
//			p += temp;
//			if (this->isRegular(p))
//				notRegular = false;
//		}
//		set[indices[i]] = p;
//		algebraicVariables.push_back(variables[indices[i]]);
//		algVars = algebraicVariables;
//		updateTriangularSetStates(p);
//	}
//	check = !stronglyNormalized;

	bool check = true;
	while (check) {
		std::vector<Symbol> names;
		RecursivePoly ex(Symbol::randomElement()),exn,temp,one;
		one.one();
		mode = TS_FIXED;
		set.clear();
		set.reserve(nVars);
		for (int i=0; i<nVars; ++i) {
			set.emplace_back(); // initialize elements of set with empty polynomials
		}
		vars.clear();
		algVars.clear();
		trcVars.clear();
		RecursivePoly p;
		std::vector<Symbol> vs,variables,algebraicVariables,transcendentalVariables;
		std::vector<int> indices;
		int nvs(nVars+nTrcVars);
		vs = Symbol::randomElements(nvs);
		variables.insert(variables.end(),vs.begin(),vs.begin()+nVars);
		transcendentalVariables.insert(transcendentalVariables.end(),vs.begin()+nVars,vs.begin()+nvs);
		// Generate random indicies for the algebraic variables
		indices = randValsInRange(0,nVars-1,nAlgVars);
		sort(indices.begin(),indices.end());
		for (int i=0; i<nAlgVars; ++i) {
			algebraicVariables.emplace_back(variables[indices[i]]);
		}
		for (size_t i=0; i<indices.size(); ++i) {
			p = TriangularSet<Field,RecursivePoly>::randomTriangularSetPolynomial(variables,indices[i],transcendentalVariables,nTerms,coefBound,pSparsity,includeNeg);
			names.clear();
			names.emplace_back(variables[indices[i]]);
			ex.setRingVariables(names);
			exn = ex^(p.leadingVariableDegree().get_si()+2);
			temp = one;
			temp *= exn;
			p += temp;
			set[indices[i]] = p;
			updateTriangularSetStates(p);
		}
		algVars = algebraicVariables;
		trcVars = transcendentalVariables;
		vars = variables;
		check = !stronglyNormalized;
	}
}

/**
 * Generate a random regular chain
 *
 * @param nVars: number of variables
 * @param nAlgVars: number of algebraic variables
 * @param nTrcVars: number of transcendental variables
 * @param maxDegs: maximum degrees among the full set of variables
 * @param coefBound: maximum coefficient size
 * @param pSparsity: sparsity of the polynomials
 * @param includeNeg: whether to include negative coefficients
 **/
template <class Field, class RecursivePoly>
void RegularChain<Field,RecursivePoly>::randomRegularChain(int nVars, int nAlgVars, int nTrcVars, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg) {
	if (size_t(nVars + nTrcVars) != maxDegs.size()) {
		std::cerr << "BPAS: error, total number of variables must equal the size of the vector of maximum degrees." << std::endl;
		exit(1);
	}
	bool check = true;
	while (check) {
		std::vector<Symbol> names;
		RecursivePoly exn,temp,one;
		one.one();
//		mode = TS_FIXED;
		set.clear();
		set.reserve(nVars);
		for (int i=0; i<nVars; ++i) {
			set.emplace_back(); // initialize elements of set with empty polynomials
		}
		vars.clear();
		algVars.clear();
		trcVars.clear();
		RecursivePoly p;
		std::vector<Symbol> vs,variables,algebraicVariables,transcendentalVariables;
		std::vector<int> indices;
		int nvs(nVars+nTrcVars);
		vs = Symbol::randomElements(nvs);
		variables.insert(variables.end(),vs.begin(),vs.begin()+nVars);
		transcendentalVariables.insert(transcendentalVariables.end(),vs.begin()+nVars,vs.begin()+nvs);
		// Generate random indicies for the algebraic variables
		indices = randValsInRange(0,nVars-1,nAlgVars);
		sort(indices.begin(),indices.end());
		for (int i=0; i<nAlgVars; ++i) {
			algebraicVariables.emplace_back(variables[indices[i]]);
		}
		for (size_t i=0; i<indices.size(); ++i) {
			p = TriangularSet<Field,RecursivePoly>::randomTriangularSetPolynomial(variables,indices[i],transcendentalVariables,maxDegs,coefBound,pSparsity,includeNeg);
			names.clear();
			names.emplace_back(variables[indices[i]]);
			set[indices[i]] = p;
			updateTriangularSetStates(p);
		}
		algVars = algebraicVariables;
		trcVars = transcendentalVariables;
		vars = variables;
		check = !stronglyNormalized;
	}
}


/// Possible instantiations of the RegularChain class ///
template class RegularChain<RationalNumber,SparseMultivariateRationalPolynomial>;
