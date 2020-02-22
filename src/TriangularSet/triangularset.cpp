#include "../../include/RegularChain/regularchain.hpp"
#include "../../include/TriangularSet/triangularset.hpp"
#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include "../../include/ring.h"
#include "../../include/Utils/RandomHelpers.hpp"
#include "../../include/Utils/SymbolHelpers.hpp"
#include <iostream>


/// Protected functions //

template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::updateTriangularSetStates(const RecursivePoly& p) {

	//Check whether ts is strongly normalized and how much of it is strongly normalized (after adding p to ts)
	RecursivePoly init = p.initial();
	Symbol v(p.leadingVariable());
	int vIndex = std::find(algVars.begin(),algVars.end(),v) - algVars.begin();
	int sNIndex = std::find(algVars.begin(),algVars.end(),sNMaxVar) - algVars.begin();
	bool vIsAtEnd,initIsConstant;
	
	if (vIndex == algVars.size()-1)
		vIsAtEnd = true;
	else
		vIsAtEnd = false;
	
	initIsConstant = isConstantPolynomial(init);
		
	if (stronglyNormalized){
		if (!initIsConstant) {	// ts no longer stronglyNormalized
			stronglyNormalized = false;
		}
		else {							// ts continues to be stronglyNormalized, update sNMaxVar
			if (vIndex < sNIndex) {
				sNMaxVar = v;
			}
			return;
		}
	}
	// determine if sNMaxVar needs updating and update it if required
	// TODO: write a version of this that takes into account known information
//	if (vIndex > sNIndex) {
//		if (vIsAtEnd) {
//			sNMaxVar = Symbol("");
//		}
//		else {
//			sNMaxVar = algVars[vIndex+1];
//		}
//	}
	updateTriangularSetStates();
	
	//Check whether or not by adding p to ts, ts remains strongly normalized
//	RecursivePoly init = p.initial();
//	std::vector<Symbol> initVars = init.variables();
////	std::cerr << "p.initial() = " << init << std::endl;
////	std::cerr << "p.initial() variables:" << std::endl;
////	printVariables(initVars);
//	if(stronglyNormalized){
//		std::vector<Symbol> commonVars;
//		// Tests to ensure that the initial of p has no transcendental variables		//
//		// TODO: This test will be replaced by a more refined one in a future version 	//
//		//		 to allow for transcendental initials, which require working over the 	//
//		//		 field of rational functions of the transcendental variables.			//
//		commonVars = setIntersection(trcVars,initVars);
//		if (!commonVars.empty())
//			stronglyNormalized = false;
//		// Tests to ensure that the initial of p has no potentially algebraic variables	//
//		// Test is robust against additions of polynomials to a ts 						//
//		commonVars = setIntersection(vars,initVars);
//		if (!commonVars.empty())
//			stronglyNormalized = false;
////		std::cerr << "TriangularSet is ";
////		if (!stronglyNormalized)
////			std::cerr << "not ";
////		std::cerr << "stronglyNormalized" << std::endl;
//	}
}

template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::updateTriangularSetStates() {

	//Check whether or not ts is strongly normalized (after polynomials have been removed)
	if (algVars.empty()) {
		stronglyNormalized = true;
		sNMaxVar = Symbol("");
		return;
	}
	RecursivePoly init;
	Symbol v;
	int i,index;
	bool maxFound = false;
	if (!stronglyNormalized) {
		for (i=set.size()-1; i>-1; --i) {
		//TODO: can be improved to start loop based on the index of the current value of sNMaxVar
			if (!set[i].isZero()) {
				if (isConstantPolynomial(set[i].initial())) {
					sNMaxVar = vars[i];
				}
				else {
					maxFound = true;
					break;
				}
			}
		}
		if (!maxFound) {
			stronglyNormalized = true;
		}
	}
	else { // check if top poly was removed, meaning that sNMaxVar is no longer an algebraic variable
		if (sNMaxVar != algVars[0]) {
			sNMaxVar = algVars[0];
		}
	}

	//Check whether or not ts is strongly normalized (after polynomials have been removed)
//	std::vector<Symbol> initVars,commonVars;
//	RecursivePoly init;
//	int i(0);
//	if (!stronglyNormalized) {
//		while(stronglyNormalized || i<size()) {
//			if (!set[i].isZero()) {
//				// Tests to ensure that the initial of p has no transcendental variables		//
//				// TODO: This test will be replaced by a more refined one in a future version 	//
//				//		 to allow for transcendental initials, which require working over the 	//
//				//		 field of rational functions of the transcendental variables.			//
//				init = set[i].initial();
//				initVars = init.variables();
//				if (!trcVars.empty()) {
//					commonVars = setIntersection(trcVars,initVars);
//					if (!commonVars.empty())
//						stronglyNormalized = false;
//				}
//				// Tests to ensure that the initial of p has no potentially algebraic variables	//
//				commonVars = setIntersection(vars,initVars);
//				if (!commonVars.empty())
//					stronglyNormalized = false;
//			}
//			++i;
//		}
//	}
}

template <class Field, class RecursivePoly>	
int TriangularSet<Field,RecursivePoly>::variableIndex(const Symbol& s) const {
	int a = std::find(vars.begin(),vars.end(),s) - vars.begin();
	return a;
}

/// Constructors ///


/**
 * Default constructor: creates an empty triangular set of variable size
 * with empty list of transcendentals
 *
 * @param
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet () : mode(TS_VARIABLE), characteristic(), stronglyNormalized(true), sNMaxVar("") {
	Field f;
	characteristic = f.characteristic;
}

/**
 * Construct an empty triangular set of fixed size in the s decreasingly ordered 
 * variables given by xs with empty list of transcendentals
 *
 * @param s: The number of polynomials
 * @param xs: The variable names
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (const std::vector<Symbol>& xs) : mode(TS_FIXED), characteristic(), stronglyNormalized(true), sNMaxVar("") {
	if (xs.size() == 0) {
			std::cerr << "BPAS: error, Fixed mode TriangularSet cannot be defined with 0 variables." << std::endl;
			exit(1);
	}
	Field f;
	characteristic = f.characteristic;
	vars = xs;
	set.reserve(xs.size());
	set.clear();
	for (int i=0; i<xs.size(); ++i) {
		set.emplace_back(); // initialize elements of set with empty polynomials
	}
}

/**
 * Construct an empty triangular set of fixed size in the s decreasingly ordered 
 * variables given by xs and list of transcendentals given by ts
 *
 * @param s: The number of polynomials
 * @param xs: The variable names
 * @param ps: The transcendental variable names
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (const std::vector<Symbol>& xs, const std::vector<Symbol>& ps) : mode(TS_FIXED), characteristic(), stronglyNormalized(true), sNMaxVar("") {
	if (xs.size() == 0) {
			std::cerr << "BPAS: error, Fixed mode TriangularSet cannot be defined with 0 variables." << std::endl;
			exit(1);
	}
	Field f;
	characteristic = f.characteristic;
	vars = xs;
	trcVars = ps;
	set.reserve(xs.size());
	set.clear();
	for (int i=0; i<xs.size(); ++i) {
		set.emplace_back(); // initialize elements of set with empty polynomials
	}
}

/**
 * Construct a variable triangular set containing p, such that the variables of p are treated as algebraic,
 * with empty list of transcendentals
 *
 * @param p: The polynomial to add
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (const RecursivePoly& p) : mode(TS_VARIABLE), characteristic(), stronglyNormalized(true), sNMaxVar("") {
	if (p.isConstant()) {
			std::cerr << "BPAS: error, cannot add a constant to a TriangularSet." << std::endl;
			exit(1);
	}
	Field f;
	characteristic = f.characteristic;
	vars = p.variables();
	set.reserve(vars.size());
	set.clear();
	int pIndex = variableIndex(p.leadingVariable());
	for (int i=0; i<vars.size(); ++i) {
		if (i != pIndex)
			set.emplace_back(); // initialize elements of set with empty polynomials
		else {
			SparseMultivariateRationalPolynomial minimizedP = p;
			minimizedP.setRingVariables(vars);
			set.emplace_back(minimizedP); // initialize pIndexth element using copy constructor on p
		}
	}
	algVars.reserve(vars.size());
	algVars.clear();
	algVars.push_back(p.leadingVariable());
//	stronglyNormalized = true;
	updateTriangularSetStates(p);
}

/**
 * Construct a variable triangular set containing p, such that the variables in ps are
 * treated as transcendental, while any remaining variables of p are treated as algebraic
 *
 * @param p: The polynomial to add
 * @param ps: The transcendental variable names
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (const RecursivePoly& p, const std::vector<Symbol>& ps) : mode(TS_VARIABLE), characteristic(), stronglyNormalized(true), sNMaxVar("") {
	if (p.isConstant()) {
			std::cerr << "BPAS: error, cannot add a constant to a TriangularSet." << std::endl;
			exit(1);
	}
	Symbol s(p.leadingVariable());
	if (std::find(ps.begin(),ps.end(),s) != ps.end()) {
			std::cerr << "BPAS: error in TriangularSet: leading variable of input polynomial cannot be transcendental." << std::endl;
			exit(1);
	}
	Field f;
	characteristic = f.characteristic;
	trcVars = ps;
	vars = orderPreservingSetDifference(p.variables(),ps);
	set.reserve(vars.size());
	set.clear();
	int pIndex = variableIndex(s);
	for (int i=0; i<vars.size(); ++i) {
		if (i != pIndex)
			set.emplace_back(); // initialize elements of set with empty polynomials
		else {
			SparseMultivariateRationalPolynomial minimizedP = p;
			minimizedP.setRingVariables(minimizedP.variables());
			set.emplace_back(minimizedP); // initialize pIndexth element using copy constructor on p
		}
	}
	algVars.reserve(vars.size());
	algVars.clear();
	algVars.push_back(s);
//	stronglyNormalized = true;
	updateTriangularSetStates(p);
}

/**
 * Copy constructor
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (const TriangularSet<Field,RecursivePoly>& a) : mode(a.mode), characteristic(a.characteristic), stronglyNormalized(a.stronglyNormalized), sNMaxVar(a.sNMaxVar) {
	unsigned long long int timerStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&timerStart);
	#endif
	if (a.vars.size() || a.trcVars.size()) {
		set = a.set;
		vars = a.vars;
		algVars = a.algVars;
		trcVars = a.trcVars;
	}
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&timerStart,&tsCopyTime);
	#endif
}

///**
// * Copy constructor
// *
// * @param a: A regular chain
// **/
//template <class Field, class RecursivePoly>
//TriangularSet<Field,RecursivePoly>::TriangularSet (const RegularChain<Field,RecursivePoly>& a) : mode(a.mode), characteristic(a.characteristic), stronglyNormalized(a.stronglyNormalized) {
//	TriangularSet<Field,RecursivePoly> ts(
//	if (a.vars.size()) {
//		set = a.set;
//		vars = a.vars;
//		algVars = a.algVars;
//		trcVars = a.trcVars;
//	}
//}

///**
// * Copy constructor
// *
// * @param a: A zero dimensional regular chain
// **/
//template <class Field, class RecursivePoly>
//TriangularSet<Field,RecursivePoly>::TriangularSet (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) : mode(a.mode), characteristic(a.characteristic), stronglyNormalized(a.stronglyNormalized) {
//	if (a.vars.size()) {
//		set = a.set;
//		vars = a.vars;
//		algVars = a.algVars;
//		trcVars = a.trcVars;
//	}
//}

/**
 * Move constructor
 *
 * @param a: An r-value reference triangular set
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (TriangularSet<Field,RecursivePoly>&& a) : mode(a.mode), characteristic(a.characteristic), stronglyNormalized(a.stronglyNormalized), sNMaxVar(a.sNMaxVar) {
	set = std::move(a.set);
	vars = std::move(a.vars);
	algVars = std::move(a.algVars);
	trcVars = std::move(a.trcVars);
	a.mode = TS_VARIABLE;
	a.set.clear();
	a.vars.clear();
	a.algVars.clear();
	a.trcVars.clear();
	a.stronglyNormalized = true;
}

///**
// * Move constructor
// *
// * @param a: An r-value reference regular chain
// **/
//template <class Field, class RecursivePoly>
//TriangularSet<Field,RecursivePoly>::TriangularSet (RegularChain<Field,RecursivePoly>&& a) : mode(a.mode), characteristic(a.characteristic), stronglyNormalized(a.stronglyNormalized) {
//	set = a.set;
//	vars = a.vars;
//	algVars = a.algVars;
//	trcVars = a.trcVars;
//	a.mode = TS_VARIABLE;
//	a.set.clear();
//	a.vars.clear();
//	a.algVars.clear();
//	a.trcVars.clear();
//	a.stronglyNormalized = true;
//}

///**
// * Move constructor
// *
// * @param a: An r-value reference regular chain
// **/
//template <class Field, class RecursivePoly>
//TriangularSet<Field,RecursivePoly>::TriangularSet (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a) : mode(a.mode), characteristic(a.characteristic), stronglyNormalized(a.stronglyNormalized) {
//	set = a.set;
//	vars = a.vars;
//	algVars = a.algVars;
//	trcVars = a.trcVars;
//	a.mode = TS_VARIABLE;
//	a.set.clear();
//	a.vars.clear();
//	a.algVars.clear();
//	a.trcVars.clear();
//	a.stronglyNormalized = true;
//}

/**
 * Computational constructor: creates a triangular set given all the data
 *
 * @param
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::TriangularSet (const std::vector<Symbol>&& vs, const std::vector<Symbol>&& avs, const std::vector<Symbol>&& tvs, const std::vector<RecursivePoly>&& ts, TriangularSetMode tsm, const mpz_class& c) : mode(tsm), characteristic(c), stronglyNormalized(false), sNMaxVar("") {
	unsigned long long int timerStart;
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&timerStart);
	#endif
	// TODO: Check that the set of mainVariables of set equals algVars
	vars = std::move(vs);
	algVars = std::move(avs);
	trcVars = std::move(tvs);
	set = std::move(ts);
	updateTriangularSetStates(); // stronglyNormalized set to false so this routine works properly
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&timerStart,&tsCopyTime);
	#endif
}

/**
 * Deconstructor
 *
 * @param
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>::~TriangularSet() {}

/// Member Functions ///
		
/**
 * Tests if the TriangularSet is empty
 *
 * @param
 **/
template <class Field, class RecursivePoly>
bool TriangularSet<Field,RecursivePoly>::isEmpty() const {
	return algVars.empty();
}
		
/**
 * Tests if the polynomial is constant relative to the TriangularSet
 *
 * @param
 **/
template <class Field, class RecursivePoly>
bool TriangularSet<Field,RecursivePoly>::isConstantPolynomial(const RecursivePoly& p) const {
	if (p.isConstant() != 0)
		return true;
	std::vector<Symbol> vs;
	vs = setDifference(p.variables(),trcVars);
	if (vs.empty())
		return true;
	else 
		return false;
}
		
/**
 * Assignment operator =
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>& TriangularSet<Field,RecursivePoly>::operator= (const TriangularSet<Field,RecursivePoly>& a) {
	if (this != &a) {
		mode = a.mode;
		characteristic = a.characteristic;
		set = a.set;
		vars = a.vars;
		algVars = a.algVars;
		trcVars = a.trcVars;
		stronglyNormalized = a.stronglyNormalized;
		sNMaxVar = a.sNMaxVar;
	}
	return *this;
}
		
/**
 * Assignment operator =
 *
 * @param a: A BPASTriangularSet
 **/
template <class Field, class RecursivePoly>
BPASTriangularSet<Field,RecursivePoly>& TriangularSet<Field,RecursivePoly>::operator= (const BPASTriangularSet<Field,RecursivePoly>& a) {
	if (dynamic_cast<const TriangularSet<Field,RecursivePoly>*>(&a))
		*this = dynamic_cast<const TriangularSet<Field,RecursivePoly>&>(a);
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to TriangularSet."));
	return *this;
}
		
/**
 * Move assignment operator =
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>& TriangularSet<Field,RecursivePoly>::operator= (TriangularSet<Field,RecursivePoly>&& a) {
	if (this != &a) {
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "moving..." << std::endl;
		#endif
		mode = a.mode;
		characteristic = a.characteristic;
		set = std::move(a.set);
		vars = std::move(a.vars);
		algVars = std::move(a.algVars);
		trcVars = std::move(a.trcVars);
		stronglyNormalized = a.stronglyNormalized;
		sNMaxVar = a.sNMaxVar;
		a.mode = TS_VARIABLE;
		a.set.clear();
		a.vars.clear();
		a.algVars.clear();
		a.trcVars.clear();
		a.stronglyNormalized = true;
		a.sNMaxVar = "";
	}
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A BPASTriangulaSet
 **/
template <class Field, class RecursivePoly>
BPASTriangularSet<Field,RecursivePoly>& TriangularSet<Field,RecursivePoly>::operator= (BPASTriangularSet<Field,RecursivePoly>&& a) {
	if (dynamic_cast<TriangularSet<Field,RecursivePoly>*>(&a)) {
		*this = dynamic_cast<TriangularSet<Field,RecursivePoly>&&>(a);
	}
	else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to TriangularSet."));
	return *this;
}
		
/**
 * Overload operator +
 * Adds a polynomial to a triangular set and returns a new triangular set
 *
 * @param p: A sparse multivariate polynomial
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly> TriangularSet<Field,RecursivePoly>::operator+ (const RecursivePoly& p) {
	TriangularSet<Field,RecursivePoly> r(*this);
	return (r += p);
}

/**
 * Overload operator +=
 * Adds a polynomial to a triangular set
 *
 * @param p: A recursively viewed polynomial
 **/
template <class Field, class RecursivePoly>
TriangularSet<Field,RecursivePoly>& TriangularSet<Field,RecursivePoly>::operator+= (const RecursivePoly& rp) {
	if (isConstantPolynomial(rp)) {
		std::cerr << "BPAS: error, cannot add a constant to a TriangularSet" << std::endl;
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "p = " << rp << std::endl;
			// TriangularSet<Field,RecursivePoly> T(*this);
			std::cerr << "T = " << *this << std::endl;
			printVariables(vars, "T's vars");
			printVariables(trcVars, "T's trcVars");
		#endif
		exit(1);
	}
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "entering +=(p) in TriangularSet:" << std::endl;
//		std::cerr << "vars:" << std::endl;
//		printVariables(vars);
//		std::cerr << "trcVars:" << std::endl;
//		printVariables(trcVars);
//		std::cerr << "p.vars:" << std::endl;
//		printVariables(rp.variables());
//		std::cerr << "p.ringVars:" << std::endl;
//		printVariables(rp.ringVariables());
	#endif
	
	SparseMultivariateRationalPolynomial p = rp;
	p.setRingVariables(p.variables());

	Symbol lv(p.leadingVariable());
	if (std::find(algVars.begin(),algVars.end(),lv) != algVars.end()) {
		std::cerr << "BPAS: error, a polynomial with the leading variable " << p.leadingVariable() << " already exists in the TriangularSet" << std::endl;
		exit(1);
	}
	if (std::find(trcVars.begin(),trcVars.end(),lv) != trcVars.end()) {
		std::cerr << "BPAS: error, cannot add a polynomial with transcendental leading variable to the TriangularSet" << std::endl;
		exit(1);
	}
	std::vector<Symbol> vs;
	vs = orderPreservingSetIntersection(vars,p.variables());
	for (int i=0; i<vs.size(); ++i) {
		if (variableIndex(vs[i]) < variableIndex(lv) && variableIndex(lv) < vs.size()) {
			std::cerr << "BPAS: error, cannot add a polynomial with leading variable " << lv << " that also contains a higher variable, i.e., " << vs[i] << ", in the order of variables in the triangular set" << std::endl;
			std::cerr << "p = " << p << std::endl;
//			std::cerr << "vars:" << std::endl;
//			printVariables(vars);
//			std::cerr << "trcVars:" << std::endl;
//			printVariables(trcVars);
//			std::cerr << "p.variables:" << std::endl;
//			printVariables(p.variables());
//			std::cerr << "variableIndex(vs[" << i << "]) = " << variableIndex(vs[i]) << std::endl;
//			std::cerr << "variableIndex(lv) = " << variableIndex(lv) << std::endl;
//			std::cerr << "vs.size = " << vs.size() << std::endl;
			exit(1);
		}
	}
	vs = orderPreservingSetDifference(p.variables(),trcVars);
	vs = orderPreservingSetDifference(vs,vars); // list of new variables in p
	if (mode == TS_FIXED) {
		if (vs.size() != 0) {
			std::cerr << "BPAS: error, cannot add a polynomial with new variables to a fixed mode TriangularSet" << std::endl;
			printVariables(vars);
			exit(1);
		}
	}
	else if (mode == TS_VARIABLE) {
		vars.reserve(vars.size()+vs.size());
		vars.insert(vars.begin(),vs.begin(),vs.end());
		set.reserve(vars.size()+vs.size());
		std::vector<RecursivePoly> temp;
		temp.reserve(vs.size());
		for (int i=0; i<vs.size(); ++i) {
			temp.emplace_back();
		}
		set.insert(set.begin(),temp.begin(),temp.end());
	}
	int pIndex = variableIndex(lv);
	set[pIndex] = p;
	vs = orderPreservingSetIntersection(allVariables(),p.variables());
	set[pIndex].setRingVariables(vs);
	algVars.push_back(lv);
	algVars = orderPreservingSetIntersection(vars,algVars);
	updateTriangularSetStates(p);
	
	return *this;
}

/**
 * Overload comparison operator ==
 *
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
bool TriangularSet<Field,RecursivePoly>::operator== (const TriangularSet<Field,RecursivePoly>& a) const {
	if(vars.size() == a.vars.size()){
		for(int i=0; i<vars.size(); i++){
			if (set[i] != a.set[i])
				return 0;
		}
		return 1;
	}
	return 0;
}

/**
 * Overload comparison operator !=
 *
 *
 * @param a: A triangular set
 **/
template <class Field, class RecursivePoly>
bool TriangularSet<Field,RecursivePoly>::operator!= (const TriangularSet<Field,RecursivePoly>& a) const {
	return !(*this==a);
}

template <class Field, class RecursivePoly>
std::vector<Symbol> TriangularSet<Field,RecursivePoly>::allVariables() const {
	std::vector<Symbol> names;
	names = vars;
	names.insert(names.end(),trcVars.begin(),trcVars.end());
	return names;
}

/**
 * Select a polynomial given the leading variable;
 * if no such polynomial, 0 is returned
 * @param x: The leading variable name
 **/
template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::select(const Symbol& s) const {
	if (std::find(algVars.begin(),algVars.end(),s) != algVars.end()) {
		int index = variableIndex(s);
		return set[index];
	}
	else
		return RecursivePoly();
}
		
/**
 * Replace each polynomial of the regular chain with its integer primitive part
 *
 * @param p: input polynomial
 **/
template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::makePrimitive() {
	long long unsigned int rcProfilingStart;
	for (int i=0; i<set.size(); ++i) {
		if (!set[i].isZero()) {
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			set[i] = set[i].primitivePart();
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
			#endif
		}
	}
}

/**
 * Returns the ts consisting of polynomials with 
 * main variable strictly less than s
 *
 * @param s: Symbol of the main variable of specified element of the triangular set
 * @param ts: The returned triangular set
 **/
template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const {
	std::vector<RecursivePoly> polys;
	bool isNotAVariable = std::find(vars.begin(),vars.end(),s) == vars.end();
	int index = variableIndex(s);
	
	if (isNotAVariable || index == vars.size()-1) {
		if (isNotAVariable) {
			std::cerr << "BPAS: warning, Symbol s in lower(s,ts) is not a variable!" << std::endl;
			std::cerr << *this << std::endl;
			std::cerr << "s: " << s << std::endl;
			// free((int*)-1);
		}
		if (mode == TS_VARIABLE) {
			ts = TriangularSet<Field,RecursivePoly>(std::vector<Symbol>(),std::vector<Symbol>(),std::vector<Symbol>(trcVars),std::move(polys),mode,characteristic);
			return;
		}
		polys.reserve(vars.size());
		for (int i=0; i<vars.size(); ++i) {
			polys.emplace_back();
		}
		ts = TriangularSet<Field,RecursivePoly>(std::vector<Symbol>(vars),std::vector<Symbol>(),std::vector<Symbol>(trcVars),std::move(polys),mode,characteristic);
		return;
	}
	
	int size = vars.size()-index-1;
	std::vector<Symbol> avs(vars.begin()+index+1,vars.end());
//	std::vector<Symbol> avs;
	avs = orderPreservingSetIntersection(avs,algVars);
	
	polys.reserve(set.size());
	for (int i=0; i<index+1; ++i)
		polys.emplace_back();
	for (int i=index+1; i<set.size(); ++i)
		polys.emplace_back(set[i]);
	// TODO: revisit the question of whether vars below should ever be replaced by vs
	ts = TriangularSet<Field,RecursivePoly>(std::vector<Symbol>(vars),std::move(avs),std::vector<Symbol>(trcVars),std::move(polys),mode,characteristic);

}

/**
 * Returns the ts consisting of polynomials with
 * main variable strictly greater than s
 *
 * @param s: Symbol of the main variable of specified element of the triangular set
 * @param ts: The returned triangular set
 **/
template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const {
	std::vector<RecursivePoly> polys;
	bool isNotAVariable = std::find(vars.begin(),vars.end(),s) == vars.end();
	int index = variableIndex(s);
	
	if (isNotAVariable || index == 0) {
		if (isNotAVariable) {
			std::cerr << "BPAS: warning, Symbol s in upper(s,ts) is not a variable!" << std::endl;
			std::cerr << *this << std::endl;
			std::cerr << "s: " << s << std::endl;
			// free((int*)-1);
		}
		if (mode == TS_VARIABLE) {
			ts = TriangularSet<Field,RecursivePoly>(std::vector<Symbol>(),std::vector<Symbol>(),std::vector<Symbol>(trcVars),std::move(polys),mode,characteristic);
			return;
		}
		polys.reserve(vars.size());
		for (int i=0; i<vars.size(); ++i) {
			polys.emplace_back();
		}
			ts = TriangularSet<Field,RecursivePoly>(std::vector<Symbol>(vars),std::vector<Symbol>(),std::vector<Symbol>(trcVars),std::move(polys),mode,characteristic);
		return;
	}
	std::vector<Symbol> avs(vars.begin(),vars.begin()+index);
//	std::vector<Symbol> avs;
	avs = orderPreservingSetIntersection(avs,algVars);
//	std::vector<Symbol> tvs(trcVars);
	
	polys.reserve(set.size());
	for (int i=0; i<index; ++i)
		polys.emplace_back(set[i]);
	for (int i=index; i<set.size(); ++i)
		polys.emplace_back();
	ts = TriangularSet<Field,RecursivePoly>(std::vector<Symbol>(vars),std::move(avs),std::vector<Symbol>(trcVars),std::move(polys),mode,characteristic);
}

template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::cutChain(const TriangularSet<Field,RecursivePoly>& T, const Symbol& v, TriangularSet<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, TriangularSet<Field,RecursivePoly>& Tgv) const {
		T.lower(v,Tlv);
		T.upper(v,Tgv);
		Tv = T.select(v);
}

template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::cutChain(const Symbol& v, RecursivePoly& Tv, TriangularSet<Field,RecursivePoly>& Tgv) const {
		Tv = select(v);
		upper(v,Tgv);
}

template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::cutChain(const Symbol& v, TriangularSet<Field,RecursivePoly>& Tlv, RecursivePoly& Tv) const {
		lower(v,Tlv);
		Tv = select(v);
}

/**
 * Pseudo division
 * Return the pseudo-remainder, the pseudo-quotients and
 * c such that c*p = ∑(q_i T_i) + r 
 * @param p: An input polynomial
 * @param quo: The array of quotients
 * @param c: The constant multiplied to the input polynomial
 **/
template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::pseudoDivide(const RecursivePoly& p, std::vector<RecursivePoly>* q, RecursivePoly* c) const {

	
	bool qflag = false;
	bool cflag = false;
	if (q == NULL) {
		q = new std::vector<RecursivePoly>;
		qflag = true;
	}
	// if (c == NULL) {
	// 	c = new RecursivePoly;
	// 	cflag = true;
	// }
	q->reserve(this->numberOfAlgebraicVariables());
	q->clear();
	// c->one();
	// Begin naive pseudoDivide code //
//	RecursivePoly r(p);
//	RecursivePoly q;
//	std::vector<RecursivePoly> Q;
//	RecursivePoly cc;
//	RecursivePoly temp;
//	std::vector<Symbol> names(r.variables());
//	int pnvars(names.size());
// 	for (int i=0; i<vars.size(); ++i) {
// 		if (!set[i].isZero()) {
// 			r = r.pseudoDivide(set[i],&q,&cc);
// 			*c *= cc;
// 			// adjust the quotients to satisfy the constraint c*p = ∑(q_i T_i) + r //
// 			if (!quo->empty()) {
// 				for (int j=0;j<quo->size();++j) {
// 					(*quo).at(j) *= cc;
// 				}
// 			}
// 			// add new quotient and compute product of initials //
// 			quo->push_back(q);
// 		}
// 	}
//	return r;
	// End naive pseudoDivide code //
	/// /// /// /// /// /// /// /// /// 

	// triangularSetPseudoDivide code //
	RecursivePoly r;
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "ts = " << *this << std::endl;
		std::cerr << "p = " << p << std::endl;
	#endif
	
	r = p.triangularSetPseudoDivide(*this, q, c);
	
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "r = " << r << std::endl;
		if (c == NULL) {
			std::cerr << "c = null" << std::endl;
		} else {
			std::cerr << "c = " << *c << std::endl;
		}
	#endif
	/// /// /// /// /// /// /// /// /// 

	if (qflag)
		delete q;
	// if (cflag)
		// delete c;
	return RecursivePoly(r);
}

/**
 * normalForm in the sense of Groebner basis
 *
 * @param p: Input polynomial
 * 
 **/

template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::normalForm(const RecursivePoly& p, std::vector<RecursivePoly>* q) const {
	if(!stronglyNormalized){
		std::cerr<< "BPAS: error, The triangular set must be strongly normalized to compute a normal form" << std::endl;
		exit(1); 
	}
	bool qflag = false;
	if (q == NULL) {
		q = new std::vector<RecursivePoly>;
		qflag = true;
	}
	q->reserve(this->numberOfAlgebraicVariables());
	q->clear();
//	std::vector<RecursivePoly> q;
	RecursivePoly r;
	
//	RecursivePoly P(p);
//	for (int i=set.size()-1; i>-1; --i)
//		s.push_back(set[i]);
	// Update allVariables based on the variables in p
//	std::vector<Symbol> vs,vs2;
//	vs = allVariables();
//	vs2 = setDifference(P.ringVariables(),vs);
//	vs.insert(vs.begin(),vs2.begin(),vs2.end());
//	// Correct variable ordering in P
//	vs2 = orderPreservingSetIntersection(vs,P.ringVariables());
//	P.setRingVariables(vs2);

	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "calling triangularSetNormalForm..." << std::endl;
	#endif

	r = p.triangularSetNormalForm(*this,q);
		
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "calling reverse..." << std::endl;
	#endif
	
	if (!qflag) {
		std::reverse(q->begin(),q->end());
	}
	else
		delete q;
	
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "leaving normalForm..." << std::endl;
	#endif
	return RecursivePoly(r);
}

template <class Field, class RecursivePoly>
bool TriangularSet<Field,RecursivePoly>::isZeroDimensionalMathematically() const {
	if (!set.empty()) {
		std::vector<Symbol> tvs,nzvs;
		for (int i=0; i<set.size(); ++i) {
			nzvs = setUnion(nzvs,set[i].variables());
		}
		tvs = this->transcendentalVariables();
		nzvs = setDifference(nzvs,tvs);
		return (this->numberOfAlgebraicVariables() == nzvs.size());
	}
	return false;
}

template <class Field, class RecursivePoly>
bool TriangularSet<Field,RecursivePoly>::canComputeInDimensionZero(const RecursivePoly& p, bool excludeMainVariable) const {
	std::vector<Symbol> thisvs(this->mainVariables()),temp;
	temp = this->transcendentalVariables();
	thisvs.insert(thisvs.end(),temp.begin(),temp.end());
	bool variablesConsistent;
	// TODO: what if there is an inconsistency and the regular chain is fixed? -> then we don't compute in zero dimension.
	if (excludeMainVariable) {
		temp.clear();
		temp.push_back(p.leadingVariable());
		temp = setDifference(p.variables(),temp);
		variablesConsistent = ((isConstantPolynomial(p)) || isSubset(temp,thisvs));
	}
	else {
		variablesConsistent = ((isConstantPolynomial(p)) || isSubset(p.variables(),thisvs));
	}
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "completed consistency check, checking is zero dimensional mathematically..." << std::endl;
	#endif
	return (variablesConsistent && isZeroDimensionalMathematically());
}
		
/**
 * returns (c, r) such that c*r = p modulo sat(T) such that
 * c has no algebraic variables
 *
 * @param p: Input polynomial
 * 
 **/
template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::reduce(const RecursivePoly& p, RecursivePoly& c, bool takeMainPrimitivePart, bool onlyInDimZero) const {
	long long unsigned int rcProfilingStart;
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "entering reduce(p,c): " << p << std::endl;
	#endif
	c.one();
	if (isConstantPolynomial(p)) {
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "leaving reduce(p,c) (constant): " << p << std::endl;
			std::cerr << "T = " << *this << std::endl;
		#endif
		return p;
	}
	
	int index;
	Symbol v = sNMaxVar;
	RecursivePoly d,q;
	std::vector<Symbol> nonAlgVars(allVariables()),mvar;
	if (!takeMainPrimitivePart)
		nonAlgVars = setDifference(nonAlgVars,mainVariables());
//	printVariables(nonAlgVars,"nonAlgVars");
	
	if (takeMainPrimitivePart) {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		q = p.mainPrimitivePart(c);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		#ifdef TRIANGULARSET_DEBUG
		std::cerr << "calling primitive part wrt nonAlgVars... " << std::endl;
		#endif
		q = p.primitivePart(nonAlgVars,c);
//		#ifdef TRIANGULARSET_DEBUG
//		std::cerr << "calling primitive part... " << std::endl;
//		#endif
//		q = p.mainPrimitivePart(c);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
	}
	#ifdef REGULARCHAIN_PROFILING
		startTimer(&rcProfilingStart);
	#endif
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "removing integer content... " << std::endl;
	#endif
	c = c.primitivePart(); // remove integer part of content
	#ifdef REGULARCHAIN_PROFILING
		stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
	#endif
	
	// Clarify the reason for having this here.
//	if (onlyInDimZero && canComputeInDimensionZero(p)) {
//		if (stronglyNormalized) {	
//			q = normalForm(q);
//		}
//		else {
//			q = pseudoDivide(q);
//		}
//		return q;
//	}
	
	index = std::find(algVars.begin(),algVars.end(),v) - algVars.begin();
	
	if (index < 2) {
		TriangularSet<Field,RecursivePoly> Tlev;
		if (index == 1) {
			lower(algVars[0],Tlev);
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			q = Tlev.normalForm(q);
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&normalFormTime);
			#endif
		}
		else {
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			q = normalForm(q);
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&normalFormTime);
			#endif
		}
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		if (takeMainPrimitivePart)
			q = q.mainPrimitivePart(d);
		else 
			q = q.primitivePart(nonAlgVars,d);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
		c *= d;
		
		if (index == 1) {
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			q = q.pseudoDivide(select(algVars[0]));
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
				startTimer(&rcProfilingStart);
			#endif
			if (takeMainPrimitivePart)
				q = q.mainPrimitivePart(d);
			else 
				q = q.primitivePart(nonAlgVars,d);
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
			#endif
			c *= d;
			#ifdef REGULARCHAIN_PROFILING
				startTimer(&rcProfilingStart);
			#endif
			q = Tlev.normalForm(q);
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&normalFormTime);
				startTimer(&rcProfilingStart);
			#endif
			if (takeMainPrimitivePart)
				q = q.mainPrimitivePart(d);
			else
				q = q.primitivePart(nonAlgVars,d);
			#ifdef REGULARCHAIN_PROFILING
				stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
			#endif
			c *= d;
			#ifdef TRIANGULARSET_DEBUG
			std::cerr << "leaving reduce(p,c) (nF -> pD -> nF): " << p << std::endl;
			#endif
			return q;
		}
		else {
			#ifdef TRIANGULARSET_DEBUG
				std::cerr << "leaving reduce(p,c) (nF): " << p << std::endl;
			#endif
			return q;
		}
	}
	else {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
		#endif
		q = pseudoDivide(q);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
			startTimer(&rcProfilingStart);
		#endif
		if (takeMainPrimitivePart)
			q = q.mainPrimitivePart(d);
		else
			q = q.primitivePart(nonAlgVars,d);
		#ifdef REGULARCHAIN_PROFILING
			stopTimerAddElapsed(&rcProfilingStart,&primitivePartTime);
		#endif
		c *= d;
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "leaving reduce(p,c) (pD): " << p << std::endl;
		#endif
		return q;
	}
}
		
/**
 * reduce the polynomial p modulo the triangular set
 *
 * @param p: Input polynomial
 * 
 **/
template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::reduce(const RecursivePoly& p) const {
	long long unsigned int rcProfilingStart;
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "entering reduce(p): " << p << std::endl;
		std::cerr << "declaring RecursivePoly variable: " << p << std::endl;
	#endif
	RecursivePoly out;
	if (isConstantPolynomial(p)) {
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "leaving reduce(p) (constant): " << p << std::endl;
			std::cerr << "T = " << *this << std::endl;
		#endif
		return p;
	}
	else if (stronglyNormalized) {
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
			out = normalForm(p);
			stopTimerAddElapsed(&rcProfilingStart,&normalFormTime);
		#else
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "calling reduce(p) (normalForm): " << p << std::endl;
			std::cerr << "T = " << *this << std::endl;
		#endif
			out = normalForm(p);
		#endif
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "leaving reduce(p) (normalForm): " << p << std::endl;
			std::cerr << "T = " << *this << std::endl;
		#endif
		return out;
	}
	else {
		// TODO: take advantage of parts of the triangular set that have constant initials
		#ifdef REGULARCHAIN_PROFILING
			startTimer(&rcProfilingStart);
			out = pseudoDivide(p);
			stopTimerAddElapsed(&rcProfilingStart,&pseudoDivideTime);
		#else
			#ifdef TRIANGULARSET_DEBUG
				std::cerr << "calling reduce(p) (pseudoDivide): " << p << std::endl;
				std::cerr << "T = " << *this << std::endl;
			#endif
			// SMQP tmp("t");
			// if (p == tmp) {
			// 	free((int*)-1);
			// 	return tmp;
			// }
			out = pseudoDivide(p);
		#endif
		#ifdef TRIANGULARSET_DEBUG
			std::cerr << "leaving reduce(p) (pseudoDivide): " << p << std::endl;
			std::cerr << "T = " << *this << std::endl;
		#endif
		return out;
	}
}

template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::randomTriangularSetPolynomial(std::vector<Symbol> variables, int algVar, std::vector<Symbol> transcendentalVariables, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg) {
	RecursivePoly out;
	int numTrcVars(0);
	std::vector<int> indices;
	std::vector<Symbol> vs,ts;
	vs.insert(vs.begin(),variables.begin()+algVar,variables.end());
	if (transcendentalVariables.size() > 0) {
		numTrcVars = randValInRange(0,transcendentalVariables.size());
		if (numTrcVars > 0) {
			indices = randValsInRange(0,transcendentalVariables.size()-1,numTrcVars);
			sort(indices.begin(),indices.end());
			for (auto i=0; i<indices.size(); ++i) {
				ts.push_back(transcendentalVariables[indices[i]]);
			}
			vs.insert(vs.end(),ts.begin(),ts.end());
		}
	}
	out.randomPolynomial(vs.size(),nTerms,coefBound,pSparsity,includeNeg);
	out.setRingVariables(vs);
	return out;
}

template <class Field, class RecursivePoly>
RecursivePoly TriangularSet<Field,RecursivePoly>::randomTriangularSetPolynomial(std::vector<Symbol> variables, int algVar, std::vector<Symbol> transcendentalVariables, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg) {
	RecursivePoly out;
	int numTrcVars(0);
	std::vector<int> indices,newMaxDegs;
	std::vector<Symbol> vs,ts;
	vs.insert(vs.end(),variables.begin()+algVar,variables.end());
	vs.insert(vs.end(),transcendentalVariables.begin(),transcendentalVariables.end());
	newMaxDegs.insert(newMaxDegs.end(),maxDegs.begin()+algVar,maxDegs.end());
	// TODO: develop randomization for which variables other than the main variable have non-zero terms
	out.randomPolynomial(newMaxDegs,coefBound,pSparsity,includeNeg);
	out.setRingVariables(vs);
	return out;
}

/**
 * Generate a random triangular set
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
void TriangularSet<Field,RecursivePoly>::randomTriangularSet(int nVars, int nAlgVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg) {
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
	for (auto i=0; i<nAlgVars; ++i) {
		algebraicVariables.push_back(variables[indices[i]]);
	}
	for (auto i=0; i<indices.size(); ++i) {
		p = randomTriangularSetPolynomial(variables,indices[i],transcendentalVariables,nTerms,coefBound,pSparsity,includeNeg);
		set[indices[i]] = p;
	}
	algVars = algebraicVariables;
	trcVars = transcendentalVariables;
	vars = variables;
	updateTriangularSetStates();
}

/**
 * Generate a random triangular set
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
void TriangularSet<Field,RecursivePoly>::randomStronglyNormalizedTriangularSet(int nVars, int nAlgVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg) {
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
		sNMaxVar = "";
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
		for (auto i=0; i<nAlgVars; ++i) {
			algebraicVariables.push_back(variables[indices[i]]);
		}
		for (auto i=0; i<indices.size(); ++i) {
			p = randomTriangularSetPolynomial(variables,indices[i],transcendentalVariables,nTerms,coefBound,pSparsity,includeNeg);
			names.clear();
			names.push_back(p.leadingVariable());
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
 * Display the triangular set
 *
 **/
template <class Field, class RecursivePoly>
void TriangularSet<Field,RecursivePoly>::display() {
	int n = algVars.size();
	if(n==0) {
		std::cout<<" /" << std::endl;
		std::cout << " | "<<std::endl;
		std::cout<<"<  "<<std::endl;
		std::cout << " | "<<std::endl;
		std::cout<<" \\" << std::endl;
	}
	else if(n==1) {
		std::cout<<" /" << std::endl;
		std::cout << " | "<<std::endl;
		std::cout<<"<  "<< set.at(n-1) << " = 0"<<std::endl;
		std::cout << " | "<<std::endl;
		std::cout<<" \\" << std::endl;
	}
	else {
		int half = n/2;
    	int rem = n%2;
		std::cout << " /" << std::endl;\
    	int index(0);
    	for (int i = 0; i < half; ++i) {
    		while (set[index].isZero()) {
    			index++;
    		}
			std::cout << " | " << set[index] << " = 0" << std::endl;
			index++;
		}
		if(rem==0) {
			std::cout<<"< "<<std::endl;
			for (int i = half; i < n; ++i) {
				while (set[index].isZero()) {
					index++;
				}
				std::cout << " | " << set[index] << " = 0" << std::endl;
				index++;
			}
		}
		else {
			while (set[index].isZero()) {
				index++;
			}
			std::cout<< "<  " << set[index] << " = 0" <<std::endl;
			index++;
			for (int i = half+1; i < n; ++i) {
				while (set[index].isZero()) {
					index++;
				}
				std::cout << " | " << set[index] << " = 0" <<std::endl;
				index++;
			}
		}
		std::cout << " \\ " << std::endl;
	}
}

/// Possible instantiations of the TriangularSet class ///
template class TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>;
