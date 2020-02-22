#include "../../include/SubResultantChain/subresultantchain.hpp"
#include "../../include/RingPolynomial/upolynomial.h"
#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include "../../include/ring.h"
#include <iostream>

/// Protected Functions ///

///**
// * Fill a sparse subresultant chain with zeros
// *
// * @param
// **/
//template <class Domain, class UnivariateDomainPoly>
//void SubResultantChain<Domain,UnivariateDomainPoly>::fillChain () {
//	UnivariateDomainPoly zero;
//	zero.zero();
//	int fullSize(chain[chain.size()-2].degree().get_ui()+2);
//	int delta;
////	std::cerr << "chain.size() = " << chain.size() << std::endl;
////	std::cerr << "fullSize = " << fullSize << std::endl;
//	if (chain.size() < fullSize) {
//		chain.reserve(fullSize);
//		for (int i=chain.size()-2; i>0; --i) {
//			if (chain[i].degree() != chain[i-1].degree()+1) {
//				delta = chain[i].degree().get_ui() - chain[i-1].degree().get_ui();
//				if (i > 1) {
//					i = i-1;
//					for (int j=0; j<delta-2; ++j)
//						chain.insert(chain.begin()+i,zero);
//				}
//				else {
//					for (int j=0; j<delta-1; ++j)
//						chain.insert(chain.begin()+i,zero);
//				}
//			}
//		}
//		if (chain[0].degree() != 0) {
//				for (int j=0; j<chain[0].degree(); ++j)
//					chain.insert(chain.begin(),zero);
//		}
//	}
////	std::cerr << "chain.size() = " << chain.size() << std::endl;
//}

///**
// * Fill a sparse subresultant chain with zeros
// *
// * @param
// **/
//template <class Domain, class UnivariateDomainPoly>
//void SubResultantChain<Domain,UnivariateDomainPoly>::fillChain (Symbol& v) {
//	UnivariateDomainPoly zero;
//	zero.zero();
//	int fullSize(chain[chain.size()-2].degree(v).get_ui()+2);
//	int delta;
////	std::cerr << "chain.size() = " << chain.size() << std::endl;
////	std::cerr << "fullSize = " << fullSize << std::endl;
//	if (chain.size() < fullSize) {
//		chain.reserve(fullSize);
//		for (int i=chain.size()-2; i>0; --i) {
//			if (chain[i].degree(v) != chain[i-1].degree(v)+1) {
//				delta = chain[i].degree(v).get_ui() - chain[i-1].degree(v).get_ui();
//				if (i > 1) {
//					i = i-1;
//					for (int j=0; j<delta-2; ++j)
//						chain.insert(chain.begin()+i,zero);
//				}
//				else {
//					for (int j=0; j<delta-1; ++j)
//						chain.insert(chain.begin()+i,zero);
//				}
//			}
//		}
//		if (chain[0].degree(v) != 0) {
//				for (int j=0; j<chain[0].degree(v); ++j)
//					chain.insert(chain.begin(),zero);
//		}
//	}
////	std::cerr << "chain.size() = " << chain.size() << std::endl;
//}

/// Constructors ///

/**
 * Default constructor: creates an empty subresultant chain of variable size
 * with empty list of transcendentals
 *
 * @param
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::SubResultantChain () : var("x") {}
		
/**
 * Constructor: creates an empty subresultant chain with identified variable
 *
 * @param v: variableName
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::SubResultantChain (const Symbol& v) : var(v) {}
		
/**
 * Default constructor: creates subresultant chain from the input polynomials
 *
 * @param a: univariate polynomial in v
 * @param b: univariate polynomial in v
 * @param v: variable of the input polynomials
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::SubResultantChain (const UnivariateDomainPoly& a, const UnivariateDomainPoly& b, const Symbol& v) : var(v) {
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "[SRC] calling a.subresultantChain(b,1):" << std::endl;
		std::cerr << "[SRC] a = " << a << std::endl;
		std::cerr << "[SRC] b = " << b << std::endl;
	#endif
	chain = a.subresultantChain(b,1);

	std::cerr << "Computing Chain Between: " << std::endl << "P: " << a << std::endl << "Q: " << b << std::endl;

	std::cerr << "The resultant: " << chain[0] << std::endl;
	printVariables(chain[0].ringVariables(), "the resultant vars");
	fprintf(stderr, "\n\n");

	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "[SRC] leaving default constructor." << std::endl;
	#endif
}
		
/**
 * Ring specification constructor: creates subresultant chain from the input polynomials in a minimal ordered subring of R
 *
 * @param a: univariate polynomial in v
 * @param b: univariate polynomial in v
 * @param v: variable of the input polynomials
 * @param R: ordered polynomial ring containing a and b
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::SubResultantChain (const UnivariateDomainPoly& a, const UnivariateDomainPoly& b, const Symbol& v, const std::vector<Symbol>& R) : var(v) {
	#ifdef TRIANGULARSET_DEBUG
//		printVariables(R,"R");
	#endif
	std::vector<Symbol> R2;
	R2 = setUnion(a.variables(),b.variables());
	#ifdef TRIANGULARSET_DEBUG
//		printVariables(R2,"AUB");
	#endif
	R2 = orderPreservingSetDifference(R,R2);
	#ifdef TRIANGULARSET_DEBUG
//		printVariables(R2,"R2-R");
	#endif
	R2 = orderPreservingSetDifference(R,R2);
	#ifdef TRIANGULARSET_DEBUG
//		printVariables(R2,"R2-R");
		std::cerr << "[SRC] calling a.subresultantChain(b,1):" << std::endl;
//		std::cerr << "[SRC] a = " << a << std::endl;
//		std::cerr << "[SRC] b = " << b << std::endl;
	#endif
	chain = a.subresultantChain(b,1);
	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "[SRC] returning to constructor." << std::endl;
	#endif
	for (int i=0; i<chain.size(); ++i)
		chain[i].setRingVariables(R);
	chain[0].setRingVariables(chain[0].variables());
	chain[chain.size()-1].setRingVariables(chain[chain.size()-1].variables());
	chain[chain.size()-2].setRingVariables(chain[chain.size()-2].variables());
	

	std::cerr << "Computing Chain Between: " << std::endl << "P: " << a << std::endl << "Q: " << b << std::endl;


	std::cerr << "The resultant: " << chain[0] << std::endl;
	printVariables(chain[0].ringVariables(), "the resultant vars");
	fprintf(stderr, "\n\n");

	#ifdef TRIANGULARSET_DEBUG
		std::cerr << "[SRC] leaving ring specification constructor." << std::endl;
	#endif
	
	
//	UnivariateDomainPoly A(a),B(b);
//	A.setRingVariables(R);
//	A.setRingVariables(A.variables());
//	B.setRingVariables(R);
//	B.setRingVariables(B.variables());
//	std::vector<Symbol> R2;
//	R2 = setUnion(A.variables(),B.variables());
//	R2 = setDifference(R,R2);
//	R2 = setDifference(R,R2);
//	std::cerr << "[SRC] calling A.subresultantChain(B,1):" << std::endl;
//	std::cerr << "[SRC] A = " << A << std::endl;
//	std::cerr << "[SRC] B = " << B << std::endl;
//	chain = A.subresultantChain(B,1);
//	std::cerr << "[SRC] returning to constructor." << std::endl;
//	for (int i=0; i<chain.size(); ++i)
//		chain[i].setRingVariables(R2);
//	chain[0].setRingVariables(chain[0].variables());
//	chain[chain.size()-1].setRingVariables(chain[chain.size()-1].variables());
//	chain[chain.size()-2].setRingVariables(chain[chain.size()-2].variables());
//	std::cerr << "[SRC] leaving ring specification constructor." << std::endl;
}

/**
 * Copy constructor
 *
 * @param a: A subresultant chain
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::SubResultantChain (const SubResultantChain<Domain,UnivariateDomainPoly>& a) : var(a.var), chain(a.chain) {}

/**
 * Move constructor
 *
 * @param a: An r-value reference subresultant chain
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::SubResultantChain (SubResultantChain<Domain,UnivariateDomainPoly>&& a) : var("x") {
	var = a.var;
	chain = std::move(a.chain);
	a.var = "x";
	a.chain.clear();
}

/**
 * Deconstructor
 *
 * @param
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>::~SubResultantChain() {}

/**
 * Assignment operator =
 *
 * @param a: A subresultant chain
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>& SubResultantChain<Domain,UnivariateDomainPoly>::operator= (const SubResultantChain<Domain,UnivariateDomainPoly>& a) {
	if (this != &a) {
		var = a.var;
		chain = a.chain;
	}
	return *this;
}

/**
 * Move assignment operator =
 *
 * @param a: A subresultant chain
 **/
template <class Domain, class UnivariateDomainPoly>
SubResultantChain<Domain,UnivariateDomainPoly>& SubResultantChain<Domain,UnivariateDomainPoly>::operator= (SubResultantChain<Domain,UnivariateDomainPoly>&& a) {
	if (this != &a) {
		var = a.var;
		chain = std::move(a.chain);
		a.var = "x";
		a.chain.clear();
	}
	return *this;
}

/**
 * Identity operator ==
 *
 *
 * @param a: A subresultant chain
 **/
template <class Domain, class UnivariateDomainPoly>
bool SubResultantChain<Domain,UnivariateDomainPoly>::operator== (SubResultantChain<Domain,UnivariateDomainPoly>& a) {
	if(chain.size() == a.chain.size()){
		for(int i=0; i<chain.size(); i++){
			if (chain[i] != a.chain[i])
				return 0;
		}
		return 1;
	}
	return 0;
}

/**
 * Negated identity operator !=
 *
 *
 * @param a: A subresultant chain
 **/
template <class Domain, class UnivariateDomainPoly>
bool SubResultantChain<Domain,UnivariateDomainPoly>::operator!= (SubResultantChain<Domain,UnivariateDomainPoly>& a) {
	return !(*this==a);
}

/**
 * Select a subresultant given the index;
 * if no such polynomial, 0 is returned
 *
 * @param i: index of the desired subresultant
 **/
template <class Domain, class UnivariateDomainPoly>
UnivariateDomainPoly SubResultantChain<Domain,UnivariateDomainPoly>::subResultantOfIndex(int i) const {
	if (i <= chain.size() && !(chain.empty()))
		return chain[i];
	else {
		std::cerr << "BPAS: error, requested subresultant does not exist!" << std::endl;
		exit(1);
	}
}
		
/**
 * Return the first polynomial in the chain
 *
 **/
template <class Domain, class UnivariateDomainPoly>
UnivariateDomainPoly SubResultantChain<Domain,UnivariateDomainPoly>::firstPolynomial() const {
	if (!chain.empty())
		return chain.back();

	UnivariateDomainPoly zero;
	zero.zero();
	return zero;
}

/**
 * Return the second polynomial in the chain
 *
 **/
template <class Domain, class UnivariateDomainPoly>
UnivariateDomainPoly SubResultantChain<Domain,UnivariateDomainPoly>::secondPolynomial() const {
	if (!chain.empty())
		return chain[chain.size()-2];

	UnivariateDomainPoly zero;
        zero.zero();
	return zero;

}
/**
 * Get the resultant of the subresultant chain
 *
 **/
template <class Domain, class UnivariateDomainPoly>
UnivariateDomainPoly SubResultantChain<Domain,UnivariateDomainPoly>::resultant(bool lazy) const {
	if (!chain.empty())
		return chain[0];

	UnivariateDomainPoly zero;
        zero.zero();
	return zero;
}
		
///**
// * Select the leading coefficient of a subresultant given the index;
// * if no such polynomial, 0 is returned
// *
// * @param i: index of the leading coefficient of the desired subresultant
// **/
//template <class Domain, class UnivariateDomainPoly>
//Domain SubResultantChain<Domain,UnivariateDomainPoly>::principleSubResultantCoefficientOfIndex(int i) const {
//	if (i <= chain.size() && !(chain.empty()))
//		return chain[i].leadingCoefficientInVariable(var);
//	else {
//		std::cerr << "BPAS: error, subresultant for requested coefficient does not exist!" << std::endl;
//		exit(1);
//	}
//}

/// Possible instantiations of the SubResultantChain class ///
//template class SubResultantChain<RationalNumber,SparseUnivariatePolynomial<RationalNumber>>;
//template class SubResultantChain<RationalNumber,SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>>;
template class SubResultantChain<SparseMultivariateRationalPolynomial,SparseMultivariateRationalPolynomial>;
//template class SubResultantChain<RationalNumber,SparseMultivariateRationalPolynomial>;
