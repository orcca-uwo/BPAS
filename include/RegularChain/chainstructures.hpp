#ifndef _CHAINSTRUCTURES_H_
#define _CHAINSTRUCTURES_H_

#include <string>
#include "../ExpressionTree/ExprTreeNode.hpp"
#include "../ExpressionTree/ExpressionTree.hpp"

/**
 * A templated data structure class designed to hold a pair consisting of a
 * recursively viewed polynomial and a regular chain.
 **/
template <class RecursivePoly, class RegularChainType>
class PolyChainPair {
	public:
		RecursivePoly poly;
		RegularChainType chain;
		
	   /**
	    * Default constructor.
	    *
	    * @param
	    **/
		PolyChainPair () {
			poly = RecursivePoly();
			chain = RegularChainType();
		}
		
	   /**
	    * A constructor that creates a PolyChainPair object consisting of copies
	    * of the supplied arguments.
	    *
	    * @param p: a recursively viewed polynomial
	    * @param rc: a regular chain
	    **/
		PolyChainPair (const RecursivePoly& p, const RegularChainType& rc) {
			poly = RecursivePoly(p);
			chain = RegularChainType(rc);
		}
		
		
	   /**
	    * Copy constructor.
	    *
	    * @param a: a PolyChainPair
	    **/
		PolyChainPair (const PolyChainPair<RecursivePoly,RegularChainType>& a) {
			poly = a.poly;
			chain = a.chain;
		}
		
		
	   /**
	    * Move constructor.
	    *
	    * @param a: a PolyChainPair
	    **/
		PolyChainPair (const PolyChainPair<RecursivePoly,RegularChainType>&& a) {
			poly = std::move(a.poly);
			chain = std::move(a.chain);
		}
		
		
	   /**
	    * Default destructor.
	    *
	    * @param
	    **/
		~PolyChainPair () {}
		
		
	   /**
	    * Assignment operator.
	    *
	    * @param a: a PolyChainPair
	    **/
		inline PolyChainPair& operator= (const PolyChainPair<RecursivePoly,RegularChainType>& a) {
			poly = a.poly;
			chain = a.chain;
			return *this;
		}
		
		
	   /**
	    * Move assignment operator =.
	    *
	    * @param a: a PolyChainPair
	    **/
		inline PolyChainPair& operator= (const PolyChainPair<RecursivePoly,RegularChainType>&& a) {
			poly = std::move(a.poly);
			chain = std::move(a.chain);
			return *this;
		}
		
//		friend bool operator==(const PolyChainPair& s1,const PolyChainPair& s2);
//		friend bool operator!=(const PolyChainPair& s1,const PolyChainPair& s2) {
//			return !(s1==s2);
//		}
		// need to extend ExpressionTree to cover PolyChainPairs
		/*ExpressionTree toExpressionTree() const {
			ExprTreeNode etn(*this);
			ExpressionTree et(&etn);
			return et;
		}*/
		
	   /**
	    * Convert the PolyChainPair to a string.
	    *
	    * @param
	    **/
		std::string toString() const {
			std::stringstream ss;
			std::string out;
			ss << "<" << poly << "," << chain << ">";
			out = ss.str();
			return out;
		}
		
	   /**
	    * Overload operator<<.
	    *
	    * @param
	    **/
		inline friend std::ostream& operator<< (std::ostream &out, const PolyChainPair<RecursivePoly,RegularChainType>& a) {
			out << a.toString();
			return out;
		}

};


/**
 * A templated data structure class designed to hold a pair consisting of a
 * boolean value and a regular chain.
 **/
template <class RegularChainType>
class BoolChainPair {
	public:
		bool isTrue;
		RegularChainType chain;
		
	   /**
	    * Default constructor.
	    *
	    * @param
	    **/
		BoolChainPair () {
			isTrue = false;
			chain = RegularChainType();
		}	
		
	   /**
	    * A constructor that creates a BoolChainPair object consisting of copies
	    * of the supplied arguments.
	    *
	    * @param p: a boolean value
	    * @param rc: a regular chain
	    **/
		BoolChainPair (bool b, const RegularChainType& rc) {
			isTrue = b;
			chain = RegularChainType(rc);
		}
		
	   /**
	    * Copy constructor.
	    *
	    * @param a: a BoolChainPair
	    **/
		BoolChainPair (const BoolChainPair<RegularChainType>& a) {
			isTrue = a.isTrue;
			chain = a.chain;
		}
		
	   /**
	    * Move constructor.
	    *
	    * @param a: a BoolChainPair
	    **/
		BoolChainPair (const BoolChainPair<RegularChainType>&& a) {
			isTrue = a.isTrue;
			chain = std::move(a.chain);
		}
		
	   /**
	    * Default destructor.
	    *
	    * @param
	    **/
		~BoolChainPair () {}
		
	   /**
	    * Assignment operator.
	    *
	    * @param
	    **/
		inline BoolChainPair& operator= (const BoolChainPair<RegularChainType>& a) {
			isTrue = a.isTrue;
			chain = a.chain;
			return *this;
		}
		
	   /**
	    * Move assignment operator.
	    *
	    * @param
	    **/
		inline BoolChainPair& operator= (const BoolChainPair<RegularChainType>&& a) {
			isTrue = a.isTrue;
			chain = std::move(a.chain);
			return *this;
		}
		
//		friend bool operator==(const BoolChainPair& s1,const BoolChainPair& s2);
//		friend bool operator!=(const BoolChainPair& s1,const BoolChainPair& s2) {
//			return !(s1==s2);
//		}
		// need to extend ExpressionTree to cover BoolChainPairs
		/*ExpressionTree toExpressionTree() const {
			ExprTreeNode etn(*this);
			ExpressionTree et(&etn);
			return et;
		}*/
		
	   /**
	    * Covert BoolChainPair object to a string.
	    *
	    * @param
	    **/
		std::string toString() const {
			std::stringstream ss;
			std::string out;
			ss << "<";
			if (isTrue)
				ss << "true";
			else
				ss << "false";
			ss << "," << chain << ">";
			out = ss.str();
			return out;
		}
		
	   /**
	    * Overload operator<<.
	    *
	    * @param
	    **/
		inline friend std::ostream& operator<< (std::ostream &out, const BoolChainPair<RegularChainType>& a) {
			out << a.toString();
			return out;
		}

};

#endif
