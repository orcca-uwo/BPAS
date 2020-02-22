#ifndef _SUBRESULTANTCHAIN_H_
#define _SUBRESULTANTCHAIN_H_

//#include "Polynomial/mpolynomial.h"
#include "../polynomial.h"
#include "../ring.h"
#include "../Ring/BPASIntegralDomain.hpp"
//#include "../Utils/TemplateHelpers.hpp"
//#include "../Utils/SymbolHelpers.hpp"
//#include "../RegularChain/regularchain.hpp"
//#include "../RegularChain/zerodimensionalregularchain.hpp"
#include <iostream>

template <class Domain, class UnivariateDomainPoly>
class SubResultantChain {
//: private Derived_from<Domain,BPASIntegralDomain>,
//						  private Derived_from<UnivariateDomainPoly,BPASUnivariatePolynomial<Domain,UnivariateDomainPoly>> {

	protected:
		std::vector<UnivariateDomainPoly> chain;
		Symbol var;
		
		void fillChain();
//		void fillChain (Symbol& v);


	public:
		/**
		 * Default constructor: creates an empty subresultant chain
		 *
		 * @param
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly> ();
		
		/**
		 * Constructor: creates an empty subresultant chain with identified variable
		 *
		 * @param v: variableName
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly> (const Symbol& v);
		
		/**
		 * Default constructor: creates subresultant chain from the input polynomials
		 *
		 * @param a: univariate polynomial in v
		 * @param b: univariate polynomial in v
		 * @param v: variable of the input polynomials
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly> (const UnivariateDomainPoly& a, const UnivariateDomainPoly& b, const Symbol& v);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A subresultant chain
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly> (const SubResultantChain<Domain,UnivariateDomainPoly>& a);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference subresultant chain
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly> (SubResultantChain<Domain,UnivariateDomainPoly>&& a);
		
		/**
		 * Deconstructor
		 *
		 * @param
		 **/
		~SubResultantChain<Domain,UnivariateDomainPoly>();
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A subresultant chain
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly>& operator= (const SubResultantChain<Domain,UnivariateDomainPoly>& a);
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A subresultant chain
		 **/
		SubResultantChain<Domain,UnivariateDomainPoly>& operator= (SubResultantChain<Domain,UnivariateDomainPoly>&& a);
		
		/**
		 * Identity operator ==
		 *
		 *
		 * @param a: A subresultant chain
		 **/
		bool operator== (SubResultantChain<Domain,UnivariateDomainPoly>& a);

		/**
		 * Negated identity operator !=
		 *
		 *
		 * @param a: A subresultant chain
		 **/
		bool operator!= (SubResultantChain<Domain,UnivariateDomainPoly>& a);
		
		/**
		 * Get the size of the subresultant chain
		 *
		 * @param
		 **/
		inline int size() const {
			return chain.size();
		}
		
		/**
		 * Get the variable name
		 *
		 * @param
		 **/
		inline Symbol variableName() const {
			return var;
		}
		
		/**
		 * Determine whether the subresultant chain is empty
		 *
		 **/
		inline bool isEmpty() const {
			return chain.empty();
		}
		
		/**
		 * Get the polynomials in the subresultant chain
		 *
		 **/
		inline std::vector<UnivariateDomainPoly> polynomials() const {
			return chain;
		}
		
		/**
		 * Select a subresultant given the index;
		 * if no such polynomial, 0 is returned
		 *
		 * @param i: index of the desired subresultant
		 **/
		UnivariateDomainPoly subResultantOfIndex(int i) const;
		
		/**
		 * Return the first polynomial in the chain
		 *
		 **/
		UnivariateDomainPoly firstPolynomial() const;
		
		/**
		 * Return the second polynomial in the chain
		 *
		 **/
		UnivariateDomainPoly secondPolynomial() const;
		
		/**
		 * Get the resultant of the subresultant chain
		 *
		 **/
		UnivariateDomainPoly resultant() const;
		
//		/**
//		 * Select the leading coefficient of a subresultant given the index;
//		 * if no such polynomial, 0 is returned
//		 *
//		 * @param i: index of the leading coefficient of the desired subresultant
//		 **/
//		Domain principleSubResultantCoefficientOfIndex(int i) const;

		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param a: A triangular set
		 **/
		inline friend std::ostream& operator<< (std::ostream& out, SubResultantChain<Domain,UnivariateDomainPoly>& a) {
			bool isNotFirst = 0;
			out << "[";
			for (int i = 0; i < a.chain.size(); ++i) {
				if (isNotFirst) { out << ", "; }
				out << a.chain[i];
				isNotFirst = 1;
			}
			out << "]";
			return out;
		}
		

		/**
		 * Convert subresultant chain to an expression tree
		 *
		 **/
		inline ExpressionTree convertToExpressionTree() {
			if (!chain.size()) {
				ExprTreeNode* node = new ExprTreeNode(EXPR_ARRAY, NULL, NULL, NULL);
				return ExpressionTree(node);
			}
			else {
				ExpressionTree t;
				t.fromVector<UnivariateDomainPoly>(chain);
				return(t);
			}
		}
};

#endif
