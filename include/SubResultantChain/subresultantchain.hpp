#ifndef _SUBRESULTANTCHAIN_H_
#define _SUBRESULTANTCHAIN_H_

#include "../Polynomial/BPASRecursivePolynomial.hpp"
#include "../Polynomial/BPASUnivarPolynomial.hpp"
#include "../RationalNumberPolynomial/mrpolynomial.h"

/**
 * A data structure class for encoded a subresultant chain between
 * two polynomials over a ring. The polynomials must be
 * either univariate or a recursively-viewed polynomial.
 * This is enforced by the otherwise unnecessary inheriance.
 *
 * Ring should be the ground ring of the polynomial ring defined
 * by Poly.
 */

template <class Ring, class Poly>
class SubResultantChain : private std::conditional<std::is_base_of<BPASUnivariatePolynomial<Ring,Poly>, Poly>::value, Derived_from<Poly, BPASUnivariatePolynomial<Ring,Poly>>, Derived_from<Poly, BPASRecursivelyViewedPolynomial<Ring,Poly>>> {

public:

	/**
	 * Default constructor: creates an empty subresultant chain.
	 *
	 **/
	SubResultantChain ();

	/**
	 * Constructor: creates an empty subresultant chain with identified variable.
	 *
	 * @param v: variableName
	 **/
	SubResultantChain (const Symbol& v);

	/**
	 * Default constructor: creates a subresultant chain from the input polynomials and their leading variable.
	 *
	 * @param a: univariate polynomial in v
	 * @param b: univariate polynomial in v
	 * @param v: variable of the input polynomials
	 **/
	SubResultantChain (const Poly& a, const Poly& b, const Symbol& v);

	//TODO remove this one, it's so awkward.
	SubResultantChain (const Poly& a, const Poly& b, const Symbol& v, std::vector<Symbol> vars);

	/**
	 * Get the size of the subresultant chain.
	 *
	 **/
	size_t size() const;

	/**
	 * Determine if the chain is empty
	 */
	bool isEmpty() const;


	/**
	 * Get the variable name.
	 *
	 **/
	Symbol variableName() const;

	/**
	 * Get the polynomials in the subresultant chain.
	 *
	 **/
	std::vector<Poly> polynomials() const;

	/**
	 * Select a subresultant given the index;
	 * if no such polynomial exists, 0 is returned.
	 *
	 * @param i: index of the desired subresultant
	 **/
	Poly subResultantOfIndex(size_t i, bool lazy = 0) const;

	/**
	 * Return the first polynomial in the chain.
	 * That is, one of the two inputs of higher degree.
	 *
	 **/
	Poly firstPolynomial() const;

	/**
	 * Return the second polynomial in the chain.
	 * That is, one of the two inputs of lower degree.
	 *
	 **/
	Poly secondPolynomial() const;

	/**
	 * Get the resultant of the subresultant chain.
	 * The lazy parameter indicates whether to compute only the subresultant and no other part
	 * of he subresultant chain. This is useful if, in the future, one never needs any subresultants.
	 *
	 * @param lazy: if true compute the resultant in a lazy way.
	 *
	 **/
	Ring resultant(bool lazy = 0) const;

	/**
	 * Select the principal subresultant coefficient of the given index;
	 * the coefficient of degree i for the i'th subresultant.
	 *
	 * @param i: index of the subresultant
	 **/
	Ring principalSubResultantCoefficientOfIndex(int i) const;

	/**
	 * Select the initial of the subresultant of the given index;
	 * the leading coefficient of the subresultant viewed as a
	 * univariate polynomial in v.
	 *
	 * @note: this is different than the principal subresultant coefficient;
	 *        this method returns the actual leading coefficient,
	 *        thus returning 0 iff the polynomial is identically 0.
	 *
	 * @param i: index of the subresultant
	 */
	Ring subResultantInitialOfIndex(size_t i) const;


	/**
	 * Given an output stream, print the subresultant chain.
	 *
	 * @param ostream: the output stream.
	 */
	void print(std::ostream& ostream) const;

	/**
	 * Overload stream operator <<.
	 *
	 * @param out: output stream object
	 * @param a: A subresultant chain.
	 **/
    friend std::ostream& operator<< (std::ostream& ostream, const SubResultantChain<Ring,Poly>& d) {
    	d.print(ostream);
        return ostream;
    }

	/**
	 * Overload stream operator <<.
	 *
	 * @param out: output stream object
	 * @param a: A subresultant chain.
	 **/
    friend std::ostream& operator<< (std::ostream& ostream, SubResultantChain<Ring,Poly>&& d) {
        d.print(ostream);
        return ostream;
    }

	/**
	 * Convert subresultant chain to an expression tree.
	 *
	 **/
	ExpressionTree convertToExpressionTree() const;
};


/* Template declares so that includers of SubResultantChain can get specializations for free */

template<> class SubResultantChain<SparseMultivariateRationalPolynomial,SparseMultivariateRationalPolynomial>;
typedef SubResultantChain<SparseMultivariateRationalPolynomial, SparseMultivariateRationalPolynomial> SMQPSubResultantChain;

#include "SubResultantChain/SRC_SMQP.hpp"

#endif
