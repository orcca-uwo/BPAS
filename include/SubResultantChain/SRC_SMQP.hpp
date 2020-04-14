#ifndef _SUBRESCHAIN_SMQP_H_
#define _SUBRESCHAIN_SMQP_H_

#include "SubResultantChain/subresultantchain.hpp"
#include "Symbol/Symbol.hpp"

#include <vector>

typedef SparseMultivariateRationalPolynomial SMQP;

/**
 * A data structure class for encoded a subresultant chain between
 * two multivariate rational number polynomials.
 */
template <>
class SubResultantChain<SMQP, SMQP> {

private:

	SMQP P, Q; //deg(p, var) > deg(q, var);
	mutable std::vector<SMQP> chain;
	mutable std::vector<bool> valid;
	mutable std::vector<SMQP> chainCoefs;
	mutable std::vector<bool> validCoefs;
	Symbol var;

	void fillChain() const;


public:

	/**
	 * Default constructor: creates an empty subresultant chain.
	 *
	 **/
	SubResultantChain() {}

	/**
	 * Constructor: creates an empty subresultant chain with identified variable.
	 *
	 * @param v: variableName
	 **/
	SubResultantChain(const Symbol& v) : var(v) {}

	/**
	 * Default constructor: creates a subresultant chain from the input polynomials and their leading variable.
	 *
	 * @param a: univariate polynomial in v
	 * @param b: univariate polynomial in v
	 * @param v: variable of the input polynomials
	 **/
	SubResultantChain(const SMQP& a, const SMQP& b, const Symbol& v);

	//TODO remove this one, it's so awkward.
	SubResultantChain(const SMQP& a, const SMQP& b, const Symbol& v, std::vector<Symbol> vars);

	/**
	 * Copy constructor.
	 *
	 * @param a: A subresultant chain
	 **/
	SubResultantChain(const SMQPSubResultantChain& a) :
		P(a.P),
		Q(a.Q),
		chain(a.chain),
		valid(a.valid),
		chainCoefs(a.chainCoefs),
		validCoefs(a.validCoefs),
		var(a.var) {

	}

	/**
	 * Move constructor.
	 *
	 * @param a: An r-value subresultant chain
	 **/
	SubResultantChain(SMQPSubResultantChain&& a) :
		P(std::move(a.P)),
		Q(std::move(a.Q)),
		chain(std::move(a.chain)),
		valid(std::move(a.valid)),
		chainCoefs(std::move(a.chainCoefs)),
		validCoefs(std::move(a.validCoefs)),
		var(std::move(a.var)) {

	}

	/**
	 * Deconstructor.
	 *
	 **/
	~SubResultantChain() {}

	/**
	 * Assignment operator =.
	 *
	 * @param a: A subresultant chain
	 **/
	SMQPSubResultantChain& operator=(const SMQPSubResultantChain& a) {
		P = a.P;
		Q = a.Q;
		chain = a.chain;
		valid = a.valid;
		chainCoefs = a.chainCoefs;
		validCoefs = a.validCoefs;
		var = a.var;
		return *this;
	}

	/**
	 * Move assignment operator =.
	 *
	 * @param a: An r-value subresultant chain
	 **/
	SMQPSubResultantChain& operator=(SMQPSubResultantChain&& a) {
		P = std::move(a.P);
		Q = std::move(a.Q);
		chain = std::move(a.chain);
		valid = std::move(a.valid);
		chainCoefs = std::move(a.chainCoefs);
		validCoefs = std::move(a.validCoefs);
		var = std::move(a.var);
		return *this;
	}

	/**
	 * Identity operator ==.
	 *
	 * @param a: A subresultant chain
	 */
	bool operator==(SMQPSubResultantChain& a);

	/**
	 * Negated identity operator !=.
	 *
	 * @param a: A subresultant chain
	 **/
	inline bool operator!=(SMQPSubResultantChain& a) {
		return !(*this == a);
	}

	/**
	 * Get the size of the subresultant chain.
	 *
	 **/
	inline size_t size() const {
		return chain.size();
	}

	inline bool isEmpty() const {
		return chain.size() == 0;
	}

	/**
	 * Get the variable name.
	 *
	 **/
	inline Symbol variableName() const {
		return var;
	}

	/**
	 * Get the polynomials in the subresultant chain.
	 *
	 **/
	std::vector<SMQP> polynomials() const;

	/**
	 * Select a subresultant given the index;
	 * if no such polynomial exists, 0 is returned.
	 *
	 * @param i: index of the desired subresultant
	 */
	SMQP subResultantOfIndex(size_t i, bool lazy = 1) const;

	/**
	 * Return the first polynomial in the chain.
	 *
	 **/
	inline SMQP firstPolynomial() const {
		return P;
	}

	/**
	 * Return the second polynomial in the chain.
	 *
	 **/
	inline SMQP secondPolynomial() const {
		return Q;
	}

	/**
	 * Get the resultant of the subresultant chain.
	 * The lazy parameter indicates whether to compute only the subresultant and no other part
	 * of he subresultant chain. This is useful if, in the future, one never needs any subresultants.
	 *
	 * @param lazy, if true compute the resultant in a lazy way.
	 *
	 **/
	SMQP resultant(bool lazy = 0) const;

	/**
	 * Select the leading coefficient of a subresultant given the index;
	 * if no such polynomial exists, 0 is returned.
	 *
	 * @param i: index of the leading coefficient of the desired subresultant
	 **/
	SMQP principleSubResultantCoefficientOfIndex(int i) const;

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
    friend std::ostream& operator<< (std::ostream& ostream, const SMQPSubResultantChain& d) {
    	d.print(ostream);
        return ostream;
    }

	/**
	 * Overload stream operator <<.
	 *
	 * @param out: output stream object
	 * @param a: A subresultant chain.
	 **/
    friend std::ostream& operator<< (std::ostream& ostream, SMQPSubResultantChain&& d) {
        d.print(ostream);
        return ostream;
    }

	/**
	 * Convert subresultant chain to an expression tree.
	 *
	 **/
	ExpressionTree convertToExpressionTree() const;
};


// typedef class SubResultantChain<RationalNumber, SparseMultivariateRationalPolynomial> SMQPSubResultantChain;

#endif
