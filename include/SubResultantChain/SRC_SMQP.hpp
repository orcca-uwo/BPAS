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

	/* The following fields are mutable because what appears to be a
	   read-only operation may actually cause things to be computed
	   and subsequently cached. */
	mutable std::vector<SMQP> chain;
	mutable std::vector<bool> valid;
	mutable std::vector<SMQP> chainCoefs;
	mutable std::vector<bool> validCoefs;
	mutable std::vector<int> chainDegs;

	mutable specSRC_AAZ* lazyInfo = NULL; //a struct cacheing data for low-level algorithms.

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
		chainDegs(a.chainDegs),
		lazyInfo(a.lazyInfo),
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
		chainDegs(std::move(a.chainDegs)),
		lazyInfo(a.lazyInfo),
		var(std::move(a.var)) {
		a.lazyInfo = NULL;
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
		chainDegs = a.chainDegs;
		var = a.var;
		lazyInfo = a.lazyInfo;
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
		chainDegs = std::move(a.chainDegs);
		lazyInfo = a.lazyInfo;
		a.lazyInfo = NULL;
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
	 * Get the variable name with respect to which this subresultant chain is computed.
	 * That is, the variable for which the input polynomials are viewed
	 * as unviariate polynomials in.
	 *
	 **/
	inline Symbol variableName() const {
		return var;
	}

	/**
	 * Get all the polynomials in the subresultant chain.
	 *
	 **/
	std::vector<SMQP> polynomials() const;

	/**
	 * Select a subresultant given the index;
	 * if no such polynomial exists, 0 is returned.
	 *
	 * @param i: index of the desired subresultant
	 */
	SMQP subResultantOfIndex(size_t i) const;

	/**
	 * Return the first polynomial in the chain, the first polynomial
	 * passed to the constructor of this object.
	 *
	 **/
	inline SMQP firstPolynomial() const {
		return P;
	}

	/**
	 * Return the second polynomial in the chain, the second polynomial
	 * passed to the constructor of this object.
	 *
	 **/
	inline SMQP secondPolynomial() const {
		return Q;
	}

	/**
	 * Get the resultant of the subresultant chain.
	 * The lazy parameter indicates whether to compute only the subresultant and no other part
	 * of the subresultant chain. This is useful if, in the future, one never needs any subresultants.
	 *
	 * @param lazy, if true compute the resultant in a lazy way.
	 *
	 **/
	SMQP resultant(bool lazy = 0) const;

	/**
	 * Select the principal subresultant coefficient of the given index;
	 * the coefficient of degree i for the i'th subresultant.
	 *
	 * @param i: index of the subresultant
	 **/
	SMQP principalSubResultantCoefficientOfIndex(size_t i) const;

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
	SMQP subResultantInitialOfIndex(size_t i) const;

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
