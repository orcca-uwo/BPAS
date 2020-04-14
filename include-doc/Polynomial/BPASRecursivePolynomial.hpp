
#ifndef _BPAS_RECURSIVE_POLY_H_
#define _BPAS_RECURSIVE_POLY_H_

#include "BPASMultivarPolynomial.hpp"

/**
 * An abstract class defining the interface of a multivariate polynomial that can be viewed recursively.
 * That is, it can be viewed as a univariate polynomial with multivariate polynomial coefficients.
 */
template <class Ring>
class BPASRecursivelyViewedPolynomial : public virtual BPASMultivariatePolynomial<Ring>
{
	public:

		/**
		 * Get the initial, or leading coefficient, with respect to the main variable of
		 * this polynomial.
		 *
		 * @return the initial of this polynomial.
		 */
	    virtual Derived initial() const = 0;

	    /**
	     * Get the main variable of this polynomial, that is, the one
	     * which is most significant under its variable ordering.
	     *
	     * @return the main variable of the polynomial.
	     */
	    virtual Symbol mainVariable() const = 0;

		/**
		 * Get the main degree of this polynomial, that is, the
		 * partial degree of the main variable.
		 *
		 * @return the main degree of the polynomial.
		 */
		virtual int mainDegree() const = 0;

		/**
		 * Get the rank of this polynomial, that is, the
		 * main variable of this polynomial raised to the main degree.
		 * This this returns a monomial.
		 *
		 * @return the rank of this polynomial.
		 */
		virtual Derived rank() const = 0;

		/**
		 * Get the tail of this polynomial, that is,
		 * this polynomial after removing its leading term.
		 * It is equivalent to the reductum of this polynomial, viewing
		 * it as a univariate polynomial in its main variable.
 		 *
 		 * @return the tail of this polynomial.
		 */
 		virtual Derived tail() const = 0;

 		/**
 		 * Get the head of this polynomial, that is,
 		 * the leading term of this polynomial.
		 * It is equivalent to the initial of this polynomial
		 * multiplied by its rank.
		 *
 		 * @return the head of this polynomial.
 		 */
 		virtual Derived head() const = 0;

 		/**
 		 * Get the separant of this polynomial, that is,
 		 * its derivative with respect to the main variable.
 		 *
 		 * @return the separant of this polynomial.
 		 */
 		virtual Derived separant() const = 0;

 		// virtual SparseUnivariatePolynomial<Derived> convertToSUP() const = 0;

};

#endif
