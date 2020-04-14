#ifndef _BPAS_POLYNOMIAL_H_
#define _BPAS_POLYNOMIAL_H_


#include "ring.h"
#include "BPASPolynomialTesters.hpp"
#include "BPASIntegralPolynomial.hpp"
#include "BPASGCDPolynomial.hpp"
#include "Utils/TemplateHelpers.hpp"
#include "Symbol/Symbol.hpp"



/**
 * An abstract class defining the interface of a polynomial over an arbitrary BPASRing.
 *
 * The actual interface of this class depends on the particular
 * specialization of the Ring template parameter
 * but will always follow that of BPASBasePolynomial at least.
 * In particular, a cascade of type checks are performed on the Ring
 * to determine if it is a GCD domain, integral domain, etc.,
 * and from this determine the interface for this polynomial.
 *
 * Derived should be the conrete class being implemented, following
 * the Curiously Recurring Template Pattern.
 *
 * Users should inherit from this class (or BPASUnivariatePolynomial,
 * or BPASMultivariatePolynomil) to get this automatic build-up of
 * the interface based on the true type of Ring,
 */
template <class Ring, class Derived>
class BPASPolynomial : public virtual BPASGCDPolynomialTester<Ring, Derived>  {};


/**
 * An abstract class defining the interface of polynomial over an arbitrary BPASRing.
 * This is class is the common interface for both univariate and multivariate
 * polynomials.
 *
 * Depending on the true type of Ring, the interface may add some additional methods.
 * @see BPASEuclideanPolynomial, BPASGCDPolynomial, BPASIntegralPolynomial.
 *
 * Since polynomials themselves form a ring, we inherit from BPASRing<Derived>,
 * meanwhile the ground ring, Ring, must also be a valid ring.
 **/
template <class Ring, class Derived>
class BPASBasePolynomial : public virtual BPASRing<Derived>,
					       private Derived_from<Ring, BPASRing<Ring>>
{
	public:

		/**
		 * Assignment operator from the ground ring.
		 *
		 * @param r: a ground ring element to assign from.
		 * @return a reference to this polynomial.
		 */
		virtual Derived& operator= (const Ring& r) = 0;

		/**
		 * Addition operator of a ground ring element.
		 *
		 * @param r: a ground ring element to add.
		 * @return the sum.
		 */
		virtual Derived operator+ (const Ring& r) const = 0;

		/**
		 * Addition assignment operator of a ground ring element.
		 *
		 * @param r: a ground ring element to add.
		 * @return a reference to this polynomial, updated after addition.
		 */
		virtual Derived& operator+= (const Ring&) = 0;

		/**
		 * Subtraction operator of a ground ring element.
		 *
		 * @param r: a ground ring element to subtract.
		 * @return the difference.
		 */
		virtual Derived operator- (const Ring&) const = 0;

		/**
		 * Subtraction assignment operator of a ground ring element.
		 *
		 * @param r: a ground ring element to subtract.
		 * @return a reference to this polynomial, updated after subtraction.
		 */
		virtual Derived& operator-= (const Ring&) = 0;

		/**
		 * Negation operation. Obtain the additive inverse of this polynomial.
		 *
		 * @return the negation of this polynomial.
		 */
		virtual Derived operator- () const = 0;

		/**
		 * Multiplication operator of a ground ring element.
		 *
		 * @param r: a ground ring element to multiply.
		 * @return the product.
		 */
		virtual Derived operator* (const Ring&) const = 0;

		/**
		 * Multiplication assignment operator of a ground ring element.
		 *
		 * @param r: a ground ring element to multiply.
		 * @return a reference to this polynomial, updated after multiplicaiton.
		 */
		virtual Derived& operator*= (const Ring&) = 0;

		/**
		 * Return the (total) degree of polynomial.
		 *
		 * @return the degree.
		 */
		virtual Integer degree() const = 0;

		/**
		 * Get the leading coefficient of this polynomial,
		 * the non-zero coefficient of the monomial with maximum degree.
		 *
		 * @return the leading coefficient as a Ring.
		 */
		virtual Ring leadingCoefficient() const = 0;

		/**
		 * Get the trailing coefficient of this polynomial,
		 * the non-zero coefficient of the monomial with minimum degree.
		 *
		 * @return the trailing coefficient as a Ring.
		 */
		virtual Ring trailingCoefficient() const = 0;

		/**
		 * Determine if the constant term of this polynomial is zero or not.
		 *
		 * @return true iff the constant term is zero.
		 */
		virtual bool isConstantTermZero() const = 0;

		/**
		 * Determine the number of non-zero terms in this polynomial.
		 *
		 * @return the number of terms.
		 */
		virtual Integer numberOfTerms() const = 0;
};




#endif