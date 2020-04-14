
#ifndef _BPAS_UPOLYNOMIAL_H_
#define _BPAS_UPOLYNOMIAL_H_

#include "BPASPolynomial.hpp"

#include "BPASEuclideanPolynomial.hpp"
#include "../Ring/BPASField.hpp"

/**
 * An abstract class defining the interface of a univariate polynomial over an arbitrary BPASRing.
 *
 * Depending on the specialization of the template Ring parameter, this class
 * might form a Euclidean domain, and if so, fulfills the Euclidean domain
 * interface also.
 *
 * Derived should be the conrete class being implemented, following
 * the Curiously Recurring Template Pattern.
 *
 */
template <class Ring>
class BPASUnivariatePolynomial : public virtual std::conditional<std::is_base_of<BPASField<Ring>, Ring>::value, BPASEuclideanPolynomial<Ring>, BPASPolynomial<Ring> >::type

{
	public:

		/**
		 * Differentiate this polynomial, setting itself to its derivative.
		 */
		virtual void differentiate() = 0; // p = dp/dx

		/**
		 * Differentiate this polynomial k times,
		 * setting itself to the final derivative.
		 * @param k, the number of times to differentiate.
		 */
		virtual void differentiate(int k) = 0;

		/**
		 * Obtain the derivative of this polynomial.
		 *
		 * @return the derivative.
		 */
		virtual Derived derivative() const = 0;

		/**
		 * Obtain the kth derivative of this polynomial.
		 *
		 * @param k, the number of times to differentiate.
		 * @return the kth derivative.
		 */
		virtual Derived derivative(int) const = 0;

		/**
		 * Evaluate this polynomial by substituting the input ring
		 * element r for the indeterminate.
		 *
		 * @param r the input ring element
		 * @return the evaluation of this polynomial at r.
		 */
		virtual Ring evaluate(const Ring& r) const = 0;

		/**
		 * Divide this polynomial by the monic polynomial d,
		 * setting this polynomial to be the remainder.
		 *
		 * @param d, the monic divison.
		 * @return the quotient.
		 */
		virtual Derived monicDivide(const Derived& d) = 0;

		/**
		 * Divide this polynomial by the monic polynomial d.
		 * The remainder is returned and the quotient returned
		 * via pointer q.
		 *
		 * @param d, the monic divisor.
		 * @param[out] q, a pointer to the quotient to return.
		 * @return the remainder.
		 */
		virtual Derived monicDivide(const Derived& d, Derived* q = NULL) const = 0;

		/**
		 * Perform a lazy pseudo-divison, also known as sparse pseudo-division.
		 * Specifically, determine q and r, as in division, but satisfying
		 * h^e * this = q*d + r, where h is the leading coefficient of d
		 * and e is a positive integer equal to the number of division steps performed.
	 	 * Return the quotient and this becomes the remainder.
		 *
		 * @param d, the divisor.
		 * @param[out] h1: The leading coefficient of b to the power e
	 	 * @param[out] h2: h to the power deg(a) - deg(b) + 1 - e
		 * @return the quotient.
		 */
		virtual Derived lazyPseudoDivide(const Derived& d, Ring* h1, Ring* h2) = 0;


		/**
		 * Perform a lazy pseudo-divison, also known as sparse pseudo-division.
		 * Specifically, determine q and r, as in division, but satisfying
		 * h^e * this = q*d + r, where h is the leading coefficient of d
		 * and e is a positive integer equal to the number of division steps performed.
	 	 * Return the quotient and this becomes the remainder.
		 *
		 * @param d, the divisor.
		 * @param[out] q:  a pointer to the quotient to output.
		 * @param[out] h1: The leading coefficient of b to the power e
	 	 * @param[out] h2: h to the power deg(this) - deg(d) + 1 - e
		 * @return the remainder.
		 */
		virtual Derived lazyPseudoDivide(const Derived& d, Derived* q, Ring* h1, Ring* h2) const = 0;

		/**
		 * Perform a pseudo-divison, specifically,
		 * determine q and r, as in division, but satisfying
		 * h^e * this = q*d + r, where h is the leading coefficient of d
		 * and e is a positive integer equal to deg(this) - deg(d) + 1
		 * Return the quotient and this becomes the remainder.
		 *
		 * @param d, the divisor.
	 	 * @param[out] he: h to the power deg(this) - deg(d) + 1
		 * @return the quotient.
		 */
		virtual Derived pseudoDivide(const Derived& d, Ring* he) = 0;

		/**
		 * Perform a pseudo-divison, specifically,
		 * determine q and r, as in division, but satisfying
		 * h^e * this = q*d + r, where h is the leading coefficient of d
		 * and e is a positive integer equal to deg(this) - deg(d) + 1
		 * Return the quotient and this becomes the remainder.
		 *
		 * @param d, the divisor.
		 * @param[out] q:  a pointer to the quotient to output.
	 	 * @param[out] he: h to the power deg(this) - deg(d) + 1
		 * @return the remainder.
		 */
		virtual Derived pseudoDivide(const Derived& d, Derived* q, Ring* he) const = 0;

		/**
		 * Get the coefficient of the monomial with degree d.
		 *
		 * @param d, the degree of the monomial.
		 * @return the coefficient of that monomial.
		 */
		virtual Ring coefficient(int d) const = 0;

		/**
		 * Set the coefficient of the monomial with degree d to
		 * be the Ring element r.
		 *
		 * @param d, the degree of the monomial.
		 * @param r, the ring element
		 */
		virtual void setCoefficient(int d, const Ring& r) = 0;

		/**
		 * Set the indeterminate of this polynomial to be the input symbol.
		 *
		 * @param sym, the new symbol for the indeterminate.
		 */
		virtual void setVariableName (const Symbol& sym) = 0;

		/**
		 * Get the indeterminate of this polynomial.
		 *
		 * @return this polynomial's indeterminate.
		 */
		virtual Symbol variable() const = 0;

		/**
		 * Shift this polynomial left i times, that is,
		 * multiply by x^i, where x is this polynomial's indeterminate.
		 *
		 * @return the shifted polynomial.
		 */
		virtual Derived operator<< (int i) const = 0; // q = p * (x^i);

		/**
		 * Shift this polynomial left i times, that is,
		 * multiply by x^i, where x is this polynomial's indeterminate,
		 * and set this polynomial to the result.
		 *
		 * @return a reference to this polynomial after shifting.
		 */
		virtual Derived& operator<<= (int i) = 0; // p = p *(x^i)

		/**
		 * Shift this polynomial right i times, that is,
		 * divide by x^i, where x is this polynomial's indeterminate,
		 *
		 * @return the shifted polynomial.
		 */
		virtual Derived operator>> (int) const = 0; // q = p / (x^i);

		/**
		 * Shift this polynomial right i times, that is,
		 * divide by x^i, where x is this polynomial's indeterminate,
		 * and set this polynomial to the result.
		 *
		 * @return a reference to this polynomial after shifting.
		 */
		virtual Derived& operator>>= (int) = 0;
};




#endif