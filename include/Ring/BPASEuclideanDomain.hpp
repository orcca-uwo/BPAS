
#ifndef _BPAS_EUCLIDEAN_DOMAIN_H_
#define _BPAS_EUCLIDEAN_DOMAIN_H_

#include "BPASGCDDomain.hpp"
#include "../Utils/TemplateHelpers.hpp"

class Integer;

/**
 * An abstract class defining the interface of a Euclidean domain.
 */
template <class Derived>
class BPASEuclideanDomain : public virtual BPASGCDDomain<Derived> {

public:

	/**
	 * Get the euclidean size of *this.
	 */
	virtual Integer euclideanSize() const = 0;

	/**
	 * Perform the eucldiean division of *this and b. Returns the
	 * remainder. If q is not NULL, then returns the quotient in q.
	 *
	 * @param b: the divisor.
	 * @param[out] q: a pointer to store the quotient in.
	 *
	 * @return the remainder.
	 */
	virtual Derived euclideanDivision(const Derived& b, Derived* q = NULL) const = 0;

	/**
	 * Perform the extended euclidean division on *this and b.
	 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
	 * @param b: the divisor.
	 * @param[out] s: the bezout coefficient of this.
	 * @param[out] t: the bezout coefficient of b.
	 * @return the remainder.
	 */
	virtual Derived extendedEuclidean(const Derived& b, Derived* s = NULL, Derived* t = NULL) const = 0;

	/**
	 * Get the quotient of *this and b.
	 *
	 * @param b: the divisor
	 * @return the quotient
	 */
	virtual Derived quotient(const Derived& b) const = 0;

	/**
	 * Get the remainder of *this and b.
	 * @param b: the divisor
	 * @return the remainder
	 */
	virtual Derived remainder(const Derived& b) const = 0;

	/**
	 * Get the remainder of *this and b;
	 * @param b: the divisor
	 * @return the remainder
	 */
	virtual Derived operator%(const Derived& b) const = 0;

	/**
	 * Assign *this to be the remainder of *this and b.
	 * @param b: the divisor
	 * @return this after assignment.
	 */
	virtual Derived& operator%=(const Derived& b) = 0;


};

#endif
