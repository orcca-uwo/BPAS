
#ifndef _BPAS_GCD_DOMAIN_H_
#define _BPAS_GCD_DOMAIN_H_

#include "BPASIntegralDomain.hpp"
#include "../Utils/TemplateHelpers.hpp"
#include "../DataStructures/Factors.hpp"


/**
 * An abstract class defining the interface of a GCD domain.
 */
template <class Derived>
class BPASGCDDomain : public virtual BPASIntegralDomain<Derived> {

public:

	/**
	 * Get GCD of *this and other.
	 *
	 * @param other: the other element to get a gcd with.
	 * @return the gcd.
	 */
	virtual Derived gcd(const Derived& other) const = 0;

	/**
	 * Compute squarefree factorization of *this.
	 *
	 * @return the square free factorization as a Factors object.
	 */
	virtual Factors<Derived> squareFree() const = 0;

};

#endif
