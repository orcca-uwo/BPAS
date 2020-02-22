
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
	 */
	virtual Derived gcd(const Derived& other) const = 0;
	
	/**
	 * Compute squarefree factorization of *this
	 */
	virtual Factors<Derived> squareFree() const = 0;

};

#endif
