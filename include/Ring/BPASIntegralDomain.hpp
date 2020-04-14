
#ifndef _BPAS_INTEGRAL_DOMAIN_H_
#define _BPAS_INTEGRAL_DOMAIN_H_


#include "BPASRing.hpp"
#include "../Utils/TemplateHelpers.hpp"


/**
 * An abstract class defining the interface of an integral domain.
 */
template <class Derived>
class BPASIntegralDomain : public virtual BPASRing<Derived> {

public:

	/**
	 * Exact division.
	 *
	 * @param d: the divisor.
	 * @return the equotient.
	 */
	virtual Derived operator/ (const Derived& d) const = 0;

	/**
	 * Exact division assignment.
	 *
	 * @param d: the divisor.
	 * @return a reference to this after assignment.
	 */
	virtual Derived& operator/= (const Derived& d) = 0;

};

#endif
