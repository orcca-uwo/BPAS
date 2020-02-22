
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
	 */
	virtual Derived operator/ (const Derived&) const = 0;

	/**
	 * Exact division assignment.
	 */
	virtual Derived& operator/= (const Derived&) = 0;

};

#endif
