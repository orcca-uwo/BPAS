


#ifndef _BPAS_GCDPOLY_H_
#define _BPAS_GCDPOLY_H_

#include "Utils/TemplateHelpers.hpp"


/**
 * An abstract class defining the interface of a polynomial ring which
 * is also an GCD domain. E.g., a polynomial over a GCD domain.
 *
 * This class is automatically determined to be a superclass of
 * a BPASPolynomial based on a template specialization of the Ring
 * parameter.
 */
template <class Ring>
class BPASGCDPolynomial : public BPASIntegralPolynomial<Ring>,
						  public BPASGCDDomain {
	  	virtual Ring content() const = 0;
		virtual Derived primitivePart() const = 0;
};





#endif
