
#ifndef _BPAS_FIELD_H_
#define _BPAS_FIELD_H_

#include "BPASEuclideanDomain.hpp"

/**
 * An abstract class defining the interface of a field.
 */
template <class Derived>
class BPASField : public virtual BPASEuclideanDomain<Derived> {

public:

	/**
	 * Get the inverse of *this.
	 *
	 * @return the inverse
	 */
	virtual Derived inverse() const = 0;
};

#endif
