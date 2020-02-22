
#ifndef _BPAS_FIELD_H_
#define _BPAS_FIELD_H_

#include "BPASEuclideanDomain.hpp"

/**
 * An abstract class defining the interface of a field.
 */
class BPASField : public virtual BPASEuclideanDomain {

public:

	/**
	 * Get the inverse of *this.
	 */
	virtual Derived inverse() const = 0;
};

#endif
