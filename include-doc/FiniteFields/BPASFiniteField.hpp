
#ifndef _BPAS_FINITE_FIELD_H_
#define _BPAS_FINITE_FIELD_H_

#include "../Ring/BPASField.hpp"

/**
 * An abstract class defining the interface of a prime field.
 * See SmallPrimeField, BigPrimeField, GeneralizedFermatPrimeField.
 */
class BPASFiniteField : public virtual BPASField {

public:

	virtual Derived findPrimitiveRootOfUnity(long int) const = 0;
};


#endif
