
/**
 * Data Structure for multivariate modulo polynomial
 * stored in a distributed dense case
 **/

#include "RingPolynomial/dmpolynomial.h"

template <class Field>
mpz_class DistributedDenseMultivariateModularPolynomial<Field>::characteristic(5);

//empty class to test compiling.
//TODO properly separate implementation out of dmpolynomial.h

DistributedDenseMultivariateModularPolynomial<RationalNumber> testMethod() {
	return DistributedDenseMultivariateModularPolynomial<RationalNumber>();
}
