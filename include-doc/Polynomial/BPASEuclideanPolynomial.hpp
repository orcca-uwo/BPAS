

#ifndef _BPAS_EUCLIDEANPOLY_H_
#define _BPAS_EUCLIDEANPOLY_H_

#include "Utils/TemplateHelpers.hpp"


/**
 * An abstract class defining the interface of a polynomial ring which
 * is also a Euclidean domain. Thus, it is a univariate polynomial
 * over a field.
 *
 * This class is automatically determined to be a superclass of
 * a BPASPolynomial based on a template specialization of the Ring
 * parameter.
 */
template <class Ring>
class BPASEuclideanPolynomial : public BPASGCDPolynomial<Ring>,
                                public BPASEuclideanDomain {}



#endif