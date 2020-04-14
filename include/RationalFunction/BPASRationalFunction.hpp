
#ifndef _BPAS_RATFUNC_H_
#define _BPAS_RATFUNC_H_

#include "../Ring/BPASFieldOfFractions.hpp"


/**
 * An abstract class defining the interface of a rational function. Domain
 * should be BPASGCDDomain.
 */
template <class Domain, class Derived>
class BPASRationalFunction : public virtual BPASFieldOfFractions<Domain, Derived> {

};


#endif
