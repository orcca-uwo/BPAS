
#ifndef _BPAS_EUCLIDEAN_DOMAIN_H_
#define _BPAS_EUCLIDEAN_DOMAIN_H_

#include "BPASGCDDomain.hpp"
#include "../Utils/TemplateHelpers.hpp"

/**
 * An abstract class defining the interface of a Euclidean domain.
 */
template <class Derived>
class BPASEuclideanDomain : public virtual BPASGCDDomain<Derived> {

public:

	/** 
	 * Get the euclidean size of *this.
	 */
	virtual Derived euclideanSize() const = 0;
	
	/**
	 * Perform the eucldiean division of *this and b. Returns the 
	 * remainder. If q is not NULL, then returns the quotient in q. 
	 */ 
	virtual Derived euclideanDivision(const Derived& b, Derived* q = NULL) const = 0;
	
	/**
	 * Perofrm the extended euclidean division on *this and b. 
	 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
	 */
	virtual Derived extendedEuclidean(const Derived& b, Derived* s = NULL, Derived* t = NULL) const = 0;
	
	/**
	 * Get the quotient of *this and b.
	 */
	virtual Derived quotient(const Derived& b) const = 0;

	/** 
	 * Get the remainder of *this and b.
	 */
	virtual Derived remainder(const Derived& b) const = 0;

	/** 
	 * Get the remainder of *this and b;
	 */
	virtual Derived operator%(const Derived& b) const = 0;

	/**
	 * Assign *this to be the remainder of *this and b.
	 */
	virtual Derived& operator%=(const Derived& b) = 0;


};

#endif
