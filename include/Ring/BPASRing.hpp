#ifndef _BPAS_RING_H_
#define _BPAS_RING_H_

#include <gmpxx.h>
#include <vector>
#include "../ExpressionTree/ExpressionTree.hpp"
#include <sstream>

/** 
 * An enumeration which describes the properties that a particular ring has.
 * This operates as a bit-wise mask such that, given an unsigned int representing
 * a particular ring's properties, one can do :
 * if (ringProp & (PRIME_FIELD | FINITE_FILE)) 
 * to test if a Ring is a (in this instance) a prime field or a finite field.   
 */
typedef enum RingProperty {
    COMMUTATIVE_RING = 0x001, 
    INTEGRAL_DOMAIN = 0x003,
    GCD_DOMAIN = 0x007,
    UNIQUE_FACTORIZATION_DOMAIN = 0x00f,
    PRINICPAL_IDEAL_DOMAIN = 0x01f,
    EUCLIDEAN_DOMAIN = 0x03f,
    FIELD = 0x07f,
    PRIME_FIELD = 0x0ff,
    FINITE_FIELD = 0x1ff,
    SMALL_PRIME_FIELD = 0x3ff,
    COMPLEX_FIELD = 0x47f
} RingProperty;

/**
 * Class enapsulates a set of RingProperty(s) that the Ring may have. 
 */
class RingProperties {
private:
    unsigned int prop;
    // int n;

public:

    /**
     * Default constructor which specifies no properties. 
     */
    RingProperties();

    /**
     * Construct a RingProperties from a RingProperty enum element.
     */
    RingProperties(RingProperty p);

    /**
     * Construct a RingProperties from a collection of RingProperty elements.
     */
    RingProperties(std::vector<RingProperty> v);

    /**
     * Determine if *this RingProperties has a RingProperty.
     */
    inline bool has(RingProperty p);

    /**
     * Determine if *this RingProperties has all properties defined
     * by the other RingProperties p.
     */
    inline bool has(const RingProperties& p);
};



/**
 * An abstract class defining the interface of a commutative ring.
 *
 * The template Derived is a concrete class derived from (implemeneting the 
 * interface of) BPASRing. This is the "curiously recurring template pattern" (CTRP).
 * This pattern is used among all sub-classes of BPASRing.
 */
template <class Derived>
class BPASRing : public virtual ExpressionTreeConvert {
public:

    /**
     * Static element describing the properties of this ring class.
     */
    static RingProperties properties;        

    /**
     * The characteristic of this ring class.
     */
    virtual mpz_class characteristic() {
        return 0;
    }

    /**
     * Determine if *this ring element is zero, that is the additive identity.
     *
     * returns true iff *this is zero.
     */
    virtual bool isZero() const = 0;

    /**
     * Make *this ring element zero.
     */
    virtual void zero() = 0;

    /**
     * Determine if *this ring element is one, that is the multiplication identity.
     *
     * returns true iff *this is one.
     */
    virtual bool isOne() const = 0;
    
    /**
     * Make *this ring element one.
     */
    virtual void one() = 0;
    
    /**
     * Obtain the unit normal (a.k.a canonical associate) of an element. 
     * If either parameters u, v, are non-NULL then the units are returned such that 
     * b = ua, v = u^-1. Where b is the unit normal of a, and is the returned value.
     */
    virtual Derived unitCanonical(Derived* u = NULL, Derived* v = NULL) const = 0;

    /**
     * Copy assignment.
     */
    virtual Derived& operator= (const Derived&) = 0;
    
    /**
     * Addition.
     */
    virtual Derived operator+ (const Derived&) const = 0;

    /**
     * Addition assignment.
     */
    virtual Derived& operator+= (const Derived&) =0;
    
    /**
     * Subtraction.
     */
    virtual Derived operator- (const Derived&) const = 0;
    
    /**
     * Subtraction assignment.
     */
    virtual Derived& operator-= (const Derived&) = 0;

    /**
     * Negation.
     */
    virtual Derived operator- () const = 0;

    /**
     * Multiplication.
     */
    virtual Derived operator* (const Derived&) const = 0;

    /**
     * Multiplication assignment.
     */
    virtual Derived& operator*= (const Derived&) = 0;
    
    /**
     * Exponentiation.
     */
    virtual Derived operator^ (long long int e) const = 0;
    
    /**
     * Exponentiation assignment.
     */
    virtual Derived& operator^= (long long int e) = 0;

    /**
     * Equality test,
     *
     * returns true iff equal
     */
    virtual bool operator== (const Derived&) const = 0;

    /**
     * Inequality test,
     *
     * returns true iff not equal.
     */
    virtual bool operator!= (const Derived&) const = 0;

    /**
     * Print the Ring element.
     *
     * Derived classes may override this to get custom printing that may
     * be more expressive (and prettier) than expression tree printing.
     */
    virtual void print(std::ostream& ostream) const {
        ostream << convertToExpressionTree().toString();
    }

    /**
     * Convert the Ring element to a string.
     *
     * Simple delegation of printing to a stringstream to obtain a string.
     * Overriding the print method is sufficient for sub-classes
     * to make use of this method. 
     *
     * returns the string representation of the Ring element. 
     */
    virtual std::string toString() const {
        std::stringstream ss;
        print(ss);
        return ss.str();
    }

    /**
     * Output operator. 
     *
     * Defines a to string conversion. 
     */
    friend std::ostream& operator<< (std::ostream& ostream, const Derived& d) {
        d.print(ostream);
        return ostream;
    }  

    friend std::ostream& operator<< (std::ostream& ostream, Derived&& d) {
        d.print(ostream);
        return ostream;
    }
    
    // virtual bool isNegativeOne() = 0;
    // virtual void negativeOne() = 0;
    // virtual int isConstant() = 0;
	// static bool isPrimeField;
	// static bool isSmallPrimeField;
	// static bool isComplexField;
};

#endif
