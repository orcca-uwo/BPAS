
#ifndef _BPAS_POLYNOMIAL_TESTERS_
#define _BPAS_POLYNOMIAL_TESTERS_

#include "../Ring/BPASIntegralDomain.hpp"
#include "../Ring/BPASGCDDomain.hpp"


//forward declares

template <class Ring, class Dervied>
class BPASBasePolynomial;

template <class Ring, class Dervied>
class BPASIntegralPolynomial;

template <class Ring, class Dervied>
class BPASGCDPolynomial;

/**
 * Via conditional inheritance, determine if the ground ring
 * template parameter Ring is an integral domain
 * or not. If so, the resulting polynomial is a BPASIntegralPolynomial
 * and forms an integral domain.
 *
 * This tester is used by BPASPolynomial for compile-time type
 * resolution.
 */
template <class Ring, class Derived>
class BPASIntegralPolynomialTester : public std::conditional<std::is_base_of<BPASIntegralDomain<Ring>, Ring>::value, BPASIntegralPolynomial<Ring, Derived>, BPASBasePolynomial<Ring, Derived>>::type {};


/**
 * Via conditional inheritance, determine if the ground ring
 * template parameter Ring is a GCD domain
 * or not. If so, the resulting polynomial is a BPASGCDPolynomial
 * and forms a GCD domain.
 *
 * This tester is used by BPASPolynomial for compile-time type
 * resolution.
 */
template <class Ring, class Derived>
class BPASGCDPolynomialTester: public std::conditional<std::is_base_of<BPASGCDDomain<Ring>, Ring>::value, BPASGCDPolynomial<Ring, Derived>, BPASIntegralPolynomialTester<Ring, Derived> >::type {};


#endif