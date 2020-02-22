
#ifndef _MAPLE_INTERFACE_HPP_
#define _MAPLE_INTERFACE_HPP_

#include <vector>
#include <string>

#include "../DataStructures/Factors.hpp"
#include "MapleWorkerThread.hpp"


/*
 * Forward declares.
 */
class SparseMultivariateIntegerPolynomial;



/**
 * A generic abstract base class for processing requests to maple
 * via the MapleWorkerThread. Different concrete implementations
 * may create their own way of communicating with that MapleWorkerThread.
 *
 * Methods provided herein take BPAS objects, converting them to 
 * sequences of strings to be executed as a sequence of maple commands.
 *
 * Concrete derived classes should simply provide their own implementation
 * of sendCommand.
 *
 */
class MapleInterface {

protected: 

	/**
	 * Send the command encoded as the vector request to the underlying MapleWorkerThread.
	 * The MapleWorkerThread returns its response in the response vector.
	 *
	 * request: vector of strings encoding the request for maple to process.
	 * response: the vector of strings encoding the result of maple's processing.
	 * 
	 */
	virtual void sendCommand(const std::vector<std::string>& request, std::vector<std::string>& response) = 0;

public: 

	/**
	 * Get the GCD of two sparse integer polynomials.
	 *
	 * @param a: the first polynomial.
	 * @param b: the second polynomial.
	 *
	 * returns the gcd of a and b.
	 */
	SparseMultivariateIntegerPolynomial gcd(const SparseMultivariateIntegerPolynomial& a, const SparseMultivariateIntegerPolynomial& b);


	/**
	 * Factor a sprase integer polynomial.
	 * @param p: the polynomial to factor.
	 * 
	 * returns the factors of p as a Factors object. See the Factors class.
	 */
	Factors<SparseMultivariateIntegerPolynomial> factor(const SparseMultivariateIntegerPolynomial& p);

	/**
     * Validate a triangularize result.
     * request should be {"triangularizeValiadte", <isLazard>, <F>, <rc>, <vars>, <our_results>}
     * response will be returned as {<boolean>}
     *
	 * returns true iff a response was generated. 
	 */
	bool TriangularizeValidation(bool isLazard, std::vector<std::string>& inputs);
		
};

#endif
