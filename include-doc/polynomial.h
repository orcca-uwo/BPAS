#ifndef _BPAS_POLYNOMIAL_H_
#define _BPAS_POLYNOMIAL_H_

/**
 * Abstract Polynomial Classes
 **/

#include "ring.h"
#include "Ring/BPASGCDDomain.hpp"
#include "Symbol/Symbol.hpp"
#include <vector>
#include "Utils/TemplateHelpers.hpp"

/**
 * An abstract class defining the interface of a polynomial over an arbitrary BPASRing.
 * This class is 
 * Polynomials themselves form a ring, BPASRing<Derived> as well as the base
 * ring, Ring, is also a valid ring. 
 **/
template <class Ring> 
class BPASBasePolynomial : public virtual BPASRing
{
	public:

		/**
		 * In addition to the ring arithmetic for Derived defined by BPASRing
		 * polynomials must do arithmetic with their base ring.
		 **/
		virtual Derived& operator= (const Ring&) = 0;
		virtual Derived operator+ (const Ring&) const = 0;
		virtual Derived& operator+= (const Ring&) = 0;
		virtual Derived operator- (const Ring&) const = 0;
		virtual Derived operator- () const = 0;
		virtual Derived& operator-= (const Ring&) = 0;
		virtual Derived operator* (const Ring&) const = 0;
		virtual Derived& operator*= (const Ring&) = 0;

		virtual Integer degree() const = 0;    // total degree
		virtual Ring leadingCoefficient() const = 0; 
		virtual Ring trailingCoefficient() const = 0; 
		virtual bool isConstantTermZero() const = 0;
		virtual Integer numberOfTerms() const = 0;
		virtual Ring content() const = 0;
		virtual Derived primitivePart() const = 0;
};

/** 
 * An abstract class defining the interface of a polynomial over an arbitrary BPASRing.
 * Users should inherit from this class rather than BPASBasePolynomial, BPASIntegralPolynomial, or BPASGCDPolynomial.
 * Inherirance of proper classes and exporting of proper functions is done automatically
 * at compile type through introspection of the template parameter.  
 */
template <class Ring>
class BPASPolynomial : public BPASBasePolynomial<Ring> {};
	
/**
 * An abstract class defining the interface of a univariate polynomial over an arbitrary BPASRing.
 */
template <class Ring>
class BPASUnivariatePolynomial : public virtual BPASPolynomial<Ring>
							
{
	public:
		virtual void differentiate() = 0; // p = dp/dx
		virtual void differentiate(int) = 0; 
		virtual Derived derivative() const = 0; // q = dp/dx
		virtual Derived derivative(int) const = 0;
		virtual Ring evaluate(const Ring& r) const = 0; 
		virtual Derived monicDivide(const Derived&) = 0;
		virtual Derived monicDivide(const Derived&, Derived*) const = 0;
		virtual Derived lazyPseudoDivide(const Derived&, Ring*, Ring*) = 0;
		virtual Derived lazyPseudoDivide(const Derived&, Derived*, Ring*, Ring*) const = 0;
		virtual Derived pseudoDivide(const Derived&, Ring*) = 0;
		virtual Derived pseudoDivide(const Derived&, Derived*, Ring*) const = 0;
		virtual Ring coefficient(int) const = 0;
		virtual void setCoefficient(int, const Ring&) = 0;
		virtual void setVariableName (const Symbol&) = 0;
		virtual Symbol variable() const = 0;
		virtual Derived operator<< (int i) const = 0; // q = p * (x^i);
		virtual Derived& operator<<= (int) = 0; // p = p *(x^i)
		virtual Derived operator>> (int) const = 0; // q = p / (x^i);
		virtual Derived& operator>>= (int) = 0;
};

/**
 * An abstract class defining the interface of a multivariate polynomial over an arbitrary BPASRing.
 */
template <class Ring>
class BPASMultivariatePolynomial : public virtual BPASPolynomial<Ring>
{
	// Monomials are the abelian free monoid generated by Symbols
	public:
		// virtual void setPolynomialRing(const std::vector<Symbol>& v) = 0;
		// virtual std::vector<Symbol> polynomialRing() = 0;

		virtual void differentiate(const Symbol&) = 0; 
		virtual void differentiate(const Symbol&, int) = 0; 
		virtual Derived derivative(const Symbol&) const = 0;
		virtual Derived derivative(const Symbol&, int) const = 0;
		
		// TODO
		// virtual void integrate(const Symbol&) = 0; 
		// virtual void integrate(const Symbol&, int) = 0; 
		// virtual Derived integral(const Symbol&) const = 0;
		// virtual Derived integral(const Symbol&, int) const = 0;

		virtual Derived evaluate(int, const Symbol*, const Ring*) const = 0; 
		virtual Derived evaluate(const std::vector<Symbol>&, const std::vector<Ring>&) const = 0;

		virtual int numberOfVariables() const = 0;
		virtual int numberOfRingVariables() const = 0;
		virtual Integer degree(const Symbol& v) const = 0; 
		// virtual unsigned int degree(const Symbol& v) const = 0; 
		virtual Ring coefficient(int, const int*) const = 0;    			
		virtual Ring coefficient(const std::vector<int>& v) const = 0;
		virtual void setCoefficient(int, const int*, const Ring& r) = 0;  
		virtual void setCoefficient(const std::vector<int>& v, const Ring& r) = 0;

		//set variables in the ring.
		virtual void setRingVariables (const std::vector<Symbol>& xs) = 0;
		//get all variables in the polynomial ring.
		virtual std::vector<Symbol> ringVariables() const = 0;

		//non-zero variables
		virtual std::vector<Symbol> variables() const = 0;

		// virtual Derived content(const std::vector<Symbol>& v) const = 0; 
		// virtual Derived primitivePart(const std::vector<Symbol>& v) const = 0;
		// virtual Derived primitivePart(const std::vector<Symbol>& v, Derived& content) const = 0;
};


/**
 * An abstract class defining the interface of a multivariate polynomial that can be viewed recursively.
 * That is, it can be viewed as a univariate polynomial with multivariate polynomial coefficients.
 */
template <class Ring>
class BPASRecursivelyViewedPolynomial : public virtual BPASMultivariatePolynomial<Ring>
{
	public:
	    virtual Derived initial() const = 0;
	    virtual Symbol mainVariable() const = 0;
		virtual int mainDegree() const = 0;
		virtual Derived rank() const = 0;
 		virtual Derived tail() const = 0;
 		virtual Derived head() const = 0;
 		virtual Derived separant() const = 0;
 		// virtual SparseUnivariatePolynomial<Derived> convertToSUP() const = 0;

};


/**
 * An abstract class defining the interface of a rational function. Domain
 * should be BPASGCDDomain.
 */
template <class Domain>
class BPASRationalFunction : public virtual BPASFieldOfFractions<Domain> {	

};

/**
 * An abstract class defining the interface of a triangular set.
 * A BPASTriangularSet is templated by a BPASRecursivelyViewedPolynomial with coefficients
 * in a BPASField.
 */
template <class Field, class RecursivePoly> 
class BPASTriangularSet
{
	public:
	
		virtual BPASTriangularSet<Field,RecursiveFieldPoly>& operator= (const BPASTriangularSet<Field,RecursiveFieldPoly>&) = 0;
		virtual BPASTriangularSet<Field,RecursiveFieldPoly>& operator= (BPASTriangularSet<Field,RecursiveFieldPoly>&&) = 0;
//		virtual void copy(const BPASTriangularSet<Field,RecursiveFieldPoly>&) = 0;
		virtual int numberOfVariables() const = 0;
		virtual std::vector<Symbol> variables() const = 0;

		virtual RecursiveFieldPoly select(const Symbol&) const = 0;
		virtual void lower(const Symbol&, BPASTriangularSet<Field,RecursiveFieldPoly>&) const = 0;
		virtual void upper(const Symbol&, BPASTriangularSet<Field,RecursiveFieldPoly>&) const = 0;
		virtual RecursiveFieldPoly pseudoDivide (const RecursiveFieldPoly&, std::vector<RecursiveFieldPoly>*, RecursiveFieldPoly*) const = 0;
		virtual RecursiveFieldPoly normalForm (const RecursiveFieldPoly&, std::vector<RecursiveFieldPoly>*) const = 0;
};

/**
 * An abstract class defining the interface of a regular chain.
 * See also BPASTriangularSet.
 */
template <class Field, class RecursivePoly>
class BPASRegularChain : public virtual BPASTriangularSet<Field,RecursivePoly>
{
	public:
	
		virtual BPASRegularChain<Field,RecursiveFieldPoly>& operator= (const BPASRegularChain<Field,RecursiveFieldPoly>&) = 0;
		virtual BPASRegularChain<Field,RecursiveFieldPoly>& operator= (BPASRegularChain<Field,RecursiveFieldPoly>&&) = 0;
//		virtual RecursiveFieldPoly normalize (const RecursiveFieldPoly&) = 0;
};

/**
 * An abstract class defining the interface of a zero-dimensional regular chain.
 */
template <class Field, class RecursivePoly>
class BPASZeroDimensionalRegularChain : public virtual BPASRegularChain<Field,RecursivePoly>
{
	public:
	
		virtual BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>& operator= (const BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>&) = 0;
		virtual BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>& operator= (BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>&&) = 0;
};

#endif
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/

