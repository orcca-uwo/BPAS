#ifndef _SMQPALTARRAY_H_ 
#define _SMQPALTARRAY_H_

#include "../polynomial.h"
#include "../Interval/interval.h"
#include "urpolynomial.h"
#include "../IntegerPolynomial/mzpolynomial.hpp"
#include "../RingPolynomial/upolynomial.h"
#include "SMQP_CppSupport-AA.hpp"
#include "SMQP_Support-AA.h"
#include "SMQP_Support_Test-AA.h"
#include "SMQP_Support_Recursive-AA.h"
#include <gmpxx.h>
#include "../ExpressionTree/ExpressionTree.hpp"
#include "../DataStructures/Factors.hpp"
#include "../TriangularSet/triangularset.hpp"
#include "../IntegerPolynomial/DUZP_Support.h"
#include "../Parser/bpas_parser.h"

class SparseMultivariateIntegerPolynomial;

/**
 * An element of the SLP of a rational number polynomial.
 */
class SLPRepresentation {

	public:
		union CoefOrInt {
			RationalNumber* c;
			int i;
		};

		int op;			// 0: '+'; 1: '*'.
		int type;		// 0: coef & variate;
						// 1: coef & result;
						// 2: variate & result;
						// 3: result & result;
						// 4: coef (constant);
		CoefOrInt a;	// variate or result or coef position
		int b;			// variate or result or coef position
		Interval res;	// Storing the result

		SLPRepresentation() {
			a.c = NULL;
		}

		SLPRepresentation(const SLPRepresentation& r) {
			op = r.op;
			type = r.type;
			b = r.b;
			res = r.res;
			if (type == 0 || type == 1 || type == 4) {
				a.c = new RationalNumber(*(r.a.c));
			} else {
				a.i = r.a.i;
			}   
		}

		~SLPRepresentation() {
			if (type == 0 || type == 1) {
				delete a.c;
			}
		}
};

/** 
 * A multivariate polynomial with RationalNumber coefficients represented sparely.
 * Only non-zero coefficients are encoded.
 */
class SparseMultivariateRationalPolynomial: public BPASRecursivelyViewedPolynomial<RationalNumber,SparseMultivariateRationalPolynomial> {

	private:
		mutable AltArr_t* poly;
		int nvar;   		  //number of variables
		Symbol* names;   //list of strings representing the variables. 
							  //names[0] is 1 if "system-generated" variables
							  //names[0] is 9 if "user-specified" variables;


		friend class SparseMultivariateIntegerPolynomial;

		/**
		 * Determine if *this and b have compatible variable orderings.
		 * xs returns the mapping of both *this and b to the compatible superset
		 * returns true iff compatible.
		 */
		bool isOrderedRing(const SparseMultivariateRationalPolynomial& b, std::vector<int>& xs) const;

		/** 
		 * Rearrange exponent vectors in place and then re-sort the polynomial.
		 */
		void reorderVarsInPlace(int varmap[]);

		/**
		 * Rearrange exponent vectors, expanding as needed, and then re-sort the polynomial.
		 */
		void expandVarsInPlace(int vars, Symbol* newvars, int varmap[]);

		/** 
		 * Returns a copy of *this under the new variable ordering supplied. 
		 * varmap is such that this.names[i] = newvars[varmap[i]]
		 * Returns an SMQP equal to *this but expended to newvars.
		 */
		SparseMultivariateRationalPolynomial expandVariables(int vars, Symbol* newvars, int varmap[]) const;

		void preparePolysForSRC(const SparseMultivariateRationalPolynomial& q, const Symbol& v, std::vector<Symbol>& superRing, bool sameRing, int* varMap, int& swapIdx, AltArrZ_t** ppZ, AltArrZ_t** qqZ) const;
    
	public:

                
		SparseMultivariateRationalPolynomial subresultantGCD (const SparseMultivariateRationalPolynomial& q) const;
    
		std::vector<SLPRepresentation> slp;

		/* Constructors */

		/**
		 * Construct a multivariate polynomial
		 *
		 **/
		SparseMultivariateRationalPolynomial();

		/**
		 * Construct a multivariate polynomial with specific number of variables
		 *
		 * @param v: Number of variables
		 **/
		SparseMultivariateRationalPolynomial(int v);

		/**
		 * Construct with a variable name such that f(x) = x;
		 *
		 * @param x: The variable name
		 **/
		SparseMultivariateRationalPolynomial (const Symbol& x);

		/**
		 * Construct a polynomial by parsing the string str.
		 */
		SparseMultivariateRationalPolynomial (const std::string& str);

        /**
		 * Construct an SMQP given the head Node*. The SMQP takes ownership
		 * of the Node list.
		 */
		SparseMultivariateRationalPolynomial(AltArr_t* aa, int vars, Symbol* varNames);

		/**
		 * Copy Constructor.
		 * 
		 * Does not reuse underlying memory allocated by b. 
		 *
		 * @param b: A sparse multivariate polynomial
		 **/
		SparseMultivariateRationalPolynomial(const SparseMultivariateRationalPolynomial& b);

		/**
		 * Move Constructor.
		 *
		 * @params b: The r-value reference polynomial.
		 */
		SparseMultivariateRationalPolynomial(SparseMultivariateRationalPolynomial&& b);

		/**
		 * Create a SMZP from an SMQP.
		 */
		SparseMultivariateRationalPolynomial(const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Create a SMQP from an Integer. 
		 */
		SparseMultivariateRationalPolynomial(const Integer& r, int nvar = 0);

		/**
		 * Create a SMQP from a RationalNumber. 
		 */
		SparseMultivariateRationalPolynomial(const RationalNumber& r, int nvar = 0);

		/**
		 * Create a SMQP from a univariate rational polynomial. 
		 * 
		 * @param p: A SUQP polynomial. 
		 **/
		SparseMultivariateRationalPolynomial (const DenseUnivariateRationalPolynomial& p);

		/**
		 * Construct from a SUP<SMQP> polynomial such that the resulting SMQP
		 * has main variable equal to the variable of s. 
		 *
		 * @param s: The SUP<SMQP> polynomial
		 **/
    	SparseMultivariateRationalPolynomial (const SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>& s);
	
		/**
		 * Destroy the polynomial and underlying node memory.
		 **/
		~SparseMultivariateRationalPolynomial();


		/* BPASRing interface */
        static mpz_class characteristic;
        static RingProperties properties;
        // static bool isPrimeField;
        // static bool isSmallPrimeField;
        // static bool isComplexField;

        /** 
         * Is this polynomial zero. 
         * returns true iff this polynomial encodes 0.
         */
		bool isZero() const;

		/**
		 * Set this polynomial to zero. Maintains existing variable ordering.
		 */ 
        void zero();

        /**
         * Is this polynomial one.
         * return true iff this polynomial encodes 1.
         */
        bool isOne() const;

        /** 
         * Sets this polynomial to one. Maintains existing variable ordering.
         */
        void one();

        /**
         * Is this polynomial negative one.
         * returns true iff this polynomial encodes -1.
         */
        bool isNegativeOne() const;

        /** 
         * Sets this polynomial to -1. Maintains existing variable ordering.
         */
        void negativeOne();

        /**
         * Determine if this polynomial encodes a constant. 
         * returns 0 if not a constant, 1 if a constant >= 0, -1 if a negative contstant. 
         */
        int isConstant() const;
    
    	/**
	     * Obtain the unit normal (a.k.a canonical associate) of an element. 
	     * If either parameters u, v, are non-NULL then the units are returned such that 
	     * b = ua, v = u^-1. Where b is the unit normal of a, and is the returned value.
	     */  
		inline SparseMultivariateRationalPolynomial unitCanonical(SparseMultivariateRationalPolynomial* u, SparseMultivariateRationalPolynomial* v) const {
        	RationalNumber lead = leadingCoefficient();
        	RationalNumber leadInv = lead.inverse();
        	if (u != NULL) {
        		*u = leadInv;
        	}
        	if (v != NULL) {
        		*v = lead;
        	}
        	return (*this * leadInv);
        }

        /* BPASPolynomial interface */

        /**
         * Assign this polynomail to equal the specified.
         */
        SparseMultivariateRationalPolynomial& operator= (const SparseMultivariateRationalPolynomial& b);

		/**
         * Movement assignment: move b to be this polynomail.
         */
        SparseMultivariateRationalPolynomial& operator= (SparseMultivariateRationalPolynomial&& b);

        /**
         * Assign this polynomial to be equal to the constant ring element.
         */
	    SparseMultivariateRationalPolynomial& operator= (const RationalNumber& r);
	    
        /**
         * Add two SMQP polynomials together, *this and the specified.
         */
		SparseMultivariateRationalPolynomial operator+ (const SparseMultivariateRationalPolynomial& b) const;
		SparseMultivariateRationalPolynomial operator+ (SparseMultivariateRationalPolynomial&& b) const;

		// /**
		//  * Add the polynomails a and b and return the sum.
		//  */
		// friend SparseMultivariateRationalPolynomial operator+ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

		/**
		 * Update *this by adding the specified polynomail to it.
		 */
		SparseMultivariateRationalPolynomial& operator+= (const SparseMultivariateRationalPolynomial& b);

		/**
		 * Subtract the specified polynomail from *this
		 */
		SparseMultivariateRationalPolynomial operator- (const SparseMultivariateRationalPolynomial& b) const;
		SparseMultivariateRationalPolynomial operator- (SparseMultivariateRationalPolynomial&& b) const;

		// /**
		//  * Subtract the polynomial b from a and return the difference.
		//  */
		// friend SparseMultivariateRationalPolynomial operator- (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

		/**
		 * Unary operation, return  *this * -1.
		 */
		SparseMultivariateRationalPolynomial operator- () const;

		/**
		 * Update *this by subtracting the specified polynomial from it.
		 */
		SparseMultivariateRationalPolynomial& operator-= (const SparseMultivariateRationalPolynomial& b);

		/**
		 * Multiply *this by the specified polynomail.
		 */
		SparseMultivariateRationalPolynomial operator* (const SparseMultivariateRationalPolynomial& b) const;
		SparseMultivariateRationalPolynomial operator* (SparseMultivariateRationalPolynomial&& b) const;

		// /** 
		//  * Multiply the polynomials a and b, returning their product.
		//  */
		// friend SparseMultivariateRationalPolynomial operator* (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

		/**
		 * Update this by multiplying by the specified polynomail.
		 */
		SparseMultivariateRationalPolynomial& operator*= (const SparseMultivariateRationalPolynomial& b);

		/**
		 * Divide *this by the specified polynomial.
		 */
		SparseMultivariateRationalPolynomial operator/ (const SparseMultivariateRationalPolynomial& b) const;
		SparseMultivariateRationalPolynomial operator/ (SparseMultivariateRationalPolynomial&& b) const;

		// friend SparseMultivariateRationalPolynomial operator/ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

		/**
		 * Update *this by dividing by the specified polynomial.
		 */
		SparseMultivariateRationalPolynomial& operator/= (const SparseMultivariateRationalPolynomial& b);

		/**
		 * Exponentiate *this by the input exponent integer.
		 * Treats negative exponents as positive.
		 */ 
		SparseMultivariateRationalPolynomial operator^ (long long int e) const;

		/**
		 * Update *this by exponentiating this to the input integer.
		 * Treats negative exponents as positive.
		 */
		SparseMultivariateRationalPolynomial& operator^= (long long int e);
		
		/**
		 * Determine if *this is equal to the specified polynomial.
		 * This takes into account the variable ordering on both poylnomials 
		 * in such a way that the same polynomial under different variable orderings
		 * are NOT equal.
		 */ 
		bool operator== (const SparseMultivariateRationalPolynomial& b) const;

		/**
		 * Determine if *this is not equal to the specified polynomial.
		 * This takes into account the variable ordering on both poylnomials 
		 * in such a way that the same polynomial under different variable orderings
		 * are NOT equal.
		 */
		bool operator!= (const SparseMultivariateRationalPolynomial& b) const;
		
		/**
		 * Output the string representation of *this to the input ostream.
		 */
		void print(std::ostream& os) const;

		/**
		 * Parse polynomial from in stream. Exactly one line is parsed. 
		 */
		friend std::istream& operator>>(std::istream& in, SparseMultivariateRationalPolynomial& p);

		/** 
		 * Parse a polynomial from string str and place in *this.
		 */
 		void fromString(const std::string& str);

		/**
		 * Given a polynomial q, compute the subresultant chain between this and q, viewing both polynomial 
		 * recursively with main variable v. 
		 *
		 * @note, if this and q do not exist in the same ambient space, the space of p will define the ordering of the subresultants. 
		 * 
		 * @param q, the other polynomial for which to compute the subresultant chain.
		 * @param v, the main variable to be used when computing the subresultant chain.
		 * @param filled, if false, degenerative cases are not returned and so indices in the returned chain do not necessarily match subresultant degrees
		 *
		 * @return a vector containing the subresultant chain whereby index 0 is resultant and subresultant degrees increase with index.
		 *
		 */
		std::vector<SparseMultivariateRationalPolynomial> subresultantChain(const SparseMultivariateRationalPolynomial& q, const Symbol& v, bool filled = 1) const;	    

		/**
		 * Subresultant Chain with respect to the leaving variable of this.
		 * Return the list of subresultants
 	     **/
 		std::vector<SparseMultivariateRationalPolynomial> subresultantChain (const SparseMultivariateRationalPolynomial& q, int filled = 1) const;

		/**
		 * Given a polynomial q, compute the subresultant of index idx between this and q, viewing both polynomial 
		 * recursively with main variable v. This function in fact computes the subresultant of index idx and idx+1
		 * and so both are returned, unless idx+1 does not exist (e.g. idx = this->degree(v)).
		 *
		 * Optionally, this function also returns the principleCoefs of the subresultants of index i to this->degree(v).
		 *
		 * @note, if the subresultant of index idx+1 is degenerative then many subresultants will be returned until the first non-degenerative one is found.
		 * @note, if this and q do not exist in the same ambient space, the space of p will define the ordering of the subresultants. 
		 * 
		 * @param q, the other polynomial for which to compute the subresultant chain.
		 * @param v, the main variable to be used when computing the subresultant chain.
		 *
		 * @return a vector containing the subresultant chain whereby index 0 is requested subresultant idx subresultant degrees increase with index.
		 *
		 */
		std::vector<SparseMultivariateRationalPolynomial> subresultantChainAtIdx (const SparseMultivariateRationalPolynomial& q, const Symbol& v, int idx = 0, std::vector<SparseMultivariateRationalPolynomial>* principleCoefs = NULL) const;
	
		/**
 		* Extended Subresultant Chain
 		* Return the list of subresultants with Besout Coefficients
 		**/
	    std::vector<std::vector<SparseMultivariateRationalPolynomial>> exSubresultantChain (const SparseMultivariateRationalPolynomial& q, const Symbol& v) const;	    
	    
	    SparseMultivariateRationalPolynomial resultant (const SparseMultivariateRationalPolynomial& q, const Symbol& v) const;
    
        /**
         * Get GCD between *this and b.
         */
        SparseMultivariateRationalPolynomial gcd(const SparseMultivariateRationalPolynomial& b) const;

        /**
         * Get the GCD between *this and b as a primitive polynomial.
         */
        SparseMultivariateRationalPolynomial primitiveGCD(const SparseMultivariateRationalPolynomial& b) const;
			
		/**
		 * Compute squarefree factorization of *this with respect to all of its variables.
		 */
		Factors<SparseMultivariateRationalPolynomial> squareFree() const;

		/**
		 * Compute squarefree factorization of *this with respect to the list of variables, vars.
		 */
		Factors<SparseMultivariateRationalPolynomial> squareFree(const std::vector<Symbol>& vars) const;

		/**
		 * Computes the square free part of this polynomail. That is, the polynomial of this
		 * divided by all square factors. This is with respect to all variables.
		 */
		SparseMultivariateRationalPolynomial squareFreePart() const;

		/**
		 * Computes the square free part of this polynomail. That is, the polynomial of this
		 * divided by all square factors. This is with respect to all variables.
		 */
		SparseMultivariateRationalPolynomial squareFreePart(std::vector<Symbol>& vars) const;

		/**
		 * Get the content with respect to all variables. That is, a single
		 * rational number. The content here is one such that this/content is
		 * an integer polynomial with content of 1.
		 *
		 * Moreover, the content is one such that the leading coefficient
		 * of the corresponding primitive part is positive.
		 */
		RationalNumber content() const;

		/**
		 * Get the content of this polynomial with respect to the variables in 
		 * the input vector v.
		 *
		 * Moreover, the content is one such that the leading coefficient
		 * of the corresponding primitive part is positive.
		 */
		SparseMultivariateRationalPolynomial content(const std::vector<Symbol>& v) const;

		/**
		 * Get the primitive part with respect to all variables, returned as an SMZP.
		 * This is equivalent to this / content();
		 */
		SparseMultivariateIntegerPolynomial primitivePartSMZP() const;

		/**
		 * Get the primitive part with respect to all variables, returned as an SMZP.
		 * This is equivalent to this / content();
		 *
		 * Simultaneously returns the rational number content in the parameter content.
		 */
		SparseMultivariateIntegerPolynomial primitivePartSMZP(RationalNumber& content) const;

		/**
		 * Get the primitive part with respect to all variables. 
		 * This is equivalent to this / content().
		 */
		SparseMultivariateRationalPolynomial primitivePart() const;

		/**
		 * Get the primitive part with respect to the variable s. 
		 * This is equivalent to this / content(s).
		 */
		SparseMultivariateRationalPolynomial primitivePart(const Symbol& s) const;

		/**
		 * Get the primitive part with respect to all variables. 
		 * This is equivalent to this / content().
		 *
		 * Simultaneously returns the rational number content in the parameter content.
		 */
		SparseMultivariateRationalPolynomial primitivePart(RationalNumber& content) const;

		/**
		 * Get the primitive part with respect to the variables in the vector v. 
		 * 
		 * returns the primitive part.
		 */
		SparseMultivariateRationalPolynomial primitivePart(const std::vector<Symbol>& v) const;

		/**
		 * Get the primitive part with respect to the variables in the vector v. 
		 * Returns the corresponding content in the content reference.
		 * returns the primitive part.
		 */
		SparseMultivariateRationalPolynomial primitivePart(const std::vector<Symbol>& v, SparseMultivariateRationalPolynomial& content) const; 

        SparseMultivariateRationalPolynomial mainPrimitivePart() const;

        SparseMultivariateRationalPolynomial mainPrimitivePart(SparseMultivariateRationalPolynomial& content) const;
	
		/**
		 * Get the leading coefficient of *this with respect to the main variable.
		 *
		 * returns the initial.
		 */
		SparseMultivariateRationalPolynomial initial() const;

		/**
		 * Get the main variable. That is, the highest-ordered variable with
		 * positive degree.
		 *
		 * returns the main variable.
		 */
		inline Symbol mainVariable() const {
			return leadingVariable();
		}

		/**
		 * Get the degree of the main variable.
		 *
		 * returns the degree.
		 */
		int mainDegree() const;

		/**
		 * Get the rank of this polynomial.
		 * That is, the main variable raised to the main degree. 
		 *
		 * returns the rank.
		 */
		SparseMultivariateRationalPolynomial rank() const;

		/**
		 * Get the head of this polynomial. That is, the initial multiplied
		 * by the rank.
		 *
		 * returns the head.
		 */
		SparseMultivariateRationalPolynomial head() const;

		/**
		 * Get the tail of this polynomial. That is, This - this.head().
		 *
		 * returns the tail.
		 */
		SparseMultivariateRationalPolynomial tail() const;

		SparseMultivariateRationalPolynomial separant() const; 

		/* BPASMultivariatePolynomial interface */

		/**
		 * Get the number of variables in this polynomial.
		 */
		int numberOfVariables() const;

		/**
		 * Get the number of variables in this polynomial ring.
		 */
		inline int numberOfRingVariables() const {
			return ringVariables().size();
		}

		/**
		 * Get the number of non-zero terms 
		 */
		Integer numberOfTerms() const;
		
		/** 
		 * Total degree.
		 */
		Integer degree() const;

		/**
		 * Get the degree of a variable 
		 */
		Integer degree(const Symbol&) const;
		
		/**
		 * Get the leading coefficient 
		 */
		RationalNumber leadingCoefficient() const;
		
		/**
		 * Get the trailing coefficient.
		 */
		RationalNumber trailingCoefficient() const;

		/**
		 * Get a coefficient, given the exponent of each variable 
		 */
		RationalNumber coefficient(int, const int*) const;

		inline RationalNumber coefficient(const std::vector<int>& v) const {
			return coefficient(v.size(), v.data());
		}
		
		/**
		 * Set a coefficient, given the exponent of each variable
		 */
		void setCoefficient(int, const int*, const RationalNumber&);

		inline void setCoefficient(const std::vector<int>& v, const RationalNumber& r) {
			setCoefficient(v.size(), v.data(), r);
		}


		/**
		 * Set variables' names.
		 *
		 * This method can be used to shrink, expand, re-order, and re-name the variables
		 * of the polynomial ring in which this polynomial belongs. 
		 *
		 * Any symbol in the input vector that matches the current variables
		 * is considered a re-order. Non-matching symbols are considered a re-name
		 * and this renames is done in the order they appear by increasing index.
		 * For example: Q[x,y,z] -> Q[y,s,t] will have y reordered to be the 
		 * first variable, x renamed to s, and z renamed to t. 
		 *
		 * If the size of the input vector is greater than the current number of variables
		 * then the polynomial ring is being expanded. Matching variables from this
		 * polynomial are re-ordered as they appear in the input vector. Current variables
		 * that have no match in the input vector are re-named in order of increasing index
		 * of unused variables from the input vector. 
		 * For example: Q[x,y,z] -> Q[s,t,z,x] will have y named to s, t a new variable, 
		 * and z and x re-ordered accordingly. 
		 *
		 * If the size of the input vector is less than the current number of variables
		 * then the polynomial ring is shrink to remove variables. Variables in the 
		 * input vector that are also found to be in the current polynomial ring 
		 * are matched and re-ordered as necessary.  
		 *
		 * It is invalid to remove a variable which is non-zero in this polynomial.
		 */
		void setRingVariables (const std::vector<Symbol>&);
		
		/**
		 * Get variable names of all variables available to this polynomial, 
		 * even those that have zero degree.
		 */
		std::vector<Symbol> ringVariables() const;

		/**
		 * Get variable names of variables with non-zero degree;
		 */
		std::vector<Symbol> variables() const;


		/* SMQP-specific */

		/**
		 * Determine if *this is equal to b.
		 * This takes into account the variable ordering on both poylnomials 
		 * in such a way that the same polynomial under different variable orderings
		 * are NOT equal.
		 */
		bool isEqual(const SparseMultivariateRationalPolynomial& b) const;

		/**
		 * Convert current object to its k-th derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the derivative, k > 0
		 **/ 
    	inline void differentiate(const Symbol& s, int k) {
    		*this = this->derivative(s, k);
    	}

		/**
		 * Convert current object to its derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/ 
    	inline void differentiate(const Symbol& s) {
    		this->differentiate(s,1);
    	}

		/**
		 * Return k-th derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the k-th derivative, k > 0
		 **/ 
    	SparseMultivariateRationalPolynomial derivative(const Symbol& s, int k) const;		
		
		/**
		 * Compute derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/ 
		inline SparseMultivariateRationalPolynomial derivative(const Symbol& s) const {
    	 	return this->derivative(s,1);
        } 

        /**
		 * Convert current object to its k-th integral
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the derivative, k > 0
		 **/ 
    	inline void integrate(const Symbol& s, int k) {
    		*this = this->integral(s, k);
    	}

		/**
		 * Convert current object to its derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/ 
    	inline void integrate(const Symbol& s) {
    		this->integrate(s,1);
    	}

		/**
		 * Return k-th integral. If Symbol s is not a symbol contained in this
		 * polynomial then it becomes the main variable of the result and
		 * integration proceeds as expected.
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the k-th derivative, k > 0
		 **/ 
    	SparseMultivariateRationalPolynomial integral(const Symbol& s, int k) const;		
		
		/**
		 * Compute integral
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/ 
		inline SparseMultivariateRationalPolynomial integral(const Symbol& s) const {
    	 	return this->integral(s,1);
        } 

		/**
		 * Evaluate f(x)
		 *
		 * @param syms: Array of Symbols to evaluate at corresponding xs
		 * @param xs: Evaluation points
		 **/
    	inline SparseMultivariateRationalPolynomial evaluate(int n, const Symbol* syms, const RationalNumber* xs) const {
    		std::vector<Symbol> vecSyms;
    		std::vector<RationalNumber> vecRats;
    		vecSyms.reserve(n);
    		vecRats.reserve(n);
    		for (int i = 0; i < n; ++i) {
    			vecSyms.push_back(syms[i]);
    			vecRats.push_back(xs[i]);
    		}
    		return evaluate(vecSyms, vecRats);
    	}

		/**
		 * Evaluate *this polynomial given the variables and their values in vars, and values.
		 * vars must be a list of variables which are a (not necessarily proper) subset.
		 * vars and values should have matching indices. i.e. the values of vars[i] is values[i].
		 *
		 * returns a new SMQP where all variables in vars have been evaluated using the values given. 
		 */
		SparseMultivariateRationalPolynomial evaluate(const std::vector<Symbol>& vars, const std::vector<RationalNumber>& values) const;

		/**
		 * Find the interpolating polynomial for the points and values specified by each vector, respectively.
		 * The points vector is a vector of vectors, where each element of the outer vector represents a multi-dimensional point. 
		 * Each multi-dimensional point should have the same number of dimensions and in the same order.
		 *
		 * returns the interpolating polynomial.
		 */
		static SparseMultivariateRationalPolynomial interpolate(const std::vector<std::vector<RationalNumber>>& points, const std::vector<RationalNumber>& vals);

		/**
		 * Divide this by polynomial b, returning the quotient and remainder in q and r, respectively.
		 *
		 * returns a boolean indicating if the division was exact.
		 */
		bool divide(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial& q, SparseMultivariateRationalPolynomial& r) const;

		/**
		 * Get the remainder of *this divided by b.
		 */
		SparseMultivariateRationalPolynomial operator% (const SparseMultivariateRationalPolynomial& b) const;

		/**
		 * Update *this by setting it to the remainder of *this divided by b.
		 */ 
		SparseMultivariateRationalPolynomial& operator%= (const SparseMultivariateRationalPolynomial& b);

		
		inline SparseMultivariateRationalPolynomial operator% (const Integer& i) const {
			SparseMultivariateRationalPolynomial ret = *this;
			ret %= i;
			return ret;
		}

		inline SparseMultivariateRationalPolynomial& operator%= (const Integer& i) {
			if (!this->isZero()) {
				applyModuloSymmetricPrimPart_AA_inp(this->poly, i.get_mpz_t());
			}
			return *this;
		}

		/** 
		 * Pseudo divide this by b. The remainder is returned.
		 * if parameter quo is not null then the quotient is returned in quo.
		 * if parameter mult is not null then the multiplier is set to the initial of b 
		 * raised to the power of degree(c, mvar(c)) - degree(b, mvar(c)) + 1.
		 *
		 * returns the pseudo remainder.
		 */
		SparseMultivariateRationalPolynomial pseudoDivide(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial* quo = NULL, SparseMultivariateRationalPolynomial* mult = NULL, bool lazy = 0) const;

		/** 
		 * Pseudo divide this by b. The remainder is returned.
		 * This function assumes that both this and b are primitive, and thus
		 * can be viewed as integer polynomials. The psuedo-division is then performed
		 * over the integers and the quotient and remainder returned also so they are primitive.
		 *
		 * If parameter quo is not null then the quotient is returned in quo.
		 * If parameter mult is not null then the multiplier is set to the initial of b 
		 * raised to the power of degree(c, mvar(c)) - degree(b, mvar(c)) + 1.
		 *
		 * @return the pseudo remainder.
		 */
		SparseMultivariateRationalPolynomial pseudoDivideBySMZP(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial* quo = NULL, SparseMultivariateRationalPolynomial* mult = NULL, bool lazy = 0) const;


		/**
		 * Add *this and a ratNum_t.
		 */		
		inline SparseMultivariateRationalPolynomial operator+ (const ratNum_t& r) const {
			return (*this + RationalNumber(r));
		}

		/**
		 * Add *this and the RationalNumber r.
		 */
		SparseMultivariateRationalPolynomial operator+ (const RationalNumber& r) const;

		/**
		 * Add ratNum_t r and SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator+ (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) {
			return (b + r);
		}

		/**
		 * Update *this by adding r
		 */
		inline SparseMultivariateRationalPolynomial& operator+= (const ratNum_t& r) {
			return (*this += RationalNumber(r));
		}

		/**
		 * Add RationalNumber r and SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator+ (const RationalNumber& r, const SparseMultivariateRationalPolynomial& b) {
			return (b + r);
		}
		
		/**
		 * Update *this by adding the RationalNumber r to it.
		 */ 
		SparseMultivariateRationalPolynomial& operator+= (const RationalNumber& r);

		/**
		 * Subtract the ratNum_t r from *this.
		 */
		inline SparseMultivariateRationalPolynomial operator- (const ratNum_t& r) const {
			return (*this - RationalNumber(r));
		}

		/**
		 * Subtract the RationalNumber r from *this.
		 */
		SparseMultivariateRationalPolynomial operator- (const RationalNumber& r) const;

		/**
		 * Subtract SMQP b from ratNum_t r.
		 */
		inline friend SparseMultivariateRationalPolynomial operator- (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) {
			return (-b + r);
		}

		/**
		 * Subtract SMQP b from ratNum_t r.
		 */
		inline friend SparseMultivariateRationalPolynomial operator- (const RationalNumber& r, const SparseMultivariateRationalPolynomial& b) {
			return (-b + r);
		}

		/**
		 * Update *this by subtracting ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial& operator-= (const ratNum_t& r) {
			return (*this -= RationalNumber(r));
		}

		/**
		 * Update *this by subtracting RationalNumber r.
		 */
		SparseMultivariateRationalPolynomial& operator-= (const RationalNumber& r);

		/**
		 * Multiply *this by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial operator* (const ratNum_t& r) const {
			return (*this * RationalNumber(r));
		}

		/**
		 * Multiply *this by RationalNumber r.
		 */
		SparseMultivariateRationalPolynomial operator* (const RationalNumber& r) const;
		
		/**
		 * Multiply ratNum_t r and SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator* (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) {
			return (b * r);
		}

		/**
		 * Update *this by multiplying by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial& operator*= (const ratNum_t& r) {
			return (*this *= RationalNumber(r));
		}

		/**
		 * Multiply RationalNumber r and SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator* (const RationalNumber& r, const SparseMultivariateRationalPolynomial& b) {
			return (b * r);
		}
		
		/**
		 * Update *this by multiplying by RationalNumber r.
		 */
		SparseMultivariateRationalPolynomial& operator*= (const RationalNumber& r);

		/**
		 * Divide *this by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial operator/ (const ratNum_t& r) const {
			return (*this / RationalNumber(r));
		}

		/**
		 * Divide *this by RationalNumber r.
		 */
		SparseMultivariateRationalPolynomial operator/ (const RationalNumber& r) const;

		/**
		 * Divide ratNum_t r by SMQP b.
		 */
		friend SparseMultivariateRationalPolynomial operator/ (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) ;

		/** 
		 * Update *this by dividing by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial& operator/= (const ratNum_t& r) {
			return (*this /= RationalNumber(r));
		}
		/**
		 * Divide RationalNumber r by SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator/ (const RationalNumber& r, const SparseMultivariateRationalPolynomial& b) {
			ratNum_t t;
			mpq_init(t);
			mpq_set(t, r.get_mpq_t());
			SparseMultivariateRationalPolynomial ret = (t / b);
			mpq_clear(t);
			return ret;
		}

		/** 
		 * Update *this by dividing by RationalNumber r.
		 */
		SparseMultivariateRationalPolynomial& operator/= (const RationalNumber& r); 

		/**
		 * Get the polynomial term at index. Returns 0 if index is beyond the 
		 * number of terms in this polynomial.
		 */
		SparseMultivariateRationalPolynomial operator[] (int index) const;

		/**
		 * Get the leading variable, that is, the highest-order variable with positive degree
		 * of this polynomial. 
		 * returns the leading variable or the empty string if this polynomial has zero variables.
		 */
		Symbol leadingVariable() const;

		/**
		 * Get the degree of this polynomial w.r.t the leading variable.
		 */
		Integer leadingVariableDegree() const;
    degree_t leadingVariableDegree_tmp() const;
		/**
		 * Is the contant term zero.
		 */
		bool isConstantTermZero() const;

		/**
		 * Get the leading coefficient w.r.t the input variable 'x'. 
		 * Returns the leading exponent as e.
		 *
		 * returns the coefficient as a new SMQP.
		 */
		SparseMultivariateRationalPolynomial leadingCoefficientInVariable (const Symbol& x, int* e = NULL) const;

		/**
		 * Convert to a SUP<SMQP> given the variable 'x'.
		 *
		 * returns the SUP<SMQP>.
		 */
		SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> convertToSUP(const Symbol& x) const;

		/**
		 * Negate all the coefficients of *this. Note, that due to the 
		 * sharing nature of underling nodes, this may alter the Nodes of
		 * other SMQP.
		 */
		void negate();

		/** 
		 * Get a copy of this such that all underlying memory is NOT shared.
		 * Note, any following arithmetic operations on the returned result
		 * will result in possibly making the underlying memory shared again. 
		 */
		SparseMultivariateRationalPolynomial deepCopy() const;

		/**
		 * Factors this polynomial into irreducible factors. 
		 * The Factors may include a common numerical (rational) factor.
		 * See also, Factors.
		 */
		Factors<SparseMultivariateRationalPolynomial> factor() const;

		/**
		 * SLP representation of the polynomial
		 **/
		void straightLineProgram();

		/**
		 * Print SLP representation
		 **/
		void printSLP(std::ostream& out = std::cout) const;

		/* Real Root Isolation */
	private:
		int isIntervalsMatchable(Intervals*, Intervals*, DenseUnivariateRationalPolynomial*, DenseUnivariateRationalPolynomial*, lfixq);
		int refineSleeveUnivariateIntervals(Intervals*, Intervals*, Intervals*, DenseUnivariateRationalPolynomial*, DenseUnivariateRationalPolynomial*, lfixq);
		void sleeveBoundURPolynomials(DenseUnivariateRationalPolynomial*, DenseUnivariateRationalPolynomial*, Intervals&, int, int);
	public: 
		/**
		 * Given one real root for x_1, .., x_{n-1},
		 * isolate positive roots for x_n
		 *
		 * @param mpIs: Roots of x_n (Output)
		 * @param apIs: A root of previous polynomials
		 * @param s: deal with 0: positive roots; 1: negative roots
		 * @param check: 1: check the leading or tail coefficient; 0: do not
		 * @param width: Interval's right - left < width
		 * @param ts: Taylor Shift option
		 *
		 * Return
		 * 1: Need to refine preious polynomials
		 * 0: Found positive real roots
		 **/		
		int positiveRealRootIsolation(Intervals* pIs, Intervals& apIs, mpq_class width, int ts=-1, bool s=0, bool check=1);

		/**
		 * Change *this to be a random non-zero polynomial.
		 * numvar: number of variables
		 * nterms: number of terms
		 * coefBound: limit on the number of bits encoding the coefficients.
		 * sparsity: succesive terms are at most sparsity away from each other
		 * includeNeg: a bool to say if coefs can be randomly negative or all positive
		 */
		void randomPolynomial(int numvar, int nterms, unsigned long int coefBound, degree_t sparsity, bool includeNeg);

		/**
		 * Change *this to be a random non-zero polynomial. The number of variables will
		 * be equal to the size of the maxDegs vector. The term whose monomial
		 * has exponents equal to those in maxDegs is guaranteed to exist in the
		 * resulting polynomial.
		 * A sparsity of 0 produces a dense polynomial. A sparsity of 1 produces only
		 * one term; one whose monomial has exponents equal to maxDegs.
		 * 
		 * coefBound: limit on the number of bits encoding the coefficients.
		 * sparsity: a percentage of zero-terms between term with maxDegs and the constant term.
		 * includeNeg: a bool to say if coefs can be randomly negative or all positive
		 */
		void randomPolynomial(std::vector<int> maxDegs, unsigned long int coefBound, float sparsity, bool includeNeg);

		/**
		 * Convert *this to an ExpressionTree.
		 *
		 * returns the new ExpressionTree.
		 */		
		ExpressionTree convertToExpressionTree() const;


/****************
* Multi-Divisor Division
*****************/	
    
/** 
 * Normal Form (or Multi-Divisor Division (MDD)) 
 * Given the dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}},
 * such that f = q_0*g_0 + ... + q_{s-1}*g_{s-1} + r.
 */
    SparseMultivariateRationalPolynomial lexNormalForm(const std::vector<Symbol>& superNames,
						       const std::vector<SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet = NULL) const;
    
  SparseMultivariateRationalPolynomial lexNormalizeDim0 (const std::vector<Symbol>& superNames, 
                                                         const std::vector<SparseMultivariateRationalPolynomial>& ts, SparseMultivariateRationalPolynomial* A) const;    

  /** 
 * Multi-Divisor Division (MDD) where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
    SparseMultivariateRationalPolynomial triangularSetNormalForm(const TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet) const;
	
/** 
 * Do the pseudo division of c by the triangular-set (divisor set) of B in the naive principle
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 */
    SparseMultivariateRationalPolynomial triangularSetPseudoDivide (const TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>& ts,
								    std::vector<SparseMultivariateRationalPolynomial>* quoSet, SparseMultivariateRationalPolynomial* h) const;
    
/** 
 * Specialized Normal Form where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s].
 */
SparseMultivariateRationalPolynomial triangularSetOnlyNormalForm (const TriangularSet<RationalNumber, SparseMultivariateRationalPolynomial>& ts) const;








/*

***/

Factors<SparseMultivariateRationalPolynomial> Factoring (SparseMultivariateRationalPolynomial& a);
Factors<SparseMultivariateRationalPolynomial> AlgebricFactorization(SparseMultivariateRationalPolynomial& minimal,SparseMultivariateRationalPolynomial*polyTotalcontent);


};





#endif //_SMQPLINKEDLIST_H_
