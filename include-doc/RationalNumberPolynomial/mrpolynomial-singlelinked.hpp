#ifndef _SMQPLINKEDLIST_H_
#define _SMQPLINKEDLIST_H_

#include "../polynomial.h"
#include "../Interval/interval.h"
#include "urpolynomial.h"
#include "../RingPolynomial/upolynomial.h"
#include "SMQP_CppSupport.hpp"
#include "SMQP_Support.h"
#include "SMQP_Support_Test.h"
#include "SMQP_Support_Recursive.h"
#include <gmpxx.h>
#include "../ExpressionTree/ExpressionTree.hpp"

/**
 * An element of the SLP of a polynomial.
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
 * Sparse Multivariate Rational Polynomial.
 * A multivariate polynomial with rational number coefficients using a 
 * sparse representation. 
 */
class SparseMultivariateRationalPolynomial: public BPASRecursivelyViewedPolynomial<RationalNumber,SparseMultivariateRationalPolynomial>,
											public BPASGCDDomain<SparseMultivariateRationalPolynomial> {

	private: 
		Node* poly; 		  //head of linked list
		mutable Node* tailNode; 		  //tail of linked list, not necessarily initialized. 
		int nvar;   		  //number of variables
		Symbol* names;   //list of strings representing the variables. 
							  //names[0] is 1 if "system-generated" variables
							  //names[0] is 9 if "user-specified" variables
		RationalNumber* constTerm; //if nvar == 0, this is the poly instead of poly.


		/**
		 * Construct an SMQP given the head Node*. The SMQP takes ownership
		 * of the Node list.
		 */
		SparseMultivariateRationalPolynomial(Node* head, int vars, Symbol* varNames);

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
		 * Returns a copy of *this under the new variable ordering supplied. 
		 * varmap is such that this.names[i] = newvars[varmap[i]]
		 * Returns an SMQP equal to *this but expended to newvars.
		 */
		SparseMultivariateRationalPolynomial expandVariables(int vars, Symbol* newvars, int varmap[]) const;

		/**
		 * Get the tail Node.
		 */
		Node* getTailNode() const;

		/** 
		 * Remove the tail Node. 
		 */
		void removeTailNode();

	public:
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

        inline SparseMultivariateRationalPolynomial unitCanonical(SparseMultivariateRationalPolynomial* u, SparseMultivariateRationalPolynomial* v) const {
        	RationalNumber lead = leadingCoefficient();
        	RationalNumber leadInv = lead.inverse();
        	if (u != NULL) {
        		*u = lead;
        	}
        	if (v != NULL) {
        		*v = leadInv;
        	}
        	return (*this * leadInv);
        }

        /**
         * Get GCD between *this and b.
         */
        SparseMultivariateRationalPolynomial gcd(const SparseMultivariateRationalPolynomial& b) const;

		
		/**
		 * Compute squarefree factorization of *this
		 */
		inline std::vector<SparseMultivariateRationalPolynomial> squareFree() const {
			std::cerr << "SparseMultivariateRationalPolynomial::squareFree NOT YET IMPLEMENTED" << std::endl;
			//TODO
			// When implementing this, remember to uncomment the corresponding test
			std::vector<SparseMultivariateRationalPolynomial> ret;
			ret.push_back(*this);
			return ret;
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

		/**
		 * Add the polynomails a and b and return the sum.
		 */
		friend SparseMultivariateRationalPolynomial operator+ (const SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

		/**
		 * Update *this by adding the specified polynomail to it.
		 */
		SparseMultivariateRationalPolynomial& operator+= (const SparseMultivariateRationalPolynomial& b);

		/**
		 * Subtract the specified polynomail from *this
		 */
		SparseMultivariateRationalPolynomial operator- (const SparseMultivariateRationalPolynomial& b) const;
		SparseMultivariateRationalPolynomial operator- (SparseMultivariateRationalPolynomial&& b) const;

		/**
		 * Subtract the polynomial b from a and return the difference.
		 */
		friend SparseMultivariateRationalPolynomial operator- (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

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

		/** 
		 * Multiply the polynomials a and b, returning their product.
		 */
		friend SparseMultivariateRationalPolynomial operator* (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

		/**
		 * Update this by multiplying by the specified polynomail.
		 */
		SparseMultivariateRationalPolynomial& operator*= (const SparseMultivariateRationalPolynomial& b);

		/**
		 * Divide *this by the specified polynomial.
		 */
		SparseMultivariateRationalPolynomial operator/ (const SparseMultivariateRationalPolynomial& b) const;
		SparseMultivariateRationalPolynomial operator/ (SparseMultivariateRationalPolynomial&& b) const;

		friend SparseMultivariateRationalPolynomial operator/ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b);

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

		RationalNumber content() const;

		SparseMultivariateRationalPolynomial primitivePart() const;

		SparseMultivariateRationalPolynomial initial() const;

		SparseMultivariateRationalPolynomial mainDegree() const;

		SparseMultivariateRationalPolynomial head() const;

		SparseMultivariateRationalPolynomial tail() const;

		SparseMultivariateRationalPolynomial separant() const; 

		/* BPASMultivariatePolynomial interface */

		/**
		 * Get the number of variables in this polynomial.
		 */
		int numberOfVariables() const;

		/**
		 * Get the number of variables in this polynomial with non-zero degree.
		 */
		int numberOfNonZeroVariables() const;

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
		RationalNumber coefficient(int, int*) const;
		
		/**
		 * Set a coefficient, given the exponent of each variable
		 */
		void setCoefficient(int, int*, const RationalNumber&);

		/**
		 * Set variables' names. It is invalid to try and set variable names
		 * with a vector whose size does not match this polynomail's nvar.
		 */
		void setVariableNames (const std::vector<Symbol>&);
		
		/**
		 * Get variable names of all variables available to this polynomial, 
		 * even those that have zero degree.
		 */
		std::vector<Symbol> variables() const;

		/**
		 * Get variable names of variables with non-zero degree;
		 */
		std::vector<Symbol> nonZeroVariables() const;


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
    	void differentiate(const Symbol& s, int k) {} //TODO

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
    	inline SparseMultivariateRationalPolynomial derivative(const Symbol& s, int k) const {
    	 	SparseMultivariateRationalPolynomial a(*this);
    	 	a.differentiate(s,k);
    	 	return a;
    	}
		
		/**
		 * Compute derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/ 
		inline SparseMultivariateRationalPolynomial derivative(const Symbol& s) const {
    	 	return this->derivative(s,0);
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
		 * Add *this and a ratNum_t.
		 */		
		inline SparseMultivariateRationalPolynomial operator+ (const ratNum_t& r) const {
			return (*this + mpq_class(r));
		}

		/**
		 * Add *this and the mpq_class r.
		 */
		SparseMultivariateRationalPolynomial operator+ (const mpq_class& r) const;

		inline SparseMultivariateRationalPolynomial operator+ (const RationalNumber& r) const {
			return (*this + r.get_mpq());
		}

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
			return (*this += mpq_class(r));
		}


		/**
		 * Add mpq_class r and SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator+ (const mpq_class& r, const SparseMultivariateRationalPolynomial& b) {
			return (b + r);
		}
		
		/**
		 * Update *this by adding the mpq_class r to it.
		 */ 
		SparseMultivariateRationalPolynomial& operator+= (const mpq_class& r);

		/**
		 * Update *this by adding the RationalNumber r to it.
		 */ 
		inline SparseMultivariateRationalPolynomial& operator+= (const RationalNumber& r) {
		    *this += r.get_mpq();
		    return *this;
		}

		/**
		 * Subtract the ratNum_t r from *this.
		 */
		inline SparseMultivariateRationalPolynomial operator- (const ratNum_t& r) const {
			return (*this - mpq_class(r));
		}

		/**
		 * Subtract the mpq_class r from *this.
		 */
		SparseMultivariateRationalPolynomial operator- (const mpq_class& r) const;

		/**
		 * Subtract the RationalNumber r from *this.
		 */
		inline SparseMultivariateRationalPolynomial operator- (const RationalNumber& r) const {
		    return (*this - r.get_mpq());
		}

		/**
		 * Subtract SMQP b from ratNum_t r.
		 */
		inline friend SparseMultivariateRationalPolynomial operator- (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) {
			return (-b + r);
		}

		/**
		 * Update *this by subtracting ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial& operator-= (const ratNum_t& r) {
			return (*this -= mpq_class(r));
		}


		/**
		 * Subtract SMQP b from mpq_class r.
		 */
		inline friend SparseMultivariateRationalPolynomial operator- (const mpq_class& r, const SparseMultivariateRationalPolynomial& b) {
			return (-b + r);
		}

		/**
		 * Update *this by subtracting mpq_class r.
		 */
		SparseMultivariateRationalPolynomial& operator-= (const mpq_class& r);

		/**
		 * Update *this by subtracting RationalNumber r.
		 */
		inline SparseMultivariateRationalPolynomial& operator-= (const RationalNumber& r) {
		    *this -= r.get_mpq();
		    return *this;
		}

		/**
		 * Multiply *this by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial operator* (const ratNum_t& r) const {
			return (*this * mpq_class(r));
		}

		/**
		 * Multiply *this by mpq_class r.
		 */
		SparseMultivariateRationalPolynomial operator* (const mpq_class& r) const;
		
		/**
		 * Multiply *this by RationalNumber r.
		 */
		inline SparseMultivariateRationalPolynomial operator* (const RationalNumber& r) const {
		    return (*this * r.get_mpq());
		}

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
			return (*this *= mpq_class(r));
		}

		/**
		 * Multiply mpq_class r and SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator* (const mpq_class& r, const SparseMultivariateRationalPolynomial& b) {
			return (b * r);
		}
		
		/**
		 * Update *this by multiplying by mpq_class r.
		 */
		SparseMultivariateRationalPolynomial& operator*= (const mpq_class& r);

		/**
		 * Update *this by multiplying by RationalNumber r.
		 */
		inline SparseMultivariateRationalPolynomial& operator*= (const RationalNumber& r) {
		    *this *= r.get_mpq();
		    return *this;
		}

		/**
		 * Divide *this by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial operator/ (const ratNum_t& r) const {
			return (*this / mpq_class(r));
		}

		/**
		 * Divide *this by mpq_class r.
		 */
		SparseMultivariateRationalPolynomial operator/ (const mpq_class& r) const;

		/**
		 * Divide *this by RationalNumber r.
		 */
		inline SparseMultivariateRationalPolynomial operator/ (const RationalNumber& r) const {
		    return (*this / r.get_mpq());
		}

		/**
		 * Divide ratNum_t r by SMQP b.
		 */
		friend SparseMultivariateRationalPolynomial operator/ (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) ;

		/** 
		 * Update *this by dividing by ratNum_t r.
		 */
		inline SparseMultivariateRationalPolynomial& operator/= (const ratNum_t& r) {
			return (*this /= mpq_class(r));
		}
		/**
		 * Divide mpq_class r by SMQP b.
		 */
		inline friend SparseMultivariateRationalPolynomial operator/ (const mpq_class& r, const SparseMultivariateRationalPolynomial& b) {
			ratNum_t t;
			mpq_init(t);
			mpq_set(t, r.get_mpq_t());
			SparseMultivariateRationalPolynomial ret = (t / b);
			mpq_clear(t);
			return ret;
		}

		/** 
		 * Update *this by dividing by mpq_class r.
		 */
		SparseMultivariateRationalPolynomial& operator/= (const mpq_class& r); 
		/** 
		 * Update *this by dividing by RationalNumber r.
		 */
		inline SparseMultivariateRationalPolynomial& operator/= (const RationalNumber& r) {
		    *this /= r.get_mpq();
		    return *this;
		}

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
		 * coefBound: upper (exclusive) limit of integer coef.
		 * sparsity: succesive terms are at most sparsity away from each other
		 * includeNeg: a bool to say if coefs can be neg or all positive
		 */
		void randomPolynomial(int numvar, int nterms, unsigned long int coefBound, degree_t sparsity, bool includeNeg);

		/**
		 * Convert *this to an ExpressionTree.
		 *
		 * returns the new ExpressionTree.
		 */		
		ExpressionTree convertToExpressionTree() const;

/****************
* SubResultantChain
*****************/

		

	/**
	 * Subresultant Chain
	 * Return the list of subresultants
	 *
	 * @param q: The other sparse univariate polynomial
	 **/
	std::vector<SparseMultivariateRationalPolynomial> subresultantChain (const SparseMultivariateRationalPolynomial& q) const;
	


	/**
	 * Subresultant Chain GCD
	 * Return the last non-zero subresultant of the current polynomial and the input polynomial if the resultant is zero and return 1 otherwise
	 *
	 * @param q: The other sparse univariate polynomial
	 **/
	SparseMultivariateRationalPolynomial subresultantGCD (const SparseMultivariateRationalPolynomial& q) const;



/****************
* Multi-Divisor Division
*****************/	

	/** Normal Form Algorithm in lexicographical polynomial ordering
	 * @param[in] superNames The vector of variable names in order
	 * @param[in] ts The triangular-set
	 * @param[out] r The remainder
	 * @param[out] quoSet The quotient-set
	 * \brief This algorithm runs by f and ts[s], and NULL quoSet[s] where s is the size of triangular-set,
	 * Computes r, quoSet[0], ... and quoSet[s] such that f = quoSet[0]*ts[0] + ... + quoSet[s-1]*ts[s-1] + r.
	 */
	SparseMultivariateRationalPolynomial lexNormalForm(const std::vector<Symbol>& superNames,
			const std::vector<SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet = NULL) const;

/*
		* Multi-Divisor Division (MDD)
		 * @param[in] dividend The dividend
		 * @param[in] divisorSet The divisor-set
		 * @param[in] t The type of MDD (type = 0 ? HeapMDD : (type = 1 ? triangularSetMDD : primitiveFactorTriangularSetMDD))
		 * @param[out] rem The remainder
		 * @param[out] quoSet The quotient-set
		 * @param[out] 0-1 Return 1, unless the C version of algorithm doesn't work
		 * \brief This algorithm computes the remainder and quoSet.
		friend bool multiDivisorDivision (const SparseMultivariateRationalPolynomial& dividend, const std::vector<SparseMultivariateRationalPolynomial>& divisorSet, std::vector<SparseMultivariateRationalPolynomial>* quoSet, SparseMultivariateRationalPolynomial* rem, int t);
*/
};

		



#endif //_SMQPLINKEDLIST_H_
