#ifndef _SMZPALTARRAY_H_
#define _SMZPALTARRAY_H_

#include "../Ring/Integer.hpp"
#include "../Polynomial/BPASRecursivePolynomial.hpp"
#include "uzpolynomial.h"
#include "../RingPolynomial/upolynomial.h"
#include "../RationalNumberPolynomial/mrpolynomial.h"
#include "SMZP_CppSupport.hpp"
#include "SMZP_Support.h"
#include "SMZP_Support_Test.h"
#include "SMZP_Support_Recursive.h"
#include <gmpxx.h>
#include "../ExpressionTree/ExpressionTree.hpp"
#include "../DataStructures/Factors.hpp"

#if defined(WITH_BLAD) && WITH_BLAD
#include "../BLADInterface/bladinterface.h"
#endif

#if defined(WITH_MAPLE) && WITH_MAPLE
#include "../MapleInterface/MapleInterfaceStream.hpp"
#endif

class SparseMultivariateRationalPolynomial;

/**
 * An element of the SLP of an integer polynomial.
 */
class SLPZRepresentation {

	public:
		union CoefOrInt {
			Integer* c;
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

		SLPZRepresentation() {
			a.c = NULL;
		}

		SLPZRepresentation(const SLPZRepresentation& r) {
			op = r.op;
			type = r.type;
			b = r.b;
			res = r.res;
			if (type == 0 || type == 1 || type == 4) {
				a.c = new Integer(*(r.a.c));
			} else {
				a.i = r.a.i;
			}
		}

		~SLPZRepresentation() {
			if (type == 0 || type == 1) {
				delete a.c;
			}
		}
};

/**
 * A multivariate polynomial with Integer coefficients using a
 * sparse representation. Only non-zero coefficients are encoded.
 */
class SparseMultivariateIntegerPolynomial: public BPASRecursivelyViewedPolynomial<Integer,SparseMultivariateIntegerPolynomial> {

	private:
		mutable AltArrZ_t* poly;
		int nvar;   		  //number of variables
		Symbol* names;   //list of strings representing the variables.
							  //names[0] is 1 if "system-generated" variables
							  //names[0] is 9 if "user-specified" variables

		friend class SparseMultivariateRationalPolynomial;
		friend class MapleInterface;

		/**
		 * Construct an SMZP given the alternating array.
		 */
		SparseMultivariateIntegerPolynomial(AltArrZ_t* aa, int vars, Symbol* varNames);

		/**
		 * Attempt to construct an SMZP given the rational number alternating array.
		 */
		SparseMultivariateIntegerPolynomial(const RationalNumber& r, AltArr_t* aa, int vars, Symbol* varNames);

		/**
		 * Determine if *this and b have compatible variable orderings.
		 * xs returns the mapping of both *this and b to the compatible superset
		 * returns true iff compatible.
		 */
		bool isOrderedRing(const SparseMultivariateIntegerPolynomial& b, std::vector<int>& xs) const;

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
		SparseMultivariateIntegerPolynomial expandVariables(int vars, Symbol* newvars, int varmap[]) const;

	public:

    SparseMultivariateIntegerPolynomial subresultantGCD (const SparseMultivariateIntegerPolynomial& q) const;

    std::vector<SparseMultivariateIntegerPolynomial> subresultantChain (const SparseMultivariateIntegerPolynomial& q, int filled = 0) const;

    SparseMultivariateIntegerPolynomial resultant (const SparseMultivariateIntegerPolynomial& q) const;

    std::vector<SLPZRepresentation> slp;

		/* Constructors */

		/**
		 * Construct a multivariate polynomial
		 *
		 **/
		SparseMultivariateIntegerPolynomial();

		/**
		 * Construct a multivariate polynomial with specific number of variables
		 *
		 * @param v: Number of variables
		 **/
		SparseMultivariateIntegerPolynomial(int v);

		/**htop

		 * Construct with a variable name such that f(x) = x;
		 *
		 * @param x: The variable name
		 **/
		SparseMultivariateIntegerPolynomial (const Symbol& x);

		/**
		 * Construct a polynomial by parsing the string str.
		 */
		SparseMultivariateIntegerPolynomial (const std::string& str);

		/**
		 * Copy Constructor.
		 *
		 * Does not reuse underlying memory allocated by b.
		 *
		 * @param b: A sparse multivariate polynomial
		 **/
		SparseMultivariateIntegerPolynomial(const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Move Constructor.
		 *
		 * @params b: The r-value reference polynomial.
		 */
		SparseMultivariateIntegerPolynomial(SparseMultivariateIntegerPolynomial&& b);

		/**
		 * Create a SMZP from an SMQP.
		 */
		SparseMultivariateIntegerPolynomial(const SparseMultivariateRationalPolynomial& b);

		/**
		 * Create a SMZP from an Integer.
		 */
		SparseMultivariateIntegerPolynomial(const Integer& r, int nvar = 0);

		/**
		 * Create a SMZP from a RationalNumber.
		 */
		SparseMultivariateIntegerPolynomial(const RationalNumber& r, int nvar = 0);

		/**
		 * Create a SMZP from a univariate rational polynomial.
		 *
		 * @param p: A SUZP polynomial.
		 **/
		SparseMultivariateIntegerPolynomial (const DenseUnivariateIntegerPolynomial& p);

		/**
		 * Construct from a SUP<SMZP> polynomial such that the resulting SMZP
		 * has main variable equal to the variable of s.
		 *
		 * @param s: The SUP<SMZP> polynomial
		 **/
    	SparseMultivariateIntegerPolynomial (const SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial>& s);

		/**
		 * Destroy the polynomial and underlying node memory.
		 **/
		~SparseMultivariateIntegerPolynomial();


		/* BPASRing interface */
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
        inline SparseMultivariateIntegerPolynomial unitCanonical(SparseMultivariateIntegerPolynomial* u, SparseMultivariateIntegerPolynomial* v) const {
        	if (isZero()) {
        		return *this;
        	}

        	if (mpz_cmp_si(this->poly->elems->coef, 0l) < 0) {
        		SparseMultivariateIntegerPolynomial ret = -(*this);
        		if (u != NULL) {
        			*u = SparseMultivariateIntegerPolynomial(Integer(-1));
        		}
        		if (v != NULL) {
        			*v = *u;
        		}
        		return ret;
        	}

        	if (u != NULL) {
        		(*u).one();
        	}
        	if (v != NULL) {
        		(*v).one();
        	}
        	return *this;
        }

        /* BPASPolynomial interface */

        /**
         * Assign this polynomail to equal the specified.
         */
        SparseMultivariateIntegerPolynomial& operator= (const SparseMultivariateIntegerPolynomial& b);

		/**
         * Movement assignment: move b to be this polynomail.
         */
        SparseMultivariateIntegerPolynomial& operator= (SparseMultivariateIntegerPolynomial&& b);

        /**
         * Assign this polynomial to be equal to the constant ring element.
         */
	    SparseMultivariateIntegerPolynomial& operator= (const Integer& r);

        /**
         * Add two SMQP polynomials together, *this and the specified.
         */
		SparseMultivariateIntegerPolynomial operator+ (const SparseMultivariateIntegerPolynomial& b) const;
		SparseMultivariateIntegerPolynomial operator+ (SparseMultivariateIntegerPolynomial&& b) const;

		/**
		 * Add the polynomails a and b and return the sum.
		 */
		// friend SparseMultivariateIntegerPolynomial operator+ (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Update *this by adding the specified polynomail to it.
		 */
		SparseMultivariateIntegerPolynomial& operator+= (const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Subtract the specified polynomail from *this
		 */
		SparseMultivariateIntegerPolynomial operator- (const SparseMultivariateIntegerPolynomial& b) const;
		SparseMultivariateIntegerPolynomial operator- (SparseMultivariateIntegerPolynomial&& b) const;

		/**
		 * Subtract the polynomial b from a and return the difference.
		 */
		// friend SparseMultivariateIntegerPolynomial operator- (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Unary operation, return  *this * -1.
		 */
		SparseMultivariateIntegerPolynomial operator- () const;

		/**
		 * Update *this by subtracting the specified polynomial from it.
		 */
		SparseMultivariateIntegerPolynomial& operator-= (const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Multiply *this by the specified polynomail.
		 */
		SparseMultivariateIntegerPolynomial operator* (const SparseMultivariateIntegerPolynomial& b) const;
		SparseMultivariateIntegerPolynomial operator* (SparseMultivariateIntegerPolynomial&& b) const;

		/**
		 * Multiply the polynomials a and b, returning their product.
		 */
		// friend SparseMultivariateIntegerPolynomial operator* (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Update this by multiplying by the specified polynomail.
		 */
		SparseMultivariateIntegerPolynomial& operator*= (const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Divide *this by the specified polynomial.
		 */
		SparseMultivariateIntegerPolynomial operator/ (const SparseMultivariateIntegerPolynomial& b) const;
		SparseMultivariateIntegerPolynomial operator/ (SparseMultivariateIntegerPolynomial&& b) const;

		// friend SparseMultivariateIntegerPolynomial operator/ (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Update *this by dividing by the specified polynomial.
		 */
		SparseMultivariateIntegerPolynomial& operator/= (const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Exponentiate *this by the input exponent integer.
		 * Treats negative exponents as positive.
		 */
		SparseMultivariateIntegerPolynomial operator^ (long long int e) const;

		/**
		 * Update *this by exponentiating this to the input integer.
		 * Treats negative exponents as positive.
		 */
		SparseMultivariateIntegerPolynomial& operator^= (long long int e);

		/**
		 * Determine if *this is equal to the specified polynomial.
		 * This takes into account the variable ordering on both poylnomials
		 * in such a way that the same polynomial under different variable orderings
		 * are NOT equal.
		 */
		bool operator== (const SparseMultivariateIntegerPolynomial& b) const;

		/**
		 * Determine if *this is not equal to the specified polynomial.
		 * This takes into account the variable ordering on both poylnomials
		 * in such a way that the same polynomial under different variable orderings
		 * are NOT equal.
		 */
		bool operator!= (const SparseMultivariateIntegerPolynomial& b) const;

		/**
		 * Output the string representation of *this to the input ostream.
		 */
		void print(std::ostream& os) const;

		/**
		 * Parse a polynomial from the in stream. Exactly one line is parsed.
		 */
		friend std::istream& operator>>(std::istream& in, SparseMultivariateIntegerPolynomial& p);

		/**
		 * Parse a polynomial from the string str and place it in *this.
		 */
		void fromString(const std::string& str);

        /**
         * Get GCD between *this and b.
         */
        SparseMultivariateIntegerPolynomial gcd(const SparseMultivariateIntegerPolynomial& b) const;

        /**
         * Get the GCD between *this and b as a primitive polynomial.
         */
        SparseMultivariateIntegerPolynomial primitiveGCD(const SparseMultivariateIntegerPolynomial& b) const;

		/**
		 * Compute squarefree factorization of *this with respect to all of its variables.
		 */
		Factors<SparseMultivariateIntegerPolynomial> squareFree() const;

		/**
		 * Compute squarefree factorization of *this with respect to the list of variables, vars.
		 */
		Factors<SparseMultivariateIntegerPolynomial> squareFree(const std::vector<Symbol>& vars) const;

		/**
		 * Computes the square free part of this polynomail. That is, the polynomial of this
		 * divided by all square factors. This is with respect to all variables.
		 */
		SparseMultivariateIntegerPolynomial squareFreePart() const;

		/**
		 * Computes the square free part of this polynomail. That is, the polynomial of this
		 * divided by all square factors. This is with respect to all variables.
		 */
		SparseMultivariateIntegerPolynomial squareFreePart(std::vector<Symbol>& vars) const;

		/**
		 * Get the content with respect to all variables. That is, a single
		 * rational number. The content here is one such that this/content is
		 * an integer polynomial with content of 1.
		 *
		 * Moreover, the content is one such that the leading coefficient
		 * of the corresponding primitive part is positive.
		 */
		Integer content() const;

		/**
		 * Get the content of this polynomial with respect to the variables in
		 * the input vector v.
		 *
		 * Moreover, the content is one such that the leading coefficient
		 * of the corresponding primitive part is positive.
		 */
		SparseMultivariateIntegerPolynomial content(const std::vector<Symbol>& v) const;

		/**
		 * Get the primitive part with respect to all variables.
		 * This is equivalent to this / content().
		 */
		SparseMultivariateIntegerPolynomial primitivePart() const;

		/**
		 * Get the primitive part with respect to all variables.
		 * This is equivalent to this / content().
		 *
		 * Simultaneously returns the rational number content in the parameter content.
		 */
		SparseMultivariateIntegerPolynomial primitivePart(Integer& content) const;

		/**
		 * Get the primitive part with respect to the variables in the vector v.
		 *
		 * returns the primitive part.
		 */
		SparseMultivariateIntegerPolynomial primitivePart(const std::vector<Symbol>& v) const;

		/**
		 * Get the primitive part with respect to the variables in the vector v.
		 * Returns the corresponding content in the content reference.
		 * returns the primitive part.
		 */
		SparseMultivariateIntegerPolynomial primitivePart(const std::vector<Symbol>& v, SparseMultivariateIntegerPolynomial& content) const;

		/**
		 * Get the leading coefficient of *this with respect to the main variable.
		 *
		 * returns the initial.
		 */
		SparseMultivariateIntegerPolynomial initial() const;

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
		SparseMultivariateIntegerPolynomial rank() const;

		/**
		 * Get the head of this polynomial. That is, the initial multiplied
		 * by the rank.
		 *
		 * returns the head.
		 */
		SparseMultivariateIntegerPolynomial head() const;

		/**
		 * Get the tail of this polynomial. That is, This - this.head().
		 *
		 * returns the tail.
		 */
		SparseMultivariateIntegerPolynomial tail() const;

		SparseMultivariateIntegerPolynomial separant() const;

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
		Integer leadingCoefficient() const;

		/**
		 * Set the leading coefficient.
		 */
		void setLeadingCoefficient(const Integer& it);

		/**
		 * Get the trailing coefficient.
		 */
		Integer trailingCoefficient() const;

		/**
		 * Get a coefficient, given the exponent of each variable
		 */
		Integer coefficient(int, const int*) const;

		inline Integer coefficient(const std::vector<int>& v) const {
			return coefficient(v.size(), v.data());
		}

		/**
		 * Set a coefficient, given the exponent of each variable
		 */
		void setCoefficient(int, const int*, const Integer&);

		inline void setCoefficient(const std::vector<int>& v, const Integer& r) {
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
		bool isEqual(const SparseMultivariateIntegerPolynomial& b) const;

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
    	SparseMultivariateIntegerPolynomial derivative(const Symbol& s, int k) const;

		/**
		 * Compute derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/
		inline SparseMultivariateIntegerPolynomial derivative(const Symbol& s) const {
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
		 * If the integral of this polynomial cannot be represented as an integer
		 * polynomial then an error occurs. See rationalIntegral().
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the k-th derivative, k > 0
		 **/
    	SparseMultivariateIntegerPolynomial integral(const Symbol& s, int k) const;

		/**
		 * Return k-th integral. If Symbol s is not a symbol contained in this
		 * polynomial then it becomes the main variable of the result and
		 * integration proceeds as expected.
		 *
		 * This method returns a rational number polynomial and therefore works
		 * regardless of the coefficients of this polynomial.
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the k-th derivative, k > 0
		 **/
		SparseMultivariateRationalPolynomial rationalIntegral(const Symbol& s, int k) const;

		/**
		 * Compute integral
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/
		inline SparseMultivariateIntegerPolynomial integral(const Symbol& s) const {
    	 	return this->integral(s,1);
        }

        /**
         * Return the integral of this, with respect to Symbol s, as a rational number polynomial.
         *
         * @param s: the symbol to integrate with respect to
         * returns the integral
         */
		SparseMultivariateRationalPolynomial rationalIntegral(const Symbol& s) const;

		/**
		 * Evaluate f(x)
		 *
		 * @param syms: Array of Symbols to evaluate at corresponding xs
		 * @param xs: Evaluation points
		 **/
    	inline SparseMultivariateIntegerPolynomial evaluate(int n, const Symbol* syms, const Integer* xs) const {
    		std::vector<Symbol> vecSyms;
    		std::vector<Integer> vecRats;
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
		SparseMultivariateIntegerPolynomial evaluate(const std::vector<Symbol>& vars, const std::vector<Integer>& values) const;

		/**
		 * Find the interpolating polynomial for the integer points and values specified by each vector, respectively.
		 * The points vector is a vector of vectors, where each element of the outer vector represents a multi-dimensional point.
		 * Each multi-dimensional point should have the same number of dimensions and in the same order.
		 *
		 * An error will occur if the resulting interpolating polynomial cannot be represented as an integer polynomial.
		 *
		 * returns the interpolating polynomial.
		 */
		static SparseMultivariateIntegerPolynomial interpolate(const std::vector<std::vector<Integer>>& points, const std::vector<Integer>& vals);

		/**
		 * Find the interpolating polynomial for the rational number points and values specified by each vector, respectively.
		 * The points vector is a vector of vectors, where each element of the outer vector represents a multi-dimensional point.
		 * Each multi-dimensional point should have the same number of dimensions and in the same order.
		 *
		 * An error will occur if the resulting interpolating polynomial cannot be represented as an integer polynomial.
		 *
		 * returns the interpolating polynomial.
		 */
		static SparseMultivariateIntegerPolynomial interpolate(const std::vector<std::vector<RationalNumber>>& points, const std::vector<RationalNumber>& vals);


		/**
		 * Divide this by polynomial b, returning the quotient and remainder in q and r, respectively.
		 *
		 * returns a boolean indicating if the division was exact.
		 */
		bool divide(const SparseMultivariateIntegerPolynomial& b, SparseMultivariateIntegerPolynomial& q, SparseMultivariateIntegerPolynomial& r) const;

		/**
		 * Divide this by polynomial b, return the quotient and remainder in q and r as rational number polynomials.
		 * This differs from the normal divide function where divisibility of coefficients does not matter as the result
		 * belongs to the polynomial ring over the rational numbers.
		 *
		 * returns a boolean indicating if the division was exact.
		 */
		bool rationalDivide(const SparseMultivariateIntegerPolynomial& b, SparseMultivariateRationalPolynomial& q, SparseMultivariateRationalPolynomial& r) const;

		/**
		 * Get the remainder of *this divided by b.
		 */
		SparseMultivariateIntegerPolynomial operator% (const SparseMultivariateIntegerPolynomial& b) const;

		/**
		 * Update *this by setting it to the remainder of *this divided by b.
		 */
		SparseMultivariateIntegerPolynomial& operator%= (const SparseMultivariateIntegerPolynomial& b);

		/**
		 * Pseudo divide this by b. The remainder is returned.
		 * if parameter quo is not null then the quotient is returned in quo.
		 * if parameter mult is not null then the multiplier is set to the initial of b
		 * raised to the power of degree(c, mvar(c)) - degree(b, mvar(c)) + 1.
		 *
		 * returns the pseudo remainder.
		 */
		SparseMultivariateIntegerPolynomial pseudoDivide(const SparseMultivariateIntegerPolynomial& b, SparseMultivariateIntegerPolynomial* quo = NULL, SparseMultivariateIntegerPolynomial* mult = NULL, bool lazy = 0) const;


		/**
		 * Add *this and a mpz_t.
		 */
		inline SparseMultivariateIntegerPolynomial operator+ (const mpz_t r) const {
			return (*this + mpz_class(r));
		}

		/**
		 * Add *this and the mpz_class r.
		 */
		SparseMultivariateIntegerPolynomial operator+ (const mpz_class& r) const;

		inline SparseMultivariateIntegerPolynomial operator+ (const Integer& r) const {
			return (*this + r.get_mpz());
		}

		/**
		 * Add mpz_t r and SMQP b.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator+ (const mpz_t r, const SparseMultivariateIntegerPolynomial& b) {
			return (b + r);
		}

		/**
		 * Update *this by adding r
		 */
		inline SparseMultivariateIntegerPolynomial& operator+= (const mpz_t r) {
			return (*this += mpz_class(r));
		}


		/**
		 * Add mpz_class r and SMQP b.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator+ (const mpz_class& r, const SparseMultivariateIntegerPolynomial& b) {
			return (b + r);
		}

		/**
		 * Update *this by adding the mpz_class r to it.
		 */
		SparseMultivariateIntegerPolynomial& operator+= (const mpz_class& r);

		/**
		 * Update *this by adding the Integer r to it.
		 */
		inline SparseMultivariateIntegerPolynomial& operator+= (const Integer& r) {
		    *this += r.get_mpz();
		    return *this;
		}

		/**
		 * Subtract the mpz_t r from *this.
		 */
		inline SparseMultivariateIntegerPolynomial operator- (const mpz_t r) const {
			return (*this - mpz_class(r));
		}

		/**
		 * Subtract the mpz_class r from *this.
		 */
		SparseMultivariateIntegerPolynomial operator- (const mpz_class& r) const;

		/**
		 * Subtract the Integer r from *this.
		 */
		inline SparseMultivariateIntegerPolynomial operator- (const Integer& r) const {
		    return (*this - r.get_mpz());
		}

		/**
		 * Subtract SMQP b from mpz_t r.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator- (const mpz_t r, const SparseMultivariateIntegerPolynomial& b) {
			return (-b + r);
		}

		/**
		 * Update *this by subtracting mpz_t r.
		 */
		inline SparseMultivariateIntegerPolynomial& operator-= (const mpz_t r) {
			return (*this -= mpz_class(r));
		}


		/**
		 * Subtract SMQP b from mpz_class r.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator- (const mpz_class& r, const SparseMultivariateIntegerPolynomial& b) {
			return (-b + r);
		}

		/**
		 * Update *this by subtracting mpz_class r.
		 */
		SparseMultivariateIntegerPolynomial& operator-= (const mpz_class& r);

		/**
		 * Update *this by subtracting Integer r.
		 */
		inline SparseMultivariateIntegerPolynomial& operator-= (const Integer& r) {
		    *this -= r.get_mpz();
		    return *this;
		}

		/**
		 * Multiply *this by mpz_t r.
		 */
		inline SparseMultivariateIntegerPolynomial operator* (const mpz_t r) const {
			return (*this * mpz_class(r));
		}

		/**
		 * Multiply *this by mpz_class r.
		 */
		SparseMultivariateIntegerPolynomial operator* (const mpz_class& r) const;

		/**
		 * Multiply *this by Integer r.
		 */
		inline SparseMultivariateIntegerPolynomial operator* (const Integer& r) const {
		    return (*this * r.get_mpz());
		}

		/**
		 * Multiply mpz_t r and SMQP b.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator* (const mpz_t r, const SparseMultivariateIntegerPolynomial& b) {
			return (b * r);
		}

		/**
		 * Update *this by multiplying by mpz_t r.
		 */
		inline SparseMultivariateIntegerPolynomial& operator*= (const mpz_t r) {
			return (*this *= mpz_class(r));
		}

		/**
		 * Multiply mpz_class r and SMQP b.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator* (const mpz_class& r, const SparseMultivariateIntegerPolynomial& b) {
			return (b * r);
		}

		/**
		 * Update *this by multiplying by mpz_class r.
		 */
		SparseMultivariateIntegerPolynomial& operator*= (const mpz_class& r);

		/**
		 * Update *this by multiplying by Integer r.
		 */
		inline SparseMultivariateIntegerPolynomial& operator*= (const Integer& r) {
		    *this *= r.get_mpz();
		    return *this;
		}

		/**
		 * Divide *this by mpz_t r.
		 */
		inline SparseMultivariateIntegerPolynomial operator/ (const mpz_t r) const {
			return (*this / mpz_class(r));
		}

		/**
		 * Divide *this by mpz_class r.
		 */
		SparseMultivariateIntegerPolynomial operator/ (const mpz_class& r) const;

		/**
		 * Divide *this by Integer r.
		 */
		inline SparseMultivariateIntegerPolynomial operator/ (const Integer& r) const {
		    return (*this / r.get_mpz());
		}

		/**
		 * Divide mpz_t r by SMQP b.
		 */
		friend SparseMultivariateIntegerPolynomial operator/ (const mpz_t r, const SparseMultivariateIntegerPolynomial& b) ;

		/**
		 * Update *this by dividing by mpz_t r.
		 */
		inline SparseMultivariateIntegerPolynomial& operator/= (const mpz_t r) {
			return (*this /= mpz_class(r));
		}
		/**
		 * Divide mpz_class r by SMQP b.
		 */
		inline friend SparseMultivariateIntegerPolynomial operator/ (const mpz_class& r, const SparseMultivariateIntegerPolynomial& b) {
			mpz_t t;
			mpz_init(t);
			mpz_set(t, r.get_mpz_t());
			SparseMultivariateIntegerPolynomial ret = (t / b);
			mpz_clear(t);
			return ret;
		}

		/**
		 * Update *this by dividing by mpz_class r.
		 */
		SparseMultivariateIntegerPolynomial& operator/= (const mpz_class& r);
		/**
		 * Update *this by dividing by Integer r.
		 */
		inline SparseMultivariateIntegerPolynomial& operator/= (const Integer& r) {
		    *this /= r.get_mpz();
		    return *this;
		}

		/**
		 * Get the polynomial term at index. Returns 0 if index is beyond the
		 * number of terms in this polynomial.
		 */
		SparseMultivariateIntegerPolynomial operator[] (int index) const;

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
		SparseMultivariateIntegerPolynomial leadingCoefficientInVariable (const Symbol& x, int* e = NULL) const;

		/**
		 * Convert to a SUP<SMQP> given the variable 'x'.
		 *
		 * returns the SUP<SMQP>.
		 */
		SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial> convertToSUP(const Symbol& x) const;

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
		SparseMultivariateIntegerPolynomial deepCopy() const;

		/**
		 * Factors this polynomial into irreducible factors.
		 * The Factors may include a common numerical (integer) factor.
		 * See also, Factors.
		 */
		Factors<SparseMultivariateIntegerPolynomial> factor() const;

		/**
		 * SLP representation of the polynomial
		 **/
		void straightLineProgram();

		/**
		 * Print SLP representation
		 **/
		void printSLP(std::ostream& out = std::cout) const;

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
};





#endif //_SMQPLINKEDLIST_H_
