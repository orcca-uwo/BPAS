#ifndef _UZPOLYNOMIAL_H_
#define _UZPOLYNOMIAL_H_


#include "../Polynomial/BPASUnivarPolynomial.hpp"
#include "modpoly.h"
#include "../ring.h"

/**
 * A univariate polynomial with Integer coefficients using a dense representation.
 * This representation stores all possible coefficients, up to a maximum degree,
 * even if they are zero.
 */
class DenseUnivariateIntegerPolynomial : public BPASUnivariatePolynomial<Integer,DenseUnivariateIntegerPolynomial> {
	private:
		Symbol name;	// Variable name
		int curd;		// Current degree
		int n;			// Maximum size of the polynomial
		mpz_class* coef;	// Coefficients

		inline void zeros() {
			for (int i = 0; i < n; ++i)
				coef[i] = 0;
		}

    	bool isEqual(const DenseUnivariateIntegerPolynomial& q) const;
    	void pomopo(const mpz_class c, const mpz_class t, const DenseUnivariateIntegerPolynomial& b);
    	void resetDegree();
		DenseUnivariateIntegerPolynomial euclideanGCD (const DenseUnivariateIntegerPolynomial& q) const;
		bool isDivide(DenseUnivariateIntegerPolynomial, DenseUnivariateIntegerPolynomial, const DenseUnivariateIntegerPolynomial&) const
		;
		DenseUnivariateIntegerPolynomial modularGCD (const DenseUnivariateIntegerPolynomial& q) const;

	public:
		static mpz_class characteristic;
		// static bool isPrimeField;
		// static bool isSmallPrimeField;
        // static bool isComplexField;
		/**
		 * Construct a polynomial
		 *
		 * @param d
		 **/
		DenseUnivariateIntegerPolynomial () : curd(0), n(1), name("%") {
			coef = new mpz_class[1];
			coef[0] = 0;
		}
		/**
		 * Construct a polynomial with degree
		 *
		 * @param d: Size of the polynomial
		 **/
		DenseUnivariateIntegerPolynomial(int s) {
			if (s < 1) { s = 1; }
			n = s;
			coef = new mpz_class[n];
			curd = 0;
			//coef[0] = 0;
			zeros();
			name = "%";
		}

		/**
		 * Construct a polynomial with a coeffient
		 * @param e: The coefficient
		 **/
		DenseUnivariateIntegerPolynomial (const Integer& e) : curd(0), n(1), name("%") {
			coef = new mpz_class[1];
			coef[0] = e.get_mpz();
		}
		DenseUnivariateIntegerPolynomial (const RationalNumber& e) : curd(0), n(1), name("%") {
			if (e.get_den() == 1) {
				coef = new mpz_class[1];
				coef[0] = e.get_num().get_mpz();
			}
			else {
				std::cout << "BPAS error, try to construct a rational number in DUZP." << std::endl;
				exit(1);
			}
		}
		/**
		 * Copy constructor
		 *
		 * @param b: A densed univariate rationl polynomial
		 **/
		DenseUnivariateIntegerPolynomial(const DenseUnivariateIntegerPolynomial& b) : curd(b.curd), name(b.name) {
			n = curd + 1;
			coef = new mpz_class[n];
			std::copy(b.coef, b.coef+n, coef);
		}
		/**
		 * Destroy the polynomial
		 *
		 * @param
		 **/
		~DenseUnivariateIntegerPolynomial() {
			delete [] coef;
		}

		/**
		 * Get degree of the polynomial
		 *
		 * @param
		 **/
		inline Integer degree() const {
			return curd;
		}

		/**
		 * Get the leading coefficient
		 *
		 * @param
		 **/
		inline Integer leadingCoefficient() const {
			return coef[curd];
		}

		inline Integer trailingCoefficient() const {
			for(int i = 0; i <= curd; ++i) {
				if (coef[i] != 0) {
					return coef[i];
				}
			}
			return 0;
		}

		inline Integer numberOfTerms() const {
			size_t c = 0;
			for (size_t i = 0; i <= curd; ++i) {
				++c;
			}
			return c;
		}

		/**
		 * Get coefficients of the polynomial, given start offset
		 *
		 * @param k: Offset
		 **/
		inline mpz_class* coefficients(int k=0) const {
#ifdef BPASDEBUG
			if (k < 0 || k >= n)
				std::cout << "BPAS: warning, try to access a non-exist coefficient " << k << " from DUZP(" << n << ")." << std::endl;
#endif
			return &coef[k];
		}
		/**
		 * Get a coefficient of the polynomial
		 *
		 * @param k: Offset
		 **/
		inline Integer coefficient(int k) const {
			if (k < 0 || k >= n) {
				mpz_class z(0);
				return Integer(z);
			}
			return Integer(coef[k]);
		}
		/**
		 * Set a coefficient of the polynomial
		 *
		 * @param k: Offset
		 * @param val: Coefficient
		 **/
		inline void setCoefficient(int k, const mpz_class value) {
			if (k >= n || k < 0) {
				std::cout << "BPAS: error, DUZP(" << n << ") but trying to access " << k << "." << std::endl;
				exit(1);
			}
			coef[k] = value;
			if (k > curd && value != 0)
				curd = k;
			resetDegree();
		}
		inline void setCoefficient(int k, const Integer& value) {
			setCoefficient(k, value.get_mpz());
		}
		inline void setCoefficient(int k, const int value) {
			setCoefficient(k, mpz_class(value));
		}


		/**
		 * Set a coefficient of the polynomial
		 *
		 * @param k: Degree of the term of which you are setting the coefficient
		 * @param val: Coefficient value
		 **/
		inline void setCoefficient(Integer k, const Integer& value);

		/**
		 * Get variable's name
		 *
		 * @param
		 **/
		inline Symbol variable() const {
			return name;
		}
		/**
		 * Set variable's name
		 *
		 * @param x: Varable's name
		 **/
		inline void setVariableName (const Symbol& x) {
			name = x;
		}
		/**
		 * Overload operator =
		 *
		 * @param b: A univariate integer polynoial
		 **/
		inline DenseUnivariateIntegerPolynomial& operator= (const DenseUnivariateIntegerPolynomial& b) {
			if (this != &b) {
				if (n) { delete [] coef; n = 0; }
				name = b.name;
				curd = b.curd;
				n = curd + 1;
				coef = new mpz_class[n];
				std::copy(b.coef, b.coef+n, coef);
			}
			return *this;
		}

		inline DenseUnivariateIntegerPolynomial& operator= (const Integer& i) {
			*this = DenseUnivariateIntegerPolynomial(i);
			return *this;
		}

		/**
		 * Overload operator !=
		 *
		 * @param b: A univariate integer polynoial
		 **/
		inline bool operator!= (const DenseUnivariateIntegerPolynomial& b) const {
			return !(isEqual(b));
		}

		/**
		 * Overload operator ==
		 *
		 * @param b: A univariate integer polynoial
		 **/
		inline bool operator== (const DenseUnivariateIntegerPolynomial& b) const {
			return isEqual(b);
		}

		/**
		 * Is zero polynomial
		 *
		 * @param
		 **/
		inline bool isZero () const {
			if (!curd)
				return (coef[0] == 0);
			return 0;
		}

		/**
		 * Zero polynomial
		 *
		 * @param
		 **/
		inline void zero() {
			curd = 0;
			zeros();
		}

		/**
		 * Is polynomial a constatn 1
		 *
		 * @param
		 **/
		inline bool isOne() const {
			if (!curd)
				return (coef[0] == 1);
			return 0;
		}

		/**
		 * Set polynomial to 1
		 *
		 * @param
		 **/
		inline void one() {
			curd = 0;
			coef[0] = 1;
			for (int i = 1; i < n; ++i)
				coef[i] = 0;
		}

		/**
		 * Is polynomial a constatn -1
		 *
		 * @param
		 **/
		inline bool isNegativeOne() const {
			if (!curd)
				return (coef[0] == -1);
			return 0;
		}

		/**
		 * Set polynomial to -1
		 *
		 * @param
		 **/
		inline void negativeOne() {
			curd = 0;
			coef[0] = -1;
			for (int i = 1; i < n; ++i)
				coef[i] = 0;
		}

		/**
		 * Is a constant
		 *
		 * @param
		 **/
		inline int isConstant() const {
			if (curd) { return 0; }
			else if (coef[0] >= 0) { return 1; }
			else { return -1; }
		}

		/**
	     * Obtain the unit normal (a.k.a canonical associate) of an element.
	     * If either parameters u, v, are non-NULL then the units are returned such that
	     * b = ua, v = u^-1. Where b is the unit normal of a, and is the returned value.
	     */
	    inline DenseUnivariateIntegerPolynomial unitCanonical(DenseUnivariateIntegerPolynomial* u = NULL, DenseUnivariateIntegerPolynomial* v = NULL) const {
	    	Integer lead = leadingCoefficient();
	    	Integer unit, uInv;
	    	lead.unitCanonical(&unit, &uInv);
	    	if (unit != 1) {
	    		if (u != NULL) {
	    			*u = unit;
	    		}
	    		if (v != NULL) {
	    			*v = uInv;
	    		}

	    		return *this * unit;

	    	}

	    	if (u != NULL) {
	    		*u = Integer(1);
	    	}
	    	if (v != NULL) {
	    		*v = Integer(1);
	    	}
	    	return *this;
	    }



		/**
		 * Content of the polynomial
		 *
		 * @param
		 **/
		inline Integer content() const {
			mpz_class c = coef[0];
			for (int i = 1; i <= curd; ++i) {
				if (coef[i] != 0) {
					mpz_gcd(c.get_mpz_t(), c.get_mpz_t(), coef[i].get_mpz_t());
					if (c == 1)
						break;
				}
			}
			return Integer(c);
		}

		inline DenseUnivariateIntegerPolynomial primitivePart() const {
			//TODO
			std::cerr << "BPAS ERROR: DUZP::primitivePart NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		/**
		 * Is the least signficant coefficient zero
		 *
		 * @param
		 **/
		inline bool isConstantTermZero() const {
			return (coef[0] == 0);
		}

		/**
		 * Overload operator ^
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
    	DenseUnivariateIntegerPolynomial operator^ (long long int e) const;

		/**
		 * Overload operator ^=
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		inline DenseUnivariateIntegerPolynomial& operator^= (long long int e) {
			*this = *this ^ e;
			return *this;
		}

		/**
		 * Overload operator <<
		 * replace by muplitying x^k
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
    	DenseUnivariateIntegerPolynomial operator<< (int k) const;

		/**
		 * Overload operator <<
		 * replace by muplitying x^k
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariateIntegerPolynomial& operator<<= (int k) {
			*this = *this << k;
			return *this;
		}

		/**
		 * Overload operator >>
		 * replace by dividing x^k, and
		 * return the quotient
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
    	DenseUnivariateIntegerPolynomial operator>> (int k) const;

		/**
		 * Overload operator >>=
		 * replace by dividing x^k, and
		 * return the quotient
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariateIntegerPolynomial& operator>>= (int k) {
			*this = *this >> k;
			return *this;
		}

		/**
		 * Overload operator +
		 *
		 * @param b: A univariate integer polynomial
		 **/
    	DenseUnivariateIntegerPolynomial operator+ (const DenseUnivariateIntegerPolynomial& b) const;

		/**
		 * Overload Operator +=
		 *
		 * @param b: A univariate integer polynomial
		 **/
		inline DenseUnivariateIntegerPolynomial& operator+= (const DenseUnivariateIntegerPolynomial& b) {
			if (curd >= b.curd)
				add(b);
			else
				*this = *this + b;
			return *this;
		}

		/**
		 * Add another polynomial to itself
		 *
		 * @param b: A univariate integer polynomial
		 **/
    	void add(const DenseUnivariateIntegerPolynomial& b);

		/**
		 * Overload Operator +
		 *
		 * @param c: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial operator+ (const Integer& c) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r += c);
		}

		inline DenseUnivariateIntegerPolynomial operator+ (const mpz_class& c) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r += c);
		}

		inline DenseUnivariateIntegerPolynomial operator+ (int c) const {
            DenseUnivariateIntegerPolynomial r (*this);
            return (r += c);
        }

		/**
		 * Overload Operator +=
		 *
		 * @param c: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial& operator+= (const Integer& c) {
			coef[0] += c.get_mpz();
			return *this;
		}

		inline DenseUnivariateIntegerPolynomial& operator+= (const mpz_class& c) {
			coef[0] += c;
			return *this;
		}

		inline DenseUnivariateIntegerPolynomial& operator+= (int c) {
            coef[0] += c;
            return *this;
        }

		inline friend DenseUnivariateIntegerPolynomial operator+ (const mpz_class& c, const DenseUnivariateIntegerPolynomial& p) {
			return (p + c);
		}

		inline friend DenseUnivariateIntegerPolynomial operator+ (int c, const DenseUnivariateIntegerPolynomial& p) {
            return (p + c);
        }

		/**
		 * Subtract another polynomial
		 *
		 * @param b: A univariate integer polynomial
		 */
    	DenseUnivariateIntegerPolynomial operator- (const DenseUnivariateIntegerPolynomial& b) const;

		/**
		 * Overload operator -=
		 *
		 * @param b: A univariate integer polynomial
		 **/
		inline DenseUnivariateIntegerPolynomial& operator-= (const DenseUnivariateIntegerPolynomial& b) {
			if (curd >= b.curd)
				subtract(b);
			else
				*this = *this - b;
			return *this;
		}

		/**
		 * Overload operator -, negate
		 *
		 * @param
		 **/
    	DenseUnivariateIntegerPolynomial operator- () const;

		/**
		 * Compute -f(x)
		 *
		 * @param
		 **/
		inline void negate() {
			for (int i = 0; i <= curd; ++i) {
				coef[i] = -coef[i];
			}
		}

		/**
		 * Subtract another polynomial from itself
		 *
		 * @param b: A univariate integer polynomial
		 **/
    	void subtract(const DenseUnivariateIntegerPolynomial& b);

		/**
		 *  Overload operator -
		 *
		 *  @param c: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial operator- (const Integer& c) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r -= c);
		}

		inline DenseUnivariateIntegerPolynomial operator- (const mpz_class& c) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r -= c);
		}

		inline DenseUnivariateIntegerPolynomial operator- (int c) const {
            DenseUnivariateIntegerPolynomial r (*this);
            return (r -= c);
        }

		/**
		 * Overload operator -=
		 *
		 * @param c: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial& operator-= (const Integer& c) {
			coef[0] -= c.get_mpz();
			return *this;
		}

		inline DenseUnivariateIntegerPolynomial& operator-= (const mpz_class& c) {
			coef[0] -= c;
			return *this;
		}

		inline DenseUnivariateIntegerPolynomial& operator-= (int c) {
            coef[0] -= c;
            return *this;
	    }

		inline friend DenseUnivariateIntegerPolynomial operator- (const mpz_class& c, const DenseUnivariateIntegerPolynomial& p) {
			return (-p + c);
		}

		inline friend DenseUnivariateIntegerPolynomial operator- (int c, const DenseUnivariateIntegerPolynomial& p) {
            return (-p + c);
        }

		/**
		 * Multiply to another polynomial
		 *
		 * @param b: A univariate integer polynomial
		 **/
    	DenseUnivariateIntegerPolynomial operator* (const DenseUnivariateIntegerPolynomial& b) const;

		/**
		 * Overload operator *=
		 *
		 * @param b: A univariate integer polynomial
		 **/
		inline DenseUnivariateIntegerPolynomial& operator*= (const DenseUnivariateIntegerPolynomial& b) {
			*this = *this * b;
			return *this;
		}

		/**
		 * Overload operator *
		 *
		 * @param e: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial operator* (const Integer& e) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r *= e);
		}

		inline DenseUnivariateIntegerPolynomial operator* (const mpz_class& e) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r *= e);
		}

		inline DenseUnivariateIntegerPolynomial operator* (int e) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r *= e);
		}

		/**
		 * Overload operator *=
		 *
		 * @param e: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial& operator*= (const Integer& e) {
			mpz_class c = e.get_mpz();
			*this *= c;
			return *this;
		}

		inline DenseUnivariateIntegerPolynomial& operator*= (const mpz_class& e) {
			if (e != 0 && e != 1) {
				for (int i = 0; i <= curd; ++i)
					coef[i] *= e;
			}
			else if (e == 0) { zero(); }
			return *this;
		}

		/**
		 * Overload operator *=
		 *
		 * @param e: A constant
		 **/
		inline DenseUnivariateIntegerPolynomial& operator*= (int e) {
			if (e != 0 && e != 1) {
				for (int i = 0; i <= curd; ++i)
					coef[i] *= e;
			}
			else if (e == 0) { zero(); }
			return *this;
		}

		inline friend DenseUnivariateIntegerPolynomial operator* (const mpz_class& e, const DenseUnivariateIntegerPolynomial& p) {
			return (p * e);
		}

		inline friend DenseUnivariateIntegerPolynomial operator* (int e, const DenseUnivariateIntegerPolynomial& p) {
            return (p * e);
        }

		/**
		 * Overload operator /
		 * ExactDivision
		 *
		 * @param b: A univariate integer polynomial
		 **/
		inline DenseUnivariateIntegerPolynomial operator/ (const DenseUnivariateIntegerPolynomial& b) const {
			DenseUnivariateIntegerPolynomial rem(*this);
			return (rem /= b);
		}

		/**
		 * Overload operator /=
		 * ExactDivision
		 *
		 * @param b: A univariate integer polynomial
		 **/
    	DenseUnivariateIntegerPolynomial& operator/= (const DenseUnivariateIntegerPolynomial& b);

		/**
		 * Overload operator /
		 *
		 * @param e: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial operator/ (const Integer& e) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r /= e);
		}

		inline DenseUnivariateIntegerPolynomial operator/ (const mpz_class& e) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r /= e);
		}

		inline DenseUnivariateIntegerPolynomial operator/ (int e) const {
			DenseUnivariateIntegerPolynomial r (*this);
			return (r /= e);
		}

		/**
		 * Overload operator /=
		 *
		 * @param e: An integer
		 **/
		inline DenseUnivariateIntegerPolynomial& operator/= (const Integer& e) {
			mpz_class c = e.get_mpz();
			return (*this /= c);
		}

    	DenseUnivariateIntegerPolynomial& operator/= (const mpz_class& e);

		inline DenseUnivariateIntegerPolynomial& operator/= (int e) {
			return (*this /= mpz_class(e));
		}

    	friend DenseUnivariateIntegerPolynomial operator/ (const mpz_class& e, const DenseUnivariateIntegerPolynomial& p);

		/**
		 * Monic division
		 * Return quotient and itself become the remainder
		 *
		 * @param b: The dividend polynomial
		 **/
    	DenseUnivariateIntegerPolynomial monicDivide(const DenseUnivariateIntegerPolynomial& b);

		/**
		 * Monic division
		 * Return quotient
		 *
		 * @param b: The dividend polynomial
		 * @param rem: The remainder polynomial
		 **/
		inline DenseUnivariateIntegerPolynomial monicDivide(const DenseUnivariateIntegerPolynomial& b, DenseUnivariateIntegerPolynomial* rem) const {
			*rem = *this;
			return rem->monicDivide(b);
		}

		/**
		 * Lazy pseudo dividsion
		 * Return the quotient and itself becomes remainder
		 * e is the exact number of division steps
		 *
		 * @param b: The dividend polynomial
		 * @param c: The leading coefficient of b to the power e
		 * @param d: That to the power deg(a) - deg(b) + 1 - e
		 **/
    	DenseUnivariateIntegerPolynomial lazyPseudoDivide (const DenseUnivariateIntegerPolynomial& b, Integer* c, Integer* d=NULL);

		/**
		 * Lazy pseudo dividsion
		 * Return the quotient
		 * e is the exact number of division steps
		 *
		 * @param b: The divident polynomial
		 * @param rem: The remainder polynomial
		 * @param c: The leading coefficient of b to the power e
		 * @param d: That to the power deg(a) - deg(b) + 1 - e
		 **/
		inline DenseUnivariateIntegerPolynomial lazyPseudoDivide (const DenseUnivariateIntegerPolynomial& b, DenseUnivariateIntegerPolynomial* rem, Integer* c, Integer* d) const {
			*rem = *this;
			return rem->lazyPseudoDivide(b, c, d);
		}

		/**
		 * Pseudo dividsion
		 * Return the quotient and itself becomes remainder
		 *
		 * @param b: The divident polynomial
		 * @param d: The leading coefficient of b
		 *           to the power deg(a) - deg(b) + 1
		 **/
    	DenseUnivariateIntegerPolynomial pseudoDivide (const DenseUnivariateIntegerPolynomial& b, Integer* d=NULL);

		/**
		 * Pseudo dividsion
		 * Return the quotient
		 *
		 * @param b: The divident polynomial
		 * @param rem: The remainder polynomial
		 * @param d: The leading coefficient of b
		 *           to the power deg(a) - deg(b) + 1
		 **/
    	DenseUnivariateIntegerPolynomial pseudoDivide (const DenseUnivariateIntegerPolynomial& b, DenseUnivariateIntegerPolynomial* rem, Integer* d) const;

		/**
		 * GCD(p, q)
		 *
		 * @param q: The other polynomial
		 **/
    	DenseUnivariateIntegerPolynomial gcd (const DenseUnivariateIntegerPolynomial& q, int type) const;

    	inline DenseUnivariateIntegerPolynomial gcd (const DenseUnivariateIntegerPolynomial& q) const {
    		return gcd(q, 0);
    	}

		/**
		 * Convert current object to its k-th derivative
		 *
		 * @param k: Order of the derivative, k > 0
		 **/
    	void differentiate(int k);

		/**
		 * Convert current object to its derivative
		 *
		 **/
    	inline void differentiate() {
    		this->differentiate(1);
    	}

		/**
		 * Return k-th derivative
		 *
		 * @param k: Order of the k-th derivative, k > 0
		 **/
    	 inline DenseUnivariateIntegerPolynomial derivative(int k) const {
    	 	DenseUnivariateIntegerPolynomial a(*this);
    	 	a.differentiate(k);
    	 	return a;
    	 }

		/**
		 * Compute derivative
		 *
		 **/
    	inline DenseUnivariateIntegerPolynomial derivative() const {
    	 	return this->derivative(1);
    	 }

		/**
		 * Convert current object to its integral with constant of integration 0
		 *
		 **/
    	void integrate();

		/**
		 * Compute integral with constant of integration 0
		 *
		 **/
    	inline DenseUnivariateIntegerPolynomial integral() {
    	 	DenseUnivariateIntegerPolynomial a(*this);
    	 	a.integrate();
    	 	return a;
    	 }

		/**
		/**
		 * Evaluate f(x)
		 *
		 * @param x: Evaluation point
		 **/
		 // THIS FUNCTION IS DEPRECATED
    	// mpz_class evaluate(const mpz_class& x) const;

		/**
		 * Evaluate f(x)
		 *
		 * @param x: Evaluation point
		 **/
    	Integer evaluate(const Integer& x) const;

		/**
		 * Square free
		 *
		 * @param
		 **/
    	Factors<DenseUnivariateIntegerPolynomial> squareFree() const;

		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param b: A univariate Integer polynoial
		 **/
    	void print(std::ostream &out) const;

    	/**
    	 * Convert this DUZP into an expression tree.
    	 */
    	ExpressionTree convertToExpressionTree() const;
};

#endif
