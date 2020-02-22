#ifndef _UNIPOLYNOMIAL_H_
#define _UNIPOLYNOMIAL_H_


#include "../DyadicRationalNumber/globals.h"
#include "../DyadicRationalNumber/Multiplication/multiplication.h"	// Taylor Shift DnC
#include "../Interval/interval.h"


void ts_modulo (lfixz*, lfixz, lfixz, int);

/**
 * A univariate polynomial with RationalNumber coefficients represented densely.
 */
class DenseUnivariateRationalPolynomial : public BPASUnivariatePolynomial<RationalNumber>,
										  public BPASEuclideanDomain{
	private:
		Symbol name;	// Variable name
		int curd;	// Current degree
		int n;		// Maximum size of the polynomial
		lfixq* coef;	// Coefficients

		inline void zeros() {
			for (int i = 0; i < n; ++i)
				coef[i] = 0;
		}
		bool isEqual(const DenseUnivariateRationalPolynomial&) const;

		void taylorShiftCVL();

		/* Taylor Shift DnC */
		/* Compute degree to the power of 2 */
		void taylorShiftDnC(int, int);

		/* Taylor Shift Dnc */
		/* Compute coefficients in each (1+x)^{2^i}, i = 0..log[2](n) */
		void taylorShiftDnC(int);
		void binomials(lfixz*, int);
		void taylorShiftBasePower2(lfixz*, int);

		/* Taylor Shift IncrementalCilkFor */
		void taylorShiftIncrementalCilkFor(lfixq*, int, int);
		void tableauBase(lfixq*, int);
		void polygonBase(lfixq*, int, int, int);
		void taylorShiftBase(lfixq*, int);
		void taylorShift3RegularCilkFor(lfixq*, int, int);

		/* Taylor Shift Tableau */
		void taylorShiftTableau(int);
		void taylorShiftBase(lfixq*, lfixq*, int);
		void tableauBaseInplace(lfixq*, lfixq*, int, int);
		void tableauConstruction(lfixq*, lfixq*, int, int, int);
		void taylorShiftGeneral(lfixq*, lfixq*, int, int);

		/* Root Bound */
		lfixz positiveRootBound();
		lfixz cauchyRootBound();
		lfixz complexRootBound();
		int rootBoundExp();

		/* Real Root Isolation */
		long taylorConstant(int, lfixz, int, lfixz);
		void genDescartes(Intervals* pIs, DenseUnivariateRationalPolynomial*, int);
		void isolateScaledUnivariatePolynomial(Intervals*, DenseUnivariateRationalPolynomial*, int);
		void isolatePositiveUnivariateRealRoots(Intervals*, DenseUnivariateRationalPolynomial*, int);
		void isolateUnivariateRealRoots(Intervals*, DenseUnivariateRationalPolynomial*, int);
		void refineUnivariateInterval(Interval*, lfixq, DenseUnivariateRationalPolynomial*, lfixq);
		void refineUnivariateIntervals(Intervals*, Intervals*, DenseUnivariateRationalPolynomial*, lfixq);
		void univariatePositiveRealRootIsolation(Intervals*, DenseUnivariateRationalPolynomial*, lfixq, int);
		void univariateRealRootIsolation(Intervals*, DenseUnivariateRationalPolynomial*, lfixq, int);


	    void pomopo(const lfixq c, const lfixq t, const DenseUnivariateRationalPolynomial& b);
    	void resetDegree();

		// gcd subroutines
		DenseUnivariateRationalPolynomial euclideanGCD (const DenseUnivariateRationalPolynomial& q) const;
		DenseUnivariateRationalPolynomial modularGCD (const DenseUnivariateRationalPolynomial& q) const;

	public:
		static mpz_class characteristic;
		static RingProperties properties;
		// static bool isPrimeField;
		// static bool isSmallPrimeField;
        // static bool isComplexField;
		/**
		 * Construct a polynomial
		 *
		 * @param d
		 **/ 
		DenseUnivariateRationalPolynomial () : curd(0), n(1), name("%") {
			coef = new lfixq[1];
			coef[0] = 0;
		}
		
		/**
		 * Construct a polynomial with degree
		 *
		 * @param d: Size of the polynomial
		 **/
		DenseUnivariateRationalPolynomial(int s) {
			if (s < 1) { s = 1; }
			n = s;
			coef = new lfixq[n];
			curd = 0;
			//coef[0] = 0;
			zeros();
			name = "%";
		}

		/**
		 * Construct a polynomial with a coefficient
		 *
		 * @param e: The coefficient
		 **/
		DenseUnivariateRationalPolynomial (const Integer& e) : curd(0), n(1), name("%")  {
			coef = new lfixq[1];
			coef[0] = mpq_class(e.get_mpz());
		}
		
		DenseUnivariateRationalPolynomial (const RationalNumber& e) : curd(0), n(1), name("%")  {
			coef = new lfixq[1];
			coef[0] = mpq_class(e.get_mpq());
		}
		
		/**
		 * Copy constructor
		 *
		 * @param b: A densed univariate rationl polynomial
		 **/ 
		DenseUnivariateRationalPolynomial(const DenseUnivariateRationalPolynomial& b) : curd(b.curd), name(b.name) {
			n = curd + 1;
			coef = new lfixq[n];
			std::copy(b.coef, b.coef+n, coef);
		}
		
		/**
		 * Destroy the polynomial
		 *
		 * @param
		 **/
		~DenseUnivariateRationalPolynomial() {
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
		inline RationalNumber leadingCoefficient() const {
			return coef[curd];
		}

		inline RationalNumber trailingCoefficient() const {
			for (size_t i = 0; i <= curd; ++i) {
				if (coef[i] != 0) {
					return coef[i];
				}
			}
			return 0;
		}

		inline Integer numberOfTerms() const {
			size_t c = 0;
			for (size_t i = 0; i <= curd; ++i) {
				if (coef[i] != 0){ 
					++c;
				}
			}
			return c;
		}
		
		/**
		 * Get coefficients of the polynomial, given start offset
		 *
		 * @param k: Offset
		 **/
		inline mpq_class* coefficients(int k=0) const {
#ifdef BPASDEBUG
			if (k < 0 || k >= n)
				std::cout << "BPAS: warning, try to access a non-exist coefficient " << k << " from DUQP(" << n << ")." << std::endl;
#endif
			return &coef[k];
		}
		
		/**
		 * Get a coefficient of the polynomial
		 *
		 * @param k: Offset
		 **/
		inline RationalNumber coefficient(int k) const {
			if (k < 0 || k >= n)
				return lfixq(0);
			return coef[k];
		}
		
		/**
		 * Set a coefficient of the polynomial
		 *
		 * @param k: Offset
		 * @param val: Coefficient
		 **/
		inline void setCoefficient(int k, const RationalNumber& value) {
			if (k >= n || k < 0) {
				std::cout << "BPAS: error, DUQP(" << n << ") but trying to access " << k << "." << std::endl;
				exit(1);
			}
			coef[k] = value.get_mpq();
			if (k > curd && value != 0)
				curd = k;
			resetDegree();
		}
		
		inline void setCoefficient(int k, double value) {
			setCoefficient(k, lfixq(value));
		}
		
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
		
		inline DenseUnivariateRationalPolynomial unitCanonical(DenseUnivariateRationalPolynomial* u = NULL, DenseUnivariateRationalPolynomial* v = NULL) const {
			RationalNumber lead = leadingCoefficient();
			RationalNumber leadInv = lead.inverse();
			DenseUnivariateRationalPolynomial ret = *this * leadInv;
			if (u != NULL) {
				*u = lead;
			}
			if (v != NULL) {
				*v = leadInv;
			}
			return ret;
		}

		/**
		 * Overload operator =
		 *
		 * @param b: A univariate rational polynoial
		 **/
		inline DenseUnivariateRationalPolynomial& operator= (const DenseUnivariateRationalPolynomial& b) {
			if (this != &b) {
				if (n) { delete [] coef; n = 0; }
				name = b.name;
				curd = b.curd;
				n = curd + 1;
				coef = new lfixq[n];
				std::copy(b.coef, b.coef+n, coef);
			}
			return *this;
		}

		inline DenseUnivariateRationalPolynomial& operator= (const RationalNumber& r) {
			*this = DenseUnivariateRationalPolynomial(r);
			return *this;
		}
		
		/**
		 * Overload operator !=
		 *
		 * @param b: A univariate rational polynoial
		 **/
		inline bool operator!= (const DenseUnivariateRationalPolynomial& b) const {
			return !(isEqual(b));
		}
		
		/**
		 * Overload operator ==
		 *
		 * @param b: A univariate rational polynoial
		 **/ 
		inline bool operator== (const DenseUnivariateRationalPolynomial& b) const {
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
			//coef[0] = 0;
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
		 * Content of the polynomial
		 *
		 * @param
		 **/ 
		inline RationalNumber content() const {
			return !isZero();
		}

		inline DenseUnivariateRationalPolynomial primitivePart() const {
			//TODO
			std::cerr << "BPAS ERROR: DUQP::primitivePart NOT YET IMPLEMENTED" << std::endl;
			return (*this);
		}
		
		/**
		 * Overload operator ^
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		DenseUnivariateRationalPolynomial operator^ (long long int e) const;
		
		/**
		 * Overload operator ^=
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		inline DenseUnivariateRationalPolynomial& operator^= (long long int e) {
			*this = *this ^ e;
			return *this;
		}
		
		/**
		 * Overload operator <<
		 * replace by muplitying x^k
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		DenseUnivariateRationalPolynomial operator<< (int k) const;
		
		/**
		 * Overload operator <<=
		 * replace by muplitying x^k
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariateRationalPolynomial& operator<<= (int k) {
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
		DenseUnivariateRationalPolynomial operator>> (int k) const;
		
		/**
		 * Overload operator >>=
		 * replace by dividing x^k, and
		 * return the quotient
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariateRationalPolynomial& operator>>= (int k) {
			*this = *this >> k;
			return *this;
		}
		
		/**
		 * Overload operator +
		 *
		 * @param b: A univariate rational polynomial
		 **/
		DenseUnivariateRationalPolynomial operator+ (const DenseUnivariateRationalPolynomial& b) const;
		
		/**
		 * Overload Operator +=
		 *
		 * @param b: A univariate rational polynomial
		 **/ 
		inline DenseUnivariateRationalPolynomial& operator+= (const DenseUnivariateRationalPolynomial& b) {
			if (curd >= b.curd)
				add(b);
			else
				*this = *this + b;
			return *this;
		}
		
		/**
		 * Add another polynomial to itself
		 *
		 * @param b: A univariate rational polynomial
		 **/ 
    	void add(const DenseUnivariateRationalPolynomial& b);
		
		/**
		 * Overload Operator +
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariateRationalPolynomial operator+ (const RationalNumber& c) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r += c);
		}
		
		inline DenseUnivariateRationalPolynomial operator+ (const mpq_class& c) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r += c);
		}
		
		/**
		 * Overload Operator +=
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariateRationalPolynomial& operator+= (const RationalNumber& c) {
			coef[0] += lfixq(c.get_mpq());
			return *this;
		}
		
		inline DenseUnivariateRationalPolynomial& operator+= (const mpq_class& c) {
			coef[0] += c;
			return *this;
		}

		inline friend DenseUnivariateRationalPolynomial operator+ (const mpq_class& c, const DenseUnivariateRationalPolynomial& p) {
			return (p + c);
		}
		
		/**
		 * Subtract another polynomial
		 *
		 * @param b: A univariate rational polynomial
		 */
    	DenseUnivariateRationalPolynomial operator- (const DenseUnivariateRationalPolynomial& b) const;
		
		/**
		 * Overload operator -=
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariateRationalPolynomial& operator-= (const DenseUnivariateRationalPolynomial& b) {
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
    	DenseUnivariateRationalPolynomial operator- () const;

		/**
		 * Subtract another polynomial from itself
		 *
		 * @param b: A univariate rational polynomial
		 **/ 
    	void subtract(const DenseUnivariateRationalPolynomial& b);

		/**
		 * Overload operator -
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariateRationalPolynomial operator- (const RationalNumber& c) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r -= c);
		}

		inline DenseUnivariateRationalPolynomial operator- (const mpq_class& c) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r -= c);
		}

		/**
		 * Overload operator -=
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariateRationalPolynomial& operator-= (const RationalNumber& c) {
			coef[0] -= lfixq(c.get_mpq());
			return *this;
		}
		
		inline DenseUnivariateRationalPolynomial& operator-= (const mpq_class& c) {
			coef[0] -= c;
			return *this;
		}

		inline friend DenseUnivariateRationalPolynomial operator- (const mpq_class& c, const DenseUnivariateRationalPolynomial& p) {
            return (-p + c);
        }

		/**
		 * Multiply to another polynomial
		 *
		 * @param b: A univariate rational polynomial
		 **/
    	DenseUnivariateRationalPolynomial operator* (const DenseUnivariateRationalPolynomial& b) const;

		/**
		 * Overload operator *=
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariateRationalPolynomial& operator*= (const DenseUnivariateRationalPolynomial& b) {
			*this = *this * b;
			return *this;
		}

		/**
		 * Overload operator *
		 *
		 * @param e: A rational number
		 **/
		inline DenseUnivariateRationalPolynomial operator* (const RationalNumber& e) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r *= e);
		}

		inline DenseUnivariateRationalPolynomial operator* (const mpq_class& e) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r *= e);
		}

		inline DenseUnivariateRationalPolynomial operator* (const sfixn& e) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r *= e);
		}
		
		/**
		 * Overload operator *=
		 *
		 * @param e: A rational number
		 **/
    	DenseUnivariateRationalPolynomial& operator*= (const RationalNumber& e);

    	DenseUnivariateRationalPolynomial& operator*= (const mpq_class& e);

		/**
		 * Overload operator *=
		 *
		 * @param e: A constant
		 **/
    	DenseUnivariateRationalPolynomial& operator*= (const sfixn& e);

		inline friend DenseUnivariateRationalPolynomial operator* (const mpq_class& c, const DenseUnivariateRationalPolynomial& p) {
            return (p * c);
        }

		inline friend DenseUnivariateRationalPolynomial operator* (const sfixn& c, const DenseUnivariateRationalPolynomial& p) {
			return (p * c);
		}

		/**
		 * Overload operator /
		 * ExactDivision
		 *
		 * @param b: A univariate rational polynomial
		 **/ 
		inline DenseUnivariateRationalPolynomial operator/ (const DenseUnivariateRationalPolynomial& b) const{
			DenseUnivariateRationalPolynomial rem(*this);
			return (rem /= b);
		}

		/**
		 * Overload operator /=
		 * ExactDivision
		 *
		 * @param b: A univariate rational polynomial
		 **/
    	DenseUnivariateRationalPolynomial& operator/= (const DenseUnivariateRationalPolynomial& b);

    	inline DenseUnivariateRationalPolynomial operator% (const DenseUnivariateRationalPolynomial& b) const {
    		DenseUnivariateRationalPolynomial ret(*this);
    		ret %= b;
    		return ret;
    	}

    	inline DenseUnivariateRationalPolynomial& operator%= (const DenseUnivariateRationalPolynomial& b) {
    		*this = this->remainder(b);
    		return *this;
    	}

		/**
		 * Overload operator /
		 *
		 * @param e: A rational number
		 **/
		inline DenseUnivariateRationalPolynomial operator/ (const RationalNumber& e) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r /= e);
		}

		inline DenseUnivariateRationalPolynomial operator/ (const mpq_class& e) const {
			DenseUnivariateRationalPolynomial r (*this);
			return (r /= e);
		}

		/**
		 * Overload operator /=
		 *
		 * @param e: A rational number
		 **/
    	DenseUnivariateRationalPolynomial& operator/= (const RationalNumber& e);

    	DenseUnivariateRationalPolynomial& operator/= (const mpq_class& e);

    	friend DenseUnivariateRationalPolynomial operator/ (const mpq_class& c, const DenseUnivariateRationalPolynomial& p);

		/**
		 * Monic division
		 * Return quotient and itself become the remainder
		 *
		 * @param b: The dividend polynomial
		 **/ 
    	DenseUnivariateRationalPolynomial monicDivide(const DenseUnivariateRationalPolynomial& b);

		/**
		 * Monic division
		 * Return quotient
		 *
		 * @param b: The dividend polynomial
		 * @param rem: The remainder polynomial
		 **/
    	DenseUnivariateRationalPolynomial monicDivide(const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* rem) const;

		/**
		 * Lazy pseudo dividsion
		 * Return the quotient and itself becomes remainder
		 * e is the exact number of division steps
		 *
		 * @param b: The dividend polynomial
		 * @param c: The leading coefficient of b to the power e
		 * @param d: That to the power deg(a) - deg(b) + 1 - e
		 **/
    	DenseUnivariateRationalPolynomial lazyPseudoDivide (const DenseUnivariateRationalPolynomial& b, RationalNumber* c, RationalNumber* d=NULL);

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
    	DenseUnivariateRationalPolynomial lazyPseudoDivide (const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* rem, RationalNumber* c, RationalNumber* d) const;

		/**
		 * Pseudo dividsion
		 * Return the quotient and itself becomes remainder
		 *
		 * @param b: The divident polynomial
		 * @param d: The leading coefficient of b
		 *           to the power deg(a) - deg(b) + 1
		 **/
    	DenseUnivariateRationalPolynomial pseudoDivide (const DenseUnivariateRationalPolynomial& b, RationalNumber* d=NULL);

		/**
		 * Pseudo dividsion
		 * Return the quotient
		 *
		 * @param b: The divident polynomial
		 * @param rem: The remainder polynomial
		 * @param d: The leading coefficient of b
		 *           to the power deg(a) - deg(b) + 1
		 **/ 
    	DenseUnivariateRationalPolynomial pseudoDivide (const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* rem, RationalNumber* d) const;

		/**
		 * s * a \equiv g (mod b), where g = gcd(a, b)
		 *
		 * @param b: A univariate polynomial
		 * @param g: The GCD of a and b
		 **/
    	DenseUnivariateRationalPolynomial halfExtendedEuclidean (const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* g) const;

		/**
		 * s*a + t*b = c, where c in the ideal (a,b)
		 *
		 * @param a: A univariate polynomial
		 * @oaran b: A univariate polynomial
		 * @param s: Either s = 0 or degree(s) < degree(b)
		 * @param t
		 **/
    	void diophantinEquationSolve(const DenseUnivariateRationalPolynomial& a, const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* s, DenseUnivariateRationalPolynomial* t) const;

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
		 * @param k: k-th derivative, k > 0
		 **/ 
    	inline DenseUnivariateRationalPolynomial derivative(int k) const {
    	 	DenseUnivariateRationalPolynomial a(*this);
    	 	a.differentiate(k);
    	 	return a;
    	}
		
		/**
		 * Compute derivative
		 *
		 **/ 
    	inline DenseUnivariateRationalPolynomial derivative() const {
    	 	return this->derivative(1);
    	}
		
		/**
		 * Compute the integral with constant of integration 0
		 *
		 * @param
		 **/
		 // THIS FUNCTION IS DEPRECATED
        //	DenseUnivariateRationalPolynomial integrate();
		
		/**
		 * Convert current object to its integral with constant of integration 0
		 *
		 **/ 
    	void integrate();
	
		/**
		 * Compute integral with constant of integration 0
		 *
		 **/ 
    	inline DenseUnivariateRationalPolynomial integral() const {
    	 	DenseUnivariateRationalPolynomial a(*this);
    	 	a.integrate();
    	 	return a;
    	 }
	
		/**
		 * Evaluate f(x)
		 *
		 * @param x: Rational evaluation point
		 **/
		 // THIS FUNCTION IS DEPRECATED
		RationalNumber evaluate(const RationalNumber& x) const;

		/**
		 * Evaluate f(x)
		 *
		 * @param x: Evaluation point
		 **/
    	Integer evaluate(const Integer& x) const;

		/**
		 * Evaluate f(x)
		 *
		 * @param x: Evaluation point in a larger ring, i.e. a ring in which the rationals can be embedded
		 **/
		template <class LargerRing>
		LargerRing evaluate(const LargerRing& x) const {
			// we might need a way of checking that this is always possible
			LargerRing a;
			if (curd) {
				LargerRing px = (LargerRing)coef[curd];
				for (int i = curd-1; i > -1; --i){
					a = (LargerRing)coef[i];
					px = (px * x) + a;
				}
				return px;
			}
			return (LargerRing)coef[0];
		}

		/**
		 * Is the least signficant coefficient zero
		 *
		 * @param
		 **/
		bool isConstantTermZero() const;

		/**
		 * GCD(p, q)
		 *
		 * @param q: The other polynomial
		 **/
    	DenseUnivariateRationalPolynomial gcd (const DenseUnivariateRationalPolynomial& q, int type) const;
		
    	inline DenseUnivariateRationalPolynomial gcd(const DenseUnivariateRationalPolynomial& q) const {
    		return gcd(q, 0);
    	}

		/**
		 * Square free
		 *
		 * @param
		 **/
    	Factors<DenseUnivariateRationalPolynomial> squareFree() const;

		/**
		 * Divide by variable if it is zero
		 *
		 * @param
		 **/
		bool divideByVariableIfCan();

		/**
		 * Number of coefficient sign variation
		 *
		 * @param
		 **/
		int numberOfSignChanges();

		/**
		 * Revert coefficients
		 *
		 * @param
		 **/

		void reciprocal();
		
		/**
		 * Homothetic operation
		 *
		 * @param k > 0: 2^(k*d) * f(2^(-k)*x);
		 **/
		void homothetic(int k=1);
		
		/**
		 * Scale transform operation
		 *
		 * @param k > 0: f(2^k*x)
		 **/
		void scaleTransform(int k);
		
		/**
		 * Compute f(-x)
		 *
		 * @param
		 **/
		void negativeVariable();

		/**
		 * Compute -f(x)
	         *
	         * @param
		 **/
		void negate();

		/**
		 * Return an integer k such that any positive root 
		   alpha of the polynomial satisfies alpha < 2^k
		 *
		 * @param
		 **/
		mpz_class rootBound();
		
		/**
		 * Taylor Shift operation by 1
		 *
		 * @param ts: Algorithm id
		 **/
		void taylorShift(int ts=-1);
		
		/**
		 * Positive real root isolation
		 * for square-free polynomials
		 *
		 * @param width: Interval's right - left < width
		 * @ts: Taylor Shift option: 0 - CMY; -1 - optimized 
		 **/ 
		inline Intervals positiveRealRootIsolate (mpq_class width, int ts=-1) {
			Intervals pIs;
			univariatePositiveRealRootIsolation(&pIs, this, width, ts);
			std::vector<Symbol> xs;
			xs.push_back(variable());
			pIs.setVariableNames(xs);
			return pIs;
		}

		/***********************************************************************
		* This function is an overload of the above function,
		* it is added to make BPAS compatable with gcc 6.0.
		* Reason: If cilk_spawn return value is a non primitive type, it results in an
		* invalid use of _Cilk_spawn error.
		*
		* Positive real root isolation
		* for square-free polynomials
		*
		* @param width: Interval's right - left < width
		* @ts: Taylor Shift option: 0 - CMY; -1 - optimized
		* @pIs: return value as a prameter
		**/

		inline void positiveRealRootIsolate (mpq_class width, Intervals pIs, int ts=-1) {
			//Intervals pIs;
			univariatePositiveRealRootIsolation(&pIs, this, width, ts);
			std::vector<Symbol> xs;
			xs.push_back(variable());
			pIs.setVariableNames(xs);
		}

		/**
		 * Real root isolation
		 * for square-free polynomials
		 *
		 * @param width: Interval's right - left < width
		 * @ts: Taylor Shift option: 0 - CMY; -1 - optimized
		 **/
		inline Intervals realRootIsolate (mpq_class width, int ts=-1) {
			Intervals pIs;
#if __GNUC__ == 4
			univariateRealRootIsolation(&pIs, this, width, ts);
#else

#endif
			return pIs;
		}
		
		/**
		 * Refine a root
		 *
		 * @param a: The root
		 * @param width: Interval's right - left < width
		 **/
		inline void refineRoot(Interval* a, mpq_class width) {
			refineUnivariateInterval(a, a->right+1, this, width);
		}
		
		/**
		 * Refine the roots
		 *
		 * @paran a: The roots
		 * @param width: Interval's right - left < width
		 **/
		inline Intervals refineRoots(Intervals& a, mpq_class width) {
			Intervals b;
			refineUnivariateIntervals(&b, &a, this, width);
			return b;
		}
		
		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param b: A univariate rational polynoial
		 **/ 
		void print(std::ostream &out) const;

		ExpressionTree convertToExpressionTree() const;



		/** BPASEuclideanDomain methods **/

		inline DenseUnivariateRationalPolynomial euclideanSize() const {
			return degree();
		}

		inline DenseUnivariateRationalPolynomial euclideanDivision(const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* q = NULL) const {
			std::cerr << "DenseUnivariateRationalPolynomial::euclideanDivision NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		inline DenseUnivariateRationalPolynomial extendedEuclidean(const DenseUnivariateRationalPolynomial& b, 
																		 DenseUnivariateRationalPolynomial* s = NULL,
																		 DenseUnivariateRationalPolynomial* t = NULL) const {
			std::cerr << "DenseUnivariateRationalPolynomial::extendedEuclidean NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		inline DenseUnivariateRationalPolynomial quotient(const DenseUnivariateRationalPolynomial& b) const {
			std::cerr << "DenseUnivariateRationalPolynomial::quotient NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}


		inline DenseUnivariateRationalPolynomial remainder(const DenseUnivariateRationalPolynomial& b) const {
			std::cerr << "DenseUnivariateRationalPolynomial::remainder NOT YET IMPLEMENTED" << std::endl;
			//TODO 
			return *this;
		}
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


