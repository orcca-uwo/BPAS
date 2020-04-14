
#ifndef _COMPLEX_RATIONAL_NUMBER_H_
#define _COMPLEX_RATIONAL_NUMBER_H_

#include "BPASField.hpp"
#include <iostream>

//forward declarations
class Integer;
class RationalNumber;
class SmallPrimeField;
class BigPrimeField;
class GeneralizedFermatPrimeField;
class DenseUnivariateIntegerPolynomial;
class DenseUnivariateRationalPolynomial;
template <class Ring>
class SparseUnivariatePolynomial;


/**
 * An arbitrary-precision complex rational number.
 */
class ComplexRationalNumber : public BPASField<ComplexRationalNumber> {

private:
	// a + i * b
	mpq_class a;
	mpq_class b;

public:



	// static bool isPrimeField;
	// static bool isSmallPrimeField;
    // static bool isComplexField;

	ComplexRationalNumber ();

	ComplexRationalNumber (const mpq_class& _a, const mpq_class& _b = mpq_class(1));

	ComplexRationalNumber (const ComplexRationalNumber& c);

	ComplexRationalNumber(int _a, int _b = 1, int _c = 0, int _d = 1);

	explicit ComplexRationalNumber (const Integer& c);

	explicit ComplexRationalNumber (const RationalNumber& c);

	explicit ComplexRationalNumber (const SmallPrimeField& c);

	explicit ComplexRationalNumber (const BigPrimeField& c);

	explicit ComplexRationalNumber (const GeneralizedFermatPrimeField& c);

	explicit ComplexRationalNumber (const DenseUnivariateIntegerPolynomial& c);

	explicit ComplexRationalNumber (const DenseUnivariateRationalPolynomial& c);

	explicit ComplexRationalNumber (const SparseUnivariatePolynomial<Integer>& c);

	explicit ComplexRationalNumber (const SparseUnivariatePolynomial<RationalNumber>& c);

	explicit ComplexRationalNumber (const SparseUnivariatePolynomial<ComplexRationalNumber>& c);

	template <class Ring>
	explicit ComplexRationalNumber (const SparseUnivariatePolynomial<Ring>& c);

	ComplexRationalNumber& operator= (const ComplexRationalNumber& c);

	ComplexRationalNumber& operator= (const mpq_class& k);

	ComplexRationalNumber& operator= (int k);

	ComplexRationalNumber& setRealPart (const RationalNumber& r);

	ComplexRationalNumber& setRealPart (const mpq_class& k);

	ComplexRationalNumber& setRealPart (int k);

	ComplexRationalNumber& setImaginaryPart (const RationalNumber& r);

	ComplexRationalNumber& setImaginaryPart (const mpq_class& k);

	ComplexRationalNumber& setImaginaryPart (int k);

	ComplexRationalNumber& set (const RationalNumber& ka, const RationalNumber& kb);

	ComplexRationalNumber& set (const mpq_class& ka, const mpq_class& kb);

	ComplexRationalNumber& set (const mpq_class& ka, int kb);

	ComplexRationalNumber& set (int ka, const mpq_class& kb);

	ComplexRationalNumber& set (int ka, int kb);

	/**
	 * Is a zero
	 *
	 * @param
	 **/
	inline bool isZero() const {
		return (a == 0 && b == 0);
	}

	/**
	 * Assign to zero
	 *
	 * @param
	 **/
	inline void zero() {
		a = 0;
		b = 0;
	}

	/**
	 * Is a 1
	 *
	 * @param
	 **/
	inline bool isOne() const {
		return (a == 1 && b == 0);
	}

	/**
	 * Assign to one
	 *
	 * @param
	 **/
	inline void one() {
		a = 1;
		b = 0;
	}

	/**
	 * Is a -1
	 *
	 * @param
	 **/
	inline bool isNegativeOne() const {
		return (a == -1 && b == 0);
	}

	/**
	 * Assign to negative one
	 *
	 * @param
	 **/
	inline void negativeOne() {
		a = -1;
		b = 0;
	}

	/**
	 * Is a constant
	 *
	 * @param
	 **/
	inline int isConstant() const {
	  if (a >= 0)
		return 1;
	  else { return -1; }
	}

	ComplexRationalNumber unitCanonical(ComplexRationalNumber* u = NULL, ComplexRationalNumber* v = NULL) const;

	inline bool operator== (const ComplexRationalNumber& c) const {
		if (a == c.a && b == c.b)
			return 1;
		else { return 0; }
	}

	inline bool operator== (const mpq_class& k) const {
		if (a == k && b == 0)
			return 1;
		else { return 0; }
	}

	inline bool operator== (int k) const {
		if (a == k && b == 0)
			return 1;
		else { return 0; }
	}

	inline bool operator!= (const ComplexRationalNumber& c) const {
		if (a == c.a && b == c.b)
			return 0;
		else { return 1; }
	}

	inline bool operator!= (const mpq_class& k) const {
		if (a == k && b == 0)
			return 1;
		else { return 0; }
	}

	inline bool operator!= (int k) const {
		if (a == k && b == 0)
			return 1;
		else { return 0; }
	}

	inline ComplexRationalNumber operator+ (const ComplexRationalNumber& c) const {
		ComplexRationalNumber r (*this);
		return (r += c);
	}

	inline ComplexRationalNumber& operator+= (const ComplexRationalNumber& c) {
		a += c.a;
		b += c.b;
		return *this;
	}

	inline ComplexRationalNumber operator- (const ComplexRationalNumber& c) const {
	 	ComplexRationalNumber r (*this);
		return (r -= c);
	}

	inline ComplexRationalNumber& operator-= (const ComplexRationalNumber& c) {
		a -= c.a;
		b -= c.b;
		return *this;
	}

	inline ComplexRationalNumber operator- () const {
		ComplexRationalNumber r (-a, -b);
		return r;
	}

	inline ComplexRationalNumber operator* (const ComplexRationalNumber& c) const {
		ComplexRationalNumber r (*this);
		return (r *= c);
	}

	inline ComplexRationalNumber& operator*= (const ComplexRationalNumber& c) {
		mpq_class t = a*c.a - b*c.b;
		mpq_class e = a*c.b + c.a*b;
		a = t;
		b = e;
		return *this;
	}

	inline ComplexRationalNumber& operator*= (const mpq_class& c) {
		a *= c;
		b *= c;
		return *this;
	}

	inline ComplexRationalNumber& operator*= (int c) {
		a *= c;
		b *= c;
		return *this;
	}

	/**
	 * Overload operator ^
	 * replace xor operation by exponentiation
	 *
	 * @param e: The exponentiation
	 **/
	inline ComplexRationalNumber operator^ (long long int e) const {
		ComplexRationalNumber r;
		if (isZero() || isOne() || e == 1)
			r = *this;
		else if (e == 2) {
			r = *this * *this;
		}
		else if (e > 2) {
			ComplexRationalNumber x (*this);
			r.one();

			while (e != 0) {
				if (e % 2)
					r *= x;
				x = x * x;
				e >>= 1;
			}
		}
		else if (e == 0) {
			r.one();
		}
		else {
			r = *this ^ (-e);
			r.inverse();
		}
		return r;
	}

	inline ComplexRationalNumber& operator^= (long long int e) {
		*this = *this ^ e;
		return *this;
	}

	inline ExpressionTree convertToExpressionTree() const {
		std::cerr << "ComplexRationalNumber::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
		exit(1);
		return ExpressionTree();
	}

	inline ComplexRationalNumber operator/ (const ComplexRationalNumber& c) const {
		ComplexRationalNumber r (*this);
		return (r /= c);
	}

	inline ComplexRationalNumber& operator/= (const ComplexRationalNumber& c) {
		if (c.isZero()) {
			std::cout << "BPAS: error, dividend is zero from ComplexRationalNumber."<< std::endl;
        	exit(1);
		}
		mpq_class r = c.a*c.a + c.b*c.b;
		mpq_class t = (a*c.a+b*c.b)/r;
		mpq_class e = (b*c.a-c.b*a)/r;
		a = t;
		b = e;
		return *this;
	}

	inline ComplexRationalNumber operator% (const ComplexRationalNumber& c) const {
		return 0;
	}

	inline ComplexRationalNumber& operator%= (const ComplexRationalNumber& c) {
		*this = 0;
		return *this;
	}

	/**
	 * GCD(a, b)
	 *
	 * @param b: The other rational number
	 **/
	inline ComplexRationalNumber gcd (const ComplexRationalNumber& c) const {
		ComplexRationalNumber e;
		if (isZero() && c.isZero())
			e.zero();
		else
			e.one();
		return e;
	}

	/**
	 * Compute squarefree factorization of *this
	 */
	inline Factors<ComplexRationalNumber> squareFree() const {
		std::vector<ComplexRationalNumber> ret;
		ret.push_back(*this);
		return ret;
	}

	/**
	 * Get the euclidean size of *this.
	 */
	inline Integer euclideanSize() const {
		return Integer(1);
	}

	/**
	 * Perform the eucldiean division of *this and b. Returns the
	 * remainder. If q is not NULL, then returns the quotient in q.
	 */
	ComplexRationalNumber euclideanDivision(const ComplexRationalNumber& b, ComplexRationalNumber* q = NULL) const;

	/**
	 * Perform the extended euclidean division on *this and b.
	 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
	 */
	ComplexRationalNumber extendedEuclidean(const ComplexRationalNumber& b, ComplexRationalNumber* s = NULL, ComplexRationalNumber* t = NULL) const;

	/**
	 * Get the quotient of *this and b.
	 */
	ComplexRationalNumber quotient(const ComplexRationalNumber& b) const;

	/**
	 * Get the remainder of *this and b.
	 */
	ComplexRationalNumber remainder(const ComplexRationalNumber& b) const;

	inline ComplexRationalNumber inverse() const {
		ComplexRationalNumber r;
		mpq_class e = a * a + b * b;
		r.a = a/e;
		r.b = -b/e;
		return r;
	}

	inline RationalNumber realPart() const {
		return RationalNumber(a);
	}

	inline RationalNumber imaginaryPart() const {
		return RationalNumber(b);
	}

	inline ComplexRationalNumber conjugate() const {
		ComplexRationalNumber r(a, -b);
		return r;
	}

	void print(std::ostream& out) const;

};


#endif
