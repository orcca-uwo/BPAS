
#ifndef _BPAS_RATIONAL_NUMBER_H_
#define _BPAS_RATIONAL_NUMBER_H_

#include "BPASField.hpp"
#include <gmpxx.h>
#include <iostream>

//forward declarations
class Integer;
class ComplexRationalNumber;
class SmallPrimeField;
class BigPrimeField;
class GeneralizedFermatPrimeField;
class DenseUnivariateIntegerPolynomial;
class DenseUnivariateRationalPolynomial;
template <class Ring>
class SparseUnivariatePolynomial;
class SparseMultivariateRationalPolynomial;

/**
 * An arbitrary-precision rational number.
 */
class RationalNumber : public BPASField<RationalNumber> {
	private:
	 	mpq_class _m;

	public:

		// static bool isPrimeField;
		// static bool isSmallPrimeField;
        // static bool isComplexField;

		RationalNumber ();

		RationalNumber (int a, int b = 1);

		RationalNumber (const std::string& digits, int base = 10);

		RationalNumber (const mpq_t& q);

		RationalNumber (const mpq_class& a);

		RationalNumber (const mpz_class& a, const mpz_class& b = mpz_class(1));

		RationalNumber (const RationalNumber& a);

		explicit RationalNumber (const Integer& a);

		explicit RationalNumber (const ComplexRationalNumber& a);

		explicit RationalNumber (const SmallPrimeField& a);

		explicit RationalNumber (const BigPrimeField& a);

		explicit RationalNumber (const GeneralizedFermatPrimeField& a);

		explicit RationalNumber (const DenseUnivariateIntegerPolynomial& a);

		explicit RationalNumber (const DenseUnivariateRationalPolynomial& a);

		explicit RationalNumber (const SparseUnivariatePolynomial<Integer>& a);

		explicit RationalNumber (const SparseUnivariatePolynomial<RationalNumber>& a);

		explicit RationalNumber (const SparseUnivariatePolynomial<ComplexRationalNumber>& a);

		explicit RationalNumber (const SparseMultivariateRationalPolynomial& a);

		template <class Ring>
		explicit RationalNumber (const SparseUnivariatePolynomial<Ring>& a);

		RationalNumber* RNpointer(RationalNumber* a);

		RationalNumber* RNpointer(SmallPrimeField* a);

		RationalNumber* RNpointer(BigPrimeField* a);

		RationalNumber* RNpointer(GeneralizedFermatPrimeField* a);

		RationalNumber& set (int a, int b);

		inline mpq_class get_mpq() const {
			return _m;
		}

		inline mpq_class& get_mpq_ref() {
			return _m;
		}

		inline const mpq_class& get_mpq_ref() const {
			return _m;
		}

		inline mpq_ptr get_mpq_t() {
			return _m.get_mpq_t();
		}

		inline mpq_srcptr get_mpq_t() const {
			return _m.get_mpq_t();
		}

		inline Integer get_num() const {
			return Integer(_m.get_num());
		}

		double get_d() const {
			return _m.get_d();
		}

		inline Integer get_den() const {
			return Integer(_m.get_den());
		}

		/**
		 * Is a zero
		 *
		 * @param
		 **/
		inline bool isZero() const {
			return (_m == 0);
		}

		/**
		 * Assign to zero
		 *
		 * @param
		 **/
		inline void zero() {
			_m = 0;
		}

		/**
		 * Is a 1
		 *
		 * @param
		 **/
		inline bool isOne() const {
			return (_m == 1);
		}

		/**
		 * Assign to one
		 *
		 * @param
		 **/
		inline void one() {
			_m = 1;
		}

		/**
		 * Is a -1
		 *
		 * @param
		 **/
		inline bool isNegativeOne() const {
			return (_m == -1);
		}

		/**
		 * Assign to negative one
		 *
		 * @param
		 **/
		inline void negativeOne() {
			_m = -1;
		}

		/**
 		 * Is a constant
 		 *
 		 * @param
 		 **/
		inline int isConstant() const {
			if (_m >= 0)
				return 1;
			else { return -1; }
		}

		/**
	     * Obtain the unit normal (a.k.a canonical associate) of an element.
	     * If either parameters u, v, are non-NULL then the units are returned such that
	     * b = ua, v = u^-1. Where b is the unit normal of a, and is the returned value.
	     */
	    RationalNumber unitCanonical(RationalNumber* u = NULL, RationalNumber* v = NULL) const;

		/**
		 * Copy assignment.
		 */
		RationalNumber& operator= (const RationalNumber& a);

	    /**
	     * Addition.
	     */
		inline RationalNumber operator+ (const RationalNumber& i) const {
			RationalNumber ret = *this;
			ret += i;
			return ret;
		}

	    /**
	     * Addition assignment.
	     */
	    inline RationalNumber& operator+= (const RationalNumber& i) {
	    	_m += i._m;
	    	return *this;
	    }

	    /**
	     * Subtraction.
	     */
	    inline RationalNumber operator- (const RationalNumber& i) const {
	    	RationalNumber ret = *this;
	    	ret -= i;
	    	return ret;
	    }

	    /**
	     * Subtraction assignment.
	     */
	    inline RationalNumber& operator-= (const RationalNumber& i) {
	    	_m -= i._m;
	    	return *this;
	    }

	    /**
	     * Negation.
	     */
	    inline RationalNumber operator- () const {
	    	RationalNumber ret;
	    	mpq_neg(ret._m.get_mpq_t(), (*this)._m.get_mpq_t());
	    	return ret;
	    }

	    /**
	     * Multiplication.
	     */
	    inline RationalNumber operator* (const RationalNumber& i) const {
	    	RationalNumber ret = *this;
	    	ret *= i;
	    	return ret;
	    }

	    /**
	     * Multiplication assignment.
	     */
	    inline RationalNumber& operator*= (const RationalNumber& i) {
	    	_m *= i._m;
	    	return *this;
	    }

	    /**
	     * Equality test,
	     *
	     * returns true iff equal
	     */
	    inline bool operator== (const RationalNumber& i) const {
	    	return (_m == i._m);
	    }

	    /**
	     * Inequality test,
	     *
	     * returns true iff not equal.
	     */
	    inline bool operator!= (const RationalNumber& i) const {
	    	return (_m != i._m);
	    }

	    inline ExpressionTree convertToExpressionTree() const {
	    	return ExpressionTree(new ExprTreeNode(_m));
	    }

	    /**
		 * Exact division.
		 */
		inline RationalNumber operator/ (const RationalNumber& i) const {
			//TODO ensure this is exact and not rounded
			RationalNumber ret = *this;
			ret /= i;
			return ret;
		}

		/**
		 * Exact division assignment.
		 */
		inline RationalNumber& operator/= (const RationalNumber& i) {
			_m /= i._m;
			return *this;
		}

		inline bool operator< (const RationalNumber& r) const {
			return _m < r._m;
		}

		inline bool operator<= (const RationalNumber& r) const {
			return _m <= r._m;
		}

		inline bool operator> (const RationalNumber& r) const {
			return _m > r._m;
		}

		inline bool operator>= (const RationalNumber& r) const {
			return _m >= r._m;
		}

		/**
		 * GCD(a, b)
		 *
		 * @param b: The other rational number
		 **/
		inline RationalNumber gcd (const RationalNumber& b) const {
			RationalNumber c;
			if (this->isZero() && b.isZero())
				c = 0;
			else
				c = 1;
			return c;
		}

		/**
		 * Compute squarefree factorization of *this
		 */
		inline Factors<RationalNumber> squareFree() const {
			std::vector<RationalNumber> ret;
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
		RationalNumber euclideanDivision(const RationalNumber& b, RationalNumber* q = NULL) const;

		/**
		 * Perform the extended euclidean division on *this and b.
		 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
		 */
		RationalNumber extendedEuclidean(const RationalNumber& b, RationalNumber* s = NULL, RationalNumber* t = NULL) const;

		/**
		 * Get the quotient of *this and b.
		 */
		RationalNumber quotient(const RationalNumber& b) const;

		/**
		 * Get the remainder of *this and b.
		 */
		RationalNumber remainder(const RationalNumber& b) const;


		/**
		 * Overload operator ^
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation
		 **/
		inline RationalNumber operator^ (long long int e) const {
			RationalNumber r;
			mpz_pow_ui(r._m.get_num_mpz_t(), _m.get_num_mpz_t(), (unsigned long int) e);
			mpz_pow_ui(r._m.get_den_mpz_t(), _m.get_den_mpz_t(), (unsigned long int) e);
			return r;
		}

		inline RationalNumber& operator^= (long long int e) {
			*this = *this ^ e;
			return *this;
		}

		inline RationalNumber operator% (const RationalNumber& r) const {
			return 0;
		}

		inline RationalNumber& operator%= (const RationalNumber& r) {
			*this = 0;
			return *this;
		}

		inline RationalNumber inverse() const {
			RationalNumber ret;
			mpz_set(ret._m.get_den_mpz_t(), _m.get_num_mpz_t());
			mpz_set(ret._m.get_num_mpz_t(), _m.get_den_mpz_t());
			ret._m.canonicalize();
			return ret;
		}



		//Friend functions

		inline friend RationalNumber operator+(int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) + r._m);
		}

		inline friend RationalNumber operator-(int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) - r._m);
		}

		inline friend RationalNumber operator*(int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) * r._m);
		}

		inline friend RationalNumber operator/(int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) / r._m);
		}

		inline friend RationalNumber operator+(long int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) + r._m);
		}

		inline friend RationalNumber operator-(long int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) - r._m);
		}

		inline friend RationalNumber operator*(long int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) * r._m);
		}

		inline friend RationalNumber operator/(long int a, const RationalNumber& r) {
			return RationalNumber(mpq_class(a) / r._m);
		}

		inline friend RationalNumber abs(const RationalNumber& i) {
			return RationalNumber(abs(i._m));
		}





};

#endif
