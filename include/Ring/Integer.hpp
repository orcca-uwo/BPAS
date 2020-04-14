
#ifndef _BPAS_INTEGER_H_
#define _BPAS_INTEGER_H_

#include "BPASEuclideanDomain.hpp"
#include <gmpxx.h>

//forward declarations
class RationalNumber;
class ComplexRationalNumber;
class SmallPrimeField;
class BigPrimeField;
class GeneralizedFermatPrimeField;
class DenseUnivariateIntegerPolynomial;
class DenseUnivariateRationalPolynomial;
template <class Ring>
class SparseUnivariatePolynomial;

/**
 * An arbitrary-precision Integer.
 */
class Integer : public BPASEuclideanDomain<Integer> {
	private:
		mpz_class _m;
	public:
		// static bool isPrimeField;
		// static bool isSmallPrimeField;
        // static bool isComplexField;

		/**
		 * Get a zero integer.
		 */
		Integer ();

		/**
		 * Construct and Integer from an int.
		 */
		Integer (int a);

		/**
		 * Construct an integer form a GMP mpz_t.
		 */
		Integer (const mpz_t& a);

		/**
		 * Construct an integer from a GMP mpz_class.
		 */
		Integer (const mpz_class& a);
		// Integer (mpz_class a);

		/**
		 * Copy assignment from another Integer.
		 */
		Integer (const Integer& a);

		/**
		 * Attempts to construct an Integer from a rational number.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const RationalNumber& a);

		/**
		 * Attempts to construct an Integer from a complex rational number.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const ComplexRationalNumber& a);

		/**
		 * Attempts to construct an Integer from a SmallPrimeField.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const SmallPrimeField& a);

		/**
		 * Attempts to construct an Integer from a BigPrimeField.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const BigPrimeField& a);

		/**
		 * Attempts to construct an Integer from a GeneralizedFermatPrimeField.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const GeneralizedFermatPrimeField& a);

		/**
		 * Attempts to construct an Integer from a DenseUnivariateIntegerPolynomial.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const DenseUnivariateIntegerPolynomial& a);

		/**
		 * Attempts to construct an Integer from a DenseUnivariateRationalPolynomial.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const DenseUnivariateRationalPolynomial& a);

		/**
		 * Attempts to construct an Integer from a SparseUnivariatePolynomial<Integer>.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const SparseUnivariatePolynomial<Integer>& a);

		/**
		 * Attempts to construct an Integer from a SparseUnivariatePolynomial<RationalNumber>.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const SparseUnivariatePolynomial<RationalNumber>& a);

		/**
		 * Attempts to construct an Integer from a SparseUnivariatePolynomial<ComplexRationalNumber>.
		 * If the conversion is not exact, causes an error.
		 */
		explicit Integer (const SparseUnivariatePolynomial<ComplexRationalNumber>& a);

		/**
		 * Attempts to construct an Integer from a SparseUniveraitePolynomial<Ring>
		 * for a generic ring. If the conversion is not exact, causes an error.
		 */
		template <class Ring>
		explicit Integer (const SparseUnivariatePolynomial<Ring>& a);

		inline mpz_class& get_mpz_ref() {
			return _m;
		}

		inline const mpz_class& get_mpz_ref() const {
			return _m;
		}

		inline mpz_class get_mpz () const {
			return _m;
		}

		inline mpz_ptr get_mpz_t() {
			return _m.get_mpz_t();
		}

		inline mpz_srcptr get_mpz_t() const {
			return _m.get_mpz_t();
		}


		inline double get_d() const {
			return _m.get_d();
		}

		inline long int get_si() const {
			return _m.get_si();
		}

		inline unsigned long int get_ui() const {
			return _m.get_ui();
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
	    Integer unitCanonical(Integer* u = NULL, Integer* v = NULL) const;

		/**
		 * Copy assignment.
		 */
		Integer& operator= (const Integer& a);

	    /**
	     * Addition.
	     */
		inline Integer operator+ (const Integer& i) const {
			Integer ret = *this;
			ret += i;
			return ret;
		}

	    /**
	     * Addition assignment.
	     */
	    inline Integer& operator+= (const Integer& i) {
	    	_m += i._m;
	    	return *this;
	    }

	    /**
	     * Subtraction.
	     */
	    inline Integer operator- (const Integer& i) const {
	    	Integer ret = *this;
	    	ret -= i;
	    	return ret;
	    }

	    /**
	     * Subtraction assignment.
	     */
	    inline Integer& operator-= (const Integer& i) {
	    	_m -= i._m;
	    	return *this;
	    }

	    /**
	     * Negation.
	     */
	    inline Integer operator- () const {
	    	Integer ret;
	    	mpz_neg(ret._m.get_mpz_t(), _m.get_mpz_t());
	    	return ret;
	    }

	    /**
	     * Multiplication.
	     */
	    inline Integer operator* (const Integer& i) const {
	    	Integer ret = *this;
	    	ret *= i;
	    	return ret;
	    }

	    /**
	     * Multiplication assignment.
	     */
	    inline Integer& operator*= (const Integer& i) {
	    	_m *= i._m;
	    	return *this;
	    }

	    /**
	     * Exponentiation.
	     */
	    inline Integer operator^ (long long int e) const {
			Integer r;
			mpz_pow_ui(r._m.get_mpz_t(), _m.get_mpz_t(), (unsigned long int) e);
			return r;
	    }

	    /**
	     * Exponentiation assignment.
	     */
	    inline Integer& operator^= (long long int e) {
	    	*this = *this ^ e;
	    	return *this;
	    }

	    /**
	     * Equality test,
	     *
	     * returns true iff equal
	     */
	    inline bool operator== (const Integer& i) const {
	    	return (mpz_cmp(_m.get_mpz_t(), i._m.get_mpz_t()) == 0);
	    }

	    /**
	     * Inequality test,
	     *
	     * returns true iff not equal.
	     */
	    inline bool operator!= (const Integer& i) const {
	    	return (mpz_cmp(_m.get_mpz_t(), i._m.get_mpz_t()) != 0);
	    }

	    inline bool operator< (const Integer& r) const {
			return _m < r._m;
		}

		inline bool operator<= (const Integer& r) const {
			return _m <= r._m;
		}

		inline bool operator> (const Integer& r) const {
			return _m > r._m;
		}

		inline bool operator>= (const Integer& r) const {
			return _m >= r._m;
		}

	    /**
	     * Convert this to an expression tree.
	     *
	     * returns an expression tree describing *this.
	     */
		inline ExpressionTree convertToExpressionTree() const {
			return ExpressionTree(new ExprTreeNode(_m));

		}

		/**
		 * Exact division.
		 */
		inline Integer operator/ (const Integer& i) const {
			//TODO ensure this is exact and not rounded
			Integer ret = *this;
			ret /= i;
			return ret;
		}

		/**
		 * Exact division assignment.
		 */
		inline Integer& operator/= (const Integer& i) {
			if (mpz_divisible_p(_m.get_mpz_t(), i._m.get_mpz_t())) {
				mpz_divexact(_m.get_mpz_t(), _m.get_mpz_t(), i._m.get_mpz_t());
			} else {
				std::cerr << "BPAS ERROR: Non-exact division in Integer: " << _m << " / " << i._m << std::endl;
				exit(1);
			}
			return *this;
		}

		inline bool divisible(unsigned long int i) {
			return mpz_divisible_ui_p(_m.get_mpz_t(), i);
		}

		inline Integer operator% (const Integer& r) const {
			Integer ret(this->_m % r._m);
			return ret;
		}

		inline Integer& operator%= (const Integer& r) {
			_m = _m % r._m;
			return *this;
		}

		/**
		 * Get GCD of *this and other.
 		 */
		inline Integer gcd(const Integer& other) const {
			Integer c;
			mpz_gcd(c._m.get_mpz_t(), _m.get_mpz_t(), other._m.get_mpz_t());
			return c;
		}
		/**
		 * Compute squarefree factorization of *this
		 */
		inline Factors<Integer> squareFree() const {
			std::vector<Integer> ret;
			ret.push_back(*this);
			return ret;
		}

		/**
		 * Get the euclidean size of *this.
		 */
		inline Integer euclideanSize() const {
			if (*this < 0) {
				return -(*this);
			}
			return *this;
		}

		/**
		 * Perform the eucldiean division of *this and b. Returns the
		 * remainder. If q is not NULL, then returns the quotient in q.
		 */
		Integer euclideanDivision(const Integer& b, Integer* q = NULL) const;

		/**
		 * Perform the extended euclidean division on *this and b.
		 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
		 */
		Integer extendedEuclidean(const Integer& b, Integer* s = NULL, Integer* t = NULL) const;

		/**
		 * Get the quotient of *this and b.
		 */
		Integer quotient(const Integer& b) const;

		/**
		 * Get the remainder of *this and b.
		 */
		Integer remainder(const Integer& b) const;


		inline friend Integer operator+(int a, const Integer& r) {
			return Integer(mpz_class(a) + r._m);
		}

		inline friend Integer operator-(int a, const Integer& r) {
			return Integer(mpz_class(a) - r._m);
		}

		inline friend Integer operator*(int a, const Integer& r) {
			return Integer(mpz_class(a) * r._m);
		}

		inline friend Integer operator/(int a, const Integer& r) {
			return Integer(mpz_class(a) / r._m);
		}

		inline friend bool operator<(int a, const Integer& r) {
			return r > a;
		}

		inline friend bool operator<=(int a, const Integer& r) {
			return r >= a;
		}

		inline friend bool operator>(int a, const Integer& r) {
			return r < a;
		}

		inline friend bool operator>=(int a, const Integer& r) {
			return r <= a;
		}

		inline friend Integer operator+(long int a, const Integer& r) {
			return Integer(mpz_class(a) + r._m);
		}

		inline friend Integer operator-(long int a, const Integer& r) {
			return Integer(mpz_class(a) - r._m);
		}

		inline friend Integer operator*(long int a, const Integer& r) {
			return Integer(mpz_class(a) * r._m);
		}

		inline friend Integer operator/(long int a, const Integer& r) {
			return Integer(mpz_class(a) / r._m);
		}

		inline friend Integer abs(const Integer& i) {
			return Integer(abs(i._m));
		}

		inline friend bool operator<(long int a, const Integer& r) {
			return r > a;
		}

		inline friend bool operator<=(long int a, const Integer& r) {
			return r >= a;
		}

		inline friend bool operator>(long int a, const Integer& r) {
			return r < a;
		}

		inline friend bool operator>=(long int a, const Integer& r) {
			return r <= a;
		}


};

#endif
