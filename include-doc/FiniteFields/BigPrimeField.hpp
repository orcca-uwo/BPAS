
#ifndef _BIG_PRIME_FIELD_H_
#define _BIG_PRIME_FIELD_H_

#include "BPASFiniteField.hpp"
#include <iostream>
#include <string>

using std::endl;
using std::cout;

//forward declarations
class Integer;
class RationalNumber;
class ComplexRationalNumber;
class SmallPrimeField;
class GeneralizedFermatPrimeField;
class DenseUnivariateIntegerPolynomial;
class DenseUnivariateRationalPolynomial;
template <class Ring>
class SparseUnivariatePolynomial;

// prime field with big number, using GMP functions
/**
 * A prime field whose prime can be arbitrarily large. 
 */
class BigPrimeField : public BPASFiniteField {

private:

	static mpz_class prime;
	mpz_class a;

public:

	static mpz_class characteristic;
    static RingProperties properties;        
	// static bool isPrimeField;
	// static bool isSmallPrimeField;
	// static bool isComplexField;

	BigPrimeField ();

	BigPrimeField (mpz_class _a);

	BigPrimeField (long int _a);

	BigPrimeField (const BigPrimeField& c);

	explicit BigPrimeField (const Integer& c);

	explicit BigPrimeField (const RationalNumber& c);

	explicit BigPrimeField (const ComplexRationalNumber& c);

	explicit BigPrimeField (const SmallPrimeField& c); //? NOT SURE

	explicit BigPrimeField (const GeneralizedFermatPrimeField& c);

	explicit BigPrimeField (const DenseUnivariateIntegerPolynomial& c);

	explicit BigPrimeField (const DenseUnivariateRationalPolynomial& c);

	explicit BigPrimeField (const SparseUnivariatePolynomial<Integer>& c);

	explicit BigPrimeField (const SparseUnivariatePolynomial<RationalNumber>& c);

	explicit BigPrimeField (const SparseUnivariatePolynomial<ComplexRationalNumber>& c);

	template <class Ring>
	explicit BigPrimeField (const SparseUnivariatePolynomial<Ring>& c);

	BigPrimeField* BPFpointer(BigPrimeField* b);

	BigPrimeField* BPFpointer(RationalNumber* a);

	BigPrimeField* BPFpointer(SmallPrimeField* a);

	BigPrimeField* BPFpointer(GeneralizedFermatPrimeField* a);

	static void setPrime(mpz_class p){
		prime = p;
	//	characteristic = p;
	}

	mpz_class Prime() const;

	mpz_class number() const;

	void whichprimefield();

	BigPrimeField findPrimitiveRootOfUnity(long int n) const {
		return BigPrimeField::findPrimitiveRootofUnity(mpz_class(n));
	}


	static BigPrimeField findPrimitiveRootofUnity(mpz_class n){
		if ( ((prime - 1) % n != 0)){
			cout << "ERROR: n does not divide prime - 1." << endl;
			return -1;
		}
		bool flag = false;
		mpz_class  q = (prime - 1) / n;
		BigPrimeField p1 (prime - 1);
		BigPrimeField c;
		int i = 0;
		mpz_class test = q * n / 2;
		srand (time(NULL));
		while(i < 20){
			c =  rand();
			if ((c^test) == p1) {
				flag = true;
				return (c^q);
			}
			i++;
		}
		if (!flag ){
			cout << "No primitive root found!"<< endl;
			return 0;
		}
	// If no primitive root found
	}

  	BigPrimeField unitCanonical(BigPrimeField* u = NULL, BigPrimeField* v = NULL) const;

	BigPrimeField& operator= (const BigPrimeField& c);

	BigPrimeField& operator= (long int k);

	BigPrimeField& operator= (const mpz_class& k);

	inline bool isZero() const {
		return (a == 0);
	}

	inline void zero() {
		a = 0;
	}

	inline bool isOne() const {
		return (a == 1);
	}

	inline void one() {
		a = 1;
	}

	inline bool isNegativeOne() const {
		return (a == (prime - 1));
	}

	inline void negativeOne() {
		a = prime - 1;
	}

	inline int isConstant() const {
		if (a >= 0)
			return 1;
		else { return -1; }
	}

	inline bool operator== (const BigPrimeField& c) const {
		if (a == c.a)
			return 1;
		else { return 0; }
	}

	inline bool operator== (const mpz_class& k) const {
		BigPrimeField r (*this);
		BigPrimeField b (k);
		if (b == r){
			return 1;
		}
		else {
			return 0;
		}
	}

	inline bool operator== (long int k) const {
		BigPrimeField r (*this);
		BigPrimeField b (k);
		if (b == r){
			return 1;
		}
		else {
			return 0;
		}
	}

	inline bool operator!= (const BigPrimeField& c) const {
		if (a == c.a)
			return 0;
		else { return 1; }
	}

	inline bool operator!= (const mpz_class& k) const {
		BigPrimeField r (*this);
		BigPrimeField b (k);
		if (b == r){
			return 0;
		}
		else {
			return 1;
		}
	}

	inline bool operator!= (long int k) const {
		BigPrimeField r (*this);
		BigPrimeField b (k);
		if (b == r){
			return 0;
		}
		else {
			return 1;
		}
	}

	inline BigPrimeField operator+ (const BigPrimeField& c) const {
		BigPrimeField r (*this);
		return (r += c);
	}

	inline BigPrimeField operator+ (long int c) const {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return (r += b);
	}

	inline BigPrimeField operator+ (const mpz_class& c) const {
		BigPrimeField r (*this);
		BigPrimeField b(c);
		return (r += b);
	}

	inline BigPrimeField& operator+= (const BigPrimeField& c) {
		a = (a + c.a);
		if(a>prime)
			a -= prime;
		return *this;
	}

	inline BigPrimeField operator+= (long int c) {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return (r += b);
	}

	inline BigPrimeField operator+= (const mpz_class& c) {
		BigPrimeField r (*this);
		BigPrimeField b(c);
		return (r += b);
	}

	inline BigPrimeField operator- (const BigPrimeField& c) const {
		BigPrimeField r (*this);
		return (r -= c);
	}
	
	inline BigPrimeField operator- (long int c) const {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return (r += b);
	}

	inline BigPrimeField operator- (const mpz_class& c) const {
		BigPrimeField r (*this);
		BigPrimeField b(c);
		return (r += b);
	}
	
	inline BigPrimeField& operator-= (const BigPrimeField& c) {
		if ((a - c.a)<0){
			a = prime+(a - c.a);
		}
		else{
			a = a - c.a;
		}
		return *this;
	}

	inline BigPrimeField operator-= (long int c) {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return (r += b);
	}
	
	inline BigPrimeField operator-= (const mpz_class& c) {
		BigPrimeField r (*this);
		BigPrimeField b(c);
		return (r += b);
	}

	inline BigPrimeField operator- () const {
		BigPrimeField ret = *this;
		ret.a = prime - a;
		return ret;
	}

	inline BigPrimeField operator* (const BigPrimeField& c) const {
		BigPrimeField r (*this);
		return (r *= c);
	}

	inline BigPrimeField operator* (const mpz_class& c) const {
		BigPrimeField r (*this);
		return (r *= c);
	}

	inline BigPrimeField operator* (long int c) const {
		BigPrimeField r (*this);
		return (r *= c);
	}

	inline BigPrimeField& operator*= (const BigPrimeField& c) {
		a = (a * c.a)%prime;
		return *this;
	}

	inline BigPrimeField& operator*= (const mpz_class& m) {
		mpz_class c(m);
		while(c<0){
			c=c+prime;
		}
		a = (a * c)%prime;
		return *this;
	}

	inline BigPrimeField& operator*= (long int c) {
		mpz_class b = c;
		while(b<0){
			b=b+prime;
		}
		a = (a * b)%prime;
		return *this;
	}

	inline BigPrimeField operator^ (long long int c) const {
		BigPrimeField r (*this);
		mpz_class b(std::to_string(c), 10);
		return (r ^ b);
	}

	inline BigPrimeField operator^ (const mpz_class& exp) const {
		BigPrimeField r;
		mpz_class e = exp;
		if (isZero() || isOne() || e == 1)
			r = *this;
		else if (e == 2) {
			r = *this * *this;
		}
		else if (e > 2) {
			BigPrimeField x (*this);
			r.one();

			while (e != 0) {
				if ((e % 2) == 1)
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
			r=r.inverse();
		}
		return r;
	}

	inline BigPrimeField& operator^= (long long int c) {
		*this = *this ^ c;
		return *this;
	}

	inline BigPrimeField& operator^= (const mpz_class& e) {
		*this = *this ^ e;
		return *this;
	}

	inline ExpressionTree convertToExpressionTree() const {
		return ExpressionTree(new ExprTreeNode(a));
	}

	inline BigPrimeField operator/ (const BigPrimeField& c) const {
		BigPrimeField r (*this);
		return (r /= c);
	}

	inline BigPrimeField operator/ (long int c) const {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return (r /= b);
	}

	inline BigPrimeField operator/ (const mpz_class& c) const {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return (r /= b);
	}

	inline BigPrimeField& operator/= (const BigPrimeField& c) {
		if (c.isZero()) {
			std::cout << "BPAS: error, dividend is zero from SmallPrimeField."<< std::endl;
			exit(1);
		}
		BigPrimeField inv = c.inverse();
		*this *= inv;
		return *this;
	}

	inline BigPrimeField& operator/= (long int c) {
		BigPrimeField b (c);
		return (*this /= b);
	}

	inline BigPrimeField& operator/= (const mpz_class& c) {
		BigPrimeField b (c);
		return (*this /= b);
	}
	
	inline BigPrimeField operator% (const BigPrimeField& c) const {
		return 0;
	}

	inline BigPrimeField& operator%= (const BigPrimeField& c) {
		*this = 0;
		return *this;
	}

	inline BigPrimeField gcd (const BigPrimeField& other) const {
		BigPrimeField q (0);
		BigPrimeField r (0);
		BigPrimeField c (other);
		BigPrimeField b (*this);
		if(b.a < c.a){
			return c.gcd(b);
		}
		while (c.a > 0) {
			q.a = b.a / c.a;
			r.a = b.a % c.a;
			b = c;
			c = r;
		}
		return b;
	}

	inline BigPrimeField gcd (long int c){
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return r.gcd(b);
	}

	inline BigPrimeField gcd (const mpz_class& c) const {
		BigPrimeField r (*this);
		BigPrimeField b (c);
		return r.gcd(b);
	}

	/**
	 * Compute squarefree factorization of *this
	 */
	inline Factors<BigPrimeField> squareFree() const {
		std::vector<BigPrimeField> ret;
		ret.push_back(*this);
		return ret;
	}

	/** 
	 * Get the euclidean size of *this.
	 */
	inline BigPrimeField euclideanSize() const {
		return (*this).number();
	}
	
	/**
	 * Perform the eucldiean division of *this and b. Returns the 
	 * remainder. If q is not NULL, then returns the quotient in q. 
	 */ 
	BigPrimeField euclideanDivision(const BigPrimeField& b, BigPrimeField* q = NULL) const;

	/**
	 * Perform the extended euclidean division on *this and b. 
	 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
	 */
	BigPrimeField extendedEuclidean(const BigPrimeField& b, BigPrimeField* s = NULL, BigPrimeField* t = NULL) const;
	
	/**
	 * Get the quotient of *this and b.
	 */
	BigPrimeField quotient(const BigPrimeField& b) const;

	/** 
	 * Get the remainder of *this and b.
	 */
	BigPrimeField remainder(const BigPrimeField& b) const;

	inline BigPrimeField inverse() const {
		BigPrimeField r(0);
		mpz_t temp;
		mpz_init(temp);
		if (!mpz_invert(temp, a.get_mpz_t(), prime.get_mpz_t()))
			mpz_set_si(temp, 0);
		mpz_class temp_class(temp);
		mpz_clear(temp);
		r = temp_class;
		return r;
	}

};

#endif 
