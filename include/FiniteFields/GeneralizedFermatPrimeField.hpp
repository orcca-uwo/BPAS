
#ifndef _GENERALIZED_FERMAT_PRIMT_FIELD_H_
#define _GENERALIZED_FERMAT_PRIMT_FIELD_H_

#include "BPASFiniteField.hpp"
#include <iostream>

typedef unsigned long long int usfixn64;

//ULMAX=2^64-1 //equivalent to -1 on 64bits machines
#define ULMAX 18446744073709551615ULL
#define SQRTR 3037000502
#define RC 9223372019674906624ULL

//forward declarations
class Integer;
class RationalNumber;
class ComplexRationalNumber;
class SmallPrimeField;
class BigPrimeField;
class DenseUnivariateIntegerPolynomial;
class DenseUnivariateRationalPolynomial;
template <class Ring>
class SparseUnivariatePolynomial;

using std::endl;
using std::cout;

// field when the prime is a generalized fermat prime:
// r = (2^w +- 2^u)
// prime = r^k
/**
 * A finite field whose prime should be a generalized fermat number.
 * That is, for r = (2^w +/- 2^u), the prime is r^k, for some k.
 */
class GeneralizedFermatPrimeField : public BPASFiniteField<GeneralizedFermatPrimeField> {

private:
	static mpz_class& prime;

public:
	usfixn64 *x;
	// number = xk-1*r^k-1 + xk-2*r^k-2 + ... + x1*r + x0
	static mpz_class characteristic;
	static RingProperties properties;

	static usfixn64 r;
	static int k;
	// static bool isPrimeField;
	// static bool isSmallPrimeField;
	// static bool isComplexField;

	GeneralizedFermatPrimeField ();

	GeneralizedFermatPrimeField (mpz_class a);

	GeneralizedFermatPrimeField (int a);

	GeneralizedFermatPrimeField (const GeneralizedFermatPrimeField& c);

	explicit GeneralizedFermatPrimeField (const Integer& c);

	explicit GeneralizedFermatPrimeField (const RationalNumber& c);

	explicit GeneralizedFermatPrimeField (const SmallPrimeField& c); //? NOT SURE

	explicit GeneralizedFermatPrimeField (const BigPrimeField& c);

	explicit GeneralizedFermatPrimeField (const DenseUnivariateIntegerPolynomial& c);

	explicit GeneralizedFermatPrimeField (const DenseUnivariateRationalPolynomial& c);

	explicit GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<Integer>& c);

	explicit GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<RationalNumber>& c);

	explicit GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<ComplexRationalNumber>& c);

	template <class Ring>
	explicit GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<Ring>& c);

	GeneralizedFermatPrimeField* GPFpointer(GeneralizedFermatPrimeField* a);

	GeneralizedFermatPrimeField* GPFpointer(RationalNumber* a);

	GeneralizedFermatPrimeField* GPFpointer(SmallPrimeField* a);

	GeneralizedFermatPrimeField* GPFpointer(BigPrimeField* a);

	static void setPrime (mpz_class p,usfixn64 R, int K){
		prime = p;
	//   characteristic = p;
		r = R;
		k = K;
	}

	void setX (mpz_class a);

	mpz_class Prime() const;

	mpz_class number() const;

	GeneralizedFermatPrimeField unitCanonical(GeneralizedFermatPrimeField* u = NULL, GeneralizedFermatPrimeField* v = NULL) const;

	static mpz_class power(mpz_class xi, mpz_class yi) {
		mpz_class res = 1;     // Initialize result

		xi = xi % prime; // Update x if it is more than or
		// equal to p

		while (yi > 0){
		// If y is odd, multiply x with result
			if ((yi % 2) == 1){
				res = (res*xi) % prime;
				yi = yi - 1;
			}

		// y must be even now
		yi = yi / 2; // y = y/2
		xi = (xi*xi) % prime;
		}
		return res;
	}

	GeneralizedFermatPrimeField findPrimitiveRootOfUnity(long int n) const {
		return GeneralizedFermatPrimeField::findPrimitiveRootofUnity(n);
	}

	static mpz_class findPrimitiveRootofUnity_plain(mpz_class n){
		if ( ((prime - 1) % n != 0)){
			cout << "ERROR: n does not divide prime - 1." << endl;
			return -1;
		}
		bool flag = false;
		mpz_class  q = (prime - 1) / n;
		mpz_class p1 = prime - 1;
		mpz_class c;
		int i = 0;
		mpz_class test = q * n / 2;


		srand (time(NULL));

		while(i < 20){   
	//	std::cout << "i " << i <<  std::endl; 	
			c =  rand();
	//	cout << "c " << c << endl;	
			if (power(c,test) == p1) {
				flag = true;
				return  power(c,q);
			}
			i++;

		}
		if (!flag ){
			cout << "No primitive root found!"<< endl;
			return 0;
		}
		return 0;
	// If no primitive root found
	}

	// find a privimitive root of unity that w^(N/2k) = r;
	static GeneralizedFermatPrimeField findPrimitiveRootofUnity(mpz_class n){
		if(n<(2*k)){
			throw std::invalid_argument( "findPrimitiveRootofUnity error: N is less than 2k");
		}
		mpz_class g = findPrimitiveRootofUnity_plain(n);
		mpz_class n2k = n / (2 * k);
		mpz_class a = power(g, n2k);
		mpz_class b = a;
		std::stringstream str;
		str << r;
		mpz_class r_mpz (str.str());
		int j = 1;
		while (b != r_mpz){
			b = (a * b)%prime;
			j++;
		}
		GeneralizedFermatPrimeField result(power(g,j));
		return result;
	}

	GeneralizedFermatPrimeField& operator= (const GeneralizedFermatPrimeField& c);

	GeneralizedFermatPrimeField& operator= (const mpz_class& c);

	GeneralizedFermatPrimeField& operator= (int c);

	inline bool isZero() const {
		for (int i = 0; i < k; i++) {
			if (x[i] != 0) {
				return 0;
			}
		}
		return 1;
	}

	inline void zero() {
		memset(x, 0x00, (k) * sizeof(usfixn64));
	}

	inline bool isOne() const {
		for (int i = k - 1;i > 0; i --){
			if (x[i]!=0){
				return 0;
			}
		}

		if (x[0] == 1){
			return 1;
		}
		else {return 0;}
	}

	inline void one() {
		for (int i = k - 1;i > 0; i --){
			x[i] = 0;
		}

		x[0] = 1;
	}

	inline bool isNegativeOne() const {
		for (int i = 0; i < (k - 1);i ++){
			if (x[i] != 0){
				return 0;
			}
		}
		if (x[k - 1] != r){
			return 0;
		}
		else {return 1;}
	}

	inline void negativeOne() {
		for (int i = 0; i < (k - 1); i ++){
			x[i] = 0;
		}
		x[k - 1] = r;
	}

	inline int isConstant() const {
		return 1;
	}

	inline bool operator== (const GeneralizedFermatPrimeField& c) const {
		for (int i = 0; i < k; i++) {
			if (x[i] != c.x[i]) {
				return 0;
			}
		}

		return 1;
	}

	inline bool operator== (const mpz_class& c) const {
		GeneralizedFermatPrimeField b (c);
		for (int i = 0; i < k; i++) {
			if (x[i] != b.x[i]) {
				return 0;
			}
		}

		return 1;
	}

	inline bool operator!= (const GeneralizedFermatPrimeField& c) const {
		for (int i = 0; i < k; i++) {
			if (x[i] != c.x[i]) {
				return 1;
			}
		}
		return 0;
	}

	inline bool operator!= (const mpz_class& c) const {
		GeneralizedFermatPrimeField b (c);
		for (int i = 0; i < k; i++) {
			if (x[i] != b.x[i]) {
				return 1;
			}
		}

		return 0;
	}

	inline GeneralizedFermatPrimeField operator+ (const GeneralizedFermatPrimeField& c) const {
		GeneralizedFermatPrimeField t (*this);
		return (t += c);
	}

	inline GeneralizedFermatPrimeField operator+ (const mpz_class& c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t += b);
	}

	inline GeneralizedFermatPrimeField& operator+= (const mpz_class& c) {
		GeneralizedFermatPrimeField b (c);
		*this += b;
		return *this;
	}

	inline GeneralizedFermatPrimeField operator+ (int c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t += b);
	}

	inline GeneralizedFermatPrimeField& operator+= (int c) {
		GeneralizedFermatPrimeField b (c);
		*this += b;
		return *this;
	}

	// computer zi = xi + yi for i in 0...k-1
	// let zk = 0
	// for i in 0...k-1, zi/r = qi*a + si
	// zi+1 = zi+1 + qi
	// if zk ==0 return
	// if zk == 1 return (r,0,0...0)
	GeneralizedFermatPrimeField& operator+= (const GeneralizedFermatPrimeField& y);

	inline GeneralizedFermatPrimeField operator- (const GeneralizedFermatPrimeField& c) const {
		GeneralizedFermatPrimeField t (*this);
		return (t -= c);
	}

	inline GeneralizedFermatPrimeField operator- (const mpz_class& c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t -= b);
	}

	//TODO this and the next -= seems wrong
	inline GeneralizedFermatPrimeField& operator-= (const mpz_class& c) {
		GeneralizedFermatPrimeField b (c);
		*this -= b;
		return *this;
	}

	inline GeneralizedFermatPrimeField operator- (int c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t -= b);
	}

	inline GeneralizedFermatPrimeField& operator-= (int c) {
		GeneralizedFermatPrimeField b (c);
		*this -= b;
		return *this;
	}

	// using similar algorithm as addition
	GeneralizedFermatPrimeField& operator-= (const GeneralizedFermatPrimeField& y);

	inline GeneralizedFermatPrimeField operator- () const {
		GeneralizedFermatPrimeField b (0);
		return (b - *this);
	}

	// part of multiplication
	void smallAdd2 (usfixn64 *xm, usfixn64* ym, short & c);

	void oneShiftRight (usfixn64 * xs);

	// x*y = s0 + s1*r + s2 * r^2
	void mulLong_2 (usfixn64 x, usfixn64 y, usfixn64 &s0,usfixn64 &s1, usfixn64 &s2);

	// special one for prime3
	void mulLong_3 (usfixn64 const &x, usfixn64 const &y, usfixn64 &s0,
		usfixn64 &s1, usfixn64 &s2);

	void multiplication (usfixn64* __restrict__ xs, const usfixn64* __restrict__ ys,
		usfixn64 permutationStride, usfixn64* lVector,
		usfixn64 *hVector, usfixn64* cVector,
		usfixn64* lVectorSub,
		usfixn64 *hVectorSub, usfixn64* cVectorSub );

	void multiplication_step2 (usfixn64* __restrict__ xs, usfixn64 permutationStride,
		usfixn64* __restrict__ lVector,
		usfixn64 * __restrict__ hVector,
		usfixn64* __restrict__ cVector);

	inline GeneralizedFermatPrimeField operator* (const GeneralizedFermatPrimeField& c) const {
		GeneralizedFermatPrimeField t (*this);
		return (t *= c);
	}

	inline GeneralizedFermatPrimeField operator* (const mpz_class& c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t *= b);
	}

	inline GeneralizedFermatPrimeField& operator*= (const mpz_class& c) {
		GeneralizedFermatPrimeField b (c);
		*this *= b;
		return *this;
	}

	inline GeneralizedFermatPrimeField operator* (int c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t *= b);
	}

	inline GeneralizedFermatPrimeField& operator*= (int c) {
		GeneralizedFermatPrimeField b (c);
		*this *= c;
		return *this;
	}

	// using GMP multiplication 
	inline GeneralizedFermatPrimeField& operator*= (const GeneralizedFermatPrimeField& c){
		mpz_class xi = number();
		mpz_class yi = c.number();
		*this = xi*yi;
		return *this;
	}

	// special multiplication for P3
	//r = 9223372054034644992
	// k = 8
	GeneralizedFermatPrimeField MultiP3 (GeneralizedFermatPrimeField ys);

	//Multiple by some power of r
	//using shift only
	GeneralizedFermatPrimeField MulPowR(int s);


	void egcd (const mpz_class& x, const mpz_class& y, mpz_class *ao, mpz_class *bo, mpz_class *vo, mpz_class P);

	inline GeneralizedFermatPrimeField inverse2(){
		mpz_class a, n, b, v;
		GeneralizedFermatPrimeField t(*this);
		a = t.number();
		egcd (a, prime, &n, &b, &v, prime);
		if (b < 0)
			b += prime;
		GeneralizedFermatPrimeField R (b);
		return R;
	}

	inline GeneralizedFermatPrimeField operator^ (long long int c) const {
		GeneralizedFermatPrimeField t (*this);
		mpz_class b (std::to_string(c), 10);
		return (t ^ b);
	}

	inline GeneralizedFermatPrimeField operator^ (const mpz_class& exp) const {
		GeneralizedFermatPrimeField t;
		mpz_class e(exp);

		if (isZero() || isOne() || e == 1)
			t = *this;
		else if (e == 2) {
			t = *this * *this;
		}
		else if (e > 2) {
			GeneralizedFermatPrimeField x (*this);
			t.one();

			while (e != 0) {
				if ((e % 2) == 1)
					t = t * x;
				x = x * x;
				e >>= 1;
			}
		}
		else if (e == 0) {
			t.one();
		}
		else {
			t = *this ^ (-e);
			t=t.inverse();
		}
		return t;
	}

	inline GeneralizedFermatPrimeField& operator^= (long long int c) {
		*this = *this ^ c;
		return *this;
	}

	inline GeneralizedFermatPrimeField& operator^= (const mpz_class& c) {
		*this = *this ^ c;
		return *this;
	}

	ExpressionTree convertToExpressionTree() const {
		return ExpressionTree(new ExprTreeNode(this->number()));
	}

	inline GeneralizedFermatPrimeField operator/ (const GeneralizedFermatPrimeField& c) const {
		GeneralizedFermatPrimeField t (*this);
		return (t /= c);
	}

	inline GeneralizedFermatPrimeField operator/ (long int c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t /= b);
	}

	inline GeneralizedFermatPrimeField operator/ (const mpz_class& c) const {
		GeneralizedFermatPrimeField t (*this);
		GeneralizedFermatPrimeField b (c);
		return (t /= b);
	}

	inline GeneralizedFermatPrimeField& operator/= (const GeneralizedFermatPrimeField& c) {
		if (c.isZero()) {
			std::cout << "BPAS: error, dividend is zero from GeneralizedFermatPrimeField."<< std::endl;
			exit(1);
		}
		GeneralizedFermatPrimeField inv = c.inverse();
		*this *= c;
		return *this;
	}

	inline GeneralizedFermatPrimeField& operator/= (long int c) {
		GeneralizedFermatPrimeField b (c);
		*this /= b;
		return *this;
	}

	inline GeneralizedFermatPrimeField& operator/= (const mpz_class& c) {
		GeneralizedFermatPrimeField b (c);
		*this /= b;
		return *this;
	}
	
	inline GeneralizedFermatPrimeField operator% (const GeneralizedFermatPrimeField& c) const {
		return 0;
	}

	inline GeneralizedFermatPrimeField& operator%= (const GeneralizedFermatPrimeField& c) {
		*this = 0;
		return *this;
	}

	inline GeneralizedFermatPrimeField gcd(const GeneralizedFermatPrimeField& a) const {
		std::cerr << "GeneralizedFermatPrimeField::gcd NOT YET IMPLEMENTED" << std::endl;
		exit(1);
		return GeneralizedFermatPrimeField();
	}
	
	/**
	 * Compute squarefree factorization of *this
	 */
	inline Factors<GeneralizedFermatPrimeField> squareFree() const {
		std::vector<GeneralizedFermatPrimeField> ret;
		ret.push_back(*this);
		return ret;
	}

	/** 
	 * Get the euclidean size of *this.
	 */
	inline GeneralizedFermatPrimeField euclideanSize() const {
		return (*this).number();
	}
	
	/**
	 * Perform the eucldiean division of *this and b. Returns the 
	 * remainder. If q is not NULL, then returns the quotient in q. 
	 */ 
	GeneralizedFermatPrimeField euclideanDivision(const GeneralizedFermatPrimeField& b, GeneralizedFermatPrimeField* q = NULL) const;

	/**
	 * Perform the extended euclidean division on *this and b. 
	 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
	 */
	GeneralizedFermatPrimeField extendedEuclidean(const GeneralizedFermatPrimeField& b, GeneralizedFermatPrimeField* s = NULL, GeneralizedFermatPrimeField* t = NULL) const;
	
	/**
	 * Get the quotient of *this and b.
	 */
	GeneralizedFermatPrimeField quotient(const GeneralizedFermatPrimeField& b) const;

	/** 
	 * Get the remainder of *this and b.
	 */
	GeneralizedFermatPrimeField remainder(const GeneralizedFermatPrimeField& b) const;

	// GMP inverse
	GeneralizedFermatPrimeField inverse() const;


};

#endif
