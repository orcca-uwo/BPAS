
#ifndef _SMALL_PRIME_FIELD_H
#define _SMALL_PRIME_FIELD_H

#include "BPASFiniteField.hpp"
#include <iostream>
#include <sstream>
#include <math.h>

using std::cout;
using std::endl;

//forward declarations
class Integer;
class RationalNumber;
class ComplexRationalNumber;
class BigPrimeField;
class GeneralizedFermatPrimeField;
class DenseUnivariateIntegerPolynomial;
class DenseUnivariateRationalPolynomial;
template <class Ring>
class SparseUnivariatePolynomial;

//NOTE the below code is not properly formatted with the corresponding
//implementation file. Only the montegomery version is setup properly.
#if ( 0 )
// simple small prime field using long long int to represent the number
class SmallPrimeField : public BPASField{

private:

 	static long int prime;
 	long long int a;

public:

	static mpz_class characteristic;
	static bool isPrimeField;
	static bool isSmallPrimeField;
	static bool isComplexField;


  	SmallPrimeField () : a(0) {}
  	SmallPrimeField (long long int _a) {
    	while (_a < 0){
      		_a = _a + prime;
    	}
    	a = _a%prime;
  	}
  	SmallPrimeField (const SmallPrimeField& c) {
    	a = c.a;
  	}
	//SmallPrimeField (Integer c);
	SmallPrimeField (RationalNumber c);
	SmallPrimeField (BigPrimeField c);
	SmallPrimeField (GeneralizedFermatPrimeField c);
	//SmallPrimeField (DenseUnivariateIntegerPolynomial c);
	//SmallPrimeField (DenseUnivariateRationalPolynomial c);
	//SmallPrimeField (SparseUnivariatePolynomial<Integer> c);
	//SmallPrimeField (SparseUnivariatePolynomial<RationalNumber> c);
	//SmallPrimeField (SparseUnivariatePolynomial<ComplexRationalNumber> c);
	//template <class Ring>
	//SmallPrimeField (SparseUnivariatePolynomial<Ring> c);

	SmallPrimeField* SPFpointer(SmallPrimeField* b) {
		return b;
	}
	SmallPrimeField* SPFpointer(RationalNumber* a) {
		std::cout << "BPAS error, try to cast pointer to Rational Number to pointer to SmallPrimeField" << std::endl;
		exit(1);
	}
	SmallPrimeField* SPFpointer(BigPrimeField* a) {
		std::cout << "BPAS error, try to cast pointer to BigPrimeField to pointer to SmallPrimeField" << std::endl;
		exit(1);
	}
	SmallPrimeField* SPFpointer(GeneralizedFermatPrimeField* a) {
		std::cout << "BPAS error, try to cast pointer to GeneralizedFermatPrimeField to pointer to SmallPrimeField" << std::endl;
		exit(1);
	}

  	static void setPrime(long int p){

  		prime = p;
  		//using string for mpz_class assignment
  		std::ostringstream ss;
		ss << p;
		std::string sp = ss.str();

		//mpz_t z;
		//mpz_init(z

		mpz_class w(sp);
	//	mpz_class w(p);
  		characteristic = w;
  	}

  	long int Prime(){
  		return prime;
  	}


  	void whichprimefield(){
		cout << "SmallPrimeField" << endl;
	}

	long long int number(){
		return a;
	}

//fing the nth primitive root of unity
	static SmallPrimeField findPrimitiveRootofUnity(long int n){
		if ( ((prime - 1) % n != 0)){
			cout << "ERROR: n does not divide prime - 1." << endl;
        	return -1;
    	}
    	bool flag = false;
    	long int  q = (prime - 1) / n;
    	SmallPrimeField p1 (prime - 1);
    	SmallPrimeField c;
    	int i = 0;
    	long int test = q * n / 2;
    	srand (time(NULL));
    	while(i < 20){
        	c =  rand();
            if ((c^test) == p1) {
                flag = true;
                return (c^q);
            }
        	i ++;
    	}
    	if (!flag ){
        	cout << "No primitive root found!"<< endl;
        	return 0;
    	}
    // If no primitive root found
	}

  	SmallPrimeField& operator= (SmallPrimeField c) {
    	if (this != &c) {
      	a = c.a;
    	}
    	return *this;
  	}
  	SmallPrimeField& operator= (long long int k) {
    	while (k < 0){
     	k = k + prime;
    	}
    	a = k%prime;
    	return *this;
  	}

  	inline bool isZero() {
    	return (a == 0);
  	}
  	inline void zero() {
    	a = 0;
  	}
  	inline bool isOne() {
    	return (a == 1);
  	}
  	inline void one() {
   		a = 1;
 	}

  	inline bool isNegativeOne() {
    	return (a == (prime - 1));
  	}
  	inline void negativeOne() {
    	a = prime - 1;
  	}

  	inline int isConstant() {
    	if (a >= 0)
    	return 1;
    	else { return -1; }
  	}


   	inline SmallPrimeField gcd (SmallPrimeField c) {
    	SmallPrimeField q (0);
     	SmallPrimeField r (0);
     	SmallPrimeField b (*this);
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
   	inline SmallPrimeField gcd (long long int c){
   		SmallPrimeField b (c);
   		SmallPrimeField r (*this);
   		return r.gcd(b);
   	}

  	inline bool operator== (SmallPrimeField& c) {
    	if (a == c.a)
      		return 1;
    	else { return 0; }
  	}

  	inline bool operator== (long long int k){
    	SmallPrimeField r (*this);
    	SmallPrimeField b (k);
    	if (b == r){
      	return 1;
    	}
    	else {
      	return 0;
    	}
 	}
  	inline bool operator!= (SmallPrimeField& c) {
    	if (a == c.a)
      		return 0;
    	else { return 1; }
  	}

  	inline bool operator!= (long long int k){
    	SmallPrimeField r (*this);
    	SmallPrimeField b (k);
    	if (b == r){
      	return 0;
    	}
    	else {
      	return 1;
    	}
  	}


  	inline SmallPrimeField operator+ (SmallPrimeField& c) {
    	SmallPrimeField r (*this);
    	return (r += c);
  	}
  	inline SmallPrimeField operator+ (long long int& c) {
    	SmallPrimeField r (*this);
    	SmallPrimeField b (c);
    	return (r += b);
  	}
  	inline SmallPrimeField operator+ (int& c) {
    	SmallPrimeField r (*this);
    	SmallPrimeField b (c);
    	return (r += b);
  	}
  	inline SmallPrimeField operator+= (long long int& c) {
    	SmallPrimeField r (*this);
    	SmallPrimeField b (c);
    	return (r += b);
  	}
  	inline SmallPrimeField& operator+= (SmallPrimeField c) {
    	a = (a + c.a)%prime;
    	return *this;
  	}

  	inline SmallPrimeField operator- (SmallPrimeField& c) {
    	SmallPrimeField r (*this);
    	return (r -= c);
  	}
  	inline SmallPrimeField operator- (long long int& c) {
    	SmallPrimeField r (*this);
    	SmallPrimeField b (c);
    	return (r -= b);
  	}
  	inline SmallPrimeField operator-= (long long int& c) {
    	SmallPrimeField r (*this);
    	SmallPrimeField b (c);
    	return (r -= b);
  	}
  	inline SmallPrimeField& operator-= (SmallPrimeField c) {
    	if ((a - c.a)<0){
      	a = prime+(a - c.a);
    	}
    	else{
      	a = a - c.a;
   		}
    	return *this;
  	}
  	inline SmallPrimeField operator- () {
    	a = prime - a;
    	return *this;
  	}
	inline SmallPrimeField operator* (SmallPrimeField& c) {
		SmallPrimeField r (*this);
		return (r *= c);
	}
  	inline SmallPrimeField operator* (long long int c) {
    	SmallPrimeField r (*this);
    	return (r *= c);
  	}

  	inline SmallPrimeField& operator*= (SmallPrimeField c) {
    	a = (a * c.a)%prime;
    	return *this;
  	}
  	inline SmallPrimeField operator*= (long long int& c) {
      	a = a * c;
      	SmallPrimeField r (a);
    	return r;
  	}
  void egcd (long long int x, long long int y, long long int *ao, long long int *bo, long long int *vo, long long int P){
  		long long int t, A, B, C, D, u, v, q;

		u = y;
		v = x;
		A = 1;
		B = 0;
		C = 0;
		D = 1;

  	do
  	{
  		q = u / v;
  		t = u;
  		u = v;
  		v = t - q * v;
  		t = A;
  		A = B;
  		B = t - q * B;
  		t = C;
  		C = D;
  		D = t - q * D;
  	}
  	while (v != 0);

  	*ao = A;
  	*bo = C;
  	*vo = u;
  }

  inline SmallPrimeField inverse2(){
    long long int n, b, v;
  	egcd (a, prime, &n, &b, &v, prime);
  	if (b < 0)
  		b += prime;
  	SmallPrimeField r (b);
  	return r;
  }


//Binary inverse

  SmallPrimeField inverse(){
  	long long int u = a;
  	long long int v = prime;
  	long long int x1 = 1;
  	long long int x2 = 0;
  	SmallPrimeField r(0);
  	while(u!=1 && v!=1){
  		while (u%2 == 0){
  			u = u >> 1;
  			if(x1%2==0)
  				x1 = x1 >> 1;
  			else
  				x1 = (x1 + prime) >> 1;
  		}
  		while (v%2 ==0){
  			v = v >> 1;
  			if(x2%2 == 0)
  				x2 = x2 >> 1;
  			else
  				x2 = (x2 + prime) >> 1;
  		}
  		if(u>=v){
  			u = u - v;
  			x1 = x1 - x2;
  		}
  		else{
  			v = v - u;
  			x2 = x2 - x1;
  		}
  	}
  	if(u == 1)
  		r = x1;
  	else
  		r = x2;
  	return r;
  }

  inline SmallPrimeField operator/ (SmallPrimeField& c) {
    SmallPrimeField r (*this);
    return (r /= c);
  }
  inline SmallPrimeField operator/ (long long int c) {
    SmallPrimeField r (*this);
    SmallPrimeField b (c);
    return (r /= b);
  }
  inline SmallPrimeField operator/= (long long int c) {
    SmallPrimeField r (*this);
    SmallPrimeField b (c);
    return (r /= b);
  }
  inline SmallPrimeField operator/= (SmallPrimeField c) {
    if (c.isZero()) {
      std::cout << "BPAS: error, dividend is zero from SmallPrimeField."<< std::endl;
            exit(1);
    }
    c=c.inverse();
    SmallPrimeField r (*this);
    r*=c;
    a = r.a;
    return *this;
  }

  inline SmallPrimeField operator^ (long long int e){
    SmallPrimeField r;
    if (isZero() || isOne() || e == 1)
      r = *this;
    else if (e == 2) {
      r = *this * *this;
    }
    else if (e > 2) {
      SmallPrimeField x (*this);
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
      r=r.inverse();
    }
    return r;
  }


  inline friend std::ostream& operator<< (std::ostream &out, SmallPrimeField c){
  	out << c.a;
  	return out;
  }
};


#else

/**
 * A prime field whose prime is 32 bits or less. Elements of this field
 * are encoded using montgomery trick.
 */
class SmallPrimeField : public BPASFiniteField<SmallPrimeField> {

private:

	long long int a;
	static long long int prime;
	//static long long int R;
	static unsigned long long int Pp;

//R = 2^Wt;
// 	static long long int Rp;

//Pp = -p^-1 mod R;
// 	static long long int N;

public:

	static mpz_class characteristic;
// static bool isPrimeField;
// static bool isSmallPrimeField;
// static bool isComplexField;

	SmallPrimeField ();

// using aR mod p to represent the number
	SmallPrimeField (long long int _a);

	SmallPrimeField (const SmallPrimeField& c);

	explicit SmallPrimeField (const Integer& c);

	explicit SmallPrimeField (const RationalNumber& c);

	explicit SmallPrimeField (const ComplexRationalNumber& c);

	explicit SmallPrimeField (const BigPrimeField& c);

	explicit SmallPrimeField (const GeneralizedFermatPrimeField& c);

	explicit SmallPrimeField (const DenseUnivariateIntegerPolynomial& c);

	explicit SmallPrimeField (const DenseUnivariateRationalPolynomial& c);

	explicit SmallPrimeField (const SparseUnivariatePolynomial<Integer>& c);

	explicit SmallPrimeField (const SparseUnivariatePolynomial<RationalNumber>& c);

	explicit SmallPrimeField (const SparseUnivariatePolynomial<ComplexRationalNumber>& c);

	template <class Ring>
	explicit SmallPrimeField (const SparseUnivariatePolynomial<Ring>& c);

	SmallPrimeField* SPFpointer(SmallPrimeField* b);

	SmallPrimeField* SPFpointer(RationalNumber* a);

	SmallPrimeField* SPFpointer(BigPrimeField* a);

	SmallPrimeField* SPFpointer(GeneralizedFermatPrimeField* a);

// return the actual number
	long long int number() const;

	void whichprimefield();

	static void setPrime(long long int p){

		prime = p;
//using string for mpz_class assignment
		std::ostringstream ss;
		ss << p;
		std::string sp = ss.str();

//mpz_t z;
//mpz_init(z

		mpz_class w(sp);
//	mpz_class w(p);
  		SmallPrimeField::characteristic = w;

		// long long int p1 = R - prime;

		// for(int i = 0; i < R; i ++){
		// 	if ((i * p1)%R == 1){
		// 		Pp = i;
		// 		break;
		// 	}
		// }
//cout << Pp << endl;
      long long int t, A, B, C, D, u, v,q;
    __asm__(                          //compute the first step of egcd where R = q*prime + v
       // "movq %%rax,%%rsi\n\t"
        "xor %%rax,%%rax\n\t"
        "movq $1,%%rdx\n\t"
        "divq %2\n\t"
       // "movq %%rsi,%%rax\n\t"
       // "movq %%rdx,%0\n\t"
       // "movq %%rax,%1\n\t"
        : "=&d" (v),"=&a" (q)
        : "b"(prime)
        :"rsi","rdi");
    A = 1;
    B = 0;
    C = 0;
    D = 1;
   // q = R/prime;
    u = prime;
    //v = R - R mod prime;
    t = A;
    A = B;
    B = t - q * B;
    t = C;
    C = D;
    D = t - q * D;

    while (v != 0){
      q = u / v;
      t = u;
      u = v;
      v = t - q * v;
      t = A;
      A = B;
      B = t - q * B;
      t = C;
      C = D;
      D = t - q * D;
   //   printf("%llu,%llu\n",A,C);
    }
     //unsigned long long int res;
    if(C < 0){
      C = 0-C;
      Pp = (unsigned long long int)C;
    }
    else{
      __asm__(                          //compute -C mod R;
       // "movq %%rax,%%rsi\n\t"
        "xor %%rax,%%rax\n\t"
        "movq $1,%%rdx\n\t"
        "sub %1,%%rax\n\t"
        "sbb $0,%%rdx\n\t"
       // "movq %%rsi,%%rax\n\t"
       // "movq %%rdx,%0\n\t"
        "movq %%rax,%%rdx\n\t"
        : "=&d" (Pp)
        : "b"((unsigned long long int)C)
        :"rsi","rdi");
    }

	}

  mpz_class getCharacteristic() const override {
    return characteristic;
  }

	long long int Prime();

	SmallPrimeField& operator= (const SmallPrimeField& c);

	SmallPrimeField& operator= (long long int k);

	SmallPrimeField findPrimitiveRootOfUnity(long int n) const {
		return SmallPrimeField::findPrimitiveRootofUnity(n);
	}


	static SmallPrimeField findPrimitiveRootofUnity(long long int n){
		if ( ((prime - 1) % n != 0)){
			cout << "ERROR: n does not divide prime - 1." << endl;
			return -1;
		}
		bool flag = false;
		long long int  q = (prime - 1) / n;
		SmallPrimeField p1 (prime - 1);
		SmallPrimeField c;
		int i = 0;
		long long int test = q * n / 2;
		srand (time(NULL));
		while(i < 20){
			c =  rand();
			if ((c^test) == p1) {
				flag = true;
				return (c^q);
			}
			i ++;
		}
		if (!flag ){
			cout << "No primitive root found!"<< endl;
			return 0;
		}
// If no primitive root found
    return 0;
	}

	inline bool isZero() const {
		return (a == 0);
	}

	inline void zero() {
		a = 0;
	}

	inline bool isOne() const {
    SmallPrimeField b(1);
		return  (a == b.a);
	}

	inline void one() {
    SmallPrimeField b(1);
		a = b.a;// R % prime;
	}

	inline bool isNegativeOne() {
    SmallPrimeField b(-1);
   	return (a == b.a);
	}
	inline void negativeOne() {
    SmallPrimeField b(-1);
 	  a = b.a;
	}

	inline int isConstant() {
		if (a >= 0)
			return 1;
		else { return -1; }
	}

  SmallPrimeField unitCanonical(SmallPrimeField* u = NULL, SmallPrimeField* v = NULL) const ;

	inline SmallPrimeField operator+ (const SmallPrimeField& c) const {
		SmallPrimeField r (*this);
		return (r += c);
	}

	inline SmallPrimeField operator+ (const long long int& c) const {
		SmallPrimeField r (*this);
		SmallPrimeField b (c);
		return (r += b);
	}
	inline SmallPrimeField operator+ (const int& c) const {
		SmallPrimeField r (*this);
		SmallPrimeField b (c);
		return (r += b);
	}
	inline SmallPrimeField operator+= (const long long int& c) {
		SmallPrimeField r (*this);
		SmallPrimeField b (c);
		return (r += b);
	}

	inline SmallPrimeField& operator+= (const SmallPrimeField& c) {
		//a = (a + c.a)%prime;
		//return *this;
    a = a + c.a;
    a -= prime;
    a += (a >> 63) & prime;
    return *this;
	}

	inline SmallPrimeField operator- (const SmallPrimeField& c) const {
		SmallPrimeField r (*this);
		return (r -= c);
	}

	inline SmallPrimeField operator- (const long long int& c) const {
		SmallPrimeField r (*this);
		SmallPrimeField b (c);
		return (r -= b);
	}
	inline SmallPrimeField operator-= (const long long int& c) {
		SmallPrimeField r (*this);
		SmallPrimeField b (c);
		return (r -= b);
	}

	inline SmallPrimeField& operator-= (const SmallPrimeField& c) {
		// if ((a - c.a)<0){
		// 	a = prime+(a - c.a);
		// }
		// else{
		// 	a = a - c.a;
		// }
		// return *this;
    a = a - c.a;
    a += (a >> 63) & prime;
    return *this;
	}
	inline SmallPrimeField operator- () const {
		SmallPrimeField ret (*this);
		ret.a = prime - a;
		return ret;
	}

	inline SmallPrimeField operator* (const SmallPrimeField& c) const {
 		SmallPrimeField r (*this);

		r*=c;
//cout << r << endl;
		return r;
	}

	inline SmallPrimeField operator* (long long int c) const {
		SmallPrimeField r (*this);
		SmallPrimeField b(c);
		return (r *= b);
	}
// Mont(a,b) = ab/R mod p
// Mont(aR,bR) = abR mod p
	inline SmallPrimeField& operator*= (const SmallPrimeField& c) {
		// long long int x = (c.a * a);
		// long long int w = x*Pp;
		// long long int y = x + prime*(w&(R-1));
		// long long int z =  y>>32 ;
		// if(z >= prime){
		// 	z = z - prime;
		// }
		// a = z;
		// return *this;
    long long int _a = a;
  __asm__(
      "mulq %2\n\t"
      "movq %%rax,%%rsi\n\t"
      "movq %%rdx,%%rdi\n\t"
      "imulq %3,%%rax\n\t"
      "mulq %4\n\t"
      "add %%rsi,%%rax\n\t"
      "adc %%rdi,%%rdx\n\t"
      "subq %4,%%rdx\n\t"
      "mov %%rdx,%%rax\n\t"
      "sar $63,%%rax\n\t"
      "andq %4,%%rax\n\t"
      "addq %%rax,%%rdx\n\t"
      : "=&d" (_a)
      : "a"(_a),"rm"(c.a),"b"((unsigned long long int)Pp),"c"(prime)
      :"rsi","rdi");
  a = _a;
    return *this;
	}

//partial Montgomery inversion
	long long int* pinverse();
//  long long int MontMul(long long int a, long long int b);

//return bc/R mod p;
	static long long int Mont(long long int b, long long int c);
  static long long int getRsquare();
//Montgomery inversion
	inline SmallPrimeField inverse() const {
		 SmallPrimeField r(*this);
	//	 SmallPrimeField rsquare(R);
     long long int rsquare = getRsquare();
		// long long int* result = r.pinverse();
		// if(result[0] < 32){
		// 	result[0] = Mont(result[0], rsquare.a);
		// 	result[1] += 32;
		// }
		// result[0] = Mont(result[0], rsquare.a);
		// result[0] = Mont(result[0], pow(2,64-result[1]));
		// r.a = result[0];
		// return r;

      if (r.a == 0)
      {
        printf("Not invertible!\n");
        return 0;
      }
   // long long int rsquare =  getRsquare(prime);
    long long int* result = r.pinverse();
    if(result[1] < 64){
      result[0] = Mont(result[0], rsquare);
      result[1] += 64;
    }
    result[0] = Mont(result[0], rsquare);
    long long int tmp;
    if(result[1] != 64){
      r.a = 1L << (128 - result[1]);
      result[0] = Mont(result[0], r.a);
    }
    r.a = result[0];
    free(result);
    return r;
	}

//Binary inverse
	SmallPrimeField inverse2();

	inline SmallPrimeField operator^ (long long int e) const {
		SmallPrimeField r;
		if (isZero() || isOne() || e == 1)
			r = *this;
		else if (e == 2) {
			r = *this * *this;
		}
		else if (e > 2) {
			SmallPrimeField x (*this);
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
			r=r.inverse();
		}
		return r;

	}

	inline SmallPrimeField& operator^= (long long int e) {
		*this = *this ^ e;
		return *this;
	}

	inline bool operator== (const SmallPrimeField& c) const {
		if ((*this).number() == c.number())
			return 1;
		else { return 0; }
	}

	inline bool operator== (long long int k) const {
		SmallPrimeField r (*this);
		SmallPrimeField b (k);
		if (b.number() == r.number()){
			return 1;
		}
		else {
			return 0;
		}
	}

	inline bool operator!= (const SmallPrimeField& c) const {
		if ((*this).number() == c.number())
			return 0;
		else { return 1; }
	}


	inline bool operator!= (long long int k) const {
		SmallPrimeField r (*this);
		SmallPrimeField b (k);
		if (b.number() == r.number()){
			return 0;
		}
		else {
			return 1;
		}
	}

	inline ExpressionTree convertToExpressionTree() const {
		return ExpressionTree(new ExprTreeNode(this->number()));
	}

	inline SmallPrimeField operator/ (const SmallPrimeField& c) const {
		SmallPrimeField r (*this);
		r /= c;
		return r;
	}

	inline SmallPrimeField operator/ (long long int c) const {
		SmallPrimeField r (*this);
		SmallPrimeField b(c);
		return (r /= b);
	}

	inline SmallPrimeField& operator/= (const SmallPrimeField& c) {
		if (c.isZero()) {
			std::cout << "BPAS: error, dividend is zero from SmallPrimeField."<< std::endl;
			exit(1);
		}
    //cout << "division" << endl;
		SmallPrimeField inv = c.inverse();
		SmallPrimeField r (*this);
		r *= inv;
    (*this).a = r.a;
		return *this;
	}

  inline SmallPrimeField operator% (const SmallPrimeField& c) const {
    return *this % c;
  }

  inline SmallPrimeField& operator%= (const SmallPrimeField& c) {
    (*this).a  = (*this).a % c.a;
    return *this;
  }

	inline SmallPrimeField gcd (const SmallPrimeField& other) const {
		SmallPrimeField q (0);
		SmallPrimeField r (0);
		SmallPrimeField b (*this);
		SmallPrimeField c (other);
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

	/**
	 * Compute squarefree factorization of *this
	 */
	inline Factors<SmallPrimeField> squareFree() const {
		std::vector<SmallPrimeField> ret;
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
	SmallPrimeField euclideanDivision(const SmallPrimeField& b, SmallPrimeField* q = NULL) const;

	/**
	 * Perform the extended euclidean division on *this and b.
	 * Returns the GCD. If s and t are not NULL, returns the bezout coefficients in them.
	 */
	SmallPrimeField extendedEuclidean(const SmallPrimeField& b, SmallPrimeField* s = NULL, SmallPrimeField* t = NULL) const;

	/**
	 * Get the quotient of *this and b.
	 */
	SmallPrimeField quotient(const SmallPrimeField& b) const;

	/**
	 * Get the remainder of *this and b.
	 */
	SmallPrimeField remainder(const SmallPrimeField& b) const;

};

#endif //montgomery switch

#endif //include guard
