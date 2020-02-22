
#include <sstream>
#include <iostream>
#include <math.h>

#include "Ring/Integer.hpp"
#include "Ring/RationalNumber.hpp"
#include "Ring/ComplexRationalNumber.hpp"
#include "FiniteFields/SmallPrimeField.hpp"
#include "FiniteFields/BigPrimeField.hpp"
#include "FiniteFields/GeneralizedFermatPrimeField.hpp"
#include "IntegerPolynomial/uzpolynomial.h"
#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h" 

using std::endl;
using std::cout;

mpz_class SmallPrimeField::characteristic(883949569);
RingProperties SmallPrimeField::properties = SMALL_PRIME_FIELD;

// bool SmallPrimeField::isPrimeField = 1;
// bool SmallPrimeField::isSmallPrimeField = 1;
// bool SmallPrimeField::isComplexField = 0;

long long int SmallPrimeField::prime = 883949569;
//long long int SmallPrimeField::R = 4294967296;
//Pp = -p^-1 mod R;
unsigned long long int SmallPrimeField::Pp = 2677397675937103871;


SmallPrimeField::SmallPrimeField () : a(0) {

}

// using aR mod p to represent the number
SmallPrimeField::SmallPrimeField (long long int _a) {
    __asm__(
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "idivq %2\n\t"
      "movq %%rdx,%%rax\n\t"
      "imulq %%rbx\n\t"
      "idivq %2\n\t"
      "movq %%rdx,%%rax\n\t"
      "sar $63,%%rax\n\t"
      "andq %2,%%rax\n\t"
      "addq %%rax,%%rdx\n\t"
      : "=&d" (_a)
      : "b"(_a),"c"(prime));
  // cout << _a << endl;
    a = _a;
}

SmallPrimeField::SmallPrimeField (const SmallPrimeField& c) {
  a = c.a;
}

SmallPrimeField::SmallPrimeField (const Integer& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}

SmallPrimeField::SmallPrimeField (const RationalNumber& c) {
	// need a case handling depending on the size of the components of c, so that they fit inside SmallPrimeField
        std::cout << "BPAS error, try to construct a SmallPrimeField from RationalNumber class." << std::endl;
        exit(1);
	/*long long int num, den;
	mpz_t z;
	mpz_init(z);
	mpz_set(z,c.get_num().get_mpz_t());
	mpz_export(&num, mpz_size(z), 1, sizeof(num), 0, 0, z);
	mpz_set(z,c.get_den().get_mpz_t());
	mpz_export(&den, &mpz_size(z), 1, sizeof(den), 0, 0, z);
	*this = (num/den);
	*/
	//*this = (long long int) c.get_num();
}

SmallPrimeField::SmallPrimeField (const ComplexRationalNumber& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}	

SmallPrimeField::SmallPrimeField (const BigPrimeField& c) {
	// need to write code to implement proper embedding
        std::cout << "BPAS error, try to construct a SmallPrimeField from BigPrimeField class." << std::endl;
        exit(1);
}

SmallPrimeField::SmallPrimeField (const GeneralizedFermatPrimeField& c) {
	// need to write code to implement proper embedding
        std::cout << "BPAS error, try to construct a SmallPrimeField from GeneralizedFermatPrimeField class." << std::endl;
        exit(1);
}

SmallPrimeField::SmallPrimeField (const DenseUnivariateIntegerPolynomial& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}

SmallPrimeField::SmallPrimeField (const DenseUnivariateRationalPolynomial& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}

SmallPrimeField::SmallPrimeField (const SparseUnivariatePolynomial<Integer>& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}

SmallPrimeField::SmallPrimeField (const SparseUnivariatePolynomial<RationalNumber>& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}

SmallPrimeField::SmallPrimeField (const SparseUnivariatePolynomial<ComplexRationalNumber>& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}

template <class Ring>
SmallPrimeField::SmallPrimeField (const SparseUnivariatePolynomial<Ring>& c) {
	std::cerr << "Cannot convert input to SmallPrimeField! " << std::endl;
	exit(1);
}	

SmallPrimeField* SmallPrimeField::SPFpointer(SmallPrimeField* b) {
	return b;
}

SmallPrimeField* SmallPrimeField::SPFpointer(RationalNumber* a) {
	std::cout << "BPAS error, try to cast pointer to Rational Number to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

SmallPrimeField* SmallPrimeField::SPFpointer(BigPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to BigPrimeField to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

SmallPrimeField* SmallPrimeField::SPFpointer(GeneralizedFermatPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to GeneralizedFermatPrimeField to pointer to SmallPrimeField" << std::endl;
	exit(1);
}


// return the actual number
long long int SmallPrimeField::number() const {
	// long long int x = a;
	// long long int w = x*Pp;
	// long long int y = x + prime*(w & (R-1));
	// long long int z =  y >> 32 ;
	// if(z >= prime){
	// 	z = z - prime;
	// }
	// return z;
	long long int z = SmallPrimeField::Mont(a,1);
	return z;
}

void SmallPrimeField::whichprimefield(){
	cout << "MontSmallPrimeField" << endl;
}
	
long long int SmallPrimeField::Prime(){
	return prime;
}

SmallPrimeField& SmallPrimeField::operator= (const SmallPrimeField& c) {
  	a = c.a;
  	return *this;
}

SmallPrimeField& SmallPrimeField::operator= (long long int k) {
	__asm__(
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "idivq %2\n\t"
      "movq %%rdx,%%rax\n\t"
      "imulq %%rbx\n\t"
      "idivq %2\n\t"
      "movq %%rdx,%%rax\n\t"
      "sar $63,%%rax\n\t"
      "andq %2,%%rax\n\t"
      "addq %%rax,%%rdx\n\t"
      : "=&d" (k)
      : "b"(k),"c"(prime));
	a = k;
	return *this;
}

SmallPrimeField SmallPrimeField::unitCanonical(SmallPrimeField* u , SmallPrimeField* v) const {
	if (isZero()) {
		if (u != NULL) {
			*u = 1;
		}
		if (v != NULL) {
			*v = 1;
		}
		return SmallPrimeField(0);
	} else {
		if (u != NULL) {
			*u = *this;
		}
		if (v != NULL) {
			*v = this->inverse();
		}
		return SmallPrimeField(1);
	}
}

//partial Montgomery inversion
long long int* SmallPrimeField::pinverse(){
    long long int u = a;
    long long int v = prime;
    long long int x2 = 0;
    long long int* r = (long long int *)malloc(2*sizeof(long long int));
    r[0] = 1;
    r[1] = 0;
    while(v > 0){
	if (v%2 == 0){
	    v = v >> 1;
	    r[0] = 2 * r[0];
	}
	else if (u%2 ==0){
	    u = u >> 1;
	    x2 = 2 * x2;
	}
	else if(v>=u){
	    v = (v - u) >> 1;
	    x2 = x2 + r[0];
	    r[0] = 2 * r[0];
	}
	else{
	    u = (u - v) >> 1;
	    r[0] = x2 + r[0];
	    x2 = 2 * x2;
	}
	r[1]++;
    }
    if(u != 1){
		printf("Not invertible!\n");
		return NULL;
    }
    if(r[0] > prime){
		r[0] -= prime;    
    }
    return r;
}

//return bc/R mod p;
long long int SmallPrimeField::Mont(long long int b, long long int c){
	// long long int x = (c * b);
	// long long int w = x*Pp;
	// long long int y = x + prime*(w&(R-1));
	// long long int z =  y>>32 ;
	// if(z >= prime){
	// 	z = z - prime;
	// }
	// return z;

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
      : "=&d" (b)
      : "a"(b),"rm"(c),"b"((unsigned long long int)Pp),"c"(prime)
      :"rsi","rdi");
  return b;
}
long long int SmallPrimeField::getRsquare(){
  
  long long int res;
   __asm__(
      "movq %%rax,%%rsi\n\t"
     // "push %%rax\n\t"
      //"push %%rdx\n\t"
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "divq %1\n\t"
      "movq %%rdx,%%rax\n\t"
      "mulq %%rdx\n\t"
      "divq %1\n\t"
      //"movq %%rdx,%0\n\t"
      //"pop %%rdx\n\t"
      //"pop %%rax\n\t"
      "movq %%rsi,%%rax\n\t"
      : "=&d" (res)
      : "rm"(prime)
      :"rsi","rax");
  return res;
}
//Binary inverse
SmallPrimeField SmallPrimeField::inverse2(){
	SmallPrimeField r(*this);
	long long int u = r.number();
	long long int v = prime;
	long long int x1 = 1;
	long long int x2 = 0;
	
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

SmallPrimeField SmallPrimeField::euclideanDivision(const SmallPrimeField& b, SmallPrimeField* q) const {
	if(q != NULL)
		*q = (*this)/b;
	return SmallPrimeField(0);
}

SmallPrimeField SmallPrimeField::extendedEuclidean(const SmallPrimeField& b, SmallPrimeField* s, SmallPrimeField* t) const {
	// Base Case
	// SmallPrimeField r(*this);
	// if(s == NULL || t == NULL){
	// 	return r.gcd(b);
	// }

 //    if (isZero()) { 
 //        (*s).Zero(); 
 //        (*t).One(); 
 //        return b; 
 //    } 
//TODO
	std::cerr << "SmallPrimeField::extendedEuclidean not yet implemented" << std::endl;
    // int x1, y1; // To store results of recursive call 
    // int gcd = gcdExtended(b%a, a, &x1, &y1); 
  
    // // Update x and y using results of recursive 
    // // call 
    // *x = y1 - (b/a) * x1; 
    // *y = x1; 
  
    // return gcd;
    return 1;
}

SmallPrimeField SmallPrimeField::quotient(const SmallPrimeField& b) const {
	//std::cerr << "SmallPrimeField::quotient not yet implemented" << std::endl;
	//exit(1);
	//return *this;
	SmallPrimeField r(*this);
	return r/=b;
}

SmallPrimeField SmallPrimeField::remainder(const SmallPrimeField& b) const {
	return SmallPrimeField(0);
	//std::cerr << "SmallPrimeField::remainder not yet implemented" << std::endl;
	//exit(1);
	//return *this;
}

