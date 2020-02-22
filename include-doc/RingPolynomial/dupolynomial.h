#ifndef _DUNIPOLYNOMIAL_H_
#define _DUNIPOLYNOMIAL_H_


#include "../polynomial.h"
#include <typeinfo>
#include "../DyadicRationalNumber/Multiplication/multiplication.h"	// Taylor Shift DnC
#include <math.h>
#include "../FFT/src/dft_general.h"
#include "../FFT/src/dft16.h"
#include "../FFT/src/dft_general_fermat.h"
#include <new>
#include "../Utils/TemplateHelpers.hpp"

/**
 * A univariate polynomial with arbitrary BPASField coefficients represented densely.
 * This class is templated by a Field which should be a BPASField.
 */
template <class Field>
class DenseUnivariatePolynomial : public BPASUnivariatePolynomial<Field>, 
							      public BPASEuclideanDomain 
{
	private:
		Symbol name;	// Variable name
		int curd;	// Current degree
		int n;		// Maximum size of the polynomial
		Field* coef;	// Coefficients

		inline void zeros() {
			for (int i = 0; i < n; ++i)
				coef[i] = 0;
		}
		inline bool isEqual(const DenseUnivariatePolynomial<Field>& q) const {
			if (curd && q.curd && (name != q.name))
				return 0;
			if (curd != q.curd)
				return 0;
			for (int i = 0; i <= curd; ++i) {
				if (coef[i] != q.coef[i])
				        return 0;
			}
			return 1;
		}


        inline void pomopo(Field c, Field t, const DenseUnivariatePolynomial<Field>& b) {
		    if (c == 1) {
			for (int i = curd, j = b.curd; j > -1; --i, --j) {
			    Field elem = t * b.coef[j];
			    elem += coef[i];
			    coef[i]  = elem;
			}
		    }
		    else {
			for (int i = curd, j = b.curd; i > -1; --i, --j) {
			    Field elem = coef[i] * c;
			    if (j > -1)
				elem += t * b.coef[j];
			    coef[i] = elem;
			}
		    }
		    resetDegree();
		}

    	inline void resetDegree() {
		    for (int i = curd; i > 0; --i) {
			if (coef[i] == 0)
			    curd = i - 1;
			else { break; }
		    }
		}

		// gcd subroutines
		inline DenseUnivariatePolynomial<Field> euclideanGCD (const DenseUnivariatePolynomial<Field>& q) const {
		    DenseUnivariatePolynomial<Field> a, b;
		    if (curd < q.curd) {
			a = q;
			b = *this;
		    }
		    else {
			a = *this;
			b = q;
		    }

		    while (!b.isZero()) {
			Field lc = b.coef[b.curd];
			b /= lc;
			DenseUnivariatePolynomial<Field> r;
			r.name = name;
			a.monicDivide(b, &r);
			a = b * lc;
			b = r;
		    }
		    return a / a.coef[a.curd];
		}


	public:
		mpz_class characteristic;
		// static bool isPrimeField;
        // static bool isComplexField;
        static mpz_class p1, p2, p3, p4, p5, p6, p7;  //static variables used to hold Generalized Fermat Numbers.

		/**
		 * Construct a polynomial
		 *
		 * @param d
		 **/
		DenseUnivariatePolynomial<Field> () : curd(0), n(1), name("%") {
			coef = new Field[1];
			coef[0] = 0;
			characteristic = coef[0].characteristic;
		}
		/**
		 * Construct a polynomial with degree
		 *
		 * @param d: Size of the polynomial
		 **/
		DenseUnivariatePolynomial<Field>(int s) {
			if (s < 1) { s = 1; }
			n = s;
			coef = new Field[n];
			curd = 0;
			//coef[0] = 0;
			zeros();
			name = "%";
			characteristic = coef[0].characteristic;
		}


		DenseUnivariatePolynomial<Field> (Field e) : curd(0), n(1), name("%")  {
			coef = new Field[1];
			coef[0] = e;
			characteristic = e.characteristic;
		}
		/**
		 * Copy constructor
		 *
		 * @param b: A dense univariate rational polynomial
		 **/
		DenseUnivariatePolynomial<Field>(const DenseUnivariatePolynomial<Field>& b) : curd(b.curd), name(b.name) {
			n = curd + 1;
			coef = new Field[n];
			std::copy(b.coef, b.coef+n, coef);
		}
		/**
		 * Destroy the polynomial
		 *
		 * @param
		 **/
		~DenseUnivariatePolynomial<Field>() {
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
		inline Field leadingCoefficient() const {
			return coef[curd];
		}

		inline Field trailingCoefficient() const {
			for (size_t i = 0; i <= curd; ++i) {
				if (coef[i] != 0) {
					return coef[i];
				}
			}
			return Field(0);
		}

		inline Integer numberOfTerms() const {
			size_t c = 0;
			for (size_t i = 0; i <= curd; ++i) {
				if (coef[i] != 0) {
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
		inline Field* coefficients(int k=0) {
#ifdef BPASDEBUG
			if (k < 0 || k >= n)
				std::cout << "BPAS: warning, try to access a non-exist coefficient " << k << " from DUFP(" << n << ")." << std::endl;
#endif
			return &coef[k];
		}
		/**
		 * Get a coefficient of the polynomial
		 *
		 * @param k: Offset
		 **/
		inline Field coefficient(int k) const {
			if (k < 0 || k >= n)
				return Field(0);
			return coef[k];
		}
		/**
		 * Set a coefficient of the polynomial
		 *
		 * @param k: Offset
		 * @param val: Coefficient
		 **/
		inline void setCoefficient(int k, const Field& value) {
			if (k >= n || k < 0) {
				std::cout << "BPAS: error, DUFP(" << n << ") but trying to access " << k << "." << std::endl;
				exit(1);
			}
			coef[k] = value;
			//std::cout << k << "   "  << curd << std::endl;
			if (k > curd && value != 0){
				curd = k;
			}
			resetDegree();
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
		/**
		 * Overload operator =
		 *
		 * @param b: A univariate rational polynoial
		 **/
		inline DenseUnivariatePolynomial<Field>& operator= (const DenseUnivariatePolynomial<Field>& b) {
			if (this != &b) {
				if (n) { delete [] coef; n = 0; }
				name = b.name;
				curd = b.curd;
				n = curd + 1;
				coef = new Field[n];
				std::copy(b.coef, b.coef+n, coef);
			}
			return *this;
		}

		inline DenseUnivariatePolynomial<Field>& operator= (const Field& f) {
			*this = DenseUnivariatePolynomial<Field>(f);
			return *this;
		}

		/**
		 * Overload operator !=
		 *
		 * @param b: A univariate rational polynoial
		 **/
		inline bool operator!= (const DenseUnivariatePolynomial<Field>& b) const {
			return !(isEqual(b));
		}
		/**
		 * Overload operator ==
		 *
		 * @param b: A univariate rational polynoial
		 **/
		inline bool operator== (const DenseUnivariatePolynomial<Field>& b) const {
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
				return (coef[0].isOne());
			return 0;
		}
		/**
		 * Set polynomial to 1
		 *
		 * @param
		 **/
		inline void one() {
			curd = 0;
			coef[0].one();
			for (int i = 1; i < n; ++i)
				coef[i].zero();
		}
		/**
		 * Is polynomial a constatn -1
		 *
		 * @param
		 **/
		inline bool isNegativeOne() const {
			if (!curd)
				return (coef[0].isNegativeOne());
			return 0;
		}
		/**
		 * Set polynomial to -1
		 *
		 * @param
		 **/
		inline void negativeOne() {
			curd = 0;
			coef[0].negativeOne();
			for (int i = 1; i < n; ++i)
				coef[i].zero();
		}
		/**
		 * Is a constant
		 *
		 * @param
		 **/
		inline int isConstant() const {
			if (curd) { return 0; }
			else if (!coef[0].isZero()) { return 1; }
			else { return -1; }
		}

		inline DenseUnivariatePolynomial<Field> unitCanonical(DenseUnivariatePolynomial<Field>* u, DenseUnivariatePolynomial<Field>* v) const {
			Field lead = this->leadingCoefficient();
			Field leadInv = lead.inverse();
			DenseUnivariatePolynomial<Field> ret = *this * leadInv;
			if (u != NULL) {
				*u = lead;
			}
			if (v != NULL) {
				*v = leadInv;
			}
			return ret;
		}


		/**
		 * Content of the polynomial
		 *
		 * @param
		 **/
		inline Field content() const {
			return Field((int)!isZero());
		}

		inline DenseUnivariatePolynomial<Field> primitivePart() const {
			std::cerr << "DenseUnivariatePolynomial<Field>::primitivePart() NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}


		/**
		 * Overload operator ^
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		inline DenseUnivariatePolynomial<Field> operator^ (long long int e) const {
		    DenseUnivariatePolynomial<Field> res;
		    res.name = name;
		    res.one();
		    unsigned long int q = e / 2, r = e % 2;
		    DenseUnivariatePolynomial<Field> power2 = *this * *this;
		    for (int i = 0; i < q; ++i)
			res *= power2;
		    if (r) { res *= *this; }
		    return res;
		}
		/**
		 * Overload operator ^=
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		inline DenseUnivariatePolynomial<Field>& operator^= (long long int e) {
			*this = *this ^ e;
			return *this;
		}
		/**
		 * Overload operator <<
		 * replace by muplitying x^k
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariatePolynomial<Field> operator<< (int k) const {
		    int s = curd + k + 1;
		    DenseUnivariatePolynomial<Field> r;
		    r.n = s;
		    r.curd = s - 1;
		    r.name = name;
		    r.coef = new Field[s];
		    for (int i = 0; i < k; ++i)
			r.coef[i] = 0;
		    for (int i = k; i < s; ++i)
			r.coef[i] = coef[i-k];
		    return r;
		}
		/**
		 * Overload operator <<=
		 * replace by muplitying x^k
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariatePolynomial<Field>& operator<<= (int k) {
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
		inline DenseUnivariatePolynomial<Field> operator>> (int k) const {
		    DenseUnivariatePolynomial<Field> r;
		    r.name = name;
		    int s = curd - k + 1;
		    if (s > 0) {
			r.n = s;
			r.curd = s - 1;
			delete [] r.coef;
			r.coef = new Field[s];
			for (int i = 0; i < s; ++i)
			    r.coef[i] = coef[i+k];
		    }
		}
		/**
		 * Overload operator >>=
		 * replace by dividing x^k, and
		 * return the quotient
		 *
		 * @param k: The exponent of variable, k > 0
		 **/
		inline DenseUnivariatePolynomial<Field>& operator>>= (int k) {
			*this = *this >> k;
			return *this;
		}
		
		/**
		 * Overload operator +
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariatePolynomial<Field> operator+ (const DenseUnivariatePolynomial<Field>& b) const {
		    if (!curd) { return (b + coef[0]); }
		    if (!b.curd) {
				return (*this + b.coef[0]);
			}
		    if (name != b.name) {
				std::cout << "BPAS: error, trying to add between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
				exit(1);
		    }

		    int size = (curd > b.curd)? curd+1 : b.curd+1;
		    DenseUnivariatePolynomial<Field> res;
		    res.n = size;
		    res.curd = size - 1;
		    res.name = name;
		    res.coef = new Field[size];
		    for (int i = 0; i < size; ++i) {
				Field elem = 0;
				if (i <= curd)
				    elem += coef[i];
				if (i <= b.curd)
				    elem += b.coef[i];
				res.coef[i] = elem;
		    }
		    res.resetDegree();
		    return res;
		}

		/**
		 * Overload Operator +=
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariatePolynomial<Field>& operator+= (const DenseUnivariatePolynomial<Field>& b) {
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
    	inline void add(const DenseUnivariatePolynomial<Field>& b) {
		    for (int i = curd; i >= 0; --i) {
			if (i <= b.curd)
			    coef[i] += b.coef[i];
		    }
		    resetDegree();
		}

		/**
		 * Overload Operator +
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field> operator+ (const Field& c) const {
			DenseUnivariatePolynomial<Field> r (*this);
			return (r += c);
		}

		/**
		 * Overload Operator +=
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field>& operator+= (const Field& c) {
			coef[0] += c;
			return *this;
		}

		inline friend DenseUnivariatePolynomial<Field> operator+ (Field c, DenseUnivariatePolynomial<Field> p) {
			return (p + c);
		}

		/**
		 * Subtract another polynomial
		 *
		 * @param b: A univariate rational polynomial
		 */
    	inline DenseUnivariatePolynomial<Field> operator- (const DenseUnivariatePolynomial<Field>& b) const {
		    if (!curd) { return (coef[0] - b); }
		    if (!b.curd) { return (*this - b.coef[0]); }
		    if (name != b.name) {
			std::cout << "BPAS: error, trying to subtract between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
			exit(1);
		    }


		    int size = (curd > b.curd)? curd+1 : b.curd+1;
		    DenseUnivariatePolynomial<Field> res;
		    res.n = size;
		    res.curd = size - 1;
		    res.name = name;
		    res.coef = new Field[size];
		    for (int i = 0; i < size; ++i) {
			Field elem = 0;
			if (i <= curd)
			    elem = coef[i];
			if (i <= b.curd)
			    elem -= b.coef[i];
			res.coef[i] = elem;
		    }
		    res.resetDegree();
		    return res;
		}
		/**
		 * Overload operator -=
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariatePolynomial<Field>& operator-= (const DenseUnivariatePolynomial<Field>& b) {
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
    	inline DenseUnivariatePolynomial<Field> operator- () const {
		    DenseUnivariatePolynomial<Field> res(curd+1);
		    res.name = name;
		    res.curd = curd;
		    for (int i = 0; i <= curd; ++i)
			res.coef[i] = -coef[i];
		    return res;
		}
		/**
		 * Subtract another polynomial from itself
		 *
		 * @param b: A univariate rational polynomial
		 **/
    	inline void subtract(const DenseUnivariatePolynomial<Field>& b) {
		    for (int i = curd; i >= 0; --i) {
			if (i <= b.curd)
			    coef[i] -= b.coef[i];
		    }
		    resetDegree();
		}
		/**
		 * Overload operator -
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field> operator- (const Field& c) const {
			DenseUnivariatePolynomial<Field> r (*this);
			return (r -= c);
		}

		/**
		 * Overload operator -=
		 *
		 * @param c: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field>& operator-= (const Field& c) {
			coef[0] -= c;
			return *this;
		}


		inline friend DenseUnivariatePolynomial<Field> operator- (const Field& c, const DenseUnivariatePolynomial<Field>& p) {
	        return (-p + c);
        }

        /**
		 * Helper function to compare a given prime number with manually computed Generalized Fermat Numbers
		 *
		 * @param p: prime number
		 **/
        inline bool SRGFNcmp(mpz_class &p){
        	if(cmp(p, p1)==0)
        		return true;
        	else if(cmp(p, p2)==0)
        		return true;
        	else if(cmp(p, p3)==0)
        		return true;
        	else if(cmp(p, p4)==0)
        		return true;
        	else if(cmp(p, p5)==0)
        		return true;
        	else if(cmp(p, p6)==0)
        		return true;
        	else if(cmp(p, p7)==0)
        		return true;
        	return false;
        }
	
		/**
		 * Multiply to another polynomial
		 *
		 * @param b: A univariate rational polynomial
		 **/
    	inline DenseUnivariatePolynomial<Field> operator* (const DenseUnivariatePolynomial<Field>& b) const {
		    if (!curd) {
				return (b * coef[0]);
			}
		    if (!b.curd) {
				return (*this * b.coef[0]);
			}
		    if (name != b.name) {
				std::cout << "BPAS: error, trying to multiply between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
				exit(1);
		    }

		    int d = curd + 1;
		    int m = b.curd + 1;
		    int size = curd + m;
		    DenseUnivariatePolynomial<Field> res;
			DenseUnivariatePolynomial<Field> *res2 = new DenseUnivariatePolynomial<Field>();
			res2->n = size;
			res2->curd = size-1;
			res2->name = name;
			res2->coef = new Field[size];

		    res.n = size;
		    res.curd = size - 1;
		    res.name = name;
		    res.coef = new Field[size];

		    if(Field::characteristic == 0){ 	//Q, R, C
		    	//TODO
				// if(Field::isPrimeField == 1){				//if Field = Q
				// 	std::cout << "we have rational number coefficients!" << std::endl;//to avoid compiler error
				// 	RationalNumber rn;
				// 	RationalNumber *coefRN = rn.RNpointer(coef);
				// 	RationalNumber *coefRNb = rn.RNpointer(b.coef);
				// 	RationalNumber *coefRNr = rn.RNpointer(res.coef);
				// 	Integer aden = 1, bden = 1;
				// 	for (int i = 0; i < d; ++i)
				// 		aden *= coefRN[i].get_den();
				// 	for (int i = 0; i < m; ++i)
				// 		bden *= coefRNb[i].get_den();
				// 	lfixz* acoef = new lfixz[d];
				// 	for (int i = 0; i < d; ++i)
				// 		acoef[i] = aden.get_mpz() * coefRN[i].get_mpq();
				// 	lfixz* bcoef = new lfixz[m];
				// 	for (int i = 0; i < m; ++i)
				// 		bcoef[i] = bden.get_mpz() * coefRNb[i].get_mpq();

				// 	lfixz* mul = new lfixz[size];
				// 	univariateMultiplication(mul, acoef, d, bcoef, m);  //need to change , comes from Multiplication.h class
				// 	lfixz den = aden * bden;
				// 	for (int i = 0; i < size; ++i) {
				// 		mpq_class elem = mul[i];
				// 		elem /= den;
				// 		coefRNr[i] = elem;
				// 	}

				// 	delete [] acoef;
				// 	delete [] bcoef;
				// 	delete [] mul;

				// 	res.resetDegree();
				// 	return res;
				// }
				//generic multiplcation complex Field
				// else{
					for(int i=0; i<d; i++){
						for(int j=0; j<m; j++){
							res.coef[i+j] += coef[i] * b.coef[j];
						}
					}
					res.resetDegree();
					return res;
				// }

		    }
		    else{  //prime charactersitcs
			  //   if(Field::isPrimeField == 1){
					// int de = this->degree() + b.degree() + 1;
					// //int de = 9999 + 4999 + 1;
					// //std::cout << "de: " << de << std::endl;
					// int  e = ceil(log2(de));
					// int N = exp2(e);
					// mpz_class ch(Field::characteristic);
		   //  		if((ch-1)%N == 0){
		   //  			if(/*Field::isSmallPrimeField*/0){ //is small prime
					// 		std::cout << "we have small prime field!" << std::endl;
					// 		SmallPrimeField pfield;

					// 		mpz_class userprime;
					// 		userprime = pfield.Prime();
					// 		SmallPrimeField nth_primitive = pfield.findPrimitiveRootofUnity(N);

					// 		SmallPrimeField spf;
					// 		SmallPrimeField *coefSPF = spf.SPFpointer(coef);
					// 		SmallPrimeField *coefSPFb = spf.SPFpointer(b.coef);
					// 		SmallPrimeField *coefSPFr = spf.SPFpointer(res.coef);

					// 		SmallPrimeField *u = new SmallPrimeField[N];
					// 		SmallPrimeField *v = new SmallPrimeField[N];
					// 		for(int j=0; j<N; j++){
					// 		  u[j].zero();
					// 		  v[j].zero();
					// 		}

					// 		std::copy(coef, coef+this->curd+1, u);
					// 		std::copy(b.coef, b.coef+b.curd+1, v);

					// 		this->FFT(u, 2, e, nth_primitive);
					// 		b.FFT(v, 2, e,  nth_primitive);

					// 		for(int i=0; i<N; i++){
					// 		  u[i] *= v[i];
					// 		}


					// 		this->IFFT(u, 2, e, nth_primitive);

					// 		std::copy(&u[0], &u[res.n], res.coef);
					// 		res.resetDegree();
					// 		return res;
					// 	}else if(0/*SRGFNcmp(ch)*/){
					// 		std::cout << "Generalized Fermat Prime Field Multiplication" << std::endl;
					// 		GeneralizedFermatPrimeField pfield;

					// 		mpz_class userprime;
					// 		userprime = pfield.Prime();
					// 		GeneralizedFermatPrimeField nth_primitive = pfield.findPrimitiveRootofUnity(N);

					// 		GeneralizedFermatPrimeField gpf;
					// 		GeneralizedFermatPrimeField *coefSPF = gpf.GPFpointer(coef);
					// 		GeneralizedFermatPrimeField *coefSPFb = gpf.GPFpointer(b.coef);
					// 		GeneralizedFermatPrimeField *coefSPFr = gpf.GPFpointer(res.coef);

					// 		GeneralizedFermatPrimeField *u = new GeneralizedFermatPrimeField[N];
					// 		GeneralizedFermatPrimeField *v = new GeneralizedFermatPrimeField[N];
					// 		for(int j=0; j<N; j++){
					// 		  u[j].zero();
					// 		  v[j].zero();
					// 		}

					// 		std::copy(coefSPF, coefSPF+this->curd+1, u);
					// 		std::copy(coefSPFb, coefSPFb+b.curd+1, v);
					// 		this->FFT(u, 16, e, nth_primitive);
					// 		b.FFT(v, 16, e,  nth_primitive);

					// 		for(int i=0; i<N; i++){
					// 		  u[i] *= v[i];
					// 		}

					// 		this->IFFT(u, 16, e, nth_primitive);
					// 		for(int j=0; j<res.n; j++){
					// 		coefSPFr[j] = u[j];
					// 		//coefSPFr[j] = u[j].number()%userprime;
					// 		}
					// 		res.resetDegree();
					// 		return res;
			  //             }else{
					// 		//std::cout << "Big Prime Field Multiplication" << std::endl;
					// 		//BigPrimeField pfield; // big prime field for using some of its function
					// 		BigPrimeField nth_primitive = BigPrimeField::findPrimitiveRootofUnity(mpz_class(N)); //cast N to MPZ CLASS

					// 		BigPrimeField *u = new BigPrimeField[N];
					// 		BigPrimeField *v = new BigPrimeField[N];

					// 		std::copy(coef, (coef)+this->curd+1, u);					
					// 		std::copy(b.coef, (b.coef)+b.curd+1, v);

					// 		this->FFT(u, 2, e, nth_primitive);
					// 		b.FFT(v, 2, e,  nth_primitive);

					// 		for(int i=0; i<N; i++){
					// 		  u[i] *= v[i];
					// 		}

					// 		this->IFFT(u, 2, e, nth_primitive);
					// 		std::copy(u, u+res2->n, res2->coef);

					// 		delete[] u;
					// 		delete[] v;
					// 		res2->resetDegree();
					// 		return *res2;
			  //             } //end of inner if ch <= 962592769
		   //  		}else{
		   //  			for(int i=0; i<curd+1; i++){
		   //                        for(int j=0; j<b.curd+1; j++){
		   //                           	 res.coef[i+j] += coef[i] * b.coef[j];
				 //              }
				 //          }
		   //  		}//end of outer if (ch-1)%N == 0
			  //   }else{
			    	for(int i=0; i<d; i++){
		                          for(int j=0; j<m; j++){
		                             	 res.coef[i+j] += coef[i] * b.coef[j];
			              }
			          }
			    // }//end of Field::isPrimeField == 1
			}//end of field::char if

		}//end of function		

		inline void FFT(SmallPrimeField* field, int K, int e, SmallPrimeField& omega){
			//std::cout  << "FFT omega "<<  omega << std::endl;
			DFT_general(field, K, e, omega);
			//DFT_16(field, omega);
		}

		inline void FFT(BigPrimeField* field, int K, int e, BigPrimeField& omega){
			DFT_general(field, K, e, omega);
		}

		inline void FFT(GeneralizedFermatPrimeField* field, int K, int e, GeneralizedFermatPrimeField& omega){
			std::cout << " i called generalized FFT " << std::endl;
			dft_general_fermat(field, K, e, omega);
		}

		inline void IFFT(SmallPrimeField* field, int K, int e, SmallPrimeField& omega){
			/*
			int N = exp2(e);
			SmallPrimeField N_inv(N);
		    	N_inv = N_inv.inverse();
		    	//std::cout  << "N inverse   " <<  N_inv << std::endl;
		    	SmallPrimeField omg = omega.inverse();
		    	//std::cout  << "IFFT omega inverse   " <<  omg << std::endl;
		    	DFT_general(field,  K,  e, omg);
		    	//DFT_16(field, omg);
		    	for(int i=0; i<N; i++){
		    		field[i] *= N_inv;
		    	}
		    	*/
		    	inverse_DFT(field, K, e, omega);
		}

		inline void IFFT(BigPrimeField* field, int K, int e, BigPrimeField& omega){
			/*
			int N = exp2(e);
			BigPrimeField N_inv(N);
			N_inv = N_inv.inverse();
		    	BigPrimeField omg = omega.inverse();

		    	DFT_general(field,  K,  e, omg);
		    	for(int i=0; i<N; i++){
		    		field[i] *= N_inv;
		    	}
		    	*/
		    	inverse_DFT(field, K, e, omega);
		}

		inline void IFFT(GeneralizedFermatPrimeField* field, int K, int e, GeneralizedFermatPrimeField& omega){
			/*
			int N = exp2(e);
			GeneralizedFermatPrimeField N_inv(N);
		    	N_inv = N_inv.inverse();
		    	GeneralizedFermatPrimeField omg = omega.inverse();

		    	DFT_general(field,  K,  e, omg);
		    	for(int i=0; i<N; i++){
		    		field[i] *= N_inv;
		    	}
		    	*/
		    	inverse_fermat_DFT(field, K, e, omega);
		}


		/**
		 * Overload operator *=
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariatePolynomial<Field>& operator*= (const DenseUnivariatePolynomial<Field>& b) {
			*this = *this * b;
			return *this;
		}
		/**
		 * Overload operator *
		 *
		 * @param e: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field> operator* (const Field& e) const {
			DenseUnivariatePolynomial<Field> r (*this);
			return (r *= e);
		}

		/**
		 * Overload operator *=
		 *
		 * @param e: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field>& operator*= (const Field& e) {
		    if (e != 0 && e != 1) {
			for (int i = 0; i <= curd; ++i)
			    coef[i] *= e;
		    }
		    else if (e == 0) { zero(); }
		    return *this;
		}

		inline friend DenseUnivariatePolynomial<Field> operator* (const Field& c, const DenseUnivariatePolynomial<Field>& p) {
            return (p * c);
        }

		/**
		 * Overload operator /
		 * ExactDivision
		 *
		 * @param b: A univariate rational polynomial
		 **/
		inline DenseUnivariatePolynomial<Field> operator/ (const DenseUnivariatePolynomial<Field>& b) const {
			DenseUnivariatePolynomial<Field> rem(*this);
			return (rem /= b);
		}
		/**
		 * Overload operator /=
		 * ExactDivision
		 *
		 * @param b: A univariate rational polynomial
		 **/
    	inline DenseUnivariatePolynomial<Field>& operator/= (const DenseUnivariatePolynomial<Field>& b) {
		    if (b.isZero()) {
			std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
			exit(1);
		    }
		    if (!b.curd)
			return (*this /= b.coef[0]);
		    if (!curd) {
			coef[0].zero();
			return *this;
		    }

		    if (name != b.name) {
			std::cout << "BPAS: error,rem trying to exact divide between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
			exit(1);
		    }

		    DenseUnivariatePolynomial<Field> rem(*this);
		    zeros();
		    while (rem.curd >= b.curd) {
			Field lc = rem.coef[rem.curd] / b.coef[b.curd];
			int diff = rem.curd - b.curd;
			rem.pomopo(1, -lc, b);
			coef[diff] = lc;
		    }
		    resetDegree();
		    if (!rem.isZero()) {
			std::cout << "BPAS: error, not exact division from DUFP." << std::endl;
			exit(1);
		    }
		    return *this;
		}

		/**
		 * Overload operator /
		 *
		 * @param e: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field> operator/ (const Field& e) const {
			DenseUnivariatePolynomial<Field> r (*this);
			return (r /= e);
		}

		/**
		 * Overload operator /=
		 *
		 * @param e: A rational number
		 **/
		inline DenseUnivariatePolynomial<Field>& operator/= (const Field& e) {
		    if (e == 0) {
			std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
			exit(1);
		    }
		    else if (e != 1) {
		    	for (int i = 0; i <= curd; ++i)
		    	    coef[i] /= e;
		    }
		    return *this;
		}

    	//friend DenseUnivariatePolynomial<Field> operator/ (Field c, DenseUnivariatePolynomial<Field> p);~
    	inline friend DenseUnivariatePolynomial<Field> operator/ (const Field& c, const DenseUnivariatePolynomial<Field>& p) {
		    if (p.isZero()) {
			std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
			exit(1);
		    }

		    DenseUnivariatePolynomial<Field> q;
		    q.name = p.name;
		    q.curd = 0;
		    q.n = 1;
		    q.coef = new Field[1];

		    if (p.isConstant())
			q.coef[0] = c / p.coef[0];
		    else
			q.coef[0].zero();
		    return q;
		}

		/**
		 * Monic division
		 * Return quotient and itself become the remainder
		 *
		 * @param b: The dividend polynomial
		 **/
    	inline DenseUnivariatePolynomial<Field> monicDivide(const DenseUnivariatePolynomial<Field>& b) {
		    if (b.isZero()) {
			std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
			exit(1);
		    }
		    else if (b.coef[b.curd] != 1) {
			std::cout << "BPAS: error, leading coefficient is not one in monicDivide() from DUFP." << std::endl;
			exit(1);
		    }
		    if (!b.curd) {
			DenseUnivariatePolynomial<Field> r (*this);
			zero();
			return r;
		    }
		    if (!curd) {
			DenseUnivariatePolynomial<Field> r;
			r.name = name;
			return r;
		    }
		    if (name != b.name) {
			std::cout << "BPAS: error, trying to monic divide between [" << name << "] and Q[" << b.name << "]." << std::endl;
			exit(1);
		    }

		    int size = curd - b.curd + 1;
		    DenseUnivariatePolynomial<Field> quo(size);
		    quo.curd = size - 1;
		    quo.name = name;
		    while (curd >= b.curd) {
			Field lc = coef[curd];
			int diff = curd - b.curd;
			pomopo(1, -lc, b);
			quo.coef[diff] = lc;
		    }
		    quo.resetDegree();
		    return quo;
		}
		/**
		 * Monic division
		 * Return quotient
		 *
		 * @param b: The dividend polynomial
		 * @param rem: The remainder polynomial
		 **/
    	inline DenseUnivariatePolynomial<Field> monicDivide(const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* rem) const {
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
    	inline DenseUnivariatePolynomial<Field> lazyPseudoDivide (const DenseUnivariatePolynomial<Field>& b, Field* c, Field* d) {
		    if (d == NULL)
			d = new Field;
		    int da = curd, db = b.curd;
		    if (b.isZero() || !db) {
			std::cout << "BPAS: error, dividend is zero or constant." << std::endl;
			exit(1);
		    }
		    c->one(), d->one();
		    if (!curd) {
			DenseUnivariatePolynomial<Field> r;
			r.name = name;
			return r;
		    }
		    if (name != b.name) {
			std::cout << "BPAS: error, trying to pseudo divide between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
			exit(1);
		    }

		    if (da < db) {
			DenseUnivariatePolynomial<Field> r;
			r.name = name;
			return r;
		    }

		    int size = curd - b.curd + 1;
		    DenseUnivariatePolynomial<Field> quo(size);
		    quo.curd = size - 1;
		    quo.name = name;
		    int e = 0, diff = da - db;
		    Field blc = b.coef[b.curd];
		    while (curd >= b.curd) {
			Field lc = coef[curd];
			int k = curd - b.curd;
			*c *= blc;
			e++;
			pomopo(blc, -coef[curd], b);
			quo.coef[k] = lc;
		    }
		    quo.resetDegree();
		    for (int i = e; i <= diff; ++i)
			*d *= blc;
		    return quo;
		}

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
    	inline DenseUnivariatePolynomial<Field> lazyPseudoDivide (const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* rem, Field* c, Field* d) const {
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
    	inline DenseUnivariatePolynomial<Field> pseudoDivide (const DenseUnivariatePolynomial<Field>& b, Field* d=NULL) {
		    Field c;
		    if (d == NULL)
			d = new Field;
		    DenseUnivariatePolynomial<Field> quo = lazyPseudoDivide(b, &c, d);
		    quo *= *d;
		    *this *= *d;
		    *d *= c;
		    return quo;
		}

		/**
		 * Pseudo dividsion
		 * Return the quotient
		 *
		 * @param b: The divident polynomial
		 * @param rem: The remainder polynomial
		 * @param d: The leading coefficient of b
		 *           to the power deg(a) - deg(b) + 1
		 **/
    	inline DenseUnivariatePolynomial<Field> pseudoDivide (const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* rem, Field* d) const {
		    Field c;
		    DenseUnivariatePolynomial<Field> quo = lazyPseudoDivide(b, rem, &c, d);
		    quo *= *d;
		    *rem *= *d;
		    *d *= c;
		    return quo;
		}

		/**
		 * s * a \equiv g (mod b), where g = gcd(a, b)
		 *
		 * @param b: A univariate polynomial
		 * @param g: The GCD of a and b
		 **/
    	inline DenseUnivariatePolynomial<Field> halfExtendedEuclidean (const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* g) {
		    if (g == NULL)
			g = new DenseUnivariatePolynomial<Field>;
		    *g = *this;

		    DenseUnivariatePolynomial<Field> a1, b1;
		    a1.name = name;
		    b1.name = b.name;
		    a1.coef[0].one();
		    while (!b.isZero()) {
			DenseUnivariatePolynomial<Field> q, r;
			q.name = r.name = name;
			Field e = b.coef[b.curd];
			if (e != 1) {
			    b /= e;
			    *g /= e;
			}
			q = g->monicDivide(b, &r);
			if (e != 1) {
			    *g = b * e;
			    b = r * e;
			}
			else {
			    *g = b;
			    b = r;
			}

			r = a1;
			r -= q * b1;
			a1 = b1;
			b1 = r;
		    }

		    a1 /= g->coef[g->curd];
		    *g /= g->coef[g->curd];

		    return a1;
		}

		/**
		 * s*a + t*b = c, where c in the ideal (a,b)
		 *
		 * @param a: A univariate polynomial
		 * @oaran b: A univariate polynomial
		 * @param s: Either s = 0 or degree(s) < degree(b)
		 * @param t
		 **/
    	inline void diophantinEquationSolve(const DenseUnivariatePolynomial<Field>& a, const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* s, DenseUnivariatePolynomial<Field>* t) {
		    DenseUnivariatePolynomial<Field> f(*this), g, q, r;
		    *s = a.halfExtendedEuclidean(b, &g);
		    if (g.coef[g.curd] != 1) {
			f /= g.coef[g.curd];
			g /= g.coef[g.curd];
		    }
		    q = f.monicDivide(g, &r);
		    if (!r.isZero()) {
			std::cout << "BPAS: error, " << *this << " is not in the ideal (" << a << ", " << b << ") from DUQP." << std::endl;
			exit(1);
		    }
		    *s *= q;

		    Field e = b.coef[b.curd];
		    if (e != 1) { b /= e; }
		    if (s->curd >= b.curd) {
			*s /= e;
			s->monicDivide(b, &r);
			*s = r * e;
		    }

		    g = *this;
		    g -= *s * a;
		    if (e != 1) { g /= e; }
		    *t = g.monicDivide(b);
		}

		/**
		 * Compute k-th differentiate
		 *
		 * @param k: k-th differentiate, k > 0
		 **/
    	inline void differentiate(int k) {
    		*this = this->derivative(k);
		}

		inline void differentiate() {
			return this->differentiate(1);
		}

		inline DenseUnivariatePolynomial<Field> derivative(int k) const {
			DenseUnivariatePolynomial<Field> ret(*this);
		    if (k <= 0) { return *this; }
		    for (int i = k; i <= ret.curd; ++i) {
			ret.coef[i-k] = ret.coef[i];
			for (int j = 0; j < k; ++j)
			    ret.coef[i-k] *= (i - j);
		    }
		    ret.curd -= k;
		    ret.resetDegree();
		    return ret;
		}

		inline DenseUnivariatePolynomial<Field> derivative() const {
			return derivative(1);
		}

		/**
		 * Compute the integral with constant of integration 0
		 *
		 * @param
		 **/
    	inline DenseUnivariatePolynomial<Field> integrate() const {
		    DenseUnivariatePolynomial<Field> b;
		    b.name = name;
		    b.n = curd+2;
		    b.coef = new Field[b.n];
		    b.coef[0].zero();
		    for (int i = 0; i <= curd; ++i)
			b.coef[i+1] = coef[i] / (i + 1);
		    b.resetDegree();
		    return b;
		}
		/**
		 * Evaluate f(x)
		 *
		 * @param x: Evaluation point
		 **/
		inline Field evaluate(const Field& x) const {
			if (curd) {
				Field px = coef[curd];
				for (int i = curd-1; i > -1; --i)
					px = px * x + coef[i];
				return px;
			}
			return coef[0];
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
		 * GCD(p, q)
		 *
		 * @param q: The other polynomial
		 **/
    	inline DenseUnivariatePolynomial<Field> gcd (const DenseUnivariatePolynomial<Field>& q, int type) const {
		    if (isZero()) { return q; }
		    if (q.isZero()) { return *this; }
		    if (!curd || !q.curd) {
			DenseUnivariatePolynomial<Field> h (1);
			h.coef[0].one();
			h.name = name;
			return h;
		    }

		    if (name != q.name) {
			std::cout << "BPAS: error, remtrying to compute GCD between Q[" << name << "] and Q[" << q.name << "]." << std::endl;
			exit(1);
		    }

			DenseUnivariatePolynomial<Field> r;
			// if (!type)
				r = euclideanGCD(q);
			//TODO modularGCD is not defined anywhere???
			// else
				// r = modularGCD(q);


			return r;
		}

		inline DenseUnivariatePolynomial<Field> gcd(const DenseUnivariatePolynomial<Field>& q) const {
			return gcd(q, 0);
		}

		/**
		 * Square free
		 *
		 * @param
		 **/
    	inline Factors<DenseUnivariatePolynomial<Field> > squareFree() const {
		    std::vector<DenseUnivariatePolynomial<Field> > sf;
		    if (!curd) {
				sf.push_back(*this);
		    }
		    else if (curd == 1) {
				DenseUnivariatePolynomial<Field> t;
				t.name = name;
				t.coef[0] = coef[curd];
				sf.push_back(t);
				t = *this / t.coef[0];
				sf.push_back(t);
		    }
		    else {
				DenseUnivariatePolynomial<Field> a (*this), b(*this);
				b.differentiate(1);
				DenseUnivariatePolynomial<Field> g = a.gcd(b);
				DenseUnivariatePolynomial<Field> x = a / g;
				DenseUnivariatePolynomial<Field> y = b / g;
				DenseUnivariatePolynomial<Field> z = -x;
				z.differentiate(1);
				z += y;

				while (!z.isZero()) {
				    g = x.gcd(z);
				    sf.push_back(g);
				    x /= g;
				    y = z / g;
				    z = -x;
				    z.differentiate(1);
				    z += y;
				}
				sf.push_back(x);

				Field e;
				e.one();
				for (int i = 0; i < sf.size(); ++i) {
				    e *= sf[i].coef[sf[i].curd];
				    sf[i] /= sf[i].coef[sf[i].curd];
				}
				DenseUnivariatePolynomial<Field> t;
				t.name = name;
				t.coef[0] = e;
				sf.insert(sf.begin(), t);
		    }

		    Factors<DenseUnivariatePolynomial<Field>> f;
		    f.setRingElement(sf[0]);
		    for (int i = 1; i < sf.size(); ++i) {
		    	f.addFactor(sf[i], i);
		    }
		    return f;
		}
		/**
		 * Divide by variable if it is zero
		 *
		 * @param
		 **/
		inline bool divideByVariableIfCan() {
			if (coef[0] != 0)
				return 0;
			else {
				curd--;
				for (int i = 0; i <= curd; ++i)
				        coef[i] = coef[i+1];
				return 1;
			}
		}
		/**
		 * Number of coefficient sign variation
		 *
		 * @param
		 **/
		//int numberOfSignChanges();

		/**
		 * Revert coefficients
		 *
		 * @param
		 **/
		inline void reciprocal() {
			for (int i = 0; i < (curd+1)/2; ++i) {
				Field elem = coef[i];
				coef[i] = coef[curd-i];
				coef[curd-i] = elem;
			}
			resetDegree();
		}
		/**
		 * Homothetic operation
		 *
		 * @param k > 0: 2^(k*d) * f(2^(-k)*x);
		 **/
		inline void homothetic(int k) {
			for (int i = 0; i <= curd; ++i)
				coef[i] <<= (curd - i) * k;
		}
		/**
		 * Scale transform operation
		 *
		 * @param k > 0: f(2^k*x)
		 **/
		inline void scaleTransform(int k) {
			for (int i = 0; i <= curd; ++i)
				coef[i] <<= k * i;
		}
		/**
		 * Compute f(-x)
		 *
		 * @param
		 **/
		inline void negativeVariable() {
			for (int i = 0; i <= curd; ++i) {
				if (i%2)
				        coef[i] = -coef[i];
			}
		}

		/**
		 * Compute -f(x)
	         *
	         * @param
		 **/
		inline void negate() {
			for (int i = 0; i <= curd; ++i)
				coef[i] = -coef[i];
		}


		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param b: A univariate rational polynoial
		 **/
		inline void print(std::ostream& out) const {		
			//std::cout << b.curd << std::endl;
			bool isFirst = 1;
			for (int i = 0; i <= this->curd; ++i) {
				if (!this->coef[i].isZero()) {
					if (!isFirst)
						out << "+";
					if (i) {
						if (!this->coef[i].isOne())
							out << this->coef[i] << "*";
						out << this->name;
						if (i > 1)
							out << "^" << i;
					}
					else { out << this->coef[i]; }
					isFirst = 0;
				}
			}
			if (isFirst) { out << "0"; }
		}

		inline ExpressionTree convertToExpressionTree() const {
			std::cerr << "DenseUnivariatePolynomial<Field>::convertToExpressionTree() NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return ExpressionTree();
		}



		/** Euclidean domain methods **/
		
		DenseUnivariatePolynomial<Field> euclideanSize() const {
			return DenseUnivariatePolynomial<Field>(degree().get_si());
		}

		DenseUnivariatePolynomial<Field> euclideanDivision(const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* q = NULL) const {
			std::cerr << "DenseUnivariatePolynomial::ExactDivision NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		DenseUnivariatePolynomial<Field> extendedEuclidean(const DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* s = NULL, DenseUnivariatePolynomial<Field>* t = NULL) const {
			std::cerr << "DenseUnivariatePolynomial::extendedEuclidean NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		DenseUnivariatePolynomial<Field> quotient(const DenseUnivariatePolynomial<Field>& b) const {
			std::cerr << "DenseUnivariatePolynomial::quotient NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;	
		}

		DenseUnivariatePolynomial<Field> remainder(const DenseUnivariatePolynomial<Field>& b) const {
			std::cerr << "DenseUnivariatePolynomial::remainder NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		DenseUnivariatePolynomial<Field> operator% (const DenseUnivariatePolynomial<Field>& b) const {
			DenseUnivariatePolynomial<Field> ret(*this);
			ret %= b;
			return ret;
		}

		DenseUnivariatePolynomial<Field>& operator%= (const DenseUnivariatePolynomial<Field>& b) {
			*this = remainder(b);
			return *this;
		}

		/**
		 * Reverse a polynomial
		 *(based on Modern Computer Algebra , Newton Iteration Chapter 9)
		 *
		 * @param
		 **/
		inline DenseUnivariatePolynomial<Field> Reverse() const {
			DenseUnivariatePolynomial<Field> y(this->degree()+1);
			y.setVariableName(this->name);
			for(int i = this->degree(); i >=0; i--){
				 y.setCoefficient(this->degree()-i, this->coefficient(i));
			}
			return y;
		}

		/**
		 * Does Newton Iteration to compute the Reverse Inverse of a polynomial
		 *(based on Modern Computer Algebra , Newton Iteration Chapter 9)
		 *
		 * @param l:  number of newton iteration (n - m + 1)
		 **/
		inline DenseUnivariatePolynomial<Field> NewtonIterationInversion(int l){
			DenseUnivariatePolynomial<Field> g(1);
			g.one(); //g0 = 1;
			DenseUnivariatePolynomial<Field> two(1);
			two = g + g; //constant 2
			//implement 2*g0 - f*g0^2
			int r = ceil(log2(l));
			for(int i=0; i<r; i++){
				g = (g*two) - ((*this)*(g*g));
				int deg = g.degree();
				int it = pow(2, i+1);
				DenseUnivariatePolynomial<Field> temp(it);
				temp.setVariableName(this->name);
				temp.zero();
				for(int j=0; j<it; j++){
					temp.setCoefficient(j, g.coefficient(j));
				}
				g = temp;
				temp.zero();
			}
			*this = g;
			return *this;
		}

		/**
		 * Fast Division (based on Modern Computer Algebra , Newton Iteration Chapter 9)
		 *
		 * @param b: divisor
		 * @param q: (out) qoutient
		 * @param r: (out) remainder
		 * @param fieldVal: prime field
		 * @param l:  number of newton iteration
		 **/
		inline void NewtonDivisionQuotient(DenseUnivariatePolynomial<Field> &b,  DenseUnivariatePolynomial<Field> &q, DenseUnivariatePolynomial<Field> &r, int l){
			//std::cout << "b = " << b << std::endl;
			DenseUnivariatePolynomial<Field> bb(b.degree());
			this->setVariableName(this->name);
			int m = this->degree() - b.degree();
			bb = b.Reverse();
			DenseUnivariatePolynomial<Field> tempQuot(m+1); //quotient temp var before mod
			this->setVariableName(this->name);
			tempQuot = this->Reverse() * bb.NewtonIterationInversion(l);
			DenseUnivariatePolynomial<Field> tempModQuot(m+1); //quotient temp var after mod
			tempModQuot.setVariableName(this->name);
			for(int j=0; j<m+1; j++){
				tempModQuot.setCoefficient(j, tempQuot.coefficient(j));
			}
			q = tempModQuot.Reverse();	        //set the qoutient value
			r = (*this) - (b*q);
		}
};

template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p1("85236826359346144956638323529482240001", 10);
template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p2("115763822272329310636028559609001827025179711501300126126825041166177555972097", 10);
template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p3("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p4("41855814947416230160905824077102044107723669195743478941314613188961602997492974631771080284897031989884853955219823218284507222180203198418365663867454557929637857661786690253443410218509127538830031125593957763469901695303351030684228602570766096943274989297544663448958866408079360000000000000001",10);
template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p5("2877263539710446421731156102715413817501924683780203686146470704479400367915116616138891969740546502029403969510010797456265193559830952776697283528940113468571280365197702193026548850188074400758303550622614323711553963509499617917475399195332219370635307438240127286098723642161848281710791141158514433179269730389996820841870847439976376347222063201890220935378395463717774388387120384211168483542746311978533439751087743189685602954995966307828869222780207146437854468321622285523442498620491234986821135672225672571988303097617848628157952640424574080207395225600000000000000000000000000000001",10);
template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p6("56616002759896959989415322854958727265928518265558894081721915112338478642152987382081777901106335963367885845155763417916565153797308803666707803977539356324269844089304564587076354736376986405278272779897162317130287620982152662704148023386047880191229049189096093126016047741788639253285397227699187048459276900664828446151952389379779872483448827148216422038130194295758324036868739847281348652718382706086985165783151098320403070769834562613385987477083199297963686734355937454052502472824036058010003416215471582444408483889852282411679533336856925440851946178882063994444491916045719853726844950223273102469496195879587497097269220475494120869128671848644970662426110042888649905237795797441676015340881806045138767954343337713455178248747715202796137422620937669272350685622178511493534367574328398416598592863781482592003147323049228191952409766232800103525201364901470397364906237602845461804081111668634216071070833972247292545003189326173445338852542217465888156430817411944210832836729444740015333742919095761005296073941774939988656174783310275846314511363344139168277105168250508129957477572900470266117913389691465601936338974995950792881631257549184987642764414960873189711746068672863213585432577",10);
template <class Field>
mpz_class DenseUnivariatePolynomial<Field>::p7("1090748133587739207496133104623209685652043340944525583769922043686052876979141040038395373676653597591795495926807366279537470195844376816521507523618991686398153905679787781478759328905566212525054647896231628355359135635798293361069723097331736557050988577986243892532466698313135887886138444662377632073755381143950221100025346676696096637756641164069467382471667320021238955178621076923370794798694345041333296701870548156651489758101607968522281523969635189200754728978386730361168951561126882965293650334538609649761122401061010853886509602124698406378631272561442673369233937827268652039439967721278901678024095092864613759484683254583736063704358175932953266749684554868685618312267884858266625457439986285293580091133149878361352074448701717504456518088365128951751801324867934946295450111407936849396611409550525348661741918663903030770894532253322458469943425510889282843922868145807907363320595378693696835191869459437314830716877183603950247193579548207267489507689973156425200759112179070584550934920498931256877825412489146370885286547409945479643242694187863206290671548410690203340958006276613646712414730669667136223216464966782734564195802691203036768257805594494102413460722901406746958482564038382712381753485759119365977102301962480617202907182258596705133279208778898957530866346338095658945650199559035759233831660684522622662629611826038362788351572091930801750440291059713981728995722747793865224424707283629868761313993150399745881815938602359160601744218529135375895347494668407421931150844787149547879919791239585278863789643049041266244553112251574272797287440972746560481170921474633126611400596027422225750981318016619477921226148845565965761129968247215480526365772473453996471123113518375556705346341505539769221542741170266489368211724332244794506899906126995079451869450859316306542657546308047179341725779510114986968935109860258498110066067846826136020542550539706591317400144521343503647697440441721640835065587388170041690014385725998766494102073306124651052994491828044263831344641288860621253508151934136220329913691383260972274841983306863608939116025856836931995877644594086124860302354708294645851839410316355926145110650494378186764629875649766904746646788333914098179144689393484186874701417787934689792155194373081952167714271117559720776224173158712008495626957577178627392043052139387289600000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001",10);

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


