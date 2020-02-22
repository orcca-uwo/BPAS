#pragma once
/**
 * Polynomial class for representing Univariate Integer Polynomial in the densed-fashion which is an array of integers using GMP Integer Objects.
 *                         
 * This class provide basic functions for polynomials ranging from arithmetic, to utils.
 *
 * @author Farnam Mansouri
 * 
 */


#include "Multiplication/Multiplication.h"
#include "Poly.h"

class DensedUnivariateIntegerPolynomial : public BPASUnivariatePolynomial{
	private:

		UnivariateIntegerPolynomial p;
		Multiplication m;

		/**
		 * Polynomial Constructor.
		 * This will instantiate the object by passing the UnivariateIntegerPolynomial object which 
		 * is an internal class and not exposed to the user.
		 * This is a private constructor which is being used in some functions.
		 *
		 * @param poly internal polynomial instance
		 */
		DensedUnivariateIntegerPolynomial(UnivariateIntegerPolynomial poly){
                        p = poly;
                }

		/**
		 * Get the pointer to the UnivariateIntegerPolynomial object.
		 */
                inline UnivariateIntegerPolynomial* getP(){
                        return &p;
                }

		inline void pomopo(mpz_class c, mpz_class t, DensedUnivariateIntegerPolynomial& b) {
                        if (c == 1) {
                                for (int i = degree(), j = b.degree(); j > -1; --i, --j) {
                                        mpz_class elem = coefficient(i) + t * b.coefficient(j);
                                        setCoefficient(i, elem);
                                }
                        }
                        else {
                                for (int i = degree(), j = b.degree(); i > -1; --i, --j) {
                                        mpz_class elem = coefficient(i) * c;
                                        if (j > -1)
                                                elem += t * b.coefficient(j);
                                        setCoefficient(i, elem);
                                }
                        }
                }
		
	public:

		/******************************* CONSTRUCTORS **************************************/

		DensedUnivariateIntegerPolynomial () : p() {}
		/**
		 * Polynomial Constructor.
		 * This will instantiate new array with (with the input size) for the coefficients...
		 * This generate the zero polynomial: P(x) = 0.
		 * The real size of the polynomial will be zero! since the coefficients are all zero.
		 *
		 * @param s size of the polynomial
		 */
		DensedUnivariateIntegerPolynomial(int s){
			p = UnivariateIntegerPolynomial(s);
		}

		/**
		 * Polynomial Constructor.
		 * Make sure that you don't delete the pointer to the coefficients when using this constructor.
		 *
		 * @param s size of the polynomial
		 * @param c the pointer to the coefficients
		 */
		DensedUnivariateIntegerPolynomial(int s, mpz_class * c){
                        p = UnivariateIntegerPolynomial(s, c);
                }

		/**
		 * Polynomial Constructor.
		 * This constructor shall be used when knowing the maximum number of the bits for the coefficients.
		 * Make sure that you don't delete the pointer to the coefficients when using this constructor.
		 *
		 * @param s size of the polynomial
		 * @param cDigits Maximum number of bits for the coefficients
		 * @param c the pointer to the coefficients
		 */
		DensedUnivariateIntegerPolynomial(int s, int cDigits, mpz_class * c){
                        p = UnivariateIntegerPolynomial(s, cDigits, c);
                }

		/**
		 * The Copy constructor
		 *
		 * @param q the polynomial which it copies from
		 */
		DensedUnivariateIntegerPolynomial(const DensedUnivariateIntegerPolynomial & q){
			p = UnivariateIntegerPolynomial(q.p);
		}

		
		/**
		 * Destroy the Polynomial object and free all of the allocated memory.
		 */
		~DensedUnivariateIntegerPolynomial() {
                        p.freeHeap();
                }


		/********************************* OVERRIDEN ***************************************/

		/**
		 * Check weather the polynomial is zero or not.
		 */
		inline bool isZero() {
			return (p.getExposedSize() == 1 && p.getCoefficient(0) == 0);
		}

		/**
		 * Set the polynomial to zero.
		 */
		inline void zero() {
			p.setToZero();
		}

		/**
		 * Check weather the polynomial is one or not.
		 */
		inline bool isOne() {
			return (p.getExposedSize() == 1 && p.getCoefficient(0) == 1);
		}

		/**
		 * Set the polynomial to the constant polynomial with the value of 1.
		 */
		inline void one() {
			p.setToOne();
		}

		/**
		 * Overloading the = operator (assignment).
		 *
		 * @param q the given polynomial which we assign from
		 */
		inline DensedUnivariateIntegerPolynomial & operator= (DensedUnivariateIntegerPolynomial q){
			if(this != &q)
				p.assign(q.getP());
			return *this;
		}

		/**
		 * Overloading the += operator. 
		 *
		 * @param q the polynomial to be added to
		 */
		inline DensedUnivariateIntegerPolynomial & operator+= (DensedUnivariateIntegerPolynomial q){
			if(size() >= q.size())	
				add(q);		
			else	
				*this = *this + q;
			return *this;
		}

		/**
		 * Overloading the + operator. Add the polynomial to another polynomial and return the result polynomial.
		 *
		 * @param q the polynomial to be added to
		 */
		inline DensedUnivariateIntegerPolynomial operator+ (DensedUnivariateIntegerPolynomial & q){
			UnivariateIntegerPolynomial r = p.add(q.getP());
			DensedUnivariateIntegerPolynomial result(r);
                        return result;
		}

		/**
		 * Add another polynomial to the polynomial in-place.
		 *
		 * @param q the polynomial which will be added from
		 */
		inline void add(DensedUnivariateIntegerPolynomial & q){
			p.addInPlace(q.getP());
                }

		/**
		 * Overloading the -= operator.
		 *
		 * @param q the polynomial which will be substracted.
		 */
		inline DensedUnivariateIntegerPolynomial & operator-= (DensedUnivariateIntegerPolynomial q){
                        if(size() >= q.size())    
                                substract(q);         
                        else    
                                *this = *this - q;
                        return *this;
                }


		/**
		 * Overloading the - operator. Subtract the polynomial from another polynomial and return the result polynomial.
		 *
		 * @param q the polynomial which will be substracted.
		 */
		inline DensedUnivariateIntegerPolynomial operator- (DensedUnivariateIntegerPolynomial & q){
                        UnivariateIntegerPolynomial r = p.substract(q.getP());
                        DensedUnivariateIntegerPolynomial result(r);
                        return result;
                }

		/**
		 * Substract the polynomial from another polynomial in-place.
		 *
		 * @param q the polynomial which will be substrated
		 */
		inline void substract(DensedUnivariateIntegerPolynomial & q){
                        p.substractInPlace(q.getP());
                }

		/**
		 * Overloading the *= operator.
		 *
		 * @param q the polynomial to be multiplied with
		 */
		inline DensedUnivariateIntegerPolynomial & operator*= (DensedUnivariateIntegerPolynomial q){
			*this = *this * q;
                        return *this;
		}

		/**
		 * Overloading the * operator. Multiply the polynomial to another polynomial and return the result polynomial.
		 *
		 * @param q the polynomial to be multiplied with
		 */
		inline DensedUnivariateIntegerPolynomial operator* (DensedUnivariateIntegerPolynomial & q){
			UnivariateIntegerPolynomial r = m.multiply(&p, q.getP());
			DensedUnivariateIntegerPolynomial result(r);
			return result;
                }

		/**
		 * Overloading the * operator. Multiply the polynomial with an integer
		 *
		 * @param factor the integer
		 */
		inline DensedUnivariateIntegerPolynomial operator*(mpz_class factor){
			UnivariateIntegerPolynomial r = p.multiply(factor);
			DensedUnivariateIntegerPolynomial result(r);
                        return result;
		}

		/**
		 * 
		 *
		 * @param b
		 */
		inline DensedUnivariateIntegerPolynomial operator/ (DensedUnivariateIntegerPolynomial& b) {
                        if (b.isZero()) {
                                std::cout << "BPAS: warning, dividend is zero." << std::endl;
                                return b;
                        }
                        if (isZero())
                                return *this;

                        int size = degree() - b.degree() + 1;
                        DensedUnivariateIntegerPolynomial quo(size);
			quo.setVariableName(variableName());
                        quo.zero();

                        DensedUnivariateIntegerPolynomial rem(*this);
                        while (rem.degree() >= b.degree()) {
                                mpz_class lc = rem.leadingCoefficient() / b.leadingCoefficient();
                                int diff = rem.degree() - b.degree();
                                rem.pomopo(1, -lc, b);
                                quo.setCoefficient(diff, lc);
                        }
                        return quo;
                }

		/**
		 * 
		 *
		 * @param
		 */
		inline DensedUnivariateIntegerPolynomial& operator/= (DensedUnivariateIntegerPolynomial b) {
                        if (b.isZero()) {
                                std::cout << "BPAS: warning, dividend is zero." << std::endl;
                                return *this;
                        }
                        if (isZero())
                                return *this;

                        DensedUnivariateIntegerPolynomial rem(*this);
                        zero();
                        while (rem.degree() >= b.degree()) {
                                mpz_class lc = rem.leadingCoefficient() / b.leadingCoefficient();
                                int diff = rem.degree() - b.degree();
                                rem.pomopo(1, -lc, b);
                                setCoefficient(diff, lc);
                        }
                        return *this;
                }

		inline void differentiate(int k) {}


		/**
		 * Overloading the = operator. Check if the polynomial is equivalant to another polynomial.
		 *
		 * @param q the other polynomial
		 */
		inline bool operator==(DensedUnivariateIntegerPolynomial & q){
			return p.isEqual(q.getP());
		}

		/**
		 * Overloading the == operator. Check if the polynomial is equivalant to a constant polynomial (with degree 0).
		 *
		 * @param constant the value of the constant polynomial
		 */
		inline bool operator==(mpz_class constant){
			return (p.getExposedSize() == 1 && p.getCoefficient(0) == constant);
		}

		/**
		 * Overloading the != operator. Check if the polynomial is not equivalant to another polynomial.
		 *
		 * @param q the other polynomial
		 */
		inline bool operator!=(DensedUnivariateIntegerPolynomial & q){
                        return !p.isEqual(q.getP());
                }

		/**
		 *  Overloading the != operator. Check if the polynomial is not equivalant to a constant polynomial (with degree 0).
		 *
		 *  @param constant the value of the constant polynomial
		 */
		inline bool operator!=(mpz_class constant){
                        return !(p.getExposedSize() == 1 && p.getCoefficient(0) == constant);
                }

		/**
		 * Copy the coefficients from another polynomial.
		 *
		 * @param q the polynomial to be copied from
		 */
		inline void copyFrom(DensedUnivariateIntegerPolynomial & q){
			p.copy(q.getP());
		}

		/**
		 * Overload stream operator <<
		 *
		 * @param out the stream object
		 * @param b the polynomial
		 */
		friend std::ostream& operator<< (std::ostream &out, DensedUnivariateIntegerPolynomial & b){
			char v = b.getP()->getVariableName();
			if(b.size() > 0)
				out << b.coefficient(0);
		        for(int i = 1; i < b.size(); ++i)
		                out << " + " << b.coefficient(i) << v << "^" << i;
			out << "\n";
			out << "FINITO\n";
		}

		/**
		 * Prints all of the coefficients of the polynomial.
		 * The variable name is 'x' by default; you may set other names.
		 */
		inline void print(){
			p.print();
		}

		/**
		 * Write all of the coefficients to a file
		 *
		 * @param fileName Name of the file
		 */
		inline void read(std::string fileName){
                        p.read(fileName);
                }

		/**
		 * Write all of the coefficients to a file
		 *
		 * @param fileName Name of the file 
		 */
		inline void write(std::string fileName){
                        p.write(fileName);
                }
		
		/**
		 * Get the degree of the polynomial.
		 */
		inline int degree() {
			return p.getExposedSize() - 1;
		}

		/**
		 * Get the size of the polynomial.
		 */
		inline int size(){
			return p.getExposedSize();
		}

		/**
		 * Get the value of the leading coefficient of the polynomial.
		 */
		inline mpz_class leadingCoefficient(){
			return p.getLeadingCoefficient();
		}

		/**
		 * Get a coefficient of the polynomial at a specified index.
		 *
		 * @param index index of the coefficient
		 */
		inline mpz_class coefficient(int index){
			return p.getCoefficient(index);
		}

		/**
		 * Set the value of a coefficient of the polynomial.
		 *
 		 * @param index index of the coefficient.
		 * @param value value of the coefficient
		 */
		inline void setCoefficient(int index, mpz_class value){
			p.setCoefficient(index, value);
		}

		/**
		 * Get the pointer to the array of the coefficients.
		 */
		inline mpz_class * coefficients(){
			return p.getCoefficients();
		}

		/**
		 * Set the variable's name for the polynomial.
		 *
		 * @param v variable's name
		 */
		inline void setVariableName(char v){
			p.setVariableName(v);
		}

		/**
		 * Return the variable name of the polynomial.
		 */
		inline char variableName(){
			return p.getVariableName();
		}
		
		/**
		 * Check if the polynomial has the root of 0.
		 */
		inline bool isTrailingCoefficientZero(){
			return p.isTrailingCoefficientZero();
		}

		/**
		 * Evaluate the polynomial at an arbitrary point.
		 *
		 * @param point the point
		 */
		inline mpz_class evaluate(mpz_class point){
			return p.evaluate(point);
		}

		/**
		 * Shrink the polynomial by reducing the size (coefficients array).
		 *
		 * @param newSize the new size for the polynomial
		 */
		inline void shrink(int newSize){
			p.shrink(newSize);
		}

		/**
		 * Grow the polynomial by increasing the size (coefficients array).
		 *
		 * @param newSize the new size of the polynomail
		 */
		inline void grow(int newSize){
			p.grow(newSize);
		}

		inline DensedUnivariateIntegerPolynomial operator-(){
			UnivariateIntegerPolynomial r = p.negateOutOfPlace();
                        DensedUnivariateIntegerPolynomial result(r);
			return result;
		}

		/**
		 * Negate the polynomial (all of the coefficients) in place.
		 */
		inline void negate(){
			p.negateInPlace();
		}

		/**
		 * Convert the Polynomial into a big integer (apply Kronecker-Substitution). 
		 * This means evaluating the polynomial at the point 2^b.
		 *
		 * @param b the log of evaluating point
		 */
		inline mpz_class convertIntoBigInteger(int b){
			if (b <= p.maxCoefficientSize())
				throw std::runtime_error("b is too small\n");
			return p.getBigIntegerSigned(b);
		}
};
