#ifndef _FRACTION_H_
#define _FRACTION_H_

#include <iostream>
#include "BPASFieldOfFractions.hpp"
#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"
#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include <algorithm>
#include <vector>

/**
 * A field of fractions templated by an arbitrary BPASGCDDomain.
 */
template <class Domain>
class Fraction: public BPASFieldOfFractions<Domain> {

	private:
		Domain den;//denominator
		Domain num;//numerator
		


	public:

		mpz_class characteristic;
		/**
		 * Construct the zero fraction function
		 *
		 * @param
		 **/ 
		Fraction<Domain>(){
			den.one();
			num.zero();
			Domain e;
			characteristic = e.characteristic;
		};
		/**
		 * Copy constructor
		 *
		 * @param b: A rational function
		 **/ 
		Fraction<Domain>(const Fraction<Domain> &b): num(b.num),den(b.den){
			characteristic = b.characteristic;
		}

		/**
		 * constructor with two parameter
		 *
		 * @param a: the numerator
		 * @param b: the denominator
		 **/
		Fraction<Domain>(Domain a, Domain b){
			if (b.isZero()) {
		    std::cout << "BPAS error: denominator is zero from Fraction<Domain>" << std::endl;
		    exit(1);
		  }
		  num = a;
		  den = b;
		  Domain e;
		  characteristic = e.characteristic;
		}

		/**
		 * Deconstructor for fraction
		 *
		 * @param
		 **/
		~Fraction<Domain> () {}


		void setNumerator(const Domain& b);
		void setDenominator(const Domain& b);

		void set(const Domain& a, const Domain& b);

		Domain numerator() const;
		Domain denominator() const;

		bool operator!=(const Fraction<Domain> &b) const;
		bool operator==(const Fraction<Domain> &b) const;

		Fraction<Domain> operator+(const Fraction<Domain> &b) const;
		Fraction<Domain>& operator+=(const Fraction<Domain> &b);
		
		Fraction<Domain> operator-(const Fraction<Domain> &b) const;
		Fraction<Domain>& operator-=(const Fraction<Domain> &b);
		
		Fraction<Domain> operator*(const Fraction<Domain> &b) const;
		Fraction<Domain>& operator*=(const Fraction<Domain> &b);

		Fraction<Domain> operator/(const Fraction<Domain> &b) const;
		Fraction<Domain>& operator/=(const Fraction<Domain> &b);

		Fraction<Domain> operator-() const;

		 /**
         * Overload operator ^
         * replace xor operation by exponentiation
         *
         * @param e: The exponentiation, e > 0
         **/
		Fraction<Domain> operator^(long long int e) const;
		/**
         * Overload operator ^=
         * replace xor operation by exponentiation
         *
         * @param e: The exponentiation, e > 0
         **/
		Fraction<Domain>& operator^=(long long int e);

		Fraction<Domain> inverse() const;

		bool isZero() const;
		void zero();
		bool isOne() const;
		void one();
		bool isNegativeOne() const;
		void negativeOne();
		int isConstant() const;
		Fraction<Domain> unitCanonical(Fraction<Domain>* u = NULL, Fraction<Domain>*v = NULL) const;

		/**
		 * Overload operator =
		 *
		 * @param b: A rational function
		 **/
		Fraction<Domain>& operator=(const Fraction<Domain>& b);

		void print(std::ostream& ostream) const;

		Fraction<Domain> quotient(const Fraction<Domain>& b) const;

		//not implemented yet
		Fraction<Domain> remainder(const Fraction<Domain>& b) const;
		//not implemented yet
		Fraction<Domain> operator%(const Fraction<Domain>& b) const;
		//not implemented yet
		Fraction<Domain>& operator%=(const Fraction<Domain>& b);


		void canonicalize();

		void differentiate();


		void normalize();
		//RefineReturn<Domain> Refine(Domain a, Domain e, Domain b, Domain f);




		ExpressionTree convertToExpressionTree() const {
		std::cerr << "UnivariateRationalFunction::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
			//TODO
		return ExpressionTree();
		}

    	Fraction<Domain> gcd(const Fraction<Domain>& b) const {
    		//std::cerr << "UnivariateRationalFunction::gcd NOT YET IMPLEMENTED" << std::endl;
    		//if both are 0 then return 0
    		Fraction<Domain> ret;
    		if(isZero() == true &&b.isZero() == true){
    			ret.zero();
    			
    		}
    		//otherwise is one
    		else{
    			ret.one();	
    		}

    		return ret;
    	}

    	Factors<Fraction<Domain>> squareFree() const {
    		std::cerr << "BPAS ERROR: Fraction::squareFree NOT YET IMPLEMENTED" << std::endl;
    		return Factors<Fraction<Domain>>(*this);
    	}

    	Fraction<Domain> euclideanSize() const {
    		//TODO
    		//if is zero, throw error exit(1)
    		//otherwise return 0

    		if(isZero() == true){
    			std::cerr << "in Fraction<Domain>, zero does not have a euclidean size" << std::endl;
    			exit(1);
    		}
    		else{
    			Fraction<Domain> ret;
    			ret.zero();
    			return ret;
    		}
    	}

    	Fraction<Domain> euclideanDivision(const Fraction<Domain>& b, Fraction<Domain>* q = NULL) const {
    		//TODO
    		//return type is remainder,
    		//q is quo,
    		//it always exact division
    		Fraction<Domain> ret;
    		ret.zero();
    		if(q!=NULL){
    			Fraction<Domain> curr(num,den);
    			*q = curr/b;
    		}

    		return ret;
    	}

    	Fraction<Domain> extendedEuclidean(const Fraction<Domain>& b, Fraction<Domain>* s = NULL,Fraction<Domain>* t = NULL) const {
    		std::cerr << "UnivariateRationalFunction::extendedEuclidean NOT YET IMPLEMENTED" << std::endl;
    		//TODO
    		//
    		return *this;
    	}

};





#endif
