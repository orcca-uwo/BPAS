
#include "Ring/Integer.hpp"
#include "Ring/RationalNumber.hpp"
#include "Ring/ComplexRationalNumber.hpp"
#include "FiniteFields/SmallPrimeField.hpp"
#include "FiniteFields/BigPrimeField.hpp"
#include "FiniteFields/GeneralizedFermatPrimeField.hpp"
#include "IntegerPolynomial/uzpolynomial.h"
#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h"

#include <iostream>

mpz_class Integer::characteristic(0);
RingProperties Integer::properties = EUCLIDEAN_DOMAIN;		

// bool Integer::isPrimeField = 0;
// bool Integer::isSmallPrimeField = 0;
// bool Integer::isComplexField = 0;

//// TODO //////////////////////////////
// Implement properly all the conversions.
////////////////////////////////////////

Integer::Integer () { 
	_m = 0; 
}

Integer::Integer (int a) {
	_m = a;
}

Integer::Integer (const mpz_t& a) : _m(a) {

}

Integer::Integer (const mpz_class& a) : _m(a) {

}

Integer::Integer (const Integer& a) {
	_m = a._m;
}

Integer::Integer (const RationalNumber& a) {
	if (a.get_den() == 1)
		*this = a.get_num();
	else {
		std::cout << "BPAS error, try to construct a rational number to Integer class." << std::endl;
		free((int*)-1);
		exit(1);
	}
}

Integer::Integer (const ComplexRationalNumber& a) {
	if (a.imaginaryPart() != 0) {
		std::cout << "BPAS error, try to construct a complex number to Integer class." << std::endl;
		exit(1);
	}
	mpq_class e = a.realPart().get_mpq();
	if (e.get_den() != 1) {
		std::cout << "BPAS error, try to construct a complex number to Integer class." << std::endl;
		exit(1);
	}
	else
		*this = e.get_num();
}

Integer::Integer (const SmallPrimeField& a) {
	mpz_t z;
	long long int n = a.number();
	mpz_import(z, 1, 1, sizeof(n), 0, 0, &n);
	mpz_class w(z);
	*this = w;
}

Integer::Integer (const BigPrimeField& a) {
	//mpz_t z;

	//mpz_import(z, 1, -1, sizeof(a.a), 0, 0, &a.a);
	mpz_class w = a.number();
	*this = w;
}

Integer::Integer (const GeneralizedFermatPrimeField& a) {
	*this = a.number();
}

Integer::Integer (const DenseUnivariateIntegerPolynomial& a) {
	if (a.degree() > 0) {
		std::cout << "BPAS error, try to construct a univariate polynomial over Z to Integer class." << std::endl;
		exit(1);
	}
	*this = a.coefficient(0);
}

Integer::Integer (const DenseUnivariateRationalPolynomial& a) {
	if (a.degree() > 0) {
		std::cout << "BPAS error, try to construct a univariate polynomial over Q to Integer class." << std::endl;
		exit(1);
	}
	mpq_class e = a.coefficient(0).get_mpq();
	if (e.get_den() == 1)
		*this = e.get_num();
	else {
		std::cout << "BPAS error, try to construct a univariate polynomila over Q to Integer class." << std::endl;
		exit(1);
	}
}

Integer::Integer (const SparseUnivariatePolynomial<Integer>& a) {
    if (a.degree() > 0) {
            std::cout << "BPAS error, try to construct a univariate polynomial to Integer class." << std::endl;
            exit(1);
    }

    *this = a.coefficient(0);
}

Integer::Integer (const SparseUnivariatePolynomial<RationalNumber>& a) {
	if (a.degree() > 0) {
		std::cout << "BPAS error, try to construct a univariate polynomial to Integer class." << std::endl;
		exit(1);
	}

	RationalNumber e = a.coefficient(0);
	*this = Integer(e);
}

Integer::Integer (const SparseUnivariatePolynomial<ComplexRationalNumber>& a) {
        if (a.degree() > 0 ) {
                std::cout << "BPAS error, try to construct a univariate polynomial to Integer class." << std::endl;
                exit(1);
        }

        ComplexRationalNumber e = a.coefficient(0);
        *this = Integer(e);
}

template <class Ring>
Integer::Integer (const SparseUnivariatePolynomial<Ring>& a) {
	std::cout << "BPAS error, try to construct a univariate polynomial to Integer class." << std::endl;
	exit(1);
}

Integer Integer::unitCanonical(Integer* u, Integer* v) const {
	if (*this >= 0) {
		if (u != NULL) {
			*u = 1;
		}
		if (v != NULL) {
			*v = 1;
		}
		return *this;
	} else {
		if (u != NULL) {
			*u = -1;
		}
		if (v != NULL) {
			*v = -1;
		}
		return -(*this);
	}
}

Integer& Integer::operator= (const Integer& a) {
    if (*this != a) {
    	_m = a._m;
    }

    return *this;
}

Integer Integer::euclideanDivision(const Integer& b, Integer* q) const {
	//TODO
	std::cerr << "Integer::euclideanDivision not yet implemented" << std::endl;
	exit(1);
	return *this;
}

Integer Integer::extendedEuclidean(const Integer& b, Integer* s, Integer* t) const {
	//TODO
	std::cerr << "Integer::extendedEuclidean not yet implemented" << std::endl;
	exit(1);
	return *this;
}

Integer Integer::quotient(const Integer& b) const {
	//TODO
	std::cerr << "Integer::quotient not yet implemented" << std::endl;
	exit(1);
	return *this;
}

Integer Integer::remainder(const Integer& b) const {
	//TODO

	std::cerr << "Integer::remainder not yet implemented" << std::endl;
	exit(1);
	return *this;
}
	
