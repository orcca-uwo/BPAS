
#include "Ring/Integer.hpp"
#include "Ring/RationalNumber.hpp"
#include "Ring/ComplexRationalNumber.hpp"
#include "FiniteFields/SmallPrimeField.hpp"
#include "FiniteFields/BigPrimeField.hpp"
#include "FiniteFields/GeneralizedFermatPrimeField.hpp"
#include "IntegerPolynomial/uzpolynomial.h"
#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h"
#include "RationalNumberPolynomial/mrpolynomial.h"

#include <iostream>

// bool RationalNumber::isPrimeField = 1;
// bool RationalNumber::isSmallPrimeField = 0;
// bool RationalNumber::isComplexField = 0;

RationalNumber::RationalNumber () {
	_m = 0;
}

RationalNumber::RationalNumber (int a, int b) {
	_m = a;
	_m /= b;
}

RationalNumber::RationalNumber (const std::string& digits, int base) : _m(digits, base) {

}


RationalNumber::RationalNumber (const mpq_t& q) : _m(q) {

}

RationalNumber::RationalNumber (const mpq_class& a) : _m(a) {

}

RationalNumber::RationalNumber (const mpz_class& a, const mpz_class& b ) : _m(a, b) {

}


RationalNumber::RationalNumber (const RationalNumber& a) : _m(a._m) {

}

RationalNumber::RationalNumber (const Integer& a) : _m(a.get_mpz()) {

}

RationalNumber::RationalNumber (const ComplexRationalNumber& a) {
	if (a.imaginaryPart() != 0) {
		std::cout << "BPAS error, try to construct a complex number to RationalNumber class." << std::endl;
		exit(1);
	}
	else {
		*this = a.realPart();
	}
}

RationalNumber::RationalNumber (const SmallPrimeField& a) {
	mpz_t z;
	long long int n = a.number();
	mpz_import(z, 1, -1, sizeof(n), 0, 0, &n);
	mpz_class w(z);
	_m = w;
}

RationalNumber::RationalNumber (const BigPrimeField& a) {
	std::cerr << "Cannot convert input to RationalNumber! " << std::endl;
	exit(1);
}

RationalNumber::RationalNumber (const GeneralizedFermatPrimeField& a) {
	std::cerr << "Cannot convert input to RationalNumber! " << std::endl;
	exit(1);
}

RationalNumber::RationalNumber (const DenseUnivariateIntegerPolynomial& a) {
	if (a.degree() > 0) {
		std::cout << "BPAS error, try to construct a univariate polynomial to RationalNumber class." << std::endl;
		exit(1);
	}
	*this = RationalNumber(a.coefficient(0));
}

RationalNumber::RationalNumber (const DenseUnivariateRationalPolynomial& a) {
	if (a.degree() > 0) {
        std::cout << "BPAS error, try to construct a univariate polynomial to RationalNumber class." << std::endl;
        exit(1);
    }
	*this = RationalNumber(a.coefficient(0));
}

RationalNumber::RationalNumber (const SparseUnivariatePolynomial<Integer>& a) {
    if (a.degree() > 0) {
        std::cout << "BPAS error, try to construct a univariate polynomial to RationalNumber class." << std::endl;
        exit(1);
    }

    *this = RationalNumber(a.coefficient(0));
}

RationalNumber::RationalNumber (const SparseUnivariatePolynomial<RationalNumber>& a) {
    if (a.degree() > 0) {
            std::cout << "BPAS error, try to construct a univariate polynomial to RationalNumber class." << std::endl;
            exit(1);
    }

    *this = RationalNumber(a.coefficient(0));
}

RationalNumber::RationalNumber (const SparseUnivariatePolynomial<ComplexRationalNumber>& a) {
    if (a.degree() > 0) {
            std::cout << "BPAS error, try to construct a univariate polynomial to RationalNumber class." << std::endl;
            exit(1);
    }

    *this = RationalNumber(a.coefficient(0));
}

template <class Ring>
RationalNumber::RationalNumber (const SparseUnivariatePolynomial<Ring>& a) {
        std::cout << "BPAS error, try to construct a univariate polynomial to RationalNumber class." << std::endl;
        exit(1);
}

RationalNumber::RationalNumber (const SparseMultivariateRationalPolynomial& a) {
    if (a.isConstant() == 0) {
            std::cout << "BPAS error, try to construct a multivariate polynomial to RationalNumber class." << std::endl;
            exit(1);
    }

    *this = RationalNumber(a.leadingCoefficient());
}

RationalNumber* RationalNumber::RNpointer(RationalNumber* a) {
	return a;
}

RationalNumber* RationalNumber::RNpointer(SmallPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to SmallPrimeField to pointer to Rational Number" << std::endl;
	exit(1);
}

RationalNumber* RationalNumber::RNpointer(BigPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to BigPrimeField to pointer to Rational Number" << std::endl;
	exit(1);
}

RationalNumber* RationalNumber::RNpointer(GeneralizedFermatPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to GeneralizedFermatPrimeField to pointer to Rational Number" << std::endl;
	exit(1);
}

RationalNumber& RationalNumber::set (int a, int b) {
	*this = a;
	*this /= b;
	return *this;
}

RationalNumber RationalNumber::unitCanonical(RationalNumber* u, RationalNumber* v) const {
	if (isZero()) {
		if (u != NULL) {
			*u = 1;
		}
		if (v != NULL) {
			*v = 1;
		}
		return 0;
	} else {
		if (u != NULL) {
			*u = *this;
		}
		if (v != NULL) {
			*v = inverse();
		}
		return 1;
	}
}


RationalNumber& RationalNumber::operator= (const RationalNumber& a) {
	if (this != &a) {
		_m = a._m;
	}
	return *this;
}


RationalNumber RationalNumber::euclideanDivision(const RationalNumber& b, RationalNumber* q) const {
	//TODO
	std::cerr << "RationalNumber::euclideanDivision not yet implemented" << std::endl;
	exit(1);
	return *this;
}

RationalNumber RationalNumber::extendedEuclidean(const RationalNumber& b, RationalNumber* s, RationalNumber* t) const {
	//TODO
	std::cerr << "RationalNumber::extendedEuclidean not yet implemented" << std::endl;
	exit(1);
	return *this;
}

RationalNumber RationalNumber::quotient(const RationalNumber& b) const {
	//TODO
	std::cerr << "RationalNumber::quotient not yet implemented" << std::endl;
	exit(1);
	return *this;
}

RationalNumber RationalNumber::remainder(const RationalNumber& b) const {
	//TODO

	std::cerr << "RationalNumber::remainder not yet implemented" << std::endl;
	exit(1);
	return *this;
}
