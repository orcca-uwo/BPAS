
#include "Ring/Integer.hpp"
#include "Ring/RationalNumber.hpp"
#include "Ring/ComplexRationalNumber.hpp"
#include "FiniteFields/SmallPrimeField.hpp"
#include "FiniteFields/BigPrimeField.hpp"
#include "FiniteFields/GeneralizedFermatPrimeField.hpp"
#include "IntegerPolynomial/uzpolynomial.h"
#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h" 

mpz_class ComplexRationalNumber::characteristic(0);
RingProperties ComplexRationalNumber::properties = COMPLEX_FIELD;

// bool ComplexRationalNumber::isPrimeField = 0;
// bool ComplexRationalNumber::isSmallPrimeField = 0;
// bool ComplexRationalNumber::isComplexField = 1;

ComplexRationalNumber::ComplexRationalNumber () : a(0), b(0) {

}

ComplexRationalNumber::ComplexRationalNumber (const mpq_class& _a, const mpq_class& _b) {
	a = _a;
	b = _b;
}

ComplexRationalNumber::ComplexRationalNumber (const ComplexRationalNumber& c) {
	a = c.a;
	b = c.b;
}

ComplexRationalNumber::ComplexRationalNumber(int _a, int _b, int _c, int _d) {
	a = mpq_class(_a, _b);
	b = mpq_class(_c, _d);
}

ComplexRationalNumber::ComplexRationalNumber (const Integer& c) {
	a = mpq_class(c.get_mpz());
	b = 0;
}

ComplexRationalNumber::ComplexRationalNumber (const RationalNumber& c) {
	a = c.get_mpq();
	b = 0;
}

ComplexRationalNumber::ComplexRationalNumber (const SmallPrimeField& c) {
	std::cerr << "Cannot convert input to ComplexRationalNumber! " << std::endl;
	exit(1);
}

ComplexRationalNumber::ComplexRationalNumber (const BigPrimeField& c) {
	std::cerr << "Cannot convert input to ComplexRationalNumber! " << std::endl;
	exit(1);
}

ComplexRationalNumber::ComplexRationalNumber (const GeneralizedFermatPrimeField& c) {
	std::cerr << "Cannot convert input to ComplexRationalNumber! " << std::endl;
	exit(1);
}

ComplexRationalNumber::ComplexRationalNumber (const DenseUnivariateIntegerPolynomial& c) {
	if (c.degree() > 0) {
	        std::cout << "BPAS error, try to construct a univariate polynomial to ComplexRationalNumber class." << std::endl;
	        exit(1);
	}
	a = c.coefficient(0).get_mpz();
	b = 0;
}

ComplexRationalNumber::ComplexRationalNumber (const DenseUnivariateRationalPolynomial& c) {
        if (c.degree() > 0) {
                std::cout << "BPAS error, try to construct a univariate polynomial to ComplexRationalNumber class." << std::endl;
                exit(1);
        }
    a = c.coefficient(0).get_mpq();
	b = 0;
}

ComplexRationalNumber::ComplexRationalNumber (const SparseUnivariatePolynomial<Integer>& c) {
        if (c.degree() > 0) {
                std::cout << "BPAS error, try to construct a univariate polynomial to ComplexRationalNumber class." << std::endl;
                exit(1);
        }

    a = c.coefficient(0).get_mpz();
	b = 0;
}

ComplexRationalNumber::ComplexRationalNumber (const SparseUnivariatePolynomial<RationalNumber>& c) {
        if (c.degree() > 0) {
                std::cout << "BPAS error, try to construct a univariate polynomial to ComplexRationalNumber class." << std::endl;
                exit(1);
        }

    a = c.coefficient(0).get_mpq();
	b = 0;
}

ComplexRationalNumber::ComplexRationalNumber (const SparseUnivariatePolynomial<ComplexRationalNumber>& c) {
        if (c.degree() > 0) {
                std::cout << "BPAS error, try to construct a univariate polynomial to ComplexRationalNumber class." << std::endl;
                exit(1);
        }

        *this = c.coefficient(0);
}

template <class Ring>
ComplexRationalNumber::ComplexRationalNumber (const SparseUnivariatePolynomial<Ring>& c) {
        std::cout << "BPAS error, try to construct a univariate polynomial to ComplexRationalNumber class." << std::endl;
        exit(1);
}

ComplexRationalNumber ComplexRationalNumber::unitCanonical(ComplexRationalNumber* u, ComplexRationalNumber* v) const {
	if(isZero()) {
		if (u != NULL) {
			*u = 1;
		}
		if (v != NULL) {
			*v = 1;
		}
		return ComplexRationalNumber(0);
	}

	if (u != NULL) {
		*u = *this;
	}
	if (v != NULL) {
		*v = inverse();
	}

	return ComplexRationalNumber(1);
}


ComplexRationalNumber& ComplexRationalNumber::operator= (const ComplexRationalNumber& c) {
	if (this != &c) {
		a = c.a;
		b = c.b;
	}
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::operator= (const mpq_class& k) {
    a = k;
    b = 0;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::operator= (int k) {
    a = k;
    b = 0;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::setRealPart (const RationalNumber& r) {
	a = r.get_mpq();
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::setRealPart (const mpq_class& k) {
    a = k;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::setRealPart (int k) {
    a = k;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::setImaginaryPart (const RationalNumber& r) {
	a = r.get_mpq();
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::setImaginaryPart (const mpq_class& k) {
    b = k;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::setImaginaryPart (int k) {
    b = k;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::set (const RationalNumber& ka, const RationalNumber& kb) {
	a = ka.get_mpq();
	b = kb.get_mpq();
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::set (const mpq_class& ka, const mpq_class& kb) {
	a = ka;
    b = kb;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::set (const mpq_class& ka, int kb) {
	a = ka;
    b = kb;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::set (int ka, const mpq_class& kb) {
	a = ka;
    b = kb;
	return *this;
}

ComplexRationalNumber& ComplexRationalNumber::set (int ka, int kb) {
	a = ka;
    b = kb;
	return *this;
}


ComplexRationalNumber ComplexRationalNumber::euclideanDivision(const ComplexRationalNumber& b, ComplexRationalNumber* q) const {
	//TODO
	std::cerr << "ComplexRationalNumber::euclideanDivision not yet implemented" << std::endl;
	exit(1);
	return *this;
}

ComplexRationalNumber ComplexRationalNumber::extendedEuclidean(const ComplexRationalNumber& b, ComplexRationalNumber* s, ComplexRationalNumber* t) const {
	//TODO
	std::cerr << "ComplexRationalNumber::extendedEuclidean not yet implemented" << std::endl;
	exit(1);
	return *this;
}

ComplexRationalNumber ComplexRationalNumber::quotient(const ComplexRationalNumber& b) const {
	//TODO
	std::cerr << "ComplexRationalNumber::quotient not yet implemented" << std::endl;
	exit(1);
	return *this;
}

ComplexRationalNumber ComplexRationalNumber::remainder(const ComplexRationalNumber& b) const {
	//TODO

	std::cerr << "SmallPrimeField::remainder not yet implemented" << std::endl;
	exit(1);
	return *this;
}

void ComplexRationalNumber::print(std::ostream &out) const {
	if (this->b == 0) {
		out << this->a;
	}
	else if (this->a == 0) {
		if (abs(this->b) == 1) {
			if (this->b == 1)
				out << "I";
			else
				out << "-I";
		}
		else {
			if (this->b < 0)
				out << "-I*" << abs(this->b);
			else
				out << "I*" << this->b;
		}
	}
	else {
		out << "(" << this->a;
		if (abs(this->b) == 1) {
			if (this->b == 1)
				out << "+I";
			else
				out << "-I";
		}
		else {
			if (this->b < 0)
				out << "-I*" << abs(this->b);
			else
				out << "+I*" << this->b;
		}
		out << ")";
	}
}
