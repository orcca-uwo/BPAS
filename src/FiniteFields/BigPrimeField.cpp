
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

using std::endl;
using std::cout;

//mpz_class BigPrimeField::prime ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
mpz_class BigPrimeField::characteristic ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);

// bool BigPrimeField::isPrimeField = 1;
// bool BigPrimeField::isSmallPrimeField = 0;
// bool BigPrimeField::isComplexField = 0;
mpz_class BigPrimeField::prime = BigPrimeField::characteristic;


BigPrimeField::BigPrimeField () : a(0) {

}

BigPrimeField::BigPrimeField (mpz_class _a) {
	_a = _a%prime;
	while (_a < 0){
		_a = _a + prime;
	}
	a = _a%prime;
}

BigPrimeField::BigPrimeField (long int _a) {
	a=_a%prime;
	while (a < 0){
		a = a + prime;
	}
//  a = a%prime;
}

BigPrimeField::BigPrimeField (const BigPrimeField& c) {
	a = c.a;
}

BigPrimeField::BigPrimeField (const Integer& c) {
	a = c.get_mpz();
	while (a < 0){
		a = a + prime;
	}
	a = a%prime;
}

BigPrimeField::BigPrimeField (const RationalNumber& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const SmallPrimeField& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const GeneralizedFermatPrimeField& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const DenseUnivariateIntegerPolynomial& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const DenseUnivariateRationalPolynomial& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const SparseUnivariatePolynomial<Integer>& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const SparseUnivariatePolynomial<RationalNumber>& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField::BigPrimeField (const SparseUnivariatePolynomial<ComplexRationalNumber>& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

template <class Ring>
BigPrimeField::BigPrimeField (const SparseUnivariatePolynomial<Ring>& c) {
	std::cerr << "Cannot convert input to BigPrimeField! " << std::endl;
	exit(1);
}

BigPrimeField* BigPrimeField::BPFpointer(BigPrimeField* b) {
	return b;
}

BigPrimeField* BigPrimeField::BPFpointer(RationalNumber* a) {
	std::cout << "BPAS error, try to cast pointer to Rational Number to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

BigPrimeField* BigPrimeField::BPFpointer(SmallPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to BigPrimeField to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

BigPrimeField* BigPrimeField::BPFpointer(GeneralizedFermatPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to GeneralizedFermatPrimeField to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

mpz_class BigPrimeField::Prime() const {
	return prime;
}

mpz_class BigPrimeField::number() const {
	return mpz_class(a);
}

void BigPrimeField::whichprimefield() {
	cout << "BigPrimeField" << endl;
}

BigPrimeField BigPrimeField::unitCanonical(BigPrimeField* u, BigPrimeField* v) const {
	if (isZero()) {
		if (u != NULL) {
			*u = 1;
		}
		if (v != NULL) {
			*v = 1;
		}
		return BigPrimeField(0);
	} else {
		if (u != NULL) {
			*u = *this;
		}
		if (v != NULL) {
			*v = this->inverse();
		}
		return BigPrimeField(1);
	}
}

BigPrimeField& BigPrimeField::operator= (const BigPrimeField& c) {
	if (this != &c) {
		a = c.a;
	}
	return *this;
}

BigPrimeField& BigPrimeField::operator= (long int k) {
	a = k;
	while (a < 0){
		a = a + prime;
	}
	a = a%prime;
	return *this;
}

BigPrimeField& BigPrimeField::operator= (const mpz_class& k) {
	a=k;
	while (a < 0){
		a = a + prime;
	}
	a = a%prime;
	return *this;
}

BigPrimeField BigPrimeField::euclideanDivision(const BigPrimeField& b, BigPrimeField* q) const {
	//TODO
	std::cerr << "BigPrimeField::euclideanDivision not yet implemented" << std::endl;
	exit(1);
	return *this;
}

BigPrimeField BigPrimeField::extendedEuclidean(const BigPrimeField& b, BigPrimeField* s, BigPrimeField* t) const {
	//TODO
	std::cerr << "BigPrimeField::extendedEuclidean not yet implemented" << std::endl;
	exit(1);
	return *this;
}

BigPrimeField BigPrimeField::quotient(const BigPrimeField& b) const {
	//TODO
	std::cerr << "BigPrimeField::quotient not yet implemented" << std::endl;
	exit(1);
	return *this;
}

BigPrimeField BigPrimeField::remainder(const BigPrimeField& b) const {
	//TODO

	std::cerr << "BigPrimeField::remainder not yet implemented" << std::endl;
	exit(1);
	return *this;
}

