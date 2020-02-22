#include "../include/IntegerPolynomial/uzpolynomial.h"
#include "../include/RationalNumberPolynomial/urpolynomial.h"
#include "../include/RationalNumberPolynomial/mrpolynomial.h"
#include "../include/ModularPolynomial/sdmpolynomial.h"
#include "../include/RingPolynomial/dmpolynomial.h"
#include "../include/RingPolynomial/upolynomial.h"
#include "../include/RingPolynomial/mpolynomial.h"
#include "../include/RationalFunction/urationalfunction.h"


mpz_class DenseUnivariateIntegerPolynomial::characteristic(0);
mpz_class DenseUnivariateRationalPolynomial::characteristic(0);
mpz_class SparseMultivariateRationalPolynomial::characteristic(0);
mpz_class SparseMultivariateIntegerPolynomial::characteristic(0);
mpz_class SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::characteristic(101);
template <class Field>
mpz_class DistributedDenseMultivariateModularPolynomial<Field>::characteristic(5);

// bool DenseUnivariateIntegerPolynomial::isPrimeField = 0;
// bool DenseUnivariateRationalPolynomial::isPrimeField = 0;
// bool SparseMultivariateRationalPolynomial::isPrimeField = 0;
// bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::isPrimeField = 0;
// template <class Field>
// bool DistributedDenseMultivariateModularPolynomial<Field>::isPrimeField = 0;
// template<>
// bool SparseUnivariatePolynomial<Integer>::isPrimeField = 0;
// template<>
// bool SparseUnivariatePolynomial<RationalNumber>::isPrimeField = 0;
// template<>
// bool SparseUnivariatePolynomial<ComplexRationalNumber>::isPrimeField = 0;
// template <class Ring>
// bool SparseUnivariatePolynomial<Ring>::isPrimeField = 0;
// template <class Ring>
// bool SparseMultivariatePolynomial<Ring>::isPrimeField = 0;
// template <class Polynomial, class Field>
// bool UnivariateRationalFunction<Polynomial, Field>::isPrimeField = 0;

// bool DenseUnivariateIntegerPolynomial::isSmallPrimeField = 0;
// bool DenseUnivariateRationalPolynomial::isSmallPrimeField = 0;
// bool SparseMultivariateRationalPolynomial::isSmallPrimeField = 0;
// bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::isSmallPrimeField = 0;
// template <class Field>
// bool DistributedDenseMultivariateModularPolynomial<Field>::isSmallPrimeField = 0;
// template<>
// bool SparseUnivariatePolynomial<Integer>::isSmallPrimeField = 0;
// template<>
// bool SparseUnivariatePolynomial<RationalNumber>::isSmallPrimeField = 0;
// template<>
// bool SparseUnivariatePolynomial<ComplexRationalNumber>::isSmallPrimeField = 0;
// template <class Ring>
// bool SparseUnivariatePolynomial<Ring>::isSmallPrimeField = 0;
// template <class Ring>
// bool SparseMultivariatePolynomial<Ring>::isSmallPrimeField = 0;
// template <class Polynomial, class Field>
// bool UnivariateRationalFunction<Polynomial, Field>::isSmallPrimeField = 0;

// bool DenseUnivariateIntegerPolynomial::isComplexField = 0;
// bool DenseUnivariateRationalPolynomial::isComplexField = 0;
// bool SparseMultivariateRationalPolynomial::isComplexField = 0;
// bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::isComplexField = 0;
// template <class Field>
// bool DistributedDenseMultivariateModularPolynomial<Field>::isComplexField = 0;
// template<>
// bool SparseUnivariatePolynomial<Integer>::isComplexField = 0;
// template<>
// bool SparseUnivariatePolynomial<RationalNumber>::isComplexField = 0;
// template<>
// bool SparseUnivariatePolynomial<ComplexRationalNumber>::isComplexField = 0;
// template <class Ring>
// bool SparseUnivariatePolynomial<Ring>::isComplexField = 0;
// template <class Ring>
// bool SparseMultivariatePolynomial<Ring>::isComplexField = 0;
// template <class Polynomial, class Field>
// bool UnivariateRationalFunction<Polynomial, Field>::isComplexField = 0;

