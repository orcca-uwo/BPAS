#ifndef _BPAS_H_
#define _BPAS_H_

#include "ModularPolynomial/sdmpolynomial.h"
#include "IntegerPolynomial/uzpolynomial.h"
#include "RingPolynomial/upolynomial.h"
#include "RingPolynomial/mpolynomial.h"
#include "RingPolynomial/dmpolynomial.h"
//#include "RingPolynomial/dupolynomial.h"
#include "Interval/interval.h"
#include "RationalNumberPolynomial/urpolynomial.h"
//#include "RationalNumberPolynomial/urpolynomial-altarr.hpp"
#include "RationalNumberPolynomial/mrpolynomial.h"
#include "SubResultantChain/subresultantchain.hpp"
#include "RegularChain/rationalregularchain.h"
#include "RegularChain/zerodimensionalregularchain.hpp"
#include "RegularChain/regularchain.hpp"
#include "TriangularSet/triangularset.hpp"
//#include "../src/TriangularSet/triangularset.cpp"
#include "RationalFunction/urationalfunction.h"
#include "RationalFunction/unumericalpolynomial.h"
#include "ExpressionTree/ExpressionTree.hpp"
#include "Symbol/Symbol.hpp"

#include "Parser/bpas_parser.h"
//#include "ring.h"

#if defined(WITH_MAPLE) && WITH_MAPLE
#include "MapleInterface/MapleInterfaceStream.hpp"
#endif

typedef SmallPrimeFieldDistributedDenseMultivariateModularPolynomial SFDDMMP;
typedef DenseUnivariateIntegerPolynomial DUZP;
typedef DenseUnivariateRationalPolynomial DUQP;
typedef SparseMultivariateRationalPolynomial SMQP;
typedef SparseMultivariateIntegerPolynomial SMZP;
typedef RationalRegularChain QRC;
typedef RationalNumber RN;
typedef ComplexRationalNumber CRN;

template <typename T>
struct Ring {
	typedef SparseUnivariatePolynomial<T> SUP;
	typedef SparseMultivariatePolynomial<T> SMP;
};

template <typename T>
struct Field {
	typedef DistributedDenseMultivariateModularPolynomial<T> DDMMP;
//	typedef DenseUnivariatePolynomial<T> DUFP;
};

std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > bivariateRationalSubresultantChain (SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial>&, SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial>&);
//long int SmallPrimeField::prime = 2147483647;
//long int SmallPrimeField::characteristic = 2147483647;
//mpz_class BigPrimeField::prime ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
//mpz_class GeneralizedFermatPrimeField::prime ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
//unsigned long long int GeneralizedFermatPrimeField::r = 9223372054034644992ULL;
//mpz_class BigPrimeField::characteristic ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
//mpz_class GeneralizedFermatPrimeField::characteristic ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
//int GeneralizedFermatPrimeField::k = 8;

#endif
