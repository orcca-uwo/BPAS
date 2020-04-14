#ifndef _UNUMERICALPOLYNOMIAL
#define _UNUMERICALPOLYNOMIAL
#include <complex>
#include "../RingPolynomial/upolynomial.h"
#include "../ring.h"
#include "urationalfunction.h"

/**
 * A univariate polynomial with numerical coefficients represented sparsely.
 * The template determines the type of the numerical coefficients.
 */
template<class NumericalType>
class SparseUnivariateDoublePolynomial
{
	private:
		std::vector<NumericalType> coef;
		std::vector<int> degs;

		template<class UnivariatePolynomialOverRealField>
		NumericalType evaluateRealPolynomial(UnivariatePolynomialOverRealField& A, NumericalType& a) {
			// need a check here that input is a real field
			NumericalType out;
			out = NumericalType(A.leadingCoefficient().get_d());
			for (int i=A.degree().get_si()-1; i>-1; i--){
				out *= a;
				if (A.coefficient(i) != 0){
					out += NumericalType(A.coefficient(i).get_d());
				}
			}
			return out;
		}
		std::complex<double> evaluateComplexPolynomial(SparseUnivariatePolynomial<ComplexRationalNumber>& A, std::complex<double>& a) {
			ComplexRationalNumber value;
			value = A.leadingCoefficient();
			std::complex<double> out(value.realPart().get_d(),value.imaginaryPart().get_d());
			for (int i=A.degree().get_si()-1; i>-1; i--){
				out *= a;
				if (!A.coefficient(i).isZero()){
					value = A.coefficient(i);
					out += std::complex<double>(value.realPart().get_d(),value.imaginaryPart().get_d());
				}
			}
			return out;
		}

	public:
		SparseUnivariateDoublePolynomial(DenseUnivariateRationalPolynomial& A) {
			// need a check here that input is a real field
			NumericalType value(0);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				value = A.coefficient(i).get_d();
				if (!(value == 0.0)){
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		SparseUnivariateDoublePolynomial(SparseUnivariatePolynomial<RationalNumber>& A) {
			// need a check here that input is a real field
			NumericalType value(0);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				value = A.coefficient(i).get_d();
				if (!(value == 0.0)){
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		// construct a SUNP by evaluating the coefficients of an SUP<DUQP> at a NumericalType
		SparseUnivariateDoublePolynomial(SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial>& A, NumericalType& a) {
			// need a check here that input is a real field
			DenseUnivariateRationalPolynomial term;
			NumericalType value(0);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				term = A.coefficient(i);
				if (!term.isZero()){
					value = evaluateRealPolynomial<DenseUnivariateRationalPolynomial>(term,a);
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		// construct a SUNP by evaluating the coefficients of an SUP<SUP<RN>> at a NumericalType
		SparseUnivariateDoublePolynomial(SparseUnivariatePolynomial< SparseUnivariatePolynomial<RationalNumber> >& A, NumericalType a) {
			// need a check here that input is a real field
			SparseUnivariatePolynomial<RationalNumber> term;
			NumericalType value(0);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				term = A.coefficient(i);
				if (!term.isZero()){
					value = evaluateRealPolynomial< SparseUnivariatePolynomial<RationalNumber> >(term,a);
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		SparseUnivariateDoublePolynomial(SparseUnivariatePolynomial<ComplexRationalNumber>& A) {
			NumericalType value(0);
			ComplexRationalNumber cTemp;
			for (int i=0; i<A.degree().get_si()+1; i++) {
				cTemp = A.coefficient(i);
				value = NumericalType(cTemp.realPart().get_d(),cTemp.imaginaryPart().get_d());
				if (!(value == std::complex<double>(0))){
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		// construct a SUNP by evaluating the coefficients of an SUP<SUP<CRN>> at a complex<double>
		SparseUnivariateDoublePolynomial(SparseUnivariatePolynomial< SparseUnivariatePolynomial<ComplexRationalNumber> >& A, std::complex<double> a) {
			SparseUnivariatePolynomial<ComplexRationalNumber> term;
			NumericalType value(0);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				term = A.coefficient(i);
				if (!term.isZero()){
					value = evaluateComplexPolynomial(term,a);
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		~SparseUnivariateDoublePolynomial() {};
		int degree() {
			int out(degs.back());
			return out;
		}
		NumericalType leadingCoefficient() {
			NumericalType out(coef.back());
			return out;
		}
		NumericalType evaluate(double t){
			NumericalType out;
			int d = coef.size() - 1;
			if (d < 0) { return out; }
			int e = degs.at(d) - 1;
			out = coef.at(d);
			d--;
			for (int i = e; i > -1; --i) {
				out *= NumericalType(t);
				if (d > -1){
					if (i == degs.at(d)) {
						out += coef.at(d);
						d--;
					}
				}
			}
			return out;
		}
		std::complex<double> evaluate(std::complex<double> t){
			std::complex<double> out;
			int d = coef.size() - 1;
			if (d < 0) { return out; }
			int e = degs.at(d) - 1;
			out = std::complex<double>(coef.at(d));
			d--;
			for (int i = e; i > -1; --i) {
				out *= t;
				if (d > -1){
					if (i == degs.at(d)) {
						out += std::complex<double>(coef.at(d));
						d--;
					}
				}
			}
			return out;
		}
};

// this should be redone with complexMPF converted into a class, so that it can have nice constructors and behave
// like a number not a struct.
/**
 * A univariate polynomial with multi-precision complex coefficients represented sparsely.
 */
class SparseUnivariateMPComplexPolynomial
{
	private:
		std::vector<complexMPF> coef;
		std::vector<int> degs;

	public:
		SparseUnivariateMPComplexPolynomial(DenseUnivariateRationalPolynomial& A, int prec) {
			complexMPF value;
			mpc_init2(value.c,4*prec);
			mpq_t zero;
			mpq_init(zero);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				mpc_set_q(value.c,A.coefficient(i).get_mpq_t(),zero);
				if (!(mpc_eq_zero(value.c))){
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}
		SparseUnivariateMPComplexPolynomial(SparseUnivariatePolynomial<RationalNumber>& A, int prec) {
			complexMPF value;
			mpc_init2(value.c,4*prec);
			mpq_t zero;
			mpq_init(zero);
			for (int i=0; i<A.degree().get_si()+1; i++) {
				mpc_set_q(value.c,A.coefficient(i).get_mpq_t(),zero);
				if (!(mpc_eq_zero(value.c))){
					coef.push_back(value);
					degs.push_back(i);
				}
			}
		}

		complexMPF evaluate(complexMPF t, int prec){
			complexMPF out;
			mpc_init2(out.c,4*prec);
			int d = coef.size() - 1;
			if (d >= 0) {
				int e = degs.at(d) - 1;
				mpc_set(out.c,coef.at(d).c);
				d--;
				for (int i = e; i > -1; --i) {
					mpc_mul(out.c,out.c,t.c);
					if (d > -1){
						if (i == degs.at(d)) {
							mpc_add(out.c,out.c,coef.at(d).c);
							d--;
						}
					}
				}
			}
			return out;
		}
};

#endif
