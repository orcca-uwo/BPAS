#include "../../include/ring.h"
#include "RationalFunction/rationalfunction_complexrationalnumberordering.h"
#include <complex>


/*bool ComplexRationalNumberOrdering::epsilonEqual(ComplexRationalNumber &a, ComplexRationalNumber &b, int prec){
	ComplexRationalNumber t;
	RationalNumber m,s;
	double d;
	t = a - b;
	m = t.realPart();
	m *= m;
	s = t.imaginaryPart();
	s *= s;
	m += s;
	s = (RationalNumber)1;
	s /= (RationalNumber)2^(prec);
	d = mpq_get_d(m.get_mpq_t());
	d = sqrt(d);
	if (d < mpq_get_d(s.get_mpq_t()))
		return true;
	else
		return false;
}*/

bool ComplexRationalNumberOrdering::complexRationalNumberOrdering(ComplexRationalNumber &A, ComplexRationalNumber &B, int &prec){
	if (!(epsilonEqual(A.realPart(),B.realPart(),prec)))
		return (A.realPart() < B.realPart());
	else {
		if (A.imaginaryPart() == 0)
			return (true);
		else if (B.imaginaryPart() == 0)
			return (false);
		else
			return (A.imaginaryPart() > B.imaginaryPart());
	}
}

// bool ComplexRationalNumberOrdering::epsilonEqual(mpq_class a, mpq_class b, int &prec){
// 	mpq_class t,s;
// 	t = a - b;
// 	t = abs(t);
// 	s = RationalNumber(mpq_class(1));
// 	s /= RationalNumber(mpq_class(2))^(prec);
// 	if (mpq_get_d(t.get_mpq_t()) < mpq_get_d(s.get_mpq_t()))
// 		return true;
// 	else
// 		return false;
// }

bool ComplexRationalNumberOrdering::epsilonEqual(RationalNumber a, RationalNumber b, int &prec){
	RationalNumber t,s;
	t = a - b;
	t = abs(t);
	s = RationalNumber(1);
	s /= RationalNumber(2)^(prec);
	if (mpq_get_d(t.get_mpq_t()) < mpq_get_d(s.get_mpq_t()))
		return true;
	else
		return false;
}

bool ComplexDoubleOrdering::complexDoubleOrdering(std::complex<double> &A, std::complex<double> &B, int &prec){
	if (!(epsilonEqual(A.real(),B.real(),prec)))
		return (A.real() < B.real());
	else {
		if (epsilonEqual(A.imag(),0.0,prec))
			return (true);
		else if (epsilonEqual(B.imag(),0.0,prec))
			return (false);
		else
			return (A.imag() > B.imag());
	}
}

bool ComplexDoubleOrdering::epsilonEqual(double a, double b, int &prec){
	double t,s;
	t = a - b;
	if (t < 0.0)
		t = -t;
	s = pow(2.0,-prec);
	if (t < s)
		return true;
	else
		return false;
}


