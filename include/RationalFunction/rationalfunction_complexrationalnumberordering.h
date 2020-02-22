#ifndef _RATIONALFUNCTION_COMPLEXRATIONALNUMBERORDERING_H
#define _RATIONALFUNCTION_COMPLEXRATIONALNUMBERORDERING_H
#include <complex>

struct CompareByRealThenReverseImaginary {
	bool operator() (ComplexRationalNumber a, ComplexRationalNumber b) {
		if (!(a.realPart() == b.realPart()))
			return (a.realPart() < b.realPart());
		else {
			if (a.imaginaryPart() == 0)
				return (true);
			else if (b.imaginaryPart() == 0)
				return (false);
			else
				return (a.imaginaryPart() > b.imaginaryPart());
		}
	}
	
	bool operator() (std::pair<ComplexRationalNumber,ComplexRationalNumber> a, std::pair<ComplexRationalNumber,ComplexRationalNumber> b) {
		if (!(a.first == b.first)) {
			if (!(a.first.realPart() == b.first.realPart()))
				return (a.first.realPart() < b.first.realPart());
			else {
				if (a.first.imaginaryPart() == 0)
					return (true);
				else if (b.first.imaginaryPart() == 0)
					return (false);
				else
					return (a.first.imaginaryPart() > b.first.imaginaryPart());
			}
		}
		else {
			if (!(a.second.realPart() == b.second.realPart()))
				return (a.second.realPart() < b.second.realPart());
			else {
				if (a.second.imaginaryPart() == 0)
					return (true);
				else if (b.second.imaginaryPart() == 0)
					return (false);
				else
					return (a.second.imaginaryPart() > b.second.imaginaryPart());
			}
		}
	}
	
	bool operator() (std::pair<ComplexRationalNumber,int> a, std::pair<ComplexRationalNumber,int> b) {
		if (!(a.first == b.first)) {
			if (!(a.first.realPart() == b.first.realPart()))
				return (a.first.realPart() < b.first.realPart());
			else {
				if (a.first.imaginaryPart() == 0)
					return (true);
				else if (b.first.imaginaryPart() == 0)
					return (false);
				else
					return (a.first.imaginaryPart() > b.first.imaginaryPart());
			}
		}
		else
			return (a.second < b.second);
	}
};

/**
 * A class to allow ordering (using sort command) of ComplexRationalNumber numbers.
 * These numbers are obtained from 
 * numerical rootfinding with precision 'prec', which determines a radius of equality.
 */
class ComplexRationalNumberOrdering {

private:
	int prec;
	
	//bool epsilonEqual(RationalNumber &a, RationalNumber &b, int prec);

	//bool epsilonEqual(ComplexRationalNumber &a, ComplexRationalNumber &b, int prec);

	bool complexRationalNumberOrdering(ComplexRationalNumber &A, ComplexRationalNumber &B, int &prec);

public:

	/**
	 * Construct a ComplexRationalNumberOrdering with precision p.
	 */
	ComplexRationalNumberOrdering(int p) : prec(p) {}
		
	static bool epsilonEqual(RationalNumber a, RationalNumber b, int &prec);
	
	bool operator()(ComplexRationalNumber A, ComplexRationalNumber B) {
		return complexRationalNumberOrdering(A,B,prec);
	}

};

/**
 * A class to allow ordering (using sort command) of complex<double> numbers.
 * Ordering is with the precision 'prec', which determines a radius of equality.
 */
class ComplexDoubleOrdering {

private:
	int prec;
	
	//bool epsilonEqual(RationalNumber &a, RationalNumber &b, int prec);

	//bool epsilonEqual(ComplexRationalNumber &a, ComplexRationalNumber &b, int prec);

	bool complexDoubleOrdering(std::complex<double> &A, std::complex<double> &B, int &prec);

public:
	ComplexDoubleOrdering(int p) {
		if (p > 52)
			prec = 52;
		else
			prec = p;
	}
		
	static bool epsilonEqual(double a, double b, int &prec);
	
	bool operator()(std::complex<double> A, std::complex<double> B) {
		return complexDoubleOrdering(A,B,prec);
	}

};
#endif
