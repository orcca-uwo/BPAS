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
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


