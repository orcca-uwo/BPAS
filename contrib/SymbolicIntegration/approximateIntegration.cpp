#include <bpas.h>
#include "RationalFunction/multiprecision_rootfinding.h"
#include "RationalFunction/rationalfunction_complexrationalnumberordering.h"

using namespace std;


//bool epsilonEqual(mpq_class a, mpq_class b, int &prec){
//	mpq_class t,s;
//	t = a - b;
//	t = abs(t);
//	s = mpq_class(1);
//	s /= mpq_class(2)^(prec);
//	if (mpq_get_d(t.get_mpq_t()) < mpq_get_d(s.get_mpq_t()))
//		return true;
//	else
//		return false;
//}

bool epsilonEqual(const RationalNumber& a, const RationalNumber& b, int &prec) {
	RationalNumber t,s;
	t = a - b;
	t = abs(t);
	s = RN(1);
	s /= RN(2)^(prec);
	if (t.get_d() < s.get_d())
		return true;
	else
		return false;
}

bool epsilonEqual(CRN &a, CRN &b, int prec){
	CRN t;
	RN m,s;
	double d;
	t = a - b;
	m = t.realPart();
	m *= m;
	s = t.imaginaryPart();
	s *= s;
	m += s;
	s = RN(1);
	s /= RN(2)^(prec);
	d = mpq_get_d(m.get_mpq_t());
	d = sqrt(d);
	if (d < mpq_get_d(s.get_mpq_t()))
		return true;
	else
		return false;
}


inline int factorial(int x) {
  return (x == 1 ? x : x * factorial(x - 1));
}

void distinctRoots(SparseUnivariatePolynomial<RN> &den, int prec, vector<CRN> *roots, vector<int> *mult){
	/* finds distinct roots of input polynomial den with multiplicities */
	/* note: a scalar version of roots_m should perhpas be written      */
	
	int mp(0);
	vector<int> m;
	vector<CRN> r;
	vector< vector<CRN> > rVec;
	vector< SparseUnivariatePolynomial<RN> > U;

	U.push_back(den);
	rVec = _rootsMultiprecision<SparseUnivariatePolynomial<RN>,RN>(U,prec);
	r = rVec.at(0);
	
	//for (int i=0;i<r.size();i++) {
	//	cout << "r[" << i << "] = " << r.at(i) << endl;
	//}
	
	sort(r.begin(),r.end(),ComplexRationalNumberOrdering(prec));
	
	// remove duplicates and count multiplicity //
	for (int i=r.size()-1; i>0; i--) {
		mp = 1;
		int j = r.size()-1-i;
		while (i > 0 && epsilonEqual(r.at(i),r.at(i-1),prec)) {
			r.erase(r.end()-j-1);
			mp++;
			i--;
		}
		m.insert(m.begin(),mp);
		if (i==1) {
			mp = 1;
			m.insert(m.begin(),mp);
		}
	}
	
	*mult = m;
	*roots = r;
}
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


