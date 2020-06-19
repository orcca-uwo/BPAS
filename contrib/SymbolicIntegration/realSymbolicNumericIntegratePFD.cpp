#include <bpas.h>
#include "RationalFunction/rationalfunction_euclideanmethods.h"
#include "RationalFunction/rationalfunction_integrationpostprocessing.h"
#include "../../include/RationalFunction/complexdouble.h"

using namespace std;

// approximateIntegration.cpp //
void distinctRoots(SparseUnivariatePolynomial<RN> &den, int prec, vector<CRN> *roots, vector<int> *mult);
//bool epsilonEqual(mpq_class a, mpq_class b, int &prec);
bool epsilonEqual(const RationalNumber& a, const RationalNumber& b, int &prec);
bool epsilonEqual(CRN &a, CRN &b, int prec);

// floatIntegrate.cpp //
SparseUnivariatePolynomial<CRN> simplifyComplexRFTerm(CRN &c, SparseUnivariatePolynomial<CRN> *A);
void normalize(SparseUnivariatePolynomial<CRN> *A, SparseUnivariatePolynomial<CRN> *B);
void normalize(SparseUnivariatePolynomial<RN> *A, SparseUnivariatePolynomial<RN> *B);

inline int factorial(int x) {
  return (x == 1 ? x : x * factorial(x - 1));
}

extern void snIntPFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING);

void realSymbolicNumericIntegratePFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING){

	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	cout << "[realSymbolicNumericIntegratePFD (snIntPFD): Symbolic-Numeric Integration using BPAS and MPSolve]" << endl;
	cout << "[Integration method: Hermite Reduction, Partial Fraction Decomposition w/ Numerical Rootfinding]" << endl;
	cout << "Starting..." << endl;
	
	if (PROFILING){	
		cout << "snIntPFD" << endl;
		cout << "--------------------------------------" << endl;
		startTimer(&start);
	}
	
	snIntPFD(f,P,g,lg,Lg,atn,Atn,prec,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "--------------------------------------" << endl;
		cout << "floatInt runtime: " << elapsed << " s" << endl;
	}
	cout << endl;
}

CRN computeResidue(CRN &root, SparseUnivariatePolynomial<CRN> &G, SparseUnivariatePolynomial<CRN> &dH) {
	
	CRN residue;
	
	residue = G.evaluate(root);
	residue /= dH.evaluate(root);
	
	return residue;
}

/*SparseUnivariatePolynomial<CRN> simplifyComplexRFTerm(CRN &c, SparseUnivariatePolynomial<CRN> *A){
	CRN cstar(c.conjugate());
	SparseUnivariatePolynomial<CRN> Astar;
	SparseUnivariatePolynomial<CRN> term;

	Astar = conjugate(*A);
	term = c*Astar;
	term += c.conjugate()* *A;
	*A *= Astar;
	
	return term;
}

void normalize(SparseUnivariatePolynomial<CRN> *A, SparseUnivariatePolynomial<CRN> *B){
	SparseUnivariatePolynomial<CRN> temp;
	temp = A->gcd(*B);
	*A /= temp;
	*B /= temp;
}

void normalize(SparseUnivariatePolynomial<RN> *A, SparseUnivariatePolynomial<RN> *B){
	SparseUnivariatePolynomial<RN> temp;
	temp = A->gcd(*B);
	*A /= temp;
	*B /= temp;
}*/

SparseUnivariatePolynomial<RN> _CRNPolyRealPart(SparseUnivariatePolynomial<CRN> &A){
	SparseUnivariatePolynomial<RN> output;
	output.setVariableName(A.variable());
	RN rTemp;
	for (int i=0; i<A.degree()+1; i++){
		rTemp = A.coefficient(i).realPart();
		if (!rTemp.isZero())
			output.setCoefficient(i,rTemp);
	}
	return output;
}

SparseUnivariatePolynomial<RN> _CRNPolyImaginaryPart(SparseUnivariatePolynomial<CRN> &A){
	SparseUnivariatePolynomial<RN> output;
	output.setVariableName(A.variable());
	RN rTemp;
	for (int i=0; i<A.degree()+1; i++){
		rTemp = A.coefficient(i).imaginaryPart();
		if (!rTemp.isZero())
			output.setCoefficient(i,rTemp);
	}
	return output;
}

void snIntPFD_BEA(SparseUnivariatePolynomial<CRN> &G, SparseUnivariatePolynomial<CRN> &H, SparseUnivariatePolynomial<CRN> &dH, vector<CRN> roots){
	
	UnivariateRationalFunction<SparseUnivariatePolynomial<CRN>,CRN> e(G,dH);
	UnivariateRationalFunction<SparseUnivariatePolynomial<CRN>,CRN> de(G,dH);
	UnivariateRationalFunction<SparseUnivariatePolynomial<CRN>,CRN> dalpha(H,dH);
	de.differentiate();
	vector<double> x;
	vector<double> xi; // to hold values of error term at a given x
	vector<double> Xi; // to hold values of max error at each x
	vector<CRN> alpha;
	vector<double> e_alpha_rp;
	vector<double> e_alpha_ip;
	vector<double> de_alpha_rp;
	vector<double> de_alpha_ip;
	vector<double> delta_alpha_rp;
	vector<double> delta_alpha_ip;
	double dTemp;
	double et1_rp;
	double et1_ip;
	double et2_rp;
	double et2_ip;
	double x_minus_alpha_rp;
	double x_minus_alpha_ip;
	CRN cTemp;
	
	// This loop evaluates delta_alpha = H(alpha)/dH(alpha), e(alpha) and de(alpha), which are
	// independent of the locations x at which the sum is evaluated
	for (int i=0;i<roots.size();i++) {
		cTemp = dalpha.evaluate(roots.at(i));
		if (!cTemp.isZero()) {
			alpha.push_back(roots.at(i));
			delta_alpha_rp.push_back(cTemp.realPart().get_d());
			delta_alpha_ip.push_back(cTemp.imaginaryPart().get_d());
			cTemp = e.evaluate(roots.at(i));
			e_alpha_rp.push_back(cTemp.realPart().get_d());
			e_alpha_ip.push_back(cTemp.imaginaryPart().get_d());
			cTemp = de.evaluate(roots.at(i));
			de_alpha_rp.push_back(cTemp.realPart().get_d());
			de_alpha_ip.push_back(cTemp.imaginaryPart().get_d());
			x.push_back(roots.at(i).realPart().get_d());
		}
	}
	
	// This loop computes the products and quotients of each term in the sum, puts the result in xi,
	// and then sums the results to compute a the summation at x, which is placed in Xi. The max
	// value of the elements in Xi is the estimate of the backward error
	for (int i=0;i<x.size();i++) {
		for (int j=0;j<alpha.size();j++) {
			x_minus_alpha_rp = x.at(i) - alpha.at(j).realPart().get_d();
			x_minus_alpha_ip = -alpha.at(j).imaginaryPart().get_d();
			if 	(x_minus_alpha_rp == 0) {
				dTemp = abs(e_alpha_rp.at(j)*delta_alpha_rp.at(j) - e_alpha_ip.at(j)*delta_alpha_ip.at(j));
				dTemp *= 1e16;
				dTemp = sqrt(dTemp);
				cout << "|x-alpha| > " << dTemp << " for alpha = " << alpha.at(j).realPart().get_d() << endl;
				x.at(i) += dTemp;
				dTemp = x.at(i) + 2*dTemp;
				//x.insert(x.begin()+i+1, dTemp);
				x_minus_alpha_rp = x.at(i) - alpha.at(j).realPart().get_d();
			}
			if (alpha.at(j).imaginaryPart() == 0) {
				et1_rp = de_alpha_rp.at(j)/x_minus_alpha_rp;
				et2_rp = e_alpha_rp.at(j)/(x_minus_alpha_rp*x_minus_alpha_rp);
				et1_rp = et1_rp + et2_rp;
				et1_rp = et1_rp*delta_alpha_rp.at(j);
				xi.push_back(et1_rp);
			}
			else {
				complexDivide(de_alpha_rp.at(j), de_alpha_ip.at(j), x_minus_alpha_rp, x_minus_alpha_ip, &et1_rp, &et1_ip);
				complexDivide(e_alpha_rp.at(j), e_alpha_ip.at(j), x_minus_alpha_rp*x_minus_alpha_rp -
					 x_minus_alpha_ip*x_minus_alpha_ip, 2*x_minus_alpha_rp*x_minus_alpha_ip, &et2_rp, &et2_ip);
				complexAdd(et1_rp, et1_ip, et2_rp, et2_ip, &et1_rp, &et1_ip);
				//complexMultiply(et1_rp, et1_ip, delta_alpha_rp, delta_alpha_ip, &et1_rp, &et1_ip);
				et1_rp = et1_rp*delta_alpha_rp.at(j) - et1_ip*delta_alpha_ip.at(j);
				xi.push_back(2*et1_rp);
			}
		}
		dTemp = 0;
		for (int j=0;j<xi.size();j++) {
		 	dTemp += xi.at(j);
		}
		Xi.push_back(dTemp);
	}
	cout << abs(*max_element(Xi.begin(),Xi.end())) << endl;

	// this needs to be rewritten to avoid computing polys except for the rfs e(x) and e'(x), and to compute
	// in terms of doubles as much as possible.
	/*double dTemp;
	RN rTemp;
	CRN one;
	CRN alpha;
	CRN deltaAlpha;
	CRN cE;
	CRN cdE;
	CRN cG;
	CRN cdG;
	CRN cH;
	CRN cdH;
	CRN cddH;
	CRN cTemp;
	SparseUnivariatePolynomial<RN> aNum;
	SparseUnivariatePolynomial<RN> bNum;
	SparseUnivariatePolynomial<RN> aDen;
	SparseUnivariatePolynomial<RN> bDen;
	SparseUnivariatePolynomial<RN> suprnTemp;
	SparseUnivariatePolynomial<CRN> linPoly; // polynomial of the form x - a
	SparseUnivariatePolynomial<CRN> cNum;
	SparseUnivariatePolynomial<CRN> cDen;
	SparseUnivariatePolynomial<CRN> dG;
	SparseUnivariatePolynomial<CRN> ddH;
	UnivariateRationalFunction<SparseUnivariatePolynomial<CRN>,CRN> urfTemp;
	dG = G;
	dG.differentiate(1);
	ddH = dH;
	ddH.differentiate(1);
	int n = roots.size();
	vector<double> Xi;
	vector<RN> realParts;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<CRN>,CRN> > xi;
	
	one.set(1,0);
	linPoly.setCoefficient(1,one);
	
	for (int i=0;i<n;i++) {
		// this needs to be edited to detect complex roots and multiply the error factor by 2 in that case
		alpha = roots.at(i);
		linPoly.setCoefficient(0,-alpha);
		cDen = linPoly*linPoly; // denominator = (x - alpha)^2
		cG = G.evaluate(alpha);
		cdG = dG.evaluate(alpha);
		cH = H.evaluate(alpha);
		cdH = dH.evaluate(alpha);
		cddH = ddH.evaluate(alpha);
		cTemp = cdG*cdH;
		cTemp -= cG*cddH;
		cTemp *= cdH*cdH;
		cNum = cTemp*linPoly; // E'(alpha)*(x - alpha)
		cTemp = cG/cdH; // E(alpha)
		cNum += cTemp;
		cTemp = -cH/cdH; // deltaAlpha
		cNum *= cTemp;
		rTemp.zero();
		if (alpha.imaginaryPart() == rTemp) {
			cNum *= 2;
			/*suprnTemp = _CRNPolyRealPart(cNum);
			urfTemp.setNumerator(suprnTemp);
			suprnTemp = _CRNPolyRealPart(cDen);
			urfTemp.setDenominator(suprnTemp);*/
		//}
		/*else {
			aNum = _CRNPolyRealPart(cNum);
			bNum = _CRNPolyRealPart(cNum);
			aDen = _CRNPolyRealPart(cDen);
			bDen = _CRNPolyRealPart(cDen);
			suprnTemp = aNum*aDen;
			suprnTemp += bNum*bDen;
			suprnTemp *= 2;
			urfTemp.setNumerator(suprnTemp);
			suprnTemp = aDen*aDen;
			suprnTemp *= bDen*bDen;
			urfTemp.setDenominator(suprnTemp);
		}*/
		//urfTemp.set(cNum,cDen);
		//xi.push_back(urfTemp);
		// and then an efficient way to graph the sum of the vector of univariate rational functions
	//}
	//cout << "xi.size() = " << xi.size() << endl;
	/*for (int i=0;i<n;i++) {
		rTemp = roots.at(i).realPart();
		realParts.push_back(rTemp);
	}*/
	/*for (int i=0;i<n;i++) {
		dTemp = 0;
		for (int j=0;j<n;j++) {
			cTemp.set(roots.at(j).realPart(),one.imaginaryPart());
			dTemp += xi.at(i).evaluate(cTemp).realPart().get_d();
		}
		Xi.push_back(dTemp);
	}
	cout << *max_element(Xi.begin(),Xi.end()) << endl;*/
}

void snIntPFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING){
// Hermite Reduction and Partial Fraction Integration of rational function with rational coefficients //
	
	Symbol variable(f.variable());
	CRN c;
	CRN cTemp;
	vector<int> mult;
	vector<CRN> roots;
	//SparseUnivariatePolynomial<RN> Q;
	SparseUnivariatePolynomial<RN> R;
	SparseUnivariatePolynomial<RN> Num(f.numerator());
	SparseUnivariatePolynomial<RN> Den(f.denominator());
	//SparseUnivariatePolynomial<RN> NewDen;
	SparseUnivariatePolynomial<CRN> cDen;
	SparseUnivariatePolynomial<CRN> cLinPoly; // poly of the form a*x + b
	SparseUnivariatePolynomial<CRN> cNum;
	SparseUnivariatePolynomial<CRN> cPolyTemp;
	SparseUnivariatePolynomial<CRN> cdH;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f2(f);
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> h(f);
	
	R.setVariableName(variable);
	//Q.setVariableName(variable);
	cNum.setVariableName(variable);
	cDen.setVariableName(variable);
	cdH.setVariableName(variable);
	cLinPoly.setVariableName(variable);
	cPolyTemp.setVariableName(variable);
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	// clear containers for results
	P->zero();
	g->clear();
	lg->clear();
	Lg->clear();
	atn->clear();
	Atn->clear();
	
	// might be an error in lazyPseudoDivide and pseudoDivide
	_euclideanDivide<SparseUnivariatePolynomial<RN>,RN>(Num,Den,P,&R);
	if (!(P->isZero()))
		*P = P->integral();
	
	normalize(&R,&Den);
	f2.set(R,Den);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	f2.hermiteReduce(g, &h);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "\thermiteReduce\t\t\t" << elapsed << " s" << endl;
	}
	
	Num = h.numerator();
	Den = h.denominator();
	
	if (PROFILING){
		startTimer(&start);
	}
	
	// precision correction handled by _rootsMultiprecision called by distinctRoots
	distinctRoots(Den,prec,&roots,&mult);
	
	//for (int i=0;i<roots.size();i++) {
	//	cout << "roots[" << i << "] = " << roots.at(i) << " has multiplicity " << mult.at(i) << endl;
	//}
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "\tdistinctRoots\t\t\t" << elapsed << " s" << endl;
	}
	
	//Num = R;
	//Num /= Den.leadingCoefficient();
	
	// remove complex conjugate from roots list //
	for (int i=roots.size()-1; i>=0; i--) {
		int j = roots.size()-1-i;
		if (roots.at(i).imaginaryPart() != 0) {
			roots.erase(roots.end()-j-1);
			mult.erase(mult.end()-j-1);
			i--;
		}
	}
	
	cNum = RNPolyToCRNPoly(Num);
	cDen = RNPolyToCRNPoly(Den);
	cdH = cDen;
	cdH.differentiate(1);
	
	RN rnTemp;
	SparseUnivariatePolynomial<RN> suprnTemp;
	suprnTemp.setVariableName(variable);
	//UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> urfTemp;
	//urfTemp.setVariableName(variable);
	//vector< SparseUnivariatePolynomial<RN> > G;
	//vector< SparseUnivariatePolynomial<RN> > Lg;
	//vector< SparseUnivariatePolynomial<RN> > Atn;
	//vector< SparseUnivariatePolynomial<RN> > Atn2;
	//vector<RN> lg;
	//vector<RN> atn;
	
	RN zero(mpq_class(0));
	
	if (PROFILING){
		startTimer(&start);
	}
	
	c.one();
	cLinPoly.setCoefficient(1,c);
		
	for (int i=0; i<roots.size(); i++){
	
		cTemp = computeResidue(roots.at(i), cNum, cdH);
	
		if (cTemp.realPart() != 0 && epsilonEqual(cTemp.realPart(),zero,prec))
			cTemp.setRealPart(zero);
		if (cTemp.imaginaryPart() != 0 && epsilonEqual(cTemp.imaginaryPart(),zero,prec))
			cTemp.setImaginaryPart(zero);
		if (roots.at(i).imaginaryPart() == 0) {
			cLinPoly.setCoefficient(0,-roots.at(i));
			lg->push_back((RN)cTemp.realPart());
			Lg->push_back(CRNPolyRealPart(cLinPoly));
		}
		else {
			if (cTemp.imaginaryPart() != 0){
				RN a;
				RN b;
				a = roots.at(i).realPart();
				b = roots.at(i).imaginaryPart();
				// consider implementing with RN polys
				c.set(-a,0);
				cLinPoly.setCoefficient(0,c);
				c.set(-b,0);
				cLinPoly /= c; // cLinPoly = (x - a)/(-b)
				rnTemp = RN(cTemp.imaginaryPart()*2);
				atn->push_back(rnTemp);
				suprnTemp = CRNPolyRealPart(cLinPoly);
				Atn->push_back(suprnTemp);
				c.one();
				cLinPoly.setCoefficient(1,c);
			}
			if (cTemp.realPart() != 0){
				cLinPoly.setCoefficient(0,-roots.at(i));
				cPolyTemp = cLinPoly;
				cLinPoly.setCoefficient(0,-roots.at(i).conjugate());
				cPolyTemp *= cLinPoly;
				rnTemp = RN(cTemp.realPart());
				lg->push_back(rnTemp);
				suprnTemp = CRNPolyRealPart(cPolyTemp);
				Lg->push_back(suprnTemp);
			}
		}
	}
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "\tIntegration Step\t\t\t" << elapsed << " s" << endl;
		startTimer(&start);
	}
	
	snIntPFD_BEA(cNum, cDen, cdH, roots);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "\tBackward Error Analysis\t\t\t" << elapsed << " s" << endl;
	}
	
	//UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> temp;
	//vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	//int i(0);
	//while (i<G.size()) {
	//	temp.set(G.at(i),G.at(i+1));
	//	g.push_back(temp);
	//	i += 2;
	//}
	//f.printIntegral(*P, g, lg, Lg, atn, Atn);

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


