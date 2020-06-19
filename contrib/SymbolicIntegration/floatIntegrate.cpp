#include <bpas.h>
#include "RationalFunction/rationalfunction_euclideanmethods.h"
#include "RationalFunction/rationalfunction_integrationpostprocessing.h"

using namespace std;

// approximateIntegration.cpp //
void distinctRoots(SparseUnivariatePolynomial<RN> &den, int prec, vector<CRN> *roots, vector<int> *mult);
//bool epsilonEqual(mpq_class a, mpq_class b, int &prec);
bool epsilonEqual(const RationalNumber& a, const RationalNumber& b, int &prec);
bool epsilonEqual(CRN &a, CRN &b, int prec);

inline int factorial(int x) {
  return (x == 1 ? x : x * factorial(x - 1));
}

extern void floatInt(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING);

void floatIntegrate(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING){

	unsigned long long start;
	float elapsed = 0;
	
	cout << "[floatIntegrate (fInt): Symbolic-Numeric Integration using BPAS and MPSolve]" << endl;
	cout << "[Integration method: Brute Force Partial Fraction Decomposition w/ Numerical Rootfinding]" << endl;
	cout << "Starting..." << endl;
	
	if (PROFILING){	
		cout << "floatInt" << endl;
		cout << "--------------------------------------" << endl;
		startTimer(&start);
	}
	
	floatInt(f,P,g,lg,Lg,atn,Atn,prec,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "--------------------------------------" << endl;
		cout << "floatInt runtime: " << elapsed << " s" << endl;
	}
	cout << endl;
}

void differentiate(SparseUnivariatePolynomial<CRN> *num, SparseUnivariatePolynomial<CRN> *den){
	/* Destructive rational function differentiation */
	// input a/d
	// output a'*d-a*d'/d^2 = a'*e-a*f/d*e; e=d/g, f=d'/g
	
	SparseUnivariatePolynomial<CRN> D(*den); // D = d
	SparseUnivariatePolynomial<CRN> dD(*den);
	SparseUnivariatePolynomial<CRN> temp;
	
	dD.differentiate(1);	// dD = d'
	temp = D.gcd(dD);		// temp = g
	D /= temp;				// D = e
	dD /= temp;				// dD = f
	temp = -*num;			// temp = -a
	temp *= dD;				// temp = -a*f
	dD = *num;				// dD = a
	dD.differentiate(1);	// dD = a'
	dD *= D;				// dD = a'*e
	temp += dD;				// temp = a'*e-a*f
	D *= *den; 				// D = d*e
	
	//temp /= D.leadingCoefficient();
	//D /= D.leadingCoefficient();
	*num = temp;
	*den = D;
	//canonicalize(num,den);
}

CRN computeResidue(CRN &root, int mult, int index, SparseUnivariatePolynomial<CRN> &R, SparseUnivariatePolynomial<CRN> &denMinor) {
	// It would be more efficient to compute all the residues for a given //
	// root at once, so that derivatives do not need to be recomputed.    //
	
	int i = mult - index;
	SparseUnivariatePolynomial<CRN> den(denMinor);
	SparseUnivariatePolynomial<CRN> num(R);
	//UnivariateRationalFunction<SparseUnivariatePolynomial<CRN>,CRN> t(R,denMinor);
	CRN residue;
	
	residue.one();
	if (i > 0) {
		for (int k=0; k<i; k++) {
			//differentiate(&num,&den);
			//t.differentiate();
		}
	residue *= factorial(i);
	residue = residue.inverse();
	}
	
	//residue *= t.numerator().evaluate(root);
	//residue /= t.denominator().evaluate(root);
	residue *= num.evaluate(root);
	residue /= den.evaluate(root);
	
	return residue;
}

SparseUnivariatePolynomial<CRN> simplifyComplexRFTerm(CRN &c, SparseUnivariatePolynomial<CRN> *A){
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
}

void floatInt(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING){
// Partial Fraction Integration of rational function with rational coefficients //
	
	Symbol variable(f.variable());
	CRN cTemp;
	vector<int> mult;
	vector<CRN> roots;
	//SparseUnivariatePolynomial<RN> Q;
	SparseUnivariatePolynomial<RN> R;
	SparseUnivariatePolynomial<RN> Num(f.numerator());
	SparseUnivariatePolynomial<RN> Den(f.denominator());
	SparseUnivariatePolynomial<RN> NewDen;
	SparseUnivariatePolynomial<CRN> cNewDen;
	SparseUnivariatePolynomial<CRN> cLinPoly; // poly of the form x + c
	SparseUnivariatePolynomial<CRN> cNum;
	SparseUnivariatePolynomial<CRN> cPolyTemp;
	
	R.setVariableName(variable);
	//Q.setVariableName(variable);
	NewDen.setVariableName(variable);
	cNum.setVariableName(variable);
	cNewDen.setVariableName(variable);
	cLinPoly.setVariableName(variable);
	cPolyTemp.setVariableName(variable);
	
	// Profiling variables
	unsigned long long start;
	unsigned long long start_is;
	float elapsed = 0;
	float elapsed_minor = 0;
	float elapsed_res = 0;
	float elapsed_pp = 0;
	
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
	
	Num = R;
	Num /= Den.leadingCoefficient();
	
	// compute the new denominator (includes numerical error, if any) //
	
	if (PROFILING){
		startTimer(&start);
	}
		
	cNewDen.one();
	cTemp.set(1,0);
	cLinPoly.setCoefficient(1,cTemp);
	for (int i=0; i<roots.size(); i++) {
		cLinPoly.setCoefficient(0,-roots.at(i));
		cNewDen *= cLinPoly^mult.at(i);
	}
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		cout << "\tnewDenominator\t\t\t" << elapsed << " s" << endl;
	}

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
	
	RN rnTemp;
	SparseUnivariatePolynomial<RN> suprnTemp;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> urfTemp;
	suprnTemp.setVariableName(variable);
	urfTemp.setVariableName(variable);
	//vector< SparseUnivariatePolynomial<RN> > G;
	//vector< SparseUnivariatePolynomial<RN> > Lg;
	//vector< SparseUnivariatePolynomial<RN> > Atn;
	//vector< SparseUnivariatePolynomial<RN> > Atn2;
	//vector<RN> lg;
	//vector<RN> atn;
	
	RN zero(mpq_class(0));
	
	if (PROFILING){
		cout << "\tintegrationStep\t\t\t" << endl;
		cout << "--------------------------------------" << endl;
		startTimer(&start_is);
	}
	for (int i=0; i<roots.size(); i++){
		for (int j=1; j<mult.at(i)+1; j++){
			// timer for computing "minors", qq(i)
	
			if (PROFILING){
				startTimer(&start);
			}

			cLinPoly.setCoefficient(0,-roots.at(i)); // cLinPoly = x - r(i)
			cPolyTemp = cNewDen; 
			// loop to make use of optimized division by linear polynomials (not clear this helps)
			//cPolyTemp /= cLinPoly^mult.at(i);
			for (int k=0; k<mult.at(i); k++)
				cPolyTemp /= cLinPoly;
	
			if (PROFILING){	
				stopTimerAddElapsed(&start,&elapsed_minor);	
				startTimer(&start);		
			}
			
			cTemp = computeResidue(roots.at(i), mult.at(i), j, cNum, cPolyTemp);
	
			if (PROFILING){	
				stopTimerAddElapsed(&start,&elapsed_res);
				startTimer(&start);
			}
			
			if (cTemp.realPart() != 0 && epsilonEqual(cTemp.realPart(),zero,prec))
				cTemp.setRealPart(zero);
			if (cTemp.imaginaryPart() != 0 && epsilonEqual(cTemp.imaginaryPart(),zero,prec))
				cTemp.setImaginaryPart(zero);
			if (j == 1)
				if (roots.at(i).imaginaryPart() == 0) {
					lg->push_back((RN)cTemp.realPart());
					Lg->push_back(CRNPolyRealPart(cLinPoly));
				}
				else {
					if (cTemp.imaginaryPart() != 0){
						RN a;
						RN b;
						CRN c;
						a = roots.at(i).realPart();
						b = roots.at(i).imaginaryPart();
						// consider implementing with RN polys
						c.set(-a,0);
						cPolyTemp.zero();
						cPolyTemp.setCoefficient(0,c);
						c.one();
						cPolyTemp.setCoefficient(1,c);
						c.set(-b,0);
						cPolyTemp /= c;
						rnTemp = RN(cTemp.imaginaryPart()*2);
						atn->push_back(rnTemp);
						suprnTemp = CRNPolyRealPart(cPolyTemp);
						Atn->push_back(suprnTemp);
					}
					if (cTemp.realPart() != 0){
						cPolyTemp = cLinPoly;
						cLinPoly.setCoefficient(0,-roots.at(i).conjugate());
						cPolyTemp *= cLinPoly;
						rnTemp = RN(cTemp.realPart());
						lg->push_back(rnTemp);
						suprnTemp = CRNPolyRealPart(cPolyTemp);
						Lg->push_back(suprnTemp);
					}
				}
			else {
				cPolyTemp = cLinPoly^(j-1);
				cPolyTemp *= 1-j;
				if (roots.at(i).imaginaryPart() == 0) {
					// rational term, real root
					//SparseUnivariatePolynomial<RN> supRN;
					suprnTemp.zero();
					suprnTemp.setCoefficient(0,cTemp.realPart());
					urfTemp.zero();
					urfTemp.setNumerator(suprnTemp);
					suprnTemp = CRNPolyRealPart(cPolyTemp);
					urfTemp.setDenominator(suprnTemp);
					g->push_back(urfTemp);
					//G.push_back(supRN);
					//G.push_back(CRNPolyRealPart(cPolyTemp));
				}
				else {
					// rational term, complex root
					SparseUnivariatePolynomial<CRN> cPolyTemp2;
					cPolyTemp2 = simplifyComplexRFTerm(cTemp,&cPolyTemp);
					urfTemp.zero();
					suprnTemp = CRNPolyRealPart(cPolyTemp2);
					urfTemp.setNumerator(suprnTemp);
					suprnTemp = CRNPolyRealPart(cPolyTemp);
					urfTemp.setDenominator(suprnTemp);
					g->push_back(urfTemp);
					//G.push_back(CRNPolyRealPart(cPolyTemp2));
					//G.push_back(CRNPolyRealPart(cPolyTemp));
				}
			}
	
			if (PROFILING){	
				stopTimerAddElapsed(&start,&elapsed_pp);
			}
		}
	}
	
	if (PROFILING){
		stopTimer(&start_is,&elapsed);
		cout << "\t\tminors (qq)\t\t" << elapsed_minor << " s" << endl;
		cout << "\t\tresidues\t\t" << elapsed_res << " s" << endl;
		cout << "\t\tpost-processing\t\t" << elapsed_pp << " s" << endl;
		cout << "--------------------------------------" << endl;
		cout << "\t\t\t\t\t" << elapsed << " s" << endl;
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


