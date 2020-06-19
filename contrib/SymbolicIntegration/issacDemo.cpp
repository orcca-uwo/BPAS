#include <bpas.h>
#ifndef SERIAL
	#include <cilktools/cilkview.h>
#endif

using namespace std;

// floatIntegrate.cpp
void floatIntegrate(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING);

SparseUnivariatePolynomial<RN> polyParse(string s);

extern void runExactExample(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runSNExampleDUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runSNExampleSUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);

int main(int argc, char *argv[]) {

	
	SparseUnivariatePolynomial<RN> A,D;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f;
	string variable = "x";
	A.setVariableName(variable);
	D.setVariableName(variable);
	
	int qsize, Asize, Dsize;
	
	f.setProfiling(false);				// default: false
	int prec = 53;						// default: 53 bit
	f.setFloatingPointPrinting(true);	// default: false
	int runCounter = 1;
	
	// Example: f = (2 + 2x - 2x^2 + 6x^3)/(-x + x^5) //
	A = polyParse("2 + 2*x - 2*x^2 + 6*x^3");
	D = polyParse("-x + x^5");
	
	f.set(A,D);
	
	runExactExample(f,prec,&runCounter);
	
	runSNExampleDUP(f,prec,&runCounter);
	
	// Example: f = (1)/(1 + x^4) //
	A = polyParse("1 + 0*x");
	D = polyParse("1 + x^4");
	
	f.set(A,D);
	
	f.setFloatingPointPrinting(false);	// default: false
	
	runSNExampleDUP(f,prec,&runCounter);
	
	f.setFloatingPointPrinting(true);	// default: false
	
	runSNExampleDUP(f,prec,&runCounter);
	
	
	// Example: f = (1)/(1+x+x^n), n < 255 //
	int n = 128;
	A.one();
	D.one();
	D.setCoefficient(1,RN(mpq_class(1)));
	D.setCoefficient(n,RN(mpq_class(1)));
	f.set(A,D);
	
	runSNExampleDUP(f,prec,&runCounter);
	
	return 0;
}




void startTimer(unsigned long long *start){
	#ifndef SERIAL
		*start = __cilkview_getticks();
	#endif
}

void endTimer(unsigned long long *start, unsigned long long *end, float *elapsed){
	#ifndef SERIAL
		*end = __cilkview_getticks();
		*elapsed = (*end - *start) / 1000.f;
	#endif
}

void printTime(float *elapsed, string name){
	#ifndef SERIAL
		cout << "--------------------------------------" << endl;
		cout << name << " runtime: " << *elapsed << " s" << endl;
		cout << "--------------------------------------" << endl;
	#endif
}

void runExactExample(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	DUQP P2;
	vector< UnivariateRationalFunction<DUQP,mpq_class> > g2;
	vector<DUQP> U;
	vector< SparseUnivariatePolynomial<DUQP> > S;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	unsigned long long end;
	float elapsed = 0;
	
	string snInt = "integrate";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	UnivariateRationalFunction<DUQP,mpq_class> j;
	j.setProfiling(PROFILING);
	j.setFloatingPointPrinting(f.isFloatingPointPrinting());
	DUQP num,den;
	num = f.numerator().convertToDUQP();
	den = f.denominator().convertToDUQP();
	j.set(num,den);
	j.integrate(&P2,&g2,&U,&S);
	
	endTimer(&start,&end,&elapsed);
	j.printIntegral(P2,g2,U,S);
	printTime(&elapsed,snInt);
	
	*counter += 1;
}

void runSNExampleDUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn;
	DUQP P2;
	vector< UnivariateRationalFunction<DUQP,mpq_class> > g2;
	vector<mpq_class> lg2,atn2;
	vector<DUQP> Lg2,Atn2;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	unsigned long long end;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrate";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	UnivariateRationalFunction<DUQP,mpq_class> j;
	j.setProfiling(PROFILING);
	j.setFloatingPointPrinting(f.isFloatingPointPrinting());
	DUQP num,den;
	num = f.numerator().convertToDUQP();
	den = f.denominator().convertToDUQP();
	j.set(num,den);
	j.realSymbolicNumericIntegrate(&P2,&g2,&lg2,&Lg2,&atn2,&Atn2,prec);
	
	endTimer(&start,&end,&elapsed);
	j.printIntegral(P2,g2,lg2,Lg2,atn2,Atn2);
	printTime(&elapsed,snInt);
	
	/*startTimer(&start);
	
	floatIntegrate(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
	
	endTimer(&start,&end,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,fInt);*/
	*counter += 1;
}

void runSNExampleSUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	unsigned long long end;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrate";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	
	endTimer(&start,&end,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,snInt);
	
	/*startTimer(&start);
	
	floatIntegrate(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
	
	endTimer(&start,&end,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,fInt);*/
	*counter += 1;
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


