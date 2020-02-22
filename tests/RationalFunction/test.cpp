#include <bpas.h>
#include <complex>
#include <mps/mps.h>
#include <mps/mpc.h>
//#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"

using namespace std;

extern void testRationalFunction();

int main(int argc, char *argv[]) {
	testRationalFunction();
	
	return 0;
}

void testRationalFunction() {
	cout << "-- RationalFunction --" << endl;
	Symbol variable("x");
	DUQP num(2);
	num.one();
	DUQP den(3);
	den.one();
	den.setCoefficient(1,2);
	den.setCoefficient(2,1);
	num.setVariableName(variable);
	den.setVariableName(variable);
	UnivariateRationalFunction<DUQP,RationalNumber> r(num,den);
	UnivariateRationalFunction<DUQP,RationalNumber> r2(r);
	DUQP dTemp,dTemp2;
	dTemp.setVariableName(variable);
	dTemp2.setVariableName(variable);
	bool testPass = true;
	r2 *= r;
	dTemp = den*den;
	dTemp2.one();
	if (!(r2.denominator() == dTemp) || !(r2.numerator() == dTemp2))
		testPass = false;
	r2 = r;
	r2 /= r;
	dTemp.one();
	if (!(r2.denominator() == dTemp) || !(r2.numerator() == dTemp2))
		testPass = false;
	r2 = r;
	r2 += r;
	dTemp = den;
	dTemp2 *= 2;
	if (!(r2.denominator() == dTemp) || !(r2.numerator() == dTemp2))
		testPass = false;
	r2 = r;
	r2 -= r;
	dTemp.one();
	dTemp2.zero();
	if (!(r2.denominator() == dTemp) || !(r2.numerator() == dTemp2))
		testPass = false;
	if (testPass)
		cout << "UnivariateRationalFunction<DUQP> arithmetic:\tpass" << endl;
	else
		cout << "UnivariateRationalFunction<DUQP> arithmetic:\tFAIL" << endl;
	
	SparseUnivariatePolynomial<RN> num2;
	SparseUnivariatePolynomial<RN> den2;
	num2.one();
	den2.one();
	den2.setCoefficient(1,RN(mpq_class(2)));
	den2.setCoefficient(2,RN(mpq_class(1)));
	num2.setVariableName(variable);
	den2.setVariableName(variable);
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> s(num2,den2);
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> s2(s);
	SparseUnivariatePolynomial<RN> sTemp,sTemp2;
	sTemp.setVariableName(variable);
	sTemp2.setVariableName(variable);
	testPass = true;
	s2 *= s;
	sTemp = den2*den2;
	sTemp2.one();
	if (!(s2.denominator() == sTemp) || !(s2.numerator() == sTemp2))
		testPass = false;
	s2 = s;
	s2 /= s;
	sTemp.one();
	if (!(s2.denominator() == sTemp) || !(s2.numerator() == sTemp2))
		testPass = false;
	s2 = s;
	s2 += s;
	sTemp = den2;
	sTemp2 *= 2;
	if (!(s2.denominator() == sTemp) || !(s2.numerator() == sTemp2))
		testPass = false;
	s2 = s;
	s2 -= s;
	sTemp.one();
	sTemp2.zero();
	if (!(s2.denominator() == sTemp) || !(s2.numerator() == sTemp2))
		testPass = false;
	if (testPass)
		cout << "UnivariateRationalFunction<SUP<RN>> arithmetic:\tpass" << endl;
	else
		cout << "UnivariateRationalFunction<SUP<RN>> arithmetic:\tFAIL" << endl;
	
	
	UnivariateRationalFunction<DUQP,RationalNumber> t(r);
	std::vector< UnivariateRationalFunction<DUQP,RationalNumber> > g;
	num.zero();
	num = (num + mpq_class(1) << 7);
	num.setCoefficient(4,-24);
	num.setCoefficient(2,-4);
	num.setCoefficient(1,8);
	num.setCoefficient(0,-8);
	den.zero();
	den = (den + mpq_class(1) << 8);
	den.setCoefficient(6,6);
	den.setCoefficient(4,12);
	den.setCoefficient(2,8);
	t.setNumerator(num);
	t.setDenominator(den);
	UnivariateRationalFunction<DUQP,RationalNumber> h;
	t.hermiteReduce(&g,&h);
	testPass = true;
	if (!(g.size() == 2))
		testPass = false;
	else {
		dTemp.one();
		dTemp2.zero();
		dTemp2 = dTemp2 + mpq_class(1) << 1;
		if (!(h.denominator() == dTemp2) || !(h.numerator() == dTemp))
			testPass = false;
		dTemp.zero();
		dTemp = dTemp + mpq_class(8) << 2;
		dTemp += mpq_class(4);
		dTemp2.zero();
		dTemp2 = dTemp2 + mpq_class(1) << 5;
		dTemp2.setCoefficient(3,4);
		dTemp2.setCoefficient(1,4);
		if (!(g.at(0).denominator() == dTemp2) || !(g.at(0).numerator() == dTemp))
			testPass = false;
		dTemp.zero();
		dTemp += mpq_class(3);
		dTemp2.zero();
		dTemp2 = dTemp2 + mpq_class(1) << 2;
		dTemp2 += mpq_class(2);
		if (!(g.at(1).denominator() == dTemp2) || !(g.at(1).numerator() == dTemp))
			testPass = false;
	}
	if (testPass)
		cout << "UnivariateRationalFunction<DUQP> hermiteReduce:\tpass" << endl;
	else
		cout << "UnivariateRationalFunction<DUQP> hermiteReduce:\tFAIL" << endl;
	
	
	std::vector< SparseUnivariatePolynomial<DUQP> > S;
	std::vector<DUQP> U;
	num.zero();
	num = num + mpq_class(1) << 4;
	num.setCoefficient(2,-3);
	num.setCoefficient(0,6);
	den.zero();
	den = den + mpq_class(1) << 6;
	den.setCoefficient(4,-5);
	den.setCoefficient(2,5);
	den.setCoefficient(0,4);
	r.set(num,den);
	cout << endl << "DUQP integration tests:" << endl;
	cout << "\tr = " << r << endl;
	r.setProfiling(false);
	r.integrateRationalFunctionLogPart(&S,&U);
	
	r.integrate(&dTemp,&g,&U,&S);
	cout << "r.integrate:" << endl;
	r.printIntegral(dTemp,g,U,S);
	
	vector<RationalNumber> lg,atn;
	vector<DUQP> Lg,Atn,Atn2;
		
	r.realSymbolicNumericIntegrate(&dTemp,&g,&lg,&Lg,&atn,&Atn,53);
	cout << "r.symbolicNumericIntegrate:" << endl;
	r.printIntegral(dTemp,g,lg,Lg,atn,Atn);
		
	den2.zero();
	den2 = (den2 + RN(mpq_class(1)) << 8);
	den2.setCoefficient(6,RN(mpq_class(6)));
	den2.setCoefficient(4,RN(mpq_class(12)));
	den2.setCoefficient(2,RN(mpq_class(8)));
	std::vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g2;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> h2;
	
	SparseUnivariatePolynomial<RN> denom,deriv;
	denom = den2;
	
	std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<RN> > > S2;
	std::vector< SparseUnivariatePolynomial<RN> > U2;
	num2.zero();
	num2 = num2 + RN(mpq_class(1)) << 4;
	num2.setCoefficient(2,RN(mpq_class(-3)));
	num2.setCoefficient(0,RN(mpq_class(6)));
	den2.zero();
	den2 = den2 + RN(mpq_class(1)) << 6;
	den2.setCoefficient(4,RN(mpq_class(-5)));
	den2.setCoefficient(2,RN(mpq_class(5)));
	den2.setCoefficient(0,RN(mpq_class(4)));
	s.set(num2,den2);
	cout << endl << "SUP<RN> integration tests:" << endl;
	cout << "\ts = " << s << endl;
	s.setProfiling(false);
	
	s.integrate(&sTemp,&g2,&U2,&S2);
	s.printIntegral(sTemp,g2,U2,S2);
	
	vector<RN> lg_2,atn_2;
	vector< SparseUnivariatePolynomial<RN> > Lg_2,Atn_2,Atn2_2;
		
	s.realSymbolicNumericIntegrate(&sTemp,&g2,&lg_2,&Lg_2,&atn_2,&Atn_2,53);
	s.setFloatingPointPrinting(true);
	s.setMatlabOutput();
	s.printIntegral(sTemp,g2,lg_2,Lg_2,atn_2,Atn_2);
		
	deriv = denom;
	deriv.differentiate(1);
//	cout << "A = " << denom << ", A' = " << deriv << ", gcd(A,A') = " << denom.gcd(deriv) << endl;
	
	cout << "SparseUnivariateDoublePolynomial tests:" << endl;
	SparseUnivariateDoublePolynomial<double> np(denom);
	SparseUnivariateDoublePolynomial<double> np2(den);
	//cout << "denom = " << denom << endl;
	//cout << "den = " << den << endl;
	denom.zero();
	denom.setCoefficient(0,RationalNumber(mpq_class(107)));
	S2.at(0) /= denom;
	SparseUnivariateDoublePolynomial<double> np3(S2.at(0),1.0);
	//cout << "S2.at(0) = " << S2.at(0) << endl;
	testPass = true;
	testPass = true;
	//cout << "np.degree() = " << np.degree() << ", np.leadingCoefficient() = " << np.leadingCoefficient() << endl;
	if (np.evaluate(0) != 0)
		testPass = false;
	if (np.evaluate(1) != 27)
		testPass = false;
	double test,value;
	value = 1.75;
	test = 8*pow(value,2);
	test += 12*pow(value,4);
	test += 6*pow(value,6);
	test += pow(value,8);
	if (np.evaluate(value) != test)
		testPass = false;
	if (np2.evaluate(0) != 4)
		testPass = false;
	if (np2.evaluate(1) != 5)
		testPass = false;
	value = 1.75;
	test = 4;
	test += 5*pow(value,2);
	test -= 5*pow(value,4);
	test += pow(value,6);
	if (np2.evaluate(value) != test)
		testPass = false;
	if (np3.evaluate(0) != -2)
		testPass = false;
	if (np3.evaluate(1) != 3)
		testPass = false;
	value = 1.75;
	test = -2;
	test += 6*pow(value,1);
	test += pow(value,2);
	test -= 2*pow(value,3);
	if (np3.evaluate(value) != test)
		testPass = false;
	if (testPass == false)
		cout << "SparseUnivariateDoublePolynomial<double> evaluate(double) test: fail" << endl;
	else
		cout << "SparseUnivariateDoublePolynomial<double> evaluate(double) test: pass" << endl;
	
	testPass = true;
	complex<double> ctest,cvalue;
	cvalue = complex<double>(1.0,1.5);
	ctest = complex<double>(1.0,0.0);
	ctest *= cvalue*cvalue;
	ctest += complex<double>(6.0,0.0);
	ctest *= cvalue*cvalue;
	ctest += complex<double>(12.0,0.0);
	ctest *= cvalue*cvalue;
	ctest += complex<double>(8.0,0.0);
	ctest *= cvalue*cvalue;
	//cout << "cvalue = " << cvalue << ", ctest = " << ctest << ", np.evaluate(cvalue) = " << np.evaluate(cvalue) << endl;
	if (np.evaluate(cvalue) != ctest)
		testPass = false;
	ctest = complex<double>(1.0,0.0);
	ctest *= cvalue*cvalue;
	ctest += complex<double>(-5.0,0.0);
	ctest *= cvalue*cvalue;
	ctest += complex<double>(5.0,0.0);
	ctest *= cvalue*cvalue;
	ctest += complex<double>(4.0,0.0);
	//cout << "cvalue = " << cvalue << ", ctest = " << ctest << ", np2.evaluate(cvalue) = " << np2.evaluate(cvalue) << endl;
	if (np2.evaluate(cvalue) != ctest)
		testPass = false;
	ctest = complex<double>(-2.0,0.0);
	ctest *= cvalue;
	ctest += complex<double>(1.0,0.0);
	ctest *= cvalue;
	ctest += complex<double>(6.0,0.0);
	ctest *= cvalue;
	ctest += complex<double>(-2.0,0.0);
	//cout << "cvalue = " << cvalue << ", ctest = " << ctest << ", np3.evaluate(cvalue) = " << np3.evaluate(cvalue) << endl;
	if (np3.evaluate(cvalue) != ctest)
		testPass = false;
	if (testPass == false)
		cout << "SparseUnivariateDoublePolynomial<double> evaluate(complex<double>) test: fail" << endl;
	else
		cout << "SparseUnivariateDoublePolynomial<double> evaluate(complex<double>) test: pass" << endl;
	
	testPass = true;
	SparseUnivariatePolynomial<CRN> cp;
	cp.setVariableName(Symbol("t"));
	CRN cTemp;
	cTemp = CRN(1, 1, 0, 1);
	cp.setCoefficient(1,cTemp);
	cTemp = CRN(0, 1, 1 ,1);

	cp.setCoefficient(3,cTemp);
	//cout << "cp = " << cp << endl;

	SparseUnivariateDoublePolynomial< complex<double> > np4(cp);
	//cout << "np4.evaluate(0) = " << np4.evaluate(0) << endl;
	//cout << "np4.evaluate(I) = " << np4.evaluate(complex<double>(0,1)) << endl;
	if (np4.evaluate(0.0) != complex<double>(0.0,0.0)) 
		testPass = false;
	if (np4.evaluate(complex<double>(0.0,1.0)) != complex<double>(1.0,1.0))
		testPass = false;
	cvalue = complex<double>(1.0,1.5);
	ctest = complex<double>(0.0,1.0);
	ctest *= cvalue;
	ctest *= cvalue;
	ctest += complex<double>(1.0,0.0);
	ctest *= cvalue;
	//cout << "np4.evaluate(" << cvalue << ") = " << np4.evaluate(cvalue) << ", ctest = " << ctest << endl;
	if (np4.evaluate(cvalue) != ctest)
		testPass = false;
	SparseUnivariatePolynomial< SparseUnivariatePolynomial<CRN> > scp;
	scp.setVariableName(Symbol("x"));
	scp.setCoefficient(0,cp);
	scp.setCoefficient(2,cp);
	//cout << "scp = " << scp << endl;
	SparseUnivariateDoublePolynomial< complex<double> > np5(scp,complex<double>(1,1));
	//cout << "np5.leadingCoefficient() = " << np5.leadingCoefficient() << endl;
	cvalue = complex<double>(1.0,1.5);
	ctest = np5.leadingCoefficient();
	ctest *= cvalue;
	ctest *= cvalue;
	ctest += np5.leadingCoefficient();
	//cout << "np5.evaluate(" << cvalue << ") = " << np5.evaluate(cvalue) << ", ctest = " << ctest << endl;
	if (np5.evaluate(cvalue) != ctest)
		testPass = false;
	if (testPass == false)
		cout << "SparseUnivariateDoublePolynomial< complex<double> > evaluate test: fail" << endl;
	else
		cout << "SparseUnivariateDoublePolynomial< complex<double> > evaluate test: pass" << endl;
	
	/*
    mpc_t x;
    mpf_t a,b;
    mpc_init2 (x, 256);
    mpc_set_d(x,-5.2,6.0);
    mpf_init_set(b,mpc_Im(x));
    cout << "x = " << x->r << "+" << x->i << "i" << endl;
    cout << "Im(x) = " << b << endl;
    mpf_init_set(a,mpc_Re(x));
    complexMPF cmpf;
    mpc_init2(cmpf.c,256);
    mpc_set(cmpf.c,x);
    SparseUnivariateMPComplexPolynomial mpcp(denom,53);
    cmpf = mpcp.evaluate(cmpf,53);
    cout << "denom(x) = (" << cmpf.c->r << "," << cmpf.c->i << ")" << endl;
    */
}


