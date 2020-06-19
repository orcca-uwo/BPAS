/* Comparison of two methods of symbolic-numeric integration of rational functions */
/*
 * This test code compares the built-in symbolic-numeric integration code of bpas
 * with brute-force partial fraction decomposition (BF-PFD) method suggested by
 * Fateman (2008).
 *
 * The built-in bpas integration is a member function of the template class
 * UnivariateRationalFunction<UnivariatePolynomialOverField,Field>,
 * where for the present context UnivariatePolynomialOverField (UPoF) and Field
 * can be either:
 *		1) UPoF = DUQP (DenseUnivariateRationalPolynomial) and
 *		   Field = RN (where RN = RationalNumber, the bpas rationals); or
 *		2) UPoF = SparseUnivariatePolynomial<RN> and
 *		   Field = RN
 * The integration algorithm is accessed through 
 * symbolicNumericIntegrate(conts,prec), where conts is a list of containers for
 * the solution (a polynomial part, rational part, log functions, arctan functions).
 * The solution can be printed by passing the containers to the printIntegral
 * member function.
 *
 * The BF-PFD method is provided in the floatIntegrate.cpp file in this directory.
 * it is accessed as the void function floatIntegrate(conts,prec), with the same
 * signature as the symbolicNumericIntegrate function.
 * 
 * This file contains two suites of test examples, accessed from "main" through 
 * the functions testExamples() and threeProblems(). In threeProblems() there are
 * three general problems that contain an integer n as a parameter: 
 * f(x) = 1/(1+x^n); f(x) = 1/(1+x+x^n); and f(x) = (x^2+2x+2)^n. Each are solved
 * by both symbolicNumericIntegrate and floatIntegrate. In testExamples() are a
 * series of textbook integration problems also solved by both methods.
 *
 * Since the integration code is polymorphic, it can be run with rational functions
 * using dense and sparsely represented polynomials, the two cases identified
 * above. Examples using the dense representation are run using the function 
 * "runExampleDUP", where "DUP" stands for "DenseUnivariatePolynomial", and 
 * examples using the sparse representation are run using the function 
 *  "runExampleSUP".
 *
 * Note that floating point output printing and whether profiling output is
 * printed are both controlled by private boolean member variables inside the
 * UnivariateRationalFunction class. These can be swiched using the member
 * functions "setFloatingPointPrinting" and "setProfiling".
 *
 * Note also that the timing of the code and the profiling of the code require the
 * the parallel version of the library, which requires an implementation of
 * CilkPlus. Under the serial version, "make serial" will make the serial version
 * of the test code, which will compute the integrals and print the results but
 * will not provide the comparison of run times for the integration methods.
 *
 */

#include <bpas.h>
#include <fstream>

using namespace std;

// floatIntegrate.cpp
void floatIntegrate(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING);

// realSymbolicNumericIntegratePFD.cpp
void realSymbolicNumericIntegratePFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> &f, SparseUnivariatePolynomial<RN> *P, vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > *g, vector<RN> *lg, vector< SparseUnivariatePolynomial<RN> > *Lg, vector<RN> *atn, vector< SparseUnivariatePolynomial<RN> > *Atn, int prec, bool PROFILING);

SparseUnivariatePolynomial<RN> polyParse(string s);


extern void testExamples();
extern void threeProblems();
extern void fibonacciTest();
extern void fall2017Talk();
extern void performanceTests();
extern void runExampleDUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runExampleSUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runExampleLRT(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runExamplePFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runExampleSimplePFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter);
extern void runExampleSUPToFile(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter, string &fname, int i);

int main(int argc, char *argv[]) {

//	fibonacciTest();
	threeProblems();
//	testExamples();
//	fall2017Talk();
//	performanceTests();
		
}

void performanceTests(){

	SparseUnivariatePolynomial<RN> A,D,T;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f;
	Symbol variable = 'x';
	A.setVariableName(variable);
	D.setVariableName(variable);
	T.setVariableName(variable);
	RN r;
	
	int qsize, Asize, Dsize;
	
	f.setProfiling(true);				// default: false
	int prec = 53;						// default: 53 bit
	f.setFloatingPointPrinting(true);	// default: false
//	f.setAnalyzingError(false);			// default: false
	int runCounter = 1;
	f.setMapleOutput();
	
	int n[13] = {12, 15, 19, 24, 31, 39, 50, 63, 81, 102, 131, 165, 212};
	remove("perftiming.txt");
	fstream fs;
	fs.open("perftiming.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	
	int iterations(0);
	int repetitions(0);
	
	f.setAnalyzingError(false);
	fs << "Example1NoEA" << std::endl; // f = (1)/(1+x^n), n < 255:	
	
	// Example 1: f = (1)/(1+x^n), n < 255 //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	iterations = 13;
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= -2;
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExampleLRT(f,prec,&runCounter);
		}
	}
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= -2;
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExamplePFD(f,prec,&runCounter);
		}
	}
		
	f.setAnalyzingError(true);
	fs << "Example1EA" << std::endl; // f = (1)/(1+x^n), n < 255:
	
	// Example 1: f = (1)/(1+x^n), n < 255 //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	iterations = 13;
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= -2;
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExampleLRT(f,prec,&runCounter);
		}
	}
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= -2;
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExamplePFD(f,prec,&runCounter);
		}
	}
	
	f.setAnalyzingError(false);	
	fs << "Example2NoEA" << std::endl; // f = (1)/(1+x+x^n), n < 255
	
	// Example 2: f = (1)/(1+x+x^n), n < 255 //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	iterations = 9;
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= 2;
			D.setCoefficient(1,RN(1));
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExampleLRT(f,prec,&runCounter);
		}
	}
	int num = 13;
	int previt = iterations;
	iterations = num - previt;
	repetitions = 1;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= 2;
			D.setCoefficient(1,RN(1));
			D.setCoefficient(n[previt+i-1],RN(1));
			f.set(A,D);
			runExampleLRT(f,prec,&runCounter);
		}
	}
	iterations = num;
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= 2;
			D.setCoefficient(1,RN(1));
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExamplePFD(f,prec,&runCounter);
		}
	}
	
	f.setAnalyzingError(true);	
	fs << "Example2EA" << std::endl; // f = (1)/(1+x+x^n), n < 255
	
	// Example 2: f = (1)/(1+x+x^n), n < 255 //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	iterations = 9;
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= 2;
			D.setCoefficient(1,RN(1));
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExampleLRT(f,prec,&runCounter);
		}
	}
	num = 13;
	previt = iterations;
	iterations = num - previt;
	repetitions = 1;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= 2;
			D.setCoefficient(1,RN(1));
			D.setCoefficient(n[previt+i-1],RN(1));
			f.set(A,D);
			runExampleLRT(f,prec,&runCounter);
		}
	}
	iterations = num;
	repetitions = 10;
	fs << iterations << "\t" << repetitions << std::endl;
	for (int j=1;j<=repetitions;++j) {
		for (int i=1;i<=iterations;++i) {
			A.one();
			D.negativeOne();
			D *= 2;
			D.setCoefficient(1,RN(1));
			D.setCoefficient(n[i-1],RN(1));
			f.set(A,D);
			runExamplePFD(f,prec,&runCounter);
		}
	}
	
	f.setAnalyzingError(true);	
	// Example 3: f = [n,n]_{exp(x)/x}(x) //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	
	// n = 8 //
	A = polyParse("-x^8-64*x^7-2016*x^6-40320*x^5-554400*x^4-5322240*x^3-34594560*x^2-138378240*x-259459200");
	D = polyParse("8*x^8-504*x^7+15120*x^6-277200*x^5+3326400*x^4-25945920*x^3+121080960*x^2-259459200*x");
	// n = 13 //
	/*A = polyParse("x^13+169*x^12+14196*x^11+780780*x^10+31231200*x^9+955674720*x^8");
	r = mpz_class(22936193280);
	A.setCoefficient(7,r);
	r = mpz_class(435787672320);
	A.setCoefficient(6,r);
	r = mpz_class(6536815084800);
	A.setCoefficient(5,r);
	r = mpz_class(76262842656000);
	A.setCoefficient(4,r);
	r = mpz_class(671113015372800);
	A.setCoefficient(3,r);
	r = mpz_class(4209708914611200);
	A.setCoefficient(2,r);
	r = mpz_class(16838835658444800);
	A.setCoefficient(1,r);
	r = mpz_class(32382376266240000);
	A.setCoefficient(0,r);
	D = polyParse("13*x^13-2184*x^12+180180*x^11-9609600*x^10+367567200*x^9");
	r = mpz_class(-10585935360);
	D.setCoefficient(8,r);
	r = mpz_class(234654900480);
	D.setCoefficient(7,r);
	r = mpz_class(-4022655436800);
	D.setCoefficient(6,r);
	r = mpz_class(52797352608000);
	D.setCoefficient(5,r);
	r = mpz_class(-516240781056000);
	D.setCoefficient(4,r);
	r = mpz_class(3562061389286400);
	D.setCoefficient(3,r);
	r = mpz_class(-15543540607795200);
	D.setCoefficient(2,r);
	r = mpz_class(32382376266240000);
	D.setCoefficient(1,r);*/
	// n = 21 //
	/*A = polyParse("x^21+441*x^20+97020*x^19+14132580*x^18+1526318640*x^17");
	r = mpz_class(129737084400);
	A.setCoefficient(16,r);
	r = mpz_class(8995104518400);
	A.setCoefficient(15,r);
	r = mpz_class(520431047136000);
	A.setCoefficient(14,r);
	r = mpz_class(25501121309664000);
	A.setCoefficient(13,r);
	r = mpz_class(1068213637082592000);
	A.setCoefficient(12,r);
	r = mpz_class("38455690934973312000");
	A.setCoefficient(11,r);
	r = mpz_class("1192126418984172672000");
	A.setCoefficient(10,r);
	r = mpz_class("31790037839577937920000");
	A.setCoefficient(9,r);
	r = mpz_class("726280095258049812480000");
	A.setCoefficient(8,r);
	r = mpz_class("14110584707870682071040000");
	A.setCoefficient(7,r);
	r = mpz_class("230472883561887807160320000");
	A.setCoefficient(6,r);
	r = mpz_class("3111383928085485396664320000");
	A.setCoefficient(5,r);
	r = mpz_class("33859178040930282257817600000");
	A.setCoefficient(4,r);
	r = mpz_class("285921947901189050177126400000");
	A.setCoefficient(3,r);
	r = mpz_class("1760677258128374677406515200000");
	A.setCoefficient(2,r);
	r = mpz_class("7042709032513498709626060800000");
	A.setCoefficient(1,r);
	r = mpz_class("13750050968240640337841356800000");
	A.setCoefficient(0,r);
	D = polyParse("21*x^21-9240*x^20+2018940*x^19-290727360*x^18");
	r = mpz_class(30889782000);
	D.setCoefficient(17,r);
	r = mpz_class(-2570029862400);
	D.setCoefficient(16,r);
	r = mpz_class(173477015712000);
	D.setCoefficient(15,r);
	r = mpz_class(-9714712879872000);
	D.setCoefficient(14,r);
	r = mpz_class(457805844463968000);
	D.setCoefficient(13,r);
	r = mpz_class("-18312233778558720000");
	D.setCoefficient(12,r);
	r = mpz_class("624447171848852352000");
	D.setCoefficient(11,r);
	r = mpz_class("-18165735908330250240000");
	D.setCoefficient(10,r);
	r = mpz_class("449601963731173693440000");
	D.setCoefficient(9,r);
	r = mpz_class("-9407056471913788047360000");
	D.setCoefficient(8,r);
	r = mpz_class("164623488258491290828800000");
	D.setCoefficient(7,r);
	r = mpz_class("-2370578230922274587934720000");
	D.setCoefficient(6,r);
	r = mpz_class("27409810795038799922995200000");
	D.setCoefficient(5,r);
	r = mpz_class("-245075955343876328723251200000");
	D.setCoefficient(4,r);
	r = mpz_class("1592993709735196136701132800000");
	D.setCoefficient(3,r);
	r = mpz_class("-6707341935727141628215296000000");
	D.setCoefficient(2,r);
	r = mpz_class("13750050968240640337841356800000");
	D.setCoefficient(1,r);*/
	f.set(A,D);
	
//	runExampleSUP(f,prec,&runCounter);
	
	// Example 4: f = 2x/(x^2 - (1 + epsilon)^2) //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	r = RationalNumber(mpq_class(1e16));
	r = RationalNumber(1)/r;
	r += RationalNumber(1);
	r *= r;
	r *= -1;
	A.zero();
	D.zero();
	A.setCoefficient(1,RationalNumber(2));
	D.setCoefficient(2,RationalNumber(1));
	D.setCoefficient(0,r);
	f.set(A,D);
	
//	runExampleSUP(f,prec,&runCounter);
	
	// I think there is a issue with the error analysis that is revealed on this problem. The code is detecting the singularities on the real axis, but I'm not sure I beleieve the accuracy assessments. Also, the error value comes out to be judged to be zero, which cannot be right. I will test out the evaluations. There is also an issue with the FE evaluation, which comes out at 10^{182}.
	
	// Example 5: f = 2x/(x^2 + epsilon^2) //
	// In this comparison we compare BPAS LRT and PFD methods using the sparse representation.  //
	r = RationalNumber(mpq_class(1e16));
	r = RationalNumber(1)/r;
	r *= r;
	A.zero();
	D.zero();
	A.setCoefficient(1,RationalNumber(2));
	D.setCoefficient(2,RationalNumber(1));
	D.setCoefficient(0,r);
	f.set(A,D);
	
//	runExampleSUP(f,prec,&runCounter);
	
	// Example 6: f = (x^2-1)/(x^4+5*x^2+7) //
	A = polyParse("x^2-1");
	D = polyParse("x^4+5*x^2+7");
	f.set(A,D);
	
	//runExampleSUP(f,prec,&runCounter);
	
	fs.close();
}

void fibonacciTest(){

	SparseUnivariatePolynomial<RN> A,D;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f;
	Symbol variable = 'x';
	A.setVariableName(variable);
	D.setVariableName(variable);
	
	int qsize, Asize, Dsize;
	
	f.setProfiling(false);				// default: false
	int prec = 53;					// default: 53 bit
	f.setFloatingPointPrinting(true);		// default: false
	f.setAnalyzingError(false);			// default: false
	int runCounter = 1;
	f.setMapleOutput();
	
	
	// Example 1: f = (1)/(1+x^n), n < 255 //
	// In this comparison we compare BF-PFD with BPAS using the sparse representation,  //
	// which exhibits faster performance than the dense representation on this problem. //

	string fname = "data.txt";
	int i(8),i_prior(5);
/*
	ofstream myfile(fname.c_str(), ios::app);
	if(!myfile.is_open()){
		cout << "file not open " << endl;
		return;
	}
	myfile.close();
*/
	for(; i<145;){
		cout << "i = " << i << endl;
		A.one();
		D.negativeOne();
		D *= 2;
		D.setCoefficient(i,RN(1));
		f.set(A,D);
	
		//runExampleSUP(f,prec,&runCounter);
		runExampleSUPToFile(f,prec,&runCounter, fname, i);
		i = i + i_prior;
		i_prior = i - i_prior;
	}
}


void fall2017Talk(){

	SparseUnivariatePolynomial<RN> A,D;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f;
	Symbol variable = 'x';
	A.setVariableName(variable);
	D.setVariableName(variable);
	
	int qsize, Asize, Dsize;
	
	f.setProfiling(false);				// default: false
	int prec = 53;					// default: 53 bit
	f.setFloatingPointPrinting(true);		// default: false
	f.setAnalyzingError(false);			// default: false
	int runCounter = 1;
	f.setMapleOutput();
	char ch;
	
	
	// Example 1: f = (1)/(1+x^n), n < 255 //
	// In this comparison we compare BF-PFD with BPAS using the sparse representation,  //
	// which exhibits faster performance than the dense representation on this problem. //
	
	A = polyParse("5 - x + x^2 + 4*x^3 - x^4 + x^5");
	D = polyParse("4 - 4*x + 5*x^2 - 2*x^3 + x^4");
	
	f.set(A,D);
	
	runExampleLRT(f,prec,&runCounter);
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	A = polyParse("-96 + 576*x + 1440*x^2 - 312*x^3 - 1280*x^4 - 192*x^5 + 256*x^6 + 72*x^7");
	D = polyParse("9 - 108*x + 288*x^2 + 468*x^3 - 78*x^4 - 252*x^5 - 32*x^6 + 36*x^7 + 9*x^8");
	
	f.set(A,D);
	
	runExampleLRT(f,prec,&runCounter);
	
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	runExamplePFD(f,prec,&runCounter);
	
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);

	int n = 128;
	A.one();
	D.negativeOne();
	D *= 2;
	D.setCoefficient(n,RN(1));
	f.set(A,D);
	
	runExampleLRT(f,prec,&runCounter);
	
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	runExampleSimplePFD(f,prec,&runCounter);
	
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	n = 128;
	A.one();
	D.one();
	D.setCoefficient(1,RN(1));
	D.setCoefficient(n,RN(1));
	f.set(A,D);
	
	runExampleLRT(f,prec,&runCounter);
	
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	runExampleSimplePFD(f,prec,&runCounter);
	
	cout << "End of Algorithm Comparison...Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	f.setProfiling(true);				// default: false
	prec = 53;					// default: 53 bit
	f.setFloatingPointPrinting(true);		// default: false
	f.setAnalyzingError(true);			// default: false
	
	A = polyParse("-7 - 2*x + 23*x^2 + 30*x^3 - 6*x^4 - 26*x^5 - 4*x^6 - 12*x^7 + x^8 + 8*x^9");
	D = polyParse("-2 - 4*x + 3*x^2 + 10*x^3 + 7*x^4 - 4*x^6 - 2*x^7 - 2*x^8 + x^10");
	f.set(A,D);
	
	runExampleLRT(f,prec,&runCounter);
	
	cout << "Press any key to continue:" << std::endl;	
	cin.get(ch);
	
	prec = 128;					// default: 53 bit
	
	runExampleLRT(f,prec,&runCounter);
	
}


void threeProblems(){

	SparseUnivariatePolynomial<RN> A,D;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f;
	Symbol variable = 'x';
	A.setVariableName(variable);
	D.setVariableName(variable);
	
	int qsize, Asize, Dsize;
	
	f.setProfiling(true);				// default: false
	int prec = 53;						// default: 53 bit
	f.setFloatingPointPrinting(true);	// default: false
	f.setAnalyzingError(false);			// default: false
	int runCounter = 1;
	f.setMapleOutput();
	
	
	// Example 1: f = (1)/(2+x^n), n < 255 //
	// In this comparison we compare BF-PFD with BPAS using the sparse representation,  //
	// which exhibits faster performance than the dense representation on this problem. //

	int n = 128;
	A.one();
	D.one();
	D *= 2;
	D.setCoefficient(n,RN(1));
	f.set(A,D);
	
//	runExampleSUP(f,prec,&runCounter);
	
	
	// Example 2: f = (1)/(2+x+x^n), n < 255 //
	// In this comparison we compare BF-PFD with BPAS using the sparse representation,  //
	// which exhibits faster performance than the dense representation on this problem. //
	
	n = 128;
	A.one();
	D.one();
	D *= 2;
	D.setCoefficient(1,RN(1));
	D.setCoefficient(n,RN(1));
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
	
	
	// Example 3: f = (1)/(1+2*x+2*x^2)^n, n <= 32 //
	// In this comparison we compare BF-PFD with BPAS using the dense representation,  //
	// which exhibits slightly better performance than the sparse representation on this problem. //
	
	n = 4;
	A.one();
	D = polyParse("1+2*x+2*x^2");
	D = D^n;
	f.set(A,D);
	
	//runExampleDUP(f,prec,&runCounter);
}

void testExamples(){
	
	SparseUnivariatePolynomial<RN> A,D;
	UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f;
	
	f.setProfiling(true);				// default: false
	int prec = 53;					// default: 53 bit
	f.setFloatingPointPrinting(true);		// default: false
	f.setAnalyzingError(true);			// default: false
	int runCounter = 1;
	
	
	A = polyParse("1 + x - x^2 + 3*x^3");
	A = 4*A;
	D = polyParse("-2*x + 2*x^5");
	
	f.set(A,D);
	
	// need to change to DUP and fix error
	runExampleSUP(f,prec,&runCounter);
	
	
	A = polyParse("6 - 3*x^2 + x^4");
	D = polyParse("4 + 5*x^2 - 5*x^4 + x^6");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);

	
	A = polyParse("-8 + 8*x - 4*x^2 - 24*x^4 + x^7");
	D = polyParse("8*x^2 + 12*x^4 + 6*x^6 + x^8");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
	
	
	A = polyParse("36 + 0*x");
	A.setVariableName('x');
	D = polyParse("-2 + x + 4*x^2 - 2*x^3 - 2*x^4 + x^5");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
	
	
	A = polyParse("5 - x + x^2 + 4*x^3 - x^4 + x^5");
	D = polyParse("4 - 4*x + 5*x^2 - 2*x^3 + x^4");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);


	A = polyParse("-7 - 2*x + 23*x^2 + 30*x^3 - 6*x^4 - 26*x^5 - 4*x^6 - 12*x^7 + x^8 + 8*x^9");
	D = polyParse("-2 - 4*x + 3*x^2 + 10*x^3 + 7*x^4 - 4*x^6 - 2*x^7 - 2*x^8 + x^10");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
	
	
	A = polyParse("1 + 0*x");
	D = polyParse("1 + x^4");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
	
	
	A = polyParse("1 + 0*x");
	D = polyParse("1 + x^6");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
	
	
	A = polyParse("-96 + 576*x + 1440*x^2 - 312*x^3 - 1280*x^4 - 192*x^5 + 256*x^6 + 72*x^7");
	D = polyParse("9 - 108*x + 288*x^2 + 468*x^3 - 78*x^4 - 252*x^5 - 32*x^6 + 36*x^7 + 9*x^8");
	
	f.set(A,D);
	
	runExampleSUP(f,prec,&runCounter);
}


void printTime(float *elapsed, string name){
	//#ifndef SERIAL
	cout << "--------------------------------------" << endl;
	cout << name << " runtime: " << *elapsed << " s" << endl;
	cout << "--------------------------------------" << endl;
	//#endif
}

void runExampleDUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn,Atn2;
	DUQP P2;
	vector< UnivariateRationalFunction<DUQP,RN> > g2;
	vector<RN> lg2,atn2;
	vector<DUQP> Lg2,Atn2a,Atn2b;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrate";
	string snIntPFD = "realSymbolicNumericIntegratePFD";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	UnivariateRationalFunction<DUQP,RN> j;
	j.setProfiling(PROFILING);
	j.setFloatingPointPrinting(f.isFloatingPointPrinting());
	j.setAnalyzingError(f.isAnalyzingError());
	DUQP num,den;
	num = f.numerator().convertToDUQP();
	den = f.denominator().convertToDUQP();
	j.set(num,den);
	
	//f.realSymbolicNumericIntegratePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,&Atn2,prec);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	//f.printIntegral(P,g,lg,Lg,atn,Atn,Atn2);
	printTime(&elapsed,snInt);
	
	startTimer(&start);
	
	f.realSymbolicNumericIntegratePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,snIntPFD);
	
	/*startTimer(&start);
	
	floatIntegrate(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,fInt);*/
//	j.realSymbolicNumericIntegratePFD(&P2,&g2,&lg2,&Lg2,&atn2,&Atn2a,prec);
//	//j.realSymbolicNumericIntegrate(&P2,&g2,&lg2,&Lg2,&atn2,&Atn2a,&Atn2b,prec);
//	
//	stopTimer(&start,&elapsed);
//	j.printIntegral(P2,g2,lg2,Lg2,atn2,Atn2a,Atn2b);
//	printTime(&elapsed,snInt);
//	
//	/*startTimer(&start);
//	
//	realSymbolicNumericIntegratePFD(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
//	
//	stopTimer(&start,&elapsed);
//	f.printIntegral(P,g,lg,Lg,atn,Atn);
//	printTime(&elapsed,snIntPFD);
//	
//	startTimer(&start);
//	
//	floatIntegrate(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
//	
//	stopTimer(&start,&elapsed);
//	f.printIntegral(P,g,lg,Lg,atn,Atn);
//	printTime(&elapsed,fInt);*/
	*counter += 1;
}

void runExampleSUP(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn,Atn2;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrate";
	string snIntPFD = "realSymbolicNumericIntegratePFD";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	//f.realSymbolicNumericIntegratePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,&Atn2,prec);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	//f.printIntegral(P,g,lg,Lg,atn,Atn,Atn2);
	printTime(&elapsed,snInt);
	
	startTimer(&start);
	
	f.realSymbolicNumericIntegratePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,snIntPFD);
	
	/*startTimer(&start);
	
	floatIntegrate(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,fInt);*/
	*counter += 1;
}

void runExampleLRT(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn,Atn2;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrateLRT";
	string snIntPFD = "realSymbolicNumericIntegratePFD";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	//f.realSymbolicNumericIntegratePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,&Atn2,prec);
	
	stopTimer(&start,&elapsed);
	//f.printIntegral(P,g,lg,Lg,atn,Atn);
	f.printIntegral(P,g,lg,Lg,atn,Atn,Atn2);
	printTime(&elapsed,snInt);
	*counter += 1;
}

void runExamplePFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn,Atn2;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrateLRT";
	string snIntPFD = "realSymbolicNumericIntegratePFD";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	f.realSymbolicNumericIntegratePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	//f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,&Atn2,prec);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	//f.printIntegral(P,g,lg,Lg,atn,Atn,Atn2);
	printTime(&elapsed,snIntPFD);
	*counter += 1;
}

void runExampleSimplePFD(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter){
	
	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn,Atn2;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrateLRT";
	string snIntPFD = "realSymbolicNumericIntegrateSimplePFD";
	
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	startTimer(&start);
	
	f.realSymbolicNumericIntegrateSimplePFD(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	//f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,&Atn2,prec);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	//f.printIntegral(P,g,lg,Lg,atn,Atn,Atn2);
	printTime(&elapsed,snIntPFD);
	*counter += 1;
}


void runExampleSUPToFile(UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> f, int prec, int *counter, string &fname, int i){
	
	ofstream myfile(fname.c_str(), ios::app);
	if(!myfile.is_open()){
		cout << "file not open " << endl;
		return;
	}	

	SparseUnivariatePolynomial<RN> P;
	vector< UnivariateRationalFunction<SparseUnivariatePolynomial<RN>,RN> > g;
	vector<RN> lg,atn;
	vector< SparseUnivariatePolynomial<RN> > Lg,Atn;
	
	bool PROFILING = f.isProfiling();
	
	// Profiling variables
	unsigned long long start;
	float elapsed = 0;
	
	string fInt = "floatIntegrate";
	string snInt = "realSymbolicNumericIntegrate";
	string snIntPFD = "realSymbolicNumericIntegratePFD";
	/*
	cout << "--------------------------------------" << endl;
	cout << "Example " << *counter << ":" << endl;
	cout << "f = " << f << endl;
	cout << "--------------------------------------" << endl;
	*/
	startTimer(&start);
	
	f.realSymbolicNumericIntegrate(&P,&g,&lg,&Lg,&atn,&Atn,prec);
	
	stopTimer(&start,&elapsed);
	//f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,snInt);
	myfile << i << ","<< elapsed << endl;

	myfile.close();
	/*
	startTimer(&start);
	
	realSymbolicNumericIntegratePFD(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
	
	stopTimer(&start,&elapsed);
	f.printIntegral(P,g,lg,Lg,atn,Atn);
	printTime(&elapsed,snIntPFD);
	
	startTimer(&start);
	
	floatIntegrate(f,&P,&g,&lg,&Lg,&atn,&Atn,prec,PROFILING);
	
	stopTimer(&start,&elapsed);
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


