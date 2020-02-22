#include "../../include/global.h"
#include <bpas.h>
#include "examples.h"

using namespace std;

extern float testTaylorShift(DUQP&, int);
extern float testUnivariate(DUQP&, lfixq, int);
extern float testMultivariate(std::vector<RationalRegularChain>, int, lfixq, int);
extern void testSLP();

int main(int argc, char *argv[]) {
	int test = 1;	// 0: Test Taylor Shift
			// 1: Test Univariate Real Root
			// 2: test Multivariate Real Root
	int choice = 1;	// Univariate or Multivariate
			// example number
	int n = 128;	// Univariate polynomial degree
	int check = 0;	// whether check the result
	if (argc > 1) {
		test = atoi(argv[1]);

		if (argc > 2 && atoi(argv[2]) > 0)
			choice = atoi(argv[2]);

		if (argc > 3 && atoi(argv[3]) > 0) {
			if (test < 2)
				n = atoi(argv[3]);
			else
				check = atoi(argv[3]);
		}
		if (test < 2 && argc > 4)
			check = atoi(argv[4]);
	}

	lfixq width(1, 20);
	if (!test) {			// Univariate Taylor Shift
		DenseUnivariateRationalPolynomial p(n);
		string filepath = "../../../examples/univariate/";
		assignPoly(&p, n, choice, filepath);	
		float elapsed = testTaylorShift(p, check);
	}
	else if (test == 1) {		// Univariate Real Root Isolation
		DUQP p(n);
		string filepath = "../../../examples/univariate/";
		assignPoly(&p, n, choice, filepath);
		float elapsed = testUnivariate(p, width, check);
	}
	else if (test == 2) {		// Multivariate Real Root Isolation
		string filepath = "../../../examples/multivariate/";
	        ostringstream convert;
	        convert << choice;
	        string filename = filepath + "multivariate_" + convert.str() + ".dat";
	        std::ifstream inputfile(filename.c_str(), std::ifstream::in);

	        int m;
	        inputfile >> m;
	        int vars;
	        inputfile >> vars;
		inputfile.close();

		char base = 'z';
		Symbol* xs = new Symbol[vars];
        	for (int i = 0; i < vars; ++i)
	                xs[i] = Symbol(string(sizeof(char), base-i));

		std::vector<RationalRegularChain> rcs (m, RationalRegularChain(vars, xs));
		assignRCs(rcs, choice, filepath, xs);
		float elapsed = testMultivariate(rcs, choice, width, check);
		delete [] xs;
	}
	else {
		testSLP();
	}

	return 0;
}

float testTaylorShift(DUQP& p, int check) {
        DUQP p2(p);
        unsigned long long ts_start;
        startTimer(&ts_start);
        p2.taylorShift();
        float elapsed;
        stopTimer(&ts_start, &elapsed);
        cout << "Taylor Shift\t" << elapsed << endl;
        if (check) {
                DUQP p3(p);
                unsigned long long o_start;
                startTimer(&o_start);
                p3.taylorShift(0);
                float elapsed2;
                stopTimer(&o_start, &elapsed2);
                cout << "CMY\t" << elapsed2 << endl;
                if (p2 == p3)
                        cout << "Pass." << endl;
                else {
                        cout << "Error." << endl;
                        exit(1);
                }
        }
        return elapsed;
}

float testUnivariate(DUQP& p, lfixq width, int check) {
        unsigned long long uni_start;
        startTimer(&uni_start);
        Intervals pIs = p.realRootIsolate(width);
        float elapsed;
        stopTimer(&uni_start, &elapsed);
        cout << "Univariate Real Root\t" << elapsed << endl;
        if (check) {
                unsigned long long o_start;
                startTimer(&o_start);
                float elapsed2;
                stopTimer(&o_start, &elapsed2);
                Intervals opIs = p.realRootIsolate(width, 0);
                cout << "CMY\t" << elapsed2 << endl;
                if (pIs.isExactEqual(opIs))
                        cout << "Pass." << endl;
                else {
                        cout << "Error." << endl;
                        exit(1);
                }
        }
        return elapsed;
}

float testMultivariate(std::vector<RationalRegularChain> rcs, int choice, lfixq width, int check) {
        unsigned long long multi_start;
        startTimer(&multi_start);
        Intervals mpIs = RealRootIsolation(rcs, width);
        float elapsed;
        stopTimer(&multi_start, &elapsed);
//cout << mpIs << endl;
        cout << "Multivariate Real Root\t" << elapsed << "\t";
        if (check) {
                ostringstream convert;
                convert << choice;
                string filename = "../../../examples/multivariate/roots_" + convert.str() + ".dat";
                ifstream inputfile(filename.c_str(), ifstream::in);
                int m;
                inputfile >> m;
                inputfile >> m;
                if (m != mpIs.numberOfIntervals()) {
                        cout << "Error: the number of roots doesn't match!" << endl;
                        exit(1);
                }
                else {
                        cout << "Pass." << endl;
                }
        }
        return elapsed;
}


/************************ Internal Test ********************************/
void testSLP() {
        ifstream inputfile("../../../examples/input.dat", ifstream::in);
        int n;
        inputfile >> n;
        int var;
        inputfile >> var;
	SMQP mp(var);

        for (int i = 0; i < n; ++i) {
		RationalNumberTerm t;
                inputfile >> t.coef;
		t.v = var;
		t.degs = new int[var];
                for (int j = 0; j < var; ++j)
                        inputfile >> t.degs[var-1-j];
		mp.setCoefficient(var, t.degs, t.coef);
        }
        inputfile.close();

        mp.straightLineProgram();
        mp.printSLP();
}
