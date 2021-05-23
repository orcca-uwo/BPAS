#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <bpas.h>
#include <vector>
#include "sysGen.hpp"
#include "../MapleTestTool/MapleTestTool.hpp"

#if defined(MAPLE_VALIDATE) && MAPLE_VALIDATE
#include "../../include/MapleInterface/MapleInterfaceStream.hpp"
#endif

using namespace std;




void testSys(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars, string sysNumber, string sysName, bool showOutput, bool isLazard, bool medianTiming) {

    // Output system name
    fstream fs;
    int nchars = 8;
    std::string name;
    if (sysName.size() > nchars)
    		name = sysName.substr(0,nchars/2) + sysName.substr(sysName.size()-nchars/2,nchars/2);
    else
    		name = sysName;
    RegularChain<RN,SMQP> rc;

    if (!medianTiming) {
	    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	    fs << name << "\t";
	    fs.close();
	    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	    fs << name << "\t";
	    fs.close();
	    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	    fs << name << "\t";
	    fs.close();

		cerr << "calling " << sysNumber << " (" << sysName << "):" << endl;
		#if REGULARCHAIN_DEBUG
		if (showOutput) {
			std::cerr << "Input polynomials:" << std::endl;
			for (size_t i = 0; i < polys.size(); ++i)
				std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
			std::cerr << "PolynomialRing:" << std::endl;
			printVariables(vars,"R");
	//		exit(0);
		}
		#endif
		bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);
	//	bool pass(true);
		#if REGULARCHAIN_DEBUG
		if (showOutput) {
			std::cerr << "Input polynomials:" << std::endl;
			for (size_t i = 0; i < polys.size(); ++i)
				std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
			std::cerr << "PolynomialRing:" << std::endl;
			printVariables(vars,"R");
		}
		#endif
		if (pass) {
			std::cerr << "(" << sysName << ")" << std::endl;
			std::cerr << "System " <<  sysNumber << " Test:\t\t\t\t\t\t PASSED" << std::endl;
		}
		else {
			std::cerr << "System " <<  sysNumber << " Test:\t\t\t\t\t\t FAILED" << std::endl;
			exit(1);
		}
	}
	else {

		// Timing variables
		long long unsigned int start(0);
		float elapsed;
		int nRuns(9);
		std::vector<float> runTimes;

		fs.open("medianTiming.txt",std::fstream::in | std::fstream::out | std::fstream::app);
		fs << name << "\t";
		// Redirect stdout
		std::fstream f("/dev/null");

		std::cerr << "Testing " << name << std::endl;

		//not thread-safe
		streambuf *oldcout,*oldcerr;

		if (!showOutput) {
		   oldcout = cout.rdbuf(); // <-- save
		   oldcerr = cerr.rdbuf(); // <-- save

		   cout.rdbuf(f.rdbuf());
		   cerr.rdbuf(f.rdbuf());
		   // cout.rdbuf (ss.rdbuf());
		   // cerr.rdbuf (ss.rdbuf());
		}
		#if REGULARCHAIN_DEBUG
		for (size_t i = 0; i < polys.size(); ++i) {
		   cerr << "F[" << i << "] = " << polys[i] << endl;
		}
		#endif

		for (int i=0; i<nRuns; ++i) {
			startTimer(&start);

			rc.triangularize(polys,isLazard);

			stopTimer(&start,&elapsed);
			runTimes.push_back(elapsed);
		}

		std::sort(runTimes.begin(),runTimes.end());
		fs << runTimes[nRuns/2] << std::endl;
		fs.close();

		if (!showOutput) {
			cout.rdbuf (oldcout);   // <-- restore
			cerr.rdbuf (oldcerr);   // <-- restore
		   f.close();
		}

	}
}


void testSys(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& algvars, std::vector<Symbol>& transVars, string sysNumber, string sysName, bool showOutput, bool isLazard, bool medianTiming) {

    // Output system name
    fstream fs;
    int nchars = 8;
    std::string name;
    if (sysName.size() > nchars)
    		name = sysName.substr(0,nchars/2) + sysName.substr(sysName.size()-nchars/2,nchars/2);
    else
    		name = sysName;
    RegularChain<RN,SMQP> rc(algvars, transVars);

    if (!medianTiming) {
	    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	    fs << name << "\t";
	    fs.close();
	    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	    fs << name << "\t";
	    fs.close();
	    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	    fs << name << "\t";
	    fs.close();

		cerr << "calling " << sysNumber << " (" << sysName << "):" << endl;
		#if REGULARCHAIN_DEBUG
		if (showOutput) {
			std::cerr << "Input polynomials:" << std::endl;
			for (size_t i = 0; i < polys.size(); ++i)
				std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
			std::cerr << "PolynomialRing:" << std::endl;
			printVariables(algvars,"R");
			std::cerr << "Transcendentals: " << std::endl;
			printVariables(transVars, "ts");
	//		exit(0);
		}
		#endif
		bool pass = triangularizeValidate(polys,rc,algvars,transVars,showOutput,isLazard);
	//	bool pass(true);
		#if REGULARCHAIN_DEBUG
		if (showOutput) {
			std::cerr << "Input polynomials:" << std::endl;
			for (size_t i = 0; i < polys.size(); ++i)
				std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
			std::cerr << "PolynomialRing:" << std::endl;
			printVariables(algvars,"R");
		}
		#endif
		if (pass) {
			std::cerr << "(" << sysName << ")" << std::endl;
			std::cerr << "System " <<  sysNumber << " Test:\t\t\t\t\t\t PASSED" << std::endl;
		}
		else {
			std::cerr << "System " <<  sysNumber << " Test:\t\t\t\t\t\t FAILED" << std::endl;
			exit(1);
		}
	}
	else {

		// Timing variables
		long long unsigned int start(0);
		float elapsed;
		int nRuns(9);
		std::vector<float> runTimes;

		fs.open("medianTiming.txt",std::fstream::in | std::fstream::out | std::fstream::app);
		fs << name << "\t";
		// Redirect stdout
		std::fstream f("/dev/null");

		std::cerr << "Testing " << name << std::endl;

		//not thread-safe
		streambuf *oldcout,*oldcerr;

		if (!showOutput) {
		   oldcout = cout.rdbuf(); // <-- save
		   oldcerr = cerr.rdbuf(); // <-- save

		   cout.rdbuf(f.rdbuf());
		   cerr.rdbuf(f.rdbuf());
		   // cout.rdbuf (ss.rdbuf());
		   // cerr.rdbuf (ss.rdbuf());
		}
		#if REGULARCHAIN_DEBUG
		for (size_t i = 0; i < polys.size(); ++i) {
		   cerr << "F[" << i << "] = " << polys[i] << endl;
		}
		#endif

		for (int i=0; i<nRuns; ++i) {
			startTimer(&start);

			rc.triangularize(polys,isLazard);

			stopTimer(&start,&elapsed);
			runTimes.push_back(elapsed);
		}

		std::sort(runTimes.begin(),runTimes.end());
		fs << runTimes[nRuns/2] << std::endl;
		fs.close();

		if (!showOutput) {
			cout.rdbuf (oldcout);   // <-- restore
			cerr.rdbuf (oldcerr);   // <-- restore
		   f.close();
		}

	}
}
//bool testSubsetAsConstructibleSetsPappus(vector<RegularChain<RN,SMQP>> lrc, bool showOutput, bool isLazard) {
//
//    stringstream ss;
//    // Redirect stdout
//    streambuf *oldcout,*oldcerr;
//
//    if (!showOutput) {
//		oldcout = cout.rdbuf(); // <-- save
//		oldcerr = cerr.rdbuf(); // <-- save
//		cout.rdbuf (ss.rdbuf());
//		cerr.rdbuf (ss.rdbuf());
//    }
//
//	std::cerr << "BPAS result:" << std::endl;
//    for (int i=0; i<lrc.size(); ++i) {
//    	cerr << "lrc[" << i << "] = " << lrc[i] << endl;
//    }
//
//    if (!showOutput) {
//		cout.rdbuf (oldcout);   // <-- restore
//		cerr.rdbuf (oldcerr);   // <-- restore
//	}

//    ExpressionTree lrcTrees;
//    lrcTrees.fromVector<RegularChain<RN,SMQP>>(lrc);

//    std::vector<std::string> inputs;
//    inputs.push_back(lrcTrees.toMapleString() + ":");

//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    char* cstr;

//    cstr = new char[inputs[0].length()+1];
//    std::strcpy(cstr, inputs[0].c_str());
//    ALGEB res = EvalMapleStatement(kv, cstr);
//    delete[] cstr;

//	std::string procStr = "SubsetAsConstructibleSetsPappus := proc (dec1::list) local n, dec2, lrc1, lrc2, lrs1, lrs2, cs1, cs2, R, i, rc; R := RegularChains:-PolynomialRing([a, b, c, d, e, t, u, v, w, x, y, z]); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; dec2 := [[y*a+(w-y)*t, b*v+c*x-v*x, (t*x-v*z)*c+(-x+z)*v*t, d*v+e*y-v*y, (u*y-v*z)*e+(-y+z)*v*u, (w*x-x*y)*t+(-y*w+x*y)*u], [a, c, e, t, u, v], [y*a+(w-y)*t, b*t+c*z-t*z, e, u, v, x], [x*a+(w-x)*u, c, d*u+e*z-u*z, t, v, y], [b*t+c*z-t*z, d*u+e*z-u*z, v, w, x, y], [a, b*v+c*x-v*x, d*v+e*y-v*y, t, u, z], [y*a+(w-y)*t, b, d*v+e*y-v*y, u, x, z], [x*a+(w-x)*u, b*v+c*x-v*x, d, t, y, z], [t, u, v, x, y, z], [b, d, w, x, y, z], [y*a+(w-y)*t, b-z, c, d-z, e, (w*x-x*y)*t+(-y*w+x*y)*u, v], [a, c, d-z, e, t, v, w-x], [a, b*v+c*x-v*x, (t*x-v*z)*c+(-x+z)*v*t, d-y, e, u, w-y], [a, b-z, c, e, u, v, w-y], [a, b*v+c*y-v*y, (t*y-v*z)*c+(-y+z)*v*t, d*v+e*y-v*y, (u*y-v*z)*e+(-y+z)*v*u, w-y, x-y], [z*a+(w-z)*u, b*v+c*z-v*z, d*v+e*y-v*y, (u*y-v*z)*e+(-y+z)*v*u, t-v, (y*w-y*z)*u+(-z*w+y*z)*v, x-z], [y*a+(w-y)*t, b, c-t, d-y, e, u, x], [y*a+(w-y)*t, b-z, c, e, u, v, x], [a-t, b, c-t, d*v+e*y-v*y, (u*y-v*z)*e+(-y+z)*v*u, w, x], [a-t, b*t+c*z-t*z, d-z, e, v, w, x], [z*a+(w-z)*t, b*v+c*x-v*x, (t*x-v*z)*c+(-x+z)*v*t, d*v+e*z-v*z, (w*x-x*z)*t+(-z*w+x*z)*v, u-v, y-z], [z*a+(w-z)*v, b*v+c*z-v*z, d*v+e*z-v*z, t-v, u-v, x-z, y-z], [x*a+(w-x)*u, b-x, c, d, e-u, t, y], [x*a+(w-x)*u, c, d-z, e, t, v, y], [a-u, b*v+c*x-v*x, (t*x-v*z)*c+(-x+z)*v*t, d, e-u, w, y], [a-u, b-z, c, d*u+e*z-u*z, v, w, y], [a-u, c, d*u+e*z-u*z, t, v, w, y], [c, e, t, u, v, x, y], [b, c-t, d, e-u, w, x, y], [b-z, c, d*u+e*z-u*z, v, w, x, y], [b*t+c*z-t*z, d-z, e, v, w, x, y], [c, d*u+e*z-u*z, t, v, w, x, y], [b*t+c*z-t*z, e, u, v, w, x, y], [a, b*v+c*x-v*x, d, e-v, t, w-x, z], [a, b, c-v, d*v+e*y-v*y, u, w-y, z], [y*a+(w-y)*t, b, d, e-v, u, x, z], [a, e, t, u, v, x, z], [a-t, b, d, e-v, w, x, z], [x*a+(w-x)*u, c, d, t, v, y, z], [a, c, t, u, v, y, z], [a-u, b, c-v, d, w, y, z], [a-u, b*v+c*x-v*x, d, t, w, y, z], [b, d, t, u, x, y, z], [d, t, v, w, x, y, z], [b, u, v, w, x, y, z], [a, b-z, c, d-z, e, v, w-y, x-y], [a-t, b-z, c, d-z, e, v, w, x], [a, b*v+c*z-v*z, d-z, e, t-v, w-z, x-z, y-z], [a, b-z, c, d*v+e*z-v*z, u-v, w-z, x-z, y-z], [a-t, b, c-t, d*v+e*z-v*z, u-v, w, x, y-z], [a-u, b-z, c, d-z, e, v, w, y], [a-u, b*v+c*z-v*z, d, e-u, t-v, w, x-z, y], [b, c, d, e, t, u, x, y], [b-z, c, d-z, e, v, w, x, y], [c, d-z, e, t, v, w, x, y], [b-z, c, e, u, v, w, x, y], [a, d, e, t, v, w, x, z], [x*a+(w-x)*v, b, c-v, d, t, u-v, y, z], [a, b, c-v, d, t, u, y, z], [a-u, c, d, t, v, w, y, z], [a, b, c, u, v, w, y, z]]; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; lrs1 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); lrs2 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs1, R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs2, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) then return true else return false end if end proc:";
//    cstr = new char[procStr.length()+1];
//    std::strcpy (cstr, procStr.c_str());
//    ALGEB testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;
//
//    ALGEB result = EvalMapleProc(kv, testProc, 1, res);
//
//	if (mapleTest->testEquality(result, "true"))
//		return true;
//	else
//		return false;
//}

bool triangularizeValidate(vector<SMQP> F, RegularChain<RN,SMQP> rc, vector<Symbol> algvars, vector<Symbol> transVars, bool showOutput, bool isLazard) {
    vector<RegularChain<RN,SMQP>> results;

    // Timing variables
    long long unsigned int start;
    float elapsed;

    // Open timing file for output
    fstream fs,fs2,fs3;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs2.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs3.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);


    // Redirect stdout
    std::fstream f("/dev/null");

    //not thread-safe
    // stringstream ss;
    streambuf *oldcout,*oldcerr;

    if (!showOutput) {
        oldcout = cout.rdbuf(); // <-- save
        oldcerr = cerr.rdbuf(); // <-- save

        cout.rdbuf(f.rdbuf());
        cerr.rdbuf(f.rdbuf());
        // cout.rdbuf (ss.rdbuf());
        // cerr.rdbuf (ss.rdbuf());
    }
    #if REGULARCHAIN_DEBUG
    for (int i=0; i<F.size(); ++i) {
        cerr << "F[" << i << "] = " << F[i] << endl;
    }
    #endif

    startTimer(&start);

	results = rc.triangularize(F,isLazard,0);

	stopTimer(&start,&elapsed);
	fs << elapsed << "\t";
	fs2 << elapsed << "\t";
	fs3 << elapsed << std::endl;
	fs3.close();

    if (!showOutput) {
		cout.rdbuf (oldcout);   // <-- restore
		cerr.rdbuf (oldcerr);   // <-- restore
        f.close();
	}

	std::cerr << "BPAS result:" << std::endl;
	std::cerr << "[";
	for (int k=0; k<results.size(); ++k) {
		std::cerr << results[k];
		if (k!=results.size()-1)
			cout << ", ";
	}
	std::cerr << "]" << std::endl;


    ExpressionTree FTree;
    FTree.fromVector<SMQP>(F);
    ExpressionTree rcTree = rc.convertToExpressionTree();
    ExpressionTree RTree;
    RTree.fromVector<Symbol>(algvars);
    ExpressionTree ParamTree;
    ParamTree.fromVector<Symbol>(transVars);
    ExpressionTree resultTrees;
    resultTrees.fromVector<RegularChain<RN,SMQP>>(results);

    std::vector<std::string> inputs;
    inputs.push_back(FTree.toMapleString() + ":");
    inputs.push_back(rcTree.toMapleString() + ":");
    inputs.push_back(RTree.toMapleString() + ":");
    inputs.push_back(ParamTree.toMapleString() + ":");
    inputs.push_back(resultTrees.toMapleString() + ":");


#if defined(MAPLE_VALIDATE) && MAPLE_VALIDATE
    //Jan23/2020:
    MapleInterfaceStream& mis = MapleInterfaceStream::instance();
    return mis.TriangularizeValidation(isLazard, inputs);
#else
    return true;
#endif
// #if defined(WITH_MAPLE) && WITH_MAPLE
//    fs << "0\t0" << std::endl;
//    fs2 << "0\t0" << std::endl;
//    fs.close();
//    fs2.close();
//    MapleInterfaceStream& mis = MapleInterfaceStream::instance();
//    bool passed = mis.TriangularizeValidation(isLazard, inputs);
//    return passed
// #elif defined(SERIAL) && SERIAL && defined(MAPLE_VALIDATE) && MAPLE_VALIDATE

#if 0
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();

    char* cstr;
    std::string evalStr;
    evalStr = "kernelopts(numcpus=1);";
    cstr = new char[evalStr.length()+1];
    std::strcpy (cstr, evalStr.c_str());
    EvalMapleStatement(kv, cstr);
    delete[] cstr;

	startTimer(&start);

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }

    //algebList: [0] = F, [1] = T, [2] = Rlist, [3] = resultChainList

//    std::string procStr = "TriangularizeValidate := proc (F::list, in_rc::list, Rlist::list, results::list) local lrc1, lrc2, n, rc, R, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(results); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, results)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc2 := RegularChains:-Triangularize(F, rc, R); if evalb(nops(lrc2) = 0) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) then return true end if; lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
//    std::string procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard'); n := nops(lrc); results := []; if evalb(nops(lrc) = 0) then results := [op(results), []] else for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do end if; results end proc:";
    std::string procStr;
    if(isLazard) {
        procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard','radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
    }
    else {
        procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R,'radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
    }
    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    ALGEB result = EvalMapleProc(kv, testProc, 3, algebList[0], algebList[1], algebList[2]);

    stopTimer(&start,&elapsed);
    fs << elapsed << "\t";
    fs2 << elapsed << "\t";

//	std::cerr << "Triangularize result: " << std::endl;
//	std::cerr << mapleTest->algebToString(kv, result) << std::endl;
//     std::cerr << "calling maple proc: \n\n";
//     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
//     std::cerr << "\n\n";

	std::vector<ALGEB> algebList2;
	algebList2.push_back(algebList[3]);
	algebList2.push_back(result);
	algebList2.push_back(algebList[2]);

	if (isLazard) {
	    procStr = "EqualAsConstructibleSets := proc (dec1::list, dec2::list, Rlist::list) local n, lrc1, lrc2, lrs1, lrs2, cs1, cs2, R, i, rc; R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; if evalb(nops(lrc1) = 1) and evalb(nops(lrc2) = 1) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc2), R) then return true end if; lrs1 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); lrs2 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs1, R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs2, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
    }
    else {
	      procStr = "KalkbrenerTest := proc(dec1::list,dec2::list,Rlist::list) local rc,R,lsi1,lsi2,li,h,J,J1,J2,i,j,n,m,pass::boolean; if evalb(nops(dec1) = 0) and evalb(nops(dec2) = 0) then return true; end if;  R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lsi1 := []; writeline(terminal,convert(ConstructSaturatedIdeals,string)); writeline(terminal,convert(SaturatedIdealsForFirstDecomposition,string)); for i to n do writeline(terminal,convert(computingSaturatedIdeal,string)); rc := dec1[i]; m := nops(rc); h := 1; for j to m-1 do  h := h * RegularChains:-Initial(rc[j],R); end do; J := PolynomialIdeals:-PolynomialIdeal(rc); J := PolynomialIdeals:-Saturate(J,h); lsi1 := [op(lsi1), J]; end do; n := nops(dec2); lsi2 := []; writeline(terminal,convert(SaturatedIdealsForSecondDecomposition,string)); for i to n do writeline(terminal,convert(computingSaturatedIdeal,string)); rc := dec2[i]; m := nops(rc); h := 1; for j to m-1 do h := h * RegularChains:-Initial(rc[j],R); end do; J := PolynomialIdeals:-PolynomialIdeal(rc); J := PolynomialIdeals:-Saturate(J,h); lsi2 := [op(lsi2), J]; end do; writeline(terminal,convert(BeginMutualInclusionOfSaturatedIdealsTest,string)); n := nops(lsi1); m := nops(lsi2); for i to n do for j to m do if PolynomialIdeals:-IdealContainment(lsi1[i],lsi2[j],lsi1[i]) then pass := true; break; else pass := false; end if; end do; if evalb(pass) then next; else break; end if; end do; if evalb(not pass) then writeline(terminal,convert(SaturatedIdealsFromFirstListNotIncludedInSecond,string)); end if; if evalb(pass) then  for i to n do for j to m do if PolynomialIdeals:-IdealContainment(lsi1[i],lsi2[j],lsi1[i]) then pass := true; break; else pass := false; end if; end do; if evalb(pass) then next; else break; end if; end do; if evalb(pass) then return true; end if; end if; writeline(terminal,convert(SaturatedIdealsFromSecondListNotIncludedInFirst,string)); writeline(terminal,convert(BeginIdentityOfIntersectionsOfSaturatedIdealsTest,string)); J1 := lsi1[1]; if n > 1 then for i from 2 to n do J1 := PolynomialIdeals:-Intersect(J1,lsi1[i]); end do; end if; J2 := lsi2[1]; if m > 1 then for i from 2 to m do J2 := PolynomialIdeals:-Intersect(J2,lsi2[i]); end do; end if; if PolynomialIdeals:-IdealContainment(J1,J2,J1) then return true; end if; for i from 1 to n do G := PolynomialIdeals:-Generators(PolynomialIdeals:-Simplify(lsi1[i])); for j from 1 to m do pass := true; for g in G do if PolynomialIdeals:-RadicalMembership(g, lsi2[j]) then else pass := false; break; end if; end do; if evalb(pass) = true then break; end if; end do; if evalb(pass) = false then break;	end if; end do; if evalb(pass) = true then for j from 1 to m do G := PolynomialIdeals:-Generators(PolynomialIdeals:-Simplify(lsi2[j])); for i from 1 to n do pass := true; for g in G do if PolynomialIdeals:-RadicalMembership(g, lsi1[i]) then else pass := false; break;	end if; end do;	if evalb(pass) = true then break; end if; end do; if evalb(pass) = false then break; end if; end do; end if; if evalb(pass) = true then	return true; end if; J1 := PolynomialIdeals:-Radical(J1); J2 := PolynomialIdeals:-Radical(J2); if PolynomialIdeals:-IdealContainment(J1,J2,J1) then	return true; else return false;	end if; end proc:";

//        procStr = "with(PolynomialIdeals):";
//        cstr = new char[procStr.length()+1];
//        std::strcpy (cstr, procStr.c_str());
//        testProc = EvalMapleStatement(kv, cstr);
//        delete[] cstr;
//        procStr = "KalkbrenerTest := proc(dec1::list,dec2::list,Rlist::list) local rc,lrc1,lrc2,lsi1,lsi2,li,h,R,item,J,J1,J2,i,j,n,m,pass::boolean; R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc]; end do; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; n := nops(lrc1); lsi1 := []; for i to n do li := RegularChains:-Inequations(lrc1[i],R); h := 1; for item in li do h := h * item; end do; J := PolynomialIdeals:-PolynomialIdeal(RegularChains:-Equations(lrc1[i],R)); J := PolynomialIdeals:-Saturate(J,h); lsi1 := [op(lsi1), J]; end do; n := nops(lrc2); lsi2 := []; for i to n do li := RegularChains:-Inequations(lrc2[i],R); h := 1; for item in li do h := h * item; end do; J := PolynomialIdeals:-PolynomialIdeal(RegularChains:-Equations(lrc2[i],R)); J := PolynomialIdeals:-Saturate(J,h); lsi2 := [op(lsi2), J]; end do; n := nops(lsi1); m := nops(lsi2); for i to n do for j to m do if (lsi1[i] subset lsi2[j] and lsi2[j] subset lsi1[i]) then pass := 'true'; break; else pass := 'false'; end if; end do; if eval(pass) then next; else break; end if; end do; if not eval(pass) then end if; if eval(pass) then for i to n do for j to m do if (lsi1[i] subset lsi2[j] and lsi2[j] subset lsi1[i]) then pass := "true"; break; else pass := "false"; end if; end do; if eval(pass) then next; else break; end if; end do; if eval(pass) then return true; end if; end if; J1 := lsi1[1]; if n > 1 then for i from 2 to n do J1 := PolynomialIdeals:-Intersect(J1,lsi1[i]); end do; end if; J2 := lsi2[1]; if m > 1 then for i from 2 to m do J2 := PolynomialIdeals:-Intersect(J2,lsi2[i]); end do; end if; if J1 subset J2 and J2 subset J1 then return true; else return false; end if: end proc:";
    }

    evalStr = "kernelopts(numcpus=24);";
    cstr = new char[evalStr.length()+1];
    std::strcpy (cstr, evalStr.c_str());
    EvalMapleStatement(kv, cstr);
    delete[] cstr;

    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    startTimer(&start);

    result = EvalMapleProc(kv, testProc, 3, algebList2[0], algebList2[1], algebList2[2]);

    stopTimer(&start,&elapsed);
    fs << elapsed << std::endl;
    fs2 << elapsed << std::endl;

    fs.close();
    fs2.close();

    if (mapleTest->testEquality(result, "true")) {
//  if (true) {
        std::cerr << "BPAS result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[0]) << "\n";
        std::cerr << "Maple result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[1]) << "\n";
        return true;
    }
    else {
        std::cerr << "BPAS result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[0]) << "\n";
        std::cerr << "Maple result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[1]) << "\n";
//      std::cerr << mapleTest->algebToString(kv, algebList2[2]) << "\n";
        std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
        return false;
    }
//#else
//    fs << "0\t0" << std::endl;
//    fs2 << "0\t0" << std::endl;
//    fs.close();
//    fs2.close();
//    std::cerr << "BPAS result:" << std::endl;
//    std::cerr << "" << "\n";
//    std::cerr << "Maple result:" << std::endl;
//    std::cerr << "" << "\n";
//    return true;
//#endif
#endif
}

bool triangularizeValidate(vector<SMQP> F, RegularChain<RN,SMQP> rc, vector<Symbol> vars, bool showOutput, bool isLazard) {
    vector<RegularChain<RN,SMQP>> results;

    // Timing variables
    long long unsigned int start;
    float elapsed;

    // Open timing file for output
    fstream fs,fs2,fs3;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs2.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs3.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);


    // Redirect stdout
    std::fstream f("/dev/null");

    //not thread-safe
    // stringstream ss;
    streambuf *oldcout,*oldcerr;

    if (!showOutput) {
        oldcout = cout.rdbuf(); // <-- save
        oldcerr = cerr.rdbuf(); // <-- save

        cout.rdbuf(f.rdbuf());
        cerr.rdbuf(f.rdbuf());
        // cout.rdbuf (ss.rdbuf());
        // cerr.rdbuf (ss.rdbuf());
    }
    #if REGULARCHAIN_DEBUG
    for (int i=0; i<F.size(); ++i) {
        cerr << "F[" << i << "] = " << F[i] << endl;
    }
    #endif

    startTimer(&start);

	results = rc.triangularize(F,isLazard,0);

	stopTimer(&start,&elapsed);
	fs << elapsed << "\t";
	fs2 << elapsed << "\t";
	fs3 << elapsed << std::endl;
	fs3.close();

    if (!showOutput) {
		cout.rdbuf (oldcout);   // <-- restore
		cerr.rdbuf (oldcerr);   // <-- restore
        f.close();
	}

	std::cerr << "BPAS result:" << std::endl;
	std::cerr << "[";
	for (int k=0; k<results.size(); ++k) {
		std::cerr << results[k];
		if (k!=results.size()-1)
			cout << ", ";
	}
	std::cerr << "]" << std::endl;

    ExpressionTree FTree;
    FTree.fromVector<SMQP>(F);
    ExpressionTree rcTree = rc.convertToExpressionTree();
    ExpressionTree RTree;
    RTree.fromVector<Symbol>(vars);
    ExpressionTree resultTrees;
    resultTrees.fromVector<RegularChain<RN,SMQP>>(results);

    std::vector<std::string> inputs;
    inputs.push_back(FTree.toMapleString() + ":");
    inputs.push_back(rcTree.toMapleString() + ":");
    inputs.push_back(RTree.toMapleString() + ":");
    inputs.push_back(resultTrees.toMapleString() + ":");

#if defined(MAPLE_VALIDATE) && MAPLE_VALIDATE
    //Jan23/2020:
    MapleInterfaceStream& mis = MapleInterfaceStream::instance();
    return mis.TriangularizeValidation(isLazard, inputs);
#else
    return true;
#endif

#if 0
// #if defined(WITH_MAPLE) && WITH_MAPLE && defined(MAPLE_VALIDATE) && MAPLE_VALIDATE
//    fs << "0\t0" << std::endl;
//    fs2 << "0\t0" << std::endl;
//    fs.close();
//    fs2.close();
//    MapleInterfaceStream& mis = MapleInterfaceStream::instance();
//    bool passed = mis.TriangularizeValidation(isLazard, inputs);
//    return passed
// #elif defined(SERIAL) && SERIAL && defined(MAPLE_VALIDATE) && MAPLE_VALIDATE

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();

    char* cstr;
    std::string evalStr;
    evalStr = "kernelopts(numcpus=1);";
    cstr = new char[evalStr.length()+1];
    std::strcpy (cstr, evalStr.c_str());
    EvalMapleStatement(kv, cstr);
    delete[] cstr;

	startTimer(&start);

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }

    //algebList: [0] = F, [1] = T, [2] = Rlist, [3] = resultChainList

//    std::string procStr = "TriangularizeValidate := proc (F::list, in_rc::list, Rlist::list, results::list) local lrc1, lrc2, n, rc, R, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(results); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, results)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc2 := RegularChains:-Triangularize(F, rc, R); if evalb(nops(lrc2) = 0) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) then return true end if; lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
//    std::string procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard'); n := nops(lrc); results := []; if evalb(nops(lrc) = 0) then results := [op(results), []] else for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do end if; results end proc:";
    std::string procStr;
    if(isLazard) {
        procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard','radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
    }
    else {
        procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R,'radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
    }
    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    ALGEB result = EvalMapleProc(kv, testProc, 3, algebList[0], algebList[1], algebList[2]);

    stopTimer(&start,&elapsed);
    fs << elapsed << "\t";
    fs2 << elapsed << "\t";

//	std::cerr << "Triangularize result: " << std::endl;
//	std::cerr << mapleTest->algebToString(kv, result) << std::endl;
//     std::cerr << "calling maple proc: \n\n";
//     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
//     std::cerr << "\n\n";

	std::vector<ALGEB> algebList2;
	algebList2.push_back(algebList[3]);
	algebList2.push_back(result);
	algebList2.push_back(algebList[2]);

	if (isLazard) {
	    procStr = "EqualAsConstructibleSets := proc (dec1::list, dec2::list, Rlist::list) local n, lrc1, lrc2, lrs1, lrs2, cs1, cs2, R, i, rc; R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; if evalb(nops(lrc1) = 1) and evalb(nops(lrc2) = 1) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc2), R) then return true end if; lrs1 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); lrs2 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs1, R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs2, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
    }
    else {
	      procStr = "KalkbrenerTest := proc(dec1::list,dec2::list,Rlist::list) local rc,R,lsi1,lsi2,li,h,J,J1,J2,i,j,n,m,pass::boolean; if evalb(nops(dec1) = 0) and evalb(nops(dec2) = 0) then return true; end if;  R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lsi1 := []; writeline(terminal,convert(ConstructSaturatedIdeals,string)); writeline(terminal,convert(SaturatedIdealsForFirstDecomposition,string)); for i to n do writeline(terminal,convert(computingSaturatedIdeal,string)); rc := dec1[i]; m := nops(rc); h := 1; for j to m-1 do  h := h * RegularChains:-Initial(rc[j],R); end do; J := PolynomialIdeals:-PolynomialIdeal(rc); J := PolynomialIdeals:-Saturate(J,h); lsi1 := [op(lsi1), J]; end do; n := nops(dec2); lsi2 := []; writeline(terminal,convert(SaturatedIdealsForSecondDecomposition,string)); for i to n do writeline(terminal,convert(computingSaturatedIdeal,string)); rc := dec2[i]; m := nops(rc); h := 1; for j to m-1 do h := h * RegularChains:-Initial(rc[j],R); end do; J := PolynomialIdeals:-PolynomialIdeal(rc); J := PolynomialIdeals:-Saturate(J,h); lsi2 := [op(lsi2), J]; end do; writeline(terminal,convert(BeginMutualInclusionOfSaturatedIdealsTest,string)); n := nops(lsi1); m := nops(lsi2); for i to n do for j to m do if PolynomialIdeals:-IdealContainment(lsi1[i],lsi2[j],lsi1[i]) then pass := true; break; else pass := false; end if; end do; if evalb(pass) then next; else break; end if; end do; if evalb(not pass) then writeline(terminal,convert(SaturatedIdealsFromFirstListNotIncludedInSecond,string)); end if; if evalb(pass) then  for i to n do for j to m do if PolynomialIdeals:-IdealContainment(lsi1[i],lsi2[j],lsi1[i]) then pass := true; break; else pass := false; end if; end do; if evalb(pass) then next; else break; end if; end do; if evalb(pass) then return true; end if; end if; writeline(terminal,convert(SaturatedIdealsFromSecondListNotIncludedInFirst,string)); writeline(terminal,convert(BeginIdentityOfIntersectionsOfSaturatedIdealsTest,string)); J1 := lsi1[1]; if n > 1 then for i from 2 to n do J1 := PolynomialIdeals:-Intersect(J1,lsi1[i]); end do; end if; J2 := lsi2[1]; if m > 1 then for i from 2 to m do J2 := PolynomialIdeals:-Intersect(J2,lsi2[i]); end do; end if; if PolynomialIdeals:-IdealContainment(J1,J2,J1) then return true; end if; for i from 1 to n do G := PolynomialIdeals:-Generators(PolynomialIdeals:-Simplify(lsi1[i])); for j from 1 to m do pass := true; for g in G do if PolynomialIdeals:-RadicalMembership(g, lsi2[j]) then else pass := false; break; end if; end do; if evalb(pass) = true then break; end if; end do; if evalb(pass) = false then break;	end if; end do; if evalb(pass) = true then for j from 1 to m do G := PolynomialIdeals:-Generators(PolynomialIdeals:-Simplify(lsi2[j])); for i from 1 to n do pass := true; for g in G do if PolynomialIdeals:-RadicalMembership(g, lsi1[i]) then else pass := false; break;	end if; end do;	if evalb(pass) = true then break; end if; end do; if evalb(pass) = false then break; end if; end do; end if; if evalb(pass) = true then	return true; end if; J1 := PolynomialIdeals:-Radical(J1); J2 := PolynomialIdeals:-Radical(J2); if PolynomialIdeals:-IdealContainment(J1,J2,J1) then	return true; else return false;	end if; end proc:";

//        procStr = "with(PolynomialIdeals):";
//        cstr = new char[procStr.length()+1];
//        std::strcpy (cstr, procStr.c_str());
//        testProc = EvalMapleStatement(kv, cstr);
//        delete[] cstr;
//        procStr = "KalkbrenerTest := proc(dec1::list,dec2::list,Rlist::list) local rc,lrc1,lrc2,lsi1,lsi2,li,h,R,item,J,J1,J2,i,j,n,m,pass::boolean; R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc]; end do; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; n := nops(lrc1); lsi1 := []; for i to n do li := RegularChains:-Inequations(lrc1[i],R); h := 1; for item in li do h := h * item; end do; J := PolynomialIdeals:-PolynomialIdeal(RegularChains:-Equations(lrc1[i],R)); J := PolynomialIdeals:-Saturate(J,h); lsi1 := [op(lsi1), J]; end do; n := nops(lrc2); lsi2 := []; for i to n do li := RegularChains:-Inequations(lrc2[i],R); h := 1; for item in li do h := h * item; end do; J := PolynomialIdeals:-PolynomialIdeal(RegularChains:-Equations(lrc2[i],R)); J := PolynomialIdeals:-Saturate(J,h); lsi2 := [op(lsi2), J]; end do; n := nops(lsi1); m := nops(lsi2); for i to n do for j to m do if (lsi1[i] subset lsi2[j] and lsi2[j] subset lsi1[i]) then pass := 'true'; break; else pass := 'false'; end if; end do; if eval(pass) then next; else break; end if; end do; if not eval(pass) then end if; if eval(pass) then for i to n do for j to m do if (lsi1[i] subset lsi2[j] and lsi2[j] subset lsi1[i]) then pass := "true"; break; else pass := "false"; end if; end do; if eval(pass) then next; else break; end if; end do; if eval(pass) then return true; end if; end if; J1 := lsi1[1]; if n > 1 then for i from 2 to n do J1 := PolynomialIdeals:-Intersect(J1,lsi1[i]); end do; end if; J2 := lsi2[1]; if m > 1 then for i from 2 to m do J2 := PolynomialIdeals:-Intersect(J2,lsi2[i]); end do; end if; if J1 subset J2 and J2 subset J1 then return true; else return false; end if: end proc:";
    }

    evalStr = "kernelopts(numcpus=24);";
    cstr = new char[evalStr.length()+1];
    std::strcpy (cstr, evalStr.c_str());
    EvalMapleStatement(kv, cstr);
    delete[] cstr;

    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    startTimer(&start);

    result = EvalMapleProc(kv, testProc, 3, algebList2[0], algebList2[1], algebList2[2]);

    stopTimer(&start,&elapsed);
    fs << elapsed << std::endl;
    fs2 << elapsed << std::endl;

    fs.close();
    fs2.close();

    if (mapleTest->testEquality(result, "true")) {
//  if (true) {
        std::cerr << "BPAS result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[0]) << "\n";
        std::cerr << "Maple result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[1]) << "\n";
        return true;
    }
    else {
        std::cerr << "BPAS result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[0]) << "\n";
        std::cerr << "Maple result:" << std::endl;
        std::cerr << mapleTest->algebToString(kv, algebList2[1]) << "\n";
//      std::cerr << mapleTest->algebToString(kv, algebList2[2]) << "\n";
        std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
        return false;
    }
//#else
//    fs << "0\t0" << std::endl;
//    fs2 << "0\t0" << std::endl;
//    fs.close();
//    fs2.close();
//    std::cerr << "BPAS result:" << std::endl;
//    std::cerr << "" << "\n";
//    std::cerr << "Maple result:" << std::endl;
//    std::cerr << "" << "\n";
//    return true;
//#endif
#endif
}

void Sys126Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(2)*(x_1^2)-(x_2^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)*(x_1^2)-(x_1)*(x_2)-(x_2)*(c_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "126" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "126" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "126" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys126Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys126Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys130Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber("-27577200")*(x_1^23)*(x_2)-RationalNumber("920137680")*(x_1^22)*(x_2)-RationalNumber("13283391240")*(x_1^21)*(x_2)-RationalNumber("106624468368")*(x_1^20)*(x_2)-RationalNumber("491746186572")*(x_1^19)*(x_2)-RationalNumber("993778599636")*(x_1^18)*(x_2)+RationalNumber("2476605785622")*(x_1^17)*(x_2)+RationalNumber("25815018153060")*(x_1^16)*(x_2)+RationalNumber("96112956427272")*(x_1^15)*(x_2)+RationalNumber("210514577654592")*(x_1^14)*(x_2)+RationalNumber("255834944199876")*(x_1^13)*(x_2)+RationalNumber("1673661361104")*(x_1^12)*(x_2)-RationalNumber("682308301717818")*(x_1^11)*(x_2)-RationalNumber("1476465320622258")*(x_1^10)*(x_2)-RationalNumber("1610410635303342")*(x_1^9)*(x_2)-RationalNumber("597224233919088")*(x_1^8)*(x_2)+RationalNumber("938901769037796")*(x_1^7)*(x_2)+RationalNumber("1697021033615460")*(x_1^6)*(x_2)+RationalNumber("1209904120743648")*(x_1^5)*(x_2)+RationalNumber("297722035302048")*(x_1^4)*(x_2)-RationalNumber("166331517404934")*(x_1^3)*(x_2)-RationalNumber("158840938909422")*(x_1^2)*(x_2)-RationalNumber("50605101044820")*(x_1)*(x_2)-RationalNumber("6044877890100")*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("1185408")*(x_1^21)*(x_2)+RationalNumber("59185728")*(x_1^20)*(x_2)+RationalNumber("1351691712")*(x_1^19)*(x_2)+RationalNumber("18700782336")*(x_1^18)*(x_2)+RationalNumber("174871916064")*(x_1^17)*(x_2)+RationalNumber("1166498753904")*(x_1^16)*(x_2)+RationalNumber("5699486912112")*(x_1^15)*(x_2)+RationalNumber("20539343601792")*(x_1^14)*(x_2)+RationalNumber("53774942979384")*(x_1^13)*(x_2)+RationalNumber("96423825757788")*(x_1^12)*(x_2)+RationalNumber("93804944335092")*(x_1^11)*(x_2)-RationalNumber("40214503702512")*(x_1^10)*(x_2)-RationalNumber("316349629027650")*(x_1^9)*(x_2)-RationalNumber("570981123783459")*(x_1^8)*(x_2)-RationalNumber("569250985467807")*(x_1^7)*(x_2)-RationalNumber("304838378563104")*(x_1^6)*(x_2)-RationalNumber("81135447444378")*(x_1^5)*(x_2)+RationalNumber("16825934075661")*(x_1^4)*(x_2)+RationalNumber("38057921178669")*(x_1^3)*(x_2)+RationalNumber("16575248060910")*(x_1^2)*(x_2)+RationalNumber("1990386622350")*(x_1)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "130" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "130" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "130" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys130Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys130Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys132Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(-1)+(x_2)-(x_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-2)-RationalNumber(4)*(x_1)*(x_2)-RationalNumber(2)*(x_1)-RationalNumber(2)*(x_2^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)-RationalNumber(2)*(x_1)*(x_2^2)-RationalNumber(4)*(x_1)-RationalNumber(4)*(x_1)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(1)-(x_1)*(x_2)+(x_2^2)-RationalNumber(3)*(x_2)-RationalNumber(2)*(x_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "132" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "132" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "132" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys132Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys132Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}


void Sys1000Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+(x_1)+RationalNumber(2)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)+RationalNumber(2)*(x_1)+RationalNumber(3)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1000" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1000" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1000" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1000Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1000Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1001Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+(x_1)+RationalNumber(3)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)+(x_1)+RationalNumber(3)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1001" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1001" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1001" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1001Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1001Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1002Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(8)+(x_1^2)-(x_2^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-16)*(x_1^2)-RationalNumber(48)-RationalNumber(16)*(x_2)-(x_1^4)-RationalNumber(4)*(x_1^2)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1002" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1002" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1002" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
}
	if (pass)
		std::cerr << "RC Sys1002Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1002Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1003Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = (x_1^2)*(x_2)+RationalNumber(3)*(x_1)+RationalNumber(2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_1^2)+(x_1)*(x_2)+RationalNumber(1)-(x_1)+(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1003" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1003" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1003" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1003Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1003Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1004Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(3)+RationalNumber(3)*(x_1)+(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(3)+RationalNumber(3)*(x_1)+RationalNumber(2)*(x_1^2)+(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1004" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1004" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1004" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1004Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1004Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1005Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(2)+RationalNumber(2)*(x_1)+(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)+(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1005" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1005" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1005" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1005Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1005Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1006Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+(x_1)+RationalNumber(5)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(4)+RationalNumber(5)*(x_1)+RationalNumber(4)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1006" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1006" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1006" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1006Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1006Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1007Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+RationalNumber(4)*(x_1)+RationalNumber(2)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)+RationalNumber(2)*(x_1)+RationalNumber(5)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1007" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1007" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1007" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1007Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1007Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1008Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = (x_1)*(x_2)+RationalNumber(2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(4)+(c_1)-(x_2)*(x_1^2)-(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1008" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1008" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1008" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1008Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1008Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1009Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+RationalNumber(2)*(x_1)+(x_1^2)+(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(3)+RationalNumber(2)*(x_1)+RationalNumber(2)*(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1009" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1009" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1009" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1009Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1009Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1010Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(2)+RationalNumber(3)*(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(3)+RationalNumber(2)*(x_1)+RationalNumber(3)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1010" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1010" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1010" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1010Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1010Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1011Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+(x_1)+RationalNumber(2)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)+(x_1)+(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1011" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1011" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1011" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1011Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1011Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1012Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(1)+RationalNumber(2)*(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(3)*(x_1^3)+(x_1^2)+RationalNumber(3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1012" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1012" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1012" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1012Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1012Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1013Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(2)*(x_1^3)+(x_1)+RationalNumber(2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)+(x_1)+(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1013" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1013" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1013" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1013Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1013Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1014Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(4)+(x_1)+RationalNumber(2)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(5)+(x_1)+RationalNumber(3)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1014" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1014" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1014" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1014Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1014Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1015Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(3)+RationalNumber(5)*(x_1)+RationalNumber(2)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(4)+(x_1)+RationalNumber(4)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1015" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1015" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1015" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1015Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1015Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1016Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(3)+RationalNumber(3)*(x_1)+(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(4)+RationalNumber(3)*(x_1)+(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1016" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1016" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1016" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1016Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1016Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1221Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly,temp;
	std::cerr << "constructing first polynomial:" << std::endl;
	poly = RationalNumber(2)*(c_1)*(x_1^2);
	poly.setRingVariables(vars);
	temp = -(c_1)*(x_1)*(x_2);
	temp.setRingVariables(vars);
	poly += temp;
	poly += RationalNumber(4)*(x_1^3);
	temp = -RationalNumber(3)*(x_1^2)*(x_2);
	temp.setRingVariables(vars);
	poly += temp;
	polys.push_back(poly);
	std::cerr << "constructing second polynomial:" << std::endl;
	poly = RationalNumber(2)*(x_1^2)-(x_2^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;
	std::cerr << "calling triangularize:" << std::endl;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1221" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1221" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1221" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1221Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1221Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1255Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("c_2");
	SMQP c_2("c_2");
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = (x_1^2)-RationalNumber(2)*(x_1)*(c_1)+RationalNumber(1)-(x_2^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-2)*(c_2)*(x_1)-RationalNumber(4)-RationalNumber(2)*(x_1)*(x_2)+RationalNumber(4)*(x_1)-(x_1^2)+RationalNumber(2)*(x_1)*(c_1)+RationalNumber(4)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1255" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1255" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1255" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1255Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1255Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1289Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("c_2");
	SMQP c_2("c_2");
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(125)*(c_1)*(c_2^2);
	poly += RationalNumber(125)*(x_2)*(c_1^2);
	poly += RationalNumber(-125)*(c_2)*(c_1^2);
	poly += RationalNumber(125)*(c_1);
	poly += RationalNumber(-15876)*(x_2);
	poly += RationalNumber(15751)*(c_2);
	poly += RationalNumber(-125)*(x_2)*(c_2^2);
	poly += RationalNumber(-125)*(c_1)*(x_2^2);
	poly += RationalNumber(125)*(x_2^2)*(c_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(15625)*(x_1^2)*(c_2);
	poly += RationalNumber(-15625)*(x_1)*(c_2^2);
	poly += RationalNumber(-15625)*(c_1)*(x_1^2);
	poly += RationalNumber(15625)*(c_1)*(c_2^2); // sign change bad
	poly += RationalNumber(15625)*(x_1)*(c_1^2);
	poly += RationalNumber(-15625)*(c_2)*(c_1^2);
	poly += RationalNumber(-1984500)*(x_1);
	poly += RationalNumber(63011844)*(c_2);
	poly += RationalNumber(-61027344)*(c_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(15625)*(x_1^2)*(x_2);
	poly += RationalNumber(-15625)*(x_1)*(x_2^2);
	poly += RationalNumber(-15625)*(c_1)*(x_1^2);
	poly += RationalNumber(15625)*(c_1)*(x_2^2);
	poly += RationalNumber(15625)*(x_1)*(c_1^2);
	poly += RationalNumber(-15625)*(x_2)*(c_1^2);
	poly += RationalNumber(-1968875)*(x_1);
	poly += RationalNumber(63011844)*(x_2);
	poly += RationalNumber(-61042969)*(c_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1289" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1289" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1289" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1289Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1289Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1302Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(6562)+RationalNumber(6562)*(x_2)-RationalNumber(7695)*(x_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-59058)+RationalNumber(10858)*(x_1)*(x_2)+RationalNumber(83748)*(x_1)-RationalNumber(52496)*(x_2)+RationalNumber(6562)*(x_2^2)-RationalNumber(23085)*(x_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-32810)-RationalNumber(79594)*(x_1)*(x_2)+RationalNumber(18553)*(x_1)*(x_2^2)-RationalNumber(6798)*(x_1^2)*(x_2)-RationalNumber(97862)*(x_1)+RationalNumber(175427)*(x_1^2)-RationalNumber(32810)*(x_2^2)-RationalNumber(65620)*(x_2)-RationalNumber(23085)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(6562)-RationalNumber(116130)*(x_1)*(x_2)-RationalNumber(58350)*(x_1)*(x_2^2)+RationalNumber(3304)*(x_1^2)*(x_2)+RationalNumber(16287)*(x_2^2)*(x_1^2)-RationalNumber(19922)*(x_1^3)*(x_2)-RationalNumber(68895)*(x_1)-RationalNumber(12413)*(x_1^2)+RationalNumber(125858)*(x_1^3)+RationalNumber(6562)*(x_2)-RationalNumber(7695)*(x_1^4);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-26248)-RationalNumber(730895)*(x_1)-RationalNumber(534495)*(x_1^2)+RationalNumber(1203563)*(x_1^3)+RationalNumber(6562)*(x_2)-RationalNumber(1459890)*(x_1)*(x_2)-RationalNumber(728995)*(x_1)*(x_2^2)-RationalNumber(340090)*(x_1^2)*(x_2)+RationalNumber(194405)*(x_2^2)*(x_1^2)-RationalNumber(109094)*(x_1^3)*(x_2)-RationalNumber(257)*(x_1^3)*(x_2^2)-RationalNumber(10348)*(x_1^4)*(x_2)+RationalNumber(95)*(x_1^4)*(x_2^2)-RationalNumber(43253)*(x_1^4);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(13124)+RationalNumber(19686)*(x_2)+RationalNumber(31023)*(x_1)+RationalNumber(285)*(x_1)*(x_2^2)+RationalNumber(4588)*(x_1)*(x_2)-RationalNumber(14635)*(x_2^2)*(x_1^2)-RationalNumber(28130)*(x_1^2)*(x_2)+RationalNumber(3163)*(x_1^3)*(x_2^2)+RationalNumber(35406)*(x_1^3)*(x_2)-RationalNumber(8828)*(x_1^4)*(x_2)-RationalNumber(41710)*(x_1^2)+RationalNumber(6562)*(x_2^2)+RationalNumber(32623)*(x_1^3)+RationalNumber(27617)*(x_1^4);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-682448)+RationalNumber(651978)*(x_1)-RationalNumber(643076)*(x_2)+RationalNumber(39372)*(x_2^2)-RationalNumber(18943)*(x_1^2)+RationalNumber(83706)*(x_1)*(x_2)+RationalNumber(15133)*(x_1)*(x_2^2)-RationalNumber(15918)*(x_1^2)*(x_2)+RationalNumber(570)*(x_2^2)*(x_1^2)+RationalNumber(29840)*(x_1^3)*(x_2)+RationalNumber(14540)*(x_1^3)*(x_2^2)+RationalNumber(5004)*(x_1^4)*(x_2)-RationalNumber(1133)*(x_1^4)*(x_2^2)-RationalNumber(11205)*(x_1^3)+RationalNumber(6232)*(x_1^4);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-360910)-RationalNumber(1300980)*(x_1)+RationalNumber(1610397)*(x_1^2)-RationalNumber(360910)*(x_2^2)-RationalNumber(1142070)*(x_1)*(x_2)+RationalNumber(157960)*(x_1)*(x_2^2)+RationalNumber(17404)*(x_1^2)*(x_2)+RationalNumber(11157)*(x_2^2)*(x_1^2)-RationalNumber(26002)*(x_1^3)*(x_2)+RationalNumber(380)*(x_1^3)*(x_2^2)+RationalNumber(7460)*(x_1^4)*(x_2)+RationalNumber(3635)*(x_1^4)*(x_2^2)-RationalNumber(721820)*(x_2)-RationalNumber(80672)*(x_1^3)-RationalNumber(4725)*(x_1^4);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1302" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1302" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1302" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1302Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1302Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1303Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("c_2");
	SMQP c_2("c_2");
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = (c_1^2)-RationalNumber(2)*(x_2^2)+(x_2)*(c_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_1^2)*(x_2)-(x_2^3)+(c_2)*(c_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(2)*(x_2^3)-(x_2)*(c_1^2)-(c_2)*(c_1^2)-RationalNumber(2)*(x_1)*(x_2^2)+(x_1)*(c_1^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1303" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1303" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1303" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1303Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1303Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1304Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(-270)*(x_1^4)*(x_2^3)-RationalNumber(314)*(x_1)*(x_2^4)-RationalNumber(689)*(x_1)*(x_2^3)+RationalNumber(1428);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(36)*(x_1^7)+RationalNumber(417)*(x_1^6)*(x_2)-RationalNumber(422)*(x_1^5)*(x_2^2)-RationalNumber(270)*(x_1^4)*(x_2^3)+RationalNumber(1428)*(x_1^3)*(x_2^4)-RationalNumber(1475)*(x_1^2)*(x_2^5)+RationalNumber(510)*(x_2^6)*(x_1)-RationalNumber(200)*(x_1^6)-RationalNumber(174)*(x_1^5)*(x_2)-RationalNumber(966)*(x_1^4)*(x_2^2)+RationalNumber(529)*(x_1^3)*(x_2^3)+RationalNumber(269)*(x_1^2)*(x_2^4)+RationalNumber(49)*(x_1)*(x_2^5)-RationalNumber(267)*(x_2^6)+RationalNumber(529)*(x_1^4)*(x_2)+RationalNumber(1303)*(x_1^2)*(x_2^3)-RationalNumber(314)*(x_1)*(x_2^4)+RationalNumber(262)*(x_2^5)+RationalNumber(36)*(x_1^4)-RationalNumber(788)*(x_2^2)*(x_1^2)-RationalNumber(689)*(x_1)*(x_2^3)+RationalNumber(177)*(x_2^4);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1304" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1304" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1304" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1304Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1304Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1305Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(3)+(x_1)+RationalNumber(2)*(x_1^2)+RationalNumber(2)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(3)+(x_1)+RationalNumber(2)*(x_1^2)+RationalNumber(3)*(x_1^3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1305" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1305" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1305" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1305Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1305Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1364Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(-1)*(x_1^5)*(x_2)-RationalNumber(8)*(x_1^4)*(x_2)-RationalNumber(20)*(x_1^3)*(x_2)-RationalNumber(17)*(x_1^2)*(x_2)-RationalNumber(6)*(x_1)*(x_2)-RationalNumber(8)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-1)*(x_1^8)*(x_2)-RationalNumber(14)*(x_1^7)*(x_2)-RationalNumber(79)*(x_1^6)*(x_2)-RationalNumber(231)*(x_1^5)*(x_2)-RationalNumber(376)*(x_1^4)*(x_2)-RationalNumber(353)*(x_1^3)*(x_2)-RationalNumber(230)*(x_1^2)*(x_2)-RationalNumber(152)*(x_1)*(x_2)-RationalNumber(64)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-2)*(x_1^6)*(x_2)-RationalNumber(26)*(x_1^5)*(x_2)-RationalNumber(134)*(x_1^4)*(x_2)-RationalNumber(350)*(x_1^3)*(x_2)-RationalNumber(488)*(x_1^2)*(x_2)-RationalNumber(344)*(x_1)*(x_2)-RationalNumber(96)*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1364" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1364" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1364" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1364Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1364Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1366Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber("12")+RationalNumber("4")*(x_2)-RationalNumber("6")*(x_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-12")-RationalNumber("4")*(x_1)*(x_2)+RationalNumber("2")*(x_2^2)-RationalNumber("3")*(x_1^2)+RationalNumber("12")*(x_1)+RationalNumber("8")*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("12")-RationalNumber("2")*(x_1)*(x_2^2)-RationalNumber("2")*(x_1^2)*(x_2)+RationalNumber("8")*(x_1)*(x_2)+RationalNumber("4")*(x_2^2)-RationalNumber("12")*(x_2)-RationalNumber("6")*(x_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-12")-(x_2^2)*(x_1^2)+RationalNumber("4")*(x_1)*(x_2^2)-RationalNumber("8")*(x_1)*(x_2)-RationalNumber("3")*(x_1^2)+RationalNumber("12")*(x_1)-RationalNumber("6")*(x_2^2)+RationalNumber("16")*(x_2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1366" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1366" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1366" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1366Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1366Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1373Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber("2")*(x_2)*(RationalNumber("-1")-RationalNumber("12")*(x_1^5)-RationalNumber("2")*(x_2^2)+RationalNumber("4")*(x_1)*(x_2^2)-RationalNumber("24")*(x_1^3)*(x_2^2)+RationalNumber("4")*(x_1)-RationalNumber("6")*(x_1^2)+RationalNumber("4")*(x_1^3)+RationalNumber("32")*(x_1^6)-RationalNumber("48")*(x_1^7)+RationalNumber("42")*(x_1^8)-RationalNumber("20")*(x_1^9)+RationalNumber("4")*(x_1^10)+RationalNumber("10")*(x_2^8)+RationalNumber("8")*(x_2^6)+RationalNumber("4")*(x_2^10)+RationalNumber("156")*(x_1^4)*(x_2^4)-RationalNumber("12")*(x_1)*(x_2^4)+RationalNumber("72")*(x_1^4)*(x_2^2)+RationalNumber("48")*(x_1^2)*(x_2^4)+RationalNumber("20")*(x_1^8)*(x_2^2)-RationalNumber("80")*(x_1^7)*(x_2^2)+RationalNumber("136")*(x_1^6)*(x_2^2)+RationalNumber("40")*(x_1^6)*(x_2^4)-RationalNumber("128")*(x_1^5)*(x_2^2)-RationalNumber("120")*(x_1^5)*(x_2^4)+RationalNumber("40")*(x_1^4)*(x_2^6)-RationalNumber("112")*(x_1^3)*(x_2^4)-RationalNumber("80")*(x_1^3)*(x_2^6)+RationalNumber("72")*(x_1^2)*(x_2^6)+RationalNumber("20")*(x_1^2)*(x_2^8)-RationalNumber("32")*(x_2^6)*(x_1)-RationalNumber("20")*(x_2^8)*(x_1));
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("24")*(x_1^5)-RationalNumber("2")*(x_2^4)+RationalNumber("6")*(x_1^4)-RationalNumber("4")*(x_1)*(x_2^2)+RationalNumber("4")*(x_2^2)*(x_1^2)+RationalNumber("48")*(x_1^3)*(x_2^2)-RationalNumber("2")*(x_1)+RationalNumber("8")*(x_1^2)-RationalNumber("12")*(x_1^3)-RationalNumber("84")*(x_1^6)+RationalNumber("144")*(x_1^7)-RationalNumber("156")*(x_1^8)+RationalNumber("108")*(x_1^9)-RationalNumber("44")*(x_1^10)+RationalNumber("8")*(x_1^11)-RationalNumber("12")*(x_2^8)-RationalNumber("12")*(x_2^6)-RationalNumber("4")*(x_2^10)-RationalNumber("440")*(x_1^4)*(x_2^4)+RationalNumber("24")*(x_1)*(x_2^4)-RationalNumber("180")*(x_1^4)*(x_2^2)-RationalNumber("108")*(x_1^2)*(x_2^4)-RationalNumber("180")*(x_1^8)*(x_2^2)+RationalNumber("368")*(x_1^7)*(x_2^2)-RationalNumber("448")*(x_1^6)*(x_2^2)-RationalNumber("280")*(x_1^6)*(x_2^4)+RationalNumber("352")*(x_1^5)*(x_2^2)+RationalNumber("456")*(x_1^5)*(x_2^4)-RationalNumber("200")*(x_1^4)*(x_2^6)+RationalNumber("272")*(x_1^3)*(x_2^4)+RationalNumber("240")*(x_1^3)*(x_2^6)-RationalNumber("160")*(x_1^2)*(x_2^6)-RationalNumber("60")*(x_1^2)*(x_2^8)+RationalNumber("64")*(x_2^6)*(x_1)+RationalNumber("44")*(x_2^8)*(x_1)+RationalNumber("40")*(x_1^9)*(x_2^2)+RationalNumber("80")*(x_1^7)*(x_2^4)+RationalNumber("80")*(x_1^5)*(x_2^6)+RationalNumber("40")*(x_1^3)*(x_2^8)+RationalNumber("8")*(x_1)*(x_2^10);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1373" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1373" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1373" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1373Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1373Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys1397Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	SMQP poly;
	poly = RationalNumber(5)*(x_1)*(x_2^4)+(x_2^4)+RationalNumber(2)*(x_2^3)+RationalNumber(2)*(x_2^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-1)*(x_2^2)-RationalNumber(2)*(x_1^5)-RationalNumber(4)*(x_2^3)-(x_1^2)*(x_2^3)+RationalNumber(2)*(x_1)*(x_2^4)+RationalNumber(3)*(x_1^3)*(x_2)+RationalNumber(64);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1397" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1397" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "1397" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys1397Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys1397Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys2852Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars,variables,params;
	vars.emplace_back("c_1");
	SMQP c_1("c_1");
	vars.emplace_back("c_2");
	SMQP c_2("c_2");
	vars.emplace_back("c_3");
	SMQP c_3("c_3");
	vars.emplace_back("c_4");
	SMQP c_4("c_4");
	vars.emplace_back("c_5");
	SMQP c_5("c_5");
	vars.emplace_back("x_10");
	SMQP x_10("x_10");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	vars.emplace_back("x_3");
	SMQP x_3("x_3");
	vars.emplace_back("x_4");
	SMQP x_4("x_4");
	vars.emplace_back("x_7");
	SMQP x_7("x_7");
	vars.emplace_back("x_8");
	SMQP x_8("x_8");
	vars.emplace_back("x_9");
	SMQP x_9("x_9");
	variables = vars;
//	variables.erase(variables.begin()+3,variables.begin()+5);
//	params.push_back("c_4");
//	params.push_back("c_5");
	SMQP poly;
	poly = c_1;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = c_2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = x_3;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = x_4;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber(-1)*(c_2)+(c_3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_7)-(c_3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_10^2)+(c_4^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_10^2)+(c_4^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_10)*(x_2);
	poly.setRingVariables(vars);
	poly += c_5;
	poly += c_4*c_1;
	poly.setRingVariables(vars);
	polys.push_back(poly);

	RegularChain<RationalNumber,SMQP> rc(vars);

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "2852" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "2852" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "2852" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys2852Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys2852Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void Sys2985Test(bool showOutput, bool isLazard) {
	std::vector<SMQP> polys;
	std::vector<Symbol> vars;
	vars.emplace_back("x_1");
	SMQP x_1("x_1");
	vars.emplace_back("x_2");
	SMQP x_2("x_2");
	vars.emplace_back("x_3");
	SMQP x_3("x_3");
	SMQP poly;
	poly = (x_1);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_3);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x_2^2)+RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	RegularChain<RationalNumber,SMQP> rc;

    // Output system name
    fstream fs;
    fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "2985" << "\t";
    fs.close();
    fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "2985" << "\t";
    fs.close();
    fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "2985" << "\t";
    fs.close();

	bool pass = triangularizeValidate(polys,rc,vars,showOutput,isLazard);

	if (showOutput) {
		std::cerr << "Input polynomials:" << std::endl;
		for (size_t i = 0; i < polys.size(); ++i)
			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
	}
	if (pass)
		std::cerr << "RC Sys2985Test:\t\t\t\t\t\t\t PASSED" << std::endl;
	else {
		std::cerr << "RC Sys2985Test:\t\t\t\t\t\t\t FAILED" << std::endl;
		exit(1);
	}
}

void getSampleSys2915(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("a");
	SparseMultivariateRationalPolynomial a("a");
	vars.emplace_back("t");
	SparseMultivariateRationalPolynomial t("t");
	vars.emplace_back("u");
	SparseMultivariateRationalPolynomial u("u");
	vars.emplace_back("v");
	SparseMultivariateRationalPolynomial v("v");
	vars.emplace_back("w");
	SparseMultivariateRationalPolynomial w("w");
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	SparseMultivariateRationalPolynomial poly;
	poly = (z)+(a);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x)+(y)-(w);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (t)+(u)+(v)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("2")*(v)*(w)+RationalNumber("2")*(u)*(a)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("3")*(v)*(w^2)+RationalNumber("3")*(u)*(a^2)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (y)*(v)*(a);
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys2920(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("b1");
	SparseMultivariateRationalPolynomial b1("b1");
	vars.emplace_back("t");
	SparseMultivariateRationalPolynomial t("t");
	vars.emplace_back("u");
	SparseMultivariateRationalPolynomial u("u");
	vars.emplace_back("v");
	SparseMultivariateRationalPolynomial v("v");
	vars.emplace_back("w");
	SparseMultivariateRationalPolynomial w("w");
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	SparseMultivariateRationalPolynomial poly;
	poly = (b1)+(y)+(z)-(t)-(w);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("2")*(z)*(u)+RationalNumber("2")*(y)*(v)+RationalNumber("2")*(t)*(w)-RationalNumber("2")*(w^2)-(w)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("3")*(z)*(u^2)+RationalNumber("3")*(y)*(v^2)-RationalNumber("3")*(t)*(w^2)+RationalNumber("3")*(w^3)+RationalNumber("3")*(w^2)-(t)+RationalNumber("4")*(w);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("6")*(x)*(z)*(v)-RationalNumber("6")*(t)*(w^2)+RationalNumber("6")*(w^3)-RationalNumber("3")*(t)*(w)+RationalNumber("6")*(w^2)-(t)+RationalNumber("4")*(w);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("4")*(z)*(u^3)+RationalNumber("4")*(y)*(v^3)+RationalNumber("4")*(t)*(w^3)-RationalNumber("4")*(w^4)-RationalNumber("6")*(w^3)+RationalNumber("4")*(t)*(w)-RationalNumber("10")*(w^2)-(w)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("8")*(x)*(z)*(u)*(v)+RationalNumber("8")*(t)*(w^3)-RationalNumber("8")*(w^4)+RationalNumber("4")*(t)*(w^2)-RationalNumber("12")*(w^3)+RationalNumber("4")*(t)*(w)-RationalNumber("14")*(w^2)-RationalNumber("3")*(w)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("12")*(x)*(z)*(v^2)+RationalNumber("12")*(t)*(w^3)-RationalNumber("12")*(w^4)+RationalNumber("12")*(t)*(w^2)-RationalNumber("18")*(w^3)+RationalNumber("8")*(t)*(w)-RationalNumber("14")*(w^2)-(w)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-24")*(t)*(w^3)+RationalNumber("24")*(w^4)-RationalNumber("24")*(t)*(w^2)+RationalNumber("36")*(w^3)-RationalNumber("8")*(t)*(w)+RationalNumber("26")*(w^2)+RationalNumber("7")*(w)+RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys2926(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	std::vector<SparseMultivariateRationalPolynomial> newPolys = {
	SMQP("2*y*w-z*w+t*w"),
	SMQP("-2*u*w^2+10*v*w^2-20*w^3+7*t*u-35*t*v+10*t*w,"),
	SMQP("2*y*w^2-2*z*w^2+6*t*w^2-7*y*t+7*z*t-21*t^2"),
	SMQP("-2*v^3+4*u*v*w+5*v^2*w-6*u*w^2-7*v*w^2+15*w^3+42*y*v-14*z*v-63*y*w+21*z*w-42*t*w+147*x"),
	SMQP("-9*u*w^3+45*v*w^3-135*w^4+14*z*v^2-14*t*v^2-28*z*u*w+70*t*u*w-14*z*v*w-196*t*v*w+28*z*w^2+602*t*w^2-294*y*z+98*z^2+294*y*t-98*z*t-147*x*u+735*x*v-2205*x*w"),
	SMQP("6*y*w^3-9*z*w^3+36*t*w^3-14*x*v^2-28*y*t*w+42*z*t*w-168*t^2*w+28*x*u*w+14*x*v*w-28*x*w^2+392*x*y-245*x*z+588*x*t"),
	SMQP("2*u*v*w-6*v^2*w-u*w^2+13*v*w^2-5*w^3+14*y*w-28*t*w"),
	SMQP("u^2*w-3*u*v*w+5*u*w^2+14*y*w-28*t*w"),
	SMQP("-2*z*u*w-2*t*u*w+4*y*v*w+6*z*v*w-2*t*v*w-16*y*w^2-10*z*w^2+22*t*w^2+42*x*w"),
	SMQP("28*y*u*w+8*z*u*w-20*t*u*w-88*y*v*w-24*z*v*w+68*t*v*w+156*y*w^2+40*z*w^2-132*t*w^2-252*x*w"),
	SMQP("-4*y*z*w+10*y*t*w+8*z*t*w-20*t^2*w+12*x*u*w-30*x*v*w+15*x*w^2"),
	SMQP("-2*y^2*w+y*z*w+2*y*t*w-2*z*t*w+4*t^2*w-6*x*u*w+12*x*v*w-6*x*w^2"),
	SMQP("8*x*y*w-4*x*z*w+8*x*t*w")};
	std::vector<Symbol> newVars = {Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w"), Symbol("x"), Symbol("y"), Symbol("z")};
	for (int i=0; i<newPolys.size(); ++i)
		newPolys[i].setRingVariables(newVars);
	vars = newVars;
	polys = newPolys;
//	vars.emplace_back("t");
//	SparseMultivariateRationalPolynomial t("t");
//	vars.emplace_back("u");
//	SparseMultivariateRationalPolynomial u("u");
//	vars.emplace_back("v");
//	SparseMultivariateRationalPolynomial v("v");
//	vars.emplace_back("w");
//	SparseMultivariateRationalPolynomial w("w");
//	vars.emplace_back("x");
//	SparseMultivariateRationalPolynomial x("x");
//	vars.emplace_back("y");
//	SparseMultivariateRationalPolynomial y("y");
//	vars.emplace_back("z");
//	SparseMultivariateRationalPolynomial z("z");
//	SparseMultivariateRationalPolynomial poly;
//	poly = RationalNumber("2")*(y)*(w);
//	poly.setRingVariables(vars);
//	poly -= -(z)*(w);
//	poly += (t)*(w);
//	polys.push_back(poly);
//	poly = RationalNumber("-2")*(u)*(w^2)+RationalNumber("10")*(v)*(w^2);
//	poly.setRingVariables(vars);
//	poly -= RationalNumber("20")*(w^3);
//	poly += RationalNumber("7")*(t)*(u);
//	poly -= RationalNumber("35")*(t)*(v);
//	poly += RationalNumber("10")*(t)*(w);
//	polys.push_back(poly);
//	poly = RationalNumber("2")*(y)*(w^2)-RationalNumber("2")*(z)*(w^2)+RationalNumber("6")*(t)*(w^2)-RationalNumber("7")*(y)*(t)+RationalNumber("7")*(z)*(t)-RationalNumber("21")*(t^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("-2")*(v^3)+RationalNumber("4")*(u)*(v)*(w);
//	poly.setRingVariables(vars);
//	poly += RationalNumber("5")*(v^2)*(w);
//	poly -= RationalNumber("6")*(u)*(w^2);
//	poly -= RationalNumber("7")*(v)*(w^2);
//	poly += RationalNumber("15")*(w^3);
//	poly += RationalNumber("42")*(y)*(v);
//	poly -= RationalNumber("14")*(z)*(v);
//	poly -= RationalNumber("63")*(y)*(w);
//	poly += RationalNumber("21")*(z)*(w);
//	poly -= RationalNumber("42")*(t)*(w);
//	poly += RationalNumber("147")*(x);
//	polys.push_back(poly);
//	poly = RationalNumber("-9")*(u)*(w^3)+RationalNumber("45")*(v)*(w^3)-RationalNumber("135")*(w^4)+RationalNumber("14")*(z)*(v^2)-RationalNumber("14")*(t)*(v^2)-RationalNumber("28")*(z)*(u)*(w)+RationalNumber("70")*(t)*(u)*(w)-RationalNumber("14")*(z)*(v)*(w)-RationalNumber("196")*(t)*(v)*(w)+RationalNumber("28")*(z)*(w^2)+RationalNumber("602")*(t)*(w^2)-RationalNumber("294")*(y)*(z)+RationalNumber("98")*(z^2)+RationalNumber("294")*(y)*(t)-RationalNumber("98")*(z)*(t)-RationalNumber("147")*(x)*(u)+RationalNumber("735")*(x)*(v)-RationalNumber("2205")*(x)*(w);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("6")*(y)*(w^3)-RationalNumber("9")*(z)*(w^3)+RationalNumber("36")*(t)*(w^3)-RationalNumber("14")*(x)*(v^2)-RationalNumber("28")*(y)*(t)*(w)+RationalNumber("42")*(z)*(t)*(w)-RationalNumber("168")*(t^2)*(w)+RationalNumber("28")*(x)*(u)*(w)+RationalNumber("14")*(x)*(v)*(w)-RationalNumber("28")*(x)*(w^2)+RationalNumber("392")*(x)*(y)-RationalNumber("245")*(x)*(z)+RationalNumber("588")*(x)*(t);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("2")*(u)*(v)*(w)-RationalNumber("6")*(v^2)*(w)-(u)*(w^2)+RationalNumber("13")*(v)*(w^2)-RationalNumber("5")*(w^3)+RationalNumber("14")*(y)*(w)-RationalNumber("28")*(t)*(w);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (u^2)*(w)-RationalNumber("3")*(u)*(v)*(w)+RationalNumber("5")*(u)*(w^2)+RationalNumber("14")*(y)*(w)-RationalNumber("28")*(t)*(w);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("-2")*(z)*(u)*(w)-RationalNumber("2")*(t)*(u)*(w)+RationalNumber("4")*(y)*(v)*(w)+RationalNumber("6")*(z)*(v)*(w)-RationalNumber("2")*(t)*(v)*(w)-RationalNumber("16")*(y)*(w^2)-RationalNumber("10")*(z)*(w^2)+RationalNumber("22")*(t)*(w^2)+RationalNumber("42")*(x)*(w);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("28")*(y)*(u)*(w)+RationalNumber("8")*(z)*(u)*(w)-RationalNumber("20")*(t)*(u)*(w)-RationalNumber("88")*(y)*(v)*(w)-RationalNumber("24")*(z)*(v)*(w)+RationalNumber("68")*(t)*(v)*(w)+RationalNumber("156")*(y)*(w^2)+RationalNumber("40")*(z)*(w^2)-RationalNumber("132")*(t)*(w^2)-RationalNumber("252")*(x)*(w);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("-4")*(y)*(z)*(w)+RationalNumber("10")*(y)*(t)*(w)+RationalNumber("8")*(z)*(t)*(w)-RationalNumber("20")*(t^2)*(w)+RationalNumber("12")*(x)*(u)*(w)-RationalNumber("30")*(x)*(v)*(w)+RationalNumber("15")*(x)*(w^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("-2")*(y^2)*(w)+(y)*(z)*(w)+RationalNumber("2")*(y)*(t)*(w)-RationalNumber("2")*(z)*(t)*(w)+RationalNumber("4")*(t^2)*(w)-RationalNumber("6")*(x)*(u)*(w)+RationalNumber("12")*(x)*(v)*(w)-RationalNumber("6")*(x)*(w^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("8")*(x)*(y)*(w)-RationalNumber("4")*(x)*(z)*(w)+RationalNumber("8")*(x)*(t)*(w);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
}

void getSampleSys2931(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("a");
	SparseMultivariateRationalPolynomial a("a");
	vars.emplace_back("b");
	SparseMultivariateRationalPolynomial b("b");
	vars.emplace_back("c");
	SparseMultivariateRationalPolynomial c("c");
	vars.emplace_back("d");
	SparseMultivariateRationalPolynomial d("d");
	vars.emplace_back("e");
	SparseMultivariateRationalPolynomial e("e");
	vars.emplace_back("t");
	SparseMultivariateRationalPolynomial t("t");
	vars.emplace_back("u");
	SparseMultivariateRationalPolynomial u("u");
	vars.emplace_back("v");
	SparseMultivariateRationalPolynomial v("v");
	vars.emplace_back("w");
	SparseMultivariateRationalPolynomial w("w");
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	SparseMultivariateRationalPolynomial poly;
	poly = (x)*(u)+RN(-1)*(u)*(w)+RN(-1)*(x)*(a);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (y)*(t)+RN(-1)*(t)*(w)+RN(-1)*(y)*(a);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x)*(v)+RN(-1)*(v)*(b)+RN(-1)*(x)*(c);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (z)*(t)+RN(-1)*(t)*(b)+RN(-1)*(z)*(c);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (y)*(v)+RN(-1)*(v)*(d)+RN(-1)*(y)*(e);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (z)*(u)+RN(-1)*(u)*(d)+RN(-1)*(z)*(e);
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys2934(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("a");
	SparseMultivariateRationalPolynomial a("a");
	vars.emplace_back("b");
	SparseMultivariateRationalPolynomial b("b");
	vars.emplace_back("c_1");
	SparseMultivariateRationalPolynomial c_1("c_1");
	vars.emplace_back("c_2");
	SparseMultivariateRationalPolynomial c_2("c_2");
	vars.emplace_back("s_1");
	SparseMultivariateRationalPolynomial s_1("s_1");
	vars.emplace_back("s_2");
	SparseMultivariateRationalPolynomial s_2("s_2");
	SparseMultivariateRationalPolynomial poly;
	poly = (c_1)*(c_2)-(s_1)*(s_2)+(c_1)-(a);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (c_1)*(s_2)+(c_2)*(s_1)+(s_1)-(b);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (c_1^2)+(s_1^2)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (c_2^2)+(s_2^2)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3105(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("p");
	SparseMultivariateRationalPolynomial p("p");
	vars.emplace_back("s");
	SparseMultivariateRationalPolynomial s("s");
	vars.emplace_back("phi");
	SparseMultivariateRationalPolynomial phi("phi");
	SparseMultivariateRationalPolynomial poly;
	poly = RationalNumber("-2")*(p^3)+RationalNumber("2")*(p^3)*(phi^3)-RationalNumber("4")*(phi^3)*(s)*(p^2)+RationalNumber("5")*(phi^3)*(s^3)*(p)-(phi^3)*(s^5);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-2")*(s)*(p^3)-RationalNumber("2")*(phi^3)*(s^2)+(phi^3)*(s^4)-RationalNumber("3")*(phi^3)*(s^2)*(p)+RationalNumber("2")*(phi^3)*(p);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-2")*(s^2)+(s^4)-RationalNumber("4")*(s^2)*(p)+(phi^2)+RationalNumber("1")+RationalNumber("4")*(p);
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3110(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	SparseMultivariateRationalPolynomial poly,poly2;
	poly = (x^2)*(y)*(z);
	poly.setRingVariables(vars);
	poly2 = (x)*(y^2)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(y)*(z^2);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(y);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x^2)*(y^2)*(z);
	poly.setRingVariables(vars);
	poly2 = (x)*(y^2)*(z^2);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x^2)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x);
	poly += poly2;
	poly2 = (z);
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x^2)*(y^2)*(z^2);
	poly.setRingVariables(vars);
	poly2 = (x^2)*(y^2)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(y^2)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (x)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (z)+RationalNumber("1");
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3111(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	SparseMultivariateRationalPolynomial poly;
	poly = RationalNumber("-1")*(x^5)+(y^5)-RationalNumber("3")*(y)-RationalNumber("1");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("5")*(y^4)-RationalNumber("3");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-20")*(x)+(y)-(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3112(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("t");
	SparseMultivariateRationalPolynomial t("t");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	SparseMultivariateRationalPolynomial poly;
	poly = (y^2)*(z)+RationalNumber("2")*(x)*(y)*(t)-RationalNumber("2")*(x)-(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-1")*(x^3)*(z)+RationalNumber("4")*(x)*(y^2)*(z)+RationalNumber("4")*(x^2)*(y)*(t)+RationalNumber("2")*(y^3)*(t)+RationalNumber("4")*(x^2)-RationalNumber("10")*(y^2)+RationalNumber("4")*(x)*(z)-RationalNumber("10")*(y)*(t)+RationalNumber("2");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("2")*(y)*(z)*(t)+(x)*(t^2)-(x)-RationalNumber("2")*(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-1")*(x)*(z^3)+RationalNumber("4")*(y)*(z^2)*(t)+RationalNumber("4")*(x)*(z)*(t^2)+RationalNumber("2")*(y)*(t^3)+RationalNumber("4")*(x)*(z)+RationalNumber("4")*(z^2)-RationalNumber("10")*(y)*(t)-RationalNumber("10")*(t^2)+RationalNumber("2");
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3113(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	vars.emplace_back("t");
	SparseMultivariateRationalPolynomial t("t");
	SparseMultivariateRationalPolynomial poly;
	poly = (y^2)*(z)+RationalNumber("2")*(x)*(y)*(t)-RationalNumber("2")*(x)-(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-1")*(x^3)*(z)+RationalNumber("4")*(x)*(y^2)*(z)+RationalNumber("4")*(x^2)*(y)*(t)+RationalNumber("2")*(y^3)*(t)+RationalNumber("4")*(x^2)-RationalNumber("10")*(y^2)+RationalNumber("4")*(x)*(z)-RationalNumber("10")*(y)*(t)+RationalNumber("2");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("2")*(y)*(z)*(t)+(x)*(t^2)-(x)-RationalNumber("2")*(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-1")*(x)*(z^3)+RationalNumber("4")*(y)*(z^2)*(t)+RationalNumber("4")*(x)*(z)*(t^2)+RationalNumber("2")*(y)*(t^3)+RationalNumber("4")*(x)*(z)+RationalNumber("4")*(z^2)-RationalNumber("10")*(y)*(t)-RationalNumber("10")*(t^2)+RationalNumber("2");
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3117(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("q");
	SparseMultivariateRationalPolynomial q("q");
	vars.emplace_back("c");
	SparseMultivariateRationalPolynomial c("c");
	vars.emplace_back("p");
	SparseMultivariateRationalPolynomial p("p");
	vars.emplace_back("d");
	SparseMultivariateRationalPolynomial d("d");
	vars.emplace_back("b");
	SparseMultivariateRationalPolynomial b("b");
	SparseMultivariateRationalPolynomial poly,poly2;
	poly = (b)-RationalNumber("2");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (c)*(c)*(q)*(q);
	poly.setRingVariables(vars);
	poly2 = -RationalNumber("2")*(d)*(c)*(q)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (d)*(d)*(q)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(c)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(d)*(p)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(b)*(p)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(d)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(d)*(d)*(p)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(b)*(p)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(p)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(d)*(d)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (c)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(d)*(p)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(d)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b)*(d)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(d)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (d)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(b)*(b)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(b)*(d)*(d)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("4")*(b)*(d)*(d)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("4")*(b)*(b)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("4")*(b)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b)*(d)*(d);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(b)*(b);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("4")*(b)+RationalNumber("2");
	poly2.setRingVariables(vars);
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-1")*(b)*(p)*(c)*(q);
	poly.setRingVariables(vars);
	poly2 = -(p)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(d)*(p)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (d)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (p)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -(b)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(b)*(d)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(d)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -(d)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b)*(d);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(d);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-2")*(p)*(q);
	poly.setRingVariables(vars);
	poly2 = -RationalNumber("2")*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -(b)*(b)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2")*(b)*(b)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -(b)*(b)+RationalNumber("2");
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("3")*(c)*(c)*(q)*(q);
	poly.setRingVariables(vars);
	poly2 = -RationalNumber("3")*(d)*(d)*(q)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("4")*(q)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("6")*(c)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("6")*(b)*(d)*(p)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("6")*(b)*(d)*(c)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("6")*(d)*(d)*(p)*(q);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("3")*(c)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("6")*(b)*(d)*(p)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("6")*(b)*(d)*(c);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("3")*(b)*(b)*(d)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("3")*(d)*(d)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b)*(p)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("6")*(b)*(b)*(d)*(d)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2")*(b)*(b)*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("4")*(p);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("3")*(b)*(b)*(d)*(d);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = (b)*(b);
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3120(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	SparseMultivariateRationalPolynomial poly;
	poly = (x^3)+(y)*(x^2)+(x)+(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x^2)*(y)+RationalNumber("3")*(x)+RationalNumber("2");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = (x^2)-(y)*(x)+(z);
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3145(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	SparseMultivariateRationalPolynomial poly,poly2;
	poly = RationalNumber("7")*(y^4)-RationalNumber("20")*(x^2);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("2160")*(x)*(x)*(z^4)+RationalNumber("1512")*(x)*(z^4)+RationalNumber("315")*(z^4)-RationalNumber("4000")*(x)*(x)-RationalNumber("2800")*(x)-RationalNumber("490");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-10080000")*(x^4)*(z^3);
	poly.setRingVariables(vars);
	poly2 = -RationalNumber("28224000")*(x^3)*(z^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("15288000")*(x)*(x)*(z^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("1978032")*(x)*(z^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("180075")*(z^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("23520000")*(x^4)*(y)*(z)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("41395200")*(x^3)*(y)*(z)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("26726560")*(x)*(x)*(y)*(z)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("7727104")*(x)*(y)*(z)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("852355")*(y)*(z)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("40320000")*(x^6)*(y)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("28800000")*(x^5)*(y)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("21168000")*(x^3)*(y)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("4939200")*(x)*(x)*(y)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("347508")*(x)*(y)*(y)*(z);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("67200000")*(x^5)*(y^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("94080000")*(x^4)*(y^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("40924800")*(x^3)*(y^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = RationalNumber("2634240")*(x)*(x)*(y^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("2300844")*(x)*(y^3);
	poly2.setRingVariables(vars);
	poly += poly2;
	poly2 = -RationalNumber("432180")*(y^3);
	poly += poly2;
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void getSampleSys3149(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars) {
	polys.clear();
	vars.clear();
	vars.emplace_back("x");
	SparseMultivariateRationalPolynomial x("x");
	vars.emplace_back("y");
	SparseMultivariateRationalPolynomial y("y");
	vars.emplace_back("t");
	SparseMultivariateRationalPolynomial t("t");
	vars.emplace_back("z");
	SparseMultivariateRationalPolynomial z("z");
	vars.emplace_back("u");
	SparseMultivariateRationalPolynomial u("u");
	vars.emplace_back("v");
	SparseMultivariateRationalPolynomial v("v");
	SparseMultivariateRationalPolynomial poly;
	poly = RationalNumber("45")*(y)+RationalNumber("35")*(u)-RationalNumber("165")*(v)-RationalNumber("36");
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("35")*(y)+RationalNumber("25")*(z)+RationalNumber("40")*(t)-RationalNumber("27")*(u);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("25")*(y)*(u)-RationalNumber("165")*(v^2)+RationalNumber("15")*(x)-RationalNumber("18")*(z)+RationalNumber("30")*(t);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("15")*(y)*(z)+RationalNumber("20")*(t)*(u)-RationalNumber("9")*(x);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-11")*(v^3)+(x)*(y)+RationalNumber("2")*(z)*(t);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("-11")*(u)*(v)+RationalNumber("3")*(v^2)+RationalNumber("99")*(x);
	poly.setRingVariables(vars);
	polys.push_back(poly);
	poly = RationalNumber("10000")*(v^3)+RationalNumber("6600")*(v)+RationalNumber("2673");
	poly.setRingVariables(vars);
	polys.push_back(poly);
}

void Sys2915Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys2915(polys,vars);
	testSys(polys,vars,"2915","Hairer-1",showOutput,isLazard);
}

void Sys2920Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys2920(polys,vars);
	testSys(polys,vars,"2920","Butcher",showOutput,isLazard);
}

void Sys2926Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys2926(polys,vars);
	testSys(polys,vars,"2926","Gerdt",showOutput,isLazard);
}

void Sys2931Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys2931(polys,vars);
	testSys(polys,vars,"2931","Pappus",showOutput,isLazard);
}

void Sys2934Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys2934(polys,vars);
	testSys(polys,vars,"2934","robotPlanoEasy",showOutput,isLazard);
}

void Sys3105Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3105(polys,vars);
	testSys(polys,vars,"3105","4-body-homog",showOutput,isLazard);
}

void Sys3110Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3110(polys,vars);
	testSys(polys,vars,"3110","Arnborg-Lazard",showOutput,isLazard);
}

void Sys3111Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3111(polys,vars);
	testSys(polys,vars,"3111","Barry",showOutput,isLazard);
}

void Sys3112Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3112(polys,vars);
	testSys(polys,vars,"3112","Caprasse-Li",showOutput,isLazard);
}

void Sys3113Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3113(polys,vars);
	testSys(polys,vars,"3113","Caprasse",showOutput,isLazard);
}

void Sys3117Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3117(polys,vars);
	testSys(polys,vars,"3117","Czapor-Geddes-Wang",showOutput,isLazard);
}

void Sys3120Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3120(polys,vars);
	testSys(polys,vars,"3120","GonzalezGonzalez",showOutput,isLazard);
}

void Sys3145Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3145(polys,vars);
	testSys(polys,vars,"3145","Rose",showOutput,isLazard);
}

void Sys3149Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys;
	std::vector<Symbol> vars;
	getSampleSys3149(polys,vars);
	testSys(polys,vars,"3149","Trinks-2",showOutput,isLazard);
}

//	for (auto t : polys)
//		cerr << t << endl;
//	for (auto t : vars)
//		cerr << t << " ";
//	cerr << endl;

// issac2005 examples //

//  pnum := 850088191; modpnp := 962592769;
void Sys2887Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("z2+z3+z4+z5"),
  		SMQP("33461*c3-94642"),
 		SMQP("c2+y2+y3+y4+y5"),
 		SMQP("9*c2-8"),
 		SMQP("36926*c1-208885"),
 		SMQP("-3*c3*y2+3*x3+3*x4+8"),
 		SMQP("-c1*y2+9*y2^2+z2"),
 		SMQP("z5^2+y5^2-c2"),
 		SMQP("3*z4*z5+3*y4*y5+x4-1"),
 		SMQP("z4^2+y4^2+x4^2-1"),
 		SMQP("3*z3*z4+3*y3*y4+3*x3*x4-1"),
 		SMQP("z3^2+y3^2+x3^2-1"),
 		SMQP("-3*c3*y2*x3+3*z2*z3+3*y2*y3+3*x3-1")};
	std::vector<Symbol> vars = {Symbol("c1"), Symbol("c2"), Symbol("c3"), Symbol("z2"), Symbol("z3"), Symbol("z4"),
		Symbol("z5"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("x3"), Symbol("x4")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2887","Chemkin",showOutput,isLazard);
}

void Sys2888Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-y^2*z^2-y^2+24*y*z-z^2-13"),
		SMQP("-x^2*z^2-x^2+24*x*z-z^2-13"),
		SMQP("-x^2*y^2-x^2+24*x*y-y^2-13")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2888","Cyclohexane",showOutput,isLazard);
}

// pnum := 358079; modpnp := 962592769;
void Sys2889Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("210*z-210"),
		SMQP("-80*d+180*y+855*z"),
		SMQP("136*d*z-114*c+152*x+720*y"),
		SMQP("112*d*y+105*c*z-144*b+126*w+595*x"),
		SMQP("90*d*x+84*c*y+78*b*z-170*a+102*v+480*w"),
		SMQP("70*d*w+65*c*x+60*b*y+55*a*z+80*u+375*v"),
		SMQP("52*d*v+48*c*w+44*b*x+40*a*y+280*u"),
		SMQP("36*d*u+33*c*v+20*b*w+27*a*x"),
		SMQP("20*c*u+18*b*v+16*a*w"),
		SMQP("8*b*u+7*a*v")};
	std::vector<Symbol> vars = {'a', 'b', 'c', 'd', 'u', 'v', 'w', 'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2889","Dessin-2",showOutput,isLazard);
}

// pnum := 105761;
void Sys2890aTest(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1 + x1*x2 + x2*x3)*x4 - 1"),
		SMQP("(x2 + x1*x3)*x4 - 2"),
		SMQP("x3*x4 - 3"),
		SMQP("x1 + x2 + x3 + 1")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2890a","Eco4",showOutput,isLazard);
}

// pnum := 105761;
void Sys2890bTest(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1 + x1*x2 + x2*x3 + x3*x4)*x5 - 1"),
		SMQP("(x2 + x1*x3 + x2*x4)*x5 - 2"),
		SMQP("(x3 + x1*x4)*x5 - 3"),
		SMQP("x4*x5 - 4"),
		SMQP("x1 + x2 + x3 + x4 + 1")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2890b","Eco5",showOutput,isLazard);
}

// pnum := 105761; // can specialize x1 = 1
// use PolynomialIdeal package in Maple
void Sys2890Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5)*x6 - 1"),
		SMQP("(x2 + x1*x3 + x2*x4 + x3*x5)*x6 - 2"),
		SMQP("(x3 + x1*x4 + x2*x5)*x6 - 3"),
		SMQP("(x4 + x1*x5)*x6 - 4"),
		SMQP("x5*x6 - 5"),
		SMQP("x1 + x2 + x3 + x4 + x5 + 1")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2890","Eco6",showOutput,isLazard);
}

// pnum := 387799; modpnp := 962592769;
void Sys2891Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6)*x7 - 1"),
		SMQP("(x2 + x1*x3 + x2*x4 + x3*x5 + x4*x6)*x7 - 2"),
		SMQP("(x3 + x1*x4 + x2*x5 + x3*x6)*x7 - 3"),
		SMQP("(x4 + x1*x5 + x2*x6)*x7 - 4"),
		SMQP("(x5 + x1*x6)*x7 - 5"),
		SMQP("x6*x7 - 6"),
		SMQP("x1 + x2 + x3 + x4 + x5 + x6 + 1")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2891","Eco7",showOutput,isLazard);
}

// pnum := 2671;
void Sys2892Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("30752*x^3 + 46128*y*x^2 - 216256*x + 46128*x*y^2 + 141980 + 30752*y^3- 216256*y"),
		SMQP("30752*x^3 + 46128*z*x^2 - 216256*x + 46128*x*z^2 + 141980 + 30752*z^3+ 216256*z"),
		SMQP("46128*y^3 + 46128*x*y^2 + 46128*z*y^2 + 46128*y*x^2 - 432512*y +46128*x*z*y + 46128*z^2*y + 46128*x^3 - 432512*x + 425940 +46128*z*x^2 - 432512*z + 46128*x*z^2 + 46128*z^3")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2892","Fabfau",showOutput,isLazard);
}

// pnum := 24499; modpnp := 962592769;
void Sys2893Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-2*q*p-2*p^2-2*q+8*p-2"),
		SMQP("-3*q*c*p+2*q*p*d+4*p^2*d+3*c*p+q*d-7*p*d"),
		SMQP("q^2*c^2-2*q^2*c*d-2*q*c*p*d+q^2*d^2+2*q*p*d^2+p^2*d^2-2*q*c^2+4*q*c*p+2*q*c*d+2*c*p*d-4*q*d^2-4*p*d^2+c^2+2*q*p+10*p^2-4*c*d+4*d^2-2*q-8*p+2"),
		SMQP("3*q^2*c^2+12*q*c*p*d-3*q^2*d^2+6*q*p*d^2-3*p^2*d^2-6*q*c^2+12*q*c*d+12*c*p*d-4*q^2+3*c^2+5*p^2-12*c*d+12*d^2-6*p+5")};
	std::vector<Symbol> vars = {'q', 'c', 'p', 'd'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2893","Fee1",showOutput,isLazard);
}

// pnum := 159223;
void Sys2894Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-3779*p2*p3*p4*p5+3708*p2*p3*p4+2425*p2*p3*p5+422*p2*p4*p5+3354*p3*p4*p5-1226*p2*p3-308*p2*p4-3093*p3*p4-719*p2*p5-2256*p3*p5-873*p4*p5-593*p2+382*p3+652*p4+897*p5+455"),
		SMQP("-98*p1*p3*p4*p5-3376*p1*p3*p4+421*p1*p3*p5+942*p1*p4*p5+778*p3*p4*p5+1209*p1*p3+899*p1*p4+1510*p3*p4-1463*p1*p5+251*p3*p5+25*p4*p5+537*p1-868*p3-697*p4-27*p5+486"),
		SMQP("1480*p1*p2*p4*p5-61*p1*p2*p4-1228*p1*p2*p5-2337*p1*p4*p5-890*p2*p4*p5+1476*p1*p2+389*p1*p4+2161*p2*p4+637*p1*p5+245*p2*p5+622*p4*p5-762*p1-848*p2-482*p4-204*p5-210"),
		SMQP("2725*p1*p2*p3*p5-2522*p1*p2*p3-2257*p1*p2*p5-1939*p1*p3*p5-1264*p2*p3*p5+1987*p1*p2-1667*p1*p3+791*p2*p3+1428*p1*p5+384*p2*p5+1110*p3*p5+183*p1-951*p2+777*p3-235*p5-152"),
		SMQP("258*p1*p2*p3*p4-20*p1*p2*p3+2387*p1*p2*p4+157*p1*p3*p4+83*p2*p3*p4-597*p1*p2-182*p1*p3+62*p2*p3-1255*p1*p4-1012*p2*p4-299*p3*p4-646*p1+600*p2+147*p3+28*p4+88")};
	std::vector<Symbol> vars = {Symbol("p1"), Symbol("p2"), Symbol("p3"), Symbol("p4"), Symbol("p5")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2894","Gametwo5",showOutput,isLazard);
}

// pnum := 116663; modpnp := 962592769;
void Sys2895Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-10*x1*x6^2+ 2*x2*x6^2-x3*x6^2+x4*x6^2+ 3*x5*x6^2+x1*x6+ 2*x2*x6+x3*x6+ 2*x4*x6+x5*x6+ 10*x1+ 2*x2-x3+ 2*x4-2*x5"),
		SMQP("2*x1*x6^2-11*x2*x6^2+ 2*x3*x6^2-2*x4*x6^2+x5*x6^2+ 2*x1*x6+x2*x6+ 2*x3*x6+x4*x6+ 3*x5*x6+ 2*x1+ 9*x2+ 3*x3-x4-2*x5"),
		SMQP("-x1*x6^2+ 2*x2*x6^2-12*x3*x6^2-x4*x6^2+x5*x6^2+x1*x6+ 2*x2*x6-2*x4*x6-2*x5*x6-x1+ 3*x2+ 10*x3+ 2*x4-x5"),
		SMQP("x1*x6^2-2*x2*x6^2-x3*x6^2-10*x4*x6^2+ 2*x5*x6^2+ 2*x1*x6+x2*x6-2*x3*x6+ 2*x4*x6+ 3*x5*x6+ 2*x1-x2+ 2*x3+ 12*x4+x5"),
		SMQP("3*x1*x6^2+x2*x6^2+x3*x6^2+ 2*x4*x6^2-11*x5*x6^2+x1*x6+ 3*x2*x6-2*x3*x6+ 3*x4*x6+ 3*x5*x6-2*x1-2*x2-x3+x4+ 10*x5"),
		SMQP("x1+x2+x3+x4+x5-1")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2895","Geneig",showOutput,isLazard);
}

// pnum := 1549;
void Sys2896Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("4*a^2 + 3*d^2 + 9*d + 9*b^2 + 7*c^2 + 7*a - 5*d*a + 2*a*b + 5*b - 3*d*b + 6*b*c + 7*c - 6*d*c - 2*a*c + 5"),
		SMQP("8*a^2 -2*d^2 + 9*a*d + 9*b*d +  6*b^2  - 4*d + 8*a + 9*a*b + 4*b + 8*c - 7*d*c - 3*a*c - 7*b*c - 6*c^2 + 2"),
		SMQP("5*a^2 + 8*d^2 + 5*a*d + 2*c*d + 3*d + 7*b^2 + 7*c^2 - 7*a + 2*a*b - 7*b - 4*d*b - 8*c - 7*a*c - 8*b*c + 8"),
		SMQP("2*a^2 + 7*d^2 + 5*a*d + 3*b*d  - 5*d + 4*a + 9*a*b + 6*b  - 4*b^2 - 9*c - 5*d*c - 7*a*c - 5*b*c - 4*c^2 + 2")};
	std::vector<Symbol> vars = {'a', 'b', 'c', 'd'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2896","ISSAC97",showOutput,isLazard);
}

// pnum := 450367; modpnp := 962592769;
void Sys2897Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x10^2-2*x5"),
		SMQP("-x6*x8+x7*x9-210*x6"),
		SMQP("-x4*x8+700000"),
		SMQP("x4*x8-x7*x9-x9*x10-410*x9"),
		SMQP("-x3*x8+x6*x8+x9*x10+210*x6+1300000"),
		SMQP("10*x2*x8+10*x3*x8+10*x6*x8+10*x7*x9-10*x2*x10-11*x7*x10-10*x9*x10-10*x10^2+1400*x6-4200*x10"),
		SMQP("-10*x1*x8-10*x2*x8-10*x3*x8-10*x4*x8-10*x6*x8+10*x2*x10+11*x7*x10+20000*x2+14*x5"),
		SMQP("-64*x2*x7-10*x7*x9-11*x7*x10+320000*x1-16*x7+7000000"),
		SMQP("-32*x2*x7-5*x2*x8-5*x2*x10+160000*x1-5000*x2"),
		SMQP("64*x2*x7-10*x1*x8+10*x7*x9+11*x7*x10-320000*x1")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2897","Methan61",showOutput,isLazard);
}

// pnum := 55313;
void Sys2898Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*x^2-2*y^2+2*z^2-2*t^2-1"),
		SMQP("2*x^3-2*y^3+2*z^3-2*t^3-1"),
		SMQP("2*x^4-2*y^4+2*z^4-2*t^4-1"),
		SMQP("2*x^5-2*y^5+2*z^5-2*t^5-1")};
	std::vector<Symbol> vars = {'x', 'y', 'z', 't'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2898","Reimer-4",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys2899Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2+x*y+y^2-2*x*z-4*y*z+3*z^2-3*x*t+2*y*t+t^2-3*x-2*y+3*z-2*t-2"),
		SMQP("2*x^2-x*y+y^2-x*z-y*z-6*z^2-x*t+y*t-5*z*t-3*t^2-5*x+y+5*z+2*t+5"),
		SMQP("-3-3*x*y+2*x*z+x*t^2-5*x*z^2-5*z^2*t-3*x*t-2*z*t+x*y*z+x*y*t-x^2*z+x^2-y^2+2*z^2+11*z-2*t-x+y+x^3+y^3-3*z^3+2*t^3-3*t^2-5*y^2*z+7*y*z^2"),
		SMQP("-15+2*x*y+11*x*t^2+5*x*z^2-z*t-4*x*y*z+6*x*y*t-x^2*z+3*x^2+2*y^2-z^2+4*z-10*t-35*x-14*y-x^3+6*y^3+15*z^3+4*t^3+5*t^2+6*y^2*z+4*y*z^2-x*z*t+6*x^2*y-12*x*y^2-7*y^2*t+2*y*t")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2899","Uteshev-Bikker",showOutput,isLazard);
}

// pnum := 7433;
void Sys2900Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("y^4+x*y^2*z+x^2-2*x*y+y^2+z^2"),
		SMQP("x*y^4+y*z^4-2*x^2*y-3"),
		SMQP("-x^3*y^2+x*y*z^3+y^4+x*y^2*z-2*x*y")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"2900","Weispfenning-94",showOutput,isLazard);
}

// ASCM09-zerodim examples //

// pnum := 7841; modpnp := 962592769;
void Sys3107Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-6*p^3+4*p^3*phi^3+15*phi^3*s^3*p-3*phi^3*s^5-12*phi^3*s*p^2-3*phi^3*s*p+phi^3*s^3"),
		SMQP("-9*phi^3*s^2*p-5*phi^3*s^2-6*s*p^3+3*phi^3*s^4+5*phi^3*p"),
		SMQP("-12*s^2*p-6*s^2+3*s^4+4*phi^2+3+12*p")};
	std::vector<Symbol> vars = {'p', 's', Symbol("phi")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3107","5-body-homog",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3109Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2*y*z + x*y^2*z + x*y*z^2 + x*y*z + x*y + x*z + y*z"),
		SMQP("x^2*y^2*z + x*y^2*z^2 + x^2*y*z + x*y*z + y*z + x + z"),
		SMQP("x^2*y^2*z ^2 + x^2*y^2*z + x*y^2*z + x*y*z + x*z + z + 1")};
	std::vector<Symbol> vars = {'z', 'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3109","Arnborg-Lazard-rev",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3114Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2 - 7*x1 + x1^2*x2 - 1/2*(x3 - x1 )"),
		SMQP("6*x1 - x1^2*x2 - 5*(x4 - x2 )"),
		SMQP("2 - 7*x3 + x3^2*x4 - 1/2*(x1 - x3 )"),
		SMQP("6*x3 - x3^2*x4 + 1 + 1/2*(x2 - x4 )")};
	std::vector<Symbol> vars = {Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3114","Chemical-reaction",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3115Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6)*(x-7)*(x-8)*(x-9)*(x-10)+1/100"),
		SMQP("((x-1)^2+(y-1)^2-2)*((x-2)^2+(y-2)^2-2)*((x-3)^2+(y-3)^2-2)*((x-4)^2+(y-4)^2-2)*((x-5)^2+(y-5)^2-2)+1/1000")};
	std::vector<Symbol> vars = {'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3115","Circles",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3116Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a + b + c + d + e"),
		SMQP("a*b + b*c + c*d + d*e + e*a"),
		SMQP("a*b*c + b*c*d + c*d*e + d*e*a + e*a*b"),
		SMQP("a*b*c*d + b*c*d*e + c*d*e*a + d*e*a*b  + e*a*b*c"),
		SMQP("a*b*c*d*e - 1")};
	std::vector<Symbol> vars = {'a', 'b', 'c', 'd', 'e'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3116","Cyclic-5",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3118Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("30752*x^3 + 46128*y*x^2 - 216256*x + 46128*x*y^2 + 141980 + 30752*y^3- 216256*y"),
		SMQP("30752*x^3 + 46128*z*x^2 - 216256*x + 46128*x*z^2 + 141980 + 30752*z^3+ 216256*z"),
		SMQP("46128*y^3 + 46128*x*y^2 + 46128*z*y^2 + 46128*y*x^2 - 432512*y +46128*x*z*y + 46128*z^2*y + 46128*x^3 - 432512*x + 425940 +46128*z*x^2 - 432512*z + 46128*x*z^2 + 46128*z^3")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3118","Fabfaux",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3119Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("1/100 - 4*s*(s - 1)*(s - b)*(s - c)"),
		SMQP("1/5 - b*c"),
		SMQP("2*s - 1 - b - c")};
	std::vector<Symbol> vars = {'c', 'b', 's'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3119","Geometric-constraints",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3121Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*x^2 + 2*y^2 + 2*z^2 + 2*t^2 + u^2 - u"),
		SMQP("2*x*y + 2*y*z + 2*z*t + 2*t*u - t"),
		SMQP("2*x*z + 2*y*t + t^2 + 2*z*u - z"),
		SMQP("2*x*t + 2*z*t + 2*y*u - y"),
		SMQP("2*x + 2*y + 2*z + 2*t + u - 1")};
	std::vector<Symbol> vars = {'x', 'y', 'z', 't', 'u'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3121","Katsura-4",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3122Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*x^7+2-3*x^2-x^4"),
		SMQP("x^2+2*y^2-5"),
		SMQP("x*z-1")};
	std::vector<Symbol> vars = {'z', 'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3122","LHLP1",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3123Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("y*(320+1600*y^4-240*y^5-471*y^6+36*y^7-48*y^2+36*y^8)"),
		SMQP("-40*y^2+3*y^3+6*y^4+8*x"),
		SMQP("-8*y*z-8-40*y^4+3*y^5+6*y^6")};
	std::vector<Symbol> vars = {'z', 'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3123","LHLP2",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3124Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2450*x^6-1241*x^4+196*x^2-49"),
		SMQP("86*x^2+35*y*x^2-77*y-14"),
		SMQP("175*x^4+70*x^2*z-76*x^2-154*z+21")};
	std::vector<Symbol> vars = {'z', 'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3124","LHLP3",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3125Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^4+y^4-1"),
		SMQP("x^5*y^2-4*x^3*y^3+x^2*y^5-1")};
	std::vector<Symbol> vars = {'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3125","LHLP4",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3126Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-7*x*y*z+6*y*z-14*x*z+9*z-3*x*y-12*y-x+1"),
		SMQP("2*x*y*z-y*z+14*z+15*x*y+14*y-15*x"),
		SMQP("-8*x*y*z+11*y*z-12*x*z-5*z+15*x*y+2*y+10*x-14")};
	std::vector<Symbol> vars = {'z', 'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3126","LHLP5",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3127Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*x1^2-x1*y1-2*x2^2+x2*y2+4"),
		SMQP("4*x1*x2-x1*y2-x2*y1"),
		SMQP("x1*y1-x2*y2-2*y1^2+2*y2^2+4"),
		SMQP("x1*y2+x2*y1-4*y1*y2")};
	std::vector<Symbol> vars = {Symbol("x2"), Symbol("y2"), Symbol("x1"), Symbol("y1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3127","LHLP6",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3128Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("1 - c*x - x*y^2 - x*z^2"),
		SMQP("1 - c*y - y*x^2 - y*z^2"),
		SMQP("1 - c*z - z*x^2 - z*y^2"),
		SMQP("8*c^6 + 378*c^3 - 27")};
	std::vector<Symbol> vars = {'x', 'y', 'c', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3128","Neural-network",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3129Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^10 + y + z - 1"),
		SMQP("x + y^10 + z - 1"),
		SMQP("x + y + z^10 - 1")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3129","NLD-10-3",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3130Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^3 + y + z + t- 1"),
		SMQP("x + y^3 + z + t -1"),
		SMQP("x + y + z^3 + t-1"),
		SMQP("x + y + z + t^3 -1")};
	std::vector<Symbol> vars = {'x', 'y', 'z', 't'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3130","NLD-3-4",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3131Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^3 + y + z + t + u- 1"),
		SMQP("x + y^3 + z + t + u-1"),
		SMQP("x + y + z^3 + t + u-1"),
		SMQP("x + y + z + t^3 + u-1"),
		SMQP("x + y + z + t + u^3 -1")};
	std::vector<Symbol> vars = {'x', 'y', 'z', 't', 'u'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3131","NLD-3-5",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3132Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^4 + y + z + t + u- 1"),
		SMQP("x + y^4 + z + t + u-1"),
		SMQP("x + y + z^4 + t + u-1"),
		SMQP("x + y + z + t^4 + u-1"),
		SMQP("x + y + z + t + u^4 -1")};
	std::vector<Symbol> vars = {'x', 'y', 'z', 't', 'u'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3132","NLD-4-5",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3133Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^7+y+z-1"),
		SMQP("y^7+x+z-1"),
		SMQP("z^7+x+y-1")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3133","NLD-7-3",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3134Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^8+y+z-1"),
		SMQP("y^8+x+z-1"),
		SMQP("z^8+x+y-1")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3134","NLD-8-3",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3135Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^9+y+z-1"),
		SMQP("y^9+x+z-1"),
		SMQP("z^9+x+y-1")};
	std::vector<Symbol> vars = {'x', 'y', 'z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3135","NLD-9-3",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3136Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x10^2+x10-x9"),
		SMQP(" x9^2-x8+x9"),
		SMQP("x8^2-x7+x8"),
		SMQP("x7^2-x6+x7"),
		SMQP("x6^2-x5+x6"),
		SMQP("x5^2-x4+x5"),
		SMQP("x4^2-x3+x4"),
		SMQP("x3^2-x2+x3"),
		SMQP("x2^2-x1+x2"),
		SMQP("x1^2-2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3136","NQL-10-2",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3137Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x10^4+x10^2-x9"),
		SMQP("x9^4+x9^2-x8"),
		SMQP("x8^4+x8^2-x7"),
		SMQP("x7^4+x7^2-x6"),
		SMQP("x6^4+x6^2-x5"),
		SMQP("x5^4+x5^2-x4"),
		SMQP("x4^4+x4^2-x3"),
		SMQP("x3^4+x3^2-x2"),
		SMQP("x2^4+x2^2-x1"),
		SMQP("x1^4-2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3137","NQL-10-4",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3138Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x15^2-x14+x15"),
		SMQP("x14^2-x13+x14"),
		SMQP("x13^2-x12+x13"),
		SMQP("x12^2-x11+x12"),
		SMQP("x11^2-x10+x11"),
		SMQP("x10^2+x10-x9"),
		SMQP("x9^2-x8+x9"),
		SMQP("x8^2-x7+x8"),
		SMQP("x7^2-x6+x7"),
		SMQP("x6^2-x5+x6"),
		SMQP("x5^2-x4+x5"),
		SMQP("x4^2-x3+x4"),
		SMQP("x3^2-x2+x3"),
		SMQP("x2^2-x1+x2"),
		SMQP("x1^2-2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10"), Symbol("x11"), Symbol("x12"), Symbol("x13"), Symbol("x14"), Symbol("x15")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3138","NQL-15-2",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3139Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x5^4+x5^2-x4"),
		SMQP("x4^4+x4^2-x3"),
		SMQP("x3^4+x3^2-x2"),
		SMQP("x2^4+x2^2-x1"),
		SMQP("x1^4-2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3139","NQL-5-4",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3140Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2 + y^2 - x*y - 1"),
		SMQP("y^2 + z^2 - y*z - a^2"),
		SMQP("z^2 + x^2 - z*x - b^2"),
		SMQP("a^2 - 1 + b - b^2"),
		SMQP("3*b^6 + 56*b^4 - 122*b^3 + 56*b^2 + 3")};
	std::vector<Symbol> vars = {'z', 'y', 'x', 'a', 'b'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3140","P3P-special",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3141Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1+5)^2+ (y1-0)^2 - 1"),
		SMQP("(x2-5)^2+ (y2-0)^2 - 1"),
		SMQP("(x3)^2  + (y3-5)^2 - 1"),
		SMQP("(x1-x2)^2 + (y1-y2)^2 - 3"),
		SMQP("(x1-x3)^2 + (y1-y3)^2 - 3"),
		SMQP("(x2-x3)^2 + (y2-y3)^2 - 3")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("y1"), Symbol("x2"), Symbol("y2"), Symbol("x3"), Symbol("y3")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3141","PlateForme2d-easy",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3142Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a^2 + a"),
		SMQP("a*b + b + a*b^2 +a"),
		SMQP("b^2*c + c + b*c^3 + b"),
		SMQP("c^3*d + d + c*d^4 + c"),
		SMQP("d^4*e + e + d*e^5 + d")};
	std::vector<Symbol> vars = {'a', 'b', 'c', 'd', 'e'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3142","R-5",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3143Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a^2 + a"),
		SMQP("a*b + b + a*b^2 +a"),
		SMQP("b^2*c + c + b*c^3 + b"),
		SMQP("c^3*d + d + c*d^4 + c"),
		SMQP("d^4*e + e + d*e^5 + d"),
		SMQP("e^5*f + f + e*f^6 + e")};
	std::vector<Symbol> vars = {'a', 'b', 'c', 'd', 'e', 'f'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3143","R-6",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3144Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x4*x13+x5*x14+x6*(1-x13-x14)"),
		SMQP("x4*x15+x5*x16-x6*(x15+x16)"),
		SMQP("x7*x13+x8*x14+x9*(1-x13-x14)"),
		SMQP("x7*x15+x8*x16-x9*(x15+x16)-1"),
		SMQP("x10*x13+x11*x14+x12*(1-x13-x14)"),
		SMQP("x10*x15+x11*x16-x12*(x15+x16)"),
		SMQP("x1*x13+x2*x14+x3*(1-x13-x14)"),
		SMQP("x1*x15+x2*x16-x3*(x15+x16)"),
		SMQP("x1*x4*x13+x2*x5*x14+x3*x6*(1-x13-x14)-1"),
		SMQP("x1*x4*x15+x2*x5*x16-x3*x6*(x15+x16)"),
		SMQP("x1*x7*x13+x2*x8*x14+x3*x9*(1-x13-x14)"),
		SMQP("x1*x7*x15+x2*x8*x16-x3*x9*(x15+x16)"),
		SMQP("x1*x10*x13+x2*x11*x14+x3*x12*(1-x13-x14)"),
		SMQP("x1*x10*x15+x2*x11*x16-x3*x12*(x15+x16)-1")};
	std::vector<Symbol> vars = {Symbol("x16"), Symbol("x15"), Symbol("x14"), Symbol("x13"), Symbol("x12"), Symbol("x11"), Symbol("x10"), Symbol("x9"), Symbol("x8"), Symbol("x7"), Symbol("x6"), Symbol("x5"), Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3144","Reif",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3146Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x20^30-x19"),
		SMQP("x19^30-x18"),
		SMQP("x18^30-x17"),
		SMQP("x17^30-x16"),
		SMQP("x16^30-x15"),
		SMQP("x15^30-x14"),
		SMQP("x14^30-x13"),
		SMQP("x13^30-x12"),
		SMQP("x12^30-x11"),
		SMQP("x11^30-x10"),
		SMQP("x10^30-x9"),
		SMQP("x9^30-x8"),
		SMQP("x8^30-x7"),
		SMQP("x7^30-x6"),
		SMQP("x6^30-x5"),
		SMQP("x5^30-x4"),
		SMQP("x4^30-x3"),
		SMQP("x3^30-x2"),
		SMQP("x2^30-x1"),
		SMQP("x1^30-2")};
	std::vector<Symbol> vars = {Symbol("x20"), Symbol("x19"), Symbol("x18"), Symbol("x17"), Symbol("x16"), Symbol("x15"), Symbol("x14"), Symbol("x13"), Symbol("x12"), Symbol("x11"), Symbol("x10"), Symbol("x9"), Symbol("x8"), Symbol("x7"), Symbol("x6"), Symbol("x5"), Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3146","Simple-NQL-20-30",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3147Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*x1*(2 - x1 - y1 ) + x2 - x1"),
		SMQP("2*x2*(2 - x2 - y2 ) + x1 - x2"),
		SMQP("2*y1*(5 - x1 - 2*y1 ) + y2 - y1"),
		SMQP("y2*(3 - 2*x2 - 4*y2 ) + y1 - y2")};
	std::vector<Symbol> vars = {Symbol("y2"), Symbol("y1"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3147","Takeuchi-Lu",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3148Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("10000*d - 323536"),
		SMQP("(x-48)^2+(y-89)^2+(28)^2-(48*a+89*b+28*c)^2-d^2"),
		SMQP("(x-77)^2+(y-3)^2+(37)^2-(77*a+3*b+37*c)^2-d^2"),
		SMQP("(x-49)^2+(y-23)^2+(57)^2-(49*a+23*b+57*c)^2-d^2"),
		SMQP("a*x+b*y"),
		SMQP("a^2+b^2+c^2-1")};
	std::vector<Symbol> vars = {'x', 'y', 'a', 'b', 'c', 'd'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3148","Themos-net",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3150Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("45*y+35*u-165*v-36"),
		SMQP("35*y+25*z+40*t-27*u"),
		SMQP("25*y*u-165*v^2+15*x-18*z+30*t"),
		SMQP("15*y*z+20*t*u-9*x"),
		SMQP("-11*v^3+x*y+2*z*t"),
		SMQP("-11*u*v+3*v^2+99*x")};
	std::vector<Symbol> vars = {'x', 'y', 't', 'z', 'u', 'v'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3150","Trinks-difficult",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3151Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2+x*y+y^2-2*x*z-4*y*z+3*z^2-3*x*t+2*y*t+t^2-3*x-2*y+3*z-2*t-2"),
		SMQP("2*x^2-x*y+y^2-x*z-y*z-6*z^2-x*t+y*t-5*z*t-3*t^2-5*x+y+5*z+2*t+5"),
		SMQP("-3-3*x*y+2*x*z+x*t^2-5*x*z^2-5*z^2*t-3*x*t-2*z*t+x*y*z+x*y*t-x^2*z+x^2-y^2+2*z^2+11*z-2*t-x+y+x^3+y^3-3*z^3+2*t^3-3*t^2-5*y^2*z+7*y*z^2"),
		SMQP("-15+2*x*y+11*x*t^2+5*x*z^2-z*t-4*x*y*z+6*x*y*t-x^2*z+3*x^2+2*y^2-z^2+4*z-10*t-35*x-14*y-x^3+6*y^3+15*z^3+4*t^3+5*t^2+6*y^2*z+4*y*z^2-x*z*t+6*x^2*y-12*x*y^2-7*y^2*t+2*y*t")};
	std::vector<Symbol> vars = {'x', 'y', 'z', 't'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3151","Uteshev-Bikker",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3152Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6)*(x-7)*(x-8)*(x-9)*(x-10)*(x-11)*(x-12)*(x-13)*(x-14)*(x-15)*(x-16)*(x-17)*(x-18)*(x-19)*(x-20)+1/100")};
	std::vector<Symbol> vars = {'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3152","Wilkinson20",showOutput,isLazard);
}

// pnum := 7841; modpnp := 962592769;
void Sys3153Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)+1/100"),
		SMQP("(y-x-1)*(y-x-2)*(y-x-3)*(y-x-4)*(y-x-5)")};
	std::vector<Symbol> vars = {'y', 'x'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3153","Wilkinsonxy",showOutput,isLazard);
}

void Sys2250Test(bool showOutput, bool isLazard) {
    std::vector<SparseMultivariateRationalPolynomial> polys = {
        SMQP("x_4*(-3+x_1)"),
        SMQP("2*x_4, -2*x_1*x_4"),
        SMQP("2*x_4*(1+x_1)"),
        SMQP("-x_1*x_3*(-1+x_1)"),
        SMQP("-x_2-2*x_1*x_5+x_1*x_2"),
        SMQP("-2*x_3+x_1*x_3-x_1^2*x_2+x_1*x_2"),
        SMQP("x_3+x_1^2*x_5-2*x_1*x_2+x_1*x_5"),
        SMQP("x_1^2*x_4+x_1*x_4+x_1*x_5+x_2")
    };
    std::vector<Symbol> vars = {Symbol("x_1"), Symbol("x_2"), Symbol("x_3"), Symbol("x_4"), Symbol("x_5")};
    for (size_t i = 0; i < polys.size(); ++i)
        polys[i].setRingVariables(vars);
    testSys(polys,vars,"2250","2250",showOutput,isLazard);
}

