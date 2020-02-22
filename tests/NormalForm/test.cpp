#include <bpas.h>
#include <vector>
#include "../MapleTestTool/MapleTestTool.hpp"


using namespace std;

///////////////////////////////// Prototypes ////////////////////////////////////
void testNormalForm();
void testMDPD();

///////////////////////////////////// MAIN ////////////////////////////////////////

int main(int argc, char** argv) {

    testNormalForm();
    testMDPD();
    return 0;
}


///////////////////////////////////// Func ////////////////////////////////////////

void testNormalForm() {

    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	SparseMultivariateRationalPolynomial p,check;
	std::vector<SparseMultivariateRationalPolynomial> set;
	std::vector<Symbol> vars,trcVars;
	int nvars(4);
	int nAlgVars(4);
	int nTrcVars(2);
	int nterms(8);
	unsigned long int coefBound(6ul);
	degree_t sparsity(25);
	bool includeNeg(1);
    ts.randomStronglyNormalizedTriangularSet(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    ts.display();
    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);


    vars = ts.variables();
    trcVars = ts.transcendentalVariables();
    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
    p.setRingVariables(vars);
    set = ts.polynomials();
    for (int i = 0; i < set.size(); ++i) {
        p *= set[i];
    }

    SparseMultivariateRationalPolynomial r;
	vector<SparseMultivariateRationalPolynomial> Q;
	r = ts.normalForm(p,&Q);

//	cout << "r = " << r << endl;
//	for (int i=0; i<Q.size(); ++i)
//		cout << "Q[" << i << "] = " << Q[i] << endl;

//	cerr << "compute constraint:" << endl;
	// set = ts.polynomials();
	// check.zero();
	// for (int i=0; i<set.size(); ++i)
	// 	check += set[i]*Q[i];
	// check += r;
	// check -= p;
//	cout << "check = " << check << endl;
    // if (!check.isZero()) {
    //     std::cerr << "TriangularSet normalForm test:\t\t\t\t\t\t\t FAILED" << std::endl;
    //     exit(1);
    // }

    // std::cerr << "TriangularSet normalForm test:\t\t\t\t\t\t\t PASSED" << std::endl;


    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree tsTree = ts.convertToExpressionTree();
    ExpressionTree remTree = r.convertToExpressionTree();
    std::vector<ExpressionTree> quoTrees;
    for (int i = 0; i < Q.size(); ++i) {
        quoTrees.push_back(Q[i].convertToExpressionTree());
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(tsTree.toMapleString() + ":");

    vars = ts.variables();
    std::string plex;
    if (vars.size() > 1) {
        plex = "plex(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            plex += vars[i].toString();
            plex += ",";
        }
        plex += vars[vars.size()-1].toString();
    } else {
        plex = "plex(" + vars[0].toString();
    }
    plex += "):";
    inputs.push_back(plex);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();
    char* cstr;

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "qList", 1));

    //algebList: [0] = dividend, [1] = divisorList, [2] = plex, [3] = qList

    std::string procStr = "Groebner:-NormalForm:";
    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    ALGEB result = EvalMapleProc(kv, testProc, 4, algebList[0], algebList[1], algebList[2], algebList[3]);

    // std::cerr << "calling maple proc: \n\n";
    // std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
    // std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
    // std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
    // std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
    // std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
    // std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
    // std::cerr << "\n\n";

    std::vector<ALGEB> resultAlgebs;
    inputs.clear();
    inputs.push_back(remTree.toMapleString() + ":");
    for (int i = 0; i < quoTrees.size(); ++i) {
        inputs.push_back(quoTrees[i].toMapleString() + ":");
    }
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        resultAlgebs.push_back(res);
        delete[] cstr;
    }

    char compareFunc[] = "verify:";
    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);

    char expandFunc[] = "expand:";
    ALGEB expandF = EvalMapleStatement(kv, expandFunc);

    result = EvalMapleProc(kv, expandF, 1, result);
    ALGEB comp = EvalMapleProc(kv, cmpF, 2, result, resultAlgebs[0]);
    M_BOOL compBool = MapleToM_BOOL(kv, comp);
    std::string mapleResStr = mapleTest->algebToString(kv, result);
    std::vector<ALGEB> mapleQuoList;
    int idx = 1;
    while (idx < resultAlgebs.size()) {
        std::string getQuoStr = "qList[" + std::to_string(idx) + "]:";
        cstr = new char[getQuoStr.length()+1];
        std::strcpy(cstr, getQuoStr.c_str());
        ALGEB mapleQuo = EvalMapleStatement(kv, cstr);
        mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo);
        comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
        compBool = compBool && MapleToM_BOOL(kv, comp);
        mapleQuoList.push_back(mapleQuo);
        ++idx;
        delete[] cstr;
    }

    if(!compBool) {
        std::cerr << "TriangularSet NormalForm test... FAIL" << std::endl;

        std::cerr << "Dividend: " << p << std::endl;
        std::cerr << "Divisors: " << tsTree.toMapleString() << std::endl << std::endl;

        std::cerr << "Got remainder: " << r << std::endl;
        std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;

        for (int i = 0; i < mapleQuoList.size(); ++i) {
            std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
            std::cerr << "Got maple quo : " << mapleTest->algebToString(kv, mapleQuoList[i]) << std::endl;
        }
        exit(1);
    }

    std::cerr << "TriangularSet NormalForm test... PASS" << std::endl;
}


void testMDPD() {
    const char* vars[] = {"x", "y", "z"}; // the variable symbols
    int nvar = 3; 
    long delta = 2;
    degree_t sparsity = 5;
    unsigned long coefBound= 5l;
    int includeNeg = 0;
    int flag = 0;
    int ntrials = 1;
    int lazard_flag = 0;
    time_t sed = time(NULL);
    // fprintf(stderr, "[seed: %lu]\n", sed);
    srand(sed);


    degrees_t maxDegs_f = 0;
    const int* offSet = getExpOffsetArray(nvar);
    for (int i = 0; i < nvar-1; ++i){
        maxDegs_f += (delta << offSet[i]);
    }
    maxDegs_f += ((nvar*nvar*delta) << offSet[nvar-1]);    
    Node* f = randomMaxDegsPoly(nvar, maxDegs_f, sparsity, coefBound, includeNeg);
    AltArr_t* f_AA = deepCopyPolynomial_AAFromNode(f, nvar);
    
    degrees_t maxDegs = 0;
    long int tmp = delta;
    for (int i = 0; i < nvar-1; ++i){
        tmp = rand() % delta;
        if (tmp == 0l){tmp = delta;}  
            maxDegs += (tmp << offSet[i]);
    }
    maxDegs += ((nvar*delta) << offSet[nvar-1]);
    // Node** G = (Node**) malloc (sizeof (Node*)*nvar);
    // randomTriangularSet(G, nvar, maxDegs, sparsity, coefBound, includeNeg);
    /// /// /// /// generate the triangularset:
    AltArr_t** G_AA = (AltArr_t**) malloc (sizeof(AltArr_t*)*nvar);
    randomTriangularSet_AA(G_AA, nvar, maxDegs, sparsity, coefBound, includeNeg, lazard_flag);
  
    /// /// /// /// pseudodivide:
    // void naiveMultiDivisorPseudoDivide_AA (AltArr_t* c, AltArr_t** B, AltArr_t** quoSet, AltArr_t** rem, AltArr_t** hPow, int nvar, int lazy, int nSet) 
    AltArr_t* Q[nvar];
    for (int i = 0; i < nvar; ++i){ Q[i] = NULL; }
    AltArr_t* r = NULL;
    AltArr_t* h = NULL;
    int e = 0;
    int lazy = 0;
    multiDivisorPseudoDivide_AA (f_AA, G_AA, Q, &r, &h, nvar, lazy, nvar);

    int rr = multiDivisorDivisionVerification_AA (f_AA, G_AA, Q, r, h, nvar, nvar);
    if (rr == 1){
        std::cout << "multiDivisorPseudoDivide test... PASS" << std::endl;
    } else {
        std::cout << "multiDivisorPseudoDivide test... FAIL" << std::endl;
        // std::cout << "f := " << polyToString (f, nvar, vars) << ";" << std::endl; 
        // for (int i = 0; i < nvar; ++i) {std::cout << "G[" << i << "] := " << polyToString (deepCopyPolynomial_NodeFromAA(G_AA[i]), nvar, vars) << ";" << std::endl; }
        // std::cout << "r := " << polyToString (deepCopyPolynomial_NodeFromAA(r), nvar,vars) << std::endl;
        // for (int i = 0; i < nvar; ++i) {std::cout << "Q[" << i << "] := " << polyToString (deepCopyPolynomial_NodeFromAA(Q[i]), nvar, vars) << ";" << std::endl; }
        // std::cout << "h := " << polyToString (deepCopyPolynomial_NodeFromAA(h), nvar,vars) << ";" << std::endl;
        // std::cout << "e := " << e << ";" << std::endl;
    }
}


// void normalForm_cppTest (int argc, char** argv)
// {
//   const char* vars[] = {"x", "y", "z"}; // the variable symbols
//   int i;
//   int nvar;
//   int nterms = 5;
//   int divisorSize = 3;
//   unsigned long int coefBound = 10ul;
//   degree_t sparsity = 3;
//   degree_t sparsityEx = 10;
//   bool includeNeg = 0;

//   //generate random polynomials:
//   SparseMultivariateRationalPolynomial p1(1),p2(2),p3(3), f(3);
//   p1.randomPolynomial(1,nterms,coefBound,sparsityEx,includeNeg);
//   p2.randomPolynomial(2,nterms,coefBound,sparsity,includeNeg);
//   p3.randomPolynomial(3,nterms,coefBound,sparsity,includeNeg);
//   f.randomPolynomial(3,20,coefBound,sparsity,includeNeg);
  
//   // set variable names:
  
//   vector<Symbol> names;
//   Symbol s = 'z';
//   names.insert(names.begin(),s);
//   p1.setRingVariables(names);

//   s = "y";
//   names.insert(names.begin(),s);
//   p2.setRingVariables(names);

//   s = "x";
//   names.insert(names.begin(),s);
//   p3.setRingVariables(names);
//   f.setRingVariables(names);


//   // // test random dividend:
//   std::cout << "==================================================" << std::endl;
//   std::cout << "[IN] " << "f := " << f << ";" << std::endl;

//   // // generate divisorSet: [p3(x,y,z), p2(y,z), p1(z)] 
//   std::vector<SparseMultivariateRationalPolynomial> divisorSet;
//   divisorSet.insert(divisorSet.begin(), p1);
//   divisorSet.insert(divisorSet.begin(), p2);
//   divisorSet.insert(divisorSet.begin(), p3);

//   // // test divisorSet:
//   for (int i = 0; i < divisorSize; ++i){
//     std::cout << "[IN] " << "G["<<i<<"] := " << divisorSet[i]  << ";" << std::endl;
//   }
  
//   // // generate quoSet and rem
//   std::vector<SparseMultivariateRationalPolynomial> quoSet;
//   SparseMultivariateRationalPolynomial rem;

//   // // test : multiDivisorDivision: 
//   rem = f.lexNormalForm(names, divisorSet, &quoSet);

//   // print quoSet:
//   for (int i = 0; i < divisorSize; ++i){
//       std::cout << "[OUT] " << "Q["<<i<<"] := " << quoSet[i] << ";" <<std::endl;  
//   }

//   // // print remainder:
//   std::cout << "[OUT] " << "r := " << rem << ";" << std::endl;

//   return;
// }


// void multiDivisorDivision_AA_cTest (int argc, char** argv)
// {
//   const int SHOW_POLY = 1;

//   time_t sed = time(NULL);
//   srand(sed);
//   int nvar = 3; // the number of variables
//   long delta = 2;
//   degree_t sparsity = 15;
//   unsigned long coefBound = 100;
//   int includeNeg = 0;

 
//  long delta_power = (long) pow ((double) delta, (double) nvar);
//   long delta_log = (long) ceil (log ((double) delta)/log (2));

//    degrees_t maxDegs = 0;
//    int* offSet = getExpOffsetArray(nvar);
//    maxDegs += (15l << offSet[0]);
//    maxDegs += (15l << offSet[1]);
//    maxDegs += (30l << offSet[2]);

//    degrees_t Tdegs = 0;
//    Tdegs += (11l << offSet[0]);
//    Tdegs += (11l << offSet[1]);
//    Tdegs += (20l << offSet[2]);


//   Node* f = randomMaxDegsPoly(nvar, maxDegs, sparsity, coefBound, includeNeg);
//   AltArr_t* f_AA = deepCopyPolynomial_AAFromNode(f,nvar);
//   // Node** G = (Node**) malloc (sizeof (Node*)*nvar);
//   // randomTriangularSet(G, nvar, Tdegs, sparsity, coefBound, includeNeg);
//   AltArr_t** G_AA = (AltArr_t**) malloc (sizeof (AltArr_t*)*nvar);
//   randomTriangularSet_AA(G_AA, nvar, Tdegs, sparsity, coefBound, includeNeg);
 

//   AltArr_t* r_AA = NULL;
//   AltArr_t* Qs_AA[nvar] = {NULL, NULL, NULL};
  
//   if (argc > 1){
//     multiDivisorDivision_AA (f_AA, G_AA, Qs_AA, &r_AA, nvar, nvar, std::atoi(argv[1]));  
//   } else{
//     multiDivisorDivision_AA (f_AA, G_AA, Qs_AA, &r_AA, nvar, nvar, 0);
//   }

//   if (SHOW_POLY != 0){
//     const char* vars[] = {"x", "y", "z"}; // the variable symbols

//     std::cout << "==================================================" << std::endl;
//     std::cout << " " << "f := " << polyToString (f,nvar,vars) << ";" << std::endl;
//     std::cout << " " << "T := [" << polyToString (deepCopyPolynomial_NodeFromAA(G_AA[0]),nvar,vars) << "," << std::endl;
//     std::cout << polyToString (deepCopyPolynomial_NodeFromAA(G_AA[1]),nvar,vars) << "," << std::endl;
//     std::cout << polyToString (deepCopyPolynomial_NodeFromAA(G_AA[2]),nvar,vars) << "];" << std::endl;
//     std::cout << " " << "Q := [" << polyToString (deepCopyPolynomial_NodeFromAA(Qs_AA[0]),nvar,vars) << "," << std::endl;
//     std::cout << polyToString (deepCopyPolynomial_NodeFromAA(Qs_AA[1]),nvar,vars) << "," << std::endl;
//     std::cout << polyToString (deepCopyPolynomial_NodeFromAA(Qs_AA[2]),nvar,vars) << "];" << std::endl;
//     std::cout << " " << "r := " << polyToString (deepCopyPolynomial_NodeFromAA(r_AA),nvar,vars) << ";" << std::endl;
//     std::cout << "==================================================" << std::endl;
//   }

//   // freePolynomial (f);
//   // freePolynomial (G[0]);
//   // freePolynomial (G[1]);
//   // freePolynomial (G[2]);
//   return;
// }

// void triangularSet_cTest (int argc, char** argv)
// {
//     const char* vars[] = {"x", "y", "z", "w", "u", "v", "k"}; // the variable symbols
//     time_t sed = time(NULL);
//     srand(sed);

//     int nvar = 7; // the number of variables
//     long delta = 5;
//     degree_t sparsity = 5;
//     unsigned long coefBound = 100;
//     int includeNeg = 0;

//     degrees_t maxDegs = (degrees_t) malloc (sizeof(degree_t)*nvar);
//     maxDegs[0] = 6;
//     maxDegs[1] = 6;
//     maxDegs[2] = 6;
//     maxDegs[3] = 10;
//     maxDegs[4] = 10;
//     maxDegs[5] = 12;
//     maxDegs[6] = 18;

//     Node** G = (Node**) malloc (sizeof (Node*)*nvar);
//     randomTriangularSet(G, nvar, maxDegs, sparsity, coefBound, includeNeg);

//     std::cout << "==================================================" << std::endl;
//     for (int i = 0; i < nvar; ++i)
//     {
//       std::cout << "T[" << i << "] := " <<
//       polyToString (G[i], nvar, vars) << ";" << std::endl;
//     }
//     std::cout << "==================================================" << std::endl;

//     freePolynomial(G[0]);
//     freePolynomial(G[1]);
//     freePolynomial(G[2]);
//     freePolynomial(G[3]);
//     freePolynomial(G[4]);
//     freePolynomial(G[5]);
//     freePolynomial(G[6]);
//     free(G);
//     return;
// }


// // argv[1] = type of MDD: heap=0, triagularSet=1, 2
// void multiDivisorDivision_cTest (int argc, char** argv)
// {
//   const int SHOW_POLY = 1;

//   time_t sed = time(NULL);
//   srand(sed);
//   int nvar = 3; // the number of variables
//   long delta = 2;
//   degree_t sparsity = 4;
//   unsigned long coefBound = 100;
//   int includeNeg = 0;

//   degrees_t degs = (degrees_t) malloc ( sizeof (degree_t)*nvar);
//   degrees_t Tdegs = (degrees_t) malloc ( sizeof (degree_t)*nvar);
//   long delta_power = (long) pow ((double) delta, (double) nvar);
//   long delta_log = (long) ceil (log ((double) delta)/log (2));
//   degs[0] = 2*delta_log;   degs[1] = 2*delta_log;   degs[2] = 2*delta_power;
//   Tdegs[0] = delta_log;    Tdegs[1] = delta_log;    Tdegs[2] = delta_power;

//   Node* f = randomMaxDegsPoly(nvar, degs, sparsity, coefBound, includeNeg);
//   Node** G = (Node**) malloc (sizeof (Node*)*nvar);
//   randomTriangularSet(G, nvar, Tdegs, sparsity, coefBound, includeNeg);
  
//   Node* r = NULL;
//   Node* Qs[3] = {NULL, NULL, NULL};
  
//   if (argc > 1){
//     multiDivisorDivision (f, G, Qs, &r, nvar, nvar, std::atoi(argv[1]));  
//   } else{
//     multiDivisorDivision (f, G, Qs, &r, nvar, nvar, 0);
//   }


//   if (SHOW_POLY != 0){
//     const char* vars[] = {"x", "y", "z"}; // the variable symbols

//     std::cout << "==================================================" << std::endl;
//     std::cout << "[IN] " << "f := " << polyToString (f,nvar,vars) << ";" << std::endl;
//     std::cout << "[IN] " << "T := [" << polyToString (G[0],nvar,vars) << "," << std::endl;
//     std::cout << polyToString (G[1],nvar,vars) << "," << std::endl;
//     std::cout << polyToString (G[2],nvar,vars) << "];" << std::endl;
//     std::cout << "[OUT] " << "Q := [" << polyToString (Qs[0],nvar,vars) << "," << std::endl;
//     std::cout << polyToString (Qs[1],nvar,vars) << "," << std::endl;
//     std::cout << polyToString (Qs[2],nvar,vars) << "];" << std::endl;
//     std::cout << "[OUT] " << "r := " << polyToString (r,nvar,vars) << ";" << std::endl;
//     std::cout << "==================================================" << std::endl;
//   }

//   freePolynomial (f);
//   freePolynomial (G[0]);
//   freePolynomial (G[1]);
//   freePolynomial (G[2]);
//   free (degs);
//   free (Tdegs); 
//   return;
// }


// void testPseudoDivide() {

// 	TriangularSet<RN, SparseMultivariateRationalPolynomial> ts;
// 	SparseMultivariateRationalPolynomial p, check;
// 	std::vector<SparseMultivariateRationalPolynomial> set;
// 	std::vector<Symbol> vars, trcVars;
// 	int nvars(4);
// 	int nAlgVars(3);
// 	int nTrcVars(1);
// 	int nterms(100);
// 	unsigned long int coefBound(3ul);
// 	degree_t sparsity(50);
// 	degree_t sparsity_div(100);
// 	bool includeNeg(1);
// 	ts.randomTriangularSet(nvars, nAlgVars, nTrcVars, nterms, coefBound,
// 			sparsity, includeNeg);
// //    ts.display();
// 	p.randomPolynomial(nvars + nTrcVars, nterms, coefBound, sparsity_div,
// 			includeNeg);
// 	vars = ts.variables();
// 	trcVars = ts.transcendentalVariables();
// 	vars.insert(vars.end(), trcVars.begin(), trcVars.end());
// 	p.setRingVariables(vars);
// 	set = ts.polynomials();
// 	std::reverse(set.begin(), set.end());
// //  std::cerr << "set.size = " << set.size() << std::endl;
// //    std::cerr << "[";
// //  for (int i=0; i<set.size(); ++i){
// //    std::cerr << set[i];
// //    if (i != set.size()-1)
// //      std::cerr << ",";
// //  }
// //  std::cerr << "]" << std::endl;
//   std::cerr << "ts = " << ts << std::endl;
//   std::cerr << "p = " << p << std::endl;
// 	vars = ts.variables();
// 	std::reverse(vars.begin(), vars.end());
// //    std::cerr << "[";
// //  for (int i=0; i<vars.size(); ++i){
// //    std::cerr << vars[i];
// //    if (i != vars.size()-1)
// //      std::cerr << ",";
// //  }
// //  std::cerr << "]" << std::endl;
// 	SparseMultivariateRationalPolynomial r, c;
// 	vector<SparseMultivariateRationalPolynomial> Q;
// //	r = ts.pseudoDivide(p, &Q, &c);
// 	r = ts.TSPseudoDivide(p, &Q, &c);
// //  cout << "c = " << c << endl;
// //  cout << "r = " << r << endl;
// //  for (int i=0; i<Q.size(); ++i)
// //    cout << "Q[" << i << "] = " << Q[i] << endl;

// //  cerr << "compute constraint:" << endl;
// //  set = ts.polynomials();
// //  check.zero();
// //  for (int i=0; i<set.size(); ++i)
// //    check += set[i]*Q[i];
// //  check += r;
// //  check -= (c*p);
// //  cout << "check = " << check << endl;
// //    if (!check.isZero()) {
// //        std::cerr << "TriangularSet pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;
// //        exit(1);
// //    }
// //
// //    std::cerr << "TriangularSet pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;

// //	ExpressionTree pTree = p.convertToExpressionTree();
// ////    ExpressionTree tsTree = ts.convertToExpressionTree();
// //	ExpressionTree tsTree;
// //	tsTree.fromVector<SparseMultivariateRationalPolynomial>(set);
// //	ExpressionTree cTree = c.convertToExpressionTree();
// //	ExpressionTree remTree = r.convertToExpressionTree();
// //	std::vector<ExpressionTree> quoTrees;
// //	for (int i = 0; i < Q.size(); ++i) {
// //		quoTrees.push_back(Q[i].convertToExpressionTree());
// //	}
// //
// //	std::vector<std::string> inputs;
// //	inputs.push_back(pTree.toMapleString() + ":");
// //	inputs.push_back(tsTree.toMapleString() + ":");
// //
// ////    vars = ts.variables();
// ////    std::reverse(vars.begin(),vars.end());
// //	std::string varList;
// //	if (vars.size() > 1) {
// //		varList = "[" + vars[0].toString() + ",";
// //		for (int i = 1; i < vars.size() - 1; ++i) {
// //			varList += vars[i].toString();//	ExpressionTree pTree = p.convertToExpressionTree();
//     ////    ExpressionTree tsTree = ts.convertToExpressionTree();
//     //	ExpressionTree tsTree;
//     //	tsTree.fromVector<SparseMultivariateRationalPolynomial>(set);
//     //	ExpressionTree cTree = c.convertToExpressionTree();
//     //	ExpressionTree remTree = r.convertToExpressionTree();
//     //	std::vector<ExpressionTree> quoTrees;
//     //	for (int i = 0; i < Q.size(); ++i) {
//     //		quoTrees.push_back(Q[i].convertToExpressionTree());
//     //	}
//     //
//     //	std::vector<std::string> inputs;
//     //	inputs.push_back(pTree.toMapleString() + ":");
//     //	inputs.push_back(tsTree.toMapleString() + ":");
//     //
//     ////    vars = ts.variables();
//     ////    std::reverse(vars.begin(),vars.end());
//     //	std::string varList;
//     //	if (vars.size() > 1) {
//     //		varList = "[" + vars[0].toString() + ",";
//     //		for (int i = 1; i < vars.size() - 1; ++i) {
//     //			varList += vars[i].toString();
//     //			varList += ",";
//     //		}
//     //		varList += vars[vars.size() - 1].toString();
//     //	} else {
//     //		varList = "[" + vars[0].toString();
//     //	}
//     //	varList += "]:";
//     //	inputs.push_back(varList);
//     //
//     //	MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//     //	mapleTest->restartMapleKernel();
//     //	MKernelVector kv = mapleTest->getMKernelVector();
//     //	char* cstr;
//     //
//     //	std::vector<ALGEB> algebList;
//     //	for (int i = 0; i < inputs.size(); ++i) {
//     //		cstr = new char[inputs[i].length() + 1];
//     //		std::strcpy(cstr, inputs[i].c_str());
//     //		ALGEB res = EvalMapleStatement(kv, cstr);
//     //		algebList.push_back(res);
//     //		delete[] cstr;
//     //	}
//     //	algebList.push_back(ToMapleName(kv, "qList", 1));
//     //	algebList.push_back(ToMapleName(kv, "m", 1));
//     //
//     //	//algebList: [0] = dividend, [1] = divisorList, [2] = varList, [3] = qList, [4] = m
//     //
//     ////    std::string procStr = "PseudoDivide := proc(p,ts,vs,q,m)\r\n local quo;\r\n quo:=[1,2,3];\r\n assign(q=quo);\r\n assign(m=2); 3;\r\n end proc:";
//     //	std::string procStr =
//     //			"multiDivisorPseudoDivision := proc (p, ts, names, q, m) local seqLoop, tsSize, r, tmpQ, Q, tmpM, totalM, i, j; seqLoop := [3, 2, 1]; tsSize := nops(ts); r := p; Q := [seq(0, i = 1 .. tsSize)]; totalM := 1; for i in seqLoop do tmpQ := 0; tmpM := 1; r := prem(r, ts[i], names[i], 'tmpM', 'tmpQ'); totalM := totalM*tmpM; for j to tsSize do Q[j] := simplify(tmpM*Q[j]) end do; Q[i] := Q[i]+tmpQ end do; assign(q = Q); assign(m = totalM); return r end proc:";
//     //	cstr = new char[procStr.length() + 1];
//     //	std::strcpy(cstr, procStr.c_str());
//     //	ALGEB testProc = EvalMapleStatement(kv, cstr);
//     //	delete[] cstr;
//     //
//     //	ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1],
//     //			algebList[2], algebList[3], algebList[4]);
//     //
//     ////     std::cerr << "calling maple proc: \n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
//     ////     std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
//     ////     std::cerr << "\n\n";
//     //
//     //	std::vector<ALGEB> resultAlgebs;
//     //	inputs.clear();
//     //	inputs.push_back(remTree.toMapleString() + ":");
//     //	inputs.push_back(cTree.toMapleString() + ":");
//     //	for (int i = 0; i < quoTrees.size(); ++i) {
//     //		inputs.push_back(quoTrees[i].toMapleString() + ":");
//     //	}
//     //	for (int i = 0; i < inputs.size(); ++i) {
//     //		cstr = new char[inputs[i].length() + 1];
//     //		std::strcpy(cstr, inputs[i].c_str());
//     //		ALGEB res = EvalMapleStatement(kv, cstr);
//     //		resultAlgebs.push_back(res);
//     //		delete[] cstr;
//     //	}
//     //
//     //	char compareFunc[] = "verify:";
//     //	ALGEB cmpF = EvalMapleStatement(kv, compareFunc);
//     //
//     //	char expandFunc[] = "expand:";
//     //	ALGEB expandF = EvalMapleStatement(kv, expandFunc);
//     //
//     //	result = EvalMapleProc(kv, expandF, 1, result);
//     //	ALGEB comp = EvalMapleProc(kv, cmpF, 2, result, resultAlgebs[0]);
//     //	M_BOOL compBool = MapleToM_BOOL(kv, comp);
//     //	if (!MapleToM_BOOL(kv, comp))
//     //		std::cerr << "remainder comp fail!" << std::endl;
//     //	std::string mapleResStr = mapleTest->algebToString(kv, result);
//     //	std::string getMStr = "m:";
//     //	cstr = new char[getMStr.length() + 1];
//     //	std::strcpy(cstr, getMStr.c_str());
//     //	ALGEB mapleM = EvalMapleStatement(kv, cstr);
//     //	mapleM = EvalMapleProc(kv, expandF, 1, mapleM);
//     //	comp = EvalMapleProc(kv, cmpF, 2, mapleM, resultAlgebs[1]);
//     //	compBool = compBool && MapleToM_BOOL(kv, comp);
//     //	if (!MapleToM_BOOL(kv, comp))
//     //		std::cerr << "multiplier comp fail!" << std::endl;
//     //	delete[] cstr;
//     //	std::vector<ALGEB> mapleQuoList;
//     //	int idx = 2;
//     ////    while (idx < resultAlgebs.size()) {
//     //	while (idx <= 4) {
//     //		std::string getQuoStr = "qList[" + std::to_string(nvars + 2 - idx)
//     //				+ "]:";
//     //		cstr = new char[getQuoStr.length() + 1];
//     //		std::strcpy(cstr, getQuoStr.c_str());
//     //		ALGEB mapleQuo = EvalMapleStatement(kv, cstr);
//     //		mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo);
//     //		comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
//     //		compBool = compBool && MapleToM_BOOL(kv, comp);
//     //		if (!MapleToM_BOOL(kv, comp)) {
//     //			std::cerr << "quotient " << idx - 2 << " comp fail!" << std::endl;
//     //		}
//     //		mapleQuoList.push_back(mapleQuo);
//     //		++idx;
//     //		delete[] cstr;
//     //	}
//     //
//     //	if (!compBool) {
//     //		std::cerr
//     //				<< "TriangularSet multiDivisorPseudoDivision test:\t\t\t\t\t FAILED"
//     //				<< std::endl;
//     //
//     //		std::cerr << "Dividend: " << p << std::endl;
//     //		std::cerr << "Divisors: " << ts << std::endl << std::endl;
//     //
//     //		std::cerr << "Got remainder: " << r << std::endl;
//     //		std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;
//     //
//     //		std::cerr << "Got multiplier: " << c << std::endl;
//     //		std::cerr << "Got maple m: " << mapleTest->algebToString(kv, mapleM)
//     //				<< std::endl << std::endl;
//     //
//     //		for (int i = 0; i < mapleQuoList.size(); ++i) {
//     //			std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
//     //			std::cerr << "Got maple quo : "
//     //					<< mapleTest->algebToString(kv, mapleQuoList[i])
//     //					<< std::endl;
//     //		}
//     //		exit(1);
//     //	}
//     //
//     //	std::cerr
//     //			<< "TriangularSet multiDivisorPseudoDivision test:\t\t\t\t\t PASSED"
//     //			<< std::endl;
// //			varList += ",";
// //		}
// //		varList += vars[vars.size() - 1].toString();
// //	} else {
// //		varList = "[" + vars[0].toString();
// //	}
// //	varList += "]:";
// //	inputs.push_back(varList);
// //
// //	MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
// //	mapleTest->restartMapleKernel();
// //	MKernelVector kv = mapleTest->getMKernelVector();
// //	char* cstr;
// //
// //	std::vector<ALGEB> algebList;
// //	for (int i = 0; i < inputs.size(); ++i) {
// //		cstr = new char[inputs[i].length() + 1];
// //		std::strcpy(cstr, inputs[i].c_str());
// //		ALGEB res = EvalMapleStatement(kv, cstr);
// //		algebList.push_back(res);
// //		delete[] cstr;
// //	}
// //	algebList.push_back(ToMapleName(kv, "qList", 1));
// //	algebList.push_back(ToMapleName(kv, "m", 1));
// //
// //	//algebList: [0] = dividend, [1] = divisorList, [2] = varList, [3] = qList, [4] = m
// //
// ////    std::string procStr = "PseudoDivide := proc(p,ts,vs,q,m)\r\n local quo;\r\n quo:=[1,2,3];\r\n assign(q=quo);\r\n assign(m=2); 3;\r\n end proc:";
// //	std::string procStr =
// //			"multiDivisorPseudoDivision := proc (p, ts, names, q, m) local seqLoop, tsSize, r, tmpQ, Q, tmpM, totalM, i, j; seqLoop := [3, 2, 1]; tsSize := nops(ts); r := p; Q := [seq(0, i = 1 .. tsSize)]; totalM := 1; for i in seqLoop do tmpQ := 0; tmpM := 1; r := prem(r, ts[i], names[i], 'tmpM', 'tmpQ'); totalM := totalM*tmpM; for j to tsSize do Q[j] := simplify(tmpM*Q[j]) end do; Q[i] := Q[i]+tmpQ end do; assign(q = Q); assign(m = totalM); return r end proc:";
// //	cstr = new char[procStr.length() + 1];
// //	std::strcpy(cstr, procStr.c_str());
// //	ALGEB testProc = EvalMapleStatement(kv, cstr);
// //	delete[] cstr;
// //
// //	ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1],
// //			algebList[2], algebList[3], algebList[4]);
// //
// ////     std::cerr << "calling maple proc: \n\n";
// ////     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
// ////     std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
// ////     std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
// ////     std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
// ////     std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
// ////     std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
// ////     std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
// ////     std::cerr << "\n\n";
// //
// //	std::vector<ALGEB> resultAlgebs;
// //	inputs.clear();
// //	inputs.push_back(remTree.toMapleString() + ":");
// //	inputs.push_back(cTree.toMapleString() + ":");
// //	for (int i = 0; i < quoTrees.size(); ++i) {
// //		inputs.push_back(quoTrees[i].toMapleString() + ":");
// //	}
// //	for (int i = 0; i < inputs.size(); ++i) {
// //		cstr = new char[inputs[i].length() + 1];
// //		std::strcpy(cstr, inputs[i].c_str());
// //		ALGEB res = EvalMapleStatement(kv, cstr);
// //		resultAlgebs.push_back(res);
// //		delete[] cstr;
// //	}
// //
// //	char compareFunc[] = "verify:";
// //	ALGEB cmpF = EvalMapleStatement(kv, compareFunc);
// //
// //	char expandFunc[] = "expand:";
// //	ALGEB expandF = EvalMapleStatement(kv, expandFunc);
// //
// //	result = EvalMapleProc(kv, expandF, 1, result);
// //	ALGEB comp = EvalMapleProc(kv, cmpF, 2, result, resultAlgebs[0]);
// //	M_BOOL compBool = MapleToM_BOOL(kv, comp);
// //	if (!MapleToM_BOOL(kv, comp))
// //		std::cerr << "remainder comp fail!" << std::endl;
// //	std::string mapleResStr = mapleTest->algebToString(kv, result);
// //	std::string getMStr = "m:";
// //	cstr = new char[getMStr.length() + 1];
// //	std::strcpy(cstr, getMStr.c_str());
// //	ALGEB mapleM = EvalMapleStatement(kv, cstr);
// //	mapleM = EvalMapleProc(kv, expandF, 1, mapleM);
// //	comp = EvalMapleProc(kv, cmpF, 2, mapleM, resultAlgebs[1]);
// //	compBool = compBool && MapleToM_BOOL(kv, comp);
// //	if (!MapleToM_BOOL(kv, comp))
// //		std::cerr << "multiplier comp fail!" << std::endl;
// //	delete[] cstr;
// //	std::vector<ALGEB> mapleQuoList;
// //	int idx = 2;
// ////    while (idx < resultAlgebs.size()) {
// //	while (idx <= 4) {
// //		std::string getQuoStr = "qList[" + std::to_string(nvars + 2 - idx)
// //				+ "]:";
// //		cstr = new char[getQuoStr.length() + 1];
// //		std::strcpy(cstr, getQuoStr.c_str());
// //		ALGEB mapleQuo = EvalMapleStatement(kv, cstr);
// //		mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo);
// //		comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
// //		compBool = compBool && MapleToM_BOOL(kv, comp);
// //		if (!MapleToM_BOOL(kv, comp)) {
// //			std::cerr << "quotient " << idx - 2 << " comp fail!" << std::endl;
// //		}
// //		mapleQuoList.push_back(mapleQuo);
// //		++idx;
// //		delete[] cstr;
// //	}
// //
// //	if (!compBool) {
// //		std::cerr
// //				<< "TriangularSet multiDivisorPseudoDivision test:\t\t\t\t\t FAILED"
// //				<< std::endl;
// //
// //		std::cerr << "Dividend: " << p << std::endl;
// //		std::cerr << "Divisors: " << ts << std::endl << std::endl;
// //
// //		std::cerr << "Got remainder: " << r << std::endl;
// //		std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;
// //
// //		std::cerr << "Got multiplier: " << c << std::endl;
// //		std::cerr << "Got maple m: " << mapleTest->algebToString(kv, mapleM)
// //				<< std::endl << std::endl;
// //
// //		for (int i = 0; i < mapleQuoList.size(); ++i) {
// //			std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
// //			std::cerr << "Got maple quo : "
// //					<< mapleTest->algebToString(kv, mapleQuoList[i])
// //					<< std::endl;
// //		}
// //		exit(1);
// //	}
// //
// //	std::cerr
// //			<< "TriangularSet multiDivisorPseudoDivision test:\t\t\t\t\t PASSED"
// //			<< std::endl;
// }
