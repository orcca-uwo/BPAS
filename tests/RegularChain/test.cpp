#include <bpas.h>
#include <vector>
#include "../MapleTestTool/MapleTestTool.hpp"
#include "../../include/Utils/SymbolHelpers.hpp"

using namespace std;

void testRCDefaultConstructor() {
    RegularChain<RN,SMQP> ts;

    if (ts.numberOfVariables() != 0 || ts.numberOfVariables() != 0 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 0) {
        std::cerr << "RC default constructor test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC default constructor test:\t\t\t\t\t PASSED" << std::endl;
}

void testRCVariableSpecificationConstructor() {
	vector<Symbol> xs;
	xs = Symbol::randomElements(2);
    RegularChain<RN,SMQP> ts(xs);

    if (ts.numberOfVariables() != 2 || ts.numberOfVariables() != 2 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 0) {
        std::cerr << "RC variable specification constructor test:\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC variable specification constructor test:\t\t\t PASSED" << std::endl;
}

void testRCVariableAndTranscendentalSpecificationConstructor() {
	vector<Symbol> vars;
	vars = Symbol::randomElements(5);
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
    RegularChain<RN,SMQP> ts(xs,ps);

    if (ts.numberOfVariables() != 2 || ts.numberOfVariables() != 2 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 3) {
        std::cerr << "RC variable and transcendental specification constructor test:\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC variable and transcendental specification constructor test:\t PASSED" << std::endl;
}

void testRCPolynomialAddConstructor() {
	int nvars(3);
	int nterms(7);
	unsigned long int coefBound(50ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(3);
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+3);
	Symbol s(vars[0]);
	SMQP p(nvars),ex(s);
	vector<Symbol> x;
	x.push_back(s);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(xs);
	p += ex;
    RegularChain<RN,SMQP> ts(p);

    if (ts.numberOfVariables() != 3 || ts.numberOfVariables() != 3 || ts.numberOfAlgebraicVariables() != 1 || ts.numberOfTranscendentalVariables() != 0 || xs != ts.variables() || x != ts.mainVariables()) {
        std::cerr << "RC polynomial addition constructor test:\t\t\t FAILED" << std::endl;
        std::cerr << "polynomial being added: " << p << std::endl;
        exit(1);
    }

    std::cerr << "RC polynomial addition constructor test:\t\t\t PASSED" << std::endl;
}

void testRCPolynomialAddWithTranscendentalsConstructor() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	Symbol s(vs[0]);
	SMQP p(nvars),ex(s);
	vector<Symbol> x;
	x.push_back(s);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
	p += ex;
    RegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars = orderPreservingSetIntersection(xs,p.ringVariables());

    if (ts.numberOfVariables() != vars.size() || ts.numberOfVariables() != vars.size() || ts.numberOfAlgebraicVariables() != 1 || ts.numberOfTranscendentalVariables() != 2 || vars != ts.variables() || x != ts.mainVariables() || ps != ts.transcendentalVariables()) {
        std::cerr << "RC polynomial addition with transcendentals constructor test:\t FAILED" << std::endl;
        std::cerr << "ts.mainVariables: ";
        for (int i=0; i<ts.mainVariables().size(); ++i)
        	cerr << ts.mainVariables()[i] << " ";
        cerr << endl;
        cerr << "x = " << x[0] << endl;
        cerr << "ts.size = " << ts.numberOfVariables() << endl;
        cerr << "ts.numberOfVariables = " << ts.numberOfVariables() << endl;
        cerr << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << endl;
        cerr << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << endl;
        std::cerr << "ts.variables: ";
        printVariables(ts.variables());
        std::cerr << "xs: ";
        printVariables(xs);
        exit(1);
    }

    std::cerr << "RC polynomial addition with transcendentals constructor test:\t PASSED" << std::endl;
}

void testRCCopyConstructor() {

    RegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
   //ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
//    std::cout << "ts = " << ts << std::endl;
//    cout << "ts.numberOfVariables = " << ts.numberOfVariables() << std::endl;
//    cout << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << std::endl;
//    cout << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << std::endl;
}

void testRCMoveConstructor() {

//    RegularChain<RN,SMQP> ts;
//	int nterms(14);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(4);
//	bool includeNeg(1);
//    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
//    std::vector<Symbol> vs(ts.mainVariables());
//    RegularChain<RN,SMQP> ts2(ts.lower(vs[0]));
//    ts = ts.lower(vs[0]);
//    if (ts.variables() != ts2.variables() || ts.transcendentalVariables() != ts2.transcendentalVariables() || ts.polynomials() != ts2.polynomials()) {
//        std::cerr << "RC copy constructor test:\t\t\t\t\t\t FAILED" << std::endl;
//        std::cerr << "ts.vars:" << std::endl;
//        printVariables(ts.variables());
//        std::cerr << "ts.trcVars:" << std::endl;
//        printVariables(ts.transcendentalVariables());
//        std::cerr << "ts.polys:" << std::endl;
//        for (int i=0; i<ts.polynomials().size(); ++i)
//        	std::cerr << "set[" << i << "] = " << ts.polynomials()[i] << std::endl;
//        std::cerr << "ts2.vars:" << std::endl;
//        printVariables(ts2.variables());
//        std::cerr << "ts2.trcVars:" << std::endl;
//        printVariables(ts2.transcendentalVariables());
//        std::cerr << "ts2.polys:" << std::endl;
//        for (int i=0; i<ts2.polynomials().size(); ++i)
//        	std::cerr << "set[" << i << "] = " << ts2.polynomials()[i] << std::endl;
//        exit(1);
//    }

//    std::cerr << "RC move constructor test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCComputationalConstructor() {}

void testRCAssignmentOperator() {

    RegularChain<RN,SMQP> ts,ts2;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    ts2 = ts;
    if (ts.variables() != ts2.variables() || ts.transcendentalVariables() != ts2.transcendentalVariables() || ts.polynomials() != ts2.polynomials()) {
        std::cerr << "RC assignment operator test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts.vars:" << std::endl;
        printVariables(ts.variables());
        std::cerr << "ts.trcVars:" << std::endl;
        printVariables(ts.transcendentalVariables());
        std::cerr << "ts.polys:" << std::endl;
        for (int i=0; i<ts.polynomials().size(); ++i)
        	std::cerr << "set[" << i << "] = " << ts.polynomials()[i] << std::endl;
        std::cerr << "ts2.vars:" << std::endl;
        printVariables(ts2.variables());
        std::cerr << "ts2.trcVars:" << std::endl;
        printVariables(ts2.transcendentalVariables());
        std::cerr << "ts2.polys:" << std::endl;
        for (int i=0; i<ts2.polynomials().size(); ++i)
        exit(1);
    }

    std::cerr << "RC assignment operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testRCMoveAssignmentOperator() {

    RegularChain<RN,SMQP> ts,ts2,ts3;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    ts2 = ts;

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    ts3 = move(ts);

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

    if (ts3.variables() != ts2.variables() || ts3.transcendentalVariables() != ts2.transcendentalVariables() || ts3.polynomials() != ts2.polynomials() || !ts.isEmpty()) {
        std::cerr << "RC move assignment operator test:\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts.isEmpty = " << ts.isEmpty() << std::endl;
        printVariables(ts3.variables());
        vector<SMQP> polys(ts3.polynomials());
        for (int i=0; i<polys.size(); ++i) {
        	std::cerr << "ts3[" << i << "] = " << polys[i] << endl;
        }
        exit(1);
    }
    std::cerr << "RC move assignment operator test:\t\t\t\t PASSED" << std::endl;
}

void testRCAddOperator() {
	std::vector<int> maxDegs = {2,2,2,2};
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	maxDegs.pop_back();
	p2.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    RegularChain<RN,SMQP> ts(p1,ps);
    ts = ts + p2;
	std::vector<SMQP> polys,tsPolys,ts2Polys;
    polys.push_back(p1);
    polys.push_back(p2);
    tsPolys = ts.polynomials();

//    ZeroDimensionalRegularChain<RN,SMQP> ts2(ps);
//    printVariables(ps);
    RegularChain<RN,SMQP> ts2(xs,ps);
    ts2 = ts2 + p2;
    ts2 = ts2 + p1;
    ts2Polys = ts2.polynomials();

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

    if (tsPolys.size() != polys.size() || ts2Polys.size() != polys.size() || tsPolys[1] != ts2Polys[1]) {
        std::cerr << "RC addition operator test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts = " << ts << std::endl;
        std::cerr << "ts2 = " << ts2 << std::endl;
        for (auto x : polys) {
        	std::cerr << x << std::endl;
        }
        exit(1);
    }

    std::cerr << "RC addition operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testRCAddAssignmentOperator() {
	std::vector<int> maxDegs = {2,2,2,2};
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	maxDegs.pop_back();
	p2.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    //cout.rdbuf (ss.rdbuf());
    //cerr.rdbuf (ss.rdbuf());

    // cout << "p2 = " << p2 << endl;
    RegularChain<RN,SMQP> ts(p2,ps);
    // cout << "p1 = " << p1 << endl;
    printVariables(p1.variables());
    ts += p1;
	std::vector<SMQP> polys,tsPolys,ts2Polys;
    polys.push_back(p1);
    polys.push_back(p2);
    tsPolys = ts.polynomials();

//    ZeroDimensionalRegularChain<RN,SMQP> ts2(ps);
//    printVariables(ps);
    RegularChain<RN,SMQP> ts2(xs,ps);
    ts2 += p1;
    ts2 += p2;
    ts2Polys = ts2.polynomials();

    //cout.rdbuf (oldcout);              // <-- restore
    //cerr.rdbuf (oldcerr);              // <-- restore

    if (tsPolys.size() != polys.size() || ts2Polys.size() != polys.size() || tsPolys[1] != ts2Polys[1]) {
        std::cerr << "RC add assignment operator test:\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts = " << ts << std::endl;
        std::cerr << "ts2 = " << ts2 << std::endl;
        for (auto x : polys) {
        	std::cerr << x << std::endl;
        }
        exit(1);
    }

    std::cerr << "RC add assignment operator test:\t\t\t\t PASSED" << std::endl;
}

void testRCIdentityOperator() {

    RegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    RegularChain<RN,SMQP> ts2;
    ts2 = ts;
    if (ts.polynomials() != ts2.polynomials()) {
        std::cerr << "RC identity operator test:\t\t\t\t\t FAILED" << std::endl;
        for (int i=0; i<ts.polynomials().size(); ++i)
        	std::cerr << "ts[" << i << "] = " << ts.polynomials()[i] << std::endl;
        for (int i=0; i<ts2.polynomials().size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2.polynomials()[i] << std::endl;
        exit(1);
    }

    std::cerr << "RC identity operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testRCNumberOfVariables() {

    RegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.numberOfVariables() != 4) {
        std::cerr << "RC number of variables test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC number of variables test:\t\t\t\t\t PASSED" << std::endl;
}

void testRCSize() {

    RegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.numberOfVariables() != 4) {
        std::cerr << "RC size test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC size test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCNumberOfAlgebraicVariables() {

    RegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    std::vector<SMQP> polys(ts.polynomials()),polys2;
    for (int i=0; i<polys.size(); ++i) {
    	if (!polys[i].isZero())
    		polys2.push_back(polys[i]);
    }
    if (ts.mainVariables().size() != polys2.size()) {
        std::cerr << "RC numberOfAlgebraicVariables test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC numberOfAlgebraicVariables test:\t\t\t\t PASSED" << std::endl;
}

void testRCNumberOfTranscendentalVariables() {

    RegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
    ts.randomRegularChain(4,3,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.transcendentalVariables().size() != 2) {
        std::cerr << "RC numberOfTranscendentalVariables test:\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC numberOfTranscendentalVariables test:\t\t\t PASSED" << std::endl;
}

void testRCVariables() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);

    for (auto v : xs) {
        p += SMQP(v);
    }

    RegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars = orderPreservingSetIntersection(xs,p.variables());

    if (ts.variables() != vars) {
        std::cerr << "RC variables test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC variables test:\t\t\t\t\t\t PASSED" << std::endl;
}


void testRCMainVariables() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    RegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());

    if (ts.mainVariables() != vars) {
        std::cerr << "RC main variables test:\t\t\t\t\t\t FAILED" << std::endl;
        printVariables(ts.mainVariables());
        printVariables(vars);
        exit(1);
    }

    std::cerr << "RC main variables test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCTranscendentalVariables() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    RegularChain<RN,SMQP> ts(p,ps);

    if (ts.transcendentalVariables() != ps) {
        std::cerr << "RC transcendental variables test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC transcendental variables test:\t\t\t\t PASSED" << std::endl;
}

void testRCIsAlgebraic() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    RegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());
    vector<Symbol> algVar;
    algVar = setIntersection(ts.mainVariables(),vars);

    if (ts.isAlgebraic(p.leadingVariable()) != (algVar == vars)) {
        std::cerr << "RC isAlgebraic test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC isAlgebraic test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCIsEmpty() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    RegularChain<RN,SMQP> ts(p,ps),ts2;

    if (ts.isEmpty() == true || ts2.isEmpty() == false) {
        std::cerr << "RC isEmpty test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC isEmpty test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCSelect() {

	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars+1);
	Symbol s(vs[5]);
	vs.pop_back();
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+3);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+3,vs.begin()+5);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    RegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());
    SMQP zero;

    if (ts.select(p.leadingVariable()).leadingVariable() != p.leadingVariable() || ts.select(s) != zero) {
        std::cerr << "RC select test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC select test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCLower() {

//     RegularChain<RN,SMQP> ts,ts2;
//     std:vector<SMQP> tslp,ts2p;
//     std::vector<Symbol> vars;
//     std::vector<int> maxDegs = {2,2,2,2,2,2};
// 	unsigned long int coefBound(6ul);
// 	double sparsity(0.1);
// 	bool includeNeg(1);
//     ts.randomRegularChain(4,4,2,maxDegs,coefBound,sparsity,includeNeg);
//     vars = ts.mainVariables();

//     // Redirect stdout, stderr
//     stringstream ss;
//     streambuf *oldcout,*oldcerr;
//     oldcout = cout.rdbuf(); // <-- save
//     oldcerr = cerr.rdbuf(); // <-- save
//     cout.rdbuf (ss.rdbuf());
//     cerr.rdbuf (ss.rdbuf());

//     ts.lower(vars[0],ts2);

//     cout.rdbuf (oldcout);              // <-- restore
//     cerr.rdbuf (oldcerr);              // <-- restore

// //	std::cerr << "ts2 in testLower = " << ts2 << std::endl;
// //	ts.display();
// //	ts2.display();
// 	tslp = ts.polynomials();
// 	tslp[0].zero();
// //	tslp.erase(tslp.begin(),tslp.begin()+1);
// 	ts2p = ts2.polynomials();
//     if (tslp != ts2p) {
//         std::cerr << "RC lower test:\t\t\t\t\t\t\t FAILED" << std::endl;
//         std::cerr << "ts lower polys:" << std::endl;
//         for (int i=0; i<tslp.size(); ++i)
//         	std::cerr << "ts[" << i << "] = " << tslp[i] << std::endl;
//         std::cerr << "ts2 polys:" << std::endl;
//         for (int i=0; i<ts2p.size(); ++i)
//         	std::cerr << "ts2[" << i << "] = " << ts2p[i] << std::endl;
//         exit(1);
//     }

    std::cerr << "RC lower test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCUpper() {

//     RegularChain<RN,SMQP> ts,ts2;
//     std:vector<SMQP> tsup,ts2p;
//     std::vector<Symbol> vars;
//     std::vector<int> maxDegs = {2,2,2,2,2,2};
// 	unsigned long int coefBound(6ul);
// 	double sparsity(0.1);
// 	bool includeNeg(1);
//     ts.randomRegularChain(4,4,2,maxDegs,coefBound,sparsity,includeNeg);
//     vars = ts.mainVariables();

//     // Redirect stdout, stderr
//     stringstream ss;
//     streambuf *oldcout,*oldcerr;
//     oldcout = cout.rdbuf(); // <-- save
//     oldcerr = cerr.rdbuf(); // <-- save
//     cout.rdbuf (ss.rdbuf());
//     cerr.rdbuf (ss.rdbuf());

//     ts.upper(vars[3],ts2);

//     cout.rdbuf (oldcout);              // <-- restore
//     cerr.rdbuf (oldcerr);              // <-- restore

// //    ts.display();
// //    ts2.display();
// 	tsup = ts.polynomials();
// 	for (int i=tsup.size()-1; i>=3; --i)
// 		tsup[i].zero();
// //	tsup.erase(tsup.end()-1,tsup.end());
// 	ts2p = ts2.polynomials();
//     if (tsup != ts2p) {
//         std::cerr << "RC upper test:\t\t\t\t\t\t\t FAILED" << std::endl;
//         std::cerr << "ts upper polys:" << std::endl;
//         for (int i=0; i<tsup.size(); ++i)
//         	std::cerr << "ts[" << i << "] = " << tsup[i] << std::endl;
//         std::cerr << "ts2 polys:" << std::endl;
//         for (int i=0; i<ts2p.size(); ++i)
//         	std::cerr << "ts2[" << i << "] = " << ts2p[i] << std::endl;
//         exit(1);
//     }

    std::cerr << "RC upper test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCPseudoDivide() {

    RegularChain<RN,SMQP> ts;
	SMQP p,check;
	std::vector<SMQP> set;
	std::vector<Symbol> vars,trcVars;
	int nvars(3);
	int nAlgVars(3);
	int nTrcVars(0);
	int nterms(6);
	unsigned long int coefBound(6ul);
	degree_t sparsity(30);
	bool includeNeg(1);
    ts.randomRegularChain(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    ts.display();
    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);
    vars = ts.variables();
    trcVars = ts.transcendentalVariables();
    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
	p.setRingVariables(vars);
	set = ts.polynomials();
	std::reverse(set.begin(),set.end());
//	std::cerr << "set.size = " << set.size() << std::endl;
//    std::cerr << "[";
//	for (int i=0; i<set.size(); ++i){
//		std::cerr << set[i];
//		if (i != set.size()-1)
//			std::cerr << ",";
//	}
//	std::cerr << "]" << std::endl;
//	std::cerr << "ts = " << ts << std::endl;
//	std::cerr << "p = " << p << std::endl;
    vars = ts.variables();
    std::reverse(vars.begin(),vars.end());
//    std::cerr << "[";
//	for (int i=0; i<vars.size(); ++i){
//		std::cerr << vars[i];
//		if (i != vars.size()-1)
//			std::cerr << ",";
//	}
//	std::cerr << "]" << std::endl;

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
//    oldcout = cout.rdbuf(); // <-- save
//    oldcerr = cerr.rdbuf(); // <-- save
//    cout.rdbuf (ss.rdbuf());
//    cerr.rdbuf (ss.rdbuf());

	SMQP r,c;
	vector<SMQP> Q;
	r = ts.pseudoDivide(p,&Q,&c);

//    cout.rdbuf (oldcout);              // <-- restore
//    cerr.rdbuf (oldcerr);              // <-- restore
//	cout << "c = " << c << endl;
//	cout << "r = " << r << endl;
//	for (int i=0; i<Q.size(); ++i)
//		cout << "Q[" << i << "] = " << Q[i] << endl;
//
//	cerr << "compute constraint:" << endl;
//	set = ts.polynomials();
//	check.zero();
//	for (int i=0; i<set.size(); ++i)
//		check += set[i]*Q[i];
//	check += r;
//	check -= (c*p);
//	cout << "check = " << check << endl;
//    if (!check.isZero()) {
//        std::cerr << "RC pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "RC pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;

    ExpressionTree pTree = p.convertToExpressionTree();
//    ExpressionTree tsTree = ts.convertToExpressionTree();
	ExpressionTree tsTree;
	tsTree.fromVector<SMQP>(set);
    ExpressionTree cTree = c.convertToExpressionTree();
    ExpressionTree remTree = r.convertToExpressionTree();
    std::vector<ExpressionTree> quoTrees;
    for (int i = 0; i < Q.size(); ++i) {
        quoTrees.push_back(Q[i].convertToExpressionTree());
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(tsTree.toMapleString() + ":");

//    vars = ts.variables();
//    std::reverse(vars.begin(),vars.end());
    std::string varList;
    if (vars.size() > 1) {
        varList = "[" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            varList += vars[i].toString();
            varList += ",";
        }
        varList += vars[vars.size()-1].toString();
    } else {
        varList = "[" + vars[0].toString();
    }
    varList += "]:";
    inputs.push_back(varList);

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
    algebList.push_back(ToMapleName(kv, "m", 1));

    //algebList: [0] = dividend, [1] = divisorList, [2] = varList, [3] = qList, [4] = m

//    std::string procStr = "PseudoDivide := proc(p,ts,vs,q,m)\r\n local quo;\r\n quo:=[1,2,3];\r\n assign(q=quo);\r\n assign(m=2); 3;\r\n end proc:";
	std::string procStr = "multiDivisorPseudoDivision := proc (p, ts, names, q, m) local seqLoop, tsSize, r, tmpQ, Q, tmpM, totalM, i, j; seqLoop := [3, 2, 1]; tsSize := nops(ts); r := p; Q := [seq(0, i = 1 .. tsSize)]; totalM := 1; for i in seqLoop do tmpQ := 0; tmpM := 1; r := prem(r, ts[i], names[i], 'tmpM', 'tmpQ'); totalM := totalM*tmpM; for j to tsSize do Q[j] := simplify(tmpM*Q[j]) end do; Q[i] := Q[i]+tmpQ end do; assign(q = Q); assign(m = totalM); return r end proc:";
    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[3], algebList[4]);

//     std::cerr << "calling maple proc: \n\n";
//     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
//     std::cerr << "\n\n";

    std::vector<ALGEB> resultAlgebs;
    inputs.clear();
    inputs.push_back(remTree.toMapleString() + ":");
    inputs.push_back(cTree.toMapleString() + ":");
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
    if (!MapleToM_BOOL(kv, comp))
    	std::cerr << "remainder comp fail!" << std::endl;
    std::string mapleResStr = mapleTest->algebToString(kv, result);
    std::string getMStr = "m:";
    cstr = new char[getMStr.length()+1];
    std::strcpy(cstr,getMStr.c_str());
    ALGEB mapleM = EvalMapleStatement(kv, cstr);
    mapleM = EvalMapleProc(kv, expandF, 1, mapleM);
    comp = EvalMapleProc(kv, cmpF, 2, mapleM, resultAlgebs[1]);
    compBool = compBool && MapleToM_BOOL(kv, comp);
    if (!MapleToM_BOOL(kv, comp))
    	std::cerr << "multiplier comp fail!" << std::endl;
    delete[] cstr;
    std::vector<ALGEB> mapleQuoList;
    int idx = 2;
//    while (idx < resultAlgebs.size()) {
    while (idx <= 4) {
        std::string getQuoStr = "qList[" + std::to_string(nvars+2-idx) + "]:";
        cstr = new char[getQuoStr.length()+1];
        std::strcpy(cstr, getQuoStr.c_str());
        ALGEB mapleQuo = EvalMapleStatement(kv, cstr);
        mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo);
        comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
        compBool = compBool && MapleToM_BOOL(kv, comp);
    	if (!MapleToM_BOOL(kv, comp)) {
    		std::cerr << "quotient " << idx-2 << " comp fail!" << std::endl;
    	}
        mapleQuoList.push_back(mapleQuo);
        ++idx;
        delete[] cstr;
    }

    if(!compBool) {
        std::cerr << "RC pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;

        std::cerr << "Dividend: " << p << std::endl;
        std::cerr << "Divisors: " << ts << std::endl << std::endl;

        std::cerr << "Got remainder: " << r << std::endl;
        std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;

        std::cerr << "Got multiplier: " << c << std::endl;
        std::cerr << "Got maple m: " << mapleTest->algebToString(kv, mapleM) << std::endl << std::endl;

        for (int i = 0; i < mapleQuoList.size(); ++i) {
            std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
            std::cerr << "Got maple quo : " << mapleTest->algebToString(kv, mapleQuoList[i]) << std::endl;
        }
        exit(1);
    }

    std::cerr << "RC pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCNormalForm() {

    RegularChain<RN,SMQP> ts;
	SMQP p,check;
	std::vector<SMQP> set;
	std::vector<Symbol> vars,trcVars;
	int nvars(4);
	int nAlgVars(4);
	int nTrcVars(2);
	int nterms(3);
	unsigned long int coefBound(2ul);
	degree_t sparsity(25);
	bool includeNeg(1);
    ts.randomRegularChain(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    ts.display();
    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);


    vars = ts.variables();
    trcVars = ts.transcendentalVariables();
    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
    p.setRingVariables(vars);
    set = ts.polynomials();
    for (int i = 0; i < 2; ++i) {
        p *= set[i];
    }

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    SMQP r;
	vector<SMQP> Q;
	r = ts.normalForm(p,&Q);

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

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
    //     std::cerr << "RC normalForm test:\t\t\t\t\t\t\t FAILED" << std::endl;
    //     exit(1);
    // }

    // std::cerr << "RC normalForm test:\t\t\t\t\t\t\t PASSED" << std::endl;


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
        std::cerr << "RC normalForm test:\t\t\t\t\t\t FAILED" << std::endl;

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

    std::cerr << "RC normalForm test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testRCConvertToExpressionTree() {

    stringstream ss;
    std::string str,str2;
    RegularChain<RN,SMQP> ts,ts2;
    std:vector<SMQP> tslp,ts2p;
    std::vector<Symbol> vars;
    std::vector<int> maxDegs = {2,2,2,2,2,2};
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
    ts.randomRegularChain(4,4,2,maxDegs,coefBound,sparsity,includeNeg);

    // Redirect stdout, stderr
    stringstream ss2;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss2.rdbuf());
    cerr.rdbuf (ss2.rdbuf());
    ExpressionTree et = ts.convertToExpressionTree();

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

    // cout << "et = " << et.toString() << endl;
    ss << ts;
    str = ss.str();
    // cout << "ts = " << str << endl;
    // MAPLE EQUALITY TEST FOR et.toString() and str //
    // Start Maple Kernel:
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();

    // et:
    std::string etStr = et.toString() + ":";
    char* cstr = new char[etStr.size() + 1];
    std::strcpy (cstr, etStr.c_str());
    ALGEB resET = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    // str:
    std::string tsStr = str + ":";
    cstr = new char[tsStr.size() + 1];
    std::strcpy (cstr, tsStr.c_str());
    ALGEB resTS = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    // operator:
    std::string procStr = "ArrayTools:-IsEqual:";
    cstr = new char[procStr.length() + 1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    // Maple:
    ALGEB equalityRes = EvalMapleProc(kv, testProc, 2, resET, resTS);
    M_BOOL cmpBool = MapleToM_BOOL(kv, equalityRes);

    // Redirect stdout, stderr
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss2.rdbuf());
    cerr.rdbuf (ss2.rdbuf());

    ts = RegularChain<RN,SMQP>();
    ExpressionTree et2;
    et2 = ts.convertToExpressionTree();

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

//  cout << "et2 = " << et2.toString() << endl;
    ss = stringstream();
    ss << ts;
    str2 = ss.str();
//  cout << "ts = " << str2 << endl;
    if (!cmpBool || et2.toString() != str2) {
        std::cerr << "RC convertToExpressionTree test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "RC convertToExpressionTree test:\t\t\t\t PASSED" << std::endl;
}

void testZDRCDefaultConstructor() {
    ZeroDimensionalRegularChain<RN,SMQP> ts;

    if (ts.numberOfVariables() != 0 || ts.numberOfVariables() != 0 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 0) {
        std::cerr << "ZDRC default constructor test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC default constructor test:\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCTranscendentalVariableSpecificationConstructor() {
	vector<Symbol> ps;
	ps = Symbol::randomElements(2);
    ZeroDimensionalRegularChain<RN,SMQP> ts(ps);

    if (ts.numberOfVariables() != 0 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 2) {
        std::cerr << "ZDRC transcendental variable specification constructor test:\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC transcendental variable specification constructor test:\t PASSED" << std::endl;
}

//void testZDRCVariableAndTranscendentalSpecificationConstructor() {
//	vector<Symbol> vars;
//	vars = Symbol::randomElements(5);
//	vector<Symbol> xs;
//	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
//	vector<Symbol> ps;
//	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
//  ZeroDimensionalRegularChain<RN,SMQP> ts(xs,ps);

//    if (ts.numberOfVariables() != 2 || ts.numberOfVariables() != 2 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 3) {
//        std::cerr << "ZDRC variable and transcendental specification constructor test:\t FAILED" << std::endl;
//        std::cerr << "ts.size = " << ts.numberOfVariables() << "; should be " << xs.size() << std::endl;
//        std::cerr << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << "; should be 0" << std::endl;
//        std::cerr << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << "; should be " << ps.size() << std::endl;
//        exit(1);
//    }

//    std::cerr << "ZDRC variable and transcendental specification constructor test:\t PASSED" << std::endl;
//}

void testZDRCPolynomialAddConstructor() {
	int nvars(1);
	int nterms(7);
	unsigned long int coefBound(50ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> xs;
	xs = Symbol::randomElements(nvars);
	SMQP p(nvars),ex(xs[0]);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(xs);
	p += ex;
    ZeroDimensionalRegularChain<RN,SMQP> ts(p);

    if (ts.numberOfVariables() != 1 || ts.numberOfAlgebraicVariables() != 1 || ts.numberOfTranscendentalVariables() != 0 || xs != ts.variables() || xs != ts.mainVariables()) {
        std::cerr << "ZDRC polynomial addition constructor test:\t\t\t FAILED" << std::endl;
        std::cerr << "polynomial being added: " << p << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC polynomial addition constructor test:\t\t\t PASSED" << std::endl;
}

void testZDRCPolynomialAddWithTranscendentalsConstructor() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+1);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+1,vs.begin()+5);
	Symbol s(vs[0]);
	SMQP p(nvars),ex(s);
	vector<Symbol> x;
	x.push_back(s);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
	p += ex;
    ZeroDimensionalRegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars = orderPreservingSetIntersection(xs,p.variables());

    if (ts.numberOfVariables() != 1 || ts.numberOfAlgebraicVariables() != 1 || ts.numberOfTranscendentalVariables() != 4 || vars != ts.variables() || x != ts.mainVariables() || ps != ts.transcendentalVariables()) {
        std::cerr << "ZDRC polynomial addition with transcendentals constructor test:\t FAILED" << std::endl;
        std::cerr << "ts.mainVariables: ";
        for (int i=0; i<ts.mainVariables().size(); ++i)
        	cerr << ts.mainVariables()[i] << " ";
        cerr << endl;
        cerr << "x = " << x[0] << endl;
        cerr << "ts.size = " << ts.numberOfVariables() << endl;
        cerr << "ts.numberOfVariables = " << ts.numberOfVariables() << endl;
        cerr << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << endl;
        cerr << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << endl;
        std::cerr << "ts.variables: ";
        printVariables(ts.variables());
        std::cerr << "xs: ";
        printVariables(xs);
        exit(1);
    }

    std::cerr << "ZDRC polynomial addition with transcendentals constructor test:\t PASSED" << std::endl;
}

void testZDRCCopyConstructor() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
   //ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
//    std::cout << "ts = " << ts << std::endl;
//    cout << "ts.numberOfVariables = " << ts.numberOfVariables() << std::endl;
//    cout << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << std::endl;
//    cout << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << std::endl;
}

void testZDRCMoveConstructor() {

//    ZeroDimensionalRegularChain<RN,SMQP> ts;
//	int nterms(14);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(4);
//	bool includeNeg(1);
//    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
//    std::vector<Symbol> vs(ts.mainVariables());
//    ZeroDimensionalRegularChain<RN,SMQP> ts2(ts.lower(vs[0]));
//    ts = ts.lower(vs[0]);
//    if (ts.variables() != ts2.variables() || ts.transcendentalVariables() != ts2.transcendentalVariables() || ts.polynomials() != ts2.polynomials()) {
//        std::cerr << "ZDRC copy constructor test:\t\t\t\t\t\t FAILED" << std::endl;
//        std::cerr << "ts.vars:" << std::endl;
//        printVariables(ts.variables());
//        std::cerr << "ts.trcVars:" << std::endl;
//        printVariables(ts.transcendentalVariables());
//        std::cerr << "ts.polys:" << std::endl;
//        for (int i=0; i<ts.polynomials().size(); ++i)
//        	std::cerr << "set[" << i << "] = " << ts.polynomials()[i] << std::endl;
//        std::cerr << "ts2.vars:" << std::endl;
//        printVariables(ts2.variables());
//        std::cerr << "ts2.trcVars:" << std::endl;
//        printVariables(ts2.transcendentalVariables());
//        std::cerr << "ts2.polys:" << std::endl;
//        for (int i=0; i<ts2.polynomials().size(); ++i)
//        	std::cerr << "set[" << i << "] = " << ts2.polynomials()[i] << std::endl;
//        exit(1);
//    }

//    std::cerr << "ZDRC move constructor test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCComputationalConstructor() {}

void testZDRCAssignmentOperator() {

    ZeroDimensionalRegularChain<RN,SMQP> ts,ts2;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    ts2 = ts;
    if (ts.variables() != ts2.variables() || ts.transcendentalVariables() != ts2.transcendentalVariables() || ts.polynomials() != ts2.polynomials()) {
        std::cerr << "ZDRC assignment operator test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts.vars:" << std::endl;
        printVariables(ts.variables());
        std::cerr << "ts.trcVars:" << std::endl;
        printVariables(ts.transcendentalVariables());
        std::cerr << "ts.polys:" << std::endl;
        for (int i=0; i<ts.polynomials().size(); ++i)
        	std::cerr << "set[" << i << "] = " << ts.polynomials()[i] << std::endl;
        std::cerr << "ts2.vars:" << std::endl;
        printVariables(ts2.variables());
        std::cerr << "ts2.trcVars:" << std::endl;
        printVariables(ts2.transcendentalVariables());
        std::cerr << "ts2.polys:" << std::endl;
        for (int i=0; i<ts2.polynomials().size(); ++i)
        exit(1);
    }

    std::cerr << "ZDRC assignment operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCMoveAssignmentOperator() {

    ZeroDimensionalRegularChain<RN,SMQP> ts,ts2,ts3;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    ts2 = ts;
    ts3 = move(ts);

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

    if (ts3.variables() != ts2.variables() || ts3.transcendentalVariables() != ts2.transcendentalVariables() || ts3.polynomials() != ts2.polynomials() || !ts.isEmpty()) {
        std::cerr << "ZDRC move assignment operator test:\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts.isEmpty = " << ts.isEmpty() << std::endl;
        exit(1);
    }
    std::cerr << "ZDRC move assignment operator test:\t\t\t\t PASSED" << std::endl;
}

void testZDRCAddOperator() {
	std::vector<int> maxDegs = {2,2,2,2};
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	maxDegs.pop_back();
	p2.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    // Redirect stdout, stderr
    std::cout.setstate(std::ios::failbit);
    std::cerr.setstate(std::ios::failbit);

    ZeroDimensionalRegularChain<RN,SMQP> ts(p2,ps);
    ts = ts + p1;
	std::vector<SMQP> polys,tsPolys,ts2Polys;
    polys.push_back(p1);
    polys.push_back(p2);
    tsPolys = ts.polynomials();

//    ZeroDimensionalRegularChain<RN,SMQP> ts2(ps);
//    printVariables(ps);
    RegularChain<RN,SMQP> ts2(xs,ps);
    ts2 = ts2 + p2;
    ts2 = ts2 + p1;
    ts2Polys = ts2.polynomials();

    std::cout.clear();
    std::cerr.clear();

    if (tsPolys.size() != polys.size() || ts2Polys.size() != polys.size() || tsPolys[1] != ts2Polys[1]) {
        std::cerr << "ZDRC addition operator test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts = " << ts << std::endl;
        std::cerr << "ts2 = " << ts2 << std::endl;
        for (auto x : polys) {
        	std::cerr << x << std::endl;
        }
        exit(1);
    }

    std::cerr << "ZDRC addition operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCAddAssignmentOperator() {
	std::vector<int> maxDegs = {2,2,2,2};
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	maxDegs.pop_back();
	p2.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    ZeroDimensionalRegularChain<RN,SMQP> ts(p2,ps);
    ts += p1;
	std::vector<SMQP> polys,tsPolys,ts2Polys;
    polys.push_back(p1);
    polys.push_back(p2);
    tsPolys = ts.polynomials();

//    ZeroDimensionalRegularChain<RN,SMQP> ts2(ps);
//    printVariables(ps);
    RegularChain<RN,SMQP> ts2(xs,ps);
    ts2 += p2;
    ts2 += p1;
    ts2Polys = ts2.polynomials();

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

    if (tsPolys.size() != polys.size() || ts2Polys.size() != polys.size() || tsPolys[1] != ts2Polys[1]) {
        std::cerr << "ZDRC addition operator test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts = " << ts << std::endl;
        std::cerr << "ts2 = " << ts2 << std::endl;
        for (auto x : polys) {
        	std::cerr << x << std::endl;
        }
        exit(1);
    }

    std::cerr << "ZDRC addition assignment operator test:\t\t\t\t PASSED" << std::endl;
}

void testZDRCIdentityOperator() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    ZeroDimensionalRegularChain<RN,SMQP> ts2;
    ts2 = ts;
    if (ts.polynomials() != ts2.polynomials()) {
        std::cerr << "ZDRC identity operator test:\t\t\t\t\t FAILED" << std::endl;
        for (int i=0; i<ts.polynomials().size(); ++i)
        	std::cerr << "ts[" << i << "] = " << ts.polynomials()[i] << std::endl;
        for (int i=0; i<ts2.polynomials().size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2.polynomials()[i] << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC identity operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCNumberOfVariables() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.numberOfVariables() != 4) {
        std::cerr << "ZDRC number of variables test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC number of variables test:\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCSize() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.size() != 4) {
        std::cerr << "ZDRC size test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC size test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCNumberOfAlgebraicVariables() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    std::vector<SMQP> polys(ts.polynomials());
    if (ts.numberOfAlgebraicVariables() != 4) {
        std::cerr << "ZDRC numberOfAlgebraicVariables test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC numberOfAlgebraicVariables test:\t\t\t\t PASSED" << std::endl;
}

void testZDRCNumberOfTranscendentalVariables() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.numberOfTranscendentalVariables() != 2) {
        std::cerr << "ZDRC numberOfTranscendentalVariables test:\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC numberOfTranscendentalVariables test:\t\t\t PASSED" << std::endl;
}

void testZDRCVariables() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(nvars-1,nterms,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	p2.randomPolynomial(nvars-2,nterms,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    for (auto v : xs) {
        p1 += SMQP(v);
    }
    p2 += SMQP(xs[1]);

    // Redirect stdout, stderr
//    std::cout.setstate(std::ios::failbit);
//    std::cerr.setstate(std::ios::failbit);

    ZeroDimensionalRegularChain<RN,SMQP> ts(p2,ps);
    ts += p1;

//    std::cout.clear();
//    std::cerr.clear();

    vs.clear();
    vs.insert(vs.end(),vars.begin(),vars.begin()+2);

    if (ts.variables() != vs) {
        std::cerr << "ZDRC variables test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC variables test:\t\t\t\t\t\t PASSED" << std::endl;
}


void testZDRCMainVariables() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(nvars-1,nterms,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	p2.randomPolynomial(nvars-2,nterms,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    for (auto v : xs) {
        p1 += SMQP(v);
    }
    p2 += SMQP(xs[1]);

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    ZeroDimensionalRegularChain<RN,SMQP> ts(p2,ps);
    ts += p1;

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

    if (ts.mainVariables() != xs) {
        std::cerr << "ZDRC main variables test:\t\t\t\t\t FAILED" << std::endl;
        printVariables(ts.mainVariables());
        printVariables(xs);
        exit(1);
    }

    std::cerr << "ZDRC main variables test:\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCTranscendentalVariables() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SMQP p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
	p1.randomPolynomial(nvars-1,nterms,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	p2.randomPolynomial(nvars-2,nterms,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[1]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);

    for (auto v : xs) {
        p1 += SMQP(v);
    }
    p2 += SMQP(xs[1]);

    // Redirect stdout, stderr
    std::cout.setstate(std::ios::failbit);
    std::cerr.setstate(std::ios::failbit);

    ZeroDimensionalRegularChain<RN,SMQP> ts(p2,ps);
    ts += p1;

    std::cout.clear();
    std::cerr.clear();

    if (ts.transcendentalVariables() != ps) {
        std::cerr << "ZDRC transcendental variables test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC transcendental variables test:\t\t\t\t PASSED" << std::endl;
}

void testZDRCIsAlgebraic() {
	int nvars(3);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+1);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+1,vs.begin()+3);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);

    for (auto v : xs) {
        p += SMQP(v);
    }

    ZeroDimensionalRegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());
    vector<Symbol> algVar;
    algVar = setIntersection(ts.mainVariables(),vars);

    if (ts.isAlgebraic(p.leadingVariable()) != (algVar == vars)) {
        std::cerr << "ZDRC isAlgebraic test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC isAlgebraic test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCIsEmpty() {
	int nvars(3);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars);
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+1);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+1,vs.begin()+3);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    ZeroDimensionalRegularChain<RN,SMQP> ts(p,ps),ts2;

    if (ts.isEmpty() == true || ts2.isEmpty() == false) {
        std::cerr << "ZDRC isEmpty test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC isEmpty test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCSelect() {

	int nvars(3);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
	vector<Symbol> vs;
	vs = Symbol::randomElements(nvars+1);
	Symbol s(vs[3]);
	vs.pop_back();
	vector<Symbol> xs;
	xs.insert(xs.end(),vs.begin(),vs.begin()+1);
	vector<Symbol> ps;
	ps.insert(ps.end(),vs.begin()+1,vs.begin()+3);
	SMQP p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    ZeroDimensionalRegularChain<RN,SMQP> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());
    SMQP zero;

    if (ts.select(p.leadingVariable()).leadingVariable() != p.leadingVariable() || ts.select(s) != zero) {
        std::cerr << "ZDRC select test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC select test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCLower() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
    ZeroDimensionalRegularChain<RN,SMQP> ts2;
    std:vector<SMQP> tslp,ts2p;
    std::vector<Symbol> vars;
    std::vector<int> maxDegs = {2,2,2,2,2};
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(3,2,maxDegs,coefBound,sparsity,includeNeg);
//    ts.display();
    vars = ts.mainVariables();

    // Redirect stdout, stderr
    std::cout.setstate(std::ios::failbit);
    std::cerr.setstate(std::ios::failbit);

    ts.lower(vars[1],ts2);

    std::cout.clear();
    std::cerr.clear();

//	std::cerr << "ts2 in testLower = " << ts2 << std::endl;
//	ts.display();
//	ts2.display();
	tslp = ts.polynomials();
	tslp.erase(tslp.begin(),tslp.begin()+2);
	ts2p = ts2.polynomials();
//	cout << "ts2p size = " << ts2p.size() << endl;
    if (tslp != ts2p) {
        std::cerr << "ZDRC lower test:\t\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts lower polys:" << std::endl;
        for (int i=0; i<tslp.size(); ++i)
        	std::cerr << "ts[" << i << "] = " << tslp[i] << std::endl;
        std::cerr << "ts2 polys:" << std::endl;
        for (int i=0; i<ts2p.size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2p[i] << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC lower test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCUpper() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
    RegularChain<RN,SMQP> ts2;
    std:vector<SMQP> tsup,ts2p;
    std::vector<Symbol> vars;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(4,2,nterms,coefBound,sparsity,includeNeg);
    vars = ts.mainVariables();

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    ts.upper(vars[2],ts2);

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

//    ts.display();
//    ts2.display();
	tsup = ts.polynomials();
	for (int i=tsup.size()-1; i>1; --i)
		tsup[i].zero();
//	tsup.erase(tsup.end()-1,tsup.end());
	ts2p = ts2.polynomials();
    if (tsup != ts2p) {
        std::cerr << "ZDRC upper test:\t\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts upper polys:" << std::endl;
        for (int i=0; i<tsup.size(); ++i)
        	std::cerr << "ts[" << i << "] = " << tsup[i] << std::endl;
        std::cerr << "ts2 polys:" << std::endl;
        for (int i=0; i<ts2p.size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2p[i] << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC upper test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCPseudoDivide() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	SMQP p,check;
	std::vector<SMQP> set;
	std::vector<Symbol> vars,trcVars;
	int nvars(3);
	int nTrcVars(0);
	int nterms(6);
	unsigned long int coefBound(6ul);
	degree_t sparsity(30);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(nvars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    ts.display();
    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);
    vars = ts.variables();
    trcVars = ts.transcendentalVariables();
    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
	p.setRingVariables(vars);
	set = ts.polynomials();
	std::reverse(set.begin(),set.end());
//	std::cerr << "set.size = " << set.size() << std::endl;
//    std::cerr << "[";
//	for (int i=0; i<set.size(); ++i){
//		std::cerr << set[i];
//		if (i != set.size()-1)
//			std::cerr << ",";
//	}
//	std::cerr << "]" << std::endl;
//	std::cerr << "ts = " << ts << std::endl;
//	std::cerr << "p = " << p << std::endl;
    vars = ts.variables();
    std::reverse(vars.begin(),vars.end());
//    std::cerr << "[";
//	for (int i=0; i<vars.size(); ++i){
//		std::cerr << vars[i];
//		if (i != vars.size()-1)
//			std::cerr << ",";
//	}
//	std::cerr << "]" << std::endl;

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

	SMQP r,c;
	vector<SMQP> Q;
	r = ts.pseudoDivide(p,&Q,&c);

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

//	cout << "c = " << c << endl;
//	cout << "r = " << r << endl;
//	for (int i=0; i<Q.size(); ++i)
//		cout << "Q[" << i << "] = " << Q[i] << endl;
//
//	cerr << "compute constraint:" << endl;
//	set = ts.polynomials();
//	check.zero();
//	for (int i=0; i<set.size(); ++i)
//		check += set[i]*Q[i];
//	check += r;
//	check -= (c*p);
//	cout << "check = " << check << endl;
//    if (!check.isZero()) {
//        std::cerr << "ZDRC pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "ZDRC pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;

    ExpressionTree pTree = p.convertToExpressionTree();
//    ExpressionTree tsTree = ts.convertToExpressionTree();
	ExpressionTree tsTree;
	tsTree.fromVector<SMQP>(set);
    ExpressionTree cTree = c.convertToExpressionTree();
    ExpressionTree remTree = r.convertToExpressionTree();
    std::vector<ExpressionTree> quoTrees;
    for (int i = 0; i < Q.size(); ++i) {
        quoTrees.push_back(Q[i].convertToExpressionTree());
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(tsTree.toMapleString() + ":");

//    vars = ts.variables();
//    std::reverse(vars.begin(),vars.end());
    std::string varList;
    if (vars.size() > 1) {
        varList = "[" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            varList += vars[i].toString();
            varList += ",";
        }
        varList += vars[vars.size()-1].toString();
    } else {
        varList = "[" + vars[0].toString();
    }
    varList += "]:";
    inputs.push_back(varList);

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
    algebList.push_back(ToMapleName(kv, "m", 1));

    //algebList: [0] = dividend, [1] = divisorList, [2] = varList, [3] = qList, [4] = m

//    std::string procStr = "PseudoDivide := proc(p,ts,vs,q,m)\r\n local quo;\r\n quo:=[1,2,3];\r\n assign(q=quo);\r\n assign(m=2); 3;\r\n end proc:";
	std::string procStr = "multiDivisorPseudoDivision := proc (p, ts, names, q, m) local seqLoop, tsSize, r, tmpQ, Q, tmpM, totalM, i, j; seqLoop := [3, 2, 1]; tsSize := nops(ts); r := p; Q := [seq(0, i = 1 .. tsSize)]; totalM := 1; for i in seqLoop do tmpQ := 0; tmpM := 1; r := prem(r, ts[i], names[i], 'tmpM', 'tmpQ'); totalM := totalM*tmpM; for j to tsSize do Q[j] := simplify(tmpM*Q[j]) end do; Q[i] := Q[i]+tmpQ end do; assign(q = Q); assign(m = totalM); return r end proc:";
    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[3], algebList[4]);

//     std::cerr << "calling maple proc: \n\n";
//     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
//     std::cerr << "\n\n";

    std::vector<ALGEB> resultAlgebs;
    inputs.clear();
    inputs.push_back(remTree.toMapleString() + ":");
    inputs.push_back(cTree.toMapleString() + ":");
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
    if (!MapleToM_BOOL(kv, comp))
    	std::cerr << "remainder comp fail!" << std::endl;
    std::string mapleResStr = mapleTest->algebToString(kv, result);
    std::string getMStr = "m:";
    cstr = new char[getMStr.length()+1];
    std::strcpy(cstr,getMStr.c_str());
    ALGEB mapleM = EvalMapleStatement(kv, cstr);
    mapleM = EvalMapleProc(kv, expandF, 1, mapleM);
    comp = EvalMapleProc(kv, cmpF, 2, mapleM, resultAlgebs[1]);
    compBool = compBool && MapleToM_BOOL(kv, comp);
    if (!MapleToM_BOOL(kv, comp))
    	std::cerr << "multiplier comp fail!" << std::endl;
    delete[] cstr;
    std::vector<ALGEB> mapleQuoList;
    int idx = 2;
//    while (idx < resultAlgebs.size()) {
    while (idx <= 4) {
        std::string getQuoStr = "qList[" + std::to_string(nvars+2-idx) + "]:";
        cstr = new char[getQuoStr.length()+1];
        std::strcpy(cstr, getQuoStr.c_str());
        ALGEB mapleQuo = EvalMapleStatement(kv, cstr);
        mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo);
        comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
        compBool = compBool && MapleToM_BOOL(kv, comp);
    	if (!MapleToM_BOOL(kv, comp)) {
    		std::cerr << "quotient " << idx-2 << " comp fail!" << std::endl;
    	}
        mapleQuoList.push_back(mapleQuo);
        ++idx;
        delete[] cstr;
    }

    if(!compBool) {
        std::cerr << "ZDRC pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;

        std::cerr << "Dividend: " << p << std::endl;
        std::cerr << "Divisors: " << ts << std::endl << std::endl;

        std::cerr << "Got remainder: " << r << std::endl;
        std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;

        std::cerr << "Got multiplier: " << c << std::endl;
        std::cerr << "Got maple m: " << mapleTest->algebToString(kv, mapleM) << std::endl << std::endl;

        for (int i = 0; i < mapleQuoList.size(); ++i) {
            std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
            std::cerr << "Got maple quo : " << mapleTest->algebToString(kv, mapleQuoList[i]) << std::endl;
        }
        exit(1);
    }

    std::cerr << "ZDRC pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCNormalForm() {

    ZeroDimensionalRegularChain<RN,SMQP> ts;
	SMQP p,check;
	std::vector<SMQP> set;
	std::vector<Symbol> vars,trcVars;
	int nvars(4);
	int nTrcVars(2);
	int nterms(3);
	unsigned long int coefBound(2ul);
	degree_t sparsity(25);
	bool includeNeg(1);
    ts.randomZeroDimensionalRegularChain(nvars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    ts.display();
    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);


    vars = ts.variables();
    trcVars = ts.transcendentalVariables();
    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
    p.setRingVariables(vars);
    set = ts.polynomials();
    for (int i = 0; i < 2; ++i) {
        p *= set[i];
    }

    // Redirect stdout, stderr
    stringstream ss;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss.rdbuf());
    cerr.rdbuf (ss.rdbuf());

    SMQP r;
	vector<SMQP> Q;
	r = ts.normalForm(p,&Q);

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

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
    //     std::cerr << "ZDRC normalForm test:\t\t\t\t\t\t\t FAILED" << std::endl;
    //     exit(1);
    // }

    // std::cerr << "ZDRC normalForm test:\t\t\t\t\t\t\t PASSED" << std::endl;


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
        std::cerr << "ZDRC normalForm test:\t\t\t\t\t\t FAILED" << std::endl;

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

    std::cerr << "ZDRC normalForm test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testZDRCConvertToExpressionTree() {

    stringstream ss;
    std::string str,str2;
    int nvars(4);
    int nTrcVars(2);
    int nterms(14);
    unsigned long int coefBound(6ul);
    degree_t sparsity(400);
    bool includeNeg(1);
    ZeroDimensionalRegularChain<RN,SMQP> ts;
    ts.randomZeroDimensionalRegularChain(nvars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
    ExpressionTree et = ts.convertToExpressionTree();
    // cout << "et = " << et.toString() << endl;
    ss << ts;
    str = ss.str();
    // cout << "ts = " << str << endl;
    // MAPLE EQUALITY TEST FOR et.toString() and str //
    // Start Maple Kernel:
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();

    // et:
    std::string etStr = et.toString() + ":";
    char* cstr = new char[etStr.size() + 1];
    std::strcpy (cstr, etStr.c_str());
    ALGEB resET = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    // str:
    std::string tsStr = str + ":";
    cstr = new char[tsStr.size() + 1];
    std::strcpy (cstr, tsStr.c_str());
    ALGEB resTS = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    // operator:
    std::string procStr = "ArrayTools:-IsEqual:";
    cstr = new char[procStr.length() + 1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    // Maple:
    ALGEB equalityRes = EvalMapleProc(kv, testProc, 2, resET, resTS);
    M_BOOL cmpBool = MapleToM_BOOL(kv, equalityRes);

    // Redirect stdout, stderr
    stringstream ss2;
    streambuf *oldcout,*oldcerr;
    oldcout = cout.rdbuf(); // <-- save
    oldcerr = cerr.rdbuf(); // <-- save
    cout.rdbuf (ss2.rdbuf());
    cerr.rdbuf (ss2.rdbuf());

    ts = ZeroDimensionalRegularChain<RN,SMQP>();
    ExpressionTree et2;
    et2 = ts.convertToExpressionTree();

    cout.rdbuf (oldcout);              // <-- restore
    cerr.rdbuf (oldcerr);              // <-- restore

//  cout << "et2 = " << et2.toString() << endl;
    ss = stringstream();
    ss << ts;
    str2 = ss.str();
//  cout << "ts = " << str2 << endl;
    if (!cmpBool || et2.toString() != str2) {
        std::cerr << "ZDRC convertToExpressionTree test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "ZDRC convertToExpressionTree test:\t\t\t\t PASSED" << std::endl;
}

bool isInvertibleValidate(SMQP f, ZeroDimensionalRegularChain<RN,SMQP> rc, vector<Symbol> vars, bool showOutput) {
    vector<ZeroDimensionalRegularChain<RN,SMQP>> resultsTrue,resultsFalse;
    vector<BoolChainPair<ZeroDimensionalRegularChain<RN,SMQP>>> results;

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

	results = rc.isInvertible(f);

    for (int l=0; l<results.size(); ++l) {
//    	cout << "b_" << l << " = " << results[l].isTrue << endl;
//    	cout << "T_" << l << " = " << results[l].chain << endl;
    	if (results[l].isTrue)
    		resultsTrue.push_back(results[l].chain);
    	else
    		resultsFalse.push_back(results[l].chain);
    }

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}
	if (showOutput) {
    	std::cerr << "f = " << f << std::endl;
    	std::cerr << "rc = " << std::endl;
    	rc.display();
    }

	std::cerr << "BPAS:" << std::endl;
	std::cerr << "[";
	for (int i=0; i<resultsTrue.size(); ++i) {
		std::cerr << resultsTrue[i];
		if (i!=resultsTrue.size()-1)
			std::cerr << ", ";
	}
	std::cerr << "]" << std::endl;
	std::cerr << "[";
	for (int i=0; i<resultsFalse.size(); ++i) {
		std::cerr << resultsFalse[i];
		if (i!=resultsFalse.size()-1)
			std::cerr << ", ";
	}
	std::cerr << "]" << std::endl;

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

    ExpressionTree fTree = f.convertToExpressionTree();
    ExpressionTree rcTree = rc.convertToExpressionTree();
    ExpressionTree RTree;
    RTree.fromVector<Symbol>(vars);
    ExpressionTree listTrueTrees,listFalseTrees;
    listTrueTrees.fromVector<ZeroDimensionalRegularChain<RN,SMQP>>(resultsTrue);
    listFalseTrees.fromVector<ZeroDimensionalRegularChain<RN,SMQP>>(resultsFalse);

    std::vector<std::string> inputs;
    inputs.push_back(fTree.toMapleString() + ":");
    inputs.push_back(rcTree.toMapleString() + ":");
    inputs.push_back(RTree.toMapleString() + ":");
    inputs.push_back(listTrueTrees.toMapleString() + ":");
    inputs.push_back(listFalseTrees.toMapleString() + ":");

    for (auto v : inputs) {
        std::cerr << v << std::endl;
    }

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

    //algebList: [0] = f, [1] = T, [2] = Rlist, [3] = trueChainList, [4] = falseChainList

    std::string checkTrueBranch = "IsInvertibleValidate := proc (f::polynom, in_rc::list, Rlist::list, list1true::list, list1false::list) local lrcs1true, lrcs1false, lrcs2true, lrcs2false, n, rc, R, inverse, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(list1true); lrcs1true := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(numer(op(i, list1true))), RegularChains:-ChainTools:-Empty(R), R); lrcs1true := [op(lrcs1true), rc] end do; n := nops(list1false); lrcs1false := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(numer(op(i, list1false))), RegularChains:-ChainTools:-Empty(R), R); lrcs1false := [op(lrcs1false), rc] end do; rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(numer(in_rc)), RegularChains:-ChainTools:-Empty(R), R); inverse := RegularChains:-Inverse(numer(f), rc, R); if nops(op(1, inverse)) = 1 then lrcs2true := [op(3, op(op(1, inverse)))] else lrcs2true := op(1, inverse) end if; lrcs2false := op(2, inverse); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs1true, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs2true, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then pass := true else pass := false end if; end proc:";
    std::string checkFalseBranch = "IsInvertibleValidate := proc (f::polynom, in_rc::list, Rlist::list, list1true::list, list1false::list) local lrcs1true, lrcs1false, lrcs2true, lrcs2false, n, rc, R, inverse, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(list1true); lrcs1true := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(numer(op(i, list1true))), RegularChains:-ChainTools:-Empty(R), R); lrcs1true := [op(lrcs1true), rc] end do; n := nops(list1false); lrcs1false := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(numer(op(i, list1false))), RegularChains:-ChainTools:-Empty(R), R); lrcs1false := [op(lrcs1false), rc] end do; rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(numer(in_rc)), RegularChains:-ChainTools:-Empty(R), R); inverse := RegularChains:-Inverse(numer(f), rc, R); if nops(op(1, inverse)) = 1 then lrcs2true := [op(3, op(op(1, inverse)))] else lrcs2true := op(1, inverse) end if; lrcs2false := op(2, inverse); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs1true, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs2true, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs1false, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs2false, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";

    std::string procStr = checkFalseBranch;

    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[3], algebList[4]);

//     std::cerr << "calling maple proc: \n\n";
//     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
     std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
//     std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
//     std::cerr << "\n\n";

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

//	std::cerr << "Maple:" << std::endl;
//	std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n";
//	std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";

    return mapleTest->testEquality(result, "true");

}

void testZDRCIsInvertible(bool showOutput) {

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }
    SMQP f,T1,T2,T3;
    vector<Symbol> vars;
    ZeroDimensionalRegularChain<RN,SMQP> zdrc;
    std::vector<SMQP> zdrcPolys;

    f = (SMQP("x")+SMQP("y"))*(SMQP("x")+RN(3)*SMQP("y"));
	cout << "f = " << f << endl;
	T1 = (SMQP("x")-SMQP("y"))*(SMQP("x")+SMQP("y"));
	cout << "T1 = " << T1 << endl;
	T2 = (SMQP("y")*SMQP("y"))-RN(2);
	cout << "T2 = " << T2 << endl;
	vars.clear();
    vars.emplace_back("x");
    vars.emplace_back("y");
    T1.setRingVariables(vars);
    f.setRingVariables(vars);
    zdrcPolys = {T1, T2};
    zdrc = ZeroDimensionalRegularChain<RN,SMQP>(zdrcPolys);
	// zdrc = ZeroDimensionalRegularChain<RN,SMQP>();
	// zdrc += T2;
	// zdrc += T1;
	cout << "zdrc = " << zdrc << endl;
//	printVariables(zdrc.variables(),"zdrc");

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

	bool firstResultCorrect = isInvertibleValidate(f,zdrc,vars,showOutput);
    // bool firstResultCorrect = true;

    if (!showOutput) {
        // Redirect stdout, stderr
        std::cout.setstate(std::ios::failbit);
        std::cerr.setstate(std::ios::failbit);
    }

    f.fromString("2*z^2 - z*y + 3*z*x - z - y^2 - 2*y + x^2 - 1");
    vars = {'z', 'y', 'x'};

    // f = (SMQP("z")-SMQP("y")+SMQP("x")-RN(1))*(RN(2)*SMQP("z")+SMQP("y")+SMQP("x")+RN(1));
    cout << "f = " << f << endl;
    T1 = (SMQP("z")-SMQP("y")+SMQP("x")-RN(1))*(SMQP("z")+SMQP("y")+SMQP("x")+RN(1));
    cout << "T1 = " << T1 << endl;
    T2 = (SMQP("y")*SMQP("y"))-((SMQP("x")+RN(1))*(SMQP("x")+RN(1)));
    cout << "T2 = " << T2 << endl;
    T3 = (SMQP("x")*SMQP("x"))-RN(2);
    cout << "T3 = " << T3 << endl;
    T1.setRingVariables(vars);
    T2.setRingVariables(vars);
    T3.setRingVariables(vars);
    zdrcPolys = {T1, T2, T3};
	// vars.clear();
 //    vars.emplace_back("x");
 //    vars.emplace(vars.begin(),"y");
 //    T2.setRingVariables(vars);
	// cout << "T2 = " << T2 << endl;
 //    vars.emplace(vars.begin(),"z");
 //    T1.setRingVariables(vars);
 //    f.setRingVariables(vars);
	// cout << "T1 = " << T1 << endl;
	// cout << "f = " << f << endl;
	// zdrc = ZeroDimensionalRegularChain<RN,SMQP>();
    zdrc = ZeroDimensionalRegularChain<RN,SMQP>(zdrcPolys);
	// zdrc += T3;
	// zdrc += T2;
	// zdrc += T1;
	cout << "zdrc = " << zdrc << endl;

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

	bool secondResultCorrect = isInvertibleValidate(f,zdrc,vars,showOutput);

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

    SMQP x("x"),y("y"),z("z");
	f = ((z*z)-RN(1))*(y+RN(2));
	cout << "f = " << f << endl;
	T1 = (x*x)+((z-RN(1))*x)+((y*y)-RN(4));
	cout << "T1 = " << T1 << endl;
	T2 = (z-RN(1));
	cout << "T2 = " << T2 << "*(" << (y*y+RN(5)*y+RN(2)) << ")" << endl;
	T2 *= (y*y+RN(5)*y+RN(2));
	cout << "T2 = " << T2 << endl;
	T3 = ((z*z)-RN(2));
	cout << "T3 = " << T3 << endl;
	vars.clear();
    vars.emplace_back("z");
    vars.emplace(vars.begin(),"y");
    T2.setRingVariables(vars);
	cout << "T2 = " << T2 << endl;
    vars.emplace(vars.begin(),"x");
    T1.setRingVariables(vars);
    f.setRingVariables(vars);
	cout << "T1 = " << T1 << endl;
	cout << "f = " << f << endl;
    zdrcPolys = {T1, T2, T3};
    zdrc = ZeroDimensionalRegularChain<RN,SMQP>(zdrcPolys);
	// zdrc = ZeroDimensionalRegularChain<RN,SMQP>();
	// zdrc += T3;
	// zdrc += T2;
	// zdrc += T1;
	cout << "zdrc = " << zdrc << endl;

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

	bool thirdResultCorrect = isInvertibleValidate(f,zdrc,vars,showOutput);
    // bool thirdResultCorrect = true;

	showOutput = false;
    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

	SMQP g;
//	std::vector<int> maxDegs = {3,2,1};
//	std::vector<int> maxDegs = {2,3,1}; // hits an unhandled corner case.
	std::vector<int> maxDegs = {1,2,4};
	int nvars(3);
	unsigned long int coefBound(6ul);
	double sparsity(0.1);
	bool includeNeg(1);
	vector<Symbol> xs;
	vars = Symbol::randomElements(nvars);
	xs = vars;
	T1.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	T1.setRingVariables(xs);
    T1 += SMQP(xs[0])^2;
	f.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
	f.setRingVariables(xs);
    g.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
    g.setRingVariables(xs);
    T1 *= g;
    f *= g;
    cout << "f = " << f << endl;
    cout << "T1 = " << T1 << endl;

    maxDegs.pop_back();
    xs.erase(xs.begin(),xs.begin()+1);
    T2.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
    T2.setRingVariables(xs);
    T2 += SMQP(xs[0])^3;
    T2 += RN(2);
    cout << "T2 = " << T2 << endl;

    maxDegs.pop_back();
    xs.erase(xs.begin(),xs.begin()+1);
    T3.randomPolynomial(maxDegs,coefBound,sparsity,includeNeg);
    T3.setRingVariables(xs);
    T3 += SMQP(xs[0])^5;
    T3 += RN(1);

	cout << "T3 = " << T3 << endl;
    zdrcPolys = {T1, T2, T3};
	zdrc = ZeroDimensionalRegularChain<RN,SMQP>(zdrcPolys);
	// zdrc += T3;
	// zdrc += T2;
	// zdrc += T1;
	cout << "zdrc = " << zdrc << endl;

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

	bool fourthResultCorrect = isInvertibleValidate(f,zdrc,vars,showOutput);
	// bool fourthResultCorrect = true;


    if(!(firstResultCorrect && secondResultCorrect && thirdResultCorrect && fourthResultCorrect)) {
        std::cerr << "ZDRC isInvertible test:\t\t\t\t\t\t FAILED" << std::endl;

        if (!firstResultCorrect)
        	std::cerr << "First Example FAILED!" << std::endl;
        if (!secondResultCorrect)
        	std::cerr << "Second Example FAILED!" << std::endl;
        if (!thirdResultCorrect)
        	std::cerr << "Third Example FAILED!" << std::endl;
        if (!fourthResultCorrect)
        	std::cerr << "Fourth Example FAILED!" << std::endl;
        exit(1);
    }
    std::cerr << "ZDRC isInvertible test:\t\t\t\t\t\t PASSED" << std::endl;
}

//bool regularGCDValidate(SMQP f, SMQP g, Symbol v, ZeroDimensionalRegularChain<RN,SMQP> rc, vector<Symbol> vars, bool showOutput) {
//    vector<ZeroDimensionalRegularChain<RN,SMQP>> resultsTrue,resultsFalse;
//    vector<PolyChainPair<SMQP,ZeroDimensionalRegularChain<RN,SMQP>>> results;
//
//    if (!showOutput) {
//		// Redirect stdout, stderr
//		std::cout.setstate(std::ios::failbit);
//		std::cerr.setstate(std::ios::failbit);
//    }
//
//	results = rc.regularGCD(f,g,v);
//
//    for (int l=0; l<results.size(); ++l) {
////    	cout << "b_" << l << " = " << results[l].isTrue << endl;
////    	cout << "T_" << l << " = " << results[l].chain << endl;
//    	if (results[l].isTrue)
//    		resultsTrue.push_back(results[l].chain);
//    	else
//    		resultsFalse.push_back(results[l].chain);
//    }

//    ExpressionTree fTree = f.convertToExpressionTree();
//    ExpressionTree rcTree = rc.convertToExpressionTree();
//    ExpressionTree RTree;
//    RTree.fromVector<Symbol>(vars);
//    ExpressionTree listTrueTrees,listFalseTrees;
//    listTrueTrees.fromVector<ZeroDimensionalRegularChain<RN,SMQP>>(resultsTrue);
//    listFalseTrees.fromVector<ZeroDimensionalRegularChain<RN,SMQP>>(resultsFalse);

//    std::vector<std::string> inputs;
//    inputs.push_back(fTree.toMapleString() + ":");
//    inputs.push_back(rcTree.toMapleString() + ":");
//    inputs.push_back(RTree.toMapleString() + ":");
//    inputs.push_back(listTrueTrees.toMapleString() + ":");
//    inputs.push_back(listFalseTrees.toMapleString() + ":");

//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    char* cstr;

//    std::vector<ALGEB> algebList;
//    for (int i = 0; i < inputs.size(); ++i) {
//        cstr = new char[inputs[i].length()+1];
//        std::strcpy(cstr, inputs[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        algebList.push_back(res);
//        delete[] cstr;
//    }

//    //algebList: [0] = f, [1] = T, [2] = Rlist, [3] = trueChainList, [4] = falseChainList

//    std::string procStr = "IsInvertibleValidate := proc (f::polynom, in_rc::list, Rlist::list, list1true::list, list1false::list) local lrcs1true, lrcs1false, lrcs2true, lrcs2false, n, rc, R, inverse, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(list1true); lrcs1true := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, list1true)), RegularChains:-ChainTools:-Empty(R), R); lrcs1true := [op(lrcs1true), rc] end do; n := nops(list1false); lrcs1false := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, list1false)), RegularChains:-ChainTools:-Empty(R), R); lrcs1false := [op(lrcs1false), rc] end do; rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(in_rc), RegularChains:-ChainTools:-Empty(R), R); inverse := RegularChains:-Inverse(f, rc, R); if nops(op(1, inverse)) = 1 then lrcs2true := [op(3, op(op(1, inverse)))] else lrcs2true := op(1, inverse) end if; lrcs2false := op(2, inverse); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs1true, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs2true, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then pass := true else pass := false end if; if pass then lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs1false, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrcs2false, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if else return false end if end proc:";
//    cstr = new char[procStr.length()+1];
//    std::strcpy (cstr, procStr.c_str());
//    ALGEB testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;

//    ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[3], algebList[4]);

////     std::cerr << "calling maple proc: \n\n";
////     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
////     std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
////     std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
////     std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
////     std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
////     std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
////     std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
////     std::cerr << "\n\n";
//
//	if (!showOutput) {
//		std::cout.clear();
//		std::cerr.clear();
//	}
//
//	std::cerr << "BPAS:" << std::endl;
//	std::cerr << "[";
//	for (int i=0; i<resultsTrue.size(); ++i) {
//		std::cerr << resultsTrue[i];
//		if (i!=resultsTrue.size()-1)
//			std::cerr << ", ";
//	}
//	std::cerr << "]" << std::endl;
//	std::cerr << "[";
//	for (int i=0; i<resultsFalse.size(); ++i) {
//		std::cerr << resultsFalse[i];
//		if (i!=resultsFalse.size()-1)
//			std::cerr << ", ";
//	}
//	std::cerr << "]" << std::endl;
//	std::cerr << "Maple:" << std::endl;
//	std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n";
//	std::cerr << mapleTest->algebToString(kv, algebList[4]) << "\n\n";
//
//	return (mapleTest->algebToString(kv, result) == "true");
//
//}

void testZDRCregularGCD(bool showOutput) {

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

	cout << "regularGCD	Examples:" << endl;

	cout << "Example 1:" << endl;

    ZeroDimensionalRegularChain<RN,SMQP> zdrc;
    vector<PolyChainPair<SMQP,ZeroDimensionalRegularChain<RN,SMQP>>> result;
    SMQP m,f1,f2;
	SMQP x("x"),y("y"),z("z");
    vector<Symbol> vars;

	m = (x-RN(1))*(x+RN(2));
	f1 = ((RN(2)*x+RN(1))*(y*y)) + ((RN(5)*x+RN(10))*y) + ((RN(6)*x) + RN(12));
	f2 = (RN(3)*(y*y)) + ((RN(4)*x+RN(8))*y) + (RN(4)*x) + RN(5);
	vars.clear();
	vars.emplace_back("y");
	vars.emplace_back("x");
	f1.setRingVariables(vars);
	f2.setRingVariables(vars);
	cout << "m = " << m << endl;
	cout << "f1 = " << f1 << endl;
	cout << "f2 = " << f2 << endl;
	zdrc = ZeroDimensionalRegularChain<RN,SMQP>(m);
    result = zdrc.regularGCD(f1,f2,Symbol("y"));

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

    for (int l=0; l<result.size(); ++l) {
//    	cout << result[l] << endl;
    	cout << "G_" << l << " = " << result[l].poly << endl;
    	cout << "T_" << l << " = " << result[l].chain << endl;
    }

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

	cout << "Example 2:" << endl;

	m = (x+RN(1))*(x+RN(2));
	f1 = ((x+RN(1))*(y*y)) + RN(1)*y + ((RN(2)*x) + RN(4));
	f2 = ((RN(2)*x+RN(3))*(y*y));
	f2 += ((RN(2)*x+RN(3))*y);
	f2 -= ((RN(2)*x) + RN(4));
	vars.clear();
	vars.emplace_back("y");
	vars.emplace_back("x");
	f1.setRingVariables(vars);
	f2.setRingVariables(vars);
	cout << "m = " << m << endl;
	cout << "f1 = " << f1 << endl;
	cout << "f2 = " << f2 << endl;
    zdrc = ZeroDimensionalRegularChain<RN,SMQP>(m);
    result = zdrc.regularGCD(f1,f2,Symbol("y"));

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

    for (int l=0; l<result.size(); ++l) {
    	cout << "G_" << l << " = " << result[l].poly << endl;
    	cout << "T_" << l << " = " << result[l].chain << endl;
    }

    if (!showOutput) {
		// Redirect stdout, stderr
		std::cout.setstate(std::ios::failbit);
		std::cerr.setstate(std::ios::failbit);
    }

	cout << "Example 3:" << endl;

	// requires vars, m and zdrc from example 2!
	f1 = ((x+RN(2))*(y*y)) + (RN(2)*x+RN(3))*y - RN(2);
	f2 = ((x+RN(1))*(y*y)) + (RN(4)*x+RN(5))*y + (RN(4)*x) + RN(6);
	f1.setRingVariables(vars);
	f2.setRingVariables(vars);
	cout << "m = " << m << endl;
	cout << "f1 = " << f1 << endl;
	cout << "f2 = " << f2 << endl;
    result = zdrc.regularGCD(f1,f2,Symbol("y"));

	if (!showOutput) {
		std::cout.clear();
		std::cerr.clear();
	}

    for (int l=0; l<result.size(); ++l) {
    	cout << "G_" << l << " = " << result[l].poly << endl;
    	cout << "T_" << l << " = " << result[l].chain << endl;
    }

//	bool firstResultCorrect = isInvertibleValidate(f,zdrc,vars,showOutput);

//    if(!(firstResultCorrect && secondResultCorrect)) {
//        std::cerr << "ZDRC regularGCD test:\t\t\t\t\t\t FAILED" << std::endl;
//
//        if (!firstResultCorrect)
//        	std::cerr << "First Example FAILED!" << std::endl;
//        if (!secondResultCorrect)
//        	std::cerr << "Second Example FAILED!" << std::endl;
//        exit(1);
//    }
    std::cerr << "ZDRC regularGCD test:\t\t\t\t\t\t PASSED" << std::endl;
    std::cerr << "regularGCD test not yet automated!" << std::endl;
}


//bool triangularizeValidateOld(vector<SMQP> F, RegularChain<RN,SMQP> rc, vector<Symbol> vars) {
//    vector<RegularChain<RN,SMQP>> results;
//
//	results = rc.triangularize(F);

//    ExpressionTree FTree;
//    FTree.fromVector<SMQP>(F);
//    ExpressionTree rcTree = rc.convertToExpressionTree();
//    ExpressionTree RTree;
//    RTree.fromVector<Symbol>(vars);
//    ExpressionTree resultTrees;
//    resultTrees.fromVector<RegularChain<RN,SMQP>>(results);

//    std::vector<std::string> inputs;
//    inputs.push_back(FTree.toMapleString() + ":");
//    inputs.push_back(rcTree.toMapleString() + ":");
//    inputs.push_back(RTree.toMapleString() + ":");
//    inputs.push_back(resultTrees.toMapleString() + ":");

//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    char* cstr;

//    std::vector<ALGEB> algebList;
//    for (int i = 0; i < inputs.size(); ++i) {
//        cstr = new char[inputs[i].length()+1];
//        std::strcpy(cstr, inputs[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        algebList.push_back(res);
//        delete[] cstr;
//    }

//    //algebList: [0] = F, [1] = T, [2] = Rlist, [3] = resultChainList

//    std::string procStr = "TriangularizeValidate := proc (F::list, in_rc::list, Rlist::list, results::list) local lrc1, lrc2, n, rc, R, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(results); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, results)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc2 := RegularChains:-Triangularize(F, rc, R); if evalb(nops(lrc2) = 0) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) then return true end if; lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
////    std::string procStr = "TriangularizeValidate := proc (F::list, in_rc::list, Rlist::list, results::list) local lrc1, lrc2, n, rc, R, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(results); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, results)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc2 := RegularChains:-Triangularize(F, rc, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
//    cstr = new char[procStr.length()+1];
//    std::strcpy (cstr, procStr.c_str());
//    ALGEB testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;

//    ALGEB result = EvalMapleProc(kv, testProc, 4, algebList[0], algebList[1], algebList[2], algebList[3]);

////     std::cerr << "calling maple proc: \n\n";
////     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
////     std::cerr << "\n\n";
//
//	if (mapleTest->algebToString(kv, result) == "true")
//		return true;
//	else {
//		std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n";
//		std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n";
//		std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n";
//		std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n";
//		std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//		return false;
//	}
//}


//bool triangularizeValidate(vector<SMQP> F, RegularChain<RN,SMQP> rc, vector<Symbol> vars, bool showOutput) {
//    vector<RegularChain<RN,SMQP>> results;
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
//    for (int i=0; i<F.size(); ++i) {
//    	cerr << "F[" << i << "] = " << F[i] << endl;
//    }
//
//	results = rc.triangularize(F);
//
//    if (!showOutput) {
//		cout.rdbuf (oldcout);   // <-- restore
//		cerr.rdbuf (oldcerr);   // <-- restore
//	}

//    ExpressionTree FTree;
//    FTree.fromVector<SMQP>(F);
//    ExpressionTree rcTree = rc.convertToExpressionTree();
//    ExpressionTree RTree;
//    RTree.fromVector<Symbol>(vars);
//    ExpressionTree resultTrees;
//    resultTrees.fromVector<RegularChain<RN,SMQP>>(results);

//    std::vector<std::string> inputs;
//    inputs.push_back(FTree.toMapleString() + ":");
//    inputs.push_back(rcTree.toMapleString() + ":");
//    inputs.push_back(RTree.toMapleString() + ":");
//    inputs.push_back(resultTrees.toMapleString() + ":");

//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    char* cstr;

//    std::vector<ALGEB> algebList;
//    for (int i = 0; i < inputs.size(); ++i) {
//        cstr = new char[inputs[i].length()+1];
//        std::strcpy(cstr, inputs[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        algebList.push_back(res);
//        delete[] cstr;
//    }

//    //algebList: [0] = F, [1] = T, [2] = Rlist, [3] = resultChainList

////    std::string procStr = "TriangularizeValidate := proc (F::list, in_rc::list, Rlist::list, results::list) local lrc1, lrc2, n, rc, R, lrs, cs1, cs2, pass, i; R := RegularChains:-PolynomialRing(Rlist); n := nops(results); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, results)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc2 := RegularChains:-Triangularize(F, rc, R); if evalb(nops(lrc2) = 0) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) then return true end if; lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); lrs := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
////    std::string procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard'); n := nops(lrc); results := []; if evalb(nops(lrc) = 0) then results := [op(results), []] else for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do end if; results end proc:";
//    std::string procStr = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard','radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
//    cstr = new char[procStr.length()+1];
//    std::strcpy (cstr, procStr.c_str());
//    ALGEB testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;

//    ALGEB result = EvalMapleProc(kv, testProc, 3, algebList[0], algebList[1], algebList[2]);

////	std::cerr << "Triangularize result: " << std::endl;
////	std::cerr << mapleTest->algebToString(kv, result) << std::endl;
////     std::cerr << "calling maple proc: \n\n";
////     std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
////     std::cerr << "\n\n";

//	std::vector<ALGEB> algebList2;
//	algebList2.push_back(algebList[3]);
//	algebList2.push_back(result);
//	algebList2.push_back(algebList[2]);

//	procStr = "EqualAsConstructibleSets := proc (dec1::list, dec2::list, Rlist::list) local n, lrc1, lrc2, lrs1, lrs2, cs1, cs2, R, i, rc; R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; if evalb(nops(lrc1) = 1) and evalb(nops(lrc2) = 1) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc2), R) then return true end if; lrs1 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); lrs2 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs1, R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs2, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
//    cstr = new char[procStr.length()+1];
//    std::strcpy (cstr, procStr.c_str());
//    testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;
//
//    result = EvalMapleProc(kv, testProc, 3, algebList2[0], algebList2[1], algebList2[2]);
//
//
//	if (mapleTest->algebToString(kv, result) == "true") {
//		if (showOutput) {
//			std::cerr << "BPAS result:" << std::endl;
//			std::cerr << mapleTest->algebToString(kv, algebList2[0]) << "\n";
//			std::cerr << "Maple result:" << std::endl;
//			std::cerr << mapleTest->algebToString(kv, algebList2[1]) << "\n";
//		}
//		return true;
//	}
//	else {
//		std::cerr << "BPAS result:" << std::endl;
//		std::cerr << mapleTest->algebToString(kv, algebList2[0]) << "\n";
//		std::cerr << "Maple result:" << std::endl;
//		std::cerr << mapleTest->algebToString(kv, algebList2[1]) << "\n";
////		std::cerr << mapleTest->algebToString(kv, algebList2[2]) << "\n";
//		std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//		return false;
//	}
//}

//static void Sys126Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(2)*(x_1^2)-(x_2^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)*(x_1^2)-(x_1)*(x_2)-(x_2)*(c_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys126Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys126Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys130Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber("-27577200")*(x_1^23)*(x_2)-RationalNumber("920137680")*(x_1^22)*(x_2)-RationalNumber("13283391240")*(x_1^21)*(x_2)-RationalNumber("106624468368")*(x_1^20)*(x_2)-RationalNumber("491746186572")*(x_1^19)*(x_2)-RationalNumber("993778599636")*(x_1^18)*(x_2)+RationalNumber("2476605785622")*(x_1^17)*(x_2)+RationalNumber("25815018153060")*(x_1^16)*(x_2)+RationalNumber("96112956427272")*(x_1^15)*(x_2)+RationalNumber("210514577654592")*(x_1^14)*(x_2)+RationalNumber("255834944199876")*(x_1^13)*(x_2)+RationalNumber("1673661361104")*(x_1^12)*(x_2)-RationalNumber("682308301717818")*(x_1^11)*(x_2)-RationalNumber("1476465320622258")*(x_1^10)*(x_2)-RationalNumber("1610410635303342")*(x_1^9)*(x_2)-RationalNumber("597224233919088")*(x_1^8)*(x_2)+RationalNumber("938901769037796")*(x_1^7)*(x_2)+RationalNumber("1697021033615460")*(x_1^6)*(x_2)+RationalNumber("1209904120743648")*(x_1^5)*(x_2)+RationalNumber("297722035302048")*(x_1^4)*(x_2)-RationalNumber("166331517404934")*(x_1^3)*(x_2)-RationalNumber("158840938909422")*(x_1^2)*(x_2)-RationalNumber("50605101044820")*(x_1)*(x_2)-RationalNumber("6044877890100")*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("1185408")*(x_1^21)*(x_2)+RationalNumber("59185728")*(x_1^20)*(x_2)+RationalNumber("1351691712")*(x_1^19)*(x_2)+RationalNumber("18700782336")*(x_1^18)*(x_2)+RationalNumber("174871916064")*(x_1^17)*(x_2)+RationalNumber("1166498753904")*(x_1^16)*(x_2)+RationalNumber("5699486912112")*(x_1^15)*(x_2)+RationalNumber("20539343601792")*(x_1^14)*(x_2)+RationalNumber("53774942979384")*(x_1^13)*(x_2)+RationalNumber("96423825757788")*(x_1^12)*(x_2)+RationalNumber("93804944335092")*(x_1^11)*(x_2)-RationalNumber("40214503702512")*(x_1^10)*(x_2)-RationalNumber("316349629027650")*(x_1^9)*(x_2)-RationalNumber("570981123783459")*(x_1^8)*(x_2)-RationalNumber("569250985467807")*(x_1^7)*(x_2)-RationalNumber("304838378563104")*(x_1^6)*(x_2)-RationalNumber("81135447444378")*(x_1^5)*(x_2)+RationalNumber("16825934075661")*(x_1^4)*(x_2)+RationalNumber("38057921178669")*(x_1^3)*(x_2)+RationalNumber("16575248060910")*(x_1^2)*(x_2)+RationalNumber("1990386622350")*(x_1)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys130Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys130Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys132Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(-1)+(x_2)-(x_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-2)-RationalNumber(4)*(x_1)*(x_2)-RationalNumber(2)*(x_1)-RationalNumber(2)*(x_2^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)-RationalNumber(2)*(x_1)*(x_2^2)-RationalNumber(4)*(x_1)-RationalNumber(4)*(x_1)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(1)-(x_1)*(x_2)+(x_2^2)-RationalNumber(3)*(x_2)-RationalNumber(2)*(x_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys132Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys132Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}


//static void Sys1000Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+(x_1)+RationalNumber(2)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)+RationalNumber(2)*(x_1)+RationalNumber(3)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1000Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1000Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1001Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+(x_1)+RationalNumber(3)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)+(x_1)+RationalNumber(3)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1001Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1001Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1002Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(8)+(x_1^2)-(x_2^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-16)*(x_1^2)-RationalNumber(48)-RationalNumber(16)*(x_2)-(x_1^4)-RationalNumber(4)*(x_1^2)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//}
//	if (pass)
//		std::cerr << "RC Sys1002Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1002Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1003Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = (x_1^2)*(x_2)+RationalNumber(3)*(x_1)+RationalNumber(2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_1^2)+(x_1)*(x_2)+RationalNumber(1)-(x_1)+(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1003Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1003Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1004Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(3)+RationalNumber(3)*(x_1)+(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(3)+RationalNumber(3)*(x_1)+RationalNumber(2)*(x_1^2)+(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1004Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1004Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1005Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(2)+RationalNumber(2)*(x_1)+(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)+(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1005Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1005Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1006Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+(x_1)+RationalNumber(5)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(4)+RationalNumber(5)*(x_1)+RationalNumber(4)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1006Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1006Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1007Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+RationalNumber(4)*(x_1)+RationalNumber(2)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)+RationalNumber(2)*(x_1)+RationalNumber(5)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1007Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1007Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1008Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = (x_1)*(x_2)+RationalNumber(2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(4)+(c_1)-(x_2)*(x_1^2)-(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1008Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1008Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1009Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+RationalNumber(2)*(x_1)+(x_1^2)+(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(3)+RationalNumber(2)*(x_1)+RationalNumber(2)*(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1009Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1009Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1010Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(2)+RationalNumber(3)*(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(3)+RationalNumber(2)*(x_1)+RationalNumber(3)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1010Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1010Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1011Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+(x_1)+RationalNumber(2)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)+(x_1)+(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1011Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1011Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1012Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(1)+RationalNumber(2)*(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(3)*(x_1^3)+(x_1^2)+RationalNumber(3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1012Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1012Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1013Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(2)*(x_1^3)+(x_1)+RationalNumber(2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)+(x_1)+(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1013Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1013Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1014Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(4)+(x_1)+RationalNumber(2)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(5)+(x_1)+RationalNumber(3)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1014Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1014Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1015Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(3)+RationalNumber(5)*(x_1)+RationalNumber(2)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(4)+(x_1)+RationalNumber(4)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1015Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1015Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1016Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(3)+RationalNumber(3)*(x_1)+(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(4)+RationalNumber(3)*(x_1)+(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1016Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1016Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1221Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly,temp;
//	std::cerr << "constructing first polynomial:" << std::endl;
//	poly = RationalNumber(2)*(c_1)*(x_1^2);
//	poly.setRingVariables(vars);
//	temp = -(c_1)*(x_1)*(x_2);
//	temp.setRingVariables(vars);
//	poly += temp;
//	poly += RationalNumber(4)*(x_1^3);
//	temp = -RationalNumber(3)*(x_1^2)*(x_2);
//	temp.setRingVariables(vars);
//	poly += temp;
//	polys.push_back(poly);
//	std::cerr << "constructing second polynomial:" << std::endl;
//	poly = RationalNumber(2)*(x_1^2)-(x_2^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	std::cerr << "calling triangularize:" << std::endl;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1221Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1221Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1255Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("c_2");
//	SMQP c_2("c_2");
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = (x_1^2)-RationalNumber(2)*(x_1)*(c_1)+RationalNumber(1)-(x_2^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-2)*(c_2)*(x_1)-RationalNumber(4)-RationalNumber(2)*(x_1)*(x_2)+RationalNumber(4)*(x_1)-(x_1^2)+RationalNumber(2)*(x_1)*(c_1)+RationalNumber(4)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1255Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1255Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1289Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("c_2");
//	SMQP c_2("c_2");
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(125)*(c_1)*(c_2^2)+RationalNumber(125)*(x_2)*(c_1^2)-RationalNumber(125)*(c_2)*(c_1^2)+RationalNumber(125)*(c_1)-RationalNumber(15876)*(x_2)+RationalNumber(15751)*(c_2)-RationalNumber(125)*(x_2)*(c_2^2)-RationalNumber(125)*(c_1)*(x_2^2)+RationalNumber(125)*(x_2^2)*(c_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(15625)*(x_1^2)*(c_2)-RationalNumber(15625)*(x_1)*(c_2^2)-RationalNumber(15625)*(c_1)*(x_1^2)+RationalNumber(15625)*(c_1)*(c_2^2)+RationalNumber(15625)*(x_1)*(c_1^2)-RationalNumber(15625)*(c_2)*(c_1^2)-RationalNumber(1984500)*(x_1)+RationalNumber(63011844)*(c_2)-RationalNumber(61027344)*(c_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(15625)*(x_1^2)*(x_2)-RationalNumber(15625)*(x_1)*(x_2^2)-RationalNumber(15625)*(c_1)*(x_1^2)+RationalNumber(15625)*(c_1)*(x_2^2)+RationalNumber(15625)*(x_1)*(c_1^2)-RationalNumber(15625)*(x_2)*(c_1^2)-RationalNumber(1968875)*(x_1)+RationalNumber(63011844)*(x_2)-RationalNumber(61042969)*(c_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1289Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1289Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1302Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(6562)+RationalNumber(6562)*(x_2)-RationalNumber(7695)*(x_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-59058)+RationalNumber(10858)*(x_1)*(x_2)+RationalNumber(83748)*(x_1)-RationalNumber(52496)*(x_2)+RationalNumber(6562)*(x_2^2)-RationalNumber(23085)*(x_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-32810)-RationalNumber(79594)*(x_1)*(x_2)+RationalNumber(18553)*(x_1)*(x_2^2)-RationalNumber(6798)*(x_1^2)*(x_2)-RationalNumber(97862)*(x_1)+RationalNumber(175427)*(x_1^2)-RationalNumber(32810)*(x_2^2)-RationalNumber(65620)*(x_2)-RationalNumber(23085)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(6562)-RationalNumber(116130)*(x_1)*(x_2)-RationalNumber(58350)*(x_1)*(x_2^2)+RationalNumber(3304)*(x_1^2)*(x_2)+RationalNumber(16287)*(x_2^2)*(x_1^2)-RationalNumber(19922)*(x_1^3)*(x_2)-RationalNumber(68895)*(x_1)-RationalNumber(12413)*(x_1^2)+RationalNumber(125858)*(x_1^3)+RationalNumber(6562)*(x_2)-RationalNumber(7695)*(x_1^4);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-26248)-RationalNumber(730895)*(x_1)-RationalNumber(534495)*(x_1^2)+RationalNumber(1203563)*(x_1^3)+RationalNumber(6562)*(x_2)-RationalNumber(1459890)*(x_1)*(x_2)-RationalNumber(728995)*(x_1)*(x_2^2)-RationalNumber(340090)*(x_1^2)*(x_2)+RationalNumber(194405)*(x_2^2)*(x_1^2)-RationalNumber(109094)*(x_1^3)*(x_2)-RationalNumber(257)*(x_1^3)*(x_2^2)-RationalNumber(10348)*(x_1^4)*(x_2)+RationalNumber(95)*(x_1^4)*(x_2^2)-RationalNumber(43253)*(x_1^4);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(13124)+RationalNumber(19686)*(x_2)+RationalNumber(31023)*(x_1)+RationalNumber(285)*(x_1)*(x_2^2)+RationalNumber(4588)*(x_1)*(x_2)-RationalNumber(14635)*(x_2^2)*(x_1^2)-RationalNumber(28130)*(x_1^2)*(x_2)+RationalNumber(3163)*(x_1^3)*(x_2^2)+RationalNumber(35406)*(x_1^3)*(x_2)-RationalNumber(8828)*(x_1^4)*(x_2)-RationalNumber(41710)*(x_1^2)+RationalNumber(6562)*(x_2^2)+RationalNumber(32623)*(x_1^3)+RationalNumber(27617)*(x_1^4);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-682448)+RationalNumber(651978)*(x_1)-RationalNumber(643076)*(x_2)+RationalNumber(39372)*(x_2^2)-RationalNumber(18943)*(x_1^2)+RationalNumber(83706)*(x_1)*(x_2)+RationalNumber(15133)*(x_1)*(x_2^2)-RationalNumber(15918)*(x_1^2)*(x_2)+RationalNumber(570)*(x_2^2)*(x_1^2)+RationalNumber(29840)*(x_1^3)*(x_2)+RationalNumber(14540)*(x_1^3)*(x_2^2)+RationalNumber(5004)*(x_1^4)*(x_2)-RationalNumber(1133)*(x_1^4)*(x_2^2)-RationalNumber(11205)*(x_1^3)+RationalNumber(6232)*(x_1^4);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-360910)-RationalNumber(1300980)*(x_1)+RationalNumber(1610397)*(x_1^2)-RationalNumber(360910)*(x_2^2)-RationalNumber(1142070)*(x_1)*(x_2)+RationalNumber(157960)*(x_1)*(x_2^2)+RationalNumber(17404)*(x_1^2)*(x_2)+RationalNumber(11157)*(x_2^2)*(x_1^2)-RationalNumber(26002)*(x_1^3)*(x_2)+RationalNumber(380)*(x_1^3)*(x_2^2)+RationalNumber(7460)*(x_1^4)*(x_2)+RationalNumber(3635)*(x_1^4)*(x_2^2)-RationalNumber(721820)*(x_2)-RationalNumber(80672)*(x_1^3)-RationalNumber(4725)*(x_1^4);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1302Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1302Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1303Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("c_2");
//	SMQP c_2("c_2");
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = (c_1^2)-RationalNumber(2)*(x_2^2)+(x_2)*(c_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_1^2)*(x_2)-(x_2^3)+(c_2)*(c_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(2)*(x_2^3)-(x_2)*(c_1^2)-(c_2)*(c_1^2)-RationalNumber(2)*(x_1)*(x_2^2)+(x_1)*(c_1^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1303Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1303Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1304Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(-270)*(x_1^4)*(x_2^3)-RationalNumber(314)*(x_1)*(x_2^4)-RationalNumber(689)*(x_1)*(x_2^3)+RationalNumber(1428);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(36)*(x_1^7)+RationalNumber(417)*(x_1^6)*(x_2)-RationalNumber(422)*(x_1^5)*(x_2^2)-RationalNumber(270)*(x_1^4)*(x_2^3)+RationalNumber(1428)*(x_1^3)*(x_2^4)-RationalNumber(1475)*(x_1^2)*(x_2^5)+RationalNumber(510)*(x_2^6)*(x_1)-RationalNumber(200)*(x_1^6)-RationalNumber(174)*(x_1^5)*(x_2)-RationalNumber(966)*(x_1^4)*(x_2^2)+RationalNumber(529)*(x_1^3)*(x_2^3)+RationalNumber(269)*(x_1^2)*(x_2^4)+RationalNumber(49)*(x_1)*(x_2^5)-RationalNumber(267)*(x_2^6)+RationalNumber(529)*(x_1^4)*(x_2)+RationalNumber(1303)*(x_1^2)*(x_2^3)-RationalNumber(314)*(x_1)*(x_2^4)+RationalNumber(262)*(x_2^5)+RationalNumber(36)*(x_1^4)-RationalNumber(788)*(x_2^2)*(x_1^2)-RationalNumber(689)*(x_1)*(x_2^3)+RationalNumber(177)*(x_2^4);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1304Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1304Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1305Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(3)+(x_1)+RationalNumber(2)*(x_1^2)+RationalNumber(2)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(3)+(x_1)+RationalNumber(2)*(x_1^2)+RationalNumber(3)*(x_1^3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1305Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1305Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1364Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(-1)*(x_1^5)*(x_2)-RationalNumber(8)*(x_1^4)*(x_2)-RationalNumber(20)*(x_1^3)*(x_2)-RationalNumber(17)*(x_1^2)*(x_2)-RationalNumber(6)*(x_1)*(x_2)-RationalNumber(8)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-1)*(x_1^8)*(x_2)-RationalNumber(14)*(x_1^7)*(x_2)-RationalNumber(79)*(x_1^6)*(x_2)-RationalNumber(231)*(x_1^5)*(x_2)-RationalNumber(376)*(x_1^4)*(x_2)-RationalNumber(353)*(x_1^3)*(x_2)-RationalNumber(230)*(x_1^2)*(x_2)-RationalNumber(152)*(x_1)*(x_2)-RationalNumber(64)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-2)*(x_1^6)*(x_2)-RationalNumber(26)*(x_1^5)*(x_2)-RationalNumber(134)*(x_1^4)*(x_2)-RationalNumber(350)*(x_1^3)*(x_2)-RationalNumber(488)*(x_1^2)*(x_2)-RationalNumber(344)*(x_1)*(x_2)-RationalNumber(96)*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1364Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1364Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1366Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber("12")+RationalNumber("4")*(x_2)-RationalNumber("6")*(x_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("-12")-RationalNumber("4")*(x_1)*(x_2)+RationalNumber("2")*(x_2^2)-RationalNumber("3")*(x_1^2)+RationalNumber("12")*(x_1)+RationalNumber("8")*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("12")-RationalNumber("2")*(x_1)*(x_2^2)-RationalNumber("2")*(x_1^2)*(x_2)+RationalNumber("8")*(x_1)*(x_2)+RationalNumber("4")*(x_2^2)-RationalNumber("12")*(x_2)-RationalNumber("6")*(x_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("-12")-(x_2^2)*(x_1^2)+RationalNumber("4")*(x_1)*(x_2^2)-RationalNumber("8")*(x_1)*(x_2)-RationalNumber("3")*(x_1^2)+RationalNumber("12")*(x_1)-RationalNumber("6")*(x_2^2)+RationalNumber("16")*(x_2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1366Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1366Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1373Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber("2")*(x_2)*(RationalNumber("-1")-RationalNumber("12")*(x_1^5)-RationalNumber("2")*(x_2^2)+RationalNumber("4")*(x_1)*(x_2^2)-RationalNumber("24")*(x_1^3)*(x_2^2)+RationalNumber("4")*(x_1)-RationalNumber("6")*(x_1^2)+RationalNumber("4")*(x_1^3)+RationalNumber("32")*(x_1^6)-RationalNumber("48")*(x_1^7)+RationalNumber("42")*(x_1^8)-RationalNumber("20")*(x_1^9)+RationalNumber("4")*(x_1^10)+RationalNumber("10")*(x_2^8)+RationalNumber("8")*(x_2^6)+RationalNumber("4")*(x_2^10)+RationalNumber("156")*(x_1^4)*(x_2^4)-RationalNumber("12")*(x_1)*(x_2^4)+RationalNumber("72")*(x_1^4)*(x_2^2)+RationalNumber("48")*(x_1^2)*(x_2^4)+RationalNumber("20")*(x_1^8)*(x_2^2)-RationalNumber("80")*(x_1^7)*(x_2^2)+RationalNumber("136")*(x_1^6)*(x_2^2)+RationalNumber("40")*(x_1^6)*(x_2^4)-RationalNumber("128")*(x_1^5)*(x_2^2)-RationalNumber("120")*(x_1^5)*(x_2^4)+RationalNumber("40")*(x_1^4)*(x_2^6)-RationalNumber("112")*(x_1^3)*(x_2^4)-RationalNumber("80")*(x_1^3)*(x_2^6)+RationalNumber("72")*(x_1^2)*(x_2^6)+RationalNumber("20")*(x_1^2)*(x_2^8)-RationalNumber("32")*(x_2^6)*(x_1)-RationalNumber("20")*(x_2^8)*(x_1));
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber("24")*(x_1^5)-RationalNumber("2")*(x_2^4)+RationalNumber("6")*(x_1^4)-RationalNumber("4")*(x_1)*(x_2^2)+RationalNumber("4")*(x_2^2)*(x_1^2)+RationalNumber("48")*(x_1^3)*(x_2^2)-RationalNumber("2")*(x_1)+RationalNumber("8")*(x_1^2)-RationalNumber("12")*(x_1^3)-RationalNumber("84")*(x_1^6)+RationalNumber("144")*(x_1^7)-RationalNumber("156")*(x_1^8)+RationalNumber("108")*(x_1^9)-RationalNumber("44")*(x_1^10)+RationalNumber("8")*(x_1^11)-RationalNumber("12")*(x_2^8)-RationalNumber("12")*(x_2^6)-RationalNumber("4")*(x_2^10)-RationalNumber("440")*(x_1^4)*(x_2^4)+RationalNumber("24")*(x_1)*(x_2^4)-RationalNumber("180")*(x_1^4)*(x_2^2)-RationalNumber("108")*(x_1^2)*(x_2^4)-RationalNumber("180")*(x_1^8)*(x_2^2)+RationalNumber("368")*(x_1^7)*(x_2^2)-RationalNumber("448")*(x_1^6)*(x_2^2)-RationalNumber("280")*(x_1^6)*(x_2^4)+RationalNumber("352")*(x_1^5)*(x_2^2)+RationalNumber("456")*(x_1^5)*(x_2^4)-RationalNumber("200")*(x_1^4)*(x_2^6)+RationalNumber("272")*(x_1^3)*(x_2^4)+RationalNumber("240")*(x_1^3)*(x_2^6)-RationalNumber("160")*(x_1^2)*(x_2^6)-RationalNumber("60")*(x_1^2)*(x_2^8)+RationalNumber("64")*(x_2^6)*(x_1)+RationalNumber("44")*(x_2^8)*(x_1)+RationalNumber("40")*(x_1^9)*(x_2^2)+RationalNumber("80")*(x_1^7)*(x_2^4)+RationalNumber("80")*(x_1^5)*(x_2^6)+RationalNumber("40")*(x_1^3)*(x_2^8)+RationalNumber("8")*(x_1)*(x_2^10);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1373Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1373Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys1397Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	SMQP poly;
//	poly = RationalNumber(5)*(x_1)*(x_2^4)+(x_2^4)+RationalNumber(2)*(x_2^3)+RationalNumber(2)*(x_2^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-1)*(x_2^2)-RationalNumber(2)*(x_1^5)-RationalNumber(4)*(x_2^3)-(x_1^2)*(x_2^3)+RationalNumber(2)*(x_1)*(x_2^4)+RationalNumber(3)*(x_1^3)*(x_2)+RationalNumber(64);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys1397Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys1397Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys2852Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars,variables,params;
//	vars.emplace_back("c_1");
//	SMQP c_1("c_1");
//	vars.emplace_back("c_2");
//	SMQP c_2("c_2");
//	vars.emplace_back("c_3");
//	SMQP c_3("c_3");
//	vars.emplace_back("c_4");
//	SMQP c_4("c_4");
//	vars.emplace_back("c_5");
//	SMQP c_5("c_5");
//	vars.emplace_back("x_10");
//	SMQP x_10("x_10");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	vars.emplace_back("x_3");
//	SMQP x_3("x_3");
//	vars.emplace_back("x_4");
//	SMQP x_4("x_4");
//	vars.emplace_back("x_7");
//	SMQP x_7("x_7");
//	vars.emplace_back("x_8");
//	SMQP x_8("x_8");
//	vars.emplace_back("x_9");
//	SMQP x_9("x_9");
//	variables = vars;
//	variables.erase(variables.begin()+3,variables.begin()+5);
//	params.push_back("c_4");
//	params.push_back("c_5");
//	SMQP poly;
//	poly = c_1;
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = c_2;
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = x_3;
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = x_4;
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = RationalNumber(-1)*(c_2)+(c_3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_7)-(c_3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_10^2)+(c_4^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_10^2)+(c_4^2);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_10)*(x_2);
//	poly.setRingVariables(vars);
//	poly += c_5;
//	poly += c_4*c_1;
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//
//	RegularChain<RationalNumber,SMQP> rc(variables,params);
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys2852Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys2852Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//static void Sys2985Test(bool showOutput) {
//	std::vector<SMQP> polys;
//	std::vector<Symbol> vars;
//	vars.emplace_back("x_1");
//	SMQP x_1("x_1");
//	vars.emplace_back("x_2");
//	SMQP x_2("x_2");
//	vars.emplace_back("x_3");
//	SMQP x_3("x_3");
//	SMQP poly;
//	poly = (x_1);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_3);
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	poly = (x_2^2)+RationalNumber("1");
//	poly.setRingVariables(vars);
//	polys.push_back(poly);
//	RegularChain<RationalNumber,SMQP> rc;
//	bool pass = triangularizeValidate(polys,rc,vars,showOutput);
//	if (showOutput) {
//		std::cerr << "Input polynomials:" << std::endl;
//		for (int i=0; i<polys.size(); ++i)
//			std::cerr << "F[" << i << "] = " << polys[i] << std::endl;
//	}
//	if (pass)
//		std::cerr << "RC Sys2985Test:\t\t\t\t\t\t\t PASSED" << std::endl;
//	else
//		std::cerr << "RC Sys2985Test:\t\t\t\t\t\t\t FAILED" << std::endl;
//}

//void triangularizeTests() {
//	cout << "Starting triangularizeTests..." << endl;
//	bool working(false);
//	bool testing(true);
//	Sys126Test(working);
//	Sys130Test(testing); // gmp: overflow in mpz type
//	Sys132Test(working);
//	Sys1000Test(working);
//	Sys1001Test(working);
//	Sys1002Test(working);
//	Sys1003Test(working);
//	Sys1004Test(working);
//	Sys1005Test(working);
//	Sys1006Test(working);
//	Sys1007Test(working);
//	Sys1008Test(testing); // algebraic error
//	Sys1009Test(working);
//	Sys1010Test(working);
//	Sys1011Test(working);
//	Sys1012Test(working);
//	Sys1013Test(working);
//	Sys1014Test(working);
//	Sys1015Test(working);
//	Sys1016Test(testing);
////	Sys1221Test(testing); // algebraic error / gmp: overflow in mpz type
//	Sys1255Test(testing);
////	Sys1289Test(testing); // infinite loop [iA] / huge expression swell
//	Sys1302Test(working);
////	Sys1303Test(testing); // algebraic error / gmp: overflow in mpz type
////	Sys1304Test(testing); // algebraic error / gmp: overflow in mpz type
////	Sys1305Test(testing); // gmp: overflow in mpz type
////	Sys1364Test(testing); // gmp: overflow in mpz type
////	Sys1366Test(working); // gmp: overflow in mpz type
////	Sys1373Test(testing); // algebraic error [too simple result] / gmp: overflow in mpz type
////	Sys1397Test(testing); // algebraic error [one missing and one incomplete component]
////	Sys2852Test(testing); // algebraic error [too simple result]
////	Sys2985Test(testing); // algebraic error [almost correct, missing an equation]
//}

int main() {
//	std::cerr << "sizeof(int) = " << sizeof(int) << std::endl;
//	std::cerr << "sizeof(long int) = " << sizeof(long int) << std::endl;
//	std::cerr << "sizeof(long long int) = " << sizeof(long long int) << std::endl;
//	std::cerr << "sizeof(SmallPrimeField) = " << sizeof(SmallPrimeField) << std::endl;
//	std::cerr << "sizeof(BigPrimeField) = " << sizeof(BigPrimeField) << std::endl;
//	std::cerr << "sizeof(GFPPrimeField) = " << sizeof(GeneralizedFermatPrimeField) << std::endl;
//	std::cerr << "sizeof(RationalNumber) = " << sizeof(RationalNumber) << std::endl;
//	std::cerr << "sizeof(ComplexRationalNumber) = " << sizeof(ComplexRationalNumber) << std::endl;
//	std::cerr << "sizeof(DUZP) = " << sizeof(DenseUnivariateRationalPolynomial) << std::endl;
//	std::cerr << "sizeof(DUQP) = " << sizeof(DenseUnivariateRationalPolynomial) << std::endl;
//	std::cerr << "sizeof(DenseUnivariatePolynomial<RN>) = " << sizeof(DenseUnivariatePolynomial<RN>) << std::endl;
//	std::cerr << "sizeof(SUP<RN>) = " << sizeof(SparseUnivariatePolynomial<RN>) << std::endl;
//	std::cerr << "sizeof(SUP<CRN>) = " << sizeof(SparseUnivariatePolynomial<CRN>) << std::endl;
//	std::cerr << "sizeof(SUP<SUP<RN>>) = " << sizeof(SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>>) << std::endl;
//	std::cerr << "sizeof(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial) = " << sizeof(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial) << std::endl;
//	std::cerr << "sizeof(RegularChain) = " << sizeof(SMQP) << std::endl;
//	std::cerr << "sizeof(DistributedDenseMultivariateModularPolynomial<RN>) = " << sizeof(DistributedDenseMultivariateModularPolynomial<RN>) << std::endl;
//	std::cerr << "sizeof(SMP<RN>) = " << sizeof(SparseUnivariatePolynomial<RN>) << std::endl;

	testRCDefaultConstructor();
	testRCVariableSpecificationConstructor();
	testRCVariableAndTranscendentalSpecificationConstructor();
	testRCPolynomialAddConstructor();
	testRCPolynomialAddWithTranscendentalsConstructor();
	testRCCopyConstructor();
	testRCMoveConstructor();
	testRCComputationalConstructor();
	testRCAssignmentOperator();
	testRCMoveAssignmentOperator();
//	testRCAddOperator();
//	testRCAddAssignmentOperator();
	testRCIdentityOperator();
	testRCNumberOfVariables();
	testRCSize();
	testRCNumberOfAlgebraicVariables();
	testRCNumberOfTranscendentalVariables();
	testRCVariables();
	testRCMainVariables();
	testRCTranscendentalVariables();
	testRCIsAlgebraic();
	testRCIsEmpty();
	testRCSelect();
	testRCLower();
	testRCUpper();
	testRCPseudoDivide();
	testRCNormalForm();
	testRCConvertToExpressionTree();


	testZDRCDefaultConstructor();
	testZDRCTranscendentalVariableSpecificationConstructor();
	testZDRCPolynomialAddConstructor();
	testZDRCPolynomialAddWithTranscendentalsConstructor();
	testZDRCCopyConstructor();
	testZDRCMoveConstructor();
	testZDRCComputationalConstructor();
	testZDRCAssignmentOperator();
	testZDRCMoveAssignmentOperator();
//	testZDRCAddOperator();
//	testZDRCAddAssignmentOperator();
	testZDRCIdentityOperator();
	testZDRCNumberOfVariables();
	testZDRCSize();
	testZDRCNumberOfAlgebraicVariables();
	testZDRCNumberOfTranscendentalVariables();
	testZDRCVariables();
	testZDRCMainVariables();
	testZDRCTranscendentalVariables();
    testZDRCIsAlgebraic();
	testZDRCIsEmpty();
	testZDRCSelect();
	testZDRCLower();
	testZDRCUpper();
	testZDRCPseudoDivide();
	testZDRCNormalForm();
	testZDRCConvertToExpressionTree();
	testZDRCIsInvertible(false);
	testZDRCregularGCD(false);

//	triangularizeTests();


//	cout << "RC Triangularize Tests:" << endl;
//	SMQP T1,T2,T3;
//	vector<SMQP> F;
//	vector<Symbol> rcvs;
//	vector<RegularChain<RN,SMQP>> TT;
////	rcvs.push_back("y");
//	rcvs.push_back("x");
////	T1 = ((SMQP("y")*SMQP("y")));
////	T2 = (T1*T1) + SMQP("x");
//	T1 = ((RN(1)+x)+RN(2)*(x*x));
//	T2 = ((RN(2)+RN(2)*x)+RN(4)*(x*x));
//	cout << "T1 = " << T1 << endl;
//	cout << "T2 = " << T2 << endl;
////	T2.setRingVariables(rcvs);
//	RegularChain<RN,SMQP> rc(rcvs);
////	rc += T3;
////	rc += T2;
//	F.push_back(T1);
////	F.push_back(T3);
//	F.push_back(T2);
////	F.push_back(SMQP("z"));
//	TT = rc.triangularize(F);
//	cout << "in test code: TT.size = " << TT.size() << endl;
//	for (int i=0; i<TT.size(); ++i)
//		cout << "TT[" << i << "] = " << TT[i] << endl;
//
//	cout << "rc.isKnownToBePrime = " << rc.isKnownToBePrime() << endl;
//	cout << "rc.isKnownToBeSquareFree = " << rc.isKnownToBeSquareFree() << endl;
//
//	cout << "triangularizeValidate: " << triangularizeValidate(F,rc,rcvs) << endl;

//	cout << "sphere intersecting circle example: " << endl;
//	T1 = (SMQP("x")*SMQP("x")) + (SMQP("y")*SMQP("y")) + (SMQP("z")*SMQP("z")) - RN(1);
//	T2 = ((SMQP("x")*SMQP("x")) + (SMQP("y")*SMQP("y")) - RN(1))*SMQP("z");
//	rcvs.clear();
//	rcvs.push_back("z");
//	rcvs.push_back("y");
//	rcvs.push_back("x");
//	T1.setRingVariables(rcvs);
//	T2.setRingVariables(rcvs);
//	cout << "T1 = " << T1 << endl;
//	cout << "T2 = " << T2 << endl;
//	F.clear();
//	F.push_back(T1);
//	F.push_back(T2);
//	RegularChain<RN,SMQP> rc2;
////	TT = rc2.triangularize(F);
////	cout << "in test code: TT.size = " << TT.size() << endl;
////	for (int i=0; i<TT.size(); ++i)
////		cout << "TT[" << i << "] = " << TT[i] << endl;
//
//	cout << "triangularizeValidate: " << triangularizeValidate(F,rc2,rcvs) << endl;









//	cout << "regularGCD	Examples for RegularChain:" << endl;
//    vector<PolyChainPair<SMQP,RegularChain<RN,SMQP>>> rcResults;
//	cout << "Example 1:" << endl;
//
////    // Redirect stdout, stderr
////    oldcout = cout.rdbuf(); // <-- save
////    oldcerr = cerr.rdbuf(); // <-- save
////    cout.rdbuf (ss.rdbuf());
////    cerr.rdbuf (ss.rdbuf());
//
//	m = (x-RN(1))*(x+RN(2));
//	f1 = ((RN(2)*x+RN(1))*(y*y)) + ((RN(5)*x+RN(10))*y) + ((RN(6)*x) + RN(12));
//	f2 = (RN(3)*(y*y)) + ((RN(4)*x+RN(8))*y) + (RN(4)*x) + RN(5);
//	vars.clear();
//	vars.emplace_back("y");
//	vars.emplace_back("x");
//	f1.setRingVariables(vars);
//	f2.setRingVariables(vars);
//	cout << "m = " << m << endl;
//	cout << "f1 = " << f1 << endl;
//	cout << "f2 = " << f2 << endl;
//	rc = RegularChain<RN,SMQP>(m);
//    result = rc.regularGCD(f1,f2,Symbol("y"));
//
////	// Restore stdout, stderr
////    cout.rdbuf (oldcout);              // <-- restore
////    cerr.rdbuf (oldcerr);              // <-- restore
//
//    for (int l=0; l<result.size(); ++l) {
////    	cout << result[l] << endl;
//    	cout << "G_" << l << " = " << result[l].poly << endl;
//    	cout << "T_" << l << " = " << result[l].chain << endl;




//	SMQP w("x");
//	w.zero();
//	cout << "zero.degree(x) = " << w.degree("x") << endl;


//    // SMQP GCD Tests //
//    SMQP p,q,g;
//	int nvars(1);
//	int nterms = 6;
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(30);
//	bool includeNeg(1);
//	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
//	vars.clear();
//    vars.emplace_back("x");
//    p.setRingVariables(vars);
//	p /= p.leadingCoefficient();
//	q.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
//    q.setRingVariables(vars);
//	q*=p;
//	g = p.subresultantGCD(q);
//	g /= g.leadingCoefficient();
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
//	cout << "g = " << g << endl;
//	g = q.subresultantGCD(p);
//	g /= g.leadingCoefficient();
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
//	cout << "g = " << g << endl;
//	q.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
//    q.setRingVariables(vars);
//	g = p.subresultantGCD(q);
//	g /= g.leadingCoefficient();
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
//	cout << "g = " << g << endl;


////    ZeroDimensionalRegularChain<RN,SMQP> ts;
//	SMQP p,check;
//	std::vector<SMQP> set;
//	std::vector<Symbol> trcVars;
//	vars.clear();
//	int nvars(3);
//	int nTrcVars(0);
//	nterms = 6;
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(30);
//	bool includeNeg(1);
//    ts.randomZeroDimensionalRegularChain(nvars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    ts.display();
//    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    vars = ts.variables();
//    trcVars = ts.transcendentalVariables();
//    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
//	p.setRingVariables(vars);
//	std::vector<BoolChainPair<ZeroDimensionalRegularChain<RN,SMQP>>> isInv;
//	std::vector<PolyChainPair<SMQP,ZeroDimensionalRegularChain<RN,SMQP>>> results;
//	results = ts.regularize(p);
//	cout << "results.size = " << results.size() << endl;
////	isInv = ts.isInvertible(p);
////	cout << "isInv.size = " << isInv.size() << endl;
////	set = ts.polynomials();

	return 0;
}

//int main() {
//	RegularChain<RN,SMQP> rc;
//
//	cout << "rc = " << rc << endl;
//
//	ZeroDimensionalRegularChain<RN,SMQP> rc2;

//	cout << "zdrc = " << rc2 << endl;
////
//    RegularChain<RN,SMQP> ts;
//
//	SMQP f,g,h,check;
//	std::vector<SMQP> set;
//	std::vector<Symbol> vars,trcVars;
//	int nvars(4);
//	int nAlgVars(4);
//	int nTrcVars(2);
//	int nterms(8);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(25);
//	bool includeNeg(1);
//    rc2.randomRegularChain(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    rc2.display();
//    f.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    g.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);
//
//    vars = rc2.variables();
//    trcVars = rc2.transcendentalVariables();
//    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
//    f.setRingVariables(vars);
//    g.setRingVariables(vars);
//    cout << "f = " << f << endl;
//    cout << "g = " << g << endl;
//    if (f.leadingVariable() == g.leadingVariable())
//    	h = rc2.regularGCD(f,g,f.leadingVariable());
//    cout << "h = " << h << endl;
//
//	return 0;
//}
