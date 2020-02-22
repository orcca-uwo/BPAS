#include <bpas.h>
#include <vector>
#include "../MapleTestTool/MapleTestTool.hpp"

using namespace std;

void testDefaultConstructor() {
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;

    if (ts.size() != 0 || ts.numberOfVariables() != 0 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 0) {
        std::cerr << "TriangularSet default constructor test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet default constructor test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testVariableSpecificationConstructor() {
	vector<Symbol> xs;
	xs = Symbol::randomElements(2);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(xs);

    if (ts.size() != 0 || ts.numberOfVariables() != 2 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 0) {
        std::cerr << "TriangularSet variable specification constructor test:\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts.size = " << ts.size() << ", should be 2" << std::endl;
        std::cerr << "ts.numberOfVariables = " << ts.numberOfVariables() << ", should be 2" << std::endl;
        std::cerr << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << ", should be 0" << std::endl;
        std::cerr << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << ", should be 0" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet variable specification constructor test:\t\t\t\t PASSED" << std::endl;
}

void testVariableAndTranscendentalSpecificationConstructor() {
	vector<Symbol> vars;
	vars = Symbol::randomElements(5);
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+5);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(xs,ps);

    if (ts.size() != 0 || ts.numberOfVariables() != 2 || ts.numberOfAlgebraicVariables() != 0 || ts.numberOfTranscendentalVariables() != 3) {
        std::cerr << "TriangularSet variable and transcendental specification constructor test:\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet variable and transcendental specification constructor test:\t PASSED" << std::endl;
}

void testPolynomialAddConstructor() {
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
	SparseMultivariateRationalPolynomial p(nvars),ex(s);
	vector<Symbol> x;
	x.push_back(s);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(xs);
	p += ex;
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p);

    if (ts.size() != 1 || ts.numberOfVariables() != 3 || ts.numberOfAlgebraicVariables() != 1 || ts.numberOfTranscendentalVariables() != 0 || xs != ts.variables() || x != ts.mainVariables()) {
        std::cerr << "TriangularSet polynomial addition constructor test:\t\t\t\t FAILED" << std::endl;
        std::cerr << "polynomial being added: " << p << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet polynomial addition constructor test:\t\t\t\t PASSED" << std::endl;
}

void testPolynomialAddWithTranscendentalsConstructor() {
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
	SparseMultivariateRationalPolynomial p(nvars),ex(s);
	vector<Symbol> x;
	x.push_back(s);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
	p += ex;
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps);
    vector<Symbol> vars;
    vars = orderPreservingSetIntersection(xs,p.variables());

    if (ts.size() != 1 || ts.numberOfVariables() != vars.size() || ts.numberOfAlgebraicVariables() != 1 || ts.numberOfTranscendentalVariables() != 2 || vars != ts.variables() || x != ts.mainVariables() || ps != ts.transcendentalVariables()) {
        std::cerr << "TriangularSet polynomial addition with transcendentals constructor test:\t FAILED" << std::endl;
        std::cerr << "ts.mainVariables: ";
        for (int i=0; i<ts.mainVariables().size(); ++i)
        	cerr << ts.mainVariables()[i] << " ";
        cerr << endl;
        cerr << "x = " << x[0] << endl;
        cerr << "ts.size = " << ts.size() << endl;
        cerr << "ts.numberOfVariables = " << ts.numberOfVariables() << endl;
        cerr << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << endl;
        cerr << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << endl;
        std::cerr << "ts.variables: ";
        printVariables(ts.variables());
        std::cerr << "xs: ";
        printVariables(xs);
        exit(1);
    }

    std::cerr << "TriangularSet polynomial addition with transcendentals constructor test:\t PASSED" << std::endl;
}

void testCopyConstructor() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
   //ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
//    std::cout << "ts = " << ts << std::endl;
//    cout << "ts.numberOfVariables = " << ts.numberOfVariables() << std::endl;
//    cout << "ts.numberOfAlgebraicVariables = " << ts.numberOfAlgebraicVariables() << std::endl;
//    cout << "ts.numberOfTranscendentalVariables = " << ts.numberOfTranscendentalVariables() << std::endl;
}

void testMoveConstructor() {
    
//    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
//	int nterms(14);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(4);
//	bool includeNeg(1);
//    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
//    std::vector<Symbol> vs(ts.mainVariables());
//    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts2(ts.lower(vs[0]));
//    ts = ts.lower(vs[0]);
//    if (ts.variables() != ts2.variables() || ts.transcendentalVariables() != ts2.transcendentalVariables() || ts.polynomials() != ts2.polynomials()) {
//        std::cerr << "TriangularSet copy constructor test:\t\t\t\t\t\t FAILED" << std::endl;
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

//    std::cerr << "TriangularSet move constructor test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testComputationalConstructor() {}

void testAssignmentOperator() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts,ts2;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    ts2 = ts;
    if (ts.variables() != ts2.variables() || ts.transcendentalVariables() != ts2.transcendentalVariables() || ts.polynomials() != ts2.polynomials()) {
        std::cerr << "TriangularSet assignment operator test:\t\t\t\t\t\t FAILED" << std::endl;
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

    std::cerr << "TriangularSet assignment operator test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testMoveAssignmentOperator() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts,ts2,ts3;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    ts2 = ts;
    ts3 = move(ts);
    if (ts3.variables() != ts2.variables() || ts3.transcendentalVariables() != ts2.transcendentalVariables() || ts3.polynomials() != ts2.polynomials() || !ts.isEmpty()) {
        std::cerr << "TriangularSet move assignment operator test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts.isEmpty = " << ts.isEmpty() << std::endl;
        exit(1);
    }
    std::cerr << "TriangularSet move assignment operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testAddOperator() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SparseMultivariateRationalPolynomial p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+4);
	p1.randomPolynomial(nvars-1,nterms,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	xs.push_back(vars[4]);
	p2.randomPolynomial(nvars-2,nterms,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[4]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);
	
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p1,ps);
    ts = ts + p2;
	std::vector<SparseMultivariateRationalPolynomial> varTSPolys,fixTSPolys;
    varTSPolys.push_back(p2);
    varTSPolys.push_back(p1);
    varTSPolys.push_back(SparseMultivariateRationalPolynomial());
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts2(xs,ps);
    ts2 = ts2 + p1;
    ts2 = ts2 + p2;
    fixTSPolys.push_back(p1);
    fixTSPolys.push_back(SparseMultivariateRationalPolynomial());
    fixTSPolys.push_back(p2);

    if (ts.polynomials() != varTSPolys || ts2.polynomials() != fixTSPolys) {
        std::cerr << "TriangularSet addition operator test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet addition operator test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testAddAssignmentOperator() {
	int nvars(5);
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
	vector<Symbol> vars;
	vars = Symbol::randomElements(nvars);
	SparseMultivariateRationalPolynomial p1,p2;
	vector<Symbol> xs;
	xs.insert(xs.end(),vars.begin(),vars.begin()+2);
	vector<Symbol> ps;
	ps.insert(ps.end(),vars.begin()+2,vars.begin()+4);
	p1.randomPolynomial(nvars-1,nterms,coefBound,sparsity,includeNeg);
	vector<Symbol> vs(xs);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p1.setRingVariables(vs);
	xs.push_back(vars[4]);
	p2.randomPolynomial(nvars-2,nterms,coefBound,sparsity,includeNeg);
	vs.clear();
	vs.push_back(vars[4]);
	vs.insert(vs.end(),ps.begin(),ps.end());
	p2.setRingVariables(vs);
	
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p1,ps);
    ts += p2;
	std::vector<SparseMultivariateRationalPolynomial> varTSPolys,fixTSPolys;
    varTSPolys.push_back(p2);
    varTSPolys.push_back(p1);
    varTSPolys.push_back(SparseMultivariateRationalPolynomial());
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts2(xs,ps);
    ts2 += p1;
    ts2 += p2;
    fixTSPolys.push_back(p1);
    fixTSPolys.push_back(SparseMultivariateRationalPolynomial());
    fixTSPolys.push_back(p2);

    if (ts.polynomials() != varTSPolys || ts2.polynomials() != fixTSPolys) {
        std::cerr << "TriangularSet addition assignment operator test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet addition assignment operator test:\t\t\t\t PASSED" << std::endl;
}

void testIdentityOperator() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts2;
    ts2 = ts;
    if (ts.polynomials() != ts2.polynomials()) {
        std::cerr << "TriangularSet identity operator test:\t\t\t\t\t\t FAILED" << std::endl;
        for (int i=0; i<ts.polynomials().size(); ++i)
        	std::cerr << "ts[" << i << "] = " << ts.polynomials()[i] << std::endl;
        for (int i=0; i<ts2.polynomials().size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2.polynomials()[i] << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet identity operator test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testNumberOfVariables() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.numberOfVariables() != 4) {
        std::cerr << "TriangularSet number of variables test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet number of variables test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testSize() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.size() != 3) {
        std::cerr << "TriangularSet size test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet size test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testNumberOfAlgebraicVariables() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(4);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    std::vector<SparseMultivariateRationalPolynomial> polys(ts.polynomials()),polys2;
    for (int i=0; i<polys.size(); ++i) {
    	if (!polys[i].isZero())
    		polys2.push_back(polys[i]);
    }
    if (ts.mainVariables().size() != polys2.size()) {
        std::cerr << "TriangularSet numberOfAlgebraicVariables test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet numberOfAlgebraicVariables test:\t\t\t\t\t PASSED" << std::endl;
}

void testNumberOfTranscendentalVariables() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
    ts.randomTriangularSet(4,3,2,nterms,coefBound,sparsity,includeNeg);
    if (ts.transcendentalVariables().size() != 2) {
        std::cerr << "TriangularSet numberOfTranscendentalVariables test:\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet numberOfTranscendentalVariables test:\t\t\t\t PASSED" << std::endl;
}

void testVariables() {
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
	SparseMultivariateRationalPolynomial p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps);
    vector<Symbol> vars;
    vars = orderPreservingSetIntersection(xs,p.variables());

    if (ts.variables() != vars) {
        std::cerr << "TriangularSet variables test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet variables test:\t\t\t\t\t\t\t PASSED" << std::endl;
}


void testMainVariables() {
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
	SparseMultivariateRationalPolynomial p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());

    if (ts.mainVariables() != vars) {
        std::cerr << "TriangularSet main variables test:\t\t\t\t\t\t FAILED" << std::endl;
        printVariables(ts.mainVariables());
        printVariables(vars);
        exit(1);
    }

    std::cerr << "TriangularSet main variables test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testTranscendentalVariables() {
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
	SparseMultivariateRationalPolynomial p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());

    if (ts.transcendentalVariables() != ps) {
        std::cerr << "TriangularSet transcendental variables test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet transcendental variables test:\t\t\t\t\t PASSED" << std::endl;
}

void testIsAlgebraic() {
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
	SparseMultivariateRationalPolynomial p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());
    vector<Symbol> algVar;
    algVar = setIntersection(ts.mainVariables(),vars);

    if (ts.isAlgebraic(p.leadingVariable()) != (algVar == vars)) {
        std::cerr << "TriangularSet isAlgebraic test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet isAlgebraic test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testIsEmpty() {
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
	SparseMultivariateRationalPolynomial p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps),ts2;

    if (ts.isEmpty() == true || ts2.isEmpty() == false) {
        std::cerr << "TriangularSet isEmpty test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet isEmpty test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testSelect() {

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
	SparseMultivariateRationalPolynomial p(nvars);

	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
	p.setRingVariables(vs);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts(p,ps);
    vector<Symbol> vars;
    vars.push_back(p.leadingVariable());
    SparseMultivariateRationalPolynomial zero;

    if (ts.select(p.leadingVariable()) != p || ts.select(s) != zero) {
        std::cerr << "TriangularSet select test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet select test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testLower() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts,ts2;
    std:vector<SparseMultivariateRationalPolynomial> tslp,ts2p;
    std::vector<Symbol> vars;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
    ts.randomTriangularSet(4,4,2,nterms,coefBound,sparsity,includeNeg);
    vars = ts.mainVariables();
    ts.lower(vars[0],ts2);
//	std::cerr << "ts2 in testLower = " << ts2 << std::endl;
//	ts.display();
//	ts2.display();
	tslp = ts.polynomials();
	tslp[0].zero();
	ts2p = ts2.polynomials();
    if (tslp != ts2p) {
        std::cerr << "TriangularSet lower test:\t\t\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts lower polys:" << std::endl;
        for (int i=0; i<tslp.size(); ++i)
        	std::cerr << "ts[" << i << "] = " << tslp[i] << std::endl;
        std::cerr << "ts2 polys:" << std::endl;
        for (int i=0; i<ts2p.size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2p[i] << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet lower test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testUpper() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts,ts2;
    std:vector<SparseMultivariateRationalPolynomial> tsup,ts2p;
    std::vector<Symbol> vars;
	int nterms(14);
	unsigned long int coefBound(6ul);
	degree_t sparsity(400);
	bool includeNeg(1);
    ts.randomTriangularSet(4,4,2,nterms,coefBound,sparsity,includeNeg);
    vars = ts.mainVariables();
    ts.upper(vars[3],ts2);
//    ts.display();
//    ts2.display();
	tsup = ts.polynomials();
	for (int i=tsup.size()-1; i>=3; --i)
		tsup[i].zero();
//	tsup.erase(tsup.end()-1,tsup.end());
	ts2p = ts2.polynomials();
    if (tsup != ts2p) {
        std::cerr << "TriangularSet upper test:\t\t\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "ts upper polys:" << std::endl;
        for (int i=0; i<tsup.size(); ++i)
        	std::cerr << "ts[" << i << "] = " << tsup[i] << std::endl;
        std::cerr << "ts2 polys:" << std::endl;
        for (int i=0; i<ts2p.size(); ++i)
        	std::cerr << "ts2[" << i << "] = " << ts2p[i] << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet upper test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testPseudoDivide() {
    
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
	SparseMultivariateRationalPolynomial p,check;
	std::vector<SparseMultivariateRationalPolynomial> set;
	std::vector<Symbol> vars,trcVars;
	int nvars(3);
	int nAlgVars(3);
	int nTrcVars(0);
	int nterms(3);
	unsigned long int coefBound(2ul);
	degree_t sparsity(50);
	bool includeNeg(1);
    ts.randomTriangularSet(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
    // ts.display();
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
	SparseMultivariateRationalPolynomial r,c;
	vector<SparseMultivariateRationalPolynomial> Q;
	r = ts.pseudoDivide(p,&Q,&c);
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
//        std::cerr << "TriangularSet pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "TriangularSet pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;

    ExpressionTree pTree = p.convertToExpressionTree();
//    ExpressionTree tsTree = ts.convertToExpressionTree();
    ExpressionTree tsTree;
    tsTree.fromVector<SparseMultivariateRationalPolynomial>(set);
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
	std::string procStr = "multiDivisorPseudoDivision := proc (p, ts, names, q, m) local seqLoop, tsSize, r, tmpQ, Q, tmpM, totalM, i, j; seqLoop := [seq(i, i = nops(ts) .. 1, -1)]; tsSize := nops(ts); r := p; Q := [seq(0, i = 1 .. tsSize)]; totalM := 1; for i in seqLoop do tmpQ := 0; tmpM := 1; r := prem(r, ts[i], names[i], 'tmpM', 'tmpQ'); totalM := totalM*tmpM; for j to tsSize do Q[j] := simplify(tmpM*Q[j]) end do; Q[i] := Q[i]+tmpQ end do; assign(q = Q); assign(m = totalM); return r end proc:";
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
        std::cerr << "TriangularSet pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;

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

    std::cerr << "TriangularSet pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;
}

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
    for (int i = 0; i < 2; ++i) {
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
        std::cerr << "TriangularSet normalForm test:\t\t\t\t\t\t\t FAILED" << std::endl;

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

    std::cerr << "TriangularSet normalForm test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testConvertToExpressionTree() {

    stringstream ss;
    std::string str,str2;
    int nvars(4);
    int nAlgVars(4);
    int nTrcVars(2);
    int nterms(14);
    unsigned long int coefBound(6ul);
    degree_t sparsity(400);
    bool includeNeg(1);
    TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
    ts.randomStronglyNormalizedTriangularSet(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
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

    ts = TriangularSet<RN,SparseMultivariateRationalPolynomial>();
    ExpressionTree et2;
    et2 = ts.convertToExpressionTree();
//  cout << "et2 = " << et2.toString() << endl;
    ss = stringstream();
    ss << ts;
    str2 = ss.str();
//  cout << "ts = " << str2 << endl;
    if (!cmpBool || et2.toString() != str2) {
        std::cerr << "TriangularSet convertToExpressionTree test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "TriangularSet convertToExpressionTree test:\t\t\t\t\t PASSED" << std::endl;
}

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
//	std::cerr << "sizeof(TriangularSet) = " << sizeof(SparseMultivariateRationalPolynomial) << std::endl;
//	std::cerr << "sizeof(DistributedDenseMultivariateModularPolynomial<RN>) = " << sizeof(DistributedDenseMultivariateModularPolynomial<RN>) << std::endl;
//	std::cerr << "sizeof(SMP<RN>) = " << sizeof(SparseUnivariatePolynomial<RN>) << std::endl;

	testDefaultConstructor();
	testVariableSpecificationConstructor();
	testVariableAndTranscendentalSpecificationConstructor();
	testPolynomialAddConstructor();
	testPolynomialAddWithTranscendentalsConstructor();
	testCopyConstructor();
	testMoveConstructor();
	testComputationalConstructor();
	testAssignmentOperator();
	testMoveAssignmentOperator();
	testAddOperator();
	testAddAssignmentOperator();
	testIdentityOperator();
	testNumberOfVariables();
	testSize();
	testNumberOfAlgebraicVariables();
	testNumberOfTranscendentalVariables();
	testVariables();
	testMainVariables();
	testTranscendentalVariables();
	testIsAlgebraic();
	testIsEmpty();
	testSelect();
	testLower();
	testUpper();
	testPseudoDivide();
	testNormalForm();
	testConvertToExpressionTree();
	
//	TriangularSet<RN,SparseMultivariateRationalPolynomial> ts;
//	int i;
//	int numvar;
//	int nterms = 7;
//	unsigned long int coefBound = 5ul;
//	degree_t sparsity = 4;
//	bool includeNeg = 1;
//	SparseMultivariateRationalPolynomial p1(3),p2(2),p3(1),ex("x"),wy("y");

//	p1.randomPolynomial(3,nterms,coefBound,sparsity,includeNeg);
//	p2.randomPolynomial(2,nterms,coefBound,sparsity,includeNeg);
//	p3.randomPolynomial(1,nterms,coefBound,sparsity,includeNeg);

//	vector<Symbol> names;

//	Symbol s = "z";
//	names.insert(names.begin(),s);
//	p3.setRingVariables(names);

//	s = "y";
//	names.insert(names.begin(),s);
//	p2.setRingVariables(names);
//	p2 *= wy;
//	
//	s = "x";
//	names.insert(names.begin(),s);
//	p1.setRingVariables(names);
//	p1 += ex;
//	
//	
//	// Checking leading variable //
//	TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial> t(names),t2;
//	t += p1;
//	cout << "t = " << t << endl;
//	t += p2;
//	cout << "t = " << t << endl;
//	t += p3;
//	cout << "t = " << t << endl;
//	
//	// Checking upper and lower //
//	t.lower("y",t2);
//	cout << "t2.variables:" << endl;
//	names = t2.variables();
//	for (i=0;i<names.size();++i) {
//		cout << "names[" << i << "] = " << names.at(i) << endl;
//	}
//	cout << "chain under y = " << t2 << endl;
//	t.upper("y",t2);
//	cout << "chain over y = " << t2 << endl;

//	// Checking pseudoDivide //
//	SparseMultivariateRationalPolynomial p4(3);
//	p4.randomPolynomial(3,nterms,coefBound,sparsity,includeNeg);
//	names.clear();
//	names.push_back("x");
//	names.push_back("y");
//	names.push_back("z");
//	p4.setRingVariables(names);
//	cout << "mvar(p4) = " << p4.leadingVariable() << ", p4 = " << p4 << endl;
//	SparseMultivariateRationalPolynomial p5(1);
//	p5.randomPolynomial(2,nterms,coefBound,sparsity,includeNeg);
//	names.clear();
//	names.push_back("x");
//	names.push_back("z");
//	p5.setRingVariables(names);
//	cout << "mvar(p5) = " << p5.leadingVariable() << ", p5 = " << p5 << endl;
////	p4 = p5;
////	names.clear();
////	names.push_back("x");
////	// One argument version //
////	SparseMultivariateRationalPolynomial rem;
////	rem = t.pseudoDivide(p4);
////	cout << "rem = " << rem << endl;
////	// Two argument version //
////	vector<SparseMultivariateRationalPolynomial> quo;
////	rem = t.pseudoDivide(p4,&quo);
////	cout << "rem = " << rem << endl;
////	for (i=0;i<quo.size();++i)
////		cout << "quo[" << i << "] = " << quo[i] << endl;
////	cout << "quo.size() = " << quo.size() << endl;

	// Three argument version //
// 	SparseMultivariateRationalPolynomial c,check;
// 	check.zero();
// 	rem = t.pseudoDivide(p4,&quo,&c);
// 	cout << "rem = " << rem << endl;
// 	cout << "c = " << c << endl;
// 	for (i=0;i<quo.size();++i)
// 		cout << "quo[" << i << "] = " << quo[i] << endl;
// 	names = t.variables();
// 	for (i=0;i<3;i++) {
// //		cout << "quo[" << i << "] vars:" << endl;
// //		printVariables(quo[i].variables());
// //		cout << "t.select(names[" << i << "]) vars:" << endl;
// //		printVariables(t.select(names[i]).variables());
// //		cout << "check vars: " << endl;
// //		printVariables(check.variables());
// 		check += quo[i]*t.select(names[i]);
// //		cout << "here?" << endl;
// 	}
	// cerr << "compute constraint:" << endl;
	// check += rem;
	// check -= (c*p4);
	// cout << "check = " << check << endl;

//	cout << "Test normalForm:" << endl;
//	SparseMultivariateRationalPolynomial one,temp;
//	one.one();
//	names.clear();
//	names.push_back("x");
//	ex.setRingVariables(names);
//	ex = ex^(p1.leadingVariableDegree().get_si());
//	temp = p1.initial();
//	temp -= one;
//	temp *= ex;
//	p1 -= temp;
//	cout << "p1 = " << p1 << endl;
////	cout << "p1.initial() = " << p1.initial() << endl;
//	wy = wy^(p2.leadingVariableDegree().get_si());
////	cout << "wy = " << wy << endl;
//	temp = p2.initial();
//	temp -= one;
//	temp *= wy;
////	cout << "temp = " << temp << endl;
//	p2 -= temp;
//	cout << "p2 = " << p2 << endl;
////	cout << "p2.initial() = " << p2.initial() << endl;
//	p3 /= p3.leadingCoefficient();
//	cout << "p3 = " << p3 << endl;
////	cout << "p3.initial() = " << p3.initial() << endl;
//	SparseMultivariateRationalPolynomial q(3);
//	q = move(p5);
//	q /= q.leadingCoefficient();
//	cout << "q = " << q << endl;
//	cout << "p5 = " << p5 << endl;
//	ts = TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>(p1);
//	cout << "ts = " << ts << endl;
//	ts += p2;
//	cout << "ts = " << ts << endl;
//	ts += p3;
//	cout << "ts = " << ts << endl;
//	SparseMultivariateRationalPolynomial r;
//	vector<SparseMultivariateRationalPolynomial> Q;
//	r = ts.normalForm(q,&Q);
//	cout << "r = " << r << endl;
//	for (int i=0; i<Q.size(); ++i)
//		cout << "Q[" << i << "] = " << Q[i] << endl;
	
//	cout << "Test variable order management:" << endl;
//	SparseMultivariateRationalPolynomial p6(5);
//	p6.randomPolynomial(5,nterms,coefBound,100,includeNeg);
//	names.clear();
//	names.push_back("x");
//	names.push_back("y");
//	names.push_back("z");
//	names.push_back("t");
//	names.push_back("u");
//	p6.setRingVariables(names);
//	SparseMultivariateRationalPolynomial p7(5);
//	p7.randomPolynomial(5,nterms,coefBound,100,includeNeg);
//	names.clear();
//	names.push_back("y");
//	names.push_back("x");
//	names.push_back("u");
//	names.push_back("t");
//	names.push_back("w");
//	p7.setRingVariables(names);
//	cout << "p6 = " << p6 << endl;
//	cout << "p7 = " << p7 << endl;
//	vector<Symbol> names2;
//	names.clear();
//	names.push_back("x");
//	names.push_back("y");
//	names.push_back("z");
//	names2.clear();
//	names2.push_back("u");
//	names2.push_back("t");
//	names2.push_back("w");
//	ts = TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>(names,names2);
//	ts += p6;
//	ts += p7;
//	ts.display();
//	cout << "variables: ";
//	for (auto i=0; i<ts.variables().size(); ++i)
//		cout << ts.variables()[i];
//	cout << endl;
//	cout << "transcendental variables: ";
//	for (auto i=0; i<ts.transcendentalVariables().size(); ++i)
//		cout << ts.transcendentalVariables()[i];
//	cout << endl;

	

	return 0;
}
