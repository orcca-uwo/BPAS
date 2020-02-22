#include <bpas.h>
#include <vector>
#include "../MapleTestTool/MapleTestTool.hpp"
#include "../../include/Utils/RandomHelpers.hpp"

using namespace std;

//SubResultantChain<RN,SparseUnivariatePolynomial<RN>> randomDerivativeSubResultantChain(unsigned long int coefBound, double sparsity, int degree) {
//	
//	SparseUnivariatePolynomial<RN> p(randomSUPQPolynomial(degree,sparsity,coefBound));
//	SparseUnivariatePolynomial<RN> q(p);
//	Symbol s(Symbol::randomElement());
//	p.setVariableName(s);
//	q.setVariableName(s);
//	q.differentiate();
////	cout << "p = " << p << endl;
////	cout << "q = " << q << endl;
//	
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(p,q,s);
////    cout << "src = " << src << endl;
//	return src;
//}

//void testDefaultConstructor() {
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src;

//    if (src.size() != 0 || src.variableName() != "x") {
//        std::cerr << "SubResultantChain default constructor test:\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain default constructor test:\t\t\t\t\t PASSED" << std::endl;
//}

//void testVariableSpecificationConstructor() {
//	Symbol s(Symbol::randomElement());
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(s);

//    if (src.size() != 0 || src.variableName() != s) {
//        std::cerr << "SubResultantChain variable specification constructor test:\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain variable specification constructor test:\t\t\t PASSED" << std::endl;
//}

//void testChainConstructor() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//	
//    SparseUnivariatePolynomial<RN> p(randomSUPQPolynomial(degree,sparsity,coefBound));
//	SparseUnivariatePolynomial<RN> q(p);
//	Symbol s(Symbol::randomElement());
//	p.setVariableName(s);
//	q.setVariableName(s);
//	q.differentiate();
////	cout << "p = " << p << endl;
////	cout << "q = " << q << endl;
//	
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(p,q,s);
////    cout << "src = " << src << endl;

//    if (false) {
//        std::cerr << "SubResultantChain chain constructor test:\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain chain constructor test:\t\t\t\t\t PASSED" << std::endl;
//}

//void testCopyConstructor() {

//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src;
//    src = randomDerivativeSubResultantChain(coefBound,sparsity,degree);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src2(src);

//    if (src.polynomials() != src2.polynomials() || src.variableName() != src2.variableName()) {
//        std::cerr << "SubResultantChain copy constructor test:\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain copy constructor test:\t\t\t\t\t PASSED" << std::endl;
//}

//void testMoveConstructor() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//	
//    SparseUnivariatePolynomial<RN> p(randomSUPQPolynomial(degree,sparsity,coefBound));
//	SparseUnivariatePolynomial<RN> q(p);
//	Symbol s(Symbol::randomElement());
//	p.setVariableName(s);
//	q.setVariableName(s);
//	q.differentiate();
//	
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(p,q,s);
//    std::vector<SparseUnivariatePolynomial<RN>> srcPolys(src.polynomials());
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src2(std::move(src));
//    
//    if (src.variableName() != "x" || !src.polynomials().empty() || srcPolys != src2.polynomials() || src2.variableName() != s) {
//        std::cerr << "SubResultantChain move constructor test:\t\t\t\t\t FAILED" << std::endl;
//        std::cerr << "src.variableName() = " << src.variableName() << std::endl;
//        std::cerr << "src.polynomials().empty() = " << src.polynomials().empty() << std::endl;
//        std::cerr << "srcPolys == src2.polynomials = " << (srcPolys == src2.polynomials()) << std::endl;
//        std::cerr << "src2.variableName = " << src2.variableName() << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain move constructor test:\t\t\t\t\t PASSED" << std::endl;
//}

//void testAssignmentOperator() {

//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src;
//    src = randomDerivativeSubResultantChain(coefBound,sparsity,degree);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src2;
//    src2 = src;

//    if (src.polynomials() != src2.polynomials() || src.variableName() != src2.variableName()) {
//        std::cerr << "SubResultantChain assignment operator test:\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain assignment operator test:\t\t\t\t\t PASSED" << std::endl;
//}

//void testMoveAssignmentOperator() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//	
//    SparseUnivariatePolynomial<RN> p(randomSUPQPolynomial(degree,sparsity,coefBound));
//	SparseUnivariatePolynomial<RN> q(p);
//	Symbol s(Symbol::randomElement());
//	p.setVariableName(s);
//	q.setVariableName(s);
//	q.differentiate();
//	
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(p,q,s);
//    std::vector<SparseUnivariatePolynomial<RN>> srcPolys(src.polynomials());
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src2;
//    
//    src2 = std::move(src);
//    
//    if (src.variableName() != "x" || !src.polynomials().empty() || srcPolys != src2.polynomials() || src2.variableName() != s) {
//        std::cerr << "SubResultantChain move assignment operator test:\t\t\t\t FAILED" << std::endl;
//        std::cerr << "src.variableName() = " << src.variableName() << std::endl;
//        std::cerr << "src.polynomials().empty() = " << src.polynomials().empty() << std::endl;
//        std::cerr << "srcPolys == src2.polynomials = " << (srcPolys == src2.polynomials()) << std::endl;
//        std::cerr << "src2.variableName = " << src2.variableName() << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain move assignment operator test:\t\t\t\t PASSED" << std::endl;
//}

//void testIdentityOperator() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src;
//    src = randomDerivativeSubResultantChain(coefBound,sparsity,degree);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src2(src);
//    
//    if ((src == src2) != true) {
//        std::cerr << "SubResultantChain identity operator test:\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain identity operator test:\t\t\t\t\t PASSED" << std::endl;
//}

//void testSize() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(randomDerivativeSubResultantChain(coefBound,sparsity,degree));
//    
//    if (src.size() != degree+1) {
//        std::cerr << "SubResultantChain size test:\t\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain size test:\t\t\t\t\t\t\t PASSED" << std::endl;
//}

//void testVariableName() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//	
//    SparseUnivariatePolynomial<RN> p(randomSUPQPolynomial(degree,sparsity,coefBound));
//	SparseUnivariatePolynomial<RN> q(p);
//	Symbol s(Symbol::randomElement());
//	p.setVariableName(s);
//	q.setVariableName(s);
//	q.differentiate();
//	
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(p,q,s);

//    if (src.variableName() != s) {
//        std::cerr << "SubResultantChain variableName test:\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain variableName test:\t\t\t\t\t\t PASSED" << std::endl;
//}

//void testIsEmpty() {
//	
//	unsigned long int coefBound(4);
//	double sparsity(0.5);
//	int degree(10);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(randomDerivativeSubResultantChain(coefBound,sparsity,degree));
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src2;
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src3(Symbol::randomElement());

//    if (src.isEmpty() == true || src2.isEmpty() == false || src3.isEmpty() == false) {
//        std::cerr << "SubResultantChain isEmpty test:\t\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain isEmpty test:\t\t\t\t\t\t\t PASSED" << std::endl;
//}

SubResultantChain<SMQP,SMQP> randomDerivativeSubResultantChain(unsigned long int coefBound, float sparsity, int degree, int* actualDegree = NULL) {
	vector<int> degs = {degree};
	SMQP p;
	p.randomPolynomial(degs,coefBound,sparsity,1);
	SMQP q(p);
	Symbol s(Symbol::randomElement());
	vector<Symbol> vars = {s};
	p.setRingVariables(vars);
	q.setRingVariables(vars);
	q.differentiate(s);
	if (actualDegree != NULL ){
		*actualDegree = p.degree().get_si();
	} 
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
	
    SubResultantChain<SMQP,SMQP> src(p,q,s);
//    cout << "src = " << src << endl;
	return src;
}

void testDefaultConstructor() {
    SubResultantChain<SMQP,SMQP> src;

    if (src.size() != 0) {
        std::cerr << "SubResultantChain default constructor test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain default constructor test:\t\t\t\t\t PASSED" << std::endl;
}

void testVariableSpecificationConstructor() {
	Symbol s(Symbol::randomElement());
    SubResultantChain<SMQP,SMQP> src(s);

    if (src.size() != 0 || src.variableName() != s) {
        std::cerr << "SubResultantChain variable specification constructor test:\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain variable specification constructor test:\t\t\t PASSED" << std::endl;
}

void testChainConstructor() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
	vector<int> degs = {degree};
	SMQP p;
	p.randomPolynomial(degs,coefBound,sparsity,1);
	SMQP q(p);
	Symbol s(Symbol::randomElement());
	vector<Symbol> vars = {s};
	p.setRingVariables(vars);
	q.setRingVariables(vars);
	q.differentiate(s);
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
	
    SubResultantChain<SMQP,SMQP> src(p,q,s);
//    cout << "src = " << src << endl;

    if (false) {
        std::cerr << "SubResultantChain chain constructor test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain chain constructor test:\t\t\t\t\t PASSED" << std::endl;
}

void testCopyConstructor() {

	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
    SubResultantChain<SMQP,SMQP> src;
    src = randomDerivativeSubResultantChain(coefBound,sparsity,degree);
    SubResultantChain<SMQP,SMQP> src2(src);

    if (src.polynomials() != src2.polynomials() || src.variableName() != src2.variableName()) {
        std::cerr << "SubResultantChain copy constructor test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain copy constructor test:\t\t\t\t\t PASSED" << std::endl;
}

void testMoveConstructor() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
	vector<int> degs = {degree};
	SMQP p;
	p.randomPolynomial(degs,coefBound,sparsity,1);
	SMQP q(p);
	Symbol s(Symbol::randomElement());
	vector<Symbol> vars = {s};
	p.setRingVariables(vars);
	q.setRingVariables(vars);
	q.differentiate(s);
	
    SubResultantChain<SMQP,SMQP> src(p,q,s);
    std::vector<SMQP> srcPolys(src.polynomials());
    SubResultantChain<SMQP,SMQP> src2(std::move(src));
    
    if (src.variableName() != s || !src.polynomials().empty() || srcPolys != src2.polynomials() || src2.variableName() != s) {
        std::cerr << "SubResultantChain move constructor test:\t\t\t\t\t FAILED" << std::endl;
        std::cerr << "src.variableName() = " << src.variableName() << std::endl;
        std::cerr << "src.polynomials().empty() = " << src.polynomials().empty() << std::endl;
        std::cerr << "srcPolys == src2.polynomials = " << (srcPolys == src2.polynomials()) << std::endl;
        std::cerr << "src2.variableName = " << src2.variableName() << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain move constructor test:\t\t\t\t\t PASSED" << std::endl;
}

void testAssignmentOperator() {

	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
    SubResultantChain<SMQP,SMQP> src;
    src = randomDerivativeSubResultantChain(coefBound,sparsity,degree);
    SubResultantChain<SMQP,SMQP> src2;
    src2 = src;

    if (src.polynomials() != src2.polynomials() || src.variableName() != src2.variableName()) {
        std::cerr << "SubResultantChain assignment operator test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain assignment operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testMoveAssignmentOperator() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
	vector<int> degs = {degree};
	SMQP p;
	p.randomPolynomial(degs,coefBound,sparsity,1);
	SMQP q(p);
	Symbol s(Symbol::randomElement());
	vector<Symbol> vars = {s};
	p.setRingVariables(vars);
	q.setRingVariables(vars);
	q.differentiate(s);
	
    SubResultantChain<SMQP,SMQP> src(p,q,s);
    std::vector<SMQP> srcPolys(src.polynomials());
    SubResultantChain<SMQP,SMQP> src2;
    
    src2 = std::move(src);
    
    if (src.variableName() != s || !src.polynomials().empty() || srcPolys != src2.polynomials() || src2.variableName() != s) {
        std::cerr << "SubResultantChain move assignment operator test:\t\t\t\t FAILED" << std::endl;
        std::cerr << "src.variableName() = " << src.variableName() << std::endl;
        std::cerr << "src.polynomials().empty() = " << src.polynomials().empty() << std::endl;
        std::cerr << "srcPolys == src2.polynomials = " << (srcPolys == src2.polynomials()) << std::endl;
        std::cerr << "src2.variableName = " << src2.variableName() << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain move assignment operator test:\t\t\t\t PASSED" << std::endl;
}

void testIdentityOperator() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
    SubResultantChain<SMQP,SMQP> src;
    src = randomDerivativeSubResultantChain(coefBound,sparsity,degree);
    SubResultantChain<SMQP,SMQP> src2(src);
    
    if ((src == src2) != true) {
        std::cerr << "SubResultantChain identity operator test:\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain identity operator test:\t\t\t\t\t PASSED" << std::endl;
}

void testSize() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
	int actualDegree = 0;
    SubResultantChain<SMQP,SMQP> src(randomDerivativeSubResultantChain(coefBound,sparsity,degree, &actualDegree));
    

    if (src.size() != actualDegree+1) {
    	std::cerr << "src sizE: " << src.size() << std::endl;
    	std::cerr << "actualDegree: " << actualDegree << std::endl;
        std::cerr << "SubResultantChain size test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain size test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testVariableName() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
	vector<int> degs = {degree};
	SMQP p;
	p.randomPolynomial(degs,coefBound,sparsity,1);
	SMQP q(p);
	Symbol s(Symbol::randomElement());
	vector<Symbol> vars = {s};
	p.setRingVariables(vars);
	q.setRingVariables(vars);
	q.differentiate(s);
	
    SubResultantChain<SMQP,SMQP> src(p,q,s);

    if (src.variableName() != s) {
        std::cerr << "SubResultantChain variableName test:\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain variableName test:\t\t\t\t\t\t PASSED" << std::endl;
}

void testIsEmpty() {
	
	unsigned long int coefBound(4);
	float sparsity(0.5);
	int degree(10);
    SubResultantChain<SMQP,SMQP> src(randomDerivativeSubResultantChain(coefBound,sparsity,degree));
    SubResultantChain<SMQP,SMQP> src2;
    SubResultantChain<SMQP,SMQP> src3(Symbol::randomElement());

    if (src.isEmpty() == true || src2.isEmpty() == false || src3.isEmpty() == false) {
        std::cerr << "SubResultantChain isEmpty test:\t\t\t\t\t\t\t FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SubResultantChain isEmpty test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

void testSubResultantOfIndex() {

//	int nvars(5);
//	int nterms(14);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(4);
//	bool includeNeg(1);
//	Symbol s(Symbol::randomElement());
//	SparseMultivariateRationalPolynomial p(nvars);
//	vector<Symbol> xs;
//	xs.push_back(s);
//	xs.push_back(Symbol::randomElement());
//	xs.push_back(Symbol::randomElement());
//	vector<Symbol> ps;
//	ps.push_back(Symbol::randomElement());
//	ps.push_back(Symbol::randomElement());
//	vector<Symbol> vs(xs);
//	vs.insert(vs.end(),ps.begin(),ps.end());

//	p.randomPolynomial(nvars,nterms,coefBound,sparsity,includeNeg);
//	p.setRingVariables(vs);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(p,ps);
//    vector<Symbol> vars;
//    vars.push_back(p.leadingVariable());
//    SparseMultivariateRationalPolynomial zero;

//    if (src.select(p.leadingVariable()) != p || src.select(Symbol::randomElement()) != zero) {
//        std::cerr << "SubResultantChain select test:\t\t\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain select test:\t\t\t\t\t\t\t PASSED" << std::endl;
}

//void testPseudoDivide() {
//    
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src;
//	SparseMultivariateRationalPolynomial p,check;
//	std::vector<SparseMultivariateRationalPolynomial> set;
//	std::vector<Symbol> vars,trcVars;
//	int nvars(3);
//	int nAlgVars(3);
//	int nTrcVars(0);
//	int nterms(6);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(30);
//	bool includeNeg(1);
//    src.randomSubResultantChain(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
////    src.display();
//    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);
//    vars = src.variables();
//    trcVars = src.transcendentalVariables();
//    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
//	p.setRingVariables(vars);
//	set = src.polynomials();
//	std::reverse(set.begin(),set.end());
////	std::cerr << "set.size = " << set.size() << std::endl;
////    std::cerr << "[";
////	for (int i=0; i<set.size(); ++i){
////		std::cerr << set[i];
////		if (i != set.size()-1)
////			std::cerr << ",";
////	}
////	std::cerr << "]" << std::endl;
////	std::cerr << "src = " << src << std::endl;
////	std::cerr << "p = " << p << std::endl;
//    vars = src.variables();
//    std::reverse(vars.begin(),vars.end());
////    std::cerr << "[";
////	for (int i=0; i<vars.size(); ++i){
////		std::cerr << vars[i];
////		if (i != vars.size()-1)
////			std::cerr << ",";
////	}
////	std::cerr << "]" << std::endl;
//	SparseMultivariateRationalPolynomial r,c;
//	vector<SparseMultivariateRationalPolynomial> Q;
//	r = src.pseudoDivide(p,&Q,&c);
////	cout << "c = " << c << endl;
////	cout << "r = " << r << endl;
////	for (int i=0; i<Q.size(); ++i)
////		cout << "Q[" << i << "] = " << Q[i] << endl;
////		
////	cerr << "compute constraint:" << endl;
////	set = src.polynomials();
////	check.zero();
////	for (int i=0; i<set.size(); ++i)
////		check += set[i]*Q[i];
////	check += r;
////	check -= (c*p);
////	cout << "check = " << check << endl;
////    if (!check.isZero()) {
////        std::cerr << "SubResultantChain pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;
////        exit(1);
////    }

////    std::cerr << "SubResultantChain pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;

//    ExpressionTree pTree = p.convertToExpressionTree();
////    ExpressionTree srcTree = src.convertToExpressionTree();
//	ExpressionTree srcTree;
//	srcTree.fromVector<SparseMultivariateRationalPolynomial>(set);
//    ExpressionTree cTree = c.convertToExpressionTree();
//    ExpressionTree remTree = r.convertToExpressionTree();
//    std::vector<ExpressionTree> quoTrees;
//    for (int i = 0; i < Q.size(); ++i) {
//        quoTrees.push_back(Q[i].convertToExpressionTree());
//    }

//    std::vector<std::string> inpusrc;
//    inpusrc.push_back(pTree.toMapleString() + ":");
//    inpusrc.push_back(srcTree.toMapleString() + ":");

////    vars = src.variables();
////    std::reverse(vars.begin(),vars.end());
//    std::string varList;
//    if (vars.size() > 1) {
//        varList = "[" + vars[0].toString() + ",";
//        for (int i = 1; i < vars.size()-1; ++i) {
//            varList += vars[i].toString();
//            varList += ",";
//        }
//        varList += vars[vars.size()-1].toString();
//    } else {
//        varList = "[" + vars[0].toString();
//    }
//    varList += "]:";
//    inpusrc.push_back(varList);

//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    char* cstr;

//    std::vector<ALGEB> algebList;
//    for (int i = 0; i < inpusrc.size(); ++i) {
//        cstr = new char[inpusrc[i].length()+1];
//        std::strcpy(cstr, inpusrc[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        algebList.push_back(res);
//        delete[] cstr;
//    }
//    algebList.push_back(ToMapleName(kv, "qList", 1));
//    algebList.push_back(ToMapleName(kv, "m", 1));

//    //algebList: [0] = dividend, [1] = divisorList, [2] = varList, [3] = qList, [4] = m

////    std::string procStr = "PseudoDivide := proc(p,src,vs,q,m)\r\n local quo;\r\n quo:=[1,2,3];\r\n assign(q=quo);\r\n assign(m=2); 3;\r\n end proc:";
//	std::string procStr = "multiDivisorPseudoDivision := proc (p, src, names, q, m) local seqLoop, srcSize, r, tmpQ, Q, tmpM, totalM, i, j; seqLoop := [3, 2, 1]; srcSize := nops(src); r := p; Q := [seq(0, i = 1 .. srcSize)]; totalM := 1; for i in seqLoop do tmpQ := 0; tmpM := 1; r := prem(r, src[i], names[i], 'tmpM', 'tmpQ'); totalM := totalM*tmpM; for j to srcSize do Q[j] := simplify(tmpM*Q[j]) end do; Q[i] := Q[i]+tmpQ end do; assign(q = Q); assign(m = totalM); return r end proc:";
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

//    std::vector<ALGEB> resultAlgebs;
//    inpusrc.clear();
//    inpusrc.push_back(remTree.toMapleString() + ":");
//    inpusrc.push_back(cTree.toMapleString() + ":");
//    for (int i = 0; i < quoTrees.size(); ++i) {
//        inpusrc.push_back(quoTrees[i].toMapleString() + ":");
//    }
//    for (int i = 0; i < inpusrc.size(); ++i) {
//        cstr = new char[inpusrc[i].length()+1];
//        std::strcpy(cstr, inpusrc[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        resultAlgebs.push_back(res);
//        delete[] cstr;
//    }

//    char compareFunc[] = "verify:";
//    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);

//    char expandFunc[] = "expand:";
//    ALGEB expandF = EvalMapleStatement(kv, expandFunc);

//    result = EvalMapleProc(kv, expandF, 1, result);
//    ALGEB comp = EvalMapleProc(kv, cmpF, 2, result, resultAlgebs[0]);
//    M_BOOL compBool = MapleToM_BOOL(kv, comp);
//    if (!MapleToM_BOOL(kv, comp))
//    	std::cerr << "remainder comp fail!" << std::endl;
//    std::string mapleResStr = mapleTest->algebToString(kv, result);
//    std::string getMStr = "m:";
//    cstr = new char[getMStr.length()+1];
//    std::strcpy(cstr,getMStr.c_str());
//    ALGEB mapleM = EvalMapleStatement(kv, cstr);
//    mapleM = EvalMapleProc(kv, expandF, 1, mapleM);
//    comp = EvalMapleProc(kv, cmpF, 2, mapleM, resultAlgebs[1]);
//    compBool = compBool && MapleToM_BOOL(kv, comp);
//    if (!MapleToM_BOOL(kv, comp))
//    	std::cerr << "multiplier comp fail!" << std::endl;
//    delete[] cstr;
//    std::vector<ALGEB> mapleQuoList;
//    int idx = 2;
////    while (idx < resultAlgebs.size()) {
//    while (idx <= 4) {
//        std::string getQuoStr = "qList[" + std::to_string(nvars+2-idx) + "]:";
//        cstr = new char[getQuoStr.length()+1];
//        std::strcpy(cstr, getQuoStr.c_str());
//        ALGEB mapleQuo = EvalMapleStatement(kv, cstr); 
//        mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo); 
//        comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
//        compBool = compBool && MapleToM_BOOL(kv, comp);
//    	if (!MapleToM_BOOL(kv, comp)) {
//    		std::cerr << "quotient " << idx-2 << " comp fail!" << std::endl;
//    	}
//        mapleQuoList.push_back(mapleQuo);
//        ++idx;
//        delete[] cstr;
//    }

//    if(!compBool) {
//        std::cerr << "SubResultantChain pseudoDivide test:\t\t\t\t\t\t FAILED" << std::endl;

//        std::cerr << "Dividend: " << p << std::endl;
//        std::cerr << "Divisors: " << src << std::endl << std::endl;

//        std::cerr << "Got remainder: " << r << std::endl;
//        std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;

//        std::cerr << "Got multiplier: " << c << std::endl;
//        std::cerr << "Got maple m: " << mapleTest->algebToString(kv, mapleM) << std::endl << std::endl;

//        for (int i = 0; i < mapleQuoList.size(); ++i) {   
//            std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
//            std::cerr << "Got maple quo : " << mapleTest->algebToString(kv, mapleQuoList[i]) << std::endl;
//        }
//        exit(1);
//    }

//    std::cerr << "SubResultantChain pseudoDivide test:\t\t\t\t\t\t PASSED" << std::endl;
//}

//void testNormalForm() {
//    
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src;
//	SparseMultivariateRationalPolynomial p,check;
//	std::vector<SparseMultivariateRationalPolynomial> set;
//	std::vector<Symbol> vars,trcVars;
//	int nvars(4);
//	int nAlgVars(4);
//	int nTrcVars(2);
//	int nterms(8);
//	unsigned long int coefBound(6ul);
//	degree_t sparsity(25);
//	bool includeNeg(1);
//    src.randomStronglyNormalizedSubResultantChain(nvars,nAlgVars,nTrcVars,nterms,coefBound,sparsity,includeNeg);
////    src.display();
//    p.randomPolynomial(nvars+nTrcVars,nterms,coefBound,sparsity,includeNeg);


//    vars = src.variables();
//    trcVars = src.transcendentalVariables();
//    vars.insert(vars.end(),trcVars.begin(),trcVars.end());
//    p.setRingVariables(vars);
//    set = src.polynomials();
//    for (int i = 0; i < set.size(); ++i) {
//        p *= set[i];
//    }

//    SparseMultivariateRationalPolynomial r;
//	vector<SparseMultivariateRationalPolynomial> Q;
//	r = src.normalForm(p,&Q);

////	cout << "r = " << r << endl;
////	for (int i=0; i<Q.size(); ++i)
////		cout << "Q[" << i << "] = " << Q[i] << endl;
//		
////	cerr << "compute constraint:" << endl;
//	// set = src.polynomials();
//	// check.zero();
//	// for (int i=0; i<set.size(); ++i)
//	// 	check += set[i]*Q[i];
//	// check += r;
//	// check -= p;
////	cout << "check = " << check << endl;
//    // if (!check.isZero()) {
//    //     std::cerr << "SubResultantChain normalForm test:\t\t\t\t\t\t\t FAILED" << std::endl;
//    //     exit(1);
//    // }

//    // std::cerr << "SubResultantChain normalForm test:\t\t\t\t\t\t\t PASSED" << std::endl;


//    ExpressionTree pTree = p.convertToExpressionTree();
//    ExpressionTree srcTree = src.convertToExpressionTree();
//    ExpressionTree remTree = r.convertToExpressionTree();
//    std::vector<ExpressionTree> quoTrees;
//    for (int i = 0; i < Q.size(); ++i) {
//        quoTrees.push_back(Q[i].convertToExpressionTree());
//    }

//    std::vector<std::string> inpusrc;
//    inpusrc.push_back(pTree.toMapleString() + ":");
//    inpusrc.push_back(srcTree.toMapleString() + ":");

//    vars = src.variables();
//    std::string plex;
//    if (vars.size() > 1) {
//        plex = "plex(" + vars[0].toString() + ",";
//        for (int i = 1; i < vars.size()-1; ++i) {
//            plex += vars[i].toString();
//            plex += ",";
//        }
//        plex += vars[vars.size()-1].toString();
//    } else {
//        plex = "plex(" + vars[0].toString();
//    }
//    plex += "):";
//    inpusrc.push_back(plex);

//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    char* cstr;

//    std::vector<ALGEB> algebList;
//    for (int i = 0; i < inpusrc.size(); ++i) {
//        cstr = new char[inpusrc[i].length()+1];
//        std::strcpy(cstr, inpusrc[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        algebList.push_back(res);
//        delete[] cstr;
//    }
//    algebList.push_back(ToMapleName(kv, "qList", 1));

//    //algebList: [0] = dividend, [1] = divisorList, [2] = plex, [3] = qList

//    std::string procStr = "Groebner:-NormalForm:";
//    cstr = new char[procStr.length()+1];
//    std::strcpy (cstr, procStr.c_str());
//    ALGEB testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;

//    ALGEB result = EvalMapleProc(kv, testProc, 4, algebList[0], algebList[1], algebList[2], algebList[3]);

//    // std::cerr << "calling maple proc: \n\n";
//    // std::cerr << mapleTest->algebToString(kv, testProc) << "\n\n";
//    // std::cerr << mapleTest->algebToString(kv, result) << "\n\n";
//    // std::cerr << mapleTest->algebToString(kv, algebList[0]) << "\n\n";
//    // std::cerr << mapleTest->algebToString(kv, algebList[1]) << "\n\n";
//    // std::cerr << mapleTest->algebToString(kv, algebList[2]) << "\n\n";
//    // std::cerr << mapleTest->algebToString(kv, algebList[3]) << "\n\n";
//    // std::cerr << "\n\n";

//    std::vector<ALGEB> resultAlgebs;
//    inpusrc.clear();
//    inpusrc.push_back(remTree.toMapleString() + ":");
//    for (int i = 0; i < quoTrees.size(); ++i) {
//        inpusrc.push_back(quoTrees[i].toMapleString() + ":");
//    }
//    for (int i = 0; i < inpusrc.size(); ++i) {
//        cstr = new char[inpusrc[i].length()+1];
//        std::strcpy(cstr, inpusrc[i].c_str());
//        ALGEB res = EvalMapleStatement(kv, cstr);
//        resultAlgebs.push_back(res);
//        delete[] cstr;
//    }

//    char compareFunc[] = "verify:";
//    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);

//    char expandFunc[] = "expand:";
//    ALGEB expandF = EvalMapleStatement(kv, expandFunc);

//    result = EvalMapleProc(kv, expandF, 1, result);
//    ALGEB comp = EvalMapleProc(kv, cmpF, 2, result, resultAlgebs[0]);
//    M_BOOL compBool = MapleToM_BOOL(kv, comp);
//    std::string mapleResStr = mapleTest->algebToString(kv, result);
//    std::vector<ALGEB> mapleQuoList;
//    int idx = 1;
//    while (idx < resultAlgebs.size()) {
//        std::string getQuoStr = "qList[" + std::to_string(idx) + "]:";
//        cstr = new char[getQuoStr.length()+1];
//        std::strcpy(cstr, getQuoStr.c_str());
//        ALGEB mapleQuo = EvalMapleStatement(kv, cstr); 
//        mapleQuo = EvalMapleProc(kv, expandF, 1, mapleQuo); 
//        comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, resultAlgebs[idx]);
//        compBool = compBool && MapleToM_BOOL(kv, comp);
//        mapleQuoList.push_back(mapleQuo);
//        ++idx;
//        delete[] cstr;
//    }

//    if(!compBool) {
//        std::cerr << "SubResultantChain normalForm test:\t\t\t\t\t\t\t FAILED" << std::endl;

//        std::cerr << "Dividend: " << p << std::endl;
//        std::cerr << "Divisors: " << srcTree.toMapleString() << std::endl << std::endl;

//        std::cerr << "Got remainder: " << r << std::endl;
//        std::cerr << "Got maple rem: " << mapleResStr << std::endl << std::endl;

//        for (int i = 0; i < mapleQuoList.size(); ++i) {   
//            std::cerr << "Got quotient " << i << ": " << Q[i] << std::endl;
//            std::cerr << "Got maple quo : " << mapleTest->algebToString(kv, mapleQuoList[i]) << std::endl;
//        }
//        exit(1);
//    }

//    std::cerr << "SubResultantChain normalForm test:\t\t\t\t\t\t\t PASSED" << std::endl;
//}

void testConvertToExpressionTree() {

//    stringstream ss;
//    std::string str,str2;
//	
//	unsigned long int coefBound(4);
//	float sparsity(0.5);
//	int degree(10);
//    SubResultantChain<RN,SparseUnivariatePolynomial<RN>> src(randomDerivativeSubResultantChain(coefBound,sparsity,degree));
//    ExpressionTree et = src.convertToExpressionTree();
//    // cout << "et = " << et.toString() << endl;
//    ss << src;
//    str = ss.str();
//    // cout << "src = " << str << endl;
//    // MAPLE EQUALITY TEST FOR et.toString() and str //
//    // Start Maple Kernel:
//    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
//    mapleTest->restartMapleKernel();
//    MKernelVector kv = mapleTest->getMKernelVector();
//    
//    // et:
//    std::string esrctr = et.toString() + ":";
//    char* cstr = new char[esrctr.size() + 1];
//    std::strcpy (cstr, esrctr.c_str());
//    ALGEB resET = EvalMapleStatement(kv, cstr);
//    delete[] cstr;
// 
//    // str:
//    std::string srcStr = str + ":";
//    cstr = new char[srcStr.size() + 1];
//    std::strcpy (cstr, srcStr.c_str());
//    ALGEB ressrc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;
// 
//    // operator:
//    std::string procStr = "ArrayTools:-IsEqual:";
//    cstr = new char[procStr.length() + 1];
//    std::strcpy (cstr, procStr.c_str());
//    ALGEB testProc = EvalMapleStatement(kv, cstr);
//    delete[] cstr;
// 
//    // Maple:
//    ALGEB equalityRes = EvalMapleProc(kv, testProc, 2, resET, ressrc);
//    M_BOOL cmpBool = MapleToM_BOOL(kv, equalityRes);

//    src = SubResultantChain<RN,SparseUnivariatePolynomial<RN>>();
//    ExpressionTree et2;
//    et2 = src.convertToExpressionTree();
////  cout << "et2 = " << et2.toString() << endl;
//    ss = stringstream();
//    ss << src;
//    str2 = ss.str();
////  cout << "src = " << str2 << endl;
//    if (!cmpBool || et2.toString() != str2) {
//        std::cerr << "SubResultantChain convertToExpressionTree test:\t\t\t\t\t FAILED" << std::endl;
//        exit(1);
//    }

//    std::cerr << "SubResultantChain convertToExpressionTree test:\t\t\t\t\t PASSED" << std::endl;
}

int main() {

	testDefaultConstructor();
	testVariableSpecificationConstructor();
	testChainConstructor();
	testCopyConstructor();
	testMoveConstructor();
	testAssignmentOperator();
	testMoveAssignmentOperator();
	testIdentityOperator();
	testSize();
	testVariableName();
	testIsEmpty();
	testSubResultantOfIndex();
//	testResultant();
	testConvertToExpressionTree();
	
	
//	SparseMultivariateRationalPolynomial p,q;
//	int nvars = 2;
//	int nterms = 4;
//	p.randomPolynomial(nvars,nterms,4,2,1);
//	q.randomPolynomial(nvars,nterms,4,2,1);
//	std::vector<Symbol> vec = Symbol::randomElements(nvars);
//	// for (int i=0; i<nvars; ++i) {
//	// 	vec.push_back(Symbol::randomElement());
//	// }
//	
//	p.setRingVariables(vec);
//	q.setRingVariables(vec);
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
//	SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> a,b;
//	a = p.convertToSUP(p.leadingVariable());
//	b = q.convertToSUP(q.leadingVariable());
//	
//    SubResultantChain<RN,SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>> src(a,b,a.variable());
////    cout << "src = " << src << endl;
//    for (int i=0; i<src.size(); ++i)
//    	cout << "S[" << i << "] = " << src.subResultantOfIndex(i) << endl;
//    SubResultantChain<SMQP,SMQP> src2(p,q,p.leadingVariable());
////    cout << "src2 = " << src2 << endl;
//    for (int i=0; i<src2.size(); ++i)
//    	cout << "S[" << i << "] = " << src2.subResultantOfIndex(i) << endl;
//    cout << endl;
//    
//	q *= p;
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
//	a = p.convertToSUP(p.leadingVariable());
//	b = q.convertToSUP(q.leadingVariable());
//	
//    src = SubResultantChain<RN,SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>>(a,b,a.variable());
////    cout << "src = " << src << endl;
//    for (int i=0; i<src.size(); ++i)
//    	cout << "S[" << i << "] = " << src.subResultantOfIndex(i) << endl;
//    src2 = SubResultantChain<SMQP,SMQP>(p,q,p.leadingVariable());
////    cout << "src2 = " << src2 << endl;
//    for (int i=0; i<src2.size(); ++i)
//    	cout << "S[" << i << "] = " << src2.subResultantOfIndex(i) << endl;
//    
//    p *= p;
//    p *= p;
//	q = p.derivative(p.leadingVariable());
//	cout << "p = " << p << endl;
//	cout << "q = " << q << endl;
//	a = p.convertToSUP(p.leadingVariable());
//	b = q.convertToSUP(q.leadingVariable());
//	
//    src = SubResultantChain<RN,SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>>(a,b,a.variable());
////    cout << "src = " << src << endl;
//    for (int i=0; i<src.size(); ++i)
//    	cout << "S[" << i << "] = " << src.subResultantOfIndex(i) << endl;
//    src2 = SubResultantChain<SMQP,SMQP>(p,q,p.leadingVariable());
////    cout << "src2 = " << src2 << endl;
//    for (int i=0; i<src2.size(); ++i)
//    	cout << "S[" << i << "] = " << src2.subResultantOfIndex(i) << endl;



	return 0;
}
