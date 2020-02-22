#include "../../include/bpas.h"
#include "../../include/Utils/Unix_Timer.h"
#include <iomanip>
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include "../MapleTestTool/MapleTestTool.hpp"

long nvar = 3;
long numTerms = 10;
long coefBound = 5ul;
degree_t sparsity = 10;
int includeNeg = 1;

///////////////////////
// - For each operator, check the following combinations: 
//      - SMQP of same nvar
//      - SMQP of different nvar
//      - SMQP that is zero
//      - SMQP that is a constant
//      - SMQP that is a constant and is nvar == 0
//
///////////////////////

void testDefaultConstructor() {
    SparseMultivariateRationalPolynomial p;

    if (!p.isZero() || p.numberOfRingVariables() != 0) {
        std::cerr << "SMQP default constructor test: FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SMQP default constructor test: PASSED" << std::endl;    
}

void testNvarConstructor(){
    SparseMultivariateRationalPolynomial p(0);

    if (!p.isZero() || p.numberOfRingVariables() != 0) {
        std::cerr << "SMQP nvar constructor test: FAILED" << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(nvar);
    if (!p.isZero() || p.numberOfRingVariables() != nvar) {
        std::cerr << "SMQP nvar constructor test: FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SMQP nvar constructor test: PASSED" << std::endl;
}

void testIdentityConstructor() {
    SparseMultivariateRationalPolynomial p("x");

    ExpressionTree pTree = p.convertToExpressionTree();
    pTree -= ExpressionTree(new ExprTreeNode(Symbol("x")));

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP identity constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial("xas_12312");
    pTree = p.convertToExpressionTree();
    pTree -= ExpressionTree(new ExprTreeNode(Symbol("xas_12312")));

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP identity constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP identity constructor test: PASSED" << std::endl;
}

void testStringConstructor() {
    std::string str = "3*x^4*y^3 + 6*x^2*y + 3*z";
    SparseMultivariateRationalPolynomial p(str);

    std::stringstream ss;
    ss << p;

    if (str != ss.str()) {
        std::cerr << "SMQP string constructor test: FAILED" << std::endl;
        std::cerr << "Expected: " << str << " but got: " << ss.str() << std::endl;
        exit(1);
    }

    std::cerr << "SMQP string constructor test: PASSED" << std::endl;
}

void testCopyConstructor() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial q = p;

    ExpressionTree pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    
    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP copy constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP copy constructor test: PASSED" << std::endl;
}

void testMoveConstructor() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();

    mpq_class val = mpq_class(rand() % coefBound + 1);
    SparseMultivariateRationalPolynomial q = p + val; 

    pTree += ExpressionTree(new ExprTreeNode(mpq_class(val)));

    pTree -= q.convertToExpressionTree();

    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP move constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP move constructor test: PASSED" << std::endl;   
}

void testSMZPConstructor() {
    SparseMultivariateIntegerPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    ExpressionTree pTree = p.convertToExpressionTree();

    SparseMultivariateRationalPolynomial q(p);

    pTree -= q.convertToExpressionTree();

    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP SMZP constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP SMZP constructor test: PASSED" << std::endl;

}

void testIntegerConstructor() {
    Integer i(47);
    SparseMultivariateRationalPolynomial p(i);

    ExpressionTree pTree = p.convertToExpressionTree();
    pTree -= i.convertToExpressionTree();

    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP integer constructor test: FAILED" << std::endl;
        std::cerr << "Expected: " << i << " but got: " << p << std::endl;
        exit(1);
    }

    std::cerr << "SMQP integer consturctor test: PASSED" << std::endl;
}

void testRationalConstructor() {
    RationalNumber i(47, 23);
    SparseMultivariateRationalPolynomial p(i);

    ExpressionTree pTree = p.convertToExpressionTree();
    pTree -= i.convertToExpressionTree();

    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP rational number constructor test: FAILED" << std::endl;
        std::cerr << "Expected: " << i << " but got: " << p << std::endl;
        exit(1);
    }

    std::cerr << "SMQP rational number consturctor test: PASSED" << std::endl;
}
ExpressionTree DUQPToExpressionTree(DenseUnivariateRationalPolynomial duqp) {

    Integer d = duqp.degree();
    Symbol var = duqp.variable();
    ExpressionTree prev;
    bool first = 1;

    mpq_class coef;
    for(int i = d.get_si(); i >=0; --i) {
        coef = duqp.coefficient(i).get_mpq();
        if (coef == 0) {
            continue;
        }

        ExpressionTree coefTree = ExpressionTree(new ExprTreeNode(coef));
        ExpressionTree monomial = ExpressionTree(new ExprTreeNode(var));
        monomial ^= ExpressionTree(new ExprTreeNode(i));
        coefTree *= monomial;

        if (first) {
            prev = coefTree;
            first = 0;
        } else {
            prev += coefTree;
        }
    }

    return prev;
}

void testDUQPConstructor() {
    unsigned long maxDegree = sparsity*nvar*numTerms;
    DenseUnivariateRationalPolynomial duqp(maxDegree);

    unsigned long curD = 0;
    for (long i = 0; i < numTerms; ++i) {
        if (nvar > 0) {
            curD += rand() % sparsity;
        }
        duqp.setCoefficient(curD, mpq_class(rand() % coefBound + 1));
    }

    SparseMultivariateRationalPolynomial p(duqp);

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree dTree = DUQPToExpressionTree(duqp);

    pTree -= dTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if(mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SMQP DUQP consturctor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP DUQP constructor test: PASSED" << std::endl;   
}

//NOTE this constructor assumed that convertToSUP() works.
ExpressionTree SUPconvertToExpressionTree(SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> sup); //foward declare
void testSUPConstructor() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> sup = p.convertToSUP(p.leadingVariable());

    SparseMultivariateRationalPolynomial q(sup);

    ExpressionTree qTree = q.convertToExpressionTree();
    qTree -= SUPconvertToExpressionTree(sup);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;

    if(mapleTest->testIfZero(qTree, &retErr) == 0) {
        std::cerr << "SMQP SUP<SMQP> constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP SUP<SMQP> constructor test: PASSED" << std::endl;   
}


void testIsZero () {
    SparseMultivariateRationalPolynomial p1;
    SparseMultivariateRationalPolynomial p2;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p2 = p1;
    if (p1.isZero() == 0 && (p1-p2).isZero() == 1) {
        std::cerr << "SMQP isZero() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP isZero() test: FAILED" << std::endl;
    exit(1);
}

void testZero() {
    SparseMultivariateRationalPolynomial p1;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    bool before = p1.isZero();
    if (nvar == 0) {
        before = (p1.leadingCoefficient() == 0);
    }
    p1.zero();
    bool after = p1.isZero();
    if (!before && after) {
        std::cerr << "SMQP zero() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP zero() test: FAILED" << std::endl;
    exit(1);
}

void testIsOne() {
    SparseMultivariateRationalPolynomial p1;
    SparseMultivariateRationalPolynomial p2;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p2 = p1;
    p1 += mpq_class(1);

    if (p1.isOne() == 0 && (p1-p2).isOne() == 1) {
        std::cerr << "SMQP isOne() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP isOne() test: FAILED" << std::endl;
    exit(1);
}

void testOne() {
    SparseMultivariateRationalPolynomial p1;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    bool before = (p1.isOne() && !(numTerms == 1 && p1.leadingCoefficient() == 1));
    p1.one();
    bool after = p1.isOne();

    if (!before && after) {
        std::cerr << "SMQP one() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP one() test: FAILED" << std::endl;
    exit(1);
}

void testIsNegativeOne() {
    SparseMultivariateRationalPolynomial p1;
    SparseMultivariateRationalPolynomial p2;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p2 = p1;
    p1 += mpq_class(-1);
    if (p1.isNegativeOne() == 0 && (p1-p2).isNegativeOne() == 1) {
        std::cerr << "SMQP isNegativeOne() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP isNegativeOne() test: FAILED" << std::endl;
    exit(1);   
}

void testNegativeOne() {
    SparseMultivariateRationalPolynomial p1;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    bool before = (p1.isNegativeOne() && !(numTerms == 1 && p1.leadingCoefficient() == -1));
    p1.negativeOne();
    bool after = p1.isNegativeOne();

    if (!before && after) {
        std::cerr << "SMQP negativeOne() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP negativeOne() test: FAILED" << std::endl;
    exit(1);
}

void testIsConstant() {
    SparseMultivariateRationalPolynomial p1;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    SparseMultivariateRationalPolynomial p2 = p1;
    p2 += mpq_class(rand() % coefBound + 1);
    if ((numTerms <= 1 || p1.isConstant() == 0) && (p1-p2).isConstant()) {
        std::cerr << "SMQP isConstant() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP isConstant() test: FAILED" << std::endl;
    exit(1);
}

void testUnitCanonical() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    SparseMultivariateRationalPolynomial lc = p.leadingCoefficient();
    while(lc.isOne()) {
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        lc = p.leadingCoefficient();
    }

    SparseMultivariateRationalPolynomial u;
    SparseMultivariateRationalPolynomial v;
    SparseMultivariateRationalPolynomial unit = p.unitCanonical(&u, &v);

    if (unit.leadingCoefficient() != RationalNumber(1)) {
        std::cerr << "SMQP unitCanonical() test: FAILED" << std::endl;
        std::cerr << "Unit canonical's leading coefficient should be 1!" << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "unitCanonical: " << unit << std::endl;
        exit(1);
    }

    if (!(lc*u).isOne()) {
        std::cerr << "SMQP unitCanonical() test: FAILED" << std::endl;
        std::cerr << "lc(p)*u should equal 1 but got: " << (lc * u) << std::endl;
        exit(1);
    }

    if (!(u*v).isOne()) {
        std::cerr << "SMQP unitCanonical() test: FAILED" << std::endl;
        std::cerr << "v should be the inverse of u but their product is: " << (u*v) << std::endl;
        exit(1); 
    }

    std::cerr << "SMQP unitCanonical() test: PASSED" << std::endl;
}


void testCopyAssignment() {
    SparseMultivariateRationalPolynomial p1;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    SparseMultivariateRationalPolynomial p2 = p1;
    if (p1 != p2) {
        std::cerr << "SMQP copy assignment test: FAILED" << std::endl;
        std::cerr << "Expected equal but got: " << p1 << " and " << p2 << std::endl;
        return;
    }

    RationalNumber it(13, 11);
    p1 = it;

    ExpressionTree pTree = p1.convertToExpressionTree();
    pTree -= it.convertToExpressionTree();
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP copy assignment from RationalNumber test: FAILED" << std::endl; 
        std::cerr << "Integer: " << it << " SMQP: " << p1;
        exit(1);
    }

    std::cerr << "SMQP copy assignment test: PASSED" << std::endl;
}


void testMoveAssignment() {
    SparseMultivariateRationalPolynomial p1;
    p1.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    SparseMultivariateRationalPolynomial p2;
    p2.one();

    ExpressionTree pTree = p1.convertToExpressionTree();
    pTree += p2.convertToExpressionTree();

    p2 = p1 + p2;

    pTree -= p2.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP move assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP move assignment test: PASSED" << std::endl; 
}

void testSMQPAddition() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    //      - SMQP of same nvar
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial sum = p+q;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree sumTree = sum.convertToExpressionTree();
    pTree += qTree;
    pTree -= sumTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition test1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP of different nvar
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    sum = p+q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition test2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    q.zero();
    sum = p+q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition test3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    sum = p + q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition test4: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    sum = p + q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition test5: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp addition test: PASSED" << std::endl;
}

void testSMQPAdditionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree qTree = q.convertToExpressionTree();

    p += q;
    ExpressionTree sumTree = p.convertToExpressionTree();
    pTree += qTree;
    pTree -= sumTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SMQP of different nvar
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    pTree = p.convertToExpressionTree();
    p+=q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    q.zero();
    pTree = p.convertToExpressionTree();
    p+=q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p += q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p += q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    std::cout << "SMQP smqp addition assignment test: PASSED" << std::endl;

}

void testUnaryNegative() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial q = -p;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();

    pTree += qTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP unary negative test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //      - SMQP that is zero
    p.zero();
    q = -p;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP unary negative test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SMQP that is a constant
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    q = -p;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP unary negative test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    p = SparseMultivariateRationalPolynomial(0);
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    q = -p;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP unary negative test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP unary negative test: PASSED" << std::endl;
}

void testSMQPSubtraction() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial diff = p-q;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree diffTree = diff.convertToExpressionTree();
    pTree -= qTree;
    pTree -= diffTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
 
    //      - SMQP of different nvar
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    diff = p-q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    q.zero();
    diff = p-q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    diff = p - q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    diff = p - q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp subtraction test: PASSED" << std::endl;
}

void testSMQPSubtractionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree qTree = q.convertToExpressionTree();

    p -= q;
    ExpressionTree diffTree = p.convertToExpressionTree();
    pTree -= qTree;
    pTree -= diffTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction assignment test 1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

        //      - SMQP of different nvar
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    pTree = p.convertToExpressionTree();
    p-=q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();



    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction assignment test 2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    q.zero();
    pTree = p.convertToExpressionTree();
    p-=q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction assignment test 3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p -= q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction assignment test 4: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p -= q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp subtraction assignment test 5: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp subtraction assignment test: PASSED" << std::endl;  
}

void testSMQPMultiplication() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial prod = p*q;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree prodTree = prod.convertToExpressionTree();
    pTree *= qTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP of different nvar
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    prod = p*q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    q.zero();
    prod= p*q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    prod = p * q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    prod = p * q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }




    std::cout << "SMQP smqp multiplication test: PASSED" << std::endl;
}

void testSMQPMultiplicationAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    ExpressionTree qTree = q.convertToExpressionTree();

    p *= q;
    ExpressionTree prodTree = p.convertToExpressionTree();
    pTree *= qTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //      - SMQP of different nvar
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    pTree = p.convertToExpressionTree();
    p*=q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    q.zero();
    pTree = p.convertToExpressionTree();
    p*=q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p *= q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is a constant and is nvar == 0    
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p *= q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp multiplication assignment test: PASSED" << std::endl;   
}

void testSMQPDivision() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial q = p;
    p *= p;
    p *= p;
    p *= p;

    SparseMultivariateRationalPolynomial quo = p/q;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTree = quo.convertToExpressionTree();

    //Due to the way the divide operation works (by returning the value of
    //interest in a new variable) we need some custom code to handle the 
    //variable naming and actually retrieving the value. Not just the proc return value.

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();


    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(qTree.toMapleString() + ":");
    inputs.push_back(quoTree.toMapleString() + ":");

    std::string procStr = "divide:";
    char* cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q1", 1));

    ALGEB result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);

    char compareFunc[] = "verify:";
    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);
    ALGEB comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);

    M_BOOL compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SMQP of different nvar
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    p *= q;
    p *= q;
    p *= q;

    quo = p/q;
    
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(quo.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q2", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SMQP that is zero
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p.zero();
    
    quo = p/q;
    
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(quo.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q3", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SMQP that is a constant
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);

    quo = p/q; 
    
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(quo.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q4", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SMQP that is a constant and is nvar == 0

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);

    quo = p/q; 

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(quo.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q5", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp division test: PASSED" << std::endl;
}

void testSMQPDivisionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    SparseMultivariateRationalPolynomial q = p;
    p *= p;
    p *= p;
    p *= p;
    ExpressionTree pTree = p.convertToExpressionTree();

    p /= q;

    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTree = p.convertToExpressionTree();

    //Due to the way the divide operation works (by returning the value of
    //interest in a new variable) we need some custom code to handle the 
    //variable naming and actually retrieving the value. Not just the proc return value.

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    MKernelVector kv = mapleTest->getMKernelVector();

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(qTree.toMapleString() + ":");
    inputs.push_back(quoTree.toMapleString() + ":");
    
    std::string procStr = "divide:";
    char* cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q6", 1));

    ALGEB result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);

    char compareFunc[] = "verify:";
    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);
    ALGEB comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);

    M_BOOL compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division assignment test: FAILED6" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }


//      - SMQP of different nvar
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    p *= q;
    p *= q;
    p *= q;

    pTree = p.convertToExpressionTree();
    p /= q;

    inputs.clear();
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q7", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division assignment test: FAILED7" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

//      - SMQP that is zero
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    p.zero();
    
    pTree = p.convertToExpressionTree();

    p/=q;

    inputs.clear();
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q8", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division assignment test: FAILED8" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

//      - SMQP that is a constant
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p /= q; 

    inputs.clear();
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q9", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division assignment test: FAILED9" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

//      - SMQP that is a constant and is nvar == 0

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p /= q; 

    inputs.clear();
    inputs.push_back(pTree.toMapleString() + ":");
    inputs.push_back(q.convertToExpressionTree().toMapleString() + ":");
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q10", 1));

    result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[3]);
    comp = EvalMapleProc(kv, cmpF, 2, algebList[3], algebList[2]);
    compBool = MapleToM_BOOL(kv, comp);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp division assignment test: FAILED10" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp division assignment test: PASSED" << std::endl;
}

void testExponentiation() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    int exp = (rand() % 10) + 1; //10 for practical purposes;
    
    SparseMultivariateRationalPolynomial prod = p ^ exp;
    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree prodTree = prod.convertToExpressionTree();

    ExpressionTree constTree = ExpressionTree(new ExprTreeNode(exp));
    pTree ^= constTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation test: FAILED1" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//      - SMQP that is zero
    p.zero();

    prod = p^exp;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation test: FAILED2" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//      - SMQP that is a constant
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    prod = p^exp;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation test: FAILED3" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }    

//      - SMQP that is a constant and is nvar == 0
    p = SparseMultivariateRationalPolynomial(0);
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    prod = p^exp;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation test: FAILED4" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//       - to exponent 0;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    prod = p^0;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    pTree ^= ExpressionTree(new ExprTreeNode(0));
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation test: FAILED5" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp exponentiation test: PASSED" << std::endl;
}

void testExponentiationAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    ExpressionTree pTree = p.convertToExpressionTree();

    int exp = (rand() % 10) + 1;
    p ^= exp;

    ExpressionTree prodTree = p.convertToExpressionTree();
    ExpressionTree constTree = ExpressionTree(new ExprTreeNode(exp));
    pTree ^= constTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation assignment test: FAILED1" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SMQP that is zero
    p.zero();
    pTree = p.convertToExpressionTree();
    p^=exp;
    prodTree = p.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation assignment test: FAILED2" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//      - SMQP that is a constant
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p^=exp;
    prodTree = p.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation assignment test: FAILED3" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }    

//      - SMQP that is a constant and is nvar == 0
    p = SparseMultivariateRationalPolynomial(0);
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p^=exp;
    prodTree = p.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation assignment test: FAILED4" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//       - to exponent 0;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = p.convertToExpressionTree();
    p^=0;
    prodTree = p.convertToExpressionTree();
    pTree ^= ExpressionTree(new ExprTreeNode(0));
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP smqp exponentiation assignment test: FAILED5" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp exponentiation assignment test: PASSED" << std::endl;
}

void testGCD() {

    SparseMultivariateRationalPolynomial p,q;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    std::vector <SparseMultivariateRationalPolynomial> subr;

    // std::cout << "Test1: " << std::endl;
    // std::cout << "p := " << p << std::endl;
    // std::cout << "q := " << q << std::endl;	
    
    SparseMultivariateRationalPolynomial g = p.gcd(q);
    
    std::vector<std::string> inputs;
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(q.convertToExpressionTree().toMapleString());
    
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if (mapleTest->testProcReturn("gcd", inputs, g.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP gcd() test: FAILED1" << std::endl;
        std::cerr << "Expected " << g << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1); 
    }

    //Test non-trivial
    SparseMultivariateRationalPolynomial r;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q = p;
    for (int i = 0; i < numTerms; ++i) {
        r.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        r = r[0];
        p *= r;

        r.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        r = r[0];
        q *= r;
    }
    
    p = p.primitivePart();
    q = q.primitivePart();
    
    // std::cout << "Test2: " << std::endl;
    // std::cout << "p := " << p << std::endl;
    // std::cout << "q := " << q << std::endl;	

    g = p.gcd(q);

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(q.convertToExpressionTree().toMapleString());

    if (mapleTest->testProcReturn("gcd", inputs, g.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP gcd() test: FAILED2" << std::endl;
        std::cerr << "Expected " << g << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1); 
    }

    //Test different nvar
    p.randomPolynomial(nvar, numTerms, coefBound, nvar == 1 ? 3 : sparsity, includeNeg);
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);

    p = p.primitivePart();
    q = q.primitivePart();

    // std::cout << "Test3: " << std::endl;
    // std::cout << "p := " << p << std::endl;
    // std::cout << "q := " << q << std::endl;	
    
    g = p.gcd(q);
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(q.convertToExpressionTree().toMapleString());
    
    if (mapleTest->testProcReturn("gcd", inputs, g.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP gcd() test: FAILED3" << std::endl;
        std::cerr << "Expected " << g << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1); 
    }

    //Test non-trivial
    p.zero();
    q.zero();
    for (int i = 0; i < numTerms; ++i) {
        r.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        r = r[0];
        p += r;

        r.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
        r = r[0];
        q += r;
    }
    
    p = p.primitivePart();
    q = q.primitivePart();

    // std::cout << "Test4: " << std::endl;
    // std::cout << "p := " << p << std::endl;
    // std::cout << "q := " << q << std::endl;	

    g = p.gcd(q);
    
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(q.convertToExpressionTree().toMapleString());

    if (mapleTest->testProcReturn("gcd", inputs, g.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP gcd() test: FAILED4" << std::endl;
        std::cerr << "Expected " << g << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1); 
    }
    
    //Test different nvar and different variables names

    std::vector<Symbol> vars = Symbol::randomElements(nvar+2);
    Symbol v1 = vars[nvar], v2 = vars[nvar+1];
    vars.erase(vars.begin()+nvar);
    vars.erase(vars.begin()+nvar+1);
    // p.setRingVariables(vars);
    vars.push_back(v1);
    vars.push_back(v2);
    q.setRingVariables(vars);



    p = p.primitivePart();
    q = q.primitivePart();
    
    // std::cout << "Test5: " << std::endl;
    // std::cout << "p := " << p << std::endl;
    // std::cout << "q := " << q << std::endl;
    
    g = p.gcd(q);
    inputs.clear();

    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(q.convertToExpressionTree().toMapleString());

    if (mapleTest->testProcReturn("gcd", inputs, g.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP gcd() test: FAILED4" << std::endl;
        std::cerr << "Expected " << g << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1); 
    }
    

    std::cerr << "SMQP gcd() test: PASSED" << std::endl;

}

void testContent() {

    SparseMultivariateRationalPolynomial p;
    
    p.zero();
    RationalNumber content = p.content();
    if (content != 0) {
        std::cerr << "SMQP content() test FAILED" << std::endl;
        std::cerr << "Content should be 0 since p is 0 but got: " << content << std::endl;
        exit(1);
    }

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    content = p.content();

    if (p.isZero() && content != 0) {
        std::cerr << "SMQP content() test FAILED" << std::endl;
        std::cerr << "Content should be 0 but got: " << content << std::endl;
        exit(1);
    }

    if (content < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SMQP content() test FAILED" << std::endl;
            std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
        }
        //change content to be positive here as maple's content is always positive.
        content *= -1;
    }



    //force p to be all integers
    // int max = pow(2,coefBound);
    // ++max;
    // for (int i = 2; i < max; ++i) {
    //     p *= RationalNumber(i);
    // }
    // content = p.content();

    std::vector<Symbol> vars = p.ringVariables();
    std::string varList;
    if (vars.size() > 0) {
        varList += "[" + vars[0].toString();
        vars.erase(vars.begin());
        for (Symbol v : vars) {
            varList += ",";
            varList += v.toString();
        }
        varList += "]";
    } else {
        varList = "x";
    }

    std::vector<std::string> inputs;
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(varList);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if (mapleTest->testProcReturn("content", inputs,content.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP content() test: FAILED" << std::endl;
        std::cerr << "Expected " << content << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1); 
    }


    if (nvar == 0) {
        std::cerr << "SMQP content() test: PASSED" << std::endl;
        return;
    }

    //test content w.r.t some variable
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    vars = p.ringVariables();
    
    std::vector<Symbol> contSyms;
    contSyms.push_back(vars[0]);
    SparseMultivariateRationalPolynomial pCont = p.content(contSyms);

    if (pCont.leadingCoefficient() < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SMQP content() test FAILED" << std::endl;
            std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
        }
        //change content to be positive here as maple's content is always positive.
        pCont *= RationalNumber(-1);
    }

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(vars[0].toString());

    if(mapleTest->testProcReturn("content", inputs, pCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP content(Symbol) test: FAILED1" << std::endl;
        std::cerr << "Expected " << pCont << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }

    //test content w.r.t to some variable with (hopefully) a non-trivial content.
    p = p.head();
    SparseMultivariateRationalPolynomial q;
    for (int i = 0; i < numTerms; ++i) {
        q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        p += q.head();
    }

    pCont = p.content(contSyms);
    if (pCont.leadingCoefficient() < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SMQP content() test FAILED" << std::endl;
            std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
        }
        //change content to be positive here as maple's content is always positive.
        pCont *= RationalNumber(-1);
    }

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(vars[0].toString());

    if(mapleTest->testProcReturn("content", inputs, pCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP content(Symbol) test: FAILED2" << std::endl;
        std::cerr << "Expected " << pCont << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }

    if (nvar > 1) {
        contSyms.clear();
        contSyms.push_back(vars[1]);
        pCont = p.content(contSyms);
        if (pCont.leadingCoefficient() < 0) {
            if (p.leadingCoefficient() > 0) {
                std::cerr << "SMQP content() test FAILED" << std::endl;
                std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
            }
            //change content to be positive here as maple's content is always positive.
            pCont *= RationalNumber(-1);
        }

        inputs.clear();
        inputs.push_back(p.convertToExpressionTree().toMapleString());
        inputs.push_back(vars[1].toString());

        if(mapleTest->testProcReturn("content", inputs, pCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP content(Symbol) test: FAILED3" << std::endl;
            std::cerr << "Expected " << pCont << " but got: " << retErr << std::endl;
            std::cerr << "p: " << p << std::endl;
            exit(1);
        }

        contSyms.push_back(vars[0]);
        pCont = p.content(contSyms);
        if (pCont.leadingCoefficient() < 0) {
            if (p.leadingCoefficient() > 0) {
                std::cerr << "SMQP content() test FAILED4" << std::endl;
                std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
            }
            //change content to be positive here as maple's content is always positive.
            pCont *= RationalNumber(-1);
        }

        inputs.clear();
        inputs.push_back(p.convertToExpressionTree().toMapleString());
        std::string varList = "[";
        varList += contSyms[0].toString() + "," + contSyms[1].toString() + "]";
        inputs.push_back(varList);

        if(mapleTest->testProcReturn("content", inputs, pCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP content(Symbol) test: FAILED5" << std::endl;
            std::cerr << "Expected " << pCont << " but got: " << retErr << std::endl;
            std::cerr << "p: " << p << std::endl;
            exit(1);
        }
    }

    std::cerr << "SMQP content() test: PASSED" << std::endl;
}

void testPrimitivePart() {

    SparseMultivariateRationalPolynomial p;
    
    p.zero();
    SparseMultivariateRationalPolynomial primPart = p.primitivePart();
    if (!primPart.isZero()) {
        std::cerr << "SMQP primitivePart() test FAILED" << std::endl;
        std::cerr << "Primitive part should be 0 since p is 0 but got: " << primPart << std::endl;
        exit(1);
    }

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    primPart = p.primitivePart();

    if (p.isZero() && !primPart.isZero()) {
        std::cerr << "SMQP primitivePart() test FAILED" << std::endl;
        std::cerr << "primPart should be 0 but got: " << primPart << std::endl;
        exit(1);
    }

    if (p.leadingCoefficient() < 0 && primPart.leadingCoefficient() < 0) {
        std::cerr << "SMQP primitivePart() test FAILED " << std::endl;
        std::cerr << "P's leadingCoefficient is negative but so is it's primitive part" << std::endl;
        exit(1);
    }


    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= RationalNumber(-1);
    }

    //force p to be all integers
    // int max = pow(2,coefBound);
    // ++max;
    // for (int i = 2; i < max; ++i) {
    //     p *= RationalNumber(i);
    // }
    // content = p.content();

    std::vector<Symbol> vars = p.ringVariables();
    std::string varList;
    if (vars.size() > 0) {
        varList += "[" + vars[0].toString();
        vars.erase(vars.begin());
        for (Symbol v : vars) {
            varList += ",";
            varList += v.toString();
        }
        varList += "]";
    } else {
        varList = "x";
    }

    std::vector<std::string> inputs;
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(varList);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if (mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart() test: FAILED" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1); 
    }


    //test primpart and content together
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    RationalNumber iCont;
    primPart = p.primitivePart(iCont);    
    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= RationalNumber(-1);
    } 
    if (iCont < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SMQP primitivePart(content) test FAILED" << std::endl;
            std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
        }
        //change content to be positive here as maple's content is always positive.
        iCont *= RationalNumber(-1);
    }
    
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(varList);
    if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart(content) test: FAILED1" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }
    if(mapleTest->testProcReturn("content", inputs, iCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart(content) test: FAILED2" << std::endl;
        std::cerr << "Expected " << iCont << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }    







    if (nvar == 0) {
        std::cerr << "SMQP primitivePart() test: PASSED" << std::endl;
        return;
    }

    //test w.r.t some variable
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    vars = p.ringVariables();
    
    std::vector<Symbol> contSyms;
    contSyms.push_back(vars[0]);
    primPart = p.primitivePart(contSyms);

    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= RationalNumber(-1);
    }

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(vars[0].toString());

    if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart(Symbol) test: FAILED1" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }

    //test content w.r.t to some variable with (hopefully) a non-trivial content.
    p = p.head();
    SparseMultivariateRationalPolynomial q;
    for (int i = 0; i < numTerms; ++i) {
        q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        p += q.head();
    }

    primPart = p.primitivePart(contSyms);

    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= RationalNumber(-1);
    }

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(vars[0].toString());

    if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart(Symbol) test: FAILED2" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }

    //test primpart and content together

    SparseMultivariateRationalPolynomial pCont;
    primPart = p.primitivePart(contSyms, pCont);
    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= RationalNumber(-1);
    }
    if (pCont.leadingCoefficient() < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SMQP content() test FAILED" << std::endl;
            std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
        }
        //change content to be positive here as maple's content is always positive.
        pCont *= RationalNumber(-1);
    }

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(vars[0].toString());
    if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart(Symbol, content) test: FAILED1" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }
    if(mapleTest->testProcReturn("content", inputs, pCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePart(Symbol, content) test: FAILED2" << std::endl;
        std::cerr << "Expected " << pCont << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }    

    if (nvar > 1) {
        contSyms.clear();
        contSyms.push_back(vars[1]);
        primPart = p.primitivePart(contSyms);
    
        if (p.leadingCoefficient() < 0) {
            //maple returns negative prim part if p is.
            primPart *= RationalNumber(-1);
        }

        inputs.clear();
        inputs.push_back(p.convertToExpressionTree().toMapleString());
        inputs.push_back(vars[1].toString());

        if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP primitivePart(Symbol) test: FAILED3" << std::endl;
            std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
            std::cerr << "p: " << p << std::endl;
            exit(1);
        }

        contSyms.push_back(vars[0]);
        primPart = p.primitivePart(contSyms);

        if (p.leadingCoefficient() < 0) {
            //maple returns negative prim part if p is.
            primPart *= RationalNumber(-1);
        }
        
        inputs.clear();
        inputs.push_back(p.convertToExpressionTree().toMapleString());
        std::string varList = "[";
        varList += contSyms[0].toString() + "," + contSyms[1].toString() + "]";
        inputs.push_back(varList);

        if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP primitivePart(Symbol) test: FAILED5" << std::endl;
            std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
            std::cerr << "p: " << p << std::endl;
            exit(1);
        }
    }

    std::cerr << "SMQP primitivePart() test: PASSED" << std::endl;
}

void testPrimitivePartSMZP() {
    SparseMultivariateRationalPolynomial p;
    SparseMultivariateIntegerPolynomial primPart;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    primPart = p.primitivePartSMZP();

    if (p.isZero() && !primPart.isZero()) {
        std::cerr << "SMQP primitivePartSMZP() test FAILED" << std::endl;
        std::cerr << "primPart should be 0 but got: " << primPart << std::endl;
        exit(1);
    }

    if (p.leadingCoefficient() < 0 && primPart.leadingCoefficient() < 0) {
        std::cerr << "SMQP primitivePartSMZP() test FAILED " << std::endl;
        std::cerr << "P's leadingCoefficient is negative but so is it's primitive part" << std::endl;
        exit(1);
    }

    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= Integer(-1);
    }

    std::vector<Symbol> vars = p.ringVariables();
    std::string varList;
    if (vars.size() > 0) {
        varList += "[" + vars[0].toString();
        vars.erase(vars.begin());
        for (Symbol v : vars) {
            varList += ",";
            varList += v.toString();
        }
        varList += "]";
    } else {
        varList = "x";
    }

    std::vector<std::string> inputs;
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(varList);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if (mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePartSMZP() test: FAILED" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1); 
    }


    //test primpart and content together
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    RationalNumber iCont;
    primPart = p.primitivePartSMZP(iCont);    
    if (p.leadingCoefficient() < 0) {
        //maple returns negative prim part if p is.
        primPart *= RationalNumber(-1);
    } 
    if (iCont < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SMQP primitivePartSMZP(content) test FAILED" << std::endl;
            std::cerr << "Content was negative but p's leading coefficient was not." << std::endl;
        }
        //change content to be positive here as maple's content is always positive.
        iCont *= RationalNumber(-1);
    }

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString());
    inputs.push_back(varList);
    if(mapleTest->testProcReturn("primpart", inputs, primPart.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePartSMZP(content) test: FAILED1" << std::endl;
        std::cerr << "Expected " << primPart << " but got: " << retErr << std::endl;
        std::cerr << "p: " << p << std::endl;
        exit(1);
    }
    if(mapleTest->testProcReturn("content", inputs, iCont.convertToExpressionTree().toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP primitivePartSMZP(content) test: FAILED2" << std::endl;
        std::cerr << "Expected " << iCont << " but got: " << retErr << std::endl;
    std::cerr << "p: " << p << std::endl;
        exit(1);
    }    

    std::cerr << "SMQP primitivePartSMZP(content) test: PASSED" << std::endl;

}

void testSquareFree() {
    
    //Test a specific example
    SparseMultivariateRationalPolynomial q(2);
    int deg[2];
    deg[0] = 3;
    deg[1] = 3;
    q.setCoefficient(2, deg, 16);
    deg[1] = 2;
    q.setCoefficient(2, deg, 16);
    deg[0] = 2;
    deg[1] = 3;
    q.setCoefficient(2, deg, -12);
    deg[1] = 2;
    q.setCoefficient(2, deg, -12);
    deg[0] = 0;
    deg[1] = 3;
    q.setCoefficient(2, deg, 1);
    deg[1] = 2;
    q.setCoefficient(2, deg, 1);

    q *= RationalNumber(4);

    Factors<SparseMultivariateRationalPolynomial> qsf = q.squareFree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();

    ALGEB qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
    std::string opStr = "op:";
    ALGEB opFunc = EvalMapleStatement(kv, const_cast<char*>(opStr.c_str()));
    std::string oneStr = "1:";
    std::string twoStr = "2:";
    ALGEB one = EvalMapleStatement(kv, const_cast<char*>(oneStr.c_str()));
    ALGEB two = EvalMapleStatement(kv, const_cast<char*>(twoStr.c_str()));
    ALGEB qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
    ALGEB qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);

    ALGEB qAlgeb = mapleTest->expressionTreeToAlgeb(q.convertToExpressionTree());
    std::string sqrStr = "sqrfree:";
    ALGEB sqrFunc = EvalMapleStatement(kv, const_cast<char*>(sqrStr.c_str()));
    ALGEB resALGEB = EvalMapleProc(kv, sqrFunc, 1, qAlgeb);

    ALGEB resContent = EvalMapleProc(kv, opFunc, 2, one, resALGEB);
    ALGEB resFactors = EvalMapleProc(kv, opFunc, 2, two, resALGEB);

    bool contentEqual = mapleTest->testEquality(resContent, qsfContent);
    bool factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");

    if (!contentEqual || !factorsEqual) {
        std::cerr << "SMQP squareFree() test: FAILED" << std::endl;
        std::cerr << "Expected factors: " << qsf << " but got " << mapleTest->algebToString(kv, resALGEB) << std::endl;
        exit(1);
    }

    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q += RationalNumber(1); //ensure constant term
    qsf = q.squareFree();
    qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
    qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
    qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);

    qAlgeb = mapleTest->expressionTreeToAlgeb(q.convertToExpressionTree());
    resALGEB = EvalMapleProc(kv, sqrFunc, 1, qAlgeb);

    resContent = EvalMapleProc(kv, opFunc, 2, one, resALGEB);
    resFactors = EvalMapleProc(kv, opFunc, 2, two, resALGEB);

    contentEqual = mapleTest->testEquality(resContent, qsfContent);
    factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");


    if (!contentEqual || !factorsEqual) {
        std::cerr << "SMQP squareFree() test: FAILED" << std::endl;
        std::cerr << "Expected factors: " << qsf << " but got " << mapleTest->algebToString(kv, resALGEB) << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1);
    }

    //non-trivial example;
    SparseMultivariateRationalPolynomial p;
    q.zero();
    for (int i = 0; i < numTerms; ++i) {
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        q += p[0];
    }

    qsf = q.squareFree();
    qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
    qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
    qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);

    qAlgeb = mapleTest->expressionTreeToAlgeb(q.convertToExpressionTree());
    resALGEB = EvalMapleProc(kv, sqrFunc, 1, qAlgeb);

    resContent = EvalMapleProc(kv, opFunc, 2, one, resALGEB);
    resFactors = EvalMapleProc(kv, opFunc, 2, two, resALGEB);

    contentEqual = mapleTest->testEquality(resContent, qsfContent);
    factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");

    if (!contentEqual || !factorsEqual) {
        std::cerr << "SMQP squareFree() test: FAILED" << std::endl;
        std::cerr << "Expected factors: " << qsf << " but got " << mapleTest->algebToString(kv, resALGEB) << std::endl;
        std::cerr << "q: " << q << std::endl;
        exit(1);
    }

    if (nvar >= 2){
        
        std::vector<Symbol> ringVars = q.ringVariables();
        std::vector<Symbol> vars;
        vars.push_back(ringVars[0]);

        qsf = q.squareFree(vars);
        qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
        qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
        qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);

        qAlgeb = mapleTest->expressionTreeToAlgeb(q.convertToExpressionTree());

        std::string varString = "[";
        varString += ringVars[0].toString() + "]:";
        ALGEB varsALGEB = EvalMapleStatement(kv, const_cast<char*>(varString.c_str()));

        resALGEB = EvalMapleProc(kv, sqrFunc, 2, qAlgeb, varsALGEB);

        resContent = EvalMapleProc(kv, opFunc, 2, one, resALGEB);
        resFactors = EvalMapleProc(kv, opFunc, 2, two, resALGEB);

        contentEqual = mapleTest->testEquality(resContent, qsfContent);
        factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");

        if (!contentEqual || !factorsEqual) {
            //try again after multiply through by negative 1
            //due to sign of contents.
            qsf[0].first *= RationalNumber(-1);
            qsf.multiplyRingElement(RationalNumber(-1));
            qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
            qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
            qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);
            contentEqual = mapleTest->testEquality(resContent, qsfContent);
            factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");
            if (!contentEqual || !factorsEqual) {
                std::cerr << "SMQP squareFree(vars) test: FAILED1" << std::endl;
                std::cerr << "Expected factors: " << qsf << " but got " << mapleTest->algebToString(kv, resALGEB) << std::endl;
                std::cerr << "q: " << q << std::endl;
                exit(1);
            }
        }

        vars.push_back(ringVars[1]);
        qsf = q.squareFree(vars);
        qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
        qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
        qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);

        qAlgeb = mapleTest->expressionTreeToAlgeb(q.convertToExpressionTree());

        varString = "[";
        varString += ringVars[0].toString() + "," + ringVars[1].toString() + "]:";
        varsALGEB = EvalMapleStatement(kv, const_cast<char*>(varString.c_str()));

        resALGEB = EvalMapleProc(kv, sqrFunc, 2, qAlgeb, varsALGEB);

        resContent = EvalMapleProc(kv, opFunc, 2, one, resALGEB);
        resFactors = EvalMapleProc(kv, opFunc, 2, two, resALGEB);

        contentEqual = mapleTest->testEquality(resContent, qsfContent);
        factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");

        if (!contentEqual || !factorsEqual) {
            //try again after multiply through by negative 1
            //due to sign of contents.
            qsf[0].first *= RationalNumber(-1);
            qsf.multiplyRingElement(RationalNumber(-1));
            qsfALGEB = mapleTest->expressionTreeToAlgeb(qsf.convertToExpressionTree());
            qsfContent = EvalMapleProc(kv, opFunc, 2, one, qsfALGEB);
            qsfFactors = EvalMapleProc(kv, opFunc, 2, two, qsfALGEB);
            contentEqual = mapleTest->testEquality(resContent, qsfContent);
            factorsEqual = mapleTest->testEquality(resFactors, qsfFactors, "as_set");
            if (!contentEqual || !factorsEqual) {
                std::cerr << "SMQP squareFree(vars) test: FAILED2" << std::endl;
                std::cerr << "Expected factors: " << qsf << " but got " << mapleTest->algebToString(kv, resALGEB) << std::endl;
                std::cerr << "q: " << q << std::endl;
                exit(1);
            }
        }
    }

    std::cerr << "SMQP squareFree() test: PASSED" << std::endl;
}

void testSquareFreePart() {

    //Test a specific example
    SparseMultivariateRationalPolynomial q(2);
    int deg[2];
    deg[0] = 3;
    deg[1] = 3;
    q.setCoefficient(2, deg, 16);
    deg[1] = 2;
    q.setCoefficient(2, deg, 16);
    deg[0] = 2;
    deg[1] = 3;
    q.setCoefficient(2, deg, -12);
    deg[1] = 2;
    q.setCoefficient(2, deg, -12);
    deg[0] = 0;
    deg[1] = 3;
    q.setCoefficient(2, deg, 1);
    deg[1] = 2;
    q.setCoefficient(2, deg, 1);

    q *= RationalNumber(4);


    Factors<SparseMultivariateRationalPolynomial> qsf = q.squareFree();
    SparseMultivariateRationalPolynomial qsfp = q.squareFreePart();
    SparseMultivariateRationalPolynomial fact;
    fact.one();
    for (int i = 0; i < qsf.size(); ++i) {
        fact *= qsf[i].first;
    }
    if (qsfp != fact) {
        std::cerr << "SMQP squareFreePart() test: FAILED1" << std::endl;
        std::cerr << "Got squareFreePart: " << qsfp << std::endl;
        std::cerr << "But got square free factorization: " << qsf << std::endl;
        exit(1);
    }


    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    qsf = q.squareFree();
    qsfp = q.squareFreePart();
    fact.one();
    for (int i = 0; i < qsf.size(); ++i) {
        fact *= qsf[i].first;
    }
    if (qsfp != fact) {
        std::cerr << "SMQP squareFreePart() test: FAILED2" << std::endl;
        std::cerr << "q: " << q << std::endl;
        std::cerr << "Got squareFreePart: " << qsfp << std::endl;
        std::cerr << "But got square free factorization: " << qsf << std::endl;
        exit(1);
    }


    //non-trivial example;
    SparseMultivariateRationalPolynomial p;
    q.zero();
    for (int i = 0; i < numTerms; ++i) {
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        q += p[0];
    }

    qsf = q.squareFree();
    qsfp = q.squareFreePart();
    fact.one();
    for (int i = 0; i < qsf.size(); ++i) {
        fact *= qsf[i].first;
    }
    if (qsfp != fact) {
        std::cerr << "SMQP squareFreePart() test: FAILED3" << std::endl;
        std::cerr << "Got squareFreePart: " << qsfp << std::endl;
        std::cerr << "But got square free factorization: " << qsf << std::endl;
        exit(1);
    }

    if (nvar >= 2){
        
        std::vector<Symbol> ringVars = q.ringVariables();
        std::vector<Symbol> vars;
        vars.push_back(ringVars[1]);

        qsf = q.squareFree(vars);
        qsfp = q.squareFreePart(vars);
        fact.one();
        for (int i = 0; i < qsf.size(); ++i) {
            fact *= qsf[i].first;
        }
        if (qsfp != fact) {
            std::cerr << "SMQP squareFreePart(vars) test: FAILED1" << std::endl;
            std::cerr << "Got squareFreePart: " << qsfp << std::endl;
            std::cerr << "But got square free factorization: " << qsf << std::endl;
            exit(1);
        }

        vars.push_back(ringVars[0]);
        qsf = q.squareFree(vars);
        qsfp = q.squareFreePart(vars);
        fact.one();
        for (int i = 0; i < qsf.size(); ++i) {
            fact *= qsf[i].first;
        }
        if (qsfp != fact) {
            std::cerr << "SMQP squareFreePart(vars) test: FAILED2" << std::endl;
            std::cerr << "Got squareFreePart: " << qsfp << std::endl;
            std::cerr << "But got square free factorization: " << qsf << std::endl;
            exit(1);
        }


    }

    std::cerr << "SMQP squareFree() test: PASSED" << std::endl;

}

void testNumberOfRingVariables() {
    SparseMultivariateRationalPolynomial p(10);
    SparseMultivariateRationalPolynomial q(24);

    if (p.numberOfRingVariables() == 10 && q.numberOfRingVariables() == 24) {
        std::cerr << "SMQP numberOfRingVariables() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SMQP numberOfRingVariables() test: FAILED" << std::endl;
    exit(1);
}

void testNumberOfVariables() {
    SparseMultivariateRationalPolynomial p(3);
    SparseMultivariateRationalPolynomial x("x");
    SparseMultivariateRationalPolynomial y("y");
    SparseMultivariateRationalPolynomial z("z");

    int zero = p.numberOfVariables();
    p += x;
    int one = p.numberOfVariables();
    p += y;
    int two = p.numberOfVariables();
    p += z;
    int three = p.numberOfVariables();

    if (zero != 0 || one != 1 || two != 2 || three != 3) {
        std::cerr << "SMQP numberOfVariables() test: FAILED" << std::endl;
        std::cerr << "zero: " << zero << " one: " << one << " two: " << two << " three: " << three << std::endl;  
        exit(1);
    }

    std::cerr << "SMQP numberOfVariables() test: PASSED" << std::endl;
}

void testNumberOfTerms() {

    if (nvar > 0) {
        SparseMultivariateRationalPolynomial p;
        p.randomPolynomial(nvar, 100, coefBound, sparsity, includeNeg);
        SparseMultivariateRationalPolynomial q;
        q.randomPolynomial(nvar, 222, coefBound, sparsity, includeNeg);
        SparseMultivariateRationalPolynomial r;
        r.one();

        if (p.numberOfTerms() == 100 && q.numberOfTerms() == 222 && r.numberOfTerms() == 1) {
            std::cerr << "SMQP numberOfTerms() test: PASSED" << std::endl;
            return;
        }
    } else {
        SparseMultivariateRationalPolynomial p;
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        SparseMultivariateRationalPolynomial q;
        q.one();

        if(p.numberOfTerms() == 1 && q.numberOfTerms() == 1) {
            std::cerr << "SMQP numberOfTerms() test: PASSED" << std::endl;
            return;
        }
    }

    std::cerr << "SMQP numberOfTerms() test: FAILED" << std::endl;
    exit(1);
}

void testDegree() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    std::vector<Symbol> vars = p.ringVariables();
    
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    ExpressionTree pTree = p.convertToExpressionTree();
    std::string retErr;
    for (int i = 0; i < vars.size(); ++i) {
        std::vector<std::string> inputs;
        inputs.push_back(pTree.toMapleString());
        inputs.push_back(vars[i].toString());
        int degree = p.degree(vars[i]).get_si();
        ExpressionTree expected = ExpressionTree(new ExprTreeNode(degree));
        if (mapleTest->testProcReturn("degree", inputs, expected.toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP degree() test: FAILED" << std::endl;
            std::cerr << "Got " << retErr << " but expected " << degree << std::endl;
            exit(1);
        }
    }

    std::cerr << "SMQP degree() test: PASSED" << std::endl;
}

void testLeadingCoefficient() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    RationalNumber leadingCoef = p.leadingCoefficient();
    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree expected = ExpressionTree(new ExprTreeNode(leadingCoef.get_mpq()));
    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testProcReturn("lcoeff", inputs, expected.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP leadingCoefficient() test: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << leadingCoef << std::endl;
        std::cerr << p << std::endl;
        exit(1); 
    }

    std::cerr << "SMQP leadingCoefficient() test: PASSED" << std::endl;
}

void testCoefficient() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    int* degs = (int*) calloc(nvar, sizeof(int));
    if (nvar > 2) {
        degs[0] = 2;
    }
    if (nvar > 4) {
        degs[1] = 3;
        degs[3] = 1;
    }
    if (nvar > 6) {
        degs[4] = 2;
    }


    RationalNumber got1 = p.coefficient(nvar, degs);

    RationalNumber expected(27, 34);
    if (got1 == expected){
        expected = RationalNumber(29, 34);
    }
    p.setCoefficient(nvar, degs, expected);

    RationalNumber got = p.coefficient(nvar, degs);

    if ( (got1 != 0 && got1 == got) || expected != got) {
        std::cerr << "SMQP coefficient() and setCoefficient() test: FAILED" << std::endl;
        std::cerr << "Got " << got << " but expected " << expected << std::endl;
        exit(1);
    }

    std::cerr << "SMQP coefficient() and setCoefficient() test: PASSED" << std::endl;
}

void testSetVariableNames() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    std::vector<Symbol> newVars; 
    for(int i = 0; i < nvar; ++i) {
        std::string var = "i_" + std::to_string(i) ;
        newVars.push_back(Symbol(var));
    }

    p.setRingVariables(newVars);
    std::vector<Symbol> gotVars = p.ringVariables();
    for (int i = 0; i < nvar; ++i) {
        if (newVars[i] != gotVars[i]) {
            std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED1" << std::endl;
            std::cerr << "Got " << gotVars[i] << " but expected " << newVars[i] << std::endl;
            exit(1);
        }
    }

    std::reverse(newVars.begin(), newVars.end());

    ExpressionTree pBefore = p.convertToExpressionTree();
    SparseMultivariateRationalPolynomial smqpBefore = p;

    p.setRingVariables(newVars);
    gotVars = p.ringVariables();

    pBefore -= p.convertToExpressionTree();

    for (int i = 0; i < nvar; ++i) {
        if (newVars[i] != gotVars[i]) {
            std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED2" << std::endl;
            std::cerr << "Got " << gotVars[i] << " but expected " << newVars[i] << std::endl;
            exit(1);
        }
    }

    //Note here that maple checks equality without regard to variable ordering.
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testIfZero(pBefore, &retErr) == 0) {
        std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED" << std::endl;
        std::cerr << "Expected polynomials to be 'equal' after changing variable ordering" << std::endl;
        std::cerr << "But their difference is: " << retErr << std::endl;
        exit(1); 
    }

    //now test expanding the ring
    newVars.push_back(Symbol("foo"));
    
    pBefore = p.convertToExpressionTree();

    p.setRingVariables(newVars);
    gotVars = p.ringVariables();
    for (int i = 0; i < newVars.size(); ++i) {
        if (newVars[i] != gotVars[i]) {
            std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED3" << std::endl;
            std::cerr << "Got " << gotVars[i] << " but expected " << newVars[i] << std::endl;
            exit(1);
        }
    }

    pBefore -= p.convertToExpressionTree();
    if (mapleTest->testIfZero(pBefore, &retErr) == 0) {
        std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED" << std::endl;
        std::cerr << "Expected polynomials to be 'equal' after expanding polynomial ring" << std::endl;
        std::cerr << "But their difference is: " << retErr << std::endl;
        exit(1); 
    }


    //now test simultaneous expanding and re-ordering;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    newVars = p.ringVariables();

    std::reverse(newVars.begin(), newVars.end());
    if (nvar > 1) {
        newVars.insert(newVars.begin()+1, Symbol("foo"));
    } else {
        newVars.insert(newVars.begin(), Symbol("foo"));
    }
    
    pBefore = p.convertToExpressionTree();

    p.setRingVariables(newVars);
    gotVars = p.ringVariables();

    for (int i = 0; i < newVars.size(); ++i) {
        if (newVars[i] != gotVars[i]) {
            std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED4" << std::endl;
            std::cerr << "Got " << gotVars[i] << " but expected " << newVars[i] << std::endl;
            exit(1);
        }
    }

    pBefore -= p.convertToExpressionTree();
    if (mapleTest->testIfZero(pBefore, &retErr) == 0) {
        std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED" << std::endl;
        std::cerr << "Expected polynomials to be 'equal' after expanding and reordered polynomial ring" << std::endl;
        std::cerr << "But their difference is: " << retErr << std::endl;
        exit(1); 
    }

    if (nvar > 0) {
        //Test shrinking the number of variables;
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        newVars = p.ringVariables();
        int rmIdx = 1;
        newVars.insert(newVars.begin()+rmIdx, Symbol("foo"));
        p.setRingVariables(newVars);

        pBefore = p.convertToExpressionTree();

        newVars.erase(newVars.begin()+rmIdx);
        p.setRingVariables(newVars);
        gotVars = p.ringVariables();

        for (int i = 0; i < newVars.size(); ++i) {
            if (newVars[i] != gotVars[i]) {
                std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED5" << std::endl;
                std::cerr << "Got " << gotVars[i] << " but expected " << newVars[i] << std::endl;
                exit(1);
            }
        }

        pBefore -= p.convertToExpressionTree();
        if (mapleTest->testIfZero(pBefore, &retErr) == 0) {
            std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED" << std::endl;
            std::cerr << "Expected polynomials to be 'equal' after shrinking polynomial ring" << std::endl;
            std::cerr << "But their difference is: " << retErr << std::endl;
            exit(1); 
        }

        //Test shrinking and reordering the number of variables;
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        newVars = p.ringVariables();
        rmIdx = 1;
        newVars.insert(newVars.begin()+rmIdx, Symbol("foo"));
        p.setRingVariables(newVars);

        pBefore = p.convertToExpressionTree();
        newVars.erase(newVars.begin() + rmIdx);
        std::reverse(newVars.begin(), newVars.end());
        p.setRingVariables(newVars);
        gotVars = p.ringVariables();

        for (int i = 0; i < newVars.size(); ++i) {
            if (newVars[i] != gotVars[i]) {
                std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED6" << std::endl;
                std::cerr << "Got " << gotVars[i] << " but expected " << newVars[i] << std::endl;
                exit(1);
            }
        }

        pBefore -= p.convertToExpressionTree();
        if (mapleTest->testIfZero(pBefore, &retErr) == 0) {
            std::cerr << "SMQP setRingVariables() and ringVariables() test: FAILED" << std::endl;
            std::cerr << "Expected polynomials to be 'equal' after shrinking and reordering polynomial ring" << std::endl;
            std::cerr << "But their difference is: " << retErr << std::endl;
            exit(1); 
        }
    }


    //Test non-zero variables:

    SparseMultivariateRationalPolynomial q(3);

    SparseMultivariateRationalPolynomial x("x");
    SparseMultivariateRationalPolynomial y("y");
    SparseMultivariateRationalPolynomial z("z");

    std::vector<Symbol> zero = q.variables();
    q += x;
    std::vector<Symbol> one = q.variables();
    q += y;
    std::vector<Symbol> two = q.variables();
    q += z;
    std::vector<Symbol> three = q.variables();

    if (zero.size() != 0) {
        std::cerr << "SMQP ringVariables() test: FAILED" << std::endl;
        std::cerr << "zero vars failed" << std::endl;
        exit(1);
    }
    if (one.size() != 1 || one[0] != "x") {
        std::cerr << "SMQP ringVariables() test: FAILED" << std::endl;
        std::cerr << "one vars failed" << std::endl;
        exit(1);
    }
    if (two.size() != 2 || std::find(two.begin(), two.end(), "x") == two.end()
                        || std::find(two.begin(), two.end(), "y") == two.end()) {
        std::cerr << "SMQP ringVariables() test: FAILED" << std::endl;
        std::cerr << "two vars failed" << std::endl;
        for (auto sym : two) {
            std::cerr << "sym: " << sym << std::endl;
        }
        exit(1);
    }
    if (three.size() != 3 || std::find(three.begin(), three.end(), "x") == three.end() 
                          || std::find(three.begin(), three.end(), "y") == three.end()
                          || std::find(three.begin(), three.end(), "z") == three.end()) {
        std::cerr << "SMQP ringVariables() test: FAILED" << std::endl;
        std::cerr << "three vars failed" << std::endl;
        for (auto sym : three) {
            std::cerr << "sym: " << sym << std::endl;
        }
        exit(1);
    } 


    std::cerr << "SMQP setRingVariables() and ringVariables() test: PASSED" << std::endl;

}

void testInitial() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial init = p.initial();
  
    if (p.isConstant() || p.isZero()) {
        if (p != init) {
            std::cerr << "SMQP initial() test: FAILED" << std::endl;
            std::cerr << "p: " << p << std::endl;
            std::cerr << "initial: " << init << std::endl;
            exit(1);
        }
        std::cerr << "SMQP initial() test: PASSED" << std::endl;
        return;
    }

    ExpressionTree pTree = p.convertToExpressionTree();

    ExpressionTree resTree = init.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    MKernelVector kv = mapleTest->getMKernelVector();
    mapleTest->restartMapleKernel();


    char polyRingCstr[] = "RegularChains:-PolynomialRing:";

    std::vector<Symbol> vars = p.ringVariables();
    std::string varList = "[" + vars[0].toString();
    vars.erase(vars.begin());
    for (Symbol v : vars) {
        varList += ",";
        varList += v.toString();
    }
    varList += "]:";

    std::string polyString = pTree.toMapleString() + ":";
    ALGEB poly = EvalMapleStatement(kv, const_cast<char*>(polyString.c_str()));

    ALGEB varListALGEB = EvalMapleStatement(kv, const_cast<char*>(varList.c_str()));

    ALGEB polyRingCmd = EvalMapleStatement(kv, polyRingCstr);
    ALGEB polyRing = EvalMapleProc(kv, polyRingCmd, 1, varListALGEB);

    char initialCstr[] = "RegularChains:-Initial:";
    ALGEB initialCmd = EvalMapleStatement(kv, initialCstr);

    ALGEB initALGEB = EvalMapleProc(kv, initialCmd, 2, poly, polyRing);
    bool equal = mapleTest->testEquality(initALGEB, resTree.toMapleString());

    if (!equal) {
        std::cerr << "SMQP initial() test: FAILED" << std::endl;
        std::cerr << "Expected: " << init << " but got: " << mapleTest->algebToString(kv, initALGEB) << std::endl;
        exit(1);
    }

    std::cerr << "SMQP initial() test: PASSED" << std::endl;
}

void testHead() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial head = p.head();
  
    if (p.isConstant() || p.isZero()) {
        if (p != head) {
            std::cerr << "SMQP head() test: FAILED" << std::endl;
            std::cerr << "p: " << p << std::endl;
            std::cerr << "head: " << head << std::endl;
            exit(1);
        }
        std::cerr << "SMQP head() test: PASSED" << std::endl;
        return;
    }

    ExpressionTree pTree = p.convertToExpressionTree();

    ExpressionTree resTree = head.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    MKernelVector kv = mapleTest->getMKernelVector();
    mapleTest->restartMapleKernel();

    char polyRingCstr[] = "RegularChains:-PolynomialRing:";

    std::vector<Symbol> vars = p.ringVariables();
    std::string varList = "[" + vars[0].toString();
    vars.erase(vars.begin());
    for (Symbol v : vars) {
        varList += ",";
        varList += v.toString();
    }
    varList += "]:";

    std::string polyString = pTree.toMapleString() + ":";
    ALGEB poly = EvalMapleStatement(kv, const_cast<char*>(polyString.c_str()));

    ALGEB varListALGEB = EvalMapleStatement(kv, const_cast<char*>(varList.c_str()));

    ALGEB polyRingCmd = EvalMapleStatement(kv, polyRingCstr);
    ALGEB polyRing = EvalMapleProc(kv, polyRingCmd, 1, varListALGEB);

    char initialCstr[] = "RegularChains:-Initial:";
    ALGEB initialCmd = EvalMapleStatement(kv, initialCstr);
    ALGEB initALGEB = EvalMapleProc(kv, initialCmd, 2, poly, polyRing);

    char rankCStr[] = "RegularChains:-Rank:";
    ALGEB rankCmd = EvalMapleStatement(kv, rankCStr);
    ALGEB rankALGEB = EvalMapleProc(kv, rankCmd, 2, poly, polyRing);

    std::string initAStr = mapleTest->algebToString(kv, initALGEB);
    std::string rankAStr = mapleTest->algebToString(kv, rankALGEB);
    std::string headExpr = "(" + initAStr + ")*(" + rankAStr + "):";
    ALGEB headALGEB = EvalMapleStatement(kv, const_cast<char*>(headExpr.c_str()));

    bool equal = mapleTest->testEquality(headALGEB, resTree.toMapleString());

    if (!equal) {
        std::cerr << "SMQP head() test: FAILED" << std::endl;
        std::cerr << "Expected: " << head << " but got: " << mapleTest->algebToString(kv, headALGEB) << std::endl;
        exit(1);
    }

    std::cerr << "SMQP head() test: PASSED" << std::endl;

}

void testTail() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial tail = p.tail();
  
    if (p.isConstant() || p.isZero()) {
        if (!tail.isZero()) {
            std::cerr << "SMQP tail() test: FAILED" << std::endl;
            std::cerr << "p: " << p << std::endl;
            std::cerr << "tail: " << tail << std::endl;
            exit(1);
        }
        std::cerr << "SMQP tail() test: PASSED" << std::endl;
        return;
    }

    ExpressionTree pTree = p.convertToExpressionTree();

    ExpressionTree resTree = tail.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    MKernelVector kv = mapleTest->getMKernelVector();
    mapleTest->restartMapleKernel();

    char polyRingCstr[] = "RegularChains:-PolynomialRing:";

    std::vector<Symbol> vars = p.ringVariables();
    std::string varList = "[" + vars[0].toString();
    vars.erase(vars.begin());
    for (Symbol v : vars) {
        varList += ",";
        varList += v.toString();
    }
    varList += "]:";

    std::string polyString = pTree.toMapleString() + ":";
    ALGEB poly = EvalMapleStatement(kv, const_cast<char*>(polyString.c_str()));
    ALGEB varListALGEB = EvalMapleStatement(kv, const_cast<char*>(varList.c_str()));
    ALGEB polyRingCmd = EvalMapleStatement(kv, polyRingCstr);
    ALGEB polyRing = EvalMapleProc(kv, polyRingCmd, 1, varListALGEB);

    char tailCstr[] = "RegularChains:-Tail:";
    ALGEB tailCmd = EvalMapleStatement(kv, tailCstr);

    ALGEB tailALGEB = EvalMapleProc(kv, tailCmd, 2, poly, polyRing);
    bool equal = mapleTest->testEquality(tailALGEB, resTree.toMapleString());

    if (!equal) {
        std::cerr << "SMQP tail() test: FAILED" << std::endl;
        std::cerr << "Expected: " << tail << " but got: " << mapleTest->algebToString(kv, tailALGEB) << std::endl;
        exit(1);
    }

    std::cerr << "SMQP tail() test: PASSED" << std::endl;
}


void testIsEqual() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial q = p.deepCopy();
    SparseMultivariateRationalPolynomial r = q;
    r += RationalNumber(2);

    if (!p.isEqual(q) || !q.isEqual(p)) {
        std::cerr << "SMQP isEqual() test: FAILED" << std::endl;
        std::cerr << "Expcted p and q to be equal: " << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;           
        exit(1);
    }
    if (p.isEqual(r) || r.isEqual(p)) {
        std::cerr << "SMQP isEqual() test: FAILED" << std::endl;
        std::cerr << "Expcted p and r to NOT be equal: " << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "r: " << r << std::endl;           
        exit(1);
    }

    std::cerr << "SMQP isEqual() test: PASSED" << std::endl;

    if (!(p == q) || !(q == p)) {
        std::cerr << "SMQP == test: FAILED" << std::endl;
        std::cerr << "Expcted p and q to be equal: " << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "q: " << q << std::endl;           
        exit(1);
    }

    std::cerr << "SMQP == test: PASSED" << std::endl;

    if (!(p != r) || !(r != p)) {
        std::cerr << "SMQP != test: FAILED" << std::endl;
        std::cerr << "Expcted p and r to NOT be equal: " << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "r: " << r << std::endl;           
        exit(1);
    }

    std::cerr << "SMQP != test: PASSED" << std::endl;
}

void testDerivative() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    if (nvar == 0) {
        SparseMultivariateRationalPolynomial dp = p.derivative(Symbol("x"));
        SparseMultivariateRationalPolynomial d3p = p.derivative(Symbol("x"), 3);
        
        SparseMultivariateRationalPolynomial q = p;
        p.differentiate(Symbol("x"));
        q.differentiate(Symbol("x"), 2);
        

        if (!dp.isZero() || !d3p.isZero() || !p.isZero() || !q.isZero()) {
            std::cerr << "SMQP derivative() test: FAILED" << std::endl;
            std::cerr << "p' should be zero but got: " << dp << std::endl;
            exit(1);
        }

        std::cerr << "SMQP derivative() test: PASSED" << std::endl;

        return;
    }

    std::vector<Symbol> syms = p.ringVariables();
    int idx = rand() % syms.size();

    SparseMultivariateRationalPolynomial dp = p.derivative(syms[idx]);
    SparseMultivariateRationalPolynomial d3p = p.derivative(syms[idx], 3);

    SparseMultivariateRationalPolynomial q = p;
    ExpressionTree pTree = p.convertToExpressionTree();

    p.differentiate(syms[idx]);
    q.differentiate(syms[idx], 3);


    if (dp != p || d3p != q) {
        std::cerr << "SMQP differentiate test() FAILED" << std::endl;
        std::cerr << "Should be equal: " << std::endl << p << std::endl << dp << std::endl;
        std::cerr << "Should be equal: " << std::endl << q << std::endl << d3p << std::endl;
        exit(1);
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString());

    ExpressionTree resTree = p.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP derivative() test: FAILED1" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString() + "$" + std::to_string(3));

    resTree = q.convertToExpressionTree();
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP derivative() test: FAILED2" << std::endl;
        std::cerr << "Expected: " << q << " but got: " << retErr << std::endl;
        exit(1);
    }

    //zero
    p.zero();
    pTree = p.convertToExpressionTree();
    p.differentiate(syms[idx]);
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString());
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP derivative() test: FAILED3" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    //constant
    p += RationalNumber(37, 13);
    pTree = p.convertToExpressionTree();
    p.differentiate(syms[idx]);
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString());
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP derivative() test: FAILED4" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    //w.r.t different variable
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = p.convertToExpressionTree();
    p.differentiate(Symbol("y"));
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back("y");
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP derivative() test: FAILED4" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP derivative() test: PASSED" << std::endl;
}

void testIntegrate() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    if (nvar == 0) {
        SparseMultivariateRationalPolynomial dp = p.integral(Symbol("x"));
        SparseMultivariateRationalPolynomial d3p = p.integral(Symbol("x"), 3);

        ExpressionTree pTree = p.convertToExpressionTree();

        SparseMultivariateRationalPolynomial q = p;
        p.integrate(Symbol("x"));
        q.integrate(Symbol("x"), 3);

        if (p != dp || q != d3p) {
            std::cerr << "SMQP integrate() test FAILED" << std::endl;
            std::cerr << "Should be equal: " << std::endl << p << std::endl << dp << std::endl;
            std::cerr << "Should be equal: " << std::endl << q << std::endl << d3p << std::endl;
            exit(1);
        }

        std::vector<std::string> inputs;
        inputs.push_back(pTree.toMapleString());
        inputs.push_back("x");

        ExpressionTree resTree = p.convertToExpressionTree();
        MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
        mapleTest->restartMapleKernel();
        std::string retErr;
        if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP integral() test: FAILED1" << std::endl;
            std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
            exit(1);
        }

        inputs.clear();
        inputs.push_back(pTree.toMapleString());
        inputs.push_back("x$" + std::to_string(3));

        resTree = q.convertToExpressionTree();
        if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP integral() test: FAILED2" << std::endl;
            std::cerr << "Expected: " << q << " but got: " << retErr << std::endl;
            exit(1);
        }
        
        std::cerr << "SMQP integrate() test: PASSED" << std::endl;

        return;
    }

    std::vector<Symbol> syms = p.ringVariables();
    int idx = rand() % syms.size();

    SparseMultivariateRationalPolynomial dp = p.integral(syms[idx]);
    SparseMultivariateRationalPolynomial d3p = p.integral(syms[idx], 3);

    SparseMultivariateRationalPolynomial q = p;
    ExpressionTree pTree = p.convertToExpressionTree();

    p.integrate(syms[idx]);
    q.integrate(syms[idx], 3);

    if (dp != p || d3p != q) {
        std::cerr << "SMQP integrate test() FAILED" << std::endl;
        std::cerr << "Should be equal: " << std::endl << p << std::endl << dp << std::endl;
        std::cerr << "Should be equal: " << std::endl << q << std::endl << d3p << std::endl;
        exit(1);
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString());

    ExpressionTree resTree = p.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP integral() test: FAILED1" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString() + "$" + std::to_string(3));

    resTree = q.convertToExpressionTree();
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP integral() test: FAILED2" << std::endl;
        std::cerr << "Expected: " << q << " but got: " << retErr << std::endl;
        exit(1);
    }

    //zero
    p.zero();
    pTree = p.convertToExpressionTree();
    p.integrate(syms[idx]);
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString());
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP integrate() test: FAILED3" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    //constant
    p += RationalNumber(37, 13);
    pTree = p.convertToExpressionTree();
    p.integrate(syms[idx]);
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms[idx].toString());
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP integrate() test: FAILED4" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    //w.r.t different variable
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = p.convertToExpressionTree();
    p.integrate(Symbol("y"));
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back("y");
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP integrate() test: FAILED" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SMQP integrate() test: PASSED" << std::endl;
}

void testEvaluate() {


    if (nvar == 0) {
        SparseMultivariateRationalPolynomial p;
        p.setCoefficient(0, NULL, RationalNumber(27,44));

        std::vector<Symbol> vars;
        std::vector<RationalNumber> vals;
        SparseMultivariateRationalPolynomial val = p.evaluate(vars, vals);
        
        if (val != RationalNumber(27,44)) {
            std::cerr << "SMQP evaluate() test: FAILED:" << std::endl;
            std::cerr << "Got: " << val << " expected: 27/44";
            exit(1);
        }

        std::cerr << "SMQP evaluate() test: PASSED" << std::endl;
        return;
    }

    //test random
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    //test all vars set
    std::vector<Symbol> vars = p.ringVariables();
    std::vector<RationalNumber> vals;
    SparseMultivariateRationalPolynomial r;
    for (int i = 0; i < vars.size(); ++i) {
        r.randomPolynomial(1, 1, coefBound, sparsity, includeNeg); //get a random rat num
        vals.push_back(r.leadingCoefficient());
    }

    SparseMultivariateRationalPolynomial res = p.evaluate(vars, vals);

    // std::ofstream ofs;
    // ofs.open("mapleTestEval.mpl");
    // ofs << "p:= " << p << ":" << std::endl << std::endl;
    // ofs << "pt:= " << vals[0] << ":" << std::endl << std::endl;
    // ofs << "time1 := time(): v := eval(p, x_1=pt): time2:= time() - time1;" << std::endl << std::endl;
    // ofs.close();

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree resTree = res.convertToExpressionTree();

    std::vector<std::string> inputs; 
    inputs.push_back(pTree.toMapleString());

    std::stringstream varss;
    varss << "[";
    for(int i = 0; i < vars.size()-1; ++i) {
        varss << vars[i] << "=" << vals[i] << ", ";
    }
    varss << vars[vars.size()-1] << "=" << vals[vars.size()-1] << "]";
    inputs.push_back(varss.str());

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testProcReturn("eval", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP evaluate() test: FAILED1" << std::endl;
        std::cerr << "Expected: " << res << " but got: " << retErr << std::endl;
        exit(1);
    }


    //test only some vars set;
    if (nvar > 1) {
        p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

        std::vector<Symbol> subVars;
        vals.clear();
        for (int i = 1; i < nvar; ++i) {
            subVars.push_back(vars[i]);
            vals.push_back(RationalNumber(rand() % coefBound));
        }        
            
        varss.str("");
        varss.clear();
        varss << "[";
        for (int i = 0; i < subVars.size()-1; ++i) {
            varss << subVars[i] << "=" << vals[i] << ", ";    
        }
        varss << subVars[subVars.size()-1] << "=" << vals[subVars.size()-1] << "]";

        res = p.evaluate(subVars, vals);

        pTree = p.convertToExpressionTree();
        resTree = res.convertToExpressionTree();
        inputs.clear();
        inputs.push_back(pTree.toMapleString());
        inputs.push_back(varss.str());

        if (mapleTest->testProcReturn("eval", inputs, resTree.toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP evaluate() test: FAILED2" << std::endl;
            std::cerr << "Expected: " << res << " but got: " << retErr << std::endl;
            exit(1);
        }
    }

    //test zero
    p.zero();

    vars = p.ringVariables();
    vals.clear();
    for (int i = 0; i < vars.size(); ++i) {
        vals.push_back(RationalNumber(rand() % coefBound));
    }
    res = p.evaluate(vars, vals);
    if (!res.isZero()) {
        std::cerr << "SMQP evaluate() test: FAILED3" << std::endl;
        std::cerr << "Expected 0 but got: " << res << std::endl;
        exit(1);
    }

    //test constant
    p.zero();
    mpq_class constTerm(rand() % coefBound);
    p += constTerm;

    vars = p.ringVariables();
    vals.clear();
    for (int i = 0; i < vars.size(); ++i) {
        vals.push_back(RationalNumber(rand() % coefBound));
    }
    res = p.evaluate(vars, vals);
    if (!res.isConstant() || !(res - constTerm).isZero()) {
        std::cerr << "SMQP evaluate() test: FAILED4" << std::endl;
        std::cerr << "Expected a constant: " << p << " but got: " << res << std::endl;
        exit(1);
    }
    
    //test constant with nvar = 0
    p = SparseMultivariateRationalPolynomial(0);
    p.zero();
    constTerm = mpq_class(rand() % coefBound);
    p += constTerm;

    vars = p.ringVariables();
    vals.clear();
    res = p.evaluate(vars, vals);
    if (!res.isConstant() || !(res - constTerm).isZero()) {
        std::cerr << "SMQP evaluate() test: FAILED5" << std::endl;
        std::cerr << "Expected a constant: " << p << " but got: " << res << std::endl;
        exit(1);
    }


    std::cerr << "SMQP evaluate() test: PASSED" << std::endl;
}

void testInterpolate() {

    int maxDeg[3] = {5, 10, 20};
    for (int idx = 0; idx < 3; ++idx) {

        std::vector<std::vector<RationalNumber>> points(maxDeg[idx]);
        std::vector<RationalNumber> vals;
        for (int i = 0; i < maxDeg[idx]; ++i) {
            SparseMultivariateRationalPolynomial t;
            t.randomPolynomial(1, 1, coefBound, 2, includeNeg); //use sparsity = 2 to ensure constant term

            RationalNumber lc = t.leadingCoefficient();
            for (int j = 0; j < i; ++j) {
                if (lc == points[j][0]) {
                    //randomize again and start search at beginning
                    t.randomPolynomial(1, 1, coefBound, 2, includeNeg); //use sparsity = 2 to ensure constant term
                    lc = t.leadingCoefficient();
                    j = -1; //reset j to beginning of loop
                }
            }

            points[i].push_back(t.leadingCoefficient());
        
            t.randomPolynomial(1, 1, coefBound, 2, includeNeg); //use sparsity = 2 to ensure constant term
            vals.push_back(t.leadingCoefficient());
        }

        std::stringstream ss;
        ss << "[";
        for (int i = 0; i < maxDeg[idx]; ++i) {
            ss << "[" << points[i][0] << "," << vals[i] << "]";
            if (i < maxDeg[idx] -1) {
                ss << ", ";
            }
        }
        ss << "]";

        SparseMultivariateRationalPolynomial p = SparseMultivariateRationalPolynomial::interpolate(points, vals);

        std::vector<std::string> inputs;
        inputs.push_back(ss.str());
        inputs.push_back(p.leadingVariable().toString());

        MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
        mapleTest->restartMapleKernel();
        std::string retErr;
        if (mapleTest->testProcReturn("CurveFitting:-PolynomialInterpolation", inputs, p.convertToExpressionTree().toMapleString(), &retErr) == 0) {
            std::cerr << "SMQP interpolate() test: FAILED " << maxDeg[idx] << std::endl;
            std::cerr << "Expected " << p << " but got: " << retErr << std::endl;
            std::cerr << "Point vals: " << ss.str() << std::endl;
            exit(1); 
        }
    }

    std::cerr << "SMQP interpolate() test: PASSED" << std::endl;
}


void testDivide() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q *= p;


    if (nvar == 0) {
        RationalNumber pCoef = p.leadingCoefficient();
        RationalNumber qCoef = q.leadingCoefficient();
        RationalNumber quo = qCoef / pCoef;
        RationalNumber rem = 0;

        SparseMultivariateRationalPolynomial quoPoly, remPoly;
        q.divide(p, quoPoly, remPoly);

        if (quoPoly != quo || remPoly != rem) {
            std::cerr << "SMQP divide() test: FAILED" << std::endl;
            exit(1);
        } 

        std::cerr << "SMQP divide() test: PASSED" << std::endl;
        return;

    }


    // q *= p;
    // q *= p;

    SparseMultivariateRationalPolynomial quo, rem;

    q.divide(p, quo, rem);

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTree = quo.convertToExpressionTree();
    ExpressionTree remTree = rem.convertToExpressionTree();

    //Due to the way the divide operation works (by returning the value of
    //interest in a new variable) we need some custom code to handle the 
    //variable naming and actually retrieving the value. Not just the proc return value.

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();


    std::vector<std::string> inputs;
    inputs.push_back(qTree.toMapleString() + ":");
    inputs.push_back("[" + pTree.toMapleString() + "]:");

    std::vector<Symbol> vars = p.ringVariables();
    std::string tdeg;
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += "):";
    inputs.push_back(tdeg);

    inputs.push_back(quoTree.toMapleString() + ":");
    inputs.push_back(remTree.toMapleString() + ":");

    std::string procStr = "Groebner:-NormalForm:";
    char* cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q2", 1));

    ALGEB result = EvalMapleProc(kv, testProc, 4, algebList[0], algebList[1], algebList[2], algebList[5]);

    char getQuoStr[] = "q2[1]:";
    ALGEB mapleQuo = EvalMapleStatement(kv, getQuoStr); 

    char compareFunc[] = "verify:";
    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);

    //algebList[3] is the quo, [4] is the rem, [5] is the namedVar.
    ALGEB comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, algebList[3]);
    ALGEB comp2 = EvalMapleProc(kv, cmpF, 2, result, algebList[4]);

    M_BOOL compBool = MapleToM_BOOL(kv, comp);
    M_BOOL comp2Bool = MapleToM_BOOL(kv, comp2);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, mapleQuo);
        std::cout << "SMQP smqp divide() test: FAILED" << std::endl;
        std::cout << "Got quotient:" << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    } 

    if (!comp2Bool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smqp divide() test: FAILED" << std::endl;
        std::cout << "Got remainder:" << retErr << std::endl;
        std::cout << "Expected: " << rem << std::endl;
        exit(1);
    }

    //test divisor with 1 term 
    p.randomPolynomial(nvar, 1, coefBound, sparsity, includeNeg);
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q *= p;

    quo.zero();
    rem.zero();
    q.divide(p, quo, rem);

    pTree = p.convertToExpressionTree();
    qTree = q.convertToExpressionTree();
    quoTree = quo.convertToExpressionTree();
    remTree = rem.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString() + ":");
    inputs.push_back("[" + pTree.toMapleString() + "]:");

    vars = p.ringVariables();
    tdeg = "";
    if (vars.size() > 1) {
        tdeg = "plex(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "plex(" + vars[0].toString();
    }
    tdeg += "):";
    inputs.push_back(tdeg);

    inputs.push_back(quoTree.toMapleString() + ":");
    inputs.push_back(remTree.toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "q3", 1));

    result = EvalMapleProc(kv, testProc, 4, algebList[0], algebList[1], algebList[2], algebList[5]);

    char getQuoStr2[] = "q3[1]:";
    mapleQuo = EvalMapleStatement(kv, getQuoStr2); 

    //algebList[3] is the quo, [4] is the rem, [5] is the namedVar.
    comp = EvalMapleProc(kv, cmpF, 2, mapleQuo, algebList[3]);
    comp2 = EvalMapleProc(kv, cmpF, 2, result, algebList[4]);

    compBool = MapleToM_BOOL(kv, comp);
    comp2Bool = MapleToM_BOOL(kv, comp2);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, mapleQuo);
        std::cout << "SMQP smzp divide() test: FAILED2" << std::endl;
        std::cout << "Got quotient:" << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    } 

    if (!comp2Bool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP smzp divide() test: FAILED2" << std::endl;
        std::cout << "Got remainder:" << retErr << std::endl;
        std::cout << "Expected: " << rem << std::endl;
        exit(1);
    }

    std::cout << "SMQP smqp divide() test: PASSED" << std::endl;
}

void testSMQPRemainder() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q *= p;
    q *= p;
    q *= p;

    SparseMultivariateRationalPolynomial rem = q%p;

    if (nvar == 0) {
        RationalNumber pCoef = p.leadingCoefficient();
        RationalNumber qCoef = q.leadingCoefficient();
        RationalNumber remCoef = qCoef % pCoef;
        if (rem != remCoef) {
            std::cerr << "SMQP remainder test: FAILED" << std::endl;
            exit(1);
        }
        std::cerr << "SMQP remainder test: PASSED" << std::endl;
        return;
    }

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree remTree = rem.convertToExpressionTree();

    std::vector<std::string> inputs;
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    std::vector<Symbol> vars = q.ringVariables();
    std::string tdeg;
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder test 1: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << rem;
        exit(1);
    }

//      - SMQP of different nvar
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    q *= p;
    q *= p;

    rem = q%p;

    pTree = p.convertToExpressionTree();
    qTree = q.convertToExpressionTree();
    remTree = rem.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = q.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder test 2: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << rem;
        exit(1);
    }

//      - SMQP that is zero
    q.zero();
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    rem = q%p;

    pTree = p.convertToExpressionTree();
    qTree = q.convertToExpressionTree();
    remTree = rem.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = p.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder test 3: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << rem;
        exit(1);
    }

//      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    p.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    rem = q%p;

    pTree = p.convertToExpressionTree();
    qTree = q.convertToExpressionTree();
    remTree = rem.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = p.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);

    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder test 4: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << rem;
        exit(1);
    }

//      - SMQP that is a constant and is nvar == 0
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    p.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    rem = q%p;

    pTree = p.convertToExpressionTree();
    qTree = q.convertToExpressionTree();
    remTree = rem.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = p.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder assignment test 5: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << q;
        exit(1);
    }


    // qTree -= remTree;
    // if (mapleTest->testIfZero(qTree, &retErr) == 0) {
    //     std::cerr << "SMQP remainder test 5: FAILED" << std::endl;
    //     std::cerr << "Got " << retErr << " but expected " << rem;
    //     exit(1);
    // }

    std::cerr << "SMQP remainder test: PASSED" << std::endl;
    return;
}

void testSMQPRemainderAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    SparseMultivariateRationalPolynomial q;
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q *= p;
    q *= p;
    q *= p;

    ExpressionTree qTree = q.convertToExpressionTree();
    
    if (nvar == 0) {
        RationalNumber pCoef = p.leadingCoefficient();
        RationalNumber qCoef = q.leadingCoefficient();
        q %= p;
        RationalNumber remCoef = qCoef % pCoef;
        if (q != remCoef) {
            std::cerr << "SMQP remainder assignment test: FAILED" << std::endl;
            exit(1);
        }
        std::cerr << "SMQP remainder assignment test: PASSED" << std::endl;
        return;
    }
 
    q %= p;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree remTree = q.convertToExpressionTree();

    std::vector<std::string> inputs;
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    std::vector<Symbol> vars = p.ringVariables();
    std::string tdeg;
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder assignment test 1: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << q;
        exit(1);
    }

    //      - SMQP of different nvar
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    q.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    q *= p;
    q *= p;

    qTree = q.convertToExpressionTree();
    q %= p;

    pTree = p.convertToExpressionTree();
    remTree = q.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = q.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder assignment test 2: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << q;
        exit(1);
    }

//      - SMQP that is zero
    q.zero();
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    qTree = q.convertToExpressionTree();
    q %= p;

    pTree = p.convertToExpressionTree();
    remTree = q.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = p.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder assignment test 3: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << q;
        exit(1);
    }

//      - SMQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    p.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    qTree = q.convertToExpressionTree();
    q%=p;

    pTree = p.convertToExpressionTree();
    remTree = q.convertToExpressionTree();

    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = p.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder assignment test 4: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << q;
        exit(1);
    }

//      - SMQP that is a constant and is nvar == 0
    q = SparseMultivariateRationalPolynomial(0);
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    p.randomPolynomial(nvar+2, numTerms, coefBound, sparsity, includeNeg);
    qTree = q.convertToExpressionTree();
    pTree = p.convertToExpressionTree();
    q%=p;

    remTree = q.convertToExpressionTree();


    inputs.clear();
    inputs.push_back(qTree.toMapleString());
    inputs.push_back("[" + pTree.toMapleString() + "]");

    vars = p.ringVariables();
    if (vars.size() > 1) {
        tdeg = "tdeg(" + vars[0].toString() + ",";
        for (int i = 1; i < vars.size()-1; ++i) {
            tdeg += vars[i].toString();
            tdeg += ",";
        }
        tdeg += vars[vars.size()-1].toString();
    } else {
        tdeg = "tdeg(" + vars[0].toString();
    }
    tdeg += ")";
    inputs.push_back(tdeg);
    if (mapleTest->testProcReturn("Groebner:-NormalForm", inputs, remTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP remainder assignment test 5: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << q;
        exit(1);
    }

    // qTree -= remTree;
    // if (mapleTest->testIfZero(qTree, &retErr) == 0) {
    //     std::cerr << "SMQP remainder assignment test 5: FAILED" << std::endl;
    //     std::cerr << "Got " << retErr << " but expected " << q;
    //     exit(1);
    // }

    std::cerr << "SMQP remainder assignment test: PASSED" << std::endl;
    return;
}

void testPseudoDivide() {

 
    if (nvar == 0) {   
        //pseudodivide is not defined then.
        return;
    }

    SparseMultivariateRationalPolynomial c;
    c.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial b;
    unsigned long numBTerms = numTerms;
    if (numTerms > 10) {
        numBTerms = numTerms/2;
    } else if (numTerms > 3) {
        numBTerms = numTerms - 2;
    }

    //Test same nvar for dividend and divisor
    b.randomPolynomial(nvar, numBTerms, coefBound, sparsity, includeNeg);
    while (b.isConstant()) {
        b.randomPolynomial(nvar, numBTerms, coefBound, sparsity, includeNeg);
    }

    SparseMultivariateRationalPolynomial quo, mult;
    SparseMultivariateRationalPolynomial rem = c.pseudoDivide(b, &quo, &mult);
   
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();
    
    ExpressionTree cTree = c.convertToExpressionTree();
    ExpressionTree bTree = b.convertToExpressionTree();
    ExpressionTree aTree = quo.convertToExpressionTree();
    ExpressionTree rTree = rem.convertToExpressionTree();
    ExpressionTree hTree = mult.convertToExpressionTree();

    Symbol var = b.leadingVariable();

    std::vector<std::string> inputs;
    inputs.push_back(cTree.toMapleString() + ":");
    inputs.push_back(bTree.toMapleString() + ":");
    inputs.push_back(var.toString() + ":"); //mvar

    inputs.push_back(aTree.toMapleString() + ":");
    inputs.push_back(rTree.toMapleString() + ":");
    inputs.push_back(hTree.toMapleString() + ":");

    std::string procStr = "prem:";
    char* cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;

    std::vector<ALGEB> algebList;
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "m2", 1));
    algebList.push_back(ToMapleName(kv, "q2", 1));

    ALGEB result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[6], algebList[7]);

    char compareFunc[] = "verify:";
    ALGEB cmpF = EvalMapleStatement(kv, compareFunc);

    char expandFunc[] = "expand:";
    ALGEB expandF = EvalMapleStatement(kv, expandFunc);

    //q vs a;
    algebList[7] = EvalMapleProc(kv, expandF, 1, algebList[7]); 
    ALGEB comp = EvalMapleProc(kv, cmpF, 2, algebList[7], algebList[3]);
    //result vs r;
    result = EvalMapleProc(kv, expandF, 1, result); 
    ALGEB comp2 = EvalMapleProc(kv, cmpF, 2, result, algebList[4]);
    //m vs e
    algebList[6] = EvalMapleProc(kv, expandF, 1, algebList[6]);
    ALGEB comp3 = EvalMapleProc(kv, cmpF, 2, algebList[6], algebList[5]);

    M_BOOL compBool = MapleToM_BOOL(kv, comp);
    M_BOOL comp2Bool = MapleToM_BOOL(kv, comp2);
    M_BOOL comp3Bool = MapleToM_BOOL(kv, comp3);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[7]);
        std::cout << "SMQP Pseudo divide() test: FAILED" << std::endl;
        std::cout << "Got quotient:" << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    } 

    if (!comp3Bool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[6]);
        std::cout << "SMQP Pesudo divide() test: FAILED" << std::endl;
        std::cout << "Got c multiplier:" << retErr << std::endl;
        std::cout << "Expected: " <<  mult  << std::endl;
        exit(1);
    }

    if (!comp2Bool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP Pesudo divide() test: FAILED" << std::endl;
        std::cout << "C: " << c << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "Got remainder:" << retErr << std::endl;
        std::cout << "Expected: " <<  rem << std::endl;
        exit(1);
    }

    //test different nvar

    c.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    int bnvar = nvar;
    if (nvar > 2) {
        bnvar -= 1;
    }
    b.randomPolynomial(bnvar, numBTerms, coefBound, sparsity, includeNeg);
    while (b.isConstant()) {
        b.randomPolynomial(bnvar, numBTerms, coefBound, sparsity, includeNeg);
    }

    std::vector<Symbol> symVec = c.ringVariables();
    if (bnvar != nvar) {
        c.setRingVariables(symVec);
        symVec.erase(symVec.begin());
        b.setRingVariables(symVec);
    }

    rem = c.pseudoDivide(b, &quo, &mult);

    cTree = c.convertToExpressionTree();
    bTree = b.convertToExpressionTree();
    aTree = quo.convertToExpressionTree();
    rTree = rem.convertToExpressionTree();
    hTree = mult.convertToExpressionTree();

    var = b.leadingVariable();

    inputs.clear();
    inputs.push_back(cTree.toMapleString() + ":");
    inputs.push_back(bTree.toMapleString() + ":");
    inputs.push_back(var.toString() + ":"); //mvar

    inputs.push_back(aTree.toMapleString() + ":");
    inputs.push_back(rTree.toMapleString() + ":");
    inputs.push_back(hTree.toMapleString() + ":");

    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "m3", 1));
    algebList.push_back(ToMapleName(kv, "q3", 1));

    //algeblist: 0 = c, 1 = b, 2 = mvar, 3 = a, 4 = r, 5 = e, 6 = 'm', 7 = 'q'

    result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[6], algebList[7]);

    //q vs a;
    algebList[7] = EvalMapleProc(kv, expandF, 1, algebList[7]); 
    comp = EvalMapleProc(kv, cmpF, 2, algebList[7], algebList[3]);
    //result vs r;
    result = EvalMapleProc(kv, expandF, 1, result); 
    comp2 = EvalMapleProc(kv, cmpF, 2, result, algebList[4]);
    //m vs e
    algebList[6] = EvalMapleProc(kv, expandF, 1, algebList[6]);
    comp3 = EvalMapleProc(kv, cmpF, 2, algebList[6], algebList[5]);

    compBool = MapleToM_BOOL(kv, comp);
    comp2Bool = MapleToM_BOOL(kv, comp2);
    comp3Bool = MapleToM_BOOL(kv, comp3);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[7]);
        std::cout << "SMQP Pseudo divide() test: FAILED2" << std::endl;
        std::cout << "Got quotient:" << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    } 

    if (!comp3Bool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[6]);
        std::cout << "SMQP Pesudo divide() test: FAILED2" << std::endl;
        std::cout << "Got c multiplier:" << retErr << std::endl;
        std::cout << "Expected: " <<  mult  << std::endl;
        exit(1);
    }

    if (!comp2Bool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP Pesudo divide() test: FAILED2" << std::endl;
        std::cout << "Got remainder:" << retErr << std::endl;
        std::cout << "Expected: " <<  rem << std::endl;
        std::cerr << "c: " << c << std::endl;
        std::cerr << "b: " << b << std::endl;
        std::cerr << "r: " << rem << std::endl;
        exit(1);
    }


    //test divisor size 1

    numBTerms = 1;
    c.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    b.randomPolynomial(nvar, numBTerms, coefBound, sparsity, includeNeg);
    while (b.isConstant()) {
        b.randomPolynomial(nvar, numBTerms, coefBound, sparsity, includeNeg);
    }

    rem = c.pseudoDivide(b, &quo, &mult);

    cTree = c.convertToExpressionTree();
    bTree = b.convertToExpressionTree();
    aTree = quo.convertToExpressionTree();
    rTree = rem.convertToExpressionTree();
    hTree = mult.convertToExpressionTree();

    var = b.leadingVariable();

    inputs.clear();
    inputs.push_back(cTree.toMapleString() + ":");
    inputs.push_back(bTree.toMapleString() + ":");
    inputs.push_back(var.toString() + ":"); //mvar

    inputs.push_back(aTree.toMapleString() + ":");
    inputs.push_back(rTree.toMapleString() + ":");
    inputs.push_back(hTree.toMapleString() + ":");


    algebList.clear();
    for (int i = 0; i < inputs.size(); ++i) {
        cstr = new char[inputs[i].length()+1];
        std::strcpy(cstr, inputs[i].c_str());
        ALGEB res = EvalMapleStatement(kv, cstr);
        algebList.push_back(res);
        delete[] cstr;
    }
    algebList.push_back(ToMapleName(kv, "m", 1));
    algebList.push_back(ToMapleName(kv, "q", 1));

    //algeblist: 0 = c, 1 = b, 2 = mvar, 3 = a, 4 = r, 5 = e, 6 = 'm', 7 = 'q'

    result = EvalMapleProc(kv, testProc, 5, algebList[0], algebList[1], algebList[2], algebList[6], algebList[7]);

    //q vs a;
    algebList[7] = EvalMapleProc(kv, expandF, 1, algebList[7]); 
    comp = EvalMapleProc(kv, cmpF, 2, algebList[7], algebList[3]);
    //result vs r;
    result = EvalMapleProc(kv, expandF, 1, result); 
    comp2 = EvalMapleProc(kv, cmpF, 2, result, algebList[4]);
    //m vs e
    algebList[6] = EvalMapleProc(kv, expandF, 1, algebList[6]);
    comp3 = EvalMapleProc(kv, cmpF, 2, algebList[6], algebList[5]);

    compBool = MapleToM_BOOL(kv, comp);
    comp2Bool = MapleToM_BOOL(kv, comp2);
    comp3Bool = MapleToM_BOOL(kv, comp3);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[7]);
        std::cout << "SMQP Pseudo divide() test: FAILED3" << std::endl;
        std::cout << "Got quotient:" << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    } 

    if (!comp3Bool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[6]);
        std::cout << "SMQP Pesudo divide() test: FAILED3" << std::endl;
        std::cout << "Got c multiplier:" << retErr << std::endl;
        std::cout << "Expected: " <<  mult  << std::endl;
        exit(1);
    }

    if (!comp2Bool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP Pesudo divide() test: FAILED3" << std::endl;
        std::cout << "Got remainder:" << retErr << std::endl;
        std::cout << "Expected: " <<  rem << std::endl;
        exit(1);
    }


    std::cout << "SMQP Pseudo divide() test: PASSED" << std::endl;
}




void testRatNum_TAddition() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    SparseMultivariateRationalPolynomial sum = p+r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree sumTree = sum.convertToExpressionTree();
    pTree += rTree;
    pTree -= sumTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    sum = p + r;

    pTree = p.convertToExpressionTree();
    pTree += rTree;
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    sum = p + r;

    SparseMultivariateRationalPolynomial sum2 = sum + r;

    pTree = sum.convertToExpressionTree();
    pTree += rTree;
    pTree -= sum2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //test ratNum + polynomial (in that order);

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = rTree;
    pTree += p.convertToExpressionTree();

    sum = r + p;

    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t addition test: PASSED" << std::endl;

}

void testRatNum_TAdditionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    ExpressionTree pTree = p.convertToExpressionTree();

    p += r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree sumTree = p.convertToExpressionTree();
    pTree += rTree;
    pTree -= sumTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 


    p.zero();
    pTree = p.convertToExpressionTree();
    
    p += r;

    pTree += rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p += r;

    pTree += rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t addition assignment test: PASSED" << std::endl;
}

void testRatNum_TSubtraction() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    SparseMultivariateRationalPolynomial diff = p-r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree diffTree = diff.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diffTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    diff = p - r;

    pTree = p.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    diff = p + r;

    SparseMultivariateRationalPolynomial diff2 = diff - r;

    pTree = diff.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diff2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //test ratNum - polynomial (in that order);
    
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = rTree;
    pTree -= p.convertToExpressionTree();

    diff = r - p;

    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t subtraction test: PASSED" << std::endl;
}

void testRatNum_TSubtractionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    ExpressionTree pTree = p.convertToExpressionTree();

    p -= r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree diffTree = p.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diffTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction assignment test 1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p -= r;

    pTree -= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction assignment test 2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p -= r;

    pTree -= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t subtraction assignment test 3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t subtraction assignment test: PASSED" << std::endl;
}

void testRatNum_TMultiplication() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    SparseMultivariateRationalPolynomial prod = p*r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree prodTree = prod.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    prod = p * r;

    pTree = p.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    prod = p + r;

    SparseMultivariateRationalPolynomial prod2 = prod * r;

    pTree = prod.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prod2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //test ratNum * polynomial (in that order);
    
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = rTree;
    pTree *= p.convertToExpressionTree();

    prod = r * p;

    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t multiplication test: PASSED" << std::endl;
}

void testRatNum_TMultiplicationAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    ExpressionTree pTree = p.convertToExpressionTree();

    p *= r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree prodTree = p.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p *= r;

    pTree *= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p *= r;

    pTree *= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t multiplication assignment test: PASSED" << std::endl;
}

void testRatNum_TDivision() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    SparseMultivariateRationalPolynomial quo = p/r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree quoTree = quo.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quoTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division test 1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    quo = p / r;

    pTree = p.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quo.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division test 2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    quo = p + r;

    SparseMultivariateRationalPolynomial quo2 = quo / r;

    pTree = quo.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quo2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division test 3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //test ratNum / polynomial (in that order);
    p = RationalNumber(r);
    mpq_mul(r, r, r);
    rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));

    pTree = rTree;
    pTree /= p.convertToExpressionTree();

    quo = r / p;

    //note, here we do pTree = quo.convertToExpressionTree on purpose! 
    //maple allows for rational functions of polynomails, so the simple
    //expression of r/poly is valid, but here we know that the quo is in fact 0
    if (p.degree() > 0) {
        pTree = quo.convertToExpressionTree();
    } else {
        pTree -= quo.convertToExpressionTree();
    }

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division test 4: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p.zero();
    p = RationalNumber(r);
    mpq_mul(r, r, r);
    rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    pTree = rTree;

    quo = r / p;

    pTree /= p.convertToExpressionTree();
    pTree -= quo.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division test 5: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
    
    std::cout << "SMQP ratNum_t division test: PASSED" << std::endl;
}

void testRatNum_TDivisionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    unsigned long val = (unsigned long) rand() % coefBound + 1;

    ratNum_t r;
    mpq_init(r);
    mpq_set_ui(r, val, 1ul);

    ExpressionTree pTree = p.convertToExpressionTree();

    p /= r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree quoTree = p.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quoTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p /= r;

    pTree /= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p /= r;

    pTree /= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP ratNum_t division assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP ratNum_t division assignment test: PASSED" << std::endl;
}

void testMPQ_ClassAddition() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    SparseMultivariateRationalPolynomial sum = p+r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree sumTree = sum.convertToExpressionTree();
    pTree += rTree;
    pTree -= sumTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    sum = p + r;

    pTree = p.convertToExpressionTree();
    pTree += rTree;
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    sum = p + r;

    SparseMultivariateRationalPolynomial sum2 = sum + r;

    pTree = sum.convertToExpressionTree();
    pTree += rTree;
    pTree -= sum2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //test ratNum + polynomial (in that order);

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = rTree;
    pTree += p.convertToExpressionTree();

    sum = r + p;

    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class addition test: PASSED" << std::endl;

}

void testMPQ_ClassAdditionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
       
    mpq_class r(rand() % coefBound + 1);

    ExpressionTree pTree = p.convertToExpressionTree();

    p += r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree sumTree = p.convertToExpressionTree();
    pTree += rTree;
    pTree -= sumTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p += r;

    pTree += rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p += r;

    pTree += rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class addition assignment test: PASSED" << std::endl;
}

void testMPQ_ClassSubtraction() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    SparseMultivariateRationalPolynomial diff = p-r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree diffTree = diff.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diffTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    diff = p - r;

    pTree = p.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    diff = p + r;

    SparseMultivariateRationalPolynomial diff2 = diff - r;

    pTree = diff.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diff2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //test ratNum - polynomial (in that order);
    
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = rTree;
    pTree -= p.convertToExpressionTree();

    diff = r - p;

    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class subtraction test: PASSED" << std::endl;
}

void testMPQ_ClassSubtractionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    ExpressionTree pTree = p.convertToExpressionTree();

    p -= r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree diffTree = p.convertToExpressionTree();
    pTree -= rTree;
    pTree -= diffTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p -= r;

    pTree -= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p -= r;

    pTree -= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class subtraction assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class subtraction assignment test: PASSED" << std::endl;
}

void testMPQ_ClassMultiplication() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    SparseMultivariateRationalPolynomial prod = p*r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree prodTree = prod.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    prod = p * r;

    pTree = p.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    prod = p + r;

    SparseMultivariateRationalPolynomial prod2 = prod * r;

    pTree = prod.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prod2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //test ratNum * polynomial (in that order);
    
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    pTree = rTree;
    pTree *= p.convertToExpressionTree();

    prod = r * p;

    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class multiplication test: PASSED" << std::endl;
}

void testMPQ_ClassMultiplicationAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    ExpressionTree pTree = p.convertToExpressionTree();

    p *= r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree prodTree = p.convertToExpressionTree();
    pTree *= rTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p *= r;

    pTree *= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p *= r;

    pTree *= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class multiplication assignment test: PASSED" << std::endl;
}

void testMPQ_ClassDivision() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    SparseMultivariateRationalPolynomial quo = p/r;

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree quoTree = quo.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quoTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division test 1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    quo = p / r;

    pTree = p.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quo.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division test 2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    quo = p + r;

    SparseMultivariateRationalPolynomial quo2 = quo / r;

    pTree = quo.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quo2.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division test 3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //test ratNum + polynomial (in that order);
    p = RationalNumber(r);
    r *= r;
    rTree = ExpressionTree(new ExprTreeNode(mpz_class(r)));

    pTree = rTree;
    pTree /= p.convertToExpressionTree();

    quo = r / p;

    //note, here we do pTree = quo.convertToExpressionTree on purpose! 
    //maple allows for rational functions of polynomails, so the simple
    //expression of r/poly is valid, but here we know that the quo is in fact 0
    if (p.degree() > 0) {
        pTree = quo.convertToExpressionTree();
    } else {
        pTree -= quo.convertToExpressionTree();
    }

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division test 4: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p.zero();
    p = RationalNumber(r);
    r *= r;
    rTree = ExpressionTree(new ExprTreeNode(r));
    pTree = rTree;

    quo = r / p;

    pTree /= p.convertToExpressionTree();
    pTree -= quo.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division test 5: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
    
    std::cout << "SMQP mpq_class division test: PASSED" << std::endl;
}

void testMPQ_ClassDivisionAssignment() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    
    mpq_class r(rand() % coefBound + 1);

    ExpressionTree pTree = p.convertToExpressionTree();

    p /= r;

    ExpressionTree rTree = ExpressionTree(new ExprTreeNode(mpq_class(r)));
    ExpressionTree quoTree = p.convertToExpressionTree();
    pTree /= rTree;
    pTree -= quoTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //now test spcial cases, likeone one the poly is 0 or a const.

    p.zero();
    pTree = p.convertToExpressionTree();
    
    p /= r;

    pTree /= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseMultivariateRationalPolynomial(0);
    p += r;
    pTree = p.convertToExpressionTree();

    p /= r;

    pTree /= rTree;
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP mpq_class division assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SMQP mpq_class division assignment test: PASSED" << std::endl;
}

void testLeadingVariable() {
    std::vector<Symbol> varArray = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t")};

    SparseMultivariateRationalPolynomial p(4);
    p.setRingVariables(varArray);

    Symbol varCheck("z");
    SparseMultivariateRationalPolynomial q(varCheck);
    p += q; //make p = 1*z;

    if (p.leadingVariable() != varCheck) {
        std::cerr << "SMQP leadingVariable() test: FAILED" << std::endl;
        std::cerr << "Expected " << varArray[0] << " but got " << p.leadingVariable() << std::endl;
        exit(1);
    }

    std::cerr << "SMQP leadingVariable() test: PASSED" << std::endl;
}

//NOTE: this test assumed that degree(std::string var) and leadingVariable() works. 
//If running tests in the order the tests are defined, then this is assured. 
void testLeadingVariableDegree() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    if (p.leadingVariableDegree() != p.degree(p.leadingVariable())) {
        std::cerr << "SMQP leadingVariableDegree() test: FAILED" << std::endl;
        std::cerr << "Expected " << p.degree(p.leadingVariable()) << " but got: " << p.leadingVariableDegree() << std::endl;
        exit(1);
    }

    std::cerr << "SMQP leadingVariableDegree() test: PASSED" << std::endl;
}

//NOTE: This test relies on coefficient(int, int*)
//This is fine if tests are performed in the order they are defined.
void testIsConstantTermZero() {
    SparseMultivariateRationalPolynomial p;
    p.zero();

    if (!p.isConstantTermZero()) {
        std::cerr << "SMQP isConstantTermZero() test: FAILED" << std::endl;
        exit(1);
    }

    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg); 
    int* degs = (int*) calloc(nvar, sizeof(int));
    RationalNumber coef = p.coefficient(nvar, degs);
    if ((coef == 0) != p.isConstantTermZero()) {
        std::cerr << "SMQP isConstantTermZero() test: FAILED" << std::endl;
        exit(1);
    }

    std::cerr << "SMQP isConstantTermZero() test: PASSED" << std::endl;
}

void testLeadingCoefficientInVariable() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    if (nvar == 0) {
        int i;
        SparseMultivariateRationalPolynomial leadingCoef = p.leadingCoefficientInVariable(Symbol("x"), &i);
        SparseMultivariateRationalPolynomial lCoef = p.leadingCoefficient();
        if (leadingCoef != lCoef || i != 0) {
            std::cerr << "SMQP leadingCoefficientInVariable() test: FAILED" << std::endl;
            std::cerr << "Expected degree 0 got: " << i << ", leading coef " << lCoef << " got: " << leadingCoef;
            exit(1);
        }
        std::cerr << "SMQP leadingCoefficientInVariable() test: PASSED" << std::endl;
        return;
    }


    int i = rand() % nvar;

    Symbol var = (p.ringVariables())[i];

    i = 0;
    SparseMultivariateRationalPolynomial leadingCoef = p.leadingCoefficientInVariable(var, &i);

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree expected = leadingCoef.convertToExpressionTree();
    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(var.toString());

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testProcReturn("lcoeff", inputs, expected.toMapleString(), &retErr) == 0) {
        std::cerr << "SMQP leadingCoefficientInVariable() test: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << leadingCoef << std::endl;
        std::cerr << p << std::endl;
        exit(1); 
    }

    if (mapleTest->testProcReturn("degree", inputs, std::to_string(i), &retErr) == 0) {
        std::cerr << "SMQP leadingCoefficientInVariable() test: FAILED" << std::endl;
        std::cerr << "Got deg(p," << var << "): " << retErr << " but expected " << i << std::endl;
        std::cerr << p << std::endl;
        exit(1);
    }

    std::cerr << "SMQP leadingCoefficientInVariable() test: PASSED" << std::endl;
}


ExpressionTree SUPconvertToExpressionTree(SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> sup) {

    Integer d = sup.degree();
    Symbol var = sup.variable();
    ExpressionTree prev;
    bool first = 1;

    SparseMultivariateRationalPolynomial coef;
    for (int i = d.get_si(); i >= 0; --i) {
        coef = sup.coefficient(i);
        if (coef.isZero()) {
            continue;
        }

        ExpressionTree coefTree = coef.convertToExpressionTree();
        ExpressionTree monomial = ExpressionTree(new ExprTreeNode(var));
        monomial ^= ExpressionTree(new ExprTreeNode(i));
        coefTree *= monomial;

        if (first) {
            prev = coefTree;
            first = 0;
        } else {
            prev += coefTree;
        }
    }

    return prev;
}


void testConvertToSUP() {

    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    Symbol var;
    if (nvar == 0) {
        var = "x";
    } else {
        int i = rand() % nvar;
        var = (p.ringVariables())[i];
    }

    SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> sup =  p.convertToSUP(var);

    ExpressionTree supTree = SUPconvertToExpressionTree(sup);

    std::vector<std::string> inputs; 
    inputs.push_back(p.convertToExpressionTree().toMapleString());

    //NOTE: we expand the SUP, because otherwise maple's verify will not properly
    //compare the two polynomials as the do not live in the same base ring.
    inputs.push_back("expand(" + supTree.toMapleString() + ")");

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testProcReturn("verify", inputs, "true", &retErr) == 0) {
        std::cerr << "SMQP convertToSUP() test: FAILED" << std::endl;
        std::cerr << "Expected equal values between: ";
        std::cerr << "p: " << p << std::endl;
        std::cerr << "and" << std::endl;
        std::cerr << "sup: " << sup << std::endl;
        exit(1);
    }

    std::cerr << "SMQP convertToSUP() test: PASSED" << std::endl;
}

void testNegate() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
    SparseMultivariateRationalPolynomial negP = p;
    negP.negate();
    ExpressionTree negPTree = negP.convertToExpressionTree();
    ExpressionTree pTree = p.convertToExpressionTree();
    pTree += negPTree;
    std::string retErr;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP negate() test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } else {
        std::cout << "SMQP negate() test: PASSED" << std::endl;
    }
}

void testDeepCopy() {
    SparseMultivariateRationalPolynomial p;
    p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

    SparseMultivariateRationalPolynomial q = p.deepCopy();

    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();

    pTree -= qTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SMQP deepCopy() test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p.zero();

    if (p == q) {
        std::cerr << "SMQP deepCopy() test: FAILED" << std::endl;
        std::cerr << "Modifying original caused change in copy!" << std::endl;
        exit(1);
    }

    std::cerr << "SMQP deepCopy() test: PASSED" << std::endl;
}

void testStraightLineProgram() { 

for(int testCase = 0; testCase < 4; ++testCase) {

    SparseMultivariateRationalPolynomial p;
    switch (testCase) {
        case 0: {
            //random poly
            p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
            break;
        }
        case 1: {
            p = SparseMultivariateRationalPolynomial(nvar);
            p.zero();
            break;            
        }
        case 2: {
            p = SparseMultivariateRationalPolynomial(nvar);
            p.zero();
            p += mpq_class(rand() % coefBound);
            break;
        }
        case 3: {
            p = SparseMultivariateRationalPolynomial(0);
            p.zero();
            p += mpq_class(rand() % coefBound);
            break;
        }
    }

    std::vector<Symbol> vars = p.ringVariables();
    int varsSize = vars.size();
    if (varsSize == 0) {
        //just so that passing to maple is easy
        vars.push_back(Symbol("x_1"));
        varsSize = 1;
    }
    
    //get its slp 
    p.straightLineProgram();

    //prepare maple proc
    std::stringstream ss;
    ss << "proc(";
    for (int i = 0; i < vars.size()-1; ++i) {
        ss << vars[i] << ", ";
    }
    ss << vars[vars.size()-1] << ") " << std::endl;
    ss << " local ";
    for (int i = 0; i < p.slp.size()-1; ++i) {
        ss << "r_" << i << ", ";
    }
    ss << "r_" << p.slp.size()-1 << ": ";
    p.printSLP(ss);
    ss << "end proc:" << std::endl;

    //get maple
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    MKernelVector kv = mapleTest->getMKernelVector();

    //get algeb for slp proc;
    std::string procStr = ss.str();
    char* procCstr = new char[procStr.length()+1];
    std::strcpy(procCstr, procStr.c_str());
    ALGEB procAlgeb = EvalMapleStatement(kv, procCstr);
    delete[] procCstr;

    //get algeb for our poly
    std::string polyStr = p.convertToExpressionTree().toMapleString() + ":";
    char* polyCstr = new char[polyStr.length()+1];
    std::strcpy(polyCstr, polyStr.c_str());
    ALGEB polyAlgeb = EvalMapleStatement(kv, polyCstr);
    delete[] polyCstr;

    //get eval and verify algebs
    char evalCStr[6] = "eval:";
    ALGEB evalAlgeb = EvalMapleStatement(kv, evalCStr);
    char verifyCStr[8] = "verify:";
    ALGEB verifyAlgeb = EvalMapleStatement(kv, verifyCStr);


    //now check SLP return vs maple's eval return for 10 random sets of values;
    int vals[varsSize];
    int maxVal = 100;

    //test 100 random sets of vals. Each iteration is quite fast.
    for (int i = 0; i < 100; ++i) {
        std::vector<ALGEB> algebList;
        for (int j = 0; j < varsSize; ++j) {
            if (i == 0) {
                vals[j] = 0;    
            } else {
                vals[j] = rand() % maxVal;
                if (rand() % 2) {
                    vals[j] = -vals[j];
                }
            }
            algebList.push_back(ToMapleInteger(kv, vals[j]));
        }
        std::stringstream varss;
        
        //build vars list for these vals
        varss << "[";
        for (int j = 0; j < varsSize-1; ++j) {
            varss << vars[j] << "=" << vals[j] << ", ";
        }
        varss << vars[varsSize-1] << "=" << vals[varsSize-1] << "]:";
        std::string varsStr = varss.str();
        char* valsCstr = new char[varsStr.length()+1];
        std::strcpy(valsCstr, varsStr.c_str());
        ALGEB varlistAlgeb = EvalMapleStatement(kv, valsCstr);
        delete[] valsCstr;

        ALGEB evalRes = EvalMapleProc(kv, evalAlgeb, 2, polyAlgeb, varlistAlgeb);

        ALGEB procRes;
        switch (algebList.size()) {
            case 1: {
                procRes = EvalMapleProc(kv, procAlgeb, 1, algebList[0]);
                break;
            }
            case 2: {
                procRes = EvalMapleProc(kv, procAlgeb, 2, algebList[0], algebList[1]);
                break;
            }
            case 3: {
                procRes = EvalMapleProc(kv, procAlgeb, 3, algebList[0], algebList[1], algebList[2]);
                break;
            }
            case 4: {
                procRes = EvalMapleProc(kv, procAlgeb, 4, algebList[0], algebList[1], algebList[2], algebList[3]);
                break;
            }
            case 5: {
                procRes = EvalMapleProc(kv, procAlgeb, 5, algebList[0], algebList[1], algebList[2], algebList[3], algebList[4]);
                break;
            }
            case 6: {
                procRes = EvalMapleProc(kv, procAlgeb, 6, algebList[0], algebList[1], algebList[2], algebList[3], algebList[4], algebList[5]);
                break;
            }
        }

        ALGEB verifyRes = EvalMapleProc(kv, verifyAlgeb, 2, evalRes, procRes);
        M_BOOL sameRes = MapleToM_BOOL(kv, verifyRes);

        if (!sameRes) {
            std::cerr << "SMQP straightLineProgram() and printSLP() test: FAILED" << std::endl;
            std::cerr << "Got SLP val  " << mapleTest->algebToString(kv, procRes) << std::endl;
            std::cerr << "Got poly val " << mapleTest->algebToString(kv, evalRes) << std::endl;
            std::cerr << "poly: " << p << std::endl;
            std::cerr << "slp: " << std::endl;
            p.printSLP(std::cerr);
            exit(1);
        }
    }//testing 100 vals for

}//outer for

    std::cerr << "SMQP straightLineProgram() test: PASSED" << std::endl;
}

void testFromString() {

    std::stringstream ss;
    SparseMultivariateRationalPolynomial p;
    
    //test "0"
    
    std::string fromString = "0";
    p.fromString(fromString);

    p.print(ss);

    std::string toString = ss.str();
    ss.str("");

    if (toString != fromString) {
        std::cerr << "SMQP: fromString() test: FAILED!" << std::endl;
        std::cerr << "Expected to get: " << fromString << " but got: " << toString << std::endl;
        exit(1);
    }

    //test univariate

    fromString = "x^10 + 2*x^4 + 12*x";
    p.fromString(fromString);

    ss.str("");
    p.print(ss);
    toString = ss.str();

    if (toString != fromString) {
        std::cerr << "SMQP: fromString() test: FAILED!" << std::endl;
        std::cerr << "Expected to get: " << fromString << " but got: " << toString << std::endl;
        exit(1);
    }

    //test multivariate

    //test multivariate with duplicated symbols in power product

    //test integers

    //test rational numbers



    std::cerr << "SMQP: fromString() test: PASSED!" << std::endl;
}

void testIstreamOperator() {

    std::stringstream ss("x^10 + 2*x^4 + 12*x");

    SparseMultivariateRationalPolynomial p;

    ss >> p;

    std::stringstream out_ss;
    p.print(out_ss);

    if (out_ss.str() != ss.str()) {
        std::cerr << "SMQP: istream operator test: FAILED!" << std::endl;
        std::cerr << "Expected: " << ss.str() << " but got: " << out_ss.str() << std::endl;
        exit(1);
    }


    std::cerr << "SMQP: istream operator test: PASSED!" << std::endl;
}

int main(int argc, char *argv[]) {
        std::cerr << "Time: " << time(NULL) << std::endl;

        if (argc > 1 && atol(argv[1]) >= 0) {
            nvar = atol(argv[1]);
        }
        if (argc > 2 && atol(argv[2]) > 0) { 
            numTerms = atol(argv[2]);
        }
        if (argc > 3 && atol(argv[3]) > 0) {
            coefBound = atol(argv[3]);
        }
        if (argc > 4 && atol(argv[4]) > 1) {
            sparsity = atol(argv[4]);
        }
        if (argc > 5 && atoi(argv[5]) >= 0) {
            includeNeg = atoi(argv[5]);
        }

        if (nvar == 0) {
            numTerms = 1;
        }

        // char* syms[] = {"x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10", "x_11", "x_12", "x_13", "x_14", "x_15", "x_16", "x_17", "x_18", "x_19", "x_20", "x_21", "x_22", "x_23", "x_24", "x_25", "x_26", "x_27", "x_28", "x_29", "x_30"};
        // AltArr_t* aa = buildRandomPoly_AA_unpk(nvar, numTerms, coefBound, sparsity, includeNeg);
        // tryPackExponentVectors_AA_inp(aa);
        // fprintf(stderr, "aa->unpacked: %d\n", aa->unpacked);
        // // printPoly_AA(stderr, aa, syms, nvar);
        // // fprintf(stderr, "\n");
        // return 0;

        // SparseMultivariateRationalPolynomial p, q, s;
        // p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        // q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        // s.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
        // p *= s;
        // q *= s;
        // RationalNumber contp = p.content();
        // RationalNumber contq = q.content();

        // timer_id id = start_timer();
        // SparseMultivariateRationalPolynomial g = p.gcd(q);
        // timer_time elapsed1 = elapsed_time(&id);
        // double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
        // std::cerr << "time: " << time << std::endl;
        // return 0;

        testDefaultConstructor();
        testNvarConstructor();
        testIdentityConstructor();
        testStringConstructor();
        testCopyConstructor();
        testMoveConstructor();
        testSMZPConstructor();
        testIntegerConstructor();
        testRationalConstructor();
        testDUQPConstructor();
        testSUPConstructor();
        testIsZero();
        testZero();
        testIsOne();
        testOne();
        testIsNegativeOne();
        testNegativeOne();
        testIsConstant();
        testUnitCanonical();
        testCopyAssignment();
        testMoveAssignment();
        testSMQPAddition();
        testSMQPAdditionAssignment();
        testUnaryNegative();
        testSMQPSubtraction();
        testSMQPSubtractionAssignment();
        testSMQPMultiplication();
        testSMQPMultiplicationAssignment();
        testSMQPDivision();
        testSMQPDivisionAssignment();
        testExponentiation();
        testExponentiationAssignment();
        testGCD(); 
        testContent();
        testPrimitivePart();
        testPrimitivePartSMZP();
        testSquareFree();
        testSquareFreePart();
        testNumberOfRingVariables();
        testNumberOfVariables();
        testNumberOfTerms();
        testDegree();
        testLeadingCoefficient();
        testCoefficient();
        testSetVariableNames();
        testIsEqual();
        testDerivative();
        testIntegrate();
        testEvaluate();
        testInterpolate();
        testDivide();
        testSMQPRemainder();
        testSMQPRemainderAssignment();
        testPseudoDivide();
        testInitial();
        testHead();
        testTail();
        testRatNum_TAddition();
        testRatNum_TAdditionAssignment();
        testRatNum_TSubtraction();
        testRatNum_TSubtractionAssignment();
        testRatNum_TMultiplication();
        testRatNum_TMultiplicationAssignment();
        testRatNum_TDivision();
        testRatNum_TDivisionAssignment();
        testMPQ_ClassAddition();
        testMPQ_ClassAdditionAssignment();
        testMPQ_ClassSubtraction();
        testMPQ_ClassSubtractionAssignment();
        testMPQ_ClassMultiplication();
        testMPQ_ClassMultiplicationAssignment();
        testMPQ_ClassDivision();
        testMPQ_ClassDivisionAssignment();
        testLeadingVariable();
        testLeadingVariableDegree();
        testIsConstantTermZero();
        testLeadingCoefficientInVariable();
        testConvertToSUP();
        testNegate();
        testDeepCopy();
        testStraightLineProgram();
        testFromString();
        testIstreamOperator();

        return 0;
}
