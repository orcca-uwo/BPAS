 #include <bpas.h>
#include "../../include/RationalNumberPolynomial/urpolynomial-altarr.hpp"
#include "../MapleTestTool/MapleTestTool.hpp"
#include <math.h>
#include <fstream>

using namespace std;

int nterms=2; unsigned long int coefBound=3; degree_t sparsity=2; int includeNeg=0;

void testIdentityConstructor() {
    SparseUnivariateRationalPolynomial p('x');

    ExpressionTree pTree = p.convertToExpressionTree();
    pTree -= ExpressionTree(new ExprTreeNode(Symbol('x')));

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SUQP identity constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    p = SparseUnivariateRationalPolynomial(Symbol("xas_12312"));
    pTree = p.convertToExpressionTree();
    pTree -= ExpressionTree(new ExprTreeNode(Symbol("xas_12312")));
    
 /*   std::cerr<<"P="<<p<<std::endl;
   std::cerr<<"ptree="<<pTree.toMapleString()<<std::endl; */
          


    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SUQP identity constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SUQP identity constructor test: PASSED" << std::endl;
    
}

void testCopyConstructor() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(  nterms,  coefBound, sparsity,  includeNeg);

    ExpressionTree pTree = p.convertToExpressionTree();

    SparseUnivariateRationalPolynomial q = p;

/*  std::cerr << "p : " << p << std::endl;
 std::cerr << "q : " << q << std::endl; */


    pTree -= q.convertToExpressionTree();
    
    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SUQP copy constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SUQP copy constructor test: PASSED" << std::endl;
}

void testMoveConstructor() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(  nterms,  coefBound, sparsity,  includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();

    mpq_class val = mpq_class(rand() % coefBound + 1);
    SparseUnivariateRationalPolynomial q = p + val; 

    pTree += ExpressionTree(new ExprTreeNode(mpq_class(val)));

    pTree -= q.convertToExpressionTree();

    std::string retErr;
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    if (mapleTest->testIfZero(pTree, &retErr) == 0) {
        std::cerr << "SUQP move constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SUQP move constructor test: PASSED" << std::endl;   
}

void testconvertToExpressionTree() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(  nterms,  coefBound, sparsity,  includeNeg);

    ExpressionTree pTree = p.convertToExpressionTree();
    SparseUnivariateRationalPolynomial q(p);

    ExpressionTree qTree = q.convertToExpressionTree();
       ExpressionTree diff=pTree- qTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;

    if(mapleTest->testIfZero(diff, &retErr) == 0) {

        std::cerr << "SUQP  constructor test: FAILED" << std::endl;
        std::cerr << "Expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cerr << "SUQP  constructor test: PASSED" << std::endl;   
}
void testIsZero () {
    SparseUnivariateRationalPolynomial p1;
    SparseUnivariateRationalPolynomial p2;
    p1.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    p2 = p1;
    if (p1.isZero() == 0 && (p1-p2).isZero() == 1) {
        std::cerr << "SUQP isZero() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SUQP isZero() test: FAILED" << std::endl;
    exit(1);
}
void testZero() {
    SparseUnivariateRationalPolynomial p1;
    p1.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    bool before = p1.isZero();
   
    p1.zero();
    bool after = p1.isZero();
    if (!before && after) {
        std::cerr << "SUQP zero() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SUQP zero() test: FAILED" << std::endl;
    exit(1);
}

void testIsOne() {
    SparseUnivariateRationalPolynomial p1;
    SparseUnivariateRationalPolynomial p2;
    p1.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    p2 = p1;
    p1 += mpq_class(1);

    if (p1.isOne() == 0 && (p1-p2).isOne() == 1) {
        std::cerr << "SUQP isOne() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SUQP isOne() test: FAILED" << std::endl;
    exit(1);
}
void testOne() {
    SparseUnivariateRationalPolynomial p1;
    p1.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    bool before = (p1.isOne() && !(nterms == 1 && p1.leadingCoefficient() == 1));
    p1.one();
    bool after = p1.isOne();

    if (!before && after) {
        std::cerr << "SUQP one() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SUQP one() test: FAILED" << std::endl;
    exit(1);
}
void testIsConstant() {
    SparseUnivariateRationalPolynomial p1;
    p1.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    SparseUnivariateRationalPolynomial p2 = p1;
    p2 += mpq_class(rand() % coefBound + 1);
    if ((nterms<= 1 || p1.isConstant() == 0) && (p1-p2).isConstant()) {
        std::cerr << "SUQP isConstant() test: PASSED" << std::endl;
        return;
    }

    std::cerr << "SUQP isConstant() test: FAILED" << std::endl;
    exit(1);
}
void testUnitCanonical() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    RationalNumber lc = p.leadingCoefficient();
    while(lc.isOne()) {
        p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
        lc = p.leadingCoefficient();
    }

    SparseUnivariateRationalPolynomial u;
    SparseUnivariateRationalPolynomial v;
    SparseUnivariateRationalPolynomial unit = p.unitCanonical(&u, &v);

    if (unit.leadingCoefficient() != RationalNumber(1)) {
        std::cerr << "SUQP unitCanonical() test: FAILED" << std::endl;
        std::cerr << "Unit canonical's leading coefficient should be 1!" << std::endl;
        std::cerr << "p: " << p << std::endl;
        std::cerr << "unitCanonical: " << unit << std::endl;
        exit(1);
    }

    if (!(u*lc).isOne()) {
        std::cerr << "SUQP unitCanonical() test: FAILED" << std::endl;
        std::cerr << "lc(p)*u should equal 1 but got: " << (u*lc) << std::endl;
        exit(1);
    }

    if (!(u*v).isOne()) {
        std::cerr << "SUQP unitCanonical() test: FAILED" << std::endl;
        std::cerr << "v should be the inverse of u but their product is: " << (u*v) << std::endl;
        exit(1); 
    }

    std::cerr << "SUQP unitCanonical() test: PASSED" << std::endl;
}
void testCopyAssignment() {
    SparseUnivariateRationalPolynomial p1;
    p1.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    SparseUnivariateRationalPolynomial p2 = p1;
    if (p1 != p2) {
        std::cerr << "SUQP copy assignment test: FAILED" << std::endl;
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
        std::cerr << "SUQP copy assignment from RationalNumber test: FAILED" << std::endl; 
        std::cerr << "Integer: " << it << " SUQP: " << p1;
        exit(1);
    }

    std::cerr << "SUQP copy assignment test: PASSED" << std::endl;
}

void testSUQPAddition() {
    SparseUnivariateRationalPolynomial p;
p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);    
    //      - SUQP of same nvar
    SparseUnivariateRationalPolynomial q;
q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    SparseUnivariateRationalPolynomial sum = p+q;

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
        std::cerr << "SUQP suqp addition test1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    sum = p+q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition test2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    q.zero();
    sum = p+q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition test3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    sum = p + q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition test4: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant and is nvar == 0    
    q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    sum = p + q;

    pTree = p.convertToExpressionTree();
    pTree += q.convertToExpressionTree();
    pTree -= sum.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition test5: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp addition test: PASSED" << std::endl;
}
void testSUQPAdditionAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();
    
    SparseUnivariateRationalPolynomial q;
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cerr << "SUQP suqp addition assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    pTree = p.convertToExpressionTree();
    p+=q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    q.zero();
    pTree = p.convertToExpressionTree();
    p+=q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p += q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant and is nvar == 0    
    q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p += q;

    pTree += q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp addition assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    std::cout << "SUQP suqp addition assignment test: PASSED" << std::endl;

}

void testSUQPSubtraction() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    SparseUnivariateRationalPolynomial q;
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    SparseUnivariateRationalPolynomial diff = p-q;

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
        std::cerr << "SUQP suqp subtraction test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
 
    //      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    diff = p-q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    q.zero();
    diff = p-q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    diff = p - q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant and is nvar == 0    
    q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    diff = p - q;

    pTree = p.convertToExpressionTree();
    pTree -= q.convertToExpressionTree();
    pTree -= diff.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp subtraction test: PASSED" << std::endl;
}

void testSUQPSubtractionAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();
    
    SparseUnivariateRationalPolynomial q;
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cerr << "SUQP suqp subtraction assignment test 1: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

        //      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    pTree = p.convertToExpressionTree();
    p-=q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();



    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction assignment test 2: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    q.zero();
    pTree = p.convertToExpressionTree();
    p-=q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction assignment test 3: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p -= q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp subtraction assignment test 4: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant and is nvar == 0    
    q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p -= q;

    pTree -= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP smup subtraction assignment test 5: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp subtraction assignment test: PASSED" << std::endl;  
}
void testSUQPMultiplication() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    SparseUnivariateRationalPolynomial q;
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    SparseUnivariateRationalPolynomial prod = p*q;

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
        std::cerr << "SUQP suqp multiplication test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    prod = p*q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    q.zero();
    prod= p*q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    prod = p * q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant and is nvar == 0    
    q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    prod = p * q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }




    std::cout << "SUQP suqp multiplication test: PASSED" << std::endl;
}
void testSUQPMultiplicationAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    ExpressionTree pTree = p.convertToExpressionTree();
    
    SparseUnivariateRationalPolynomial q;
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cerr << "SUQP suqp multiplication assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    } 

    //      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    pTree = p.convertToExpressionTree();
    p*=q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    q.zero();
    pTree = p.convertToExpressionTree();
    p*=q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p *= q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is a constant and is nvar == 0    
    q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p *= q;

    pTree *= q.convertToExpressionTree();
    pTree -= p.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication assignment test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp multiplication assignment test: PASSED" << std::endl;   
}
void testSUQPDivision() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    SparseUnivariateRationalPolynomial q = p;
    p *= p;
    p *= p;
    p *= p;

    SparseUnivariateRationalPolynomial quo = p/q;

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
        std::cout << "SUQP suqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

 //      - SUQP that is zero
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }
 
//      - SUQP that is a constant
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SUQP that is a constant and is nvar == 0

    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    q = SparseUnivariateRationalPolynomial();
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
        std::cout << "SUQP suqp division test: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp division test: PASSED" << std::endl;
}
void testSMQPDivisionAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    SparseUnivariateRationalPolynomial q = p;
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
        std::cout << "SUQP suqp division assignment test: FAILED6" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }


//      - SUQP of different nvar
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division assignment test: FAILED7" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

//      - SUQP that is zero
   /*  q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division assignment test: FAILED8" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

//      - SUQP that is a constant
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division assignment test: FAILED9" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    } */

//      - SUQP that is a constant and is nvar == 0

    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    q = SparseUnivariateRationalPolynomial();
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
        std::cout << "SUQP suqp division assignment test: FAILED10" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << p << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp division assignment test: PASSED" << std::endl;
}
void testExponentiation() {
    SparseUnivariateRationalPolynomial   p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    int exp = (rand() % 10) + 1; //10 for practical purposes;
    
    SparseUnivariateRationalPolynomial prod = p ^ exp;
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
        std::cerr << "SUQP suqp exponentiation test: FAILED1" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//      - SUQP that is zero
    p.zero();

    prod = p^exp;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation test: FAILED2" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//      - SUQP that is a constant
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    prod = p^exp;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation test: FAILED3" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }    

//      - SUQP that is a constant and is nvar == 0
    p = SparseUnivariateRationalPolynomial();
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    prod = p^exp;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation test: FAILED4" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//       - to exponent 0;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    prod = p^0;
    pTree = p.convertToExpressionTree();
    prodTree = prod.convertToExpressionTree();
    pTree ^= ExpressionTree(new ExprTreeNode(0));
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation test: FAILED5" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp exponentiation test: PASSED" << std::endl;
}
void testExponentiationAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

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
        std::cerr << "SUQP sUqp exponentiation assignment test: FAILED1" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    //      - SUQP that is zero
    p.zero();
    pTree = p.convertToExpressionTree();
    p^=exp;
    prodTree = p.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation assignment test: FAILED2" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//      - SUQP that is a constant
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p^=exp;
    prodTree = p.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation assignment test: FAILED3" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }    

//      - SUQP that is a constant and is nvar == 0
    p = SparseUnivariateRationalPolynomial();
    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    pTree = p.convertToExpressionTree();
    p^=exp;
    prodTree = p.convertToExpressionTree();
    pTree ^= constTree;
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation assignment test: FAILED4" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

//       - to exponent 0;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    pTree = p.convertToExpressionTree();
    p^=0;
    prodTree = p.convertToExpressionTree();
    pTree ^= ExpressionTree(new ExprTreeNode(0));
    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp exponentiation assignment test: FAILED5" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp exponentiation assignment test: PASSED" << std::endl;
}
/* 
void testGCD() {

    SparseUnivariateRationalPolynomial p,q;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    SparseUnivariateRationalPolynomial g = p.gcd(q);
    
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
    SparseUnivariateRationalPolynomial r;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    q = p;
    for (int i = 0; i < nterms; ++i) {
        r.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
        r = r[0];
        p *= r;

        r.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
        r = r[0];
        q *= r;
    }

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

    

    //Test non-trivial
    p.zero();
    q.zero();
    for (int i = 0; i < nterms; ++i) {
        r.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
        r = r[0];
        p += r;

        r.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
        r = r[0];
        q += r;
    }
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
 */
/* void testContent() {

    SparseUnivariateRationalPolynomial p;
    
    p.zero();
    RationalNumber content = p.content();
    if (content != 0) {
        std::cerr << "SUQP content() test FAILED" << std::endl;
        std::cerr << "Content should be 0 since p is 0 but got: " << content << std::endl;
        exit(1);
    }

    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    content = p.content();

    if (p.isZero() && content != 0) {
        std::cerr << "SUQP content() test FAILED" << std::endl;
        std::cerr << "Content should be 0 but got: " << content << std::endl;
        exit(1);
    }

    if (content < 0) {
        if (p.leadingCoefficient() > 0) {
            std::cerr << "SUQP content() test FAILED" << std::endl;
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
        varList = 'x';
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
 */
//////////////////////////////////////////////////////////
/* void testDegree() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

        int degree = p.degree().get_si();
        ExpressionTree expected = ExpressionTree(new ExprTreeNode(degree));
        if (mapleTest->testProcReturn("degree", inputs, expected.toMapleString(), &retErr) == 0) {
            std::cerr << "SUQP degree() test: FAILED" << std::endl;
            std::cerr << "Got " << retErr << " but expected " << degree << std::endl;
            exit(1);
        
    }

    std::cerr << "SUQP degree() test: PASSED" << std::endl;
} */
//////////////////////////////////////////////////////////////////
void testLeadingCoefficient() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    RationalNumber leadingCoef = p.leadingCoefficient();
    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree expected = ExpressionTree(new ExprTreeNode(leadingCoef.get_mpq()));
    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    if (mapleTest->testProcReturn("lcoeff", inputs, expected.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP leadingCoefficient() test: FAILED" << std::endl;
        std::cerr << "Got " << retErr << " but expected " << leadingCoef << std::endl;
        std::cerr << p << std::endl;
        exit(1); 
    }

    std::cerr << "SUQP leadingCoefficient() test: PASSED" << std::endl;
}

void testCoefficient() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    int degs=2;

    RationalNumber got1 = p.coefficient( degs);

    RationalNumber expected(27, 34);
    if (got1 == expected){
        expected = RationalNumber(29, 34);
    }
    p.setCoefficient(degs, expected);

    RationalNumber got = p.coefficient(degs);

    if ( (got1 != 0 && got1 == got) || expected != got) {
        std::cerr << "SUQP coefficient() and setCoefficient() test: FAILED" << std::endl;
        std::cerr << "Got " << got << " but expected " << expected << std::endl;
        exit(1);
    }

    std::cerr << "SUQP coefficient() and setCoefficient() test: PASSED" << std::endl;
}
void testDerivative() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    

    Symbol  syms = p.variable();
  

    SparseUnivariateRationalPolynomial dp = p.derivative();
    SparseUnivariateRationalPolynomial d3p = p.derivative( 3);

    SparseUnivariateRationalPolynomial q = p;
    ExpressionTree pTree = p.convertToExpressionTree();

    p.differentiate();
    q.differentiate( 3);


    if (dp != p || d3p != q) {
        std::cerr << "SUQP differentiate test() FAILED" << std::endl;
        std::cerr << "Should be equal: " << std::endl << p << std::endl << dp << std::endl;
        std::cerr << "Should be equal: " << std::endl << q << std::endl << d3p << std::endl;
        exit(1);
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString());

    ExpressionTree resTree = p.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP derivative() test: FAILED1" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString() + "$" + std::to_string(3));

    resTree = q.convertToExpressionTree();
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP derivative() test: FAILED2" << std::endl;
        std::cerr << "Expected: " << q << " but got: " << retErr << std::endl;
        exit(1);
    }

    //zero
    p.zero();
    pTree = p.convertToExpressionTree();
    p.differentiate();
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString());
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP derivative() test: FAILED3" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    //constant
    p += RationalNumber(37, 13);
    pTree = p.convertToExpressionTree();
    p.differentiate();
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString());
    if(mapleTest->testProcReturn("diff", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP derivative() test: FAILED4" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    

    std::cerr << "SUQP derivative() test: PASSED" << std::endl;
}

void testIntegrate() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

       

   Symbol  syms = p.variable();
  

    SparseUnivariateRationalPolynomial dp = p.Integral(1);
    SparseUnivariateRationalPolynomial d3p = p.Integral(3);

    SparseUnivariateRationalPolynomial q = p;
    ExpressionTree pTree = p.convertToExpressionTree();

    p.integrate(1);
    q.integrate( 3);

    if (dp != p || d3p != q) {
        std::cerr << "SUQP integrate test() FAILED" << std::endl;
        std::cerr << "Should be equal: " << std::endl << p << std::endl << dp << std::endl;
        std::cerr << "Should be equal: " << std::endl << q << std::endl << d3p << std::endl;
        exit(1);
    }

    std::vector<std::string> inputs;
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString());

    ExpressionTree resTree = p.convertToExpressionTree();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP integral() test: FAILED1" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString() + "$" + std::to_string(3));

    resTree = q.convertToExpressionTree();
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP integral() test: FAILED2" << std::endl;
        std::cerr << "Expected: " << q << " but got: " << retErr << std::endl;
        exit(1);
    }

    //zero
    p.zero();
    pTree = p.convertToExpressionTree();
    p.integrate(1);
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString());
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP integrate() test: FAILED3" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    }

    //constant
     p += RationalNumber(37, 13);
    pTree = p.convertToExpressionTree();
    p.integrate(1);
    resTree = p.convertToExpressionTree();
    inputs.clear();
    inputs.push_back(pTree.toMapleString());
    inputs.push_back(syms.toString());
    if(mapleTest->testProcReturn("int", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP integrate() test: FAILED4" << std::endl;
        std::cerr << "Expected: " << p << " but got: " << retErr << std::endl;
        exit(1);
    } 

    

    std::cerr << "SUQP integrate() test: PASSED" << std::endl;
}
void testEvaluate() {
   
    //test random
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    //test all vars set
    Symbol vars = p.variable();
    RationalNumber vals;
    SparseUnivariateRationalPolynomial r;
    
        r.randomPolynomialU( 1, coefBound, sparsity, includeNeg);         //get a random rat num
       
        vals=r.leadingCoefficient();
    

    RationalNumber res = p.evaluate( vals);

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

        varss << vars << "=" << vals<< ", ";
    
    varss << vars << "=" << vals<< "]";
    inputs.push_back(varss.str());

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testProcReturn("eval", inputs, resTree.toMapleString(), &retErr) == 0) {
        std::cerr << "SUQP evaluate() test: FAILED1" << std::endl;
        std::cerr << "Expected: " << res << " but got: " << retErr << std::endl;
        exit(1);
    }

    //test zero
    p.zero();

    vars = p.variable();
    //vals.clear();
   
        vals=(RationalNumber(rand() % coefBound));
    
    res = p.evaluate( vals);
    if (!res.isZero()) {
        std::cerr << "SUQP evaluate() test: FAILED3" << std::endl;
        std::cerr << "Expected 0 but got: " << res << std::endl;
        exit(1);
    }

    //test constant
    p.zero();
    mpq_class constTerm(rand() % coefBound);
    p += constTerm;

    vars = p.variable();

  
        vals=RationalNumber(rand() % coefBound);
    
    res = p.evaluate( vals);
    if (!res.isConstant() || !(res - constTerm).isZero()) {
        std::cerr << "SMQP evaluate() test: FAILED4" << std::endl;
        std::cerr << "Expected a constant: " << p << " but got: " << res << std::endl;
        exit(1);
    }
    
    //test constant with nvar = 0
    p = SparseUnivariateRationalPolynomial();
    p.zero();
    constTerm = mpq_class(rand() % coefBound);
    p += constTerm;

    vars = p.variable();

    res = p.evaluate( vals);
    if (!res.isConstant() || !(res - constTerm).isZero()) {
        std::cerr << "SUQP evaluate() test: FAILED5" << std::endl;
        std::cerr << "Expected a constant: " << p << " but got: " << res << std::endl;
        exit(1);
    }


    std::cerr << "SUQP evaluate() test: PASSED" << std::endl;
}


/* void testPseudoDivide() {
 
    SparseUnivariateRationalPolynomial c;
    c.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    SparseUnivariateRationalPolynomial b;
    unsigned long numBTerms = nterms;
    if (nterms > 10) {
        numBTerms = nterms/2;
    } else if (nterms > 3) {
        numBTerms = nterms - 2;
    }

    //Test same nvar for dividend and divisor
    b.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    while (b.isConstant()) {
        b.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    }

    SparseUnivariateRationalPolynomial quo;
    RationalNumber mult;
    SparseUnivariateRationalPolynomial rem = c.pseudoDivide(b, &quo, &mult);
   
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();
    
    
    ExpressionTree cTree = c.convertToExpressionTree();
    ExpressionTree bTree = b.convertToExpressionTree();
    ExpressionTree aTree = quo.convertToExpressionTree();
    ExpressionTree rTree = rem.convertToExpressionTree();
    ExpressionTree hTree = mult.convertToExpressionTree();
  

    Symbol var = b.variable();//leadingvariable

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
    //m vs    h
    algebList[6] = EvalMapleProc(kv, expandF, 1, algebList[6]);
    ALGEB comp3 = EvalMapleProc(kv, cmpF, 2, algebList[6], algebList[5]);

    M_BOOL compBool = MapleToM_BOOL(kv, comp);
    M_BOOL comp2Bool = MapleToM_BOOL(kv, comp2);
    M_BOOL comp3Bool = MapleToM_BOOL(kv, comp3);

    if (!compBool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[7]);
        std::cout << "SMQP Pseudo divide() test: FAILEDfirst" << std::endl;
        std::cout << "Got quotient:" << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    } 

    if (!comp3Bool) {
        std::string retErr = mapleTest->algebToString(kv, algebList[6]);
        std::cout << "SMQP Pesudo divide() test: FAILEDsecond" << std::endl;
        std::cout << "Got c multiplier:" << retErr << std::endl;
        std::cout << "Expected: " <<  mult  << std::endl;
        exit(1);
    }

    if (!comp2Bool) {
        std::string retErr = mapleTest->algebToString(kv, result);
        std::cout << "SMQP Pesudo divide() test: FAILEDthird" << std::endl;
        std::cout << "Got remainder:" << retErr << std::endl;
        std::cout << "Expected: " <<  rem << std::endl;
        exit(1);
    }
   
    //test divisor size 1

    numBTerms = 1;
    c.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    b.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    while (b.isConstant()) {
        b.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    }

    rem = c.pseudoDivide(b, &quo, &mult);

    cTree = c.convertToExpressionTree();
    bTree = b.convertToExpressionTree();
    aTree = quo.convertToExpressionTree();
    rTree = rem.convertToExpressionTree();
    hTree = mult.convertToExpressionTree();

    var = b.variable();

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
 */


void testSUQPShiftRight() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    
   
    //q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    SparseUnivariateRationalPolynomial q('x');


    int n=rand() % coefBound;
    q=q^(n);


    SparseUnivariateRationalPolynomial prod = p*q;


    p=p<<n;

    ExpressionTree pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree prodTree = prod.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift right test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

// SUQP zero
  p.zero();
  
   prod = p*q;
   p=p<<n;

   pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
   prodTree = prod.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= prodTree;

  
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift right test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
 // SUQP constant 

    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    prod = p * q;
 p=p<<n;

   pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
   prodTree = prod.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= prodTree;

  
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift right test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
    std::cout << "SUQP suqp multiplication test: PASSED" << std::endl;
}
void testSUQPShiftRightAssignment() {

    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);

    
   
    //q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    SparseUnivariateRationalPolynomial q('x');


    int n=rand() % coefBound;
    q=q^(n);


    SparseUnivariateRationalPolynomial prod = p*q;


    p<<=n;

    ExpressionTree pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree prodTree = prod.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift right  assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

// SUQP zero
  p.zero();
  
   prod = p*q;
    p<<=n;

   pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
   prodTree = prod.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= prodTree;

  
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift right  assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
 // SUQP constant 

    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    prod = p * q;
 p<<=n;

   pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
   prodTree = prod.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= prodTree;

  
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift right  assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
    std::cout << "SUQP shift right  assignment test: PASSED" << std::endl;
}

void testSUQPShiftLeft() {

 SparseUnivariateRationalPolynomial p;
 p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
   
    //q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    SparseUnivariateRationalPolynomial q('x'),quo,rem(p);

  

    int n=rand() % coefBound;
    q=q^(n);

quo=rem.monicDivide(q);

p>>=n;

    ExpressionTree pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTree = quo.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= quoTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift left assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    

    //      - SUQP that is zero
    p.zero();
    p>>=n;
    ExpressionTree pTreezero = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTreezero = quo.convertToExpressionTree();
    //pTree *= qTree;
    pTreezero -= quoTreezero;
 

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift left assignment test : FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

 // SUQP constant 

    p.zero();
    p += mpq_class(rand() % coefBound + 1);
    p>>=n;
 ExpressionTree pTreeconst = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTreeconst = quo.convertToExpressionTree();
    //pTree *= qTree;
    pTreeconst -= quoTreeconst;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift left assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

   
    std::cout << "SUQP suqp shift left test: PASSED" << std::endl;
}
void testSUQPShiftLeftAssignment() {

 SparseUnivariateRationalPolynomial p;
 p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
   
    //q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    SparseUnivariateRationalPolynomial q('x'),quo,rem(p);

  

    int n=rand() % coefBound;
    q=q^(n);

quo=rem.monicDivide(q);



p>>=n;

    ExpressionTree pTree = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTree = quo.convertToExpressionTree();
    //pTree *= qTree;
    pTree -= quoTree;
  
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift left test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    

    //      - SUQP that is zero
    p.zero();
p>>=n;

    ExpressionTree pTreezero = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTreezero = quo.convertToExpressionTree();
    //pTree *= qTree;
    pTreezero -= quoTreezero;

  

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift left assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

 // SUQP constant 

    p.zero();
p>>=n;

    p += mpq_class(rand() % coefBound + 1);
    ExpressionTree pTreeconst = p.convertToExpressionTree();
    //ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree quoTreeconst = quo.convertToExpressionTree();
    //pTree *= qTree;
    pTreeconst -= quoTreeconst;


    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP shift left assignment test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

   
    std::cout << "SUQP suqp shift left assignment  test: PASSED" << std::endl;
}
void testSUQPDivisionRationalNumber() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
mpq_class val = mpq_class(rand() % coefBound + 1);
RationalNumber q=RationalNumber(val);
   
    p *= p;
    p *= p;
    p *= p;

    SparseUnivariateRationalPolynomial quo = p/q;

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
        std::string retErr = mapleTest->algebToString(kv, algebList[3]);
        std::cout << "SUQP suqp division rational number test1: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }



 //      - SUQP that is zero
   // q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division rational number test2: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }
 
//      - SUQP that is a constant
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
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
        std::cout << "SUQP suqp division rational number  test3: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SUQP that is a constant and is nvar == 0

    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
   // q = SparseUnivariateRationalPolynomial();
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
        std::cout << "SUQP suqp division test rational number4 : FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp division rational number test: PASSED" << std::endl;
}
void testSUQPDivisionRationalNumberAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
  mpq_class val = mpq_class(rand() % coefBound + 1);
 RationalNumber q=RationalNumber(val);
   
   /*  p *= p;
    p *= p;
    p *= p; */
    ExpressionTree pTree = p.convertToExpressionTree();

   p/=q;

    SparseUnivariateRationalPolynomial quo =p;

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
        std::cout << "SUQP suqp division assignment rational number test1: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }



 //      - SUQP that is zero
   // q.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    p.zero();
    inputs.clear();

    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");
    
      p/=q;
      quo=p;
    
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
        std::cout << "SUQP suqp division assignment rational number test2: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }
 
//      - SUQP that is a constant
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");


    q.zero();
    q += mpq_class(rand() % coefBound + 1);

    p/=q; 
    quo=p;
     
    
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
        std::string retErr = mapleTest->algebToString(kv, algebList[3]);
        std::cout << "SUQP suqp division assignment rational number  test3: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

//      - SUQP that is a constant and is nvar == 0

    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
   // q = SparseUnivariateRationalPolynomial();
    q.zero();
    q += mpq_class(rand() % coefBound + 1);

    inputs.clear();
    inputs.push_back(p.convertToExpressionTree().toMapleString() + ":");
     p/=q; 
     quo=p;
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
        std::cout << "SUQP suqp division assignment test rational number4: FAILED" << std::endl;
        std::cout << "Got " << retErr << std::endl;
        std::cout << "Expected: " << quo << std::endl;
        exit(1);
    }

    std::cout << "SUQP suqp division rational number test: PASSED" << std::endl;
}
/* 
void testSquareFree() {
    q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
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
    std::cerr << "SUQP squareFree() test: PASSED" << std::endl;

}
 */
void testSUQPMultiplicationRationalNumber() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    mpq_class val = mpq_class(rand() % coefBound );
 RationalNumber q=RationalNumber(val);

    SparseUnivariateRationalPolynomial prod = p*q;

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
        std::cerr << "SUQP suqp multiplication rational number test: FAILED" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    
    //      - SUQP that is zero
    q.zero();
    prod= p*q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication rational number test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP that is a constant
    q.zero();
    q += mpq_class(rand() % coefBound + 1);
    prod = p * q;

    pTree = p.convertToExpressionTree();
    pTree *= q.convertToExpressionTree();
    pTree -= prod.convertToExpressionTree();

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication rational number test: FAILED" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
    std::cout << "SUQP suqp multiplication rational number  test: PASSED" << std::endl;
}
 void testSUQPMultiplicationRationalNumberAssignment() {
    SparseUnivariateRationalPolynomial p;
    p.randomPolynomialU(nterms, coefBound, sparsity, includeNeg);
    
    mpq_class val = mpq_class(rand() % coefBound );
    RationalNumber q=RationalNumber(val);

    SparseUnivariateRationalPolynomial prod = p;
    p*=q;
    ExpressionTree pTree = p.convertToExpressionTree();
    ExpressionTree qTree = q.convertToExpressionTree();
    ExpressionTree prodTree = prod.convertToExpressionTree();
    prodTree *= q.convertToExpressionTree();

    pTree -= prodTree;

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();

    std::string retErr;
    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication rational number asignment test: FAILED3" << std::endl; 
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }

    
    //      - SUQP that is zero
    p.zero();
    prod = p;

     p*=q;
     pTree = p.convertToExpressionTree();
     qTree = q.convertToExpressionTree();
     prodTree = prod.convertToExpressionTree();
    prodTree *= q.convertToExpressionTree();

    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication rational number assignment test: FAILED2" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }


    //      - SUQP that is a constant
    p.zero();
    p += mpq_class(rand() % coefBound );
        std::cerr << "p===" << p<< std::endl;

    prod = p;
        std::cerr << "prod===" << prod<< std::endl;

     p*=q;
        std::cerr << p*q<<"===" << p<< std::endl;

     pTree = p.convertToExpressionTree();
     qTree = q.convertToExpressionTree();
     prodTree = prod.convertToExpressionTree();
    prodTree *= q.convertToExpressionTree();

    pTree -= prodTree;

    mapleTest->testIfZero(pTree, &retErr);
    if (retErr != "") {
        std::cerr << "SUQP suqp multiplication rational number assignment test: FAILED2" << std::endl;
        std::cerr << ">>>> expected 0 but got: " << retErr << std::endl;
        exit(1);
    }
    std::cout << "SUQP suqp multiplication rational number assignment test: PASSED" << std::endl;
}
/////////////////////////////////////////////////////////////////////////////////////
int main(void)
{
//testIdentityConstructor();
//testCopyConstructor();
//testMoveConstructor();
 //testconvertToExpressionTree();
//testIsZero () ;
 //testZero();
//  testIsOne() ;
//testOne();
//testIsConstant();
//testOne();
//testCopyAssignment() ;
//testSUQPAdditionAssignment();
//testSUQPSubtraction();
//testSUQPSubtractionAssignment();
//testSUQPMultiplication();
//testSUQPMultiplicationAssignment();
//testSUQPDivision(); 
//testSMQPDivisionAssignment();
//testExponentiation();
//testExponentiationAssignment();
 //testLeadingCoefficient();
 //testCoefficient();
//testDerivative();
//testIntegrate();
//testEvaluate();
//testPseudoDivide();  ?? poly to string
 //testSUQPShiftRight();
// testSUQPShiftLeft();
 //testSUQPShiftRightAssignment();
//testSUQPShiftLeftAssignment();
//testSUQPDivisionRationalNumber();

 //testSUQPDivisionRationalNumberAssignment();

// testSUQPDivisionRationalNumberAssignment();


 
//testSUQPMultiplicationRationalNumber();
//testSUQPMultiplicationRationalNumberAssignment();

//exit(0);


   std::cerr<<std::endl;
   std::cerr<<std::endl;
   std::cerr<<std::endl;
   std::cerr<<std::endl;
   Symbol sym ='x';
   SparseUnivariateRationalPolynomial ss('y'),gg('y');


   SparseUnivariateRationalPolynomial q(sym);
   SparseUnivariateRationalPolynomial p;
   SparseUnivariateRationalPolynomial x('y');
   SparseUnivariateRationalPolynomial o; 
///////////////////////////////////////////////////////

ss=ss*ss; 
ss=(ss+1);
 
 ss=ss*2;
ss = ss + 5;
//////////////////////////////////////////////
gg=gg*gg*2;
// gg==(gg+1);


  gg=gg*ss;
      

   gg.poly=  sortPolynomial_AAU(gg.poly);

    
          gg.poly=  sortPolynomial_AAU(gg.poly);
   

SparseUnivariateRationalPolynomial pol;
   
 /////////////////////////////////////////////////////////////////////////////////////maple test of iedentity constructor

            ofstream myfile;
            myfile.open ("programtime1.txt");
         SparseUnivariateRationalPolynomial pol_g;
         SparseUnivariateRationalPolynomial pol_f;
         SparseUnivariateRationalPolynomial pol_a;
                  
         pol_g.setVariableName('y');
         pol_f.setVariableName('y');
         pol_a.setVariableName('y');
         long long unsigned int start;
         long long unsigned int start2;

         float elapsed;
         float Time[3];
         float total;

          float elapsed2;
         float Time2[3];
         float total2;
  std::cerr << "//////////////////////////////////////////////////" << std::endl;
  
        float sparsity=0.5; int includeNeg=0;
        for(int coef_size=1;coef_size<=11;coef_size=coef_size +1){
         unsigned long int coefBound=coef_size; 
          fprintf(stderr, "coefbound=%ld\n",coefBound);
           for (int deg_fg=0;deg_fg<=50;deg_fg=deg_fg + 5){
        for(int deg_a=0;deg_a<=22;deg_a=deg_a + 5){
         for(int i=0; i<1;i++)
         {pol_g.buildRandomPolyFromMax_seededU(deg_fg, coefBound,  sparsity,  includeNeg);
         pol_f.buildRandomPolyFromMax_seededU(deg_fg, coefBound,  sparsity,  includeNeg);
         pol_a.buildRandomPolyFromMax_seededU(deg_a, coefBound,  sparsity,  includeNeg);

        // std::cerr << "\npol_g ==" << pol_g<< std::endl;
       //  std::cerr << "\npol_f ==" << pol_f<< std::endl;
       //  std::cerr << "\npol_a ==                                       " << pol_a<< std::endl;
         gg=pol_g  * pol_a;
         pol=pol_f  * pol_a;
         cerr << "f = " << pol_f << endl;
         cerr << "g = " << pol_g << endl;
         cerr << "a = " << pol_a << endl;
         cerr << "f*a =" << pol << endl;
         cerr << "g*a = " << gg << endl;
        
                  ///////////////////  Modular method       

         startTimer(&start);
         //x = pol.ModularGCD_AAU(gg);
         x = pol.ModularGCD_AAU(gg);

         stopTimer(&start,&elapsed);
         cerr << "gcdMod(pol,gg) = " << x << endl << endl;
         cerr << endl << "Modular method time: " << elapsed << endl;
        // std::cerr << "\n\n\n\n x ==" << x<< std::endl;

         Time[i]=elapsed; 
 
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////// euclidean Algorithems 


        startTimer(&start2);
         x = pol.gcd(gg);
         stopTimer(&start2,&elapsed2);
         cerr << "gcd(pol,gg) = " << x << endl << endl;
         cerr << endl << "Euclidean GCD time: " << elapsed2 << endl;
        // std::cerr << "\n\n\n\n x ==" << x<< std::endl;

         Time2[i]=elapsed2;





///////////////////////////////////////////////////////////////////////////////////////////////////////

         }
           total=Time[0]+Time[1]+Time[2];
           total=total/3;
         ////////////////////////////////////////////////////// Modular Time 
            total2=Time2[0]+Time2[1]+Time2[2];
           total2=total2/3;



           /////////////////////////////////////////////// Time

            myfile <<coef_size<<"\t \t  "<< deg_fg<<"\t \t "<<deg_a<< "\t \t "<<total<<"\t \t" <<total2<<"\n";
                          
         }
     }
 }
     myfile.close();


std::fprintf(stderr,"\n++++++++++++++++++++++++++++++++++++++\n");
	return 0;
}

