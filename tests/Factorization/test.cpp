#include "../../include/bpas.h"
#include <iomanip>
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include "../MapleTestTool/MapleTestTool.hpp"

#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include "../../include/RationalNumberPolynomial/SMQP_Support-AA.h"

long nvar = 3;
long numTerms = 10;
long coefBound = 5ul;
degree_t sparsity = 10;
int includeNeg = 1;

int main() {

    AltArr_t*** factors;
    int **exponet;
   

    AltArr_t *minimalPoly = makePolynomial_AA(3, 2);// 3 term and 2 variable
    mpq_init(minimalPoly->elems[0].coef);
    mpq_set_si(minimalPoly->elems[0].coef, -3, 1);
   // partialDegreeTerm_AA(minimalPoly, 1, 1); return second degree of the second term
    

   //printAA(minimalPoly);
const degree_t degsList[2] = {1, 1};
setDegrees_AA_inp(minimalPoly, 0,degsList, 2);// first term ,nvar=2

    //    void Factorization_AA ( minimalPoly,  minimalPoly,  factors, exponet );*/
 
minimalPoly->size=1;
const char* syms[]={"x","y"};
//printPoly_AA(stderr, minimalPoly,syms,2);

//Factorization_AA ( minimalPoly,  minimalPoly,  factors, exponet );
//////////////////////////////////////
SparseMultivariateRationalPolynomial polyTotalcontent;
polyTotalcontent=SparseMultivariateRationalPolynomial();
SparseMultivariateRationalPolynomial polyfirst;
polyfirst=SparseMultivariateRationalPolynomial();

polyfirst=SMQP("x^4-3");//"x^4-3"   x^2-2
std::cerr << "\n polyfirst=" <<polyfirst<<"\n";
SparseMultivariateRationalPolynomial Irrudipoly;
Irrudipoly=SMQP("(z-1)^2*a^2"); // Irrudipoly=SMQP("");z^4*a^0-2
std::cerr<<"\n Irrudipoly="<<Irrudipoly<<"\n";

  Factors<SparseMultivariateRationalPolynomial> factor = Irrudipoly.AlgebricFactorization (polyfirst,&polyTotalcontent);
std::cerr<<"\n polyTotalcontent="<<polyTotalcontent<<"\n";
    return 0;
    
}


