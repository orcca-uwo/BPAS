    
#include <gmpxx.h>

#include "../../include/IntegerPolynomial/DUZP_Support.h"
#include "../../include/IntegerPolynomial/SMZP_Support_Test.h"

#include "../../include/ExpressionTree/ExpressionTree.hpp"
#include "../../tests/MapleTestTool/MapleTestTool.hpp"

#include "../../include/Utils/Unix_Timer.h"


polysize_t maxDeg = 10;
int coefBound = 5;
float sparsity = 0.9;
int includeNeg = 0;

int main(int argc, char** argv) {

    if (argc > 1 && atol(argv[1]) > 0) { 
        maxDeg = atol(argv[1]);
    }
    if (argc > 2 && atol(argv[2]) > 0) {
        coefBound = atol(argv[2]);
    }
    if (argc > 3 && atof(argv[3])) {
        sparsity = atof(argv[3]);
    }
    if (argc > 4 && atoi(argv[4]) >= 0) {
        includeNeg = atoi(argv[4]);
    }

	// DUZP_t* p = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
	// printPoly_DUZP(p, "x");

 //    Prime_ptr* pptr = (Prime_ptr*) malloc(sizeof(Prime_ptr));
 //    *pptr = prime64_ptr[0];
 //    duspoly_t* pp = convertToPrimeField_DUZP(p, pptr);
 //    printPolynomialConvertingOut_spX(pp, pptr);

 //    DUZP_t* ppp = deepCopyPolynomial_DUZPFromspX(pp, pptr);
 //    printPoly_DUZP(ppp, "x");

 //    DUZP_t* pppp = deepCopyPolynomial_DUZP(ppp);
 //    printPoly_DUZP(pppp, "x");

 //    DUZP_t* prim_part = primitivePart_DUZP(pppp);
 //    printPoly_DUZP(prim_part, "x");

 //    mpz_t cont;
 //    mpz_init(cont);
 //    content_DUZP(p, cont);
 //    gmp_fprintf(stderr, "Content: %Zd\n", cont);

 //    mpz_clear(cont);
 //    mpz_init(cont);

 //    DUZP_t* prim_part2 = primitivePartAndContent_DUZP( pppp, cont);
 //    printPoly_DUZP(prim_part2, "x");
 //    gmp_fprintf(stderr, "Content: %Zd\n", cont);

 //    primitivePart_DUZP_inp(pppp);
 //    printPoly_DUZP(pppp, "x");

    DUZP_t* p = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    // mpz_set_si(p->coefs[0], 123l);
    DUZP_t* q = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    // mpz_set_si(q->coefs[0], 17l);a
    DUZP_t* s = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    DUZP_t* prod = multiplyPolynomials_DUZP(p,q);
    DUZP_t* prod2 = multiplyPolynomials_DUZP(p,s);
    // printPoly_DUZP_maple(p, "x", "p");
    // printPoly_DUZP_maple(q, "x", "q");
    // printPoly_DUZP_maple(s, "x", "s");
    printPoly_DUZP_maple(prod, "x", "p");
    printPoly_DUZP_maple(prod2, "x", "q");

    // DUZP_t* quo = NULL;
    // int divisible = divideTest_DUZP(prod, p, &quo);
    // printPoly_DUZP_maple(quo, "x", "tempquo");

    // fprintf(stderr, "divisible: %d\n", divisible);

    // if (!divisible || !isEqual_DUZP(q, quo) ) {
    //     fprintf(stderr, "divideTest_DUZP FAIELD!\n");
    //     exit(1);
    // }

    timer_id timeID1 = start_timer();
    DUZP_t* g = GCD_DUZP(prod, prod2);
    timer_time time1 = elapsed_time(&timeID1);
    double time1D = (time1.tv_sec + ((double)time1.tv_usec / 1000000));
    fprintf(stderr, "time: %0.5f\n", time1D);

    // printPoly_DUZP_maple(g, "x", "g");
    // if (mpz_sgn(p->coefs[p->lt]) < 0) {
    //     mpz_t z;
    //     mpz_init_set_si(z, -1l);
    //     multiplyByInteger_DUZP_inp(p, z);
    //     mpz_clear(z);
    // }
    // if (!isEqual_DUZP(g, p)) {
    //     fprintf(stderr, "GCD_DUZP FAILED!\n");
    //     exit(1);
    // }

    // DUZP_t* g2 = GCD_DUZP(p, q);
    // if (g2->lt != 0) {
    //     fprintf(stderr, "GCD_DUZP FAILED!\n");
    //     exit(1);
    // }



}
