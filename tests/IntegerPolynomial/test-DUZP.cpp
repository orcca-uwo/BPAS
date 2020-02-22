    
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

    if (argc > 1 && atol(argv[1]) >= 0) { 
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


/****************
 * Test conversion
 ****************/

    // DUZP_t* p = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
//     AltArr_t* qAA = convertToAltArr_DUZP(p);
//     AltArrZ_t* zAA = convertToAltArrZ_DUZP(p);
// //    printPoly_DUZP(p, "x");
//     char* syms[] = {"x"};
// //    printPoly_AA(stderr, qAA, syms, 1);
// //    fprintf(stderr, "\n");
// //    printPoly_AAZ(stderr, zAA, syms, 1);
// //    fprintf(stderr, "\n");
// //    fprintf(stderr, "\n");

//     DUZP_t* pq = convertFromAltArr_DUZP(qAA);
//     DUZP_t* pz = convertFromAltArrZ_DUZP(zAA);
// //    printPoly_DUZP(pq, "x");
// //    printPoly_DUZP(pz, "x");

//     if(!isEqual_DUZP(p, pq)) {
// 	fprintf(stderr, "Failed to convert to/from Q!\n");
// 	printPoly_DUZP(p, "x");
// 	printPoly_AA(stderr, qAA, syms, 1);
// 	fprintf(stderr, "\n");
// 	printPoly_DUZP(pq, "x");
//     	return 1;
//     }

//     if(!isEqual_DUZP(p, pz)) {
//         fprintf(stderr, "Failed to convert to/from Z!\n");
//         printPoly_DUZP(p, "x");
//         printPoly_AAZ(stderr, zAA, syms, 1);
//         fprintf(stderr, "\n");
//         printPoly_DUZP(pz, "x");
//         return 1;
//     }
	

//     return 0;


/****************
 * Test AAZ univar GCD via DUZP
 ****************/

    DUZP_t* p = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    mpz_set_si(p->coefs[0], 123l); //ensure polys have constant term so GCD check can exactly equal p
    DUZP_t* q = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    mpz_set_si(q->coefs[0], 17l); //ensure polys have constant term so GCD check can exactly equal p
    DUZP_t* s = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
        
    mpz_set_si(p->coefs[p->lt], 1l);
    mpz_set_si(q->coefs[q->lt], 1l);
    mpz_set_si(s->coefs[s->lt], 1l);    

    p = multiplyPolynomials_DUZP(p,s);
    q = multiplyPolynomials_DUZP(q,s);
    AltArrZ_t* pAAZ = convertToAltArrZ_DUZP(p);
    AltArrZ_t* qAAZ = convertToAltArrZ_DUZP(q);

    timer_id id = start_timer();
    AltArrZ_t* gAAZ = univariateGCD_AAZ(pAAZ, qAAZ);
    timer_time elapsed1 = elapsed_time(&id);
    double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));

    DUZP_t* g = convertFromAltArrZ_DUZP(gAAZ);
    if (!isEqual_DUZP(g, s)) {
        fprintf(stderr, "univariateGCD_AAZ FAILED!\n");
        exit(1);
    }

    fprintf(stderr, "univariateGCD_AAZ PASSED!\n");

    AltArr_t* pAA = convertToAltArr_DUZP(p);
    AltArr_t* qAA = convertToAltArr_DUZP(q);
    AltArr_t* gAA = univariateGCD_AA(pAA, qAA);

    // char* syms[] = {"x"};
    // fprintf(stderr, "\np:\n");
    // printPoly_AA(stderr, pAA, syms, 1);
    // fprintf(stderr, "\nq:\n");
    // printPoly_AA(stderr, qAA, syms, 1);
    // fprintf(stderr, "\ng:\n");
    // printPoly_AA(stderr, gAA, syms, 1);
    // fprintf(stderr, "\n\n");

    g = convertFromAltArr_DUZP(gAA);
    // printPoly_DUZP_maple(g, "x", "g");
    // printPoly_DUZP_maple(s, "x", "s");
    if (!isEqual_DUZP(g, s)) {
        fprintf(stderr, "univariateGCD_AA FAILED!\n");
        exit(1);
    }

    fprintf(stderr, "univariateGCD_AA PASSED!\n");

    std::cerr << "gcd time: " << time << std::endl;

    return 0;

/****************
 * Test DUZP GCD
 ****************/

    // DUZP_t* p = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    // mpz_set_si(p->coefs[0], 123l); //ensure polys have constant term so GCD check can exactly equal p
    // DUZP_t* q = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    // mpz_set_si(q->coefs[0], 17l);
    // DUZP_t* s = buildRandomPoly_DUZP(maxDeg, coefBound, sparsity, includeNeg);
    // DUZP_t* prod = multiplyPolynomials_DUZP(p,q);
    // DUZP_t* prod2 = multiplyPolynomials_DUZP(p,s);
    // printPoly_DUZP_maple(p, "x", "p");
    // printPoly_DUZP_maple(q, "x", "q");
    // printPoly_DUZP_maple(s, "x", "s");
    // printPoly_DUZP_maple(prod, "x", "prod");
    // printPoly_DUZP_maple(prod2, "x", "prod2");

    // // DUZP_t* quo = NULL;
    // // int divisible = divideTest_DUZP(prod, p, &quo);
    // // printPoly_DUZP_maple(quo, "x", "tempquo");

    // // fprintf(stderr, "divisible: %d\n", divisible);

    // // if (!divisible || !isEqual_DUZP(q, quo) ) {
    // //     fprintf(stderr, "divideTest_DUZP FAIELD!\n");
    // //     exit(1);
    // // }

    // DUZP_t* g = GCD_DUZP(prod, prod2);
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
