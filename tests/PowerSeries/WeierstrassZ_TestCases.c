

#include "WeierstrassZ_TestCases.h"
#include "../../include/Parser/bpas_parser.h"



/**
 * f : a upops
 * bound : the maximum integer for which homogPart_PSZ gets called
 * sym : a list of variables
 */



void TestingWeierstrassPreparation(UpopsZ_t* f, int bound, const char** sym) {

    UpopsZ_t* p_out = NULL;
    UpopsZ_t* alpha_out = NULL;
    weierstrassPreparation_UPOPSZ(f, &p_out, &alpha_out);
    UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ(p_out, alpha_out);
    UpopsZ_t* sub =  subUnivariatePolyOverPowerSeries_UPOPSZ(prod,f);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "TestingWeierstrassPreparation\n");
#endif
    for (int d = 1; d <= bound; ++d) {
        homogPart_PSZ(d, p_out->data[0]);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
        fprintf(stderr, "In degree: %d\n ", d);
            fprintf(stderr, "\n");
            fprintf(stderr, "current f = \n");
            print_UPOPSZ(stderr,  f, sym);
            fprintf(stderr, "\n");
        if (d % 2 == 0) {
            fprintf(stderr, "polyPart Alpha of degree d = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPSZ(stderr,  alpha_out, sym);
            fprintf(stderr, "\n");
        } else {
            fprintf(stderr, "\n");
            fprintf(stderr, "polyPart P of degree d = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPSZ(stderr,  p_out, sym);
            fprintf(stderr, "\n");
        }
#endif

        // fprintf(stderr, "\n\nComputing poly part of prod:\n");
        // PolyZ_ptr polyPartProd = polynomialPart_UPOPSZ(prod, d);
        // fprintf(stderr, "\n Done computing poly part of prod:\n");
        // printPoly_AAZ(stderr, polyPartProd, sym, polyPartProd == NULL ? 0 : polyPartProd->nvar);
        // fprintf(stderr, "\n" );
        // fprintf(stderr, "\n\nComputing poly part of sub:\n");

        PolyZ_ptr polyPart = polynomialPart_UPOPSZ(d, sub);
        if (!isZero_AAZ(polyPart)) {
            printPoly_AAZ(stderr, polyPart, sym, polyPart == NULL ? 0 : polyPart->nvar);
            fprintf(stderr, "ERROR IN WEIERSTRASS PREPARATION!\n");
            exit(1);
        }
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
        if (d % 2 != 0) {
            fprintf(stderr, "\n");
            fprintf(stderr, "polyPart Alpha of degree d = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPSZ(stderr,  alpha_out, sym);
            fprintf(stderr, "\n");

        } else {
            fprintf(stderr, "\n");
            fprintf(stderr, "polyPart P of degree d  = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPSZ(stderr,  p_out, sym);
            fprintf(stderr, "\n");

        }
        fprintf(stderr, "\n");
#endif
    }

    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);
    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
}





/**
 *
 */
void test1(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArrZ_t* a = generate_altarrZ_var_defined("1+ z", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "****************The Input Upops****************:\n");
	printPoly_AAZ (stderr,a, sym, 3);
	fprintf(stderr, "\n");
#endif
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	UpopsZ_t* p_out = NULL;
	UpopsZ_t* alpha_out = NULL;
	weierstrassPreparation_UPOPSZ(upop, &p_out, &alpha_out);
    updateToDeg_UPOPSZ(7, upop);
	PowerSeriesZ_t* b = alpha_out->data[0];
	homogPart_PSZ(7, b);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "P is:\n");
	print_UPOPSZ(stderr,  p_out, sym);
	fprintf(stderr, "\n");
	fprintf(stderr, "alpha is:\n");
	print_UPOPSZ(stderr,  alpha_out, sym);
    fprintf(stderr, "\n");
#endif

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 1: \t\t\t\t PASSED \n");

}


/**
 *
 */
void test2(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("x+ y*z", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "****************The Input Upops****************:\n");
    printPoly_AAZ (stderr,a, sym, 3);
    fprintf(stderr, "\n");
#endif
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	UpopsZ_t* p_out = NULL;
	UpopsZ_t* alpha_out = NULL;
	weierstrassPreparation_UPOPSZ(upop, &p_out, &alpha_out);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);


}

/**
 *
 */
void test3(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArrZ_t* a = generate_altarrZ_var_defined("y+x^4*z+x*z^2", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "****************The Input Upops****************:\n");
	printPoly_AAZ (stderr,a, sym, 3);
	fprintf(stderr, "\n");
#endif
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	UpopsZ_t* p_out = NULL;
	UpopsZ_t* alpha_out = NULL;
	weierstrassPreparation_UPOPSZ(upop, &p_out, &alpha_out);

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	printPoly_AAZ (stderr, p_out->data[0]->polys[0], sym, 3);
    fprintf(stderr, "\n");
#endif

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

}


/**
 *
 */
void test4(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArrZ_t* a = generate_altarrZ_var_defined("y+z+x*z^2 ", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "****************The Input Upops****************:\n");
	printPoly_AAZ (stderr,a, sym, 3);
	fprintf(stderr, "\n");
#endif
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 4: \t\t\t\t PASSED \n");
}


/**
 * results in rational numbers
 */
// void test5(int testDeg)
// {
// 	const char* sym[] = {"z", "y", "x"};
// 	AltArrZ_t* a = generate_altarrZ_var_defined("y + 2*z + (x + 1 )*z^2", sym, 3);
// 	// fprintf(stderr, "****************The Input Upops****************:\n");
// 	// printPoly_AAZ (stderr,a, sym, 3);
// 	// fprintf(stderr, "\n");
// 	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
// 	TestingWeierstrassPreparation(upop, testDeg, sym);

//     destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
//     freePolynomial_AAZ(a);

//     fprintf(stderr, "Weierstrass Test 5: \t\t\t\t PASSED \n");
// }




/**
 * Because of this example I had to define c's before F's
 */
void test6(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArrZ_t* a = generate_altarrZ_var_defined("x^2 + y^2 + z+ (y + 1)*z^2 + z^3", sym, 3);
	// fprintf(stderr, "****************The Input Upops****************:\n");
	// printPoly_AAZ (stderr,a, sym, 3);
	// fprintf(stderr, "\n");
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 6: \t\t\t\t PASSED \n");
}

/**
 *
*/
void test7(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("x^2 + y^2 + z + z^2", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AAZ (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 7: \t\t\t\t PASSED \n");
}

/**
 * Maple example
 */
void test8(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("x + z", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AAZ (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 8: \t\t\t\t PASSED \n");
}

/**
 * Maple example
 */
void test9(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("1 + x*z", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AAZ (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 9: \t\t\t\t PASSED \n");
}

/**
 * Maple example
 */
void test10(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("y+z+ (x+1)*z^2", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AAZ (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 10: \t\t\t\t PASSED \n");
}



/**
 *
 */
void test11(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("y +z^2 +z^3", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AAZ (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 11: \t\t\t\t PASSED \n");
}



/**
 *
 */
void test12(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArrZ_t* a = generate_altarrZ_var_defined("x^2 + y^2 + (y + 1)*z^2 + z^3", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AAZ (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(a);

    fprintf(stderr, "Weierstrass Test 12: \t\t\t\t PASSED \n");
}



/**
 * Maple example : limited resources on my machine
 */
void test13(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1-y", sym, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("x", sym, 3);
	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
	PowerSeriesZ_t* ps_array[2];
	ps_array[0] = ps_3;
	ps_array[1] = quo_1;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 13: \t\t\t\t PASSED \n");
}


/**
 *
 */
void test14(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};

	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("x", sym, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y-x^2", sym, 3);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* ps_array[2];
	ps_array[0] =  ps_2;
	ps_array[1] = ps_3;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);

    fprintf(stderr, "Weierstrass Test 14: \t\t\t\t PASSED \n");
}



/**
 * Memory consumption
 */
void test15(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("x*y+x^3", sym, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y^2+x^2", sym, 3);
	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_3);
	PowerSeriesZ_t* ps_array[2];
	ps_array[0] = ps_2;
	ps_array[1] = quo_1;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 15: \t\t\t\t PASSED \n");
}

/**
 * Memory consumption
 */
void test16(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("x^2+y^2", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y^2+2*x^3", sym, 3);
    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_3);
    PowerSeriesZ_t* ps_array[4];
    ps_array[0] = ps_2;
    ps_array[1] = zeroPowerSeries_PSZ();
    ps_array[2] = onePowerSeries_PSZ(3);
    ps_array[3] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,4);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_array[1]);
    destroyPowerSeries_PSZ(ps_array[2]);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 16: \t\t\t\t PASSED \n");

}



 /**
 * It seems Maple is faster
 */
void test17(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1-x", sym, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y-x^2", sym, 3);
	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
	//PowerSeriesZ_t* sub_1 = subPowerSeries_PSZ(quo_1,onePowerSeries_PSZ());
	homogPart_PSZ(2, quo_1);
	PowerSeriesZ_t* ps_array[2];
	ps_array[0] = ps_3;
	ps_array[1] = quo_1;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);

	// print_UPOPSZ(stderr, upop , sym);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 17: \t\t\t\t PASSED \n");


}


/**
 * Maple example
 */
void test18(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1-x", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y", sym, 3);
    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
    homogPart_PSZ(6, quo_1);
    PowerSeriesZ_t* quo_2 = dividePowerSeries_PSZ(ps_1, ps_3);
    homogPart_PSZ(7, quo_2);
    PowerSeriesZ_t* ps_array[2];
    ps_array[0] = quo_1;
    ps_array[1] = quo_2;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);

    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(quo_1);
    destroyPowerSeries_PSZ(quo_2);

    fprintf(stderr, "Weierstrass Test 18: \t\t\t\t PASSED \n");
}


/**
 * Paper example1
 */
void test19(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1+x+y", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y", sym, 3);
    PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("x", sym, 3);

    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);

    PowerSeriesZ_t* ps_array[4];
    ps_array[0] = ps_4;
    ps_array[1] = ps_3;
    ps_array[2] = ps_1;
    ps_array[3] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,4);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);
    freePolynomial_AAZ(poly_4);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(ps_4);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 19: \t\t\t\t PASSED \n");
}

/**
 * Paper example2
 */
void test20(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1+x+y", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y", sym, 3);
    PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("x", sym, 3);

    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);

    PowerSeriesZ_t* ps_array[4];
    ps_array[0] = ps_4;
    ps_array[1] = ps_1;
    ps_array[2] = ps_3;
    ps_array[3] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,4);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);
    freePolynomial_AAZ(poly_4);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(ps_4);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 20: \t\t\t\t PASSED \n");
}


/**
 * Paper example2
 */
void test21(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1+x+y", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y", sym, 3);
    PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("y+1", sym, 3);
    PolyZ_ptr  poly_5 = generate_altarrZ_var_defined("x", sym, 3);

    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
    PowerSeriesZ_t* ps_5 = convertPolyToPowerSeries_PSZ(poly_5);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);

    PowerSeriesZ_t* ps_array[5];
    ps_array[0] = ps_5;
    ps_array[1] = ps_4;
    ps_array[2] = zeroPowerSeries_PSZ();
    ps_array[3] = ps_3;
    ps_array[4] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,5);
    TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);
    freePolynomial_AAZ(poly_4);
    freePolynomial_AAZ(poly_5);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_array[2]);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(ps_4);
    destroyPowerSeries_PSZ(ps_5);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 21: \t\t\t\t PASSED \n");
}


void test22(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1+x+y", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y", sym, 3);
    // PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("y+1", sym, 3);
    PolyZ_ptr  poly_5 = generate_altarrZ_var_defined("x", sym, 3);
    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    // PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
    PowerSeriesZ_t* ps_5 = convertPolyToPowerSeries_PSZ(poly_5);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
    PowerSeriesZ_t* ps_array[5];
    ps_array[0] = ps_5;
    ps_array[1] = ps_3;
    ps_array[2] = onePowerSeries_PSZ(3);
    ps_array[3] = onePowerSeries_PSZ(3);
    ps_array[4] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,5);
    TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);
    freePolynomial_AAZ(poly_5);

    destroyPowerSeries_PSZ(ps_1);
    destroyPowerSeries_PSZ(ps_array[2]);
    destroyPowerSeries_PSZ(ps_array[3]);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(ps_5);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 22: \t\t\t\t PASSED \n");

}


void test23(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1+x+y", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y", sym, 3);
    // PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("y+1", sym, 3);
    PolyZ_ptr  poly_5 = generate_altarrZ_var_defined("x", sym, 3);
    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    // PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
    PowerSeriesZ_t* ps_5 = convertPolyToPowerSeries_PSZ(poly_5);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
    PowerSeriesZ_t* ps_array[5];
    ps_array[0] = ps_5;
    ps_array[1] = ps_3;
    ps_array[2] = ps_3;//eroPowerSeries_PSZ();
    ps_array[3] = ps_3;//zeroPowerSeries_PSZ();
    ps_array[4] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,5);
    TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);
    freePolynomial_AAZ(poly_5);

    destroyPowerSeries_PSZ(ps_1);
    //destroyPowerSeries_PSZ(ps_array[2]);
    //destroyPowerSeries_PSZ(ps_array[3]);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(ps_5);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 23: \t\t\t\t PASSED \n");

}


void test24(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", sym, 3);
    PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1+x+y", sym, 3);
    PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("y", sym, 3);
    // PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("y+1", sym, 3);
    PolyZ_ptr  poly_5 = generate_altarrZ_var_defined("x", sym, 3);
    PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
    PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
    PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
    // PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
    PowerSeriesZ_t* ps_5 = convertPolyToPowerSeries_PSZ(poly_5);
    PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
    PowerSeriesZ_t* ps_array[5];
    ps_array[0] = ps_5;
    ps_array[1] = ps_3;
    ps_array[2] = ps_3; //onePowerSeries_PSZ(3);
    ps_array[3] = onePowerSeries_PSZ(3);//onePowerSeries_PSZ(3);
    ps_array[4] = quo_1;
    UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,5);
    TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);
    freePolynomial_AAZ(poly_5);

    destroyPowerSeries_PSZ(ps_1);
    //destroyPowerSeries_PSZ(ps_array[2]);
    destroyPowerSeries_PSZ(ps_array[3]);
    destroyPowerSeries_PSZ(ps_2);
    destroyPowerSeries_PSZ(ps_3);
    destroyPowerSeries_PSZ(ps_5);
    destroyPowerSeries_PSZ(quo_1);

    fprintf(stderr, "Weierstrass Test 24: \t\t\t\t PASSED \n");

}

