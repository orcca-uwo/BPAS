

#include "Weierstrass_TestCases.h"
#include "../../include/Parser/bpas_parser.h"


/**
 * f : a upops
 * bound : the maximum integer for which homogPart_PS gets called
 * sym : a list of variables
 */

void TestingWeierstrassPreparation(Upops_t* f, int bound, const char** sym) {

    Upops_t* p_out = NULL;
    Upops_t* alpha_out = NULL;
    weierstrassPreparation_UPOPS(f, &p_out, &alpha_out);
    Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS(p_out, alpha_out);
    Upops_t* sub =  subUnivariatePolyOverPowerSeries_UPOPS(prod,f);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "TestingWeierstrassPreparation\n");
#endif
    for (int d = 1; d <= bound; ++d) {
        homogPart_PS(d, p_out->data[0]);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
        fprintf(stderr, "In degree: %d\n ", d);
            fprintf(stderr, "\n");
            fprintf(stderr, "current f = \n");
            print_UPOPS(stderr,  f, sym);
            fprintf(stderr, "\n");
        if (d % 2 == 0) {
            fprintf(stderr, "polyPart Alpha of degree d = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPS(stderr,  alpha_out, sym);
            fprintf(stderr, "\n");
        } else {
            fprintf(stderr, "\n");
            fprintf(stderr, "polyPart P of degree d = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPS(stderr,  p_out, sym);
            fprintf(stderr, "\n");
        }
#endif

        // fprintf(stderr, "\n\nComputing poly part of prod:\n");
        // Poly_ptr polyPartProd = polynomialPart_UPOPS(prod, d);
        // fprintf(stderr, "\n Done computing poly part of prod:\n");
        // printPoly_AA(stderr, polyPartProd, sym, polyPartProd == NULL ? 0 : polyPartProd->nvar);
        // fprintf(stderr, "\n" );
        // fprintf(stderr, "\n\nComputing poly part of sub:\n");

        Poly_ptr polyPart = polynomialPart_UPOPS(d, sub);
        if (!isZero_AA(polyPart)) {
            printPoly_AA(stderr, polyPart, sym, polyPart == NULL ? 0 : polyPart->nvar);
            fprintf(stderr, "ERROR IN WEIERSTRASS PREPARATION!\n");
            exit(1);
        }
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
        if (d % 2 != 0) {
            fprintf(stderr, "\n");
            fprintf(stderr, "polyPart Alpha of degree d = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPS(stderr,  alpha_out, sym);
            fprintf(stderr, "\n");

        } else {
            fprintf(stderr, "\n");
            fprintf(stderr, "polyPart P of degree d  = %d\n", d);
            fprintf(stderr, "\n");
            print_UPOPS(stderr,  p_out, sym);
            fprintf(stderr, "\n");

        }
        fprintf(stderr, "\n");
#endif
    }

    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
    // fprintf(stderr, "p->refcount: %d, alpha->refCount: %d\n", p_out->refCount, alpha_out->refCount);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
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
	AltArr_t* a = generate_altarr_var_defined("1+ z", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "****************The Input Upops****************:\n");
	printPoly_AA (stderr,a, sym, 3);
	fprintf(stderr, "\n");
#endif
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Upops_t* p_out = NULL;
	Upops_t* alpha_out = NULL;
	weierstrassPreparation_UPOPS(upop, &p_out, &alpha_out);
    updateToDeg_UPOPS(7, upop);
	PowerSeries_t* b = alpha_out->data[0];
	homogPart_PS(7, b);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "P is:\n");
	print_UPOPS(stderr,  p_out, sym);
	fprintf(stderr, "\n");
	fprintf(stderr, "alpha is:\n");
	print_UPOPS(stderr,  alpha_out, sym);
    fprintf(stderr, "\n");
#endif

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 1: \t\t\t\t PASSED \n");

}


/**
 *
 */
void test2(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("x+ y*z", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "****************The Input Upops****************:\n");
    printPoly_AA (stderr,a, sym, 3);
    fprintf(stderr, "\n");
#endif
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Upops_t* p_out = NULL;
	Upops_t* alpha_out = NULL;
	weierstrassPreparation_UPOPS(upop, &p_out, &alpha_out);


    destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);


}

/**
 *
 */
void test3(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y+x^4*z+x*z^2", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "****************The Input Upops****************:\n");
	printPoly_AA (stderr,a, sym, 3);
	fprintf(stderr, "\n");
#endif
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Upops_t* p_out = NULL;
	Upops_t* alpha_out = NULL;
	weierstrassPreparation_UPOPS(upop, &p_out, &alpha_out);

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	printPoly_AA (stderr, p_out->data[0]->polys[0], sym, 3);
    fprintf(stderr, "\n");
#endif

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(p_out);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

}


/**
 *
 */
void test4(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y+z+x*z^2 ", sym, 3);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "****************The Input Upops****************:\n");
	printPoly_AA (stderr,a, sym, 3);
	fprintf(stderr, "\n");
#endif
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 4: \t\t\t\t PASSED \n");
}


/**
 *
 */
void test5(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y + 2*z + (x + 1 )*z^2", sym, 3);
	// fprintf(stderr, "****************The Input Upops****************:\n");
	// printPoly_AA (stderr,a, sym, 3);
	// fprintf(stderr, "\n");
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 5: \t\t\t\t PASSED \n");
}




/**
 * Because of this example I had to define c's before F's
 */
void test6(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("x^2 + y^2 + z+ (y + 1)*z^2 + z^3", sym, 3);
	// fprintf(stderr, "****************The Input Upops****************:\n");
	// printPoly_AA (stderr,a, sym, 3);
	// fprintf(stderr, "\n");
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 6: \t\t\t\t PASSED \n");
}

/**
 *
*/
void test7(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("x^2 + y^2 + z + z^2", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AA (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 7: \t\t\t\t PASSED \n");
}

/**
 * Maple example
 */
void test8(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("x + z", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AA (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 8: \t\t\t\t PASSED \n");
}

/**
 * Maple example
 */
void test9(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("1 + x*z", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AA (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 9: \t\t\t\t PASSED \n");
}

/**
 * Maple example
 */
void test10(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("y+z+ (x+1)*z^2", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AA (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 10: \t\t\t\t PASSED \n");
}



/**
 *
 */
void test11(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("y +z^2 +z^3", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AA (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 11: \t\t\t\t PASSED \n");
}



/**
 *
 */
void test12(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    AltArr_t* a = generate_altarr_var_defined("x^2 + y^2 + (y + 1)*z^2 + z^3", sym, 3);
    // fprintf(stderr, "****************The Input Upops****************:\n");
    // printPoly_AA (stderr,a, sym, 3);
    // fprintf(stderr, "\n");
    Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(a);

    fprintf(stderr, "Weierstrass Test 12: \t\t\t\t PASSED \n");
}



/**
 * Maple example : limited resources on my machine
 */
void test13(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("1-y", sym, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("x", sym, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
	PowerSeries_t* ps_array[2];
	ps_array[0] = ps_3;
	ps_array[1] = quo_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 13: \t\t\t\t PASSED \n");
}


/**
 *
 */
void test14(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};

	Poly_ptr  poly_2 = generate_altarr_var_defined("x", sym, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("1-y-x^2", sym, 3);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* ps_array[2];
	ps_array[0] =  ps_2;
	ps_array[1] = ps_3;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);

    fprintf(stderr, "Weierstrass Test 14: \t\t\t\t PASSED \n");
}



/**
 * Memory consumption
 */
void test15(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
	Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("x*y+x^3", sym, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("1-y^2+x^2", sym, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_3);
	PowerSeries_t* ps_array[2];
	ps_array[0] = ps_2;
	ps_array[1] = quo_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 15: \t\t\t\t PASSED \n");
}

/**
 * Memory consumption
 */
void test16(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("x^2+y^2", sym, 3);
    Poly_ptr  poly_3 = generate_altarr_var_defined("1-y^2+2*x^3", sym, 3);
    PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
    PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
    PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
    PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_3);
    PowerSeries_t* ps_array[4];
    ps_array[0] = ps_2;
    ps_array[1] = zeroPowerSeries_PS();
    ps_array[2] = onePowerSeries_PS(3);
    ps_array[3] = quo_1;
    Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,4);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_array[1]);
    destroyPowerSeries_PS(ps_array[2]);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 16: \t\t\t\t PASSED \n");

}



 /**
 * It seems Maple is faster
 */
void test17(int testDeg)
{
	const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("1-x", sym, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("y-x^2", sym, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
	//PowerSeries_t* sub_1 = subPowerSeries_PS(quo_1,onePowerSeries_PS());
	homogPart_PS(2, quo_1);
	PowerSeries_t* ps_array[2];
	ps_array[0] = ps_3;
	ps_array[1] = quo_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);

	// print_UPOPS(stderr, upop , sym);
	TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 17: \t\t\t\t PASSED \n");


}


/**
 * Maple example
 */
void test18(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("1-x", sym, 3);
    Poly_ptr  poly_3 = generate_altarr_var_defined("1-y", sym, 3);
    PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
    PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
    PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
    PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
    homogPart_PS(6, quo_1);
    PowerSeries_t* quo_2 = dividePowerSeries_PS(ps_1, ps_3);
    homogPart_PS(7, quo_2);
    PowerSeries_t* ps_array[2];
    ps_array[0] = quo_1;
    ps_array[1] = quo_2;
    Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);

    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);
    destroyPowerSeries_PS(quo_2);

    fprintf(stderr, "Weierstrass Test 18: \t\t\t\t PASSED \n");
}


/**
 * Paper example1
 */
void test19(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("1+x+y", sym, 3);
    Poly_ptr  poly_3 = generate_altarr_var_defined("y", sym, 3);
    Poly_ptr  poly_4 = generate_altarr_var_defined("x", sym, 3);

    PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
    PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
    PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
    PowerSeries_t* ps_4 = convertPolyToPowerSeries_PS(poly_4);
    PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);

    PowerSeries_t* ps_array[4];
    ps_array[0] = ps_4;
    ps_array[1] = ps_3;
    ps_array[2] = ps_1;
    ps_array[3] = quo_1;
    Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,4);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);
    freePolynomial_AA(poly_4);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(ps_4);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 19: \t\t\t\t PASSED \n");
}

/**
 * Paper example2
 */
void test20(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("1+x+y", sym, 3);
    Poly_ptr  poly_3 = generate_altarr_var_defined("y", sym, 3);
    Poly_ptr  poly_4 = generate_altarr_var_defined("x", sym, 3);

    PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
    PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
    PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
    PowerSeries_t* ps_4 = convertPolyToPowerSeries_PS(poly_4);
    PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);

    PowerSeries_t* ps_array[4];
    ps_array[0] = ps_4;
    ps_array[1] = ps_1;
    ps_array[2] = ps_3;
    ps_array[3] = quo_1;
    Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,4);
    TestingWeierstrassPreparation(upop, testDeg, sym);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);
    freePolynomial_AA(poly_4);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(ps_4);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 20: \t\t\t\t PASSED \n");
}


/**
 * Paper example2
 */
void test21(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("1+x+y", sym, 3);
    Poly_ptr  poly_3 = generate_altarr_var_defined("y", sym, 3);
    Poly_ptr  poly_4 = generate_altarr_var_defined("y+1", sym, 3);
    Poly_ptr  poly_5 = generate_altarr_var_defined("x", sym, 3);

    PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
    PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
    PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
    PowerSeries_t* ps_4 = convertPolyToPowerSeries_PS(poly_4);
    PowerSeries_t* ps_5 = convertPolyToPowerSeries_PS(poly_5);
    PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);

    PowerSeries_t* ps_array[5];
    ps_array[0] = ps_5;
    ps_array[1] = ps_4;
    ps_array[2] = zeroPowerSeries_PS();
    ps_array[3] = ps_3;
    ps_array[4] = quo_1;
    Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,5);
    TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);
    freePolynomial_AA(poly_4);
    freePolynomial_AA(poly_5);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_array[2]);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(ps_4);
    destroyPowerSeries_PS(ps_5);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 21: \t\t\t\t PASSED \n");
}


void test22(int testDeg)
{
    const char* sym[] = {"z", "y", "x"};
    Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("1+x+y", sym, 3);
    Poly_ptr  poly_3 = generate_altarr_var_defined("y", sym, 3);
    // Poly_ptr  poly_4 = generate_altarr_var_defined("y+1", sym, 3);
    Poly_ptr  poly_5 = generate_altarr_var_defined("x", sym, 3);
    PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
    PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
    PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
    // PowerSeries_t* ps_4 = convertPolyToPowerSeries_PS(poly_4);
    PowerSeries_t* ps_5 = convertPolyToPowerSeries_PS(poly_5);
    PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
    PowerSeries_t* ps_array[5];
    ps_array[0] = ps_5;
    ps_array[1] = ps_3;
    ps_array[2] = onePowerSeries_PS(3);
    ps_array[3] = onePowerSeries_PS(3);
    ps_array[4] = quo_1;
    Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,5);
    TestingWeierstrassPreparation(upop, testDeg, sym);


    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);
    freePolynomial_AA(poly_5);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_array[2]);
    destroyPowerSeries_PS(ps_array[3]);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(ps_5);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Weierstrass Test 22: \t\t\t\t PASSED \n");

}
