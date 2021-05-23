#include "../../include/PowerSeries/PowerSeries.h"
#include "../../include/PowerSeries/UnivariatePolynomialOverPowerSeries.h"
#include "../../include/Parser/bpas_parser.h"

const char* sym[] = {"x", "y"};
const char* syms[] = {"x", "y", "z", "t"};
const char* vars[] = {"x", "y", "z"};
const char* Vars[] = {"z", "x", "y"};
 /**
* Upops example 1
*/
int test1()
{


	Poly_ptr a = generate_altarr_var_defined("x^4*t+z^5-2*y^5+x*t^2*z+1", syms, 4);
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("t", syms, 4);
	Poly_ptr mainVarPoly_1 = upop->data[4]->polys[1];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("t^2*z", syms, 4);
	Poly_ptr mainVarPoly_2 = upop->data[1]->polys[3];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("z^5-2*y^5", syms, 4);
	Poly_ptr mainVarPoly_3 = upop->data[0]->polys[5];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_4 = upop->data[0]->polys[0];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);

	if (!isZero_AA(subpoly_1) || !isZero_AA(subpoly_2) || !isZero_AA(subpoly_3) || !isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 1: \t\t\t\t\t FAILED \n");
		exit(1);
	} else {
	    fprintf(stderr, "UPOPS Test 1: \t\t\t\t\t PASSED \n");
	}

	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);

	return 0;
}


/**
 * Upops example 2
 */
int test2()
{

	Poly_ptr a = generate_altarr_var_defined("2", vars, 3);
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Poly_ptr testpoly = generate_altarr_var_defined("2", vars, 3);
	Poly_ptr mainVarPoly = upop->data[0]->polys[0];
	testpoly = subPolynomials_AA_inp(testpoly,mainVarPoly,4);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "UPOPS Test 2: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 2: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	freePolynomial_AA(a);
	freePolynomial_AA(testpoly);
    return 0;
}

/**
 * Upops example 3
 */
int test3()
{

	Poly_ptr a = generate_altarr_var_defined("2+x", vars, 3);
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("1", vars, 3);
	Poly_ptr mainVarPoly_1  =  upop->data[1]->polys[0];
	testpoly_1 = subPolynomials_AA_inp(testpoly_1,mainVarPoly_1,3);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("2", vars, 3);
	Poly_ptr mainVarPoly_2  =  upop->data[0]->polys[0];
	testpoly_2 = subPolynomials_AA_inp(testpoly_2,mainVarPoly_2,3);

	if (isZero_AA(testpoly_2) && isZero_AA(testpoly_2)) {
	    fprintf(stderr, "UPOPS Test 3: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 3: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    return 0;
}

/**
 * Upops example 4
 */
int test4()
{

	Poly_ptr a = generate_altarr_var_defined("4+3*z*x^2+y*t*x^4+y", syms, 4);
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("y*t", syms, 4);
	Poly_ptr mainVarPoly_1 = upop->data[4]->polys[2];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("3*z", syms, 4);
	Poly_ptr mainVarPoly_2 = upop->data[2]->polys[1];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("4", syms, 4);
	Poly_ptr mainVarPoly_3 = upop->data[0]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("y", syms, 4);
	Poly_ptr mainVarPoly_4 = upop->data[0]->polys[1];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 4: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 4: \t\t\t\t\t FAILED \n");
		exit(1);
	}


	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);

    return 0;
}



/**
 * Upops example 5
 */
int test5()
{

	Poly_ptr a = generate_altarr_var_defined("0", syms, 4);
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	Poly_ptr mainVarPoly = upop->data[0]->polys[0];
	if (isZero_AA( mainVarPoly))  {
	    fprintf(stderr, "UPOPS Test 5: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 5: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    return 0;
}


/**
 * zero Upops
 */
int test6()
{

	Upops_t* zero = zero_UPOPS();
	Poly_ptr mainVarPoly = zero->data[0]->polys[0];
	if (isZero_AA( mainVarPoly)) {
	    fprintf(stderr, "UPOPS Test 6: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 6: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
    return 0;
}


/**
 * one Upops
 */
int test7()
{

	Upops_t* one = one_UPOPS(4);
	Poly_ptr testpoly = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly = one->data[0]->polys[0];
	Poly_ptr subpoly = subPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "UPOPS Test 7: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 7: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(testpoly);
	freePolynomial_AA(subpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
    return 0;
}


/**
 * Upops addition example 1 : addition of two zero Upops
 */
int test8()
{

	Upops_t* zero = zero_UPOPS();
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (zero, zero);
	Poly_ptr testpoly = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly = sum->data[0]->polys[0];
	Poly_ptr subpoly = subPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "UPOPS Test 8: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 8: \t\t\t\t\t FAILED \n");
		exit(1);
    }

    freePolynomial_AA(subpoly);
    freePolynomial_AA(testpoly);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);
	return 0;
}


/**
 * Upops addition example 2 : addition of two one Upops
 */
int test9()
{

	Upops_t* one = one_UPOPS(4);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (one, one);
	Poly_ptr testpoly = generate_altarr_var_defined("2", syms, 4);
	Poly_ptr mainVarPoly =  sum->data[0]->polys[0];
	Poly_ptr subpoly = subPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "UPOPS Test 9: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 9: \t\t\t\t\t FAILED \n");
		exit(1);
    }

    freePolynomial_AA(subpoly);
    freePolynomial_AA(testpoly);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
    destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);
	return 0;
}


/**
 * Upops addition example 3 : addition of zero Upops and one Upops
 */
int test10()
{
	Upops_t* zero = zero_UPOPS();
	Upops_t* one = one_UPOPS(4);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (zero, one);
	Poly_ptr testpoly = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly = sum->data[0]->polys[0];
	Poly_ptr subpoly = subPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "UPOPS Test 10: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 10: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly);
	freePolynomial_AA(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);
    return 0;
}

/**
 * Upops addition example 4 : addition of zero Upops and non-constant Upops
 */
int test11()
{
	Upops_t*  zero = zero_UPOPS();
	Poly_ptr poly = generate_altarr_var_defined("4+3*z*x^2+y*t*x^5+y", syms, 4);
	Upops_t*  upops = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (zero, upops);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("3*z", syms, 4);
	Poly_ptr mainVarPoly_1 = sum->data[2]->polys[1];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("y*t", syms, 4);
	Poly_ptr mainVarPoly_2 = sum->data[5]->polys[2];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("4", syms, 4);
	Poly_ptr mainVarPoly_3 = sum->data[0]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("y", syms, 4);
	Poly_ptr mainVarPoly_4 = sum->data[0]->polys[1];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 11: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 11: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(poly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);


    return 0;
}



/**
 * Upops addition example 5 : addition of zero Upops and non-constant Upops
 */
int test12()
{
	Poly_ptr poly1 = generate_altarr_var_defined("z", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+y", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("y*t", syms, 4);
	updateToDeg_UPOPS(2, sum);
	Poly_ptr mainVarPoly_1 = homogPart_PS(2,  sum->data[6]);

	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_2 = sum->data[6]->polys[0];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("3*z", syms, 4);
	Poly_ptr mainVarPoly_3 = sum->data[2]->polys[1];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("y+z", syms, 4);
	Poly_ptr mainVarPoly_4 = sum->data[0]->polys[1];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 12: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 12: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);
    return 0;
}

/**
 * Upops addition example 6 : addition of two Upops
 */
int test13()
{
	Poly_ptr poly1 = generate_altarr_var_defined("3*z*x^2+(1+y*t^6)*x^6+y+x", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("z+2*y*x^2", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(8, sum);

	Poly_ptr testpoly_1 = generate_altarr_var_defined("y+z", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(1,  sum->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_2 = sum->data[1]->polys[0];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 =  NULL;
	Poly_ptr mainVarPoly_3 = sum->data[3]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("2*y+3*z", syms, 4);
	Poly_ptr mainVarPoly_4 = sum->data[2]->polys[1];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);
	Poly_ptr testpoly_5 = generate_altarr_var_defined("y*t^6", syms, 4);
	Poly_ptr mainVarPoly_5 = sum->data[6]->polys[7];
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,4);
	Poly_ptr testpoly_6 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_6 = sum->data[6]->polys[0];
	Poly_ptr subpoly_6 = subPolynomials_AA(mainVarPoly_6,testpoly_6,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5) && isZero_AA(subpoly_6)) {
	    fprintf(stderr, "UPOPS Test 13: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 13: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(subpoly_6);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(testpoly_6);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);

    return 0;
}


/**
 * Upops addition example 7 : addition of two Upops
 */
int test14()
{
	Poly_ptr poly1 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+y", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("z^6", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(2, sum);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(5,  sum->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("y", syms, 4);
	Poly_ptr mainVarPoly_2 = sum->data[0]->polys[1];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("3*z", syms, 4);
	Poly_ptr mainVarPoly_3 = sum->data[2]->polys[1];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("y*t", syms, 4);
	Poly_ptr mainVarPoly_4 = sum->data[6]->polys[2];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);
	Poly_ptr testpoly_5 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_5 = sum->data[6]->polys[0];
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5)) {
	    fprintf(stderr, "UPOPS Test 14: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 14: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);

    return 0;
}

/**
 * Upops subtraction example 1
 */
int test15()
{

	Upops_t* zero = zero_UPOPS();
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (zero, zero);
	Poly_ptr testpoly = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly = sub->data[0]->polys[0];
	Poly_ptr subpoly = subPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "UPOPS Test 15: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 15: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly);
	freePolynomial_AA(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);

    return 0;
}




/**
 * Upops subtraction example 2
 */
int test16()
{

	Upops_t* one = one_UPOPS(4);
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (one, one);
	updateToDeg_UPOPS(3, sub);
	Poly_ptr mainVarPoly =  sub->data[0]->polys[0];
	if (isZero_AA(mainVarPoly)) {
	    fprintf(stderr, "UPOPS Test 16: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 16: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
    return 0;
}


/**
 * Upops subtraction example 3
 */
int test17()
{
	Upops_t* zero = zero_UPOPS();
	Upops_t* one = one_UPOPS(4);
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (zero, one);
	updateToDeg_UPOPS(4, sub);
	Poly_ptr mainVarPoly = sub->data[0]->polys[0];
	if (isNegativeOne_AA(mainVarPoly)) {
	    fprintf(stderr, "UPOPS Test 17: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 17: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
    return 0;
}

/**
 * Upops subtraction example 4
 */
int test18()
{
	Upops_t* zero = zero_UPOPS();
	Poly_ptr poly = generate_altarr_var_defined("4+3*z*x^2+y*t*x^4+y", syms, 4);
	Upops_t*  upops = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly);
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (zero, upops);
	updateToDeg_UPOPS(2, sub);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("-3*z", syms, 4);
	Poly_ptr mainVarPoly_1 = sub->data[2]->polys[1];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("-4", syms, 4);
	Poly_ptr mainVarPoly_2 = sub->data[0]->polys[0];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("-y", syms, 4);
	Poly_ptr mainVarPoly_3 = sub->data[0]->polys[1];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("-y*t", syms, 4);
	Poly_ptr mainVarPoly_4 = sub->data[4]->polys[2];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 18: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 18: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(poly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
    return 0;
}



/**
 * Upops subtraction example 5: issue with printPolyClean
 */
int test19()
{
	Poly_ptr poly1 = generate_altarr_var_defined("z", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("3*z*x^2+(1+y^3*t)*x^6+y", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(2, sub);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("-y^3*t", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(4,  sub->data[6]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("z-y", syms, 4);
	Poly_ptr mainVarPoly_2 = sub->data[0]->polys[1];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("-3*z", syms, 4);
	Poly_ptr mainVarPoly_3 = sub->data[2]->polys[1];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("-1", syms, 4);
	Poly_ptr mainVarPoly_4 = sub->data[6]->polys[0];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);


	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 19: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 19: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
    return 0;
}

/**
 * Upops subtraction example 6
 */
int test20()
{
	Poly_ptr poly1 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+y", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("z", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(2, sub);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("y-z", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(1,  sub->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("y*t", syms, 4);
	Poly_ptr mainVarPoly_2 = sub->data[6]->polys[2];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_3 = sub->data[6]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("3*z", syms, 4);
	Poly_ptr mainVarPoly_4 = sub->data[2]->polys[1];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 20: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 20: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);

    return 0;
}


/**
 * Upops subtraction example 7
 */
int test21()
{
	Poly_ptr poly1 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+t^8", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("z^6", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(6, sub);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(5,  sub->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("-z^6", syms, 4);
	Poly_ptr mainVarPoly_2 = sub->data[0]->polys[6];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("y*t", syms, 4);
	Poly_ptr mainVarPoly_3 = sub->data[6]->polys[2];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_4 = sub->data[6]->polys[0];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);
	Poly_ptr testpoly_5 = generate_altarr_var_defined("3*z", syms, 4);
	Poly_ptr mainVarPoly_5 = sub->data[2]->polys[1];
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5)) {
	    fprintf(stderr, "UPOPS Test 21: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 21: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);

    return 0;
}


/**
 * Upops multiplication example 1
 */
int test22()
{

	Upops_t* zero = zero_UPOPS();
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (zero, zero);
	Poly_ptr testpoly = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly = prod->data[0]->polys[0];
	Poly_ptr addpoly = addPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(addpoly)) {
	    fprintf(stderr, "UPOPS Test 22: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 22: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(testpoly);
	freePolynomial_AA(addpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	return 0;

}

/**
 * Upops multiplication example 2
 */
int test23()
{

	Upops_t* one = one_UPOPS(4);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (one, one);
	Poly_ptr testpoly = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly =  prod->data[0]->polys[0];
	Poly_ptr subpoly = subPolynomials_AA(mainVarPoly,testpoly,4);
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "UPOPS Test 23: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 23: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly);
	freePolynomial_AA(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
    return 0;
}


/**
 * Upops multiplication  example 3
 */
int test24()
{
	Upops_t* zero = zero_UPOPS();
	Upops_t* one = one_UPOPS(4);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (zero, one);
	Poly_ptr testpoly = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly = prod->data[0]->polys[0];
	Poly_ptr prodpoly = subPolynomials_AA(mainVarPoly,testpoly,4);

	if (isZero_AA(prodpoly)) {
	    fprintf(stderr, "UPOPS Test 24: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 24: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(prodpoly);
	freePolynomial_AA(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);
    return 0;
}


/**
 * Upops multiplication  example 4
 */
int test25()
{
	Upops_t* zero = zero_UPOPS();
	Poly_ptr poly = generate_altarr_var_defined("4+3*z*x^2+y*t*x^4+y", syms, 4);
	Upops_t*  upops = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (zero, upops);
	updateToDeg_UPOPS(1, prod);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("0", syms, 4);
	Poly_ptr mainVarPoly_1 = prod->data[2]->polys[0];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr mainVarPoly_2 = prod->data[0]->polys[0];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_1,4);
	Poly_ptr mainVarPoly_3 = prod->data[1]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_1,4);
	Poly_ptr mainVarPoly_4 = prod->data[3]->polys[0];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_1,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 25: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 25: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(poly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(zero);

    return 0;
}


/**
 * Upops multiplication example 5
 */
int test26()
{
	Poly_ptr poly1 = generate_altarr_var_defined("z^2+x^2", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+y^3+1", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(2, prod);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("y^3+3*z^3", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(3,  prod->data[2]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("z^2", syms, 4);
	Poly_ptr mainVarPoly_2 = prod->data[0]->polys[2];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = NULL;
	Poly_ptr mainVarPoly_3 = prod->data[1]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr mainVarPoly_4 = prod->data[2]->polys[0];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,4);
	Poly_ptr mainVarPoly_5 = prod->data[3]->polys[0];
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_3,4);
	Poly_ptr mainVarPoly_6 = prod->data[4]->polys[0];
	Poly_ptr subpoly_6 = subPolynomials_AA(mainVarPoly_6,testpoly_3,4);
	Poly_ptr mainVarPoly_7 = prod->data[5]->polys[0];
	Poly_ptr subpoly_7 = subPolynomials_AA(mainVarPoly_7,testpoly_3,4);
	Poly_ptr mainVarPoly_8 = prod->data[8]->polys[0];
	Poly_ptr subpoly_8 = subPolynomials_AA
	(mainVarPoly_8,testpoly_4,4);
	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5) && isZero_AA(subpoly_6) && isZero_AA(subpoly_7) && isZero_AA(subpoly_8)) {
	    fprintf(stderr, "UPOPS Test 26: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 26: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(subpoly_6);
	freePolynomial_AA(subpoly_7);
	freePolynomial_AA(subpoly_8);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);

	return 0;

}


/**
 * Upops multiplication example 6
 */
int test27()
{
	Poly_ptr poly1 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+y+1", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("t", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(2, prod);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("t*y", syms, 4);
	Poly_ptr mainVarPoly_1 = homogPart_PS(2,  prod->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("t", syms, 4);
	Poly_ptr mainVarPoly_2 = prod->data[0]->polys[1];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("3*t*z", syms, 4);
	Poly_ptr mainVarPoly_3 = homogPart_PS(2,  prod->data[2]);
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);
	Poly_ptr mainVarPoly_4 = prod->data[6]->polys[1];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_2,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 27: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 27: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
    return 0;
}


/**
 * Upops multiplication example 7
 */
int test28()
{
	Poly_ptr poly1 = generate_altarr_var_defined("3*z*x^2+(1+y*t)*x^6+t+z^3", syms, 4);
	Poly_ptr poly2 = generate_altarr_var_defined("z^6+t^2", syms, 4);
	Upops_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly1);
	Upops_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly2);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (upops1, upops2);
	updateToDeg_UPOPS(3, prod);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("t^3", syms, 4);
	Poly_ptr mainVarPoly_1 = prod->data[0]->polys[3];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,4);
	Poly_ptr testpoly_2 = NULL;
	Poly_ptr mainVarPoly_2 = prod->data[2]->polys[0];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,4);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("t^2", syms, 4);
	Poly_ptr mainVarPoly_3 = prod->data[6]->polys[2];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,4);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3)) {
	    fprintf(stderr, "UPOPS Test 28: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 28: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(poly1);
	freePolynomial_AA(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
    return 0;
}

/**
 *  Upops addition example 8: using power series as coefficients of the upops
 */
int test29()
{
	Poly_ptr  poly_1 = generate_altarr_var_defined("1", Vars, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("1-x", Vars, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("1-y", Vars, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
	homogPart_PS(2, quo_1);
	PowerSeries_t* quo_2 = dividePowerSeries_PS(ps_1, ps_3);
	homogPart_PS(2, quo_2);

	PowerSeries_t* ps_array[2];
	ps_array[0] = quo_1;
	ps_array[1] = quo_2;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);
	Upops_t* sum = addUnivariatePolyOverPowerSeries_UPOPS (upop, upop);
	updateToDeg_UPOPS(2, sum);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("2", Vars, 3);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("2*x", Vars, 3);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("2*x^2", Vars, 3);
	Poly_ptr mainVarPoly_1 = sum->data[0]->polys[0];
	Poly_ptr mainVarPoly_2 = sum->data[0]->polys[1];
	Poly_ptr mainVarPoly_3 = sum->data[0]->polys[2];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,3);
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,3);
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,3);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("2", Vars, 3);
	Poly_ptr testpoly_5 = generate_altarr_var_defined("2*y", Vars, 3);
	Poly_ptr testpoly_6 = generate_altarr_var_defined("2*y^2", Vars, 3);
	Poly_ptr mainVarPoly_4 = sum->data[1]->polys[0];
	Poly_ptr mainVarPoly_5 = sum->data[1]->polys[1];
	Poly_ptr mainVarPoly_6 = sum->data[1]->polys[2];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,3);
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,3);
	Poly_ptr subpoly_6 = subPolynomials_AA(mainVarPoly_6,testpoly_6,3);



	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5) && isZero_AA(subpoly_6)) {
	    fprintf(stderr, "UPOPS Test 29: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 29: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(subpoly_6);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(testpoly_6);
	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(quo_1);
	destroyPowerSeries_PS(quo_2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sum);
    return 0;

}



/**
 *  Upops multiplication example 8: using power series as coefficients of the upops
 */
int test30()
{
	Poly_ptr  poly_1 = generate_altarr_var_defined("1", Vars, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("1-x", Vars, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("1-y", Vars, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
	homogPart_PS(2, quo_1);
	PowerSeries_t* quo_2 = dividePowerSeries_PS(ps_1, ps_3);
	homogPart_PS(2, quo_2);
	PowerSeries_t* ps_array[2];
	ps_array[0] = quo_1;
	ps_array[1] = quo_2;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (upop, upop);
	updateToDeg_UPOPS(2, prod);
	Poly_ptr testpoly_1 = generate_altarr_var_defined("1", Vars, 3);
	Poly_ptr testpoly_2 = generate_altarr_var_defined("2*x", Vars, 3);
	Poly_ptr testpoly_3 = generate_altarr_var_defined("3*x^2", Vars, 3);
	Poly_ptr testpoly_4 = generate_altarr_var_defined("4*x^3", Vars, 3);
	Poly_ptr testpoly_5 = generate_altarr_var_defined("5*x^4", Vars, 3);
	Poly_ptr mainVarPoly_1 = prod->data[0]->polys[0];
	Poly_ptr mainVarPoly_2 = prod->data[0]->polys[1];
	Poly_ptr mainVarPoly_3 = prod->data[0]->polys[2];
	Poly_ptr mainVarPoly_4 = homogPart_PS(3,  prod->data[0]);
	Poly_ptr mainVarPoly_5 = homogPart_PS(4,  prod->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,3);
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,3);
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,3);
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,3);
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,3);

	Poly_ptr testpoly_6 = generate_altarr_var_defined("2", Vars, 3);
	Poly_ptr testpoly_7 = generate_altarr_var_defined("2*x+2*y", Vars, 3);
	Poly_ptr testpoly_8 = generate_altarr_var_defined("2*x^2+2*x*y+2*y^2", Vars, 3);
	Poly_ptr testpoly_9 = generate_altarr_var_defined("2*x^3+2*y^3+2*x^2*y+2*x*y^2", Vars, 3);
	Poly_ptr testpoly_10 = generate_altarr_var_defined("2*x^4+2*y^4+2*x*y^3+2*x^3*y+2*x^2*y^2", Vars, 3);
	Poly_ptr mainVarPoly_6 = prod->data[1]->polys[0];
	Poly_ptr mainVarPoly_7 = prod->data[1]->polys[1];
	Poly_ptr mainVarPoly_8 = prod->data[1]->polys[2];
	Poly_ptr mainVarPoly_9 = homogPart_PS(3,  prod->data[1]);
	Poly_ptr mainVarPoly_10 = homogPart_PS(4,  prod->data[1]);
	Poly_ptr subpoly_6 = subPolynomials_AA(mainVarPoly_6,testpoly_6,3);
	Poly_ptr subpoly_7 = subPolynomials_AA(mainVarPoly_7,testpoly_7,3);
	Poly_ptr subpoly_8 = subPolynomials_AA(mainVarPoly_8,testpoly_8,3);
	Poly_ptr subpoly_9 = subPolynomials_AA(mainVarPoly_9,testpoly_9,3);
	Poly_ptr subpoly_10 = subPolynomials_AA(mainVarPoly_10,testpoly_10,3);


	Poly_ptr testpoly_11 = generate_altarr_var_defined("1", Vars, 3);
	Poly_ptr testpoly_12 = generate_altarr_var_defined("2*y", Vars, 3);
	Poly_ptr testpoly_13 = generate_altarr_var_defined("3*y^2", Vars, 3);
	Poly_ptr testpoly_14 = generate_altarr_var_defined("4*y^3", Vars, 3);
	Poly_ptr testpoly_15 = generate_altarr_var_defined("5*y^4", Vars, 3);
	Poly_ptr mainVarPoly_11 = prod->data[2]->polys[0];
	Poly_ptr mainVarPoly_12 = prod->data[2]->polys[1];
	Poly_ptr mainVarPoly_13 = prod->data[2]->polys[2];
	Poly_ptr mainVarPoly_14 = homogPart_PS(3,  prod->data[2]);
	Poly_ptr mainVarPoly_15 = homogPart_PS(4,  prod->data[2]);
	Poly_ptr subpoly_11 = subPolynomials_AA(mainVarPoly_11,testpoly_11,3);
	Poly_ptr subpoly_12 = subPolynomials_AA(mainVarPoly_12,testpoly_12,3);
	Poly_ptr subpoly_13 = subPolynomials_AA(mainVarPoly_13,testpoly_13,3);
	Poly_ptr subpoly_14 = subPolynomials_AA(mainVarPoly_14,testpoly_14,3);
	Poly_ptr subpoly_15 = subPolynomials_AA(mainVarPoly_15,testpoly_15,3);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5) && isZero_AA(subpoly_6) && isZero_AA(subpoly_7) && isZero_AA(subpoly_8)
		 && isZero_AA(subpoly_9) && isZero_AA(subpoly_10) && isZero_AA(subpoly_11) && isZero_AA(subpoly_12)
		  && isZero_AA(subpoly_13) && isZero_AA(subpoly_14) && isZero_AA(subpoly_15)) {
	    fprintf(stderr, "UPOPS Test 30: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 30: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(subpoly_6);
	freePolynomial_AA(subpoly_7);
	freePolynomial_AA(subpoly_8);
	freePolynomial_AA(subpoly_9);
	freePolynomial_AA(subpoly_10);
	freePolynomial_AA(subpoly_11);
	freePolynomial_AA(subpoly_12);
	freePolynomial_AA(subpoly_13);
	freePolynomial_AA(subpoly_14);
	freePolynomial_AA(subpoly_15);

	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(testpoly_6);
	freePolynomial_AA(testpoly_7);
	freePolynomial_AA(testpoly_8);
	freePolynomial_AA(testpoly_9);
	freePolynomial_AA(testpoly_10);
	freePolynomial_AA(testpoly_11);
	freePolynomial_AA(testpoly_12);
	freePolynomial_AA(testpoly_13);
	freePolynomial_AA(testpoly_14);
	freePolynomial_AA(testpoly_15);

	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(quo_1);
	destroyPowerSeries_PS(quo_2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
    return 0;

}


/**
 *  Upops multiplication example 9: using power series as coefficients of the upops
 */
int test31()

{

	Poly_ptr  poly_1 = generate_altarr_var_defined("1", Vars, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("x*y+x^4", Vars, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("1-y^2+x^2", Vars, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_3);
	homogPart_PS(6, quo_1);
	PowerSeries_t* ps_array[2];
	ps_array[0] = ps_2;
	ps_array[1] = quo_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,2);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (upop, upop);
	updateToDeg_UPOPS(6, prod);

	Poly_ptr testpoly_1 = generate_altarr_var_defined("x^2*y^2", Vars, 3);
	Poly_ptr mainVarPoly_1 = prod->data[0]->polys[4];
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,3);

	Poly_ptr testpoly_2 = generate_altarr_var_defined("2*x^4+2*x*y^3-2*x^3*y", Vars, 3);
	Poly_ptr mainVarPoly_2 = prod->data[1]->polys[4];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,3);

	Poly_ptr testpoly_3 = generate_altarr_var_defined("2*x*y", Vars, 3);
	Poly_ptr mainVarPoly_3 = prod->data[1]->polys[2];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,3);

	Poly_ptr testpoly_4 = generate_altarr_var_defined("1", Vars, 3);
	Poly_ptr mainVarPoly_4 = prod->data[2]->polys[0];
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,3);

	Poly_ptr testpoly_5 = generate_altarr_var_defined("2*y^2-2*x^2", Vars, 3);
	Poly_ptr mainVarPoly_5 = prod->data[2]->polys[2];
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,3);

	Poly_ptr testpoly_6 = generate_altarr_var_defined("3*y^4-6*x^2*y^2+3*x^4", Vars, 3);
	Poly_ptr mainVarPoly_6 = prod->data[2]->polys[4];
	Poly_ptr subpoly_6 = subPolynomials_AA(mainVarPoly_6,testpoly_6,3);

	Poly_ptr testpoly_7 = generate_altarr_var_defined("-4*x^6+12*x^4*y^2-12*x^2*y^4+4*y^6", Vars, 3);
	Poly_ptr mainVarPoly_7 = prod->data[2]->polys[6];
	Poly_ptr subpoly_7 = subPolynomials_AA(mainVarPoly_7,testpoly_7,3);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5) && isZero_AA(subpoly_6) && isZero_AA(subpoly_7) ) {
	    fprintf(stderr, "UPOPS Test 31: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 31: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(subpoly_6);
	freePolynomial_AA(subpoly_7);

	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(testpoly_6);
	freePolynomial_AA(testpoly_7);

	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(quo_1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
    return 0;

}


/**
 *  Upops multiplication example 10: using power series as coefficients of the upops
 */
int test32()

{
	Poly_ptr  poly_1 = generate_altarr_var_defined("x+1", Vars, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("x^2+y^2", Vars, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("1-y^2+2*x^3", Vars, 3);
	Poly_ptr  poly_4 = generate_altarr_var_defined("x", Vars, 3);
	Poly_ptr  poly_5 = generate_altarr_var_defined("x+y", Vars, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_3);
	homogPart_PS(4, quo_1);
	PowerSeries_t* ps_4 = convertPolyToPowerSeries_PS(poly_4);
	PowerSeries_t* ps_5 = convertPolyToPowerSeries_PS(poly_5);
	PowerSeries_t* ps_array[4];
	ps_array[0] = ps_2;
	ps_array[1] =  ps_4;
	ps_array[2] = ps_5;
	ps_array[3] = quo_1;

	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,4);
	Upops_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPS (upop, upop);
	updateToDeg_UPOPS(4, prod);

	Poly_ptr testpoly_1 = generate_altarr_var_defined("2*x^2*y^2+x^4+y^4", Vars, 3);
	Poly_ptr mainVarPoly_1 = homogPart_PS(4, prod->data[0]);
	Poly_ptr subpoly_1 = subPolynomials_AA(mainVarPoly_1,testpoly_1,3);

	Poly_ptr testpoly_2 = NULL;
	Poly_ptr mainVarPoly_2 = prod->data[1]->polys[0];
	Poly_ptr subpoly_2 = subPolynomials_AA(mainVarPoly_2,testpoly_2,3);


	Poly_ptr testpoly_3 =  NULL;
	Poly_ptr mainVarPoly_3 = prod->data[2]->polys[0];
	Poly_ptr subpoly_3 = subPolynomials_AA(mainVarPoly_3,testpoly_3,3);

	Poly_ptr testpoly_4 = generate_altarr_var_defined("4*x^2+2*y^2+2*x*y", Vars, 3);;
	Poly_ptr mainVarPoly_4 = homogPart_PS(2,prod->data[3]);
	Poly_ptr subpoly_4 = subPolynomials_AA(mainVarPoly_4,testpoly_4,3);

	Poly_ptr testpoly_5 = generate_altarr_var_defined("2*x", Vars, 3);;
	Poly_ptr mainVarPoly_5 = prod->data[4]->polys[1];
	Poly_ptr subpoly_5 = subPolynomials_AA(mainVarPoly_5,testpoly_5,3);


	Poly_ptr testpoly_6 = generate_altarr_var_defined("2*x+2*y", Vars, 3);;
	Poly_ptr mainVarPoly_6 = prod->data[5]->polys[1];
	Poly_ptr subpoly_6 = subPolynomials_AA(mainVarPoly_6,testpoly_6,3);

	Poly_ptr testpoly_7 = generate_altarr_var_defined("1", Vars, 3);;
	Poly_ptr mainVarPoly_7 = prod->data[6]->polys[0];
	Poly_ptr subpoly_7 = subPolynomials_AA(mainVarPoly_7,testpoly_7,3);

	Poly_ptr testpoly_8 = generate_altarr_var_defined("2*x", Vars, 3);;
	Poly_ptr mainVarPoly_8 = prod->data[6]->polys[1];
	Poly_ptr subpoly_8 = subPolynomials_AA(mainVarPoly_8,testpoly_8,3);

	Poly_ptr testpoly_9 = generate_altarr_var_defined("2*y^2+x^2", Vars, 3);;
	Poly_ptr mainVarPoly_9 = prod->data[6]->polys[2];
	Poly_ptr subpoly_9 = subPolynomials_AA(mainVarPoly_9,testpoly_9,3);

	Poly_ptr testpoly_10 = generate_altarr_var_defined("4*x*y^2-4*x^3", Vars, 3);;
	Poly_ptr mainVarPoly_10 = prod->data[6]->polys[3];
	Poly_ptr subpoly_10 = subPolynomials_AA(mainVarPoly_10,testpoly_10,3);

	Poly_ptr testpoly_11 = generate_altarr_var_defined("3*y^4+2*x^2*y^2-8*x^4", Vars, 3);;
	Poly_ptr mainVarPoly_11 = prod->data[6]->polys[4];
	Poly_ptr subpoly_11 = subPolynomials_AA(mainVarPoly_11,testpoly_11,3);

	if (isZero_AA(subpoly_1) && isZero_AA(subpoly_2) && isZero_AA(subpoly_3) && isZero_AA(subpoly_4)
		&& isZero_AA(subpoly_5) && isZero_AA(subpoly_6) && isZero_AA(subpoly_7) && isZero_AA(subpoly_8)
		&& isZero_AA(subpoly_9) && isZero_AA(subpoly_10) && isZero_AA(subpoly_11)) {
		fprintf(stderr, "UPOPS Test 32: \t\t\t\t\t PASSED \n");
	} else {
		fprintf(stderr, "UPOPS Test 32: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(subpoly_1);
	freePolynomial_AA(subpoly_2);
	freePolynomial_AA(subpoly_3);
	freePolynomial_AA(subpoly_4);
	freePolynomial_AA(subpoly_5);
	freePolynomial_AA(subpoly_6);
	freePolynomial_AA(subpoly_7);
	freePolynomial_AA(subpoly_8);
	freePolynomial_AA(subpoly_9);
	freePolynomial_AA(subpoly_10);
	freePolynomial_AA(subpoly_11);

	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);
	freePolynomial_AA(testpoly_3);
	freePolynomial_AA(testpoly_4);
	freePolynomial_AA(testpoly_5);
	freePolynomial_AA(testpoly_6);
	freePolynomial_AA(testpoly_7);
	freePolynomial_AA(testpoly_8);
	freePolynomial_AA(testpoly_9);
	freePolynomial_AA(testpoly_10);
	freePolynomial_AA(testpoly_11);

	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);
	freePolynomial_AA(poly_4);
	freePolynomial_AA(poly_5);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(ps_4);
	destroyPowerSeries_PS(ps_5);
	destroyPowerSeries_PS(quo_1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	return 0;
}


/**
 *
 */
int test33()
{
	 Poly_ptr  poly_1 = generate_altarr_var_defined("x+1", Vars, 3);
	 Poly_ptr  poly_2 = generate_altarr_var_defined("x^2+y^2", Vars, 3);
	 Poly_ptr  poly_3 = generate_altarr_var_defined("1-y^2+2*x^3", Vars, 3);
	 Poly_ptr  poly_4 = generate_altarr_var_defined("x", Vars, 3);
	 Poly_ptr  poly_5 = generate_altarr_var_defined("x+y", Vars, 3);
	 PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	 PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	 PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	 PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_3);
	 homogPart_PS(4, quo_1);
	 PowerSeries_t* ps_4 = convertPolyToPowerSeries_PS(poly_4);
	 PowerSeries_t* ps_5 = convertPolyToPowerSeries_PS(poly_5);
	 PowerSeries_t* ps_array[4];
	 ps_array[0] = ps_2;
	 ps_array[1] =  ps_4;
	 ps_array[2] = ps_5;
	 ps_array[3] = quo_1;
	 Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,4);
	 // print_UPOPS(stderr,  upop, Vars);
	 // fprintf(stderr, "\n");
	 // fprintf(stderr, "Up to degree 7 \n");
	 updateToDeg_UPOPS(7, upop);
	 // print_UPOPS(stderr,  upop, Vars);
	 // fprintf(stderr, "\n");
	 // fprintf(stderr, "Up to degree 17 \n");
	 updateToDeg_UPOPS(17, upop);
	 // print_UPOPS(stderr,  upop, Vars);


	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);
	freePolynomial_AA(poly_4);
	freePolynomial_AA(poly_5);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(ps_4);
	destroyPowerSeries_PS(ps_5);
	destroyPowerSeries_PS(quo_1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);

	//Same as previous example, just with higher degree
    fprintf(stderr, "UPOPS Test 33: \t\t\t\t\t PASSED \n");

	return 0;
}


// int test34()
// {
// 	PowerSeries_t* a = geometricSeries_PS(3);
//     //Poly_ptr s = homogPart_PS(9,a);
// 	Poly_ptr  poly_1 = generate_altarr_var_defined("1", vars, 3);
// 	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
// 	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, a);
// 	Poly_ptr s = homogPart_PS(1,quo_1);
// 	//print_PS(stderr, a,sym);
// 	//printPoly_AA(stderr, s, vars,3);
// 	//fprintf(stderr, "\n");
//     return 0;

// }
