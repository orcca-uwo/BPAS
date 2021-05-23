#include "../../include/PowerSeries/PowerSeriesZ.h"
#include "../../include/PowerSeries/UPOPS_Z.h"
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


	PolyZ_ptr a = generate_altarrZ_var_defined("x^4*t+z^5-2*y^5+x*t^2*z+1", syms, 4);
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("t", syms, 4);
	PolyZ_ptr mainVarPoly_1 = upop->data[4]->polys[1];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("t^2*z", syms, 4);
	PolyZ_ptr mainVarPoly_2 = upop->data[1]->polys[3];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("z^5-2*y^5", syms, 4);
	PolyZ_ptr mainVarPoly_3 = upop->data[0]->polys[5];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_4 = upop->data[0]->polys[0];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);

	if (!isZero_AAZ(subpoly_1) || !isZero_AAZ(subpoly_2) || !isZero_AAZ(subpoly_3) || !isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 1: \t\t\t\t\t FAILED \n");
		exit(1);
	} else {
	    fprintf(stderr, "UPOPS Test 1: \t\t\t\t\t PASSED \n");
	}

	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);

	return 0;
}


/**
 * Upops example 2
 */
int test2()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("2", vars, 3);
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("2", vars, 3);
	PolyZ_ptr mainVarPoly = upop->data[0]->polys[0];
	testpoly = subPolynomials_AAZ_inp(testpoly,mainVarPoly,4);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "UPOPS Test 2: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 2: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(testpoly);
    return 0;
}

/**
 * Upops example 3
 */
int test3()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("2+x", vars, 3);
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("1", vars, 3);
	PolyZ_ptr mainVarPoly_1  =  upop->data[1]->polys[0];
	testpoly_1 = subPolynomials_AAZ_inp(testpoly_1,mainVarPoly_1,3);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("2", vars, 3);
	PolyZ_ptr mainVarPoly_2  =  upop->data[0]->polys[0];
	testpoly_2 = subPolynomials_AAZ_inp(testpoly_2,mainVarPoly_2,3);

	if (isZero_AAZ(testpoly_2) && isZero_AAZ(testpoly_2)) {
	    fprintf(stderr, "UPOPS Test 3: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 3: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    return 0;
}

/**
 * Upops example 4
 */
int test4()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("4+3*z*x^2+y*t*x^4+y", syms, 4);
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("y*t", syms, 4);
	PolyZ_ptr mainVarPoly_1 = upop->data[4]->polys[2];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("3*z", syms, 4);
	PolyZ_ptr mainVarPoly_2 = upop->data[2]->polys[1];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("4", syms, 4);
	PolyZ_ptr mainVarPoly_3 = upop->data[0]->polys[0];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("y", syms, 4);
	PolyZ_ptr mainVarPoly_4 = upop->data[0]->polys[1];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 4: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 4: \t\t\t\t\t FAILED \n");
		exit(1);
	}


	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);

    return 0;
}



/**
 * Upops example 5
 */
int test5()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("0", syms, 4);
	UpopsZ_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(a);
	PolyZ_ptr mainVarPoly = upop->data[0]->polys[0];
	if (isZero_AAZ( mainVarPoly))  {
	    fprintf(stderr, "UPOPS Test 5: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 5: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
    return 0;
}


/**
 * zero Upops
 */
int test6()
{

	UpopsZ_t* zero = zero_UPOPSZ();
	PolyZ_ptr mainVarPoly = zero->data[0]->polys[0];
	if (isZero_AAZ( mainVarPoly)) {
	    fprintf(stderr, "UPOPS Test 6: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 6: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
    return 0;
}


/**
 * one Upops
 */
int test7()
{

	UpopsZ_t* one = one_UPOPSZ(4);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly = one->data[0]->polys[0];
	PolyZ_ptr subpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "UPOPS Test 7: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 7: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(testpoly);
	freePolynomial_AAZ(subpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
    return 0;
}


/**
 * Upops addition example 1 : addition of two zero Upops
 */
int test8()
{

	UpopsZ_t* zero = zero_UPOPSZ();
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (zero, zero);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly = sum->data[0]->polys[0];
	PolyZ_ptr subpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "UPOPS Test 8: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 8: \t\t\t\t\t FAILED \n");
		exit(1);
    }

    freePolynomial_AAZ(subpoly);
    freePolynomial_AAZ(testpoly);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);
	return 0;
}


/**
 * Upops addition example 2 : addition of two one Upops
 */
int test9()
{

	UpopsZ_t* one = one_UPOPSZ(4);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (one, one);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("2", syms, 4);
	PolyZ_ptr mainVarPoly =  sum->data[0]->polys[0];
	PolyZ_ptr subpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "UPOPS Test 9: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 9: \t\t\t\t\t FAILED \n");
		exit(1);
    }

    freePolynomial_AAZ(subpoly);
    freePolynomial_AAZ(testpoly);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
    destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);
	return 0;
}


/**
 * Upops addition example 3 : addition of zero Upops and one Upops
 */
int test10()
{
	UpopsZ_t* zero = zero_UPOPSZ();
	UpopsZ_t* one = one_UPOPSZ(4);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (zero, one);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly = sum->data[0]->polys[0];
	PolyZ_ptr subpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "UPOPS Test 10: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 10: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly);
	freePolynomial_AAZ(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);
    return 0;
}

/**
 * Upops addition example 4 : addition of zero Upops and non-constant Upops
 */
int test11()
{
	UpopsZ_t*  zero = zero_UPOPSZ();
	PolyZ_ptr poly = generate_altarrZ_var_defined("4+3*z*x^2+y*t*x^5+y", syms, 4);
	UpopsZ_t*  upops = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (zero, upops);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("3*z", syms, 4);
	PolyZ_ptr mainVarPoly_1 = sum->data[2]->polys[1];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("y*t", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sum->data[5]->polys[2];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("4", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sum->data[0]->polys[0];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("y", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sum->data[0]->polys[1];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 11: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 11: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(poly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);


    return 0;
}



/**
 * Upops addition example 5 : addition of zero Upops and non-constant Upops
 */
int test12()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("z", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+y", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("y*t", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(2,  sum->data[6]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sum->data[6]->polys[0];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("3*z", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sum->data[2]->polys[1];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("y+z", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sum->data[0]->polys[1];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 12: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 12: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);
    return 0;
}

/**
 * Upops addition example 6 : addition of two Upops
 */
int test13()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t^6)*x^6+y+x", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("z+2*y*x^2", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("y+z", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(1,  sum->data[0]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sum->data[1]->polys[0];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 =  NULL;
	PolyZ_ptr mainVarPoly_3 = sum->data[3]->polys[0];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("2*y+3*z", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sum->data[2]->polys[1];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);
	PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("y*t^6", syms, 4);
	PolyZ_ptr mainVarPoly_5 = sum->data[6]->polys[7];
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,4);
	PolyZ_ptr testpoly_6 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_6 = sum->data[6]->polys[0];
	PolyZ_ptr subpoly_6 = subPolynomials_AAZ(mainVarPoly_6,testpoly_6,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5) && isZero_AAZ(subpoly_6)) {
	    fprintf(stderr, "UPOPS Test 13: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 13: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(subpoly_6);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(testpoly_6);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);

    return 0;
}


/**
 * Upops addition example 7 : addition of two Upops
 */
int test14()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+y", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("z^6", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(5,  sum->data[0]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("y", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sum->data[0]->polys[1];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("3*z", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sum->data[2]->polys[1];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("y*t", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sum->data[6]->polys[2];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);
	PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_5 = sum->data[6]->polys[0];
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5)) {
	    fprintf(stderr, "UPOPS Test 14: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 14: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);

    return 0;
}

/**
 * Upops subtraction example 1
 */
int test15()
{

	UpopsZ_t* zero = zero_UPOPSZ();
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (zero, zero);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly = sub->data[0]->polys[0];
	PolyZ_ptr subpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "UPOPS Test 15: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 15: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly);
	freePolynomial_AAZ(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);

    return 0;
}




/**
 * Upops subtraction example 2
 */
int test16()
{

	UpopsZ_t* one = one_UPOPSZ(4);
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (one, one);
	updateToDeg_UPOPSZ(3, sub);
	PolyZ_ptr mainVarPoly =  sub->data[0]->polys[0];
	if (isZero_AAZ(mainVarPoly)) {
	    fprintf(stderr, "UPOPS Test 16: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 16: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
    return 0;
}


/**
 * Upops subtraction example 3
 */
int test17()
{
	UpopsZ_t* zero = zero_UPOPSZ();
	UpopsZ_t* one = one_UPOPSZ(4);
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (zero, one);
	updateToDeg_UPOPSZ(4, sub);
	PolyZ_ptr mainVarPoly = sub->data[0]->polys[0];
	if (isNegativeOne_AAZ(mainVarPoly)) {
	    fprintf(stderr, "UPOPS Test 17: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 17: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
    return 0;
}

/**
 * Upops subtraction example 4
 */
int test18()
{
	UpopsZ_t* zero = zero_UPOPSZ();
	PolyZ_ptr poly = generate_altarrZ_var_defined("4+3*z*x^2+y*t*x^4+y", syms, 4);
	UpopsZ_t*  upops = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly);
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (zero, upops);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("-3*z", syms, 4);
	PolyZ_ptr mainVarPoly_1 = sub->data[2]->polys[1];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("-4", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sub->data[0]->polys[0];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("-y", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sub->data[0]->polys[1];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("-y*t", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sub->data[4]->polys[2];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 18: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 18: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(poly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);
    return 0;
}



/**
 * Upops subtraction example 5: issue with printPolyClean
 */
int test19()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("z", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("3*z*x^2+(1+y^3*t)*x^6+y", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("-y^3*t", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(4,  sub->data[6]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("z-y", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sub->data[0]->polys[1];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("-3*z", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sub->data[2]->polys[1];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("-1", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sub->data[6]->polys[0];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);


	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 19: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 19: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);
    return 0;
}

/**
 * Upops subtraction example 6
 */
int test20()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+y", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("z", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("y-z", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(1,  sub->data[0]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("y*t", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sub->data[6]->polys[2];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sub->data[6]->polys[0];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("3*z", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sub->data[2]->polys[1];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 20: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 20: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);

    return 0;
}


/**
 * Upops subtraction example 7
 */
int test21()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+t^8", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("z^6", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* sub = subUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(5,  sub->data[0]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("-z^6", syms, 4);
	PolyZ_ptr mainVarPoly_2 = sub->data[0]->polys[6];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("y*t", syms, 4);
	PolyZ_ptr mainVarPoly_3 = sub->data[6]->polys[2];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_4 = sub->data[6]->polys[0];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);
	PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("3*z", syms, 4);
	PolyZ_ptr mainVarPoly_5 = sub->data[2]->polys[1];
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5)) {
	    fprintf(stderr, "UPOPS Test 21: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 21: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sub);

    return 0;
}


/**
 * Upops multiplication example 1
 */
int test22()
{

	UpopsZ_t* zero = zero_UPOPSZ();
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (zero, zero);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly = prod->data[0]->polys[0];
	PolyZ_ptr addpoly = addPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(addpoly)) {
	    fprintf(stderr, "UPOPS Test 22: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 22: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(testpoly);
	freePolynomial_AAZ(addpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
	return 0;

}

/**
 * Upops multiplication example 2
 */
int test23()
{

	UpopsZ_t* one = one_UPOPSZ(4);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (one, one);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly =  prod->data[0]->polys[0];
	PolyZ_ptr subpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "UPOPS Test 23: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 23: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly);
	freePolynomial_AAZ(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    return 0;
}


/**
 * Upops multiplication  example 3
 */
int test24()
{
	UpopsZ_t* zero = zero_UPOPSZ();
	UpopsZ_t* one = one_UPOPSZ(4);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (zero, one);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly = prod->data[0]->polys[0];
	PolyZ_ptr prodpoly = subPolynomials_AAZ(mainVarPoly,testpoly,4);

	if (isZero_AAZ(prodpoly)) {
	    fprintf(stderr, "UPOPS Test 24: \t\t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "UPOPS Test 24: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(prodpoly);
	freePolynomial_AAZ(testpoly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(one);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);
    return 0;
}


/**
 * Upops multiplication  example 4
 */
int test25()
{
	UpopsZ_t* zero = zero_UPOPSZ();
	PolyZ_ptr poly = generate_altarrZ_var_defined("4+3*z*x^2+y*t*x^4+y", syms, 4);
	UpopsZ_t*  upops = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (zero, upops);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("0", syms, 4);
	PolyZ_ptr mainVarPoly_1 = prod->data[2]->polys[0];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr mainVarPoly_2 = prod->data[0]->polys[0];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_1,4);
	PolyZ_ptr mainVarPoly_3 = prod->data[1]->polys[0];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_1,4);
	PolyZ_ptr mainVarPoly_4 = prod->data[3]->polys[0];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_1,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 25: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 25: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(poly);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(zero);

    return 0;
}


/**
 * Upops multiplication example 5
 */
int test26()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("z^2+x^2", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+y^3+1", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("y^3+3*z^3", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(3,  prod->data[2]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("z^2", syms, 4);
	PolyZ_ptr mainVarPoly_2 = prod->data[0]->polys[2];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = NULL;
	PolyZ_ptr mainVarPoly_3 = prod->data[1]->polys[0];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr mainVarPoly_4 = prod->data[2]->polys[0];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,4);
	PolyZ_ptr mainVarPoly_5 = prod->data[3]->polys[0];
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_3,4);
	PolyZ_ptr mainVarPoly_6 = prod->data[4]->polys[0];
	PolyZ_ptr subpoly_6 = subPolynomials_AAZ(mainVarPoly_6,testpoly_3,4);
	PolyZ_ptr mainVarPoly_7 = prod->data[5]->polys[0];
	PolyZ_ptr subpoly_7 = subPolynomials_AAZ(mainVarPoly_7,testpoly_3,4);
	PolyZ_ptr mainVarPoly_8 = prod->data[8]->polys[0];
	PolyZ_ptr subpoly_8 = subPolynomials_AAZ
	(mainVarPoly_8,testpoly_4,4);
	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5) && isZero_AAZ(subpoly_6) && isZero_AAZ(subpoly_7) && isZero_AAZ(subpoly_8)) {
	    fprintf(stderr, "UPOPS Test 26: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 26: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(subpoly_6);
	freePolynomial_AAZ(subpoly_7);
	freePolynomial_AAZ(subpoly_8);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);

	return 0;

}


/**
 * Upops multiplication example 6
 */
int test27()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+y+1", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("t", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("t*y", syms, 4);
	PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(2,  prod->data[0]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("t", syms, 4);
	PolyZ_ptr mainVarPoly_2 = prod->data[0]->polys[1];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("3*t*z", syms, 4);
	PolyZ_ptr mainVarPoly_3 = homogPart_PSZ(2,  prod->data[2]);
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);
	PolyZ_ptr mainVarPoly_4 = prod->data[6]->polys[1];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_2,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)) {
	    fprintf(stderr, "UPOPS Test 27: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 27: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    return 0;
}


/**
 * Upops multiplication example 7
 */
int test28()
{
	PolyZ_ptr poly1 = generate_altarrZ_var_defined("3*z*x^2+(1+y*t)*x^6+t+z^3", syms, 4);
	PolyZ_ptr poly2 = generate_altarrZ_var_defined("z^6+t^2", syms, 4);
	UpopsZ_t*  upops1 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly1);
	UpopsZ_t*  upops2 = convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(poly2);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (upops1, upops2);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("t^3", syms, 4);
	PolyZ_ptr mainVarPoly_1 = prod->data[0]->polys[3];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,4);
	PolyZ_ptr testpoly_2 = NULL;
	PolyZ_ptr mainVarPoly_2 = prod->data[2]->polys[0];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,4);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("t^2", syms, 4);
	PolyZ_ptr mainVarPoly_3 = prod->data[6]->polys[2];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,4);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3)) {
	    fprintf(stderr, "UPOPS Test 28: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 28: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(poly1);
	freePolynomial_AAZ(poly2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upops1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    return 0;
}

/**
 *  Upops addition example 8: using power series as coefficients of the upops
 */
int test29()
{
	PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", Vars, 3);
	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1-x", Vars, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y", Vars, 3);
	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
	homogPart_PSZ(2, quo_1);
	PowerSeriesZ_t* quo_2 = dividePowerSeries_PSZ(ps_1, ps_3);
	homogPart_PSZ(2, quo_2);

	PowerSeriesZ_t* ps_array[2];
	ps_array[0] = quo_1;
	ps_array[1] = quo_2;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);
	UpopsZ_t* sum = addUnivariatePolyOverPowerSeries_UPOPSZ (upop, upop);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("2", Vars, 3);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("2*x", Vars, 3);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("2*x^2", Vars, 3);
	PolyZ_ptr mainVarPoly_1 = sum->data[0]->polys[0];
	PolyZ_ptr mainVarPoly_2 = sum->data[0]->polys[1];
	PolyZ_ptr mainVarPoly_3 = sum->data[0]->polys[2];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,3);
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,3);
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,3);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("2", Vars, 3);
	PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("2*y", Vars, 3);
	PolyZ_ptr testpoly_6 = generate_altarrZ_var_defined("2*y^2", Vars, 3);
	PolyZ_ptr mainVarPoly_4 = sum->data[1]->polys[0];
	PolyZ_ptr mainVarPoly_5 = sum->data[1]->polys[1];
	PolyZ_ptr mainVarPoly_6 = sum->data[1]->polys[2];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,3);
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,3);
	PolyZ_ptr subpoly_6 = subPolynomials_AAZ(mainVarPoly_6,testpoly_6,3);



	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5) && isZero_AAZ(subpoly_6)) {
	    fprintf(stderr, "UPOPS Test 29: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 29: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(subpoly_6);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(testpoly_6);
	freePolynomial_AAZ(poly_1);
	freePolynomial_AAZ(poly_2);
	freePolynomial_AAZ(poly_3);
	destroyPowerSeries_PSZ(ps_1);
	destroyPowerSeries_PSZ(ps_2);
	destroyPowerSeries_PSZ(ps_3);
	destroyPowerSeries_PSZ(quo_1);
	destroyPowerSeries_PSZ(quo_2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(sum);
    return 0;

}



/**
 *  Upops multiplication example 8: using power series as coefficients of the upops
 */
int test30()
{
	PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", Vars, 3);
	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("1-x", Vars, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y", Vars, 3);
	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_2);
	homogPart_PSZ(2, quo_1);
	PowerSeriesZ_t* quo_2 = dividePowerSeries_PSZ(ps_1, ps_3);
	homogPart_PSZ(2, quo_2);
	PowerSeriesZ_t* ps_array[2];
	ps_array[0] = quo_1;
	ps_array[1] = quo_2;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (upop, upop);
	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("1", Vars, 3);
	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("2*x", Vars, 3);
	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("3*x^2", Vars, 3);
	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("4*x^3", Vars, 3);
	PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("5*x^4", Vars, 3);
	PolyZ_ptr mainVarPoly_1 = prod->data[0]->polys[0];
	PolyZ_ptr mainVarPoly_2 = prod->data[0]->polys[1];
	PolyZ_ptr mainVarPoly_3 = prod->data[0]->polys[2];
	PolyZ_ptr mainVarPoly_4 = homogPart_PSZ(3,  prod->data[0]);
	PolyZ_ptr mainVarPoly_5 = homogPart_PSZ(4,  prod->data[0]);
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,3);
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,3);
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,3);
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,3);
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,3);

	PolyZ_ptr testpoly_6 = generate_altarrZ_var_defined("2", Vars, 3);
	PolyZ_ptr testpoly_7 = generate_altarrZ_var_defined("2*x+2*y", Vars, 3);
	PolyZ_ptr testpoly_8 = generate_altarrZ_var_defined("2*x^2+2*x*y+2*y^2", Vars, 3);
	PolyZ_ptr testpoly_9 = generate_altarrZ_var_defined("2*x^3+2*y^3+2*x^2*y+2*x*y^2", Vars, 3);
	PolyZ_ptr testpoly_10 = generate_altarrZ_var_defined("2*x^4+2*y^4+2*x*y^3+2*x^3*y+2*x^2*y^2", Vars, 3);
	PolyZ_ptr mainVarPoly_6 = prod->data[1]->polys[0];
	PolyZ_ptr mainVarPoly_7 = prod->data[1]->polys[1];
	PolyZ_ptr mainVarPoly_8 = prod->data[1]->polys[2];
	PolyZ_ptr mainVarPoly_9 = homogPart_PSZ(3,  prod->data[1]);
	PolyZ_ptr mainVarPoly_10 = homogPart_PSZ(4,  prod->data[1]);
	PolyZ_ptr subpoly_6 = subPolynomials_AAZ(mainVarPoly_6,testpoly_6,3);
	PolyZ_ptr subpoly_7 = subPolynomials_AAZ(mainVarPoly_7,testpoly_7,3);
	PolyZ_ptr subpoly_8 = subPolynomials_AAZ(mainVarPoly_8,testpoly_8,3);
	PolyZ_ptr subpoly_9 = subPolynomials_AAZ(mainVarPoly_9,testpoly_9,3);
	PolyZ_ptr subpoly_10 = subPolynomials_AAZ(mainVarPoly_10,testpoly_10,3);


	PolyZ_ptr testpoly_11 = generate_altarrZ_var_defined("1", Vars, 3);
	PolyZ_ptr testpoly_12 = generate_altarrZ_var_defined("2*y", Vars, 3);
	PolyZ_ptr testpoly_13 = generate_altarrZ_var_defined("3*y^2", Vars, 3);
	PolyZ_ptr testpoly_14 = generate_altarrZ_var_defined("4*y^3", Vars, 3);
	PolyZ_ptr testpoly_15 = generate_altarrZ_var_defined("5*y^4", Vars, 3);
	PolyZ_ptr mainVarPoly_11 = prod->data[2]->polys[0];
	PolyZ_ptr mainVarPoly_12 = prod->data[2]->polys[1];
	PolyZ_ptr mainVarPoly_13 = prod->data[2]->polys[2];
	PolyZ_ptr mainVarPoly_14 = homogPart_PSZ(3,  prod->data[2]);
	PolyZ_ptr mainVarPoly_15 = homogPart_PSZ(4,  prod->data[2]);
	PolyZ_ptr subpoly_11 = subPolynomials_AAZ(mainVarPoly_11,testpoly_11,3);
	PolyZ_ptr subpoly_12 = subPolynomials_AAZ(mainVarPoly_12,testpoly_12,3);
	PolyZ_ptr subpoly_13 = subPolynomials_AAZ(mainVarPoly_13,testpoly_13,3);
	PolyZ_ptr subpoly_14 = subPolynomials_AAZ(mainVarPoly_14,testpoly_14,3);
	PolyZ_ptr subpoly_15 = subPolynomials_AAZ(mainVarPoly_15,testpoly_15,3);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5) && isZero_AAZ(subpoly_6) && isZero_AAZ(subpoly_7) && isZero_AAZ(subpoly_8)
		 && isZero_AAZ(subpoly_9) && isZero_AAZ(subpoly_10) && isZero_AAZ(subpoly_11) && isZero_AAZ(subpoly_12)
		  && isZero_AAZ(subpoly_13) && isZero_AAZ(subpoly_14) && isZero_AAZ(subpoly_15)) {
	    fprintf(stderr, "UPOPS Test 30: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 30: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(subpoly_6);
	freePolynomial_AAZ(subpoly_7);
	freePolynomial_AAZ(subpoly_8);
	freePolynomial_AAZ(subpoly_9);
	freePolynomial_AAZ(subpoly_10);
	freePolynomial_AAZ(subpoly_11);
	freePolynomial_AAZ(subpoly_12);
	freePolynomial_AAZ(subpoly_13);
	freePolynomial_AAZ(subpoly_14);
	freePolynomial_AAZ(subpoly_15);

	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(testpoly_6);
	freePolynomial_AAZ(testpoly_7);
	freePolynomial_AAZ(testpoly_8);
	freePolynomial_AAZ(testpoly_9);
	freePolynomial_AAZ(testpoly_10);
	freePolynomial_AAZ(testpoly_11);
	freePolynomial_AAZ(testpoly_12);
	freePolynomial_AAZ(testpoly_13);
	freePolynomial_AAZ(testpoly_14);
	freePolynomial_AAZ(testpoly_15);

	freePolynomial_AAZ(poly_1);
	freePolynomial_AAZ(poly_2);
	freePolynomial_AAZ(poly_3);
	destroyPowerSeries_PSZ(ps_1);
	destroyPowerSeries_PSZ(ps_2);
	destroyPowerSeries_PSZ(ps_3);
	destroyPowerSeries_PSZ(quo_1);
	destroyPowerSeries_PSZ(quo_2);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    return 0;

}


/**
 *  Upops multiplication example 9: using power series as coefficients of the upops
 */
int test31()

{

	PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", Vars, 3);
	PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("x*y+x^4", Vars, 3);
	PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y^2+x^2", Vars, 3);
	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_3);
	homogPart_PSZ(6, quo_1);
	PowerSeriesZ_t* ps_array[2];
	ps_array[0] = ps_2;
	ps_array[1] = quo_1;
	UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,2);
	UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (upop, upop);

	PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("x^2*y^2", Vars, 3);
	PolyZ_ptr mainVarPoly_1 = prod->data[0]->polys[4];
	PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,3);

	PolyZ_ptr testpoly_2 = generate_altarrZ_var_defined("2*x^4+2*x*y^3-2*x^3*y", Vars, 3);
	PolyZ_ptr mainVarPoly_2 = prod->data[1]->polys[4];
	PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,3);

	PolyZ_ptr testpoly_3 = generate_altarrZ_var_defined("2*x*y", Vars, 3);
	PolyZ_ptr mainVarPoly_3 = prod->data[1]->polys[2];
	PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,3);

	PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("1", Vars, 3);
	PolyZ_ptr mainVarPoly_4 = prod->data[2]->polys[0];
	PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,3);

	PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("2*y^2-2*x^2", Vars, 3);
	PolyZ_ptr mainVarPoly_5 = prod->data[2]->polys[2];
	PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,3);

	PolyZ_ptr testpoly_6 = generate_altarrZ_var_defined("3*y^4-6*x^2*y^2+3*x^4", Vars, 3);
	PolyZ_ptr mainVarPoly_6 = prod->data[2]->polys[4];
	PolyZ_ptr subpoly_6 = subPolynomials_AAZ(mainVarPoly_6,testpoly_6,3);

	PolyZ_ptr testpoly_7 = generate_altarrZ_var_defined("-4*x^6+12*x^4*y^2-12*x^2*y^4+4*y^6", Vars, 3);
	PolyZ_ptr mainVarPoly_7 = prod->data[2]->polys[6];
	PolyZ_ptr subpoly_7 = subPolynomials_AAZ(mainVarPoly_7,testpoly_7,3);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5) && isZero_AAZ(subpoly_6) && isZero_AAZ(subpoly_7) ) {
	    fprintf(stderr, "UPOPS Test 31: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 31: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(subpoly_6);
	freePolynomial_AAZ(subpoly_7);

	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(testpoly_6);
	freePolynomial_AAZ(testpoly_7);

	freePolynomial_AAZ(poly_1);
	freePolynomial_AAZ(poly_2);
	freePolynomial_AAZ(poly_3);
	destroyPowerSeries_PSZ(ps_1);
	destroyPowerSeries_PSZ(ps_2);
	destroyPowerSeries_PSZ(ps_3);
	destroyPowerSeries_PSZ(quo_1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    return 0;

}


/**
 *  Upops multiplication example 10: using power series as coefficients of the upops
 */
int test32()

{
	 PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("x+1", Vars, 3);
	 PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("x^2+y^2", Vars, 3);
	 PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y^2+2*x^3", Vars, 3);
	 PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("x", Vars, 3);
	 PolyZ_ptr  poly_5 = generate_altarrZ_var_defined("x+y", Vars, 3);
	 PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	 PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	 PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	 PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_3);
	 homogPart_PSZ(4, quo_1);
	 PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
	 PowerSeriesZ_t* ps_5 = convertPolyToPowerSeries_PSZ(poly_5);
	 PowerSeriesZ_t* ps_array[4];
	 ps_array[0] = ps_2;
	 ps_array[1] =  ps_4;
	 ps_array[2] = ps_5;
	 ps_array[3] = quo_1;

	 UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,4);
	 UpopsZ_t* prod = multiplyUnivariatePolyOverPowerSeries_UPOPSZ (upop, upop);

	 PolyZ_ptr testpoly_1 = generate_altarrZ_var_defined("2*x^2*y^2+x^4+y^4", Vars, 3);
	 PolyZ_ptr mainVarPoly_1 = homogPart_PSZ(4, prod->data[0]);
	 PolyZ_ptr subpoly_1 = subPolynomials_AAZ(mainVarPoly_1,testpoly_1,3);

	 PolyZ_ptr testpoly_2 = NULL;
	 PolyZ_ptr mainVarPoly_2 = prod->data[1]->polys[0];
	 PolyZ_ptr subpoly_2 = subPolynomials_AAZ(mainVarPoly_2,testpoly_2,3);


	 PolyZ_ptr testpoly_3 =  NULL;
	 PolyZ_ptr mainVarPoly_3 = prod->data[2]->polys[0];
	 PolyZ_ptr subpoly_3 = subPolynomials_AAZ(mainVarPoly_3,testpoly_3,3);

	 PolyZ_ptr testpoly_4 = generate_altarrZ_var_defined("4*x^2+2*y^2+2*x*y", Vars, 3);;
	 PolyZ_ptr mainVarPoly_4 = homogPart_PSZ(2,prod->data[3]);
	 PolyZ_ptr subpoly_4 = subPolynomials_AAZ(mainVarPoly_4,testpoly_4,3);

	 PolyZ_ptr testpoly_5 = generate_altarrZ_var_defined("2*x", Vars, 3);;
	 PolyZ_ptr mainVarPoly_5 = prod->data[4]->polys[1];
	 PolyZ_ptr subpoly_5 = subPolynomials_AAZ(mainVarPoly_5,testpoly_5,3);


	 PolyZ_ptr testpoly_6 = generate_altarrZ_var_defined("2*x+2*y", Vars, 3);;
	 PolyZ_ptr mainVarPoly_6 = prod->data[5]->polys[1];
	 PolyZ_ptr subpoly_6 = subPolynomials_AAZ(mainVarPoly_6,testpoly_6,3);

	 PolyZ_ptr testpoly_7 = generate_altarrZ_var_defined("1", Vars, 3);;
	 PolyZ_ptr mainVarPoly_7 = prod->data[6]->polys[0];
	 PolyZ_ptr subpoly_7 = subPolynomials_AAZ(mainVarPoly_7,testpoly_7,3);

	 PolyZ_ptr testpoly_8 = generate_altarrZ_var_defined("2*x", Vars, 3);;
	 PolyZ_ptr mainVarPoly_8 = prod->data[6]->polys[1];
	 PolyZ_ptr subpoly_8 = subPolynomials_AAZ(mainVarPoly_8,testpoly_8,3);

	 PolyZ_ptr testpoly_9 = generate_altarrZ_var_defined("2*y^2+x^2", Vars, 3);;
	 PolyZ_ptr mainVarPoly_9 = prod->data[6]->polys[2];
	 PolyZ_ptr subpoly_9 = subPolynomials_AAZ(mainVarPoly_9,testpoly_9,3);

	 PolyZ_ptr testpoly_10 = generate_altarrZ_var_defined("4*x*y^2-4*x^3", Vars, 3);;
	 PolyZ_ptr mainVarPoly_10 = prod->data[6]->polys[3];
	 PolyZ_ptr subpoly_10 = subPolynomials_AAZ(mainVarPoly_10,testpoly_10,3);

	 PolyZ_ptr testpoly_11 = generate_altarrZ_var_defined("3*y^4+2*x^2*y^2-8*x^4", Vars, 3);;
	 PolyZ_ptr mainVarPoly_11 = prod->data[6]->polys[4];
	 PolyZ_ptr subpoly_11 = subPolynomials_AAZ(mainVarPoly_11,testpoly_11,3);

	if (isZero_AAZ(subpoly_1) && isZero_AAZ(subpoly_2) && isZero_AAZ(subpoly_3) && isZero_AAZ(subpoly_4)
		&& isZero_AAZ(subpoly_5) && isZero_AAZ(subpoly_6) && isZero_AAZ(subpoly_7) && isZero_AAZ(subpoly_8)
		 && isZero_AAZ(subpoly_9) && isZero_AAZ(subpoly_10) && isZero_AAZ(subpoly_11)) {
	    fprintf(stderr, "UPOPS Test 32: \t\t\t\t\t PASSED \n");
	} else {
	    fprintf(stderr, "UPOPS Test 32: \t\t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(subpoly_1);
	freePolynomial_AAZ(subpoly_2);
	freePolynomial_AAZ(subpoly_3);
	freePolynomial_AAZ(subpoly_4);
	freePolynomial_AAZ(subpoly_5);
	freePolynomial_AAZ(subpoly_6);
	freePolynomial_AAZ(subpoly_7);
	freePolynomial_AAZ(subpoly_8);
	freePolynomial_AAZ(subpoly_9);
	freePolynomial_AAZ(subpoly_10);
	freePolynomial_AAZ(subpoly_11);

	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);
	freePolynomial_AAZ(testpoly_3);
	freePolynomial_AAZ(testpoly_4);
	freePolynomial_AAZ(testpoly_5);
	freePolynomial_AAZ(testpoly_6);
	freePolynomial_AAZ(testpoly_7);
	freePolynomial_AAZ(testpoly_8);
	freePolynomial_AAZ(testpoly_9);
	freePolynomial_AAZ(testpoly_10);
	freePolynomial_AAZ(testpoly_11);

	freePolynomial_AAZ(poly_1);
	freePolynomial_AAZ(poly_2);
	freePolynomial_AAZ(poly_3);
	freePolynomial_AAZ(poly_4);
	freePolynomial_AAZ(poly_5);
	destroyPowerSeries_PSZ(ps_1);
	destroyPowerSeries_PSZ(ps_2);
	destroyPowerSeries_PSZ(ps_3);
	destroyPowerSeries_PSZ(ps_4);
	destroyPowerSeries_PSZ(ps_5);
	destroyPowerSeries_PSZ(quo_1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(prod);
    return 0;
}


/**
 *
 */
int test33()
{
	 PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("x+1", Vars, 3);
	 PolyZ_ptr  poly_2 = generate_altarrZ_var_defined("x^2+y^2", Vars, 3);
	 PolyZ_ptr  poly_3 = generate_altarrZ_var_defined("1-y^2+2*x^3", Vars, 3);
	 PolyZ_ptr  poly_4 = generate_altarrZ_var_defined("x", Vars, 3);
	 PolyZ_ptr  poly_5 = generate_altarrZ_var_defined("x+y", Vars, 3);
	 PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
	 PowerSeriesZ_t* ps_2 = convertPolyToPowerSeries_PSZ(poly_2);
	 PowerSeriesZ_t* ps_3 = convertPolyToPowerSeries_PSZ(poly_3);
	 PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, ps_3);
	 homogPart_PSZ(4, quo_1);
	 PowerSeriesZ_t* ps_4 = convertPolyToPowerSeries_PSZ(poly_4);
	 PowerSeriesZ_t* ps_5 = convertPolyToPowerSeries_PSZ(poly_5);
	 PowerSeriesZ_t* ps_array[4];
	 ps_array[0] = ps_2;
	 ps_array[1] =  ps_4;
	 ps_array[2] = ps_5;
	 ps_array[3] = quo_1;
	 UpopsZ_t* upop = convertArrayOfPSToUPOPS_UPOPSZ(ps_array,4);
	 // print_UPOPSZ(stderr,  upop, Vars);
	 // fprintf(stderr, "\n");
	 // fprintf(stderr, "Up to degree 7 \n");
	 updateToDeg_UPOPSZ(7, upop);
	 // print_UPOPSZ(stderr,  upop, Vars);
	 // fprintf(stderr, "\n");
	 // fprintf(stderr, "Up to degree 17 \n");
	 updateToDeg_UPOPSZ(17, upop);
	 // print_UPOPSZ(stderr,  upop, Vars);


	freePolynomial_AAZ(poly_1);
	freePolynomial_AAZ(poly_2);
	freePolynomial_AAZ(poly_3);
	freePolynomial_AAZ(poly_4);
	freePolynomial_AAZ(poly_5);
	destroyPowerSeries_PSZ(ps_1);
	destroyPowerSeries_PSZ(ps_2);
	destroyPowerSeries_PSZ(ps_3);
	destroyPowerSeries_PSZ(ps_4);
	destroyPowerSeries_PSZ(ps_5);
	destroyPowerSeries_PSZ(quo_1);
	destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(upop);

	//Same as previous example, just with higher degree
    fprintf(stderr, "UPOPS Test 33: \t\t\t\t\t PASSED \n");

	return 0;
}


// int test34()
// {
// 	PowerSeriesZ_t* a = geometricSeries_PSZ(3);
//     //PolyZ_ptr s = homogPart_PSZ(9,a);
// 	PolyZ_ptr  poly_1 = generate_altarrZ_var_defined("1", vars, 3);
// 	PowerSeriesZ_t* ps_1 = convertPolyToPowerSeries_PSZ(poly_1);
// 	PowerSeriesZ_t* quo_1 = dividePowerSeries_PSZ(ps_1, a);
// 	PolyZ_ptr s = homogPart_PSZ(1,quo_1);
// 	//print_PSZ(stderr, a,sym);
// 	//printPoly_AAZ(stderr, s, vars,3);
// 	//fprintf(stderr, "\n");
//     return 0;

// }
