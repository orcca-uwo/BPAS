#include "../../include/PowerSeries/PowerSeriesZ.h"
#include "../../include/Utils/Unix_Timer.h"
#include "../../include/Parser/bpas_parser.h"
#include <time.h>


const char* sym[] = {"x", "y"};
const char* syms[] = {"x", "y", "z", "t"};
const char* vars[] = {"x", "y", "z"};



void TestingPowerSeriesDivision(PowerSeriesZ_t* f, PowerSeriesZ_t* g, int bound, const char** sym, int n) {

	PowerSeriesZ_t* h = dividePowerSeries_PSZ(f,g);
	PowerSeriesZ_t* prod = multiplyPowerSeries_PSZ(g,h);
	PowerSeriesZ_t* sub = subPowerSeries_PSZ(f,prod);

	for (int d = 1; d <= bound; ++d) {
		//fprintf(stderr, "In degree: %d\n ", d);
		//fprintf(stderr, "\n");
		//fprintf(stderr, "polyPart division of degree d = %d\n", d);
	    //fprintf(stderr, "\n");
		//printPolyClean_AAZ(stderr,polynomialPart_PSZ(d, h), sym, n);
		//fprintf(stderr, "\n");
		PolyZ_ptr polyPart = polynomialPart_PSZ(d,  sub);

		if (!isZero_AAZ(polyPart)) {
			 fprintf(stderr, "ERROR IN TEST POWER SERIES DIVISION\n");
			 exit(1);
		}

	}

	destroyPowerSeries_PSZ(h);
	destroyPowerSeries_PSZ(prod);
	destroyPowerSeries_PSZ(sub);
}

void TestingPowerSeriesInversion(PowerSeriesZ_t* f, int bound, const char** sym, int n) {

	PowerSeriesZ_t* h = inversePowerSeries_PSZ(f);
	PowerSeriesZ_t* prod = multiplyPowerSeries_PSZ(f, h);
	for (int d = 1; d <= bound; ++d) {
		//fprintf(stderr, "In degree: %d\n ", d);
		//fprintf(stderr, "\n");
		//fprintf(stderr, "polyPart division of degree d = %d\n", d);
	    //fprintf(stderr, "\n");
		//printPolyClean_AAZ(stderr,polynomialPart_PSZ(d, h), sym, n);
		//fprintf(stderr, "\n");
		PolyZ_ptr polyPart = polynomialPart_PSZ(d,  prod);

		if (!isOne_AAZ(polyPart)) {
			 fprintf(stderr, "ERROR IN TEST POWER SERIES DIVISION\n");
			 exit(1);
		}

	}

	destroyPowerSeries_PSZ(h);
	destroyPowerSeries_PSZ(prod);
}




int test1()
{
	PolyZ_ptr  polynom = generate_altarrZ_var_defined("x^10+1", syms, 4);
	PowerSeriesZ_t*  a = convertPolyToPowerSeries_PSZ(polynom);
	PowerSeriesZ_t* b = multiplyPowerSeries_PSZ(a,a);
	PolyZ_ptr homog = homogPart_PSZ(20, b);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("x^20", syms, 4);
	testpoly = subPolynomials_AAZ(testpoly,homog,4);

	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 1: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 1: \t\t\t\t FAILED \n");
        exit(1);
	}
	freePolynomial_AAZ(polynom);
	freePolynomial_AAZ(testpoly);

	destroyPowerSeries_PSZ (a);
	destroyPowerSeries_PSZ (b);
    return 0;
}


int test2()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("x+100", syms, 4);
	PolyZ_ptr  b = generate_altarrZ_var_defined("z^5+x^4*t+1", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);
	//printf("%d\n",c->genOrder);
	//printPoly_AAZ (stderr, homogPart_PSZ(0, c), syms, 4);
	PolyZ_ptr homog = homogPart_PSZ(0, c);
	if (isConstant_AAZ(homog))  {
	    fprintf(stderr, "PowerSeriesZ Test 2: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 2: \t\t\t\t FAILED \n");
	    exit(1);
	}

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);

	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	destroyPowerSeries_PSZ (c);
    return 0;
}

int test3()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("x*z*t^4+y+x^2+1", syms, 4);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("x*z*t^4", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PolyZ_ptr subpoly = subPolynomials_AAZ(a_ps->polys[6],testpoly,4);
	//for (int i = 0; i <=a_ps->deg; i++) {
	//		fprintf (stderr, "a_ps->polys[%d]=", i);
	//		printPoly_AAZ (stderr, a_ps->polys[i], syms, 4);
	//		fprintf(stderr, "\n");
	//	}
	if (isZero_AAZ(subpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 3: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 3: \t\t\t\t FAILED \n");
	    exit(1);
	}

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(subpoly);
	freePolynomial_AAZ(testpoly);
	destroyPowerSeries_PSZ (a_ps);
    return 0;
}


int test4()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("1", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PolyZ_ptr  b = generate_altarrZ_var_defined("0", syms, 4);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);
	PolyZ_ptr homog = homogPart_PSZ(3, c);
	if (isZero_AAZ(homog)) {
	    fprintf(stderr, "PowerSeriesZ Test 4: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 4: \t\t\t\t FAILED \n");
	    exit(1);
	}
	//for (int i = 0; i <=a_ps->deg; i++) {
	//		fprintf (stderr, "a_ps->polys[%d]=", i);
	//		printPoly_AAZ (stderr, a_ps->polys[i], syms, 4);
	//		fprintf(stderr, "\n");
	//	}

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	destroyPowerSeries_PSZ (c);
    return 0;
}

int test5()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("1", syms, 4);
	PolyZ_ptr  b = generate_altarrZ_var_defined("2*x^2 + 2*y^2 + 2*x + 2*y +2", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);
	PolyZ_ptr homog = homogPart_PSZ(2, c);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("2*x^2 + 2*y^2", syms, 4);
	testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);
	//for (int i = 0; i <=a_ps->deg; i++) {
	//		fprintf (stderr, "a_ps->polys[%d]=", i);
	//		printPoly_AAZ (stderr, a_ps->polys[i], syms, 4);
	//		fprintf(stderr, "\n");
	//	}
	if (isZero_AAZ(testpoly))  {
	    fprintf(stderr, "PowerSeriesZ Test 5: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 5: \t\t\t\t FAILED \n");
	    exit(1);
	}
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	destroyPowerSeries_PSZ (c);
    return 0;
}



int test6()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("z+1", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,a_ps);
	PolyZ_ptr homog = homogPart_PSZ(2, c);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("z^2", syms, 4);
	testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 6: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 6: \t\t\t\t FAILED \n");
		exit(1);
	}
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(testpoly);
	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (a_ps);
	a = NULL;

    return 0;
}

int test7()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("x+1", syms, 4);
    PolyZ_ptr  b = generate_altarrZ_var_defined("x^2+1", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);
	PolyZ_ptr homog = homogPart_PSZ(3, c);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("x^3", syms, 4);
	testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 7: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 7: \t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly);
	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
    return 0;

}

/**
 * Note that multiplication process truncates the result up to the min degree of ancestors
 * in the following the truncation degree is 3. Thus the array of polys does not contain terms
 * of degree at least 4. But when we call homogPart_PSZ and ask for terms of degree 5 and look at
 * polys again we can see that this time the array polys has been updated and contains terms of degree 5
 */
int test8()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("x^3 + x + z", vars, 3);
	PolyZ_ptr  b = generate_altarrZ_var_defined("z^2*y*x + z*y + 1", vars, 3);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);
	PolyZ_ptr homog = homogPart_PSZ(5, c);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("x^3*z*y+x^2*y*z^2+z^3*y*x", vars, 3);
	testpoly = subPolynomials_AAZ_inp(testpoly,homog,3);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 8: \t\t\t\t PASSED \n");
	}
	else  {
	    fprintf(stderr, "PowerSeriesZ Test 8: \t\t\t\t FAILED \n");
		exit(1);
	}

	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (b_ps);
	destroyPowerSeries_PSZ (a_ps);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly);
	c = NULL;
	return 0;
}


int test9()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("x^2*t+z^7*y*x-2*x^4+1", syms, 4);
	PolyZ_ptr  b = generate_altarrZ_var_defined("x^20+x*t*y*z+z^4+2", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);

	PolyZ_ptr d = homogPart_PSZ(8, c);
	if (isZero_AAZ(d)|| isZero_AAZ(c->polys[8])) {
	    fprintf(stderr, "PowerSeriesZ Test 9: \t\t\t\t FAILED \n");
		exit(1);
	} else {
	    fprintf(stderr, "PowerSeriesZ Test 9: \t\t\t\t PASSED \n");
	}

	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);

	return 0;
}


int test10()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("x^2*t^3+z^7*y*x-2*x^4+y^12+x", syms, 4);
	PolyZ_ptr  b = generate_altarrZ_var_defined("x^2+x*t*y*z+z^4*y*x*t+2", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = multiplyPowerSeries_PSZ(a_ps,b_ps);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("-2*x^6", syms, 4);
	PolyZ_ptr homog = homogPart_PSZ(6, c);
	testpoly = subPolynomials_AAZ(testpoly, homog,4);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 10: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 10: \t\t\t\t FAILED \n");
	    exit(1);
	}

	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(testpoly);
	freePolynomial_AAZ(b);
    return 0;
}


int test11()
{
	PolyZ_ptr a = generate_altarrZ_var_defined("-2*x+x^2+1",syms, 4);
	PolyZ_ptr b = generate_altarrZ_var_defined("1",syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = dividePowerSeries_PSZ(b_ps,a_ps);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("5*x^4", syms, 4);
	PolyZ_ptr homog = homogPart_PSZ(4,c);
	testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 11: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 11: \t\t\t\t FAILED \n");
	    exit(1);
	}

	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly);
    return 0;

}


int test12()
{
	PolyZ_ptr a = generate_altarrZ_var_defined("-2*x*t+x^2+y*z+t^4+1",syms, 4);
	PolyZ_ptr b = generate_altarrZ_var_defined("1-t-x",syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* c = dividePowerSeries_PSZ(b_ps,a_ps);
	PolyZ_ptr testpoly = generate_altarrZ_var_defined("t^4-4*t^2*x^2+4*t*x^3+4*t*x*y*z-x^4-2*x^2*y*z-y^2*z^2", syms, 4);
	PolyZ_ptr homog = homogPart_PSZ(4,c);
	testpoly = addPolynomials_AAZ_inp(testpoly, homog,4);
	if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 12: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeriesZ Test 12: \t\t\t\t FAILED \n");
	    exit(1);
	}

    // for (int i = 0; i <= 8; ++i) {
    //     printPoly_AAZ (stderr, homogPart_PSZ(i,c), syms, 4);
    //     fprintf(stderr, "\n");
    // }

	destroyPowerSeries_PSZ (c);
	destroyPowerSeries_PSZ (a_ps);
	destroyPowerSeries_PSZ (b_ps);
	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly);
    return 0;

}


int test13()
{
    PolyZ_ptr  a = generate_altarrZ_var_defined("x^3 + x + z", vars, 3);
    PolyZ_ptr  b = generate_altarrZ_var_defined("z^2*y*x + z*y + 1", vars, 3);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = addPowerSeries_PSZ(a_ps,b_ps);
    PolyZ_ptr homog = homogPart_PSZ(4, c);
    PolyZ_ptr testpoly = generate_altarrZ_var_defined("z^2*y*x", vars, 3);
    testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);
    if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 13: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeriesZ Test 13: \t\t\t\t FAILED \n");
        exit(1);
    }

    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    freePolynomial_AAZ(a);
    freePolynomial_AAZ(b);
    freePolynomial_AAZ(testpoly);
    return 0;
 }



int test14()
{
    PolyZ_ptr  a = generate_altarrZ_var_defined("1+x+y+z+x^2+y^2+z^2", vars, 3);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t* c = addPowerSeries_PSZ(a_ps,a_ps);
    PolyZ_ptr poly = c->polys[0];
    if (isConstant_AAZ(poly)) {
	    fprintf(stderr, "PowerSeriesZ Test 14: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeriesZ Test 14: \t\t\t\t FAILED \n");
        exit(1);
    }

    freePolynomial_AAZ(a);
    destroyPowerSeries_PSZ(a_ps);
    destroyPowerSeries_PSZ(c);
    return 0;
}


int test15()
{
	PolyZ_ptr  a = generate_altarrZ_var_defined("1+x-y+z+2*x^2-4*y^2+z^6", vars, 3);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t* c = subPowerSeries_PSZ(a_ps,a_ps);
    PolyZ_ptr test = polynomialPart_PSZ(6, c);
    if (isZero_AAZ(test)) {
	    fprintf(stderr, "PowerSeriesZ Test 15: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeriesZ Test 15: \t\t\t\t FAILED \n");
        exit(1);
    }

    freePolynomial_AAZ(a);
    freePolynomial_AAZ(test);
    destroyPowerSeries_PSZ(a_ps);
    destroyPowerSeries_PSZ(c);
    return 0;
 }


int test16()
{
    PolyZ_ptr  a = generate_altarrZ_var_defined("x+100", syms, 4);
    PolyZ_ptr  b = generate_altarrZ_var_defined("z^5-y^5+x*t^2*z+x^4*t+1", syms, 4);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = subPowerSeries_PSZ(a_ps,b_ps);
    PolyZ_ptr homog = homogPart_PSZ(5, c);
    PolyZ_ptr testpoly = generate_altarrZ_var_defined("-x^4*t+y^5-z^5", syms, 4);
    testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);
    if (isZero_AAZ(testpoly)) {
	    fprintf(stderr, "PowerSeriesZ Test 16: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeriesZ Test 16: \t\t\t\t FAILED \n");
    	exit(1);
    }

    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    freePolynomial_AAZ(a);
    freePolynomial_AAZ(b);
    freePolynomial_AAZ(testpoly);
    return 0;


}


/**
 *
 */
int test17()
{

	PolyZ_ptr  a = generate_altarrZ_var_defined("x*y^6+z^3-x+100", syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	int x = isUnit_PSZ(a_ps);
	if (x == 0) {
		fprintf(stderr, "PowerSeriesZ Test 17: \t\t\t\t PASSED \n");
		 // printf("It is a unit\n");
	} else {
		fprintf(stderr, "PowerSeriesZ Test 17: \t\t\t\t FAILED \n");
		exit(1);
		 // printf("It is not a unit\n");
	}
	// printf("%d\n", x);

	freePolynomial_AAZ(a);
	destroyPowerSeries_PSZ(a_ps);
	return 0;
}


/**
 *
 */
int test18()
{

	 PolyZ_ptr  a = generate_altarrZ_var_defined("x-y+z^3*t+3*x^2-2*x*y", syms, 4);
	 PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	 int x = isUnit_PSZ(a_ps);
	 if (!x) {
		fprintf(stderr, "PowerSeriesZ Test 18: \t\t\t\t PASSED \n");
		 // printf("It is not a unit\n");
	 } else {
		fprintf(stderr, "PowerSeriesZ Test 18: \t\t\t\t FAILED \n");
		 // printf("It is a unit\n");
		exit(1);
	 }

	 freePolynomial_AAZ(a);
	 destroyPowerSeries_PSZ(a_ps);
     return 0;
}


/**
 * construct constant power series
 */
int test19()
{
	PolyZ_ptr a = generate_altarrZ_var_defined("2",sym, 2);
	PowerSeriesZ_t* constant = constPowerSeries_PSZ(a->elems->coef, 2);
	PolyZ_ptr homog = homogPart_PSZ(0, constant);
	if (isExactlyEqual_AAZ(a, homog)) {
		fprintf(stderr, "PowerSeriesZ Test 19: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 19: \t\t\t\t FAILED \n");
	    exit(1);
	}

	freePolynomial_AAZ(a);
	destroyPowerSeries_PSZ(constant);


	return 0;
}


/**
 *
 */
int test20()
{

	PowerSeriesZ_t* ps = zeroPowerSeries_PSZ();
	PowerSeriesZ_t* sub = subPowerSeries_PSZ(ps,ps);
	PolyZ_ptr homog = homogPart_PSZ(0, sub);
	if (isZero_AAZ(homog)) {
		fprintf(stderr, "PowerSeriesZ Test 20: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 20: \t\t\t\t FAILED \n");
		exit(1);
	}

	destroyPowerSeries_PSZ(ps);
	destroyPowerSeries_PSZ(sub);
    return 0;
}



/**
 *
 */
int test21()
{
	PowerSeriesZ_t* psOne = onePowerSeries_PSZ(1);
	PowerSeriesZ_t* psZero = zeroPowerSeries_PSZ();
	PowerSeriesZ_t* multiply = multiplyPowerSeries_PSZ(psOne,psZero);
	PolyZ_ptr homog = homogPart_PSZ(0, multiply);
	if (isZero_AAZ(homog)) {
		fprintf(stderr, "PowerSeriesZ Test 21: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 21: \t\t\t\t FAILED \n");
		exit(1);
	}

	destroyPowerSeries_PSZ(psOne);
	destroyPowerSeries_PSZ(psZero);
	destroyPowerSeries_PSZ(multiply);

    return 0;
}


/**
*
*/
int test22()
{
    PolyZ_ptr a = generate_altarrZ_var_defined("-2*x+x^2+1",vars, 3);
    PolyZ_ptr b = generate_altarrZ_var_defined("0",vars, 3);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = dividePowerSeries_PSZ(b_ps,a_ps);
    PolyZ_ptr homog = homogPart_PSZ(0,c);
    if (isZero_AAZ(homog)) {
		fprintf(stderr, "PowerSeriesZ Test 22: \t\t\t\t PASSED \n");
	}
    else {
		fprintf(stderr, "PowerSeriesZ Test 22: \t\t\t\t FAILED \n");
    	exit(1);
    }

    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    freePolynomial_AAZ(a);
    freePolynomial_AAZ(b);
    return 0;
}


/**
* Subtract zero power series from zero power series
*/
int test23()
{
    PolyZ_ptr a = generate_altarrZ_var_defined("0",vars, 3);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t* c = subPowerSeries_PSZ(b_ps,a_ps);
    PolyZ_ptr  homog = homogPart_PSZ(0,c);
    if (isZero_AAZ(homog)) {
		fprintf(stderr, "PowerSeriesZ Test 23: \t\t\t\t PASSED \n");
    }
    else {
		fprintf(stderr, "PowerSeriesZ Test 23: \t\t\t\t FAILED \n");
		exit(1);
    }
    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    freePolynomial_AAZ(a);
    return 0;
}


/**
* Subtract one power series from zero power series
*/
int test24()
{
    PolyZ_ptr a = generate_altarrZ_var_defined("1",vars, 3);
    PolyZ_ptr b = generate_altarrZ_var_defined("0",vars, 3);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = subPowerSeries_PSZ(b_ps,a_ps);
    PolyZ_ptr  homog = homogPart_PSZ(0,c);
    if (isNegativeOne_AAZ(homog)) {
		fprintf(stderr, "PowerSeriesZ Test 24: \t\t\t\t PASSED \n");
    }
    else {
		fprintf(stderr, "PowerSeriesZ Test 24: \t\t\t\t FAILED \n");
        exit(1);
    }

    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    freePolynomial_AAZ(a);
    freePolynomial_AAZ(b);
    return 0;
}


/**
*
*/
int test25()
{
    PolyZ_ptr a = generate_altarrZ_var_defined("1",sym, 2);
    PolyZ_ptr b = generate_altarrZ_var_defined("1 + x + y + x^2 + y^2",sym, 2);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = addPowerSeries_PSZ(b_ps,a_ps);
    PolyZ_ptr  homog = homogPart_PSZ(2,c);
    PolyZ_ptr  testpoly = generate_altarrZ_var_defined("x^2 + y^2", sym, 2);
    testpoly = subPolynomials_AAZ_inp(testpoly, homog,2);
    if (isZero_AAZ(testpoly)) {
		fprintf(stderr, "PowerSeriesZ Test 25: \t\t\t\t PASSED \n");
    }

    else {
		fprintf(stderr, "PowerSeriesZ Test 25: \t\t\t\t FAILED \n");
        exit(1);
    }

    freePolynomial_AAZ(a);
    freePolynomial_AAZ(b);
    freePolynomial_AAZ(testpoly);
    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    return 0;
}


/**
* Construct a constant power series
*/
int test26()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("2",sym, 2);
	PowerSeriesZ_t* constant = constPowerSeries_PSZ(a->elems->coef, 2);
	PolyZ_ptr b = generate_altarrZ_var_defined("1 + x + y + x^2 + y^2",sym, 2);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	PowerSeriesZ_t* sum = addPowerSeries_PSZ(b_ps,constant);
	PowerSeriesZ_t* sub = subPowerSeries_PSZ(b_ps,constant);
	PolyZ_ptr  testpoly_1 = generate_altarrZ_var_defined("3", sym, 2);
	PolyZ_ptr  homog_1 = homogPart_PSZ(0,sum);
	testpoly_1 = subPolynomials_AAZ_inp(testpoly_1, homog_1,2);
	PolyZ_ptr  testpoly_2 = generate_altarrZ_var_defined("-1", sym, 2);
	PolyZ_ptr  homog_2 = homogPart_PSZ(0,sub);
	testpoly_2 = subPolynomials_AAZ_inp(testpoly_2,homog_2,2);
	if (isZero_AAZ(testpoly_1) && isZero_AAZ(testpoly_2)) {
		fprintf(stderr, "PowerSeriesZ Test 26: \t\t\t\t PASSED \n");
	} else {
		fprintf(stderr, "PowerSeriesZ Test 26: \t\t\t\t FAILED \n");
		exit(1);
	}


	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly_1);
	freePolynomial_AAZ(testpoly_2);

	destroyPowerSeries_PSZ(constant);
	destroyPowerSeries_PSZ(sum);
	destroyPowerSeries_PSZ(sub);
	return 0;

}


int test27()
{

	PowerSeriesZ_t* ps = zeroPowerSeries_PSZ();
	PowerSeriesZ_t* sub = subPowerSeries_PSZ(ps,ps);
	PowerSeriesZ_t* sum = addPowerSeries_PSZ(ps,ps);
	PolyZ_ptr  homogSum = homogPart_PSZ(0,sum);
	PolyZ_ptr  homogSub = homogPart_PSZ(0,sub);
	PowerSeriesZ_t* multiply = multiplyPowerSeries_PSZ(ps,ps);
	PolyZ_ptr  homogMultiply = homogPart_PSZ(0,multiply);
	if ((isZero_AAZ(homogSum)) && (isZero_AAZ(homogSub)) && ((isZero_AAZ(homogMultiply)))) {
		fprintf(stderr, "PowerSeriesZ Test 27: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 27: \t\t\t\t FAILED \n");
		exit(1);
	}


	destroyPowerSeries_PSZ(multiply);
	destroyPowerSeries_PSZ(sum);
	destroyPowerSeries_PSZ(sub);
	destroyPowerSeries_PSZ(ps);
	return 0;
}



int test28()
{

	PowerSeriesZ_t* ps = zeroPowerSeries_PSZ();
	PowerSeriesZ_t* sub = subPowerSeries_PSZ(ps,ps);
	PowerSeriesZ_t* sum = addPowerSeries_PSZ(ps,ps);
	PolyZ_ptr  homogSum = homogPart_PSZ(10,sum);
	PolyZ_ptr  homogSub = homogPart_PSZ(20,sub);
	PowerSeriesZ_t* multiply = multiplyPowerSeries_PSZ(ps,ps);
	PolyZ_ptr  homogMultiply = homogPart_PSZ(10,multiply);
	if ((isZero_AAZ(homogSum)) && (isZero_AAZ(homogSub)) && ((isZero_AAZ(homogMultiply)))) {
		fprintf(stderr, "PowerSeriesZ Test 28: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 28: \t\t\t\t FAILED \n");
		exit(1);
	}


	destroyPowerSeries_PSZ(multiply);
	destroyPowerSeries_PSZ(sum);
	destroyPowerSeries_PSZ(sub);
	destroyPowerSeries_PSZ(ps);
    return 0;
}


/**
*
*/
int test29()
{
    PolyZ_ptr a = generate_altarrZ_var_defined("1",sym, 2);
    PolyZ_ptr b = generate_altarrZ_var_defined("1 + x + y + x^2 + y^2",sym, 2);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = subPowerSeries_PSZ(a_ps,b_ps);
    PolyZ_ptr  homog = homogPart_PSZ(2,c);
    PolyZ_ptr  testpoly = generate_altarrZ_var_defined("-x^2 - y^2", sym, 2);
    testpoly = subPolynomials_AAZ_inp(testpoly,homog,2);
    if (isZero_AAZ(testpoly)) {
		fprintf(stderr, "PowerSeriesZ Test 29: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 29: \t\t\t\t FAILED \n");
		exit(1);
	}

    freePolynomial_AAZ(a);
    freePolynomial_AAZ(b);
    freePolynomial_AAZ(testpoly);

    destroyPowerSeries_PSZ (c);
    destroyPowerSeries_PSZ (a_ps);
    destroyPowerSeries_PSZ (b_ps);
    return 0;
}


int test30()
{
	PolyZ_ptr  a = NULL;
	PolyZ_ptr  b = generate_altarrZ_var_defined("z^5-y^5+x*t^2*z+x^4*t+1", syms, 4);
    PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t* c = subPowerSeries_PSZ(a_ps,b_ps);
    PolyZ_ptr homog = homogPart_PSZ(5, c);
    PolyZ_ptr testpoly = generate_altarrZ_var_defined("-z^5+y^5-x^4*t", syms, 4);
    testpoly = subPolynomials_AAZ_inp(testpoly,homog,4);

    if (isZero_AAZ(testpoly)) {
		fprintf(stderr, "PowerSeriesZ Test 30: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeriesZ Test 30: \t\t\t\t FAILED \n");
		exit(1);
	}


	freePolynomial_AAZ(b);
	freePolynomial_AAZ(testpoly);

	destroyPowerSeries_PSZ(a_ps);
	destroyPowerSeries_PSZ(b_ps);
	destroyPowerSeries_PSZ(c);

	return 0;
}

int test31()
{

	PolyZ_ptr  b = generate_altarrZ_var_defined("z^5-y^5+x*t^2*z+x^4*t+1", syms, 4);

    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t*  c_ps = negatePowerSeries_PSZ(b_ps);

    homogPart_PSZ(5, c_ps);

    PolyZ_ptr poly_1 = generate_altarrZ_var_defined("z^5-y^5+x^4*t", syms, 4);
    poly_1 = addPolynomials_AAZ_inp(poly_1, c_ps->polys[5],4);
    PolyZ_ptr poly_2 = generate_altarrZ_var_defined("x*t^2*z", syms, 4);
    poly_2 = addPolynomials_AAZ_inp(poly_2,c_ps->polys[4],4);
    PolyZ_ptr poly_3 = generate_altarrZ_var_defined("1", syms, 4);
    poly_3 = addPolynomials_AAZ_inp(poly_3,c_ps->polys[0],4);

    if (isZero_AAZ(poly_1) && isZero_AAZ(poly_2) && isZero_AAZ(poly_3)) {
		fprintf(stderr, "PowerSeriesZ Test 31: \t\t\t\t PASSED \n");
    } else {
		fprintf(stderr, "PowerSeriesZ Test 31: \t\t\t\t FAILED \n");
		exit(1);
    }



    destroyPowerSeries_PSZ(c_ps);
    destroyPowerSeries_PSZ(b_ps);
    freePolynomial_AAZ(b);
    freePolynomial_AAZ(poly_1);
    freePolynomial_AAZ(poly_2);
    freePolynomial_AAZ(poly_3);

    return 0;

}

int test32()
{

	PolyZ_ptr  b = NULL;
    PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
    PowerSeriesZ_t*  c_ps = negatePowerSeries_PSZ(b_ps);
    b = addPolynomials_AAZ_inp(b, c_ps->polys[0],4);

    if (isZero_AAZ(b)) {
		fprintf(stderr, "PowerSeriesZ Test 32: \t\t\t\t PASSED \n");
    } else {
		fprintf(stderr, "PowerSeriesZ Test 32: \t\t\t\t FAILED \n");
    	exit(1);
    }

    destroyPowerSeries_PSZ(c_ps);
    destroyPowerSeries_PSZ(b_ps);
    freePolynomial_AAZ(b);
    return 0;

}


int test33()
{
	// PolyZ_ptr a = generate_altarrZ_var_defined("1",sym, 2);
	PolyZ_ptr b = generate_altarrZ_var_defined("1-x",sym, 2);
	// PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	TestingPowerSeriesInversion(b_ps, 30, sym, 2);

	// freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	// destroyPowerSeries_PSZ(a_ps);
	destroyPowerSeries_PSZ(b_ps);

	fprintf(stderr, "PowerSeriesZ Test 33: \t\t\t\t PASSED \n");
	return 0;
}


int test34()
{
	PolyZ_ptr a = generate_altarrZ_var_defined("y",sym, 2);
	PolyZ_ptr b = generate_altarrZ_var_defined("1-x",sym, 2);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	TestingPowerSeriesDivision(a_ps, b_ps,10, sym, 2);

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	destroyPowerSeries_PSZ(a_ps);
	destroyPowerSeries_PSZ(b_ps);

	fprintf(stderr, "PowerSeriesZ Test 34: \t\t\t\t PASSED \n");

	return 0;
}


int test35()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("-2*x*t+x^2+y*z+t^4+1",syms, 4);
	PolyZ_ptr b = generate_altarrZ_var_defined("1-t-x",syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	TestingPowerSeriesDivision(b_ps, a_ps,20, syms, 4);

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	destroyPowerSeries_PSZ(a_ps);
	destroyPowerSeries_PSZ(b_ps);

	fprintf(stderr, "PowerSeriesZ Test 35: \t\t\t\t PASSED \n");

	return 0;
}


int test36()
{

	PolyZ_ptr a = generate_altarrZ_var_defined("-2*x*t+x^2+y*z+t^4+1",syms, 4);
	PolyZ_ptr b = generate_altarrZ_var_defined("1-t*x^3-x^6+2*y",syms, 4);
	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
	TestingPowerSeriesDivision(b_ps, a_ps,20, syms, 4);

	freePolynomial_AAZ(a);
	freePolynomial_AAZ(b);
	destroyPowerSeries_PSZ(a_ps);
	destroyPowerSeries_PSZ(b_ps);

	fprintf(stderr, "PowerSeriesZ Test 36: \t\t\t\t PASSED \n");
	return 0;
}


// int test37()
// {
// 	PolyZ_ptr a = generate_altarrZ_var_defined("1+x+y+z",vars, 3);
// 	PolyZ_ptr b = generate_altarrZ_var_defined("1",vars, 3);
// 	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
// 	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
// 	PowerSeriesZ_t* c = dividePowerSeries_PSZ(b_ps,a_ps);
// 	unsigned long long start;
// 	_startTimer(&start);
// 	PolyZ_ptr homog = homogPart_PSZ(600,c);
// 	float elapsedTime;
// 	_stopTimer(&start, &elapsedTime);
// 	fprintf(stderr, "Time sys37: %f\n", elapsedTime);
//     return 0;

// }


// int test38()
// {
// 	PolyZ_ptr a = generate_altarrZ_var_defined("1+x+y",sym, 2);
// 	PolyZ_ptr b = generate_altarrZ_var_defined("1",sym, 2);
// 	PowerSeriesZ_t*  a_ps = convertPolyToPowerSeries_PSZ(a);
// 	PowerSeriesZ_t*  b_ps = convertPolyToPowerSeries_PSZ(b);
// 	PowerSeriesZ_t* c = dividePowerSeries_PSZ(b_ps,a_ps);
// 	for (int indx = 100; indx <= 2000; indx = indx +100)
// 	TestingPowerSeriesDivision(b_ps, a_ps, indx, sym, 2);
//     return 0;

// }
