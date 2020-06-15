#include "../../include/PowerSeries/PowerSeries.h"
#include "../../include/Utils/Unix_Timer.h"
#include "../../include/Parser/bpas_parser.h"
#include <time.h>


const char* sym[] = {"x", "y"};
const char* syms[] = {"x", "y", "z", "t"};
const char* vars[] = {"x", "y", "z"};



void TestingPowerSeriesDivision(PowerSeries_t* f, PowerSeries_t* g, int bound, const char** sym, int n) {

	PowerSeries_t* h = dividePowerSeries_PS(f,g);
	PowerSeries_t* prod = multiplyPowerSeries_PS(g,h);
	PowerSeries_t* sub = subPowerSeries_PS(f,prod);

	for (int d = 1; d <= bound; ++d) {
		//fprintf(stderr, "In degree: %d\n ", d);
		//fprintf(stderr, "\n");
		//fprintf(stderr, "polyPart division of degree d = %d\n", d);
	    //fprintf(stderr, "\n");
		//printPolyClean_AA(stderr,polynomialPart_PS(d, h), sym, n);
		//fprintf(stderr, "\n");
		Poly_ptr polyPart = polynomialPart_PS(d,  sub);

		if (!isZero_AA(polyPart)) {
			 fprintf(stderr, "ERROR IN TEST POWER SERIES DIVISION\n");
			 exit(1);
		}

	}

	destroyPowerSeries_PS(h);
	destroyPowerSeries_PS(prod);
	destroyPowerSeries_PS(sub);
}

void TestingPowerSeriesInversion(PowerSeries_t* f, int bound, const char** sym, int n) {

	PowerSeries_t* h = inversePowerSeries_PS(f);
	PowerSeries_t* prod = multiplyPowerSeries_PS(f, h);
	for (int d = 1; d <= bound; ++d) {
		//fprintf(stderr, "In degree: %d\n ", d);
		//fprintf(stderr, "\n");
		//fprintf(stderr, "polyPart division of degree d = %d\n", d);
	    //fprintf(stderr, "\n");
		//printPolyClean_AA(stderr,polynomialPart_PS(d, h), sym, n);
		//fprintf(stderr, "\n");
		Poly_ptr polyPart = polynomialPart_PS(d,  prod);

		if (!isOne_AA(polyPart)) {
			 fprintf(stderr, "ERROR IN TEST POWER SERIES DIVISION\n");
			 exit(1);
		}

	}

	destroyPowerSeries_PS(h);
	destroyPowerSeries_PS(prod);
}




int test1()
{
	Poly_ptr  polynom = generate_altarr_var_defined("x^10+1", syms, 4);
	PowerSeries_t*  a = convertPolyToPowerSeries_PS(polynom);
	PowerSeries_t* b = multiplyPowerSeries_PS(a,a);
	Poly_ptr homog = homogPart_PS(20, b);
	Poly_ptr testpoly = generate_altarr_var_defined("x^20", syms, 4);
	testpoly = subPolynomials_AA(testpoly,homog,4);

	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 1: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 1: \t\t\t\t FAILED \n");
        exit(1);
	}
	freePolynomial_AA(polynom);
	freePolynomial_AA(testpoly);

	destroyPowerSeries_PS (a);
	destroyPowerSeries_PS (b);
    return 0;
}


int test2()
{
	Poly_ptr  a = generate_altarr_var_defined("x+100", syms, 4);
	Poly_ptr  b = generate_altarr_var_defined("z^5+x^4*t+1", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);
	//printf("%d\n",c->genOrder);
	//printPoly_AA (stderr, homogPart_PS(0, c), syms, 4);
	Poly_ptr homog = homogPart_PS(0, c);
	if (isConstant_AA(homog))  {
	    fprintf(stderr, "PowerSeries Test 2: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 2: \t\t\t\t FAILED \n");
	    exit(1);
	}

	freePolynomial_AA(a);
	freePolynomial_AA(b);

	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	destroyPowerSeries_PS (c);
    return 0;
}

int test3()
{
	Poly_ptr  a = generate_altarr_var_defined("x*z*t^4+y+x^2+1", syms, 4);
	Poly_ptr testpoly = generate_altarr_var_defined("x*z*t^4", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	Poly_ptr subpoly = subPolynomials_AA(a_ps->polys[6],testpoly,4);
	//for (int i = 0; i <=a_ps->deg; i++) {
	//		fprintf (stderr, "a_ps->polys[%d]=", i);
	//		printPoly_AA (stderr, a_ps->polys[i], syms, 4);
	//		fprintf(stderr, "\n");
	//	}
	if (isZero_AA(subpoly)) {
	    fprintf(stderr, "PowerSeries Test 3: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 3: \t\t\t\t FAILED \n");
	    exit(1);
	}

	freePolynomial_AA(a);
	freePolynomial_AA(subpoly);
	freePolynomial_AA(testpoly);
	destroyPowerSeries_PS (a_ps);
    return 0;
}


int test4()
{
	Poly_ptr  a = generate_altarr_var_defined("1", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	Poly_ptr  b = generate_altarr_var_defined("0", syms, 4);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);
	Poly_ptr homog = homogPart_PS(3, c);
	if (isZero_AA(homog)) {
	    fprintf(stderr, "PowerSeries Test 4: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 4: \t\t\t\t FAILED \n");
	    exit(1);
	}
	//for (int i = 0; i <=a_ps->deg; i++) {
	//		fprintf (stderr, "a_ps->polys[%d]=", i);
	//		printPoly_AA (stderr, a_ps->polys[i], syms, 4);
	//		fprintf(stderr, "\n");
	//	}

	freePolynomial_AA(a);
	freePolynomial_AA(b);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	destroyPowerSeries_PS (c);
    return 0;
}

int test5()
{
	Poly_ptr  a = generate_altarr_var_defined("1", syms, 4);
	Poly_ptr  b = generate_altarr_var_defined("2*x^2 + 2*y^2 + 2*x + 2*y +2", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);
	Poly_ptr homog = homogPart_PS(2, c);
	Poly_ptr testpoly = generate_altarr_var_defined("2*x^2 + 2*y^2", syms, 4);
	testpoly = subPolynomials_AA_inp(testpoly,homog,4);
	//for (int i = 0; i <=a_ps->deg; i++) {
	//		fprintf (stderr, "a_ps->polys[%d]=", i);
	//		printPoly_AA (stderr, a_ps->polys[i], syms, 4);
	//		fprintf(stderr, "\n");
	//	}
	if (isZero_AA(testpoly))  {
	    fprintf(stderr, "PowerSeries Test 5: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 5: \t\t\t\t FAILED \n");
	    exit(1);
	}
	freePolynomial_AA(a);
	freePolynomial_AA(b);
	freePolynomial_AA(testpoly);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	destroyPowerSeries_PS (c);
    return 0;
}



int test6()
{
	Poly_ptr  a = generate_altarr_var_defined("z+1", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,a_ps);
	Poly_ptr homog = homogPart_PS(2, c);
	Poly_ptr testpoly = generate_altarr_var_defined("z^2", syms, 4);
	testpoly = subPolynomials_AA_inp(testpoly,homog,4);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 6: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 6: \t\t\t\t FAILED \n");
		exit(1);
	}
	freePolynomial_AA(a);
	freePolynomial_AA(testpoly);
	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (a_ps);
	a = NULL;

    return 0;
}

int test7()
{
	Poly_ptr  a = generate_altarr_var_defined("x+1", syms, 4);
    Poly_ptr  b = generate_altarr_var_defined("x^2+1", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);
	Poly_ptr homog = homogPart_PS(3, c);
	Poly_ptr testpoly = generate_altarr_var_defined("x^3", syms, 4);
	testpoly = subPolynomials_AA_inp(testpoly,homog,4);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 7: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 7: \t\t\t\t FAILED \n");
		exit(1);
	}

	freePolynomial_AA(a);
	freePolynomial_AA(b);
	freePolynomial_AA(testpoly);
	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
    return 0;

}

/**
 * Note that multiplication process truncates the result up to the min degree of ancestors
 * in the following the truncation degree is 3. Thus the array of polys does not contain terms
 * of degree at least 4. But when we call homogPart_PS and ask for terms of degree 5 and look at
 * polys again we can see that this time the array polys has been updated and contains terms of degree 5
 */
int test8()
{
	Poly_ptr  a = generate_altarr_var_defined("x^3 + x + z", vars, 3);
	Poly_ptr  b = generate_altarr_var_defined("z^2*y*x + z*y + 1", vars, 3);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);
	Poly_ptr homog = homogPart_PS(5, c);
	Poly_ptr testpoly = generate_altarr_var_defined("x^3*z*y+x^2*y*z^2+z^3*y*x", vars, 3);
	testpoly = subPolynomials_AA_inp(testpoly,homog,3);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 8: \t\t\t\t PASSED \n");
	}
	else  {
	    fprintf(stderr, "PowerSeries Test 8: \t\t\t\t FAILED \n");
		exit(1);
	}

	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (b_ps);
	destroyPowerSeries_PS (a_ps);
	freePolynomial_AA(a);
	freePolynomial_AA(b);
	freePolynomial_AA(testpoly);
	c = NULL;
	return 0;
}


int test9()
{
	Poly_ptr  a = generate_altarr_var_defined("x^2*t+z^7*y*x-2*x^4+1", syms, 4);
	Poly_ptr  b = generate_altarr_var_defined("x^20+x*t*y*z+z^4+2", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);

	Poly_ptr d = homogPart_PS(8, c);
	if (isZero_AA(d)|| isZero_AA(c->polys[8])) {
	    fprintf(stderr, "PowerSeries Test 9: \t\t\t\t FAILED \n");
		exit(1);
	} else {
	    fprintf(stderr, "PowerSeries Test 9: \t\t\t\t PASSED \n");
	}

	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	freePolynomial_AA(a);
	freePolynomial_AA(b);

	return 0;
}


int test10()
{
	Poly_ptr  a = generate_altarr_var_defined("x^2*t^3+z^7*y*x-2*x^4+y^12+x", syms, 4);
	Poly_ptr  b = generate_altarr_var_defined("x^2+x*t*y*z+z^4*y*x*t+2", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = multiplyPowerSeries_PS(a_ps,b_ps);
	Poly_ptr testpoly = generate_altarr_var_defined("-2*x^6", syms, 4);
	Poly_ptr homog = homogPart_PS(6, c);
	testpoly = subPolynomials_AA(testpoly, homog,4);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 10: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 10: \t\t\t\t FAILED \n");
	    exit(1);
	}

	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	freePolynomial_AA(a);
	freePolynomial_AA(testpoly);
	freePolynomial_AA(b);
    return 0;
}


int test11()
{
	Poly_ptr a = generate_altarr_var_defined("-2*x+x^2+1",syms, 4);
	Poly_ptr b = generate_altarr_var_defined("1",syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = dividePowerSeries_PS(b_ps,a_ps);
	Poly_ptr testpoly = generate_altarr_var_defined("5*x^4", syms, 4);
	Poly_ptr homog = homogPart_PS(4,c);
	testpoly = subPolynomials_AA_inp(testpoly,homog,4);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 11: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 11: \t\t\t\t FAILED \n");
	    exit(1);
	}

	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	freePolynomial_AA(a);
	freePolynomial_AA(b);
	freePolynomial_AA(testpoly);
    return 0;

}


int test12()
{
	Poly_ptr a = generate_altarr_var_defined("-2*x*t+x^2+y*z+t^4+1",syms, 4);
	Poly_ptr b = generate_altarr_var_defined("1-t-x",syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* c = dividePowerSeries_PS(b_ps,a_ps);
	Poly_ptr testpoly = generate_altarr_var_defined("t^4-4*t^2*x^2+4*t*x^3+4*t*x*y*z-x^4-2*x^2*y*z-y^2*z^2", syms, 4);
	Poly_ptr homog = homogPart_PS(4,c);
	testpoly = addPolynomials_AA_inp(testpoly, homog,4);
	if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 12: \t\t\t\t PASSED \n");
	}
	else {
	    fprintf(stderr, "PowerSeries Test 12: \t\t\t\t FAILED \n");
	    exit(1);
	}

    // for (int i = 0; i <= 8; ++i) {
    //     printPoly_AA (stderr, homogPart_PS(i,c), syms, 4);
    //     fprintf(stderr, "\n");
    // }

	destroyPowerSeries_PS (c);
	destroyPowerSeries_PS (a_ps);
	destroyPowerSeries_PS (b_ps);
	freePolynomial_AA(a);
	freePolynomial_AA(b);
	freePolynomial_AA(testpoly);
    return 0;

}


int test13()
{
    Poly_ptr  a = generate_altarr_var_defined("x^3 + x + z", vars, 3);
    Poly_ptr  b = generate_altarr_var_defined("z^2*y*x + z*y + 1", vars, 3);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = addPowerSeries_PS(a_ps,b_ps);
    Poly_ptr homog = homogPart_PS(4, c);
    Poly_ptr testpoly = generate_altarr_var_defined("z^2*y*x", vars, 3);
    testpoly = subPolynomials_AA_inp(testpoly,homog,4);
    if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 13: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeries Test 13: \t\t\t\t FAILED \n");
        exit(1);
    }

    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    freePolynomial_AA(a);
    freePolynomial_AA(b);
    freePolynomial_AA(testpoly);
    return 0;
 }



int test14()
{
    Poly_ptr  a = generate_altarr_var_defined("1+x+y+z+x^2+y^2+z^2", vars, 3);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t* c = addPowerSeries_PS(a_ps,a_ps);
    Poly_ptr poly = c->polys[0];
    if (isConstant_AA(poly)) {
	    fprintf(stderr, "PowerSeries Test 14: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeries Test 14: \t\t\t\t FAILED \n");
        exit(1);
    }

    freePolynomial_AA(a);
    destroyPowerSeries_PS(a_ps);
    destroyPowerSeries_PS(c);
    return 0;
}


int test15()
{
	Poly_ptr  a = generate_altarr_var_defined("1+x-y+z+2*x^2-4*y^2+z^6", vars, 3);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t* c = subPowerSeries_PS(a_ps,a_ps);
    Poly_ptr test = polynomialPart_PS(6, c);
    if (isZero_AA(test)) {
	    fprintf(stderr, "PowerSeries Test 15: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeries Test 15: \t\t\t\t FAILED \n");
        exit(1);
    }

    freePolynomial_AA(a);
    freePolynomial_AA(test);
    destroyPowerSeries_PS(a_ps);
    destroyPowerSeries_PS(c);
    return 0;
 }


int test16()
{
    Poly_ptr  a = generate_altarr_var_defined("x+100", syms, 4);
    Poly_ptr  b = generate_altarr_var_defined("z^5-y^5+x*t^2*z+x^4*t+1", syms, 4);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = subPowerSeries_PS(a_ps,b_ps);
    Poly_ptr homog = homogPart_PS(5, c);
    Poly_ptr testpoly = generate_altarr_var_defined("-x^4*t+y^5-z^5", syms, 4);
    testpoly = subPolynomials_AA_inp(testpoly,homog,4);
    if (isZero_AA(testpoly)) {
	    fprintf(stderr, "PowerSeries Test 16: \t\t\t\t PASSED \n");
    }
    else {
	    fprintf(stderr, "PowerSeries Test 16: \t\t\t\t FAILED \n");
    	exit(1);
    }

    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    freePolynomial_AA(a);
    freePolynomial_AA(b);
    freePolynomial_AA(testpoly);
    return 0;


}


/**
 *
 */
int test17()
{

	Poly_ptr  a = generate_altarr_var_defined("x*y^6+z^3-x+100", syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	int x = isUnit_PS(a_ps);
	if (x) {
		fprintf(stderr, "PowerSeries Test 17: \t\t\t\t PASSED \n");
		 // printf("It is a unit\n");
	} else {
		fprintf(stderr, "PowerSeries Test 17: \t\t\t\t FAILED \n");
		exit(1);
		 // printf("It is not a unit\n");
	}
	// printf("%d\n", x);

	freePolynomial_AA(a);
	destroyPowerSeries_PS(a_ps);
	return 0;
}


/**
 *
 */
int test18()
{

	 Poly_ptr  a = generate_altarr_var_defined("x-y+z^3*t+3*x^2-2*x*y", syms, 4);
	 PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	 int x = isUnit_PS(a_ps);
	 if (!x) {
		fprintf(stderr, "PowerSeries Test 18: \t\t\t\t PASSED \n");
		 // printf("It is not a unit\n");
	 } else {
		fprintf(stderr, "PowerSeries Test 18: \t\t\t\t FAILED \n");
		 // printf("It is a unit\n");
		exit(1);
	 }

	 freePolynomial_AA(a);
	 destroyPowerSeries_PS(a_ps);
     return 0;
}


/**
 * construct constant power series
 */
int test19()
{
	Poly_ptr a = generate_altarr_var_defined("2",sym, 2);
	PowerSeries_t* constant = constPowerSeries_PS(a->elems->coef, 2);
	Poly_ptr homog = homogPart_PS(0, constant);
	if (isExactlyEqual_AA(a, homog)) {
		fprintf(stderr, "PowerSeries Test 19: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 19: \t\t\t\t FAILED \n");
	    exit(1);
	}

	freePolynomial_AA(a);
	destroyPowerSeries_PS(constant);


	return 0;
}


/**
 *
 */
int test20()
{

	PowerSeries_t* ps = zeroPowerSeries_PS();
	PowerSeries_t* sub = subPowerSeries_PS(ps,ps);
	Poly_ptr homog = homogPart_PS(0, sub);
	if (isZero_AA(homog)) {
		fprintf(stderr, "PowerSeries Test 20: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 20: \t\t\t\t FAILED \n");
		exit(1);
	}

	destroyPowerSeries_PS(ps);
	destroyPowerSeries_PS(sub);
    return 0;
}



/**
 *
 */
int test21()
{
	PowerSeries_t* psOne = onePowerSeries_PS(1);
	PowerSeries_t* psZero = zeroPowerSeries_PS();
	PowerSeries_t* multiply = multiplyPowerSeries_PS(psOne,psZero);
	Poly_ptr homog = homogPart_PS(0, multiply);
	if (isZero_AA(homog)) {
		fprintf(stderr, "PowerSeries Test 21: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 21: \t\t\t\t FAILED \n");
		exit(1);
	}

	destroyPowerSeries_PS(psOne);
	destroyPowerSeries_PS(psZero);
	destroyPowerSeries_PS(multiply);

    return 0;
}


/**
*
*/
int test22()
{
    Poly_ptr a = generate_altarr_var_defined("-2*x+x^2+1",vars, 3);
    Poly_ptr b = generate_altarr_var_defined("0",vars, 3);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = dividePowerSeries_PS(b_ps,a_ps);
    Poly_ptr homog = homogPart_PS(0,c);
    if (isZero_AA(homog)) {
		fprintf(stderr, "PowerSeries Test 22: \t\t\t\t PASSED \n");
	}
    else {
		fprintf(stderr, "PowerSeries Test 22: \t\t\t\t FAILED \n");
    	exit(1);
    }

    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    freePolynomial_AA(a);
    freePolynomial_AA(b);
    return 0;
}


/**
* Subtract zero power series from zero power series
*/
int test23()
{
    Poly_ptr a = generate_altarr_var_defined("0",vars, 3);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t* c = subPowerSeries_PS(b_ps,a_ps);
    Poly_ptr  homog = homogPart_PS(0,c);
    if (isZero_AA(homog)) {
		fprintf(stderr, "PowerSeries Test 23: \t\t\t\t PASSED \n");
    }
    else {
		fprintf(stderr, "PowerSeries Test 23: \t\t\t\t FAILED \n");
		exit(1);
    }
    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    freePolynomial_AA(a);
    return 0;
}


/**
* Subtract one power series from zero power series
*/
int test24()
{
    Poly_ptr a = generate_altarr_var_defined("1",vars, 3);
    Poly_ptr b = generate_altarr_var_defined("0",vars, 3);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = subPowerSeries_PS(b_ps,a_ps);
    Poly_ptr  homog = homogPart_PS(0,c);
    if (isNegativeOne_AA(homog)) {
		fprintf(stderr, "PowerSeries Test 24: \t\t\t\t PASSED \n");
    }
    else {
		fprintf(stderr, "PowerSeries Test 24: \t\t\t\t FAILED \n");
        exit(1);
    }

    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    freePolynomial_AA(a);
    freePolynomial_AA(b);
    return 0;
}


/**
*
*/
int test25()
{
    Poly_ptr a = generate_altarr_var_defined("1",sym, 2);
    Poly_ptr b = generate_altarr_var_defined("1 + x + y + x^2 + y^2",sym, 2);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = addPowerSeries_PS(b_ps,a_ps);
    Poly_ptr  homog = homogPart_PS(2,c);
    Poly_ptr  testpoly = generate_altarr_var_defined("x^2 + y^2", sym, 2);
    testpoly = subPolynomials_AA_inp(testpoly, homog,2);
    if (isZero_AA(testpoly)) {
		fprintf(stderr, "PowerSeries Test 25: \t\t\t\t PASSED \n");
    }

    else {
		fprintf(stderr, "PowerSeries Test 25: \t\t\t\t FAILED \n");
        exit(1);
    }

    freePolynomial_AA(a);
    freePolynomial_AA(b);
    freePolynomial_AA(testpoly);
    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    return 0;
}


/**
* Construct a constant power series
*/
int test26()
{

	Poly_ptr a = generate_altarr_var_defined("2",sym, 2);
	PowerSeries_t* constant = constPowerSeries_PS(a->elems->coef, 2);
	Poly_ptr b = generate_altarr_var_defined("1 + x + y + x^2 + y^2",sym, 2);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	PowerSeries_t* sum = addPowerSeries_PS(b_ps,constant);
	PowerSeries_t* sub = subPowerSeries_PS(b_ps,constant);
	Poly_ptr  testpoly_1 = generate_altarr_var_defined("3", sym, 2);
	Poly_ptr  homog_1 = homogPart_PS(0,sum);
	testpoly_1 = subPolynomials_AA_inp(testpoly_1, homog_1,2);
	Poly_ptr  testpoly_2 = generate_altarr_var_defined("-1", sym, 2);
	Poly_ptr  homog_2 = homogPart_PS(0,sub);
	testpoly_2 = subPolynomials_AA_inp(testpoly_2,homog_2,2);
	if (isZero_AA(testpoly_1) && isZero_AA(testpoly_2)) {
		fprintf(stderr, "PowerSeries Test 26: \t\t\t\t PASSED \n");
	} else {
		fprintf(stderr, "PowerSeries Test 26: \t\t\t\t FAILED \n");
		exit(1);
	}


	freePolynomial_AA(a);
	freePolynomial_AA(b);
	freePolynomial_AA(testpoly_1);
	freePolynomial_AA(testpoly_2);

	destroyPowerSeries_PS(constant);
	destroyPowerSeries_PS(sum);
	destroyPowerSeries_PS(sub);
	return 0;

}


int test27()
{

	PowerSeries_t* ps = zeroPowerSeries_PS();
	PowerSeries_t* sub = subPowerSeries_PS(ps,ps);
	PowerSeries_t* sum = addPowerSeries_PS(ps,ps);
	Poly_ptr  homogSum = homogPart_PS(0,sum);
	Poly_ptr  homogSub = homogPart_PS(0,sub);
	PowerSeries_t* multiply = multiplyPowerSeries_PS(ps,ps);
	Poly_ptr  homogMultiply = homogPart_PS(0,multiply);
	if ((isZero_AA(homogSum)) && (isZero_AA(homogSub)) && ((isZero_AA(homogMultiply)))) {
		fprintf(stderr, "PowerSeries Test 27: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 27: \t\t\t\t FAILED \n");
		exit(1);
	}


	destroyPowerSeries_PS(multiply);
	destroyPowerSeries_PS(sum);
	destroyPowerSeries_PS(sub);
	destroyPowerSeries_PS(ps);
	return 0;
}



int test28()
{

	PowerSeries_t* ps = zeroPowerSeries_PS();
	PowerSeries_t* sub = subPowerSeries_PS(ps,ps);
	PowerSeries_t* sum = addPowerSeries_PS(ps,ps);
	Poly_ptr  homogSum = homogPart_PS(10,sum);
	Poly_ptr  homogSub = homogPart_PS(20,sub);
	PowerSeries_t* multiply = multiplyPowerSeries_PS(ps,ps);
	Poly_ptr  homogMultiply = homogPart_PS(10,multiply);
	if ((isZero_AA(homogSum)) && (isZero_AA(homogSub)) && ((isZero_AA(homogMultiply)))) {
		fprintf(stderr, "PowerSeries Test 28: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 28: \t\t\t\t FAILED \n");
		exit(1);
	}


	destroyPowerSeries_PS(multiply);
	destroyPowerSeries_PS(sum);
	destroyPowerSeries_PS(sub);
	destroyPowerSeries_PS(ps);
    return 0;
}


/**
*
*/
int test29()
{
    Poly_ptr a = generate_altarr_var_defined("1",sym, 2);
    Poly_ptr b = generate_altarr_var_defined("1 + x + y + x^2 + y^2",sym, 2);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = subPowerSeries_PS(a_ps,b_ps);
    Poly_ptr  homog = homogPart_PS(2,c);
    Poly_ptr  testpoly = generate_altarr_var_defined("-x^2 - y^2", sym, 2);
    testpoly = subPolynomials_AA_inp(testpoly,homog,2);
    if (isZero_AA(testpoly)) {
		fprintf(stderr, "PowerSeries Test 29: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 29: \t\t\t\t FAILED \n");
		exit(1);
	}

    freePolynomial_AA(a);
    freePolynomial_AA(b);
    freePolynomial_AA(testpoly);

    destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a_ps);
    destroyPowerSeries_PS (b_ps);
    return 0;
}


int test30()
{
	Poly_ptr  a = NULL;
	Poly_ptr  b = generate_altarr_var_defined("z^5-y^5+x*t^2*z+x^4*t+1", syms, 4);
    PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t* c = subPowerSeries_PS(a_ps,b_ps);
    Poly_ptr homog = homogPart_PS(5, c);
    Poly_ptr testpoly = generate_altarr_var_defined("-z^5+y^5-x^4*t", syms, 4);
    testpoly = subPolynomials_AA_inp(testpoly,homog,4);

    if (isZero_AA(testpoly)) {
		fprintf(stderr, "PowerSeries Test 30: \t\t\t\t PASSED \n");
	}
	else {
		fprintf(stderr, "PowerSeries Test 30: \t\t\t\t FAILED \n");
		exit(1);
	}


	freePolynomial_AA(b);
	freePolynomial_AA(testpoly);

	destroyPowerSeries_PS(a_ps);
	destroyPowerSeries_PS(b_ps);
	destroyPowerSeries_PS(c);

	return 0;
}

int test31()
{

	Poly_ptr  b = generate_altarr_var_defined("z^5-y^5+x*t^2*z+x^4*t+1", syms, 4);

    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t*  c_ps = negatePowerSeries_PS(b_ps);

    homogPart_PS(5, c_ps);

    Poly_ptr poly_1 = generate_altarr_var_defined("z^5-y^5+x^4*t", syms, 4);
    poly_1 = addPolynomials_AA_inp(poly_1, c_ps->polys[5],4);
    Poly_ptr poly_2 = generate_altarr_var_defined("x*t^2*z", syms, 4);
    poly_2 = addPolynomials_AA_inp(poly_2,c_ps->polys[4],4);
    Poly_ptr poly_3 = generate_altarr_var_defined("1", syms, 4);
    poly_3 = addPolynomials_AA_inp(poly_3,c_ps->polys[0],4);

    if (isZero_AA(poly_1) && isZero_AA(poly_2) && isZero_AA(poly_3)) {
		fprintf(stderr, "PowerSeries Test 31: \t\t\t\t PASSED \n");
    } else {
		fprintf(stderr, "PowerSeries Test 31: \t\t\t\t FAILED \n");
		exit(1);
    }



    destroyPowerSeries_PS(c_ps);
    destroyPowerSeries_PS(b_ps);
    freePolynomial_AA(b);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    return 0;

}

int test32()
{

	Poly_ptr  b = NULL;
    PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
    PowerSeries_t*  c_ps = negatePowerSeries_PS(b_ps);
    b = addPolynomials_AA_inp(b, c_ps->polys[0],4);

    if (isZero_AA(b)) {
		fprintf(stderr, "PowerSeries Test 32: \t\t\t\t PASSED \n");
    } else {
		fprintf(stderr, "PowerSeries Test 32: \t\t\t\t FAILED \n");
    	exit(1);
    }

    destroyPowerSeries_PS(c_ps);
    destroyPowerSeries_PS(b_ps);
    freePolynomial_AA(b);
    return 0;

}


int test33()
{
	// Poly_ptr a = generate_altarr_var_defined("1",sym, 2);
	Poly_ptr b = generate_altarr_var_defined("1-x",sym, 2);
	// PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	TestingPowerSeriesInversion(b_ps, 30, sym, 2);

	// freePolynomial_AA(a);
	freePolynomial_AA(b);
	// destroyPowerSeries_PS(a_ps);
	destroyPowerSeries_PS(b_ps);

	fprintf(stderr, "PowerSeries Test 33: \t\t\t\t PASSED \n");
	return 0;
}


int test34()
{
	Poly_ptr a = generate_altarr_var_defined("y",sym, 2);
	Poly_ptr b = generate_altarr_var_defined("1-x",sym, 2);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	TestingPowerSeriesDivision(a_ps, b_ps,10, sym, 2);

	freePolynomial_AA(a);
	freePolynomial_AA(b);
	destroyPowerSeries_PS(a_ps);
	destroyPowerSeries_PS(b_ps);

	fprintf(stderr, "PowerSeries Test 34: \t\t\t\t PASSED \n");

	return 0;
}


int test35()
{

	Poly_ptr a = generate_altarr_var_defined("-2*x*t+x^2+y*z+t^4+1",syms, 4);
	Poly_ptr b = generate_altarr_var_defined("1-t-x",syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	TestingPowerSeriesDivision(b_ps, a_ps,20, syms, 4);

	freePolynomial_AA(a);
	freePolynomial_AA(b);
	destroyPowerSeries_PS(a_ps);
	destroyPowerSeries_PS(b_ps);

	fprintf(stderr, "PowerSeries Test 35: \t\t\t\t PASSED \n");

	return 0;
}


int test36()
{

	Poly_ptr a = generate_altarr_var_defined("-2*x*t+x^2+y*z+t^4+1",syms, 4);
	Poly_ptr b = generate_altarr_var_defined("1-t*x^3-x^6+2*y",syms, 4);
	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
	TestingPowerSeriesDivision(b_ps, a_ps,20, syms, 4);

	freePolynomial_AA(a);
	freePolynomial_AA(b);
	destroyPowerSeries_PS(a_ps);
	destroyPowerSeries_PS(b_ps);

	fprintf(stderr, "PowerSeries Test 36: \t\t\t\t PASSED \n");
	return 0;
}


// int test37()
// {
// 	Poly_ptr a = generate_altarr_var_defined("1+x+y+z",vars, 3);
// 	Poly_ptr b = generate_altarr_var_defined("1",vars, 3);
// 	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
// 	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
// 	PowerSeries_t* c = dividePowerSeries_PS(b_ps,a_ps);
// 	unsigned long long start;
// 	_startTimer(&start);
// 	Poly_ptr homog = homogPart_PS(600,c);
// 	float elapsedTime;
// 	_stopTimer(&start, &elapsedTime);
// 	fprintf(stderr, "Time sys37: %f\n", elapsedTime);
//     return 0;

// }


// int test38()
// {
// 	Poly_ptr a = generate_altarr_var_defined("1+x+y",sym, 2);
// 	Poly_ptr b = generate_altarr_var_defined("1",sym, 2);
// 	PowerSeries_t*  a_ps = convertPolyToPowerSeries_PS(a);
// 	PowerSeries_t*  b_ps = convertPolyToPowerSeries_PS(b);
// 	PowerSeries_t* c = dividePowerSeries_PS(b_ps,a_ps);
// 	for (int indx = 100; indx <= 2000; indx = indx +100)
// 	TestingPowerSeriesDivision(b_ps, a_ps, indx, sym, 2);
//     return 0;

// }
