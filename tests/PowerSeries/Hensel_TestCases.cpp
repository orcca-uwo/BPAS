
#include <gmpxx.h>

#include "Hensel_TestCases.hpp"

#include "../../include/PowerSeries/PowerSeries.h"
#include "../../include/PowerSeries/UnivariatePolynomialOverPowerSeries.h"
#include "../../include/PowerSeries/UPOPS_Weierstrass.h"
#include "../../include/PowerSeries/UPOPS_Hensel.hpp"
#include "../../include/Utils/Unix_Timer.h"
#include "../../include/Parser/bpas_parser.h"
#include <time.h>


/**
 * f : the input upop
 * C : the given roots
 * bound : the upper bound for the accuracy
 */
void TestingFactorizationViaHensel(Upops_t* f, mpq_t* C, int n, int bound, const char** sym, int nvar ) {

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "F:\n");
	print_UPOPS(stderr, f , sym);
	fprintf(stderr, "\n");
#endif

	Upops_t** facts_out = NULL;
	HenselFactorization_UPOPS(f, C, n, &facts_out);

	Upops_t* prod = facts_out[0];
	reserve_UPOPS(prod);
	for (int i = 1; i <= n-1; ++i) {
		Upops_t* tmp = multiplyUnivariatePolyOverPowerSeries_UPOPS(prod, facts_out[i]);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
		prod = tmp;
	}
	Upops_t* sub =  subUnivariatePolyOverPowerSeries_UPOPS(prod,f);

	for (int d = 1; d <= bound; ++d) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
		fprintf(stderr, "\n");
		fprintf(stderr, "In degree: %d\n ", d);
		fprintf(stderr, "\n");
#endif
		for (int j = 0; j <= n-1; ++j) {
			updateToDeg_UPOPS(d, facts_out[j]);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
			fprintf(stderr, "\n");
			fprintf(stderr, "The upops at index: %d\n",j);
			fprintf(stderr, "\n");
			fprintf(stderr, "\n");
			print_UPOPS(stderr, facts_out[j] , sym);
			fprintf(stderr, "\n");
#endif
		}

		for (int k = 0; k <= sub->deg; ++k) {
			Poly_ptr polyPart = polynomialPart_UPOPS(d, sub);
			if (!isZero_AA(polyPart)) {
				fprintf(stderr, "Testing Factorization Via Hensel: ERROR\n");
				exit(1);
			}
			freePolynomial_AA(polyPart);
		}

	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	for (int i = 0; i < n; ++i) {
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(facts_out[i]);
	}
	free(facts_out);

}


#if defined(WITH_NTL) && WITH_NTL

/**
 * f : the input upop
 * C : the given roots
 * bound : the upper bound for the accuracy
 */
void TestingFactorizationViaHensel(Upops_t* f, int bound, const char** sym, int nvar) {

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "F:\n");
	print_UPOPS(stderr, f , sym);
	fprintf(stderr, "\n");
#endif

	int n;
	Upops_t** facts_out = NULL;
	HenselFactorization_UPOPS(f, &facts_out, &n);

	Upops_t* prod = facts_out[0];
	reserve_UPOPS(prod);
	for (int i = 1; i <= n-1; ++i) {
		Upops_t* tmp = multiplyUnivariatePolyOverPowerSeries_UPOPS(prod, facts_out[i]);
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
		prod = tmp;
	}
	Upops_t* sub =  subUnivariatePolyOverPowerSeries_UPOPS(prod,f);

	for (int d = 1; d <= bound; ++d) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
		fprintf(stderr, "\n");
		fprintf(stderr, "In degree: %d\n ", d);
		fprintf(stderr, "\n");
#endif
		for (int j = 0; j <= n-1; ++j) {
			updateToDeg_UPOPS(d, facts_out[j]);
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
			fprintf(stderr, "\n");
			fprintf(stderr, "The upops at index: %d\n",j);
			fprintf(stderr, "\n");
			fprintf(stderr, "\n");
			print_UPOPS(stderr, facts_out[j] , sym);
			fprintf(stderr, "\n");
#endif
		}

		for (int k = 0; k <= sub->deg; ++k) {
			Poly_ptr polyPart = polynomialPart_UPOPS(d, sub);
			if (!isZero_AA(polyPart)) {
				fprintf(stderr, "Testing Factorization Via Hensel: ERROR\n");
				exit(1);
			}
			freePolynomial_AA(polyPart);
		}
	}

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	for (int i = 0; i < n; ++i) {
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(facts_out[i]);
	}
	free(facts_out);
}

#endif

/**
 * c_r : the root
 * f : the input upop
 * sym :
 */
void TestingTaylorShift(mpq_t c_r, Upops_t* f, int n, const char** sym) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
    fprintf(stderr, "The input upops:\n");
	print_UPOPS(stderr, f, sym);
	fprintf(stderr, "\n");
#endif
	Upops_t* W_out = taylorShift_UPOPS(f, c_r);
	mpq_t negateVal;
	mpq_init(negateVal);
	mpq_neg(negateVal, c_r);
	Upops_t* W_new_out = taylorShift_UPOPS(W_out, negateVal);
	mpq_clear(negateVal);

	Upops_t* sub = subUnivariatePolyOverPowerSeries_UPOPS(f, W_new_out);
	for (int i = 0; i < n; ++i) {
        Poly_ptr polyPart = polynomialPart_UPOPS(i, sub);
        // Poly_ptr outPoly = polynomialPart_UPOPS(W_out, i);
    	// printPoly_AA(stderr, outPoly, sym, outPoly->nvar);
        if (!isZero_AA(polyPart)) {
        	fprintf(stderr, "Testing Taylor Shift ERROR\n");
        	printPoly_AA(stderr, polyPart, sym, polyPart->nvar);
        	exit(1);
        }
        freePolynomial_AA(polyPart);
	}


	destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(W_new_out);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(W_out);
	// fprintf(stderr, "After shift and shift back:\n");
	// updateToDeg_UPOPS(n, W_new_out);
	// print_UPOPS(stderr, W_new_out, sym);
	// fprintf(stderr, "\n");

}


void test1() {
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

    mpq_t c;
    mpq_init(c);
    mpq_set_si(c, -3, 7);
	TestingTaylorShift(c, upop, 100,sym);
	mpq_clear(c);

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

    fprintf(stderr, "Hensel Test  1: \t\t\t\t PASSED \n");
}


/**
 */
void test2()
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

	mpq_t c;
	mpq_init(c);
	mpq_set_si(c, 45, 1);
	TestingTaylorShift(c, upop, 100, sym);
	mpq_clear(c);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Hensel Test  2: \t\t\t\t PASSED \n");

}


/**
 */
void test3()
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

    mpq_t c;
    mpq_init(c);
    mpq_set_ui(c, 5, 1);
	TestingTaylorShift(c, upop, 30,sym);

	mpq_clear(c);

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

    fprintf(stderr, "Hensel Test  3: \t\t\t\t PASSED \n");
}


/**
 */
void test4()
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

    mpq_t c;
    mpq_init(c);
    mpq_set_si(c, 17, 2);

    TestingTaylorShift(c, upop, 50, sym);

	mpq_clear(c);

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

    fprintf(stderr, "Hensel Test  4: \t\t\t\t PASSED \n");
}


/**
 *  Computes the homogPart of degree 6
 */
void test5()
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

	mpq_t c;
	mpq_init(c);
	mpq_set_si(c, 43, 7);
	TestingTaylorShift(c, upop, 50, sym);
	mpq_clear(c);

    destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
    freePolynomial_AA(poly_1);
    freePolynomial_AA(poly_2);
    freePolynomial_AA(poly_3);

    destroyPowerSeries_PS(ps_1);
    destroyPowerSeries_PS(ps_2);
    destroyPowerSeries_PS(ps_3);
    destroyPowerSeries_PS(quo_1);

    fprintf(stderr, "Hensel Test  5: \t\t\t\t PASSED \n");

}

/**
 *  Taylor shift initialization
 */
void test6()
{
	const char* sym[] = {"z", "y", "x"};

	Poly_ptr  poly_0 = generate_altarr_var_defined("x+1", sym, 3);
	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(poly_0);

	mpq_t c_r;
	mpq_init(c_r);
	mpz_set_ui(mpq_numref(c_r), 3ul);
	mpz_set_ui(mpq_denref(c_r), 2ul);

	TestingTaylorShift(c_r, upop, 50, sym);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	freePolynomial_AA(poly_0);
	mpq_clear(c_r);

    fprintf(stderr, "Hensel Test  6: \t\t\t\t PASSED \n");
}


/**
 *  Taylor shift initialization
 */
void test7()
{
	const char* sym[] = {"z", "y", "x"};

	Poly_ptr  poly_0 = generate_altarr_var_defined("x+1", sym, 3);
	Poly_ptr  poly_1 = generate_altarr_var_defined("1-y-x^2", sym, 3);
	PowerSeries_t* ps_0 = convertPolyToPowerSeries_PS(poly_0);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* psArr[2];
	psArr[0] = ps_0;
	psArr[1] = ps_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(psArr,2);

	mpq_t c_r;
	mpq_init(c_r);
	mpz_set_ui(mpq_numref(c_r), 3ul);
	mpz_set_ui(mpq_denref(c_r), 2ul);

	TestingTaylorShift(c_r, upop, 50, sym);

	mpq_clear(c_r);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	freePolynomial_AA(poly_0);
	freePolynomial_AA(poly_1);
	destroyPowerSeries_PS(ps_0);
	destroyPowerSeries_PS(ps_1);

    fprintf(stderr, "Hensel Test  7: \t\t\t\t PASSED \n");

}


/**
 *  Taylor shift initialization
 */
void test8()
{
	const char* sym[] = {"z", "y", "x"};

	Poly_ptr  poly_0 = generate_altarr_var_defined("x+1+z^3", sym, 3);
	Poly_ptr  poly_1 = generate_altarr_var_defined("1-y-x^2+x^3+y*z^2+z^3", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("-2+x*y+x^3", sym, 3);
	PowerSeries_t* ps_0 = convertPolyToPowerSeries_PS(poly_0);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* psArr[3];
	psArr[0] = ps_0;
	psArr[1] = ps_1;
	psArr[2] = ps_2;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(psArr,3);
	mpq_t c_r;
	mpq_init(c_r);
	mpz_set_ui(mpq_numref(c_r), 3ul);
	mpz_set_ui(mpq_denref(c_r), 2ul);

	TestingTaylorShift(c_r, upop, 50, sym);

	mpq_clear(c_r);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyPowerSeries_PS(ps_0);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	freePolynomial_AA(poly_0);
	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);

    fprintf(stderr, "Hensel Test  8: \t\t\t\t PASSED \n");
}


/**
 *  Taylor shift asking for more terms via calling the generator
 */
void test9()
{
	const char* sym[] = {"z", "y", "x"};

	Poly_ptr  poly_0 = generate_altarr_var_defined("x^6+1", sym, 3);
	Poly_ptr  poly_1 = generate_altarr_var_defined("1-y-x^2+x^3+y*x^5", sym, 3);
    Poly_ptr  poly_2 = generate_altarr_var_defined("-2+x^3*y^3+x^3", sym, 3);
	PowerSeries_t* ps_0 = convertPolyToPowerSeries_PS(poly_0);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* psArr[3];
	psArr[0] = ps_0;
	psArr[1] = ps_1;
	psArr[2] = ps_2;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(psArr,3);

	mpq_t c_r;
	mpq_init(c_r);
	mpz_set_ui(mpq_numref(c_r), 3ul);
	mpz_set_ui(mpq_denref(c_r), 2ul);

	TestingTaylorShift(c_r, upop, 50, sym);

	mpq_clear(c_r);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyPowerSeries_PS(ps_0);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	freePolynomial_AA(poly_0);
	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);

    fprintf(stderr, "Hensel Test  9: \t\t\t\t PASSED \n");


}

void test10()
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4+(8-x)*z^3+18*z^2-27", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t c;
	mpq_init(c);
	mpz_set_ui(mpq_numref(c), 1ul);
	mpz_set_ui(mpq_denref(c), 1ul);
	TestingTaylorShift(c,f,20,sym);
	mpq_clear(c);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

    fprintf(stderr, "Hensel Test 10: \t\t\t\t PASSED \n");
}


void test11()
{
	const char* sym[] = {"z", "y", "x"};

	Poly_ptr  poly_0 = generate_altarr_var_defined("x^6+1", sym, 3);
	Poly_ptr  poly_1 = generate_altarr_var_defined("1-y-x^2+x^3+y*x^5", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("-2+x^3*y^3+x^3", sym, 3);
	PowerSeries_t* ps_0 = convertPolyToPowerSeries_PS(poly_0);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* psArr[3];
	psArr[0] = ps_0;
	psArr[1] = ps_1;
	psArr[2] = ps_2;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(psArr,3);
	mpq_t c_r;
	mpq_init(c_r);
	mpz_set_ui(mpq_numref(c_r), 3ul);
	mpz_set_ui(mpq_denref(c_r), 2ul);
	TestingTaylorShift(c_r,upop,6,sym);

	mpq_clear(c_r);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyPowerSeries_PS(ps_0);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	freePolynomial_AA(poly_0);
	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);

    fprintf(stderr, "Hensel Test 11: \t\t\t\t PASSED \n");
}

/**
 * sys1: Factorization via Hensel 1000
 */
void test12(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^2-x*z-1", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);

	mpq_t C[2];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 1ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -1, 1);

	// unsigned long long start;
	// _startTimer(&start);
	// TestingFactorizationViaHensel(f, degree, sym, 3);

	TestingFactorizationViaHensel(f, C, 2, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys1: %f\n", elapsedTime);


	mpq_clear(C[0]);
	mpq_clear(C[1]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 12: \t\t\t\t PASSED \n");
}

/**
 * sys2: Factorization via Hensel 1000
 */
void test13(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^2+(2-x)*z-x", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[2];
    mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
    mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -2, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 2,degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys2: %f\n", elapsedTime);


	mpq_clear(C[0]);
	mpq_clear(C[1]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 13: \t\t\t\t PASSED \n");

}



/**
 * sys3 : Factorization via Hensel 1000
 */
void test14(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("2+(-1-x)*z+(-2+x)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[3];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 1ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpz_set_ui(mpq_numref(C[1]), 2ul);
	mpz_set_ui(mpq_denref(C[1]), 1ul);
	mpq_init(C[2]);
	mpq_set_si(C[2], -1, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 3, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys3: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);
	mpq_clear(C[2]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 14: \t\t\t\t PASSED \n");
}

/**
 * Factorization via Hensel
 */
void test15()
{
	const char* sym[] = {"z", "y", "x"};

	AltArr_t* a = generate_altarr_var_defined("z^3-2*z^2 + 3*x^2*z^2 + x*z - z + 2 + x", sym, 3);


	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[3];
	mpq_init(C[0]);
	mpq_init(C[1]);
	mpq_init(C[2]);
	mpz_set_ui(mpq_numref(C[0]), 1ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpz_set_si(mpq_numref(C[1]), -1l);
	mpz_set_si(mpq_numref(C[2]), 2l);


	TestingFactorizationViaHensel(f, C, 3, 50, sym, 3);


	mpq_clear(C[0]);
	mpq_clear(C[1]);
	mpq_clear(C[2]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 15: \t\t\t\t PASSED \n");
}


/**
 * Factorization via Hensel
 */
void test16(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("-1+y", sym, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("x", sym, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
	PowerSeries_t* ps_array[3];
	ps_array[0] = quo_1;
	ps_array[1] = ps_3;
	ps_array[2] =  ps_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,3);
	mpq_t C[2];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 1ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -1, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(upop, C, 2, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys16: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(quo_1);
	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);

	fprintf(stderr, "Hensel Test 16: \t\t\t\t PASSED \n");
}

/**
 * Factorization via Hensel (Maple is very slow on this example (because it computes complex roots!
 * but we only compute one root) (Checcking the correctness)
 */
// void test17()
// {
// 	const char* sym[] = {"z", "y", "x"};
// 	AltArr_t* a = generate_altarr_var_defined("-1-x*z+y*z^2-y*z^3+(x+y)*z^4+z^5", sym, 3);
// 	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
// 	mpq_t* C = (mpq_t*) malloc(sizeof(mpq_t)*1);
// 	mpq_init(C[0]);
// 	mpz_set_ui(mpq_numref(C[0]), 1ul);
// 	mpz_set_ui(mpq_denref(C[0]), 1ul);
// 	//TestingFactorizationViaHensel(f, C, 1, 10, sym, 3);
// 	Upops_t** facts_out = NULL;
// 	HenselFactorization_UPOPS(f, C, 1, &facts_out);
// 	homogPart_PS(22, facts_out[0]->data[0]);
// 	fprintf(stderr, "First factor\n");
// 	print_UPOPS(stderr, facts_out[0], sym);
// 	fprintf(stderr, "\n");

// }

/**
 * sys4 : Factorization via Hensel 1000
 */
void test18(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y+z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[2];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -1, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 2, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys4: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 18: \t\t\t\t PASSED \n");
}


/**
 * sys6 : Factorization via Hensel 1000
 */
void test19(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("8+x+(12+y^2)*z+(y+6)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[1];
	mpq_init(C[0]);
	mpq_set_si(C[0], -2, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 1, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys6: %f\n", elapsedTime);


	mpq_clear(C[0]);
	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 19: \t\t\t\t PASSED \n");
}


/**
 * sys10 : Factorization via Hensel:qubic 1000
 */
void test20(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y^2+x^2+(y+1)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[2];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -1, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 2,degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys10: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 20: \t\t\t\t PASSED \n");
}


/**
 * sys7 : Factorization via Hensel 1000
 */
void test21(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y+(1/2)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[2];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -1, 2);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 2, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys7: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);
	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 21: \t\t\t\t PASSED \n");


}

/**
 * sys8 : Factorization via Hensel:qubic  1000
 */
void test22(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^3+(3+x)*z^2-(6+y^2)*z-8*(x+y+1)", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[3];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 2ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -1, 1);
	mpq_init(C[2]);
	mpq_set_si(C[2], -4, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 3, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys8: %f\n", elapsedTime);


	mpq_clear(C[0]);
	mpq_clear(C[1]);
	mpq_clear(C[2]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 22: \t\t\t\t PASSED \n");

}


/**
 * sys9 : Factorization via Hensel: quartic  1000
 */
void test23(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4-(26+y^2)*z^2+25*(x+1)", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[4];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 1ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpz_set_ui(mpq_numref(C[1]), 5ul);
	mpz_set_ui(mpq_denref(C[1]), 1ul);
	mpq_init(C[2]);
	mpq_set_si(C[2], -1, 1);
	mpq_init(C[3]);
	mpq_set_si(C[3], -5, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 4, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys9: %f\n", elapsedTime);



	mpq_clear(C[0]);
	mpq_clear(C[1]);
	mpq_clear(C[2]);
	mpq_clear(C[3]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 23: \t\t\t\t PASSED \n");

}


/**
 * sys5 : Factorization via Hensel: degree 6 1000
 */
void test24(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4+(8-x)*z^3+18*z^2-27", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[2] ;
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 1ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpq_set_si(C[1], -3, 1);
	// unsigned long long start;
	// _startTimer(&start);

	TestingFactorizationViaHensel(f, C, 2, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys5: %f\n", elapsedTime);


	mpq_clear(C[0]);
	mpq_clear(C[1]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 24: \t\t\t\t PASSED \n");

}


void test25(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^3-5*z^2+(6+x)*z", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[3];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpz_set_ui(mpq_numref(C[1]), 2ul);
	mpz_set_ui(mpq_denref(C[1]), 1ul);
	mpq_init(C[2]);
	mpz_set_ui(mpq_numref(C[2]), 3ul);
	mpz_set_ui(mpq_denref(C[2]), 1ul);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 3, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys11: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);
	mpq_clear(C[2]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 25: \t\t\t\t PASSED \n");
}


void test26(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4+(4+x)*z^3-7*z^2-10*z", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[4];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpz_set_ui(mpq_numref(C[1]), 2ul);
	mpz_set_ui(mpq_denref(C[1]), 1ul);
	mpq_init(C[2]);
	mpq_set_si(C[2], -5, 1);
	mpq_init(C[3]);
	mpq_set_si(C[3], -1, 1);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 4, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys12: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);
	mpq_clear(C[2]);
	mpq_clear(C[3]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 26: \t\t\t\t PASSED \n");

}

void test27(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^2-9*z+x", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	mpq_t C[2];
	mpq_init(C[0]);
	mpz_set_ui(mpq_numref(C[0]), 0ul);
	mpz_set_ui(mpq_denref(C[0]), 1ul);
	mpq_init(C[1]);
	mpz_set_ui(mpq_numref(C[1]), 9ul);
	mpz_set_ui(mpq_denref(C[1]), 1ul);
	// unsigned long long start;
   // _startTimer(&start);
	TestingFactorizationViaHensel(f, C, 2, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys13: %f\n", elapsedTime);

	mpq_clear(C[0]);
	mpq_clear(C[1]);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 27: \t\t\t\t PASSED \n");
}



#if defined(WITH_NTL) && WITH_NTL

/**
 * sys1: Factorization via Hensel 1000
 */
void test28(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^2-x*z-1", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);


	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);

	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys1: %f\n", elapsedTime);


	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 28: \t\t\t\t PASSED \n");
}

/**
 * sys2: Factorization via Hensel 1000
 */
void test29(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^2+(2-x)*z-x", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys2: %f\n", elapsedTime);


	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 29: \t\t\t\t PASSED \n");

}



/**
 * sys3 : Factorization via Hensel 1000
 */
void test30(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("2+(-1-x)*z+(-2+x)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys3: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 30: \t\t\t\t PASSED \n");
}

/**
 * Factorization via Hensel
 */
void test31(int degree)
{
	const char* sym[] = {"z", "y", "x"};

	AltArr_t* a = generate_altarr_var_defined("z^3-2*z^2 + 3*x^2*z^2 + x*z - z + 2 + x", sym, 3);


	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);

	TestingFactorizationViaHensel(f, degree, sym, 3);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 31: \t\t\t\t PASSED \n");
}


/**
 * Factorization via Hensel
 */
void test32(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	Poly_ptr  poly_1 = generate_altarr_var_defined("1", sym, 3);
	Poly_ptr  poly_2 = generate_altarr_var_defined("-1+y", sym, 3);
	Poly_ptr  poly_3 = generate_altarr_var_defined("x", sym, 3);
	PowerSeries_t* ps_1 = convertPolyToPowerSeries_PS(poly_1);
	PowerSeries_t* ps_2 = convertPolyToPowerSeries_PS(poly_2);
	PowerSeries_t* ps_3 = convertPolyToPowerSeries_PS(poly_3);
	PowerSeries_t* quo_1 = dividePowerSeries_PS(ps_1, ps_2);
	PowerSeries_t* ps_array[3];
	ps_array[0] = quo_1;
	ps_array[1] = ps_3;
	ps_array[2] =  ps_1;
	Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,3);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(upop, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys16: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(upop);
	destroyPowerSeries_PS(ps_1);
	destroyPowerSeries_PS(ps_2);
	destroyPowerSeries_PS(ps_3);
	destroyPowerSeries_PS(quo_1);
	freePolynomial_AA(poly_1);
	freePolynomial_AA(poly_2);
	freePolynomial_AA(poly_3);

	fprintf(stderr, "Hensel Test 32: \t\t\t\t PASSED \n");
}

/**
 * Factorization via Hensel (Maple is very slow on this example (because it computes complex roots!
 * but we only compute one root) (Checcking the correctness)
 */
// void test17()
// {
// 	const char* sym[] = {"z", "y", "x"};
// 	AltArr_t* a = generate_altarr_var_defined("-1-x*z+y*z^2-y*z^3+(x+y)*z^4+z^5", sym, 3);
// 	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
// 	mpq_t* C = (mpq_t*) malloc(sizeof(mpq_t)*1);
// 	mpq_init(C[0]);
// 	mpz_set_ui(mpq_numref(C[0]), 1ul);
// 	mpz_set_ui(mpq_denref(C[0]), 1ul);
// 	//TestingFactorizationViaHensel(f, C, 1, 10, sym, 3);
// 	Upops_t** facts_out = NULL;
// 	HenselFactorization_UPOPS(f, C, 1, &facts_out);
// 	homogPart_PS(22, facts_out[0]->data[0]);
// 	fprintf(stderr, "First factor\n");
// 	print_UPOPS(stderr, facts_out[0], sym);
// 	fprintf(stderr, "\n");

// }

/**
 * sys4 : Factorization via Hensel 1000
 */
void test33(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y+z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys4: %f\n", elapsedTime);

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 33: \t\t\t\t PASSED \n");
}


/**
 * sys6 : Factorization via Hensel 1000
 */
void test34(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("8+x+(12+y^2)*z+(y+6)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys6: %f\n", elapsedTime);

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 34: \t\t\t\t PASSED \n");
}


/**
 * sys10 : Factorization via Hensel:qubic 1000
 */
void test35(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y^2+x^2+(y+1)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f,degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys10: %f\n", elapsedTime);

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 35: \t\t\t\t PASSED \n");
}


/**
 * sys7 : Factorization via Hensel 1000
 */
void test36(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("y+(1/2)*z^2+z^3", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys7: %f\n", elapsedTime);

	freePolynomial_AA(a);
	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);

	fprintf(stderr, "Hensel Test 36: \t\t\t\t PASSED \n");


}

/**
 * sys8 : Factorization via Hensel:qubic  1000
 */
void test37(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^3+(3+x)*z^2-(6+y^2)*z-8*(x+y+1)", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys8: %f\n", elapsedTime);


	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 37: \t\t\t\t PASSED \n");

}


/**
 * sys9 : Factorization via Hensel: quartic  1000
 */
void test38(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4-(26+y^2)*z^2+25*(x+1)", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys9: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 38: \t\t\t\t PASSED \n");

}


/**
 * sys5 : Factorization via Hensel: degree 6 1000
 */
void test39(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4+(8-x)*z^3+18*z^2-27", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);

	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys5: %f\n", elapsedTime);


	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 39: \t\t\t\t PASSED \n");

}


void test40(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^3-5*z^2+(6+x)*z", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f,  degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys11: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 40: \t\t\t\t PASSED \n");
}


void test41(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^4+(4+x)*z^3-7*z^2-10*z", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
	// _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys12: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 41: \t\t\t\t PASSED \n");

}

void test42(int degree)
{
	const char* sym[] = {"z", "y", "x"};
	AltArr_t* a = generate_altarr_var_defined("z^2-9*z+x", sym, 3);
	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	// unsigned long long start;
   // _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, 3);
	// float elapsedTime;
	// _stopTimer(&start, &elapsedTime);
	// fprintf(stderr, "Time sys13: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	fprintf(stderr, "Hensel Test 42: \t\t\t\t PASSED \n");
}

#endif