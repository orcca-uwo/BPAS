

#include "PowerSeriesZ_TestCases.h"
#include "../../include/Parser/bpas_parser.h"
#include "../../include/Utils/Unix_Timer.h"

int main(int argc, char** argv) {

	// const char* vars[] = {"x", "y", "z"};

 //    PolyZ_ptr  polynom = generate_altarrZ_var_defined("1+x+y+z", vars, 3);
 //    PowerSeriesZ_t*  a = convertPolyToPowerSeries_PSZ(polynom);
 //    PowerSeriesZ_t*  b = inversePowerSeries_PSZ(a);

 //    unsigned long long start;
 //    _startTimer(&start);
 //    PolyZ_ptr p = homogPart_PSZ(400, b);
 //    float elapsed;
 //    _stopTimer(&start, &elapsed);

 //    // printPoly_AA(stderr, p, vars, 2);
 //    fprintf(stderr, "\np->size: %d\n", p->size);
 //    fprintf(stderr, "\ntime: %9.4f\n", elapsed);
 //    // int longCount = 0;
 //    // for (int i = 0; i < p->size; ++i) {
 //    // 	if (mpz_sizeinbase(mpq_numref(p->elems[i].coef), 2) > 63) {
 //    // 		++longCount;
 //    // 	}
 //    // }
 //    // fprintf(stderr, "\nlongCounf: %d\n", longCount);


 //    freePolynomial_AAZ(polynom);

 //    destroyPowerSeries_PSZ (a);
 //    destroyPowerSeries_PSZ (b);

	// return 0;

	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	test8();
	test9();
	test10();
	test11();
	test12();
	test13();
	test14();
	test15();
	test16();
	test17();
	test18();
	test19();
	test20();
	test21();
	test22();
	test23();
	test24();
	test25();
	test26();
	test27();
	test28();
	test29();
	test30();
	test31();
	test32();
	test33();
	test34();
	test35();
	test36();

	return 0;
}

