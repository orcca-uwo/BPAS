

#ifndef _POWERSERIES_TESTCASES_
#define _POWERSERIES_TESTCASES_

#include "../../include/PowerSeries/PowerSeries.h"

#ifdef __cplusplus
extern "C" {
#endif

void TestingPowerSeriesDivision(PowerSeries_t* f, PowerSeries_t* g, int bound, const char** sym, int n);


void test1();

void test2();

void test3();

void test4();

void test5();

void test6();

void test7();

/**
 * Note that multiplication process truncates the result up to the min degree of ancestors
 * in the following the truncation degree is 3. Thus the array of polys does not contain terms
 * of degree at least 4. But when we call homogPart_PS and ask for terms of degree 5 and look at
 * polys again we can see that this time the array polys has been updated and contains terms of degree 5
 */
void test8();

void test9();

void test10();

void test11();

void test12();

void test13();

void test14();

void test15();

void test16();

void test17();

void test18();

void test19();

void test20();

void test21();

void test22();

void test23();

void test24();

void test25();

void test26();

void test27();

void test28();

void test29();

void test30();

void test31();

void test32();

void test33();

void test34();

void test35();

void test36();

void test37();

void test38();

void test39();

#ifdef __cplusplus
}
#endif


#endif