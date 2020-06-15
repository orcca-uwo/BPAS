#ifndef _WEIERSTRASS_TESTCASES_
#define _WEIERSTRASS_TESTCASES_

#include "../../include/PowerSeries/UPOPS_Weierstrass.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * f : a upops
 * bound : the maximum integer for which homogPart_PS gets called
 * sym : a list of variables
 * nvar : the number of variables
 */
void TestingWeiestrassPreparation(Upops_t* f, int bound, const char** sym);

void test1(int testDeg);

void test2(int testDeg);

void test3(int testDeg);

void test4(int testDeg);

void test5(int testDeg);

void test6(int testDeg);

void test7(int testDeg);

void test8(int testDeg);

void test9(int testDeg);

void test10(int testDeg);

void test11(int testDeg);

void test12(int testDeg);

void test13(int testDeg);

void test14(int testDeg);

void test15(int testDeg);

void test16(int testDeg);

void test17(int testDeg);

void test18(int testDeg);

void test19(int testDeg);

void test20(int testDeg);

void test21(int testDeg);

void test22(int testDeg);

#ifdef __cplusplus
}
#endif


#endif
