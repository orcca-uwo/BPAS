
#include <stdlib.h>
#include "Weierstrass_TestCases.h"

int main(int argc, char** argv) {
	int testDeg = 50;
	if (argc > 1 && atoi(argv[1])) {
		testDeg = atoi(argv[1]);
	}


	test1(testDeg);
	// test2(testDeg); //tests for invalid input on purpose and fails
	// test3(testDeg); //tests for invalid input on purpose and fails
	test4(testDeg);
	test5(testDeg);
	test6(testDeg);
	test7(testDeg);
	test8(testDeg);
	test9(testDeg);
	test10(testDeg);
	test11(testDeg);
	test12(testDeg);
	test13(testDeg);
	test14(testDeg);
	test15(testDeg);
	test16(testDeg);
	test17(testDeg);
	test18(testDeg);
	test19(testDeg);
	test20(testDeg);
	test21(testDeg);
	test22(testDeg);
	return 0;

}
