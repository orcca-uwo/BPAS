
#include <stdlib.h>
#include "WeierstrassZ_TestCases.h"

int main(int argc, char** argv) {
	int testDeg = 50;
	if (argc > 1 && atoi(argv[1])) {
		testDeg = atoi(argv[1]);
	}

//	test22(10);



	// return 0;

	test1(testDeg);
	// test2(testDeg); //tests for invalid input on purpose and fails
	// test3(testDeg); //tests for invalid input on purpose and fails
	// test5(testDeg); //results in rational numbers
	test4(testDeg);
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
	test23(testDeg);
	test24(testDeg);


	return 0;

}
