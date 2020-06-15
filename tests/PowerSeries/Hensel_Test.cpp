
#include "Hensel_TestCases.hpp"

#include <stdlib.h>

int main(int argc, char** argv) {
	int degree =50;
	if (argc > 1) {
        degree = atoi(argv[1]);
    }


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
	test12(degree);
	test13(degree);
	test14(degree);
	test15();
	test16(degree);
	test18(degree);
	test19(degree);
	test20(degree);
	test21(degree);
	test22(degree);
	test23(degree);
	test24(degree);
	test25(degree);
	test26(degree);
	test27(degree);

#if defined(WITH_NTL) && WITH_NTL
	test28(degree);
	test29(degree);
	test30(degree);
	test31(degree);
	test32(degree);
	test33(degree);
	test34(degree);
	test35(degree);
	test36(degree);
	test37(degree);
	test38(degree);
	test39(degree);
	test40(degree);
	test41(degree);
	test42(degree);
#endif

	return 0;

}