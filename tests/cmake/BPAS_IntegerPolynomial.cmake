

BPAS_ADD_SANITY_TEST(duz_mul_test1 duz_mul_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/main64.cpp
	ARGUMENTS 512 512 1 0
)
BPAS_ADD_SANITY_TEST(duz_mul_test2 duz_mul_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/main64.cpp
	ARGUMENTS 512 512 3 0
)
BPAS_ADD_SANITY_TEST(duz_mul_test3 duz_mul_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/main64.cpp
	ARGUMENTS 512 512 4 0
)
BPAS_ADD_SANITY_TEST(duz_mul_test4 duz_mul_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/main64.cpp
	ARGUMENTS 512 512 5 0
)
BPAS_ADD_SANITY_TEST(duz_mul_test5 duz_mul_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/main64.cpp
	ARGUMENTS 512 51200 6 1
)


if (BPAS_WITH_NTL) 
BPAS_ADD_VALIDATE_TEST(duz_factor_test duz_factor_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/test-DUZP-Factoring.cpp
)
BPAS_ADD_SANITY_TEST(smz_factor_test smz_factor_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/IntegerPolynomial/test-SMZP-Factoring.cpp
)
endif()


