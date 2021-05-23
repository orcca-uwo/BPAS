
BPAS_ADD_SANITY_TEST(ModularPolynomial_test ModularPolynomial_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/ModularPolynomial/test.cpp
)

BPAS_ADD_SANITY_TEST(DUSP_test DUSP_test_exe FILES
	${CMAKE_CURRENT_SOURCE_DIR}/ModularPolynomial/dusp_test/test.cpp
)

BPAS_ADD_SANITY_TEST(DUSP_mul_time DUSP_mul_time_exe FILES
	${CMAKE_CURRENT_SOURCE_DIR}/ModularPolynomial/dusp_test/dusp_mul.cpp
)

BPAS_ADD_SANITY_TEST(DUSP_division_time DUSP_division_time_exe FILES
	${CMAKE_CURRENT_SOURCE_DIR}/ModularPolynomial/dusp_test/dusp_division.cpp
)
