
if( (NOT BPAS_WITH_MAPLE) OR BPAS_BUILD_SERIAL ) 

BPAS_ADD_SANITY_TEST(TriangularSet_test TriangularSet_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/TriangularSet/test.cpp
)


BPAS_ADD_VALIDATE_TEST(TriangularSet_test TriangularSet_vtest_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/TriangularSet/test.cpp
)

endif()