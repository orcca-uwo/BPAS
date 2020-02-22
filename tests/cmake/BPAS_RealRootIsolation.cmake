
BPAS_ADD_SANITY_TEST(RealRootIsolation_test1 RealRootIsolation_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/RealRootIsolation/test.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RealRootIsolation/examples.cpp
	ARGUMENTS 1 1 4096 1
)

BPAS_ADD_SANITY_TEST(RealRootIsolation_test2 RealRootIsolation_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/RealRootIsolation/test.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RealRootIsolation/examples.cpp
	ARGUMENTS 1 2 4096 1
)