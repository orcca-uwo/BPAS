


BPAS_ADD_SANITY_TEST(SubResultantChain_test SubResultantChain_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/SubResultantChain/test.cpp
)

BPAS_ADD_VALIDATE_TEST(SubResultantChain_test SubResultantChain_vtest_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/SubResultantChain/test.cpp
)