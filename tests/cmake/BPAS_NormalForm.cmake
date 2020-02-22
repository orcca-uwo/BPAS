
if ((NOT BPAS_WITH_MAPLE) OR BPAS_BUILD_SERIAL)

BPAS_ADD_SANITY_TEST(NormalForm_test NormalForm_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/NormalForm/test.cpp
)

BPAS_ADD_VALIDATE_TEST(NormalForm_test NormalForm_vtest_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/NormalForm/test.cpp
)

endif()
