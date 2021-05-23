
BPAS_ADD_SANITY_TEST(RegularChain_test1 RegularChain_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/test.cpp
)

BPAS_ADD_SANITY_TEST(RegularChain_test2 RegularChain_test_sys_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/sysTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/sysGen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/sysGen2.cpp
)

BPAS_ADD_VALIDATE_TEST(RegularChain_test1 RegularChain_vtest_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/test.cpp
)

if (BPAS_WITH_MAPLE)

BPAS_ADD_VALIDATE_TEST(RegularChain_test2 RegularChain_vtest_sys_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/sysTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/sysGen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RegularChain/sysGen2.cpp
)

endif()
