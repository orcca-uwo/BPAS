

BPAS_ADD_SANITY_TEST(Parser_test Parser_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/Parser/myparser.cpp
	ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/Parser/data.txt
)

BPAS_ADD_VALIDATE_TEST(Parser_test Parser_vtest_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/Parser/myparser.cpp
	ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/Parser/data.txt
)