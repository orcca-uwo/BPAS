
BPAS_ADD_SANITY_TEST(PolyhedralSets_test1 PolyhedralSets_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralSets/test.c
	ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralSets/testCase/t1.txt 5 10
)

BPAS_ADD_SANITY_TEST(PolyhedralSets_test2 PolyhedralSets_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralSets/test.c
	ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralSets/testCase/t3.txt 4 8
)