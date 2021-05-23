
BPAS_ADD_SANITY_TEST(PowerSeries_test_ps PowerSeries_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeries_TestCases.c
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeries_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeriesZ_test_ps PowerSeriesZ_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeriesZ_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeriesZ_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeries_test_upops PowerSeries_upops_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Upops_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Upops_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeriesZ_test_upops PowerSeriesZ_upops_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/UpopsZ_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/UpopsZ_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeries_test_weierstrass PowerSeries_weierstrass_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Weierstrass_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Weierstrass_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeriesZ_test_weierstrass PowerSeriesZ_weierstrass_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/WeierstrassZ_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/WeierstrassZ_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeries_test_hensel PowerSeries_hensel_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Hensel_TestCases.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Hensel_Test.cpp
)

BPAS_ADD_SANITY_TEST(PowerSeries_parallel_test PowerSeries_parallel_test_exe FILES
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeries_Parallel_Test.cpp
)

BPAS_ADD_SANITY_TEST(PowerSeries_parallel_test_weierstrass PowerSeries_weierstrass_parallel_test_exe FILES
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Weierstrass_Parallel_Test.cpp
	ARGUMENTS 50 2 0 6
)

BPAS_ADD_SANITY_TEST(PowerSeries_parallel_test_hensel PowerSeries_hensel_parallel_test_exe FILES
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Hensel_Parallel_Test.cpp
)
