
BPAS_ADD_SANITY_TEST(PowerSeries_test_ps PowerSeries_test_exe FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeries_TestCases.c
	${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/PowerSeries_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeries_test_upops PowerSeries_upops_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Upops_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Upops_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeries_test_weierstrass PowerSeries_weierstrass_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Weierstrass_TestCases.c
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Weierstrass_Test.c
)

BPAS_ADD_SANITY_TEST(PowerSeries_test_hensel PowerSeries_hensel_test_exe FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Hensel_TestCases.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSeries/Hensel_Test.cpp
)



