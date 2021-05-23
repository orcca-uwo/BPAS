target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/PowerSeries.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/PowerSeriesZ.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UnivariatePolynomialOverPowerSeries.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Z.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Weierstrass.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Z_Weierstrass.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Hensel.cpp
)

if(NOT BPAS_BUILD_SERIAL) 

target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/PowerSeries_Parallel.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Parallel.cpp
)

endif()
