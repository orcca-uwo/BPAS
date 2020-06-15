target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/PowerSeries.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UnivariatePolynomialOverPowerSeries.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Weierstrass.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/PowerSeries/UPOPS_Hensel.cpp
)
