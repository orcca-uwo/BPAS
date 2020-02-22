target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/ModularPolynomial/sdmpolynomial.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/ModularPolynomial/DUSP_Support_Test.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/ModularPolynomial/DUSP_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/ModularPolynomial/DUSP_FFT_Support.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/ModularPolynomial/DUSP_Support_Factoring.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/ModularPolynomial/DUSP_NTL_Support.cpp
)
