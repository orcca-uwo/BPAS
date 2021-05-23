target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Hensel.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/DUZP_Hensel.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/DUZP_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/DBZP_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Support_Unpacked.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Support_Test.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Support_Recursive.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Support_Recursive_Unpacked.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Factoring_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/mzpolynomial.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/uzpolynomial.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Poly.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/BivariatePoly.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/Mul.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/MulNaive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/MulKS.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/MulDnC.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/MulToom4.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/MulToom8.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/Multiplication/MulSSA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/modpoly.cpp
)

if( BPAS_WITH_NTL ) 
target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/DUZP_Support_Factoring.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/NTL_mpz_DUZP_conversion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/IntegerPolynomial/SMZP_Support_Factoring.cpp
)
endif()
