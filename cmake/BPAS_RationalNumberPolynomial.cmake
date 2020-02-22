target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/urpolynomial_realroot.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/urpolynomial_taylorshift.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/mrpolynomial.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_Support_Unpacked.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_Support_Recursive_Unpacked.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_Support-AA.c
	# ${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_SparseInterpolation-AA.c
	# ${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_Interpolator-AA.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/urpolynomial.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SUQP_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_Support_Test-AA.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/SMQP_Support_Recursive-AA.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalNumberPolynomial/mrpolynomial_mdd-altarr.cpp
)
