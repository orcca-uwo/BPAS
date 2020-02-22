target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/complexdouble.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/rationalfunction_euclideanmethods.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/rationalfunction_symbolicintegration.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/multiprecision_rootfinding.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/rationalfunction_integrationpostprocessing.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/rationalfunction_complexrationalnumberordering.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/rationalfunction_symbolicnumericintegration.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RationalFunction/rationalfunction_integrationprinting.cpp
)
