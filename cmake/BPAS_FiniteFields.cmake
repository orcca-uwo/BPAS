target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/SmallPrimeField.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/BigPrimeField.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/GeneralizedFermatPrimeField.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/BigPrimeField_Support.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/SmallPrimeField_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/GFPF_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/cpu_timer.c
)


set_source_files_properties(
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/GFPF_Support.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/FiniteFields/cpu_timer.c
	PROPERTIES COMPILE_FLAGS "-Ofast"
)