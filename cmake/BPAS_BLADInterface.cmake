

if(BPAS_WITH_BLAD)
	target_sources(${BPAS_LIB_TARGET} PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src/BLADInterface/bladinterface.c
	)
endif()