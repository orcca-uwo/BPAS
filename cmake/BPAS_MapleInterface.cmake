

# ----------------------------------------
# locate Maple libraries
set( MAPLEC_LIBRARIES "" CACHE STRING "Libraries for MAPLEC, to manually override search" )
set( MAPLEC_INCLUDES "" CACHE STRING "Include directory for MAPLEC, to manually override search" )
if ( "${MAPLEC_LIBRARIES}" STREQUAL "" OR "${MAPLEC_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for MPSolve. To override, set MAPLEC_LIBRARIES and MAPLEC_INCLUDES using ccmake." )
    unset( MAPLEC_LIBRARIES CACHE )
    unset( MAPLEC_INCLUDES CACHE )
    find_package( MAPLEC REQUIRED )
    if (${MAPLEC_FOUND})
		set(MAPLE_HOME "${MAPLEC_INCLUDES}/../.." CACHE STRING "Path to installation of Maple on this machine." FORCE)
	endif()
else()
    message( STATUS "User set MAPLEC_LIBRARIES. To change, edit MAPLEC_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${MAPLEC_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
	if ("${MAPLEC_INCLUDES}" STREQUAL "")
		set(MAPLE_HOME "" CACHE STRING "Path to installation of Maple on this machine.")
	else()
		if ("${MAPLEC_INCLUDES}" STREQUAL "MAPLEC_INCLUDES-NOTFOUND")
			set(MAPLE_HOME "" CACHE STRING "Path to installation of Maple on this machine.")
		else ()
			set(MAPLE_HOME "${MAPLEC_INCLUDES}/../.." CACHE STRING "Path to installation of Maple on this machine." FORCE)
		endif()
	endif()
endif()
message( STATUS "    MAPLEC_LIBRARIES:      ${MAPLEC_LIBRARIES}"      )
message( STATUS "    MAPLEC_INCLUDES:      ${MAPLEC_INCLUDES}"      )
set( EXTERNAL_LIBS ${MAPLEC_LIBRARIES} ${EXTERNAL_LIBS} )
include_directories( ${MAPLEC_INCLUDES} )


if("${MAPLE_HOME}" STREQUAL "") 
	message(WARNING "MAPLE_HOME is empty. Cannot run validate-tests.")
	add_custom_command(OUTPUT "NEED_MAPLE_WARNING.txt" COMMAND echo "\\nERROR: Cannot run validate-tests! MAPLE_HOME is not set.\\n" COMMAND exit 1 VERBATIM)
	add_custom_target(BPAS_NEED_MAPLE_WARNING DEPENDS "NEED_MAPLE_WARNING.txt")
	add_dependencies(${BPAS_LIB_TARGET} BPAS_NEED_MAPLE_WARNING)
endif()

target_sources(${BPAS_LIB_TARGET} PRIVATE
	${CMAKE_SOURCE_DIR}/src/MapleInterface/MapleInterface.c
	${CMAKE_SOURCE_DIR}/src/MapleInterface/MapleInterface.cpp
	${CMAKE_SOURCE_DIR}/src/MapleInterface/MapleWorkerThread.cpp
	${CMAKE_SOURCE_DIR}/src/MapleInterface/MapleInterfaceStream.cpp
	${CMAKE_SOURCE_DIR}/src/MapleInterface/MapleInterfacePipe.cpp
	${CMAKE_SOURCE_DIR}/src/MapleInterface/MaplePipePool.cpp
)
