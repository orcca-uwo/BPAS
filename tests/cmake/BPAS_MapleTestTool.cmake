

set (MAPLE_SANITY_LIB mapleTestToolSanity)
add_library(${MAPLE_SANITY_LIB} "")
target_sources(${MAPLE_SANITY_LIB} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/MapleTestTool/MapleTestTool.cpp
)
target_link_libraries(${MAPLE_SANITY_LIB} ${BPAS_SANITY_TEST_LIBS})
add_dependencies(${BPAS_SANITY_TEST_TARGET} ${MAPLE_SANITY_LIB})
set( BPAS_SANITY_TEST_LIBS ${BPAS_SANITY_TEST_LIBS} ${MAPLE_SANITY_LIB} )

set (MAPLE_VALIDATE_LIB mapleTestToolValidate)
add_library(${MAPLE_VALIDATE_LIB} "")
target_sources(${MAPLE_VALIDATE_LIB} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/MapleTestTool/MapleTestTool.cpp
)
target_compile_options(${MAPLE_VALIDATE_LIB} PUBLIC 
	${BPAS_VALIDATE_TEST_ARGS}
)
target_link_libraries(${MAPLE_VALIDATE_LIB} ${BPAS_VALIDATE_TEST_LIBS})
add_dependencies(${BPAS_VALIDATE_TEST_TARGET} ${MAPLE_VALIDATE_LIB})
set( BPAS_VALIDATE_TEST_LIBS ${BPAS_VALIDATE_TEST_LIBS} ${MAPLE_VALIDATE_LIB} )


# ----------------------------------------
# locate Maple libraries
set( MAPLEC_LIBRARIES "" CACHE STRING "Libraries for MAPLEC, to manually override search" )
set( MAPLEC_INCLUDES "" CACHE STRING "Include directory for MAPLEC, to manually override search" )
if ( "${MAPLEC_LIBRARIES}" STREQUAL "" OR "${MAPLEC_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for MPSolve. To override, set MAPLEC_LIBRARIES and MAPLEC_INCLUDES using ccmake." )
    unset( MAPLEC_LIBRARIES CACHE )
    unset( MAPLEC_INCLUDES CACHE )
    find_package( MAPLEC )
    if (${MAPLEC_FOUND})
		set(MAPLE_HOME "${MAPLEC_INCLUDES}/../.." CACHE STRING "Path to installation of Maple on this machine." FORCE)
	endif()
else()
    message( STATUS "User set MAPLEC_LIBRARIES. To change, edit MAPLEC_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${MAPLEC_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "Warning: file ${LIB} does not exist.\n" )
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


if("${MAPLE_HOME}" STREQUAL "") 
	message(WARNING "MAPLE_HOME is empty. Cannot run validate-tests.")
	add_custom_command(OUTPUT "NEED_MAPLE_WARNING.txt" COMMAND echo "\\nERROR: Cannot run validate-tests! MAPLE_HOME is not set.\\n" COMMAND exit 1 VERBATIM)
	add_custom_target(BPAS_NEED_MAPLE_WARNING DEPENDS "NEED_MAPLE_WARNING.txt")
	add_dependencies(${BPAS_VALIDATE_TEST_TARGET} BPAS_NEED_MAPLE_WARNING)
else()
	include_directories( ${MAPLEC_INCLUDES} )
	set( BPAS_VALIDATE_TEST_LIBS ${MAPLEC_LIBRARIES} ${BPAS_VALIDATE_TEST_LIBS} )
endif()
