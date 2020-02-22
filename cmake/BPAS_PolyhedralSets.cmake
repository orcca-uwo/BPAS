

# ----------------------------------------
# locate MPFR libraries
set( MPFR_LIBRARIES "" CACHE STRING "Libraries for MPFR, to manually override search" )
set( MPFR_INCLUDES "" CACHE STRING "Include directory for MPFR, to manually override search" )
if ( "${MPFR_LIBRARIES}" STREQUAL "" OR "${MPFR_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for MPSolve. To override, set MPFR_LIBRARIES and MPFR_INCLUDES using ccmake." )
    unset( MPFR_LIBRARIES CACHE )
    unset( MPFR_INCLUDES CACHE )
    find_package( MPFR REQUIRED )
else()
    message( STATUS "User set MPFR_LIBRARIES. To change, edit MPFR_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${MPFR_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
endif()
message( STATUS "    MPFR_LIBRARIES:      ${MPFR_LIBRARIES}"      )
message( STATUS "    MPFR_INCLUDES:      ${MPFR_INCLUDES}"      )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${MPFR_LIBRARIES} )
include_directories( ${MPFR_INCLUDES} )

# ----------------------------------------
# locate FLINT libraries
set( FLINT_LIBRARIES "" CACHE STRING "Libraries for FLINT, to manually override search" )
set( FLINT_INCLUDES "" CACHE STRING "Include directory for FLINT, to manually override search" )
if ( "${FLINT_LIBRARIES}" STREQUAL "" OR "${FLINT_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for MPSolve. To override, set FLINT_LIBRARIES and FLINT_INCLUDES using ccmake." )
    unset( FLINT_LIBRARIES CACHE )
    unset( FLINT_INCLUDES CACHE )
    find_package( FLINT REQUIRED )
else()
    message( STATUS "User set FLINT_LIBRARIES. To change, edit FLINT_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${FLINT_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
endif()
message( STATUS "    FLINT_LIBRARIES:      ${FLINT_LIBRARIES}"      )
message( STATUS "    FLINT_INCLUDES:      ${FLINT_INCLUDES}"      )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${FLINT_LIBRARIES} )
include_directories( ${FLINT_INCLUDES} )

# ----------------------------------------
# locate CDDGMP libraries
set( CDDGMP_LIBRARIES "" CACHE STRING "Libraries for CDDGMP, to manually override search" )
set( CDDGMP_INCLUDES "" CACHE STRING "Include directory for CDDGMP, to manually override search" )
if ( "${CDDGMP_LIBRARIES}" STREQUAL "" OR "${CDDGMP_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for MPSolve. To override, set CDDGMP_LIBRARIES and CDDGMP_INCLUDES using ccmake." )
    unset( CDDGMP_LIBRARIES CACHE )
    unset( CDDGMP_INCLUDES CACHE )
    find_package( CDDGMP REQUIRED )
else()
    message( STATUS "User set CDDGMP_LIBRARIES. To change, edit CDDGMP_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${CDDGMP_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
endif()
message( STATUS "    CDDGMP_LIBRARIES:      ${CDDGMP_LIBRARIES}"      )
message( STATUS "    CDDGMP_INCLUDES:      ${CDDGMP_INCLUDES}"      )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${CDDGMP_LIBRARIES} )
include_directories( ${CDDGMP_INCLUDES} )




set(BPAS_POLYHEDRAL_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_inequality.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_kohler.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_projection.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_fme.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_linAlg.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_balas.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_unrolledll.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_blockElimination.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyhedralSets/FME_Support_extreme.c
)
mark_as_advanced(BPAS_POLYHEDRAL_SRCS)

set_source_files_properties(
    ${BPAS_POLYHEDRAL_SRCS}
    PROPERTIES COMPILE_FLAGS "-lmpfr -lcddgmp -lgmp -lflint -DGMPRATIONAL=1"
)

target_sources(${BPAS_LIB_TARGET} PRIVATE
    ${BPAS_POLYHEDRAL_SRCS}
)


