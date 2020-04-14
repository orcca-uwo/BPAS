cmake_minimum_required( VERSION 3.1 FATAL_ERROR)

project( BPAS C CXX )


# ---------------------------------------
# SETUP
# ---------------------------------------

# Setup cmake paths, variables, etc.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#Verbose makefiles
set(CMAKE_VERBOSE_MAKEFILE off)
#don't regenerate files during make
set(CMAKE_SUPPRESS_REGENERATION on)

set(CMAKE_INCLUDE_PATH /usr/include;/usr/include/x86_64-linux-gnu CACHE FILEPATH "cmake include path")

set(CMAKE_INSTALL_PREFIX /usr/local CACHE STRING "cmake install prefix")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# ---------------------------------------
# OUT OF SOURCE BUILD
# ---------------------------------------
# force an out-of-source build, to not overwrite the existing Makefiles
# (out-of-source is cleaner, too)
string( COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" BPAS_COMPILE_INPLACE )
if(  BPAS_COMPILE_INPLACE )
	message( FATAL_ERROR "Compiling BPAS should be done in an out-of-source directory. To proceed:
    rm -rf CMakeCache.txt CMakeFiles/   # delete files in ${CMAKE_SOURCE_DIR}
    mkdir build
    cd build
    cmake ..
    make" )
endif()



#########################################
# FAIL IF NOT LINUX X86... temporary
#########################################
option(BPAS_BUILD_LINUXINTEL64 "Build BPAS for an intel x86_64 chip on Linux" ON)
if (NOT BPAS_BUILD_LINUXINTEL64)
    message(FATAL_ERROR "BPAS can only be built on intel x86_64 machines running LINUX... sorry.")
endif()



# --------------------------------------
# Declare the library
set(BPAS_LIB_TARGET "bpas")
add_library( ${BPAS_LIB_TARGET} "" )

# -----------------------------------------
# include dir
include_directories( ${CMAKE_SOURCE_DIR}/include )



# ---------------------------------------
# BUILD TYPES, DEFAULT FLAGS.
# ---------------------------------------

# Set a default build type if none was specified
set(default_build_type "Release")
# if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
#   set(default_build_type "Debug")
# endif()
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo")
    # "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()



set(CMAKE_C_FLAGS "")
set(CMAKE_C_FLAGS_RELEASE "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "")
set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
set(CMAKE_CXX_FLAGS_DEBUG "")

# ----------------------------------------
# Ensure 6 <= GCC < 8 for cilk
set(GCC_MIN_VERSION 6.0)
set(GCC_MAX_VERSION 8.0)

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_MIN_VERSION OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER GCC_MAX_VERSION OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL GCC_MAX_VERSION)
    message(FATAL_ERROR "ERROR: BPAS requires GCC version at least ${GCC_MIN_VERSION} and less than ${GCC_MAX_VERSION} to build but found ${CMAKE_CXX_COMPILER_VERSION}")
endif()


# ----------------------------------------
# use C++11 and C99
# see http://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message( WARNING "The compiler ${CMAKE_CXX_COMPILER} doesn't support the -std=c++11 flag. Some code may not compile.")
endif()

CHECK_C_COMPILER_FLAG("-std=gnu11" COMPILER_SUPPORTS_GNU11)
if(COMPILER_SUPPORTS_GNU11)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11")
else()
    message( WARNING "The compiler ${CMAKE_C_COMPILER} doesn't support the -std=gnu11 flag. Some code may not compile.")
endif()

CHECK_C_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
CHECK_C_COMPILER_FLAG("-march=x86-64" COMPILER_SUPPORTS_MARCH_X86)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    option(BPAS_OPTIMIZE_FOR_NATIVE "Build with -march=native; overrides OPTIMIZE_FOR_X86" ON)
endif()

if(COMPILER_SUPPORTS_MARCH_X86)
    if (COMPILER_SUPPORTS_MARCH_NATIVE) 
        option(BPAS_OPTIMIZE_FOR_X86 "Build with -march=x86-64" OFF)
		if (BPAS_OPTIMIZE_FOR_NATIVE)
   			set(BPAS_OPTIMIZE_FOR_X86 OFF CACHE BOOL "Build with -march=x86-64" FORCE)
		endif()	
	else()
        option(BPAS_OPTIMIZE_FOR_X86 "Build with -march=x86-64" ON)
    endif()
endif()

if(NOT COMPILER_SUPPORTS_MARCH_X86 AND NOT COMPILER_SUPPORTS_MARCH_NATIVE)
    message( WARNING "The compiler ${CMAKE_C_COMPILER} doesn't support the -march flag. Some code may not compile properly.")
endif()


if(BPAS_OPTIMIZE_FOR_NATIVE) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
elseif(BPAS_OPTIMIZE_FOR_X86)
    set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64")
endif()

# Note: Compiler flags are called in this order: "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_{BUILD_TYPE}}" 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic -fcilkplus -DLINUXINTEL64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -fcilkplus -DLINUXINTEL64")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2 -fomit-frame-pointer -DNDEBUG=1")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g -Wall -Wno-unused-function -Werror='return-type'")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -pg -Wall -Wno-unused-function -Werror='return-type'")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -fomit-frame-pointer -DNDEBUG=1")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g -Wall -Wno-unused-function -Wno-unused-variable -Werror='return-type'")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -pg -Wall -Wno-unused-function -Wno-unused-variable -Werror='return-type' -Werror='discarded-qualifiers' -Werror='implicit'")




# ---------------------------------------
# STANDARD LIBRARIES
# ---------------------------------------
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} "m" "pthread")

# ---------------------------------------
# LOCATE LIBRARIES
# ---------------------------------------

# ----------------------------------------
# locate GMP libraries
set( GMP_LIBRARIES "" CACHE STRING "Libraries for GMP, to manually override search" )
set( GMP_INCLUDES "" CACHE STRING "Include directory for GMP, to manually override search" )
if ( "${GMP_LIBRARIES}" STREQUAL "" OR "${GMP_INCLUDES}" STREQUAL "")
	message( STATUS "Searching for GMP. To override, set GMP_LIBRARIES and GMP_INCLUDES using ccmake." )
	unset( GMP_LIBRARIES CACHE )
	unset( GMP_INCLUDES CACHE )
	find_package( GMP REQUIRED )
else()
	message( STATUS "User set GMP_LIBRARIES. To change, edit GMP_LIBRARIES using ccmake (set to empty to enable search)." )
        # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
	foreach( LIB ${GMP_LIBRARIES} )
                if ( NOT EXISTS ${LIB} )
                        message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
                endif()
        endforeach()
endif()
message( STATUS "    GMP_LIBRARIES:      ${GMP_LIBRARIES}"      )
message( STATUS "    GMP_INCLUDES :      ${GMP_INCLUDES}"       )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${GMP_LIBRARIES} )
include_directories( ${GMP_INCLUDES} )

# -----------------------------------------
# locate GMPXX libraries
set( GMPXX_LIBRARIES "" CACHE STRING "Libraries for GMPXX, to manually override search" )
set( GMPXX_INCLUDES "" CACHE STRING "Include directory for GMPXX, to manually override search" )
if ( "${GMPXX_LIBRARIES}" STREQUAL "" OR "${GMPXX_INCLUDES}" STREQUAL "")
        message( STATUS "Searching for GMPXX. To override, set GMPXX_LIBRARIES and GMPXX_INCLUDES using ccmake." )
        unset( GMPXX_LIBRARIES CACHE )
        unset( GMPXX_INCLUDES CACHE )
        find_package( GMPXX REQUIRED )
else()
        message( STATUS "User set GMPXX_LIBRARIES. To change, edit GMPXX_LIBRARIES using ccmake (set to empty to enable search)." )
        # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
        foreach( LIB ${GMPXX_LIBRARIES} )
                if ( NOT EXISTS ${LIB} )
                        message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
                endif()
        endforeach()
endif()
message( STATUS "    GMPXX_LIBRARIES:      ${GMPXX_LIBRARIES}"      )
message( STATUS "    GMPXX_INCLUDES :      ${GMPXX_INCLUDES}"       )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${GMPXX_LIBRARIES} )
include_directories( ${GMPXX_INCLUDES} )

# ----------------------------------------
# locate MPSolve libraries
set( MPSOLVE_LIBRARIES "" CACHE STRING "Libraries for MPSOLVE, to manually override search" )
set( MPSOLVE_INCLUDES "" CACHE STRING "Include directory for MPSOLVE, to manually override search" )
if ( "${MPSOLVE_LIBRARIES}" STREQUAL "" OR "${MPSOLVE_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for MPSolve. To override, set MPSOLVE_LIBRARIES and MPSOLVE_INCLUDES using ccmake." )
    unset( MPSOLVE_LIBRARIES CACHE )
    unset( MPSOLVE_INCLUDES CACHE )
    find_package( MPSOLVE 2.2 REQUIRED )
    #So actually we need version 3.0+, but the mps.h wasn't updated to say version 3.0 
    #until 3.1.5....
else()
    message( STATUS "User set MPSOLVE_LIBRARIES. To change, edit MPSOLVE_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${MPSOLVE_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
endif()
message( STATUS "    MPSOLVE_LIBRARIES:      ${MPSOLVE_LIBRARIES}"      )
message( STATUS "    MPSOLVE_INCLUDES:      ${MPSOLVE_INCLUDES}"      )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${MPSOLVE_LIBRARIES} )
include_directories( ${MPSOLVE_INCLUDES} )


# ----------------------------------------
# modpn
set( MODPN_LIBRARIES "${CMAKE_SOURCE_DIR}/src/libmodpnLINUXINTEL64.so" CACHE STRING "Libraries for modpn. Included with BPAS distribution.")
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${MODPN_LIBRARIES})
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libmodpnLINUXINTEL64.so COMMAND cp ${CMAKE_SOURCE_DIR}/src/libmodpnLINUXINTEL64.so ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Copying libmodpn.so to current build directory...")
add_custom_target(COPY_MODPN_LIB DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libmodpnLINUXINTEL64.so)
add_dependencies(${BPAS_LIB_TARGET} COPY_MODPN_LIB)







# ----------------------------------------
# BPAS Options and Custom Defines
# ----------------------------------------

option(BPAS_WITH_BLAD "Build with BLAD" OFF)
if(BPAS_WITH_BLAD)
# ----------------------------------------
# locate BLAD libraries
set( BLAD_LIBRARIES "" CACHE STRING "Libraries for BLAD, to manually override search" )
set( BLAD_INCLUDES "" CACHE STRING "Include directory for BLAD, to manually override search" )
if ( "${BLAD_LIBRARIES}" STREQUAL "" OR "${BLAD_INCLUDES}" STREQUAL "")
    message( STATUS "Searching for BLAD. To override, set BLAD_LIBRARIES and BLAD_INCLUDES using ccmake." )
    unset( BLAD_LIBRARIES CACHE )
    unset( BLAD_INCLUDES CACHE )
    find_package( BLAD 3.0 REQUIRED )
else()
    message( STATUS "User set BLAD_LIBRARIES. To change, edit BLAD_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${BLAD_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
endif()
message( STATUS "    BLAD_LIBRARIES:      ${BLAD_LIBRARIES}"      )
message( STATUS "    BLAD_INCLUDES:      ${BLAD_INCLUDES}"      )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${BLAD_LIBRARIES} )
include_directories( ${BLAD_INCLUDES} )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWITH_BLAD=1 -DBA0_OLDDEF=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_BLAD=1")
endif()



option(BPAS_WITH_NTL "Build with NTL" OFF)
if(BPAS_WITH_NTL)
# ----------------------------------------
# locate NTL libraries
set( NTL_LIBRARIES "" CACHE STRING "Libraries for NTL, to manually override search" )
set( NTL_INCLUDES "" CACHE STRING "Include directory for NTL, to manually override search" )
if ( "${NTL_LIBRARIES}" STREQUAL "" OR "${NTL_LIBRARIES}" STREQUAL "NTL_LIBRARIES-NOTFOUND"
    OR "${NTL_INCLUDES}" STREQUAL "" OR "${NTL_INCLUDES}" STREQUAL "NTL_INCLUDES-NOTFOUND")
    message( STATUS "Searching for NTL. To override, set NTL_LIBRARIES and NTL_INCLUDES using ccmake." )
    unset( NTL_LIBRARIES CACHE )
    unset( NTL_INCLUDES CACHE )
    find_package( NTL 8.0 REQUIRED )
else()
    message( STATUS "User set NTL_LIBRARIES. To change, edit NTL_LIBRARIES using ccmake (set to empty to enable search)." )
    # Check existence -- but this may be okay, if the user entered, e.g., -lgmp instead of /path/to/gmp.a
    foreach( LIB ${NTL_LIBRARIES} )
        if ( NOT EXISTS ${LIB} )
            message( WARNING "\n      Warning: file ${LIB} does not exist.\n" )
        endif()
    endforeach()
endif()
message( STATUS "    NTL_LIBRARIES:      ${NTL_LIBRARIES}"      )
message( STATUS "    NTL_INCLUDES:      ${NTL_INCLUDES}"      )
set( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${NTL_LIBRARIES} )
include_directories( ${NTL_INCLUDES} )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWITH_NTL=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_NTL=1")
endif()

option(BPAS_WITH_MAPLE "Build with Maple" ON)
if(BPAS_WITH_MAPLE)
    include(${CMAKE_SOURCE_DIR}/cmake/BPAS_MapleInterface.cmake)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWITH_MAPLE=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MAPLE=1")
    message(WARNING "When linking to BPAS when built with Maple it is required to set your LD_LIBRARY_PATH and MAPLE environment variables to point to your Maple installation.")
endif()

option(BPAS_BUILD_SERIAL "Build without parallelism or cilk" OFF)
if(BPAS_BUILD_SERIAL)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSERIAL=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSERIAL=1")
endif()

option(BPAS_BUILD_PROFILING "Build with profiling" OFF)
if(BPAS_BUILD_PROFILING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPROFILING_ENABLED=1 -DREGULARCHAIN_PROFILING=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROFILING_ENABLED=1 -DREGULARCHAIN_PROFILING=1")
endif()

option(BPAS_BUILD_GMP_PROFILING "Build with profiling of gmp" OFF)
if(BPAS_BUILD_GMP_PROFILING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPROFILING_GMP_ENABLED=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROFILING_GMP_ENABLED=1")
endif()

option(BPAS_BUILD_RC_PROFILING "Build with regular chain profiling" OFF)
if(BPAS_BUILD_RC_PROFILING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DREGULARCHAIN_PROFILING=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DREGULARCHAIN_PROFILING=1")
endif()

# Needs cilkview. If this debug option is resurrected need to find cilkview
# option(BPAS_BUILD_TDEBUG "Build with more verbose debugging in FieldPolynomial" OFF)
# if(BPAS_BUILD_TDEBUG)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DTDEBUG=1")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTDEBUG=1")
# endif()

option(BPAS_BUILD_TS_DEBUG "Build with TriangularSet debugging" OFF)
if(BPAS_BUILD_TS_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DTRIANGULARSET_DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRIANGULARSET_DEBUG")
endif()

option(BPAS_BUILD_ZDRC_DEBUG "Build with ZeroDimensionalRegularChain debugging" OFF)
if(BPAS_BUILD_ZDRC_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DZERODIMENSIONALREGULARCHAIN_DEBUG=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DZERODIMENSIONALREGULARCHAIN_DEBUG=1")
endif()

option(BPAS_BUILD_RC_DEBUG "Build with RegularChain debugging" OFF)
if(BPAS_BUILD_RC_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DREGULARCHAIN_DEBUG=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DREGULARCHAIN_DEBUG=1")
endif()

option(BPAS_BUILD_RC_DEBUG_RRC "Build with RegularChain debugging of remove redundant components" OFF)
if(BPAS_BUILD_RC_DEBUG_RRC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DREGULARCHAIN_DEBUGII=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DREGULARCHAIN_DEBUGII=1")
endif()


option(BPAS_BUILD_RC_GENERATORS "Build with RegularChain parallel generators" OFF)
if(BPAS_BUILD_RC_GENERATORS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DRC_WITH_GENERATORS=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRC_WITH_GENERATORS=1")
endif()

option(BPAS_BUILD_FACTORING_DEBUG "Build with factoring debugging" OFF)
if (BPAS_BUILD_FACTORING_DEBUG) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSMZP_FACTORING_DEBUG=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSMZP_FACTORING_DEBUG=1")
endif()

#-----------------------------------------
# Finally set cache compile flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "C compiler flags" FORCE)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "Release C flags" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "Release with deb info C flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "Debug CXX flags" FORCE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "CXX compuiler flags" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "Release CXX flags" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "Release with deb info CXX flags" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "Debug CXX flags" FORCE)


# -----------------------------------------
# Link with externals
target_link_libraries( ${BPAS_LIB_TARGET} ${EXTERNAL_LIBS} )



# ----------------------------------------
# BPAS SOURCES & pre-build scripts
# ----------------------------------------

# --------------------------------------
# Global 
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Src.cmake/ )

# --------------------------------------
# BLADInterface 
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_BLADInterface.cmake/ )

# --------------------------------------
# FFT 
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_FFT.cmake/ )

# --------------------------------------
# Ring
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Ring.cmake/ )

# --------------------------------------
# DataStructures
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_DataStructures.cmake/ )

# --------------------------------------
# FiniteFields
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_FiniteFields.cmake/ )

# --------------------------------------
# Symbol
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Symbol.cmake/ )

# --------------------------------------
# ExpressionTree
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_ExpressionTree.cmake/ )

# --------------------------------------
# RingPolynomial
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_RingPolynomial.cmake/ )

# --------------------------------------
# ModularPolynomial
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_ModularPolynomial.cmake/ )

# --------------------------------------
# IntegerPolynomial
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_IntegerPolynomial.cmake/ )

# --------------------------------------
# DyadicRationalNumber
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_DyadicRationalNumber.cmake/ )

# --------------------------------------
# Interval
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Interval.cmake/ )

# --------------------------------------
# RationalNumberPolynomial
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_RationalNumberPolynomial.cmake/ )

# --------------------------------------
# Subresultant
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Subresultant.cmake/ )

# --------------------------------------
# SubResultantChain
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_SubResultantChain.cmake/ )

# --------------------------------------
# TriangularSet
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_TriangularSet.cmake/ )

# --------------------------------------
# RegularChain
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_RegularChain.cmake/ )

# --------------------------------------
# RealRootIsolation
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_RealRootIsolation.cmake/ )

# --------------------------------------
# RationalFunction
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_RationalFunction.cmake/ )

# --------------------------------------
# Parser
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Parser.cmake/ )

# --------------------------------------
# Utils
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_Utils.cmake/ )

# --------------------------------------
# LinearAlgebra
include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_LinearAlgebra.cmake/ )



# --------------------------------------
# PolyhedralSets
option(BPAS_WITH_POLYHEDRA "Build with PolyhedralSets" OFF)
if(BPAS_WITH_POLYHEDRA)
    include( ${CMAKE_SOURCE_DIR}/cmake/BPAS_PolyhedralSets.cmake/ )
endif()


#file(GLOB_RECURSE all_c_files ${CMAKE_SOURCE_DIR}/src/*/*.c) 
#set_source_files_properties(${all_c_files} PROPERTIES LANGUAGE CXX)


# --------------------------------------
# TESTING
# --------------------------------------

include(CTest)
add_subdirectory( ${CMAKE_SOURCE_DIR}/tests ${CMAKE_BINARY_DIR}/tests EXCLUDE_FROM_ALL)



# --------------------------------------
# DOC
# --------------------------------------

set (DOC_TARGET doc)
add_custom_target(${DOC_TARGET})
add_custom_command(
	TARGET ${DOC_TARGET}
	COMMAND rm -rf help/html help-simple/html
	COMMAND doxygen Doxyfile
	COMMAND doxygen Doxyfile-simple
	COMMAND tclsh doc/insertSimpleDiagram.tcl help-simple/html help/html
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_custom_command(
	TARGET ${DOC_TARGET}
	COMMAND make
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/latex
)
add_custom_command(
	TARGET ${DOC_TARGET}
	COMMAND tclsh insertClassTex.tcl latex ../help/html
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
)
add_custom_command(
	TARGET ${DOC_TARGET}
	COMMAND rm -rf help-simple
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
