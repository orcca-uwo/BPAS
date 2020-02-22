
cmake_minimum_required( VERSION 3.1 FATAL_ERROR)

project( BPAS_TESTS C CXX )


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

set(CTEST_OUTPUT_ON_FAILURE ON)

# ---------------------------------------
# Magic to find shared libraries for tests
# ---------------------------------------
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")



# ---------------------------------------
# TEST TARGETS
# ---------------------------------------

# --------------------------------------
# Declare target to build tests
set( BPAS_BUILD_TESTS_TARGET "check" )
add_custom_target(${BPAS_BUILD_TESTS_TARGET} DEPENDS ${BPAS_LIB_TARGET} )

# --------------------------------------
# Declare the sanity test target
set( BPAS_SANITY_TEST_TARGET "sanity-tests")
set( BPAS_SANITY_TESTS_MARK "BPAS_stest_")
set( BPAS_SANITY_TEST_ARGS )
set( BPAS_SANITY_TEST_LIBS ${BPAS_LIB_TARGET} ${EXTERNAL_LIBS})
add_custom_target(${BPAS_SANITY_TEST_TARGET} DEPENDS ${BPAS_LIB_TARGET} COMMAND ${CMAKE_CTEST_COMMAND} -R ${BPAS_SANITY_TESTS_MARK} --output-on-failure)

# --------------------------------------
# Declare the validate-tests target
set(BPAS_VALIDATE_TEST_TARGET "validate-tests")
set( BPAS_VALIDATE_TESTS_MARK "BPAS_vtest_")
set( BPAS_VALIDATE_TEST_LIBS ${BPAS_LIB_TARGET} ${EXTERNAL_LIBS})
set( BPAS_VALIDATE_TEST_ARGS "-DMAPLE_VALIDATE=1")
add_custom_target(${BPAS_VALIDATE_TEST_TARGET} DEPENDS ${BPAS_LIB_TARGET} COMMAND ${CMAKE_CTEST_COMMAND} -R ${BPAS_VALIDATE_TESTS_MARK} --output-on-failure)

# --------------------------------------
# Declare the verbose sanity test target
set( BPAS_V_SANITY_TEST_TARGET "sanity-tests-verbose")
add_custom_target(${BPAS_V_SANITY_TEST_TARGET} DEPENDS ${BPAS_LIB_TARGET} COMMAND ${CMAKE_CTEST_COMMAND} -R ${BPAS_SANITY_TESTS_MARK} --verbose)

# --------------------------------------
# Declare the verbose validate test target
set( BPAS_V_SANITY_TEST_TARGET "validate-tests-verbose")
add_custom_target(${BPAS_V_VALIDATE_TEST_TARGET} DEPENDS ${BPAS_LIB_TARGET} COMMAND ${CMAKE_CTEST_COMMAND} -R ${BPAS_VALIDATE_TESTS_MARK} --verbose)

# -----------------------------------------
# include dir
include_directories( ${CMAKE_SOURCE_DIR}/include )




# ----------------------------------------
# Give warnings based on specific compilation configurations.
# ----------------------------------------
if (BPAS_WITH_MAPLE AND !BPAS_BUILD_SERIAL) 
	message(WARNING "Cannot compile some tests if BUILD_WITH_MAPLE is ON and BPAS_BUILD_SERIAL is OFF.")
endif()



# ----------------------------------------
# BPAS Test Macros
# ----------------------------------------
include (${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_TestMacros.cmake )



# ----------------------------------------
# BPAS MapleTestTool
# ----------------------------------------
include (${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_MapleTestTool.cmake )



# ----------------------------------------
# BPAS SOURCES & pre-build scripts
# ----------------------------------------

# --------------------------------------
# DataStructures
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_DataStructures.cmake/ )

# --------------------------------------
# ExpressionTree
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_ExpressionTree.cmake/ )

# --------------------------------------
# FFT
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_FFT.cmake/ )

# --------------------------------------
# ModularPolynomial
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_ModularPolynomial.cmake/ )

# --------------------------------------
# IntegerPolynomial
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_IntegerPolynomial.cmake/ )

# --------------------------------------
# RealRootIsolation
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_RealRootIsolation.cmake/ )

# --------------------------------------
# SUPOperations
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_SUPOperations.cmake/ )

# --------------------------------------
# SMPOperations
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_SMPOperations.cmake/ )

# --------------------------------------
# SquareFree
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_SquareFree.cmake/ )

# --------------------------------------
# Ring
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_Ring.cmake/ )

# --------------------------------------
# SubResultantChain
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_SubResultantChain.cmake/ )

# --------------------------------------
# TriangularSet
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_TriangularSet.cmake/ )

# --------------------------------------
# NormalForm
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_NormalForm.cmake/ )

# --------------------------------------
# RegularChain
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_RegularChain.cmake/ )

# --------------------------------------
# Parser
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_Parser.cmake/ )

# --------------------------------------
# RationalFunction
include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_RationalFunction.cmake/ )

# --------------------------------------
# PolyhedralSets
option(BPAS_WITH_POLYHEDRA "Build with PolyhedralSets" OFF)
if(BPAS_WITH_POLYHEDRA)
    include( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BPAS_PolyhedralSets.cmake/ )
endif()





