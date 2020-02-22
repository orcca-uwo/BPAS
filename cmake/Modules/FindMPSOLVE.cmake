# Try to find the MPSolve library
# https://github.com/robol/MPSolve
#
#
# Once done this will define
#
#  MPSOLVE_FOUND - system has MPSOLVE lib with correct version
#  MPSOLVE_INCLUDES - the MPSOLVE include directory
#  MPSOLVE_LIBRARIES - the MPSOLVE library
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.

set(MPSOLVE_POSSIBLE_INCLUDE_PATHS
  /usr/include
  /usr/include/mps
  /usr/local/include
  /usr/local/include/mps
)

set(MPSOLVE_POSSIBLE_LIB_PATHS
  /usr/lib
  /usr/local/lib
)

find_path(MPSOLVE_INCLUDES NAMES mps.h PATHS ${INCLUDE_INSTALL_DIR} HINTS ${MPSOLVE_POSSIBLE_INCLUDE_PATHS})
find_path(MPSOLVE_INCLUDES NAMES mps.c PATHS ${INCLUDE_INSTALL_DIR} HINTS ${MPSOLVE_POSSIBLE_INCLUDE_PATHS})

if(NOT MPSOLVE_FIND_VERSION)
  if(NOT MPSOLVE_FIND_VERSION_MAJOR)
    set(MPSOLVE_FIND_VERSION_MAJOR 3)
  endif()
  if(NOT MPSOLVE_FIND_VERSION_MINOR)
    set(MPSOLVE_FIND_VERSION_MINOR 0)
  endif()
  if(NOT MPSOLVE_FIND_VERSION_PATCH)
    set(MPSOLVE_FIND_VERSION_PATCH 0)
  endif()
  set(MPSOLVE_FIND_VERSION
    "${MPSOLVE_FIND_VERSION_MAJOR}.${MPSOLVE_FIND_VERSION_MINOR}.${MPSOLVE_FIND_VERSION_PATCH}")
endif()

if(MPSOLVE_INCLUDES)
  file(GLOB MPSOLVE_HEADERS "${MPSOLVE_INCLUDES}/mps.h")
  foreach(mps_header_filename ${MPSOLVE_HEADERS})
    file(READ "${mps_header_filename}" _mps_version_header)
    string(REGEX MATCH
      "This file is part of MPSolve ([0-9]+\\.[0-9]+\\.?[0-9]*)" _mps_version_match
      "${_mps_version_header}")
    if(_mps_version_match)
      set(MPSOLVE_VERSION "${CMAKE_MATCH_1}")
      #message("MPSOLVE_VERSION=${MPSOLVE_VERSION}")
      break()
    endif()
    string(REGEX MATCH
      "\\*\\*.*Version ([0-9]+\\.[0-9]+)" _mps_version_match
      "${_mps_version_header}")
    if(_mps_version_match)
      set(MPSOLVE_VERSION "${CMAKE_MATCH_1}")
      #message("MPSOLVE_VERSION=${MPSOLVE_VERSION}")
      break()
    endif()
  endforeach()

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT MPSOLVE_VERSION)
    set(MPSOLVE_VERSION_OK FALSE)
    message(STATUS "MPSOLVE version was not detected")
  elseif(${MPSOLVE_VERSION} VERSION_LESS ${MPSOLVE_FIND_VERSION})
    set(MPSOLVE_VERSION_OK FALSE)
    message(STATUS "MPSOLVE version ${MPSOLVE_VERSION} found in ${MPSOLVE_INCLUDES}, "
                   "but at least version ${MPSOLVE_FIND_VERSION} is required")
  else()
    set(MPSOLVE_VERSION_OK TRUE)
  endif()
endif()


find_library(MPSOLVE_LIBRARIES mps PATHS ${LIB_INSTALL_DIR} HINTS ${MPSOLVE_POSSIBLE_LIB_PATHS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPSOLVE DEFAULT_MSG
                                  MPSOLVE_INCLUDES MPSOLVE_LIBRARIES MPSOLVE_VERSION_OK)

#mark_as_advanced(MPSOLVE_INCLUDES MPSOLVE_LIBRARIES)
