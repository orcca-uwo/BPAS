# Try to find the NTL library
# https://www.shoup.net/ntl/
#
# Once done this will define
#
#  NTL_FOUND - system has NTL lib with correct version
#  NTL_INCLUDES - the NTL include directory
#  NTL_LIBRARIES - the NTL library
#  NTL_VERSION - NTL version
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.

set(NTL_POSSIBLE_INCLUDE_PATHS
  /usr/include
  /usr/include/NTL
  /usr/include/ntl
  /usr/local/include
  /usr/local/include/NTL
  /usr/local/include/ntl
)

set(NTL_POSSIBLE_LIB_PATHS
  /usr/lib
  /usr/local/lib
)

find_path(NTL_INCLUDES NAMES version.h PATHS ${INCLUDE_INSTALL_DIR} HINTS ${NTL_POSSIBLE_INCLUDE_PATHS})


# Set NTL_FIND_VERSION to 3.0.0 if no minimum version is specified
if(NOT NTL_FIND_VERSION)
  if(NOT NTL_FIND_VERSION_MAJOR)
    set(NTL_FIND_VERSION_MAJOR 8)
  endif()
  if(NOT NTL_FIND_VERSION_MINOR)
    set(NTL_FIND_VERSION_MINOR 8)
  endif()
  if(NOT NTL_FIND_VERSION_PATCH)
    set(NTL_FIND_VERSION_PATCH 0)
  endif()
  set(NTL_FIND_VERSION
    "${NTL_FIND_VERSION_MAJOR}.${NTL_FIND_VERSION_MINOR}.${NTL_FIND_VERSION_PATCH}")
endif()

if(NTL_INCLUDES)
  file(GLOB NTL_HEADERS "${NTL_INCLUDES}/version.h")
  foreach(ntl_header_filename ${NTL_HEADERS})
    file(READ "${ntl_header_filename}" _ntl_version_header)
    string(REGEX MATCH
      "define[ \t]*NTL_VERSION[ \t]*\"*([0-9]+\\.[0-9]+\\.[0-9]+)\"*" _ntl_version_match
      "${_ntl_version_header}")
    if(_ntl_version_match)
      set(NTL_VERSION "${CMAKE_MATCH_1}")
      break()
    endif()
  endforeach()

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT NTL_VERSION)
    set(NTL_VERSION_OK FALSE)
    message(STATUS "NTL version was not detected")
  elseif(${NTL_VERSION} VERSION_LESS ${NTL_FIND_VERSION})
    set(NTL_VERSION_OK FALSE)
    message(STATUS "NTL version ${NTL_VERSION} found in ${NTL_INCLUDES}, "
                   "but at least version ${NTL_FIND_VERSION} is required")
  else()
    set(NTL_VERSION_OK TRUE)
  endif()
endif()


find_library(NTL_LIBRARIES ntl PATHS ${LIB_INSTALL_DIR} HINTS ${NTL_POSSIBLE_LIB_PATHS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NTL DEFAULT_MSG
                                  NTL_INCLUDES NTL_LIBRARIES NTL_VERSION_OK)

#mark_as_advanced(NTL_INCLUDES NTL_LIBRARIES)
