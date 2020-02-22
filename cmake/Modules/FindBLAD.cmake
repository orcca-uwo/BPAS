# Try to find the BLAD library
# http://cristal.univ-lille.fr/~boulier/pmwiki/pmwiki.php?n=Main.BLAD
#
#
# Once done this will define
#
#  BLAD_FOUND - system has BLAD lib with correct version
#  BLAD_INCLUDES - the BLAD include directory
#  BLAD_LIBRARIES - the BLAD library
#  BLAD_VERSION - BLAD version
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.

set(BLAD_POSSIBLE_INCLUDE_PATHS
  /usr/include
)

set(BLAD_POSSIBLE_LIB_PATHS
  /usr/lib
  /usr/local/lib
)

find_path(BLAD_INCLUDES NAMES blad.h `PATHS ${INCLUDE_INSTALL_DIR} HINTS ${BLAD_POSSIBLE_INCLUDE_PATHS})


# Set BLAD_FIND_VERSION to 3.0.0 if no minimum version is specified
if(NOT BLAD_FIND_VERSION)
  if(NOT BLAD_FIND_VERSION_MAJOR)
    set(BLAD_FIND_VERSION_MAJOR 3)
  endif()
  if(NOT BLAD_FIND_VERSION_MINOR)
    set(BLAD_FIND_VERSION_MINOR 0)
  endif()
  if(NOT BLAD_FIND_VERSION_PATCH)
    set(BLAD_FIND_VERSION_PATCH 0)
  endif()
  set(BLAD_FIND_VERSION
    "${BLAD_FIND_VERSION_MAJOR}.${BLAD_FIND_VERSION_MINOR}.${BLAD_FIND_VERSION_PATCH}")
endif()

if(BLAD_INCLUDES)
  file(GLOB BLAD_HEADERS "${BLAD_INCLUDES}/ba0.h")
  foreach(blad_header_filename ${BLAD_HEADERS})
    file(READ "${blad_header_filename}" _blad_version_header)
    string(REGEX MATCH
      "define[ \t]*BLAD_VERSION[ \t]*\"*([0-9]+\\.[0-9]+\\.[0-9]+)\"*" _blad_version_match
      "${_blad_version_header}")
    if(_blad_version_match)
      set(BLAD_VERSION "${CMAKE_MATCH_1}")
      break()
    endif()
  endforeach()

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT BLAD_VERSION)
    set(BLAD_VERSION_OK FALSE)
    message(STATUS "BLAD version was not detected")
  elseif(${BLAD_VERSION} VERSION_LESS ${BLAD_FIND_VERSION})
    set(BLAD_VERSION_OK FALSE)
    message(STATUS "BLAD version ${BLAD_VERSION} found in ${BLAD_INCLUDES}, "
                   "but at least version ${BLAD_FIND_VERSION} is required")
  else()
    set(BLAD_VERSION_OK TRUE)
  endif()
endif()


find_library(BLAD_LIBRARIES blad PATHS ${LIB_INSTALL_DIR} HINTS ${BLAD_POSSIBLE_LIB_PATHS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLAD DEFAULT_MSG
                                  BLAD_INCLUDES BLAD_LIBRARIES BLAD_VERSION_OK)

#mark_as_advanced(BLAD_INCLUDES BLAD_LIBRARIES)
