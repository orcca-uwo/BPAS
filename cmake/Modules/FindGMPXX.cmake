# Try to find the GMP library
# https://gmplib.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(GMPXX 6.0.0)
# to require version 6.0.0 to newer of GMP.
#
# Once done this will define
#
#  GMPXX_FOUND - system has GMP lib with correct version
#  GMPXX_INCLUDES - the GMP include directory
#  GMPXX_LIBRARIES - the GMP library
#  GMPXX_VERSION - GMP version
#
# Copyright (c) 2018 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.


find_path(GMPXX_INCLUDES NAMES gmpxx.h PATHS $ENV{GMPDIR} ${INCLUDE_INSTALL_DIR})

# Set GMP_FIND_VERSION to 5.1.0 if no minimum version is specified
if(NOT GMPXX_FIND_VERSION)
  if(NOT GMPXX_FIND_VERSION_MAJOR)
    set(GMPXX_FIND_VERSION_MAJOR 5)
  endif()
  if(NOT GMPXX_FIND_VERSION_MINOR)
    set(GMPXX_FIND_VERSION_MINOR 1)
  endif()
  if(NOT GMPXX_FIND_VERSION_PATCH)
    set(GMPXX_FIND_VERSION_PATCH 0)
  endif()
  set(GMPXX_FIND_VERSION
    "${GMPXX_FIND_VERSION_MAJOR}.${GMPXX_FIND_VERSION_MINOR}.${GMPXX_FIND_VERSION_PATCH}")
endif()

#gmpxx relies on the underlying gmp version, so lets assume its okay
set(GMPXX_VERSION_OK TRUE)

find_library(GMPXX_LIBRARIES gmpxx PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMPXX DEFAULT_MSG
                                  GMPXX_INCLUDES GMPXX_LIBRARIES GMPXX_VERSION_OK)

#mark_as_advanced(GMP_INCLUDES GMP_LIBRARIES)
