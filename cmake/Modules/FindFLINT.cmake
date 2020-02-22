# Try to find the Flint library
# https://github.com/wbhart/flint2
#
#
# Once done this will define
#
#  FLINT_FOUND - system has FLINT lib with correct version
#  FLINT_INCLUDES - the FLINT include directory
#  FLINT_LIBRARIES - the FLINT library
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.

set(FLINT_POSSIBLE_INCLUDE_PATHS
  /usr/include
  /usr/include/flint
  /usr/local/include
  /usr/local/include/flint
)

set(FLINT_POSSIBLE_LIB_PATHS
  /usr/lib
  /usr/local/lib
)

find_path(FLINT_INCLUDES NAMES flint.h PATHS ${INCLUDE_INSTALL_DIR} HINTS ${FLINT_POSSIBLE_INCLUDE_PATHS})
find_path(FLINT_INCLUDES NAMES fmpq.h PATHS ${INCLUDE_INSTALL_DIR} HINTS ${FLINT_POSSIBLE_INCLUDE_PATHS})

#lets assume its okay if it has fmpq.h
set(FLINT_VERSION_OK TRUE)

find_library(FLINT_LIBRARIES flint PATHS ${LIB_INSTALL_DIR} HINTS ${FLINT_POSSIBLE_LIB_PATHS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG
                                  FLINT_INCLUDES FLINT_LIBRARIES FLINT_VERSION_OK)

# mark_as_advanced(FLINT_INCLUDES FLINT_LIBRARIES)
