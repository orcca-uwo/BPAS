# Try to find the CDDGMP library
# https://github.com/cddlib/cddlib
#
#
# Once done this will define
#
#  CDDGMP_FOUND - system has CDDGMP lib with correct version
#  CDDGMP_INCLUDES - the CDDGMP include directory
#  CDDGMP_LIBRARIES - the CDDGMP library
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.

set(CDDGMP_POSSIBLE_INCLUDE_PATHS
  /usr/include
  /usr/local/include
  /usr/include/cdd
)

set(CDDGMP_POSSIBLE_LIB_PATHS
  /usr/lib
  /usr/local/lib
)

find_path(CDDGMP_INCLUDES NAMES cddmp.h PATHS ${INCLUDE_INSTALL_DIR} HINTS ${CDDGMP_POSSIBLE_INCLUDE_PATHS})

#lets assume its okay
set(CDDGMP_VERSION_OK TRUE)


find_library(CDDGMP_LIBRARIES cddgmp PATHS ${LIB_INSTALL_DIR} HINTS ${CDDGMP_POSSIBLE_LIB_PATHS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CDDGMP DEFAULT_MSG
                                  CDDGMP_INCLUDES CDDGMP_LIBRARIES CDDGMP_VERSION_OK)

#mark_as_advanced(CDDGMP_INCLUDES CDDGMP_LIBRARIES)
