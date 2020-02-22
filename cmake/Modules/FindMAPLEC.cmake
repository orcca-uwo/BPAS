# Try to find the Maple library
# https://maplesoft.com/products/maple/
#
#
# Once done this will define
#
#  MAPLEC_FOUND - system has MAPLEC lib with correct version
#  MAPLEC_INCLUDES - the MAPLEC include directory
#  MAPLEC_LIBRARIES - the MAPLEC library
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Redistribution and use is allowed according to the terms of the BSD license.

set(MAPLEC_POSSIBLE_INCLUDE_PATHS
  /usr/include
  /usr/local/include
  $ENV{HOME}/maple2019
  $ENV{HOME}/maple2019/extern/include
  /opt/maple2019
  /opt/maple2019/extern/include
  $ENV{HOME}/maple2018
  $ENV{HOME}/maple2018/extern/include
  /opt/maple2018
  /opt/maple2018/extern/include
  $ENV{HOME}/maple2017
  $ENV{HOME}/maple2017/extern/include
  /opt/maple2017
  /opt/maple2017/extern/include
  $ENV{HOME}/maple2016
  $ENV{HOME}/maple2016/extern/include
  /opt/maple2016
  /opt/maple2016/extern/include
  $ENV{HOME}/maple2015
  $ENV{HOME}/maple2015/extern/include
  /opt/maple2015
  /opt/maple2015/extern/include
)

set(MAPLEC_POSSIBLE_LIB_PATHS
  /usr/lib
  /usr/local/lib
  $ENV{HOME}/maple2019
  $ENV{HOME}/maple2019/bin.X86_64_LINUX
  /opt/maple2019
  /opt/maple2019/bin.X86_64_LINUX
  $ENV{HOME}/maple2018
  $ENV{HOME}/maple2018/bin.X86_64_LINUX
  /opt/maple2018
  /opt/maple2018/bin.X86_64_LINUX
  $ENV{HOME}/maple2017
  $ENV{HOME}/maple2017/bin.X86_64_LINUX
  /opt/maple2017
  /opt/maple2017/bin.X86_64_LINUX
  $ENV{HOME}/maple2016
  $ENV{HOME}/maple2016/bin.X86_64_LINUX
  /opt/maple2016
  /opt/maple2016/bin.X86_64_LINUX
  $ENV{HOME}/maple2015
  $ENV{HOME}/maple2015/bin.X86_64_LINUX
  /opt/maple2015
  /opt/maple2015/bin.X86_64_LINUX
)

find_path(MAPLEC_INCLUDES NAMES maplec.h PATHS ${INCLUDE_INSTALL_DIR} HINTS ${MAPLEC_POSSIBLE_INCLUDE_PATHS})

#lets assume its okay
set(MAPLEC_VERSION_OK TRUE)

find_library(MAPLEC_LIBRARIES maplec PATHS ${LIB_INSTALL_DIR} HINTS ${MAPLEC_POSSIBLE_LIB_PATHS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MAPLEC DEFAULT_MSG
                                  MAPLEC_INCLUDES MAPLEC_LIBRARIES MAPLEC_VERSION_OK)

#mark_as_advanced(MAPLEC_INCLUDES MAPLEC_LIBRARIES)
