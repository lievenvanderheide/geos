################################################################################
# Part of CMake configuration for GEOS
#
# Copyright (C) 2020 Martin Davis
#
# This is free software; you can redistribute and/or modify it under
# the terms of the GNU Lesser General Public Licence as published
# by the Free Software Foundation.
# See the COPYING file for more information.
################################################################################

add_executable(geosop
  GeosOp.cpp
  GeometryOp.cpp
  )

target_link_libraries(geosop PRIVATE geos geos_c)

install(TARGETS geosop
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)
