#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************

include_guard (GLOBAL)

get_filename_component (prog ${CMAKE_CURRENT_LIST_DIR} NAME)

file (GLOB sources *.f *.f90 *.F *.F90 *.c)

add_executable (${prog}.exe ${sources})

if (DEFINED ${prog}_deplibs)
  target_link_libraries (${prog}.exe ${${prog}_deplibs})
endif ()

# Info for the main project
set (${prog}_src ${CMAKE_CURRENT_LIST_DIR})
set (${prog}_sources ${sources})
