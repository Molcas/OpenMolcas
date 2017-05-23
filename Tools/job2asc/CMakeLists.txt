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
if (NOT DEFINED BASE_DIR)
        set (BASE_DIR ${PROJECT_SOURCE_DIR})
endif ()
file (RELATIVE_PATH tooldir ${BASE_DIR} ${CMAKE_CURRENT_LIST_DIR})
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${tooldir})

foreach (prog IN ITEMS job2asc asc2job)
        add_executable(${prog}.exe ${prog}.f)
        target_link_libraries (${prog}.exe molcas)
        target_link_libraries (${prog}.exe ${EXTERNAL_LIBRARIES})
endforeach (prog)
