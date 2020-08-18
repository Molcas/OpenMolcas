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

file (GLOB data *.DAT)

foreach (datafile ${data})
        file (COPY ${datafile} DESTINATION ${EXECUTABLE_OUTPUT_PATH})
endforeach (datafile)

file (GLOB sources *.f90)

add_executable (caspt2_mini.exe ${sources})
target_link_libraries (caspt2_mini.exe molcas)
target_link_libraries (caspt2_mini.exe ${EXTERNAL_LIBRARIES})

