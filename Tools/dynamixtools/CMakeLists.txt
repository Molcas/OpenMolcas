#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2015,2016,2019, Ignacio Fdez. Galván                   *
#***********************************************************************
if (NOT DEFINED BASE_DIR)
        set (BASE_DIR ${PROJECT_SOURCE_DIR})
endif ()
file (RELATIVE_PATH tooldir ${BASE_DIR} ${CMAKE_CURRENT_LIST_DIR})
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${tooldir})

if (PYTHONINTERP_FOUND)

        foreach (prog IN ITEMS dynamixtools)
                file (READ ${CMAKE_CURRENT_LIST_DIR}/${prog}.py content)
                file (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog}.py "#!/usr/bin/env ${PYTHON_EXECUTABLE}\n")
                file (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog}.py "${content}")
                file (COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog}.py
                      DESTINATION ${EXECUTABLE_OUTPUT_PATH}
                      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                )
        endforeach (prog)

else ()

        message (WARNING "No appropriate python interpreter found, dynamixtools will be disabled")

endif ()
