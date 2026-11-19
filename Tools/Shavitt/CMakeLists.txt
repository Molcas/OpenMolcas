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
# Copyright (C) 2026, Roland Lindh                                     *
#***********************************************************************

# Shavitt: the split graph of RASSCF wave functions (Python, no compilation).
# The script, the package shavittlib and the checks are copied to the build
# directory and installed under Tools/Shavitt.

if (NOT DEFINED BASE_DIR)
  set (BASE_DIR ${PROJECT_SOURCE_DIR})
endif ()
file (RELATIVE_PATH tooldir ${BASE_DIR} ${CMAKE_CURRENT_LIST_DIR})
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${tooldir})

if (Python_FOUND AND Python_VERSION VERSION_GREATER_EQUAL 3.8)

  # the script, with the shebang of the Python found by CMake
  file (STRINGS ${CMAKE_CURRENT_LIST_DIR}/shavitt.py first LIMIT_COUNT 1)
  file (READ ${CMAKE_CURRENT_LIST_DIR}/shavitt.py content)
  string (LENGTH "${first}\n" skip)
  string (SUBSTRING "${content}" ${skip} -1 content)
  file (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/shavitt.py "#!/usr/bin/env ${Python_EXECUTABLE}\n")
  file (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/shavitt.py "${content}")
  file (COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/shavitt.py
        DESTINATION ${EXECUTABLE_OUTPUT_PATH}
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
  install (PROGRAMS
           ${EXECUTABLE_OUTPUT_PATH}/shavitt.py
           DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}
  )

  # the package (with the logo) and the checks
  foreach (dir IN ITEMS shavittlib tests)
    file (COPY ${CMAKE_CURRENT_LIST_DIR}/${dir}
          DESTINATION ${EXECUTABLE_OUTPUT_PATH}
          PATTERN "__pycache__" EXCLUDE
    )
    install (DIRECTORY
             ${EXECUTABLE_OUTPUT_PATH}/${dir}
             DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}
             PATTERN "__pycache__" EXCLUDE
    )
  endforeach (dir)
  install (FILES
           ${CMAKE_CURRENT_LIST_DIR}/README.md
           DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}
  )

else ()

  message (WARNING "No appropriate Python interpreter found, Shavitt will be disabled (Python >= 3.8 required)")

endif ()
