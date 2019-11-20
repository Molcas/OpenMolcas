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
# Copyright (C) 2017, Ignacio Fdez. Galván                             *
#***********************************************************************
if (NOT DEFINED BASE_DIR)
        set (BASE_DIR ${PROJECT_SOURCE_DIR})
endif ()
file (RELATIVE_PATH tooldir ${BASE_DIR} ${CMAKE_CURRENT_LIST_DIR})
set (OUTPUT_PATH ${PROJECT_BINARY_DIR}/${tooldir})

file (GLOB py_sources "*.py")
list (REMOVE_ITEM py_sources "${CMAKE_CURRENT_LIST_DIR}/pack.py")

if (NOT PYTHONINTERP_FOUND)
        set (dir_error "disabled (no appropriate python interpreter found)" PARENT_SCOPE)
        return()
endif ()

set (PYTHONINTERP_FOUND "${PYTHONINTERP_FOUND}" PARENT_SCOPE)

execute_process (
        COMMAND ${PYTHON_EXECUTABLE} imports.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE PYMOLCAS_IMPORT_OUTPUT
        RESULT_VARIABLE PYMOLCAS_IMPORT_RC
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

if (PYMOLCAS_IMPORT_RC)
        message (WARNING "Some python modules are not available: ${PYMOLCAS_IMPORT_OUTPUT}")
        set (dir_error "disabled" PARENT_SCOPE)
        return()
endif ()

set (PYMOLCAS_TARGET "${OUTPUT_PATH}/pymolcas_")
set (PYMOLCAS_TARGET "${PYMOLCAS_TARGET}" PARENT_SCOPE)

add_custom_command (OUTPUT ${PYMOLCAS_TARGET}
COMMAND ${PYTHON_EXECUTABLE} export.py ${PYMOLCAS_TARGET}
                WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        DEPENDS ${py_sources}
        )

add_custom_target (pymolcas ALL
        DEPENDS ${PYMOLCAS_TARGET}
        )

add_custom_command (TARGET pymolcas
POST_BUILD
        COMMAND ${OPENMOLCAS_DIR}/sbin/install_pymolcas.sh ${PYMOLCAS_TARGET}
        )
