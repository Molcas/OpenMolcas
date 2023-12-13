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


if (Python_FOUND AND Python_VERSION VERSION_GREATER_EQUAL 3.6)

  foreach (prog IN ITEMS auger_main.py)
    file (READ ${CMAKE_CURRENT_LIST_DIR}/${prog} content)
    file (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog} "#!/usr/bin/env ${Python_EXECUTABLE}\n")
    file (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog} "${content}")
    file (COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog}
          DESTINATION ${EXECUTABLE_OUTPUT_PATH}
          FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    install (PROGRAMS
             ${EXECUTABLE_OUTPUT_PATH}/${prog}
             DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}
    )

  endforeach (prog)

  set (AUGER_SRC ${CMAKE_CURRENT_LIST_DIR}/auger_oca)
  set (AUGER_SRC_OUTDIR ${EXECUTABLE_OUTPUT_PATH}/auger_oca)
  foreach (progs IN ITEMS auger_driver.py basis_id.py initi.py __init__.py input_parse.py mbs.py my_variables.py OCA.dat oca_integrals.py orbtransf_sz.py rt2mzz.py symm.py)
    file (READ ${AUGER_SRC}/${progs} content)
    file (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/auger_oca/${progs} "#!/usr/bin/env ${Python_EXECUTABLE}\n")
    file (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/auger_oca/${progs} "${content}")
    file (COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/auger_oca/${progs}
          DESTINATION ${AUGER_SRC_OUTDIR}
          FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    install (PROGRAMS
             ${AUGER_SRC_OUTDIR}/${progs}
             DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}/auger_oca
    )
  endforeach (progs)

else ()

  message (WARNING "No appropriate Python interpreter found, Auger OCA will be disabled (Python >= 3.6 required)")

endif ()
