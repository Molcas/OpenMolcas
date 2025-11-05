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


if (Python_FOUND AND Python_VERSION VERSION_GREATER_EQUAL 3.9)

  # Install ga_driver.py as an executable script
  foreach (prog IN ITEMS ga_driver.py)
    file (COPY ${CMAKE_CURRENT_LIST_DIR}/${prog}
          DESTINATION ${EXECUTABLE_OUTPUT_PATH}
          FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    install (PROGRAMS
             ${EXECUTABLE_OUTPUT_PATH}/${prog}
             DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}
    )

  endforeach (prog)

  # Install GA_mod module
  set (GA_MOD_SRC ${CMAKE_CURRENT_LIST_DIR}/src/GA_mod)
  set (GA_MOD_OUTDIR ${EXECUTABLE_OUTPUT_PATH}/src/GA_mod)
  foreach (progs IN ITEMS __init__.py auxiliary_functions.py config.py crossover.py extend_ordering.py gen_ref_dicts.py GUGA_diag.py measure_fitness.py population.py run_GA.py sampling.py)
    file (COPY ${GA_MOD_SRC}/${progs}
          DESTINATION ${GA_MOD_OUTDIR}
          FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    install (PROGRAMS
             ${GA_MOD_OUTDIR}/${progs}
             DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}/src/GA_mod
    )
  endforeach (progs)

  # Install FCIDUMP_tools module
  set (FCIDUMP_SRC ${CMAKE_CURRENT_LIST_DIR}/src/FCIDUMP_tools)
  set (FCIDUMP_OUTDIR ${EXECUTABLE_OUTPUT_PATH}/src/FCIDUMP_tools)
  foreach (progs IN ITEMS __init__.py IntegralClass.py xyz2heisenberg.py)
    file (COPY ${FCIDUMP_SRC}/${progs}
          DESTINATION ${FCIDUMP_OUTDIR}
          FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    install (PROGRAMS
             ${FCIDUMP_OUTDIR}/${progs}
             DESTINATION ${CMAKE_INSTALL_PREFIX}/${tooldir}/src/FCIDUMP_tools
    )
  endforeach (progs)

else ()

  message (WARNING "No appropriate Python interpreter found, GA4GUGA will be disabled (Python >= 3.9 required)")

endif ()
