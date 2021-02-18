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

set (sources
     main.c
     parnell_base.c
     parnell.c
     parnell_cmd.c
     parnell_collect.c
     parnell_copy.c
     parnell_exec.c
     parnell_init.c
     parnell_reduce.c
     parnell_remove.c
     parnell_replica.c
     parnell_rmlist.c
     parnell_scatter.c
     parnell_translate.c
     parnell_unlink.c
     parnell_wipe.c
)

get_filename_component (prog ${CMAKE_CURRENT_LIST_DIR} NAME)

add_executable (${prog}.exe ${sources})

if (DEFINED ${prog}_deplibs)
  target_link_libraries (${prog}.exe ${${prog}_deplibs})
endif ()

# Info for the main project
set (${prog}_src     ${CMAKE_CURRENT_LIST_DIR})
set (${prog}_sources ${sources})
