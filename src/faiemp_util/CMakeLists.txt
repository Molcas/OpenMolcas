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

set (sources
     addfragdens.F90
     drv2el_faiemp.F90
     drvg_faiemp.F90
     fragexpand.F90
     fragpcont.F90
     fragpgrd.F90
     fragpint.F90
     fragpmem.F90
     fragpmmg.F90
     getfragment.F90
     makedens.F90
     prepp_faiemp.F90
     reindexfrag.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
