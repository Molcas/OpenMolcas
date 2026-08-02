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
     bctwoham.F90
     eafmod.F90
     ga_replace_ma.c
     gadgop.F90
     gadgop_scal.F90
     gadupl.F90
     gaigop.F90
     gaigop_scal.F90
     gainit.F90
     gatskl.F90
     gatskl_zero.F90
     gasync.F90
     gaterminate.F90
     nap_time.c
     par_range.F90
     pardiv.F90
     pfget_ascii_f.F90
     rsvtsk.F90
     scalapack_mod.F90
     slaveschdir.c
     sync_data.F90
     sync_th.F90
     task_manager.F90
     tlist_mod.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     eafmod.F90
     scalapack_mod.F90
     task_manager.F90
     tlist_mod.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
