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
     main.F90
     crystalfield.F90
     crystalfield_1.F90
     fetch_data_runfile_all.F90
     fetch_data_runfile_init.F90
     individual_ranks.F90
     magnetization.F90
     newcf.F90
     readin_single.F90
     recover_cf.F90
     restart_check.F90
     restart_sa.F90
     single_aniso2.F90
     single_aniso_open.F90
     sort_kq.F90
     susceptibility.F90
     termcf.F90
     torque.F90
     xt_dmoverdh_single.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
