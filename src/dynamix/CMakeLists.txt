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
     cre_dyn.F90
     dynamix.F90
     dynamix_lib.F90
     dyn.fh
     getmassdx.F90
     grom.F90
     init_dynamix.F90
     MD.fh
     nhcThermo.F90
     prgm.fh
     project_in_for.F90
     project_in_vel.F90
     project_out_for.F90
     project_out_vel.F90
     randomGauss.F90
     rdinp_dynamix.F90
     readin_dynamix.F90
     velver_dump.F90
     velver_first.F90
     velver_second.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
