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
     dminv.F90
     pot.F90
     prwf.F90
     simpsn.F90
     sort_pot.F90
     spectc.F90
     spline.F90
     vibinp.F90
     vibmat.F90
     vibrot.F90
     vibrotmain.F90
     vibrot_globals.F90
     vibtrm.F90
     wr_vibrot_info1.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
