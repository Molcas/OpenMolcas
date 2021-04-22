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
     charge_grid_it.F90
     drvmo.F90
     dumpm2msi.F90
     gridexpandselect.F90
     grid_it.F90
     grid_it_globals.F90
     input_grid_it.F90
     linenr.F90
     makepab.F90
     moevaldel.F90
     moevalder.F90
     moeval.F90
     mycoor.F90
     mygetkey.F90
     noblanks.F90
     one_charge.F90
     opengrid.F90
     outmo.F90
     pickorb.F90
     power.F90
     print_ddnorm.F90
     printheader.F90
     printline.F90
     printtitles.F90
     prt_lusc.c
     prtlusendgrid.F90
     save_ddnorm.F90
     soadpt.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
