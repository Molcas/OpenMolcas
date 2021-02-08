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
     checkquater.F90
     cross.F90
     generategeoms.F90
     getkandc.F90
     modangle.F90
     normalizevec.F90
     printgeom.F90
     quater.F90
     quaterfinish.F90
     quater_globals.F90
     quaterinit.F90
     quaterreport.F90
     quaterrotation.F90
     quatersetup.F90
     quatersolve.F90
     quater_sub.F90
     rdinput.F90
     readgeo.F90
     rotategeom.F90
     rotategeoms.F90
     setmatrix.F90
     setvect.F90
     setvecttrans.F90
     translategeom.F90
     translategeoms.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
