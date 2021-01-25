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
     checkquater.f90
     cross.f90
     generategeoms.f90
     getkandc.f90
     modangle.f90
     normalizevec.f90
     printgeom.f90
     quater.f90
     quaterfinish.f90
     quater_globals.f90
     quaterinit.f90
     quaterreport.f90
     quaterrotation.f90
     quatersetup.f90
     quatersolve.f90
     quater_sub.F90
     rdinput.f90
     readgeo.f90
     rotategeom.f90
     rotategeoms.f90
     setmatrix.f90
     setvect.f90
     setvecttrans.f90
     translategeom.f90
     translategeoms.f90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
