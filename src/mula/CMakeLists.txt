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
     calc_r00.F90
     fcval.F90
     fgenerator.F90
     inoutmod.F90
     intensitymod.F90
     intforcefield.F90
     iotools.F90
     iscd_routines.F90
     isc_maketab2.F90
     isc_routines.F90
     linalg.F90
     logevec.F90
     matelmod.F90
     mk_nincdec.F90
     mkvibwind2.F90
     mula.F90
     mula_global.F90
     optmod.F90
     potdist.F90
     potkin.F90
     randommod.F90
     rottranrem.F90
     setupdipmat2.F90
     setupharmdip.F90
     setuphmat2.F90
     solvesort.F90
     tabmod.F90
     vibmod.F90
     writecartcoord.F90
     writefreq.F90
     writeintcoord.F90
     writeint.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
