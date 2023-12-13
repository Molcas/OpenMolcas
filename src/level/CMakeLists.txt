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
     af3x3lev.F90
     alfas.F90
     cdjoelas.F90
     dampf.F90
     genint.F90
     level_common.F90
     level.F90
     level_rdinp.F90
     level_spline.F90
     levqad.F90
     levxpc.F90
     main.F90
     masses.F90
     matxel.F90
     plyintrp.F90
     potgen.F90
     prepot.F90
     qboundas.F90
     scattlendet.F90
     scecoras.F90
     schrqas.F90
     splint.F90
     widthas.F90
     zheevj3.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
