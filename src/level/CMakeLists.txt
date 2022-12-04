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
     level.f
     read_input.F90
     AF3X3LEV.f
     alfAS.f
     level.f
     masses.f
     potgen.f
     prepot.f
     scattlenDER.f
     scecorAS.f
     schrqAS.f
     widthAS.f
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
