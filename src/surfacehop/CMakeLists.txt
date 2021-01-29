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
     getIphInfo.F90
     initial_surfacehop.F90
     rdinp_surfacehop.F90
     restart_surfacehop.F90
     surfacehop.F90
     surfacehop_globals.F90
     tully.F90
     tully_variables.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
