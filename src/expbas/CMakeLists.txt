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
     cre_expwfn.F90
     desym.F90
     driverbas.F90
     expandbas.F90
     expbas.F90
     info_expbas.F90
     readinp_expbas.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
