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
     ffpt.F90
     ffpt_global.F90
     mkcom.F90
     prdiop.F90
     prinp.F90
     ptadd.F90
     ptdipo.F90
     ptefgr.F90
     ptefld.F90
     ptglbl.F90
     ptokt0.F90
     ptokt1.F90
     ptquad.F90
     ptrela.F90
     rd1int.F90
     rdinp.F90
     selectloc.F90
     transnow.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
