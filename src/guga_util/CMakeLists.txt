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
     fmmm.F90
     fmul2.F90
     guga_util_global.F90
     icunp.F90
     jacscf.F90
     jsunp.F90
     mtrans.F90
     order.F90
     siadd.F90
     squar2.F90
     squar.F90
     squarm.F90
     squarn.F90
     tradd.F90
     vneg.F90
     wr_guga_info.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
