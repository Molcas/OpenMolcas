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
     sg_citra.F90
     sg_cxinit.F90
     sg_d1mat.F90
     sg_epq_psi.F90
     sg_mstow.F90
     sg_pkwlk.F90
     sg_prwf.F90
     sg_print.F90
     sg_upkwlk.F90
     sguga.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     sguga.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
