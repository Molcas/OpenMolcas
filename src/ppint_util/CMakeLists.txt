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
     ang1.F90
     ang2.F90
     assemble_ppgrd.F90
     cortab.F90
     facab.F90
     hermit.F90
     hrecur.F90
     hroot.F90
     lmnvgn.F90
     molcas_order.F90
     ppgrd.F90
     ppint_arrays.F90
     ppint.F90
     ppmem.F90
     ppmmg.F90
     pseud1.F90
     pseud2.F90
     pseudo.F90
     ptwt.F90
     qbess.F90
     qcomp.F90
     qpasy.F90
     rad1.F90
     rad2.F90
     ssibfn.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
