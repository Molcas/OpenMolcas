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
     drv1el.F90
     drv2el.F90
     drv_fck.F90
     input_seward.F90
     nemo_opt1.F90
     output1_seward.F90
     seward_banner.F90
     seward.F90
     tr_prm_cnt.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
