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
     active_space_solver_cfg.F90
     qcmaquis_info.F90
     qcmaquis_interface_wrapper.F90
     refwfn.F90
)

set (modfile_list
     active_space_solver_cfg.F90
     qcmaquis_info.F90
     qcmaquis_interface_wrapper.F90
     refwfn.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.txt)
