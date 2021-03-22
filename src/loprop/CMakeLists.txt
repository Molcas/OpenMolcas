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
     build_amatrix.F90
     comp_f.F90
     dynamic_properties.F90
     fluctuating.F90
     make_fluctuating_charges.F90
     move_polar.F90
     polar.F90
     print_local.F90
     print_mpprop.F90
     readin_polar.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
