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
     gen_prim_density_matrix.F90
     get_density_matrix_mpprop.F90
     get_mpprop.F90
     get_mpprop_input.F90
     get_ocof.F90
     get_orbcen.F90
     get_orbe.F90
     get_polar.F90
     get_prim_atom_tab.F90
     get_prim_density_matrix.F90
     laurapol.F90
     mpprop.F90
     mpprop_globals.F90
     noverp.F90
     storempaslop.F90
     wr_cord.F90
     wr_files.F90
     wr_mpprop.F90
     wr_prop.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
