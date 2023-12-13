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
     cre_out.F90
     cre_prep.F90
     cut_matrices.F90
     equation.F90
     equation_sph.F90
     get_dipole.F90
     get_dm0.F90
     get_hcsf.F90
     get_vsoc.F90
     hamdens.F90
     integrators.F90
     kab.F90
     k_external.F90
     pop.F90
     prepare_decay.F90
     propagate.F90
     propagate_sph.F90
     pulse.F90
     read_input.F90
     read_prep.F90
     read_rassd.F90
     read_rassisd.F90
     rhodyn_data.F90
     rhodyn.F90
     rhodyn_init.F90
     rhodyn_utils.F90
     soci.F90
     uci.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
