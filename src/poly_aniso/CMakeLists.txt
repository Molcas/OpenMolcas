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
     aniso_lines_exchange3.F90
     aniso_lines_exchange9.F90
     chi_sum.F90
     dipol_exchange.F90
     dirvect.F90
     distance.F90
     dzyaloshinsky_morya_exchange.F90
     exchctl.F90
     fetch_init_const.F90
     fetch_neq.F90
     find_aniso_format.F90
     input_process.F90
     jcoeff.F90
     jito_exchange_int.F90
     jkqpar_naoya.F90
     ke_covalent.F90
     ke_exchange.F90
     kinetic_exchange.F90
     lines_exchange.F90
     magnetization_pa.F90
     momloc2.F90
     msum.F90
     pa_diagham.F90
     pa_preigen.F90
     pa_prham.F90
     poly_aniso.F90
     poly_aniso_open.F90
     popanalysis.F90
     prep_mom_exchange.F90
     pr_ito_int.F90
     project_exch.F90
     readin_poly.F90
     set_defaults.F90
     susceptibility_pa.F90
     torque_pa.F90
     xt_dmoverdh.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
