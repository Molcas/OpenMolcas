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
     ao2mo_1e.F90
     clsfls_mcpdft.F90
     compute_mcpdft_energy.F90
     energy_mcwfn.F90
     fock_update.F90
     get_charge.F90
     get_coulomb.F90
     get_hcore.F90
     inppri.F90
     mcpdft.F90
     mcpdft_init.F90
     mcpdft_input.F90
     mcpdft_output.F90
     mspdft.F90
     mspdftgrad.F90
     mspdftgrad_util.F90
     mspdft_util.F90
     ontop_functional.F90
     opnfls_mcpdft.F90
     p2_contraction.F90
     print_mcpdft_2.F90
     proc_inpx.F90
     readvc.F90
     ref_energy.F90
     savefock_mspdft.F90
     savefock_pdft.F90
     validate_wfn.F90
     write_pdft_job.F90
    )

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
