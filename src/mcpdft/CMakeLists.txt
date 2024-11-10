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
     ao2mo_1e.F90
     chkinp.f
     clsfls_mcpdft.F90
     energy_mcwfn.F90
     fock_update.f
     hybridpdft.F90
     inppri.F90
     main.F90
     mcpdft.f
     mcpdft_init.F90
     mcpdft_input.F90
     mcpdft_output.F90
     msctl.f
     mspdft.F90
     mspdft_util.F90
     mspdftgrad.F90
     mspdftgrad_util.F90
     ontop_functional.F90
     opnfls_mcpdft.F90
     print_mcpdft_2.F90
     proc_inpX.f
     readvc.F90
     ref_energy.F90
     savefock_mspdft.F90
     write_pdft_job.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
