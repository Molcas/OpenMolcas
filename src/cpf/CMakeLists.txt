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
     abcd_cpf.F90
     abci_cpf.F90
     ab_cpf.F90
     ai_cpf.F90
     alloc_cpf.F90
     apprim.F90
     count_cpf.F90
     cpfctl.F90
     cpf.F90
     cpf_global.F90
     cupdate.F90
     decomp.F90
     dens_cpf.F90
     densct_cpf.F90
     diagc_cpf.F90
     diag_cpf.F90
     diagct_cpf.F90
     diis_cpf.F90
     dsq2.F90
     epsbis.F90
     epsprim.F90
     faibj_cpf.F90
     fij_cpf.F90
     ifock.F90
     iijj_cpf.F90
     ijij_cpf.F90
     ijkl_cpf.F90
     indmat_cpf.F90
     ipo_cpf.F90
     mabcd.F90
     mabci.F90
     mab.F90
     mai.F90
     mdiagc.F90
     mdsq2.F90
     mfaibj.F90
     mfij.F90
     mijkl.F90
     mpsq2.F90
     natct.F90
     natorb_cpf.F90
     next.F90
     npset.F90
     onect.F90
     prwf_cpf.F90
     psq2.F90
     readin_cpf.F90
     restart_cpfmcpf.F90
     scatter.F90
     sdci_cpf.F90
     secord.F90
     sing.F90
     solve.F90
     sorta_cpf.F90
     sortb_cpf.F90
     sort_cpf.F90
     start_cpf.F90
     thetset.F90
     tpqset.F90
     twoct.F90
     vecsum_cpfmcpf.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
