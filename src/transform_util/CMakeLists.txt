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
     cho_genc.F90
     cho_gene.F90
     cho_motra.F90
     cho_motra_inner.F90
     chomp2_gene.F90
     chomp2_traa.F90
     chomp2_tractl.F90
     chomp2_tras.F90
     chomp2_twoel.F90
     cho_traa.F90
     cho_tractl.F90
     cho_tra.F90
     cho_tras.F90
     cho_tr_drv.F90
     cho_twoel.F90
     coul.F90
     def_subblocke.F90
     def_tcvx.F90
     exch.F90
     filler.F90
     get_nvnode.F90
     lenint.F90
     local_triang.F90
     mem_est.F90
     mkl1.F90
     mkl2.F90
     mkl3.F90
     mklij.F90
     mksubs.F90
     printsquaremat.F90
     printtriangmat.F90
     prodsa_1.F90
     prodsa_2.F90
     prodsa_2t.F90
     prodss_1.F90
     put_tra_comm.F90
     rbuff_tra2.F90
     rbufo_tra2.F90
     rbuf_tra2.F90
     rdint2.F90
     setup_caspt2_tra.F90
     tr2nsa1.F90
     tr2nsa2.F90
     tr2nsa3.F90
     tr2nsb.F90
     tr2sq.F90
     tractl.F90
     transp_MOs.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
