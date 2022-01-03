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
     cho_cas_drv.F90
     cho_caspt2_getbasenm.F90
     cho_eval_waxy.F90
     cho_fmcscf.F90
     cho_fmo_red.F90
     cho_fock_dft_red.F90
     cho_fock_rassi.F90
     cho_fock_rassi_x.F90
     cho_focktwo.F90
     cho_focktwo_red.F90
     cho_fscf.F90
     cho_ftwo_mo.F90
     cho_geth1.F90
     cho_get_mo.F90
     cho_getshfull.F90
     cho_lk_casscf.F90
     cho_lk_maxvecperbatch.F90
     cho_lk_mclr.F90
     cho_lk_rassi.F90
     cho_lk_rassi_x.F90
     cho_lk_scf.F90
     cho_lk_screeningthreshold.F90
     #cho_lr_mos.F90    These two files are #include'd in
     #choras_drv.F90    the fock_util_interface module
     cho_rassi_twxy.F90
     choscf_mem.F90
     cho_sum.F90
     cho_x_quit.F90
     coul_dmb.F90
     findmax.F90
     fmmfck.F90
     focktwo_drv.F90
     focktwo.F90
     fock_util_global.F90
     fock_util_interface.F90
     fold.F90
     fold_mat.F90
     fold_tmat.F90
     ftwo_drv.F90
     ftwo.F90
     get_tuvx.F90
     getumat_t1.F90
     mk_ishp_rs.F90
     order_arrays.F90
     rdchovec.F90
     set_nna.F90
     set_nnbsf.F90
     swap_full2rs.F90
     swap_rs2full.F90
     swap_tosqrt.F90
     thouless_t1.F90
     tra2a.F90
     tra2b.F90
     tra2c.F90
     tra_ctl2.F90
     tractl2.F90
     tractl_drv.F90
     tradrv.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
