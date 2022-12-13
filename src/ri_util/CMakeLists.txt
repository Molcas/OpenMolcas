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
     a_3c_qv_s.F90
     cd_ainv.F90
     cd_ainv_inner.F90
     cho_factor.F90
     cho_get_grad.F90
     cho_reorder_ri.F90
     cho_ri_final.F90
     cho_ri_putinfo.F90
     cho_ri_swapvecunit.F90
     cho_x_setab.F90
     compute_a_jk.F90
     compute_auxvec.F90
     compute_b.F90
     compute_txy.F90
     compute_v12.F90
     contract_zpk_tpxy.F90
     create_chunk.F90
     destroy_chunk.F90
     drv2el_3center_ri.F90
     drv2el_ri_diag.F90
     drvg1_2center_ri.F90
     drvg1_3center_ri.F90
     drvg1_ri.F90
     fix_coeff.F90
     free_tsk2.F90
     gen_qvec.F90
     get_auxiliary_shells.F90
     get_chunk.F90
     get_maxdg.F90
     get_mxos.F90
     get_pivot_idx.F90
     get_pivot_idx_w.F90
     indsft_ri_2.F90
     indsft_ri_3.F90
     inicho_ri.F90
     inicho_ri_xtras.F90
     init_tsk2.F90
     in_place_diag.F90
     in_place_square.F90
     integral_ri_2.F90
     integral_ri_3.F90
     integral_ricd.F90
     inv_cho_factor.F90
     iramax.F90
     mk_acd_accd_shells.F90
     mk_coeffs.F90
     mk_dummy_shell.F90
     mk_indkl.F90
     mk_iso2ind.F90
     mk_list2.F90
     mk_ricd_shells.F90
     mk_ri_shells.F90
     mk_tint_p.F90
     mk_tvt.F90
     mk_tvtf.F90
     modify_tint_p.F90
     mult_3c_qv_s.F90
     mult_rijk_qkl.F90
     mult_vk_qv_s.F90
     mult_with_q_caspt2.F90
     mult_with_q_mp2.F90
     mult_zp_qv_s.F90
     nmemam.F90
     nsize_3c.F90
     nsize_rv.F90
     o2n.F90
     ofembed_dmat.F90
     off_diagonal.F90
     pget1_cd2.F90
     pget1_cd3.F90
     pget1_ri2.F90
     pget1_ri3.F90
     pget2_cd2.F90
     pget2_cd3.F90
     pget2_ri2.F90
     pget2_ri3.F90
     pivot_mat.F90
     plf_ri_2.F90
     plf_ri_3.F90
     plf_ricd.F90
     post_2center_ri.F90
     put_chunk.F90
     remap_u_k.F90
     remap_v_k.F90
     remove_high_exponents.F90
     renorm2.F90
     renorm2_inner.F90
     reord_vk.F90
     restore_mat.F90
     ri_glob.F90
     ri_procedures.F90
     ri_xdiag.F90
     rm_auxshell.F90
     rsv_tsk2.F90
     set_cho_adrvec.F90
     setchoindx_ri.F90
     setup_aux.F90
     setup_aux_inner.F90
     sort_mat.F90
     square_a.F90
     termcho_ri.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     ri_glob.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
