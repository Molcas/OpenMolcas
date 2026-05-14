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
     abc_axes.F90
     atens.F90
     barrier.F90
     calcmagn1.F90
     calcmagn2.F90
     cart2tensor.F90
     cart2tensor_minus.F90
     cdet.F90
     chi.F90
     complex_1_trace2.F90
     decomp_aniso.F90
     dev.F90
     diag_c2.F90
     diag_r2.F90
     dmatrix.F90
     finddetr.F90
     g_high.F90
     g_high_1.F90
     generate_isotrop_site.F90
     hdir.F90
     hdir2.F90
     io_data.F90
     jkqpar.F90
     magn.F90
     magn_no_mf.F90
     magn_zj_par.F90
     mean_field.F90
     mean_field_all.F90
     mean_field_exch.F90
     moments.F90
     mu_order.F90
     newjkqpar.F90
     norder.F90
     operators.F90
     pa_pseudo.F90
     pa_sort.F90
     plot_barrier.F90
     plot_mh_no_exp.F90
     plot_mh_with_exp.F90
     plot_xt_no_exp.F90
     plot_xt_with_exp.F90
     prbar.F90
     pseudospin.F90
     read_aniso_old_exch.F90
     read_binary_aniso.F90
     read_formatted_aniso.F90
     read_formatted_aniso_old.F90
     read_formatted_aniso_poly.F90
     read_formatted_aniso_poly_new.F90
     read_formatted_new_aniso.F90
     read_hdf5_all.F90
     read_hdf5_init.F90
     read_hdf5_poly.F90
     real_1_trace2.F90
     recover_exch_ham_from_liviu_ito.F90
     recover_exch_ham_from_naoya_ito.F90
     recover_exch_ham_from_stevens_eso.F90
     reverse.F90
     rotmom.F90
     rotmom2.F90
     rtrace.F90
     set_an.F90
     set_knm.F90
     set_nm.F90
     set_t.F90
     spin.F90
     spin_phase.F90
     spin_phase2.F90
     tensor2cart.F90
     tensor2cart_minus.F90
     trace.F90
     trace2.F90
     trace_exch.F90
     trace_exch2.F90
     transham.F90
     utils.F90
     utmu.F90
     utmu2.F90
     utmul.F90
     write_binary_aniso.F90
     write_formatted_aniso.F90
     write_formatted_aniso_poly.F90
     write_new_formatted_aniso.F90
     zeem_sa.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
