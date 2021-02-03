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
     fmm_2darray_sort.F90
     fmm_aux_qlm_builder.F90
     fmm_boundary.F90
     fmm_box_builder.F90
     fmm_box_packer.F90
     fmm_box_utils.F90
     fmm_car_to_sph.F90
     fmm_driver.F90
     fmm_global_paras.F90
     fmm_integral_utils.F90
     fmm_interface.F90
     fmm_J_builder.F90
     fmm_local_search.F90
     fmm_multiple_T_worker.F90
     fmm_multipole_ints.F90
     fmm_multi_T_buffer.F90
     fmm_proc_selector.c
     fmm_qlm_builder.F90
     fmm_qlm_utils.F90
     fmm_scale_T_buffer.F90
     fmm_scheme_builder.F90
     fmm_shell_pairs.F90
     fmm_sort_paras.F90
     fmm_sort_T_pairs.F90
     fmm_stats.F90
     fmm_T_buffer.F90
     fmm_T_contractors.F90
     fmm_T_pair_builder.F90
     fmm_T_pair_mould.F90
     fmm_T_pair_tests.F90
     fmm_tree_buffer.F90
     fmm_T_worker.F90
     fmm_utils.F90
     fmm_Vff_driver.F90
     fmm_W_buffer.F90
     fmm_W_contractors.F90
     fmm_W_pair_builder.F90
     fmm_W_worker.F90
)

set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.txt)
