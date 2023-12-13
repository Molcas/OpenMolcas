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
     fmm_j_builder.F90
     fmm_local_search.F90
     fmm_multiple_t_worker.F90
     fmm_multipole_ints.F90
     fmm_multi_t_buffer.F90
     fmm_proc_selector.c
     fmm_qlm_builder.F90
     fmm_qlm_utils.F90
     fmm_scale_t_buffer.F90
     fmm_scheme_builder.F90
     fmm_shell_pairs.F90
     fmm_sort_paras.F90
     fmm_sort_t_pairs.F90
     fmm_stats.F90
     fmm_t_buffer.F90
     fmm_t_contractors.F90
     fmm_t_pair_builder.F90
     fmm_t_pair_mould.F90
     fmm_t_pair_tests.F90
     fmm_tree_buffer.F90
     fmm_t_worker.F90
     fmm_utils.F90
     fmm_vff_driver.F90
     fmm_w_buffer.F90
     fmm_w_contractors.F90
     fmm_w_pair_builder.F90
     fmm_w_worker.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
