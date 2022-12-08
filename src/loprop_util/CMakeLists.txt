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
     check_bond.F90
     coretopoint.F90
     cutoff_error.F90
     desymmetrize.F90
     diff_aux1.F90
     diff_mothergoose.F90
     diff_numerical.F90
     diff_thrsmul.F90
     elpointpot.F90
     elpot.F90
     epotpoint.F90
     error_for_t.F90
     find_dipole_center.F90
     get_density_matrix.F90
     golden2.F90
     golden.F90
     gramschmidt.F90
     infotomp.F90
     init_loprop.F90
     levmarquart.F90
     localize_loprop_drv.F90
     localize_loprop.F90
     local_properties.F90
     loprop_arrays.F90
     loprop.F90
     loprop_print.F90
     lowdin_lp.F90
     min_mult_error.F90
     mnbrak2.F90
     mnbrak.F90
     move_ec.F90
     move_prop.F90
     multipole_e.F90
     multipole_expansion.F90
     pickpoints.F90
     print_t_values.F90
     read_h0.F90
     read_multipole_int.F90
     reexpand.F90
     rotate_dipole.F90
     solvea.F90
     step1.F90
     step2.F90
     step3.F90
     step4.F90
     transmu.F90
     triangularize.F90
     ttotal.F90
     warnings_lp.F90
     wegotthis.F90
     yougetthis.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     loprop_arrays.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
