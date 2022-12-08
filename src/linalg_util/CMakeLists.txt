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
     blockdiagonal_matrices.F90
     dgeadd.F90
     dgesub.F90
     dgetmi.F90
     dgetmo.F90
     dnaxpy.F90
     dndot.F90
     dtrans.F90
     dyax.F90
     dzaxpy.F90
     eigen_molcas.F90
     exp_schur.F90
     fzero.F90
     idmin.F90
     ilasrt.F90
     iunitmat.F90
     izero.F90
     linalg_mod.F90
     matinvert.F90
     nidiag.F90
     nidiag_new.F90
     not_dgeev.F90
     schur_skew.F90
     unitmat.F90
     vecphase.F90
     xeigen.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     blockdiagonal_matrices.F90
     linalg_mod.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
