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
     blas_mod.F
     blas_mod_f90.F90
     close_blas.F90
     delayed_blas.F90
     delayed_lapack.F90
     delayed_legacy.F90
     dlfcn.F90
     initialize_blas.F90
     iso_c_utilities.F90
     lapack_mod.F
     lapack_mod_f90.F90
     legacy_mod.F90
     link_blas.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
