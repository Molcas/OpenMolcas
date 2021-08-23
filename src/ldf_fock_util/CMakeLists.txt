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
     ldf_auxbasvectornorm.F90
     ldf_blockmatrixnorm.F90
     ldf_computecoulombintermediates0.F90
     ldf_computecoulombintermediates.F90
     ldf_computeu.F90
     ldf_fci.F90
     ldf_fock_coulomberroranalysis.F90
     ldf_fock_coulomberror.F90
     ldf_fock_coulombonly0_1.F90
     ldf_fock_coulombonly0_2.F90
     ldf_fock_coulombonly0_3.F90
     ldf_fock_coulombonly0_.F90
     ldf_fock_coulombonly0.F90
     ldf_fock_coulombonly_.F90
     ldf_fock_coulombonly.F90
     ldf_fock_coulombonly_xidi.F90
     ldf_fock_coulombupperbound.F90
     ldf_fock_coulombupperbound_full.F90
     ldf_fock_coulombupperboundnorm.F90
     ldf_fock_coulombupperboundnorm_full.F90
     ldf_fock_cub.F90
     ldf_ftst.F90
     ldf_fvifc.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
