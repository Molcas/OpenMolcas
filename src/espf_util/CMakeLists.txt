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
     bdvgrd.F90
     calcdb.F90
     calcdt.F90
     drv1_pot.F90
     drvespf.F90
     drvpot.F90
     espf_analysis.F90
     espf_energy.F90
     espf.F90
     espf_global.F90
     espf_grad.F90
     espf_mltp.F90
     espf_write.F90
     extnuc.F90
     fetch_qmmm.F90
     gmx_wrapper.c
     h1_espf.F90
     initb.F90
     initdb.F90
     ipl_espf.F90
     la_morok.F90
     mkgrid.F90
     mmcount.F90
     no_espf.F90
     opt_mmo.F90
     pcm_cavity.F90
     pnt.F90
     potind.F90
     prepare.F90
     readin_espf.F90
     rungromacs.F90
     runtinker.F90
     selfen.F90
     write_qmmm.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
