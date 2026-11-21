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
     cidia_sa.F90
     cidia_td.F90
     cgs_mod.F90
     clsfls_mclr.F90
     dwder.F90
     exphinvv.F90
     e2.F90
     e2_td.F90
     get_cmo_mclr.F90
     h0.F90
     input_mclr.F90
     ippage.F90
     isrot.F90
     mclr_data.F90
     mclr_procedures.F90
     mkcipre.F90
     mkp1.F90
     mkp1inv.F90
     oitd.F90
     prinp_mclr.F90
     pcm_grad.F90
     reload.F90
     rddens.F90
     sa_prec.F90
     sa_prec2.F90
     sg_setup_mclr.F90
     tcmo.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     input_mclr.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
