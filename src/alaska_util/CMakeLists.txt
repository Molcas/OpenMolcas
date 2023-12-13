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
     cmbnm2.F90
     cmbnmlt1.F90
     cmbnrf1.F90
     cmbns1.F90
     cmbnt1.F90
     cmbnw1.F90
     cosgrd.F90
     create_grads.F90
     csfgrad.F90
     drvh1.F90
     finfld.F90
     grd_interface.F90
     int_parm_g.F90
     int_prep_g.F90
     knegrd.F90
     knemmg.F90
     m1grd.F90
     m1mmg.F90
     m2grd.F90
     m2mmg.F90
     mltgrd.F90
     mltgrdnuc.F90
     mltmmg.F90
     nagrd.F90
     nammg.F90
     oneel_g.F90
     ovrgrd.F90
     ovrmmg.F90
     pcmgrd.F90
     pcmmmg.F90
     prepre_g.F90
     prgrad.F90
     prjgrd.F90
     prjmmg.F90
     psoao1.F90
     query_grads.F90
     read_grad.F90
     rfgrd.F90
     rfmmg.F90
     screen_g.F90
     size_sob_g.F90
     soao_g.F90
     srogrd.F90
     srommg.F90
     store_grad.F90
     store_not_grad.F90
     twoel_g.F90
     welgrd.F90
     welmmg.F90
     xfdgrd.F90
     xfdmmg.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     grd_interface.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
