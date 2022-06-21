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
     main.F90
     arrange_orbital_molcas.F90
     gugadrt_active_drt.F90
     gugadrt_check_rcas3.F90
     gugadrt_dbl_downwalk.F90
     gugadrt_dbl_upwalk.F90
     gugadrt_ext_downwalk.F90
     gugadrt.F90
     gugadrt_global.F90
     gugadrt_gugafinalize.F90
     gugadrt_mole_inf.F90
     gugadrt_njexcit.F90
     gugadrt_paras_calculate.F90
     gugadrt_rcas.F90
     gugadrt_ref_gfs.F90
     gugadrt_rst.F90
     gugainit.F90
     packnod.F90
     redabkm.F90
     upacknod.F90
     writedrt.F90
     wrtabkm.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
