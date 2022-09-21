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
     alaska_banner.F90
     alaska.F90
     alaska_info.F90
     alaska_super_driver.F90
     annihil_rho.F90
     caspt2_grad_btamp.F90
     caspt2_grad_fwdcnt.F90
     chk_numerical.F90
     cho_alaska_rdinp.F90
     drvdftg.F90
     drvemb_.F90
     drvembg.F90
     drvg1.F90
     drvh1_emb.F90
     drvn1_emb.F90
     drvn1.F90
     inputg.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
