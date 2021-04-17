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
     autocut.F90
     cho_fock_motra.F90
     donei.F90
     fcin.F90
     ftwoi.F90
     init_motra.F90
     inpctl_motra.F90
     motra.F90
     motra_global.F90
     ortho_motra.F90
     orthox_motra.F90
     prinp.F90
     rd1int_motra.F90
     rdcmo_motra.F90
     rdinp_motra.F90
     rdrfld.F90
     tr1ctl.F90
     tr2ctl.F90
     tramo.F90
     traone_motra.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
