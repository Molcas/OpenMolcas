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
     chorpa_drv.F90
     chorpa_motra.F90
     chorpa_motra_reordercmo.F90
     rpa_checkintegralrepresentation.F90
     rpa_cleanup.F90
     rpa.F90
     rpa_freezer.F90
     rpa_frz.F90
     rpa_globals.F90
     rpa_iuhf.F90
     rpa_ppinp.F90
     rpa_prinp.F90
     rpa_rdinp.F90
     rpa_rdorb.F90
     rpa_rdorb_frominporb.F90
     rpa_rdorb_fromrunfile.F90
     rpa_rdrun.F90
     rpa_readintegerinput.F90
     rpa_setinc.F90
     rpa_setintegralrepresentation.F90
     rpa_setup_add_info.F90
     rpa_setup.F90
     rpa_warn.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
