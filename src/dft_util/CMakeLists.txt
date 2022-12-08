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
     check_fthaw.F90
     cwrap_drvnq.F90
     driver.F90
     drvdft.F90
     drvemb.F90
     fexp.F90
     functionals.F90
     get_denergy.F90
     get_exfac.F90
     ksdft_info.F90
     libxc_parameters.F90
     ndsd_ts.F90
     nucatt.F90
     ofembed.F90
     ofe_print.F90
     overlap.F90
     vemb_exc_states.F90
     vt_lim.F90
     wrap_drvnq.F90
     xlambda.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     ksdft_info.F90
     ofembed.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
