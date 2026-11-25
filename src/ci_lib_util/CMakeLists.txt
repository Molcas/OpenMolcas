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
     ci_input.F90
     ci_lib.F90
     ci_lib_checks.F90
     ci_lib_core.F90
     ci_lib_defs.F90
     ci_lib_info.F90
     ci_lib_states.F90
     ci_sd_csf_codec.F90
     ci_sg_lucia_bridge.F90
     ci_sg_ug2sg.F90
     ci_verify_options.F90
     orbital_spaces.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     ci_input.F90
     ci_lib.F90
     ci_lib_checks.F90
     ci_lib_core.F90
     ci_lib_defs.F90
     ci_lib_info.F90
     ci_lib_states.F90
     ci_sd_csf_codec.F90
     ci_sg_lucia_bridge.F90
     ci_verify_options.F90
     orbital_spaces.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
