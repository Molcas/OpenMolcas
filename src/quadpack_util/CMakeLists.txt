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
     d1mach.F90
     dqage.F90
     dqag.F90
     dqagie.F90
     dqagi.F90
     dqelg.F90
     dqk15.F90
     dqk15i.F90
     dqk21.F90
     dqk31.F90
     dqk41.F90
     dqk51.F90
     dqk61.F90
     dqpsrt.F90
     xerror.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
