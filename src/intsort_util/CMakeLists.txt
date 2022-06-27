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
     mkord.F90
     mksrt0.F90
     mksrt1.F90
     mksrt2.F90
     savebin.F90
     sort0.F90
     sort1a.F90
     sort1b.F90
     sort1c.F90
     sort2a.F90
     sort2b.F90
     sort2.F90
     sort3.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     sort_data.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
