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
     cgetmem.c
     c_null_alloc.c
     getmem.F90
     inimem.F90
     mma_module.F90
     stdalloc.F90
    )

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     mma_module.F90
     stdalloc.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
