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
     embedding_global.F90
     embpotene.F90
     embpotfreemem.F90
     embpotinit.F90
     embpotkernel.F90
     embpotmem.F90
     embpotoutput.F90
     embpotrdrun.F90
     embpotwrrun.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     embedding_global.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
