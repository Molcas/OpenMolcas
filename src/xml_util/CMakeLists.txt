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
     xml_cdump.F90
     xml_cdumpa.c
     xml_cdumpb.c
     xml_cdumpc.c
     xml_cdumps.c
     xml_close.F90
     xml_closec.c
     xml_ddump.F90
     xml_ddumpc.c
     xml_idump.F90
     xml_idumpc.c
     xml_open.F90
     xml_openc.c
     xml_prspec.c
     xml_sdump.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
