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
     basisconsistency.F90
     basisreader.F90
     basistbl.F90
     basistype.F90
     basistype_mod.F90
     crprod.F90
     foundatomicnumber.F90
     get_basistype.F90
     get_dnumber.F90
     get_inumber.F90
     pick_words.F90
     read_xbas.F90
     stdsewinput.F90
     transtbl.F90
     vec.F90
     xmatreader.F90
     xmatrixconverter.F90
     zmatconv.F90
     zmatconv_mod.F90
     zmatreader.F90
     zmatrixconverter.F90
     zmatrixconverter_gw.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
