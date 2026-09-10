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

# Shared grid utilities: evaluate MOs/densities on a set of grid points.
# Used by grid_it and by embedding_util (FDE), so they live in their own
# *_util library rather than being referenced across directories.
set (sources
     moeval.F90
     outmo.F90
     power.F90
     soadpt.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
