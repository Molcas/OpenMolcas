#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2026, Okan Koeksal                                     *
#***********************************************************************

set (sources
     main.F90
     eda_data.F90
     eda_parallel.F90
     eda_input.F90
     eda_reference.F90
     eda_integrals.F90
     eda_energy.F90
     eda.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
