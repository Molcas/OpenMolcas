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

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)

# Workaround for gcc spurious warning about large variables,
# versions are approximate
# (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98411)
# This could also be removed if large arrays are made allocatable
if (BIGOT AND NOT OPENMP AND
    CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "10.4")
  target_compile_options (${prog}_obj PRIVATE "-fmax-stack-var-size=8388608")
endif ()
