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
     main.F90
     cre_gsswfn.F90
     fckbyint.F90
     fmod1n.F90
     fmod1s.F90
     fockoper.F90
     golowdin.F90
     gopickup.F90
     gopop.F90
     gosort.F90
     guessorb.F90
     guessorb_global.F90
     initgo.F90
     inpctl_guessorb.F90
     virt_space.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
