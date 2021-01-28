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
     check_genano.f90
     cpovlp.f90
     figcls.f90
     figopn.f90
     figprt.f90
     free_genano.f90
     genano.f90
     genano_globals.F90
     init_genano.f90
     inpctl_genano.f90
     mkano.f90
     mktype.f90
     nophase.f90
     prhead.f90
     proj1.f90
     proj2.f90
     rdcmo.f90
     sort_genano.f90
     sphave.f90
     updens.f90
     uporb.f90
     zlcase.f90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
