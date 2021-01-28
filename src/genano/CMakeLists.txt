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
     check_genano.F90
     cpovlp.F90
     figcls.F90
     figopn.F90
     figprt.F90
     free_genano.F90
     genano.F90
     genano_globals.F90
     init_genano.F90
     inpctl_genano.F90
     mkano.F90
     mktype.F90
     nophase.F90
     prhead.F90
     proj1.F90
     proj2.F90
     rdcmo.F90
     sort_genano.F90
     sphave.F90
     updens.F90
     uporb.F90
     zlcase.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.txt)
