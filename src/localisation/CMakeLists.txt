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
     bitmap_localisation_atom.F90
     bitmap_localisation.F90
     chk_input.F90
     domain_localisation.F90
     genbmp_localisation.F90
     gengnu_localisation.F90
     get_cnos.F90
     get_etwo_act.F90
     getinfo_localisation_0.F90
     get_nat_lorb.F90
     localisation.F90
     localisation_globals.F90
     localise_iterative.F90
     localise_noniterative.F90
     loc_nat_orb.F90
     pao_analysis.F90
     primo_localisation.F90
     rdvec_localisation.F90
     readinp_localisation.F90
     testloc.F90
     wrvec_localisation.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
