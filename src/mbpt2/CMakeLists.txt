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
     bjai.F90
     build_mp2dens_old.F90
     check_amp2.F90
     cho_sosmp2_col.F90
     cho_sosmp2_decchk.F90
     cho_sosmp2_decdrv.F90
     cho_sosmp2_drv.F90
     cho_sosmp2_energy.F90
     cho_sosmp2_setup.F90
     cho_sosmp2_setup_prt.F90
     compute_shanks.F90
     delghost_mbpt.F90
     finish_wdensity.F90
     fnomp2_drv.F90
     fno_mp2.F90
     fnomp2_putinf.F90
     freezer.F90
     frzdel.F90
     gamma_new.F90
     get_tr_dab.F90
     lovmp2_drv.F90
     lovmp2_putinf.F90
     mbpt2_global.F90
     mp2ap.F90
     mp2dens_drv.F90
     mp2diag.F90
     mp2_driver.F90
     mp2gdens_setup.F90
     prinp_mbpt2.F90
     rdinp.F90
     rdint2_mp2.F90
     rdmbpt.F90
     rhs_mp2.F90
     rhs_mp2_help1.F90
     rhs_mp2_help2.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
