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
     add_infoc.c
     add_info.F90
     bragg_slater.F90
     charge.F90
     chkvec.F90
     chk_vec_uhf.F90
     CMakeLists.txt
     cmp_hfc.F90
     dens_if.F90
     dens_if_scf.F90
     dnuclearmass.F90
     ecp_shells.F90
     fileorb.F90
     freeze_default.F90
     freq_i.F90
     freq_molden.F90
     getdens.F90
     gfprnt_i.F90
     imostabundantisotope.F90
     inporbfmt.F90
     interf.F90
     inuclearchargefromsymbol.F90
     isoloop.F90
     isotopes.F90
     molden_dysorb.F90
     molden_interface.F90
     nat_bond_order.F90
     orb2tpidx.F90
     orb2tpidx_sym.F90
     orb2tpstr.F90
     orb2tpstr_sym.F90
     orbtype.F90
     primo.F90
     prop.F90
     prout.F90
     prppnt.F90
     prpt_.F90
     prpt.F90
     prpt_old.F90
     rdvec_.F90
     rdvec.F90
     rdvec_hdf5.F90
     s2calc.F90
     seek_n_destroy.F90
     spinav.F90
     tpidx2orb.F90
     tpidx2orb_sym.F90
     tpidx2tpstr.F90
     tpstr2orb.F90
     tpstr2orb_sym.F90
     tpstr2tpidx.F90
     vdwrad.F90
     vecsort.F90
     wrvec_.F90
     wrvec.F90
     xprop.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     isotopes.F90
     prppnt.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
