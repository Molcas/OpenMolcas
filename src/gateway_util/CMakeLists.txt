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
     agin.F90
     amatrix.F90
     amfi_info.F90
     basis2run.F90
     basis_info.F90
     calcamt.F90
     calc_lstvec.F90
     center_info.F90
     clean_bname.F90
     cow.F90
     decode.F90
     dmp_efp.F90
     dmpinf.F90
     dmwdot.F90
     drvn0.F90
     dstchk.F90
     external_centers.F90
     fill_rinfo1.F90
     find_basis_set.F90
     fix_fockop.F90
     gateway_global.F90
     gateway_info.F90
     gen_geolist.F90
     geolist.F90
     geonew.F90
     geonew_pc.F90
     getbs.F90
     getecp.F90
     getpam.F90
     index_center.F90
     index_nosym.F90
     info2runfile.F90
     lbl2nr.F90
     ldf_checkconfig.F90
     ldf_checkthrs.F90
     ldf_setoptionflag.F90
     mergebs.F90
     merge_lists.F90
     misc.F90
     mk_chdisp.F90
     name_to_lm.F90
     nrmlx.F90
     nrmlz.F90
     oeisg.F90
     one_int.F90
     order_axis.F90
     ordexp1.F90
     ordexp.F90
     ovlmp.F90
     period.F90
     phase_info.F90
     print_basis2.F90
     print_basis.F90
     print_geometry.F90
     print_isotopes.F90
     print_opinfo.F90
     print_symmetry.F90
     processrp.F90
     process_weights.F90
     pso_stuff.F90
     rdbsl.F90
     rdctl_seward.F90
     rdmx.F90
     reorder_gw.F90
     ricd_info.F90
     rigrot.F90
     saddle.F90
     set_fake_eris.F90
     settargetaccuracy_ldf.F90
     setup_offao.F90
     sizes_of_seward.F90
     soao_info.F90
     soctl_seward.F90
     symmetry_info.F90
     vexch.F90
     vqr.F90
     xyz.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     basis_info.F90
     center_info.F90
     external_centers.F90
     gateway_global.F90
     gateway_info.F90
     phase_info.F90
     pso_stuff.F90
     ricd_info.F90
     sizes_of_seward.F90
     soao_info.F90
     symmetry_info.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
