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
     check_use.F90
     clrruncache.F90
     crdrun.F90
     cwrrun.F90
     cxrdrun.F90
     cxwrrun.F90
     drdrun.F90
     dumprun.F90
     dwrrun.F90
     dxrdrun.F90
     dxwrrun.F90
     ffrun.F90
     ffxrun.F90
     fin_run_use.F90
     get_bastype.F90
     get_carray.F90
     get_cmo.F90
     get_coord_all.F90
     get_coord_all_.F90
     get_coord_full.F90
     get_coord_new_all.F90
     get_d1ao_var.F90
     get_d1sao_var.F90
     get_darray.F90
     get_darray_chk.F90
     get_dscalar.F90
     get_dscalar_.F90
     get_grad_full.F90
     get_iarray.F90
     get_iscalar.F90
     get_iscalar_.F90
     get_lscalar.F90
     get_mass_all.F90
     get_mass.F90
     get_nadc.F90
     get_name_all.F90
     get_name.F90
     get_name_full.F90
     get_natoms_all.F90
     get_natoms_full.F90
     get_nmode_all.F90
     get_nuc_charge_all.F90
     get_temp.F90
     get_zarray.F90
     gxrdrun.F90
     gxwrrun.F90
     gzrwrun.F90
     init_run_use.F90
     irdrun.F90
     iwrrun.F90
     ixrdrun.F90
     ixwrrun.F90
     mkrun.F90
     namerun.F90
     opnrun.F90
     put_analhess.F90
     put_carray.F90
     put_coord_full.F90
     put_coord_new.F90
     put_darray.F90
     put_dscalar.F90
     put_iarray.F90
     put_iscalar.F90
     put_lscalar.F90
     put_nadc.F90
     put_nucattr.F90
     put_pc_coord_new.F90
     put_temp.F90
     put_zarray.F90
     qpg_carray.F90
     qpg_darray.F90
     qpg_dscalar.F90
     qpg_iarray.F90
     qpg_iscalar.F90
     qpg_zarray.F90
     run2h5.F90
     runfile_data.F90
     runfile_procedures.F90
     store_energies.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     runfile_procedures.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
