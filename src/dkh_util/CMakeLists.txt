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
     aodkhexp.F90
     bss_ts1e.F90
     copy_mag_ints.F90
     dkh_cofu.F90
     dkh_cofu_spec.F90
     dkh_geneu.F90
     dkh_ham.F90
     dkh_info.F90
     dkh_prop.F90
     dkh_ts1e.F90
     dkh_wgene.F90
     dkh_woplft.F90
     dkh_woprig.F90
     dkh_wspec.F90
     dkh_xpx.F90
     dkrelint_dp.F90
     dmxma.F90
     merge_mag_ints.F90
     x2c_makx.F90
     x2c_ts1e.F90
     xdr_dmatinv.F90
     xdr_dmatsqrt.F90
     xdr_fpfw.F90
     xdr_fpfwprop.F90
     xdr_ham.F90
     xdr_indx.F90
     xdr_info_local.F90
     xdr_local_ham.F90
     xdr_mkutls.F90
     xdr_prop.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     dkh_info.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
