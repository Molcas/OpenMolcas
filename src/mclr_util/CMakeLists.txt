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
     adads1.F90
     adadst.F90
     ads1.F90
     adst.F90
     adtor2_mclr.F90
     annstr.F90
     cidia4.F90
     cidia_mclr.F90
     cidia_sa.F90
     cidia_td.F90
     cisize.F90
     cmp2cn.F90
     cndet_mclr.F90
     cnford.F90
     cnfstr_mclr.F90
     cnhcn2.F90
     cnhcnm.F90
     cntost.F90
     confg2.F90
     crestr.F90
     csdiag_mclr.F90
     csdtmt.F90
     csdtvc_mclr_2.F90
     csfdet_mclr.F90
     csfinf.F90
     degvec.F90
     densi2_mclr.F90
     detctl.F90
     detstr_mclr.F90
     dihdj2_mclr.F90
     dxtyp.F90
     frmdsc_mclr.F90
     gasdn2_mclr.F90
     gatrmt.F90
     genstr_mclr.F90
     getcnf_mclr.F90
     geth1i_mclr.F90
     getinc_abs.F90
     getinc_abs_td.F90
     getinc_abt.F90
     getint_mclr.F90
     getint_td.F90
     gsaxpy.F90
     gsbbd1_mclr.F90
     gsbbd2a_mclr.F90
     gsbbd2b_mclr.F90
     gsdnbb2_mclr.F90
     gsttbl_mclr.F90
     gt1dia_mclr.F90
     gth1es_mclr.F90
     gtijkl_mclr.F90
     gtjk_mclr.F90
     gtsttp_1.F90
     gtsttp_2.F90
     gtsttps_1.F90
     gtsttps_2.F90
     h0.F90
     h0mat_mclr.F90
     h0csf.F90
     iaibcm_mclr.F90
     icispc.F90
     icisps.F90
     iel13.F90
     incoos.F90
     intcsf.F90
     intdia.F90
     ioctp2_mclr.F90
     ippage.F90
     isymst_mclr.F90
     memstr.F90
     mkcipre.F90
     mkp1.F90
     mkp1inv.F90
     mxresc.F90
     ncastr_mclr.F90
     newtyp_mclr.F90
     ngeth1.F90
     nrasdt.F90
     nstrso_mclr.F90
     numst3.F90
     numst4_mclr.F90
     nxtblk_mclr.F90
     orbinf_mclr.F90
     orbord.F90
     ordstr_mclr.F90
     pmplfm.F90
     psttbl_mclr.F90
     rassg4.F90
     rsbb1e_mclr.F90
     rsbb2a_mclr.F90
     rsbb2bn_mclr.F90
     rsmxmn_mclr.F90
     rssbcbn_mclr.F90
     rssbcbn_td.F90
     sa_prec.F90
     sa_prec2.F90
     scarmt.F90
     scdtc2_mclr.F90
     sdcmrf_mclr_1.F90
     sdcmrf_mclr_2.F90
     sdcmrf_pre.F90
     sigmavec.F90
     skickj_mclr.F90
     spncom_mclr.F90
     strinf.F90
     strtyp.F90
     sxtyp.F90
     todsc_mclr.F90
     tripk2_1.F90
     tripk2_2.F90
     weight_mclr.F90
     zbase.F90
     zoos.F90
     zsmcl.F90
     zstinf_mclr.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     ippage.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
