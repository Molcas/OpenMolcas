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
     acore.F90
     clr2.F90
     clrbuf.F90
     cmbn2dc.F90
     cmbnacb1.F90
     cmbnacb2.F90
     cmbneldot.F90
     cmbnel.F90
     cmbns1_mck.F90
     cmbns2a.F90
     cmbns2b.F90
     cmbns2.F90
     cmbnt1_mck.F90
     cmbnt2.F90
     cnt1el2.F90
     cnt1el.F90
     cnthlf_mck.F90
     cntrct_mck.F90
     coreb.F90
     crsph_mck.F90
     ctldns.F90
     ctrlmo.F90
     dan.F90
     dede_mck.F90
     derctr.F90
     din.F90
     dot1el2.F90
     dot1el.F90
     drvel1.F90
     drvetc.F90
     drvg2.F90
     drvh1_mck.F90
     drvh2.F90
     drvk2_mck.F90
     drvn1_mck.F90
     drvn2.F90
     elgrddot.F90
     elgrd.F90
     elmem.F90
     fckacc_mck.F90
     hssprt.F90
     inputh.F90
     k2loop_mck.F90
     knegrd_mck.F90
     knehss.F90
     knemem_mck.F90
     knemmh.F90
     ltocore.F90
     ltosph.F90
     lu2lu.F90
     m1grd_mck.F90
     m1hss.F90
     m1kernel.F90
     m1mm1.F90
     m1mmh.F90
     main.F90
     makemo.F90
     mckinley_banner.F90
     mckinley.F90
     mckinley_global.F90
     mck_interface.F90
     mkfck.F90
     moacc.F90
     mult_sro.F90
     na2mem.F90
     nagrd_mck.F90
     nahss.F90
     namem_mck.F90
     nammh.F90
     nonatwo.F90
     nucind.F90
     opnfls_mckinley.F90
     ovrgrd_mck.F90
     ovrhss.F90
     ovrmem_mck.F90
     ovrmmh.F90
     pckint_mck.F90
     pckmo2.F90
     pickmo.F90
     prjgrd_mck.F90
     prjhss.F90
     prjmm1.F90
     prjmmh.F90
     psoao0_h.F90
     psoao2.F90
     request_mclr_run.F90
     rtocore.F90
     rtosph.F90
     schint_mck.F90
     screen_mck.F90
     smadna.F90
     sort_mck.F90
     srogrd_mck.F90
     srohss.F90
     sromm1.F90
     srommh.F90
     sttstc.F90
     supermac.F90
     symado_mck2.F90
     symado_mck.F90
     tcrtnc_h.F90
     tnchlf_h.F90
     translation.F90
     twodns.F90
     twoel_mck.F90
     wrdisk.F90
     wrhdsk.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
