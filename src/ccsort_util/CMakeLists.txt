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
     abpack.F90
     action_ccsort.F90
     addintabc1.F90
     addintab.F90
     addinta.F90
     addpqij.F90
     ampack.F90
     ccsort_global.F90
     ccsort_grc0.F90
     ccsort_t3grc0.F90
     chkinp_ccsort.F90
     daopen.F90
     daread.F90
     dawri.F90
     dawrite.F90
     dawrtmap.F90
     dawrtmediate.F90
     deflength.F90
     deft3par.F90
     esb_ic_1.F90
     esb_ic_2.F90
     esb_ic_3.F90
     expandfok.F90
     expmpq.F90
     exppqij.F90
     exppsb.F90
     fokupdate1.F90
     fokupdate2.F90
     initintabc1.F90
     inittemp.F90
     initwrk.F90
     ireorg1.F90
     ireorg2.F90
     ireorg3.F90
     ireorg.F90
     mkabpqmap.F90
     mkaddress.F90
     mkampq.F90
     mkampqmap.F90
     mkintsta.F90
     mkmapampq.F90
     mkmappqij.F90
     mktempanam.F90
     mod1.F90
     mod2.F90
     mreorg1.F90
     mreorg.F90
     prinppn.F90
     rdinppn.F90
     rdtraone.F90
     reorg.F90
     t3intpck1.F90
     t3intpck2.F90
     t3reorg.F90
     unpackk.F90
     unpackk_ic_1.F90
     unpackk_ic_2.F90
     unpackk_ic_3.F90
     unpackk_pck.F90
     unpackk_zr.F90
     vf.F90
     zasun.F90
     zasun_pck.F90
     zasun_zr.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
