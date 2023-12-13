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
     amfi.F90
     amfi_global.F90
     angular.F90
     buildcoul.F90
     cartonex.F90
     cartoney.F90
     cartonez.F90
     chngcont.F90
     contandmult.F90
     contcasaoo.F90
     contcasaso.F90
     contcasb1oo.F90
     contcasb1so.F90
     contcasb2oo.F90
     contcasb2so.F90
     contcascoo.F90
     contcascso.F90
     cont.F90
     contone.F90
     contract.F90
     couple3j.F90
     daxpint.F90
     drv_amfi.F90
     finite.F90
     gen1overr3.F90
     gencouldim.F90
     gencoul.F90
     genovlp.F90
     genpowers.F90
     genprexyz13.F90
     genprexyz14.F90
     genprexyz15a.F90
     gentkin.F90
     getaos2.F90
     getaos.F90
     getcg.F90
     getlimit.F90
     getocc_ao.F90
     inidf.F90
     initfrac.F90
     initired.F90
     kindiag.F90
     kinemat.F90
     lmdepang.F90
     mcheckxy.F90
     mcheckz.F90
     mkangl0.F90
     mkanglmin.F90
     prefac.F90
     readbas.F90
     regge3j.F90
     symtrafo.F90
     tkinet.F90
     tosigx.F90
     tosigy.F90
     tosigz.F90
     trans_amfi.F90
     transcon.F90
     two2mean12a.F90
     two2mean12b.F90
     two2mean13.F90
     two2mean34a.F90
     two2mean34b.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
