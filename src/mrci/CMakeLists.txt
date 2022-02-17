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
     abcd_mrci.F90
     abci_mrci.F90
     abd.F90
     ab_mrci.F90
     abtd.F90
     aid.F90
     ai_mrci.F90
     aitd.F90
     alloc_mrci.F90
     ci_select_mrci.F90
     count_mrci.F90
     cscale.F90
     csftra.F90
     cstart.F90
     dcorr.F90
     densct.F90
     diagc.F90
     diagct.F90
     faibj2.F90
     faibj3.F90
     faibj5.F90
     faibj.F90
     fijd.F90
     fij_mrci.F90
     fijtd.F90
     fmul2.F90
     hz.F90
     hzlp1.F90
     hzlp2.F90
     iijj.F90
     ijij.F90
     ijkl.F90
     indmat.F90
     ipo.F90
     jsunp.F90
     loop70.F90
     mkdao.F90
     mkhref.F90
     mktdao.F90
     mqct.F90
     mrci.F90
     mrci_global.F90
     mtrans.F90
     natorb_mrci.F90
     order.F90
     pkvec.F90
     pmatel.F90
     propct.F90
     prorb.F90
     prwf_mrci.F90
     readin_mrci.F90
     refci.F90
     sdci_mrci.F90
     seceq.F90
     secne.F90
     secular.F90
     siadd.F90
     sigma.F90
     sorta.F90
     sortb.F90
     sort_mrci.F90
     squar2.F90
     squar.F90
     squarm.F90
     squarn.F90
     tradd.F90
     upkvec.F90
     vneg.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
