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
     abdata.F90
     ass1a.F90
     ass1b.F90
     ass1.F90
     ass2a.F90
     ass2b.F90
     ass2.F90
     ass3a.F90
     ass3b.F90
     ass3.F90
     assg1.F90
     assg1_mck.F90
     assg2.F90
     cff2d.F90
     cff2dq.F90
     cff2ds.F90
     closer.F90
     distg1.F90
     distg1x.F90
     distg2.F90
     drvrys.F90
     exp_1.F90
     exp_2.F90
     fake.F90
     her_rw.F90
     hrr2da.F90
     hrr2da_mck.F90
     hrr2db.F90
     hrr2db_mck.F90
     hrrctl.F90
     hrrctl_mck.F90
     memrg1.F90
     memrg2.F90
     memrys.F90
     memrys_g.F90
     modu2.F90
     pppp.F90
     ppps.F90
     ppss.F90
     pr2d.F90
     psps.F90
     psss.F90
     read_rysrw.F90
     rs2dgh.F90
     rs2dmm.F90
     rtswgh.F90
     rys01.F90
     rys11.F90
     rys22.F90
     rys2d.F90
     rys2dg.F90
     rys33.F90
     rys44.F90
     rys55.F90
     rys66.F90
     rys77.F90
     rys88.F90
     rys99.F90
     rysef0.F90
     rysef1.F90
     rysef2.F90
     rysef3.F90
     rysef4.F90
     rysef.F90
     rys.F90
     rysg1.F90
     rysg2.F90
     rysscratch.F90
     sether.F90
     setupr.F90
     setup_rw.F90
     sppp.F90
     spsp.F90
     sspp.F90
     sssp.F90
     ssss.F90
     teri1.F90
     teri.F90
     teris.F90
     terisq.F90
     tnai1.F90
     tnai.F90
     vcff2d.F90
     vrys2d.F90
     vrys2dm.F90
     vrys_rw.F90
     vrysrw.F90
     xcff2d.F90
     xrys2d.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     her_rw.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
