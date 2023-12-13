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
     c1diis.F90
     c2diis.F90
     chkpnt.F90
     con_opt.F90
     convrg.F90
     dfp.F90
     drvuph.F90
     dstinf.F90
     eu.F90
     find_distance.F90
     fixhess.F90
     free_updmask.F90
     g_nrm.F90
     init_updmask.F90
     kriging_update.F90
     main.F90
     mep_dir.F90
     msp.F90
     mxlbls.F90
     newcar_kriging.F90
     newh.F90
     newh_mod.F90
     newq.F90
     outzmat.F90
     prepare_kriging.F90
     qnr.F90
     reset_thrgrd.F90
     rlxctl.F90
     rotate_bp.F90
     rs_i_rfo.F90
     rs_p_rfo.F90
     slapaf_banner.F90
     slapaf.F90
     slstatus.F90
     ts_bfgs.F90
     update_h.F90
     update_inner.F90
     update_kriging.F90
     update_sl.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
