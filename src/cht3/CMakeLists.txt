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
     block_interf.F90
     calc_mp2.F90
     check_create_klvab_t3_mem.F90
     check_loops.F90
     check_mat.F90
     cht3.F90
     cht3_global.F90
     cht3_procedures.F90
     create_klvab_t3.F90
     decomp2ind.F90
     defcommon.F90
     defparreord_t3.F90
     delf.F90
     ex23.F90
     exmap3_231.F90
     expa1_uhf.F90
     expa2_uhf.F90
     expand4_12.F90
     ext_o_32.F90
     gather_t2anti_blocked.F90
     gather_t2_blocked.F90
     gather_t2.F90
     gather_t2_fblocked.F90
     generate_juzekoe.F90
     gen_oovo.F90
     gen_vvoo_blocked.F90
     gen_vvoo.F90
     gen_vvvo.F90
     getrest_t3.F90
     getx_t3.F90
     grow_l1.F90
     grow_l2.F90
     grow_t2anti_blocked1.F90
     grow_t2anti_blocked2.F90
     grow_t2_blocked.F90
     grow_t2_fblocked1.F90
     grow_t2_fblocked2.F90
     grow_t2neq.F90
     grow_vvoo_blocked.F90
     grow_vvoo.F90
     grow_w3.F90
     inireord_t3.F90
     klvaa_vvv.F90
     main.F90
     map2_21_t3.F90
     map3_132_t3.F90
     map3_231_t3.F90
     map3_321_t3.F90
     map4_3412_t3.F90
     multi_opendir.F90
     multi_readir.F90
     multi_wridir.F90
     my_block.F90
     t3ampl_bti.F90
     t3_bta_aac.F90
     t3_bt_aaa.F90
     t3_bta_abc.F90
     t3_bt_aac.F90
     t3_bt_abc.F90
     t3_bt_acc.F90
     t3loopa.F90
     t3loopb.F90
     v_size_t3.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
