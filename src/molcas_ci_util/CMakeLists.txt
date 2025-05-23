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
     alloc.F90
     cands.F90
     caspt2_global.F90
     check_invec.F90
     chksplit.F90
     cidia.F90
     ciovlp.F90
     ciselect.F90
     cisx.F90
     citrans.F90
     citrans_free.F90
     citrans_init.F90
     cndet.F90
     cnfstr.F90
     cnhcn.F90
     csdiag.F90
     csdtvc.F90
     csfbas.F90
     csfcount.F90
     cstart.F90
     cxinit.F90
     davctl.F90
     davctl_mod.F90
     david5.F90
     dblock.F90
     detstr2.F90
     dgmm2.F90
     diagord.F90
     dihdj_molcas.F90
     dmatdmat.F90
     dmrginfo.F90
     expcsf.F90
     explh2.F90
     faroald.F90
     faroald_free.F90
     faroald_init.F90
     fciqmc_interface.F90
     fndmnx.F90
     fockoc.F90
     gas_data.F90
     gatvcs.F90
     general_data.F90
     get_cm.F90
     get_cm_.F90
     get_d1a_rasscf.F90
     get_d1i_rasscf.F90
     get_umn.F90
     getcnf.F90
     geth2a.F90
     getstepvector.F90
     gugactl.F90
     gugx.F90
     hcsce.F90
     ini_david.F90
     init_lex.F90
     input_mclr.F90
     ipcsfsplit.F90
     iphase.F90
     isgnum.F90
     lex_iter.F90
     load_ci_vec.F90
     load_h_diag.F90
     load_sig_vec.F90
     load_tmp_ci_vec.F90
     load_tmp_sig_vec.F90
     lucia2molcas.F90
     lucia2molcas_free.F90
     matml4.F90
     mkclist.F90
     mkcot.F90
     mkcoup.F90
     mkguga_free.F90
     mkguga_mod.F90
     mkmaw.F90
     mknsm.F90
     mkseg.F90
     mksgnum.F90
     motrac.F90
     newgs.F90
     nrcoup.F90
     ordstr.F90
     output_ras.F90
     page_in.F90
     page_out.F90
     pageno.F90
     phpcsf.F90
     pmat_rasscf.F90
     prttim.F90
     rasdef.F90
     rasscf_global.F90
     rassi_data.F90
     rd1int.F90
     rd2int_rasscf.F90
     recno.F90
     reord2.F90
     restr.F90
     rmvert.F90
     save_ci_vec.F90
     save_h_diag.F90
     save_sig_vec.F90
     save_tmp_ci_vec.F90
     save_tmp_sig_vec.F90
     scavcs.F90
     schmidt.F90
     second_quantization.F90
     segtab.F90
     setsxci.F90
     sgprwf.F90
     sigma1.F90
     sigmadet_cvb.F90
     simplify.F90
     sort_cdet.F90
     spinfo.F90
     splitcas_data.F90
     splitctl.F90
     stepvec.F90
     stepvector_next.F90
     sxci.F90
     term_david.F90
     timers.F90
     transadd.F90
     wadr.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     cands.F90
     caspt2_global.F90
     csfbas.F90
     dmrginfo.F90
     general_data.F90
     gugx.F90
     input_mclr.F90
     rasscf_global.F90
     rassi_data.F90
     spinfo.F90
     timers.F90
     wadr.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
