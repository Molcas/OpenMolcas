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
     alloc.f
     check_invec.f
     chksplit.F90
     cidia_ci_util.F90
     ciovlp.F90
     ciselect.F90
     citrans.F90
     citrans_free.F90
     citrans_init.F90
     cisx.f
     cndet.F90
     cnfstr.F90
     cnhcn.F90
     csdiag_ci_util.F90
     csdtvc.F90
     csfbas.F90
     cstart_ci_util.F90
     davctl.F90
     davctl_mod.F90
     david5.F90
     dblock.F90
     detstr2.F90
     dgmm2_molcas.F90
     diagord.F90
     dihdj_molcas.F90
     dmatdmat.f
     expcsf.F90
     explh2.F90
     faroald.F90
     faroald_free.F90
     faroald_init.F90
     fndmnx.F90
     fockoc.f
     gatvcs.F90
     general_data.F90
     get_cm_.F90
     get_cm.F90
     get_d1a.f
     get_d1i.f
     getcnf_lucia.F90
     geth2a.F90
     getstepvector.F90
     get_umn.F90
     gugactl.F90
     gugx.F90
     hcsce.F90
     ini_david.F90
     init_lex.F90
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
     mkdaw.F90
     mkdrt.F90
     mkdrt0.F90
     mkguga.F90
     mkmid.F90
     mknsm.f
     mkraw.F90
     mksgnum.f
     motrac.F90
     newgs.F90
     ordstr.F90
     page_in.F90
     pageno.F90
     page_out.F90
     phpcsf.F90
     pmat_rasscf.f
     prdrt.F90
     prttim.f
     rd1int.f
     rd2int_rasscf.f
     recno.F90
     restr.F90
     reord2.F90
     save_ci_vec.F90
     save_h_diag.F90
     save_sig_vec.F90
     save_tmp_ci_vec.F90
     save_tmp_sig_vec.F90
     scavcs.F90
     schmidt.F90
     second_quantization.F90
     setsxci.f
     sigmadet_cvb.F90
     simplify.F90
     sort_cdet.F90
     splitctl.F90
     splitcas_data.F90
     stepvec.F90
     stepvector_next.F90
     sxci.F90
     term_david.F90
     transadd.F90
     wadr.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     csfbas.F90
     general_data.F90
     gugx.F90
     splitcas_data.F90
     sxci.F90
     wadr.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
