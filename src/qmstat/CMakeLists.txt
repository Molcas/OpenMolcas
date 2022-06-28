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
     abboth.F90
     abnone.F90
     abone.F90
     allenginsberg.F90
     analyze_q.F90
     aointegrate.F90
     avermep.F90
     boostrep.F90
     bornmayerbk.F90
     chk_oneham.F90
     ciselector.F90
     clasclas.F90
     contractovl.F90
     contrasbas.F90
     cooout.F90
     coult0_1.F90
     coult0_2.F90
     coult0_4.F90
     coult0_5.F90
     coultn_1.F90
     coultn_2.F90
     coultn_4.F90
     coultn_5.F90
     dcorrcorr.F90
     densi_mo.F90
     densist.F90
     dispenergy.F90
     editstart.F90
     eqras.F90
     eqscf.F90
     exnemo.F90
     expectus.F90
     exras.F90
     exscf.F90
     extract.F90
     extractr.F90
     extracts.F90
     fetchtdm.F90
     ffactor.F90
     geogen.F90
     get8.F90
     get9.F90
     get_centers.F90
     get_qmstat_input.F90
     get_slater.F90
     haveweconv.F90
     hel.F90
     helstate.F90
     idubfac.F90
     isitvalid.F90
     m2trans.F90
     mandatoryinp.F90
     mbpt2corr.F90
     moldendump.F90
     momentmod.F90
     moreduce.F90
     multinew.F90
     niceoutput.F90
     noverp_q.F90
     nypart.F90
     offatom.F90
     oneoverr.F90
     oneoverr_sl.F90
     orbrot2.F90
     overlq.F90
     pararoot.F90
     placeit9.F90
     placeit.F90
     planevectors.F90
     polink.F90
     polins.F90
     polprep.F90
     polras.F90
     polscf.F90
     polsolv.F90
     put8.F90
     put9.F90
     qfread.F90
     qmposition.F90
     qmstat.F90
     qmstat_global.F90
     qmstat_init.F90
     qmstat_procedures.F90
     rash0.F90
     rasrastrans.F90
     rassihandm.F90
     reaind.F90
     revolution.F90
     rotation_qmstat.F90
     scfh0.F90
     scfhandm.F90
     singp.F90
     sl_grad.F90
     spherical.F90
     sqtotri_q.F90
     statemmeao.F90
     statemme.F90
     statemmemo.F90
     tdmtrans.F90
     tkp.F90
     transrot.F90
     wrrdsim.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
