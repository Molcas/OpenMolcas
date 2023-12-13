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
     analysis_domain.F90
     anasize_localisation.F90
     basfun_atom_.F90
     basfun_atom.F90
     bestmatch.F90
     boys.F90
     boys_iter.F90
     checkdomain.F90
     chk_unitary.F90
     choloc.F90
     choloc_p.F90
     choloc_xp.F90
     cho_ov_loc.F90
     computefuncb2.F90
     computefuncer.F90
     computefunc.F90
     definedomain.F90
     definepairdomain.F90
     delete_ghosts.F90
     diag_localisation.F90
     domain_histogram.F90
     edmistonruedenberg.F90
     edmistonruedenberg_iter.F90
     erchk_localisation.F90
     erloc.F90
     evalerfun.F90
     fwt_haar.F90
     genbmp.F90
     genbmp_loc.F90
     generateb.F90
     generatep.F90
     getat_localisation.F90
     get_can_lorb.F90
     getdens_localisation.F90
     getgrad_boys.F90
     getgrad_er.F90
     getgrad_pm.F90
     get_orb_select.F90
     getovlp_localisation.F90
     getrawpaos.F90
     getsh_localisation.F90
     getu_er.F90
     getumat_localisation.F90
     get_vir_select.F90
     inv_fwt_haar.F90
     irnge.F90
     log2.F90
     makedomaincomplete.F90
     ortho_orb.F90
     orthopao_localisation.F90
     paoloc.F90
     pipekmezey.F90
     pipekmezey_iter.F90
     pmloc.F90
     rotateorbb.F90
     rotateorb_er.F90
     rotateorb.F90
     rot_st.F90
     sort_localisation_1.F90
     sort_localisation.F90
     sq2tri.F90
     sqrtmt.F90
     tri2rec.F90
     unfold.F90
     updateb.F90
     updatep.F90
     wavelet_transform.F90
     xdiaxt.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
