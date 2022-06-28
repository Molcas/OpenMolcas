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
     at34r.F90
     auxc.F90
     bssint.F90
     cplab.F90
     cplabr.F90
     cp_prop_int.F90
     crelop.F90
     dcopiv.F90
     der.F90
     diag_dkh.F90
     diagr.F90
     dkre1r.F90
     even2.F90
     even2r.F90
     even3.F90
     even3r.F90
     extc.F90
     gam.F90
     jacob_rel.F90
     matinv.F90
     minvd.F90
     nprop_int.F90
     phi.F90
     piknam.F90
     prmat.F90
     prsq.F90
     rdname.F90
     relop.F90
     repmat.F90
     scfcli2.F90
     scfcli4.F90
     scfcli.F90
     simplm.F90
     sog.F90
     sogr.F90
     theta.F90
     trsm_dkh.F90
     trsmr2.F90
     trsmr.F90
     trsmrn.F90
     trsmt.F90
     trsmtr.F90
     vpbmbpv.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
