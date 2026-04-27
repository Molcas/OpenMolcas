!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module RASSIWfn

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: wfn_cmo, wfn_cmo_or, wfn_detcoeff, wfn_detcoeff_or, wfn_detocc, wfn_detocc_or, wfn_overlap, wfn_sfs_amfi, &
                     wfn_sfs_angmom, wfn_sfs_coef, wfn_sfs_edipmom, wfn_sfs_energy, wfn_sfs_tdm, wfn_sfs_tm, wfn_sfs_tsdm, &
                     wfn_sfs_wetdm, wfn_sos_angmomi, wfn_sos_angmomr, wfn_sos_coefi, wfn_sos_coefr, wfn_sos_dys, wfn_sos_edipmomi, &
                     wfn_sos_edipmomr, wfn_sos_energy, wfn_sos_hsoi, wfn_sos_hsor, wfn_sos_spini, wfn_sos_spinr, wfn_sos_tm, &
                     wfn_sos_vsoi, wfn_sos_vsor

public :: wfn_cmo, wfn_cmo_or, wfn_detcoeff, wfn_detcoeff_or, wfn_detocc, wfn_detocc_or, wfn_overlap, wfn_sfs_amfi, &
          wfn_sfs_angmom, wfn_sfs_coef, wfn_sfs_edipmom, wfn_sfs_energy, wfn_sfs_tdm, wfn_sfs_tm, wfn_sfs_tsdm, wfn_sfs_wetdm, &
          wfn_sos_angmomi, wfn_sos_angmomr, wfn_sos_coefi, wfn_sos_coefr, wfn_sos_dys, wfn_sos_edipmomi, wfn_sos_edipmomr, &
          wfn_sos_energy, wfn_sos_hsoi, wfn_sos_hsor, wfn_sos_spini, wfn_sos_spinr, wfn_sos_tm, wfn_sos_vsoi, wfn_sos_vsor

end module RASSIWfn
