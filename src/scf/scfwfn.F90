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

module SCFWfn

integer wfn_fileid
integer wfn_energy
integer wfn_mocoef, wfn_occnum, wfn_dens, wfn_spindens
integer wfn_supsym, wfn_basfns, wfn_ovlmat, wfn_orbene
integer wfn_tpidx
integer wfn_mocoef_a, wfn_occnum_a, wfn_orbene_a, wfn_tpidx_a
integer wfn_mocoef_b, wfn_occnum_b, wfn_orbene_b, wfn_tpidx_b

end module SCFWfn

