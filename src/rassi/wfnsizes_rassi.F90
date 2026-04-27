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

subroutine wfnsizes_rassi()

use rasdef, only: NRS1, NRS1T, NRS2, NRS2T, NRS3, NRS3T
use rassi_aux, only: nasht_save
use Symmetry_Info, only: nIrrep
use rassi_data, only: NAES, NASH, NASHT, NBASF, NBST, NISH, NISHT, NOSH, NOSHT, NSSH, NSSHT
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: isym

! The structure of the orbital space:
NISHT = sum(NISH(1:nIrrep))
NASHT = 0
NRS1T = sum(NRS1(1:nIrrep))
NRS2T = sum(NRS2(1:nIrrep))
NRS3T = sum(NRS3(1:nIrrep))
NOSHT = sum(NOSH(1:nIrrep))
NSSHT = sum(NSSH(1:nIrrep))
NBST = sum(NBASF(1:nIrrep))
do ISYM=1,nIrrep
  NAES(ISYM) = NASHT
  NASHT = NASHT+NASH(ISYM)
end do
nasht_save = nasht

end subroutine wfnsizes_rassi
