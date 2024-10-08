!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

subroutine TrimCMO(CMO,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine trim CMO's from nBas x nBas to nBas x nOrb.             *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp) :: iFrom, iSym, iTo, nData

if (all(nOrb(:nSym-1) == nBas(:nSym-1))) return ! difference in last irrep is irrelevant
if (any(nOrb(:) > nBas(:))) then
  write(u6,*) 'Error in TrimCMO'
  call Abend()
end if
!----------------------------------------------------------------------*
! Transfer orbitals.                                                   *
!----------------------------------------------------------------------*
iFrom = 0
iTo = 0
do iSym=1,nSym
  nData = nBas(iSym)*nOrb(iSym)
  if (iFrom /= iTo) CMO(iTo+1:iTo+nData) = CMO(iFrom+1:iFrom+nData)
  iFrom = iFrom+nBas(iSym)**2
  iTo = iTo+nData
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine TrimCMO
