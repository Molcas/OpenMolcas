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

subroutine TrimEor(E_or,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine trim orbital energy vectors.                            *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: E_or(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp) :: iFrom, iSym, iTo, nData

if (all(nOrb(:nSym-1) == nBas(:nSym-1))) return ! difference in last irrep is irrelevant
if (any(nOrb(:) > nBas(:))) then
  write(u6,*) 'Error in TrimEor'
  call Abend()
end if
!----------------------------------------------------------------------*
! Transfer orbital energies.                                           *
!----------------------------------------------------------------------*
iFrom = 0
iTo = 0
do iSym=nSym,1,-1
  nData = nOrb(iSym)
  if (iFrom /= iTo) E_or(iTo+1:iTo+nData) = E_or(iFrom+1:iFrom+nData)
  iFrom = iFrom+nBas(iSym)
  iTo = iTo+nData
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine TrimEor
