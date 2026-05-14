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

subroutine PadEor(E_or,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine pads orbital energy vectors.                            *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: E_or(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp) :: iFrom(8), iPtr, iSym, iTo(8), ndata

!----------------------------------------------------------------------*
! Transfer orbital energies.                                           *
!----------------------------------------------------------------------*
iFrom(1) = nOrb(1)
iTo(1) = nOrb(1)
do iSym=1,nSym-1
  iFrom(iSym+1) = iFrom(iSym)+nOrb(iSym+1)
  iTo(iSym+1) = iTo(iSym)+nOrb(iSym+1)+nBas(iSym)-nOrb(iSym)
end do
do iSym=nSym,1,-1
  ndata = nOrb(iSym)
  E_or(iTo(iSym)-ndata+1:iTo(iSym)) = E_or(iFrom(iSym)-ndata+1:iFrom(iSym))
  if (nBas(iSym) > nOrb(iSym)) then
    ndata = nBas(iSym)-nOrb(iSym)
    iPtr = iTo(iSym)+1
    E_or(iPtr:iPtr+ndata-1) = Zero
  end if
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine PadEor
