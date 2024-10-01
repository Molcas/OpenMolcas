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

subroutine PadEor(Eor1,Eor2,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine pads orbital energy vectors.                            *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Eor1(*), Eor2(*)
integer(kind=iwp) :: nSym, nBas(nSym), nOrb(nSym)
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
  Eor2(iTo(iSym)-ndata+1:iTo(iSym)) = Eor1(iFrom(iSym)-ndata+1:iFrom(iSym))
  if (nBas(iSym) > nOrb(iSym)) then
    ndata = nBas(iSym)-nOrb(iSym)
    iPtr = iTo(iSym)+1
    Eor2(iPtr:iPtr+ndata-1) = Zero
  end if
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine PadEor
