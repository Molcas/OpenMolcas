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

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
real*8 Eor1(*)
real*8 Eor2(*)
integer nSym
integer nBas(*)
integer nOrb(*)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer iFrom(8)
integer iTo(8)
integer iPtr
integer iSym
integer ndata
integer i

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
  do i=1,ndata
    Eor2(iTo(iSym)+1-i) = Eor1(iFrom(iSym)+1-i)
  end do
  if (nBas(iSym) > nOrb(iSym)) then
    ndata = nBas(iSym)-nOrb(iSym)
    iPtr = iTo(iSym)+1
    call dCopy_(ndata,[Zero],0,Eor2(iPtr),1)
  end if
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine PadEor
