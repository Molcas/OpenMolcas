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

subroutine TrimCMO(CMO1,CMO2,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine trim CMO's from nBas x nBas to nBas x nOrb.             *
!                                                                      *
!***********************************************************************

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
real*8 CMO1(*)
real*8 CMO2(*)
integer nSym
integer nBas(*)
integer nOrb(*)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer iFrom(8)
integer iTo(8)
integer iSym
integer ndata, i

!----------------------------------------------------------------------*
! Transfer orbitals.                                                   *
!----------------------------------------------------------------------*
iFrom(1) = 1
iTo(1) = 1
do iSym=1,nSym-1
  iFrom(iSym+1) = iFrom(iSym)+nBas(iSym)*nBas(iSym)
  iTo(iSym+1) = iTo(iSym)+nBas(iSym)*nOrb(iSym)
  if (iTo(iSym+1) > iFrom(iSym+1)) then
    write(6,*) 'Error in TrimCMO'
    call Abend()
  end if
end do
do iSym=1,nSym
  ndata = nBas(iSym)*nOrb(iSym)

  ! Note that CMO1 and CMO2 might overlap. Hence, we cannot use
  ! an ordinary call to DCopy!

  if (iFrom(iSym) /= iTo(iSym)) then
    do i=0,nData-1
      CMO2(iTo(iSym)+i) = CMO1(iFrom(iSym)+i)
    end do
  end if
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine TrimCMO
