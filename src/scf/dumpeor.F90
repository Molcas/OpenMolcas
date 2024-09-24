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
!***********************************************************************
!                                                                      *
! This routine dumps orbital energies on the runfile, expands them if  *
! necessary.                                                           *
!                                                                      *
!***********************************************************************

subroutine DumpEor(Label,Eor,nSym,nBas,nOrb)

use stdalloc, only: mma_allocate, mma_deallocate

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character*(*) Label
real*8 Eor(*)
integer nSym
integer nBas(*)
integer nOrb(*)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer npDump
integer iFrom(8)
integer iTo(8)
integer iSym
integer ndata
real*8, dimension(:), allocatable :: Dump

!----------------------------------------------------------------------*
! Preliminaries                                                        *
!----------------------------------------------------------------------*
npDump = 0
do iSym=1,nSym
  npDump = npDump+nBas(iSym)
end do
call mma_allocate(Dump,npDump,Label='DumpOE')
!----------------------------------------------------------------------*
! Dump orbital energies                                                *
!----------------------------------------------------------------------*
iFrom(1) = 1
iTo(1) = 1
do iSym=1,nSym-1
  iFrom(iSym+1) = iFrom(iSym)+nOrb(iSym)
  iTo(iSym+1) = iTo(iSym)+nBas(iSym)
end do
do iSym=nSym,1,-1
  ndata = nOrb(iSym)
  call dcopy_(ndata,Eor(iFrom(iSym)),1,Dump(iTo(iSym)),1)
end do
call Put_dArray(Label,Dump,npDump)
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
call mma_deallocate(Dump)

return

end subroutine DumpEor
