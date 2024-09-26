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
! This routine dumps CMO's onto runfile, expands them if necessary.    *
!                                                                      *
!***********************************************************************

subroutine DumpCMO(Label,CMO,nSym,nBas,nOrb)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
character(len=*) :: Label
real(kind=wp) :: CMO(*)
integer(kind=iwp) :: nSym, nBas(*), nOrb(*)
integer(kind=iwp) :: iFrom(8), iSym, iTo(8), ndata, npDump
real(kind=wp), allocatable :: Dump(:)

!----------------------------------------------------------------------*
! Preliminaries                                                        *
!----------------------------------------------------------------------*
npDump = 0
do iSym=1,nSym
  npDump = npDump+nBas(iSym)*nBas(iSym)
end do
call mma_allocate(Dump,npDump,Label='Dump')
!----------------------------------------------------------------------*
! Dump orbitals                                                        *
!----------------------------------------------------------------------*
iFrom(1) = 1
iTo(1) = 1
do iSym=1,nSym-1
  iFrom(iSym+1) = iFrom(iSym)+nBas(iSym)*nOrb(iSym)
  iTo(iSym+1) = iTo(iSym)+nBas(iSym)*nBas(iSym)
end do
do iSym=nSym,1,-1
  ndata = nBas(iSym)*nOrb(iSym)
  call DCopy_(ndata,CMO(iFrom(iSym)),1,Dump(iTo(iSym)),1)
end do
call Put_dArray(Label,Dump,npDump)
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
call mma_deallocate(Dump)

return

end subroutine DumpCMO
