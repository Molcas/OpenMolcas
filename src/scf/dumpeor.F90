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

subroutine DumpEor(Label,E_or,nSym,nBas,nOrb)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: E_or(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp) :: iFrom(8), iSym, iTo(8), ndata, npDump
real(kind=wp), allocatable :: Dump(:)

!----------------------------------------------------------------------*
! Preliminaries                                                        *
!----------------------------------------------------------------------*
npDump = sum(nBas(:))
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
  Dump(iTo(iSym):iTo(iSym)+ndata-1) = E_or(iFrom(iSym):iFrom(iSym)+ndata-1)
end do
call Put_dArray(Label,Dump,npDump)
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
call mma_deallocate(Dump)

return

end subroutine DumpEor
