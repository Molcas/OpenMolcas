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
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************

subroutine MkExSB32(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSBt)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,2) (p,q secondary,active) of  *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

use Cho_Tra, only: nAsh, nSsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), allocatable, intent(out) :: AddSB(:)
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
real(kind=wp), intent(in) :: AddSBt(*)
integer(kind=iwp) :: iIx, LenSB, LxType
logical(kind=iwp) :: SameLx
real(kind=wp), allocatable :: Lx0(:), Ly0(:)

! SubBlock 3 2
LenSB = nSsh(iSymA)*nAsh(iSymB)
call mma_allocate(AddSB,LenSB,Label='AddSB')
if ((iSymA == iSymB) .and. (iSymI == iSymJ) .and. (iI == iJ)) then
  ! SB 3,2 = (SB 2,3)+
  call Trnsps(nSsh(iSymB),nAsh(iSymA),AddSBt,AddSB)
  return
end if

! Build Lx
call mma_allocate(Lx0,nSsh(iSymA)*numV,Label='Lx0')
LxType = 0
iIx = 0
SameLx = .false.
call MkL3(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

! Build Ly
call mma_allocate(Ly0,nAsh(iSymB)*numV,Label='Ly0')
call MkL2(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

! Generate the SubBlock
call DGEMM_('N','T',nAsh(iSymB),nSsh(iSymA),numV,One,Ly0,nAsh(iSymB),Lx0,nSsh(iSymA),Zero,AddSB,nAsh(iSymB))

call mma_deallocate(Ly0)
call mma_deallocate(Lx0)

return

end subroutine MkExSB32
