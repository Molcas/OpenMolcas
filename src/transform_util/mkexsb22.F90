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

subroutine MkExSB22(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,2) (p,q active) of the        *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

use Cho_Tra, only: nAsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), allocatable, intent(out) :: AddSB(:)
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
integer(kind=iwp) :: iIx, LenSB, LxType
logical(kind=iwp) :: SameLx
real(kind=wp), allocatable :: Lx0(:), Ly0(:)

! SubBlock 2 2
LenSB = nAsh(iSymA)*nAsh(iSymB)
call mma_allocate(AddSB,LenSB,Label='AddSB')

! Build Lx
call mma_allocate(Lx0,nAsh(iSymA)*numV,Label='Lx0')
LxType = 0
iIx = 0
SameLx = .false.
call MkL2(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

! Build Ly
call mma_allocate(Ly0,nAsh(iSymB)*numV,Label='Ly0')
if (iSymA == iSymB) SameLx = .true.
call MkL2(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

! Generate the SubBlock
if (.not. SameLx) then
  call DGEMM_('N','T',nAsh(iSymB),nAsh(iSymA),numV,One,Ly0,nAsh(iSymB),Lx0,nAsh(iSymA),Zero,AddSB,nAsh(iSymB))
else
  call DGEMM_('N','T',nAsh(iSymA),nAsh(iSymA),numV,One,Lx0,nAsh(iSymA),Lx0,nAsh(iSymA),Zero,AddSB,nAsh(iSymA))
end if

call mma_deallocate(Ly0)
call mma_deallocate(Lx0)

return

end subroutine MkExSB22
