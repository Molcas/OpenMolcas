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

subroutine ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddEx,LenE)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           January-February 2005                                      *
! Modified for Cholesky-MP2 May 2005                                   *
!***********************************************************************

use Cho_Tra, only: nSsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV, LenE
real(kind=wp), intent(inout) :: AddEx(LenE)
integer(kind=iwp) :: iIx, LxType
logical(kind=iwp) :: SameLx
real(kind=wp), allocatable :: Lx0(:), Ly0(:)

! SubBlock 3 3

! Build Lx
call mma_allocate(Lx0,nSsh(iSymA)*numV,Label='Lx0')
LxType = 0
iIx = 0
SameLx = .false.
call MkL3(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

! Build Ly
call mma_allocate(Ly0,nSsh(iSymB)*numV,Label='Ly0')
if (iSymA == iSymB) SameLx = .true.
call MkL3(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

! Generate the SubBlock
if (.not. SameLx) then
  call DGEMM_('N','T',nSsh(iSymB),nSsh(iSymA),numV,One,Ly0,nSsh(iSymB),Lx0,nSsh(iSymA),One,AddEx,nSsh(iSymB))
else
  call DGEMM_('N','T',nSsh(iSymA),nSsh(iSymA),numV,One,Lx0,nSsh(iSymA),Lx0,nSsh(iSymA),One,AddEx,nSsh(iSymA))
end if

call mma_deallocate(Ly0)
call mma_deallocate(Lx0)

return

end subroutine ChoMP2_GenE
