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

subroutine MkCouSB22(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,2) (p,q active) of the        *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

use Cho_Tra, only: nAsh, TCVX
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), allocatable, intent(out) :: AddSB(:)
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
integer(kind=iwp) :: LenAB, LenSB
real(kind=wp), allocatable :: Lij(:)

! SubBlock 2 2
LenSB = nAsh(iSymA)*nAsh(iSymB)
call mma_allocate(AddSB,LenSB,Label='AddSB')

! Define Lab
LenAB = LenSB
!-----------------------------------------------------------------------
!if (IfTest) then
!  write(u6,*) '     MkCouSB22: TCVD(',iSymA,',',iSymB,')'
!  write(u6,'(8F10.6)') TCVX(4,iSymA,iSymB)%A(:,:)
!  call XFlush(u6)
!end if
!-----------------------------------------------------------------------

! Build Lij
call mma_allocate(Lij,NumV,Label='Lij')
call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

! Generate the SubBlock
call DGEMM_('N','N',LenAB,1,numV,One,TCVX(4,iSymA,iSymB)%A,LenAB,Lij,NumV,Zero,AddSB,LenSB)

call mma_deallocate(Lij)

return

end subroutine MkCouSB22
