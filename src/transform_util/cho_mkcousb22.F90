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

use Cho_Tra

implicit real*8(a-h,o-z)
implicit integer(i-n)
real*8, allocatable :: AddSB(:)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
real*8, allocatable :: Lij(:)

! SubBlock 2 2
LenSB = nAsh(iSymA)*nAsh(iSymB)
call mma_allocate(AddSB,LenSB,Label='AddSB')

! Define Lab
LenAB = LenSB
!-----------------------------------------------------------------------
!if (IfTest) then
!  write(6,*) '     MkCouSB22: TCVD(',iSymA,',',iSymB,')'
!  write(6,'(8F10.6)')TCVX(4,iSymA,iSymB)%A(:,:)
!  call XFlush(6)
!end if
!-----------------------------------------------------------------------

! Build Lij
call mma_allocate(Lij,NumV,Label='Lij')
call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

! Generate the SubBlock
call DGEMM_('N','N',LenAB,1,numV,1.0d0,TCVX(4,iSymA,iSymB)%A,LenAB,Lij,NumV,0.0d0,AddSB,LenSB)

call mma_deallocate(Lij)

return

end subroutine MkCouSB22
