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

subroutine MkExSB11(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(1,1) (p,q inactive) of the      *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

use Cho_Tra

implicit real*8(a-h,o-z)
implicit integer(i-n)
real*8, allocatable :: AddSB(:)
integer iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
logical SameLx

real*8, allocatable :: Lx0(:), Ly0(:)

! SubBlock 1 1
LenSB = nIsh(iSymA)*nIsh(iSymB)
call mma_allocate(AddSB,LenSB,Label='LenSB')

! Build Lx
call mma_allocate(Lx0,nIsh(iSymA)*numV,Label='Lx0')
LxType = 0
iIx = 0
SameLx = .false.
call MkL1(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

! Build Ly
call mma_allocate(Ly0,nIsh(iSymB)*numV,Label='Ly0')
if (iSymA == iSymB) SameLx = .true.
call MkL1(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

! Generate the SubBlock (Ly*Lx)
if (.not. SameLx) then
  call DGEMM_('N','T',nIsh(iSymB),nIsh(iSymA),numV,1.0d0,Ly0,nIsh(iSymB),Lx0,nIsh(iSymA),0.0d0,AddSB,nIsh(iSymB))
else
  call DGEMM_('N','T',nIsh(iSymA),nIsh(iSymA),numV,1.0d0,Lx0,nIsh(iSymA),Lx0,nIsh(iSymA),0.0d0,AddSB,nIsh(iSymA))
end if

call mma_deallocate(Ly0)
call mma_deallocate(Lx0)

return

end subroutine MkExSB11
