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

subroutine MkExSB33(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           Jannuary February 2005                                     *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,3) (p,q secondary) of the     *
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

! SubBlock 3 3
LenSB = nSsh(iSymA)*nSsh(iSymB)
call mma_allocate(AddSB,LenSB,Label='AddSB')

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
  call DGEMM_('N','T',nSsh(iSymB),nSsh(iSymA),numV,1.0d0,Ly0,nSsh(iSymB),Lx0,nSsh(iSymA),0.0d0,AddSB,nSsh(iSymB))
else
  call DGEMM_('N','T',nSsh(iSymA),nSsh(iSymA),numV,1.0d0,Lx0,nSsh(iSymA),Lx0,nSsh(iSymA),0.0d0,AddSB,nSsh(iSymA))
end if

call mma_deallocate(Ly0)
call mma_deallocate(Lx0)

return

end subroutine MkExSB33
