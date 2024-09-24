!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! Compute the dx parameters value

subroutine dX()

use LnkLst, only: SCF_V, LLx, LLDelt, LstPtr, GetNod, iVPtr, PutVec
use InfSCF, only: Iter, Iter_Start, mOV
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer jpgrd, inode, i
real*8, dimension(:), allocatable :: Scr

call mma_allocate(Scr,mOV,Label='Scr')

! Loop over all iterations starting at Iter_Start

do i=Iter_Start,Iter-1

  !dX(i) = X(i+1)-X(i)

  jpgrd = LstPtr(i+1,LLx)   ! Pointer to X(i+1)

  call GetNod(i,LLx,inode) ! X(i)
  if (inode == 0) then
    write(6,*) 'inode == 0'
    call Abend()
  end if
  call iVPtr(Scr,mOV,inode)

  Scr(:) = SCF_V(jpgrd)%A(:)-Scr(:)

  call PutVec(Scr,mOV,i,'OVWR',LLDelt)
end do

call mma_deallocate(Scr)

return

end subroutine dX
