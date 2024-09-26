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

! Compute the difference between consecutive gradients

subroutine dGrd()

use LnkLst, only: GetNod, iVPTr, LLdGrd, LLGrad, LstPtr, PutVec, SCF_V
use InfSCF, only: iter, Iter_Start, mOV
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, inode, jpgrd
real(kind=wp), allocatable :: Scr(:)

call mma_allocate(Scr,mOV,Label='Scr')

! Loop over all iterations starting at Iter_Start+1

do i=Iter_Start+1,Iter

  !dg(i-1) = g(i)-g(i-1)

  jpgrd = LstPtr(i,LLGrad)
  call GetNod(i-1,LLGrad,inode)
  if (inode == 0) then
    write(u6,*) 'inode == 0'
    call Abend()
  end if
  call iVPtr(Scr,mOV,inode)

  Scr(:) = SCF_V(jpgrd)%A(:)-Scr(:)

  call PutVec(Scr,mOV,i-1,'OVWR',LLdGrd)
end do

call mma_deallocate(Scr)

return

end subroutine dGrd
