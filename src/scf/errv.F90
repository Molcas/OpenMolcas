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
! Copyright (C) Martin Schuetz                                         *
!               2017, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ErrV(lvec,ivec,QNRstp,ErrVec)
!***********************************************************************
!                                                                      *
!     computes error vector for DIIS, in 1st order case, this          *
!     is just grad(ivec): in 2nd order case (QNRstp == .true.)         *
!     this is -H(iterso)*grad(ivec)                                    *
!     the pointer to the proper error vector is returned as function   *
!     val.                                                             *
!                                                                      *
!***********************************************************************

use LnkLst, only: LLGrad, GetNod, iVPtr
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer lVec, iVec
real*8 ErrVec(lVec)
logical QNRstp
! local vars
integer inode
real*8, allocatable :: Grad(:)

call GetNod(ivec,LLGrad,inode)
if (inode == 0) then
  ! Hmmm, no entry found in LList, that's strange
  write(6,*) 'ErrV: no entry found in LList!'
  call Abend()
end if

if (QNRstp) then

  ! for qNR step compute delta = - H^{-1}g

  call mma_allocate(Grad,lvec,Label='Grad')
  call iVPtr(Grad,lvec,inode)
# ifdef _DEBUGPRINT_
  write(6,*) 'ErrV, iVec=',iVec
  call NrmClc(Grad,lVec,'ErrV','Grad')
# endif
  call SOrUpV(Grad,lvec,ErrVec,'DISP','BFGS')
# ifdef _DEBUGPRINT_
  call NrmClc(ErrVec,lVec,'ErrV','ErrVec')
# endif
  call mma_deallocate(Grad)

else

  ! Pick up the gradient

  call iVPtr(ErrVec,lvec,inode)

end if

return

end subroutine ErrV
