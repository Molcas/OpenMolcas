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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine PotEnergy(A,nMat,iCre,iAnn,energy,grad,Hess,D3,D4,max_term,W,max_Ord,nOsc,nOscOld)
!  Purpose:
!    Calculate matrix elements of potential energy terms.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMat(0:ndim1,ndim2), iAnn(0:ndim1,ndim2), iCre(0:ndim1,ndim2), max_term, max_Ord, nOsc, nOscOld
real(kind=wp), intent(out) :: A(0:max_Ord,0:max_Ord)
real(kind=wp), intent(in) :: energy, grad(nOscOld), Hess(nOscOld,nOscOld), D3(nOscOld,nOscOld,nOscOld), &
                             D4(nOscOld,nOscOld,nOscOld,nOscOld), W(nOscOld,nOsc)
integer(kind=iwp) :: i
real(kind=wp) :: rdx(4)
real(kind=wp), allocatable :: D3_2(:,:,:), D4_2(:,:,:,:), grad_2(:), Hess_2(:,:), Temp(:)

! Zeroth order term.
do i=0,max_Ord
  A(i,i) = Energy
end do
rdx(1) = One
rdx(2) = One
rdx(3) = One
rdx(4) = One
call mma_allocate(Temp,nOscOld**4,label='Temp')

! First order terms.
if (max_term > 0) then
  call mma_allocate(grad_2,nOsc,label='grad_2')
  call DGEMM_('T','N',1,nOsc,nOscOld,One,grad,nOscOld,W,nOscOld,Zero,grad_2,1)
  call Mul1(nMat,A,icre,iann,grad_2,max_Ord,nOsc,rdx)
  call mma_deallocate(grad_2)
end if

! Second order terms.
if (max_term > 1) then
  call mma_allocate(Hess_2,nOsc,nOsc,label='Hess_2')
  call DGEMM_('T','N',nOscOld,nOsc,nOscOld,One,Hess,nOscOld,W,nOscOld,Zero,Temp,nOscOld)
  call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,Hess_2,nOsc)
  call Mul2(nMat,A,icre,iann,Hess_2,max_Ord,nOsc,rdx)
  call mma_deallocate(Hess_2)
end if

! Third order terms.
if (max_term > 2) then
  call mma_allocate(D3_2,nOsc,nOsc,nOsc,label='D3_2')
  call DGEMM_('T','N',nOscOld**2,nOsc,nOscOld,One,D3,nOscOld,W,nOscOld,Zero,Temp,nOscOld**2)
  call DGEMM_('T','N',nOsc*nOscOld,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,D3_2,nOsc*nOscOld)
  call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,D3_2,nOscOld,W,nOscOld,Zero,Temp,nOsc**2)
  call dcopy_(nOsc**3,Temp,1,D3_2,1)
  call Mul3(nMat,A,icre,iann,D3_2,max_Ord,nOsc,rdx)
  call mma_deallocate(D3_2)
end if

! Fourth order terms.
if (max_term > 3) then
  call mma_allocate(D4_2,nOsc,nOsc,nOsc,nOsc,label='D4_2')
  call DGEMM_('T','N',nOscOld**3,nOsc,nOscOld,One,D4,nOscOld,W,nOscOld,Zero,Temp,nOscOld**3)
  call DGEMM_('T','N',nOsc*nOscOld**2,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,D4_2,nOsc*nOscOld**2)
  call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,D4_2,nOscOld,W,nOscOld,Zero,Temp,nOsc**2*nOscOld)
  call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,D4_2,nOsc**3)
  call Mul4(nMat,A,icre,iann,D4_2,max_Ord,nOsc,rdx)
  call mma_deallocate(D4_2)
end if

call mma_deallocate(Temp)

end subroutine PotEnergy
