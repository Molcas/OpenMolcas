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

subroutine XDR_dmatsqrt(a,n)
! Compute the inverse square root of a real symmetric matrix : A**(-1/2)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n,n)
integer(kind=iwp) :: i, info
real(kind=wp) :: dia
real(kind=wp), allocatable :: Tmp(:), Cr(:,:), W(:)

call mma_allocate(Tmp,8*n,label='tmp')
call mma_allocate(Cr,n,n,label='Cr')
call mma_allocate(W,n,label='Eig')
Cr(:,:) = a(:,:)
call dsyev_('V','L',n,Cr,n,W,Tmp,8*n,info)
do i=1,n
  dia = One/sqrt(sqrt(W(i)))
  Cr(:,i) = Cr(:,i)*dia
end do
call dmxma(n,'N','T',Cr,Cr,a,One)
call mma_deallocate(Tmp)
call mma_deallocate(Cr)
call mma_deallocate(W)

return

end subroutine XDR_dmatsqrt
