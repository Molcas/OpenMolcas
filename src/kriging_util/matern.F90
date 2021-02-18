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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine matern(dh,m,d1,d2)

use kriging_mod, only: pAI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Five, Seven
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d1, d2
real(kind=wp), intent(in) :: dh(d1,d2)
real(kind=wp), intent(out) :: m(d1,d2)
integer(kind=iwp) :: i
real(kind=wp) :: a, d
real(kind=wp), allocatable :: d0(:,:)

call mma_Allocate(d0,d1,d2,label='d0')

d0(:,:) = sqrt(dh)
select case (pAI)
  case (0) ! v = 1/2
    m(:,:) = exp(-d0(:,:))
  case (1) ! v = 3/2
    m(:,:) = exp(-sqrt(Three)*d0(:,:))*(sqrt(Three)*d0(:,:)+One)
  case (2) ! v = 5/2
    m(:,:) = exp(-sqrt(Five)*d0(:,:))*(Five/Three*d0(:,:)**2+sqrt(Five)*d0(:,:)+One)
  case (3) ! v = 7/2
    m(:,:) = exp(-sqrt(Seven)*d0(:,:))*(Seven/15.0_wp*sqrt(Seven)*d0(:,:)**3+14.0_wp/Five*d0(:,:)**2+sqrt(Seven)*d0(:,:)+One)
  case default
    ! For this expresion you can check https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
    ! and equations (11) and (12) on ref.
    a = gamma(pAI+One)/gamma(Two*pAI+One)
    m(:,:) = Zero
    do i=0,pAI
      d = real(i,kind=wp)
      m(:,:) = m(:,:)+(gamma(pAI+One+d)/(gamma(d+One)*gamma(pAI+One-d)))*(Two*sqrt(Two*pAI+One)*d0(:,:))**(pAI-i)
    end do
    m(:,:) = a*m(:,:)*exp(-sqrt(Two*pAI+One)*d0(:,:))
end select

call mma_deallocate(d0)

end subroutine matern
