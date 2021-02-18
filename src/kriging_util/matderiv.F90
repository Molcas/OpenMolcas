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

subroutine matderiv(diff_Order,d,m,d1,d2)

use kriging_mod, only: anMd, h, pAI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Five, Seven, Eight, Twelve
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: diff_Order, d1, d2
real(kind=wp), intent(in) :: d(d1,d2)
real(kind=wp), intent(out) :: m(d1,d2)
! Local variables
integer(kind=iwp) :: k
real(kind=wp) :: nr, kr, a, t
real(kind=wp), allocatable :: b(:,:), dh(:,:), c(:,:)

call mma_Allocate(b,d1,d2,label='b')
call mma_Allocate(dh,d1,d2,label='dh')
call mma_Allocate(c,d1,d2,label='c')

m(:,:) = 0

! Analytic or Numerical Derivatives

if (anMd) then
  t = sqrt(Two*pAI+One)
  dh(:,:) = sqrt(d)
  c(:,:) = (Two*pAI+One)/(Two*pAI-One)*exp(-t*dh)
  if (pAI > 3 .or. pAI < 1) then
    write(u6,*) 'Analytical Matern derivatives (anamat=.True.)'
    write(u6,*) 'is only valid for pAI = 1, 2 and 3 (v = 3/2, 5/2 and 7/2)'
  else
    select case (pAI)
      case (1)
        select case (diff_Order)
          case (1)
            m(:,:) = -c(:,:)/Two
          case (2)
            m(:,:) = c(:,:)*merge(0.75_wp*t/dh,dh,dh /= 0)/Three
          case (3)
            m(:,:) = -(Two*t-Three*dh)*c(:,:)
        end select
      case (2)
        select case (diff_Order)
          case (1)
            m(:,:) = -c(:,:)*(One+t*dh)/Two
          case (2)
            m(:,:) = c(:,:)*Five/Four
          case (3)
            m(:,:) = merge(-Five/Eight*t/dh,dh,dh /= 0)*c(:,:)
        end select
      case (3)
        select case (diff_Order)
          case (1)
            m(:,:) = -c(:,:)*(One+t*dh+dh**2)/Two
          case (2)
            m(:,:) = c(:,:)*Seven*(One+t*dh)/Twelve
          case (3)
            m(:,:) = -c(:,:)*49.0_wp/24.0_wp
        end select
    end select
  end if
else
  ! write (u6,*) 'Numerical Matern derivatives num',diff_Order
  nr = real(diff_Order,kind=wp)
  a = gamma(nr+One)/h**diff_Order
  b = Zero
  do k=0,diff_Order
    kr = real(k,kind=wp)
    dh(:,:) = d(:,:)+kr*h
    call matern(dh,m,size(dh,1),size(dh,2))
    b(:,:) = b(:,:)+One**(k+1)/(gamma(nr-kr+One)*gamma(kr+One))*m(:,:)
  end do
  m(:,:) = a*b*(-One)**(nr+1)
end if

call mma_deAllocate(b)
call mma_deAllocate(dh)
call mma_deAllocate(c)

end subroutine matderiv
