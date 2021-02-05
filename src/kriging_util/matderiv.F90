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

use kriging_mod

implicit none
#include "stdalloc.fh"
integer diff_Order, d1, d2
real*8 d(d1,d2), m(d1,d2)
! Local variables
integer k
real*8 nr, kr, a, t
real*8, allocatable :: b(:,:), dh(:,:), c(:,:)

call mma_Allocate(b,d1,d2,label='b')
call mma_Allocate(dh,d1,d2,label='dh')
call mma_Allocate(c,d1,d2,label='c')

m = 0

! Analytic or Numerical Derivatives

if (anMd) then
  t = sqrt(2.0d0*pAI+1.0d0)
  dh(:,:) = sqrt(d)
  c(:,:) = (2.0d0*pAI+1.0d0)/(2.0d0*pAI-1.0d0)*exp(-t*dh)
  if (pAI > 3 .or. pAI < 1) then
    write(6,*) 'Analytical Matern derivatives (anamat=.True.)'
    write(6,*) 'is only valid for pAI = 1, 2 and 3 (v = 3/2, 5/2 and 7/2)'
  else
    select case (pAI)
      case (1)
        select case (diff_Order)
          case (1)
            m = -c/2.0d0
          case (2)
            m = c*merge(0.75d0*t/dh,dh,dh /= 0)/3.0d0
          case (3)
            m = -(2.0d0*t-3.0d0*dh)*c
        end select
      case (2)
        select case (diff_Order)
          case (1)
            m = -c*(1.0d0+t*dh)/2.0d0
          case (2)
            m = c*5.0d0/4.0d0
          case (3)
            m = merge(-5.0d0/8.0d0*t/dh,dh,dh /= 0)*c
        end select
      case (3)
        select case (diff_Order)
          case (1)
            m = -c*(1.0d0+t*dh+dh**2)/2.0d0
          case (2)
            m = c*7.0d0*(1.0d0+t*dh)/12.0d0
          case (3)
            m = -c*49.0d0/24.0d0
        end select
    end select
  end if
else
  ! write (6,*) 'Numerical Matern derivatives num',diff_Order
  nr = dble(diff_Order)
  a = gamma(nr+1.0d0)/h**diff_Order
  b = 0.0d0
  do k=0,diff_Order
    kr = dble(k)
    dh(:,:) = d(:,:)+kr*h
    call matern(dh,m,size(dh,1),size(dh,2))
    b(:,:) = b(:,:)+dble((-1)**(k+1))/(gamma(nr-kr+1.0d0)*gamma(kr+1.0d0))*m
  end do
  m = a*b*dble((-1)**(nr+1))
end if

call mma_deAllocate(b)
call mma_deAllocate(dh)
call mma_deAllocate(c)

end subroutine matderiv
