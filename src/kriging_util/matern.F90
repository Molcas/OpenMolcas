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

use kriging_mod

implicit none
#include "stdalloc.fh"
integer d1, d2, i
real*8 a, d, dh(d1,d2), m(d1,d2)
real*8, allocatable :: d0(:,:)

call mma_Allocate(d0,d1,d2,label='d0')

d0(:,:) = sqrt(dh)
select case (pAI)
  case (0) ! v = 1/2
    m = exp(-d0)
  case (1) ! v = 3/2
    m = exp(-sqrt(3.0d0)*d0)*(sqrt(3.0d0)*d0+1.0d0)
  case (2) ! v = 5/2
    m = exp(-sqrt(5.0d0)*d0)*(5.0d0/3.0d0*d0**2+sqrt(5.0d0)*d0+1.0d0)
  case (3) ! v = 7/2
    m = exp(-sqrt(7.0d0)*d0)*(7.0d0/15.0d0*sqrt(7.0d0)*d0**3+14.0d0/5.0d0*d0**2+sqrt(7.0d0)*d0+1.0d0)
  case default
    ! For this expresion you can check https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
    ! and equations (11) and (12) on ref.
    a = gamma(pAI+1.0d0)/gamma(2.0d0*pAI+1.0d0)
    m = 0.0d0
    do i=0,pAI
      d = dble(i)
      m = m+(gamma(pAI+1.0d0+d)/(gamma(d+1.0d0)*gamma(pAI+1.0d0-d)))*(2.0d0*sqrt(2.0d0*pAI+1.0d0)*d0)**(pAI-i)
    end do
    m = a*m*exp(-sqrt(2.0d0*pAI+1.0d0)*d0)
end select

call mma_deallocate(d0)

end subroutine matern
