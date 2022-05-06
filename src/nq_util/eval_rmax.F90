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

function Eval_RMax(alpha,m,R_L)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Eval_RMax

!write(6,*) 'alpha,m,R_L=',alpha,m,R_L

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute r_k_H as a function of m and R_L
!
! Eq(19) R. Lindh, P.-A. Malmqvist, L. Gagliardi,
! TCA, 106:178-187 (2001)

if (mod(m+3,2) == 0) then
  Gamma = One
  do i=2,(m+3)/2,1
    Gamma = Gamma*dble(i-1)
  end do
else
  Gamma = sqrt(Pi)
  do i=5,m+3,2
    Gamma = Gamma*dble(i-1)/Two
  end do
end if

!x = Alpha*(r_k_H)**2

x = 10.0d0 ! Start value
do
  x_new = log((Gamma/R_L)*x**(Half*(dble(m)+One)))
  !write(6,*) 'x,x_new=',x,x_new
  if (abs(x-x_new) <= 1.0D-8) exit
  x = x_new
end do

Eval_RMax = sqrt(x/alpha)
!write(6,*) 'Eval_RMax=',Eval_RMax
!                                                                      *
!***********************************************************************
!                                                                      *

return

end function Eval_RMax
