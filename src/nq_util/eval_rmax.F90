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

use Constants, only: One, Ten, Half, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Eval_RMax
real(kind=wp), intent(in) :: alpha, R_L
integer(kind=iwp), intent(in) :: m
real(kind=wp) :: Gmma, x, x_new
integer(kind=iwp) :: i

!write(u6,*) 'alpha,m,R_L=',alpha,m,R_L

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute r_k_H as a function of m and R_L
!
! Eq(19) R. Lindh, P.-A. Malmqvist, L. Gagliardi,
! TCA, 106:178-187 (2001)

if (mod(m+3,2) == 0) then
  Gmma = One
  do i=2,(m+3)/2,1
    Gmma = Gmma*real(i-1,kind=wp)
  end do
else
  Gmma = sqrt(Pi)
  do i=5,m+3,2
    Gmma = Gmma*real(i-1,kind=wp)*Half
  end do
end if

!x = Alpha*(r_k_H)**2

x = Ten ! Start value
do
  x_new = log((Gmma/R_L)*x**(Half*(real(m,kind=wp)+One)))
  !write(u6,*) 'x,x_new=',x,x_new
  if (abs(x-x_new) <= 1.0e-8_wp) exit
  x = x_new
end do

Eval_RMax = sqrt(x/alpha)
!write(u6,*) 'Eval_RMax=',Eval_RMax
!                                                                      *
!***********************************************************************
!                                                                      *

return

end function Eval_RMax
