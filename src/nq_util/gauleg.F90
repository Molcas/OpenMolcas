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

subroutine gauleg(x1,x2,xw,n)

use Constants, only: Zero, One, Two, Half, Quart, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: x1, x2
real(kind=wp), intent(out) :: xw(2,n)
integer(kind=iwp) :: i, j, m
real(kind=wp) :: p1, p2, p3, pp, xl, xm, z, z1
real(kind=wp), parameter :: EPS = 3.e-14_wp

m = (n+1)/2
xm = Half*(x2+x1)
xl = Half*(x2-x1)
!write(u6,*) 'm=',m
do i=1,m
  z = cos(Pi*(real(i,kind=wp)-Quart)/(real(n,kind=wp)+Half))
  do
    p1 = One
    p2 = Zero
    do j=1,n
      p3 = p2
      p2 = p1
      p1 = ((Two*real(j,kind=wp)-One)*z*p2-(real(j,kind=wp)-One)*p3)/real(j,kind=wp)
    end do
    pp = real(n,kind=wp)*(z*p1-p2)/(z*z-One)
    z1 = z
    z = z1-p1/pp
    if (abs(z-z1) <= EPS) exit
  end do
  xw(1,i) = xm-xl*z
  xw(1,n+1-i) = xm+xl*z
  xw(2,i) = Two*xl/((One-z*z)*pp*pp)
  xw(2,n+1-i) = xw(2,i)
  if (abs(xw(1,i)) < EPS) xw(1,i) = Zero
  if (abs(xw(1,n+1-i)) < EPS) xw(1,n+1-i) = Zero
  if (abs(xw(2,i)) < EPS) xw(1,i) = Zero
  if (abs(xw(2,n+1-i)) < EPS) xw(1,n+1-i) = Zero
end do

return

end subroutine gauleg
