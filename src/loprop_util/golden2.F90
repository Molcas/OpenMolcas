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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

function Golden2(ax,bx,cx,f,tol_x,tol_f,xmin,q_a,q_b,dipole_a,dipole_b,r_a,r_b)

use Constants, only: One, Three, Five, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Golden2
real(kind=wp), intent(in) :: ax, bx, cx, tol_x, tol_f
real(kind=wp), intent(out) :: xmin
! External function f and its arguments
real(kind=wp), external :: f
real(kind=wp), intent(in) :: q_a, q_b, dipole_a, dipole_b, r_a, r_b
real(kind=wp) :: f2, f3, x1, x2, x3, x4
real(kind=wp), parameter :: Ratio = Half*(Three-sqrt(Five)), RM = One-Ratio
logical(kind=iwp), parameter :: Absolute = .true.

x1 = ax
x4 = cx
if (abs(cx-bx) > abs(bx-ax)) then
  x2 = bx
  x3 = RM*bx+Ratio*cx
else
  x2 = RM*bx+Ratio*ax
  x3 = bx
end if
f2 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x2,Absolute)
f3 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x3,Absolute)
do while ((abs(x4-x1) > tol_x*(abs(x1)+abs(x2))) .and. (abs(f3-f2) > tol_f*(abs(f2)+abs(f3))))
  if (f2 < f3) then
    x4 = x3
    x3 = x2
    f3 = f2
    x2 = RM*x3+Ratio*x1
    f2 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x2,Absolute)
  else
    x1 = x2
    x2 = x3
    f2 = f3
    x3 = RM*x2+Ratio*x4
    f3 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x3,Absolute)
  end if
end do
if (f2 < f3) then
  xmin = x2
  Golden2 = f2
else
  xmin = x3
  Golden2 = f3
end if

end function Golden2
