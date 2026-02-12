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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine CMSFitTrigonometric(x,y)

implicit none
real*8, dimension(4) :: x, y
real*8 s12, s23, c12, c23, d12, d23, k, a, b, c, phi, psi1, psi2, val1, val2
real*8 atan1

s12 = sin(4.0d0*x(1))-sin(4.0d0*x(2))
s23 = sin(4.0d0*x(2))-sin(4.0d0*x(3))
c12 = cos(4.0d0*x(1))-cos(4.0d0*x(2))
c23 = cos(4.0d0*x(2))-cos(4.0d0*x(3))
d12 = y(1)-y(2)
d23 = y(2)-y(3)
k = s12/s23
c = (d12-k*d23)/(c12-k*c23)
b = (d12-c*c12)/s12
a = y(1)-b*sin(4.0d0*x(1))-c*cos(4.0d0*x(1))
phi = atan(b/c)
atan1 = atan(1.0d0)
psi1 = phi/4.0d0
if (psi1 > atan1) then
  psi2 = psi1-atan1
else
  psi2 = psi1+atan1
end if
val1 = b*sin(4.0d0*psi1)+c*cos(4.0d0*psi1)
val2 = b*sin(4.0d0*psi2)+c*cos(4.0d0*psi2)
if (val1 > val2) then
  x(4) = psi1
  !y(4) = val1
else
  x(4) = psi2
  !y(4) = val2
end if
y(4) = a+sqrt(b**2+c**2)
!write(6,*) a,b,c,x(4),y(4)

end subroutine CMSFitTrigonometric
