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

real*8 function D1Mach(I)
! DOUBLE-PRECISION MACHINE CONSTANTS
!
! D1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!
! D1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.

integer I
real*8 eps, xmin, xmax
! PAM 2008 These declarations should be used with machar
!real*8 epsneg
!integer ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp

! PAM 2008 Do not use machar any more. Try hardcoded values:
!call machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,maxexp,eps,epsneg,xmin,xmax)
EPS = 1.0D-15
XMIN = 1.0D-300
XMAX = 1.0D+300

if (I == 1) then
  D1Mach = xmin
else if (I == 2) then
  D1Mach = xmax
else if (I == 4) then
  D1Mach = eps
else
  D1Mach = -1d0
end if

return

end function D1Mach
