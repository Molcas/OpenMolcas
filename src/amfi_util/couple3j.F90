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

function couple3J(l1,l2,l3,m1,m2,m3)
!bs this routine calculates the coupling of three angular momenta to zero
!bs
!bs   Int dOmega i^(l1+l2+l3) Y^l1_m1 (Omega) Y^l2_m2 (Omega) Y^l3_m3 (Omega) =
!bs   sqrt( (2l1+1)(2l2+1)(2l2+3)/ 4Pi)  * 3J(l1,l2,l3,0,0,0) *
!bs   3J(l1,l2,l3,m1,m2,m3)

use Constants, only: Zero, Quart, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: couple3J
integer(kind=iwp), intent(in) :: l1, l2, l3, m1, m2, m3
integer(kind=iwp) :: l1d, l2d, l3d, m1d, m2d, m3d
real(kind=wp) :: fac1, fac2, fac3
real(kind=wp), parameter :: inv4pi = Quart/Pi
real(kind=wp), external :: regge3j

!bs initialize couple3J-coefficient
couple3J = Zero
!bs quick check
if (m1+m2+m3 /= 0) return
!bs double all values for regge3j
l1d = l1+l1
l2d = l2+l2
l3d = l3+l3
m1d = m1+m1
m2d = m2+m2
m3d = m3+m3
fac1 = sqrt(real(l1d+1,kind=wp)*real(l2d+1,kind=wp)*real(l3d+1,kind=wp)*inv4pi)
fac2 = regge3j(l1d,l2d,l3d,0,0,0)
fac3 = regge3j(l1d,l2d,l3d,m1d,m2d,m3d)
couple3J = fac1*fac2*fac3

return

end function couple3J
