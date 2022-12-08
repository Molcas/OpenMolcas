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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  ElPointPot
!
!> @brief
!>   Evaluate the electric potential from point multipoles in a given point
!> @author A. Ohrn
!>
!> @param[in] rinv Inverse distance between multipole and point where the electric potential is to be computed
!> @param[in] x    \f$ x \f$-coordinate for vector between the two relevant points
!> @param[in] y    \f$ y \f$-coordinate of the same vector
!> @param[in] z    \f$ z \f$-coordinate of the same vector
!> @param[in] L    Angular momentum of the multipole, \p L = ``0`` charge, \p L = ``1`` dipole, etc.
!> @param[in] D    The components of the multipole, ordered in the usual Molcas fashion. Should not be in Buckingham form,
!>                 rather as pure moments
!>
!> @return The electric potential
!***********************************************************************

function ElPointPot(rinv,x,y,z,L,D)

use Constants, only: Zero, Two, Three, Four, Six, Nine, Twelve
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: ElPointPot
integer(kind=iwp), intent(in) :: L
real(kind=wp), intent(in) :: rinv, x, y, z, D((L+1)*(L+2)/2)
real(kind=wp) :: Pot, r3, r5, r7, r9, r11, x2, x3, x4, x5, y2, y3, y4, y5, z2, z3, z4, z5
real(kind=wp), parameter :: r15 = 15.0_wp, r24 = 24.0_wp, r45 = 45.0_wp, r90 = 90.0_wp, r105 = 105.0_wp, r120 = 120.0_wp, &
                            r225 = 225.0_wp, r315 = 315.0_wp, r630 = 630.0_wp, r945 = 945.0_wp, r1050 = 1050.0_wp

Pot = Zero
r3 = rinv**3
r5 = rinv**5
x2 = x*x
y2 = y*y
z2 = z*z
r7 = rinv**7
x3 = x**3
y3 = y**3
z3 = z**3
r9 = rinv**9
x4 = x**4
y4 = y**4
z4 = z**4
r11 = rinv**11
x5 = x**5
y5 = y**5
z5 = z**5
if (L == 0) then
  Pot = Pot+D(1)*rinv
else if (L == 1) then
  Pot = Pot+D(1)*x*r3
  Pot = Pot+D(2)*y*r3
  Pot = Pot+D(3)*z*r3
else if (L == 2) then
  Pot = Pot+D(1)*(Three*x2*r5-r3)
  Pot = Pot+Two*D(2)*(Three*x*y*r5)
  Pot = Pot+Two*D(3)*(Three*x*z*r5)
  Pot = Pot+D(4)*(Three*y2*r5-r3)
  Pot = Pot+Two*D(5)*(Three*y*z*r5)
  Pot = Pot+D(6)*(Three*z2*r5-r3)
  Pot = Pot/Two
else if (L == 3) then
  Pot = Pot+D(1)*(r15*x3*r7-Nine*x*r5)
  Pot = Pot+Three*D(2)*(r15*x2*y*r7-Three*y*r5)
  Pot = Pot+Three*D(3)*(r15*x2*z*r7-Three*z*r5)
  Pot = Pot+Three*D(4)*(r15*x*y2*r7-Three*x*r5)
  Pot = Pot+Six*D(5)*(r15*x*y*z*r7)
  Pot = Pot+Three*D(6)*(r15*x*z2*r7-Three*x*r5)
  Pot = Pot+D(7)*(r15*y3*r7-Nine*y*r5)
  Pot = Pot+Three*D(8)*(r15*y2*z*r7-Three*z*r5)
  Pot = Pot+Three*D(9)*(r15*y*z2*r7-Three*y*r5)
  Pot = Pot+D(10)*(r15*z3*r7-Nine*z*r5)
  Pot = Pot/Six
else if (L == 4) then
  Pot = Pot+D(1)*(r105*x4*r9-r90*x2*r7+Nine*r5)
  Pot = Pot+D(2)*Four*(r105*x3*y*r9-r45*x*y*r7)
  Pot = Pot+D(3)*Four*(r105*x3*z*r9-r45*x*z*r7)
  Pot = Pot+D(4)*Six*(r105*x2*y2*r9-r15*x2*r7-r15*y2*r7+Three*r5)
  Pot = Pot+D(5)*Twelve*(r105*x2*y*z*r9-r15*y*z*r7)
  Pot = Pot+D(6)*Six*(r105*x2*z2*r9-r15*x2*r7-r15*z2*r7+Three*r5)
  Pot = Pot+D(7)*Four*(r105*x*y3*r9-r45*x*y*r7)
  Pot = Pot+D(8)*Twelve*(r105*x*y2*z*r9-r15*x*z*r7)
  Pot = Pot+D(9)*Twelve*(r105*x*y*z2*r9-r15*x*y*r7)
  Pot = Pot+D(10)*Four*(r105*x*z3*r9-r45*x*z*r7)
  Pot = Pot+D(11)*(r105*y4*r9-r90*y2*r7+Nine*r5)
  Pot = Pot+D(12)*Four*(r105*y3*z*r9-r45*y*z*r7)
  Pot = Pot+D(13)*Six*(r105*y2*z2*r9-r15*z2*r7-r15*y2*r7+Three*r5)
  Pot = Pot+D(14)*Four*(r105*y*z3*r9-r45*y*z*r7)
  Pot = Pot+D(15)*(r105*z4*r9-r90*z2*r7+Nine*r5)
  Pot = Pot/r24
else if (L == 5) then
  Pot = Pot+D(1)*(r945*x5*r11-r1050*x3*r9+r225*x*r7)
  Pot = Pot+D(2)*(r945*x4*y*r11-r630*x2*y*r9+r45*y*r7)
  Pot = Pot+D(3)*(r945*x4*z*r11-r630*x2*z*r9+r45*z*r7)
  Pot = Pot+D(4)*(r945*x3*y2*r11-r315*x*y2*r9-r105*x3*r9+r45*x*r7)
  Pot = Pot+D(5)*(r945*x3*y*z*r11-r315*x*y*z*r9)
  Pot = Pot+D(6)*(r945*x3*z2*r11-r315*x*z2*r9-r105*x3*r9+r45*x*r7)
  Pot = Pot+D(7)*(r945*x2*y3*r11-r315*x2*y*r9-r105*y3*r9+r45*y*r7)
  Pot = Pot+D(8)*(r945*x2*y2*z*r11-r105*y2*z*r9-r105*x2*z*r9+r15*z*r7)
  Pot = Pot+D(9)*(r945*x2*y*z2*r11-r105*z2*y*r9-r105*x2*y*r9+r15*y*r7)
  Pot = Pot+D(10)*(r945*x2*z3*r11-r315*x2*z*r9-r105*z3*r9+r45*z*r7)
  Pot = Pot+D(11)*(r945*x*y4*r11-r630*x*y2*r9+r45*x*r7)
  Pot = Pot+D(12)*(r945*x*y3*z*r11-r315*x*y*z*r9)
  Pot = Pot+D(13)*(r945*x*y2*z2*r11-r105*x*y2*r9-r105*x*z2*r9+r15*x*r7)
  Pot = Pot+D(14)*(r945*x*y*z3*r11-r315*x*y*z*r9)
  Pot = Pot+D(15)*(r945*x*z4*r11-r630*x*z2*r9+r45*x*r7)
  Pot = Pot+D(16)*(r945*y5*r11-r1050*y3*r9+r225*y*r7)
  Pot = Pot+D(17)*(r945*y4*z*r11-r630*y2*z+r45*z*r7)
  Pot = Pot+D(18)*(r945*y3*z2*r11-r315*y*z2*r9-r105*y3*r9+r45*y*r9)
  Pot = Pot+D(19)*(r945*y2*z3*r11-r315*y2*z*r9-r105*z3*r9+r45*z*r9)
  Pot = Pot+D(20)*(r945*y*z4*r11-r630*y*z2*r9+r45*y*r7)
  Pot = Pot+D(21)*(r945*z5*r11-r1050*z3*r9+r225*z*r7)
  Pot = Pot/r120
end if

ElPointPot = Pot

return

end function ElPointPot
