************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Anders Ohrn                                            *
************************************************************************
*  ElPointPot
*
*> @brief
*>   Evaluate the electric potential from point multipoles in a given point
*> @author A. Ohrn
*>
*> @param[in] rinv Inverse distance between multipole and point where the electric potential is to be computed
*> @param[in] x    \f$ x \f$-coordinate for vector between the two relevant points
*> @param[in] y    \f$ y \f$-coordinate of the same vector
*> @param[in] z    \f$ z \f$-coordinate of the same vector
*> @param[in] L    Angular momentum of the multipole, \p L = ``0`` charge, \p L = ``1`` dipole, etc.
*> @param[in] D    The components of the multipole, ordered in the usual Molcas fashion. Should not be in Buckingham form, rather as pure moments
*>
*> @return The electric potential
************************************************************************
      REAL*8 Function ElPointPot(rinv,x,y,z,L,D)
      Implicit REAL*8 (a-h,o-z)

      Dimension D((L+1)*(L+2)/2)

      Pot=0.0d0
      r3=rinv**3
      r5=rinv**5
      x2=x*x
      y2=y*y
      z2=z*z
      r7=rinv**7
      x3=x2*x
      y3=y2*y
      z3=z2*z
      r9=rinv**9
      x4=x3*x
      y4=y3*y
      z4=z3*z
      r11=rinv**11
      x5=x4*x
      y5=y4*y
      z5=z4*z
      If(L.eq.0) then
        Pot=Pot+D(1)*rinv
      Elseif(L.eq.1) then
        Pot=Pot+D(1)*x*r3
        Pot=Pot+D(2)*y*r3
        Pot=Pot+D(3)*z*r3
      Elseif(L.eq.2) then
        Pot=Pot+D(1)*(3.0d0*x2*r5-1.0d0*r3)
        Pot=Pot+2.0d0*D(2)*(3.0d0*x*y*r5)
        Pot=Pot+2.0d0*D(3)*(3.0d0*x*z*r5)
        Pot=Pot+D(4)*(3.0d0*y2*r5-1.0d0*r3)
        Pot=Pot+2.0d0*D(5)*(3.0d0*y*z*r5)
        Pot=Pot+D(6)*(3.0d0*z2*r5-1.0d0*r3)
        Pot=Pot/2.0d0
      Elseif(L.eq.3) then
        Pot=Pot+D(1)*(15.0d0*x3*r7-9.0d0*x*r5)
        Pot=Pot+3.0d0*D(2)*(15.0d0*x2*y*r7-3.0d0*y*r5)
        Pot=Pot+3.0d0*D(3)*(15.0d0*x2*z*r7-3.0d0*z*r5)
        Pot=Pot+3.0d0*D(4)*(15.0d0*x*y2*r7-3.0d0*x*r5)
        Pot=Pot+6.0d0*D(5)*(15.0d0*x*y*z*r7)
        Pot=Pot+3.0d0*D(6)*(15.0d0*x*z2*r7-3.0d0*x*r5)
        Pot=Pot+D(7)*(15.0d0*y3*r7-9.0d0*y*r5)
        Pot=Pot+3.0d0*D(8)*(15.0d0*y2*z*r7-3.0d0*z*r5)
        Pot=Pot+3.0d0*D(9)*(15.0d0*y*z2*r7-3.0d0*y*r5)
        Pot=Pot+D(10)*(15.0d0*z3*r7-9.0d0*z*r5)
        Pot=Pot/6.0d0
      Elseif(L.eq.4) then
        Pot=Pot+D(1)*(105.0d0*x4*r9-90.0d0*x2*r7
     &               +9.0d0*r5)
        Pot=Pot+D(2)*4.0d0*(105.0d0*x3*y*r9-45.0d0*x*y*r7)
        Pot=Pot+D(3)*4.0d0*(105.0d0*x3*z*r9-45.0d0*x*z*r7)
        Pot=Pot+D(4)*6.0d0*(105.0d0*x2*y2*r9-15.0d0*x2*r7
     &               -15.0d0*y2*r7+3.0d0*r5)
        Pot=Pot+D(5)*12.0d0*(105.0d0*x2*y*z*r9-15.0d0*y*z*r7)
        Pot=Pot+D(6)*6.0d0*(105.0d0*x2*z2*r9-15.0d0*x2*r7
     &               -15.0d0*z2*r7+3.0d0*r5)
        Pot=Pot+D(7)*4.0d0*(105.0d0*x*y3*r9-45.0d0*x*y*r7)
        Pot=Pot+D(8)*12.0d0*(105.0d0*x*y2*z*r9-15.0d0*x*z*r7)
        Pot=Pot+D(9)*12.0d0*(105.0d0*x*y*z2*r9-15.0d0*x*y*r7)
        Pot=Pot+D(10)*4.0d0*(105.0d0*x*z3*r9-45.0d0*x*z*r7)
        Pot=Pot+D(11)*(105.0d0*y4*r9-90.0d0*y2*r7
     &                +9.0d0*r5)
        Pot=Pot+D(12)*4.0d0*(105.0d0*y3*z*r9-45.0d0*y*z*r7)
        Pot=Pot+D(13)*6.0d0*(105.0d0*y2*z2*r9-15.0d0*z2*r7
     &                -15.0d0*y2*r7+3.0d0*r5)
        Pot=Pot+D(14)*4.0d0*(105.0d0*y*z3*r9-45.0d0*y*z*r7)
        Pot=Pot+D(15)*(105.0d0*z4*r9-90.0d0*z2*r7
     &                +9.0d0*r5)
        Pot=Pot/24.0d0
      Elseif(L.eq.5) then
        Pot=Pot+D(1)*(945.0d0*x5*r11-1050.0d0*x3*r9+225.0d0*x*r7)
        Pot=Pot+D(2)*(945.0d0*x4*y*r11-630.0d0*x2*y*r9+45.0d0*y*r7)
        Pot=Pot+D(3)*(945.0d0*x4*z*r11-630.0d0*x2*z*r9+45.0d0*z*r7)
        Pot=Pot+D(4)*(945.0d0*x3*y2*r11-315.0d0*x*y2*r9-105.0d0*x3*r9
     &               +45.0d0*x*r7)
        Pot=Pot+D(5)*(945.0d0*x3*y*z*r11-315.0d0*x*y*z*r9)
        Pot=Pot+D(6)*(945.0d0*x3*z2*r11-315.0d0*x*z2*r9-105.0d0*x3*r9
     &               +45.0d0*x*r7)
        Pot=Pot+D(7)*(945.0d0*x2*y3*r11-315.0d0*x2*y*r9-105.0d0*y3*r9
     &               +45.0d0*y*r7)
        Pot=Pot+D(8)*(945.0d0*x2*y2*z*r11-105.0d0*y2*z*r9
     &               -105.0d0*x2*z*r9+15.0d0*z*r7)
        Pot=Pot+D(9)*(945.0d0*x2*y*z2*r11-105.0d0*z2*y*r9
     &               -105.0d0*x2*y*r9+15.0d0*y*r7)
        Pot=Pot+D(10)*(945.0d0*x2*z3*r11-315.0d0*x2*z*r9-105.0d0*z3*r9
     &                +45.0d0*z*r7)
        Pot=Pot+D(11)*(945.0d0*x*y4*r11-630.0d0*x*y2*r9+45.0d0*x*r7)
        Pot=Pot+D(12)*(945.0d0*x*y3*z*r11-315.0d0*x*y*z*r9)
        Pot=Pot+D(13)*(945.0d0*x*y2*z2*r11-105.0d0*x*y2*r9
     &                -105.0d0*x*z2*r9+15.0d0*x*r7)
        Pot=Pot+D(14)*(945.0d0*x*y*z3*r11-315.0d0*x*y*z*r9)
        Pot=Pot+D(15)*(945.0d0*x*z4*r11-630.0d0*x*z2*r9+45.0d0*x*r7)
        Pot=Pot+D(16)*(945.0d0*y5*r11-1050.0d0*y3*r9+225.0d0*y*r7)
        Pot=Pot+D(17)*(945.0d0*y4*z*r11-630.0d0*y2*z+45.0d0*z*r7)
        Pot=Pot+D(18)*(945.0d0*y3*z2*r11-315.0d0*y*z2*r9-105.0d0*y3*r9
     &                +45.0d0*y*r9)
        Pot=Pot+D(19)*(945.0d0*y2*z3*r11-315.0d0*y2*z*r9-105.0d0*z3*r9
     &                +45.0d0*z*r9)
        Pot=Pot+D(20)*(945.0d0*y*z4*r11-630.0d0*y*z2*r9+45.0d0*y*r7)
        Pot=Pot+D(21)*(9450.0d0*z5*r11-1050.0d0*z3*r9+225.0d0*z*r7)
        Pot=Pot/120.0d0
      Endif

      ElPointPot=Pot

      Return
      End
