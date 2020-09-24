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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine VWN_V(mGrid,Rho,nRho,iSpin,
     &                 F_xc,dF_dRho,ndF_dRho,Coeff,T_X)
************************************************************************
*                                                                      *
* Object: To compute functional V from the VWN80 paper which fits the  *
*         Ceperley-Alder solution to the uniform electron gas.         *
*         This is the functional to use in LDA and LSDA.               *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
*  LDA Stuff
      Real*8 A(3),b(3),c(3),x0(3),Q(3),Xx0(3),e(3),d_e(3)
      data A  / 0.0621814D+00,  0.0310907D+00, -0.0337740D+00/
      data b  / 3.7274400D+00,  7.0604200D+00,  1.1310700D+00/
      data c  /12.9352000D+00, 18.0578000D+00, 13.0045000D+00/
      data x0 /-0.1049800D+00, -0.3250000D+00, -0.0047584D+00/
*     data Xx0/12.5549000D+00, 15.8687880D+00, 12.9991400D+00/
      data f_d_2/1.709920D+00/
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
      THIRD =One/Three
      FTHIRD=Four/Three
C     CVX   =(two**THIRD)*((three/Pi)**THIRD)
      Fact = One / (Two**FTHIRD-Two)
      Q(1)=Sqrt(Four*c(1)-b(1)**2)
      Q(2)=Sqrt(Four*c(2)-b(2)**2)
      Q(3)=Sqrt(Four*c(3)-b(3)**2)
      Xx0(1)=(x0(1)**2+b(1)*x0(1)+c(1))
      Xx0(2)=(x0(2)**2+b(2)*x0(2)+c(2))
      Xx0(3)=(x0(3)**2+b(3)*x0(3)+c(3))
*
      z     = Zero ! dummy initialize
      d_z_a = Zero ! dummy initialize
      d_z_b = Zero ! dummy initialize
      fz    = Zero ! dummy initialize
*
      Rho_min=T_X*1.0D-2
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute value of energy and integrad on the grid
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
         d_alpha =Rho(ipR,iGrid)
         DTot=Two*d_alpha
         If (DTot.lt.T_X) Go To 100
*
*------- Accumulate Correlation contributions to the AO intergrals
*
         r_s=(Three/(DTot*Pi*Four))**THIRD
         x=Sqrt(r_s)
         d_x=-x/(Six*DTot)
*
         Xx   =  x**2     + b(1)*x + c(1)
         d_Xx = Two*x     + b(1)
         d_Xx0= Two*x0(1) + b(1)
*
         e(1)     = A(1) * (
     &                LOG( x**2/Xx )
     &            + Two * b(1)/Q(1)
     &              * ( One - x0(1)*d_Xx0/Xx0(1) )
     &              * ATAN( Q(1)/d_Xx )
     &            - b(1)*x0(1)/Xx0(1)
     &              * LOG( (x-x0(1))**2/Xx )
     &                      )
*
         xtemp =-(Two * d_x) / (Q(1)**2+d_Xx**2)
*
         d_e(1)   = A(1) * (
     &                (Two*Xx-x*d_Xx)*d_x/(x*Xx)
     &            + Two*b(1)
     &              * ( One - x0(1)*d_Xx0/Xx0(1) )
     &              * xtemp
     &            - b(1)*x0(1)
     &              * (Two*Xx-(x-x0(1))*d_Xx)
     &              * d_x / (Xx*xx0(1)*(x-x0(1)))
     &                      )
*
*------- original formula in Rydbergs -> 0.5 converts to hartree
*
         half_func = e(1)
*
         F_xc(iGrid)=F_xc(iGrid)+Coeff*Half*half_func*DTot
*
         func_d_rho_c = DTot*d_e(1)
*
         dF_dRho(ipR,iGrid) = dF_dRho(ipR,iGrid)
     &                       + Coeff*Half*(func_d_rho_c
     &                                            +half_func)
*
 100     Continue
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/=1
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
         d_alpha =Max(Rho_min,Rho(ipRa,iGrid))
         d_beta  =Max(Rho_min,Rho(ipRb,iGrid))
         DTot=d_alpha+d_beta
         If (DTot.lt.T_X) Go To 200
         DSpn=d_alpha-d_beta
         nSpin=3
*
*------- Accumulate Correlation contributions to the AO intergrals
*
         r_s=(Three/(DTot*Pi*Four))**THIRD
         x=Sqrt(r_s)
         d_x=-x/(Six*DTot)
*
         Do j=1,nSpin
            Xx   =  x**2     + b(j)*x + c(j)
            d_Xx = Two*x     + b(j)
            d_Xx0= Two*x0(j) + b(j)
*
            e(j)     = A(j) * (
     &                   LOG( x**2/Xx )
     &               + Two * b(j)/Q(j)
     &                 * ( One - x0(j)*d_Xx0/Xx0(j) )
     &                 * ATAN( Q(j)/d_Xx )
     &               - b(j)*x0(j)/Xx0(j)
     &                 * LOG( (x-x0(j))**2/Xx )
     &                         )
*
            xtemp =-(Two * d_x) / (Q(j)**2+d_Xx**2)
*
            d_e(j)   = A(j) * (
     &                   (Two*Xx-x*d_Xx)*d_x/(x*Xx)
     &               + Two*b(j)
     &                 * ( One - x0(j)*d_Xx0/Xx0(j) )
     &                 * xtemp
     &               - b(j)*x0(j)
     &                 * (Two*Xx-(x-x0(j))*d_Xx)
     &                 * d_x / (Xx*xx0(j)*(x-x0(j)))
     &                         )
         End Do
*
*------- original formula in Rydbergs -> 0.5 converts to hartree
*
         z = DSpn/DTot
         d_z_a= Two* d_beta /DTot**2
         d_z_b=-Two* d_alpha/DTot**2
         fz=( (One+z)**FTHIRD
     &       +(One-z)**FTHIRD-Two ) * Fact
         half_func = e(1)
     &             + ((e(3)*fz) / f_d_2) * (One-z**4)
     &             + fz*( e(2)- e(1)) * z**4
*
         F_xc(iGrid)=F_xc(iGrid)+Coeff*Half*half_func*DTot
*
         d_fz= Fact * FTHIRD * d_z_a
     &       *( (One+z)**THIRD - (One-z)**THIRD )
         d_z4=Four*z**3*d_z_a
*
         d_half_func_a = d_e(1)
     &                 + (d_e(3)*fz+e(3)*d_fz) /f_d_2* (One-z**4)
     &                 -       (e(3)*fz)       /f_d_2*    d_z4
     &                 + d_fz*(  e(2)-  e(1))*z**4
     &                 +   fz*(d_e(2)-d_e(1))*z**4
     &                 +   fz*(  e(2)-  e(1))*d_z4
*
c+++  original formula in Rydbergs -> 0.5 converts to hartree
*
         func_d_rho_alpha_c =  DTot*d_half_func_a
*
         d_fz= Fact * FTHIRD * d_z_b
     &       *( (One+z)**THIRD - (One-z)**THIRD )
         d_z4=Four*z**3*d_z_b
*
         d_half_func_b = d_e(1)
     &                 + (d_e(3)*fz+e(3)*d_fz) /f_d_2* (One-z**4)
     &                 -       (e(3)*fz)       /f_d_2*    d_z4
     &                 + d_fz*(  e(2)-  e(1))*z**4
     &                 +   fz*(d_e(2)-d_e(1))*z**4
     &                 +   fz*(  e(2)-  e(1))*d_z4
         func_d_rho_beta_c  = DTot*d_half_func_b
*
         dF_dRho(ipRa,iGrid) = dF_dRho(ipRa,iGrid)
     &                       + Coeff*Half*(func_d_rho_alpha_c
     &                                            +half_func)
         dF_dRho(ipRb,iGrid) = dF_dRho(ipRb,iGrid)
     &                       + Coeff*Half*(func_d_rho_beta_c
     &                                            +half_func)
*
 200     Continue
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
