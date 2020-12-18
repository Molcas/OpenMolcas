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
*               2001, Laura Gagliardi                                  *
************************************************************************
      Subroutine VWN_III(mGrid,Rho,nRho,iSpin,F_xc,
     &                   dF_dRho,ndF_dRho,Coeff,T_X)
************************************************************************
*                                                                      *
* Object: To compute functional III from the VWN80 paper which fits the*
*         Ceperley-Alder solution to the uniform electron gas.         *
*         This is the functional to use in B3LYP.                      *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
*             Laura Gagliardi, Dipartimento di Chimica G. Ciamician,   *
*             University of Bologna, ITALY. October 2001               *
************************************************************************
      use KSDFT_Info, only: tmpB
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
*  LDA Stuff
      Real*8 A(2),b(2),c(2),x0(2),Q(2),xx0(2),e(2),d_e(2)
      data A  / 0.0621814D+00,  0.0310907D+00/
      data x0 /-0.4092860D+00, -0.7432940D+00/
      data b  /13.0720000D+00, 20.1231000D+00/
      data c  /42.7198000D+00,101.5780000D+00/
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
      THIRD =One/Three
      FTHIRD=Four/Three
      Fact = One / (Two**FTHIRD-Two)
      Q(1)=Sqrt(Four*c(1)-b(1)**2)
      Q(2)=Sqrt(Four*c(2)-b(2)**2)
      Xx0(1)=(x0(1)**2+b(1)*x0(1)+c(1))
      Xx0(2)=(x0(2)**2+b(2)*x0(2)+c(2))
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
         If (DTot.le.T_X) Go To 100
*
*------- Accumulate Correlation contributions to the AO intergrals
*
         r_s=(Three/(DTot*Pi*Four))**THIRD
         x=Sqrt(r_s)
         d_x=-x/(Six*DTot)
*
         j=1
         Xx           =  x*x     + b(1)*x + c(1)
         d_Xx         = Two*x     + b(1)
         d_Xx0        = Two*x0(1) + b(1)
         x0j_Xx0j     = x0(1)/Xx0(1)
         xmx0j        = x-x0(1)
*
         e(1)     = A(1) * (
     &                LOG( x*x/Xx )
     &            + Two * b(1)/Q(1)
     &              * ( One - x0j_Xx0j*d_Xx0 )
     &              * ATAN( Q(1)/d_Xx )
     &            - b(1)*x0j_Xx0j
     &              * LOG( xmx0j*xmx0j/Xx )
     &                      )
*
         xtemp =-(Two * d_x) / (Q(1)*Q(1)+d_Xx*d_Xx)
*
         d_e(1)   = A(1) * (
     &                (Two*Xx-x*d_Xx)*d_x/(x*Xx)
     &            + Two*b(1)
     &              * ( One - x0(1)*d_Xx0/Xx0(1) )
     &              * xtemp
     &            - b(1)*x0(1)
     &              * (Two*Xx-xmx0j*d_Xx)
     &              * d_x / (Xx*xx0(1)*xmx0j)
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
     &                               + Coeff*Half*(func_d_rho_c
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
         rho_a=d_alpha
         rho_b=d_beta
         DTot=d_alpha+d_beta
         If (DTot.le.T_X) Go To 200
         DSpn=d_alpha-d_beta
      If(.True.) Then
c      If(.False.) Then
      t1 = rho_a+rho_b
      t2 = 1/t1
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908D0*t3
      t5 = t2**(1.D0/6.D0)
      t7 = t4+0.1029581201D2*t5+0.427198D2
      t8 = 1/t7
      t11 = log(0.6203504908D0*t3*t8)
      t12 = 0.621814D-1*t11
      t13 = 0.1575246636D1*t5
      t14 = t13+0.13072D2
      t17 = atan(0.4489988864D-1/t14)
      t18 = 0.4104394588D2*t17
      t19 = 0.787623318D0*t5
      t20 = t19+0.409286D0
      t21 = t20**2
      t23 = log(t21*t8)
      t24 = 0.8862747535D-2*t23
      t26 = t4+0.1584942279D2*t5+0.101578D3
      t27 = 1/t26
      t30 = log(0.6203504908D0*t3*t27)
      t31 = 0.310907D-1*t30
      t32 = t13+0.201231D2
      t35 = atan(0.1171685282D1/t32)
      t36 = 0.1237636055D1*t35
      t37 = t19+0.743294D0
      t38 = t37**2
      t40 = log(t38*t27)
      t41 = 0.5334620013D-2*t40
      t43 = (rho_a-rho_b)*t2

      t44 = 1+t43
      t45 = t44**(1.D0/3.D0)
      t47 = 1-t43
      t48 = t47**(1.D0/3.D0)
      t51 = 2**(1.D0/3.D0)
      t56 = rho_a-1.D0*rho_b
      t57 = t56*t2
      t58 = 1.D0+t57
      t59 = t58**(1.D0/3.D0)
      t62 = 1.D0-1.D0*t57
      t63 = t62**(1.D0/3.D0)
      t65 = t59*t58+t63*t62-2.D0
      t66 = t31+t36+t41-t12-t18-t24
      t67 = t65*t66
      vwn = 0.5D0*(t12+t18+t24+0.192366105D1*t67)*t1
      t71 = t3**2
      t72 = 1/t71
      t74 = t1**2
      t75 = 1/t74
      t78 = t7**2
      t79 = 1/t78
      t82 = 0.2067834969D0*t72*t75
      t83 = t5**2

      t84 = t83**2
      t86 = 1/t84/t5
      t87 = t86*t75
      t89 = -t82-0.1715968668D1*t87
      t93 = 1/t3
      t96 = 0.1002359165D0*(-0.2067834969D0*t72*t8*t75-0.6203504908D0*t3
     &*t79*t89)*t93*t7
      t97 = t14**2
      t98 = 1/t97
      t105 = 0.4838287602D0*t98*t86*t75/(1+0.2016D-2*t98)
      t115 = 0.8862747535D-2*(-0.262541106D0*t20*t8*t87-t21*t79*t89)/t21
     &*t7
      t116 = t56*t75
      t119 = 1.D0*t2
      t120 = 1.D0*t116
      t130 = t26**2
      t131 = 1/t130
      t134 = -t82-0.2641570465D1*t87
      t141 = t32**2
      t142 = 1/t141
      t162 = 0.192366105D1*t65*(0.5011795824D-1*(-0.2067834969D0*t72*t27

     &*t75-0.6203504908D0*t3*t131*t134)*t93*t26+0.3807160955D0*t142*t86*
     &t75/(1+0.13728464D1*t142)+0.5334620013D-2*(-0.262541106D0*t37*t27*
     &t87-t38*t131*t134)/t38*t26-t96-t105-t115)
      t166 = 0.310907D-1*t11
      t167 = 0.2052197294D2*t17
      t168 = 0.4431373768D-2*t23
      t169 = 0.961830525D0*t67
      dvwndra = 0.5D0*(t96+t105+t115+0.192366105D1*(4.D0/3.D0*t59*(t2-t1
     &16)+4.D0/3.D0*t63*(-t119+t120))*t66+t162)*t1+t166+t167+t168+t169
      dvwndrb = 0.5D0*(t96+t105+t115+0.192366105D1*(4.D0/3.D0*t59*(-t119
     &-t116)+4.D0/3.D0*t63*(t119+t120))*t66+t162)*t1+t166+t167+t168+t169

*
         F_xc(iGrid)=F_xc(iGrid)+Coeff*vwn
         tmpB(iGrid)=tmpB(iGrid)+Coeff*vwn
*
         dF_dRho(ipRa,iGrid) = dF_dRho(ipRa,iGrid)
     &                               + Coeff*dvwndra
         dF_dRho(ipRb,iGrid) = dF_dRho(ipRb,iGrid)
     &                               + Coeff*dvwndrb
      Else

*
*------- Accumulate Correlation contributions to the AO intergrals
*
         r_s=(Three/(DTot*Pi*Four))**THIRD
         x=Sqrt(r_s)
         d_x=-x/(Six*DTot)
*
         Do j=1,2
            Xx           =  x*x     + b(j)*x + c(j)
            d_Xx         = Two*x     + b(j)
            d_Xx0        = Two*x0(j) + b(j)
            x0j_Xx0j     = x0(j)/Xx0(j)
            xmx0j        = x-x0(j)
*
            e(j)     = A(j) * (
     &                   LOG( x*x/Xx )
     &               + Two * b(j)/Q(j)
     &                 * ( One - x0j_Xx0j*d_Xx0 )
     &                 * ATAN( Q(j)/d_Xx )
     &               - b(j)*x0j_Xx0j
     &                 * LOG( xmx0j*xmx0j/Xx )
     &                         )
*
            xtemp =-(Two * d_x) / (Q(j)*Q(j)+d_Xx*d_Xx)
*
            d_e(j)   = A(j) * (
     &                   (Two*Xx-x*d_Xx)*d_x/(x*Xx)
     &               + Two*b(j)
     &                 * ( One - x0(j)*d_Xx0/Xx0(j) )
     &                 * xtemp
     &               - b(j)*x0(j)
     &                 * (Two*Xx-xmx0j*d_Xx)
     &                 * d_x / (Xx*xx0(j)*xmx0j)
     &                         )
         End Do
*
*------- original formula in Rydbergs -> 0.5 converts to hartree
*
         e2me1    =    e(2)- e(1)
         de2mde1  =  d_e(2)-d_e(1)
         z = DSpn/DTot
         d_z_a= Two* d_beta /DTot*DTot
         d_z_b=-Two* d_alpha/DTot*DTot
         fz=( (One+z)**FTHIRD
     &       +(One-z)**FTHIRD-Two ) * Fact
         half_func = e(1) + fz * e2me1 !* z**4
*
         F_xc(iGrid)=F_xc(iGrid)+Coeff*Half*half_func*DTot
*
         d_fz= Fact * FTHIRD * d_z_a
     &       *( (One+z)**THIRD - (One-z)**THIRD )
*
         d_half_func_a = d_e(1)
     &                 + d_fz*e2me1
     &                 +   fz*de2mde1
*
c+++  original formula in Rydbergs -> 0.5 converts to hartree
*
         func_d_rho_alpha_c =  DTot*d_half_func_a

         d_fz= Fact * FTHIRD * d_z_b
     &       *( (One+z)**THIRD - (One-z)**THIRD )
*
         d_half_func_b = d_e(1)
     &                 + d_fz*e2me1
     &                 +   fz*de2mde1
         func_d_rho_beta_c  = DTot*d_half_func_b
*
         dF_dRho(ipRa,iGrid) = dF_dRho(ipRa,iGrid)
     &                               + Coeff*Half*(func_d_rho_alpha_c
     &                                            +half_func)
         dF_dRho(ipRb,iGrid) = dF_dRho(ipRb,iGrid)
     &                               + Coeff*Half*(func_d_rho_beta_c
     &                                            +half_func)

      End If
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
