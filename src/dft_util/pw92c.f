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
* Copyright (C) Sergey Gusarov                                         *
************************************************************************
      Subroutine PW92C(Rho,nRho,mGrid,
     &                   dF_dRho,ndF_dRho,
     &                   Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object:  PW92C Functional(Formula taken from Molpro Manual)          *
*                                                                      *
* Called from: Do_batch                                                *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Authors: Gusarov S,,Department of Theoretical Chemistry         *
*              University of LUnd, SWEDEN                              *
*              Maple 8.0                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "nq_index.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),Rho(nRho,mGrid),F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
      Rho_Min=T_X*1.0D-2
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*                                                                      *
************************************************************************
*                                                                      *
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
      rho_a=Rho(1,iGrid)
      rho_b=rho_a
      rho_tot = 2*rho_a
      if(rho_tot.lt.T_X) Go To 100
*
      THREE13 = 3**(1.D0/3.D0)
      FOUR13 = 4**(1.D0/3.D0)
      t1 = FOUR13**2
      t2 = THREE13*t1
      t3 = 2**(1.D0/3.D0)
      t4 = t3**2
      t6 = 1/0.3141592653589793D1
      t8 = (1/rho_a*t6)**(1.D0/3.D0)
      t10 = t2*t4*t8
      t12 = 1+0.267125D-1*t10
      t14 = sqrt(4.D0)
      t15 = sqrt(2.D0)
      t16 = t14*t15
      t17 = sqrt(t10)
      t24 = THREE13**2
      t25 = t24*FOUR13
      t26 = t8**2
      t30 = 0.9494625D0*t16*t17+0.44845D0*t10+0.25596875D-1*t16*t17*t10+
     &0.616175D-1*t25*t3*t26
      t33 = 1+0.1608182432D2/t30
      t34 = log(t33)
      F_term = -0.124364D0*rho_a*t12*t34
      t39 = t2*t4
      t40 = 1/t26
      t41 = rho_a**2
      t42 = 1/t41
      t43 = t40*t42
      t48 = t30**2
      t51 = 4**(1.D0/6.D0)
      t52 = 2**(1.D0/6.D0)
      t53 = t51*t52
      t58 = THREE13*t40*t42*t6
      dFdRhoa = -0.62182D-1*t12*t34+2*rho_a*(0.2768394458D-3*t39*t43*t6*
     &t34+t12/t48*(-0.632975D0*t53/t17*t58-0.7474166668D-
     &1*t39*t43*t6-0.5119375D-1*t53*t17*t58-0.2053916667D-1*t25*t3/t8*t4
     &2*t6)/t33)

      Functional= F_term
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
      dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+Coeff*dFdRhoa

100   Continue
      End Do

*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/= 1
      Else
*                                                                      *
************************************************************************
*                                                                      *
*
      Do iGrid = 1, mGrid
*
      rho_a=Max(Rho_min,Rho(1,iGrid))
      rho_b=Max(Rho_min,Rho(2,iGrid))
      rho_tot = rho_a+rho_b
*
      if(rho_tot.lt.T_X) Go To 200
*
      THREE13 = 3**(1.D0/3.D0)
      FOUR13 = 4**(1.D0/3.D0)
      t1 = FOUR13**2
      t2 = THREE13*t1
      t3 = 1/rho_tot
      t4 = 1/0.3141592653589793D1
      t6 = (t3*t4)**(1.D0/3.D0)
      t7 = t2*t6
      t9 = 1+0.53425D-1*t7
      t10 = sqrt(4.D0)
      t11 = sqrt(t7)
      t12 = t10*t11
      t16 = t10*t11*t7
      t18 = THREE13**2
      t19 = t18*FOUR13
      t20 = t6**2
      t21 = t19*t20
      t23 = 0.1898925D1*t12+0.8969D0*t7+0.1023875D0*t16+0.123235D0*t21
      t26 = 1+0.1608182432D2/t23
      t27 = log(t26)
      t29 = 0.62182D-1*t9*t27
      t31 = 1+0.278125D-1*t7
      t36 = 0.258925D1*t12+
     &      0.905775D0*t7+0.5501625D-1*t16+0.1241775D0*t21
      t39 = 1+0.2960857464D2/t36
      t40 = log(t39)
      t41 = t31*t40
      t42 = rho_a-rho_b
      t43 = t42*t3
      t44 = 1+t43
      t45 = t44**(1.D0/3.D0)
      t47 = 1-t43
      t48 = t47**(1.D0/3.D0)
      t50 = t45*t44+t48*t47-2
      t51 = 2**(1.D0/3.D0)
      t54 = 1/(2*t51-2)
      t55 = t50*t54
      t56 = t42**2
      t57 = t56**2
      t58 = rho_tot**2
      t59 = t58**2
      t60 = 1/t59
      t62 = 1-t57*t60
      t63 = t55*t62
      t65 = 0.197517897D-1*t41*t63
      t67 = 1+0.5137D-1*t7
      t72 = 0.3529725D1*t12+0.1549425D1*t7+
     &      0.2103875D0*t16+0.1562925D0*t21
      t75 = 1+0.3216468318D2/t72
      t76 = log(t75)
      t79 = -0.3109D-1*t67*t76+t29
      t80 = t79*t50
      t81 = t54*t57
      t82 = t81*t60
      t83 = t80*t82
      F_term = rho_tot*(-t29+t65+t83)
      t85 = 1/t20
      t86 = t2*t85
      t87 = 1/t58
      t88 = t87*t4
      t91 = 0.1107357783D-2*t86*t88*t27
      t92 = t23**2
      t95 = 4**(1.D0/6.D0)
      t99 = t85*t87
      t100 = t99*t4
      t101 = t95/t11*THREE13*t100
      t103 = t2*t100
      t107 = t95*t11*THREE13*t100
      t112 = t19/t6*t87*t4
      t118 = t9/t92*(-0.126595D1*t101-0.2989666667D0*t103
     &-0.204775D0*t107-0.8215666667D-1*t112)/t26
      t123 = 0.1831155503D-3*t2*t99*t4*t40*t63
      t124 = t36**2
      t138 = 0.5848223396D0*t31/t124*(-0.1726166667D1*t101-0.301925D0*t1
     &03-0.1100325D0*t107-0.82785D-1*t112)/t39*t50*t54*t62
      t139 = t42*t87
      t140 = t3-t139
      t145 = 4.D0/3.D0*t45*t140-4.D0/3.D0*t48*t140
      t150 = t56*t42
      t151 = t150*t60
      t153 = 1/t59/rho_tot
      t154 = t57*t153
      t163 = t72**2
      t177 = (0.5323644332D-3*t86*t88*t76+0.1D1*t67/t163*(-0.235315D1*t1
     &01-0.516475D0*t103-0.420775D0*t107-0.104195D0*t112)/t75-t91-t118)*
     &t50*t82
      t183 = 4*t80*t54*t150*t60
      t186 = 4*t80*t81*t153
      dFdRhoa = -t29+t65+t83+rho_tot*(t91+t118-t123-t138+0.197517897D-1*
     &t41*t145*t54*t62+0.197517897D-1*t41*t55*(-4*t151+4*t154)+t177+t79*
     &t145*t82+t183-t186)
      t189 = -t3-t139
      t194 = 4.D0/3.D0*t45*t189-4.D0/3.D0*t48*t189
      dFdRhob = -t29+t65+t83+rho_tot*(t91+t118-t123-t138+0.197517897D-1*
     &t41*t194*t54*t62+0.197517897D-1*t41*t55*(4*t151+4*t154)+t177+t79*t
     &194*t82-t183-t186)

      Functional= F_term
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*                                                                      *
************************************************************************
*                                                                      *
*
      dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+Coeff*dFdRhoa
      dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+Coeff*dFdRhob

200   Continue
      End Do
      Endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
