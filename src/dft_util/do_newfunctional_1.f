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
* Copyright (C) Per Ake Malmqvist                                      *
*               2004, Giovanni Ghigo                                   *
************************************************************************
      Subroutine Do_NewFunctional_1(Rho,nRho,mGrid,
     &                   dF_dRho,ndF_dRho,
     &                   Coeff,nD,F_xc,
     &                   P2_ontop,nP2_ontop,
     &                   dF_dP2ontop,ndF_dP2ontop,T_X)
************************************************************************
*                                                                      *
* Object:  Lyp Functional(Formula taken from Molpro Manual)            *
*                                                                      *
*      Author: Per-AAke Malmquist,Department of Theoretical Chemistry  *
*              University of LUnd, SWEDEN                              *
*      Modified by G. Ghigo, Department of Theoretical Chemistry,      *
*                  University of Lund, SWEDEN. June 2004               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
#include "pamint.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),
     &       Rho(nRho,mGrid),P2_ontop(nP2_ontop,mGrid),
     &       F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
      Common /points/ GGrid(3,6000)
      Common /ipoints/ iiGrid
*
      aConst=0.049180d0*Acoef
      bConst=0.1320d0*Bcoef
      cConst=0.25330d0*Ccoef
      dConst=0.3490d0*Dcoef
CGG   Cf=2.8712340001881918D0
      Rho_min=T_X*1.0D-2
*
      Do iGrid=1,mGrid
*
       If(nD.eq.1) Then
         Rho_tot=2.0d0*Rho(1,iGrid)
         If (Rho_tot.lt.T_X) Go To 199
         Rhox=2.0d0*Rho(2,iGrid)
         Rhoy=2.0d0*Rho(3,iGrid)
         Rhoz=2.0d0*Rho(4,iGrid)
         gradRho2=Rhox*Rhox+Rhoy*Rhoy+Rhoz*Rhoz
         grad2Rho=2.0d0*Rho(5,iGrid)
         P2=P2_ontop(1,iGrid)
         P2x=P2_ontop(2,iGrid)
         P2y=P2_ontop(3,iGrid)
         P2z=P2_ontop(4,iGrid)
         gradP22=P2x*P2x+P2y*P2y+P2z*P2z
         grad2P2=P2_ontop(5,iGrid)
         P2sp=P2_ontop(6,iGrid)
         gradP2gradRho=Rhox*P2x+Rhoy*P2y+Rhoz*P2z
       Else
         Rho_tot=Max(Rho_Min,Rho(1,iGrid))+Max(Rho_Min,Rho(2,iGrid))
         If (Rho_tot.lt.T_X) Go To 199
         Rhox=Rho(3,iGrid)+Rho(6,iGrid)
         Rhoy=Rho(4,iGrid)+Rho(7,iGrid)
         Rhoz=Rho(5,iGrid)+Rho(8,iGrid)
         gradRho2=Rhox*Rhox+Rhoy*Rhoy+Rhoz*Rhoz
         grad2Rho=Rho(9,iGrid)+Rho(10,iGrid)
         P2 =P2_ontop(1,iGrid)
         P2x=P2_ontop(2,iGrid)
         P2y=P2_ontop(3,iGrid)
         P2z=P2_ontop(4,iGrid)
         gradP22=P2x*P2x+P2y*P2y+P2z*P2z
         grad2P2=P2_ontop(5,iGrid)
         P2sp=P2_ontop(6,iGrid)
         gradP2gradRho=Rhox*P2x+Rhoy*P2y+Rhoz*P2z
       End If
*
*    T1 term:
*
      t1 = Rho_tot**(1.D0/3.D0)
      t4 = 1+dConst/t1
      t5 = 1/t4
      denom = 1.D0*t5
      t6 = aConst*P2
      t7 = 1/Rho_tot
      T1_term = -4*t6*t5*t7
      t11 = t4**2
      t13 = Rho_tot**2
      dT1_dRho = -4.D0/3.D0*t6/t11/t1/t13*dConst+4*t6*t5/t13
      dT1_dP2 = -4*aConst*t5*t7
*
*    T2 term:
*
c      t1 = Rho_tot**(1.D0/3.D0)
c      t2 = 1/t1
c      t4 = 1+dConst*t2
c      t5 = 1/t4
c      denom = 1.D0*t5
c      T2_term = -aConst*Rho_tot*t5
c      t10 = t4**2
c      dT2_dRho = -aConst*t5-aConst*t2/t10*dConst/3
*
*   T3 term:
*
      t1 = aConst*bConst
      t2 = Rho_tot**(1.D0/3.D0)
      t3 = 1/t2
      t5 = exp(-cConst*t3)
      t7 = 1+dConst*t3
      t8 = 1/t7
      t9 = t5*t8
      t10 = t2**2
      t12 = 1/t10/Rho_tot
      t14 = t1*t9*t12
      Z_term = -t14
      t15 = t1*cConst
      t16 = Rho_tot**2
      t17 = t16*Rho_tot
      t18 = 1/t17
      t23 = t1*t5
      t24 = t7**2
      t25 = 1/t24
      t35 = -t15*t18*t5*t8/3-t23*t25*t18*dConst/3+5.D0/3.D0*t1*
     &         t9/t10/t16
      T3_term = -t35*gradP2gradRho/4+t23*t8*t12*P2sp
      dT3_dRho = t35*(grad2P2/4-P2sp)
      t43 = t16**2
      t44 = 1/t43
      t49 = cConst**2
      t52 = 1/t2/t43
      t53 = t52*t5
      t64 = dConst**2
      dT3_dP2 = (14.D0/9.D0*t15*t44*t5*t8-t1*t49*t53*t8/9-2.D0/9.D0*t15*
     &t53*t25*dConst-2.D0/9.D0*t23/t24/t7*t52*t64+14.D0/9.D0*t23*t25*t44
     &*dConst-40.D0/9.D0*t1*t9/t10/t17)*gradRho2/4+t35*grad2Rho/4
      dT3_dP2sp = t14
*
      Functional= T1_term   +   T3_term
      F_xc(iGrid)=F_xc(iGrid)+Coeff*Functional
*
      dF_dRho(ipR,iGrid)=     dT1_dRho + dT3_dRho
      dF_dP2ontop(1,iGrid)= dT1_dP2  + dT3_dP2
      dF_dP2ontop(2,iGrid)= 0.0D0
      dF_dP2ontop(3,iGrid)= 0.0D0
      dF_dP2ontop(4,iGrid)= 0.0D0
      dF_dP2ontop(5,iGrid)= 0.0D0
      dF_dP2ontop(6,iGrid)= dT3_dP2sp
*
 199  Continue
      End Do  ! iGrid
c
      Return
      End
