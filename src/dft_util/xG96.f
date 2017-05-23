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
* Copyright (C) 2006, Per Ake Malmqvist                                *
*               2009, Grigory A. Shamov                                *
************************************************************************
      Subroutine xG96(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object Gill-1996 Exchange: Peter Gill, Mol Phys 1996 vol 89, 433     *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. June 2006  (B88 code)        *
*             Grigory A. Shamov, U of Manitoba, Dec 2009               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)

* IDORD=Order of derivatives to request from XPBE:
      idord=1
*
      Rho_Min=T_X*1.0D-2
*
      if (ispin.eq.1) then
* ispin=1 means spin zero.

* T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=rho(ipR,iGrid)
         if(rhoa.lt.Ta) goto 110
         grdrhoa_x=rho(ipdRx,iGrid)
         grdrhoa_y=rho(ipdRy,iGrid)
         grdrhoa_z=rho(ipdRz,iGrid)
         sigmaaa=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2

         call xG96_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(2.0D0*Fa)
         dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+Coeff*dFdrhoa
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+Coeff*dFdgammaaa
* Note: For xpbe, dFdgammaab is zero.
 110     continue
        end do

      else
* ispin .ne. 1, use both alpha and beta components.

        do iGrid=1,mgrid
         rhoa=Max(Rho_Min,rho(ipRa,iGrid))
         rhob=Max(Rho_Min,rho(ipRb,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         grdrhoa_x=rho(ipdRxa,iGrid)
         grdrhoa_y=rho(ipdRya,iGrid)
         grdrhoa_z=rho(ipdRza,iGrid)
         sigmaaa=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
         call xG96_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=rho(ipdRxb,iGrid)
         grdrhob_y=rho(ipdRyb,iGrid)
         grdrhob_z=rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call xG96_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(Fa+Fb)
         dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+Coeff*dFdrhoa
         dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+Coeff*dFdrhob
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
* Note: For xpbe, dFdgammaab is zero.
         dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+Coeff*dFdgammaaa
         dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+Coeff*dFdgammabb
 210     continue
        end do

      end if

      Return
      End

      subroutine xG96_(idord,rho_s,gamma_s,B88,dB88dr,dB88dg,d2B88dr2,
     &                d2B88drdg,d2B88dg2)
      implicit real*8 (a-h,o-z)
      parameter(third=1.0d0/3.0d0)
      parameter(four3=4.0d0/3.0d0)
      parameter(seven3=7.0d0/3.0d0)
C      parameter(bcoef=0.0042d0)
C     parameter(xldacff=0.930525736349100025D0)

      parameter(b=1.0d0/137.0d0)

C      rho = min(rho_s , 1.0D-16 )
         rho=rho_s
C      gamma = min(gamma_s, 1.0D-16 )
        gamma=gamma_s
      r43 = rho**four3
      rhoinv=1.0d0/rho
* lda part:
C     xlda=-xldacff*r43
* Note: Use x=sqrt(gamma)/rho**four3
      x = sqrt(gamma)/r43

      hgi = 0.50d0 / gamma

C Gill-1996 Exchange: Peter Gill, Mol Phys 1996 vol 89, 433
C      B88 = - b * x ** (3.0d0 / 2.0d0)  * r43

      B88 =  - b * x * sqrt(x)  * r43

      if(idord.lt.1) goto 99

      dB88Dr = 2.d+0*b*rho**((-5.d+0)/3.d+0)*gamma**(3.d+0/4.d+0)/3.d+0
      dB88Dg = (-3.d+0)*b*rho**((-2.d+0)/3.d+0)*gamma**((-1.d+0)/4.d+0)/
     1   4.d+0


      if(idord.lt.2) goto 99


      d2B88Dr2 = (-1.d+1)*b*rho**((-8.d+0)/3.d+0)*gamma**(3.d+0/4.d+0)/9
     1   .d+0
      d2B88DrDg = b*rho**((-5.d+0)/3.d+0)*gamma**((-1.d+0)/4.d+0
     2   )/2.d+0
      d2B88Dg2 = 3.d+0*b*rho**((-2.d+0)/3.d+0)*gamma**((-5.d+0
     3   )/4.d+0)/1.6d+1


  99  continue

      return
      end
