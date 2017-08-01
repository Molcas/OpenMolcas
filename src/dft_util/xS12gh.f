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
      Subroutine xS12gh(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X,gh_switch)

!      Subroutine xS12g(mGrid,Rho,nRho,P2_ontop,
!     &                nP2_ontop,iSpin,F_xc,
!     &                dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
!     &                T_X)

************************************************************************
*                                                                      *
* Object S12g from Marcel Swart
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
      integer gh_switch
!#include "real.fh"
!#include "nq_index.fh"
!#include "WrkSpc.fh"
!#include "print.fh"
!      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
!     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
!     &       dF_dP2ontop(ndF_dP2ontop,mGrid)


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

         call xS12g_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2,gh_switch)
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

         call xS12g_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2,gh_switch)

         grdrhob_x=rho(ipdRxb,iGrid)
         grdrhob_y=rho(ipdRyb,iGrid)
         grdrhob_z=rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call xS12g_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2,gh_switch)

!         write(6,*) rhoa,rhob,sigmaaa,sigmabb,Fa,Fb

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

!      If (.False.) Then
!         Call Unused_real_array(P2_ontop)
!         Call Unused_real_array(dF_dP2ontop)
!      End If

      Return
      End

      subroutine xS12g_(idord,rho_s,gamma_s,B88,dB88dr,dB88dg,d2B88dr2,
     &                d2B88drdg,d2B88dg2,gh_switch)
      implicit real*8 (a-h,o-z)
      parameter(third=1.0d0/3.0d0)
      parameter(four3=4.0d0/3.0d0)
      parameter(seven3=7.0d0/3.0d0)
C      parameter(bcoef=0.0042d0)
C     parameter(xldacff=0.930525736349100025D0)
      integer gh_switch

      parameter(b=1.0d0/137.0d0)


      if(gh_switch.eq.1) then
* GGA non-hybrid parameter set
        rA = 1.03842032d0
        rK = 0.757d0
        rB = 1d0 + rK - rA
        rA = rA
        rC = 0.00403198d0
        rD = 0.00104596d0
        rE = 0.00594635d0
      elseif(gh_switch.eq.2) then
* GGA hybrid parameter set
        rA = 1.02543951d0
        rK = 0.757d0
        rB = 1d0 + rK - rA
        rA = rA - 0.25d0
        rC = 0.00761554d0
        rD = 0.00211063d0
        rE = 0.00604672d0
      endif


      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)


C      rho = min(rho_s , 1.0D-16 )
         rho=rho_s
C      gamma = min(gamma_s, 1.0D-16 )
      gamma=gamma_s
      rho13 = rho**(1.d0/3.d0)
      rho43 = rho**four3
      rhoinv=1.0d0/rho
* lda part:
C     xlda=-xldacff*rho43
* Note: Use x=sqrt(gamma)/rho**four3
      x = sqrt(gamma)/rho43
      x2 = x*x
      dxdr = -4.d0/3.d0*x/rho


      hgi = 0.50d0 / gamma

      gdenom = 1d0 + rC*x2 + rD*x2*x2
      hdenom = 1d0 + rE*x2
      ums = 1d0 - 1d0 / gdenom
      vms = 1d0 - 1d0 / hdenom
      g = C*rB*ums*vms
c
      dudx = (2d0*rC*x + 4d0*rD*x2*x)/(gdenom**2)
      dvdx = 2d0*rE*x/(hdenom**2)
      dg = C*rB*(dudx*vms + ums*dvdx)



      B88 =  rA*rho43*C
      B88 =  B88 + rho43*g

      if(idord.lt.1) goto 99

      dB88Dr = rA*(4d0/3d0)*rho13*C + (4d0/3d0)*rho13*(g-x*dg)

      t = dg / dsqrt(gamma)
      dB88Dg = t * 0.5d0

      if(idord.lt.2) goto 99

      write(6,*) 'S12g 2nd derivs not programmed'
      Call Abend()
!      d2B88Dr2 = (-1.d+1)*b*rho**((-8.d+0)/3.d+0)*gamma**(3.d+0/4.d+0)/9
!     1   .d+0
!      d2B88DrDg = b*rho**((-5.d+0)/3.d+0)*gamma**((-1.d+0)/4.d+0
!     2   )/2.d+0
!      d2B88Dg2 = 3.d+0*b*rho**((-2.d+0)/3.d+0)*gamma**((-5.d+0
!     3   )/4.d+0)/1.6d+1


  99  continue

      return
      end
