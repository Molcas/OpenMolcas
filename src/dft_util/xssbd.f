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
*               2010, Grigory A. Shamov                                *
*               2017, G. Li Manni & A. Cohen                           *
************************************************************************
      Subroutine xSSBD(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object:  SSB-D exchange functional by Swart                          *
*                                                                      *
*      Author: G. Li Manni, A. Cohen, Max Planck Institute Stuttgart   *
*              Summer 2017, edited in Cambridge (UK) & Palermo (Sicily)*
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

         call xSSBD_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
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
         call xSSBD_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=rho(ipdRxb,iGrid)
         grdrhob_y=rho(ipdRyb,iGrid)
         grdrhob_z=rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call xSSBD_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
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

      subroutine XSSBD_(idord,rho_s,gamma_s,F,dFdr,dFdg,
     & d2Fdr2,  d2Fdrdg, d2Fdg2)
      implicit none

      integer idord
      real*8 rho_s,gamma_s,F,dFdr,dFdg
      real*8  d2Fdr2,  d2Fdrdg, d2Fdg2
      real*8 rho, gamma, pi
      parameter(pi=3.14159265358979d0)
      real*8 rho43,rrho,rho13,gam12,s,d1s1,d1s2
      real*8 g,gp,d1g1,d1g2
      real*8 rA, rB, rC, rD, rE, rU, RF, fac
      real*8 C, Cs, Cx
      parameter (rA=1.079966d0, rB=0.197465d0, rC=0.272729d0)
      parameter (rE=5.873645d0, rU=-0.749940d0)
      parameter (rD=rB*(1.0d0-rU), rF=0.949488d0)

      rho=rho_s+1.0D-16
      gamma=gamma_s+1.0D-16
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**(1.d0/3.d0)
      Cs = 0.5d0/(3d0*pi*pi)**(1.d0/3.d0)
      Cs = Cs * C               ! account for including C in rho43

      rho43 = C*(2d0*rho)**(4.d0/3.d0)
      rrho = 0.5d0/rho
      rho13 = 4.d0/3.d0*rho43*rrho
      gam12 = 2d0*dsqrt(gamma)
      s = Cs*gam12/rho43
      d1s1 = -4.d0/3.d0*s*rrho
      d1s2 = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c

      g= rB*s*s/(1d0+rC*s*s)
     +               - rD*s*s/(1d0+rE*s**4)
      gp= 2d0*rB*s/(1d0+rC*s*s)**2 +
     +         (2d0*rD*rE*s**5 - 2d0*rD*s)/(1d0+rE*s**4)**2

!            g=gssb0(s)
!            gp=gssb1(s)
c
      Cx = 3.d0/4.d0*(3.d0/pi)**(1.d0/3.d0)*2.d0**(1.d0/3.d0)
      fac = rU*rF*Cx*rB/(2*(3*pi**2)**(1.d0/3.d0))**2/2**(2/3.d0)
      d1g1 = gp*d1s1
      d1g2 = gp*d1s2
      F = rho43*g*0.5d0+rho43*rA*0.5d0
     & +fac*gamma/((rho)**(4.d0/3.d0)+0.1d0)

      if(idord.lt.1) goto 99
            dFdr = (rho13*g+rho43*d1g1)+rho13*rA
     &  - fac*gamma/(rho**(4.d0/3.d0)+0.1d0)**2*4.d0/3.d0*rho**(1/3.d0)
            dFdg =  0.5d0*rho43*d1g2
     &  + fac/((rho)**(4.d0/3.d0)+0.1d0)

      if(idord.lt.2) goto 99

      write(6,*) '2nd derivatives not programmed ssb1'
      Call Abend()

 99   continue

      return
      end

