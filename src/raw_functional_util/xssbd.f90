#define _NEWCODE_
#ifdef _NEWCODE_
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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
      Subroutine xSSBD(mGrid,Coeff,nD,F_xc)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vSigma
      use libxc
      implicit none
      integer :: mGrid, nD, nRho
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! pbe exchange
      integer*4, parameter :: func_id = 92

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=2.00D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=4.00D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)
      End If

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=0.50D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=0.25D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
      End If

      Return

    End Subroutine xSSBD
#else
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
*               2017, Giovanni Li Manni                                *
*               2017, Aron Cohen                                       *
************************************************************************
      Subroutine xSSBD(mGrid,
     &                Coeff,iSpin,F_xc)
************************************************************************
*                                                                      *
* Object:  SSB-D exchange functional by Swart                          *
*                                                                      *
*      Author: G. Li Manni, A. Cohen, Max Planck Institute Stuttgart   *
*              Summer 2017, edited in Cambridge (UK) & Palermo (Sicily)*
************************************************************************
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20

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
         rhoa=rho(1,iGrid)
         if(rhoa.lt.Ta) goto 110
         sigmaaa=Sigma(1,iGrid)

         call xSSBD_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(2.0D0*Fa)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
* Note: For xpbe, dFdgammaab is zero.
 110     continue
        end do

      else
* ispin .ne. 1, use both alpha and beta components.

        do iGrid=1,mgrid
         rhoa=Max(Rho_Min,rho(1,iGrid))
         rhob=Max(Rho_Min,rho(2,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         sigmaaa=Sigma(1,iGrid)
         call xSSBD_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         sigmabb=Sigma(3,iGrid)
         call xSSBD_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(Fa+Fb)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
         vRho(2,iGrid)=vRho(2,iGrid)+Coeff*dFdrhob
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
* Note: For xpbe, dFdgammaab is zero.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
         vSigma(3,iGrid)=vSigma(3,iGrid)+Coeff*dFdgammabb
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
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(d2Fdr2)
         Call Unused_real_array(d2Fdrdg)
         Call Unused_real_array(d2Fdg2)
      End If
      end
#endif
