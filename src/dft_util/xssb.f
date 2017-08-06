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
************************************************************************
      Subroutine xSSBSW(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object:  SSB-sw exchange functional by Swart, Sola and Bickelhaupt   *
*           this is the first one, "switching between OPTX and PBE"    *
*          derived with Maxima                                         *
*                                                                      *
*                                                                      *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author: template by Per Ake Malmqvist,                          *
*             University of Lund, SWEDEN. June 2006                    *
*             Grigory A Shamov, 2010                                   *
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

         call xSSBSW_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
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
         call xSSBSW_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=rho(ipdRxb,iGrid)
         grdrhob_y=rho(ipdRyb,iGrid)
         grdrhob_z=rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call xSSBSW_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
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

      subroutine XSSBSW_(idord,rho_s,gamma_s,F,dFdr,dFdg,
     & d2Fdr2,  d2Fdrdg, d2Fdg2)
      implicit none

      integer idord
      real*8 rho_s,gamma_s,F,dFdr,dFdg
      real*8  d2Fdr2,  d2Fdrdg, d2Fdg2
      real*8 rho, gamma, pi
      parameter(pi=3.14159265358979d0)

      real*8 T1, T2, T3, T4, T5, T6, T7,  T8,  T9,  T10
      real*8 T11, T12, T13, T14, T15, T16,  T17, T18
      real*8 T19, T20, T21, T22, T23

      rho=rho_s+1.0D-16
      gamma=gamma_s+1.0D-16


      T1=rho**(4.d+0/3.d+0)
      T2=1/rho**(8.d+0/3.d+0)
      T3=4.1867733411865504d-3*T2*gamma+1.d+0
      T4=1/T3
      T5=1/rho**(1.6d+1/3.d+0)
      T6=gamma**2
      T7=1.0930391066596375d-3*T5*T6+1.d+0
      T8=1/T7
      T9=-2.973605770238684d-3*T2*gamma*T8+3.150500329583
     5   405d-3*T2*gamma*T4+1.05151d+0
      T10=1/rho**(1.9d+1/3.d+0)
      T11=1/T3**2
      T12=1/rho**(1.1000000000000001d+1/3.d+0)
      T13=1/rho**9
      T14= gamma**3
      T15=1/T7**2
      T16=7.929615387303156d-3*T12*gamma*T8-1.73
     8   347594381847d-5*T13*T14*T15-8.401334212222413d-3*T12*gamma*T4+3
     9   .5174482110131294d-5*T10*T6*T11
      T17=rho**(1.d+0/3.d+0)
      T18=1/rho**8
      T19=-2.973605770238684d-3*T2*T8+6.500534789319262d-6*T18*T
     ;   6*T15+3.150500329583405d-3*T2*T4-1.3190430791299237d-5*T5*gamma
     <   *T11
      T20=1/rho**10
      T21=1/T3**3
      T22=1/rho**(1.3999999999999999d+ 1/3.d+0)
      T23=1/T7**3

      F = -9.305257363491001d-1*T1*T9

C dont forget to take sqrt of gamma

      if(idord.lt.1) goto 99

      dFdr = -1
     >   .2407009817988002d+0*T17*T9-9.305257363491001d-1*T1*T16

      dFdg =
     ?   -9.305257363491001d-1*T1*T19

      if(idord.lt.2) goto 99

      d2Fdr2 = -4.135669939329333d-1*T9/
     @   rho**(2.d+0/3.d+0)-2.4814019635976003d+0*T17*T16-9.305257363491
     1   001d-1*T1*(-2.907525642011157d-2*T22*gamma*T8+2.022388601121548
     2   3d-4*T20*T14*T15-2.0210741301837976d-7*gamma**5*T23/rho**(4.599
     3   9999999999996d+1/3.d+0)+3.080489211148218d-2*T22*gamma*T4-3.165
     4   7033899118164d-4*T6*T11/rho**(2.2000000000000003d+1/3.d+0)+7.85
     5   4271146066177d-7*T20*T14*T21)

      d2Fdrdg = -1.2407009817988002d+0*
     6   T17*T19-9.305257363491001d-1*T1*(7.929615387303156d-3*T12*T8-6.
     7   93390377527388d-5*T13*T6*T15+7.579027988189241d-8*gamma**4*T23/
     8   rho**(4.3d+1/3.d+0)-8.401334212222413d-3*T12*T4+1.0552344633039
     9   39d-4*T10*gamma*T11-2.9453516797748164d-7*T13*T6*T21)

      d2Fdg2 =
     =   -9.305257363491001d-1*T1*(1.9501604367957787d-5*T18*gamma*T15-2
     ;   .842135495570965d-8*T14*T23/rho**(4.d+1/3.d+0)-2.63808615825984
     <   75d-5*T5*T11+1.1045068799155562d-7*T18*gamma*T21)

 99   continue

      return
      end
