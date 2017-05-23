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
      Subroutine xOPT(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object: To compute the GGA term of the OPTX exchange functional.     *
*     Note that it is GGA part only, without the LDA part.             *
*     Expressions derived with Maxima  from the formula taken from     *
*     ref: Handy, Cohen J. Mol.Phys 2001, 99, 403-412                  *
*                                                                      *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. June 2006 (B88 code)         *
*             Grigory Shamov, U of Manitoba, Jan 2009                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),
     &       dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
cGLM     &       F_xca(mGrid),F_xcb(mGrid),tmpB(mGrid)

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

         call OPTX_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
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
         call OPTX_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=rho(ipdRxb,iGrid)
         grdrhob_y=rho(ipdRyb,iGrid)
         grdrhob_z=rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call OPTX_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(Fa+Fb)
cGLM         F_xca(iGrid)=F_xca(iGrid)+Coeff*(Fa)
cGLM         F_xcb(iGrid)=F_xcb(iGrid)+Coeff*(Fb)
         dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+Coeff*dFdrhoa
         dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+Coeff*dFdrhob
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
* Note: For xpbe, dFdgammaab is zero.
         dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+Coeff*dFdgammaaa
         dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+Coeff*dFdgammabb
 210     continue
        end do
cGLM        call DCopy_(mGrid,F_xc,1,tmpB,1)

      end if

      Return
      End



      subroutine OPTX_(idord,rho_s,gamma_s, B88, dB88dr, dB88dg,
     & d2B88dr2, d2B88drdg, d2B88dg2)
C      needs 1.05151 coefficient for the Dirac part!
      implicit None
      integer idord
      real*8 rho_s, gamma_s, B88, dB88dr, dB88dg
      real*8 d2B88dr2, d2B88drdg, d2B88dg2
C
      real*8 rho, gamma, y, one, four
      parameter(one=1.0d0)
      parameter(four=4.0d0)
      parameter(y=0.006d0)
C      parameter(bcoef=1.43169d0) for
C      OPTX coeff moved outside, to allow for O3LYP and O2PLYP and KT3
      real*8 T1,T2,T3,T4,T5,T6,T7,T8,T9,T10

      rho=rho_s+1.0D-16
      gamma=gamma_s+1.0D-16

      T1=one/rho**4
      T2=gamma**2
      T3=   6.0d-3*gamma/rho**(8.0d+0/3.0d+0)+1
      T4=one/T3**2
      T5=one/rho**(2.3d+1/3.0d+0)
      T6=gamma**3
      T7=one/T3**3
      T8=one/rho**5
      T9=one/rho**(2.0d+1/3.0d+0)
      T10=one/T3**4

      B88 = -3.60d-5*T1*T2*T4

      if(idord.lt.1) goto 99

      dB88Dr = 1.44d-4*T8*T2*T4-1.1520d-6*T5*T6*T7
      dB88Dg = 4.32d-7*T9*T2*T7-7.20d-5*T1*gamma*T4


      if(idord.lt.2) goto 99

      d2B88dr2 = -7.20d-4*T2*T4/rho**6+1.34400
     7  00d-5*T6*T7/rho**(2.6d+1/3.0d+0)-5.52960d-8*gamma*
     8  *4*T10/rho**(3.4d+1/3.0d+0)

      d2B88drdg = 2.8800d-4*T
     9  8*gamma*T4-5.18400d-6*T5*T2*T7+2.073600d-8
     =  *T6*T10/rho**(3.1d+1/3.0d+0)

      d2B88dg2 = -7.200d-5*
     ;  T1*T4+1.728d-6*T9*gamma*T7-7.77600d-9*T2*T10/rho**(2
     <  .8d+1/3.0d+0)


 99   continue

      return
      end
