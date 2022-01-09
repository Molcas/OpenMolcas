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
      Subroutine xG96(mGrid,Coeff,nD,F_xc)
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

      ! Gill 96 exchange
      integer*4, parameter :: func_id = 107

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

    End Subroutine xG96
#else
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
! Copyright (C) 2006, Per Ake Malmqvist                                *
!               2009, Grigory A. Shamov                                *
!***********************************************************************
      Subroutine xG96(mGrid,Coeff,iSpin,F_xc)
!***********************************************************************
!                                                                      *
! Object Gill-1996 Exchange: Peter Gill, Mol Phys 1996 vol 89, 433     *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
!             University of Lund, SWEDEN. June 2006  (B88 code)        *
!             Grigory A. Shamov, U of Manitoba, Dec 2009               *
!***********************************************************************
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20

       Call DiracX(mGrid,iSpin,F_xc,Coeff)
! IDORD=Order of derivatives to request from XPBE:
      idord=1
!
      Rho_Min=T_X*1.0D-2
!
      if (ispin.eq.1) then
! ispin=1 means spin zero.

! T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=rho(1,iGrid)
         if(rhoa.lt.Ta) goto 110
         sigmaaa=Sigma(1,iGrid)

         call xG96_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(2.0D0*Fa)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
! Note: For xpbe, dFdgammaab is zero.
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
         call xG96_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         sigmabb=Sigma(3,iGrid)
         call xG96_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(Fa+Fb)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
         vRho(2,iGrid)=vRho(2,iGrid)+Coeff*dFdrhob
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
! Note: For xpbe, dFdgammaab is zero.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
         vSigma(3,iGrid)=vSigma(3,iGrid)+Coeff*dFdgammabb
 210     continue
        end do

      end if

      Return
      End

      subroutine xG96_(idord,rho_s,gamma_s,B88,dB88dr,dB88dg,d2B88dr2,d2B88drdg,d2B88dg2)
      implicit real*8 (a-h,o-z)
      parameter(third=1.0d0/3.0d0)
      parameter(four3=4.0d0/3.0d0)
      parameter(seven3=7.0d0/3.0d0)
!      parameter(bcoef=0.0042d0)
!     parameter(xldacff=0.930525736349100025D0)

      parameter(b=1.0d0/137.0d0)

!      rho = min(rho_s , 1.0D-16 )
         rho=rho_s
!      gamma = min(gamma_s, 1.0D-16 )
        gamma=gamma_s
      r43 = rho**four3
! lda part:
!     xlda=-xldacff*r43
! Note: Use x=sqrt(gamma)/rho**four3
      x = sqrt(gamma)/r43

! Gill-1996 Exchange: Peter Gill, Mol Phys 1996 vol 89, 433
!      B88 = - b * x ** (3.0d0 / 2.0d0)  * r43

      B88 =  - b * x * sqrt(x)  * r43

      if(idord.lt.1) goto 99

      dB88Dr = 2.d+0*b*rho**((-5.d+0)/3.d+0)*gamma**(3.d+0/4.d+0)/3.d+0
      dB88Dg = (-3.d+0)*b*rho**((-2.d+0)/3.d+0)*gamma**((-1.d+0)/4.d+0)/4.d+0


      if(idord.lt.2) goto 99


      d2B88Dr2 = (-1.d+1)*b*rho**((-8.d+0)/3.d+0)*gamma**(3.d+0/4.d+0)/9.d+0
      d2B88DrDg = b*rho**((-5.d+0)/3.d+0)*gamma**((-1.d+0)/4.d+0)/2.d+0
      d2B88Dg2 = 3.d+0*b*rho**((-2.d+0)/3.d+0)*gamma**((-5.d+0)/4.d+0)/1.6d+1


  99  continue

      return
      end
#endif
