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
      Subroutine xPBE(mGrid,Coeff,nD,F_xc)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vSigma
      use nq_Grid, only: vRho ! to be removed
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
      integer*4, parameter :: func_id = 101

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=2.00D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=4.00D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)

!        F_xc(1:mGrid)=0.0D0
!         Rho(:,1:mGrid)=0.2D0
!        vRho(:,1:mGrid)=0.0D0
!         Sigma(:,1:mGrid)=0.04D0
!        vSigma(:,1:mGrid)=0.0D0
      End If

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=0.50D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=0.25D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
      End If

!     Call RecPrt('F_xc(new)',' ',F_xc,1,mGrid)
!     Call RecPrt('vRho(new)',' ',vRho,SIZE(vRho,1),mGrid)
!     Call RecPrt('vSigma(new)',' ',vSigma,SIZE(vSigma,1),mGrid)
!     Stop 123
      Return

    End Subroutine xPBE
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
!***********************************************************************
      Subroutine XPBE(mGrid,Coeff,iSpin,F_xc)
!***********************************************************************
!                                                                      *
! Object: To compute the functional called x_pbe in the Density        *
! Functional Repository (http://www.cse.clrc.ac.uk/qcg/dft)            *
! Original reference article: J.P. Perdew, K. Burke, and M. Ernzerhof, *
!  "Generalized gradient approximation made simple,"                   *
!      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
!             University of Lund, SWEDEN. June 2006                    *
!***********************************************************************
      use KSDFT_Info, only: F_xca, F_xcb
      use nq_Grid, only: Rho, Sigma, l_casdft
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
! Call arguments:
! Weights(mGrid) (input) integration weights.
! Rho(nRho,mGrid) (input) Density and density derivative values,
!   Rho(1,iGrid) is rho_alpha values, Rho(2,iGrid) is rho_beta values
!   Rho(i,iGrid) is grad_rho_alpha (i=3..5 for d/dx, d/dy, d/dz)
!   Rho(i,iGrid) is grad_rho_beta  (i=6..8 for d/dx, d/dy, d/dz)
! dF_dRho (inout) are (I believe) values of derivatives of the
!   DFT functional (*NOT* derivatives of Fock matrix contributions).
! F_xc is values of the DFT energy density functional (surprised?)

! IDORD=Order of derivatives to request from XPBE:
      idord=1


      if (ispin.eq.1) then
! ispin=1 means spin zero.
! T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(1,iGrid))
         if(rhoa.lt.Ta) goto 110
         sigmaaa=Sigma(1,iGrid)

         call xpbe_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(2.0D0*Fa)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
! Note: For xpbe, dFdgammaab is zero.
 110     continue
        end do
      else
! ispin .ne. 1, use both alpha and beta components.
        If (l_casdft) Then
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(1,iGrid))
         rhob=max(1.0D-24,Rho(2,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) Cycle
         sigmaaa=Sigma(1,iGrid)
         call xpbe_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         sigmabb=Sigma(3,iGrid)
         call xpbe_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc (iGrid)=F_xc (iGrid)+Coeff*(Fa+Fb)
         F_xca(iGrid)=F_xca(iGrid)+Coeff*(Fa)
         F_xcb(iGrid)=F_xcb(iGrid)+Coeff*(   Fb)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
         vRho(2,iGrid)=vRho(2,iGrid)+Coeff*dFdrhob
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
! Note: For xpbe, dFdgammaab is zero.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
         vSigma(3,iGrid)=vSigma(3,iGrid)+Coeff*dFdgammabb
        end do
        Else
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(1,iGrid))
         rhob=max(1.0D-24,Rho(2,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) Cycle
         sigmaaa=Sigma(1,iGrid)
         call xpbe_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         sigmabb=Sigma(3,iGrid)
         call xpbe_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc (iGrid)=F_xc (iGrid)+Coeff*(Fa+Fb)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
         vRho(2,iGrid)=vRho(2,iGrid)+Coeff*dFdrhob
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
! Note: For xpbe, dFdgammaab is zero.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
         vSigma(3,iGrid)=vSigma(3,iGrid)+Coeff*dFdgammabb
        end do
        End If

      end if
*     Call RecPrt('F_xc',' ',F_xc,1,mGrid)
*     Call RecPrt('vRho',' ',vRho,SIZE(vRho,1),mGrid)
*     Call RecPrt('vSigma',' ',vSigma,SIZE(vSigma,1),mGrid)
*     Stop 123

      Return
      End
#endif
