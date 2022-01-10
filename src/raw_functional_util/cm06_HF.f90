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
      Subroutine cM06_HF(mGrid,Coeff,nD,F_xc)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, Sigma, Tau, Lapl
      use nq_Grid, only: vSigma
      use nq_Grid, only: vTau
      use libxc
      implicit none
      integer :: mGrid, nD, nRho
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! PBE correlation
      integer*4, parameter :: func_id = 234

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=2.00D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=4.00D0*Sigma(:,1:mGrid)
         If (Allocated(Lapl)) Lapl(:,1:mGrid)=2.00D0*Lapl(:,1:mGrid)
         vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)
         vTau(:,1:mGrid)=2.00D0*vTau(:,1:mGrid)
      Else
         If (Allocated(Tau)) Tau(:,1:mGrid)=0.50D0*Tau(:,1:mGrid)
         vTau(:,1:mGrid)=2.00D0*vTau(:,1:mGrid)
      End If

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=0.50D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=0.25D0*Sigma(:,1:mGrid)
         If (Allocated(Lapl)) Lapl(:,1:mGrid)=0.50D0*Lapl(:,1:mGrid)
         vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
         vTau(:,1:mGrid)=0.50D0*vTau(:,1:mGrid)
      Else
         If (Allocated(Tau)) Tau(:,1:mGrid)=2.00D0*Tau(:,1:mGrid)
         vTau(:,1:mGrid)=0.50D0*vTau(:,1:mGrid)
      End If
      Return

    End Subroutine cM06_HF
#else
     Subroutine CM06_HF(mGrid,CoeffA,iSpin,F_xc)
      Implicit None
      Integer mGrid, iSpin
      Real*8 CoeffA, F_xc(mGrid)
      Integer, parameter :: ijzy=2
      Call CM06(mGrid,CoeffA,iSpin,F_xc,ijzy)
      End Subroutine CM06_HF

#endif
