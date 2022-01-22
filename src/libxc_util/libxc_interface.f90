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
! Copyright (C) 2022, Susi Lehtola                                     *
!                     Roland Lindh                                     *
!***********************************************************************
Subroutine libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)
use xc_f03_lib_m
use nq_Grid, only: Rho,  Sigma,  Tau,  Lapl
use nq_Grid, only:vRho, vSigma, vTau, vLapl
use nq_Grid, only: l_casdft
use nq_Grid, only: F_xca, F_xcb
use libxc
use Definitions, only: wp, iwp, LibxcReal, LibxcSize
implicit none
integer(kind=iwp) :: mGrid, nD, iGrid
Real(kind=wp) :: F_xc(mGrid)
Real(kind=wp) :: Coeff

TYPE(xc_f03_func_t) :: xc_func      ! xc functional
TYPE(xc_f03_func_info_t) :: xc_info ! xc functional info

if ((LibxcSize /= iwp) .or. (LibxcReal /= wp)) then
  write(6,*) 'Libxc type mismatch!'
  call abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Evaluate energy depending on the family
select case (xc_f03_func_info_get_family(xc_info))
!                                                                      *
!***********************************************************************
!                                                                      *
case(XC_FAMILY_LDA)
!                                                                      *
!***********************************************************************
!                                                                      *
   func(1:mGrid) = 0.0D0! Initialize memory
   dfunc_drho(:,1:mGrid) = 0.0D0

   If (Only_exc) Then
      call xc_f03_lda_exc(xc_func, mGrid, Rho(1,1), func(1))
   Else
      call xc_f03_lda_exc_vxc(xc_func, mGrid, Rho(1,1), func(1), dfunc_drho(1,1))
   End If

   ! Libxc evaluates energy density per particle; multiply by
   ! density to get out what we really want
   ! Collect the potential

   If (nD.eq.1) Then
      If (Only_exc) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
         End Do
      End If
   Else
      If (Only_exc) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
            vRho(2,iGrid) = vRho(2,iGrid) + Coeff*dfunc_drho(2, iGrid)
         End Do
      End If

      If (l_casdft) Then
         select case(xc_f03_func_info_get_kind(xc_info))
            case (XC_EXCHANGE);
               dFunc_dRho(:,1:mGrid)=Rho(:,1:mGrid)
               Rho(2,1:mGrid)=0.0D0
               func(1:mGrid)=0.0D0
               call xc_f03_lda_exc(xc_func, mGrid, Rho(1,1), func(1))
               Do iGrid = 1, mGrid
                  F_xca(iGrid) = F_xca(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
               End Do
               Rho(1,1:mGrid)=0.0D0
               Rho(2,1:mGrid)=dFunc_dRho(2,1:mGrid)
               func(:)=0.0D0
               call xc_f03_lda_exc(xc_func, mGrid, Rho(1,1), func(1))
               Do iGrid = 1, mGrid
                  F_xcb(iGrid) = F_xcb(iGrid) + Coeff*func(iGrid)*Rho(2, iGrid)
               End Do
               Rho(:,1:mGrid)=dFunc_dRho(:,1:mGrid)
         end Select
      End If
   End If
!                                                                      *
!***********************************************************************
!                                                                      *
case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
!                                                                      *
!***********************************************************************
!                                                                      *
   func(1:mGrid) = 0.0D0 ! Initialize memory
   dfunc_drho(:,1:mGrid) = 0.0D0
   dfunc_dSigma(:,1:mGrid) = 0.0D0

   If (Only_exc) Then
      call xc_f03_gga_exc(xc_func, mGrid, Rho(1,1), Sigma(1,1), func(1))
   Else
      call xc_f03_gga_exc_vxc(xc_func, mGrid, Rho(1,1), Sigma(1,1), func(1), dfunc_dRho(1,1), dfunc_dSigma(1,1))
   End If

   ! Libxc evaluates energy density per particle; multiply by
   ! density to get out what we really want
   ! Collect the potential

   If (nD.eq.1) Then
      If (Only_exc) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
            vSigma(1,iGrid) = vSigma(1,iGrid) + Coeff*dfunc_dSigma(1, iGrid)
         End Do
      End If
   Else
      If (Only_exc) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
            vRho(2,iGrid) = vRho(2,iGrid) + Coeff*dfunc_drho(2, iGrid)
            vSigma(1,iGrid) = vSigma(1,iGrid) + Coeff*dfunc_dSigma(1, iGrid)
            vSigma(2,iGrid) = vSigma(2,iGrid) + Coeff*dfunc_dSigma(2, iGrid)
            vSigma(3,iGrid) = vSigma(3,iGrid) + Coeff*dfunc_dSigma(3, iGrid)
         End Do
      End If

      If (l_casdft) Then
         select case(xc_f03_func_info_get_kind(xc_info))
            case (XC_EXCHANGE);
               dFunc_dRho(:,1:mGrid)=Rho(:,1:mGrid)
               Rho(2,1:mGrid)=0.0D0
               func(1:mGrid)=0.0D0
               call xc_f03_gga_exc(xc_func, mGrid, Rho(1,1), Sigma(1,1), func(1))
               Do iGrid = 1, mGrid
                  F_xca(iGrid) = F_xca(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
               End Do
               Rho(1,1:mGrid)=0.0D0
               Rho(2,1:mGrid)=dFunc_dRho(2,:)
               func(1:mGrid)=0.0D0
               call xc_f03_gga_exc(xc_func, mGrid, Rho(1,1), Sigma(1,1), func(1))
               Do iGrid = 1, mGrid
                  F_xcb(iGrid) = F_xcb(iGrid) + Coeff*func(iGrid)*Rho(2, iGrid)
               End Do
               Rho(:,1:mGrid)=dFunc_dRho(:,1:mGrid)
         end Select
      End If
   End If
!                                                                      *
!***********************************************************************
!                                                                      *
case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
!                                                                      *
!***********************************************************************
!                                                                      *
   func(1:mGrid) = 0.0D0 ! Initialize memory
   dfunc_drho(:,1:mGrid) = 0.0D0
   dfunc_dSigma(:,1:mGrid) = 0.0D0
   if (Allocated(Tau)) dfunc_dTau(:,1:mGrid) = 0.0D0
   if (Allocated(Lapl)) dfunc_dLapl(:,1:mGrid) = 0.0D0

   If (Only_exc) Then
      call xc_f03_mgga_exc(xc_func, mGrid, Rho(1,1), Sigma(1,1), Lapl(1,1), Tau(1,1), func(1) )
   Else
      call xc_f03_mgga_exc_vxc(xc_func, mGrid, Rho(1,1), Sigma(1,1), Lapl(1,1), Tau(1,1),       &
                               func(1), dfunc_dRho(1,1), dfunc_dSigma(1,1), dfunc_dLapl(1,1), dfunc_dTau(1,1))
   End If

   ! Libxc evaluates energy density per particle; multiply by
   ! density to get out what we really want
   ! Collect the potential

   If (nD.eq.1) Then
      If (Only_exc) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
            vSigma(1,iGrid) = vSigma(1,iGrid) + Coeff*dfunc_dSigma(1, iGrid)
         End Do
         If (Allocated(Tau)) Then
            Do iGrid = 1, mGrid
               vTau(1,iGrid) = vTau(1,iGrid) + Coeff*dfunc_dTau(1, iGrid)
            End Do
         End If
         If (Allocated(Lapl)) Then
            Do iGrid = 1, mGrid
               vLapl(1,iGrid) = vLapl(1,iGrid) + Coeff*dfunc_dLapl(1, iGrid)
            End Do
         End If
      End If
   Else
      If (Only_exc) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
            vRho(2,iGrid) = vRho(2,iGrid) + Coeff*dfunc_drho(2, iGrid)
            vSigma(1,iGrid) = vSigma(1,iGrid) + Coeff*dfunc_dSigma(1, iGrid)
            vSigma(2,iGrid) = vSigma(2,iGrid) + Coeff*dfunc_dSigma(2, iGrid)
            vSigma(3,iGrid) = vSigma(3,iGrid) + Coeff*dfunc_dSigma(3, iGrid)
         End Do
         If (Allocated(Tau)) Then
            Do iGrid = 1, mGrid
               vTau(1,iGrid) = vTau(1,iGrid) + Coeff*dfunc_dTau(1, iGrid)
               vTau(2,iGrid) = vTau(2,iGrid) + Coeff*dfunc_dTau(2, iGrid)
            End Do
         End If
         If (Allocated(Lapl)) Then
            Do iGrid = 1, mGrid
               vLapl(1,iGrid) = vLapl(1,iGrid) + Coeff*dfunc_dLapl(1, iGrid)
               vLapl(2,iGrid) = vLapl(2,iGrid) + Coeff*dfunc_dLapl(2, iGrid)
            End Do
         End If
      End If

      If (l_casdft) Then
         Write (6,*)  "Uncharted territory!"
         Call Abend()
         select case(xc_f03_func_info_get_kind(xc_info))
            case (XC_EXCHANGE);
               dFunc_dRho(:,1:mGrid)=Rho(:,1:mGrid)
               Rho(2,1:mGrid)=0.0D0
               func(1:mGrid)=0.0D0
               call xc_f03_mgga_exc(xc_func, mGrid, Rho(1,1), Sigma(1,1), Lapl(1,1), Tau(1,1), func(1))
               Do iGrid = 1, mGrid
                  F_xca(iGrid) = F_xca(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
               End Do
               Rho(1,1:mGrid)=0.0D0
               Rho(2,1:mGrid)=dFunc_dRho(2,:)
               func(1:mGrid)=0.0D0
               call xc_f03_mgga_exc(xc_func, mGrid, Rho(1,1), Sigma(1,1), Lapl(1,1), Tau(1,1), func(1))
               Do iGrid = 1, mGrid
                  F_xcb(iGrid) = F_xcb(iGrid) + Coeff*func(iGrid)*Rho(2, iGrid)
               End Do
               Rho(:,1:mGrid)=dFunc_dRho(:,1:mGrid)
         end Select
      End If
   End If
!                                                                      *
!***********************************************************************
!                                                                      *
   Case Default
     Write (6,*) "Libxc family not properly identified."
     Call Abend()
!                                                                      *
!***********************************************************************
!                                                                      *
end select
!                                                                      *
!***********************************************************************
!                                                                      *

Return

End Subroutine libxc_interface
