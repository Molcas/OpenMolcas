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
!               2022, Roland Lindh                                     *
!***********************************************************************

subroutine libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

use xc_f03_lib_m
use nq_Grid, only: Rho, Sigma, Tau, Lapl
use nq_Grid, only: vRho, vSigma, vTau, vLapl
use nq_Grid, only: l_casdft
use nq_Grid, only: F_xca, F_xcb
use libxc
use Definitions, only: wp, iwp, LibxcReal, LibxcSize

implicit none
integer(kind=iwp) :: mGrid, nD, iGrid
real(kind=wp) :: F_xc(mGrid)
real(kind=wp) :: Coeff
type(xc_f03_func_t) :: xc_func      ! xc functional
type(xc_f03_func_info_t) :: xc_info ! xc functional info

if ((LibxcSize /= iwp) .or. (LibxcReal /= wp)) then
  write(6,*) 'Libxc type mismatch!'
  call abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Evaluate energy depending on the family
select case (xc_f03_func_info_get_family(xc_info))
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case (XC_FAMILY_LDA)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    func(1:mGrid) = 0.0d0 ! Initialize memory
    dfunc_drho(:,1:mGrid) = 0.0d0

    if (Only_exc) then
      call xc_f03_lda_exc(xc_func,mGrid,Rho(1,1),func(1))
    else
      call xc_f03_lda_exc_vxc(xc_func,mGrid,Rho(1,1),func(1),dfunc_drho(1,1))
    end if

    ! Libxc evaluates energy density per particle; multiply by
    ! density to get out what we really want
    ! Collect the potential

    if (nD == 1) then
      if (Only_exc) then
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
        end do
      else
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
          vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_drho(1,iGrid)
        end do
      end if
    else
      if (Only_exc) then
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*(Rho(1,iGrid)+Rho(2,iGrid))
        end do
      else
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*(Rho(1,iGrid)+Rho(2,iGrid))
          vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_drho(1,iGrid)
          vRho(2,iGrid) = vRho(2,iGrid)+Coeff*dfunc_drho(2,iGrid)
        end do
      end if

      if (l_casdft) then
        select case (xc_f03_func_info_get_kind(xc_info))
          case (XC_EXCHANGE)
            dFunc_dRho(:,1:mGrid) = Rho(:,1:mGrid)
            Rho(2,1:mGrid) = 0.0d0
            func(1:mGrid) = 0.0d0
            call xc_f03_lda_exc(xc_func,mGrid,Rho(1,1),func(1))
            do iGrid=1,mGrid
              F_xca(iGrid) = F_xca(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
            end do
            Rho(1,1:mGrid) = 0.0d0
            Rho(2,1:mGrid) = dFunc_dRho(2,1:mGrid)
            func(:) = 0.0d0
            call xc_f03_lda_exc(xc_func,mGrid,Rho(1,1),func(1))
            do iGrid=1,mGrid
              F_xcb(iGrid) = F_xcb(iGrid)+Coeff*func(iGrid)*Rho(2,iGrid)
            end do
            Rho(:,1:mGrid) = dFunc_dRho(:,1:mGrid)
        end select
      end if
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    func(1:mGrid) = 0.0d0 ! Initialize memory
    dfunc_drho(:,1:mGrid) = 0.0d0
    dfunc_dSigma(:,1:mGrid) = 0.0d0

    if (Only_exc) then
      call xc_f03_gga_exc(xc_func,mGrid,Rho(1,1),Sigma(1,1),func(1))
    else
      call xc_f03_gga_exc_vxc(xc_func,mGrid,Rho(1,1),Sigma(1,1),func(1),dfunc_dRho(1,1),dfunc_dSigma(1,1))
    end if

    ! Libxc evaluates energy density per particle; multiply by
    ! density to get out what we really want
    ! Collect the potential

    if (nD == 1) then
      if (Only_exc) then
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
        end do
      else
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
          vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_drho(1,iGrid)
          vSigma(1,iGrid) = vSigma(1,iGrid)+Coeff*dfunc_dSigma(1,iGrid)
        end do
      end if
    else
      if (Only_exc) then
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*(Rho(1,iGrid)+Rho(2,iGrid))
        end do
      else
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*(Rho(1,iGrid)+Rho(2,iGrid))
          vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_drho(1,iGrid)
          vRho(2,iGrid) = vRho(2,iGrid)+Coeff*dfunc_drho(2,iGrid)
          vSigma(1,iGrid) = vSigma(1,iGrid)+Coeff*dfunc_dSigma(1,iGrid)
          vSigma(2,iGrid) = vSigma(2,iGrid)+Coeff*dfunc_dSigma(2,iGrid)
          vSigma(3,iGrid) = vSigma(3,iGrid)+Coeff*dfunc_dSigma(3,iGrid)
        end do
      end if

      if (l_casdft) then
        select case (xc_f03_func_info_get_kind(xc_info))
          case (XC_EXCHANGE)
            dFunc_dRho(:,1:mGrid) = Rho(:,1:mGrid)
            Rho(2,1:mGrid) = 0.0d0
            func(1:mGrid) = 0.0d0
            call xc_f03_gga_exc(xc_func,mGrid,Rho(1,1),Sigma(1,1),func(1))
            do iGrid=1,mGrid
              F_xca(iGrid) = F_xca(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
            end do
            Rho(1,1:mGrid) = 0.0d0
            Rho(2,1:mGrid) = dFunc_dRho(2,:)
            func(1:mGrid) = 0.0d0
            call xc_f03_gga_exc(xc_func,mGrid,Rho(1,1),Sigma(1,1),func(1))
            do iGrid=1,mGrid
              F_xcb(iGrid) = F_xcb(iGrid)+Coeff*func(iGrid)*Rho(2,iGrid)
            end do
            Rho(:,1:mGrid) = dFunc_dRho(:,1:mGrid)
        end select
      end if
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (XC_FAMILY_MGGA,XC_FAMILY_HYB_MGGA)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    func(1:mGrid) = 0.0d0 ! Initialize memory
    dfunc_drho(:,1:mGrid) = 0.0d0
    dfunc_dSigma(:,1:mGrid) = 0.0d0
    if (allocated(Tau)) dfunc_dTau(:,1:mGrid) = 0.0d0
    if (allocated(Lapl)) dfunc_dLapl(:,1:mGrid) = 0.0d0

    if (Only_exc) then
      call xc_f03_mgga_exc(xc_func,mGrid,Rho(1,1),Sigma(1,1),Lapl(1,1),Tau(1,1),func(1))
    else
      call xc_f03_mgga_exc_vxc(xc_func,mGrid,Rho(1,1),Sigma(1,1),Lapl(1,1),Tau(1,1),func(1),dfunc_dRho(1,1),dfunc_dSigma(1,1), &
                               dfunc_dLapl(1,1),dfunc_dTau(1,1))
    end if

    ! Libxc evaluates energy density per particle; multiply by
    ! density to get out what we really want
    ! Collect the potential

    if (nD == 1) then
      if (Only_exc) then
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
        end do
      else
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
          vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_drho(1,iGrid)
          vSigma(1,iGrid) = vSigma(1,iGrid)+Coeff*dfunc_dSigma(1,iGrid)
        end do
        if (allocated(Tau)) then
          do iGrid=1,mGrid
            vTau(1,iGrid) = vTau(1,iGrid)+Coeff*dfunc_dTau(1,iGrid)
          end do
        end if
        if (allocated(Lapl)) then
          do iGrid=1,mGrid
            vLapl(1,iGrid) = vLapl(1,iGrid)+Coeff*dfunc_dLapl(1,iGrid)
          end do
        end if
      end if
    else
      if (Only_exc) then
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*(Rho(1,iGrid)+Rho(2,iGrid))
        end do
      else
        do iGrid=1,mGrid
          F_xc(iGrid) = F_xc(iGrid)+Coeff*func(iGrid)*(Rho(1,iGrid)+Rho(2,iGrid))
          vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_drho(1,iGrid)
          vRho(2,iGrid) = vRho(2,iGrid)+Coeff*dfunc_drho(2,iGrid)
          vSigma(1,iGrid) = vSigma(1,iGrid)+Coeff*dfunc_dSigma(1,iGrid)
          vSigma(2,iGrid) = vSigma(2,iGrid)+Coeff*dfunc_dSigma(2,iGrid)
          vSigma(3,iGrid) = vSigma(3,iGrid)+Coeff*dfunc_dSigma(3,iGrid)
        end do
        if (allocated(Tau)) then
          do iGrid=1,mGrid
            vTau(1,iGrid) = vTau(1,iGrid)+Coeff*dfunc_dTau(1,iGrid)
            vTau(2,iGrid) = vTau(2,iGrid)+Coeff*dfunc_dTau(2,iGrid)
          end do
        end if
        if (allocated(Lapl)) then
          do iGrid=1,mGrid
            vLapl(1,iGrid) = vLapl(1,iGrid)+Coeff*dfunc_dLapl(1,iGrid)
            vLapl(2,iGrid) = vLapl(2,iGrid)+Coeff*dfunc_dLapl(2,iGrid)
          end do
        end if
      end if

      if (l_casdft) then
        write(6,*) 'Uncharted territory!'
        call Abend()
        select case (xc_f03_func_info_get_kind(xc_info))
          case (XC_EXCHANGE)
            dFunc_dRho(:,1:mGrid) = Rho(:,1:mGrid)
            Rho(2,1:mGrid) = 0.0d0
            func(1:mGrid) = 0.0d0
            call xc_f03_mgga_exc(xc_func,mGrid,Rho(1,1),Sigma(1,1),Lapl(1,1),Tau(1,1),func(1))
            do iGrid=1,mGrid
              F_xca(iGrid) = F_xca(iGrid)+Coeff*func(iGrid)*Rho(1,iGrid)
            end do
            Rho(1,1:mGrid) = 0.0d0
            Rho(2,1:mGrid) = dFunc_dRho(2,:)
            func(1:mGrid) = 0.0d0
            call xc_f03_mgga_exc(xc_func,mGrid,Rho(1,1),Sigma(1,1),Lapl(1,1),Tau(1,1),func(1))
            do iGrid=1,mGrid
              F_xcb(iGrid) = F_xcb(iGrid)+Coeff*func(iGrid)*Rho(2,iGrid)
            end do
            Rho(:,1:mGrid) = dFunc_dRho(:,1:mGrid)
        end select
      end if
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case default
    write(6,*) 'Libxc family not properly identified.'
    call Abend()
    !                                                                  *
    !*******************************************************************
    !                                                                  *
end select
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine libxc_interface
