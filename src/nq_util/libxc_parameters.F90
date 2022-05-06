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
! Copyright (C) 2000,2022, Roland Lindh                                *
!               2022, Susi Lehtola                                     *
!***********************************************************************

module libxc_parameters

use xc_f03_lib_m
use Definitions, only: LibxcInt

implicit none
private

public :: nFuncs_max, nFuncs, Coeffs, func_id, xc_func, xc_info, Initiate_Libxc_functionals, Remove_Libxc_functionals, &
          libxc_functionals

integer, parameter :: nFuncs_max = 4
integer :: i
integer :: nFuncs = 0
real*8 :: Coeffs(nFuncs_Max) = [(0.0d0,i=1,nFuncs_Max)]
integer(kind=LibxcInt) :: func_id(nFuncs_Max) = [(0_LibxcInt,i=1,nFuncs_Max)]
type(xc_f03_func_t) :: xc_func(nFuncs_Max) ! xc functional
type(xc_f03_func_info_t) :: xc_info(nFuncs_Max) ! xc functional info

!                                                                      *
!***********************************************************************
!                                                                      *
contains
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Initiate_Libxc_functionals(nD)

  use nq_Grid, only: l_casdft
  use KSDFT_Info, only: CoefR, CoefX

  implicit none
  integer nD, iFunc
  real*8 :: Coeff

  ! if it is a mixed functional and we do MC-PDFT split it up in the components for
  ! further analysis.
  if ((nFuncs == 1) .and. l_casdft) then
    call xc_f03_func_init(xc_func(1),func_id(1),int(nD,kind=LibxcInt))
    nFuncs = max(1,int(xc_f03_num_aux_funcs(xc_func(1))))

    if (nFuncs /= 1) then
      call xc_f03_aux_func_ids(xc_func(1),func_id)
      call xc_f03_aux_func_weights(xc_func(1),Coeffs)
    end if
    call xc_f03_func_end(xc_func(1))

  end if
  do iFunc=1,nFuncs
    ! Initialize libxc functional: nD = 2 means spin-polarized
    call xc_f03_func_init(xc_func(iFunc),func_id(iFunc),int(nD,kind=LibxcInt))
    ! Get the functional's information
    xc_info(iFunc) = xc_f03_func_get_info(xc_func(iFunc))

    ! Reset coefficients according to input

    Coeff = Coeffs(iFunc)
    select case (xc_f03_func_info_get_kind(xc_info(iFunc)))
      case (XC_EXCHANGE)
        Coeff = Coeff*CoefX
      case (XC_CORRELATION)
        Coeff = Coeff*CoefR
    end select
    Coeffs(iFunc) = Coeff

  end do

end subroutine Initiate_Libxc_functionals
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Remove_Libxc_functionals()

  implicit none
  integer iFunc

  do iFunc=1,nFuncs
    call xc_f03_func_end(xc_func(iFunc))
  end do
  Coeffs(:) = 0.0d0
  func_id(:) = 0

end subroutine Remove_Libxc_functionals
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine libxc_functionals(mGrid,nD)

  use nq_Grid, only: F_xc, F_xca, F_xcb, l_casdft
  use nq_Grid, only: vRho, vSigma, vTau, vLapl

  implicit none
  integer mGrid, nD, iFunc
  real*8 Coeff
  real*8, parameter :: Zero = 0.0d0

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  vRho(:,1:mGrid) = Zero
  if (allocated(vSigma)) vSigma(:,1:mGrid) = Zero
  if (allocated(vTau)) vTau(:,1:mGrid) = Zero
  if (allocated(vLapl)) vLapl(:,1:mGrid) = Zero
  F_xc(1:mGrid) = Zero
  if (l_casdft) then
    F_xca(1:mGrid) = Zero
    F_xcb(1:mGrid) = Zero
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  do iFunc=1,nFuncs
    Coeff = Coeffs(iFunc)
    call libxc_interface(xc_func(iFunc),xc_info(iFunc),mGrid,nD,F_xc,Coeff)
  end do

  return

end subroutine libxc_functionals
!                                                                      *
!***********************************************************************
!                                                                      *
end module libxc_parameters
