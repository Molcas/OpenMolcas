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

use xc_f03_lib_m, only: XC_CORRELATION, XC_EXCHANGE, xc_f03_aux_func_ids, xc_f03_aux_func_weights, xc_f03_func_end, &
                        xc_f03_func_get_info, xc_f03_func_info_get_kind, xc_f03_func_info_get_n_ext_params, xc_f03_func_info_t, &
                        xc_f03_func_init, xc_f03_func_set_ext_params, xc_f03_func_t, xc_f03_num_aux_funcs
use Constants, only: Zero
use Definitions, only: wp, iwp, LibxcInt

implicit none
private

integer(kind=iwp), parameter :: nFuncs_max = 4
integer(kind=iwp) :: nFuncs = 0
integer(kind=LibxcInt) :: func_id(nFuncs_Max) = 0_LibxcInt
real(kind=wp) :: Coeffs(nFuncs_Max) = Zero
logical(kind=iwp) :: lExtParams = .false.
type(xc_f03_func_t) :: xc_func(nFuncs_Max)      ! xc functional
type(xc_f03_func_info_t) :: xc_info(nFuncs_Max) ! xc functional info
real(kind=wp), allocatable :: FuncExtParams(:,:)

public :: Clear_External_Params, Coeffs, Free_External_Params, func_id, FuncExtParams, Get_External_Params, &
          Initiate_Libxc_Functionals, lExtParams, libxc_functionals, nFuncs, nFuncs_max, Put_External_Params, &
          Remove_Libxc_Functionals, Set_External_Params

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

  integer(kind=iwp), intent(in) :: nD
  integer(kind=iwp) :: iFunc
  real(kind=wp) :: Coeff

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

  integer(kind=iwp) :: iFunc

  do iFunc=1,nFuncs
    call xc_f03_func_end(xc_func(iFunc))
  end do
  Coeffs(:) = Zero
  func_id(:) = 0

end subroutine Remove_Libxc_functionals
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine libxc_functionals(mGrid,nD)

  use nq_Grid, only: F_xc, F_xca, F_xcb, l_casdft, vLapl, vRho, vSigma, vTau

  integer(kind=iwp), intent(in) :: mGrid, nD
  integer(kind=iwp) :: iFunc
  real(kind=wp) :: Coeff

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
subroutine Set_External_Params()

  integer(kind=iwp) :: iFunc

  if (allocated(FuncExtParams)) then
    ! Parameters stored on the RunFile may belong to another functional (e.g. left by a previous run)
    if (size(FuncExtParams,2) < nFuncs) then
      call WarningMessage(2,'External parameters do not match the functional (number of functions)!')
      call Quit_OnUserError()
    end if
    do iFunc=1,nFuncs
      if (xc_f03_func_info_get_n_ext_params(xc_info(iFunc)) > size(FuncExtParams,1)) then
        call WarningMessage(2,'External parameters do not match the functional (number of parameters)!')
        call Quit_OnUserError()
      end if
    end do
    do iFunc=1,nFuncs
      call xc_f03_func_set_ext_params(xc_func(iFunc),FuncExtParams(:,iFunc))
    end do
    !call mma_deallocate(FuncExtParams)
  else
    call WarningMessage(2,'External Parameter Arrays Not Initialized!')
    call Quit_OnUserError()
  end if

end subroutine Set_External_Params
!                                                                      *
!***********************************************************************
!                                                                      *
! The external parameters are read in the MCPDFT module. To make them available to the
! modules that run afterwards (e.g. the grid contribution to the gradient in ALASKA),
! they are kept on the RunFile as [number of functions, parameters per function, values].
subroutine Put_External_Params()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: nF, nP
  real(kind=wp), allocatable :: Flat(:)

  if (.not. allocated(FuncExtParams)) return
  nP = size(FuncExtParams,1)
  nF = size(FuncExtParams,2)
  call mma_allocate(Flat,2+nP*nF,Label='ExtFuncParams')
  Flat(1) = real(nF,kind=wp)
  Flat(2) = real(nP,kind=wp)
  Flat(3:) = reshape(FuncExtParams,[nP*nF])
  call Put_dArray('ExtFuncParams',Flat,size(Flat))
  call mma_deallocate(Flat)

end subroutine Put_External_Params
!                                                                      *
!***********************************************************************
!                                                                      *
! Flag on the RunFile that no external parameters are in use (avoids picking up stale ones)
subroutine Clear_External_Params()

  call Put_dArray('ExtFuncParams',[Zero,Zero],2)

end subroutine Clear_External_Params
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Get_External_Params()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: nF, nP, nFlat
  logical(kind=iwp) :: Found
  real(kind=wp), allocatable :: Flat(:)

  if (lExtParams) return
  call Qpg_dArray('ExtFuncParams',Found,nFlat)
  if ((.not. Found) .or. (nFlat < 2)) return
  call mma_allocate(Flat,nFlat,Label='ExtFuncParams')
  call Get_dArray('ExtFuncParams',Flat,nFlat)
  nF = nint(Flat(1))
  nP = nint(Flat(2))
  if ((nF > 0) .and. (nP > 0)) then
    if (nFlat /= 2+nP*nF) then
      call WarningMessage(2,'Corrupted external functional parameters on the RunFile!')
      call Quit_OnUserError()
    end if
    if (allocated(FuncExtParams)) call mma_deallocate(FuncExtParams)
    call mma_allocate(FuncExtParams,nP,nF,Label='FuncExtParams')
    FuncExtParams(:,:) = reshape(Flat(3:),[nP,nF])
    lExtParams = .true.
  end if
  call mma_deallocate(Flat)

end subroutine Get_External_Params
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Free_External_Params()

  use stdalloc, only: mma_deallocate

  if (allocated(FuncExtParams)) call mma_deallocate(FuncExtParams)
  lExtParams = .false.

end subroutine Free_External_Params
!                                                                      *
!***********************************************************************
!                                                                      *
end module libxc_parameters
