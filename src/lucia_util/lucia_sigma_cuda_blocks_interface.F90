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
! Copyright (C) 2026, Meng Wang                                        *
!***********************************************************************


module LUCIA_SIGMA_CUDA_BLOCKS_INTERFACE

  use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

  implicit none
  private

  public :: LUCIA_SIGMA_CUDA_BLOCKS_BEGIN, LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED, &
            LUCIA_SIGMA_CUDA_BLOCKS_END, LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN, &
            LUCIA_SIGMA_CUDA_BLOCKS_HOST_END, LUCIA_SIGMA_CUDA_BLOCKS_RELEASE

  interface
    function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN(SB,CB,NSB,NCB) result(Status) &
        bind(C,name='lucia_sigma_cuda_blocks_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*)
      integer(c_int64_t), value, intent(in) :: NSB, NCB
    end function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN

    function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED(SB,CB,NSB,NCB) result(Status) &
        bind(C,name='lucia_sigma_cuda_blocks_begin_zeroed')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*)
      integer(c_int64_t), value, intent(in) :: NSB, NCB
    end function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED

    function LUCIA_SIGMA_CUDA_BLOCKS_END() result(Status) &
        bind(C,name='lucia_sigma_cuda_blocks_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_SIGMA_CUDA_BLOCKS_END

    function LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN(SB,CB,NSB,NCB) result(Status) &
        bind(C,name='lucia_sigma_cuda_blocks_host_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*)
      integer(c_int64_t), value, intent(in) :: NSB, NCB
    end function LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN

    function LUCIA_SIGMA_CUDA_BLOCKS_HOST_END() result(Status) &
        bind(C,name='lucia_sigma_cuda_blocks_host_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_SIGMA_CUDA_BLOCKS_HOST_END

    subroutine LUCIA_SIGMA_CUDA_BLOCKS_RELEASE() &
        bind(C,name='lucia_sigma_cuda_blocks_release')
    end subroutine LUCIA_SIGMA_CUDA_BLOCKS_RELEASE
  end interface

end module LUCIA_SIGMA_CUDA_BLOCKS_INTERFACE
