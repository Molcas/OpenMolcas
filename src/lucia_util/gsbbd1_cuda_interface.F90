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

module GSBBD1_CUDA_INTERFACE

  use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

  implicit none
  private

  public :: LUCIA_GSBBD1_CUDA_BEGIN, LUCIA_GSBBD1_CUDA_END, LUCIA_GSBBD1_CUDA_MAPS_BEGIN, &
            LUCIA_GSBBD1_CUDA_MAPS_END, LUCIA_GSBBD1_CUDA_ROUTE, LUCIA_GSBBD1_CUDA_RELEASE

  interface
    function LUCIA_GSBBD1_CUDA_BEGIN(RHO1,SRHO1,SB,CB,NACOB,NROW,NSB,NCB,IDOSRHO1) result(Status) &
        bind(C,name='lucia_gsbbd1_cuda_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: RHO1(*), SRHO1(*)
      real(c_double), intent(in) :: SB(*), CB(*)
      integer(c_int64_t), value, intent(in) :: NACOB, NROW, NSB, NCB, IDOSRHO1
    end function LUCIA_GSBBD1_CUDA_BEGIN

    function LUCIA_GSBBD1_CUDA_END() result(Status) bind(C,name='lucia_gsbbd1_cuda_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_GSBBD1_CUDA_END

    function LUCIA_GSBBD1_CUDA_MAPS_BEGIN(I1,XI1S,I2,XI2S,NKASTR,D1,D2) result(Status) &
        bind(C,name='lucia_gsbbd1_cuda_maps_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      integer(c_int64_t), intent(in) :: I1(*), I2(*)
      real(c_double), intent(in) :: XI1S(*), XI2S(*)
      integer(c_int64_t), value, intent(in) :: NKASTR, D1, D2
    end function LUCIA_GSBBD1_CUDA_MAPS_BEGIN

    function LUCIA_GSBBD1_CUDA_MAPS_END() result(Status) bind(C,name='lucia_gsbbd1_cuda_maps_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_GSBBD1_CUDA_MAPS_END

    function LUCIA_GSBBD1_CUDA_ROUTE(RHO1,SRHO1,SB,CB,I1,XI1S,I2,XI2S,NACOB,NROW,NSB,NCB,IBOT,NIBTC,NKASTR,KBOT,LKABTC,D1,D2, &
        OFF1,OFF2,IDOSRHO1,XAB) result(Status) bind(C,name='lucia_gsbbd1_cuda_route')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: RHO1(*), SRHO1(*)
      real(c_double), intent(in) :: SB(*), CB(*), XI1S(*), XI2S(*)
      integer(c_int64_t), intent(in) :: I1(*), I2(*)
      integer(c_int64_t), value, intent(in) :: NACOB, NROW, NSB, NCB, IBOT, NIBTC, NKASTR, KBOT, LKABTC, D1, D2, OFF1, OFF2, &
                                               IDOSRHO1
      real(c_double), value, intent(in) :: XAB
    end function LUCIA_GSBBD1_CUDA_ROUTE

    subroutine LUCIA_GSBBD1_CUDA_RELEASE() bind(C,name='lucia_gsbbd1_cuda_release')
    end subroutine LUCIA_GSBBD1_CUDA_RELEASE
  end interface

end module GSBBD1_CUDA_INTERFACE
