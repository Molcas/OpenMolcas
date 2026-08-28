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

#include "compiler_features.h"
#ifdef _CUDA_BLAS_

module GSBBD2A_CUDA_INTERFACE

  use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

  implicit none
  private

  public :: LUCIA_GSBBD2A_CUDA_BEGIN, LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN, LUCIA_GSBBD2A_CUDA_DENSITY_END, &
            LUCIA_GSBBD2A_CUDA_BLOCK_END, LUCIA_GSBBD2A_CUDA_END, LUCIA_GSBBD2A_CUDA_ROUTE, LUCIA_GSBBD2A_CUDA_RELEASE

  interface
    function LUCIA_GSBBD2A_CUDA_BEGIN(X,SB,CB,NROW,NSB,NCB,XCOUNT) result(Status) &
        bind(C,name='lucia_gsbbd2a_cuda_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: X(*)
      real(c_double), intent(in) :: SB(*), CB(*)
      integer(c_int64_t), value, intent(in) :: NROW, NSB, NCB, XCOUNT
    end function LUCIA_GSBBD2A_CUDA_BEGIN

    function LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN(RHO2,RHO2S,RHO2A,NACOB,IPACK) result(Status) &
        bind(C,name='lucia_gsbbd2a_cuda_density_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*)
      integer(c_int64_t), value, intent(in) :: NACOB, IPACK
    end function LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN

    function LUCIA_GSBBD2A_CUDA_DENSITY_END() result(Status) &
        bind(C,name='lucia_gsbbd2a_cuda_density_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_GSBBD2A_CUDA_DENSITY_END

    function LUCIA_GSBBD2A_CUDA_BLOCK_END(NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NACOB,IPACK) result(Status) &
        bind(C,name='lucia_gsbbd2a_cuda_block_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
      integer(c_int64_t), value, intent(in) :: NI, IOFF, NJ, JOFF, NK, KOFF, NL, LOFF, NACOB, IPACK
    end function LUCIA_GSBBD2A_CUDA_BLOCK_END

    function LUCIA_GSBBD2A_CUDA_END() result(Status) bind(C,name='lucia_gsbbd2a_cuda_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_GSBBD2A_CUDA_END

    function LUCIA_GSBBD2A_CUDA_ROUTE(X,SB,CB,I1,XI1S,I2,XI2S,NROW,NSB,NCB,IBOT,NIBTC,MAXK,NKBTC,NI,NK,NJ,NL,NIK,NJL,IKSM,JLSM, &
        FACTOR) &
        result(Status) bind(C,name='lucia_gsbbd2a_cuda_route')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: X(*)
      real(c_double), intent(in) :: SB(*), CB(*), XI1S(*), XI2S(*)
      integer(c_int64_t), intent(in) :: I1(*), I2(*)
      integer(c_int64_t), value, intent(in) :: NROW, NSB, NCB, IBOT, NIBTC, MAXK, NKBTC, NI, NK, NJ, NL, NIK, NJL, IKSM, JLSM
      real(c_double), value, intent(in) :: FACTOR
    end function LUCIA_GSBBD2A_CUDA_ROUTE

    subroutine LUCIA_GSBBD2A_CUDA_RELEASE() bind(C,name='lucia_gsbbd2a_cuda_release')
    end subroutine LUCIA_GSBBD2A_CUDA_RELEASE
end interface

end module GSBBD2A_CUDA_INTERFACE

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GSBBD2A_CUDA_INTERFACE)

#endif
