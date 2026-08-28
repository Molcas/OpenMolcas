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

module GSBBD2B_CUDA_INTERFACE

  use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

  implicit none
  private

  public :: LUCIA_GSBBD2B_CUDA_BEGIN, LUCIA_GSBBD2B_CUDA_END, LUCIA_GSBBD2B_CUDA_BEGIN_MAPS, &
            LUCIA_GSBBD2B_CUDA_END_MAPS, LUCIA_GSBBD2B_CUDA_ROUTE, LUCIA_GSBBD2B_CUDA_RELEASE

  interface
    function LUCIA_GSBBD2B_CUDA_BEGIN(SB,CB,RHO2,RHO2S,RHO2A,S2_TERM1,NIA,NIB,NJA,NJB,NORB,IPACK) result(Status) &
        bind(C,name='lucia_gsbbd2b_cuda_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(in) :: SB(*), CB(*)
      real(c_double), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*), S2_TERM1
      integer(c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB, NORB, IPACK
    end function LUCIA_GSBBD2B_CUDA_BEGIN

    function LUCIA_GSBBD2B_CUDA_BEGIN_MAPS(I1,XI1,I3,XI3,NKASTR,NI,NJ) result(Status) &
        bind(C,name='lucia_gsbbd2b_cuda_begin_maps')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      integer(c_int64_t), intent(in) :: I1(*), I3(*)
      real(c_double), intent(in) :: XI1(*), XI3(*)
      integer(c_int64_t), value, intent(in) :: NKASTR, NI, NJ
    end function LUCIA_GSBBD2B_CUDA_BEGIN_MAPS

    function LUCIA_GSBBD2B_CUDA_END_MAPS() result(Status) bind(C,name='lucia_gsbbd2b_cuda_end_maps')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_GSBBD2B_CUDA_END_MAPS

    function LUCIA_GSBBD2B_CUDA_END() result(Status) bind(C,name='lucia_gsbbd2b_cuda_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_GSBBD2B_CUDA_END

    function LUCIA_GSBBD2B_CUDA_ROUTE(X,SB,CB,I1,XI1,I3,XI3,I4,XI4,I2,XI2,NIA,NIB,NJA,NJB,NKASTR,KABOT,LKABTC,NKBSTR, &
        NI,NJ,NK,NL,IKORD,IOFF,JOFF,KOFF,LOFF,NORB,IPACK,S2_ACTIVE) result(Status) bind(C,name='lucia_gsbbd2b_cuda_route')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(out) :: X(*)
      real(c_double), intent(in) :: SB(*), CB(*), XI1(*), XI3(*), XI4(*), XI2(*)
      integer(c_int64_t), intent(in) :: I1(*), I3(*), I4(*), I2(*)
      integer(c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB, NKASTR, KABOT, LKABTC, NKBSTR, NI, NJ, NK, NL, IKORD, &
                                               IOFF, JOFF, KOFF, LOFF, NORB, IPACK, S2_ACTIVE
    end function LUCIA_GSBBD2B_CUDA_ROUTE

    subroutine LUCIA_GSBBD2B_CUDA_RELEASE() bind(C,name='lucia_gsbbd2b_cuda_release')
    end subroutine LUCIA_GSBBD2B_CUDA_RELEASE
end interface

end module GSBBD2B_CUDA_INTERFACE

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GSBBD2B_CUDA_INTERFACE)

#endif
