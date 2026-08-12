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

module RSBB2BN_CUDA_INTERFACE

  use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

  implicit none
  private

  public :: LUCIA_RSBB2BN_CUDA_BEGIN, LUCIA_RSBB2BN_CUDA_FLUSH, LUCIA_RSBB2BN_CUDA_END, LUCIA_RSBB2BN_CUDA_ROUTE, &
            LUCIA_RSBB2BN_CUDA_RELEASE

  interface
    function LUCIA_RSBB2BN_CUDA_BEGIN(SB,CB,NIA,NIB,NJA,NJB) result(Status) bind(C,name='lucia_rsbb2bn_cuda_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*)
      integer(c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB
    end function LUCIA_RSBB2BN_CUDA_BEGIN

    function LUCIA_RSBB2BN_CUDA_FLUSH() result(Status) bind(C,name='lucia_rsbb2bn_cuda_flush')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_RSBB2BN_CUDA_FLUSH

    function LUCIA_RSBB2BN_CUDA_END() result(Status) bind(C,name='lucia_rsbb2bn_cuda_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_RSBB2BN_CUDA_END

    function LUCIA_RSBB2BN_CUDA_ROUTE(SB,CB,XINT,I1,XI1S,I3,XI3S,I4,XI4S,I2,XI2S,NIA,NIB,NJA,NJB,NKASTR,KABOT,LKABTC, &
        NKBSTR,NI,NJ,NK,NL,IKORD) result(Status) bind(C,name='lucia_rsbb2bn_cuda_route')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*), XINT(*), XI1S(*), XI3S(*), XI4S(*), XI2S(*)
      integer(c_int64_t), intent(in) :: I1(*), I3(*), I4(*), I2(*)
      integer(c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB, NKASTR, KABOT, LKABTC, NKBSTR, NI, NJ, NK, NL, IKORD
    end function LUCIA_RSBB2BN_CUDA_ROUTE

    subroutine LUCIA_RSBB2BN_CUDA_RELEASE() bind(C,name='lucia_rsbb2bn_cuda_release')
    end subroutine LUCIA_RSBB2BN_CUDA_RELEASE
  end interface

end module RSBB2BN_CUDA_INTERFACE
