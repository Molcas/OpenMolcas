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

module SKICKJ_CUDA_INTERFACE

use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

implicit none
private

public :: LUCIA_SKICKJ_CUDA_ROUTE3, LUCIA_SKICKJ_CUDA_RELEASE

interface
  function LUCIA_SKICKJ_CUDA_ROUTE3(SKII,CKJJ,XIJKL,NKA,NKB,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,NIB,NJB,FACS) &
      result(Status) bind(C,name='lucia_skickj_cuda_route3')
    import :: c_double, c_int64_t
    integer(c_int64_t) :: Status
    real(c_double), intent(inout) :: SKII(*)
    real(c_double), intent(in) :: CKJJ(*), XIJKL(*), XKBIB(*), XKBJB(*)
    integer(c_int64_t), intent(in) :: KBIB(*), KBJB(*)
    integer(c_int64_t), value, intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, IKORD, NIB, NJB
    real(c_double), value, intent(in) :: FACS
  end function LUCIA_SKICKJ_CUDA_ROUTE3

  subroutine LUCIA_SKICKJ_CUDA_RELEASE() bind(C,name='lucia_skickj_cuda_release')
  end subroutine LUCIA_SKICKJ_CUDA_RELEASE
end interface

end module SKICKJ_CUDA_INTERFACE
