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

module ABTOR2_CUDA_INTERFACE

use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

implicit none
private

public :: LUCIA_ABTOR2_CUDA_ROUTE, LUCIA_ABTOR2_CUDA_RELEASE

interface
  function LUCIA_ABTOR2_CUDA_ROUTE(RHO2B,SKII,CKJJ,NKA,NKB,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,NIB,NJB) &
      result(Status) bind(C,name='lucia_abtor2_cuda_route')
    import :: c_double, c_int64_t
    integer(c_int64_t) :: Status
    real(c_double), intent(inout) :: RHO2B(*)
    real(c_double), intent(in) :: SKII(*), CKJJ(*), XKBIB(*), XKBJB(*)
    integer(c_int64_t), intent(in) :: KBIB(*), KBJB(*)
    integer(c_int64_t), value, intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, IKORD, NIB, NJB
  end function LUCIA_ABTOR2_CUDA_ROUTE

  subroutine LUCIA_ABTOR2_CUDA_RELEASE() bind(C,name='lucia_abtor2_cuda_release')
  end subroutine LUCIA_ABTOR2_CUDA_RELEASE
end interface

end module ABTOR2_CUDA_INTERFACE

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(ABTOR2_CUDA_INTERFACE)

#endif
