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

module RSBB1E_CUDA_INTERFACE

  use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

  implicit none
  private

  public :: LUCIA_RSBB1E_CUDA_BEGIN, LUCIA_RSBB1E_CUDA_END, LUCIA_RSBB1E_CUDA_ROUTE, LUCIA_RSBB1E_CUDA_RELEASE

  interface
    function LUCIA_RSBB1E_CUDA_BEGIN(SB,CB,NROW) result(Status) bind(C,name='lucia_rsbb1e_cuda_begin')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*)
      integer(c_int64_t), value, intent(in) :: NROW
    end function LUCIA_RSBB1E_CUDA_BEGIN

    function LUCIA_RSBB1E_CUDA_END() result(Status) bind(C,name='lucia_rsbb1e_cuda_end')
      import :: c_int64_t
      integer(c_int64_t) :: Status
    end function LUCIA_RSBB1E_CUDA_END

    function LUCIA_RSBB1E_CUDA_ROUTE(SB,CB,H,I1,XI1,I2,XI2,NROW,NCB,NSB,NKASTR,NKAEFF,D1,D2,MAXK) &
        result(Status) bind(C,name='lucia_rsbb1e_cuda_route')
      import :: c_double, c_int64_t
      integer(c_int64_t) :: Status
      real(c_double), intent(inout) :: SB(*)
      real(c_double), intent(in) :: CB(*), H(*), XI1(*), XI2(*)
      integer(c_int64_t), intent(in) :: I1(*), I2(*)
      integer(c_int64_t), value, intent(in) :: NROW, NCB, NSB, NKASTR, NKAEFF, D1, D2, MAXK
    end function LUCIA_RSBB1E_CUDA_ROUTE

    subroutine LUCIA_RSBB1E_CUDA_RELEASE() bind(C,name='lucia_rsbb1e_cuda_release')
    end subroutine LUCIA_RSBB1E_CUDA_RELEASE
end interface

end module RSBB1E_CUDA_INTERFACE

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(RSBB1E_CUDA_INTERFACE)

#endif
