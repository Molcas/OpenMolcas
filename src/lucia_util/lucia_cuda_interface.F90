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

module LUCIA_CUDA_INTERFACE

use, intrinsic :: iso_c_binding, only: c_double, c_int64_t

implicit none
private

public :: LUCIA_ABTOR2_CUDA_RELEASE, LUCIA_ABTOR2_CUDA_ROUTE, LUCIA_GSBBD1_CUDA_BEGIN, LUCIA_GSBBD1_CUDA_END, &
          LUCIA_GSBBD1_CUDA_MAPS_BEGIN, LUCIA_GSBBD1_CUDA_MAPS_END, LUCIA_GSBBD1_CUDA_RELEASE, LUCIA_GSBBD1_CUDA_ROUTE, &
          LUCIA_GSBBD2A_CUDA_BEGIN, LUCIA_GSBBD2A_CUDA_BLOCK_END, LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN, &
          LUCIA_GSBBD2A_CUDA_DENSITY_END, LUCIA_GSBBD2A_CUDA_END, LUCIA_GSBBD2A_CUDA_RELEASE, LUCIA_GSBBD2A_CUDA_ROUTE, &
          LUCIA_GSBBD2B_CUDA_BEGIN, LUCIA_GSBBD2B_CUDA_BEGIN_MAPS, LUCIA_GSBBD2B_CUDA_END, LUCIA_GSBBD2B_CUDA_END_MAPS, &
          LUCIA_GSBBD2B_CUDA_RELEASE, LUCIA_GSBBD2B_CUDA_ROUTE, LUCIA_RSBB1E_CUDA_BEGIN, LUCIA_RSBB1E_CUDA_END, &
          LUCIA_RSBB1E_CUDA_RELEASE, LUCIA_RSBB1E_CUDA_ROUTE, LUCIA_RSBB2A_CUDA_BEGIN, LUCIA_RSBB2A_CUDA_END, &
          LUCIA_RSBB2A_CUDA_RELEASE, LUCIA_RSBB2A_CUDA_ROUTE, LUCIA_RSBB2A_CUDA_XINT_BEGIN, LUCIA_RSBB2A_CUDA_XINT_END, &
          LUCIA_RSBB2BN_CUDA_BEGIN, LUCIA_RSBB2BN_CUDA_END, LUCIA_RSBB2BN_CUDA_FLUSH, LUCIA_RSBB2BN_CUDA_RELEASE, &
          LUCIA_RSBB2BN_CUDA_ROUTE, LUCIA_SIGMA_CUDA_BLOCKS_BEGIN, LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED, &
          LUCIA_SIGMA_CUDA_BLOCKS_END, LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN, LUCIA_SIGMA_CUDA_BLOCKS_HOST_END, &
          LUCIA_SIGMA_CUDA_BLOCKS_RELEASE, LUCIA_SKICKJ_CUDA_RELEASE, LUCIA_SKICKJ_CUDA_ROUTE3

interface
  function LUCIA_ABTOR2_CUDA_ROUTE(RHO2B,SKII,CKJJ,NKA,NKB,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,NIB,NJB) &
    bind(C,name='lucia_abtor2_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_ABTOR2_CUDA_ROUTE
    real(kind=c_double), intent(inout) :: RHO2B(*)
    real(kind=c_double), intent(in) :: SKII(*), CKJJ(*), XKBIB(*), XKBJB(*)
    integer(kind=c_int64_t), intent(in) :: KBIB(*), KBJB(*)
    integer(kind=c_int64_t), value, intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, IKORD, NIB, NJB
  end function LUCIA_ABTOR2_CUDA_ROUTE

  subroutine LUCIA_ABTOR2_CUDA_RELEASE() bind(C,name='lucia_abtor2_cuda_release')
  end subroutine LUCIA_ABTOR2_CUDA_RELEASE

  function LUCIA_GSBBD1_CUDA_BEGIN(RHO1,SRHO1,SB,CB,NACOB,NROW,NSB,NCB,IDOSRHO1) bind(C,name='lucia_gsbbd1_cuda_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD1_CUDA_BEGIN
    real(kind=c_double), intent(inout) :: RHO1(*), SRHO1(*)
    real(kind=c_double), intent(in) :: SB(*), CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NACOB, NROW, NSB, NCB, IDOSRHO1
  end function LUCIA_GSBBD1_CUDA_BEGIN

  function LUCIA_GSBBD1_CUDA_END() bind(C,name='lucia_gsbbd1_cuda_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD1_CUDA_END
  end function LUCIA_GSBBD1_CUDA_END

  function LUCIA_GSBBD1_CUDA_MAPS_BEGIN(I1,XI1S,I2,XI2S,NKASTR,D1,D2) bind(C,name='lucia_gsbbd1_cuda_maps_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD1_CUDA_MAPS_BEGIN
    integer(kind=c_int64_t), intent(in) :: I1(*), I2(*)
    real(kind=c_double), intent(in) :: XI1S(*), XI2S(*)
    integer(kind=c_int64_t), value, intent(in) :: NKASTR, D1, D2
  end function LUCIA_GSBBD1_CUDA_MAPS_BEGIN

  function LUCIA_GSBBD1_CUDA_MAPS_END() bind(C,name='lucia_gsbbd1_cuda_maps_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD1_CUDA_MAPS_END
  end function LUCIA_GSBBD1_CUDA_MAPS_END

  function LUCIA_GSBBD1_CUDA_ROUTE(RHO1,SRHO1,SB,CB,I1,XI1S,I2,XI2S,NACOB,NROW,NSB,NCB,IBOT,NIBTC,NKASTR,KBOT,LKABTC,D1,D2, &
                                   OFF1,OFF2,IDOSRHO1,XAB) bind(C,name='lucia_gsbbd1_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD1_CUDA_ROUTE
    real(kind=c_double), intent(inout) :: RHO1(*), SRHO1(*)
    real(kind=c_double), intent(in) :: SB(*), CB(*), XI1S(*), XI2S(*)
    integer(kind=c_int64_t), intent(in) :: I1(*), I2(*)
    integer(kind=c_int64_t), value, intent(in) :: NACOB, NROW, NSB, NCB, IBOT, NIBTC, NKASTR, KBOT, LKABTC, D1, D2, OFF1, OFF2, &
                                                  IDOSRHO1
    real(kind=c_double), value, intent(in) :: XAB
  end function LUCIA_GSBBD1_CUDA_ROUTE

  subroutine LUCIA_GSBBD1_CUDA_RELEASE() bind(C,name='lucia_gsbbd1_cuda_release')
  end subroutine LUCIA_GSBBD1_CUDA_RELEASE

  function LUCIA_GSBBD2A_CUDA_BEGIN(X,SB,CB,NROW,NSB,NCB,XCOUNT) bind(C,name='lucia_gsbbd2a_cuda_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2A_CUDA_BEGIN
    real(kind=c_double), intent(inout) :: X(*)
    real(kind=c_double), intent(in) :: SB(*), CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NROW, NSB, NCB, XCOUNT
  end function LUCIA_GSBBD2A_CUDA_BEGIN

  function LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN(RHO2,RHO2S,RHO2A,NACOB,IPACK) bind(C,name='lucia_gsbbd2a_cuda_density_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN
    real(kind=c_double), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*)
    integer(kind=c_int64_t), value, intent(in) :: NACOB, IPACK
  end function LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN

  function LUCIA_GSBBD2A_CUDA_DENSITY_END() bind(C,name='lucia_gsbbd2a_cuda_density_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2A_CUDA_DENSITY_END
  end function LUCIA_GSBBD2A_CUDA_DENSITY_END

  function LUCIA_GSBBD2A_CUDA_BLOCK_END(NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NACOB,IPACK) bind(C,name='lucia_gsbbd2a_cuda_block_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2A_CUDA_BLOCK_END
    integer(kind=c_int64_t), value, intent(in) :: NI, IOFF, NJ, JOFF, NK, KOFF, NL, LOFF, NACOB, IPACK
  end function LUCIA_GSBBD2A_CUDA_BLOCK_END

  function LUCIA_GSBBD2A_CUDA_END() bind(C,name='lucia_gsbbd2a_cuda_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2A_CUDA_END
  end function LUCIA_GSBBD2A_CUDA_END

  function LUCIA_GSBBD2A_CUDA_ROUTE(X,SB,CB,I1,XI1S,I2,XI2S,NROW,NSB,NCB,IBOT,NIBTC,MAXK,NKBTC,NI,NK,NJ,NL,NIK,NJL,IKSM,JLSM, &
                                    FACTOR) bind(C,name='lucia_gsbbd2a_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2A_CUDA_ROUTE
    real(kind=c_double), intent(inout) :: X(*)
    real(kind=c_double), intent(in) :: SB(*), CB(*), XI1S(*), XI2S(*)
    integer(kind=c_int64_t), intent(in) :: I1(*), I2(*)
    integer(kind=c_int64_t), value, intent(in) :: NROW, NSB, NCB, IBOT, NIBTC, MAXK, NKBTC, NI, NK, NJ, NL, NIK, NJL, IKSM, JLSM
    real(kind=c_double), value, intent(in) :: FACTOR
  end function LUCIA_GSBBD2A_CUDA_ROUTE

  subroutine LUCIA_GSBBD2A_CUDA_RELEASE() bind(C,name='lucia_gsbbd2a_cuda_release')
  end subroutine LUCIA_GSBBD2A_CUDA_RELEASE

  function LUCIA_GSBBD2B_CUDA_BEGIN(SB,CB,RHO2,RHO2S,RHO2A,S2_TERM1,NIA,NIB,NJA,NJB,NORB,IPACK) &
    bind(C,name='lucia_gsbbd2b_cuda_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2B_CUDA_BEGIN
    real(kind=c_double), intent(in) :: SB(*), CB(*)
    real(kind=c_double), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*), S2_TERM1
    integer(kind=c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB, NORB, IPACK
  end function LUCIA_GSBBD2B_CUDA_BEGIN

  function LUCIA_GSBBD2B_CUDA_BEGIN_MAPS(I1,XI1,I3,XI3,NKASTR,NI,NJ) bind(C,name='lucia_gsbbd2b_cuda_begin_maps')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2B_CUDA_BEGIN_MAPS
    integer(kind=c_int64_t), intent(in) :: I1(*), I3(*)
    real(kind=c_double), intent(in) :: XI1(*), XI3(*)
    integer(kind=c_int64_t), value, intent(in) :: NKASTR, NI, NJ
  end function LUCIA_GSBBD2B_CUDA_BEGIN_MAPS

  function LUCIA_GSBBD2B_CUDA_END_MAPS() bind(C,name='lucia_gsbbd2b_cuda_end_maps')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2B_CUDA_END_MAPS
  end function LUCIA_GSBBD2B_CUDA_END_MAPS

  function LUCIA_GSBBD2B_CUDA_END() bind(C,name='lucia_gsbbd2b_cuda_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2B_CUDA_END
  end function LUCIA_GSBBD2B_CUDA_END

  function LUCIA_GSBBD2B_CUDA_ROUTE(X,SB,CB,I1,XI1,I3,XI3,I4,XI4,I2,XI2,NIA,NIB,NJA,NJB,NKASTR,KABOT,LKABTC,NKBSTR,NI,NJ,NK,NL, &
                                    IKORD,IOFF,JOFF,KOFF,LOFF,NORB,IPACK,S2_ACTIVE) bind(C,name='lucia_gsbbd2b_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_GSBBD2B_CUDA_ROUTE
    real(kind=c_double), intent(out) :: X(*)
    real(kind=c_double), intent(in) :: SB(*), CB(*), XI1(*), XI3(*), XI4(*), XI2(*)
    integer(kind=c_int64_t), intent(in) :: I1(*), I3(*), I4(*), I2(*)
    integer(kind=c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB, NKASTR, KABOT, LKABTC, NKBSTR, NI, NJ, NK, NL, IKORD, IOFF, &
                                                  JOFF, KOFF, LOFF, NORB, IPACK, S2_ACTIVE
  end function LUCIA_GSBBD2B_CUDA_ROUTE

  subroutine LUCIA_GSBBD2B_CUDA_RELEASE() bind(C,name='lucia_gsbbd2b_cuda_release')
  end subroutine LUCIA_GSBBD2B_CUDA_RELEASE

  function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN(SB,CB,NSB,NCB) bind(C,name='lucia_sigma_cuda_blocks_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_SIGMA_CUDA_BLOCKS_BEGIN
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NSB, NCB
  end function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN

  function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED(SB,CB,NSB,NCB) bind(C,name='lucia_sigma_cuda_blocks_begin_zeroed')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NSB, NCB
  end function LUCIA_SIGMA_CUDA_BLOCKS_BEGIN_ZEROED

  function LUCIA_SIGMA_CUDA_BLOCKS_END() bind(C,name='lucia_sigma_cuda_blocks_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_SIGMA_CUDA_BLOCKS_END
  end function LUCIA_SIGMA_CUDA_BLOCKS_END

  function LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN(SB,CB,NSB,NCB) bind(C,name='lucia_sigma_cuda_blocks_host_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NSB, NCB
  end function LUCIA_SIGMA_CUDA_BLOCKS_HOST_BEGIN

  function LUCIA_SIGMA_CUDA_BLOCKS_HOST_END() bind(C,name='lucia_sigma_cuda_blocks_host_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_SIGMA_CUDA_BLOCKS_HOST_END
  end function LUCIA_SIGMA_CUDA_BLOCKS_HOST_END

  subroutine LUCIA_SIGMA_CUDA_BLOCKS_RELEASE() bind(C,name='lucia_sigma_cuda_blocks_release')
  end subroutine LUCIA_SIGMA_CUDA_BLOCKS_RELEASE

  function LUCIA_RSBB1E_CUDA_BEGIN(SB,CB,NROW) bind(C,name='lucia_rsbb1e_cuda_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB1E_CUDA_BEGIN
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NROW
  end function LUCIA_RSBB1E_CUDA_BEGIN

  function LUCIA_RSBB1E_CUDA_END() bind(C,name='lucia_rsbb1e_cuda_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB1E_CUDA_END
  end function LUCIA_RSBB1E_CUDA_END

  function LUCIA_RSBB1E_CUDA_ROUTE(SB,CB,H,I1,XI1,I2,XI2,NROW,NCB,NSB,NKASTR,NKAEFF,D1,D2,MAXK) &
    bind(C,name='lucia_rsbb1e_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB1E_CUDA_ROUTE
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*), H(*), XI1(*), XI2(*)
    integer(kind=c_int64_t), intent(in) :: I1(*), I2(*)
    integer(kind=c_int64_t), value, intent(in) :: NROW, NCB, NSB, NKASTR, NKAEFF, D1, D2, MAXK
  end function LUCIA_RSBB1E_CUDA_ROUTE

  subroutine LUCIA_RSBB1E_CUDA_RELEASE() bind(C,name='lucia_rsbb1e_cuda_release')
  end subroutine LUCIA_RSBB1E_CUDA_RELEASE

  function LUCIA_RSBB2A_CUDA_BEGIN(SB,CB,NROW,NSB,NCB) bind(C,name='lucia_rsbb2a_cuda_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2A_CUDA_BEGIN
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NROW, NSB, NCB
  end function LUCIA_RSBB2A_CUDA_BEGIN

  function LUCIA_RSBB2A_CUDA_END() bind(C,name='lucia_rsbb2a_cuda_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2A_CUDA_END
  end function LUCIA_RSBB2A_CUDA_END

  function LUCIA_RSBB2A_CUDA_XINT_BEGIN(XINT,NIK,NJL) bind(C,name='lucia_rsbb2a_cuda_xint_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2A_CUDA_XINT_BEGIN
    real(kind=c_double), intent(in) :: XINT(*)
    integer(kind=c_int64_t), value, intent(in) :: NIK, NJL
  end function LUCIA_RSBB2A_CUDA_XINT_BEGIN

  function LUCIA_RSBB2A_CUDA_XINT_END() bind(C,name='lucia_rsbb2a_cuda_xint_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2A_CUDA_XINT_END
  end function LUCIA_RSBB2A_CUDA_XINT_END

  function LUCIA_RSBB2A_CUDA_ROUTE(SB,CB,XINT,CMAP,CSIGN,SMAP,SSIGN,NROW,NSB,NCB,IBOT,NIBTC,NKBTC,NIK,NJL,FACTOR) &
    bind(C,name='lucia_rsbb2a_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2A_CUDA_ROUTE
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*), XINT(*), CSIGN(*), SSIGN(*)
    integer(kind=c_int64_t), intent(in) :: CMAP(*), SMAP(*)
    integer(kind=c_int64_t), value, intent(in) :: NROW, NSB, NCB, IBOT, NIBTC, NKBTC, NIK, NJL
    real(kind=c_double), value, intent(in) :: FACTOR
  end function LUCIA_RSBB2A_CUDA_ROUTE

  subroutine LUCIA_RSBB2A_CUDA_RELEASE() bind(C,name='lucia_rsbb2a_cuda_release')
  end subroutine LUCIA_RSBB2A_CUDA_RELEASE

  function LUCIA_RSBB2BN_CUDA_BEGIN(SB,CB,NIA,NIB,NJA,NJB) bind(C,name='lucia_rsbb2bn_cuda_begin')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2BN_CUDA_BEGIN
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*)
    integer(kind=c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB
  end function LUCIA_RSBB2BN_CUDA_BEGIN

  function LUCIA_RSBB2BN_CUDA_FLUSH() bind(C,name='lucia_rsbb2bn_cuda_flush')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2BN_CUDA_FLUSH
  end function LUCIA_RSBB2BN_CUDA_FLUSH

  function LUCIA_RSBB2BN_CUDA_END() bind(C,name='lucia_rsbb2bn_cuda_end')
    import :: c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2BN_CUDA_END
  end function LUCIA_RSBB2BN_CUDA_END

  function LUCIA_RSBB2BN_CUDA_ROUTE(SB,CB,XINT,I1,XI1S,I3,XI3S,I4,XI4S,I2,XI2S,NIA,NIB,NJA,NJB,NKASTR,KABOT,LKABTC,NKBSTR,NI,NJ, &
                                    NK,NL,IKORD) bind(C,name='lucia_rsbb2bn_cuda_route')
    import :: c_double, c_int64_t
    integer(kind=c_int64_t) :: LUCIA_RSBB2BN_CUDA_ROUTE
    real(kind=c_double), intent(inout) :: SB(*)
    real(kind=c_double), intent(in) :: CB(*), XINT(*), XI1S(*), XI3S(*), XI4S(*), XI2S(*)
    integer(kind=c_int64_t), intent(in) :: I1(*), I3(*), I4(*), I2(*)
    integer(kind=c_int64_t), value, intent(in) :: NIA, NIB, NJA, NJB, NKASTR, KABOT, LKABTC, NKBSTR, NI, NJ, NK, NL, IKORD
  end function LUCIA_RSBB2BN_CUDA_ROUTE

  subroutine LUCIA_RSBB2BN_CUDA_RELEASE() bind(C,name='lucia_rsbb2bn_cuda_release')
  end subroutine LUCIA_RSBB2BN_CUDA_RELEASE

  function LUCIA_SKICKJ_CUDA_ROUTE3(SKII,CKJJ,XIJKL,NKA,NKB,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,NIB,NJB,FACS) &
    bind(C,name='lucia_skickj_cuda_route3')
    import :: c_double, c_int64_t
    integer(c_int64_t) :: LUCIA_SKICKJ_CUDA_ROUTE3
    real(kind=c_double), intent(inout) :: SKII(*)
    real(kind=c_double), intent(in) :: CKJJ(*), XIJKL(*), XKBIB(*), XKBJB(*)
    integer(kind=c_int64_t), intent(in) :: KBIB(*), KBJB(*)
    integer(kind=c_int64_t), value, intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, IKORD, NIB, NJB
    real(kind=c_double), value, intent(in) :: FACS
  end function LUCIA_SKICKJ_CUDA_ROUTE3

  subroutine LUCIA_SKICKJ_CUDA_RELEASE() bind(C,name='lucia_skickj_cuda_release')
  end subroutine LUCIA_SKICKJ_CUDA_RELEASE
end interface

end module LUCIA_CUDA_INTERFACE

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(LUCIA_CUDA_INTERFACE)

#endif
