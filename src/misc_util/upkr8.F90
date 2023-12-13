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
! Copyright (C) 1993,2000, Markus P. Fuelscher                         *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine UPKR8(iOpt,nData,nByte,InBuf,OutBuf)
!***********************************************************************
!                                                                      *
!     purpose: decode pack double precision floating point numbers     *
!                                                                      *
!    Calling parameters:                                               *
!    iOpt  : Bit switch                                                *
!            bits 0-3 used for packing algorithm                       *
!    nData : number of unpacked double precision floating point numbers*
!    nByte : length of pack array in bytes                             *
!    OutBuf: contains on output unpacked numbers                       *
!    InBuf : contains packed integral on input                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P. O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     replacement of packing routines                                  *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 2000                                 *
!                                                                      *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_loc
use Pack_mod, only: Init_do_setup_d, isPack, PkThrs
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp), intent(inout) :: nData
integer(kind=iwp), intent(out) :: nByte
real(kind=wp), intent(in) :: InBuf(*)
real(kind=wp), intent(_OUT_) :: OutBuf(*)
integer(kind=iwp) :: Kase, nComp
interface
  subroutine rld_r8(in_,n_in,out_,n_out,thr) bind(C,name='rld_r8_')
    use Definitions, only: MOLCAS_C_INT, MOLCAS_C_REAL
    real(kind=MOLCAS_C_REAL) :: in_(*), out_(*), thr
    integer(kind=MOLCAS_C_INT) :: n_in, n_out
  end subroutine rld_r8
end interface

!----------------------------------------------------------------------*
! If data have not been packed copy them from                          *
! the input buffer to the output buffer                                *
!----------------------------------------------------------------------*

if (.not. isPack) then
  OutBuf(1:nData) = InBuf(1:nData)
  nByte = 8*nData
  return
end if

!----------------------------------------------------------------------*
! Unpack the data and transfer the result from                         *
! the input buffer to the output buffer                                *
!----------------------------------------------------------------------*

!call upkERI(nData,PkThrs,nByte,InBuf,OutBuf)
Kase = ibits(iOpt,0,4)
if (Kase == 0) then
  call tcd_r8_wrapper(InBuf,nComp,OutBuf,nData,PkThrs,Init_do_setup_d)
  Init_do_setup_d = 0
  nByte = nComp
else
  call rld_r8(InBuf,nComp,OutBuf,nData,PkThrs)
  nByte = 8*nComp
end if

!----------------------------------------------------------------------*

return

contains

subroutine tcd_r8_wrapper(in_,n_in,out_,n_out,thr,Init_do_setup_d)

  real(kind=wp), target, intent(in) :: in_(*)
  integer(kind=iwp), intent(out) :: n_in
  real(kind=wp), intent(_OUT_) :: out_(*)
  integer(kind=iwp), intent(in) :: n_out, Init_do_setup_d
  real(kind=wp), intent(in) :: thr
  interface
    subroutine tcd_r8(in_,n_in,out_,n_out,thr,Init_do_setup_d) bind(C,name='tcd_r8_')
      use, intrinsic :: iso_c_binding, only: c_ptr
      use Definitions, only: MOLCAS_C_INT, MOLCAS_C_REAL
      type(c_ptr), value :: in_
      integer(kind=MOLCAS_C_INT) :: n_in, n_out, Init_do_setup_d
      real(kind=MOLCAS_C_REAL) :: out_(*), thr
    end subroutine
  end interface

  call tcd_r8(c_loc(in_),n_in,out_,n_out,thr,Init_do_setup_d)

end subroutine tcd_r8_wrapper

end subroutine UPKR8
