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

subroutine PKR8(iOpt,nData,nByte,InBuf,OutBuf)
!***********************************************************************
!                                                                      *
!     purpose: pack the incoming buffer of double precision floating   *
!              point number into the outgoing buffer.                  *
!                                                                      *
!    Calling parameters:                                               *
!    iOpt  : Bit switch                                                *
!            bits 0-3 used for packing algorithm                       *
!    nData : number of unpacked integrals                              *
!    nByte : length of pack array in bytes                             *
!    InBuf : on input unpacked double precision floating point numbers *
!            (the content of this buffer is destroyed on output!)      *
!    OutBuf: contains packed numbers  on output                        *
!                                                                      *
!    Local data declarations:                                          *
!    Round : Table of constant used for rounding                       *
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
use Pack_mod, only: Init_do_setup_e, isPack, PkThrs
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iOpt, nData
integer(kind=iwp), intent(out) :: nByte
real(kind=wp), intent(in) :: InBuf(*)
real(kind=wp), intent(_OUT_) :: OutBuf(*)
integer(kind=iwp) :: Kase, nComp

!----------------------------------------------------------------------*
! If no packing is desired copy the data from                          *
! the input buffer to the output buffer                                *
!----------------------------------------------------------------------*

if (.not. isPack) then
  OutBuf(1:nData) = InBuf(1:nData)
  nByte = 8*nData
  return
end if

!----------------------------------------------------------------------*
! Pack the data and transfer the result from                           *
! the input buffer to the output buffer                                *
!----------------------------------------------------------------------*

!call pkERI(nData,PkThrs,nByte,InBuf,OutBuf)
Kase = ibits(iOpt,0,4)
if (Kase == 0) then
  call tce_r8_wrapper(InBuf,nData,OutBuf,nComp,PkThrs,Init_do_setup_e)
  Init_do_setup_e = 0
  nByte = nComp
else
  call rle_r8(InBuf,nData,OutBuf,nComp,PkThrs)
  nByte = 8*nComp
end if

!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

return

contains

subroutine tce_r8_wrapper(in_,n_in,out_,n_out,thr,Init_do_setup_e)

  real(kind=wp), intent(in) :: in_(*), thr
  integer(kind=iwp), intent(in) :: n_in, Init_do_setup_e
  real(kind=wp), target, intent(_OUT_) :: out_(*)
  integer(kind=iwp), intent(out) :: n_out
  interface
    subroutine tce_r8(in_,n_in,out_,n_out,thr,Init_do_setup_e) bind(C,name='tce_r8_')
      use, intrinsic :: iso_c_binding, only: c_ptr
      use Definitions, only: MOLCAS_C_INT, MOLCAS_C_REAL
      real(kind=MOLCAS_C_REAL) :: in_(*), thr
      integer(kind=MOLCAS_C_INT) :: n_in, n_out, Init_do_setup_e
      type(c_ptr), value :: out_
    end subroutine tce_r8
  end interface

  call tce_r8(in_,n_in,c_loc(out_),n_out,thr,Init_do_setup_e)

end subroutine tce_r8_wrapper

end subroutine PKR8
