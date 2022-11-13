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
!***********************************************************************

subroutine R8LEN(iOpt,nData,Buf,iLen)
!***********************************************************************
!                                                                      *
!     Subroutine R8Len                                                 *
!                                                                      *
!     purpose: compute the packed size of a double precision           *
!              floating point number                                   *
!                                                                      *
!    Calling parameters:                                               *
!    iOpt  : Bit switch                                                *
!            bits 0-3 used for packing algorithm                       *
!    nData : number of data elements                                   *
!    Buf   : buffer containing the data                                *
!    iLen  : buffer containing the packed length                       *
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

use Pack_mod, only: Init_do_setup_l, isPack, PkThrs
use Definitions, only: wp, iwp, RtoB

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iOpt, nData
real(kind=wp), intent(in) :: Buf(*)
integer(kind=iwp), intent(_OUT_) :: iLen(*)
integer(kind=iwp) :: i, iZero, Kase

!----------------------------------------------------------------------*

if (.not. isPack) then
  call ICOPY(nData,[RtoB],0,iLen,1)
else
  !call ERIlen(nData,PkThrs,Buf,iLen)
  Kase = ibits(iOpt,0,4)
  if (Kase == 0) then
    call tcl_r8(Buf,nData,iLen,PkThrs,Init_do_setup_l)
    Init_do_setup_l = 0
  else
    iZero = 8
    call iCopy(nData,[8],0,iLen,1)
    do i=1,nData
      if (abs(Buf(i)) < PkThrs) then
        iLen(i) = iZero
        iZero = 0
      else
        !iLen(i) = 8
        iZero = 8
      end if
    end do
  end if
end if

!----------------------------------------------------------------------*

return

end subroutine R8LEN
