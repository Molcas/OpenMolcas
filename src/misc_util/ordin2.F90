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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine OrdIn2(iOpt,Buf,lBuf,iBatch)
!***********************************************************************
!                                                                      *
!    Purpose: Read a buffer of ordered two electron integrals          *
!                                                                      *
!    Calling arguments:                                                *
!    Buf    : contains on output the integrals                         *
!    lBuf   : number of integrals to be transfered                     *
!    iOpt   : option code (iOpt=1:start reading at first integral)     *
!                         (iOpt=2:continue reading)                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use TwoDat, only: RAMD
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iOpt, lBuf, iBatch
real(kind=wp), intent(out) :: Buf(lBuf)
integer(kind=iwp) :: iOff

!----------------------------------------------------------------------*
! If this is the first block of a symmetry batch                       *
! get the disk disk start address and load the buffer                  *
!----------------------------------------------------------------------*
if (iOpt == 1) RAMD%next = RAMD%adr(iBatch)
!----------------------------------------------------------------------*
! If the number of requested integrals is larger than                  *
! the current buffer first drain the current buffer and                *
! read as many subsequent buffers as needed                            *
!----------------------------------------------------------------------*
iOff = RAMD%next
Buf(:) = RAMD%ints(iOff:iOff+lBuf-1)
RAMD%next = RAMD%next+lBuf

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine OrdIn2
