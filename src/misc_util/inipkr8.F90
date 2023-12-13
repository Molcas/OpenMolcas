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

subroutine INIPKR8(PkAcc,PkMode)
!***********************************************************************
!                                                                      *
!     purpose: initial step for the packing and unpacking              *
!              routines. It has to be called prior to any calls to     *
!              the subroutines PKR8 and UPKR8. This routine determines *
!              the number of bytes which can be stripped of from a     *
!              given double precision number,i.e. the least significant*
!              bytes are replaced by zero) such that the resulting     *
!              absolute error is smaller than a given threshold.       *
!                                                                      *
!    Calling parameters:                                               *
!    PkAcc  : desired packing accuracy                                 *
!    PkMode : This is a flag to turn packing on and off                *
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

use Pack_mod, only: Init_do_setup_d, Init_do_setup_e, Init_do_setup_l, isPack, PkThrs
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: PkAcc
logical(kind=iwp), intent(inout) :: PkMode

!----------------------------------------------------------------------*
! Packing is disabled on the C90 due to the CRAY specific data         *
! representation.                                                      *
!----------------------------------------------------------------------*

#ifdef _CRAY_C90_
PkMode = .false.
#endif

!----------------------------------------------------------------------*
! initialize the packing table                                         *
!----------------------------------------------------------------------*

PkThrs = PkAcc
isPack = PkMode

Init_do_setup_e = 1
Init_do_setup_d = 1
Init_do_setup_l = 1

return

end subroutine INIPKR8
