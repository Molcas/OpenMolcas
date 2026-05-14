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
! Copyright (C) 1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine SetTim()
!***********************************************************************
!                                                                      *
!     This subroutine is used to start a global timer.                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark                                                     *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use UnixInfo, only: Heuer
use Definitions, only: wp

implicit none
real(kind=wp) :: elapse, usercpu, syscpu

call timingcinit()
call timingc(elapse,usercpu,syscpu)
Heuer(1) = usercpu
Heuer(2) = usercpu
Heuer(3) = elapse
Heuer(4) = elapse

end subroutine SetTim
