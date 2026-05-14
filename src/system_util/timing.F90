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

subroutine Timing(CPUA,CPUE,TIOA,TIOE)
!***********************************************************************
!                                                                      *
!     This subroutine returns the cpu time used since the timer has    *
!     been set and, in addition the cpu time used since the timer has  *
!     been called the previous time.                                   *
!                                                                      *
!     calling arguments:                                               *
!     CPUA     : Type double precision real, output                    *
!                Cpu time since start of the timer                     *
!     CPUE     : Type double precision real, output                    *
!                CPU time since last call to Timing                    *
!     TIOA     : Type double precision real, output                    *
!                Wall clock time since start of timer                  *
!     TIOE     : Type double precision real, output                    *
!                Wall clock time since last call to timing             *
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
real(kind=wp), intent(out) :: CPUA, CPUE, TIOA, TIOE
real(kind=wp) :: elapse, usercpu, syscpu

call timingc(elapse,usercpu,syscpu)
! times relative to saved values
CPUA = usercpu-Heuer(1)
CPUE = usercpu-Heuer(2)
TIOA = elapse-Heuer(3)
TIOE = elapse-Heuer(4)
! save current times
Heuer(2) = usercpu
Heuer(4) = elapse

end subroutine Timing
