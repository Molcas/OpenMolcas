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
!***********************************************************************
!                                                                      *
!   This file contains Fortran front-ends for the C timing routines:   *
!                                                                      *
!   (void)timingcinit()    : initialization                            *
!   (void)timingc(double elapse,double usercpu,double syscpu)          *
!                          : returns elapse, user cpu                  *
!                            and system cpu times                      *
!                                                                      *
!   (INT)iclock()          : returns (user) cpu clock ticks            *
!   (INT)inc_clock         : returns number of clock ticks per second  *
!                                                                      *
!***********************************************************************

subroutine SetTim()
!***********************************************************************
!                                                                      *
!     This subroutine has two entry points. The first, called          *
!     SetTim is used to start a global timer. The second, called       *
!     Timing, returns the cpu time used since the timer has been set   *
!     and, in addition the cpu time used since the timer has been      *
!     called the previous time.                                        *
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

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: elapse, usercpu, syscpu
integer(kind=iwp), external :: inc_clock
#include "SysCtl.fh"

call timingcinit()
call timingc(elapse,usercpu,syscpu)
Heuer(1) = usercpu
Heuer(2) = usercpu
Heuer(3) = elapse
Heuer(4) = elapse
ClkInc = inc_clock()

return

end subroutine SetTim
