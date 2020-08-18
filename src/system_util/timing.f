************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
*   This file contains Fortran front-ends for the C timing routines:   *
*                                                                      *
*   (void)timingcinit()    : initialization                            *
*   (void)timingc(double elapse,double usercpu,double syscpu)          *
*                          : returns elapse, user cpu                  *
*                            and system cpu times                      *
*                                                                      *
*   (INT)iclock()          : returns (user) cpu clock ticks            *
*   (INT)inc_clock         : returns number of clock ticks per second  *
*                                                                      *
************************************************************************
      Subroutine Timing(CPUA,CPUE,TIOA,TIOE)
************************************************************************
*                                                                      *
*     This subroutine has two entry points. The first, called          *
*     SetTim is used to start a global timer. The second, called       *
*     Timing, returns the cpu time used since the timer has been       *
*     and, in addition the cpu time used since the timer has been      *
*     called the previous time.                                        *
*                                                                      *
*     calling arguments:                                               *
*     CPUA     : Type double precision real, output                    *
*                Cpu time since start of the timer                     *
*     CPUE     : Type double precision real, output                    *
*                CPU time since last call to Timing                    *
*     TIOA     : Type double precision real, output                    *
*                Wall clock time since start of timer                  *
*     TIOE     : Type double precision real, output                    *
*                Wall clock time since last call to timing             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark                                                     *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

#include "SysCtl.fh"
*
      Real*8 CPUA,CPUE,TIOA,TIOE
      REal*8 elapse,usercpu,syscpu

      call timingc(elapse,usercpu,syscpu)
* --- times relative to values in Common /$SysBuf/ --------------------*
      CPUA = usercpu-Heuer(1)
      CPUE = usercpu-Heuer(2)
      TIOA = elapse-Heuer(3)
      TIOE = elapse-Heuer(4)
* --- save current times in the Common /$SysBuf/ ----------------------*
      Heuer(2) = usercpu
      Heuer(4) = elapse
      Return
c
      entry SetTim
      call timingcinit()
      call timingc(elapse,usercpu,syscpu)
      Heuer(1) = usercpu
      Heuer(2) = usercpu
      Heuer(3) = elapse
      Heuer(4) = elapse
      ClkInc=inc_clock()
      Return
      End
c
      Function WallCl()
      Real*8 WallCl,elapse,usercpu,syscpu

      call timingc(elapse,usercpu,syscpu)
      WallCl=elapse
      Return
      End
c
      Function Seconds()
      Real*8 Seconds,elapse,usercpu,syscpu

      call timingc(elapse,usercpu,syscpu)
      Seconds=usercpu
      Return
      End
c
      SubRoutine CWTime(usercpu,elapse)
      Real*8 elapse,usercpu,syscpu

      call timingc(elapse,usercpu,syscpu)
      Return
      End
