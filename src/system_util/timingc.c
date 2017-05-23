/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Markus P. Fuelscher                              *
***********************************************************************/
/**********************************************************************/
/*                                                                    */
/*    Extract elapse and cpu times                                    */
/*                                                                    */
/*    first version by:                                               */
/*    M.P. Fuelscher                                                  */
/*    University of Lund, Sweden, 1994                                */
/*                                                                    */
/*   (void)timingcinit() : initialization                             */
/*   (void)timingc(double elapse,double usercpu,double syscpu)        */
/*                       : returns elapse, user cpu                   */
/*                            and system cpu times                    */
/*                                                                    */
/*   (INT)iclock()       : returns (user) cpu clock ticks             */
/*   (INT)inc_clock      : returns number of clock ticks per second   */
/*--------------------------------------------------------------------*/
/*                                                                    */
/*    history: none                                                   */
/*                                                                    */
/*--------------------------------------------------------------------*/
#ifdef OLD_LINUX
#include <confname.h>
#else
#ifndef _WIN32_
#include <unistd.h>
#endif
#endif

#ifndef _WIN32_
#include <sys/times.h>
#endif
#include <time.h>
/*
 * #ifdef _CRAY_C90_
 * #include <time.h>
 * #else
 * #include <sys/time.h>
 * #endif
 */

#include "molcastype.h"

#ifdef _CAPITALS_
#define timingcinit TIMINGCINIT
#define timingc TIMINGC
#define iclock ICLOCK
#define inc_clock INC_CLOCK
#else
#ifndef ADD_
#define timingcinit timingcinit_
#define timingc timingc_
#define iclock iclock_
#define inc_clock inc_clock_
#endif
#endif

static double clockticks;

void timingcinit()
{
#ifndef _WIN32_
  clockticks=(double)sysconf((int)_SC_CLK_TCK);
#endif

}
void timingc(double *elapse, double *usercpu, double *syscpu)
{
#ifndef _WIN32_
  struct tms buf;
  *elapse=(double)times(&buf)/clockticks;
  *usercpu=(double)buf.tms_utime/clockticks;
  *syscpu=(double)buf.tms_stime/clockticks;
#endif
  return;
}
INT iclock()
{
#ifndef _WIN32_
  struct tms buf;
  (void)times(&buf);
  return((INT)buf.tms_utime);
#else
  return 0;
#endif
}
INT inc_clock()
{
#ifndef _WIN32_
  return((INT)sysconf(_SC_CLK_TCK));
#else
  return 0;
#endif
}
