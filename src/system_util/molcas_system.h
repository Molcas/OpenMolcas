/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#ifndef _MOLCAS_SYSTEM_H_
#define _MOLCAS_SYSTEM_H_

#include "molcastype.h"

/* functions which interface with Fortran */

#ifdef _CAPITALS_
#  define write_rc WRITE_RC
#  define write_cc WRITE_CC
#  define write_pid WRITE_PID
#  define set_sighandlers  SET_SIGHANDLERS
#else
#  ifndef ADD_
#    define write_rc write_rc_
#    define write_cc write_cc_
#    define write_pid write_pid_
#    define set_sighandlers  set_sighandlers_
#  endif
#endif

void write_rc(INT* code);
void write_cc(INT* code);
void write_pid();
void set_sighandlers (INT* rank);

/* functions only interfacing with C */

void molcas_sighandler (int signo);

#endif
