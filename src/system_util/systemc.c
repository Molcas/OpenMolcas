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
#include <stdio.h>
#include <sys/types.h>
#ifndef _WIN32_
#include <sys/wait.h>
#include <unistd.h>
#endif
#include <signal.h>
#include "molcastype.h"

#ifdef _CAPITALS_
#define systemc SYSTEMC
#define systemc2 SYSTEMC2
#else
#ifndef ADD_
#define systemc systemc_
#define systemc2 systemc2_
#endif
#endif


int system(const char *);

void systemc(char *c, INT *lc, INT *RCC)
   {
#ifdef _WIN32_
 system(c);
#else
   int rc;
   pid_t pid;
   void (*signal_save)();
   c[*lc]='\0';
   signal_save=signal(SIGCHLD,SIG_DFL);
/* puts(c); */
/* puts("here we go"); */
   if (!(pid=fork()))
     execl("/bin/sh", "sh","-c",c,(char *)0);
/* puts("we forked"); */
   waitpid(pid,&rc,0);
/* puts("we get pid"); */
   *RCC=rc;
   signal(SIGCHLD,signal_save);
#endif
/* puts("we getting out"); */
   }
void systemc2(char *c, INT *lc, INT *RCC)
{
   c[*lc]='\0';
   *RCC=system(c);
}
