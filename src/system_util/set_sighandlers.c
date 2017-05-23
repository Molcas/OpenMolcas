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
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>

#include "molcas_system.h"

char* getenvc(const char*);

void set_sighandlers (INT* rank) {
        /* SIGALRM handler for timeouts */
        signal(SIGALRM, molcas_sighandler);

        /* install an alarm */
        char* molcas_timelim = getenvc("MOLCAS_TIMELIM"); /* timeout in seconds */
        if(molcas_timelim != NULL) { /* By default, no timelim */
#ifdef _DEMO_
                int timeout_sec = 60;
#else
                int timeout_sec = atoi(molcas_timelim);
#endif
                alarm(timeout_sec);
                if (*rank == 0) {
                        printf("The total execution time is limited to %d seconds.\n", timeout_sec);
                }
                free(molcas_timelim);
        }

        /* interrupt handler needed, as MPI doesn't always respond to
           SIGINT with an immediate quit */
        signal(SIGINT, molcas_sighandler);

#ifdef _GA_
        /* GA overwrites some signal handlers inside
           armci/src/common/signaltrap.c, which is annoying because it
           doesn't usually tell us much about the actual problem. So,
           I reset them here to default handler */
        signal(SIGCHLD, SIG_DFL);
        signal(SIGBUS, SIG_DFL);
        signal(SIGFPE, SIG_DFL);
        signal(SIGILL, SIG_DFL);
        signal(SIGSEGV, SIG_DFL);
        signal(SIGSYS, SIG_DFL);
        signal(SIGTRAP, SIG_DFL);
        signal(SIGHUP, SIG_DFL);
        signal(SIGTERM, SIG_DFL);
        signal(SIGIOT, SIG_DFL);
        signal(SIGCONT, SIG_DFL);
        signal(SIGXCPU, SIG_DFL);
#endif

        return;
}
