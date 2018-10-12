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

#include "warnings.h"

#include "molcas_system.h"

void molcas_sighandler (int signo) {
        INT rc = signo;

#ifdef _PAUSE_ON_SIGNAL_
        pause();
#endif

        switch (signo) {
        case SIGALRM:
                rc = _RC_TIMEOUT_;
                write_rc(&rc);
                printf("Maximum execution time reached\n");
                exit(signo);
                break;
        case SIGINT:
                write_rc(&rc);
                exit(signo);
                break;
        default:
                /* if we are set to catch any other signals,
                   just reset the default handler and reraise */
                write_rc(&rc);
                signal(signo, SIG_DFL);
                raise(signo);
        }
}
