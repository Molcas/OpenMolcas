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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "molcas_system.h"

/* WARNING: this function might be called from inside a signal
   handler, so use async-signal safe functions only (e.g. no
   (f)printf!) */

void write_rc(INT* code) {
        int rc;
        int negative;
        int wrc;

        char to_ascii [] = "0123456789";

        char stack [5];
        int stack_pos;

        int fd;

        rc = *code;

        /* we don't play games */
        if (rc > 999) rc = 999;
        if (rc < -99) rc = -99;

        if (rc < 0) {
                negative = 1;
                rc = -rc;
        } else {
                negative = 0;
        }

        stack_pos = sizeof(stack) - 1;

        stack[stack_pos--] = '\n';
        do {
                stack[stack_pos--] = to_ascii[rc % 10];
        } while ((rc = rc / 10) && stack_pos > 0);

        if (negative) {
                stack[stack_pos] = '-';
        } else {
                stack_pos++;
        }

        /* write rc string to file */
        fd = open("rc.local", O_CREAT|O_WRONLY|O_TRUNC|O_SYNC, S_IRUSR|S_IWUSR);
        /* we won't check if write/close fails,
           because there is nothing to do about it */
        /* except for keeping the compiler happy */
        wrc = write(fd, stack + stack_pos, sizeof(stack) - stack_pos);
        (void)wrc;
        close(fd);
}
