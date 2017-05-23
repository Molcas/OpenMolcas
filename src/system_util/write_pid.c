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
#include <unistd.h>

#include "molcas_system.h"

/* write pid to file */

void write_pid() {
        FILE *stream = fopen("pid", "w");
        fprintf(stream, "%ld\n", (long int)getpid());
        int rc = fclose(stream);
        if (rc) perror("write_pid()");
}
