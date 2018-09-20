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
* Copyright (C) 2010, Steven Vancoillie                                *
***********************************************************************/

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:42:14 stevenv"
*
*       File:         parnell_exec.c
*       Author:       Steven Vancoillie
*       Date:         Spring 2010
*
*       parnell_exec - executes a command by calling system() on all nodes
*
*       WARNING: this is actually not guaranteed to work in an MPI environment
*
*/

#include "parnell.h"
#include <sys/wait.h>

parnell_status_t
parnell_exec (int argc, char ** argv)
{
        (void)argc;

        pid_t pid;
        parnell_status_t status = PARNELL_OK;
        if (MyRank == 0 && nProcs > 1) {
                fprintf (stdout, "==> WARNING <==\npossible unsafe operation\n==> WARNING <==\n");
        }
        pid = fork();
        if (pid == 0) {
                int rc = execvp(*argv, argv);
                perror("while calling execvp");
                fprintf(stderr, "%d parnell: failed to execute command, rc = %d!\n", MyRank, rc);
                status = PARNELL_ERROR;
        } else {
                waitpid(pid, NULL, 0);
        }
        return status;
}
