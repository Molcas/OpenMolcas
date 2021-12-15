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
*               2013, Victor P. Vysotskiy                              *
***********************************************************************/

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:38:55 stevenv"
 *
 *       File:         parnell.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell - a program for handling copying and deleting of files
 *                 in a parallel environment
 *
 *       History:
 *       V.P. Vysotskiy, University of Lund, Sweden, 2013
 *       Modified for a bunch of subsequent calls
 *
 */

#include "parnell.h"
/* initialize global variables */
char  MyWorkDir[FILENAME_MAX];
char  WorkDir[FILENAME_MAX];

parnell_status_t
parnell (int argc, char * argv[])
{
        /* IMPORTANT NOTE: argc and argv are passed from parnell_cmd,
         * and as such you cannot rely on argv being NULL-terminated
         * at position argc, i.e. argv[argc] != NULL. */

        parnell_status_t status = PARNELL_START;

        char task;

        /* The MPI specification does not require a specific implementation
         * to distribute the command line arguments to the slaves. However,
         * this is almost invariably the case so we rely on that here. The
         * older code copied the arguments into a buffer first and then
         * broadcasted the result to the slaves. This new code completely
         * relies on argc/argv being available to all MPI processes. */

        if (argc < 2) {
                fprintf (stderr, "parnell: no arguments, exiting");
                status = PARNELL_ERROR;
        } else {
                task = argv[1][0];
                argc -= 2; argv += 2;
                status = PARNELL_CONTINUE;
        }

        /* decide to continue or stop */
        if (status != PARNELL_CONTINUE) goto error;

        /* when debugging, print what task is going to be executed */
#ifdef _DEBUGPRINT_
        if (MyRank == 0) {
                printf("parnell: %c", task);
                for (int i=0; i<argc; i++) printf(" %s", argv[i]);
                printf("\n");
        }
        fflush(NULL);
#endif

        if (task == 'b') {
                /* create work directory and subdirectories */
                status = parnell_base (argc, argv);
        } else {
                /* first initialize WorkDir, MyWorkDir */
                if (parnell_init() != PARNELL_OK) goto error;

                /* decide on which function to execute */
                switch (task) {
                case 'c' :
                        status = parnell_copy (argc, argv);
                        break;
                case 'x' :
                        parnell_rmlist(*argv);
                        /* failures during removal are ignored */
                        status = PARNELL_OK;
                        break;
                case 'w' :
                        status = parnell_wipe ();
                        break;
                case '!' :
                        status = parnell_exec(argc, argv);
                        break;
                default :
                        fprintf(stderr,"%d parnell: unknown task character '%c'\n", MyRank, task);
                        status = PARNELL_ERROR;
                        break;
                }
        }

 error:
        fflush(NULL);
        return status;
}
