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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:02:37 stevenv"
 *
 *       File:         parnell_init.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_init - initialize the work directory name as the current
 *                       working directory, and get the subdirectory name for
 *                       each slave called tmp_$n, where $n is the rank number
 *                       of the process, and then change directory to the
 *                       subdirectory.
 *
 *       History:
 *       V.P. Vysotskiy, University of Lund, Sweden, 2013
 *       Modified for a bunch of subsequent calls
 */

#include "parnell.h"

parnell_status_t
parnell_init (void)
{
        char tmpWorkDir[FILENAME_MAX+7];

        /* all process should have been started in the main work directory */
        if (WorkDir[0]==0 && getcwd(WorkDir, FILENAME_MAX)==NULL) {
                perror ("while calling getcwd");
                fprintf(stderr,"%d parnell_init: fatal error, could not determine current working directory\n", MyRank);
                return PARNELL_ERROR;
        }
        /* set MyWorkDir and switch to it */
        if(MyWorkDir[0]==0) {
                if (MyRank == 0) {
                        strncpy (MyWorkDir, WorkDir, FILENAME_MAX);
                } else {
                        snprintf (tmpWorkDir, FILENAME_MAX+7, "%s/tmp_%d", WorkDir, MyRank);
                        strncpy (MyWorkDir, tmpWorkDir, FILENAME_MAX-1);
                        MyWorkDir[FILENAME_MAX-1] = 0;
                        if (chdir (MyWorkDir) != 0) {
                                perror ("cannot change directory");
                                fprintf(stderr,"%d parnell_init: fatal error, could not switch to directory %s\n", MyRank, MyWorkDir);
                                return PARNELL_ERROR;
                        }
                }
        }
        return PARNELL_OK;
}
