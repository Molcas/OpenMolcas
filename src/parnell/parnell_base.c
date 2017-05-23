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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:00:16 stevenv"
 *
 *       File:         parnell_base.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_base - create a work directory and then set up a subdirectory
 *                       for each slave under the main work directory called
 *                       tmp_$n, where $n is the rank number of the process.
 *
 */

#include "parnell.h"

parnell_status_t
parnell_base (int argc, char **argv)
{
        struct stat info;

        if (argc != 1)  {
                fprintf(stderr,"%d parnell_base: expecting 1 argument, received %d\n", MyRank, argc);
                return PARNELL_ERROR;
        }

        /* all processes need to create the parent directory */
        strncpy (WorkDir, *argv, FILENAME_MAX);
        if (stat (WorkDir, &info) != 0) {
                if (errno != ENOENT) {
                        perror ("unexpected error while accessing directory");
                        fprintf(stderr,"%d parnell_base: cannot handle problem with %s\n", MyRank, WorkDir);
                        return PARNELL_ERROR;
                } else {
                        if (mkdir (WorkDir, S_IRWXU | S_IRWXG | S_IRWXO ) != 0) {
                                if (errno != EEXIST) {
                                        perror ("while calling mkdir");
                                        fprintf(stderr,"%d parnell_base: cannot make directory %s\n", MyRank, WorkDir);
                                        return PARNELL_ERROR;
                                }
                        } else {
                                /* when debugging, print success */
#ifdef _DEBUG_
                                fprintf(stdout,"> %d parnell_base: successfully created %s\n", MyRank, WorkDir);
#endif
                        }
                }
        } else if (!S_ISDIR(info.st_mode)) {
                fprintf(stderr,"%d parnell_base: %s is not a directory\n", MyRank, WorkDir);
                return PARNELL_ERROR;
        }

        /* only slave processes need to create a subdirectory */
        if (MyRank == 0) {
                strncpy(MyWorkDir, WorkDir, FILENAME_MAX);
        } else {
                snprintf(MyWorkDir, FILENAME_MAX, "%s/tmp_%d", WorkDir, MyRank);
                if (stat (MyWorkDir, &info) != 0) {
                        if (errno != ENOENT) {
                                perror ("unexpected error while accessing directory");
                                fprintf(stderr,"%d parnell_base: cannot handle problem with %s\n", MyRank, MyWorkDir);
                                return PARNELL_ERROR;
                        } else {
                                if (mkdir (MyWorkDir, S_IRWXU | S_IRWXG | S_IRWXO ) != 0) {
                                        if (errno != EEXIST) {
                                                perror ("while calling mkdir");
                                                fprintf(stderr,"%d parnell_base: cannot make directory %s\n", MyRank, MyWorkDir);
                                                return PARNELL_ERROR;
                                        }
                                } else {
#ifdef _DEBUG_
                                        fprintf(stdout,"> %d parnell_base: successfully created %s\n", MyRank, MyWorkDir);
#endif
                                }
                        }
                } else if (!S_ISDIR (info.st_mode)) {
                        fprintf(stderr,"%d parnell_base: %s is not a directory\n", MyRank, MyWorkDir);
                        return PARNELL_ERROR;
                }
        }

        return PARNELL_OK;
}
