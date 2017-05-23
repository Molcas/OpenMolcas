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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:38:16 stevenv"
 *
 *       File:         parnell_unlink.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_unlink - delete files specified in a colon-separated list from
 *                       the main work directory and all subdirectories
 *
 */

#include "parnell.h"

parnell_status_t
parnell_unlink (char *fpath)
{
        struct stat info, wrk_info;
        parnell_status_t status = PARNELL_OK;

        if (stat (MyWorkDir, &wrk_info) != 0) {
                perror ("cannot stat directory");
                fprintf(stderr,"%d parnell_unlink: cannot get status of work directory %s\n", MyRank, MyWorkDir);
                return PARNELL_ERROR;
        }

        /* check if file is located in the work directory */
        if (stat (dirname(fpath), &info) == 0) {
                if (!S_ISDIR (info.st_mode)) {
                        status = PARNELL_ERROR;
                } else if (info.st_ino != wrk_info.st_ino) {
                        status = PARNELL_ERROR;
                }
        } else {
                perror ("cannot stat directory");
                status = PARNELL_ERROR;
        }

        if (status == PARNELL_ERROR) {
                fprintf(stderr,"%d parnell_unlink: file not in work directory %s\n", MyRank, fpath);
                goto exit;
        }

        /* try to delete file and catch errors but don't act on them */
        if (stat (fpath, &info)) {
                /* if error other than "No such file or directory", report it */
                if (errno != ENOENT) {
                        perror("parnell_unlink: error while calling stat on file");
                        status = PARNELL_ERROR;
                }
        } else {
                if (unlink (fpath)) {
                        perror("parnell_unlink: error trying to delete file");
                        status = PARNELL_ERROR;
                }
        }

 exit:
        return status;
}
