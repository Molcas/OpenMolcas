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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:36:15 stevenv"
 *
 *       File:         parnell_translate.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_translate - get correct relative source and destination names
 *                            for a subsequent copy operation
 *
 */

#include "parnell.h"

/* set proper file names for copy operations, relative to MyWorkDir */
parnell_status_t
parnell_translate (char * src_name, char * dst_name)
{
        parnell_status_t status = PARNELL_START;

        struct stat wrk_info;
        if (stat (MyWorkDir, &wrk_info) != 0) {
                perror ("cannot stat directory");
                fprintf(stderr,"%d parnell_translate: cannot get status of work directory %s\n", MyRank, MyWorkDir);
                return PARNELL_ERROR;
        }

        struct stat src_info, dst_info;
        /* if source exists, it must be a file */
        if (stat (src_name, &src_info) == 0) {
                if (!S_ISREG (src_info.st_mode)) {
                        fprintf(stderr,"%d parnell_translate: not a regular source file %s\n", MyRank, src_name);
                        return PARNELL_ERROR;
                }
        }

        char * dir_name;
        char * tmp_name = (char*)malloc (FILENAME_MAX);

        /* check if source directory is work directory */
        tmp_name[FILENAME_MAX-1]='\0';
        strncpy (tmp_name, src_name, FILENAME_MAX-1);
        dir_name = dirname (tmp_name);
        if (stat (dir_name, &src_info) == 0) {
                if (!S_ISDIR (src_info.st_mode)) {
                        fprintf(stderr,"%d parnell_translate: not a source directory %s\n", MyRank, dir_name);
                        status = PARNELL_ERROR; goto exit;
                } else if (src_info.st_ino == wrk_info.st_ino) {
                        strncpy (tmp_name, src_name, FILENAME_MAX-1);
                        src_name[FILENAME_MAX-1]='\0';
                        strncpy (src_name, basename (tmp_name), FILENAME_MAX-1);
                }
        } else {
                perror ("cannot stat directory");
                fprintf(stderr,"%d parnell_translate: cannot get status of source directory %s\n", MyRank, dir_name);
                status = PARNELL_ERROR; goto exit;
        }

        /* check if destination directory is work directory */
        strncpy (tmp_name, dst_name, FILENAME_MAX-1);
        if (stat (dst_name, &dst_info) == 0) {
                if (S_ISDIR (dst_info.st_mode)) {
                        strncpy (tmp_name, src_name, FILENAME_MAX-1);
                        if (dst_info.st_ino == wrk_info.st_ino) {
                                dst_name[FILENAME_MAX-1]='\0';
                                strncpy (dst_name, basename (tmp_name), FILENAME_MAX-1);
                        } else {
                                strcat (dst_name, "/");
                                strcat (dst_name, basename (tmp_name));
                        }
                } else {
                        strncpy (tmp_name, dst_name, FILENAME_MAX-1);
                        dir_name = dirname (tmp_name);
                        if (stat (dir_name, &dst_info) == 0) {
                                if (S_ISDIR (dst_info.st_mode)) {
                                        if (dst_info.st_ino == wrk_info.st_ino) {
                                                strncpy (tmp_name, dst_name, FILENAME_MAX-1);
                                                dst_name[FILENAME_MAX-1]='\0';
                                                strncpy (dst_name, basename (tmp_name), FILENAME_MAX-1);
                                        }
                                } else {
                                        fprintf(stderr,"%d parnell_translate: no valid destination directory %s\n", MyRank, dir_name);
                                        status = PARNELL_ERROR; goto exit;
                                }
                        } else {
                                perror ("cannot stat directory");
                                fprintf(stderr,"%d parnell_translate: cannot get status of destination directory %s\n", MyRank, dir_name);
                                status = PARNELL_ERROR; goto exit;
                        }
                }
        } else {
                strncpy (tmp_name, dst_name, FILENAME_MAX-1);
                dir_name = dirname (tmp_name);
                if (stat (dir_name, &dst_info) == 0) {
                        if (S_ISDIR (dst_info.st_mode)) {
                                if (dst_info.st_ino == wrk_info.st_ino) {
                                        strncpy (tmp_name, dst_name, FILENAME_MAX-1);
                                        dst_name[FILENAME_MAX-1]='\0';
                                        strncpy (dst_name, basename (tmp_name), FILENAME_MAX-1);
                                }
                        } else {
                                fprintf(stderr,"%d parnell_translate: no valid destination directory %s\n", MyRank, dir_name);
                                status = PARNELL_ERROR; goto exit;
                        }
                } else {
                        perror ("cannot stat directory");
                        fprintf(stderr,"%d parnell_translate: cannot get status of destination directory %s\n", MyRank, dir_name);
                        status = PARNELL_ERROR; goto exit;
                }
        }
        status = PARNELL_CONTINUE;

 exit:
        free(tmp_name);

        return status;
}
