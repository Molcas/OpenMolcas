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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:49:45 stevenv"
 *
 *       File:         parnell_replica.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_replica - performs a simple copy operation
 */

#include "parnell.h"

parnell_status_t
parnell_replica (char * src_name, char * dst_name)
{
        parnell_status_t status = PARNELL_START;
        FILE* src_file = NULL;
        FILE* dst_file = NULL;
        struct stat src_info, dst_info;
        size_t bytes_read, bytes_written;

        /* check if source and destination file are different */
        if (stat (src_name, &src_info) != 0) {
                if (errno != EOVERFLOW) {
                        perror("cannot stat file");
                        fprintf(stderr,"%d parnell_replica: cannot get status of source %s\n", MyRank, src_name);
                        return PARNELL_ERROR;
                }
        } else if (!S_ISREG (src_info.st_mode)) {
                fprintf(stderr,"%d parnell_replica: not a regular source file %s\n", MyRank, src_name);
                return PARNELL_ERROR;
        }
        /* source is OK, proceed */
        if(stat (dst_name, &dst_info) != 0) {
                if (errno != ENOENT && errno != EOVERFLOW) {
                        perror("cannot stat file");
                        fprintf(stderr,"%d parnell_replica: cannot handle status of destination %s\n", MyRank, dst_name);
                        return PARNELL_ERROR;
                }
        } else if (!S_ISREG (dst_info.st_mode)) {
                fprintf(stderr,"%d parnell_replica: not a regular destination file %s\n", MyRank, dst_name);
                return PARNELL_ERROR;
        } else {
                if (src_info.st_ino == dst_info.st_ino) return PARNELL_OK; /* src = dst, no work to do */
        }

        /* open source and destination */
        if ((src_file = fopen (src_name, "r")) == NULL) {
                perror("cannot open file for reading");
                fprintf(stderr,"%d parnell_replica: error opening source file %s\n", MyRank, src_name);
                return PARNELL_ERROR;
        }
        if ((dst_file = fopen (dst_name, "w")) == NULL) {
                perror("cannot open file for writing");
                fprintf(stderr,"%d parnell_replica: error opening destination file %s\n", MyRank, dst_name);
                fclose (src_file);
                return PARNELL_ERROR;
        }

        /* start copy operation */
        void * buffer = (void*)malloc (PARNELL_BUFSIZE);
        bytes_read = PARNELL_BUFSIZE;
        while (bytes_read) {
                bytes_read = fread (buffer, 1, PARNELL_BUFSIZE, src_file);
                if (bytes_read != PARNELL_BUFSIZE && feof(src_file) == 0) {
                        perror("premature end while reading");
                        fprintf(stderr,"%d parnell_replica: cannot read from source file %s\n", MyRank, src_name);
                        status = PARNELL_ERROR; goto exit;
                } else {
                        bytes_written = fwrite (buffer, 1, bytes_read, dst_file);
                        if (bytes_written != bytes_read) {
                                perror("premature end while writing");
                                fprintf(stderr,"%d parnell_replica: cannot write to destination file %s\n", MyRank, dst_name);
                                status = PARNELL_ERROR; goto exit;
                        }
                }
        }
        status = PARNELL_OK;

 exit:
        free (buffer);

        fclose (src_file);
        fclose (dst_file);

        return status;
}
