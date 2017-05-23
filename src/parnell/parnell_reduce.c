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

/* -*- mode: C -*- Time-stamp: "2015-07-28 17:56:43 stevenv"
 *
 *       File:         parnell_reduce.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_reduce - reads a number from the file with name "src_name"
 *                        on each slave and write the global maximum to
 *                         "dst_name" on the master (use e.g. for return code)
 *
 */

#include "parnell.h"

parnell_status_t
parnell_reduce (char * src_name, char * dst_name)
{
        int status = PARNELL_START;
#ifdef _MOLCAS_MPP_
        int status_send;
#endif

        void* buffer;
        buffer = (void*)malloc (PARNELL_BUFSIZE);

        FILE* src_file = NULL;
        FILE* dst_file = NULL;
        size_t bytes_read;
        long int lrc;
#ifdef _MOLCAS_MPP_
        long int lrc_send;
#endif
        char* endptr;

        /* try to open source file */
        if ((src_file = fopen (src_name, "r")) == NULL) {
                perror("cannot open file for reading");
                fprintf(stderr,"%d parnell_reduce: cannot open source file %s\n", MyRank, src_name);
                status = PARNELL_ERROR;
        } else {
                /* start reading from disk and transmit to the master */
                bytes_read = fread (buffer, 1, PARNELL_BUFSIZE, src_file);
                if (bytes_read != PARNELL_BUFSIZE && feof(src_file) == 0) {
                        perror("premature end while reading");
                        fprintf(stderr,"%d parnell_reduce: cannot read from source file %s\n", MyRank, src_name);
                        status = PARNELL_ERROR;
                } else if (bytes_read == PARNELL_BUFSIZE) {
                        fprintf(stderr,"%d parnell_reduce: max buffer size reached wile reading %s\n", MyRank, src_name);
                        status = PARNELL_ERROR;
                } else {
                        /* strtol converts a string to an integer _with_ error checking (atoi does no error checking) */
                        errno = 0;
                        lrc = strtol(buffer, &endptr, 10);
                        if ( (errno == ERANGE && (lrc == LONG_MIN || lrc == LONG_MAX)) \
                             || (errno != 0 && lrc == 0) ) {
                                perror("strtol");
                                fprintf(stderr,"%d parnell_reduce: invalid digits in source %s\n", MyRank, src_name);
                                status = PARNELL_ERROR;
                        } else if (endptr == buffer) {
                                fprintf(stderr,"%d parnell_reduce: no digits in source %s\n", MyRank, src_name);
                                status = PARNELL_ERROR;
                        } else {
                                status = PARNELL_CONTINUE;
                        }
                }
        }

#ifdef _MOLCAS_MPP_
        status_send = status;
        if (nProcs > 1) MPI_Allreduce (&status_send, &status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

        if (status != PARNELL_CONTINUE) goto exit;

#ifdef _MOLCAS_MPP_
        lrc_send = lrc;
        if (nProcs > 1) MPI_Reduce (&lrc_send, &lrc, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
#endif

        /* write rc to file on master */
        if (MyRank == 0) {
                if ((dst_file = fopen (dst_name, "w")) == NULL) {
                        perror("cannot open file for writing");
                        fprintf(stderr,"%d parnell_reduce: cannot open destination file %s\n", MyRank, dst_name);
                        status = PARNELL_ERROR;
                        goto exit;
                }
                fprintf(dst_file, "%ld\n", lrc);
        }

        status = PARNELL_OK;

 exit:
        free(buffer);

        return (parnell_status_t)status;
}
