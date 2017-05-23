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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:44:05 stevenv"
 *
 *       File:         parnell_scatter.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_scatter - takes a file with "src_name" and copies it to
 *                          "dst_name" in the main work directory "dst_name" on
 *                          the master and to the subdirectories of all the
 *                          slaves, thus distributing the file.
 *
 */

#include "parnell.h"

parnell_status_t
parnell_scatter (char * src_name, char * dst_name)
{
        parnell_status_t status = PARNELL_START;

        void* buffer;
        buffer = (void*)malloc (PARNELL_BUFSIZE);

        /* start MPI copy operation: master reads source file and transmits it,
         * slaves receive transmitted data and write to their own subdirectory */
        if (MyRank == 0) {
                FILE* src_file = NULL;
                size_t bytes_read;
                /* make local copy on the master */
                status = parnell_replica (src_name, dst_name);

                if (nProcs == 1) goto exit;

                /* open the source file */
                if (status == PARNELL_OK) {
                        if ((src_file = fopen (src_name, "r")) == NULL) {
                                perror("cannot open file for reading");
                                fprintf(stderr,"%d parnell_scatter: error opening source file %s\n", MyRank, src_name);
                                status = PARNELL_ERROR;
                        } else {
                                status = PARNELL_CONTINUE;
                        }
                }
#ifdef _MOLCAS_MPP_
                /* communicate to slaves to continue or stop */
                MPI_Bcast (&status, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
                if (status != PARNELL_CONTINUE) goto exit;
                /* read source and broadcast to slaves */
                bytes_read = 1;
                while (bytes_read) {
                        bytes_read = fread (buffer, 1, PARNELL_BUFSIZE, src_file);
                        if (bytes_read != PARNELL_BUFSIZE && feof(src_file) == 0) {
                                perror("premature end while reading");
                                fprintf(stderr,"%d parnell_scatter: error reading source file %s\n", MyRank, src_name);
                                status = PARNELL_ERROR; goto exit;
                        }
#ifdef _MOLCAS_MPP_
                        MPI_Bcast (&bytes_read, 1, MPI_LONG, 0, MPI_COMM_WORLD);
                        if (bytes_read) {
                                MPI_Bcast (buffer, bytes_read, MPI_BYTE, 0, MPI_COMM_WORLD);
                        }
#endif
                }
                fclose (src_file);
                status = PARNELL_OK;
        } else {
#ifdef _MOLCAS_MPP_
                FILE* dst_file = NULL;
                size_t bytes_received;
                size_t bytes_written;
                MPI_Bcast (&status, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
                if (status != PARNELL_CONTINUE) goto exit;
                /* open destination file */
                if ((dst_file = fopen (dst_name, "w")) == NULL) {
                        perror("cannot open file for writing");
                        fprintf(stderr,"%d parnell_scatter: error opening destination file %s\n", MyRank, dst_name);
                        status = PARNELL_ERROR; goto exit;
                }
                /* receive buffer and write to disk */
                bytes_received = 1;
                while (bytes_received) {
                        MPI_Bcast (&bytes_received, 1, MPI_LONG, 0, MPI_COMM_WORLD);
                        if (bytes_received) {
                                MPI_Bcast (buffer, bytes_received, MPI_BYTE, 0, MPI_COMM_WORLD);
                                bytes_written = fwrite (buffer, 1, bytes_received, dst_file);
                                if (bytes_written != bytes_received) {
                                        perror("premature end while writing");
                                        fprintf(stderr,"%d parnell_scatter: error writing destination file %s\n", MyRank, dst_name);
                                        status = PARNELL_ERROR; goto exit;
                                }
                        }
                }
                status = PARNELL_OK;
#endif
        }

 exit:
        free(buffer);

        return status;
}
