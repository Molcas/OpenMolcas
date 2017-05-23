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

/* -*- mode: C -*- Time-stamp: "2010-07-03 16:38:14 stevenv"
 *
 *       File:         parnell_collect.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_collect - copies a file with name "src_name" on each slave to
 *                         "dst_name.$n" on the master, where $n is the rank of
 *                         the process, thereby collecting them - needs rigorous
 *                         error-checking
 *
 */

#include "parnell.h"

parnell_status_t
parnell_collect (char * src_name, char * dst_name)
{
        parnell_status_t status = PARNELL_START;

        void* buffer;
        buffer = (void*)malloc (PARNELL_BUFSIZE);

        if (MyRank != 0) {
                FILE* src_file = NULL;
                size_t bytes_read;
                /* try to open source file */
                if ((src_file = fopen (src_name, "r")) == NULL) {
                        perror("cannot open file for reading");
                        fprintf(stderr,"%d parnell_collect: cannot open source file %s\n", MyRank, src_name);
                        status = PARNELL_ERROR;
                } else {
                        status = PARNELL_CONTINUE;
                }
#ifdef _MOLCAS_MPP_
                MPI_Send (&status, sizeof(parnell_status_t), MPI_BYTE, 0, MyRank, MPI_COMM_WORLD);
#endif
                if (status == PARNELL_CONTINUE) {
                        /* start reading from disk and transmit to the master */
                        bytes_read = PARNELL_BUFSIZE;
                        while (bytes_read) {
                                bytes_read = fread (buffer, 1, PARNELL_BUFSIZE, src_file);
                                if (bytes_read != PARNELL_BUFSIZE && feof(src_file) == 0) {
                                        perror("premature end while reading");
                                        fprintf(stderr,"%d parnell_collect: cannot read from source file %s\n", MyRank, src_name);
                                        status = PARNELL_ERROR;
                                } else {
                                        status = PARNELL_CONTINUE;
                                }
#ifdef _MOLCAS_MPP_
                                MPI_Send (&status, sizeof(parnell_status_t), MPI_BYTE, 0, MyRank, MPI_COMM_WORLD);
#endif
                                if (status == PARNELL_CONTINUE) {
#ifdef _MOLCAS_MPP_
                                        MPI_Send (&bytes_read, sizeof(size_t), MPI_BYTE, 0, MyRank, MPI_COMM_WORLD);
                                        if (bytes_read) {
                                                MPI_Send (buffer, bytes_read, MPI_BYTE, 0, MyRank, MPI_COMM_WORLD);
                                        }
#endif
                                } else {
                                        break;
                                }
                        }
                        fclose(src_file);
                        if (status == PARNELL_CONTINUE) status = PARNELL_OK;
                } else {
                        goto exit;
                }
        } else {
                /* take master copy */
                status = parnell_replica (src_name, dst_name);

#ifdef _MOLCAS_MPP_
                if (nProcs == 1) goto exit;

                FILE* dst_file = NULL;
                size_t bytes_received;
                size_t bytes_written;
                /* collect from slaves */
                char * file_name = (char*)malloc (FILENAME_MAX);

                int iRank, iAlive;
                int n_slaves = nProcs - 1;
                int n_seconds = 0;
                parnell_status_t * status_slave = (parnell_status_t*)malloc (sizeof(parnell_status_t)*nProcs);

                MPI_Status Stat;
                MPI_Request * Req = (MPI_Request*)malloc (sizeof(MPI_Request)*nProcs);

                for (iRank=1; iRank<nProcs; ++iRank) {
                        status_slave[iRank] = PARNELL_START;
                        MPI_Irecv (&status_slave[iRank], sizeof(int), MPI_BYTE, iRank, iRank, MPI_COMM_WORLD, &Req[iRank]);
                }

                /* wait for data from slaves */
                while (n_slaves > 0 && n_seconds < PARNELL_TIMEOUT) {
                        /* check all slaves for what to do */
                        for (iRank=1; iRank<nProcs; ++iRank) {
                                MPI_Test (&Req[iRank], &iAlive, &Stat);
                                if (iAlive && status_slave[iRank] == PARNELL_CONTINUE) {
                                        /* open destination file */
                                        snprintf (file_name, FILENAME_MAX, "%s.%d", dst_name, iRank);
                                        if ((dst_file = fopen (file_name, "w")) == NULL) {
                                                perror("cannot open file for writing");
                                                fprintf(stderr,"%d parnell_collect: cannot open destination file %s\n", MyRank, file_name);
                                                status_slave[iRank] = PARNELL_ERROR;
                                        } else {
                                                /* start receiving from slave iRank and write to disk */
                                                bytes_received = PARNELL_BUFSIZE;
                                                while (bytes_received && status_slave[iRank] == PARNELL_CONTINUE) {
                                                        MPI_Recv ((status_slave + iRank), sizeof(parnell_status_t), MPI_BYTE, iRank, iRank, MPI_COMM_WORLD, &Stat);
                                                        if (status_slave[iRank] == PARNELL_CONTINUE) {
                                                                MPI_Recv (&bytes_received, 1, MPI_LONG, iRank, iRank, MPI_COMM_WORLD, &Stat);
                                                                if (bytes_received) {
                                                                        MPI_Recv (buffer, bytes_received, MPI_BYTE, iRank, iRank, MPI_COMM_WORLD, &Stat);
                                                                        bytes_written = fwrite (buffer, 1, bytes_received, dst_file);
                                                                        if (bytes_written != bytes_received) {
                                                                                perror("premature end while writing");
                                                                                fprintf(stderr,"%d parnell_collect: cannot write to destination file %s\n", MyRank, file_name);
                                                                                status_slave[iRank] = PARNELL_ERROR;
                                                                        }
                                                                }
                                                        }
                                                }
                                                fclose(dst_file);
                                                if (status_slave[iRank] == PARNELL_CONTINUE) status_slave[iRank] = PARNELL_FINISHED;
                                        }
                                        --n_slaves;
                                }
                        }
                        sleep (1); ++n_seconds;
                }
                for (iRank=1; iRank<nProcs; ++iRank) {
                        if (status_slave[iRank] != PARNELL_FINISHED) {
                                fprintf(stderr,"%d parnell_collect: problem with process %d\n", MyRank, iRank);
                        }
                }
                free(status_slave);
                free(Req);
                free(file_name);
#endif
        }

 exit:
        free(buffer);

        return status;
}
