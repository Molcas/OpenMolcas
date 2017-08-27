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
* Copyright (C) 2012, Steven Vancoillie                                *
***********************************************************************/

/* -*- mode: C -*- Time-stamp: "2010-07-02 14:52:33 stevenv"
 *
 *       File:         parnell_copy.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_copy - handle copying of files in parallel according to the mode:
 *
 * 0 = replica: regular copy, but only executed on the master node, typical
 *              for output files (emil SAVE)
 * 1 = scatter: source file on the master is distributed to the work directory
 *              of all the slaves, typical for input files (emil COPY)
 * 2 = collect: source files on all slaves are gathered on the master with
 *              with a rank suffix (emil COLLECT)
 * 3 = replica: regular copy, executed on all the nodes (emil CLONE)
 *
 */

#include "parnell.h"

parnell_status_t
parnell_copy (int argc, char **argv)
{
        parnell_status_t status = PARNELL_START;

        int mode;
        char src_name[FILENAME_MAX];
        char dst_name[FILENAME_MAX];

        /* read the arguments and get the proper relative file names */
        if (argc != 3) {
                fprintf (stderr, "parnell_copy: expecting 3 arguments (mode source dest):\n");
                for (int i=0; i<argc; i++) fprintf(stderr, " %s", argv[i]);
                fprintf(stderr, "\n");
                status = PARNELL_ERROR;
        } else {
                mode = argv[0][0];
                if (MyRank == 0) {
                        strncpy (src_name, argv[1], FILENAME_MAX);
                        strncpy (dst_name, argv[2], FILENAME_MAX);
                        status = parnell_translate (src_name, dst_name);
                }
#ifdef _MOLCAS_MPP_
                if (nProcs > 1) MPI_Bcast (&status, sizeof(status), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
        }

        /* decide to continue or stop */
        if (status != PARNELL_CONTINUE) goto exit;

#ifdef _MOLCAS_MPP_
        if (nProcs > 1) {
                MPI_Bcast (src_name, FILENAME_MAX, MPI_BYTE, 0, MPI_COMM_WORLD);
                MPI_Bcast (dst_name, FILENAME_MAX, MPI_BYTE, 0, MPI_COMM_WORLD);
        }
#endif
        /* perform the actual copy action depending on the mode */
        switch (mode) {
        case '0':
                if (MyRank == 0) {
                        status = parnell_replica (src_name, dst_name);
                } else {
                        status = PARNELL_OK;
                }
                break;
        case '1':
                status = parnell_scatter (src_name, dst_name);
                break;
        case '2':
                status = parnell_collect (src_name, dst_name);
                break;
        case '3':
                status = parnell_replica (src_name, dst_name);
                break;
        case '4':
                status = parnell_reduce  (src_name, dst_name);
                break;
        default:
                fprintf (stderr, "%d parnell_copy: invalid mode number: %d\n", MyRank, mode);
                status = PARNELL_ERROR;
                break;
        }

 exit:
        return status;
}
