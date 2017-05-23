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

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:52:02 stevenv"
 *
 *       File:         main.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell - a program for handling copying and deleting of files
 *                 in a parallel environment - main
 *
 */

#include "parnell.h"

int MyRank = 0;
int nProcs = 1;

int main (int argc, char * argv[])
{
        parnell_status_t status = PARNELL_START;

#ifdef _MOLCAS_MPP_
        char *molcas_nprocs = getenv("MOLCAS_NPROCS");
        if ((molcas_nprocs == NULL) || (atoi(molcas_nprocs) != 1)) {
                MPI_Init(&argc, &argv);
                MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
                MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
        }
#endif

        status = parnell_cmd(argc, argv);

        if (status != PARNELL_OK) {
                fprintf (stderr, "%d parnell: ABORTING\n", MyRank); fflush (NULL);
#ifdef _MOLCAS_MPP_
                if (nProcs > 1) MPI_Abort (MPI_COMM_WORLD, 1);
#endif
                return 1;
        } else {
#ifdef _MOLCAS_MPP_
                if (nProcs > 1) MPI_Finalize();
#endif
                return 0;
        }
}
