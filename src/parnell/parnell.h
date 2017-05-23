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

/* -*- mode: C -*- Time-stamp: "2010-07-03 11:49:39 stevenv"
 *
 *       File:         parnell.h
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell -- program for handling copying and deleting of files in a
 *                  parallel environment - header file
 *
 *       History:
 *       V.P. Vysotskiy, University of Lund, Sweden, 2013
 *       Modified for a bunch of subsequent calls
 */

#ifndef PARNELL_H
#define PARNELL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libgen.h>
#include <limits.h>
#include <dirent.h>

#include <sys/stat.h>
#include <sys/types.h>
#ifndef _WIN32_
#  include <unistd.h>
#endif

#ifdef _MOLCAS_MPP_
#  include <mpi.h>
#endif

#include "parnell_status.h"

/* macro definitions */

#ifndef FILENAME_MAX
#  define FILENAME_MAX 4096
#endif

#define PARNELL_BUFSIZE 4096
#define PARNELL_TIMEOUT 5

/* global variables */

extern int MyRank;
extern char MyWorkDir[FILENAME_MAX];
extern int nProcs;
extern char WorkDir[FILENAME_MAX];

/* function prototypes */
parnell_status_t
parnell_cmd (int argc, char* argv[]);

parnell_status_t
parnell (int argc, char* argv[]);

parnell_status_t
parnell_base (int argc, char **argv);

parnell_status_t
parnell_init (void);

parnell_status_t
parnell_copy (int argc, char **argv);

parnell_status_t
parnell_translate (char * src_name, char * dst_name);

parnell_status_t
parnell_scatter (char * src_name, char * dst_name);

parnell_status_t
parnell_collect (char * src_name, char * dst_name);

parnell_status_t
parnell_replica (char * src_name, char * dst_name);

parnell_status_t
parnell_reduce (char * src_name, char * dst_name);

parnell_status_t
parnell_rmlist(char * rmlist);

parnell_status_t
parnell_wipe (void);

parnell_status_t
parnell_unlink (char *fpath);

parnell_status_t
parnell_exec (int argc, char ** argv);

#endif /* PARNELL_H */
