/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molcastype.h"

#ifndef _WIN32_
  #include <unistd.h>
#endif

#ifdef _CAPITALS_
  #define slaveschdir SLAVESCHDIR
#else
  #ifndef ADD_
    #define slaveschdir slaveschdir_
  #endif
#endif

/* This should be enough to cover, if not increase :) */
#define SRANK_SIZE 8

void slaveschdir(INT *MyRank, INT *iRC)
{
  *iRC = 0;
  if (*MyRank != 0) {
    char s_rank[SRANK_SIZE];
    char * p_WorkDir = (char*)malloc (sizeof(char)*FILENAME_MAX);
    char * p_subdir = (char*)malloc (sizeof(char)*FILENAME_MAX);
    if (snprintf (s_rank, SRANK_SIZE, INT_FORMAT, *MyRank) == -1) {
      perror ("while printing to string");
      fprintf(stderr,"%s slaveschdir: fatal error, could not write rank to string\n", s_rank);
      *iRC = 99; goto exit;
    }
    if (getcwd (p_WorkDir, FILENAME_MAX) == NULL) {
      perror ("while calling getcwd");
      fprintf(stderr,"%s slaveschdir: fatal error, could not determine current working directory\n", s_rank);
      *iRC = 99; goto exit;
    }
    sprintf (p_subdir, "%s/tmp_%s", p_WorkDir, s_rank);
    if (chdir (p_subdir) != 0) {
      perror ("cannot change directory");
      fprintf(stderr,"%s slaveschdir: fatal error, could not switch to directory %s\n", s_rank, p_subdir);
      *iRC = 99; goto exit;
    }
exit:
    free (p_WorkDir);
    free (p_subdir);
    return;
  }
  return;
}
