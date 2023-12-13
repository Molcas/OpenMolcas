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
*               2023, Ignacio Fdez. Galvan                             *
***********************************************************************/

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:38:16 stevenv"
 *
 *       File:         parnell_wipe.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_wipe - delete all files from the work directory
 *
 */

#include "parnell.h"

parnell_status_t parnell_wipe(int argc, char **argv) {
  struct stat info;
  struct dirent *entry;
  DIR *cwd_ptr;
  char tmpdir[16];
  int skip;
  char *flag;

  if (argc > 0) {
    flag = argv[0];
  } else {
    flag = "none";
  }

  if (!(cwd_ptr = opendir(MyWorkDir))) {
    perror("parnell_wipe: error trying to open work directory");
    fprintf(stderr, "%d parnell_wipe: could not wipe %s\n", MyRank, MyWorkDir);
    return PARNELL_ERROR;
  }

  while ((entry = readdir(cwd_ptr))) {
    /* skip dot entries */
    if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, ".."))
      continue;
    if (lstat(entry->d_name, &info) == 0) {
      /* skip slave directories in the master */
      skip = 1;
      if (MyRank == 0 && S_ISDIR(info.st_mode)) {
        for (int i = 1; i < nProcs; i++) {
          snprintf(tmpdir, 15, "tmp_%d", i);
          if (!strcmp(entry->d_name, tmpdir)) {
            skip = 0;
            break;
          }
        }
      }
      if (!skip)
        continue;
      if (S_ISREG(info.st_mode) || S_ISLNK(info.st_mode)||S_ISDIR(info.st_mode))
        parnell_unlink(entry->d_name);
    } else {
      /* if error other than "No such file or directory", report it */
      if (errno != ENOENT)
        perror("parnell_wipe: error while calling stat on file");
    }
  }
  closedir(cwd_ptr);

  /* remove the WorkDir */
  if (!strcmp(flag, "remove")) {
    /* tmp_n subdir in slaves */
    if (MyRank > 1 && lstat(MyWorkDir, &info) == 0) {
      if (S_ISDIR(info.st_mode))
        remove(MyWorkDir);
    }
    if (lstat(WorkDir, &info) == 0) {
      if (S_ISDIR(info.st_mode))
        remove(MyWorkDir);
    }
  }

  return PARNELL_OK;
}
