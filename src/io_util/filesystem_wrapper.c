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
* Copyright (C) 2019, Oskar Weser                                      *
***********************************************************************/

#define _XOPEN_SOURCE 500
#include <ftw.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "molcastype.h"

/* C_SIZE_T (or in general unsigned ints) is not supported by FORTRAN */
/* Implicit cast to c_int */
INT strlen_wrapper(const char*const* str)
{
  return strlen(*str);
}

void getcwd_wrapper(char* path, const INT* n, INT* err)
{
  if (getcwd(path, *n) == path) {
    *err = 0;
    INT i = -1;
/* This is necessary for FORTRAN trim() to work correctly.*/
    while (path[++i] != '\0');
    for (; i < *n; i++) {
      path[i] = ' ';
    }
  } else {
    *err = 1;
  }
}

void chdir_wrapper(const char* path, INT *err)
{
  *err = chdir(path);
}

void symlink_wrapper(const char* to, const char* from, INT* err)
{
  *err = symlink(to, from);
}

/* MODE_T (or in general unsigned ints) is not supported by FORTRAN */
/* Implicit cast to c_int */
void mkdir_wrapper(const char*  path, const INT* mode, INT* err)
{
    *err = mkdir(path, *mode);
}

INT get_errno() {
  return errno;
}


/* Method to recursively remove a directory and its contents
   (from a stackoverflow.com answer) */

/* private method to be used as argument for rmrf */
static int unlink_cb(const char* fpath, const struct stat* sb,
                     int typeflag, struct FTW* ftwbuf)
{
    int rv = remove(fpath);

    (void)sb;
    (void)typeflag;
    (void)ftwbuf;

    if (rv)
        perror(fpath);

    return rv;
}

void remove_wrapper(const char* path, INT* err)
{
    *err = (INT) nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}
