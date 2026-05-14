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
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "molcastype.h"

/* C_SIZE_T (or in general unsigned ints) is not supported by FORTRAN */
/* Implicit cast to c_int */
INT strlen_wrapper(const char *const *str) {
  return strlen(*str);
}

INT access_wrapper(const char *path) {
  /* https://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c */
  return access(path, F_OK);
}

void getcwd_wrapper(char *path, const INT *n, INT *err) {
  if (getcwd(path, *n) == path) {
    *err = 0;
    INT i = -1;
    /* This is necessary for FORTRAN trim() to work correctly.*/
    while (path[++i] != '\0') {
    }
    for (; i < *n; i++) {
      path[i] = ' ';
    }
  } else {
    *err = 1;
  }
}

void chdir_wrapper(const char *path, INT *err) {
  *err = chdir(path);
}

void symlink_wrapper(const char *to, const char *from, INT *err) {
  *err = symlink(to, from);
}

/* MODE_T (or in general unsigned ints) is not supported by FORTRAN */
/* Implicit cast to c_int */
void mkdir_wrapper(const char *path, const INT *mode, INT *err) {
  *err = mkdir(path, *mode);
}

INT get_errno() {
  return errno;
}

/* Method to recursively remove a directory and its contents
(from a stackoverflow.com answer) */

/* private method to be used as argument for nftw */
static int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf) {
  int rv = remove(fpath);

  (void)sb;
  (void)typeflag;
  (void)ftwbuf;

  if (rv)
    perror(fpath);

  return rv;
}

void remove_wrapper(const char *path, INT *err) {
  *err = (INT)nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}

void copy(const char *src, const char *dst, INT *err) {
  char buf[BUFSIZ];
  size_t size;

  *err = 0;
  FILE *source = fopen(src, "rb");

  if (!source) {
    *err = 1;
    return;
  }

  FILE *dest = fopen(dst, "wb");

  // feof(FILE* stream) returns non-zero if the end of file indicator for stream is set

  while ((size = fread(buf, 1, BUFSIZ, source))) {
    fwrite(buf, 1, size, dest);
  }

  fclose(source);
  fclose(dest);
}
