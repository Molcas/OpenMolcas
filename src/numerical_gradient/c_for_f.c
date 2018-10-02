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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
***********************************************************************/
/* C wrappers to be called from Fortran (newdir.f) */

#define _XOPEN_SOURCE 500
#include <ftw.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "molcastype.h"

/* Method to create a directory
   (uses int* instead of mode_t as argument, because mode_t not portable with Fortran) */
INT mkdir_for_f(const char *pathname, INT *mode)
{
    INT rv = mkdir(pathname, *mode);
    return rv;
}

/* Method to get the current directory
   (uses int* instead of size_t as argument, because size_t not portable with Fortran) */
char *getcwd_for_f(char *buf, INT *size)
{
    char *cwd = getcwd(buf, *size);
    return cwd;
}

/* Method to change directory
   (for consistency) */
INT chdir_for_f(const char *path)
{
    INT rv = chdir(path);
    return rv;
}

/* Method to recursively remove a directory and its contents
   (from a stackoverflow.com answer) */

/* private method to be used as argument for rmrf */
static int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
    int rv = remove(fpath);

    (void)sb;
    (void)typeflag;
    (void)ftwbuf;

    if (rv)
        perror(fpath);

    return rv;
}

INT rmrf(char *path)
{
    return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}
