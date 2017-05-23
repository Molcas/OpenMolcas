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

#ifdef _SOLARIS_
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#endif
#include "molcastype.h"

#ifdef _CAPITALS_
#define fndlnk FNDLNK
#else
#ifndef ADD_
#define fndlnk fndlnk_
#endif
#endif


     void fndlnk(irc,fnin,fnout)
     char *fnin,*fnout;
     INT *irc;
    {
#ifdef _SOLARIS_
     char fntmp[256];
     struct stat tmp;
     char *ptr;
     ptr=strchr(fnin,' ');
      if(ptr!=NULL) *ptr=0;
     *irc=lstat(fnin,&tmp);
     strcpy(fnout,fnin);

     if ((tmp.st_mode & S_IFMT) == S_IFLNK )
     {

      INT length=100;
      INT fnlng=readlink(fnin,fnout,length);
      fnout[fnlng]='\0';
      strcpy(fntmp,fnout);
      fndlnk(irc,fntmp,fnout);
     }
     *irc=lstat(fnout,&tmp);
#endif

    }
