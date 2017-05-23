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
#include <stdlib.h>
#ifndef _WIN32_
#include <unistd.h>
#endif
#include <time.h>
#ifndef _WIN32_
#include <pwd.h>
#endif
#include <sys/types.h>
#include <string.h>

#include "molcastype.h"

char* getenvc(const char*);

#ifdef _CAPITALS_
#define unixinfoc UNIXINFOC
#else
#ifndef ADD_
#define unixinfoc unixinfoc_
#endif
#endif

void unixinfoc(INT *pid, INT *ppid,
   INT *sec, INT *min, INT *hour, INT *mday, INT *mon, INT *year, INT *wday, INT *yday, INT *isdst,
   char *username, char *realname, char *homedir, char *shell,
   char *molcasdir)
   {
#ifndef _WIN32_
   time_t timep;
   struct tm *structtimep;
#ifdef _PLEASE_DELETE_ME_
   uid_t  uId;
   struct passwd *user;
#endif
   int i;
   char *envvar;
/* Process ID : */
   *pid=(INT)getpid();
   *ppid=(INT)getppid();
/* Information on date and time : */
   timep=time(NULL);
   structtimep=localtime(&timep);
   *sec=(INT)structtimep->tm_sec;
   *min=(INT)structtimep->tm_min;
   *hour=(INT)structtimep->tm_hour;
   *mday=(INT)structtimep->tm_mday;
   *mon=(INT)structtimep->tm_mon;
   *year=(INT)structtimep->tm_year;
   *wday=(INT)structtimep->tm_wday;
   *yday=(INT)structtimep->tm_yday;
   *isdst=(INT)structtimep->tm_isdst;
/* Login name */
#ifdef _PLEASE_DELETE_ME_
   uId=geteuid();
   user=getpwuid(uId);
   if(user) {
     for(i=0;i<(int)strlen(user->pw_name);i++)  username[i]=user->pw_name[i];
     for(i=0;i<(int)strlen(user->pw_gecos);i++) realname[i]=user->pw_gecos[i];
     for(i=0;i<(int)strlen(user->pw_dir);i++)   homedir[i]=user->pw_dir[i];
     for(i=0;i<(int)strlen(user->pw_shell);i++) shell[i]=user->pw_shell[i];
   }
#endif
/* Environment variables */
/* Molcas directory */
   if((envvar=getenvc("MOLCAS"))!=NULL)
      {
      for(i=0;i<(int)strlen(envvar);i++)
         molcasdir[i]=envvar[i];
      free(envvar);
      }
#endif
   }
