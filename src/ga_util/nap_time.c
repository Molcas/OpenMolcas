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
#ifndef _WIN32_
#include <unistd.h>
#endif

#include "molcastype.h"

#ifdef _MOLCAS_MPP_
#include "ga.h"
#endif

#ifdef _CAPITALS_
#define nap_time NAP_TIME
#else
#ifndef ADD_
#define nap_time nap_time_
#endif
#endif

char* getenvc(char*);

void nap_time()
   {
#ifndef _MOLCAS_MPP_
   return;
#else
    char hostname[256];
    int MyRank;
    char *mysleep;
    int sleep_time;
    if((mysleep=getenvc("MOLCAS_NAP"))==NULL) return;
    sscanf(mysleep,"%d",&sleep_time);
    if(sleep_time>0)
    {
    gethostname(hostname, sizeof(hostname));
    MyRank=GA_Nodeid();
    printf("PID %d on %s with rank %d ready for attach\n", getpid(), hostname,MyRank);
    fflush(stdout);
    sleep(sleep_time);
    }
#endif
   }
