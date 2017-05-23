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
#include <time.h>
#include <assert.h>
#include "molcastype.h"
#ifndef _WIN32_
#include <unistd.h>
#endif

#ifdef _CAPITALS_
#define get_genome GET_GENOME
#else
#ifndef ADD_
#define get_genome get_genome_
#endif
#endif

#define LEN_DNA 256
#define LEN_HOST 64
#define LEN_TIME 64

void get_genome(char * cDNA, INT *nDNA)
   {
    char hostname[LEN_HOST];
    gethostname(hostname, sizeof(hostname));
    int pid = getpid();
    time_t t = time(NULL);
    char sTime[LEN_TIME];
    int i;

    i=strftime(sTime, LEN_TIME, "%c", localtime(&t));
    assert (i < LEN_TIME);
    for (i=0;i<LEN_DNA;i++) cDNA[i] = ' ';
    *nDNA=snprintf (cDNA, LEN_DNA, "HOST %s PID %d DATE %s", hostname, pid, sTime);
    assert (*nDNA < LEN_DNA);
    *nDNA=8*((strlen(cDNA)+7)/8);
    assert (*nDNA < LEN_DNA);
    cDNA[strlen(cDNA)]=' ';
   }
