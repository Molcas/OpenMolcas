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
#include <ctype.h>
#include <molcastype.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#ifdef _CAPITALS_
#define open_molcas_info OPEN_MOLCAS_INFO
#else
#ifndef ADD_
#define open_molcas_info open_molcas_info_
#endif
#endif

#ifdef _CAPITALS_
#define close_molcas_info CLOSE_MOLCAS_INFO
#else
#ifndef ADD_
#define close_molcas_info close_molcas_info_
#endif
#endif

#ifdef _CAPITALS_
#define add_molcas_info ADD_MOLCAS_INFO
#else
#ifndef ADD_
#define add_molcas_info add_molcas_info_
#endif
#endif

FILE *f;

int open_molcas_info()
{
struct stat buffer;
if(stat("molcas_info",&buffer)!=0)
{
f=fopen("molcas_info","w");
fprintf(f,"%s","###########\n# MOLCAS-Info_File Vers.No. 1.2\n###########\n");
}
else
{
f=fopen("molcas_info","a");
}
return 0;
}

int close_molcas_info()
{
fclose(f);
return 0;
}

int add_molcas_info(char *str, INT *n)
{
str[*n]=0;
fprintf(f,"%s\n",str);
return 0;
}
