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
* Copyright (C) 2001-2005, Valera Veryazov                             *
***********************************************************************/
/*
 *  getenv
 *
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
#ifdef _MOLCAS_MPP_
#  include <ga.h>
#endif
*/

#include "molcastype.h"

char* getenvc(const char*);

#ifdef _CAPITALS_
#define getenvf2c GETENVF2C
#else
#ifndef ADD_
#define getenvf2c getenvf2c_
#endif
#endif

#ifdef _CAPITALS_
#define getenvinit GETENVINIT
#else
#ifndef ADD_
#define getenvinit getenvinit_
#endif
#endif

/* MAXSTR = max len of ENV value */
/* MAXENV = max len of MOLCAS related ENV in total */
/* I do NOT care about memory leak in this code! */

#define MAXSTR 256
#define MAXENV 4096
static char MOLCAS_ENV[MAXENV];


void getenvf2c(char *name, INT *ilen, char *value, INT *maxlen, INT *irl)
   {
   char *envvar;
   INT len=0;
   char *ptr;
   char *name0;
    name0=(char*) malloc( (*ilen+1)*sizeof(char));

     if(!name0) return;
    strncpy(name0,name,*ilen);
    name0[*ilen]=0;

    ptr=strchr(name0,' ');
    if(ptr)  *ptr=0;
   if((envvar=getenvc(name0))!=NULL)
      {
      len=strlen(envvar);
      if (len>=*maxlen) len=*maxlen-1;
      strncpy(value,envvar,len);
      value[len]=0;
      free(envvar);
      }
     *irl=len;
     free(name0);
  return;
   }
char *getenvc(const char *name)
   {

char Name[MAXSTR];
int i;
char *value;
char *ptr;
char *ptr2;
 Name[0]='\n';
 Name[1]=0;
 i=strlen(name);
   if(i>MAXSTR-2)
   {
        fprintf(stderr,"Environment variable %s is too long!\n",name);
        return NULL;
   }
 strcat(Name,name);
 strcat(Name,"=");
 /* debug
 puts("env=");
 puts(MOLCAS_ENV);

 */
   ptr=strstr(MOLCAS_ENV,Name);

   if(ptr==NULL)
   {
/* trying to get env */
       ptr=getenv(name);
       if(ptr!=NULL)
       {
           value=(char *)malloc((strlen(ptr)+1) *sizeof(char));
           strcpy(value,ptr);
           return value;
       }
       return NULL;
   }
   ptr=ptr+i+2;
   ptr2=strchr(ptr,'\n');
   if(ptr2==NULL)
   {
        fprintf(stderr,"Environment variable %s is not terminated!\n",name);
        return NULL;
   }
   i=ptr2-ptr;
   if(i > MAXSTR)
   {
        fprintf(stderr,"Environment value for %s is too long!\n",name);
        return NULL;
   }
   value=(char *)malloc((i+1) *sizeof(char));
   if(value==NULL) return NULL;
    strncpy(value,ptr,i);
    value[i]=0;
    return value;
}

/* a routine reads molcas.env file and returns 0/-1 */
INT getenvinit(void)
{
FILE *MYENV;
char line[MAXSTR];
int i=0;
   /* system("pwd"); */
   MYENV=fopen("molcas.env","r");
   if (MYENV==NULL)
           {
             fprintf(stderr,"Unable to open molcas.env file\n");

           return -1;
           }

        MOLCAS_ENV[0]='\n';
        MOLCAS_ENV[1]=0;
   while(!feof(MYENV))
   {
           if(fgets(line,MAXSTR,MYENV)==NULL) continue;
           if(line[0]!='#')
           {
                   line[MAXSTR-1]=0;
                   i=i+strlen(line);
                   if(i>MAXENV) return -1;
                   strcat(MOLCAS_ENV,line);
           }
   }
   fclose(MYENV);
/*  puts(MOLCAS_ENV);  */
 return 0;

}
