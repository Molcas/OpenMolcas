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
* Copyright (C) 1992, Markus P. Fuelscher                              *
***********************************************************************/
/**********************************************************************/
/*                                                                    */
/*    Extract the current date and time.                              */
/*                                                                    */
/*    Note:                                                           */
/*    The VS-FORTRAN subrotuines, Datim and Datimx, are replaced      */
/*    by this routine.                                                */
/*                                                                    */
/*--------------------------------------------------------------------*/
/*                                                                    */
/*    written by:                                                     */
/*    M.P. Fuelscher                                                  */
/*    University of Lund, Sweden, 1992                                */
/*                                                                    */
/*--------------------------------------------------------------------*/
/*                                                                    */
/*    history: none                                                   */
/*                                                                    */
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#ifndef _WIN32_
#include <sys/time.h>
#endif
#include <string.h>

#include "molcastype.h"

#ifdef _CAPITALS_
#define datimx DATIMX
#else
#ifndef ADD_
#define datimx datimx_
#endif
#endif


void datimx(char *TimeStamp)
{
#ifdef _WIN32_
strcpy(TimeStamp,"Once upon a time..");
#else
  static int CTIME_RES_LENGTH = 24;
  char *ptr;
  time_t x;
  struct timeval t;
  struct timezone tz;
  if ( gettimeofday(&t,&tz) != 0 ) {
       printf(" *** Error in procedure datimx: %s\n",strerror(errno));
       exit(20);
     }
     else
     {
       x=(time_t)t.tv_sec;
       ptr=ctime(&x);
       if(ptr!=NULL) {
        strncpy(TimeStamp,ptr, CTIME_RES_LENGTH);
        TimeStamp[CTIME_RES_LENGTH+1]=0;
        }
     }
#endif
}
