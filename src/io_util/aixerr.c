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
* Copyright (C) 1993, Per-Olof Widmark                                 *
***********************************************************************/
/******************************************************************************/
/*                                                                            */
/*                               A I X - I / O                                */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* rc=AixErr(string)                                                          */
/*                                                                            */
/* This routine returns the last error return code to the caller. String also */
/* contains the error message. The string is of fixed length 80 for easy      */
/* compatibillity with fortran.                                               */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Author:  Per-Olof Widmark                                                  */
/*          IBM Sweden                                                        */
/* Written: December 1993                                                     */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* History: none                                                              */
/*                                                                            */
/******************************************************************************/

#include <errno.h>
#include <string.h>
#include "molcastype.h"

#ifdef _CAPITALS_
#define aixerr AIXERR
#else
#ifndef ADD_
#define aixerr aixerr_
#endif
#endif

INT aixerr(s)
 char s[];

{
   char *p;
   INT k;
   k=0;
   if( errno>0 ) {
      p=strerror(errno);
      while ( (k<80) && (*p!=0) ) s[k++]=*(p++);
   } else {
      strcpy(s,"Unknown error");
/*    to avoid conflicts with the FORTRAN function of the same name           */
/*    replace the following line by an explicit number:                       */
/*    k=strlen("Unknown error");                                              */
      k=13;
   }
   while ( k<80 ) s[k++]=' ';
   return(errno);
}
