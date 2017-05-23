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
* Copyright (C) Markus P. Fuelscher                                    *
*               Luis Serrano-Andres                                    *
***********************************************************************/
/******************************************************************************/
/*                                                                            */
/*  objective:                                                                */
/*  This routine searches for the environment variable MOLCASDISK which       */
/*  determines the largest file size for a multiple file dataset. Multiple    */
/*  file datasets are used to cut huge files into several files of smaller    */
/*  size and thus allows to overcome the 2 GByte limit common to many         */
/*  variants of UNIX.                                                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*  The max size of a file is defined by the environment variable             */
/*  MOLCASDISK, and the size is expresses in units of megabytes.              */
/*  The following three cases may ocurr                                       */
/*  1. MOLCASDISK is not defined -- The multiple file I/O option is           */
/*     ignored (default).                                                     */
/*  2. MOLCASDISK is set to 0 -- A max. file size of 2048 MB is assumed.      */
/*  3. MOLCASDISK is set to a nonzero value -- The max file size is set to    */
/*     this value.                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*  Usage: filelen=allocdisk()                                                */
/*  filelen - Max file size in megabytes.                                     */
/*            Type: integer.                                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Authors: M. P. Fuelscher and Luis Serrano-Andres                           */
/*          University of Lund                                                */
/*          Sweden                                                            */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "molcastype.h"

#ifdef _CAPITALS_
#define allocdisk ALLOCDISK
#else
#ifndef ADD_
#define allocdisk allocdisk_
#endif
#endif

char* getenvc(const char*);


INT allocdisk()

{

   char *disksize;
   char *ptr;
   INT MBytes;

   disksize=getenvc("MOLCAS_DISK");
   if(disksize==NULL) {
      fprintf(stdout,"grabit: MOLCAS_DISK is not defined \n");
      MBytes=0;
   } else {
      MBytes=0;
      ptr=disksize;
      while(*disksize!='\0') {
         if(isdigit(*disksize)) MBytes=10*MBytes+*disksize-'0';
         disksize++;
      }
      free(ptr);
#ifdef _I8_
      if(MBytes==0) MBytes=204700;
#else
      if(MBytes==0) MBytes=2047;
#endif
#ifdef _CRAY_C90_
      MBytes=0;
#endif
/*
      fprintf(stdout,"grabit: MOLCASDISK is ");
      fprintf(stdout,INT_FORMAT,MBytes);
      fprintf(stdout," MBytes \n");
*/
   }

   return(MBytes);

}
