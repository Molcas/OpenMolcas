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
*               2012,2013, Victor P. Vysotskiy                         *
*               2020, Ignacio Fdez. Galvan                             *
***********************************************************************/
/******************************************************************************/
/*                                                                            */
/*                               A I X - I / O                                */
/*                                                                            */
/*  The fast I/O calls the following C-langue primitives:                     */
/*  open, close, read, write, lseek, remove and fsync                         */
/*  This file includes the FORTRAN to C-language interfaces.                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    written by:                                                             */
/*    M.P. Fuelscher                                                          */
/*    University of Lund, Sweden, 1992                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    history: Victor P. Vysotskiy, University of Lund, Sweden, 2012-2013     */
/*            - collect all declarations in the 'cio.h' header file           */
/*            - thread-safe pwrite/pread IO                                   */
/*            - filesize information via 'c_stat'                             */
/*            Ignacio Fdez. Galvan, 2020                                      */
/*            - add c_rename                                                  */
/*                                                                            */
/******************************************************************************/

#include <fcntl.h>
#ifndef _WIN32_
#include <unistd.h>
#else
#include <windows.h>
#define open  _lopen
#define close _lclose
#define read  _lread
#define write _lwrite
#define lseek _llseek
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/uio.h>
#include "molcastype.h"
#include "cio.h"

#define MIN(x,y) (x<y? x : y)
/*--- c_open -----------------------------------------------------------------*/
INT c_open(Path)
 char *Path;

{

 INT rc;
 INT oFlag;
 INT oMode;
#ifdef _CRAY_C90_
 char fn[256];
 oFlag=O_CREAT|O_RDWR;
 (void)strcpy(fn,Path);
 rc = open(fn,oFlag,0644);
#else
#ifndef _WIN32_
 oFlag=O_CREAT|O_RDWR;
 oMode=S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR;
 rc = open(Path,oFlag,oMode);
#else
 oFlag=OF_READWRITE;
 rc=open(Path,oFlag);
#endif
#endif
 if(rc<0) {
   oFlag=O_RDONLY;
#ifdef _CRAY_C90_
   rc = open(fn,oFlag);
#else
#ifndef _WIN32_
   rc = open(Path,oFlag);
#else
   oFlag=OF_READ;
   rc=open(Path,oFlag);
#endif
#endif
 }

 return rc;

}
/*--- c_open_w -----------------------------------------------------------------*/
INT c_openw(Path)
 char *Path;

{

 INT rc;
 INT oFlag;
 INT oMode;
#ifdef _CRAY_C90_
 char fn[256];
 oFlag=O_CREAT|O_RDWR|O_TRUNC;
 (void)strcpy(fn,Path);
 rc = open(fn,oFlag,0644);
#else
#ifndef _WIN32_
 oFlag=O_CREAT|O_RDWR|O_TRUNC;
 oMode=S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR;
 rc = open(Path,oFlag,oMode);
#else
 oFlag=OF_READWRITE;
 rc=open(Path,oFlag);
#endif
#endif
 return rc;

}

/*--- c_close ----------------------------------------------------------------*/
INT c_close(FileDescriptor)
 INT *FileDescriptor;

{
 INT rc;
 rc = close(*FileDescriptor);
 return rc;
}

/*--- c_read -----------------------------------------------------------------*/
INT c_read(FileDescriptor,Buffer,nBytes)
 INT *FileDescriptor;
 char *Buffer;
 INT *nBytes;

{
 INT rc=0;
 INT bfrblk=1024*1024;
 INT i=0;
 INT remains;
 INT readlength;
 remains=*nBytes;
 while (remains > 0){
      readlength = MIN(bfrblk,remains);
      rc = (INT)read(*FileDescriptor,(void *)(Buffer+i),(size_t)(readlength));
      if ( rc == readlength ) { i = i+readlength; rc = i; remains = remains - bfrblk;}
      else { rc = 0; return rc ;}
 }
 return rc;
}


/*--- c_pread -----------------------------------------------------------------*/
INT c_pread(INT *FileDescriptor,char *Buffer,INT *nBytes,INT *Offset) {
 INT rc=0;
 rc=(INT) pread((int) *FileDescriptor,(void *) Buffer, (size_t) *nBytes, (off_t) *Offset);
 return rc;
}



/*--- c_write ----------------------------------------------------------------*/
INT c_write(FileDescriptor,Buffer,nBytes)
 INT *FileDescriptor;
 char *Buffer;
 INT *nBytes;

{
 INT rc=0;
 INT bfrblk=1024*1024;
 INT i=0;
 INT remains;
 INT writelength;
 remains=*nBytes;
 while (remains > 0){
      writelength = MIN(bfrblk,remains);
      rc = (INT)write(*FileDescriptor,(void *)(Buffer+i),(size_t)(writelength));
      if ( rc == writelength ) { i = i+writelength; rc = i; remains = remains - bfrblk;}
      else { rc = 0; return rc ;}
 }
 return rc;
}


/*--- c_pwrite ----------------------------------------------------------------*/
INT c_pwrite(INT *FileDescriptor,char *Buffer,INT *nBytes, INT *Offset) {
INT rc=0;

 rc = (INT) pwrite((int) *FileDescriptor,(void *)(Buffer),(size_t)(*nBytes), (off_t) *Offset);
 return rc;
}


/*--- c_lseek ----------------------------------------------------------------*/
INT c_lseek(FileDescriptor,Offset)
 INT *FileDescriptor;
 INT *Offset;

{
#ifdef _WIN32_
typedef long off_t;
#endif
 INT rc;
 rc = (INT)lseek(*FileDescriptor,(off_t)(*Offset),SEEK_SET);
 return rc;
}

/*--- c_remove ---------------------------------------------------------------*/
INT c_remove(FileName)
 char *FileName;

{
 INT rc;
#ifdef _CAPITALS_
 char fn[256];
#endif

#ifdef _CAPITALS_
 (void)strcpy(fn,FileName);
 rc = remove(fn);
#else
#ifndef _WIN32_
 rc = remove(FileName);
#else
 rc = DeleteFile(FileName);
#endif
#endif
 return rc;
}

/*--- c_fsync ----------------------------------------------------------------*/
INT c_fsync(FileDescriptor)
 INT *FileDescriptor;

{
 INT rc;
#ifndef _WIN32_
 rc = fsync(*FileDescriptor);
#else
 rc=0;
#endif
 return rc;
}

/*--- c_copy ----------------------------------------------------------------*/
INT c_copy(FileDescriptor1, FileDescriptor2)
 INT *FileDescriptor1, *FileDescriptor2;
{
 INT rc;
 char *Buffer;
 struct stat stat;
 size_t rce;

 rc=fstat(*FileDescriptor1, &stat);

 rce=stat.st_size;
 Buffer=(char*) malloc(sizeof(char)*(rce+1));
      rc = (INT)read(*FileDescriptor1,(void *)(Buffer),(size_t)(rce));
      rc = (INT)write(*FileDescriptor2,(void *)(Buffer),(size_t)(rce));
 free(Buffer);
 return rc;
}

/*--- c_stat ----------------------------------------------------------------*/
INT c_stat(FileDescriptor)
 INT *FileDescriptor;
{
 INT rc;
 struct stat flstat;
 off_t fsize;

 rc=fstat(*FileDescriptor, &flstat);
 (void)rc;
 fsize=flstat.st_size;
 return fsize;
}

/*--- c_rename --------------------------------------------------------------*/
INT c_rename(FileName,NewName)
 char *FileName;
 char *NewName;

{
 INT rc;
#ifdef _CAPITALS_
 char fn[256];
 char nn[256];
#endif

#ifdef _CAPITALS_
 (void)strcpy(fn,FileName);
 (void)strcpy(nn,FileName);
 rc = rename(fn,nn);
#else
#ifndef _WIN32_
 rc = rename(FileName,NewName);
#else
 rc = MoveFile(FileName,NewName);
#endif
#endif
 return rc;
}
