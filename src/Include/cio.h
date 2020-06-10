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
* Copyright (C) 2012-2013, Victor P. Vysotskiy                         *
*               2020, Ignacio Fdez. Galvan                             *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/*              THE NATIVE MOLCAS'S AIX IO LAYER                          */
/*                                                                        */
/* Just a header file                                                     */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Victor P. Vysotskiy                                           */
/*          Lund University, Sweden                                       */
/* Written: 2012-2013                                                     */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* History:                                                               */
/*                                                                        */
/**************************************************************************/
#ifdef _CAPITALS_
#define c_open   C_OPEN
#define c_openw  C_OPENW
#define c_close  C_CLOSE
#define c_read   C_READ
#define c_pread  C_PREAD
#define c_write  C_WRITE
#define c_pwrite C_PWRITE
#define c_lseek  C_LSEEK
#define c_remove C_REMOVE
#define c_fsync  C_FSYNC
#define c_copy   C_COPY
#define c_stat   C_STAT
#define c_rename C_RENAME
#else
#ifndef ADD_
#define c_open   c_open_
#define c_openw  c_openw_
#define c_close  c_close_
#define c_read   c_read_
#define c_pread  c_pread_
#define c_write  c_write_
#define c_pwrite c_pwrite_
#define c_lseek  c_lseek_
#define c_remove c_remove_
#define c_fsync  c_fsync_
#define c_copy   c_copy_
#define c_stat   c_stat_
#define c_rename c_rename_
#endif
#endif

INT c_open(char *Path);
INT c_openw(char *Path);
INT c_close(INT *FileDescriptor);
INT c_read(INT *FileDescriptor,char *Buffer,INT *nBytes);
INT c_pread(INT *FileDescriptor,char *Buffer,INT *nBytes,INT *Offset);
INT c_write(INT *FileDescriptor,char *Buffer,INT *nBytes);
INT c_pwrite(INT *FileDescriptor,char *Buffer,INT *nBytes, INT *Offset);
INT c_lseek(INT *FileDescriptor,INT *Offset);
INT c_remove(char *FileName);
INT c_fsync(INT *FileDescriptor);
INT c_copy(INT *FileDescriptor1, INT *FileDescriptor2);
INT c_stat(INT *FileDescriptor);
INT c_rename(char *FileName,char *NewName);
