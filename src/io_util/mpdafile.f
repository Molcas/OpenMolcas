************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Markus P. Fuelscher                              *
*               1996, Luis Serrano-Andres                              *
*               2012, Victor P. Vysotskiy                              *
*               2016, Steven Vancoillie                                *
************************************************************************
      Subroutine MpDaFile(Lu,MaxFileSizel,iOpt,Buf,lBuf,iDisk)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     File systems or disks may be to small to keep a given large      *
*     dataset. Therefore, it may be split into subsets which then      *
*     may be linked to different physical units. This process is       *
*     controlled by the environment variable MOLCASDISK                *
*     (c.f. allocdisk). This routine computes the proper splitting     *
*     and calls the IO routines.                                       *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
*     MaxFileSize : integer, input                                     *
*               max length of a file.                                  *
*     iOpt    : integer, input                                         *
*               option code                                            *
*               iOpt= 0 dummy write                                    *
*               iOpt= 1 synchronous write                              *
*               iOpt= 2 synchronous read                               *
*               iOpt= 3 synchronous gather and write                   *
*               iOpt= 4 synchronous read and scatter                   *
*               iOpt= 5 synchronous rewind                             *
*               iOpt= 6 asynchronous write                             *
*               iOpt= 7 asynchronous read                              *
*               iOpt= 8 asynchronous gather and write                  *
*               iOpt= 9 asynchronous read and scatter                  *
*               iOpt=10 asynchronous rewind                            *
*               Note: At present the asynchronous modes are not        *
*                     supported and work identically the synchronous   *
*                     modes                                            *
*     Buf     : array of characters, input/output                      *
*               Buffer carrying/receiving the data to write/read       *
*     lBuf    : integer, input                                         *
*               length of the buffer Buf in bytes!!                    *
*     iDisk   : integer, input/output                                  *
*               disk address as byte offset!!!                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher, University of Lund, Sweden, 1996                 *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
*     Steven Vancoillie: use of byte lengths/offests (2016)            *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Character Buf(*)

#include "SysDef.fh"
#include "blksize.fh"
#include "filesize.fh"
#include "fio.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#endif
      Character*8 Stdnam, ext
      Character*100 My_Progname,Get_Progname
      External Get_Progname
      Integer Length_Progname
      Character*256 tmp
      Character*80 Text
      Character*16 TheName
      Data TheName/'MpDaFile'/

      My_Progname=Get_Progname()
      Length_Progname=StrnLn(My_Progname)
      If ( Query ) Call qEnter(TheName)
*
      max_File_Size = MaxFileSizel*10**6
      max_Bytes     = MIN(max_File_Length,max_File_Size)
      n_Bytes       = lBuf
      offset        = iDisk/max_Bytes
*  set correct char to append: 0-9, A-Z
           start_char=48

         if(offset.le.9) then
           start_char=48
         else
           start_char=55
         endif
c
      If (offset.lt.0 .or. offset.gt.MaxSplitFile-1 ) Then
      StdNam=LuName(Lu)
      Call PrgmTranslate(StdNam,tmp,ltmp)
         Write (6,*) '          Current I/O Status as follows'
         Write (6,*)
         iDum=iPrintLevel(3)
         Call FASTIO('STATUS')
       Call SysAbendFileMsg(TheName,StdNam,
     &   'Extensions out of range!',
     &   'increase MOLCAS_DISK value or MaxSplitFile in fio.fh ')

         Call Abend()
      End If
      LU_mod    = MPUnit(offset,LU)
      iDisk_mod = iDisk-offset*max_Bytes
*
*---- If extension not opened before: open and initiate here!
*
      StdNam=LuName(Lu)
      Call PrgmTranslate(StdNam,tmp,ltmp)
      lStdNam=ltmp
      If (Lu_Mod.lt.0) Then
         Lu_Mod=isfreeunit(Lu)
         MPUnit(offset,LU)=Lu_Mod
         temp = 0
         tmp(1+lStdNam:1+lStdNam)=Char(start_char+offset)
         iext=strnln(StdNam)
         ext=StdNam
         if(offset.le.9) then
          ext(1+iext:1+iext)=Char(start_char+offset)
         else
          ext(1+iext:1+iext)=Char(start_char+offset/10)
          ext(2+iext:2+iext)=Char(start_char+offset-offset/10*10)
         endif
*
         iRc = AixOpn(temp,tmp,.false.)
         If ( iRc.ne.0 ) then
            iRc = AixErr(Text)
            Call SysFileMsg(TheName,'MSG: open', Lu_Mod,Text)
         End If
         isOpen(Lu_Mod)    = 1
         FSCB(Lu_Mod)    = temp
*
cvv         LuName(Lu_Mod)    = tmp
         LuName(Lu_Mod)    = ext
         Call SetLuMark(Lu_Mod)
         Addr(Lu_Mod)    = 0
         Multi_File(Lu_Mod)=.True.
         MPUnit(0,LU_Mod)=Lu
#ifdef _OLD_IO_STAT_
         MxAddr(Lu_Mod)  = 0
         Count(1,Lu_Mod) = 0
         Count(2,Lu_Mod) = 0
         Count(3,Lu_Mod) = 0
         Count(4,Lu_Mod) = 0
#endif
         MBL(Lu_Mod)     = MBL(Lu)
*
      End If
*
      If ( (iDisk_mod+n_Bytes).le.max_Bytes ) then
*
*------- The whole record in contained on a single file.
*
         p_Buf     = 1
         n_Copy    = lBuf
         CALL CHDAFILE(LU_mod,iOpt,Buf(p_Buf),n_Copy,iDisk_mod)
      Else
*
*------- The record must be split over several files.
*
         p_Buf     = 1
         n_Save    = lBuf
         n_Copy    = (max_Bytes-iDisk_mod)
         Do while( n_Save.gt.0 )
            If (LU_MOD.lt.0) Then
               Lu_Mod=isfreeunit(Lu)
               MPUnit(offset,LU)=Lu_Mod
               temp = 0
           start_char=48
         if(offset.le.9) then
           start_char=48
         else
           start_char=55
         endif

               tmp(1+lStdNam:1+lStdNam)=Char(start_char+offset)
         iext=strnln(StdNam)
         ext=StdNam
         if(offset.le.9) then
           ext(1+iext:1+iext)=Char(start_char+offset)
         else
          ext(1+iext:1+iext)=Char(start_char+offset/10)
          ext(2+iext:2+iext)=Char(start_char+offset-offset/10*10)
         endif
*
               iRc = AixOpn(temp,tmp,.false.)
               If ( iRc.ne.0 ) then
                  iRc = AixErr(Text)
                  Call SysFileMsg(TheName,'MSG: open', Lu_Mod,Text)
               End If
               isOpen(Lu_mod)    = 1
               FSCB(Lu_mod)    = temp
*
cvv               LuName(Lu_Mod)    = tmp
               LuName(Lu_Mod)    = ext
               Call SetLuMark(Lu_Mod)

               Addr(Lu_Mod)    = 0
               Multi_File(Lu_Mod)=.True.
               MPUnit(0,LU_Mod)=Lu
#ifdef _OLD_IO_STAT_
               MxAddr(Lu_Mod)  = 0
               Count(1,Lu_Mod) = 0
               Count(2,Lu_Mod) = 0
               Count(3,Lu_Mod) = 0
               Count(4,Lu_Mod) = 0
#endif
               MBL(Lu_Mod)     = MBL(Lu)
*
            End If
            CALL CHDAFILE(LU_mod,iOpt,Buf(p_Buf),n_Copy,iDisk_mod)
            p_Buf     = p_Buf+n_Copy
            n_Save    = n_Save-n_Copy
            n_Copy    = Min(max_Bytes,n_Save)
            offset = offset + 1
      If (offset.gt.MaxSplitFile-1 ) Then
         Write (6,*) '          Current I/O Status as follows'
         Write (6,*)
         iDum=iPrintLevel(3)
         Call FASTIO('STATUS')
       Call SysAbendFileMsg(TheName,StdNam,
     &   'Extensions out of range!',
     &   'increase MOLCAS_DISK value or MaxSplitFile in fio.fh ')

         Call Abend()
      End If

            Lu_Mod=MPUnit(offset,LU)
            iDisk_mod = 0
         End Do
      End If

      If ( Query ) Call qExit(TheName)

      Return
      End
